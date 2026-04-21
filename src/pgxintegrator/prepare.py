from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
import gzip
import shutil
import subprocess

from .config import AppConfig, load_config
from .gene_registry import ensure_gene_coordinates, normalize_build
from .manifests import (
    DatasetRecord,
    SampleRecord,
    load_dataset_manifest,
    load_gene_list,
    load_sample_manifest,
    list_gene_names,
    summarize_chrom_datasets,
    summarize_gene_datasets,
)


INTERNAL_PADDING_BP = 10000


class PreparePlanError(ValueError):
    """Raised when a prepare plan cannot be built."""


class PrepareExecutionError(RuntimeError):
    """Raised when the prepare stage fails during execution."""


@dataclass
class PlatformInputPlan:
    platform: str
    source_path: Path
    source_build: str
    layout: str
    liftover_needed: bool
    chrom: Optional[str] = None
    gene: Optional[str] = None
    dataset_id: Optional[str] = None
    extract_chrom: Optional[str] = None
    extract_start: Optional[int] = None
    extract_end: Optional[int] = None


@dataclass
class GenePreparePlan:
    gene: str
    chrom: str
    start: int
    end: int
    extract_start: int
    extract_end: int
    target_build: str
    platform_inputs: Dict[str, PlatformInputPlan] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)


@dataclass
class PreparePlan:
    project_name: str
    target_build: str
    registry_tsv: Path
    output_dir: Path
    sample_count: int
    gene_count: int
    gene_plans: List[GenePreparePlan]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def build_prepare_plan(
    cfg: AppConfig,
    registry_tsv: Optional[Union[str, Path]] = None,
) -> PreparePlan:
    """
    Build a prepare plan.

    This function does not run bcftools, liftover, or integration. It only:
      - loads manifests,
      - resolves gene names to coordinates,
      - matches genes to candidate input files by platform,
      - resolves source-build extraction coordinates per platform.
    """
    samples = load_sample_manifest(cfg.inputs.sample_manifest)
    datasets = load_dataset_manifest(cfg.inputs.dataset_manifest)
    gene_list_records = load_gene_list(cfg.inputs.gene_list)

    requested_genes = list_gene_names(gene_list_records)
    target_build = normalize_build(cfg.target_build)
    registry_path = _resolve_registry_tsv(cfg, registry_tsv)

    target_gene_rows = ensure_gene_coordinates(
        genes=requested_genes,
        build=target_build,
        registry_tsv=str(registry_path),
    )
    target_gene_map = _rows_to_gene_lookup(target_gene_rows)

    source_builds = sorted(set(normalize_build(record.build) for record in datasets))
    per_build_gene_maps = {target_build: target_gene_map}

    for source_build in source_builds:
        if source_build in per_build_gene_maps:
            continue
        source_rows = ensure_gene_coordinates(
            genes=requested_genes,
            build=source_build,
            registry_tsv=str(registry_path),
        )
        per_build_gene_maps[source_build] = _rows_to_gene_lookup(source_rows)

    chrom_datasets = summarize_chrom_datasets(datasets)
    gene_datasets = summarize_gene_datasets(datasets)

    gene_plans: List[GenePreparePlan] = []
    for gene in requested_genes:
        target_row = target_gene_map.get(gene)
        if target_row is None:
            raise PreparePlanError(
                "Could not resolve target-build coordinates for gene {}.".format(gene)
            )

        chrom = _normalize_chrom(str(target_row["chrom"]))
        start = int(target_row["start"])
        end = int(target_row["end"])

        extract_start = max(1, start - INTERNAL_PADDING_BP)
        extract_end = end + INTERNAL_PADDING_BP

        gene_plan = GenePreparePlan(
            gene=gene,
            chrom=chrom,
            start=start,
            end=end,
            extract_start=extract_start,
            extract_end=extract_end,
            target_build=target_build,
        )

        for platform in _platforms_in_order(datasets):
            selected = _select_dataset_for_gene(
                platform=platform,
                gene=gene,
                chrom=chrom,
                gene_datasets=gene_datasets,
                chrom_datasets=chrom_datasets,
            )

            if selected is None:
                gene_plan.warnings.append(
                    "No {} dataset found for gene {}".format(platform, gene)
                )
                continue

            source_build = normalize_build(selected.build)
            source_gene_map = per_build_gene_maps.get(source_build)
            if source_gene_map is None or gene not in source_gene_map:
                raise PreparePlanError(
                    "Could not resolve {}-build coordinates for gene {}.".format(
                        source_build,
                        gene,
                    )
                )

            source_row = source_gene_map[gene]
            source_chrom = _normalize_chrom(str(source_row["chrom"]))
            source_start = int(source_row["start"])
            source_end = int(source_row["end"])

            platform_extract_start = max(1, source_start - INTERNAL_PADDING_BP)
            platform_extract_end = source_end + INTERNAL_PADDING_BP

            gene_plan.platform_inputs[platform] = PlatformInputPlan(
                platform=platform,
                source_path=selected.path,
                source_build=source_build,
                layout=selected.layout,
                liftover_needed=(source_build != target_build),
                chrom=selected.chrom,
                gene=selected.gene,
                dataset_id=selected.dataset_id,
                extract_chrom=source_chrom,
                extract_start=platform_extract_start,
                extract_end=platform_extract_end,
            )

        if not gene_plan.platform_inputs:
            gene_plan.warnings.append(
                "No usable datasets found for gene {} on chromosome {}".format(
                    gene,
                    chrom,
                )
            )

        gene_plans.append(gene_plan)

    return PreparePlan(
        project_name=cfg.project_name,
        target_build=target_build,
        registry_tsv=registry_path,
        output_dir=cfg.prepare.output_dir,
        sample_count=len(samples),
        gene_count=len(gene_plans),
        gene_plans=gene_plans,
    )


def run_prepare_dry_run(
    config_path: Union[str, Path],
    registry_tsv: Optional[Union[str, Path]] = None,
) -> PreparePlan:
    cfg = load_config(config_path)
    return build_prepare_plan(cfg, registry_tsv=registry_tsv)


def run_prepare(
    config_path: Union[str, Path],
    registry_tsv: Optional[Union[str, Path]] = None,
    dry_run: bool = False,
) -> PreparePlan:
    """
    Execute the prepare stage.

    Workflow per gene:
      1. subset padded gene region and samples from each selected platform dataset
         using coordinates in that platform's source build
      2. rename samples to canonical sample_id values
      3. liftover with UCSC liftOver if needed
      4. harmonize chromosome naming to match the target reference FASTA
      5. normalize
      6. integrate platforms into one per-gene VCF
    """
    cfg = load_config(config_path)
    plan = build_prepare_plan(cfg, registry_tsv=registry_tsv)

    if dry_run:
        return plan

    samples = load_sample_manifest(cfg.inputs.sample_manifest)
    _execute_prepare_plan(cfg, plan, samples)
    return plan


def format_prepare_plan(plan: PreparePlan) -> str:
    lines: List[str] = []
    lines.append("Project: {}".format(plan.project_name))
    lines.append("Target build: {}".format(plan.target_build))
    lines.append("Sample count: {}".format(plan.sample_count))
    lines.append("Gene count: {}".format(plan.gene_count))
    lines.append("Registry TSV: {}".format(plan.registry_tsv))
    lines.append("Prepare output dir: {}".format(plan.output_dir))
    lines.append("Internal padding: {} bp".format(INTERNAL_PADDING_BP))
    lines.append("")

    for gene_plan in plan.gene_plans:
        lines.append(
            "[{}] raw={}:{}-{} | extract={}:{}-{}".format(
                gene_plan.gene,
                gene_plan.chrom,
                gene_plan.start,
                gene_plan.end,
                gene_plan.chrom,
                gene_plan.extract_start,
                gene_plan.extract_end,
            )
        )

        if gene_plan.platform_inputs:
            for platform in sorted(gene_plan.platform_inputs.keys()):
                item = gene_plan.platform_inputs[platform]
                lines.append(
                    "  - {} | {} | build={} | liftover={} | source={} | extract={}:{}-{}".format(
                        item.platform,
                        item.layout,
                        item.source_build,
                        "yes" if item.liftover_needed else "no",
                        item.source_path,
                        item.extract_chrom,
                        item.extract_start,
                        item.extract_end,
                    )
                )
        else:
            lines.append("  - no platform inputs")

        for warning in gene_plan.warnings:
            lines.append("  ! {}".format(warning))

        lines.append("")

    return "\n".join(lines).rstrip() + "\n"


# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------


def _execute_prepare_plan(
    cfg: AppConfig,
    plan: PreparePlan,
    samples: List[SampleRecord],
) -> None:
    output_dir = cfg.prepare.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    work_root = output_dir / "_work"
    work_root.mkdir(parents=True, exist_ok=True)

    sample_root = output_dir / "_samples"
    sample_root.mkdir(parents=True, exist_ok=True)

    for gene_plan in plan.gene_plans:
        final_vcf = output_dir / "{}.vcf.gz".format(gene_plan.gene)
        final_index = Path(str(final_vcf) + ".csi")

        if cfg.prepare.resume and final_vcf.exists() and final_index.exists():
            continue

        if not gene_plan.platform_inputs:
            raise PrepareExecutionError(
                "Cannot prepare gene {} because no usable platform inputs were found.".format(
                    gene_plan.gene
                )
            )

        gene_work_dir = work_root / gene_plan.gene
        gene_work_dir.mkdir(parents=True, exist_ok=True)

        ready_inputs: Dict[str, Path] = {}
        for platform, platform_plan in sorted(gene_plan.platform_inputs.items()):
            sample_ids_file, rename_map_file = _write_sample_files(
                samples=samples,
                platform=platform,
                output_dir=sample_root,
            )

            region_bed = gene_work_dir / "{}.{}.bed".format(gene_plan.gene, platform)
            extracted_vcf = gene_work_dir / "{}.{}.extract.vcf.gz".format(gene_plan.gene, platform)
            renamed_vcf = gene_work_dir / "{}.{}.renamed.vcf.gz".format(gene_plan.gene, platform)
            lifted_vcf = gene_work_dir / "{}.{}.lifted.vcf.gz".format(gene_plan.gene, platform)
            normalized_vcf = gene_work_dir / "{}.{}.norm.vcf.gz".format(gene_plan.gene, platform)

            _write_platform_bed_file(platform_plan=platform_plan, bed_path=region_bed)

            _extract_platform_region(
                platform_plan=platform_plan,
                sample_ids_file=sample_ids_file,
                region_bed=region_bed,
                output_vcf=extracted_vcf,
            )
            _index_vcf(extracted_vcf)

            _rename_samples(
                input_vcf=extracted_vcf,
                rename_map_file=rename_map_file,
                output_vcf=renamed_vcf,
            )
            _index_vcf(renamed_vcf)

            current_vcf = renamed_vcf
            if platform_plan.liftover_needed:
                if not cfg.liftover.enabled:
                    raise PrepareExecutionError(
                        "Gene {} platform {} requires liftover from {} to {}, but liftover is disabled.".format(
                            gene_plan.gene,
                            platform,
                            platform_plan.source_build,
                            gene_plan.target_build,
                        )
                    )
                _liftover_vcf_with_ucsc(
                    input_vcf=current_vcf,
                    output_vcf=lifted_vcf,
                    chain_file=cfg.liftover.chain_file,
                )
                _index_vcf(lifted_vcf)
                current_vcf = lifted_vcf

            _normalize_vcf(
                input_vcf=current_vcf,
                output_vcf=normalized_vcf,
                reference_fasta=cfg.reference.fasta,
            )
            _index_vcf(normalized_vcf)
            ready_inputs[platform] = normalized_vcf

        if "wes" in ready_inputs and "array" in ready_inputs:
            _integrate_gene_vcfs(
                wes_vcf=ready_inputs["wes"],
                array_vcf=ready_inputs["array"],
                output_vcf=final_vcf,
                temp_dir=gene_work_dir,
            )
        elif "wes" in ready_inputs:
            _copy_with_index(ready_inputs["wes"], final_vcf)
        elif "array" in ready_inputs:
            _copy_with_index(ready_inputs["array"], final_vcf)
        else:
            raise PrepareExecutionError(
                "No prepared platform VCFs available for gene {}.".format(gene_plan.gene)
            )

        _index_vcf(final_vcf)

        if not cfg.prepare.keep_intermediate:
            _cleanup_gene_work_dir(gene_work_dir)


# ---------------------------------------------------------------------------
# VCF processing helpers
# ---------------------------------------------------------------------------


def _write_platform_bed_file(platform_plan: PlatformInputPlan, bed_path: Path) -> None:
    """
    Write a BED file for source-build extraction for one platform.

    BED is 0-based, half-open. Chromosome naming is adjusted to match the source VCF.
    """
    if (
        platform_plan.extract_chrom is None
        or platform_plan.extract_start is None
        or platform_plan.extract_end is None
    ):
        raise PrepareExecutionError(
            "Platform extraction coordinates are missing for {}.".format(
                platform_plan.platform
            )
        )

    bed_start = max(0, platform_plan.extract_start - 1)
    bed_end = platform_plan.extract_end

    source_has_chr = _vcf_uses_chr_prefix(platform_plan.source_path)
    bed_chrom = _format_chrom_for_source_vcf(
        chrom=platform_plan.extract_chrom,
        add_chr=source_has_chr,
    )

    with bed_path.open("w", encoding="utf-8") as handle:
        handle.write(
            "{}\t{}\t{}\t{}\n".format(
                bed_chrom,
                bed_start,
                bed_end,
                platform_plan.platform,
            )
        )


def _extract_platform_region(
    platform_plan: PlatformInputPlan,
    sample_ids_file: Path,
    region_bed: Path,
    output_vcf: Path,
) -> None:
    cmd = [
        "bcftools",
        "view",
        "--force-samples",
        "-S",
        str(sample_ids_file),
        "-R",
        str(region_bed),
        "-Oz",
        "-o",
        str(output_vcf),
        str(platform_plan.source_path),
    ]
    _run_command(cmd)


def _rename_samples(
    input_vcf: Path,
    rename_map_file: Path,
    output_vcf: Path,
) -> None:
    cmd = [
        "bcftools",
        "reheader",
        "-s",
        str(rename_map_file),
        "-o",
        str(output_vcf),
        str(input_vcf),
    ]
    _run_command(cmd)


def _liftover_vcf_with_ucsc(
    input_vcf: Path,
    output_vcf: Path,
    chain_file: Optional[Path],
) -> None:
    """
    Liftover VCF coordinates using UCSC liftOver, following the study pipeline style:

      VCF -> BED with variant tag -> liftOver -> rewrite VCF coordinates -> bgzip
    """
    if chain_file is None:
        raise PrepareExecutionError("liftover.chain_file is required for UCSC liftOver.")

    lift_over_exe = shutil.which("liftOver")
    if lift_over_exe is None:
        raise PrepareExecutionError(
            "UCSC liftOver executable was not found in PATH."
        )

    work_dir = output_vcf.parent
    plain_input_vcf = work_dir / output_vcf.name.replace(".vcf.gz", ".prelift.vcf")
    plain_output_vcf = work_dir / output_vcf.name.replace(".vcf.gz", ".postlift.vcf")
    bed_path = work_dir / output_vcf.name.replace(".vcf.gz", ".variants.bed")
    lifted_bed = work_dir / output_vcf.name.replace(".vcf.gz", ".lifted.bed")
    unmapped_bed = work_dir / output_vcf.name.replace(".vcf.gz", ".unmapped.bed")

    _decompress_vcf_to_plain(input_vcf, plain_input_vcf)
    _write_variant_bed_from_vcf(plain_input_vcf, bed_path)

    cmd = [
        lift_over_exe,
        str(bed_path),
        str(chain_file),
        str(lifted_bed),
        str(unmapped_bed),
    ]
    _run_command(cmd)

    coord_map = _read_lifted_bed_map(lifted_bed)
    _rewrite_vcf_with_lifted_coords(
        input_vcf=plain_input_vcf,
        output_vcf=plain_output_vcf,
        coord_map=coord_map,
    )
    _bgzip_file(plain_output_vcf, output_vcf)

    for tmp_path in [plain_input_vcf, plain_output_vcf, bed_path, lifted_bed]:
        if tmp_path.exists():
            tmp_path.unlink()

    if unmapped_bed.exists():
        target_unmapped = Path(str(output_vcf) + ".unmapped.bed")
        if target_unmapped.exists():
            target_unmapped.unlink()
        shutil.move(str(unmapped_bed), str(target_unmapped))


def _decompress_vcf_to_plain(input_vcf: Path, output_vcf: Path) -> None:
    with gzip.open(input_vcf, "rt", encoding="utf-8") as f_in, output_vcf.open(
        "w",
        encoding="utf-8",
    ) as f_out:
        shutil.copyfileobj(f_in, f_out)


def _write_variant_bed_from_vcf(input_vcf: Path, bed_path: Path) -> None:
    with input_vcf.open("r", encoding="utf-8") as handle, bed_path.open(
        "w",
        encoding="utf-8",
    ) as out:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue

            chrom = _normalize_chrom(parts[0])
            ucsc_chrom = _to_ucsc_chrom(chrom)
            pos = int(parts[1])
            ref = parts[3]
            alt = parts[4]
            variant_id = "{}:{}:{}:{}".format(chrom, pos, ref, alt)

            out.write(
                "{}\t{}\t{}\t{}\n".format(
                    ucsc_chrom,
                    max(0, pos - 1),
                    pos,
                    variant_id,
                )
            )


def _read_lifted_bed_map(lifted_bed: Path) -> Dict[str, Tuple[str, str]]:
    coord_map: Dict[str, Tuple[str, str]] = {}
    if not lifted_bed.exists():
        return coord_map

    with lifted_bed.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            new_chrom = _normalize_chrom(parts[0])
            new_pos = parts[2]
            variant_id = parts[3]
            coord_map[variant_id] = (new_chrom, new_pos)

    return coord_map


def _rewrite_vcf_with_lifted_coords(
    input_vcf: Path,
    output_vcf: Path,
    coord_map: Dict[str, Tuple[str, str]],
) -> None:
    """
    Rewrite VCF coordinates after liftOver.

    Variants that fail to lift are dropped, which is safer than leaving them in the
    old coordinate system.
    """
    kept = 0

    with input_vcf.open("r", encoding="utf-8") as f_in, output_vcf.open(
        "w",
        encoding="utf-8",
    ) as f_out:
        for line in f_in:
            if line.startswith("#"):
                f_out.write(line)
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue

            chrom = _normalize_chrom(parts[0])
            pos = parts[1]
            ref = parts[3]
            alt = parts[4]
            key = "{}:{}:{}:{}".format(chrom, pos, ref, alt)

            lifted = coord_map.get(key)
            if lifted is None:
                continue

            new_chrom, new_pos = lifted
            parts[0] = new_chrom
            parts[1] = new_pos
            f_out.write("\t".join(parts) + "\n")
            kept += 1

    if kept == 0:
        raise PrepareExecutionError(
            "All variants failed liftover while rewriting VCF coordinates."
        )


def _normalize_vcf(
    input_vcf: Path,
    output_vcf: Path,
    reference_fasta: Path,
) -> None:
    current_input = input_vcf

    ref_has_chr = _reference_uses_chr_prefix(reference_fasta)
    vcf_has_chr = _vcf_uses_chr_prefix(input_vcf)

    if ref_has_chr != vcf_has_chr:
        rename_map = output_vcf.parent / output_vcf.name.replace(
            ".norm.vcf.gz",
            ".rename_chroms.txt",
        )
        renamed_vcf = output_vcf.parent / output_vcf.name.replace(
            ".norm.vcf.gz",
            ".prechr.vcf.gz",
        )

        _write_chr_rename_map(rename_map, add_chr=ref_has_chr)
        _rename_vcf_contigs(input_vcf, rename_map, renamed_vcf)
        _index_vcf(renamed_vcf)
        current_input = renamed_vcf

    cmd = [
        "bcftools",
        "norm",
        "-f",
        str(reference_fasta),
        "-m",
        "-both",
        "-Oz",
        "-o",
        str(output_vcf),
        str(current_input),
    ]
    _run_command(cmd)


def _reference_uses_chr_prefix(reference_fasta: Path) -> bool:
    fai_path = Path(str(reference_fasta) + ".fai")
    if fai_path.exists():
        with fai_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                if not line.strip():
                    continue
                contig = line.split("\t", 1)[0]
                return contig.startswith("chr")

    with reference_fasta.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                contig = line[1:].strip().split()[0]
                return contig.startswith("chr")

    raise PrepareExecutionError(
        "Could not determine chromosome naming style from reference FASTA: {}".format(
            reference_fasta
        )
    )


def _vcf_uses_chr_prefix(vcf_path: Path) -> bool:
    with gzip.open(vcf_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##contig=<ID="):
                contig = line.split("##contig=<ID=", 1)[1].split(",", 1)[0].rstrip(">")
                return contig.startswith("chr")
            if line.startswith("#CHROM"):
                continue
            if line.startswith("#"):
                continue
            chrom = line.split("\t", 1)[0]
            return chrom.startswith("chr")

    raise PrepareExecutionError(
        "Could not determine chromosome naming style from VCF: {}".format(vcf_path)
    )


def _write_chr_rename_map(rename_map_path: Path, add_chr: bool) -> None:
    with rename_map_path.open("w", encoding="utf-8") as handle:
        for chrom in [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]:
            if add_chr:
                new_name = "chrM" if chrom == "MT" else "chr{}".format(chrom)
                handle.write("{}\t{}\n".format(chrom, new_name))
            else:
                old_name = "chrM" if chrom == "MT" else "chr{}".format(chrom)
                new_name = "MT" if chrom == "MT" else chrom
                handle.write("{}\t{}\n".format(old_name, new_name))


def _rename_vcf_contigs(
    input_vcf: Path,
    rename_map_file: Path,
    output_vcf: Path,
) -> None:
    cmd = [
        "bcftools",
        "annotate",
        "--rename-chrs",
        str(rename_map_file),
        "-Oz",
        "-o",
        str(output_vcf),
        str(input_vcf),
    ]
    _run_command(cmd)


def _integrate_gene_vcfs(
    wes_vcf: Path,
    array_vcf: Path,
    output_vcf: Path,
    temp_dir: Path,
) -> None:
    """
    Integrate normalized gene-level VCFs.

    Rule:
      - keep all WES records
      - append only array records whose normalized (CHROM, POS, REF, ALT)
        keys are not already present in WES
    """
    wes_samples = _read_vcf_samples(wes_vcf)
    array_samples = _read_vcf_samples(array_vcf)
    if wes_samples != array_samples:
        raise PrepareExecutionError(
            "Sample names/order mismatch between WES and array VCFs for integration."
        )

    unsorted_vcf = temp_dir / "{}.unsorted.vcf".format(output_vcf.stem.replace(".vcf", ""))

    wes_keys = set()
    with gzip.open(wes_vcf, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            wes_keys.add((fields[0], fields[1], fields[3], fields[4]))

    with open(unsorted_vcf, "w", encoding="utf-8") as out_handle:
        with gzip.open(wes_vcf, "rt", encoding="utf-8") as handle:
            for line in handle:
                out_handle.write(line)

        with gzip.open(array_vcf, "rt", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 5:
                    continue
                key = (fields[0], fields[1], fields[3], fields[4])
                if key in wes_keys:
                    continue
                out_handle.write(line)

    cmd = [
        "bcftools",
        "sort",
        "-Oz",
        "-o",
        str(output_vcf),
        str(unsorted_vcf),
    ]
    _run_command(cmd)

    if unsorted_vcf.exists():
        unsorted_vcf.unlink()


# ---------------------------------------------------------------------------
# Sample helpers
# ---------------------------------------------------------------------------


def _write_sample_files(
    samples: List[SampleRecord],
    platform: str,
    output_dir: Path,
) -> Tuple[Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)

    sample_ids_file = output_dir / "{}.samples.txt".format(platform)
    rename_map_file = output_dir / "{}.rename.txt".format(platform)

    with sample_ids_file.open("w", encoding="utf-8") as sample_handle, rename_map_file.open(
        "w",
        encoding="utf-8",
    ) as rename_handle:
        count = 0
        for record in samples:
            source_name = _sample_source_name(record, platform)
            if source_name is None:
                continue
            sample_handle.write(source_name + "\n")
            rename_handle.write("{}\t{}\n".format(source_name, record.sample_id))
            count += 1

    if count == 0:
        raise PrepareExecutionError(
            "No samples could be resolved for platform {} from sample manifest.".format(platform)
        )

    return sample_ids_file, rename_map_file


def _sample_source_name(record: SampleRecord, platform: str) -> Optional[str]:
    if platform == "wes":
        return record.wes_id or record.sample_id
    if platform == "array":
        return record.array_id or record.sample_id
    return record.sample_id


# ---------------------------------------------------------------------------
# Generic command helpers
# ---------------------------------------------------------------------------


def _copy_with_index(source_vcf: Path, dest_vcf: Path) -> None:
    shutil.copy2(str(source_vcf), str(dest_vcf))
    source_index = _find_index_file(source_vcf)
    if source_index is not None:
        dest_index = Path(str(dest_vcf) + source_index.suffix)
        shutil.copy2(str(source_index), str(dest_index))


def _find_index_file(vcf_path: Path) -> Optional[Path]:
    csi = Path(str(vcf_path) + ".csi")
    if csi.exists():
        return csi
    tbi = Path(str(vcf_path) + ".tbi")
    if tbi.exists():
        return tbi
    return None


def _index_vcf(vcf_path: Path) -> None:
    cmd = ["bcftools", "index", "-f", str(vcf_path)]
    _run_command(cmd)


def _bgzip_file(input_path: Path, output_path: Path) -> None:
    bgzip_exe = shutil.which("bgzip")
    if bgzip_exe is None:
        raise PrepareExecutionError("bgzip executable was not found in PATH.")

    with output_path.open("wb") as out_handle:
        cmd = [bgzip_exe, "-c", str(input_path)]
        result = subprocess.run(cmd, stdout=out_handle, stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            raise PrepareExecutionError(
                "Command failed: {}\n{}".format(" ".join(cmd), result.stderr)
            )


def _read_vcf_samples(vcf_path: Path) -> List[str]:
    cmd = ["bcftools", "query", "-l", str(vcf_path)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise PrepareExecutionError(
            "Failed to read VCF samples from {}:\n{}".format(vcf_path, result.stderr)
        )
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def _run_command(cmd: List[str]) -> None:
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise PrepareExecutionError(
            "Command failed: {}\nSTDOUT:\n{}\nSTDERR:\n{}".format(
                " ".join(cmd),
                result.stdout,
                result.stderr,
            )
        )


def _cleanup_gene_work_dir(gene_work_dir: Path) -> None:
    if not gene_work_dir.exists():
        return

    for path in gene_work_dir.iterdir():
        if path.is_file():
            path.unlink()
        elif path.is_dir():
            shutil.rmtree(path)

    try:
        gene_work_dir.rmdir()
    except OSError:
        pass


# ---------------------------------------------------------------------------
# Internal plan helpers
# ---------------------------------------------------------------------------


def _resolve_registry_tsv(
    cfg: AppConfig,
    registry_tsv: Optional[Union[str, Path]] = None,
) -> Path:
    if registry_tsv is not None:
        return Path(registry_tsv).expanduser().resolve()

    if cfg.config_path is None:
        raise PreparePlanError("Config path is missing; cannot derive default registry TSV path.")

    default_path = cfg.config_path.parent / "resources" / "gene_registry_{}.tsv".format(
        normalize_build(cfg.target_build)
    )
    return default_path.resolve()


def _rows_to_gene_lookup(rows: List[Dict[str, str]]) -> Dict[str, Dict[str, str]]:
    lookup: Dict[str, Dict[str, str]] = {}
    for row in rows:
        gene = str(row["gene"]).upper()
        lookup[gene] = row
    return lookup


def _platforms_in_order(datasets: List[DatasetRecord]) -> List[str]:
    preferred = ["wes", "array"]
    seen = set(record.platform for record in datasets)
    ordered = [platform for platform in preferred if platform in seen]
    extras = sorted(platform for platform in seen if platform not in preferred)
    return ordered + extras


def _select_dataset_for_gene(
    platform: str,
    gene: str,
    chrom: str,
    gene_datasets: Dict[str, Dict[str, DatasetRecord]],
    chrom_datasets: Dict[str, Dict[str, DatasetRecord]],
) -> Optional[DatasetRecord]:
    gene_bucket = gene_datasets.get(platform, {})
    if gene in gene_bucket:
        return gene_bucket[gene]

    chrom_bucket = chrom_datasets.get(platform, {})
    if chrom in chrom_bucket:
        return chrom_bucket[chrom]

    return None


def _normalize_chrom(chrom: str) -> str:
    value = chrom.strip()
    if not value:
        raise PreparePlanError("Empty chromosome value encountered while building plan.")
    if value.lower().startswith("chr"):
        value = value[3:]
    value = value.upper()
    if value == "M":
        value = "MT"
    return value


def _to_ucsc_chrom(chrom: str) -> str:
    value = _normalize_chrom(chrom)
    if value == "MT":
        return "chrM"
    return "chr{}".format(value)


def _format_chrom_for_source_vcf(chrom: str, add_chr: bool) -> str:
    value = _normalize_chrom(chrom)
    if add_chr:
        if value == "MT":
            return "chrM"
        return "chr{}".format(value)
    return value