from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import yaml

from .aldy import (
    format_aldy_plan,
    run_aldy,
    run_aldy_dry_run,
    run_aldy_simple,
    run_aldy_simple_dry_run,
)
from .gene_registry import normalize_build
from .prepare import format_prepare_plan, run_prepare, run_prepare_dry_run
from .split import (
    format_split_plan,
    run_split,
    run_split_dry_run,
    run_split_simple,
    run_split_simple_dry_run,
)


def main(argv: Optional[List[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if not getattr(args, "command", None):
        parser.print_help()
        return 1

    try:
        if args.command == "prepare":
            return _handle_prepare(args)
        if args.command == "split":
            return _handle_split(args)
        if args.command == "aldy":
            return _handle_aldy(args)

        parser.error("Unknown command: {}".format(args.command))
        return 1

    except Exception as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 2


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pgx",
        description="PGxIntegrator command-line interface",
    )
    subparsers = parser.add_subparsers(dest="command")

    # ------------------------------------------------------------------
    # prepare
    # ------------------------------------------------------------------
    prepare_parser = subparsers.add_parser(
        "prepare",
        help="Prepare per-gene integrated VCFs",
        description=(
            "Prepare supports two modes:\n"
            "  1) Config mode:  pgx prepare --config config.yaml\n"
            "  2) Simple mode:  pgx prepare --wes <file_or_dir> --array <file_or_dir> --genes ...\n\n"
            "Simple mode accepts either a single .vcf.gz file or a directory.\n"
            "Directory mode only works when filenames clearly indicate chromosome\n"
            "(e.g. chr10.vcf.gz) or requested gene (e.g. CYP3A5.vcf.gz).\n"
            "Ambiguous cases should use --config."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    prepare_parser.add_argument("--config", default=None, help="Path to YAML config file (advanced mode)")
    prepare_parser.add_argument("--dry-run", action="store_true", help="Print the prepare plan without executing commands")
    prepare_parser.add_argument("--registry-tsv", default=None, help="Optional path to gene registry TSV")

    prepare_parser.add_argument("--project-name", default="pgx_prepare_run", help="Project name for internally generated config")
    prepare_parser.add_argument("--reference-fasta", default=None, help="Reference FASTA path (required in simple mode)")
    prepare_parser.add_argument("--reference-build", default="hg38", help="Target/reference build (default: hg38)")
    prepare_parser.add_argument("--wes", default=None, help="WES input path: a .vcf.gz file or a directory containing .vcf.gz files")
    prepare_parser.add_argument("--wes-build", default=None, help="Build of WES input(s) (default: same as --reference-build)")
    prepare_parser.add_argument("--array", default=None, help="Array input path: a .vcf.gz file or a directory containing .vcf.gz files")
    prepare_parser.add_argument("--array-build", default=None, help="Build of array input(s) (default: hg19)")
    prepare_parser.add_argument(
        "--sample-manifest",
        default=None,
        help=(
            "Optional sample manifest TSV. If omitted, an identity mapping is generated "
            "from VCF sample names. If WES/array sample IDs differ, this file is required."
        ),
    )
    prepare_parser.add_argument("--genes", action="append", default=None, help="Comma-separated genes, e.g. --genes CYP2C19,CYP3A5")
    prepare_parser.add_argument("--output-dir", default="results/prepared", help="Prepare output directory (default: results/prepared)")
    prepare_parser.add_argument("--split-output-dir", default=None, help="Internal split output directory to write into temporary config")
    prepare_parser.add_argument("--aldy-output-dir", default=None, help="Internal Aldy output directory to write into temporary config")
    prepare_parser.add_argument("--keep-intermediate", action="store_true", help="Keep prepare intermediate files")
    prepare_parser.add_argument("--no-resume", action="store_true", help="Disable resume for prepare")

    # ------------------------------------------------------------------
    # split
    # ------------------------------------------------------------------
    split_parser = subparsers.add_parser(
        "split",
        help="Split prepared per-gene VCFs into single-sample VCFs",
        description=(
            "Split supports two modes:\n"
            "  1) Config mode:  pgx split --config config.yaml\n"
            "  2) Simple mode:  pgx split --input results/prepared --output-dir results/single_sample_vcfs\n\n"
            "Simple mode reads canonical sample IDs directly from prepared VCFs."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    split_parser.add_argument("--config", default=None, help="Path to YAML config file (advanced mode)")
    split_parser.add_argument("--dry-run", action="store_true", help="Print the split plan without executing commands")
    split_parser.add_argument("--gene", action="append", default=None, help="Restrict split to one gene. Repeat for multiple genes.")

    split_parser.add_argument("--project-name", default="pgx_split_run", help="Project name in simple mode")
    split_parser.add_argument("--input", "--input-dir", dest="input_path", default=None, help="Prepared input path: a per-gene VCF or a directory of per-gene VCFs")
    split_parser.add_argument("--genes", action="append", default=None, help="Comma-separated genes, e.g. --genes CYP2C19,CYP3A5")
    split_parser.add_argument("--output-dir", default="results/single_sample_vcfs", help="Split output directory (default: results/single_sample_vcfs)")
    split_parser.add_argument("--samples-per-batch", type=int, default=1000, help="Samples per batch directory (default: 1000)")
    split_parser.add_argument("--jobs", type=int, default=1, help="Number of parallel split jobs (default: 1)")
    split_parser.add_argument("--no-resume", action="store_true", help="Disable resume for split")

    # ------------------------------------------------------------------
    # aldy
    # ------------------------------------------------------------------
    aldy_parser = subparsers.add_parser(
        "aldy",
        help="Run Aldy on split per-sample VCFs",
        description=(
            "Aldy supports two modes:\n"
            "  1) Config mode:  pgx aldy --config config.yaml\n"
            "  2) Simple mode:  pgx aldy --input results/single_sample_vcfs --output-dir results/aldy_results --reference-fasta ref.fa\n\n"
            "Simple mode is useful as a parallel Aldy runner."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    aldy_parser.add_argument("--config", default=None, help="Path to YAML config file (advanced mode)")
    aldy_parser.add_argument("--dry-run", action="store_true", help="Print the Aldy plan without executing commands")
    aldy_parser.add_argument("--gene", action="append", default=None, help="Restrict Aldy to one gene. Repeat for multiple genes.")

    aldy_parser.add_argument("--project-name", default="pgx_aldy_run", help="Project name in simple mode")
    aldy_parser.add_argument("--input", "--input-dir", dest="input_path", default=None, help="Split input path: a split sample VCF or a directory of split inputs")
    aldy_parser.add_argument("--genes", action="append", default=None, help="Comma-separated genes, e.g. --genes CYP2C19,CYP3A5")
    aldy_parser.add_argument("--output-dir", default="results/aldy_results", help="Aldy output directory (default: results/aldy_results)")
    aldy_parser.add_argument("--reference-fasta", default=None, help="Reference FASTA path (required in simple mode)")
    aldy_parser.add_argument("--genome", default="hg38", help="Aldy genome option (default: hg38)")
    aldy_parser.add_argument("--profile", default="wgs", help="Aldy profile option (default: wgs)")
    aldy_parser.add_argument("--jobs", type=int, default=1, help="Number of parallel Aldy jobs (default: 1)")
    aldy_parser.add_argument("--no-resume", action="store_true", help="Disable resume for Aldy")

    return parser


# ---------------------------------------------------------------------------
# Handlers
# ---------------------------------------------------------------------------


def _handle_prepare(args: argparse.Namespace) -> int:
    if args.config:
        config_path = Path(args.config)
        registry_tsv = Path(args.registry_tsv) if args.registry_tsv else None

        if args.dry_run:
            plan = run_prepare_dry_run(config_path=config_path, registry_tsv=registry_tsv)
            print(format_prepare_plan(plan))
            return 0

        plan = run_prepare(config_path=config_path, registry_tsv=registry_tsv, dry_run=False)
        print(format_prepare_plan(plan))
        print("Prepare finished.")
        return 0

    return _handle_prepare_simple_mode(args)


def _handle_split(args: argparse.Namespace) -> int:
    if args.config:
        config_path = Path(args.config)
        genes = _collect_gene_args(repeated=args.gene, comma_genes=None)

        if args.dry_run:
            plan = run_split_dry_run(config_path=config_path, genes=genes)
            print(format_split_plan(plan))
            return 0

        plan = run_split(config_path=config_path, genes=genes, dry_run=False)
        print(format_split_plan(plan))
        print("Split finished.")
        return 0

    return _handle_split_simple_mode(args)


def _handle_aldy(args: argparse.Namespace) -> int:
    if args.config:
        config_path = Path(args.config)
        genes = _collect_gene_args(repeated=args.gene, comma_genes=None)

        if args.dry_run:
            plan = run_aldy_dry_run(config_path=config_path, genes=genes)
            print(format_aldy_plan(plan))
            return 0

        plan = run_aldy(config_path=config_path, genes=genes, dry_run=False)
        print(format_aldy_plan(plan))
        print("Aldy finished.")
        return 0

    return _handle_aldy_simple_mode(args)


# ---------------------------------------------------------------------------
# Prepare simple mode
# ---------------------------------------------------------------------------


def _handle_prepare_simple_mode(args: argparse.Namespace) -> int:
    genes = _parse_genes_argument(args.genes)
    if not genes:
        raise ValueError("Simple mode requires --genes, e.g. --genes CYP2C19,CYP3A5")

    if not args.reference_fasta:
        raise ValueError("Simple mode requires --reference-fasta")

    if not args.wes and not args.array:
        raise ValueError("Simple mode requires at least one of --wes or --array")

    reference_fasta = Path(args.reference_fasta).expanduser().resolve()
    reference_build = normalize_build(args.reference_build)

    if not reference_fasta.exists():
        raise ValueError("Reference FASTA not found: {}".format(reference_fasta))

    wes_input = Path(args.wes).expanduser().resolve() if args.wes else None
    array_input = Path(args.array).expanduser().resolve() if args.array else None
    wes_build = normalize_build(args.wes_build) if args.wes_build else reference_build
    array_build = normalize_build(args.array_build) if args.array_build else "hg19"

    output_dir = Path(args.output_dir).expanduser().resolve()
    split_output_dir = (
        Path(args.split_output_dir).expanduser().resolve()
        if args.split_output_dir
        else output_dir.parent / "single_sample_vcfs"
    )
    aldy_output_dir = (
        Path(args.aldy_output_dir).expanduser().resolve()
        if args.aldy_output_dir
        else output_dir.parent / "aldy_results"
    )

    registry_tsv = (
        Path(args.registry_tsv).expanduser().resolve()
        if args.registry_tsv
        else (output_dir.parent / "resources" / "gene_registry_{}.tsv".format(reference_build))
    )

    with tempfile.TemporaryDirectory(prefix="pgx_prepare_") as tmp:
        tmpdir = Path(tmp)
        input_dir = tmpdir / "inputs"
        input_dir.mkdir(parents=True, exist_ok=True)

        dataset_rows, wes_first_vcf, array_first_vcf = _build_dataset_rows(
            wes_input=wes_input,
            wes_build=wes_build,
            array_input=array_input,
            array_build=array_build,
            requested_genes=genes,
        )

        sample_manifest_path = _resolve_or_build_sample_manifest(
            sample_manifest_arg=args.sample_manifest,
            tmpdir=input_dir,
            wes_first_vcf=wes_first_vcf,
            array_first_vcf=array_first_vcf,
        )

        gene_list_path = _write_gene_list(
            genes=genes,
            output_path=input_dir / "gene_list.tsv",
        )

        dataset_manifest_path = _write_dataset_manifest(
            rows=dataset_rows,
            output_path=input_dir / "dataset_manifest.tsv",
        )

        config_path = _write_temp_prepare_config(
            output_path=tmpdir / "config.yaml",
            project_name=args.project_name,
            reference_fasta=reference_fasta,
            reference_build=reference_build,
            sample_manifest=sample_manifest_path,
            dataset_manifest=dataset_manifest_path,
            gene_list=gene_list_path,
            prepare_output_dir=output_dir,
            split_output_dir=split_output_dir,
            aldy_output_dir=aldy_output_dir,
            keep_intermediate=bool(args.keep_intermediate),
            resume=(not args.no_resume),
            liftover_enabled=_should_enable_liftover(
                reference_build=reference_build,
                wes_build=wes_build if wes_input else None,
                array_build=array_build if array_input else None,
            ),
        )

        if args.dry_run:
            plan = run_prepare_dry_run(config_path=config_path, registry_tsv=registry_tsv)
            print(format_prepare_plan(plan))
            return 0

        plan = run_prepare(config_path=config_path, registry_tsv=registry_tsv, dry_run=False)
        print(format_prepare_plan(plan))
        print("Prepare finished.")
        return 0


# ---------------------------------------------------------------------------
# Split simple mode
# ---------------------------------------------------------------------------


def _handle_split_simple_mode(args: argparse.Namespace) -> int:
    if not args.input_path:
        raise ValueError("Simple split mode requires --input or --input-dir")

    genes = _collect_gene_args(repeated=args.gene, comma_genes=args.genes)
    input_path = Path(args.input_path).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()

    if args.dry_run:
        plan = run_split_simple_dry_run(
            input_path=input_path,
            output_dir=output_dir,
            genes=genes,
            samples_per_batch=args.samples_per_batch,
            project_name=args.project_name,
        )
        print(format_split_plan(plan))
        return 0

    plan = run_split_simple(
        input_path=input_path,
        output_dir=output_dir,
        genes=genes,
        samples_per_batch=args.samples_per_batch,
        jobs=args.jobs,
        resume=(not args.no_resume),
        project_name=args.project_name,
        dry_run=False,
    )
    print(format_split_plan(plan))
    print("Split finished.")
    return 0


# ---------------------------------------------------------------------------
# Aldy simple mode
# ---------------------------------------------------------------------------


def _handle_aldy_simple_mode(args: argparse.Namespace) -> int:
    if not args.input_path:
        raise ValueError("Simple Aldy mode requires --input or --input-dir")

    if not args.reference_fasta:
        raise ValueError("Simple Aldy mode requires --reference-fasta")

    genes = _collect_gene_args(repeated=args.gene, comma_genes=args.genes)
    input_path = Path(args.input_path).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    reference_fasta = Path(args.reference_fasta).expanduser().resolve()

    if args.dry_run:
        plan = run_aldy_simple_dry_run(
            input_path=input_path,
            output_dir=output_dir,
            reference_fasta=reference_fasta,
            genome=args.genome,
            profile=args.profile,
            genes=genes,
            project_name=args.project_name,
        )
        print(format_aldy_plan(plan))
        return 0

    plan = run_aldy_simple(
        input_path=input_path,
        output_dir=output_dir,
        reference_fasta=reference_fasta,
        genome=args.genome,
        profile=args.profile,
        genes=genes,
        jobs=args.jobs,
        resume=(not args.no_resume),
        project_name=args.project_name,
        dry_run=False,
    )
    print(format_aldy_plan(plan))
    print("Aldy finished.")
    return 0


# ---------------------------------------------------------------------------
# Shared temp-config builders for prepare
# ---------------------------------------------------------------------------


def _build_dataset_rows(
    wes_input: Optional[Path],
    wes_build: str,
    array_input: Optional[Path],
    array_build: str,
    requested_genes: List[str],
) -> Tuple[List[Dict[str, str]], Optional[Path], Optional[Path]]:
    rows: List[Dict[str, str]] = []
    wes_first_vcf = None
    array_first_vcf = None

    if wes_input is not None:
        wes_rows, wes_first_vcf = _scan_dataset_input(
            dataset_input=wes_input,
            platform="wes",
            build=wes_build,
            requested_genes=requested_genes,
        )
        rows.extend(wes_rows)

    if array_input is not None:
        array_rows, array_first_vcf = _scan_dataset_input(
            dataset_input=array_input,
            platform="array",
            build=array_build,
            requested_genes=requested_genes,
        )
        rows.extend(array_rows)

    if not rows:
        raise ValueError("No input VCFs were found from the supplied inputs.")

    return rows, wes_first_vcf, array_first_vcf


def _scan_dataset_input(
    dataset_input: Path,
    platform: str,
    build: str,
    requested_genes: List[str],
) -> Tuple[List[Dict[str, str]], Path]:
    if not dataset_input.exists():
        raise ValueError("Input path not found: {}".format(dataset_input))

    if dataset_input.is_file():
        if not dataset_input.name.endswith(".vcf.gz"):
            raise ValueError("Input file must end with .vcf.gz: {}".format(dataset_input))
        vcfs = [dataset_input.resolve()]
    elif dataset_input.is_dir():
        vcfs = sorted(
            path.resolve()
            for path in dataset_input.glob("*.vcf.gz")
            if not path.name.endswith(".vcf.gz.tbi")
            and not path.name.endswith(".vcf.gz.csi")
        )
    else:
        raise ValueError("Unsupported input path: {}".format(dataset_input))

    if not vcfs:
        raise ValueError(
            "No .vcf.gz files found in {} for platform {}".format(dataset_input, platform)
        )

    rows: List[Dict[str, str]] = []
    seen_chroms = set()
    seen_genes = set()

    for path in vcfs:
        stem = path.name[:-7]
        gene = _infer_gene_from_name(stem, requested_genes)
        chrom = _infer_chrom_from_name(stem)

        if gene is not None:
            layout = "gene_by_file"
            if gene in seen_genes:
                raise ValueError(
                    "Ambiguous {} input: more than one file matched gene {}. Use --config.".format(
                        platform, gene
                    )
                )
            seen_genes.add(gene)
            chrom_value = ""
            gene_value = gene
        elif chrom is not None:
            layout = "cohort_by_chrom"
            if chrom in seen_chroms:
                raise ValueError(
                    "Ambiguous {} input: more than one file matched chromosome {}. Use --config.".format(
                        platform, chrom
                    )
                )
            seen_chroms.add(chrom)
            chrom_value = chrom
            gene_value = ""
        else:
            raise ValueError(
                "Could not infer chromosome or requested gene from file name '{}'. "
                "Use clearer names or switch to --config.".format(path.name)
            )

        rows.append(
            {
                "dataset_id": platform,
                "platform": platform,
                "build": build,
                "layout": layout,
                "chrom": chrom_value,
                "gene": gene_value,
                "path": str(path),
            }
        )

    return rows, vcfs[0]


def _resolve_or_build_sample_manifest(
    sample_manifest_arg: Optional[str],
    tmpdir: Path,
    wes_first_vcf: Optional[Path],
    array_first_vcf: Optional[Path],
) -> Path:
    if sample_manifest_arg:
        path = Path(sample_manifest_arg).expanduser().resolve()
        if not path.exists():
            raise ValueError("Sample manifest not found: {}".format(path))
        return path

    wes_samples = _read_vcf_samples(wes_first_vcf) if wes_first_vcf else []
    array_samples = _read_vcf_samples(array_first_vcf) if array_first_vcf else []

    if not wes_samples and not array_samples:
        raise ValueError("Could not infer sample names because no input VCFs were available.")

    if wes_samples and array_samples:
        if set(wes_samples) != set(array_samples):
            raise ValueError(
                "WES and array sample IDs differ. Please provide --sample-manifest "
                "for explicit ID mapping."
            )
        canonical_samples = list(wes_samples)
    elif wes_samples:
        canonical_samples = list(wes_samples)
    else:
        canonical_samples = list(array_samples)

    output_path = tmpdir / "sample_manifest.tsv"
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("sample_id\tarray_id\twes_id\n")
        for sample_id in canonical_samples:
            array_id = sample_id if array_samples else ""
            wes_id = sample_id if wes_samples else ""
            handle.write("{}\t{}\t{}\n".format(sample_id, array_id, wes_id))

    return output_path.resolve()


def _write_gene_list(genes: List[str], output_path: Path) -> Path:
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("gene\n")
        for gene in genes:
            handle.write("{}\n".format(gene))
    return output_path.resolve()


def _write_dataset_manifest(rows: List[Dict[str, str]], output_path: Path) -> Path:
    rows = sorted(
        rows,
        key=lambda row: (
            row["platform"],
            row["layout"],
            row["chrom"],
            row["gene"],
            row["path"],
        ),
    )
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("dataset_id\tplatform\tbuild\tlayout\tchrom\tgene\tpath\n")
        for row in rows:
            handle.write(
                "{dataset_id}\t{platform}\t{build}\t{layout}\t{chrom}\t{gene}\t{path}\n".format(
                    **row
                )
            )
    return output_path.resolve()


def _write_temp_prepare_config(
    output_path: Path,
    project_name: str,
    reference_fasta: Path,
    reference_build: str,
    sample_manifest: Path,
    dataset_manifest: Path,
    gene_list: Path,
    prepare_output_dir: Path,
    split_output_dir: Path,
    aldy_output_dir: Path,
    keep_intermediate: bool,
    resume: bool,
    liftover_enabled: bool,
) -> Path:
    cfg = {
        "project_name": project_name,
        "reference": {
            "fasta": str(reference_fasta.resolve()),
            "build": reference_build,
        },
        "inputs": {
            "sample_manifest": str(sample_manifest.resolve()),
            "dataset_manifest": str(dataset_manifest.resolve()),
            "gene_list": str(gene_list.resolve()),
        },
        "prepare": {
            "output_dir": str(prepare_output_dir.resolve()),
            "normalize": True,
            "harmonize": True,
            "sample_missingness_threshold": 0.05,
            "keep_intermediate": bool(keep_intermediate),
            "resume": bool(resume),
        },
        "liftover": {
            "enabled": bool(liftover_enabled),
            "target_build": reference_build,
        },
        "split": {
            "output_dir": str(split_output_dir.resolve()),
            "input_dir": str(prepare_output_dir.resolve()),
            "samples_per_batch": 1000,
            "jobs": 1,
            "resume": True,
        },
        "aldy": {
            "output_dir": str(aldy_output_dir.resolve()),
            "input_dir": str(split_output_dir.resolve()),
            "jobs": 1,
            "genome": reference_build,
            "profile": "wgs",
            "reference_fasta": str(reference_fasta.resolve()),
            "resume": True,
        },
    }

    with output_path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(cfg, handle, sort_keys=False)

    return output_path.resolve()


def _should_enable_liftover(
    reference_build: str,
    wes_build: Optional[str],
    array_build: Optional[str],
) -> bool:
    builds = []
    if wes_build is not None:
        builds.append(normalize_build(wes_build))
    if array_build is not None:
        builds.append(normalize_build(array_build))
    return any(build != reference_build for build in builds)


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------


def _parse_genes_argument(values: Optional[List[str]]) -> List[str]:
    if not values:
        return []

    genes: List[str] = []
    seen = set()

    for chunk in values:
        for item in chunk.split(","):
            gene = item.strip().upper()
            if not gene:
                continue
            if gene not in seen:
                genes.append(gene)
                seen.add(gene)

    return genes


def _collect_gene_args(
    repeated: Optional[List[str]],
    comma_genes: Optional[List[str]],
) -> Optional[List[str]]:
    merged: List[str] = []
    seen = set()

    for gene in repeated or []:
        value = gene.strip().upper()
        if value and value not in seen:
            merged.append(value)
            seen.add(value)

    for gene in _parse_genes_argument(comma_genes):
        if gene not in seen:
            merged.append(gene)
            seen.add(gene)

    return merged or None


def _infer_gene_from_name(stem: str, requested_genes: List[str]) -> Optional[str]:
    token = stem.upper()
    for gene in requested_genes:
        pattern = r"(^|[^A-Z0-9]){}([^A-Z0-9]|$)".format(re.escape(gene))
        if re.search(pattern, token):
            return gene
        if token == gene:
            return gene
    return None


def _infer_chrom_from_name(stem: str) -> Optional[str]:
    token = stem.upper()
    match = re.search(
        r"(^|[^A-Z0-9])(?:CHR)?(?P<chrom>(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT))([^A-Z0-9]|$)",
        token,
    )
    if not match:
        return None

    chrom = match.group("chrom")
    if chrom == "M":
        chrom = "MT"
    return chrom


def _read_vcf_samples(vcf_path: Optional[Path]) -> List[str]:
    if vcf_path is None:
        return []

    cmd = ["bcftools", "query", "-l", str(vcf_path)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise ValueError(
            "Failed to read VCF samples from {}:\n{}".format(vcf_path, result.stderr)
        )

    return [line.strip() for line in result.stdout.splitlines() if line.strip()]