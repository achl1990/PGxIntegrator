from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple, Union
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

from .config import AppConfig, load_config
from .manifests import load_gene_list, load_sample_manifest, list_gene_names


class SplitPlanError(ValueError):
    """Raised when a split plan cannot be built."""


class SplitExecutionError(RuntimeError):
    """Raised when the split stage fails during execution."""


@dataclass
class SplitSamplePlan:
    sample_id: str
    batch_name: str
    output_vcf: Path


@dataclass
class SplitGenePlan:
    gene: str
    input_vcf: Path
    sample_count: int
    batch_count: int
    sample_plans: List[SplitSamplePlan] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)


@dataclass
class SplitPlan:
    project_name: str
    input_dir: Path
    output_dir: Path
    gene_count: int
    total_sample_jobs: int
    gene_plans: List[SplitGenePlan]


# ---------------------------------------------------------------------------
# Public API: config mode
# ---------------------------------------------------------------------------


def build_split_plan(
    cfg: AppConfig,
    genes: Optional[List[str]] = None,
) -> SplitPlan:
    """
    Build a split plan from YAML config mode.
    """
    samples = load_sample_manifest(cfg.inputs.sample_manifest)
    gene_list_records = load_gene_list(cfg.inputs.gene_list)

    configured_genes = list_gene_names(gene_list_records)
    selected_genes = _select_genes(configured_genes, genes)
    input_path = _resolve_split_input_path(cfg)
    preferred_sample_order = [record.sample_id for record in samples]

    gene_inputs = _discover_prepared_gene_inputs(
        input_path=input_path,
        requested_genes=selected_genes,
    )

    return _build_split_plan_from_inputs(
        project_name=cfg.project_name,
        input_root=input_path,
        output_dir=cfg.split.output_dir,
        gene_inputs=gene_inputs,
        samples_per_batch=cfg.split.samples_per_batch,
        preferred_sample_order=preferred_sample_order,
    )


def run_split_dry_run(
    config_path: Union[str, Path],
    genes: Optional[List[str]] = None,
) -> SplitPlan:
    cfg = load_config(config_path)
    return build_split_plan(cfg, genes=genes)


def run_split(
    config_path: Union[str, Path],
    genes: Optional[List[str]] = None,
    dry_run: bool = False,
) -> SplitPlan:
    cfg = load_config(config_path)
    plan = build_split_plan(cfg, genes=genes)

    if dry_run:
        return plan

    _execute_split_jobs(
        plan=plan,
        jobs=max(1, cfg.split.jobs),
        resume=cfg.split.resume,
    )
    return plan


# ---------------------------------------------------------------------------
# Public API: simple mode
# ---------------------------------------------------------------------------


def build_split_plan_simple(
    input_path: Union[str, Path],
    output_dir: Union[str, Path],
    genes: Optional[List[str]] = None,
    samples_per_batch: int = 1000,
    project_name: str = "pgx_split_run",
) -> SplitPlan:
    """
    Build a split plan directly from prepared per-gene VCF input(s).

    Simple mode does not require:
      - YAML config
      - sample manifest
      - reference FASTA

    It infers:
      - genes from prepared VCF file names like GENE.vcf.gz
      - sample IDs directly from the prepared VCF sample columns
    """
    resolved_input = Path(input_path).expanduser().resolve()
    resolved_output = Path(output_dir).expanduser().resolve()
    selected_genes = _normalize_genes(genes)

    gene_inputs = _discover_prepared_gene_inputs(
        input_path=resolved_input,
        requested_genes=selected_genes,
    )

    return _build_split_plan_from_inputs(
        project_name=project_name,
        input_root=resolved_input,
        output_dir=resolved_output,
        gene_inputs=gene_inputs,
        samples_per_batch=samples_per_batch,
        preferred_sample_order=None,
    )


def run_split_simple_dry_run(
    input_path: Union[str, Path],
    output_dir: Union[str, Path],
    genes: Optional[List[str]] = None,
    samples_per_batch: int = 1000,
    project_name: str = "pgx_split_run",
) -> SplitPlan:
    return build_split_plan_simple(
        input_path=input_path,
        output_dir=output_dir,
        genes=genes,
        samples_per_batch=samples_per_batch,
        project_name=project_name,
    )


def run_split_simple(
    input_path: Union[str, Path],
    output_dir: Union[str, Path],
    genes: Optional[List[str]] = None,
    samples_per_batch: int = 1000,
    jobs: int = 1,
    resume: bool = True,
    project_name: str = "pgx_split_run",
    dry_run: bool = False,
) -> SplitPlan:
    plan = build_split_plan_simple(
        input_path=input_path,
        output_dir=output_dir,
        genes=genes,
        samples_per_batch=samples_per_batch,
        project_name=project_name,
    )

    if dry_run:
        return plan

    _execute_split_jobs(
        plan=plan,
        jobs=max(1, jobs),
        resume=resume,
    )
    return plan


def format_split_plan(plan: SplitPlan) -> str:
    lines: List[str] = []
    lines.append("Project: {}".format(plan.project_name))
    lines.append("Split input dir: {}".format(plan.input_dir))
    lines.append("Split output dir: {}".format(plan.output_dir))
    lines.append("Gene count: {}".format(plan.gene_count))
    lines.append("Total sample jobs: {}".format(plan.total_sample_jobs))
    lines.append("")

    for gene_plan in plan.gene_plans:
        lines.append(
            "[{}] input={} | samples={} | batches={}".format(
                gene_plan.gene,
                gene_plan.input_vcf,
                gene_plan.sample_count,
                gene_plan.batch_count,
            )
        )

        preview = gene_plan.sample_plans[:5]
        for item in preview:
            lines.append("  - {} -> {}".format(item.sample_id, item.output_vcf))

        if len(gene_plan.sample_plans) > len(preview):
            lines.append(
                "  - ... {} more sample jobs".format(
                    len(gene_plan.sample_plans) - len(preview)
                )
            )

        for warning in gene_plan.warnings:
            lines.append("  ! {}".format(warning))

        lines.append("")

    return "\n".join(lines).rstrip() + "\n"


# ---------------------------------------------------------------------------
# Shared plan builder
# ---------------------------------------------------------------------------


def _build_split_plan_from_inputs(
    project_name: str,
    input_root: Path,
    output_dir: Path,
    gene_inputs: List[Tuple[str, Path]],
    samples_per_batch: int,
    preferred_sample_order: Optional[List[str]],
) -> SplitPlan:
    if samples_per_batch < 1:
        raise SplitPlanError("samples_per_batch must be at least 1.")

    gene_plans: List[SplitGenePlan] = []
    total_sample_jobs = 0

    for gene, input_vcf in gene_inputs:
        if not input_vcf.exists():
            raise SplitPlanError(
                "Prepared VCF not found for gene {}: {}".format(gene, input_vcf)
            )

        vcf_samples = _read_vcf_samples(input_vcf)
        if preferred_sample_order is None:
            selected_samples = list(vcf_samples)
        else:
            selected_samples = [sample_id for sample_id in preferred_sample_order if sample_id in vcf_samples]

        if not selected_samples:
            raise SplitPlanError(
                "No matching samples found in prepared VCF for gene {}.".format(gene)
            )

        gene_output_dir = output_dir / gene
        sample_plans: List[SplitSamplePlan] = []

        for idx, sample_id in enumerate(selected_samples):
            batch_idx = idx // samples_per_batch
            batch_name = "batch_{:03d}".format(batch_idx)
            output_vcf = gene_output_dir / batch_name / "{}.vcf.gz".format(sample_id)
            sample_plans.append(
                SplitSamplePlan(
                    sample_id=sample_id,
                    batch_name=batch_name,
                    output_vcf=output_vcf,
                )
            )

        batch_names = sorted(set(item.batch_name for item in sample_plans))
        gene_plans.append(
            SplitGenePlan(
                gene=gene,
                input_vcf=input_vcf,
                sample_count=len(sample_plans),
                batch_count=len(batch_names),
                sample_plans=sample_plans,
            )
        )
        total_sample_jobs += len(sample_plans)

    return SplitPlan(
        project_name=project_name,
        input_dir=input_root,
        output_dir=output_dir,
        gene_count=len(gene_plans),
        total_sample_jobs=total_sample_jobs,
        gene_plans=gene_plans,
    )


# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------


def _execute_split_jobs(
    plan: SplitPlan,
    jobs: int,
    resume: bool,
) -> None:
    plan.output_dir.mkdir(parents=True, exist_ok=True)

    split_jobs = _collect_split_jobs(plan=plan, resume=resume)
    if not split_jobs:
        return

    if jobs == 1:
        for job in split_jobs:
            _run_single_split_job(job)
        return

    errors: List[str] = []
    with ThreadPoolExecutor(max_workers=jobs) as executor:
        future_to_job = {
            executor.submit(_run_single_split_job, job): job
            for job in split_jobs
        }

        for future in as_completed(future_to_job):
            job = future_to_job[future]
            try:
                future.result()
            except Exception as exc:
                errors.append(
                    "[gene={} sample={}] {}".format(
                        job["gene"],
                        job["sample_id"],
                        exc,
                    )
                )

    if errors:
        raise SplitExecutionError(
            "One or more split jobs failed:\n{}".format("\n".join(errors))
        )


def _collect_split_jobs(plan: SplitPlan, resume: bool) -> List[dict]:
    jobs: List[dict] = []

    for gene_plan in plan.gene_plans:
        for sample_plan in gene_plan.sample_plans:
            output_vcf = sample_plan.output_vcf
            output_index = Path(str(output_vcf) + ".csi")
            output_vcf.parent.mkdir(parents=True, exist_ok=True)

            if resume and output_vcf.exists() and output_index.exists():
                continue

            jobs.append(
                {
                    "gene": gene_plan.gene,
                    "input_vcf": gene_plan.input_vcf,
                    "sample_id": sample_plan.sample_id,
                    "output_vcf": output_vcf,
                }
            )

    return jobs


def _run_single_split_job(job: dict) -> None:
    input_vcf = job["input_vcf"]
    sample_id = job["sample_id"]
    output_vcf = job["output_vcf"]

    _split_single_sample(
        input_vcf=input_vcf,
        sample_id=sample_id,
        output_vcf=output_vcf,
    )
    _index_vcf(output_vcf)


# ---------------------------------------------------------------------------
# Split helpers
# ---------------------------------------------------------------------------


def _split_single_sample(input_vcf: Path, sample_id: str, output_vcf: Path) -> None:
    cmd = [
        "bcftools",
        "view",
        "-s",
        sample_id,
        "-Oz",
        "-o",
        str(output_vcf),
        str(input_vcf),
    ]
    _run_command(cmd)


# ---------------------------------------------------------------------------
# Discovery helpers
# ---------------------------------------------------------------------------


def _discover_prepared_gene_inputs(
    input_path: Path,
    requested_genes: Optional[List[str]] = None,
) -> List[Tuple[str, Path]]:
    if not input_path.exists():
        raise SplitPlanError("Input path not found: {}".format(input_path))

    requested = _normalize_genes(requested_genes)

    if input_path.is_file():
        if not input_path.name.endswith(".vcf.gz"):
            raise SplitPlanError(
                "Prepared input file must end with .vcf.gz: {}".format(input_path)
            )

        if requested and len(requested) > 1:
            raise SplitPlanError(
                "A single prepared VCF file cannot satisfy multiple genes."
            )

        gene = requested[0] if requested else _infer_gene_from_prepared_vcf(input_path)
        return [(gene, input_path.resolve())]

    if not input_path.is_dir():
        raise SplitPlanError("Input path is not a file or directory: {}".format(input_path))

    vcfs = sorted(
        path.resolve()
        for path in input_path.glob("*.vcf.gz")
        if not path.name.endswith(".vcf.gz.csi")
        and not path.name.endswith(".vcf.gz.tbi")
    )

    if not vcfs:
        raise SplitPlanError(
            "No prepared per-gene VCFs found in {}".format(input_path)
        )

    gene_to_vcf = {}
    for path in vcfs:
        gene = _infer_gene_from_prepared_vcf(path)
        if gene in gene_to_vcf:
            raise SplitPlanError(
                "Ambiguous prepared inputs: more than one VCF matched gene {}.".format(gene)
            )
        gene_to_vcf[gene] = path

    if requested:
        missing = [gene for gene in requested if gene not in gene_to_vcf]
        if missing:
            raise SplitPlanError(
                "Prepared VCF(s) not found for requested gene(s): {}".format(
                    ", ".join(missing)
                )
            )
        return [(gene, gene_to_vcf[gene]) for gene in requested]

    discovered = sorted(gene_to_vcf.keys())
    return [(gene, gene_to_vcf[gene]) for gene in discovered]


def _infer_gene_from_prepared_vcf(vcf_path: Path) -> str:
    if not vcf_path.name.endswith(".vcf.gz"):
        raise SplitPlanError(
            "Prepared VCF must end with .vcf.gz: {}".format(vcf_path)
        )
    gene = vcf_path.name[:-7].strip().upper()
    if not gene:
        raise SplitPlanError(
            "Could not infer gene name from prepared VCF file: {}".format(vcf_path)
        )
    return gene


# ---------------------------------------------------------------------------
# Generic command helpers
# ---------------------------------------------------------------------------


def _index_vcf(vcf_path: Path) -> None:
    cmd = ["bcftools", "index", "-f", str(vcf_path)]
    _run_command(cmd)


def _read_vcf_samples(vcf_path: Path) -> List[str]:
    cmd = ["bcftools", "query", "-l", str(vcf_path)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise SplitPlanError(
            "Failed to read VCF samples from {}:\n{}".format(vcf_path, result.stderr)
        )
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def _run_command(cmd: List[str]) -> None:
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise SplitExecutionError(
            "Command failed: {}\nSTDOUT:\n{}\nSTDERR:\n{}".format(
                " ".join(cmd),
                result.stdout,
                result.stderr,
            )
        )


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _resolve_split_input_path(cfg: AppConfig) -> Path:
    if cfg.split.input_dir is not None:
        return cfg.split.input_dir
    return cfg.prepare.output_dir


def _select_genes(configured_genes: List[str], genes: Optional[List[str]]) -> List[str]:
    if not genes:
        return configured_genes

    selected: List[str] = []
    configured_set = set(configured_genes)
    for gene in genes:
        normalized = gene.strip().upper()
        if not normalized:
            continue
        if normalized not in configured_set:
            raise SplitPlanError(
                "Requested gene is not present in the configured gene list: {}".format(
                    normalized
                )
            )
        if normalized not in selected:
            selected.append(normalized)

    if not selected:
        raise SplitPlanError("No valid genes were selected for splitting.")

    return selected


def _normalize_genes(genes: Optional[List[str]]) -> Optional[List[str]]:
    if not genes:
        return None

    normalized: List[str] = []
    seen = set()

    for gene in genes:
        value = gene.strip().upper()
        if not value:
            continue
        if value not in seen:
            normalized.append(value)
            seen.add(value)

    return normalized or None