from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

from .config import AppConfig, load_config
from .manifests import load_gene_list, list_gene_names


class AldyPlanError(ValueError):
    """Raised when an Aldy plan cannot be built."""


class AldyExecutionError(RuntimeError):
    """Raised when the Aldy stage fails during execution."""


@dataclass
class AldySamplePlan:
    sample_id: str
    batch_name: str
    input_vcf: Path
    output_file: Path
    log_file: Path
    command: List[str] = field(default_factory=list)


@dataclass
class AldyGenePlan:
    gene: str
    sample_count: int
    batch_count: int
    sample_plans: List[AldySamplePlan] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)


@dataclass
class AldyPlan:
    project_name: str
    input_dir: Path
    output_dir: Path
    profile: str
    genome: str
    reference_fasta: Path
    gene_count: int
    total_jobs: int
    gene_plans: List[AldyGenePlan]


# ---------------------------------------------------------------------------
# Public API: config mode
# ---------------------------------------------------------------------------


def build_aldy_plan(
    cfg: AppConfig,
    genes: Optional[List[str]] = None,
) -> AldyPlan:
    configured_genes = list_gene_names(load_gene_list(cfg.inputs.gene_list))
    selected_genes = _select_genes(configured_genes, genes)
    input_path = _resolve_aldy_input_path(cfg)

    discovered = _discover_aldy_inputs(
        input_path=input_path,
        requested_genes=selected_genes,
    )

    return _build_aldy_plan_from_inputs(
        project_name=cfg.project_name,
        input_root=input_path,
        output_dir=cfg.aldy.output_dir,
        discovered_inputs=discovered,
        reference_fasta=_resolve_aldy_reference_fasta(cfg),
        genome=cfg.aldy.genome,
        profile=cfg.aldy.profile,
    )


def run_aldy_dry_run(
    config_path: Union[str, Path],
    genes: Optional[List[str]] = None,
) -> AldyPlan:
    cfg = load_config(config_path)
    return build_aldy_plan(cfg, genes=genes)


def run_aldy(
    config_path: Union[str, Path],
    genes: Optional[List[str]] = None,
    dry_run: bool = False,
) -> AldyPlan:
    cfg = load_config(config_path)
    plan = build_aldy_plan(cfg, genes=genes)

    if dry_run:
        return plan

    _execute_aldy_jobs(
        plan=plan,
        jobs=max(1, cfg.aldy.jobs),
        resume=cfg.aldy.resume,
    )
    return plan


# ---------------------------------------------------------------------------
# Public API: simple mode
# ---------------------------------------------------------------------------


def build_aldy_plan_simple(
    input_path: Union[str, Path],
    output_dir: Union[str, Path],
    reference_fasta: Union[str, Path],
    genome: str = "hg38",
    profile: str = "wgs",
    genes: Optional[List[str]] = None,
    project_name: str = "pgx_aldy_run",
) -> AldyPlan:
    resolved_input = Path(input_path).expanduser().resolve()
    resolved_output = Path(output_dir).expanduser().resolve()
    resolved_reference = Path(reference_fasta).expanduser().resolve()

    if not resolved_reference.exists():
        raise AldyPlanError("Reference FASTA not found: {}".format(resolved_reference))

    selected_genes = _normalize_genes(genes)
    discovered = _discover_aldy_inputs(
        input_path=resolved_input,
        requested_genes=selected_genes,
    )

    return _build_aldy_plan_from_inputs(
        project_name=project_name,
        input_root=resolved_input,
        output_dir=resolved_output,
        discovered_inputs=discovered,
        reference_fasta=resolved_reference,
        genome=genome,
        profile=profile,
    )


def run_aldy_simple_dry_run(
    input_path: Union[str, Path],
    output_dir: Union[str, Path],
    reference_fasta: Union[str, Path],
    genome: str = "hg38",
    profile: str = "wgs",
    genes: Optional[List[str]] = None,
    project_name: str = "pgx_aldy_run",
) -> AldyPlan:
    return build_aldy_plan_simple(
        input_path=input_path,
        output_dir=output_dir,
        reference_fasta=reference_fasta,
        genome=genome,
        profile=profile,
        genes=genes,
        project_name=project_name,
    )


def run_aldy_simple(
    input_path: Union[str, Path],
    output_dir: Union[str, Path],
    reference_fasta: Union[str, Path],
    genome: str = "hg38",
    profile: str = "wgs",
    genes: Optional[List[str]] = None,
    jobs: int = 1,
    resume: bool = True,
    project_name: str = "pgx_aldy_run",
    dry_run: bool = False,
) -> AldyPlan:
    plan = build_aldy_plan_simple(
        input_path=input_path,
        output_dir=output_dir,
        reference_fasta=reference_fasta,
        genome=genome,
        profile=profile,
        genes=genes,
        project_name=project_name,
    )

    if dry_run:
        return plan

    _execute_aldy_jobs(
        plan=plan,
        jobs=max(1, jobs),
        resume=resume,
    )
    return plan


def format_aldy_plan(plan: AldyPlan) -> str:
    lines: List[str] = []
    lines.append("Project: {}".format(plan.project_name))
    lines.append("Aldy input dir: {}".format(plan.input_dir))
    lines.append("Aldy output dir: {}".format(plan.output_dir))
    lines.append("Profile: {}".format(plan.profile))
    lines.append("Genome: {}".format(plan.genome))
    lines.append("Reference FASTA: {}".format(plan.reference_fasta))
    lines.append("Gene count: {}".format(plan.gene_count))
    lines.append("Total jobs: {}".format(plan.total_jobs))
    lines.append("")

    for gene_plan in plan.gene_plans:
        lines.append(
            "[{}] samples={} | batches={}".format(
                gene_plan.gene,
                gene_plan.sample_count,
                gene_plan.batch_count,
            )
        )

        preview = gene_plan.sample_plans[:5]
        for item in preview:
            lines.append(
                "  - {} | input={} | out={} | log={}".format(
                    item.sample_id,
                    item.input_vcf,
                    item.output_file,
                    item.log_file,
                )
            )

        if len(gene_plan.sample_plans) > len(preview):
            lines.append(
                "  - ... {} more Aldy jobs".format(
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


def _build_aldy_plan_from_inputs(
    project_name: str,
    input_root: Path,
    output_dir: Path,
    discovered_inputs: List[Tuple[str, str, str, Path]],
    reference_fasta: Path,
    genome: str,
    profile: str,
) -> AldyPlan:
    gene_to_items: Dict[str, List[Tuple[str, str, Path]]] = {}

    for gene, batch_name, sample_id, input_vcf in discovered_inputs:
        gene_to_items.setdefault(gene, []).append((batch_name, sample_id, input_vcf))

    gene_plans: List[AldyGenePlan] = []
    total_jobs = 0

    for gene in sorted(gene_to_items.keys()):
        entries = sorted(gene_to_items[gene], key=lambda x: (x[0], x[1]))
        sample_plans: List[AldySamplePlan] = []
        batch_names = set()

        for batch_name, sample_id, input_vcf in entries:
            batch_names.add(batch_name)
            gene_output_dir = output_dir / gene / batch_name
            output_file = gene_output_dir / "{}_aldy.out".format(sample_id)
            log_file = gene_output_dir / "{}_aldy.log".format(sample_id)

            command = _build_aldy_command(
                gene=gene,
                input_vcf=input_vcf,
                output_file=output_file,
                log_file=log_file,
                reference_fasta=reference_fasta,
                genome=genome,
                profile=profile,
            )

            sample_plans.append(
                AldySamplePlan(
                    sample_id=sample_id,
                    batch_name=batch_name,
                    input_vcf=input_vcf,
                    output_file=output_file,
                    log_file=log_file,
                    command=command,
                )
            )

        gene_plans.append(
            AldyGenePlan(
                gene=gene,
                sample_count=len(sample_plans),
                batch_count=len(batch_names),
                sample_plans=sample_plans,
            )
        )
        total_jobs += len(sample_plans)

    return AldyPlan(
        project_name=project_name,
        input_dir=input_root,
        output_dir=output_dir,
        profile=profile,
        genome=genome,
        reference_fasta=reference_fasta,
        gene_count=len(gene_plans),
        total_jobs=total_jobs,
        gene_plans=gene_plans,
    )


# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------


def _execute_aldy_jobs(
    plan: AldyPlan,
    jobs: int,
    resume: bool,
) -> None:
    _require_aldy_executable()
    plan.output_dir.mkdir(parents=True, exist_ok=True)

    aldy_jobs = _collect_aldy_jobs(plan=plan, resume=resume)
    if not aldy_jobs:
        return

    if jobs == 1:
        for job in aldy_jobs:
            _run_single_aldy_job(job)
        return

    errors: List[str] = []
    with ThreadPoolExecutor(max_workers=jobs) as executor:
        future_to_job = {
            executor.submit(_run_single_aldy_job, job): job
            for job in aldy_jobs
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
        raise AldyExecutionError(
            "One or more Aldy jobs failed:\n{}".format("\n".join(errors))
        )


def _collect_aldy_jobs(plan: AldyPlan, resume: bool) -> List[dict]:
    jobs: List[dict] = []

    for gene_plan in plan.gene_plans:
        for sample_plan in gene_plan.sample_plans:
            sample_plan.output_file.parent.mkdir(parents=True, exist_ok=True)

            if resume and sample_plan.output_file.exists() and sample_plan.log_file.exists():
                continue

            jobs.append(
                {
                    "gene": gene_plan.gene,
                    "sample_id": sample_plan.sample_id,
                    "input_vcf": sample_plan.input_vcf,
                    "output_file": sample_plan.output_file,
                    "log_file": sample_plan.log_file,
                    "command": sample_plan.command,
                }
            )

    return jobs


def _run_single_aldy_job(job: dict) -> None:
    result = subprocess.run(job["command"], capture_output=True, text=True)
    if result.returncode != 0:
        raise AldyExecutionError(
            "Command failed: {}\nSTDOUT:\n{}\nSTDERR:\n{}".format(
                " ".join(job["command"]),
                result.stdout,
                result.stderr,
            )
        )


# ---------------------------------------------------------------------------
# Discovery helpers
# ---------------------------------------------------------------------------


def _discover_aldy_inputs(
    input_path: Path,
    requested_genes: Optional[List[str]] = None,
) -> List[Tuple[str, str, str, Path]]:
    requested = _normalize_genes(requested_genes)

    if not input_path.exists():
        raise AldyPlanError("Input path not found: {}".format(input_path))

    if input_path.is_file():
        if not input_path.name.endswith(".vcf.gz"):
            raise AldyPlanError("Input file must end with .vcf.gz: {}".format(input_path))

        gene, batch_name = _infer_gene_and_batch_for_single_vcf(
            input_vcf=input_path.resolve(),
            requested_genes=requested,
        )
        sample_id = input_path.name[:-7]
        return [(gene, batch_name, sample_id, input_path.resolve())]

    if not input_path.is_dir():
        raise AldyPlanError("Input path is not a file or directory: {}".format(input_path))

    vcfs = sorted(
        path.resolve()
        for path in input_path.rglob("*.vcf.gz")
        if not path.name.endswith(".vcf.gz.csi")
        and not path.name.endswith(".vcf.gz.tbi")
    )

    if not vcfs:
        raise AldyPlanError("No split input VCFs found under {}".format(input_path))

    discovered: List[Tuple[str, str, str, Path]] = []
    seen = set()

    for vcf_path in vcfs:
        rel = vcf_path.relative_to(input_path)
        gene, batch_name = _infer_gene_and_batch_from_relative(
            root=input_path,
            rel=rel,
            requested_genes=requested,
        )

        if requested and gene not in requested:
            continue

        sample_id = vcf_path.name[:-7]
        key = (gene, batch_name, sample_id)
        if key in seen:
            raise AldyPlanError(
                "Ambiguous split inputs: duplicate sample job for gene={}, batch={}, sample={}".format(
                    gene, batch_name, sample_id
                )
            )
        seen.add(key)
        discovered.append((gene, batch_name, sample_id, vcf_path))

    if requested:
        found_genes = set(item[0] for item in discovered)
        missing = [gene for gene in requested if gene not in found_genes]
        if missing:
            raise AldyPlanError(
                "Split input(s) not found for requested gene(s): {}".format(
                    ", ".join(missing)
                )
            )

    if not discovered:
        raise AldyPlanError("No Aldy input jobs were discovered under {}".format(input_path))

    return discovered


def _infer_gene_and_batch_for_single_vcf(
    input_vcf: Path,
    requested_genes: Optional[List[str]],
) -> Tuple[str, str]:
    parent_name = input_vcf.parent.name.upper()
    grandparent_name = input_vcf.parent.parent.name.upper() if input_vcf.parent.parent else ""

    batch_name = "batch_000"
    gene = None

    if parent_name.startswith("BATCH_"):
        batch_name = input_vcf.parent.name
        if grandparent_name:
            gene = grandparent_name

    if gene is None and requested_genes and len(requested_genes) == 1:
        gene = requested_genes[0]

    if gene is None:
        raise AldyPlanError(
            "Could not infer gene for single input VCF {}. Supply --gene/--genes or use a standard gene/batch directory structure.".format(
                input_vcf
            )
        )

    return gene, batch_name


def _infer_gene_and_batch_from_relative(
    root: Path,
    rel: Path,
    requested_genes: Optional[List[str]],
) -> Tuple[str, str]:
    parts = rel.parts

    if len(parts) >= 3:
        gene = parts[0].upper()
        batch_name = parts[1]
        return gene, batch_name

    if len(parts) == 2:
        root_name = root.name.upper()
        batch_name = parts[0]

        if requested_genes and len(requested_genes) == 1:
            gene = requested_genes[0]
        else:
            gene = root_name

        return gene, batch_name

    if len(parts) == 1:
        if requested_genes and len(requested_genes) == 1:
            return requested_genes[0], "batch_000"

        raise AldyPlanError(
            "Could not infer gene for input file {} under directory {}. Supply --gene/--genes or use a standard gene/batch directory structure.".format(
                rel,
                root,
            )
        )

    raise AldyPlanError(
        "Could not infer gene/batch structure from input path {}".format(rel)
    )


# ---------------------------------------------------------------------------
# Command helpers
# ---------------------------------------------------------------------------


def _build_aldy_command(
    gene: str,
    input_vcf: Path,
    output_file: Path,
    log_file: Path,
    reference_fasta: Path,
    genome: str,
    profile: str,
) -> List[str]:
    return [
        "aldy",
        "genotype",
        "--gene",
        gene,
        "--genome",
        genome,
        "--profile",
        profile,
        "--reference",
        str(reference_fasta),
        "--output",
        str(output_file),
        "--log",
        str(log_file),
        str(input_vcf),
    ]


def _resolve_aldy_input_path(cfg: AppConfig) -> Path:
    if cfg.aldy.input_dir is not None:
        return cfg.aldy.input_dir
    return cfg.split.output_dir


def _resolve_aldy_reference_fasta(cfg: AppConfig) -> Path:
    if cfg.aldy.reference_fasta is not None:
        return cfg.aldy.reference_fasta
    return cfg.reference.fasta


def _require_aldy_executable() -> None:
    if shutil.which("aldy") is None:
        raise AldyExecutionError(
            "Aldy executable was not found in PATH. Install Aldy or activate the correct environment."
        )


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


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
            raise AldyPlanError(
                "Requested gene is not present in the configured gene list: {}".format(
                    normalized
                )
            )
        if normalized not in selected:
            selected.append(normalized)

    if not selected:
        raise AldyPlanError("No valid genes were selected for Aldy.")

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