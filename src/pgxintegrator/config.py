from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Optional, Dict, Union

import yaml


class ConfigError(ValueError):
    """Raised when the user configuration is missing required values or is invalid."""


DEFAULT_HG19_TO_HG38_CHAIN = (
    Path(__file__).resolve().parent
    / "resources"
    / "chains"
    / "hg19ToHg38.over.chain.gz"
)


@dataclass
class ReferenceConfig:
    fasta: Path
    build: str


@dataclass
class InputsConfig:
    sample_manifest: Path
    dataset_manifest: Path
    gene_list: Path


@dataclass
class PrepareConfig:
    output_dir: Path
    normalize: bool = True
    harmonize: bool = True
    sample_missingness_threshold: float = 0.05
    keep_intermediate: bool = False
    resume: bool = True
    temp_dir: Optional[Path] = None


@dataclass
class LiftoverConfig:
    enabled: bool = False
    target_build: Optional[str] = None
    chain_file: Optional[Path] = None


@dataclass
class SplitConfig:
    output_dir: Path
    input_dir: Optional[Path] = None
    samples_per_batch: int = 1000
    jobs: int = 1
    resume: bool = True


@dataclass
class AldyConfig:
    output_dir: Path
    input_dir: Optional[Path] = None
    jobs: int = 1
    genome: str = "hg38"
    profile: str = "wgs"
    reference_fasta: Optional[Path] = None
    resume: bool = True


@dataclass
class AppConfig:
    project_name: str
    reference: ReferenceConfig
    inputs: InputsConfig
    prepare: PrepareConfig
    liftover: LiftoverConfig
    split: SplitConfig
    aldy: AldyConfig
    config_path: Optional[Path] = field(default=None, repr=False)

    @property
    def target_build(self) -> str:
        """Return the unified target build used across the pipeline."""
        return self.liftover.target_build or self.reference.build


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def load_config(
    config_path: Union[str, Path],
    overrides: Optional[Mapping[str, Any]] = None,
) -> AppConfig:
    """
    Load an application config from YAML and optionally apply nested overrides.

    Parameters
    ----------
    config_path:
        Path to the YAML config file.
    overrides:
        Nested mapping with values that should replace config values. Example:
        {"aldy": {"jobs": 32}, "prepare": {"resume": False}}
    """
    path = Path(config_path).expanduser().resolve()
    if not path.exists():
        raise ConfigError("Config file does not exist: {}".format(path))

    with path.open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle) or {}

    if not isinstance(raw, dict):
        raise ConfigError("Top-level YAML structure must be a mapping.")

    merged = _deep_merge(raw, dict(overrides or {}))
    cfg = _build_app_config(merged, config_path=path)
    _validate_config(cfg)
    return cfg


# ---------------------------------------------------------------------------
# Builders
# ---------------------------------------------------------------------------


def _build_app_config(raw: Mapping[str, Any], config_path: Path) -> AppConfig:
    project_name = _require_str(raw, "project_name")

    reference_raw = _require_section(raw, "reference")
    inputs_raw = _require_section(raw, "inputs")
    prepare_raw = _require_section(raw, "prepare")
    split_raw = _require_section(raw, "split")
    aldy_raw = _require_section(raw, "aldy")
    liftover_raw = raw.get("liftover", {}) or {}

    reference = ReferenceConfig(
        fasta=_to_path(_require_str(reference_raw, "fasta"), base_dir=config_path.parent),
        build=_require_str(reference_raw, "build"),
    )

    inputs = InputsConfig(
        sample_manifest=_to_path(
            _require_str(inputs_raw, "sample_manifest"),
            base_dir=config_path.parent,
        ),
        dataset_manifest=_to_path(
            _require_str(inputs_raw, "dataset_manifest"),
            base_dir=config_path.parent,
        ),
        gene_list=_to_path(
            _require_str(inputs_raw, "gene_list"),
            base_dir=config_path.parent,
        ),
    )

    prepare = PrepareConfig(
        output_dir=_to_path(_require_str(prepare_raw, "output_dir"), base_dir=config_path.parent),
        normalize=bool(prepare_raw.get("normalize", True)),
        harmonize=bool(prepare_raw.get("harmonize", True)),
        sample_missingness_threshold=float(
            prepare_raw.get("sample_missingness_threshold", 0.05)
        ),
        keep_intermediate=bool(prepare_raw.get("keep_intermediate", False)),
        resume=bool(prepare_raw.get("resume", True)),
        temp_dir=_optional_path(prepare_raw.get("temp_dir"), base_dir=config_path.parent),
    )

    requested_target_build = _optional_str(liftover_raw.get("target_build")) or reference.build
    chain_file = _optional_path(liftover_raw.get("chain_file"), base_dir=config_path.parent)
    chain_file = _resolve_default_chain_file(
        chain_file=chain_file,
        target_build=requested_target_build,
    )

    liftover = LiftoverConfig(
        enabled=bool(liftover_raw.get("enabled", False)),
        target_build=requested_target_build,
        chain_file=chain_file,
    )

    split = SplitConfig(
        output_dir=_to_path(_require_str(split_raw, "output_dir"), base_dir=config_path.parent),
        input_dir=_optional_path(split_raw.get("input_dir"), base_dir=config_path.parent),
        samples_per_batch=int(split_raw.get("samples_per_batch", 1000)),
        jobs=int(split_raw.get("jobs", 1)),
        resume=bool(split_raw.get("resume", True)),
    )

    aldy = AldyConfig(
        output_dir=_to_path(_require_str(aldy_raw, "output_dir"), base_dir=config_path.parent),
        input_dir=_optional_path(aldy_raw.get("input_dir"), base_dir=config_path.parent),
        jobs=int(aldy_raw.get("jobs", 1)),
        genome=_optional_str(aldy_raw.get("genome")) or reference.build,
        profile=_optional_str(aldy_raw.get("profile")) or "wgs",
        reference_fasta=_optional_path(
            aldy_raw.get("reference_fasta"),
            base_dir=config_path.parent,
        ),
        resume=bool(aldy_raw.get("resume", True)),
    )

    return AppConfig(
        project_name=project_name,
        reference=reference,
        inputs=inputs,
        prepare=prepare,
        liftover=liftover,
        split=split,
        aldy=aldy,
        config_path=config_path,
    )


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def _validate_config(cfg: AppConfig) -> None:
    if cfg.prepare.sample_missingness_threshold < 0 or cfg.prepare.sample_missingness_threshold > 1:
        raise ConfigError("prepare.sample_missingness_threshold must be between 0 and 1.")

    if cfg.split.samples_per_batch < 1:
        raise ConfigError("split.samples_per_batch must be at least 1.")

    if cfg.split.jobs < 1:
        raise ConfigError("split.jobs must be at least 1.")

    if cfg.aldy.jobs < 1:
        raise ConfigError("aldy.jobs must be at least 1.")

    if not cfg.reference.fasta.exists():
        raise ConfigError("Reference FASTA not found: {}".format(cfg.reference.fasta))

    if not cfg.inputs.sample_manifest.exists():
        raise ConfigError("Sample manifest not found: {}".format(cfg.inputs.sample_manifest))

    if not cfg.inputs.dataset_manifest.exists():
        raise ConfigError("Dataset manifest not found: {}".format(cfg.inputs.dataset_manifest))

    if not cfg.inputs.gene_list.exists():
        raise ConfigError("Gene list not found: {}".format(cfg.inputs.gene_list))

    if cfg.aldy.reference_fasta is not None and not cfg.aldy.reference_fasta.exists():
        raise ConfigError("Aldy reference FASTA not found: {}".format(cfg.aldy.reference_fasta))

    if cfg.liftover.target_build and cfg.liftover.target_build != cfg.reference.build:
        raise ConfigError(
            "liftover.target_build must match reference.build in the current release design."
        )

    if cfg.liftover.enabled:
        if cfg.liftover.chain_file is None:
            raise ConfigError(
                "No chain file could be resolved for liftover. Provide liftover.chain_file "
                "or place the default hg19ToHg38 chain at {}".format(DEFAULT_HG19_TO_HG38_CHAIN)
            )
        if not cfg.liftover.chain_file.exists():
            raise ConfigError(
                "Liftover chain file not found: {}".format(cfg.liftover.chain_file)
            )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _resolve_default_chain_file(
    chain_file: Optional[Path],
    target_build: str,
) -> Optional[Path]:
    if chain_file is not None:
        return chain_file

    target = target_build.strip().lower()
    if target in {"hg38", "grch38"}:
        return DEFAULT_HG19_TO_HG38_CHAIN.resolve()

    return None


def _require_section(raw: Mapping[str, Any], key: str) -> Mapping[str, Any]:
    value = raw.get(key)
    if not isinstance(value, dict):
        raise ConfigError("Missing or invalid section: {}".format(key))
    return value


def _require_str(raw: Mapping[str, Any], key: str) -> str:
    value = raw.get(key)
    if value is None:
        raise ConfigError("Missing required field: {}".format(key))
    if not isinstance(value, str) or not value.strip():
        raise ConfigError("Field must be a non-empty string: {}".format(key))
    return value.strip()


def _optional_str(value: Any) -> Optional[str]:
    if value is None:
        return None
    if not isinstance(value, str):
        raise ConfigError("Expected a string value.")
    stripped = value.strip()
    return stripped or None


def _to_path(value: str, base_dir: Path) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = (base_dir / path).resolve()
    else:
        path = path.resolve()
    return path


def _optional_path(value: Any, base_dir: Path) -> Optional[Path]:
    if value is None:
        return None
    if not isinstance(value, str) or not value.strip():
        raise ConfigError("Expected a non-empty path string.")
    return _to_path(value, base_dir=base_dir)


def _deep_merge(base: Dict[str, Any], updates: Dict[str, Any]) -> Dict[str, Any]:
    merged: Dict[str, Any] = dict(base)
    for key, value in updates.items():
        if (
            key in merged
            and isinstance(merged[key], dict)
            and isinstance(value, Mapping)
        ):
            merged[key] = _deep_merge(dict(merged[key]), dict(value))
        else:
            merged[key] = value
    return merged