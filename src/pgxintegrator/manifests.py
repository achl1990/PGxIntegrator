from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Tuple
import csv


class ManifestError(ValueError):
    """Raised when a manifest file is missing required fields or contains invalid rows."""


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------


@dataclass
class SampleRecord:
    sample_id: str
    array_id: Optional[str] = None
    wes_id: Optional[str] = None


@dataclass
class DatasetRecord:
    dataset_id: str
    platform: str
    build: str
    layout: str
    path: Path
    chrom: Optional[str] = None
    gene: Optional[str] = None


@dataclass
class GeneNameRecord:
    gene: str


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def load_sample_manifest(path: Path) -> List[SampleRecord]:
    """
    Load and validate a sample manifest.

    Required columns:
      - sample_id

    Optional columns:
      - array_id
      - wes_id
    """
    rows = _read_tsv(path)
    _require_columns(rows.header, required=["sample_id"])

    records: List[SampleRecord] = []
    seen_sample_ids = set()

    for row_num, row in rows.records:
        sample_id = _clean_required(row, "sample_id", row_num)
        array_id = _clean_optional(row, "array_id")
        wes_id = _clean_optional(row, "wes_id")

        if sample_id in seen_sample_ids:
            raise ManifestError(
                "Duplicate sample_id in sample manifest at row {}: {}".format(
                    row_num, sample_id
                )
            )
        seen_sample_ids.add(sample_id)

        records.append(
            SampleRecord(
                sample_id=sample_id,
                array_id=array_id,
                wes_id=wes_id,
            )
        )

    if not records:
        raise ManifestError("Sample manifest contains no data rows: {}".format(path))

    return records


def load_dataset_manifest(path: Path) -> List[DatasetRecord]:
    """
    Load and validate a dataset manifest.

    Required columns:
      - dataset_id
      - platform
      - build
      - layout
      - path

    Additional required columns depend on layout:

    For layout = cohort_by_chrom:
      - chrom

    For layout = gene_by_file:
      - gene

    Supported platform values:
      - wes
      - array

    Supported layout values:
      - cohort_by_chrom
      - gene_by_file
    """
    rows = _read_tsv(path)
    _require_columns(
        rows.header,
        required=["dataset_id", "platform", "build", "layout", "path"],
    )

    records: List[DatasetRecord] = []
    seen_keys = set()
    manifest_dir = path.expanduser().resolve().parent

    for row_num, row in rows.records:
        dataset_id = _clean_required(row, "dataset_id", row_num)
        platform = _normalize_platform(_clean_required(row, "platform", row_num), row_num)
        build = _clean_required(row, "build", row_num)
        layout = _normalize_layout(_clean_required(row, "layout", row_num), row_num)
        input_path = _resolve_path(_clean_required(row, "path", row_num), manifest_dir)

        if not input_path.exists():
            raise ManifestError(
                "Input file not found in dataset manifest at row {}: {}".format(
                    row_num, input_path
                )
            )

        chrom: Optional[str] = None
        gene: Optional[str] = None

        if layout == "cohort_by_chrom":
            chrom = _normalize_chrom(_clean_required(row, "chrom", row_num), row_num)
        elif layout == "gene_by_file":
            gene = _clean_required(row, "gene", row_num).upper()

        key = (dataset_id, platform, build, layout, chrom, gene, str(input_path))
        if key in seen_keys:
            raise ManifestError(
                "Duplicate dataset entry in dataset manifest at row {}: {}".format(
                    row_num, key
                )
            )
        seen_keys.add(key)

        records.append(
            DatasetRecord(
                dataset_id=dataset_id,
                platform=platform,
                build=build,
                layout=layout,
                path=input_path,
                chrom=chrom,
                gene=gene,
            )
        )

    if not records:
        raise ManifestError("Dataset manifest contains no data rows: {}".format(path))

    return records


def load_gene_list(path: Path) -> List[GeneNameRecord]:
    """
    Load and validate a user-facing gene list.

    Required columns:
      - gene

    Example:
        gene
        CYP2C19
        CYP2D6
        CYP3A5
    """
    rows = _read_tsv(path)
    _require_columns(rows.header, required=["gene"])

    records: List[GeneNameRecord] = []
    seen_genes = set()

    for row_num, row in rows.records:
        gene = _clean_required(row, "gene", row_num).upper()

        if gene in seen_genes:
            raise ManifestError(
                "Duplicate gene in gene list at row {}: {}".format(row_num, gene)
            )
        seen_genes.add(gene)

        records.append(GeneNameRecord(gene=gene))

    if not records:
        raise ManifestError("Gene list contains no data rows: {}".format(path))

    return records


def group_dataset_records_by_platform(
    records: List[DatasetRecord],
) -> Dict[str, List[DatasetRecord]]:
    """
    Group dataset records by platform.
    """
    grouped: Dict[str, List[DatasetRecord]] = {}
    for record in records:
        grouped.setdefault(record.platform, []).append(record)
    return grouped


def summarize_chrom_datasets(
    records: List[DatasetRecord],
) -> Dict[str, Dict[str, DatasetRecord]]:
    """
    Build a nested mapping for chromosome-based datasets:

        summary[platform][chrom] = DatasetRecord

    Only records with layout == 'cohort_by_chrom' are included.
    """
    summary: Dict[str, Dict[str, DatasetRecord]] = {}

    for record in records:
        if record.layout != "cohort_by_chrom":
            continue

        if record.chrom is None:
            raise ManifestError(
                "Chromosome-based dataset record is missing chrom: {}".format(record)
            )

        platform_bucket = summary.setdefault(record.platform, {})
        if record.chrom in platform_bucket:
            raise ManifestError(
                "Duplicate chromosome dataset for platform {} chrom {}".format(
                    record.platform, record.chrom
                )
            )
        platform_bucket[record.chrom] = record

    return summary


def summarize_gene_datasets(
    records: List[DatasetRecord],
) -> Dict[str, Dict[str, DatasetRecord]]:
    """
    Build a nested mapping for gene-based datasets:

        summary[platform][gene] = DatasetRecord

    Only records with layout == 'gene_by_file' are included.
    """
    summary: Dict[str, Dict[str, DatasetRecord]] = {}

    for record in records:
        if record.layout != "gene_by_file":
            continue

        if record.gene is None:
            raise ManifestError(
                "Gene-based dataset record is missing gene: {}".format(record)
            )

        platform_bucket = summary.setdefault(record.platform, {})
        if record.gene in platform_bucket:
            raise ManifestError(
                "Duplicate gene dataset for platform {} gene {}".format(
                    record.platform, record.gene
                )
            )
        platform_bucket[record.gene] = record

    return summary


def list_gene_names(records: List[GeneNameRecord]) -> List[str]:
    """
    Return gene names as an ordered list.
    """
    return [record.gene for record in records]


# ---------------------------------------------------------------------------
# Internal TSV reader
# ---------------------------------------------------------------------------


@dataclass
class _TSVRows:
    header: List[str]
    records: List[Tuple[int, Mapping[str, str]]]


def _read_tsv(path: Path) -> _TSVRows:
    resolved = path.expanduser().resolve()
    if not resolved.exists():
        raise ManifestError("Manifest file not found: {}".format(resolved))

    with resolved.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ManifestError("Manifest file has no header: {}".format(resolved))

        header = [field.strip() for field in reader.fieldnames]
        records: List[Tuple[int, Mapping[str, str]]] = []

        for row_index, row in enumerate(reader, start=2):
            normalized_row = {}
            for key, value in row.items():
                if key is None:
                    continue
                normalized_row[key.strip()] = value
            records.append((row_index, normalized_row))

    return _TSVRows(header=header, records=records)


# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------


def _require_columns(header: List[str], required: List[str]) -> None:
    header_set = set(header)
    missing = [column for column in required if column not in header_set]
    if missing:
        raise ManifestError(
            "Manifest is missing required column(s): {}".format(", ".join(missing))
        )


def _clean_required(row: Mapping[str, str], key: str, row_num: int) -> str:
    value = row.get(key)
    if value is None or not str(value).strip():
        raise ManifestError(
            "Missing required value for '{}' at row {}".format(key, row_num)
        )
    return str(value).strip()


def _clean_optional(row: Mapping[str, str], key: str) -> Optional[str]:
    value = row.get(key)
    if value is None:
        return None
    cleaned = str(value).strip()
    return cleaned or None


def _normalize_platform(value: str, row_num: int) -> str:
    platform = value.strip().lower()
    allowed = {"wes", "array"}
    if platform not in allowed:
        raise ManifestError(
            "Unsupported platform at row {}: {}. Allowed values: {}".format(
                row_num, value, ", ".join(sorted(allowed))
            )
        )
    return platform


def _normalize_layout(value: str, row_num: int) -> str:
    layout = value.strip().lower()
    allowed = {"cohort_by_chrom", "gene_by_file"}
    if layout not in allowed:
        raise ManifestError(
            "Unsupported layout at row {}: {}. Allowed values: {}".format(
                row_num, value, ", ".join(sorted(allowed))
            )
        )
    return layout


def _normalize_chrom(value: str, row_num: int) -> str:
    chrom = value.strip()
    if not chrom:
        raise ManifestError("Empty chromosome value at row {}".format(row_num))

    chrom_lower = chrom.lower()
    if chrom_lower.startswith("chr"):
        chrom = chrom[3:]

    chrom = chrom.upper()
    if chrom == "M":
        chrom = "MT"

    return chrom


def _resolve_path(value: str, base_dir: Path) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = (base_dir / path).resolve()
    else:
        path = path.resolve()
    return path