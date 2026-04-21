from pathlib import Path
import csv
import requests


def normalize_build(build: str) -> str:
    build = build.strip().lower()
    if build in {"hg38", "grch38"}:
        return "hg38"
    if build in {"hg19", "grch37"}:
        return "hg19"
    raise ValueError("Unsupported build: {}".format(build))


def ensembl_base_url(build: str) -> str:
    build = normalize_build(build)
    if build == "hg38":
        return "https://rest.ensembl.org"
    if build == "hg19":
        return "https://grch37.rest.ensembl.org"
    raise ValueError("Unsupported build: {}".format(build))


def load_registry_tsv(registry_tsv: str) -> list:
    path = Path(registry_tsv)
    if not path.exists():
        return []

    rows = []
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def write_registry_tsv(rows: list, registry_tsv: str) -> None:
    path = Path(registry_tsv)
    path.parent.mkdir(parents=True, exist_ok=True)

    rows = sorted(
        rows,
        key=lambda x: (
            x["build"],
            x["chrom"],
            int(x["start"]),
            int(x["end"]),
            x["gene"],
        ),
    )

    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["gene", "build", "chrom", "start", "end"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def fetch_gene_coordinate(gene: str, build: str) -> dict:
    gene = gene.strip().upper()
    build = normalize_build(build)
    url = "{}/lookup/symbol/homo_sapiens/{}?content-type=application/json".format(
        ensembl_base_url(build), gene
    )

    r = requests.get(url, timeout=30)
    if r.status_code == 404:
        raise ValueError("Gene not found in Ensembl: {}".format(gene))
    r.raise_for_status()

    data = r.json()

    chrom = str(data["seq_region_name"]).replace("chr", "").upper()
    if chrom == "M":
        chrom = "MT"

    return {
        "gene": gene,
        "build": build,
        "chrom": chrom,
        "start": str(data["start"]),
        "end": str(data["end"]),
    }


def ensure_gene_coordinates(genes: list, build: str, registry_tsv: str) -> list:
    """
    Make sure requested genes exist in the local TSV for the requested build.
    Missing genes are fetched from Ensembl once and appended to the TSV.

    Returns a list of rows for the requested genes/build.
    """
    build = normalize_build(build)
    genes = [g.strip().upper() for g in genes if g.strip()]
    genes = list(dict.fromkeys(genes))  # remove duplicates, keep order

    existing_rows = load_registry_tsv(registry_tsv)

    existing_lookup = {
        (row["gene"].upper(), normalize_build(row["build"])): row
        for row in existing_rows
    }

    missing = [
        gene for gene in genes
        if (gene, build) not in existing_lookup
    ]

    new_rows = []
    for gene in missing:
        row = fetch_gene_coordinate(gene, build)
        new_rows.append(row)
        existing_lookup[(gene, build)] = row

    if new_rows:
        all_rows = existing_rows + new_rows
        write_registry_tsv(all_rows, registry_tsv)

    resolved = [existing_lookup[(gene, build)] for gene in genes]
    return resolved