#!/usr/bin/env python3
from __future__ import annotations

import csv
import os
from pathlib import Path

import pandas as pd
import snaf


def required_env(name: str) -> str:
    value = os.environ.get(name, "").strip()
    if not value:
        raise RuntimeError(f"{name} is required")
    return value


def read_hla(path: Path) -> list[str]:
    values: list[str] = []
    for token in path.read_text(encoding="utf-8").replace(",", "\n").splitlines():
        allele = token.strip()
        if allele and allele.startswith("HLA-"):
            values.append(allele)
    if not values:
        raise RuntimeError(f"no HLA class-I alleles found in {path}")
    return values


def main() -> int:
    outdir = Path(required_env("NEOAG_SNAF_OUTDIR"))
    matrix_path = Path(required_env("NEOAG_SNAF_MATRIX"))
    db_dir = Path(required_env("NEOAG_SNAF_DB"))
    hla_file = Path(required_env("NEOAG_SNAF_HLA_FILE"))
    sample_id = required_env("NEOAG_SNAF_SAMPLE_ID")
    cores = int(os.environ.get("NEOAG_SNAF_CORES", "8"))
    netmhcpan = os.environ.get("NEOAG_NETMHCPAN_BIN") or None
    binding_method = os.environ.get(
        "NEOAG_SNAF_BINDING_METHOD", "netMHCpan" if netmhcpan else "MHCflurry"
    )

    expected = [
        db_dir / "Alt91_db/Hs_Ensembl_exon_add_col.txt",
        db_dir / "Alt91_db/mRNA-ExonIDs.txt",
        db_dir / "Alt91_db/Hs_gene-seq-2000_flank.fa",
        db_dir / "controls/GTEx_junction_counts.h5ad",
    ]
    missing = [str(path) for path in expected if not path.is_file()]
    if missing:
        raise RuntimeError("incomplete SNAF reference: " + ", ".join(missing))

    outdir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(matrix_path, sep="\t", index_col=0)
    candidates = [column for column in df.columns if "replicate" not in str(column).lower()]
    if not candidates:
        raise RuntimeError(f"no primary sample column found in {matrix_path}")
    df = df.loc[:, [candidates[0]]]
    df.columns = [sample_id]
    df = df.loc[~df.index.duplicated(keep="first")]
    df.to_csv(outdir / "snaf_junction_count_matrix.tsv", sep="\t")

    hlas = read_hla(hla_file)
    if binding_method.lower() == "netmhcpan" and not netmhcpan:
        raise RuntimeError("NEOAG_NETMHCPAN_BIN is required for SNAF NetMHCpan mode")
    snaf.initialize(
        df=df,
        db_dir=str(db_dir),
        gtex_mode="count",
        binding_method=binding_method,
        software_path=netmhcpan,
    )
    query = snaf.JunctionCountMatrixQuery(
        junction_count_matrix=df,
        cores=cores,
        outdir=str(outdir),
        filter_mode="maxmin",
    )
    query.run(hlas=[hlas], outdir=str(outdir))
    snaf.JunctionCountMatrixQuery.generate_results(
        path=str(outdir / "after_prediction.p"), outdir=str(outdir)
    )

    source = outdir / "T_candidates" / f"T_antigen_candidates_{sample_id}.txt"
    if not source.is_file():
        source = outdir / "T_candidates" / "T_antigen_candidates_all.txt"
    if not source.is_file():
        raise RuntimeError("SNAF completed without a T-antigen candidate table")
    with source.open(encoding="utf-8", errors="replace", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    fields = [
        "sample_id", "event_id", "gene", "chrom", "start", "end", "strand",
        "junction_reads", "peptide", "hla_allele", "binding_rank",
        "immunogenicity", "tumor_specificity_mean", "tumor_specificity_mle",
        "source_tool", "evidence_status",
    ]
    output = outdir / "snaf_candidates.tsv"
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            coord = str(row.get("coord", ""))
            chrom = coord.split(":", 1)[0] if ":" in coord else ""
            span = coord.split(":", 1)[1].split("(", 1)[0] if ":" in coord else ""
            strand = coord.rsplit("(", 1)[1].rstrip(")") if "(" in coord else ""
            start, end = (span.split("-", 1) + [""])[:2] if span else ("", "")
            writer.writerow({
                "sample_id": sample_id,
                "event_id": row.get("uid", ""),
                "gene": row.get("symbol", ""),
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "junction_reads": row.get("junction_count", ""),
                "peptide": row.get("peptide", ""),
                "hla_allele": row.get("hla", ""),
                "binding_rank": row.get("binding_affinity", ""),
                "immunogenicity": row.get("immunogenicity", ""),
                "tumor_specificity_mean": row.get("tumor_specificity_mean", ""),
                "tumor_specificity_mle": row.get("tumor_specificity_mle", ""),
                "source_tool": "SNAF",
                "evidence_status": "SNAF_GTEX_SUPPORTED",
            })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
