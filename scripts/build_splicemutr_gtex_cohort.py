#!/usr/bin/env python3
"""Build a traceable SpliceMutr cohort from a GTEx junction-count matrix."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any


DEFAULT_API = "https://gtexportal.org/api/v2/dataset/sample"


def get_json(url: str, retries: int = 5) -> dict[str, Any]:
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "open-neo/SpliceMutr"})
            with urllib.request.urlopen(request, timeout=120) as response:
                return json.load(response)
        except Exception as error:  # Network failures are retried and then surfaced.
            last_error = error
            if attempt + 1 < retries:
                time.sleep(min(30, 2 ** attempt))
    raise RuntimeError(f"GTEx metadata request failed after {retries} attempts: {last_error}")


def fetch_metadata(api: str, dataset: str, tissue: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    page = 0
    while True:
        query = urllib.parse.urlencode(
            {"datasetId": dataset, "tissueSiteDetailId": tissue, "page": page, "pageSize": 250}
        )
        payload = get_json(f"{api}?{query}")
        rows.extend(payload.get("data", []))
        paging = payload.get("paging_info", {})
        if page + 1 >= int(paging.get("numberOfPages", 1)):
            return rows
        page += 1


def load_matrix_samples(gct: Path) -> list[str]:
    with gzip.open(gct, "rt", encoding="utf-8") as handle:
        handle.readline()
        handle.readline()
        return handle.readline().rstrip("\n").split("\t")[2:]


def select_samples(
    metadata: list[dict[str, Any]], matrix_samples: set[str], count: int, min_rin: float
) -> list[dict[str, Any]]:
    eligible = [
        row for row in metadata
        if row.get("dataType") == "RNASEQ"
        and float(row.get("rin") or 0) >= min_rin
        and row.get("sampleId") in matrix_samples
    ]
    eligible.sort(key=lambda row: (-float(row.get("rin") or 0), row.get("sampleId", "")))
    selected: list[dict[str, Any]] = []
    donors: set[str] = set()
    for row in eligible:
        donor = str(row.get("subjectId") or "")
        if not donor or donor in donors:
            continue
        selected.append(row)
        donors.add(donor)
        if len(selected) == count:
            return selected
    raise ValueError(
        f"only {len(selected)} distinct-donor RNASEQ samples with RIN >= {min_rin} "
        f"were present in the junction matrix; requested {count}"
    )


def split_junction_name(name: str) -> tuple[str, int, int, str]:
    try:
        locus, strand = name.rsplit(":", 1)
        chrom, interval = locus.split(":", 1)
        start, end = interval.split("-", 1)
        return chrom, int(start), int(end), strand
    except (ValueError, TypeError) as error:
        raise ValueError(f"invalid GTEx junction identifier: {name!r}") from error


def export_junctions(
    gct: Path, selected: list[dict[str, Any]], outdir: Path, min_reads: int
) -> dict[str, int]:
    outdir.mkdir(parents=True, exist_ok=True)
    sample_ids = [str(row["sampleId"]) for row in selected]
    outputs = {sample: (outdir / f"{sample}.junc").open("w", encoding="utf-8") for sample in sample_ids}
    retained = {sample: 0 for sample in sample_ids}
    try:
        with gzip.open(gct, "rt", encoding="utf-8") as handle:
            handle.readline()
            handle.readline()
            header = handle.readline().rstrip("\n").split("\t")
            indices = {sample: header.index(sample) for sample in sample_ids}
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                chrom, start, end, strand = split_junction_name(fields[0])
                for sample, index in indices.items():
                    reads = int(float(fields[index]))
                    if reads < min_reads:
                        continue
                    retained[sample] += 1
                    values = [chrom, str(start - 1), str(end + 1),
                              f"JUNC{retained[sample]}", str(reads), strand]
                    outputs[sample].write("\t".join(values) + "\n")
    finally:
        for handle in outputs.values():
            handle.close()
    return retained


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gct", required=True, type=Path)
    parser.add_argument("--tumor-sample-id", required=True)
    parser.add_argument("--tumor-star-sj", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--tissue-id", default="Adipose_Visceral_Omentum")
    parser.add_argument("--normal-count", type=int, default=12)
    parser.add_argument("--min-rin", type=float, default=7.0)
    parser.add_argument("--min-junction-reads", type=int, default=3)
    parser.add_argument("--metadata-dataset", default="gtex_v10")
    parser.add_argument("--junction-dataset", default="GTEx_V11")
    parser.add_argument("--api-url", default=DEFAULT_API)
    parser.add_argument("--metadata-json", type=Path)
    args = parser.parse_args()

    gct = args.gct.expanduser().resolve()
    tumor_sj = args.tumor_star_sj.expanduser().resolve()
    if not gct.is_file() or not tumor_sj.is_file():
        raise FileNotFoundError("the GTEx GCT and tumor STAR junction file must both exist")
    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    if args.metadata_json:
        metadata = json.loads(args.metadata_json.read_text(encoding="utf-8"))
        if isinstance(metadata, dict):
            metadata = metadata.get("data", [])
    else:
        metadata = fetch_metadata(args.api_url, args.metadata_dataset, args.tissue_id)
        (outdir / "gtex_metadata.raw.json").write_text(
            json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
        )

    matrix_samples = load_matrix_samples(gct)
    selected = select_samples(metadata, set(matrix_samples), args.normal_count, args.min_rin)
    junc_dir = outdir / "normal_junctions"
    retained = export_junctions(gct, selected, junc_dir, args.min_junction_reads)

    selection_fields = ["sample_id", "subject_id", "tissue", "rin", "sex", "age_bracket",
                        "retained_junctions", "normal_source", "normal_match"]
    selection_rows = [
        {"sample_id": row["sampleId"], "subject_id": row.get("subjectId", ""),
         "tissue": args.tissue_id, "rin": row.get("rin", ""), "sex": row.get("sex", ""),
         "age_bracket": row.get("ageBracket", ""),
         "retained_junctions": retained[str(row["sampleId"])],
         "normal_source": args.junction_dataset, "normal_match": "PUBLIC_PROXY"}
        for row in selected
    ]
    with (outdir / "gtex_normal_selection.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=selection_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(selection_rows)

    sample_fields = ["sample_id", "role", "star_sj", "junction_file", "normal_source",
                     "normal_match", "tissue", "dataset"]
    with (outdir / "splicemutr.samples.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=sample_fields, delimiter="\t")
        writer.writeheader()
        for row in selected:
            sample = str(row["sampleId"])
            writer.writerow({"sample_id": sample, "role": "normal",
                             "junction_file": junc_dir / f"{sample}.junc",
                             "normal_source": args.junction_dataset,
                             "normal_match": "PUBLIC_PROXY", "tissue": args.tissue_id,
                             "dataset": args.junction_dataset})
        writer.writerow({"sample_id": args.tumor_sample_id, "role": "tumor",
                         "star_sj": tumor_sj})

    manifest = {
        "status": "READY", "normal_source": args.junction_dataset,
        "metadata_source": args.metadata_dataset, "metadata_proxy_for_junction_dataset": True,
        "normal_match": "PUBLIC_PROXY", "tissue": args.tissue_id,
        "distinct_donors": len(selected), "min_rin": args.min_rin,
        "min_junction_reads": args.min_junction_reads,
        "missing_public_normal_junction_is_not_matched_normal_negative": True,
        "gct": str(gct), "tumor_star_sj": str(tumor_sj),
        "sample_sheet": str(outdir / "splicemutr.samples.tsv"),
    }
    (outdir / "splicemutr_gtex_cohort.manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
