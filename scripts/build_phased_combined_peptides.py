#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from neoag.haplotype import combined_haplotype_peptides
from neoag.schemas import EVENT_FIELDS, PEPTIDE_FIELDS
from neoag.utils import read_tsv, safe_id, write_tsv


def main() -> None:
    parser = argparse.ArgumentParser(description="Replace phased cis SNV components with combined-mutant peptides")
    parser.add_argument("--raw-events", required=True)
    parser.add_argument("--raw-peptides", required=True)
    parser.add_argument("--variant-peptides", required=True)
    parser.add_argument("--phasing-summary", required=True)
    parser.add_argument("--variant-key", action="append", required=True)
    parser.add_argument("--hla", action="append", required=True)
    parser.add_argument("--out-events", required=True)
    parser.add_argument("--out-peptides", required=True)
    parser.add_argument("--combined-only", help="Optional TSV containing only combined peptide-HLA rows")
    args = parser.parse_args()

    summary = json.loads(Path(args.phasing_summary).read_text())
    if summary.get("haplotype_status") != "PHASED_CIS":
        raise SystemExit("Combined reconstruction requires PHASED_CIS evidence")

    keys = set(args.variant_key)
    catalog = [row for row in read_tsv(args.variant_peptides) if (row.get("variant_key") or row.get("event_id")) in keys]
    component_ids, protein_change, windows = combined_haplotype_peptides(catalog)
    events = read_tsv(args.raw_events)
    peptides = read_tsv(args.raw_peptides)
    component_events = [row for row in events if row.get("event_id") in keys]
    if len(component_events) != len(keys):
        raise SystemExit("Not all component event IDs were found in raw_events")

    representative = dict(component_events[0])
    positions = sorted(int(row["pos"]) for row in component_events)
    phase_group = safe_id(f"phase_{representative.get('chrom')}_{positions[0]}_{positions[-1]}_{representative.get('transcript_id')}")
    component_labels = [key.split("|", 1)[-1] for key in sorted(keys)]
    combined_event_id = f"{representative.get('gene')}|{';'.join(component_labels)}|PHASED_CIS"
    combined_event = dict(representative)
    combined_event.update({
        "event_id": combined_event_id,
        "event_name": combined_event_id,
        "event_type": "SNV",
        "consequence": "phased_multi_substitution",
        "pos": str(positions[0]),
        "ref": ";".join(row.get("ref", "") for row in sorted(component_events, key=lambda r: int(r["pos"]))),
        "alt": ";".join(row.get("alt", "") for row in sorted(component_events, key=lambda r: int(r["pos"]))),
        "tumor_depth": str(summary.get("informative_fragments", "")),
        "tumor_alt_count": str(summary.get("cis_alt_fragments", "")),
        "tumor_vaf": f"{summary.get('cis_alt_fragments', 0) / max(summary.get('informative_fragments', 1), 1):.6f}",
        "phase_group_id": phase_group,
        "haplotype_status": "PHASED_CIS_COMBINED",
        "phase_support_reads": str(summary.get("cis_alt_fragments", "")),
        "phase_total_informative_reads": str(summary.get("informative_fragments", "")),
        "phase_confidence": str(summary.get("phase_confidence", "")),
        "component_event_ids": component_ids,
        "combined_protein_change": protein_change,
        "redundancy_group": phase_group,
        "source": "read_backed_phased_haplotype",
    })
    retained_events = [row for row in events if row.get("event_id") not in keys]
    retained_events.append(combined_event)

    retained_peptides = [row for row in peptides if row.get("event_id") not in keys]
    combined_rows = []
    for window in windows:
        for hla in args.hla:
            peptide_id = safe_id(f"{combined_event_id}_{window['peptide']}_{hla}")
            row = {field: "" for field in PEPTIDE_FIELDS}
            row.update({
                "peptide_id": peptide_id, "event_id": combined_event_id,
                "sample_id": representative.get("sample_id", ""), "event_type": "SNV",
                "mutation_source": representative.get("mutation_source", "SNV"),
                "peptide_consequence": "phased_multi_substitution",
                "gene": representative.get("gene", ""), "peptide": window["peptide"],
                "wildtype_peptide": window["wildtype_peptide"], "hla_allele": hla,
                "mhc_class": "I", "source_tool": "read_backed_phased_haplotype",
                "phase_group_id": phase_group, "haplotype_status": "PHASED_CIS_COMBINED",
                "phase_support_reads": str(summary.get("cis_alt_fragments", "")),
                "phase_total_informative_reads": str(summary.get("informative_fragments", "")),
                "phase_confidence": str(summary.get("phase_confidence", "")),
                "component_event_ids": component_ids, "combined_protein_change": protein_change,
                "redundancy_group": phase_group,
                "mutation_positions_in_peptide": window["mutation_positions_in_peptide"],
            })
            combined_rows.append(row)
    retained_peptides.extend(combined_rows)
    write_tsv(args.out_events, retained_events, EVENT_FIELDS)
    write_tsv(args.out_peptides, retained_peptides, PEPTIDE_FIELDS)
    if args.combined_only:
        write_tsv(args.combined_only, combined_rows, PEPTIDE_FIELDS)
    print(json.dumps({
        "combined_event_id": combined_event_id, "phase_group_id": phase_group,
        "component_events_removed": len(component_events),
        "component_peptide_hla_rows_removed": len(peptides) - len(retained_peptides) + len(combined_rows),
        "combined_unique_peptides": len(windows), "combined_peptide_hla_rows": len(combined_rows),
        "combined_protein_change": protein_change,
    }, indent=2))


if __name__ == "__main__":
    main()
