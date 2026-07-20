#!/usr/bin/env python3
"""Score + rank peptides directly from a wide `variant_peptides*annotated.tsv`
table (extract-variant-peptides output already annotated with NetMHCpan /
MHCflurry / PRIME / BigMHC / IEDB / NetMHCstabPan columns) plus a matching
`raw_events.tsv` (EVENT_FIELDS schema).

This bypasses the multi-file CLI orchestration (`run-v03`, which expects
separate raw netmhcpan.xls / mhcflurry.csv tool outputs) and instead calls
the same underlying scoring functions directly, using the upstream tool
columns that are already embedded in the wide table -- so no evidence is
silently downgraded to placeholder values.

Usage:
    python score_from_annotated_table.py \
        --annotated variant_peptides_annotated.tsv \
        --raw-events raw_events.tsv \
        --profile default \
        --outdir out/

Outputs (written to --outdir):
    ranked_peptides.tsv  -- one row per peptide x HLA, sorted by efficacy_score desc
    ranked_events.tsv    -- one row per variant event, sorted by event_score desc
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path


def _locate_neoag_src(explicit: str | None) -> Path:
    """Find the neo repo's `src/` directory so `import neoag_v03` works,
    without assuming this script lives at a fixed relative path inside the
    repo (this script is bundled with the SKILL and may be invoked from
    outside the repo checkout).

    Resolution order: --neoag-root arg > $NEOAG_ROOT env var > walking up
    from the current working directory looking for a `src/neoag_v03` dir.
    """
    candidates = []
    if explicit:
        candidates.append(Path(explicit))
    if os.environ.get("NEOAG_ROOT"):
        candidates.append(Path(os.environ["NEOAG_ROOT"]))
    cwd = Path.cwd()
    candidates.extend([cwd, *cwd.parents])
    for c in candidates:
        src = Path(c) / "src"
        if (src / "neoag_v03").is_dir():
            return src
    raise SystemExit(
        "Could not locate the neo repo's src/neoag_v03 package. "
        "Pass --neoag-root /path/to/neo, or set NEOAG_ROOT=/path/to/neo, "
        "or run this script from inside (or below) the neo repo checkout."
    )


_ap_pre = argparse.ArgumentParser(add_help=False)
_ap_pre.add_argument("--neoag-root")
_pre_args, _ = _ap_pre.parse_known_args()
sys.path.insert(0, str(_locate_neoag_src(_pre_args.neoag_root)))

from neoag_v03.config import load_profile
from neoag_v03.model_layers import enrich_event_layers, enrich_peptide_layers
from neoag_v03.safety import apply_event_safety, apply_peptide_safety
from neoag_v03.scoring_v03 import score_event, score_peptide
from neoag_v03.presentation import grade as pres_grade
from neoag_v03.utils import read_tsv, write_tsv, norm_rank, clamp, to_float, safe_id


def build_presentation_dict(row: dict, profile: dict) -> dict:
    """Compute binding_evidence_score / presentation_evidence_score / grade
    from the wide table's real upstream-tool columns, using the exact same
    formula as presentation.py::build_presentation_evidence (norm_rank +
    profile-configurable weighted average), instead of the CLI's file-based
    netmhcpan_evidence.tsv/mhcflurry_evidence.tsv merge step.
    """
    w = profile.get("presentation_weights", {})
    w_ba = float(w.get("netmhcpan_ba", 0.25))
    w_el = float(w.get("netmhcpan_el", 0.35))
    w_mhcf = float(w.get("mhcflurry_presentation", 0.30))
    w_proc = float(w.get("mhcflurry_processing", 0.10))

    ba = row.get("netmhcpan_mt_rank_ba", "")
    el = row.get("netmhcpan_mt_rank_el", "")
    pct = row.get("mhcflurry_mt_affinity_percentile", "")
    proc = row.get("mhcflurry_mt_processing_score", "")
    pres = row.get("mhcflurry_mt_presentation_score", "")

    ba_s = norm_rank(ba) if ba not in ("", None) else None
    el_s = norm_rank(el) if el not in ("", None) else None
    pct_s = norm_rank(pct) if pct not in ("", None) else None
    proc_s = clamp(to_float(proc, -1)) if proc not in ("", None) else None
    pres_s = clamp(to_float(pres, -1)) if pres not in ("", None) else None

    binding_parts = [x for x in [ba_s, pct_s] if x is not None]
    binding = max(binding_parts) if binding_parts else 0.0

    num = den = 0.0
    for val, wt in [(ba_s, w_ba), (el_s, w_el), (pres_s, w_mhcf), (proc_s, w_proc)]:
        if val is not None:
            num += val * wt
            den += wt
    presentation = (num / den) if den else 0.0
    complete = min(1.0, den / (w_ba + w_el + w_mhcf + w_proc)) if den else 0.0

    return {
        "netmhcpan_ba_rank": str(to_float(ba, 99.0)),
        "netmhcpan_el_rank": str(to_float(el, 99.0)),
        "netmhcstabpan_score": row.get("netmhcstabpan_score", ""),
        "netmhcstabpan_rank": row.get("netmhcstabpan_rank", ""),
        "mhcflurry_affinity_percentile": str(to_float(pct, 99.0)),
        "mhcflurry_processing_score": str(to_float(proc, 0.0)),
        "mhcflurry_presentation_score": str(to_float(pres, 0.0)),
        "binding_evidence_score": f"{binding:.4f}",
        "presentation_evidence_score": f"{presentation:.4f}",
        "presentation_evidence_grade": pres_grade(binding, presentation, complete),
        # BigMHC_IM is the immunogenicity source used here; PRIME/IEDB are
        # also carried through in case the active profile weights them.
        "bigmhc_im_score": row.get("bigmhc_im_score", ""),
        "prime_score": row.get("prime_score", ""),
        "iedb_immunogenicity_score": row.get("iedb_immunogenicity_score", ""),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotated", required=True, help="variant_peptides_annotated.tsv")
    ap.add_argument("--raw-events", required=True, help="raw_events.tsv (EVENT_FIELDS schema)")
    ap.add_argument("--profile", default="default")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--neoag-root", help="Path to the neo repo checkout (containing src/neoag_v03). "
                     "Also settable via NEOAG_ROOT env var; auto-detected if omitted and this "
                     "script is run from inside the repo.")
    args = ap.parse_args()

    profile = load_profile(args.profile)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    events_raw = {r["event_id"]: r for r in read_tsv(args.raw_events)}
    annotated_rows = read_tsv(args.annotated)

    # --- events: enrich + safety + score, once per event_id ---
    scored_events: dict[str, dict] = {}
    for eid, ev in events_raw.items():
        e = dict(ev)
        e = enrich_event_layers(e)
        e = apply_event_safety(e, profile, {})
        e = score_event(e, profile)
        scored_events[eid] = e

    # --- peptides: one row per (peptide, hla) pair, already 1:1 in the wide table ---
    normal_ligands = {
        (r.get("mutant_peptide") or "").upper()
        for r in annotated_rows
        if str(r.get("in_normal_proteome", "")).strip().lower() == "yes"
    }

    scored_peptides = []
    skipped = 0
    for row in annotated_rows:
        eid = row.get("variant_key", "")
        event = scored_events.get(eid)
        if event is None:
            skipped += 1
            continue
        peptide = {
            "peptide_id": row.get("peptide_id") or safe_id(f"{eid}_{row.get('hla_allele','')}_{row.get('mutant_peptide','')}"),
            "event_id": eid,
            "sample_id": event.get("sample_id", ""),
            "event_type": event.get("event_type", ""),
            "gene": row.get("gene", ""),
            "peptide": row.get("mutant_peptide", ""),
            "wildtype_peptide": row.get("wildtype_peptide", ""),
            "crosses_junction": row.get("crosses_junction", ""),
            "contains_novel_aa": row.get("contains_novel_aa", ""),
            "rna_junction_reads": row.get("rna_junction_reads", ""),
            "hla_allele": row.get("hla_allele", ""),
            "mhc_class": "II" if any(x in (row.get("hla_allele") or "").upper() for x in ("DR", "DQ", "DP")) else "I",
            "wildtype_binding_rank": row.get("netmhcpan_wt_rank_ba", "99"),
            "self_similarity_score": "0.0",
        }
        peptide = enrich_peptide_layers(peptide, event)
        peptide = apply_peptide_safety(peptide, event, profile, normal_ligands)

        presentation = build_presentation_dict(row, profile)
        summary = {
            "mhc_i_integrity_score": event.get("appm_mhc_i_integrity", "") or "1.0",
            "mhc_ii_integrity_score": event.get("appm_mhc_ii_integrity", "") or "1.0",
            "hla_loh_alleles": "",
        }
        peptide = score_peptide(peptide, event, profile, presentation, summary)
        scored_peptides.append(peptide)

    scored_peptides.sort(key=lambda p: to_float(p.get("efficacy_score"), 0.0), reverse=True)
    for i, p in enumerate(scored_peptides, 1):
        p["pathogenic_rank"] = str(i)

    events_list = sorted(scored_events.values(), key=lambda e: to_float(e.get("event_score"), 0.0), reverse=True)

    peptide_fields = list(scored_peptides[0].keys()) if scored_peptides else []
    event_fields = list(events_list[0].keys()) if events_list else []
    write_tsv(outdir / "ranked_peptides.tsv", scored_peptides, peptide_fields)
    write_tsv(outdir / "ranked_events.tsv", events_list, event_fields)

    print(f"Scored {len(scored_peptides)} peptide x HLA rows across {len(events_list)} events.")
    if skipped:
        print(f"WARNING: skipped {skipped} rows whose variant_key had no matching event_id in --raw-events.")
    print(f"Wrote: {outdir / 'ranked_peptides.tsv'}")
    print(f"Wrote: {outdir / 'ranked_events.tsv'}")


if __name__ == "__main__":
    main()
