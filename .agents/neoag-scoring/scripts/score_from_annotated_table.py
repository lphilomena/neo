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


def _first_existing(*candidates: Path) -> Path | None:
    for c in candidates:
        if c and c.exists():
            return c
    return None


def discover_run_dir(run_dir: Path) -> dict[str, Path | None]:
    """Inventory a real neoag pipeline output directory (see the example
    tree this skill ships with) and return the paths it can find for each
    piece of evidence, or None if that piece is missing. Callers use this to
    tell the user up front what can be scored with real data vs what will
    fall back to defaults -- never silently.
    """
    return {
        "ranked_peptides": _first_existing(run_dir / "scoring" / "ranked_peptides.v03.tsv"),
        "ranked_events": _first_existing(run_dir / "scoring" / "ranked_events.v03.tsv"),
        "raw_events": _first_existing(run_dir / "parsed" / "raw_events.tsv", run_dir / "upstream" / "parsed" / "raw_events.tsv"),
        "raw_peptides": _first_existing(run_dir / "parsed" / "raw_peptides.tsv", run_dir / "upstream" / "parsed" / "raw_peptides.tsv"),
        "annotated_wide": _first_existing(
            run_dir / "upstream" / "tools" / "variant_peptides.annotated.tsv",
            run_dir / "upstream" / "tools" / "variant_peptides.tsv",
        ),
        "presentation_evidence": _first_existing(run_dir / "presentation" / "presentation_evidence.tsv"),
        "ccf": _first_existing(run_dir / "clonality" / "ccf_2.tsv", run_dir / "clonality" / "ccf_lite.tsv"),
        "escape_flags": _first_existing(run_dir / "immune_escape" / "peptide_escape_flags.tsv"),
        "appm_summary": _first_existing(run_dir / "appm" / "appm_summary.tsv"),
        "safety_event": _first_existing(run_dir / "safety" / "event_safety.tsv"),
        "safety_peptide": _first_existing(run_dir / "safety" / "peptide_safety.tsv"),
    }


def print_inventory_report(found: dict[str, Path | None]) -> None:
    labels = {
        "ranked_peptides": "已有排序表 scoring/ranked_peptides.v03.tsv",
        "ranked_events": "已有事件排序表 scoring/ranked_events.v03.tsv",
        "raw_events": "raw_events.tsv (打分必需)",
        "raw_peptides": "raw_peptides.tsv (若无宽表annotated表则需要)",
        "annotated_wide": "variant_peptides.annotated.tsv 宽表 (若无raw_peptides.tsv则需要)",
        "presentation_evidence": "真实呈递证据 presentation/presentation_evidence.tsv",
        "ccf": "真实克隆性证据 clonality/ccf_2.tsv 或 ccf_lite.tsv",
        "escape_flags": "真实免疫逃逸证据 immune_escape/peptide_escape_flags.tsv",
        "appm_summary": "APPM通路完整性汇总 appm/appm_summary.tsv",
        "safety_event": "事件层安全性证据 safety/event_safety.tsv",
        "safety_peptide": "肽段层安全性证据 safety/peptide_safety.tsv",
    }
    print("=== 目录清点结果 ===")
    for key, label in labels.items():
        mark = "[found]" if found.get(key) else "[missing]"
        print(f"  {mark:10s} {label}")


def load_presentation_evidence(path: Path) -> dict[str, dict]:
    """Key by peptide_hla_key if present, else safe_id(peptide+hla_allele) --
    matches presentation.py::by_key so this is a drop-in replacement for the
    wide-table recompute when a real presentation_evidence.tsv exists."""
    out = {}
    for r in read_tsv(path):
        key = r.get("peptide_hla_key") or safe_id(f"{r.get('peptide','')}_{r.get('hla_allele','')}")
        out[key] = r
    return out


def load_ccf_overrides(path: Path) -> dict[str, dict]:
    """Key by event_id. Works for both ccf_2.tsv (CCF_FIELDS) and
    ccf_lite.tsv (CCF_LITE_FIELDS) -- both carry ccf_status/clonality_multiplier."""
    return {r.get("event_id", ""): r for r in read_tsv(path) if r.get("event_id")}


def load_escape_overrides(path: Path) -> dict[str, dict]:
    """Key by peptide_id. Reads whatever columns are present rather than a
    fixed schema, since peptide_escape_flags.tsv's exact column set has
    varied across neoag_v03 versions; only escape_multiplier/escape_severity
    (falling back to escape_status if escape_severity isn't present) and
    priority_cap are actually consumed downstream."""
    return {r.get("peptide_id", ""): r for r in read_tsv(path) if r.get("peptide_id")}


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
    ap.add_argument("--run-dir", help="Root of a real neoag pipeline output directory "
                     "(the one with scoring/, presentation/, appm/, clonality/, immune_escape/, "
                     "parsed/, upstream/ subfolders). If given, --annotated/--raw-events/"
                     "--presentation-evidence/--ccf/--escape-flags are auto-discovered from it "
                     "and don't need to be passed separately.")
    ap.add_argument("--force", action="store_true",
                     help="With --run-dir: proceed with a fresh scoring run even if "
                     "scoring/ranked_peptides.v03.tsv already exists there.")
    ap.add_argument("--annotated", help="variant_peptides_annotated.tsv (wide table). "
                     "Required unless --run-dir finds one.")
    ap.add_argument("--raw-events", help="raw_events.tsv (EVENT_FIELDS schema). "
                     "Required unless --run-dir finds one.")
    ap.add_argument("--presentation-evidence", help="Optional real presentation/presentation_evidence.tsv "
                     "-- used instead of recomputing binding/presentation scores from the wide table.")
    ap.add_argument("--ccf", help="Optional real clonality/ccf_2.tsv or ccf_lite.tsv "
                     "-- overrides ccf_status/clonality_multiplier already in raw_events.tsv.")
    ap.add_argument("--escape-flags", help="Optional real immune_escape/peptide_escape_flags.tsv "
                     "-- without this, escape_multiplier stays at 1.0 (not evaluated) for every peptide.")
    ap.add_argument("--profile", default="default")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--neoag-root", help="Path to the neo repo checkout (containing src/neoag_v03). "
                     "Also settable via NEOAG_ROOT env var; auto-detected if omitted and this "
                     "script is run from inside the repo.")
    args = ap.parse_args()

    annotated_path = args.annotated
    raw_events_path = args.raw_events
    presentation_evidence_path = args.presentation_evidence
    ccf_path = args.ccf
    escape_flags_path = args.escape_flags

    if args.run_dir:
        run_dir = Path(args.run_dir)
        found = discover_run_dir(run_dir)
        print_inventory_report(found)
        if found["ranked_peptides"] and not args.force:
            print(f"\nscoring/ranked_peptides.v03.tsv 已经存在于 {found['ranked_peptides']}，"
                  "这次没有重新打分（直接解读这份现成的表即可）。"
                  "如果确实要重新跑，加 --force。")
            return
        annotated_path = annotated_path or (str(found["annotated_wide"]) if found["annotated_wide"] else None)
        raw_events_path = raw_events_path or (str(found["raw_events"]) if found["raw_events"] else None)
        presentation_evidence_path = presentation_evidence_path or (str(found["presentation_evidence"]) if found["presentation_evidence"] else None)
        ccf_path = ccf_path or (str(found["ccf"]) if found["ccf"] else None)
        escape_flags_path = escape_flags_path or (str(found["escape_flags"]) if found["escape_flags"] else None)
        print()

    if not raw_events_path or not annotated_path:
        raise SystemExit("Need --raw-events and --annotated (or --run-dir pointing at a directory "
                          "that contains parsed/raw_events.tsv and an upstream/tools/variant_peptides*.tsv).")

    profile = load_profile(args.profile)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    events_raw = {r["event_id"]: r for r in read_tsv(raw_events_path)}
    annotated_rows = read_tsv(annotated_path)

    ccf_overrides = load_ccf_overrides(Path(ccf_path)) if ccf_path else {}
    presentation_overrides = load_presentation_evidence(Path(presentation_evidence_path)) if presentation_evidence_path else {}
    escape_overrides = load_escape_overrides(Path(escape_flags_path)) if escape_flags_path else {}
    if ccf_path:
        print(f"Using real clonality evidence from {ccf_path} (overrides raw_events.tsv's own ccf columns).")
    if presentation_evidence_path:
        print(f"Using real presentation evidence from {presentation_evidence_path} (skipping wide-table recompute where matched).")
    if escape_flags_path:
        print(f"Using real immune-escape evidence from {escape_flags_path}.")
    else:
        print("No --escape-flags / immune_escape/peptide_escape_flags.tsv found -- "
              "escape_multiplier will stay at 1.0 (not evaluated) for every peptide.")

    # --- events: enrich + safety + score, once per event_id ---
    scored_events: dict[str, dict] = {}
    for eid, ev in events_raw.items():
        e = dict(ev)
        if eid in ccf_overrides:
            ov = ccf_overrides[eid]
            for k in ("ccf_estimate", "ccf_status", "clonality_multiplier"):
                if ov.get(k) not in (None, ""):
                    e[k] = ov[k]
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
    presentation_from_sidecar = 0
    escape_applied = 0
    for row in annotated_rows:
        eid = row.get("variant_key", "")
        event = scored_events.get(eid)
        if event is None:
            skipped += 1
            continue
        peptide_id = row.get("peptide_id") or safe_id(f"{eid}_{row.get('hla_allele','')}_{row.get('mutant_peptide','')}")
        peptide = {
            "peptide_id": peptide_id,
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

        pres_key = safe_id(f"{peptide['peptide']}_{peptide['hla_allele']}")
        if pres_key in presentation_overrides:
            presentation = presentation_overrides[pres_key]
            presentation_from_sidecar += 1
        else:
            presentation = build_presentation_dict(row, profile)
        summary = {
            "mhc_i_integrity_score": event.get("appm_mhc_i_integrity", "") or "1.0",
            "mhc_ii_integrity_score": event.get("appm_mhc_ii_integrity", "") or "1.0",
            "hla_loh_alleles": "",
        }
        peptide = score_peptide(peptide, event, profile, presentation, summary)

        if peptide_id in escape_overrides:
            ov = escape_overrides[peptide_id]
            mult = to_float(ov.get("escape_multiplier"), 1.0)
            peptide["escape_multiplier"] = f"{mult:.4f}"
            severity = ov.get("escape_severity") or (
                "ESCAPE_REJECT" if mult <= 0.0 else ("ESCAPE_CAUTION" if mult < 1.0 else "ESCAPE_PASS")
            )
            peptide["escape_severity"] = severity
            peptide["escape_status"] = ov.get("escape_status", peptide.get("escape_status", ""))
            if severity == "ESCAPE_REJECT":
                peptide["safety_status"] = "FAIL"
            elif severity == "ESCAPE_CAUTION" and peptide.get("safety_status") != "FAIL":
                peptide["safety_status"] = "CAUTION"
            peptide["efficacy_score"] = f"{to_float(peptide.get('efficacy_score'), 0.0) * mult:.4f}"
            escape_applied += 1

        scored_peptides.append(peptide)

    scored_peptides.sort(key=lambda p: to_float(p.get("efficacy_score"), 0.0), reverse=True)
    for i, p in enumerate(scored_peptides, 1):
        p["pathogenic_rank"] = str(i)

    events_list = sorted(scored_events.values(), key=lambda e: to_float(e.get("event_score"), 0.0), reverse=True)

    peptide_fields = list(dict.fromkeys(k for p in scored_peptides for k in p.keys()))
    event_fields = list(dict.fromkeys(k for e in events_list for k in e.keys()))
    write_tsv(outdir / "ranked_peptides.tsv", scored_peptides, peptide_fields)
    write_tsv(outdir / "ranked_events.tsv", events_list, event_fields)

    print(f"\nScored {len(scored_peptides)} peptide x HLA rows across {len(events_list)} events.")
    if presentation_evidence_path:
        print(f"  {presentation_from_sidecar}/{len(scored_peptides)} peptides matched real presentation evidence "
              f"(rest fell back to wide-table recompute).")
    if escape_flags_path:
        print(f"  {escape_applied}/{len(scored_peptides)} peptides matched real immune-escape evidence.")
    if skipped:
        print(f"WARNING: skipped {skipped} rows whose variant_key had no matching event_id in --raw-events.")
    print(f"Wrote: {outdir / 'ranked_peptides.tsv'}")
    print(f"Wrote: {outdir / 'ranked_events.tsv'}")


if __name__ == "__main__":
    main()
