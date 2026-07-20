from __future__ import annotations

import argparse
import subprocess
from pathlib import Path
from .common import ensure_dir


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Plan or run Project B evidence scoring")
    ap.add_argument("--raw-events", required=True)
    ap.add_argument("--raw-peptides", required=True)
    ap.add_argument("--presentation", required=True)
    ap.add_argument("--appm-summary")
    ap.add_argument("--appm-peptide-modifiers")
    ap.add_argument("--ccf")
    ap.add_argument("--peptide-safety")
    ap.add_argument("--peptide-escape-flags")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--execute", action="store_true")
    args = ap.parse_args(argv)
    outdir = ensure_dir(args.outdir)
    score_dir = ensure_dir(outdir / "scoring")
    out_events = score_dir / "ranked_events.tsv"
    out_peptides = score_dir / "ranked_peptides.tsv"
    cmd = ["neoag", "score", "--raw-events", args.raw_events, "--raw-peptides", args.raw_peptides, "--presentation", args.presentation, "--out-events", str(out_events), "--out-peptides", str(out_peptides)]
    optional = {
        "--appm-summary": args.appm_summary,
        "--appm-peptide-modifiers": args.appm_peptide_modifiers,
        "--ccf": args.ccf,
        "--peptide-safety": args.peptide_safety,
        "--peptide-escape-flags": args.peptide_escape_flags,
    }
    for flag, val in optional.items():
        if val:
            cmd += [flag, val]
    plan = "# Evidence scoring plan\n\n```bash\n" + " ".join(cmd) + "\n```\n"
    (outdir / "evidence_scoring_plan.md").write_text(plan, encoding="utf-8")
    if not args.execute:
        return 0
    proc = subprocess.run(cmd, text=True, capture_output=True)
    (outdir / "evidence_scoring.stdout.txt").write_text(proc.stdout, encoding="utf-8")
    (outdir / "evidence_scoring.stderr.txt").write_text(proc.stderr, encoding="utf-8")
    return proc.returncode

if __name__ == "__main__":
    raise SystemExit(main())
