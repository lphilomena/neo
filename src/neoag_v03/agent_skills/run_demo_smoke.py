from __future__ import annotations

import argparse
import subprocess
from pathlib import Path
from .common import ensure_dir


def run(cmd: list[str], cwd: Path, out_prefix: Path, execute: bool) -> int:
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    (out_prefix.with_suffix(".cmd.txt")).write_text(" ".join(cmd) + "\n", encoding="utf-8")
    if not execute:
        return 0
    proc = subprocess.run(cmd, cwd=str(cwd), text=True, capture_output=True)
    (out_prefix.with_suffix(".stdout.txt")).write_text(proc.stdout, encoding="utf-8")
    (out_prefix.with_suffix(".stderr.txt")).write_text(proc.stderr, encoding="utf-8")
    return proc.returncode


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Run Project B demo and smoke checks")
    ap.add_argument("--project-root", default=".")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--execute", action="store_true")
    ap.add_argument("--skip-pytest", action="store_true")
    ap.add_argument("--skip-demo", action="store_true")
    ap.add_argument("--run-nextflow", action="store_true")
    args = ap.parse_args(argv)
    root = Path(args.project_root).resolve()
    outdir = ensure_dir(args.outdir)
    steps = []
    rc = 0
    if not args.skip_pytest:
        cmd = ["python", "-m", "pytest", "-q"]
        steps.append("pytest -q")
        rc = rc or run(cmd, root, outdir / "pytest", args.execute)
    if not args.skip_demo:
        demo_out = outdir / "demo_v043"
        cmd = ["python", "-m", "neoag_v03.cli", "run-demo", "--outdir", str(demo_out), "--sample-id", "DEMO001"]
        steps.append("run-demo")
        rc = rc or run(cmd, root, outdir / "run_demo", args.execute)
    if args.run_nextflow:
        cmd = ["bin/neoag-nextflow", "run", "workflows/main.nf", "--pvac_files", "data/fixtures/pvacseq_aggregated.tsv", "--outdir", str(outdir / "demo_nf")]
        steps.append("nextflow fixture")
        rc = rc or run(cmd, root, outdir / "nextflow", args.execute)
    md = ["# Demo and smoke plan", "", f"Project root: `{root}`", f"Execute: {args.execute}", "", "## Steps"] + [f"- {s}" for s in steps]
    (outdir / "smoke_test_report.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    return rc

if __name__ == "__main__":
    raise SystemExit(main())
