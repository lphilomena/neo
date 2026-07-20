from __future__ import annotations

import argparse
from pathlib import Path
from .common import ensure_dir, executable_status, run_command, write_json, write_tsv

DEFAULT_TOOLS = [
    "neoag", "nextflow", "java", "gatk", "vep", "netMHCpan", "mhcflurry-predict",
    "netMHCstabpan", "PRIME", "bigmhc_predict", "deepimmuno-cnn.py", "LOHHLA",
    "runFACETS.R", "ascat.R", "star-fusion-neoag", "arriba", "fusioncatcher-neoag", "easyfuse-neoag", "pyclone",
]

REF_KEYS = [
    "genome_fasta", "genome_fai", "gencode_gtf", "vep_cache_dir", "normal_expression",
    "normal_hla_ligands", "reference_proteome", "normal_junctions", "population_sv_vcf",
]


def parse_simple_manifest(path: Path) -> dict[str, str]:
    vals: dict[str, str] = {}
    if not path.exists():
        return vals
    for line in path.read_text(encoding="utf-8").splitlines():
        s = line.strip()
        if not s or s.startswith("#") or ":" not in s:
            continue
        k, v = s.split(":", 1)
        vals[k.strip()] = v.strip().strip('"\'')
    return vals


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Project B P0 skill: tool and reference QC")
    ap.add_argument("--project-root", default=".")
    ap.add_argument("--reference-manifest")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--mini-smoke", action="store_true", help="Run lightweight command smoke tests when safe")
    args = ap.parse_args(argv)
    outdir = ensure_dir(args.outdir)
    root = Path(args.project_root).resolve()

    tool_rows = []
    for exe in DEFAULT_TOOLS:
        row = executable_status(exe)
        tool_rows.append(row)
    write_tsv(outdir / "tool_qc_report.tsv", tool_rows, ["tool", "status", "path", "message"])

    ref_rows = []
    if args.reference_manifest:
        refs = parse_simple_manifest(Path(args.reference_manifest))
        for key in REF_KEYS:
            val = refs.get(key, "")
            exists = bool(val) and Path(val).exists()
            ref_rows.append({"reference": key, "path": val, "status": "OK" if exists else "MISSING", "message": "exists" if exists else "not found or unset"})
    write_tsv(outdir / "reference_qc_report.tsv", ref_rows, ["reference", "path", "status", "message"])

    smoke = []
    if args.mini_smoke:
        commands = [
            ["neoag", "--help"],
            ["vep", "--help"],
            ["mhcflurry-predict", "--help"],
            ["gatk", "--version"],
            ["netMHCpan", "-h"],
        ]
        for cmd in commands:
            if executable_status(cmd[0])["status"] == "OK":
                smoke.append(run_command(cmd, cwd=root, timeout=30))
            else:
                smoke.append({"cmd": " ".join(cmd), "ok": False, "returncode": 127, "stdout": "", "stderr": "entrypoint not found"})
    write_json(outdir / "tool_smoke_report.json", smoke)

    missing = [r for r in tool_rows if r["status"] != "OK"]
    md = ["# Tool and reference QC", "", f"Project root: `{root}`", "", "## Tool status"]
    for r in tool_rows:
        md.append(f"- {r['tool']}: **{r['status']}** {r['message']} {r['path']}")
    if ref_rows:
        md.append("\n## Reference status")
        for r in ref_rows:
            md.append(f"- {r['reference']}: **{r['status']}** `{r['path']}`")
    if missing:
        md.append("\n## Notes")
        md.append("Missing tools should be interpreted as missing evidence capability, not negative biological evidence.")
    (outdir / "tool_reference_qc.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
