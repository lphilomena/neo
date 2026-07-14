from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from .io import ensure_dir, write_json, write_tsv
from .registry import CATEGORY_LABELS, SKILLS_BY_NAME, list_specs, registry_dict
from .runner import run_skill, validate_skill_dirs


def _parse_kv(items: list[str]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for item in items:
        if "=" not in item:
            raise SystemExit(f"--arg requires key=value, got {item!r}")
        k, v = item.split("=", 1)
        # basic JSON parsing for booleans/numbers/lists
        try:
            out[k] = json.loads(v)
        except Exception:
            out[k] = v
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="NeoAg A/B/C/D skill taxonomy CLI")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list", help="List all skills")
    desc = sub.add_parser("describe", help="Describe one skill")
    desc.add_argument("skill")
    exp = sub.add_parser("export-registry", help="Export registry JSON")
    exp.add_argument("--out", required=True)
    val = sub.add_parser("validate", help="Validate .agents/skills directories")
    val.add_argument("--root", default=".")
    val.add_argument("--outdir", default="work/skill_taxonomy_validate")
    run = sub.add_parser("run", help="Run a skill handler")
    run.add_argument("skill")
    run.add_argument("--outdir", required=True)
    run.add_argument("--dry-run", action="store_true")
    run.add_argument("--arg", action="append", default=[], help="key=value arguments passed to skill")
    args = ap.parse_args(argv)

    if args.cmd == "list":
        print("# NeoAg skills")
        for cat, label in CATEGORY_LABELS.items():
            print(f"\n[{cat}] {label}")
            for spec in list_specs(cat):
                print(f"- {spec.name}: {spec.purpose} (risk={spec.risk_level})")
        return 0
    if args.cmd == "describe":
        spec = SKILLS_BY_NAME[args.skill]
        print(json.dumps(spec.to_dict(), ensure_ascii=False, indent=2))
        return 0
    if args.cmd == "export-registry":
        write_json(args.out, registry_dict())
        print(args.out)
        return 0
    if args.cmd == "validate":
        outdir = ensure_dir(args.outdir)
        res = validate_skill_dirs(args.root)
        write_json(outdir / "skill_taxonomy_validation.json", res)
        write_tsv(outdir / "skill_taxonomy_validation.tsv", res["rows"])
        print(json.dumps({"status": res["status"], "missing": res["missing"]}, ensure_ascii=False))
        return 0 if res["status"] == "PASS" else 1
    if args.cmd == "run":
        kv = _parse_kv(args.arg)
        kv["outdir"] = args.outdir
        res = run_skill(args.skill, kv, dry_run=args.dry_run)
        print(json.dumps(res, ensure_ascii=False, indent=2))
        return 0 if res.get("status") not in {"FAIL"} else 2
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
