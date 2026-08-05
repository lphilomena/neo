from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

from .contracts import InstallCheckInput, ReviewInput, RunInput, validate_json_schema
from .errors import exit_code_for_result
from .install_check import run_install_check
from .run import run_open_neo
from .review import run_review


def _add_common_output(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--case-id")


def _json_list(value: str) -> list[str]:
    try:
        loaded = json.loads(value)
        if isinstance(loaded, list):
            return [str(x) for x in loaded]
    except Exception:
        pass
    return [x for x in value.replace(";", ",").split(",") if x]


def _tool_result(value: str) -> tuple[str, str, str]:
    if "=" not in value or ":" not in value.split("=", 1)[0]:
        raise argparse.ArgumentTypeError("tool result must use DOMAIN:TOOL=PATH")
    left, path = value.split("=", 1)
    domain, tool = left.split(":", 1)
    return domain.strip(), tool.strip(), path.strip()


def _claude_code_channel(value: str) -> str:
    if value in {"stable", "latest"} or re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+", value):
        return value
    raise argparse.ArgumentTypeError("Claude Code channel must be stable, latest, or X.Y.Z")


def _validate_public_input(command: str, args: dict[str, Any]) -> list[str]:
    skill = {"install-check": "open-neo-install-check", "run": "open-neo-run", "review": "open-neo-review"}[command]
    contract = {"install-check": InstallCheckInput, "run": RunInput, "review": ReviewInput}[command]
    normalized = contract.from_mapping(args).to_mapping()
    root = Path(__file__).resolve().parents[3]
    schema_path = root / ".agents" / "skills" / skill / "references" / "INPUT_SCHEMA.json"
    if not schema_path.is_file():
        return []
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    return validate_json_schema(normalized, schema)


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="Open-Neo public macro Skills CLI")
    sub = ap.add_subparsers(dest="command", required=True)

    install = sub.add_parser("install-check", help="Install/environment/reference check for a new machine")
    _add_common_output(install)
    install.add_argument("--project-root", default=".")
    install.add_argument("--release-tarball")
    install.add_argument("--sha256")
    install.add_argument("--deployment-tier", choices=["review", "core", "prediction", "full"], default="core")
    install.add_argument("--mode", choices=["plan", "verify", "repair", "install", "resume"], default="verify")
    install.add_argument("--tools-manifest")
    install.add_argument("--reference-manifest")
    install.add_argument("--sample-manifest")
    install.add_argument(
        "--conda-base",
        help="Existing Conda/Miniforge root; use this NAS/site installation instead of auto-installing Miniforge",
    )
    install.add_argument("--profile", default="local")
    install.add_argument("--run-demo", action="store_true")
    install.add_argument("--run-pytest", action="store_true")
    install.add_argument("--run-nextflow", action="store_true")
    install.add_argument("--mini-smoke", action="store_true")
    install.add_argument("--no-release-audit", action="store_true")
    install.add_argument("--approved", action="store_true")
    install.add_argument("--deploy-root", default="/opt/neoag")
    install.add_argument("--tools-root")
    install.add_argument("--reference-root")
    install.add_argument("--licensed-root")
    install.add_argument("--asset-source-host")
    install.add_argument("--asset-source-root")
    install.add_argument("--asset-manifest")
    install.add_argument("--deployment-reference-manifest")
    install.add_argument("--install-timeout", type=int, default=7200)
    install.add_argument("--allow-download", action="store_true")
    install.add_argument(
        "--install-claude-code", action="store_true",
        help="Install Claude Code; requires approved mutating mode and --allow-download",
    )
    install.add_argument(
        "--claude-code-channel", type=_claude_code_channel, default="stable",
        help="Claude Code release channel or exact version (default: stable)",
    )
    install.add_argument(
        "--installer-profile", choices=["minimal", "standard", "all-open", "all"], default=None,
        help="Installer scope; defaults to minimal for review/core and standard for prediction/full",
    )
    install.add_argument("--no-sync-assets", action="store_true")

    run = sub.add_parser("run", help="Detect inputs, route, run the pipeline, and emit weighted plus evidence-consensus rankings")
    _add_common_output(run)
    run.add_argument("--project-root", default=".")
    run.add_argument("--sample-manifest")
    run.add_argument("--mode", choices=["plan", "dry-run", "execute", "resume", "ranking-only"])
    run.add_argument("--approved", action="store_true")
    run.add_argument("--allow-partial", action="store_true")
    run.add_argument("--doctor", action=argparse.BooleanOptionalAction, default=True)
    run.add_argument("--tools-manifest")
    run.add_argument("--reference-manifest")
    run.add_argument("--execution-profile", default="local")
    run.add_argument("--automatic-tool-policy", choices=["all-available", "balanced", "minimal"], default="all-available")
    run.add_argument("--mini-smoke", action="store_true")
    run.add_argument("--release-audit", action="store_true")
    run.add_argument("--stub", action="store_true")
    run.add_argument("--profile", default="default")
    run.add_argument("--genome-build", default="GRCh38")
    run.add_argument("--sample-id")
    run.add_argument("--input-dir")
    run.add_argument("--tumor-dna-bam")
    run.add_argument("--normal-dna-bam")
    run.add_argument("--tumor-rna-bam")
    run.add_argument("--tumor-dna-fastq", action="append")
    run.add_argument("--normal-dna-fastq", action="append")
    run.add_argument("--tumor-rna-fastq", action="append")
    run.add_argument("--tumor-sample-id")
    run.add_argument("--normal-sample-id")
    run.add_argument("--assay-type", choices=["WGS", "WES", "PANEL", "CAPTURE", "RNA"])
    run.add_argument("--somatic-vcf")
    run.add_argument("--fusion-tsv")
    run.add_argument("--splice-junction-tsv")
    run.add_argument("--sv-vcf", action="append")
    run.add_argument("--capture-bed")
    run.add_argument("--peptide-csv")
    run.add_argument("--raw-events")
    run.add_argument("--raw-peptides")
    run.add_argument("--sv-raw-events")
    run.add_argument("--sv-raw-peptides")
    run.add_argument("--hla-file")
    run.add_argument("--hla-alleles", type=_json_list)
    run.add_argument("--expression-tsv")
    run.add_argument("--transcript-expression-tsv")
    run.add_argument("--rna-evidence-tsv")
    run.add_argument("--rna-quant-method", choices=["auto", "salmon", "rsem"], default="auto")
    run.add_argument("--salmon-index")
    run.add_argument("--tx2gene")
    run.add_argument("--rsem-reference")
    run.add_argument("--star-index")
    run.add_argument("--ctat-genome-lib")
    run.add_argument("--easyfuse-ref")
    run.add_argument("--normal-readthrough")
    run.add_argument("--snaf-workflow")
    run.add_argument("--snaf-db")
    run.add_argument("--snaf-python")
    run.add_argument("--altanalyze-image")
    run.add_argument("--splicemutr-workflow")
    run.add_argument("--rna-threads", type=int, default=16)
    run.add_argument("--purity-tsv")
    run.add_argument("--cnv-tsv")
    run.add_argument("--hla-loh-tsv")
    run.add_argument("--normal-expression")
    run.add_argument("--normal-hla-ligands")
    run.add_argument("--reference-proteome")
    run.add_argument("--normal-junctions")
    run.add_argument("--reference-fasta")
    run.add_argument("--gencode-gtf")
    run.add_argument("--vep-cache")
    run.add_argument("--production-manifest")
    run.add_argument("--result-dir")
    run.add_argument("--comprehensive-evidence")
    run.add_argument("--weighted-baseline")
    run.add_argument("--rules")
    run.add_argument("--provenance")
    run.add_argument("--force", action="store_true")
    run.add_argument("--timeout", type=int, default=7200)
    run.add_argument("--tool-result", action="append", type=_tool_result, default=[], metavar="DOMAIN:TOOL=PATH")
    run.add_argument("--gateway-url", help="Required for execute/resume when not already invoked by NeoAg Gateway")
    run.add_argument("--gateway-wait", action="store_true", help="Poll Gateway until the submitted job finishes")

    review = sub.add_parser("review", help="Review event-level consensus, design experiments, and generate reports")
    _add_common_output(review)
    review.add_argument("--result-dir", required=True)
    review.add_argument("--top-n", type=int, default=12)
    review.add_argument("--clinical-context")
    review.add_argument("--disease-profile")
    review.add_argument("--therapy-context", default="research")
    review.add_argument(
        "--reports", type=_json_list, default=["patient", "technical", "onepage"],
        help="Comma-separated report outputs: patient,technical,onepage; use none for review tables only",
    )
    return ap


def main(argv: list[str] | None = None) -> int:
    args = vars(build_parser().parse_args(argv))
    command = args.pop("command")
    if command == "install-check":
        args["release_audit"] = not args.pop("no_release_audit", False)
    elif command == "run":
        tool_results: dict[str, dict[str, str]] = {}
        for domain, tool, path in args.pop("tool_result", []):
            tool_results.setdefault(domain, {})[tool] = path
        if tool_results:
            args["tool_results"] = tool_results
    validation_errors = _validate_public_input(command, args)
    if validation_errors:
        result = {
            "schema_version": "open-neo-macro-skill-v1",
            "skill": f"open-neo-{command}",
            "status": "BLOCKED",
            "blocking_issues": validation_errors,
            "steps": [],
            "outputs": {},
        }
    elif command == "install-check":
        result = run_install_check(args)
    elif command == "run":
        result = run_open_neo(args)
    else:
        result = run_review(args)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return exit_code_for_result(result)


if __name__ == "__main__":
    raise SystemExit(main())
