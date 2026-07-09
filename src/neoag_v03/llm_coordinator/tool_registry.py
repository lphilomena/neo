from __future__ import annotations

import json
from pathlib import Path
from typing import Any

DEFAULT_TOOL_REGISTRY: dict[str, dict[str, Any]] = {
    "neoag-sliding-run": {
        "kind": "skill",
        "intents": ["sliding_run", "workflow_run_request"],
        "required_inputs": ["somatic_vcf", "hla"],
        "optional_inputs": ["expression", "fusion", "splice", "purity", "hla_loh"],
        "outputs": ["ranked_peptides.v03.tsv", "ranked_events.v03.tsv", "evidence_report.v03.html"],
        "risk": "medium",
        "description": "Run the SNV/InDel sliding-window neoantigen workflow.",
    },
    "spechla-from-bam": {
        "kind": "external_workflow",
        "intents": ["workflow_run_request", "setup_tool", "check_status"],
        "required_inputs": ["bam"],
        "optional_inputs": ["bai", "reference_build"],
        "outputs": ["HLA typing result files", "realign BAM/VCF", "SpecHLA logs"],
        "risk": "medium",
        "description": "Extract HLA reads from an aligned BAM and run SpecHLA.",
    },
    "facets": {
        "kind": "external_workflow",
        "intents": ["workflow_run_request", "inspect_results"],
        "required_inputs": ["tumor_bam", "normal_bam", "snp_reference"],
        "outputs": ["purity", "ploidy", "CNV segments", "FACETS plots"],
        "risk": "medium",
        "description": "Tumor purity, ploidy, and allele-specific CNV analysis.",
    },
    "purple": {
        "kind": "external_workflow",
        "intents": ["workflow_run_request", "inspect_results"],
        "required_inputs": ["tumor_bam", "normal_bam", "amber", "cobalt", "reference"],
        "outputs": ["purity", "ploidy", "copy-number", "driver summary"],
        "risk": "medium",
        "description": "PURPLE/AMBER/COBALT purity and copy-number workflow.",
    },
    "project-overview": {
        "kind": "explanation",
        "intents": ["project_overview", "general_explanation"],
        "required_inputs": [],
        "outputs": ["final_response.md"],
        "risk": "low",
        "description": "Explain project scope, inputs, outputs, and workflow capabilities.",
    },
}


def load_tool_registry(project_root: str | Path = ".", skill_registry: dict[str, Any] | None = None) -> dict[str, Any]:
    root = Path(project_root)
    registry = dict(DEFAULT_TOOL_REGISTRY)
    extra = root / ".agents" / "config" / "tools_registry.json"
    if extra.exists():
        try:
            data = json.loads(extra.read_text(encoding="utf-8"))
            if isinstance(data, dict):
                registry.update(data)
        except Exception:
            pass
    for name, meta in (skill_registry or {}).items():
        item = registry.setdefault(name, {})
        item.setdefault("kind", "skill")
        item.setdefault("intents", meta.get("intents", []))
        item.setdefault("outputs", meta.get("outputs", []))
        item.setdefault("risk", "low")
        item.setdefault("description", meta.get("description", f"Registered skill {name}"))
    return registry


def registry_for_prompt(registry: dict[str, Any], max_chars: int = 24000) -> str:
    compact = {
        name: {
            "kind": meta.get("kind"),
            "intents": meta.get("intents", []),
            "required_inputs": meta.get("required_inputs", []),
            "outputs": meta.get("outputs", []),
            "risk": meta.get("risk", "low"),
            "description": meta.get("description", ""),
        }
        for name, meta in registry.items()
    }
    return json.dumps(compact, ensure_ascii=False, indent=2)[:max_chars]
