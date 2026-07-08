from __future__ import annotations

INTENT_TO_SKILLS: dict[str, list[str]] = {
    "input_check": ["neoag-input-qc"],
    "tool_check": ["neoag-tool-and-reference-qc"],
    "demo_smoke": ["neoag-run-demo-and-smoke"],
    "run_scoring": ["neoag-input-qc", "neoag-evidence-scoring"],
    "appm_escape_review": ["neoag-input-qc", "neoag-hla-loh-appm-review"],
    "ccf_review": ["neoag-input-qc", "neoag-ccf-clonality-review"],
    "ranking_compare": ["neoag-input-qc", "neoag-ranking-compare"],
    "patient_report_update": ["neoag-input-qc", "neoag-ranking-compare", "neoag-hla-loh-appm-review", "neoag-ccf-clonality-review", "neoag-patient-report"],
    "experiment_design": ["neoag-input-qc", "neoag-ranking-compare", "neoag-hla-loh-appm-review", "neoag-ccf-clonality-review"],
    "general_explanation": [],
}

INTENT_DESCRIPTIONS = {
    "input_check": "Inspect provided sample files/result directory and recommend the appropriate Project B workflow.",
    "tool_check": "Check external tool and reference readiness.",
    "demo_smoke": "Run or plan Project B pytest/run-demo/Nextflow smoke tests.",
    "run_scoring": "Plan or run evidence scoring from raw_events/raw_peptides/presentation inputs.",
    "appm_escape_review": "Review APPM, HLA LOH, and immune escape evidence.",
    "ccf_review": "Review CCF, purity, clonality, and CCF modifiers.",
    "ranking_compare": "Compare NetMHCpan42 ranking with recommendation ranking.",
    "patient_report_update": "Generate or update the patient-facing report.",
    "experiment_design": "Create an experimental validation candidate grouping plan.",
    "general_explanation": "Explain concepts without executing workflows.",
}
