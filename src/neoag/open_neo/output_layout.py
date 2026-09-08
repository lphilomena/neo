from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Mapping

from neoag.controlled_execution.io_utils import ensure_dir, write_tsv


OUTPUT_SECTIONS = (
    ("01_tools", "Raw caller outputs and tool-specific working results."),
    ("02_evidence", "Normalized evidence, expression, safety and QC tables."),
    ("03_consensus", "Cross-tool consensus and conflict-resolution outputs."),
    ("04_ranking", "Final event and peptide rankings plus score tables."),
    ("05_reports", "Patient, technical and review reports."),
    ("06_logs", "Execution logs, audit trails and issue records."),
    ("07_manifests", "Run manifests, input inventories and reproducibility metadata."),
)


def _safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_.") or "artifact"


def _section_for_artifact(name: str, path: Path) -> str:
    text = f"{name} {path}".lower()
    if "report" in text or path.suffix.lower() in {".html", ".pdf", ".docx", ".pptx"}:
        return "05_reports"
    if any(token in text for token in ("log", "audit", "issue")):
        return "06_logs"
    if any(token in text for token in ("manifest", "route", "inventory", "config", "plan")):
        return "07_manifests"
    if any(token in text for token in ("rank", "score", "weighted_baseline", "all_tool_results")):
        return "04_ranking"
    # The final consensus peptide/event tables are the published rankings;
    # reserve the consensus section for tool-level reconciliation outputs.
    if "consensus_peptide" in text or "consensus_event" in text:
        return "04_ranking"
    if any(token in text for token in ("consensus", "conflict", "hla_loh")):
        return "03_consensus"
    if any(token in text for token in ("tool", "caller", "branch")):
        return "01_tools"
    return "02_evidence"


def _link(destination: Path, target: Path) -> str:
    """Create/refresh only a link owned by this view; never move raw output."""
    if destination.is_symlink():
        try:
            if destination.resolve() == target.resolve():
                return "REUSED"
        except OSError:
            pass
        destination.unlink()
    elif destination.exists():
        return "PRESERVED_EXISTING"
    destination.symlink_to(target, target_is_directory=target.is_dir())
    return "LINKED"


def materialize_output_view(
    output_root: str | Path,
    *,
    source_root: str | Path | None = None,
    artifacts: Mapping[str, Any] | None = None,
    layout_paths: Mapping[str, str | Path] | None = None,
    producer: str = "open-neo",
) -> dict[str, str]:
    """Build a canonical, non-destructive navigational view of a run's output.

    The production runner and macro skills retain their native trees for safe
    resume.  This function adds only symlinks and a TSV index under
    ``deliverables/`` so all entrypoints expose the same user-facing layout.
    """
    output_root = Path(output_root).resolve()
    source_root = Path(source_root or output_root).resolve()
    view = ensure_dir(output_root / "deliverables")
    for section, _ in OUTPUT_SECTIONS:
        ensure_dir(view / section)

    rows: list[dict[str, str]] = []

    def add(section: str, label: str, value: str | Path | None) -> None:
        if not value:
            return
        target = Path(value)
        if not target.exists():
            return
        target = target.resolve()
        destination = view / section / _safe_name(label)
        status = _link(destination, target)
        rows.append({
            "section": section,
            "label": label,
            "source_path": str(target),
            "view_path": str(destination),
            "kind": "directory" if target.is_dir() else "file",
            "status": status,
        })

    native_groups = {
        "01_tools": ("tools", "upstream/tools", "pipeline/production/branches", "production/branches"),
        "02_evidence": ("parsed", "upstream/parsed", "presentation", "safety", "appm", "clonality", "immune_escape", "evidence"),
        "03_consensus": ("tool_consensus", "hla_loh", "hla_loh_consensus", "consensus"),
        "04_ranking": ("scoring", "ranking"),
        "05_reports": ("reports",),
        "06_logs": ("logs",),
        "07_manifests": ("manifest", "manifests"),
    }
    for section, relatives in native_groups.items():
        for relative in relatives:
            add(section, relative.replace("/", "_"), source_root / relative)

    for name, value in (layout_paths or {}).items():
        add(_section_for_artifact(name, Path(value)), name, value)
    for name, value in (artifacts or {}).items():
        if isinstance(value, (str, Path)):
            add(_section_for_artifact(str(name), Path(value)), str(name), value)

    index = view / "output_index.tsv"
    write_tsv(index, rows)
    readme = view / "README.md"
    section_lines = "\n".join(f"- `{section}/`: {description}" for section, description in OUTPUT_SECTIONS)
    readme.write_text(
        "# Open-Neo deliverables\n\n"
        f"This view was generated by `{producer}`. It is an index of the original results, not a copy. "
        "All tool outputs remain in their native locations so resume and provenance remain valid.\n\n"
        f"{section_lines}\n\n"
        "`output_index.tsv` records every linked item and its original path.\n",
        encoding="utf-8",
    )
    return {
        "deliverables_dir": str(view),
        "deliverables_index": str(index),
        "deliverables_readme": str(readme),
    }
