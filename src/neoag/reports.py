"""Legacy v0.3 report entrypoint — delegates to dual reports."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from .reports_dual import ReportBundle, load_report_bundle, make_technical_report


def make_report(
    path,
    profile: Mapping[str, Any],
    events,
    peptides,
    appm_summary,
    validation_rows,
    *,
    outdir: str | Path | None = None,
    provenance: Mapping[str, Any] | None = None,
    sample_id: str = "",
    entry_mode: str = "",
):
    """Backward-compatible single technical report writer."""
    bundle = load_report_bundle(
        profile=profile,
        events=events,
        peptides=peptides,
        appm_summary=appm_summary,
        validation_rows=validation_rows,
        outdir=outdir,
        provenance=provenance,
        sample_id=sample_id,
        entry_mode=entry_mode,
    )
    make_technical_report(path, bundle)
