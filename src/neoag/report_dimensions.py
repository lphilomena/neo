from __future__ import annotations

"""Machine-readable alignment between report language and ranking fields.

The report may describe more evidence than the automatic EC ranking uses. This
module makes that boundary explicit and auditable:

* ``automatic``: directly used by hard-fail/cap/R/Pareto logic.
* ``track_specific``: used only for selected event classes.
* ``review_overlay``: reviewed after pipeline ranking; must not silently change it.
* ``post_ranking_validation``: experimental outcomes, not an input to EC-v2.x.
"""

from pathlib import Path
import tomllib
from typing import Any, Mapping

from .utils import read_tsv, write_tsv


DEFAULT_MAP = Path(__file__).resolve().parents[2] / "configs" / "report" / "report_dimension_map_v1.toml"


def load_report_dimension_map(path: str | Path | None = None) -> dict[str, Any]:
    target = Path(path) if path else DEFAULT_MAP
    if not target.is_file():
        raise FileNotFoundError(f"Report dimension map not found: {target}")
    with target.open("rb") as handle:
        payload = tomllib.load(handle)
    payload["_path"] = str(target)
    return payload


def _as_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return [str(item) for item in value]
    return [str(value)]


def audit_report_dimensions(
    input_tsv: str | Path,
    output_tsv: str | Path,
    map_path: str | Path | None = None,
) -> dict[str, Any]:
    """Audit whether report-declared dimensions have corresponding data fields.

    This is a schema/traceability audit only. It does not reinterpret biological
    calls and does not alter R grades, Pareto fronts, or event ranks.
    """

    rows = read_tsv(input_tsv)
    available = set(rows[0]) if rows else set()
    mapping = load_report_dimension_map(map_path)
    dimensions = mapping.get("dimensions", {})
    output: list[dict[str, str]] = []
    blocking = 0
    for key, spec_any in dimensions.items():
        spec = spec_any if isinstance(spec_any, Mapping) else {}
        fields = _as_list(spec.get("fields"))
        present = [field for field in fields if field in available]
        missing = [field for field in fields if field not in available]
        layer = str(spec.get("layer", ""))
        required = bool(spec.get("required_for_report", False))
        status = "PASS"
        if not fields:
            status = "INFO"
        elif not present:
            status = "MISSING"
        elif missing:
            status = "PARTIAL"
        if required and status in {"MISSING", "PARTIAL"}:
            blocking += 1
        output.append({
            "dimension_key": str(key),
            "label_zh": str(spec.get("label_zh", key)),
            "layer": layer,
            "ranking_role": str(spec.get("ranking_role", "")),
            "applies_to": ",".join(_as_list(spec.get("applies_to"))),
            "declared_fields": ",".join(fields),
            "present_fields": ",".join(present),
            "missing_fields": ",".join(missing),
            "required_for_report": "yes" if required else "no",
            "status": status,
            "interpretation_boundary": str(spec.get("interpretation_boundary", "")),
        })
    write_tsv(output_tsv, output)
    return {
        "input": str(input_tsv),
        "output": str(output_tsv),
        "map": str(mapping.get("_path", map_path or DEFAULT_MAP)),
        "dimensions": len(output),
        "blocking_dimensions": blocking,
        "status": "PASS" if blocking == 0 else "PARTIAL",
    }


__all__ = ["DEFAULT_MAP", "load_report_dimension_map", "audit_report_dimensions"]
