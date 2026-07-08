"""Evidence provenance metadata for standard intermediate TSVs (spec §6.3)."""

from __future__ import annotations

import os
import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

from .schemas import EVIDENCE_PROVENANCE_FIELDS, TOOL_PROVENANCE_TOOLS

# Generic provenance columns on every evidence table row.
STATUS_REAL = "real"
STATUS_STUB = "stub"
STATUS_MISSING = "missing"
STATUS_NOT_USED = "not_used"
STATUS_INVALID_FOR_PRODUCTION = "invalid_for_production"

SOURCE_LOCAL = "local"
SOURCE_FIXTURE = "fixture"
SOURCE_EXTERNAL = "external"
SOURCE_DERIVED = "derived"
SOURCE_MISSING = "missing"

MODE_TOOL_RUN = "tool_run"
MODE_PASSTHROUGH = "passthrough"
MODE_CONVERTED = "converted"
MODE_DERIVED = "derived"
MODE_STUB = "stub"

IMMUNOGENICITY_PRODUCTION_INVALID = frozenset({"prime", "bigmhc_im", "deepimmuno", "iedb"})


@dataclass
class ProvenanceRecord:
    tool: str
    source: str = SOURCE_LOCAL
    status: str = STATUS_REAL
    mode: str = MODE_TOOL_RUN
    file: str = ""
    version: str = ""

    def as_fields(self) -> dict[str, str]:
        return {
            "evidence_source": self.source,
            "evidence_tool": self.tool,
            "evidence_tool_version": self.version or detect_tool_version(self.tool),
            "evidence_mode": self.mode,
            "evidence_file": self.file,
            "evidence_status": self.status,
        }

    def tool_summary(self) -> dict[str, str]:
        ver = self.version or detect_tool_version(self.tool)
        return {
            f"{self.tool}_source": self.source,
            f"{self.tool}_version": ver,
            f"{self.tool}_status": self.status,
        }


def detect_tool_version(tool: str) -> str:
    """Best-effort version string for registered upstream tools."""
    from .tools.registry import ROOT, TOOL_REGISTRY

    if tool == "netmhcpan":
        ver_file = ROOT / "tools" / "netMHCpan" / "data" / "version"
        if ver_file.is_file():
            return ver_file.read_text(encoding="utf-8", errors="ignore").strip()
    if tool == "netmhcstabpan":
        ver_file = ROOT / "tools" / "netMHCstabpan" / "data" / "version"
        if ver_file.is_file():
            return ver_file.read_text(encoding="utf-8", errors="ignore").strip()
    if tool == "mhcflurry":
        try:
            out = subprocess.run(
                ["python3", "-c", "import mhcflurry; print(mhcflurry.__version__)"],
                capture_output=True,
                text=True,
                timeout=15,
                check=False,
            )
            if out.returncode == 0 and out.stdout.strip():
                return out.stdout.strip()
        except (OSError, subprocess.TimeoutExpired):
            pass
    if tool == "facets":
        return "mskcc-facets"
    if tool == "lohhla":
        return "LOHHLA"
    if tool == "vep":
        return os.environ.get("NEOAG_VEP_ENV", "vep")
    spec = TOOL_REGISTRY.get(tool)
    if spec:
        return spec.executable
    return "unknown"


def provenance_for_netmhcpan_run(
    path: str | Path,
    backend: str,
    *,
    stub: bool = False,
) -> ProvenanceRecord:
    """Provenance for NetMHCpan outputs produced by neoag tool runner."""
    if stub:
        return provenance_stub("netmhcpan")
    if backend in {"iedb", "api", "remote"}:
        return ProvenanceRecord(
            tool="netmhcpan",
            source=SOURCE_EXTERNAL,
            status=STATUS_REAL,
            mode=MODE_TOOL_RUN,
            file=str(path),
            version="iedb-netmhcpan",
        )
    return ProvenanceRecord(
        tool="netmhcpan",
        source=SOURCE_LOCAL,
        status=STATUS_REAL,
        mode=MODE_TOOL_RUN,
        file=str(path),
        version=detect_tool_version("netmhcpan"),
    )


def provenance_from_file(
    tool: str,
    path: str | Path | None,
    *,
    mode: str = MODE_PASSTHROUGH,
    source: str = SOURCE_LOCAL,
    status: str = STATUS_REAL,
    version: str = "",
) -> ProvenanceRecord:
    file_str = str(path) if path else ""
    if file_str and "fixtures" in file_str.replace("\\", "/"):
        source = SOURCE_FIXTURE
    return ProvenanceRecord(
        tool=tool,
        source=source,
        status=status,
        mode=mode,
        file=file_str,
        version=version,
    )


def provenance_stub(tool: str, *, production_invalid: bool = False) -> ProvenanceRecord:
    status = STATUS_INVALID_FOR_PRODUCTION if production_invalid else STATUS_STUB
    return ProvenanceRecord(
        tool=tool,
        source=SOURCE_FIXTURE,
        status=status,
        mode=MODE_STUB,
        file="stub",
    )


def provenance_not_used(tool: str) -> ProvenanceRecord:
    return ProvenanceRecord(
        tool=tool,
        source=SOURCE_MISSING,
        status=STATUS_NOT_USED,
        mode=MODE_DERIVED,
        file="",
    )


def provenance_missing(tool: str) -> ProvenanceRecord:
    return ProvenanceRecord(
        tool=tool,
        source=SOURCE_MISSING,
        status=STATUS_MISSING,
        mode=MODE_DERIVED,
        file="",
    )


def provenance_derived(tool: str, path: str | Path | None, *, upstream: str = "") -> ProvenanceRecord:
    return ProvenanceRecord(
        tool=tool,
        source=SOURCE_DERIVED,
        status=STATUS_REAL,
        mode=MODE_DERIVED,
        file=str(path) if path else upstream,
    )


def attach_provenance(
    rows: list[dict[str, str]],
    record: ProvenanceRecord | Mapping[str, str],
) -> list[dict[str, str]]:
    if isinstance(record, ProvenanceRecord):
        fields = record.as_fields()
    else:
        fields = dict(record)
    if not rows:
        return [fields]
    return [{**row, **fields} for row in rows]


def write_evidence_tsv(
    path: str | Path,
    rows: list[dict[str, str]],
    data_fields: list[str],
    provenance: ProvenanceRecord | Mapping[str, str],
) -> None:
    from .utils import write_tsv

    fields = list(data_fields) + EVIDENCE_PROVENANCE_FIELDS
    write_tsv(path, attach_provenance(rows, provenance), fields)


@dataclass
class ProvenanceRegistry:
    """Track per-tool provenance for composite tables and provenance.v03.json."""

    _records: dict[str, ProvenanceRecord] = field(default_factory=dict)

    def set(self, record: ProvenanceRecord) -> None:
        self._records[record.tool] = record

    def register_real(
        self,
        tool: str,
        path: str | Path | None,
        *,
        mode: str = MODE_TOOL_RUN,
        source: str = SOURCE_LOCAL,
        version: str = "",
    ) -> ProvenanceRecord:
        rec = provenance_from_file(tool, path, mode=mode, source=source, status=STATUS_REAL, version=version)
        self.set(rec)
        return rec

    def register_passthrough(self, tool: str, path: str | Path | None) -> ProvenanceRecord:
        rec = provenance_from_file(tool, path, mode=MODE_PASSTHROUGH)
        self.set(rec)
        return rec

    def register_converted(self, tool: str, path: str | Path | None) -> ProvenanceRecord:
        rec = provenance_from_file(tool, path, mode=MODE_CONVERTED)
        self.set(rec)
        return rec

    def register_stub(self, tool: str, *, production_invalid: bool | None = None) -> ProvenanceRecord:
        invalid = (
            production_invalid
            if production_invalid is not None
            else tool in IMMUNOGENICITY_PRODUCTION_INVALID
        )
        rec = provenance_stub(tool, production_invalid=invalid)
        self.set(rec)
        return rec

    def register_not_used(self, tool: str) -> ProvenanceRecord:
        rec = provenance_not_used(tool)
        self.set(rec)
        return rec

    def register_missing(self, tool: str) -> ProvenanceRecord:
        rec = provenance_missing(tool)
        self.set(rec)
        return rec

    def register_derived(self, tool: str, path: str | Path | None, *, upstream: str = "") -> ProvenanceRecord:
        rec = provenance_derived(tool, path, upstream=upstream)
        self.set(rec)
        return rec

    def get(self, tool: str) -> ProvenanceRecord | None:
        return self._records.get(tool)

    def has(self, tool: str) -> bool:
        return tool in self._records

    def fields_for(self, tool: str, default: ProvenanceRecord | None = None) -> dict[str, str]:
        rec = self._records.get(tool) or default
        if rec is None:
            rec = provenance_missing(tool)
        return rec.as_fields()

    def tool_summary_fields(self, tools: tuple[str, ...] = TOOL_PROVENANCE_TOOLS) -> dict[str, str]:
        out: dict[str, str] = {}
        for tool in tools:
            rec = self._records.get(tool) or provenance_not_used(tool)
            out.update(rec.tool_summary())
        return out

    def to_json(self) -> dict[str, Any]:
        return {
            tool: {
                "source": rec.source,
                "version": rec.version or detect_tool_version(tool),
                "status": rec.status,
                "mode": rec.mode,
                "file": rec.file,
            }
            for tool, rec in sorted(self._records.items())
        }


def without_provenance(fields: list[str]) -> list[str]:
    return [f for f in fields if f not in EVIDENCE_PROVENANCE_FIELDS]


def infer_passthrough_source(path: str | Path | None) -> str:
    from .tools.registry import ROOT

    if not path:
        return SOURCE_MISSING
    p = str(path).replace("\\", "/")
    if "/fixtures/" in p or p.startswith("data/fixtures"):
        return SOURCE_FIXTURE
    if str(ROOT) in p or "/tools/" in p:
        return SOURCE_LOCAL
    return SOURCE_EXTERNAL
