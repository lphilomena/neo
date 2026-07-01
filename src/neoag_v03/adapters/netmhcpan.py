from __future__ import annotations

from pathlib import Path

from ..utils import split_ws_or_tab, to_float, safe_id
from ..schemas import NETMHCPAN_EVIDENCE_FIELDS
from ..evidence_provenance import (
    ProvenanceRecord,
    provenance_from_file,
    without_provenance,
    write_evidence_tsv,
    MODE_PASSTHROUGH,
)

STANDARD_XLS_HEADER = [
    "Pos", "Peptide", "HLA", "Score_EL", "%Rank_EL", "Score_BA", "%Rank_BA", "BindLevel",
]


def _hla_with_star(allele: str) -> str:
    raw = str(allele or "").strip().upper()
    if not raw:
        return ""
    if not raw.startswith("HLA-"):
        raw = f"HLA-{raw}"
    body = raw[4:]
    if "*" in body:
        locus, digits = body.split("*", 1)
        compact = digits.replace(":", "")
        if len(compact) == 4 and compact.isdigit():
            return f"HLA-{locus}*{compact[:2]}:{compact[2:]}"
        return raw
    if len(body) >= 2:
        compact = body[1:].replace(":", "")
        if len(compact) == 4 and compact.isdigit():
            return f"HLA-{body[0]}*{compact[:2]}:{compact[2:]}"
        return f"HLA-{body[0]}*{body[1:]}"
    return raw


def parse_netmhcpan_local_stdout(stdout: str, source: str = "") -> list[dict[str, str]]:
    """Parse NetMHCpan 4.2 PEPTIDEMHC stdout into standard binding rows."""
    rows: list[dict[str, str]] = []
    for line in stdout.splitlines():
        if "PEPLIST" not in line:
            continue
        parts = line.split()
        try:
            peplist_idx = parts.index("PEPLIST")
        except ValueError:
            continue
        if peplist_idx < 3 or len(parts) < peplist_idx + 6:
            continue
        hla = _hla_with_star(parts[1])
        peptide = parts[2]
        nums = parts[peplist_idx + 1 :]
        if len(nums) < 5:
            continue
        el_score, el_rank, _mid, ba_rank, ic50 = nums[0], nums[1], nums[2], nums[3], nums[4]
        rows.append({
            "Peptide": peptide,
            "HLA": hla,
            "Score_EL": el_score,
            "%Rank_EL": el_rank,
            "Score_BA": ic50,
            "%Rank_BA": ba_rank,
            "BindLevel": "",
        })
    return rows


def write_netmhcpan_standard_xls(path: str | Path, rows: list[dict[str, str]]) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(STANDARD_XLS_HEADER)]
    for i, row in enumerate(rows):
        lines.append("\t".join([
            str(i),
            row.get("Peptide", ""),
            row.get("HLA", ""),
            row.get("Score_EL", ""),
            row.get("%Rank_EL", ""),
            row.get("Score_BA", ""),
            row.get("%Rank_BA", ""),
            row.get("BindLevel", ""),
        ]))
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_netmhcpan(path: str | Path, sample_id: str = "") -> list[dict[str, str]]:
    p = Path(path)
    rows = []
    header = None
    for line in p.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = split_ws_or_tab(line)
        low = [x.lower() for x in parts]
        if "peptide" in low and any(x in low for x in ["hla", "allele"]):
            header = parts
            continue
        if header:
            rec = {header[i]: parts[i] for i in range(min(len(header), len(parts)))}
            peptide = rec.get("Peptide") or rec.get("peptide") or ""
            allele = rec.get("HLA") or rec.get("Allele") or rec.get("MHC") or rec.get("allele") or ""
            allele = _hla_with_star(allele)
            if not peptide or not allele:
                continue
            ba_rank = rec.get("%Rank_BA") or rec.get("Rank_BA") or rec.get("BA_Rank") or rec.get("BA_rank") or rec.get("Rank") or "99"
            el_rank = rec.get("%Rank_EL") or rec.get("Rank_EL") or rec.get("EL_Rank") or rec.get("EL_rank") or rec.get("Rank") or "99"
            ba_score = rec.get("Score_BA") or rec.get("BA_score") or rec.get("Score") or "0"
            rows.append({
                "sample_id": sample_id,
                "peptide": peptide,
                "hla_allele": allele,
                "peptide_hla_key": safe_id(f"{peptide}_{allele}"),
                "netmhcpan_ba_score": str(to_float(ba_score, 0.0)),
                "netmhcpan_ba_rank": str(to_float(ba_rank, 99.0)),
                "netmhcpan_el_score": str(to_float(rec.get("Score_EL") or rec.get("Score"), 0.0)),
                "netmhcpan_el_rank": str(to_float(el_rank, 99.0)),
                "source_file": str(p),
            })
    return rows

def write_netmhcpan_evidence(path, rows, provenance: ProvenanceRecord | None = None):
    src = rows[0].get("source_file") if rows else path
    prov = provenance or provenance_from_file("netmhcpan", src, mode=MODE_PASSTHROUGH)
    write_evidence_tsv(path, rows, without_provenance(NETMHCPAN_EVIDENCE_FIELDS), prov)
