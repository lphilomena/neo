from __future__ import annotations

from pathlib import Path

from ..utils import safe_id, split_ws_or_tab, to_float
from ..evidence_provenance import ProvenanceRecord, provenance_from_file, without_provenance, write_evidence_tsv
from ..schemas import NETMHCSTABPAN_EVIDENCE_FIELDS


def parse_netmhcstabpan(path: str | Path, sample_id: str = "") -> list[dict[str, str]]:
    """Parse NetMHCstabpan output (IEDB API TSV or consolidated benchmark TSV)."""
    p = Path(path)
    rows: list[dict[str, str]] = []
    header: list[str] | None = None
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
        elif len(parts) >= 8 and parts[5].isalpha():
            rec = {
                "allele": parts[0],
                "peptide": parts[5],
                "score": parts[6],
                "percentile_rank": parts[7],
            }
        elif len(parts) >= 4:
            rec = {
                "Peptide": parts[0],
                "HLA": parts[1],
                "score": parts[2],
                "percentile_rank": parts[3],
            }
        else:
            continue
        peptide = rec.get("Peptide") or rec.get("peptide") or ""
        allele = rec.get("HLA") or rec.get("allele") or rec.get("Allele") or ""
        if not peptide or not allele:
            continue
        rank = (
            rec.get("percentile_rank")
            or rec.get("Percentile_rank")
            or rec.get("%Rank")
            or rec.get("rank")
            or "99"
        )
        score = rec.get("score") or rec.get("Score") or "0"
        rows.append({
            "sample_id": sample_id,
            "peptide": peptide,
            "hla_allele": allele,
            "peptide_hla_key": safe_id(f"{peptide}_{allele}"),
            "netmhcstabpan_score": str(to_float(score, 0.0)),
            "netmhcstabpan_rank": str(to_float(rank, 99.0)),
            "source_file": str(p),
        })
    return rows


def write_netmhcstabpan_evidence(path, rows, provenance: ProvenanceRecord | None = None):
    src = rows[0].get("source_file") if rows else path
    prov = provenance or provenance_from_file("netmhcstabpan", src, mode="passthrough")
    write_evidence_tsv(path, rows, without_provenance(NETMHCSTABPAN_EVIDENCE_FIELDS), prov)
