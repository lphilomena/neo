from __future__ import annotations
from pathlib import Path
from ..utils import first, read_csv, read_tsv, safe_id, to_float
from ..evidence_provenance import ProvenanceRecord, provenance_from_file, without_provenance, write_evidence_tsv
from ..schemas import MHCFLURRY_EVIDENCE_FIELDS

def parse_mhcflurry(path: str | Path, sample_id: str = "") -> list[dict[str, str]]:
    p = Path(path)
    rows_in = read_csv(p) if p.suffix.lower() == ".csv" else read_tsv(p)
    out = []
    for r in rows_in:
        pep = first(r, ["peptide", "Peptide", "sequence"], "")
        allele = first(r, ["allele", "Allele", "hla", "HLA"], "")
        if pep and allele:
            out.append({
                "sample_id": sample_id,
                "peptide": pep,
                "hla_allele": allele,
                "peptide_hla_key": safe_id(f"{pep}_{allele}"),
                "mhcflurry_affinity": str(to_float(first(r, ["mhcflurry_affinity", "affinity", "prediction", "ic50"], "0"), 0.0)),
                "mhcflurry_affinity_percentile": str(to_float(first(r, ["mhcflurry_affinity_percentile", "affinity_percentile", "percentile_rank", "percentile"], "99"), 99.0)),
                "mhcflurry_processing_score": str(to_float(first(r, ["mhcflurry_processing_score", "processing_score"], "0"), 0.0)),
                "mhcflurry_presentation_score": str(to_float(first(r, ["mhcflurry_presentation_score", "presentation_score"], "0"), 0.0)),
                "source_file": str(p),
            })
    return out

def write_mhcflurry_evidence(path, rows, provenance: ProvenanceRecord | None = None):
    src = rows[0].get("source_file") if rows else path
    prov = provenance or provenance_from_file("mhcflurry", src, mode="passthrough")
    write_evidence_tsv(path, rows, without_provenance(MHCFLURRY_EVIDENCE_FIELDS), prov)
