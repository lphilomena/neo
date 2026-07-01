from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .protein_reconstruct import ProteinReconstruction


def read_hla_alleles(path_or_values: str | Path | list[str] | None) -> list[str]:
    if path_or_values is None:
        return []
    if isinstance(path_or_values, list):
        vals = path_or_values
    else:
        s = str(path_or_values)
        p = Path(s)
        if p.exists():
            vals = []
            for line in p.read_text(encoding="utf-8").splitlines():
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                vals.extend([x.strip() for x in line.replace(",", " ").split() if x.strip()])
        else:
            vals = [x.strip() for x in s.replace(",", " ").split() if x.strip()]
    out: list[str] = []
    for v in vals:
        if not v:
            continue
        if v not in out:
            out.append(v)
    return out


def _int_or_none(s: str) -> int | None:
    try:
        if s in {"", ".", "NA", "None"}:
            return None
        return int(float(s))
    except Exception:
        return None


@dataclass
class PeptideCandidate:
    peptide_id: str
    event_id: str
    sample_id: str
    protein_sequence_id: str
    peptide: str
    wildtype_peptide: str
    hla_allele: str
    mhc_class: str
    peptide_length: int
    peptide_start_aa: int
    peptide_end_aa: int
    crosses_junction: str
    contains_novel_aa: str
    novel_aa_count: int
    wildtype_match_status: str
    reference_proteome_match: str = "not_assessed"
    normal_hla_ligand_overlap: str = "no"

    def to_sidecar_row(self) -> dict[str, str]:
        return {
            "peptide_id": self.peptide_id,
            "event_id": self.event_id,
            "sample_id": self.sample_id,
            "protein_sequence_id": self.protein_sequence_id,
            "peptide": self.peptide,
            "wildtype_peptide": self.wildtype_peptide,
            "hla_allele": self.hla_allele,
            "mhc_class": self.mhc_class,
            "peptide_length": str(self.peptide_length),
            "peptide_start_aa": str(self.peptide_start_aa),
            "peptide_end_aa": str(self.peptide_end_aa),
            "crosses_junction": self.crosses_junction,
            "contains_novel_aa": self.contains_novel_aa,
            "novel_aa_count": str(self.novel_aa_count),
            "wildtype_match_status": self.wildtype_match_status,
            "reference_proteome_match": self.reference_proteome_match,
            "normal_hla_ligand_overlap": self.normal_hla_ligand_overlap,
        }


def build_mhc1_peptides(
    proteins: list[ProteinReconstruction],
    hla_alleles: list[str],
    *,
    lengths: tuple[int, ...] = (8, 9, 10, 11),
    normal_ligands: set[str] | None = None,
    reference_peptides: set[str] | None = None,
) -> list[PeptideCandidate]:
    normal_ligands = normal_ligands or set()
    reference_peptides = reference_peptides or set()
    out: list[PeptideCandidate] = []
    seen: set[tuple[str, str, str, int]] = set()
    for prot in proteins:
        seq = (prot.protein_sequence or "").strip().upper().replace("*", "")
        wt = (prot.wt_protein_sequence or "").strip().upper().replace("*", "")
        if not seq:
            continue
        junction = _int_or_none(prot.junction_aa_position)
        novel_start = _int_or_none(prot.novel_start_aa)
        if novel_start is None:
            novel_start = junction if junction is not None else 0
        for L in lengths:
            if len(seq) < L:
                continue
            for start0 in range(0, len(seq) - L + 1):
                end0 = start0 + L  # half-open
                peptide = seq[start0:end0]
                if "X" in peptide or "*" in peptide:
                    continue
                crosses = junction is not None and start0 < junction < end0
                contains_novel = start0 <= novel_start < end0 or start0 >= novel_start
                # For fusion, require crossing the join. For frameshift, require novel sequence.
                if prot.protein_type == "SV_Fusion" and not crosses:
                    continue
                if prot.protein_type != "SV_Fusion" and not contains_novel:
                    continue
                wt_pep = wt[start0:end0] if wt and end0 <= len(wt) else ""
                wt_status = "same_as_wt" if wt_pep == peptide and wt_pep else ("wt_window_available" if wt_pep else "no_wt_equivalent")
                if wt_status == "same_as_wt":
                    continue
                ref_match = "yes" if peptide in reference_peptides else "no" if reference_peptides else "not_assessed"
                if ref_match == "yes":
                    continue
                normal_overlap = "yes" if peptide in normal_ligands else "no"
                for hla in hla_alleles:
                    key = (prot.event_id, peptide, hla, start0)
                    if key in seen:
                        continue
                    seen.add(key)
                    pid = f"SVPEP_{prot.event_id}_{hla.replace('*','').replace(':','')}_{start0+1}_{L}_{abs(hash(peptide)) % 1000000}"
                    out.append(PeptideCandidate(
                        peptide_id=pid,
                        event_id=prot.event_id,
                        sample_id=prot.sample_id,
                        protein_sequence_id=prot.protein_sequence_id,
                        peptide=peptide,
                        wildtype_peptide=wt_pep,
                        hla_allele=hla,
                        mhc_class="I",
                        peptide_length=L,
                        peptide_start_aa=start0 + 1,
                        peptide_end_aa=end0,
                        crosses_junction="yes" if crosses else "no",
                        contains_novel_aa="yes" if contains_novel else "no",
                        novel_aa_count=max(0, end0 - max(start0, novel_start)),
                        wildtype_match_status=wt_status,
                        reference_proteome_match=ref_match,
                        normal_hla_ligand_overlap=normal_overlap,
                    ))
    return out
