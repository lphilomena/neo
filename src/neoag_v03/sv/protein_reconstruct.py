from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .reference import FastaReference, GtfAnnotation, TranscriptModel, revcomp
from .sv_merge import SVCluster

CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}


def translate_dna(seq: str, stop_at_stop: bool = True) -> str:
    seq = "".join(ch for ch in seq.upper() if ch in "ACGTN")
    aa: list[str] = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        x = CODON_TABLE.get(codon, "X")
        if x == "*" and stop_at_stop:
            break
        aa.append(x)
    return "".join(aa)


def _closest_index(pos_map: list[int], genomic_pos: int) -> int:
    if not pos_map:
        return -1
    try:
        return pos_map.index(genomic_pos)
    except ValueError:
        return min(range(len(pos_map)), key=lambda i: abs(pos_map[i] - genomic_pos))


def _prefix_suffix(tx: TranscriptModel, ref: FastaReference, genomic_pos: int) -> tuple[str, str, int, str, list[int]]:
    seq, pos_map = tx.cds_sequence_and_map(ref)
    if not seq:
        return "", "", -1, "no_cds", []
    idx = _closest_index(pos_map, genomic_pos)
    if idx < 0:
        return "", "", -1, "no_position_map", pos_map
    return seq[: idx + 1], seq[idx + 1 :], idx, "ok", pos_map


def first_difference(wt: str, mut: str) -> int:
    for i, (a, b) in enumerate(zip(wt, mut)):
        if a != b:
            return i
    if len(mut) > len(wt):
        return len(wt)
    return max(0, len(mut) - 1)


@dataclass
class ProteinReconstruction:
    protein_sequence_id: str
    event_id: str
    sample_id: str
    gene: str
    transcript_id: str
    protein_type: str
    protein_sequence: str
    wt_protein_sequence: str = ""
    wt_prefix_aa: str = ""
    novel_aa: str = ""
    junction_aa_position: str = ""
    novel_start_aa: str = ""
    frameshift_start_aa: str = ""
    stop_gain_position: str = ""
    in_frame: str = "unknown"
    reconstruction_method: str = ""
    reconstruction_confidence: str = "low"
    reconstruction_reason: str = ""

    def to_row(self) -> dict[str, str]:
        return {k: str(getattr(self, k)) for k in self.__dataclass_fields__.keys()}


def reconstruct_cluster_protein(
    cluster: SVCluster,
    event_id: str,
    sample_id: str,
    annotation: GtfAnnotation,
    ref: FastaReference,
) -> tuple[ProteinReconstruction | None, dict[str, Any]]:
    rep = cluster.representative
    tx1, tx1_distance, tx1_basis = annotation.best_coding_transcript_near(rep.chrom1, rep.pos1)
    tx2, tx2_distance, tx2_basis = annotation.best_coding_transcript_near(rep.chrom2, rep.pos2)
    gene1 = tx1.gene_name if tx1 else annotation.gene_at(rep.chrom1, rep.pos1)
    gene2 = tx2.gene_name if tx2 else annotation.gene_at(rep.chrom2, rep.pos2)
    meta = {
        "gene1": gene1,
        "gene2": gene2,
        "transcript1": tx1.transcript_id if tx1 else "",
        "transcript2": tx2.transcript_id if tx2 else "",
        "effect_class": "SV_Noncoding",
        "fusion_in_frame": "unknown",
        "frameshift": "unknown",
        "junction_aa_position": "",
        "reconstruction_status": "failed",
        "reconstruction_reason": "",
    }
    if not tx1:
        meta["reconstruction_reason"] = f"breakpoint1_no_coding_transcript:{tx1_basis}"
        return None, meta

    same_gene = bool(gene1 and gene1 == gene2 and rep.chrom1 == rep.chrom2)
    svtype = rep.svtype.upper()

    try:
        if svtype in {"BND", "TRA", "CTX"} or (gene2 and gene1 and gene1 != gene2):
            if not tx2:
                meta["reconstruction_reason"] = f"fusion_partner_no_coding_transcript:{tx2_basis}"
                return None, meta
            left_prefix, _left_suffix, idx1, status1, _map1 = _prefix_suffix(tx1, ref, rep.pos1)
            _right_prefix, right_suffix, idx2, status2, _map2 = _prefix_suffix(tx2, ref, rep.pos2)
            if status1 != "ok" or status2 != "ok" or not left_prefix or not right_suffix:
                meta["reconstruction_reason"] = f"fusion_sequence_failed:{status1},{status2}"
                return None, meta
            mut_cds = left_prefix + right_suffix
            mut_prot = translate_dna(mut_cds, stop_at_stop=True)
            wt1 = translate_dna(tx1.cds_sequence_and_map(ref)[0], stop_at_stop=True)
            junction_aa = max(0, len(left_prefix) // 3)
            novel_start = min(junction_aa, max(0, len(mut_prot) - 1))
            novel = mut_prot[novel_start: novel_start + 80]
            in_frame = "yes" if len(left_prefix) % 3 == 0 else "no"
            meta.update({
                "effect_class": "SV_Fusion",
                "fusion_in_frame": in_frame,
                "frameshift": "no" if in_frame == "yes" else "yes",
                "junction_aa_position": str(junction_aa),
                "reconstruction_status": "ok",
                "reconstruction_reason": f"heuristic_fusion_prefix_suffix:{tx1_basis}:{tx1_distance}bp,{tx2_basis}:{tx2_distance}bp",
            })
            protein = ProteinReconstruction(
                protein_sequence_id=f"PROT_{event_id}",
                event_id=event_id,
                sample_id=sample_id,
                gene=f"{gene1}::{gene2}",
                transcript_id=f"{tx1.transcript_id}::{tx2.transcript_id}",
                protein_type="SV_Fusion",
                protein_sequence=mut_prot,
                wt_protein_sequence=wt1,
                wt_prefix_aa=wt1[:junction_aa],
                novel_aa=novel,
                junction_aa_position=str(junction_aa),
                novel_start_aa=str(novel_start),
                frameshift_start_aa="" if in_frame == "yes" else str(junction_aa),
                in_frame=in_frame,
                reconstruction_method="heuristic_cds_prefix_suffix",
                reconstruction_confidence="medium" if mut_prot else "low",
                reconstruction_reason=(
                    "DNA breakpoint to CDS prefix/suffix; "
                    f"lookup={tx1_basis}:{tx1_distance}bp,{tx2_basis}:{tx2_distance}bp; "
                    "validate with RNA junction when available"
                ),
            )
            return protein, meta

        # Same-gene coding SVs: first-pass deletion/insertion/duplication model.
        cds, pos_map = tx1.cds_sequence_and_map(ref)
        if not cds:
            meta["reconstruction_reason"] = "no_cds_sequence"
            return None, meta
        wt_prot = translate_dna(cds, stop_at_stop=True)
        idx1 = _closest_index(pos_map, rep.pos1)
        idx2 = _closest_index(pos_map, rep.pos2)
        if idx1 < 0:
            meta["reconstruction_reason"] = "breakpoint_not_mapped_to_cds"
            return None, meta
        i, j = sorted([idx1, idx2 if idx2 >= 0 else idx1])
        if svtype == "DEL":
            mut_cds = cds[:i] + cds[j+1:]
            deleted_len = max(0, j - i + 1)
            frameshift = deleted_len % 3 != 0
            effect = "SV_Frameshift" if frameshift else "SV_Junction"
        elif svtype == "DUP":
            duplicated = cds[i:j+1]
            mut_cds = cds[:j+1] + duplicated + cds[j+1:]
            frameshift = len(duplicated) % 3 != 0
            effect = "SV_Frameshift" if frameshift else "SV_Junction"
        elif svtype in {"INS", "BND"} and rep.inserted_sequence:
            ins = rep.inserted_sequence.upper()
            mut_cds = cds[:i+1] + ins + cds[i+1:]
            frameshift = len(ins) % 3 != 0
            effect = "SV_Insertion" if not frameshift else "SV_Frameshift"
        else:
            # Inversion/other local rearrangements are represented as a junction by reversing the affected CDS piece.
            segment = cds[i:j+1]
            mut_cds = cds[:i] + revcomp(segment) + cds[j+1:]
            frameshift = False
            effect = "SV_Junction"
        mut_prot = translate_dna(mut_cds, stop_at_stop=True)
        novel_start = first_difference(wt_prot, mut_prot)
        novel = mut_prot[novel_start: novel_start + 80]
        meta.update({
            "gene1": gene1,
            "gene2": gene2 or gene1,
            "effect_class": effect,
            "fusion_in_frame": "no",
            "frameshift": "yes" if frameshift else "no",
            "junction_aa_position": str(max(0, i // 3)),
            "reconstruction_status": "ok",
            "reconstruction_reason": f"heuristic_same_transcript_cds_edit:{tx1_basis}:{tx1_distance}bp",
        })
        protein = ProteinReconstruction(
            protein_sequence_id=f"PROT_{event_id}",
            event_id=event_id,
            sample_id=sample_id,
            gene=gene1,
            transcript_id=tx1.transcript_id,
            protein_type=effect,
            protein_sequence=mut_prot,
            wt_protein_sequence=wt_prot,
            wt_prefix_aa=wt_prot[:novel_start],
            novel_aa=novel,
            junction_aa_position=str(max(0, i // 3)),
            novel_start_aa=str(novel_start),
            frameshift_start_aa=str(novel_start) if frameshift else "",
            in_frame="no" if frameshift else "yes",
            reconstruction_method="heuristic_same_transcript_cds_edit",
            reconstruction_confidence="medium" if mut_prot else "low",
            reconstruction_reason=(
                "DNA breakpoint to transcript CDS edit; "
                f"lookup={tx1_basis}:{tx1_distance}bp; validate with RNA when available"
            ),
        )
        return protein, meta
    except Exception as exc:
        meta["reconstruction_status"] = "failed"
        meta["reconstruction_reason"] = f"exception:{type(exc).__name__}:{exc}"
        return None, meta
