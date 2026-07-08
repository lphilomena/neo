"""WES SV Phase 1.5 adapter — exome-capture-limited SV neoantigen triage."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

from .phase1 import build_sv_phase1_raw


def build_sv_wes_phase1_5_raw(
    *,
    sample_id: str,
    sv_vcfs: Iterable[str | Path],
    reference_fasta: str | Path,
    gencode_gtf: str | Path,
    hla: str | Path | list[str],
    outdir: str | Path,
    profile_name: str = "sv_wes_phase1_5",
    callers: Iterable[str] | None = None,
    tumor_sample_name: str | None = None,
    normal_sample_name: str | None = None,
    expression_tsv: str | Path | None = None,
    rna_junction_tsv: str | Path | None = None,
    normal_expression_tsv: str | Path | None = None,
    normal_hla_ligands_tsv: str | Path | None = None,
    merge_distance_bp: int = 200,
    allow_tier2: bool = True,
    capture_bed: str | Path | None = None,
    capture_near_bp: int = 250,
    capture_slop_bp: int = 1000,
) -> dict[str, str]:
    """Build raw_events/raw_peptides using WES Phase 1.5 filters and profile."""
    return build_sv_phase1_raw(
        sample_id=sample_id,
        sv_vcfs=sv_vcfs,
        reference_fasta=reference_fasta,
        gencode_gtf=gencode_gtf,
        hla=hla,
        outdir=outdir,
        profile_name=profile_name,
        callers=callers,
        tumor_sample_name=tumor_sample_name,
        normal_sample_name=normal_sample_name,
        expression_tsv=expression_tsv,
        rna_junction_tsv=rna_junction_tsv,
        normal_expression_tsv=normal_expression_tsv,
        normal_hla_ligands_tsv=normal_hla_ligands_tsv,
        merge_distance_bp=merge_distance_bp,
        allow_tier2=allow_tier2,
        wes_mode=True,
        capture_bed=capture_bed,
        capture_near_bp=capture_near_bp,
        capture_slop_bp=capture_slop_bp,
    )


class WESAdapter:
    """WES SV Phase 1.5 adapter for Project B / neoag-v03."""

    def build_raw(self, **kwargs) -> dict[str, str]:
        return build_sv_wes_phase1_5_raw(**kwargs)
