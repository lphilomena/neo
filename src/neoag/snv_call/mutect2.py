"""GATK Mutect2 paired calling for WES Phase 1."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path


def _resolve_gatk(explicit: str | None = None) -> str:
    if explicit:
        return explicit
    found = shutil.which("gatk")
    if not found:
        raise FileNotFoundError("gatk not found on PATH; install GATK4 or set tools.executables.gatk")
    return found


def ensure_ref_index(reference_fasta: str | Path, *, gatk: str | None = None) -> Path:
    """Ensure FASTA .fai and .dict exist (creates them if missing)."""
    ref = Path(reference_fasta).resolve()
    if not ref.is_file():
        raise FileNotFoundError(f"Reference FASTA not found: {ref}")
    fai = ref.with_suffix(ref.suffix + ".fai") if ref.suffix else Path(str(ref) + ".fai")
    if not fai.is_file():
        samtools = shutil.which("samtools")
        if not samtools:
            raise FileNotFoundError(f"Missing {fai.name}; install samtools or provide pre-indexed reference")
        subprocess.run([samtools, "faidx", str(ref)], check=True)
    dict_path = ref.with_suffix(".dict")
    if not dict_path.is_file():
        gatk_bin = _resolve_gatk(gatk)
        subprocess.run([gatk_bin, "CreateSequenceDictionary", "-R", str(ref), "-O", str(dict_path)], check=True)
    return ref


def run_mutect2_paired(
    *,
    tumor_bam: str | Path,
    normal_bam: str | Path,
    reference_fasta: str | Path,
    intervals_bed: str | Path,
    tumor_sample_name: str,
    normal_sample_name: str,
    out_vcf: str | Path,
    workdir: str | Path | None = None,
    gatk: str | None = None,
) -> Path:
    """Run GATK Mutect2 on tumor-normal WES BAMs."""
    ref = ensure_ref_index(reference_fasta, gatk=gatk)
    tumor_bam = Path(tumor_bam).resolve()
    normal_bam = Path(normal_bam).resolve()
    intervals_bed = Path(intervals_bed).resolve()
    out_vcf = Path(out_vcf)
    work = Path(workdir or out_vcf.parent)
    work.mkdir(parents=True, exist_ok=True)
    out_vcf.parent.mkdir(parents=True, exist_ok=True)

    for label, bam in (("tumor", tumor_bam), ("normal", normal_bam)):
        if not bam.is_file():
            raise FileNotFoundError(f"{label} BAM not found: {bam}")
        bai = Path(f"{bam}.bai")
        if not bai.is_file():
            raise FileNotFoundError(f"Missing BAM index: {bai}")

    if not intervals_bed.is_file():
        raise FileNotFoundError(f"WES intervals BED not found: {intervals_bed}")

    gatk_bin = _resolve_gatk(gatk)
    cmd = [
        gatk_bin,
        "Mutect2",
        "-R",
        str(ref),
        "-I",
        str(tumor_bam),
        "-I",
        str(normal_bam),
        "-normal",
        normal_sample_name,
        "-L",
        str(intervals_bed),
        "-O",
        str(out_vcf),
    ]
    subprocess.run(cmd, cwd=work, check=True)
    if not out_vcf.is_file():
        raise FileNotFoundError(f"Mutect2 did not produce: {out_vcf}")
    return out_vcf


def run_filter_mutect_calls(
    *,
    reference_fasta: str | Path,
    unfiltered_vcf: str | Path,
    out_vcf: str | Path,
    workdir: str | Path | None = None,
    gatk: str | None = None,
    germline_resource: str | Path | None = None,
    panel_of_normals: str | Path | None = None,
) -> Path:
    """Run GATK FilterMutectCalls on Mutect2 output."""
    ref = ensure_ref_index(reference_fasta, gatk=gatk)
    unfiltered_vcf = Path(unfiltered_vcf).resolve()
    out_vcf = Path(out_vcf)
    work = Path(workdir or out_vcf.parent)
    work.mkdir(parents=True, exist_ok=True)
    out_vcf.parent.mkdir(parents=True, exist_ok=True)

    if not unfiltered_vcf.is_file():
        raise FileNotFoundError(f"Unfiltered VCF not found: {unfiltered_vcf}")

    gatk_bin = _resolve_gatk(gatk)
    cmd = [
        gatk_bin,
        "FilterMutectCalls",
        "-R",
        str(ref),
        "-V",
        str(unfiltered_vcf),
        "-O",
        str(out_vcf),
    ]
    if germline_resource:
        gr = Path(germline_resource)
        if gr.is_file():
            cmd.extend(["--germline-resource", str(gr)])
    if panel_of_normals:
        pon = Path(panel_of_normals)
        if pon.is_file():
            cmd.extend(["--panel-of-normals", str(pon)])

    subprocess.run(cmd, cwd=work, check=True)
    if not out_vcf.is_file():
        raise FileNotFoundError(f"FilterMutectCalls did not produce: {out_vcf}")
    return out_vcf
