from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from typing import Any

from neoag.controlled_execution.io_utils import write_json, write_tsv


def _fastq_values(value: Any) -> list[str]:
    values = value if isinstance(value, list) else ([value] if value else [])
    return [str(item) for item in values if item]


def _read_number(path: str) -> int | None:
    name = Path(path).name.lower()
    read1_tokens = ("_r1", ".r1", "-r1", "_1.f", ".1.f", "-1.f")
    read2_tokens = ("_r2", ".r2", "-r2", "_2.f", ".2.f", "-2.f")
    if any(token in name for token in read1_tokens):
        return 1
    if any(token in name for token in read2_tokens):
        return 2
    return None


def _fastq_batches(value: Any) -> list[tuple[str, str]]:
    values = _fastq_values(value)
    if len(values) < 2 or len(values) % 2:
        return []
    read1 = [item for item in values if _read_number(item) == 1]
    read2 = [item for item in values if _read_number(item) == 2]
    if len(read1) == len(read2) and len(read1) + len(read2) == len(values):
        return list(zip(read1, read2))
    return [(values[i], values[i + 1]) for i in range(0, len(values), 2)]


def _safe_id(value: Any) -> str:
    safe = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in str(value or "sample"))
    return safe.strip("._-") or "sample"


def prepare_rna_evidence(
    inputs: dict[str, Any],
    *,
    project_root: str | Path,
    outdir: str | Path,
    execute: bool,
    method: str = "auto",
    timeout: int = 7200,
) -> dict[str, Any]:
    """Plan or run RNA quantification and allele counting behind the Gateway boundary."""
    root = Path(project_root).resolve()
    od = Path(outdir).resolve()
    od.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []
    outputs: dict[str, str] = {}

    for label, key in (("gene_tpm", "expression_tsv"), ("transcript_tpm", "transcript_expression_tsv"), ("rna_alt_vaf", "rna_evidence_tsv")):
        path = str(inputs.get(key) or "")
        if path and Path(path).is_file():
            outputs[label] = path
            rows.append({"stage": label, "status": "REUSED", "command_preview": "", "output": path, "message": "declared evidence reused"})

    batches = _fastq_batches(inputs.get("tumor_rna_fastq"))
    pair = batches[0] if len(batches) == 1 else None
    if len(batches) > 1:
        merge_dir = od / "merged_fastq"
        sample = _safe_id(inputs.get("sample_id") or "sample")
        pair = (str(merge_dir / f"{sample}_R1.fq.gz"), str(merge_dir / f"{sample}_R2.fq.gz"))
        merge_command = ["bash", str(root / "scripts/run_merge_paired_fastq.sh"), "--outdir", str(merge_dir), "--sample-id", sample]
        for r1, r2 in batches:
            merge_command.extend(["--fastq1", r1, "--fastq2", r2])
        if not execute:
            rows.append({"stage": "rna_fastq_merge", "status": "PLANNED", "command_preview": " ".join(merge_command), "output": str(merge_dir), "message": "merge multi-batch paired RNA FASTQ before quantification"})
        else:
            proc = subprocess.run(merge_command, cwd=root, text=True, capture_output=True, timeout=timeout)
            (od / "rna_fastq_merge.log").write_text(proc.stdout + ("\n--- STDERR ---\n" if proc.stderr else "") + proc.stderr, encoding="utf-8")
            status = "PASS" if proc.returncode == 0 and Path(pair[0]).is_file() and Path(pair[1]).is_file() else "FAILED"
            rows.append({"stage": "rna_fastq_merge", "status": status, "command_preview": " ".join(merge_command), "output": str(merge_dir), "message": "merged R1/R2 FASTQ" if status == "PASS" else f"returncode={proc.returncode}"})
            if status != "PASS":
                pair = None
    quant_method = method
    if quant_method == "auto":
        quant_method = "salmon" if inputs.get("salmon_index") and inputs.get("tx2gene") else ("rsem" if inputs.get("rsem_reference") else "")
    if pair and not (outputs.get("gene_tpm") and outputs.get("transcript_tpm")):
        quant_dir = od / (quant_method or "quantification_unconfigured")
        if quant_method == "salmon":
            command = ["bash", str(root / "scripts/run_salmon_fastq_to_tpm.sh"), "--fastq1", pair[0], "--fastq2", pair[1], "--outdir", str(quant_dir), "--sample-id", str(inputs.get("sample_id") or "sample"), "--salmon-index", str(inputs.get("salmon_index") or ""), "--tx2gene", str(inputs.get("tx2gene") or "")]
        elif quant_method == "rsem":
            command = ["bash", str(root / "scripts/run_rsem_fastq_to_tpm.sh"), "--fastq1", pair[0], "--fastq2", pair[1], "--outdir", str(quant_dir), "--sample-id", str(inputs.get("sample_id") or "sample"), "--rsem-reference", str(inputs.get("rsem_reference") or "")]
        else:
            command = []
        if not command:
            rows.append({"stage": "rna_quantification", "status": "UNASSESSED", "command_preview": "", "output": "", "message": "RNA FASTQ supplied but Salmon/RSEM reference was not declared"})
        elif not execute:
            rows.append({"stage": "rna_quantification", "status": "PLANNED", "command_preview": " ".join(command), "output": str(quant_dir), "message": "Gateway approval required for execution"})
        else:
            proc = subprocess.run(command, cwd=root, text=True, capture_output=True, timeout=timeout)
            (od / "rna_quantification.log").write_text(proc.stdout + ("\n--- STDERR ---\n" if proc.stderr else "") + proc.stderr, encoding="utf-8")
            gene = quant_dir / "gene_tpm.tsv"
            transcript = quant_dir / "transcript_tpm.tsv"
            status = "PASS" if proc.returncode == 0 and gene.is_file() and transcript.is_file() else "FAILED"
            rows.append({"stage": "rna_quantification", "status": status, "command_preview": " ".join(command), "output": str(quant_dir), "message": "gene and transcript TPM" if status == "PASS" else f"returncode={proc.returncode}"})
            if status == "PASS":
                outputs.update({"gene_tpm": str(gene), "transcript_tpm": str(transcript)})

    rna_bam = str(inputs.get("tumor_rna_bam") or "")
    somatic_vcf = str(inputs.get("somatic_vcf") or "")
    if rna_bam and somatic_vcf and not outputs.get("rna_alt_vaf"):
        alt_out = od / "rna_alt_vaf.tsv"
        command = [sys.executable, str(root / "scripts/rna_allele_counts_pysam.py"), "--somatic-vcf", somatic_vcf, "--rna-bam", rna_bam, "--output-tsv", str(alt_out)]
        if not execute:
            rows.append({"stage": "rna_alt_vaf", "status": "PLANNED", "command_preview": " ".join(command), "output": str(alt_out), "message": "Gateway approval required for execution"})
        else:
            proc = subprocess.run(command, cwd=root, text=True, capture_output=True, timeout=timeout)
            (od / "rna_alt_vaf.log").write_text(proc.stdout + ("\n--- STDERR ---\n" if proc.stderr else "") + proc.stderr, encoding="utf-8")
            status = "PASS" if proc.returncode == 0 and alt_out.is_file() else "FAILED"
            rows.append({"stage": "rna_alt_vaf", "status": status, "command_preview": " ".join(command), "output": str(alt_out), "message": "RNA ref/alt reads, depth and VAF" if status == "PASS" else f"returncode={proc.returncode}"})
            if status == "PASS":
                outputs["rna_alt_vaf"] = str(alt_out)

    for label, key in (("fusion_junction_reads", "fusion_tsv"), ("splice_junction_reads", "splice_junction_tsv")):
        path = str(inputs.get(key) or "")
        if path and Path(path).is_file():
            outputs[label] = path
            rows.append({"stage": label, "status": "REUSED", "command_preview": "", "output": path, "message": "junction evidence retained from caller table"})

    write_tsv(od / "rna_preprocessing_status.tsv", rows)
    status = "SKIPPED" if not rows else ("FAILED" if any(row["status"] == "FAILED" for row in rows) else ("PARTIAL" if any(row["status"] == "UNASSESSED" for row in rows) else ("PASS" if execute else "PLANNED")))
    summary = {"status": status, "execute": execute, "method": quant_method or "UNCONFIGURED", "outputs": outputs, "stages": rows}
    write_json(od / "rna_preprocessing_summary.json", summary)
    return summary
