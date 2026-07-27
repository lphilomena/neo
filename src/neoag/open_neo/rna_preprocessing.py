from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from typing import Any

from neoag.controlled_execution.io_utils import write_json, write_tsv


def _fastq_pair(value: Any) -> tuple[str, str] | None:
    values = value if isinstance(value, list) else ([value] if value else [])
    values = [str(item) for item in values if item]
    return (values[0], values[1]) if len(values) >= 2 else None


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

    pair = _fastq_pair(inputs.get("tumor_rna_fastq"))
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
