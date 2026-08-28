from __future__ import annotations

import json
import os
import re
import shlex
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any


PROFILE_NAME = "rna_fusion_splice_v1"


@dataclass(frozen=True)
class ProfileRequirement:
    field: str
    status: str
    required: bool
    detail: str


def _fastq_values(value: Any) -> list[str]:
    values = value if isinstance(value, list) else ([value] if value else [])
    return [str(item) for item in values if str(item)]


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


def _pair(value: Any) -> tuple[str, str] | None:
    batches = _fastq_batches(value)
    return batches[0] if len(batches) == 1 else None


def _safe_id(value: Any) -> str:
    safe = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in str(value or "sample"))
    return safe.strip("._-") or "sample"


def is_rna_fastq_profile_candidate(inputs: dict[str, Any]) -> bool:
    batches = _fastq_batches(inputs.get("tumor_rna_fastq"))
    dna_inputs = any(inputs.get(key) for key in (
        "tumor_dna_bam", "normal_dna_bam", "tumor_dna_fastq", "normal_dna_fastq"
    ))
    return bool(batches) and not dna_inputs and not inputs.get("production_manifest")


def profile_requirements(inputs: dict[str, Any]) -> list[ProfileRequirement]:
    batches = _fastq_batches(inputs.get("tumor_rna_fastq"))
    rows: list[ProfileRequirement] = []

    def add(field: str, value: Any, *, required: bool, detail: str) -> None:
        if isinstance(value, list):
            present = bool(value) and all(Path(str(item)).exists() for item in value)
        elif field == "hla_alleles_or_hla_file":
            hla_file = str(inputs.get("hla_file") or "")
            present = bool(inputs.get("hla_alleles")) or bool(hla_file and Path(hla_file).is_file())
        else:
            present = bool(value) and Path(str(value)).exists()
        rows.append(ProfileRequirement(field, "READY" if present else ("MISSING" if required else "UNASSESSED"), required, detail))

    fastqs = _fastq_values(inputs.get("tumor_rna_fastq"))
    fastqs_present = bool(batches) and all(Path(item).exists() for item in fastqs)
    rows.append(ProfileRequirement(
        "tumor_rna_fastq_pair",
        "READY" if fastqs_present else "MISSING",
        True,
        "one or more paired-end tumor RNA FASTQ batches, supplied as R1/R2 pairs",
    ))
    add("hla_alleles_or_hla_file", None, required=True, detail="required for peptide-HLA prediction")
    add("reference_fasta", inputs.get("reference_fasta"), required=True, detail="GRCh38 FASTA used by STAR and downstream reconstruction")
    add("gencode_gtf", inputs.get("gencode_gtf"), required=True, detail="GENCODE GTF matching the FASTA and indexes")
    add("star_index", inputs.get("star_index"), required=True, detail="STAR genome index matching FASTA/GTF")
    add("easyfuse_ref", inputs.get("easyfuse_ref"), required=True, detail="EasyFuse reference bundle")
    salmon_ready = bool(
        inputs.get("salmon_index") and Path(str(inputs["salmon_index"])).is_dir()
        and inputs.get("tx2gene") and Path(str(inputs["tx2gene"])).is_file()
    )
    rsem_prefix = str(inputs.get("rsem_reference") or "")
    rsem_ready = bool(
        rsem_prefix and Path(rsem_prefix).parent.is_dir()
        and list(Path(rsem_prefix).parent.glob(Path(rsem_prefix).name + ".*"))
    )
    rows.append(ProfileRequirement(
        "expression_reference", "READY" if salmon_ready or rsem_ready else "MISSING", True,
        "Salmon index plus tx2gene, or a matching RSEM reference prefix",
    ))
    add("salmon_index", inputs.get("salmon_index"), required=False, detail="Salmon transcriptome index; paired with tx2gene")
    add("tx2gene", inputs.get("tx2gene"), required=False, detail="transcript-to-gene mapping; paired with Salmon index")
    rows.append(ProfileRequirement("rsem_reference", "READY" if rsem_ready else "UNASSESSED", False, "RSEM reference prefix matching FASTA/GTF"))
    add("normal_junctions", inputs.get("normal_junctions"), required=False, detail="normal-junction background; missing keeps splice safety partial")
    add("normal_readthrough", inputs.get("normal_readthrough"), required=False, detail="normal/read-through fusion background")
    add("normal_expression", inputs.get("normal_expression"), required=False, detail="normal tissue/HSPC expression background")
    add("normal_hla_ligands", inputs.get("normal_hla_ligands"), required=False, detail="normal HLA ligandome background")
    add("reference_proteome", inputs.get("reference_proteome"), required=False, detail="reference proteome exact-match safety check")
    add("snaf_db", inputs.get("snaf_db"), required=True, detail="official SNAF reference root containing Alt91_db and controls")
    add("snaf_workflow", inputs.get("snaf_workflow"), required=False, detail="optional site override; built-in workflow is used when snaf_db is configured")
    custom_splicemutr = bool(inputs.get("splicemutr_workflow"))
    add("splicemutr_workflow", inputs.get("splicemutr_workflow"), required=True, detail="NeoAg production SpliceMutr workflow")
    add("splicemutr_config", inputs.get("splicemutr_config"), required=not custom_splicemutr,
        detail="required cohort/normal-control configuration for the bundled workflow")
    return rows


def _q(value: Any) -> str:
    return shlex.quote(str(value))


def _toml(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, list):
        return "[" + ", ".join(json.dumps(str(item)) for item in value) + "]"
    return json.dumps(str(value))


def _stage(lines: list[str], name: str, *, required: bool, command: str,
           outputs: dict[str, str], depends_on: list[str] | None = None,
           source: str = "", data_row_outputs: list[str] | None = None) -> None:
    thread_match = re.search(r"--(?:threads|cores)(?:=|\s+)(\d+)", command)
    ceiling = int(thread_match.group(1)) if thread_match else 16
    key = name.lower()
    if "easyfuse" in key:
        cpus, memory_gb = 8, 64.0
    elif any(token in key for token in ("rna_alignment", "star_fusion", "fusioncatcher", "arriba")):
        cpus, memory_gb = 8, 48.0
    elif "rsem" in key:
        cpus, memory_gb = 4, 32.0
    elif any(token in key for token in ("rna_expression", "salmon")):
        cpus, memory_gb = 4, 16.0
    elif any(token in key for token in ("snaf", "splicemutr")):
        cpus, memory_gb = 6, 48.0
    elif "fastq_qc" in key:
        cpus, memory_gb = 4, 8.0
    elif "fastq_merge" in key:
        cpus, memory_gb = 2, 8.0
    else:
        cpus, memory_gb = 2, 8.0
    cpus = max(1, min(cpus, ceiling))
    command = re.sub(r"(--threads(?:=|\s+))\d+", lambda match: match.group(1) + str(cpus), command)
    command = re.sub(r"(--cores(?:=|\s+))\d+", lambda match: match.group(1) + str(cpus), command)
    lines += ["", f"[stages.{name}]", f"required = {_toml(required)}"]
    if source:
        lines.append(f"source = {_toml(source)}")
    if depends_on:
        lines.append(f"depends_on = {_toml(depends_on)}")
    if data_row_outputs:
        lines.append(f"data_row_outputs = {_toml(data_row_outputs)}")
    lines.append(f"cpus = {cpus}")
    lines.append(f"memory_gb = {memory_gb}")
    lines.append(f"command = {_toml(command)}")
    lines.append(f"[stages.{name}.outputs]")
    lines.extend(f"{key} = {_toml(value)}" for key, value in outputs.items())


def generate_rna_fusion_splice_manifest(
    inputs: dict[str, Any],
    path: str | Path,
    *,
    project_root: str | Path,
    outdir: str | Path,
) -> dict[str, Any]:
    batches = _fastq_batches(inputs.get("tumor_rna_fastq"))
    if not batches:
        raise ValueError("rna_fusion_splice profile requires paired tumor RNA FASTQ files supplied as R1/R2 pairs")
    requirements = profile_requirements(inputs)
    missing_required = [row.field for row in requirements if row.required and row.status != "READY"]
    sample_id = str(inputs.get("sample_id") or "SAMPLE001")
    safe_sample_id = _safe_id(sample_id)
    merge_fastq = len(batches) > 1
    pair = (
        (f"{{outdir}}/rna/merged_fastq/{safe_sample_id}_R1.fq.gz", f"{{outdir}}/rna/merged_fastq/{safe_sample_id}_R2.fq.gz")
        if merge_fastq else batches[0]
    )
    hla_file = str(inputs.get("hla_file") or "")
    hla_alleles = [str(value) for value in inputs.get("hla_alleles") or []]
    threads = int(inputs.get("rna_threads") or 16)
    root = Path(project_root).resolve()
    run_out = Path(outdir).resolve()
    manifest_path = Path(path)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    if not hla_file and hla_alleles:
        generated_hla = manifest_path.with_name("rna_fusion_splice.hla.txt")
        generated_hla.write_text("\n".join(hla_alleles) + "\n", encoding="utf-8")
        hla_file = str(generated_hla.resolve())

    def script(name: str) -> str:
        return _q(root / "scripts" / name)

    lines = [
        "# Auto-generated by open-neo-run. Review before execute.",
        "# SNAF and SpliceMutr are default splice stages; missing assets block execution.",
        "[run]",
        f"sample_id = {_toml(sample_id)}",
        f"profile = {_toml(PROFILE_NAME)}",
        f"outdir = {_toml(str(run_out))}",
        f"hla_file = {_toml(hla_file)}",
        f"hla_alleles = {_toml(hla_alleles)}",
        "tools_stub = false",
        "immunogenicity_stub = false",
        f"expected_peptide_sources = {_toml(['EasyFuse', 'SNAF', 'SpliceMutr'])}",
        'presentation_predictors = ["netmhcpan", "mhcflurry", "netmhcstabpan", "netchop"]',
        'required_presentation_predictors = ["netmhcpan", "mhcflurry", "netmhcstabpan", "netchop"]',
        'reports = "patient,technical"',
        f"total_cpus = {int(inputs.get('total_cpus') or threads)}",
        f"total_memory_gb = {float(inputs.get('total_memory_gb') or 0)}",
        f"max_parallel_stages = {int(inputs.get('max_parallel_stages') or 3)}",
        "default_stage_cpus = 1",
        "default_stage_memory_gb = 4.0",
    ]

    if merge_fastq:
        merge_args = " ".join(
            f"--fastq1 {_q(r1)} --fastq2 {_q(r2)}" for r1, r2 in batches
        )
        merge_command = (
            f"bash {script('run_merge_paired_fastq.sh')} --sample-id {_q(safe_sample_id)} "
            f"--outdir {{outdir}}/rna/merged_fastq {merge_args}"
        )
        _stage(lines, "fastq_merge", required=True, command=merge_command,
               outputs={
                   "fastq1": f"{{outdir}}/rna/merged_fastq/{safe_sample_id}_R1.fq.gz",
                   "fastq2": f"{{outdir}}/rna/merged_fastq/{safe_sample_id}_R2.fq.gz",
                   "summary": "{outdir}/rna/merged_fastq/merge_fastq.summary.tsv",
               })

    qc_command = (
        f"bash {script('run_rna_fastq_qc.sh')} --fastq1 {_q(pair[0])} --fastq2 {_q(pair[1])} "
        f"--sample-id {_q(sample_id)} --threads {threads} --outdir {{outdir}}/rna/qc"
    )
    _stage(lines, "fastq_qc", required=True, command=qc_command,
           outputs={"qc_complete": "{outdir}/rna/qc/qc.complete.json"},
           depends_on=["fastq_merge"] if merge_fastq else None)

    star_command = ""
    if inputs.get("star_index") and inputs.get("gencode_gtf"):
        star_command = (
            f"bash {script('run_star_rna_fastq.sh')} --fastq1 {_q(pair[0])} --fastq2 {_q(pair[1])} "
            f"--star-index {_q(inputs['star_index'])} --gtf {_q(inputs['gencode_gtf'])} "
            f"--sample-id {_q(sample_id)} --threads {threads} --outdir {{outdir}}/rna/star"
        )
    _stage(lines, "rna_alignment", required=True, command=star_command,
           outputs={
               "rna_bam": "{outdir}/rna/star/Aligned.sortedByCoord.out.bam",
               "chimeric_junction": "{outdir}/rna/star/Chimeric.out.junction",
               "star_junctions": "{outdir}/rna/star/SJ.out.tab",
           }, depends_on=["fastq_qc"])

    expression_command = ""
    quant_method = str(inputs.get("rna_quant_method") or "auto")
    salmon_primary = quant_method in {"auto", "salmon"} and inputs.get("salmon_index") and inputs.get("tx2gene")
    rsem_available = bool(inputs.get("rsem_reference"))
    if salmon_primary:
        expression_command = (
            f"bash {script('run_salmon_fastq_to_tpm.sh')} --fastq1 {_q(pair[0])} --fastq2 {_q(pair[1])} "
            f"--sample-id {_q(sample_id)} --threads {threads} --salmon-index {_q(inputs['salmon_index'])} "
            f"--tx2gene {_q(inputs['tx2gene'])} --outdir {{outdir}}/rna/expression"
        )
    elif quant_method in {"auto", "rsem"} and rsem_available:
        expression_command = (
            f"bash {script('run_rsem_fastq_to_tpm.sh')} --fastq1 {_q(pair[0])} --fastq2 {_q(pair[1])} "
            f"--sample-id {_q(sample_id)} --threads {threads} --rsem-reference {_q(inputs['rsem_reference'])} "
            "--outdir {outdir}/rna/expression"
        )
    _stage(lines, "rna_expression", required=True, command=expression_command,
           outputs={
               "expression": "{outdir}/rna/expression/gene_tpm.tsv",
               "transcript_expression": "{outdir}/rna/expression/transcript_tpm.tsv",
           }, depends_on=["fastq_qc"])

    if salmon_primary and rsem_available and quant_method == "auto":
        rsem_crosscheck_command = (
            f"bash {script('run_rsem_fastq_to_tpm.sh')} --fastq1 {_q(pair[0])} --fastq2 {_q(pair[1])} "
            f"--sample-id {_q(sample_id)} --threads {threads} --rsem-reference {_q(inputs['rsem_reference'])} "
            "--outdir {outdir}/rna/rsem_expression"
        )
        _stage(lines, "rsem_expression_crosscheck", required=False, command=rsem_crosscheck_command,
               outputs={
                   "rsem_expression": "{outdir}/rna/rsem_expression/gene_tpm.tsv",
                   "rsem_transcript_expression": "{outdir}/rna/rsem_expression/transcript_tpm.tsv",
               }, depends_on=["fastq_qc"])

    easyfuse_ref = str(inputs.get("easyfuse_ref") or "")
    easyfuse_command = (
        f"EASYFUSE_SAMPLE_ID={_q(sample_id)} EASYFUSE_FQ1={_q(pair[0])} EASYFUSE_FQ2={_q(pair[1])} "
        f"NEOAG_EASYFUSE_REF={_q(easyfuse_ref)} OUTDIR={{outdir}}/branches/fusion/easyfuse "
        f"bash {script('run_easyfuse_sample.sh')}"
    ) if easyfuse_ref else ""
    _stage(lines, "easyfuse_discovery", required=True, command=easyfuse_command,
           outputs={"fusion_tsv": f"{{outdir}}/branches/fusion/easyfuse/{sample_id}/fusions.pass.csv"},
           depends_on=["fastq_qc"])

    fusion_union_command = (
        f"PYTHONPATH={_q(root / 'src')} {_q(Path(sys.executable))} "
        f"{script('build_easyfuse_candidate_union.py')} "
        f"--pass-csv {{outdir}}/branches/fusion/easyfuse/{_q(sample_id)}/fusions.pass.csv "
        f"--all-csv {{outdir}}/branches/fusion/easyfuse/{_q(sample_id)}/fusions.csv "
        f"--output {{outdir}}/branches/fusion/easyfuse/{_q(sample_id)}/fusions.openneo_union.csv "
        f"--summary {{outdir}}/branches/fusion/easyfuse/{_q(sample_id)}/fusions.openneo_union.summary.json"
    )
    _stage(lines, "fusion_candidate_union", required=True, command=fusion_union_command,
           outputs={
               "fusion_tsv": f"{{outdir}}/branches/fusion/easyfuse/{sample_id}/fusions.openneo_union.csv",
               "summary": f"{{outdir}}/branches/fusion/easyfuse/{sample_id}/fusions.openneo_union.summary.json",
           }, depends_on=["easyfuse_discovery"])

    fusion_consensus_command = (
        f"{_q(Path(sys.executable))} {script('review_rna_fusions.py')} "
        f"--easyfuse {{outdir}}/branches/fusion/easyfuse/{_q(sample_id)}/fusions.openneo_union.csv "
        f"--normal-readthrough {_q(inputs.get('normal_readthrough') or '')} "
        f"--outdir {{outdir}}/branches/fusion/consensus"
    )
    _stage(lines, "fusion_cross_validation", required=True, command=fusion_consensus_command,
           outputs={
               "fusion_consensus": "{outdir}/branches/fusion/consensus/fusion_consensus.tsv",
               "fusion_background_review": "{outdir}/branches/fusion/consensus/fusion_background_review.tsv",
           },
           depends_on=["fusion_candidate_union"])

    regtools_command = (
        f"bash {script('run_regtools_junctions.sh')} --bam {{outdir}}/rna/star/Aligned.sortedByCoord.out.bam "
        f"--sample-id {_q(sample_id)} --out {{outdir}}/branches/splice/regtools_junctions.tsv"
    )
    _stage(lines, "junction_extraction", required=True, command=regtools_command,
           outputs={"junctions": "{outdir}/branches/splice/regtools_junctions.tsv"},
           depends_on=["rna_alignment"])

    snaf_workflow = str(inputs.get("snaf_workflow") or os.environ.get("SNAF_WORKFLOW") or "")
    snaf_db = str(inputs.get("snaf_db") or os.environ.get("NEOAG_SNAF_DB") or "")
    snaf_command = ""
    if snaf_workflow:
        snaf_command = (
            f"bash {script('run_snaf_sample.sh')} --workflow {_q(snaf_workflow)} "
            f"--bam-dir {{outdir}}/rna/star --hla-file {_q(hla_file)} --sample-id {_q(sample_id)} "
            f"--outdir {{outdir}}/branches/splice/snaf && "
            f"test -s {{outdir}}/branches/splice/snaf/snaf_candidates.tsv"
        )
    elif snaf_db:
        snaf_command = (
            f"SNAF_PYTHON={_q(inputs.get('snaf_python') or 'python')} "
            f"NEOAG_ALTANALYZE_IMAGE={_q(inputs.get('altanalyze_image') or 'neoag-altanalyze:snaf')} "
            f"bash {script('run_snaf_pipeline.sh')} "
            f"--bam {{outdir}}/rna/star/Aligned.sortedByCoord.out.bam "
            f"--hla-file {_q(hla_file)} --sample-id {_q(sample_id)} --db-dir {_q(snaf_db)} "
            f"--threads {threads} --outdir {{outdir}}/branches/splice/snaf"
        )
    _stage(lines, "snaf_discovery", required=bool(snaf_command), command=snaf_command,
           outputs={"candidate_table": "{outdir}/branches/splice/snaf/snaf_candidates.tsv"},
           depends_on=["junction_extraction", "rna_expression"])

    bundled_splicemutr = root / "workflows" / "splicemutr" / "SpliceMutr.smk"
    splicemutr_config = str(inputs.get("splicemutr_config") or os.environ.get("SPLICEMUTR_CONFIG") or "")
    splicemutr_samples = str(inputs.get("splicemutr_samples") or os.environ.get("SPLICEMUTR_SAMPLES") or "")
    splicemutr_workflow = str(
        inputs.get("splicemutr_workflow") or os.environ.get("SPLICEMUTR_WORKFLOW")
        or (bundled_splicemutr if splicemutr_config else "")
    )
    splicemutr_command = ""
    if splicemutr_workflow and splicemutr_config:
        samples_override = f" --config sample_sheet={_q(splicemutr_samples)}" if splicemutr_samples else ""
        splicemutr_command = (
            f"NEOAG_SPLICEMUTR_OUTDIR={{outdir}}/branches/splice/splicemutr "
            f"NEOAG_SPLICEMUTR_BAM={{outdir}}/rna/star/Aligned.sortedByCoord.out.bam "
            f"NEOAG_SPLICEMUTR_JUNCTIONS={{outdir}}/branches/splice/regtools_junctions.tsv "
            f"splicemutr-neoag workflow {_q(splicemutr_workflow)} --cores {threads} "
            f"--configfile {_q(splicemutr_config)} --config "
            f"output_dir={{outdir}}/branches/splice/splicemutr target_sample_id={_q(sample_id)}"
            f"{samples_override} && "
            f"test -s {{outdir}}/branches/splice/splicemutr/splicemutr_candidates.tsv"
        )
    elif splicemutr_workflow:
        # Backward-compatible site workflows may carry their own configfile.
        splicemutr_command = (
            f"NEOAG_SPLICEMUTR_OUTDIR={{outdir}}/branches/splice/splicemutr "
            f"NEOAG_SPLICEMUTR_BAM={{outdir}}/rna/star/Aligned.sortedByCoord.out.bam "
            f"NEOAG_SPLICEMUTR_JUNCTIONS={{outdir}}/branches/splice/regtools_junctions.tsv "
            f"splicemutr-neoag workflow {_q(splicemutr_workflow)} --cores {threads} && "
            f"test -s {{outdir}}/branches/splice/splicemutr/splicemutr_candidates.tsv"
        )
    _stage(lines, "splicemutr_discovery", required=bool(splicemutr_command), command=splicemutr_command,
           outputs={"candidate_table": "{outdir}/branches/splice/splicemutr/splicemutr_candidates.tsv"},
           depends_on=["junction_extraction", "rna_expression"])

    fusion_norm = (
        f"PYTHONPATH={_q(root / 'src')} {_q(Path(sys.executable))} -m neoag.cli build-intermediates "
        f"--entry-mode fusion --easyfuse-tsv {{outdir}}/branches/fusion/easyfuse/{_q(sample_id)}/fusions.openneo_union.csv "
        f"--hla-file {_q(hla_file)} --require-nonempty-peptides "
        f"--sample-id {_q(sample_id)} --outdir {{outdir}}/branches/fusion/intermediates"
    )
    _stage(lines, "fusion_peptide_generation", required=True, command=fusion_norm, source="EasyFuse",
           outputs={
               "raw_events": "{outdir}/branches/fusion/intermediates/parsed/raw_events.tsv",
               "raw_peptides": "{outdir}/branches/fusion/intermediates/parsed/raw_peptides.tsv",
           }, depends_on=["fusion_candidate_union"], data_row_outputs=["raw_peptides"])

    splice_norm = (
        f"PYTHONPATH={_q(root / 'src')} {_q(Path(sys.executable))} "
        f"{script('normalize_rna_fusion_splice.py')} --sample-id {_q(sample_id)} "
        f"--profile {_q(PROFILE_NAME)} --genome-build {_q(inputs.get('genome_build') or 'GRCh38')} "
        f"--junctions {{outdir}}/branches/splice/regtools_junctions.tsv "
        f"--annotation-gtf {_q(inputs.get('gencode_gtf') or '')} "
        f"--snaf {{outdir}}/branches/splice/snaf/snaf_candidates.tsv "
        f"--splicemutr {{outdir}}/branches/splice/splicemutr/splicemutr_candidates.tsv "
        f"--normal-junctions {_q(inputs.get('normal_junctions') or '')} "
        f"--outdir {{outdir}}/branches/splice/intermediates"
    )
    splice_norm += (
        f" && PYTHONPATH={_q(root / 'src')} {_q(Path(sys.executable))} "
        f"{script('rebuild_splice_origins_from_splicemutr.py')} --sample-id {_q(sample_id)} "
        f"--genome-build {_q(inputs.get('genome_build') or 'GRCh38')} "
        f"--candidates {{outdir}}/branches/splice/intermediates/raw_peptides.tsv "
        f"--splicemutr-glob '{{outdir}}/branches/splice/splicemutr/formed_transcripts/**/*_data_splicemutr_cp_corrected.txt' "
        f"--outdir {{outdir}}/branches/splice/intermediates/formal_origins"
    )
    _stage(lines, "splice_candidate_normalization", required=bool(snaf_command and splicemutr_command), command=splice_norm, source="splice_consensus",
           outputs={
               "raw_events": "{outdir}/branches/splice/intermediates/raw_events.tsv",
               "raw_peptides": "{outdir}/branches/splice/intermediates/formal_origins/raw_peptides.formal_origins.tsv",
               "rna_junction_tsv": "{outdir}/branches/splice/intermediates/rna_junction_evidence.tsv",
               "splice_junctions": "{outdir}/branches/splice/intermediates/splice_junctions.tsv",
               "splice_tool_evidence": "{outdir}/branches/splice/intermediates/splice_tool_evidence.long.tsv",
               "splice_peptide_provenance": "{outdir}/branches/splice/intermediates/splice_peptide_provenance.tsv",
               "formal_origin_summary": "{outdir}/branches/splice/intermediates/formal_origins/rebuild_summary.json",
               "splice_event_merge_provenance": "{outdir}/branches/splice/intermediates/splice_event_merge_provenance.tsv",
               "splice_peptide_merge_provenance": "{outdir}/branches/splice/intermediates/splice_peptide_merge_provenance.tsv",
               "splice_merge_conflicts": "{outdir}/branches/splice/intermediates/splice_merge_conflicts.tsv",
               "splice_consensus": "{outdir}/branches/splice/intermediates/splice_consensus.tsv",
               "splice_consensus_provenance": "{outdir}/branches/splice/intermediates/splice_consensus_provenance.tsv",
               "splice_consensus_conflicts": "{outdir}/branches/splice/intermediates/splice_consensus_conflicts.tsv",
               "junction_aliases": "{outdir}/branches/splice/intermediates/junction_aliases.tsv",
               "splice_evidence_conflicts": "{outdir}/branches/splice/intermediates/evidence_conflicts.tsv",
               "splice_qc": "{outdir}/branches/splice/intermediates/splice_qc.tsv",
               "provenance_manifest": "{outdir}/branches/splice/intermediates/provenance_manifest.json",
           }, depends_on=["junction_extraction", "snaf_discovery", "splicemutr_discovery"])

    lines += [
        "",
        "[evidence]",
        'expression = "{outdir}/rna/expression/gene_tpm.tsv"',
    ]
    if salmon_primary and rsem_available and quant_method == "auto":
        lines += [
            'rsem_expression = "{outdir}/rna/rsem_expression/gene_tpm.tsv"',
            'rsem_transcript_expression = "{outdir}/rna/rsem_expression/transcript_tpm.tsv"',
        ]
    lines += [
        'rna_junction_tsv = "{outdir}/branches/splice/intermediates/rna_junction_evidence.tsv"',
        f"normal_junctions = {_toml(str(inputs.get('normal_junctions') or ''))}",
        f"normal_expression = {_toml(str(inputs.get('normal_expression') or ''))}",
        f"normal_hla_ligands = {_toml(str(inputs.get('normal_hla_ligands') or ''))}",
        f"reference_proteome = {_toml(str(inputs.get('reference_proteome') or ''))}",
    ]
    manifest_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return {
        "profile": PROFILE_NAME,
        "manifest": str(manifest_path),
        "requirements": [row.__dict__ for row in requirements],
        "missing_required": missing_required,
        "ready_for_execute": not missing_required,
    }
