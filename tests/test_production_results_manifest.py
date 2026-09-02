from __future__ import annotations

import json
import subprocess
import sys
import tomllib
from pathlib import Path

from neoag.utils import read_tsv
from neoag.production_runner import _materialize_reusable_immunogenicity_outputs, _outputs_ready


ROOT = Path(__file__).resolve().parents[1]


def write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def make_star_index(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    for name in ("Genome", "SA", "SAindex"):
        write(path / name, "fixture\n")
    write(path / "genomeParameters.txt", "versionGenome 2.7.11b\nsjdbOverhang 149\n")
    return path


def test_outputs_ready_rejects_existing_table_with_blank_required_fields(tmp_path):
    purity = write(
        tmp_path / "recommended_purity.tsv",
        "sample_id\tpurity\tploidy\nS1\t\t\n",
    )
    outputs = {"purity": str(purity)}
    assert not _outputs_ready(
        outputs,
        data_row_outputs=["purity"],
        required_output_fields=["purity:purity", "purity:ploidy"],
    )
    purity.write_text("sample_id\tpurity\tploidy\nS1\t0.61\t2.0\n", encoding="utf-8")
    assert _outputs_ready(
        outputs,
        data_row_outputs=["purity"],
        required_output_fields=["purity:purity", "purity:ploidy"],
    )


def test_generator_builds_all_three_upstream_consensus_stages(tmp_path):
    optitype = write(tmp_path / "optitype/result.tsv", "A1\tA2\tB1\tB2\tC1\tC2\nA*02:01\tA*11:01\tB*15:01\tB*40:01\tC*03:04\tC*08:01\n").parent
    spechla = write(tmp_path / "spechla/hla.result.txt", "A\tA*02:01\tA*11:01\nB\tB*15:01\tB*40:01\nC\tC*03:04\tC*08:01\n").parent
    facets = write(tmp_path / "facets/facets_purity.txt", "purity\t0.60\n").parent
    sequenza = write(tmp_path / "sequenza/result.tsv", "purity\tploidy\n0.62\t2.1\n").parent
    purple = write(tmp_path / "purple/purity.tsv", "purity\tploidy\n0.61\t2.0\n").parent
    lohhla = write(tmp_path / "lohhla/hla_loh.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    spechla_loh = write(tmp_path / "spechla_loh/hla_loh.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    star = write(tmp_path / "star/star-fusion.fusion_predictions.tsv", "#FusionName\tJunctionReadCount\nEWSR1--WT1\t20\n")
    output = tmp_path / "production.toml"
    command = [sys.executable, str(ROOT / "scripts/generate_production_from_results_manifest.py"), "--project-root", str(ROOT), "--sample-id", "S1", "--outdir", str(tmp_path / "run"), "--output", str(output), "--optitype", str(optitype), "--spechla-typing", str(spechla), "--facets", str(facets), "--sequenza", str(sequenza), "--purple", str(purple), "--lohhla", str(lohhla), "--spechla-loh", str(spechla_loh), "--star-fusion", str(star)]
    subprocess.run(command, check=True)
    manifest = tomllib.loads(output.read_text(encoding="utf-8"))
    assert manifest["run"]["hla_file"].endswith("evidence/hla_typing/recommended_hla.txt")
    assert {"hla_typing_consensus", "purity_cnv_consensus", "hla_loh_consensus", "fusion_candidates"} <= set(manifest["stages"])
    assert "--jaffal" not in manifest["stages"]["fusion_candidates"]["command"]
    assert manifest["evidence"]["purity"].endswith("evidence/purity_cnv/recommended_purity.tsv")
    assert manifest["evidence"]["hla_loh"].endswith("evidence/hla_loh/hla_loh_consensus.tsv")
    assert manifest["run"]["profile"].endswith("profiles/sarcoma_rna_supported_v2_provisional.toml")
    assert manifest["evidence"]["evidence_consensus_rules"].endswith(
        "configs/ranking/sarcoma_evidence_consensus_v3_source_chain.toml"
    )


def test_generator_and_runner_reuse_immunogenicity_evidence_independently(tmp_path):
    hla = write(tmp_path / "hla.txt", "HLA-A*02:01\n")
    purity = write(tmp_path / "purity.tsv", "sample_id\tpurity\tploidy\nS1\t0.60\t2.0\n")
    cnv = write(tmp_path / "cnv.tsv", "chrom\tstart\tend\ttotal_cn\nchr1\t1\t1000\t2\n")
    lohhla = write(tmp_path / "lohhla.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    spechla_loh = write(tmp_path / "spechla_loh.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    vcf = write(tmp_path / "somatic.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\t.\tA\tT\t.\tPASS\t.\n")
    prime = write(tmp_path / "existing/prime_evidence.tsv", "peptide\thla_allele\tprime_score\nSIINFEKL\tHLA-A*02:01\t0.8\n")
    bigmhc = write(tmp_path / "existing/bigmhc_im_evidence.tsv", "peptide\thla_allele\tbigmhc_im_score\nSIINFEKL\tHLA-A*02:01\t0.7\n")
    deep = write(tmp_path / "existing/deepimmuno_evidence.tsv", "peptide\thla_allele\tdeepimmuno_score\nSIINFEKL\tHLA-A*02:01\t0.6\n")
    output = tmp_path / "production.toml"
    subprocess.run([
        sys.executable, str(ROOT / "scripts/generate_production_from_results_manifest.py"),
        "--project-root", str(ROOT), "--sample-id", "S1", "--outdir", str(tmp_path / "run"),
        "--output", str(output), "--hla-file", str(hla), "--purity", str(purity), "--cnv", str(cnv),
        "--lohhla", str(lohhla), "--spechla-loh", str(spechla_loh), "--somatic-vcf", str(vcf),
        "--prime-evidence", str(prime), "--bigmhc-evidence", str(bigmhc),
        "--deepimmuno-evidence", str(deep),
    ], check=True)
    manifest = tomllib.loads(output.read_text(encoding="utf-8"))
    assert manifest["evidence"]["prime_evidence"] == str(prime)
    assert manifest["evidence"]["bigmhc_evidence"] == str(bigmhc)
    assert manifest["evidence"]["deepimmuno_evidence"] == str(deep)

    final = tmp_path / "run/final"
    reused = _materialize_reusable_immunogenicity_outputs(final, manifest["evidence"])
    assert reused == ["prime", "bigmhc_im", "deepimmuno"]
    assert (final / "presentation/prime_evidence.tsv").read_text() == prime.read_text()
    assert (final / "presentation/bigmhc_im_evidence.tsv").read_text() == bigmhc.read_text()
    assert (final / "presentation/deepimmuno_evidence.tsv").read_text() == deep.read_text()


def test_generator_propagates_paired_vcf_sample_roles(tmp_path):
    hla = write(tmp_path / "hla.txt", "HLA-A*02:01\n")
    purity = write(tmp_path / "purity.tsv", "sample_id\tpurity\tploidy\nS1\t0.60\t2.0\n")
    cnv = write(tmp_path / "cnv.tsv", "chrom\tstart\tend\ttotal_cn\nchr1\t1\t1000\t2\n")
    lohhla = write(tmp_path / "lohhla.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    spechla_loh = write(tmp_path / "spechla_loh.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    vcf = write(tmp_path / "somatic.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tblood\ttumor\n")
    output = tmp_path / "production.toml"
    subprocess.run([
        sys.executable, str(ROOT / "scripts/generate_production_from_results_manifest.py"),
        "--project-root", str(ROOT), "--sample-id", "S1", "--outdir", str(tmp_path / "run"),
        "--output", str(output), "--hla-file", str(hla), "--purity", str(purity), "--cnv", str(cnv),
        "--lohhla", str(lohhla), "--spechla-loh", str(spechla_loh), "--somatic-vcf", str(vcf),
        "--tumor-sample-name", "tumor", "--normal-sample-name", "blood",
    ], check=True)
    command = tomllib.loads(output.read_text(encoding="utf-8"))["stages"]["snv_indel_candidates"]["command"]
    assert '--tumor-sample-name "tumor"' in command
    assert '--normal-sample-name "blood"' in command


def test_generator_filters_splice_candidates_before_production_scoring(tmp_path):
    hla = write(tmp_path / "hla.txt", "HLA-A*02:01\n")
    facets = write(tmp_path / "facets/facets_purity.txt", "purity\t0.60\n").parent
    ascat = write(tmp_path / "ascat/ascat_summary.tsv", "sample_id\tpurity\tploidy\nS1\t0.64\t2.1\n").parent
    lohhla = write(tmp_path / "lohhla/hla_loh.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    spechla_loh = write(tmp_path / "spechla_loh/hla_loh.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    junctions = write(tmp_path / "junctions.tsv", "chrom\tstart\tend\tstrand\tjunction_reads\nchr1\t10\t20\t+\t12\n")
    snaf = write(tmp_path / "snaf.tsv", "event_id\tpeptide\thla_allele\tbinding_rank\nE1\tACDEFGHIK\tHLA-A*02:01\t0.8\n")
    output = tmp_path / "production.toml"
    subprocess.run([sys.executable, str(ROOT / "scripts/generate_production_from_results_manifest.py"), "--project-root", str(ROOT), "--sample-id", "S1", "--outdir", str(tmp_path / "run"), "--output", str(output), "--hla-file", str(hla), "--facets", str(facets), "--ascat", str(ascat), "--lohhla", str(lohhla), "--spechla-loh", str(spechla_loh), "--junctions", str(junctions), "--snaf", str(snaf)], check=True)
    manifest = tomllib.loads(output.read_text(encoding="utf-8"))
    splice = manifest["stages"]["splice_candidates"]
    assert "filter_splice_production_candidates.py" in splice["command"]
    assert splice["outputs"]["raw_peptides"].endswith("production_selected/raw_peptides.tsv")


def test_generator_enriches_splice_qc_from_star_and_rna_bam(tmp_path):
    hla = write(tmp_path / "hla.txt", "HLA-A*02:01\n")
    facets = write(tmp_path / "facets/facets_purity.txt", "purity\t0.60\n").parent
    ascat = write(tmp_path / "ascat/ascat_summary.tsv", "sample_id\tpurity\tploidy\nS1\t0.64\t2.1\n").parent
    lohhla = write(tmp_path / "lohhla/hla_loh.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    spechla_loh = write(tmp_path / "spechla_loh/hla_loh.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    star_sj = write(tmp_path / "star/SJ.out.tab", "chr1\t10\t20\t1\t1\t0\t12\t0\t25\n")
    rna_bam = write(tmp_path / "star/Aligned.sortedByCoord.out.bam", "fixture\n")
    write(Path(str(rna_bam) + ".bai"), "fixture\n")
    rna_vaf = write(tmp_path / "rna_alt_vaf.tsv", "chrom\tpos\trna_depth\nchr1\t10\t20\n")
    snaf = write(tmp_path / "snaf.tsv", "event_id\tpeptide\thla_allele\tbinding_rank\nE1\tACDEFGHIK\tHLA-A*02:01\t0.8\n")
    normal = write(tmp_path / "normal_junctions.tsv", "junction_id\nchr1:10-20:+\n")
    write(Path(str(normal) + ".sqlite"), "fixture\n")
    output = tmp_path / "production.toml"
    subprocess.run([
        sys.executable, str(ROOT / "scripts/generate_production_from_results_manifest.py"),
        "--project-root", str(ROOT), "--sample-id", "S1", "--outdir", str(tmp_path / "run"),
        "--output", str(output), "--hla-file", str(hla), "--facets", str(facets),
        "--ascat", str(ascat), "--lohhla", str(lohhla), "--spechla-loh", str(spechla_loh),
        "--star-sj", str(star_sj), "--snaf", str(snaf), "--rna-bam", str(rna_bam),
        "--rna-vaf", str(rna_vaf), "--normal-junctions", str(normal),
    ], check=True)
    manifest = tomllib.loads(output.read_text(encoding="utf-8"))
    splice = manifest["stages"]["splice_candidates"]
    assert "build_splice_junction_qc_from_star_bam.py" in splice["command"]
    assert "--splice-consensus" in splice["command"]
    assert "star_bam_qc/raw_events.enriched.tsv" in splice["command"]
    assert splice["outputs"]["junction_read_qc"].endswith("splice_junction_qc.enriched.tsv")
    assert "rna_bam_input" in splice["depends_on"]

    output_without_normal = tmp_path / "production-without-normal.toml"
    subprocess.run([
        sys.executable, str(ROOT / "scripts/generate_production_from_results_manifest.py"),
        "--project-root", str(ROOT), "--sample-id", "S1", "--outdir", str(tmp_path / "run-no-normal"),
        "--output", str(output_without_normal), "--hla-file", str(hla), "--facets", str(facets),
        "--ascat", str(ascat), "--lohhla", str(lohhla), "--spechla-loh", str(spechla_loh),
        "--star-sj", str(star_sj), "--snaf", str(snaf), "--rna-bam", str(rna_bam),
        "--rna-vaf", str(rna_vaf),
    ], check=True)
    command_without_normal = tomllib.loads(output_without_normal.read_text(encoding="utf-8"))["stages"]["splice_candidates"]["command"]
    assert "build_splice_junction_qc_from_star_bam.py" in command_without_normal
    assert "--normal-junction-sqlite" not in command_without_normal


def test_generator_accepts_facets_and_ascat_when_other_purity_tools_failed(tmp_path):
    hla = write(tmp_path / "hla.txt", "HLA-A*02:01\n")
    fasta = write(tmp_path / "GRCh38.fa", ">chr1\nACGT\n")
    facets = write(tmp_path / "facets/facets_purity.txt", "purity\t0.60\n").parent
    ascat = write(tmp_path / "ascat/ascat_summary.tsv", "sample_id\tpurity\tploidy\nS1\t0.64\t2.1\n").parent
    lohhla = write(tmp_path / "lohhla/hla_loh.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    spechla_loh = write(tmp_path / "spechla_loh/hla_loh.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    star = write(tmp_path / "star/star-fusion.fusion_predictions.tsv", "#FusionName\tJunctionReadCount\nEWSR1--WT1\t20\n")
    output = tmp_path / "production.toml"
    command = [sys.executable, str(ROOT / "scripts/generate_production_from_results_manifest.py"), "--project-root", str(ROOT), "--sample-id", "S1", "--outdir", str(tmp_path / "run"), "--output", str(output), "--hla-file", str(hla), "--reference-fasta", str(fasta), "--facets", str(facets), "--ascat", str(ascat), "--lohhla", str(lohhla), "--spechla-loh", str(spechla_loh), "--somatic-vcf", str(star), "--star-fusion", str(star)]
    subprocess.run(command, check=True)
    manifest = tomllib.loads(output.read_text(encoding="utf-8"))
    assert manifest["run"]["required_tool_groups"]["purity_cnv"]["tools"] == ["facets", "ascat"]
    assert "purity_facets" in manifest["stages"]
    assert "purity_ascat" in manifest["stages"]
    assert "purity_sequenza" not in manifest["stages"]
    assert "purity_purple" not in manifest["stages"]
    assert manifest["run"]["profile"].endswith("profiles/sarcoma_rna_supported_v2_provisional.toml")
    assert manifest["evidence"]["evidence_consensus_rules"].endswith(
        "configs/ranking/sarcoma_evidence_consensus_v3_source_chain.toml"
    )
    assert f'NEOAG_REFERENCE_FASTA="{fasta}"' in manifest["stages"]["snv_indel_candidates"]["command"]


def test_generator_builds_star_and_rna_allele_count_stages_from_fastq(tmp_path):
    hla = write(tmp_path / "hla.txt", "HLA-A*02:01\n")
    purity = write(tmp_path / "purity.tsv", "sample_id\tpurity\tploidy\nS1\t0.60\t2.0\n")
    cnv = write(tmp_path / "cnv.tsv", "chrom\tstart\tend\ttotal_cn\nchr1\t1\t1000\t2\n")
    lohhla = write(tmp_path / "lohhla.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    spechla_loh = write(tmp_path / "spechla_loh.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    vcf = write(tmp_path / "somatic.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\t.\tA\tT\t.\tPASS\t.\n")
    fastq1 = write(tmp_path / "rna_R1.fastq", "@r1\nACGT\n+\nFFFF\n")
    fastq2 = write(tmp_path / "rna_R2.fastq", "@r1\nTGCA\n+\nFFFF\n")
    star_index = make_star_index(tmp_path / "star_index")
    gtf = write(tmp_path / "gencode.gtf", 'chr1\ttest\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG1";\n')
    output = tmp_path / "production.toml"
    subprocess.run([
        sys.executable, str(ROOT / "scripts/generate_production_from_results_manifest.py"),
        "--project-root", str(ROOT), "--sample-id", "S1",
        "--outdir", str(tmp_path / "run"), "--output", str(output),
        "--hla-file", str(hla), "--purity", str(purity), "--cnv", str(cnv),
        "--lohhla", str(lohhla), "--spechla-loh", str(spechla_loh),
        "--somatic-vcf", str(vcf), "--rna-fastq1", str(fastq1),
        "--rna-fastq2", str(fastq2), "--star-index", str(star_index),
        "--gencode-gtf", str(gtf), "--rna-threads", "8",
    ], check=True)
    manifest = tomllib.loads(output.read_text(encoding="utf-8"))
    assert "rna_star_alignment" in manifest["stages"]
    assert manifest["stages"]["rna_star_index"]["source"] == "STAR_INDEX_REUSE_EXPLICIT"
    assert manifest["stages"]["rna_star_alignment"]["depends_on"] == ["rna_star_index"]
    assert "run_star_rna_fastq.sh" in manifest["stages"]["rna_star_alignment"]["command"]
    assert manifest["stages"]["rna_alt_vaf"]["depends_on"] == ["rna_star_alignment"]
    assert "rna_allele_counts_pysam.py" in manifest["stages"]["rna_alt_vaf"]["command"]
    assert manifest["evidence"]["rna_vaf"] == "{outdir}/rna/rna_alt_vaf.tsv"


def test_generator_reuses_valid_easyfuse_star_index(tmp_path):
    hla = write(tmp_path / "hla.txt", "HLA-A*02:01\n")
    purity = write(tmp_path / "purity.tsv", "sample_id\tpurity\tploidy\nS1\t0.60\t2.0\n")
    cnv = write(tmp_path / "cnv.tsv", "chrom\tstart\tend\ttotal_cn\nchr1\t1\t1000\t2\n")
    lohhla = write(tmp_path / "lohhla.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    spechla_loh = write(tmp_path / "spechla_loh.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    vcf = write(tmp_path / "somatic.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    fastq1 = write(tmp_path / "rna_R1.fastq", "@r1\nACGT\n+\nFFFF\n")
    fastq2 = write(tmp_path / "rna_R2.fastq", "@r1\nTGCA\n+\nFFFF\n")
    easyfuse_index = make_star_index(tmp_path / "easyfuse/ref_genome.fa.star.idx")
    gtf = write(tmp_path / "gencode.gtf", 'chr1\ttest\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG1";\n')
    output = tmp_path / "manifest/production.toml"

    subprocess.run([
        sys.executable, str(ROOT / "scripts/generate_production_from_results_manifest.py"),
        "--project-root", str(ROOT), "--sample-id", "S1", "--outdir", str(tmp_path / "run"),
        "--output", str(output), "--hla-file", str(hla), "--purity", str(purity), "--cnv", str(cnv),
        "--lohhla", str(lohhla), "--spechla-loh", str(spechla_loh), "--somatic-vcf", str(vcf),
        "--rna-fastq1", str(fastq1), "--rna-fastq2", str(fastq2),
        "--easyfuse-star-index", str(easyfuse_index), "--gencode-gtf", str(gtf),
    ], check=True)

    manifest = tomllib.loads(output.read_text(encoding="utf-8"))
    index_stage = manifest["stages"]["rna_star_index"]
    assert index_stage["source"] == "STAR_INDEX_REUSE_EASYFUSE"
    assert index_stage["outputs"]["star_index"] == str(easyfuse_index.resolve())
    validation = json.loads((output.parent / "star_index_validation.json").read_text())
    assert validation["status"] == "VALIDATED_REUSE"
    assert validation["selected_source"] == "EASYFUSE"


def test_generator_plans_star_index_rebuild_when_reuse_is_invalid(tmp_path):
    hla = write(tmp_path / "hla.txt", "HLA-A*02:01\n")
    purity = write(tmp_path / "purity.tsv", "sample_id\tpurity\tploidy\nS1\t0.60\t2.0\n")
    cnv = write(tmp_path / "cnv.tsv", "chrom\tstart\tend\ttotal_cn\nchr1\t1\t1000\t2\n")
    lohhla = write(tmp_path / "lohhla.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    spechla_loh = write(tmp_path / "spechla_loh.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    vcf = write(tmp_path / "somatic.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    fastq1 = write(tmp_path / "rna_R1.fastq", "@r1\nACGT\n+\nFFFF\n")
    fastq2 = write(tmp_path / "rna_R2.fastq", "@r1\nTGCA\n+\nFFFF\n")
    invalid_index = tmp_path / "easyfuse/incomplete"
    invalid_index.mkdir(parents=True)
    fasta = write(tmp_path / "GRCh38.fa", ">chr1\nACGT\n")
    gtf = write(tmp_path / "gencode.gtf", 'chr1\ttest\tgene\t1\t4\t.\t+\t.\tgene_id "ENSG1";\n')
    output = tmp_path / "manifest/production.toml"

    subprocess.run([
        sys.executable, str(ROOT / "scripts/generate_production_from_results_manifest.py"),
        "--project-root", str(ROOT), "--sample-id", "S1", "--outdir", str(tmp_path / "run"),
        "--output", str(output), "--hla-file", str(hla), "--purity", str(purity), "--cnv", str(cnv),
        "--lohhla", str(lohhla), "--spechla-loh", str(spechla_loh), "--somatic-vcf", str(vcf),
        "--rna-fastq1", str(fastq1), "--rna-fastq2", str(fastq2), "--reference-fasta", str(fasta),
        "--easyfuse-star-index", str(invalid_index), "--gencode-gtf", str(gtf),
    ], check=True)

    manifest = tomllib.loads(output.read_text(encoding="utf-8"))
    index_stage = manifest["stages"]["rna_star_index"]
    assert index_stage["source"] == "STAR_INDEX_BUILD"
    assert "build_star_index.sh" in index_stage["command"]
    validation = json.loads((output.parent / "star_index_validation.json").read_text())
    assert validation["status"] == "PLANNED_BUILD"


def test_generator_reuses_existing_rna_vaf_without_star(tmp_path):
    hla = write(tmp_path / "hla.txt", "HLA-A*02:01\n")
    purity = write(tmp_path / "purity.tsv", "sample_id\tpurity\tploidy\nS1\t0.60\t2.0\n")
    cnv = write(tmp_path / "cnv.tsv", "chrom\tstart\tend\ttotal_cn\nchr1\t1\t1000\t2\n")
    lohhla = write(tmp_path / "lohhla.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    spechla_loh = write(tmp_path / "spechla_loh.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    vcf = write(tmp_path / "somatic.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    rna_vaf = write(tmp_path / "rna_alt_vaf.tsv", "chrom\tpos\tref\talt\trna_depth\trna_alt_reads\trna_vaf\nchr1\t10\tA\tT\t20\t5\t0.25\n")
    output = tmp_path / "production.toml"
    subprocess.run([
        sys.executable, str(ROOT / "scripts/generate_production_from_results_manifest.py"),
        "--project-root", str(ROOT), "--sample-id", "S1",
        "--outdir", str(tmp_path / "run"), "--output", str(output),
        "--hla-file", str(hla), "--purity", str(purity), "--cnv", str(cnv),
        "--lohhla", str(lohhla), "--spechla-loh", str(spechla_loh),
        "--somatic-vcf", str(vcf), "--rna-vaf", str(rna_vaf),
    ], check=True)
    manifest = tomllib.loads(output.read_text(encoding="utf-8"))
    assert "rna_alt_vaf_input" in manifest["stages"]
    assert "rna_star_alignment" not in manifest["stages"]
    assert manifest["evidence"]["rna_vaf"] == str(rna_vaf.resolve())


def test_fusion_union_keeps_event_only_callers_and_provided_peptides(tmp_path):
    hla = write(tmp_path / "hla.txt", "HLA-A*02:01\n")
    star = write(tmp_path / "star.tsv", "#FusionName\tLeftBreakpoint\tRightBreakpoint\tJunctionReadCount\nEWSR1--WT1\tchr22:1:+\tchr11:2:-\t20\n")
    arriba = write(tmp_path / "arriba.tsv", "gene1\tgene2\tbreakpoint1\tbreakpoint2\tjunction_peptide\tsplit_reads1\nCPSF6\tRARG\tchr12:3:+\tchr12:4:-\tQYYSTPWTF\t12\n")
    outdir = tmp_path / "union"
    subprocess.run([sys.executable, str(ROOT / "scripts/build_fusion_caller_union.py"), "--sample-id", "S1", "--hla-file", str(hla), "--star-fusion", str(star), "--arriba", str(arriba), "--outdir", str(outdir)], cwd=ROOT, check=True)
    events = read_tsv(outdir / "raw_events.tsv")
    peptides = read_tsv(outdir / "raw_peptides.tsv")
    audit = read_tsv(outdir / "fusion_caller_union.tsv")
    assert {row["gene"] for row in events} == {"EWSR1::WT1", "CPSF6::RARG"}
    assert any(row["peptide"] == "QYYSTPWTF" and row["hla_allele"] == "HLA-A*02:01" for row in peptides)
    assert any(row["gene_pair"] == "EWSR1::WT1" and row["peptide_status"] == "ORF_PEPTIDE_UNAVAILABLE_REVIEW_ONLY" for row in audit)


def test_fusion_union_rescues_only_exact_diagnostic_whitelist_from_unfiltered_easyfuse(tmp_path):
    hla = write(tmp_path / "hla.txt", "HLA-A*02:01\n")
    header = (
        "BPID;Fusion_Gene;Breakpoint1;Breakpoint2;FTID;prediction_class;prediction_prob;"
        "fusioncatcher_detected;star_detected;arriba_detected;tool_count;fusioncatcher_junc;"
        "fusioncatcher_span;fusioncatcher_anch;frame;type;neo_peptide_sequence;neo_peptide_sequence_bp"
    )
    unfiltered = write(
        tmp_path / "fusions.csv",
        header + "\n"
        "ews;EWSR1_WT1;22:100:+;11:200:-;tx1;;0;1;0;0;1;20;10;18;in_frame;trans;ACDEFGHIKLMNPQRSTVWY;10\n"
        "noise;NOISE_BACKGROUND;1:10:+;2:20:-;tx2;;0;1;0;0;1;30;20;18;in_frame;trans;ACDEFGHIKLMNPQRSTVWY;10\n",
    )
    outdir = tmp_path / "union-rescue"
    subprocess.run([
        sys.executable, str(ROOT / "scripts/build_fusion_caller_union.py"),
        "--sample-id", "S1", "--hla-file", str(hla),
        "--easyfuse-unfiltered", str(unfiltered), "--outdir", str(outdir),
    ], cwd=ROOT, check=True)
    events = read_tsv(outdir / "raw_events.tsv")
    peptides = read_tsv(outdir / "raw_peptides.tsv")
    audit = read_tsv(outdir / "fusion_caller_union.tsv")
    rescue = read_tsv(outdir / "diagnostic_fusion_rescue.tsv")
    assert {row["gene"] for row in events} == {"EWSR1::WT1"}
    assert peptides and all(row["crosses_junction"] == "yes" for row in peptides)
    assert all(row["gene"] == "EWSR1::WT1" for row in peptides)
    assert audit[0]["admission_policy"] == "DIAGNOSTIC_WHITELIST_RESCUE"
    assert audit[0]["rescue_reason"] == "diagnostic_whitelist_not_in_easyfuse_pass"
    assert {row["fusion_gene"] for row in rescue} == {"EWSR1::WT1"}


def test_easyfuse_audit_joins_source_by_event_id_not_filtered_row_order(tmp_path):
    hla = write(tmp_path / "hla.txt", "HLA-A*02:01\n")
    header = (
        "BPID;Fusion_Gene;Breakpoint1;Breakpoint2;FTID;prediction_class;prediction_prob;"
        "fusioncatcher_detected;star_detected;arriba_detected;fusioncatcher_junc;"
        "fusioncatcher_span;fusioncatcher_anch;frame;type;neo_peptide_sequence;neo_peptide_sequence_bp"
    )
    easyfuse = write(
        tmp_path / "fusions.pass.csv",
        header + "\n"
        "noise;NOISE_EVENT;10:10:+;14:20:-;noise_tx;negative;0.1;1;0;0;20;10;30;in_frame;trans;ACDEFGHIKLMNPQRST;8\n"
        "ews;EWSR1_WT1;22:100:+;11:200:-;ews_tx;positive;0.9;1;1;0;20;10;30;in_frame;trans;ACDEFGHIKLMNPQRST;8\n",
    )
    outdir = tmp_path / "union-order"
    subprocess.run([
        sys.executable, str(ROOT / "scripts/build_fusion_caller_union.py"),
        "--sample-id", "S1", "--hla-file", str(hla), "--easyfuse", str(easyfuse),
        "--disable-diagnostic-fusion-rescue", "--no-targeted-fusion-rescue",
        "--outdir", str(outdir),
    ], cwd=ROOT, check=True)
    events = read_tsv(outdir / "raw_events.tsv")
    assert len(events) == 1
    assert events[0]["gene"] == "EWSR1::WT1"
    assert events[0]["breakpoint1"] == "22:100:+"
    assert events[0]["breakpoint2"] == "11:200:-"
    assert "10:10" not in events[0]["adjacency_key"]


def test_confirmed_expressed_product_generates_closed_fusion_origin_chain(tmp_path):
    hla = write(tmp_path / "hla.txt", "HLA-A*02:01\n")
    header = (
        "BPID;Fusion_Gene;Breakpoint1;Breakpoint2;FTID;prediction_class;prediction_prob;"
        "fusioncatcher_detected;star_detected;arriba_detected;fusioncatcher_junc;"
        "fusioncatcher_span;fusioncatcher_anch;frame;type;neo_peptide_sequence;neo_peptide_sequence_bp"
    )
    easyfuse = write(
        tmp_path / "fusions.pass.csv",
        header + "\n"
        "ews;EWSR1_WT1;22:100:+;11:200:-;caller_tx;positive;0.9;1;1;1;20;10;30;in_frame;trans;ACDEFGHIKLMNPQRST;8\n",
    )
    junctions = write(
        tmp_path / "Chimeric.out.junction",
        "chr22\t100\t+\tchr11\t200\t-\t1\t0\t0\tread1\n"
        "chr22\t100\t+\tchr11\t200\t-\t1\t0\t0\tread2\n"
        "chr22\t100\t+\tchr11\t200\t-\t1\t0\t0\tread3\n",
    )
    products = write(
        tmp_path / "expressed_products.tsv",
        "genome_build\tchrom1\tpos1\tstrand1\tchrom2\tpos2\tstrand2\tgene1\tgene2\ttranscript1\ttranscript2\tprotein_sequence\twildtype_protein_sequence\tjunction_aa_position\tin_frame\torf_status\tnmd_status\tsource_tool\tsource_record_id\n"
        "GRCh38\t22\t100\t+\t11\t200\t-\tEWSR1\tWT1\tENST1\tENST2\tMACDEFGHIKLMNPQRSTVWY\t\t10\tyes\tCONFIRMED\tNOT_AT_RISK\tAGFusion\tAGF1\n",
    )
    outdir = tmp_path / "union-confirmed-product"
    subprocess.run([
        sys.executable, str(ROOT / "scripts/build_fusion_caller_union.py"),
        "--sample-id", "S1", "--hla-file", str(hla), "--easyfuse", str(easyfuse),
        "--star-chimeric", str(junctions), "--fusion-expressed-products", str(products),
        "--disable-diagnostic-fusion-rescue", "--no-targeted-fusion-rescue",
        "--outdir", str(outdir),
    ], cwd=ROOT, check=True)
    origin = read_tsv(outdir / "fusion_peptide_origin_chain.tsv")
    closed = [row for row in origin if row["source_chain_status"] == "CLOSED_ORF_AND_EXACT_JUNCTION"]
    assert closed
    assert all(row["orf_id"] and "|" in row["fusion_junction_display"] for row in closed)
    assert all(row["exact_verified_junction_reads"] == "3" for row in closed)
    queue = read_tsv(outdir / "fusion_orf_completion_queue.tsv")
    assert queue[0]["orf_completion_status"] == "ORF_SOURCE_CHAIN_COMPLETE"


def test_generator_unions_existing_fusion_callers_including_jaffal(tmp_path):
    hla = write(tmp_path / "hla.txt", "HLA-A*02:01\n")
    purity = write(tmp_path / "purity.tsv", "sample_id\tpurity\tploidy\nS1\t0.60\t2.0\n")
    cnv = write(tmp_path / "cnv.tsv", "chrom\tstart\tend\ttotal_cn\nchr1\t1\t1000\t2\n")
    lohhla = write(tmp_path / "lohhla.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    spechla_loh = write(tmp_path / "spechla_loh.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    jaffal = write(tmp_path / "long-rna/jaffal/output/jaffa_results.csv", "fusion genes,chrom1,base1,chrom2,base2,spanning reads,classification\nEWSR1::WT1,chr22,100,chr11,200,15,HighConfidence\nNOISE::EVENT,chr1,10,chr2,20,2,PotentialTransSplicing\n")
    caller_root = tmp_path / "long-rna"
    stabpan = write(
        tmp_path / "netmhcstabpan_evidence.tsv",
        "sample_id\tpeptide\thla_allele\tnetmhcstabpan_score\tnetmhcstabpan_rank\n"
        "S1\tAAAAAAAAA\tHLA-A*02:01\t0.8\t0.7\n",
    )
    output = tmp_path / "production.toml"
    subprocess.run([
        sys.executable, str(ROOT / "scripts/generate_production_from_results_manifest.py"),
        "--project-root", str(ROOT), "--sample-id", "S1", "--outdir", str(tmp_path / "run"),
        "--output", str(output), "--hla-file", str(hla), "--purity", str(purity), "--cnv", str(cnv),
        "--lohhla", str(lohhla), "--spechla-loh", str(spechla_loh),
        "--jaffal", str(jaffal), "--fusion-caller-root", str(caller_root),
        "--skip-netmhcstabpan", "--netmhcstabpan-evidence", str(stabpan),
    ], check=True)
    manifest = tomllib.loads(output.read_text(encoding="utf-8"))
    assert manifest["evidence"]["netmhcstabpan"] == str(stabpan)
    assert "NETMHCSTABPAN_LOCAL_UNAVAILABLE" not in manifest["run"].get("production_limitations", [])
    fusion = manifest["stages"]["fusion_candidates"]
    assert "build_fusion_caller_union.py" in fusion["command"]
    assert "--jaffal" in fusion["command"]
    assert "--caller-root" in fusion["command"]
    assert fusion["outputs"]["fusion_consensus"].endswith("fusion_consensus.tsv")

    union_dir = tmp_path / "union"
    subprocess.run([
        sys.executable, str(ROOT / "scripts/build_fusion_caller_union.py"),
        "--sample-id", "S1", "--hla-file", str(hla), "--jaffal", str(jaffal),
        "--outdir", str(union_dir),
    ], cwd=ROOT, check=True)
    events = read_tsv(union_dir / "raw_events.tsv")
    assert len(events) == 1
    assert events[0]["gene"] == "EWSR1::WT1"
    assert events[0]["provided_rna_junction_reads"] == "15"
    assert events[0]["rna_junction_reads"] == "0"
    assert events[0]["junction_match_status"] == "NO_EXACT_JUNCTION_MATCH"


def test_generator_passes_explicit_star_chimeric_for_fusion_read_backlink(tmp_path):
    hla = write(tmp_path / "hla.txt", "HLA-A*02:01\n")
    purity = write(tmp_path / "purity.tsv", "sample_id\tpurity\tploidy\nS1\t0.60\t2.0\n")
    cnv = write(tmp_path / "cnv.tsv", "chrom\tstart\tend\ttotal_cn\nchr1\t1\t1000\t2\n")
    lohhla = write(tmp_path / "lohhla.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    spechla_loh = write(tmp_path / "spechla_loh.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    easyfuse = write(tmp_path / "easyfuse/fusions.pass.csv", "gene5,gene3,Breakpoint1,Breakpoint2\nEWSR1,WT1,chr22:100,chr11:200\n")
    chimeric = write(tmp_path / "star/Chimeric.out.junction", "chr22\t100\t+\tchr11\t200\t-\t1\t0\t0\tread1\n")
    output = tmp_path / "production.toml"
    subprocess.run([
        sys.executable, str(ROOT / "scripts/generate_production_from_results_manifest.py"),
        "--project-root", str(ROOT), "--sample-id", "S1", "--outdir", str(tmp_path / "run"),
        "--output", str(output), "--hla-file", str(hla), "--purity", str(purity), "--cnv", str(cnv),
        "--lohhla", str(lohhla), "--spechla-loh", str(spechla_loh),
        "--easyfuse", str(easyfuse), "--star-chimeric", str(chimeric),
    ], check=True)
    command = tomllib.loads(output.read_text(encoding="utf-8"))["stages"]["fusion_candidates"]["command"]
    assert "--star-chimeric" in command
    assert str(chimeric) in command


def test_generator_adds_star_bam_splice_qc_before_prefilter(tmp_path):
    hla = write(tmp_path / "hla.txt", "HLA-A*02:01\n")
    purity = write(tmp_path / "purity.tsv", "sample_id\tpurity\tploidy\nS1\t0.60\t2.0\n")
    cnv = write(tmp_path / "cnv.tsv", "chrom\tstart\tend\ttotal_cn\nchr1\t1\t1000\t2\n")
    lohhla = write(tmp_path / "lohhla.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    spechla_loh = write(tmp_path / "spechla_loh.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    junctions = write(tmp_path / "junctions.tsv", "chrom\tstart\tend\tstrand\nchr1\t100\t200\t+\n")
    snaf = write(tmp_path / "snaf.tsv", "uid\nchr1:100-200:+\n")
    splicemutr = tmp_path / "splicemutr"
    splicemutr.mkdir()
    rna_bam = write(tmp_path / "rna.bam", "fixture\n")
    star_sj = write(tmp_path / "SJ.out.tab", "1\t100\t200\t1\t0\t0\t5\t0\t20\n")
    normal_sqlite = write(tmp_path / "normal.sqlite", "fixture\n")
    output = tmp_path / "production.toml"

    subprocess.run([
        sys.executable, str(ROOT / "scripts/generate_production_from_results_manifest.py"),
        "--project-root", str(ROOT), "--sample-id", "S1", "--outdir", str(tmp_path / "run"),
        "--output", str(output), "--hla-file", str(hla), "--purity", str(purity), "--cnv", str(cnv),
        "--lohhla", str(lohhla), "--spechla-loh", str(spechla_loh),
        "--junctions", str(junctions), "--snaf", str(snaf), "--splicemutr", str(splicemutr),
        "--splice-rna-bam", str(rna_bam), "--splice-star-sj", str(star_sj),
        "--normal-junction-sqlite", str(normal_sqlite),
    ], check=True)

    command = tomllib.loads(output.read_text(encoding="utf-8"))["stages"]["splice_candidates"]["command"]
    assert "build_splice_junction_qc_from_star_bam.py" in command
    assert str(rna_bam) in command and str(star_sj) in command
    assert str(normal_sqlite) in command
    assert command.index("build_splice_junction_qc_from_star_bam.py") < command.index("filter_splice_production_candidates.py")


def test_existing_upstream_results_produce_hla_purity_and_loh_consensus(tmp_path):
    optitype = write(tmp_path / "optitype/result.tsv", "A1\tA2\tB1\tB2\tC1\tC2\nA*02:01\tA*11:01\tB*15:01\tB*40:01\tC*03:04\tC*08:01\n").parent
    spechla = write(tmp_path / "spechla/hla.result.txt", "A\tA*02:01\tA*11:01\nB\tB*15:01\tB*40:01\nC\tC*03:04\tC*08:01\n").parent
    hla_out = tmp_path / "hla_consensus"
    subprocess.run([sys.executable, "-m", "neoag.agent_skills.hla_typing_compare", "--result-dir", str(optitype), "--result-dir", str(spechla), "--sample-id", "S1", "--outdir", str(hla_out)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    subprocess.run([sys.executable, str(ROOT / "scripts/hla_consensus_to_file.py"), "--consensus", str(hla_out / "hla_typing_consensus.tsv"), "--output", str(hla_out / "recommended_hla.txt")], cwd=ROOT, check=True)
    assert "HLA-A*02:01" in (hla_out / "recommended_hla.txt").read_text(encoding="utf-8")

    facets = write(tmp_path / "facets/facets_purity.txt", "purity\t0.60\n").parent
    sequenza = write(tmp_path / "sequenza/result.tsv", "purity\tploidy\n0.62\t2.1\n").parent
    purple = write(tmp_path / "purple/purity.tsv", "purity\tploidy\n0.61\t2.0\n").parent
    purity_out = tmp_path / "purity_consensus"
    subprocess.run([sys.executable, "-m", "neoag.agent_skills.purity_cnv_review", "--result-dir", str(facets), "--result-dir", str(sequenza), "--result-dir", str(purple), "--sample-id", "S1", "--outdir", str(purity_out)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    purity_rows = read_tsv(purity_out / "purity_cnv_consensus.tsv")
    assert purity_rows and int(purity_rows[0]["n_tools"]) >= 2

    lohhla = write(tmp_path / "lohhla.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    spechla_loh = write(tmp_path / "spechla_loh.tsv", "hla_allele\tloh_status\nHLA-A*02:01\tretained\n")
    loh_out = tmp_path / "loh_consensus"
    subprocess.run([sys.executable, str(ROOT / "scripts/build_hla_loh_consensus.py"), "--sample-id", "S1", "--lohhla", str(lohhla), "--spechla", str(spechla_loh), "--outdir", str(loh_out)], cwd=ROOT, check=True)
    assert read_tsv(loh_out / "hla_loh_consensus.tsv")[0]["consensus_status"] == "CONSENSUS_RETAINED"
