import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
WRAPPER = ROOT / "scripts" / "run_production_case.sh"


def make_star_index(path: Path) -> Path:
    path.mkdir(parents=True)
    for name in ("Genome", "SA", "SAindex", "genomeParameters.txt"):
        (path / name).write_text("fixture\n")
    return path


def test_wrapper_uses_v2_baseline_and_v3_consensus_and_passes_production_inputs(tmp_path):
    project = tmp_path / "project"
    profile = project / "profiles/sarcoma_rna_supported_v2_provisional.toml"
    profile.parent.mkdir(parents=True)
    profile.write_text("[metadata]\nname='test'\n")
    consensus_rules = project / "configs/ranking/sarcoma_evidence_consensus_v3_source_chain.toml"
    consensus_rules.parent.mkdir(parents=True)
    consensus_rules.write_text("[metadata]\nname='consensus-test'\n")
    (project / "scripts").mkdir()

    case_root = tmp_path / "case"
    case_root.mkdir()
    outdir = tmp_path / "out"
    stabpan = tmp_path / "stabpan"
    (stabpan / "Linux_x86_64/bin").mkdir(parents=True)
    (stabpan / "data").mkdir()
    stabpan_bin = stabpan / "Linux_x86_64/bin/netMHCstabpan"
    stabpan_bin.write_text("#!/bin/sh\nexit 0\n")
    stabpan_bin.chmod(0o755)

    paths = {}
    for name in (
        "somatic.vcf.gz", "reference.fa", "gencode.gtf", "rna_R1.fq.gz",
        "rna_R2.fq.gz", "sequenza.tsv", "purple.tsv",
    ):
        path = tmp_path / name
        path.write_text("fixture\n")
        paths[name] = path
    star_index = make_star_index(tmp_path / "star_index")
    star_chimeric = tmp_path / "Chimeric.out.junction"
    star_chimeric.write_text("chr22\t100\t+\tchr11\t200\t-\t1\t0\t0\tread1\n")
    star = tmp_path / "STAR"
    samtools = tmp_path / "samtools"
    for executable in (star, samtools):
        executable.write_text("#!/bin/sh\nexit 0\n")
        executable.chmod(0o755)

    invocation_log = tmp_path / "python.invocations"
    fake_python = tmp_path / "python"
    fake_python.write_text(
        "#!/bin/sh\n"
        f"printf '%s\\n' \"$*\" >> {invocation_log}\n"
        "exit 0\n"
    )
    fake_python.chmod(0o755)

    subprocess.run([
        "bash", str(WRAPPER),
        "--sample-id", "S1",
        "--case-root", str(case_root),
        "--outdir", str(outdir),
        "--somatic-vcf", str(paths["somatic.vcf.gz"]),
        "--project-root", str(project),
        "--python", str(fake_python),
        "--netmhcstabpan-home", str(stabpan),
        "--reference-fasta", str(paths["reference.fa"]),
        "--gencode-gtf", str(paths["gencode.gtf"]),
        "--sequenza", str(paths["sequenza.tsv"]),
        "--purple", str(paths["purple.tsv"]),
        "--rna-fastq1", str(paths["rna_R1.fq.gz"]),
        "--rna-fastq2", str(paths["rna_R2.fq.gz"]),
        "--star-index", str(star_index),
        "--star-chimeric", str(star_chimeric),
        "--star-executable", str(star),
        "--samtools-executable", str(samtools),
        "--rna-threads", "7",
    ], check=True)

    invocation_text = invocation_log.read_text()
    invocations = invocation_text.splitlines()
    assert "event_track" in invocation_text and "splice_junction" in invocation_text
    assert "WT_STRONG_BINDING_REVIEW" in invocation_text
    assert "PUTATIVE_TCR_FACING" in invocation_text
    assert "STRUCTURAL_ROLE_UNCERTAIN" in invocation_text
    generator_call = next(line for line in invocations if "generate_production_from_results_manifest.py" in line)
    assert "--profile" in generator_call
    assert "sarcoma_rna_supported_v2_provisional.toml" in generator_call
    assert "--evidence-consensus-rules" in generator_call
    assert "sarcoma_evidence_consensus_v3_source_chain.toml" in generator_call
    for flag in (
        "--reference-fasta", "--gencode-gtf", "--sequenza", "--purple",
        "--rna-fastq1", "--rna-fastq2", "--star-index", "--star-executable",
        "--star-chimeric", "--samtools-executable", "--rna-threads",
    ):
        assert flag in generator_call

    wrapper_text = WRAPPER.read_text()
    assert "verify_mtwt_output_fields" in wrapper_text
    assert "ranked_peptides.evidence_consensus.tsv" in wrapper_text
    assert "mt_wt_interpretation_caution" in wrapper_text
    assert "verify_splice_prefilter_outputs" in wrapper_text
    assert "splice_prefilter_funnel.tsv" in wrapper_text
    assert "UNIQUE_JUNCTION_READS" in wrapper_text
    assert "TOTAL_JUNCTION_COVERAGE" in wrapper_text
    assert "PSI" in wrapper_text


def test_wrapper_rejects_multiple_rna_input_modes(tmp_path):
    result = subprocess.run([
        "bash", str(WRAPPER), "--sample-id", "S1", "--case-root", str(tmp_path),
        "--outdir", str(tmp_path / "out"), "--somatic-vcf", str(tmp_path / "x.vcf"),
        "--rna-fastq1", "r1.fq.gz", "--rna-fastq2", "r2.fq.gz",
        "--rna-bam", "rna.bam", "--star-index", "index", "--gencode-gtf", "g.gtf",
    ], text=True, capture_output=True)
    assert result.returncode == 2
    assert "only one RNA allele-evidence input mode" in result.stderr
