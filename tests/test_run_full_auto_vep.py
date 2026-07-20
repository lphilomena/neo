from pathlib import Path

from neoag_v03.tools import upstream
from neoag_v03.tools import runner
from neoag_v03.tools.registry import RunContext


def test_vcf_has_csq_annotations_plain_file_with_gz_suffix(tmp_path):
    vcf = tmp_path / "annotated.vcf.gz"
    vcf.write_text(
        '##fileformat=VCFv4.2\n'
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Allele|Consequence">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        encoding="utf-8",
    )

    assert upstream.vcf_has_csq_annotations(vcf)


def test_auto_annotates_unannotated_variant_peptide_vcf(tmp_path, monkeypatch):
    input_vcf = tmp_path / "input.vcf"
    input_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t10\t.\tA\tT\t.\tPASS\t.\n",
        encoding="utf-8",
    )
    ref_fasta = tmp_path / "ref.fa"
    ref_fasta.write_text(">chr1\nAAAA\n", encoding="utf-8")
    tools_dir = tmp_path / "tools"
    tools_dir.mkdir()

    def fake_run_vep(**kwargs):
        out = Path(kwargs["output_vcf"])
        out.write_text(
            '##fileformat=VCFv4.2\n'
            '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Allele|Consequence">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t10\t.\tA\tT\t.\tPASS\tCSQ=T|missense_variant\n",
            encoding="utf-8",
        )
        return {"annotated_vcf": str(out), "command": "vep mock"}

    monkeypatch.setattr("neoag_v03.vep.annotate.run_vep_pvacseq_annotate", fake_run_vep)

    annotated, outputs = upstream._auto_annotate_variants_vcf(
        {
            "inputs": {
                "variant_peptide_extraction": True,
                "reference_fasta": str(ref_fasta),
            }
        },
        variants_vcf=input_vcf,
        tools_dir=tools_dir,
        sample_id="S1",
    )

    assert annotated == tools_dir / "S1.vep.annotated.vcf"
    assert upstream.vcf_has_csq_annotations(annotated)
    assert outputs["vep_annotated_vcf"] == str(annotated)


def test_run_upstream_auto_annotates_before_variant_peptides(tmp_path, monkeypatch):
    input_vcf = tmp_path / "input.vcf"
    input_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t10\t.\tA\tT\t.\tPASS\t.\n",
        encoding="utf-8",
    )
    ref_fasta = tmp_path / "ref.fa"
    ref_fasta.write_text(">chr1\nAAAA\n", encoding="utf-8")
    cfg = tmp_path / "run.toml"
    cfg.write_text(
        f"""
[sample]
id = "S1"
profile = "default"

[tools]
enabled = []

[inputs]
entry_mode = "snv_indel"
variant_peptide_extraction = true
variants_vcf = "{input_vcf}"
reference_fasta = "{ref_fasta}"
hla_alleles = ["HLA-A*02:01"]
extract_appm_from_vcf = false
""",
        encoding="utf-8",
    )

    def fake_run_vep(**kwargs):
        out = Path(kwargs["output_vcf"])
        out.write_text(
            '##fileformat=VCFv4.2\n'
            '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Allele|Consequence">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t10\t.\tA\tT\t.\tPASS\tCSQ=T|missense_variant\n",
            encoding="utf-8",
        )
        return {"annotated_vcf": str(out), "command": "vep mock"}

    def fake_variant_upstream(_cfg, *, variants_vcf, parsed_dir, **_kwargs):
        assert variants_vcf.name == "S1.vep.annotated.vcf"
        assert upstream.vcf_has_csq_annotations(variants_vcf)
        tools_dir = _kwargs["tools_dir"]
        variant_peptides = tools_dir / "variant_peptides.tsv"
        raw_peptides = parsed_dir / "raw_peptides.tsv"
        raw_events = parsed_dir / "raw_events.tsv"
        variant_peptides.write_text("peptide_id\tmutant_peptide\thla_allele\n", encoding="utf-8")
        raw_peptides.write_text("peptide_id\tpeptide\thla_allele\n", encoding="utf-8")
        raw_events.write_text("event_id\tgene\n", encoding="utf-8")
        return {
            "variant_peptides": str(variant_peptides),
            "raw_peptides": str(raw_peptides),
            "raw_events": str(raw_events),
        }

    monkeypatch.setattr("neoag_v03.vep.annotate.run_vep_pvacseq_annotate", fake_run_vep)
    monkeypatch.setattr(upstream, "run_variant_peptide_upstream", fake_variant_upstream)

    outputs = upstream.run_upstream(cfg, tmp_path / "out")

    assert outputs["vep_annotated_vcf"].endswith("S1.vep.annotated.vcf")
    assert outputs["raw_peptides"].endswith("raw_peptides.tsv")


def test_vep_appm_prefers_config_reference_fasta(tmp_path, monkeypatch):
    variants_vcf = tmp_path / "variants.vcf"
    variants_vcf.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
    config_fasta = tmp_path / "config.fa"
    config_fasta.write_text(">chr1\nAAAA\n", encoding="utf-8")
    env_fasta = tmp_path / "environment.fa"
    env_fasta.write_text(">chr1\nCCCC\n", encoding="utf-8")
    captured = {}

    def fake_run_cmd(cmd, _workdir, **_kwargs):
        captured["cmd"] = cmd

    def fake_vep_to_appm(_raw, out_tsv):
        out_tsv.write_text("event_id\tgene\n", encoding="utf-8")

    monkeypatch.setenv("NEOAG_REFERENCE_FASTA", str(env_fasta))
    monkeypatch.setenv("NEOAG_RUNNER_MODE", "conda")
    monkeypatch.setattr(runner, "_run_cmd", fake_run_cmd)
    monkeypatch.setattr(runner, "vep_to_appm_tsv", fake_vep_to_appm)

    ctx = RunContext(
        sample_id="S1",
        outdir=tmp_path / "out",
        variants_vcf=variants_vcf,
        reference_fasta=config_fasta,
    )
    out_tsv = tmp_path / "vep_appm.tsv"
    runner.run_vep_appm(ctx, out_tsv)

    cmd = captured["cmd"]
    fasta_index = cmd.index("--fasta")
    assert cmd[fasta_index + 1] == str(config_fasta)
    assert str(env_fasta) not in cmd
