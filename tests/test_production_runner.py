from pathlib import Path

from neoag_v03.production_runner import run_production
from neoag_v03.tools.registry import ROOT
from neoag_v03.utils import read_tsv


def _manifest(tmp_path: Path, *, include_fusion: bool = True, include_splice: bool = True) -> Path:
    hla = tmp_path / "hla.txt"
    hla.write_text("HLA-A*02:01\nHLA-B*07:02\n", encoding="utf-8")
    stages = [
        f'''[stages.snv]
required = true
source = "pVACseq"
[stages.snv.outputs]
pvac_file = "{ROOT / 'data/fixtures/pvacseq_aggregated.tsv'}"
'''
    ]
    if include_fusion:
        stages.append(
            f'''[stages.fusion]
required = false
source = "pVACfuse"
depends_on = ["snv"]
[stages.fusion.outputs]
pvac_file = "{ROOT / 'data/fixtures/pvacfuse_aggregated.tsv'}"
'''
        )
    if include_splice:
        stages.append(
            f'''[stages.splice]
required = false
source = "pVACsplice"
depends_on = ["snv"]
[stages.splice.outputs]
pvac_file = "{ROOT / 'data/fixtures/pvacsplice_aggregated.tsv'}"
'''
        )
    path = tmp_path / "production.toml"
    path.write_text(
        f'''[run]
sample_id = "PROD1"
profile = "default"
hla_file = "{hla}"
tools_stub = true
immunogenicity_stub = true
expected_peptide_sources = ["pVACseq", "pVACfuse", "pVACsplice"]
presentation_predictors = ["netmhcpan", "mhcflurry"]
required_presentation_predictors = ["netmhcpan", "mhcflurry"]

{''.join(stages)}
''',
        encoding="utf-8",
    )
    return path


def test_production_runner_merges_all_sources_and_ranks(tmp_path):
    result = run_production(
        _manifest(tmp_path),
        outdir=tmp_path / "run",
        project_root=ROOT,
        execute=True,
    )

    assert result.status == "PASS"
    assert result.source_status == "COMPLETE"
    assert set(result.detected_sources) == {"pVACseq", "pVACfuse", "pVACsplice"}
    assert (tmp_path / "run/final/scoring/ranked_peptides.v03.tsv").is_file()
    coverage = read_tsv(tmp_path / "run/peptide_source_coverage.tsv")[0]
    assert coverage["status"] == "COMPLETE"


def test_production_runner_missing_optional_sources_is_low_confidence(tmp_path):
    result = run_production(
        _manifest(tmp_path, include_fusion=False, include_splice=False),
        outdir=tmp_path / "run",
        project_root=ROOT,
        execute=True,
    )

    assert result.status == "LOW_CONFIDENCE"
    assert result.source_status == "LOW_CONFIDENCE"
    assert result.missing_sources == ["pVACfuse", "pVACsplice"]
    assert (tmp_path / "run/final/scoring/ranked_peptides.v03.tsv").is_file()


def test_production_runner_dry_run_plans_commands(tmp_path):
    manifest = tmp_path / "dry.toml"
    manifest.write_text(
        '''[run]
sample_id = "DRY"

[stages.hla]
required = true
command = "touch {outdir}/hla.txt"
[stages.hla.outputs]
hla_file = "{outdir}/hla.txt"

[stages.fusion]
required = false
source = "pVACfuse"
command = "false"
[stages.fusion.outputs]
pvac_file = "{outdir}/fusion.tsv"
''',
        encoding="utf-8",
    )

    result = run_production(manifest, outdir=tmp_path / "run", project_root=ROOT)

    assert result.status == "DRY_RUN"
    assert result.dry_run
    assert [stage.status for stage in result.stages] == ["PLANNED", "PLANNED"]
    assert not (tmp_path / "run/hla.txt").exists()


def test_production_runner_reuses_declared_outputs(tmp_path):
    manifest = _manifest(tmp_path, include_fusion=False, include_splice=False)
    result = run_production(
        manifest,
        outdir=tmp_path / "run",
        project_root=ROOT,
        execute=True,
        skip_ranking=True,
    )

    assert result.stages[0].status == "REUSED"
    assert result.stages[-1].status == "SKIPPED"
    assert (tmp_path / "run/merged/raw_peptides.tsv").is_file()
