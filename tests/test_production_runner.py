from pathlib import Path
from datetime import datetime

from neoag.production_runner import run_production
from neoag.tools.registry import ROOT
from neoag.utils import read_tsv


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

[evidence]
evidence_consensus_rules = "{ROOT / 'configs/ranking/sarcoma_evidence_consensus_v3_source_chain.toml'}"

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

    assert result.status == "FAILED"
    assert result.source_status == "COMPLETE"
    assert set(result.detected_sources) == {"pVACseq", "pVACfuse", "pVACsplice"}
    assert (tmp_path / "run/final/scoring/ranked_peptides.tsv").is_file()
    assert (tmp_path / "run/final/reports/evidence_report.patient.html").is_file()
    assert (tmp_path / "run/final/reports/evidence_report.technical.html").is_file()
    coverage = read_tsv(tmp_path / "run/peptide_source_coverage.tsv")[0]
    assert coverage["status"] == "COMPLETE"
    assert any(stage.name == "mtwt_assessment_gate" and stage.status == "FAILED" for stage in result.stages)


def test_production_runner_missing_optional_sources_is_low_confidence(tmp_path):
    result = run_production(
        _manifest(tmp_path, include_fusion=False, include_splice=False),
        outdir=tmp_path / "run",
        project_root=ROOT,
        execute=True,
    )

    assert result.status == "FAILED"
    assert result.source_status == "LOW_CONFIDENCE"
    assert result.missing_sources == ["pVACfuse", "pVACsplice"]
    assert (tmp_path / "run/final/scoring/ranked_peptides.tsv").is_file()
    assert any(stage.name == "mtwt_assessment_gate" and stage.status == "FAILED" for stage in result.stages)


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
    generated = (tmp_path / "run/run.production.generated.toml").read_text(encoding="utf-8")
    assert "evidence_consensus_rules" in generated
    assert "sarcoma_evidence_consensus_v3_source_chain.toml" in generated


def test_production_runner_blocks_incomplete_cross_tool_group(tmp_path):
    hla = tmp_path / "hla.txt"; hla.write_text("HLA-A*02:01\n")
    result_file = tmp_path / "result.tsv"; result_file.write_text("status\nPASS\n")
    manifest = tmp_path / "gate.toml"
    manifest.write_text(f'''[run]
sample_id = "GATE"
hla_file = "{hla}"

[run.required_tool_groups.purity_cnv]
tools = ["facets", "sequenza", "purple"]
min_successful = 2
require_all_declared = true

[stages.purity_facets]
required = false
[stages.purity_facets.outputs]
facets_result = "{result_file}"

[stages.purity_sequenza]
required = false
[stages.purity_sequenza.outputs]
sequenza_result = "{result_file}"
''')
    result = run_production(manifest, outdir=tmp_path / "run", project_root=ROOT, execute=True)
    assert result.status == "BLOCKED"
    gate = read_tsv(tmp_path / "run/production_release_gate.tsv")[0]
    assert gate["status"] == "FAIL"
    assert "purple" in gate["reason"]


def test_production_runner_runs_independent_stages_within_resource_budget(tmp_path):
    manifest = tmp_path / "parallel.toml"
    manifest.write_text('''[run]
sample_id = "PARALLEL"
total_cpus = 2
total_memory_gb = 4
max_parallel_stages = 2

[stages.a]
required = true
cpus = 1
memory_gb = 2
command = "mkdir -p {outdir}; sleep 0.4; touch {outdir}/a.done"
[stages.a.outputs]
done = "{outdir}/a.done"

[stages.b]
required = true
cpus = 1
memory_gb = 2
command = "mkdir -p {outdir}; sleep 0.4; touch {outdir}/b.done"
[stages.b.outputs]
done = "{outdir}/b.done"
''', encoding="utf-8")

    result = run_production(
        manifest, outdir=tmp_path / "run", project_root=ROOT,
        execute=True, skip_ranking=True,
    )

    stages = {stage.name: stage for stage in result.stages}
    assert stages["a"].status == "PASS"
    assert stages["b"].status == "PASS"
    starts = [datetime.fromisoformat(stages[name].started_at) for name in ("a", "b")]
    assert abs((starts[0] - starts[1]).total_seconds()) < 0.3
    schedule = read_tsv(tmp_path / "run/production_resource_schedule.tsv")
    assert {row["cpus"] for row in schedule} == {"1"}


def test_production_runner_queues_stages_when_memory_budget_is_full(tmp_path):
    manifest = tmp_path / "memory_queue.toml"
    manifest.write_text('''[run]
sample_id = "MEMORY_QUEUE"
total_cpus = 2
total_memory_gb = 4
max_parallel_stages = 2

[stages.a]
required = true
cpus = 1
memory_gb = 3
command = "mkdir -p {outdir}; sleep 0.35; touch {outdir}/a.done"
[stages.a.outputs]
done = "{outdir}/a.done"

[stages.b]
required = true
cpus = 1
memory_gb = 3
command = "mkdir -p {outdir}; sleep 0.35; touch {outdir}/b.done"
[stages.b.outputs]
done = "{outdir}/b.done"
''', encoding="utf-8")

    result = run_production(
        manifest, outdir=tmp_path / "run", project_root=ROOT,
        execute=True, skip_ranking=True,
    )

    stages = {stage.name: stage for stage in result.stages}
    first_finished = datetime.fromisoformat(stages["a"].finished_at)
    second_started = datetime.fromisoformat(stages["b"].started_at)
    assert second_started >= first_finished


def test_production_runner_blocks_stage_larger_than_global_budget(tmp_path):
    manifest = tmp_path / "oversized.toml"
    manifest.write_text('''[run]
sample_id = "OVERSIZED"
total_cpus = 2
total_memory_gb = 4
max_parallel_stages = 2

[stages.heavy]
required = true
cpus = 3
memory_gb = 2
command = "touch {outdir}/should-not-run"
[stages.heavy.outputs]
done = "{outdir}/should-not-run"
''', encoding="utf-8")

    result = run_production(manifest, outdir=tmp_path / "run", project_root=ROOT, execute=True)

    assert result.status == "BLOCKED"
    assert result.stages[0].status == "BLOCKED"
    assert "exceeds global budget" in result.stages[0].message
    assert not (tmp_path / "run/should-not-run").exists()


def test_production_runner_does_not_reuse_header_only_required_table(tmp_path):
    events = tmp_path / "raw_events.tsv"
    peptides = tmp_path / "raw_peptides.tsv"
    events.write_text("event_id\nE1\n", encoding="utf-8")
    peptides.write_text("peptide_id\n", encoding="utf-8")
    manifest = tmp_path / "header_only.toml"
    manifest.write_text(
        f'''[run]
sample_id = "HEADER_ONLY"

[stages.fusion]
required = true
source = "EasyFuse"
data_row_outputs = ["raw_peptides"]
command = "true"
[stages.fusion.outputs]
raw_events = "{events}"
raw_peptides = "{peptides}"
''',
        encoding="utf-8",
    )

    result = run_production(manifest, outdir=tmp_path / "run", project_root=ROOT)

    assert result.stages[0].status == "PLANNED"
