from __future__ import annotations

import csv
import io
import json
import os
import subprocess
import sys
import tarfile
from pathlib import Path

from neoag.open_neo.contracts import MacroResult, MacroStep, RunInput, validate_json_schema
from neoag.open_neo.capability_planner import build_automatic_production_plan
from neoag.open_neo.errors import FailureCode, exit_code_for_result
from neoag.controlled_execution.doctor import CheckRow
from neoag.open_neo.install_check import (
    DEFAULT_ASSET_SOURCE_HOST,
    DEFAULT_ASSET_SOURCE_ROOT,
    _apply_default_asset_source,
    _apply_production_runtime_fixes,
    _asset_manifest_uses_placeholder_source,
    _assess_tier,
    _deployment_command,
    _collect_claude_code_status,
    _required_asset_sources_missing,
    _rewrite_asset_manifest,
    _run_deployment,
    _contig_style,
    _vep_cache_root,
    _stage_release,
    _write_local_manifest_templates,
    run_install_check,
)
from neoag.open_neo.auto_config import configure_machine
from neoag.open_neo.cli import build_parser
from neoag.open_neo.review import _event_kind, _merge_review_context, build_review_rows, run_review, select_first_batch
from neoag.open_neo.review_integrity import audit_review_inputs
from neoag.open_neo.routing import inspect_manifest
from neoag.open_neo.run import run_open_neo
from neoag.open_neo.public_assets import _download, _extract_split_archive, _normalize_hf_endpoint, sync_public_assets
from neoag.controlled_execution.pipeline_runner import _manifest_file_hashes
from neoag.open_neo.state import RunLayout, load_run_state, resume_step_decision
from neoag.open_neo.rna_preprocessing import prepare_rna_evidence
from neoag.open_neo.rna_fusion_splice_profile import (
    generate_rna_fusion_splice_manifest,
    is_rna_fastq_profile_candidate,
)
from neoag.open_neo.tool_consensus import build_tool_consensus
from neoag.production_runner import load_production_manifest, run_production
from neoag.sample_identity.bam_matcher import parse_bam_matcher_short
from neoag.skill_taxonomy.registry import SKILLS_BY_NAME
from neoag.skill_taxonomy.runner import run_skill
from neoag.skill_taxonomy.review_skills import _validation_route


def test_run_input_contract_and_public_schema_accept_rna_fastq():
    schema = json.loads(
        (Path.cwd() / ".agents/skills/open-neo-run/references/INPUT_SCHEMA.json").read_text(encoding="utf-8")
    )
    request = RunInput.from_mapping({
        "outdir": "work/run",
        "tumor_rna_fastq": ["tumor_R1.fastq.gz", "tumor_R2.fastq.gz"],
        "star_index": "/ref/star",
        "ctat_genome_lib": "/ref/ctat",
        "rna_threads": 8,
    }).to_mapping()
    assert validate_json_schema(request, schema) == []
    assert validate_json_schema({"outdir": "work/run"}, schema)
    invalid = dict(request, rna_threads=0, tumor_rna_fastq=[1])
    assert "BELOW_MINIMUM:rna_threads:1" in validate_json_schema(invalid, schema)
    assert "INVALID_ITEM_TYPE:tumor_rna_fastq" in validate_json_schema(invalid, schema)


def test_public_asset_cli_defaults_and_opt_out():
    install = build_parser().parse_args([
        "install-check", "--outdir", "work/install", "--no-sync-public-assets",
        "--hf-endpoint", "https://hf-mirror.com",
    ])
    run = build_parser().parse_args([
        "run", "--outdir", "work/run", "--input-dir", "inputs",
        "--public-asset-root", "/srv/open-neo/refs",
    ])
    assert install.sync_public_assets is False
    assert install.hf_endpoint == "https://hf-mirror.com"
    assert run.sync_public_assets is True
    assert run.public_asset_repo == "open-neo/open-neo-public-assets"
    assert run.public_asset_root == "/srv/open-neo/refs"


def test_install_readiness_normalizes_vep_cache_and_contig_style(tmp_path: Path):
    cache = tmp_path / "vep"
    version = cache / "homo_sapiens/112_GRCh38"
    version.mkdir(parents=True)
    assert _vep_cache_root(str(cache)) == cache
    assert _vep_cache_root(str(version)) == cache
    assert _vep_cache_root(str(tmp_path / "empty")) is None
    assert _contig_style("chr1") == "CHR"
    assert _contig_style("1") == "NO_CHR"
    assert _contig_style("") == "UNKNOWN"


def test_public_asset_plan_is_offline_and_detects_marker(tmp_path: Path):
    root = tmp_path / "refs"
    planned = sync_public_assets(root, cache_dir=tmp_path / "cache", execute=False)
    assert planned["status"] == "PLANNED"
    for relative in ("data/ref", "data/normal", "data/easyfuse"):
        (root / relative).mkdir(parents=True, exist_ok=True)
    (root / ".open_neo_public_assets.json").write_text("{}\n", encoding="utf-8")
    reused = sync_public_assets(root, cache_dir=tmp_path / "cache", execute=False)
    assert reused["status"] == "REUSED"


def test_public_asset_split_archive_extracts_safely(tmp_path: Path):
    payload = io.BytesIO()
    with tarfile.open(fileobj=payload, mode="w:gz") as archive:
        info = tarfile.TarInfo("ref/example.txt")
        content = b"open-neo\n"
        info.size = len(content)
        archive.addfile(info, io.BytesIO(content))
    raw = payload.getvalue()
    split = len(raw) // 2
    parts = [tmp_path / "part-000", tmp_path / "part-001"]
    parts[0].write_bytes(raw[:split])
    parts[1].write_bytes(raw[split:])
    destination = tmp_path / "data"
    _extract_split_archive(parts, destination)
    assert (destination / "ref/example.txt").read_text(encoding="utf-8") == "open-neo\n"


def test_public_asset_download_uses_official_hf_cli_and_local_cache(
    tmp_path: Path, monkeypatch,
):
    monkeypatch.delenv("HF_ENDPOINT", raising=False)
    fake_hf = tmp_path / "bin/hf"
    fake_hf.parent.mkdir()
    fake_hf.write_text("#!/bin/sh\n", encoding="utf-8")
    fake_hf.chmod(0o755)
    calls: list[list[str]] = []

    environments: list[dict[str, str]] = []

    def fake_run(command, **kwargs):
        calls.append(command)
        environments.append(kwargs["env"])
        local_dir = Path(command[command.index("--local-dir") + 1])
        remote_path = command[3]
        output = local_dir / remote_path
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_bytes(b"asset-data")
        return subprocess.CompletedProcess(command, 0, stdout=f"{output}\n", stderr="")

    monkeypatch.setattr("neoag.open_neo.public_assets.shutil.which", lambda name: str(fake_hf) if name == "hf" else None)
    monkeypatch.setattr("neoag.open_neo.public_assets.subprocess.run", fake_run)
    destination = tmp_path / "cache/archive.part-000"
    action = _download(
        "open-neo/open-neo-public-assets", "main", "archive.part-000", destination,
        expected_size=len(b"asset-data"),
    )
    assert action == "DOWNLOADED"
    assert destination.read_bytes() == b"asset-data"
    assert calls[0][1:4] == ["download", "open-neo/open-neo-public-assets", "archive.part-000"]
    assert "--repo-type" in calls[0]
    assert "--local-dir" in calls[0]
    assert environments[0]["HF_ENDPOINT"] == "https://huggingface.co"


def test_public_asset_endpoint_prefers_argument_then_environment(monkeypatch):
    monkeypatch.setenv("HF_ENDPOINT", "https://hf-mirror.com/")
    assert _normalize_hf_endpoint() == "https://hf-mirror.com"
    assert _normalize_hf_endpoint("https://huggingface.co/") == "https://huggingface.co"


def test_run_state_requires_matching_output_signature_for_reuse(tmp_path: Path):
    source = tmp_path / "source.tsv"
    source.write_text("input\n", encoding="utf-8")
    output = tmp_path / "result.tsv"
    output.write_text("a\n", encoding="utf-8")
    result_path = tmp_path / "skill_result.json"
    result = MacroResult(skill="open-neo-run", run_id="RUN1", case_id="CASE1", mode="execute")
    result.steps.append(MacroStep("ranking", "ranking", "PASS", inputs={"source": str(source)}, outputs={"table": str(output)}))
    result.finish("PASS").write(result_path)
    state = load_run_state(tmp_path / "run_state.json")
    decision = resume_step_decision(state, "ranking")
    assert decision["decision"] == "REUSE"
    output.write_text("changed\n", encoding="utf-8")
    decision = resume_step_decision(state, "ranking")
    assert decision["decision"] == "RUN"
    assert decision["reason"].startswith("OUTPUT_HASH_CHANGED:")


def test_failure_codes_have_stable_cli_exit_mapping():
    assert exit_code_for_result({"status": "PASS"}) == 0
    assert exit_code_for_result({"status": "BLOCKED", "blocking_issues": [FailureCode.APPROVAL_REQUIRED.value]}) == 3
    assert exit_code_for_result({"status": "NEEDS_RANKING", "blocking_issues": [FailureCode.NEEDS_RANKING.value]}) == 4


def test_review_event_kind_prioritizes_explicit_dna_event_type():
    assert _event_kind("SNV", "splice_junction") == "MISSENSE"
    assert _event_kind("InDel", "splice_acceptor_variant") == "FRAMESHIFT"
    assert _event_kind("", "splice_junction") == "SPLICE"


def test_validation_route_prioritizes_explicit_event_type_over_splice_consequence():
    snv = {"event_type": "SNV", "peptide_consequence": "splice_junction", "wildtype_peptide": "AAAAAAAAA"}
    indel = {"event_type": "InDel", "peptide_consequence": "splice_junction", "wildtype_peptide": "AAAAAAAAA"}
    splice = {"event_type": "Splice", "peptide_consequence": "splice_junction"}
    assert _validation_route(snv) == "short_peptide_plus_wt_control"
    assert _validation_route(indel) == "short_peptide_plus_wt_control"
    assert _validation_route(splice) == "targeted_rna_then_junction_long_peptide_or_minigene"


def test_review_cli_accepts_report_selection():
    args = build_parser().parse_args([
        "review", "--result-dir", "results/case", "--outdir", "reviews/case",
        "--reports", "patient,technical",
    ])
    assert args.reports == ["patient", "technical"]
    schema = json.loads(
        (Path.cwd() / ".agents/skills/open-neo-review/references/INPUT_SCHEMA.json").read_text(encoding="utf-8")
    )
    assert validate_json_schema(vars(args), schema) == []
    invalid = dict(vars(args), reports=["clinical"])
    assert "INVALID_ITEM_ENUM:reports" in validate_json_schema(invalid, schema)


def test_review_clinical_context_overrides_stale_run_metadata():
    merged = _merge_review_context(
        {"disease": "sarcoma_profile", "sample_id": "CASE001"},
        {"disease": "DSRCT", "analysis_goal": "neoantigen review", "notes": ""},
    )
    assert merged["disease"] == "DSRCT"
    assert merged["analysis_goal"] == "neoantigen review"
    assert merged["sample_id"] == "CASE001"
    assert "notes" not in merged


def test_public_macro_registry_and_internal_composition():
    for name in ("open-neo-install-check", "open-neo-run", "open-neo-review"):
        spec = SKILLS_BY_NAME[name]
        assert spec.category == "M"
        assert spec.visibility == "public"
        assert spec.entrypoint.startswith("open-neo ")
        assert spec.composes
    assert SKILLS_BY_NAME["neoag-vcf"].category == "A"
    assert SKILLS_BY_NAME["neoag-ranking"].category == "B"
    assert SKILLS_BY_NAME["neoag-experiment-design"].category == "C"
    assert SKILLS_BY_NAME["neoag-doctor"].category == "D"


def test_input_detection_multi_route(tmp_path: Path):
    project = Path.cwd()
    manifest = tmp_path / "sample.yaml"
    manifest.write_text(
        f"""schema_version: open-neo-sample-v1
case_id: CASE_ROUTE
sample_id: SAMPLE_ROUTE
genome_build: GRCh38
inputs:
  fusion_tsv: {project / 'data/fixtures/easyfuse_fusions.pass.tsv'}
  splice_junction_tsv: {project / 'data/fixtures/regtools_splice_junctions.tsv'}
  hla_alleles:
    - HLA-A*02:01
    - HLA-B*07:02
execution:
  profile: default
""",
        encoding="utf-8",
    )
    result = inspect_manifest(manifest)
    assert result.status == "PASS"
    assert {route.route for route in result.routes} == {"fusion", "splice"}
    assert not result.missing


def test_install_check_review_tier_source_checkout(tmp_path: Path):
    result = run_install_check({
        "project_root": str(Path.cwd()),
        "deployment_tier": "review",
        "mode": "plan",
        "release_audit": False,
        "outdir": str(tmp_path / "install"),
    })
    assert result["status"] in {"READY", "PARTIAL"}
    assert Path(result["outputs"]["deployment_report"]).is_file()


def test_install_check_release_requires_checksum(tmp_path: Path):
    tarball = tmp_path / "release.tar"
    tarball.write_bytes(b"release")
    result = run_install_check({
        "project_root": str(Path.cwd()),
        "release_tarball": str(tarball),
        "mode": "plan",
        "outdir": str(tmp_path / "install"),
    })
    assert result["status"] == "BLOCKED"
    assert "CHECKSUM_REQUIRED" in result["blocking_issues"]


def test_install_check_safely_stages_release_tarball(tmp_path: Path):
    source = tmp_path / "source/release"
    source.mkdir(parents=True)
    (source / "pyproject.toml").write_text("[project]\nname='release'\n", encoding="utf-8")
    (source / "README.md").write_text("release\n", encoding="utf-8")
    archive = tmp_path / "release.tar.gz"
    with tarfile.open(archive, "w:gz") as handle:
        handle.add(source, arcname="release")
    import hashlib
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    layout = RunLayout.create(tmp_path / "install")
    project_root, outputs, reused = _stage_release(archive, digest, layout)
    assert (project_root / "pyproject.toml").is_file()
    assert Path(outputs["release_staging"]).is_file()
    assert reused is False
    project_root_2, _, reused_2 = _stage_release(archive, digest, layout)
    assert project_root_2 == project_root
    assert reused_2 is True


def test_install_check_rejects_unsafe_release_member(tmp_path: Path):
    archive = tmp_path / "unsafe.tar"
    payload = b"unsafe\n"
    with tarfile.open(archive, "w") as handle:
        member = tarfile.TarInfo("../escape.txt")
        member.size = len(payload)
        handle.addfile(member, io.BytesIO(payload))
    import hashlib
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    layout = RunLayout.create(tmp_path / "install")
    try:
        _stage_release(archive, digest, layout)
        assert False, "unsafe archive should be rejected"
    except ValueError as exc:
        assert "unsafe archive path" in str(exc)
    assert not (tmp_path / "escape.txt").exists()


def test_prediction_tier_requires_complete_reference_manifest(tmp_path: Path):
    tools = [
        "python", "neoag", "neoag-skill", "vep", "netmhcpan", "mhcflurry",
        "netmhcstabpan", "netchop", "pvacseq", "prime",
    ]
    rows = [CheckRow("tool", name, "OK") for name in tools]
    status, requirements = _assess_tier("prediction", rows, None)
    assert status == "PARTIAL"
    assert any(row["kind"] == "reference" and row["status"] == "MISSING" for row in requirements)

    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\nA\n", encoding="utf-8")
    Path(str(fasta) + ".fai").write_text("chr1\t1\t6\t1\t2\n", encoding="utf-8")
    fasta.with_suffix(".dict").write_text("@HD\tVN:1.6\n", encoding="utf-8")
    gtf = tmp_path / "gencode.gtf"; gtf.write_text("# gtf\n", encoding="utf-8")
    vep = tmp_path / "vep"; vep.mkdir()
    proteome = tmp_path / "normal.fa"; proteome.write_text(">P\nA\n", encoding="utf-8")
    normal_expression = tmp_path / "normal_expression.tsv"
    normal_expression.write_text("gene\tmax_tpm\nGENE1\t0\n", encoding="utf-8")
    manifest = tmp_path / "references.json"
    manifest.write_text(json.dumps({"references": {
        "reference_fasta": {"path": str(fasta)}, "gencode_gtf": {"path": str(gtf)},
        "vep_cache": {"path": str(vep)}, "normal_proteome": {"path": str(proteome)},
        "normal_expression": {"path": str(normal_expression)},
    }}), encoding="utf-8")
    status, requirements = _assess_tier("prediction", rows, manifest)
    assert status == "READY"
    Path(str(fasta) + ".fai").unlink()
    status, _ = _assess_tier("prediction", rows, manifest)
    assert status == "PARTIAL"


def test_prediction_tier_requires_all_production_presentation_tools():
    base_tools = ["python", "neoag", "neoag-skill", "vep", "netmhcpan", "mhcflurry", "pvacseq", "prime"]
    rows = [CheckRow("tool", name, "OK") for name in base_tools]
    status, requirements = _assess_tier("prediction", rows, None)
    missing_tools = {row["requirement"] for row in requirements if row["kind"] == "tool" and row["status"] == "MISSING"}
    assert status == "PARTIAL"
    assert {"netmhcstabpan", "netchop"} <= missing_tools

def test_prediction_tier_immunogenicity_support_is_advisory(tmp_path: Path):
    tools = [
        "python", "neoag", "neoag-skill", "vep", "netmhcpan", "mhcflurry",
        "netmhcstabpan", "netchop", "pvacseq",
    ]
    rows = [CheckRow("tool", name, "OK") for name in tools]
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\nA\n", encoding="utf-8")
    Path(str(fasta) + ".fai").write_text("chr1\t1\t6\t1\t2\n", encoding="utf-8")
    fasta.with_suffix(".dict").write_text("@HD\tVN:1.6\n", encoding="utf-8")
    gtf = tmp_path / "gencode.gtf"; gtf.write_text("# gtf\n", encoding="utf-8")
    vep = tmp_path / "vep"; vep.mkdir()
    proteome = tmp_path / "normal.fa"; proteome.write_text(">P\nA\n", encoding="utf-8")
    normal_expression = tmp_path / "normal_expression.tsv"
    normal_expression.write_text("gene\tmax_tpm\nGENE1\t0\n", encoding="utf-8")
    manifest = tmp_path / "references.json"
    manifest.write_text(json.dumps({"references": {
        "reference_fasta": {"path": str(fasta)}, "gencode_gtf": {"path": str(gtf)},
        "vep_cache": {"path": str(vep)}, "normal_proteome": {"path": str(proteome)},
        "normal_expression": {"path": str(normal_expression)},
    }}), encoding="utf-8")
    status, requirements = _assess_tier("prediction", rows, manifest)
    advisory = [row for row in requirements if row["requirement"] == "immunogenicity_support"]
    assert status == "READY"
    assert advisory and advisory[0]["required"] == "false"
    assert advisory[0]["status"] == "WARN"

def test_install_check_writes_comprehensive_local_manifests(tmp_path: Path):
    layout = RunLayout.create(tmp_path / "install")
    outputs = _write_local_manifest_templates(layout, Path.cwd(), {"deploy_root": tmp_path / "neoag"})
    tools = Path(outputs["tools_manifest_template"]).read_text(encoding="utf-8")
    refs = Path(outputs["reference_manifest_template"]).read_text(encoding="utf-8")
    assert all(name in tools for name in ["netmhcpan", "spechla", "purple", "easyfuse", "splicemutr", "bam_matcher"])
    assert all(name in refs for name in ["normal_ligandome", "normal_junctions", "ctat_genome_lib", "sequenza_gc_wiggle", "bam_matcher_loci"])
    assert str(tmp_path / "neoag/refs/data/ref/hg38") in refs
    # The local template must inherit the formal manifest's canonical paths
    # and directory markers rather than the legacy fallback paths.
    assert "normal_junctions.GRCh38.tsv.gz" in refs
    assert "/data/hla/spechla/db" in refs
    assert "marker: versionInfo.json" in refs
    assert "salmon_tx2gene:" in refs
    assert "/data/rna/rsem_reference/gencode_v49/gencode_v49_gene_first" in refs
    assert "marker: .grp" in refs
    assert "G1000_loci_hg38.txt" not in refs
    assert "G1000_alleles_hg38.txt" not in refs
    assert "snaf.workflow.yaml" not in refs
    assert "splicemutr.workflow.yaml" not in refs
    asset_manifest = Path(outputs["asset_manifest_local"])
    asset_text = asset_manifest.read_text(encoding="utf-8")
    assert "reference_fasta_fai" in asset_text
    assert "1e74081a49ceb9739cc14c812fbb8b3db978eb80ba8e5350beb80d8ad8dfef3b" in asset_text
    target_paths = [line.split("\t")[2] for line in asset_text.splitlines() if "\t" in line and not line.startswith("#")]
    assert all(not path.startswith("/srv/") for path in target_paths[1:])


def test_install_tier_assessment_marks_inaccessible_reference_missing(tmp_path: Path, monkeypatch):
    inaccessible = "/root/other-machine/data/ref/hg38/Homo_sapiens_assembly38.fasta"
    manifest = tmp_path / "references.json"
    manifest.write_text(json.dumps({"references": {
        "reference_fasta": {"path": inaccessible},
    }}), encoding="utf-8")
    original_exists = Path.exists

    def guarded_exists(path: Path) -> bool:
        if str(path) == inaccessible:
            raise PermissionError(inaccessible)
        return original_exists(path)

    monkeypatch.setattr(Path, "exists", guarded_exists)
    status, rows = _assess_tier("prediction", [], manifest)
    by_name = {row["requirement"]: row for row in rows}
    assert status == "BLOCKED"
    assert by_name["reference_fasta"]["status"] == "MISSING"


def test_auto_config_discovers_tools_references_and_templates(tmp_path: Path):
    project = tmp_path / "project"
    tools_root = tmp_path / "tools"
    refs_root = tmp_path / "refs"
    (project / "scripts").mkdir(parents=True)
    (project / "scripts/run_bam_matcher_pair.sh").write_text("#!/bin/sh\n", encoding="utf-8")
    (tools_root / "bin").mkdir(parents=True)
    for executable in ("HLA-LA.pl", "bam-matcher"):
        path = tools_root / "bin" / executable
        path.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
        path.chmod(0o755)
    graph = refs_root / "data/hla/PRG_MHC_GRCh38_withIMGT"
    graph.mkdir(parents=True)
    fasta = refs_root / "data/ref/hg38/Homo_sapiens_assembly38.fasta"
    fasta.parent.mkdir(parents=True)
    fasta.write_text(">chr1\nA\n", encoding="utf-8")
    loci = refs_root / "data/sample_identity/bam_matcher.common_snps.hg38.vcf"
    loci.parent.mkdir(parents=True)
    loci.write_text("##reference=GRCh38\n#CHROM\tPOS\tID\tREF\tALT\n", encoding="utf-8")
    tools_manifest = tmp_path / "tools.json"
    tools_manifest.write_text(json.dumps({"tools": {
        "hla_la": {"executable": "HLA-LA.pl"},
        "bam_matcher": {"executable": "bam-matcher"},
        "snaf": {"executable": "snaf"},
    }}), encoding="utf-8")
    refs_manifest = tmp_path / "references.json"
    refs_manifest.write_text(json.dumps({"genome_build": "GRCh38", "references": {
        "reference_fasta": {"path": "missing.fa"},
        "hla_la_graph": {"path": "missing-graph"},
        "bam_matcher_loci": {"path": "missing.vcf"},
        "snaf_workflow": {"path": "missing.yaml"},
    }}), encoding="utf-8")
    result = configure_machine(
        project_root=project, tools_manifest=tools_manifest, reference_manifest=refs_manifest,
        outdir=tmp_path / "configured", tools_root=tools_root, reference_root=refs_root,
        licensed_root=tmp_path / "licensed", publish_local=True,
    )
    rows = {(row["component_type"], row["component"]): row for row in result.rows}
    assert rows[("tool", "hla_la")]["status"] == "CONFIGURED"
    assert rows[("tool", "bam_matcher")]["status"] == "CONFIGURED"
    assert rows[("tool", "snaf")]["status"] == "UNAVAILABLE"
    templates = Path(result.outputs["command_templates"]).read_text(encoding="utf-8")
    assert "hla_la" in templates and "bam_matcher" in templates
    assert (project / "configs/local/tools_manifest.generated.yaml").is_file()
    assert Path(result.outputs["configuration_status"]).is_file()


def test_auto_config_skips_inaccessible_paths_from_another_machine(tmp_path: Path, monkeypatch):
    project = tmp_path / "project"; project.mkdir()
    tools_root = tmp_path / "tools"; tools_root.mkdir()
    refs_root = tmp_path / "refs"
    portable_fasta = refs_root / "data/ref/hg38/Homo_sapiens_assembly38.fasta"
    portable_fasta.parent.mkdir(parents=True)
    portable_fasta.write_text(">chr1\nA\n", encoding="utf-8")
    tools_manifest = tmp_path / "tools.json"
    tools_manifest.write_text(json.dumps({"tools": {}}), encoding="utf-8")
    refs_manifest = tmp_path / "references.json"
    inaccessible = "/root/other-machine/data/ref/hg38/Homo_sapiens_assembly38.fasta"
    refs_manifest.write_text(json.dumps({"genome_build": "GRCh38", "references": {
        "reference_fasta": {"path": inaccessible},
    }}), encoding="utf-8")

    original_exists = Path.exists

    def guarded_exists(path: Path) -> bool:
        if str(path) == inaccessible:
            raise PermissionError(inaccessible)
        return original_exists(path)

    monkeypatch.setattr(Path, "exists", guarded_exists)
    result = configure_machine(
        project_root=project, tools_manifest=tools_manifest, reference_manifest=refs_manifest,
        outdir=tmp_path / "configured", tools_root=tools_root, reference_root=refs_root,
        licensed_root=tmp_path / "licensed",
    )
    row = next(row for row in result.rows if row["component"] == "reference_fasta")
    assert row["status"] == "CONFIGURED"
    assert row["resolved_path"] == str(portable_fasta.resolve())


def test_auto_config_discovers_project_wrappers_conda_envs_and_portable_assets(tmp_path: Path):
    project = tmp_path / "project"
    tools_root = tmp_path / "miniforge"
    refs_root = tmp_path / "neodata4git"
    (project / "bin").mkdir(parents=True)
    wrapper = project / "bin/neoag-nextflow"
    wrapper.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    wrapper.chmod(0o755)
    for executable in ("STAR", "regtools", "salmon"):
        path = tools_root / "envs/neoag-fusion/bin" / executable
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
        path.chmod(0o755)
    prime = refs_root / "data/predictors/prime/PRIME"
    prime.parent.mkdir(parents=True)
    prime.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    prime.chmod(0o755)
    portable_assets = {
        "normal_proteome": "data/normal/proteome/Homo_sapiens.GRCh38.pep.all.fa",
        "normal_ligandome": "data/normal/ligandome/normal_ms_ligands.tsv",
        "normal_junctions": "data/normal/junctions/normal_junctions.GRCh38.tsv.gz",
        "easyfuse_ref": "data/ref/ctat/current",
    }
    for relative in portable_assets.values():
        path = refs_root / relative
        if path.suffix:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("fixture\n", encoding="utf-8")
        else:
            path.mkdir(parents=True, exist_ok=True)
    tools_manifest = tmp_path / "tools.json"
    tools_manifest.write_text(json.dumps({"tools": {
        "nextflow": {"executable": "/root/old/bin/nextflow"},
        "star": {"executable": "star"},
        "regtools": {"executable": "regtools-neoag"},
        "salmon": {"executable": "salmon"},
        "prime": {"executable": "PRIME"},
    }}), encoding="utf-8")
    refs_manifest = tmp_path / "references.json"
    refs_manifest.write_text(json.dumps({"genome_build": "GRCh38", "references": {
        name: {"path": f"/root/old/{name}"} for name in portable_assets
    }}), encoding="utf-8")
    result = configure_machine(
        project_root=project, tools_manifest=tools_manifest, reference_manifest=refs_manifest,
        outdir=tmp_path / "configured", tools_root=tools_root, reference_root=refs_root,
        licensed_root=tmp_path / "licensed",
    )
    rows = {(row["component_type"], row["component"]): row for row in result.rows}
    for tool in ("nextflow", "star", "regtools", "salmon", "prime"):
        assert rows[("tool", tool)]["status"] == "CONFIGURED"
    for reference in portable_assets:
        assert rows[("reference", reference)]["status"] == "CONFIGURED"
    assert rows[("tool", "nextflow")]["resolved_path"] == str(wrapper.resolve())
    assert rows[("tool", "star")]["resolved_path"].endswith("neoag-fusion/bin/STAR")
    assert rows[("tool", "prime")]["resolved_path"] == str(prime.resolve())


def test_auto_config_rejects_wrong_build_bam_matcher_panel(tmp_path: Path):
    project = tmp_path / "project"; project.mkdir()
    tools_root = tmp_path / "tools"; (tools_root / "bin").mkdir(parents=True)
    executable = tools_root / "bin/bam-matcher"
    executable.write_text("#!/bin/sh\n", encoding="utf-8"); executable.chmod(0o755)
    fasta = tmp_path / "GRCh38.fa"; fasta.write_text(">chr1\nA\n", encoding="utf-8")
    loci = tmp_path / "bam_matcher.hg19.vcf"
    loci.write_text("##reference=GRCh37\n#CHROM\tPOS\tID\tREF\tALT\n", encoding="utf-8")
    tools_manifest = tmp_path / "tools.json"
    tools_manifest.write_text(json.dumps({"tools": {"bam_matcher": {"executable": "bam-matcher"}}}), encoding="utf-8")
    refs_manifest = tmp_path / "refs.json"
    refs_manifest.write_text(json.dumps({"genome_build": "GRCh38", "references": {
        "reference_fasta": {"path": str(fasta)}, "bam_matcher_loci": {"path": str(loci)},
    }}), encoding="utf-8")
    result = configure_machine(
        project_root=project, tools_manifest=tools_manifest, reference_manifest=refs_manifest,
        outdir=tmp_path / "configured", tools_root=tools_root, reference_root=tmp_path,
        licensed_root=tmp_path / "licensed",
    )
    rows = {(row["component_type"], row["component"]): row for row in result.rows}
    assert rows[("reference", "bam_matcher_loci")]["status"] == "BUILD_MISMATCH"
    assert rows[("tool", "bam_matcher")]["status"] == "PARTIAL"


def test_local_asset_manifest_rewrites_source_and_target_roots(tmp_path: Path):
    source = tmp_path / "assets.tsv"
    source.write_text(
        "asset_name\tsource_path\ttarget_path\tkind\trequired\tsha256\tmarker\n"
        "ref\t/srv/neoag-assets/source/data/ref.fa\t/srv/neoag-assets/install/data/ref.fa\tfile\t1\t-\t-\n"
        "tool\t/srv/neoag-assets/source/tools/tool\t/srv/neoag-tools/tools/tool\tdir\t1\t-\tbin/tool\n"
        "licensed\t/srv/neoag-assets/source/licensed/tool\t/srv/neoag-licensed/tool\tdir\t0\t-\ttool\n",
        encoding="utf-8",
    )
    source_root = tmp_path / "source"
    (source_root / "data").mkdir(parents=True)
    (source_root / "data/ref.fa").write_text(">chr1\nA\n", encoding="utf-8")
    output = _rewrite_asset_manifest(
        source, tmp_path / "local.tsv", tools_root=tmp_path / "tools",
        reference_root=tmp_path / "refs", licensed_root=tmp_path / "licensed",
        source_root=source_root,
    )
    text = output.read_text(encoding="utf-8")
    assert str(source_root / "data/ref.fa") in text
    assert str(tmp_path / "refs/data/ref.fa") in text
    assert str(tmp_path / "tools/tools/tool") in text
    assert str(tmp_path / "licensed/tool") in text
    assert _required_asset_sources_missing(output, "") == [f"tool={source_root}/tools/tool"]
    assert _required_asset_sources_missing(output, "asset-host") == []


def test_deployment_profile_defaults_follow_requested_tier(tmp_path: Path):
    project = tmp_path / "project"
    script = project / ".agents/skills/neoag-remote-deploy/scripts/16_install_new_machine.sh"
    script.parent.mkdir(parents=True)
    script.write_text("#!/bin/sh\n", encoding="utf-8")
    layout = RunLayout.create(tmp_path / "run")
    full = _deployment_command({"deployment_tier": "full"}, project, layout, execute=False)
    core = _deployment_command({"deployment_tier": "core"}, project, layout, execute=False)
    explicit = _deployment_command(
        {"deployment_tier": "full", "installer_profile": "all-open"}, project, layout, execute=False,
    )
    assert "--all-open" in full
    assert "--minimal" in core
    assert "--all-open" in explicit


def test_standard_installer_profile_covers_open_neo_run_groups():
    """Skill1 full→--standard must install fusion + bam-matcher (not side-path only)."""
    script = Path(__file__).resolve().parents[1] / (
        ".agents/skills/neoag-remote-deploy/scripts/16_install_new_machine.sh"
    )
    text = script.read_text(encoding="utf-8")
    # Extract the --standard) INSTALL_TOOL_GROUPS=(...) line
    match = None
    for line in text.splitlines():
        if line.strip().startswith("--standard)"):
            match = line
            break
    assert match is not None
    for flag in (
        "--core-env", "--vep", "--gatk", "--immunogenicity", "--optitype",
        "--facets", "--splice", "--lohhla", "--fusion", "--bam-matcher", "--install-torch",
    ):
        assert flag in match, f"missing {flag} in --standard groups"
    assert "--ascat-pyclone" not in match


def test_install_skill_defaults_to_central_asset_server(monkeypatch):
    monkeypatch.delenv("OPEN_NEO_ASSET_SOURCE_HOST", raising=False)
    monkeypatch.delenv("OPEN_NEO_ASSET_SOURCE_ROOT", raising=False)
    defaults = _apply_default_asset_source({})
    assert defaults["asset_source_host"] == DEFAULT_ASSET_SOURCE_HOST
    assert defaults["asset_source_root"] == DEFAULT_ASSET_SOURCE_ROOT
    explicit = _apply_default_asset_source({
        "asset_source_host": "user@other-host",
        "asset_source_root": "/srv/other-assets",
    })
    assert explicit["asset_source_host"] == "user@other-host"
    assert explicit["asset_source_root"] == "/srv/other-assets"


def test_asset_manifest_placeholder_detection_does_not_require_conda(tmp_path: Path):
    manifest = tmp_path / "assets.tsv"
    manifest.write_text(
        "asset_name\tsource_path\ttarget_path\n"
        "reference\t/srv/neoag-assets/source/data/ref\t/srv/neoag-assets/install/data/ref\n",
        encoding="utf-8",
    )
    assert _asset_manifest_uses_placeholder_source(manifest) is True


def test_runtime_fixes_select_facets_capable_r_environment(tmp_path: Path):
    project = tmp_path / "project"
    conda_base = tmp_path / "miniforge3"
    rscript = conda_base / "envs/neoag-r/bin/Rscript"
    rscript.parent.mkdir(parents=True)
    rscript.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    rscript.chmod(0o755)

    outputs = _apply_production_runtime_fixes(
        {
            "deploy_root": str(tmp_path / "deploy"),
            "tools_root": str(tmp_path / "tools"),
            "reference_root": str(tmp_path / "refs"),
            "conda_base": str(conda_base),
        },
        project,
    )
    local_env = Path(outputs["runtime_env_overrides"])
    text = local_env.read_text(encoding="utf-8")
    assert f"FACETS_R_ENV_PREFIX={conda_base / 'envs/neoag-r'}" in text


def test_install_skill_propagates_explicit_conda_base(tmp_path):
    project = tmp_path / "project"
    script = project / ".agents/skills/neoag-remote-deploy/scripts/16_install_new_machine.sh"
    script.parent.mkdir(parents=True)
    script.write_text("#!/bin/sh\n", encoding="utf-8")
    layout = RunLayout.create(tmp_path / "run")
    command = _deployment_command(
        {"deployment_tier": "full", "conda_base": "/nas/apps/miniforge3"},
        project,
        layout,
        execute=False,
    )
    assert command[command.index("--conda-base") + 1] == "/nas/apps/miniforge3"


def test_install_skill_propagates_conda_package_cache_source(tmp_path):
    project = tmp_path / "project"
    script = project / ".agents/skills/neoag-remote-deploy/scripts/16_install_new_machine.sh"
    script.parent.mkdir(parents=True)
    script.write_text("#!/bin/sh\n", encoding="utf-8")
    layout = RunLayout.create(tmp_path / "run")
    command = _deployment_command(
        {"deployment_tier": "full", "conda_pkgs_source": "/mnt/assets/conda_pkgs"},
        project,
        layout,
        execute=False,
    )
    assert command[command.index("--conda-pkgs-source") + 1] == "/mnt/assets/conda_pkgs"


def test_install_cli_accepts_conda_package_cache_source():
    args = build_parser().parse_args([
        "install-check", "--outdir", "work/install",
        "--conda-pkgs-source", "/mnt/assets/conda_pkgs",
    ])
    assert args.conda_pkgs_source == "/mnt/assets/conda_pkgs"


def test_install_skill_propagates_claude_code_bootstrap(tmp_path):
    project = tmp_path / "project"
    script = project / ".agents/skills/neoag-remote-deploy/scripts/16_install_new_machine.sh"
    script.parent.mkdir(parents=True)
    script.write_text("#!/bin/sh\n", encoding="utf-8")
    layout = RunLayout.create(tmp_path / "run")
    command = _deployment_command(
        {
            "deployment_tier": "core",
            "install_claude_code": True,
            "claude_code_channel": "2.1.170",
            "allow_download": True,
        },
        project,
        layout,
        execute=True,
    )
    assert "--claude-code" in command
    assert command[command.index("--claude-code-channel") + 1] == "2.1.170"
    assert "--allow-download" in command


def test_install_cli_accepts_claude_code_options():
    args = build_parser().parse_args([
        "install-check", "--project-root", ".", "--outdir", "work/install",
        "--install-claude-code", "--claude-code-channel", "stable",
    ])
    assert args.install_claude_code is True
    assert args.claude_code_channel == "stable"


def test_install_skill_collects_claude_code_status(tmp_path):
    layout = RunLayout.create(tmp_path / "run")
    report = layout.root / "deployment/readme_tools/claude_code/claude_code_install_report.md"
    report.parent.mkdir(parents=True)
    report.write_text(
        "# Claude Code install report\n\nBinary: `/opt/user/.local/bin/claude`\n"
        "Version: `2.1.170 (Claude Code)`\n",
        encoding="utf-8",
    )
    status, outputs = _collect_claude_code_status(layout, requested=True, executed=True)
    assert status["status"] == "READY"
    assert status["version"] == "2.1.170 (Claude Code)"
    assert Path(outputs["claude_code_status"]).is_file()
    assert outputs["claude_code_install_report"] == str(report)


def test_install_skill_launcher_rejects_explicit_old_python(tmp_path):
    launcher = Path.cwd() / ".agents/skills/open-neo-install-check/scripts/run.sh"
    fake_python = tmp_path / "python"
    fake_python.write_text("#!/bin/sh\nexit 1\n", encoding="utf-8")
    fake_python.chmod(0o755)
    proc = subprocess.run(
        ["bash", str(launcher), "--help"],
        cwd=Path.cwd(),
        env={"PATH": os.environ.get("PATH", ""), "NEOAG_PYTHON": str(fake_python)},
        text=True,
        capture_output=True,
        check=False,
    )
    assert proc.returncode == 31
    assert "OPEN_NEO_PYTHON_UNSUPPORTED" in proc.stderr


def test_install_skill_asset_server_defaults_can_be_overridden_by_environment(monkeypatch):
    monkeypatch.setenv("OPEN_NEO_ASSET_SOURCE_HOST", "asset-user@asset-host")
    monkeypatch.setenv("OPEN_NEO_ASSET_SOURCE_ROOT", "/data/neoag-assets")
    configured = _apply_default_asset_source({})
    assert configured["asset_source_host"] == "asset-user@asset-host"
    assert configured["asset_source_root"] == "/data/neoag-assets"


def test_install_checkpoint_can_be_reused(tmp_path: Path):
    command = ["bash", "-c", "printf done"]
    checkpoint = tmp_path / "checkpoint.json"
    log = tmp_path / "install.log"
    status, _, failure = _run_deployment(command, tmp_path, log, checkpoint, timeout=30, resume=False)
    assert status == "PASS" and failure == ""
    status, _, failure = _run_deployment(command, tmp_path, log, checkpoint, timeout=30, resume=True)
    assert status == "REUSED" and failure == ""


def test_install_checkpoint_records_timeout(tmp_path: Path, monkeypatch):
    class TimeoutProcess:
        pid = -1
        returncode = None

        def __init__(self, *args, **kwargs):
            self.communicate_calls = 0

        def communicate(self, timeout=None):
            self.communicate_calls += 1
            if self.communicate_calls == 1:
                raise subprocess.TimeoutExpired(cmd=["bash", "installer.sh"], timeout=timeout)
            self.returncode = -15
            return "partial output", None

        def poll(self):
            return self.returncode

        def terminate(self):
            self.returncode = -15

        def kill(self):
            self.returncode = -9

        def wait(self, timeout=None):
            if self.returncode is None:
                self.returncode = -15
            return self.returncode

    monkeypatch.setattr(subprocess, "Popen", TimeoutProcess)
    checkpoint = tmp_path / "checkpoint.json"
    status, log, failure = _run_deployment(
        ["bash", "installer.sh"], tmp_path, tmp_path / "install.log", checkpoint,
        timeout=60, resume=False,
    )
    assert status == "FAILED"
    assert "timed out after 60 seconds" in failure
    assert Path(log).read_text(encoding="utf-8") == "partial output"
    assert json.loads(checkpoint.read_text(encoding="utf-8"))["status"] == "TIMEOUT"


def test_install_resume_requires_approval(tmp_path: Path):
    result = run_install_check({
        "project_root": str(Path.cwd()), "deployment_tier": "review",
        "mode": "resume", "release_audit": False, "outdir": str(tmp_path / "install"),
    })
    assert result["status"] == "APPROVAL_REQUIRED"
    assert "APPROVAL_REQUIRED" in result["blocking_issues"]


def test_open_neo_run_plan_writes_route_and_config(tmp_path: Path):
    result = run_open_neo({
        "mode": "plan",
        "project_root": str(Path.cwd()),
        "sample_id": "PLAN",
        "fusion_tsv": str(Path.cwd() / "data/fixtures/easyfuse_fusions.pass.tsv"),
        "splice_junction_tsv": str(Path.cwd() / "data/fixtures/regtools_splice_junctions.tsv"),
        "hla_alleles": ["HLA-A*02:01", "HLA-B*07:02"],
        "doctor": False,
        "outdir": str(tmp_path / "plan"),
    })
    assert result["status"] == "PASS"
    assert Path(result["outputs"]["route_plan"]).is_file()
    config = Path(result["outputs"]["generated_run_config"])
    assert config.is_file()
    assert 'entry_mode = "e2e"' in config.read_text(encoding="utf-8")
    assert Path(result["outputs"]["pipeline_plan"]).is_file()
    assert Path(result["outputs"]["consensus_tool_run_status"]).is_file()
    assert Path(result["outputs"]["consensus_tool_consensus_summary"]).is_file()


def test_open_neo_run_missing_hla_is_blocked(tmp_path: Path):
    result = run_open_neo({
        "mode": "plan",
        "project_root": str(Path.cwd()),
        "fusion_tsv": str(Path.cwd() / "data/fixtures/easyfuse_fusions.pass.tsv"),
        "doctor": False,
        "outdir": str(tmp_path / "plan"),
    })
    assert result["status"] == "BLOCKED"
    assert "HLA_MISSING" in result["blocking_issues"]


def test_input_directory_header_detection_and_bam_index_qc(tmp_path: Path):
    data_dir = tmp_path / "inputs"
    data_dir.mkdir()
    fusion = data_dir / "caller.tsv"
    fusion.write_text("gene1\tgene2\tjunction_reads\nA\tB\t5\n", encoding="utf-8")
    hla = data_dir / "typing.txt"
    hla.write_text("HLA-A*02:01 HLA-B*07:02\n", encoding="utf-8")
    manifest = tmp_path / "sample.yaml"
    manifest.write_text(f"case_id: AUTO\nsample_id: AUTO\ninputs:\n  input_dir: {data_dir}\n", encoding="utf-8")
    result = inspect_manifest(manifest, input_dir=data_dir, output_dir=tmp_path / "out")
    assert result.status == "PASS"
    assert result.inputs["fusion_tsv"] == str(fusion)
    assert result.inputs["hla_file"] == str(hla)


def test_bam_without_index_blocks_production_route(tmp_path: Path):
    bam = tmp_path / "tumor.bam"
    bam.write_bytes(b"BAM")
    manifest = tmp_path / "sample.yaml"
    manifest.write_text(f"case_id: BAM\nsample_id: BAM\ninputs:\n  tumor_dna_bam: {bam}\n", encoding="utf-8")
    result = inspect_manifest(manifest, output_dir=tmp_path / "out")
    assert result.status == "BLOCKED"
    assert any(item["field"] == "tumor_dna_bam_index" for item in result.missing)


def test_tool_consensus_emits_domain_outputs(tmp_path: Path):
    optitype = tmp_path / "sample_optitype_result.txt"
    spechla = tmp_path / "sample_spechla_result.txt"
    optitype.write_text("HLA-A*02:01 HLA-B*07:02 HLA-C*07:02\n", encoding="utf-8")
    spechla.write_text("HLA-A*02:01 HLA-B*07:02 HLA-C*07:02\n", encoding="utf-8")
    outputs = build_tool_consensus({
        "sample_id": "S1",
        "tool_results": {"hla_typing": {"optitype": str(optitype), "spechla": str(spechla)}},
    }, tmp_path / "consensus")
    assert Path(outputs["hla_typing_consensus.tsv"]).is_file()
    text = Path(outputs["tool_consensus_summary.tsv"]).read_text(encoding="utf-8")
    assert "hla_typing\tCONSISTENT" in text
    assert Path(outputs["tool_evidence.long.tsv"]).is_file()


def test_rna_preprocessing_plans_gene_transcript_tpm_and_alt_vaf(tmp_path: Path):
    fastq1 = tmp_path / "tumor_R1.fastq.gz"
    fastq2 = tmp_path / "tumor_R2.fastq.gz"
    bam = tmp_path / "tumor.rna.bam"
    vcf = tmp_path / "somatic.vcf"
    tx2gene = tmp_path / "tx2gene.tsv"
    salmon_index = tmp_path / "salmon_index"
    salmon_index.mkdir()
    for path in (fastq1, fastq2, bam):
        path.write_bytes(b"fixture")
    tx2gene.write_text("transcript_id\tgene_id\nTX1\tG1\n", encoding="utf-8")
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t1\t.\tA\tT\t.\tPASS\t.\n", encoding="utf-8")
    result = prepare_rna_evidence(
        {
            "sample_id": "RNA1",
            "tumor_rna_fastq": [str(fastq1), str(fastq2)],
            "tumor_rna_bam": str(bam),
            "somatic_vcf": str(vcf),
            "salmon_index": str(salmon_index),
            "tx2gene": str(tx2gene),
        },
        project_root=Path.cwd(),
        outdir=tmp_path / "rna",
        execute=False,
    )
    stages = {row["stage"]: row for row in result["stages"]}
    assert stages["rna_quantification"]["status"] == "PLANNED"
    assert "run_salmon_fastq_to_tpm.sh" in stages["rna_quantification"]["command_preview"]
    assert stages["rna_alt_vaf"]["status"] == "PLANNED"
    assert "rna_allele_counts_pysam.py" in stages["rna_alt_vaf"]["command_preview"]
    assert Path(tmp_path / "rna/rna_preprocessing_status.tsv").is_file()


def _rna_profile_inputs(tmp_path: Path) -> dict[str, object]:
    files: dict[str, Path] = {}
    for name in ("tumor_R1.fastq.gz", "tumor_R2.fastq.gz", "GRCh38.fa", "gencode.gtf", "tx2gene.tsv"):
        path = tmp_path / name
        path.write_bytes(b"fixture")
        files[name] = path
    for name in ("star_index", "easyfuse_ref", "ctat", "salmon_index"):
        path = tmp_path / name
        path.mkdir()
        files[name] = path
    (files["easyfuse_ref"] / "fusioncatcher_index").mkdir()
    rsem_prefix = tmp_path / "rsem" / "reference"
    rsem_prefix.parent.mkdir()
    Path(str(rsem_prefix) + ".grp").write_text("fixture\n", encoding="utf-8")
    snaf_db = tmp_path / "snaf_db"
    (snaf_db / "controls").mkdir(parents=True, exist_ok=True)
    (snaf_db / "Alt91_db").mkdir(exist_ok=True)
    (snaf_db / "controls/GTEx_junction_counts.h5ad").write_bytes(b"fixture")
    splicemutr_workflow = tmp_path / "run_splicemutr.smk"
    splicemutr_workflow.write_text("rule all:\n    input: []\n", encoding="utf-8")
    return {
        "sample_id": "RNA_PROFILE",
        "tumor_rna_fastq": [str(files["tumor_R1.fastq.gz"]), str(files["tumor_R2.fastq.gz"])],
        "hla_alleles": ["HLA-A*02:01", "HLA-B*07:02"],
        "reference_fasta": str(files["GRCh38.fa"]),
        "gencode_gtf": str(files["gencode.gtf"]),
        "star_index": str(files["star_index"]),
        "easyfuse_ref": str(files["easyfuse_ref"]),
        "ctat_genome_lib": str(files["ctat"]),
        "salmon_index": str(files["salmon_index"]),
        "rsem_reference": str(rsem_prefix),
        "snaf_db": str(snaf_db),
        "splicemutr_workflow": str(splicemutr_workflow),
        "tx2gene": str(files["tx2gene.tsv"]),
        "rna_threads": 8,
    }


def test_rna_fastq_profile_generator_emits_full_dag(tmp_path: Path):
    inputs = _rna_profile_inputs(tmp_path)
    assert is_rna_fastq_profile_candidate(inputs)
    result = generate_rna_fusion_splice_manifest(
        inputs,
        tmp_path / "profile.toml",
        project_root=Path.cwd(),
        outdir=tmp_path / "run",
    )
    assert result["ready_for_execute"] is True
    text = Path(result["manifest"]).read_text(encoding="utf-8")
    for stage in (
        "fastq_qc", "rna_alignment", "rna_expression", "rsem_expression_crosscheck", "easyfuse_discovery",
        "junction_extraction", "fusion_cross_validation",
        "snaf_discovery", "splicemutr_discovery", "fusion_peptide_generation",
        "splice_candidate_normalization",
    ):
        assert f"[stages.{stage}]" in text
    assert "run_star_rna_fastq.sh" in text
    assert "run_rsem_fastq_to_tpm.sh" in text
    assert "run_fusioncatcher_sample.sh" not in text
    assert "--fusioncatcher" not in text
    assert "[stages.star_fusion_discovery]" not in text
    assert "[stages.fusioncatcher_discovery]" not in text
    assert "[stages.arriba_discovery]" not in text
    assert "--star-fusion" not in text
    assert "--arriba" not in text
    assert "--easyfuse" in text
    assert "--outdir {outdir}/rna/rsem_expression" in text
    assert "normalize_rna_fusion_splice.py" in text
    assert (tmp_path / "rna_fusion_splice.hla.txt").is_file()
    parsed = load_production_manifest(result["manifest"])
    assert parsed["run"]["profile"] == "rna_fusion_splice_v1"
    planned = run_production(
        result["manifest"], project_root=Path.cwd(), outdir=tmp_path / "production-plan"
    )
    assert {stage.name for stage in planned.stages} >= {
        "rna_alignment", "easyfuse_discovery", "splice_candidate_normalization"
    }



def test_rna_fastq_profile_merges_multi_batch_fastq(tmp_path: Path):
    inputs = _rna_profile_inputs(tmp_path)
    for name in ("tumor_batch2_R1.fastq.gz", "tumor_batch2_R2.fastq.gz"):
        path = tmp_path / name
        path.write_bytes(b"fixture")
    inputs["tumor_rna_fastq"] = [
        str(tmp_path / "tumor_R1.fastq.gz"),
        str(tmp_path / "tumor_batch2_R1.fastq.gz"),
        str(tmp_path / "tumor_R2.fastq.gz"),
        str(tmp_path / "tumor_batch2_R2.fastq.gz"),
    ]
    result = generate_rna_fusion_splice_manifest(
        inputs, tmp_path / "multi-batch-profile.toml", project_root=Path.cwd(), outdir=tmp_path / "run",
    )
    assert result["ready_for_execute"] is True
    text = Path(result["manifest"]).read_text(encoding="utf-8")
    assert "[stages.fastq_merge]" in text
    assert "run_merge_paired_fastq.sh" in text
    assert text.count("{outdir}/rna/merged_fastq/RNA_PROFILE_R1.fq.gz") >= 4
    assert text.count("{outdir}/rna/merged_fastq/RNA_PROFILE_R2.fq.gz") >= 4
    assert "depends_on = [\"fastq_merge\"]" in text


def test_rna_fastq_profile_accepts_rsem_expression_reference(tmp_path: Path):
    inputs = _rna_profile_inputs(tmp_path)
    inputs.pop("salmon_index")
    inputs.pop("tx2gene")
    prefix = Path(str(inputs["rsem_reference"]))
    inputs["rna_quant_method"] = "rsem"
    result = generate_rna_fusion_splice_manifest(
        inputs, tmp_path / "rsem-profile.toml", project_root=Path.cwd(), outdir=tmp_path / "run",
    )
    assert result["ready_for_execute"] is True
    text = Path(result["manifest"]).read_text(encoding="utf-8")
    assert "run_rsem_fastq_to_tpm.sh" in text
    assert "run_fusioncatcher_sample.sh" not in text
    assert "--fusioncatcher" not in text
    assert "[stages.star_fusion_discovery]" not in text
    assert "[stages.fusioncatcher_discovery]" not in text
    assert "[stages.arriba_discovery]" not in text
    assert "--star-fusion" not in text
    assert "--arriba" not in text
    assert "--easyfuse" in text
    assert "run_salmon_fastq_to_tpm.sh" not in text


def test_rna_fastq_profile_uses_builtin_snaf_when_reference_is_configured(tmp_path: Path):
    inputs = _rna_profile_inputs(tmp_path)
    snaf_db = tmp_path / "snaf_db"
    (snaf_db / "controls").mkdir(parents=True, exist_ok=True)
    (snaf_db / "Alt91_db").mkdir(exist_ok=True)
    for relative in (
        "controls/GTEx_junction_counts.h5ad",
        "Alt91_db/Hs_Ensembl_exon_add_col.txt",
        "Alt91_db/mRNA-ExonIDs.txt",
        "Alt91_db/Hs_gene-seq-2000_flank.fa",
    ):
        (snaf_db / relative).write_bytes(b"fixture")
    inputs["snaf_db"] = str(snaf_db)
    inputs["snaf_python"] = "/opt/snaf/bin/python"
    inputs["altanalyze_image"] = "neoag-altanalyze:snaf"
    result = generate_rna_fusion_splice_manifest(
        inputs, tmp_path / "snaf-profile.toml", project_root=Path.cwd(), outdir=tmp_path / "run",
    )
    text = Path(result["manifest"]).read_text(encoding="utf-8")
    assert "run_snaf_pipeline.sh" in text
    assert str(snaf_db) in text
    assert "neoag-altanalyze:snaf" in text


def test_open_neo_run_auto_generates_rna_profile_in_plan_mode(tmp_path: Path):
    inputs = _rna_profile_inputs(tmp_path)
    result = run_open_neo({
        **inputs,
        "mode": "plan",
        "doctor": False,
        "project_root": str(Path.cwd()),
        "outdir": str(tmp_path / "openneo"),
    })
    assert result["status"] == "PASS"
    manifest = Path(result["outputs"]["generated_production_manifest"])
    assert manifest.is_file()
    requirements = Path(result["outputs"]["rna_fusion_splice_requirements"])
    assert requirements.is_file()
    assert "rna_fusion_splice_v1" in manifest.read_text(encoding="utf-8")


def test_open_neo_run_blocks_multiple_rna_allele_evidence_modes(tmp_path: Path):
    rna_bam = tmp_path / "rna.bam"
    rna_bam.write_bytes(b"bam")
    rna_vaf = tmp_path / "rna_vaf.tsv"
    rna_vaf.write_text("event_id\trna_vaf\nE1\t0.1\n", encoding="utf-8")
    result = run_open_neo({
        "sample_id": "RNA_CONFLICT",
        "tumor_rna_bam": str(rna_bam),
        "rna_evidence_tsv": str(rna_vaf),
        "mode": "plan",
        "doctor": False,
        "project_root": str(Path.cwd()),
        "outdir": str(tmp_path / "openneo"),
    })
    assert result["status"] == "BLOCKED"
    assert FailureCode.AMBIGUOUS_INPUT.value in result["blocking_issues"]
    assert "RNA allele-evidence input mode" in result["steps"][0]["detail"]


def _automatic_plan_inputs(tmp_path: Path) -> tuple[dict[str, object], Path, Path]:
    tumor = tmp_path / "tumor.bam"; tumor.write_bytes(b"bam")
    normal = tmp_path / "normal.bam"; normal.write_bytes(b"bam")
    vcf = tmp_path / "somatic.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1\t.\tA\tT\t.\tPASS\t.\n", encoding="utf-8")
    hla = tmp_path / "hla.txt"; hla.write_text("HLA-A*02:01\nHLA-B*07:02\n", encoding="utf-8")
    fasta = tmp_path / "GRCh38.fa"; fasta.write_text(">chr1\nA\n", encoding="utf-8")
    facets = tmp_path / "common.vcf.gz"; facets.write_bytes(b"fixture")
    sequenza_gc = tmp_path / "gc50.wig.gz"; sequenza_gc.write_bytes(b"fixture")
    spechla_db = tmp_path / "spechla_db"; spechla_db.mkdir()
    hla_la_graph = tmp_path / "hla_la_graph"; hla_la_graph.mkdir()
    purple_ref = tmp_path / "purple_reference"; purple_ref.mkdir()
    tools = tmp_path / "tools.json"
    tools.write_text(json.dumps({"tools": {
        "samtools": {"executable": "/bin/true"},
        "facets": {"executable": "/bin/true"}, "sequenza": {"executable": "/bin/true"},
        "spechla": {"executable": "/bin/true"}, "purple": {"executable": "/bin/true"},
        "optitype": {"executable": "/bin/true"}, "hla_la": {"executable": "/bin/true"},
        "lohhla": {"executable": "/bin/true"}, "netmhcpan": {"executable": "/bin/true"},
        "mhcflurry": {"executable": "/bin/true"}, "netmhcstabpan": {"executable": "/bin/true"},
        "netchop": {"executable": "/bin/true"},
    }}), encoding="utf-8")
    refs = tmp_path / "refs.json"
    refs.write_text(json.dumps({"references": {
        "reference_fasta": {"path": str(fasta)}, "facets_snp_vcf": {"path": str(facets)},
        "sequenza_gc_wiggle": {"path": str(sequenza_gc)},
        "spechla_db": {"path": str(spechla_db)}, "purple_reference": {"path": str(purple_ref)},
        "hla_la_graph": {"path": str(hla_la_graph)},
    }}), encoding="utf-8")
    return {
        "sample_id": "AUTO1", "tumor_dna_bam": str(tumor), "normal_dna_bam": str(normal),
        "somatic_vcf": str(vcf), "hla_file": str(hla),
    }, tools, refs


def test_pipeline_plan_does_not_hash_large_bam_contents(tmp_path: Path):
    bam = tmp_path / "large.bam"
    with bam.open("wb") as handle:
        handle.truncate(60 * 1024 * 1024)
    rows = _manifest_file_hashes({"inputs": {"tumor_dna_bam": str(bam)}})
    assert rows[0]["sha256"] == "not_computed_large_file"
    assert rows[0]["size_bytes"] == str(60 * 1024 * 1024)
    assert rows[0]["mtime_ns"]


def test_capability_planner_builds_dna_hla_purity_loh_and_ranking_dag(tmp_path: Path):
    inputs, tools, refs = _automatic_plan_inputs(tmp_path)
    plan = build_automatic_production_plan(
        inputs, tmp_path / "manifests/automatic.toml", project_root=Path.cwd(),
        outdir=tmp_path / "run", tools_manifest=tools, reference_manifest=refs,
    )
    assert plan.status in {"READY", "PARTIAL"}
    text = Path(plan.manifest).read_text(encoding="utf-8")
    for stage in ("snv_indel_candidates", "purity_facets", "purity_sequenza", "purity_purple", "purity_consensus", "hla_loh_multi_tool"):
        assert f"[stages.{stage}]" in text
    assert "run_hla_loh_multi_tool.sh" in text
    assert "FACETS_CVAL_PRE=50" in text
    assert "FACETS_CVAL_PROC=300" in text
    assert "BIN_WINDOW=${SEQUENZA_BIN_WINDOW:-500}" in text
    assert 'hla_loh = "{outdir}/hla_loh/recommended_hla_loh.tsv"' in text
    assert "profiles/sarcoma_rna_supported_v2_provisional.toml" in text
    assert "configs/ranking/sarcoma_evidence_consensus_v3_source_chain.toml" in text
    assert 'required_presentation_predictors = ["netmhcpan", "mhcflurry", "netmhcstabpan", "netchop"]' in text
    assert set(["facets", "sequenza", "purple", "lohhla", "spechla", "netmhcpan", "mhcflurry"]) <= set(plan.selected_tools)
    rows = list(csv.DictReader(Path(plan.outputs["capability_decisions"]).open(), delimiter="\t"))
    assert any(row["domain"] == "purity_cnv" and row["status"] == "SELECTED" for row in rows)


def test_capability_planner_runs_three_hla_callers_from_normal_bam(tmp_path: Path):
    inputs, tools, refs = _automatic_plan_inputs(tmp_path)
    inputs["hla_file"] = str(tmp_path / "missing_hla.txt")
    plan = build_automatic_production_plan(
        inputs, tmp_path / "three-hla.toml", project_root=Path.cwd(), outdir=tmp_path / "run",
        tools_manifest=tools, reference_manifest=refs,
    )
    config = load_production_manifest(plan.manifest)
    assert {"hla_optitype", "hla_hla_la", "hla_spechla", "hla_consensus"} <= set(config["stages"])
    assert "run_optitype_sample.sh" in config["stages"]["hla_optitype"]["command"]
    assert "run_hla_la_sample.sh" in config["stages"]["hla_hla_la"]["command"]
    assert "run_spechla_sample.sh" in config["stages"]["hla_spechla"]["command"]
    assert any(row.tool == "provided_hla" and row.status == "INVALID" for row in plan.decisions)


def test_bam_matcher_parser_normalizes_match_mismatch_and_inconclusive(tmp_path: Path):
    header = "# BAM1\t BAM2\t DP_thresh\t FracCommon\t Same\t Same_hom\t Same_het\t Different\t 1het-2het\t 1het-2hom\t 1het-2sub\t 1hom-2het\t 1hom-2hom\t 1sub-2het\t Conclusion\n"
    cases = [("SAME", "MATCH", 0.99, 240, 2), ("DIFFERENT", "MISMATCH", 0.42, 100, 140), ("INCONCLUSIVE", "INSUFFICIENT_DATA", 0.0, 5, 4)]
    for index, (conclusion, expected, fraction, same, different) in enumerate(cases):
        path = tmp_path / f"identity{index}.tsv"
        path.write_text(header + f"normal.bam\ttumor.bam\t15\t{fraction}\t{same}\t1\t1\t{different}\t0\t0\t0\t0\t0\t0\t{conclusion}\n", encoding="utf-8")
        parsed = parse_bam_matcher_short(path)
        assert parsed["sample_identity_status"] == expected
        assert parsed["sites_compared"] == same + different


def test_capability_planner_gates_paired_analyses_on_bam_matcher(tmp_path: Path):
    inputs, tools, refs = _automatic_plan_inputs(tmp_path)
    tools_data = json.loads(tools.read_text(encoding="utf-8"))
    tools_data["tools"]["bam_matcher"] = {"executable": "/bin/true"}
    tools.write_text(json.dumps(tools_data), encoding="utf-8")
    loci = tmp_path / "identity.hg38.vcf"
    loci.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1\t.\tA\tT\t.\tPASS\t.\n", encoding="utf-8")
    refs_data = json.loads(refs.read_text(encoding="utf-8"))
    refs_data["references"]["bam_matcher_loci"] = {"path": str(loci)}
    refs.write_text(json.dumps(refs_data), encoding="utf-8")
    plan = build_automatic_production_plan(inputs, tmp_path / "auto.toml", project_root=Path.cwd(), outdir=tmp_path / "run", tools_manifest=tools, reference_manifest=refs)
    config = load_production_manifest(plan.manifest)
    assert "sample_identity_bam_matcher" in config["stages"]
    for stage in ("purity_facets", "purity_sequenza", "purity_purple", "hla_loh_multi_tool"):
        assert "sample_identity_bam_matcher" in config["stages"][stage]["depends_on"]
    assert any(row.domain == "sample_identity" and row.status == "SELECTED" for row in plan.decisions)


def test_tool_consensus_emits_sample_identity_status(tmp_path: Path):
    identity = tmp_path / "sample_identity.tsv"
    identity.write_text("sample_identity_status\tofficial_conclusion\tconfidence\tfraction_common\tsites_compared\nMATCH\tSAME\thigh\t0.99\t200\n", encoding="utf-8")
    outputs = build_tool_consensus({"tool_results": {"sample_identity": {"bam-matcher": str(identity)}}}, tmp_path / "consensus")
    rows = list(csv.DictReader(Path(outputs["sample_identity_consensus.tsv"]).open(), delimiter="\t"))
    assert rows[0]["sample_identity_status"] == "MATCH"


def test_open_neo_run_raw_bam_plan_uses_capability_aware_manifest(tmp_path: Path):
    inputs, tools, refs = _automatic_plan_inputs(tmp_path)
    result = run_open_neo({
        **inputs, "tools_manifest": str(tools), "reference_manifest": str(refs),
        "project_root": str(Path.cwd()), "doctor": False, "mode": "plan",
        "outdir": str(tmp_path / "openneo"),
    })
    assert result["status"] == "PASS"
    manifest = Path(result["outputs"]["automatic_production_manifest"])
    assert manifest.is_file()
    assert "[stages.purity_consensus]" in manifest.read_text(encoding="utf-8")
    assert Path(result["outputs"]["capability_plan"]).is_file()


def test_open_neo_run_clears_hla_missing_when_auto_typing_is_selected(tmp_path: Path):
    inputs, tools, refs = _automatic_plan_inputs(tmp_path)
    inputs["hla_file"] = str(tmp_path / "missing_hla.txt")
    result = run_open_neo({
        **inputs, "tools_manifest": str(tools), "reference_manifest": str(refs),
        "project_root": str(Path.cwd()), "doctor": False, "mode": "plan",
        "outdir": str(tmp_path / "openneo-auto-hla"),
    })
    assert result["status"] == "PASS"
    assert "hla_file" not in result["missing_evidence"]
    assert "hla_alleles_or_hla_file" not in result["missing_evidence"]


def test_capability_planner_combines_dna_and_rna_fastq_routes(tmp_path: Path):
    inputs, tools, refs = _automatic_plan_inputs(tmp_path)
    r1 = tmp_path / "rna_R1.fastq.gz"; r1.write_bytes(b"fixture")
    r2 = tmp_path / "rna_R2.fastq.gz"; r2.write_bytes(b"fixture")
    extra_refs = {
        "star_index": tmp_path / "star", "ctat_genome_lib": tmp_path / "ctat",
        "easyfuse_ref": tmp_path / "easyfuse", "salmon_index": tmp_path / "salmon",
    }
    for path in extra_refs.values(): path.mkdir()
    tx2gene = tmp_path / "tx2gene.tsv"; tx2gene.write_text("tx\tgene\n", encoding="utf-8")
    ref_data = json.loads(refs.read_text(encoding="utf-8"))
    for key, path in extra_refs.items(): ref_data["references"][key] = {"path": str(path)}
    ref_data["references"]["tx2gene"] = {"path": str(tx2gene)}
    refs.write_text(json.dumps(ref_data), encoding="utf-8")
    inputs["tumor_rna_fastq"] = [str(r1), str(r2)]
    plan = build_automatic_production_plan(inputs, tmp_path / "auto.toml", project_root=Path.cwd(), outdir=tmp_path / "run", tools_manifest=tools, reference_manifest=refs)
    text = Path(plan.manifest).read_text(encoding="utf-8")
    assert "[stages.snv_indel_candidates]" in text
    assert "[stages.rna_expression]" in text
    assert "[stages.easyfuse_discovery]" in text
    assert "[stages.rna_alt_vaf]" in text
    parsed = load_production_manifest(plan.manifest)
    assert parsed["stages"]["rna_alt_vaf"]["depends_on"] == ["rna_alignment"]
    assert parsed["evidence"]["rna_vaf"] == "{outdir}/rna/rna_alt_vaf.tsv"
    assert {"snv_indel", "rna_expression", "fusion", "splice"} <= set(plan.routes)


def test_rna_fastq_profile_reports_missing_required_assets(tmp_path: Path):
    fastq1 = tmp_path / "R1.fastq.gz"
    fastq2 = tmp_path / "R2.fastq.gz"
    fastq1.write_bytes(b"fixture")
    fastq2.write_bytes(b"fixture")
    result = generate_rna_fusion_splice_manifest(
        {"sample_id": "RNA_MISSING", "tumor_rna_fastq": [str(fastq1), str(fastq2)]},
        tmp_path / "profile.toml",
        project_root=Path.cwd(),
        outdir=tmp_path / "run",
    )
    assert result["ready_for_execute"] is False
    assert "hla_alleles_or_hla_file" in result["missing_required"]
    assert "star_index" in result["missing_required"]


def test_rna_fusion_cross_validation_marks_normal_background(tmp_path: Path):
    easyfuse = tmp_path / "fusions.pass.csv"
    easyfuse.write_text(
        "Fusion_Gene,junction_reads,frame\nEWSR1::WT1,20,in-frame\n",
        encoding="utf-8",
    )
    arriba = tmp_path / "arriba.tsv"
    arriba.write_text(
        "gene1\tgene2\tsplit_reads1\treading_frame\nEWSR1\tWT1\t12\tin-frame\n",
        encoding="utf-8",
    )
    normal = tmp_path / "normal.tsv"
    normal.write_text("gene1\tgene2\nEWSR1\tWT1\n", encoding="utf-8")
    outdir = tmp_path / "fusion-review"
    subprocess.run([
        sys.executable, str(Path.cwd() / "scripts/review_rna_fusions.py"),
        "--easyfuse", str(easyfuse), "--arriba", str(arriba),
        "--normal-readthrough", str(normal), "--outdir", str(outdir),
    ], check=True)
    rows = list(csv.DictReader((outdir / "fusion_consensus.tsv").open(), delimiter="\t"))
    assert rows[0]["n_tools"] == "2"
    assert rows[0]["normal_readthrough_status"] == "DETECTED"
    assert rows[0]["status"] == "NORMAL_BACKGROUND_REVIEW"


def test_direct_execute_requires_gateway(tmp_path: Path):
    comprehensive, weighted = _write_minimal_evidence(tmp_path)
    result = run_open_neo({
        "mode": "execute", "approved": True, "doctor": False,
        "comprehensive_evidence": str(comprehensive), "weighted_baseline": str(weighted),
        "outdir": str(tmp_path / "execute"),
    })
    assert result["status"] == "BLOCKED"
    assert "GATEWAY_REQUIRED" in result["blocking_issues"]


def _write_minimal_evidence(tmp_path: Path) -> tuple[Path, Path]:
    comprehensive = tmp_path / "comprehensive_peptide_evidence.tsv"
    comprehensive.write_text(
        "peptide_id\tevent_id\tevent_type\tpeptide\thla_allele\t"
        "cross_platform_status\trna_alt_reads\trna_vaf\t"
        "netmhcpan_mt_rank_el\tmhcflurry_presentation_score\t"
        "mutant_specificity_status\tccf_estimate\tccf_confidence\t"
        "appm_multiplier\thla_loh_status\trestricting_hla_lost\t"
        "safety_status\tsafety_evidence_completeness\n"
        "P1\tE1\tSNV\tSYFPEITHI\tHLA-A*02:01\t"
        "CROSS_PLATFORM_PASS_CONCORDANT\t8\t0.20\t0.2\t0.8\t"
        "MT_SPECIFIC\t0.9\thigh\t1.0\tRETAINED\tfalse\tPASS\tCOMPLETE\n",
        encoding="utf-8",
    )
    weighted = tmp_path / "ranked_peptides.tsv"
    weighted.write_text("peptide_id\tefficacy_score\tfinal_priority\nP1\t0.8\tB\n", encoding="utf-8")
    return comprehensive, weighted


def test_open_neo_run_ranking_only_keeps_baseline_and_emits_consensus(tmp_path: Path):
    comprehensive, weighted = _write_minimal_evidence(tmp_path)
    outdir = tmp_path / "ranking"
    result = run_open_neo({
        "mode": "ranking-only",
        "project_root": str(Path.cwd()),
        "comprehensive_evidence": str(comprehensive),
        "weighted_baseline": str(weighted),
        "doctor": False,
        "outdir": str(outdir),
    })
    assert result["status"] == "PASS"
    assert result["production_command"] == "neoag evidence-rank"
    assert weighted.read_text(encoding="utf-8").startswith("peptide_id")
    for name in (
        "ranked_peptides.weighted_baseline.tsv",
        "ranked_peptides.evidence_consensus.tsv",
        "ranked_events.evidence_consensus.tsv",
        "all_tool_results.tsv",
    ):
        assert (outdir / name).is_file()
    all_tool_header = (outdir / "all_tool_results.tsv").read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert "legacy_weighted_rank" in all_tool_header
    assert "evidence_rank" in all_tool_header
    assert (outdir / "pipeline/tool_consensus/presentation_consensus.tsv").is_file()


def _write_review_fixture(root: Path) -> None:
    scoring = root / "scoring"
    scoring.mkdir(parents=True)
    (scoring / "ranked_peptides.weighted_baseline.tsv").write_text(
        "peptide_id\tefficacy_score\tfinal_priority\nP1\t0.8\tB\nP2\t0.7\tC\n", encoding="utf-8"
    )
    (scoring / "ranked_peptides.evidence_consensus.tsv").write_text(
        "evidence_rank\tpeptide_id\tevent_id\tevent_type\tgene\tpeptide\thla_allele\tevidence_grade\tpareto_front\trna_support_state\tsafety_state\tpresentation_consensus_state\tmutant_specificity_state\tclonality_state\tccf_confidence\thla_appm_state\tevidence_completeness_state\thard_failure\n"
        "1\tP1\tE1\tSNV\tGENE1\tSYFPEITHI\tHLA-A*02:01\tR1\t1\tRNA_CONFIRMED\tSAFETY_PASS\tPRESENTATION_CONSISTENT_STRONG\tMT_SPECIFIC\tCLONAL_LIKE\thigh\tHLA_APPM_RETAINED\tCOMPLETE\tno\n"
        "2\tP2\tE2\tfusion\tGENE2::GENE3\tABCDEFGHI\tHLA-B*07:02\tR3\t1\tRNA_UNASSESSED\tSAFETY_PARTIAL\tPRESENTATION_SINGLE_TOOL\tNOVEL_JUNCTION\tunresolved\tlow\tHLA_LOH_UNASSESSED\tPARTIAL\tno\n",
        encoding="utf-8",
    )
    (scoring / "all_tool_results.tsv").write_text((scoring / "ranked_peptides.evidence_consensus.tsv").read_text(encoding="utf-8"), encoding="utf-8")
    (scoring / "validation_plan.tsv").write_text("event_id\trecommended_validation\nE1\tMT/WT pair\nE2\ttargeted RNA\n", encoding="utf-8")
    (scoring / "run_manifest.json").write_text(json.dumps({"schema_version": "1.0", "algorithm": "discrete_state_grade_track_pareto_v2", "status": "PROVISIONAL_RESEARCH_ONLY", "outputs": {}}) + "\n", encoding="utf-8")
    (scoring / "ranked_events.evidence_consensus.tsv").write_text(
        "event_evidence_rank\tevent_group_id\tevent_id\tgene\tevent_type\tbest_evidence_grade\tbest_pareto_front\tmanual_review_required\trecommended_next_steps\tevent_consensus_trace\t"
        "representative_1_peptide_id\trepresentative_1_peptide\trepresentative_1_hla_allele\n"
        "1\tEVENT:E1\tE1\tGENE1\tSNV\tR1\t1\tno\tMT/WT peptide assay\tcomplete\tP1\tSYFPEITHI\tHLA-A*02:01\n"
        "2\tEVENT:E2\tE2\tGENE2::GENE3\tFusion\tR3\t1\tyes\ttargeted RNA and second caller\tRNA-only\tP2\tABCDEFGHI\tHLA-B*07:02\n",
        encoding="utf-8",
    )


def test_open_neo_review_is_event_level_and_non_mutating(tmp_path: Path):
    result_dir = tmp_path / "result"
    _write_review_fixture(result_dir)
    weighted = result_dir / "scoring/ranked_peptides.weighted_baseline.tsv"
    source_hashes = {path: path.read_bytes() for path in (result_dir / "scoring").iterdir() if path.is_file()}
    outdir = tmp_path / "review"
    result = run_review({"result_dir": str(result_dir), "top_n": 2, "outdir": str(outdir)})
    assert result["status"] == "PASS_WITH_WARNINGS"
    assert weighted.read_bytes() == source_hashes[weighted]
    assert all(path.read_bytes() == content for path, content in source_hashes.items())
    review_rows = list(csv.DictReader((outdir / "review/candidate_review.tsv").open(), delimiter="\t"))
    assert [row["gene"] for row in review_rows] == ["GENE1", "GENE2::GENE3"]
    assert (outdir / "review/first_batch_experiment_set.tsv").is_file()
    first_batch = list(csv.DictReader((outdir / "review/first_batch_experiment_set.tsv").open(), delimiter="\t"))
    assert [row["gene"] for row in first_batch] == ["GENE1", "GENE2::GENE3"]
    assert first_batch[1]["experiment_priority"] == "TARGETED_RNA_FIRST"
    completion = list(csv.DictReader((outdir / "review/evidence_completion_queue.tsv").open(), delimiter="\t"))
    assert [row["gene"] for row in completion] == ["GENE2::GENE3"]
    patient_report = outdir / "reports/patient_report.html"
    assert patient_report.is_file()
    assert "患者沟通版" in patient_report.read_text(encoding="utf-8")
    assert not (outdir / "reports/patient_report.md").exists()
    assert not (outdir / "reports/patient_report.docx").exists()
    assert not (outdir / "reports/production_patient").exists()
    assert (outdir / "review/integrity/review_integrity.json").is_file()
    assert (outdir / "review/experiment_design/short_peptide_pool.tsv").is_file()
    assert (outdir / "review/experiment_design/targeted_rna_validation_plan.tsv").is_file()
    assert (outdir / "review/hla_loh_appm_review/appm_escape_review.md").is_file()
    assert (outdir / "review/ccf_clonality_review/ccf_clonality_review.md").is_file()
    assert (outdir / "reports/technical_report.md").is_file()
    assert review_rows[0]["pipeline_r_grade"] == "R1"
    assert review_rows[0]["experiment_priority"] == "EXPERIMENT_PRIORITY_HIGH"


def test_open_neo_review_can_skip_document_generation(tmp_path: Path):
    result_dir = tmp_path / "result"
    _write_review_fixture(result_dir)
    outdir = tmp_path / "review"
    result = run_review({"result_dir": str(result_dir), "reports": ["none"], "outdir": str(outdir)})
    assert result["status"] == "PASS_WITH_WARNINGS"
    assert not (outdir / "reports/patient_report.html").exists()
    assert not (outdir / "reports/technical_report.md").exists()
    assert next(step for step in result["steps"] if step["step_id"] == "06")["status"] == "SKIPPED"


def test_open_neo_review_rejects_unknown_report_type(tmp_path: Path):
    result = run_review({"result_dir": str(tmp_path / "result"), "reports": ["clinical"], "outdir": str(tmp_path / "review")})
    assert result["status"] == "BLOCKED"
    assert "INVALID_REPORT_SELECTION" in result["blocking_issues"]


def test_open_neo_review_blocks_incomplete_required_result_set(tmp_path: Path):
    result_dir = tmp_path / "result"
    _write_review_fixture(result_dir)
    (result_dir / "scoring/ranked_peptides.weighted_baseline.tsv").unlink()
    result = run_review({"result_dir": str(result_dir), "outdir": str(tmp_path / "review")})
    assert result["status"] == "BLOCKED"
    assert "REVIEW_INTEGRITY_BLOCKED" in result["blocking_issues"]


def test_open_neo_review_returns_needs_ranking_without_event_consensus(tmp_path: Path):
    result_dir = tmp_path / "result"
    _write_review_fixture(result_dir)
    (result_dir / "scoring/ranked_events.evidence_consensus.tsv").unlink()
    result = run_review({"result_dir": str(result_dir), "outdir": str(tmp_path / "review")})
    assert result["status"] == "NEEDS_RANKING"
    assert "NEEDS_RANKING" in result["blocking_issues"]


def test_open_neo_review_blocks_promoted_hard_failure(tmp_path: Path):
    result_dir = tmp_path / "result"
    _write_review_fixture(result_dir)
    peptide = result_dir / "scoring/ranked_peptides.evidence_consensus.tsv"
    text = peptide.read_text(encoding="utf-8").replace("\tCOMPLETE\tno\n", "\tCOMPLETE\tyes\n", 1)
    peptide.write_text(text, encoding="utf-8")
    (result_dir / "scoring/all_tool_results.tsv").write_text(text, encoding="utf-8")
    result = run_review({"result_dir": str(result_dir), "outdir": str(tmp_path / "review")})
    assert result["status"] == "BLOCKED"


def test_review_integrity_blocks_mass_r4_with_missing_core_presentation(tmp_path: Path):
    result_dir = tmp_path / "result"
    _write_review_fixture(result_dir)
    scoring = result_dir / "scoring"
    header = (
        "evidence_rank\tpeptide_id\tevent_id\tevent_type\tgene\tpeptide\thla_allele\t"
        "evidence_grade\tpareto_front\trna_support_state\tsafety_state\t"
        "presentation_consensus_state\tmutant_specificity_state\tclonality_state\t"
        "ccf_confidence\thla_appm_state\tevidence_completeness_state\thard_failure\t"
        "netmhcpan_mt_rank_el\tmhcflurry_presentation_score\n"
    )
    rows = "".join(
        f"{index + 1}\tP{index}\tE1\tSNV\tGENE1\tAAAAAAAAA\tHLA-A*02:01\t"
        "R4\t1\tRNA_UNASSESSED\tSAFETY_PARTIAL\tPRESENTATION_WEAK\tUNASSESSED\t"
        "UNRESOLVED\tlow\tHLA_LOH_UNASSESSED\tPARTIAL\tno\t99\t0\n"
        for index in range(20)
    )
    peptide_path = scoring / "ranked_peptides.evidence_consensus.tsv"
    peptide_path.write_text(header + rows, encoding="utf-8")
    (scoring / "all_tool_results.tsv").write_text(header + rows, encoding="utf-8")
    artifacts = {
        "run_manifest": str(scoring / "run_manifest.json"),
        "consensus_events": str(scoring / "ranked_events.evidence_consensus.tsv"),
        "consensus_peptides": str(peptide_path),
        "weighted_baseline": str(scoring / "ranked_peptides.weighted_baseline.tsv"),
        "all_tool_results": str(scoring / "all_tool_results.tsv"),
        "validation_plan": str(scoring / "validation_plan.tsv"),
    }
    result = audit_review_inputs(artifacts, tmp_path / "audit")
    assert result["status"] == "BLOCKED"
    checks = {row["check"]: row for row in result["checks"]}
    assert checks["evidence_grade_distribution"]["status"] == "WARN"
    assert checks["core_presentation_coverage"]["status"] == "FAIL"


def test_review_priority_uses_evidence_reason_codes_and_phase_deduplication():
    events = [
        {"event_evidence_rank": "1", "event_group_id": "G1", "event_id": "E1", "gene": "FUS", "event_type": "Fusion", "best_evidence_grade": "R3", "representative_1_peptide_id": "P1", "representative_1_peptide": "AAAA", "representative_1_hla_allele": "HLA-A*02:01"},
        {"event_evidence_rank": "2", "event_group_id": "G2", "event_id": "E2", "gene": "SPL", "event_type": "Splice", "best_evidence_grade": "R3", "phase_group_id": "PH1", "representative_1_peptide_id": "P2", "representative_1_peptide": "BBBB", "representative_1_hla_allele": "HLA-B*07:02"},
        {"event_evidence_rank": "3", "event_group_id": "G3", "event_id": "E3", "gene": "DRV", "event_type": "SNV", "best_evidence_grade": "R4", "manual_review_required": "yes", "phase_group_id": "PH1", "representative_1_peptide_id": "P3", "representative_1_peptide": "CCCC", "representative_1_hla_allele": "HLA-C*07:02"},
    ]
    peptides = [
        {"peptide_id": "P1", "rna_support_state": "RNA_UNASSESSED", "safety_state": "SAFETY_PARTIAL"},
        {"peptide_id": "P2", "rna_support_state": "RNA_UNASSESSED", "safety_state": "SAFETY_PARTIAL"},
        {"peptide_id": "P3", "rna_support_state": "RNA_CONFIRMED", "safety_state": "SAFETY_PASS"},
    ]
    all_tool = [{"peptide_id": "P1", "fusion_consensus_status": "SINGLE_CALLER_SUPPORTED"}]
    rows = build_review_rows(events, peptides, all_tool)
    assert rows[0]["experiment_priority"] == "FUSION_CONFIRMATION_FIRST"
    assert rows[0]["review_reason"] == "REVIEW_FUSION_SINGLE_CALLER"
    assert rows[1]["experiment_priority"] == "TARGETED_RNA_FIRST"
    assert rows[2]["experiment_priority"] == "MANUAL_REVIEW_ONLY"
    selected = select_first_batch(rows, 12)
    assert len([row for row in selected if row.get("phase_group_id") == "PH1"]) <= 1


def test_open_neo_review_blocks_source_output_mutation(tmp_path: Path):
    result_dir = tmp_path / "result"
    _write_review_fixture(result_dir)
    result = run_review({"result_dir": str(result_dir), "outdir": str(result_dir)})
    assert result["status"] == "BLOCKED"
    assert "REPORT_BOUNDARY_VIOLATION" in result["blocking_issues"]


def test_macro_skills_run_through_skill_runner(tmp_path: Path):
    comprehensive, weighted = _write_minimal_evidence(tmp_path)
    result = run_skill("open-neo-run", {
        "comprehensive_evidence": str(comprehensive),
        "weighted_baseline": str(weighted),
        "doctor": False,
        "outdir": str(tmp_path / "skill"),
    })
    assert result["status"] == "PASS"
    assert result["algorithm_owner"] == "src/neoag/evidence_consensus.py"
