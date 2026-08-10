from __future__ import annotations

import json
from hashlib import sha256
from pathlib import Path

from neoag.controlled_execution.doctor import run_doctor
from neoag.controlled_execution.pipeline_runner import run_pipeline_full, plan_pipeline
from neoag.controlled_execution.release_audit import scan_release_boundary


def test_doctor_minimal_dry_run(tmp_path):
    sample = tmp_path / "sample.json"
    sample.write_text(json.dumps({"sample_id": "T001", "inputs": {"ranked_peptides_recommendation": "missing.tsv"}}), encoding="utf-8")
    out = tmp_path / "doctor"
    res = run_doctor(project_root=Path.cwd(), outdir=out, sample_manifest=sample, release_audit=False, allow_execute=False)
    assert res.status in {"READY", "PARTIAL", "BLOCKED"}
    assert (out / "doctor_status.json").is_file()
    assert (out / "doctor_checks.tsv").is_file()


def test_pipeline_full_dry_run_outputs(tmp_path):
    sample = tmp_path / "sample.json"
    (tmp_path / "raw_events.tsv").write_text("event_id\tevent_type\nE1\tSNV\n", encoding="utf-8")
    (tmp_path / "raw_peptides.tsv").write_text("peptide_id\tevent_id\tpeptide\nP1\tE1\tAAAAAAAAA\n", encoding="utf-8")
    (tmp_path / "presentation_evidence.tsv").write_text("peptide_id\thla_allele\nP1\tHLA-A*02:01\n", encoding="utf-8")
    sample.write_text(json.dumps({
        "sample_id": "T002",
        "inputs": {
            "raw_events": str(tmp_path / "raw_events.tsv"),
            "raw_peptides": str(tmp_path / "raw_peptides.tsv"),
            "presentation_evidence": str(tmp_path / "presentation_evidence.tsv")
        }
    }), encoding="utf-8")
    out = tmp_path / "pipeline"
    run = run_pipeline_full(sample_manifest=sample, outdir=out, project_root=Path.cwd(), dry_run=True, run_doctor_first=False)
    assert run.status in {"DRY_RUN", "PARTIAL", "PASS"}
    assert (out / "run_manifest.json").is_file()
    assert (out / "pipeline_plan.md").is_file()
    assert any(s.name == "recommendation-ranking" for s in run.steps)


def test_release_audit_detects_private_path(tmp_path):
    root = tmp_path / "pkg"
    root.mkdir()
    (root / "x.py").write_text("PATH='/home/na/miniforge3/bin'\n", encoding="utf-8")
    res = scan_release_boundary(root)
    assert res["summary"]["private_path_hits"] >= 1
    assert res["status"] == "UNSAFE"


def test_release_audit_allows_deploy_prefixes_and_skips_vendor_trees(tmp_path):
    """Skill1 machine audit: scan deploy_root, allowlist install paths, skip vendor clones."""
    deploy = tmp_path / "neo-agent"
    (deploy / "env_tool" / "bin").mkdir(parents=True)
    (deploy / "env_tool" / "tools" / "upstream").mkdir(parents=True)
    wrapper = deploy / "env_tool" / "bin" / "vep"
    wrapper.write_text(
        "#!/usr/bin/env bash\n"
        f"ROOT={deploy / 'env_tool'}\n"
        "exec \"$ROOT/miniforge3/envs/neoag-vep/bin/vep\" \"$@\"\n",
        encoding="utf-8",
    )
    # Upstream docs under tools/ must not fail the deploy audit.
    (deploy / "env_tool" / "tools" / "upstream" / "README.md").write_text(
        "Built at /home/paul/localwork/bam-matcher\npatient sample demo\n",
        encoding="utf-8",
    )
    skip = {
        ".git", "work", "results", "tools", "conda_envs", "refs",
        "licensed_tools", "sources", "conda_pkgs", "pkgs",
    }
    res_ok = scan_release_boundary(
        deploy,
        skip_dirs=skip,
        allowed_path_prefixes=[str(deploy)],
    )
    assert res_ok["summary"]["private_path_hits"] == 0
    assert res_ok["summary"]["patient_hint_hits"] == 0
    assert res_ok["status"] in {"PASS", "PARTIAL"}

    (deploy / "env_tool" / "leak.py").write_text("x='/home/other/secret'\n", encoding="utf-8")
    res_leak = scan_release_boundary(
        deploy,
        skip_dirs=skip,
        allowed_path_prefixes=[str(deploy)],
    )
    assert any(h["path"].endswith("leak.py") for h in res_leak["private_path_hits"])


def test_plan_pipeline_ranks_result_review():
    steps = plan_pipeline({"sample_id": "T003", "inputs": {"ranked_peptides_recommendation": "a.tsv", "ranked_peptides_netmhcpan42": "b.tsv"}})
    names = [s.name for s in steps]
    assert "input-qc" in names
    assert "reports" in names



def test_release_audit_records_and_skips_generated_dirs(tmp_path):
    root = tmp_path / "pkg"
    (root / "results" / "nested").mkdir(parents=True)
    (root / "results" / "nested" / "x.py").write_text("PATH='/home/na/should_not_be_scanned_when_skipped'\n", encoding="utf-8")
    (root / "src").mkdir(parents=True)
    (root / "src" / "ok.py").write_text("print('ok')\n", encoding="utf-8")
    res = scan_release_boundary(root)
    assert res["summary"]["cache_hits"] >= 1
    assert res["summary"]["skipped_dirs"] >= 1
    assert res["summary"]["private_path_hits"] == 0


def test_doctor_mini_smoke_and_specific_checks(tmp_path):
    ref = tmp_path / "reference.json"
    ref.write_text(json.dumps({
        "reference_fasta": str(tmp_path / "missing.fa"),
        "gencode_gtf": str(tmp_path / "missing.gtf"),
        "vep_cache": str(tmp_path / "vep" / "homo_sapiens" / "115_GRCh38"),
        "facets_snp_vcf": str(tmp_path / "common_snp.vcf.gz"),
    }), encoding="utf-8")
    out = tmp_path / "doctor"
    res = run_doctor(project_root=Path.cwd(), outdir=out, reference_manifest=ref, release_audit=False, allow_execute=False, mini_smoke=True)
    names = {(r.category, r.name) for r in res.rows}
    assert ("mini_smoke", "neoag.cli") in names
    assert ("reference_specific", "reference_fasta.fai") in names
    assert ("reference_specific", "vep_cache_layout") in names
    assert ("entrypoint", "neoag-doctor") in names


def test_doctor_treats_top_level_reference_sha256_as_metadata(tmp_path):
    gtf = tmp_path / "gencode.gtf.gz"
    gtf.write_bytes(b"reference")
    ref = tmp_path / "reference.json"
    ref.write_text(json.dumps({
        "gencode_gtf": {
            "path": str(gtf),
            "sha256": sha256(b"reference").hexdigest(),
        },
    }), encoding="utf-8")
    result = run_doctor(
        project_root=Path.cwd(), outdir=tmp_path / "doctor-hash",
        reference_manifest=ref, release_audit=False, allow_execute=False,
    )
    checks = {(row.category, row.name): row.status for row in result.rows}
    assert checks[("reference", "gencode_gtf.path")] == "OK"
    assert checks[("reference_hash", "gencode_gtf.sha256")] == "OK"
    assert ("reference", "gencode_gtf.sha256") not in checks


def test_doctor_ignores_reference_descriptive_metadata_and_hash_placeholders(tmp_path):
    fasta = tmp_path / "GRCh38.fa"
    fasta.write_text(">1\nA\n", encoding="utf-8")
    ref = tmp_path / "reference.json"
    ref.write_text(json.dumps({"references": {"reference_fasta": {
        "path": str(fasta),
        "source": "Ensembl reference bundle",
        "version": "GRCh38",
        "marker": ".fai",
        "sha256": "-",
        "required": True,
    }}}), encoding="utf-8")
    result = run_doctor(
        project_root=Path.cwd(), outdir=tmp_path / "doctor-metadata",
        reference_manifest=ref, release_audit=False, allow_execute=False,
    )
    checks = {(row.category, row.name): row.status for row in result.rows}
    assert checks[("reference", "references.reference_fasta.path")] == "OK"
    assert not any(name.endswith((".source", ".version", ".marker")) for category, name in checks if category == "reference")
    assert ("reference_hash", "reference_fasta.sha256") not in checks


def test_doctor_uses_configured_java_executable_for_version_check(tmp_path):
    java = tmp_path / "java"
    java.write_text("#!/bin/sh\necho openjdk-test >&2\n", encoding="utf-8")
    java.chmod(0o755)
    tools = tmp_path / "tools.json"
    tools.write_text(json.dumps({"tools": {"java": {"executable": str(java)}}}), encoding="utf-8")
    result = run_doctor(
        project_root=Path.cwd(), outdir=tmp_path / "doctor-java",
        tools_manifest=tools, release_audit=False, allow_execute=True,
    )
    checks = {(row.category, row.name): row for row in result.rows}
    assert checks[("workflow", "java_version")].status == "OK"
    assert str(java) in checks[("workflow", "java_version")].path


def test_doctor_reports_inaccessible_cross_machine_paths_without_crashing(tmp_path, monkeypatch):
    inaccessible_tool = "/root/other-machine/licensed/netMHCpan"
    inaccessible_ref = "/root/other-machine/ref/GRCh38.fa"
    tools = tmp_path / "tools.json"
    tools.write_text(json.dumps({"tools": {"netmhcpan": {"executable": inaccessible_tool}}}), encoding="utf-8")
    refs = tmp_path / "refs.json"
    refs.write_text(json.dumps({"references": {"reference_fasta": {"path": inaccessible_ref}}}), encoding="utf-8")
    original_exists = Path.exists

    def guarded_exists(path: Path) -> bool:
        if str(path) in {inaccessible_tool, inaccessible_ref}:
            raise PermissionError(str(path))
        return original_exists(path)

    monkeypatch.setattr(Path, "exists", guarded_exists)
    result = run_doctor(
        project_root=Path.cwd(), outdir=tmp_path / "doctor-inaccessible",
        tools_manifest=tools, reference_manifest=refs, release_audit=False, allow_execute=False,
    )
    checks = {(row.category, row.name): row for row in result.rows}
    assert checks[("tool", "netmhcpan")].status in {"LICENSE_REQUIRED", "MISSING"}
    assert checks[("reference", "references.reference_fasta.path")].status == "MISSING"



def test_gateway_validation_and_risk_policy(tmp_path):
    from neoag.controlled_execution.gateway import _approval_required, _risk, _validate_request

    roots = [tmp_path]
    try:
        _validate_request("/ranking-compare", {"recommendation": "a.tsv"}, roots)
        assert False, "missing netmhcpan42 should fail"
    except Exception as exc:
        assert getattr(exc, "status", "") == "BAD_REQUEST"

    assert _risk("/patient-report", {}) == "MEDIUM"
    assert _risk("/pipeline-full", {"sample_manifest": "sample.json", "execute": True}) == "HIGH"
    assert _approval_required("/pipeline-full", {"sample_manifest": "sample.json", "execute": True}) is True
    assert _approval_required("/pipeline-full", {"sample_manifest": "sample.json", "execute": True, "approved": True}) is False

    try:
        _validate_request("/doctor", {"outdir": "/not/allowed/out"}, roots)
        assert False, "unsafe output path should fail"
    except Exception as exc:
        assert getattr(exc, "status", "") == "PATH_NOT_ALLOWED"


def test_gateway_job_store_persists_status(tmp_path):
    from neoag.controlled_execution.gateway import JobStore

    store = JobStore(tmp_path / "jobs")
    job_id = store.add({"route": "/health", "status": "QUEUED", "risk_level": "LOW", "request": {"case_id": "C1"}})
    store.update(job_id, status="PASS", result={"ok": True})
    loaded = store.get(job_id)
    assert loaded["status"] == "PASS"
    assert loaded["result"]["ok"] is True
    paths = store.paths(job_id)
    assert json.loads(paths["request"].read_text(encoding="utf-8"))["case_id"] == "C1"
    assert json.loads(paths["approval"].read_text(encoding="utf-8"))["approval_status"] == "not_required"
    assert json.loads(paths["status"].read_text(encoding="utf-8"))["status"] == "PASS"
    assert len(paths["audit"].read_text(encoding="utf-8").splitlines()) >= 2


def test_gateway_route_registry_includes_expected_routes():
    from neoag.controlled_execution.gateway import ROUTE_SPECS

    for route in ["/doctor", "/input-qc", "/ranking-compare", "/appm-review", "/ccf-review", "/patient-report", "/pipeline-full"]:
        assert route in ROUTE_SPECS
    assert "sample_manifest" in ROUTE_SPECS["/pipeline-full"].required
    assert "netmhcpan42" in ROUTE_SPECS["/ranking-compare"].required



def test_pipeline_full_declared_ranked_outputs_do_not_block(tmp_path):
    from neoag.controlled_execution.pipeline_runner import run_pipeline_full

    sample = tmp_path / "sample.json"
    sample.write_text(json.dumps({
        "sample_id": "T004",
        "inputs": {
            "ranked_peptides_recommendation": str(tmp_path / "missing.recommendation.tsv"),
            "ranked_peptides_netmhcpan42": str(tmp_path / "missing.netmhcpan42.tsv"),
            "evidence_report": str(tmp_path / "missing.html"),
        },
    }), encoding="utf-8")
    out = tmp_path / "pipeline"
    run = run_pipeline_full(sample_manifest=sample, outdir=out, project_root=Path.cwd(), dry_run=True, run_doctor_first=False)
    by_id = {s.step_id: s for s in run.steps}
    assert run.status == "DRY_RUN"
    assert by_id["05"].status == "UNASSESSED"
    assert by_id["05"].mode == "unassessed_existing_results"
    assert by_id["06"].status == "UNASSESSED"
    assert by_id["11"].status == "UNASSESSED"
    header = (out / "pipeline_status.tsv").read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert "requires_execution" in header


def test_pipeline_full_requires_execution_flag_for_vcf(tmp_path):
    from neoag.controlled_execution.pipeline_runner import plan_pipeline

    vcf = tmp_path / "tumor.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
    steps = plan_pipeline({"sample_id": "T005", "inputs": {"somatic_vcf": str(vcf)}})
    by_id = {s.step_id: s for s in steps}
    assert by_id["05"].mode == "requires_execution"
    assert by_id["05"].requires_execution is True



def test_gateway_request_metadata_has_mvp_audit_fields():
    from neoag.controlled_execution.gateway import _request_metadata

    meta = _request_metadata(
        "/pipeline-full",
        {"sample_manifest": "sample.yaml", "sample_id": "S1", "case_id": "C1", "user": "tester", "approved": True},
        job_id="job-1",
        risk="HIGH",
        status="QUEUED",
        output_dir="work/out",
    )
    for key in ["request_id", "case_id", "sample_id", "user", "timestamp", "input_manifest", "skill_or_pipeline", "command_preview", "risk_level", "approval_status", "job_id", "output_dir", "status", "failure_reason"]:
        assert key in meta
    assert meta["approval_status"] == "approved"
    assert meta["input_manifest"] == "sample.yaml"

    open_meta = _request_metadata(
        "/open/run", {"result_dir": "results/S1"}, job_id="job-2",
        risk="MEDIUM", status="QUEUED", output_dir="work/out",
    )
    assert open_meta["task_type"] == "open_run"
    assert open_meta["skill_or_pipeline"] == "open_run"
    assert open_meta["command_preview"] == "open-neo run"



def test_manifest_first_run_manifest_schema(tmp_path):
    from neoag.controlled_execution.pipeline_runner import run_pipeline_full

    sample = tmp_path / "sample.json"
    ref = tmp_path / "reference.json"
    tools = tmp_path / "tools.json"
    sample.write_text(json.dumps({
        "sample_id": "T006",
        "tumor": {"wgs_bam": str(tmp_path / "tumor.bam")},
        "normal": {"wgs_bam": str(tmp_path / "normal.bam")},
        "inputs": {"ranked_peptides_recommendation": str(tmp_path / "missing.tsv")},
    }), encoding="utf-8")
    ref.write_text(json.dumps({"genome_build": "GRCh38", "reference_fasta": str(tmp_path / "ref.fa"), "gencode_gtf": str(tmp_path / "g.gtf"), "vep_cache": str(tmp_path / "vep")}), encoding="utf-8")
    tools.write_text(json.dumps({"bwa": {"mode": "container", "image": "bwa@sha256:abc", "executable": "bwa"}, "netmhcpan": {"mode": "local_license", "path": "/opt/netMHCpan", "license_required": True}}), encoding="utf-8")
    out = tmp_path / "run"
    run = run_pipeline_full(sample_manifest=sample, reference_manifest=ref, tools_manifest=tools, outdir=out, project_root=Path.cwd(), dry_run=True, run_doctor_first=False)
    data = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
    for key in ["run_id", "sample_id", "pipeline_version", "git_sha", "input_hashes", "reference_hashes", "tool_versions", "container_digests", "start_time", "end_time", "status", "output_manifest"]:
        assert key in data
    assert data["sample_id"] == "T006"
    assert any(x["tool"] == "bwa" and x["container_digest"] == "abc" for x in data["container_digests"])
    assert (out / "manifest_validation.tsv").is_file()
    assert (out / "reference_hashes.tsv").is_file()
    assert (out / "tool_versions.tsv").is_file()
