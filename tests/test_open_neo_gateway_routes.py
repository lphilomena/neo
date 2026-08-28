import json

from neoag.controlled_execution.gateway import ROUTE_SPECS, _resolve_project_root, _risk


def test_open_neo_gateway_routes_registered():
    assert "/open/install-check" in ROUTE_SPECS
    assert "/open/run" in ROUTE_SPECS
    assert "/open/review" in ROUTE_SPECS


def test_open_neo_gateway_mode_risk():
    assert _risk("/open/run", {"mode": "plan"}) == "LOW"
    assert _risk("/open/run", {"mode": "execute", "result_dir": "results/case"}) == "MEDIUM"
    assert _risk("/open/run", {"mode": "execute", "tumor_dna_bam": "tumor.bam"}) == "HIGH"
    assert _risk("/open/run", {"mode": "resume", "production_manifest": "production.toml"}) == "HIGH"
    assert _risk("/open/install-check", {"mode": "verify"}) == "LOW"
    assert _risk("/open/install-check", {"mode": "install"}) == "HIGH"
    assert _risk("/open/install-check", {"mode": "resume"}) == "HIGH"
    assert _risk("/open/review", {"reports": ["none"]}) == "LOW"
    assert _risk("/open/review", {"reports": ["patient", "technical"]}) == "MEDIUM"


def test_open_neo_gateway_accepts_rna_fusion_splice_profile_fields():
    optional = set(ROUTE_SPECS["/open/run"].optional)
    assert {
        "star_index", "ctat_genome_lib", "easyfuse_ref", "normal_readthrough",
        "snaf_workflow", "splicemutr_workflow", "rna_threads",
        "star_sjdb_overhang", "fusion_caller_root",
        "event_top_n", "candidate_top_n",
    } <= optional


def test_open_neo_gateway_accepts_all_open_run_local_input_and_tool_fields():
    optional = set(ROUTE_SPECS["/open/run"].optional)
    assert {
        "rna_fastq1", "rna_fastq2", "rna_bam", "rna_vaf",
        "evidence_consensus_rules", "predictor_deps", "netmhcpan_home",
        "netmhcstabpan_home", "samtools_executable", "star_executable",
        "star_index_build_dir", "easyfuse_star_index", "prime_evidence",
        "bigmhc_evidence", "deepimmuno_evidence",
    } <= optional


def test_open_neo_gateway_accepts_review_ranking_limits():
    optional = set(ROUTE_SPECS["/open/review"].optional)
    assert {"event_top_n", "candidate_top_n"} <= optional


def test_gateway_resolves_relative_project_root_against_configured_root(tmp_path):
    gateway_root = tmp_path / "checkout"
    assert _resolve_project_root(".", gateway_root) == gateway_root.resolve()
    assert _resolve_project_root("nested", gateway_root) == (gateway_root / "nested").resolve()


def test_open_neo_gateway_manifest_risk_uses_declared_inputs(tmp_path):
    precomputed = tmp_path / "precomputed.json"
    precomputed.write_text(json.dumps({"inputs": {"somatic_vcf": "sample.vcf.gz"}}), encoding="utf-8")
    raw = tmp_path / "raw.json"
    raw.write_text(json.dumps({"tumor": {"wgs_bam": "tumor.bam"}}), encoding="utf-8")
    assert _risk("/open/run", {"mode": "execute", "sample_manifest": str(precomputed)}) == "MEDIUM"
    assert _risk("/open/run", {"mode": "execute", "sample_manifest": str(raw)}) == "HIGH"
