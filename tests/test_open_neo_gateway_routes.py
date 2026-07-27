import json

from neoag.controlled_execution.gateway import ROUTE_SPECS, _risk


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
    assert _risk("/open/review", {"reports": ["none"]}) == "LOW"
    assert _risk("/open/review", {"reports": ["patient", "technical"]}) == "MEDIUM"


def test_open_neo_gateway_accepts_rna_fusion_splice_profile_fields():
    optional = set(ROUTE_SPECS["/open/run"].optional)
    assert {
        "star_index", "ctat_genome_lib", "easyfuse_ref", "normal_readthrough",
        "snaf_workflow", "splicemutr_workflow", "rna_threads",
    } <= optional


def test_open_neo_gateway_manifest_risk_uses_declared_inputs(tmp_path):
    precomputed = tmp_path / "precomputed.json"
    precomputed.write_text(json.dumps({"inputs": {"somatic_vcf": "sample.vcf.gz"}}), encoding="utf-8")
    raw = tmp_path / "raw.json"
    raw.write_text(json.dumps({"tumor": {"wgs_bam": "tumor.bam"}}), encoding="utf-8")
    assert _risk("/open/run", {"mode": "execute", "sample_manifest": str(precomputed)}) == "MEDIUM"
    assert _risk("/open/run", {"mode": "execute", "sample_manifest": str(raw)}) == "HIGH"
