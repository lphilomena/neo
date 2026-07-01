from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

BENCHMARK_TESTS = {
    "test_benchmark_improve.py",
    "test_benchmark_system.py",
}

BENCHMARK_TEST_ITEMS = {
    "test_benchmark_v042_external_required_outputs",
}

EXTERNAL_TESTS = {
    "test_check_tools_env.py",
    "test_facets_lohhla.py",
    "test_immunogenicity_tools.py",
    "test_netmhcpan_local.py",
    "test_prime_parallel.py",
    "test_tools.py",
}

EXTERNAL_TEST_ITEMS = {
    "test_netmhcpan_tool_run_provenance_columns",
}

TIMEOUT_BY_CATEGORY = {
    "unit": 60,
    "integration": 120,
    "benchmark": 600,
    "external": 900,
}

TIMEOUT_EXPLANATION = {
    "unit": "fast unit contract",
    "integration": "focused fixture integration smoke",
    "benchmark": "synthetic/sensitivity benchmark smoke",
    "external": "external tool or local installation smoke",
}

INTEGRATION_TESTS = {
    "test_appm_escape_consistency.py",
    "test_intermediates.py",
    "test_rc_cli.py",
    "test_snv_phase1_wes.py",
    "test_splice_junction.py",
    "test_sv_phase1.py",
    "test_sv_score_v03.py",
    "test_sv_wes_phase1_5.py",
    "test_v03.py",
    "test_v041_appm_ccf_escape.py",
    "test_v04_evidence_safety_escape.py",
    "test_variant_peptide_upstream.py",
}


def pytest_addoption(parser: pytest.Parser) -> None:
    group = parser.getgroup("neoag")
    group.addoption("--run-integration", action="store_true", help="run integration tests")
    group.addoption("--run-benchmark", action="store_true", help="run benchmark tests")
    group.addoption("--run-external", action="store_true", help="run tests requiring external tools or installs")
    group.addoption("--run-all", action="store_true", help="run unit, integration, benchmark, and external tests")


def _category(path: Path, item_name: str = "") -> str:
    name = path.name
    if item_name in EXTERNAL_TEST_ITEMS:
        return "external"
    if item_name in BENCHMARK_TEST_ITEMS:
        return "benchmark"
    if name in BENCHMARK_TESTS:
        return "benchmark"
    if name in EXTERNAL_TESTS:
        return "external"
    if name in INTEGRATION_TESTS:
        return "integration"
    return "unit"


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    run_all = config.getoption("--run-all")
    timeout_plugin = config.pluginmanager.hasplugin("timeout")
    enabled = {
        "unit": True,
        "integration": run_all or config.getoption("--run-integration"),
        "benchmark": run_all or config.getoption("--run-benchmark"),
        "external": run_all or config.getoption("--run-external"),
    }
    skip_reasons = {
        "integration": (
            "integration test skipped by default: may run cross-module CLI/workflow smoke; "
            "enable with --run-integration or --run-all; timeout 120s when pytest-timeout is installed"
        ),
        "benchmark": (
            "benchmark test skipped by default: may run synthetic/sensitivity/external-required checks; "
            "enable with --run-benchmark or --run-all; timeout 600s when pytest-timeout is installed"
        ),
        "external": (
            "external-tool test skipped by default: requires installed/licensed tools, network, or local data; "
            "enable with --run-external or --run-all; timeout 900s when pytest-timeout is installed"
        ),
    }
    for item in items:
        category = _category(Path(str(item.fspath)), item.name)
        item.add_marker(getattr(pytest.mark, category))
        if timeout_plugin and not any(mark.name == "timeout" for mark in item.iter_markers()):
            item.add_marker(pytest.mark.timeout(TIMEOUT_BY_CATEGORY[category]))
        if not enabled.get(category, False):
            item.add_marker(pytest.mark.skip(reason=skip_reasons[category]))
