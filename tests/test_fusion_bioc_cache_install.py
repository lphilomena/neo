from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_easyfuse_native_is_the_default_without_duplicate_callers() -> None:
    script = (ROOT / "scripts/install_fusion_tools.sh").read_text(encoding="utf-8")

    assert 'INSTALL_MODE="${NEOAG_FUSION_INSTALL_MODE:-easyfuse}"' in script
    assert 'NEOAG_EASYFUSE_VERSION:-2.2.1' in script
    assert "install_easyfuse_native" in script
    assert "modules/arriba/environment.yml" in script
    assert "modules/starfusion/starfusion/environment.yml" in script
    assert "modules/fusioncatcher/environment.yml" in script
    assert 'case "${INSTALL_MODE}"' in script
    assert "standalone)" in script


def test_standalone_fallback_uses_resumable_bioconductor_cache() -> None:
    script = (ROOT / "scripts/install_fusion_tools.sh").read_text(encoding="utf-8")

    assert "with_bioc_data_cache.sh" in script
    assert "NEOAG_FUSION_BIOC_PACKAGE_KEY" in script
    assert "genomeinfodbdata-1.2.11" in script
    assert '"${transaction[@]}" env create' in script
    assert '"${transaction[@]}" env update' in script
    assert "--override-channels" not in script
    assert '[[ -x "${CONDA_BASE}/bin/conda" ]]' in script


def test_deploy_entrypoint_checks_easyfuse_instead_of_arriba() -> None:
    script = (ROOT / "scripts/deploy_external_tools.sh").read_text(encoding="utf-8")

    assert "easyfuse_installed()" in script
    assert "Install EasyFuse-native fusion stack" in script
    assert 'FUSION_SKIP="${SKIP_FUSION:-${SKIP_ARRIBA:-0}}"' in script
    assert "arriba_installed()" not in script


def test_remote_installer_exposes_explicit_standalone_fallback() -> None:
    script = (
        ROOT
        / ".agents/skills/neoag-remote-deploy/scripts/13_install_readme_tools.sh"
    ).read_text(encoding="utf-8")

    assert "--standalone-fusion" in script
    assert 'FUSION_INSTALL_MODE="${NEOAG_FUSION_INSTALL_MODE:-easyfuse}"' in script
    assert 'NEOAG_FUSION_INSTALL_MODE="${FUSION_INSTALL_MODE}"' in script


def test_easyfuse_runner_supports_module_native_layout() -> None:
    script = (ROOT / "scripts/run_easyfuse_sample.sh").read_text(encoding="utf-8")

    assert 'EASYFUSE_LAYOUT="module-native"' in script
    assert "sample\\tfastq_1\\tfastq_2" in script
    assert "EasyFuse module-native layout" in script
    assert "NEOAG_BIOC_CACHE_HELPER" in script
    assert "NEOAG_EASYFUSE_BIOC_PACKAGE_KEY" in script
    assert 'if [[ "${EASYFUSE_LAYOUT}" == "module-native" ]]' in script


def test_bioconductor_cache_download_resumes_and_rotates_mirrors() -> None:
    helper = (
        ROOT
        / ".agents/skills/neoag-remote-deploy/scripts/with_bioc_data_cache.sh"
    ).read_text(encoding="utf-8")

    assert "--continue=true" in helper
    assert "--lowest-speed-limit=32K" in helper
    assert 'for pass in 1 2; do' in helper
    assert 'rm -f "$TARBALL.aria2"' in helper
