from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_splice_installer_enables_pinned_snaf_by_default():
    script = (ROOT / "scripts/install_splice_tools.sh").read_text(encoding="utf-8")
    assert 'INSTALL_SNAF="${NEOAG_INSTALL_SNAF:-1}"' in script
    assert 'SNAF_ENV_NAME="${NEOAG_SNAF_ENV:-neoag-snaf}"' in script
    assert 'create -n "${SNAF_ENV_NAME}" -c conda-forge python=3.8 pip -y' in script
    assert "https://github.com/frankligy/SNAF.git" in script
    assert "e23ce39512a1a7f58c74e59b4b7cedc89248b908" in script
    assert "SNAF/archive/${SNAF_GIT_REF}.tar.gz" in script
    assert "https://gh-proxy.com/https://github.com/frankligy/SNAF/archive/" in script
    assert "https://mirrors.aliyun.com/pypi/simple" in script
    assert 'SNAF_ARCHIVE_CACHE="${SNAF_ARCHIVE_CACHE:-${NEOAG_SNAF_ARCHIVE_CACHE:-' in script
    assert 'curl -fL --retry 3 --connect-timeout 20' in script
    assert '"protobuf==3.20.3"' in script
    assert "SNAF import OK; TensorFlow" in script


def test_remote_install_skill_has_explicit_snaf_opt_out():
    script = (
        ROOT / ".agents/skills/neoag-remote-deploy/scripts/13_install_readme_tools.sh"
    ).read_text(encoding="utf-8")
    assert "INSTALL_SNAF=1" in script
    assert "--skip-snaf" in script
    assert 'NEOAG_INSTALL_SNAF="$INSTALL_SNAF"' in script


def test_splicemutr_is_pinned_and_installed_by_default():
    group = (ROOT / "scripts/install_splice_tools.sh").read_text(encoding="utf-8")
    installer = (ROOT / "scripts/install_splicemutr.sh").read_text(encoding="utf-8")
    env_yml = (ROOT / "conda/env.neoag-splicemutr.yml").read_text(encoding="utf-8")
    assert 'INSTALL_SPLICEMUTR="${NEOAG_INSTALL_SPLICEMUTR:-1}"' in group
    assert "scripts/install_splicemutr.sh" in group
    assert "ac0d17005cb37810bc1e6c9a50d7707f8bd3ae66" in installer
    assert "FertigLab/splicemutr/archive/${REF}.tar.gz" in installer
    assert "splicemutr-neoag" in installer
    assert "doctor" in installer
    assert "retry-all-errors" in installer
    assert "continue-at = -" in installer
    assert "GenomeInfoDbData_1.2.11.tar.gz" in installer
    assert "md5sum -c" in installer
    assert 'export BASH_ENV="${CURL_RETRY_HOME}/bash_env"' in installer
    assert "snakemake-minimal=7.32.4" in env_yml
    assert "bioconductor-genomicfeatures" in env_yml


def test_remote_install_skill_has_explicit_splicemutr_opt_out():
    script = (
        ROOT / ".agents/skills/neoag-remote-deploy/scripts/13_install_readme_tools.sh"
    ).read_text(encoding="utf-8")
    assert "INSTALL_SPLICEMUTR=1" in script
    assert "--skip-splicemutr" in script
    assert 'NEOAG_INSTALL_SPLICEMUTR="$INSTALL_SPLICEMUTR"' in script
