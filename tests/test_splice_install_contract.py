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
