from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_extended_splice_installer_is_portable_and_pinned():
    script = (ROOT / "scripts/install_extended_splice_tools.sh").read_text()
    for value in ("3.1.1", "2.0.0", "7.1.1", "2.0.1"):
        assert value in script
    assert "NEOAG_ENV_ROOT" in script
    assert "NEOAG_SPLADDER_ENV" in script
    assert "NEOAG_IMMUNOPEPPER_ENV" in script
    assert "NEOAG_SPLICE_TOOLS_ROOT" in script
    assert "/mnt/zjl-bgi-zzb" not in script
    assert "/home/na" not in script
    assert 'rm -rf "${source_dir}"' not in script


def test_extended_splice_verifier_covers_all_tools():
    script = (ROOT / "scripts/verify_extended_splice_tools.sh").read_text()
    for tool in ("SplAdder", "ImmunoPepper", "pVACbind", "IRFinder-S"):
        assert tool in script


def test_extended_splice_conda_specs_exist():
    for filename in (
        "env.neoag-spladder.yml",
        "env.neoag-immunopepper.yml",
    ):
        assert (ROOT / "conda" / filename).is_file()
