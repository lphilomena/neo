from __future__ import annotations

import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SKILL = ROOT / ".agents" / "skills" / "neoag-remote-deploy"
SCRIPTS = SKILL / "scripts"
INSTALLER = SCRIPTS / "17_install_claude_code.sh"


def _run(*args: object) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(arg) for arg in args],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )


def test_claude_code_installer_dry_run_is_network_free(tmp_path: Path) -> None:
    proc = _run("bash", INSTALLER, "--outdir", tmp_path, "--channel", "stable")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    report = tmp_path / "claude_code_install_report.md"
    assert report.is_file()
    assert "DRY_RUN" in report.read_text(encoding="utf-8")
    assert "https://claude.ai/install.sh" in proc.stdout


def test_claude_code_execute_requires_download_approval(tmp_path: Path) -> None:
    proc = _run("bash", INSTALLER, "--outdir", tmp_path, "--execute")
    assert proc.returncode == 23
    assert "DOWNLOAD_NOT_APPROVED" in proc.stderr


def test_claude_code_channel_rejects_shell_input(tmp_path: Path) -> None:
    proc = _run("bash", INSTALLER, "--outdir", tmp_path, "--channel", "stable;false")
    assert proc.returncode == 2
    assert "CLAUDE_CODE_CHANNEL_INVALID" in proc.stderr


def test_claude_code_is_wired_into_both_macro_installers() -> None:
    readme_installer = (SCRIPTS / "13_install_readme_tools.sh").read_text(encoding="utf-8")
    machine_installer = (SCRIPTS / "16_install_new_machine.sh").read_text(encoding="utf-8")
    skill_text = (SKILL / "SKILL.md").read_text(encoding="utf-8")
    for text in (readme_installer, machine_installer, skill_text):
        assert "--claude-code" in text
        assert "--claude-code-channel" in text
    assert "17_install_claude_code.sh" in readme_installer
    assert "never performs login or stores credentials" in skill_text
