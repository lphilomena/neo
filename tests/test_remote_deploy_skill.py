from __future__ import annotations

import csv
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SKILL = ROOT / ".agents" / "skills" / "neoag-remote-deploy"
SCRIPTS = SKILL / "scripts"


def _run(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(arg) for arg in args],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )


def test_remote_deploy_shell_scripts_have_valid_syntax() -> None:
    for script in sorted(SCRIPTS.glob("*.sh")):
        proc = _run("bash", "-n", script)
        assert proc.returncode == 0, f"{script}: {proc.stderr}"


def test_remote_deploy_skill_has_no_private_machine_defaults() -> None:
    prohibited = (
        "/home/na",
        "/mnt/zjl-bgi-zzb",
        "/root/neo",
        "10.200.50.134",
        "M1ML150017383",
        "chenxiaoliang",
        "xiaoliang",
    )
    for path in sorted(SKILL.rglob("*")):
        if not path.is_file():
            continue
        text = path.read_text(encoding="utf-8")
        for token in prohibited:
            assert token not in text, f"private default {token!r} found in {path}"


def test_new_machine_entrypoint_dry_run_needs_no_download_approval(tmp_path: Path) -> None:
    outdir = tmp_path / "run"
    proc = _run(
        "bash",
        SCRIPTS / "16_install_new_machine.sh",
        "--project-root",
        ROOT,
        "--tools-root",
        tmp_path / "tools",
        "--reference-root",
        tmp_path / "refs",
        "--licensed-root",
        tmp_path / "licensed",
        "--outdir",
        outdir,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    report = outdir / "new_machine_install_report.md"
    assert report.is_file()
    assert "DRY_RUN" in report.read_text(encoding="utf-8")
    assert "DOWNLOAD_NOT_APPROVED" not in proc.stdout + proc.stderr


def test_real_vcf_smoke_requires_explicit_inputs(tmp_path: Path) -> None:
    proc = _run(
        "bash",
        SCRIPTS / "16_install_new_machine.sh",
        "--project-root",
        ROOT,
        "--outdir",
        tmp_path / "run",
        "--no-sync-assets",
        "--no-runtime-validate",
        "--no-verify",
        "--run-real-vcf-smoke",
    )
    assert proc.returncode == 46
    assert "REAL_VCF_REQUIRED" in proc.stderr


def test_production_asset_manifest_keeps_pinned_reference_assets() -> None:
    manifest = ROOT / "configs" / "assets" / "production_assets.tsv"
    with manifest.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t"))
    by_name = {row["asset_name"]: row for row in rows}
    required = {
        "hla_ligand_atlas_peptides_raw",
        "iedb_human_ms_ligands",
        "iedb_mhc_ligand_manifest",
        "cancer_gene_list",
        "normal_proteome",
        "reference_fasta",
    }
    assert required <= by_name.keys()
    assert sum(row["sha256"] != "-" for row in rows) >= 11
    assert all("/mnt/zjl-bgi-zzb" not in row["source_path"] for row in rows)
