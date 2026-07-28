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
    assert {"netmhcpan_container_image", "netmhcstabpan_container_image"} <= by_name.keys()
    assert by_name["spechla_db"]["marker"] == "HLA/hla.ref.extend.fa"


def test_shared_asset_mode_creates_links_without_replacing_targets(tmp_path: Path) -> None:
    shared = tmp_path / "shared"
    source = shared / "data" / "ref" / "tiny.fa"
    source.parent.mkdir(parents=True)
    source.write_text(">tiny\nACGT\n", encoding="utf-8")
    manifest = tmp_path / "assets.tsv"
    manifest.write_text(
        "asset_name\tsource_path\ttarget_path\tkind\trequired\tsha256\tmarker\n"
        "tiny\t/srv/neoag-assets/source/data/ref/tiny.fa\t"
        "/srv/neoag-assets/install/data/ref/tiny.fa\tfile\t1\t-\t-\n",
        encoding="utf-8",
    )
    refs = tmp_path / "refs"
    proc = _run(
        "bash", SCRIPTS / "15_sync_asset_manifest.sh",
        "--project-root", ROOT,
        "--asset-manifest", manifest,
        "--shared-asset-root", shared,
        "--reference-root", refs,
        "--tools-root", tmp_path / "tools",
        "--licensed-root", tmp_path / "licensed",
        "--outdir", tmp_path / "out",
        "--execute",
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    target = refs / "data" / "ref" / "tiny.fa"
    assert target.is_symlink()
    assert target.resolve() == source

    proc = _run(
        "bash", SCRIPTS / "15_sync_asset_manifest.sh",
        "--project-root", ROOT,
        "--asset-manifest", manifest,
        "--shared-asset-root", shared,
        "--reference-root", refs,
        "--tools-root", tmp_path / "tools",
        "--licensed-root", tmp_path / "licensed",
        "--outdir", tmp_path / "out2",
        "--execute",
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "present" in proc.stdout


def test_bioconductor_cache_helper_is_wired_into_sequenza_and_ascat() -> None:
    helper = SCRIPTS / "with_bioc_data_cache.sh"
    assert helper.is_file()
    assert _run("bash", "-n", helper).returncode == 0
    readme_installer = (SCRIPTS / "13_install_readme_tools.sh").read_text(encoding="utf-8")
    ascat_installer = (ROOT / "scripts" / "install_ascat_pyclone.sh").read_text(encoding="utf-8")
    assert "genomeinfodbdata-1.2.9" in readme_installer
    assert "genomeinfodbdata-1.2.13" in ascat_installer


def test_shared_netmhcpan_asset_is_not_repaired_in_place() -> None:
    installer = (SCRIPTS / "13_install_readme_tools.sh").read_text(encoding="utf-8")
    assert '[[ -L "$LICENSED_ROOT/netMHCpan" ]]' in installer
    assert "skipping in-place native repair" in installer
