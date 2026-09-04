from __future__ import annotations

import os
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
HELPER = ROOT / "scripts" / "lib" / "resolve_spechla_db.sh"


def _write_marker(db: Path) -> Path:
    marker = db / "ref" / "hla.ref.extend.fa"
    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.write_text(">hla\nACGT\n", encoding="utf-8")
    return db


def _resolve(env: dict[str, str], extra_root: Path | None = None) -> subprocess.CompletedProcess[str]:
    cmd = ["bash", "-lc", f"source {HELPER}; neoag_resolve_spechla_db {extra_root or ''}"]
    merged = os.environ.copy()
    for key in (
        "SPECHLA_DB",
        "SPECHLA_HOME",
        "NEOAG_REF_BUNDLE",
        "NEOAG_REFERENCE_ROOT",
        "OPEN_NEO_REFERENCE_ROOT",
        "REFERENCE_ROOT",
        "NEOAG_TOOLS_ROOT",
        "OPEN_NEO_TOOLS_ROOT",
        "NEOAG_PROJECT_ROOT",
    ):
        merged.pop(key, None)
    merged.update(env)
    return subprocess.run(cmd, text=True, capture_output=True, check=False, env=merged, cwd=ROOT)


def test_resolve_spechla_db_prefers_upstream_layout(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    db = _write_marker(bundle / "data" / "hla" / "spechla" / "db")
    proc = _resolve({"NEOAG_REFERENCE_ROOT": str(bundle)})
    assert proc.returncode == 0, proc.stderr
    assert Path(proc.stdout.strip()) == db.resolve()


def test_resolve_spechla_db_accepts_legacy_alias(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    db = _write_marker(bundle / "data" / "hla" / "spechla_db")
    proc = _resolve({"NEOAG_REFERENCE_ROOT": str(bundle)})
    assert proc.returncode == 0, proc.stderr
    assert Path(proc.stdout.strip()) == db.resolve()


def test_resolve_spechla_db_ignores_stale_env_and_finds_bundle(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    db = _write_marker(bundle / "data" / "hla" / "spechla" / "db")
    proc = _resolve({
        "SPECHLA_DB": str(tmp_path / "missing" / "spechla_db"),
        "NEOAG_REFERENCE_ROOT": str(bundle),
    })
    assert proc.returncode == 0, proc.stderr
    assert Path(proc.stdout.strip()) == db.resolve()
