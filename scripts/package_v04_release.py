#!/usr/bin/env python3
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTDIR = ROOT / 'work' / 'releases'

INCLUDE_ROOTS = [
    'assets', 'bin', 'conda', 'conf', 'docs', 'modules', 'profiles', 'resources',
    'scripts', 'src', 'tests', 'workflows',
]
INCLUDE_FILES = [
    'pyproject.toml', 'README.md', 'RELEASE.md', 'RELEASE_V04_20260616.md',
    'RELEASE_REFRESH_20260615.md', 'CHANGELOG_V04_EVIDENCE_SAFETY_ESCAPE.md',
    'CHANGELOG_V041_APPM_CCF_IMMUNE_ESCAPE.md', 'CHANGELOG_V042.md',
    'CHANGELOG_V043.md', 'README_PATCH_V043.md', 'README_PATCH_V043_CCF21.md',
]
INCLUDE_DATA_ROOTS = [
    'data/fixtures', 'data/fixtures_snv', 'data/fixtures_sv', 'data/improve',
]
INCLUDE_OPTIONAL_ROOTS = ['vendor']
EXCLUDE_TOP_LEVEL = {
    '.git', '.venv', '.venv.local', '.pytest_cache', '.nextflow',
    'tools', 'results', 'work', 'dist', 'conda_packs',
}
EXCLUDE_PARTS = {'__pycache__'}
EXCLUDE_SUFFIXES = {'.pyc', '.pyo'}
EXCLUDE_PREFIXES = ('.nextflow.log',)
EXCLUDE_NAMES = {'human_v102.tar.gz.ad', 'neoantigen2-main.zip', 'migrate_export.log', 'wget_ad.log', 'tools.env.local.sh'}
EXCLUDE_REL_PREFIXES = (
    'resources/arriba_v2.5.1_official/',
)


def excluded(path: Path) -> bool:
    rel = path.relative_to(ROOT)
    parts = set(rel.parts)
    if rel.parts and rel.parts[0] in EXCLUDE_TOP_LEVEL:
        return True
    if parts & EXCLUDE_PARTS:
        return True
    if path.name in EXCLUDE_NAMES:
        return True
    if any(path.name.startswith(prefix) for prefix in EXCLUDE_PREFIXES):
        return True
    if path.suffix in EXCLUDE_SUFFIXES:
        return True
    if path.name.startswith('._'):
        return True
    if '.bak_' in path.name:
        return True
    rel_posix = rel.as_posix()
    if any(rel_posix.startswith(prefix) for prefix in EXCLUDE_REL_PREFIXES):
        return True
    return False


def iter_files() -> list[Path]:
    files: set[Path] = set()
    for root in INCLUDE_ROOTS + INCLUDE_DATA_ROOTS + INCLUDE_OPTIONAL_ROOTS:
        base = ROOT / root
        if not base.exists():
            continue
        for p in base.rglob('*'):
            if p.is_file() and not excluded(p):
                files.add(p)
    for rel in INCLUDE_FILES:
        p = ROOT / rel
        if p.is_file() and not excluded(p):
            files.add(p)
    return sorted(files, key=lambda x: x.relative_to(ROOT).as_posix())


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description='Build lightweight v04 release archive.')
    parser.add_argument('--outdir', default=str(DEFAULT_OUTDIR))
    parser.add_argument('--name', default=None)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    stamp = dt.datetime.now(dt.timezone.utc).strftime('%Y%m%d')
    name = args.name or f'neoag_event_pipeline_v043_lite_{stamp}.tar.gz'
    tar_path = outdir / name
    files = iter_files()
    root_name = tar_path.name[:-7] if tar_path.name.endswith('.tar.gz') else tar_path.stem
    with tarfile.open(tar_path, 'w:gz') as tf:
        for p in files:
            tf.add(p, arcname=str(Path(root_name) / p.relative_to(ROOT)))
    digest = sha256(tar_path)
    sha_path = tar_path.with_suffix(tar_path.suffix + '.sha256')
    sha_path.write_text(f'{digest}  {tar_path.name}\n', encoding='utf-8')
    print(tar_path)
    print(f'files={len(files)} size={tar_path.stat().st_size} sha256={digest}')


if __name__ == '__main__':
    main()
