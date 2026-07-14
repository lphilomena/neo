#!/usr/bin/env python3
import argparse
import shutil
from pathlib import Path

parser = argparse.ArgumentParser(description='Create target-machine local manifests from templates.')
parser.add_argument('--project-root', default='.')
parser.add_argument('--outdir', default='configs/local')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--deploy-root', default='/opt/neoag')
parser.add_argument('--rewrite-deploy-root', action='store_true', help='Rewrite /opt/neoag in templates to --deploy-root')
args = parser.parse_args()
root = Path(args.project_root).resolve()
outdir = (root / args.outdir).resolve()
outdir.mkdir(parents=True, exist_ok=True)
source_dir = root / 'configs' / 'controlled_execution'
items = {
    'tools_manifest.example.yaml': 'tools_manifest.yaml',
    'reference_manifest.example.yaml': 'reference_manifest.yaml',
    'sample_manifest.example.yaml': 'sample_manifest.yaml',
}
for src_name, dst_name in items.items():
    src = source_dir / src_name
    dst = outdir / dst_name
    if not src.exists():
        raise SystemExit(f'MANIFEST_TEMPLATE_MISSING: {src}')
    if dst.exists() and not args.overwrite:
        print(f'keep_existing={dst}')
    else:
        content = src.read_text()
        if args.rewrite_deploy_root:
            content = content.replace('/opt/neoag', str(Path(args.deploy_root)))
        dst.write_text(content)
        print(f'created={dst}')
print('Edit these manifests for target-machine paths. Do not commit private local paths.')
