#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

parser = argparse.ArgumentParser(description='Write deployment_report.md from deployment outputs.')
parser.add_argument('--project-root', default='.')
parser.add_argument('--workdir', default='work/remote_deploy')
args = parser.parse_args()
root = Path(args.project_root).resolve()
workdir = (root / args.workdir).resolve()
workdir.mkdir(parents=True, exist_ok=True)
report = workdir / 'deployment_report.md'
doctor_status = workdir / 'doctor' / 'doctor_status.json'
status = 'UNKNOWN'
if doctor_status.exists():
    try:
        data = json.loads(doctor_status.read_text())
        status = data.get('status') or data.get('overall_status') or 'UNKNOWN'
    except Exception:
        status = 'UNKNOWN'
lines = [
    '# Deployment report',
    '',
    f'Project root: `{root}`',
    f'Deployment workdir: `{workdir}`',
    f'Doctor status: `{status}`',
    '',
    '## Outputs',
    '',
]
for rel in [
    'preflight_report.md',
    'smoke_test_report.md',
    'doctor/doctor_summary.md',
    'doctor/blocking_issues.tsv',
    'doctor/recommended_fixes.md',
]:
    path = workdir / rel
    if path.exists():
        lines.append(f'- `{path}`')
lines += [
    '',
    '## Safety boundary',
    '',
    'This deployment report does not certify clinical use. Missing tools or references are deployment gaps, not biological negative evidence.',
]
report.write_text('\n'.join(lines) + '\n')
print(f'deployment_report={report}')
