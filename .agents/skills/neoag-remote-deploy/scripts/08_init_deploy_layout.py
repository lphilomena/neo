#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

parser = argparse.ArgumentParser(description="Initialize a four-way NeoAg deployment layout.")
parser.add_argument("--deploy-root", default="/opt/neoag")
parser.add_argument("--project-root", default=".")
parser.add_argument("--outdir", default="work/remote_deploy")
parser.add_argument("--create", action="store_true", help="Create directories. Without this flag, only write a plan.")
args = parser.parse_args()

deploy_root = Path(args.deploy_root)
project_root = Path(args.project_root).resolve()
outdir = (project_root / args.outdir).resolve()
outdir.mkdir(parents=True, exist_ok=True)

dirs = [
    "code",
    "conf/profiles",
    "envs/conda",
    "envs/containers/docker",
    "envs/containers/apptainer",
    "refs/GRCh38/fasta",
    "refs/GRCh38/gencode",
    "refs/GRCh38/vep",
    "refs/GRCh38/indices/bwa",
    "refs/GRCh38/indices/star",
    "refs/hla/optitype",
    "refs/hla/hla-la/graph",
    "refs/hla/spechla/reference",
    "refs/hla/lohhla/hla_reference",
    "refs/hla/hla_allele_metadata",
    "refs/mhc/netmhcpan/4.2",
    "refs/mhc/netmhcstabpan",
    "refs/mhc/mhcflurry/models",
    "refs/cnv/facets",
    "refs/cnv/ascat",
    "refs/cnv/purple",
    "refs/cnv/sequenza",
    "refs/fusion/arriba",
    "refs/fusion/star_fusion",
    "refs/fusion/easyfuse",
    "refs/fusion/normal_readthrough",
    "refs/safety/normal_proteome",
    "refs/safety/normal_ligandome",
    "refs/safety/normal_expression",
    "refs/safety/normal_junctions",
    "refs/validated/assay_results",
    "refs/_staging",
    "runs",
]

rows = []
for rel in dirs:
    path = deploy_root / rel
    status = "PLANNED"
    if args.create:
        try:
            path.mkdir(parents=True, exist_ok=True)
            status = "CREATED_OR_EXISTS"
        except Exception as exc:
            status = f"FAILED: {exc}"
    rows.append({"relative_path": rel, "path": str(path), "status": status})

(outdir / "deployment_layout.json").write_text(json.dumps({"deploy_root": str(deploy_root), "directories": rows}, indent=2) + "\n")
(outdir / "deployment_layout.md").write_text(
    "# NeoAg deployment layout\n\n"
    f"Deploy root: `{deploy_root}`\n\n"
    + "\n".join(f"- `{r['path']}`: {r['status']}" for r in rows)
    + "\n\nBy default refs should be read-only during pipeline execution. Use refs/_staging for new derived assets.\n"
)
print(f"deployment_layout={outdir / 'deployment_layout.md'}")
