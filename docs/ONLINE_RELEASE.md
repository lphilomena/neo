# Online Release Guide

This document defines the public online release boundary for NeoAg Event Pipeline v0.4.3. The shape follows the same practical principles used by OpenRare-style repositories: keep the repository installable, auditable, and runnable with lightweight fixtures, while leaving heavyweight tools and private deployment state outside the source release.

## Release Boundary

Included in the online package:

- Python source under `src/`
- CLI wrappers under `bin/`
- Nextflow workflows and modules
- portable configs and profiles
- setup scripts and release scripts
- tests and lightweight fixtures
- Markdown documentation, license, notice, citation, release notes, manifests, and conda environment manifests including `conda/env.neoag-ascat-v3.yml`

Excluded from the online package:

- `.git`, virtual environments, caches, and bytecode
- `tools/`, `results/`, `work/`, `dist/`, `conda_packs/`
- `.nextflow*` metadata and logs
- large references under `data/ref`, `data/vep`, `data/external`, `data/examples`
- local/private deployment files such as `conf/tools.env.local.sh`, `conf/site.config`, `conf/private/*`, `conf/*.private.toml`, `conf/*.local.toml`
- patient-specific scripts, sample identifiers, and site-local absolute paths
- editable office artifacts such as `.docx` and `.pptx`

## Build

```bash
python scripts/package_online_release.py --outdir work/releases
```

The script writes:

- `neoag_event_pipeline_v043_online_<date>.tar.gz`
- `neoag_event_pipeline_v043_online_<date>.tar.gz.sha256` external checksum file
- `neoag_event_pipeline_v043_online_<date>.manifest.json`

## Smoke Test

```bash
tmpdir=$(mktemp -d)
tar -xzf work/releases/neoag_event_pipeline_v043_online_<date>.tar.gz -C "$tmpdir"
cd "$tmpdir"/neoag_event_pipeline_v043_online_<date>
python -m pip install -e '.[test]'
pytest -q
neoag-v03 run-demo --outdir work/demo_v043 --sample-id DEMO001
```

For Nextflow fixture smoke:

The online package includes the lightweight `bin/nextflow` launcher but does not bundle the Nextflow dependency cache. Online mode may download runtime dependencies into `NXF_HOME` on first use. Offline mode must pre-stage a populated `NXF_HOME` cache, Java, and any conda/container assets before launch. Reuse a shared writable cache when available, for example `export NXF_HOME=/path/to/nextflow_cache`.

```bash
bin/neoag-nextflow run workflows/main.nf   -w /tmp/neoag_nf_work   --pvac_files data/fixtures/pvacseq_aggregated.tsv   --outdir /tmp/neoag_nf_demo   --sample_id NF_DEMO
```

## Deployment Notes

For complete environment, tool, and reference-data setup, read `docs/INSTALL_AND_DATA.md`. Set `NEOAG_TOOLS_ROOT` to a shared artifact installation before real-data runs. The fixture demo and default test suite do not require licensed or heavyweight tools.

```bash
export NEOAG_TOOLS_ROOT=/path/to/neoag_artifacts
source conf/tools.env.sh
neoag-v03 check-tools
```

## Verification Performed For This Release

The online release should pass, at minimum:

- `pytest -q`
- `neoag-v03 run-demo --outdir work/demo_v043 --sample-id DEMO001`
- `bin/neoag-nextflow run workflows/main.nf` with lightweight fixture input
- `scripts/check_release_boundary.sh`

## Interpretation Boundary

NeoAg Event Pipeline is a research triage and validation-planning tool. It does not make clinical treatment recommendations.
