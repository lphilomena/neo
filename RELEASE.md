# NeoAg Event Pipeline v0.4.3 Online Release

- Release name: `v043_online_20260629`
- Version: `v0.4.3`
- Release date: 2026-06-29
- Package: `neoag_event_pipeline_v043_online_20260629.tar.gz`
- External checksum: `neoag_event_pipeline_v043_online_20260629.tar.gz.sha256`
- External package manifest: `neoag_event_pipeline_v043_online_20260629.manifest.json`
- Boundary: online source/config/test/docs/fixture release; heavyweight tools, references, results, work directories, virtual environments, patient-specific scripts, site-local absolute paths, and Nextflow caches are external artifacts.

## What Is Included

This release contains the v0.4.3 neoantigen prioritization pipeline with:

- schema-compatible event and peptide parsing/scoring outputs. The schema suffix is a stable table-schema label, not the software version.
- APPM 2.0 gene/module/peptide evidence and input-status tracking.
- CCF 2.0 clonality estimates from purity, copy-number, and VAF context.
- Peptide safety gates for normal expression, normal ligandome, normal junction, matched-normal, and reference proteome context.
- Immune escape evidence from HLA LOH, APPM, CCF, B2M/JAK/APM context.
- WGS/WES SV support, including WES capture-aware Phase 1.5 evidence caps.
- Dual-audience HTML reports: patient communication report and research/technical report with dashboard/provenance.
- CLI demo, validation plan, tests, profiles, workflows, and lightweight fixtures.
- Installation and data setup guide: `docs/INSTALL_AND_DATA.md`.

## Package Boundary

Included:

- source code under `src/`
- CLI wrappers under `bin/`
- portable configs, profiles, modules, workflows, scripts
- Markdown documentation, license, notice, citation, and release notes
- lightweight fixtures under `data/fixtures`, `data/fixtures_snv`, `data/fixtures_sv`, and `data/improve`
- tests and release-safe QA checks

Excluded:

- `.git`, `.venv`, `.venv.local`, `.pytest_cache`
- `.nextflow*`, Nextflow dependency caches, and local work directories
- `tools/`, `results/`, `work/`, `dist/`, `conda_packs/`
- large references such as `data/ref`, `data/vep`, `data/external`, and `data/examples`
- local/private deployment files such as `conf/tools.env.local.sh`, `conf/site.config`, `conf/private/*`, `conf/*.private.toml`, and `conf/*.local.toml`
- patient-specific scripts, sample identifiers, and site-local absolute paths
- editable office artifacts such as `.docx` and `.pptx`
- licensed external tool binaries

## Verification Summary

The online release was verified on the source tree and the unpacked online package.

| Check | Scope | Status |
| --- | --- | --- |
| `pytest -q` | default release-safe tests | PASS |
| `neoag run-demo --outdir ...` | CLI fixture demo | PASS |
| unpacked package `pytest -q` | package smoke | PASS |
| unpacked package demo | package smoke | PASS |
| Nextflow fixture with pre-populated `NXF_HOME` | workflow smoke | PASS |
| package forbidden-path scan | release boundary | PASS |

## v043_online_20260629 RNA / Junction / HLA LOH Additions

This refresh completes three production-facing evidence gaps:

- RNA allele support: `build-evidence-layer --rna-vaf` now materializes `rna_alt_reads`, `rna_ref_reads`, `rna_depth`, `rna_vaf`, `rna_vaf_source`, and `rna_support_status` in `parsed/rna_junction_evidence.tsv`.
- Fusion/splice targeted RNA validation: fusion evidence and splice junction evidence now produce `targeted_validation_status`, `targeted_validation_source`, and `targeted_validation_method` alongside junction reads.
- HLA LOH cross-validation: `neoag crosscheck-hla-loh` compares normalized LOHHLA and SpecHLA calls, writes a detailed cross-check table, and can emit a downstream-compatible consensus `hla_loh.tsv`.

## Default Pytest vs Run-All

Default pytest is intentionally lightweight:

```bash
pytest -q
```

It runs fast unit tests and release-safe contract checks, while skipping integration, benchmark, external-tool, and long workflow tests.

`--run-all` is a maintainer/release verification mode:

```bash
pytest -q --run-all
```

It opts into broader checks and may require external tools, references, network access, writable Nextflow cache, and substantially more runtime. Use `--run-integration`, `--run-benchmark`, and `--run-external` when you only want one class of expanded checks.

## Nextflow Online/Offline Boundary

The online package includes workflow source and the lightweight `bin/nextflow` launcher, but it does not bundle the Nextflow runtime dependency cache.

Online mode:

```bash
export NXF_HOME=/path/to/writable/nextflow_cache
bin/neoag-nextflow -version
```

The first run may download Nextflow runtime dependencies into `NXF_HOME`.

Offline mode:

- pre-stage Java 11+
- pre-stage a populated writable `NXF_HOME` cache
- pre-stage conda/container assets required by the selected workflow
- pre-stage all external tools and reference data
- launch via `bin/neoag-nextflow`, not a root-owned `.nextflow` directory

## Tool And Data Setup

External tools and large references are intentionally not bundled. Read the full deployment guide before real-data runs:

- `docs/INSTALL_AND_DATA.md`
- `docs/TOOLS_SETUP.md`
- `docs/SITE_CONFIG_BOUNDARY.md`

Typical setup:

```bash
export NEOAG_TOOLS_ROOT=/path/to/neoag_artifacts
source conf/tools.env.sh
neoag check-tools
```

The fixture demo and default tests do not require licensed or heavyweight tools.

## Quick Smoke

```bash
python -m pip install -e '.[test]'
pytest -q
neoag run-demo --outdir work/demo_v043 --sample-id DEMO001
```

## Interpretation Boundary

This is a research triage and validation-planning pipeline. It does not make clinical treatment recommendations. Candidate rankings require assay validation, disease context, HLA typing, tumor purity, expression/protein support, and clinical governance.
