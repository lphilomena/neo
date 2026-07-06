# NeoAg Event Pipeline v0.4.3

NeoAg Event Pipeline is a research-oriented neoantigen prioritization pipeline. It converts SNV/InDel,
fusion, splice, structural-variant, and peptide-only candidates into standardized event and peptide-HLA
tables, then layers presentation, APPM, CCF, safety, immune-escape, validation-plan, and report evidence.

**Interpretation boundary**: the pipeline produces computational triage and validation-planning outputs.
It does not make clinical diagnoses, clinical resistance calls, or validated treatment recommendations.

## What It Does

- Parses pVACtools-like SNV/fusion/splice outputs, or generates sliding-window variant peptides directly
  from VEP-annotated VCFs (with optional automatic VEP annotation).
- Scores MHC presentation evidence (NetMHCpan, MHCflurry, optional stability/immunogenicity tools).
- Builds APPM 2.0 (antigen presentation machinery), CCF 2.1 (clonality), peptide safety, and immune-escape
  evidence layers.
- Produces ranked event/peptide tables, a validation plan, and dual-audience (patient/technical) HTML reports.
- Runs via the `neoag-v03` CLI directly, or through the included Nextflow wrappers.

The `.v03.tsv` suffix in ranked outputs is a schema-compatibility label, not the software version.

## Six Input Entries

The pipeline is organized around 6 independent input entries. Each one runs from its own raw input through
a shared scoring/report tail. See **[`docs/USAGE_GUIDE.md`](docs/USAGE_GUIDE.md)** for the full command-by-command
tutorial (module chains, intermediate files, final outputs, parameters, and required environment per module).

| Entry | Input | AI skill (`.agent/skills/`) |
|---|---|---|
| SNV/InDel | Somatic VCF (+ optional pVACseq) | `neoag-vcf` |
| Fusion | EasyFuse `fusions.pass.csv` (+ optional pVACfuse) | `neoag-fusion` |
| Splice junction | VCF + RegTools junction TSV | `neoag-splice` |
| SV (WGS) | Structural-variant VCF + GTF + FASTA | `neoag-sv-wgs` |
| SV (WES) | Same, + capture BED | `neoag-sv-wes` |
| Peptide-only | Peptide + HLA CSV/TSV | `neoag-peptide-csv` |

If you're using an AI coding/agent tool that supports `.agent/skills/`, start with `pipeline-get` — it checks
the environment, lists the 6 entries above, and routes to the right one based on your input file. Each entry
skill and `docs/USAGE_GUIDE.md` describe the exact same commands; keep both in sync when either changes.

## Quick Start

```bash
python -m pip install -e '.[test]'
neoag-v03 run-demo --entry-mode snv_indel --outdir work/demo_snv --sample-id DEMO001
```

`run-demo --entry-mode {snv_indel,fusion,splice_junction,sv_wgs,sv_wes,peptide_only}` runs a full,
fixture-backed smoke test for a single entry and prints only the tools relevant to that entry — no need
to install every optional tool before trying the pipeline.

## Tests

```bash
pytest -q                       # fast unit tests + release-safe contract checks (default)
pytest -q --run-integration
pytest -q --run-benchmark
pytest -q --run-external
pytest -q --run-all             # everything; needs external tools/network/Nextflow cache
```

## Documentation

| Doc | Covers |
|---|---|
| [`docs/USAGE_GUIDE.md`](docs/USAGE_GUIDE.md) | Full tutorial: per-entry command chains + module parameter/output/environment reference |
| [`docs/INSTALL_AND_DATA.md`](docs/INSTALL_AND_DATA.md) | Base environment, Python/conda, Nextflow, external tool install table, reference data, acceptance checks |
| [`docs/CHANGELOG.md`](docs/CHANGELOG.md) | Version history (v0.4 → v0.4.3, skills refactor) |
| [`docs/RELEASE_BOUNDARY.md`](docs/RELEASE_BOUNDARY.md) | What ships in the lightweight online release vs. what stays external/local |
| [`docs/SITE_CONFIG_BOUNDARY.md`](docs/SITE_CONFIG_BOUNDARY.md) | Which config files are committed vs. kept site-local |

## License And Citation

See [`LICENSE`](LICENSE), [`NOTICE`](NOTICE), and [`CITATION.cff`](CITATION.cff).
