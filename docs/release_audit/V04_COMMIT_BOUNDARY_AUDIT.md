# V04 Commit Boundary Audit

## Branch Status

- Requested branch: `codex/v04-release-boundary`
- Actual branch: `main`
- Blocker: `.git/` is owned by `root`, so user `na` cannot create branches, update `.git/info/exclude`, or stage/commit files in this checkout.

## Recommended Commit 1: v04 evidence/safety/escape + release boundary

Stage these files only after Git metadata ownership is fixed:

- `tests/test_appm_escape_consistency.py`

- `tests/test_v041_appm_ccf_escape.py`

- `src/neoag/ccf_v2.py`

- `src/neoag/appm_v2.py`

- `docs/V041_APPM_CCF_IMMUNE_ESCAPE.md`

- `CHANGELOG_V041_APPM_CCF_IMMUNE_ESCAPE.md`

- `CHANGELOG_V04_EVIDENCE_SAFETY_ESCAPE.md`
- `RELEASE_REFRESH_20260615.md`
- `RELEASE_V04_20260616.md`
- `docs/RELEASE_BOUNDARY_V04.md`
- `docs/V04_EVIDENCE_SAFETY_ESCAPE.md`
- `pyproject.toml`
- `profiles/default.toml`
- `profiles/sv_wes_phase1_5.toml`
- `conf/tools.env.sh`
- `conf/easyfuse.nextflow.config`
- `conf/run.snv_wes.example.toml`
- `conf/run.snv_wes_fixture.toml`
- `conf/snv_wes_demo.config`
- `bin/neoag-nextflow`
- `scripts/write_release_manifest.py`
- `scripts/package_v04_release.py`
- `docs/release_audit/V04_COMMIT_BOUNDARY_AUDIT.md`
- `scripts/stage_v04_release_boundary.sh`
- `src/neoag/adapters/bigmhc_im.py`
- `src/neoag/adapters/deepimmuno.py`
- `src/neoag/adapters/easyfuse_adapter.py`
- `src/neoag/adapters/event_catalog.py`
- `src/neoag/adapters/facets.py`
- `src/neoag/adapters/iedb_immunogenicity.py`
- `src/neoag/adapters/lohhla.py`
- `src/neoag/adapters/mhcflurry.py`
- `src/neoag/adapters/netmhcpan.py`
- `src/neoag/adapters/netmhcstabpan.py`
- `src/neoag/adapters/peptide_input.py`
- `src/neoag/adapters/prime.py`
- `src/neoag/appm_lite.py`
- `src/neoag/ccf_lite.py`
- `src/neoag/cli.py`
- `src/neoag/config.py`
- `src/neoag/evidence_layer.py`
- `src/neoag/evidence_provenance.py`
- `src/neoag/immune_escape.py`
- `src/neoag/immune_escape_resistance.py`
- `src/neoag/immunogenicity_composite.py`
- `src/neoag/input_router.py`
- `src/neoag/model_layers.py`
- `src/neoag/peptide_safety_gate.py`
- `src/neoag/pipeline.py`
- `src/neoag/preflight.py`
- `src/neoag/presentation.py`
- `src/neoag/schemas.py`
- `src/neoag/scoring.py`
- `src/neoag/snv_call/__init__.py`
- `src/neoag/snv_call/mutect2.py`
- `src/neoag/snv_call/pipeline.py`
- `src/neoag/sv/phase1.py`
- `src/neoag/sv/schemas_sv.py`
- `src/neoag/sv/score_pipeline.py`
- `src/neoag/sv/wes_adapter.py`
- `src/neoag/sv/wes_capture.py`
- `src/neoag/sv/wes_filter.py`
- `src/neoag/tools/registry.py`
- `src/neoag/tools/runner.py`
- `src/neoag/vep/__init__.py`
- `src/neoag/vep/annotate.py`
- `src/neoag/vep/extract_peptides.py`
- `data/fixtures/easyfuse_fusions.pass.tsv`
- `data/fixtures/hla_loh_lost.tsv`
- `data/fixtures/vep_immune_escape.tsv`
- `data/fixtures_snv/hla.txt`
- `data/fixtures_snv/mini_ref.dict`
- `data/fixtures_snv/mini_ref.fa`
- `data/fixtures_snv/mini_ref.fa.fai`
- `data/fixtures_snv/mini_somatic.vcf`
- `data/fixtures_snv/wes_capture.bed`
- `tests/conftest.py`
- `tests/test_easyfuse_adapter.py`
- `tests/test_evidence_provenance.py`
- `tests/test_immune_escape_resistance.py`
- `tests/test_intermediates.py`
- `tests/test_peptide_input.py`
- `tests/test_preflight.py`
- `tests/test_snv_phase1_wes.py`
- `tests/test_tools.py`
- `tests/test_pipeline.py`
- `tests/test_v04_evidence_safety_escape.py`
- `modules/gatk_filter_mutect_calls/main.nf`
- `modules/gatk_mutect2/main.nf`
- `modules/snv_write_run_config/main.nf`
- `workflows/snv_phase1_wes.nf`
- `workflows/snv_phase1_wes_fixture.nf`

## Explicitly Excluded From Commit Boundary

- README.md (tracked deletion; do not include until intentionally replacing docs)
- docs/TOOLS_SETUP.md (tracked deletion; do not include until docs migration is reviewed)
- scripts/install_lohhla_facets_nextflow.sh (tracked deletion; do not include without tool-install migration decision)
- scripts/install_netmhcpan.sh (tracked deletion; do not include without tool-install migration decision)
- scripts/package_lite_release.sh (tracked deletion; replaced by package_v04_release.py but deletion should be reviewed)
- results/** (local run outputs)
- work/** (local run outputs and v04 archive artifacts)
- tools/** (external installed tools)
- data/ref/** and data/vep/** (large external references/caches)
- dist/** (root-owned historical release archives)
- conda_packs/** (large local environment packs)
- human_v102.tar.gz.ad, neoantigen2-main.zip, migrate_export.log, wget_ad.log (local downloads/logs)

## Verification Already Completed

- Cross-module APPM/immune_escape/scoring consistency tests added for B2M, JAK/IFNG, CIITA/MHC-II, and no-input scenarios.

- APPM2 native sidecars and peptide-specific modifiers added; scoring consumes `appm_peptide_modifiers.tsv`.

- v0.4.1 APPM2/CCF2/Immune Escape2 patch merged; fixture pipeline run completed.

- WES SV priority caps now read from `[wes_confidence_caps]` profile settings; custom profile override is covered by regression test.

- `.venv/bin/python -m pytest -q` -> `131 passed`
- `bin/neoag run-demo --sample-id V04DEMO --profile default --outdir work/release_demo_v04_basic` -> completed
- fixture `run` with safety/immune escape inputs -> completed
- v04 lightweight archive generated in `work/releases/` but intentionally excluded from Git commit.

## Permission Fix Needed Before Actual Branch/Commit

```bash
sudo chown -R na:na .git
git switch -c codex/v04-release-boundary
scripts/stage_v04_release_boundary.sh
git commit -m "Add v04 evidence safety escape release boundary"
```
