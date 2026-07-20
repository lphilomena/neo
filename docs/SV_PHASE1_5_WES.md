# WES SV Phase 1.5

Exome-capture-limited SV neoantigen triage for paired tumor-normal **WES** samples.

Phase 1 (`sv_wgs_phase1`) targets WGS with higher split/pair support thresholds. Phase 1.5 relaxes DNA support requirements and elevates **RNA junction evidence** to define `WES_Tier1` / `WES_Tier2` tiers.

## WES tier rules

| Tier | Criteria |
|------|----------|
| `WES_Tier1` | RNA junction reads ≥ 3, **or** base WGS Tier1 confidence |
| `WES_Tier2` | RNA junction reads ≥ 1, **or** base WGS Tier2 confidence |
| `WES_Tier3` | Rejected unless `--tier1-only` is off and tier is WES_Tier2 |

Evidence scope: `EXOME_CAPTURE_LIMITED` — intergenic and off-target SVs may be missed.

## CLI

```bash
source conf/tools.env.sh

# Adapter only
neoag sv-build-raw-wes \
  --sample-id SVMINI \
  --profile sv_wes_phase1_5 \
  --sv-vcf data/fixtures_sv/mini_sv.vcf \
  --callers GRIDSS2 \
  --reference-fasta data/fixtures_sv/mini_ref.fa \
  --gencode-gtf data/fixtures_sv/mini.gtf \
  --hla data/fixtures_sv/hla.txt \
  --expression data/fixtures_sv/expression.tsv \
  --rna-junctions data/fixtures_sv/rna_junctions.tsv \
  --outdir results/SVMINI_sv_wes_adapter

# End-to-end (adapter + NetMHCpan/MHCflurry + score)
neoag sv-run-full-wes \
  --sample-id SVMINI \
  --profile sv_wes_phase1_5 \
  --sv-vcf data/fixtures_sv/mini_sv.vcf \
  --callers GRIDSS2 \
  --reference-fasta data/fixtures_sv/mini_ref.fa \
  --gencode-gtf data/fixtures_sv/mini.gtf \
  --hla data/fixtures_sv/hla.txt \
  --rna-junctions data/fixtures_sv/rna_junctions.tsv \
  --binding-stub \
  --outdir results/SVMINI_sv_wes_e2e
```

## Nextflow

```bash
nextflow run workflows/sv_phase1_5_wes.nf -c conf/sv_wes_demo.config
```

Uses `params.wes_mode = true` so `SV_BUILD_RAW` invokes `sv-build-raw-wes`, then runs `NEOAG_SV_SCORE`.

## Outputs

Same layout as WGS Phase 1, with WES-specific provenance:

- `parsed/raw_events.tsv`, `parsed/raw_peptides.tsv`
- `sv/sv_events.full.tsv` (tiers labelled `WES_Tier1/2/3`)
- `provenance.sv_wes_phase1_5.json`

Downstream scoring uses profile `sv_wes_phase1_5` (see `profiles/sv_wes_phase1_5.toml`).

## Limitations

- Exome capture limits breakpoint discovery; provide RNA junction TSV when available.
- Same heuristic protein reconstruction caveats as WGS Phase 1 apply.
- Not a clinical-grade SV caller; validate candidates experimentally.

See also: [SV_PHASE1_WGS.md](SV_PHASE1_WGS.md)
