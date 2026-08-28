# Sarcoma RNA-supported v2 provisional profile

`sarcoma_rna_supported_v2_provisional` is a research-only scoring profile for
sarcoma samples with matched tumor RNA evidence. It is not a disease-trained
model and is not calibrated for clinical decision-making.

## Changes from v1

- Keeps normal critical-tissue gene expression and HSPC expression as
  `CAUTION` evidence instead of using either signal alone as a hard rejection.
- Distinguishes HPA single-cell HSPC `nCPM` from bulk-expression `TPM` and uses
  separate provisional thresholds.
- Makes safety multipliers profile-configurable. The v2 caution multiplier is
  `0.85`; legacy profiles retain their previous behavior.
- Requires at least 1 TPM event expression, 3 RNA alternate reads, 0.02 RNA
  VAF, and 0.5 allele expression (`gene TPM * RNA VAF`) for non-junction RNA
  support.
- Reduces the fusion source multiplier from 1.55 to 1.20.
- Applies a `C_CAUTION` cap to RNA-only fusion candidates while clonality is
  unresolved or the normal-junction reference has not been assessed.
- Retains mutant-versus-wild-type specificity and phasing priority caps.

## Usage

Pass the profile name to commands that accept `--profile`:

```bash
neoag score \
  --profile sarcoma_rna_supported_v2_provisional \
  ...
```

Use the same profile for peptide safety, scoring, validation-plan generation,
and report generation within a synchronized run.

## Interpretation

`A` or `B` should not be forced by relaxing missing-evidence controls. A
high-scoring RNA-only fusion remains `C_CAUTION` until DNA clonality and normal
junction evidence are resolved. A high-scoring SNV/InDel with weak or marginal
mutant-versus-wild-type specificity also remains `C_CAUTION` and requires a
paired mutant/wild-type validation design.

The HPA nCPM thresholds and all score weights remain provisional. Recalibrate
them against an experimentally validated sarcoma cohort before production use.
