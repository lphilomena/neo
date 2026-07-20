# Production Neoantigen Workflow

This document defines the complete production flow and its machine-readable handoff points. The discovery branches run in parallel; they are not one long serial command.

## 1. Required sample inputs

| Input | Main use |
| --- | --- |
| Tumor and matched-normal DNA BAM/CRAM | Somatic variants, purity, CNV, CCF, and HLA LOH |
| Somatic PASS VCF | SNV/InDel peptide generation and VEP/APPM evidence |
| Tumor RNA FASTQ/BAM | Gene/transcript expression, RNA variant support, fusion, and aberrant splice evidence |
| Patient HLA alleles | All peptide-HLA binding and presentation predictions |
| GRCh38 FASTA, GTF, and VEP cache | Consistent variant, transcript, fusion, and splice interpretation |

## 2. Parallel evidence branches

| Branch | Recommended tools | Required handoff |
| --- | --- | --- |
| HLA typing | OptiType plus SpecHLA and/or HLA-LA cross-check | Normalized HLA allele list |
| Purity/CNV/CCF | FACETS; Sequenza/ASCAT/PURPLE cross-check where available | `purity.tsv`, `cnv_segments.tsv`, CCF evidence |
| HLA LOH/APPM | LOHHLA plus HLA/CNV cross-check | `hla_loh.tsv` and APPM completeness state |
| RNA expression | Salmon gene and transcript TPM; RSEM cross-check where required | Gene expression TSV; transcript abundance retained in branch outputs |
| SNV/InDel | VEP plus pVACseq or sliding-window extraction | `pvacseq_aggregated.tsv` or standard raw tables |
| Fusion | EasyFuse/Arriba/STAR-Fusion discovery, frame review, then pVACfuse | `pvacfuse_aggregated.tsv` with junction-crossing peptides |
| Aberrant splice | RegTools plus pVACsplice; optional SNAF cross-check | `pvacsplice_aggregated.tsv` with junction support |

HLA typing must finish before peptide-HLA prediction. Purity, CNV, HLA LOH, and RNA evidence may run in parallel and are merged during ranking.

## 3. Candidate peptide contract

Every peptide source must be normalized into `raw_events.tsv` and `raw_peptides.tsv`. Fusion and splice candidates are not complete merely because a caller found an event.

A production fusion peptide must:

- cross the fusion breakpoint or contain fusion-derived novel amino acids;
- retain caller, frame, junction-read, and transcript-expression evidence;
- enter NetMHCpan and MHCflurry using the patient's HLA alleles;
- enter PRIME/BigMHC or another configured immunogenicity layer;
- pass normal-proteome, normal-expression, and normal-ligand safety review.

A production splice peptide must similarly cross the abnormal exon junction, retain junction support, and pass the same presentation and safety layers.

## 4. Unified presentation and ranking

Use `conf/run.production_multisource.example.toml` after the three peptide branches finish:

```bash
source conf/tools.env.sh
neoag-v03 run-full \
  --config conf/run.production_multisource.example.toml \
  --outdir results/SAMPLE001_multisource
```

The production example requires all of these peptide sources:

- `pVACseq`
- `pVACfuse`
- `pVACsplice`

It also requires NetMHCpan and MHCflurry outputs. Missing peptide sources do not stop the run: the output is marked `LOW_CONFIDENCE` with `missing_peptide_sources`. Missing required predictor outputs still fail because presentation ranking would be invalid. NetMHCstabpan remains optional and skipped by default.

The merged order is:

Source coverage is persisted in `upstream/tools/peptide_source_coverage.tsv`, including expected, detected, and missing sources plus the final completeness status.

1. Normalize SNV/InDel, fusion, and splice candidates.
2. Merge all peptide-HLA pairs.
3. Run NetMHCpan and MHCflurry on the merged set.
4. Add PRIME/BigMHC and configured immunogenicity evidence.
5. Add RNA expression and event-level RNA support.
6. Add purity, CNV, CCF, HLA LOH, and APPM evidence.
7. Apply normal-proteome/ligandome safety and immune-escape modifiers.
8. Produce unified event, peptide, validation-plan, and report outputs.

## 5. Resume without repeating prediction

To reuse completed NetMHCpan and MHCflurry results:

1. Remove `netmhcpan` and `mhcflurry` from `[tools].enabled`.
2. Set `[inputs].netmhcpan` and `[inputs].mhcflurry` to the existing output files.
3. Keep `required_presentation_predictors = ["netmhcpan", "mhcflurry"]`.
4. Rerun `run-full`; cached outputs are loaded and upstream prediction is not repeated.

The same approach applies to cached `vep_appm`, expression, purity, CNV, and HLA LOH sidecars.

## 6. Production acceptance checks

A complete run must contain:

| Check | Expected result |
| --- | --- |
| HLA typing | Consensus alleles and conflict flags recorded |
| Source coverage | Detected and missing sources recorded; missing fusion/splice marks `LOW_CONFIDENCE` |
| Presentation | NetMHCpan and MHCflurry outputs present |
| Fusion quality | Breakpoint/frame/junction evidence retained |
| Splice quality | Abnormal junction and RNA read support retained |
| DNA context | Purity, CNV, and CCF assessed or explicitly unassessed |
| Immune escape | HLA LOH/APPM assessed or explicitly unassessed |
| RNA context | Gene expression plus transcript/junction evidence available |
| Safety | Normal expression/proteome/ligand evidence assessed |
| Final outputs | Ranked events, ranked peptides, validation plan, and reports present |

Key outputs are under `scoring/`, `presentation/`, `evidence/`, and `reports/`. RNA-only fusion/splice events should use RNA abundance and junction support rather than an invented DNA CCF.
