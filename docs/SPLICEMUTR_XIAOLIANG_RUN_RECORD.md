# SpliceMutr run record: Xiaoliang liver tumor RNA

## Scope

- RNA sample: `FP500004780_L01_203`
- Patient analysis ID: `M1ML150017383`
- Genome build: GRCh38
- Completed comparison date: 2026-07-23
- Purpose: reconstruct splice-junction peptides with SpliceMutr and compare them
  with SNAF by exact junction, peptide, HLA allele, and NetMHCpan EL-rank class.

## Runtime and assets

- SpliceMutr environment: `work/envs/neoag-splicemutr`
- Python: 3.11.0
- R: 4.3.1
- Source snapshot: `work/splice_tools/SpliceMutr-main`
- The historical source directory does not retain an independently verifiable
  upstream SpliceMutr tag or commit. Future runs must record this in the run
  manifest.
- Reference assets include a GRCh38/Gencode STAR index and TxDb. The install
  skill now verifies the SpliceMutr runtime and its BSgenome/annotation assets.

## Inputs and coordinate policy

The strict production input is:

```text
/mnt/zjl-bgi-zzb/peixunban/gl/data/chenxiaoliang_data/work/snaf/
FP500004780_L01_203_20260721/snaf_result_8_11/
complete_comparison_splicemutr/snaf_splicemutr_9_10_strict_candidates.tsv
```

Its SHA-256 is:

```text
3555d4256027b5d70c4ddea15cb418ba090325bccca3eaaa3684b910a8477716
```

SNAF and SpliceMutr source coordinates are exon splice boundaries. They are
converted to GRCh38 1-based closed intron coordinates with `start + 1` and
`end - 1`. Matching never falls back to gene, nearest locus, or peptide alone.

## Cross-validation rules

An exact SNAF-SpliceMutr match requires all of the following:

1. Same normalized junction coordinate and strand.
2. Same peptide sequence.
3. Same HLA allele.
4. Same NetMHCpan EL-rank class.

EL-rank classes are strong (`<= 0.5%`), weak (`> 0.5%` and `<= 2%`), and
non-binder (`> 2%`). `strong_cross_validated` additionally requires the
strong-binder class.

## Results

### 9/10-mer production set

| Metric | Count |
|---|---:|
| SNAF candidate rows | 73,675 |
| Unique SNAF events | 14,993 |
| SpliceMutr reconstructed rows | 37,603 |
| Events with exact structure | 8,975 |
| Events with a junction-spanning reconstructed peptide | 1,831 |
| Strict peptide-HLA rows | 7,962 |
| Strong strict rows | 2,208 |
| Events with at least one strong strict result | 1,076 |

The 7,962 strict rows cover 1,831 splice events. Downstream prediction covered
NetMHCpan, MHCflurry, PRIME, and BigMHC for all 7,962 rows. The resulting
evidence grades were 3,033 R3 and 4,929 R4; all 4,929 R4 rows had an exact
reference-proteome hard failure.

### 8/11-mer reconstruction comparison

| Length | SNAF events | SpliceMutr events | Shared events | Exact shared junction-peptides | Exact shared binders |
|---:|---:|---:|---:|---:|---:|
| 8 | 19,262 | 4,919 | 4,592 | 15,266 | 851 |
| 11 | 19,115 | 4,918 | 4,563 | 21,013 | 1,194 |

The 2,045 exact shared binders represent 1,070 junctions and 1,613 unique
peptides. For identical junction-peptide-HLA combinations, NetMHCpan EL-rank
class concordance was 100%.

## Canonical outputs

Primary cross-validation directory:

```text
/mnt/zjl-bgi-zzb/peixunban/gl/data/chenxiaoliang_data/work/snaf/
FP500004780_L01_203_20260721/snaf_result_8_11/complete_comparison_splicemutr
```

Integrated production directory:

```text
/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata4git/work/analysis_results/
chenxiaoliang_wgs_wes_easyfuse_joint_sequenza012_rna_validation/
phased_run_splice_crossvalidated_20260723
```

Important files:

- `snaf_splicemutr_9_10_strict_candidates.tsv`: strict 9/10-mer input.
- `snaf_splicemutr_8_11_strict_binders_unique.tsv`: strict 8/11-mer binders.
- `SNAF_vs_SpliceMutr_complete_report_8_9_10_11.md`: complete comparison.
- `splice_snaf/crossvalidated_splice_prepare_manifest.json`: normalized inputs,
  hashes, coordinate policy, and prepared outputs.
- `qc/splice_crossvalidated_merge_manifest.json`: downstream counts, tool
  coverage, safety, grades, and output hashes.
- `scoring/splice_crossvalidated_evaluation_summary.tsv`: compact result table.

## Interpretation limits

- SpliceMutr reconstruction was seeded from SNAF-filtered junctions. Event-level
  overlap is therefore not independent discovery validation.
- The informative independent layer is transcript/ORF/peptide reconstruction.
- Identical NetMHCpan results are reproducibility evidence, not independent
  biological confirmation, because both branches used the same model.
- RNA-only splice events do not receive an inferred CCF. Their CCF remains
  `RNA_ONLY_UNRESOLVED` unless supported by an appropriate DNA event model.
- A shared junction is not automatically a confirmed neoantigen. It still
  requires frame/ORF support, presentation, safety, and experimental review.
