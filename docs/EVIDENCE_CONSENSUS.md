# Evidence-consensus parallel ranking

## Purpose

The evidence-consensus layer is a parallel review view. It does not modify or
replace the existing weighted neoantigen ranking, `efficacy_score`, or
`final_priority`. The full pipeline first builds
`comprehensive_peptide_evidence.tsv`; consensus code reads that table rather
than calling or changing `score()`.

## CLI

The protected first-phase entry point is:

```bash
neoag evidence-rank \
  --comprehensive-evidence results/sample/scoring/comprehensive_peptide_evidence.tsv \
  --weighted-baseline results/sample/scoring/ranked_peptides.tsv \
  --rules configs/ranking/sarcoma_evidence_consensus_v1.toml \
  --provenance results/sample/provenance.json \
  --outdir results/sample/scoring/evidence_consensus \
  --mode parallel \
  --track all \
  --emit-event-ranking \
  --compare-weighted \
  --deterministic
```

Defaults are `--mode parallel`, `--track all`, event ranking, weighted
comparison, and deterministic ordering. Optional tracks are `missense`,
`frameshift`, `fusion`, `splice`, `dna_sv`, and `manual_review`.

The first-phase parser does not expose `--replace-primary-ranking`; passing it
is an error. The older `evidence-consensus-rank` command remains available only
as a compatibility entry point.

## Data flow

```text
raw evidence tables -> score() -> ranked_peptides.tsv
                   -> comprehensive_peptide_evidence.tsv
                   -> evidence_states.tsv
                   -> hard-failure and priority-cap constraints
                   -> evidence grade R1-R4
                   -> within-track Pareto fronts
                   -> deterministic tie-break
                   -> peptide/event consensus rankings and comparison
```

Candidates are first assigned R1-R4 and then routed to `MISSENSE`,
`FRAMESHIFT`, `FUSION`, `SPLICE`, `DNA_SV`, or `MANUAL_REVIEW`. Pareto fronts
are calculated only within the same evidence grade and track, so a conventional
SNV never numerically dominates an RNA-only fusion. Manual-review routing does
not change `biological_event_track`, which remains the basis for evidence-state
interpretation.

Track-specific dimensions are:

- Missense: event authenticity, RNA support, presentation consensus, mutant
  specificity, safety, HLA/APM, clonality, and evidence completeness.
- Frameshift: event authenticity, RNA support, novel-tail evidence,
  presentation consensus, safety, HLA/APM, and clonality.
- Fusion/splice: junction authenticity, junction reads, frame evidence,
  presentation consensus, normal-junction safety, HLA/APM, and completeness.
- DNA SV: event authenticity, RNA support, novel-sequence evidence,
  presentation consensus, safety, HLA/APM, clonality, and completeness.
- Manual review: the full base evidence vector; driver status cannot upgrade
  the R grade.

Identical vectors are collapsed before non-dominated sorting, so complexity
follows the number of unique state combinations rather than the number of
peptide-HLA rows.

## Deterministic tie-break

Within one R grade, track, and Pareto front, candidates are ordered by this
fixed sequence:

1. safety evidence completeness;
2. event authenticity;
3. measured RNA support;
4. presentation consensus;
5. mutant specificity;
6. restricting-HLA/APM integrity;
7. CCF confidence/clonality;
8. overall evidence completeness;
9. NetMHCpan mutant EL rank, lower first and missing last;
10. MHCflurry presentation score, higher first and missing last;
11. `peptide_id` lexical order.

`evidence_rank_key` serializes the grade, track, Pareto front, state labels,
normalized predictor values, and peptide identifier. The key is intended for
audit and comparison; the structured columns remain authoritative. Legacy
weighted rank and input row order are not tie-break criteria, so row-order
changes cannot alter the consensus ordering.

CCF confidence is normalized independently as `CCF_HIGH_CONFIDENCE`,
`CCF_MEDIUM_CONFIDENCE`, `CCF_LOW_CONFIDENCE`, or `CCF_UNASSESSED`; the
tie-break uses `ccf_confidence_grade`, while the Pareto vector may still use the
broader clonality state. Missing predictor values and unassessed CCF confidence
sort after assessed values.

Implementation is split across:

- `src/neoag/evidence_states.py`: provisional threshold-to-state derivation.
- `src/neoag/pareto.py`: deterministic Pareto fronts over unique vectors.
- `src/neoag/evidence_consensus.py`: orchestration, hard failures, caps, R1-R4,
  outputs, comparison, and provenance.
- `configs/ranking/sarcoma_evidence_consensus_v1.toml`: independent rules.

It uses the current schema fields for seven layers:

| Layer | Preferred numeric field | Availability/status examples |
| --- | --- | --- |
| Presentation | `l3_hla_presentation_score` | `presentation_gate_status`, presentation grade |
| Binding | `l3_hla_binding_score` | presentation gate/grade |
| Expression | `l3_expression_score` | `expression_evidence_status` |
| Clonality | `l3_clonality_score`, then `ccf_estimate` | CCF status/resolution/confidence |
| Tumor specificity | `l3_tumor_specificity_score` | cross-platform status |
| Safety | `l3_normal_tissue_safety_score` | safety completeness/status |
| APPM | `l3_apm_integrity_score` | APPM completeness/integrity |

## Evidence states

- `SUPPORTED`: evidence was assessed. A `FAIL` or rejection is assessed
  negative evidence, not a data conflict.
- `PARTIAL`: incomplete, caution, review, or low-confidence evidence.
- `MISSING`: unavailable or unresolved evidence. It is not interpreted as a
  negative result.
- `CONFLICT`: discordant or internally inconsistent evidence requiring review.

## Output fields

- `evidence_consensus_score`: weighted score across all configured layers;
  missing layers contribute no evidence and remain listed explicitly.
- `evidence_consensus_quality_score`: score normalized over assessed layers.
- `evidence_completeness_score`: fraction of configured evidence weight that
  was assessed.
- `evidence_consensus_status`: `COMPLETE`, `PARTIAL_EVIDENCE`, `LOW_EVIDENCE`,
  `CONFLICT_REVIEW`, or `UNASSESSED`.
- `evidence_supporting_layers`, `evidence_assessed_layers`,
  `evidence_missing_layers`, `evidence_conflict_layers`.
- `evidence_layer_states`: auditable per-layer state summary.
- `evidence_rank`: deterministic rank in the parallel output.
- `hard_failure`, `hard_failure_codes`, `hard_failure_reasons`,
  `consensus_priority_cap`: explicit
  constraints applied only to the consensus branch.
- `evidence_grade_uncapped`, `evidence_grade`: R1-R4 before and after hard
  failure/priority-cap constraints.
- `biological_event_track`: biological event type used for state derivation.
- `evidence_track`, `pareto_dimensions`, `pareto_front`, `track_rank`:
  grade-and-track-specific Pareto placement and its auditable dimensions.
- `consensus_trace`: row-level explanation of missing/conflicting evidence and
  any grade constraint.
- `*_reason_code`: stable machine-readable reason code for every normalized
  evidence domain; `*_reason` retains its human-readable explanation.
- `evidence_reason_codes`: compact `domain:code` trace across all domains.
  Priority-cap codes remain in `evidence_grade_cap_reasons`.

The generated bundle contains:

- `ranked_peptides.tsv`: unchanged weighted production ranking and stable
  compatibility interface.
- `ranked_peptides.weighted_baseline.tsv`: explicit hard-link/copy alias of the
  unchanged weighted ranking.
- `comprehensive_peptide_evidence.tsv`: authoritative merged evidence table.
- `all_tool_results.tsv`: canonical user-facing evidence table with a stable
  schema version, record type, and deterministic record identifier.
- `all_tool_results.manifest.json`: input/output checksums, source manifest,
  dimensions, required fields, and conflict/missing-layer QC counts.
- `comprehensive_evidence_manifest.json`: merged source paths, checksums, row
  counts, precedence version, output checksum, and conflict summary.
- `evidence_states.tsv`: normalized state/value/source for every evidence layer.
- `ranked_peptides.evidence_consensus.tsv`: independent peptide ranking.
- `ranked_events.evidence_consensus.tsv`: event aggregation based on each
  event's best consensus-ranked peptide.
- `evidence_consensus_summary.tsv`: compact counts by grade, track, hard failure,
  manual review, event, and conflict status.
- `ranking_compare_weighted_vs_consensus.tsv`: old/new rank, rank shift,
  constraints, and deterministic difference reason for each peptide.
- `ranking_compare_weighted_vs_consensus.md`: human-readable comparison
  summary and the largest absolute rank shifts.
- `evidence_consensus_run.json`: input/output checksums, rules identity,
  algorithm version, counts, and output manifest.
- `evidence_conflicts.tsv`: discordant evidence layers requiring manual
  reconciliation.

`weighted_vs_consensus_comparison.tsv` remains as a compatibility alias for
older callers.

For historical result directories lacking the standard `ranked_peptides.tsv`
name, the builder recognizes `ranked_peptides.cancer_annotated.tsv` and then
`ranked_peptides.v03.tsv`, materializing the standard interface and weighted
baseline aliases without changing the source file.

## Event-level deduplication

`ranked_peptides.evidence_consensus.tsv` intentionally retains every
peptide-HLA row. `ranked_events.evidence_consensus.tsv` is the primary table
for experimental-priority review and applies these deterministic rules:

1. group component events by `phase_group_id` when present, otherwise by
   `event_id`;
2. retain only the best-ranked peptide for each `event_id + HLA` pair;
3. remove repeated overlapping windows sharing `redundancy_group`;
4. retain at most two representative peptide-HLA combinations per event or
   phased event group.

The event table records `member_event_ids`, total peptide count, post-HLA
candidate count, representative count, and complete `representative_1_*` /
`representative_2_*` fields. No peptide rows are deleted from the detailed
peptide table.

All thresholds in `sarcoma_evidence_consensus_v1.toml` are
`PROVISIONAL_RESEARCH_ONLY`. Phase 1 compares algorithm behavior; it does not
claim that these thresholds are clinically validated.

## Authoritative field precedence

`comprehensive_peptide_evidence.tsv` uses field-specific precedence version
`1.0`; it no longer fills only blank cells. The main authorities are:

| Field family | Authority |
| --- | --- |
| Peptide/HLA/event identifiers | annotated peptides, then raw peptides |
| Event/gene/transcript/variant | raw events |
| Presentation predictors | presentation evidence, then annotated tool output |
| Gene/transcript expression | expression evidence |
| RNA ALT/VAF/junction | RNA junction evidence |
| CCF/clonality | CCF 2 |
| APPM | APPM peptide modifiers |
| HLA LOH/escape | peptide escape flags |
| Safety | peptide safety plus namespaced event safety |
| Legacy score/priority | ranked peptides |
| Experimental method | validation plan |

Every row records `evidence_source_precedence_version`, a complete JSON
`evidence_field_sources` map, `comprehensive_evidence_schema_version`,
`evidence_conflict_fields`, and compact conflict details. Different non-empty
values are written to `evidence_source_conflicts.tsv`; the final
`evidence_conflicts.tsv` combines these source conflicts with derived-state
conflicts. Numerically equivalent formatting is not reported as a conflict.

## Normalized states

- Event: `EVENT_CONFIRMED`, `EVENT_STRONG`, `EVENT_PARTIAL`,
  `EVENT_SAMPLE_SPECIFIC`, `EVENT_CONFLICT`, `EVENT_ARTIFACT_RISK`,
  `EVENT_UNASSESSED`.
- RNA: `RNA_CONFIRMED`, `RNA_LOW_SUPPORT`, `GENE_EXPRESSION_ONLY`,
  `RNA_NEGATIVE`, `RNA_UNASSESSED`. Gene TPM alone never proves mutant-allele
  expression.
- Presentation: `PRESENTATION_CONSISTENT_STRONG`, `PRESENTATION_MODERATE`,
  `PRESENTATION_DISCORDANT`, `PRESENTATION_SINGLE_TOOL`, `PRESENTATION_WEAK`,
  `PRESENTATION_UNASSESSED`. NetMHCpan and MHCflurry are the two core groups;
  stability is auxiliary, while PRIME/BigMHC/DeepImmuno form one correlated
  immunogenicity-like group.
- Mutant specificity reuses `MT_SPECIFIC`, `MARGINAL_MT_ADVANTAGE`,
  `MT_WT_SIMILAR`, `WT_BETTER`, `NON_MUTANT_SEQUENCE`, and `UNASSESSED`.
- HLA/APPM: `HLA_APPM_RETAINED`, `HLA_APPM_CAUTION`, `HLA_LOH_UNASSESSED`,
  `RESTRICTING_HLA_LOST`, `MAJOR_APPM_DEFECT`.
- Safety: `SAFETY_PASS`, `SAFETY_PARTIAL`, `SAFETY_REVIEW`,
  `SAFETY_HIGH_RISK`, `SAFETY_REJECT`. Missing normal proteome, ligandome, or
  normal-junction evidence is partial, never pass.

## Hard failures and caps

Hard failures go directly to R4 and retain stable reason codes, including
`HARD_REFERENCE_PROTEOME_MATCH`, `HARD_NORMAL_JUNCTION`,
`HARD_RESTRICTING_HLA_LOST`, `HARD_MATCHED_NORMAL_SUPPORT`,
`HARD_EVENT_ARTIFACT`, `HARD_NON_MUTANT_SEQUENCE`, and
`HARD_SAFETY_REJECT`.

State-driven caps are applied before Pareto ranking. R3 caps include RNA-only
fusion, unassessed normal junction, partial safety, required phasing,
capture-limited WES SV, MT/WT similarity, single-tool presentation, and an
unreproduced source call. R2 caps include unassessed HLA LOH, low APPM evidence,
low-confidence CCF, marginal MT advantage, and unassessed mutant specificity.
`WT_BETTER` and major APPM defects are capped at R4.

## R1-R4 interpretation

The consensus branch uses only R labels; legacy A/B/C/D values remain in
explicitly named legacy columns and never name the new ranking.

- `R1` - first-batch experimental priority. Requires no hard failure; strong or
  confirmed event authenticity; confirmed RNA ALT/junction support; consistent
  NetMHCpan plus MHCflurry support; MT specificity or an explicit novel
  junction/tail; retained restricting HLA without major APPM defect; complete
  safety pass; and supported CCF unless the RNA event track does not require DNA
  CCF.
- `R2` - advance with caution. Event, RNA, and core presentation must still be
  established, with at most one caution such as low-confidence CCF, partial
  APPM evidence, or incomplete safety background.
- `R3` - evidence completion first. Typical triggers are unconfirmed RNA,
  single-tool presentation, RNA-only fusion, weak MT/WT differential,
  unassessed normal junction, required phasing, or partial safety. Suggested
  next steps are targeted RNA, IGV, RT-PCR/Sanger, a second caller/tool group,
  normal-tissue junction review, and paired MT/WT prediction. It does not
  recommend direct ELISpot before those gaps are addressed.
- `R4` - do not advance automatically. Includes hard failure, high safety risk,
  restricting-HLA loss, WT-better evidence, high artifact risk, or consistently
  weak core presentation.

`manual_review_required=yes` may retain mechanism-relevant KRAS, TP53, or
EWSR1::WT1 events for expert review, but it never upgrades their R grade.

The default weights are presentation 0.25, binding 0.20, expression 0.15, and
0.10 each for clonality, tumor specificity, safety, and APPM.

## Run

```bash
neoag evidence-consensus-rank \
  --input results/sample/scoring/comprehensive_peptide_evidence.tsv \
  --output results/sample/scoring/ranked_peptides.evidence_consensus.tsv \
  --rules configs/ranking/sarcoma_evidence_consensus_v1.toml
```

The full pipeline creates this output automatically and records both ranking
paths in provenance. Candidate interpretation should compare the two rankings;
the consensus rank must not be presented as an independently validated model.
