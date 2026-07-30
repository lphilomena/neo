# NeoAg v0.5.0 Formal Splice Provenance Layer

## 1. Scope

v0.5.0 converts splice-derived candidates from loosely connected caller tables into an exact five-level provenance graph:

```text
canonical junction
  → biological splice event
  → transcript hypothesis
  → ORF / translated segment
  → peptide origin
  → peptide-HLA presentation
```

The authoritative tables are the formal `splice_*.tsv` outputs. `raw_events.tsv`, `raw_peptides.tsv`, and `rna_junction_evidence.tsv` remain compatibility projections for the existing NeoAg ranking and report layers.

This is a research pipeline. A computationally generated ORF or peptide is not an experimentally validated transcript, translated product, HLA ligand, or T-cell target.

## 2. Core invariants

1. Junction identity is exact: `genome build + chromosome + 1-based closed intron start + intron end + strand`.
2. No gene-name, nearest-coordinate, or largest-read fallback is permitted.
3. Every foreign-key relation uses a stable formal identifier.
4. pVACbind results map only through the exact generated FASTA `Index` and an exact epitope occurrence in that ORF.
5. Normal non-detection becomes supportive only when locus coverage is explicitly adequate.
6. A partial ImmunoPepper translation is labelled as a partial translated segment, not a confirmed full-length ORF.
7. Conflicts are materialized in `splice_conflicts.tsv`; they are never silently overwritten.
8. Multi-tool agreement is evaluated by independent evidence groups, not raw tool count.

## 3. Stable identifiers

```text
SJ|GRCh38|chr1|151|200|+      canonical junction
SEV|<digest>                   biological splice event
STH|<digest>                   transcript hypothesis
ORF|<digest>                   ORF / translated segment
PEP|<digest>                   unique amino-acid sequence
POR|<digest>                   peptide sequence in a specific ORF/event origin
PRE|<digest>                   peptide-HLA presentation result
```

IDs are deterministic from biological identity fields. Caller-local names remain provenance only.

## 4. Supported inputs

### Exact RNA junction evidence

- RegTools annotated tables;
- RegTools BED12 `junctions extract` output;
- STAR `SJ.out.tab`.

### Splice-event graph

- SplAdder confirmed-event GFF3;
- SplAdder build-mode headerless event TXT;
- SplAdder headered test-mode TSV/TXT.

SplAdder HDF5 is deliberately not parsed in v0.5.0. Export GFF3 or TXT first.

### Intron retention

- IRFinder-S result tables, with an explicit coordinate-system declaration.

### Translation hypotheses

- ImmunoPepper `*_sample_peptides_meta`;
- ImmunoPepper junction-kmer outputs.

Official semicolon-separated `modifiedExonsCoord` and site-exported `chr:start-end` forms are supported. `isJunctionList` is retained as source context and is not misinterpreted as “crosses junction”; crossing status derives from exon structure and `isIsolated`.

### Presentation

- pVACbind `*.all_epitopes.tsv` for MHC-I, MHC-II, or combined runs.

The parser accepts pVACtools 6.x/7.x header aliases, but the producing pVACtools version must be recorded and locked in a production deployment.

### Normal background

- exact normal-junction files;
- explicit coverage-aware normal assessment tables.

## 5. Authoritative outputs

```text
splice_junctions.tsv
splice_events.tsv
splice_event_junction_links.tsv
splice_transcript_hypotheses.tsv
splice_orfs.tsv
splice_peptide_origins.tsv
splice_peptide_origin_links.tsv
splice_pvacbind_predictions.tsv
splice_normal_background.tsv
splice_tool_evidence.long.tsv
splice_consensus.tsv
splice_conflicts.tsv
splice_qc.tsv
splice_pvacbind_input.fasta
splice_pvacbind_fasta_map.tsv
provenance_manifest.json
```

Compatibility outputs:

```text
raw_events.tsv
raw_peptides.tsv
rna_junction_evidence.tsv
```

## 6. Two-pass production run

The recommended driver is:

```bash
bash scripts/run_splice_provenance_v050.sh \
  --sample-id PATIENT_001 \
  --outdir results/PATIENT_001/splice_v050 \
  --genome-build GRCh38 \
  --junctions regtools_junctions.tsv \
  --junction-coordinate-system bed12 \
  --spladder-gff3 merge_graphs_exon_skip_C3.confirmed.gff3 \
  --irfinder IRFinder-IR-nondir.txt \
  --irfinder-coordinate-system intron_1based_closed \
  --immunopepper-meta ref_sample_peptides_meta.tsv \
  --immunopepper-kmers ref_graph_kmer_JuncExpr.tsv \
  --normal-junctions protocol_matched_normal_junctions.tsv \
  --normal-coverage normal_coverage.tsv \
  --hla-file hla.txt \
  --pvacbind-algorithms MHCflurry,MHCflurryEL \
  --tool-version STAR=2.7.11b \
  --tool-version RegTools=1.0.0 \
  --tool-version SplAdder=3.1.1 \
  --tool-version IRFinder-S=2.0.1 \
  --tool-version ImmunoPepper=<locked-commit> \
  --tool-version pVACbind=7.1.1 \
  --strict
```

The driver performs:

1. formal layer build without presentation;
2. exact ORF FASTA and FASTA-map generation;
3. optional pVACbind execution;
4. final rebuild with exact FASTA-index mapping;
5. manifest/hash and referential-integrity validation.

To build from already generated pVACbind files, use the Python CLI directly:

```bash
neoag-splice-layer build \
  --sample-id PATIENT_001 \
  --outdir results/PATIENT_001/splice_layer \
  --junctions regtools_junctions.tsv \
  --spladder-gff3 merge_graphs_exon_skip_C3.confirmed.gff3 \
  --immunopepper-meta ref_sample_peptides_meta.tsv \
  --pvacbind-fasta-map pre_pvacbind/splice_pvacbind_fasta_map.tsv \
  --pvacbind pvacbind/MHC_Class_I/PATIENT_001.MHC_I.all_epitopes.tsv \
  --strict
```

## 7. Evidence grades

### Event evidence

- `E0`: no exact RNA confirmation;
- `E1`: exact split-read support;
- `E2`: exact RNA support plus an independent event model such as SplAdder or IRFinder-S;
- `E3`: DNA-causal, long-read, or equivalent higher-order evidence.

### ORF evidence

- `O0`: invalid or unresolved;
- `O1`: one valid generator;
- `O2`: independent generators agree on event, frame, and protein sequence;
- `O3`: long-read/protein/ligandome validation.

### Normal safety

- `N0`: detected in normal background;
- `N1`: incomplete assessment;
- `N2`: one adequate-coverage negative source;
- `N3`: at least two independent adequate-coverage negative sources.

Final `R1–R4` tiers are constrained by caps and hard failures. Unknown peptide novelty is a review cap, not an automatic assertion of novelty. Explicitly non-junction-crossing and non-novel candidates can receive `HARD_NO_NOVEL_AMINO_ACID`.

## 8. pVACbind mapping safeguards

A presentation row is accepted only when:

1. `Index` maps to exactly one generated FASTA record;
2. that FASTA record maps to exactly one ORF/event chain;
3. the epitope occurs exactly in the mapped ORF;
4. the reported position matches, or a unique exact sequence occurrence allows a recorded correction.

Unmapped, ambiguous, or sequence-inconsistent results are written to the conflict table and contribute no presentation evidence.

## 9. Known boundaries

- ImmunoPepper output represents local two-/three-exon translations; v0.5.0 does not relabel these as full-length transcripts.
- SplAdder HDF5, long-read isoform callers, moPepGen, splice2neo, EasyQuant, pVACsplice, and k4neo are scheduled for later layers/releases or remain external evidence providers.
- IRFinder coordinate semantics must be declared by the site; do not rely on filename inference.
- Thresholds and evidence tiers are engineering defaults requiring calibration on validated datasets.
- External tools, reference assets, HLA predictors, and licensed software are not bundled in the lightweight source release.
