# SV hardened implementation

## 1. Intended use

This implementation separates two products:

1. `DNA_SV_event`: a PASS-filtered, canonical genomic adjacency detected from
   tumor/normal DNA. It is a technical event hypothesis.
2. `expressed_rearrangement_product`: an exact-breakpoint RNA-supported,
   externally reconstructed translatable protein product. Only this product may
   generate peptide candidates.

The implementation is research software. It does not make a clinical claim and
does not replace orthogonal validation.

## 2. Fail-closed rules

The hardened defaults are deliberate:

- Multi-sample VCFs require explicit tumor and normal sample names.
- Only records whose VCF FILTER is exactly `PASS` are admitted by default.
- RNA evidence must contain both oriented breakends and a genome build.
- Gene-pair-only RNA evidence is rejected.
- RNA QC requires at least 3 split reads, 2 unique fragment starts, 10 bp
  minimum anchor and MAPQ 20.
- WES mode requires a capture BED.
- DNA-only heuristic protein reconstructions are labelled
  `ORF_HYPOTHESIS_UNRESOLVED` and cannot generate peptides.
- Peptide identifiers use SHA-256 and are stable across processes.

## 3. Exact RNA junction schema

Required columns:

```text
genome_build chrom1 pos1 strand1 chrom2 pos2 strand2
split_reads unique_start_count min_anchor_bp min_mapq
```

Recommended provenance columns:

```text
spanning_reads source_tool source_record_id
```

Coordinates must use the same reference build as the SV VCF. Breakend order is
irrelevant because the pipeline canonicalizes the pair. Strand/orientation is
part of the identity and is not optional.

See `data/fixtures_sv/rna_junctions_exact.tsv`.

## 4. Confirmed expressed-product schema

Required columns:

```text
genome_build chrom1 pos1 strand1 chrom2 pos2 strand2
gene1 gene2 transcript1 transcript2 protein_sequence
junction_aa_position in_frame orf_status
```

`orf_status` must be `CONFIRMED`, `COMPLETE`, or `TRANSLATABLE`.
`junction_aa_position` is a zero-based boundary in the supplied protein: the
left product ends immediately before that position and the right/novel product
begins at that position.

Recommended columns:

```text
wildtype_protein_sequence source_tool source_record_id
```

The table should be produced by a versioned Arriba/STAR-Fusion plus AGFusion,
EasyFuse, local transcript assembly, long-read transcript, or equivalent
adapter. The upstream adapter remains responsible for validating the exact
exon structure, translation direction, CDS phase, inserted bases and ORF.

See `data/fixtures_sv/expressed_products.tsv`.

## 5. WGS command

```bash
neoag sv-build-raw \
  --sample-id PATIENT01 \
  --sv-vcf manta.vcf.gz svaba.vcf.gz gridss.vcf.gz \
  --callers Manta SvABA GRIDSS2 \
  --tumor-sample-name TUMOR \
  --normal-sample-name NORMAL \
  --genome-build GRCh38 \
  --reference-fasta GRCh38.fa \
  --gencode-gtf gencode.gtf \
  --hla hla.txt \
  --rna-junctions exact_rna_junctions.tsv \
  --expressed-products expressed_products.tsv \
  --expression gene_tpm.tsv \
  --normal-expression normal_expression.tsv \
  --normal-hla-ligands normal_ligands.tsv \
  --outdir results/PATIENT01_sv
```

Omit `--rna-junctions` and `--expressed-products` to build a DNA event catalogue
only. In that mode, `raw_peptides.tsv` is intentionally empty.

## 6. WES command

```bash
neoag sv-build-raw-wes \
  --sample-id PATIENT01 \
  --sv-vcf wes_sv.vcf.gz \
  --callers Manta \
  --tumor-sample-name TUMOR \
  --normal-sample-name NORMAL \
  --genome-build GRCh38 \
  --capture-bed capture_targets.bed \
  --reference-fasta GRCh38.fa \
  --gencode-gtf gencode.gtf \
  --hla hla.txt \
  --rna-junctions exact_rna_junctions.tsv \
  --expressed-products expressed_products.tsv \
  --outdir results/PATIENT01_wes_sv
```

WES is a rescue/triage route. Capture context and RNA evidence do not turn WES
into genome-wide SV detection.

## 7. Outputs and interpretation

- `sv/sv_events.full.tsv`: all admitted DNA event hypotheses, including
  `adjacency_key`, exact RNA match/QC and reconstruction state.
- `sv/sv_protein_reconstruction.tsv`: only confirmed expressed products.
- `sv/sv_event_to_peptide.tsv`: residue-window provenance for generated peptides.
- `parsed/raw_events.tsv` and `parsed/raw_peptides.tsv`: downstream scoring inputs.
- `provenance.sv_phase1.json` or `provenance.sv_wes_phase1_5.json`: run settings.

Important states:

- `hypothesis_only`: DNA event retained; no validated product and no peptides.
- `confirmed_expressed_product`: exact product supplied and eligible for the RNA
  QC gate.
- `RNA_JUNCTION_LOW_QC`: matching row exists but fails split-read/uniqueness/
  anchor/MAPQ thresholds.
- `NO_EXACT_MATCH`: no evidence with the exact canonical adjacency.

## 8. What remains intentionally unsolved

The code does not claim to reconstruct derivative chromosomes, chromothripsis,
templated insertions, chained rearrangements or arbitrary inversion products.
Those require a graph/assembly layer and must not be approximated by nearest-CDS
splicing. The integration point is the confirmed expressed-product TSV.

Before production research use, add caller-version-specific filters, PoN,
population-SV, repeat/mappability and blacklist annotations; benchmark against
truth sets and representative real caller VCFs; and perform RT-PCR/Sanger or
another orthogonal assay for selected products.

## 9. Validation performed for this delivery

- Python compile check for all `src` and `tests` modules.
- WGS adapter smoke run through event, product and peptide generation.
- WES adapter smoke run with capture BED.
- Full scoring smoke run with the repository's explicit binding stubs.
- Direct assertions for sample identity, VCF FILTER, swapped breakend
  canonicalization, gene-pair evidence rejection, deterministic IDs and
  DNA-only peptide suppression.

The execution environment did not contain pytest and network installation was
unavailable, so the complete pytest suite was not executed here. The new tests
are included in `tests/test_sv_hardening.py` for execution in the normal test
environment.

