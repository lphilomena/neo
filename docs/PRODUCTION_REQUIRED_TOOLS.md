# Production required tools

This document distinguishes installation from production evidence. A tool is
not considered complete merely because its executable passes a smoke test; the
sample-specific result must exist, be parseable, and be recorded in the run
manifest.

| Evidence domain | Required run tools | Cross-check / policy | Applicability |
| --- | --- | --- | --- |
| Input identity and DNA QC | samtools, bam-matcher | Both must pass for paired tumor/normal DNA | Paired DNA |
| HLA typing | OptiType, SpecHLA | At least two valid tools; HLA-LA is the preferred third check | No validated external HLA consensus |
| Purity and CNV | FACETS, Sequenza, PURPLE | Run all three; at least two valid estimates. Use the median only when concordant; preserve conflicts and mark purity/CCF low confidence | Paired DNA |
| HLA-I LOH | LOHHLA, SpecHLA | Both must produce comparable allele-level calls. Preserve LOST, RETAINED, UNASSESSED, and DISCORDANT | Paired DNA |
| Variant annotation | Ensembl VEP | Pinned VEP/cache/build must match | SNV/InDel |
| RNA expression | Salmon | RSEM is an optional cross-check | Short-read RNA-seq |
| Fusion discovery | EasyFuse, STAR-Fusion | At least two valid callers; FusionCatcher/Arriba may add support | Short-read RNA-seq |
| Aberrant splice discovery | SNAF, SpliceMutr | Both workflows must run; RegTools/pVACsplice provide supporting evidence | Short-read RNA-seq |
| MHC-I presentation | NetMHCpan, MHCflurry | Both core predictors are required for patient-facing priority candidates | Candidate peptides |
| Immunogenicity | PRIME, BigMHC-IM | Both required; DeepImmuno is a cross-check | Candidate peptides |
| Stability/processing | NetMHCstabpan, NetChop 3.1d | Sample-specific stability and cleavage results are both required; missing or skipped evidence blocks production release | Candidate peptides |

## Hard consensus rules

1. Purity/CNV runs attempt FACETS, Sequenza, and PURPLE. At least two must
   produce valid purity values before CCF can receive reliable positive weight.
2. Concordant purity values use their median. Material disagreement keeps all
   values and marks purity and downstream CCF as `LOW_CONFIDENCE`.
3. HLA-I LOH requires both LOHHLA and SpecHLA. A single-tool result is
   `SINGLE_TOOL_ONLY`, not a confirmed retained/lost call.
4. LOHHLA/SpecHLA disagreement is `DISCORDANT`; the workflow must not silently
   prefer one tool.
5. NetMHCstabpan and NetChop 3.1d are mandatory for production candidate
   ranking. Their authorized packages or approved compatible services must be
   configured before a production run; skipping is allowed only for debugging.
6. Conditional branches that lack their required input are `NOT_APPLICABLE`,
   not `PASS` and not biological negatives.

The authoritative machine-readable policy is
`configs/tools/tools_manifest.yaml:production_requirements`.
