# Tool cross-validation rules

- HLA typing: compare allele calls; preserve disagreement.
- HLA LOH: distinguish class I from class II and apply loss only to the restricting allele.
- Fusion: use one authoritative adjacency-level consensus. The fixed short-read panel is Arriba + STAR-Fusion + FusionCatcher; EasyFuse PASS is an aggregator status, while long-read and DNA-SV are orthogonal confirmation columns. Require exact orientation-compatible breakpoints before combining caller support, preserve caller availability separately from zero calls, and keep targeted rescue outside ordinary PASS consensus. Junction reads, frame, ORF, normal read-through and tissue background remain required biological evidence beyond caller count.
- Splice: distinguish DNA splice-effect prediction, RNA junction observation, and neoepitope reconstruction.
- Presentation: NetMHCpan and MHCflurry are the core evidence group; NetMHCstabpan and NetChop are required production support. PRIME, MixMHCpred, BigMHC and DeepImmuno are advisory immunogenicity support and must not block the main workflow when unavailable.
- CNV/purity/CCF: propagate disagreement as lower confidence rather than forcing one answer.
