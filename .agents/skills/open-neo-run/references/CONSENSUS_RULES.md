# Tool cross-validation rules

- HLA typing: compare allele calls; preserve disagreement.
- HLA LOH: distinguish class I from class II and apply loss only to the restricting allele.
- Fusion: use caller support, junction reads, frame, normal read-through, and tissue background rather than caller count alone.
- Splice: distinguish DNA splice-effect prediction, RNA junction observation, and neoepitope reconstruction.
- Presentation: NetMHCpan and MHCflurry are the core evidence group; NetMHCstabpan and NetChop are required production support. PRIME, MixMHCpred, BigMHC and DeepImmuno are advisory immunogenicity support and must not block the main workflow when unavailable.
- CNV/purity/CCF: propagate disagreement as lower confidence rather than forcing one answer.
