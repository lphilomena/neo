# Open-Neo review rules

1. Review is read-only and must use an output directory different from the
   source result directory.
2. Required run manifest, event consensus, peptide consensus, weighted
   baseline, all-tool evidence and validation plan must pass integrity checks.
3. Missing event consensus returns `NEEDS_RANKING`; it never silently falls
   back to a weighted Top20.
4. Preserve `pipeline_r_grade` and `pipeline_event_rank`. Skill3 adds review
   fields but does not rewrite pipeline ranking.
5. Review event-level consensus first. Keep at most one or two representative
   peptide-HLA pairs per event and deduplicate phase/redundancy groups.
6. Hard failures cannot enter R1/R2 or a direct experimental set.
7. Missing RNA, safety, HLA LOH, APPM or CCF evidence is not negative evidence. Fusion/splice caller absence must be reported from Skill2 consensus/provenance fields rather than inferred by rescanning raw caller directories.
8. R1/R2 with adequate evidence may receive direct experiment priority. R3 may
   enter a clearly labelled targeted-RNA, fusion-confirmation or phasing queue.
9. R4 drivers remain manual-review-only and are not upgraded for mechanism
   importance.
10. Outputs describe computational candidates and research validation routes,
    not confirmed neoantigens, treatment benefit, drug advice or a finalized
    vaccine design.
11. Before rendering, inspect the R-grade distribution. An all-R4 result is a
    review warning and becomes blocking when NetMHCpan or MHCflurry coverage is
    absent, because mass fallback values must not be presented as biology.
12. Tool coverage denominators use unique applicable peptide-HLA combinations,
    not all duplicated rows. Report assessed, missing, unsupported and failed
    separately.
13. Show only event tracks present in the formal event ranking. Do not create
    empty DNA-SV, Fusion or Splice sections from consequence text alone.
14. When Splice is present, read `splice_prefilter_funnel.tsv` and distinguish
    strict PASS from prioritized REVIEW. Summarize every UNASSESSED stage; do
    not describe a top-N review selection as full funnel validation.
