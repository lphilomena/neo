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
7. Missing RNA, safety, HLA LOH, APPM or CCF evidence is not negative evidence.
8. R1/R2 with adequate evidence may receive direct experiment priority. R3 may
   enter a clearly labelled targeted-RNA, fusion-confirmation or phasing queue.
9. R4 drivers remain manual-review-only and are not upgraded for mechanism
   importance.
10. Outputs describe computational candidates and research validation routes,
    not confirmed neoantigens, treatment benefit, drug advice or a finalized
    vaccine design.
