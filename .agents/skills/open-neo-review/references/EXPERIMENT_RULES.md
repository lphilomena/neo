# Experimental routes

- Missense SNV: paired mutant/wild-type short peptide, ELISpot or tetramer.
- Frameshift: novel-tail long peptide and/or minigene.
- Fusion: breakpoint confirmation, junction-covering long peptide and/or fusion minigene.
- Splice: targeted RNA, abnormal-junction long peptide and/or splice minigene.
- DNA SV: DNA breakpoint plus RNA transcript confirmation before peptide testing.
- Driver/manual-review events remain a separate review lane and are not automatically promoted.

Review priority is evidence-driven:

- R1/R2 with RNA confirmation and safety support may enter the experiment lane.
- Missing RNA evidence routes to `TARGETED_RNA_FIRST`.
- Single-caller fusion routes to `FUSION_CONFIRMATION_FIRST`.
- Unresolved phase routes to `PHASING_FIRST`.
- Partial safety routes to `SAFETY_REVIEW_FIRST`.
- Hard failures remain excluded even when the weighted score is high.

The first-batch set deduplicates event, phase, and redundancy groups and limits HLA/event-type concentration. It is a deterministic research heuristic, not a mathematical treatment optimizer.
