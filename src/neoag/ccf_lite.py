from __future__ import annotations
from .ccf_v2 import build_ccf_2, load_purity, find_cn, estimate_ccf, clonality_status as status, clonality_multiplier as mult


def build_ccf_lite(events_tsv, purity_tsv, cnv_tsv, profile, out):
    """Backward-compatible CCF-lite entry point backed by CCF 2.0.

    Existing score still consumes ccf_estimate/ccf_status/clonality_multiplier.
    New outputs additionally include ccf_best/min/max, copy-number context, method
    and confidence fields.
    """
    return build_ccf_2(events_tsv, purity_tsv, cnv_tsv, profile, out)
