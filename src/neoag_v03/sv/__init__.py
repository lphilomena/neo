"""SV Phase 1 / 1.5 support for neoag-v03.

This package converts tumor-normal WGS/WES structural-variant calls into the
existing v0.3 event/peptide schema. It deliberately separates SV-specific
provenance (sv_events.full.tsv, sv_event_to_peptide.tsv) from the generic
raw_events.tsv/raw_peptides.tsv consumed by scoring_v03.
"""

from .phase1 import build_sv_phase1_raw
from .wes_adapter import WESAdapter, build_sv_wes_phase1_5_raw

__all__ = ["build_sv_phase1_raw", "build_sv_wes_phase1_5_raw", "WESAdapter"]
