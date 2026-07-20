from __future__ import annotations

import re
from dataclasses import dataclass

_BND_RE = re.compile(r"([\[\]])([^:\[\]]+):(\d+)([\[\]])")


@dataclass(frozen=True)
class BndPartner:
    mate_chrom: str
    mate_pos: int
    local_piece: str
    remote_orientation: str
    alt_form: str
    strand1: str
    strand2: str


def parse_bnd_alt(alt: str) -> BndPartner | None:
    """Parse a VCF breakend ALT allele.

    The four canonical BND patterns encode both mate location and orientation:
      * N]chr2:123]    local sequence followed by reverse-strand remote join
      * N[chr2:123[    local sequence followed by forward-strand remote join
      * ]chr2:123]N    reverse-strand remote join followed by local sequence
      * [chr2:123[N    forward-strand remote join followed by local sequence

    This function does not claim to fully reconstruct graph topology; it returns
    a deterministic orientation approximation used for first-pass annotation.
    """
    if not alt or ("[" not in alt and "]" not in alt):
        return None
    m = _BND_RE.search(alt)
    if not m:
        return None
    left_bracket, chrom, pos, right_bracket = m.groups()
    before = alt[: m.start()]
    after = alt[m.end() :]
    remote_orientation = left_bracket + right_bracket
    # If the mate bracket precedes the local bases, the local retained side is 5'.
    mate_first = m.start() == 0
    local_piece = after if mate_first else before
    # Practical approximation for event classification. Precise reconstruction
    # remains flagged as heuristic downstream.
    if mate_first:
        strand1 = "-" if left_bracket == "]" else "+"
        strand2 = "+" if left_bracket == "]" else "-"
        alt_form = f"mate_first_{left_bracket}{right_bracket}"
    else:
        strand1 = "+" if left_bracket == "]" else "-"
        strand2 = "-" if left_bracket == "]" else "+"
        alt_form = f"local_first_{left_bracket}{right_bracket}"
    return BndPartner(
        mate_chrom=chrom,
        mate_pos=int(pos),
        local_piece=local_piece,
        remote_orientation=remote_orientation,
        alt_form=alt_form,
        strand1=strand1,
        strand2=strand2,
    )
