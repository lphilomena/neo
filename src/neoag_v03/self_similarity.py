"""Self-similarity scoring for autoimmune-risk peptide safety checks
(scoring audit fix -- self_similarity_score).

Previously ``self_similarity_score`` was hardcoded to ``"0.0"`` in every
adapter with zero computation behind it anywhere in the codebase. Because
``safety.py`` only flags a candidate as high autoimmune risk when
``wildtype_binding_rank <= 0.5 AND self_similarity_score >= 0.85`` (see
``[safety] self_similarity_caution`` in profiles/*.toml), a permanently-zero
score means that check can never fire -- every candidate silently passes it,
regardless of how similar the mutant peptide actually is to its own
wildtype counterpart.

This module computes a real similarity score between the mutant peptide and
its paired wildtype peptide (both already carried in the PEPTIDE_FIELDS /
PEPTIDE_SAFETY_FIELDS schemas as ``peptide`` and ``wildtype_peptide`` --
for a missense-derived candidate they are the same length, differing at the
mutated position(s)) using the standard BLOSUM62 substitution matrix.

Scope note: this answers "does this mutant peptide look like the *specific*
self-peptide it was derived from" -- the narrow, tractable question that
the existing safety-gate threshold is designed around. It does NOT search
the whole human proteome for the closest-matching *unrelated* self-peptide
(that broader question is what the still-unpopulated
``closest_self_peptide``/``closest_self_gene`` columns in
peptide_safety_gate.py are meant for, and requires a proteome-wide index
that is out of scope here). Peptides with no wildtype counterpart at all
(frameshift, fusion, splice-junction -- genuinely novel junction sequences)
get a score of 0.0, which is the correct answer for *this specific* check:
there is no corresponding normal self-peptide to be similar to.
"""
from __future__ import annotations

from typing import Mapping

# Standard BLOSUM62 substitution matrix (Henikoff & Henikoff 1992), the same
# matrix TESLA's foreignness/agretopicity-style analyses use for peptide
# alignment. Values are log-odds substitution scores; higher = more similar.
_BLOSUM62_RAW = """
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
""".strip().splitlines()


def _build_blosum62() -> dict[tuple[str, str], int]:
    header = _BLOSUM62_RAW[0].split()
    matrix: dict[tuple[str, str], int] = {}
    for line in _BLOSUM62_RAW[1:]:
        parts = line.split()
        row_aa = parts[0]
        for col_aa, val in zip(header, parts[1:]):
            matrix[(row_aa, col_aa)] = int(val)
    return matrix


BLOSUM62 = _build_blosum62()
_DEFAULT_MISMATCH_PENALTY = -4  # score for any pair involving a non-standard/unknown residue


def _pair_score(a: str, b: str) -> int:
    return BLOSUM62.get((a, b), _DEFAULT_MISMATCH_PENALTY)


def _alignment_score(seq_a: str, seq_b: str) -> int:
    """Ungapped, position-by-position BLOSUM62 score. Mutant/wildtype
    peptide pairs from a missense substitution are the same length and
    already positionally aligned (same sliding-window frame), so a full
    gapped alignment (Needleman-Wunsch) isn't needed for this comparison --
    an ungapped sum is both simpler and exactly what's being asked: how
    much does this specific window of sequence still resemble self."""
    return sum(_pair_score(a, b) for a, b in zip(seq_a.upper(), seq_b.upper()))


def compute_self_similarity(
    mutant_peptide: str,
    wildtype_peptide: str | None,
    profile: Mapping | None = None,
) -> float:
    """Return self_similarity_score in [0, 1].

    similarity = alignment_score(mutant, wildtype) / alignment_score(wildtype, wildtype)

    Dividing by the wildtype's self-alignment score (its own maximum
    possible BLOSUM62 score) turns the unbounded log-odds score into an
    intuitive "what fraction of the wildtype's own identity does the mutant
    still retain" value: 1.0 for a peptide identical to its wildtype
    counterpart (e.g. an anchor-distal mutation with no effect on this
    9-mer window), decreasing toward 0 as substitutions accumulate or land
    on residues BLOSUM62 treats as highly dissimilar.

    Returns 0.0 (no autoimmune-risk signal from *this* check) when there is
    no wildtype counterpart at all -- frameshift/fusion/splice-junction
    peptides have no corresponding normal self-peptide by construction.
    """
    mt = (mutant_peptide or "").strip().upper()
    wt = (wildtype_peptide or "").strip().upper()
    if not mt or not wt:
        return 0.0
    if len(mt) != len(wt):
        # Length mismatch (e.g. in-frame indel where the safety-relevant
        # notion of a "paired wildtype peptide" breaks down) -- can't do a
        # meaningful positional BLOSUM62 comparison, so fall back to a
        # coarse identity-fraction over the shared prefix rather than
        # silently returning 0 (which would misleadingly read as "definitely
        # not self-similar").
        shared = min(len(mt), len(wt))
        if shared == 0:
            return 0.0
        matches = sum(1 for a, b in zip(mt[:shared], wt[:shared]) if a == b)
        return round(matches / shared, 4)

    self_score = _alignment_score(wt, wt)
    if self_score <= 0:
        return 0.0
    mt_score = _alignment_score(mt, wt)
    similarity = mt_score / self_score
    return round(max(0.0, min(1.0, similarity)), 4)
