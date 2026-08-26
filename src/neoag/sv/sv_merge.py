from __future__ import annotations

from dataclasses import dataclass, field
from .sv_callset import SVRecord
from .identity import Breakend


def _canonical_endpoints(r: SVRecord) -> tuple[Breakend, Breakend]:
    return tuple(sorted((
        Breakend.create(r.chrom1, r.pos1, r.strand1),
        Breakend.create(r.chrom2, r.pos2, r.strand2),
    )))  # type: ignore[return-value]


def _norm_chrom_pair(r: SVRecord) -> tuple[str, str, int, int, str, str]:
    a, b = _canonical_endpoints(r)
    return (a.chrom, b.chrom, a.pos, b.pos, a.strand, b.strand)


def _compatible(a: SVRecord, b: SVRecord, distance: int = 200, consider_type: bool = True, consider_strand: bool = True) -> bool:
    if consider_type and a.svtype != b.svtype:
        # BND and TRA frequently represent the same event.
        if {a.svtype, b.svtype} <= {"BND", "TRA", "CTX"}:
            pass
        else:
            return False
    a1, a2 = _canonical_endpoints(a)
    b1, b2 = _canonical_endpoints(b)
    if a1.chrom != b1.chrom or a2.chrom != b2.chrom:
        return False
    if abs(a1.pos - b1.pos) > distance or abs(a2.pos - b2.pos) > distance:
        return False
    if consider_strand:
        for x, y in ((a1.strand, b1.strand), (a2.strand, b2.strand)):
            if x in "+-" and y in "+-" and x != y:
                return False
    return True


@dataclass
class SVCluster:
    records: list[SVRecord] = field(default_factory=list)

    @property
    def representative(self) -> SVRecord:
        # Prefer assembly-style callers for breakpoint precision; otherwise max support.
        priority = {"GRIDSS": 3, "GRIDSS2": 3, "SvABA": 2, "svaba": 2, "Manta": 1, "manta": 1}
        return sorted(
            self.records,
            key=lambda r: (priority.get(r.caller, 0), r.tumor_alt_support, -r.breakpoint_precision_bp),
            reverse=True,
        )[0]

    @property
    def callers(self) -> list[str]:
        return sorted({r.caller for r in self.records})

    @property
    def record_ids(self) -> list[str]:
        return [r.record_id for r in self.records]

    def max_int(self, attr: str) -> int:
        return max((int(getattr(r, attr, 0) or 0) for r in self.records), default=0)

    def best_precision(self) -> int:
        vals = [r.breakpoint_precision_bp for r in self.records if r.breakpoint_precision_bp > 0]
        return min(vals) if vals else 0


def cluster_sv_records(records: list[SVRecord], distance: int = 200, consider_type: bool = True, consider_strand: bool = True) -> list[SVCluster]:
    clusters: list[SVCluster] = []
    for rec in sorted(records, key=_norm_chrom_pair):
        placed = False
        for cl in clusters:
            if _compatible(cl.representative, rec, distance, consider_type, consider_strand):
                cl.records.append(rec)
                placed = True
                break
        if not placed:
            clusters.append(SVCluster([rec]))
    return clusters
