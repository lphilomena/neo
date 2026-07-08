
"""Capture-region utilities for WES SV Phase 1.5.

This module is intentionally dependency-light: it parses BED files, creates
expanded BED sidecars, and classifies SV breakends as ON_TARGET / NEAR_TARGET /
OFF_TARGET. It does not replace caller-specific targeted settings; it provides
the capture-aware evidence layer consumed by the WES adapter.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .sv_merge import SVCluster


@dataclass(frozen=True)
class BedInterval:
    chrom: str
    start0: int
    end0: int

    def contains_1based(self, pos: int) -> bool:
        return self.start0 <= pos - 1 < self.end0

    def distance_1based(self, pos: int) -> int:
        p0 = pos - 1
        if self.start0 <= p0 < self.end0:
            return 0
        if p0 < self.start0:
            return self.start0 - p0
        return p0 - self.end0 + 1


def _parse_bed_line(line: str) -> BedInterval | None:
    if not line.strip() or line.startswith('#'):
        return None
    parts = line.rstrip('\n').split('\t')
    if len(parts) < 3:
        parts = line.split()
    if len(parts) < 3:
        return None
    try:
        return BedInterval(parts[0], max(0, int(float(parts[1]))), max(0, int(float(parts[2]))))
    except Exception:
        return None


class CaptureRegions:
    """In-memory BED intervals with simple distance classification."""

    def __init__(self, intervals: Iterable[BedInterval] = ()):  # noqa: D401
        by_chrom: dict[str, list[BedInterval]] = {}
        for iv in intervals:
            if iv.end0 <= iv.start0:
                continue
            by_chrom.setdefault(iv.chrom, []).append(iv)
        for chrom in list(by_chrom):
            by_chrom[chrom] = sorted(by_chrom[chrom], key=lambda x: (x.start0, x.end0))
        self.by_chrom = by_chrom

    @classmethod
    def from_bed(cls, bed: str | Path | None) -> "CaptureRegions | None":
        if not bed:
            return None
        p = Path(bed)
        if not p.exists():
            raise FileNotFoundError(f"Missing capture BED: {p}")
        intervals = []
        with p.open('r', encoding='utf-8', errors='ignore') as fh:
            for line in fh:
                iv = _parse_bed_line(line)
                if iv:
                    intervals.append(iv)
        return cls(intervals)

    def nearest_distance(self, chrom: str, pos: int) -> int | None:
        intervals = self.by_chrom.get(chrom) or []
        if not intervals:
            # Try alternate chr-prefix spelling; this catches chr1 vs 1 in toy fixtures.
            alt = chrom[3:] if chrom.startswith('chr') else 'chr' + chrom
            intervals = self.by_chrom.get(alt) or []
        if not intervals:
            return None
        best: int | None = None
        for iv in intervals:
            d = iv.distance_1based(pos)
            if best is None or d < best:
                best = d
                if best == 0:
                    return 0
            # Intervals sorted; once we are far beyond pos, no later interval will improve.
            if iv.start0 > pos - 1 and best is not None and iv.start0 - (pos - 1) > best:
                break
        return best

    @staticmethod
    def status_from_distance(distance: int | None, near_bp: int = 250, slop_bp: int = 1000) -> str:
        if distance is None:
            return 'NOT_ASSESSED'
        if distance == 0:
            return 'ON_TARGET'
        if distance <= near_bp:
            return f'NEAR_TARGET_{near_bp}BP'
        if distance <= slop_bp:
            return f'NEAR_TARGET_{slop_bp}BP'
        return 'OFF_TARGET'

    def classify_breakend(self, chrom: str, pos: int, near_bp: int = 250, slop_bp: int = 1000) -> tuple[str, str]:
        d = self.nearest_distance(chrom, pos)
        return self.status_from_distance(d, near_bp, slop_bp), '' if d is None else str(d)

    def classify_cluster(self, cluster: SVCluster, near_bp: int = 250, slop_bp: int = 1000) -> dict[str, str]:
        rep = cluster.representative
        s1, d1 = self.classify_breakend(rep.chrom1, rep.pos1, near_bp, slop_bp)
        s2, d2 = self.classify_breakend(rep.chrom2, rep.pos2, near_bp, slop_bp)
        statuses = {s1, s2}
        if 'ON_TARGET' in statuses or f'NEAR_TARGET_{near_bp}BP' in statuses:
            interpret = 'HIGH'
        elif f'NEAR_TARGET_{slop_bp}BP' in statuses:
            interpret = 'MEDIUM'
        elif statuses == {'OFF_TARGET'} or (s1 == 'OFF_TARGET' and s2 == 'OFF_TARGET'):
            interpret = 'UNINTERPRETABLE'
        elif 'NOT_ASSESSED' in statuses:
            interpret = 'NOT_ASSESSED'
        else:
            interpret = 'LOW'
        return {
            'breakend1_capture_status': s1,
            'breakend2_capture_status': s2,
            'breakend1_capture_distance_bp': d1,
            'breakend2_capture_distance_bp': d2,
            'capture_interpretability': interpret,
        }


def write_expanded_beds(capture_bed: str | Path, outdir: str | Path, *, near_bp: int = 250, slop_bp: int = 1000) -> dict[str, str]:
    """Write sorted/merged-ish target and expanded BED files.

    The implementation avoids external bedtools for portability. It merges only
    overlapping intervals on identical contigs and writes unclipped expanded
    intervals unless a genome file is supplied elsewhere by workflow wrappers.
    """
    regions = CaptureRegions.from_bed(capture_bed)
    if regions is None:
        return {}
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    paths = {
        'targets': out / 'capture.targets.bed',
        f'expanded_{near_bp}bp': out / f'capture.expanded_{near_bp}bp.bed',
        f'expanded_{slop_bp}bp': out / f'capture.expanded_{slop_bp}bp.bed',
    }
    merged: list[BedInterval] = []
    for chrom, ivs in regions.by_chrom.items():
        cur: BedInterval | None = None
        for iv in ivs:
            if cur is None:
                cur = iv
            elif iv.start0 <= cur.end0:
                cur = BedInterval(chrom, cur.start0, max(cur.end0, iv.end0))
            else:
                merged.append(cur); cur = iv
        if cur is not None:
            merged.append(cur)
    def write(path: Path, pad: int = 0):
        with path.open('w', encoding='utf-8') as fh:
            for iv in merged:
                fh.write(f"{iv.chrom}\t{max(0, iv.start0-pad)}\t{iv.end0+pad}\n")
    write(paths['targets'], 0)
    write(paths[f'expanded_{near_bp}bp'], near_bp)
    write(paths[f'expanded_{slop_bp}bp'], slop_bp)
    return {k: str(v) for k, v in paths.items()}


def annotate_cluster_capture(cluster: SVCluster, capture: CaptureRegions | None, *, near_bp: int = 250, slop_bp: int = 1000) -> dict[str, str]:
    """Attach capture annotation to an SVCluster and return it."""
    if capture is None:
        ann = {
            'evidence_scope': 'EXOME_CAPTURE_LIMITED',
            'breakend1_capture_status': 'CAPTURE_NOT_PROVIDED',
            'breakend2_capture_status': 'CAPTURE_NOT_PROVIDED',
            'breakend1_capture_distance_bp': '',
            'breakend2_capture_distance_bp': '',
            'capture_interpretability': 'NOT_ASSESSED',
            'capture_filter_status': 'REVIEW',
            'capture_filter_reason': 'capture_bed_not_provided',
        }
    else:
        ann = capture.classify_cluster(cluster, near_bp=near_bp, slop_bp=slop_bp)
        interpret = ann.get('capture_interpretability', 'NOT_ASSESSED')
        ann['evidence_scope'] = 'EXOME_CAPTURE_LIMITED'
        ann['capture_filter_status'] = 'PASS' if interpret in {'HIGH','MEDIUM'} else ('REVIEW' if interpret == 'LOW' else 'FAIL')
        ann['capture_filter_reason'] = 'capture_near_target' if interpret in {'HIGH','MEDIUM'} else ('low_capture_interpretability' if interpret == 'LOW' else 'outside_capture_scope')
    setattr(cluster, 'capture_annotation', ann)
    return ann


def cluster_capture_annotation(cluster: SVCluster) -> dict[str, str]:
    return getattr(cluster, 'capture_annotation', {}) or {}
