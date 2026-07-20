from __future__ import annotations

import gzip
from dataclasses import dataclass, field
from pathlib import Path


def open_text(path: str | Path):
    p = Path(path)
    if str(p).endswith(".gz"):
        return gzip.open(p, "rt", encoding="utf-8", errors="ignore")
    return p.open("r", encoding="utf-8", errors="ignore")


RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(RC_TABLE)[::-1].upper()


class FastaReference:
    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.seqs = self._load(self.path)

    @staticmethod
    def _load(path: Path) -> dict[str, str]:
        seqs: dict[str, list[str]] = {}
        name: str | None = None
        with open_text(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    name = line[1:].strip().split()[0]
                    seqs.setdefault(name, [])
                elif name:
                    seqs[name].append(line.strip())
        return {k: "".join(v).upper() for k, v in seqs.items()}

    def fetch(self, chrom: str, start: int, end: int, strand: str = "+") -> str:
        if chrom not in self.seqs:
            alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
            if alt in self.seqs:
                chrom = alt
            else:
                raise KeyError(f"Contig not found in FASTA: {chrom}")
        s = max(1, int(start))
        e = min(len(self.seqs[chrom]), int(end))
        if e < s:
            return ""
        seq = self.seqs[chrom][s - 1 : e]
        return revcomp(seq) if strand == "-" else seq


@dataclass
class Feature:
    chrom: str
    start: int
    end: int
    strand: str
    feature_type: str
    gene_id: str
    gene_name: str
    transcript_id: str = ""
    exon_number: str = ""
    phase: str = "."


@dataclass
class TranscriptModel:
    transcript_id: str
    gene_id: str
    gene_name: str
    chrom: str
    strand: str
    cds: list[Feature] = field(default_factory=list)
    exons: list[Feature] = field(default_factory=list)

    def sorted_cds(self) -> list[Feature]:
        if self.strand == "-":
            return sorted(self.cds, key=lambda f: f.start, reverse=True)
        return sorted(self.cds, key=lambda f: f.start)

    def sorted_exons(self) -> list[Feature]:
        if self.strand == "-":
            return sorted(self.exons, key=lambda f: f.start, reverse=True)
        return sorted(self.exons, key=lambda f: f.start)

    def cds_sequence_and_map(self, ref: FastaReference) -> tuple[str, list[int]]:
        seq_parts: list[str] = []
        pos_map: list[int] = []
        for f in self.sorted_cds():
            if self.strand == "-":
                seq_parts.append(ref.fetch(f.chrom, f.start, f.end, "-"))
                pos_map.extend(list(range(f.end, f.start - 1, -1)))
            else:
                seq_parts.append(ref.fetch(f.chrom, f.start, f.end, "+"))
                pos_map.extend(list(range(f.start, f.end + 1)))
        return "".join(seq_parts).upper(), pos_map


def parse_gtf_attributes(attr: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for raw in attr.strip().rstrip(";").split(";"):
        raw = raw.strip()
        if not raw:
            continue
        if " " in raw:
            k, v = raw.split(" ", 1)
            out[k] = v.strip().strip('"')
        elif "=" in raw:
            k, v = raw.split("=", 1)
            out[k] = v.strip().strip('"')
    return out


class GtfAnnotation:
    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.features: list[Feature] = []
        self.transcripts: dict[str, TranscriptModel] = {}
        self.gene_features: list[Feature] = []
        self._parse()

    def _parse(self) -> None:
        with open_text(self.path) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                chrom, _source, ftype, start, end, _score, strand, phase, attrs = parts
                ad = parse_gtf_attributes(attrs)
                gene_id = ad.get("gene_id", "")
                gene_name = ad.get("gene_name", gene_id)
                transcript_id = ad.get("transcript_id", "")
                feat = Feature(
                    chrom=chrom,
                    start=int(start),
                    end=int(end),
                    strand=strand,
                    feature_type=ftype,
                    gene_id=gene_id,
                    gene_name=gene_name,
                    transcript_id=transcript_id,
                    exon_number=ad.get("exon_number", ""),
                    phase=phase,
                )
                if ftype == "gene":
                    self.gene_features.append(feat)
                if ftype in {"CDS", "exon"} and transcript_id:
                    tx = self.transcripts.get(transcript_id)
                    if not tx:
                        tx = TranscriptModel(transcript_id, gene_id, gene_name, chrom, strand)
                        self.transcripts[transcript_id] = tx
                    if ftype == "CDS":
                        tx.cds.append(feat)
                    elif ftype == "exon":
                        tx.exons.append(feat)
                    self.features.append(feat)

    def features_at(self, chrom: str, pos: int, feature_types: set[str] | None = None) -> list[Feature]:
        chroms = {chrom, chrom[3:] if chrom.startswith("chr") else "chr" + chrom}
        out = []
        for f in self.features:
            if f.chrom in chroms and f.start <= pos <= f.end and (feature_types is None or f.feature_type in feature_types):
                out.append(f)
        return out

    def gene_at(self, chrom: str, pos: int) -> str:
        feats = self.features_at(chrom, pos, {"CDS", "exon"})
        if feats:
            cds = [f for f in feats if f.feature_type == "CDS"]
            chosen = (cds or feats)[0]
            return chosen.gene_name or chosen.gene_id
        chroms = {chrom, chrom[3:] if chrom.startswith("chr") else "chr" + chrom}
        for f in self.gene_features:
            if f.chrom in chroms and f.start <= pos <= f.end:
                return f.gene_name or f.gene_id
        return ""

    def best_transcript_at(self, chrom: str, pos: int) -> TranscriptModel | None:
        feats = self.features_at(chrom, pos, {"CDS"}) or self.features_at(chrom, pos, {"exon"})
        if not feats:
            return None
        tids = {f.transcript_id for f in feats if f.transcript_id in self.transcripts}
        if not tids:
            return None
        # Prefer the transcript with the longest CDS.
        return sorted((self.transcripts[t] for t in tids), key=lambda tx: sum(f.end - f.start + 1 for f in tx.cds), reverse=True)[0]

    @staticmethod
    def _tx_cds_len(tx: TranscriptModel) -> int:
        return sum(f.end - f.start + 1 for f in tx.cds)

    @staticmethod
    def _tx_distance(tx: TranscriptModel, pos: int) -> int:
        spans = tx.exons or tx.cds
        if not spans:
            return 10**12
        if any(f.start <= pos <= f.end for f in spans):
            return 0
        return min(abs(pos - f.start) if pos < f.start else abs(pos - f.end) for f in spans)

    def best_coding_transcript_near(self, chrom: str, pos: int, max_distance_bp: int = 200_000) -> tuple[TranscriptModel | None, int, str]:
        """Return a coding transcript at, within the same gene as, or near a breakpoint.

        Long-read SV breakpoints often fall in introns rather than CDS/exons. For
        first-pass neoantigen triage, mapping to the closest coding transcript of
        the same gene preserves candidate fusions while keeping the reconstruction
        explicitly heuristic.
        """
        exact = self.best_transcript_at(chrom, pos)
        if exact:
            return exact, 0, "overlapping_cds_or_exon"
        chroms = {chrom, chrom[3:] if chrom.startswith("chr") else "chr" + chrom}
        gene = self.gene_at(chrom, pos)
        candidates = [
            tx for tx in self.transcripts.values()
            if tx.chrom in chroms and tx.cds and (not gene or tx.gene_name == gene or tx.gene_id == gene)
        ]
        if candidates:
            best = sorted(candidates, key=lambda tx: (self._tx_distance(tx, pos), -self._tx_cds_len(tx)))[0]
            return best, self._tx_distance(best, pos), "same_gene_nearest_coding_transcript"
        nearby = [
            tx for tx in self.transcripts.values()
            if tx.chrom in chroms and tx.cds and self._tx_distance(tx, pos) <= max_distance_bp
        ]
        if not nearby:
            return None, 0, "no_nearby_coding_transcript"
        best = sorted(nearby, key=lambda tx: (self._tx_distance(tx, pos), -self._tx_cds_len(tx)))[0]
        return best, self._tx_distance(best, pos), "nearby_coding_transcript"

    def transcript_by_id(self, transcript_id: str) -> TranscriptModel | None:
        return self.transcripts.get(transcript_id)
