from __future__ import annotations

import gzip
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Any

from .bnd_parser import parse_bnd_alt


def open_text(path: str | Path):
    p = Path(path)
    if str(p).endswith(".gz"):
        return gzip.open(p, "rt", encoding="utf-8", errors="ignore")
    return p.open("r", encoding="utf-8", errors="ignore")


def parse_info(info: str) -> dict[str, str | bool]:
    d: dict[str, str | bool] = {}
    if not info or info == ".":
        return d
    for part in info.split(";"):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
        else:
            d[part] = True
    return d


def _as_int(value: Any, default: int = 0) -> int:
    try:
        if value is None:
            return default
        s = str(value).strip()
        if s in {"", ".", "NA", "nan"}:
            return default
        if "," in s:
            return int(float(s.split(",")[-1]))
        return int(float(s))
    except Exception:
        return default


def _last_int(value: str | None) -> int:
    if not value:
        return 0
    if value in {".", ""}:
        return 0
    if "," in value:
        return _as_int(value.split(",")[-1], 0)
    return _as_int(value, 0)


def _sum_ints(value: str | None) -> int:
    if not value or value in {".", ""}:
        return 0
    total = 0
    for part in str(value).replace("|", ",").split(","):
        total += _as_int(part, 0)
    return total


def _count_list(value: Any) -> int:
    if value in {None, "", ".", "NA"}:
        return 0
    return len([x for x in str(value).replace("|", ",").split(",") if x.strip()])


def parse_format(fmt: str, sample: str) -> dict[str, str]:
    if not fmt or fmt == "." or not sample:
        return {}
    keys = fmt.split(":")
    vals = sample.split(":")
    return {k: vals[i] if i < len(vals) else "" for i, k in enumerate(keys)}


def evidence_from_sample(sample_map: dict[str, str]) -> tuple[int, int, int, int]:
    """Return (sr, pe, alt_support, depth) from heterogeneous SV FORMAT fields."""
    sr = 0
    pe = 0
    # Manta commonly uses PR/SR as ref,alt pairs. DELLY may use DV/RV.
    if "SR" in sample_map:
        sr = _last_int(sample_map.get("SR"))
    elif "RV" in sample_map:
        sr = _as_int(sample_map.get("RV"), 0)
    elif "SRQ" in sample_map:
        sr = 0

    if "PR" in sample_map:
        pe = _last_int(sample_map.get("PR"))
    elif "PE" in sample_map:
        pe = _last_int(sample_map.get("PE"))
    elif "DV" in sample_map:
        pe = _as_int(sample_map.get("DV"), 0)

    alt = sr + pe
    if alt <= 0:
        for key in ("AD", "AO", "TIR", "TAR", "DR", "DV"):
            if key in sample_map:
                if key == "AD":
                    alt = _last_int(sample_map.get(key))
                elif key in {"AO", "DV"}:
                    alt = _sum_ints(sample_map.get(key))
                break
    depth = _as_int(sample_map.get("DP"), 0)
    if depth <= 0:
        # Approximate depth from ref/alt evidence where possible.
        if "AD" in sample_map:
            depth = _sum_ints(sample_map.get("AD"))
        elif "PR" in sample_map or "SR" in sample_map:
            depth = _sum_ints(sample_map.get("PR")) + _sum_ints(sample_map.get("SR"))
        elif "DR" in sample_map or "DV" in sample_map:
            depth = _as_int(sample_map.get("DR"), 0) + _as_int(sample_map.get("DV"), 0)
    if alt <= 0:
        alt = sr + pe
    return sr, pe, alt, depth


@dataclass
class SVRecord:
    record_id: str
    caller: str
    chrom1: str
    pos1: int
    chrom2: str
    pos2: int
    svtype: str
    alt: str = ""
    ref: str = "N"
    qual: str = "."
    filt: str = "."
    strand1: str = "."
    strand2: str = "."
    cipos: str = ""
    ciend: str = ""
    svlen: str = ""
    inserted_sequence: str = ""
    sniffles_support: int = 0
    sniffles_coverage: int = 0
    sniffles_precise: str = ""
    sniffles_rnames_count: int = 0
    tumor_sr: int = 0
    tumor_pe: int = 0
    tumor_alt_support: int = 0
    tumor_local_depth: int = 0
    normal_sr: int = 0
    normal_pe: int = 0
    normal_alt_support: int = 0
    normal_local_depth: int = 0
    info: dict[str, str | bool] = field(default_factory=dict)
    source_path: str = ""

    @property
    def breakpoint_precision_bp(self) -> int:
        vals: list[int] = []
        for raw in (self.cipos, self.ciend):
            if not raw:
                continue
            try:
                parts = [int(float(x)) for x in str(raw).split(",") if x not in {"", "."}]
                if len(parts) >= 2:
                    vals.append(abs(parts[1] - parts[0]))
            except Exception:
                pass
        return max(vals) if vals else 0

    @property
    def vaf_like(self) -> float:
        if self.tumor_local_depth > 0:
            return max(0.0, min(1.0, self.tumor_alt_support / self.tumor_local_depth))
        denom = self.tumor_alt_support + self.normal_alt_support
        if denom > 0:
            return max(0.0, min(1.0, self.tumor_alt_support / denom))
        return 0.0


def _select_sample_indices(sample_names: list[str], tumor_sample_name: str | None, normal_sample_name: str | None) -> tuple[int | None, int | None]:
    if not sample_names:
        return None, None
    tumor_i = sample_names.index(tumor_sample_name) if tumor_sample_name in sample_names else None
    normal_i = sample_names.index(normal_sample_name) if normal_sample_name in sample_names else None
    if tumor_i is None:
        # Conservative default for Manta tumor-normal VCFs is often NORMAL,TUMOR.
        tumor_i = len(sample_names) - 1
    if normal_i is None and len(sample_names) >= 2:
        normal_i = 0 if tumor_i != 0 else 1
    return tumor_i, normal_i


def parse_vcf_records(
    path: str | Path,
    caller: str,
    *,
    tumor_sample_name: str | None = None,
    normal_sample_name: str | None = None,
) -> list[SVRecord]:
    records: list[SVRecord] = []
    sample_names: list[str] = []
    tumor_i: int | None = None
    normal_i: int | None = None
    with open_text(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                sample_names = header[9:]
                tumor_i, normal_i = _select_sample_indices(sample_names, tumor_sample_name, normal_sample_name)
                continue
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, rid, ref, alt, qual, filt, info_raw = parts[:8]
            info = parse_info(info_raw)
            svtype = str(info.get("SVTYPE") or "")
            bnd = parse_bnd_alt(alt)
            chrom2 = str(info.get("CHR2") or info.get("MATE_CHR") or "")
            pos2 = int(float(str(info.get("END") or pos)))
            strand1 = "."
            strand2 = "."
            if bnd:
                chrom2 = chrom2 or bnd.mate_chrom
                pos2 = bnd.mate_pos
                strand1, strand2 = bnd.strand1, bnd.strand2
                svtype = svtype or "BND"
            if not chrom2:
                chrom2 = chrom
            if not svtype:
                if alt.startswith("<") and alt.endswith(">"):
                    svtype = alt.strip("<>")
                elif bnd:
                    svtype = "BND"
                else:
                    svtype = "SV"
            strands = str(info.get("STRANDS") or info.get("CT") or "")
            if strands and strand1 == ".":
                if len(strands) >= 2 and strands[0] in "+-" and strands[1] in "+-":
                    strand1, strand2 = strands[0], strands[1]
            inserted = str(info.get("INSSEQ") or info.get("INSERTED_SEQUENCE") or info.get("SEQ") or "")
            sniffles_support = _as_int(info.get("SUPPORT"), 0)
            sniffles_coverage = _as_int(info.get("COVERAGE"), 0)
            sniffles_precise = "yes" if "PRECISE" in info else ("no" if "IMPRECISE" in info else "")
            sniffles_rnames_count = _count_list(info.get("RNAMES"))
            fmt_map_t: dict[str, str] = {}
            fmt_map_n: dict[str, str] = {}
            if len(parts) >= 10:
                fmt = parts[8]
                samples = parts[9:]
                if tumor_i is not None and tumor_i < len(samples):
                    fmt_map_t = parse_format(fmt, samples[tumor_i])
                if normal_i is not None and normal_i < len(samples):
                    fmt_map_n = parse_format(fmt, samples[normal_i])
            t_sr, t_pe, t_alt, t_depth = evidence_from_sample(fmt_map_t)
            n_sr, n_pe, n_alt, n_depth = evidence_from_sample(fmt_map_n)
            # Fallback to INFO fields used by some SV callers.
            if t_alt <= 0:
                t_sr = _as_int(info.get("SR"), t_sr)
                t_pe = _as_int(info.get("PE"), t_pe)
                t_alt = t_sr + t_pe
            if t_alt <= 0 and sniffles_support > 0:
                # Sniffles2 stores long-read support in INFO/SUPPORT rather than SR/PE.
                t_pe = sniffles_support
                t_alt = sniffles_support
            if t_depth <= 0 and sniffles_coverage > 0:
                t_depth = sniffles_coverage
            record_id = rid if rid and rid != "." else f"{caller}_{chrom}_{pos}_{svtype}_{len(records)+1}"
            records.append(SVRecord(
                record_id=record_id,
                caller=caller,
                chrom1=chrom,
                pos1=int(float(pos)),
                chrom2=chrom2,
                pos2=pos2,
                svtype=str(svtype),
                alt=alt,
                ref=ref,
                qual=qual,
                filt=filt,
                strand1=strand1,
                strand2=strand2,
                cipos=str(info.get("CIPOS") or ""),
                ciend=str(info.get("CIEND") or ""),
                svlen=str(info.get("SVLEN") or ""),
                inserted_sequence=inserted if inserted not in {".", "NA"} else "",
                sniffles_support=sniffles_support,
                sniffles_coverage=sniffles_coverage,
                sniffles_precise=sniffles_precise,
                sniffles_rnames_count=sniffles_rnames_count,
                tumor_sr=t_sr,
                tumor_pe=t_pe,
                tumor_alt_support=t_alt,
                tumor_local_depth=t_depth,
                normal_sr=n_sr,
                normal_pe=n_pe,
                normal_alt_support=n_alt,
                normal_local_depth=n_depth,
                info=info,
                source_path=str(path),
            ))
    return records


def read_sv_inputs(
    sv_vcfs: Iterable[str | Path],
    callers: Iterable[str] | None = None,
    *,
    tumor_sample_name: str | None = None,
    normal_sample_name: str | None = None,
) -> list[SVRecord]:
    paths = [Path(p) for p in sv_vcfs]
    caller_names = list(callers or [])
    out: list[SVRecord] = []
    for i, p in enumerate(paths):
        caller = caller_names[i] if i < len(caller_names) and caller_names[i] else p.stem.replace(".vcf", "")
        out.extend(parse_vcf_records(p, caller, tumor_sample_name=tumor_sample_name, normal_sample_name=normal_sample_name))
    return out
