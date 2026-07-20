
"""Peptide-level safety evidence resolver.

The gate is deliberately conservative and explainable. It writes a sidecar
`peptide_safety.tsv` that score can consume, instead of hiding safety logic
inside a single scalar ranking score.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from .utils import read_tsv, write_tsv, to_float, first
from .safety import load_normal_expression
from .adapters.peptide_input import normalize_hla_allele

PEPTIDE_SAFETY_FIELDS = [
    'peptide_id','event_id','sample_id','event_type','mutation_source','peptide_consequence','gene',
    'peptide','hla_allele','mhc_class','matched_normal_status','normal_alt_reads','normal_vaf','tumor_only_flag',
    'reference_proteome_exact_match','reference_proteome_match_type','reference_match_gene','reference_match_protein',
    'reference_match_position','reference_match_length','normal_hla_ligand_exact_match','normal_ligand_tissue',
    'normal_ligand_hla','normal_ligand_source_protein','normal_ligand_class','normal_tissue_max_tpm',
    'normal_tissue_max_tissue','normal_hspc_tpm','critical_tissue_hit','critical_tissue_name','normal_junction_seen',
    'normal_junction_source','normal_junction_max_reads','normal_junction_tissue','wildtype_peptide','mt_binding_rank',
    'wt_binding_rank','mt_wt_fold_change','mutation_position_in_peptide','mutation_anchor_only','anchor_risk_status',
    'closest_self_peptide','closest_self_gene','closest_self_similarity','closest_self_hla_binding_rank',
    'closest_self_normal_expression_tpm','safety_tier','safety_status','safety_reason','safety_multiplier','review_required',
]

EVENT_SAFETY_FIELDS = [
    'event_id','sample_id','gene','event_type','mutation_source','normal_expression_status','normal_junction_status',
    'matched_normal_status','event_safety_status','event_safety_reason'
]


def _read_fasta_sequences(path: str | Path | None) -> dict[str, str]:
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f'Missing reference proteome FASTA: {p}')
    seqs: dict[str, str] = {}
    name = None; chunks: list[str] = []
    with p.open('r', encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(chunks).upper().replace('*','')
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            seqs[name] = ''.join(chunks).upper().replace('*','')
    return seqs


def build_reference_index(path: str | Path | None, lengths: tuple[int, ...] = (8,9,10,11)) -> dict[str, list[tuple[str,int]]]:
    seqs = _read_fasta_sequences(path)
    idx: dict[str, list[tuple[str,int]]] = {}
    for prot, seq in seqs.items():
        for L in lengths:
            if len(seq) < L:
                continue
            for i in range(0, len(seq)-L+1):
                pep = seq[i:i+L]
                if 'X' in pep:
                    continue
                idx.setdefault(pep, []).append((prot, i+1))
    return idx


def _normal_ligand_rows(path: str | Path | None) -> list[dict[str, str]]:
    if not path:
        return []
    p = Path(path)
    if not p.exists():
        return []
    out=[]
    for r in read_tsv(p):
        pep = first(r, ['peptide','Peptide','sequence','Sequence'], '').upper().strip()
        if not pep:
            continue
        out.append({
            'peptide': pep,
            'hla_allele': first(r, ['hla_allele','HLA','allele'], ''),
            'tissue': first(r, ['tissue','source_tissue','normal_tissue','source'], ''),
            'source_protein': first(r, ['source_protein','protein','gene','source_gene'], ''),
            'hla_class': first(r, ['hla_class','mhc_class','class'], ''),
        })
    return out


def _normal_junctions(path: str | Path | None) -> dict[str, dict[str, str]]:
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        return {}
    out={}
    for r in read_tsv(p):
        keys=[]
        for k in ('event_id','junction_id','gene_pair','fusion'):
            v = r.get(k)
            if v:
                keys.append(v)
        if r.get('gene1') and r.get('gene2'):
            keys.extend([f"{r.get('gene1')}::{r.get('gene2')}", f"{r.get('gene2')}::{r.get('gene1')}"])
        for k in keys:
            out[str(k)] = r
    return out


def _split_genes(gene: str) -> list[str]:
    raw = str(gene or '')
    for sep in ('::',';','|',','):
        raw = raw.replace(sep, ' ')
    return [x for x in raw.split() if x]


def _normal_expr_for_genes(gene_value: str, normal_expr: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    nt=0.0; nh=0.0; crit=False; found=False
    for g in _split_genes(gene_value):
        row = normal_expr.get(g) or {}
        if row:
            found=True
            nt=max(nt, to_float(row.get('normal_tissue_max_tpm'), 0.0))
            nh=max(nh, to_float(row.get('normal_hspc_tpm'), 0.0))
            crit=crit or to_float(row.get('critical_tissue_hit'), 0.0) > 0
    return {'found': found, 'normal_tissue_max_tpm': nt, 'normal_hspc_tpm': nh, 'critical_tissue_hit': crit}


def _find_normal_ligand(peptide: str, hla: str, ligs: list[dict[str, str]]) -> dict[str, str] | None:
    h_norm = normalize_hla_allele(hla)
    exact = [r for r in ligs if r['peptide'] == peptide]
    if not exact:
        return None
    for r in exact:
        if normalize_hla_allele(r.get('hla_allele','')) == h_norm and h_norm:
            return r
    return exact[0]


def _event_normal_support(event: Mapping[str, Any]) -> tuple[bool, str, str, str]:
    n_alt = to_float(event.get('normal_alt_count') or event.get('normal_alt_support'), 0.0)
    n_vaf = to_float(event.get('normal_vaf'), 0.0)
    if n_alt >= 3 and (n_vaf >= 0.02 or n_vaf == 0.0):
        return True, str(int(n_alt)), f'{n_vaf:.4f}', 'normal_alt_support'
    return False, str(int(n_alt)), f'{n_vaf:.4f}', 'no_clear_normal_support'


def _anchor_only(peptide: str, wt_peptide: str, hla: str) -> tuple[str, str, str]:
    if not peptide or not wt_peptide or len(peptide) != len(wt_peptide):
        return 'unknown', '', 'unassessed'
    diffs = [i+1 for i,(a,b) in enumerate(zip(peptide, wt_peptide)) if a != b]
    if not diffs:
        return 'no', '', 'no_difference'
    L = len(peptide)
    # Approximate MHC-I anchor positions; class/allele-specific anchor maps can replace this later.
    anchors = {2, L}
    only = all(d in anchors for d in diffs)
    return 'yes' if only else 'no', ','.join(map(str,diffs)), 'anchor_only' if only else 'tcr_facing_change_present'


def safety_multiplier_for(status: str) -> float:
    s = str(status or '').upper()
    if s in {'SAFETY_PASS','PASS'}:
        return 1.0
    if s in {'SAFETY_REVIEW','REVIEW','CAUTION'}:
        return 0.55
    if s in {'SAFETY_HIGH_RISK','HIGH_RISK'}:
        return 0.20
    if s in {'SAFETY_REJECT','REJECT','FAIL'}:
        return 0.0
    return 0.60


def build_peptide_safety_gate(
    *,
    raw_events: str | Path,
    raw_peptides: str | Path,
    out_peptide_safety: str | Path,
    out_event_safety: str | Path | None = None,
    profile: Mapping[str, Any] | None = None,
    normal_expression: str | Path | None = None,
    normal_hla_ligands: str | Path | None = None,
    reference_proteome: str | Path | None = None,
    normal_junctions: str | Path | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    events = {r.get('event_id',''): r for r in read_tsv(raw_events)}
    peptides = read_tsv(raw_peptides)
    norm_expr = load_normal_expression(normal_expression)
    ligs = _normal_ligand_rows(normal_hla_ligands)
    ref_idx = build_reference_index(reference_proteome)
    junctions = _normal_junctions(normal_junctions)
    safety_cfg = (profile or {}).get('safety', {})
    ref_reject = str(safety_cfg.get('reference_proteome_exact_match_action','FAIL')).upper() in {'FAIL','REJECT'}
    normal_lig_action = str(safety_cfg.get('normal_hla_ligand_overlap_action','CAUTION')).upper()
    rows=[]; event_rows=[]
    for p in peptides:
        e = events.get(p.get('event_id',''), {})
        peptide = (p.get('peptide') or '').upper().strip()
        hla = p.get('hla_allele','')
        reasons=[]; status='SAFETY_PASS'; review='no'
        normal_support, normal_alt, normal_vaf, normal_status = _event_normal_support(e)
        if normal_support:
            status='SAFETY_REJECT'; reasons.append('matched_normal_support')
        ref_hits = ref_idx.get(peptide, []) if peptide else []
        if ref_hits and ref_reject:
            status='SAFETY_REJECT'; reasons.append('reference_proteome_exact_match')
        lig = _find_normal_ligand(peptide, hla, ligs) if peptide else None
        if lig:
            if normal_lig_action in {'FAIL','REJECT'} or normalize_hla_allele(lig.get('hla_allele','')) == normalize_hla_allele(hla):
                if status != 'SAFETY_REJECT':
                    status='SAFETY_HIGH_RISK'
                reasons.append('normal_hla_ligand_exact_match')
            else:
                if status == 'SAFETY_PASS': status='SAFETY_REVIEW'
                reasons.append('normal_hla_ligand_overlap')
        event_key = e.get('event_id','')
        gene_pair = e.get('gene','')
        jrow = junctions.get(event_key) or junctions.get(gene_pair)
        is_junction_type = str(p.get('peptide_consequence') or p.get('event_type') or '').lower() in {'fusion','splice_junction','sv_junction','exon_deletion_junction'} or p.get('crosses_junction') == 'yes'
        if jrow and is_junction_type:
            max_reads = to_float(first(jrow, ['normal_reads','junction_reads','reads'], '0'), 0.0)
            if max_reads >= float(safety_cfg.get('normal_junction_reject_reads', 3)):
                status='SAFETY_REJECT'; reasons.append('normal_junction_seen')
            else:
                if status == 'SAFETY_PASS': status='SAFETY_REVIEW'
                reasons.append('low_level_normal_junction_seen')
        expr = _normal_expr_for_genes(p.get('gene') or e.get('gene',''), norm_expr)
        wt_rank = to_float(p.get('wildtype_binding_rank') or p.get('netmhcpan_wt_rank_el') or p.get('netmhcpan_wt_rank_ba'), 99.0)
        mt_rank = to_float(p.get('netmhcpan_el_rank') or p.get('netmhcpan_mt_rank_el') or p.get('binding_rank'), 99.0)
        anchor_only, mut_pos, anchor_status = _anchor_only(peptide, p.get('wildtype_peptide',''), hla)
        if anchor_only == 'yes' and wt_rank <= float(safety_cfg.get('wildtype_strong_binding_rank', 0.5)):
            if status == 'SAFETY_PASS': status='SAFETY_REVIEW'
            reasons.append('anchor_only_mutation_with_wt_binding')
        sim = to_float(p.get('self_similarity_score'), 0.0)
        if sim >= float(safety_cfg.get('self_similarity_high_risk', safety_cfg.get('self_similarity_caution', 0.85))) and expr.get('critical_tissue_hit'):
            if status not in {'SAFETY_REJECT'}: status='SAFETY_HIGH_RISK'
            reasons.append('high_self_similarity_critical_tissue')
        if status in {'SAFETY_REVIEW','SAFETY_HIGH_RISK'}:
            review='yes'
        if not reasons:
            reasons.append('no_major_signal')
        row={
            'peptide_id': p.get('peptide_id',''), 'event_id': p.get('event_id',''), 'sample_id': p.get('sample_id',''),
            'event_type': p.get('event_type',''), 'mutation_source': p.get('mutation_source',''),
            'peptide_consequence': p.get('peptide_consequence',''), 'gene': p.get('gene',''), 'peptide': peptide,
            'hla_allele': hla, 'mhc_class': p.get('mhc_class',''), 'matched_normal_status': normal_status,
            'normal_alt_reads': normal_alt, 'normal_vaf': normal_vaf, 'tumor_only_flag': 'yes' if not e else 'no',
            'reference_proteome_exact_match': 'yes' if ref_hits else 'no' if ref_idx else 'not_assessed',
            'reference_proteome_match_type': 'exact' if ref_hits else '',
            'reference_match_gene': '', 'reference_match_protein': ref_hits[0][0] if ref_hits else '',
            'reference_match_position': str(ref_hits[0][1]) if ref_hits else '', 'reference_match_length': str(len(peptide)) if ref_hits else '',
            'normal_hla_ligand_exact_match': 'yes' if lig else 'no',
            'normal_ligand_tissue': (lig or {}).get('tissue',''), 'normal_ligand_hla': (lig or {}).get('hla_allele',''),
            'normal_ligand_source_protein': (lig or {}).get('source_protein',''), 'normal_ligand_class': (lig or {}).get('hla_class',''),
            'normal_tissue_max_tpm': f"{expr.get('normal_tissue_max_tpm',0.0):.4f}",
            'normal_tissue_max_tissue': '', 'normal_hspc_tpm': f"{expr.get('normal_hspc_tpm',0.0):.4f}",
            'critical_tissue_hit': 'yes' if expr.get('critical_tissue_hit') else 'no', 'critical_tissue_name': '',
            'normal_junction_seen': 'yes' if jrow else 'no', 'normal_junction_source': first(jrow or {}, ['source','dataset'], ''),
            'normal_junction_max_reads': first(jrow or {}, ['normal_reads','junction_reads','reads'], ''),
            'normal_junction_tissue': first(jrow or {}, ['tissue','normal_tissue'], ''),
            'wildtype_peptide': p.get('wildtype_peptide',''), 'mt_binding_rank': f'{mt_rank:.4f}', 'wt_binding_rank': f'{wt_rank:.4f}',
            'mt_wt_fold_change': '', 'mutation_position_in_peptide': mut_pos, 'mutation_anchor_only': anchor_only,
            'anchor_risk_status': anchor_status, 'closest_self_peptide': '', 'closest_self_gene': '',
            'closest_self_similarity': p.get('self_similarity_score',''), 'closest_self_hla_binding_rank': '',
            'closest_self_normal_expression_tpm': '', 'safety_tier': status, 'safety_status': 'FAIL' if status=='SAFETY_REJECT' else ('CAUTION' if status in {'SAFETY_REVIEW','SAFETY_HIGH_RISK'} else 'PASS'),
            'safety_reason': ';'.join(dict.fromkeys(reasons)), 'safety_multiplier': f'{safety_multiplier_for(status):.4f}',
            'review_required': review,
        }
        rows.append(row)
    # event-level rollup
    by_event: dict[str, list[dict[str, Any]]] = {}
    for r in rows:
        by_event.setdefault(r['event_id'], []).append(r)
    for eid, rs in by_event.items():
        e = events.get(eid, {})
        worst = 'PASS'
        reasons=[]
        if any(r['safety_status']=='FAIL' for r in rs): worst='FAIL'
        elif any(r['safety_tier']=='SAFETY_HIGH_RISK' for r in rs): worst='CAUTION'
        elif any(r['safety_status']=='CAUTION' for r in rs): worst='CAUTION'
        for r in rs:
            if r['safety_reason'] != 'no_major_signal': reasons.append(r['safety_reason'])
        event_rows.append({
            'event_id': eid, 'sample_id': e.get('sample_id',''), 'gene': e.get('gene',''), 'event_type': e.get('event_type',''),
            'mutation_source': e.get('mutation_source',''), 'normal_expression_status': 'assessed' if norm_expr else 'not_assessed',
            'normal_junction_status': 'assessed' if junctions else 'not_assessed', 'matched_normal_status': '',
            'event_safety_status': worst, 'event_safety_reason': ';'.join(dict.fromkeys(reasons)) if reasons else 'no_major_signal'
        })
    write_tsv(out_peptide_safety, rows, PEPTIDE_SAFETY_FIELDS)
    if out_event_safety:
        write_tsv(out_event_safety, event_rows, EVENT_SAFETY_FIELDS)
    return rows, event_rows


def load_peptide_safety(path: str | Path | None) -> dict[str, dict[str, str]]:
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        return {}
    return {r.get('peptide_id',''): r for r in read_tsv(p) if r.get('peptide_id')}
