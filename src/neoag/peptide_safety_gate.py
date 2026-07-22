
"""Peptide-level safety evidence resolver.

The gate is deliberately conservative and explainable. It writes a sidecar
`peptide_safety.tsv` that score can consume, instead of hiding safety logic
inside a single scalar ranking score.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from .utils import read_tsv, write_tsv, to_float, first
from .safety import apply_event_safety, load_normal_expression
from .safety import combine_safety_reasons, worst_safety_status
from .adapters.peptide_input import normalize_hla_allele

PEPTIDE_SAFETY_FIELDS = [
    'peptide_id','event_id','sample_id','event_type','mutation_source','peptide_consequence','gene',
    'peptide','hla_allele','mhc_class','matched_normal_status','normal_alt_reads','normal_vaf','tumor_only_flag',
    'reference_proteome_exact_match','reference_proteome_match_type','reference_match_gene','reference_match_protein',
    'reference_match_position','reference_match_length','normal_hla_ligand_exact_match','normal_ligand_tissue',
    'normal_ligand_hla','normal_ligand_source_protein','normal_ligand_class','normal_tissue_max_tpm',
    'normal_tissue_max_tissue','critical_tissue_max_tpm','normal_hspc_tpm','normal_hspc_unit','critical_tissue_hit','critical_tissue_name','normal_junction_seen',
    'normal_junction_source','normal_junction_max_reads','normal_junction_tissue','wildtype_peptide','mt_binding_rank',
    'wt_binding_rank','mt_wt_fold_change','mutation_position_in_peptide','mutation_anchor_only','anchor_risk_status',
    'closest_self_peptide','closest_self_gene','closest_self_similarity','closest_self_hla_binding_rank',
    'closest_self_normal_expression_tpm','safety_tier','safety_status','safety_reason','safety_multiplier','review_required',
    'normal_expression_status','normal_hspc_status','reference_proteome_status','normal_ligandome_status',
    'anchor_assessment_status','normal_junction_assessment_status','safety_evidence_completeness',
    'safety_missing_layers','safety_priority_cap',
]

EVENT_SAFETY_FIELDS = [
    'event_id','sample_id','gene','event_type','mutation_source','normal_expression_status','normal_junction_status',
    'matched_normal_status','event_safety_status','event_safety_reason',
    'normal_hspc_status','reference_proteome_status','normal_ligandome_status','anchor_assessment_status',
    'normal_tissue_max_tpm','normal_tissue_max_tissue','critical_tissue_max_tpm','critical_tissue_name',
    'normal_hspc_tpm','normal_hspc_unit','safety_evidence_completeness','safety_missing_layers'
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


def build_reference_index(
    path: str | Path | None,
    lengths: tuple[int, ...] = (8,9,10,11),
    target_peptides: set[str] | None = None,
) -> dict[str, list[tuple[str,int]]]:
    seqs = _read_fasta_sequences(path)
    idx: dict[str, list[tuple[str,int]]] = {}
    for prot, seq in seqs.items():
        for L in lengths:
            if len(seq) < L:
                continue
            for i in range(0, len(seq)-L+1):
                pep = seq[i:i+L]
                if 'X' in pep or (target_peptides is not None and pep not in target_peptides):
                    continue
                idx.setdefault(pep, []).append((prot, i+1))
    return idx


def _normal_ligand_rows(
    path: str | Path | None,
    target_peptides: set[str] | None = None,
) -> dict[str, list[dict[str, str]]]:
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        return {}
    out: dict[str, list[dict[str, str]]] = {}
    for r in read_tsv(p):
        pep = first(r, ['peptide','Peptide','sequence','Sequence'], '').upper().strip()
        if not pep or (target_peptides is not None and pep not in target_peptides):
            continue
        out.setdefault(pep, []).append({
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
        normalized=[]
        for key in keys:
            value = str(key).strip()
            if not value:
                continue
            normalized.extend([value, value.replace('--', '::')])
            if '::' in value.replace('--', '::'):
                left, right = value.replace('--', '::').split('::', 1)
                normalized.append(f'{right}::{left}')
        for key in normalized:
            out[key] = r
    return out


def _normal_junction_scopes(path: str | Path | None) -> set[str]:
    if not path or not Path(path).is_file():
        return set()
    scopes = {
        str(row.get('junction_class') or '').strip()
        for row in read_tsv(path)
        if str(row.get('junction_class') or '').strip()
    }
    # Backward compatibility for user-provided files predating scoped refs.
    return scopes or {'all_junctions'}


def _normal_junction_kind(
    peptide: Mapping[str, Any],
    event: Mapping[str, Any],
) -> str:
    consequence = str(peptide.get('peptide_consequence') or '').lower()
    event_type = str(peptide.get('event_type') or event.get('event_type') or '').lower()
    mutation_source = str(peptide.get('mutation_source') or event.get('mutation_source') or '').lower()
    if consequence == 'exon_deletion_junction':
        return 'not_applicable'
    # VEP can label an SNV/InDel as splice_region_variant while the peptide
    # catalog labels its protein product as splice_junction. Without an
    # observed RNA exon-exon junction this is still a sequence variant, not a
    # junction that can be looked up in GTEx.
    if event_type in {'snv', 'indel'} and mutation_source in {'snv', 'indel', 'wes', 'wgs'}:
        return 'not_applicable'
    if consequence == 'splice_junction' or event_type in {'splice', 'splice_junction'}:
        return 'splice'
    if consequence in {'fusion', 'sv_junction'} or event_type in {'fusion', 'sv_junction'}:
        return 'fusion'
    # crosses_junction also marks some variant-generated InDel peptides. Those
    # are not normal RNA splice/fusion junctions and must not require this layer.
    if peptide.get('crosses_junction') == 'yes' and mutation_source in {'sv', 'fusion'}:
        return 'fusion'
    if peptide.get('crosses_junction') == 'yes' and mutation_source == 'splice':
        return 'splice'
    return 'not_applicable'


def _candidate_junction_keys(
    peptide: Mapping[str, Any],
    event: Mapping[str, Any],
) -> list[str]:
    keys: list[str] = []
    for row in (peptide, event):
        for field in ('junction_id', 'event_id', 'event_name', 'gene_pair', 'fusion'):
            value = str(row.get(field) or '').strip()
            if value and value not in keys:
                keys.append(value)
        chrom = str(row.get('chrom') or '').strip()
        start = str(row.get('junction_start') or row.get('start') or '').strip()
        end = str(row.get('junction_end') or row.get('end') or '').strip()
        strand = str(row.get('strand') or '').strip()
        if chrom and start and end:
            coordinate = f'{chrom}:{start}-{end}'
            for value in (coordinate, f'{coordinate}:{strand}' if strand else ''):
                if value and value not in keys:
                    keys.append(value)
    return keys


def _split_genes(gene: str) -> list[str]:
    raw = str(gene or '')
    for sep in ('::',';','|',','):
        raw = raw.replace(sep, ' ')
    return [x for x in raw.split() if x]


def _normal_expr_for_genes(gene_value: str, normal_expr: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    nt=0.0; nh=0.0; crit=False; found=False; critical_tpm=0.0
    nt_tissue=''; critical_tissue=''; hspc_unit=''; expr_assessed=True; hspc_assessed=True
    genes = _split_genes(gene_value)
    for g in genes:
        row = normal_expr.get(g) or {}
        if row:
            found=True
            row_nt=to_float(row.get('normal_tissue_max_tpm'), 0.0)
            if row_nt >= nt: nt=row_nt; nt_tissue=row.get('normal_tissue_max_tissue','')
            row_hspc=to_float(row.get('normal_hspc_tpm'), 0.0)
            if row_hspc >= nh: nh=row_hspc; hspc_unit=row.get('normal_hspc_unit','TPM')
            row_critical=to_float(row.get('critical_tissue_max_tpm'), row_nt)
            if row_critical >= critical_tpm: critical_tpm=row_critical; critical_tissue=row.get('critical_tissue_name','')
            crit=crit or to_float(row.get('critical_tissue_hit'), 0.0) > 0
            expr_assessed=expr_assessed and row.get('normal_expression_status','ASSESSED') == 'ASSESSED'
            hspc_assessed=hspc_assessed and row.get('normal_hspc_status','ASSESSED') == 'ASSESSED'
        else:
            expr_assessed=False; hspc_assessed=False
    return {
        'found': found, 'normal_tissue_max_tpm': nt, 'normal_tissue_max_tissue': nt_tissue,
        'critical_tissue_max_tpm': critical_tpm, 'critical_tissue_name': critical_tissue,
        'normal_hspc_tpm': nh, 'normal_hspc_unit': hspc_unit or 'TPM', 'critical_tissue_hit': crit,
        'normal_expression_status': 'ASSESSED' if found and expr_assessed else 'UNASSESSED',
        'normal_hspc_status': 'ASSESSED' if found and hspc_assessed else 'UNASSESSED',
    }


def _find_normal_ligand(
    peptide: str,
    hla: str,
    ligs: Mapping[str, list[dict[str, str]]],
) -> dict[str, str] | None:
    h_norm = normalize_hla_allele(hla)
    exact = ligs.get(peptide, [])
    if not exact:
        return None
    for r in exact:
        if normalize_hla_allele(r.get('hla_allele','')) == h_norm and h_norm:
            return r
    return exact[0]


def _event_normal_support(event: Mapping[str, Any]) -> tuple[bool, str, str, str]:
    has_evidence = any(str(event.get(key, '')).strip() for key in ('normal_alt_count','normal_alt_support','normal_vaf','normal_depth'))
    if not has_evidence:
        return False, '', '', 'not_assessed'
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


def safety_multiplier_for(status: str, profile: Mapping[str, Any] | None = None) -> float:
    cfg = (profile or {}).get('safety', {})
    s = str(status or '').upper()
    if s in {'SAFETY_PASS','PASS'}:
        return 1.0
    if s in {'SAFETY_REVIEW','REVIEW','CAUTION'}:
        return float(cfg.get('caution_multiplier', 0.55))
    if s in {'SAFETY_PARTIAL','PARTIAL'}:
        return float(cfg.get('partial_multiplier', 0.75))
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
    target_peptides = {(row.get('peptide') or '').upper().strip() for row in peptides if row.get('peptide')}
    ligs = _normal_ligand_rows(normal_hla_ligands, target_peptides=target_peptides)
    ref_idx = build_reference_index(reference_proteome, target_peptides=target_peptides)
    junctions = _normal_junctions(normal_junctions)
    junction_scopes = _normal_junction_scopes(normal_junctions)
    ref_assessed = bool(reference_proteome and Path(reference_proteome).is_file())
    ligand_assessed = bool(normal_hla_ligands and Path(normal_hla_ligands).is_file())
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
        gene_pair = e.get('gene','')
        candidate_junction_keys = _candidate_junction_keys(p, e)
        jrow = next((junctions.get(key) for key in candidate_junction_keys if junctions.get(key)), None)
        if not jrow:
            jrow = junctions.get(gene_pair)
        junction_kind = _normal_junction_kind(p, e)
        is_junction_type = junction_kind != 'not_applicable'
        if junction_kind == 'fusion':
            junction_assessed = bool({'normal_recurrent_fusion', 'all_junctions'} & junction_scopes)
        elif junction_kind == 'splice':
            junction_assessed = bool({'normal_splice_junction', 'all_junctions'} & junction_scopes)
        else:
            junction_assessed = False
        if jrow and is_junction_type:
            max_reads = to_float(first(jrow, ['normal_reads','normal_sample_count','junction_reads','reads'], '0'), 0.0)
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
        expression_status = expr.get('normal_expression_status', 'UNASSESSED')
        hspc_status = expr.get('normal_hspc_status', 'UNASSESSED')
        reference_status = 'ASSESSED' if ref_assessed else 'UNASSESSED'
        ligandome_status = 'ASSESSED' if ligand_assessed else 'UNASSESSED'
        junction_type = is_junction_type
        junction_status = ('ASSESSED' if junction_assessed else 'UNASSESSED') if junction_type else 'NOT_APPLICABLE'
        if anchor_status != 'unassessed':
            anchor_assessment = 'ASSESSED'
        elif junction_type or str(p.get('peptide_consequence','')).lower() in {'frameshift','insertion'}:
            anchor_assessment = 'NOT_APPLICABLE'
        else:
            anchor_assessment = 'UNASSESSED'
        required_layers = {
            'normal_expression': expression_status,
            'normal_hspc': hspc_status,
            'reference_proteome': reference_status,
            'normal_ligandome': ligandome_status,
            'anchor_only': anchor_assessment,
            'normal_junction': junction_status,
        }
        assessed = [name for name, layer_status in required_layers.items() if layer_status != 'NOT_APPLICABLE']
        missing = [name for name in assessed if required_layers[name] != 'ASSESSED']
        completeness = (len(assessed) - len(missing)) / max(len(assessed), 1)
        if missing:
            status = worst_safety_status(status, 'SAFETY_PARTIAL')
            reasons.append('safety_evidence_incomplete')
        if status in {'SAFETY_REVIEW','SAFETY_HIGH_RISK','SAFETY_PARTIAL'}:
            review='yes'
        if not reasons:
            reasons.append('no_major_signal')
        row={
            'peptide_id': p.get('peptide_id',''), 'event_id': p.get('event_id',''), 'sample_id': p.get('sample_id',''),
            'event_type': p.get('event_type',''), 'mutation_source': p.get('mutation_source',''),
            'peptide_consequence': p.get('peptide_consequence',''), 'gene': p.get('gene',''), 'peptide': peptide,
            'hla_allele': hla, 'mhc_class': p.get('mhc_class',''), 'matched_normal_status': normal_status,
            'normal_alt_reads': normal_alt, 'normal_vaf': normal_vaf, 'tumor_only_flag': 'yes' if not e else 'no',
            'reference_proteome_exact_match': 'yes' if ref_hits else 'no' if ref_assessed else 'not_assessed',
            'reference_proteome_match_type': 'exact' if ref_hits else '',
            'reference_match_gene': '', 'reference_match_protein': ref_hits[0][0] if ref_hits else '',
            'reference_match_position': str(ref_hits[0][1]) if ref_hits else '', 'reference_match_length': str(len(peptide)) if ref_hits else '',
            'normal_hla_ligand_exact_match': 'yes' if lig else 'no' if ligand_assessed else 'not_assessed',
            'normal_ligand_tissue': (lig or {}).get('tissue',''), 'normal_ligand_hla': (lig or {}).get('hla_allele',''),
            'normal_ligand_source_protein': (lig or {}).get('source_protein',''), 'normal_ligand_class': (lig or {}).get('hla_class',''),
            'normal_tissue_max_tpm': f"{expr.get('normal_tissue_max_tpm',0.0):.4f}" if expr.get('found') else '',
            'normal_tissue_max_tissue': expr.get('normal_tissue_max_tissue',''),
            'critical_tissue_max_tpm': f"{expr.get('critical_tissue_max_tpm',0.0):.4f}" if expr.get('found') else '',
            'normal_hspc_tpm': f"{expr.get('normal_hspc_tpm',0.0):.4f}" if expr.get('found') else '',
            'normal_hspc_unit': expr.get('normal_hspc_unit','') if expr.get('found') else '',
            'critical_tissue_hit': ('yes' if expr.get('critical_tissue_hit') else 'no') if expr.get('found') else 'not_assessed',
            'critical_tissue_name': expr.get('critical_tissue_name',''),
            'normal_junction_seen': (
                'yes' if jrow and is_junction_type else
                'no' if junction_assessed else
                'not_applicable' if junction_kind == 'not_applicable' else
                'not_assessed'
            ),
            'normal_junction_source': first(jrow or {}, ['source','dataset'], '') if is_junction_type else '',
            'normal_junction_max_reads': first(jrow or {}, ['normal_reads','normal_sample_count','junction_reads','reads'], ''),
            'normal_junction_tissue': first(jrow or {}, ['tissue','normal_tissue'], ''),
            'wildtype_peptide': p.get('wildtype_peptide',''), 'mt_binding_rank': f'{mt_rank:.4f}', 'wt_binding_rank': f'{wt_rank:.4f}',
            'mt_wt_fold_change': '', 'mutation_position_in_peptide': mut_pos, 'mutation_anchor_only': anchor_only,
            'anchor_risk_status': anchor_status, 'closest_self_peptide': '', 'closest_self_gene': '',
            'closest_self_similarity': p.get('self_similarity_score',''), 'closest_self_hla_binding_rank': '',
            'closest_self_normal_expression_tpm': '', 'safety_tier': status, 'safety_status': 'FAIL' if status=='SAFETY_REJECT' else ('CAUTION' if status in {'SAFETY_REVIEW','SAFETY_HIGH_RISK'} else ('SAFETY_PARTIAL' if status=='SAFETY_PARTIAL' else 'PASS')),
            'safety_reason': ';'.join(dict.fromkeys(reasons)), 'safety_multiplier': f'{safety_multiplier_for(status, profile):.4f}',
            'review_required': review,
            'normal_expression_status': expression_status, 'normal_hspc_status': hspc_status,
            'reference_proteome_status': reference_status, 'normal_ligandome_status': ligandome_status,
            'anchor_assessment_status': anchor_assessment, 'normal_junction_assessment_status': junction_status,
            'safety_evidence_completeness': f'{completeness:.4f}', 'safety_missing_layers': ';'.join(missing),
            'safety_priority_cap': str(safety_cfg.get('partial_priority_cap','C_CAUTION')) if missing else '',
        }
        rows.append(row)
    # event-level rollup
    by_event: dict[str, list[dict[str, Any]]] = {}
    for r in rows:
        by_event.setdefault(r['event_id'], []).append(r)
    for eid, rs in by_event.items():
        e = events.get(eid, {})
        event_assessment = apply_event_safety(dict(e), profile or {}, norm_expr)
        worst = event_assessment.get('safety_status', 'SAFETY_PARTIAL')
        reasons = [event_assessment.get('safety_reason', '')]
        matched_normal = any(r.get('matched_normal_status') == 'normal_alt_support' for r in rs)
        if matched_normal:
            worst = worst_safety_status(worst, 'FAIL')
            reasons.append('matched_normal_support')
        junction_rows = [r for r in rs if r.get('normal_junction_seen') == 'yes']
        if junction_rows:
            max_normal_junction_reads = max(to_float(r.get('normal_junction_max_reads'), 0.0) for r in junction_rows)
            if max_normal_junction_reads >= float(safety_cfg.get('normal_junction_reject_reads', 3)):
                worst = worst_safety_status(worst, 'FAIL')
                reasons.append('normal_junction_seen')
            else:
                worst = worst_safety_status(worst, 'CAUTION')
                reasons.append('low_level_normal_junction_seen')
        missing_layers = sorted({layer for r in rs for layer in r.get('safety_missing_layers','').split(';') if layer})
        event_rows.append({
            'event_id': eid, 'sample_id': e.get('sample_id',''), 'gene': e.get('gene',''), 'event_type': e.get('event_type',''),
            'mutation_source': e.get('mutation_source',''), 'normal_expression_status': rs[0].get('normal_expression_status','UNASSESSED'),
            'normal_junction_status': rs[0].get('normal_junction_assessment_status','UNASSESSED'),
            'matched_normal_status': 'normal_alt_support' if matched_normal else rs[0].get('matched_normal_status',''),
            'normal_hspc_status': rs[0].get('normal_hspc_status','UNASSESSED'),
            'reference_proteome_status': rs[0].get('reference_proteome_status','UNASSESSED'),
            'normal_ligandome_status': rs[0].get('normal_ligandome_status','UNASSESSED'),
            'anchor_assessment_status': 'PARTIAL' if len({r.get('anchor_assessment_status') for r in rs}) > 1 else rs[0].get('anchor_assessment_status','UNASSESSED'),
            'normal_tissue_max_tpm': event_assessment.get('normal_tissue_max_tpm',''),
            'normal_tissue_max_tissue': event_assessment.get('normal_tissue_max_tissue',''),
            'critical_tissue_max_tpm': event_assessment.get('critical_tissue_max_tpm',''),
            'critical_tissue_name': event_assessment.get('critical_tissue_name',''),
            'normal_hspc_tpm': event_assessment.get('normal_hspc_tpm',''),
            'normal_hspc_unit': event_assessment.get('normal_hspc_unit',''),
            'safety_evidence_completeness': f"{min(to_float(r.get('safety_evidence_completeness'),0.0) for r in rs):.4f}",
            'safety_missing_layers': ';'.join(missing_layers),
            'event_safety_status': worst, 'event_safety_reason': combine_safety_reasons(*reasons, 'safety_evidence_incomplete' if missing_layers else '')
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
