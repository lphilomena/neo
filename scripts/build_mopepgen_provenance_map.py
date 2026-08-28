#!/usr/bin/env python3
"""Build an exact moPepGen→NeoAg provenance map without gene-level fallback.

The input is a moPepGen peptide table (or a reviewed table derived from it).  A
second exact map may supply canonical splice_event_id/junction_id per FASTA
header or variant_id.  Records that cannot resolve uniquely are retained with
resolution_status=UNRESOLVED and cause --strict to fail.
"""
from __future__ import annotations
import argparse, csv, hashlib, json
from collections import defaultdict
from pathlib import Path

ALIASES = {
    "header": ("mopepgen_header", "fasta_header", "header", "record_id", "sequence_name"),
    "variant": ("variant_id", "mopepgen_variant_id", "variant_ids"),
    "peptide": ("peptide_sequence", "sequence", "peptide"),
    "event": ("splice_event_id", "event_id"),
    "junction": ("junction_id", "canonical_junction_id", "junc_id"),
    "gene": ("gene", "gene_name", "gene_symbol"),
    "transcript": ("transcript_id", "reference_transcript_id"),
}
def get(row, key):
    lower={str(k).lower():str(v).strip() for k,v in row.items()}
    for name in ALIASES[key]:
        if lower.get(name.lower()): return lower[name.lower()]
    return ""
def read(path):
    with Path(path).open(newline='',encoding='utf-8-sig') as f:
        sample=f.read(4096); f.seek(0)
        try: dialect=csv.Sniffer().sniff(sample,delimiters='\t,')
        except csv.Error: dialect=csv.excel_tab
        return list(csv.DictReader(f,dialect=dialect))
def split_tokens(value):
    return [x.strip() for x in value.replace(',', ';').split(';') if x.strip()]
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--peptide-table',required=True)
    ap.add_argument('--exact-map',action='append',default=[])
    ap.add_argument('--output',required=True)
    ap.add_argument('--strict',action='store_true')
    args=ap.parse_args()
    by_header=defaultdict(list); by_variant=defaultdict(list)
    for path in args.exact_map:
        for row in read(path):
            h=get(row,'header'); variants=split_tokens(get(row,'variant'))
            if h: by_header[h].append(row)
            for v in variants: by_variant[v].append(row)
    out=[]; unresolved=0
    for i,row in enumerate(read(args.peptide_table),1):
        header=get(row,'header'); variants=split_tokens(get(row,'variant'))
        candidates=list(by_header.get(header,[]))
        for v in variants: candidates.extend(by_variant.get(v,[]))
        unique={json.dumps(r,sort_keys=True):r for r in candidates}.values()
        events={get(r,'event') for r in unique if get(r,'event')}
        junctions={get(r,'junction') for r in unique if get(r,'junction')}
        genes={get(r,'gene') for r in unique if get(r,'gene')}
        transcripts={get(r,'transcript') for r in unique if get(r,'transcript')}
        status='RESOLVED_EXACT' if len(events)==1 and len(junctions)>=1 else 'UNRESOLVED'
        if status=='UNRESOLVED': unresolved+=1
        out.append({
          'mopepgen_header':header,'variant_id':';'.join(variants),
          'splice_event_id':next(iter(events)) if len(events)==1 else '',
          'junction_id':';'.join(sorted(junctions)),
          'gene':next(iter(genes)) if len(genes)==1 else get(row,'gene'),
          'transcript_id':next(iter(transcripts)) if len(transcripts)==1 else get(row,'transcript'),
          'peptide_sequence':get(row,'peptide'),'crosses_junction':row.get('crosses_junction','UNASSESSED'),
          'contains_novel_aa':row.get('contains_novel_aa','true'),
          'source_row_number':str(i),'resolution_status':status,
        })
    fields=list(out[0]) if out else ['mopepgen_header','variant_id','splice_event_id','junction_id','gene','transcript_id','peptide_sequence','crosses_junction','contains_novel_aa','source_row_number','resolution_status']
    path=Path(args.output); path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(out)
    h=hashlib.sha256(path.read_bytes()).hexdigest()
    print(json.dumps({'output':str(path),'sha256':h,'rows':len(out),'unresolved':unresolved},indent=2))
    if args.strict and unresolved: raise SystemExit(f'{unresolved} moPepGen records lack exact event+junction provenance')
if __name__=='__main__': main()
