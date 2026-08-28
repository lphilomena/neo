#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def load_tx2gene(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with path.open('r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                parts = line.rstrip('\n').split(',')
            if len(parts) >= 2:
                mapping[parts[0].split('.')[0]] = parts[1]
                mapping[parts[0]] = parts[1]
    return mapping


def main() -> int:
    ap = argparse.ArgumentParser(description='Convert Salmon quant.sf transcript TPM to gene TPM')
    ap.add_argument('--quant-sf', required=True)
    ap.add_argument('--tx2gene', required=True, help='TSV/CSV: transcript_id gene_id [gene_name]')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    tx2gene = load_tx2gene(Path(args.tx2gene))
    gene_tpm: dict[str, float] = defaultdict(float)
    missing = 0
    with Path(args.quant_sf).open('r', encoding='utf-8', errors='replace') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            tx = (row.get('Name') or row.get('target_id') or '').strip()
            if not tx:
                continue
            gene = tx2gene.get(tx) or tx2gene.get(tx.split('.')[0])
            if not gene:
                missing += 1
                continue
            try:
                tpm = float(row.get('TPM') or row.get('tpm') or 0.0)
            except Exception:
                tpm = 0.0
            gene_tpm[gene] += tpm
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open('w', encoding='utf-8', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=['gene_id', 'tpm'], delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for gene, tpm in sorted(gene_tpm.items(), key=lambda x: (-x[1], x[0])):
            writer.writerow({'gene_id': gene, 'tpm': f'{tpm:.6f}'})
    meta = out.with_suffix(out.suffix + '.summary.txt')
    meta.write_text(f'genes={len(gene_tpm)}\nmissing_transcripts={missing}\n', encoding='utf-8')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
