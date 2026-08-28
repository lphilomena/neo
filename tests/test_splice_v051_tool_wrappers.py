from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _write(path: Path, text: str, executable: bool = False) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    if executable:
        path.chmod(0o755)
    return path


def _run(*args: str, env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    merged = os.environ.copy()
    if env:
        merged.update(env)
    return subprocess.run(args, cwd=ROOT, env=merged, text=True, capture_output=True, check=True)


def test_easyquant_wrapper_emits_list_and_manifest(tmp_path: Path):
    query = _write(tmp_path / "queries.tsv", "name\tsequence\tposition\nSEQ1\tAAAACCCCGGGG\t5\n")
    bam = _write(tmp_path / "rna.bam", "fixture\n")
    fake = _write(
        tmp_path / "bin" / "bp_quant",
        """#!/usr/bin/env bash
set -euo pipefail
if [[ ${1:-} == --version ]]; then echo 'bp-quant 0.test'; exit 0; fi
out=''
while [[ $# -gt 0 ]]; do [[ $1 == -o ]] && { out=$2; shift 2; continue; }; shift; done
mkdir -p "$out"
printf 'name\\tpos\\tjunc\\tspan\\tanch\\nSEQ1\\t5\\t7\\t2\\t12\\n' > "$out/quantification.tsv"
""",
        executable=True,
    )
    out = tmp_path / "easy"
    cp = _run(
        "bash", "scripts/run_easyquant_sample.sh", "--query-table", str(query),
        "--bam", str(bam), "--outdir", str(out), "--bp-quant-bin", str(fake),
    )
    list_path = Path(cp.stdout.strip().splitlines()[-1])
    assert list_path.is_file()
    assert "quantification.tsv" in list_path.read_text()
    manifest = json.loads((out / "easyquant_run_manifest.json").read_text())
    assert manifest["tool"] == "EasyQuant/bp-quant"
    assert manifest["quantification"]


def test_k4neo_wrapper_requires_license_and_emits_categorized_lists(tmp_path: Path):
    query = _write(tmp_path / "k4.tsv", "cts_id\tcts_seq\nSEQ1\tAAAACCCCGGGG\n")
    db = tmp_path / "db"; db.mkdir()
    index = _write(tmp_path / "index.tsv", "sample\tindex\nS1\tidx\n")
    fake = _write(
        tmp_path / "bin" / "k4neo-annotator",
        """#!/usr/bin/env bash
set -euo pipefail
if [[ ${1:-} == --version ]]; then echo 'k4neo 0.test'; exit 0; fi
prefix=''
while [[ $# -gt 0 ]]; do [[ $1 == --output ]] && { prefix=$2; shift 2; continue; }; shift; done
mkdir -p "$(dirname "$prefix")"
printf 'cts_id\\ttissue\\tsample_rate\\nSEQ1\\tblood\\t0\\n' > "${prefix}_healthy_sample_rate_raptor.tsv"
printf 'cts_id\\tannotation\\nSEQ1\\tnone\\n' > "${prefix}_annotated_raptor.tsv"
printf 'cts_id\\tuniqueness_rate\\nSEQ1\\t1\\n' > "${prefix}_uniqueness.tsv"
""",
        executable=True,
    )
    out = tmp_path / "k4out"
    rejected = subprocess.run(
        ["bash", "scripts/run_k4neo_sample.sh", "--query-table", str(query), "--database", str(db),
         "--index", str(index), "--outdir", str(out), "--k4neo-bin", str(fake)],
        cwd=ROOT, text=True, capture_output=True,
    )
    assert rejected.returncode != 0
    cp = _run(
        "bash", "scripts/run_k4neo_sample.sh", "--query-table", str(query), "--database", str(db),
        "--index", str(index), "--outdir", str(out), "--k4neo-bin", str(fake), "--license-accepted",
    )
    paths = [Path(x) for x in cp.stdout.splitlines() if x.strip()]
    assert len(paths) == 3
    assert all(path.is_file() for path in paths)
    manifest = json.loads((out / "k4neo_run_manifest.json").read_text())
    assert manifest["license_acknowledged"] is True
    assert manifest["healthy_sample_rate"]


def test_pvacsplice_wrapper_emits_all_epitopes_list(tmp_path: Path):
    junctions = _write(tmp_path / "cis.tsv", "chrom\tstart\tend\nchr1\t150\t200\n")
    vcf = _write(tmp_path / "vep.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    ref = _write(tmp_path / "ref.fa", ">chr1\nACGT\n")
    gtf = _write(tmp_path / "genes.gtf", "chr1\ttest\tgene\t1\t4\t.\t+\t.\tgene_id \"G1\";\n")
    fake = _write(
        tmp_path / "bin" / "pvacsplice",
        """#!/usr/bin/env bash
set -euo pipefail
if [[ ${1:-} == --version ]]; then echo 'pVACtools 7.test'; exit 0; fi
# positional run: junctions sample hla algorithms outdir vcf fasta gtf
out=${6}; sample=${3}
mkdir -p "$out/MHC_Class_I"
printf 'Index\\tHLA Allele\\tEpitope Seq\\n1\\tHLA-A*02:01\\tARNDCEQGH\\n' > "$out/MHC_Class_I/${sample}.all_epitopes.tsv"
""",
        executable=True,
    )
    out = tmp_path / "pvacsplice"
    cp = _run(
        "bash", "scripts/run_pvacsplice_sample.sh", "--junctions", str(junctions),
        "--annotated-vcf", str(vcf), "--sample-id", "S1", "--hla", "HLA-A*02:01",
        "--algorithms", "MHCflurry", "--outdir", str(out), "--ref-fasta", str(ref),
        "--gtf", str(gtf), "--pvacsplice-bin", str(fake),
    )
    list_path = Path(cp.stdout.strip().splitlines()[-1])
    assert list_path.is_file()
    assert "S1.all_epitopes.tsv" in list_path.read_text()
    manifest = json.loads((out / "pvacsplice_run_manifest.json").read_text())
    assert manifest["all_epitopes"]


def test_v051_driver_help_and_shell_syntax():
    _run("bash", "-n", "scripts/run_splice_provenance_v051.sh")
    cp = _run("bash", "scripts/run_splice_provenance_v051.sh", "--help")
    assert "Three-pass contract" in cp.stdout
