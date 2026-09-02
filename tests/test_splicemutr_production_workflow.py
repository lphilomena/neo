from __future__ import annotations

import csv
import gzip
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def write_star(path: Path, reads: int = 12) -> None:
    path.write_text(f"chr1\t101\t200\t1\t1\t0\t{reads}\t0\t20\n", encoding="utf-8")


def test_prepare_splicemutr_cohort_orders_normal_reference_first(tmp_path: Path) -> None:
    tumor, normal1, normal2 = (tmp_path / name for name in ("tumor.SJ", "n1.SJ", "n2.SJ"))
    for path in (tumor, normal1, normal2):
        write_star(path)
    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "sample_id\trole\tstar_sj\n"
        f"T\ttumor\t{tumor}\nN1\tnormal\t{normal1}\nN2\tnormal\t{normal2}\n",
        encoding="utf-8",
    )
    out = tmp_path / "out"
    subprocess.run([sys.executable, str(ROOT / "scripts/prepare_splicemutr_cohort.py"),
                    "--samples", str(samples), "--outdir", str(out)], check=True)
    groups = (out / "groups_file.txt").read_text(encoding="utf-8").splitlines()
    assert groups == ["N1\t1", "N2\t1", "T\t2"]
    summary = json.loads((out / "cohort_summary.json").read_text(encoding="utf-8"))
    assert summary["status"] == "READY"
    assert summary["normal_is_reference_group"] is True


def test_prepare_splicemutr_cohort_blocks_underpowered_normal_group(tmp_path: Path) -> None:
    tumor, normal = tmp_path / "tumor.SJ", tmp_path / "normal.SJ"
    write_star(tumor); write_star(normal)
    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "sample_id\trole\tstar_sj\n" f"T\ttumor\t{tumor}\nN\tnormal\t{normal}\n",
        encoding="utf-8",
    )
    result = subprocess.run([sys.executable, str(ROOT / "scripts/prepare_splicemutr_cohort.py"),
                             "--samples", str(samples), "--outdir", str(tmp_path / "out")],
                            text=True, capture_output=True)
    assert result.returncode != 0
    assert "requires >= 2 normal" in result.stderr


def test_prepare_splicemutr_cohort_accepts_public_proxy_junctions(tmp_path: Path) -> None:
    tumor = tmp_path / "tumor.SJ"
    write_star(tumor)
    normal1, normal2 = tmp_path / "n1.junc", tmp_path / "n2.junc"
    for path in (normal1, normal2):
        path.write_text("chr1\t100\t201\tJUNC1\t12\t+\n", encoding="utf-8")
    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "sample_id\trole\tstar_sj\tjunction_file\tnormal_source\tnormal_match\ttissue\n"
        f"N1\tnormal\t\t{normal1}\tGTEx_V11\tPUBLIC_PROXY\tAdipose_Visceral_Omentum\n"
        f"N2\tnormal\t\t{normal2}\tGTEx_V11\tPUBLIC_PROXY\tAdipose_Visceral_Omentum\n"
        f"T\ttumor\t{tumor}\t\t\t\t\n",
        encoding="utf-8",
    )
    out = tmp_path / "out"
    subprocess.run([sys.executable, str(ROOT / "scripts/prepare_splicemutr_cohort.py"),
                    "--samples", str(samples), "--outdir", str(out)], check=True)
    summary = json.loads((out / "cohort_summary.json").read_text(encoding="utf-8"))
    assert summary["public_proxy_normal_samples"] == 2
    normalized = list(csv.DictReader(
        (out / "cohort_samples.normalized.tsv").open(encoding="utf-8"), delimiter="\t"
    ))
    assert normalized[0]["source_type"] == "LEAFCUTTER_JUNC"
    assert normalized[-1]["source_type"] == "STAR_SJ"


def test_splicemutr_workflow_and_open_neo_contract() -> None:
    workflow = (ROOT / "workflows/splicemutr/SpliceMutr.smk").read_text(encoding="utf-8")
    assert "rule leafcutter_differential" in workflow
    assert 'ANALYSIS_MODE == "outlier"' in workflow
    assert "extract_leafcutter_outlier_introns.R" in workflow
    assert "bsgenome_package" in workflow
    assert "run_splicemutr_form_parallel.sh" in workflow
    assert "leafcutter_outlier_adjust_method" in workflow
    assert "rule form_transcripts" in workflow
    assert "rule export_candidates" in workflow
    cli = (ROOT / "src/neoag/open_neo/cli.py").read_text(encoding="utf-8")
    assert 'run.add_argument("--splicemutr-config")' in cli
    assert 'run.add_argument("--splicemutr-samples")' in cli


def test_splicemutr_export_streams_gzip_normal_catalog(tmp_path: Path) -> None:
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "chr\tstart\tend\tstrand\tjuncs\tpeptide\tpep_junc_loc\ttx_id\tgene\terror\n"
        "chr1\t100\t200\t+\tchr1:100-200:+\tACDEFGHIKLMN\t5\tTX1\tGENE1\ttx\n",
        encoding="utf-8",
    )
    normal = tmp_path / "normal.tsv.gz"
    with gzip.open(normal, "wt", encoding="utf-8") as handle:
        handle.write("junction_id\tchromosome\nchr1:100-200:+\tchr1\n")
    output = tmp_path / "candidates.tsv"
    subprocess.run([
        sys.executable, str(ROOT / "scripts/export_splicemutr_candidates.py"),
        "--metadata", str(metadata), "--sample-id", "T", "--out", str(output),
        "--normal-junctions", str(normal),
    ], check=True)
    rows = list(csv.DictReader(output.open(encoding="utf-8"), delimiter="\t"))
    assert rows
    assert {row["normal_junction_status"] for row in rows} == {"DETECTED_BROAD_NORMAL"}


def test_splicemutr_snakefile_dry_run(tmp_path: Path) -> None:
    snakemake = shutil.which("snakemake")
    if not snakemake:
        return
    sj_files = []
    for sample in ("N1", "N2", "T"):
        path = tmp_path / f"{sample}.SJ.out.tab"
        write_star(path)
        sj_files.append(path)
    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "sample_id\trole\tstar_sj\n"
        f"N1\tnormal\t{sj_files[0]}\nN2\tnormal\t{sj_files[1]}\nT\ttumor\t{sj_files[2]}\n",
        encoding="utf-8",
    )
    txdb = tmp_path / "txdb.sqlite"; txdb.touch()
    config = tmp_path / "config.yaml"
    config.write_text(
        f"output_dir: {tmp_path / 'out'}\nsample_sheet: {samples}\ntarget_sample_id: T\n"
        f"splicemutr_home: {tmp_path / 'SpliceMutr'}\npython: {sys.executable}\nrscript: Rscript\n"
        f"leafcutter_cluster_script: cluster.py\nleafcutter_ds_script: leafcutter_ds.R\n"
        f"leafviz_prepare_results_script: prepare_results.R\nleafcutter_exon_file: exons.txt\n"
        f"leafcutter_annotation_prefix: annotation\ntxdb: {txdb}\nbsgenome_name: TestGenome\n",
        encoding="utf-8",
    )
    result = subprocess.run(
        [snakemake, "-n", "-s", str(ROOT / "workflows/splicemutr/SpliceMutr.smk"),
         "--configfile", str(config), "--cores", "1"],
        cwd=tmp_path, text=True, capture_output=True,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
    )
    assert result.returncode == 0, result.stderr + result.stdout
