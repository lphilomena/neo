import sqlite3
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "build_normal_junction_index.py"


def test_normal_junction_index_preserves_frequency_metadata(tmp_path):
    source = tmp_path / "normal.tsv"
    source.write_text(
        "junction_id\tnormal_samples\tnormal_reads\tnormal_total_reads\tnormal_tissues\ttissue\tsource\tdataset\n"
        "chr1:101-200:+\t1\t2\t2\t1\tLIVER\tGTEx\tGTEx_v8\n"
        "chr1:101-200:+\t3\t8\t8\t2\tBLOOD\tGTEx\tGTEx_v8\n",
        encoding="utf-8",
    )
    output = tmp_path / "normal.sqlite"
    subprocess.run([sys.executable, str(SCRIPT), "--input", str(source), "--output", str(output)], check=True)
    subprocess.run([
        sys.executable, str(SCRIPT), "--input", str(source), "--output", str(output), "--check",
    ], check=True)
    with sqlite3.connect(output) as connection:
        row = connection.execute(
            "select normal_samples,normal_reads,normal_tissues from junction_ids where junction_id=?",
            ("chr1:101-200:+",),
        ).fetchone()
        version = connection.execute("select value from meta where key='schema_version'").fetchone()[0]
    assert row == (3, 8, 2)
    assert version == "2"
