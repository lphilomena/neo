import os
from pathlib import Path
from unittest.mock import patch

import pytest

from neoag.adapters.prime import prime_parallel_jobs, split_peptide_chunks
from neoag.tools.registry import RunContext
from neoag.tools.runner import _run_prime_batch, _run_prime_external


def test_prime_parallel_jobs_default(monkeypatch):
    monkeypatch.delenv("NEOAG_PRIME_JOBS", raising=False)
    assert prime_parallel_jobs() == 4


def test_prime_parallel_jobs_env(monkeypatch):
    monkeypatch.setenv("NEOAG_PRIME_JOBS", "8")
    assert prime_parallel_jobs() == 8


def test_prime_parallel_jobs_invalid(monkeypatch):
    monkeypatch.setenv("NEOAG_PRIME_JOBS", "bad")
    assert prime_parallel_jobs(default=3) == 3


def test_prime_parallel_jobs_minimum_one(monkeypatch):
    monkeypatch.setenv("NEOAG_PRIME_JOBS", "0")
    assert prime_parallel_jobs() == 1


def test_split_peptide_chunks_single():
    peptides = ["AAA", "BBB", "CCC"]
    assert split_peptide_chunks(peptides, 1) == [peptides]


def test_split_peptide_chunks_even():
    peptides = ["A", "B", "C", "D", "E", "F"]
    chunks = split_peptide_chunks(peptides, 3)
    assert chunks == [["A", "B"], ["C", "D"], ["E", "F"]]


def test_split_peptide_chunks_more_jobs_than_peptides():
    peptides = ["A", "B"]
    chunks = split_peptide_chunks(peptides, 8)
    assert chunks == [["A"], ["B"]]


def test_run_prime_external_parallel_batches(tmp_path, monkeypatch):
    monkeypatch.setenv("NEOAG_PRIME_JOBS", "2")
    calls: list[int] = []

    def fake_batch(batch_dir, peptides, prime_alleles, prime_exe, mixmhcpred):
        calls.append(len(peptides))
        batch_dir = Path(batch_dir)
        batch_dir.mkdir(parents=True, exist_ok=True)
        raw_out = batch_dir / "prime_out.tsv"
        raw_out.write_text(
            "Peptide\tScore_A0201\t%Rank_A0201\n"
            + "\n".join(f"{p}\t0.1\t1.0" for p in peptides)
            + "\n",
            encoding="utf-8",
        )
        return raw_out

    pairs = [("PEP1", "HLA-A*02:01"), ("PEP2", "HLA-A*02:01"), ("PEP3", "HLA-A*02:01")]
    ctx = RunContext(sample_id="S1", outdir=tmp_path, stub=False, raw_peptides=tmp_path / "raw.tsv")
    out_tsv = tmp_path / "prime_evidence.tsv"

    with patch("neoag.tools.runner._resolve_prime_exe", return_value=Path("/fake/PRIME")), patch(
        "neoag.tools.runner._run_prime_batch", side_effect=fake_batch
    ):
        _run_prime_external(pairs, out_tsv, ctx)

    assert len(calls) == 2
    assert sum(calls) == 3
    assert out_tsv.is_file()
    text = out_tsv.read_text(encoding="utf-8")
    assert "PEP1" in text
    assert "PEP3" in text


def test_run_prime_external_groups_requested_pairs_by_allele(tmp_path, monkeypatch):
    monkeypatch.setenv("NEOAG_PRIME_JOBS", "2")
    calls: list[tuple[tuple[str, ...], tuple[str, ...]]] = []

    def fake_batch(batch_dir, peptides, prime_alleles, prime_exe, mixmhcpred):
        calls.append((tuple(peptides), tuple(prime_alleles)))
        batch_dir = Path(batch_dir)
        batch_dir.mkdir(parents=True, exist_ok=True)
        allele = prime_alleles[0]
        raw_out = batch_dir / "prime_out.tsv"
        raw_out.write_text(
            f"Peptide\tScore_{allele}\t%Rank_{allele}\n"
            + "\n".join(f"{p}\t0.1\t1.0" for p in peptides)
            + "\n",
            encoding="utf-8",
        )
        return raw_out

    pairs = [
        ("SHAREDPEP", "HLA-A*02:01"),
        ("SHAREDPEP", "HLA-B*07:02"),
        ("AONLYPEP", "HLA-A*02:01"),
    ]
    ctx = RunContext(sample_id="S1", outdir=tmp_path, stub=False, raw_peptides=tmp_path / "raw.tsv")
    out_tsv = tmp_path / "prime_evidence.tsv"

    with patch("neoag.tools.runner._resolve_prime_exe", return_value=Path("/fake/PRIME")), patch(
        "neoag.tools.runner._run_prime_batch", side_effect=fake_batch
    ):
        _run_prime_external(pairs, out_tsv, ctx)

    assert len(calls) == 2
    assert all(len(alleles) == 1 for _, alleles in calls)
    assert sum(len(peptides) for peptides, _ in calls) == 3
    assert {alleles[0] for _, alleles in calls} == {"A0201", "B0702"}
    text = out_tsv.read_text(encoding="utf-8")
    assert "HLA-A*02:01" in text
    assert "HLA-B*07:02" in text


def test_run_prime_batch_rejects_empty_success_output(tmp_path, monkeypatch):
    python = tmp_path / "prime-env" / "bin" / "python"
    python.parent.mkdir(parents=True)
    python.write_text("", encoding="utf-8")
    monkeypatch.setenv("NEOAG_PRIME_PYTHON", str(python))

    def fake_run(cmd, check, cwd, env):
        (Path(cwd) / "prime_out.tsv").write_text("\n", encoding="utf-8")

    with patch("subprocess.run", side_effect=fake_run):
        with pytest.raises(RuntimeError, match="empty output"):
            _run_prime_batch(
                tmp_path / "batch",
                ["SLYNTVATL"],
                ["A0201"],
                Path("/tools/PRIME"),
                "/tools/MixMHCpred",
            )


def test_run_prime_batch_prepends_configured_python(tmp_path, monkeypatch):
    python = tmp_path / "prime-env" / "bin" / "python"
    python.parent.mkdir(parents=True)
    python.write_text("", encoding="utf-8")
    monkeypatch.setenv("NEOAG_PRIME_PYTHON", str(python))
    monkeypatch.setenv("PATH", "/unrelated/optitype/bin:/usr/bin")

    def fake_run(cmd, check, cwd, env):
        (Path(cwd) / "prime_out.tsv").write_text(
            "Peptide\tScore_A0201\t%Rank_A0201\nSLYNTVATL\t0.1\t1.0\n",
            encoding="utf-8",
        )

    with patch("subprocess.run", side_effect=fake_run) as run:
        _run_prime_batch(
            tmp_path / "batch",
            ["SLYNTVATL"],
            ["A0201"],
            Path("/tools/PRIME"),
            "/tools/MixMHCpred",
        )

    env = run.call_args.kwargs["env"]
    assert env["MIXMHCPRED_PYTHON"] == str(python.resolve())
    assert env["PATH"].split(os.pathsep)[0] == str(python.resolve().parent)
