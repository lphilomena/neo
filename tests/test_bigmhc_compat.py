from __future__ import annotations

import sys
import csv
from types import SimpleNamespace

from neoag.tools import bigmhc_compat
from neoag.tools import runner
from neoag.utils import read_tsv


def test_bigmhc_compat_disables_mkldnn_and_forwards_arguments(tmp_path, monkeypatch):
    observed = tmp_path / "observed.txt"
    script = tmp_path / "predict.py"
    script.write_text(
        "import os, sys\n"
        f"open({str(observed)!r}, 'w').write('|'.join(sys.argv) + '\\n' + os.environ['OMP_NUM_THREADS'])\n",
        encoding="utf-8",
    )
    fake_torch = SimpleNamespace(
        backends=SimpleNamespace(mkldnn=SimpleNamespace(enabled=True)),
        set_num_threads=lambda value: setattr(fake_torch, "threads", value),
    )
    monkeypatch.setitem(sys.modules, "torch", fake_torch)
    monkeypatch.setenv("NEOAG_BIGMHC_CPU_THREADS", "2")
    monkeypatch.delenv("OMP_NUM_THREADS", raising=False)
    monkeypatch.setattr(sys, "argv", ["bigmhc_compat", str(script), "-m=im", "-d=cpu"])

    assert bigmhc_compat.main() == 0
    assert fake_torch.backends.mkldnn.enabled is False
    assert fake_torch.threads == 2
    assert observed.read_text(encoding="utf-8").splitlines() == [
        f"{script}|-m=im|-d=cpu",
        "2",
    ]


def test_bigmhc_runner_chunks_and_reuses_completed_chunks(tmp_path, monkeypatch):
    bigmhc_dir = tmp_path / "bigmhc"
    (bigmhc_dir / "src").mkdir(parents=True)
    (bigmhc_dir / "src/predict.py").write_text("# fixture\n", encoding="utf-8")
    ctx = SimpleNamespace(
        executables={"bigmhc_dir": str(bigmhc_dir), "bigmhc_python": sys.executable},
        sample_id="S1",
    )
    pairs = [(f"PEPTIDE{i}", "HLA-A*02:01") for i in range(5)]
    output = tmp_path / "bigmhc_im.tsv"
    calls = []

    def fake_run(cmd, **kwargs):
        calls.append(cmd)
        input_path = next(arg[3:] for arg in cmd if arg.startswith("-i="))
        output_path = next(arg[3:] for arg in cmd if arg.startswith("-o="))
        with open(input_path, encoding="utf-8", newline="") as source, open(
            output_path, "w", encoding="utf-8", newline=""
        ) as target:
            rows = list(csv.DictReader(source))
            writer = csv.DictWriter(target, fieldnames=["mhc", "pep", "BigMHC_IM"])
            writer.writeheader()
            for row in rows:
                writer.writerow({**row, "BigMHC_IM": "0.75"})
        return SimpleNamespace(returncode=0)

    monkeypatch.setenv("NEOAG_BIGMHC_CHUNK_SIZE", "2")
    monkeypatch.setattr(runner.subprocess if hasattr(runner, "subprocess") else __import__("subprocess"), "run", fake_run)
    runner._run_bigmhc_im_external(pairs, output, ctx)
    assert len(calls) == 3
    assert all(cmd[1:4] == ["-m", "neoag.tools.bigmhc_compat", str(bigmhc_dir / "src/predict.py")] for cmd in calls)
    assert len(read_tsv(output)) == 5

    runner._run_bigmhc_im_external(pairs, output, ctx)
    assert len(calls) == 3
