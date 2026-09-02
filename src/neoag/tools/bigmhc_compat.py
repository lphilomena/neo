"""Run the external BigMHC predictor with a conservative CPU backend."""

from __future__ import annotations

import os
import runpy
import sys
from pathlib import Path


def main() -> int:
    if len(sys.argv) < 2:
        raise SystemExit("usage: python -m neoag.tools.bigmhc_compat PREDICT_PY [ARGS...]")
    predict_py = Path(sys.argv[1]).resolve()
    args = sys.argv[2:]
    default_threads = str(min(8, os.cpu_count() or 1))
    os.environ.setdefault(
        "OMP_NUM_THREADS", os.environ.get("NEOAG_BIGMHC_CPU_THREADS", default_threads)
    )
    os.environ.setdefault(
        "MKL_NUM_THREADS", os.environ.get("NEOAG_BIGMHC_CPU_THREADS", default_threads)
    )
    import torch

    if os.environ.get("NEOAG_BIGMHC_ENABLE_MKLDNN", "").strip().lower() not in {
        "1",
        "true",
        "yes",
    }:
        torch.backends.mkldnn.enabled = False
    torch.set_num_threads(
        max(1, int(os.environ.get("NEOAG_BIGMHC_CPU_THREADS", default_threads)))
    )
    sys.argv = [str(predict_py), *args]
    runpy.run_path(str(predict_py), run_name="__main__")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
