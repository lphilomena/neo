#!/usr/bin/env python3
"""Resumable parallel HTTP range downloader for large public references."""

from __future__ import annotations

import argparse
import concurrent.futures
import os
import time
import urllib.request
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--url", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--size", type=int, required=True)
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--chunk-mib", type=int, default=128)
    parser.add_argument("--retries", type=int, default=30)
    return parser.parse_args()


def download_part(url: str, part: Path, start: int, end: int, retries: int) -> tuple[Path, int]:
    expected = end - start + 1
    part.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(retries):
        current = part.stat().st_size if part.exists() else 0
        if current == expected:
            return part, expected
        if current > expected:
            part.unlink()
            current = 0
        request = urllib.request.Request(url, headers={"Range": f"bytes={start + current}-{end}"})
        try:
            with urllib.request.urlopen(request, timeout=120) as response, part.open("ab") as handle:
                while True:
                    block = response.read(4 * 1024 * 1024)
                    if not block:
                        break
                    handle.write(block)
        except Exception as exc:
            print(f"retry part={part.name} attempt={attempt + 1}: {exc}", flush=True)
            time.sleep(min(2 ** min(attempt, 5), 30))
            continue
    current = part.stat().st_size if part.exists() else 0
    if current != expected:
        raise RuntimeError(f"incomplete {part}: {current}/{expected}")
    return part, expected


def main() -> int:
    args = parse_args()
    if args.out.is_file() and args.out.stat().st_size == args.size:
        print(f"already complete: {args.out}")
        return 0
    chunk = args.chunk_mib * 1024 * 1024
    parts_dir = args.out.with_name(args.out.name + ".parts")
    tasks = []
    for index, start in enumerate(range(0, args.size, chunk)):
        end = min(start + chunk - 1, args.size - 1)
        tasks.append((args.url, parts_dir / f"part.{index:05d}", start, end, args.retries))
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(download_part, *task) for task in tasks]
        for done, future in enumerate(concurrent.futures.as_completed(futures), 1):
            part, size = future.result()
            print(f"completed {done}/{len(tasks)} {part.name} {size}", flush=True)
    assembling = args.out.with_name(args.out.name + ".assembling")
    with assembling.open("wb") as target:
        for _, part, _, _, _ in tasks:
            with part.open("rb") as source:
                while True:
                    block = source.read(8 * 1024 * 1024)
                    if not block:
                        break
                    target.write(block)
    if assembling.stat().st_size != args.size:
        raise RuntimeError(f"assembled size mismatch: {assembling.stat().st_size}/{args.size}")
    os.replace(assembling, args.out)
    print(f"downloaded: {args.out} ({args.size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
