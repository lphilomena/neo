#!/usr/bin/env bash
# Install NetMHCstabpan (DTU tarball) or IEDB API shim (default when no tarball).
#
# Usage:
#   bash scripts/install_netmhcstabpan.sh [--iedb]
#   bash scripts/install_netmhcstabpan.sh [path/to/netMHCstabpan-1.0a.Linux.tar.gz]
#
# Default install dir: <project>/tools/netMHCstabpan
# Pipeline upstream uses IEDB API (see runner.py); local binary is optional for check-tools.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
NETMHCSTABPAN_HOME="${NETMHCSTABPAN_HOME:-${ROOT}/tools/netMHCstabpan}"
DATA_URL="https://services.healthtech.dtu.dk/services/NetMHCstabpan-1.0/data.tar.gz"
MODE="auto"
TARBALL=""

for arg in "$@"; do
  case "${arg}" in
    --iedb) MODE="iedb" ;;
    -h|--help)
      sed -n '2,10p' "$0"
      exit 0
      ;;
    *)
      if [[ -f "${arg}" ]]; then
        TARBALL="${arg}"
        MODE="dtu"
      fi
      ;;
  esac
done

if [[ "${MODE}" == "auto" && -z "${TARBALL}" ]]; then
  for candidate in \
    "${ROOT}/vendor/netMHCstabpan-1.0a.Linux.tar.gz" \
    "${ROOT}/vendor/netMHCstabpan-1.0cstatic.Linux.tar.gz" \
    "${ROOT}/vendor/netMHCstabpan-1.0.Linux.tar.gz"; do
    if [[ -f "${candidate}" ]]; then
      TARBALL="${candidate}"
      MODE="dtu"
      break
    fi
  done
fi

install_iedb_shim() {
  echo "==> Installing NetMHCstabpan IEDB shim to ${NETMHCSTABPAN_HOME} ..."
  mkdir -p "${NETMHCSTABPAN_HOME}"
  cat > "${NETMHCSTABPAN_HOME}/netMHCstabpan" <<'PY'
#!/usr/bin/env python3
"""NetMHCstabpan CLI shim — calls IEDB MHCI API (same backend as neoag upstream)."""
from __future__ import annotations

import argparse
import csv
import ssl
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

IEDB_URL = "https://tools-cluster-interface.iedb.org/tools_api/mhci/"
VERSION = "neoag-iedb-shim-1.0"


def iedb_row(peptide: str, allele: str) -> dict[str, str]:
    data = urllib.parse.urlencode({
        "method": "netmhcstabpan",
        "sequence_text": peptide,
        "allele": allele,
        "length": str(len(peptide)),
    }).encode()
    req = urllib.request.Request(IEDB_URL, data=data, method="POST")
    ctx = ssl.create_default_context()
    last_err = ""
    for attempt in range(6):
        try:
            with urllib.request.urlopen(req, context=ctx, timeout=300) as resp:
                text = resp.read().decode("utf-8", errors="ignore")
            if text.startswith("allele\t"):
                lines = [ln for ln in text.splitlines() if ln.strip()]
                if len(lines) >= 2:
                    header = lines[0].split("\t")
                    vals = lines[1].split("\t")
                    row = dict(zip(header, vals))
                    return {
                        "Peptide": peptide,
                        "HLA": allele,
                        "score": row.get("score", "0"),
                        "percentile_rank": row.get("percentile_rank", row.get("rank", "99")),
                    }
            last_err = text[:200]
        except urllib.error.HTTPError as exc:
            last_err = f"HTTP {exc.code}"
            if exc.code in {429, 500, 502, 503, 504} and attempt < 5:
                time.sleep(5 * (attempt + 1))
                continue
            raise
        except (urllib.error.URLError, TimeoutError) as exc:
            last_err = str(exc)
            if attempt < 5:
                time.sleep(5 * (attempt + 1))
                continue
            raise
        time.sleep(2)
    raise RuntimeError(f"IEDB netmhcstabpan failed for {peptide}/{allele}: {last_err}")


def load_pairs(path: str) -> list[tuple[str, str]]:
    pairs: list[tuple[str, str]] = []
    with open(path, newline="", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t") if "\t" in line else line.split()
            if len(parts) < 2:
                continue
            low = [p.lower() for p in parts]
            if "peptide" in low or "hla" in low or "allele" in low:
                continue
            pairs.append((parts[0], parts[1]))
    return pairs


def main() -> int:
    parser = argparse.ArgumentParser(
        description="NetMHCstabpan (IEDB API shim for neoag; DTU binary not required)."
    )
    parser.add_argument("--version", action="store_true", help="Print version")
    parser.add_argument("-f", "--file", help="TSV with peptide and HLA columns (no header or Peptide/HLA header)")
    parser.add_argument("-a", "--allele", action="append", default=[], help="HLA allele (repeatable)")
    parser.add_argument("-p", "--peptide", action="append", default=[], help="Peptide sequence (repeatable)")
    parser.add_argument("-o", "--out", default="-", help="Output TSV path (default stdout)")
    args, _rest = parser.parse_known_args()

    if args.version:
        print(VERSION)
        return 0

    pairs: list[tuple[str, str]] = []
    if args.file:
        pairs.extend(load_pairs(args.file))
    if args.peptide and args.allele:
        for pep in args.peptide:
            for hla in args.allele:
                pairs.append((pep, hla))
    if not pairs and not sys.stdin.isatty():
        pairs.extend(load_pairs("/dev/stdin"))

    if not pairs:
        print(
            "NetMHCstabpan IEDB shim — neoag upstream uses this API by default.\n"
            "Usage: netMHCstabpan -f pairs.tsv -o out.tsv\n"
            "       netMHCstabpan -p SIINFEKL -a HLA-A*02:01",
            file=sys.stderr,
        )
        return 0

    header = ["Peptide", "HLA", "score", "percentile_rank"]
    out_fh = sys.stdout if args.out == "-" else open(args.out, "w", encoding="utf-8", newline="")
    try:
        writer = csv.DictWriter(out_fh, fieldnames=header, delimiter="\t")
        writer.writeheader()
        for i, (peptide, allele) in enumerate(pairs, 1):
            writer.writerow(iedb_row(peptide, allele))
            if i % 25 == 0:
                time.sleep(2.0)
            time.sleep(0.3)
    finally:
        if out_fh is not sys.stdout:
            out_fh.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
PY
  chmod +x "${NETMHCSTABPAN_HOME}/netMHCstabpan"
  echo "==> IEDB shim installed. neoag upstream already calls IEDB in runner.py."
}

install_dtu() {
  WORK="$(mktemp -d)"
  trap 'rm -rf "${WORK}"' EXIT

  echo "==> Extracting ${TARBALL} ..."
  tar -xzf "${TARBALL}" -C "${WORK}"

  SRC=""
  if [[ -f "${WORK}/netMHCstabpan" ]]; then
    SRC="${WORK}"
  elif [[ -d "${WORK}/netMHCstabpan-1.0" ]]; then
    SRC="${WORK}/netMHCstabpan-1.0"
  else
    SRC="$(find "${WORK}" -maxdepth 3 -type f -name 'netMHCstabpan' | head -1 | xargs dirname 2>/dev/null || true)"
  fi

  if [[ -z "${SRC}" || ! -f "${SRC}/netMHCstabpan" ]]; then
    echo "ERROR: netMHCstabpan binary not found after extract. Contents:" >&2
    find "${WORK}" -maxdepth 3 | head -30 >&2
    exit 1
  fi

  echo "==> Installing DTU NetMHCstabpan to ${NETMHCSTABPAN_HOME} ..."
  mkdir -p "${NETMHCSTABPAN_HOME}"
  rsync -a "${SRC}/" "${NETMHCSTABPAN_HOME}/"
  chmod +x "${NETMHCSTABPAN_HOME}/netMHCstabpan"

  if [[ ! -d "${NETMHCSTABPAN_HOME}/data" ]]; then
    echo "==> Downloading data.tar.gz from DTU ..."
    curl -fsSL "${DATA_URL}" -o "${NETMHCSTABPAN_HOME}/data.tar.gz"
    tar -xzf "${NETMHCSTABPAN_HOME}/data.tar.gz" -C "${NETMHCSTABPAN_HOME}"
    rm -f "${NETMHCSTABPAN_HOME}/data.tar.gz"
  fi

  mkdir -p "${NETMHCSTABPAN_HOME}/tmp"

  if grep -q 'setenv[[:space:]]*NMHOME' "${NETMHCSTABPAN_HOME}/netMHCstabpan" 2>/dev/null; then
    sed -i "s|^setenv[[:space:]]*NMHOME.*|setenv\tNMHOME\t${NETMHCSTABPAN_HOME}|" \
      "${NETMHCSTABPAN_HOME}/netMHCstabpan"
    sed -i "s|^setenv[[:space:]]*TMPDIR[[:space:]]*/tmp|setenv  TMPDIR  ${NETMHCSTABPAN_HOME}/tmp|" \
      "${NETMHCSTABPAN_HOME}/netMHCstabpan"
  fi

  # Fix legacy DTU paths if present
  if grep -q '/net/sund-nas.win.dtu.dk' "${NETMHCSTABPAN_HOME}/netMHCstabpan" 2>/dev/null; then
    sed -i "s|/net/sund-nas.win.dtu.dk/storage/services/www/packages/netMHCstabpan/1.0/netMHCstabpan-1.0/data|${NETMHCSTABPAN_HOME}/data|g" \
      "${NETMHCSTABPAN_HOME}/netMHCstabpan"
  fi
}

if [[ "${MODE}" == "iedb" || ( "${MODE}" == "auto" && -z "${TARBALL}" ) ]]; then
  install_iedb_shim
elif [[ "${MODE}" == "dtu" ]]; then
  install_dtu
else
  cat <<EOF
NetMHCstabpan 安装包未找到。

选项 A — IEDB API（与 neoag 上游一致，无需 DTU 许可证）：
  bash scripts/install_netmhcstabpan.sh --iedb

选项 B — DTU 本地二进制（学术用户）：
1. 注册：https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=netMHCstabpan&version=1.0&packageversion=1.0cstatic&platform=Linux
2. 下载 Linux 包，例如 netMHCstabpan-1.0a.Linux.tar.gz
3. 安装：
   cp /path/to/netMHCstabpan-1.0a.Linux.tar.gz ${ROOT}/vendor/
   bash scripts/install_netmhcstabpan.sh ${ROOT}/vendor/netMHCstabpan-1.0a.Linux.tar.gz

安装目录（默认）: ${NETMHCSTABPAN_HOME}
EOF
  exit 1
fi

TOOLS_ENV="${ROOT}/conf/tools.env.sh"
if [[ -f "${TOOLS_ENV}" ]]; then
  if ! grep -q '^export NETMHCSTABPAN_HOME=' "${TOOLS_ENV}"; then
    cat >> "${TOOLS_ENV}" <<EOF

# NetMHCstabpan — install: bash scripts/install_netmhcstabpan.sh [--iedb]
export NETMHCSTABPAN_HOME="\${NEOAG_TOOLS_ROOT}/tools/netMHCstabpan"
if [[ -x "\${NETMHCSTABPAN_HOME}/netMHCstabpan" ]]; then
  export PATH="\${NETMHCSTABPAN_HOME}:\${PATH}"
fi
EOF
  fi
fi

echo "==> Testing ..."
export PATH="${NETMHCSTABPAN_HOME}:${PATH}"
if netMHCstabpan --version 2>&1 | head -1; then
  echo "==> NetMHCstabpan ready at ${NETMHCSTABPAN_HOME}"
  echo "    Add to shell: source ${TOOLS_ENV}"
else
  echo "WARN: netMHCstabpan --version failed" >&2
  exit 1
fi
