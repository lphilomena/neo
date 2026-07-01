#!/usr/bin/env bash
# Install or repair NetMHCpan 4.x from DTU tarball (academic license required).
#
# Usage:
#   bash scripts/install_netmhcpan.sh [path/to/netMHCpan-4.2c.Linux.tar.gz]
#   bash scripts/install_netmhcpan.sh --repair
#
# Default install dir: <project>/tools/netMHCpan

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
NETMHCPAN_HOME="${NETMHCPAN_HOME:-${ROOT}/tools/netMHCpan}"
DATA_URL_42="https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/data.tar.gz"
DATA_URL_41="https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/data.tar.gz"
CONDA_BASE="${CONDA_BASE:-${NEOAG_CONDA_BASE:-$(conda info --base 2>/dev/null || echo ${HOME}/miniconda3)}}"
SYSROOT="${CONDA_BASE}/envs/neoag-tools/x86_64-conda-linux-gnu/sysroot"
LD_LINUX="${SYSROOT}/lib/ld-linux-x86-64.so.2"
PATCHELF="${CONDA_BASE}/envs/neoag-tools/bin/patchelf"

repair_netmhcpan_layout() {
  local home="$1"
  local bin_dir="${home}/Linux_x86_64/bin"
  local wrap_dir="${home}/.wrapper-bin"
  local link_dir="${home}/bin"

  if [[ ! -d "${bin_dir}" ]]; then
    echo "ERROR: ${bin_dir} not found; install NetMHCpan tarball first." >&2
    return 1
  fi

  mkdir -p "${home}/tmp" "${wrap_dir}" "${link_dir}"

  if [[ ! -x "${LD_LINUX}" ]]; then
    echo "WARN: conda sysroot loader missing at ${LD_LINUX}" >&2
    echo "      Run: conda install -n neoag-tools -c conda-forge sysroot_linux-64=2.34" >&2
  fi

  if [[ -x "${PATCHELF}" && -f "${LD_LINUX}" ]]; then
    echo "==> Patching NetMHCpan ELF binaries (patchelf) ..."
    for bin in "${bin_dir}"/*; do
      file "${bin}" 2>/dev/null | grep -q ELF || continue
      "${PATCHELF}" --set-interpreter "${LD_LINUX}" "${bin}" 2>/dev/null || true
      "${PATCHELF}" --set-rpath "${SYSROOT}/lib64:${SYSROOT}/lib" "${bin}" 2>/dev/null || true
    done
  fi

  echo "==> Creating ${wrap_dir} launchers ..."
  for real_bin in "${bin_dir}"/*; do
    [[ -f "${real_bin}" && -x "${real_bin}" ]] || continue
    file "${real_bin}" 2>/dev/null | grep -q ELF || continue
    name="$(basename "${real_bin}")"
    wrap="${wrap_dir}/${name}"
    cat > "${wrap}" <<EOF
#!/usr/bin/env bash
set -euo pipefail
if [[ -x "${LD_LINUX}" ]]; then
  exec "${LD_LINUX}" --library-path "${SYSROOT}/lib64:${SYSROOT}/lib" "${real_bin}" "\$@"
else
  exec "${real_bin}" "\$@"
fi
EOF
    chmod +x "${wrap}"
  done

  echo "==> Linking ${link_dir}/ -> .wrapper-bin ..."
  for wrap in "${wrap_dir}"/*; do
    [[ -x "${wrap}" ]] || continue
    name="$(basename "${wrap}")"
    ln -sf "../.wrapper-bin/${name}" "${link_dir}/${name}"
  done

  # Bash frontend (replaces tcsh launcher when present).
  if [[ ! -f "${home}/netMHCpan" ]] || grep -q 'tcsh' "${home}/netMHCpan" 2>/dev/null; then
    cat > "${home}/netMHCpan" <<'LAUNCHER'
#!/usr/bin/env bash
# NetMHCpan 4.2 frontend (conda sysroot for bundled binaries).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
export NETMHCPAN_HOME="${NETMHCPAN_HOME:-${ROOT}}"
export NETMHCpan="${NETMHCpan:-${NETMHCPAN_HOME}}"
export TMPDIR="${NEOAG_NETMHCPAN_TMPDIR:-${NETMHCPAN_HOME}/tmp}"
mkdir -p "${TMPDIR}"

WRAPPER="${NETMHCPAN_HOME}/.wrapper-bin/netMHCpan-4.2"
if [[ -x "${WRAPPER}" ]]; then
  exec "${WRAPPER}" "$@"
fi

PLATFORM_DIR="${NETMHCPAN_HOME}/Linux_x86_64"
BIN="${PLATFORM_DIR}/bin/netMHCpan-4.2"
if [[ ! -x "${BIN}" ]]; then
  echo "netMHCpan binary not found under ${NETMHCPAN_HOME}" >&2
  exit 127
fi

CONDA_BASE="${NEOAG_CONDA_BASE:-$(conda info --base 2>/dev/null || echo ${HOME}/miniconda3)}"
SYSROOT="${CONDA_BASE}/envs/neoag-tools/x86_64-conda-linux-gnu/sysroot"
LD_LINUX="${SYSROOT}/lib/ld-linux-x86-64.so.2"
if [[ -x "${LD_LINUX}" ]]; then
  exec "${LD_LINUX}" --library-path "${SYSROOT}/lib64:${SYSROOT}/lib" "${BIN}" "$@"
else
  exec "${BIN}" "$@"
fi
LAUNCHER
    chmod +x "${home}/netMHCpan"
  fi
}

MODE="install"
TARBALL=""
for arg in "$@"; do
  case "${arg}" in
    --repair) MODE="repair" ;;
    -h|--help)
      sed -n '2,10p' "$0"
      exit 0
      ;;
    *)
      if [[ -f "${arg}" ]]; then
        TARBALL="${arg}"
      fi
      ;;
  esac
done

if [[ "${MODE}" == "repair" ]]; then
  echo "==> Repairing NetMHCpan layout at ${NETMHCPAN_HOME} ..."
  repair_netmhcpan_layout "${NETMHCPAN_HOME}"
  export NETMHCPAN_HOME NETMHCpan="${NETMHCPAN_HOME}"
  export PATH="${NETMHCPAN_HOME}:${PATH}"
  echo "==> Smoke test ..."
  tmp="$(mktemp -d)"
  printf 'SIINFEKL HLA-A02:01 0\n' > "${tmp}/pep.pmhc"
  if netMHCpan -pmhc -BA -f "${tmp}/pep.pmhc" -t -99.9 2>&1 | grep -q PEPLIST; then
    echo "==> NetMHCpan repair OK at ${NETMHCPAN_HOME}"
  else
    echo "WARN: netMHCpan repair smoke test failed" >&2
    exit 1
  fi
  rm -rf "${tmp}"
  exit 0
fi

if [[ -z "${TARBALL}" ]]; then
  for candidate in \
    "${ROOT}/vendor/netMHCpan-4.2c.Linux.tar.gz" \
    "${ROOT}/vendor/netMHCpan-4.2.Linux.tar.gz" \
    "${ROOT}/vendor/netMHCpan-4.1b.linux.tar.gz" \
    "${ROOT}/vendor/netMHCpan-4.1.Linux.tar.gz" \
    "${ROOT}/vendor/netMHCpan-4.1b.Linux.tar.gz"; do
    if [[ -f "${candidate}" ]]; then
      TARBALL="${candidate}"
      break
    fi
  done
fi

if [[ -z "${TARBALL}" || ! -f "${TARBALL}" ]]; then
  cat <<EOF
NetMHCpan 安装包未找到。

请完成以下步骤后重试：

1. 在 DTU 注册（学术用户）：
   https://services.healthtech.dtu.dk/cgi-bin/request.cgi?tool_id=NetMHCpan

2. 从邮件链接下载 Linux 压缩包，例如 netMHCpan-4.2c.Linux.tar.gz

3. 复制到项目并安装：
   cp /path/to/netMHCpan-4.2c.Linux.tar.gz ${ROOT}/vendor/
   bash scripts/install_netmhcpan.sh ${ROOT}/vendor/netMHCpan-4.2c.Linux.tar.gz

已安装仅修复布局（wrapper + bin/ 链接）：
   bash scripts/install_netmhcpan.sh --repair

安装目录（默认）: ${NETMHCPAN_HOME}
EOF
  exit 1
fi

WORK="$(mktemp -d)"
trap 'rm -rf "${WORK}"' EXIT

echo "==> Extracting ${TARBALL} ..."
tar -xzf "${TARBALL}" -C "${WORK}"

SRC=""
if [[ -f "${WORK}/netMHCpan" ]]; then
  SRC="${WORK}"
elif [[ -d "${WORK}/netMHCpan-4.2" ]]; then
  SRC="${WORK}/netMHCpan-4.2"
elif [[ -d "${WORK}/netMHCpan-4.1" ]]; then
  SRC="${WORK}/netMHCpan-4.1"
else
  SRC="$(find "${WORK}" -maxdepth 3 -type f -name 'netMHCpan' | head -1 | xargs dirname 2>/dev/null || true)"
fi

if [[ -z "${SRC}" || ! -f "${SRC}/netMHCpan" ]]; then
  echo "ERROR: netMHCpan binary not found after extract. Contents:" >&2
  find "${WORK}" -maxdepth 3 | head -30 >&2
  exit 1
fi

echo "==> Installing to ${NETMHCPAN_HOME} ..."
mkdir -p "${NETMHCPAN_HOME}"
rsync -a "${SRC}/" "${NETMHCPAN_HOME}/"
chmod +x "${NETMHCPAN_HOME}/netMHCpan" 2>/dev/null || true

if [[ ! -d "${NETMHCPAN_HOME}/data" ]]; then
  echo "==> Downloading data.tar.gz from DTU ..."
  if curl -fsSL "${DATA_URL_42}" -o "${NETMHCPAN_HOME}/data.tar.gz"; then
    :
  else
    curl -fsSL "${DATA_URL_41}" -o "${NETMHCPAN_HOME}/data.tar.gz"
  fi
  tar -xzf "${NETMHCPAN_HOME}/data.tar.gz" -C "${NETMHCPAN_HOME}"
  rm -f "${NETMHCPAN_HOME}/data.tar.gz"
fi

repair_netmhcpan_layout "${NETMHCPAN_HOME}"

echo "==> Testing ..."
export NETMHCPAN_HOME NETMHCpan="${NETMHCPAN_HOME}"
export PATH="${NETMHCPAN_HOME}:${PATH}"
if netMHCpan -h 2>&1 | head -3; then
  echo "==> NetMHCpan installed at ${NETMHCPAN_HOME}"
  TOOLS_ENV="${ROOT}/conf/tools.env.sh"
  if [[ -f "${TOOLS_ENV}" ]] && ! grep -q 'NetMHCpan — installed via scripts/install_netmhcpan.sh' "${TOOLS_ENV}"; then
    cat >> "${TOOLS_ENV}" <<EOF

# NetMHCpan — installed via scripts/install_netmhcpan.sh
export NETMHCPAN_HOME="${NETMHCPAN_HOME}"
export NETMHCpan="${NETMHCPAN_HOME}"
export NEOAG_NETMHCPAN_BIN="${NETMHCPAN_HOME}/netMHCpan"
export NEOAG_NETMHCPAN_TMPDIR="${NEOAG_NETMHCPAN_TMPDIR:-${NETMHCPAN_HOME}/tmp}"
export PATH="${NETMHCPAN_HOME}:${NETMHCPAN_HOME}/bin:\${PATH}"
EOF
  fi
  echo "    source ${ROOT}/conf/tools.env.sh"
else
  echo "WARN: netMHCpan -h failed; try: bash scripts/install_netmhcpan.sh --repair" >&2
  exit 1
fi
