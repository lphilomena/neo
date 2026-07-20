from __future__ import annotations
import os
import shutil
import subprocess
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

from .prep import (
    unique_peptide_hla_pairs,
    write_peptide_fasta,
    write_mhcflurry_peptides_csv,
    write_netmhcpan_pmhc_input,
    netmhcpan_allele_string,
    mhcflurry_allele_list,
)
from .postprocess import vep_to_appm_tsv, lohhla_to_hla_loh_tsv, facets_to_purity_tsv
from .registry import TOOL_REGISTRY, RunContext, ROOT, resolve_runner_mode, RunnerMode, tool_docker_image


def strict_production_enabled() -> bool:
    return os.environ.get("NEOAG_STRICT_MODE", "").strip().lower() in {"1", "true", "yes", "production", "strict"}


def require_non_strict(reason: str) -> None:
    if strict_production_enabled():
        raise RuntimeError(f"Strict production mode forbids {reason}")


@dataclass
class ToolStatus:
    name: str
    executable: str
    available: bool
    resolved_path: str
    message: str


def _resolve_executable(name: str, override: str | None = None) -> tuple[bool, str, str]:
    import os
    if name == "vep":
        override = os.environ.get("NEOAG_VEP_BIN") or override
    if name == "bigmhc_im":
        if override and Path(override).is_file() and Path(override).name == "predict.py":
            return True, str(Path(override).resolve()), "BigMHC predict.py found"
        for candidate in (
            os.environ.get("BIGMHC_DIR", ""),
            f"{os.environ.get('NEOAG_TOOLS_ROOT', '')}/tools/bigmhc",
        ):
            if not candidate:
                continue
            predict_py = Path(candidate) / "src" / "predict.py"
            if predict_py.is_file():
                return True, str(predict_py.resolve()), "BigMHC predict.py found"
    if name == "prime":
        if override and Path(override).is_file():
            return True, str(Path(override).resolve()), "PRIME found"
        prime_path = os.environ.get("NEOAG_PRIME_BIN") or ""
        if not prime_path:
            prime_home = os.environ.get("PRIME_HOME", "")
            if prime_home:
                prime_path = str(Path(prime_home) / "PRIME")
        if prime_path and Path(prime_path).is_file():
            return True, str(Path(prime_path).resolve()), "PRIME found"
    if name == "netmhcpan":
        tools_root = os.environ.get("NEOAG_TOOLS_ROOT", "")
        home = os.environ.get("NETMHCPAN_HOME") or (f"{tools_root}/tools/netMHCpan" if tools_root else "")
        for candidate in (
            override,
            f"{home}/netMHCpan" if home else None,
            f"{home}/.wrapper-bin/netMHCpan-4.2" if home else None,
        ):
            if candidate and Path(candidate).is_file():
                return True, str(Path(candidate).resolve()), "netMHCpan found"
    if name == "star_fusion":
        tools_root = os.environ.get("NEOAG_TOOLS_ROOT", "")
        home = os.environ.get("NEOAG_STAR_FUSION_HOME") or (f"{tools_root}/tools/STAR-Fusion" if tools_root else "")
        for candidate in (
            override,
            f"{tools_root}/bin/star-fusion-neoag" if tools_root else None,
            f"{home}/STAR-Fusion" if home else None,
        ):
            if candidate and Path(candidate).is_file():
                return True, str(Path(candidate).resolve()), "STAR-Fusion found"
    if name == "fusioncatcher":
        tools_root = os.environ.get("NEOAG_TOOLS_ROOT", "")
        home = os.environ.get("NEOAG_FUSIONCATCHER_HOME") or (f"{tools_root}/tools/fusioncatcher" if tools_root else "")
        for candidate in (
            override,
            f"{tools_root}/bin/fusioncatcher-neoag" if tools_root else None,
            f"{home}/bin/fusioncatcher" if home else None,
        ):
            if candidate and Path(candidate).is_file():
                return True, str(Path(candidate).resolve()), "FusionCatcher found"
    if name == "easyfuse":
        tools_root = os.environ.get("NEOAG_TOOLS_ROOT", "")
        for candidate in (
            override,
            f"{tools_root}/bin/easyfuse-neoag" if tools_root else None,
        ):
            if candidate and Path(candidate).is_file():
                return True, str(Path(candidate).resolve()), "EasyFuse found"
    if name == "pyclone":
        for candidate in (
            override,
            os.environ.get("NEOAG_PYCLONE_BIN", ""),
            f"{os.environ.get('NEOAG_TOOLS_ROOT', '')}/bin/pyclone" if os.environ.get("NEOAG_TOOLS_ROOT") else None,
        ):
            if candidate and Path(candidate).is_file():
                return True, str(Path(candidate).resolve()), "PyClone-VI found"
    if name == "ascat":
        tools_root = os.environ.get("NEOAG_TOOLS_ROOT", "")
        home = os.environ.get("ASCAT_HOME") or (f"{tools_root}/tools/ascat" if tools_root else "")
        for candidate in (
            override,
            f"{home}/ascat.R" if home else None,
            f"{tools_root}/bin/ascat.R" if tools_root else None,
        ):
            if candidate and Path(candidate).is_file():
                return True, str(Path(candidate).resolve()), "ASCAT found"
    if name == "deepimmuno":
        try:
            from ..adapters.deepimmuno import resolve_deepimmuno_dir

            custom = override or os.environ.get("DEEPIMMUNO_DIR", "")
            deep_dir = resolve_deepimmuno_dir(custom or None)
            script = deep_dir / "deepimmuno-cnn.py"
            if script.is_file():
                return True, str(script.resolve()), "DeepImmuno-CNN found"
        except FileNotFoundError as exc:
            return False, override or "deepimmuno-cnn.py", str(exc)
    exe = override or TOOL_REGISTRY[name].executable
    path = shutil.which(exe)
    if path:
        return True, path, "found on PATH"
    if Path(exe).is_file():
        return True, str(Path(exe).resolve()), "found at explicit path"
    return False, exe, "not found on PATH"


def check_tool(name: str, executable: str | None = None) -> ToolStatus:
    if name not in TOOL_REGISTRY:
        return ToolStatus(name, executable or "", False, "", f"unknown tool: {name}")
    spec = TOOL_REGISTRY[name]
    ok, path, msg = _resolve_executable(name, executable or spec.executable)
    if ok and spec.version_args:
        try:
            subprocess.run(
                [path, *spec.version_args],
                capture_output=True,
                text=True,
                timeout=30,
                check=False,
            )
        except (OSError, subprocess.TimeoutExpired) as exc:
            ok = False
            msg = f"failed version check: {exc}"
    return ToolStatus(name, executable or spec.executable, ok, path if ok else "", msg)


def check_all_tools(executables: dict[str, str] | None = None) -> list[ToolStatus]:
    ex = executables or {}
    return [check_tool(n, ex.get(n)) for n in TOOL_REGISTRY]


def _stub_copy(fixture: str, dest: Path) -> None:
    src = Path(fixture)
    if not src.is_file():
        raise FileNotFoundError(f"Stub fixture missing: {src}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dest)


# ---------------------------------------------------------------------------
# Runner mode dispatch: conda (subprocess) vs docker (container)
# ---------------------------------------------------------------------------

_current_tool: str | None = None


@contextmanager
def _tool_context(tool_name: str):
    """Set the current tool name for the duration of one runner invocation."""
    global _current_tool
    old = _current_tool
    _current_tool = tool_name
    try:
        yield
    finally:
        _current_tool = old


def _effective_mode(tool_name: str) -> RunnerMode:
    """Return the effective runner mode for *tool_name*.

    Docker mode is only used when:
    1. NEOAG_RUNNER_MODE is ``docker``, AND
    2. The tool has a known Docker image in DOCKER_IMAGES.
    Otherwise conda mode is used (silent fallback).
    """
    if resolve_runner_mode() != RunnerMode.DOCKER:
        return RunnerMode.CONDA
    if tool_docker_image(tool_name) is None:
        return RunnerMode.CONDA
    return RunnerMode.DOCKER


def _collect_docker_mounts(cmd: list[str], workdir: Path) -> list[str]:
    """Collect ``-v host:container`` mount arguments for a Docker run.

    Strategy: mount every filesystem root that the command touches
    (project root, tools root, workdir tree, and any external data mounts)
    at the identical path inside the container so all host paths resolve.

    Additionally, for well-known tools, remap host paths to the container's
    expected default locations (e.g. VEP cache → /opt/vep/.vep).
    """
    mounts: dict[str, str] = {}  # host_path → container_path (not always identical)

    def _add(host: Path, container: Path | None = None) -> None:
        try:
            resolved = host.resolve()
        except (OSError, RuntimeError):
            return
        for p in [resolved] + list(resolved.parents):
            if p.exists():
                target = str(container) if container else str(p)
                mounts[str(p)] = target
                return

    # Always mount the project root and tools root
    project_root = os.environ.get("NEOAG_PROJECT_ROOT", "")
    if project_root and Path(project_root).exists():
        mounts[project_root] = project_root
    tools_root = os.environ.get("NEOAG_TOOLS_ROOT", "")
    if tools_root and Path(tools_root).exists() and tools_root != project_root:
        mounts[tools_root] = tools_root

    # Mount workdir tree
    _add(workdir)

    # Mount directories referenced in command arguments
    for arg in cmd:
        if arg.startswith("-"):
            continue
        p = Path(arg)
        if p.exists():
            mounts[str(p.resolve())] = str(p.resolve())
        elif p.parent.exists():
            mounts[str(p.parent.resolve())] = str(p.parent.resolve())

    # ---- Tool-specific remaps ----
    # VEP: the official container expects cache at /opt/vep/.vep
    vep_cache = os.environ.get("NEOAG_VEP_CACHE", "").strip()
    if vep_cache and Path(vep_cache).exists():
        mounts[vep_cache] = "/opt/vep/.vep"

    # pVACtools: the official container has IEDB at /opt/iedb
    # No action needed — pvacseq auto-detects /opt/iedb

    # Mount external data directories from environment (host=container)
    for env_var in ("NEOAG_REFERENCE_FASTA",
                    "NETMHCPAN_HOME", "MHCFLURRY_DOWNLOADS_DIR",
                    "BIGMHC_DIR", "PRIME_HOME", "MIXMHCPRED_HOME",
                    "CTAT_GENOME_LIB", "NEOAG_EASYFUSE_REF",
                    "NEOAG_CTAT_LIB_DIR", "NEOAG_SHARED_REF_DIR"):
        val = os.environ.get(env_var, "").strip()
        if val:
            p = Path(val)
            if p.exists():
                mounts[str(p.resolve())] = str(p.resolve())
            else:
                # Mount the nearest existing parent (file may not exist yet on dev)
                _add(p)

    # Mount conda base so tools that still use conda envs (NetMHCpan etc.) can find shared libs
    conda_base = os.environ.get("NEOAG_CONDA_BASE", "")
    if conda_base and Path(conda_base).exists():
        mounts[conda_base] = conda_base

    return [arg for pair in (("-v", f"{h}:{c}") for h, c in mounts.items()) for arg in pair]


def _run_docker(image: str, cmd: list[str], workdir: Path, *, env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    """Execute *cmd* inside a Docker container from *image*.

    Returns the completed process so callers can inspect stdout/stderr.
    """
    workdir.mkdir(parents=True, exist_ok=True)

    mounts = _collect_docker_mounts(cmd, workdir)

    # Build env-var passthrough: always pass core vars, then add caller's env dict.
    #
    # PATH handling:
    #   - neoag-* custom images ship their own tool binaries.  Don't pass PATH
    #     at all — let the container's built-in ENV PATH take effect.  Passing
    #     the host PATH would shadow container binaries with conda Docker wrapper
    #     scripts, causing recursive "docker: not found" failures.
    #   - Other images (ensembl-vep, pvacseq, etc.) may need host conda env
    #     paths for shebangs like python3.11 to resolve correctly.
    if image.startswith("neoag-"):
        system_path = None  # let container use its own PATH
    else:
        system_path = (
            f"{os.environ.get('NEOAG_CONDA_BASE', '')}/envs/neoag-tools/bin:"
            f"{os.environ.get('NEOAG_CONDA_BASE', '')}/envs/neoag-vep/bin:"
            f"{os.environ.get('PATH', '/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin')}"
        )

    passthrough = {
        "NEOAG_PROJECT_ROOT": os.environ.get("NEOAG_PROJECT_ROOT", ""),
        "NEOAG_TOOLS_ROOT": os.environ.get("NEOAG_TOOLS_ROOT", ""),
        "TMPDIR": os.environ.get("TMPDIR", "/tmp"),
        # Conda env libs (TensorFlow, CUDA, etc.)
        "LD_LIBRARY_PATH": os.environ.get("LD_LIBRARY_PATH", ""),
        # MHCflurry needs these
        "TF_USE_LEGACY_KERAS": os.environ.get("TF_USE_LEGACY_KERAS", "1"),
        "MHCFLURRY_DOWNLOADS_DIR": os.environ.get("MHCFLURRY_DOWNLOADS_DIR", ""),
    }
    if system_path is not None:
        passthrough["PATH"] = system_path
    if env:
        passthrough.update(env)
    # For neoag-* images, ensure the container uses its own environment.
    # _netmhcpan_subprocess_env() returns a copy of the full host environ
    # which includes conda PATH and host tool paths — these would shadow
    # the container's own tool binaries with host Docker wrapper scripts.
    if system_path is None:
        passthrough.pop("PATH", None)
        passthrough.pop("NETMHCPAN_HOME", None)
        passthrough.pop("NETMHCpan", None)

    # Custom (neoag-*) images use ENTRYPOINT ["/bin/bash", "-lc"] which
    # wraps the command; "exec" inside wrapper scripts fails under that.
    # Bypass the entrypoint so the command runs directly.
    entrypoint_args = ["--entrypoint", ""] if image.startswith("neoag-") else []

    # Resolve bare executable names to full host paths so they are
    # accessible inside the container via volume mounts.
    # Skip resolution for neoag-* images: entrypoint bypass means the
    # command runs directly inside the container.  Resolving to a host
    # conda wrapper would cause recursive "docker: not found" errors.
    if cmd and not os.path.isabs(cmd[0]) and "/" not in cmd[0] and not entrypoint_args:
        resolved = shutil.which(cmd[0])
        if resolved:
            cmd = [resolved, *cmd[1:]]

    docker_cmd = [
        "/usr/bin/docker", "run", "--rm",
        *entrypoint_args,
        "--workdir", str(workdir),
        *[arg for pair in (("-e", f"{k}={v}") for k, v in passthrough.items()) for arg in pair],
        *mounts,
        image,
        *cmd,
    ]

    proc = subprocess.run(
        docker_cmd,
        cwd=workdir,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"Docker command failed ({proc.returncode}): {' '.join(cmd)}\n"
            f"Image: {image}\n"
            f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
        )
    return proc


def _run_cmd(cmd: list[str], workdir: Path, *, env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    """Execute *cmd* using the currently active runner mode.

    In Docker mode the command runs inside a container with volume mounts;
    in conda mode it runs as a local subprocess.  *env* supplies extra
    environment variables (merged on top of ``os.environ`` for conda,
    passed as ``-e`` flags for Docker).

    Returns the completed process so callers can inspect stdout/stderr.
    Raises RuntimeError if the command exits non-zero.
    """
    workdir.mkdir(parents=True, exist_ok=True)

    # Determine mode based on current tool context
    tool = _current_tool
    if tool and _effective_mode(tool) == RunnerMode.DOCKER:
        image = tool_docker_image(tool)
        if image:
            return _run_docker(image, cmd, workdir, env=env)

    # Default: conda mode (subprocess)
    subprocess_env = os.environ.copy()
    if env:
        subprocess_env.update(env)
    proc = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True, env=subprocess_env)
    if proc.returncode != 0:
        raise RuntimeError(
            f"Command failed ({proc.returncode}): {' '.join(cmd)}\n"
            f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
        )
    return proc


def _iedb_mhci_row(method: str, peptide: str, allele: str) -> dict[str, str]:
    import ssl
    import time
    import urllib.error
    import urllib.parse
    import urllib.request

    data = urllib.parse.urlencode({
        "method": method,
        "sequence_text": peptide,
        "allele": allele,
        "length": str(len(peptide)),
    }).encode()
    req = urllib.request.Request(
        "https://tools-cluster-interface.iedb.org/tools_api/mhci/",
        data=data,
        method="POST",
    )
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
                    return dict(zip(header, vals))
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
    raise RuntimeError(f"IEDB {method} failed for {peptide}/{allele}: {last_err}")


def _run_netmhcpan_iedb(pairs: list[tuple[str, str]], out_xls: Path) -> None:
    import time

    work = out_xls.parent / "netmhcpan"
    work.mkdir(parents=True, exist_ok=True)
    progress = work / "iedb_progress.tsv"
    done: set[tuple[str, str]] = set()
    rows: list[dict[str, str]] = []
    if progress.is_file():
        import csv
        with progress.open(newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                key = (row.get("Peptide", ""), row.get("HLA", ""))
                if key[0] and key[1]:
                    done.add(key)
                    rows.append(row)

    total = len(pairs)
    header = ["Pos", "Peptide", "HLA", "Score_EL", "%Rank_EL", "Score_BA", "%Rank_BA", "BindLevel"]
    with progress.open("a" if done else "w", encoding="utf-8", newline="") as fh:
        if not done:
            fh.write("\t".join(header) + "\n")
        writer = None
        if done:
            import csv
            writer = csv.DictWriter(fh, fieldnames=header, delimiter="\t", extrasaction="ignore")

        for i, (peptide, allele) in enumerate(pairs):
            if (peptide, allele) in done:
                continue
            ba = _iedb_mhci_row("netmhcpan_ba", peptide, allele)
            time.sleep(0.3)
            el = _iedb_mhci_row("netmhcpan_el", peptide, allele)
            row = {
                "Pos": "0",
                "Peptide": peptide,
                "HLA": allele,
                "Score_EL": el.get("score", el.get("Score_EL", "0")),
                "%Rank_EL": el.get("percentile_rank", el.get("%Rank_EL", "99")),
                "Score_BA": ba.get("ic50", ba.get("Score_BA", "0")),
                "%Rank_BA": ba.get("percentile_rank", ba.get("%Rank_BA", "99")),
                "BindLevel": "",
            }
            rows.append(row)
            done.add((peptide, allele))
            if writer is None:
                fh.write("\t".join(row.get(col, "") for col in header) + "\n")
            else:
                writer.writerow(row)
            fh.flush()
            completed = len(done)
            if completed % 10 == 0:
                print(f"netmhcpan IEDB: {completed}/{total}", flush=True)
            if completed % 25 == 0:
                time.sleep(2.0)

    out_xls.parent.mkdir(parents=True, exist_ok=True)
    with out_xls.open("w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(row.get(col, "") for col in header) + "\n")


def _write_stabpan_stub(pairs: list[tuple[str, str]], out_tsv: Path) -> None:
    header = ["Peptide", "HLA", "score", "percentile_rank"]
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(header) + "\n")
        for peptide, allele in pairs:
            h = abs(hash(f"{peptide}_{allele}"))
            score = f"{5.0 + (h % 500) / 100.0:.2f}"
            rank = f"{0.1 + (h % 990) / 10.0:.2f}"
            fh.write("\t".join([peptide, allele, score, rank]) + "\n")


def _run_netmhcstabpan_iedb(pairs: list[tuple[str, str]], out_tsv: Path) -> None:
    import time

    work = out_tsv.parent / "netmhcstabpan"
    work.mkdir(parents=True, exist_ok=True)
    progress = work / "iedb_progress.tsv"
    done: set[tuple[str, str]] = set()
    rows: list[dict[str, str]] = []
    if progress.is_file():
        import csv
        with progress.open(newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                key = (row.get("Peptide", ""), row.get("HLA", ""))
                if key[0] and key[1]:
                    done.add(key)
                    rows.append(row)

    total = len(pairs)
    header = ["Peptide", "HLA", "score", "percentile_rank"]
    with progress.open("a" if done else "w", encoding="utf-8", newline="") as fh:
        if not done:
            fh.write("\t".join(header) + "\n")
        writer = None
        if done:
            import csv
            writer = csv.DictWriter(fh, fieldnames=header, delimiter="\t", extrasaction="ignore")

        for peptide, allele in pairs:
            if (peptide, allele) in done:
                continue
            resp = _iedb_mhci_row("netmhcstabpan", peptide, allele)
            row = {
                "Peptide": peptide,
                "HLA": allele,
                "score": resp.get("score", "0"),
                "percentile_rank": resp.get("percentile_rank", resp.get("rank", "99")),
            }
            rows.append(row)
            done.add((peptide, allele))
            if writer is None:
                fh.write("\t".join(row.get(col, "") for col in header) + "\n")
            else:
                writer.writerow(row)
            fh.flush()
            completed = len(done)
            if completed % 10 == 0:
                print(f"netmhcstabpan IEDB: {completed}/{total}", flush=True)
            if completed % 25 == 0:
                time.sleep(2.0)
            time.sleep(0.3)

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(row.get(col, "") for col in header) + "\n")


def run_netmhcstabpan(ctx: RunContext, out_tsv: Path) -> Path:
    spec = TOOL_REGISTRY["netmhcstabpan"]
    if not ctx.raw_peptides:
        raise ValueError("netmhcstabpan requires raw_peptides.tsv")
    pairs = unique_peptide_hla_pairs(ctx.raw_peptides)
    if not pairs:
        raise ValueError(f"netmhcstabpan: no peptide/HLA pairs in {ctx.raw_peptides}")
    out_tsv = out_tsv.resolve()
    if ctx.stub:
        _write_stabpan_stub(pairs, out_tsv)
        return out_tsv
    _run_netmhcstabpan_iedb(pairs, out_tsv)
    return out_tsv


def _netmhcpan_subprocess_env() -> dict[str, str]:
    env = os.environ.copy()
    home = env.get("NETMHCPAN_HOME") or env.get("NETMHCpan") or ""
    if not home:
        tools_root = env.get("NEOAG_TOOLS_ROOT", "")
        if tools_root:
            home = str(Path(tools_root) / "tools" / "netMHCpan")
    if home:
        env["NETMHCPAN_HOME"] = home
        env["NETMHCpan"] = home
        tmpdir = env.get("NEOAG_NETMHCPAN_TMPDIR") or str(Path(home) / "tmp")
        Path(tmpdir).mkdir(parents=True, exist_ok=True)
        env["TMPDIR"] = tmpdir
    return env


def resolve_netmhcpan_backend() -> str:
    """local = project NetMHCpan 4.2; iedb = IEDB API only."""
    mode = os.environ.get("NEOAG_NETMHCPAN_BACKEND", "local").strip().lower()
    return "iedb" if mode in {"iedb", "api", "remote"} else "local"


def netmhcpan_iedb_fallback_enabled() -> bool:
    return os.environ.get("NEOAG_NETMHCPAN_ALLOW_IEDB_FALLBACK", "").strip().lower() in {"1", "true", "yes"}


def _run_netmhcpan_local(
    pairs: list[tuple[str, str]],
    out_xls: Path,
    *,
    netmhcpan_exe: str,
) -> None:
    from ..adapters.netmhcpan import parse_netmhcpan_local_stdout, write_netmhcpan_standard_xls

    work = out_xls.parent / "netmhcpan"
    work.mkdir(parents=True, exist_ok=True)
    chunk_size = int(os.environ.get("NEOAG_NETMHCPAN_LOCAL_CHUNK_SIZE", "500") or "500")
    chunk_size = max(chunk_size, 1)
    all_rows: list[dict[str, str]] = []
    for chunk_idx, start_idx in enumerate(range(0, len(pairs), chunk_size), start=1):
        chunk = pairs[start_idx:start_idx + chunk_size]
        pep_mhc = work / f"peptides.chunk{chunk_idx:04d}.pmhc"
        write_netmhcpan_pmhc_input(chunk, pep_mhc)
        cmd = [
            netmhcpan_exe,
            "-pmhc",
            "-BA",
            "-f",
            str(pep_mhc),
            "-t",
            "-99.9",
        ]
        proc = _run_cmd(cmd, work, env=_netmhcpan_subprocess_env())
        if proc.returncode != 0:
            raise RuntimeError(
                f"NetMHCpan 4.2 failed ({proc.returncode}) on chunk {chunk_idx} "
                f"({len(chunk)} pairs, {start_idx + len(chunk)}/{len(pairs)} total): {' '.join(cmd)}\n"
                f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
            )
        rows = parse_netmhcpan_local_stdout(proc.stdout, source=str(out_xls))
        if not rows:
            raise RuntimeError(
                f"NetMHCpan 4.2 produced no binding rows for chunk {chunk_idx} "
                f"({len(chunk)} peptide/HLA pairs).\nstdout tail:\n{proc.stdout[-2000:]}"
            )
        all_rows.extend(rows)
        print(f"netmhcpan local: {min(start_idx + chunk_size, len(pairs))}/{len(pairs)}", flush=True)
    write_netmhcpan_standard_xls(out_xls, all_rows)


def _run_netmhcpan_local_by_allele(
    pairs: list[tuple[str, str]],
    out_xls: Path,
    *,
    netmhcpan_exe: str,
) -> None:
    """Run NetMHCpan locally without -pmhc.

    Some NetMHCpan 4.2 builds crash in peptide-MHC input mode on larger chunks.
    The standard peptide input mode is slower but more robust: group by allele,
    pass unique peptides with -p, and parse stdout into the same standard table.
    """
    from ..adapters.netmhcpan import parse_netmhcpan_local_stdout, write_netmhcpan_standard_xls

    work = out_xls.parent / "netmhcpan_by_allele"
    work.mkdir(parents=True, exist_ok=True)
    by_allele: dict[str, list[str]] = {}
    seen: dict[str, set[str]] = {}
    for peptide, hla in pairs:
        if not peptide or not hla:
            continue
        if hla not in seen:
            seen[hla] = set()
            by_allele[hla] = []
        if peptide in seen[hla]:
            continue
        seen[hla].add(peptide)
        by_allele[hla].append(peptide)

    all_rows: list[dict[str, str]] = []
    completed = 0
    total = sum(len(v) for v in by_allele.values())
    for hla, peptides in sorted(by_allele.items()):
        allele = netmhcpan_allele_string([hla])
        pep_file = work / f"{allele.replace('*', '').replace(':', '').replace(',', '_')}.pep"
        pep_file.write_text("\n".join(peptides) + "\n", encoding="utf-8")
        cmd = [
            netmhcpan_exe,
            "-p",
            "-BA",
            "-a",
            allele,
            "-f",
            str(pep_file),
            "-t",
            "-99.9",
        ]
        proc = _run_cmd(cmd, work, env=_netmhcpan_subprocess_env())
        if proc.returncode != 0:
            raise RuntimeError(
                f"NetMHCpan allele-mode failed ({proc.returncode}) for {hla}: {' '.join(cmd)}\n"
                f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
            )
        rows = parse_netmhcpan_local_stdout(proc.stdout, source=str(out_xls))
        if not rows:
            raise RuntimeError(
                f"NetMHCpan allele-mode produced no rows for {hla}; stdout tail:\n{proc.stdout[-2000:]}"
            )
        all_rows.extend(rows)
        completed += len(peptides)
        print(f"netmhcpan local allele-mode: {completed}/{total}", flush=True)
    write_netmhcpan_standard_xls(out_xls, all_rows)

def _set_netmhcpan_provenance(ctx: RunContext, out_xls: Path, backend: str, *, stub: bool = False) -> None:
    from ..evidence_provenance import provenance_for_netmhcpan_run

    ctx.tool_provenance["netmhcpan"] = provenance_for_netmhcpan_run(out_xls, backend, stub=stub)


def run_netmhcpan(ctx: RunContext, out_xls: Path) -> Path:
    spec = TOOL_REGISTRY["netmhcpan"]
    if ctx.stub:
        require_non_strict("stub NetMHCpan output")
        _stub_copy(spec.fixture_outputs["xls"], out_xls)
        _set_netmhcpan_provenance(ctx, out_xls, "local", stub=True)
        return out_xls
    if not ctx.raw_peptides:
        raise ValueError("netmhcpan requires raw_peptides.tsv")
    pairs = unique_peptide_hla_pairs(ctx.raw_peptides)
    if not pairs:
        raise ValueError(f"netmhcpan: no peptide/HLA pairs in {ctx.raw_peptides}")
    out_xls = out_xls.resolve()
    if resolve_netmhcpan_backend() != "iedb":
        try:
            _run_netmhcpan_local(pairs, out_xls, netmhcpan_exe=ctx.exe("netmhcpan"))
            _set_netmhcpan_provenance(ctx, out_xls, "local")
            return out_xls
        except Exception as exc:
            try:
                print(f"netmhcpan: -pmhc local failed ({exc}); retrying local allele mode", flush=True)
                _run_netmhcpan_local_by_allele(pairs, out_xls, netmhcpan_exe=ctx.exe("netmhcpan"))
                _set_netmhcpan_provenance(ctx, out_xls, "local")
                return out_xls
            except Exception as allele_exc:
                exc = allele_exc
            if not netmhcpan_iedb_fallback_enabled():
                raise RuntimeError(
                    "NetMHCpan local mode failed and IEDB fallback is disabled by default. "
                    "Set NEOAG_NETMHCPAN_ALLOW_IEDB_FALLBACK=1 to allow remote fallback."
                ) from exc
            if strict_production_enabled():
                raise RuntimeError(f"Strict production mode forbids NetMHCpan silent fallback after local failure: {exc}") from exc
            print(f"netmhcpan: local failed ({exc}); falling back to IEDB API", flush=True)
    _run_netmhcpan_iedb(pairs, out_xls)
    _set_netmhcpan_provenance(ctx, out_xls, "iedb")
    return out_xls


def run_mhcflurry(ctx: RunContext, out_csv: Path) -> Path:
    spec = TOOL_REGISTRY["mhcflurry"]
    if ctx.stub:
        require_non_strict("stub MHCflurry output")
        _stub_copy(spec.fixture_outputs["csv"], out_csv)
        return out_csv
    if not ctx.raw_peptides:
        raise ValueError("mhcflurry requires raw_peptides.tsv")
    pairs = unique_peptide_hla_pairs(ctx.raw_peptides)
    if not pairs:
        raise ValueError(f"mhcflurry: no peptide/HLA pairs in {ctx.raw_peptides}")
    work = (ctx.outdir / "tools" / "mhcflurry").resolve()
    pep_csv = work / "peptides.csv"
    out_csv = out_csv.resolve()
    write_mhcflurry_peptides_csv(pairs, pep_csv)
    cmd = [
        ctx.exe("mhcflurry"),
        str(pep_csv),
        "--out", str(out_csv),
    ]
    _run_cmd(cmd, work)
    return out_csv


def run_pvacseq(ctx: RunContext, out_tsv: Path) -> Path:
    spec = TOOL_REGISTRY["pvacseq"]
    if ctx.stub:
        _stub_copy(spec.fixture_outputs["aggregated"], out_tsv)
        return out_tsv
    if not ctx.tumor_vcf:
        raise ValueError("pvacseq requires tumor_vcf (or tools.stub=true)")
    if not ctx.hla_alleles:
        raise ValueError("pvacseq requires hla_alleles")
    if ctx.normal_vcf and ctx.normal_vcf != ctx.tumor_vcf:
        raise ValueError(
            "pvacseq expects a single combined VCF in tumor_vcf; "
            "set normal_sample_name for the matched-normal column"
        )
    work = (ctx.outdir / "tools" / "pvacseq").resolve()
    hla = ",".join(ctx.hla_alleles)
    tumor_sample = ctx.tumor_sample_name or ctx.sample_id
    algos = ctx.prediction_algorithms.replace(",", " ").split()
    cmd = [
        ctx.exe("pvacseq"), "run",
        str(ctx.tumor_vcf.resolve()),
        tumor_sample,
        hla,
        *algos,
        str(work),
        "-e1", "8,9,10,11",
        "-e2", "12,13,14,15,16,17,18,19,20,21",
    ]
    if ctx.pass_only:
        cmd.append("--pass-only")
    if ctx.normal_sample_name:
        cmd.extend(["--normal-sample-name", ctx.normal_sample_name])
    if ctx.phased_vcf:
        cmd.extend(["-p", str(ctx.phased_vcf.resolve())])
    _run_cmd(cmd, work)
    candidates = sorted(work.rglob("*aggregated*.tsv"))
    if not candidates:
        candidates = sorted(work.rglob("*.tsv"))
    if not candidates:
        raise FileNotFoundError(f"pvacseq finished but no aggregated TSV in {work}")
    shutil.copy2(candidates[0], out_tsv)
    return out_tsv


def run_pvacfuse(ctx: RunContext, out_tsv: Path) -> Path:
    spec = TOOL_REGISTRY["pvacfuse"]
    if ctx.stub:
        _stub_copy(spec.fixture_outputs["aggregated"], out_tsv)
        return out_tsv
    if not ctx.fusion_tsv:
        raise ValueError("pvacfuse requires fusion_tsv (AGFusion dir or Arriba TSV) or tools.stub=true")
    if not ctx.hla_alleles:
        raise ValueError("pvacfuse requires hla_alleles")
    work = (ctx.outdir / "tools" / "pvacfuse").resolve()
    hla = ",".join(ctx.hla_alleles)
    algos = ctx.prediction_algorithms.replace(",", " ").split()
    fusion = ctx.fusion_tsv.resolve() if ctx.fusion_tsv else ctx.fusion_tsv
    cmd = [
        ctx.exe("pvacfuse"), "run",
        str(fusion),
        ctx.sample_id,
        hla,
        *algos,
        str(work),
        "-e1", "8,9,10,11",
    ]
    _run_cmd(cmd, work)
    candidates = sorted(work.rglob("*aggregated*.tsv"))
    if not candidates:
        candidates = sorted(work.rglob("*.tsv"))
    if not candidates:
        raise FileNotFoundError(f"pvacfuse finished but no aggregated TSV in {work}")
    shutil.copy2(candidates[0], out_tsv)
    return out_tsv


def run_pvacsplice(ctx: RunContext, out_tsv: Path) -> Path:
    spec = TOOL_REGISTRY["pvacsplice"]
    if ctx.stub:
        _stub_copy(spec.fixture_outputs["aggregated"], out_tsv)
        return out_tsv
    if not ctx.splice_junction_tsv:
        raise ValueError("pvacsplice requires splice_junction_tsv (RegTools annotated junctions)")
    vcf = ctx.variants_vcf or ctx.tumor_vcf
    if not vcf:
        raise ValueError("pvacsplice requires variants_vcf/tumor_vcf or tools.stub=true")
    ref_fasta = ctx.reference_fasta
    if not ref_fasta:
        import os
        ref_fasta = Path(os.environ.get("NEOAG_REFERENCE_FASTA", "")) if os.environ.get("NEOAG_REFERENCE_FASTA") else None
    gtf = ctx.gencode_gtf
    if not ref_fasta or not Path(ref_fasta).is_file():
        raise ValueError("pvacsplice requires reference_fasta (inputs.reference_fasta or NEOAG_REFERENCE_FASTA)")
    if not gtf or not Path(gtf).is_file():
        raise ValueError("pvacsplice requires gencode_gtf (inputs.gencode_gtf)")
    if not ctx.hla_alleles:
        raise ValueError("pvacsplice requires hla_alleles")
    work = (ctx.outdir / "tools" / "pvacsplice").resolve()
    hla = ",".join(ctx.hla_alleles)
    sample_name = ctx.tumor_sample_name or ctx.sample_id
    algos = ctx.prediction_algorithms.replace(",", " ").split()
    cmd = [
        ctx.exe("pvacsplice"), "run",
        str(ctx.splice_junction_tsv.resolve()),
        sample_name,
        hla,
        *algos,
        str(work),
        str(Path(vcf).resolve()),
        str(Path(ref_fasta).resolve()),
        str(Path(gtf).resolve()),
        "-e1", "8,9,10,11",
    ]
    if ctx.pass_only:
        cmd.append("--pass-only")
    if ctx.normal_sample_name:
        cmd.extend(["--normal-sample-name", ctx.normal_sample_name])
    junction_score = getattr(ctx, "junction_score", None)
    if junction_score is not None:
        cmd.extend(["-j", str(junction_score)])
    _run_cmd(cmd, work)
    candidates = sorted(work.rglob("*aggregated*.tsv"))
    if not candidates:
        candidates = sorted(work.rglob("*.tsv"))
    if not candidates:
        raise FileNotFoundError(f"pvacsplice finished but no aggregated TSV in {work}")
    shutil.copy2(candidates[0], out_tsv)
    return out_tsv


def run_vep_appm(ctx: RunContext, out_tsv: Path) -> Path:
    spec = TOOL_REGISTRY["vep"]
    if ctx.stub:
        _stub_copy(spec.fixture_outputs["appm"], out_tsv)
        return out_tsv
    if not ctx.variants_vcf:
        raise ValueError("vep requires variants_vcf or tools.stub=true")
    work = ctx.outdir / "tools" / "vep"
    raw = work / "vep_raw.tsv"
    import os
    docker_mode = _effective_mode("vep") == RunnerMode.DOCKER
    # In Docker mode use the bare executable name (inside the container);
    # in conda mode resolve to the full host path (may include conda env prefix).
    vep_bin = "vep" if docker_mode else ctx.exe("vep")
    cmd = [
        *vep_bin.split(),
        "-i", str(ctx.variants_vcf),
        "--symbol", "--tab", "--force_overwrite",
        "-o", str(raw),
    ]
    if os.environ.get("NEOAG_VEP_ONLINE", "").lower() in {"1", "true", "yes"}:
        cmd.extend(["--database", "--species", "homo_sapiens"])
    else:
        cmd.extend(["--cache", "--offline"])
        cache_dir = os.environ.get("NEOAG_VEP_CACHE", "").strip()
        if cache_dir:
            if docker_mode:
                # VEP official container expects cache at /opt/vep/.vep
                cmd.extend(["--dir_cache", "/opt/vep/.vep"])
            else:
                cmd.extend(["--dir_cache", cache_dir])
        cache_version = os.environ.get("NEOAG_VEP_CACHE_VERSION", "").strip()
        if cache_version:
            cmd.extend(["--cache_version", cache_version])
        reference_fasta = os.environ.get("NEOAG_REFERENCE_FASTA", "").strip()
        if reference_fasta and Path(reference_fasta).is_file():
            cmd.extend(["--fasta", reference_fasta])
    _run_cmd(cmd, work)
    vep_to_appm_tsv(raw, out_tsv)
    return out_tsv


def run_lohhla(ctx: RunContext, out_tsv: Path) -> Path:
    import os

    spec = TOOL_REGISTRY["lohhla"]
    if ctx.stub:
        _stub_copy(spec.fixture_outputs["hla_loh"], out_tsv)
        return out_tsv
    pred_path = (
        ctx.lohhla_prediction
        or _path_or_none(ctx.executables.get("lohhla_prediction"))
        or _path_or_none(os.environ.get("NEOAG_LOHHLA_PREDICTION", ""))
    )
    if pred_path and pred_path.is_file():
        lohhla_to_hla_loh_tsv(pred_path, out_tsv)
        return out_tsv
    raise NotImplementedError(
        "LOHHLA requires lohhla_prediction file (HLAlossPrediction_CI) or tools.stub=true. "
        "Set inputs.lohhla_prediction in run config or NEOAG_LOHHLA_PREDICTION."
    )


def _path_or_none(val: str | Path | None) -> Path | None:
    if not val:
        return None
    p = Path(val)
    return p if p.is_file() else None


def run_facets(ctx: RunContext, out_tsv: Path) -> Path:
    spec = TOOL_REGISTRY["facets"]
    if ctx.stub:
        _stub_copy(spec.fixture_outputs["purity"], out_tsv)
        return out_tsv
    if not ctx.facets_rds:
        raise ValueError("facets requires facets_rds output or tools.stub=true")
    work = ctx.outdir / "tools" / "facets"
    work.mkdir(parents=True, exist_ok=True)
    raw = work / "facets_purity.txt"
    cncf = work / "facets_cncf.tsv"
    facets_r = Path(ctx.exe("facets"))
    if not facets_r.is_absolute():
        facets_r = (ROOT / facets_r).resolve()
    cmd = ["Rscript", str(facets_r), str(ctx.facets_rds.resolve()), str(raw.resolve()), str(cncf.resolve())]
    _run_cmd(cmd, work)
    facets_to_purity_tsv(raw, ctx.sample_id, out_tsv)
    return out_tsv


def _write_prime_stub(pairs: list[tuple[str, str]], out_tsv: Path, sample_id: str) -> None:
    from ..adapters.prime import predict_pairs, write_prime_evidence

    write_prime_evidence(out_tsv, [
        {**row, "sample_id": sample_id, "source_file": "stub"}
        for row in predict_pairs(pairs)
    ])


def _resolve_prime_exe(ctx: RunContext) -> Path:
    import os

    prime_exe = Path(ctx.exe("prime"))
    if not prime_exe.is_absolute():
        found = shutil.which(str(prime_exe))
        if found:
            prime_exe = Path(found)
        elif (ROOT / prime_exe).is_file():
            prime_exe = (ROOT / prime_exe).resolve()
        else:
            env_prime = os.environ.get("NEOAG_PRIME_BIN", "")
            if env_prime and Path(env_prime).is_file():
                prime_exe = Path(env_prime)
    return prime_exe


def _run_prime_batch(
    batch_dir: Path,
    peptides: list[str],
    prime_alleles: list[str],
    prime_exe: Path,
    mixmhcpred: str | None,
) -> Path:
    import os
    import subprocess

    batch_dir.mkdir(parents=True, exist_ok=True)
    pep_file = batch_dir / "peptides.txt"
    with pep_file.open("w", encoding="utf-8") as fh:
        for peptide in peptides:
            fh.write(peptide + "\n")
    raw_out = batch_dir / "prime_out.tsv"
    cmd = [
        str(prime_exe),
        "-i",
        str(pep_file),
        "-o",
        str(raw_out),
        "-a",
        ",".join(prime_alleles),
    ]
    if mixmhcpred:
        cmd.extend(["-mix", mixmhcpred])
    _run_cmd(cmd, batch_dir, env=os.environ.copy())
    return raw_out


def _prime_rows_by_peptide(raw_paths: list[Path]) -> dict[str, dict[str, str]]:
    from ..adapters.prime import read_prime_wide_rows

    by_pep: dict[str, dict[str, str]] = {}
    for raw_out in raw_paths:
        for row in read_prime_wide_rows(raw_out):
            peptide = (row.get("Peptide") or row.get("peptide") or "").strip().upper()
            if peptide:
                by_pep[peptide] = row
    return by_pep


def _prime_evidence_rows(
    pairs: list[tuple[str, str]],
    by_pep: dict[str, dict[str, str]],
    sample_id: str,
    source_file: str,
) -> list[dict[str, str]]:
    from ..adapters.prime import extract_prime_pair, predict_pair_stub

    merged: list[dict[str, str]] = []
    for peptide, allele in pairs:
        row = by_pep.get(peptide.upper(), {})
        score, rank = extract_prime_pair(row, peptide, allele)
        if not score and not rank:
            score, rank = predict_pair_stub(peptide, allele)
        merged.append({
            "sample_id": sample_id,
            "peptide": peptide,
            "hla_allele": allele,
            "prime_score": score,
            "prime_rank": rank,
            "source_file": source_file,
        })
    return merged


def _run_prime_external(pairs: list[tuple[str, str]], out_tsv: Path, ctx: RunContext) -> None:
    import os
    from concurrent.futures import ThreadPoolExecutor, as_completed

    from ..adapters.prime import (
        prime_allele_tag,
        prime_parallel_jobs,
        split_peptide_chunks,
        write_prime_evidence,
    )

    work = out_tsv.parent / "prime"
    work.mkdir(parents=True, exist_ok=True)
    unique_peptides = sorted({peptide.upper() for peptide, _ in pairs})
    prime_alleles = sorted({prime_allele_tag(allele) for _, allele in pairs})
    prime_exe = _resolve_prime_exe(ctx)
    mixmhcpred = ctx.executables.get("mixmhcpred") or os.environ.get("MIXMHCPRED_BIN")
    n_jobs = prime_parallel_jobs()
    chunks = split_peptide_chunks(unique_peptides, n_jobs)

    if len(chunks) == 1:
        raw_out = _run_prime_batch(work, chunks[0], prime_alleles, prime_exe, mixmhcpred)
        by_pep = _prime_rows_by_peptide([raw_out])
        source_file = str(raw_out)
    else:
        print(
            f"prime: running {len(chunks)} batches in parallel "
            f"(NEOAG_PRIME_JOBS={n_jobs}, {len(unique_peptides)} peptides)",
            flush=True,
        )
        batch_dirs = [work / f"batch_{i:03d}" for i in range(len(chunks))]
        raw_paths: list[Path] = []
        with ThreadPoolExecutor(max_workers=len(chunks)) as pool:
            futures = {
                pool.submit(
                    _run_prime_batch,
                    batch_dir,
                    chunk,
                    prime_alleles,
                    prime_exe,
                    mixmhcpred,
                ): batch_dir
                for batch_dir, chunk in zip(batch_dirs, chunks)
            }
            for fut in as_completed(futures):
                raw_paths.append(fut.result())
        by_pep = _prime_rows_by_peptide(raw_paths)
        source_file = str(work)

    write_prime_evidence(out_tsv, _prime_evidence_rows(pairs, by_pep, ctx.sample_id, source_file))


def run_prime(ctx: RunContext, out_tsv: Path) -> Path:
    from ..tools.prep import unique_peptide_hla_pairs

    if not ctx.raw_peptides:
        raise ValueError("prime requires raw_peptides.tsv")
    pairs = unique_peptide_hla_pairs(ctx.raw_peptides)
    if not pairs:
        raise ValueError(f"prime: no peptide/HLA pairs in {ctx.raw_peptides}")
    out_tsv = out_tsv.resolve()
    if ctx.stub:
        require_non_strict("stub PRIME output")
        _write_prime_stub(pairs, out_tsv, ctx.sample_id)
        return out_tsv
    import shutil

    exe = ctx.exe("prime")
    if shutil.which(exe) or Path(exe).is_file():
        _run_prime_external(pairs, out_tsv, ctx)
        return out_tsv
    require_non_strict("PRIME fallback to stub when executable is unavailable")
    _write_prime_stub(pairs, out_tsv, ctx.sample_id)
    return out_tsv


def _write_bigmhc_im_stub(pairs: list[tuple[str, str]], out_tsv: Path, sample_id: str) -> None:
    from ..adapters.bigmhc_im import predict_pairs, write_bigmhc_im_evidence

    write_bigmhc_im_evidence(out_tsv, [
        {**row, "sample_id": sample_id, "source_file": "stub"}
        for row in predict_pairs(pairs)
    ])


def _run_bigmhc_im_external(pairs: list[tuple[str, str]], out_tsv: Path, ctx: RunContext) -> None:
    import os
    import subprocess

    work = out_tsv.parent / "bigmhc_im"
    work.mkdir(parents=True, exist_ok=True)
    in_csv = work / "input.csv"
    with in_csv.open("w", encoding="utf-8", newline="") as fh:
        fh.write("mhc,pep\n")
        for peptide, allele in pairs:
            fh.write(f"{allele},{peptide}\n")
    bigmhc_dir = Path(ctx.executables.get("bigmhc_dir") or os.environ.get("BIGMHC_DIR", ""))
    predict_py = bigmhc_dir / "src" / "predict.py"
    out_prd = work / "input.csv.prd"
    cmd = [
        "python3",
        str(predict_py),
        f"-i={in_csv}",
        "-m=im",
        "-a=0",
        "-p=1",
        "-c=1",
        f"-o={out_prd}",
        "-d=cpu",
    ]
    _run_cmd(cmd, bigmhc_dir / "src", env=os.environ.copy())
    from ..adapters.bigmhc_im import parse_bigmhc_im, write_bigmhc_im_evidence

    write_bigmhc_im_evidence(out_tsv, parse_bigmhc_im(out_prd, ctx.sample_id))


def run_bigmhc_im(ctx: RunContext, out_tsv: Path) -> Path:
    import os
    from ..tools.prep import unique_peptide_hla_pairs

    if not ctx.raw_peptides:
        raise ValueError("bigmhc_im requires raw_peptides.tsv")
    pairs = unique_peptide_hla_pairs(ctx.raw_peptides)
    if not pairs:
        raise ValueError(f"bigmhc_im: no peptide/HLA pairs in {ctx.raw_peptides}")
    out_tsv = out_tsv.resolve()
    if ctx.stub:
        require_non_strict("stub BigMHC_IM output")
        _write_bigmhc_im_stub(pairs, out_tsv, ctx.sample_id)
        return out_tsv
    bigmhc_dir = Path(ctx.executables.get("bigmhc_dir") or os.environ.get("BIGMHC_DIR", ""))
    predict_py = bigmhc_dir / "src" / "predict.py"
    if predict_py.is_file():
        _run_bigmhc_im_external(pairs, out_tsv, ctx)
        return out_tsv
    require_non_strict("BigMHC_IM fallback to stub when predict.py is unavailable")
    _write_bigmhc_im_stub(pairs, out_tsv, ctx.sample_id)
    return out_tsv


def _write_deepimmuno_stub(pairs: list[tuple[str, str]], out_tsv: Path, sample_id: str) -> None:
    from ..adapters.deepimmuno import predict_pairs, write_deepimmuno_evidence

    write_deepimmuno_evidence(
        out_tsv,
        [{**row, "sample_id": sample_id, "source_file": "stub"} for row in predict_pairs(pairs)],
    )


def _run_deepimmuno_external(pairs: list[tuple[str, str]], out_tsv: Path, ctx: RunContext) -> None:
    import os

    from ..adapters.deepimmuno import resolve_deepimmuno_dir, run_deepimmuno_batch, write_deepimmuno_evidence

    custom = ctx.executables.get("deepimmuno_dir") or os.environ.get("DEEPIMMUNO_DIR", "")
    deep_dir = resolve_deepimmuno_dir(custom or None)
    rows = run_deepimmuno_batch(pairs, deep_dir, sample_id=ctx.sample_id)
    write_deepimmuno_evidence(out_tsv, rows)


def run_deepimmuno(ctx: RunContext, out_tsv: Path) -> Path:
    import os

    from ..adapters.deepimmuno import resolve_deepimmuno_dir
    from ..tools.prep import unique_peptide_hla_pairs

    if not ctx.raw_peptides:
        raise ValueError("deepimmuno requires raw_peptides.tsv")
    pairs = unique_peptide_hla_pairs(ctx.raw_peptides)
    if not pairs:
        raise ValueError(f"deepimmuno: no peptide/HLA pairs in {ctx.raw_peptides}")
    out_tsv = out_tsv.resolve()
    if ctx.stub:
        require_non_strict("stub DeepImmuno output")
        _write_deepimmuno_stub(pairs, out_tsv, ctx.sample_id)
        return out_tsv
    try:
        custom = ctx.executables.get("deepimmuno_dir") or os.environ.get("DEEPIMMUNO_DIR", "")
        resolve_deepimmuno_dir(custom or None)
        _run_deepimmuno_external(pairs, out_tsv, ctx)
        return out_tsv
    except FileNotFoundError:
        require_non_strict("DeepImmuno fallback to stub when model directory is unavailable")
        _write_deepimmuno_stub(pairs, out_tsv, ctx.sample_id)
        return out_tsv


def run_optitype(ctx: RunContext, out_csv: Path) -> Path:
    """Run OptiType HLA typing from BAM or FASTQ input.

    OptiType produces a *_result.tsv with 4-digit HLA typing (A1,A2,B1,B2,C1,C2).
    """
    import shutil

    spec = TOOL_REGISTRY["optitype"]
    if ctx.stub:
        _stub_copy(spec.fixture_outputs["hla_csv"], out_csv)
        return out_csv

    # Determine input files: BAM takes precedence, then FASTQ pairs
    input_bam = ctx.variants_vcf  # not correct for BAM, but ctx doesn't have bam_path yet
    # Use environment variables or executables config for BAM paths
    import os
    bam = (os.environ.get("NEOAG_OPTITYPE_BAM") or ""
           or str(ctx.executables.get("optitype_bam", "")))
    fq1 = (os.environ.get("NEOAG_OPTITYPE_FQ1") or ""
           or str(ctx.executables.get("optitype_fq1", "")))
    fq2 = (os.environ.get("NEOAG_OPTITYPE_FQ2") or ""
           or str(ctx.executables.get("optitype_fq2", "")))
    seq_type = os.environ.get("NEOAG_OPTITYPE_SEQ_TYPE", "dna").strip().lower()

    if not bam and not fq1:
        raise RuntimeError(
            "OptiType requires input BAM or FASTQ. "
            "Set NEOAG_OPTITYPE_BAM or NEOAG_OPTITYPE_FQ1/FQ2 env vars, "
            "or use tools.stub=true."
        )

    work = ctx.outdir / "tools" / "optitype_work"
    work.mkdir(parents=True, exist_ok=True)

    cmd = [ctx.exe("optitype"), "run", "-o", str(work.resolve()), "-p", ctx.sample_id]

    if bam:
        cmd.extend(["-i", bam])
    else:
        cmd.extend(["-i", fq1])
        if fq2:
            cmd.extend(["-i", fq2])

    if seq_type == "rna":
        cmd.append("--rna")
    else:
        cmd.append("--dna")

    threads = os.environ.get("NEOAG_OPTITYPE_THREADS", "4")
    cmd.extend(["--threads", threads])

    _run_cmd(cmd, work)

    # Find the result CSV
    result_files = list(work.glob(f"{ctx.sample_id}*_result.tsv"))
    if not result_files:
        result_files = list(work.glob("*_result.tsv"))
    if not result_files:
        raise FileNotFoundError(f"OptiType did not produce a _result.tsv in {work}")

    result_csv = result_files[0]
    shutil.copy2(result_csv, out_csv)
    return out_csv


RUNNERS = {
    "netmhcpan": lambda ctx, p: run_netmhcpan(ctx, p),
    "netmhcstabpan": lambda ctx, p: run_netmhcstabpan(ctx, p),
    "prime": lambda ctx, p: run_prime(ctx, p),
    "bigmhc_im": lambda ctx, p: run_bigmhc_im(ctx, p),
    "deepimmuno": lambda ctx, p: run_deepimmuno(ctx, p),
    "mhcflurry": lambda ctx, p: run_mhcflurry(ctx, p),
    "pvacseq": lambda ctx, p: run_pvacseq(ctx, p),
    "pvacfuse": lambda ctx, p: run_pvacfuse(ctx, p),
    "pvacsplice": lambda ctx, p: run_pvacsplice(ctx, p),
    "vep": lambda ctx, p: run_vep_appm(ctx, p),
    "lohhla": lambda ctx, p: run_lohhla(ctx, p),
    "facets": lambda ctx, p: run_facets(ctx, p),
    "optitype": lambda ctx, p: run_optitype(ctx, p),
}


def run_tool(name: str, ctx: RunContext, output: str | Path) -> Path:
    if name not in RUNNERS:
        raise KeyError(f"No runner for tool '{name}'. Known: {sorted(RUNNERS)}")
    if ctx.stub:
        require_non_strict(f"stub mode for tool {name}")

    mode = _effective_mode(name)
    if mode == RunnerMode.DOCKER:
        # In Docker mode the executable lives inside the container;
        # skip the host-PATH existence check.
        image = tool_docker_image(name)
        print(f"[runner] {name}: docker mode ({image})", flush=True)
    else:
        st = check_tool(name, ctx.exe(name))
        if not st.available and not ctx.stub:
            raise RuntimeError(
                f"{name} ({ctx.exe(name)}) not available: {st.message}. "
                "Install the tool, set executables.<name> in run config, or enable tools.stub."
            )
        if resolve_runner_mode() == RunnerMode.DOCKER and not tool_docker_image(name):
            print(f"[runner] {name}: no docker image — falling back to conda", flush=True)

    with _tool_context(name):
        return RUNNERS[name](ctx, Path(output))
