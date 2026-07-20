from __future__ import annotations

import argparse
import json
import os
import threading
import uuid
from dataclasses import dataclass, field
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Callable
from urllib.parse import urlparse

from .audit import AuditLogger
from .doctor import run_doctor
from .io_utils import ensure_dir, now_iso, read_json, write_json
from .pipeline_runner import run_pipeline_full

try:
    from neoag.agent_skills.appm_review import main as appm_review_main
    from neoag.agent_skills.ccf_review import main as ccf_review_main
    from neoag.agent_skills.input_qc import main as input_qc_main
    from neoag.agent_skills.patient_report import main as patient_report_main
    from neoag.agent_skills.ranking_compare import main as ranking_compare_main
except Exception:  # pragma: no cover - import safety
    input_qc_main = ranking_compare_main = appm_review_main = ccf_review_main = patient_report_main = None

HIGH_RISK_FLAGS = {"execute", "hpc_submit", "install_tools", "overwrite", "allow_overwrite", "delete", "delete_files"}
WRITE_PATH_FIELDS = {"outdir", "output", "output_dir", "result_dir"}
READ_PATH_HINTS = ("file", "path", "dir", "manifest", "recommendation", "netmhcpan42", "evidence_report", "ranked_peptides", "ccf", "purity")


@dataclass(frozen=True)
class RouteSpec:
    route: str
    risk_level: str
    required: tuple[str, ...] = ()
    optional: tuple[str, ...] = ()
    high_risk_when: tuple[str, ...] = ()
    description: str = ""


ROUTE_SPECS: dict[str, RouteSpec] = {
    "/doctor": RouteSpec(
        route="/doctor",
        risk_level="LOW",
        optional=("project_root", "outdir", "tools_manifest", "reference_manifest", "sample_manifest", "profile", "run_demo", "run_pytest", "run_nextflow", "mini_smoke", "skip_release_audit", "execute", "approved"),
        high_risk_when=("execute",),
        description="Read-only health checks by default; execute=true may run optional smoke tests.",
    ),
    "/input-qc": RouteSpec(
        route="/input-qc",
        risk_level="LOW",
        optional=("outdir", "manifest", "result_dir"),
        description="Input and output directory QC.",
    ),
    "/ranking-compare": RouteSpec(
        route="/ranking-compare",
        risk_level="LOW",
        required=("recommendation", "netmhcpan42"),
        optional=("outdir",),
        description="Compare recommendation ranking with NetMHCpan ranking.",
    ),
    "/appm-review": RouteSpec(
        route="/appm-review",
        risk_level="LOW",
        optional=("outdir", "evidence_report", "appm_summary", "appm_gene_status", "appm_submodule_scores", "peptide_escape_flags", "ranked_peptides", "hla_loh"),
        description="APPM and immune-escape review.",
    ),
    "/ccf-review": RouteSpec(
        route="/ccf-review",
        risk_level="LOW",
        optional=("outdir", "ccf", "purity", "ranked_peptides"),
        description="CCF and clonality review.",
    ),
    "/patient-report": RouteSpec(
        route="/patient-report",
        risk_level="LOW",
        optional=("outdir", "evidence_report", "recommendation", "netmhcpan42", "ranking_compare_report", "appm_review", "ccf_review"),
        description="Patient-facing draft report output.",
    ),
    "/pipeline-full": RouteSpec(
        route="/pipeline-full",
        risk_level="LOW",
        required=("sample_manifest",),
        optional=("outdir", "tools_manifest", "reference_manifest", "project_root", "profile", "execute", "strict", "approved", "overwrite", "hpc_submit"),
        high_risk_when=("execute", "overwrite", "hpc_submit"),
        description="Manifest-driven full pipeline; dry-run by default.",
    ),
}


class GatewayError(Exception):
    def __init__(self, status: str, message: str, http_code: int = 400, details: dict[str, Any] | None = None):
        super().__init__(message)
        self.status = status
        self.message = message
        self.http_code = http_code
        self.details = details or {}


class JobStore:
    def __init__(self, outdir: Path):
        self.outdir = ensure_dir(outdir)
        self.jobs: dict[str, dict[str, Any]] = {}
        self.lock = threading.Lock()

    def add(self, payload: dict[str, Any]) -> str:
        job_id = payload.get("job_id") or f"job-{uuid.uuid4().hex[:10]}"
        with self.lock:
            self.jobs[job_id] = {"job_id": job_id, "created_at": now_iso(), **payload}
        write_json(self.outdir / f"{job_id}.json", self.jobs[job_id])
        return job_id

    def update(self, job_id: str, **kwargs: Any) -> None:
        with self.lock:
            self.jobs.setdefault(job_id, {"job_id": job_id})
            self.jobs[job_id].update(kwargs)
            self.jobs[job_id]["updated_at"] = now_iso()
        write_json(self.outdir / f"{job_id}.json", self.jobs[job_id])

    def get(self, job_id: str) -> dict[str, Any] | None:
        with self.lock:
            if job_id in self.jobs:
                return self.jobs[job_id]
        p = self.outdir / f"{job_id}.json"
        if p.is_file():
            try:
                return read_json(p)
            except Exception:
                return None
        return None


def _truthy(data: dict[str, Any], key: str) -> bool:
    value = data.get(key)
    if isinstance(value, str):
        return value.lower() in {"1", "true", "yes", "y"}
    return bool(value)


def _risk(route: str, data: dict[str, Any]) -> str:
    spec = ROUTE_SPECS.get(route)
    if not spec:
        return "LOW"
    if any(_truthy(data, flag) for flag in set(spec.high_risk_when) | HIGH_RISK_FLAGS):
        return "HIGH"
    return spec.risk_level


def _approval_required(route: str, data: dict[str, Any]) -> bool:
    return _risk(route, data) == "HIGH" and not _truthy(data, "approved")


def _is_relative_to(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.resolve())
        return True
    except Exception:
        return False


def _safe_write_path(path: str | Path, allowed_roots: list[Path]) -> bool:
    p = Path(path).expanduser()
    if not p.is_absolute():
        return True
    return any(_is_relative_to(p, root) for root in allowed_roots)


def _looks_like_path_key(key: str) -> bool:
    low = key.lower()
    return any(h in low for h in READ_PATH_HINTS)


def _validate_request(route: str, data: dict[str, Any], allowed_roots: list[Path]) -> None:
    spec = ROUTE_SPECS.get(route)
    if spec is None:
        raise GatewayError("NOT_FOUND", f"unknown route: {route}", 404)
    if not isinstance(data, dict):
        raise GatewayError("BAD_REQUEST", "request body must be a JSON object", 400)
    missing = [k for k in spec.required if data.get(k) in {None, ""}]
    if missing:
        raise GatewayError("BAD_REQUEST", "missing required fields", 400, {"missing": missing})
    allowed = set(spec.required) | set(spec.optional) | {"job_id", "approved", "wait", "request_id", "case_id", "sample_id", "user", "input_manifest"} | HIGH_RISK_FLAGS
    unknown = sorted(k for k in data if k not in allowed)
    if unknown:
        raise GatewayError("BAD_REQUEST", "unknown fields", 400, {"unknown": unknown, "allowed": sorted(allowed)})
    if any(_truthy(data, flag) for flag in {"install_tools", "delete", "delete_files"}):
        raise GatewayError("UNSUPPORTED_HIGH_RISK_ACTION", "install/delete actions are not implemented by this gateway", 400)
    for key, value in data.items():
        if key in WRITE_PATH_FIELDS and isinstance(value, str) and not _safe_write_path(value, allowed_roots):
            raise GatewayError("PATH_NOT_ALLOWED", f"write path is outside allowed roots: {key}", 400, {"path": value})
        if _looks_like_path_key(key) and isinstance(value, str) and "\x00" in value:
            raise GatewayError("BAD_REQUEST", f"invalid path value: {key}", 400)


def _run_skill_func(func: Callable[[list[str]], Any] | None, argv: list[str]) -> tuple[str, int]:
    if func is None:
        return "skill entrypoint unavailable", 2
    try:
        rc = func(argv)
        return "", int(rc or 0)
    except SystemExit as e:
        return "", int(e.code or 0)
    except Exception as exc:
        return str(exc), 2


def _request_metadata(route: str, data: dict[str, Any], *, job_id: str, risk: str, status: str, output_dir: str = "", failure_reason: str = "") -> dict[str, Any]:
    input_manifest = data.get("sample_manifest") or data.get("manifest") or data.get("input_manifest") or ""
    return {
        "request_id": data.get("request_id") or job_id,
        "case_id": data.get("case_id") or data.get("sample_id") or "",
        "sample_id": data.get("sample_id") or data.get("case_id") or "",
        "user": data.get("user") or os.environ.get("USER", ""),
        "timestamp": now_iso(),
        "input_manifest": input_manifest,
        "skill_or_pipeline": route.lstrip("/"),
        "command_preview": f"neoag-gateway {route}",
        "risk_level": risk,
        "approval_status": "approved" if _truthy(data, "approved") else ("required" if risk == "HIGH" else "not_required"),
        "job_id": job_id,
        "output_dir": output_dir,
        "status": status,
        "failure_reason": failure_reason,
    }


class GatewayHandler(BaseHTTPRequestHandler):
    server_version = "NeoAgGateway/0.2"

    def _json(self, code: int, payload: dict[str, Any]) -> None:
        body = json.dumps(payload, ensure_ascii=False, indent=2).encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _read_body(self) -> dict[str, Any]:
        length = int(self.headers.get("Content-Length", "0") or "0")
        if length <= 0:
            return {}
        raw = self.rfile.read(length).decode("utf-8")
        try:
            data = json.loads(raw) if raw else {}
        except Exception:
            raise GatewayError("BAD_JSON", "request body is not valid JSON", 400)
        if not isinstance(data, dict):
            raise GatewayError("BAD_REQUEST", "request body must be a JSON object", 400)
        return data

    @property
    def store(self) -> JobStore:
        return self.server.job_store  # type: ignore[attr-defined]

    @property
    def logger(self) -> AuditLogger:
        return self.server.audit_logger  # type: ignore[attr-defined]

    @property
    def project_root(self) -> Path:
        return self.server.project_root  # type: ignore[attr-defined]

    @property
    def allowed_roots(self) -> list[Path]:
        return self.server.allowed_roots  # type: ignore[attr-defined]

    def log_message(self, fmt: str, *args: Any) -> None:  # keep stdout clean for wrappers
        return

    def do_GET(self) -> None:
        parsed = urlparse(self.path)
        if parsed.path == "/health":
            self._json(200, {"status": "OK", "service": "neoag-gateway", "time": now_iso(), "routes": sorted(ROUTE_SPECS)})
            return
        if parsed.path.startswith("/job-status/"):
            job_id = parsed.path.rstrip("/").split("/")[-1]
            job = self.store.get(job_id)
            if job:
                self._json(200, job)
            else:
                self._json(404, {"status": "NOT_FOUND", "job_id": job_id})
            return
        self._json(404, {"status": "NOT_FOUND", "path": parsed.path})

    def do_POST(self) -> None:
        parsed = urlparse(self.path)
        route = parsed.path.rstrip("/") or "/"
        try:
            data = self._read_body()
            _validate_request(route, data, self.allowed_roots)
        except GatewayError as exc:
            self._json(exc.http_code, {"status": exc.status, "message": exc.message, **exc.details})
            return
        except Exception as exc:
            self._json(500, {"status": "FAIL", "message": str(exc)})
            return

        try:
            self._enqueue_or_reject(route, data)
        except GatewayError as exc:
            self._json(exc.http_code, {"status": exc.status, "message": exc.message, **exc.details})
        except Exception as exc:
            self._json(500, {"status": "FAIL", "message": str(exc)})

    def _enqueue_or_reject(self, route: str, data: dict[str, Any]) -> None:
        risk = _risk(route, data)
        if _approval_required(route, data):
            job_id = self.store.add({"route": route, "status": "APPROVAL_REQUIRED", "risk_level": risk, "request": data})
            meta = _request_metadata(route, data, job_id=job_id, risk=risk, status="APPROVAL_REQUIRED", failure_reason="approved=true required for HIGH risk request")
            self.store.update(job_id, **{k: v for k, v in meta.items() if k != "job_id"})
            self.logger.log("gateway.request", "APPROVAL_REQUIRED", route, case_id=meta["case_id"], request_id=meta["request_id"], risk_level=risk, approval_status=meta["approval_status"], command_preview=meta["command_preview"], metadata={"request": data, **meta})
            self._json(403, {"job_id": job_id, "status": "APPROVAL_REQUIRED", "route": route, "risk_level": risk, "message": "Set approved=true after reviewing dry-run plan."})
            return

        job_id = self.store.add({"route": route, "status": "QUEUED", "risk_level": risk, "request": data})
        output_dir = str(self._planned_outdir(data, job_id))
        meta = _request_metadata(route, data, job_id=job_id, risk=risk, status="QUEUED", output_dir=output_dir)
        self.store.update(job_id, **{k: v for k, v in meta.items() if k != "job_id"})
        self.logger.log("gateway.request", "QUEUED", route, case_id=meta["case_id"], request_id=meta["request_id"], risk_level=risk, approval_status=meta["approval_status"], command_preview=meta["command_preview"], metadata={"request": data, **meta})
        thread = threading.Thread(target=self._run_job, args=(route, data, job_id, risk), daemon=True)
        thread.start()
        self._json(202, {"job_id": job_id, "status": "QUEUED", "route": route, "risk_level": risk, "status_url": f"/job-status/{job_id}"})

    def _run_job(self, route: str, data: dict[str, Any], job_id: str, risk: str) -> None:
        start_meta = _request_metadata(route, data, job_id=job_id, risk=risk, status="RUNNING", output_dir=str(self._planned_outdir(data, job_id)))
        self.store.update(job_id, started_at=now_iso(), **{k: v for k, v in start_meta.items() if k != "job_id"})
        self.logger.log("gateway.job", "START", route, case_id=start_meta["case_id"], request_id=start_meta["request_id"], risk_level=risk, approval_status=start_meta["approval_status"], command_preview=start_meta["command_preview"], metadata=start_meta)
        try:
            result = self._dispatch(route, data, job_id)
            status = result.get("status", "PASS")
            finish_meta = _request_metadata(route, data, job_id=job_id, risk=risk, status=status, output_dir=str(self._planned_outdir(data, job_id)), failure_reason=result.get("error", "") or result.get("failure_reason", ""))
            self.store.update(job_id, finished_at=now_iso(), result=result, **{k: v for k, v in finish_meta.items() if k != "job_id"})
            self.logger.log("gateway.job", status, route, case_id=finish_meta["case_id"], request_id=finish_meta["request_id"], risk_level=risk, approval_status=finish_meta["approval_status"], command_preview=finish_meta["command_preview"], metadata={"result": result, **finish_meta})
        except Exception as exc:
            payload = {"status": "FAIL", "error": str(exc)}
            fail_meta = _request_metadata(route, data, job_id=job_id, risk=risk, status="FAIL", output_dir=str(self._planned_outdir(data, job_id)), failure_reason=str(exc))
            self.store.update(job_id, finished_at=now_iso(), result=payload, **{k: v for k, v in fail_meta.items() if k != "job_id"})
            self.logger.log("gateway.job", "FAIL", route, case_id=fail_meta["case_id"], request_id=fail_meta["request_id"], risk_level=risk, approval_status=fail_meta["approval_status"], command_preview=fail_meta["command_preview"], metadata={"result": payload, **fail_meta})

    def _planned_outdir(self, data: dict[str, Any], job_id: str) -> Path:
        raw = data.get("outdir") or str(self.store.outdir / job_id)
        p = Path(raw)
        if not p.is_absolute():
            p = self.project_root / p
        return p

    def _job_outdir(self, data: dict[str, Any], job_id: str) -> Path:
        raw = data.get("outdir") or str(self.store.outdir / job_id)
        p = Path(raw)
        if not p.is_absolute():
            p = self.project_root / p
        ensure_dir(p)
        return p

    def _dispatch(self, route: str, data: dict[str, Any], job_id: str) -> dict[str, Any]:
        outdir = self._job_outdir(data, job_id)
        if route == "/doctor":
            res = run_doctor(
                project_root=data.get("project_root") or self.project_root,
                outdir=outdir,
                tools_manifest=data.get("tools_manifest"),
                reference_manifest=data.get("reference_manifest"),
                sample_manifest=data.get("sample_manifest"),
                profile=data.get("profile", "local"),
                run_demo=bool(data.get("run_demo")),
                run_pytest=bool(data.get("run_pytest")),
                run_nextflow=bool(data.get("run_nextflow")),
                mini_smoke=bool(data.get("mini_smoke")),
                release_audit=not bool(data.get("skip_release_audit", False)),
                allow_execute=bool(data.get("execute", False)),
            )
            return {"status": res.status, "outputs": res.outputs}
        if route == "/pipeline-full":
            run = run_pipeline_full(
                sample_manifest=data["sample_manifest"],
                tools_manifest=data.get("tools_manifest"),
                reference_manifest=data.get("reference_manifest"),
                project_root=data.get("project_root") or self.project_root,
                outdir=outdir,
                profile=data.get("profile", "local"),
                dry_run=not bool(data.get("execute", False)),
                allow_partial=not bool(data.get("strict", False)),
            )
            return {"status": run.status, "run_id": run.run_id, "outputs": {"run_manifest": str(Path(run.output_dir) / "run_manifest.json"), "pipeline_plan": str(Path(run.output_dir) / "pipeline_plan.md")}}
        if route == "/input-qc":
            argv = ["--outdir", str(outdir)]
            if data.get("manifest"):
                argv += ["--manifest", data["manifest"]]
            if data.get("result_dir"):
                argv += ["--result-dir", data["result_dir"]]
            err, rc = _run_skill_func(input_qc_main, argv)
            return {"status": "PASS" if rc == 0 else "FAIL", "error": err, "outputs": {"input_status": str(outdir / "input_status.json")}}
        if route == "/ranking-compare":
            argv = ["--recommendation", data["recommendation"], "--netmhcpan42", data["netmhcpan42"], "--outdir", str(outdir)]
            err, rc = _run_skill_func(ranking_compare_main, argv)
            return {"status": "PASS" if rc == 0 else "FAIL", "error": err, "outputs": {"report": str(outdir / "ranking_compare_report.md")}}
        if route == "/appm-review":
            argv = ["--outdir", str(outdir)]
            for key in ["evidence_report", "appm_summary", "appm_gene_status", "appm_submodule_scores", "peptide_escape_flags", "ranked_peptides", "hla_loh"]:
                if data.get(key):
                    argv += ["--" + key.replace("_", "-"), data[key]]
            err, rc = _run_skill_func(appm_review_main, argv)
            return {"status": "PASS" if rc == 0 else "FAIL", "error": err, "outputs": {"report": str(outdir / "appm_escape_review.md")}}
        if route == "/ccf-review":
            argv = ["--outdir", str(outdir)]
            for key in ["ccf", "purity", "ranked_peptides"]:
                if data.get(key):
                    argv += ["--" + key.replace("_", "-"), data[key]]
            err, rc = _run_skill_func(ccf_review_main, argv)
            return {"status": "PASS" if rc == 0 else "FAIL", "error": err, "outputs": {"report": str(outdir / "ccf_clonality_review.md")}}
        if route == "/patient-report":
            argv = ["--outdir", str(outdir)]
            for key in ["evidence_report", "recommendation", "netmhcpan42", "ranking_compare_report", "appm_review", "ccf_review"]:
                if data.get(key):
                    argv += ["--" + key.replace("_", "-"), data[key]]
            err, rc = _run_skill_func(patient_report_main, argv)
            return {"status": "PASS" if rc == 0 else "FAIL", "error": err, "outputs": {"report_md": str(outdir / "patient_report.md")}}
        return {"status": "FAIL", "error": f"unhandled route {route}"}


def _allowed_roots(project_root: Path, outdir: Path, extra_roots: list[str] | None = None) -> list[Path]:
    roots = [project_root.resolve(), outdir.resolve()]
    env_roots = os.environ.get("NEOAG_GATEWAY_ALLOWED_ROOTS", "")
    for raw in list(extra_roots or []) + [x for x in env_roots.split(os.pathsep) if x]:
        roots.append(Path(raw).expanduser().resolve())
    dedup: list[Path] = []
    for root in roots:
        if root not in dedup:
            dedup.append(root)
    return dedup


def run_gateway(*, host: str, port: int, project_root: str | Path, outdir: str | Path, allowed_root: list[str] | None = None) -> None:
    od = ensure_dir(outdir)
    server = ThreadingHTTPServer((host, port), GatewayHandler)
    server.project_root = Path(project_root).resolve()  # type: ignore[attr-defined]
    server.job_store = JobStore(od / "jobs")  # type: ignore[attr-defined]
    server.audit_logger = AuditLogger(od / "audit_log.jsonl")  # type: ignore[attr-defined]
    server.allowed_roots = _allowed_roots(server.project_root, od, allowed_root)  # type: ignore[attr-defined]
    print(f"NeoAg Gateway listening on http://{host}:{port}")
    print(f"  outdir: {od}")
    print("  allowed_roots:")
    for root in server.allowed_roots:  # type: ignore[attr-defined]
        print(f"    - {root}")
    server.serve_forever()


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="NeoAg Gateway: controlled execution gateway")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8000)
    ap.add_argument("--project-root", default=".")
    ap.add_argument("--outdir", default="work/neoag_gateway")
    ap.add_argument("--allowed-root", action="append", default=[], help="Additional root allowed for write/output paths; can be repeated")
    args = ap.parse_args(argv)
    run_gateway(host=args.host, port=args.port, project_root=args.project_root, outdir=args.outdir, allowed_root=args.allowed_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
