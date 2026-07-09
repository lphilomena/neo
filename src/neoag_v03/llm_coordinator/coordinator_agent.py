from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any

from neoag_v03.agents.skill_router import find_named_files
from .case_state import append_audit, new_case_id, write_case_state
from .context_builder import build_context
from .guardrails import computational_boundary_note
from .input_state import build_input_state
from .model_provider import make_provider
from .planner import classify_with_llm, plan_with_llm
from .task_parser import parse_task_with_llm
from .tool_registry import load_tool_registry
from .state_tracker import append_task_event
from .result_collector import flatten_output_links
from .schemas import CaseState, ExecutionMode
from .skill_executor import execute_plan, load_registry
from .summarizer import summarize_with_llm


def _mode(value: str) -> ExecutionMode:
    if value not in {"plan", "execute-safe", "execute-with-approval"}:
        raise argparse.ArgumentTypeError("mode must be plan, execute-safe, or execute-with-approval")
    return value  # type: ignore[return-value]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="LLM-assisted Coordinator Agent for Project B Skills")
    ap.add_argument("--message", required=True, help="User request in natural language")
    ap.add_argument("--file", action="append", default=[], help="Input/result file path; may be repeated")
    ap.add_argument("--result-dir", help="Existing Project B result directory")
    ap.add_argument("--project-root", default=".")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--mode", type=_mode, default="plan", help="plan | execute-safe | execute-with-approval")
    ap.add_argument("--allow-high-risk", action="store_true", help="Allow approval-gated steps in execute-with-approval mode")
    ap.add_argument("--llm-provider", default="rule", help="rule | litellm | vllm | qwen | deepseek")
    ap.add_argument("--model", help="Model name for LiteLLM/vLLM/Qwen/DeepSeek")
    ap.add_argument("--api-base", help="OpenAI-compatible API base, e.g. local vLLM http://localhost:8000/v1")
    ap.add_argument("--api-key-env", default="OPENAI_API_KEY")
    ap.add_argument("--temperature", type=float, default=0.0)
    ap.add_argument("--case-id")
    ap.add_argument("--sample-id", help="Optional sample ID used by result-review skills to filter outputs")
    ap.add_argument("--execute-input-qc", action="store_true", help="Run input-qc before planning to build deterministic input state")
    ap.add_argument("--no-llm-summary", action="store_true", help="Use deterministic final summary")
    ap.add_argument("--min-intent-confidence", type=float, default=0.80)
    args = ap.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    case_id = args.case_id or new_case_id("neoag_llm")
    audit_path = outdir / "audit_log.jsonl"
    append_audit(audit_path, {"ts": time.time(), "event": "start", "case_id": case_id, "message": args.message, "mode": args.mode})

    context = build_context(args.message, files=args.file, result_dir=args.result_dir, project_root=args.project_root)
    (outdir / "context.json").write_text(json.dumps(context, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    input_state = build_input_state(args.file, args.result_dir, outdir=outdir, execute_input_qc=args.execute_input_qc)
    provider = make_provider(args.llm_provider, model=args.model, api_base=args.api_base, api_key_env=args.api_key_env)
    registry = load_registry(args.project_root)
    tool_registry = load_tool_registry(args.project_root, registry)
    parsed_task = parse_task_with_llm(args.message, provider, context.get("available_files", []), tool_registry)
    if args.sample_id:
        parsed_task.parameters["sample_id"] = args.sample_id
    (outdir / "task_spec.json").write_text(json.dumps(parsed_task.to_dict(), ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    append_task_event(args.project_root, outdir, {"event": "task_parsed", "case_id": case_id, "task": parsed_task.to_dict(), "outdir": str(outdir)})
    file_kinds = [str(f.get("kind")) for f in context.get("available_files", [])]
    intent = classify_with_llm(args.message, provider, file_kinds=file_kinds, min_confidence=args.min_intent_confidence)
    if parsed_task.intent != "unknown":
        intent.intent = parsed_task.intent
        intent.user_goal = parsed_task.user_visible_summary or intent.user_goal
        intent.requires_execution = parsed_task.action in {"run", "assist"} and parsed_task.risk != "low"
        intent.requires_human_approval = parsed_task.risk == "high"
        intent.source = parsed_task.source
    append_audit(audit_path, {"ts": time.time(), "event": "intent", "intent": intent.to_dict(), "task": parsed_task.to_dict()})

    plan = plan_with_llm(args.message, intent, input_state, tool_registry, provider)
    append_audit(audit_path, {"ts": time.time(), "event": "plan", "plan": plan.to_dict()})

    # Always write plan before any execution.
    (outdir / "coordinator_plan.json").write_text(json.dumps(plan.to_dict(), ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    plan_lines = ["# LLM-assisted Coordinator Plan", "", f"Case ID: `{case_id}`", f"Intent: **{intent.intent}**", f"Intent source: **{intent.source}**", f"Mode: **{args.mode}**", "", "## Steps"]
    for step in plan.steps:
        plan_lines.append(f"- `{step.step_id}` **{step.skill}** ({step.mode}) — {step.reason or 'registered skill step'}")
    if plan.approval_required:
        plan_lines += ["", "## Approval required", plan.approval_reason or "This plan contains approval-gated operations."]
    if plan.missing_inputs:
        plan_lines += ["", "## Missing inputs"]
        for m in plan.missing_inputs:
            plan_lines.append(f"- {m}")
    plan_lines += ["", "## Boundary", computational_boundary_note()]
    (outdir / "coordinator_plan.md").write_text("\n".join(plan_lines) + "\n", encoding="utf-8")

    files = find_named_files(args.result_dir, args.file)
    if args.mode == "plan":
        results = execute_plan(plan, files, args.result_dir, outdir, args.project_root, mode="plan", allow_high_risk=False, task_parameters=parsed_task.parameters)
    else:
        allow = args.mode == "execute-with-approval" and args.allow_high_risk
        results = execute_plan(plan, files, args.result_dir, outdir, args.project_root, mode=args.mode, allow_high_risk=allow, task_parameters=parsed_task.parameters)
    append_audit(audit_path, {"ts": time.time(), "event": "skill_results", "results": [r.to_dict() for r in results]})
    append_task_event(args.project_root, outdir, {"event": "skill_results", "case_id": case_id, "status": "PASS" if all(r.status in {"PASS", "PLANNED", "SKIPPED", "APPROVAL_REQUIRED"} for r in results) else "FAIL", "outdir": str(outdir), "results": [r.to_dict() for r in results]})

    outputs = flatten_output_links(results)
    final_text, validation_warnings = summarize_with_llm(args.message, intent, plan, results, provider, use_llm=not args.no_llm_summary and args.llm_provider != "rule")
    if validation_warnings:
        append_audit(audit_path, {"ts": time.time(), "event": "output_guardrail_warnings", "warnings": validation_warnings})
    (outdir / "final_response.md").write_text(final_text, encoding="utf-8")

    state = CaseState(
        case_id=case_id,
        message=args.message,
        mode=args.mode,
        intent=intent,
        input_state=input_state,
        plan=plan,
        skill_results=results,
        warnings=validation_warnings,
        outputs=outputs,
        boundary_note=computational_boundary_note(),
        audit_log=[{"audit_log": str(audit_path)}],
    )
    write_case_state(outdir / "case_state.json", state)
    append_task_event(args.project_root, outdir, {"event": "finished", "case_id": case_id, "status": "PASS" if all(r.status in {"PASS", "PLANNED", "SKIPPED", "APPROVAL_REQUIRED"} for r in results) else "FAIL", "outdir": str(outdir), "outputs": outputs})
    print(json.dumps({"case_id": case_id, "status": "PASS" if all(r.status in {"PASS", "PLANNED", "SKIPPED", "APPROVAL_REQUIRED"} for r in results) else "FAIL", "outdir": str(outdir), "plan": str(outdir / "coordinator_plan.md"), "case_state": str(outdir / "case_state.json"), "final_response": str(outdir / "final_response.md"), "outputs": outputs}, ensure_ascii=False, indent=2))
    return 0 if all(r.status in {"PASS", "PLANNED", "SKIPPED", "APPROVAL_REQUIRED"} for r in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
