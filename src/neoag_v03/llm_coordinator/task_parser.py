from __future__ import annotations

import json
from typing import Any

from neoag_v03.agents.skill_router import classify_intent
from .json_utils import extract_json_object
from .model_provider import BaseModelProvider
from .prompts import COORDINATOR_SYSTEM_PROMPT, TASK_JSON_PROMPT_TEMPLATE
from .schemas import ParsedTask, validate_intent
from .tool_registry import registry_for_prompt


def _risk_for_message(message: str, intent: str) -> str:
    m = message.lower()
    if any(x in m for x in ["删除", "覆盖", "push", "上传", "安装", "rm ", "force"]):
        return "high"
    if intent in {"sliding_run", "workflow_run_request", "setup_tool"} or any(x in m for x in ["运行", "跑", "run", "执行"]):
        return "medium"
    return "low"


def _parameters_for_message(message: str) -> dict[str, Any]:
    params: dict[str, Any] = {}
    m = message.strip()
    match = __import__("re").search(r"(?:sample[_ -]?id|样本(?:编号|ID)?)[：:= ]+([A-Za-z0-9_.-]+)", m, flags=__import__("re").I)
    if match:
        params["sample_id"] = match.group(1)
    lower = m.lower()
    if "salmon" in lower:
        params["method"] = "salmon"
    elif "rsem" in lower:
        params["method"] = "rsem"
    return params


def _workflow_for_message(message: str) -> str | None:
    m = message.lower()
    for key in ["rna_tpm", "fusion_rna", "hla_typing", "purity_cnv", "spechla", "hla-la", "optitype", "facets", "purple", "sequenza", "ascat", "easyfuse", "lohhla", "netmhcpan", "sliding"]:
        if key in m:
            return key
    if any(x in m for x in ["tpm", "表达量", "表达矩阵", "rsem", "salmon", "kallisto", "fastq生成tpm", "fastq 生成 tpm"]):
        return "rna_tpm"
    if any(x in m for x in ["easyfuse", "star-fusion", "star_fusion", "arriba", "fusioncatcher", "fusion", "融合"]):
        return "fusion_rna"
    if any(x in m for x in ["hla分型", "hla typing", "分型结果"]):
        return "hla_typing"
    if any(x in m for x in ["纯度", "cnv", "copy number", "拷贝数"]):
        return "purity_cnv"
    return None


def rule_parse_task(message: str, available_files: list[dict[str, Any]] | None = None) -> ParsedTask:
    intent = validate_intent(classify_intent(message))
    workflow = _workflow_for_message(message)
    wants_run = any(x in message.lower() for x in ["运行", "跑", "run", "执行", "生成", "量化", "quant"])
    # Keep specific workflow intents when available so the planner can select the right skill.
    if workflow and wants_run and intent == "input_check":
        intent = "workflow_run_request"
    inputs = [str(f.get("path")) for f in (available_files or []) if f.get("path")]
    action = "explain"
    if intent in {"sliding_run", "workflow_run_request", "setup_tool"} or wants_run:
        action = "run"
    elif intent == "check_status":
        action = "status"
    elif intent in {"inspect_results", "ranking_compare"}:
        action = "inspect"
    elif intent in {"git_release", "update_docs", "data_transfer", "debug_error"}:
        action = "assist"
    return ParsedTask(
        intent=intent,
        task_type=intent,
        action=action,
        workflow=workflow,
        target=workflow,
        inputs=inputs,
        parameters=_parameters_for_message(message),
        risk=_risk_for_message(message, intent),
        missing_fields=[],
        needs_status_tracking=action in {"run", "status", "inspect", "assist"},
        user_visible_summary=message[:240],
        source="rule",
    )


def _coerce_task(obj: dict[str, Any], fallback: ParsedTask, source: str) -> ParsedTask:
    intent = validate_intent(str(obj.get("intent") or fallback.intent))
    risk = str(obj.get("risk") or fallback.risk).lower()
    if risk not in {"low", "medium", "high"}:
        risk = fallback.risk
    raw_params = obj.get("parameters") if isinstance(obj.get("parameters"), dict) else {}
    raw_inputs = obj.get("inputs") if isinstance(obj.get("inputs"), list) else fallback.inputs
    return ParsedTask(
        intent=intent,
        task_type=str(obj.get("task_type") or intent),
        action=str(obj.get("action") or fallback.action),
        workflow=str(obj.get("workflow") or fallback.workflow or "") or None,
        target=str(obj.get("target") or fallback.target or "") or None,
        inputs=[str(x) for x in raw_inputs],
        parameters={str(k): v for k, v in raw_params.items()},
        risk=risk,
        missing_fields=[str(x) for x in (obj.get("missing_fields") or [])],
        needs_status_tracking=bool(obj.get("needs_status_tracking", fallback.needs_status_tracking)),
        user_visible_summary=str(obj.get("user_visible_summary") or fallback.user_visible_summary),
        source=source,
    )


def parse_task_with_llm(message: str, provider: BaseModelProvider, available_files: list[dict[str, Any]], tool_registry: dict[str, Any]) -> ParsedTask:
    fallback = rule_parse_task(message, available_files)
    if getattr(provider, "name", "") == "rule":
        return fallback
    prompt = TASK_JSON_PROMPT_TEMPLATE.format(
        message=message,
        available_files=json.dumps(available_files, ensure_ascii=False, indent=2)[:5000],
        tool_registry=registry_for_prompt(tool_registry),
    )
    try:
        resp = provider.complete([
            {"role": "system", "content": COORDINATOR_SYSTEM_PROMPT},
            {"role": "user", "content": prompt},
        ], temperature=0.0, response_format="json_object")
        obj = extract_json_object(resp.text)
        if obj:
            return _coerce_task(obj, fallback, resp.provider)
    except Exception:
        pass
    return fallback
