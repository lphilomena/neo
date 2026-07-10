from __future__ import annotations

import json
import os
from dataclasses import dataclass
from typing import Any

from .json_utils import extract_json_object
from .schema_validation import validate_intent_payload
from .schemas import IntentResult, validate_intent
from neoag_v03.agents.skill_router import classify_intent


@dataclass
class LLMResponse:
    text: str
    provider: str
    raw: Any | None = None


class BaseModelProvider:
    name = "base"

    def complete(self, messages: list[dict[str, str]], temperature: float = 0.0, response_format: str | None = None) -> LLMResponse:
        raise NotImplementedError


class RuleBasedProvider(BaseModelProvider):
    """Offline deterministic provider.

    This is not a real LLM. It lets tests and secured deployments run the same
    Coordinator path without network or model dependencies. It is also used as a
    fallback when LLM JSON is invalid or confidence is low.
    """

    name = "rule"

    def complete(self, messages: list[dict[str, str]], temperature: float = 0.0, response_format: str | None = None) -> LLMResponse:
        user = "\n".join(m.get("content", "") for m in messages if m.get("role") == "user")
        intent = classify_intent(user)
        text = json.dumps(
            {
                "intent": intent,
                "confidence": 0.76,
                "user_goal": user[:200],
                "requires_execution": intent not in {"general_explanation", "unknown"},
                "requires_human_approval": any(k in user.lower() for k in ["hpc", "slurm", "sbatch", "安装", "删除", "覆盖"]),
                "missing_information_questions": [],
            },
            ensure_ascii=False,
        )
        return LLMResponse(text=text, provider=self.name)


class LiteLLMProvider(BaseModelProvider):
    """LiteLLM-backed provider.

    Set model names through --model or NEOAG_LLM_MODEL. For local vLLM, configure
    LiteLLM with an OpenAI-compatible endpoint, e.g. --api-base http://localhost:8000/v1.
    """

    name = "litellm"

    def __init__(self, model: str, api_base: str | None = None, api_key_env: str = "OPENAI_API_KEY") -> None:
        self.model = model
        self.api_base = api_base
        self.api_key_env = api_key_env

    def complete(self, messages: list[dict[str, str]], temperature: float = 0.0, response_format: str | None = None) -> LLMResponse:
        try:
            from litellm import completion  # type: ignore
        except Exception as e:  # pragma: no cover - optional dependency
            raise RuntimeError("LiteLLM is not installed. Install with project extra [agent-llm] or use --llm-provider rule.") from e
        kwargs: dict[str, Any] = {
            "model": self.model,
            "messages": messages,
            "temperature": temperature,
        }
        if self.api_base:
            kwargs["api_base"] = self.api_base
        if os.environ.get(self.api_key_env):
            kwargs["api_key"] = os.environ.get(self.api_key_env)
        if response_format == "json_object":
            kwargs["response_format"] = {"type": "json_object"}
        resp = completion(**kwargs)
        try:
            text = resp.choices[0].message.content or ""
        except Exception:
            text = str(resp)
        return LLMResponse(text=text, provider=self.name, raw=resp)


def make_provider(kind: str = "rule", model: str | None = None, api_base: str | None = None, api_key_env: str = "OPENAI_API_KEY") -> BaseModelProvider:
    kind = (kind or "rule").lower()
    if kind in {"rule", "mock", "offline"}:
        return RuleBasedProvider()
    if kind in {"litellm", "openai-compatible", "vllm", "qwen", "deepseek"}:
        return LiteLLMProvider(model=model or os.environ.get("NEOAG_LLM_MODEL", "qwen/qwen3-32b"), api_base=api_base or os.environ.get("NEOAG_LLM_API_BASE"), api_key_env=api_key_env)
    raise ValueError(f"Unknown provider: {kind}")


def parse_intent_response(text: str, source: str = "llm") -> IntentResult | None:
    obj = extract_json_object(text)
    if not obj:
        return None
    validate_intent_payload(obj)
    intent = validate_intent(str(obj.get("intent", "unknown")))
    try:
        confidence = float(obj.get("confidence", 0.0))
    except Exception:
        confidence = 0.0
    questions = obj.get("missing_information_questions") or []
    if not isinstance(questions, list):
        questions = [str(questions)]
    return IntentResult(
        intent=intent,
        confidence=confidence,
        user_goal=str(obj.get("user_goal", "")),
        requires_execution=bool(obj.get("requires_execution", intent not in {"unknown", "general_explanation"})),
        requires_human_approval=bool(obj.get("requires_human_approval", False)),
        missing_information_questions=[str(q) for q in questions],
        source=source,
    )
