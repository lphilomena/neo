from __future__ import annotations

from typing import Any, TypedDict

class NeoagGraphState(TypedDict, total=False):
    message: str
    files: list[str]
    result_dir: str | None
    project_root: str
    outdir: str
    context: dict[str, Any]
    input_state: dict[str, Any]
    intent: dict[str, Any]
    plan: dict[str, Any]
    skill_results: list[dict[str, Any]]
    final_response: str


def build_langgraph_app():
    """Build a LangGraph app when langgraph is installed.

    This optional adapter is intentionally light. The production CLI uses the
    deterministic coordinator modules directly; this function gives teams a
    starting point for deploying the same nodes in a LangGraph runtime.
    """
    try:
        from langgraph.graph import StateGraph, END  # type: ignore
    except Exception as e:  # pragma: no cover - optional dependency
        raise RuntimeError("LangGraph is not installed. Install optional extra [agent-llm] to use this adapter.") from e

    def build_context_node(state: NeoagGraphState) -> NeoagGraphState:
        from .context_builder import build_context
        state["context"] = build_context(state.get("message", ""), state.get("files", []), state.get("result_dir"), state.get("project_root", "."))
        return state

    def input_state_node(state: NeoagGraphState) -> NeoagGraphState:
        from .input_state import build_input_state
        istate = build_input_state(state.get("files", []), state.get("result_dir"), outdir=state.get("outdir"), execute_input_qc=False)
        state["input_state"] = istate.to_dict()
        return state

    graph = StateGraph(NeoagGraphState)
    graph.add_node("build_context", build_context_node)
    graph.add_node("input_state", input_state_node)
    graph.set_entry_point("build_context")
    graph.add_edge("build_context", "input_state")
    graph.add_edge("input_state", END)
    return graph.compile()
