# Changelog: LLM-assisted Coordinator P1

Added an optional LLM-assisted Coordinator layer for Project B:

- `neoag-llm-agent` CLI.
- `src/neoag_v03/llm_coordinator/` package.
- Rule-based/offline provider plus optional LiteLLM provider.
- Context builder, input-state builder, intent classifier, Skill planner, guardrails, Skill executor, result collector, summarizer, case-state and audit logging.
- Optional LangGraph adapter in `llm_coordinator/langgraph_app.py`.
- Example LiteLLM/vLLM/Qwen configuration files under `configs/llm/`.
- Example startup scripts for local vLLM and LiteLLM Gateway.
- Documentation for model API and agent framework selection.

Default mode has no external LLM dependency. Install optional dependencies with:

```bash
python -m pip install -e '.[agent-llm]'
```
