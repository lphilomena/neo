# Model API and Agent Framework Selection for Project B

## Recommended stack

```text
Agent orchestration: LangGraph, optional adapter
LLM gateway: LiteLLM
Self-hosted inference: vLLM OpenAI-compatible server
Primary local model: Qwen-family open-weight model
Development fallback: DeepSeek API for de-identified tasks only
Skills: .agents/skills/neoag-* P0 Skills Pack
Execution: Project B CLI / Nextflow / HPC wrappers
State: case_state.json, audit_log.jsonl
```

## Rationale

Project B is not a simple chatbot. It requires controlled multi-step workflows, file-state inspection, Skill invocation, human approval gates, and strict research-use boundaries. A local model behind vLLM protects patient data, LiteLLM provides a provider-neutral API, and LangGraph can be adopted for durable graph orchestration without changing Project B Skills.

## Model profile policy

- Real patient/sample data should stay on local/approved infrastructure.
- Cloud APIs can be used only with de-identified summaries or development fixtures.
- Tool-calling and planning should run at low temperature with structured JSON.
- Patient reports must pass boundary checks before delivery.

## Framework policy

- Do not let an LLM directly execute arbitrary shell commands.
- Use Skills for SOPs and Project B CLI/Nextflow for computation.
- Use an execution gateway or Skill executor to enforce dry-run, approval, and provenance.
- Add specialist agents only after Coordinator + Skills are stable.
