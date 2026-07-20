# Project B LLM-assisted Coordinator P1

This document describes the second-step agent layer for Project B: a **LLM-assisted Coordinator** on top of the existing P0 Skills Pack. It does not replace Project B CLI, Nextflow, or external bioinformatics tools. It routes user requests to registered Skills, applies guardrails, collects results, and summarizes next steps.

## Why this layer exists

Project B already has most of the neoantigen analysis backbone: HLA typing, SNV/InDel, fusion, peptide generation, binding/presentation prediction, safety, APPM/immune escape, CCF/clonality, recommendation ranking, and reports. The missing layer is interactive task planning: natural-language intent recognition, input-state checking, Skill routing, clarification, and safe report generation.

## Architecture

```text
User message + uploaded files
        ↓
LLM-assisted Coordinator
        ↓
Context Builder → Input State Builder → Intent Router → Skill Planner
        ↓
Guardrails + Human Approval Gate
        ↓
Skill Executor
        ↓
Result Collector → LLM/Deterministic Summarizer
        ↓
case_state.json + coordinator_plan.md + final_response.md
```

## CLI

Plan-only mode, no Skill execution:

```bash
neoag-llm-agent \
  --message "比较 recommendation 和 NetMHCpan42 排序差异" \
  --file ranked_peptides.recommendation.tsv \
  --file ranked_peptides.netmhcpan42.tsv \
  --outdir work/llm_plan \
  --mode plan
```

Execute safe Skills:

```bash
neoag-llm-agent \
  --message "比较 recommendation 和 NetMHCpan42 排序差异" \
  --file ranked_peptides.recommendation.tsv \
  --file ranked_peptides.netmhcpan42.tsv \
  --outdir work/llm_execute \
  --mode execute-safe
```

Use LiteLLM / local vLLM / Qwen-compatible endpoint:

```bash
neoag-llm-agent \
  --message "请根据最新结果更新患者沟通版 Word" \
  --file evidence_report.v04x_latest.html \
  --file ranked_peptides.recommendation.tsv \
  --file ranked_peptides.netmhcpan42.tsv \
  --outdir work/llm_patient_report \
  --mode execute-safe \
  --llm-provider litellm \
  --model openai/qwen3-32b \
  --api-base http://localhost:8000/v1 \
  --api-key-env LOCAL_VLLM_API_KEY
```

## Modes

- `plan`: write `coordinator_plan.*` and `case_state.json`; do not execute Skills.
- `execute-safe`: execute read-only/low-risk Skills such as input QC, ranking compare, APPM review, CCF review, and draft patient report generation.
- `execute-with-approval`: can execute approval-gated steps only with `--allow-high-risk`.

High-risk operations include HPC submission, tool installation, deletion, overwriting result directories, downloading large references, and generating formal patient-facing reports without review.

## Outputs

- `context.json`: discovered files and brief context.
- `coordinator_plan.json`: structured Skill plan.
- `coordinator_plan.md`: human-readable plan.
- `case_state.json`: persistent case state.
- `audit_log.jsonl`: execution/audit events.
- `final_response.md`: final user-facing summary.
- Per-Skill subdirectories, for example `neoag-ranking-compare/`.

## Model API recommendation

- Production/private patient data: local Qwen-family open-weight model served by vLLM OpenAI-compatible API.
- Unified gateway: LiteLLM.
- De-identified development fallback: DeepSeek API or another OpenAI-compatible provider.
- Small classification tasks can use a smaller local model; reporting/planning can use a larger local model.

## Framework recommendation

The code is dependency-light by default. Optional LangGraph integration is provided in `neoag.llm_coordinator.langgraph_app`. Use LangGraph when you want durable stateful orchestration, human-in-the-loop nodes, and deployment-grade graph runtime.

## Clinical/research boundary

The Coordinator must not state that a computational candidate is a confirmed neoantigen, must not promise clinical benefit, and must not diagnose resistance yes/no. Missing evidence must be reported as missing/unassessed, not negative.
