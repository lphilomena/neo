"""LLM-assisted Coordinator for Project B.

This package implements the P1.2 coordinator layer: LLM-assisted intent
classification, skill planning, guardrails, safe skill execution, result
collection, and case-state/audit logging. It is intentionally optional:
all code works in rule-based/mock mode without installing LangGraph, LiteLLM,
or any cloud model SDK.
"""

__all__ = ["__version__"]
__version__ = "0.1.0"
