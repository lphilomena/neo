"""Controlled execution layer for Project B.

This package implements Phase 0-2 enhancement code:
- Phase 0: release boundary audit and NeoAg Doctor.
- Phase 1: stdlib Gateway with safe JSON endpoints.
- Phase 2: manifest-driven pipeline-full runner.

All modules are intentionally dependency-light and safe-by-default. Heavy
bioinformatics tools are checked via manifests and smoke hooks but are not
bundled or executed unless explicitly requested.
"""

from .version import ENHANCEMENT_VERSION

__all__ = ["ENHANCEMENT_VERSION"]
