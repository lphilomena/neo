# Open-Neo macro Skills v1

- Added public category M and three public macro Skills.
- Added `open-neo` CLI with `install-check`, `run`, and `review` subcommands.
- Added deterministic input inspection and route planning.
- Added plan/dry-run/execute/resume/ranking-only modes.
- Added event-level review, first-batch experiment heuristic, and patient/technical reports.
- Preserved all existing A/B/C/D Skills and the weighted ranking baseline.
- Installation/repair now delegates to the portable new-machine installer after explicit approval; downloads remain separately opt-in.
- Missing HLA, required references, or supplied files now block execution instead of becoming optional evidence warnings.
- Resume reuses completed outputs by default; force reruns require an explicit flag.
- Review accepts consensus-only results, keeps R3 in an evidence-completion queue, uses clinical context, and prevents writing into the source result directory.
- Preserved the existing `neoag-workflow-select` CLI while adding the `open-neo` CLI.
