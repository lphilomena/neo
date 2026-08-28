from __future__ import annotations

from collections import defaultdict, deque
from itertools import combinations
from typing import Any


def _identifier(row: dict[str, Any], index: int) -> str:
    return str(row.get("_pareto_id") or row.get("peptide_id") or f"row:{index}")


def _vector(row: dict[str, Any], dimensions: list[str]) -> tuple[float, ...]:
    values: list[float] = []
    for dimension in dimensions:
        try:
            values.append(float(row.get(dimension, 0)))
        except (TypeError, ValueError):
            values.append(0.0)
    return tuple(values)


def _dominates(left: tuple[float, ...], right: tuple[float, ...]) -> bool:
    return all(a >= b for a, b in zip(left, right)) and any(a > b for a, b in zip(left, right))


def nondominated_fronts(
    rows: list[dict[str, Any]],
    dimensions: list[str],
) -> dict[str, int]:
    """Return deterministic Pareto fronts after collapsing duplicate vectors.

    Larger values are preferred for every dimension. Pairwise comparisons are
    performed only for unique discrete vectors, so repeated peptide-HLA rows do
    not make complexity quadratic in the full candidate count.
    """

    if not dimensions:
        raise ValueError("Pareto ranking requires at least one dimension")
    vector_members: dict[tuple[float, ...], list[str]] = defaultdict(list)
    for index, row in enumerate(rows):
        vector_members[_vector(row, dimensions)].append(_identifier(row, index))
    vectors = sorted(vector_members, reverse=True)
    outgoing: dict[tuple[float, ...], list[tuple[float, ...]]] = {vector: [] for vector in vectors}
    indegree = {vector: 0 for vector in vectors}
    for left, right in combinations(vectors, 2):
        if _dominates(left, right):
            outgoing[left].append(right)
            indegree[right] += 1
        elif _dominates(right, left):
            outgoing[right].append(left)
            indegree[left] += 1

    queue = deque(sorted((vector for vector in vectors if indegree[vector] == 0), reverse=True))
    vector_front: dict[tuple[float, ...], int] = {}
    current_front = 1
    while queue:
        next_queue: list[tuple[float, ...]] = []
        while queue:
            vector = queue.popleft()
            vector_front[vector] = current_front
            for child in outgoing[vector]:
                indegree[child] -= 1
                if indegree[child] == 0:
                    next_queue.append(child)
        queue = deque(sorted(next_queue, reverse=True))
        current_front += 1

    result: dict[str, int] = {}
    for vector, members in vector_members.items():
        for member in members:
            result[member] = vector_front[vector]
    return result
