"""Hardware coupling graphs used by the routing pass."""

from __future__ import annotations

import math


def all_to_all(n: int) -> list[tuple[int, int]]:
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def line(n: int) -> list[tuple[int, int]]:
    return [(i, i + 1) for i in range(n - 1)]


def ring(n: int) -> list[tuple[int, int]]:
    return line(n) + ([(n - 1, 0)] if n > 2 else [])


def square_grid(n: int) -> list[tuple[int, int]]:
    """Row-major grid on ``ceil(sqrt(n)) x ceil(n / ceil(sqrt(n)))``."""
    w = int(math.ceil(math.sqrt(n)))
    h = int(math.ceil(n / w))
    edges = []
    for r in range(h):
        for c in range(w):
            q = r * w + c
            if q >= n:
                continue
            if c + 1 < w and r * w + c + 1 < n:
                edges.append((q, r * w + c + 1))
            if r + 1 < h and (r + 1) * w + c < n:
                edges.append((q, (r + 1) * w + c))
    return edges


def heavy_hex(n: int) -> list[tuple[int, int]]:
    """A heavy-hex-like subgraph: a line with degree-3 'flag' qubits attached.

    Not an exact device layout; a degree-limited connected subgraph with the
    same routing character (max degree 3).  Declared as such rather than
    labelled with a vendor name.
    """
    edges = line(n)
    for q in range(2, n - 2, 4):
        target = q + 3
        if target < n and len([e for e in edges if q in e]) < 3:
            edges.append((q, target))
    return sorted({tuple(sorted(e)) for e in edges})


def modular_two_cluster(n: int) -> list[tuple[int, int]]:
    """Two all-to-all clusters joined by a single inter-module link."""
    half = n // 2
    edges = [(i, j) for i in range(half) for j in range(i + 1, half)]
    edges += [(i, j) for i in range(half, n) for j in range(i + 1, n)]
    if half < n:
        edges.append((half - 1, half))
    return sorted({tuple(sorted(e)) for e in edges})


def neutral_atom_grid(n: int, radius: float = 1.5) -> list[tuple[int, int]]:
    """Programmable geometry: unit-square lattice with a blockade radius."""
    w = int(math.ceil(math.sqrt(n)))
    pos = {q: (q % w, q // w) for q in range(n)}
    edges = []
    for i in range(n):
        for j in range(i + 1, n):
            dx = pos[i][0] - pos[j][0]
            dy = pos[i][1] - pos[j][1]
            if math.hypot(dx, dy) <= radius + 1e-9:
                edges.append((i, j))
    return edges


TOPOLOGIES = {
    "all_to_all": all_to_all,
    "line": line,
    "ring": ring,
    "square_grid": square_grid,
    "heavy_hex": heavy_hex,
    "modular_two_cluster": modular_two_cluster,
    "neutral_atom_grid": neutral_atom_grid,
}


def build(name: str, n: int) -> list[tuple[int, int]]:
    if name not in TOPOLOGIES:
        raise KeyError(f"unknown topology {name!r}; have {sorted(TOPOLOGIES)}")
    return TOPOLOGIES[name](n)
