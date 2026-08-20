"""Framework-neutral circuit IR, lowering to {1q, CX}, and hardware routing.

Conventions (declared once, used everywhere):

* Qubit ordering is little-endian: basis index ``b`` assigns qubit ``q`` the bit
  ``(b >> q) & 1``.  Qubit 0 is the least significant bit.
* A :class:`Circuit` may declare trailing *clean* ancillas.  Every construction
  in this package returns its ancillas to |0> so that the data register is
  recovered exactly by projecting onto the all-zero ancilla slice.
* Resource counts are always taken from an *emitted* circuit after
  :func:`lower_to_cx` (and after :func:`route` when a topology is declared).
  No count in this package is produced by a closed-form formula.
"""

from __future__ import annotations

import math
from collections import deque
from dataclasses import dataclass, field
from typing import Iterable, Sequence

ONE_QUBIT_GATES = frozenset({"h", "x", "y", "z", "s", "sdg", "t", "tdg", "ry", "rz"})
TWO_QUBIT_GATES = frozenset({"cx", "cz", "swap"})
MULTI_GATES = frozenset({"mcx", "mcry"})
NATIVE_GATES = ONE_QUBIT_GATES | {"cx"}


@dataclass(frozen=True)
class Gate:
    """A single instruction.

    For ``mcx``/``mcry`` the final entry of ``qubits`` is the target and the
    preceding entries are the controls.
    """

    name: str
    qubits: tuple[int, ...]
    param: float = 0.0

    @property
    def target(self) -> int:
        return self.qubits[-1]

    @property
    def controls(self) -> tuple[int, ...]:
        return self.qubits[:-1]


@dataclass
class Circuit:
    """An ordered instruction list over ``n_qubits`` wires."""

    n_qubits: int
    gates: list[Gate] = field(default_factory=list)
    n_ancilla: int = 0
    metadata: dict = field(default_factory=dict)

    # -- construction ---------------------------------------------------
    @property
    def n_data(self) -> int:
        return self.n_qubits - self.n_ancilla

    def add(self, name: str, qubits: Sequence[int], param: float = 0.0) -> "Circuit":
        qs = tuple(int(q) for q in qubits)
        if len(set(qs)) != len(qs):
            raise ValueError(f"repeated qubit in {name}{qs}")
        for q in qs:
            if not 0 <= q < self.n_qubits:
                raise ValueError(f"qubit {q} out of range for width {self.n_qubits}")
        self.gates.append(Gate(name, qs, float(param)))
        return self

    def h(self, q: int) -> "Circuit":
        return self.add("h", (q,))

    def x(self, q: int) -> "Circuit":
        return self.add("x", (q,))

    def z(self, q: int) -> "Circuit":
        return self.add("z", (q,))

    def s(self, q: int) -> "Circuit":
        return self.add("s", (q,))

    def sdg(self, q: int) -> "Circuit":
        return self.add("sdg", (q,))

    def ry(self, q: int, theta: float) -> "Circuit":
        return self.add("ry", (q,), theta)

    def rz(self, q: int, theta: float) -> "Circuit":
        return self.add("rz", (q,), theta)

    def cx(self, c: int, t: int) -> "Circuit":
        return self.add("cx", (c, t))

    def cz(self, a: int, b: int) -> "Circuit":
        return self.add("cz", (a, b))

    def swap(self, a: int, b: int) -> "Circuit":
        return self.add("swap", (a, b))

    def mcx(self, controls: Sequence[int], target: int) -> "Circuit":
        return self.add("mcx", (*controls, target))

    def mcry(self, controls: Sequence[int], target: int, theta: float) -> "Circuit":
        return self.add("mcry", (*controls, target), theta)

    # -- combinators ----------------------------------------------------
    def copy(self) -> "Circuit":
        return Circuit(self.n_qubits, list(self.gates), self.n_ancilla, dict(self.metadata))

    def compose(self, other: "Circuit", qubit_map: Sequence[int] | None = None) -> "Circuit":
        """Append ``other`` onto this circuit, relabelling its wires."""
        if qubit_map is None:
            if other.n_qubits > self.n_qubits:
                raise ValueError("cannot compose a wider circuit without a qubit map")
            qubit_map = range(other.n_qubits)
        qubit_map = list(qubit_map)
        for g in other.gates:
            self.add(g.name, tuple(qubit_map[q] for q in g.qubits), g.param)
        return self

    def inverse(self) -> "Circuit":
        """Inverse circuit.  Only gates used by this package are supported."""
        inv = Circuit(self.n_qubits, [], self.n_ancilla, dict(self.metadata))
        flip = {"s": "sdg", "sdg": "s", "t": "tdg", "tdg": "t"}
        for g in reversed(self.gates):
            if g.name in {"ry", "rz", "mcry"}:
                inv.add(g.name, g.qubits, -g.param)
            elif g.name in flip:
                inv.add(flip[g.name], g.qubits, g.param)
            else:
                inv.add(g.name, g.qubits, g.param)
        return inv


# ----------------------------------------------------------------------
# Lowering
# ----------------------------------------------------------------------

def _toffoli(out: Circuit, a: int, b: int, t: int) -> None:
    """Standard ancilla-free Toffoli: exactly 6 CX."""
    out.h(t)
    out.cx(b, t)
    out.add("tdg", (t,))
    out.cx(a, t)
    out.add("t", (t,))
    out.cx(b, t)
    out.add("tdg", (t,))
    out.cx(a, t)
    out.add("t", (b,))
    out.add("t", (t,))
    out.h(t)
    out.cx(a, b)
    out.add("t", (a,))
    out.add("tdg", (b,))
    out.cx(a, b)


def _mcx_vchain(out: Circuit, controls: Sequence[int], target: int, anc: Sequence[int]) -> None:
    """Multi-controlled X via a clean-ancilla V-chain.

    Uses ``len(controls) - 2`` clean ancillas and ``2*len(controls) - 3``
    Toffolis; every ancilla is returned to |0>.
    """
    c = list(controls)
    if len(c) == 0:
        out.x(target)
        return
    if len(c) == 1:
        out.cx(c[0], target)
        return
    if len(c) == 2:
        _toffoli(out, c[0], c[1], target)
        return
    k = len(c) - 2
    a = list(anc[:k])
    # compute
    _toffoli(out, c[0], c[1], a[0])
    for j in range(k - 1):
        _toffoli(out, c[j + 2], a[j], a[j + 1])
    _toffoli(out, c[-1], a[-1], target)
    # uncompute
    for j in reversed(range(k - 1)):
        _toffoli(out, c[j + 2], a[j], a[j + 1])
    _toffoli(out, c[0], c[1], a[0])


def lowering_ancillas_needed(circ: Circuit) -> int:
    """Peak clean-ancilla demand of :func:`lower_to_cx` on ``circ``."""
    need = 0
    for g in circ.gates:
        if g.name in MULTI_GATES:
            need = max(need, len(g.controls) - 2)
    return max(need, 0)


def lower_to_cx(circ: Circuit) -> Circuit:
    """Rewrite ``circ`` over the native gate set {1-qubit rotations, CX}.

    The returned circuit may be wider than the input: lowering allocates clean
    ancillas for multi-controlled gates.  ``n_ancilla`` of the result counts the
    input's declared ancillas plus the lowering ancillas.
    """
    extra = lowering_ancillas_needed(circ)
    out = Circuit(circ.n_qubits + extra, [], circ.n_ancilla + extra, dict(circ.metadata))
    anc = list(range(circ.n_qubits, circ.n_qubits + extra))
    for g in circ.gates:
        if g.name in ONE_QUBIT_GATES:
            out.add(g.name, g.qubits, g.param)
        elif g.name == "cx":
            out.cx(*g.qubits)
        elif g.name == "cz":
            a, b = g.qubits
            out.h(b)
            out.cx(a, b)
            out.h(b)
        elif g.name == "swap":
            a, b = g.qubits
            out.cx(a, b)
            out.cx(b, a)
            out.cx(a, b)
        elif g.name == "mcx":
            _mcx_vchain(out, g.controls, g.target, anc)
        elif g.name == "mcry":
            ctrls, t, th = g.controls, g.target, g.param
            if len(ctrls) == 0:
                out.ry(t, th)
            else:
                out.ry(t, th / 2.0)
                _mcx_vchain(out, ctrls, t, anc)
                out.ry(t, -th / 2.0)
                _mcx_vchain(out, ctrls, t, anc)
        else:
            raise ValueError(f"cannot lower gate {g.name!r}")
    out.metadata["lowering_ancillas"] = extra
    return out


# ----------------------------------------------------------------------
# Routing
# ----------------------------------------------------------------------

def _bfs_path(adj: dict[int, set[int]], src: int, dst: int) -> list[int]:
    if src == dst:
        return [src]
    prev = {src: src}
    q = deque([src])
    while q:
        u = q.popleft()
        for v in adj.get(u, ()):  # deterministic: adjacency built from sorted edges
            if v not in prev:
                prev[v] = u
                if v == dst:
                    path = [v]
                    while path[-1] != src:
                        path.append(prev[path[-1]])
                    return list(reversed(path))
                q.append(v)
    raise ValueError(f"coupling graph is disconnected: no path {src} -> {dst}")


def route(circ: Circuit, edges: Iterable[tuple[int, int]]) -> tuple[Circuit, int, list[int]]:
    """Insert SWAPs so every CX acts on a physical edge.

    Deliberately simple and deliberately honest: virtual qubit ``v`` starts on
    physical qubit ``v``; for each illegal CX we walk the control along a
    shortest path.  Returns ``(routed, n_swaps, final_layout)`` where
    ``final_layout[v]`` is the physical wire holding virtual qubit ``v``.
    """
    adj: dict[int, set[int]] = {}
    edge_set: set[frozenset[int]] = set()
    for a, b in edges:
        a, b = int(a), int(b)
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)
        edge_set.add(frozenset((a, b)))
    # sort adjacency for reproducible BFS
    adj = {k: set(sorted(v)) for k, v in adj.items()}

    layout = list(range(circ.n_qubits))  # virtual -> physical
    out = Circuit(circ.n_qubits, [], circ.n_ancilla, dict(circ.metadata))
    n_swaps = 0
    for g in circ.gates:
        if g.name in ONE_QUBIT_GATES:
            out.add(g.name, (layout[g.qubits[0]],), g.param)
            continue
        if g.name != "cx":
            raise ValueError("route() expects a circuit lowered to {1q, cx}")
        cv, tv = g.qubits
        cp, tp = layout[cv], layout[tv]
        if frozenset((cp, tp)) not in edge_set:
            path = _bfs_path(adj, cp, tp)
            for i in range(len(path) - 2):
                u, v = path[i], path[i + 1]
                out.add("swap", (u, v))
                n_swaps += 1
                iu = layout.index(u)
                iv = layout.index(v)
                layout[iu], layout[iv] = layout[iv], layout[iu]
            cp = layout[cv]
        out.cx(cp, layout[tv])
    routed = Circuit(out.n_qubits, [], out.n_ancilla, dict(out.metadata))
    for g in out.gates:  # expand SWAPs last so swap_count is recorded separately
        if g.name == "swap":
            a, b = g.qubits
            routed.cx(a, b)
            routed.cx(b, a)
            routed.cx(a, b)
        else:
            routed.add(g.name, g.qubits, g.param)
    routed.metadata["swap_count"] = n_swaps
    routed.metadata["final_layout"] = list(layout)
    return routed, n_swaps, layout


def is_hardware_legal(circ: Circuit, edges: Iterable[tuple[int, int]]) -> bool:
    edge_set = {frozenset((int(a), int(b))) for a, b in edges}
    for g in circ.gates:
        if g.name in TWO_QUBIT_GATES:
            if frozenset(g.qubits) not in edge_set:
                return False
        elif g.name in MULTI_GATES:
            return False
    return True


# ----------------------------------------------------------------------
# Resource accounting (always from an emitted circuit)
# ----------------------------------------------------------------------

def resources(circ: Circuit) -> dict:
    """Resource metrics read off the instruction list."""
    two_q = [g for g in circ.gates if g.name in TWO_QUBIT_GATES]
    if any(g.name in MULTI_GATES for g in circ.gates):
        raise ValueError("resources() requires a lowered circuit")
    depth_layer = [0] * circ.n_qubits
    depth2_layer = [0] * circ.n_qubits
    for g in circ.gates:
        lvl = max(depth_layer[q] for q in g.qubits) + 1
        for q in g.qubits:
            depth_layer[q] = lvl
        if g.name in TWO_QUBIT_GATES:
            lvl2 = max(depth2_layer[q] for q in g.qubits) + 1
            for q in g.qubits:
                depth2_layer[q] = lvl2
    return {
        "n_qubits": circ.n_qubits,
        "n_ancilla": circ.n_ancilla,
        "gate_count": len(circ.gates),
        "two_qubit_gate_count": len(two_q),
        "cx_count": sum(1 for g in two_q if g.name == "cx"),
        "swap_count": int(circ.metadata.get("swap_count", 0)),
        "depth": max(depth_layer) if depth_layer else 0,
        "two_qubit_depth": max(depth2_layer) if depth2_layer else 0,
    }


def to_qasm(circ: Circuit) -> str:
    """OpenQASM 2.0 export for a lowered circuit (interchange only)."""
    lines = ["OPENQASM 2.0;", 'include "qelib1.inc";', f"qreg q[{circ.n_qubits}];"]
    for g in circ.gates:
        qs = ",".join(f"q[{q}]" for q in g.qubits)
        if g.name in {"ry", "rz"}:
            lines.append(f"{g.name}({g.param:.17g}) {qs};")
        elif g.name == "swap":
            lines.append(f"swap {qs};")
        else:
            lines.append(f"{g.name} {qs};")
    return "\n".join(lines) + "\n"
