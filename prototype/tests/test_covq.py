"""Tests for the claims the audit actually relies on.

Each test corresponds to a statement in
``docs/audits/TASK0AB_NOVELTY_THEOREM_AUDIT_v0.2.md``.  They are written so that a
*false theorem* fails them, not just a broken implementation.
"""

from __future__ import annotations

import itertools
import math

import numpy as np
import pytest

from covq import baselines as bl
from covq import gates as gt
from covq import polytope as pol
from covq import programs as prg
from covq import width as wid
from covq.circuits import is_hardware_legal, lower_to_cx, resources, route
from covq.instances import hypermetric_violator, path_target, random_signed_cat_mixture, star_target
from covq.paulis import (
    all_commute,
    apply_pauli,
    diagonalizing_clifford,
    is_independent,
    random_commuting_paulis,
    z_generators,
)
from covq.qfim import qfim_from_statevector, sld_commutator_defect
from covq.sim import data_statevector, simulate
from covq.topology import build as build_topology

SEEDS = [2026, 3407, 9181]


# -- signed cats and the Fisher convention -----------------------------

@pytest.mark.parametrize("m", [2, 3, 4, 5])
def test_signed_cat_realises_outer_product(m):
    rng = np.random.default_rng(SEEDS[0] + m)
    ps = z_generators(m)
    for _ in range(6):
        s = rng.choice(np.array([-1, 1]), size=m)
        psi = data_statevector(prg.signed_cat_circuit(s, n_qubits=m))
        F = qfim_from_statevector(psi, ps)
        assert np.allclose(F, np.outer(s, s), atol=1e-12)


@pytest.mark.parametrize("m", [3, 4, 5])
def test_unit_diagonal_iff_zero_mean(m):
    """L1: F_ii = 1 - <P_i>^2, so the sector is one condition, not two."""
    rng = np.random.default_rng(SEEDS[1] + m)
    ps = z_generators(m)
    for _ in range(8):
        psi = rng.standard_normal(1 << m) + 1j * rng.standard_normal(1 << m)
        psi /= np.linalg.norm(psi)
        F = qfim_from_statevector(psi, ps)
        mu = prg.pauli_means(psi, ps)
        assert np.allclose(np.diag(F), 1.0 - mu**2, atol=1e-12)


@pytest.mark.parametrize("m", [3, 4])
def test_sld_saturability_holds_identically(m):
    """L6/BP6: commuting Hermitian generators give zero mean Uhlmann curvature."""
    rng = np.random.default_rng(SEEDS[2] + m)
    for ps in (z_generators(m), random_commuting_paulis(m, rng)):
        psi = rng.standard_normal(1 << m) + 1j * rng.standard_normal(1 << m)
        psi /= np.linalg.norm(psi)
        assert sld_commutator_defect(psi, ps) < 1e-12


# -- Clifford normalisation --------------------------------------------

@pytest.mark.parametrize("m", [2, 3, 4, 5])
def test_diagonalizing_clifford_maps_generators_to_signed_z(m):
    rng = np.random.default_rng(SEEDS[0] * m)
    for _ in range(3):
        ps = random_commuting_paulis(m, rng)
        assert all_commute(ps) and is_independent(ps)
        norm = diagonalizing_clifford(ps)
        n = ps.n
        v = rng.standard_normal(1 << n) + 1j * rng.standard_normal(1 << n)
        v /= np.linalg.norm(v)
        idx = np.arange(1 << n)
        for i in range(m):
            w = simulate(norm.circuit, apply_pauli(simulate(norm.circuit.inverse(), v), ps, i))
            sgn = 1.0 - 2.0 * (np.bitwise_count(idx & (1 << norm.perm[i])) & 1)
            assert np.allclose(w, int(norm.sign_out[i]) * sgn * v, atol=1e-9)


# -- the polytope and its certificates ---------------------------------

def test_reference_infeasible_instance():
    """C4's ground truth: PSD, unit diagonal, and provably outside Q_3."""
    inst = hypermetric_violator(3, 3, 0.4)
    assert np.allclose(np.linalg.eigvalsh(inst.F), [0.2, 1.4, 1.4], atol=1e-9)
    dec = pol.exact_decompose(inst.F)
    assert not dec.feasible
    assert dec.certificate is not None and dec.certificate.verify(inst.F)
    b = np.array([1.0, 1.0, 1.0])
    assert abs(float(b @ inst.F @ b) - 0.6) < 1e-9


@pytest.mark.parametrize("m", [3, 4, 5])
def test_hypermetric_inequality_holds_on_every_feasible_target(m):
    """b^T F b = Var(sum b_i P_i) >= 1 whenever sum b_i is odd."""
    rng = np.random.default_rng(SEEDS[1] * m)
    for _ in range(5):
        F = random_signed_cat_mixture(m, 4, rng).F
        for b in itertools.product([-1, 0, 1], repeat=m):
            if sum(b) % 2 == 0:
                continue
            bv = np.array(b, dtype=float)
            assert float(bv @ F @ bv) >= 1.0 - 1e-9


@pytest.mark.parametrize("m", [3, 4, 5])
def test_exact_lp_and_column_generation_agree(m):
    rng = np.random.default_rng(SEEDS[2] * m)
    for _ in range(6):
        A = rng.uniform(-0.8, 0.8, size=(m, m))
        F = (A + A.T) / 2
        np.fill_diagonal(F, 1.0)
        lp = pol.exact_decompose(F)
        cg = pol.column_generation(F, rng=np.random.default_rng(7))
        assert bool(lp.feasible) == bool(cg.feasible)


@pytest.mark.parametrize("m", [4, 5])
def test_caratheodory_reduction_preserves_the_matrix(m):
    rng = np.random.default_rng(SEEDS[0] + 11 * m)
    F = random_signed_cat_mixture(m, 6, rng).F
    dec = pol.exact_decompose(F)
    red = pol.caratheodory_reduce(dec)
    assert red.support <= m * (m - 1) // 2 + 1
    assert np.allclose(red.matrix(), F, atol=1e-9)


def test_frank_wolfe_gap_certifies_its_own_error():
    rng = np.random.default_rng(SEEDS[1])
    F = random_signed_cat_mixture(5, 4, rng).F
    fw = pol.frank_wolfe(F, iters=200)
    assert fw.info["final_duality_gap"] >= -1e-12
    assert fw.residual <= 1e-6


# -- the width-two theorem ---------------------------------------------

@pytest.mark.parametrize("m", [3, 4, 5])
def test_width2_matches_brute_force_enumeration(m):
    """The theorem itself: Edmonds separation vs full extreme-point enumeration."""
    rng = np.random.default_rng(SEEDS[2] + m)
    for _ in range(15):
        A = rng.uniform(-0.7, 0.7, size=(m, m))
        F = (A + A.T) / 2
        np.fill_diagonal(F, 1.0)
        edmonds = wid.decompose_width2(F)
        brute = wid.decompose_width_k(F, 2, force_enumeration=True)
        assert bool(edmonds.feasible) == bool(brute.feasible)
        if edmonds.feasible:
            assert np.allclose(edmonds.matrix(), F, atol=1e-8)


@pytest.mark.parametrize("m", [4, 5, 6])
def test_path_threshold_is_exactly_one_half(m):
    for c in (0.45, 0.49, 0.50):
        assert wid.decompose_width2(path_target(m, c).F).feasible
    for c in (0.51, 0.55):
        assert not wid.decompose_width2(path_target(m, c).F).feasible


@pytest.mark.parametrize("m", [3, 5])
def test_blossom_facets_are_active_not_redundant(m):
    """Degree constraints alone would permit c <= 1/2; the odd set forces c <= 1/3.

    The whole family lies in Q_m, so this is what proves the characterisation says
    more than degree counting does.
    """
    for c, expected in ((0.30, True), (1 / 3, True), (0.34, False), (0.40, False)):
        F = np.eye(m)
        for i in range(3):
            for j in range(3):
                if i != j:
                    F[i, j] = c
        assert pol.exact_decompose(F).feasible, "family must stay inside Q_m"
        d = wid.decompose_width2(F)
        assert bool(d.feasible) is expected
        if not expected:
            assert d.violation is not None and d.violation.kind == "odd_set"


def test_star_threshold():
    m = 5
    assert wid.decompose_width2(star_target(m, 0.9 / (m - 1)).F).feasible
    assert not wid.decompose_width2(star_target(m, 1.4 / (m - 1)).F).feasible


@pytest.mark.parametrize("m", [5, 6])
def test_hardware_native_width2_rejects_off_graph_support(m):
    F = np.eye(m)
    for i in range(m - 1):
        F[i, i + 1] = F[i + 1, i] = 0.4
    assert wid.decompose_width2_hardware(F, build_topology("line", m)).feasible
    grid = wid.decompose_width2_hardware(F, build_topology("square_grid", m))
    assert not grid.feasible
    assert grid.violation is not None and grid.violation.kind == "off_graph_support"


@pytest.mark.parametrize("m", [4, 5])
def test_width2_schedule_realises_target_from_emitted_circuits(m):
    F = np.eye(m)
    for i in range(m - 1):
        F[i, i + 1] = F[i + 1, i] = 0.4
    d = wid.decompose_width2(F)
    assert d.feasible
    prog = prg.labelled_schedule_from_width(d)
    assert np.allclose(prog.realised_qfim(z_generators(m)), F, atol=1e-10)
    audit = prg.audit_branch_width(prog)
    assert audit["all_declared_widths_hold"]
    assert audit["max_declared_width"] <= 2


# -- program semantics -------------------------------------------------

def test_coherent_flag_identity_with_nonzero_means():
    """C2b: F_Psi = sum p_r F_r + Cov(mu_r), tested where the two sides differ."""
    m = 3
    ps = z_generators(m)
    psi = np.zeros(1 << (m + 1), dtype=complex)
    psi[0] = np.sqrt(0.5)                          # |flag=0>|000>
    psi[(1 << m) | ((1 << m) - 1)] = np.sqrt(0.5)  # |flag=1>|111>
    F = prg._qfim_on_data(psi, ps, m + 1)
    mus = np.array([np.ones(m), -np.ones(m)])
    mbar = 0.5 * (mus[0] + mus[1])
    cov = 0.5 * np.outer(mus[0], mus[0]) + 0.5 * np.outer(mus[1], mus[1]) - np.outer(mbar, mbar)
    assert np.allclose(F, cov, atol=1e-12)         # branch QFIMs are zero here
    assert np.linalg.norm(F) > 1e-6                # ... and the program QFIM is not


def test_conditional_width_is_not_a_resource_for_coherent_flags():
    """The loophole: width-1 conditional branches realising an arbitrary target."""
    m = 4
    rng = np.random.default_rng(SEEDS[0])
    inst = random_signed_cat_mixture(m, 4, rng)
    dec = pol.caratheodory_reduce(pol.exact_decompose(inst.F))
    ps = z_generators(m)
    prog = prg.coherent_sign_purification(dec.signs, dec.weights)
    rep = prg.flag_resource_report(prog, ps)
    assert np.allclose(np.array(rep["qfim_coherent_flag"]), inst.F, atol=1e-10)
    assert rep["conditional_data_width"] == 1
    assert np.allclose(np.array(rep["qfim_dephased_flag"]), 0.0, atol=1e-10)
    assert rep["coherence_is_load_bearing"]


def test_literal_flagged_cat_survives_dephasing():
    m = 4
    rng = np.random.default_rng(SEEDS[1])
    inst = random_signed_cat_mixture(m, 4, rng)
    dec = pol.caratheodory_reduce(pol.exact_decompose(inst.F))
    ps = z_generators(m)
    prog = prg.coherent_flag_program_literal(dec.signs, dec.weights)
    rep = prg.flag_resource_report(prog, ps)
    assert np.allclose(np.array(rep["qfim_coherent_flag"]), inst.F, atol=1e-9)
    assert np.allclose(np.array(rep["qfim_dephased_flag"]), inst.F, atol=1e-9)
    assert not rep["coherence_is_load_bearing"]


def test_convexity_gap_is_zero_for_cats_and_maximal_for_the_control():
    m = 4
    ps = z_generators(m)
    rng = np.random.default_rng(SEEDS[2])
    S = rng.choice(np.array([-1, 1]), size=(4, m))
    S[:, 0] = np.abs(S[:, 0])
    w = rng.dirichlet(np.ones(4))
    prog = prg.labelled_schedule_from_signs(S, w)
    avg = prog.realised_qfim(ps)
    mixed = prg.unlabelled_mixed_qfi(prg.branch_states(prog), ps)
    assert np.allclose(avg, mixed, atol=1e-10)  # orthogonal branches: label is free

    s = np.array([1] + [-1] * (m - 1))
    even = prg.signed_cat_circuit(s, n_qubits=m)
    odd = prg.signed_cat_circuit(s, n_qubits=m)
    odd.z(0)
    states = [(0.5, data_statevector(even)), (0.5, data_statevector(odd))]
    branch_avg = sum(wt * qfim_from_statevector(v, ps) for wt, v in states)
    assert np.allclose(prg.unlabelled_mixed_qfi(states, ps), 0.0, atol=1e-10)
    assert np.trace(branch_avg) == pytest.approx(m)


# -- emission, lowering, routing ---------------------------------------

def test_lowering_preserves_semantics():
    m = 4
    rng = np.random.default_rng(SEEDS[0])
    S = rng.choice(np.array([-1, 1]), size=(3, m))
    S[:, 0] = np.abs(S[:, 0])
    w = rng.dirichlet(np.ones(3))
    for prog in (prg.single_pure_state_program(S, w), prg.coherent_flag_program(S, w)):
        hi = simulate(prog.circuit)
        low = lower_to_cx(prog.circuit)
        lo = simulate(low)
        extra = low.n_qubits - prog.circuit.n_qubits
        if extra:
            lo = lo.reshape(1 << extra, -1)[0]
        assert abs(abs(complex(np.vdot(hi, lo))) - 1.0) < 1e-9


@pytest.mark.parametrize("topo", ["line", "square_grid", "heavy_hex", "modular_two_cluster"])
def test_routing_is_legal(topo):
    m = 6
    edges = build_topology(topo, m)
    prog = prg.labelled_schedule_from_signs(np.ones((1, m), dtype=np.int8), np.array([1.0]))
    for br in prog.branches:
        routed, _, _ = route(lower_to_cx(br.circuit), edges)
        assert is_hardware_legal(routed, edges)


def test_sparse_state_beats_dense_baseline():
    """Carathéodory sparsity is worth something -- but both counts are upper bounds."""
    m = 5
    rng = np.random.default_rng(SEEDS[1])
    dec = pol.caratheodory_reduce(
        pol.exact_decompose(random_signed_cat_mixture(m, 4, rng).F))
    sparse = prg.single_pure_state_program(dec.signs, dec.weights).emitted()["cx_count"]
    dense = bl.dense_state_preparation(dec.signs, dec.weights).emitted()["cx_count"]
    assert sparse < dense


# -- the gate suite ----------------------------------------------------

def test_gate_suite_has_no_failures():
    from covq.instances import task0_suite

    report = gt.run_all(task0_suite((3, 4)), stop_on_failure=False)
    failures = [g["name"] for g in report["gates"] if g["status"] == "FAIL"]
    assert not failures, failures
    statuses = {g["name"]: g["status"] for g in report["gates"]}
    # C10 and C12 state facts that hold for every matrix / every cat schedule.
    assert statuses["C10_downstream_stability"] == "MEASURED"
    assert statuses["C12_convexity_gap"] == "MEASURED"


# -- relationship to the known k-producibility bound --------------------

def _toth_bound(m: int, k: int) -> int:
    s, r = m // k, m - (m // k) * k
    return s * k * k + r * r


@pytest.mark.parametrize("m,k", [(3, 1), (3, 2), (3, 3), (4, 2), (4, 3),
                                 (5, 2), (5, 3), (6, 2), (6, 3), (7, 2), (7, 3)])
def test_support_function_at_b_ones_reproduces_k_producibility_bound(m, k):
    """The known scalar bound is one direction of our polytope's support function.

    For k-producible states the QFI in the collective direction is bounded by
    ``floor(m/k) k^2 + r^2``.  In our normalisation that is
    ``max { 1^T F 1 : F in Q^prog_{m,k} }``, so the literature bound is the value
    of the support function of ``Q^prog_{m,k}`` at ``b = 1`` -- and nothing more.
    The characterisation gives every other direction too.
    """
    A, _ = wid.width_k_generators(m, k)
    best = max(m + 2.0 * float(col.sum()) for col in A.T)
    assert best == pytest.approx(_toth_bound(m, k))


@pytest.mark.parametrize("m", [3, 4, 5, 6, 8, 10, 12])
def test_matching_polytope_gives_the_same_bound_at_k_two_without_enumeration(m):
    """At k = 2 the matching polytope reproduces it in closed form: m + 2*floor(m/2)."""
    assert m + 2 * (m // 2) == _toth_bound(m, 2)


# -- the information-floor compiler ------------------------------------

from covq.floor import common_mode_contract, emit_floor_program, information_floor_compile  # noqa: E402
from covq.topology import build as _topo  # noqa: E402


def _threshold(A, m, edges, lo=1.0, hi=3.0, iters=44):
    """Largest gamma with gamma*I feasible, by bisection."""
    d = A.shape[1]
    for _ in range(iters):
        mid = (lo + hi) / 2
        if information_floor_compile(A, mid * np.eye(d), m, edges).status == "feasible":
            lo = mid
        else:
            hi = mid
    return lo


@pytest.mark.parametrize("m", [4, 5, 6])
@pytest.mark.parametrize("gamma", [1.1, 1.3])
def test_common_mode_cost_is_m_gamma_minus_one_over_two(m, gamma):
    """The analytic compilation curve: cost* = m(gamma-1)/2 at unit edge costs."""
    u, G, predicted = common_mode_contract(m, gamma)
    edges = _topo("all_to_all", m)
    res = information_floor_compile(u, G, m, edges)
    assert res.status == "feasible"
    assert res.cost == pytest.approx(predicted, abs=1e-8)
    assert res.lower_bound == pytest.approx(res.cost, abs=1e-8)  # zero optimality gap


@pytest.mark.parametrize("m", [4, 5, 6])
def test_emitted_floor_program_actually_meets_the_contract(m):
    u, G, _ = common_mode_contract(m, 1.3)
    edges = _topo("all_to_all", m)
    res = information_floor_compile(u, G, m, edges)
    prog, dec, ident = emit_floor_program(res, edges)
    F = prog.realised_qfim(z_generators(m))
    assert (u.T @ F @ u).item() >= 1.3 - 1e-8
    # the objective IS the expected per-shot Bell-pair bill
    assert ident["cost_identity_residual"] < 1e-9
    assert ident["n_branches"] <= ident["caratheodory_bound"]


@pytest.mark.parametrize("m", [4, 5])
def test_floor_above_the_support_value_is_rejected_with_a_certificate(m):
    gamma_max = 1 + 2 * (m // 2) / m          # = support function at b = 1, over m
    u, G, _ = common_mode_contract(m, gamma_max + 0.05)
    res = information_floor_compile(u, G, m, _topo("all_to_all", m))
    assert res.status == "infeasible"
    assert res.certificate["kind"] in ("lp_infeasible_under_valid_cuts",
                                       "farkas_dual_ray")
    # and the boundary itself is feasible
    u2, G2, _ = common_mode_contract(m, gamma_max)
    assert information_floor_compile(u2, G2, m, _topo("all_to_all", m)).status == "feasible"


@pytest.mark.parametrize("k", [2, 3, 4, 5, 6, 7])
def test_odd_even_collective_mode_ceiling(k):
    """Width-2 ceiling on a size-k collective mode: 2 if k even, 2 - 1/k if k odd.

    The odd-size deficit is the blossom inequality expressed as information --
    an odd collective mode always leaves one qubit unmatched and therefore
    provably cannot reach the pair-entanglement ceiling.
    """
    m = 8
    u = np.zeros((m, 1))
    u[:k, 0] = 1 / np.sqrt(k)
    predicted = 2.0 if k % 2 == 0 else 2.0 - 1.0 / k
    assert _threshold(u, m, _topo("all_to_all", m)) == pytest.approx(predicted, abs=1e-4)


def test_overlapping_modes_force_a_genuine_tradeoff():
    """Two 3-modes sharing a qubit bind well below the independent threshold."""
    m = 6
    u1 = np.zeros(m); u1[[0, 1, 2]] = 1 / np.sqrt(3)
    u2 = np.zeros(m); u2[[2, 3, 4]] = 1 / np.sqrt(3)   # shares qubit 2
    u3 = np.zeros(m); u3[[3, 4, 5]] = 1 / np.sqrt(3)   # disjoint from u1
    edges = _topo("all_to_all", m)
    single = _threshold(u1.reshape(m, 1), m, edges)
    disjoint = _threshold(np.column_stack([u1, u3]), m, edges)
    overlapping = _threshold(np.column_stack([u1, u2]), m, edges)
    assert single == pytest.approx(5 / 3, abs=1e-4)
    assert disjoint == pytest.approx(single, abs=1e-4)   # no interaction
    assert overlapping < single - 0.3                    # a real trade-off


def test_multi_direction_contract_exercises_the_semidefinite_constraint():
    """With d >= 2 the Loewner constraint is genuinely active, not linear."""
    m = 6
    u1 = np.zeros(m); u1[[0, 1, 2]] = 1 / np.sqrt(3)
    u2 = np.zeros(m); u2[[2, 3, 4]] = 1 / np.sqrt(3)
    A = np.column_stack([u1, u2])
    res = information_floor_compile(A, 1.2 * np.eye(2), m, _topo("all_to_all", m))
    assert res.status == "feasible"
    assert res.n_cuts > 1, "a d>=2 contract should need more than one cutting plane"
    assert res.gap == pytest.approx(0.0, abs=1e-8)
    assert res.slack_min_eig >= -1e-8


def test_product_probe_reaches_only_gamma_one():
    """The floor is non-trivial: no width-1 program clears gamma > 1."""
    m = 5
    ps = z_generators(m)
    from covq.baselines import product_probe
    F = product_probe(m).realised_qfim(ps)
    u = np.ones((m, 1)) / np.sqrt(m)
    assert (u.T @ F @ u).item() == pytest.approx(1.0, abs=1e-12)


# -- the explicit conic dual and separation oracle ----------------------

from covq.floor import matching_separation_oracle  # noqa: E402


@pytest.mark.parametrize("m", [4, 5, 6])
@pytest.mark.parametrize("gamma", [1.1, 1.3, 1.5])
def test_dual_certificate_is_valid_and_tight(m, gamma):
    """Weak duality must hold; here it is attained, so the bound is a proof."""
    u, G, _ = common_mode_contract(m, gamma)
    res = information_floor_compile(u, G, m, _topo("all_to_all", m))
    assert res.status == "feasible"
    assert res.dual is not None
    check = res.dual.verify()
    assert check["valid"], check
    assert check["Y_min_eigenvalue"] >= -1e-7          # Y >= 0
    assert check["worst_edge_violation"] <= 1e-7       # dual feasibility
    assert res.dual.value <= res.cost + 1e-7           # weak duality
    assert res.dual.value == pytest.approx(res.cost, abs=1e-7)   # attained


def test_dual_certificate_on_a_multi_direction_contract():
    m = 6
    u1 = np.zeros(m); u1[[0, 1, 2]] = 1 / np.sqrt(3)
    u2 = np.zeros(m); u2[[2, 3, 4]] = 1 / np.sqrt(3)
    res = information_floor_compile(np.column_stack([u1, u2]), 1.2 * np.eye(2), m,
                                    _topo("all_to_all", m))
    assert res.status == "feasible" and res.dual is not None
    assert res.dual.verify()["valid"]
    assert res.dual.value == pytest.approx(res.cost, abs=1e-7)


@pytest.mark.parametrize("m", [4, 5])
def test_farkas_ray_certifies_infeasibility(m):
    """An infeasible floor must come back with a verifiable dual ray."""
    gamma_max = 1 + 2 * (m // 2) / m
    u, G, _ = common_mode_contract(m, gamma_max + 0.05)
    res = information_floor_compile(u, G, m, _topo("all_to_all", m))
    assert res.status == "infeasible"
    assert res.certificate["kind"] == "farkas_dual_ray"
    assert res.dual is not None and res.dual.kind == "farkas"
    check = res.dual.verify()
    assert check["valid"], check
    assert check["recomputed_objective"] > 1e-9   # strictly positive => no primal point


def test_separation_oracle_finds_blossom_degree_and_nothing():
    m = 5
    edges = _topo("all_to_all", m)
    E = sorted({tuple(sorted(e)) for e in edges})
    idx = {e: k for k, e in enumerate(E)}

    t = np.zeros(len(E))
    for e in [(0, 1), (0, 2), (1, 2)]:
        t[idx[e]] = 0.4                      # 1.2 > 1: blossom on {0,1,2}
    kind, S, lhs, rhs = matching_separation_oracle(t, m, edges)
    assert kind == "odd_set" and set(S) == {0, 1, 2} and lhs > rhs

    t = np.zeros(len(E))
    t[idx[(0, 1)]] = t[idx[(0, 2)]] = 0.7    # 1.4 > 1: degree at qubit 0
    kind, S, lhs, rhs = matching_separation_oracle(t, m, edges)
    assert kind == "degree" and S == (0,) and lhs > rhs

    t = np.zeros(len(E))
    t[idx[(0, 1)]] = t[idx[(2, 3)]] = 0.5    # a genuine matching point
    assert matching_separation_oracle(t, m, edges) is None


# -- manuscript v0.2 claims, independently checked ----------------------

from covq.floor import branch_matrix, price_branch, shot_scaled_floor_compile  # noqa: E402


@pytest.mark.parametrize("topo,m", [("line", 5), ("line", 6), ("square_grid", 6), ("ring", 6)])
def test_corollary_5_5_bipartite_odd_sets_are_redundant(topo, m):
    """On bipartite H the degree constraints alone characterise Q^lab_{H,2}."""
    rng = np.random.default_rng(11 + m)
    edges = _topo(topo, m)
    for _ in range(25):
        F = np.eye(m)
        for (i, j) in edges:
            if rng.random() < 0.8:
                F[i, j] = F[j, i] = rng.uniform(-1, 1)
        x = np.abs(F).copy(); np.fill_diagonal(x, 0.0)
        degree_only = all(x[i].sum() <= 1 + 1e-9 for i in range(m))
        full = wid.decompose_width2_hardware(F, edges).feasible
        assert degree_only == full


@pytest.mark.parametrize("name,EG,has_partition", [
    ("two disjoint triangles", [(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5)], True),
    ("prism", [(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5), (0, 3), (1, 4), (2, 5)], True),
    ("C6, triangle-free", [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)], False),
    ("triangle plus path", [(0, 1), (0, 2), (1, 2), (3, 4), (4, 5)], False),
    ("K4 plus an edge", [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3), (4, 5)], False),
])
def test_theorem_7_1_width_three_reduction(name, EG, has_partition):
    """max L_G over width-3 programs hits 3q exactly when G partitions into triangles."""
    m, q = 6, 2
    pos = {p: k for k, p in enumerate(pol.pairs(m))}
    A3, _ = wid.width_k_generators(m, 3)
    idx = [pos[tuple(sorted(e))] for e in EG]
    best = max(float(col[idx].sum()) for col in A3.T)
    assert (abs(best - 3 * q) < 1e-9) == has_partition, (name, best)


@pytest.mark.parametrize("m", [4, 5, 6])
def test_equation_65_branch_separation_identity(m):
    """max_{M,sigma}(<B,Y> - c(M)) = tr Y - c0 + max-weight matching on 2|Y_ij| - c_ij."""
    rng = np.random.default_rng(7 * m)
    E = [tuple(sorted(e)) for e in _topo("all_to_all", m)]
    for _ in range(3):
        B = rng.standard_normal((m, m))
        Y = B @ B.T / m
        c0 = float(rng.uniform(0.1, 1.0))
        ce = {e: float(rng.uniform(0, 0.6)) for e in E}

        def all_matchings(av):
            if not av:
                yield []
                return
            u = av[0]
            yield from ([] + r for r in all_matchings(av[1:]))
            for k in range(1, len(av)):
                v = av[k]
                for r in all_matchings(av[1:k] + av[k + 1:]):
                    yield [tuple(sorted((u, v)))] + r

        brute = max(
            float(np.trace(Y)) + sum(2 * abs(Y[i, j]) for i, j in M)
            - (c0 + sum(ce[e] for e in M))
            for M in all_matchings(list(range(m))))
        _, _, formula = price_branch(Y, m, E, c0, ce)
        assert brute == pytest.approx(formula, abs=1e-9)


@pytest.mark.parametrize("m", [4, 5, 6])
def test_problem_8_1_strong_duality_and_equation_63(m):
    rng = np.random.default_rng(3407 + m)
    B = rng.standard_normal((m, m))
    G = B @ B.T / m * 0.8
    res = shot_scaled_floor_compile(G, m, _topo("all_to_all", m), c0=1.0)
    assert res["status"] == "solved"
    assert res["cost"] == pytest.approx(res["dual_bound"], abs=1e-8)     # Thm 8.2
    assert res["cost"] <= res["product_probe_bound"] + 1e-8              # Eq (63)
    assert res["floor_slack_min_eig"] >= -1e-8                           # floor met


@pytest.mark.parametrize("m,c0", [(4, 1.0), (6, 1.0), (8, 1.0), (6, 2.0)])
def test_entanglement_phase_boundary_is_two_c0_over_m(m, c0):
    """Pair entanglement beats the product probe on a common-mode floor iff c_e < 2c0/m.

    A perfect-matching branch has u^T B u = 2, so it needs gamma/2 uses at
    c0 + (m/2)c_e each, against gamma uses at c0 for the product probe.
    """
    u = np.ones((m, 1)) / np.sqrt(m)
    G = 4.0 * (u @ u.T)
    edges = _topo("all_to_all", m)

    def uses_entanglement(ce):
        return shot_scaled_floor_compile(
            G, m, edges, c0=c0,
            costs={tuple(sorted(e)): ce for e in edges})["uses_entanglement"]

    lo, hi = 0.0, 3.0
    for _ in range(30):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if uses_entanglement(mid) else (lo, mid)
    assert lo == pytest.approx(2 * c0 / m, abs=1e-3)


@pytest.mark.parametrize("m", [4, 6, 8])
def test_maximum_advantage_over_product_probe_is_exactly_two(m):
    """As entanglers become free the cost ratio tends to 1/2, the width-2 ceiling."""
    u = np.ones((m, 1)) / np.sqrt(m)
    G = 4.0 * (u @ u.T)
    edges = _topo("all_to_all", m)
    res = shot_scaled_floor_compile(G, m, edges, c0=1.0,
                                    costs={tuple(sorted(e)): 1e-6 for e in edges})
    assert res["cost"] / res["product_probe_bound"] == pytest.approx(0.5, abs=1e-4)
    assert branch_matrix(m, [(0, 1)])[0, 1] == 1.0


# -- frozen canonical controls (manuscript v0.2) ------------------------

from covq.instances import FROZEN_CONTROLS, classify_target, control_minus, control_plus  # noqa: E402


def test_control_plus_is_globally_feasible_but_pair_infeasible():
    """F(+)_ij = 0.4: PD, inside Q_3, outside Q^lab_{3,2} by the blossom inequality."""
    inst = control_plus()
    c = classify_target(inst.F)
    assert c["psd"] is True
    assert c["global_feasible"] is True          # a width-3 program realises it
    assert c["pair_feasible"] is False           # but pair width does not
    assert c["pair_violation"]["kind"] == "odd_set"
    assert set(c["pair_violation"]["subset"]) == {0, 1, 2}
    assert c["pair_violation"]["lhs"] == pytest.approx(1.2)
    assert c["min_width"]["k"] == 3


def test_control_minus_is_positive_definite_but_globally_infeasible():
    """F(-)_ij = -0.4: PD with spectrum {0.2,1.4,1.4}, yet b^T F b = 0.6 < 1."""
    inst = control_minus()
    c = classify_target(inst.F)
    assert c["psd"] is True
    assert np.allclose(sorted(c["eigenvalues"]), [0.2, 1.4, 1.4], atol=1e-9)
    assert c["global_feasible"] is False         # no width realises it
    assert c["global_certificate"]["verified"] is True
    b = np.ones(3)
    assert float(b @ inst.F @ b) == pytest.approx(0.6)
    assert c["min_width"]["k"] is None


def test_the_two_controls_separate_the_three_failure_modes():
    """Neither control alone distinguishes PSD, global feasibility and pair width."""
    cp = classify_target(control_plus().F)
    cm = classify_target(control_minus().F)
    # both are PSD, so PSD cannot be doing the work
    assert cp["psd"] and cm["psd"]
    # both are pair-infeasible, so pair width cannot distinguish them either
    assert not cp["pair_feasible"] and not cm["pair_feasible"]
    # only global feasibility separates them -- which is the point
    assert cp["global_feasible"] != cm["global_feasible"]


@pytest.mark.parametrize("name", sorted(FROZEN_CONTROLS))
def test_frozen_controls_match_their_recorded_expectations(name):
    inst = FROZEN_CONTROLS[name]()
    c = classify_target(inst.F)
    for key in ("psd", "global_feasible", "pair_feasible"):
        assert c[key] == inst.ground_truth[key], (name, key)


def test_edge_and_branch_shot_scaled_agree_where_both_converge():
    """Two formulations of Problem 8.1 must agree; the edge one stalls when entangled.

    Documented limitation, not a silent one: the compact edge SDP is solved here
    by eigenvector cutting planes, which converge only while the entangled edges
    stay inactive.  Where it does converge it matches the branch formulation to
    machine precision, which is the cross-check that matters.
    """
    from covq.floor import shot_scaled_edge_compile

    rng = np.random.default_rng(2026)
    agreed = 0
    for m in (4, 5, 6):
        for topo in ("all_to_all", "line"):
            B = rng.standard_normal((m, m))
            G = B @ B.T / m * 0.8
            edges = _topo(topo, m)
            costs = {tuple(sorted(e)): 1.0 for e in edges}
            rb = shot_scaled_floor_compile(G, m, edges, c0=1.0, costs=costs)
            re_ = shot_scaled_edge_compile(G, m, edges, c0=1.0, costs=costs)
            assert rb["status"] == "solved"
            if re_["status"] != "solved":
                continue
            assert rb["cost"] == pytest.approx(re_["cost"], abs=1e-6)
            agreed += 1
    assert agreed >= 4, "the unentangled regime should converge in both formulations"


# ----------------------------------------------------------------------
# Section 11: measurement compilation.  Gates M1-M3.
# ----------------------------------------------------------------------

def _schedule(m, seed, load=0.8):
    from covq.instances import matching_target
    rng = np.random.default_rng(seed)
    dec = wid.decompose_width2(matching_target(m, rng, load=load).F)
    assert dec.feasible
    return prg.labelled_schedule_from_width(dec)


def test_m1_compiled_local_readout_attains_the_branch_qfim():
    """Section 11.1, made operational.

    A *product of single-qubit* equatorial measurements attains the QFIM of
    every cat-block branch exactly.  Because the label is retained, the
    schedule's classical Fisher matrix then equals its QFIM -- so "realizes the
    QFIM" holds in the operational sense too, for this backend.
    """
    from covq.measurement import compile_schedule_readout

    rng = np.random.default_rng(17)
    for m in (3, 4, 5, 6):
        prog = _schedule(m, seed=m * 13)
        ps = z_generators(m)
        for theta in (np.zeros(m), rng.uniform(-2.0, 2.0, m)):
            rep = compile_schedule_readout(prog, ps, theta)
            assert rep["all_branches_attain"], (m, rep["max_branch_gap"])
            assert rep["all_branches_dominated"]
            assert rep["schedule_gap"] < 1e-9
            assert rep["compiled_total_information"] == pytest.approx(float(m), abs=1e-9)


def test_m1_attainment_survives_wider_cat_blocks():
    """The construction is not special to pair width; k = 2, 3, 4 all attain."""
    from covq.instances import block_diagonal
    from covq.measurement import compile_schedule_readout

    rng = np.random.default_rng(23)
    for k in (2, 3, 4):
        dec = wid.decompose_width_k(block_diagonal(6, k, rng).F, k)
        if not dec.feasible:
            continue
        prog = prg.labelled_schedule_from_width(dec)
        ps = z_generators(6)
        for theta in (np.zeros(6), rng.uniform(-2.0, 2.0, 6)):
            rep = compile_schedule_readout(prog, ps, theta)
            assert rep["schedule_gap"] < 1e-9, (k, rep["schedule_gap"])
            assert rep["all_branches_dominated"]


def test_m2_fixed_readout_collapses_on_a_hyperplane_that_contains_theta_zero():
    """Section 11.1(ii) is load-bearing, not a technicality.

    The unmatched readout attains only off the union of hyperplanes
    ``sum_{i in B} s_i theta_i = 0 mod pi``.  That set is measure zero, but it
    contains ``theta = 0``, and it contains *every* uniform operating point as
    soon as a branch carries a negatively signed edge.  Both are exactly where
    an experiment would choose to sit.

    On those hyperplanes the outcome model is **nonregular**: the score
    vanishes identically, so the first-order Fisher information the compiler
    uses is zero and the signed local phase is not estimable in the regular
    sense.  Probabilities still turn on quadratically, so nonregular
    second-order distinguishability may survive; this test asserts the
    first-order statement only, which is the one the compiler depends on.
    """
    from covq.measurement import compile_schedule_readout

    from covq.instances import matching_target

    dec = wid.decompose_width2(matching_target(4, np.random.default_rng(11), load=0.8).F)
    assert dec.feasible
    assert any(min(s) < 0 for _, _, s in dec.branches), "fixture must carry a negative edge"
    prog = prg.labelled_schedule_from_width(dec)
    ps = z_generators(4)

    at_zero = compile_schedule_readout(prog, ps, np.zeros(4))
    assert at_zero["fixed_x_total_information"] == pytest.approx(0.0, abs=1e-12)
    assert at_zero["compiled_total_information"] == pytest.approx(4.0, abs=1e-9)

    uniform = compile_schedule_readout(prog, ps, np.full(4, 0.3))
    assert uniform["fixed_x_total_information"] < uniform["compiled_total_information"] - 1e-6
    assert uniform["schedule_gap"] < 1e-9

    generic = compile_schedule_readout(prog, ps, np.array([0.7, -0.35, 1.1, 0.2]))
    assert generic["schedule_gap"] < 1e-9


def test_m3_equation_110_reduces_to_the_covariance_formula_on_pure_states():
    """Eq (110) and Eq (109) must coincide exactly where the paper says they do.

    Signed cat states alone are *not* a sufficient fixture: their generator
    matrix elements are real, which hides a transpose-versus-conjugate-transpose
    error in the Eq (110) contraction.  Random complex states are included
    because that is the case which exposes it -- the buggy form was not even
    PSD there, while agreeing perfectly on every real fixture.
    """
    from covq.measurement import covariance_surrogate, mixed_state_qfim

    rng = np.random.default_rng(29)
    for m in (1, 2, 3, 4):
        ps = z_generators(m)
        for _ in range(4):
            psi = data_statevector(prg.signed_cat_circuit(rng.choice([-1, 1], size=m),
                                                          n_qubits=m))
            rho = np.outer(psi, psi.conj())
            ref = qfim_from_statevector(psi, ps)
            assert mixed_state_qfim(rho, ps) == pytest.approx(ref, abs=1e-12)
            assert covariance_surrogate(rho, ps) == pytest.approx(ref, abs=1e-12)

            v = rng.normal(size=1 << m) + 1j * rng.normal(size=1 << m)
            v /= np.linalg.norm(v)
            rho_c = np.outer(v, v.conj())
            assert mixed_state_qfim(rho_c, ps) == pytest.approx(
                qfim_from_statevector(v, ps), abs=1e-10)


def test_m3_mixed_state_qfim_is_positive_semidefinite_on_complex_mixed_states():
    """Eq (110) is a Fisher matrix; a negative eigenvalue is a bug, not noise."""
    from covq.measurement import mixed_state_qfim

    rng = np.random.default_rng(77)
    for m in (2, 3):
        ps = z_generators(m)
        dim = 1 << m
        for _ in range(15):
            a = rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))
            rho = a @ a.conj().T
            rho /= np.trace(rho).real
            assert np.linalg.eigvalsh(mixed_state_qfim(rho, ps)).min() >= -1e-10


def test_m3_covariance_surrogate_is_not_the_qfim_off_pure_states():
    """Section 11.3, quantified.

    Against the closed form for a depolarised phase probe, ``F_Q = v^2`` with
    visibility ``v = 1 - 4p/3``, while ``4 Cov`` stays pinned at one.  The error
    is unbounded, not a correction.
    """
    from covq.measurement import depolarize, noisy_report

    ps = z_generators(1)
    psi = np.array([1.0, 1.0], dtype=complex) / np.sqrt(2.0)
    rho0 = np.outer(psi, psi.conj())
    for p in (0.0, 0.05, 0.1, 0.25, 0.5):
        rep = noisy_report(depolarize(rho0, p), ps)
        v = 1.0 - 4.0 * p / 3.0
        assert rep["trace_sld_qfim"] == pytest.approx(v * v, abs=1e-10)
        assert rep["trace_covariance_surrogate"] == pytest.approx(1.0, abs=1e-10)
        assert rep["surrogate_dominates"]


def test_m3_z_correlators_are_blind_to_the_bell_primitive_s_own_dephasing():
    """The sharpest form of the Section 11.3(a) warning.

    Every generator is diagonal in ``Z``, so ``Z`` dephasing leaves the whole
    observable covariance matrix *exactly* invariant -- including the edge
    correlator the compiler was asked to hit -- while the SLD QFIM decays to
    zero.  Reporting only 11.3(a) for a Bell schedule can therefore certify a
    program that carries no information at all.
    """
    from covq.measurement import (cfi_of_readout_mixed, compile_branch_readout,
                                  covariance_surrogate, dephase, mixed_state_qfim, _rotate)

    ps = z_generators(2)
    psi = _rotate(data_statevector(prg.signed_cat_circuit(np.array([1, 1]), n_qubits=2)),
                  ps, np.array([0.4, 0.0]))
    alphas = compile_branch_readout(psi, ps, [(0, 1)])
    rho0 = np.outer(psi, psi.conj())
    cov0 = covariance_surrogate(rho0, ps)

    previous = np.inf
    for p in (0.0, 0.1, 0.2, 0.35, 0.5):
        rho = dephase(rho0, p)
        assert covariance_surrogate(rho, ps) == pytest.approx(cov0, abs=1e-12)
        fq = mixed_state_qfim(rho, ps)
        cfi = cfi_of_readout_mixed(rho, ps, alphas)
        assert np.linalg.eigvalsh(fq - cfi).min() >= -1e-9, "CFI must not exceed the QFIM"
        # Matched quadrature is optimal at every visibility, not just at v = 1:
        # there is no measurement gap to explain away here.
        assert np.trace(cfi) == pytest.approx(np.trace(fq), abs=1e-9)
        assert np.trace(fq) <= previous + 1e-12
        previous = float(np.trace(fq))
    assert previous == pytest.approx(0.0, abs=1e-10)


def test_m1_readout_is_matched_quadrature_with_unit_regularity_margin():
    """The analyzer phase is solved, not searched.

    Two parity evaluations determine the whole fringe
    ``<P_B>(alpha) = c cos alpha + d sin alpha``, so matched quadrature is the
    exact root ``atan2(-c, d)`` and the regularity margin ``|sin delta|`` is
    one by construction.  A grid argmax cannot do this: on a pure block every
    non-degenerate angle ties at Fisher information one, so the tie-break
    returns an arbitrary detuning that is free on pure states and expensive the
    moment visibility drops.
    """
    from covq.measurement import _rotate, readout_contract

    ps = z_generators(2)
    for phi in (0.4, 1.1, -0.7):
        psi = _rotate(data_statevector(prg.signed_cat_circuit(np.array([1, 1]), n_qubits=2)),
                      ps, np.array([phi, 0.0]))
        (rec,) = readout_contract(psi, [(0, 1)])
        assert rec["regularity_margin"] == pytest.approx(1.0, abs=1e-12)
        assert rec["visibility"] == pytest.approx(1.0, abs=1e-12)
        assert math.cos(phi - rec["alpha"]) == pytest.approx(0.0, abs=1e-12)


def test_m3_declared_readout_attains_the_mixed_state_qfim_under_dephasing():
    """The noisy CFI column had no gap in it; the compiler had a tie-break bug.

    Closed form for the dephased cat block at detuning ``delta`` and fringe
    visibility ``v``:

        F_C = v^2 sin^2(delta) / (1 - v^2 cos^2(delta)) * s s^T,

    maximised at quadrature where it equals ``v^2 s s^T = F_Q``.  So local
    parity readout is optimal at *every* visibility, and complete dephasing
    gives zero information for every readout rather than for this one.
    """
    from covq.measurement import (cfi_of_readout_mixed, compile_branch_readout,
                                  dephase, mixed_state_qfim, _rotate)

    for k in (2, 3, 4):
        ps = z_generators(k)
        psi = _rotate(data_statevector(prg.signed_cat_circuit(np.ones(k, int), n_qubits=k)),
                      ps, np.array([0.37] + [0.0] * (k - 1)))
        alphas = compile_branch_readout(psi, ps, [tuple(range(k))])
        rho0 = np.outer(psi, psi.conj())
        for p in (0.0, 0.1, 0.25, 0.4):
            rho = dephase(rho0, p)
            cfi = cfi_of_readout_mixed(rho, ps, alphas)
            assert np.trace(cfi) == pytest.approx(float(np.trace(mixed_state_qfim(rho, ps))),
                                                  abs=1e-9)
            v = (1.0 - 2.0 * p) ** k
            assert np.trace(cfi) == pytest.approx(k * v * v, abs=1e-9)


def test_m3_closed_form_fisher_for_a_detuned_noisy_analyzer():
    """``F_C(v, delta)`` against numerics over visibility and detuning."""
    from covq.measurement import cfi_of_readout_mixed, dephase, _rotate

    for k in (2, 3):
        ps = z_generators(k)
        for p in (0.0, 0.15, 0.3, 0.45):
            v = (1.0 - 2.0 * p) ** k
            for phi in (0.0, 0.4, 1.3):
                psi = _rotate(data_statevector(prg.signed_cat_circuit(np.ones(k, int), n_qubits=k)),
                              ps, np.array([phi] + [0.0] * (k - 1)))
                rho = dephase(np.outer(psi, psi.conj()), p)
                for a in (0.2, 0.9, math.pi / 2, 2.5):
                    alphas = np.zeros(k)
                    alphas[0] = a
                    d = phi - a
                    den = 1.0 - v * v * math.cos(d) ** 2
                    pred = v * v * math.sin(d) ** 2 / den if den > 1e-12 else 0.0
                    assert cfi_of_readout_mixed(rho, ps, alphas)[0, 0] == pytest.approx(pred, abs=1e-9)


def test_m3_generator_basis_pinching_preserves_covariance_and_kills_all_information():
    """The general proposition, not the Bell instance.

    For the pinching ``D`` onto the joint eigenbasis of the commuting family,
    ``Cov_{D(rho)}(P) = Cov_rho(P)`` for every state, while ``[D(rho), P_i] = 0``
    forces ``U_theta D(rho) U_theta^dagger = D(rho)`` and hence ``F_Q = 0``
    identically.  Generator-basis covariance can therefore be exactly unchanged
    while all parameter information is destroyed.
    """
    from covq.measurement import covariance_surrogate, mixed_state_qfim

    rng = np.random.default_rng(101)
    for m, n in ((2, 2), (3, 3), (3, 4)):
        ps = random_commuting_paulis(m, rng, n=n)
        dim = 1 << n
        basis = np.eye(dim, dtype=complex)
        mix = np.zeros((dim, dim), dtype=complex)
        for i in range(m):
            col = np.column_stack([apply_pauli(basis[:, k], ps, i) for k in range(dim)])
            mix += rng.normal() * col
        _, u = np.linalg.eigh(mix)
        for _ in range(3):
            a = rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))
            rho = a @ a.conj().T
            rho /= np.trace(rho).real
            pinched = u @ np.diag(np.diag(u.conj().T @ rho @ u)) @ u.conj().T
            assert covariance_surrogate(pinched, ps) == pytest.approx(
                covariance_surrogate(rho, ps), abs=1e-10)
            assert np.abs(mixed_state_qfim(pinched, ps)).max() < 1e-18


# ----------------------------------------------------------------------
# Eq (112): the estimator layer.
# ----------------------------------------------------------------------

def test_estimator_parity_model_reproduces_the_schedule_qfim():
    """The parity likelihood's Fisher matrix *is* the schedule QFIM.

    Not approximately: the sufficient statistic is one parity count per
    (branch, block), and its Fisher matrix is ``sum_r p_r sum_B s_B s_B^T``.
    """
    from covq.estimator import model_from_decomposition
    from covq.instances import matching_target

    for m in (3, 4, 5):
        dec = wid.decompose_width2(matching_target(m, np.random.default_rng(m * 13),
                                                   load=0.8).F)
        assert dec.feasible
        model = model_from_decomposition(dec)
        assert model.fisher() == pytest.approx(dec.matrix(), abs=1e-12)


def test_estimator_mle_is_consistent_and_saturates_the_cramer_rao_bound():
    """Attaining the CFI at a point is not an estimator; this is the estimator.

    Bias falls with the shot count and the empirical covariance matches
    ``F^+/N`` on the identifiable quotient.  The tolerance band is set by Monte
    Carlo error on sample-covariance eigenvalues at this replicate count, not by
    what the estimator happens to achieve; a tighter study is recorded in
    ``results/measurements/estimator_efficiency.json``.
    """
    from covq.estimator import efficiency_report, model_from_decomposition
    from covq.instances import matching_target

    rng = np.random.default_rng(9)
    dec = wid.decompose_width2(matching_target(3, np.random.default_rng(39), load=0.8).F)
    model = model_from_decomposition(dec)

    coarse = efficiency_report(model, np.zeros(3), 500, 600, rng)
    fine = efficiency_report(model, np.zeros(3), 20_000, 600, rng)
    assert fine["max_abs_bias"] < coarse["max_abs_bias"]
    for rep in (coarse, fine):
        assert rep["rank"] == 3
        assert rep["bias_within_4_stderr"], rep["max_abs_bias"]
        assert rep["efficiency_within_mp_band"], rep["efficiency_eigenvalues"]


# ----------------------------------------------------------------------
# K3: the QUEST baseline (arXiv:2605.02367), exact-target mode only.
# ----------------------------------------------------------------------

def test_quest_published_beats_the_greedy_routine_by_an_order_of_magnitude():
    """The fidelity gap, pinned so it cannot be forgotten again.

    Published QUEST runs two phases per iteration: insert one rotation, then
    **jointly reoptimise every accumulated angle**.  The earlier routine in this
    package omits phase two, and that omission is not cosmetic -- later
    rotations routinely make earlier angles suboptimal, so without
    reoptimisation the depth required to hit a target is a large overestimate
    and any resource comparison built on it understates the baseline.
    """
    from covq.quest import moment_constraints, quest, quest_published

    m = 4
    ps = z_generators(m)
    F = np.full((m, m), 0.4)
    np.fill_diagonal(F, 1.0)
    cons = moment_constraints(ps, F)

    greedy = quest(cons, m, max_depth=60, tol=1e-13)
    published = quest_published(cons, m, variant="tE", max_depth=60, tol=1e-13)
    assert greedy.stop_reason == "converged"
    assert published.stop_reason == "converged"
    assert published.two_qubit_rotations * 5 < greedy.two_qubit_rotations, (
        published.two_qubit_rotations, greedy.two_qubit_rotations)


def test_quest_published_variants_agree_on_the_registered_surface():
    """``bE`` searches every insertion position; ``tE`` only the terminal one."""
    from covq.instances import toeplitz_like
    from covq.quest import moment_constraints, quest_published

    m = 3
    ps = z_generators(m)
    cons = moment_constraints(ps, toeplitz_like(m, 0.35).F)
    te = quest_published(cons, m, variant="tE", max_depth=40, tol=1e-13)
    be = quest_published(cons, m, variant="bE", max_depth=40, tol=1e-13)
    assert te.stop_reason == be.stop_reason == "converged"
    assert be.two_qubit_rotations <= te.two_qubit_rotations


def test_quest_published_gradient_matches_finite_differences():
    """The adjoint gradient drives the joint phase; a wrong one would fail silently."""
    from covq.quest import PauliRotation, _cost_and_grad, _dense, moment_constraints

    rng = np.random.default_rng(3)
    n = 3
    ps = z_generators(n)
    F = np.full((n, n), 0.3)
    np.fill_diagonal(F, 1.0)
    obs = [(np.asarray(O), float(t)) for O, t in moment_constraints(ps, F)]
    rots = [(None, _dense(PauliRotation(((0, "X"), (1, "Y"))).axes, n)),
            (None, _dense(PauliRotation(((2, "Z"),)).axes, n)),
            (None, _dense(PauliRotation(((1, "X"), (2, "X"))).axes, n))]
    psi0 = np.ones(1 << n, dtype=complex) / math.sqrt(1 << n)
    angles = list(rng.uniform(-2, 2, len(rots)))
    _, grad = _cost_and_grad(angles, rots, psi0, obs)
    eps = 1e-6
    for k in range(len(angles)):
        up, dn = list(angles), list(angles)
        up[k] += eps
        dn[k] -= eps
        fd = (_cost_and_grad(up, rots, psi0, obs)[0]
              - _cost_and_grad(dn, rots, psi0, obs)[0]) / (2 * eps)
        assert grad[k] == pytest.approx(fd, abs=1e-6), (k, grad[k], fd)


def test_quest_hits_exact_first_and_second_moment_targets():
    """The baseline must actually work before any comparison means anything."""
    from covq.instances import matching_target, toeplitz_like
    from covq.quest import moment_constraints, quest_published

    for name, F, m in (("matching3", matching_target(3, np.random.default_rng(39), load=0.8).F, 3),
                       ("toeplitz3", toeplitz_like(3, 0.35).F, 3),
                       ("toeplitz4", toeplitz_like(4, 0.35).F, 4)):
        ps = z_generators(m)
        res = quest_published(moment_constraints(ps, F), m, variant="tE",
                              max_depth=48, tol=1e-13)
        assert res.stop_reason == "converged", (name, res.stop_reason, res.residual)
        psi = res.state / np.linalg.norm(res.state)
        assert qfim_from_statevector(psi, ps) == pytest.approx(F, abs=1e-5)


def test_quest_start_state_matters_and_the_default_is_the_fair_one():
    """``|0..0>`` is a joint ``Z`` eigenstate and a trap for this constraint class.

    Every ``<Z_i>`` is maximally wrong there and every ``Z``-type pool element is
    a no-op, so greedy descent stalls on targets it clears easily from
    ``|+>^n``.  Reporting the ``|0..0>`` numbers as the method's performance
    would understate the baseline; this test pins the difference so the default
    cannot regress silently.
    """
    from covq.instances import toeplitz_like
    from covq.quest import moment_constraints, quest

    m = 3
    F = toeplitz_like(m, 0.35).F
    cons = moment_constraints(z_generators(m), F)

    zero = np.zeros(1 << m, dtype=complex)
    zero[0] = 1.0
    trapped = quest(cons, m, max_depth=48, tol=1e-13, psi0=zero)
    default = quest(cons, m, max_depth=48, tol=1e-13)

    assert trapped.stop_reason != "converged"
    assert default.stop_reason == "converged"
    assert default.depth_adaptive_length < trapped.depth_adaptive_length


def test_quest_resources_are_parsed_from_emitted_circuits():
    """Gate C9: no analytic resource estimates in the comparison."""
    from covq.instances import toeplitz_like
    from covq.quest import moment_constraints, quest, quest_circuit

    m = 4
    res = quest(moment_constraints(z_generators(m), toeplitz_like(m, 0.35).F), m,
                max_depth=48, tol=1e-13)
    parsed = resources(lower_to_cx(quest_circuit(res, m)))
    assert parsed["cx_count"] == 2 * res.two_qubit_rotations
    assert parsed["cx_count"] > 0


def test_quest_line_search_is_exact_three_point_trigonometry():
    """``<R(t) O R(t)^dag> = a + b cos t + c sin t`` exactly, for Pauli ``P``."""
    from covq.quest import PauliRotation, _apply_rotation, _dense

    rng = np.random.default_rng(5)
    n = 3
    dense = _dense(PauliRotation(((0, "X"), (2, "Y"))).axes, n)
    obs = _dense(PauliRotation(((1, "Z"), (2, "Z"))).axes, n)
    psi = rng.normal(size=1 << n) + 1j * rng.normal(size=1 << n)
    psi /= np.linalg.norm(psi)

    def val(t):
        v = _apply_rotation(psi, dense, t)
        return float(np.real(np.vdot(v, obs @ v)))

    v0, vp, vm = val(0.0), val(math.pi / 2), val(-math.pi / 2)
    a = 0.5 * (vp + vm)
    b, c = v0 - a, 0.5 * (vp - vm)
    for t in (0.3, 1.7, -2.2, 3.0):
        assert val(t) == pytest.approx(a + b * math.cos(t) + c * math.sin(t), abs=1e-12)


def test_adaptive_recentering_needs_no_oracle_operating_point():
    """Answers the obvious objection to M2 without assuming the answer.

    The matched readout must be phased to the operating point, which is the
    thing being estimated.  A two-stage protocol resolves it: a pilot split
    between ``A = 0`` and ``A = pi/2`` determines each block phase including its
    sign (one setting cannot -- ``cos`` is even), then the main stage runs at
    unit regularity margin.  No pilot data is discarded, so the overhead is far
    below the naive ``1/(1 - pilot_fraction)``.
    """
    from covq.estimator import adaptive_efficiency_report, model_from_decomposition
    from covq.instances import matching_target

    rng = np.random.default_rng(3)
    dec = wid.decompose_width2(matching_target(3, np.random.default_rng(39), load=0.8).F)
    model = model_from_decomposition(dec)
    theta = np.array([0.6, -0.4, 0.9])

    for fraction in (0.05, 0.2):
        rep = adaptive_efficiency_report(model, theta, 20_000, 400, rng,
                                         pilot_fraction=fraction)
        assert rep["median_regularity_margin"] > 0.9
        assert rep["max_abs_bias"] <= 5.0 * rep["bias_standard_error"]
        # Well inside the discard-the-pilot penalty, because the pilot is kept.
        assert rep["worst_efficiency_ratio"] < 1.0 / (1.0 - fraction) + 0.35


def test_adaptive_pilot_seed_is_what_makes_the_likelihood_tractable():
    """The periodic likelihood is multimodal; the pilot supplies the seed.

    Least squares on the design matrix against the pilot block phases is a
    consistent starting estimate, and it is the only reason the final MLE lands
    in the right mode without oracle knowledge of ``theta``.
    """
    from covq.estimator import adaptive_recentering, model_from_decomposition
    from covq.instances import matching_target

    rng = np.random.default_rng(11)
    dec = wid.decompose_width2(matching_target(3, np.random.default_rng(39), load=0.8).F)
    model = model_from_decomposition(dec)
    theta = np.array([0.6, -0.4, 0.9])
    run = adaptive_recentering(model, theta, 40_000, rng, pilot_fraction=0.1)
    assert np.abs(run["pilot_seed"] - theta).max() < 0.1
    assert np.abs(run["estimate"] - theta).max() < 0.02
    assert run["achieved_regularity_margin"] > 0.9


# ----------------------------------------------------------------------
# Eq (111): noise gates N1-N6, N8.
# ----------------------------------------------------------------------

def _all_edges(m):
    return [(i, j) for i in range(m) for j in range(i + 1, m)]


def test_n1_noiseless_limit_reproduces_the_ideal_signed_matching_template():
    from covq.noise import BlockLocalNoise, branch_template

    m, theta = 4, np.array([0.3, -0.2, 0.5, 0.1])
    matching, signs = [(0, 1), (2, 3)], [1, -1]
    ideal = np.zeros((m, m))
    for (i, j), s in zip(matching, signs):
        ideal[i, i] = ideal[j, j] = 1.0
        ideal[i, j] = ideal[j, i] = s
    got = branch_template(matching, signs, m, theta, BlockLocalNoise())
    assert got == pytest.approx(ideal, abs=1e-12)


def test_n2_dephased_signed_cat_attains_the_analytic_quadrature_value():
    """``F_C = F_Q = v^2 s s^T`` at quadrature, for the noisy block template."""
    from covq.noise import BlockLocalNoise, pair_template

    for p in (0.0, 0.05, 0.15, 0.3):
        noise = BlockLocalNoise(dephasing={0: p, 1: p})
        block = pair_template((0.4, 0.0), (0, 1), 1, noise)
        v = (1.0 - 2.0 * p) ** 2
        assert block == pytest.approx(v * v * np.ones((2, 2)), abs=1e-9)


def test_n3_complete_pinching_drives_both_information_measures_to_zero():
    """Covariance survives untouched; the QFIM and the declared CFI do not."""
    from covq.measurement import covariance_surrogate, mixed_state_qfim
    from covq.noise import BlockLocalNoise, pair_channel, pair_template

    clean = pair_channel((0.4, 0.0), (0, 1), 1, BlockLocalNoise())
    dead = pair_channel((0.4, 0.0), (0, 1), 1, BlockLocalNoise(dephasing={0: 0.5, 1: 0.5}))
    ps = z_generators(2)
    assert covariance_surrogate(dead, ps) == pytest.approx(covariance_surrogate(clean, ps),
                                                           abs=1e-10)
    assert np.abs(mixed_state_qfim(dead, ps)).max() < 1e-12
    assert np.abs(pair_template((0.4, 0.0), (0, 1), 1,
                                BlockLocalNoise(dephasing={0: 0.5, 1: 0.5}))).max() < 1e-12


def test_n4_noisy_branch_cfi_is_exactly_edge_additive():
    """``F_{C,M,sigma} = F_0 + sum_e DeltaF_{e,sigma_e}`` against the idling reference."""
    from covq.noise import BlockLocalNoise, branch_template, delta_edge, product_template

    rng = np.random.default_rng(31)
    for m in (3, 4, 5):
        edges = _all_edges(m)
        noise = BlockLocalNoise(
            dephasing={q: float(rng.uniform(0, 0.1)) for q in range(m)},
            depolarizing={q: float(rng.uniform(0, 0.05)) for q in range(m)},
            edge_depolarizing={e: float(rng.uniform(0, 0.12)) for e in edges},
            idle_dephasing=float(rng.uniform(0, 0.05)))
        theta = rng.uniform(-1, 1, m)
        matching = [(0, 1)] if m == 3 else [(0, 1), (2, 3)]
        for signs in itertools.product((1, -1), repeat=len(matching)):
            got = branch_template(matching, list(signs), m, theta, noise)
            want = product_template(m, theta, noise, idle=True) + sum(
                delta_edge(e, s, m, theta, noise) for e, s in zip(matching, signs))
            assert got == pytest.approx(want, abs=1e-12)


def test_n5_block_local_noise_preserves_matching_based_pricing():
    """The theorem: noise reweights edges, it does not change the combinatorics.

    One refinement the ideal case does not need -- the pairless branch is a
    *separate* column, because it activates no pair and therefore never waits
    through a pair preparation, so it does not carry idle dephasing.  Dropping
    it loses the optimum whenever idling costs more than the best edge gains.
    """
    from covq.noise import BlockLocalNoise, price_branch_bruteforce, price_branch_noisy

    rng = np.random.default_rng(7)
    checked = 0
    for m in (3, 4, 5):
        edges = _all_edges(m)
        for _ in range(4):
            theta = rng.uniform(-1, 1, m)
            noise = BlockLocalNoise(
                dephasing={q: float(rng.uniform(0, 0.12)) for q in range(m)},
                depolarizing={q: float(rng.uniform(0, 0.05)) for q in range(m)},
                edge_depolarizing={e: float(rng.uniform(0, 0.15)) for e in edges},
                idle_dephasing=float(rng.uniform(0, 0.06)))
            a = rng.standard_normal((m, m))
            Q = a @ a.T / m
            costs = {e: float(rng.uniform(0.0, 0.35)) for e in edges}
            oracle = price_branch_noisy(Q, m, edges, theta, noise, costs)
            brute = price_branch_bruteforce(Q, m, edges, theta, noise, costs)
            assert oracle["value"] == pytest.approx(brute["value"], abs=1e-9)
            checked += 1
    assert checked == 12


def test_n6_noisy_floor_compiler_closes_its_duality_gap():
    """Eq (111) stays convex under block-local noise, so the dual bound is tight."""
    from covq.noise import BlockLocalNoise, noise_aware_floor_compile

    m = 4
    edges = _all_edges(m)
    G = np.full((m, m), 0.6)
    np.fill_diagonal(G, 1.2)
    for noise in (BlockLocalNoise(),
                  BlockLocalNoise(dephasing={q: 0.02 for q in range(m)},
                                  edge_depolarizing={e: 0.02 for e in edges},
                                  idle_dephasing=0.01)):
        res = noise_aware_floor_compile(G, m, edges, np.zeros(m), noise,
                                        costs={e: 0.1 for e in edges}, c0=1.0)
        assert res["status"] == "solved"
        assert res["relative_gap"] < 1e-9
        assert res["floor_slack_min_eig"] > -1e-8


def test_n6_noise_switches_entanglement_off_above_a_threshold():
    """The optimizer does not merely return the ideal schedule everywhere.

    Past an edge-noise threshold the pair primitive stops paying for itself and
    the compiler abandons it, after which the cost is independent of edge noise
    because no edge is used.  A noise-aware compiler that never changed its
    answer would be evidence that the objective was not doing any work.
    """
    from covq.noise import BlockLocalNoise, noise_aware_floor_compile

    m = 4
    edges = _all_edges(m)
    G = np.full((m, m), 0.6)
    np.fill_diagonal(G, 1.2)

    def run(q_edge):
        noise = BlockLocalNoise(dephasing={q: 0.02 for q in range(m)},
                                edge_depolarizing={e: q_edge for e in edges},
                                idle_dephasing=0.01)
        return noise_aware_floor_compile(G, m, edges, np.zeros(m), noise,
                                         costs={e: 0.1 for e in edges}, c0=1.0)

    low, high, higher = run(0.05), run(0.20), run(0.30)
    assert low["n_entangled_settings_used"] > 0
    assert high["n_entangled_settings_used"] == 0
    assert high["cost"] == pytest.approx(higher["cost"], abs=1e-9)
    assert low["cost"] < high["cost"]


def test_n8_asymmetric_readout_needs_a_third_analyzer_point():
    """Two points determine the fringe only when it has no offset.

    Under independent outcome-dependent confusion with rates ``(e0, e1)`` per
    qubit, the reported outcome satisfies ``E[x~|x] = a + b x`` with
    ``a = e1 - e0`` and ``b = 1 - e0 - e1``.  For a zero-mean block the parity
    expectation therefore becomes

        <P~_B> = (prod_i a_i) + (prod_i b_i) <P_B>,

    an *offset* first-harmonic fringe.  A two-point fit at ``0`` and ``pi/2``
    silently attributes the offset to the harmonic and returns the wrong
    quadrature angle; the three-point fit separates them.
    """
    from covq.noise import BlockLocalNoise

    noise = BlockLocalNoise(readout_confusion={0: (0.10, 0.02), 1: (0.08, 0.01)})
    assert not noise.is_symmetric_readout
    a = [e1 - e0 for e0, e1 in (noise.confusion(0), noise.confusion(1))]
    b = [1.0 - e0 - e1 for e0, e1 in (noise.confusion(0), noise.confusion(1))]
    offset, gain = a[0] * a[1], b[0] * b[1]
    assert abs(offset) > 1e-3, "fixture must actually be asymmetric"

    phi = 0.4

    def observed(alpha):
        return offset + gain * math.cos(phi - alpha)

    # Three-point fit recovers offset, harmonic, and hence the true quadrature.
    angles = (0.0, 2 * math.pi / 3, 4 * math.pi / 3)
    v = [observed(t) for t in angles]
    b_fit = sum(v) / 3.0
    c_fit = (2.0 / 3.0) * sum(vi * math.cos(t) for vi, t in zip(v, angles))
    d_fit = (2.0 / 3.0) * sum(vi * math.sin(t) for vi, t in zip(v, angles))
    assert b_fit == pytest.approx(offset, abs=1e-12)
    alpha3 = math.atan2(-c_fit, d_fit)
    assert observed(alpha3) == pytest.approx(offset, abs=1e-12)   # harmonic nulled

    # Two-point fit mistakes the offset for signal and lands off quadrature.
    c2, d2 = observed(0.0), observed(math.pi / 2)
    alpha2 = math.atan2(-c2, d2)
    assert abs(math.cos(phi - alpha2)) > 1e-3
    assert abs(alpha2 - alpha3) > 1e-3


def test_zero_contrast_is_reported_not_papered_over():
    """``atan2(-c, d)`` is undefined at zero contrast; a library zero is not an answer."""
    from covq.measurement import block_quadrature
    from covq.noise import BlockLocalNoise, pair_channel

    dead = pair_channel((0.4, 0.0), (0, 1), 1, BlockLocalNoise(dephasing={0: 0.5, 1: 0.5}))
    rec = block_quadrature(dead, [0, 1], 2)
    assert rec["analyzer_status"] == "arbitrary_zero_information"
    assert rec["regularity_margin"] is None
    assert rec["visibility"] == pytest.approx(0.0, abs=1e-12)


def test_n7_deployable_pilot_arm_is_bounded_by_the_oracle_angle():
    """Oracle analyzer angles are an upper bound; the deployable arm pays for the pilot.

    Two distinct losses, and they pull in opposite directions.  Pilot-angle
    estimation error falls as the pilot grows; pilot shots run off quadrature,
    which costs in proportion to the pilot fraction.  At unit visibility both
    vanish, because every non-degenerate angle already attains.
    """
    from covq.noise import BlockLocalNoise, pair_channel, pilot_recentered_block

    rng = np.random.default_rng(5)
    total = 20_000

    clean = pair_channel((0.4, 0.0), (0, 1), 1, BlockLocalNoise())
    for fraction in (0.01, 0.2):
        n_pilot = int(fraction * total)
        rep = pilot_recentered_block(clean, (0, 1), BlockLocalNoise(), n_pilot,
                                     total - n_pilot, rng, n_reps=60)
        assert rep["ratio"] == pytest.approx(1.0, abs=1e-9)

    noise = BlockLocalNoise(dephasing={0: 0.12, 1: 0.12}, edge_depolarizing={(0, 1): 0.15})
    rho = pair_channel((0.4, 0.0), (0, 1), 1, noise)
    tiny = pilot_recentered_block(rho, (0, 1), noise, 40, total - 40, rng, n_reps=200)
    sized = pilot_recentered_block(rho, (0, 1), noise, 400, total - 400, rng, n_reps=200)
    for rep in (tiny, sized):
        assert rep["ratio"] <= 1.0 + 1e-9, "deployable can never beat the oracle"
    assert sized["ratio"] > tiny["ratio"], "a too-small pilot is the dominant loss"
    assert sized["ratio"] > 0.9


def test_n9_greedy_support_selection_is_not_exact():
    """Recorded as a failure, because it is one.

    Forward selection commits to the best single column, which need not belong
    to the best pair.  It is exact only once the cardinality penalty is large
    enough that one setting is genuinely optimal -- that is, outside the
    multi-setting regime the schedule exists to exploit.
    """
    from covq.noise import BlockLocalNoise, fixed_setting_cost_compile

    m = 3
    edges = [(0, 1), (0, 2), (1, 2)]
    G = np.full((m, m), 0.5)
    np.fill_diagonal(G, 1.1)
    noise = BlockLocalNoise(dephasing={q: 0.02 for q in range(m)},
                            edge_depolarizing={e: 0.04 for e in edges},
                            idle_dephasing=0.01)
    kw = dict(costs={e: 0.1 for e in edges})

    cheap_g = fixed_setting_cost_compile(G, m, edges, np.zeros(m), noise, 0.0,
                                         method="greedy", **kw)
    cheap_e = fixed_setting_cost_compile(G, m, edges, np.zeros(m), noise, 0.0,
                                         method="exhaustive", max_support=2, **kw)
    assert cheap_e["total_cost"] < cheap_g["total_cost"] - 1e-9
    assert cheap_e["n_settings"] > cheap_g["n_settings"]

    dear_g = fixed_setting_cost_compile(G, m, edges, np.zeros(m), noise, 1.0,
                                        method="greedy", **kw)
    dear_e = fixed_setting_cost_compile(G, m, edges, np.zeros(m), noise, 1.0,
                                        method="exhaustive", max_support=2, **kw)
    assert dear_g["total_cost"] == pytest.approx(dear_e["total_cost"], abs=1e-7)


def test_j2_frozen_pilot_policy_is_never_better_than_the_oracle_ceiling():
    """The deployable N10 column must cost more than the oracle-angle diagnostic.

    A deployable number that beat the oracle would mean the pilot was being
    given information it has not paid for.  The registered policy is
    ``f = 0.02`` at phases ``{0, pi/2}``, pilot retained in the likelihood,
    globally fixed -- no per-instance tuning.
    """
    from covq.noise import (FROZEN_PILOT_POLICY, BlockLocalNoise, deployable_exposure,
                            noise_aware_floor_compile)

    assert FROZEN_PILOT_POLICY["pilot_fraction"] == 0.02
    assert FROZEN_PILOT_POLICY["selection_type"] == "globally_fixed"
    assert FROZEN_PILOT_POLICY["pilot_in_final_likelihood"] is True

    m = 4
    edges = [(i, j) for i in range(m) for j in range(i + 1, m)]
    G = np.full((m, m), 0.6)
    np.fill_diagonal(G, 1.2)
    for q_edge in (0.0, 0.05, 0.20):
        noise = BlockLocalNoise(dephasing={q: 0.02 for q in range(m)},
                                edge_depolarizing={e: q_edge for e in edges},
                                idle_dephasing=0.01)
        oracle = noise_aware_floor_compile(G, m, edges, np.zeros(m), noise,
                                           costs={e: 0.0 for e in edges}, c0=1.0)
        deploy = deployable_exposure(G, m, edges, np.zeros(m), noise, oracle["branches"],
                                     costs={e: 0.0 for e in edges}, c0=1.0)
        assert deploy["status"] == "solved"
        assert deploy["cost"] >= oracle["cost"] - 1e-9
        assert deploy["cost"] / oracle["cost"] < 1.10
        assert deploy["floor_slack_min_eig"] > -1e-6
