"""Tests for the claims the audit actually relies on.

Each test corresponds to a statement in
``docs/audits/TASK0AB_NOVELTY_THEOREM_AUDIT_v0.2.md``.  They are written so that a
*false theorem* fails them, not just a broken implementation.
"""

from __future__ import annotations

import itertools

import numpy as np
import pytest

from covq import baselines as bl
from covq import gates as gt
from covq import polytope as pol
from covq import programs as prg
from covq import width as wid
from covq.circuits import is_hardware_legal, lower_to_cx, route
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
