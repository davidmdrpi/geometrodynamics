"""Is the v4 CKM realization a prediction, or a locally flexible fit?

THE ANSWER
──────────
**rank K = 4. The CKM realization is a fit, not a prediction.**

The decisive object is not another pseudoinverse. Build the two response maps
over one common parameter chart ``x``,

    J_M = ∂y_M/∂x ,   y_M = ln(m_{s,c,b,t} / m_d)
    J_F = ∂y_F/∂x ,   y_F = (θ₁₂, θ₂₃, θ₁₃, δ)

take the mass-preserving tangent space ``N_M = ker J_M``, and form

    K = J_F N_M .

``rank K`` counts the physically independent CKM directions reachable **without
disturbing the masses**. If it equals 4 — the full dimension of the physical
flavor space of a unitary 3×3 matrix — then the mass-preserving parameter
freedom already spans everything the CKM can be, and fitting it predicts
nothing.

It equals 4, robustly: stable across finite-difference steps ``1e-5…1e-7`` and
null-space cutoffs ``1e-6…1e-10``, with the four singular values spread over
only ``379×`` and no near-degeneracy. Confirmed by direct construction — an
arbitrary target ``δy_F`` is realised to ``1e-14`` while the masses hold to
``1e-14``.

There is therefore **no left-null vector**, hence **no first-order invariant
relation** ``wᵀδy_F = 0`` to compare against experiment. The round was set up
to extract one if it existed. It does not.

The φ_h A/B test
────────────────
The library treats the Hopf holonomy ``φ_h = π/k₅`` as *derived* and as the
source of CP structure. If holding it fixed dropped the flavor rank from 4 to
3, and the missing direction were the observed CP relation, that would be
evidence that topology removes one calibration freedom.

    φ_h RELEASED   rank K = 4   σ = [2.611, 0.468, 0.00827, 0.00689]
    φ_h FIXED      rank K = 4   σ = [0.542, 0.366, 0.00793, 0.00357]

**Fixing it does not lower the rank.** The other fitted matrix elements absorb
arbitrary CKM data on their own, so the derived phase is not, by itself,
producing a flavor prediction. It is the single most *efficient* CP handle —
releasing it multiplies the leading singular value by ``4.8`` — but efficiency
is not identifiability.

This confirms and sharpens PR #173, which found that *adding* ``φ_h`` as an
input left the observable rank unchanged. That was a statement about spanning;
this is the stronger one: with ``φ_h`` held at its derived value, the
mass-preserving freedom still covers all four physical flavor coordinates.

The calibration-DOF census, measured rather than counted
────────────────────────────────────────────────────────
"New symbol count" and "calibration degree-of-freedom count" are different
numbers. Measuring the second as ``rank J_F`` restricted to each knob group:

    group                          symbols   measured flavor rank
    v4 targeted etas                  3              3
    eta_k3k5_minus retune             1              1
    diag_shift_plus                   3              2
    diag_shift_minus                  3              2
    all diagonal shifts               6              3
    phi_h                             1              1
    v4 additions, phi_h fixed         9              4
    v4 additions, phi_h released     10              4

Two structural facts fall out. **The trace direction of each diagonal-shift
triple is an exact CKM gauge** — ``|J_F·1| ≈ 2e-10`` against ``|J_M·1| = 12.5``
— so a uniform shift within a block moves masses and cannot have been selected
using flavor data. And the realised values are **traceless to ~1e-10** in both
blocks, with ``diag_shift_plus`` further collapsing to ``(+d, −d, 0)``.

So the nine v4 additions with ``φ_h`` fixed supply **4** independent flavor
directions — exactly the dimension of the physical flavor space.

The "+3 parameters for +5 independent observables" claim
────────────────────────────────────────────────────────
It cannot hold, for a reason that needs no computation: **a unitary 3×3 CKM has
exactly four physical parameters.** The nine quoted flavor-CP observables
(``|V_us|, |V_cb|, |V_ub|, |V_td|, |V_ts|, J, β, γ, sin δ``) are all functions
of those four, so at most four of them are independent. "+5 independent
observables" exceeds the ceiling.

Against that ceiling the measured calibration dimension is 4. **Net predictive
surplus ≤ 0**, not ``+2``.

What this does not say
──────────────────────
Nothing here says the v4 numbers are wrong, or that the CKM agreement is
accidental in the sense of being lucky. The lock reproduces nine observables at
``≤ 1%`` and that is a real property of the realisation. What the rank shows is
that the agreement is **not evidence for the Hamiltonian**, because the same
Hamiltonian could have reproduced any *locally neighbouring* CKM equally well at
the same masses.

Two scope limits worth stating. This is a **local** statement at the v4 lock: a
rank is a first-order object and says nothing about how far the mass-preserving
surface extends before nonlinearity or positivity bites. And the ``|δx|``
required is not small — an arbitrary unit ``δy_F`` needs ``|δx| ≈ 50–80`` in
these coordinates — so "reachable" is a statement about directions, not about
comfortable distances.
"""

from __future__ import annotations

import math
from dataclasses import replace
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.linalg import null_space

from geometrodynamics.qcd.quark_spectrum import (LOCKED_QUARK_PARAMS,
                                                 LOCKED_QUARK_PARAMS_V4,
                                                 extract_ckm_matrix,
                                                 extract_physical_spectrum)

__all__ = [
    "MASS_OBSERVABLES",
    "FLAVOR_OBSERVABLES",
    "SCALAR_COORDS",
    "TUPLE_COORDS",
    "COORD_NAMES",
    "V3_KNOBS",
    "V4_TARGETED_ETAS",
    "DIAG_SHIFT_PLUS",
    "DIAG_SHIFT_MINUS",
    "PHYSICAL_FLAVOR_DIMENSION",
    "mass_observables",
    "flavor_observables",
    "jacobian",
    "mass_preserving_tangent_space",
    "flavor_response_on_mass_preserving_directions",
    "measure_the_observables_are_four_and_rephasing_invariant",
    "measure_both_jacobians_converge",
    "measure_the_calibration_dof_census",
    "measure_the_mass_preserving_flavor_rank",
    "measure_the_reachability_is_nonlinearly_true",
    "measure_the_restricted_v4_calibration_rank",
    "measure_the_phi_h_ab_test",
    "measure_the_counting_claim",
    "measure_the_flavor_ledger",
]

#: The four independently scored mass ratios (`u` is zero by construction,
#: `d` is the anchor).
MASS_OBSERVABLES: Tuple[str, ...] = ("s", "c", "b", "t")

#: Four genuinely independent flavor coordinates — NOT ten redundant
#: moduli/angles/invariants. `J`, `α`, `β`, `γ`, `|V_td|`, `|V_ts|` are
#: functions of these and are validation outputs, not extra rows.
FLAVOR_OBSERVABLES: Tuple[str, ...] = ("theta12", "theta23", "theta13", "delta")

#: A unitary 3×3 matrix modulo rephasing has exactly this many real parameters.
PHYSICAL_FLAVOR_DIMENSION: int = 4

V3_KNOBS: Tuple[str, ...] = ("beta", "gamma_q", "transport", "pinhole",
                             "resistance", "uplift_asymmetry", "chi_q_k3",
                             "eta_k3k5_minus")
V4_TARGETED_ETAS: Tuple[str, ...] = ("eta_k1k3_plus", "eta_k1k3_minus",
                                     "eta_k1k5_minus")
SCALAR_COORDS: Tuple[str, ...] = V3_KNOBS + V4_TARGETED_ETAS + ("phi_h",)
TUPLE_COORDS: Tuple[Tuple[str, int], ...] = tuple(
    [("diag_shift_plus", i) for i in range(3)]
    + [("diag_shift_minus", i) for i in range(3)])
COORD_NAMES: Tuple[str, ...] = tuple(
    list(SCALAR_COORDS) + [f"{n}[{i}]" for n, i in TUPLE_COORDS])
DIAG_SHIFT_PLUS: Tuple[str, ...] = tuple(f"diag_shift_plus[{i}]"
                                         for i in range(3))
DIAG_SHIFT_MINUS: Tuple[str, ...] = tuple(f"diag_shift_minus[{i}]"
                                          for i in range(3))

_KEYS: List = list(SCALAR_COORDS) + list(TUPLE_COORDS)
_DEFAULT_STEP = 1e-6


# ── the two observable maps ─────────────────────────────────────────────────

def _bump(params, key, h: float):
    if isinstance(key, tuple):
        name, index = key
        values = list(getattr(params, name))
        values[index] += h
        return replace(params, **{name: tuple(values)})
    return replace(params, **{key: getattr(params, key) + h})


def mass_observables(params=None) -> np.ndarray:
    """``ln(m_X / m_d)`` for the four scored species."""
    spectrum = extract_physical_spectrum(
        params if params is not None else LOCKED_QUARK_PARAMS_V4)
    return np.array([math.log(spectrum[s] / spectrum["d"])
                     for s in MASS_OBSERVABLES])


def flavor_observables(params=None) -> np.ndarray:
    """``(θ₁₂, θ₂₃, θ₁₃, δ)`` in the PDG convention.

    Built only from ``|V_ij|`` and the Jarlskog invariant, both invariant under
    the arbitrary per-eigenvector phases ``eigh`` returns. Using the raw complex
    entries would make the Jacobian a function of LAPACK's phase convention.
    """
    V = extract_ckm_matrix(params if params is not None
                           else LOCKED_QUARK_PARAMS_V4)
    mod = np.abs(V)
    s13 = float(mod[0, 2])
    c13 = math.sqrt(max(0.0, 1.0 - s13 * s13))
    theta12 = math.asin(min(1.0, float(mod[0, 1]) / c13))
    theta23 = math.asin(min(1.0, float(mod[1, 2]) / c13))
    theta13 = math.asin(min(1.0, s13))
    jarlskog = float(np.imag(V[0, 1] * V[1, 2]
                             * np.conj(V[0, 2]) * np.conj(V[1, 1])))
    denominator = (math.cos(theta12) * math.sin(theta12)
                   * math.cos(theta23) * math.sin(theta23)
                   * c13 * c13 * s13)
    sin_delta = jarlskog / denominator if abs(denominator) > 1e-300 else 0.0
    delta = math.asin(max(-1.0, min(1.0, sin_delta)))
    return np.array([theta12, theta23, theta13, delta])


def jarlskog(params=None) -> float:
    """The rephasing invariant, kept as a *validation output* not an observable row."""
    V = extract_ckm_matrix(params if params is not None
                           else LOCKED_QUARK_PARAMS_V4)
    return float(np.imag(V[0, 1] * V[1, 2]
                         * np.conj(V[0, 2]) * np.conj(V[1, 1])))


# ── the two Jacobians, on one common chart ──────────────────────────────────

def jacobian(which: str, step: float = _DEFAULT_STEP,
             coords: Optional[Sequence[str]] = None,
             params=None) -> np.ndarray:
    """``∂y/∂x`` by central difference, in **linear** (absolute) coordinates.

    Linear rather than log because several v4 knobs sit at or near zero, and
    because the round's conclusion is a **rank** — invariant under any
    invertible reparameterisation, so no chart choice is being smuggled in.
    """
    base = params if params is not None else LOCKED_QUARK_PARAMS_V4
    fn = mass_observables if which == "mass" else flavor_observables
    names = list(COORD_NAMES) if coords is None else list(coords)
    columns = []
    for name in names:
        key = _KEYS[COORD_NAMES.index(name)]
        columns.append((fn(_bump(base, key, step))
                        - fn(_bump(base, key, -step))) / (2.0 * step))
    return np.column_stack(columns)


def mass_preserving_tangent_space(step: float = _DEFAULT_STEP,
                                  coords: Optional[Sequence[str]] = None,
                                  rcond: float = 1e-8) -> np.ndarray:
    """``N_M = ker J_M`` — the directions that move no mass at first order."""
    return null_space(jacobian("mass", step, coords), rcond=rcond)


def flavor_response_on_mass_preserving_directions(
        step: float = _DEFAULT_STEP, coords: Optional[Sequence[str]] = None,
        rcond: float = 1e-8) -> Tuple[np.ndarray, np.ndarray]:
    """``K = J_F N_M``, and its singular values."""
    N = mass_preserving_tangent_space(step, coords, rcond)
    if N.shape[1] == 0:
        empty = np.zeros((len(FLAVOR_OBSERVABLES), 0))
        return empty, np.zeros(0)
    K = jacobian("flavor", step, coords) @ N
    return K, np.linalg.svd(K, compute_uv=False)


# ── the measurements ────────────────────────────────────────────────────────

def measure_the_observables_are_four_and_rephasing_invariant() -> dict:
    """F0 — four independent coordinates, not ten redundant ones.

    The ceiling matters more than the numbers: a unitary ``3×3`` matrix modulo
    rephasing carries exactly four real parameters, so no choice of flavor
    observables can supply more than four independent rows.
    """
    V = extract_ckm_matrix(LOCKED_QUARK_PARAMS_V4)
    y = flavor_observables()
    unitarity = float(np.max(np.abs(V.conj().T @ V - np.eye(3))))

    rng = np.random.default_rng(0)
    drift = 0.0
    for _ in range(5):
        left = np.diag(np.exp(1j * rng.uniform(0, 2 * math.pi, 3)))
        right = np.diag(np.exp(1j * rng.uniform(0, 2 * math.pi, 3)))
        W = left @ V @ right
        j_w = float(np.imag(W[0, 1] * W[1, 2]
                            * np.conj(W[0, 2]) * np.conj(W[1, 1])))
        drift = max(drift, abs(j_w - jarlskog()),
                    float(np.max(np.abs(np.abs(W) - np.abs(V)))))

    return {
        "observables": list(FLAVOR_OBSERVABLES),
        "values_rad": y.tolist(),
        "values_deg": [math.degrees(v) for v in y],
        "jarlskog": jarlskog(),
        "ckm_unitarity_defect": unitarity,
        "max_drift_under_random_rephasing": drift,
        "physical_flavor_dimension": PHYSICAL_FLAVOR_DIMENSION,
        "validation_outputs_not_rows": (
            "J, alpha, beta, gamma, |V_td|, |V_ts| and the remaining moduli are "
            "functions of these four and are NOT independent observable rows"),
        "why_this_matters": (
            "A unitary 3x3 CKM modulo rephasing has exactly four real "
            "parameters, so no observable set can supply more than four "
            "independent flavor rows"),
    }


def measure_both_jacobians_converge(
        steps: Sequence[float] = (1e-4, 1e-5, 1e-6, 1e-7)) -> dict:
    """F1 — both maps run through eigen-decompositions; check they differentiate."""
    rows = []
    for step in steps:
        JM = jacobian("mass", step)
        JF = jacobian("flavor", step)
        rows.append({"step": step,
                     "mass_frobenius": float(np.linalg.norm(JM)),
                     "flavor_frobenius": float(np.linalg.norm(JF)),
                     "rank_mass": int(np.linalg.matrix_rank(JM, tol=1e-8)),
                     "rank_flavor": int(np.linalg.matrix_rank(JF, tol=1e-8))})
    m = [r["mass_frobenius"] for r in rows]
    f = [r["flavor_frobenius"] for r in rows]
    return {
        "rows": rows,
        "mass_relative_spread": (max(m) - min(m)) / max(m),
        "flavor_relative_spread": (max(f) - min(f)) / max(f),
        "both_converged": ((max(m) - min(m)) / max(m) < 1e-5
                           and (max(f) - min(f)) / max(f) < 1e-5),
        "ranks_stable": len({r["rank_mass"] for r in rows}) == 1
                        and len({r["rank_flavor"] for r in rows}) == 1,
    }


def measure_the_calibration_dof_census() -> dict:
    """F2 — count what was *selected using flavor data*, by measuring, not naming.

    The measured calibration dimension of a knob group is ``rank J_F``
    restricted to it: how many independent physical flavor directions that
    group can actually move.
    """
    JF = jacobian("flavor")
    JM = jacobian("mass")

    groups = {
        "v3 knobs": list(V3_KNOBS),
        "v4 targeted etas": list(V4_TARGETED_ETAS),
        "eta_k3k5_minus retune": ["eta_k3k5_minus"],
        "diag_shift_plus": list(DIAG_SHIFT_PLUS),
        "diag_shift_minus": list(DIAG_SHIFT_MINUS),
        "all diagonal shifts": list(DIAG_SHIFT_PLUS) + list(DIAG_SHIFT_MINUS),
        "phi_h": ["phi_h"],
        "v4 additions, phi_h fixed": (list(V4_TARGETED_ETAS)
                                      + list(DIAG_SHIFT_PLUS)
                                      + list(DIAG_SHIFT_MINUS)),
        "v4 additions, phi_h released": (list(V4_TARGETED_ETAS)
                                         + list(DIAG_SHIFT_PLUS)
                                         + list(DIAG_SHIFT_MINUS) + ["phi_h"]),
        "everything": list(COORD_NAMES),
    }
    rows = []
    for name, members in groups.items():
        cols = [COORD_NAMES.index(m) for m in members]
        sub = JF[:, cols]
        sv = np.linalg.svd(sub, compute_uv=False)
        rows.append({"group": name, "symbols": len(members),
                     "measured_flavor_rank": int(
                         np.linalg.matrix_rank(sub, tol=1e-8)),
                     "singular_values": sv.tolist()})

    # the trace of each diagonal-shift triple: an exact CKM gauge?
    traces = {}
    for name, members in (("diag_shift_plus", DIAG_SHIFT_PLUS),
                          ("diag_shift_minus", DIAG_SHIFT_MINUS)):
        direction = np.zeros(len(COORD_NAMES))
        for m in members:
            direction[COORD_NAMES.index(m)] = 1.0
        realised = getattr(LOCKED_QUARK_PARAMS_V4, name)
        traces[name] = {
            "flavor_response_norm": float(np.linalg.norm(JF @ direction)),
            "mass_response_norm": float(np.linalg.norm(JM @ direction)),
            "is_an_exact_ckm_gauge": bool(
                np.linalg.norm(JF @ direction) < 1e-8),
            "realised_value": list(realised),
            "realised_trace": float(sum(realised)),
            "realised_is_traceless": bool(abs(sum(realised)) < 1e-9),
        }

    return {
        "rows": rows,
        "diagonal_shift_traces": traces,
        "the_trace_direction_is_flavor_blind": all(
            v["is_an_exact_ckm_gauge"] for v in traces.values()),
        "both_realised_triples_are_traceless": all(
            v["realised_is_traceless"] for v in traces.values()),
        "diag_shift_plus_collapses_further": (
            "the realised (+d, -d, 0) sits in a ONE-parameter family inside the "
            "two-dimensional traceless subspace"),
        "symbol_count_is_not_dof_count": (
            "`diag_shift_plus` and diag_shift_minus carry three symbols each and "
            "measure rank 2 each, because the trace direction moves masses but "
            "no flavor observable and so cannot have been selected on flavor "
            "data"),
    }


def measure_the_mass_preserving_flavor_rank(
        steps: Sequence[float] = (1e-5, 1e-6, 1e-7),
        rconds: Sequence[float] = (1e-6, 1e-8, 1e-10)) -> dict:
    """F3 — the decisive object. ``rank K`` where ``K = J_F ker(J_M)``.

    Rank is invariant under invertible reparameterisation of ``x``, which is
    exactly why it is the right object: no metric, no pseudoinverse, no chart.
    """
    grid = []
    for step in steps:
        for rcond in rconds:
            N = mass_preserving_tangent_space(step, rcond=rcond)
            K, sv = flavor_response_on_mass_preserving_directions(step,
                                                                  rcond=rcond)
            grid.append({"step": step, "rcond": rcond,
                         "kernel_dimension": int(N.shape[1]),
                         "rank_K": int(np.linalg.matrix_rank(K, tol=1e-8)),
                         "singular_values": sv.tolist()})

    K, sv = flavor_response_on_mass_preserving_directions()
    N = mass_preserving_tangent_space()
    JM = jacobian("mass")

    # A TANGENT-SPACE construction -- algebra on the same first-order matrices
    # whose rank was just computed. It is NOT independent evidence, and the
    # first draft presented it as though it were. Kept, correctly labelled.
    rng = np.random.default_rng(1)
    trials = []
    for _ in range(3):
        target = rng.normal(size=len(FLAVOR_OBSERVABLES))
        target /= np.linalg.norm(target)
        coefficients, *_ = np.linalg.lstsq(K, target, rcond=None)
        dx = N @ coefficients
        trials.append({
            "flavor_miss": float(np.linalg.norm(K @ coefficients - target)),
            "mass_disturbance": float(np.max(np.abs(JM @ dx))),
            "parameter_excursion": float(np.linalg.norm(dx))})

    ranks = {row["rank_K"] for row in grid}
    return {
        "grid": grid,
        "rank_K": int(np.linalg.matrix_rank(K, tol=1e-8)),
        "kernel_dimension": int(N.shape[1]),
        "singular_values": sv.tolist(),
        "singular_value_spread": float(sv[0] / sv[-1]) if sv.size else 0.0,
        "kernel_is_a_kernel": float(np.max(np.abs(JM @ N))),
        "rank_stable_across_the_grid": len(ranks) == 1,
        "tangent_space_construction": trials,
        "the_tangent_space_construction_is_not_independent_evidence": (
            "it solves K c = dy_F and reports J_M N c = 0, which is algebra on "
            "the same first-order matrices whose rank was just computed; the "
            "nonlinear check is measure_the_reachability_is_nonlinearly_true"),
        "spans_the_whole_physical_flavor_space": bool(
            int(np.linalg.matrix_rank(K, tol=1e-8))
            == PHYSICAL_FLAVOR_DIMENSION),
        "left_null_vectors": max(
            0, len(FLAVOR_OBSERVABLES)
            - int(np.linalg.matrix_rank(K, tol=1e-8))),
        "so_there_is_no_invariant_relation": (
            "Rank K = 4 leaves no w with w^T K = 0, so the model predicts no "
            "first-order flavor relation w^T dy_F = 0. The round was built to "
            "extract one if it existed"),
        "verdict": ("The mass-preserving parameter freedom spans the entire "
                    "physical flavor space. Fitting the CKM is a successful "
                    "realisation, not a prediction"),
    }


def measure_the_reachability_is_nonlinearly_true(
        epsilons: Sequence[float] = (1e-4, 5e-5, 2.5e-5, 1.25e-5)) -> dict:
    """F3c — the independent nonlinear check the first draft was missing.

    The tangent-space construction solves ``K c = δy_F`` and reports
    ``J_M N c ≈ 0``. That is algebra on the matrices whose rank was just
    computed, not evidence that the *Hamiltonian* does it.

    The real test: take one of those directions, scale it by ``ε``, re-run the
    actual Hamiltonian, and check that the mass error and the flavor miss both
    fall as ``O(ε²)``. If the linearisation is honest, halving ``ε`` quarters
    both.

    It does — ratios ``4.00`` throughout. So the reachability claim is true of
    the model and not only of its Jacobian. What it licenses remains a
    statement about **infinitesimal** CKM directions: the excursion needed for
    a unit ``δy_F`` is ``|δx| ≈ 50–80``, far outside the regime this scaling
    tests.
    """
    JM = jacobian("mass")
    JF = jacobian("flavor")
    N = null_space(JM, rcond=1e-8)
    K = JF @ N
    base_flavor = flavor_observables()
    base_mass = mass_observables()

    rng = np.random.default_rng(2)
    target = rng.normal(size=len(FLAVOR_OBSERVABLES))
    target /= np.linalg.norm(target)
    coefficients, *_ = np.linalg.lstsq(K, target, rcond=None)
    dx = N @ coefficients

    keys = list(SCALAR_COORDS) + list(TUPLE_COORDS)
    rows = []
    for eps in epsilons:
        params = LOCKED_QUARK_PARAMS_V4
        for key, step in zip(keys, eps * dx):
            params = _bump(params, key, float(step))
        mass_error = float(np.max(np.abs(mass_observables(params) - base_mass)))
        miss = float(np.linalg.norm(
            (flavor_observables(params) - base_flavor) - eps * target))
        rows.append({"epsilon": eps, "mass_error": mass_error,
                     "flavor_miss": miss})

    mass_ratios = [a["mass_error"] / b["mass_error"]
                   for a, b in zip(rows[:-1], rows[1:])]
    miss_ratios = [a["flavor_miss"] / b["flavor_miss"]
                   for a, b in zip(rows[:-1], rows[1:])]
    return {
        "rows": rows,
        "mass_error_ratios": mass_ratios,
        "flavor_miss_ratios": miss_ratios,
        "expected_ratio_for_second_order": 4.0,
        "both_are_second_order": bool(
            all(3.6 < r < 4.4 for r in mass_ratios)
            and all(3.6 < r < 4.4 for r in miss_ratios)),
        "parameter_excursion_for_a_unit_target": float(np.linalg.norm(dx)),
        "what_it_licenses": (
            "The claim is about INFINITESIMAL CKM directions reachable locally, "
            "not about arbitrary finite CKM matrices; the excursion needed for "
            "a unit dy_F is far outside the regime this scaling tests"),
    }


def measure_the_restricted_v4_calibration_rank() -> dict:
    """F3b — **the decisive test for predictive surplus**, corrected after review.

    The headline ``rank K`` lets *every* coordinate move, including the eight
    v3 knobs. But the v4 construction was explicitly **additive over a frozen
    v3 lock**, so the question that actually bears on predictive surplus is the
    restricted one:

        ``K_v4 = J_F[:,G] · ker(J_M[:,G])``   with the v3 knobs held fixed.

    The first draft reported ``rank J_F[:,G]`` instead — which lets the masses
    move — and separately reported the unrestricted ``K``. Neither establishes
    that the *actual* v4 calibration freedoms span all four CKM directions
    **while preserving the frozen masses**.

    They do. ``rank K_v4 = 4`` with the v3 knobs fixed and ``φ_h`` fixed at its
    derived value, so the headline is **stronger** than the first draft claimed,
    not weaker.

    ``G`` includes the ``eta_k3k5_minus`` **retune** (``5.0 → 5.586``). The
    round's own rule was to count every numerical quantity selected using flavor
    data regardless of whether the symbol already existed, and the first draft's
    "9 symbols" group broke that rule by omitting it. Including it also improves
    the conditioning markedly — the fourth singular value rises from ``8.0e-6``
    to ``3.1e-3``.
    """
    JM = jacobian("mass")
    JF = jacobian("flavor")
    etas = list(V4_TARGETED_ETAS)
    diag = list(DIAG_SHIFT_PLUS) + list(DIAG_SHIFT_MINUS)
    groups = {
        "v4 additions only, phi_h fixed": etas + diag,
        "v4 additions only, phi_h released": etas + diag + ["phi_h"],
        "with the eta_k3k5 retune, phi_h fixed": (etas + ["eta_k3k5_minus"]
                                                  + diag),
        "with the eta_k3k5 retune, phi_h released": (etas + ["eta_k3k5_minus"]
                                                     + diag + ["phi_h"]),
    }
    rows = []
    for name, members in groups.items():
        cols = [COORD_NAMES.index(m) for m in members]
        N = null_space(JM[:, cols], rcond=1e-8)
        if N.shape[1] == 0:
            rows.append({"group": name, "coordinates": len(members),
                         "kernel_dimension": 0, "rank_K_v4": 0,
                         "singular_values": []})
            continue
        K = JF[:, cols] @ N
        sv = np.linalg.svd(K, compute_uv=False)
        rows.append({"group": name, "coordinates": len(members),
                     "kernel_dimension": int(N.shape[1]),
                     "rank_K_v4": int(np.linalg.matrix_rank(K, tol=1e-8)),
                     "singular_values": sv.tolist(),
                     "smallest_singular_value": float(sv[-1])})

    canonical = next(r for r in rows
                     if r["group"] == "with the eta_k3k5 retune, phi_h fixed")
    bare = next(r for r in rows
                if r["group"] == "v4 additions only, phi_h fixed")
    return {
        "rows": rows,
        "canonical_group": canonical["group"],
        "canonical_coordinates": canonical["coordinates"],
        "rank_K_v4": canonical["rank_K_v4"],
        "every_group_reaches_the_full_flavor_space": all(
            r["rank_K_v4"] == PHYSICAL_FLAVOR_DIMENSION for r in rows),
        "including_the_retune_improves_conditioning": bool(
            canonical["smallest_singular_value"]
            > bare["smallest_singular_value"]),
        "conditioning_improvement": [bare["smallest_singular_value"],
                                     canonical["smallest_singular_value"]],
        "verdict": (
            "With the frozen v3 lock held fixed AND phi_h at its derived value, "
            f"the {canonical['coordinates']} numerical quantities selected in "
            "the v4 flavor re-lock still span all four physical CKM directions "
            "at fixed masses. The surplus claim fails on the restricted "
            "question, not merely on the permissive one"),
        "what_the_first_draft_reported_instead": (
            "It reported rank J_F[:,G], which lets the masses move, alongside an "
            "unrestricted K over all 18 coordinates; neither is the restricted "
            "mass-preserving question"),
    }


def measure_the_phi_h_ab_test() -> dict:
    """F4 — does holding the derived Hopf phase fixed remove a flavor freedom?

    The test the library's own framing demands: ``φ_h = π/k₅`` is called derived
    and is called the source of CP structure. If fixing it dropped the rank to
    3, topology would have removed one calibration freedom.
    """
    released = list(COORD_NAMES)
    fixed = [c for c in COORD_NAMES if c != "phi_h"]
    out = {}
    for label, coords in (("phi_h_released", released), ("phi_h_fixed", fixed)):
        N = mass_preserving_tangent_space(coords=coords)
        K, sv = flavor_response_on_mass_preserving_directions(coords=coords)
        out[label] = {"n_coordinates": len(coords),
                      "kernel_dimension": int(N.shape[1]),
                      "rank_K": int(np.linalg.matrix_rank(K, tol=1e-8)),
                      "singular_values": sv.tolist()}
    lowered = out["phi_h_fixed"]["rank_K"] < out["phi_h_released"]["rank_K"]
    ratio = (out["phi_h_released"]["singular_values"][0]
             / out["phi_h_fixed"]["singular_values"][0])

    # Is it actually a CP direction? Chart-independent, because rescaling the
    # coordinate scales the column without rotating it.
    column = jacobian("flavor")[:, COORD_NAMES.index("phi_h")]
    unit = column / np.linalg.norm(column)
    phi_h_direction = {
        "column": column.tolist(),
        "unit": unit.tolist(),
        "by_observable": dict(zip(FLAVOR_OBSERVABLES, unit.tolist())),
        "delta_share": float(abs(unit[FLAVOR_OBSERVABLES.index("delta")])),
    }
    return {
        "per_case": out,
        "fixing_phi_h_lowers_the_rank": bool(lowered),
        "leading_singular_value_ratio": float(ratio),
        "the_ratio_is_chart_dependent": (
            "singular-value MAGNITUDE is Euclidean in the chosen linear "
            "absolute parameter coordinates; rescaling the phi_h column changes "
            "it. Only the RANK is invariant. Kept as a scoped numerical "
            "diagnostic, not a physical efficiency statement"),
        "phi_h_response_direction": phi_h_direction,
        "phi_h_is_delta_dominated": bool(phi_h_direction["delta_share"] > 0.99),
        "why_the_direction_claim_survives_the_chart": (
            "rescaling the phi_h coordinate scales its J_F column but does not "
            "rotate it, so the share of the response lying along delta is "
            "invariant -- unlike the singular-value ratio"),
        "verdict": (
            "Fixing the derived phase does NOT lower the flavor rank: with "
            "phi_h held at pi/k_5 the other fitted matrix elements still span "
            "all four physical flavor directions at fixed masses. The derived "
            "phase is the single most EFFICIENT CP handle — releasing it "
            f"multiplies the leading singular value by {ratio:.1f} — but "
            "efficiency is not identifiability, and it produces no flavor "
            "prediction by itself"),
        "relation_to_pr_173": (
            "PR #173 found that ADDING phi_h as an input left the observable "
            "rank unchanged. This is the stronger statement: with phi_h held "
            "FIXED, the mass-preserving freedom still covers the whole "
            "physical flavor space"),
    }


def measure_the_counting_claim() -> dict:
    """F5 — test "+3 parameters for +5 independent observables, net +2"."""
    census = measure_the_calibration_dof_census()
    fixed = next(r for r in census["rows"]
                 if r["group"] == "v4 additions, phi_h fixed")
    released = next(r for r in census["rows"]
                    if r["group"] == "v4 additions, phi_h released")
    restricted = measure_the_restricted_v4_calibration_rank()
    canonical = next(r for r in restricted["rows"]
                     if r["group"] == "with the eta_k3k5 retune, phi_h fixed")
    return {
        "restricted_mass_preserving_rank": restricted["rank_K_v4"],
        "restricted_group_coordinates": canonical["coordinates"],
        "the_restricted_question_is_the_decisive_one": (
            "the v4 construction is additive over a FROZEN v3 lock, so the "
            "surplus question is whether the v4 calibration freedoms alone span "
            "the CKM at fixed masses -- rank K_v4, not rank J_F[:,G]"),
        "the_group_includes_the_eta_k3k5_retune": True,
        "claimed_new_parameters": 3,
        "claimed_new_independent_observables": 5,
        "claimed_net_surplus": 2,
        "physical_flavor_dimension": PHYSICAL_FLAVOR_DIMENSION,
        "the_observable_claim_exceeds_the_ceiling": (
            5 > PHYSICAL_FLAVOR_DIMENSION),
        "measured_calibration_dimension_phi_h_fixed":
            fixed["measured_flavor_rank"],
        "measured_calibration_dimension_phi_h_released":
            released["measured_flavor_rank"],
        "symbols_added_phi_h_fixed": fixed["symbols"],
        "measured_net_surplus_phi_h_fixed": (
            PHYSICAL_FLAVOR_DIMENSION - fixed["measured_flavor_rank"]),
        "verdict": (
            "'+5 independent observables' exceeds the ceiling: a unitary 3x3 "
            "CKM has exactly four physical parameters, so at most four of the "
            "nine quoted flavor-CP observables are independent. Against that "
            f"ceiling the v4 calibration group -- the three targeted etas, the "
            f"eta_k3k5 RETUNE, and the six diagonal shifts, "
            f"{canonical['coordinates']} numerical quantities with phi_h fixed "
            f"-- spans rank {restricted['rank_K_v4']} at FIXED v3 masses, so "
            "the net predictive surplus is at most zero, not +2"),
    }


def measure_the_flavor_ledger() -> dict:
    """F6 — what the flavor identifiability audit settles."""
    rank = measure_the_mass_preserving_flavor_rank()
    restricted = measure_the_restricted_v4_calibration_rank()
    ab = measure_the_phi_h_ab_test()
    counting = measure_the_counting_claim()
    census = measure_the_calibration_dof_census()

    entries = [
        {"claim": "the v4 CKM realization is a prediction of the Hamiltonian",
         "verdict": "WITHDRAWN",
         "evidence": f"rank K = {rank['rank_K']} over all coordinates AND "
                     f"rank K_v4 = {restricted['rank_K_v4']} over the v4 "
                     "calibration group alone with the frozen v3 lock and "
                     "phi_h both held fixed -- the restricted question, which "
                     "is the decisive one"},
        {"claim": "the first draft's headline used the right restricted map",
         "verdict": "NO -- THIS ROUND'S OWN GAP",
         "evidence": "it reported rank J_F[:,G] (masses free) beside an "
                     "unrestricted K over all 18 coordinates; neither is the "
                     "mass-preserving restricted map. Computing it gives the "
                     "same rank 4, so the headline strengthens"},
        {"claim": "the tangent-space construction is independent evidence",
         "verdict": "NO -- RELABELLED",
         "evidence": "it is algebra on the same first-order matrices; the "
                     "independent check is an epsilon-scaling re-run of the "
                     "Hamiltonian, which gives clean O(eps^2) (ratios 4.00)"},
        {"claim": "there is a first-order flavor relation to compare with data",
         "verdict": "NONE EXISTS",
         "evidence": f"rank K = 4 leaves {rank['left_null_vectors']} left-null "
                     "vectors, so no w^T dy_F = 0 invariant is predicted"},
        {"claim": "the derived phi_h = pi/k_5 produces a flavor prediction",
         "verdict": "NOT BY ITSELF",
         "evidence": "holding it fixed leaves rank K = 4; it is the most "
                     f"efficient CP handle (leading sv x"
                     f"{ab['leading_singular_value_ratio']:.1f}) but not an "
                     "identifying one"},
        {"claim": "+3 parameters bought +5 independent observables, net +2",
         "verdict": "REFUTED",
         "evidence": "a unitary 3x3 CKM has exactly 4 physical parameters, so "
                     "'+5 independent' exceeds the ceiling; the measured "
                     "calibration dimension of the v4 additions is "
                     f"{counting['measured_calibration_dimension_phi_h_fixed']} "
                     "with phi_h fixed, giving net surplus <= 0"},
        {"claim": "the six diagonal-shift numbers are six fit freedoms",
         "verdict": "MISCOUNTED",
         "evidence": "each triple measures flavor rank 2, not 3: the trace "
                     "direction is an exact CKM gauge that moves masses only, "
                     "and both realised triples are traceless to ~1e-10"},
        {"claim": "the masses and the CKM constrain the parameters jointly",
         "verdict": "NOT AT FIRST ORDER",
         "evidence": f"ker J_M is {rank['kernel_dimension']}-dimensional and "
                     "J_F restricted to it still reaches rank 4; the two "
                     "sectors do not intersect in a constraining way here"},
    ]
    return {
        "entries": entries,
        "headline": (
            "rank K = 4: the mass-preserving parameter freedom spans the whole "
            "physical flavor space, so the v4 CKM agreement is a successful "
            "but locally flexible realisation, not a prediction — and holding "
            "the derived phi_h fixed does not change that"),
        "what_this_does_not_say": (
            "The v4 numbers are not wrong and the <= 1% agreement across nine "
            "observables is real. What the rank shows is that the agreement is "
            "not evidence FOR the Hamiltonian, because the same Hamiltonian "
            "could have reproduced any LOCALLY NEIGHBOURING CKM equally well "
            "at the same masses. The claim is first-order: the epsilon-scaling "
            "check licenses infinitesimal directions, not arbitrary finite "
            "CKM matrices"),
        "scope": (
            "A local, first-order statement at the v4 lock. Rank says nothing "
            "about how far the mass-preserving surface extends before "
            "nonlinearity or positivity bites, and the parameter excursions "
            "required are not small"),
        "the_recommendation": (
            "Per the round's own terms: rank K = 4, so downgrade the CKM "
            "realization to a successful but locally flexible realisation, "
            "stop the quark parameter archaeology, and return to the trunk"),
    }
