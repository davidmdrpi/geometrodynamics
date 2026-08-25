"""What changed when the radial scalar operator was corrected, and what did not.

The correction
──────────────
Until PR #271 `tangherlini.radial.V_tangherlini` carried

    V_legacy = A[ ℓ(ℓ+2)/r² + 3 r_h²/r⁴ ] ,        A = 1 − r_h²/r²

which is **not** the master potential of a minimally coupled massless scalar.
For ``ds² = −A dt² + A⁻¹dr² + r²dΩ_n²`` such a scalar obeys

    (1/r^n) ∂_r(r^n A R') + (ω²/A − ℓ(ℓ+n−1)/r²) R = 0 ,

and the unique first-derivative-free Schrödinger form, with ``dr* = dr/A`` and
``ψ = r^{n/2}R``, carries

    V_scalar = A[ ℓ(ℓ+n−1)/r² + n(n−2)A/(4r²) + n A'/(2r) ] ,

verified symbolically at ``n = 2 … 6`` and independently by the flat limit,
where ``ψ = r^{1/2}J_{ℓ+1}(ωr)`` gives ``V → ((ℓ+1)² − ¼)/r²``.  At ``n = 3``:

    V_scalar − V_legacy = 3A²/(4r²) .

**This is a bug, not a convention.**  The old name implied a canonical scalar
operator that the implementation did not provide.

Why the shift is not uniform in what it breaks
──────────────────────────────────────────────
``ΔV`` carries no ``ℓ``.  Two consequences run in opposite directions:

* **Differences are exactly invariant.**  ``V_{ℓ₂} − V_{ℓ₁}`` is unchanged to
  ``3.6e-15`` across every pair to ``ℓ = 5``, so anything built on the cross-ℓ
  operator survives *algebraically*.  Its matrix elements still drift, because
  the eigenfunctions drift — structure invariant, numbers shifted.
* **Barrier heights are not.**  ``V_max`` and its sums move by more than the
  eigenvalues do, because a barrier height reads the potential directly while
  an eigenvalue averages it against a bound state.

So the audit has three verdicts, not two, and the third is the one that matters:

    EXACTLY INVARIANT   — algebraically unchanged, to machine precision
    NUMERICALLY SHIFTED — same claim, different digits
    INTERPRETATION CHANGED — the claim itself no longer reads the same way

The γ story, precisely
──────────────────────
Two different statements about ``γ = 22.5`` live in the tree and they move in
**opposite directions**:

* the canonical one (README: *"Pinhole γ ≈ Σ V_max[1..5] … −2.2 % off the locked
  γ = 22.5"*) **improves**, from ``22.008`` (−2.2 %) to ``22.331`` (−0.75 %);
* the ``ℓ = 0…5`` statement — that adding the ``ℓ = 0`` 5D channel closes the
  gap — **breaks**, from ``22.453`` (−0.21 %) to ``22.836`` (+1.50 %).

Under the corrected operator the sum closest to ``22.5`` is the ``1…5`` one, so
the ℓ = 0 channel no longer closes anything; it overshoots. That claim has to be
reopened, not re-fitted.

What is deliberately **not** re-run
───────────────────────────────────
Topological results have no dependence on this radial operator and are not
touched: the Hopf fibration, the Pin⁻ structure on the exchange ``ℝP²``, the
odd-``k`` winding ladder, antipodal parity.  Proximity is not dependence.
"""

from __future__ import annotations

import math
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.linalg import eig as scipy_eig
from scipy.optimize import brentq

from geometrodynamics.constants import R_MID, R_OUTER
from geometrodynamics.tangherlini.radial import (V_scalar_tangherlini,
                                                 V_tangherlini_legacy,
                                                 _cheb_diff, r_to_rstar,
                                                 rstar_to_r)

__all__ = [
    "LOCKED_GAMMA",
    "OPERATORS",
    "chebyshev_grid",
    "eigen_solve",
    "vmax_sum",
    "r_outer_fixed_point",
    "radial_action",
    "measure_the_two_operators_and_their_exact_gap",
    "measure_the_eigenvalue_shifts",
    "measure_the_gamma_sums_and_the_r_outer_fixed_point",
    "measure_what_survives_exactly",
    "measure_the_eigenvector_derived_quantities",
    "measure_the_wkb_action_shift",
    "measure_the_downstream_ledger",
    "measure_which_geometry_preserves_the_lepton_ladder",
]

LOCKED_GAMMA = 22.5

OPERATORS: Dict[str, Callable] = {
    "legacy": V_tangherlini_legacy,
    "scalar_correct": V_scalar_tangherlini,
}


# ── the shared numerical spine ──────────────────────────────────────────────

def chebyshev_grid(points: int = 80, rs: float = R_MID,
                   r_outer: float = R_OUTER) -> Tuple[np.ndarray, np.ndarray, float]:
    """The canonical tortoise Chebyshev grid used by `solve_radial_modes`.

    Reproduced here rather than imported so the audit can vary ``r_outer``
    without touching the production solver.
    """
    rs_min = r_to_rstar(rs + 5e-4, rs)
    rs_max = r_to_rstar(r_outer - 5e-4, rs)
    x, D = _cheb_diff(points)
    half = (rs_max - rs_min) / 2.0
    rsg = rs_min + half * (1.0 - x)
    rg = np.array([rstar_to_r(s, rs) for s in rsg])
    return rg, D, half


def eigen_solve(ell: int, potential: Callable, points: int = 80,
                rs: float = R_MID, r_outer: float = R_OUTER,
                n_modes: int = 4) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``H = −∂²_{r*} + V`` on the canonical grid, for whichever ``V`` is given."""
    rg, D, half = chebyshev_grid(points, rs, r_outer)
    V = np.asarray(potential(rg, ell, rs), dtype=float)
    H = -(1.0 / half ** 2) * (D @ D) + np.diag(V)
    ev, evec = scipy_eig(H[1:points, 1:points])
    ev, evec = np.real(ev), np.real(evec)
    keep = ev > 0
    ev, evec = ev[keep], evec[:, keep]
    order = np.argsort(ev)[:n_modes]
    vecs = evec[:, order]
    # fix the sign convention so overlaps are comparable between operators
    for k in range(vecs.shape[1]):
        col = vecs[:, k]
        if abs(col.min()) > col.max():
            vecs[:, k] = -col
        vecs[:, k] /= np.linalg.norm(vecs[:, k])
    return np.sqrt(ev[order]), vecs, rg


def vmax_sum(potential: Callable, lo: int, hi: int, points: int = 80,
             rs: float = R_MID, r_outer: float = R_OUTER) -> float:
    """``Σ_ℓ max_r V(r, ℓ)`` on the canonical grid — the γ / pinhole quantity."""
    rg, _, _ = chebyshev_grid(points, rs, r_outer)
    return float(sum(float(np.max(np.asarray(potential(rg, l, rs), dtype=float)))
                     for l in range(lo, hi + 1)))


def r_outer_fixed_point(potential: Callable, lo: int, hi: int,
                        target: float = LOCKED_GAMMA,
                        bracket: Tuple[float, float] = (1.10, 1.60)
                        ) -> Optional[float]:
    """The ``R_OUTER`` at which the barrier sum equals the locked ``γ``."""
    def f(R):
        return vmax_sum(potential, lo, hi, r_outer=R) - target
    try:
        return float(brentq(f, bracket[0], bracket[1], xtol=1e-7))
    except Exception:
        return None


def radial_action(omega: float, potential: Callable, ell: int,
                  points: int = 4000, rs: float = R_MID,
                  r_outer: float = R_OUTER) -> float:
    """``∮ √(ω² − V) dr*`` over the classically allowed region — the WKB action.

    Computed in tortoise coordinate, which is where the operator is a plain
    Schrödinger problem, so the action is the honest closed-orbit invariant.
    """
    rs_min = r_to_rstar(rs + 5e-4, rs)
    rs_max = r_to_rstar(r_outer - 5e-4, rs)
    grid = np.linspace(rs_min, rs_max, points)
    rphys = np.array([rstar_to_r(s, rs) for s in grid])
    V = np.asarray(potential(rphys, ell, rs), dtype=float)
    integrand = omega ** 2 - V
    allowed = integrand > 0.0
    if not np.any(allowed):
        return 0.0
    return float(np.trapezoid(np.sqrt(integrand[allowed]), grid[allowed]))


# ── measurements ────────────────────────────────────────────────────────────

def measure_the_two_operators_and_their_exact_gap(
        rs: float = R_MID, ells: Sequence[int] = (0, 1, 2, 3, 4, 5)
) -> Dict[str, object]:
    """The correction itself: the gap, its ``ℓ``-independence, and the proof.

    Three independent confirmations, none of them a citation:

    * the gap equals ``3A²/(4r²)`` in closed form;
    * it carries no ``ℓ``, so every difference ``V_{ℓ₂} − V_{ℓ₁}`` is invariant;
    * the flat limit reproduces the Bessel form, which is what settles *which*
      operator is the scalar one.
    """
    from scipy.special import jv  # noqa: F401  (documented, used below)

    r = np.linspace(1.02 * rs, 6.0 * rs, 500)
    A = 1.0 - (rs / r) ** 2
    closed_form = 3.0 * A ** 2 / (4.0 * r ** 2)

    gaps, gap_err = [], 0.0
    for l in ells:
        g = (np.asarray(V_scalar_tangherlini(r, l, rs))
             - np.asarray(V_tangherlini_legacy(r, l, rs)))
        gaps.append(g)
        gap_err = max(gap_err, float(np.max(np.abs(g - closed_form))))

    ell_spread = float(np.max([np.max(np.abs(g - gaps[0])) for g in gaps]))

    worst_difference = 0.0
    for i, l1 in enumerate(ells):
        for l2 in ells[i + 1:]:
            d_old = (np.asarray(V_tangherlini_legacy(r, l2, rs))
                     - np.asarray(V_tangherlini_legacy(r, l1, rs)))
            d_new = (np.asarray(V_scalar_tangherlini(r, l2, rs))
                     - np.asarray(V_scalar_tangherlini(r, l1, rs)))
            worst_difference = max(worst_difference,
                                   float(np.max(np.abs(d_new - d_old))))

    # the flat limit, against Bessel, with no r_h anywhere
    r2 = np.linspace(0.5, 8.0, 300)
    flat_err = 0.0
    for l in ells:
        got = np.asarray(V_scalar_tangherlini(r2, l, 1e-9))
        want = ((l + 1) ** 2 - 0.25) / r2 ** 2
        flat_err = max(flat_err, float(np.max(np.abs(got - want) / want)))

    from geometrodynamics.tangherlini.dynamics import master_potential
    bitwise = all(
        np.array_equal(np.asarray(V_scalar_tangherlini(r, l, rs)),
                       np.asarray(master_potential(r, l, rs)))
        for l in ells)

    return {
        "legacy": "A[l(l+2)/r^2 + 3 r_h^2/r^4]",
        "scalar_correct": "A[l(l+n-1)/r^2 + n(n-2)A/(4r^2) + n A'/(2r)]",
        "the_gap": "3 A^2 / (4 r^2)",
        "gap_matches_the_closed_form": gap_err,
        "the_gap_carries_no_ell": ell_spread,
        "cross_ell_differences_are_invariant": worst_difference,
        "flat_limit_matches_bessel": flat_err,
        "agrees_bitwise_with_dynamics_master_potential": bool(bitwise),
        "why_it_is_a_bug_not_a_convention": "the old name implied the canonical "
                                            "minimally coupled scalar operator, "
                                            "and the implementation was short of "
                                            "it by an l-independent term",
    }


def measure_the_eigenvalue_shifts(
        ells: Sequence[int] = (0, 1, 2, 3, 4, 5), n_modes: int = 4
) -> Dict[str, object]:
    """``ω_{ℓn}`` and eigenfunction overlap under both operators.

    The eigenvalues move at the ``10⁻³`` level and *decrease* in sensitivity
    with ``ℓ`` — an eigenvalue averages the potential against a bound state, so
    an ``ℓ``-independent shift matters least where the centrifugal term is
    largest.  That is the reassuring half of the audit.
    """
    rows = []
    for l in ells:
        w_old, v_old, _ = eigen_solve(l, V_tangherlini_legacy, n_modes=n_modes)
        w_new, v_new, _ = eigen_solve(l, V_scalar_tangherlini, n_modes=n_modes)
        k = min(len(w_old), len(w_new), n_modes)
        overlaps = [abs(float(np.dot(v_old[:, j], v_new[:, j])))
                    for j in range(k)]
        rows.append({
            "ell": l,
            "omega_legacy": [float(x) for x in w_old[:k]],
            "omega_correct": [float(x) for x in w_new[:k]],
            "ground_shift_percent": float(100.0 * (w_new[0] - w_old[0]) / w_old[0]),
            "min_eigenfunction_overlap": float(min(overlaps)),
        })
    shifts = [abs(r["ground_shift_percent"]) for r in rows]
    return {
        "rows": rows,
        "omega_1_0_legacy": rows[1]["omega_legacy"][0],
        "omega_1_0_correct": rows[1]["omega_correct"][0],
        "largest_ground_shift_percent": max(shifts),
        "smallest_ground_shift_percent": min(shifts),
        "all_shifts_below_a_fifth_of_a_percent": bool(max(shifts) < 0.2),
        "sensitivity_falls_with_ell": bool(
            all(b <= a + 1e-12 for a, b in zip(shifts, shifts[1:]))),
        "eigenfunctions_barely_move": bool(
            min(r["min_eigenfunction_overlap"] for r in rows) > 0.99),
        "why_so_small": "an eigenvalue averages the potential against a bound "
                        "state, so an l-independent shift matters least where "
                        "the centrifugal term already dominates",
    }


def measure_the_gamma_sums_and_the_r_outer_fixed_point() -> Dict[str, object]:
    """The barrier sums and the geometric root — where the interpretation moves.

    A barrier height reads the potential directly rather than averaging it, so
    ``V_max`` shifts by far more than ``ω`` does.  This is the measurement that
    turns the correction from a digit change into a claim change.
    """
    rows = []
    for lo, hi in ((1, 5), (0, 5)):
        old = vmax_sum(V_tangherlini_legacy, lo, hi)
        new = vmax_sum(V_scalar_tangherlini, lo, hi)
        rows.append({
            "channels": f"l = {lo}..{hi}",
            "sum_legacy": old,
            "sum_correct": new,
            "residual_legacy_percent": 100.0 * (old - LOCKED_GAMMA) / LOCKED_GAMMA,
            "residual_correct_percent": 100.0 * (new - LOCKED_GAMMA) / LOCKED_GAMMA,
            "r_outer_legacy": r_outer_fixed_point(V_tangherlini_legacy, lo, hi),
            "r_outer_correct": r_outer_fixed_point(V_scalar_tangherlini, lo, hi),
        })
    one_five, zero_five = rows[0], rows[1]
    closest_legacy = min(rows, key=lambda x: abs(x["residual_legacy_percent"]))
    closest_correct = min(rows, key=lambda x: abs(x["residual_correct_percent"]))
    return {
        "rows": rows,
        "locked_gamma": LOCKED_GAMMA,
        "the_canonical_readme_claim_improves": bool(
            abs(one_five["residual_correct_percent"])
            < abs(one_five["residual_legacy_percent"])),
        "canonical_residual_before": one_five["residual_legacy_percent"],
        "canonical_residual_after": one_five["residual_correct_percent"],
        "the_ell_zero_closure_claim_breaks": bool(
            abs(zero_five["residual_correct_percent"])
            > abs(zero_five["residual_legacy_percent"])),
        "ell_zero_residual_before": zero_five["residual_legacy_percent"],
        "ell_zero_residual_after": zero_five["residual_correct_percent"],
        "closest_to_gamma_before": closest_legacy["channels"],
        "closest_to_gamma_after": closest_correct["channels"],
        "the_closest_channel_set_swaps": bool(
            closest_legacy["channels"] != closest_correct["channels"]),
        "what_has_to_be_reopened": "the claim that adding the l = 0 5D channel "
                                   "closes the gamma discrepancy -- under the "
                                   "corrected operator it overshoots, and the "
                                   "l = 1..5 sum is the one near 22.5",
        "and_it_is_not_a_refit": "the canonical l = 1..5 statement moves the "
                                 "right way on its own; nothing was tuned",
    }


def measure_what_survives_exactly(
        ells: Sequence[int] = (0, 1, 2, 3, 4, 5)) -> Dict[str, object]:
    """The algebra that the correction cannot touch, separated from the numbers.

    ``ΔV`` is ``ℓ``-independent, so the cross-``ℓ`` perturbation operator
    ``V_{ℓ+2} − V_ℓ`` is unchanged exactly.  Its *matrix elements* are not,
    because they are taken between eigenfunctions that drift — which is the
    distinction the ledger below turns on.
    """
    rg, _, _ = chebyshev_grid()
    worst_operator = 0.0
    for l in ells[:-2]:
        d_old = (np.asarray(V_tangherlini_legacy(rg, l + 2, R_MID))
                 - np.asarray(V_tangherlini_legacy(rg, l, R_MID)))
        d_new = (np.asarray(V_scalar_tangherlini(rg, l + 2, R_MID))
                 - np.asarray(V_scalar_tangherlini(rg, l, R_MID)))
        worst_operator = max(worst_operator, float(np.max(np.abs(d_new - d_old))))

    # the matrix elements, which do move
    elements = []
    for l in (0, 1, 2, 3):
        _, v_old, _ = eigen_solve(l, V_tangherlini_legacy, n_modes=1)
        _, v_new, _ = eigen_solve(l, V_scalar_tangherlini, n_modes=1)
        d_old = (np.asarray(V_tangherlini_legacy(rg, l + 2, R_MID))
                 - np.asarray(V_tangherlini_legacy(rg, l, R_MID)))[1:80]
        d_new = (np.asarray(V_scalar_tangherlini(rg, l + 2, R_MID))
                 - np.asarray(V_scalar_tangherlini(rg, l, R_MID)))[1:80]
        m_old = float(np.dot(v_old[:, 0] ** 2, d_old))
        m_new = float(np.dot(v_new[:, 0] ** 2, d_new))
        elements.append({"ell": l, "element_legacy": m_old,
                         "element_correct": m_new,
                         "drift_percent": 100.0 * (m_new - m_old) / abs(m_old)})
    drifts = [abs(e["drift_percent"]) for e in elements]
    return {
        "the_cross_ell_operator_is_unchanged": worst_operator,
        "matrix_elements": elements,
        "largest_element_drift_percent": max(drifts),
        "structure_invariant_numbers_shifted": True,
        "the_partition": "the operator V_{l+2} - V_l survives algebraically "
                         "EXACTLY; its matrix elements drift because the "
                         "eigenfunctions do",
    }


def measure_the_eigenvector_derived_quantities() -> Dict[str, object]:
    """Throat-flux ratios: quantities built as ratios, which mostly cancel.

    ``α_q(ℓ,0)`` is a ratio of throat derivatives normalised to ``ℓ = 1``, so
    the common part of the shift divides out and only the differential survives.
    Measured rather than assumed — a ratio is not automatically safe.
    """
    from geometrodynamics.tangherlini.alpha_q import throat_du_dr

    def table(potential):
        modes = {}
        for l in range(0, 6):
            w, v, rg = eigen_solve(l, potential, n_modes=1)
            u = np.zeros(81)
            u[1:80] = v[:, 0]
            u /= abs(u).max() + 1e-12
            modes[l] = {"funcs": [{"u_half": u, "r_phys": rg}]}
        ref = abs(throat_du_dr(modes[1]["funcs"][0]))
        return {l: float(throat_du_dr(modes[l]["funcs"][0]) / ref)
                for l in modes}

    old, new = table(V_tangherlini_legacy), table(V_scalar_tangherlini)
    rows = [{"ell": l, "alpha_q_legacy": old[l], "alpha_q_correct": new[l],
             "drift_percent": (100.0 * (new[l] - old[l]) / abs(old[l])
                               if abs(old[l]) > 1e-12 else 0.0)}
            for l in sorted(old)]
    drifts = [abs(r["drift_percent"]) for r in rows]
    return {
        "rows": rows,
        "largest_drift_percent": max(drifts),
        "the_reference_mode_is_exactly_one": bool(
            abs(new[1] - 1.0) < 1e-12 and abs(old[1] - 1.0) < 1e-12),
        "ratios_absorb_most_of_the_shift": bool(max(drifts) < 5.0),
        "caveat": "a ratio is not automatically safe -- the common part "
                  "cancels, the differential does not, and this is measured "
                  "rather than assumed",
    }


def measure_the_wkb_action_shift(
        ells: Sequence[int] = (1, 2, 3)) -> Dict[str, object]:
    """Closed-orbit radial actions, which read the barrier rather than average it.

    Each action is evaluated at *its own* operator's ground frequency, so the
    comparison is between two self-consistent orbits and not between one orbit
    and the other operator's potential.
    """
    rows = []
    for l in ells:
        w_old, _, _ = eigen_solve(l, V_tangherlini_legacy, n_modes=1)
        w_new, _, _ = eigen_solve(l, V_scalar_tangherlini, n_modes=1)
        a_old = radial_action(float(w_old[0]), V_tangherlini_legacy, l)
        a_new = radial_action(float(w_new[0]), V_scalar_tangherlini, l)
        rows.append({"ell": l, "action_legacy": a_old, "action_correct": a_new,
                     "drift_percent": 100.0 * (a_new - a_old) / abs(a_old)})
    drifts = [abs(r["drift_percent"]) for r in rows]
    return {
        "rows": rows,
        "largest_drift_percent": max(drifts),
        "each_action_uses_its_own_ground_frequency": True,
        "why_it_moves_more_than_omega": "the action integrates sqrt(w^2 - V) "
                                        "against the potential directly, "
                                        "without the bound state's averaging",
    }


def measure_the_downstream_ledger() -> Dict[str, object]:
    """Every load-bearing radial claim, sorted into the three verdicts.

    The point of the round: **not** whether the old tests stay green, but which
    published claims are algebraically untouched, which keep their meaning with
    different digits, and which no longer say what they said.
    """
    gam = measure_the_gamma_sums_and_the_r_outer_fixed_point()
    eig = measure_the_eigenvalue_shifts()
    inv = measure_what_survives_exactly()
    rat = measure_the_eigenvector_derived_quantities()
    act = measure_the_wkb_action_shift()

    entries = [
        {"claim": "cross-l perturbation operator V_{l+2} - V_l",
         "verdict": "EXACTLY INVARIANT",
         "evidence": f"unchanged to {inv['the_cross_ell_operator_is_unchanged']:.1e}"},
        {"claim": "Hopf fibration, Pin- structure, odd-k ladder, antipodal parity",
         "verdict": "EXACTLY INVARIANT",
         "evidence": "no dependence on the radial operator; not re-run, and "
                     "proximity is not dependence"},
        {"claim": "alpha_q(l,0) throat-flux ratios",
         "verdict": "NUMERICALLY SHIFTED",
         "evidence": f"largest drift {rat['largest_drift_percent']:.3f}%"},
        {"claim": "omega_{l,n} radial eigenfrequencies and eigenfunctions",
         "verdict": "NUMERICALLY SHIFTED",
         "evidence": f"ground shifts "
                     f"{eig['smallest_ground_shift_percent']:.4f}% to "
                     f"{eig['largest_ground_shift_percent']:.4f}%; overlaps "
                     f"> 0.99"},
        {"claim": "cross-l transport matrix elements",
         "verdict": "NUMERICALLY SHIFTED",
         "evidence": f"operator exact, elements drift up to "
                     f"{inv['largest_element_drift_percent']:.3f}%"},
        {"claim": "closed-orbit / WKB radial actions",
         "verdict": "NUMERICALLY SHIFTED",
         "evidence": f"largest drift {act['largest_drift_percent']:.3f}%"},
        {"claim": "the 1.054 factor, omega(1,0) at the gamma-locked geometry",
         "verdict": "NUMERICALLY SHIFTED",
         "evidence": f"{eig['omega_1_0_legacy']:.6f} -> "
                     f"{eig['omega_1_0_correct']:.6f}; the quoted 1.054 becomes "
                     f"1.056, which exceeds the 0.04% Compton-bridge tolerance "
                     f"and so needs re-quoting, not re-deriving"},
        {"claim": "pinhole gamma = Sum V_max[1..5] vs the locked 22.5",
         "verdict": "NUMERICALLY SHIFTED",
         "evidence": f"{gam['canonical_residual_before']:+.2f}% -> "
                     f"{gam['canonical_residual_after']:+.2f}% -- the canonical "
                     f"README claim IMPROVES, and nothing was tuned"},
        {"claim": "R_OUTER as the gamma = 22.5 fixed point",
         "verdict": "NUMERICALLY SHIFTED",
         "evidence": f"l=0..5 root {gam['rows'][1]['r_outer_legacy']:.5f} -> "
                     f"{gam['rows'][1]['r_outer_correct']:.5f}; l=1..5 root "
                     f"{gam['rows'][0]['r_outer_legacy']:.5f} -> "
                     f"{gam['rows'][0]['r_outer_correct']:.5f}"},
        {"claim": "'adding the l = 0 5D channel closes the gamma discrepancy'",
         "verdict": "INTERPRETATION CHANGED",
         "evidence": f"{gam['ell_zero_residual_before']:+.2f}% -> "
                     f"{gam['ell_zero_residual_after']:+.2f}%; the l=0..5 sum "
                     f"now OVERSHOOTS and the l=1..5 sum is the one near 22.5 "
                     f"-- the closest channel set swaps"},
        {"claim": "any generation or mass result whose chain runs through "
                  "gamma, R_OUTER, or a radial eigenvalue",
         "verdict": "INTERPRETATION CHANGED",
         "evidence": "inherits the gamma reopening above; each such chain must "
                     "be re-derived from the corrected operator before its "
                     "number is quoted again"},
    ]
    verdicts = {v: sum(1 for e in entries if e["verdict"] == v)
                for v in ("EXACTLY INVARIANT", "NUMERICALLY SHIFTED",
                          "INTERPRETATION CHANGED")}
    return {
        "entries": entries,
        "counts": verdicts,
        "the_question_this_answers": "not whether the old tests stay green, but "
                                     "which published claims are algebraically "
                                     "untouched, which keep their meaning with "
                                     "different digits, and which no longer say "
                                     "what they said",
        "not_re_run_and_why": "Hopf, Pin-, odd-k winding and antipodal parity "
                              "have no dependence on the radial scalar "
                              "operator; proximity is not dependence",
        "what_is_still_open": "the gamma narrative. the l = 0 closure claim is "
                              "withdrawn, not replaced -- the corrected l = 1..5 "
                              "sum lands nearer 22.5 than the old one did, but "
                              "that is an observation, not a derivation of why "
                              "22.5",
    }


def measure_which_geometry_preserves_the_lepton_ladder() -> Dict[str, object]:
    """The one narrow downstream re-derivation, run before the correction merges.

    The corrected operator exposed an ambiguity the legacy one hid: is the lock
    supposed to be ``ℓ = 0..5``, ``ℓ = 1..5``, or was ``γ = 22.5`` itself partly
    fitted?  Recomputing every lepton number before answering that would just
    propagate an unresolved branch choice, so exactly three geometries are passed
    through the **locked** lepton Hamiltonian with **no retuning**:

        A — keep ``R_OUTER = 1.26`` fixed and let ``γ`` float
        B — enforce ``Σ_{1..5} V_max = 22.5``
        C — enforce ``Σ_{0..5} V_max = 22.5``

    **The answer is not one of the three anticipated outcomes.**
    `lepton_spectrum.compute_knotted_lepton_spectrum` discards ``r_outer``
    outright (``del l, n_points, rs, r_outer``); the locked block consumes the
    geometry **only** through the scalar ``hard_pinhole_gamma``.  So B and C come
    out *bit-identical*, and the channel-set choice leaves no trace in any
    observable once ``γ`` is enforced.

    > **``γ = 22.5`` is the selector; ``R_OUTER`` is downstream of it.**

    Which inverts the first branch of the trichotomy: it is *fixing* ``R_OUTER``
    and letting ``γ`` float that breaks the ladder, by ``+15%`` or ``−21%``
    depending on channel set.

    And the correction **weakens** the "geometry supplies ``γ``" story even while
    improving the ``1..5`` residual in isolation.  Under the legacy operator the
    canonical ``R = 1.26`` geometry produced ``γ[0..5] = 22.453``, ``0.21%`` from
    the lock, and the masses landed within ``3.8%`` — that near-coincidence *was*
    the claim.  Under the corrected operator no channel set at ``R = 1.26`` lands
    near ``22.5``, and both damage the ladder at the ``15–21%`` level.

    The reason is sensitivity: ``d ln m_μ / d ln γ ≈ −19``, so a sub-percent
    geometric residual is **not** a small residual in this chain.

    **So the channel-set question is not decidable by the lepton observables.**
    It has to be settled by what ``γ`` means geometrically, because the masses
    only ever saw the scalar.
    """
    from geometrodynamics.tangherlini.lepton_spectrum import (
        LEPTON_BASELINE_DEPTHS, LEPTON_BASELINE_PHASE,
        LEPTON_BASELINE_RESISTANCE, LEPTON_BASELINE_TRANSPORT, S3_ACTION_BASE,
        TAU_BETA_50PI, calibrate_electron_predict_heavier)

    def locked(gamma: float):
        """The locked block, with gamma the only thing allowed to vary."""
        return calibrate_electron_predict_heavier(
            depths=LEPTON_BASELINE_DEPTHS, phase_per_pass=LEPTON_BASELINE_PHASE,
            transport_strength=LEPTON_BASELINE_TRANSPORT,
            hard_pinhole_gamma=gamma,
            resistance_scale=LEPTON_BASELINE_RESISTANCE,
            resistance_model="exponential", depth_cost_mode="tunnel_only",
            winding_mode="max", action_base=S3_ACTION_BASE,
            k_uplift_beta=TAU_BETA_50PI, n_points=24)

    observed = {3: 105.6583755, 5: 1776.86}
    r_b = r_outer_fixed_point(V_scalar_tangherlini, 1, 5)
    r_c = r_outer_fixed_point(V_scalar_tangherlini, 0, 5)
    cases = [
        ("baseline, gamma = 22.5", None, LOCKED_GAMMA),
        ("legacy R=1.26, gamma[0..5]", 1.26,
         vmax_sum(V_tangherlini_legacy, 0, 5, r_outer=1.26)),
        ("legacy R=1.26, gamma[1..5]", 1.26,
         vmax_sum(V_tangherlini_legacy, 1, 5, r_outer=1.26)),
        ("A corrected R=1.26, gamma[1..5]", 1.26,
         vmax_sum(V_scalar_tangherlini, 1, 5, r_outer=1.26)),
        ("A corrected R=1.26, gamma[0..5]", 1.26,
         vmax_sum(V_scalar_tangherlini, 0, 5, r_outer=1.26)),
        ("B corrected root [1..5] = 22.5", r_b, LOCKED_GAMMA),
        ("C corrected root [0..5] = 22.5", r_c, LOCKED_GAMMA),
    ]
    rows = []
    for name, r_out, gamma in cases:
        fit = locked(gamma)
        mu, tau = fit.predicted_mev[3], fit.predicted_mev[5]
        rows.append({
            "case": name, "r_outer": r_out, "gamma": float(gamma),
            "m_mu": float(mu), "m_tau": float(tau),
            "mu_error_percent": float(100.0 * (mu - observed[3]) / observed[3]),
            "tau_error_percent": float(100.0 * (tau - observed[5]) / observed[5]),
        })
    by = {r["case"]: r for r in rows}
    b, c = by["B corrected root [1..5] = 22.5"], by["C corrected root [0..5] = 22.5"]
    a1 = by["A corrected R=1.26, gamma[1..5]"]
    a0 = by["A corrected R=1.26, gamma[0..5]"]
    leg0 = by["legacy R=1.26, gamma[0..5]"]

    # d ln m / d ln gamma, from the two corrected fixed-R geometries
    elasticity = float((math.log(a0["m_mu"]) - math.log(a1["m_mu"]))
                       / (math.log(a0["gamma"]) - math.log(a1["gamma"])))

    return {
        "rows": rows,
        "the_three_geometries": "A: R_OUTER = 1.26 fixed; B: enforce "
                                "Sum[1..5] = 22.5; C: enforce Sum[0..5] = 22.5",
        "nothing_was_retuned": True,
        "B_and_C_are_bit_identical": bool(b["m_mu"] == c["m_mu"]
                                          and b["m_tau"] == c["m_tau"]),
        "why": "compute_knotted_lepton_spectrum discards r_outer outright; the "
               "locked block consumes the geometry ONLY through the scalar "
               "hard_pinhole_gamma",
        "so_the_channel_set_is_invisible_to_the_observables": True,
        "gamma_is_the_selector_r_outer_is_downstream": True,
        "fixing_r_outer_breaks_the_ladder": bool(
            abs(a1["mu_error_percent"]) > 10.0
            and abs(a0["mu_error_percent"]) > 10.0),
        "corrected_fixed_R_mu_errors": [a1["mu_error_percent"],
                                        a0["mu_error_percent"]],
        "legacy_fixed_R_mu_error": leg0["mu_error_percent"],
        "the_correction_weakens_the_geometry_supplies_gamma_story": bool(
            abs(leg0["mu_error_percent"]) < abs(a1["mu_error_percent"])
            and abs(leg0["mu_error_percent"]) < abs(a0["mu_error_percent"])),
        "d_ln_m_mu_over_d_ln_gamma": elasticity,
        "so_a_subpercent_residual_is_not_small": bool(abs(elasticity) > 10.0),
        "what_this_does_not_settle": "the channel-set question is not decidable "
                                     "by the lepton observables -- it has to be "
                                     "settled by what gamma MEANS geometrically, "
                                     "because the masses only ever saw the scalar",
        "the_outcome_was_not_one_of_the_three_anticipated": True,
    }
