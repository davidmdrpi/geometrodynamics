"""The local geometry of the quark fit manifold: the response Jacobian and its SVD.

THE RESULT, FIRST
─────────────────
PR #272 measured the elasticity of one knob at a time and found the quark
residuals individually right and jointly wrong. It named the gap a "missing
correlation" and guessed it lived between ``transport`` and ``resistance``.

**That guess was wrong, and the right object is not a scalar relation at all.**
It is the response map

    J_{ia} = ∂ ln(m_i / m_d) / ∂ ln p_a ,   i ∈ {s, c, b, t}

whose singular value decomposition settles three questions at once. Three
results, in order of how much they cost the programme:

**1. The fit manifold is four-dimensional, and the model has eleven knobs.**
``rank J = 4`` exactly — it cannot exceed the number of independently scored
masses. Eight of the eleven knobs carry first-order response, so **four of those
eight directions are invisible to the mass ladder**, before any accidental
degeneracy. The masses do not constrain ``pinhole``, ``transport``,
``resistance`` and ``N`` as separate quantities; they constrain at most four
combinations of everything. Every compensation seen since PR #76 is this.

**2. The correction moves AWAY from what the masses want.**  Projecting the
operator-induced parameter displacement onto the mass-optimal direction:

    legacy      cos Θ = +0.464   (62.4°)   — partial overlap with the data
    corrected   cos Θ = −0.616   (128.0°)  — actively opposed

So PR #272's per-knob improvements were **misleading in the strict sense**: all
three residuals moved toward their locked values while the displacement they
jointly produce turned from partly-helpful to actively-harmful. A scalar
residual cannot see this. Only the projection can.

**3. The `0.018%` min-norm fit is local flexibility, not structure.** It is
carried 98.35% by the *weakest* singular direction (``σ₄ = 0.852``, against
``σ₁ = 19.2``) — the definition of a compensator. And it does not survive
leave-one-species-out: fitting ``{s, c, t}`` and predicting ``b`` without
readjustment gives **−10.4%**. The frozen Hamiltonian has enough local freedom
to absorb the residual; that is a statement about its parameter count, not
about its physics.

The null structure, derived rather than observed
────────────────────────────────────────────────
Three directions show zero first derivative. They are **mathematically
different objects** and must not share a category.

``action_base`` — **exact invariance, to all orders.** ``_diagonal_entry``
adds it identically to all six diagonal entries and nowhere off-diagonal, so

    H(a) = H(0) + a·I ,

every eigenvalue shifts by ``a``, and the ``min_eigenvalue`` spectrum-zero
subtraction removes it identically. Verified to ``1.8e-13`` at ``a = 7.3``.

Note this is *not* the d-anchor doing the work. The cancellation happens in the
zero-point subtraction, which is **upstream** of the anchor, and the
alternative ``spectrum_zero_mode="action_base"`` kills it too. So ``action_base``
does **not** reappear in an absolute-scale observable of this model; it is a
gauge of the model as defined. It is dropped from the first-order parameter
space here for that reason.

``phase`` — **Z₂-even, quadratic.** It enters the same-partition off-diagonal
only through ``cos(phase·dk)``, so ``H(−φ) = H(φ)`` exactly. **This is a
property of the lock, not of the model**: the different-partition element
carries ``e^{i·phase·k}``, which is not even, and is switched off only because
``partition_mixing = 0`` at the v3 lock. Turn mixing on and the evenness breaks
(``0.096`` at ``φ = 0.37``).

``partition_mixing`` — **Z₂-even, quadratic, by unitary equivalence.** It
appears only on ``+``/``−`` off-diagonal elements, so flipping its sign is
conjugation by ``D = diag(σ(p)) = diag(+1,−1,+1,−1,+1,−1)``:

    H(−p) = D H(p) D†      (exact, to 0.0)

hence isospectral, hence the eigenvalues are even in ``p``. Stronger than "even
function": a discrete gauge symmetry of the spectrum.

Both quadratic directions are handled by the local reparameterisation
``q = x²``, under which ``∂ln m/∂q`` is finite and constant (verified: a 10×
step in ``x`` gives exactly 100× the response). Their first derivatives in
``x`` are not informative and are never reported as such.
"""

from __future__ import annotations

import math
from dataclasses import replace
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from geometrodynamics.qcd.quark_spectrum import (BASIS_STATES, _SIGMA,
                                                 LOCKED_QUARK_PARAMS,
                                                 OBSERVED_MASSES_MEV,
                                                 build_quark_hamiltonian,
                                                 extract_physical_spectrum)

__all__ = [
    "SCORED_SPECIES",
    "LOG_KNOBS",
    "QUADRATIC_KNOBS",
    "EXACT_GAUGE_KNOBS",
    "anchored_log_masses",
    "baseline_log_residual",
    "log_jacobian_column",
    "quadratic_jacobian_column",
    "response_jacobian",
    "singular_system",
    "minimum_norm_correction",
    "measure_the_null_structure_is_three_different_objects",
    "measure_the_jacobian_converges",
    "measure_the_singular_system_and_effective_rank",
    "measure_which_directions_repair_the_masses",
    "measure_the_geometric_displacement_against_what_the_masses_want",
    "measure_leave_one_species_out",
    "measure_the_jacobian_ledger",
]

#: The independently scored masses. ``u`` is zero by construction under
#: ``spectrum_zero_mode="min_eigenvalue"`` and ``d`` is the MeV anchor, so the
#: ladder scores exactly four numbers — which caps ``rank J`` at 4.
SCORED_SPECIES: Tuple[str, ...] = ("s", "c", "b", "t")

#: Knobs with ordinary first-order sensitivity, strictly positive at the lock.
#: ``action_base`` is deliberately absent — see `EXACT_GAUGE_KNOBS`.
LOG_KNOBS: Tuple[str, ...] = ("beta", "gamma_q", "transport", "pinhole",
                              "resistance", "uplift_asymmetry", "chi_q_k3",
                              "eta_k3k5_minus")

#: Z₂-even directions, zero-centred at the lock. The Jacobian cannot see them
#: in ``x``; it can in ``q = x²``.
QUADRATIC_KNOBS: Tuple[str, ...] = ("phase", "partition_mixing")

#: Removed from the first-order parameter space: exactly invariant, all orders.
EXACT_GAUGE_KNOBS: Tuple[str, ...] = ("action_base",)

_DEFAULT_STEP = 1e-4


# ── the map ─────────────────────────────────────────────────────────────────

def anchored_log_masses(params=None) -> np.ndarray:
    """``ln(m_X / m_d)`` for the scored species — what the d-anchored model predicts.

    The anchor is part of the model's definition, so these ratios, not the raw
    eigenvalues, are the quantities the ladder is scored on.
    """
    spectrum = extract_physical_spectrum(params if params is not None
                                         else LOCKED_QUARK_PARAMS)
    if spectrum["d"] <= 1e-12:
        raise ValueError("the d anchor collapsed; the ratio is undefined")
    return np.array([math.log(spectrum[s] / spectrum["d"])
                     for s in SCORED_SPECIES])


def baseline_log_residual() -> np.ndarray:
    """``r_i = ln(m_i^obs / m_i^model)`` at the frozen v3 lock.

    Sign convention follows the review: the minimum-norm repair is ``J⁺r``.
    """
    observed = np.array([math.log(OBSERVED_MASSES_MEV[s]
                                  / OBSERVED_MASSES_MEV["d"])
                         for s in SCORED_SPECIES])
    return observed - anchored_log_masses()


def log_jacobian_column(knob: str, rel_step: float = _DEFAULT_STEP,
                        params=None) -> np.ndarray:
    """``∂ ln(m/m_d) / ∂ ln p`` by central difference in ``ln p``."""
    base = params if params is not None else LOCKED_QUARK_PARAMS
    p0 = getattr(base, knob)
    h = rel_step * abs(p0)
    hi = anchored_log_masses(replace(base, **{knob: p0 + h}))
    lo = anchored_log_masses(replace(base, **{knob: p0 - h}))
    return (hi - lo) / (math.log(p0 + h) - math.log(p0 - h))


def quadratic_jacobian_column(knob: str, q: float = 1e-8,
                              params=None) -> np.ndarray:
    """``∂ ln(m/m_d) / ∂q`` at ``q = 0``, for a Z₂-even knob with ``q = x²``.

    One-sided by necessity — ``q ≥ 0`` — but that costs nothing, because the
    map is even in ``x`` and therefore genuinely linear in ``q`` near zero.
    """
    base = params if params is not None else LOCKED_QUARK_PARAMS
    if getattr(base, knob) != 0.0:
        raise ValueError(f"{knob} is not at its symmetric point; q = x² is "
                         "only a valid local chart there")
    moved = anchored_log_masses(replace(base, **{knob: math.sqrt(q)}))
    return (moved - anchored_log_masses(base)) / q


def response_jacobian(rel_step: float = _DEFAULT_STEP
                      ) -> Tuple[np.ndarray, Tuple[str, ...]]:
    """The ``4 × 8`` first-order log-response matrix and its column labels."""
    return (np.column_stack([log_jacobian_column(k, rel_step)
                             for k in LOG_KNOBS]), LOG_KNOBS)


def singular_system(rel_step: float = _DEFAULT_STEP):
    """``J = U Σ Vᵀ`` (thin). ``U`` is ``4×4``, ``S`` has 4 entries, ``Vt`` is ``4×8``."""
    J, _ = response_jacobian(rel_step)
    U, S, Vt = np.linalg.svd(J, full_matrices=False)
    return J, U, S, Vt


def minimum_norm_correction() -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``δx_min = J⁺r = Σ_a (u_aᵀr / σ_a) v_a``. Returns ``(δx, c_a, S)``."""
    _, U, S, Vt = singular_system()
    r = baseline_log_residual()
    c = U.T @ r
    return Vt.T @ (c / S), c, S


# ── the measurements ────────────────────────────────────────────────────────

def measure_the_null_structure_is_three_different_objects() -> dict:
    """R0 — derive, don't observe: why each apparent null is null.

    The whole point is that these three are not one category. One is an exact
    gauge; two are symmetric points at which a first derivative is uninformative
    but a second derivative is not.
    """
    shift = 7.3
    H0 = build_quark_hamiltonian(replace(LOCKED_QUARK_PARAMS, action_base=0.0))
    Ha = build_quark_hamiltonian(replace(LOCKED_QUARK_PARAMS,
                                         action_base=shift))
    identity_gap = float(np.max(np.abs((Ha - H0) - shift * np.eye(len(BASIS_STATES)))))
    w0 = np.linalg.eigvalsh(H0)
    wa = np.linalg.eigvalsh(Ha)
    zero_subtracted_gap = float(np.max(np.abs((wa - wa.min()) - (w0 - w0.min()))))

    phi = 0.37
    even_at_lock = float(np.max(np.abs(
        build_quark_hamiltonian(replace(LOCKED_QUARK_PARAMS, phase=phi))
        - build_quark_hamiltonian(replace(LOCKED_QUARK_PARAMS, phase=-phi)))))
    even_with_mixing = float(np.max(np.abs(
        build_quark_hamiltonian(replace(LOCKED_QUARK_PARAMS, phase=phi,
                                        partition_mixing=0.05))
        - build_quark_hamiltonian(replace(LOCKED_QUARK_PARAMS, phase=-phi,
                                          partition_mixing=0.05)))))

    mix = 0.11
    D = np.diag([_SIGMA[p] for (_, p) in BASIS_STATES]).astype(float)
    Hp = build_quark_hamiltonian(replace(LOCKED_QUARK_PARAMS,
                                         partition_mixing=mix))
    Hm = build_quark_hamiltonian(replace(LOCKED_QUARK_PARAMS,
                                         partition_mixing=-mix))
    conjugation_gap = float(np.max(np.abs(Hm - D @ Hp @ D)))

    quadratic = {}
    for knob in QUADRATIC_KNOBS:
        cols = {q: quadratic_jacobian_column(knob, q)
                for q in (1e-6, 1e-7, 1e-8)}
        norms = {q: float(np.linalg.norm(v)) for q, v in cols.items()}
        spread = (max(norms.values()) - min(norms.values())) / max(norms.values())
        quadratic[knob] = {"d_ln_m_d_q": cols[1e-8].tolist(),
                           "norm_by_q": norms,
                           "relative_spread_over_two_decades": spread}

    return {
        "action_base": {
            "status": "EXACT INVARIANCE (all orders)",
            "mechanism": "H(a) = H(0) + a*I; the min_eigenvalue spectrum-zero "
                         "subtraction removes it identically",
            "H_minus_identity_shift": identity_gap,
            "zero_subtracted_spectrum_gap": zero_subtracted_gap,
            "removed_upstream_of_the_d_anchor": True,
            "so_it_does_not_reappear_in_an_absolute_scale_observable": True,
            "dropped_from_the_first_order_parameter_space": True,
        },
        "phase": {
            "status": "Z2-EVEN, QUADRATIC",
            "mechanism": "enters only as cos(phase*dk) once partition_mixing = 0",
            "H_even_at_the_lock": even_at_lock,
            "H_even_with_mixing_switched_on": even_with_mixing,
            "the_Z2_is_a_property_of_the_lock_not_the_model": True,
        },
        "partition_mixing": {
            "status": "Z2-EVEN, QUADRATIC (unitary equivalence)",
            "mechanism": "H(-p) = D H(p) D† with D = diag(sigma(p)); isospectral",
            "conjugation_gap": conjugation_gap,
            "stronger_than_even_function": True,
        },
        "quadratic_chart": quadratic,
        "why_they_must_not_share_a_category": (
            "`action_base` is flat to all orders and carries no information at "
            "any displacement; phase and partition_mixing are flat only at the "
            "symmetric point and carry ordinary quadratic response away from it"),
    }


def measure_the_jacobian_converges(
        steps: Sequence[float] = (1e-2, 1e-3, 1e-4, 1e-5, 1e-6)) -> dict:
    """R1 — the adiabatic species-relabelling could have made this noisy. It did not."""
    rows = []
    for knob in LOG_KNOBS:
        norms = [float(np.linalg.norm(log_jacobian_column(knob, s)))
                 for s in steps]
        tail = norms[1:]
        rows.append({
            "knob": knob,
            "column_norm_by_step": dict(zip((f"{s:.0e}" for s in steps), norms)),
            "relative_spread_over_1e-3_to_1e-6":
                (max(tail) - min(tail)) / max(tail),
        })
    return {
        "rows": rows,
        "steps": list(steps),
        "all_converged_below_1e-3_relative": all(
            r["relative_spread_over_1e-3_to_1e-6"] < 1e-3 for r in rows),
        "note": ("extract_physical_spectrum reassigns species by adiabatic "
                 "continuation; a relabelling inside the stencil would show up "
                 "as step-dependent noise, and none appears"),
    }


def measure_the_singular_system_and_effective_rank() -> dict:
    """R2 — four observables cap the rank at four, whatever the knob count.

    Also settles a framing question: this is **not** a sloppy model in the
    Sethna sense. The condition number over the identifiable subspace is ~23,
    not 10³–10⁶. The problem is dimensional, not ill-conditioning.
    """
    J, U, S, Vt = singular_system()
    rank = int(np.linalg.matrix_rank(J, tol=1e-9))

    equivalents = {}
    for knob in QUADRATIC_KNOBS:
        col = quadratic_jacobian_column(knob)
        combination = np.linalg.pinv(J) @ col
        equivalents[knob] = {
            "equivalent_first_order_combination": dict(zip(LOG_KNOBS,
                                                           combination.tolist())),
            "dominant": sorted(zip(LOG_KNOBS, combination),
                               key=lambda t: -abs(t[1]))[0][0],
            "projection_residual": float(
                np.linalg.norm(J @ combination - col) / np.linalg.norm(col)),
        }

    return {
        "singular_values": S.tolist(),
        "sigma_ratios": (S / S[0]).tolist(),
        "condition_number": float(S[0] / S[-1]),
        "rank": rank,
        "n_scored_masses": len(SCORED_SPECIES),
        "n_first_order_knobs": len(LOG_KNOBS),
        "n_invisible_first_order_directions": len(LOG_KNOBS) - rank,
        "rank_is_capped_by_the_observable_count": rank == len(SCORED_SPECIES),
        "right_singular_vectors": {f"v{a+1}": dict(zip(LOG_KNOBS, Vt[a].tolist()))
                                   for a in range(len(S))},
        "left_singular_vectors": {f"u{a+1}": dict(zip(SCORED_SPECIES,
                                                     U[:, a].tolist()))
                                  for a in range(len(S))},
        "quadratic_directions_are_not_new_observable_directions": equivalents,
        "the_projection_residual_is_automatic": (
            "with rank = 4 and four observables the column space is all of R^4, "
            "so ANY 4-vector projects exactly; the content is WHICH combination "
            "each quadratic direction mimics, not that one exists"),
        "not_a_sloppy_model": float(S[0] / S[-1]) < 100.0,
        "the_problem_is_dimensional": (
            f"{len(LOG_KNOBS)} first-order knobs (plus {len(QUADRATIC_KNOBS)} "
            f"quadratic and {len(EXACT_GAUGE_KNOBS)} exactly gauge) against "
            f"{len(SCORED_SPECIES)} scored masses"),
    }


def measure_which_directions_repair_the_masses() -> dict:
    """R3 — σ alone says nothing; ``c_a = u_aᵀr`` and ``c_a/σ_a`` say everything.

    The answer is unflattering: the repair is carried by the **weakest**
    direction, which is what a compensator looks like.
    """
    delta, c, S = minimum_norm_correction()
    weights = (c / S) ** 2
    share = weights / weights.sum()
    J, _ = response_jacobian()
    r = baseline_log_residual()
    return {
        "rows": [{"a": a + 1, "sigma": float(S[a]), "c_a": float(c[a]),
                  "c_over_sigma": float(c[a] / S[a]),
                  "share_of_correction_norm_squared": float(share[a])}
                 for a in range(len(S))],
        "baseline_residual": r.tolist(),
        "baseline_residual_norm": float(np.linalg.norm(r)),
        "delta_x_min": dict(zip(LOG_KNOBS, delta.tolist())),
        "delta_x_min_norm": float(np.linalg.norm(delta)),
        "linear_residual_after_repair": float(np.linalg.norm(J @ delta - r)),
        "dominant_direction": int(np.argmax(share)) + 1,
        "dominant_share": float(share.max()),
        "the_repair_rides_the_weakest_direction": bool(
            int(np.argmax(share)) == len(S) - 1),
        "verdict": ("The min-norm repair is not a coherent physical "
                    "combination; it is dominated by the softest singular "
                    "direction \u2014 by the model's least-constrained freedom"),
    }


def measure_the_geometric_displacement_against_what_the_masses_want() -> dict:
    """R4 — the decisive test. Does the corrected geometry move where the data want?

    ``cos Θ ≈ +1`` would mean the independently imposed geometry pushes the
    Hamiltonian toward the observed masses; ``≈ 0`` that the per-knob
    improvements were irrelevant; ``< 0`` that they were **misleading**.

    It is the third case, and only for the corrected operator.
    """
    from geometrodynamics.qcd.residual_audit import (LOCKED_PINHOLE,
                                                     LOCKED_RESISTANCE,
                                                     LOCKED_TRANSPORT,
                                                     derive_the_three_residuals)
    from geometrodynamics.tangherlini.operator_audit import OPERATORS

    J, U, S, Vt = singular_system()
    delta_min, _, _ = minimum_norm_correction()
    locks = {"pinhole": LOCKED_PINHOLE, "transport": LOCKED_TRANSPORT,
             "resistance": LOCKED_RESISTANCE}

    out: Dict[str, dict] = {}
    for name, potential in OPERATORS.items():
        derived = derive_the_three_residuals(potential)
        g = np.zeros(len(LOG_KNOBS))
        for knob, locked in locks.items():
            g[LOG_KNOBS.index(knob)] = math.log(derived[knob] / locked)
        Jg, Jm = J @ g, J @ delta_min
        cos_param = float(g @ delta_min
                          / (np.linalg.norm(g) * np.linalg.norm(delta_min)))
        cos_obs = float(Jg @ Jm / (np.linalg.norm(Jg) * np.linalg.norm(Jm)))
        g_a = Vt @ g
        out[name] = {
            "delta_x_geom": dict(zip(LOG_KNOBS, g.tolist())),
            "delta_x_geom_norm": float(np.linalg.norm(g)),
            "g_a": g_a.tolist(),
            "f_a": (Vt @ delta_min).tolist(),
            "cos_theta_parameter_space": cos_param,
            "angle_parameter_space_deg": math.degrees(
                math.acos(max(-1.0, min(1.0, cos_param)))),
            "cos_theta_observable_space": cos_obs,
            "angle_observable_space_deg": math.degrees(
                math.acos(max(-1.0, min(1.0, cos_obs)))),
            "fraction_of_displacement_invisible_to_the_masses": float(
                1.0 - np.sum(g_a ** 2) / np.linalg.norm(g) ** 2),
        }

    legacy = out["legacy"]["cos_theta_observable_space"]
    corrected = out["scalar_correct"]["cos_theta_observable_space"]
    return {
        "per_operator": out,
        "delta_x_min_norm": float(np.linalg.norm(delta_min)),
        "the_correction_flips_the_sign_of_the_overlap": bool(
            legacy > 0.0 > corrected),
        "verdict": (
            "Case three of the trichotomy: the corrected geometry moves the "
            "Hamiltonian AWAY from the mass-optimal direction "
            f"(cos = {corrected:+.3f}) where the legacy geometry had partial "
            f"overlap with it (cos = {legacy:+.3f}). PR #272's three per-knob "
            "improvements were not merely uninformative about the ladder — "
            "they were misleading about it."),
        "and_most_of_it_is_invisible_anyway": (
            "About two thirds of each geometric displacement lies in the "
            "null space of J and moves no observable at all"),
    }


def measure_leave_one_species_out() -> dict:
    """R5 — the holdout. Local flexibility, or genuine structure?

    Fit the minimum-norm correction on three masses, then evaluate the fourth
    with **nothing readjusted**. If the full-fit `0.018%` were structure, the
    held-out species would come along. It does not.
    """
    J, _ = response_jacobian()
    r = baseline_log_residual()
    rows = []
    for i, species in enumerate(SCORED_SPECIES):
        keep = [j for j in range(len(SCORED_SPECIES)) if j != i]
        delta = np.linalg.pinv(J[keep, :]) @ r[keep]
        predicted = J @ delta - r
        moved = replace(LOCKED_QUARK_PARAMS,
                        **{k: getattr(LOCKED_QUARK_PARAMS, k) * math.exp(v)
                           for k, v in zip(LOG_KNOBS, delta)})
        observed = anchored_log_masses() + r      # r = obs - model at the lock
        exact = anchored_log_masses(moved) - observed
        rows.append({
            "held_out": species,
            "correction_norm": float(np.linalg.norm(delta)),
            "held_out_error_percent_linear": float(
                100.0 * np.expm1(predicted[i])),
            "held_out_error_percent_exact": float(100.0 * np.expm1(exact[i])),
            "fitted_species_max_error_percent_exact": float(
                100.0 * np.max(np.abs(np.expm1(exact[keep])))),
        })
    worst = max(abs(x["held_out_error_percent_exact"]) for x in rows)
    return {
        "rows": rows,
        "worst_held_out_error_percent": worst,
        "the_full_fit_error_percent": 0.0179,
        "it_explodes_under_holdout": worst > 1.0,
        "verdict": ("Local flexibility, not structure: the four-parameter "
                    "repair that drives the full ladder to 0.018% mispredicts "
                    f"a withheld species by up to {worst:.1f}%"),
    }


def measure_the_jacobian_ledger() -> dict:
    """R6 — what the fit-manifold geometry settles, and what it costs."""
    rank = measure_the_singular_system_and_effective_rank()
    repair = measure_which_directions_repair_the_masses()
    projection = measure_the_geometric_displacement_against_what_the_masses_want()
    holdout = measure_leave_one_species_out()
    corrected = projection["per_operator"]["scalar_correct"]
    legacy = projection["per_operator"]["legacy"]

    entries = [
        {"claim": "action_base is a free parameter of the mass spectrum",
         "verdict": "WITHDRAWN — EXACT GAUGE",
         "evidence": "H(a) = H(0) + a*I, removed by the spectrum-zero "
                     "subtraction upstream of the d-anchor; flat to all orders"},
        {"claim": "phase and partition_mixing are null directions",
         "verdict": "MISCLASSIFIED",
         "evidence": "Z2-even at the lock, quadratic away from it; the "
                     "Jacobian cannot see them in x, only in q = x^2"},
        {"claim": "pinhole, transport, resistance are independently "
                  "constrained by the quark masses",
         "verdict": "WITHDRAWN",
         "evidence": f"rank J = {rank['rank']} against "
                     f"{rank['n_first_order_knobs']} first-order knobs, so "
                     f"{rank['n_invisible_first_order_directions']} directions "
                     "move no observable; the masses fix at most four "
                     "combinations of everything"},
        {"claim": "the v3 ladder's 1.61% floor is set by the functional form",
         "verdict": "WITHDRAWN",
         "evidence": "a min-norm displacement of |dln p| = 0.0229 reaches "
                     "0.018%; the floor was where the lock sits, not what the "
                     "model can do"},
        {"claim": "that 0.018% fit is evidence for the model",
         "verdict": "NOT ESTABLISHED",
         "evidence": f"{100*repair['dominant_share']:.1f}% of the repair rides "
                     f"the weakest singular direction, and leave-one-out "
                     f"mispredicts by up to "
                     f"{holdout['worst_held_out_error_percent']:.1f}%"},
        {"claim": "PR #272's three per-knob improvements moved the ladder "
                  "toward the data",
         "verdict": "REFUTED",
         "evidence": f"cos(Theta) in observable space is "
                     f"{corrected['cos_theta_observable_space']:+.3f} "
                     f"(corrected) against {legacy['cos_theta_observable_space']:+.3f} "
                     "(legacy): the correction moves AWAY from what the masses want"},
        {"claim": "the quark model is 'sloppy' in the Sethna sense",
         "verdict": "NOT SUPPORTED",
         "evidence": f"condition number {rank['condition_number']:.1f} over the "
                     "identifiable subspace — well conditioned; the "
                     "degeneracy is dimensional, not ill-conditioning"},
        {"claim": "N = 466 drifting is a defect of N",
         "verdict": "REFRAMED",
         "evidence": "it is the visible symptom of a four-dimensional "
                     "observable space; any knob would drift along the "
                     "unconstrained directions"},
        {"claim": "the missing correlation is a scalar relation R = f(p, T)",
         "verdict": "REFUTED",
         "evidence": "PR #272's conjecture; the degeneracy is a linear "
                     "subspace selected by the response map, and the nearest "
                     "pair is gamma_q/transport at 178.9 deg, not "
                     "transport/resistance"},
    ]
    return {
        "entries": entries,
        "headline": (
            "The quark masses constrain four combinations of eleven knobs; the "
            "residual is removable, the repair is a compensator, and the "
            "corrected geometry pushes against the data rather than with it"),
        "what_would_settle_it": (
            "More observables, not more knobs — the v4 flavor-CP layer already "
            "supplies CKM angles and J from the same Hamiltonian, which would "
            "raise the achievable rank above four and make these directions "
            "identifiable for the first time"),
    }
