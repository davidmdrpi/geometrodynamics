"""The local geometry of the quark fit manifold: the response Jacobian and its SVD.

THE DURABLE RESULT, FIRST
─────────────────────────
PR #272's residual round measured the elasticity of one knob at a time, found
the quark residuals individually right and jointly wrong, and guessed the gap
was a scalar relation between ``transport`` and ``resistance``. That guess was
wrong, and so was the object: a collection of scalars cannot express what is
happening. The right object is the response map

    J_{ia} = ∂ ln(m_i / m_d) / ∂ ln p_a ,   i ∈ {s, c, b, t}

Three statements survive review, and they are the ones to build on because none
of them depends on a choice of metric in parameter space:

**1. Four scored masses cannot identify the current parameterisation.**
``rank J = 4`` exactly — capped by the observable count — while eight knobs
carry first-order response and three more are gauge or quadratic. ``J^T J`` has
four exact zero eigenvalues. So ``pinhole``, ``transport``, ``resistance`` and
``N`` are not separately constrained by the mass ladder; at most four
combinations are.

**2. The positive-knob Jacobian is numerically full row rank**, converged to
five digits across three decades of step size, with no adiabatic-relabelling
noise. The response map is a real derivative, not a difference artefact.

**3. Per-knob proximity to a fitted value does not determine the effect on the
spectrum.** The three radial residuals do not compose linearly into the ladder,
and their individual agreement carries no information about the sign or size of
their joint effect.

What the review corrected, and how
──────────────────────────────────
Four things published in the first draft of this module were wrong or
overclaimed. They are recorded here because the programme's own rule is that
wrong first drafts get written down.

**The ``phase`` symmetry is model-wide, not lock-only.** The first draft said
the ``Z₂`` was "a property of the lock", because the *matrix* stops satisfying
``H(−φ) = H(φ)`` once ``partition_mixing ≠ 0``. But the Hamiltonian satisfies

    H(−φ, p) = H(φ, p)*        for arbitrary p        (exact, to 0.0)

— same-partition entries are real ``cos(φ·dk)``, different-partition entries
``e^{iφk}`` conjugate — and since ``H`` is Hermitian, ``H*`` = ``Hᵀ`` has
identical eigenvalues. The **spectrum** is even in ``φ`` for every ``p``,
verified to ``0.0`` on the eigenvalues and on the anchored masses at
``p = 0.05, 0.3``. The first draft mistook "matrix not equal" for "spectrum not
even". The correct distinction is between two kinds of ``Z₂``:

    phase             antiunitary (complex-conjugation) Z₂ of the spectrum
    partition_mixing  unitary-conjugation Z₂,  H(−p) = D H(p) D†

both quadratic at zero, for different reasons.

**The headline projection used the wrong displacement.** ``g_legacy`` and
``g_corrected`` are two *candidate displacements from the fitted lock*, not the
move the operator correction makes. That move is

    Δg = g_corrected − g_legacy ,

and the question "did correcting the operator move the geometry toward what the
masses still needed?" compares ``J·Δg`` against ``r − J·g_legacy``. It gives
**cos = +0.873**: the correction points strongly *toward* the residual the
legacy triple leaves. The first draft's ``cos = −0.616`` is a true statement
about ``g_corrected`` as a displacement from the lock, and a false basis for
"the correction moves away from the data".

The two objectives also disagree, which is the substance rather than a
technicality:

    L2 log-residual  |r − Jg|   0.0548 (legacy) → 0.0433 (corrected)   IMPROVES
    max relative error          3.44%  (legacy) → 3.78%  (corrected)   WORSENS

Both are reported. Neither is privileged, and no claim is made that rests on
one without saying so.

**The min-norm geometry is metric-dependent.** ``δx_min``, the ``98.35%``
share, the right singular vectors, the parameter-space cosine and the
"invisible fraction" are all Euclidean in the chosen eight log coordinates.
Rescaling any column changes them. They are scoped throughout to *the eight
positive knobs in unit log coordinates with the quadratic coordinates held
fixed*, and the physical conclusions are carried by the observable-space
quantities ``J·g`` and ``r``, which are chart-independent.

**Leave-one-species-out was not a predictive holdout.** For each holdout
``J_keep`` is ``3×8`` with a five-dimensional kernel, and because the full ``J``
has rank 4 the held-out row is *not* in the span of the other three — so there
is always a ``z ∈ ker(J_keep)`` moving the withheld species. Explicitly:
``δ + λz`` fits the held-out mass to ``~1e-15`` while holding the other three
exact, at ``1.02–4.77×`` the norm. The ``−10.4%`` therefore measures the
**minimum-log-norm regulariser**, not the Hamiltonian. Rank deficiency already
established the local flexibility; the holdout added an overclaim on top and is
now reported as a regulariser diagnostic.

The null structure, derived rather than observed
────────────────────────────────────────────────
``action_base`` — **exact invariance, all orders.** ``_diagonal_entry`` adds it
identically to all six diagonal entries and nowhere off-diagonal, so
``H(a) = H(0) + a·I``; every eigenvalue shifts by ``a`` and the
``min_eigenvalue`` spectrum-zero subtraction removes it. A ``3×`` change moves
the anchored masses by ``1.6e-12``. This is *not* the d-anchor doing the work —
the cancellation is upstream of it, and ``spectrum_zero_mode="action_base"``
kills it too — so it does not reappear in an absolute-scale observable of this
model. Dropped from the first-order parameter space.

``phase``, ``partition_mixing`` — the two ``Z₂``s above, handled by the local
chart ``q = x²`` under which ``∂ln m/∂q`` is finite and constant (a 10× step in
``x`` gives 100× the response, to 2%). Their first derivatives in ``x`` are
uninformative and are never reported as such.

On "sloppiness"
───────────────
The ratio of the four **nonzero** singular values is ``22.6``, so there is no
long hierarchy among the identifiable directions. That is not the same as "not
a sloppy model": ``J^T J`` has four exact zeros, and the full eleven-coordinate
model has more local non-identifiability still. The dominant pathology is
**structural non-identifiability**, not ill-conditioning — and the null-space
content varies enormously by coordinate, from ``uplift_asymmetry`` at
``1.0000`` identifiable to ``eta_k3k5_minus`` at ``0.0003``. "Any knob would
drift" is false; which knobs drift is itself structure.
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
    "measure_the_minimum_norm_regulariser_does_not_predict_a_held_out_mass",
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

    Three directions, three mechanisms, and the ``Z₂``s are of two different
    kinds — one antiunitary, one unitary.
    """
    shift = 7.3
    n = len(BASIS_STATES)
    H0 = build_quark_hamiltonian(replace(LOCKED_QUARK_PARAMS, action_base=0.0))
    Ha = build_quark_hamiltonian(replace(LOCKED_QUARK_PARAMS,
                                         action_base=shift))
    identity_gap = float(np.max(np.abs((Ha - H0) - shift * np.eye(n))))
    w0, wa = np.linalg.eigvalsh(H0), np.linalg.eigvalsh(Ha)
    zero_subtracted_gap = float(np.max(np.abs((wa - wa.min())
                                              - (w0 - w0.min()))))

    # phase: H(-phi, p) = conj(H(phi, p)) for ARBITRARY p, so the SPECTRUM is
    # even everywhere -- the matrix equality that fails at p != 0 is not the
    # symmetry that matters.
    conj_rows, matrix_rows, spectrum_rows = [], [], []
    for p in (0.0, 0.05, 0.3):
        for phi in (0.1, 0.37):
            Hp = build_quark_hamiltonian(replace(LOCKED_QUARK_PARAMS,
                                                 phase=phi,
                                                 partition_mixing=p))
            Hm = build_quark_hamiltonian(replace(LOCKED_QUARK_PARAMS,
                                                 phase=-phi,
                                                 partition_mixing=p))
            conj_rows.append(float(np.max(np.abs(Hm - np.conj(Hp)))))
            matrix_rows.append(float(np.max(np.abs(Hm - Hp))))
            spectrum_rows.append(float(np.max(np.abs(
                np.linalg.eigvalsh(Hp) - np.linalg.eigvalsh(Hm)))))
            if p != 0.0:
                spectrum_rows.append(float(np.max(np.abs(
                    anchored_log_masses(replace(LOCKED_QUARK_PARAMS, phase=phi,
                                                partition_mixing=p))
                    - anchored_log_masses(replace(LOCKED_QUARK_PARAMS,
                                                  phase=-phi,
                                                  partition_mixing=p))))))

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
            "status": "ANTIUNITARY Z2 OF THE SPECTRUM, QUADRATIC",
            "mechanism": "H(-phi, p) = conj(H(phi, p)) for ARBITRARY p; H "
                         "Hermitian so conj(H) = H^T is isospectral",
            "max_conjugation_gap_over_all_mixings": max(conj_rows),
            "max_matrix_gap_over_all_mixings": max(matrix_rows),
            "max_spectrum_and_anchored_mass_gap": max(spectrum_rows),
            "the_symmetry_is_model_wide_not_lock_only": True,
            "corrected_from_the_first_draft": (
                "the first draft called this a property of the lock because "
                "the MATRIX equality fails at partition_mixing != 0; the "
                "SPECTRUM is even for every mixing, which is what the masses "
                "see"),
        },
        "partition_mixing": {
            "status": "UNITARY-CONJUGATION Z2, QUADRATIC",
            "mechanism": "H(-p) = D H(p) D† with D = diag(sigma(p)); isospectral",
            "conjugation_gap": conjugation_gap,
            "stronger_than_even_function": True,
        },
        "the_two_Z2s_are_different_kinds": (
            "phase is an antiunitary (complex-conjugation) symmetry of the "
            "spectrum; partition_mixing is a unitary-conjugation symmetry. "
            "Both give a vanishing first derivative and a live second one, for "
            "different reasons"),
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

    **Scope.** Every right-space quantity below is Euclidean in the eight
    positive knobs in unit log coordinates, with the quadratic coordinates held
    fixed. Rescaling a column changes the singular vectors and the min-norm
    solution. The chart-independent statements are the rank, the zero count of
    ``J^T J``, and anything computed in observable space.
    """
    J, U, S, Vt = singular_system()
    rank = int(np.linalg.matrix_rank(J, tol=1e-9))
    gram_eigs = np.linalg.eigvalsh(J.T @ J)
    exact_zeros = int(np.sum(np.abs(gram_eigs) < 1e-12))

    projector = Vt.T @ Vt
    identifiable_share = {}
    for j, knob in enumerate(LOG_KNOBS):
        e = np.zeros(len(LOG_KNOBS))
        e[j] = 1.0
        identifiable_share[knob] = float(e @ projector @ e)

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
        "nonzero_singular_value_ratio": float(S[0] / S[-1]),
        "gram_eigenvalues": gram_eigs.tolist(),
        "gram_exact_zeros": exact_zeros,
        "rank": rank,
        "n_scored_masses": len(SCORED_SPECIES),
        "n_first_order_knobs": len(LOG_KNOBS),
        "n_first_order_null_directions": len(LOG_KNOBS) - rank,
        "rank_is_capped_by_the_observable_count": rank == len(SCORED_SPECIES),
        "identifiable_share_by_knob": identifiable_share,
        "most_identifiable_knob": max(identifiable_share,
                                      key=identifiable_share.get),
        "least_identifiable_knob": min(identifiable_share,
                                       key=identifiable_share.get),
        "right_singular_vectors": {f"v{a+1}": dict(zip(LOG_KNOBS, Vt[a].tolist()))
                                   for a in range(len(S))},
        "left_singular_vectors": {f"u{a+1}": dict(zip(SCORED_SPECIES,
                                                     U[:, a].tolist()))
                                  for a in range(len(S))},
        "quadratic_directions_add_no_observable_direction": equivalents,
        "the_projection_residual_is_automatic": (
            "with rank = 4 and four observables the column space is all of "
            "R^4, so ANY 4-vector projects exactly; this does NOT justify "
            "omitting the quadratic coordinates from a claim about the full "
            "11-knob fit manifold, and no such claim is made"),
        "no_long_hierarchy_among_nonzero_singular_values": float(S[0] / S[-1]) < 100.0,
        "the_dominant_pathology_is_structural_non_identifiability": (
            f"the four nonzero singular values span only "
            f"{float(S[0]/S[-1]):.1f}x, but J^T J carries {exact_zeros} exact "
            f"zeros in these {len(LOG_KNOBS)} coordinates and the full "
            f"11-coordinate model carries more; the problem is rank, not "
            f"conditioning"),
        "which_knobs_drift_is_itself_structure": (
            f"null-space content varies from {identifiable_share['uplift_asymmetry']:.4f} "
            f"(uplift_asymmetry, fully identifiable) to "
            f"{identifiable_share['eta_k3k5_minus']:.4f} (eta_k3k5_minus, "
            f"essentially entirely null) -- 'any knob would drift' is false"),
        "metric_scope": ("Euclidean in the 8 positive knobs in unit log "
                         "coordinates, quadratic coordinates held fixed"),
    }


def measure_which_directions_repair_the_masses() -> dict:
    """R3 — ``c_a = u_aᵀr`` and ``c_a/σ_a``, plus the exact nonlinear re-run.

    **Scope.** The share decomposition is metric-dependent (see
    `measure_the_singular_system_and_effective_rank`). The exact re-run at the
    end is not: it applies ``exp(δ)`` to the lock and re-solves.
    """
    delta, c, S = minimum_norm_correction()
    weights = (c / S) ** 2
    share = weights / weights.sum()
    J, _ = response_jacobian()
    r = baseline_log_residual()

    # the exact nonlinear consequence -- computed, never quoted from prose
    moved = replace(LOCKED_QUARK_PARAMS,
                    **{k: getattr(LOCKED_QUARK_PARAMS, k) * math.exp(v)
                       for k, v in zip(LOG_KNOBS, delta)})
    observed = anchored_log_masses() + r
    exact_residual = anchored_log_masses(moved) - observed
    exact_max = float(100.0 * np.max(np.abs(np.expm1(exact_residual))))
    locked_max = float(100.0 * np.max(np.abs(np.expm1(r))))
    largest_move = float(100.0 * max(abs(math.expm1(v)) for v in delta))

    return {
        "rows": [{"a": a + 1, "sigma": float(S[a]), "c_a": float(c[a]),
                  "c_over_sigma": float(c[a] / S[a]),
                  "share_of_correction_norm_squared": float(share[a])}
                 for a in range(len(S))],
        "baseline_residual": r.tolist(),
        "baseline_residual_norm": float(np.linalg.norm(r)),
        "locked_max_rel_err_percent": locked_max,
        "delta_x_min": dict(zip(LOG_KNOBS, delta.tolist())),
        "delta_x_min_norm": float(np.linalg.norm(delta)),
        "largest_single_knob_move_percent": largest_move,
        "linear_residual_after_repair": float(np.linalg.norm(J @ delta - r)),
        "exact_max_rel_err_percent_after_repair": exact_max,
        "exact_signed_residual_after_repair": exact_residual.tolist(),
        "dominant_direction": int(np.argmax(share)) + 1,
        "dominant_share": float(share.max()),
        "the_repair_rides_the_weakest_direction": bool(
            int(np.argmax(share)) == len(S) - 1),
        "verdict": ("The frozen Hamiltonian has enough local freedom to absorb "
                    f"essentially the whole residual: {locked_max:.2f}% -> "
                    f"{exact_max:.4f}% under a displacement whose largest "
                    f"single knob moves {largest_move:.2f}%. That is a "
                    "statement about parameter count, not about the physics"),
        "metric_caveat": ("the share decomposition and delta_x_min are "
                          "Euclidean in unit log coordinates; rescaling a "
                          "column changes them. The exact re-run does not "
                          "depend on that choice"),
    }


def measure_the_geometric_displacement_against_what_the_masses_want() -> dict:
    """R4 — does the operator correction move the geometry toward the data?

    Two different questions live here and the first draft ran them together.

    * ``g_legacy`` and ``g_corrected`` are two **candidate displacements from
      the fitted lock**. Comparing ``J·g`` to ``r`` asks which candidate lands
      nearer the data.
    * ``Δg = g_corrected − g_legacy`` is **the move the correction actually
      makes**. Comparing ``J·Δg`` to ``r − J·g_legacy`` asks whether correcting
      the operator pushed toward the residual the legacy triple left.

    They give opposite signs, and the second is the one that answers "did the
    correction help": **cos = +0.873**.
    """
    from geometrodynamics.qcd.residual_audit import (LOCKED_PINHOLE,
                                                     LOCKED_RESISTANCE,
                                                     LOCKED_TRANSPORT,
                                                     derive_the_three_residuals)
    from geometrodynamics.tangherlini.operator_audit import OPERATORS

    J, U, S, Vt = singular_system()
    r = baseline_log_residual()
    delta_min, _, _ = minimum_norm_correction()
    locks = {"pinhole": LOCKED_PINHOLE, "transport": LOCKED_TRANSPORT,
             "resistance": LOCKED_RESISTANCE}

    def cosine(a, b):
        return float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)))

    g: Dict[str, np.ndarray] = {}
    per_operator: Dict[str, dict] = {}
    for name, potential in OPERATORS.items():
        derived = derive_the_three_residuals(potential)
        v = np.zeros(len(LOG_KNOBS))
        for knob, locked in locks.items():
            v[LOG_KNOBS.index(knob)] = math.log(derived[knob] / locked)
        g[name] = v
        residual = r - J @ v
        per_operator[name] = {
            "delta_x_geom": dict(zip(LOG_KNOBS, v.tolist())),
            "delta_x_geom_norm": float(np.linalg.norm(v)),
            "cos_to_r_observable_space": cosine(J @ v, r),
            "cos_to_delta_x_min_parameter_space": cosine(v, delta_min),
            "l2_log_residual": float(np.linalg.norm(residual)),
            "max_rel_err_percent_linear": float(
                100.0 * np.max(np.abs(np.expm1(-residual)))),
            "fraction_of_displacement_in_the_null_space": float(
                1.0 - np.sum((Vt @ v) ** 2) / np.linalg.norm(v) ** 2),
        }

    move = g["scalar_correct"] - g["legacy"]
    residual_after_legacy = r - J @ g["legacy"]
    cos_move = cosine(J @ move, residual_after_legacy)

    l2 = {k: per_operator[k]["l2_log_residual"] for k in per_operator}
    mx = {k: per_operator[k]["max_rel_err_percent_linear"] for k in per_operator}
    return {
        "per_operator": per_operator,
        "the_operator_induced_move": {
            "delta_g": dict(zip(LOG_KNOBS, move.tolist())),
            "delta_g_norm": float(np.linalg.norm(move)),
            "residual_left_by_the_legacy_triple_norm": float(
                np.linalg.norm(residual_after_legacy)),
            "cos_move_vs_residual_after_legacy": cos_move,
            "angle_deg": math.degrees(math.acos(max(-1.0, min(1.0, cos_move)))),
        },
        "the_two_objectives_disagree": {
            "l2_log_residual": l2,
            "max_rel_err_percent_linear": mx,
            "l2_improves": l2["scalar_correct"] < l2["legacy"],
            "max_rel_err_worsens": mx["scalar_correct"] > mx["legacy"],
        },
        "verdict": (
            "The correction, taken as the move it actually makes, points "
            f"toward the residual the legacy triple leaves (cos = {cos_move:+.3f}), "
            "and lowers the L2 log-residual "
            f"({l2['legacy']:.4f} -> {l2['scalar_correct']:.4f}). It "
            "simultaneously worsens the repository's max-relative-error score "
            f"({mx['legacy']:.2f}% -> {mx['scalar_correct']:.2f}%). Both are "
            "true; the direction of 'improvement' is objective-dependent and "
            "no metric-free claim is available"),
        "withdrawn": (
            "The first draft's headline -- 'the corrected geometry moves AWAY "
            "from what the masses want', from cos(J g_corrected, r) = -0.616 "
            "-- is a true statement about g_corrected as a displacement FROM "
            "THE LOCK, and an invalid basis for a claim about what correcting "
            "the operator does. Withdrawn"),
        "what_survives": (
            "Per-knob proximity to a fitted value carries no information about "
            "the sign or size of the joint effect on the spectrum -- which was "
            "the point, and does not need the sign of any cosine"),
    }


def measure_the_minimum_norm_regulariser_does_not_predict_a_held_out_mass() -> dict:
    """R5 — a diagnostic of the regulariser, **not** a predictive holdout.

    For each holdout ``J_keep`` is ``3×8`` with a five-dimensional kernel, and
    the full ``J`` has rank 4, so the held-out row is not in the span of the
    other three: there is always a ``z ∈ ker(J_keep)`` that moves the withheld
    species. ``δ + λz`` fits it to machine precision while keeping the other
    three exact. So the pseudoinverse's miss measures the **minimum-log-norm
    choice**, not the model's predictive content.

    Reported because the first draft read it the other way.
    """
    from scipy.linalg import null_space

    J, _ = response_jacobian()
    r = baseline_log_residual()
    rows = []
    for i, species in enumerate(SCORED_SPECIES):
        keep = [j for j in range(len(SCORED_SPECIES)) if j != i]
        J_keep = J[keep, :]
        delta = np.linalg.pinv(J_keep) @ r[keep]
        kernel = null_space(J_keep)
        reach = J[i, :] @ kernel
        miss = float((J @ delta - r)[i])

        repaired = None
        if np.linalg.norm(reach) > 1e-12:
            lam = -miss / float(reach @ reach) * reach
            delta2 = delta + kernel @ lam
            resid2 = J @ delta2 - r
            repaired = {
                "held_out_error_percent": float(100.0 * np.expm1(resid2[i])),
                "fitted_three_max_error_percent": float(
                    100.0 * np.max(np.abs(np.expm1(resid2[keep])))),
                "norm_ratio_to_min_norm": float(
                    np.linalg.norm(delta2) / np.linalg.norm(delta)),
            }
        rows.append({
            "held_out": species,
            "kernel_dimension": int(kernel.shape[1]),
            "held_out_row_reachable_from_kernel": float(np.linalg.norm(reach)),
            "min_norm_correction_norm": float(np.linalg.norm(delta)),
            "min_norm_held_out_error_percent": float(100.0 * np.expm1(miss)),
            "after_a_kernel_shift": repaired,
        })

    return {
        "rows": rows,
        "every_held_out_species_is_reachable": all(
            r_["held_out_row_reachable_from_kernel"] > 1e-9 for r_ in rows),
        "a_kernel_shift_fits_the_held_out_mass_exactly": all(
            r_["after_a_kernel_shift"] is not None
            and abs(r_["after_a_kernel_shift"]["held_out_error_percent"]) < 1e-9
            for r_ in rows),
        "at_norm_cost_between": [
            min(r_["after_a_kernel_shift"]["norm_ratio_to_min_norm"]
                for r_ in rows),
            max(r_["after_a_kernel_shift"]["norm_ratio_to_min_norm"]
                for r_ in rows)],
        "verdict": ("This is a property of the minimum-log-norm regulariser, "
                    "not of the Hamiltonian. The rank deficiency already "
                    "establishes the local flexibility; reading these misses "
                    "as failed predictions was an overclaim and is withdrawn"),
    }


def measure_the_jacobian_ledger() -> dict:
    """R6 — what the fit-manifold geometry settles, what it costs, and what
    this round had to withdraw."""
    rank = measure_the_singular_system_and_effective_rank()
    repair = measure_which_directions_repair_the_masses()
    projection = measure_the_geometric_displacement_against_what_the_masses_want()
    holdout = measure_the_minimum_norm_regulariser_does_not_predict_a_held_out_mass()
    move = projection["the_operator_induced_move"]
    objectives = projection["the_two_objectives_disagree"]

    entries = [
        {"claim": "four scored masses can identify the current quark "
                  "parameterisation",
         "verdict": "WITHDRAWN",
         "evidence": f"rank J = {rank['rank']} against "
                     f"{rank['n_first_order_knobs']} first-order knobs; "
                     f"J^T J carries {rank['gram_exact_zeros']} exact zeros in "
                     "these coordinates and the full 11-knob model more"},
        {"claim": "action_base is a free parameter of the mass spectrum",
         "verdict": "WITHDRAWN — EXACT GAUGE",
         "evidence": "H(a) = H(0) + a*I, removed by the spectrum-zero "
                     "subtraction upstream of the d-anchor; flat to all orders"},
        {"claim": "phase and partition_mixing are null directions",
         "verdict": "MISCLASSIFIED",
         "evidence": "two different Z2s -- phase antiunitary on the spectrum, "
                     "partition_mixing unitary-conjugation -- both quadratic; "
                     "the Jacobian sees them in q = x^2, not in x"},
        {"claim": "the phase Z2 is a property of the lock, not the model",
         "verdict": "WITHDRAWN — THIS ROUND'S OWN ERROR",
         "evidence": "H(-phi, p) = conj(H(phi, p)) for arbitrary p, so the "
                     "SPECTRUM is even at every mixing (0.0 on eigenvalues and "
                     "anchored masses); the first draft mistook matrix "
                     "inequality for spectral asymmetry"},
        {"claim": "per-knob proximity to a fitted value determines the effect "
                  "on the spectrum",
         "verdict": "REFUTED",
         "evidence": "the three radial residuals do not compose linearly into "
                     "the ladder; individual agreement carries no information "
                     "about the sign or size of the joint effect"},
        {"claim": "the corrected geometry moves AWAY from what the masses want",
         "verdict": "WITHDRAWN — THIS ROUND'S OWN ERROR",
         "evidence": "that used g_corrected as a displacement from the lock. "
                     "The move the correction makes is Delta g = g_corrected - "
                     f"g_legacy, and cos(J.Delta g, r - J.g_legacy) = "
                     f"{move['cos_move_vs_residual_after_legacy']:+.3f} -- "
                     "toward the residual, not away"},
        {"claim": "one objective settles whether the correction helped",
         "verdict": "NOT AVAILABLE",
         "evidence": f"L2 log-residual improves "
                     f"({objectives['l2_log_residual']['legacy']:.4f} -> "
                     f"{objectives['l2_log_residual']['scalar_correct']:.4f}) "
                     f"while max relative error worsens "
                     f"({objectives['max_rel_err_percent_linear']['legacy']:.2f}% "
                     f"-> {objectives['max_rel_err_percent_linear']['scalar_correct']:.2f}%); "
                     "the direction of improvement is objective-dependent"},
        {"claim": "the v3 ladder's 1.61% floor is set by the functional form",
         "verdict": "WITHDRAWN",
         "evidence": f"a displacement whose largest single knob moves "
                     f"{repair['largest_single_knob_move_percent']:.2f}% "
                     f"reaches {repair['exact_max_rel_err_percent_after_repair']:.4f}% "
                     "on an exact nonlinear re-run"},
        {"claim": "that repair is evidence for the model",
         "verdict": "NOT ESTABLISHED",
         "evidence": f"{100*repair['dominant_share']:.1f}% of it rides the "
                     "weakest singular direction in the chosen log chart, and "
                     "the rank deficiency makes some such displacement "
                     "guaranteed to exist"},
        {"claim": "leave-one-species-out shows the Hamiltonian fails to predict",
         "verdict": "WITHDRAWN — THIS ROUND'S OWN OVERCLAIM",
         "evidence": "ker(J_keep) is 5-dimensional and the held-out row is "
                     "reachable from it in every case; a kernel shift fits the "
                     "withheld mass to ~1e-15 at "
                     f"{holdout['at_norm_cost_between'][0]:.2f}-"
                     f"{holdout['at_norm_cost_between'][1]:.2f}x the norm. It "
                     "measures the regulariser, not the model"},
        {"claim": "the quark model is 'sloppy' in the Sethna sense",
         "verdict": "NOT THE RIGHT DIAGNOSIS",
         "evidence": f"the four nonzero singular values span only "
                     f"{rank['nonzero_singular_value_ratio']:.1f}x -- no long "
                     "hierarchy -- but the rank deficiency is exact; the "
                     "pathology is structural non-identifiability, not "
                     "ill-conditioning"},
        {"claim": "N = 466 drifting is a defect of N, and any knob would drift",
         "verdict": "REFRAMED, AND THE SECOND HALF IS FALSE",
         "evidence": f"identifiable share runs from "
                     f"{rank['identifiable_share_by_knob']['uplift_asymmetry']:.4f} "
                     f"(uplift_asymmetry) to "
                     f"{rank['identifiable_share_by_knob']['eta_k3k5_minus']:.4f} "
                     "(eta_k3k5_minus); which knobs drift is itself structure"},
        {"claim": "the missing correlation is a scalar relation R = f(p, T)",
         "verdict": "REFUTED",
         "evidence": "the degeneracy is a linear subspace selected by the "
                     "response map; the nearest pair is gamma_q/transport at "
                     "178.9 deg, not transport/resistance"},
    ]
    return {
        "entries": entries,
        "headline": (
            "Four scored masses cannot identify the current quark "
            "parameterisation. The positive-knob Jacobian is numerically full "
            "row rank, and the radial residual triple's individual proximity "
            "to fitted knobs is not enough to infer its effect on the spectrum"),
        "what_would_settle_it": (
            "More observables, not more knobs. The v4 flavor-CP layer supplies "
            "CKM angles and J from the same Hamiltonian and inherits the v3 "
            "eigenvalues, so it adds observable rows -- but it is NOT "
            "automatically column-free: v4 QuarkParams carries phi_h, targeted "
            "eta couplings and per-shell diagonal shifts. A joint "
            "identifiability audit has to decide which of those are externally "
            "derived and which belong in the response columns"),
        "metric_scope": rank["metric_scope"],
    }
