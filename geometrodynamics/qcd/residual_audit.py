"""The downstream re-derivation PR #271 deferred: the quark residual sector.

THE RESULT, FIRST
─────────────────
PR #271 corrected the 5D radial scalar operator and reopened the *lepton*
sector's geometric closure — ``γ = 22.5`` survived as a requirement of the
locked surrogate but lost its derivation from the canonical ``R_OUTER = 1.26``
barrier geometry.  The obvious expectation was that the quark sector, whose
three residual knobs are derived from the *same* eigensolver, would break the
same way.

**It does not. All three quark residuals move toward their locked values.**

    residual                     locked      legacy         corrected
    pinhole   Σ V_max[1..5]      22.25       22.008 (−1.09%)  22.331 (+0.36%)
    transport ⟨u₁|ΔV|u₂⟩ mean     0.54        0.5447 (+0.88%)   0.5438 (+0.70%)
    resistance transport·ln(α₅/α₁) 0.14       0.1407 (+0.49%)   0.1400 (−0.02%)

The two sectors were reading the same barrier and disagreeing about it.  Under
the legacy operator the lepton target (22.5) and the quark target (22.25) both
sat *above* the computed sum (22.008); the correction raises the sum by 1.47 %
and carries it past the quark target while leaving it short of the lepton one.

Why the sectors part company: elasticity
────────────────────────────────────────
Neither sector's agreement is worth what its percentage suggests, and they are
not worth the same amount.  At their respective locks

    lepton    d ln m_μ / d ln γ        = −17.471
    quark     d ln m_s / d ln pinhole  =  +4.798

so the quark ladder is ~3.6× less stiff.  The lepton's improved −0.75 % residual
is still a 9 % muon error; the quark's +0.36 % is a ~1.7 % strange error.  Same
operator, same barrier, different consequence — because the two Hamiltonians
amplify their geometric scalar by different factors.

**A percentage agreement between a geometric quantity and a fitted knob is
meaningless until multiplied by the elasticity of the thing it feeds.**  That is
the transferable lesson of this round, and it is why #271's headline and this
round's headline point in opposite directions without contradicting each other.

What the ladder actually wants
──────────────────────────────
Scanning the locked quark Hamiltonian over the pinhole alone, the ladder is
minimised at ``22.228`` (max rel err 1.348 %), not at the fitted ``22.25``
(1.610 %).  Measured against *that* — the ladder's own optimum, not a fitted
knob — the correction roughly halves the miss:

    legacy    22.008  →  −0.989 %   (ladder error 3.613 %)
    corrected 22.331  →  +0.464 %   (ladder error 3.397 %)

The composite is not the sum of its parts
─────────────────────────────────────────
Substituting all three derived residuals at once, with **no retuning**, the
corrected set scores *worse* than the legacy set (3.78 % vs 3.44 %) even though
each individual residual is closer to its lock.  The legacy triple enjoys a
partial cancellation: its pinhole error (3.61 % alone) is walked back to 3.44 %
by its transport and resistance errors pushing the other way.  The corrected
triple has no such luck — 3.40 % alone becomes 3.78 % together.

Neither composite reaches the fitted lock's 1.61 %.  **"Each knob is derived to
within 1 %" and "the derived knobs reproduce the ladder" are different claims,
and only the first was ever established.**  That is true of the legacy numbers
too; the correction did not create this gap, it exposed it.

The one thing that genuinely improves: cross-sector R_OUTER
───────────────────────────────────────────────────────────
Each sector, read as a demand on ``R_OUTER``, gives

                    legacy               corrected
    lepton   γ=22.50   1.28737 (+2.17 %)   1.26788 (+0.63 %)
    quark    pin=22.25  1.27229 (+0.98 %)   1.25645 (−0.28 %)
    split               1.179 %             0.906 %

Under the legacy operator both sectors demanded ``R_OUTER > 1.26``, putting the
canonical value **outside** the bracket they define — 0.81 bracket-widths below
it.  Under the corrected operator they **straddle** it: the quark sector wants
less, the lepton sector wants more, and 1.26 sits at 31 % of the way across.

This is evidence *for* ``R_OUTER = 1.26``, and it is not the evidence #271
reopened — that was a single-sector fixed point, this is a two-sector bracket.
It is also weak: two numbers 0.9 % apart bracket any value between them, and
nothing here derives 1.26 rather than admitting it.  Recorded as suggestive,
not as a restored derivation.

What is deliberately **not** re-run
───────────────────────────────────
The same exclusion as #271, for the same reason: the four quark shell-index
axioms (ε, η, χ, phase) are expressible in ``k₅ = 5`` alone, the Hopf transport
phase is topological, and the CKM/flavor-CP layer inherits the v3 eigenvalues by
construction.  None of them reads the radial scalar operator.  Proximity is not
dependence.
"""

from __future__ import annotations

import math
from dataclasses import replace
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
from scipy.optimize import minimize_scalar

from geometrodynamics.constants import R_INNER, R_MID, R_OUTER
from geometrodynamics.qcd.quark_spectrum import (LOCKED_QUARK_PARAMS,
                                                 OBSERVED_MASSES_MEV,
                                                 QUARK_SPECIES,
                                                 extract_physical_spectrum)
from geometrodynamics.tangherlini.alpha_q import throat_du_dr
from geometrodynamics.tangherlini.lepton_spectrum import (
    LEPTON_BASELINE_PINHOLE, LEPTON_BASELINE_RESISTANCE)
from geometrodynamics.tangherlini.operator_audit import (OPERATORS,
                                                         chebyshev_grid,
                                                         eigen_solve,
                                                         r_outer_fixed_point,
                                                         vmax_sum)
from geometrodynamics.tangherlini.radial import r_to_rstar

__all__ = [
    "LOCKED_PINHOLE",
    "LOCKED_TRANSPORT",
    "LOCKED_RESISTANCE",
    "FITTED_SPECIES",
    "TRANSPORT_PAIRS",
    "mode_profiles",
    "derive_pinhole",
    "derive_transport",
    "derive_alpha_q_table",
    "derive_resistance",
    "derive_the_three_residuals",
    "ladder_error",
    "measure_the_three_residuals_under_both_operators",
    "measure_the_ladder_without_retuning",
    "measure_the_pinhole_elasticity",
    "measure_the_cross_sector_r_outer_bracket",
    "measure_the_lepton_resistance_selector",
    "measure_the_n_stability",
    "measure_the_two_sectors_side_by_side",
    "measure_the_quark_ledger",
]

# The three residual knobs of the frozen v3 quark lock (quark_spectrum.py).
LOCKED_PINHOLE = 22.25
LOCKED_TRANSPORT = 0.54
LOCKED_RESISTANCE = 0.14

# ``u`` is zero by construction under ``min_eigenvalue``; ``d`` is the MeV
# anchor.  The quoted "1.6 % max rel err" is over the remaining four.
FITTED_SPECIES: Tuple[str, ...] = ("s", "c", "b", "t")

# The same-partition off-diagonals the model actually carries.
TRANSPORT_PAIRS: Tuple[Tuple[int, int], ...] = ((1, 3), (3, 5), (1, 5))

_LMAX = 5


# ── the three derivations, with the operator injectable ─────────────────────

def mode_profiles(potential: Callable, points: int = 80, lmax: int = _LMAX
                  ) -> Tuple[Dict[int, dict], np.ndarray]:
    """``solve_radial_modes``'s ``funcs`` dict, for whichever ``V`` is given.

    Reproduced from the production solver rather than imported, because
    ``solve_radial_modes`` hard-wires the canonical ``V_tangherlini``; the audit
    has to run both operators through numerics that are otherwise identical.
    """
    grid, _, _ = chebyshev_grid(points, R_MID, R_OUTER)
    out: Dict[int, dict] = {}
    for ell in range(1, lmax + 1):
        omegas, vecs, _ = eigen_solve(ell, potential, points=points, n_modes=4)
        funcs = []
        for k in range(vecs.shape[1]):
            u = np.zeros(points + 1)
            u[1:points] = vecs[:, k]
            if abs(u.min()) > u.max():
                u = -u
            u /= abs(u).max() + 1e-12
            funcs.append({"u_half": u, "r_phys": grid})
        out[ell] = {"omega": omegas, "funcs": funcs}
    return out, grid


def derive_pinhole(potential: Callable, points: int = 80,
                   r_outer: float = R_OUTER) -> float:
    """``Σ_{ℓ=1..5} max_r V(r, ℓ)`` on the eigensolver's tortoise grid.

    The tortoise grid is the canonical evaluation domain — the same one that
    defines the bound modes.  A raw-``r`` linspace gives a different number
    (see `derive_pinhole_raw_grid`); the difference is quadrature, not physics,
    and the tortoise value is the one the README quotes.
    """
    return vmax_sum(potential, 1, _LMAX, points=points, r_outer=r_outer)


def derive_pinhole_raw_grid(potential: Callable, points: int = 400) -> float:
    """The same sum on a uniform ``r`` linspace — kept because
    ``scripts/experiment_residuals_from_geometry.py`` uses it and reports a
    visibly different number."""
    grid = np.linspace(R_INNER + 0.01, R_OUTER - 0.01, points)
    return float(sum(float(np.max(np.asarray(potential(grid, ell, R_MID),
                                             dtype=float)))
                     for ell in range(1, _LMAX + 1)))


def derive_transport(potential: Callable, points: int = 80
                     ) -> Tuple[float, Dict[Tuple[int, int], float]]:
    """``mean_pairs |⟨u_{ℓ₁}| V_{ℓ₂} − V_{ℓ₁} |u_{ℓ₂}⟩|``, L2-normalised in ``r*``.

    The genuine off-diagonal element — both eigenvectors, not one squared.  The
    operator ``V_{ℓ₂} − V_{ℓ₁}`` is *exactly* invariant under the correction
    (``ΔV = 3A²/4r²`` carries no ``ℓ``), so every part of the drift here comes
    from the eigenvectors.
    """
    profiles, grid = mode_profiles(potential, points)
    rstar = np.array([r_to_rstar(float(v), R_MID) for v in grid])
    order = np.argsort(rstar)
    xs = rstar[order]

    normed: Dict[int, np.ndarray] = {}
    for ell in range(1, _LMAX + 1):
        u = profiles[ell]["funcs"][0]["u_half"][order]
        normed[ell] = u / math.sqrt(float(np.trapezoid(u * u, xs)))

    per_pair: Dict[Tuple[int, int], float] = {}
    for (l1, l2) in TRANSPORT_PAIRS:
        delta = (np.asarray(potential(grid, l2, R_MID), dtype=float)
                 - np.asarray(potential(grid, l1, R_MID), dtype=float))[order]
        per_pair[(l1, l2)] = abs(float(
            np.trapezoid(normed[l1] * normed[l2] * delta, xs)))
    return float(np.mean(list(per_pair.values()))), per_pair


def derive_alpha_q_table(potential: Callable, points: int = 80
                         ) -> Dict[int, float]:
    """``α_q(ℓ,0) = (du/dr)|_throat / |(du/dr)|_throat at ℓ=1|`` for ℓ = 1..5."""
    profiles, _ = mode_profiles(potential, points)
    ref = abs(throat_du_dr(profiles[1]["funcs"][0], R_MID))
    return {ell: float(throat_du_dr(profiles[ell]["funcs"][0], R_MID) / ref)
            for ell in range(1, _LMAX + 1)}


def derive_resistance(potential: Callable, points: int = 80,
                      transport: Optional[float] = None) -> float:
    """``transport · ln(α_q(5,0) / α_q(1,0))``.

    ``transport`` defaults to the value derived from the *same* operator.  The
    README's ``−0.43 %`` was computed with the **locked** ``0.54`` instead; both
    conventions are reported by `measure_the_three_residuals_under_both_operators`
    because they disagree in sign about whether the legacy value overshoots.
    """
    table = derive_alpha_q_table(potential, points)
    if transport is None:
        transport, _ = derive_transport(potential, points)
    return float(transport * math.log(table[5] / table[1]))


def derive_the_three_residuals(potential: Callable, points: int = 80) -> dict:
    """All three residual knobs from one operator, in one pass."""
    transport, per_pair = derive_transport(potential, points)
    table = derive_alpha_q_table(potential, points)
    log_ratio = math.log(table[5] / table[1])
    return {
        "pinhole": derive_pinhole(potential, points),
        "pinhole_raw_grid": derive_pinhole_raw_grid(potential),
        "transport": transport,
        "transport_pairs": {f"{a},{b}": v for (a, b), v in per_pair.items()},
        "alpha_q": table,
        "alpha_q_log_ratio": log_ratio,
        "resistance": float(transport * log_ratio),
        "resistance_at_locked_transport": float(LOCKED_TRANSPORT * log_ratio),
    }


# ── passing them through the locked ladder ──────────────────────────────────

def ladder_error(**overrides) -> Tuple[Dict[str, float], Dict[str, float], float]:
    """The locked quark Hamiltonian, ``d``-anchored, with **no retuning**.

    Returns ``(predicted_mev, signed_relative_error, max_abs_relative_error)``
    over `FITTED_SPECIES`.
    """
    params = replace(LOCKED_QUARK_PARAMS, **overrides) if overrides \
        else LOCKED_QUARK_PARAMS
    spectrum = extract_physical_spectrum(params)
    scale = OBSERVED_MASSES_MEV["d"] / spectrum["d"]
    predicted = {s: float(spectrum[s] * scale) for s in QUARK_SPECIES}
    signed = {s: float((predicted[s] - OBSERVED_MASSES_MEV[s])
                       / OBSERVED_MASSES_MEV[s]) for s in FITTED_SPECIES}
    return predicted, signed, max(abs(v) for v in signed.values())


# ── the measurements ────────────────────────────────────────────────────────

def measure_the_three_residuals_under_both_operators() -> dict:
    """Q1 — the headline table: every residual moves toward its locked value."""
    derived = {name: derive_the_three_residuals(V)
               for name, V in OPERATORS.items()}
    locked = {"pinhole": LOCKED_PINHOLE, "transport": LOCKED_TRANSPORT,
              "resistance": LOCKED_RESISTANCE}
    rows = []
    for knob, lock in locked.items():
        legacy = derived["legacy"][knob]
        corrected = derived["scalar_correct"][knob]
        rows.append({
            "residual": knob,
            "locked": lock,
            "legacy": legacy,
            "corrected": corrected,
            "legacy_rel_pct": 100.0 * (legacy - lock) / lock,
            "corrected_rel_pct": 100.0 * (corrected - lock) / lock,
            "drift_pct": 100.0 * (corrected - legacy) / legacy,
            "moves_toward_lock": abs(corrected - lock) < abs(legacy - lock),
        })
    return {
        "rows": rows,
        "derived": derived,
        "all_three_move_toward_the_lock": all(r["moves_toward_lock"]
                                              for r in rows),
        "note": ("the README's resistance residual (-0.43%) used the LOCKED "
                 "transport 0.54, not the derived one; both conventions are "
                 "carried in `derived[*]['resistance*']`"),
    }


def measure_the_ladder_without_retuning() -> dict:
    """Q2 — substitute the derived residuals into the locked ladder, retune nothing.

    The composite ordering **reverses** the individual ordering: each corrected
    residual is closer to its lock, yet the corrected triple scores worse,
    because the legacy triple's errors partially cancel.
    """
    baseline = ladder_error()
    derived = {name: derive_the_three_residuals(V)
               for name, V in OPERATORS.items()}
    knobs = ("pinhole", "transport", "resistance")

    composite, one_at_a_time = {}, {}
    for name, vals in derived.items():
        composite[name] = ladder_error(**{k: vals[k] for k in knobs})[2]
        one_at_a_time[name] = {k: ladder_error(**{k: vals[k]})[2] for k in knobs}

    return {
        "locked_max_rel_err": baseline[2],
        "locked_signed": baseline[1],
        "composite_max_rel_err": composite,
        "one_at_a_time_max_rel_err": one_at_a_time,
        "corrected_composite_is_worse": (composite["scalar_correct"]
                                         > composite["legacy"]),
        "corrected_pinhole_alone_is_better": (
            one_at_a_time["scalar_correct"]["pinhole"]
            < one_at_a_time["legacy"]["pinhole"]),
        "legacy_composite_beats_its_own_pinhole": (
            composite["legacy"] < one_at_a_time["legacy"]["pinhole"]),
        "neither_composite_reaches_the_lock": all(
            v > baseline[2] for v in composite.values()),
    }


def measure_the_pinhole_elasticity() -> dict:
    """Q3 — how much mass error one percent of pinhole error buys.

    The quark analogue of #271's ``d ln m_μ/d ln γ = −17.471``.  Reported as a
    local central derivative at the lock *and* as a secant across the two
    derived pinholes — they are different quantities and #271's review was
    right to insist they be labelled separately.
    """
    def masses(pinhole: float) -> Dict[str, float]:
        return ladder_error(pinhole=pinhole)[0]

    def worst(pinhole: float) -> float:
        return ladder_error(pinhole=pinhole)[2]

    eps = 1e-3 * LOCKED_PINHOLE
    hi, lo = masses(LOCKED_PINHOLE + eps), masses(LOCKED_PINHOLE - eps)
    denom = math.log(LOCKED_PINHOLE + eps) - math.log(LOCKED_PINHOLE - eps)
    local = {s: float((math.log(hi[s]) - math.log(lo[s])) / denom)
             for s in FITTED_SPECIES}

    legacy = derive_pinhole(OPERATORS["legacy"])
    corrected = derive_pinhole(OPERATORS["scalar_correct"])
    ma, mb = masses(legacy), masses(corrected)
    span = math.log(corrected) - math.log(legacy)
    secant = {s: float((math.log(mb[s]) - math.log(ma[s])) / span)
              for s in FITTED_SPECIES}

    opt = minimize_scalar(worst, bracket=(21.9, LOCKED_PINHOLE, 22.5),
                          method="brent", options={"xtol": 1e-8})
    wanted = float(opt.x)
    return {
        "local_d_ln_m_d_ln_pinhole_at_the_lock": local,
        "secant_across_the_two_derived_pinholes": secant,
        "stiffest_species": max(local, key=lambda s: abs(local[s])),
        "pinhole_the_ladder_wants": wanted,
        "ladder_error_at_what_it_wants": float(opt.fun),
        "ladder_error_at_the_v3_lock": worst(LOCKED_PINHOLE),
        "distance_from_what_it_wants_pct": {
            "v3_lock": 100.0 * (LOCKED_PINHOLE - wanted) / wanted,
            "legacy": 100.0 * (legacy - wanted) / wanted,
            "corrected": 100.0 * (corrected - wanted) / wanted,
        },
        "correction_halves_the_miss": (abs(corrected - wanted)
                                       < abs(legacy - wanted)),
    }


def measure_the_cross_sector_r_outer_bracket() -> dict:
    """Q4 — read as demands on ``R_OUTER``, do the two sectors bracket 1.26?

    Legacy: both demand more than 1.26, so the canonical value falls outside.
    Corrected: they straddle it.  This is genuinely new evidence for 1.26 and
    it is *not* the single-sector fixed point #271 reopened — but two numbers
    0.9 % apart bracket anything between them, so it is suggestive only.
    """
    out: Dict[str, dict] = {}
    for name, V in OPERATORS.items():
        lepton = r_outer_fixed_point(V, 1, _LMAX, target=LEPTON_BASELINE_PINHOLE,
                                     bracket=(1.05, 1.70))
        quark = r_outer_fixed_point(V, 1, _LMAX, target=LOCKED_PINHOLE,
                                    bracket=(1.05, 1.70))
        lo, hi = min(lepton, quark), max(lepton, quark)
        width = hi - lo
        out[name] = {
            "lepton_r_outer": lepton,
            "quark_r_outer": quark,
            "lepton_vs_canonical_pct": 100.0 * (lepton - R_OUTER) / R_OUTER,
            "quark_vs_canonical_pct": 100.0 * (quark - R_OUTER) / R_OUTER,
            "split_pct_of_mean": 100.0 * width / (0.5 * (lepton + quark)),
            "brackets_canonical": bool(lo <= R_OUTER <= hi),
            "canonical_position_in_bracket": float((R_OUTER - lo) / width),
        }
    return {
        "per_operator": out,
        "correction_narrows_the_split": (
            out["scalar_correct"]["split_pct_of_mean"]
            < out["legacy"]["split_pct_of_mean"]),
        "only_the_corrected_operator_brackets_1_26": (
            out["scalar_correct"]["brackets_canonical"]
            and not out["legacy"]["brackets_canonical"]),
        "caveat": ("a 0.9% bracket admits any value inside it; this supports "
                   "R_OUTER = 1.26 but does not derive it"),
    }


def measure_the_lepton_resistance_selector() -> dict:
    """Q5 — the one lepton-sector claim that still reads a radial eigenvalue.

    README: *"Resistance = 7π/100 — selected over 4·(ω−1) by R_OUTER bisection."*
    The selector is #271-reopened, but under the corrected operator the rejected
    competitor degrades from ``+0.48 %`` to ``+2.50 %``, so ``7π/100`` now wins
    on proximity as well.  **Conclusion survives; its stated reason does not.**
    """
    closed_form = 7.0 * math.pi / 100.0
    rows = {}
    for name, V in OPERATORS.items():
        omega = float(eigen_solve(1, V, n_modes=1)[0][0])
        competitor = 4.0 * (omega - 1.0)
        rows[name] = {
            "omega_1_0": omega,
            "four_omega_minus_one": competitor,
            "competitor_rel_pct": 100.0 * (competitor - LEPTON_BASELINE_RESISTANCE)
            / LEPTON_BASELINE_RESISTANCE,
        }
    return {
        "locked_resistance": LEPTON_BASELINE_RESISTANCE,
        "seven_pi_over_100": closed_form,
        "closed_form_rel_pct": 100.0 * (closed_form - LEPTON_BASELINE_RESISTANCE)
        / LEPTON_BASELINE_RESISTANCE,
        "per_operator": rows,
        "competitor_beat_the_closed_form_under_legacy": (
            abs(rows["legacy"]["competitor_rel_pct"])
            < abs(100.0 * (closed_form - LEPTON_BASELINE_RESISTANCE)
                  / LEPTON_BASELINE_RESISTANCE)),
        "closed_form_wins_under_the_correction": (
            abs(100.0 * (closed_form - LEPTON_BASELINE_RESISTANCE)
                / LEPTON_BASELINE_RESISTANCE)
            < abs(rows["scalar_correct"]["competitor_rel_pct"])),
        "verdict": ("conclusion survives, stated reason does not: the "
                    "R_OUTER-bisection selector is reopened by PR #271, and "
                    "under the legacy operator it had selected the WORSE-fitting "
                    "of the two candidates"),
    }


def _coordinate_descent(base, err: Callable, axes: List[Tuple], rounds: int = 6):
    current, current_err = base, err(base)
    for _ in range(rounds):
        moved = False
        for values, apply in axes:
            best_v, best_e = None, current_err
            for v in values:
                candidate = err(apply(current, float(v)))
                if candidate < best_e:
                    best_e, best_v = candidate, float(v)
            if best_v is not None:
                current, current_err, moved = apply(current, best_v), best_e, True
        if not moved:
            break
    return current, current_err


def measure_the_n_stability(rounds: int = 6) -> dict:
    """Q6 — does pinning the residuals to *corrected* geometry stabilise ``N``?

    It does not, under either operator.  The v3 lock's own conclusion — ``β``
    is the model's last fit knob and ``N = 466`` is a compensator, not a
    topological invariant — is unchanged by the correction.  Only the digits
    move, including ``N`` at PDG masses itself (466 → 458).
    """
    beta_grid = np.arange(380, 561, 2)
    grids = {"transport": np.arange(0.40, 0.80, 0.01),
             "resistance": np.arange(0.06, 0.22, 0.005),
             "pinhole": np.arange(18.0, 27.0, 0.25)}
    perturbations = {
        "PDG": OBSERVED_MASSES_MEV,
        "c_x_1.10": {**OBSERVED_MASSES_MEV, "c": OBSERVED_MASSES_MEV["c"] * 1.10},
        "b_x_1.10": {**OBSERVED_MASSES_MEV, "b": OBSERVED_MASSES_MEV["b"] * 1.10},
        "t_x_1.10": {**OBSERVED_MASSES_MEV, "t": OBSERVED_MASSES_MEV["t"] * 1.10},
    }

    def make_err(observed):
        target = np.array([observed[s] for s in QUARK_SPECIES], dtype=float)
        skip = QUARK_SPECIES.index("u")

        def err(params) -> float:
            try:
                spectrum = extract_physical_spectrum(params)
            except Exception:
                return float("inf")
            if spectrum["d"] <= 1e-6:
                return float("inf")
            pred = np.array([spectrum[s] * observed["d"] / spectrum["d"]
                             for s in QUARK_SPECIES], dtype=float)
            rel = np.abs(pred - target) / target
            return float(np.max([rel[i] for i in range(len(QUARK_SPECIES))
                                 if i != skip]))
        return err

    derived = {name: derive_the_three_residuals(V)
               for name, V in OPERATORS.items()}
    modes = {"baseline": None, "legacy": derived["legacy"],
             "corrected": derived["scalar_correct"]}

    out: Dict[str, dict] = {}
    for mode, vals in modes.items():
        if vals is None:
            base = LOCKED_QUARK_PARAMS
            axes = [(beta_grid, lambda p, v: replace(p, beta=v * math.pi / 2.0))]
            axes += [(grids[k], (lambda key: lambda p, v: replace(p, **{key: v}))(k))
                     for k in grids]
        else:
            base = replace(LOCKED_QUARK_PARAMS,
                           **{k: vals[k] for k in grids})
            axes = [(beta_grid, lambda p, v: replace(p, beta=v * math.pi / 2.0))]
        found = {}
        for name, observed in perturbations.items():
            lock, _ = _coordinate_descent(base, make_err(observed), axes, rounds)
            found[name] = int(round(lock.beta * 2.0 / math.pi))
        span = max(found.values()) - min(found.values())
        out[mode] = {"N": found, "range": [min(found.values()),
                                           max(found.values())],
                     "width": span}
    return {
        "per_mode": out,
        "N_drifts_under_every_mode": all(v["width"] > 0 for v in out.values()),
        "verdict": ("unchanged by the correction: beta remains the model's "
                    "last fit knob"),
    }


def measure_the_two_sectors_side_by_side() -> dict:
    """Q7 — the same barrier, the same correction, two opposite verdicts.

    Both sectors are put in the same frame: take the residual each one has
    against its own locked scalar under the corrected operator, and report the
    **measured** heaviest-species error that residual actually produces — not
    the linearisation.

    This also corrects a number PR #271 left in the README.  That row read
    *"the residual improves to −0.75 %, but d ln m_μ/d ln γ = −17.5 at the lock
    makes that a 9 % muon error"*.  The elasticity is right and the residual is
    right, but ``9 %`` follows from neither: linearising gives ``+14.0 %`` and
    the locked block actually returns ``+15.2 %`` at ``γ = 22.331``.  The audit's
    own A/B/C table carried the correct value the whole time.
    """
    from geometrodynamics.tangherlini.operator_audit import \
        measure_which_geometry_preserves_the_lepton_ladder

    lepton = measure_which_geometry_preserves_the_lepton_ladder()
    corrected_at_canonical = next(
        c for c in lepton["rows"] if c["case"] == "A corrected R=1.26, gamma[1..5]")
    lepton_gamma = float(corrected_at_canonical["gamma"])
    lepton_residual = 100.0 * (lepton_gamma - LEPTON_BASELINE_PINHOLE) \
        / LEPTON_BASELINE_PINHOLE
    lepton_elasticity = float(
        lepton["local_d_ln_m_mu_over_d_ln_gamma_at_22p5"])
    lepton_measured = float(corrected_at_canonical["mu_error_percent"])
    lepton_linear = 100.0 * (math.exp(lepton_elasticity * lepton_residual / 100.0)
                             - 1.0)

    quark_pinhole = derive_pinhole(OPERATORS["scalar_correct"])
    quark_residual = 100.0 * (quark_pinhole - LOCKED_PINHOLE) / LOCKED_PINHOLE
    elasticity = measure_the_pinhole_elasticity()
    quark_elasticity = elasticity["local_d_ln_m_d_ln_pinhole_at_the_lock"]["s"]
    at_lock = ladder_error()[1]["s"]
    at_derived = ladder_error(pinhole=quark_pinhole)[1]["s"]
    quark_measured = 100.0 * (at_derived - at_lock)
    quark_linear = 100.0 * (math.exp(quark_elasticity * quark_residual / 100.0)
                            - 1.0)

    return {
        "lepton": {
            "scalar": "gamma = Sum V_max[1..5] vs the locked 22.5",
            "corrected_value": lepton_gamma,
            "residual_pct": lepton_residual,
            "elasticity": lepton_elasticity,
            "linearised_error_pct": lepton_linear,
            "measured_error_pct": lepton_measured,
            "species": "mu",
        },
        "quark": {
            "scalar": "pinhole = Sum V_max[1..5] vs the fitted 22.25",
            "corrected_value": quark_pinhole,
            "residual_pct": quark_residual,
            "elasticity": quark_elasticity,
            "linearised_error_pct": quark_linear,
            "measured_error_pct": quark_measured,
            "species": "s",
        },
        "elasticity_ratio": abs(lepton_elasticity / quark_elasticity),
        "the_readme_9_percent_is_wrong": True,
        "readme_said_pct": 9.0,
        "corrected_to_pct": lepton_measured,
        "why_the_verdicts_differ": (
            "the same barrier feeds both sectors and the correction moves it "
            "the same way; the lepton Hamiltonian amplifies its residual by "
            f"{abs(lepton_elasticity):.1f} and the quark ladder by "
            f"{abs(quark_elasticity):.1f}, so a comparable geometric agreement "
            "is a 15% miss in one sector and a 2% miss in the other"),
    }


def measure_the_quark_ledger() -> dict:
    """Q7 — the verdict table, in #271's three categories."""
    residuals = measure_the_three_residuals_under_both_operators()
    ladder = measure_the_ladder_without_retuning()
    elasticity = measure_the_pinhole_elasticity()
    bracket = measure_the_cross_sector_r_outer_bracket()
    selector = measure_the_lepton_resistance_selector()

    def pct(knob: str, key: str) -> float:
        return next(r[key] for r in residuals["rows"] if r["residual"] == knob)

    entries = [
        {"claim": "quark shell-index axioms (eps, eta, chi, phase) in k_5 = 5",
         "verdict": "EXACTLY INVARIANT",
         "evidence": "expressible in k_5 alone; reads no radial operator"},
        {"claim": "transport operator V_l2 - V_l1 as an algebraic object",
         "verdict": "EXACTLY INVARIANT",
         "evidence": "the correction 3A^2/4r^2 carries no l"},
        {"claim": "pinhole = Sum V_max[1..5] vs the fitted 22.25",
         "verdict": "NUMERICALLY SHIFTED, AND IMPROVED",
         "evidence": f"{pct('pinhole', 'legacy_rel_pct'):+.2f}% -> "
                     f"{pct('pinhole', 'corrected_rel_pct'):+.2f}%"},
        {"claim": "transport = mean off-diagonal vs the fitted 0.54",
         "verdict": "NUMERICALLY SHIFTED, AND IMPROVED",
         "evidence": f"{pct('transport', 'legacy_rel_pct'):+.2f}% -> "
                     f"{pct('transport', 'corrected_rel_pct'):+.2f}%"},
        {"claim": "resistance = transport * ln(alpha_5/alpha_1) vs the fitted 0.14",
         "verdict": "NUMERICALLY SHIFTED, AND IMPROVED",
         "evidence": f"{pct('resistance', 'legacy_rel_pct'):+.2f}% -> "
                     f"{pct('resistance', 'corrected_rel_pct'):+.2f}%"},
        {"claim": "the derived residuals reproduce the quark ladder",
         "verdict": "INTERPRETATION CHANGED",
         "evidence": "never established under either operator: composite "
                     f"{100*ladder['composite_max_rel_err']['legacy']:.2f}% "
                     f"(legacy) and "
                     f"{100*ladder['composite_max_rel_err']['scalar_correct']:.2f}% "
                     f"(corrected) against the fitted lock's "
                     f"{100*ladder['locked_max_rel_err']:.2f}%; per-knob "
                     "agreement is not ladder agreement"},
        {"claim": "N = 466 as a stable integer of the residual-pinned fit",
         "verdict": "NUMERICALLY SHIFTED",
         "evidence": "still a compensator; N at PDG masses moves under the "
                     "correction and the drift width stays O(50)"},
        {"claim": "lepton resistance = 7pi/100 selected over 4(omega-1)",
         "verdict": "INTERPRETATION CHANGED",
         "evidence": "conclusion survives on proximity "
                     f"({selector['closed_form_rel_pct']:+.2f}% vs "
                     f"{selector['per_operator']['scalar_correct']['competitor_rel_pct']:+.2f}%) "
                     "but the R_OUTER-bisection selector that chose it is "
                     "reopened by PR #271"},
        {"claim": "R_OUTER = 1.26 from a single-sector fixed point",
         "verdict": "INTERPRETATION CHANGED",
         "evidence": "reopened by PR #271 and not restored here"},
        {"claim": "R_OUTER = 1.26 bracketed by two independent sectors",
         "verdict": "NEW, AND WEAK",
         "evidence": "legacy puts 1.26 outside the lepton/quark bracket; "
                     "corrected straddles it at "
                     f"{bracket['per_operator']['scalar_correct']['canonical_position_in_bracket']:.2f} "
                     "of the way across a "
                     f"{bracket['per_operator']['scalar_correct']['split_pct_of_mean']:.2f}% "
                     "window — suggestive, not a derivation"},
    ]
    return {
        "entries": entries,
        "headline": ("the quark residual sector survives the operator "
                     "correction and improves; the lepton sector's closure did "
                     "not, and the difference is elasticity: "
                     f"d ln m_s/d ln pinhole = "
                     f"{elasticity['local_d_ln_m_d_ln_pinhole_at_the_lock']['s']:+.3f} "
                     "against the lepton's -17.471"),
    }
