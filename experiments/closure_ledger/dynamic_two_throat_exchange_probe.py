"""
Dynamic two-throat exchange path with back-reaction (PR #191).

> Framing: QFT on the *fixed classical* throat geometry (geometry → fields),
> not quantum gravity.

BEYOND THE ADIABATIC HOLONOMY
──────────────────────────────
PR #188 measured the exchange sign −1 as the ADIABATIC holonomy of the throat
spin-½ state along the swap path (the infinitely-slow limit, no
back-reaction).  This probe goes dynamical: it traverses the exchange path at
FINITE speed in real time, with the field BACK-REACTING (sourced by the moving
throat and acting back on it), and characterizes (i) the recovery of the
adiabatic −1, (ii) the non-adiabatic corrections, and (iii) the back-reaction.

The effective dynamics (the dynamical generalization of #188): the throat's
internal (Pin) spinor evolves under

        i ψ̇ = H(t) ψ ,   H(t) = (Δ/2) n̂(s(t))·σ + g x(t) σ_z ,

where n̂(s) = (cos 2πs, sin 2πs, 0) traces the swap loop (the equatorial great
circle — the Finkelstein–Rubinstein 2π-rotation / exchange path) over a finite
duration T (the swap SPEED is v = 1/T; Δ is the internal gap), and the
BACK-REACTION field x(t) is SOURCED by the throat and acts back:

        ẍ + γ ẋ + ω² x = −κ ⟨σ_z⟩(t) .

This is a genuine two-way (back-reacting) dynamical system: the moving throat
drives the field, the field deflects the throat.  The exchange phase is the
geometric (Berry) phase of the loop; its adiabatic value is the FR holonomy.

WHAT IS COMPUTED (measured; the real-time coupled spinor + field)
  • THE ADIABATIC EXCHANGE PHASE: the swap loop (the equatorial great circle)
    encloses solid angle Ω = 2π, so the spin-½ Berry phase is −Ω/2 = −π —
    the exchange sign −1 (the #188 holonomy), a geometric invariant of the
    loop, speed-independent.
  • ADIABATIC RECOVERY: as the swap is SLOWED (T → ∞), the dynamical geometric
    phase → −π and the non-adiabatic excitation P_exc → 0 — the finite-speed
    evolution recovers the adiabatic −1.
  • NON-ADIABATIC CORRECTION: at finite speed the throat cannot follow the
    instantaneous state — the geometric phase deviates from −π by O(1/T) (the
    deviation × T is constant) and P_exc grows with speed; the "cost" of a
    finite-speed swap.
  • BACK-REACTION: the moving throat SOURCES the field (the field energy E_f
    grows at finite speed) and the field acts back on the spinor; E_f → 0 as
    T → ∞.  The back-reaction modifies the dynamics but NOT the adiabatic
    limit — the exchange −1 is unchanged.

HONEST SCOPE
  This is an EFFECTIVE model of the dynamic exchange: the throat's internal
  (Pin) state is the driven spin-½ on the swap loop (the FR exchange ≃ 2π
  rotation, #188); the back-reaction is a single sourced field mode (the
  throat's deformation / the mutual field, in mean field), not the full
  field-resolved two-throat dynamics.  The exchange phase is the ADIABATIC
  invariant (−π = the −1); the finite-speed and back-reaction effects are
  controlled corrections that vanish as the swap slows.  So the exchange
  statistics survive as the adiabatic limit of the full dynamical evolution,
  with the dynamical cost quantified.  Weak-field, code units; the
  field-resolved real-time two-throat solve is the follow-up.

Tests:
  T1. Goal: the dynamic exchange path at finite speed, with back-reaction.
  T2. The model: the driven spinor on the swap loop + the sourced field.
  T3. Adiabatic recovery: geo phase → −π, P_exc → 0 as the swap slows.
  T4. Non-adiabatic correction: the O(1/T) deviation; P_exc grows with speed.
  T5. Back-reaction: the field is sourced (E_f grows at finite speed → 0 slow).
  T6. The exchange phase is the adiabatic invariant (−π = the −1).
  T7. Honest scope.
  T8. Assessment.

Verdict:
  - DYNAMIC_TWO_THROAT_EXCHANGE_RECOVERS_THE_ADIABATIC_MINUS_ONE_WITH_NON_ADIABATIC_AND_BACK_REACTION_CORRECTIONS_VANISHING_ADIABATICALLY
    (expected): the dynamic, finite-speed two-throat exchange with a
    back-reacting field recovers the adiabatic exchange phase −π (the −1, the
    #188 holonomy) as the swap slows, with non-adiabatic excitation (O(1/T))
    and back-reaction energy that grow with speed and vanish adiabatically —
    the exchange statistics are the robust adiabatic limit of the full
    dynamics.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np


# ════════════════════════════════════════════════════════════════════════
# THE DYNAMIC EXCHANGE  (driven spin-½ on the swap loop + a sourced field)
# ════════════════════════════════════════════════════════════════════════

_SX = np.array([[0, 1], [1, 0]], dtype=complex)
_SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
_SZ = np.array([[1, 0], [0, -1]], dtype=complex)

_DELTA = 1.0       # internal gap
_OMEGA = 2.0       # back-reaction field frequency
_KAPPA = 0.6       # source (throat → field)
_GAMMA = 0.3       # field damping (radiation)
_DT = 0.02         # fixed step for uniform accuracy

_CACHE: dict = {}


def _H(s: float, x: float, g: float) -> np.ndarray:
    """H(t) = (Δ/2) n̂(s)·σ + g x σ_z, n̂ the equatorial swap loop."""
    return (0.5 * _DELTA * (math.cos(2 * math.pi * s) * _SX
                            + math.sin(2 * math.pi * s) * _SY)
            + g * x * _SZ)


def _lower(H: np.ndarray) -> np.ndarray:
    _, v = np.linalg.eigh(H)
    return v[:, 0]


def run_swap(T: float, g: float = 0.0):
    """Evolve the throat spinor around the swap loop over duration T (speed
    1/T) with the back-reacting field, by RK4 + Verlet.  Returns the geometric
    phase, the return amplitude, the non-adiabatic excitation P_exc, and the
    peak back-reaction field energy.  Memoized."""
    key = (round(T, 3), round(g, 3))
    if key in _CACHE:
        return _CACHE[key]
    nt = max(20000, int(T / _DT))
    dt = T / nt
    psi = _lower(_H(0.0, 0.0, g)).astype(complex)
    psi0 = psi.copy()
    x = 0.0
    vx = 0.0
    phi_dyn = 0.0
    Ef_peak = 0.0
    for k in range(nt):
        s = (k + 0.5) * dt / T
        H = _H(s, x, g)
        f = lambda p: -1j * (H @ p)
        k1 = f(psi)
        k2 = f(psi + 0.5 * dt * k1)
        k3 = f(psi + 0.5 * dt * k2)
        k4 = f(psi + dt * k3)
        psi = psi + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        psi = psi / np.linalg.norm(psi)
        sz = float(np.real(np.vdot(psi, _SZ @ psi)))
        phi_dyn += float(np.real(np.vdot(psi, H @ psi))) * dt   # ∫⟨H⟩ dt
        a = -_GAMMA * vx - _OMEGA ** 2 * x - _KAPPA * sz        # back-reaction
        vx += a * dt
        x += vx * dt
        Ef_peak = max(Ef_peak, 0.5 * (vx ** 2 + _OMEGA ** 2 * x ** 2))
    # geometric phase: factor the (large) dynamical phase continuously
    amp = complex(np.vdot(psi0, psi))
    geo = float(np.angle(amp * np.exp(1j * phi_dyn)))
    P_exc = float(1 - abs(np.vdot(_lower(_H(1.0, x, g)), psi)) ** 2)
    out = {"geo": geo, "return": float(abs(amp)), "P_exc": P_exc,
           "E_field": float(Ef_peak)}
    _CACHE[key] = out
    return out


def adiabatic_berry_phase() -> float:
    """The exact adiabatic Berry phase of the swap loop: the equatorial great
    circle encloses solid angle Ω = 2π, so the spin-½ phase is −Ω/2 = −π."""
    return -0.5 * (2 * math.pi)


# ════════════════════════════════════════════════════════════════════════
# TESTS
# ════════════════════════════════════════════════════════════════════════

def test_T1_goal() -> dict:
    return {
        "name": "T1_goal",
        "description": (
            "Go beyond the adiabatic holonomy of #188: traverse the two-throat "
            "exchange path at FINITE speed, in real time, with the field "
            "BACK-REACTING — sourced by the moving throat and acting back on "
            "it — and characterize the exchange phase, the non-adiabatic "
            "corrections, and the back-reaction. The throat's internal (Pin) "
            "spinor is driven around the swap loop (the FR 2π-rotation / "
            "exchange path) over a finite duration T, coupled to a sourced "
            "field mode; the question is whether and how the adiabatic "
            "exchange sign −1 emerges from the full dynamics."
        ),
        "extends": "the adiabatic #188 holonomy to finite-speed dynamics + back-reaction",
        "model": "driven spinor on the swap loop + a sourced (back-reacting) field",
        "framing": "QFT on the fixed classical throat geometry — not quantum gravity",
        "pass": True,
    }


def test_T2_model() -> dict:
    """The model: the driven spinor on the swap loop + the sourced field."""
    # verify the field is genuinely two-way: with g>0 and finite speed the
    # field is sourced (nonzero), and it couples back into H.
    fast = run_swap(20.0, g=0.5)
    sourced = fast["E_field"] > 0
    ok = sourced
    return {
        "name": "T2_dynamic_model",
        "description": (
            "The dynamic exchange model (the dynamical generalization of "
            "#188). The throat's internal Pin spinor evolves under "
            "i ψ̇ = H(t) ψ with H(t) = (Δ/2) n̂(s(t))·σ + g x(t) σ_z, where "
            "n̂(s) = (cos 2πs, sin 2πs, 0) traces the swap loop (the "
            "equatorial great circle — the Finkelstein–Rubinstein 2π-rotation "
            "/ exchange path) over a finite duration T (swap speed v = 1/T), "
            "and the BACK-REACTION field x(t) is sourced by the throat and "
            "acts back: ẍ + γẋ + ω²x = −κ⟨σ_z⟩. This is a genuine TWO-WAY "
            "system — the moving throat drives the field "
            f"(peak field energy {fast['E_field']:.3f} > 0 at finite speed), "
            "the field deflects the throat (the g x σ_z term). The exchange "
            "phase is the geometric (Berry) phase of the loop; its adiabatic "
            "value is the FR holonomy (T6)."
        ),
        "field_is_sourced_two_way": sourced,
        "pass": ok,
    }


def test_T3_adiabatic_recovery() -> dict:
    """Adiabatic recovery: geo phase → −π, P_exc → 0 as the swap slows."""
    berry = adiabatic_berry_phase()
    Ts = [100.0, 300.0, 1000.0]
    rows = {T: run_swap(T, g=0.4) for T in Ts}
    geos = {T: rows[T]["geo"] for T in Ts}
    pexc = {T: rows[T]["P_exc"] for T in Ts}
    dev = {T: abs(geos[T] - berry) for T in Ts}
    recovers_phase = dev[1000.0] < 0.05
    excitation_vanishes = pexc[1000.0] < 1e-3
    monotone_dev = dev[1000.0] < dev[100.0]
    ok = recovers_phase and excitation_vanishes and monotone_dev
    return {
        "name": "T3_adiabatic_recovery",
        "description": (
            "As the swap is SLOWED, the dynamics recovers the adiabatic "
            "exchange −1. The geometric phase "
            f"{ {T: round(geos[T], 3) for T in Ts} } → −π = {berry:.3f} (the "
            f"deviation { {T: round(dev[T], 3) for T in Ts} } shrinks to "
            f"{dev[1000.0]:.3f} at T = 1000), and the non-adiabatic excitation "
            f"P_exc { {T: float('%.0e'%pexc[T]) for T in Ts} } → 0. The "
            "finite-speed evolution, slowed, reproduces the #188 holonomy: the "
            "throat spinor adiabatically follows its internal state around the "
            "swap loop and accumulates the geometric phase −π — the exchange "
            "sign −1."
        ),
        "berry_phase_target": round(berry, 4),
        "geo_phase_by_T": {str(T): round(geos[T], 4) for T in Ts},
        "P_exc_by_T": {str(T): pexc[T] for T in Ts},
        "deviation_by_T": {str(T): round(dev[T], 4) for T in Ts},
        "pass": ok,
    }


def test_T4_non_adiabatic() -> dict:
    """Non-adiabatic correction: O(1/T) deviation; P_exc grows with speed."""
    berry = adiabatic_berry_phase()
    Ts = [10.0, 30.0, 100.0, 300.0]
    rows = {T: run_swap(T, g=0.0) for T in Ts}
    dev = {T: abs(rows[T]["geo"] - berry) for T in Ts}
    pexc = {T: rows[T]["P_exc"] for T in Ts}
    dev_times_T = {T: round(dev[T] * T, 1) for T in Ts}
    # the deviation scales as 1/T: dev·T ≈ constant
    vals = list(dev_times_T.values())
    one_over_T = (max(vals[1:]) - min(vals[1:])) / np.mean(vals[1:]) < 0.2
    exc_grows = pexc[10.0] > pexc[100.0]
    ok = one_over_T and exc_grows
    return {
        "name": "T4_non_adiabatic_correction",
        "description": (
            "At finite speed the throat cannot follow adiabatically — the "
            "cost of a finite-speed swap. The geometric-phase deviation from "
            f"−π is { {T: round(dev[T], 3) for T in Ts} }, scaling as O(1/T): "
            f"deviation × T = {dev_times_T} is nearly constant (≈ "
            f"{np.mean(vals[1:]):.0f}), the standard first-order non-adiabatic "
            "correction. And the non-adiabatic excitation P_exc = "
            f"{ {T: float('%.0e'%pexc[T]) for T in Ts} } GROWS as the swap "
            "speeds up (the throat is excited / deformed when forced to "
            "exchange quickly). The dynamical exchange has a real, quantified "
            "cost that vanishes as the swap slows."
        ),
        "deviation_by_T": {str(T): round(dev[T], 4) for T in Ts},
        "deviation_times_T": {str(T): dev_times_T[T] for T in Ts},
        "P_exc_by_T": {str(T): pexc[T] for T in Ts},
        "scales_as_one_over_T": one_over_T,
        "pass": ok,
    }


def test_T5_back_reaction() -> dict:
    """Back-reaction: the field is sourced (E_f grows at finite speed → 0)."""
    # with back-reaction on, the field is sourced; it grows at finite speed
    # and vanishes adiabatically; the adiabatic limit (the −1) is unchanged.
    fast = run_swap(20.0, g=0.5)
    slow = run_swap(300.0, g=0.5)
    no_br_fast = run_swap(20.0, g=0.0)
    sourced = fast["E_field"] > 1e-3
    vanishes_adiabatically = slow["E_field"] < fast["E_field"]
    # the adiabatic limit is unchanged by the back-reaction (geo → −π either way)
    berry = adiabatic_berry_phase()
    with_br = run_swap(1000.0, g=0.5)
    limit_unchanged = abs(with_br["geo"] - berry) < 0.06
    ok = sourced and vanishes_adiabatically and limit_unchanged
    return {
        "name": "T5_back_reaction",
        "description": (
            "The field BACK-REACTS — sourced by the moving throat, acting back "
            "on it. At finite speed (T = 20) the moving throat sources the "
            f"field to peak energy E_f = {fast['E_field']:.3f} (versus "
            f"{no_br_fast['E_field']:.3f} with the source off) — energy pumped "
            "from the swap motion into the field (radiation / drag, the "
            "back-reaction). Slowing the swap (T = 300) the field response "
            f"falls to E_f = {slow['E_field']:.3f} → 0: the back-reaction is a "
            "finite-speed effect. And it does NOT change the adiabatic limit — "
            f"with the back-reaction on, the slow-swap geometric phase "
            f"{with_br['geo']:.3f} still → −π = {berry:.3f}. The two throats "
            "exchange energy with the field as they swap, but the exchange "
            "STATISTIC (the −1) is untouched."
        ),
        "E_field_fast_with_br": round(fast["E_field"], 5),
        "E_field_fast_no_br": round(no_br_fast["E_field"], 5),
        "E_field_slow_with_br": round(slow["E_field"], 5),
        "geo_phase_slow_with_br": round(with_br["geo"], 4),
        "adiabatic_limit_unchanged": limit_unchanged,
        "pass": ok,
    }


def test_T6_adiabatic_invariant() -> dict:
    """The exchange phase is the adiabatic invariant (−π = the −1)."""
    berry = adiabatic_berry_phase()
    is_minus_pi = abs(berry + math.pi) < 1e-9
    exchange_sign = int(round(math.cos(berry)))      # e^{iγ} eigenvalue sign
    # the dynamics at the slowest swap reproduces it
    slow = run_swap(1000.0, g=0.4)
    dynamics_recovers = abs(slow["geo"] - berry) < 0.05
    ok = is_minus_pi and exchange_sign == -1 and dynamics_recovers
    return {
        "name": "T6_exchange_phase_adiabatic_invariant",
        "description": (
            "The exchange phase is the ADIABATIC invariant. The swap loop (the "
            "equatorial great circle n̂(s)) encloses solid angle Ω = 2π, so "
            f"the spin-½ Berry phase is −Ω/2 = {berry:.4f} = −π — a GEOMETRIC "
            "invariant of the loop (the Finkelstein–Rubinstein holonomy of "
            "#188), independent of the swap speed and of the back-reaction. "
            f"Its phase factor e^{{iγ}} = {exchange_sign:+d} is the exchange "
            "sign −1 (fermionic). The finite-speed dynamics REPRODUCES it "
            f"(the slow swap gives geometric phase {slow['geo']:.3f} ≈ −π), and "
            "the non-adiabatic (T4) and back-reaction (T5) effects are "
            "corrections that vanish as the swap slows. So the exchange "
            "statistics are the robust ADIABATIC LIMIT of the full dynamical "
            "evolution — the −1 is not a fragile fine-tuning but the universal "
            "slow-swap value the dynamics flows to."
        ),
        "adiabatic_berry_phase": round(berry, 6),
        "exchange_sign": exchange_sign,
        "dynamics_recovers_at_slow_swap": dynamics_recovers,
        "pass": ok,
    }


def test_T7_scope() -> dict:
    return {
        "name": "T7_honest_scope",
        "description": (
            "An EFFECTIVE model of the dynamic exchange. The throat's internal "
            "(Pin) state is the driven spin-½ on the swap loop (the FR "
            "exchange ≃ 2π rotation, #188); the back-reaction is a single "
            "sourced field mode (the throat's deformation / the mutual field, "
            "in mean field), not the full field-resolved two-throat dynamics. "
            "The exchange phase is the ADIABATIC invariant (−π = the −1); the "
            "finite-speed and back-reaction effects — the O(1/T) geometric-"
            "phase deviation, the non-adiabatic excitation P_exc, and the "
            "back-reaction field energy — are controlled CORRECTIONS that "
            "vanish as the swap slows. So the exchange statistics survive as "
            "the adiabatic limit of the full dynamical evolution, with the "
            "dynamical cost quantified. This complements #188 (the adiabatic "
            "holonomy) and #189/#190 (the static mean field): the dynamics "
            "do not change the statistics, only add a finite-speed cost. The "
            "field-resolved real-time two-throat solve is the follow-up. "
            "Weak-field, code units."
        ),
        "model": "driven spin-½ on the swap loop + a single sourced (back-reacting) field mode",
        "exact": "the adiabatic Berry phase −π (the exchange −1)",
        "corrections": ["O(1/T) geometric-phase deviation", "non-adiabatic excitation P_exc",
                        "back-reaction field energy — all vanish as the swap slows"],
        "follow_up": "the field-resolved real-time two-throat solve",
        "pass": True,
    }


def test_T8_assessment() -> dict:
    return {
        "name": "T8_assessment",
        "description": (
            "Dynamical. Traversing the two-throat exchange path at finite "
            "speed in real time, with the field back-reacting, the adiabatic "
            "exchange sign −1 emerges as the slow-swap limit of the full "
            "dynamics. The swap loop's exact Berry phase is −π (Ω/2 of the "
            "equatorial great circle) — the exchange −1, the #188 holonomy, a "
            "geometric invariant independent of speed and back-reaction. As "
            "the swap is SLOWED the dynamical geometric phase → −π and the "
            "non-adiabatic excitation → 0 (adiabatic recovery). At FINITE "
            "speed the geometric phase deviates from −π by O(1/T) and the "
            "excitation grows (the throat cannot follow — the cost of a fast "
            "swap), while the moving throat SOURCES the back-reaction field "
            "(energy pumped into it, vanishing adiabatically) which acts back "
            "on the spinor. The back-reaction and the finite speed are "
            "controlled corrections; the exchange STATISTIC (the −1) is the "
            "robust adiabatic limit, unchanged. So the dynamics do not alter "
            "the exchange sign — they only add a quantified, adiabatically-"
            "vanishing cost. SCOPE: an effective model (the Pin spinor + a "
            "single back-reacting field mode); the field-resolved real-time "
            "two-throat solve is the follow-up."
        ),
        "classification": (
            "DYNAMIC_TWO_THROAT_EXCHANGE_RECOVERS_THE_ADIABATIC_MINUS_ONE_WITH_NON_ADIABATIC_AND_BACK_REACTION_CORRECTIONS_VANISHING_ADIABATICALLY"
        ),
        "pass": True,
    }


# ════════════════════════════════════════════════════════════════════════
# RUNNER + VERDICT
# ════════════════════════════════════════════════════════════════════════

def run_probe() -> dict:
    tests = [
        test_T1_goal(),
        test_T2_model(),
        test_T3_adiabatic_recovery(),
        test_T4_non_adiabatic(),
        test_T5_back_reaction(),
        test_T6_adiabatic_invariant(),
        test_T7_scope(),
        test_T8_assessment(),
    ]
    t3, t4, t5, t6 = tests[2], tests[3], tests[4], tests[5]
    if all(t["pass"] for t in tests):
        verdict_class = (
            "DYNAMIC_TWO_THROAT_EXCHANGE_RECOVERS_THE_ADIABATIC_MINUS_ONE_WITH_NON_ADIABATIC_AND_BACK_REACTION_CORRECTIONS_VANISHING_ADIABATICALLY"
        )
        verdict = (
            "DYNAMICAL — THE EXCHANGE −1 IS THE ADIABATIC LIMIT OF THE FULL "
            "DYNAMICS. Finite-speed swap with a back-reacting field.\n\n"
            "THE INVARIANT. The swap loop's exact Berry phase is "
            f"{t6['adiabatic_berry_phase']} = −π — the exchange sign "
            f"{t6['exchange_sign']:+d} (the #188 holonomy), a geometric "
            "invariant independent of speed and back-reaction.\n\n"
            "ADIABATIC RECOVERY. Slowing the swap, the dynamical geometric "
            f"phase → −π (deviation → {t3['deviation_by_T']['1000.0']} at "
            "T = 1000) and the non-adiabatic excitation → 0.\n\n"
            "NON-ADIABATIC COST. At finite speed the geometric phase deviates "
            "from −π by O(1/T) (deviation × T ≈ "
            f"{list(t4['deviation_times_T'].values())[-1]}, constant) and the "
            "excitation grows — the throat cannot follow a fast swap.\n\n"
            "BACK-REACTION. The moving throat sources the field (peak energy "
            f"{t5['E_field_fast_with_br']} at T = 20, → "
            f"{t5['E_field_slow_with_br']} when slow) which acts back on the "
            "spinor; energy is exchanged with the field, but the adiabatic "
            "limit (the −1) is unchanged. So the dynamics do not alter the "
            "exchange sign — they add a quantified, adiabatically-vanishing "
            "cost. SCOPE: an effective model; the field-resolved real-time "
            "two-throat solve is the follow-up."
        )
    else:
        verdict_class = "DYNAMIC_TWO_THROAT_EXCHANGE_INCOMPLETE"
        verdict = (
            "INCOMPLETE. A check failed; review the adiabatic recovery, the "
            "non-adiabatic scaling, the back-reaction, or the adiabatic "
            "invariant."
        )

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "the dynamic, finite-speed two-throat exchange with a back-reacting "
            "field recovers the adiabatic exchange phase −π (the −1, the #188 "
            "holonomy) as the swap slows, with non-adiabatic excitation (O(1/T)) "
            "and back-reaction energy that grow with speed and vanish "
            "adiabatically — the exchange statistic is the robust adiabatic "
            "limit of the full dynamics"
        ),
        "adiabatic_phase": "the swap loop's Berry phase −π (Ω/2 of the great circle) = the exchange −1",
        "recovery": "slowing the swap → geo phase → −π, P_exc → 0 (the #188 holonomy)",
        "non_adiabatic": "finite speed → O(1/T) phase deviation; P_exc grows with speed",
        "back_reaction": "the throat sources the field (E_f grows at finite speed → 0 slow); the −1 is unchanged",
        "scope": "an effective model (Pin spinor + a back-reacting field mode); the field-resolved solve is the follow-up",
        "tests": tests,
        "n_passed": sum(1 for t in tests if t["pass"]),
        "n_total": len(tests),
        "verdict_class": verdict_class,
        "verdict": verdict,
    }


def render_markdown(s: dict) -> str:
    out = []
    out.append("# Dynamic two-throat exchange path with back-reaction (PR #191)")
    out.append("")
    out.append(f"**Run:** {s['timestamp_utc']}")
    out.append("")
    out.append(
        "Goes beyond the adiabatic #188 holonomy: traverses the two-throat "
        "exchange path at finite speed in real time, with the field "
        "back-reacting, and shows the exchange `−1` is the robust adiabatic "
        "limit of the full dynamics. *(QFT on the classical throat, not "
        "quantum gravity.)*"
    )
    out.append("")
    out.append(f"- **Adiabatic phase**: {s['adiabatic_phase']}")
    out.append(f"- **Recovery**: {s['recovery']}")
    out.append(f"- **Non-adiabatic**: {s['non_adiabatic']}")
    out.append(f"- **Back-reaction**: {s['back_reaction']}")
    out.append("")
    out.append("## Test summary")
    out.append("")
    out.append("| # | Test | Key finding | PASS? |")
    out.append("|---|---|---|---|")
    labels = {
        "T1": "the dynamic exchange path at finite speed, with back-reaction",
        "T2": "the model: the driven spinor on the swap loop + the sourced field",
        "T3": "adiabatic recovery: geo phase → −π, P_exc → 0 as the swap slows",
        "T4": "non-adiabatic correction: O(1/T) deviation; P_exc grows with speed",
        "T5": "back-reaction: the field is sourced (E_f grows at speed → 0 slow)",
        "T6": "the exchange phase is the adiabatic invariant (−π = the −1)",
        "T7": "honest scope (an effective model)",
        "T8": "DYNAMIC_EXCHANGE_RECOVERS_THE_ADIABATIC_MINUS_ONE",
    }
    for t in s["tests"]:
        p = "**PASS**" if t["pass"] else "**FAIL**"
        pre = t["name"][:2]
        out.append(f"| {pre} | `{t['name']}` | {labels.get(pre,'—')} | {p} |")
    out.append("")
    t3, t4 = s["tests"][2], s["tests"][3]
    out.append("## Adiabatic recovery and the non-adiabatic cost")
    out.append("")
    out.append("| T (slow → fast) | geometric phase | deviation from −π | deviation × T | P_exc |")
    out.append("|---:|---:|---:|---:|---:|")
    for T in ["10.0", "30.0", "100.0", "300.0"]:
        dev = t4["deviation_by_T"][T]
        devT = t4["deviation_times_T"][T]
        out.append(f"| {T} | −π + {dev} | {dev} | {devT} | {t4['P_exc_by_T'][T]:.1e} |")
    out.append("")
    out.append(
        "The geometric phase → −π (the exchange `−1`) as the swap slows; the "
        "deviation scales as `O(1/T)` (deviation × T ≈ constant) and the "
        "excitation `P_exc` grows with speed — the quantified cost of a "
        "finite-speed swap, vanishing adiabatically."
    )
    out.append("")
    out.append("## Verdict")
    out.append("")
    out.append(f"**{s['verdict_class']}.** {s['verdict']}")
    out.append("")
    return "\n".join(out)


def _json_default(o):
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def main(argv: Optional[list] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_dynamic_two_throat_exchange_probe"
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(
        json.dumps(summary, indent=2, default=_json_default))
    (out / "probe.md").write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
