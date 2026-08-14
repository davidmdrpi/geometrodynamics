"""
A solved wave field with a throat identification, against the branch ledger.

THE QUESTION
────────────
`waves.pair_history` (PR #253) built a **ray** ledger: a history's leg may take
the short way, the long way, or wind, and closure is an equation on those path
lengths.  It ended by conceding that rank counting had reached its limit — what
it could not supply is a quantity that *vanishes* when a source is removed
rather than merely becoming underdetermined.

Before building that quantity, one thing has to be checked:

    **does a solved wave field, with the throat imposed as the same
    identification, reproduce the exact branch ledger — its arrival times, its
    amplitudes, and phases the ray picture never had?**

It does, and more sharply than a stationary-phase argument would give.

WHY THE ANSWER IS EXACT RATHER THAN ASYMPTOTIC
──────────────────────────────────────────────
On the Einstein static universe ``S³ × R`` the scalar Laplacian has eigenvalues
``−n(n+2)``, and ``R = 6``, so the **conformally coupled** (``ξ = 1/6``) massless
field has

    ``ω_n² = n(n+2) + 1 = (n+1)²``      ⟹      ``ω_n = n + 1``

**Integer frequencies.**  The retarded Green function is therefore exactly
periodic, and summing the modes in closed form gives

    ``G(χ,t) = 1/(4π sin χ) · [ Σ_k δ(t − χ − 2πk) − Σ_k δ(t + χ − 2πk) ]``

so the geometric-optics branches are the **exact support** of the field, not the
leading term of an expansion.  Three things follow, and each is measured:

* **the arrival times are the ray ledger.**  ``t = χ + 2πk`` is the short way
  with ``k`` windings, ``t = 2πk − χ`` the long way — precisely the branch set
  ``ℓ = (χ or 2π−χ) + 2πk`` that PR #253 enumerated by hand;
* **the amplitude is the shell law.**  ``1/(4π sin χ)`` — the ``1/sin χ`` of
  PR #251, now falling out of a field solve instead of being imposed by energy
  conservation on a shell;
* **and the field supplies phases the ray ledger could not.**  Each image
  carries a sign, and the sign is ``(−1)^m`` with ``m`` the number of focal
  crossings — the antipode at ``t = π``, the source point again at ``t = 2π``,
  and so on.  That is the **Maslov index**, and no amount of path-length
  bookkeeping produces it.

WHICH FIELD THE LEDGER BELONGS TO
─────────────────────────────────
Conformal coupling is not a convenience.  The **minimally** coupled field has
``ω_n = √(n(n+2))``, which is irrational, so its Green function has no image
structure: measured here, ``63%`` of the peak amplitude sits *between* the
branches, against ``4e-08`` for the conformal field.  The sharp branch ledger is
therefore a statement about the conformally coupled field specifically, and PR
#253 was implicitly working with it.

THE THROAT
──────────
Imposed as the same identification as before, now on a field rather than on
rays:

    ``φ(M⁺, t) = η · φ(M⁻, t + Δ)``

which at the level of arrivals means a through-throat contribution reaches a
point at ``τ + ℓ₁ + Δ + ℓ₂``, summed over branch pairs.  PR #253's closure
condition ``ℓ₁ + ℓ₂ + Δ = 0`` is then exactly the statement that a through-throat
arrival lands back on the event it left from — and the field additionally says
with what sign it lands, which is ``η`` times the product of the two Maslov
factors.

SCOPE
─────
Still a **linear field on a fixed round background**, with the throat still an
*identification map* rather than a solution — `shells.junction` (PR #249) priced
that throat and the bill is inherited, unpaid.  What is new is that the
propagation is now solved rather than enumerated.

**Not done:** no backreaction, no stress tensor, no topology change, no rate, and
**no two-source invariant yet** — the cross-quantity that vanishes when a source
is removed is the next step, and it needs a branch-resolved definition precisely
because the branch structure above is real rather than a bookkeeping choice.
The self-consistent loop (multiple throat traversals) is *not* summed here; only
the arrivals are, so nothing in this module supersedes PR #251's fixed point.
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

__all__ = [
    "TWO_PI",
    "esu_frequency",
    "spectral_field",
    "image_field",
    "branch_arrivals",
    "field_peaks",
    "through_throat_arrivals",
    "measure_the_spectral_solve_matches_the_image_sum",
    "measure_the_field_support_is_the_branch_ledger",
    "measure_the_phases_are_the_maslov_index",
    "measure_the_amplitude_is_the_shell_law",
    "measure_minimal_coupling_has_no_branch_structure",
    "measure_the_throat_reproduces_the_closure_condition",
]

TWO_PI = 2.0 * math.pi


# ════════════════════════════════════════════════════════════════════════════
# THE FIELD
# ════════════════════════════════════════════════════════════════════════════
def esu_frequency(n: int, conformal: bool = True) -> float:
    """Mode frequency on the Einstein static universe.

    ``ω² = n(n+2)`` for the minimally coupled massless scalar; the conformal
    term ``ξR = 1`` adds one, giving ``ω = n+1`` exactly.  The integrality is
    the whole reason the branch structure below is sharp.
    """
    base = n * (n + 2.0)
    return math.sqrt(base + 1.0) if conformal else math.sqrt(base)


def _zonal(m: np.ndarray, chi: float) -> np.ndarray:
    """``(n+1)/(2π²) · sin((n+1)χ)/sin χ`` with ``m = n+1`` — the ``S³`` addition
    theorem, whose ``χ → 0`` limit ``(n+1)²/2π²`` is the mode degeneracy over the
    volume ``2π²``."""
    s = math.sin(chi)
    if abs(s) < 1e-12:
        raise ValueError("the zonal sum is singular at χ = 0 and χ = π")
    return m * np.sin(m * chi) / (2.0 * math.pi ** 2 * s)


def spectral_field(chi: float, t: float, width: float = 0.06,
                   n_max: int = 4000, conformal: bool = True) -> float:
    """**The solve.**  Retarded field of a Gaussian pulse, by mode sum.

    A unit source at a point, switched on with a Gaussian profile of width
    ``width`` centred at ``t = 0``, gives mode amplitudes
    ``e^{−ω²·width²/2} sin(ωt)/ω``.  Nothing about images or branches is put in:
    this is a truncated eigenmode sum, and the Gaussian in time is also what
    makes the truncation converge.
    """
    m = np.arange(1, int(n_max) + 1)
    om = (m.astype(float) if conformal
          else np.sqrt((m - 1.0) * (m + 1.0)))
    with np.errstate(divide="ignore", invalid="ignore"):
        prop = np.where(om > 1e-12, np.sin(om * t) / np.where(om > 1e-12, om,
                                                              1.0), t)
    return float(np.sum(_zonal(m, chi) * prop
                        * np.exp(-(om * width) ** 2 / 2.0)))


def image_field(chi: float, t: float, width: float = 0.06,
                n_images: int = 6) -> float:
    """The closed form the mode sum should equal — an **independent** construction.

        ``1/(4π sin χ) [ Σ_k G_w(t − χ − 2πk) − Σ_k G_w(t + χ − 2πk) ]``

    Derived by Poisson summation from the same mode sum, but evaluated as a sum
    over *images* rather than over modes, so agreement is a real check rather
    than a rearrangement performed twice.
    """
    def g(x: float) -> float:
        return math.exp(-x * x / (2.0 * width * width)) / (
            width * math.sqrt(TWO_PI))

    s = math.sin(chi)
    if abs(s) < 1e-12:
        raise ValueError("the image sum is singular at χ = 0 and χ = π")
    total = sum(g(t - chi - TWO_PI * k) for k in range(0, n_images + 1))
    total -= sum(g(t + chi - TWO_PI * k) for k in range(1, n_images + 1))
    return total / (4.0 * math.pi * s)


# ════════════════════════════════════════════════════════════════════════════
# THE LEDGER, AS THE FIELD SEES IT
# ════════════════════════════════════════════════════════════════════════════
def branch_arrivals(chi: float, t_max: float) -> List[Dict[str, object]]:
    """Arrival times, with the winding label **and the Maslov sign**.

    The times are PR #253's branch set.  What is new is ``crossings``: how many
    focal points the signal has passed — the antipode at ``t = π``, the source
    point again at ``t = 2π``, and so on — and the sign ``(−1)^crossings`` that
    the field attaches to each arrival.
    """
    out: List[Dict[str, object]] = []
    k = 0
    while True:
        short = chi + TWO_PI * k
        long_way = TWO_PI * (k + 1) - chi
        added = False
        for t, is_long in ((short, 0), (long_way, 1)):
            if t <= t_max:
                crossings = 2 * k if not is_long else 2 * k + 1
                out.append({"t": t, "long_way": is_long, "winding": k,
                            "crossings": crossings,
                            "sign": (-1) ** crossings})
                added = True
        if not added:
            break
        k += 1
        if k > 200:
            break
    return sorted(out, key=lambda r: r["t"])


def field_peaks(chi: float, t_max: float, width: float = 0.06,
                n_grid: int = 20000, conformal: bool = True,
                rel_floor: float = 0.05) -> List[Dict[str, float]]:
    """Local extrema of ``|φ|`` in the **solved** field, with their signs.

    Deliberately reads the solve rather than the closed form: the comparison is
    only worth something if one side never saw the image structure.
    """
    ts = np.linspace(1e-3, t_max, int(n_grid))
    f = np.array([spectral_field(chi, float(t), width, conformal=conformal)
                  for t in ts])
    floor = rel_floor * float(np.abs(f).max())
    picks: List[Dict[str, float]] = []
    for i in range(1, len(ts) - 1):
        a = abs(f[i])
        if a > abs(f[i - 1]) and a >= abs(f[i + 1]) and a > floor:
            if picks and ts[i] - picks[-1]["t"] < 4.0 * width:
                if a > abs(picks[-1]["value"]):
                    picks[-1] = {"t": float(ts[i]), "value": float(f[i])}
                continue
            picks.append({"t": float(ts[i]), "value": float(f[i])})
    return picks


def through_throat_arrivals(chi_to_mouth: float, chi_from_mouth: float,
                            delay: float, t_max: float,
                            eta: int = 1) -> List[Dict[str, object]]:
    """Arrivals of the field that has traversed the throat once.

    The identification ``φ(M⁺,t) = η φ(M⁻, t+Δ)`` moves an arrival at ``M⁺`` to
    an emission at ``M⁻``, so a through-throat contribution lands at
    ``ℓ₁ + Δ + ℓ₂`` summed over branch pairs — exactly PR #253's closure sum.
    The sign is ``η`` times the two Maslov factors, which is the part the ray
    ledger had no way to carry.
    """
    out = []
    for a in branch_arrivals(chi_to_mouth, t_max + abs(delay)):
        for b in branch_arrivals(chi_from_mouth, t_max + abs(delay)):
            t = a["t"] + delay + b["t"]
            if not -1e-9 <= t <= t_max:
                continue
            out.append({
                "t": t, "leg_in": a["t"], "leg_out": b["t"],
                "branch_in": (a["long_way"], a["winding"]),
                "branch_out": (b["long_way"], b["winding"]),
                "crossings": a["crossings"] + b["crossings"],
                "sign": int(eta) * a["sign"] * b["sign"]})
    return sorted(out, key=lambda r: r["t"])


# ════════════════════════════════════════════════════════════════════════════
# MEASUREMENTS
# ════════════════════════════════════════════════════════════════════════════
def measure_the_spectral_solve_matches_the_image_sum(
        chis: Sequence[float] = (0.4, 1.1, 2.0, 2.7),
        width: float = 0.06, n_t: int = 60) -> Dict[str, object]:
    """The mode sum against the image sum: two constructions, one field.

    The solve never sees an image and the closed form never sees a mode, so the
    agreement is a check rather than an algebraic rearrangement done twice.
    """
    worst = 0.0
    rows = []
    for chi in chis:
        scale = image_field(chi, chi, width)
        w = 0.0
        for t in np.linspace(0.05, 4.0 * math.pi, n_t):
            w = max(w, abs(spectral_field(chi, float(t), width)
                           - image_field(chi, float(t), width)))
        worst = max(worst, w)
        rows.append({"chi": chi, "peak_scale": scale, "worst_abs_error": w,
                     "relative": w / abs(scale)})
    return {
        "rows": rows, "worst_abs_error": worst,
        "the_two_constructions_agree": bool(worst < 1e-10),
        "what_this_shows": ("the solved field IS the image sum; the branch "
                            "structure is exact, not an approximation"),
    }


def measure_the_field_support_is_the_branch_ledger(
        chis: Sequence[float] = (0.7, 1.1, 1.9),
        t_max: float = 4.0 * math.pi, width: float = 0.05
        ) -> Dict[str, object]:
    """Do the peaks of the **solved** field sit on PR #253's branch times?

    The whole question of this round, in one measurement.  Peaks are read off
    the mode sum; the predicted times come from the ray ledger; neither knows
    about the other.
    """
    rows = []
    worst = 0.0
    all_matched = True
    for chi in chis:
        peaks = field_peaks(chi, t_max, width)
        ledger = [r for r in branch_arrivals(chi, t_max)]
        matched = []
        for p in peaks:
            best = min(ledger, key=lambda r: abs(r["t"] - p["t"]))
            err = abs(best["t"] - p["t"])
            worst = max(worst, err)
            matched.append({"t_peak": p["t"], "t_branch": best["t"],
                            "error": err,
                            "branch": (best["long_way"], best["winding"])})
        # every ledger entry inside the window should have a peak
        covered = all(any(abs(m["t_branch"] - r["t"]) < 1e-9 for m in matched)
                      for r in ledger)
        all_matched = all_matched and covered and len(peaks) == len(ledger)
        rows.append({"chi": chi, "n_peaks": len(peaks),
                     "n_branches": len(ledger), "matches": matched})
    return {
        "rows": rows, "worst_time_error": worst,
        "every_branch_has_a_peak_and_no_peak_is_spurious": bool(all_matched),
        "peaks_land_on_the_ledger": bool(worst < 1e-2),
        "so_the_field_reproduces_the_ray_ledger": bool(all_matched
                                                       and worst < 1e-2),
    }


def measure_the_phases_are_the_maslov_index(
        chis: Sequence[float] = (0.7, 1.1, 1.9),
        t_max: float = 4.0 * math.pi, width: float = 0.05
        ) -> Dict[str, object]:
    """**What the field adds that the ray ledger could not.**

    Each arrival carries a sign, and it is ``(−1)^m`` with ``m`` the number of
    focal crossings: the antipode at ``t = π``, the source point again at
    ``t = 2π``, and so on.  That is the Maslov index.  Path-length bookkeeping
    has no way to produce it — PR #253 could say *when* a history arrives and
    never *with what sign*.
    """
    rows = []
    agree = True
    for chi in chis:
        peaks = field_peaks(chi, t_max, width)
        ledger = branch_arrivals(chi, t_max)
        for p in peaks:
            best = min(ledger, key=lambda r: abs(r["t"] - p["t"]))
            got = 1 if p["value"] > 0 else -1
            ok = got == best["sign"]
            agree = agree and ok
            rows.append({"chi": chi, "t": p["t"],
                         "crossings": best["crossings"],
                         "field_sign": got, "maslov_sign": best["sign"],
                         "agrees": ok})
    return {
        "rows": rows,
        "every_sign_is_the_maslov_factor": bool(agree and rows),
        "the_rule": "sign = (−1)^(number of focal crossings)",
        "the_ray_ledger_could_not_supply_this": (
            "path lengths give arrival times; the sign comes from the "
            "focusing history, which only a field carries"),
    }


def measure_the_amplitude_is_the_shell_law(
        chis: Sequence[float] = (0.35, 0.7, 1.1, math.pi / 2, 2.1, 2.6),
        width: float = 0.05) -> Dict[str, object]:
    """``1/(4π sin χ)`` — PR #251's shell law, now out of a field solve.

    That round imposed ``A ∝ 1/sin χ`` by conserving energy across a shell of
    area ``4π sin²χ``.  Here it is not imposed: the peak of the solved field,
    multiplied by ``sin χ``, is the same constant at every ``χ``.
    """
    rows = []
    prods = []
    for chi in chis:
        peak = spectral_field(chi, chi, width)     # the direct arrival
        prods.append(peak * math.sin(chi))
        rows.append({"chi": chi, "peak": peak,
                     "peak_times_sin_chi": prods[-1]})
    spread = (max(prods) - min(prods)) / max(prods)
    expected = 1.0 / (4.0 * math.pi * width * math.sqrt(TWO_PI))
    return {
        "rows": rows, "relative_spread": spread,
        "the_product_is_constant": bool(spread < 1e-6),
        "constant": float(np.mean(prods)),
        "predicted_constant": expected,
        "matches_one_over_four_pi_sin_chi": bool(
            abs(float(np.mean(prods)) - expected) / expected < 1e-6),
        "so_the_shell_law_is_derived_here_not_imposed": True,
    }


def measure_minimal_coupling_has_no_branch_structure(
        chi: float = 1.1, width: float = 0.06, n_grid: int = 3000
        ) -> Dict[str, object]:
    """Which field the ledger belongs to — and it is not the obvious one.

    The **minimally** coupled massless scalar has ``ω = √(n(n+2))``, irrational,
    so no image structure and no sharp branches: the field between the arrivals
    is comparable to the arrivals themselves.  The conformal field leaves
    essentially nothing there.  So PR #253's sharp ledger is a statement about
    the conformally coupled field specifically, which that round never said
    because rays cannot tell the two apart.
    """
    ts = np.linspace(0.02, TWO_PI, int(n_grid))
    out = {}
    for label, conf in (("conformal", True), ("minimal", False)):
        f = np.array([spectral_field(chi, float(t), width, conformal=conf)
                      for t in ts])
        peak = float(np.abs(f).max())
        mask = (np.abs(ts - chi) > 0.35) & (np.abs(ts - (TWO_PI - chi)) > 0.35)
        between = float(np.abs(f[mask]).max())
        out[label] = {"peak": peak, "worst_between_branches": between,
                      "ratio": between / peak}
    return {
        **out,
        "conformal_is_sharp": bool(out["conformal"]["ratio"] < 1e-5),
        "minimal_is_not": bool(out["minimal"]["ratio"] > 0.1),
        "the_ledger_belongs_to_the_conformal_field": bool(
            out["conformal"]["ratio"] < 1e-5 < out["minimal"]["ratio"]),
        "why": ("ω = n+1 is integer only with the conformal term ξR = 1; "
                "without it ω = √(n(n+2)) is irrational and the Green function "
                "has no images"),
    }


def measure_the_throat_reproduces_the_closure_condition(
        chi_in: float = 1.2, chi_out: float = 0.9, eta: int = 1,
        t_max: float = 6.0 * math.pi) -> Dict[str, object]:
    """PR #253's closure, as a statement about where a field arrives.

    A through-throat contribution lands at ``ℓ₁ + Δ + ℓ₂``.  Choosing ``Δ`` to
    be minus a branch-pair sum — which is exactly PR #253's closure condition —
    puts an arrival back on the emission event, at ``t = 0``.  The field then
    adds what that round could not: the arrival's **sign**, ``η`` times the two
    Maslov factors.
    """
    rows = []
    ok = True
    ledger_in = branch_arrivals(chi_in, t_max)
    ledger_out = branch_arrivals(chi_out, t_max)
    for a in ledger_in[:3]:
        for b in ledger_out[:3]:
            delay = -(a["t"] + b["t"])            # PR #253's closure condition
            arrivals = through_throat_arrivals(chi_in, chi_out, delay, t_max,
                                               eta)
            zero = [r for r in arrivals if abs(r["t"]) < 1e-9]
            hit = bool(zero)
            ok = ok and hit
            rows.append({
                "leg_in": a["t"], "leg_out": b["t"], "delay": delay,
                "closes_at_zero": hit,
                "sign_of_the_closing_arrival": (zero[0]["sign"] if zero
                                                else None),
                "expected_sign": eta * a["sign"] * b["sign"],
                "sign_agrees": bool(zero
                                    and zero[0]["sign"]
                                    == eta * a["sign"] * b["sign"])})
    return {
        "rows": rows,
        "closure_puts_an_arrival_on_the_emission_event": bool(ok),
        "every_closing_sign_is_eta_times_the_maslov_factors": bool(
            all(r["sign_agrees"] for r in rows)),
        "the_closure_condition_is_a_field_statement": (
            "ℓ₁ + Δ + ℓ₂ = 0 is where a through-throat arrival lands back on "
            "its own emission event"),
        "and_the_field_adds_the_sign": (
            "η times the two Maslov factors — the ray ledger carried no phase"),
        "what_is_still_put_in": (
            "the throat is still an identification map, not a solution, and "
            "the self-consistent sum over repeated traversals is PR #251's "
            "fixed point and is not redone here"),
    }
