"""
What a drawn wave has to obey, and what a height field cannot.

Three complaints motivated this module, and they do not all land the same way.

CAUSALITY — THE FRONT IS FINE, THE OFFSET IS NOT
─────────────────────────────────────────────────
"Why is the antipode moving before the wavefront arrives?"  The front is not the
culprit: the amplitude at the antipode is ``5e-133`` at launch, ``6e-40`` at
``t = 1.4`` and ``7e-17`` at ``t = 2.0``, rising to ``O(1)`` only at ``t ≈ π``.
Signal ahead of the front sits at the leapfrog's own noise floor, ``~1e-07``.
Nothing propagates faster than the front.

What *does* move instantly is a **constant**.  On a closed surface
``d²/dt² ∫u dA = ∫Δu dA = 0``, so the mean field is at worst linear in time — and
``throat_wavefront._outgoing_velocity`` already kills the linear term by keeping
``∫u_t dA = 0`` with a correction that stays inside the pulse.  What it cannot
kill is the *constant*: a Gaussian of width ``w`` carries a monopole ``w²/4``
(measured ``0.008056`` against ``0.008100`` predicted), and ``ℓ = 0`` has
``ω = 0``, so that offset never propagates, never radiates and never decays.  It
displaces the whole sphere — antipode included — from ``t = 0`` forever.

But — and this correction matters — the monopole is the sphere's *average*, not
the antipode's value.  Ahead of the front the higher-``ℓ`` modes cancel it
exactly, which is why the far side really is dark.  What the offset does instead
is **permanent**: it is a DC displacement that never radiates away, so the
surface never returns to undeformed.  The time-averaged displacement is
``0.0081`` at every point of the sphere, and the quietest instant of a whole
run still leaves ``max|u| = 0.094``.

So the complaint lands, but one step over from where it was aimed: nothing moves
early, and nothing ever stops moving.

THE FIX HAS TO STAY INSIDE THE PULSE
────────────────────────────────────
Subtracting the mean from ``f`` is the obvious repair and it is worse than the
disease: it leaves ``−w²/4`` at the antipode, exactly the displacement it was
meant to remove, and spreads it over the entire far side.  The correction has to
be *localised*, which a difference of two bumps is:

    ``f = e^{−(d/w)²} − c · e^{−(d/kw)²}`` ,   ``c`` fixed by ``∫f dA = 0``.

Monopole ``4e-19``, antipode ``2e-34``, and still localised.  This is the same
argument the velocity correction already makes, applied to the position.

AND WHY THE MODE IS UNPHYSICAL ANYWAY
─────────────────────────────────────
Electromagnetism has no monopole radiation and gravity has none at ``ℓ = 0`` or
``ℓ = 1``.  So the mode that produces the instant global response is precisely
the one real radiation is forbidden to have — it is not a small blemish on the
analogy, it is outside it.  The spin-2 field of ``spin2_tidal`` is already clean
here: ``h = sin²d·q`` kills ``ℓ = 0`` and ``ℓ = 1`` identically.

THE RING DOES CONCENTRATE — BUT NOT WHERE THE EYE IS LOOKING
────────────────────────────────────────────────────────────
"Height should increase as the ring compresses."  It does, and correctly:
amplitude on the converging front follows ``A ∝ 1/√(sin d)``, which is exactly
``A²·(circumference) = const``.  Two things hide it.

* **The growth is nearly all in the last stretch.**  ``1/√(sin d)`` is flat
  across the middle of the trip and only diverges within a pulse width of the
  antipode, so for most of the run there is nothing to see.
* **On a compact surface the launch is itself a focus.**  The source is a point
  and so is the antipode, so the focal height cannot much exceed the launch
  height — measured ``1.05×``.  The ``5.97×`` growth is relative to the
  *equator*, not the launch.  An expectation of unbounded growth belongs to an
  open geometry, not this one.

HEIGHT IS NOT ENERGY
────────────────────
The sharpest statement of why a deforming surface misleads.  The energy density
is ``u̇² + |∇u|²``, so a *constant* offset contributes **zero** energy while
displacing every point of the surface.  The monopole is the extreme case: it is
the most global feature of the drawn shape and carries none of the energy.  The
eye is drawn to height; the physics is in the gradient.
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from geometrodynamics.viz.antipodal_focusing import SphereWaveSim
from geometrodynamics.viz.throat_wavefront import BareSphereSim

ANTIPODAL_TIME = math.pi
RETURN_TIME = 2.0 * math.pi


# ════════════════════════════════════════════════════════════════════════════
# A SOURCE WITH NO MONOPOLE
# ════════════════════════════════════════════════════════════════════════════
def monopole_free_profile(d: np.ndarray, w: float, k: float = 2.0,
                          with_derivative: bool = False):
    """A localised profile whose mean over the sphere is zero.

    A difference of two concentric bumps: the wide one is subtracted with just
    enough weight to cancel the narrow one's monopole.  Both are localised, so
    unlike "subtract the mean" the correction never touches the far side.
    """
    if k <= 1.0:
        raise ValueError("k must exceed 1 so the two bumps differ")
    sd = np.sin(d)
    narrow = np.exp(-((d / w) ** 2))
    wide = np.exp(-((d / (k * w)) ** 2))
    c = (float(np.trapezoid(narrow * sd, d))
         / float(np.trapezoid(wide * sd, d)))
    if not with_derivative:
        return narrow - c * wide, c
    # −∂f/∂d analytically.  A numerical gradient here is not a rounding
    # detail: it seeds high-wavenumber grid modes that the leapfrog propagates
    # at the wrong speed, and they show up as precursor signal on the far side
    # five orders of magnitude above the analytic launch.
    outgoing = (2.0 * d / (w * w) * narrow
                - c * 2.0 * d / ((k * w) ** 2) * wide)
    return narrow - c * wide, c, outgoing


def compact_monopole_free_profile(d: np.ndarray, a: float = 0.30,
                                  k: float = 2.0):
    """Monopole-free **and** compactly supported — exactly zero past ``k·a``.

    A Gaussian never vanishes, so "nothing ahead of the front" is only ever true
    to the size of its tail, and the wide bump needed to cancel a monopole has a
    tail a hundred orders of magnitude fatter than the narrow one it corrects.
    Building both from a smooth finite-support bump ``(1 − x²)⁴`` removes the
    question: the far side is not small, it is zero, and stays zero until the
    front arrives.

    Returns ``(f, c, -∂f/∂d)`` with the derivative analytic.
    """
    if k <= 1.0:
        raise ValueError("k must exceed 1 so the two bumps differ")
    if not 0.0 < k * a < math.pi:
        raise ValueError("the support k·a must fit inside the sphere")

    def bump(x):
        y = np.zeros_like(x)
        m = np.abs(x) < 1.0
        y[m] = (1.0 - x[m] ** 2) ** 4
        return y

    def dbump(x, scale):
        y = np.zeros_like(x)
        m = np.abs(x) < 1.0
        y[m] = -8.0 * x[m] / scale * (1.0 - x[m] ** 2) ** 3
        return y

    sd = np.sin(d)
    narrow, wide = bump(d / a), bump(d / (k * a))
    c = (float(np.trapezoid(narrow * sd, d))
         / float(np.trapezoid(wide * sd, d)))
    f = narrow - c * wide
    outgoing = -(dbump(d / a, a) - c * dbump(d / (k * a), k * a))
    return f, c, outgoing


def sphere_monopole(u: np.ndarray, d: np.ndarray) -> float:
    """``∫u dA / ∫dA`` on the sphere — the mean the closed geometry conserves."""
    sd = np.sin(d)
    return float(np.trapezoid(u * sd, d) / np.trapezoid(sd, d))


class CleanSphereSim:
    """A scalar pulse launched one-way **and** with no monopole in either slot.

    ``BareSphereSim`` already launches one-way with ``∫u_t dA = 0``.  This adds
    the missing condition ``∫u dA = 0``, using a localised profile so the far
    side is left at zero rather than displaced by the correction.

    ``monopole_free=False`` reproduces the ordinary launch exactly, so the two
    can be compared on one clock.
    """

    def __init__(self, n_radial: int = 2400, pulse_width: float = 0.18,
                 k: float = 2.0, monopole_free: bool = True,
                 compact: bool = True) -> None:
        self.pulse_width = float(pulse_width)
        self.k = float(k)
        self.monopole_free = bool(monopole_free)
        self.compact = bool(compact)
        self._sim = SphereWaveSim(n=n_radial)
        self.reset()

    # ── launch ──────────────────────────────────────────────────────────────
    def reset(self) -> None:
        s = self._sim
        d, w = s.theta, self.pulse_width
        sd = np.sin(d)
        if self.monopole_free and self.compact:
            f, _, grad = compact_monopole_free_profile(d, a=w * 5.0 / 3.0,
                                                       k=self.k)
        elif self.monopole_free:
            f, _, grad = monopole_free_profile(d, w, self.k,
                                               with_derivative=True)
        else:
            f = np.exp(-((d / w) ** 2))
            grad = 2.0 * d / (w * w) * f
        # One-way: u_t = −∂f/∂d, corrected inside the pulse so ∫u_t dA = 0.
        # The corrector cannot be ``f`` itself once ``f`` is monopole-free —
        # that is exactly the profile with no mean to correct with, and using
        # it leaves the mean free to ramp.  A narrow Gaussian is localised and
        # has a mean, so it does both jobs.
        corrector = np.exp(-((d / w) ** 2))
        denom = float(np.sum(corrector * sd))
        c = -float(np.sum(grad * sd)) / denom
        v = grad + c * corrector
        s.u = f.copy()
        s.u_prev = f - s.dt * v
        s.t = 0.0

    @property
    def t(self) -> float:
        return self._sim.t

    @property
    def d(self) -> np.ndarray:
        return self._sim.theta

    @property
    def u(self) -> np.ndarray:
        return self._sim.u

    def advance_to(self, t: float) -> None:
        self._sim.advance_to(t)

    def field_at_distance(self, d) -> np.ndarray:
        return self._sim.sample(np.asarray(d, dtype=float))

    def monopole(self) -> float:
        return sphere_monopole(self._sim.u, self._sim.theta)

    def energy_density(self) -> np.ndarray:
        s = self._sim
        u_t = (s.u - s.u_prev) / s.dt
        return u_t ** 2 + np.gradient(s.u, s.dth) ** 2


# ════════════════════════════════════════════════════════════════════════════
# C1 — THE FRONT IS CAUSAL
# ════════════════════════════════════════════════════════════════════════════
def measure_the_front_is_causal(pulse_width: float = 0.18,
                                n_radial: int = 2400,
                                margin: float = 4.0) -> Dict[str, object]:
    """Nothing outruns the front — checked against the scheme's own noise floor.

    The amplitude at the antipode is tracked from launch: if the front were
    leaking ahead of itself the far side would light up early.  It does not; it
    climbs through 130 orders of magnitude in step with ``t → π``.
    """
    sim = BareSphereSim(n_theta=8, n_phi=8, pulse_width=pulse_width,
                        n_radial=n_radial)
    sim.reset()
    d = np.linspace(0.0, math.pi, 3000)
    rows, ahead_max = [], 0.0
    for t in (0.0, 0.4, 0.8, 1.4, 2.0, 2.4, 2.8, 3.0):
        if t > 0.0:
            sim.advance_to(t)
        u = sim.field_at_distance(d)
        beyond = d > (sim.t + margin * pulse_width)
        ahead = float(np.max(np.abs(u[beyond]))) if beyond.any() else float("nan")
        if beyond.any():
            ahead_max = max(ahead_max, ahead)
        rows.append({"t": sim.t, "antipode": float(u[-1]),
                     "ahead_of_the_front": ahead})
    early = [r for r in rows if r["t"] < 2.2]
    return {
        "rows": rows,
        "worst_ahead_of_the_front": ahead_max,
        "antipode_at_launch": rows[0]["antipode"],
        "antipode_at_two": [r["antipode"] for r in rows if abs(r["t"] - 2.0) < 0.1][0],
        "antipode_at_pi": rows[-1]["antipode"],
        "nothing_outruns_the_front": bool(ahead_max < 1e-5),
        "the_antipode_stays_dark_until_the_front": bool(
            all(abs(r["antipode"]) < 1e-10 for r in early)),
    }


# ════════════════════════════════════════════════════════════════════════════
# C2 — THE CONSTANT THAT MOVES EVERYTHING AT ONCE
# ════════════════════════════════════════════════════════════════════════════
def measure_the_constant_monopole(pulse_width: float = 0.18,
                                  n_radial: int = 2400,
                                  frames: int = 60) -> Dict[str, object]:
    """A Gaussian carries ``w²/4`` of monopole, and a closed surface keeps it.

    ``ℓ = 0`` has ``ω = 0``, so the offset neither propagates nor decays.  It is
    *not* an early response — ahead of the front the higher modes cancel it
    exactly, which ``measure_the_front_is_causal`` confirms across 130 orders of
    magnitude.  It is a **permanent** rigid displacement: a mode with nowhere to
    go, so the surface never returns to its undeformed state.
    """
    sim = BareSphereSim(n_theta=8, n_phi=8, pulse_width=pulse_width,
                        n_radial=n_radial)
    sim.reset()
    d = np.linspace(0.0, math.pi, 4000)
    monos = []
    for i in range(frames):
        sim.advance_to((i + 1) * RETURN_TIME / frames)
        monos.append(sphere_monopole(sim.field_at_distance(d), d))
    sim.reset()
    m0 = sphere_monopole(sim.field_at_distance(d), d)
    predicted = pulse_width ** 2 / 4.0
    drift = (max(monos) - min(monos)) / max(abs(m0), 1e-30)
    return {
        "pulse_width": pulse_width,
        "monopole_at_launch": m0,
        "predicted_w_squared_over_four": predicted,
        "relative_error_against_prediction": abs(m0 - predicted) / predicted,
        "monopole_drift_over_a_full_return": drift,
        "peak_amplitude": 1.0,
        "offset_as_a_fraction_of_peak": abs(m0),
        "the_monopole_is_w_squared_over_four": bool(
            abs(m0 - predicted) < 0.02 * predicted),
        "it_never_changes": bool(drift < 1e-3),
    }


def measure_the_monopole_carries_no_energy(pulse_width: float = 0.18
                                           ) -> Dict[str, object]:
    """The most global feature of the drawn surface carries none of the energy.

    The energy density is ``u̇² + |∇u|²``, and a constant has neither.  So the
    offset that displaces every point of the sphere contributes exactly zero —
    which is the cleanest statement of why drawing a field as a height puts the
    eye on the wrong quantity.
    """
    n = 2000
    d = np.linspace(0.0, math.pi, n)
    dd = float(d[1] - d[0])
    offset = pulse_width ** 2 / 4.0
    flat = np.full_like(d, offset)
    energy_of_offset = float(np.sum((np.gradient(flat, dd) ** 2) * np.sin(d))) * dd

    pulse = np.exp(-((d / pulse_width) ** 2))
    energy_of_pulse = float(np.sum((np.gradient(pulse, dd) ** 2) * np.sin(d))) * dd
    return {
        "offset": offset,
        "energy_of_the_offset": energy_of_offset,
        "energy_of_the_pulse": energy_of_pulse,
        "offset_energy_fraction": energy_of_offset / max(energy_of_pulse, 1e-30),
        "offset_displaces_every_point": True,
        "it_carries_no_energy": bool(energy_of_offset < 1e-20 * energy_of_pulse),
    }


def measure_a_localised_correction_is_the_only_one_that_works(
        pulse_width: float = 0.18, n: int = 6000) -> Dict[str, object]:
    """Three candidate sources, and only one leaves the far side alone.

    The plain Gaussian carries the monopole.  Subtracting its mean removes the
    monopole and displaces the antipode by exactly the amount it removed — the
    disease moved, not cured.  A difference of two localised bumps removes it
    and leaves the antipode at zero, because the correction never leaves the
    pulse.
    """
    d = np.linspace(0.0, math.pi, n)
    rows = []

    g = np.exp(-((d / pulse_width) ** 2))
    rows.append({"source": "gaussian", "monopole": sphere_monopole(g, d),
                 "antipode": float(g[-1]), "peak": float(np.max(np.abs(g)))})

    shifted = g - sphere_monopole(g, d)
    rows.append({"source": "gaussian minus its mean",
                 "monopole": sphere_monopole(shifted, d),
                 "antipode": float(shifted[-1]),
                 "peak": float(np.max(np.abs(shifted)))})

    f, c = monopole_free_profile(d, pulse_width, k=2.0)
    rows.append({"source": "difference of two bumps",
                 "monopole": sphere_monopole(f, d), "antipode": float(f[-1]),
                 "peak": float(np.max(np.abs(f))), "coefficient": c})

    by = {r["source"]: r for r in rows}
    return {
        "rows": rows,
        "gaussian_monopole": by["gaussian"]["monopole"],
        "mean_subtracted_antipode": by["gaussian minus its mean"]["antipode"],
        "localised_monopole": by["difference of two bumps"]["monopole"],
        "localised_antipode": by["difference of two bumps"]["antipode"],
        "subtracting_the_mean_just_moves_it": bool(
            abs(by["gaussian minus its mean"]["antipode"])
            > 0.5 * abs(by["gaussian"]["monopole"])),
        "the_localised_source_is_clean": bool(
            abs(by["difference of two bumps"]["monopole"]) < 1e-12
            and abs(by["difference of two bumps"]["antipode"]) < 1e-12),
    }


def measure_the_clean_source_leaves_the_antipode_alone(
        pulse_width: float = 0.18, n_radial: int = 2400,
        frames: int = 60) -> Dict[str, object]:
    """Three launches on one clock: the monopole, and the far side before arrival.

    The ordinary Gaussian is beautifully quiet on the far side and carries the
    monopole.  Cancelling the monopole with a *wider Gaussian* fixes the mean and
    costs four orders of magnitude of far-side quiet, because the corrector's
    tail is fatter than the pulse it corrects.  Built from compactly supported
    bumps instead, both conditions hold at once — and the far side is quieter
    than the original, because it is exactly zero rather than merely small.
    """
    rows = []
    for label, kw in (("gaussian", dict(monopole_free=False)),
                      ("gaussian, monopole-free",
                       dict(monopole_free=True, compact=False)),
                      ("compact, monopole-free",
                       dict(monopole_free=True, compact=True))):
        sim = CleanSphereSim(n_radial=n_radial, pulse_width=pulse_width, **kw)
        worst, monos = 0.0, []
        for i in range(frames):
            sim.advance_to((i + 1) * 0.6 * ANTIPODAL_TIME / frames)
            worst = max(worst, abs(float(sim.u[-1])))
            monos.append(sim.monopole())
        rows.append({"source": label, "antipode_before_arrival": worst,
                     "monopole": float(np.mean(monos))})
    by = {r["source"]: r for r in rows}
    plain, wide, compact = (by["gaussian"], by["gaussian, monopole-free"],
                            by["compact, monopole-free"])
    return {
        "rows": rows,
        "gaussian_monopole": plain["monopole"],
        "compact_monopole": compact["monopole"],
        "gaussian_far_side": plain["antipode_before_arrival"],
        "compact_far_side": compact["antipode_before_arrival"],
        "wide_gaussian_far_side": wide["antipode_before_arrival"],
        "a_wider_corrector_costs_far_side_quiet": bool(
            wide["antipode_before_arrival"]
            > 1e3 * plain["antipode_before_arrival"]),
        "the_compact_source_wins_on_both": bool(
            abs(compact["monopole"]) < 1e-3 * abs(plain["monopole"])
            and compact["antipode_before_arrival"]
            <= plain["antipode_before_arrival"]),
    }


# ════════════════════════════════════════════════════════════════════════════
# C3 — THE RING DOES CONCENTRATE, WHERE YOU CANNOT SEE IT
# ════════════════════════════════════════════════════════════════════════════
def measure_the_ring_bookkeeping(pulse_width: float = 0.18,
                                 n_radial: int = 2400, frames: int = 300
                                 ) -> Dict[str, object]:
    """``A²·(circumference)`` is what a converging ring conserves.

    Equivalent to ``A ∝ 1/√(sin d)``, and reported as the spread of
    ``A²·sin d`` along the converging leg rather than as a fitted exponent.
    """
    sim = BareSphereSim(n_theta=8, n_phi=8, pulse_width=pulse_width,
                        n_radial=n_radial)
    sim.reset()
    d = np.linspace(0.0, math.pi, 2000)
    rows = []
    for i in range(frames):
        sim.advance_to((i + 1) * ANTIPODAL_TIME / frames)
        u = np.abs(sim.field_at_distance(d))
        j = int(np.argmax(u))
        rows.append({"t": sim.t, "front": float(d[j]), "height": float(u[j])})
    conv = [r for r in rows
            if 0.5 * math.pi < r["front"] < math.pi - 3.0 * pulse_width]
    invariant = [r["height"] ** 2 * math.sin(r["front"]) for r in conv]
    mean = sum(invariant) / len(invariant)
    return {
        "n_samples": len(conv),
        "mean_A2_sin_d": mean,
        "spread": (max(invariant) - min(invariant)) / mean,
        "first_front": conv[0]["front"], "last_front": conv[-1]["front"],
        "energy_on_the_ring_is_conserved": bool(
            (max(invariant) - min(invariant)) / mean < 0.25),
    }


def measure_the_growth_is_invisible_for_most_of_the_trip(
        pulse_width: float = 0.18, n_radial: int = 2400, frames: int = 300,
        tolerance: float = 0.5) -> Dict[str, object]:
    """``1/√(sin d)`` is flat across the middle and only diverges at the end.

    So "the ring compresses and the height grows" is true and nearly invisible:
    the fraction of the trip over which the height stays within ``tolerance`` of
    its equatorial minimum is reported, along with where the growth actually is.
    """
    sim = BareSphereSim(n_theta=8, n_phi=8, pulse_width=pulse_width,
                        n_radial=n_radial)
    sim.reset()
    d = np.linspace(0.0, math.pi, 2000)
    rows = []
    for i in range(frames):
        sim.advance_to((i + 1) * ANTIPODAL_TIME / frames)
        u = np.abs(sim.field_at_distance(d))
        j = int(np.argmax(u))
        rows.append({"t": sim.t, "front": float(d[j]), "height": float(u[j])})
    outbound = [r for r in rows if 0.2 < r["t"] < 0.85 * ANTIPODAL_TIME]
    floor = min(r["height"] for r in outbound)
    flat = [r for r in rows if r["height"] <= (1.0 + tolerance) * floor]
    peak = max(rows, key=lambda r: r["height"])
    return {
        "equator_height": floor,
        "peak_height": peak["height"],
        "growth_from_the_equator": peak["height"] / floor,
        "tolerance": tolerance,
        "flat_fraction_of_the_trip": len(flat) / len(rows),
        "growth_happens_after_t": max(r["t"] for r in flat),
        "growth_window_fraction": 1.0 - len(flat) / len(rows),
        "most_of_the_trip_shows_nothing": bool(len(flat) / len(rows) > 0.6),
    }


def measure_launch_and_focus_are_the_same(
        widths: Sequence[float] = (0.30, 0.18, 0.10),
        n_radial: int = 2400, frames: int = 300) -> Dict[str, object]:
    """On a compact surface the source is a focus too, so the focus cannot win.

    The antipodal height is compared with the *launch* height rather than with
    the equator.  The ratio stays near one for every pulse width, which is why
    an expectation of unbounded growth does not belong on a sphere — the two
    ends of the trip are geometrically the same situation.
    """
    rows = []
    for w in widths:
        sim = BareSphereSim(n_theta=8, n_phi=8, pulse_width=float(w),
                            n_radial=n_radial)
        sim.reset()
        d = np.linspace(0.0, math.pi, 2000)
        launch = float(np.max(np.abs(sim.field_at_distance(d))))
        peak = 0.0
        for i in range(frames):
            sim.advance_to((i + 1) * ANTIPODAL_TIME / frames)
            peak = max(peak, float(np.max(np.abs(sim.field_at_distance(d)))))
        rows.append({"pulse_width": float(w), "launch_height": launch,
                     "focal_height": peak, "ratio": peak / launch})
    ratios = [r["ratio"] for r in rows]
    return {
        "rows": rows,
        "mean_ratio": sum(ratios) / len(ratios),
        "worst_ratio": max(ratios),
        "the_focus_does_not_beat_the_launch": bool(max(ratios) < 1.6),
        "unbounded_growth_needs_an_open_geometry": True,
    }


def measure_the_offset_never_leaves(pulse_width: float = 0.18,
                                    n_radial: int = 2400,
                                    returns: float = 2.0,
                                    frames: int = 800) -> Dict[str, object]:
    """The surface never goes home: a DC displacement with nowhere to radiate.

    Two readings of the same fact.  The time-averaged displacement at every
    point converges on ``w²/4`` — the ``ℓ = 0`` amplitude, which oscillating
    modes cannot contribute to.  And the *quietest instant* of a whole run still
    leaves the sphere visibly deformed, because the offset is always there.

    This is the honest form of the complaint about the far side: nothing moves
    before the front, and nothing ever settles after it.
    """
    sim = BareSphereSim(n_theta=8, n_phi=8, pulse_width=pulse_width,
                        n_radial=n_radial)
    sim.reset()
    probes = np.array([0.3, 1.0, 1.6, 2.4, math.pi - 1e-6])
    d = np.linspace(0.0, math.pi, 1500)
    acc = np.zeros_like(probes)
    quietest = math.inf
    t_end = returns * RETURN_TIME
    for i in range(frames):
        sim.advance_to((i + 1) * t_end / frames)
        acc += sim.field_at_distance(probes)
        quietest = min(quietest, float(np.max(np.abs(sim.field_at_distance(d)))))
    averages = acc / frames
    predicted = pulse_width ** 2 / 4.0
    return {
        "probe_distances": probes.tolist(),
        "time_averaged_displacement": averages.tolist(),
        "predicted_offset": predicted,
        "far_side_average": float(averages[-1]),
        "relative_error_at_the_antipode": abs(float(averages[-1]) - predicted)
        / predicted,
        "quietest_moment": quietest,
        "every_point_is_offset": bool(
            all(a > 0.5 * predicted for a in averages)),
        "the_surface_never_returns_home": bool(quietest > predicted),
    }
