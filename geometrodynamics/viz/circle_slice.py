"""
A great-circle slice of the scalar wave, and a bulk it can wrap through.

THE SLICE
─────────
Take the scalar field on ``S²`` and cut it with the great circle through the
source and its antipode.  Parametrised by ``σ ∈ [−π, π)``, the geodesic
distance from the source along that circle is simply ``d = |σ|`` — so the
circle carries **both** halves of the wave at once: two lobes running in
opposite directions, meeting head-on at ``σ = ±π``, which is the antipode.

Nothing is re-solved.  The field is the 2-D solve sampled at ``d(σ)``, so the
slice inherits the real caustic — including the ``1/√(sin d)`` focusing that a
naive 1-D wave equation on a circle would not have.  This matters: on a circle
two counter-propagating pulses would merely superpose to ``2×``; on the sphere
they arrive as a collapsing ring and the amplitude does something else.

THE BULK
────────
The slice is drawn inside the vacuole — an annulus between ``R_inner`` and
``R_outer``.  The crossing rule is the obvious one: a radius that would pass
``R_outer`` re-enters at ``R_inner``, so the wave that reaches up into the bulk
comes back **inside** the circle.  That glues the two boundaries and makes the
radial coordinate a circle of circumference ``gap = R_outer − R_inner``; the
drawn plane is an annulus, but the space the curve really lives on is a
**torus** ``S¹_σ × S¹_ρ``.

WHAT THAT BUYS, AND WHAT IT DOES NOT
────────────────────────────────────
On a torus a closed curve has a winding number in each direction, and winding
numbers are integers that continuous motion cannot change — exactly the sort of
stable object a crossing rule is supposed to produce.

It does not produce one here, and the reason is worth stating plainly.  The
curve is a **graph** ``r = f(σ)`` over the circle with ``f`` single-valued, so
its radial winding number is identically zero at every amplitude and every
time.  Every outward crossing of the seam is matched by an inward one; the
signed count is exactly ``0``, no matter how many unsigned crossings the
amplitude buys.

So the seam can be crossed as often as you like without ever accumulating
topological charge.  A height field cannot wind.  That is a constraint on any
bulk crossing rule of this kind, and it says where a stable topological object
would have to come from instead: not from the amplitude of a scalar height,
but from a curve that is free to stop being a graph.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from geometrodynamics.viz.throat_wavefront import BareSphereSim
from geometrodynamics.viz.warped_sphere import NestedShells

TWO_PI = 2.0 * math.pi
ANTIPODAL_TIME = math.pi
RETURN_TIME = 2.0 * math.pi


# ════════════════════════════════════════════════════════════════════════════
# THE BULK
# ════════════════════════════════════════════════════════════════════════════
@dataclass
class BulkAnnulus:
    """The vacuole, with its outer boundary glued to its inner one.

    The gluing is the whole content: it makes the radial direction periodic
    with period ``gap``, so ``(σ, ρ)`` is a torus rather than an annulus.  A
    radius driven past ``R_outer`` reappears just inside ``R_inner`` — the wave
    that reaches into the bulk wraps around and shows up *inside* the circle.
    """

    shells: NestedShells = None  # type: ignore[assignment]

    def __post_init__(self) -> None:
        if self.shells is None:
            self.shells = NestedShells()

    @property
    def r_inner(self) -> float:
        return self.shells.r_inner

    @property
    def r_outer(self) -> float:
        return self.shells.r_outer

    @property
    def r_mid(self) -> float:
        return self.shells.r_mid

    @property
    def gap(self) -> float:
        return self.r_outer - self.r_inner

    # ── the rule ────────────────────────────────────────────────────────────
    def wrap(self, r) -> Tuple[np.ndarray, np.ndarray]:
        """Fold a radius into the annulus; report how many times it wrapped.

        Returns ``(drawn, sheet)`` where ``drawn ∈ [R_inner, R_outer)`` is
        where the point is painted and ``sheet`` is the integer count of
        boundary crossings it took to get there — positive outward.
        """
        a = np.asarray(r, dtype=float)
        offset = (a - self.r_inner) / self.gap
        sheet = np.floor(offset)
        return self.r_inner + (offset - sheet) * self.gap, sheet.astype(int)

    def reaches_the_seam(self, r) -> bool:
        a = np.asarray(r, dtype=float)
        return bool(np.any(a >= self.r_outer) or np.any(a < self.r_inner))


# ════════════════════════════════════════════════════════════════════════════
# THE SLICE
# ════════════════════════════════════════════════════════════════════════════
class CircleSlice:
    """The great circle through both poles, carrying the scalar wave.

    Parameters
    ----------
    sim
        The scalar solve.  Defaults to a ``BareSphereSim`` — the field is a
        function of geodesic distance alone, so the slice samples it exactly.
    gain
        Display amplitude ``ε``.  ``None`` calibrates from the run's own peak
        so the wave uses a fixed fraction of the half-gap.
    fill
        With ``gain=None``, the fraction of the half-gap the peak should use.
        Above ``1.0`` the wave crosses the seam and wraps.
    """

    def __init__(self, sim: Optional[BareSphereSim] = None,
                 bulk: Optional[BulkAnnulus] = None,
                 n_sigma: int = 721, gain: Optional[float] = None,
                 fill: float = 0.72, pulse_width: float = 0.18,
                 n_radial: int = 900) -> None:
        self.sim = sim if sim is not None else BareSphereSim(
            n_theta=8, n_phi=8, pulse_width=pulse_width, n_radial=n_radial)
        self.bulk = bulk or BulkAnnulus()
        self.n_sigma = int(n_sigma)
        # closed: σ = −π and σ = +π are the same point, both kept so the drawn
        # curve closes and the winding count sees the whole loop
        self.sigma = np.linspace(-math.pi, math.pi, self.n_sigma)
        self.fill = float(fill)
        self.gain = float(gain) if gain is not None else self._calibrate()

    # ── clock ───────────────────────────────────────────────────────────────
    @property
    def t(self) -> float:
        return self.sim.t

    def reset(self) -> None:
        self.sim.reset()

    def advance_to(self, t_target: float) -> None:
        self.sim.advance_to(t_target)

    def _calibrate(self, samples: int = 200,
                   t_end: float = RETURN_TIME) -> float:
        self.sim.reset()
        peak = 0.0
        for i in range(samples):
            self.sim.advance_to((i + 1) * t_end / samples)
            peak = max(peak, float(np.max(np.abs(self.field()))))
        self.sim.reset()
        half = 0.5 * self.bulk.gap
        return self.fill * half / max(peak, 1e-30)

    # ── the field on the circle ─────────────────────────────────────────────
    def distance(self) -> np.ndarray:
        """``d = |σ|`` — the geodesic distance from the source along the slice."""
        return np.abs(self.sigma)

    def field(self) -> np.ndarray:
        """The 2-D solve, sampled on the slice.  Not a separate 1-D wave."""
        return self.sim.field_at_distance(self.distance())

    def radius(self, gain: Optional[float] = None) -> np.ndarray:
        """The unwrapped radius — before the crossing rule is applied."""
        eps = self.gain if gain is None else float(gain)
        return self.bulk.r_mid + eps * self.field()

    def points(self, gain: Optional[float] = None) -> np.ndarray:
        """The drawn curve in ℝ², after wrapping.  Breaks at the seam."""
        drawn, _ = self.bulk.wrap(self.radius(gain=gain))
        return np.stack([drawn * np.cos(self.sigma),
                         drawn * np.sin(self.sigma)], axis=-1)

    def segments(self, gain: Optional[float] = None) -> List[np.ndarray]:
        """The drawn curve, split wherever it crosses the seam.

        The curve is continuous on the torus and discontinuous in the drawn
        plane; splitting is a drawing decision, not a physical one.
        """
        drawn, sheet = self.bulk.wrap(self.radius(gain=gain))
        xy = np.stack([drawn * np.cos(self.sigma),
                       drawn * np.sin(self.sigma)], axis=-1)
        cuts = np.nonzero(np.diff(sheet) != 0)[0]
        out, start = [], 0
        for c in cuts:
            if c + 1 > start:
                out.append(xy[start:c + 1])
            start = c + 1
        out.append(xy[start:])
        return [s for s in out if len(s) > 1]

    def sheet(self, gain: Optional[float] = None) -> np.ndarray:
        return self.bulk.wrap(self.radius(gain=gain))[1]

    # ── topology ────────────────────────────────────────────────────────────
    def seam_crossings(self, gain: Optional[float] = None) -> Dict[str, object]:
        """Signed and unsigned crossings of the glued boundary.

        The signed total is the radial winding number.  It is zero for a graph
        — which is the point.
        """
        sheet = self.sheet(gain=gain)
        steps = np.diff(sheet)
        outward = int(np.sum(steps[steps > 0]))
        inward = int(-np.sum(steps[steps < 0]))
        return {
            "outward": outward,
            "inward": inward,
            "unsigned": outward + inward,
            "signed": outward - inward,
            "sheets_visited": int(np.ptp(sheet)) + 1,
        }

    def winding_number(self, gain: Optional[float] = None) -> int:
        """Degree of ``σ ↦ ρ`` on the torus, from unwrapped increments.

        Computed independently of ``seam_crossings`` — by summing shortest-path
        steps around the loop — so the two are a cross-check rather than one
        number reported twice.
        """
        drawn, _ = self.bulk.wrap(self.radius(gain=gain))
        rho = (drawn - self.bulk.r_inner) / self.bulk.gap
        step = np.diff(rho)
        step -= np.round(step)                       # shortest way round
        return int(round(float(np.sum(step))))

    def excursion(self, gain: Optional[float] = None) -> Dict[str, float]:
        r = self.radius(gain=gain)
        b = self.bulk
        return {
            "max_radius": float(np.max(r)),
            "min_radius": float(np.min(r)),
            "outward_fraction": float((np.max(r) - b.r_mid)
                                      / (b.r_outer - b.r_mid)),
            "inward_fraction": float((b.r_mid - np.min(r))
                                     / (b.r_mid - b.r_inner)),
            "wraps": bool(b.reaches_the_seam(r)),
        }

    # ── the unrolled picture ────────────────────────────────────────────────
    def unrolled(self, gain: Optional[float] = None
                 ) -> Tuple[np.ndarray, np.ndarray]:
        """``(σ, ρ)`` on the torus, with ``ρ ∈ [0, 1)`` in units of the gap."""
        drawn, _ = self.bulk.wrap(self.radius(gain=gain))
        return self.sigma, (drawn - self.bulk.r_inner) / self.bulk.gap


# ════════════════════════════════════════════════════════════════════════════
# MEASUREMENTS
# ════════════════════════════════════════════════════════════════════════════
def measure_the_slice_is_the_sphere(slice_: Optional[CircleSlice] = None,
                                    t: float = 1.4) -> Dict[str, object]:
    """The slice is a cut of the 2-D solve, not a 1-D wave of its own.

    Checked two ways: the field at ``σ`` equals the sphere's field at geodesic
    distance ``|σ|`` measured through the full ``(θ, φ)`` machinery, and the
    slice is exactly symmetric under ``σ → −σ`` because both arms are the same
    distance from a point source.
    """
    s = slice_ or CircleSlice()
    s.reset()
    s.advance_to(t)
    u = s.field()

    # ...the same points, reached through the (θ, φ) route instead
    sim = s.sim
    n_hat = np.array([math.sin(sim.source_theta) * math.cos(sim.source_phi),
                      math.sin(sim.source_theta) * math.sin(sim.source_phi),
                      math.cos(sim.source_theta)])
    tmp = np.array([0.0, 0.0, 1.0])
    if abs(float(tmp @ n_hat)) > 0.9:
        tmp = np.array([1.0, 0.0, 0.0])
    e = tmp - float(tmp @ n_hat) * n_hat
    e = e / np.linalg.norm(e)
    pts = (np.cos(s.sigma)[:, None] * n_hat[None, :]
           + np.sin(s.sigma)[:, None] * e[None, :])
    theta = np.arccos(np.clip(pts[:, 2], -1.0, 1.0))
    phi = np.arctan2(pts[:, 1], pts[:, 0])
    via_sphere = sim.field_at(theta, phi)

    peak = float(np.max(np.abs(u))) or 1.0
    mismatch = float(np.max(np.abs(u - via_sphere))) / peak
    mirror = float(np.max(np.abs(u - u[::-1]))) / peak
    return {
        "time": s.t,
        "peak": peak,
        "worst_mismatch_against_the_sphere": mismatch,
        "worst_mirror_asymmetry": mirror,
        "is_a_cut_of_the_sphere": bool(mismatch < 1e-9),
        "both_arms_are_the_same": bool(mirror < 1e-9),
    }


def measure_the_lobes_meet_at_the_antipode(slice_: Optional[CircleSlice] = None,
                                           frames: int = 400
                                           ) -> Dict[str, object]:
    """Two lobes leave the source in opposite directions and arrive together.

    On the slice the antipode is a single point ``σ = ±π``, so the meeting is
    a head-on collision — and the amplitude there is the sphere's caustic, not
    a 1-D superposition, because the field being sampled is the 2-D solve.
    """
    s = slice_ or CircleSlice()
    s.reset()
    best = {"amp": -math.inf, "time": 0.0}
    launch = float(np.max(np.abs(s.field())))
    half = s.n_sigma // 2
    for i in range(frames):
        s.advance_to((i + 1) * 1.15 * ANTIPODAL_TIME / frames)
        u = s.field()
        # the lobe positions: peak of each arm, excluding the source cell
        amp_at_antipode = abs(float(u[0]))       # σ = −π ≡ +π
        if amp_at_antipode > best["amp"]:
            best = {"amp": amp_at_antipode, "time": s.t}
        left = float(s.sigma[int(np.argmax(np.abs(u[:half])))])
        right = float(s.sigma[half + int(np.argmax(np.abs(u[half:])))])
    return {
        "launch_amplitude": launch,
        "antipodal_amplitude": best["amp"],
        "arrival_time": best["time"],
        "antipodal_gain": best["amp"] / max(launch, 1e-30),
        "final_left_lobe": left,
        "final_right_lobe": right,
        "lobes_are_symmetric": bool(abs(left + right) < 0.05),
        "arrives_on_time": bool(abs(best["time"] - ANTIPODAL_TIME) < 0.25),
    }


def measure_the_wrap_threshold(slice_: Optional[CircleSlice] = None,
                               frames: int = 200) -> Dict[str, object]:
    """The gain at which the wave first reaches the seam, predicted and found.

    ``ε_crit = (R_outer − R_mid) / max|u|`` — and the measured threshold is
    searched for rather than assumed.
    """
    s = slice_ or CircleSlice()
    s.reset()
    peak = 0.0
    for i in range(frames):
        s.advance_to((i + 1) * RETURN_TIME / frames)
        peak = max(peak, float(np.max(np.abs(s.field()))))
    predicted = (s.bulk.r_outer - s.bulk.r_mid) / max(peak, 1e-30)

    lo, hi = 0.0, 4.0 * predicted
    for _ in range(60):                     # bisection on "does it wrap?"
        mid = 0.5 * (lo + hi)
        wraps = False
        s.reset()
        for i in range(frames):
            s.advance_to((i + 1) * RETURN_TIME / frames)
            if s.bulk.reaches_the_seam(s.radius(gain=mid)):
                wraps = True
                break
        if wraps:
            hi = mid
        else:
            lo = mid
    measured = 0.5 * (lo + hi)
    return {
        "run_peak": peak,
        "predicted_threshold": predicted,
        "measured_threshold": measured,
        "relative_error": abs(measured - predicted) / max(predicted, 1e-30),
        "display_gain": s.gain,
        "display_gain_is_below_threshold": bool(s.gain < predicted),
        "threshold_is_the_half_gap_over_the_peak": bool(
            abs(measured - predicted) < 5e-3 * predicted),
    }


def measure_the_seam_is_crossed_in_pairs(
        slice_: Optional[CircleSlice] = None,
        gains: Sequence[float] = (1.0, 1.6, 2.4, 3.6, 5.0),
        frames: int = 160) -> Dict[str, object]:
    """However often the wave wraps, the signed count is exactly zero.

    This is the load-bearing negative result.  The drawn curve is a graph
    ``r = f(σ)`` over the circle with ``f`` single-valued, so on the torus its
    radial winding number is identically zero: every outward crossing of the
    glued boundary is paid for by an inward one.  Amplitude buys unsigned
    crossings and never buys topological charge.

    Both the crossing ledger and the independently-computed degree are
    reported, so the claim is cross-checked rather than restated.
    """
    s = slice_ or CircleSlice()
    thr = measure_the_wrap_threshold(s, frames=frames)["predicted_threshold"]
    rows = []
    worst_signed, worst_winding, most_unsigned = 0, 0, 0
    for g in gains:
        eps = float(g) * thr
        s.reset()
        best = {"unsigned": 0, "signed": 0, "winding": 0, "time": 0.0,
                "sheets": 1}
        for i in range(frames):
            s.advance_to((i + 1) * RETURN_TIME / frames)
            c = s.seam_crossings(gain=eps)
            w = s.winding_number(gain=eps)
            worst_signed = max(worst_signed, abs(c["signed"]))
            worst_winding = max(worst_winding, abs(w))
            if c["unsigned"] > best["unsigned"]:
                best = {"unsigned": c["unsigned"], "signed": c["signed"],
                        "winding": w, "time": s.t,
                        "sheets": c["sheets_visited"]}
        most_unsigned = max(most_unsigned, best["unsigned"])
        rows.append({"gain_over_threshold": float(g), "gain": eps, **best})
    return {
        "threshold": thr,
        "rows": rows,
        "most_unsigned_crossings": most_unsigned,
        "worst_signed_total": worst_signed,
        "worst_winding_number": worst_winding,
        "amplitude_buys_crossings": bool(most_unsigned > 0),
        "the_signed_count_is_always_zero": bool(worst_signed == 0),
        "a_height_field_cannot_wind": bool(worst_winding == 0),
    }


def measure_different_waves_wrap_differently(
        widths: Sequence[float] = (0.36, 0.24, 0.14, 0.08),
        gain: float = 0.80, frames: int = 200) -> Dict[str, object]:
    """Different scalar waves meet the same bulk at the **same** gain.

    Normalising each pulse by its own wrap threshold would make them identical
    by construction — every wave would reach the same multiple of the gap — so
    they are all driven at one fixed gain instead, and the bulk is the same
    annulus for all of them.

    What separates them is not the crossing *count*.  Every pulse of a given
    sign crosses the seam in and out again, so counts move in twos and are
    nearly width-independent; the launch amplitude is ``1`` for all of them and
    their peaks barely differ.  What does separate them is **how much of the
    circle is on the far sheet** — a wide pulse carries a broad arc through the
    bulk, a narrow one pokes a thin spike through — and that is the quantity
    reported here.

    The winding column is zero for every one of them, which is why they are
    worth putting side by side.
    """
    rows = []
    for w in widths:
        s = CircleSlice(pulse_width=float(w), gain=gain)
        s.reset()
        peak_arc, peak_reach, most, sheets, at_time = 0.0, 0.0, 0, 1, 0.0
        winding, signed = 0, 0
        for i in range(frames):
            s.advance_to((i + 1) * RETURN_TIME / frames)
            sheet = s.sheet(gain=gain)
            arc = float(np.mean(sheet != 0))       # fraction of the circle away
            c = s.seam_crossings(gain=gain)
            r = s.radius(gain=gain)
            winding = max(winding, abs(s.winding_number(gain=gain)))
            signed = max(signed, abs(c["signed"]))
            most = max(most, c["unsigned"])
            sheets = max(sheets, c["sheets_visited"])
            peak_reach = max(peak_reach,
                             float((np.max(r) - s.bulk.r_mid) / s.bulk.gap))
            if arc > peak_arc:
                peak_arc, at_time = arc, s.t
        rows.append({"pulse_width": float(w), "max_crossings": most,
                     "sheets": sheets, "arc_on_the_far_sheet": peak_arc,
                     "at_time": at_time, "reach_in_gaps": peak_reach,
                     "winding": winding, "signed": signed})
    arcs = [r["arc_on_the_far_sheet"] for r in rows]
    # the arc is the pulse: (arc fraction × 2π) / pulse width, near-constant
    per_width = [r["arc_on_the_far_sheet"] * TWO_PI / r["pulse_width"]
                 for r in rows]
    for row, k in zip(rows, per_width):
        row["arc_over_pulse_width"] = k
    return {
        "gain": float(gain),
        "rows": rows,
        "crossing_counts": [r["max_crossings"] for r in rows],
        "arc_range": (min(arcs), max(arcs)),
        "arc_ratio": max(arcs) / max(min(arcs), 1e-30),
        "mean_arc_over_pulse_width": sum(per_width) / len(per_width),
        "arc_over_pulse_width_spread": max(per_width) - min(per_width),
        "the_far_sheet_arc_is_the_pulse": bool(
            max(per_width) - min(per_width) < 0.35 * (
                sum(per_width) / len(per_width))),
        "wide_pulses_carry_a_broader_arc": bool(arcs[0] > arcs[-1]),
        "they_differ_in_arc_not_in_count": bool(
            max(arcs) > 1.5 * min(arcs)
            and len(set(r["max_crossings"] for r in rows)) <= 2),
        "none_of_them_winds": bool(all(r["winding"] == 0 for r in rows)
                                   and all(r["signed"] == 0 for r in rows)),
    }


def measure_the_curve_closes(slice_: Optional[CircleSlice] = None,
                             t: float = 1.4) -> Dict[str, object]:
    """``σ = −π`` and ``σ = +π`` are the same point of the sphere.

    So the drawn curve closes exactly, wrapped or not — which is what makes
    the winding number a well-posed integer in the first place.
    """
    s = slice_ or CircleSlice()
    s.reset()
    s.advance_to(t)
    thr = (s.bulk.r_outer - s.bulk.r_mid) / max(
        float(np.max(np.abs(s.field()))), 1e-30)
    out = {}
    for label, eps in (("below", 0.5 * thr), ("above", 3.0 * thr)):
        r = s.radius(gain=eps)
        drawn, _ = s.bulk.wrap(r)
        out[f"{label}_endpoint_gap"] = abs(float(r[0] - r[-1]))
        out[f"{label}_drawn_gap"] = abs(float(drawn[0] - drawn[-1]))
        out[f"{label}_winding"] = s.winding_number(gain=eps)
    out["closes_exactly"] = bool(out["below_endpoint_gap"] < 1e-12
                                 and out["above_endpoint_gap"] < 1e-12)
    return out
