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

    The gluing makes the radial direction periodic, so ``(σ, ρ)`` is a torus
    rather than an annulus: a radius driven past ``R_outer`` reappears inside
    ``R_inner``, and the wave that reaches into the bulk shows up *inside* the
    circle.

    **How** it is glued is a choice, and not a cosmetic one.

    ``mode="translate"``
        Shift the radius by ``gap``.  Radial *offsets* are carried across
        unchanged.  But the two boundary circles have different circumferences,
        so a feature emerging at ``R_inner`` keeps its radial height while
        sitting on a shorter arc: its shape is distorted by ``R_outer/R_inner``.
        Worse, the sheets are arithmetically spaced, so going inward they march
        straight through ``r = 0`` into negative radius.

    ``mode="conformal"``
        Multiply the radius by ``R_inner/R_outer`` — a translation in
        ``ln r`` instead of in ``r``.  Radial offsets scale with the boundary
        they cross, so an emerging feature is a faithfully **scaled copy**:
        aspect ratio preserved exactly.  The sheets are geometrically spaced,
        accumulating at the origin without ever reaching it, so no radius is
        ever negative.

    The two agree to first order in the excursion and diverge as it grows.
    ``translate`` is kept as the default so earlier results are reproducible;
    ``conformal`` is the one that respects the geometry it is folding.
    """

    shells: NestedShells = None  # type: ignore[assignment]
    mode: str = "translate"

    def __post_init__(self) -> None:
        if self.shells is None:
            self.shells = NestedShells()
        if self.mode not in ("translate", "conformal"):
            raise ValueError("mode must be 'translate' or 'conformal'")

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

    @property
    def scale_per_crossing(self) -> float:
        """Radius factor picked up crossing outward through ``R_outer``.

        ``1`` for the translate rule — it shifts rather than scales — and
        ``R_inner/R_outer`` for the conformal one.
        """
        return 1.0 if self.mode == "translate" else self.r_inner / self.r_outer

    # ── the coordinate the seam translates in ───────────────────────────────
    def to_w(self, r) -> np.ndarray:
        a = np.asarray(r, dtype=float)
        if self.mode == "translate":
            return a
        if np.any(a <= 0.0):
            raise ValueError("the conformal seam needs a positive radius; "
                             "pair it with the multiplicative radial law")
        return np.log(a)

    def from_w(self, w) -> np.ndarray:
        a = np.asarray(w, dtype=float)
        return a if self.mode == "translate" else np.exp(a)

    @property
    def period(self) -> float:
        return (self.gap if self.mode == "translate"
                else math.log(self.r_outer / self.r_inner))

    # ── the rule ────────────────────────────────────────────────────────────
    def wrap(self, r) -> Tuple[np.ndarray, np.ndarray]:
        """Fold a radius into the annulus; report how many times it wrapped.

        Returns ``(drawn, sheet)`` where ``drawn ∈ [R_inner, R_outer)`` is
        where the point is painted and ``sheet`` is the integer count of
        boundary crossings it took to get there — positive outward.
        """
        w0 = float(self.to_w(np.array([self.r_inner]))[0])
        offset = (self.to_w(r) - w0) / self.period
        sheet = np.floor(offset)
        return self.from_w(w0 + (offset - sheet) * self.period), sheet.astype(int)

    def sheet_edges(self, n: int = 4) -> Dict[str, List[float]]:
        """Where the boundaries of the neighbouring sheets sit, in real radius.

        Arithmetic for the translate rule — which walks through ``r = 0`` — and
        geometric for the conformal one, which cannot.
        """
        w0 = float(self.to_w(np.array([self.r_inner]))[0])
        out = [float(self.from_w(np.array([w0 + (k + 1) * self.period]))[0])
               for k in range(n)]
        inward = []
        for k in range(n):
            w = w0 - (k + 1) * self.period
            inward.append(float(self.from_w(np.array([w]))[0])
                          if self.mode == "conformal"
                          else self.r_inner - (k + 1) * self.period)
        return {"outward": out, "inward": inward}

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
                 n_radial: int = 900,
                 radial_law: str = "additive") -> None:
        self.sim = sim if sim is not None else BareSphereSim(
            n_theta=8, n_phi=8, pulse_width=pulse_width, n_radial=n_radial)
        self.bulk = bulk or BulkAnnulus()
        if radial_law not in ("additive", "multiplicative"):
            raise ValueError("radial_law must be 'additive' or 'multiplicative'")
        self.radial_law = radial_law
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
        """The unwrapped radius — before the crossing rule is applied.

        ``additive`` is ``R_mid + εu``: the obvious reading, and the one that
        can be driven to a **negative** radius.  ``multiplicative`` is
        ``R_mid·exp(εu)``, which is manifestly positive and is the law that
        pairs with a conformal seam — the two agree to first order in ``ε``.
        """
        eps = self.gain if gain is None else float(gain)
        u = self.field()
        if self.radial_law == "additive":
            return self.bulk.r_mid + eps * u
        return self.bulk.r_mid * np.exp(eps * u)

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
        _, rho = self.unrolled(gain=gain)
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
        """``(σ, ρ)`` on the torus, ``ρ ∈ [0, 1)`` in units of the seam period.

        The period is the gap for a translate seam and ``ln(R_outer/R_inner)``
        for a conformal one, so the strip is the right flat picture either way.
        """
        b = self.bulk
        drawn, _ = b.wrap(self.radius(gain=gain))
        w0 = float(b.to_w(np.array([b.r_inner]))[0])
        return self.sigma, (b.to_w(drawn) - w0) / b.period


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


# ════════════════════════════════════════════════════════════════════════════
# HOW THE SEAM IS GLUED — the scaling is a choice, and it shows
# ════════════════════════════════════════════════════════════════════════════
def measure_the_emerging_feature_keeps_its_shape(
        shells: Optional[NestedShells] = None, height: float = 0.10,
        half_width: float = 0.05) -> Dict[str, object]:
    """A feature crossing the seam: does it come back a scaled copy, or squashed?

    The two boundary circles have different circumferences, so carrying a
    radial *offset* across unchanged — the translate rule — lands the feature
    on a shorter arc with its full height intact.  Its aspect ratio is
    multiplied by ``R_outer/R_inner``: the emerging wave is not the same wave.

    The conformal rule scales the offset by the same factor as the boundary, so
    height and arc length shrink together and the shape is preserved exactly.
    That is the difference between folding the annulus and folding the geometry
    the annulus has.
    """
    sh = shells or NestedShells()
    ri, ro = sh.r_inner, sh.r_outer
    before = height / (ro * half_width)
    rows = []
    for mode in ("translate", "conformal"):
        b = BulkAnnulus(sh, mode=mode)
        drawn, sheet = b.wrap(np.array([ro + height]))
        emerged_height = float(drawn[0]) - ri
        after = emerged_height / (ri * half_width)
        rows.append({
            "mode": mode,
            "height_before": height,
            "height_after": emerged_height,
            "arc_before": ro * half_width,
            "arc_after": ri * half_width,
            "aspect_before": before,
            "aspect_after": after,
            "aspect_distortion": after / before,
            "sheet": int(sheet[0]),
        })
    by = {r["mode"]: r for r in rows}
    return {
        "circumference_ratio": ro / ri,
        "rows": rows,
        "translate_distortion": by["translate"]["aspect_distortion"],
        "conformal_distortion": by["conformal"]["aspect_distortion"],
        "translate_squashes_the_feature": bool(
            abs(by["translate"]["aspect_distortion"] - ro / ri) < 1e-9),
        "conformal_preserves_the_shape": bool(
            abs(by["conformal"]["aspect_distortion"] - 1.0) < 1e-9),
    }


def measure_the_translate_rule_runs_off_the_origin(
        shells: Optional[NestedShells] = None, n: int = 4) -> Dict[str, object]:
    """Going inward, arithmetic sheets reach ``r = 0`` and then go negative.

    That is not a rounding problem, it is the rule: subtracting a fixed ``gap``
    from a radius has nothing to stop it.  Geometric sheets accumulate at the
    origin and never arrive, which is what a radius should do.
    """
    sh = shells or NestedShells()
    out = {}
    for mode in ("translate", "conformal"):
        edges = BulkAnnulus(sh, mode=mode).sheet_edges(n=n)
        out[mode] = edges
        out[f"{mode}_min_inward"] = min(edges["inward"])
        out[f"{mode}_goes_nonpositive"] = bool(min(edges["inward"]) <= 0.0)
    return {
        **out,
        "translate_reaches_negative_radius": out["translate_goes_nonpositive"],
        "conformal_stays_positive": bool(not out["conformal_goes_nonpositive"]),
    }


def measure_the_rule_does_not_rescue_the_winding(
        gains: Sequence[float] = (0.4, 1.0, 3.0, 12.0, 60.0),
        frames: int = 120) -> Dict[str, object]:
    """The zero is not an artefact of how the seam was glued.

    This is the check the scaling objection deserves: rebuild the whole thing on
    a conformal seam with a multiplicative radial law — a genuinely different
    identification, with a different sheet structure and a different notion of
    size — and ask again.  The winding number is still identically zero, for
    the same reason as before: ``ρ(σ)`` comes from a single-valued function on
    the circle whichever coordinate the seam translates in.
    """
    rows = []
    worst = 0
    for mode, law in (("translate", "additive"), ("conformal", "multiplicative"),
                      ("conformal", "additive")):
        s = CircleSlice(bulk=BulkAnnulus(mode=mode), radial_law=law,
                        n_sigma=721)
        s.reset()
        best = {"unsigned": 0, "signed": 0, "winding": 0}
        blew_up = False
        for i in range(frames):
            s.advance_to((i + 1) * RETURN_TIME / frames)
            for g in gains:
                try:
                    c = s.seam_crossings(gain=g)
                    w = s.winding_number(gain=g)
                except ValueError:
                    blew_up = True          # a negative radius, honestly reported
                    continue
                worst = max(worst, abs(c["signed"]), abs(w))
                if c["unsigned"] > best["unsigned"]:
                    best = {"unsigned": c["unsigned"], "signed": c["signed"],
                            "winding": w}
        rows.append({"seam": mode, "radial_law": law, **best,
                     "hit_a_nonpositive_radius": blew_up})
    return {
        "gains": [float(g) for g in gains],
        "rows": rows,
        "worst_signed_or_winding": worst,
        "every_rule_gives_zero": bool(worst == 0),
        "the_negative_result_survives_the_choice": bool(worst == 0),
    }


def winding_curve(bulk: BulkAnnulus, turns: int = 1, n: int = 1441,
                  r_start: Optional[float] = None
                  ) -> Tuple[np.ndarray, np.ndarray]:
    """A curve that *does* wind — for comparison, since the wave never does.

    In the seam's own coordinate it is a straight ramp advancing exactly
    ``turns`` periods as ``σ`` goes once around.  On a conformal seam that is a
    logarithmic spiral: it comes back to the same point of the quotient having
    been **magnified** by ``(R_outer/R_inner)^turns``.  On a translate seam it
    comes back shifted by ``turns × gap`` with no magnification at all.

    This is what the winding number looks like when it is not zero, and it is
    not a graph: no single-valued ``r = f(σ)`` can do it.
    """
    sigma = np.linspace(-math.pi, math.pi, int(n))
    w0 = float(bulk.to_w(np.array([r_start or bulk.r_mid]))[0])
    w = w0 + turns * bulk.period * (sigma + math.pi) / TWO_PI
    return sigma, np.asarray(bulk.from_w(w), dtype=float)


def measure_winding_is_a_magnification(turns: int = 1) -> Dict[str, object]:
    """On a conformal seam the winding number is visible as a scale factor.

    Which is the real answer to "the scaling is a choice": the choice decides
    what a winding number would *look* like.  Glue by translation and a wound
    curve returns displaced; glue conformally and it returns magnified by
    ``(R_outer/R_inner)^w``.  Either way the physical wave has ``w = 0``, so it
    returns to exactly itself — and that is not a coincidence, it is what lets a
    scale-changing seam be consistent for a single-valued height at all.

    Whether the factor is a genuine *scale* is decided by starting the curve at
    several radii.  On a conformal seam the magnification is the same for all of
    them; on a translate seam the ratio drifts with the starting radius, because
    a fixed shift is not a scale.
    """
    starts = (0.80, 1.00, 1.20)
    rows = []
    for mode in ("translate", "conformal"):
        b = BulkAnnulus(mode=mode)
        mags = []
        for r0 in starts:
            _, r = winding_curve(b, turns=turns, r_start=r0)
            _, sheet = b.wrap(r)
            mags.append(float(r[-1] / r[0]))
            climbed = int(sheet[-1] - sheet[0])
        rows.append({
            "seam": mode,
            "start_radii": list(starts),
            "magnifications": mags,
            "magnification_spread": max(mags) - min(mags),
            "displacement": float(b.period * turns) if mode == "translate" else None,
            "sheets_climbed": climbed,
            "predicted_magnification": (
                None if mode == "translate"
                else (b.r_outer / b.r_inner) ** turns),
        })
    by = {r["seam"]: r for r in rows}
    conf, tran = by["conformal"], by["translate"]
    return {
        "turns": int(turns),
        "rows": rows,
        "conformal_magnification": conf["magnifications"][1],
        "conformal_spread": conf["magnification_spread"],
        "translate_spread": tran["magnification_spread"],
        "a_wound_curve_climbs_the_sheets": bool(
            all(r["sheets_climbed"] == turns for r in rows)),
        "conformal_turns_winding_into_a_scale": bool(
            conf["magnification_spread"] < 1e-9
            and abs(conf["magnifications"][1]
                    - conf["predicted_magnification"]) < 1e-9),
        "translate_gives_no_scale_at_all": bool(
            tran["magnification_spread"] > 0.1),
    }
