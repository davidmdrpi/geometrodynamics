"""
Can the converging ring reach across the gap — and can it ever intersect?

THE RING IS REAL
────────────────
The scalar wave is not only its focal pulse.  A ring leaves the source, thins to
a minimum at the equator, and then **grows** as it converges on the antipode —
classic ``1/√(sin d)`` focusing.  Measured here: ``0.158`` at the equator rising
to ``0.925`` at the focus, a factor of ``5.9``, all of it *before* the final
pulse.  So it is fair to ask whether that ring, rather than the pulse, is what
first reaches across the vacuole.

TWO QUESTIONS THAT LOOK LIKE ONE
────────────────────────────────
**Can it reach across?**  Yes, and easily.  Spanning is a race between the
excursion ``ε·u`` and the half-gap ``δ``, so it is bought by raising the gain or
shrinking ``δ`` — either works, and the threshold is exactly ``δ / max|u|``.
Shrinking the gap also buys *dwell*: the fraction of the wave's life spent
spanning grows as the gap closes, so the reaching stops being an instant at the
focus and becomes a sustained state.

**Can it intersect?**  No — and no amount of either knob changes that.  A curve
drawn as ``r = f(σ)`` with ``f`` single-valued assigns exactly one radius to
each direction, so it is *embedded* by construction.  Two of its points can
never occupy the same place.  This is the same obstruction that made the winding
number zero, seen from the side: a graph cannot cross itself any more than it
can wind.  It is checked here by an actual segment-intersection sweep over gap
and gain rather than by assertion.

WHAT DOES INTERSECT
───────────────────
Let the material move **sideways**.  If each point carries a tangential
displacement as well as a radial one,

    σ(σ₀) = σ₀ + λ ε ∂_σ u ,      r(σ₀) = R_mid + ε u ,

then the map ``σ₀ ↦ σ`` can fold.  It folds where

    ∂σ/∂σ₀ = 1 + λ ε ∂²_σ u  <  0 ,

so the threshold is ``λ ε = 1 / max(−∂²_σ u)`` — set by the **curvature of the
front**, not by the gap at all.  Past it the curve is multivalued in ``σ``, stops
being a graph, and self-intersects.

That is the honest reconciliation.  The ring is exactly the right place to look:
``∂²_σ u`` is largest at the steep converging front, so the fold appears there
first and at the moment of tightest focus.  What the ring cannot do is fold
while it is still a height — the freedom it needs is tangential, and the gap is
the wrong knob for buying it.
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from geometrodynamics.viz.circle_slice import (
    ANTIPODAL_TIME,
    RETURN_TIME,
    TWO_PI,
    BulkAnnulus,
    CircleSlice,
)
from geometrodynamics.viz.throat_wavefront import BareSphereSim
from geometrodynamics.viz.warped_sphere import NestedShells


# ════════════════════════════════════════════════════════════════════════════
# AN HONEST INTERSECTION TEST
# ════════════════════════════════════════════════════════════════════════════
def _orient(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
    return ((b[..., 0] - a[..., 0]) * (c[..., 1] - a[..., 1])
            - (b[..., 1] - a[..., 1]) * (c[..., 0] - a[..., 0]))


def self_intersections(xy: np.ndarray, closed: bool = True) -> int:
    """Count crossing pairs of segments of a polyline — properly, not by proxy.

    All pairs are tested with orientation predicates, skipping neighbours (which
    share an endpoint by construction).  This is deliberately the blunt
    ``O(n²)`` version: the claim being checked is a negative one, so the test
    should have no cleverness in it to be wrong about.
    """
    p = np.asarray(xy, dtype=float)
    if closed and not np.allclose(p[0], p[-1]):
        p = np.vstack([p, p[:1]])
    a, b = p[:-1], p[1:]
    n = len(a)
    i, j = np.triu_indices(n, k=2)
    if closed:                       # first and last segments are neighbours too
        keep = ~((i == 0) & (j == n - 1))
        i, j = i[keep], j[keep]
    a1, b1, a2, b2 = a[i], b[i], a[j], b[j]
    d1 = _orient(a1, b1, a2) * _orient(a1, b1, b2)
    d2 = _orient(a2, b2, a1) * _orient(a2, b2, b1)
    return int(np.count_nonzero((d1 < 0.0) & (d2 < 0.0)))


# ════════════════════════════════════════════════════════════════════════════
# THE RING
# ════════════════════════════════════════════════════════════════════════════
def measure_the_ring_grows_as_it_converges(pulse_width: float = 0.18,
                                           n_radial: int = 1800,
                                           frames: int = 300
                                           ) -> Dict[str, object]:
    """The ring thins to the equator and then grows — before any focal pulse.

    Its amplitude tracks ``1/√(sin d)``, the geometric factor for a ring whose
    circumference is closing.  Reported against that law rather than merely
    asserted to follow it.
    """
    sim = BareSphereSim(n_theta=8, n_phi=8, pulse_width=pulse_width,
                        n_radial=n_radial)
    sim.reset()
    d = np.linspace(0.0, math.pi, 1500)
    rows = []
    for i in range(frames):
        sim.advance_to((i + 1) * ANTIPODAL_TIME / frames)
        u = np.abs(sim.field_at_distance(d))
        j = int(np.argmax(u))
        rows.append({"t": sim.t, "front": float(d[j]), "height": float(u[j])})

    # The minimum is taken on the OUTBOUND leg only.  Sweeping the whole
    # history instead finds the momentary null just past the focus, which is
    # the Gouy phase passing through zero and has nothing to do with the ring.
    outbound = [r for r in rows if 0.2 < r["t"] < 0.85 * ANTIPODAL_TIME]
    equator = min(outbound, key=lambda r: r["height"])
    peak = max(rows, key=lambda r: r["height"])

    # ...against the 1/√(sin d) law, normalised at the equatorial minimum.
    # Restricted to sin d > 0.25: inside a pulse width of the antipode the
    # ring is no longer thin compared with its own radius and the law stops
    # applying, which is reported rather than fitted away.
    checks = []
    for r in rows:
        s = math.sin(r["front"])
        if 0.25 < s < 0.98 and r["front"] > equator["front"] and r["t"] < ANTIPODAL_TIME:
            predicted = equator["height"] * math.sqrt(
                math.sin(equator["front"]) / s)
            checks.append(r["height"] / predicted)
    return {
        "pulse_width": pulse_width,
        "equator_distance": equator["front"],
        "equator_height": equator["height"],
        "peak_height": peak["height"],
        "peak_time": peak["t"],
        "peak_distance": peak["front"],
        "growth_from_the_equator": peak["height"] / equator["height"],
        "law_ratio_mean": float(np.mean(checks)) if checks else float("nan"),
        "law_ratio_spread": (float(np.max(checks) - np.min(checks))
                             if checks else float("nan")),
        "launch_height": rows[0]["height"],
        "the_ring_thins_then_grows": bool(
            equator["height"] < rows[0]["height"]
            and peak["height"] > 2.0 * equator["height"]),
        "follows_one_over_root_sin": bool(
            checks and abs(float(np.mean(checks)) - 1.0) < 0.35),
    }


# ════════════════════════════════════════════════════════════════════════════
# REACHING ACROSS
# ════════════════════════════════════════════════════════════════════════════
def measure_when_the_ring_spans_the_gap(
        deltas: Sequence[float] = (0.26, 0.16, 0.09, 0.05),
        gains: Sequence[float] = (0.05, 0.10, 0.20, 0.40),
        pulse_width: float = 0.18, frames: int = 240) -> Dict[str, object]:
    """Reaching across is easy, and both knobs buy it.

    The threshold is exactly ``δ / max|u|``: raise the gain or shrink the gap and
    the excursion spans.  Shrinking the gap buys something extra — **dwell**, the
    fraction of the wave's life spent spanning — so the reaching stops being an
    instant at the focus and becomes a sustained state.

    Also recorded: whether the *first* spanning happens before the focal pulse,
    i.e. whether the converging ring gets there first.  It does, once the gap is
    small enough that spanning starts well before ``t = π``.
    """
    sim = BareSphereSim(n_theta=8, n_phi=8, pulse_width=pulse_width,
                        n_radial=1200)
    d = np.linspace(0.0, math.pi, 1200)

    sim.reset()
    peak = 0.0
    for i in range(frames):
        sim.advance_to((i + 1) * RETURN_TIME / frames)
        peak = max(peak, float(np.max(np.abs(sim.field_at_distance(d)))))

    rows = []
    for delta in deltas:
        for gain in gains:
            sim.reset()
            spanning, first_t, first_d = 0, None, None
            conv_t, conv_d = None, None
            for i in range(frames):
                t = (i + 1) * RETURN_TIME / frames
                sim.advance_to(t)
                u = sim.field_at_distance(d)
                j = int(np.argmax(np.abs(u)))
                if gain * abs(float(u[j])) >= delta:
                    spanning += 1
                    if first_t is None:
                        first_t, first_d = sim.t, float(d[j])
                    # ...and separately, on the CONVERGING leg only.  The
                    # launch spans trivially — the source is the tallest thing
                    # in the run — so it has to be excluded to answer the
                    # question that was actually asked about the ring.
                    if conv_t is None and 0.55 * ANTIPODAL_TIME < sim.t < ANTIPODAL_TIME:
                        conv_t, conv_d = sim.t, float(d[j])
            rows.append({
                "delta": float(delta), "gain": float(gain),
                "threshold_gain": float(delta) / max(peak, 1e-30),
                "spans": bool(spanning > 0),
                "dwell_fraction": spanning / frames,
                "first_spanning_time": first_t,
                "first_spanning_distance": first_d,
                "converging_span_time": conv_t,
                "converging_span_distance": conv_d,
                "lead_before_the_focus": (None if conv_t is None
                                          else ANTIPODAL_TIME - conv_t),
            })
    spanning_rows = [r for r in rows if r["spans"]]
    leads = [r["lead_before_the_focus"] for r in rows
             if r["lead_before_the_focus"] is not None]
    return {
        "run_peak": peak,
        "rows": rows,
        "n_spanning": len(spanning_rows),
        "max_dwell_fraction": max((r["dwell_fraction"] for r in rows),
                                  default=0.0),
        "longest_lead_before_the_focus": max(leads, default=0.0),
        "threshold_is_delta_over_peak": bool(all(
            r["spans"] == (r["gain"] >= r["threshold_gain"] - 1e-9)
            for r in rows)),
        "shrinking_the_gap_buys_dwell": bool(
            max((r["dwell_fraction"] for r in rows
                 if r["delta"] == min(deltas)), default=0.0)
            > max((r["dwell_fraction"] for r in rows
                   if r["delta"] == max(deltas)), default=0.0)),
        "shrinking_the_gap_buys_lead": bool(
            max((r["lead_before_the_focus"] or 0.0 for r in rows
                 if r["delta"] == min(deltas)), default=0.0)
            > max((r["lead_before_the_focus"] or 0.0 for r in rows
                   if r["delta"] == max(deltas)), default=0.0)),
        "the_converging_ring_spans_well_before_the_pulse": bool(
            max(leads, default=0.0) > 0.5),
    }


# ════════════════════════════════════════════════════════════════════════════
# ...AND STILL NEVER INTERSECTING
# ════════════════════════════════════════════════════════════════════════════
def measure_a_graph_never_intersects(
        deltas: Sequence[float] = (0.26, 0.12, 0.05),
        gains: Sequence[float] = (0.2, 1.0, 5.0, 40.0),
        pulse_width: float = 0.18, frames: int = 90,
        n_sigma: int = 1201) -> Dict[str, object]:
    """Swept over gap and gain, with a real intersection test, the answer is no.

    A single-valued ``f`` puts exactly one radius at each direction, so the drawn
    curve is embedded and two of its points can never coincide.  Shrinking the
    gap and raising the gain buy spanning and seam crossings; they never buy a
    crossing of the curve with itself.
    """
    rows = []
    worst = 0
    for delta in deltas:
        shells = NestedShells(r_mid=1.0, delta=float(delta))
        for gain in gains:
            s = CircleSlice(bulk=BulkAnnulus(shells, mode="conformal"),
                            radial_law="multiplicative",
                            pulse_width=pulse_width, n_sigma=n_sigma,
                            gain=float(gain))
            s.reset()
            hits, most_cross = 0, 0
            for i in range(frames):
                s.advance_to((i + 1) * RETURN_TIME / frames)
                hits = max(hits, self_intersections(s.points(gain=gain)))
                most_cross = max(most_cross,
                                 s.seam_crossings(gain=gain)["unsigned"])
            worst = max(worst, hits)
            rows.append({"delta": float(delta), "gain": float(gain),
                         "self_intersections": hits,
                         "seam_crossings": most_cross})
    return {
        "rows": rows,
        "worst_self_intersections": worst,
        "most_seam_crossings": max(r["seam_crossings"] for r in rows),
        "a_graph_is_embedded": bool(worst == 0),
        "and_it_was_genuinely_driven": bool(
            max(r["seam_crossings"] for r in rows) > 0),
    }


def measure_the_intersection_test_can_see_one() -> Dict[str, object]:
    """A control, so the negative result above is not just a broken detector.

    A limaçon ``r = a + b cos σ`` with ``b > a`` has a genuine inner loop and
    crosses itself; a circle does not.  Both go through the same counter.
    """
    sig = np.linspace(0.0, TWO_PI, 1441)

    def polar(r):
        return np.stack([r * np.cos(sig), r * np.sin(sig)], axis=-1)

    circle = self_intersections(polar(np.ones_like(sig)))
    limacon = self_intersections(polar(0.4 + 1.0 * np.cos(sig)))
    # ...offset off the sample grid: the lemniscate's crossing sits exactly on a
    # node at σ = 0, π, and the orientation test is strict, so an exactly
    # endpoint-on-segment touch is deliberately not counted as a crossing
    off = sig + 0.5 * float(sig[1] - sig[0])
    figure8 = self_intersections(np.stack(
        [np.sin(2.0 * off), np.sin(off)], axis=-1))
    # a fold of the kind this module is actually hunting for
    folded = self_intersections(polar(1.0 + 0.0 * sig)
                                + 0.35 * np.stack([np.cos(3.0 * sig),
                                                   np.zeros_like(sig)], axis=-1))
    return {
        "circle": circle,
        "limacon_with_an_inner_loop": limacon,
        "figure_eight": figure8,
        "a_folded_loop": folded,
        "the_detector_works": bool(circle == 0 and limacon > 0
                                   and figure8 > 0 and folded > 0),
    }


# ════════════════════════════════════════════════════════════════════════════
# WHAT DOES INTERSECT — LETTING THE MATERIAL MOVE SIDEWAYS
# ════════════════════════════════════════════════════════════════════════════
class LagrangianSlice:
    """The slice as *material points*, free to move tangentially as well.

    Each point is labelled by its rest position ``σ₀`` and carries

        ``σ = σ₀ + λ ε ∂_σ u``     ``r = R_mid + ε u``   (or ``R_mid·e^{εu}``)

    ``λ = 0`` is the height field of v46 — a graph, embedded forever.  ``λ > 0``
    lets the map ``σ₀ ↦ σ`` fold, and a folded map is not a graph.

    The mixing ``λ`` is a modelling choice and is reported as one: it says how
    much of the wave's displacement is along the surface rather than across it.
    Nothing here derives it from the scalar equation.
    """

    def __init__(self, sim: Optional[BareSphereSim] = None,
                 shells: Optional[NestedShells] = None,
                 mixing: float = 1.0, gain: float = 0.2,
                 n_sigma: int = 2001, pulse_width: float = 0.18,
                 n_radial: int = 2400,
                 radial_law: str = "multiplicative") -> None:
        self.sim = sim if sim is not None else BareSphereSim(
            n_theta=8, n_phi=8, pulse_width=pulse_width, n_radial=n_radial)
        self.shells = shells or NestedShells()
        self.mixing = float(mixing)
        self.gain = float(gain)
        self.radial_law = radial_law
        self.n_sigma = int(n_sigma)
        # endpoint EXCLUDED: σ = −π and σ = +π are the same point, so the
        # sample set must not contain both.  Keeping both and differencing
        # one-sidedly at the ends makes ∂_σu non-periodic, which closes the
        # drawn curve with a chord that crosses it — a spurious intersection
        # reported at times when the map has not folded at all.
        self.sigma0 = np.linspace(-math.pi, math.pi, self.n_sigma,
                                  endpoint=False)
        self.dsigma = float(TWO_PI / self.n_sigma)

    @property
    def t(self) -> float:
        return self.sim.t

    def reset(self) -> None:
        self.sim.reset()

    def advance_to(self, t: float) -> None:
        self.sim.advance_to(t)

    # ── the field and its σ-derivatives ─────────────────────────────────────
    def field(self) -> np.ndarray:
        return self.sim.field_at_distance(np.abs(self.sigma0))

    def _du(self) -> np.ndarray:
        """Periodic central difference — exactly odd, so the curve closes."""
        u = self.field()
        return (np.roll(u, -1) - np.roll(u, 1)) / (2.0 * self.dsigma)

    def _d2u(self) -> np.ndarray:
        u = self.field()
        return (np.roll(u, -1) - 2.0 * u + np.roll(u, 1)) / self.dsigma ** 2

    # ── the drawn curve ─────────────────────────────────────────────────────
    def sigma(self, gain: Optional[float] = None,
              mixing: Optional[float] = None) -> np.ndarray:
        eps = self.gain if gain is None else float(gain)
        lam = self.mixing if mixing is None else float(mixing)
        return self.sigma0 + lam * eps * self._du()

    def radius(self, gain: Optional[float] = None) -> np.ndarray:
        eps = self.gain if gain is None else float(gain)
        u = self.field()
        if self.radial_law == "additive":
            return self.shells.r_mid + eps * u
        return self.shells.r_mid * np.exp(eps * u)

    def points(self, gain: Optional[float] = None,
               mixing: Optional[float] = None) -> np.ndarray:
        sig = self.sigma(gain=gain, mixing=mixing)
        r = self.radius(gain=gain)
        return np.stack([r * np.cos(sig), r * np.sin(sig)], axis=-1)

    # ── the fold ────────────────────────────────────────────────────────────
    def jacobian(self, gain: Optional[float] = None,
                 mixing: Optional[float] = None) -> np.ndarray:
        """``∂σ/∂σ₀`` — the map folds where this goes negative."""
        eps = self.gain if gain is None else float(gain)
        lam = self.mixing if mixing is None else float(mixing)
        return 1.0 + lam * eps * self._d2u()

    def is_folded(self, gain: Optional[float] = None,
                  mixing: Optional[float] = None) -> bool:
        return bool(np.any(self.jacobian(gain=gain, mixing=mixing) < 0.0))

    def fold_margin(self, gain: Optional[float] = None,
                    mixing: Optional[float] = None) -> float:
        return float(np.min(self.jacobian(gain=gain, mixing=mixing)))

    def predicted_fold_product(self) -> float:
        """The ``λε`` at which this instant first folds: ``1 / max(−∂²_σ u)``."""
        worst = float(np.max(-self._d2u()))
        return float("inf") if worst <= 0.0 else 1.0 / worst


def measure_the_fold_threshold(mixing: float = 1.0, pulse_width: float = 0.18,
                               frames: int = 200, n_sigma: int = 2001,
                               n_radial: int = 2400) -> Dict[str, object]:
    """Where the curve stops being a graph, predicted and then found.

    The prediction is ``λ ε = 1 / max(−∂²_σ u)`` — pure front curvature.  The
    measurement bisects on it and then checks the drawn curve with the same
    segment test used for the negative result, so the two are comparable.

    Folding is **necessary but not sufficient** for a crossing, and that is
    reported rather than glossed: past the threshold the map is multivalued in
    ``σ``, but the two branches only intersect *in the plane* if their radii
    also overlap.  What never happens is the converse — a crossing without a
    fold — because that is the embedded-graph theorem again.
    """
    s = LagrangianSlice(mixing=mixing, pulse_width=pulse_width,
                        n_sigma=n_sigma, n_radial=n_radial)

    s.reset()
    best = {"product": math.inf, "time": 0.0, "distance": 0.0}
    for i in range(frames):
        s.advance_to((i + 1) * ANTIPODAL_TIME / frames)
        p = s.predicted_fold_product()
        if p < best["product"]:
            j = int(np.argmax(-s._d2u()))
            best = {"product": p, "time": s.t,
                    "distance": float(abs(s.sigma0[j]))}
    predicted_gain = best["product"] / mixing

    def folds_somewhere(gain: float) -> Tuple[bool, int]:
        s.reset()
        hits = 0
        folded = False
        for i in range(frames):
            s.advance_to((i + 1) * ANTIPODAL_TIME / frames)
            if s.is_folded(gain=gain):
                folded = True
                hits = max(hits, self_intersections(s.points(gain=gain),
                                                    closed=False))
        return folded, hits

    lo, hi = 0.0, 4.0 * predicted_gain
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if folds_somewhere(mid)[0]:
            hi = mid
        else:
            lo = mid
    measured_gain = 0.5 * (lo + hi)
    _, hits_above = folds_somewhere(1.6 * predicted_gain)
    _, hits_below = folds_somewhere(0.6 * predicted_gain)

    return {
        "mixing": float(mixing),
        "predicted_fold_product": best["product"],
        "predicted_gain": predicted_gain,
        "measured_gain": measured_gain,
        "relative_error": abs(measured_gain - predicted_gain)
        / max(predicted_gain, 1e-30),
        "folds_first_at_time": best["time"],
        "folds_first_at_distance": best["distance"],
        "distance_to_the_antipode": math.pi - best["distance"],
        "self_intersections_above": hits_above,
        "self_intersections_below": hits_below,
        "the_threshold_is_front_curvature": bool(
            abs(measured_gain - predicted_gain) < 0.05 * predicted_gain),
        "past_it_the_curve_self_intersects": bool(hits_above > 0),
        "below_it_the_curve_does_not": bool(hits_below == 0),
        "it_folds_at_the_converging_ring": bool(
            best["distance"] > 0.5 * math.pi),
    }


def measure_folding_ignores_the_gap(
        deltas: Sequence[float] = (0.26, 0.12, 0.05),
        mixing: float = 1.0, pulse_width: float = 0.18,
        frames: int = 120) -> Dict[str, object]:
    """The two knobs are not the same knob.

    The gap sets when the wave *reaches across*; the front's curvature sets when
    it *folds*.  Changing the gap moves the first and leaves the second exactly
    where it was — which is why shrinking the vacuole can never buy an
    intersection.
    """
    rows = []
    for delta in deltas:
        s = LagrangianSlice(shells=NestedShells(r_mid=1.0, delta=float(delta)),
                            mixing=mixing, pulse_width=pulse_width,
                            n_sigma=1601)
        s.reset()
        product = math.inf
        for i in range(frames):
            s.advance_to((i + 1) * ANTIPODAL_TIME / frames)
            product = min(product, s.predicted_fold_product())
        rows.append({"delta": float(delta),
                     "span_threshold_gain": float(delta),   # ÷ peak, below
                     "fold_threshold_product": product})
    products = [r["fold_threshold_product"] for r in rows]
    return {
        "mixing": float(mixing),
        "rows": rows,
        "fold_threshold_spread": max(products) - min(products),
        "fold_threshold_is_gap_independent": bool(
            max(products) - min(products) < 1e-6 * max(products)),
        "span_threshold_scales_with_the_gap": bool(
            rows[0]["span_threshold_gain"] > rows[-1]["span_threshold_gain"]),
    }


def measure_the_fold_threshold_converges(
        pulse_width: float = 0.18,
        grids: Sequence[Tuple[int, int]] = ((1001, 1200), (2001, 2400),
                                            (4001, 4800)),
        frames: int = 120) -> Dict[str, object]:
    """The threshold is a number, not a resolution — checked before it is used.

    ``∂²_σ u`` at a focusing front is exactly the quantity a coarse grid gets
    wrong, so this refines the angular sampling and the radial solve **together**
    and watches the threshold settle.  Refining ``σ`` alone does not converge:
    it just samples interpolation noise in a solve that has not kept up, which
    is why the two are tied here.
    """
    rows = []
    for n_sigma, n_radial in grids:
        s = LagrangianSlice(mixing=1.0, pulse_width=pulse_width,
                            n_sigma=n_sigma, n_radial=n_radial)
        s.reset()
        best = math.inf
        for i in range(frames):
            s.advance_to((i + 1) * ANTIPODAL_TIME / frames)
            best = min(best, s.predicted_fold_product())
        rows.append({"n_sigma": n_sigma, "n_radial": n_radial,
                     "max_negative_curvature": 1.0 / best,
                     "fold_product": best})
    products = [r["fold_product"] for r in rows]
    drift = abs(products[-1] - products[-2]) / max(products[-1], 1e-30)
    return {
        "pulse_width": pulse_width,
        "rows": rows,
        "last_step_drift": drift,
        "converged_fold_product": products[-1],
        "it_converges": bool(drift < 0.02),
    }


def measure_the_fold_threshold_scales_with_the_pulse(
        widths: Sequence[float] = (0.36, 0.24, 0.18, 0.12),
        frames: int = 120, n_sigma: int = 2001,
        n_radial: int = 2400) -> Dict[str, object]:
    """Narrower fronts fold sooner, as the square of their width.

    ``∂²_σ u`` at the focus goes like ``1/w²`` — the front's height is set by the
    focusing and its curvature by its own width — so the threshold ``λε`` scales
    as ``w²``.  The constant is reported so the law can be checked rather than
    admired.
    """
    rows = []
    for w in widths:
        s = LagrangianSlice(mixing=1.0, pulse_width=float(w), n_sigma=n_sigma,
                            n_radial=n_radial)
        s.reset()
        best = math.inf
        for i in range(frames):
            s.advance_to((i + 1) * ANTIPODAL_TIME / frames)
            best = min(best, s.predicted_fold_product())
        rows.append({"pulse_width": float(w), "fold_product": best,
                     "product_over_width_squared": best / (float(w) ** 2)})
    ks = [r["product_over_width_squared"] for r in rows]
    mean = sum(ks) / len(ks)
    return {
        "rows": rows,
        "mean_constant": mean,
        "spread": max(ks) - min(ks),
        "narrower_folds_sooner": bool(rows[-1]["fold_product"]
                                      < rows[0]["fold_product"]),
        "scales_as_the_width_squared": bool(
            max(ks) - min(ks) < 0.12 * mean),
    }


def measure_folding_is_necessary_for_crossing(
        mixing: float = 1.0, pulse_width: float = 0.18,
        multiples: Sequence[float] = (0.6, 1.5, 3.0),
        frames: int = 60, n_sigma: int = 801,
        n_radial: int = 1200) -> Dict[str, object]:
    """A crossing always comes with a fold; a fold does not always cross.

    The asymmetry is the point.  ``crossing without fold`` is the embedded-graph
    theorem restated, so it must be zero at every drive.  ``fold without
    crossing`` is allowed and does happen: the map can be multivalued in ``σ``
    while the two branches stay radially apart, so they never meet in the plane.
    """
    s = LagrangianSlice(mixing=mixing, pulse_width=pulse_width,
                        n_sigma=n_sigma, n_radial=n_radial)
    s.reset()
    threshold = math.inf
    for i in range(frames):
        s.advance_to((i + 1) * ANTIPODAL_TIME / frames)
        threshold = min(threshold, s.predicted_fold_product())

    rows = []
    fold_no_cross = cross_no_fold = 0
    for m in multiples:
        gain = float(m) * threshold / mixing
        s.reset()
        folded = crossed = 0
        for i in range(frames):
            s.advance_to((i + 1) * 1.12 * ANTIPODAL_TIME / frames)
            f = s.is_folded(gain=gain)
            h = self_intersections(s.points(gain=gain))
            folded += int(f)
            crossed += int(h > 0)
            fold_no_cross += int(f and h == 0)
            cross_no_fold += int(h > 0 and not f)
        rows.append({"multiple": float(m), "gain": gain,
                     "folded_frames": folded, "crossing_frames": crossed})
    return {
        "threshold_product": threshold,
        "rows": rows,
        "fold_without_crossing": fold_no_cross,
        "crossing_without_fold": cross_no_fold,
        "a_crossing_always_folds": bool(cross_no_fold == 0),
        "a_fold_need_not_cross": bool(fold_no_cross > 0),
        "nothing_crosses_below_threshold": bool(rows[0]["crossing_frames"] == 0),
    }
