"""
Draw the wave as vectors normal to the surface, and they intersect.

THE ABSTRACTION THAT WAS MISSING
────────────────────────────────
Every earlier round drew the field as the *tip* of a radial displacement — the
graph ``r = f(σ)``.  A graph is embedded by construction, which is why it could
not wind (``circle_slice``), could not self-intersect at any gap or gain
(``slice_folding``), and needed an invented tangential freedom before anything
would fold.

Draw the **vectors** instead of their tips and the obstruction evaporates, for a
reason that is entirely classical: *neighbouring normals to a curve meet at its
centre of curvature*.  The family of normals has an envelope — the evolute — and
a normal of length ``L`` crosses its neighbours as soon as ``L`` exceeds the
local radius of curvature ``ρ = 1/κ``.  Nothing has to be added by hand.

WHAT THAT BUYS
──────────────
* **Intersections, at last.**  On the same wave, at the same time, the graph of
  the tips gives ``0`` self-intersections while the normal field gives
  hundreds.
* **A threshold with meaning.**  It is the radius of curvature of the deformed
  surface, and the converging ring drives it *down*: measured ``0.1408`` at
  mid-flight, ``0.1087`` converging, ``0.0540`` at the focus.  The wave sharpens
  its own surface until its normals cross.
* **The gap matters again.**  ``slice_folding`` established that shrinking the
  vacuole could never buy an intersection, because the fold threshold did not
  know the gap existed.  For normals it does: the vector length is what spans
  the gap, so ``L`` and ``δ`` are the same knob.

AND THE RESET IS A SECOND MECHANISM
───────────────────────────────────
A normal long enough to leave through ``R_outer`` re-enters at ``R_inner`` — at
the angle where it exited, continuing in the same direction.  That re-entered
stub starts deep inside the annulus and shoots outward across everything, so it
crosses vectors it could never have reached.  Measured at the focus with
``L = 0.35``: ``306`` crossings among the normals themselves, ``398`` once the
reset is included.  The two mechanisms are separable and both are real.
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from geometrodynamics.viz.circle_slice import (
    ANTIPODAL_TIME,
    RETURN_TIME,
    BulkAnnulus,
    CircleSlice,
)
from geometrodynamics.viz.warped_sphere import NestedShells


# ════════════════════════════════════════════════════════════════════════════
# SEGMENT CROSSINGS
# ════════════════════════════════════════════════════════════════════════════
def _cross(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]


def segment_crossings(a: np.ndarray, b: np.ndarray) -> int:
    """Count transversal crossings among a bundle of open segments.

    Unlike the polyline test in ``slice_folding`` these segments share no
    endpoints, so every pair is a fair candidate and nothing has to be skipped.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    total = 0
    for i in range(len(a) - 1):
        u = b[i] - a[i]
        v = b[i + 1:] - a[i + 1:]
        d1 = _cross(u, a[i + 1:] - a[i]) * _cross(u, b[i + 1:] - a[i])
        d2 = _cross(v, a[i] - a[i + 1:]) * _cross(v, b[i] - a[i + 1:])
        total += int(np.count_nonzero((d1 < 0.0) & (d2 < 0.0)))
    return total


# ════════════════════════════════════════════════════════════════════════════
# THE NORMAL FIELD OF A DEFORMED SLICE
# ════════════════════════════════════════════════════════════════════════════
class NormalField:
    """The slice's outward normals, drawn as vectors rather than as a graph.

    ``length`` is the drawn vector length.  Nothing else here is a display
    choice: the directions and the curvature are the surface's own.
    """

    def __init__(self, slice_: Optional[CircleSlice] = None,
                 delta: float = 0.26, pulse_width: float = 0.18,
                 n_sigma: int = 4001, n_radial: int = 2400,
                 gain: float = 0.30, stride: int = 8) -> None:
        shells = NestedShells(r_mid=1.0, delta=float(delta))
        self.shells = shells
        # n_radial is matched to n_sigma deliberately.  κ needs a second
        # derivative of the field, so refining the angular sampling against a
        # coarse radial solve does not converge — it just samples interpolation
        # noise, the same trap ``slice_folding`` documents for the fold
        # threshold.  ``measure_the_curvature_converges`` holds this down.
        self.slice = slice_ or CircleSlice(
            bulk=BulkAnnulus(shells, mode="conformal"),
            radial_law="multiplicative", pulse_width=pulse_width,
            n_sigma=n_sigma, n_radial=n_radial, gain=gain)
        self.gain = float(gain)
        self.stride = int(stride)

    # ── clock ───────────────────────────────────────────────────────────────
    @property
    def t(self) -> float:
        return self.slice.t

    def reset(self) -> None:
        self.slice.reset()

    def advance_to(self, t: float) -> None:
        self.slice.advance_to(t)

    # ── geometry of the drawn curve ─────────────────────────────────────────
    def _radius(self) -> np.ndarray:
        return self.slice.radius(gain=self.gain)

    def frame(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """``(X, n̂, κ)`` — position, outward unit normal, signed curvature.

        Derivatives are periodic (``np.roll``), because ``σ = ±π`` is one point.
        """
        sig = self.slice.sigma
        r = self._radius()
        ds = float(sig[1] - sig[0])
        rp = (np.roll(r, -1) - np.roll(r, 1)) / (2.0 * ds)
        rpp = (np.roll(r, -1) - 2.0 * r + np.roll(r, 1)) / ds ** 2
        cos, sin = np.cos(sig), np.sin(sig)
        X = np.stack([r * cos, r * sin], axis=-1)
        Xp = np.stack([rp * cos - r * sin, rp * sin + r * cos], axis=-1)
        T = Xp / np.linalg.norm(Xp, axis=1)[:, None]
        N = np.stack([T[:, 1], -T[:, 0]], axis=-1)
        if float(np.mean(np.sum(N * X, axis=1))) < 0.0:
            N = -N                              # keep it pointing outward
        kappa = ((r ** 2 + 2.0 * rp ** 2 - r * rpp)
                 / (r ** 2 + rp ** 2) ** 1.5)
        return X, N, kappa

    def radius_of_curvature(self) -> np.ndarray:
        """``1/κ`` where the normals converge, ``inf`` where they diverge."""
        _, _, kappa = self.frame()
        return np.where(kappa > 1e-9, 1.0 / np.maximum(kappa, 1e-9), np.inf)

    def envelope_distance(self) -> float:
        """The smallest ``L`` at which *some* pair of normals must have met."""
        return float(np.min(self.radius_of_curvature()))

    # ── the drawn vectors ───────────────────────────────────────────────────
    def _sample(self) -> np.ndarray:
        """Strided indices with the duplicated endpoint dropped.

        ``CircleSlice.sigma`` deliberately keeps both ``σ = −π`` and ``σ = +π``
        so the drawn curve closes, but they are the *same point*.  Sampling both
        puts two coincident normals in the bundle, and the orientation test
        scores that as a crossing at every length, including ``1e-04`` — a
        length at which nothing can possibly have met.
        """
        return np.arange(0, self.slice.n_sigma - 1, self.stride)

    def vectors(self, length: float) -> Tuple[np.ndarray, np.ndarray]:
        X, N, _ = self.frame()
        i = self._sample()
        return X[i], X[i] + float(length) * N[i]

    def _exit_parameter(self, P: np.ndarray, D: np.ndarray,
                        radius: float) -> np.ndarray:
        """Where ``P + tD`` first leaves the circle of the given radius."""
        a = np.sum(D * D, axis=1)
        b = 2.0 * np.sum(P * D, axis=1)
        c = np.sum(P * P, axis=1) - radius ** 2
        disc = b * b - 4.0 * a * c
        t = np.full(len(P), np.inf)
        ok = disc > 0.0
        t[ok] = (-b[ok] + np.sqrt(disc[ok])) / (2.0 * a[ok])
        return t

    def vectors_with_reset(self, length: float
                           ) -> Tuple[np.ndarray, np.ndarray, np.ndarray,
                                      np.ndarray]:
        """Vectors clipped at ``R_outer``, plus the stubs that re-enter.

        A vector that leaves through the outer shell comes back at ``R_inner``
        *at the angle where it left*, continuing in the same direction — the
        same crossing rule the slice itself uses, applied to the vector rather
        than to the surface.
        """
        X, N, _ = self.frame()
        i = self._sample()
        P, D = X[i], N[i]
        L = float(length)
        t_exit = self._exit_parameter(P, D, self.shells.r_outer)
        over = t_exit < L

        tips = P + L * D
        clipped = np.where(over[:, None], P + t_exit[:, None] * D, tips)

        exit_pt = P[over] + t_exit[over][:, None] * D[over]
        ang = np.arctan2(exit_pt[:, 1], exit_pt[:, 0])
        base = self.shells.r_inner * np.stack([np.cos(ang), np.sin(ang)],
                                              axis=-1)
        stub = base + (L - t_exit[over])[:, None] * D[over]
        return P, clipped, base, stub

    # ── counting ────────────────────────────────────────────────────────────
    def crossings(self, length: float, with_reset: bool = False) -> int:
        if not with_reset:
            a, b = self.vectors(length)
            return segment_crossings(a, b)
        a, b, a2, b2 = self.vectors_with_reset(length)
        return segment_crossings(np.vstack([a, a2]), np.vstack([b, b2]))

    def n_wrapped(self, length: float) -> int:
        _, _, base, _ = self.vectors_with_reset(length)
        return int(len(base))


# ════════════════════════════════════════════════════════════════════════════
# MEASUREMENTS
# ════════════════════════════════════════════════════════════════════════════
def measure_normals_intersect_where_a_graph_cannot(
        t: float = 3.0, delta: float = 0.26, gain: float = 0.30,
        lengths: Sequence[float] = (0.10, 0.20, 0.26, 0.35, 0.50)
        ) -> Dict[str, object]:
    """The same wave, the same instant, drawn two ways.

    The graph of the tips is embedded — that is the whole of the earlier
    negative result.  The normals whose tips those are cross each other freely.
    Nothing about the field changed; only which object was drawn.
    """
    from geometrodynamics.viz.slice_folding import self_intersections

    nf = NormalField(delta=delta, gain=gain)
    nf.reset()
    nf.advance_to(t)
    graph = self_intersections(nf.slice.points(gain=gain))
    rows = [{"length": float(L), "normal_crossings": nf.crossings(L)}
            for L in lengths]
    return {
        "time": t,
        "graph_self_intersections": graph,
        "rows": rows,
        "most_normal_crossings": max(r["normal_crossings"] for r in rows),
        "envelope_distance": nf.envelope_distance(),
        "the_graph_never_crosses": bool(graph == 0),
        "the_normals_do": bool(max(r["normal_crossings"] for r in rows) > 0),
    }


def measure_the_threshold_is_the_radius_of_curvature(
        times: Sequence[float] = (1.2, 2.0, 2.6, 3.0), delta: float = 0.26,
        gain: float = 0.30) -> Dict[str, object]:
    """Neighbouring normals meet at the centre of curvature — so ``L* ≈ ρ_min``.

    The drawn bundle is sampled at a finite stride, so the first *drawn*
    crossing lags the continuous envelope: infinitesimally-neighbouring normals
    meet at ``ρ``, sampled ones need a little more.  Both are reported, and the
    point is that they move together as the wave sharpens.
    """
    rows = []
    for t in times:
        nf = NormalField(delta=delta, gain=gain)
        nf.reset()
        nf.advance_to(float(t))
        rho = nf.envelope_distance()
        lo, hi = 0.0, 1.2
        for _ in range(34):
            mid = 0.5 * (lo + hi)
            if nf.crossings(mid) > 0:
                hi = mid
            else:
                lo = mid
        rows.append({"t": float(t), "rho_min": rho,
                     "first_drawn_crossing": 0.5 * (lo + hi)})
    rhos = [r["rho_min"] for r in rows]
    return {
        "rows": rows,
        "rho_at_the_start": rhos[0],
        "rho_at_the_focus": rhos[-1],
        "sharpening_factor": rhos[0] / max(rhos[-1], 1e-30),
        "the_focus_sharpens_the_surface": bool(rhos[-1] < rhos[0]),
        "the_envelope_bounds_the_drawn_crossing": bool(
            all(r["first_drawn_crossing"] >= r["rho_min"] - 1e-9
                for r in rows)),
    }


def measure_the_reset_adds_intersections(
        t: float = 3.0, delta: float = 0.26, gain: float = 0.30,
        lengths: Sequence[float] = (0.20, 0.26, 0.35, 0.50)
        ) -> Dict[str, object]:
    """A second, independent mechanism: the stub that comes back inside.

    A vector leaving through ``R_outer`` re-enters at ``R_inner`` and shoots
    outward from deep inside the annulus, so it crosses vectors it could never
    otherwise have reached.  Counted against the same bundle without the reset,
    so the two mechanisms are separated rather than conflated.
    """
    nf = NormalField(delta=delta, gain=gain)
    nf.reset()
    nf.advance_to(t)
    rows = []
    for L in lengths:
        plain = nf.crossings(L)
        both = nf.crossings(L, with_reset=True)
        rows.append({"length": float(L), "wrapped": nf.n_wrapped(L),
                     "without_reset": plain, "with_reset": both,
                     "added": both - plain})
    return {
        "time": t,
        "rows": rows,
        "most_added": max(r["added"] for r in rows),
        "biggest_ratio": max(r["with_reset"] / max(r["without_reset"], 1)
                             for r in rows),
        "the_reset_adds_crossings": bool(max(r["added"] for r in rows) > 0),
        "it_grows_with_the_stub": bool(
            rows[-1]["added"] > rows[0]["added"]),
    }


def measure_the_gap_matters_again(
        deltas: Sequence[float] = (0.40, 0.26, 0.16, 0.09),
        t: float = 3.0, gain: float = 0.30) -> Dict[str, object]:
    """What the height representation had severed, the normals restore.

    ``slice_folding`` showed the fold threshold was *independent* of the gap, so
    shrinking the vacuole could never buy an intersection.  Here the vector
    length is what spans the gap, so the two are the same knob: each bundle is
    drawn at ``L = δ``, the length that just reaches across, and the crossing
    count is read off.
    """
    rows = []
    for delta in deltas:
        nf = NormalField(delta=float(delta), gain=gain)
        nf.reset()
        nf.advance_to(t)
        L = float(delta)
        rows.append({"delta": float(delta), "length": L,
                     "rho_min": nf.envelope_distance(),
                     "wrapped": nf.n_wrapped(L),
                     "crossings": nf.crossings(L),
                     "with_reset": nf.crossings(L, with_reset=True)})
    return {
        "time": t,
        "rows": rows,
        "the_gap_changes_the_count": bool(
            len({r["with_reset"] for r in rows}) > 1),
        "every_gap_crosses": bool(all(r["with_reset"] > 0 for r in rows)),
    }


def segment_crossing_points(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Where the segments actually meet — the caustic, as points to draw.

    Same predicate as ``segment_crossings``, but returning the intersections
    rather than counting them, so a picture can show the envelope instead of
    asserting it.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    out: List[np.ndarray] = []
    for i in range(len(a) - 1):
        u = b[i] - a[i]
        v = b[i + 1:] - a[i + 1:]
        d1 = _cross(u, a[i + 1:] - a[i]) * _cross(u, b[i + 1:] - a[i])
        d2 = _cross(v, a[i] - a[i + 1:]) * _cross(v, b[i] - a[i + 1:])
        hit = np.nonzero((d1 < 0.0) & (d2 < 0.0))[0]
        if not len(hit):
            continue
        j = i + 1 + hit
        denom = _cross(u[None, :], b[j] - a[j])
        ok = np.abs(denom) > 1e-14
        if not np.any(ok):
            continue
        t = _cross(a[j][ok] - a[i], b[j][ok] - a[j]) / denom[ok]
        out.append(a[i] + t[:, None] * u[None, :])
    return np.vstack(out) if out else np.zeros((0, 2))


def measure_the_curvature_converges(
        t: float = 3.06, delta: float = 0.26, gain: float = 0.30,
        grids: Sequence[Tuple[int, int]] = ((1001, 900), (2001, 1200),
                                            (4001, 2400), (8001, 4800))
        ) -> Dict[str, object]:
    """``ρ_min`` is a number, not a resolution — checked before it is quoted.

    The curvature needs a second derivative of the solved field, so the angular
    sampling and the radial solve have to be refined **together**.  Refining
    ``σ`` against a coarse solve reports a ``ρ_min`` that is mostly
    interpolation noise, and reports it confidently.
    """
    rows = []
    for n_sigma, n_radial in grids:
        nf = NormalField(delta=delta, gain=gain, n_sigma=n_sigma,
                         n_radial=n_radial, stride=8)
        nf.reset()
        nf.advance_to(t)
        rows.append({"n_sigma": n_sigma, "n_radial": n_radial,
                     "rho_min": nf.envelope_distance()})
    rhos = [r["rho_min"] for r in rows]
    drift = abs(rhos[-1] - rhos[-2]) / max(rhos[-1], 1e-30)

    # ...and what a mismatched grid claims instead
    bad = NormalField(delta=delta, gain=gain, n_sigma=8001, n_radial=900,
                      stride=8)
    bad.reset()
    bad.advance_to(t)
    return {
        "time": t,
        "rows": rows,
        "converged_rho_min": rhos[-1],
        "last_step_drift": drift,
        "mismatched_grid_rho_min": bad.envelope_distance(),
        "mismatched_error": abs(bad.envelope_distance() - rhos[-1])
        / max(rhos[-1], 1e-30),
        "it_converges": bool(drift < 0.06),
        "a_mismatched_grid_gets_it_wrong": bool(
            abs(bad.envelope_distance() - rhos[-1]) > 0.1 * rhos[-1]),
    }
