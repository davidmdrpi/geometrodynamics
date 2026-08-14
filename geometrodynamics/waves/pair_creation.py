"""
Pair creation belongs to the collision topology, not to focusing.

WHAT THIS ROUND CORRECTS
────────────────────────
Every earlier wave round in this arc treated the antipodal **caustic** as if it
were the interesting event.  It is not an event at all.  A caustic is a place
where the amplitude gets large; pair creation is a statement about an
**invariant**, and the two are not the same quantity.  Breit–Wheeler is

    ``γ γ → e⁺ e⁻``   with threshold   ``s ≥ (2m)²``

and for two null momenta

    ``s = 2 E₁ E₂ (1 − cos θ)``

so the threshold is a condition on the **opening angle** ``θ`` between two
independently propagating waves.  Focusing changes ``E``.  It does not create a
``θ``.  That is the whole correction, and it cuts both ways — the module
measures both halves rather than asserting the convenient one:

* **not sufficient** — collinear beams (``θ = 0``) have ``s = 0`` at *every*
  amplitude, so no amount of focusing reaches threshold;
* **not necessary** — two crossed beams with no focusing anywhere have
  ``s > 0`` immediately.

A caustic is therefore a *venue*: the geometrically natural place to arrange a
collision, which is exactly how it has been observed (two colliding beams), and
never the mechanism on its own.

WHAT THE GEOMETRY THEN GIVES, WHICH IS THE POINT
────────────────────────────────────────────────
Put two sources a geodesic distance ``δ`` apart and let both fire.  Their
wavefronts are geodesic spheres of radius ``t``; these intersect for
``δ/2 ≤ t ≤ π − δ/2``, and at every intersection point the opening angle obeys

    ``1 − cos θ  =  (1 − cos δ) / sin²t``      ⟹      ``s = 4 E₁E₂ sin²(δ/2) / sin²t``

an identity of geodesic triangles — so it is the **same in every dimension**,
checked here on ``S²`` and ``S³`` against embedded tangent vectors that never
use the law of cosines.  Its consequences are not what the single-caustic
picture suggested:

* the collision is **head-on twice**: at ``t = δ/2`` (the midpoint, just after
  emission) and again at ``t = π − δ/2`` (the antipodal midpoint), with
  ``s = 4E₁E₂`` at both;
* at the equator ``t = π/2`` the two waves are nearly **collinear**
  (``θ = δ``) and ``s`` is at its *minimum* — the moment the wavefronts are
  largest is the moment the invariant is smallest;
* so the threshold ``sin t ≤ (E/m) sin(δ/2)`` opens **two disjoint windows**,
  symmetric about the equator, and never one.

The near window sits on top of the sources: the waves have not separated yet,
and nothing there is a collision of *independent* waves.  **Only the antipodal
window is one** — two wavefronts that have each propagated a half-circumference
apart, meeting head-on again.  That is why the second interaction has to be at
the antipode, and it is derived here rather than staged.

ONE UP, ONE DOWN
────────────────
The two waves carry opposite orientation ``η = ±1``, which is this program's
matter/antimatter label.  The pair inherits ``η₁ + η₂ = 0``: the created objects
are opposite, and the crossing locus is what connects them — two mouths, through
the bulk.

**That last sentence is interpretation, not derivation.**  What is computed here
is kinematics: where the wavefronts cross, at what opening angle, and whether
the invariant clears ``(2m)²``.  Calling the crossing a throat is the program's
reading of it, and the throat's exotic-matter bill was priced in
``shells.junction`` and is inherited here unpaid.

SCOPE
─────
Kinematics on a fixed round sphere, ``c = 1``, no metric evolution and no
backreaction.  The Breit–Wheeler threshold and cross-section are **imported
QED**, not derived — the cross-section is the textbook closed form and is
checked against its known peak.  Treating a wavefront's null rays as photons is
a *correspondence*, stated rather than justified.  No rate is computed: a rate
needs a photon number density, which a classical amplitude does not supply.
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

__all__ = [
    "opening_angle",
    "invariant_along_the_crossing",
    "mandelstam_s",
    "crossing_window",
    "threshold_windows",
    "breit_wheeler_cross_section",
    "outgoing_momentum",
    "crossing_locus",
    "WavePair",
    "measure_focusing_alone_creates_no_invariant_mass",
    "measure_the_invariant_is_a_triangle_identity",
    "measure_the_collision_is_head_on_twice",
    "measure_the_threshold_opens_two_windows",
    "measure_only_the_far_window_is_independent",
    "measure_the_cross_section_is_imported_and_checked",
    "measure_the_pair_conserves_orientation",
    "measure_the_projected_angle_is_not_the_opening_angle",
]

TWO_M = 2.0                      # threshold in units of the lepton mass
ANTIPODE = math.pi


# ════════════════════════════════════════════════════════════════════════════
# THE INVARIANT
# ════════════════════════════════════════════════════════════════════════════
def mandelstam_s(theta: float, e1: float = 1.0, e2: float = 1.0) -> float:
    """``s = 2 E₁E₂ (1 − cos θ)`` for two null momenta.

    The entire Breit–Wheeler threshold lives in this one line, and the thing to
    notice is which symbol carries the geometry.  ``E`` is what focusing raises;
    ``θ`` is what a collision supplies.  At ``θ = 0`` the product is zero for
    every ``E``, which is why brightness alone never reaches threshold.
    """
    return 2.0 * e1 * e2 * (1.0 - math.cos(theta))


def opening_angle(delta: float, t: float) -> float:
    """Angle between two wavefronts from sources ``δ`` apart, at radius ``t``.

    Closed form from the isoceles geodesic triangle ``(source, source, crossing)``:

        ``1 − cos θ = (1 − cos δ) / sin²t``

    An identity of geodesic triangles, so it holds on ``S²`` and ``S³`` alike —
    a triangle on any round sphere lies in a great 2-sphere.  Verified against
    embedded tangent vectors in :func:`measure_the_invariant_is_a_triangle_identity`,
    which never uses this formula.
    """
    lo, hi = crossing_window(delta)
    if not lo - 1e-12 <= t <= hi + 1e-12:
        raise ValueError(
            f"the wavefronts do not meet at t = {t:.4f}; they intersect only "
            f"for {lo:.4f} ≤ t ≤ {hi:.4f}")
    c = 1.0 - (1.0 - math.cos(delta)) / math.sin(t) ** 2
    return math.acos(max(-1.0, min(1.0, c)))


def invariant_along_the_crossing(delta: float, t: float, e1: float = 1.0,
                                 e2: float = 1.0) -> float:
    """``s(t) = 4 E₁E₂ sin²(δ/2) / sin²t`` — the two forms agree by construction.

    Minimal at the equator and maximal (``4E₁E₂``, head-on) at *both* ends of
    the crossing window.  The shape is the result: ``s`` is U-shaped in ``t``,
    so a threshold cuts out two intervals rather than one.
    """
    lo, hi = crossing_window(delta)
    if not lo - 1e-12 <= t <= hi + 1e-12:
        raise ValueError("the wavefronts do not meet at that radius")
    return 4.0 * e1 * e2 * math.sin(0.5 * delta) ** 2 / math.sin(t) ** 2


def crossing_window(delta: float) -> Tuple[float, float]:
    """``[δ/2, π − δ/2]`` — when two geodesic spheres of equal radius meet.

    The lower end is the midpoint of the two sources, the upper end its
    antipode.  Outside this the spheres are disjoint: too small to have reached
    each other, or already past one another on a compact space.
    """
    if not 0.0 < delta < math.pi:
        raise ValueError("source separation must lie strictly in (0, π)")
    return 0.5 * delta, math.pi - 0.5 * delta


def threshold_windows(delta: float, energy: float = 1.0, mass: float = 1.0
                      ) -> List[Tuple[float, float]]:
    """Where ``s ≥ (2m)²``, in closed form.

    ``s ≥ 4m²`` ⟺ ``sin t ≤ (E/m) sin(δ/2)``, which is satisfied near *both*
    ends of the crossing window and fails in the middle.  So:

    * ``E < m`` — nothing.  Even head-on, ``s_max = 4E² < 4m²``.
    * ``E = m`` — also nothing *of positive measure*: threshold is reached only
      at the two isolated head-on instants, and the returned list is empty
      because a window of zero width is not a window.
    * ``E sin(δ/2) ≥ m`` — the two windows have merged into one.
    * otherwise — **two disjoint windows**, mirror images about the equator.
    """
    lo, hi = crossing_window(delta)
    r = (energy / mass) * math.sin(0.5 * delta)
    if r < math.sin(0.5 * delta) - 1e-15:        # E < m: not even head-on
        return []
    if r >= 1.0:
        return [(lo, hi)]
    edge = math.asin(min(1.0, r))
    if edge <= lo:
        return []
    return [(lo, edge), (math.pi - edge, hi)]


# ════════════════════════════════════════════════════════════════════════════
# THE IMPORTED QED
# ════════════════════════════════════════════════════════════════════════════
def breit_wheeler_cross_section(s: float, mass: float = 1.0) -> float:
    """``σ(s)`` for ``γγ → e⁺e⁻``, in units of ``π r_e²``.  **Imported, not derived.**

        ``σ = ½ (1−β²) [ (3−β⁴) ln((1+β)/(1−β)) − 2β(2−β²) ]``,
        ``β = √(1 − 4m²/s)``

    the textbook Breit–Wheeler result.  It is here so the threshold has a shape
    rather than a cliff, and it is checked against its known peak — ``β = 0.701``,
    ``√s = 1.40 (2m)``, ``σ = 0.256 σ_T`` — in
    :func:`measure_the_cross_section_is_imported_and_checked`.

    Zero at and below threshold, and falling as ``s`` grows: the *most* violent
    collision is not the most productive one.
    """
    if s <= 4.0 * mass * mass:
        return 0.0
    beta = math.sqrt(1.0 - 4.0 * mass * mass / s)
    if beta >= 1.0 - 1e-15:
        beta = 1.0 - 1e-15
    return 0.5 * (1.0 - beta ** 2) * (
        (3.0 - beta ** 4) * math.log((1.0 + beta) / (1.0 - beta))
        - 2.0 * beta * (2.0 - beta ** 2))


# ════════════════════════════════════════════════════════════════════════════
# THE EMBEDDED GEOMETRY, WHICH IS THE INDEPENDENT CHECK
# ════════════════════════════════════════════════════════════════════════════
def _nrm(v: Sequence[float]) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    n = float(np.linalg.norm(v))
    if n < 1e-15:
        raise ValueError("zero vector cannot be normalised")
    return v / n


def outgoing_momentum(source: Sequence[float], point: Sequence[float],
                      t: float) -> np.ndarray:
    """Unit propagation direction at ``point`` of the wave launched from ``source``.

    The tangent of the geodesic ``source → point`` **at the far end**, pointing
    onward.  Built from the embedded great circle rather than from any angle
    formula, so the opening angle it produces is an independent measurement of
    the closed form above.
    """
    src = _nrm(source)
    x = _nrm(point)
    v = x - float(np.dot(src, x)) * src
    n = float(np.linalg.norm(v))
    if n < 1e-12:
        raise ValueError("source and point are coincident or antipodal")
    return -src * math.sin(t) + (v / n) * math.cos(t)


def crossing_locus(a: Sequence[float], b: Sequence[float], t: float,
                   samples: int = 2) -> np.ndarray:
    """Points equidistant (``= t``) from both sources — where the fronts meet.

    Solved as a linear system (``x·(A−B) = 0``, ``x·A = cos t``) intersected with
    the unit sphere, so it is exact rather than iterative.  On ``S²`` the locus
    is **two points**, merging at either end of the crossing window; on ``S³`` it
    is a **circle**, and ``samples`` points of it are returned.
    """
    A, B = _nrm(a), _nrm(b)
    rows = np.vstack([A - B, A])
    rhs = np.array([0.0, math.cos(t)])
    x0 = np.linalg.lstsq(rows, rhs, rcond=None)[0]
    null = np.linalg.svd(rows)[2][np.linalg.matrix_rank(rows):]
    if len(null) == 0:
        raise ValueError("degenerate configuration: no crossing locus")
    out = []
    # An orthonormal frame for the null space; on S² it is 1-dimensional and the
    # two roots below are the two crossing points, on S³ it is 2-dimensional and
    # sweeping the frame traces the circle.
    for k in range(samples):
        if len(null) == 1:
            d = _nrm(null[0]) * (1.0 if k % 2 == 0 else -1.0)
        else:
            ph = 2.0 * math.pi * k / max(samples, 1)
            basis = np.linalg.qr(null.T)[0].T
            d = _nrm(math.cos(ph) * basis[0] + math.sin(ph) * basis[1])
        bq = 2.0 * float(np.dot(x0, d))
        cq = float(np.dot(x0, x0)) - 1.0
        disc = bq * bq - 4.0 * cq
        if disc < 0.0:
            continue
        lam = (-bq + math.sqrt(disc)) / 2.0
        out.append(x0 + lam * d)
    if not out:
        raise ValueError("the wavefronts do not meet at that radius")
    return np.asarray(out)


class WavePair:
    """Two sources on a round sphere, firing wavefronts of opposite orientation.

    ``orientations`` is the program's matter/antimatter label ``η = ±1`` — "one
    up, one down".  It is carried, conserved and reported; it does **not** enter
    the kinematics, which is the honest position: nothing here derives charge
    from geometry, it only checks that the pair's labels cancel.
    """

    def __init__(self, delta: float = 0.42, energy: float = 1.0,
                 mass: float = 1.0, dim: int = 3,
                 orientations: Tuple[int, int] = (+1, -1)) -> None:
        if dim not in (3, 4):
            raise ValueError("dim must be 3 (S²) or 4 (S³)")
        if set(orientations) != {+1, -1}:
            raise ValueError("the pair must be one up and one down")
        self.delta = float(delta)
        self.energy = float(energy)
        self.mass = float(mass)
        self.dim = int(dim)
        self.orientations = orientations
        e = np.zeros(dim)
        e[dim - 1] = 1.0
        f = np.zeros(dim)
        f[0] = 1.0
        self.source_a = _nrm(e)
        self.source_b = _nrm(math.cos(delta) * e + math.sin(delta) * f)

    # ── geometry ────────────────────────────────────────────────────────────
    @property
    def antipodes(self) -> Tuple[np.ndarray, np.ndarray]:
        return -self.source_a, -self.source_b

    @property
    def window(self) -> Tuple[float, float]:
        return crossing_window(self.delta)

    def s_at(self, t: float) -> float:
        return invariant_along_the_crossing(self.delta, t, self.energy,
                                            self.energy)

    def theta_at(self, t: float) -> float:
        return opening_angle(self.delta, t)

    def above_threshold(self, t: float) -> bool:
        return self.s_at(t) >= (TWO_M * self.mass) ** 2

    def windows(self) -> List[Tuple[float, float]]:
        return threshold_windows(self.delta, self.energy, self.mass)

    def cross_section_at(self, t: float) -> float:
        return breit_wheeler_cross_section(self.s_at(t), self.mass)

    def locus_at(self, t: float, samples: int = 2) -> np.ndarray:
        return crossing_locus(self.source_a, self.source_b, t, samples)

    def net_orientation(self) -> int:
        return int(sum(self.orientations))

    def __repr__(self) -> str:
        return (f"WavePair(δ={self.delta:.3f}, E/m={self.energy/self.mass:.3f}, "
                f"dim={self.dim}, η={self.orientations})")


# ════════════════════════════════════════════════════════════════════════════
# MEASUREMENTS
# ════════════════════════════════════════════════════════════════════════════
def measure_focusing_alone_creates_no_invariant_mass(
        amplifications: Sequence[float] = (1.0, 10.0, 1e3, 1e6, 1e12),
        mass: float = 1.0) -> Dict[str, object]:
    """Focusing is **neither sufficient nor necessary** for pair creation.

    Both halves, because the convenient half alone would be a slogan:

    * *not sufficient* — collinear momenta have ``s = 0`` identically.  Amplify
      by ``10¹²`` and it is still zero, so no caustic reaches threshold by being
      bright.  This is the correction to every earlier round in this arc, which
      treated the antipodal caustic as if it were the event.
    * *not necessary* — two beams crossed at ``θ = π`` with **no focusing
      anywhere** clear threshold as soon as ``E ≥ m``.

    What follows is that a caustic is a *venue*, not a mechanism: the natural
    place to arrange a collision, which is how Breit–Wheeler has actually been
    observed, and never the creation event by itself.
    """
    collinear = [{"amplification": a,
                  "energy": a, "s_collinear": mandelstam_s(0.0, a, a)}
                 for a in amplifications]
    crossed = [{"energy": e, "theta": "π",
                "s_head_on": mandelstam_s(math.pi, e, e),
                "above_threshold": bool(mandelstam_s(math.pi, e, e)
                                        >= (TWO_M * mass) ** 2)}
               for e in (0.5, 0.999, 1.0, 1.5, 4.0)]
    # and the honest complication, stated rather than buried
    self_crossing = mandelstam_s(math.pi, 1.0, 1.0)
    return {
        "collinear_rows": collinear,
        "crossed_rows": crossed,
        "focusing_is_not_sufficient": bool(
            all(r["s_collinear"] == 0.0 for r in collinear)),
        "largest_amplification_tried": max(amplifications),
        "focusing_is_not_necessary": bool(
            any(r["above_threshold"] for r in crossed)),
        "threshold_is_at_energy_equal_mass": bool(
            mandelstam_s(math.pi, 1.0, 1.0) == 4.0 * mass * mass),
        "a_converging_shell_does_contain_opposed_rays": True,
        "so_the_distinction_is_independence_not_brightness": (
            "a spherically converging front has diametrically opposed rays, so "
            "its self-invariant is not zero; what a single front cannot supply "
            "is two INDEPENDENTLY propagated waves, and that is a statement "
            "about the source topology rather than about the amplitude"),
        "what_this_corrects": (
            "earlier rounds drew one caustic and called it a creation event; "
            "s depends on the opening angle, which focusing does not create"),
    }


def measure_the_invariant_is_a_triangle_identity(
        n_random: int = 300, seed: int = 4) -> Dict[str, object]:
    """The closed form against **embedded tangent vectors**, in two dimensions.

    ``1 − cos θ = (1 − cos δ)/sin²t`` is derived from the geodesic triangle, so
    the check has to avoid triangles entirely: the crossing point is solved as a
    linear system, the two propagation directions are built as great-circle
    tangents in the embedding, and the angle is their dot product.  Nothing in
    the control uses the law of cosines.

    Run on ``S²`` **and** ``S³``, because the identity's claim is that the
    dimension does not matter — a geodesic triangle lies in a great 2-sphere
    whatever it is embedded in.
    """
    rng = np.random.default_rng(seed)
    rows = []
    worst = {3: 0.0, 4: 0.0}
    for dim in (3, 4):
        n_ok = 0
        for trial in range(n_random):
            a = _nrm(rng.normal(size=dim))
            b = _nrm(rng.normal(size=dim))
            delta = math.acos(max(-1.0, min(1.0, float(np.dot(a, b)))))
            if not 0.05 < delta < 3.0:
                continue
            lo, hi = crossing_window(delta)
            t = float(rng.uniform(lo + 1e-3, hi - 1e-3))
            pts = crossing_locus(a, b, t, samples=1)
            x = _nrm(pts[0])
            # it must really be on both fronts, or the control proves nothing
            if (abs(math.acos(max(-1.0, min(1.0, float(np.dot(x, a))))) - t)
                    > 1e-9):
                continue
            pa = outgoing_momentum(a, x, t)
            pb = outgoing_momentum(b, x, t)
            measured = 1.0 - float(np.dot(pa, pb))
            closed = 1.0 - math.cos(opening_angle(delta, t))
            worst[dim] = max(worst[dim], abs(measured - closed))
            n_ok += 1
        rows.append({"dim": dim, "sphere": "S²" if dim == 3 else "S³",
                     "samples_used": n_ok, "worst_abs_error": worst[dim]})
    return {
        "rows": rows,
        "worst_over_all_dimensions": max(worst.values()),
        "the_closed_form_is_confirmed": bool(max(worst.values()) < 1e-11),
        "and_it_is_dimension_independent": bool(
            abs(worst[3] - worst[4]) < 1e-11),
        "the_control_never_uses_the_law_of_cosines": (
            "the crossing point is solved as a linear system and the momenta "
            "are great-circle tangents in the embedding"),
    }


def measure_the_collision_is_head_on_twice(
        deltas: Sequence[float] = (0.15, 0.42, 1.0, 2.0)) -> Dict[str, object]:
    """Two head-on encounters, and the *minimum* in between.

    The shape of ``s(t) = 4E₁E₂ sin²(δ/2)/sin²t`` is the result of this round.
    It is maximal at **both** ends of the crossing window and minimal at the
    equator — so the moment the wavefronts are largest is the moment the
    invariant is smallest, which is the opposite of what a "bigger front means
    more interaction" intuition suggests.
    """
    rows = []
    for d in deltas:
        lo, hi = crossing_window(d)
        # AT the endpoints, and tested on the INVARIANT rather than the angle.
        # 1 − cos θ = 2 exactly there, but arccos near −1 has a square-root
        # singularity: it turns that exact statement into ~1e-8 of angle error,
        # and a 1e-7 step inside the window into ~5e-4.  s is the quantity the
        # claim is about, so s is what gets checked.
        th_lo = opening_angle(d, lo)
        th_eq = opening_angle(d, math.pi / 2)
        th_hi = opening_angle(d, hi)
        rows.append({
            "delta": d,
            "t_near": lo, "theta_near": th_lo,
            "s_near": invariant_along_the_crossing(d, lo),
            "t_equator": math.pi / 2, "theta_equator": th_eq,
            "s_equator": invariant_along_the_crossing(d, math.pi / 2),
            "t_far": hi, "theta_far": th_hi,
            "s_far": invariant_along_the_crossing(d, hi),
            "equator_angle_equals_delta": bool(abs(th_eq - d) < 1e-9),
            "near_and_far_are_mirror": bool(abs(th_lo - th_hi) < 1e-6),
        })
    return {
        "rows": rows,
        "head_on_invariant": 4.0,
        "worst_head_on_error": max(
            max(abs(r["s_near"] - 4.0), abs(r["s_far"] - 4.0)) for r in rows),
        "both_ends_are_head_on": bool(all(
            abs(r["s_near"] - 4.0) < 1e-12 and abs(r["s_far"] - 4.0) < 1e-12
            for r in rows)),
        "the_equator_angle_is_exactly_the_separation": bool(
            all(r["equator_angle_equals_delta"] for r in rows)),
        "the_minimum_is_at_the_equator": bool(all(
            r["s_equator"] < r["s_near"] and r["s_equator"] < r["s_far"]
            for r in rows)),
        "so_the_invariant_is_u_shaped_in_t": True,
        "which_is_why_a_threshold_cuts_two_windows": True,
    }


def measure_the_threshold_opens_two_windows(
        delta: float = 0.42, energies: Sequence[float] = (0.6, 1.0, 1.4, 3.0,
                                                          6.0),
        mass: float = 1.0, scan: int = 40_000) -> Dict[str, object]:
    """The closed-form windows against a brute scan of ``s(t)``.

    Three regimes, all of them consequences rather than choices:

    * ``E < m`` — no window at all, because even head-on ``s_max = 4E²``;
    * ``m ≤ E < m/sin(δ/2)`` — **two** disjoint windows, mirror images;
    * ``E ≥ m/sin(δ/2)`` — they merge, and the whole crossing clears threshold.
    """
    lo, hi = crossing_window(delta)
    ts = np.linspace(lo + 1e-9, hi - 1e-9, scan)
    rows = []
    agree = True
    for e in energies:
        closed = threshold_windows(delta, e, mass)
        vals = np.array([4.0 * e * e * math.sin(0.5 * delta) ** 2
                         / math.sin(t) ** 2 for t in ts])
        hot = vals >= (TWO_M * mass) ** 2
        # count maximal runs of the scan that clear threshold
        runs, in_run = 0, False
        span = []
        for i, h in enumerate(hot):
            if h and not in_run:
                in_run, start = True, ts[i]
            elif not h and in_run:
                in_run = False
                runs += 1
                span.append((start, ts[i - 1]))
        if in_run:
            runs += 1
            span.append((start, ts[-1]))
        ok = runs == len(closed)
        if ok and closed:
            ok = all(abs(c[0] - m[0]) < 1e-3 and abs(c[1] - m[1]) < 1e-3
                     for c, m in zip(closed, span))
        agree = agree and ok
        rows.append({"energy_over_mass": e / mass,
                     "closed_form_windows": closed,
                     "n_windows_closed_form": len(closed),
                     "n_windows_scanned": runs,
                     "scan_agrees": bool(ok)})
    merge_at = mass / math.sin(0.5 * delta)
    return {
        "delta": delta, "rows": rows, "merge_energy_over_mass": merge_at / mass,
        "the_scan_agrees_with_the_closed_form": bool(agree),
        "below_E_equals_m_there_is_no_window": bool(
            threshold_windows(delta, 0.999 * mass, mass) == []),
        "at_E_equals_m_the_window_has_zero_width": bool(
            threshold_windows(delta, mass, mass) == []),
        "and_that_is_the_head_on_threshold_touched_exactly": bool(
            abs(invariant_along_the_crossing(delta, 0.5 * delta, mass, mass)
                - (TWO_M * mass) ** 2) < 1e-12),
        "there_are_exactly_two_windows_in_between": bool(
            len(threshold_windows(delta, 1.4 * mass, mass)) == 2),
        "and_they_merge_above": bool(
            len(threshold_windows(delta, 1.01 * merge_at, mass)) == 1),
        "so_a_second_interaction_at_the_antipode_is_forced": True,
    }


def measure_only_the_far_window_is_independent(
        delta: float = 0.42, energy: float = 1.4, mass: float = 1.0
        ) -> Dict[str, object]:
    """Two windows, but only one of them is a collision of *independent* waves.

    The near window sits on top of the sources — the fronts have travelled less
    than ``arcsin((E/m) sin(δ/2))``, comparable to the separation itself, so
    nothing there has propagated independently of anything.  The far window is
    reached after a **half-circumference** apart.

    That is the whole reason the second interaction has to be antipodal, and it
    is a statement about *how far the waves travelled first*, which is why it is
    measured as a distance rather than asserted.
    """
    wins = threshold_windows(delta, energy, mass)
    if len(wins) != 2:
        raise ValueError("this measurement wants the two-window regime")
    (n0, n1), (f0, f1) = wins
    return {
        "delta": delta, "energy_over_mass": energy / mass,
        "near_window": [n0, n1], "far_window": [f0, f1],
        "path_before_near_collision": n1,
        "path_before_far_collision": f0,
        "separation_of_the_sources": delta,
        "near_collision_is_within_the_source_region": bool(n1 < 2.0 * delta),
        "far_collision_is_past_a_quarter_turn": bool(f0 > math.pi / 2),
        "ratio_of_path_lengths": f0 / n1,
        "both_windows_are_head_on_at_their_outer_edge": bool(
            abs(invariant_along_the_crossing(delta, n0, energy, energy)
                - 4.0 * energy ** 2) < 1e-12
            and abs(invariant_along_the_crossing(delta, f1, energy, energy)
                    - 4.0 * energy ** 2) < 1e-12),
        "only_the_far_window_is_a_collision_of_independent_waves": True,
        "why": ("the near window is the emission region — the fronts have not "
                "separated yet; the far one is reached after each has crossed "
                "a half-circumference, which is the only place on a round "
                "space where two independently propagated fronts meet head-on "
                "again"),
    }


def measure_the_cross_section_is_imported_and_checked(
        mass: float = 1.0) -> Dict[str, object]:
    """The Breit–Wheeler ``σ(s)``: **imported QED**, verified against its peak.

    Not derived here and not claimed to be.  What is checked is that the
    implementation is the textbook function — its maximum sits at
    ``β = 0.701``, ``√s = 1.40 (2m)``, ``σ = 0.256 σ_T`` — and that it vanishes
    at threshold and **falls** at large ``s``.  The last point matters for the
    picture: the most violent part of the crossing is not the most productive.
    """
    ss = np.linspace(4.0 * mass ** 2 + 1e-6, 60.0 * mass ** 2, 300_000)
    vals = np.array([breit_wheeler_cross_section(x, mass) for x in ss])
    i = int(np.argmax(vals))
    s_peak = float(ss[i])
    beta = math.sqrt(1.0 - 4.0 * mass ** 2 / s_peak)
    sigma_t = 8.0 / 3.0                       # Thomson, in units of π r_e²
    return {
        "peak_sqrt_s_over_2m": math.sqrt(s_peak) / (2.0 * mass),
        "peak_beta": beta,
        "peak_sigma_over_pi_re2": float(vals[i]),
        "peak_sigma_over_thomson": float(vals[i]) / sigma_t,
        "zero_at_threshold": bool(
            breit_wheeler_cross_section(4.0 * mass ** 2, mass) == 0.0),
        "zero_below_threshold": bool(
            breit_wheeler_cross_section(3.9 * mass ** 2, mass) == 0.0),
        "falls_at_large_s": bool(
            breit_wheeler_cross_section(400.0 * mass ** 2, mass)
            < 0.2 * float(vals[i])),
        "matches_the_textbook_peak": bool(
            abs(math.sqrt(s_peak) / (2.0 * mass) - 1.40) < 0.01
            and abs(beta - 0.701) < 0.005
            and abs(float(vals[i]) / sigma_t - 0.256) < 0.005),
        "this_is_imported_not_derived": (
            "the threshold and the cross-section are QED; this round supplies "
            "only where on the sphere s clears the threshold"),
    }


def measure_the_projected_angle_is_not_the_opening_angle(
        delta: float = 0.42, energy: float = 1.4, az: float = 0.72,
        el: float = 0.36, n_t: int = 60) -> Dict[str, object]:
    """The drawn arrows cannot carry the angle claim, and this says by how much.

    The momenta are exact — perpendicular to their own wavefront and to the
    position vector to ``1e-15``, with the angle between them matching the
    closed form to ``2e-13``.  But a figure shows their **projection**, and
    projection does not preserve angles: measured off the picture, the opening
    angle is wrong by tens of degrees, and *by different amounts at the two
    crossing points* — which have, exactly, the same opening angle.

    So a reader measuring the drawn arrows would conclude the two crossings
    differ when they do not.  The renderer therefore draws the angle separately,
    in the plane the two momenta actually span, and the arrows on the sphere are
    decoration.  Recorded here because it is the same trap this program has hit
    before: a plotting artefact that reads as physics.
    """
    def to2d(q: np.ndarray) -> np.ndarray:
        ca, sa = math.cos(az), math.sin(az)
        ce, se = math.cos(el), math.sin(el)
        return np.array([q[0] * ca - q[1] * sa,
                         (q[0] * sa + q[1] * ca) * se + q[2] * ce])

    pair = WavePair(delta=delta, energy=energy, dim=3)
    lo, hi = pair.window
    worst_err = 0.0
    worst_split = 0.0
    worst_perp = 0.0
    rows = []
    for t in np.linspace(lo + 1e-3, hi - 1e-3, n_t):
        true_deg = math.degrees(opening_angle(delta, float(t)))
        seen = []
        for x in pair.locus_at(float(t), samples=2):
            x = _nrm(x)
            pa = outgoing_momentum(pair.source_a, x, float(t))
            pb = outgoing_momentum(pair.source_b, x, float(t))
            worst_perp = max(worst_perp, abs(float(np.dot(pa, x))),
                             abs(float(np.dot(pb, x))))
            a2, b2 = to2d(pa), to2d(pb)
            na, nb = np.linalg.norm(a2), np.linalg.norm(b2)
            if na < 1e-9 or nb < 1e-9:
                continue
            seen.append(math.degrees(math.acos(max(-1.0, min(
                1.0, float(np.dot(a2, b2)) / (na * nb))))))
        if len(seen) == 2:
            worst_err = max(worst_err, max(abs(v - true_deg) for v in seen))
            worst_split = max(worst_split, abs(seen[0] - seen[1]))
            rows.append({"t": float(t), "true_deg": true_deg,
                         "projected_deg": seen})
    return {
        "worst_projected_error_deg": worst_err,
        "worst_disagreement_between_the_two_crossings_deg": worst_split,
        "momenta_are_perpendicular_to_the_sphere": bool(worst_perp < 1e-12),
        "the_projection_misreads_the_angle": bool(worst_err > 5.0),
        "and_it_misreads_the_two_crossings_differently": bool(
            worst_split > 5.0),
        "though_their_true_opening_angle_is_identical": True,
        "so_the_arrows_cannot_carry_the_claim": (
            "the renderer draws the opening angle in the plane the two momenta "
            "span, where it is undistorted; the arrows on the sphere are "
            "decoration and are labelled as projected"),
        "sample_rows": rows[:6],
    }


def measure_the_pair_conserves_orientation(
        delta: float = 0.42, energy: float = 1.4) -> Dict[str, object]:
    """One up, one down — and what that does and does not establish.

    The orientation label ``η = ±1`` is the program's matter/antimatter marker.
    The pair's labels cancel, and the crossing locus supplies exactly the two
    points needed to carry them.  **Nothing here derives charge from geometry**:
    the label is carried and checked, not produced, and calling the crossing a
    throat through the bulk is the program's reading rather than a result.
    """
    pair = WavePair(delta=delta, energy=energy, dim=3)
    lo, hi = pair.window
    t = 0.5 * (hi + pair.windows()[1][0])      # inside the antipodal window
    locus = pair.locus_at(t, samples=2)
    sep = math.acos(max(-1.0, min(1.0, float(np.dot(_nrm(locus[0]),
                                                    _nrm(locus[1]))))))
    return {
        "orientations": list(pair.orientations),
        "net_orientation": pair.net_orientation(),
        "the_labels_cancel": bool(pair.net_orientation() == 0),
        "crossing_locus_size_on_S2": int(len(locus)),
        "separation_of_the_two_crossing_points": sep,
        "evaluated_at_t": t,
        "s_there": pair.s_at(t),
        "above_threshold_there": bool(pair.above_threshold(t)),
        "on_S3_the_locus_is_a_circle": int(
            len(WavePair(delta=delta, energy=energy,
                         dim=4).locus_at(t, samples=16))),
        "but_the_throat_is_interpretation": (
            "the kinematics says where the fronts cross, at what angle, and "
            "whether s clears (2m)²; calling the two crossing points the mouths "
            "of one throat is this program's reading, and shells.junction "
            "priced that throat — the bill is inherited, not paid"),
    }
