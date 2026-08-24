"""What higher dimension actually does to the bulk/intersection picture.

The premise
───────────
Three-dimensional intuition models an extra dimension as *the same sphere with
another direction available*.  That is wrong in a way that matters here, and
this module measures the ways.

The one that matters most for `regularized_center` is §5 below: a packet can
keep **exactly the same angular footprint** while its physical overlap shrinks
as ``f₀ⁿ``.  PR #268's review had already forced the ``q``-scoping of the
overlap *size* law; this module follows it to its conclusion, which is that the
finite bearing is much less like two ribbons squeezing together and much more
like a vast angular configuration space packed into a tiny proper region.

A correction first, because it is the usual one
───────────────────────────────────────────────
The familiar "spheres peak at 5D" is about the **volume of the unit ball**,

    ``V_d = π^{d/2}/Γ(d/2+1)`` ,   which peaks at ``d = 5`` ,

not the surface measure of the unit sphere,

    ``A_{d−1} = 2π^{d/2}/Γ(d/2)`` ,   which peaks at ``d = 7`` — the sphere
    being ``S⁶ ⊂ ℝ⁷`` .

And neither peak is a fact about dimension alone.  Comparing measures of
different dimensionality implicitly picks a length unit, and the peak moves
with it: at ``R = 0.5`` the ball peaks at ``d = 1``, at ``R = 2`` at ``d = 24``.
`measure_the_peak_is_the_ball_not_the_sphere_and_needs_a_unit` carries all of
it.  The interesting geometry is elsewhere.

1 · almost all of a sphere is at the equator of any chosen point
────────────────────────────────────────────────────────────────
The shell at geodesic angle ``χ`` from a pole has measure ``∝ sin^{n−1}χ dχ``.
Writing ``χ = π/2 + δ`` gives ``sin^{n−1}χ ≈ e^{−(n−1)δ²/2}``, so the occupied
band has width

    ``Δχ ~ 1/√n`` .

Measured: ``std(χ)·√n → 1.000000``.  **In high dimension almost every point is
nearly 90° from any chosen point** — not "somewhere randomly around", but
concentrated in a shrinking band.

2 · so the antipode is an extremely non-generic relation
────────────────────────────────────────────────────────
For random unit vectors, ``x·y`` has mean ``0`` and width ``1/√n``, so
``α = arccos(x·y) → π/2``.  Random points do **not** tend toward ``α = π``.
Measured, the fraction with ``α > 0.99π`` is ``3.2e-04`` on ``S²`` and
indistinguishable from zero by ``n = 10``.

That sharpens what this repository's antipodal identification is claiming.
Selecting ``x ↔ −x`` is not "pairing a point with the far one".  It picks a
vanishing-measure relation out of an enormous set of nearly orthogonal
alternatives — and it gets *more* non-generic as dimension rises.

3 · the wavefront is a point, then an equatorial shell, then a point
────────────────────────────────────────────────────────────────────
A wave launched at a pole of ``S^n_R`` has transverse measure
``𝒜(χ) = |S^{n−1}|·(R sin χ)^{n−1}``: zero at the pole, maximal at ``χ = π/2``,
zero again at the antipode.  Higher ``n`` makes the middle overwhelmingly
dominant while the two poles get relatively tinier — so the dimensional-descent
picture *strengthens* with dimension, and with it the contrast between where
the propagation measure lives and where the congruence reconverges.

4 · the embedding centre is a whole direction space collapsed to a point
───────────────────────────────────────────────────────────────────────
An intrinsic ``S^n`` has no centre; embedded in ``ℝ^{n+1}`` it has exactly one,
and that point is the unique fixed point of the whole rotation group
``SO(n+1)``.  Every direction ``n̂ ∈ S^n`` gives a radial line through it, so the
origin is where an ``S^n``'s worth of directions coincide.

**That is a reading of what `regularized_center` does.**  Replacing ``f = 0``
with ``f₀ > 0`` does not merely avoid a singularity — it restores the direction
space the origin had crushed, at finite size.  The bearing *is* the blown-up
direction space, and `measure_the_bearing_is_the_blown_up_direction_space` puts
a number on how much of it comes back.

5 · and the collapse is ``f ⁿ``, not ``f``
──────────────────────────────────────────
For ``dℓ² = ds² + f(s)²dΩ_n²`` the transverse measure of an angular patch is
``f ⁿ dΩ_n``, so shrinking from ``F`` to ``f₀`` costs

    ``(f₀/F)ⁿ`` ,

which at ``f₀/F = 1e-03`` is ``1e-03`` on a circle, ``1e-06`` on ``S²``,
``1e-09`` on ``S³``.  The 2-D drawing's ``ℓ ∝ f`` badly understates it.

> **The angular overlap can stay finite while the physical overlap collapses
> as ``f₀ⁿ``.**

6 · the antipodal quotient flips orientability with parity — and the
    repository uses two quotients that are always on opposite sides
────────────────────────────────────────────────────────────────────
``ℝP^n`` is orientable **iff ``n`` is odd**, because the antipodal map is the
restriction of ``−I`` on ``ℝ^{n+1}`` and ``det(−I) = (−1)^{n+1}``.  So ``ℝP²``
is non-orientable, ``ℝP³`` orientable, ``ℝP⁴`` non-orientable: a parity effect,
not a size effect.

The refinement that matters here is that this repository carries **two**
antipodal quotients, at consecutive ``n``:

* the **spatial** slice, ``S^d/± = ℝP^d`` — at ``d = 3``, ``ℝP³``, *orientable*;
* the **two-body exchange** configuration space,
  ``(ℝ^d∖0)/± ≃ ℝP^{d−1} × ℝ₊`` — at ``d = 3``, ``ℝP²``, *non-orientable*, and
  it is the ``ℝP²`` whose Pin⁻ structure makes the throat a spinor.

Being one apart, they are **always of opposite parity**, so raising the spatial
dimension by one *swaps which of them is non-orientable*.  At ``d = 4`` the
exchange space would be the orientable ``ℝP³``, and the structure the repo's
spin-statistics result rests on would have to be re-derived rather than carried
over.  ``d = 2`` is different again: the exchange space is ``ℝP¹ ≃ S¹`` with
``π₁ = ℤ``, the braid group, not ``ℤ₂``.

7 · and ``S³`` is not a generic sphere, so nothing here extrapolates smoothly
────────────────────────────────────────────────────────────────────────────
``S³ ≅ SU(2)``, the unit quaternions; it is **parallelizable** (a global
orthonormal frame ``q·i, q·j, q·k``, verified to ``6.7e-16``), which only
``S¹``, ``S³`` and ``S⁷`` are; and it carries the Hopf fibration
``S¹ → S³ → S²``.  ``S²`` admits no nowhere-zero vector field at all.

So dimensions do not simply add room.  Some are algebraically special, and the
one this repository is built on is one of them — which is a standing argument
against reading any of §1–§6 as a smooth trend to be extrapolated.

Scope
─────
Every measurement here is about **measure, orientation and frames on round
spheres**.  Nothing solves a field equation, nothing evolves, and nothing here
claims the model should move to a different dimension — §6 in particular is a
statement about what *would* have to be redone, not a recommendation.
"""

from __future__ import annotations

import math
from typing import Dict, List, Tuple

import numpy as np
from scipy.integrate import quad

__all__ = [
    "ball_volume",
    "sphere_area",
    "shell_measure",
    "front_measure",
    "patch_collapse",
    "projective_space_is_orientable",
    "quaternion_frame",
    "hopf_map",
    "measure_the_peak_is_the_ball_not_the_sphere_and_needs_a_unit",
    "measure_the_shell_measure_concentrates_at_the_equator",
    "measure_the_antipode_is_a_vanishing_measure_relation",
    "measure_the_front_is_a_point_then_a_shell_then_a_point",
    "measure_the_bearing_is_the_blown_up_direction_space",
    "measure_the_patch_collapses_as_f_to_the_n",
    "measure_the_quotient_flips_orientability_with_parity",
    "measure_s3_is_not_a_generic_sphere",
]


# ════════════════════════════════════════════════════════════════════════════
# CLOSED FORMS
# ════════════════════════════════════════════════════════════════════════════
def ball_volume(dim: int, radius: float = 1.0) -> float:
    """``V_d = π^{d/2} R^d / Γ(d/2 + 1)`` — the ``d``-**ball**, peaking at ``d=5``."""
    d = int(dim)
    return math.pi ** (d / 2.0) * float(radius) ** d / math.gamma(d / 2.0 + 1.0)


def sphere_area(dim: int, radius: float = 1.0) -> float:
    """``A_{d−1} = 2π^{d/2} R^{d−1} / Γ(d/2)`` — the **sphere** ``S^{d−1} ⊂ ℝ^d``.

    Peaks at ``d = 7``, i.e. at ``S⁶``, two dimensions later than the ball.
    """
    d = int(dim)
    return 2.0 * math.pi ** (d / 2.0) * float(radius) ** (d - 1) \
        / math.gamma(d / 2.0)


def shell_measure(chi: float, sphere_dim: int) -> float:
    """``sin^{n−1}χ`` — the unnormalised measure of the shell at angle ``χ``
    from a pole of ``S^n``."""
    return math.sin(float(chi)) ** (int(sphere_dim) - 1)


def front_measure(chi: float, sphere_dim: int, radius: float = 1.0) -> float:
    """``|S^{n−1}|·(R sin χ)^{n−1}`` — a wavefront's transverse measure.

    Zero at the pole, maximal at ``χ = π/2``, zero again at the antipode.
    """
    n = int(sphere_dim)
    return sphere_area(n) * (float(radius) * math.sin(float(chi))) ** (n - 1)


def patch_collapse(neck: float, scale: float, angular_dim: int) -> float:
    """``(f₀/F)^n`` — what a fixed angular patch's *physical* measure does.

    The transverse measure of an angular patch in ``ds² + f²dΩ_n²`` is
    ``fⁿ dΩ_n``, so the angular footprint is untouched while the physical one
    falls as the ``n``-th power.  This is the same ``f^q`` weighting that
    `regularized_center.monopole_resistance` carries, followed through to the
    overlap region rather than to the resistance.
    """
    return (float(neck) / float(scale)) ** int(angular_dim)


def projective_space_is_orientable(dim: int) -> bool:
    """``ℝP^n`` is orientable iff ``n`` is odd.

    The antipodal map on ``S^n`` is the restriction of ``−I`` on ``ℝ^{n+1}``;
    orienting the sphere by (outward normal, tangent frame), the induced frame
    is ``−I`` applied to the whole ``(n+1)``-frame, so the orientation factor
    is ``det(−I_{n+1}) = (−1)^{n+1}``.
    """
    return int(dim) % 2 == 1


def quaternion_frame(q: np.ndarray) -> np.ndarray:
    """``(q·i, q·j, q·k)`` — a global orthonormal tangent frame on ``S³``.

    Nowhere zero, which is what "parallelizable" means and what ``S²`` cannot
    do.  Returned as ``(N, 3, 4)``: three tangent vectors per point.
    """
    q = np.atleast_2d(np.asarray(q, float))

    def mul(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        w1, x1, y1, z1 = a.T
        w2, x2, y2, z2 = b.T
        return np.stack([w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
                         w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
                         w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
                         w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2], axis=1)

    units = [np.tile(u, (len(q), 1)) for u in
             ([0., 1., 0., 0.], [0., 0., 1., 0.], [0., 0., 0., 1.])]
    return np.stack([mul(q, u) for u in units], axis=1)


def hopf_map(q: np.ndarray) -> np.ndarray:
    """``S³ → S²`` with ``S¹`` fibres, in complex coordinates.

    Writing ``q = (w,x,y,z)`` as ``z₁ = w + ix``, ``z₂ = y + iz``, the map is
    ``(2Re z₁z̄₂, 2Im z₁z̄₂, |z₁|² − |z₂|²)``, and multiplying *both* complex
    components by a common phase leaves all three components alone — so the
    fibre is exactly the circle.

    Worth writing in this form rather than as a quaternion sandwich: the
    obvious quaternionic expression is easy to get wrong by pairing the wrong
    components, and the first draft of this module did (the image missed the
    sphere by ``0.96``).
    """
    q = np.atleast_2d(np.asarray(q, float))
    z1 = q[:, 0] + 1j * q[:, 1]
    z2 = q[:, 2] + 1j * q[:, 3]
    p = z1 * np.conj(z2)
    return np.stack([2.0 * p.real, 2.0 * p.imag,
                     np.abs(z1) ** 2 - np.abs(z2) ** 2], axis=1)


# ════════════════════════════════════════════════════════════════════════════
# MEASUREMENTS
# ════════════════════════════════════════════════════════════════════════════
def measure_the_peak_is_the_ball_not_the_sphere_and_needs_a_unit(
) -> Dict[str, object]:
    """The usual "spheres peak at 5D" is two claims, and both need care.

    It is the **ball** whose volume peaks at ``d = 5``.  The **sphere**'s
    surface measure peaks two dimensions later, at ``d = 7`` — the sphere being
    ``S⁶ ⊂ ℝ⁷``.

    And neither is a fact about dimension by itself.  ``V_d`` and ``A_{d−1}``
    have different physical units at different ``d``, so comparing them across
    ``d`` picks a length scale, and the peak follows it: at ``R = 0.5`` the
    ball peaks at ``d = 1``, at ``R = 2`` at ``d = 24``.  The "peak at 5" is a
    statement about the unit ball and nothing more.
    """
    rows = [{"dim": d, "ball_volume": ball_volume(d),
             "sphere_area": sphere_area(d), "sphere": f"S^{d - 1}"}
            for d in range(1, 13)]
    # In log space, and over a range wide enough that the answer is a peak
    # rather than the end of the search: at R = 4 the ball peaks near d = 130,
    # so a range of 80 would have reported its own upper bound as the answer.
    def log_ball(d: int, r: float) -> float:
        return (d / 2.0) * math.log(math.pi) + d * math.log(r) \
            - math.lgamma(d / 2.0 + 1.0)

    def log_sphere(d: int, r: float) -> float:
        return math.log(2.0) + (d / 2.0) * math.log(math.pi) \
            + (d - 1) * math.log(r) - math.lgamma(d / 2.0)

    span = range(1, 2000)

    def argmax_with_ties(f, tol: float = 1e-12) -> Tuple[int, List[int]]:
        """The peak, and every dimension tied with it.

        Reported rather than resolved: at ``R = 1/2`` the sphere measure is
        *exactly* tied between ``d = 2`` and ``d = 3`` (``2πR = 4πR²`` at
        ``R = 1/2``), and a bare ``max`` picks a side by floating-point
        accident.  This arc has been bitten by a silent maximiser degeneracy
        before.
        """
        best = max(f(d) for d in span)
        tied = [d for d in span if abs(f(d) - best) <= tol * max(abs(best), 1.0)]
        return tied[0], tied

    scans = []
    for r in (0.5, 1.0, 2.0, 4.0):
        bp, bt = argmax_with_ties(lambda d: log_ball(d, r))
        sp, st = argmax_with_ties(lambda d: log_sphere(d, r))
        scans.append({"radius": r, "ball_peak_dim": bp, "sphere_peak_dim": sp,
                      "ball_tied_dims": bt, "sphere_tied_dims": st,
                      "hit_the_search_ceiling": bool(max(bp, sp) >= 1999)})
    unit = next(s for s in scans if s["radius"] == 1.0)
    return {
        "rows": rows,
        "radius_scan": scans,
        "ball_peaks_at": unit["ball_peak_dim"],
        "sphere_peaks_at": unit["sphere_peak_dim"],
        "the_peaking_sphere_is": f"S^{unit['sphere_peak_dim'] - 1} in R^"
                                 f"{unit['sphere_peak_dim']}",
        "the_ball_peaks_at_five": bool(unit["ball_peak_dim"] == 5),
        "the_sphere_peaks_two_later": bool(
            unit["sphere_peak_dim"] == unit["ball_peak_dim"] + 2),
        "the_unit_peaks_are_unique": bool(len(unit["ball_tied_dims"]) == 1
                                          and len(unit["sphere_tied_dims"]) == 1),
        "the_peak_moves_with_the_unit": bool(
            len({s["ball_peak_dim"] for s in scans}) > 1),
        "no_peak_is_at_the_search_ceiling": bool(
            not any(s["hit_the_search_ceiling"] for s in scans)),
        "ties_are_reported_not_resolved": "at R = 1/2 the sphere measure is "
                                          "exactly tied between d = 2 and "
                                          "d = 3, since 2 pi R = 4 pi R^2 "
                                          "there; a bare max picks a side by "
                                          "floating-point accident",
        "why": "V_d and A_{d-1} carry different units at different d, so "
               "comparing them across d implicitly chooses a length scale",
        "so_the_familiar_claim_is": "about the UNIT ball, and is not a "
                                    "dimension-invariant fact",
    }


def measure_the_shell_measure_concentrates_at_the_equator(
) -> Dict[str, object]:
    """``sin^{n−1}χ`` collapses onto ``χ = π/2`` with width ``1/√n``.

    Putting ``χ = π/2 + δ`` gives ``sin^{n−1}χ ≈ e^{−(n−1)δ²/2}``, a Gaussian
    of width ``1/√(n−1)``.  The measured ``std(χ)·√n`` runs
    ``0.967 → 0.999 → 1.000000`` over ``n = 2 … 1000``, which is the claim
    rather than an illustration of it.

    **Almost every point of a high-dimensional sphere is nearly 90° from any
    chosen point.**  Not scattered — concentrated in a shrinking band.
    """
    rows = []
    for n in (2, 3, 4, 7, 15, 50, 200, 1000):
        w = lambda c, n=n: shell_measure(c, n)
        z = quad(w, 0.0, math.pi, limit=400)[0]
        mean = quad(lambda c: c * w(c), 0.0, math.pi, limit=400)[0] / z
        var = quad(lambda c: (c - mean) ** 2 * w(c), 0.0, math.pi,
                   limit=400)[0] / z
        sd = math.sqrt(var)
        band = quad(w, math.pi / 2 - 1.0 / math.sqrt(n),
                    math.pi / 2 + 1.0 / math.sqrt(n), limit=400)[0] / z
        rows.append({"sphere_dim": n, "mean_chi": mean, "std_chi": sd,
                     "std_times_sqrt_n": sd * math.sqrt(n),
                     "mass_within_one_over_sqrt_n": band})
    tail = [r for r in rows if r["sphere_dim"] >= 50]
    return {
        "rows": rows,
        "law": "shell measure ~ sin^(n-1)(chi); chi = pi/2 + delta gives "
               "exp(-(n-1) delta^2 / 2), so the band width is ~1/sqrt(n)",
        "the_mean_is_always_the_equator": bool(
            all(abs(r["mean_chi"] - math.pi / 2) < 1e-9 for r in rows)),
        "std_times_sqrt_n_tends_to_one": bool(
            all(abs(r["std_times_sqrt_n"] - 1.0) < 1e-3 for r in tail)),
        "the_band_narrows_as_one_over_sqrt_n": bool(
            all(a["std_chi"] > b["std_chi"] for a, b in zip(rows, rows[1:]))),
        "what_it_means": "almost every point is nearly 90 degrees from any "
                         "chosen point -- concentrated, not scattered",
    }


def measure_the_antipode_is_a_vanishing_measure_relation(
        samples: int = 200000, seed: int = 7) -> Dict[str, object]:
    """Random directions tend to ``π/2``, so ``x ↔ −x`` is extremely special.

    ``x·y`` has mean ``0`` and width ``1/√n`` for random unit vectors, so the
    typical separation is a right angle and the near-antipodal fraction
    collapses: ``3.2e-04`` on ``S²``, and nothing at all in ``2e+05`` samples
    by ``n = 10``.

    **That is a point in favour of the identification being interesting.**
    Selecting ``x ↔ −x`` is not "pairing a point with the far one" — the far
    ones are a vanishing fraction of an enormous nearly-orthogonal majority,
    and the fraction shrinks as dimension grows.
    """
    rng = np.random.default_rng(int(seed))
    rows = []
    for n in (3, 4, 10, 50, 500):
        a = rng.normal(size=(int(samples), n))
        a /= np.linalg.norm(a, axis=1, keepdims=True)
        b = rng.normal(size=(int(samples), n))
        b /= np.linalg.norm(b, axis=1, keepdims=True)
        dot = np.einsum("ij,ij->i", a, b)
        ang = np.arccos(np.clip(dot, -1.0, 1.0))
        rows.append({
            "ambient_dim": n,
            "std_dot": float(dot.std()),
            "std_dot_times_sqrt_n": float(dot.std() * math.sqrt(n)),
            "mean_angle_over_pi": float(ang.mean() / math.pi),
            "fraction_far_from_right_angle": float(
                np.mean(np.abs(ang - math.pi / 2) > 0.2)),
            "fraction_near_antipodal": float(np.mean(ang > 0.99 * math.pi)),
        })
    return {
        "rows": rows,
        "law": "x.y has mean 0 and width 1/sqrt(n), so alpha -> pi/2",
        "std_scales_as_one_over_sqrt_n": bool(
            all(abs(r["std_dot_times_sqrt_n"] - 1.0) < 0.02 for r in rows)),
        "the_mean_angle_is_a_right_angle": bool(
            all(abs(r["mean_angle_over_pi"] - 0.5) < 0.01 for r in rows)),
        "near_antipodal_fraction_collapses": bool(
            rows[0]["fraction_near_antipodal"]
            > 10.0 * rows[2]["fraction_near_antipodal"]
            or rows[2]["fraction_near_antipodal"] == 0.0),
        "so_the_antipodal_identification_is_non_generic":
            "it selects a vanishing-measure relation out of an enormous set "
            "of nearly orthogonal alternatives, and gets MORE non-generic as "
            "dimension rises",
        "what_is_not_claimed": "that this makes the identification correct -- "
                               "only that it is not the bland 'pair it with "
                               "the far point' it can sound like",
    }


def measure_the_front_is_a_point_then_a_shell_then_a_point(
) -> Dict[str, object]:
    """``𝒜(χ) ∝ sin^{n−1}χ``: the descent picture strengthens with dimension.

    A wave launched at a pole has zero transverse measure there, maximum at
    the equator, and zero again at the antipode.  Raising ``n`` makes the
    middle overwhelmingly dominant — the share of the path-integrated measure
    within ``0.1`` rad of the equator runs ``10.0% → 12.7% → 20.2%`` from
    ``S²`` to ``S⁷`` — while the poles get relatively tinier.

    So the contrast the model turns on gets sharper, not softer: **almost all
    of the propagation measure sits around the equator, yet the whole radial
    congruence still reconverges at one antipodal locus.**
    """
    rows = []
    grid = np.linspace(1e-9, math.pi - 1e-9, 20001)
    for n in (2, 3, 4, 7, 15):
        f = lambda c, n=n: shell_measure(c, n)
        z = quad(f, 0.0, math.pi, limit=400)[0]
        peak = float(max(grid, key=f))
        core = quad(f, math.pi / 2 - 0.1, math.pi / 2 + 0.1, limit=400)[0] / z
        rows.append({
            "sphere_dim": n,
            "front_exponent": n - 1,
            "peak_chi_over_pi": peak / math.pi,
            "measure_within_a_tenth_of_the_equator": core,
            "front_at_quarter_turn_over_peak": front_measure(math.pi / 4, n)
            / front_measure(math.pi / 2, n),
        })
    return {
        "rows": rows,
        "law": "front measure = |S^(n-1)| (R sin chi)^(n-1)",
        "the_peak_is_always_the_equator": bool(
            all(abs(r["peak_chi_over_pi"] - 0.5) < 1e-4 for r in rows)),
        "the_poles_vanish": bool(front_measure(0.0, 3) == 0.0
                                 and front_measure(math.pi, 3) < 1e-30),
        "the_middle_dominates_more_with_dimension": bool(
            all(a["measure_within_a_tenth_of_the_equator"]
                < b["measure_within_a_tenth_of_the_equator"]
                for a, b in zip(rows, rows[1:]))),
        "the_contrast_the_model_turns_on": "almost all propagation measure is "
                                           "equatorial, yet the whole radial "
                                           "congruence reconverges at one "
                                           "antipodal locus",
    }


def measure_the_bearing_is_the_blown_up_direction_space(
) -> Dict[str, object]:
    """A reading of PR #268: the bearing *is* the direction space, restored.

    An embedded ``S^n_R ⊂ ℝ^{n+1}`` has exactly one centre, and it is the
    unique fixed point of ``SO(n+1)``.  Every direction ``n̂ ∈ S^n`` runs a
    radial line through it, so the origin is where an entire ``S^n``'s worth of
    directions coincide — *the direction space crushed to zero size*.

    Regularising ``f = 0 → f₀ > 0`` restores it at finite size.  The bearing's
    own measure is ``f₀ⁿ·|S^n|``, so how much comes back is a strong function
    of ``n``: at ``f₀ = 1e-03`` it is ``6.3e-03`` for a circle and ``2.0e-08``
    for an ``S³`` — nonzero in every case, which is the whole point, and
    vanishingly small in a way the 2-D drawing cannot suggest.

    This is offered as an *interpretation* of a construction that stands on
    its own geometry, not as a new result.  Nothing below is used by
    `regularized_center`.
    """
    rows = []
    for n in (1, 2, 3, 4):
        for f0 in (1e-2, 1e-3):
            rows.append({
                "angular_dim": n, "neck": f0,
                "unit_direction_space_measure": sphere_area(n + 1),
                "bearing_measure": sphere_area(n + 1) * f0 ** n,
                "fraction_of_unit_scale": f0 ** n,
            })
    return {
        "rows": rows,
        "the_centre_is_unique": "|O - X| = R for every surface point, and O is "
                                "the unique fixed point of SO(n+1)",
        "what_the_origin_loses": "an entire S^n of directions collapses to one "
                                 "invariant point",
        "what_the_bearing_restores": "that direction space at finite size, "
                                     "with measure f0^n |S^n|",
        "every_bearing_measure_is_positive": bool(
            all(r["bearing_measure"] > 0.0 for r in rows)),
        "and_shrinks_faster_in_higher_dimension": bool(
            all(a["fraction_of_unit_scale"] > b["fraction_of_unit_scale"]
                for a, b in zip(
                    [r for r in rows if r["neck"] == 1e-3],
                    [r for r in rows if r["neck"] == 1e-3][1:]))),
        "this_is_an_interpretation_not_a_result":
            "PR #268's bearing stands on its own geometry; this section only "
            "says what the construction can be read as doing, and nothing in "
            "regularized_center depends on it",
    }


def measure_the_patch_collapses_as_f_to_the_n() -> Dict[str, object]:
    """**The correction that matters for the picture.**  ``(f₀/F)ⁿ``, not ``f₀/F``.

    For ``dℓ² = ds² + f²dΩ_n²`` the transverse measure of an angular patch is
    ``fⁿ dΩ_n``.  So a packet keeps **exactly the same angular footprint**
    while its physical overlap falls as the ``n``-th power of the squeeze:

        ``f₀/F = 1e-03``  →  ``1e-03`` (``S¹``), ``1e-06`` (``S²``),
        ``1e-09`` (``S³``), ``1e-12`` (``S⁴``).

    PR #268's review had already forced the *size* law to be scoped by
    dimension; this is that scoping followed to its conclusion.  It changes
    what the finite-bearing picture is a picture *of*: much less like two
    ribbons squeezing together, much more like a vast angular configuration
    space packed into an extremely small proper region.

    The yes/no overlap criterion is untouched — it was always angular.
    """
    rows = []
    for ratio in (1e-1, 1e-2, 1e-3):
        row = {"squeeze": ratio}
        for n in (1, 2, 3, 4):
            row[f"n_{n}"] = patch_collapse(ratio, 1.0, n)
        rows.append(row)
    drawn, physical = [], []
    for n in (1, 2, 3, 4):
        drawn.append(patch_collapse(1e-3, 1.0, 1))
        physical.append(patch_collapse(1e-3, 1.0, n))
    return {
        "rows": rows,
        "law": "transverse measure of an angular patch = f^n dOmega_n, so a "
               "squeeze from F to f0 costs (f0/F)^n",
        "the_angular_footprint_is_untouched": True,
        "the_drawn_2d_length_understates_it": bool(
            physical[-1] < 1e-6 * drawn[-1]),
        "understatement_factor_at_n_three": float(
            patch_collapse(1e-3, 1.0, 1) / patch_collapse(1e-3, 1.0, 3)),
        "the_collapse_steepens_with_dimension": bool(
            all(a > b for a, b in zip(physical, physical[1:]))),
        "what_the_picture_is_of": "not two ribbons squeezing together, but a "
                                  "vast angular configuration space packed "
                                  "into an extremely small proper region",
        "the_yes_no_criterion_is_unaffected": "whether two fronts meet was "
                                              "always an angular question, "
                                              "and stays one",
    }


def measure_the_quotient_flips_orientability_with_parity(
) -> Dict[str, object]:
    """``ℝP^n`` orientable iff ``n`` odd — and this repo uses two, one apart.

    ``det(−I_{n+1}) = (−1)^{n+1}``, so the antipodal map preserves orientation
    exactly when ``n`` is odd.  A **parity** effect, not a size effect, and one
    that no amount of "the next sphere up" intuition would suggest.

    The consequence for this repository is specific.  It carries two antipodal
    quotients at *consecutive* ``n``:

    * the **spatial** slice ``S^d/± = ℝP^d`` — at ``d = 3``, orientable;
    * the **two-body exchange** space ``(ℝ^d∖0)/± ≃ ℝP^{d−1} × ℝ₊`` — at
      ``d = 3``, the non-orientable ``ℝP²`` whose Pin⁻ structure is what makes
      the throat a spinor.

    One apart means **always opposite parity**, so raising the spatial
    dimension by one swaps which is which: at ``d = 4`` the exchange space is
    the orientable ``ℝP³``, and the spin-statistics mechanism would have to be
    re-derived rather than carried across.  ``d = 2`` fails differently —
    ``ℝP¹ ≃ S¹`` has ``π₁ = ℤ``, the braid group.

    Stated as a consequence *if* the dimension moved.  It is not an argument
    that it should.
    """
    rows = []
    for d in (2, 3, 4, 5, 6):
        rows.append({
            "spatial_dim": d,
            "spatial_quotient": f"RP^{d}",
            "spatial_orientable": projective_space_is_orientable(d),
            "exchange_quotient": f"RP^{d - 1}",
            "exchange_orientable": projective_space_is_orientable(d - 1),
            "exchange_pi1": "Z (braid)" if d == 2 else "Z_2",
            "opposite_parity": bool(projective_space_is_orientable(d)
                                    != projective_space_is_orientable(d - 1)),
        })
    dets = [{"n": n, "det_minus_identity": float(
        round(np.linalg.det(-np.eye(n + 1)))),
        "orientable": projective_space_is_orientable(n)} for n in range(1, 8)]
    here = next(r for r in rows if r["spatial_dim"] == 3)
    return {
        "rows": rows,
        "determinants": dets,
        "why": "the antipodal map is the restriction of -I on R^(n+1), so the "
               "orientation factor is det(-I) = (-1)^(n+1)",
        "orientable_iff_odd": bool(
            all(d["orientable"] == (d["det_minus_identity"] > 0)
                for d in dets)),
        "the_two_quotients_are_always_opposite": bool(
            all(r["opposite_parity"] for r in rows)),
        "at_the_repos_dimension": {
            "spatial": f"{here['spatial_quotient']} orientable",
            "exchange": f"{here['exchange_quotient']} NON-orientable "
                        f"-- this is where the Pin- structure lives"},
        "raising_the_dimension_swaps_them": bool(
            rows[1]["spatial_orientable"] != rows[2]["spatial_orientable"]),
        "it_is_parity_not_size": True,
        "what_is_not_claimed": "that the model should move dimension -- this "
                               "says what would have to be re-derived if it "
                               "did, nothing more",
    }


def measure_s3_is_not_a_generic_sphere(samples: int = 20000,
                                       seed: int = 3) -> Dict[str, object]:
    """A standing argument against extrapolating any of the above smoothly.

    ``S³`` is not "the next sphere".  It is ``SU(2)``, the unit quaternions;
    it is **parallelizable**, and the frame ``(q·i, q·j, q·k)`` is verified
    orthonormal and tangent to ``6.7e-16`` at every sampled point; and it
    carries the Hopf fibration ``S¹ → S³ → S²``, checked here by confirming
    that the image lands on ``S²`` and that a common phase on both complex
    components moves the point but not its image.

    Only ``S¹``, ``S³`` and ``S⁷`` are parallelizable — ``S²`` admits no
    nowhere-zero vector field at all.  So dimension does not simply add room:
    some dimensions are algebraically special, and the one this repository is
    built on is one of them.
    """
    rng = np.random.default_rng(int(seed))
    q = rng.normal(size=(int(samples), 4))
    q /= np.linalg.norm(q, axis=1, keepdims=True)

    frame = quaternion_frame(q)
    full = np.concatenate([q[:, None, :], frame], axis=1)
    gram = np.einsum("nik,njk->nij", full, full)
    gram_err = float(np.abs(gram - np.eye(4)).max())
    norms = np.linalg.norm(frame, axis=2)

    h = hopf_map(q)
    on_sphere = float(np.abs(np.linalg.norm(h, axis=1) - 1.0).max())
    phase = np.exp(1j * rng.uniform(0.0, 2 * math.pi, len(q)))
    z1, z2 = (q[:, 0] + 1j * q[:, 1]) * phase, (q[:, 2] + 1j * q[:, 3]) * phase
    moved = np.stack([z1.real, z1.imag, z2.real, z2.imag], axis=1)
    fibre_err = float(np.abs(hopf_map(moved) - h).max())
    fibre_dist = float(np.median(np.linalg.norm(moved - q, axis=1)))

    return {
        "identifications": ["S^3 = SU(2) = the unit quaternions",
                            "Hopf fibration S^1 -> S^3 -> S^2",
                            "parallelizable (only S^1, S^3, S^7 are)"],
        "frame_orthonormality_error": gram_err,
        "smallest_frame_vector_norm": float(norms.min()),
        "the_frame_is_global_and_nowhere_zero": bool(
            gram_err < 1e-12 and norms.min() > 0.999),
        "hopf_image_on_the_sphere_error": on_sphere,
        "hopf_fibre_error": fibre_err,
        "median_distance_moved_along_a_fibre": fibre_dist,
        "the_fibre_is_a_circle_not_a_point": bool(
            fibre_err < 1e-12 and fibre_dist > 0.1),
        "s2_has_no_such_frame": "the hairy-ball theorem: S^2 admits no "
                                "nowhere-zero tangent vector field",
        "so_nothing_here_extrapolates_smoothly":
            "dimensions do not simply add room; S^1, S^3 and S^7 carry "
            "division-algebra structure that their neighbours do not, and "
            "this repository is built on one of them",
    }
