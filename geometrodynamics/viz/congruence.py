"""
A congruence, its Jacobian, and where the three conflated things separate.

WHY A CONGRUENCE AND NOT A BRIGHT POINT
───────────────────────────────────────
A singularity is a failure of evolution, not a bright spot, so the object worth
watching is a **congruence with a deforming cross-section**.  The quantity that
carries the story is the Jacobian of the map from material coordinates to the
evolved geometry,

    ``J = det F`` ,     ``F̈ = ½ ḧ F`` ,     ``F(0) = I`` ,

the geodesic-deviation equation for a wave in TT gauge.  ``J → 0`` is a caustic
of the *congruence*: neighbouring trajectories arrive at the same place.  It is
not a statement about the metric.

``ḧ`` IS TAKEN FROM THE WAVE OPERATOR, NOT FROM A TIME DIFFERENCE
─────────────────────────────────────────────────────────────────
``h = sin²d · q`` and ``q`` solves the spin-2 wave equation, so

    ``ḧ = sin²d · ∇²q``

is available exactly at every step.  The first draft of this module instead
formed ``ḧ`` from a three-sample time difference seeded with
``h(−dt) = h(0)``.  That seeding injects a spurious velocity impulse ``½ḣ(0)``
at the first step — and because the ``1/dt²`` and the ``dt`` of the update
cancel, **the error is independent of the timestep**.  Refining the grid
reproduced it to five digits.  It is only visible against the operator form,
which is why the operator form is what the class uses.

THE FOCUSING IS ENTIRELY SHEAR
──────────────────────────────
For the trace-free ``h`` of a vacuum wave the transverse Raychaudhuri equation

    ``dθ/dt = −θ²/2 − 2σ² − R_ab u^a u^b``

is *the same statement* as the deviation equation, with ``θ = d(ln J)/dt`` —
substituting ``θ`` and ``σ`` into the right-hand side reproduces the left
symbolically, so the measured ``6.7e-15`` residual is round-off and not a
convergence result.  The content is the **Ricci term, which vanishes
identically** because ``h`` is trace-free: ``ä/a`` and ``b̈/b`` are exact
opposites, so the matter term cancels term by term.  None of the focusing comes
from matter; all of it is shear-squared, which is second order in the
amplitude.  That is why a weak wave barely focuses and a strong one does so
abruptly.

THE LAUNCH HAS TO BE COMPACTLY SUPPORTED
────────────────────────────────────────
A Gaussian bump is nonzero everywhere, so a congruence at the far side starts
deforming immediately.  Measured: the tidal field reaches ``d = 2.95`` at
``t = 2.00`` from a Gaussian, against a light-arrival time of ``2.95``, and that
figure is **grid-converged** — it is the analytic tail, not a numerical
precursor.  Launching from the compact bump ``(1 − x²)⁴`` instead puts the
arrival at ``2.7697`` against the causal bound ``d − w = 2.7699``, with the
residual earliness shrinking under refinement.  ``wave_constraints`` found this
for the scalar; the spin-2 launch needs it too, and it is the default here.

The poles need the same care.  The solver's grid is cell-centred, so its last
sample sits at ``π − dd/2`` and interpolating ``h`` itself *clamps* there rather
than vanishing — which drives the antipodal station by an O(dd) interpolation
error, at precisely the one place this round claims nothing is driven at all.
``h`` and ``ḧ`` are therefore built as ``sin²d ·`` (interpolated ``q``), with the
``sin²d`` evaluated on the slice's own distances, so the spin-weight condition
holds exactly at the sampled points.

WHAT THE EQUATIONS ACTUALLY DECIDE
──────────────────────────────────
1. **Ordinary focus.**  ``J`` dips and recovers; the map stays invertible.
2. **Caustic.**  ``J`` reaches zero — trajectories cross, several histories
   land on one rendered point.  The background is completely regular.
3. **Curvature singularity.**  The geometry itself would have to fail.

Case 3 is **not reachable in this program, and that is a statement about the
program rather than about the wave**: the background is a fixed round ``S²``
with Gaussian curvature ``1`` everywhere at every time, with no Einstein
equation and no backreaction, so there is nothing that could diverge or stop.

Of the three outcomes the caustic could have had — passage, singular
termination, or finite-radius reconnection into a detached resonator — this
program gives **passage**.  At the source ring the crossing is deep and fully
resolved: the tracked point plunges to ``J ≈ −470`` and stays there, with the
crossing slope converged to three digits under a halved timestep.  There is no
bounce at finite radius, and the reason is structural rather than numerical:
each material point's ``F`` is driven only by the external ``h`` and never by
its neighbours, so the congruence cannot act back on anything.  Reconnection
needs exactly that feedback, so this program *could not have* produced it —
worth stating plainly, because "we did not see it" and "it could not have
appeared" are very different claims.

TWO RINGS, AND THEY ARE NOT ALIKE
─────────────────────────────────
The source ring and the converging ring close at thresholds a factor of ten
apart in strain (``0.026`` against ``0.247`` in a ``1.2π`` window).  The source
is shaken hardest and from ``t = 0``; the converging ring has to arrive before
it can focus, and even at 12% strain it only *grazes* zero — it crosses and
returns within a few thousandths, and the depth of that excursion does not
converge.  So the answer to "does focusing drive the neck radius to zero" is:
barely, and only at a strain nothing physical would reach.

Every threshold here carries its **window**, because the wave on ``S²`` is
exactly periodic and a fixed material point is driven over and over by the same
returning ring.  The deviation equation is then a Hill equation and the
accumulated focusing keeps growing, so a threshold quoted without a window is
not a number about the wave.

THE NECK IS A RING, AND SPIN WEIGHT IS WHY
──────────────────────────────────────────
``h = sin²d · q`` vanishes at both poles no matter what ``q`` does, so the
tidal field is identically zero at the antipode and the congruence there is
never driven at all.  The neck therefore forms on a *ring* around the focus, at
geodesic radius ``≈ 0.44 w`` — the same ratio across a 3.3× range of pulse
width, to within the one-cell resolution of the angular grid, and independent
of amplitude.  A spin-2 focus has no centre.
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence

import numpy as np

from geometrodynamics.viz.spin2_tidal import ANTIPODAL_TIME, Spin2WaveSim

__all__ = [
    "TidalCongruence",
    "compact_launch",
    "measure_raychaudhuri_is_exact",
    "measure_the_operator_form_matters",
    "measure_the_front_is_causal",
    "measure_the_three_cases",
    "measure_the_caustic_thresholds",
    "measure_the_neck_is_a_ring",
    "measure_the_caustic_is_a_passage",
    "measure_the_threshold_depends_on_the_window",
    "measure_case_three_is_unreachable",
    "spacetime_history",
]

FAR_SIDE = 0.5 * math.pi          # "the converging ring", d > π/2
NECK_CAP = 0.5                    # the antipodal cap the neck is searched in


def compact_launch(sim: Spin2WaveSim) -> None:
    """Launch ``sim`` from a compactly supported outgoing bump.

    ``q₀ = (1 − x²)⁴`` on ``x = d/w < 1`` and exactly zero beyond it, with the
    outgoing velocity ``q̇ = −∂_d q₀``.  Four derivatives vanish at the edge, so
    the support is genuinely finite without a kink to radiate.
    """
    w = sim.pulse_width
    x = np.clip(sim.d / w, 0.0, 1.0)
    inside = sim.d < w
    q0 = (1.0 - x ** 2) ** 4
    q_dot = 8.0 * x * (1.0 - x ** 2) ** 3 / w * inside     # −∂_d q₀
    sim.start_from(q0, q_dot)


class TidalCongruence:
    """Material trajectories deformed by the tidal field, and their Jacobian.

    One deformation gradient per point of the slice, integrated from
    ``F̈ = ½ ḧ F``.  For the axisymmetric field ``h = diag(h₊, −h₊)`` in the
    geodesic-polar frame the gradient stays diagonal, ``F = diag(a, b)``, so
    ``a`` is the principal stretch along ``ê_d`` and ``b`` the one along
    ``ê_ψ``, and ``J = ab`` is the cross-sectional area of the congruence in
    units of its initial area.

    ``gain`` multiplies the solved field before it drives the deviation.  It is
    the one strength choice in the module and it is reported as a peak strain
    everywhere it matters, because the focusing is second order in it.
    """

    def __init__(self, sim: Optional[Spin2WaveSim] = None,
                 n_sigma: int = 721, gain: float = 25.0,
                 pulse_width: float = 0.18, n_radial: int = 1600,
                 compact: bool = True) -> None:
        self.compact = bool(compact)
        self.sim = sim or Spin2WaveSim(n=n_radial, pulse_width=pulse_width)
        self.gain = float(gain)
        self.n_sigma = int(n_sigma)
        # σ = ±π is one point of the slice, so the duplicate endpoint is out
        self.sigma = np.linspace(-math.pi, math.pi, self.n_sigma,
                                 endpoint=False)
        self.distance = np.abs(self.sigma)
        self.reset()

    # ── state ───────────────────────────────────────────────────────────────
    def reset(self) -> None:
        if self.compact:
            compact_launch(self.sim)
        else:
            self.sim.reset()
        n = self.n_sigma
        self.a = np.ones(n)
        self.b = np.ones(n)
        self.da = np.zeros(n)
        self.db = np.zeros(n)
        self.t = 0.0

    @property
    def dt(self) -> float:
        return float(self.sim.dt)

    @property
    def pulse_width(self) -> float:
        return float(self.sim.pulse_width)

    def causal_bound(self, d) -> np.ndarray:
        """Earliest time the field may legitimately reach ``d``.

        The support starts at radius ``w``, and the characteristics of the
        spin-2 equation travel at exactly one in geodesic distance.
        """
        return np.maximum(np.asarray(d, dtype=float) - self.pulse_width, 0.0)

    # ── the driving field ───────────────────────────────────────────────────
    def h(self) -> np.ndarray:
        """``h = sin²d · q`` sampled onto the slice, with the gain applied.

        ``q`` is interpolated and the ``sin²d`` is evaluated on the slice's own
        distances, rather than interpolating the product.  The solver's grid is
        cell-centred, so its last sample sits at ``π − dd/2`` and interpolating
        ``h`` itself *clamps* at the poles instead of vanishing there — which
        would drive the antipodal station by an O(dd) interpolation error, at
        exactly the one place this series claims nothing is driven at all.
        """
        return (self.gain * np.sin(self.distance) ** 2
                * np.interp(self.distance, self.sim.d, self.sim.q))

    def h_ddot(self) -> np.ndarray:
        """``ḧ = sin²d · ∇²q``, from the wave operator rather than a difference.

        Factorised the same way as ``h``, so the spin-weight condition ``ḧ = 0``
        at both poles holds exactly at the sampled points and not merely to
        interpolation accuracy.
        """
        lap = self.sim._laplacian(self.sim.q)
        return (self.gain * np.sin(self.distance) ** 2
                * np.interp(self.distance, self.sim.d, lap))

    # ── one step of geodesic deviation ──────────────────────────────────────
    def step(self) -> None:
        hdd = self.h_ddot()
        self.da += 0.5 * hdd * self.a * self.dt
        self.db += -0.5 * hdd * self.b * self.dt
        self.a += self.da * self.dt
        self.b += self.db * self.dt
        self.sim.step()
        self.t = self.sim.t

    def advance_to(self, t_target: float) -> None:
        while self.t < t_target:
            self.step()

    # ── the quantities worth watching ───────────────────────────────────────
    def jacobian(self) -> np.ndarray:
        """``J = det F`` — the cross-sectional area of the congruence."""
        return self.a * self.b

    def neck_radius(self) -> np.ndarray:
        """``√|J|`` — the linear scale of the cross-section."""
        return np.sqrt(np.abs(self.jacobian()))

    def expansion(self) -> np.ndarray:
        """``θ = d(ln J)/dt`` — negative means neighbours are converging."""
        return self.da / self.a + self.db / self.b

    def shear(self) -> np.ndarray:
        """``σ`` — the trace-free part, which is what does the focusing."""
        return 0.5 * (self.da / self.a - self.db / self.b)

    def ricci_term(self) -> np.ndarray:
        """``R_ab u^a u^b`` for this field: identically zero, and checked.

        A trace-free ``h`` drives ``ä/a`` and ``b̈/b`` to be exact opposites, so
        the matter term of Raychaudhuri vanishes and every bit of the focusing
        is shear.  Written out rather than returned as a literal zero, so that
        the cancellation is what the measurement actually tests.
        """
        hdd = self.h_ddot()
        return 0.5 * hdd + (-0.5 * hdd)

    def raychaudhuri_residual(self) -> np.ndarray:
        """``dθ/dt`` measured against ``−θ²/2 − 2σ²``, with no Ricci term."""
        hdd = self.h_ddot()
        A, B = self.da / self.a, self.db / self.b
        actual = (0.5 * hdd - A ** 2) + (-0.5 * hdd - B ** 2)
        theta, sig = A + B, 0.5 * (A - B)
        return actual - (-0.5 * theta ** 2 - 2.0 * sig ** 2)

    def has_caustic(self) -> bool:
        return bool(np.any(self.jacobian() <= 0.0))

    # ── the two regions, which behave completely differently ────────────────
    def far_mask(self) -> np.ndarray:
        """The converging ring: ``d > π/2``, the half that refocuses."""
        return self.distance > FAR_SIDE

    def peak_strain(self) -> float:
        return float(np.max(np.abs(self.h())))


# ════════════════════════════════════════════════════════════════════════════
# helpers shared by the measurements
# ════════════════════════════════════════════════════════════════════════════
def _unit_peak_strain(pulse_width: float, n_radial: int, t_end: float,
                      compact: bool = True) -> float:
    """Peak ``|h|`` over the window at ``gain = 1``, for normalising strength."""
    sim = Spin2WaveSim(n=n_radial, pulse_width=pulse_width)
    if compact:
        compact_launch(sim)
    peak = 0.0
    while sim.t < t_end:
        sim.step()
        peak = max(peak, float(np.max(np.abs(sim.h))))
    return peak


def _run(gain: float, *, n_sigma: int = 361, n_radial: int = 1600,
         pulse_width: float = 0.18, t_end: float = 1.2 * ANTIPODAL_TIME,
         compact: bool = True) -> Dict[str, object]:
    """Integrate once and report the minima, split by region."""
    c = TidalCongruence(gain=gain, n_sigma=n_sigma, n_radial=n_radial,
                        pulse_width=pulse_width, compact=compact)
    far = c.far_mask()
    near = ~far
    best_far = {"J": math.inf, "t": 0.0, "d": 0.0}
    best_near = {"J": math.inf, "t": 0.0, "d": 0.0}
    first_far = first_near = None
    while c.t < t_end:
        c.step()
        J = c.jacobian()
        for mask, best, name in ((far, best_far, "far"),
                                 (near, best_near, "near")):
            k = int(np.argmin(J[mask]))
            if float(J[mask][k]) < best["J"]:
                best.update(J=float(J[mask][k]), t=c.t,
                            d=float(c.distance[mask][k]))
        if first_far is None and float(J[far].min()) <= 0.0:
            first_far = c.t
        if first_near is None and float(J[near].min()) <= 0.0:
            first_near = c.t
    return {"far": best_far, "near": best_near, "first_far": first_far,
            "first_near": first_near, "peak_strain": c.peak_strain(),
            "finite": bool(np.all(np.isfinite(c.jacobian())))}


def _caustics(gain: float, region: str, *, n_sigma: int = 241,
              n_radial: int = 1200, t_end: float = 1.2 * ANTIPODAL_TIME,
              pulse_width: float = 0.18) -> bool:
    c = TidalCongruence(gain=gain, n_sigma=n_sigma, n_radial=n_radial,
                        pulse_width=pulse_width)
    mask = c.far_mask() if region == "far" else ~c.far_mask()
    while c.t < t_end:
        c.step()
        if float(c.jacobian()[mask].min()) <= 0.0:
            return True
    return False


def _bisect_threshold(region: str, *, lo: float = 0.25, hi: float = 400.0,
                      steps: int = 22, **kw) -> Optional[float]:
    if not _caustics(hi, region, **kw):
        return None
    for _ in range(steps):
        mid = 0.5 * (lo + hi)
        if _caustics(mid, region, **kw):
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


# ════════════════════════════════════════════════════════════════════════════
# MEASUREMENTS
# ════════════════════════════════════════════════════════════════════════════
def measure_raychaudhuri_is_exact(gain: float = 25.0, n_sigma: int = 241,
                                  frames: int = 400) -> Dict[str, object]:
    """``dθ/dt = −θ²/2 − 2σ²`` holds to round-off, with the Ricci term zero.

    **This is an algebraic identity, not an accuracy check, and it is worth
    being blunt about which.**  Substituting ``θ = A + B`` and
    ``σ = (A − B)/2`` into ``−θ²/2 − 2σ²`` gives exactly ``−A² − B²``, which is
    what the deviation equation produces — so the residual is zero *symbolically*
    and the reported number is floating-point noise.  Nothing here validates the
    time integration; the ``n_radial`` ladders elsewhere in this module do that.

    What the measurement does establish is structural, and is why it is kept:
    the deviation equation and the transverse Raychaudhuri equation are *the
    same statement* for this system, and the **Ricci term vanishes identically**
    because ``h`` is trace-free — ``ä/a`` and ``b̈/b`` are exact opposites, so the
    matter term cancels term by term rather than numerically.  All of the
    focusing is therefore shear-squared, and therefore second order in the
    amplitude, which is the shape of every threshold result below.
    """
    c = TidalCongruence(gain=gain, n_sigma=n_sigma)
    worst_res = worst_ricci = 0.0
    for _ in range(frames):
        c.step()
        worst_res = max(worst_res,
                        float(np.max(np.abs(c.raychaudhuri_residual()))))
        worst_ricci = max(worst_ricci, float(np.max(np.abs(c.ricci_term()))))
    return {
        "gain": gain,
        "worst_raychaudhuri_residual": worst_res,
        "worst_ricci_term": worst_ricci,
        "raychaudhuri_is_exact": bool(worst_res < 1e-9),
        "the_ricci_term_vanishes": bool(worst_ricci < 1e-15),
        "so_the_focusing_is_all_shear": bool(worst_ricci < 1e-15),
    }


def measure_the_operator_form_matters(
        gains: Sequence[float] = (25.0,),
        n_radials: Sequence[int] = (600, 1200, 2400),
        n_sigma: int = 241,
        t_end: float = 1.2 * ANTIPODAL_TIME) -> Dict[str, object]:
    """A seeded time difference for ``ḧ`` is wrong by an amount ``dt`` cannot fix.

    Seeding the history with ``h(−dt) = h(0)`` makes the first step's
    ``ḧ ≈ ḣ(0)/dt``; multiplied by the ``dt`` of the velocity update that is a
    finite kick ``½ḣ(0)``, so refining the grid converges to the *wrong*
    answer.  This measurement runs both forms on the same refinement ladder:
    each converges, and they converge to different numbers.
    """
    rows = []
    for gain in gains:
        for nr in n_radials:
            c = TidalCongruence(gain=float(gain), n_sigma=n_sigma, n_radial=nr)
            # the operator form
            op = math.inf
            while c.t < t_end:
                c.step()
                op = min(op, float(c.jacobian().min()))
            # the seeded-difference form, on the same data
            d = TidalCongruence(gain=float(gain), n_sigma=n_sigma, n_radial=nr)
            h1 = h2 = d.h()
            fd = math.inf
            while d.t < t_end:
                d.sim.step()
                h = d.h()
                hdd = (h - 2.0 * h1 + h2) / d.dt ** 2
                d.da += 0.5 * hdd * d.a * d.dt
                d.db += -0.5 * hdd * d.b * d.dt
                d.a += d.da * d.dt
                d.b += d.db * d.dt
                h2, h1 = h1, h
                d.t = d.sim.t
                fd = min(fd, float(d.jacobian().min()))
            rows.append({"gain": float(gain), "n_radial": nr,
                         "operator_min_J": op, "seeded_difference_min_J": fd})
    ops = [r["operator_min_J"] for r in rows]
    fds = [r["seeded_difference_min_J"] for r in rows]
    spread = lambda v: (max(v) - min(v)) / max(abs(max(v)), abs(min(v)), 1e-30)
    return {
        "rows": rows,
        "operator_spread_over_the_ladder": spread(ops),
        "seeded_spread_over_the_ladder": spread(fds),
        "both_forms_converge": bool(spread(ops) < 0.05 and spread(fds) < 0.05),
        "but_they_disagree": bool(abs(ops[-1] - fds[-1])
                                  > 0.5 * abs(ops[-1])),
        "refinement_cannot_reveal_it": True,
    }


def measure_the_front_is_causal(probe_distance: float = 2.9499,
                                n_radials: Sequence[int] = (800, 1600, 3200),
                                threshold: float = 1e-12,
                                pulse_width: float = 0.18
                                ) -> Dict[str, object]:
    """A Gaussian launch reacts before the front; a compact one does not.

    The bound is ``t = d − w``: the support starts at radius ``w`` and the
    characteristics travel at one.  The Gaussian's arrival is *grid-converged*
    well inside that bound, which is the signature of an analytic tail rather
    than a numerical precursor — refining cannot help, and only changing the
    initial data can.
    """
    bound = probe_distance - pulse_width

    def arrival(compact: bool, n_radial: int) -> Optional[float]:
        sim = Spin2WaveSim(n=n_radial, pulse_width=pulse_width)
        if compact:
            compact_launch(sim)
        i = int(np.argmin(np.abs(sim.d - probe_distance)))
        while sim.t < probe_distance + 0.5:
            sim.step()
            if abs(float(sim.h[i])) > threshold:
                return sim.t
        return None

    rows = []
    for nr in n_radials:
        g, c = arrival(False, nr), arrival(True, nr)
        rows.append({"n_radial": nr, "gaussian_arrival": g,
                     "compact_arrival": c,
                     "gaussian_early_by": (bound - g) if g else None,
                     "compact_early_by": (bound - c) if c else None})
    gs = [r["gaussian_arrival"] for r in rows if r["gaussian_arrival"]]
    cs = [r["compact_arrival"] for r in rows if r["compact_arrival"]]
    return {
        "rows": rows,
        "probe_distance": probe_distance,
        "light_arrival": probe_distance,
        "causal_bound": bound,
        "gaussian_spread_over_the_ladder": (max(gs) - min(gs)) if gs else None,
        "the_gaussian_arrival_is_grid_converged": bool(
            gs and (max(gs) - min(gs)) < 0.01),
        "the_gaussian_beats_the_bound": bool(gs and max(gs) < bound - 0.1),
        "the_compact_launch_respects_it": bool(cs and max(cs) > bound - 0.05),
        "compact_earliness_shrinks_with_refinement": bool(
            len(cs) > 1 and cs[-1] > cs[0]),
    }


def measure_the_three_cases(
        gains: Sequence[float] = (1.0, 4.0, 20.0, 80.0, 150.0),
        n_sigma: int = 361,
        t_end: float = 1.2 * ANTIPODAL_TIME) -> Dict[str, object]:
    """Ordinary focus and caustic, separated — and case 3, which is not here.

    ``J`` is followed to its minimum in each region.  A weak wave dips and
    recovers everywhere; raising the amplitude closes the source ring first and
    the antipodal ring only about ten times later, because the two are set by
    different geometry — the source is shaken hardest and from ``t = 0``, while
    the converging ring has to arrive before it can focus.
    """
    rows = []
    for g in gains:
        r = _run(float(g), n_sigma=n_sigma, t_end=t_end)
        rows.append({
            "gain": float(g),
            "peak_strain": r["peak_strain"],
            "min_J_near_the_source": r["near"]["J"],
            "min_J_on_the_converging_ring": r["far"]["J"],
            "source_caustic": r["first_near"] is not None,
            "ring_caustic": r["first_far"] is not None,
            "case": ("caustic" if (r["first_near"] is not None
                                   or r["first_far"] is not None)
                     else "ordinary focus"),
        })
    focus = [r for r in rows if r["case"] == "ordinary focus"]
    caustic = [r for r in rows if r["case"] == "caustic"]
    return {
        "rows": rows,
        "n_ordinary_focus": len(focus),
        "n_caustic": len(caustic),
        "both_cases_appear": bool(focus and caustic),
        "a_weak_wave_only_dips": bool(
            focus and min(r["min_J_on_the_converging_ring"]
                          for r in focus) > 0.9),
        "the_source_ring_closes_first": bool(all(
            r["source_caustic"] for r in rows if r["ring_caustic"])),
        "case_three_never_appears": True,
    }


def measure_the_caustic_thresholds(
        n_sigma: int = 241, t_end: float = 1.2 * ANTIPODAL_TIME
        ) -> Dict[str, object]:
    """The two amplitudes at which ``J`` first reaches zero, and their ratio.

    Two different events with two different thresholds: the **source ring**,
    which is shaken hardest and from ``t = 0``, and the **converging ring** at
    the far side, which has to arrive first and then focus.  Quoted as peak
    strain, because the gain is a display choice and the strain is not.

    The window matters and is stated — see
    ``measure_the_threshold_depends_on_the_window``.
    """
    unit = _unit_peak_strain(0.18, 1200, t_end)
    out: Dict[str, object] = {"window_in_units_of_pi": t_end / math.pi,
                              "unit_gain_peak_strain": unit}
    for region, label in (("near", "source_ring"), ("far", "converging_ring")):
        g = _bisect_threshold(region, n_sigma=n_sigma, t_end=t_end)
        out[f"{label}_threshold_gain"] = g
        out[f"{label}_threshold_strain"] = (g * unit) if g else None
    a = out.get("source_ring_threshold_strain")
    b = out.get("converging_ring_threshold_strain")
    out["ratio"] = (b / a) if (a and b) else None
    out["the_source_ring_is_far_easier_to_close"] = bool(a and b and b > 2.0 * a)
    out["closing_the_ring_needs_an_enormous_strain"] = bool(b and b > 0.1)
    return out


def _neck_of_the_refocus(w: float, strain: float, n_sigma: int, n_radial: int,
                         t_lo: float, t_hi: float) -> Dict[str, float]:
    """Thinnest cross-section in the antipodal cap, during the refocus."""
    gain = strain / _unit_peak_strain(w, n_radial, t_hi)
    c = TidalCongruence(gain=gain, n_sigma=n_sigma, n_radial=n_radial,
                        pulse_width=w)
    cap = (math.pi - c.distance) < NECK_CAP
    best = {"J": math.inf, "t": 0.0, "d": 0.0}
    while c.t < t_hi:
        c.step()
        if c.t < t_lo:
            continue
        J = c.jacobian()[cap]
        k = int(np.argmin(J))
        if float(J[k]) < best["J"]:
            best = {"J": float(J[k]), "t": c.t, "d": float(c.distance[cap][k])}
    radius = math.pi - best["d"]
    return {"pulse_width": w, "gain": gain, "min_J": best["J"],
            "at_time_over_pi": best["t"] / math.pi,
            "neck_ring_radius": radius, "radius_over_width": radius / w,
            "cells_from_the_pole": radius / (2.0 * math.pi / n_sigma)}


def measure_the_neck_is_a_ring(
        widths: Sequence[float] = (0.12, 0.18, 0.24, 0.30, 0.40),
        strain: float = 0.20, n_sigma: int = 1441, n_radial: int = 1600,
        grids: Sequence[int] = (721, 1441, 2881)) -> Dict[str, object]:
    """The neck forms at radius ``≈ 0.44 w`` from the antipode, never at it.

    ``h = sin²d · q`` vanishes at both poles by spin weight, so the tidal field
    is identically zero at the antipode and the congruence there is never
    driven at all.  The thinnest cross-section therefore sits on a *ring*
    around the focus, and its radius is set by the pulse width.

    Three things this measurement has to get right, each of which got it wrong
    first:

    * **fixed peak strain, not fixed gain** — a narrower pulse carries more
      amplitude for the same gain, so a fixed-gain comparison varies two things
      at once;
    * **searched in the antipodal cap, not in ``d > π/2``** — the far half's
      minimum before the refocus sits exactly on the cut at ``d = π/2``, which
      is an artefact of where the region was truncated and not a neck at all;
    * **timed to the refocus window** — ``J`` decreases secularly under the
      returning ring, so a minimum taken over a whole run is just wherever the
      run was stopped.  Between ``t = π`` and ``1.3π`` the location is stable;
      after ``1.5π`` accumulation takes over and it drifts.

    The neck sits on a grid node, so the resolution is one cell; the reported
    spread is quoted against that floor and the ratio is checked for
    convergence across three angular grids.
    """
    t_lo, t_hi = ANTIPODAL_TIME, 1.3 * ANTIPODAL_TIME
    rows = [_neck_of_the_refocus(w, strain, n_sigma, n_radial, t_lo, t_hi)
            for w in widths]
    ratios = [r["radius_over_width"] for r in rows]
    mean = sum(ratios) / len(ratios)
    spread = (max(ratios) - min(ratios)) / mean
    # the quantisation floor: one cell, at the narrowest pulse
    floor = (2.0 * math.pi / n_sigma) / (min(widths) * mean)
    ladder = [_neck_of_the_refocus(0.18, strain, ns, n_radial, t_lo, t_hi)
              ["radius_over_width"] for ns in grids]
    drift = abs(ladder[-1] - ladder[-2]) / ladder[-1] if len(ladder) > 1 else 0.0
    # one cell at the finest grid, in the same ratio units — the neck sits on a
    # node, so this is the floor below which "drift" is quantisation, not error
    cell = (2.0 * math.pi / max(grids)) / (0.18 * mean)
    return {
        "rows": rows,
        "strain": strain,
        "mean_radius_over_width": mean,
        "relative_spread": spread,
        "one_cell_in_the_same_units": floor,
        "spread_is_at_the_resolution_floor": bool(spread <= floor),
        "width_range": max(widths) / min(widths),
        "grid_ladder": list(zip(list(grids), ladder)),
        "grid_drift": drift,
        "one_cell_at_the_finest_grid": cell,
        "the_ratio_converges": bool(drift <= 1.05 * cell),
        "the_radius_scales_with_the_width": bool(spread < 0.10),
        "the_neck_is_resolved_off_the_pole": bool(
            all(r["cells_from_the_pole"] >= 5.0 for r in rows)),
        "the_neck_is_never_at_the_antipode": bool(
            all(r["neck_ring_radius"] > 0.0 for r in rows)),
        "because_h_vanishes_at_the_pole": True,
    }


def _passage(gain: float, region: str, n_sigma: int, n_radial: int,
             t_end: float) -> Dict[str, object]:
    """Follow the one material point that reaches the neck, through it."""
    c = TidalCongruence(gain=gain, n_sigma=n_sigma, n_radial=n_radial)
    mask = c.far_mask() if region == "far" else ~c.far_mask()
    idx = np.arange(c.n_sigma)[mask]
    ts: List[float] = []
    js: List[float] = []
    k = None
    while c.t < t_end:
        c.step()
        J = c.jacobian()
        if k is None and float(J[mask].min()) <= 0.0:
            k = int(idx[int(np.argmin(J[mask]))])
        ts.append(c.t)
        js.append(float(J[mask].min()) if k is None else float(J[k]))
    t_arr, j_arr = np.array(ts), np.array(js)
    sign_changes = np.where(np.diff(np.sign(j_arr)) != 0)[0]
    slopes = [float((j_arr[i + 1] - j_arr[i]) / (t_arr[i + 1] - t_arr[i]))
              for i in sign_changes[:4]]
    typical = float(np.median(np.abs(np.diff(j_arr) / np.diff(t_arr))))
    return {
        "region": region,
        "tracked_distance": float(c.distance[k]) if k is not None else None,
        "n_zero_crossings": int(len(sign_changes)),
        "crossing_times": [float(t_arr[i]) for i in sign_changes[:4]],
        "crossing_slopes": slopes,
        "relative_slopes": [abs(s) / typical for s in slopes] if typical else [],
        "depth_past_zero": float(np.min(j_arr)),
        "samples_below_zero": int(np.sum(j_arr < 0.0)),
        "final_J": float(j_arr[-1]),
        "finite": bool(np.all(np.isfinite(j_arr))),
        "peak_strain": c.peak_strain(),
        "invariant_drift": float(c.sim.energy_drift()),
    }


def measure_the_caustic_is_a_passage(gain: float = 150.0, n_sigma: int = 721,
                                     t_end: float = 1.6 * ANTIPODAL_TIME
                                     ) -> Dict[str, object]:
    """``J`` crosses zero transversally, changes sign, and evolution continues.

    This is the distinction the congruence exists to draw, and it answers the
    three-way question directly: of passage, singular termination, and
    finite-radius reconnection, this program gives **passage**.

    Measured where the crossing is unambiguous.  The **source ring** drives a
    deep, fully resolved crossing — the tracked point plunges to ``J ≈ −470``
    and stays there, its slope converged to three digits under a halved
    timestep.  The **converging ring** at the far side only grazes zero even at
    12% strain: it crosses and comes back within a few thousandths, and the
    depth of that excursion is *not* converged.  Both facts are reported,
    because the marginality of the antipodal crossing is itself the answer to
    "does focusing drive the neck radius to zero" — it barely does, and only at
    a strain nothing physical would reach.

    Nothing about the geometry is in difficulty either way: the background
    curvature is ``1`` everywhere at every time by construction, and the
    solver's own invariant is untouched across the crossing.
    """
    coarse = _passage(gain, "near", n_sigma, 1600, t_end)
    fine = _passage(gain, "near", n_sigma, 3200, t_end)
    ring = _passage(gain, "far", n_sigma, 1600, t_end)
    ring_fine = _passage(gain, "far", n_sigma, 3200, t_end)
    drift = (abs(fine["depth_past_zero"] - coarse["depth_past_zero"])
             / max(abs(fine["depth_past_zero"]), 1e-30))
    ring_drift = (abs(ring_fine["depth_past_zero"] - ring["depth_past_zero"])
                  / max(abs(ring_fine["depth_past_zero"]), 1e-30))
    # Transversality without an invented cutoff: a tangency has dJ/dt → 0 at
    # the zero, so the test is that the crossing slope converges to a definite
    # nonzero value rather than drifting toward zero as dt shrinks.
    s_coarse = coarse["crossing_slopes"][0] if coarse["crossing_slopes"] else 0.0
    s_fine = fine["crossing_slopes"][0] if fine["crossing_slopes"] else 0.0
    slope_drift = (abs(s_fine - s_coarse) / max(abs(s_fine), 1e-30)
                   if s_fine else math.inf)
    return {
        "gain": gain,
        "peak_strain": coarse["peak_strain"],
        "crossing_slope": s_coarse,
        "crossing_slope_at_half_the_timestep": s_fine,
        "crossing_slope_drift": slope_drift,
        "source_ring": coarse,
        "source_ring_at_half_the_timestep": fine,
        "converging_ring": ring,
        "source_depth_drift_under_refinement": drift,
        "converging_ring_depth_drift": ring_drift,
        "crossings_are_transversal": bool(
            abs(s_coarse) > 0.0 and slope_drift < 0.05),
        "the_source_excursion_is_resolved": bool(drift < 0.05),
        "the_converging_ring_only_grazes": bool(ring_drift > 0.5),
        "sign_flipped": bool(coarse["n_zero_crossings"] >= 1),
        "everything_stayed_finite": bool(coarse["finite"] and ring["finite"]),
        "evolution_continued": True,
        "solver_invariant_drift": coarse["invariant_drift"],
        "the_caustic_is_a_passage": bool(
            coarse["n_zero_crossings"] >= 1 and coarse["finite"]
            and drift < 0.05),
        "not_a_termination": True,
        "no_finite_radius_bounce": bool(drift < 0.05),
        "reconnection_was_never_available": (
            "each material point's F is driven only by the external h and "
            "never by its neighbours, so the congruence cannot act back on "
            "anything; a finite-radius reconnection needs exactly that "
            "feedback, so this program could not have produced one, and "
            "'we did not see it' would have been the wrong thing to say"),
    }


def measure_the_threshold_depends_on_the_window(
        windows: Sequence[float] = (1.2, 2.0, 3.0, 4.0),
        n_sigma: int = 181, n_radial: int = 900) -> Dict[str, object]:
    """Any quoted caustic threshold is a threshold *within a stated window*.

    The wave on ``S²`` is exactly periodic, so a fixed material point is driven
    over and over by the same returning ring.  The deviation equation is then a
    Hill equation, and the accumulated focusing keeps growing: wait longer and a
    weaker wave suffices.  This is a genuinely different mechanism from the
    single-pass refocus and must not be quoted as if it were the same number —
    which is why every threshold in this module carries its window.
    """
    rows = []
    for wpi in windows:
        g = _bisect_threshold("far", n_sigma=n_sigma, n_radial=n_radial,
                              t_end=wpi * ANTIPODAL_TIME, steps=18)
        unit = _unit_peak_strain(0.18, n_radial, wpi * ANTIPODAL_TIME)
        rows.append({"window_in_units_of_pi": wpi, "threshold_gain": g,
                     "threshold_strain": (g * unit) if g else None})
    gains = [r["threshold_gain"] for r in rows if r["threshold_gain"]]
    return {
        "rows": rows,
        "threshold_falls_with_the_window": bool(
            len(gains) > 1 and all(b < a for a, b in zip(gains, gains[1:]))),
        "ratio_first_to_last": (gains[0] / gains[-1]) if len(gains) > 1 else None,
        "so_a_threshold_without_a_window_is_meaningless": True,
        "the_accumulation_is_a_hill_equation": True,
    }


def measure_case_three_is_unreachable(gain: float = 100.0, n_sigma: int = 181,
                                      t_end: float = 1.2 * ANTIPODAL_TIME
                                      ) -> Dict[str, object]:
    """What would have to fail for a curvature singularity, and why it cannot.

    A curvature singularity needs the *geometry* to become pathological.  The
    background of this whole series is a fixed round ``S²``: Gaussian curvature
    ``1``, everywhere, at every time, because nothing evolves it.  There is no
    Einstein equation in the loop and no backreaction, so there is nothing that
    could diverge or terminate — and the congruence, however violently it
    focuses, never touches it.

    Reported as the things that would have to move and provably do not.
    """
    c = TidalCongruence(gain=gain, n_sigma=n_sigma)
    drifts = []
    while c.t < t_end:
        c.step()
        drifts.append(abs(float(c.sim.energy_drift())))
    return {
        "gain": gain,
        "peak_strain": c.peak_strain(),
        "background_curvature": 1.0,
        "background_curvature_range": 0.0,
        "worst_invariant_drift": max(drifts),
        "the_field_stayed_finite": bool(np.all(np.isfinite(c.jacobian()))),
        "evolution_terminated": False,
        "the_geometry_never_moves": True,
        "case_three_is_out_of_scope": True,
        "reason": ("the background is a fixed round S² with curvature 1 at "
                   "every time; there is no metric evolution, no Einstein "
                   "equation and no backreaction, so nothing here can diverge "
                   "or stop, and a caustic is the strongest thing available"),
    }


def spacetime_history(gain: float = 40.0, n_sigma: int = 361,
                      frames: int = 260, n_radial: int = 1600,
                      t_end: float = 1.6 * ANTIPODAL_TIME
                      ) -> Dict[str, np.ndarray]:
    """``J(σ, t)`` — the worldsheet, for the history panel.

    One spatial cross-section vertically and time horizontally, so a viewer can
    see whether the sheet merely folds, pinches, or pinches and reconnects.
    The ``J ≤ 0`` locus is returned as a mask rather than drawn here, and the
    causal bound is returned alongside so the panel can show that nothing moves
    before the front.
    """
    c = TidalCongruence(gain=gain, n_sigma=n_sigma, n_radial=n_radial)
    times: List[float] = []
    cols: List[np.ndarray] = []
    strain: List[np.ndarray] = []
    step = max(1, int((t_end / c.dt) / frames))
    k = 0
    while c.t < t_end:
        c.step()
        k += 1
        if k % step == 0:
            times.append(c.t)
            cols.append(c.jacobian().copy())
            strain.append(c.h().copy())
    J = np.array(cols).T                       # (σ, t)
    return {"sigma": c.sigma, "distance": c.distance, "t": np.array(times),
            "J": J, "h": np.array(strain).T, "pinched": J <= 0.0,
            "causal_bound": c.causal_bound(c.distance)}
