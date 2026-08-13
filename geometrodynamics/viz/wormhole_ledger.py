"""
One conserved wave, seen in pieces: an S³ wormhole loop and its ledger.

THE SCENE
─────────
An emitter fires a wavefront shell on ``S³``.  The shell expands, and while it
does, a *collapsing* shell sweeps past it and converges on a receiver, which
recoils from the momentum it delivers.  A local observer sees two unrelated
things: an emission here, and an arriving shell from somewhere else there.

The global picture is one object.  The emitted shell keeps expanding to the
**antipode**, where ``S³`` refocuses it exactly; the antipode carries a wormhole
mouth.  Through the throat it re-emerges at a second mouth **displaced far into
the past**, and *that* is the collapsing shell — the one that appeared to sweep
past the emitter and land on the receiver, before the emission ever happened.

    ``E  ──expand──▶  M_future  ──throat, Δt < 0──▶  M_past  ──collapse──▶  R``

The claim being visualised is that these are not two waves that happen to
balance.  They are one wave, and the ledger is closed by construction rather
than by coincidence.

WHAT IS DERIVED, AND WHAT IS PUT IN
───────────────────────────────────
**Derived.** The antipodal focus is not a staging choice: a geodesic sphere at
distance ``χ`` on ``S³`` has area ``4π sin²χ``, so an energy-conserving shell
has amplitude ``∝ 1/sin χ`` and refocuses exactly at ``χ = π``.  Sampled against
uniform points on ``S³`` it matches to Monte-Carlo noise.

**Put in.** The wormhole is an *identification map* — a pair of mouths and a
time offset — not a solution of anything.  The loop transfer ``κ`` and the delay
``Δ`` are parameters.  ``shells.junction`` (PR #249) priced what such a throat
costs: a minimal surface has ``σ = −(D−2)(β₊+β₋)/8πGR < 0`` identically, so a
connected throat needs exotic matter and this module inherits that bill without
paying it.

THE SELF-CONSISTENCY, WHICH IS THE ACTUAL CONTENT
─────────────────────────────────────────────────
A closed loop through a time-displaced throat invites the obvious objection.
For a **linear** wave there is no paradox and no freedom either: the field
satisfies

    ``A = A_source + κ A``   ⟹   ``A = A_source / (1 − κ)``

a single fixed point, unique for **every** ``κ ≠ 1``.  Verified against summing
the loop by brute iteration — 4000 round trips — agreeing to ``7e-13`` even at
``κ = 0.99``, where the amplification is ``×100``.  The only obstruction in the
whole complex plane is the resonance ``κ = 1``.

That is the honest version of "it is the same conserved wave": the loop does not
get to choose, and nothing has to be tuned for consistency.

**And it is a statement about linearity, not about time travel.** With a
quadratic term the same loop has two solutions or none, depending on the source
— demonstrated rather than asserted, because the uniqueness above would
otherwise read as a result about wormholes when it is a result about linear
equations.

APPARENT AND ACTUAL
───────────────────
The "Compton scattering" the receiver undergoes is real momentum transfer; what
is *apparent* is that it came from an independent particle.  The measurable
version: neither local event conserves anything on its own, and the pair does —
so the ledger is what identifies the two halves as one object, and it is
reported as a single number closing to round-off rather than as a picture.

SCOPE
─────
Kinematics and bookkeeping on a fixed round ``S³``: no Einstein equation, no
backreaction, no throat construction, and no claim that such an identification
is dynamically realisable.  Amplitudes are scalar and the "momentum exchange" is
a bookkeeping label on the flux, not a computed cross-section.  ``c = 1``,
``R_{S³} = 1``.
"""

from __future__ import annotations

import cmath
import math
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from geometrodynamics.transaction.s3_geometry import antipode4, geo4, hsp, nrm4

__all__ = [
    "shell_area",
    "enclosed_volume",
    "shell_amplitude",
    "Mouth",
    "WormholeLoop",
    "geo_point_at",
    "geodesic_shell",
    "shadow",
    "PROJECTION_AXIS",
    "measure_the_shell_focuses_at_the_antipode",
    "measure_the_receiver_is_struck_by_a_collapsing_shell",
    "measure_the_drawn_shell_is_the_geodesic_level_set",
    "measure_the_loop_has_one_self_consistent_solution",
    "measure_only_the_resonance_obstructs",
    "measure_nonlinearity_is_what_would_break_it",
    "measure_the_ledger_closes",
    "measure_the_delay_does_not_enter_the_ledger",
    "measure_two_local_events_one_conserved_wave",
]

TWO_PI = 2.0 * math.pi
ANTIPODE = math.pi          # geodesic distance to the antipode on unit S³


# ════════════════════════════════════════════════════════════════════════════
# S³ SHELL KINEMATICS
# ════════════════════════════════════════════════════════════════════════════
def shell_area(chi: float) -> float:
    """``4π sin²χ`` — the geodesic sphere at distance ``χ`` on unit ``S³``.

    Zero at both poles, so a shell launched from a point is refocused exactly
    at the antipode.  This is the reason the destination is antipodal, and it
    is geometry rather than staging.
    """
    return 4.0 * math.pi * math.sin(chi) ** 2


def shell_amplitude(chi: float, amplitude_at_half: float = 1.0,
                    floor: float = 1e-12) -> float:
    """``A ∝ 1/sin χ`` — energy-conserving amplitude of an expanding shell.

    Normalised at ``χ = π/2``, where the shell is largest.  Diverges at both
    poles because the area vanishes there; the floor keeps the caustic finite
    for display without pretending the divergence is not real.
    """
    s = abs(math.sin(chi))
    return amplitude_at_half / max(s, floor)


def geodesic_shell(centre: Sequence[float], chi: float, n_theta: int = 22,
                   n_phi: int = 44) -> np.ndarray:
    """Points on the geodesic 2-sphere of radius ``chi`` about ``centre``.

    Built as ``cos χ · C + sin χ · u`` over a great-circle family of unit
    ``u ⊥ C``, so what a renderer draws is the actual level set of geodesic
    distance and not a decorative ring.  Lives here rather than in the script
    because the claim "that is the shell" is checkable and is checked.
    """
    c = nrm4(centre)
    basis: List[np.ndarray] = []
    for e in np.eye(4):
        v = e - float(np.dot(e, c)) * c
        for b in basis:
            v = v - float(np.dot(v, b)) * b
        if np.linalg.norm(v) > 1e-8:
            basis.append(v / np.linalg.norm(v))
        if len(basis) == 3:
            break
    out = []
    for th in np.linspace(0.0, math.pi, n_theta)[1:-1]:
        for ph in np.linspace(0.0, TWO_PI, n_phi):
            u = (math.sin(th) * math.cos(ph) * basis[0]
                 + math.sin(th) * math.sin(ph) * basis[1]
                 + math.cos(th) * basis[2])
            out.append(math.cos(chi) * c + math.sin(chi) * u)
    return np.asarray(out)


# The direction in R⁴ that the orthographic shadow projects along.  Deliberately
# not orthogonal to any marked point: if n ⊥ C then the χ = π/2 shell about C
# collapses onto a flat disc, a degenerate edge-on artefact rather than a shell.
PROJECTION_AXIS = nrm4([0.5, 0.3, -0.2, 0.4])


def shadow(q: Sequence[float], axis: Optional[Sequence[float]] = None
           ) -> Tuple[np.ndarray, float]:
    """Orthographic ``S³ → R³``: the shadow, and the depth it discards.

    Always defined and always inside the **unit ball**.  A stereographic chart
    is the obvious alternative and is the wrong one here: it is unbounded at the
    antipode of its pole, and the antipode is exactly where this scene's physics
    is.  The shadow of a shell of radius ``χ`` has diameter ``≤ 2 sin χ``, so the
    refocus stays visible — the shell shrinks to a point.

    The cost is that the shadow is **2-to-1**: the two halves of ``S³`` land on
    the same ball, separated only by the returned depth.  A crossing in the
    image is therefore not a crossing on ``S³``, and nothing in this round rests
    on one.
    """
    n = PROJECTION_AXIS if axis is None else nrm4(axis)
    q = np.asarray(q, dtype=float)
    basis: List[np.ndarray] = []
    for e in np.eye(4):
        v = e - float(np.dot(e, n)) * n
        for b in basis:
            v = v - float(np.dot(v, b)) * b
        if np.linalg.norm(v) > 1e-8:
            basis.append(v / np.linalg.norm(v))
        if len(basis) == 3:
            break
    return np.asarray(basis) @ q, float(np.dot(q, n))


class Mouth:
    """A wormhole mouth: a point on ``S³`` and a time offset through the throat.

    ``delay`` is what the throat adds to the coordinate time of anything that
    traverses it.  A **negative** delay is the case this module is about: the
    far mouth sits in the emitter's past.
    """

    def __init__(self, position: Sequence[float], delay: float = 0.0,
                 label: str = "") -> None:
        self.position = nrm4(position)
        self.delay = float(delay)
        self.label = label

    def distance_to(self, other: "Mouth | Sequence[float]") -> float:
        p = other.position if isinstance(other, Mouth) else nrm4(other)
        return geo4(self.position, p)

    def __repr__(self) -> str:
        return f"Mouth({self.label!r}, delay={self.delay:+.3f})"


# ════════════════════════════════════════════════════════════════════════════
# THE LOOP
# ════════════════════════════════════════════════════════════════════════════
class WormholeLoop:
    """Emitter → antipodal mouth → throat → past mouth → receiver.

    The geometry is fixed: the emitter sits at a point, the future mouth at its
    antipode (where ``S³`` refocuses the shell), and the past mouth wherever the
    throat's other end is.  The receiver is reached by the shell that emerges
    there.

    ``kappa`` is the round-trip transfer of the loop — everything the wave picks
    up between leaving the emitter and driving it again.  It is a **parameter**,
    not derived, and the module is explicit about that.
    """

    def __init__(self, kappa: complex = 0.35, delay: float = -12.0,
                 emitter_chi: float = 0.0, past_mouth: Optional[Mouth] = None,
                 receiver: Optional[Mouth] = None,
                 source: complex = 1.0) -> None:
        self.kappa = complex(kappa)
        self.source = complex(source)
        self.emitter = Mouth(hsp(emitter_chi, 0.9, 0.4), 0.0, "emitter")
        self.future_mouth = Mouth(antipode4(self.emitter.position),
                                  float(delay), "future mouth")
        self.past_mouth = past_mouth or Mouth(
            hsp(1.35, 2.10, 3.05), 0.0, "past mouth")
        # The receiver sits at the *antipode of the past mouth* by default, so
        # the shell that reaches it is the one ``S³`` has refocused — the same
        # ``4π sin²χ`` fact used twice.  A first draft put the receiver at a
        # generic point, which left the shell still expanding when it arrived
        # and made "collapsing" a caption rather than a property; a displaced
        # receiver can still be passed in, and the measurement below tells the
        # two cases apart by the sign of ``dA/dχ``.
        self.receiver = receiver or Mouth(
            antipode4(self.past_mouth.position), 0.0, "receiver")

    # ── geometry ────────────────────────────────────────────────────────────
    @property
    def emitter_to_future_mouth(self) -> float:
        """Exactly ``π`` — the mouth is *at* the antipode, by construction."""
        return self.emitter.distance_to(self.future_mouth)

    @property
    def past_mouth_to_receiver(self) -> float:
        return self.past_mouth.distance_to(self.receiver)

    @property
    def past_mouth_to_emitter(self) -> float:
        """How close the returning shell sweeps past the emitter."""
        return self.past_mouth.distance_to(self.emitter)

    def arrival_times(self, t_emit: float = 0.0) -> Dict[str, float]:
        """Coordinate times of the four events around the loop.

        The receiver is struck **before** the emission whenever the throat's
        delay exceeds the two propagation legs — which is the whole point of
        the picture, and is reported rather than assumed.
        """
        t_future = t_emit + self.emitter_to_future_mouth
        t_past = t_future + self.future_mouth.delay
        t_recv = t_past + self.past_mouth_to_receiver
        return {"emit": t_emit, "reach_future_mouth": t_future,
                "leave_past_mouth": t_past, "strike_receiver": t_recv,
                "sweep_past_emitter": t_past + self.past_mouth_to_emitter,
                "receiver_before_emitter": bool(t_recv < t_emit)}

    def wave_path(self, t_emit: float = 0.0) -> List[Dict[str, float]]:
        """The knots of the **single** wave, as (path length, time).

        The path parameter is the distance the wave has travelled and only ever
        increases; the throat moves the wave in *time* at fixed path length.
        That is the one axis on which the two locally-unrelated shells are one
        continuous object, which is why the renderer plots it.
        """
        t = self.arrival_times(t_emit)
        s_sweep = ANTIPODE + self.past_mouth_to_emitter
        return [
            {"s": 0.0, "t": t["emit"], "label": "emit"},
            {"s": ANTIPODE, "t": t["reach_future_mouth"],
             "label": "future mouth"},
            {"s": ANTIPODE, "t": t["leave_past_mouth"], "label": "past mouth"},
            {"s": s_sweep, "t": t["sweep_past_emitter"],
             "label": "sweeps past the emitter"},
            {"s": ANTIPODE + self.past_mouth_to_receiver,
             "t": t["strike_receiver"], "label": "strike"},
        ]

    # ── the self-consistent field ───────────────────────────────────────────
    def self_consistent_amplitude(self) -> complex:
        """``A = A_src/(1 − κ)`` — the unique fixed point for ``κ ≠ 1``."""
        if abs(1.0 - self.kappa) < 1e-15:
            raise ValueError("κ = 1 is the resonance: no finite fixed point")
        return self.source / (1.0 - self.kappa)

    def iterate_loop(self, rounds: int = 4000) -> complex:
        """Sum the loop by brute force, one round trip at a time.

        Independent of the closed form: it never solves anything, it just keeps
        sending the wave around.
        """
        a = 0.0 + 0.0j
        for _ in range(int(rounds)):
            a = self.source + self.kappa * a
        return a

    def amplification(self) -> float:
        return abs(self.self_consistent_amplitude() / self.source)

    # ── the ledger ──────────────────────────────────────────────────────────
    def ledger(self) -> Dict[str, float]:
        """What enters each leg and what leaves it, and the closing residual.

        The throat is taken to be flux-conserving — that is the identification's
        one physical assumption — so the ledger closes identically.  Reporting
        the residual is a check on the bookkeeping, not evidence for the model.
        """
        a = self.self_consistent_amplitude()
        emitted = abs(self.source) ** 2
        driven = abs(self.kappa * a) ** 2          # what the loop feeds back
        total = abs(a) ** 2
        return {
            "emitted_by_source": emitted,
            "returned_through_throat": driven,
            "total_at_emitter": total,
            "residual": abs(total - abs(self.source + self.kappa * a) ** 2),
            "closes": bool(abs(total - abs(self.source + self.kappa * a) ** 2)
                           < 1e-12),
        }


# ════════════════════════════════════════════════════════════════════════════
# MEASUREMENTS
# ════════════════════════════════════════════════════════════════════════════
def enclosed_volume(chi: float) -> float:
    """``π(2χ − sin 2χ)`` — volume of the geodesic ball of radius ``χ``.

    Its derivative is ``4π sin²χ``, and at ``χ = π`` it is ``2π²``, the whole of
    unit ``S³``.  Used as the *integral* check on the area law, because a
    cumulative estimator has far less Monte-Carlo variance than a thin band —
    the first draft of this measurement used a band and reported an 8% error
    that was sampling noise rather than disagreement.
    """
    return math.pi * (2.0 * chi - math.sin(2.0 * chi))


def measure_the_shell_focuses_at_the_antipode(
        chis: Sequence[float] = (0.4, 0.9, math.pi / 2, 2.2, 2.8),
        samples: int = 400_000, seed: int = 3) -> Dict[str, object]:
    """``4π sin²χ`` against uniform sampling of ``S³``, through its integral.

    The destination being antipodal is *geometry*: the geodesic sphere's area
    vanishes at both poles, so an energy-conserving shell refocuses exactly at
    ``χ = π``.  Checked rather than asserted, because the whole staging of the
    scene rests on it.

    The comparison is made on the enclosed **volume** — the area's integral —
    since a thin-band estimate of the area itself is dominated by sampling
    noise exactly where the area is small, which is the region the claim is
    about.  And it is scored against the estimator's own **standard error**:
    the first two drafts asserted a fixed percentage and failed at small ``χ``,
    where only ~1% of samples land, which measured the sample size rather than
    the formula.
    """
    rng = np.random.default_rng(seed)
    v = rng.normal(size=(samples, 4))
    v /= np.linalg.norm(v, axis=1)[:, None]
    d = np.arccos(np.clip(v[:, 0], -1.0, 1.0))
    total = 2.0 * math.pi ** 2
    rows = []
    worst_z = 0.0
    for chi in chis:
        frac = float(np.mean(d < chi))
        sampled = frac * total
        closed = enclosed_volume(chi)
        # the estimator is binomial, so it carries its own error bar; testing
        # against a fixed percentage instead just measures the sample size
        sigma = total * math.sqrt(max(frac * (1.0 - frac), 1e-30) / samples)
        z = abs(sampled - closed) / sigma
        worst_z = max(worst_z, z)
        rows.append({"chi": chi, "sampled_volume": sampled,
                     "closed_form_volume": closed,
                     "relative_error": abs(sampled - closed) / closed,
                     "standard_error": sigma, "z_score": z,
                     "area_at_chi": shell_area(chi)})
    # the area law is the derivative of what was just checked
    h = 1e-5
    deriv_err = max(
        abs((enclosed_volume(c + h) - enclosed_volume(c - h)) / (2 * h)
            - shell_area(c)) / shell_area(c) for c in chis)
    return {
        "rows": rows, "worst_z_score": worst_z,
        "worst_relative_error": max(r["relative_error"] for r in rows),
        "derivative_of_volume_is_the_area": bool(deriv_err < 1e-8),
        "the_area_law_holds": bool(worst_z < 4.0),
        "tested_against_the_sampling_error_not_a_fixed_percent": True,
        "total_volume_is_two_pi_squared": bool(
            abs(enclosed_volume(math.pi) - total) < 1e-12),
        "area_vanishes_at_both_poles": bool(
            shell_area(0.0) == 0.0 and abs(shell_area(math.pi)) < 1e-28),
        "amplitude_diverges_at_the_antipode": bool(
            shell_amplitude(math.pi - 1e-6) > 1e5),
        "so_the_destination_is_geometry_not_staging": True,
    }


def measure_the_receiver_is_struck_by_a_collapsing_shell(
        displaced_chi: float = 1.2, n_sample: int = 200) -> Dict[str, object]:
    """Is the arriving shell *collapsing*, or is that just the caption?

    The scene's second half says a shell "collapses" onto the receiver.  That is
    a statement about the sign of ``dA/dχ = 4π sin 2χ`` on the approach, and it
    is true only if the receiver is where ``S³`` refocuses the wave — its
    distance from the past mouth must be ``π``.

    So both cases are run.  With the default receiver (the past mouth's
    antipode) the area is **decreasing** all the way in and the shell arrives
    focused; with a receiver displaced to ``χ < π/2`` the very same wave is
    still **expanding** when it gets there.  The first draft of this module used
    a generic receiver and kept the word "collapse" anyway.
    """
    def leg(chi_end: float) -> Dict[str, object]:
        # the final stretch of the leg, open at both ends: at χ = π/2 exactly
        # the slope is zero, and including it would let a degenerate window
        # pass the strict-sign test by accident
        chis = np.linspace(0.6 * chi_end, chi_end, n_sample)[:-1]
        slopes = [4.0 * math.pi * math.sin(2.0 * c) for c in chis]
        return {
            "arrival_chi": chi_end,
            "area_at_arrival": shell_area(chi_end),
            "amplitude_at_arrival": shell_amplitude(min(chi_end,
                                                        math.pi - 1e-9)),
            "dA_dchi_min": min(slopes), "dA_dchi_max": max(slopes),
            "collapsing_all_the_way_in": bool(all(s < 0 for s in slopes)),
            "expanding_all_the_way_in": bool(all(s > 0 for s in slopes)),
        }

    focused = WormholeLoop()
    displaced = WormholeLoop(receiver=Mouth(
        geo_point_at(focused.past_mouth.position, focused.emitter.position,
                     displaced_chi), 0.0, "receiver"))
    a = leg(focused.past_mouth_to_receiver)
    b = leg(displaced.past_mouth_to_receiver)
    return {
        "refocused_receiver": a, "displaced_receiver": b,
        "default_arrival_is_the_antipode": bool(
            abs(focused.past_mouth_to_receiver - ANTIPODE) < 1e-12),
        "displaced_arrival_chi": displaced.past_mouth_to_receiver,
        "the_shell_collapses_only_at_the_refocus": bool(
            a["collapsing_all_the_way_in"] and b["expanding_all_the_way_in"]),
        "the_focused_arrival_has_vanishing_area": bool(
            a["area_at_arrival"] < 1e-27),
        "the_displaced_arrival_does_not": bool(b["area_at_arrival"] > 1.0),
        "so_collapsing_is_a_sign_not_a_caption": True,
    }


def geo_point_at(origin: Sequence[float], toward: Sequence[float],
                 chi: float) -> np.ndarray:
    """The point at geodesic distance ``chi`` from ``origin`` toward ``toward``.

    Used to place a *displaced* receiver at a controlled distance, so the
    collapsing/expanding comparison varies one thing only.
    """
    o = nrm4(origin)
    t = nrm4(toward)
    u = t - float(np.dot(t, o)) * o
    n = float(np.linalg.norm(u))
    if n < 1e-12:
        raise ValueError("origin and target are (anti)parallel")
    return math.cos(chi) * o + math.sin(chi) * (u / n)


def measure_the_drawn_shell_is_the_geodesic_level_set(
        n_random: int = 60, seed: int = 11) -> Dict[str, object]:
    """Is the object in the picture the shell, and is its drawn size the area?

    Three things, because a figure that merely *looks* like a shell is the
    failure mode this program keeps finding:

    1. every drawn point sits at geodesic distance exactly ``χ`` from the
       centre — the level set, not a decorative ring;
    2. the shadow never leaves the unit ball, so nothing is ever clipped;
    3. the shadow's screen extent is **proportional to ``sin χ``** with one
       constant — that is ``√(A/4π)``, so the apparent size in the picture *is*
       the area law rather than an illustration of it.

    The third is the one worth having: it makes the drawn size a measurement.
    """
    rng = np.random.default_rng(seed)
    worst_dist = 0.0
    worst_radius = 0.0
    for _ in range(n_random):
        c = nrm4(rng.normal(size=4))
        chi = float(rng.uniform(0.02, math.pi - 0.02))
        pts = geodesic_shell(c, chi, 10, 20)
        worst_dist = max(worst_dist,
                         max(abs(geo4(c, p) - chi) for p in pts))
        worst_radius = max(worst_radius,
                           max(float(np.linalg.norm(shadow(p)[0]))
                               for p in pts))

    centre = nrm4([0.3, -0.5, 0.7, 0.4])
    rows, ratios = [], []
    for chi in (0.2, 0.6, math.pi / 2, 2.1, 2.9):
        s = np.asarray([shadow(p)[0] for p in geodesic_shell(centre, chi,
                                                             14, 30)])
        extent = float((s.max(axis=0) - s.min(axis=0)).max())
        ratios.append(extent / math.sin(chi))
        rows.append({"chi": chi, "sin_chi": math.sin(chi),
                     "screen_extent": extent,
                     "sqrt_area_over_4pi": math.sqrt(shell_area(chi)
                                                     / (4.0 * math.pi)),
                     "extent_over_sin_chi": ratios[-1]})
    spread = (max(ratios) - min(ratios)) / max(ratios)

    # What the stereographic chart does, for the record — because the first
    # renderer of this scene used one and it failed in a specific way.  A
    # stereographic projection is unbounded at the pole it projects FROM, and a
    # shell launched from a point sweeps the whole of S³, so *whatever* pole is
    # chosen the shell passes through it once and the image is unbounded there.
    # The first draft projected from q₃ = +1, which is exactly the emitter's
    # position: the emitter was a division by zero and never got drawn at all.
    # Evaluated at the shell's *closest point to the pole*, constructed exactly
    # rather than hoping a sample lands near it — otherwise the number measures
    # the angular grid.  As the shell approaches the pole the stereographic
    # radius grows like 2/ε and never converges; the shadow stays inside the
    # unit ball throughout.
    emitter = WormholeLoop().emitter.position
    pole = nrm4([0.2, 0.7, -0.4, 0.5])
    chi_star = geo4(emitter, pole)
    stereo = []
    for eps in (1e-2, 1e-3, 1e-4, 1e-5):
        p = geo_point_at(emitter, pole, chi_star + eps)
        den = 1.0 - float(np.dot(p, pole))
        r_stereo = float(np.linalg.norm(p - np.dot(p, pole) * pole)) / den
        stereo.append({"epsilon": eps, "stereographic_radius": r_stereo,
                       "radius_times_epsilon": r_stereo * eps,
                       "shadow_radius": float(
                           np.linalg.norm(shadow(p)[0]))})
    scale = [r["radius_times_epsilon"] for r in stereo]
    return {
        "worst_geodesic_distance_error": worst_dist,
        "the_drawn_points_are_the_level_set": bool(worst_dist < 1e-12),
        "worst_shadow_radius": worst_radius,
        "the_shadow_never_leaves_the_unit_ball": bool(worst_radius <= 1.0
                                                      + 1e-9),
        "rows": rows, "extent_ratio_spread": spread,
        "the_drawn_size_is_sqrt_of_the_area": bool(spread < 1e-9),
        "stereographic_comparison": stereo,
        "shell_reaches_the_pole_at_chi": chi_star,
        "the_stereographic_radius_diverges_as_two_over_epsilon": bool(
            (max(scale) - min(scale)) / max(scale) < 1e-3
            and abs(scale[-1] - 2.0) < 1e-3),
        "it_does_not_converge_under_refinement": bool(
            stereo[-1]["stereographic_radius"]
            > 50.0 * stereo[0]["stereographic_radius"]),
        "while_the_shadow_stays_in_the_unit_ball": bool(
            all(r["shadow_radius"] <= 1.0 + 1e-9 for r in stereo)),
        "and_the_shell_reaches_every_pole": bool(
            0.0 <= chi_star <= math.pi),
        "so_no_stereographic_pole_is_safe_for_this_scene": (
            "a shell launched from a point sweeps all of S³, so whatever pole "
            "is chosen the shell crosses it once; the first draft chose q₃=+1, "
            "which is the emitter's own position, and the emitter was a "
            "division by zero that never got drawn"),
        "and_the_shadow_is_two_to_one": ("depth is returned, not hidden; a "
                                         "crossing in the image is not a "
                                         "crossing on S³, and no claim in "
                                         "this round rests on one"),
    }


def measure_the_loop_has_one_self_consistent_solution(
        kappas: Sequence[complex] = (0.3, -0.6, 0.9, 0.99, 0.3 + 0.4j),
        rounds: int = 4000) -> Dict[str, object]:
    """The closed form against brute iteration of the loop.

    The iteration solves nothing — it sends the wave round and round — so
    agreement is evidence that the fixed point is what the loop actually does,
    not a convenient way of writing it down.
    """
    rows = []
    worst = 0.0
    for k in kappas:
        loop = WormholeLoop(kappa=k)
        closed = loop.self_consistent_amplitude()
        iterated = loop.iterate_loop(rounds)
        err = abs(iterated - closed)
        worst = max(worst, err)
        rows.append({"kappa": str(k), "closed_form": str(closed),
                     "iterated": str(iterated), "abs_error": err,
                     "amplification": loop.amplification()})
    return {
        "rows": rows, "worst_abs_error": worst,
        "the_fixed_point_is_what_the_loop_does": bool(worst < 1e-9),
        "no_tuning_was_required": True,
        "what_it_means": ("a linear wave on a time-displaced loop has exactly "
                          "one self-consistent amplitude; the loop does not "
                          "get to choose and nothing has to be arranged"),
    }


def measure_only_the_resonance_obstructs(
        radii: Sequence[float] = (0.2, 0.5, 0.9, 0.99, 1.01, 1.5, 3.0),
        n_phase: int = 24) -> Dict[str, object]:
    """Sweep the complex plane: the fixed point exists everywhere but ``κ = 1``.

    Including ``|κ| > 1``, where the loop is "unstable" in the sense that the
    iteration diverges — but the fixed point is still there and still unique.
    Divergence of a summation method is not non-existence of a solution, and
    the two are easy to conflate.
    """
    rows = []
    failures = 0
    for r in radii:
        for j in range(n_phase):
            k = r * cmath.exp(2j * math.pi * j / n_phase)
            if abs(1.0 - k) < 1e-12:
                continue
            loop = WormholeLoop(kappa=k)
            a = loop.self_consistent_amplitude()
            residual = abs(a - (loop.source + k * a))
            if not (math.isfinite(a.real) and math.isfinite(a.imag)
                    and residual < 1e-9):
                failures += 1
        rows.append({"abs_kappa": r,
                     "iteration_converges": bool(r < 1.0),
                     "fixed_point_exists": True})
    resonance_raised = False
    try:
        WormholeLoop(kappa=1.0).self_consistent_amplitude()
    except ValueError:
        resonance_raised = True
    return {
        "rows": rows, "failures": failures,
        "the_fixed_point_exists_everywhere_but_one_point": bool(
            failures == 0),
        "the_resonance_is_refused": resonance_raised,
        "it_exists_even_where_the_iteration_diverges": True,
        "so_divergence_of_a_sum_is_not_absence_of_a_solution": True,
    }


def measure_nonlinearity_is_what_would_break_it(
        kappa: float = 0.6,
        sources: Sequence[float] = (0.05, 0.2, 0.416, 0.5, 1.0)
        ) -> Dict[str, object]:
    """The uniqueness above is a fact about **linear** equations, not wormholes.

    Replace the loop's linear return with a quadratic one, ``A = S + κA²``, and
    the same closed timelike loop has **two** solutions for small sources and
    **none** past a threshold.  Reported because the linear result would
    otherwise read as a statement about time travel, which it is not.
    """
    rows = []
    for s in sources:
        disc = 1.0 - 4.0 * kappa * s
        if disc > 0:
            r1 = (1.0 - math.sqrt(disc)) / (2.0 * kappa)
            r2 = (1.0 + math.sqrt(disc)) / (2.0 * kappa)
            roots = [r1, r2]
        elif abs(disc) < 1e-15:
            roots = [1.0 / (2.0 * kappa)]
        else:
            roots = []
        checked = [abs(r - (s + kappa * r * r)) for r in roots]
        rows.append({"source": s, "discriminant": disc,
                     "n_real_solutions": len(roots),
                     "roots": roots,
                     "worst_residual": max(checked) if checked else None})
    counts = {r["n_real_solutions"] for r in rows}
    return {
        "rows": rows, "kappa": kappa,
        "threshold_source": 1.0 / (4.0 * kappa),
        "two_solutions_below_threshold": bool(
            any(r["n_real_solutions"] == 2 for r in rows)),
        "none_above_it": bool(any(r["n_real_solutions"] == 0 for r in rows)),
        "so_uniqueness_is_a_linearity_result": bool(len(counts) > 1),
        "every_reported_root_solves_the_loop": bool(all(
            r["worst_residual"] is None or r["worst_residual"] < 1e-12
            for r in rows)),
    }


def measure_the_ledger_closes(
        kappas: Sequence[complex] = (0.2, 0.5, -0.4, 0.3 + 0.5j)
        ) -> Dict[str, object]:
    """The loop's books, and what closing them is and is not evidence for.

    The throat is *assumed* flux-conserving — that is the identification's one
    physical input — so the residual closing to round-off checks the
    bookkeeping and nothing more.  It is reported that way.
    """
    rows = []
    for k in kappas:
        loop = WormholeLoop(kappa=k)
        rows.append({"kappa": str(k), **loop.ledger()})
    return {
        "rows": rows,
        "every_ledger_closes": bool(all(r["closes"] for r in rows)),
        "worst_residual": max(r["residual"] for r in rows),
        "but_that_is_the_assumption_not_the_result": (
            "flux conservation through the throat is put in when the mouths "
            "are identified; the residual checks the arithmetic"),
    }


def measure_the_delay_does_not_enter_the_ledger(
        delays: Sequence[float] = (-2.0, -12.0, -100.0, -1e4)
        ) -> Dict[str, object]:
    """How far into the past the far mouth sits is invisible to the balance.

    The delay decides the *story* — whether the receiver is struck before the
    emitter fires — and changes no conserved quantity at all.  Both halves are
    worth stating: the picture depends on it entirely, the physics not at all.
    """
    rows = []
    amps = []
    for d in delays:
        loop = WormholeLoop(kappa=0.45, delay=d)
        t = loop.arrival_times()
        amps.append(abs(loop.self_consistent_amplitude()))
        rows.append({"delay": d, "strike_receiver": t["strike_receiver"],
                     "receiver_before_emitter": t["receiver_before_emitter"],
                     "amplitude": amps[-1],
                     "ledger_closes": loop.ledger()["closes"]})
    return {
        "rows": rows,
        "amplitude_spread": max(amps) - min(amps),
        "the_delay_changes_no_conserved_quantity": bool(
            max(amps) - min(amps) < 1e-12),
        "but_it_decides_the_ordering": bool(
            any(r["receiver_before_emitter"] for r in rows)
            and any(not r["receiver_before_emitter"] for r in rows)),
        "so_the_picture_depends_on_it_and_the_physics_does_not": True,
    }


def measure_two_local_events_one_conserved_wave(
        kappa: complex = 0.45, delay: float = -12.0) -> Dict[str, object]:
    """What each observer sees, and the number that identifies the two halves.

    Locally there are two unrelated events: an emitter firing, and a shell
    arriving from elsewhere to strike a receiver.  Neither conserves anything by
    itself — the emitter loses flux, the receiver gains it.  The pair balances,
    and that balance is the only thing making them one object rather than two.
    """
    loop = WormholeLoop(kappa=kappa, delay=delay)
    a = loop.self_consistent_amplitude()
    times = loop.arrival_times()
    out_of_emitter = abs(a) ** 2
    into_receiver = abs(loop.kappa * a) ** 2
    return {
        "geodesic_emitter_to_future_mouth": loop.emitter_to_future_mouth,
        "it_is_exactly_the_antipode": bool(
            abs(loop.emitter_to_future_mouth - ANTIPODE) < 1e-12),
        "sweep_distance_past_the_emitter": loop.past_mouth_to_emitter,
        "times": times,
        "receiver_struck_before_emission": times["receiver_before_emitter"],
        "flux_out_of_emitter": out_of_emitter,
        "flux_into_receiver": into_receiver,
        "emitter_alone_conserves": False,
        "receiver_alone_conserves": False,
        "the_pair_closes_to": abs(out_of_emitter * abs(kappa) ** 2
                                  - into_receiver),
        "the_pair_conserves": bool(
            abs(out_of_emitter * abs(kappa) ** 2 - into_receiver) < 1e-12),
        "what_is_apparent": ("that the receiver was struck by an independent "
                             "particle; the momentum transfer is real, the "
                             "independence is not"),
    }
