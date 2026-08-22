"""Two waves on the circle slice: one driven outward, one driven inward.

Where this sits
───────────────
`circle_slice` (PR #246's v46) put one scalar wave on the great circle through a
source and its antipode, drew it as a radial height ``r = R_mid + εu(|σ|)`` in
the vacuole, and glued ``R_outer`` to ``R_inner`` so the radial direction is a
circle and the picture lives on a torus.  Its result was a negative one, and a
sharp one:

    the curve is a **graph** ``r = f(σ)``, so its radial winding number is
    identically zero.  Every outward crossing of the seam is paid for by an
    inward one.  **A height field cannot wind**, and a single wave running to
    its antipode never meets itself.

That is a statement about *one* wave.  The arc since then has needed *two* —
two mouths, two sources, an interference term that is bilinear and vanishes
unless both are present.  So the question this module asks is the one v46 left
open:

    one wave pulsing **outward** and one pulsing **inward**, both refocusing at
    the antipode — do they connect, inner to outer, and where?

The construction
────────────────
Two curves over the same circle, mirror images about the mid-radius::

    r_A(σ) = R_mid + ε u(d_A)        (driven outward)
    r_B(σ) = R_mid − ε u(d_B)        (driven inward)

with ``d_A = |σ|`` and ``d_B`` measured from ``B``'s own source, which may sit
on top of ``A``'s (``offset = 0``) or at the far side (``offset = π``).

On the torus two points at the same ``σ`` are joined by two radial arcs.  The
one that matters here is the one **through the seam** — out past ``R_outer``,
round, and back in at ``R_inner`` — of length ``gap − |δ|`` with
``δ = r_A − r_B = ε(u_A + u_B)``.  It closes when

    ``|δ| = gap`` ,   i.e.   ``ε u = gap/2``   for the co-located pair.

The answer, and the identity underneath it
──────────────────────────────────────────
**Yes, and exactly where one wave would have wrapped.**  ``ε u = gap/2`` is
*also* the condition for a single wave to cross the seam at all, because
crossing means ``R_mid + εu > R_outer``.  So:

    the amplitude at which one wave first crosses the seam and the amplitude at
    which the outward and inward pair first touch through it are **the same
    number**, at the same time and the same place.

v46's "the wave comes back inside the circle" and this module's "the two pulses
connect inner to outer" are one event described twice.  That is an identity, not
a coincidence, and it is checked as one.

Three details worth having
──────────────────────────
**The contact is on the seam.**  At threshold one curve sits at ``R_outer`` and
the other at ``R_inner``, and after gluing those are the same point.  Which is
which is not the obvious way round: the antipodal refocus is a **rarefaction**
(``u < 0``), so it is the *inward*-driven wave that bulges out to ``R_outer``
and the outward-driven one that dips to ``R_inner``.

**The contact set is an arc, not a pair of points.**  At threshold it is a
single tangency at ``σ = ±π``.  Above threshold it opens into one arc, bounded
by two genuine crossings of the two curves, on which the band between them
covers the *entire* radial circle.  That is the thing a single graph cannot do:
one wave crossing the seam still leaves every radius outside itself, while two
waves past threshold leave none.

**Meeting mid-flight is harder, not easier.**  With the sources at opposite ends
(``offset = π``) the two travelling pulses cross at the quarter points at
``t = π/2``, which looks like the natural place for a connection and is the
worst one: they *partially cancel*, and the threshold there is ``7.4×`` the
cheapest for a co-located pair and ``9.0×`` for an antipodal one.  (A first
draft of this line guessed ``4×`` from a partial scan; the measured penalties
are those.)  The cheapest connection is always at a **refocus** — and a
co-located pair reaches it at half the amplitude an antipodally-sourced one
needs, because at a refocus *both* of a co-located pair are at peak while only
one of an antipodal pair is.

The degenerate case, and what moving off it shows
─────────────────────────────────────────────────
Everything above is the **co-located** pair, and a co-located pair is the most
degenerate configuration available: both wave histories hang off one antipodal
axis, so bringing them together at one pole invites either exact overlap or
exact cancellation and tests neither.  The interesting question is the one that
needs the sources apart:

    can an inner-going branch launched from one axis meet an outer-going branch
    that has crossed the identification and re-entered at the inner boundary on
    a *different* axis — and does such a pair reach places that a like-signed
    pair cannot?

Two knobs answer it.  ``offset`` is the angle ``α`` between the two sources, and
``signs`` is the radial sense each wave is driven in::

    r_A(σ) = R_mid + s_A ε u(d_A)          d_A = |σ|
    r_B(σ) = R_mid + s_B ε u(d_B)          d_B = |σ − α|   (wrapped)

    δ = r_A − r_B = ε (s_A u_A − s_B u_B)

so **opposed** signs give the *sum* of the two fields and **like** signs give the
*difference*.  That one line is the whole asymmetry.

**There are two cases, not four.**  ``(out, out)`` and ``(in, in)`` give
identically the same ``|δ|``, and so do ``(out, in)`` and ``(in, out)`` —
exactly, not nearly.  Flipping both signs is a reflection about ``R_mid``, which
is an isometry of the glued radial circle.  So this representation *cannot*
distinguish inner–inner from outer–outer, and it is worth being clear why: the
radial direction here carries the field's **amplitude**, not its direction of
propagation.  A picture in which those two differ would have to encode
propagation in the curve, and this one does not.

**The bisector.**  ``σ = α/2`` is equidistant from both sources, so
``u_A = u_B`` there identically — at every time, at every amplitude.  For a
like-signed pair that makes ``δ ≡ 0``: the two curves are *the same curve* on
that axis and are never separated by so much as a hair, so no gain however large
carries them through the seam there.  For an opposed pair the same equality
makes ``δ = 2εu(α/2)`` — the largest it can be.  The bisector is where one pair
is maximally connected and the other is identically not.

There are **two** such axes, ``α/2`` and ``α/2 − π``, and the far one is the
cheaper of the two because it sits nearer the antipodal caustic.

**So the answer is yes, and the offset is what produces it.**  Above threshold
the opposed pair's through-seam contact opens into an arc centred on the
bisector to machine precision, on which the like-signed pair's contact set is
empty at every offset tested.  At ``α = 0`` the bisector collapses onto the
source axis and there is nothing to see — which is exactly the degeneracy, found
by measurement rather than assumed.  Turning ``α`` up slides the exclusive
connection continuously from the source axis to the quarter point, at a cost
rising from ``0.220`` to ``1.66`` — monotonically but for a ``0.02%`` turn-over
into the symmetric endpoint — and reached at the time the two outgoing pulses
cross, ``t = α/2``, to ``0.003π``.

It is paid for.  From ``α = 0.125π`` up, the globally cheapest connection sits
exactly on one of the **four axes** — either source or either antipode — is
available to *both* pairs, and costs ``1.7–3.7×`` less.  Exclusive is not the
same as cheap, and both numbers are reported.

Honest scope
────────────
Everything v46 put in is still put in.  The crossing rule is a *representation*
choice and not a derived boundary condition; nothing makes either wave
dynamically aware of the seam or of the other wave; the field is a linear scalar
on a fixed round background, so the two waves do not interact at all — they are
drawn on the same torus and the question is purely whether their images meet.
The gain is a **display** amplitude, as it was in v46, so every threshold quoted
here is in those units and is a statement about the picture rather than about a
measured field strength.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .circle_slice import ANTIPODAL_TIME, RETURN_TIME, TWO_PI, CircleSlice

__all__ = [
    "TwoWaveSlice",
    "CO_LOCATED",
    "ANTIPODAL_SOURCES",
    "OUTWARD_INWARD",
    "BOTH_OUTWARD",
    "BOTH_INWARD",
    "INWARD_OUTWARD",
    "measure_the_pair_touches_exactly_where_one_wave_wraps",
    "measure_where_the_two_pulses_connect",
    "measure_the_contact_is_an_arc_the_band_covers",
    "measure_meeting_mid_flight_is_harder",
    "measure_like_signs_are_one_case_not_two",
    "measure_the_bisector_is_degenerate_for_like_signs",
    "measure_only_the_opposed_pair_connects_on_the_bisector",
    "measure_the_offset_slides_the_connection",
]

CO_LOCATED = 0.0
ANTIPODAL_SOURCES = math.pi

#: Radial sense each wave is driven in.  ``+1`` outward, ``−1`` inward.  There
#: are only **two** distinct cases here, not four: flipping both signs is a
#: reflection about ``R_mid``, an isometry of the glued radial circle, so
#: ``BOTH_OUTWARD`` and ``BOTH_INWARD`` give identically the same ``|δ|`` — and
#: so do the two opposed orderings.  Checked at zero, not asserted.
OUTWARD_INWARD = (1.0, -1.0)
INWARD_OUTWARD = (-1.0, 1.0)
BOTH_OUTWARD = (1.0, 1.0)
BOTH_INWARD = (-1.0, -1.0)


@dataclass
class TwoWaveSlice:
    """Two scalar waves on one circle slice, driven in opposite radial senses.

    ``offset`` is the angle ``α`` between the two sources — continuous, with
    ``0`` the co-located pair (both refocusing at the same antipode) and ``π``
    sources at opposite ends.  ``signs`` is the radial sense each is driven in,
    ``+1`` outward and ``−1`` inward; the default is the opposed pair.
    """

    slice_: CircleSlice = field(default_factory=CircleSlice)
    offset: float = CO_LOCATED
    signs: Tuple[float, float] = OUTWARD_INWARD

    def __post_init__(self) -> None:
        s_a, s_b = self.signs
        if not (abs(abs(float(s_a)) - 1.0) < 1e-12
                and abs(abs(float(s_b)) - 1.0) < 1e-12):
            raise ValueError("signs must each be +1 (outward) or -1 (inward)")
        self.signs = (float(s_a), float(s_b))

    @property
    def opposed(self) -> bool:
        """Whether the two waves are driven in opposite radial senses.

        This is the only thing about ``signs`` that changes ``|δ|``: opposed
        gives ``ε|u_A + u_B|`` and like-signed gives ``ε|u_A − u_B|``.
        """
        return self.signs[0] * self.signs[1] < 0.0

    # ── geometry of the base circle ─────────────────────────────────────────
    @property
    def sigma(self) -> np.ndarray:
        """The drawn parameterisation, with ``σ = ±π`` both present."""
        return self.slice_.sigma

    # ── the two axes equidistant from both sources ──────────────────────────
    @property
    def bisector(self) -> float:
        """``σ = α/2`` — equidistant from both sources going the short way.

        ``u_A = u_B`` there identically, so a like-signed pair has ``δ ≡ 0``
        and an opposed one has ``δ = 2εu(α/2)``: the axis where one pair is
        maximally connected and the other is identically not.
        """
        return 0.5 * float(self.offset)

    @property
    def far_bisector(self) -> float:
        """``σ = α/2 − π`` — the other one, equidistant going the long way.

        There are always two, and this is the cheaper of the pair for an
        opposed configuration because it sits nearer the antipodal caustic.
        """
        return 0.5 * float(self.offset) - math.pi

    def fields_at(self, sigma) -> Tuple[np.ndarray, np.ndarray]:
        """``(u_A, u_B)`` at arbitrary ``σ`` — off the drawing grid.

        The bisector is only a grid point by accident, and the degeneracy
        claimed for it is exact, so it has to be evaluated where it actually is
        rather than at the nearest sample.
        """
        s = np.atleast_1d(np.asarray(sigma, dtype=float))
        d_a, d_b = self._distances(s)
        at = self.slice_.sim.field_at_distance
        return at(d_a), at(d_b)

    @property
    def sigma_closed(self) -> np.ndarray:
        """The circle with the duplicated endpoint removed.

        Counting anything on the circle has to use this one.  A first draft
        counted contact regions on the drawn array and reported ``2`` at the
        tangency, which was the single point ``σ = ±π`` seen twice.
        """
        return self.slice_.sigma[:-1]

    @property
    def gap(self) -> float:
        return self.slice_.bulk.gap

    def _distances(self, sigma: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Geodesic distance from each source, both wrapped.

        ``A`` sits at ``σ = 0``, so ``|σ|`` is right on the drawing grid, but
        the bisectors are asked for off it and one of them is outside
        ``[−π, π]``; both distances are wrapped so the two sources are treated
        by exactly the same expression.
        """
        d_a = np.abs(np.mod(sigma + math.pi, TWO_PI) - math.pi)
        d_b = np.abs(np.mod(sigma - self.offset + math.pi, TWO_PI) - math.pi)
        return d_a, d_b

    # ── the two fields, and the two radii ───────────────────────────────────
    def fields(self, closed: bool = False) -> Tuple[np.ndarray, np.ndarray]:
        s = self.sigma_closed if closed else self.sigma
        d_a, d_b = self._distances(s)
        at = self.slice_.sim.field_at_distance
        return at(d_a), at(d_b)

    def radii(self, gain: Optional[float] = None,
              closed: bool = False) -> Tuple[np.ndarray, np.ndarray]:
        """``(r_A, r_B)`` before wrapping, each driven in its own sense."""
        eps = self.slice_.gain if gain is None else float(gain)
        u_a, u_b = self.fields(closed=closed)
        mid = self.slice_.bulk.r_mid
        s_a, s_b = self.signs
        return mid + s_a * eps * u_a, mid + s_b * eps * u_b

    def separation(self, gain: Optional[float] = None,
                   closed: bool = False) -> np.ndarray:
        """``δ = r_A − r_B``.  The pair spans this much of the radial circle."""
        r_a, r_b = self.radii(gain=gain, closed=closed)
        return r_a - r_b

    def outward_gap(self, gain: Optional[float] = None,
                    closed: bool = False) -> np.ndarray:
        """``gap − |δ|`` — the arc from one curve to the other **through the
        seam**.  This is the quantity that closes when they connect."""
        return self.gap - np.abs(self.separation(gain=gain, closed=closed))

    def covered_fraction(self, gain: Optional[float] = None) -> float:
        """``max|δ| / gap``.  Reaches ``1`` exactly at contact."""
        return float(np.max(np.abs(self.separation(gain=gain, closed=True)))
                     / self.gap)

    # ── wrapped curves, for drawing ─────────────────────────────────────────
    def segments(self, gain: Optional[float] = None
                 ) -> Tuple[List[np.ndarray], List[np.ndarray]]:
        out = []
        for r in self.radii(gain=gain):
            drawn, sheet = self.slice_.bulk.wrap(r)
            xy = np.stack([drawn * np.cos(self.sigma),
                           drawn * np.sin(self.sigma)], axis=-1)
            cuts = np.nonzero(np.diff(sheet) != 0)[0]
            pieces, start = [], 0
            for c in cuts:
                if c + 1 > start:
                    pieces.append(xy[start:c + 1])
                start = c + 1
            pieces.append(xy[start:])
            out.append([p for p in pieces if len(p) > 1])
        return out[0], out[1]

    # ── contact ─────────────────────────────────────────────────────────────
    def contact_gain(self) -> float:
        """The display gain at which the pair first touches through the seam.

        ``|δ| = gap`` with ``δ = ε(s_A u_A − s_B u_B)``, so this is
        ``gap / max|s_A u_A − s_B u_B|`` at the current time — the *sum* of the
        two fields for an opposed pair and their *difference* for a like-signed
        one.
        """
        u_a, u_b = self.fields(closed=True)
        s_a, s_b = self.signs
        peak = float(np.max(np.abs(s_a * u_a - s_b * u_b)))
        return float("inf") if peak <= 0.0 else self.gap / peak

    def separation_at(self, sigma, gain: Optional[float] = None) -> np.ndarray:
        """``δ`` at arbitrary ``σ``, off the drawing grid."""
        eps = self.slice_.gain if gain is None else float(gain)
        u_a, u_b = self.fields_at(sigma)
        s_a, s_b = self.signs
        return eps * (s_a * u_a - s_b * u_b)

    def contact_gain_at(self, sigma) -> float:
        """The gain at which the pair spans the whole radial circle at ``σ``.

        ``inf`` when the two curves coincide there — which is what happens on a
        bisector for a like-signed pair, at every time and every amplitude.
        """
        d = float(np.max(np.abs(self.separation_at(sigma, gain=1.0))))
        return float("inf") if d <= 0.0 else self.gap / d

    def contact(self, gain: Optional[float] = None) -> Dict[str, object]:
        """Whether, where, and at what radius the two curves meet."""
        s = self.sigma_closed
        delta = self.separation(gain=gain, closed=True)
        over = np.abs(delta) >= self.gap * (1.0 - 1e-12)
        wrapped = np.concatenate([over, over[:1]])
        flips = int(np.count_nonzero(np.diff(wrapped.astype(int))))
        k = int(np.argmax(np.abs(delta)))
        r_a, r_b = self.radii(gain=gain, closed=True)
        return {
            "connected": bool(over.any()),
            "covered_fraction": float(np.max(np.abs(delta)) / self.gap),
            "arcs": flips // 2,
            "crossings": flips,
            "sigma_of_closest_approach": float(s[k]),
            "sigma_over_pi": float(s[k] / math.pi),
            "outward_gap_at_closest": float(self.gap - abs(delta[k])),
            "radius_a": float(r_a[k]),
            "radius_b": float(r_b[k]),
            "r_inner": self.slice_.bulk.r_inner,
            "r_outer": self.slice_.bulk.r_outer,
            "the_refocus_is_a_rarefaction": bool(delta[k] < 0.0),
        }

    # ── time ────────────────────────────────────────────────────────────────
    def scan(self, samples: int = 240,
             t_end: float = RETURN_TIME) -> List[Dict[str, float]]:
        """``contact_gain`` over a whole return period."""
        rows = []
        self.slice_.reset()
        for i in range(int(samples)):
            t = (i + 1) * t_end / int(samples)
            self.slice_.advance_to(t)
            u_a, u_b = self.fields(closed=True)
            tot = u_a + u_b
            k = int(np.argmax(np.abs(tot)))
            rows.append({"t": float(t), "peak": float(abs(tot[k])),
                         "sigma": float(self.sigma_closed[k]),
                         "contact_gain": self.contact_gain()})
        self.slice_.reset()
        return rows


# ════════════════════════════════════════════════════════════════════════════
# MEASUREMENTS
# ════════════════════════════════════════════════════════════════════════════
def measure_the_pair_touches_exactly_where_one_wave_wraps(
        pair: Optional[TwoWaveSlice] = None,
        samples: int = 400) -> Dict[str, object]:
    """**The identity.**  One wave wrapping and two waves touching are one event.

    A single wave crosses the seam when ``R_mid + εu > R_outer``, i.e.
    ``εu > gap/2``.  The co-located pair spans ``|δ| = 2ε|u|`` of the radial
    circle and touches through the seam when that reaches ``gap`` — i.e. when
    ``εu = gap/2``.  Same condition, so the same threshold gain, at the same
    time and the same ``σ``.

    Checked against the run's own peak rather than asserted: both thresholds are
    computed independently and compared.
    """
    pair = pair or TwoWaveSlice()
    sl = pair.slice_
    sl.reset()
    best, when, where = 0.0, 0.0, 0.0
    for i in range(int(samples)):
        t = (i + 1) * RETURN_TIME / int(samples)
        sl.advance_to(t)
        u = sl.sim.field_at_distance(np.abs(pair.sigma_closed))
        k = int(np.argmax(np.abs(u)))
        if abs(u[k]) > best:
            best, when, where = float(abs(u[k])), float(t), float(
                pair.sigma_closed[k])
    sl.reset()
    half = 0.5 * pair.gap
    wrap_gain = half / best                       # one wave reaches R_outer
    sl.advance_to(when)
    touch_gain = pair.contact_gain()              # two waves span the gap
    sl.reset()
    return {
        "run_peak": best,
        "at_time": when,
        "at_sigma": where,
        "single_wave_wrap_gain": float(wrap_gain),
        "pair_contact_gain": float(touch_gain),
        "relative_difference": float(abs(touch_gain / wrap_gain - 1.0)),
        "they_are_the_same_threshold": bool(
            abs(touch_gain / wrap_gain - 1.0) < 1e-12),
        "why": "both are eps*u = gap/2 -- crossing the seam and spanning the "
               "radial period are the same inequality",
    }


def measure_where_the_two_pulses_connect(
        pair: Optional[TwoWaveSlice] = None) -> Dict[str, object]:
    """**The answer.**  At the antipode, on the seam, at the refocus.

    Driven to the contact threshold at ``t = π``, one curve sits exactly at
    ``R_outer`` and the other exactly at ``R_inner`` — the same point once the
    boundaries are glued — at ``σ = ±π``.

    And the assignment is the other way round from the obvious guess.  The
    antipodal refocus is a **rarefaction**, so it is the *inward*-driven wave
    that reaches ``R_outer`` and the outward-driven one that dips to
    ``R_inner``.
    """
    pair = pair or TwoWaveSlice()
    sl = pair.slice_
    sl.reset()
    sl.advance_to(ANTIPODAL_TIME)
    gain = pair.contact_gain()
    got = pair.contact(gain=gain)
    u_a, _ = pair.fields(closed=True)
    k = int(np.argmax(np.abs(u_a)))
    radii = sorted((got["radius_a"], got["radius_b"]))
    sl.reset()
    return {
        "time": ANTIPODAL_TIME,
        "contact_gain": float(gain),
        "sigma": got["sigma_of_closest_approach"],
        "sigma_over_pi": got["sigma_over_pi"],
        "refocus_amplitude": float(u_a[k]),
        "it_is_a_rarefaction": bool(u_a[k] < 0.0),
        "radius_of_the_outward_wave": got["radius_a"],
        "radius_of_the_inward_wave": got["radius_b"],
        "r_inner": got["r_inner"],
        "r_outer": got["r_outer"],
        "distance_to_the_seam": float(
            max(abs(radii[0] - got["r_inner"]),
                abs(radii[1] - got["r_outer"]))),
        "connected": bool(got["connected"]),
        "the_contact_is_on_the_seam": bool(
            max(abs(radii[0] - got["r_inner"]),
                abs(radii[1] - got["r_outer"])) < 1e-9),
        "the_contact_is_at_the_antipode": bool(
            abs(abs(got["sigma_over_pi"]) - 1.0) < 1e-9),
        "the_inward_wave_is_the_one_that_reaches_r_outer": bool(
            got["radius_b"] > got["radius_a"]),
    }


def measure_the_contact_is_an_arc_the_band_covers(
        pair: Optional[TwoWaveSlice] = None,
        multipliers: Sequence[float] = (0.98, 0.999, 1.0, 1.001, 1.05, 1.3)
        ) -> Dict[str, object]:
    """What two waves can do that one cannot.

    Below threshold the two curves never meet.  **At** threshold they touch once
    — a tangency at ``σ = ±π``.  **Above** it that point opens into a single arc,
    bounded by two genuine crossings, on which the band between the two curves
    covers the whole radial circle: at those ``σ`` there is no radius left
    outside the pair.

    A single wave past its own wrap threshold does nothing of the kind.  It
    crosses the seam, reappears inside, and still leaves every radius outside
    itself — its winding number is zero because it is a graph.  Two graphs
    bound a band, and a band can be radially surjective.
    """
    pair = pair or TwoWaveSlice()
    sl = pair.slice_
    sl.reset()
    sl.advance_to(ANTIPODAL_TIME)
    base = pair.contact_gain()
    rows = []
    for m in multipliers:
        got = pair.contact(gain=base * float(m))
        rows.append({"gain_over_threshold": float(m),
                     "covered_fraction": got["covered_fraction"],
                     "connected": got["connected"],
                     "arcs": got["arcs"],
                     "crossings": got["crossings"]})
    sl.reset()
    below = [r for r in rows if r["gain_over_threshold"] < 1.0]
    above = [r for r in rows if r["gain_over_threshold"] > 1.0]
    at = [r for r in rows if r["gain_over_threshold"] == 1.0]
    return {
        "rows": rows,
        "nothing_connects_below_threshold": bool(
            all(not r["connected"] for r in below)),
        "it_touches_at_threshold": bool(all(r["connected"] for r in at)),
        "one_arc_above": bool(all(r["arcs"] == 1 for r in above)),
        "the_covered_fraction_tracks_the_gain": float(max(
            abs(r["covered_fraction"] / r["gain_over_threshold"] - 1.0)
            for r in rows)),
        "what_one_wave_cannot_do": "a single wave is a graph, so its radial "
                                   "winding is zero and every radius stays "
                                   "outside it however hard it is driven; two "
                                   "graphs bound a band, and past threshold "
                                   "that band leaves no radius out",
    }


def measure_meeting_mid_flight_is_harder(
        samples: int = 240) -> Dict[str, object]:
    """Antipodal sources meet at the quarter points — and partly cancel there.

    With ``B`` launched from the far side, the two travelling pulses cross at
    ``σ = ±π/2`` at ``t = π/2``.  That looks like the natural place for an
    inner-to-outer connection and it is the **worst** one: the pulses partially
    cancel, so ``|u_A + u_B|`` is *smaller* there than either pulse alone, and
    the contact threshold is about ``4×`` what a refocus needs.

    The cheapest connection is always at a refocus, in both configurations.
    """
    out: Dict[str, object] = {}
    for name, offset in (("co_located", CO_LOCATED),
                         ("antipodal_sources", ANTIPODAL_SOURCES)):
        pair = TwoWaveSlice(offset=offset)
        rows = pair.scan(samples=samples)
        best = min(rows, key=lambda r: r["contact_gain"])
        mid = min(rows, key=lambda r: abs(r["t"] - 0.5 * ANTIPODAL_TIME))
        out[name] = {
            "cheapest_contact_gain": best["contact_gain"],
            "at_time": best["t"],
            "at_time_over_pi": best["t"] / math.pi,
            "at_sigma_over_pi": best["sigma"] / math.pi,
            "mid_flight_contact_gain": mid["contact_gain"],
            "mid_flight_penalty": mid["contact_gain"] / best["contact_gain"],
        }
    co = out["co_located"]
    anti = out["antipodal_sources"]
    return {
        **out,
        "both_are_cheapest_at_a_refocus": bool(
            abs(co["at_time_over_pi"] - round(co["at_time_over_pi"])) < 0.05
            and abs(anti["at_time_over_pi"]
                    - round(anti["at_time_over_pi"])) < 0.05),
        "mid_flight_is_harder": bool(co["mid_flight_penalty"] > 1.0
                                     and anti["mid_flight_penalty"] > 1.0),
        "worst_penalty": float(max(co["mid_flight_penalty"],
                                   anti["mid_flight_penalty"])),
        "antipodal_over_co_located": float(anti["cheapest_contact_gain"]
                                           / co["cheapest_contact_gain"]),
        "co_located_is_about_twice_as_cheap": bool(
            1.8 < anti["cheapest_contact_gain"]
            / co["cheapest_contact_gain"] < 2.2),
    }


# ════════════════════════════════════════════════════════════════════════════
# OFF THE DEGENERATE AXIS — what the offset and the signs buy
# ════════════════════════════════════════════════════════════════════════════
def measure_like_signs_are_one_case_not_two(
        offsets: Sequence[float] = (0.0, 0.15, 0.30, 0.50, 0.70, 0.85, 1.0),
        times: Sequence[float] = (0.5, 1.0, 1.5)) -> Dict[str, object]:
    """There are **two** configurations here, not four — and that is a limit.

    ``(out, out)`` and ``(in, in)`` give identically the same ``|δ|``, and so do
    the two opposed orderings.  Flipping both signs is a reflection about
    ``R_mid``, which is an isometry of the glued radial circle, so the picture
    cannot tell inner–inner from outer–outer.

    Worth stating as a **limitation of the representation** rather than as a
    result about waves: the radial direction here carries the field's
    *amplitude*, not its direction of propagation.  Distinguishing an
    inner-going branch from an outer-going one would need a curve that encodes
    propagation, and this one does not.

    Two residues are reported, because they are different claims.  As a
    *difference of fields* the agreement is exact — the same bits.  Through the
    *drawn radii* it is one ulp of ``R_mid``, because ``(R_mid + εu_A) −
    (R_mid + εu_B)`` and ``(R_mid − εu_A) − (R_mid − εu_B)`` round differently.
    A first draft asserted the drawn residue at zero and it failed; the identity
    was never in doubt, only which quantity it was being asserted about.
    """
    worst_like, worst_opposed = 0.0, 0.0
    worst_like_field, worst_opposed_field = 0.0, 0.0
    rows = []
    for f in offsets:
        alpha = float(f) * math.pi
        oo = TwoWaveSlice(offset=alpha, signs=BOTH_OUTWARD)
        ii = TwoWaveSlice(offset=alpha, signs=BOTH_INWARD)
        oi = TwoWaveSlice(offset=alpha, signs=OUTWARD_INWARD)
        io = TwoWaveSlice(offset=alpha, signs=INWARD_OUTWARD)
        # one clock drives them all: they share nothing but the sim they read
        sl = oo.slice_
        for p in (ii, oi, io):
            p.slice_ = sl
        sl.reset()
        d_like = d_opp = f_like = f_opp = 0.0
        for tf in times:
            sl.advance_to(float(tf) * math.pi)
            u_a, u_b = oo.fields(closed=True)
            # the identity, as a difference of fields: the same bits
            f_like = max(f_like, float(np.max(np.abs(
                np.abs(u_a - u_b) - np.abs(-u_a + u_b)))))
            f_opp = max(f_opp, float(np.max(np.abs(
                np.abs(u_a + u_b) - np.abs(-u_a - u_b)))))
            # ...and through the drawn radii, where R_mid costs an ulp
            d_like = max(d_like, float(np.max(np.abs(
                np.abs(oo.separation(gain=1.0, closed=True))
                - np.abs(ii.separation(gain=1.0, closed=True))))))
            d_opp = max(d_opp, float(np.max(np.abs(
                np.abs(oi.separation(gain=1.0, closed=True))
                - np.abs(io.separation(gain=1.0, closed=True))))))
        sl.reset()
        worst_like = max(worst_like, d_like)
        worst_opposed = max(worst_opposed, d_opp)
        worst_like_field = max(worst_like_field, f_like)
        worst_opposed_field = max(worst_opposed_field, f_opp)
        rows.append({"offset_over_pi": float(f),
                     "out_out_vs_in_in": d_like,
                     "out_in_vs_in_out": d_opp,
                     "out_out_vs_in_in_as_fields": f_like,
                     "out_in_vs_in_out_as_fields": f_opp})
    ulp = float(np.spacing(TwoWaveSlice().slice_.bulk.r_mid))
    return {
        "rows": rows,
        "worst_out_out_vs_in_in": worst_like,
        "worst_out_in_vs_in_out": worst_opposed,
        "worst_as_fields": max(worst_like_field, worst_opposed_field),
        "one_ulp_of_r_mid": ulp,
        "drawn_residue_in_ulps": max(worst_like, worst_opposed) / ulp,
        "out_out_is_in_in": bool(worst_like_field == 0.0),
        "the_two_opposed_orderings_agree": bool(worst_opposed_field == 0.0),
        "there_are_two_cases_not_four": bool(worst_like_field == 0.0
                                             and worst_opposed_field == 0.0),
        "the_drawn_residue_is_one_ulp_of_r_mid": bool(
            max(worst_like, worst_opposed) <= 2.0 * ulp),
        "why": "flipping both signs reflects the picture about R_mid, an "
               "isometry of the glued radial circle",
        "the_limitation": "the radial direction carries field amplitude, not "
                          "direction of propagation, so this representation "
                          "cannot distinguish inner-inner from outer-outer",
    }


def measure_the_bisector_is_degenerate_for_like_signs(
        offsets: Sequence[float] = (0.0, 0.15, 0.30, 0.50, 0.70, 0.85, 1.0),
        samples: int = 240) -> Dict[str, object]:
    """``σ = α/2`` is equidistant from both sources, so ``u_A = u_B`` on it.

    For a **like-signed** pair that makes ``δ ≡ 0``: the two curves are the same
    curve on that axis, at every time and every amplitude, so no gain carries
    them through the seam there.  For an **opposed** pair the same equality
    makes ``δ = 2εu(α/2)`` — as large as it can be.

    There are two such axes, ``α/2`` and ``α/2 − π``, and both are checked.  The
    residue is at the floating-point floor of the wrapped distance rather than
    exactly zero, so the tolerance is stated relative to the field's own scale;
    that residue is the ``mod``, not the field.
    """
    rows = []
    worst_like, floor = 0.0, 0.0
    for f in offsets:
        alpha = float(f) * math.pi
        row = {"offset_over_pi": float(f)}
        for name, sign in (("near", +1), ("far", -1)):
            opp = TwoWaveSlice(offset=alpha, signs=OUTWARD_INWARD)
            like = TwoWaveSlice(offset=alpha, signs=BOTH_OUTWARD)
            like.slice_ = opp.slice_
            axis = opp.bisector if sign > 0 else opp.far_bisector
            sl = opp.slice_
            sl.reset()
            d_like, d_opp, scale, at = 0.0, 0.0, 0.0, 0.0
            for i in range(int(samples)):
                t = (i + 1) * RETURN_TIME / int(samples)
                sl.advance_to(t)
                u_a, _ = opp.fields_at(axis)
                d_like = max(d_like, float(np.max(np.abs(
                    like.separation_at(axis, gain=1.0)))))
                v = float(np.max(np.abs(opp.separation_at(axis, gain=1.0))))
                scale = max(scale, float(np.max(np.abs(u_a))))
                if v > d_opp:
                    d_opp, at = v, t
            sl.reset()
            worst_like = max(worst_like, d_like)
            floor = max(floor, d_like / max(scale, 1e-30))
            row[f"{name}_bisector_over_pi"] = axis / math.pi
            row[f"{name}_like_signed_separation"] = d_like
            row[f"{name}_opposed_threshold"] = (
                opp.gap / d_opp if d_opp > 0 else float("inf"))
            row[f"{name}_at_time_over_pi"] = at / math.pi
        rows.append(row)
    far_cheaper = [r for r in rows
                   if r["offset_over_pi"] not in (0.0, 1.0)]
    return {
        "rows": rows,
        "worst_like_signed_separation": worst_like,
        "worst_relative_residue": floor,
        "the_like_signed_pair_never_separates_on_a_bisector": bool(
            floor < 1e-13),
        "the_opposed_pair_always_does": bool(
            all(math.isfinite(r[f"{n}_opposed_threshold"])
                for r in rows for n in ("near", "far"))),
        "the_far_bisector_is_the_cheaper_one": bool(
            all(r["far_opposed_threshold"] < r["near_opposed_threshold"]
                for r in far_cheaper)),
        "why_the_far_one_is_cheaper": "it sits nearer the antipodal caustic, "
                                      "where the field is amplified",
        "why": "sigma = alpha/2 is equidistant from both sources, so u_A = u_B "
               "there identically; like signs subtract that to zero and "
               "opposed signs add it to twice one field",
    }


def measure_only_the_opposed_pair_connects_on_the_bisector(
        offsets: Sequence[float] = (0.15, 0.30, 0.50, 0.70, 1.0),
        drive: float = 1.15, samples: int = 400) -> Dict[str, object]:
    """**The answer to the question the offset was added for.**

    Drive each configuration to ``1.15×`` the opposed pair's own bisector
    threshold, at the time that threshold is reached, and look at where the two
    curves span the whole radial circle.

    The opposed pair's contact set is a single arc **centred on the bisector**
    — to machine zero, at every offset — and the like-signed pair's contact set
    is *empty on that arc*.  So yes: an opposed pair possesses off-antipodal
    connections that a like-signed pair does not possess at any amplitude, and
    the offset is what produces them.  At ``α = 0`` the bisector collapses onto
    the source axis and there is nothing off-axis to find, which is the
    degeneracy, measured rather than assumed.
    """
    rows = []
    worst_centre, like_hits = 0.0, 0
    for f in offsets:
        alpha = float(f) * math.pi
        opp = TwoWaveSlice(offset=alpha, signs=OUTWARD_INWARD)
        like = TwoWaveSlice(offset=alpha, signs=BOTH_OUTWARD)
        like.slice_ = opp.slice_
        sl, axis = opp.slice_, opp.bisector
        sl.reset()
        best, at = 0.0, 0.0
        for i in range(int(samples)):
            t = (i + 1) * RETURN_TIME / int(samples)
            sl.advance_to(t)
            v = float(np.max(np.abs(opp.separation_at(axis, gain=1.0))))
            if v > best:
                best, at = v, t
        thr = opp.gap / best
        sl.reset()                 # the clock only runs forwards; `at` is behind
        sl.advance_to(at)
        gain = float(drive) * thr
        s = opp.sigma_closed
        d_sigma = TWO_PI / len(s)
        on_o = np.abs(opp.separation(gain=gain, closed=True)) >= opp.gap
        on_l = np.abs(like.separation(gain=gain, closed=True)) >= opp.gap
        near = np.nonzero(on_o & (np.abs(s - axis) < 0.5))[0]
        centre = (0.5 * (float(s[near].min()) + float(s[near].max()))
                  if len(near) else float("nan"))
        sl.reset()
        worst_centre = max(worst_centre, abs(centre - axis))
        like_hits += int(on_l[near].sum()) if len(near) else 0
        rows.append({
            "offset_over_pi": float(f),
            "bisector_over_pi": axis / math.pi,
            "threshold": thr,
            "at_time_over_pi": at / math.pi,
            "opposed_arc": len(near) * d_sigma,
            "arc_centre_over_pi": centre / math.pi,
            "centre_minus_bisector": centre - axis,
            "like_signed_samples_on_that_arc": (
                int(on_l[near].sum()) if len(near) else 0),
            "distance_to_the_nearest_source": min(
                abs(axis), abs(axis - alpha)) / math.pi,
            "distance_to_the_nearest_antipode": min(
                abs(abs(axis) - math.pi),
                abs(abs(axis - alpha) - math.pi)) / math.pi,
        })
    return {
        "drive_over_threshold": float(drive),
        "rows": rows,
        "worst_centre_offset": worst_centre,
        "the_arc_is_centred_on_the_bisector": bool(worst_centre < 1e-12),
        "every_offset_opens_an_arc": bool(all(r["opposed_arc"] > 0.0
                                              for r in rows)),
        "the_like_signed_pair_reaches_none_of_it": bool(like_hits == 0),
        "it_is_off_both_axes": bool(
            all(r["distance_to_the_nearest_source"] > 0.0
                and r["distance_to_the_nearest_antipode"] > 0.0
                for r in rows)),
        "answer": "yes -- an opposed (inner-outer) pair connects through the "
                  "seam on an arc centred on the bisector between the two "
                  "source axes, which is off both the sources and their "
                  "antipodes, and on which a like-signed pair is identically "
                  "coincident and so connects at no amplitude at all",
    }


def measure_the_offset_slides_the_connection(
        n: int = 49, samples: int = 400) -> Dict[str, object]:
    """What the slider does: the exclusive connection moves, and gets dearer.

    Sweeping ``α`` slides the bisector continuously from the source axis to the
    quarter point and raises the threshold there from ``0.220`` to ``1.66``,
    monotonically — with one measured exception, a turn-over of about ``0.02%``
    into the symmetric endpoint ``α = π``, visible only on a sweep fine enough
    to sample it.

    That number has a history worth keeping.  A first pass put the dip at
    ``0.08%``, four times larger, and "confirmed" it by refining the *time* grid
    fourfold and watching it fail to shrink.  It failed to shrink because time
    was never the problem: that pass evaluated the bisector at the nearest point
    of the ``σ`` grid, and at ``α = 0.958π`` the bisector falls exactly halfway
    between two samples.  Refining the axis a discrepancy does not live on will
    always leave it standing, and leaving it standing is not evidence.  The
    bisector is evaluated off-grid here, at the angle it actually has, and what
    survives that is one fifth the size and still real.

    The timing is the pulse-crossing time.  The two outgoing pulses reach the
    bisector together at ``t = α/2``; the measured peak leads that by a constant
    ``0.107`` rad, which is where the pulse's own extremum sits relative to its
    geodesic arrival — checked here against a single wave rather than fitted.

    And exclusivity is not cheapness.  From ``α = 0.125π`` upward the globally
    cheapest connection sits **exactly on one of the four axes** — either
    source, or either antipode — is available to *both* pairs, and costs
    ``1.7–3.7×`` less than the bisector one.  Below that it drifts off axis, by
    at most ``0.03π``, while the two return focuses are still close enough to
    blend.  The pinning offset is *measured*, not derived:
    a first draft attributed the drift to the pulse width and the numbers do not
    support that, so what is reported is where it pins and by how much it
    strays.
    """
    # where a single pulse's extremum sits relative to its geodesic arrival
    probe = TwoWaveSlice()
    sl = probe.slice_
    sl.reset()
    leads = {d: (0.0, 0.0) for d in (0.25, 0.5, 0.75, 1.0, 1.5)}
    for i in range(2 * int(samples)):
        t = (i + 1) * ANTIPODAL_TIME / (2 * int(samples))
        sl.advance_to(t)
        for d in leads:
            v = abs(float(sl.sim.field_at_distance(np.array([d]))[0]))
            if v > leads[d][0]:
                leads[d] = (v, t)
    sl.reset()
    lead = [t - d for d, (_, t) in leads.items()]
    mean_lead = sum(lead) / len(lead)

    rows, drops = [], []
    prev = -math.inf
    for i in range(int(n)):
        f = i / (int(n) - 1)
        alpha = f * math.pi
        opp = TwoWaveSlice(offset=alpha, signs=OUTWARD_INWARD)
        sl, axis = opp.slice_, opp.bisector
        sl.reset()
        best, at, gbest, gsig = 0.0, 0.0, math.inf, 0.0
        for j in range(int(samples)):
            t = (j + 1) * RETURN_TIME / int(samples)
            sl.advance_to(t)
            v = float(np.max(np.abs(opp.separation_at(axis, gain=1.0))))
            if v > best:
                best, at = v, t
            g = opp.contact_gain()
            if g < gbest:
                gbest = g
                k = int(np.argmax(np.abs(
                    opp.separation(gain=1.0, closed=True))))
                gsig = float(opp.sigma_closed[k])
        sl.reset()
        thr = opp.gap / best
        if thr < prev:
            drops.append({"offset_over_pi": f, "drop": prev - thr,
                          "relative": (prev - thr) / prev})
        prev = thr
        rows.append({
            "offset_over_pi": f,
            "bisector_over_pi": axis / math.pi,
            "bisector_threshold": thr,
            "at_time_over_pi": at / math.pi,
            "predicted_time_over_pi": (0.5 * alpha + mean_lead) / math.pi,
            "cheapest_threshold_anywhere": gbest,
            "cheapest_at_sigma_over_pi": gsig / math.pi,
            "price_of_exclusivity": thr / gbest,
        })
    # The crossing law only means anything once the crossing is clear of the
    # launch: at alpha = 0.104 pi it predicts t = 0.018 pi, inside the pulse's
    # own width, and the measured peak is the first frame.  That row is named
    # rather than absorbed into the tolerance.
    width = float(probe.slice_.sim.pulse_width)
    timed = [r for r in rows
             if r["predicted_time_over_pi"] * math.pi > width]
    too_early = [r["offset_over_pi"] for r in rows
                 if r["offset_over_pi"] > 0.1 and r not in timed]
    timing = max(abs(r["at_time_over_pi"] - r["predicted_time_over_pi"])
                 for r in timed)
    # the price is only a comparison once the cheapest point has pinned to an
    # axis; below that the bisector *is* roughly the cheapest point and the
    # ratio is near 1 by construction rather than by result
    prices = [r["price_of_exclusivity"] for r in rows
              if r["offset_over_pi"] >= 0.125]
    # "an axis" means one of the four: either source, or either antipode.  A
    # first pass compared only against A's, which the sweep promptly falsified
    # -- above alpha = 0.66 pi the winner alternates between A's axis and B's,
    # the two being degenerate by symmetry and separated only by which grid
    # index argmax happens to reach first.
    grid = TWO_PI / len(probe.sigma_closed)

    def off_axis(r):
        alpha = r["offset_over_pi"] * math.pi
        s = r["cheapest_at_sigma_over_pi"] * math.pi
        return min(abs((s - a + math.pi) % TWO_PI - math.pi)
                   for a in (0.0, math.pi, alpha, alpha - math.pi))

    for r in rows:
        r["cheapest_is_off_axis_by_over_pi"] = off_axis(r) / math.pi
    strays = [r for r in rows if off_axis(r) > 0.5 * grid]
    pins_from = (min(r["offset_over_pi"] for r in rows
                     if all(off_axis(q) <= 0.5 * grid for q in rows
                            if q["offset_over_pi"] >= r["offset_over_pi"]))
                 if any(off_axis(r) <= 0.5 * grid for r in rows) else None)
    apart = [r for r in rows
             if pins_from is not None and r["offset_over_pi"] >= pins_from]
    return {
        "rows": rows,
        "pulse_lead": mean_lead,
        "pulse_width": float(probe.slice_.sim.pulse_width),
        "pulse_lead_spread": max(lead) - min(lead),
        "threshold_at_zero_offset": rows[0]["bisector_threshold"],
        "threshold_at_pi": rows[-1]["bisector_threshold"],
        "threshold_range": (rows[-1]["bisector_threshold"]
                            / rows[0]["bisector_threshold"]),
        "worst_timing_error_over_pi": timing,
        "offsets_whose_crossing_is_inside_the_launch": too_early,
        "the_timing_is_the_pulse_crossing": bool(timing < 0.01),
        "non_monotone_steps": drops,
        "it_rises_monotonically_except_at_the_endpoint": bool(
            all(d["offset_over_pi"] > 0.99 for d in drops)
            and all(d["relative"] < 0.01 for d in drops)),
        "price_of_exclusivity_range": (min(prices), max(prices)),
        "it_pins_to_an_axis_from_offset_over_pi": pins_from,
        "offsets_at_or_above_that": len(apart),
        "worst_drift_below_it": max((off_axis(r) / math.pi for r in strays),
                                    default=0.0),
        "offsets_where_it_drifts": [r["offset_over_pi"] for r in strays],
        "the_cheapest_connection_sits_on_one_of_the_four_axes": bool(
            apart and all(off_axis(r) <= 0.5 * grid for r in apart)),
        "it_drifts_off_axis_only_at_small_offset": bool(
            all(r["offset_over_pi"] < 0.125 for r in strays)),
        "the_drift_is_small": bool(
            max((off_axis(r) for r in strays), default=0.0) < 0.04 * math.pi),
        "exclusive_is_not_cheap": bool(min(prices) > 1.5),
    }
