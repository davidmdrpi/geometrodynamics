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
    "measure_the_pair_touches_exactly_where_one_wave_wraps",
    "measure_where_the_two_pulses_connect",
    "measure_the_contact_is_an_arc_the_band_covers",
    "measure_meeting_mid_flight_is_harder",
]

CO_LOCATED = 0.0
ANTIPODAL_SOURCES = math.pi


@dataclass
class TwoWaveSlice:
    """Two scalar waves on one circle slice, driven in opposite radial senses.

    ``offset`` is where ``B``'s source sits relative to ``A``'s: ``0`` for the
    co-located pair (both refocusing at the same antipode) and ``π`` for sources
    at opposite ends of the circle.
    """

    slice_: CircleSlice = field(default_factory=CircleSlice)
    offset: float = CO_LOCATED

    # ── geometry of the base circle ─────────────────────────────────────────
    @property
    def sigma(self) -> np.ndarray:
        """The drawn parameterisation, with ``σ = ±π`` both present."""
        return self.slice_.sigma

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
        d_a = np.abs(sigma)
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
        """``(r_A, r_B)`` before wrapping — outward and inward."""
        eps = self.slice_.gain if gain is None else float(gain)
        u_a, u_b = self.fields(closed=closed)
        mid = self.slice_.bulk.r_mid
        return mid + eps * u_a, mid - eps * u_b

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

        ``|δ| = gap`` with ``δ = ε(u_A + u_B)``, so this is
        ``gap / max|u_A + u_B|`` at the current time.
        """
        u_a, u_b = self.fields(closed=True)
        peak = float(np.max(np.abs(u_a + u_b)))
        return float("inf") if peak <= 0.0 else self.gap / peak

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
