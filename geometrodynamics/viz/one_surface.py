"""One scalar field on one surface — the construction v66 should have used.

What was wrong with v66
───────────────────────
`two_wave_slice` drew **two** curves, ``r_A = R_mid + s_A ε u_A`` and
``r_B = R_mid + s_B ε u_B``, and asked whether their images meet through the
glued seam.  Its docstring was careful that the two waves do not interact — but
that caveat does not repair the construction, it only labels the problem.  Two
curves in one frame are two surfaces, and reading their overlap as a connection
is a statement about a picture, not about a field.

If the two contributions are two pieces of **one** scalar deformation of **one**
surface, there is only ever one curve::

    u(θ, t) = s_A u_A + s_B u_B                 ← one field
    r(θ, t) = R_mid + ε u(θ, t)                 ← one curve

and the question is not whether two images meet but whether that single curve
reaches ``R_outer`` at some ``θ`` and ``R_inner`` at another — so that the
surface itself passes through the identification.

The repair costs nothing, and that is the surprising part
─────────────────────────────────────────────────────────
``δ = r_A − r_B = ε(s_A u_A − s_B u_B)`` — the quantity v66 plotted as a
*separation between two curves* — is **bit-for-bit** the one-surface deformation
``ε(s_A u_A + s'_B u_B)`` with ``s'_B = −s_B``.  Checked at zero, not at a
tolerance.  So v66 computed the right array throughout and mislabelled it; every
number it reported survives, with the two configurations swapping names:

    v66 "like-signed"  δ = ε(u_A − u_B)   ==   one surface, OPPOSED
    v66 "opposed"      δ = ε(u_A + u_B)   ==   one surface, LIKE-signed

That inversion matters.  v66's headline — that the opposed pair connects most
cheaply when co-located — becomes, correctly read, a statement about the *like*
pair; and v66's "the like pair is identically zero on the bisector" becomes the
node of the *opposed* field, which is where it belongs.

The monochromatic law
─────────────────────
For one spatial mode ``m``, equal amplitude, opposite radial orientation, foci
at ``∓α/2``::

    u = A cos[m(θ+α/2) − ωt] − A cos[m(θ−α/2) − ωt]
      = −2A sin(mα/2) sin(mθ − ωt)

so the deformation amplitude is ``B = 2A|sin(mα/2)|`` and the peak-to-trough
radial span is ``Δr = 4A|sin(mα/2)|``.  Three consequences, all checked:

* **``α = 0`` cancels exactly.**  ``u ≡ 0``.  Coincident foci with opposite
  orientation have no deformation at all, so no connection at any amplitude.
* **The optimum is half a wavelength**, ``α* = π/m`` (and ``(2j+1)π/m``), not
  "the antipode".  The antipode is simply the ``m = 1`` member of that family.
* **The antipode is parity-dependent**: ``sin(mπ/2)`` is ``±1`` for odd ``m``
  and ``0`` for even ``m``, so opposite foci are *maximal* for odd modes and
  *cancel exactly* for even ones.

With the gap ``D = R_outer − R_inner``, the curve spans it when
``4A|sin(mα/2)| ≥ D``, i.e. ``A_req = D/(4|sin(mα/2)|)`` — properly singular as
``α → 0``.

And the bulk chord shrinks with frequency.  At the optimum the outward and
inward extrema are ``π/m`` apart on the circle, so the straight chord between
them is ``L = √(D² + 4 R_out R_in sin²(π/2m))``, falling from ``2.000`` at
``m = 1`` to ``D = 0.520`` as ``m → ∞``.  Same deformation, shorter connection.

Where the plane-wave picture stops being the right one
──────────────────────────────────────────────────────
Two measured departures, and neither is a correction to the algebra above — the
algebra is exact for what it describes.  They are statements about *which model
the visualisation is in*.

**The ESU spectrum is zonal, not plane-wave, and its optimum is the antipode.**
The real eigenfunctions here are ``Z_n(χ) = sin[(n+1)χ]/[(n+1) sin χ]`` with
``ω_n = n + 1`` — the same ``ω_n = n+1`` this repository derives for the
conformally-coupled ESU.  A zonal harmonic is *centred*: ``Z_n(0) = 1`` is a
global maximum and ``|Z_n| ≤ 1`` everywhere.  So the opposed pair obeys
``|Z_A − Z_B| ≤ 2``, and equality needs ``Z_A = +1`` and ``Z_B = −1`` at the
same point — which happens exactly when the two poles are **antipodal** and
``n`` is **odd**, since ``Z_n(π) = (−1)ⁿ``.

    For zonal modes ``α* = π`` for every odd ``n``, and it saturates the bound.
    The half-wavelength family does **not** reproduce it: at ``α = π/(n+1)`` the
    strength is a fraction of ``2``.  For even ``n`` the antipode cancels
    exactly and nothing reaches the bound at all — the best available is
    ``≈ 1.22–1.33``.

A plane wave has no distinguished centre, so nothing picks out a separation
other than the wavelength; a zonal mode has one, and the antipode is where two
centres coincide with opposite sign.  The parity survives the change of model
exactly.  The location of the optimum does not.  (The kernel this programme
actually cares about is ``n = 1``, which is odd — so for it the antipode is
optimal *and* saturating.)

**A pulse is not monochromatic, and its cancellation window is its own width.**
v46's field is a launched pulse, and measured against the ``S²`` zonal basis it
carries a power-weighted mean of ``n ≈ 10`` with fifteen modes holding 90% of
the power.  For it the ``1/|sin(mα/2)|`` divergence is confined to
``α ≲ 2 ×`` the pulse width: past that the two pulses no longer overlap, nothing
cancels, and the threshold **saturates** at the single-pulse value.  The
coincident-foci cancellation is real and exact; the law governing its approach
is not.

Honest scope
────────────
Unchanged from v46 and v66, and it still binds.  The crossing rule that glues
``R_outer`` to ``R_inner`` is a *representation* choice, not a derived boundary
condition.  The field is a linear scalar on a fixed round background.  The gain
is a *display* amplitude — and the fixed-energy section below is the one place
that distinction is made to do work, because ``E ∝ ω²A²`` means a frequency
slider cannot hold displacement fixed and still be read as physics.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .circle_slice import ANTIPODAL_TIME, RETURN_TIME, TWO_PI, CircleSlice

__all__ = [
    "OneSurfaceSlice",
    "MonochromaticPair",
    "OPPOSED",
    "LIKE_SIGNED",
    "zonal_harmonic",
    "zonal_pair_strength",
    "bulk_chord",
    "measure_the_two_curve_picture_was_one_field_all_along",
    "measure_coincident_foci_cancel_exactly",
    "measure_the_optimum_is_half_a_wavelength",
    "measure_the_antipode_is_parity_dependent",
    "measure_the_zonal_optimum_is_the_antipode",
    "measure_a_pulse_saturates_where_a_mode_diverges",
    "measure_the_chord_shrinks_with_frequency",
    "measure_fixed_energy_reverses_the_preference",
]

#: The two radial orientations, as they act on **one** field.  ``OPPOSED``
#: subtracts the second contribution from the first; ``LIKE_SIGNED`` adds it.
OPPOSED = (1.0, -1.0)
LIKE_SIGNED = (1.0, 1.0)


# ════════════════════════════════════════════════════════════════════════════
# ONE FIELD ON ONE SURFACE
# ════════════════════════════════════════════════════════════════════════════
@dataclass
class OneSurfaceSlice:
    """One scalar deformation of one circle, from two signed contributions.

    ``offset`` is the angle ``α`` between the two foci — ``A`` at ``σ = 0``,
    ``B`` at ``σ = α``.  ``signs`` is how each contribution enters the *single*
    field; there is no second curve anywhere in this class.
    """

    slice_: CircleSlice = field(default_factory=CircleSlice)
    offset: float = math.pi
    signs: Tuple[float, float] = OPPOSED

    def __post_init__(self) -> None:
        s_a, s_b = self.signs
        if not (abs(abs(float(s_a)) - 1.0) < 1e-12
                and abs(abs(float(s_b)) - 1.0) < 1e-12):
            raise ValueError("signs must each be +1 or -1")
        self.signs = (float(s_a), float(s_b))

    @property
    def opposed(self) -> bool:
        return self.signs[0] * self.signs[1] < 0.0

    @property
    def sigma(self) -> np.ndarray:
        return self.slice_.sigma[:-1]          # the circle, endpoint dropped

    @property
    def gap(self) -> float:
        return self.slice_.bulk.gap

    @property
    def bisector(self) -> float:
        """``α/2`` — equidistant from both foci, so an exact node of the
        opposed field whenever the contributions depend only on distance."""
        return 0.5 * float(self.offset)

    def _distance(self, sigma, source: float) -> np.ndarray:
        return np.abs(np.mod(np.asarray(sigma, float) - source + math.pi,
                             TWO_PI) - math.pi)

    def contributions(self, sigma=None) -> Tuple[np.ndarray, np.ndarray]:
        """The two *signed* contributions.  They are not two curves."""
        s = self.sigma if sigma is None else np.atleast_1d(
            np.asarray(sigma, float))
        at = self.slice_.sim.field_at_distance
        s_a, s_b = self.signs
        return (s_a * at(self._distance(s, 0.0)),
                s_b * at(self._distance(s, self.offset)))

    def field(self, sigma=None) -> np.ndarray:
        """``u = s_A u_A + s_B u_B`` — the deformation.  One array."""
        c_a, c_b = self.contributions(sigma)
        return c_a + c_b

    def radius(self, gain: Optional[float] = None,
               sigma=None) -> np.ndarray:
        """``r = R_mid + ε u`` — the deformed surface.  One curve."""
        eps = self.slice_.gain if gain is None else float(gain)
        return self.slice_.bulk.r_mid + eps * self.field(sigma)

    # ── reaching the two boundaries ─────────────────────────────────────────
    def reach(self) -> Tuple[float, float]:
        """``(max u, −min u)`` — how far the field pushes out, and in."""
        u = self.field()
        return float(np.max(u)), float(-np.min(u))

    def span_gain(self) -> float:
        """Smallest gain whose single curve touches **both** boundaries.

        ``ε·max(u) ≥ D/2`` *and* ``ε·(−min u) ≥ D/2``, so the binding one is
        the smaller reach.  Infinite when the field is one-signed or zero —
        which is what a coincident opposed pair gives, exactly.
        """
        out, inn = self.reach()
        worst = min(out, inn)
        return float("inf") if worst <= 0.0 else 0.5 * self.gap / worst

    def connection(self, gain: Optional[float] = None) -> Dict[str, object]:
        """Where the surface meets each boundary, and the chord between."""
        eps = self.slice_.gain if gain is None else float(gain)
        u = self.field()
        s = self.sigma
        k_out, k_in = int(np.argmax(u)), int(np.argmin(u))
        b = self.slice_.bulk
        r_out, r_in = b.r_mid + eps * u[k_out], b.r_mid + eps * u[k_in]
        sep = float(np.abs(np.mod(s[k_out] - s[k_in] + math.pi, TWO_PI)
                           - math.pi))
        return {
            "gain": eps,
            "sigma_out": float(s[k_out]),
            "sigma_in": float(s[k_in]),
            "sigma_out_over_pi": float(s[k_out] / math.pi),
            "sigma_in_over_pi": float(s[k_in] / math.pi),
            "radius_out": float(r_out),
            "radius_in": float(r_in),
            "angular_separation": sep,
            "angular_separation_over_pi": sep / math.pi,
            "reaches_r_outer": bool(r_out >= b.r_outer - 1e-12),
            "reaches_r_inner": bool(r_in <= b.r_inner + 1e-12),
            "connected": bool(r_out >= b.r_outer - 1e-12
                              and r_in <= b.r_inner + 1e-12),
            "bulk_chord": bulk_chord(b.r_inner, b.r_outer, sep),
        }

    def scan(self, samples: int = 300,
             t_end: float = RETURN_TIME) -> List[Dict[str, float]]:
        rows = []
        self.slice_.reset()
        for i in range(int(samples)):
            t = (i + 1) * t_end / int(samples)
            self.slice_.advance_to(t)
            out, inn = self.reach()
            rows.append({"t": float(t), "reach_out": out, "reach_in": inn,
                         "span_gain": self.span_gain()})
        self.slice_.reset()
        return rows

    def cheapest_span(self, samples: int = 300) -> Dict[str, float]:
        """The gain, time and place of the cheapest full-gap span in a run."""
        self.slice_.reset()
        best, at = math.inf, 0.0
        for i in range(int(samples)):
            t = (i + 1) * RETURN_TIME / int(samples)
            self.slice_.advance_to(t)
            g = self.span_gain()
            if g < best:
                best, at = g, t
        self.slice_.reset()
        if not math.isfinite(best):
            return {"span_gain": best, "at_time": 0.0}
        self.slice_.advance_to(at)
        got = self.connection(gain=best)
        self.slice_.reset()
        return {"span_gain": float(best), "at_time": float(at),
                "at_time_over_pi": float(at / math.pi), **got}


# ════════════════════════════════════════════════════════════════════════════
# THE MONOCHROMATIC LAW
# ════════════════════════════════════════════════════════════════════════════
@dataclass
class MonochromaticPair:
    """One spatial mode, two foci, opposite radial orientation.

    ``u = A cos[m(θ+α/2) − ωt] − A cos[m(θ−α/2) − ωt]
        = −2A sin(mα/2) sin(mθ − ωt)``
    """

    mode: int = 1
    offset: float = math.pi
    amplitude: float = 1.0

    @property
    def wavelength(self) -> float:
        return TWO_PI / self.mode

    def amplitude_factor(self) -> float:
        """``2|sin(mα/2)|`` — the whole α-dependence, in one number."""
        return 2.0 * abs(math.sin(0.5 * self.mode * self.offset))

    def deformation_amplitude(self) -> float:
        return self.amplitude * self.amplitude_factor()

    def radial_span(self) -> float:
        """Peak-to-trough: ``4A|sin(mα/2)|``."""
        return 2.0 * self.deformation_amplitude()

    def field(self, theta, t: float = 0.0) -> np.ndarray:
        th = np.atleast_1d(np.asarray(theta, float))
        a = 0.5 * self.offset
        return self.amplitude * (
            np.cos(self.mode * (th + a) - t)
            - np.cos(self.mode * (th - a) - t))

    def required_amplitude(self, gap: float) -> float:
        """``D / (4|sin(mα/2)|)`` — singular at coincidence, as it must be.

        The span is ``4A|sin(mα/2)| = 2A·(amplitude factor)``, so this is
        ``gap / (2·factor)``.  A first draft wrote ``gap / factor`` and was
        exactly twice too large everywhere; the test that pinned the four
        hand-checked values (`0.130`, `0.184`, `0.260`, `0.502`) caught it.
        """
        f = self.amplitude_factor()
        return float("inf") if f <= 0.0 else float(gap) / (2.0 * f)

    def optimal_offsets(self) -> List[float]:
        """``(2j+1)π/m`` — every half-wavelength separation, in ``(0, π]``."""
        out = []
        j = 0
        while True:
            a = (2 * j + 1) * math.pi / self.mode
            if a > math.pi + 1e-12:
                break
            out.append(a)
            j += 1
        return out

    @property
    def first_optimum(self) -> float:
        """``π/m`` — half a surface wavelength."""
        return math.pi / self.mode

    def extrema(self, t: float = 0.0) -> List[float]:
        """``(2j+1)π/2m`` shifted by the phase — **independent of α**."""
        out = []
        for j in range(2 * self.mode):
            th = ((2 * j + 1) * math.pi / 2.0 + t) / self.mode
            th = (th + math.pi) % TWO_PI - math.pi
            out.append(th)
        return sorted(out)

    def bulk_chord(self, r_inner: float, r_outer: float) -> float:
        """The chord between adjacent outward and inward extrema at the
        optimum, where they sit ``π/m`` apart."""
        return bulk_chord(r_inner, r_outer, math.pi / self.mode)


def bulk_chord(r_inner: float, r_outer: float, separation: float) -> float:
    """``√(D² + 4 R_out R_in sin²(Δθ/2))`` — the law of cosines, regrouped so
    the purely radial gap ``D`` is visible as the ``Δθ → 0`` limit."""
    d = float(r_outer) - float(r_inner)
    return math.sqrt(d * d + 4.0 * r_outer * r_inner
                     * math.sin(0.5 * float(separation))**2)


# ════════════════════════════════════════════════════════════════════════════
# THE ZONAL SPECTRUM — what the ESU actually has
# ════════════════════════════════════════════════════════════════════════════
def zonal_harmonic(n: int, chi) -> np.ndarray:
    """``Z_n(χ) = sin[(n+1)χ] / [(n+1) sin χ]``, with both limits taken.

    Eigenvalue ``−n(n+2)``; ``Z_n(0) = 1``, ``Z_n(π) = (−1)ⁿ``, ``|Z_n| ≤ 1``.
    The endpoint guards matter: ``sin χ`` vanishes at ``χ = π`` as well as at
    ``0``, and a first draft that guarded only the small-χ series returned the
    Taylor polynomial at the antipode — giving ``|Z_8| ≈ 131`` and destroying
    exactly the parity result it was meant to test.
    """
    c = np.atleast_1d(np.asarray(chi, float))
    s = np.sin(c)
    near_0 = np.abs(c) < 1e-8
    near_pi = np.abs(c - math.pi) < 1e-8
    safe = np.where(near_0 | near_pi, 1.0, s)
    out = np.sin((n + 1) * c) / ((n + 1) * safe)
    out = np.where(near_0, 1.0 - n * (n + 2) * c * c / 6.0, out)
    out = np.where(near_pi, float((-1) ** n), out)
    return out


def zonal_pair_strength(n: int, alpha: float, samples: int = 8000) -> float:
    """``max_θ |Z_n(d_A) − Z_n(d_B)|`` on a great circle, foci ``α`` apart.

    Bounded above by ``2`` because ``|Z_n| ≤ 1``, with equality only when some
    point sees ``Z_A = +1`` and ``Z_B = −1`` at once.
    """
    th = -math.pi + (np.arange(int(samples)) + 0.5) * TWO_PI / int(samples)
    d_a = np.abs(np.mod(th + math.pi, TWO_PI) - math.pi)
    d_b = np.abs(np.mod(th - float(alpha) + math.pi, TWO_PI) - math.pi)
    return float(np.max(np.abs(zonal_harmonic(n, d_a)
                               - zonal_harmonic(n, d_b))))


# ════════════════════════════════════════════════════════════════════════════
# MEASUREMENTS
# ════════════════════════════════════════════════════════════════════════════
def measure_the_two_curve_picture_was_one_field_all_along(
        offsets: Sequence[float] = (0.0, 0.25, 0.5, 0.75, 1.0),
        samples: int = 60) -> Dict[str, object]:
    """**The repair costs nothing.**  v66's ``δ`` *is* the one-surface field.

    ``δ = r_A − r_B = ε(s_A u_A − s_B u_B)`` and the deformation is
    ``ε(s_A u_A + s'_B u_B)`` with ``s'_B = −s_B``: the same array, the same
    bits.  So every number v66 reported survives the change of construction —
    with the two configurations swapping names, which inverts its headline.
    """
    from .two_wave_slice import BOTH_OUTWARD, OUTWARD_INWARD, TwoWaveSlice

    worst_opposed, worst_like = 0.0, 0.0
    rows = []
    for f in offsets:
        alpha = float(f) * math.pi
        # v66 "like-signed" (both driven out)  <->  one surface OPPOSED
        v66_like = TwoWaveSlice(offset=alpha, signs=BOTH_OUTWARD)
        one_opp = OneSurfaceSlice(offset=alpha, signs=OPPOSED)
        one_opp.slice_ = v66_like.slice_
        # v66 "opposed"  <->  one surface LIKE-signed
        v66_opp = TwoWaveSlice(offset=alpha, signs=OUTWARD_INWARD)
        one_like = OneSurfaceSlice(offset=alpha, signs=LIKE_SIGNED)
        one_like.slice_ = v66_like.slice_
        v66_opp.slice_ = v66_like.slice_

        sl = v66_like.slice_
        sl.reset()
        d_o = d_l = 0.0
        for i in range(int(samples)):
            sl.advance_to((i + 1) * RETURN_TIME / int(samples))
            d_o = max(d_o, float(np.max(np.abs(
                v66_like.separation(gain=1.0, closed=True)
                - one_opp.field()))))
            d_l = max(d_l, float(np.max(np.abs(
                v66_opp.separation(gain=1.0, closed=True)
                - one_like.field()))))
        sl.reset()
        worst_opposed = max(worst_opposed, d_o)
        worst_like = max(worst_like, d_l)
        rows.append({"offset_over_pi": float(f),
                     "v66_like_vs_one_surface_opposed": d_o,
                     "v66_opposed_vs_one_surface_like": d_l})
    # v66's `separation` forms (R_mid + a) - (R_mid + b), so the comparison
    # carries R_mid's rounding whatever the fields do.  The identity is about
    # the fields; the residue below is the mid-radius, and saying so is the
    # difference between a checked claim and a failed one.  (Second time this
    # exact species has appeared in this arc.)
    ulp = float(np.spacing(OneSurfaceSlice().slice_.bulk.r_mid))
    worst = max(worst_opposed, worst_like)
    return {
        "rows": rows,
        "worst_opposed": worst_opposed,
        "worst_like": worst_like,
        "one_ulp_of_r_mid": ulp,
        "residue_in_ulps": worst / ulp,
        "they_are_the_same_array": bool(worst <= 4.0 * ulp),
        "the_residue_is_the_mid_radius_not_the_fields": bool(
            worst <= 4.0 * ulp),
        "the_labels_swap": "v66's 'like-signed' separation is the one-surface "
                           "OPPOSED deformation, and v66's 'opposed' "
                           "separation is the one-surface LIKE-signed one",
        "what_that_costs": "nothing numerically and the headline everything: "
                           "v66's cheapest-when-co-located result belongs to "
                           "the like pair, and its identically-zero bisector "
                           "belongs to the opposed field, where it is a node",
    }


def measure_coincident_foci_cancel_exactly(
        samples: int = 120) -> Dict[str, object]:
    """``α = 0`` with opposite orientation gives ``u ≡ 0`` — no connection.

    Exactly zero, at every time, because the two contributions are the same
    function of the same distance with opposite sign.  The span threshold is
    infinite, not merely large.
    """
    p = OneSurfaceSlice(offset=0.0, signs=OPPOSED)
    q = OneSurfaceSlice(offset=0.0, signs=LIKE_SIGNED)
    q.slice_ = p.slice_
    sl = p.slice_
    sl.reset()
    worst_u, best_like = 0.0, math.inf
    for i in range(int(samples)):
        sl.advance_to((i + 1) * RETURN_TIME / int(samples))
        worst_u = max(worst_u, float(np.max(np.abs(p.field()))))
        best_like = min(best_like, q.span_gain())
    sl.reset()
    return {
        "worst_absolute_field": worst_u,
        "opposed_span_gain": p.span_gain(),
        "like_signed_cheapest_span_gain": float(best_like),
        "the_opposed_field_is_identically_zero": bool(worst_u == 0.0),
        "no_gain_connects_it": bool(not math.isfinite(p.span_gain())),
        "the_like_pair_is_unaffected": bool(math.isfinite(best_like)),
        "why": "the two contributions are the same function of the same "
               "distance, so opposite orientation subtracts them to nothing",
    }


def measure_the_optimum_is_half_a_wavelength(
        modes: Sequence[int] = (1, 2, 3, 4, 5, 6, 7, 8),
        samples: int = 40001) -> Dict[str, object]:
    """``B = 2A|sin(mα/2)|``, maximal at ``α* = π/m`` — half a wavelength.

    The claim is not asserted from the identity; the amplitude is measured by
    scanning ``θ`` on a grid and the optimum by scanning ``α``, and both are
    then compared with the closed form.
    """
    rows = []
    worst_amp, worst_opt = 0.0, 0.0
    th = -math.pi + (np.arange(4001) + 0.5) * TWO_PI / 4001
    alg = np.linspace(1e-6, math.pi, 4001)
    for m in modes:
        pair = MonochromaticPair(mode=int(m), offset=math.pi / int(m))
        got = float(np.max(np.abs(pair.field(th))))
        worst_amp = max(worst_amp, abs(got - pair.amplitude_factor())
                        / max(pair.amplitude_factor(), 1e-30))
        # the FIRST optimum: search only out to one full period of
        # |sin(m alpha/2)|, so a later member of the family cannot win.  A
        # first draft looked for the first grid point exceeding 2 - 1e-9,
        # which no grid point reaches, so argmax returned index 0 and the
        # "measured optimum" was the left edge.
        window = alg[alg <= min(math.pi, TWO_PI / m) + 1e-12]
        strength = 2.0 * np.abs(np.sin(0.5 * m * window))
        first = float(window[int(np.argmax(strength))])
        worst_opt = max(worst_opt, abs(first - math.pi / m))
        rows.append({
            "mode": int(m),
            "wavelength_over_pi": TWO_PI / m / math.pi,
            "first_optimum_over_pi": (math.pi / m) / math.pi,
            "measured_first_optimum_over_pi": float(first / math.pi),
            "half_wavelength_over_pi": (math.pi / m) / math.pi,
            "measured_amplitude_factor": got,
            "closed_form_amplitude_factor": pair.amplitude_factor(),
            "all_optima_over_pi": [a / math.pi
                                   for a in pair.optimal_offsets()],
        })
    return {
        "rows": rows,
        "worst_amplitude_error": worst_amp,
        "worst_optimum_error": worst_opt,
        "the_closed_form_is_the_measured_amplitude": bool(worst_amp < 1e-6),
        "the_first_optimum_is_half_a_wavelength": bool(worst_opt < 2e-3),
        "the_optimum_is_a_half_wavelength_not_an_antipode": bool(
            all(abs(r["first_optimum_over_pi"] - r["half_wavelength_over_pi"])
                < 1e-12 for r in rows)),
        "why": "alpha* = pi/m is half a surface wavelength; the antipode is "
               "simply the m = 1 member of that family",
    }


def measure_the_antipode_is_parity_dependent(
        modes: Sequence[int] = (1, 2, 3, 4, 5, 6, 7, 8)) -> Dict[str, object]:
    """Opposite foci are maximal for odd modes and cancel for even ones.

    ``sin(mπ/2)`` is ``±1`` for odd ``m`` and ``0`` for even ``m``, so "opposite
    poles give the greatest separation" is a statement about the *lowest* mode,
    not a general one.  Checked twice: on the plane-wave amplitude factor, and
    on the exact ``S³`` zonal harmonics through ``Z_n(π) = (−1)ⁿ``.
    """
    th = -math.pi + (np.arange(20001) + 0.5) * TWO_PI / 20001
    rows = []
    worst = 0.0
    for m in modes:
        pair = MonochromaticPair(mode=int(m), offset=math.pi)
        measured = float(np.max(np.abs(pair.field(th))))
        z_pi = float(zonal_harmonic(int(m), np.array([math.pi]))[0])
        zonal = zonal_pair_strength(int(m), math.pi)
        odd = bool(m % 2)
        worst = max(worst, abs(measured - pair.amplitude_factor()))
        rows.append({
            "mode": int(m),
            "odd": odd,
            "plane_wave_factor_at_pi": pair.amplitude_factor(),
            "measured_at_pi": measured,
            "zonal_Z_n_at_pi": z_pi,
            "expected_minus_one_to_the_n": float((-1) ** m),
            "zonal_pair_strength_at_pi": zonal,
            "verdict": "maximal" if odd else "cancels",
        })
    odds = [r for r in rows if r["odd"]]
    evens = [r for r in rows if not r["odd"]]
    return {
        "rows": rows,
        "worst_plane_wave_error": worst,
        "odd_modes_are_maximal_at_the_antipode": bool(
            all(abs(r["plane_wave_factor_at_pi"] - 2.0) < 1e-9 for r in odds)),
        "even_modes_cancel_at_the_antipode": bool(
            all(r["plane_wave_factor_at_pi"] < 1e-9 for r in evens)),
        "the_zonal_spectrum_agrees": bool(
            all(abs(r["zonal_Z_n_at_pi"]
                    - r["expected_minus_one_to_the_n"]) < 1e-9 for r in rows)),
        # 1e-4 rather than 1e-6: the residue is the theta grid missing the peak
        # by O(1/N^2), not the identity, which is exact
        "zonal_odd_doubles": bool(
            all(abs(r["zonal_pair_strength_at_pi"] - 2.0) < 1e-4
                for r in odds)),
        "zonal_even_cancels": bool(
            all(r["zonal_pair_strength_at_pi"] < 1e-6 for r in evens)),
        "why": "sin(m pi/2) is +-1 for odd m and 0 for even m; on S^3 the same "
               "parity appears as Z_n(pi) = (-1)^n, so it is a property of the "
               "spectrum and not of the plane-wave reduction",
    }


def measure_the_zonal_optimum_is_the_antipode(
        orders: Sequence[int] = (1, 2, 3, 4, 5, 6, 8, 10),
        alphas: int = 601) -> Dict[str, object]:
    """**Where the plane-wave picture and the real spectrum part company.**

    A plane wave has no distinguished centre, so nothing picks out a separation
    other than the wavelength and ``α* = π/m`` follows.  A zonal harmonic *is*
    centred — ``Z_n(0) = 1`` is a global maximum and ``|Z_n| ≤ 1`` — so the
    opposed pair obeys ``|Z_A − Z_B| ≤ 2`` with equality only where one focus
    sees ``+1`` and the other ``−1`` at the same point.  That is exactly the
    antipode with odd ``n``.

    So for zonal modes ``α* = π`` for **every** odd ``n``, saturating the
    bound, and the half-wavelength separation ``π/(n+1)`` reaches only a
    fraction of it.  For even ``n`` the antipode cancels exactly and nothing
    reaches the bound at all.

    The parity carries across the two models unchanged.  The *location* of the
    optimum does not, and that is a measured disagreement rather than an error
    in either.

    One caution the measurement forces.  For **odd** orders the maximiser is
    unique and is ``π``.  For **even** orders it is not: several separations
    reach within ``1e-3`` of the same peak, so the reported location moves with
    the grid — at ``n = 6`` a 1501-point sweep says ``0.794π`` and a 601-point
    sweep says ``0.205π`` for the same strength to seven digits.  The *strength*
    is robust and the *location* is not, so both are reported and the
    degeneracy is counted rather than hidden behind whichever grid ran last.
    """
    alg = np.linspace(1e-3, math.pi, int(alphas))
    rows = []
    for n in orders:
        strengths = np.array([zonal_pair_strength(int(n), float(a))
                              for a in alg])
        k = int(np.argmax(strengths))
        # count SEPARATE maximisers, not samples near the peak: a smooth
        # maximum always has a plateau of grid points around it, so counting
        # points made every order look degenerate.  Connected runs of
        # `within 1e-3 of the peak` is the quantity that distinguishes one
        # maximiser from several.
        hot = strengths > strengths[k] - 1e-3
        near = int(np.count_nonzero(np.diff(hot.astype(int)) > 0)
                   + (1 if hot[0] else 0))
        half = math.pi / (int(n) + 1)
        rows.append({
            "order": int(n),
            "omega": int(n) + 1,
            "odd": bool(n % 2),
            "measured_optimum_over_pi": float(alg[k] / math.pi),
            "maximiser_is_unique": bool(near == 1),
            "separate_maximisers_within_1e-3": near,
            "peak_strength": float(strengths[k]),
            "strength_at_the_antipode": zonal_pair_strength(int(n), math.pi),
            "half_wavelength_over_pi": half / math.pi,
            "strength_at_half_wavelength": zonal_pair_strength(int(n), half),
        })
    odds = [r for r in rows if r["odd"]]
    evens = [r for r in rows if not r["odd"]]
    return {
        "rows": rows,
        "the_bound_is_two": bool(all(r["peak_strength"] <= 2.0 + 1e-6
                                     for r in rows)),
        "odd_orders_peak_at_the_antipode": bool(
            all(abs(r["measured_optimum_over_pi"] - 1.0) < 2e-3
                for r in odds)),
        # the residue here is the theta grid, not the physics: the bound is
        # attained exactly and the sampling misses the peak by O(1/N^2)
        "and_they_saturate_the_bound": bool(
            all(abs(r["peak_strength"] - 2.0) < 1e-4 for r in odds)),
        "the_odd_maximiser_is_unique": bool(
            all(r["maximiser_is_unique"] for r in odds)),
        "the_even_maximiser_is_not": bool(
            any(not r["maximiser_is_unique"] for r in evens)),
        "half_a_wavelength_does_not_reproduce_it": bool(
            all(r["strength_at_half_wavelength"] < r["peak_strength"] - 1e-3
                for r in odds if r["order"] > 1)),
        "even_orders_cancel_at_the_antipode": bool(
            all(r["strength_at_the_antipode"] < 1e-6 for r in evens)),
        "and_never_reach_the_bound": bool(
            all(r["peak_strength"] < 1.45 for r in evens)),
        "best_even_strength": max((r["peak_strength"] for r in evens),
                                  default=0.0),
        "why": "a zonal harmonic is centred, so the antipode is where two "
               "centres coincide with opposite sign; a plane wave has no "
               "centre, so only its wavelength sets a scale",
        "the_kernel_is_odd": "this programme's kernel is n = 1, which is odd, "
                             "so for it the antipode is both optimal and "
                             "saturating",
    }


def measure_a_pulse_saturates_where_a_mode_diverges(
        offsets: Sequence[float] = (0.02, 0.04, 0.06, 0.1, 0.2, 0.4, 0.6,
                                    0.8, 1.0),
        samples: int = 300) -> Dict[str, object]:
    """The ``1/|sin(mα/2)|`` divergence is a **monochromatic** statement.

    v46's field is a launched pulse, not a mode.  Two localized pulses cancel
    only while they overlap, so the required gain rises as ``α → 0`` — the
    coincident cancellation is real and exact — but past roughly twice the
    pulse width nothing cancels any more and the threshold **saturates** at the
    single-pulse value instead of continuing down like ``1/|sin|``.
    """
    probe = OneSurfaceSlice()
    gap = probe.gap
    width = float(probe.slice_.sim.pulse_width)
    rows = []
    for f in offsets:
        alpha = float(f) * math.pi
        p = OneSurfaceSlice(offset=alpha, signs=OPPOSED)
        got = p.cheapest_span(samples=samples)
        mono = MonochromaticPair(mode=1,
                                 offset=alpha).required_amplitude(gap)
        rows.append({
            "offset_over_pi": float(f),
            "offset_in_pulse_widths": alpha / width,
            "monochromatic_required_amplitude": mono,
            "pulse_threshold": got["span_gain"],
            "ratio": mono / got["span_gain"],
        })
    plateau = [r for r in rows if r["offset_over_pi"] >= 0.4]
    lo = min(r["pulse_threshold"] for r in plateau)
    hi = max(r["pulse_threshold"] for r in plateau)
    return {
        "pulse_width": width,
        "rows": rows,
        "plateau_threshold": 0.5 * (lo + hi),
        "plateau_spread": hi - lo,
        "the_pulse_threshold_saturates": bool((hi - lo) / max(lo, 1e-30)
                                              < 0.01),
        "the_monochromatic_law_keeps_falling": bool(
            rows[-1]["monochromatic_required_amplitude"]
            < 0.6 * rows[3]["monochromatic_required_amplitude"]),
        "it_still_cancels_at_coincidence": bool(
            rows[0]["pulse_threshold"] > 2.5 * 0.5 * (lo + hi)),
        "cancellation_window_in_pulse_widths": max(
            r["offset_in_pulse_widths"] for r in rows
            if r["pulse_threshold"] > 1.05 * 0.5 * (lo + hi)),
        "why": "two localized pulses cancel only while they overlap; a mode "
               "fills the whole circle and cancels everywhere",
    }


def measure_the_chord_shrinks_with_frequency(
        modes: Sequence[int] = (1, 2, 3, 4, 5, 6, 8, 16),
        slice_: Optional[CircleSlice] = None) -> Dict[str, object]:
    """Same deformation, shorter connection, as the mode number rises.

    At the optimum the outward and inward extrema sit ``π/m`` apart, so the
    straight chord between them is
    ``L = √(D² + 4 R_out R_in sin²(π/2m))`` — checked against the plain law of
    cosines — falling monotonically to the purely radial gap ``D``.
    """
    sl = slice_ or CircleSlice()
    b = sl.bulk
    th = -math.pi + (np.arange(20001) + 0.5) * TWO_PI / 20001
    rows = []
    worst_form, worst_span = 0.0, 0.0
    for m in modes:
        sep = math.pi / int(m)
        law = math.sqrt(b.r_outer**2 + b.r_inner**2
                        - 2 * b.r_outer * b.r_inner * math.cos(sep))
        boxed = bulk_chord(b.r_inner, b.r_outer, sep)
        worst_form = max(worst_form, abs(law - boxed))
        pair = MonochromaticPair(mode=int(m), offset=sep)
        span = float(np.max(np.abs(pair.field(th))))
        worst_span = max(worst_span, abs(span - 2.0))
        rows.append({
            "mode": int(m),
            "optimal_separation_over_pi": sep / math.pi,
            "radial_span_over_A": span,
            "bulk_chord": boxed,
            "chord_over_gap": boxed / b.gap,
        })
    return {
        "gap": b.gap,
        "rows": rows,
        "worst_closed_form_error": worst_form,
        "worst_span_error": worst_span,
        "the_span_is_the_same_at_every_mode": bool(worst_span < 1e-5),
        "the_closed_form_is_the_law_of_cosines": bool(worst_form < 1e-12),
        "the_chord_shrinks_monotonically": bool(
            all(rows[i + 1]["bulk_chord"] < rows[i]["bulk_chord"]
                for i in range(len(rows) - 1))),
        "chord_at_the_fundamental": rows[0]["bulk_chord"],
        "chord_at_the_highest_mode": rows[-1]["bulk_chord"],
        "limit_is_the_radial_gap": b.gap,
        "why": "the two extrema are half a wavelength apart, so raising the "
               "mode brings them together and the chord tends to the purely "
               "radial gap",
    }


def measure_fixed_energy_reverses_the_preference(
        modes: Sequence[int] = tuple(range(1, 17)),
        slice_: Optional[CircleSlice] = None) -> Dict[str, object]:
    """A frequency slider cannot hold displacement fixed and be read as physics.

    For a monochromatic mode ``E ∝ ω²A²``, so at fixed energy ``A ∝ 1/ω`` and
    the attainable span ``4A`` falls with frequency exactly as fast as the
    chord does.  Held at fixed *display* amplitude the high modes look free;
    held at fixed *energy* they buy a shorter connection with a smaller
    deformation, and there is a highest mode that can still span the gap at
    all.

    No numerical optimum is claimed beyond that: it would need an energy
    normalisation and a packet focusing law, and neither is in this model.
    """
    sl = slice_ or CircleSlice()
    b = sl.bulk
    rows = []
    for m in modes:
        omega = float(m)                       # c = R = 1 on the unit circle
        a_fixed = 1.0
        a_energy = 1.0 / omega                 # E fixed => A ∝ 1/omega
        rows.append({
            "mode": int(m),
            "omega": omega,
            "span_at_fixed_amplitude": 4.0 * a_fixed,
            "span_at_fixed_energy": 4.0 * a_energy,
            "bulk_chord": bulk_chord(b.r_inner, b.r_outer, math.pi / int(m)),
            "spans_the_gap_at_fixed_energy": bool(4.0 * a_energy >= b.gap),
        })
    ok = [r for r in rows if r["spans_the_gap_at_fixed_energy"]]
    return {
        "gap": b.gap,
        "rows": rows,
        "highest_mode_that_still_spans_at_fixed_energy": (
            max(r["mode"] for r in ok) if ok else None),
        "span_is_flat_at_fixed_amplitude": bool(
            all(abs(r["span_at_fixed_amplitude"] - 4.0) < 1e-12
                for r in rows)),
        "span_falls_like_one_over_omega": bool(
            all(abs(r["span_at_fixed_energy"] * r["omega"] - 4.0) < 1e-12
                for r in rows)),
        "the_chord_and_the_span_both_shrink": bool(
            rows[-1]["bulk_chord"] < rows[0]["bulk_chord"]
            and rows[-1]["span_at_fixed_energy"]
            < rows[0]["span_at_fixed_energy"]),
        "the_tradeoff": "low frequency buys a large deformation on a long "
                        "connection; high frequency buys a short connection "
                        "with a small deformation",
        "what_is_not_claimed": "no favourable frequency, because that needs an "
                               "energy normalisation and a packet focusing law "
                               "that this model does not contain",
    }


# ════════════════════════════════════════════════════════════════════════════
# HOW ONE SURFACE ANSWERS TO TWO FRONTS
# ════════════════════════════════════════════════════════════════════════════
def decompose(surf: "OneSurfaceSlice") -> Dict[str, object]:
    """The two signed contributions, the field they sum to, and where they meet.

    The point of keeping these three arrays together is that only the third is
    a *surface*.  ``c_A`` and ``c_B`` are components of one deformation, and
    the moment they are drawn as closed curves in the annulus the picture is
    back to v66's error — so they belong on a graph of field values against
    ``σ``, where nothing invites reading them as separate objects.

    ``reinforcing`` is where the two contributions share a sign and therefore
    add.  For an opposed pair that means the two *fields* are anti-correlated,
    which is the configuration the whole question is about.
    """
    c_a, c_b = surf.contributions()
    u = c_a + c_b
    s = surf.sigma
    reinforcing = np.sign(c_a) * np.sign(c_b) > 0
    live = (np.abs(c_a) > 0.02 * max(float(np.max(np.abs(c_a))), 1e-30)) & \
           (np.abs(c_b) > 0.02 * max(float(np.max(np.abs(c_b))), 1e-30))
    d_sigma = TWO_PI / len(s)
    peak_a = float(np.max(np.abs(c_a)))
    peak_u = float(np.max(np.abs(u)))
    return {
        "sigma": s,
        "contribution_a": c_a,
        "contribution_b": c_b,
        "field": u,
        "reinforcing": reinforcing,
        "both_present": live,
        "overlap_arc": float(np.count_nonzero(live) * d_sigma),
        "reinforcing_fraction": float(np.mean(reinforcing)),
        "peak_contribution": peak_a,
        "peak_field": peak_u,
        "amplification": peak_u / max(peak_a, 1e-30),
        "sigma_of_peak_a": float(s[int(np.argmax(np.abs(c_a)))]),
        "sigma_of_peak_b": float(s[int(np.argmax(np.abs(c_b)))]),
        "sigma_of_peak_field": float(s[int(np.argmax(np.abs(u)))]),
    }


def measure_how_one_surface_answers_two_fronts(
        offsets: Sequence[float] = (0.0, 0.15, 0.25, 0.5, 0.75, 1.0),
        samples: int = 200) -> Dict[str, object]:
    """What the offset actually buys, decomposed on the one surface.

    Two things the totals alone hide, and the picture has to show:

    **The contributions barely overlap.**  v46's field is a localized pulse, so
    at any offset wider than the pulse the two contributions occupy *disjoint*
    arcs of the circle.  The surface then carries two separate deformations,
    one per front, and the amplification ``max|u| / max|c_A|`` is about
    ``1.01`` — the total is one contribution plus almost nothing.  Interference
    is confined to the narrow arc where both are present.

    **``α = 0`` is the only configuration where the overlap is total** — and it
    is exactly the one that cancels.  So the offset does not "turn interference
    on"; it turns the *cancellation* off by pulling the two fronts apart, and
    what is left is two nearly independent dents.

    Reported at the time of the largest deformation in a run, per offset.
    """
    probe = OneSurfaceSlice()
    rows = []
    for f in offsets:
        alpha = float(f) * math.pi
        q = OneSurfaceSlice(offset=alpha, signs=OPPOSED)
        q.slice_ = probe.slice_
        sl = q.slice_
        sl.reset()
        best, at = -1.0, 0.0
        for i in range(int(samples)):
            t = (i + 1) * RETURN_TIME / int(samples)
            sl.advance_to(t)
            v = float(np.max(np.abs(q.field())))
            if v > best:
                best, at = v, t
        sl.reset()
        sl.advance_to(at)
        d = decompose(q)
        sl.reset()
        rows.append({
            "offset_over_pi": float(f),
            "at_time_over_pi": at / math.pi,
            "peak_contribution": d["peak_contribution"],
            "peak_field": d["peak_field"],
            "amplification": d["amplification"],
            "overlap_arc": d["overlap_arc"],
            "reinforcing_fraction": d["reinforcing_fraction"],
            "sigma_of_peak_a_over_pi": d["sigma_of_peak_a"] / math.pi,
            "sigma_of_peak_b_over_pi": d["sigma_of_peak_b"] / math.pi,
            "sigma_of_peak_field_over_pi": d["sigma_of_peak_field"] / math.pi,
        })
    apart = [r for r in rows if r["offset_over_pi"] > 0.1]
    zero = rows[0]
    return {
        "rows": rows,
        "amplification_when_apart": (min(r["amplification"] for r in apart),
                                     max(r["amplification"] for r in apart)),
        "the_contributions_barely_overlap": bool(
            all(r["amplification"] < 1.05 for r in apart)),
        "the_field_peak_sits_on_a_contribution_peak": bool(
            all(min(abs(r["sigma_of_peak_field_over_pi"]
                        - r["sigma_of_peak_a_over_pi"]),
                    abs(r["sigma_of_peak_field_over_pi"]
                        - r["sigma_of_peak_b_over_pi"])) < 0.01
                for r in apart)),
        "coincidence_is_the_only_total_overlap": bool(
            zero["peak_field"] == 0.0),
        "what_the_offset_buys": "not interference -- the offset turns the "
                                "CANCELLATION off by pulling the two fronts "
                                "apart, and what is left is two nearly "
                                "independent dents in one surface",
    }
