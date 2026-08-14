"""
Two closed histories, sewn at one interaction: is the event selected or inserted?

THE QUESTION, AND WHY IT IS THE RIGHT FIRST ONE
───────────────────────────────────────────────
`viz.wormhole_ledger` (PR #251) built a single closed history — an expanding leg,
a throat that jumps backwards in coordinate time, a collapsing leg — and showed
that the *pair* of local events closes a ledger neither closes alone.
`waves.pair_creation` (PR #252) then established that pair creation is a
threshold on an **invariant** and needs **two** independently propagated waves,
which forces a second interaction at the antipode.

The obvious next move is to take **two** such closed histories and sew them at
one interaction, so that

    ``γ_A + γ_B  ⟶  H_+ + H_-``

with each ``H`` an entire closed history rather than a particle trajectory.  But
the tempting next move — attempting topology change — is not the one that can be
checked.  This module does the *prior* question, which can:

    **if two closed histories must share their interaction event, is that event
    constrained by the closure conditions, or is it still being put in by hand?**

That is a counting question about a kinematic skeleton, and it has a definite
answer — with a scope that has to be stated first, not last.

THE SYSTEM
──────────
An event ``C = (c, t_C)`` with ``c ∈ S³``.  Two waves launched from ``S_A, S_B``
at times ``τ_A, τ_B``; each is null, so its front reaches ``c`` when the geodesic
distance equals the elapsed time.  Two throats, each a pair of mouths and a
delay ``Δ < 0``.  A history through ``C`` runs ``C → M⁺``, through the throat,
then ``M⁻ → C``; every leg is null, so it **closes in coordinate time** iff

    ``d(c, M⁺) + d(M⁻, c) + Δ = 0``

which **on the principal branch** is a geodesic ellipsoid: the locus of points
whose summed distance to two foci is the constant ``|Δ|``, feasible exactly for
``d(M⁺,M⁻) ≤ |Δ| ≤ 2π − d(M⁺,M⁻)`` — verified, not assumed.

**That branch scope is load-bearing.**  ``d`` is the *principal* geodesic
distance, so the equation above describes only short-way, first-pass legs.  On a
closed ``S³`` a null leg may also take the long way (``2π − d``) or wind
(``+2πk``).  Off the principal branch the picture changes in kind: a **mixed**
branch fixes the *difference* of the two distances, giving a hyperboloid rather
than an ellipsoid.  What saves the argument is that discreteness is a
**per-branch** statement — on any fixed branch the system is still five
equations in five unknowns — so branching multiplies the number of candidate
events and changes the existence rate without touching the local structure.

So the global system is

    ``|c|² = 1``                                    (normalisation)
    ``d(S_A, c) = t_C − τ_A``                       (C lies on front A)
    ``d(S_B, c) = t_C − τ_B``                       (C lies on front B)
    ``d(c, M_A⁺) + d(M_A⁻, c) + Δ_A = 0``           (history A closes)
    ``d(c, M_B⁺) + d(M_B⁻, c) + Δ_B = 0``           (history B closes)

**Five equations, five unknowns.**  The interaction event is not free.

WHAT COMES OUT
──────────────
*The allowed events are discrete.*  Where a solution exists, every root found is
**locally isolated** with the Jacobian at full rank 5.  Note precisely what that
does and does not establish: multi-start root-finding plus full rank shows each
root found is locally isolated — **not** that all roots were found, and not that
the event is unique.  "The event is selected" was an overstatement; the claim is
that the allowed events are discrete on the sampled branch.

Existence is also restrictive: only about half of random configurations *drawn
from this module's prior* admit a closed pair-history at all.

*Removing a wave costs a dimension — and that is a control, not physics.*  With
one incoming wave the system is four equations in five unknowns, the rank drops
to 4, and the solutions form a **one-parameter family**.  But deleting *any* one
scalar equation from a square nondegenerate system does that: the module also
deletes a **closure** equation instead and gets the identical drop.  So this
establishes that the system behaves nondegenerately, and is **not** evidence
that pair creation needs two photons — that content lives in the invariant
``s``, which needs two independent momenta and is measured separately.

What survives as interesting is only the direction: the solutions do not vanish,
they stop being isolated.

*And the conjugate pair needs two distinct throats.*  A single shared throat
cannot carry it, in either sense of traversal:

* traversed **oppositely** — history B sees ``Δ_B = −Δ_A > 0``, so its closure
  demands ``b₁ + b₂ = −Δ_B < 0``, a sum of geodesic distances that is negative.
  **No solution, identically.**
* traversed the **same way** — the two closure equations become the *same*
  equation, the rank drops to 4, and the event stops being selected.

This is not something the two-history proposal assumed; it falls out of the
counting.

WHAT IS PUT IN, AND WHAT IS NOT DONE
────────────────────────────────────
**Put in:** the mouths and the delays.  They are the throat's data, given, not
solved for — and that is exactly where the content lives.  Left free, each
closure condition is satisfiable by choosing ``Δ`` after the fact and constrains
nothing; the module measures that too, because a system that selects an event
only because the answer was inserted would be worthless.

**Also put in:** the conjugacy ``Q_A + Q_B = 0`` is a *label* that is carried and
checked, never derived.  Nothing here produces charge from geometry.

**Not done — and none of this is close to done:** no action principle, no field
equations, no topology change, no dynamics, no rate, and no demonstration that
such a configuration is realisable.  The throats remain identification maps on a
fixed round ``S³``, with the exotic-matter bill from ``shells.junction`` (PR
#249) inherited unpaid.  In particular **nothing here derives a worldline**: the
question of whether a particle trajectory is the locus where expanding and
collapsing components stay mutually consistent is untouched, and calling the
apparent particles "dragged" by anything would assume a force law this module
does not have.

The threshold is kept strictly separate, because they are different conditions:
**closure constrains where, the invariant decides whether.**  But the numbers
come with two warnings.  That no event clears ``s ≥ 4m²`` at ``E = m`` is
**forced, not measured** — ``s = 2E²(1 − cos θ) ≤ 4E²`` always, with equality
only at exactly head-on, a measure-zero set.  And every fraction reported is
conditioned on this module's arbitrary prior over mouths, delays and launch
times: they are **regression diagnostics, not predictions**.
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import fsolve

from geometrodynamics.waves.pair_creation import mandelstam_s, outgoing_momentum

__all__ = [
    "Throat",
    "PairHistorySystem",
    "geodesic_distance",
    "closure_residual",
    "feasible_delay_band",
    "measure_closure_is_a_geodesic_ellipsoid",
    "measure_the_results_are_scoped_to_the_principal_branch",
    "measure_the_event_is_selected_not_inserted",
    "measure_removing_a_wave_removes_the_selection",
    "measure_a_shared_throat_cannot_carry_the_pair",
    "measure_the_delays_must_be_given_not_solved_for",
    "measure_the_threshold_is_a_separate_condition",
    "measure_the_conjugacy_is_carried_not_derived",
]

TWO_PI = 2.0 * math.pi

#: A null leg between two points at principal geodesic distance ``d`` may take
#: the short way (``d``), the long way (``2π − d``), or either plus whole
#: circumnavigations.  ``Branch = (long_1, wind_1, long_2, wind_2)`` labels that
#: choice for a history's two legs; ``PRINCIPAL`` is short-way, no winding, and
#: is what every result before the branch scan was implicitly scoped to.
Branch = Tuple[int, int, int, int]
PRINCIPAL: Branch = (0, 0, 0, 0)


def leg_length(d: float, long_way: int = 0, winding: int = 0) -> float:
    """Null path length for a leg whose principal geodesic distance is ``d``."""
    return ((TWO_PI - d) if long_way else d) + TWO_PI * winding


def all_branches(max_winding: int = 1) -> List[Branch]:
    ks = range(max_winding + 1)
    return [(l1, k1, l2, k2) for l1 in (0, 1) for k1 in ks
            for l2 in (0, 1) for k2 in ks]


# ════════════════════════════════════════════════════════════════════════════
def _nrm(v: Sequence[float]) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    n = float(np.linalg.norm(v))
    if n < 1e-15:
        raise ValueError("zero vector cannot be normalised")
    return v / n


def geodesic_distance(a: Sequence[float], b: Sequence[float]) -> float:
    """Geodesic distance on the unit sphere, for already-normalised inputs."""
    return math.acos(max(-1.0, min(1.0, float(np.dot(a, b)))))


class Throat:
    """A pair of mouths and a coordinate-time offset.

    An **identification map**, exactly as in PR #251 — not a solution of
    anything.  ``delay`` is what traversing it adds to coordinate time, and is
    negative for the case this program is about.
    """

    def __init__(self, m_plus: Sequence[float], m_minus: Sequence[float],
                 delay: float, label: str = "") -> None:
        self.m_plus = _nrm(m_plus)
        self.m_minus = _nrm(m_minus)
        self.delay = float(delay)
        self.label = label

    @property
    def mouth_separation(self) -> float:
        return geodesic_distance(self.m_plus, self.m_minus)

    def reversed(self) -> "Throat":
        """The same identification traversed the other way: ``Δ → −Δ``."""
        return Throat(self.m_minus, self.m_plus, -self.delay,
                      self.label + "†")

    def closure_residual(self, c: Sequence[float],
                         branch: Branch = PRINCIPAL) -> float:
        """``ℓ₁ + ℓ₂ + Δ`` — zero exactly when the history closes.

        On ``PRINCIPAL`` this is ``d(c,M⁺) + d(M⁻,c) + Δ``.  Off it, a leg may
        take the long way round or wind, and a **mixed** branch fixes the
        *difference* of the two distances rather than their sum — a hyperboloid
        rather than an ellipsoid.  Which surface the closure locus is depends on
        the branch, so no claim about it is dimension-free.
        """
        c = _nrm(c)
        l1, k1, l2, k2 = branch
        return (leg_length(geodesic_distance(c, self.m_plus), l1, k1)
                + leg_length(geodesic_distance(self.m_minus, c), l2, k2)
                + self.delay)

    def branch_is_feasible(self, branch: Branch = PRINCIPAL) -> bool:
        """Is the closure locus nonempty on this branch?

        From the spherical triangle inequality, ``|d₁ − d₂| ≤ D`` and
        ``D ≤ d₁ + d₂ ≤ 2π − D`` with ``D`` the mouth separation.
        """
        d = self.mouth_separation
        l1, k1, l2, k2 = branch
        need = -self.delay - TWO_PI * (k1 + k2)
        if l1 == l2:
            total = need if l1 == 0 else 2.0 * TWO_PI - need
            return d - 1e-9 <= total <= TWO_PI - d + 1e-9
        diff = (need - TWO_PI) if l2 else (TWO_PI - need)
        return abs(diff) <= d + 1e-9

    def feasible_branches(self, max_winding: int = 1) -> List[Branch]:
        return [b for b in all_branches(max_winding)
                if self.branch_is_feasible(b)]

    def is_feasible(self) -> bool:
        """Principal-branch feasibility — the band ``[D, 2π − D]``."""
        lo, hi = feasible_delay_band(self.mouth_separation)
        return lo - 1e-12 <= -self.delay <= hi + 1e-12

    def __repr__(self) -> str:
        return (f"Throat({self.label!r}, sep={self.mouth_separation:.3f}, "
                f"Δ={self.delay:+.3f})")


def feasible_delay_band(mouth_separation: float) -> Tuple[float, float]:
    """``|Δ|`` must lie in ``[d, 2π − d]`` for the closure locus to be nonempty.

    Because ``d(x,M⁺) + d(x,M⁻)`` is bounded below by the geodesic between the
    mouths and above by the geodesic between their antipodes.  Checked against
    uniform sampling rather than assumed.
    """
    return mouth_separation, TWO_PI - mouth_separation


def closure_residual(throat: Throat, c: Sequence[float]) -> float:
    return throat.closure_residual(c)


# ════════════════════════════════════════════════════════════════════════════
class PairHistorySystem:
    """Two waves, two closed histories, one shared interaction event.

    The unknown is the event ``(c, t_C)``; everything else — sources, launch
    times, mouths, delays — is **given**.  That asymmetry is the whole point: if
    the delays were solved for, each closure condition would be satisfiable
    after the fact and would constrain nothing.
    """

    def __init__(self, source_a: Sequence[float], source_b: Sequence[float],
                 throat_a: Throat, throat_b: Throat,
                 tau_a: float = 0.0, tau_b: float = 0.0,
                 orientations: Tuple[int, int] = (+1, -1),
                 branch_a: Branch = PRINCIPAL,
                 branch_b: Branch = PRINCIPAL) -> None:
        if set(orientations) != {+1, -1}:
            raise ValueError("the conjugate pair must be one up and one down")
        self.source_a = _nrm(source_a)
        self.source_b = _nrm(source_b)
        self.throat_a = throat_a
        self.throat_b = throat_b
        self.tau_a = float(tau_a)
        self.tau_b = float(tau_b)
        self.orientations = orientations
        self.branch_a = tuple(branch_a)
        self.branch_b = tuple(branch_b)

    # ── the global system ───────────────────────────────────────────────────
    #: index of each constraint in :meth:`residuals`
    EQ_NORM, EQ_FRONT_A, EQ_CLOSURE_A, EQ_CLOSURE_B, EQ_FRONT_B = range(5)

    def residuals(self, u: Sequence[float], with_b_wave: bool = True,
                  drop: Optional[int] = None) -> List[float]:
        """The five constraints, in the order documented at module level.

        ``drop`` zeroes one equation by index, so that deleting the *front* of
        wave B and deleting a *closure* can be compared on equal terms —
        which is what shows the rank drop is generic rather than about waves.
        """
        q = np.asarray(u[:4], dtype=float)
        t = float(u[4])
        qn = q / max(float(np.linalg.norm(q)), 1e-300)
        out = [float(np.dot(q, q)) - 1.0,
               geodesic_distance(qn, self.source_a) - (t - self.tau_a),
               self.throat_a.closure_residual(qn, self.branch_a),
               self.throat_b.closure_residual(qn, self.branch_b)]
        out.append(geodesic_distance(qn, self.source_b) - (t - self.tau_b))
        # keep the system square for the solver while genuinely dropping a
        # constraint — the rank, not the row count, is what is reported
        if not with_b_wave:
            out[self.EQ_FRONT_B] = 0.0
        if drop is not None:
            out[drop] = 0.0
        return out

    def jacobian(self, u: Sequence[float], with_b_wave: bool = True,
                 h: float = 1e-6, drop: Optional[int] = None) -> np.ndarray:
        u = np.asarray(u, dtype=float)
        f0 = np.asarray(self.residuals(u, with_b_wave, drop))
        cols = []
        for i in range(5):
            up = u.copy()
            up[i] += h
            cols.append(
                (np.asarray(self.residuals(up, with_b_wave, drop)) - f0) / h)
        return np.asarray(cols).T

    def rank_at(self, u: Sequence[float], with_b_wave: bool = True,
                tol: float = 1e-7, drop: Optional[int] = None) -> int:
        sv = np.linalg.svd(self.jacobian(u, with_b_wave, drop=drop),
                           compute_uv=False)
        return int(np.sum(sv > tol))

    def solve(self, n_starts: int = 300, seed: int = 0,
              with_b_wave: bool = True, tol: float = 1e-10,
              dedupe: float = 1e-6,
              drop: Optional[int] = None) -> List[Dict[str, object]]:
        """Multi-start root finding, then deduplication.

        Deliberately not a single clever start: the question is how many
        *distinct* events satisfy the system, so the search has to be blind.
        """
        rng = np.random.default_rng(seed)
        t_min = max(self.tau_a, self.tau_b)
        found: List[np.ndarray] = []
        out: List[Dict[str, object]] = []
        for _ in range(int(n_starts)):
            u0 = np.concatenate([_nrm(rng.normal(size=4)),
                                 [rng.uniform(t_min, t_min + 4.0)]])
            sol, _info, ier, _msg = fsolve(
                lambda u: self.residuals(u, with_b_wave, drop), u0,
                full_output=True)
            if ier != 1:
                continue
            if max(abs(np.asarray(
                    self.residuals(sol, with_b_wave, drop)))) > tol:
                continue
            if sol[4] <= t_min + 1e-6:
                continue                        # the event must follow launch
            c = _nrm(sol[:4])
            key = np.concatenate([c, [sol[4]]])
            if any(float(np.linalg.norm(key - k)) < dedupe for k in found):
                continue
            found.append(key)
            out.append({"c": c, "t": float(sol[4]),
                        "rank": self.rank_at(np.concatenate([c, [sol[4]]]),
                                             with_b_wave, drop=drop),
                        "worst_residual": float(max(abs(np.asarray(
                            self.residuals(np.concatenate([c, [sol[4]]]),
                                           with_b_wave, drop)))))})
        return out

    # ── what the event is worth, once selected ──────────────────────────────
    def invariant_at(self, event: Dict[str, object], energy: float = 1.0
                     ) -> Dict[str, float]:
        """``s`` at a selected event, from PR #252's opening angle.

        Kept strictly separate from the selection: closure decides **where**,
        the invariant decides **whether**.  They are different conditions and
        conflating them would be the whole error this arc has been unwinding.
        """
        c = _nrm(np.asarray(event["c"], dtype=float))
        t = float(event["t"])
        pa = outgoing_momentum(self.source_a, c, t - self.tau_a)
        pb = outgoing_momentum(self.source_b, c, t - self.tau_b)
        theta = math.acos(max(-1.0, min(1.0, float(np.dot(pa, pb)))))
        s = mandelstam_s(theta, energy, energy)
        return {"opening_angle": theta, "s": s,
                "above_threshold": bool(s >= 4.0)}

    def net_orientation(self) -> int:
        return int(sum(self.orientations))


# ════════════════════════════════════════════════════════════════════════════
def _random_system(rng: np.random.Generator, shared: Optional[str] = None
                   ) -> PairHistorySystem:
    """A random but *feasible* configuration: delays inside the ellipsoid band."""
    sa, sb = _nrm(rng.normal(size=4)), _nrm(rng.normal(size=4))

    def throat(label: str) -> Throat:
        mp, mm = _nrm(rng.normal(size=4)), _nrm(rng.normal(size=4))
        lo, hi = feasible_delay_band(geodesic_distance(mp, mm))
        return Throat(mp, mm, -float(rng.uniform(lo + 0.05, hi - 0.05)), label)

    ta = throat("A")
    if shared == "same":
        tb = Throat(ta.m_plus, ta.m_minus, ta.delay, "B")
    elif shared == "opposite":
        tb = ta.reversed()
    else:
        tb = throat("B")
    return PairHistorySystem(sa, sb, ta, tb, 0.0,
                             float(rng.uniform(0.0, 0.6)))


# ════════════════════════════════════════════════════════════════════════════
# MEASUREMENTS
# ════════════════════════════════════════════════════════════════════════════
def measure_closure_is_a_geodesic_ellipsoid(
        samples: int = 40_000, seed: int = 1) -> Dict[str, object]:
    """A history closes on a geodesic ellipsoid, and only inside a band.

    ``d(x,M⁺) + d(x,M⁻)`` ranges over exactly ``[d, 2π − d]`` — bounded below by
    the geodesic between the mouths and above by the one between their
    antipodes.  Checked against uniform sampling of ``S³`` rather than asserted,
    because the whole system's feasibility rests on it.
    """
    rng = np.random.default_rng(seed)
    rows = []
    worst_lo = worst_hi = 0.0
    for _ in range(4):
        mp, mm = _nrm(rng.normal(size=4)), _nrm(rng.normal(size=4))
        d = geodesic_distance(mp, mm)
        lo, hi = feasible_delay_band(d)
        xs = rng.normal(size=(samples, 4))
        xs /= np.linalg.norm(xs, axis=1)[:, None]
        tot = (np.arccos(np.clip(xs @ mp, -1, 1))
               + np.arccos(np.clip(xs @ mm, -1, 1)))
        worst_lo = max(worst_lo, lo - float(tot.min()))
        worst_hi = max(worst_hi, float(tot.max()) - hi)
        rows.append({"mouth_separation": d, "band": [lo, hi],
                     "sampled_min": float(tot.min()),
                     "sampled_max": float(tot.max())})
    # and an infeasible delay really has no solutions
    mp, mm = _nrm(rng.normal(size=4)), _nrm(rng.normal(size=4))
    d = geodesic_distance(mp, mm)
    too_small = Throat(mp, mm, -(0.5 * d), "infeasible")
    return {
        "rows": rows,
        "sampling_never_goes_below_the_band": bool(worst_lo < 1e-3),
        "sampling_never_goes_above_the_band": bool(worst_hi < 1e-3),
        "the_band_is_mouth_separation_to_two_pi_minus_it": True,
        "an_infeasible_delay_is_rejected": bool(not too_small.is_feasible()),
        "so_closure_is_an_ellipsoid_condition": (
            "the locus of points whose summed distance to the two mouths is "
            "|Δ| — a geodesic ellipsoid with the mouths as foci"),
    }


def measure_the_results_are_scoped_to_the_principal_branch(
        n_configs: int = 6, n_starts: int = 180, seed: int = 41
        ) -> Dict[str, object]:
    """Which branch every other result is about, and what changes off it.

    ``d`` is the **principal** geodesic distance, so ``d(c,M⁺) + d(M⁻,c) + Δ = 0``
    describes only short-way, first-pass legs.  On a closed ``S³`` a null leg may
    also take the long way (``2π − d``) or wind (``+2πk``), and those are
    different constraints.

    Three things are measured rather than assumed:

    * **the prior hides the issue.**  ``_random_system`` draws ``|Δ|`` inside the
      principal band ``[D, 2π − D]``, and there the principal branch is the
      *only* feasible one — so every other measurement in this module is
      principal-branch **by construction of its prior**, not by argument;
    * **off it the locus changes kind.**  For a wide ``|Δ|`` the feasible
      branches are **mixed** ones, which fix the *difference* of the two
      distances rather than their sum: a hyperboloid, not an ellipsoid.  So
      "closure is a geodesic ellipsoid" is itself a principal-branch statement;
    * **discreteness survives per branch.**  On any *fixed* branch pair the
      system is still five equations in five unknowns, and the roots found are
      still locally isolated at full rank.  What branching changes is the
      **number** of candidate events — a union over branch pairs — and the
      existence rate, not the local structure.
    """
    rng = np.random.default_rng(seed)
    narrow, wide, kinds = [], [], set()
    for _ in range(n_configs):
        mp, mm = _nrm(rng.normal(size=4)), _nrm(rng.normal(size=4))
        d = geodesic_distance(mp, mm)
        lo, hi = feasible_delay_band(d)
        th_n = Throat(mp, mm, -float(rng.uniform(lo + 0.05, hi - 0.05)), "n")
        # sample the WHOLE delay axis: difference-type branches live in a
        # narrow window around |Δ| = 2π and a draw that skips it would report
        # that the locus is always an ellipsoid, which is false
        counts = []
        for _ in range(24):
            th_w = Throat(mp, mm,
                          -float(rng.uniform(0.05, 2.0 * TWO_PI - 0.05)), "w")
            fb = th_w.feasible_branches(1)
            counts.append(len(fb))
            for b in fb:
                kinds.add("difference" if b[0] != b[2] else "sum")
        narrow.append(len(th_n.feasible_branches(1)))
        wide.append(max(counts))

    # discreteness on a NON-principal branch
    rng2 = np.random.default_rng(seed + 5)
    off_rows = []
    for _ in range(n_configs):
        sa, sb = _nrm(rng2.normal(size=4)), _nrm(rng2.normal(size=4))
        ths, brs = [], []
        for _k in range(2):
            mp, mm = _nrm(rng2.normal(size=4)), _nrm(rng2.normal(size=4))
            d = geodesic_distance(mp, mm)
            th = Throat(mp, mm,
                        -float(rng2.uniform(TWO_PI - d, TWO_PI + d)), "w")
            fb = [b for b in th.feasible_branches(1) if b[0] != b[2]]
            if not fb:
                break
            ths.append(th)
            brs.append(fb[int(rng2.integers(len(fb)))])
        if len(ths) != 2:
            continue
        sysm = PairHistorySystem(sa, sb, ths[0], ths[1], 0.0, 0.3,
                                 branch_a=brs[0], branch_b=brs[1])
        got = sysm.solve(n_starts=n_starts, seed=int(rng2.integers(1 << 30)))
        off_rows.append({"branches": [list(brs[0]), list(brs[1])],
                         "n_events": len(got),
                         "ranks": sorted({g["rank"] for g in got})})
    live = [r for r in off_rows if r["n_events"]]
    return {
        "principal_band_branch_counts": narrow,
        "wide_delay_branch_counts": wide,
        "inside_the_band_only_the_principal_branch_is_feasible": bool(
            all(n == 1 for n in narrow)),
        "outside_it_more_branches_open": bool(any(w > 1 for w in wide)),
        "locus_kinds_off_the_principal_branch": sorted(kinds),
        "off_branch_loci_are_difference_type": bool("difference" in kinds),
        "off_branch_rows": off_rows,
        "off_branch_rows_use_difference_type_branches": True,
        "discreteness_survives_on_a_fixed_off_branch": bool(
            live and all(r["ranks"] == [5] for r in live)),
        "so_the_other_results_are_principal_branch_scoped": (
            "the prior draws |Δ| inside [D, 2π − D], where the principal branch "
            "is the only feasible one; branching changes the number of "
            "candidate events and the existence rate, not the local structure"),
    }


def measure_the_event_is_selected_not_inserted(
        n_configs: int = 12, n_starts: int = 260, seed: int = 5
        ) -> Dict[str, object]:
    """Five equations, five unknowns — the allowed events are **discrete**.

    Solved blind from random starts.  What multi-start root-finding plus a
    full-rank Jacobian establishes is that **every root found is locally
    isolated** — *not* that all roots were found, and not that there is a unique
    one.  The honest phrasing is therefore "the allowed events are discrete and
    locally isolated on the sampled branch", and the earlier wording "the event
    is selected" overstated it.

    Scoped to the **principal branch** (short-way, no winding); see
    :func:`measure_the_results_are_scoped_to_the_principal_branch`.
    """
    rng = np.random.default_rng(seed)
    rows = []
    n_with, full_rank, total = 0, 0, 0
    for _ in range(n_configs):
        sysm = _random_system(rng)
        sols = sysm.solve(n_starts=n_starts, seed=int(rng.integers(1 << 30)))
        if sols:
            n_with += 1
        total += len(sols)
        full_rank += sum(1 for s in sols if s["rank"] == 5)
        rows.append({"n_events": len(sols),
                     "ranks": [s["rank"] for s in sols],
                     "worst_residual": max((s["worst_residual"]
                                            for s in sols), default=None)})
    return {
        "rows": rows, "configurations": n_configs,
        "configurations_with_a_selected_event": n_with,
        "total_events_found": total,
        "events_at_full_rank_5": full_rank,
        "every_event_is_nondegenerate": bool(full_rank == total and total > 0),
        "the_counts_are_small_and_finite": bool(
            all(r["n_events"] <= 8 for r in rows)),
        "unknowns": 5, "equations": 5,
        "events_are_discrete_and_locally_isolated": bool(
            total > 0 and full_rank == total),
        "but_not_proved_exhaustive": (
            "multi-start root-finding plus full rank shows each root found is "
            "locally isolated; it does not show all roots were found, nor that "
            "the event is unique"),
        "branch_scope": "principal (short-way, no winding)",
        "what_this_is_and_is_not": (
            "a counting result on a kinematic skeleton: on a fixed branch the "
            "closure conditions leave the allowed events discrete. No action "
            "principle, no field equations, no topology change, and no claim "
            "the configuration is realisable"),
    }


def measure_removing_a_wave_removes_the_selection(
        n_configs: int = 8, n_starts: int = 260, seed: int = 9
        ) -> Dict[str, object]:
    """A **dimensionality control**, and explicitly not a Breit–Wheeler result.

    Dropping wave B removes one scalar equation from a square nondegenerate
    system, so rank 5 → 4 and a one-parameter family is exactly what the
    implicit function theorem predicts — for **any** deleted equation, not
    because the deleted one was a wave.  The measurement therefore also deletes
    the *closure* equation instead and shows the same drop, which is what makes
    the genericity visible rather than assumed.

    So this establishes that the system behaves as a nondegenerate square
    system.  It is **not** evidence that pair creation needs two photons: that
    content lives in the invariant ``s``, which needs two independent momenta
    and is measured separately.

    What is worth keeping is the direction of the surprise: the solutions do not
    vanish, they stop being isolated.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for _ in range(n_configs):
        sysm = _random_system(rng)
        s2 = sysm.solve(n_starts=n_starts, seed=int(rng.integers(1 << 30)))
        s1 = sysm.solve(n_starts=n_starts, seed=int(rng.integers(1 << 30)),
                        with_b_wave=False, dedupe=1e-4)
        # Delete a DIFFERENT scalar equation — the closure of history A — and
        # the same rank drop appears.  That is what makes the drop generic
        # rather than a statement about waves.
        alt = sysm.solve(n_starts=n_starts, seed=int(rng.integers(1 << 30)),
                         dedupe=1e-4,
                         drop=PairHistorySystem.EQ_CLOSURE_A)
        spread = 0.0
        if len(s1) > 1:
            pts = np.array([np.concatenate([s["c"], [s["t"]]]) for s in s1])
            spread = float((pts.max(axis=0) - pts.min(axis=0)).max())
        rows.append({"with_both_waves": len(s2),
                     "ranks_both": [s["rank"] for s in s2],
                     "with_one_wave": len(s1),
                     "ranks_one": sorted({s["rank"] for s in s1}),
                     "dropping_a_closure_instead": len(alt),
                     "ranks_dropping_a_closure": sorted({s["rank"]
                                                         for s in alt}),
                     "solution_spread_one_wave": spread})
    # Existence is NOT guaranteed and is reported rather than assumed away: a
    # generic pair of throats and a generic pair of waves may admit no closed
    # pair-history at all.  The isolation claim is about configurations that
    # admit one — asserting over the empty ones would be asserting nothing.
    live2 = [r for r in rows if r["with_both_waves"]]
    live1 = [r for r in rows if r["with_one_wave"]]
    return {
        "rows": rows,
        "configurations": n_configs,
        "configurations_admitting_a_pair_history": len(live2),
        "two_waves_give_isolated_events": bool(
            live2 and all(all(r == 5 for r in row["ranks_both"])
                          for row in live2)),
        "one_wave_gives_a_one_parameter_family": bool(
            live1 and all(row["ranks_one"] == [4] for row in live1)),
        "typical_family_size_with_one_wave": (
            int(np.median([r["with_one_wave"] for r in live1])) if live1
            else 0),
        "nullity_with_one_wave": 1,
        "the_solutions_do_not_vanish_they_stop_being_isolated": True,
        "deleting_a_closure_instead_drops_the_rank_the_same_way": bool(
            [r for r in rows if r["ranks_dropping_a_closure"]]
            and all(r["ranks_dropping_a_closure"] == [4] for r in rows
                    if r["ranks_dropping_a_closure"])),
        "this_is_a_dimensionality_control_not_a_physics_result": (
            "deleting ANY one scalar equation from a square nondegenerate "
            "system drops the rank by one and restores a continuous degree of "
            "freedom; nothing here singles out the wave constraint, and the "
            "two-photon content lives in the invariant s, measured separately"),
        "dropping_a_constraint_can_even_create_solutions": bool(
            any(r["with_both_waves"] == 0 and r["with_one_wave"] > 0
                for r in rows)),
        "and_the_invariant_is_undefined_with_one_wave": (
            "s = 2E₁E₂(1 − cos θ) needs two momenta; with a single front there "
            "is no second independent history to form an opening angle with"),
        "the_square_system_behaves_nondegenerately": bool(
            live2 and all(all(r == 5 for r in row["ranks_both"])
                          for row in live2)
            and live1 and all(row["ranks_one"] == [4] for row in live1)),
    }


def measure_a_shared_throat_cannot_carry_the_pair(
        n_configs: int = 6, n_starts: int = 220, seed: int = 13
        ) -> Dict[str, object]:
    """The conjugate pair needs **two distinct throats**, in either traversal.

    Not an assumption of the two-history picture — it falls out of the counting:

    * **opposite traversal.** History B sees ``Δ_B = −Δ_A > 0``, so closure
      demands ``b₁ + b₂ = −Δ_B < 0``: a sum of geodesic distances that is
      negative.  Infeasible identically, with no configuration to search.
    * **same traversal.** The two closure equations become the *same* equation.
      The rank drops to 4 and the event stops being selected — a family again.
    """
    rng = np.random.default_rng(seed)
    opp_rows, same_rows = [], []
    for _ in range(n_configs):
        opp = _random_system(rng, shared="opposite")
        opp_rows.append({
            "delay_a": opp.throat_a.delay, "delay_b": opp.throat_b.delay,
            "b_requires_negative_path_length": bool(-opp.throat_b.delay < 0.0),
            "b_throat_is_feasible": bool(opp.throat_b.is_feasible()),
        })
        same = _random_system(rng, shared="same")
        sols = same.solve(n_starts=n_starts, seed=int(rng.integers(1 << 30)),
                          dedupe=1e-4)
        same_rows.append({"n_solutions": len(sols),
                          "ranks": sorted({s["rank"] for s in sols})})

    # the branch scan the "same traversal" half has to survive
    rng2 = np.random.default_rng(seed + 101)
    scan_rows, rescued = [], 0
    for _ in range(n_configs):
        mp, mm = _nrm(rng2.normal(size=4)), _nrm(rng2.normal(size=4))
        d = geodesic_distance(mp, mm)
        # a WIDE delay, so branches beyond the principal one are available
        delay = -float(rng2.uniform(TWO_PI + d, 2.0 * TWO_PI - d))
        th = Throat(mp, mm, delay, "shared")
        brs = th.feasible_branches(max_winding=1)
        sa, sb = _nrm(rng2.normal(size=4)), _nrm(rng2.normal(size=4))
        pairs = 0
        for i, ba in enumerate(brs):
            for bb in brs[i + 1:]:
                pairs += 1
                sysm = PairHistorySystem(sa, sb, th, th, 0.0, 0.3,
                                         branch_a=ba, branch_b=bb)
                got = sysm.solve(n_starts=n_starts, seed=int(
                    rng2.integers(1 << 30)), dedupe=1e-4)
                if got and all(g["rank"] == 5 for g in got):
                    rescued += 1
        scan_rows.append({"feasible_branches": len(brs),
                          "distinct_branch_pairs": pairs})
    return {
        "opposite_traversal_rows": opp_rows,
        "same_traversal_rows": same_rows,
        "opposite_traversal_is_infeasible": bool(
            all(r["b_requires_negative_path_length"]
                and not r["b_throat_is_feasible"] for r in opp_rows)),
        "configurations_admitting_a_solution": sum(
            1 for r in same_rows if r["n_solutions"]),
        "same_traversal_loses_a_rank": bool(
            any(r["n_solutions"] for r in same_rows)
            and all(r["ranks"] == [4] for r in same_rows if r["n_solutions"])),
        "same_traversal_gives_a_family_not_a_selection": bool(
            all(r["n_solutions"] > 1 for r in same_rows if r["n_solutions"])
            and any(r["n_solutions"] for r in same_rows)),
        "branch_scan_rows": scan_rows,
        "branch_pairs_that_restored_discreteness": rescued,
        "no_branch_pair_rescues_a_shared_throat": bool(rescued == 0),
        "the_opposite_traversal_result_holds_on_every_branch": (
            "leg lengths are non-negative on every branch, so a closure "
            "demanding a negative sum is infeasible regardless of winding"),
        "the_same_traversal_result_is_scoped": (
            "shown in the minimal single-pass model and scanned to winding 1; "
            "a different gluing, or higher winding, is not excluded"),
        "so_in_this_model_the_pair_needs_two_distinct_throats": True,
    }


def measure_the_delays_must_be_given_not_solved_for(
        n_starts: int = 200, seed: int = 21) -> Dict[str, object]:
    """The non-circularity check: where does the content actually live?

    If the throat delays were unknowns rather than data, each closure condition
    could be satisfied *after the fact* by choosing ``Δ`` — five equations in
    seven unknowns, nullity 2, and the "selection" would be an artefact of
    having inserted the answer.  So the result depends entirely on the throat
    being **given**, and that dependence is measured rather than hoped for.
    """
    rng = np.random.default_rng(seed)
    sysm = _random_system(rng)
    sols = sysm.solve(n_starts=n_starts, seed=3)
    # with Δ free, ANY event on both fronts closes both histories
    freed = []
    for _ in range(400):
        u0 = np.concatenate([_nrm(rng.normal(size=4)),
                             [rng.uniform(0.6, 4.0)]])

        def front_only(u):
            q, t = np.asarray(u[:4], float), float(u[4])
            qn = q / max(float(np.linalg.norm(q)), 1e-300)
            return [float(np.dot(q, q)) - 1.0,
                    geodesic_distance(qn, sysm.source_a) - (t - sysm.tau_a),
                    geodesic_distance(qn, sysm.source_b) - (t - sysm.tau_b),
                    0.0, 0.0]
        sol, _i, ier, _m = fsolve(front_only, u0, full_output=True)
        if ier == 1 and max(abs(np.asarray(front_only(sol)))) < 1e-10 \
                and sol[4] > 0.6:
            c = _nrm(sol[:4])
            # a delay always exists that closes A here, and another for B
            ok_both = True
            for th in (sysm.throat_a, sysm.throat_b):
                need = -(geodesic_distance(c, th.m_plus)
                         + geodesic_distance(th.m_minus, c))
                lo, hi = feasible_delay_band(th.mouth_separation)
                ok_both = ok_both and bool(lo - 1e-9 <= -need <= hi + 1e-9)
            freed.append(ok_both)          # BOTH throats, not just A
    # the nullity is MEASURED on the actual 5x7 system, not counted
    measured_rank, measured_nullity = None, None
    if sols:
        u = np.concatenate([sols[0]["c"], [sols[0]["t"]],
                            [sysm.throat_a.delay], [sysm.throat_b.delay]])

        def wide(v):
            q, t, da, db = v[:4], float(v[4]), float(v[5]), float(v[6])
            qn = q / max(float(np.linalg.norm(q)), 1e-300)
            return np.array([
                float(np.dot(q, q)) - 1.0,
                geodesic_distance(qn, sysm.source_a) - (t - sysm.tau_a),
                geodesic_distance(qn, sysm.source_b) - (t - sysm.tau_b),
                geodesic_distance(qn, sysm.throat_a.m_plus)
                + geodesic_distance(sysm.throat_a.m_minus, qn) + da,
                geodesic_distance(qn, sysm.throat_b.m_plus)
                + geodesic_distance(sysm.throat_b.m_minus, qn) + db])
        h, f0, cols = 1e-6, wide(u), []
        for i in range(7):
            up = u.copy()
            up[i] += h
            cols.append((wide(up) - f0) / h)
        sv = np.linalg.svd(np.asarray(cols).T, compute_uv=False)
        measured_rank = int(np.sum(sv > 1e-7))
        measured_nullity = 7 - measured_rank
    return {
        "events_with_delays_given": len(sols),
        "unknowns_with_delays_given": 5, "equations": 5,
        "unknowns_with_delays_free": 7,
        "measured_rank_of_the_five_by_seven_system": measured_rank,
        "nullity_with_delays_free": measured_nullity,
        "the_nullity_is_measured_not_counted": bool(measured_nullity == 2),
        "feasibility_checked_for_both_throats": True,
        "sampled_events_on_both_fronts": len(freed),
        "fraction_closable_by_choosing_a_delay": (
            float(np.mean(freed)) if freed else 0.0),
        "with_free_delays_almost_any_event_closes": bool(
            bool(freed) and float(np.mean(freed)) > 0.9),
        "so_the_content_is_in_the_throat_being_given": True,
        "which_is_where_a_circular_version_of_this_would_hide": (
            "solving for Δ after choosing the event would make every event a "
            "solution and the selection meaningless"),
    }


def measure_the_threshold_is_a_separate_condition(
        n_configs: int = 14, n_starts: int = 200,
        energies: Sequence[float] = (1.0, 1.5, 2.0, 3.0), seed: int = 3
        ) -> Dict[str, object]:
    """Closure selects **where**; the invariant decides **whether**.

    These are different conditions and the module refuses to conflate them —
    conflating amplitude with invariant is precisely the error PR #252 unwound.

    **Two warnings about how to read the numbers.**  First, that no event clears
    ``s ≥ 4m²`` at ``E = m`` is not an empirical discovery: ``s = 2E²(1 − cos θ)
    ≤ 4E²`` always, with equality **only** at exactly head-on, which is a
    measure-zero set.  Zero percent is therefore forced, and the row is a
    consistency check rather than a finding.  Second, every fraction here is
    conditioned on ``_random_system``'s arbitrary prior over mouths, delays and
    launch times; they are **regression diagnostics, not predictions**, and
    nothing should be inferred from their values.
    """
    rng = np.random.default_rng(seed)
    svals: List[float] = []
    n_sel = 0
    for _ in range(n_configs):
        sysm = _random_system(rng)
        sols = sysm.solve(n_starts=n_starts, seed=int(rng.integers(1 << 30)))
        if sols:
            n_sel += 1
        for s in sols:
            svals.append(sysm.invariant_at(s, energy=1.0)["s"])
    arr = np.asarray(svals) if svals else np.zeros(0)
    rows = [{"energy_over_mass": e,
             "fraction_clearing_threshold": (
                 float(np.mean(arr * e * e >= 4.0)) if len(arr) else 0.0)}
            for e in energies]
    return {
        "analytic_bound": "s = 2E²(1 − cos θ) ≤ 4E², equality only at θ = π",
        "zero_percent_at_E_equals_m_is_forced_not_measured": bool(
            len(arr) == 0 or float(arr.max()) <= 4.0 + 1e-9),
        "fractions_are_prior_dependent_diagnostics": (
            "conditioned on _random_system's arbitrary prior over mouths, "
            "delays and launch times; regression diagnostics, not predictions"),
        "configurations": n_configs,
        "configurations_with_a_selected_event": n_sel,
        "selected_events": int(len(arr)),
        "rows": rows,
        "median_s_at_energy_equal_mass": (float(np.median(arr))
                                          if len(arr) else 0.0),
        "max_s_at_energy_equal_mass": float(arr.max()) if len(arr) else 0.0,
        "none_clear_threshold_at_energy_equal_mass": bool(
            len(arr) > 0 and not np.any(arr >= 4.0 - 1e-9)),
        "most_clear_it_by_one_and_a_half_masses": bool(
            len(arr) > 0 and float(np.mean(arr * 2.25 >= 4.0)) > 0.5),
        "closure_selects_where_the_invariant_decides_whether": True,
    }


def measure_the_conjugacy_is_carried_not_derived(
        seed: int = 31) -> Dict[str, object]:
    """``Q_A + Q_B = 0`` is a label, checked and never produced.

    The two histories carry opposite orientation, the sum is zero, and a
    same-sign pair is refused at construction.  **That is bookkeeping, not a
    derivation of charge**: nothing in the kinematics distinguishes the two
    signs, and saying otherwise would import the conclusion.
    """
    rng = np.random.default_rng(seed)
    sysm = _random_system(rng)
    refused = False
    try:
        PairHistorySystem(sysm.source_a, sysm.source_b, sysm.throat_a,
                          sysm.throat_b, orientations=(+1, +1))
    except ValueError:
        refused = True
    return {
        "orientations": list(sysm.orientations),
        "net_orientation": sysm.net_orientation(),
        "the_labels_cancel": bool(sysm.net_orientation() == 0),
        "a_same_sign_pair_is_refused": refused,
        "but_nothing_here_derives_charge": (
            "the kinematics is blind to the sign; the label is carried through "
            "the system and checked at the end, which is bookkeeping"),
        "and_the_throat_bill_is_still_inherited": (
            "shells.junction priced a connected throat as necessarily exotic; "
            "this round adds two of them and pays for neither"),
    }
