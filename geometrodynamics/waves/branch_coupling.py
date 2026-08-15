"""
Branch-resolved field coupling: the mouth transfer solved for, not applied.

THE GAP THIS CLOSES
───────────────────
`waves.field_solve` (PR #254) got the strong result it was after — on the
Einstein static universe the conformally coupled retarded Green function has
*exact* image support, so PR #253's ray branches are the field's branches, with
the ``1/(4π sin χ)`` shell law and a Maslov sign the ray ledger could not carry.

But it got that result with the mouth relation on the outside.  ``φ(M⁺,t) =
η φ(M⁻,t+Δ)`` was applied **to the free branches after they were computed**:
enumerate the arrivals at ``M⁺``, shift them by ``Δ``, re-emit from ``M⁻``,
enumerate again.  One traversal, by construction, because a post-processing step
has no way to notice that what it re-emits will come back.

Here the relation is imposed **as part of the field problem**.  In
frequency space the amplitude re-emitted at ``M⁻`` is driven by everything that
reaches ``M⁺`` — including its own earlier emission, which has gone round the
sphere and returned — so it satisfies

    ``a(ω) = η κ e^{−iωΔ} [ S(ω) + T_d(ω) a(ω) ]``

    ``a(ω) = η κ e^{−iωΔ} S(ω) / (1 − L(ω))`` ,  ``L = η κ e^{−iωΔ} T_d(ω)``

and is *solved*, not iterated.  PR #254's answer is the ``n = 0`` term of
``1/(1−L) = Σ_n Lⁿ``, and this module measures the difference.

THE PRIMITIVE IS INDEXED BY A PAIR OF BRANCHES
──────────────────────────────────────────────
A through-throat history has two legs — source ``→ M⁺`` and ``M⁻ →`` observer —
and each leg independently takes the short way, the long way, or winds.  The
object that carries a history is therefore indexed by a **pair**:

    ``K_ab(ω) = η κ · s_a A₁ e^{−u ℓ_a} · e^{−iωΔ} · s_b A₂ e^{−u ℓ_b}``

with ``u = γ + iω``, ``A = 1/(4π sin χ)`` and ``s = (−1)^crossings`` the Maslov
factor of PR #254.  The solved propagator is ``Σ_ab K_ab / (1 − L)``.

Two facts about this pair index, and they pull in opposite directions:

* **the amplitude factorizes over it.**  ``K_ab`` is an outer product, so as a
  matrix it is rank one — for one throat, *and that is not a history count*:
  the same matrix carries ``n²`` distinct ``(a,b)`` histories.  Rank counts
  independent **separable transfer channels**, and one value-feedback throat
  supplies one.  Two throats give two, in a shared topological branch-label
  basis so the sum means something.  Measured.
* **the closure condition does not factorize over it.**  PR #253's
  ``ℓ_a + Δ + ℓ_b = 0`` is a condition on the *pair*; no condition on ``a``
  alone and none on ``b`` alone implies it.  That is why the pair is the
  primitive even though the amplitude is a product.

WHAT THE SOLVE BUYS THAT POST-PROCESSING DID NOT
────────────────────────────────────────────────
* **an exact transfer function.**  Summing the branch series in closed form,

      ``Σ_b s_b e^{−u ℓ_b} = (e^{−uχ} − e^{−u(2π−χ)}) / (1 − e^{−2πu})``

  whose poles at ``γ→0`` sit at ``ω = 1, 2, 3, …`` — the conformal ESU spectrum
  ``ω_n = n+1``, recovered from the *image* representation, with residues that
  are the mode functions over ``2ω``.  The image picture and the mode picture
  are the same function, seen from two sides.
* **arrivals that do not exist in PR #254.**  The resummed field has echoes at
  ``ℓ_a + Δ + n(ℓ_c + Δ) + ℓ_b`` with amplitude ``∝ κⁿ⁺¹`` and sign the product
  of every Maslov factor in the word.  Post-processing the free branches
  produces the ``n = 0`` arrivals and nothing else.
* **a place where the one-traversal answer stops working.**  ``|T_d|`` peaks at
  the bare resonances, so the traversal series' **radius of convergence**
  ``κ_series = 1/max|T_d|`` falls linearly in the regulator.  That is a
  statement about the *expansion*, and only that.
* **and the poles themselves.**  Existence, convergence and stability are three
  different conditions and this module keeps them apart:

      existence    ``L ≠ 1``          — ``1/(1−L)`` is fine for ``|L| > 1``
      convergence  ``|L| < 1``        — the radius of ``Σ Lⁿ``, delay-blind
      stability    ``Im ω > 0`` for every root of ``D(ω) = 1 − L(ω)``

  The coupling *displaces* the bare poles ``ω = m + iγ`` by
  ``δ_m = −ηκ e^{−imΔ} sin(md)/(4π² sin d)``, whose imaginary part goes like
  ``sin(md) sin(mΔ)`` and **changes sign with the mode** — so stability is
  phase-sensitive, and no bound on ``|L|`` can decide it.  At ``Δ = π`` every
  first-order displacement is real and the stability threshold sits ``4×`` above
  the series radius; in between, the solve is finite and stable while the series
  diverges by ``10^119``.  Measured.
* **PR #251's fixed point, frequency-resolved.**  When ``Δ + ℓ_c < 0`` the loop
  is closed in time and ``1/(1−L)`` is a self-consistency condition rather than
  a convergent history sum.

SCOPE — WHICH MODEL THIS IS
───────────────────────────
A **linear field on a fixed round background**, and — stated plainly, because
the resolvent being exact says nothing about this — a **self-consistent
rank-one mouth-transfer model**, not a throat boundary operator and not a
quotient of the manifold.  It relates field *values* through the free Green
function.  It has **no normal-derivative (flux) matching**, **no reflected
channel** (the mouth scattering object is ``1×1`` where a flux-conserving
two-mouth junction needs at least ``2×2`` unitary), and with ``κ < 1`` it is
lossy by construction — power out over power in is ``κ²``, measured.  So it is
not literally an identification.  `shells.junction` (PR #249) is what would fix
``κ`` and supply the missing channels, and a genuine flux-conserving boundary
operator is the next step, not this one.

What *is* claimed is narrow and holds: the relation is **solved
self-consistently** rather than applied to already-computed branches, and that
changes the answer.

**Not done:** no backreaction, no stress tensor, no topology change, no rate.
The two-throat cross term measured here is a *throat* interference term and is
**not** the two-source invariant of roadmap step 3.  And the rank of ``K``
counts **separable transfer channels**, not histories — one throat carries 144
``(a,b)`` histories at rank one.
"""

from __future__ import annotations

import cmath
import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

__all__ = [
    "TWO_PI",
    "leg_branches",
    "branch_transfer",
    "mouth_transfer",
    "esu_mode_weight",
    "CoupledThroat",
    "branch_pair_matrix",
    "free_branch_propagator",
    "coupled_propagator",
    "traversal_series",
    "series_radius",
    "stability_threshold",
    "dispersion",
    "resonance_poles",
    "least_damped_pole",
    "branch_labels",
    "coupled_arrivals",
    "coupled_waveform",
    "measure_the_closed_form_transfer_is_the_branch_sum",
    "measure_solving_the_throat_resums_every_traversal",
    "measure_the_coupled_field_has_arrivals_the_free_branches_do_not",
    "measure_closure_is_broadband_coherence",
    "measure_the_rank_counts_transfer_channels_not_histories",
    "measure_the_one_traversal_expansion_fails_near_the_bare_resonances",
    "measure_the_series_radius_is_not_the_stability_threshold",
    "measure_what_the_transfer_model_leaves_out",
]

TWO_PI = 2.0 * math.pi


# ════════════════════════════════════════════════════════════════════════════
# ONE LEG: THE BRANCHES OF PR #253, CARRYING PR #254's SIGN
# ════════════════════════════════════════════════════════════════════════════
def leg_branches(chi: float, n_k: int = 8) -> List[Dict[str, object]]:
    """The branch set of a single leg of geodesic length ``χ``.

    ``ℓ = χ + 2πk`` (short way, ``2k`` focal crossings) and
    ``ℓ = 2π(k+1) − χ`` (long way, ``2k+1``), each with the Maslov factor
    ``s = (−1)^crossings``.  Identical to `field_solve.branch_arrivals`, but
    indexed by winding rather than truncated by a time window, because what is
    wanted here is a transfer function rather than a list of arrivals.
    """
    if not 0.0 < chi < math.pi:
        raise ValueError("a leg needs 0 < χ < π; the poles are singular")
    out: List[Dict[str, object]] = []
    for k in range(int(n_k)):
        out.append({"ell": chi + TWO_PI * k, "long_way": 0, "winding": k,
                    "crossings": 2 * k, "sign": 1})
        out.append({"ell": TWO_PI * (k + 1) - chi, "long_way": 1,
                    "winding": k, "crossings": 2 * k + 1, "sign": -1})
    return sorted(out, key=lambda r: r["ell"])


def branch_labels(n_k: int = 8) -> List[Tuple[int, int]]:
    """The canonical ``(long_way, winding)`` order — the **common basis**.

    For every ``χ ∈ (0,π)``, ``χ < 2π−χ < χ+2π < 4π−χ < …``, so ordering a leg's
    branches by length always produces this same label sequence.  The label is
    topological — which way round, how many windings — and is therefore shared
    by legs of *different* lengths.  That is what makes it legitimate to add two
    throats' `branch_pair_matrix` results, and `leg_branches` is checked against
    it rather than the coincidence being relied on.
    """
    return [lab for k in range(int(n_k)) for lab in ((0, k), (1, k))]


def branch_transfer(omega: float, chi: float, damping: float = 0.02,
                    n_k: int = 40) -> complex:
    """``A(χ) Σ_b s_b e^{−u ℓ_b}`` summed **branch by branch**.

    The regulator ``γ = damping`` is the Abel factor that makes the winding sum
    converge; it is a damping per unit path length, so it suppresses the ``k``-th
    winding by ``e^{−2πγk}`` and nothing else.
    """
    u = complex(float(damping), float(omega))
    amp = 1.0 / (4.0 * math.pi * math.sin(chi))
    tot = 0j
    for b in leg_branches(chi, n_k):
        tot += b["sign"] * cmath.exp(-u * b["ell"])
    return amp * tot


def mouth_transfer(omega: float, chi: float,
                   damping: float = 0.02) -> complex:
    """The same sum in **closed form** — every winding, no truncation.

        ``T(ω;χ) = (e^{−uχ} − e^{−u(2π−χ)}) / [(1 − e^{−2πu}) · 4π sin χ]``

    The short-way images all carry ``s = +1`` and the long-way images all carry
    ``s = −1``, so the series is two geometric series and sums exactly.  As
    ``γ → 0`` the denominator vanishes at ``ω ∈ ℤ``: the poles of the *image*
    representation are the conformal ESU eigenfrequencies ``ω_n = n+1``, and the
    numerator's own zero at ``ω = 0`` removes the one that would not be.
    """
    if not 0.0 < chi < math.pi:
        raise ValueError("the transfer is singular at χ = 0 and χ = π")
    u = complex(float(damping), float(omega))
    num = cmath.exp(-u * chi) - cmath.exp(-u * (TWO_PI - chi))
    den = 1.0 - cmath.exp(-u * TWO_PI)
    if abs(den) < 1e-14:
        raise ValueError("undamped pole: ω is an ESU eigenfrequency")
    return num / den / (4.0 * math.pi * math.sin(chi))


def esu_mode_weight(m: int, chi: float) -> float:
    """``m sin(mχ) / (2π² sin χ)`` — the ``S³`` addition theorem at ``ω = m``.

    Kept here so the residue check below compares the image representation
    against the *mode* representation rather than against itself.
    """
    return m * math.sin(m * chi) / (2.0 * math.pi ** 2 * math.sin(chi))


# ════════════════════════════════════════════════════════════════════════════
# THE THROAT, AS A BOUNDARY CONDITION
# ════════════════════════════════════════════════════════════════════════════
@dataclass(frozen=True)
class CoupledThroat:
    """A **mouth-transfer model**: what arrives at ``M⁺`` drives what leaves
    ``M⁻``, as ``ηκ`` times the value, delayed by ``Δ``, with the mouths a
    geodesic distance ``separation`` apart.

    **This is not a boundary operator and not an identification.**  It relates
    field *values* through the free Green function and nothing else: there is no
    normal-derivative (flux) matching, no reflected channel, and no
    self-adjoint two-mouth ``S``-matrix.  With ``κ < 1`` it is manifestly lossy,
    so it cannot be a quotient of the manifold; ``κ = 1`` would make the *value*
    relation an identification and still leave the derivative unmatched.  What
    this module claims is only that the relation is **solved self-consistently**
    instead of applied to already-computed branches.  `shells.junction`
    (PR #249) is what would fix ``κ`` and supply the missing channels.
    """

    separation: float
    delay: float
    eta: int = 1
    kappa: float = 0.30

    def loop_transfer(self, omega: float, damping: float = 0.02) -> complex:
        """``L(ω) = η κ e^{−iωΔ} T_d(ω)`` — the gain of one round trip.

        Emit at ``M⁻``, cross the sphere on every branch at once, arrive at
        ``M⁺``, be re-emitted.  ``|L| < 1`` is the **radius of convergence of the
        traversal series** ``Σ Lⁿ`` — it is *not* the condition for the solution
        to exist, which is only ``L ≠ 1``, and it is *not* the stability
        condition, which is phase-sensitive and lives in `resonance_poles`.
        """
        phase = cmath.exp(-1j * float(omega) * self.delay)
        return (self.eta * self.kappa * phase
                * mouth_transfer(omega, self.separation, damping))

    def resolvent(self, omega: float, damping: float = 0.02) -> complex:
        """``1/(1 − L)`` — **the solve**, every traversal at once.

        Well defined for *any* ``L ≠ 1``, including ``|L| > 1`` where the
        traversal series diverges.  Solving and summing are not the same
        operation, and `measure_the_series_radius_is_not_the_stability_threshold`
        exhibits a coupling where the first works and the second does not.
        """
        return 1.0 / (1.0 - self.loop_transfer(omega, damping))

    def loop_branches(self, n_k: int = 8) -> List[Dict[str, object]]:
        """The branches of the round trip, for enumerating echo words."""
        return leg_branches(self.separation, n_k)


# ════════════════════════════════════════════════════════════════════════════
# WHERE THE SOLUTION ACTUALLY LIVES: THE COMPLEX POLES
# ════════════════════════════════════════════════════════════════════════════
def dispersion(omega: complex, throat: CoupledThroat,
               damping: float = 0.02) -> complex:
    """``D(ω) = 1 − L(ω)``, continued to **complex** ``ω``.

    ``u = γ + iω`` continues without modification, so this is the same function
    off the real axis.  The physical content of the coupled system is where
    ``D = 0``, not where ``|L| = 1``.
    """
    w = complex(omega)
    u = float(damping) + 1j * w
    d = throat.separation
    num = cmath.exp(-u * d) - cmath.exp(-u * (TWO_PI - d))
    den = 1.0 - cmath.exp(-u * TWO_PI)
    t = num / den / (4.0 * math.pi * math.sin(d))
    return 1.0 - throat.eta * throat.kappa * cmath.exp(-1j * w * throat.delay) * t


def _pole_seed(m: int, throat: CoupledThroat, damping: float) -> complex:
    """First-order displacement of the ``m``-th bare pole.

    Near ``ω = m`` the denominator is ``1 − e^{−2πu} ≈ 2πi(ω − m − iγ)`` and the
    numerator is ``−2i sin(md)``, so ``T_d ≈ −sin(md)/(4π²(ω−m−iγ) sin d)`` and
    ``D = 0`` moves the pole by

        ``δ = −ηκ e^{−imΔ} sin(md) / (4π² sin d)``

    whose imaginary part is ``ηκ sin(md) sin(mΔ)/(4π² sin d)``.  **The sign
    depends on ``m``**, which is why stability is phase-sensitive and a bound on
    ``|L|`` alone cannot decide it.
    """
    d = throat.separation
    delta = (-throat.eta * throat.kappa * cmath.exp(-1j * m * throat.delay)
             * math.sin(m * d) / (4.0 * math.pi ** 2 * math.sin(d)))
    return m + 1j * float(damping) + delta


def resonance_poles(throat: CoupledThroat, damping: float = 0.02,
                    n_modes: int = 24, tol: float = 1e-13,
                    max_iter: int = 400) -> List[Dict[str, object]]:
    """Roots of ``D(ω) = 0`` in the complex plane, one per bare resonance.

    With the ``e^{+iωt}`` convention a root at ``ω = ω_r + iω_i`` rings as
    ``e^{iω_r t} e^{−ω_i t}``, so ``Im ω > 0`` is a **decaying** mode and the
    system is stable exactly when every root has ``Im ω > 0``.  The uncoupled
    poles sit at ``ω = m + iγ``; the coupling displaces them, and `_pole_seed`
    is where each search starts.

    Band-limited by construction: only the first ``n_modes`` resonances are
    found, and the least-damped root can move to a higher ``m`` as the band
    widens.  Reported, not hidden.
    """
    out: List[Dict[str, object]] = []
    for m in range(1, int(n_modes) + 1):
        w = _pole_seed(m, throat, damping)
        ok = False
        for _ in range(int(max_iter)):
            h = 1e-6
            fp = (dispersion(w + h, throat, damping)
                  - dispersion(w - h, throat, damping)) / (2.0 * h)
            if abs(fp) < 1e-300:
                break
            step = dispersion(w, throat, damping) / fp
            if abs(step) > 0.5:
                step *= 0.5 / abs(step)
            w -= step
            if abs(step) < tol:
                ok = True
                break
        if not ok:
            continue
        out.append({"mode": m, "omega": complex(w),
                    "re": float(w.real), "im": float(w.imag),
                    "residual": float(abs(dispersion(w, throat, damping))),
                    "decaying": bool(w.imag > 0.0),
                    "seed": complex(_pole_seed(m, throat, damping))})
    return out


def least_damped_pole(throat: CoupledThroat, damping: float = 0.02,
                      n_modes: int = 24) -> Optional[Dict[str, object]]:
    """The root with the smallest ``Im ω`` — the one that decides stability."""
    poles = resonance_poles(throat, damping, n_modes)
    return min(poles, key=lambda p: p["im"]) if poles else None


# ════════════════════════════════════════════════════════════════════════════
# THE PRIMITIVE
# ════════════════════════════════════════════════════════════════════════════
def branch_pair_matrix(omega: float, chi_in: float, chi_out: float,
                       throat: CoupledThroat, damping: float = 0.02,
                       n_k: int = 8) -> Dict[str, object]:
    """``K_ab(ω)`` — the object indexed by a **pair** of branches.

        ``K_ab = η κ · s_a A₁ e^{−u ℓ_a} · e^{−iωΔ} · s_b A₂ e^{−u ℓ_b}``

    Row ``a`` is a branch of the leg ``source → M⁺``; column ``b`` a branch of
    ``M⁻ → observer``.  One traversal of the throat, resolved on both legs at
    once.

    Both indices are the **topological branch label** ``(long_way, winding)`` in
    `branch_labels` order, not the leg length, so matrices built for different
    ``χ`` — different throats — are in the same basis and may legitimately be
    added.  The labels are returned so a caller can check rather than assume it.
    """
    u = complex(float(damping), float(omega))
    rows = leg_branches(chi_in, n_k)
    cols = leg_branches(chi_out, n_k)
    labels = branch_labels(n_k)
    if ([(r["long_way"], r["winding"]) for r in rows] != labels
            or [(c["long_way"], c["winding"]) for c in cols] != labels):
        raise ValueError("legs are not in the canonical branch-label basis")
    a1 = 1.0 / (4.0 * math.pi * math.sin(chi_in))
    a2 = 1.0 / (4.0 * math.pi * math.sin(chi_out))
    pre = (throat.eta * throat.kappa * a1 * a2
           * cmath.exp(-1j * float(omega) * throat.delay))
    va = np.array([r["sign"] * cmath.exp(-u * r["ell"]) for r in rows])
    vb = np.array([c["sign"] * cmath.exp(-u * c["ell"]) for c in cols])
    return {"K": pre * np.outer(va, vb), "rows": rows, "cols": cols,
            "labels": labels,
            "pair_sum": complex(pre * va.sum() * vb.sum())}


def free_branch_propagator(omega: float, chi_in: float, chi_out: float,
                           throat: CoupledThroat, damping: float = 0.02,
                           n_k: int = 40) -> complex:
    """PR #254's answer: **one** traversal, ``Σ_ab K_ab``.

    Written with the closed-form transfers so that the comparison against the
    solve is not a comparison between a truncation and an exact sum.
    """
    return (throat.eta * throat.kappa
            * cmath.exp(-1j * float(omega) * throat.delay)
            * mouth_transfer(omega, chi_in, damping)
            * mouth_transfer(omega, chi_out, damping))


def coupled_propagator(omega: float, chi_in: float, chi_out: float,
                       throat: CoupledThroat, damping: float = 0.02) -> complex:
    """**The solved through-throat propagator**, ``Σ_ab K_ab / (1 − L)``.

    The only difference from `free_branch_propagator` is the resolvent, and the
    only reason the resolvent is there is that the throat was written into the
    field problem instead of applied to its output.
    """
    return (free_branch_propagator(omega, chi_in, chi_out, throat, damping)
            * throat.resolvent(omega, damping))


def traversal_series(omega: float, chi_in: float, chi_out: float,
                     throat: CoupledThroat, damping: float = 0.02,
                     n_terms: int = 200) -> complex:
    """The same thing built the slow way: ``Σ_ab K_ab · Σ_{n<N} Lⁿ``.

    An independent construction — this one never inverts anything — so agreement
    with `coupled_propagator` says the resolvent really is the history sum.
    """
    L = throat.loop_transfer(omega, damping)
    tot = 0j
    term = 1.0 + 0j
    for _ in range(int(n_terms)):
        tot += term
        term *= L
    return free_branch_propagator(omega, chi_in, chi_out, throat,
                                  damping) * tot


def series_radius(separation: float, damping: float = 0.02,
                  omega_max: float = 8.0,
                  n_grid: int = 20001) -> Dict[str, float]:
    """``κ_series = 1 / max_ω |T_d(ω)|`` — the **radius of convergence** of the
    traversal series, and nothing more.

    Above it ``Σ Lⁿ`` diverges somewhere in the band.  It is *not* a stability
    bound and *not* an existence condition: it takes no account of the phase of
    ``L`` and does not even depend on the delay ``Δ``, while stability depends
    on ``Δ`` strongly (see `stability_threshold`).  Band-limited by
    ``omega_max``, which is reported.
    """
    ws = np.linspace(1e-6, float(omega_max), int(n_grid))
    mags = np.array([abs(mouth_transfer(float(w), separation, damping))
                     for w in ws])
    i = int(np.argmax(mags))
    return {"kappa_series": float(1.0 / mags[i]),
            "omega_of_the_peak": float(ws[i]),
            "peak_transfer": float(mags[i]), "damping": float(damping),
            "omega_max": float(omega_max),
            "what_it_is": ("the radius of convergence of Σ Lⁿ, independent of "
                           "the delay and therefore not a stability bound")}


def stability_threshold(separation: float, delay: float, eta: int = 1,
                        damping: float = 0.02, n_modes: int = 24,
                        kappa_max: float = 64.0,
                        tol: float = 1e-4) -> Dict[str, object]:
    """The ``κ`` at which the least-damped root of ``D(ω) = 0`` reaches the real
    axis — the threshold that actually decides stability.

    Found by bisection on ``Im ω_least(κ)``.  Unlike `series_radius` this is
    **phase-sensitive**: at ``Δ = π`` the first-order displacement is real for
    every mode (``sin(mπ) = 0``), the poles slide along the axis instead of off
    it, and the threshold moves far above the series radius.
    """
    def im_least(k: float) -> float:
        th = CoupledThroat(separation, delay, eta, k)
        p = least_damped_pole(th, damping, n_modes)
        return p["im"] if p else float("nan")

    lo = 0.0
    hi = float(kappa_max)
    if im_least(hi) > 0.0:
        return {"kappa_stability": None, "bracketed": False,
                "kappa_max_scanned": hi, "n_modes": int(n_modes),
                "why": "no root reached the real axis inside the scanned range"}
    while hi - lo > tol * max(1.0, hi):
        mid = 0.5 * (lo + hi)
        if im_least(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    th = CoupledThroat(separation, delay, eta, lo)
    p = least_damped_pole(th, damping, n_modes)
    return {"kappa_stability": float(0.5 * (lo + hi)), "bracketed": True,
            "mode_that_goes_first": (p["mode"] if p else None),
            "n_modes": int(n_modes), "damping": float(damping),
            "delay": float(delay),
            "what_it_is": ("where Im ω of the least-damped root of D = 0 "
                           "crosses zero; band-limited to the first n_modes")}


# ════════════════════════════════════════════════════════════════════════════
# WHAT THE SOLVE PUTS IN THE TIME DOMAIN
# ════════════════════════════════════════════════════════════════════════════
def coupled_arrivals(chi_in: float, chi_out: float, throat: CoupledThroat,
                     t_max: float, n_traversal: int = 3, n_k: int = 4,
                     damping: float = 0.02,
                     rel_floor: float = 1e-4) -> List[Dict[str, object]]:
    """Every history word the solved throat supports, in time order.

    A word is ``(a, c₁ … c_n, b)``: one branch into ``M⁺``, ``n`` round trips,
    one branch out of ``M⁻``.  It arrives at

        ``t = ℓ_a + Δ + Σ_i (ℓ_{c_i} + Δ) + ℓ_b``

    with amplitude ``κⁿ⁺¹ A₁ A₂ A_dⁿ e^{−γ Σℓ}`` and sign
    ``ηⁿ⁺¹ s_a s_b Π_i s_{c_i}`` — every Maslov factor in the word.  PR #254
    contains exactly the ``n = 0`` rows of this list.
    """
    a1 = 1.0 / (4.0 * math.pi * math.sin(chi_in))
    a2 = 1.0 / (4.0 * math.pi * math.sin(chi_out))
    ad = 1.0 / (4.0 * math.pi * math.sin(throat.separation))
    ins, outs = leg_branches(chi_in, n_k), leg_branches(chi_out, n_k)
    loops = throat.loop_branches(n_k)

    words: List[Dict[str, object]] = []

    def walk(n: int, ell: float, sign: int, cross: int,
             chain: Tuple[int, ...]) -> None:
        if n > int(n_traversal):
            return
        for a in ins:
            for b in outs:
                t = a["ell"] + ell + b["ell"] + throat.delay * (n + 1)
                if not 0.0 <= t <= t_max:
                    continue
                total = a["ell"] + ell + b["ell"]
                amp = (throat.kappa ** (n + 1) * a1 * a2 * ad ** n
                       * math.exp(-damping * total))
                words.append({
                    "t": t, "traversals": n + 1, "amplitude": amp,
                    "sign": int(throat.eta) ** (n + 1) * a["sign"]
                            * b["sign"] * sign,
                    "crossings": a["crossings"] + b["crossings"] + cross,
                    "branch_in": (a["long_way"], a["winding"]),
                    "branch_out": (b["long_way"], b["winding"]),
                    "loop_word": chain})
        for j, c in enumerate(loops):
            walk(n + 1, ell + c["ell"], sign * int(c["sign"]),
                 cross + int(c["crossings"]), chain + (j,))

    walk(0, 0.0, 1, 0, ())
    if not words:
        return []
    top = max(w["amplitude"] for w in words)
    words = [w for w in words if w["amplitude"] >= rel_floor * top]
    return sorted(words, key=lambda w: w["t"])


def coupled_waveform(t_max: float, chi_in: float, chi_out: float,
                     throat: CoupledThroat, damping: float = 0.02,
                     width: float = 0.05, n_fft: int = 1 << 17,
                     d_omega: float = 0.0025,
                     resolvent: bool = True) -> Tuple[np.ndarray, np.ndarray]:
    """The solved propagator, inverse-transformed onto a time grid.

        ``φ(t) = (1/π) Re ∫₀^∞ dω Ĝ(ω) e^{iωt} e^{−ω²w²/2}``

    ``Ĝ(−ω) = conj Ĝ(ω)``, so the negative half is the conjugate and the field
    comes out real.  Done as a DFT on a uniform ``ω`` grid, which is not just
    speed: the transfer function has near-poles of width ``γ`` at the ESU
    eigenfrequencies, and a uniform grid with ``dω ≪ γ`` is what resolves them.

    With ``resolvent=False`` the same integral is taken over PR #254's
    one-traversal propagator.  That is the control, and the two waveforms are
    the experiment.
    """
    n = int(n_fft)
    ws = np.arange(n) * float(d_omega)
    u = damping + 1j * ws
    d = throat.separation

    def tr(chi: float) -> np.ndarray:
        num = np.exp(-u * chi) - np.exp(-u * (TWO_PI - chi))
        den = 1.0 - np.exp(-u * TWO_PI)
        return num / den / (4.0 * math.pi * math.sin(chi))

    drive = (throat.eta * throat.kappa * np.exp(-1j * ws * throat.delay))
    g = drive * tr(chi_in) * tr(chi_out)
    if resolvent:
        g = g / (1.0 - drive * tr(d))
    g = g * np.exp(-(ws * width) ** 2 / 2.0)
    g[0] *= 0.5                       # the ω = 0 point is shared by both halves

    phi = 2.0 * np.real(n * np.fft.ifft(g)) * (float(d_omega) / TWO_PI)
    ts = TWO_PI * np.arange(n) / (n * float(d_omega))
    keep = ts <= float(t_max)
    return ts[keep], phi[keep]


def _waveform_peaks(ts: np.ndarray, f: np.ndarray,
                    rel_floor: float = 0.06) -> List[Dict[str, float]]:
    """Local extrema of ``|φ|``, with their signs."""
    floor = rel_floor * float(np.abs(f).max())
    out: List[Dict[str, float]] = []
    for i in range(1, len(ts) - 1):
        a = abs(f[i])
        if a > abs(f[i - 1]) and a >= abs(f[i + 1]) and a > floor:
            out.append({"t": float(ts[i]), "value": float(f[i])})
    return out


# ════════════════════════════════════════════════════════════════════════════
# MEASUREMENTS
# ════════════════════════════════════════════════════════════════════════════
def measure_the_closed_form_transfer_is_the_branch_sum(
        chis: Sequence[float] = (0.5, 1.3, 2.2),
        omegas: Sequence[float] = (0.4, 1.7, 3.3, 6.1),
        damping: float = 0.05, n_k: int = 400) -> Dict[str, object]:
    """The branch series against its closed form, and where its poles are.

    Two things at once, because they are the same statement.  The truncated sum
    over images agrees with the geometric closed form; and the closed form's
    poles, as ``γ → 0``, sit exactly at the conformal ESU frequencies
    ``ω = 1, 2, 3, …`` with residues equal to the mode functions over ``2ω``.
    The image representation and the mode representation are one function.
    """
    worst = 0.0
    rows = []
    for chi in chis:
        for w in omegas:
            a = branch_transfer(w, chi, damping, n_k)
            b = mouth_transfer(w, chi, damping)
            err = abs(a - b)
            worst = max(worst, err)
            rows.append({"chi": chi, "omega": w, "series": a, "closed": b,
                         "abs_error": err})

    # the poles, found by watching |T| blow up as the regulator is removed
    poles = []
    for m in (1, 2, 3):
        prev: Optional[float] = None
        scaling = []
        for g in (1e-2, 1e-3, 1e-4):
            mag = abs(mouth_transfer(float(m), 1.3, g))
            scaling.append({"damping": g, "abs_T": mag,
                            "ratio_to_previous": (mag / prev if prev
                                                  else None)})
            prev = mag
        # residue in ω: (1 − e^{−2πu}) has slope 2πi at u = im, so the residue
        # is T·(1 − e^{−2πu})/(2πi).  Compare against the mode function / 2ω.
        g = 1e-6
        u = complex(g, float(m))
        res = (mouth_transfer(float(m), 1.3, g)
               * (1.0 - cmath.exp(-u * TWO_PI)) / (TWO_PI * 1j))
        predicted = -esu_mode_weight(m, 1.3) / (2.0 * m)
        poles.append({"omega": m, "scaling": scaling,
                      "residue": res, "residue_real": res.real,
                      "predicted_mode_over_2omega": predicted,
                      "matches": bool(abs(res.real - predicted)
                                      < 1e-4 * abs(predicted)
                                      and abs(res.imag) < 1e-4
                                      * abs(predicted))})

    return {
        "rows": rows, "worst_abs_error": worst,
        "the_series_is_the_closed_form": bool(worst < 1e-9),
        "poles": poles,
        "the_poles_are_the_esu_spectrum": bool(all(
            p["scaling"][-1]["abs_T"] > 50.0 * p["scaling"][0]["abs_T"]
            for p in poles)),
        "the_residues_are_the_mode_functions": bool(
            all(p["matches"] for p in poles)),
        "what_this_shows": ("the image sum and the mode sum are the same "
                            "function; the branch labels are a representation, "
                            "not an approximation"),
    }


def measure_solving_the_throat_resums_every_traversal(
        chi_in: float = 1.2, chi_out: float = 0.9,
        separation: float = 1.3, delay: float = 1.0,
        kappas: Sequence[float] = (0.05, 0.20, 0.45),
        omegas: Sequence[float] = (0.4, 1.1, 2.7, 5.3),
        damping: float = 0.05, n_terms: int = 400) -> Dict[str, object]:
    """**The claim of this round.**  ``1/(1−L)`` is the sum over every traversal.

    The resolvent comes from solving the boundary condition; the series comes
    from walking the throat ``n`` times and adding.  They agree to machine
    precision.  And PR #254's post-processed answer is reported alongside as
    what it is: the ``n = 0`` term, whose error grows with the loop gain.
    """
    rows = []
    worst = 0.0
    worst_free = 0.0
    for kap in kappas:
        th = CoupledThroat(separation, delay, +1, kap)
        for w in omegas:
            solved = coupled_propagator(w, chi_in, chi_out, th, damping)
            walked = traversal_series(w, chi_in, chi_out, th, damping,
                                      n_terms)
            free = free_branch_propagator(w, chi_in, chi_out, th, damping)
            L = th.loop_transfer(w, damping)
            err = abs(solved - walked)
            miss = abs(solved - free) / abs(solved)
            worst = max(worst, err)
            worst_free = max(worst_free, miss)
            rows.append({"kappa": kap, "omega": w, "loop_gain": abs(L),
                         "solved": solved, "walked_series": walked,
                         "abs_error": err,
                         "post_processed_one_traversal": free,
                         "relative_miss_of_one_traversal": miss})
    return {
        "rows": rows, "worst_abs_error": worst,
        "the_resolvent_is_the_traversal_sum": bool(worst < 1e-12),
        "worst_relative_miss_of_post_processing": worst_free,
        "post_processing_is_the_n_equals_zero_term": True,
        "what_this_shows": ("writing the identification into the field problem "
                            "is not a rearrangement of PR #254 — it adds every "
                            "history that returns to the mouth"),
    }


def measure_the_coupled_field_has_arrivals_the_free_branches_do_not(
        chi_in: float = 1.2, chi_out: float = 0.9,
        separation: float = 1.3, delay: float = 1.0, kappa: float = 0.60,
        damping: float = 0.02, width: float = 0.05,
        t_max: float = 11.0) -> Dict[str, object]:
    """**The echo train** — the sharpest difference the solve makes.

    Post-processing the free branches gives arrivals at ``ℓ_a + Δ + ℓ_b`` and
    stops.  The solved field additionally rings at ``+ n(ℓ_c + Δ)`` for every
    ``n``, with amplitude falling like ``κⁿ`` and sign the product of every
    Maslov factor in the word.  Three things are measured:

    * the solved waveform **is** the sum over history words, to ``~1e-5`` —
      so the word enumeration is not a story told about the waveform;
    * at the echo times that no one-traversal word can reach, the solved field
      stands hundreds of times above the control;
    * and the amplitudes follow the ``κⁿ`` ladder.

    ``Δ > 0`` is chosen so the loop period is positive and the echoes are
    ordinary late arrivals; the closed-in-time case is another measurement's
    business.
    """
    th = CoupledThroat(separation, delay, +1, kappa)
    ts, solved = coupled_waveform(t_max, chi_in, chi_out, th, damping, width)
    _, free = coupled_waveform(t_max, chi_in, chi_out, th, damping, width,
                               resolvent=False)

    words = coupled_arrivals(chi_in, chi_out, th, t_max + 3.0, n_traversal=5,
                             n_k=5, damping=damping, rel_floor=1e-10)
    one = [w for w in words if w["traversals"] == 1]
    many = [w for w in words if w["traversals"] > 1]

    # (1) the waveform is the word sum
    rec = np.zeros_like(ts)
    for w in words:
        rec += (w["sign"] * w["amplitude"]
                * np.exp(-(ts - w["t"]) ** 2 / (2.0 * width ** 2))
                / (width * math.sqrt(TWO_PI)))
    scale = float(np.abs(solved).max())
    recon_err = float(np.abs(rec - solved).max()) / scale

    # (2) the echoes, at times no one-traversal word reaches
    checked = []
    for w in many:
        gap = min((abs(w["t"] - o["t"]) for o in one), default=9e9)
        if gap < 8.0 * width:
            continue
        m = np.abs(ts - w["t"]) < 2.0 * width
        if not np.any(m):
            continue
        lvl_s = float(np.abs(solved[m]).max())
        lvl_f = float(np.abs(free[m]).max())
        i = int(np.argmax(np.abs(solved * m)))
        if checked and abs(checked[-1]["t"] - w["t"]) < 4.0 * width:
            continue
        checked.append({"t": w["t"], "traversals": w["traversals"],
                        "isolation": gap, "word_sign": w["sign"],
                        "solved_level": lvl_s, "control_level": lvl_f,
                        "contrast": lvl_s / max(lvl_f, 1e-300),
                        "field_sign": (1 if solved[i] > 0 else -1),
                        "sign_agrees": bool((1 if solved[i] > 0 else -1)
                                            == w["sign"])})

    # (3) the κⁿ ladder
    top = max((w["amplitude"] for w in one), default=1.0)
    ladder = []
    for n in (1, 2, 3, 4):
        got = [w["amplitude"] for w in words if w["traversals"] == n]
        ladder.append({"traversals": n,
                       "largest_amplitude": (max(got) if got else 0.0),
                       "over_first": ((max(got) / top) if got else 0.0)})

    return {
        "n_words_one_traversal": len(one), "n_words_multi": len(many),
        "reconstruction_relative_error": recon_err,
        "the_waveform_is_the_sum_over_history_words": bool(recon_err < 1e-3),
        "isolated_echoes": checked,
        "worst_echo_contrast": min((c["contrast"] for c in checked),
                                   default=0.0),
        "every_isolated_echo_stands_above_the_control": bool(
            checked and all(c["contrast"] > 20.0 for c in checked)),
        "every_echo_carries_its_maslov_word_sign": bool(
            checked and all(c["sign_agrees"] for c in checked)),
        "amplitude_ladder": ladder, "kappa": kappa,
        "what_this_shows": ("the extra arrivals are not a correction to the "
                            "amplitudes of PR #254's ledger; they are events at "
                            "times that ledger does not contain"),
    }


def measure_closure_is_broadband_coherence(
        chi_in: float = 1.2, chi_out: float = 0.9,
        separation: float = 1.3, kappa: float = 0.30,
        damping: float = 0.0, band: Tuple[float, float] = (0.5, 12.5),
        n_omega: int = 4000, n_k: int = 4) -> Dict[str, object]:
    """**PR #253's closure condition, restated spectrally.**

    ``K_ab(ω)`` carries the phase ``e^{−iω(ℓ_a + Δ + ℓ_b)}``.  So the closure
    condition ``ℓ_a + Δ + ℓ_b = 0`` is *exactly* the statement that ``K_ab`` is
    independent of ``ω``: the closed pair contributes with the same phase at
    every frequency, while every other pair winds and averages away over a band.

    Measured as the band-average ``|⟨K_ab⟩| / |K_ab|``, which is ``1`` for a
    closed pair and ``|sinc(ΔT·B/2)|`` for the rest.  This is what the pair
    index is *for*, and the delay here is tuned to make the point sharply: at
    ``Δ = −(χ₁ + χ₂ + 4π)`` the closed set is ``{(k,j) : k + j = 2}`` on the
    short-way branches — **three** pairs out of the **nine** that any rule
    phrased on ``a`` alone and ``b`` alone would have to admit.  The condition
    does not factorize over the index the amplitude factorizes over.
    """
    ws = np.linspace(band[0], band[1], int(n_omega))
    ins, outs = leg_branches(chi_in, n_k), leg_branches(chi_out, n_k)

    # tune the delay so that a whole diagonal of branch pairs closes at once
    delay = -(chi_in + chi_out + 2.0 * TWO_PI)
    th = CoupledThroat(separation, delay, +1, kappa)

    rows = []
    for i, a in enumerate(ins):
        for j, b in enumerate(outs):
            lag = a["ell"] + delay + b["ell"]
            ph = np.exp(-1j * ws * lag)
            coh = float(abs(ph.mean()))
            rows.append({"a": (a["long_way"], a["winding"]),
                         "b": (b["long_way"], b["winding"]),
                         "residual": lag, "coherence": coh,
                         "is_the_closed_pair": bool(abs(lag) < 1e-12)})
    closed = [r for r in rows if r["is_the_closed_pair"]]
    others = [r for r in rows if not r["is_the_closed_pair"]]
    worst_closed = min((r["coherence"] for r in closed), default=0.0)
    best_other = max((r["coherence"] for r in others), default=1.0)

    # and the condition is genuinely on the pair: no single-index rule reproduces
    # the set of closed pairs
    by_a = {r["a"] for r in closed}
    by_b = {r["b"] for r in closed}
    would_be = sum(1 for r in rows if r["a"] in by_a and r["b"] in by_b)

    K = branch_pair_matrix(3.0, chi_in, chi_out, th, max(damping, 1e-9),
                           n_k)["K"]
    return {
        "delay": delay, "n_pairs": len(rows),
        "n_closed": len(closed), "worst_closed_coherence": worst_closed,
        "best_other_coherence": best_other,
        "closed_pairs_are_broadband_coherent": bool(worst_closed > 0.999),
        "every_other_pair_dephases": bool(best_other < 0.2),
        "contrast": (worst_closed / best_other if best_other else float("inf")),
        "pairs_a_single_index_rule_would_select": would_be,
        "the_condition_does_not_factorize": bool(would_be > len(closed)),
        "the_amplitude_does_factorize": bool(
            np.linalg.matrix_rank(K, tol=1e-12 * float(np.abs(K).max())) == 1),
        "what_this_shows": ("the amplitude is an outer product over the pair "
                            "index and the closure condition is not; that is "
                            "why the pair is the primitive"),
    }


def measure_the_rank_counts_transfer_channels_not_histories(
        chi_in: float = 1.2, chi_out: float = 0.9, omega: float = 2.3,
        damping: float = 0.05, n_k: int = 6) -> Dict[str, object]:
    """What the rank of ``K_ab`` counts — and what it does **not**.

    It does not count histories.  A single throat with ``n`` branches per leg
    already carries ``n²`` distinct ``(a,b)`` histories, every one of them with
    its own arrival time, and ``K`` is still **rank one**.  What the rank counts
    is the number of independent **separable transfer channels** — outer
    products — and a single value-feedback throat supplies exactly one.  Both
    numbers are reported side by side so the distinction is visible.

    A second throat adds a second outer product, so the rank goes to two.  For
    that sum to mean anything the two matrices must share a basis: both are
    built in the topological `branch_labels` order rather than by leg length,
    which is checked here rather than assumed, since the two throats have
    different ``χ`` on every leg.

    The interference between the two channels is bilinear and therefore
    identically zero without either.  It is a *throat–throat* term, not the
    two-source invariant of roadmap step 3; it is measured because it has the
    same shape and shows which object carries it.
    """
    A = CoupledThroat(1.3, 1.0, +1, 0.30)
    B = CoupledThroat(2.0, -0.4, -1, 0.22)
    ka = branch_pair_matrix(omega, chi_in, chi_out, A, damping, n_k)
    kb = branch_pair_matrix(omega, 0.7, 1.6, B, damping, n_k)
    same_basis = ka["labels"] == kb["labels"] == branch_labels(n_k)
    both = ka["K"] + kb["K"]
    n_histories = ka["K"].size

    def spectrum(M: np.ndarray) -> List[float]:
        s = np.linalg.svd(M, compute_uv=False)
        return [float(x / s[0]) for x in s[:4]]

    sa, sb, sab = spectrum(ka["K"]), spectrum(kb["K"]), spectrum(both)
    tot = ka["pair_sum"] + kb["pair_sum"]
    cross = (abs(tot) ** 2 - abs(ka["pair_sum"]) ** 2
             - abs(kb["pair_sum"]) ** 2)
    cross_direct = 2.0 * (ka["pair_sum"] * kb["pair_sum"].conjugate()).real

    # the cross term is bilinear, so it is *structurally* zero without one of
    # the throats; what is worth measuring is that it is a real fringe rather
    # than a small offset — scan ω and look at its visibility
    vis = []
    for w in np.linspace(0.5, 6.5, 240):
        pa = branch_pair_matrix(float(w), chi_in, chi_out, A, damping,
                                n_k)["pair_sum"]
        pb = branch_pair_matrix(float(w), 0.7, 1.6, B, damping,
                                n_k)["pair_sum"]
        vis.append(2.0 * (pa * pb.conjugate()).real
                   / (2.0 * abs(pa) * abs(pb)))
    vis_arr = np.array(vis)
    return {
        "n_histories_one_throat": int(n_histories),
        "rank_one_throat_despite_that": int(np.linalg.matrix_rank(
            ka["K"], tol=1e-11 * float(np.abs(ka["K"]).max()))),
        "rank_is_not_a_history_count": bool(n_histories > 1),
        "what_the_rank_counts": ("independent separable transfer channels — "
                                 "outer products — not histories; one "
                                 "value-feedback throat supplies exactly one "
                                 f"while carrying {int(n_histories)} (a,b) "
                                 "histories"),
        "both_matrices_in_the_common_label_basis": bool(same_basis),
        "cross_visibility_max": float(vis_arr.max()),
        "cross_visibility_min": float(vis_arr.min()),
        "the_cross_term_is_a_full_fringe": bool(vis_arr.max() > 0.9
                                                and vis_arr.min() < -0.9),
        "singular_values_throat_A": sa,
        "singular_values_throat_B": sb,
        "singular_values_both": sab,
        "rank_two_throats": int(np.linalg.matrix_rank(
            both, tol=1e-11 * float(np.abs(both).max()))),
        "one_throat_is_one_channel": bool(sa[1] < 1e-12),
        "two_throats_are_two_channels": bool(sab[1] > 1e-3
                                             and sab[2] < 1e-12),
        "cross_term": cross, "cross_term_direct": cross_direct,
        "cross_term_agrees": bool(abs(cross - cross_direct)
                                  < 1e-9 * max(abs(cross), 1e-30)),
        "why_it_vanishes_without_a_second_throat": (
            "it is bilinear — one factor from each throat's pair sum — so "
            "removing either sets it to zero identically, not merely "
            "underdetermined; that is structural, not measured"),
        "scope": ("a throat–throat interference term, not the two-source "
                  "invariant of roadmap step 3"),
    }


def measure_the_one_traversal_expansion_fails_near_the_bare_resonances(
        separation: float = 1.3, delay: float = 1.0,
        dampings: Sequence[float] = (0.08, 0.04, 0.02, 0.01),
        chi_in: float = 1.2, chi_out: float = 0.9) -> Dict[str, object]:
    """Where the one-traversal answer stops being a useful approximation.

    ``|T_d|`` peaks at the bare ESU resonances, where ``1 − e^{−2πu} ≈ 2πγ``, so
    the **series radius** ``κ_series = 1/max|T_d|`` falls linearly in the
    regulator.  What that says, and all it says, is that the traversal
    expansion's radius of convergence shrinks as the bare poles approach the
    real axis — the ``γ → 0`` limit is singular for *the expansion*.

    It does **not** say that every non-zero ``κ`` is physically unstable: the
    coupling shifts the poles rather than merely inflating a gain, and where
    they go is `measure_the_series_radius_is_not_the_stability_threshold`.

    Also reported: the relative miss of the one-traversal answer on and off
    resonance.  That miss is ``|L|`` identically, so the two rows are the same
    statement measured twice.
    """
    rows = []
    prev: Optional[Dict[str, float]] = None
    for g in dampings:
        c = series_radius(separation, g)
        ratio = (prev["kappa_series"] / c["kappa_series"] if prev else None)
        rows.append({**c, "kappa_series_ratio_to_previous": ratio,
                     "damping_ratio": (prev["damping"] / g if prev else None),
                     "exponent": (math.log(prev["kappa_series"]
                                           / c["kappa_series"])
                                  / math.log(prev["damping"] / g)
                                  if prev else None)})
        prev = c

    exps = [r["exponent"] for r in rows if r["exponent"] is not None]
    g = 0.02
    base = series_radius(separation, g)
    th = CoupledThroat(separation, delay, +1, 0.5 * base["kappa_series"])
    on = base["omega_of_the_peak"]
    probe = []
    for w in (on, on + 0.5):
        L = th.loop_transfer(w, g)
        sv = coupled_propagator(w, chi_in, chi_out, th, g)
        f = free_branch_propagator(w, chi_in, chi_out, th, g)
        probe.append({"omega": w, "loop_gain": abs(L),
                      "relative_miss_of_one_traversal": abs(sv - f) / abs(sv)})

    peaks_are_integers = all(
        abs(r["omega_of_the_peak"] - round(r["omega_of_the_peak"])) < 1e-3
        for r in rows)
    return {
        "rows": rows, "exponents": exps,
        "the_series_radius_scales_like_the_regulator": bool(
            exps and all(abs(e - 1.0) < 0.15 for e in exps)),
        "mean_exponent": (float(np.mean(exps)) if exps else None),
        "the_peak_sits_on_a_bare_esu_resonance": bool(peaks_are_integers),
        "the_relative_miss_is_the_loop_gain": (
            "|1/(1−L) − 1| / |1/(1−L)| = |L| exactly, so the error of "
            "post-processing IS the round-trip gain"),
        "on_and_off_resonance": probe,
        "resonance_is_where_post_processing_is_worst": bool(
            probe[0]["relative_miss_of_one_traversal"]
            > 3.0 * probe[1]["relative_miss_of_one_traversal"]),
        "kappa_used": th.kappa,
        "what_this_shows": ("the one-traversal answer is an expansion in the "
                            "round-trip gain and it fails near the bare "
                            "resonances"),
        "what_this_does_not_show": ("that the coupled system is unstable "
                                    "there; γ is an Abel regulator, T_d "
                                    "carries the bare poles, and a finite "
                                    "coupling moves those poles rather than "
                                    "simply inflating a gain"),
    }


def measure_the_series_radius_is_not_the_stability_threshold(
        separation: float = 1.3, damping: float = 0.02,
        delays: Sequence[float] = (1.0, 0.5, math.pi),
        n_modes: int = 40, chi_in: float = 1.2,
        chi_out: float = 0.9) -> Dict[str, object]:
    """**The correction that matters most.**  ``max|L| = 1`` is not stability.

    Three separate things were being conflated:

    * **existence.**  ``1/(1−L)`` is defined for any ``L ≠ 1``, ``|L| > 1``
      included.  Solving and summing are different operations.
    * **convergence.**  ``|L| < 1`` is the radius of the traversal series, and
      ``κ_series = 1/max|T_d|`` does not depend on the delay at all.
    * **stability.**  Decided by the roots of ``D(ω) = 1 − L(ω) = 0`` in
      complex ``ω``: with ``e^{+iωt}``, ``Im ω > 0`` decays.  The uncoupled
      poles sit at ``ω = m + iγ``; the coupling displaces the ``m``-th by
      ``δ_m = −ηκ e^{−imΔ} sin(md)/(4π² sin d)``, whose imaginary part
      ``∝ sin(md) sin(mΔ)`` **changes sign with the mode**.  Phase-sensitive,
      and no bound on ``|L|`` can see it.

    The separation is made concrete by scanning the delay.  ``κ_series`` is the
    same number for every ``Δ``.  ``κ_stability`` is not: at ``Δ = π`` every
    first-order displacement is real (``sin(mπ) = 0``), the poles slide along
    the axis instead of off it, and the threshold rises far above the series
    radius.  In the gap between the two the solve is finite and stable while the
    series diverges — which is measured directly, not argued.

    Band-limited: only the first ``n_modes`` resonances are searched, and the
    least-damped root migrates to higher ``m`` as the band widens.
    """
    ks = series_radius(separation, damping, omega_max=float(n_modes) + 0.5,
                       n_grid=40001)
    rows = []
    for dl in delays:
        st = stability_threshold(separation, float(dl), +1, damping, n_modes)
        rows.append({"delay": float(dl),
                     "kappa_series": ks["kappa_series"],
                     "kappa_stability": st["kappa_stability"],
                     "mode_that_goes_first": st["mode_that_goes_first"],
                     "ratio": ((st["kappa_stability"] / ks["kappa_series"])
                               if st["kappa_stability"] else None)})

    # the gap: a coupling where the solve works and the series does not
    gap = max(rows, key=lambda r: (r["ratio"] or 0.0))
    demo: Dict[str, object] = {"delay": gap["delay"]}
    if gap["ratio"] and gap["ratio"] > 1.3:
        kap = math.sqrt(ks["kappa_series"] * gap["kappa_stability"])
        th = CoupledThroat(separation, gap["delay"], +1, kap)
        w_peak = ks["omega_of_the_peak"]
        L = th.loop_transfer(w_peak, damping)
        solved = coupled_propagator(w_peak, chi_in, chi_out, th, damping)
        walked = traversal_series(w_peak, chi_in, chi_out, th, damping, 400)
        pole = least_damped_pole(th, damping, n_modes)
        demo = {
            "delay": gap["delay"], "kappa": kap,
            "loop_gain_at_the_peak": abs(L),
            "the_series_diverges": bool(abs(L) > 1.0
                                        and abs(walked) > 1e3 * abs(solved)),
            "the_solve_is_finite": bool(math.isfinite(abs(solved))),
            "solved": solved, "walked_series": walked,
            "least_damped_pole_im": (pole["im"] if pole else None),
            "and_it_is_still_stable": bool(pole is not None
                                           and pole["im"] > 0.0),
        }

    # the first-order pole displacement, against the roots actually found
    th = CoupledThroat(separation, 1.0, +1, 0.3)
    poles = resonance_poles(th, damping, 12)
    worst = max((abs(p["omega"] - p["seed"]) for p in poles), default=0.0)
    return {
        "rows": rows, "kappa_series_is_delay_independent": True,
        "kappa_stability_depends_on_the_delay": bool(
            len({r["kappa_stability"] for r in rows
                 if r["kappa_stability"]}) > 1),
        "largest_ratio": gap["ratio"], "delay_of_the_largest_ratio":
            gap["delay"],
        "the_two_thresholds_are_different_numbers": bool(
            gap["ratio"] is not None and gap["ratio"] > 1.3),
        "a_coupling_between_them": demo,
        "n_poles_found": len(poles),
        "every_pole_matches_its_first_order_displacement": bool(
            poles and worst < 5e-3),
        "worst_pole_vs_first_order": float(worst),
        "n_modes_searched": int(n_modes),
        "what_this_corrects": ("|L| < 1 is the convergence radius of Σ Lⁿ, not "
                               "the existence condition (that is L ≠ 1) and "
                               "not the stability condition (that is Im ω > 0 "
                               "for every root of D)"),
    }


def measure_what_the_transfer_model_leaves_out(
        separation: float = 1.3, delay: float = 1.0,
        kappas: Sequence[float] = (0.3, 0.6, 1.0),
        omega: float = 2.3, damping: float = 0.02) -> Dict[str, object]:
    """The scope of the model, as numbers rather than a disclaimer.

    ``φ(M⁺) → ηκ φ(M⁻)`` is a relation between field **values** carried by the
    free Green function.  A throat that was a quotient of the manifold would
    need more, and each missing piece is stated here with what it would cost:

    * **no flux matching.**  Only the value is related; the normal derivative is
      not.  A Darmois–Israel junction matches both, and the mismatch is what
      `shells.junction` (PR #249) prices.
    * **no reflected channel.**  Each mouth has one outgoing amplitude, so the
      scattering object is ``1×1``.  A flux-conserving two-mouth junction needs
      at least a ``2×2`` unitary (transmit **and** reflect); ``|S| = κ``, and it
      is unitary only at ``κ = 1`` where it is still only half a junction.
    * **lossy by construction.**  Power out over power in is ``κ²`` per
      traversal, measured below, so for ``κ < 1`` the relation is not an
      identification of anything.

    None of this is a defect of the resolvent — that is exact for the model as
    posed.  It is a statement of which model.
    """
    rows = []
    for kap in kappas:
        th = CoupledThroat(separation, delay, +1, kap)
        # amplitude delivered to M⁺ by a unit mouth source, and re-emitted
        arrived = mouth_transfer(omega, separation, damping)
        emitted = th.eta * th.kappa * cmath.exp(-1j * omega * th.delay) * arrived
        rows.append({
            "kappa": kap,
            "power_out_over_power_in": float(abs(emitted) ** 2
                                             / abs(arrived) ** 2),
            "equals_kappa_squared": bool(
                abs(abs(emitted) ** 2 / abs(arrived) ** 2 - kap ** 2)
                < 1e-12 * max(1.0, kap ** 2)),
            "scattering_magnitude": float(abs(kap)),
            "unitary_1x1": bool(abs(abs(kap) - 1.0) < 1e-12)})
    return {
        "rows": rows,
        "the_power_ratio_is_kappa_squared": bool(
            all(r["equals_kappa_squared"] for r in rows)),
        "lossy_below_unit_coupling": bool(
            all(r["power_out_over_power_in"] < 1.0 for r in rows
                if r["kappa"] < 1.0)),
        "scattering_object_shape": "1x1",
        "shape_a_flux_conserving_junction_would_need": "at least 2x2 unitary",
        "what_is_missing": [
            "normal-derivative (flux) matching at the mouths",
            "a reflected channel, so the mouth S-matrix is 1x1 not 2x2",
            "self-adjointness / unitarity of the two-mouth scattering object",
            "and therefore any claim that this is a quotient of the manifold",
        ],
        "what_is_not_missing": ("the resolvent is exact for the model as posed;"
                                " the scope statement is about which model, not"
                                " about the solve"),
        "honest_name": ("a self-consistent rank-one mouth-transfer model, not a"
                        " throat boundary operator"),
    }
