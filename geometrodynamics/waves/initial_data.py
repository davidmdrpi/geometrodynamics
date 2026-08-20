"""Initial-data geometry on the resolved neck: the signed areal response.

Where this sits
───────────────
PR #263 asked whether ``A + B`` produces a metric response rescaling ``A`` or
``B`` cannot, measured it in the transverse-traceless sector, and then — asked
for a *geometric* verdict — proved that sector cannot give one:

    ``δA/A = −½⟨h_nn⟩``, and that average vanishes **identically** for any
    transverse-traceless field.  Tracelessness kills the isotropic part of
    ``⟨n_i n_j f(k·n̂)⟩`` and transversality kills the rest.  Building more
    harmonics adds more exact zeros.

So a signed ``δA/A`` has to come from the **scalar** sector.  This module does
that as an **initial-data** problem rather than an evolution, which is what
removes the two things PR #263 avoided the scalar sector for.

The constraint
──────────────
On a slice of the Einstein static universe the background has ``K̄_ij = 0``, and
the ``K`` terms in the Hamiltonian constraint are *quadratic*, so at the order
the field's stress tensor enters,

    ``δR⁽³⁾ = 16πG δρ``

with no time derivatives, **no sound speed, and no Eddington mode** — those are
properties of the *evolution*, and a constraint does not have them.  The fluid
is held rigid (``δρ_fluid = 0``), which is consistent precisely because the
scalar's stress tensor is separately conserved.

Writing the spatial metric as ``g_ij = ψ⁴ĝ_ij`` with ``ψ = 1 + u``:

    ``δR⁽³⁾ = −4u R̂ − 8∇̂²u``      and      ``δA/A = 4u``

so with ``R̂ = 6`` on the unit sphere the constraint is

    ``∇²u + 3u = −2πG δρ`` .

The conformal ansatz is **not** a restriction here.  A transverse-traceless
piece contributes nothing to ``δR⁽³⁾`` (both terms of ``−∇²(tr h) + ∇^i∇^j
h_ij`` vanish), and a longitudinal piece is a diffeomorphism of a constant-``R̂``
background, so it contributes nothing either.  The conformal factor is the whole
of what this equation sees, and ``4u`` is the whole of the areal response.

Why the resolved neck is required, twice
────────────────────────────────────────
Neither reason is a refinement.  Both are the difference between a finite answer
and no answer.

**The source diverges at the mouths.**  The traceless projection of PR #263
survived the ``1/χ⁴`` singularities because the angular average killed them;
an *energy* density has no such cancellation.  Measured: the interference
density ``ΔT₀₀`` is bounded at a **source** — it tends to ``≈ 0.1`` as
``r → 0``, so the radial integrand vanishes — but at a **mouth** it goes as
``1/χ⁴`` with the integrand *growing* into the singularity, ``1.1e-03`` at
``r = 0.32`` rising to ``1.4e-01`` at ``r = 0.01``.  Both waves drive the same
mouths, so their mouth terms multiply.  Removing the balls removes the
singularity from the domain.

**And the operator is exactly degenerate on the whole sphere.**  ``∇²u = −3u``
is ``−∇²u + u = 4u``, and ``4 = (n+1)²`` at ``n = 1``: the constraint operator
sits **exactly on** the ESU's dipole level.  On ``S³`` it has a genuine kernel —
the four harmonics ``x^A`` — and PR #261's fixed-ambient Green function has a
literal pole there, which is not a numerical difficulty but the statement that
the problem is unsolvable as posed.

Removing the balls lifts it, and the lifting has **exact closed forms — four of
them**, because the ``k = 1`` multiplet has four members.  About a mouth's
centre it splits into one ``ℓ = 0`` member ``ψ = cos χ`` and three ``ℓ = 1``
members ``ψ = sin χ``, and both are degenerate at ``λ = 4``:

    ``N₀(a, 4)  =  +4π sin²a · tan a  ⟶  +4π a³``
    ``N₁(a, 4)  =  −4π sin²a · cot a  =  −2π sin 2a  ⟶  −4π a``

verified against the shooting solve to ``1e-09`` and ``2e-10`` respectively.
**The signs differ and so do the orders.**  The monopole stiffness is positive
and ``O(a³)``; the three dipole stiffnesses are **negative** and ``O(a)`` — in
those channels the exterior *assists* rather than resists.  Only ``ℓ ≥ 2`` is
genuinely non-degenerate, staying ``O(1)`` and positive.

A first version of this module reported the ``ℓ = 1`` channels as positive and
``O(1)``, i.e. not degenerate at all.  That was an artefact of a bug in
`neck`: the radial equation carried ``ℓ(ℓ+2)`` — the ``S³`` harmonic eigenvalue
— where the angular Laplacian on the **two**-sphere of directions gives
``ℓ(ℓ+1)``.  ``ℓ = 0`` is untouched by it, which is why every closed-form check
passed and the error survived a merge.  **The three dipole partners of the
degeneracy were hidden by the wrong centrifugal term.**

This matters for the solve rather than only for bookkeeping.  The solvability
condition ``∫δρ Y₁ dV = 0`` has **four** components, one per kernel member, and
a reduced monopole model addresses exactly one of them.  The other three have to
be absorbed by the ``ℓ = 1`` channels at the mouths — whose stiffness is
negative, so they absorb with the opposite sign.

What is still put in
────────────────────
The fluid is **rigid**: ``δρ_fluid = 0``.  That is one named assumption on the
front page rather than an unnamed sound speed, and it is consistent because the
scalar's stress tensor is separately conserved.  A responsive fluid would change
the number; how much is a question for the evolution problem, not this one.

The **single-wave** pieces still need resolved sources — ``T₀₀[A]`` alone goes
as ``1/χ⁴`` at its own source and its integral diverges, exactly as PR #260's
point mouth did.  The interference piece, which is what ``Δh`` is built from,
does not.
"""

from __future__ import annotations

import math
from typing import Dict, List, Sequence, Tuple

import numpy as np

from .backreaction import (
    BackreactionSetup, _direction_rule, _tangent_basis, solve_batch,
    stress_series,
)
from .finite_mouth import FOUR_PI
from .neck import exterior_dtn, exterior_dtn_monopole, exterior_log_derivative
from .two_wave import orthonormal_frame

__all__ = [
    "CONSTRAINT_EIGENVALUE",
    "KERNEL_PROJECTOR",
    "areal_response_from_conformal_factor",
    "kernel_component",
    "regularised_green",
    "constraint_operator_eigenvalue",
    "measure_the_constraint_operator_is_degenerate_on_the_sphere",
    "measure_removing_the_balls_lifts_the_degeneracy",
    "measure_the_interference_source_needs_the_resolved_neck",
    "radial_energy_profile",
]

#: ``∇²u = −3u`` is ``−∇²u + u = λu`` at ``λ = 4`` — which is ``(n+1)²`` at
#: ``n = 1``, the ESU's dipole level.  The coincidence is the whole reason this
#: calculation needs a resolved neck.
CONSTRAINT_EIGENVALUE = 4.0


def constraint_operator_eigenvalue(scalar_curvature: float = 6.0) -> float:
    """``λ`` such that the constraint operator is `neck`'s ``−∇² + 1 − λ``.

    The Hamiltonian constraint gives ``∇²u + (R̂/2)u = −2πG δρ``, so
    ``∇²u = −(R̂/2)u`` and ``λ = 1 + R̂/2``.  On the unit sphere ``R̂ = 6`` and
    ``λ = 4``, which is exactly ``(n+1)²`` at ``n = 1``.
    """
    return 1.0 + 0.5 * float(scalar_curvature)


def areal_response_from_conformal_factor(u) -> np.ndarray:
    """``δA/A = 4u`` — areas scale as ``ψ⁴``.

    Sign convention, since the whole round turns on it: ``u < 0`` means the
    conformal factor shrinks, so the mouth's area falls and the geometry moves
    **toward** a neck.  ``u > 0`` moves away from one.
    """
    return 4.0 * np.asarray(u, dtype=float)


KERNEL_PROJECTOR = 2.0 / math.pi ** 2


def regularised_green(chi) -> np.ndarray:
    """The constraint operator's Green function, with the ``k = 1`` pole removed.

    ``G(χ,ω) = sin(ω(π−χ))/(4π sin χ sin πω)`` has a **literal pole** at
    ``ω = 2`` — that is the degeneracy, not a numerical difficulty.  Expanding
    about it with ``ω = 2 + ε`` and ``u = π − χ``:

        ``sin(ωu) ≈ −sin 2χ + ε u cos 2χ`` ,   ``sin πω = sin πε ≈ πε``

    so the pole part is ``−cos χ / (2π²ε)`` — residue exactly the ``k = 1``
    harmonic about the source point, as it has to be — and the finite remainder
    is

        ``G_⊥(χ) = (π − χ) cos(2χ) / (4π² sin χ)``

    which satisfies

        ``(∇² + 3) G_⊥ = −δ³(x−y) + (2/π²) cos χ`` .

    The inhomogeneous term is the **normalised kernel projector**: the ``k = 1``
    harmonics are ``x^A`` with ``‖x^A‖² = π²/2`` and ``Σ_A x^A y^A = cos χ``, so
    the projector kernel is ``2 cos χ/π²`` — measured as a constant ``0.202642``
    across ``χ ∈ [0.3, 2.9]``.  So ``G_⊥`` inverts the operator on the kernel's
    orthogonal complement, which is the most that can be asked of it.

    Two limits check the construction: ``G_⊥ · 4πχ → 1`` at short distance, the
    correct three-dimensional normalisation; and ``G_⊥(π) = 1/(4π²)`` exactly,
    finite at the antipode, which is the ``k = 1`` removal doing its work.
    """
    x = np.asarray(chi, dtype=float)
    return (math.pi - x) * np.cos(2.0 * x) / (4.0 * math.pi ** 2 * np.sin(x))


def kernel_component(values, weights, points) -> np.ndarray:
    """``∫ f x^A dV`` — the source's overlap with the operator's kernel.

    Nonzero means the free problem is **unsolvable**: the four numbers here are
    the obstruction that the mouths have to absorb.  Reported rather than
    assumed away, because a solvability condition that is quietly violated
    produces a plausible answer to a question that was never well posed.
    """
    f = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    return np.einsum("p,p,pA->A", w, f, np.asarray(points, dtype=float))


# ════════════════════════════════════════════════════════════════════════════
# MEASUREMENTS
# ════════════════════════════════════════════════════════════════════════════
def measure_the_constraint_operator_is_degenerate_on_the_sphere(
        k_max: int = 6) -> Dict[str, object]:
    """``∇² + 3`` has an exact kernel on ``S³``, and it is not an accident.

    Scalar harmonics on the unit ``S³`` obey ``∇²Y_k = −k(k+2)Y_k``, so the
    constraint operator has eigenvalues ``3 − k(k+2)``:

        ``k = 0`` → ``+3``,  ``k = 1`` → **``0``**,  ``k = 2`` → ``−5``, …

    Two things follow, and both matter.  The operator is **indefinite** — so
    this is not a nice positive elliptic solve — and it has a **kernel** at
    ``k = 1``, which imposes a solvability condition ``∫δρ Y₁ dV = 0`` on the
    source.

    The degeneracy is exactly the ESU's own dipole level: ``λ = 1 + R̂/2 = 4``
    is ``(n+1)²`` at ``n = 1``.  PR #261's fixed-ambient Green function has a
    literal pole there — the model does not merely lose accuracy at the
    constraint eigenvalue, it is undefined at it.
    """
    rows = []
    for k in range(int(k_max) + 1):
        lap = -k * (k + 2)
        rows.append({"k": int(k), "laplacian_eigenvalue": float(lap),
                     "operator_eigenvalue": float(lap + 3.0),
                     "free_esu_lambda": float((k + 1) ** 2)})
    kernel = [r["k"] for r in rows if abs(r["operator_eigenvalue"]) < 1e-12]
    negative = [r["k"] for r in rows if r["operator_eigenvalue"] < -1e-12]
    return {
        "rows": rows,
        "constraint_eigenvalue": CONSTRAINT_EIGENVALUE,
        "derived_eigenvalue": constraint_operator_eigenvalue(),
        "kernel_modes": kernel,
        "negative_modes": negative,
        "the_operator_has_a_kernel": bool(kernel == [1]),
        "the_operator_is_indefinite": bool(len(negative) > 0),
        "it_sits_on_the_esu_dipole_level": bool(
            abs(constraint_operator_eigenvalue() - 4.0) < 1e-12
            and abs((1 + 1) ** 2 - 4.0) < 1e-12),
        "why_it_matters": ("the kernel imposes a solvability condition on the "
                           "source, and the fixed-ambient Green function has a "
                           "literal pole at this eigenvalue -- the point-mouth "
                           "model is not inaccurate here, it is undefined"),
    }


def measure_removing_the_balls_lifts_the_degeneracy(
        radii: Sequence[float] = (0.01, 0.02, 0.05, 0.1, 0.2, 0.35),
        ells: Sequence[int] = (0, 1, 2, 3)) -> Dict[str, object]:
    """And the lifting has a closed form: ``N₀(a, 4) = 4π a³``.

    On ``S³ ∖ B_a`` the exterior stiffness at the degenerate eigenvalue is not
    zero, so the operator is invertible and the solve is well posed.  The
    monopole member of the ``k = 1`` multiplet about the mouth's centre is
    ``cos χ``, whose logarithmic derivative at ``χ = a`` is ``−tan a``; hence

        ``N₀ = −4π sin²a · (ψ'/ψ) = 4π sin²a · tan a  ⟶  4π a³``

    which is exact, and reproduces the measured stiffness to ``1e-09`` at every
    radius including where the ``4πa³`` form has started to drift.  The small-``a``
    coefficient is a constant ``12.56637 = 4π`` across a decade and a half in
    ``a``, against ``4πa`` at ``λ = 0`` — so the ratio is ``a²``, and that extra
    ``a²`` is the degeneracy itself.

    **The multiplet has four members and they do not behave alike.**  The three
    ``ℓ = 1`` partners are ``ψ = sin χ``, giving
    ``N₁ = −4π sin²a · cot a = −2π sin 2a ⟶ −4π a``: **negative**, and only
    ``O(a)``.  So the dipole channels are softer than the monopole by two orders
    in ``a`` *and* of the opposite sign — the exterior assists there instead of
    resisting.  Only ``ℓ ≥ 2`` is genuinely non-degenerate.  It vanishes as ``a → 0``: **the
    degeneracy returns in the point limit**, which is the precise sense in which
    this calculation needs a resolved neck rather than merely preferring one.

    The higher multipoles are not degenerate — they stay ``O(1)`` — so the
    softness is confined to the one channel the coincidence produces.
    """
    rows = []
    for a in radii:
        n4 = exterior_dtn_monopole(float(a), CONSTRAINT_EIGENVALUE)
        n0 = exterior_dtn_monopole(float(a), 0.0)
        rows.append({"radius": float(a), "N0_at_the_eigenvalue": float(n4),
                     "over_a_cubed": float(n4 / float(a) ** 3),
                     "N0_at_zero": float(n0),
                     "ratio": float(n4 / n0)})
    multipoles = []
    for a in (0.05, 0.2):
        multipoles.append({
            "radius": float(a),
            "by_ell": [float(exterior_dtn(float(a), CONSTRAINT_EIGENVALUE,
                                          int(l))) for l in ells]})
    coeffs = [r["over_a_cubed"] for r in rows]
    # the closed form is exact; what is measured is the integrator converging to
    # it, and the dipole is the harder solve -- v = sin^2 chi starts as e^2 at
    # the antipode against a singular centrifugal term, so the error grows as the
    # ball shrinks.  Reported as a trend rather than pinned to a magic tolerance.
    dipole_errors = []
    for a in radii:
        exact = -2.0 * math.pi * math.sin(2.0 * float(a))
        got = exterior_dtn(float(a), CONSTRAINT_EIGENVALUE, 1)
        dipole_errors.append(abs(got - exact) / abs(exact))
    shot = [-(FOUR_PI * math.sin(a) ** 2)
            * exterior_log_derivative(a, CONSTRAINT_EIGENVALUE, 0)
            for a in radii]
    worst_shoot = max(abs(s - r["N0_at_the_eigenvalue"])
                      / abs(r["N0_at_the_eigenvalue"])
                      for s, r in zip(shot, rows))
    return {
        "rows": rows,
        "multipoles": multipoles,
        "four_pi": FOUR_PI,
        "coefficient_spread": float(max(coeffs) - min(coeffs)),
        "worst_shooting_vs_closed_form": float(worst_shoot),
        "it_is_nonzero_everywhere": bool(
            all(r["N0_at_the_eigenvalue"] > 0.0 for r in rows)),
        "it_is_four_pi_a_cubed": bool(
            all(abs(c - FOUR_PI) / FOUR_PI < 2e-3 for c in coeffs)),
        "it_vanishes_in_the_point_limit": bool(
            rows[0]["N0_at_the_eigenvalue"] < 1e-4),
        "the_dipole_partners_are_also_degenerate": bool(
            all(m["by_ell"][1] < 0.0 for m in multipoles)),
        "worst_dipole_closed_form_error": float(max(dipole_errors)),
        "dipole_error_falls_with_radius": bool(
            dipole_errors[0] > 10.0 * dipole_errors[-1]),
        "the_dipole_closed_form_holds": bool(max(dipole_errors) < 1e-5),
        "only_ell_two_and_above_is_free": bool(
            all(m["by_ell"][2] > 0.5 and m["by_ell"][3] > 0.5
                for m in multipoles)),
        "what_it_settles": ("the resolved neck makes the constraint operator "
                            "invertible and the point limit makes it singular "
                            "again, so this is a requirement rather than a "
                            "refinement"),
    }


def radial_energy_profile(centre: Sequence[float],
                          radii: Sequence[float] = (0.01, 0.02, 0.04, 0.08,
                                                    0.16, 0.32),
                          setup: BackreactionSetup = None,
                          window: Tuple[float, float] = (4.0, 30.0),
                          n_theta: int = 10, n_phi: int = 20
                          ) -> List[Dict[str, float]]:
    """Angular-averaged ``T₀₀`` at a set of radii about a point.

    The radial integrand of ``∫T₀₀ dV`` is ``⟨T₀₀⟩ sin²r``: constant means the
    integral converges, growing into the singularity means it does not.  This is
    the diagnostic that decides whether a source can be used at all.
    """
    s = setup or BackreactionSetup()
    c = np.asarray(centre, dtype=float)
    c = c / np.linalg.norm(c)
    basis = _tangent_basis(c)
    dirs, wd = _direction_rule(int(n_theta), int(n_phi))
    w = wd / wd.sum()
    times = s.grid.times
    sl = (times > window[0]) & (times < window[1])
    out = []
    for r in radii:
        pts = math.cos(r) * c[None, :] + math.sin(r) * (dirs @ basis)
        frames = np.array([orthonormal_frame(p) for p in pts])
        base = s._setup(pts[0])
        a = solve_batch(base, s.source_a, base.pulse_a, pts, frames)
        b = solve_batch(base, s.source_b, base.pulse_b, pts, frames)
        total = {k: a[k] + b[k] for k in a}
        ta, tb = stress_series(a), stress_series(b)
        tt = stress_series(total)
        row = {"radius": float(r)}
        for name, x in (("single", ta), ("cross", tt - ta - tb)):
            e = np.einsum("p,pt->t", w, x[..., 0, 0])
            rms = float(np.sqrt(np.mean(e[sl] ** 2)))
            row[name] = rms
            row[f"{name}_integrand"] = rms * math.sin(r) ** 2
        out.append(row)
    return out


def measure_the_interference_source_needs_the_resolved_neck(
        setup: BackreactionSetup = None) -> Dict[str, object]:
    """The source converges at a source and diverges at a mouth.

    PR #263's traceless projection survived the ``1/χ⁴`` singularities because
    the angular average killed them exactly.  An **energy** density has no such
    cancellation, so the diagnostic has to be run again for this round's source
    and it gives a different answer in the two places.

    * at a **source**, ``ΔT₀₀`` is bounded — it tends to ``≈ 0.1`` as
      ``r → 0``, so the radial integrand ``⟨ΔT₀₀⟩ sin²r`` vanishes and the
      integral converges.  (``T₀₀`` for a *single* wave does not: it goes as
      ``1/χ⁴`` and its integral diverges, which is the point-source self-energy
      and is why only the interference piece is computable here.)
    * at a **mouth**, ``ΔT₀₀`` goes as ``1/χ⁴`` and the integrand **grows** into
      the singularity.  Both waves drive the same mouths, so their mouth terms
      multiply and the cancellation that saved the source is not available.

    So the resolved neck is not an improvement on the point-mouth model for this
    calculation.  It is the difference between an integral that exists and one
    that does not.
    """
    s = setup or BackreactionSetup()
    c1, _ = s.mouths()
    at_source = radial_energy_profile(s.source_a, setup=s)
    at_mouth = radial_energy_profile(c1, setup=s)

    def inward_growth(rows, key):
        """How much the radial integrand grows going **inward** from its
        minimum.

        Comparing the innermost radius to the outermost is the wrong
        statistic and it passed for the wrong reason in a first draft: the
        outermost shell is large simply because its volume is, so a divergent
        integrand can still look smaller there.  What divergence means is that
        the integrand **turns around and grows as the singularity is
        approached**, so the comparison has to be against the minimum.
        """
        vals = [r[key] for r in rows]
        return float(vals[0] / max(min(vals), 1e-300))

    return {
        "at_the_source": at_source,
        "at_the_mouth": at_mouth,
        "cross_inward_growth_at_the_source": inward_growth(at_source,
                                                           "cross_integrand"),
        "cross_inward_growth_at_the_mouth": inward_growth(at_mouth,
                                                          "cross_integrand"),
        "single_inward_growth_at_the_source": inward_growth(at_source,
                                                            "single_integrand"),
        "the_cross_term_converges_at_the_source": bool(
            inward_growth(at_source, "cross_integrand") < 1.5),
        "the_cross_term_diverges_at_the_mouth": bool(
            inward_growth(at_mouth, "cross_integrand") > 10.0),
        "the_single_wave_term_diverges_at_its_own_source": bool(
            inward_growth(at_source, "single_integrand") > 3.0),
        "why_the_neck_is_required": ("both waves drive the same mouths, so "
                                     "their mouth terms multiply in the cross "
                                     "term and the angular cancellation that "
                                     "saved PR #263's traceless projection is "
                                     "not available for an energy density"),
    }
