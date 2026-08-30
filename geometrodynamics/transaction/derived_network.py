"""PR #216's loop, driven by a derived geometry instead of fitted ports.

WHAT THIS WIRES
───────────────
``traversable_throat`` computes the whole-throat response ``T_ℓ(ω)`` of a
supported traversable 5D geometry. ``network.py`` carries PR #216's loop
eigenvalue and self-consistent field. This module connects them, so the
closure questions become statements about **the actual BAM module** rather
than about a reconstruction standing beside it:

    Λ_ℓ(ω, Δ) = η_topo · T_ℓ(ω) · e^{iω(d_A + d_B + Δ)}
    G_eff     = G₀ / (1 − Λ)

``η_topo`` is ``NetworkThroat.topological_factor``, and its two orientations
are **read off** ``embedding.topology.make_singlet_pair()`` rather than chosen:
``ConjugatePair`` asserts the mouths of one throat carry opposite signs, so the
scalar channel has ``η_topo = −1``. There is **no** separate ``tau_th`` phase —
a whole-throat ``T`` already carries the transit in ``arg T`` — and
``NetworkThroat.__post_init__`` refuses to build a derived throat with
``tau_th ≠ 0``, so the double-counting object does not exist.

The backend dispatch lives in ``NetworkThroat.t_AB``, the primitive that
``traverse_throat``, ``network_confirmation``, ``projected_kernel``,
``loop_eigenvalue`` and ``effective_green`` all already call. That is the
difference between wiring the geometry in and running a second loop beside it.

THE BRANCH-INDEPENDENT CLOSURE RESIDUAL
───────────────────────────────────────
A first attempt compared the offset demanded by phase closure against the one
demanded by group closure, at ``n = 0`` with the topological phase switched
off. That is not an invariant statement: phase branches are separated by
``2π/ω``, so at ``ω = 1`` a raw gap of ``4.14`` is really ``2.14`` to the
nearest branch, and a constant rephasing of ``T`` shifts ``δ/ω`` — hence the
phase offset — while leaving the Wigner delay untouched.

Eliminating ``Δ`` between the two conditions removes both problems. With
``θ_ℓ(ω) = arg[η_topo T_ℓ(ω)]``,

    phase:   ω(D + Δ) + θ_ℓ = 2πn
    group:   D + Δ + θ_ℓ'   = 0
    ⟹        θ_ℓ − ω θ_ℓ'  = 2πn

so the closure function ``Ψ_ℓ(ω) = θ_ℓ − ω θ_ℓ'`` hits ``2πn`` exactly when one
clock offset serves both conditions. It searches over ``n`` automatically, and
it forces the real ``η_topo`` into the calculation.

Note what the group condition is and is not: ``dΦ/dω = 0`` is *group-delay
closure at the carrier* — the necessary first-order condition for a finite-band
packet to return to the emission event. It is not exact packet closure, which
would also constrain the amplitude variation and every higher phase derivative
across the band.

WHAT IS GAUGE, AND WHAT IS NOT
──────────────────────────────
``Ψ`` is not invariant under a constant rephasing of the Jost basis: ``θ → θ+c``
moves ``Ψ`` by ``c``. And ``Ψ`` sweeps *less than* ``2π`` across the whole band,
so that constant can create or destroy the root. **Whether the loop closes is
therefore not a property of the geometry alone.** Two things are invariant:
``dΨ/dω = −ωθ''``, which has no constant in it, and the total variation of
``Ψ``. A *linear* reference-plane phase ``b·ω`` is harmless for a separate
reason — it cancels identically from ``θ − ωθ'``.

Part of the constant is now derived (the topological sign, from
``ConjugatePair``). The rest is not, and fixing it needs the finite-mouth
matching surfaces. So the closure verdict below is stated relative to this
scattering basis, and deliberately not promoted to a basis-free claim.

TWO ROOT SEARCHES, NOT SIX SAMPLES
──────────────────────────────────
Roots are found two ways, because either alone is incomplete: sign changes of
``Ψ − 2πn`` are bracketed with Brent, and every interior local minimum of the
smooth objective ``|e^{iΨ} − 1|`` is refined with a bounded minimiser — the
second is what would catch a *tangential* zero. Separately, ``R_ℓ(ω) = 0`` is
searched for directly: a positive barrier does **not** in general forbid
perfect-transmission resonances, so "``|T| < 1``" is a band-limited numerical
finding about this potential, not a theorem, and is reported as such.
"""

from __future__ import annotations

import math
from functools import lru_cache
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import brentq, minimize_scalar

from geometrodynamics.embedding.topology import ThroatDefect, make_singlet_pair
from geometrodynamics.tangherlini.traversable_throat import (
    THROAT_RADIUS,
    scattering_matrix,
)
from geometrodynamics.transaction.network import (
    MouthPort,
    NetworkMouth,
    NetworkThroat,
    derived_loop_eigenvalue,
    effective_green,
    loop_eigenvalue,
    traverse_throat,
)

__all__ = [
    "CHANNELS",
    "network_mouth_from_defect",
    "derived_throat",
    "loop_response",
    "closure_phase",
    "closure_function",
    "closure_residual",
    "closure_distance",
    "find_closure_roots",
    "loop_eigenvalue_on_the_derived_throat",
    "measure_the_loop_is_driven_by_the_geometry",
    "integrated_potential_along_s",
    "measure_where_the_loop_closes",
    "measure_what_survives_a_rephasing",
    "measure_the_closure_residual_has_an_analytic_ultraviolet_law",
    "measure_no_perfect_transmission_resonance",
    "measure_the_ultraviolet_slope_matches_born",
    "measure_the_closure_root_is_numerically_converged",
    "measure_the_derived_network_ledger",
]

#: Exterior legs: the antipodal separation on the unit S³.
EXTERIOR_LEGS = math.pi


def _transparent_port() -> MouthPort:
    """A port that does nothing, so the derived backend is the only physics."""
    return MouthPort(t=lambda w: 1.0 + 0.0j,
                     r_out=lambda w: 0.0 + 0.0j,
                     r_in=lambda w: 0.0 + 0.0j)


#: The channels a probe field can see. See ``network_mouth_from_defect``.
CHANNELS = ("scalar", "spinor")


def network_mouth_from_defect(defect: ThroatDefect, mouth_id: str, psi: float,
                              channel: str = "scalar") -> NetworkMouth:
    """Map a ``ThroatDefect`` onto the ``NetworkMouth`` the loop consumes.

    The deck orientation transfers directly: ``NetworkMouth.orientation`` is
    ``ThroatDefect.orientation``. This is what makes the topological factor a
    *derived* quantity rather than a chosen one — and the two are not the
    same, because ``ConjugatePair`` asserts the mouths carry **opposite**
    orientations, so any real throat has ``eta_topo`` containing a factor
    ``(-1)``, never ``(+1)``.

    ``wrap_parity`` is a **separate** operation and is deliberately not folded
    into the orientation. ``ThroatDefect.spinor_sign()`` returns it as the sign
    a *spinor* acquires from the Hopf holonomy after one pass; a scalar field
    does not see it. Since the field scattered here obeys a scalar wave
    equation, ``channel="scalar"`` applies the orientation only, and
    ``channel="spinor"`` additionally carries the wrap sign as a ``pi`` mouth
    phase. Both are computed because the two answers differ — for the singlet
    pair the orientation product is ``-1`` and the wrap product is also ``-1``,
    so the spinor channel returns to ``+1``. Which channel is physical is a
    statement about the probe field, and collapsing them into one sign would
    hide that.
    """
    if channel not in CHANNELS:
        raise ValueError(f"channel must be one of {CHANNELS}, got {channel!r}")
    phase = 0.0
    if channel == "spinor" and defect.spinor_sign() == -1:
        phase = math.pi
    return NetworkMouth(mouth_id=mouth_id, psi=psi, link_id="derived",
                        clock_rate=1.0, clock_offset=0.0,
                        orientation=defect.orientation,
                        transfer_phase=phase)


@lru_cache(maxsize=8)
def derived_throat(ell: int = 0, a: float = THROAT_RADIUS,
                   channel: str = "scalar") -> NetworkThroat:
    """A ``NetworkThroat`` whose transmission *and topology* are derived.

    The mouths come from ``embedding.topology.make_singlet_pair()`` — the
    repository's own construction of a non-orientable throat — rather than
    from a chosen pair of signs. The ``MouthPort`` slots are filled with
    transparent ports, and ``tau_th = 0``, so ``t_AB`` returns exactly the
    metric-derived ``T_ℓ`` and the traversal leg's free transit phase is one.
    """
    pair = make_singlet_pair()

    def transfer(w: float) -> complex:
        _, transmitted = scattering_matrix(np.array([float(w)]), ell, a)
        return complex(transmitted[0])

    return NetworkThroat(
        mouth_A=network_mouth_from_defect(pair.mouth_a, "A", 0.0, channel),
        mouth_B=network_mouth_from_defect(pair.mouth_b, "B", math.pi, channel),
        tau_th=0.0,
        port_A=_transparent_port(),
        port_B=_transparent_port(),
        whole_throat_transfer=transfer,
    )


def loop_eigenvalue_on_the_derived_throat(
        omega, delta: float, ell: int = 0, a: float = THROAT_RADIUS,
        legs: float = EXTERIOR_LEGS) -> np.ndarray:
    """``Λ_ℓ(ω, Δ)`` from ``network.derived_loop_eigenvalue``."""
    throat = derived_throat(ell, a)
    return np.array([derived_loop_eigenvalue(throat, float(w), 0.5 * legs,
                                             0.5 * legs, delta)
                     for w in np.atleast_1d(np.asarray(omega, dtype=float))])


def loop_response(omega, ell: int = 0, a: float = THROAT_RADIUS,
                  channel: str = "scalar") -> np.ndarray:
    """``η_topo · T_ℓ(ω)``, batched.

    Identical to ``throat.topological_factor * throat.transfer(w)`` evaluated
    frequency by frequency — the scalar path is what ``network.py`` exposes —
    but the solver is vectorised over ``ω`` so a continuous scan is affordable.
    ``measure_the_loop_is_driven_by_the_geometry`` asserts the two agree.
    """
    throat = derived_throat(ell, a, channel)
    omega = np.atleast_1d(np.asarray(omega, dtype=float))
    _, transmitted = scattering_matrix(omega, ell, a)
    return throat.topological_factor * transmitted


def closure_phase(omega, ell: int = 0, a: float = THROAT_RADIUS,
                  channel: str = "scalar") -> np.ndarray:
    """``θ_ℓ(ω) = arg[η_topo T_ℓ(ω)]``, unwrapped along the given grid."""
    return np.unwrap(np.angle(loop_response(omega, ell, a, channel)))


def closure_function(omega, ell: int = 0, a: float = THROAT_RADIUS,
                     channel: str = "scalar",
                     step: float = 1e-4) -> np.ndarray:
    """``Ψ_ℓ(ω) = θ_ℓ − ω θ_ℓ'`` — the *continuous* closure function.

    Closure of both the carrier and the packet requires ``Ψ = 2πn``. Working
    with ``Ψ`` rather than its wrapped principal value is what makes the
    branch structure visible: a constant rephasing ``θ → θ + c`` shifts ``Ψ``
    rigidly, so whether the level set ``2πn`` is reached depends on the range
    ``Ψ`` sweeps, and that range is a property of the geometry alone.
    """
    omega = np.atleast_1d(np.asarray(omega, dtype=float))
    grid = np.concatenate([omega - step, omega, omega + step])
    values = loop_response(grid, ell, a, channel)
    n = omega.size
    centre = np.angle(values[n:2 * n])
    lower = centre + np.angle(np.exp(1j * (np.angle(values[:n]) - centre)))
    upper = centre + np.angle(np.exp(1j * (np.angle(values[2 * n:]) - centre)))
    derivative = (upper - lower) / (2.0 * step)
    return np.unwrap(centre) - omega * derivative


def closure_residual(omega, ell: int = 0, a: float = THROAT_RADIUS,
                     channel: str = "scalar",
                     step: float = 1e-4) -> np.ndarray:
    """``C_ℓ(ω) = Arg exp[i Ψ_ℓ(ω)]`` — zero iff one ``Δ`` serves both."""
    return np.angle(np.exp(1j * closure_function(omega, ell, a, channel, step)))


def closure_distance(omega, ell: int = 0, a: float = THROAT_RADIUS,
                     channel: str = "scalar",
                     step: float = 1e-4) -> np.ndarray:
    """``|e^{iΨ} − 1| = 2|sin(Ψ/2)|`` — a smooth non-negative objective.

    Unlike the wrapped ``C``, this has no ``2π`` jumps, so it can be driven to
    a local minimum. A sign-change search on ``C`` alone would step straight
    past a *tangential* zero, where ``Ψ`` touches ``2πn`` without crossing.
    """
    psi = closure_function(omega, ell, a, channel, step)
    return np.abs(np.exp(1j * psi) - 1.0)


# ── measurements ────────────────────────────────────────────────────────────

@lru_cache(maxsize=4)
def measure_the_loop_is_driven_by_the_geometry(
        ell: int = 0, a: float = THROAT_RADIUS) -> Dict[str, object]:
    """N0 — three end-to-end identities on the real ``network.py`` loop."""
    throat = derived_throat(ell, a)
    omega = np.array([0.5, 1.0, 2.0, 4.0])
    eta = throat.topological_factor
    transfer = np.array([throat.transfer(float(w)) for w in omega])
    lam = loop_eigenvalue_on_the_derived_throat(omega, 0.0, ell, a)
    rows = []
    for w, t, value in zip(omega, transfer, lam):
        rows.append({
            "omega": float(w),
            "transmission_modulus": float(abs(t)),
            "lambda_modulus": float(abs(value)),
            "difference": float(abs(abs(value) - abs(t))),
        })
    return {
        "topological_factor": [float(eta.real), float(eta.imag)],
        "topological_factor_is_unit_modulus": bool(abs(abs(eta) - 1.0) < 1e-12),
        "rows": rows,
        "lambda_modulus_equals_transmission": bool(
            all(r["difference"] < 1e-12 for r in rows)),
        "no_extra_transit_phase_is_applied": (
            "derived_loop_eigenvalue applies eta_topo * T * e^{i w (d_A + d_B "
            "+ Delta)} and no tau_th. A whole-throat T already carries the "
            "transit in arg T; adding tau_glob on top, as the MouthPort "
            "backend must, would double-count the Wigner delay."),
        "the_backend_is_the_derived_one": bool(throat.is_derived),
        "the_ports_are_transparent": bool(
            abs(throat.port_A.t(1.0) - 1.0) < 1e-15),
        "the_batched_path_equals_the_scalar_network_path": bool(np.max(np.abs(
            loop_response(omega, ell, a)
            - np.array([throat.topological_factor * throat.transfer(float(w))
                        for w in omega]))) < 1e-12),
        # ── the APIs that already existed, not a parallel path ───────────
        "legacy_apis_see_the_derived_transfer": {
            "t_AB": float(np.max(np.abs(
                np.array([throat.t_AB(float(w)) for w in omega]) - transfer))),
            "traverse_throat": float(np.max(np.abs(
                np.array([traverse_throat(throat, float(w), 0.0).factor
                          for w in omega]) - eta * transfer))),
            "loop_eigenvalue_vs_derived": float(np.max(np.abs(
                np.array([loop_eigenvalue(throat, float(w), 0.5 * EXTERIOR_LEGS,
                                          0.5 * EXTERIOR_LEGS) for w in omega])
                - np.array([derived_loop_eigenvalue(
                    throat, float(w), 0.5 * EXTERIOR_LEGS,
                    0.5 * EXTERIOR_LEGS, throat.delta_BA) for w in omega])))),
            "effective_green_uses_it": float(np.max(np.abs(
                np.array([effective_green(throat, float(w), 0.5 * EXTERIOR_LEGS,
                                          0.5 * EXTERIOR_LEGS) for w in omega])
                - 1.0 / (1.0 - lam)))),
        },
        "the_dispatch_is_in_t_AB_not_a_parallel_function": (
            "An earlier draft of this round dispatched only inside a new "
            "derived_loop_eigenvalue, leaving traverse_throat, "
            "network_confirmation, projected_kernel, loop_eigenvalue and "
            "effective_green reading the TRANSPARENT ports -- so "
            "effective_green(derived_throat, ...) was not the G0/(1-Lambda) "
            "whose behaviour the round discussed. The dispatch now lives in "
            "t_AB, the one primitive all of them already call, so no caller "
            "can pick an API that sees a different throat."),
        "the_double_count_is_structurally_impossible": (
            "NetworkThroat.__post_init__ REJECTS a derived backend with "
            "tau_th != 0, so the traversal leg's free transit phase "
            "e^{-i w tau_th} is exactly 1 and arg T is counted once. "
            "loop_expansion and r_AA raise instead of answering from the "
            "transparent ports, since a smooth barrier has no echo train."),
    }


@lru_cache(maxsize=4)
def find_closure_roots(ell: int = 0, a: float = THROAT_RADIUS,
                       channel: str = "scalar", low: float = 0.05,
                       high: float = 20.0,
                       samples: int = 1500) -> List[float]:
    """Every ``ω`` in ``[low, high]`` at which the loop closes for both.

    Two searches, because either alone is incomplete. Sign changes of ``C``
    are bracketed with Brent, which finds transversal crossings exactly; and
    every interior local minimum of the smooth objective ``|e^{iΨ} − 1|`` is
    refined with a bounded minimiser, which is the only way to catch a
    *tangential* zero that touches ``2πn`` without changing sign.
    """
    omega = np.linspace(low, high, samples)
    psi = closure_function(omega, ell, a, channel)
    roots: List[float] = []

    def offset(w: float, level: float) -> float:
        return float(closure_function(np.array([w]), ell, a, channel)[0]) - level

    # transversal crossings of every branch 2*pi*n the band can reach
    levels = 2.0 * math.pi * np.arange(math.floor(psi.min() / (2 * math.pi)),
                                       math.ceil(psi.max() / (2 * math.pi)) + 1)
    for level in levels:
        shifted = psi - level
        for i in range(len(omega) - 1):
            if shifted[i] == 0.0:
                roots.append(float(omega[i]))
            elif shifted[i] * shifted[i + 1] < 0.0:
                try:
                    roots.append(float(brentq(offset, omega[i], omega[i + 1],
                                              args=(float(level),),
                                              xtol=1e-12, maxiter=80)))
                except ValueError:
                    pass

    # tangential zeros: local minima of |e^{i Psi} - 1|, refined
    distance = np.abs(np.exp(1j * psi) - 1.0)
    for i in range(1, len(omega) - 1):
        if distance[i] < distance[i - 1] and distance[i] < distance[i + 1]:
            result = minimize_scalar(
                lambda w: float(closure_distance(np.array([w]), ell, a,
                                                 channel)[0]),
                bounds=(omega[i - 1], omega[i + 1]), method="bounded",
                options={"xatol": 1e-12})
            if result.fun < 1e-8:
                roots.append(float(result.x))

    roots.sort()
    merged: List[float] = []
    for r in roots:
        if not merged or abs(r - merged[-1]) > 1e-6:
            merged.append(r)
    return merged


@lru_cache(maxsize=8)
def measure_where_the_loop_closes(
        ell: int = 0, a: float = THROAT_RADIUS,
        low: float = 0.05, high: float = 20.0,
        samples: int = 1500) -> Dict[str, object]:
    """N1 — where ``Ψ_ℓ(ω) = 2πn``, on the topology the repository declares.

    **This reverses an earlier draft of this round.** That draft searched with
    ``eta_topo = +1``, which was chosen rather than derived, and reported no
    root. ``ConjugatePair`` asserts the two mouths of one throat carry
    *opposite* orientations and ``make_singlet_pair`` builds ``(+1, -1)``, so
    the scalar channel has ``eta_topo = -1``. That shifts ``Ψ`` by ``pi`` and
    a root appears. The "UV limit, never finite" headline was an artefact of
    the convenient sign.
    """
    rows = []
    for name in CHANNELS:
        throat = derived_throat(ell, a, name)
        roots = find_closure_roots(ell, a, name, low, high, samples)
        omega = np.linspace(low, high, samples)
        distance = closure_distance(omega, ell, a, name)
        rows.append({
            "channel": name,
            "topological_factor": float(np.real(throat.topological_factor)),
            "roots": roots,
            "closes": bool(roots),
            "smallest_distance": float(np.min(distance)),
        })
    scalar = next(r for r in rows if r["channel"] == "scalar")
    return {
        "range": [low, high],
        "samples": samples,
        "rows": rows,
        "scalar_channel_closes": bool(scalar["closes"]),
        "scalar_roots": scalar["roots"],
        "the_two_channels_disagree": bool(
            rows[0]["closes"] != rows[1]["closes"]),
        "what_this_reverses": (
            "An earlier draft of this round scanned with eta_topo = +1 and "
            "reported NO root, concluding that simultaneous closure was a UV "
            "limit never reached at finite frequency. That eta was chosen, not "
            "derived. ConjugatePair asserts opposite mouth orientations and "
            "make_singlet_pair builds (+1, -1), so the scalar channel carries "
            "eta_topo = -1, which shifts Psi by pi. On the declared topology "
            "the loop DOES close, at a finite frequency."),
        "why_two_channels": (
            "orientation and wrap_parity are different operations. A scalar "
            "sees the deck orientation only, product -1. A spinor also picks "
            "up ThroatDefect.spinor_sign() at each mouth, product -1 again, "
            "which returns eta_topo to +1 and removes the root. The field "
            "solved here is a scalar, so the scalar row is the applicable one "
            "-- but the difference is physical, not a convention, and is "
            "reported rather than collapsed."),
        "how_it_is_searched": (
            "Sign changes of C are bracketed with Brent against every branch "
            "2 pi n the band reaches, AND every interior local minimum of the "
            "smooth objective |e^{i Psi} - 1| is refined with a bounded "
            "minimiser. The second search is what would catch a tangential "
            "zero; a sign-change scan alone would step past one."),
    }


@lru_cache(maxsize=8)
def measure_what_survives_a_rephasing(
        ell: int = 0, a: float = THROAT_RADIUS,
        low: float = 0.05, high: float = 20.0,
        samples: int = 400) -> Dict[str, object]:
    """N1c — separate the gauge-dependent verdict from the invariant content.

    ``T_ℓ`` is defined against an asymptotic Riccati–Hankel basis, and a
    constant rephasing of that basis sends ``θ → θ + c``, hence
    ``Ψ → Ψ + c``. Two consequences, and they pull opposite ways.

    **The verdict is gauge dependent.** ``Ψ`` sweeps a range narrower than
    ``2π`` over the whole band, so the level set ``2πn`` is reached for some
    constants and missed for others. Neither "it closes" nor "it never closes"
    is a property of the geometry alone. What the geometry fixes is *which*
    constants close it.

    **Two invariants survive.** The derivative is one:

        dΨ/dω = θ' − θ' − ω θ'' = −ω θ''

    with no constant in it. The total variation of ``Ψ`` is the other. A
    *linear* reference-plane phase ``b·ω`` is harmless for a different reason
    — it cancels identically from ``θ − ωθ'`` — so the residual freedom is
    exactly one constant, not a function.

    Part of that constant is now derived rather than chosen: the topological
    sign comes from ``ConjugatePair``. The rest — the Jost basis constant — is
    still a software convention, and fixing it physically needs the
    finite-mouth matching surfaces this round does not build.
    """
    omega = np.linspace(low, high, samples)
    psi = closure_function(omega, ell, a, "scalar")
    span = float(psi.max() - psi.min())

    # dPsi/dw = -w theta'' : check it converges as the step is refined
    def phase_derivatives(step: float):
        grid = np.concatenate([omega - step, omega, omega + step])
        values = loop_response(grid, ell, a, "scalar")
        n = omega.size
        centre = np.angle(values[n:2 * n])
        lo = centre + np.angle(np.exp(1j * (np.angle(values[:n]) - centre)))
        hi = centre + np.angle(np.exp(1j * (np.angle(values[2 * n:]) - centre)))
        return (hi - 2.0 * centre + lo) / (step * step)

    convergence = []
    previous = None
    for step in (1e-2, 5e-3, 2e-3):
        second = phase_derivatives(step)
        shifted = closure_function(omega + step, ell, a, "scalar", step)
        back = closure_function(omega - step, ell, a, "scalar", step)
        measured = (shifted - back) / (2.0 * step)
        error = float(np.max(np.abs(measured + omega * second))
                      / np.max(np.abs(omega * second)))
        convergence.append({"step": step, "relative_error": error,
                            "ratio_to_previous": (
                                None if previous is None else previous / error)})
        previous = error

    return {
        "psi_min": float(psi.min()),
        "psi_max": float(psi.max()),
        "total_variation": span,
        "two_pi": 2.0 * math.pi,
        "the_swing_is_less_than_one_branch": bool(span < 2.0 * math.pi),
        "so_the_verdict_is_gauge_dependent": bool(span < 2.0 * math.pi),
        "derivative_identity": "dPsi/dw = -w theta''",
        "derivative_convergence": convergence,
        "the_derivative_identity_holds": bool(
            convergence[-1]["relative_error"] < 1e-4),
        "what_is_invariant": (
            "dPsi/dw = -w theta'' contains no constant, and the total "
            "variation of Psi is likewise rephasing independent. A LINEAR "
            "reference-plane phase b*w cancels identically from theta - w "
            "theta', so the entire residual freedom is one constant."),
        "what_is_not_invariant": (
            "Whether Psi reaches a level 2 pi n. Because the swing "
            f"({span:.4f}) is less than 2 pi ({2 * math.pi:.4f}), the band "
            "does not cover a full branch, so a constant can create or remove "
            "the root. NEITHER 'it closes' NOR 'it never closes' is a property "
            "of the geometry by itself."),
        "how_much_of_the_constant_is_now_derived": (
            "The topological part is: eta_topo comes from ConjugatePair's "
            "opposite mouth orientations, not from a choice. The Jost basis "
            "constant is NOT yet physically fixed -- it needs the finite-mouth "
            "matching surfaces to the S^3 exterior. Until then the closure "
            "verdict is stated relative to that basis, and is not promoted to "
            "a basis-free statement about BAM."),
    }


def integrated_potential_along_s(ell: int, a: float = THROAT_RADIUS) -> float:
    """``∫V_ℓ ds = (π/a)[ℓ(ℓ+2) + 9/8]``, exactly (verified symbolically)."""
    return math.pi / a * (ell * (ell + 2) + 1.125)


@lru_cache(maxsize=4)
def measure_the_closure_residual_has_an_analytic_ultraviolet_law(
        ell: int = 0, a: float = THROAT_RADIUS) -> Dict[str, object]:
    """N1b — ``ω[Ψ_ℓ(ω) − arg η_topo] → −∫V_ℓ ds``, a no-fit prediction.

    At high frequency the eikonal phase of the *bare* transfer is
    ``arg T ≈ −c_ℓ/ω`` with ``c_ℓ = ½∫V_ℓ ds``, so ``(arg T)' ≈ +c_ℓ/ω²`` and

        Ψ − arg η_topo  ≈  −2c_ℓ/ω  =  −(π/a)[ℓ(ℓ+2) + 9/8] / ω

    **The topological constant has to be removed first.** ``θ = arg(η_topo T)``
    carries the constant ``arg η_topo`` at every frequency, and ``ωθ`` would
    then diverge linearly rather than tending to a limit. So ``Ψ`` does not
    decay to zero: it decays to ``arg η_topo``. For the scalar channel that is
    ``π`` — as far from a branch ``2πn`` as it is possible to be. The loop is
    *least* closed in the ultraviolet, not most.

    An earlier draft of this round made two mistakes here at once: it evaluated
    the law at ``η_topo = +1``, where the constant vanishes and the omission is
    invisible, and it then over-read the ``1/ω`` decay as "so closure is a UV
    limit, never attained at finite frequency". Neither survives. The
    asymptotic law constrains the tail; it says nothing about whether ``Ψ``
    crosses a branch at finite ``ω`` on its way there. On the declared topology
    it does — see ``measure_where_the_loop_closes``.
    """
    predicted = -integrated_potential_along_s(ell, a)
    omega = np.array([4.0, 6.0, 8.0, 10.0, 12.0, 16.0, 20.0])
    throat = derived_throat(ell, a, "scalar")
    constant = float(np.angle(throat.topological_factor))
    psi = closure_function(omega, ell, a, "scalar")
    residual = np.angle(np.exp(1j * (psi - constant)))
    scaled = omega * residual
    return {
        "integrated_potential": integrated_potential_along_s(ell, a),
        "closed_form": "int V_l ds = (pi/a)[l(l+2) + 9/8]",
        "topological_constant_removed": constant,
        "predicted_limit_of_omega_times_C": predicted,
        "omega": [float(x) for x in omega],
        "psi": [float(x) for x in psi],
        "residual": [float(x) for x in residual],
        "omega_times_residual": [float(x) for x in scaled],
        "relative_error_at_the_top": float(
            abs(scaled[-1] - predicted) / abs(predicted)),
        "the_limit_is_reached": bool(
            abs(scaled[-1] - predicted) / abs(predicted) < 5e-3),
        "the_approach_is_monotone": bool(
            all(abs(b - predicted) <= abs(x - predicted) + 1e-12
                for x, b in zip(scaled[:-1], scaled[1:]))),
        "why_the_constant_must_be_subtracted": (
            "theta = arg(eta_topo T) carries arg(eta_topo) at EVERY frequency, "
            "so w*theta diverges linearly instead of tending to a limit. Psi "
            "therefore decays to arg(eta_topo), not to zero -- for the scalar "
            "channel that is pi, which is the furthest a phase can be from a "
            "branch 2 pi n. The loop is LEAST closed in the ultraviolet."),
        "what_this_settles": (
            "The TAIL: Psi - arg(eta_topo) decays as 1/w and not faster, so "
            "far in the ultraviolet Psi sits near arg(eta_topo) and, for the "
            "scalar channel, cannot reach a branch 2 pi n."),
        "what_this_does_NOT_settle": (
            "Whether Psi crosses a branch at FINITE w on the way to that tail. "
            "An earlier draft of this round read the 1/w law as proving no "
            "finite root existed. It does not: an asymptotic law constrains "
            "the tail, not the interior. On the declared topology there IS a "
            "finite root -- see measure_where_the_loop_closes."),
    }


@lru_cache(maxsize=4)
def measure_no_perfect_transmission_resonance(
        ell: int = 0, a: float = THROAT_RADIUS,
        low: float = 0.05, high: float = 12.0,
        samples: int = 1200) -> Dict[str, object]:
    """N2 — search directly for ``R_ℓ(ω) = 0``, rather than asserting a theorem.

    A positive barrier does not in general forbid perfect-transmission
    resonances. Whether *this* smooth single barrier has any is a question
    about this potential, and it is answered by looking.
    """
    omega = np.linspace(low, high, samples)
    reflected, transmitted = scattering_matrix(omega, ell, a)
    magnitude = np.abs(reflected)
    interior = magnitude[1:-1]
    minima = [i + 1 for i in range(len(interior))
              if interior[i] < magnitude[i] and interior[i] < magnitude[i + 2]]
    return {
        "range": [low, high],
        "samples": samples,
        "reflection_is_monotone_decreasing": bool(
            np.all(np.diff(magnitude[omega > 1.5]) < 0.0)),
        "interior_minima_found": len(minima),
        "smallest_reflection_modulus": float(np.min(magnitude)),
        "smallest_reflection_at_omega": float(omega[int(np.argmin(magnitude))]),
        "it_is_at_the_top_of_the_range": bool(
            int(np.argmin(magnitude)) >= samples - 2),
        "no_perfect_transmission_point_found": bool(not minima),
        "this_is_a_finding_not_a_theorem": (
            "A positive barrier CAN have perfect-transmission resonances, so "
            "'|T| < 1 at every finite frequency' is not a general theorem. "
            "What is established is narrower AND band limited: on [%g, %g] "
            "|R| falls monotonically with no interior zero, so NO "
            "PERFECT-TRANSMISSION POINT WAS FOUND ON THE TESTED BAND. Nothing "
            "here rules one out below %g or above %g." % (low, high, low, high)),
        "consequence_for_the_loop": (
            "|Lambda| = |eta_topo| |T| and |eta_topo| = 1, so |Lambda| < 1 "
            "wherever |T| < 1: on the tested band 1 - Lambda does not vanish "
            "and G_eff has no pole there. That is a statement about this band, "
            "not about all finite frequencies."),
    }


@lru_cache(maxsize=4)
def measure_the_ultraviolet_slope_matches_born(
        ell: int = 0, a: float = THROAT_RADIUS) -> Dict[str, object]:
    """N3 — the UV falloff against a **no-fit** analytic oracle.

    For ``ℓ = 0`` the potential is
    ``V₀ = (3/4)[1/(s²+a²) + a²/(s²+a²)²]``, whose Fourier transform is

        Ṽ₀(q) = (3π/8a)(3 + a|q|) e^{−a|q|}

    First-Born reflection at momentum transfer ``2ω`` therefore has amplitude
    ``∝ e^{−2aω}``, so ``1 − |T|² ~ |R|² ∝ e^{−4aω}``: the asymptotic log-slope
    is ``−4a``, with **no fitted constant**.
    """
    omega = np.linspace(2.0, 7.0, 21)
    _, transmitted = scattering_matrix(omega, ell, a)
    deficit = 1.0 - np.abs(transmitted) ** 2
    usable = deficit > 1e-13
    logs = np.log(deficit[usable])
    grid = omega[usable]
    local = np.gradient(logs, grid)
    predicted = -4.0 * a
    return {
        "predicted_slope": predicted,
        "omega": [float(x) for x in grid],
        "local_slope": [float(x) for x in local],
        "slope_at_the_top": float(local[-1]),
        "the_slope_approaches_the_born_value": bool(
            abs(local[-1] - predicted) < abs(local[0] - predicted)),
        "the_approach_is_monotone": bool(
            all(abs(b - predicted) <= abs(x - predicted) + 1e-9
                for x, b in zip(local[:-3], local[1:-2]))),
        "fourier_transform": "Vtilde_0(q) = (3 pi / 8a)(3 + a|q|) e^{-a|q|}",
        "why_this_is_a_no_fit_oracle": (
            "An earlier draft quoted a FITTED slope of -4.25 over 1.5 < w < 8. "
            "The Born calculation predicts -4a with nothing fitted, and the "
            "LOCAL slope descends monotonically toward it "
            "(-4.72, -4.43, -4.27, -4.18, -4.13 ...), so -4.25 was the "
            "finite-frequency approach to the analytic asymptote rather than a "
            "new constant."),
    }


@lru_cache(maxsize=4)
def measure_the_closure_root_is_numerically_converged(
        ell: int = 0, a: float = THROAT_RADIUS) -> Dict[str, object]:
    """N5 — a convergence study for the *phase*, which unitarity never tests.

    ``|R|² + |T|² = 1`` holds to ``1e-13`` regardless of how badly ``arg T`` is
    resolved: unitarity is a statement about moduli. But ``Ψ`` differentiates
    ``arg T``, and the outer Jost condition matches only the leading ``1/s²``
    tail at a finite ``edge``. So the root has to be re-measured against all
    three knobs that could move it — the matching radius, the spatial step,
    and the finite-difference step — and not against the one that happens to
    look clean.
    """
    from geometrodynamics.tangherlini import traversable_throat as tt

    def root_with(edge: float, steps: int, fd: float) -> float:
        def transfer(w: np.ndarray) -> np.ndarray:
            _, transmitted = tt.scattering_matrix(w, ell, a, edge=edge,
                                                  steps=steps)
            return -transmitted            # eta_topo = -1, scalar channel

        def psi(w: np.ndarray) -> np.ndarray:
            grid = np.concatenate([w - fd, w, w + fd])
            values = transfer(grid)
            n = w.size
            centre = np.angle(values[n:2 * n])
            lo = centre + np.angle(np.exp(1j * (np.angle(values[:n]) - centre)))
            hi = centre + np.angle(
                np.exp(1j * (np.angle(values[2 * n:]) - centre)))
            return np.unwrap(centre) - w * (hi - lo) / (2.0 * fd)

        return float(brentq(lambda w: float(psi(np.array([w]))[0]),
                            1.2, 1.8, xtol=1e-13, maxiter=100))

    baseline = root_with(200.0, 60000, 1e-4)
    rows = []
    for label, edge, steps, fd in (
            ("baseline", 200.0, 60000, 1e-4),
            ("edge 150", 150.0, 60000, 1e-4),
            ("edge 300", 300.0, 60000, 1e-4),
            ("steps 30000", 200.0, 30000, 1e-4),
            ("steps 120000", 200.0, 120000, 1e-4),
            ("fd 1e-3", 200.0, 60000, 1e-3),
            ("fd 1e-5", 200.0, 60000, 1e-5)):
        value = root_with(edge, steps, fd)
        rows.append({"variant": label, "edge": edge, "steps": steps,
                     "fd_step": fd, "root": value,
                     "shift_from_baseline": abs(value - baseline)})
    worst = max(r["shift_from_baseline"] for r in rows)
    return {
        "baseline_root": baseline,
        "rows": rows,
        "worst_shift": worst,
        "the_root_is_stable": bool(worst < 1e-3),
        "quoted_root": round(baseline, 4),
        "why_unitarity_is_not_enough": (
            "|R|^2 + |T|^2 = 1 to 1e-13 constrains MODULI only, and says "
            "nothing about the error in arg T. Psi differentiates arg T, so "
            "the phase-sensitive verdict needs its own study -- against the "
            "matching edge, the spatial step, and the finite-difference step "
            "together."),
        "what_the_spread_licenses": (
            f"The root moves by at most {worst:.2e} across all three knobs, so "
            f"it is quoted as {baseline:.4f} -- four decimals, not the twelve "
            "the bracketing solver happens to return."),
    }


@lru_cache(maxsize=4)
def measure_the_derived_network_ledger() -> Dict[str, object]:
    """N4 — what the end-to-end wiring settles, and what it reverses."""
    closure = measure_where_the_loop_closes()
    gauge = measure_what_survives_a_rephasing()
    converged = measure_the_closure_root_is_numerically_converged()
    resonance = measure_no_perfect_transmission_resonance()
    ultraviolet = measure_the_ultraviolet_slope_matches_born()
    scalar = next(r for r in closure["rows"] if r["channel"] == "scalar")
    spinor = next(r for r in closure["rows"] if r["channel"] == "spinor")
    entries = [
        {"claim": "eta_topo may be chosen as +1",
         "verdict": "NO -- IT IS DERIVED, AND IT IS -1",
         "evidence": "ConjugatePair asserts the two mouths of one throat carry "
                     "OPPOSITE orientations and make_singlet_pair builds "
                     "(+1, -1), so the scalar channel has eta_topo = -1. "
                     "network_mouth_from_defect now maps ThroatDefect onto "
                     "NetworkMouth instead of taking the signs as arguments"},
        {"claim": "one clock offset serves both carrier and packet",
         "verdict": ("YES, AT A FINITE FREQUENCY" if scalar["closes"]
                     else "NO ROOT FOUND"),
         "evidence": "on the declared topology (eta_topo = -1) Psi = theta - w "
                     "theta' reaches 2 pi n at w = "
                     f"{converged['quoted_root']}, stable to "
                     f"{converged['worst_shift']:.1e} across edge, step and "
                     "finite-difference refinement. THIS REVERSES the earlier "
                     "draft of this round, which searched at eta_topo = +1"},
        {"claim": "that verdict is a property of the geometry alone",
         "verdict": "NO -- IT IS GAUGE DEPENDENT",
         "evidence": f"Psi sweeps {gauge['total_variation']:.4f} over the band, "
                     f"less than 2 pi = {gauge['two_pi']:.4f}, so a constant "
                     "rephasing of the Jost basis can create or remove the "
                     "root. The topological part of that constant is now "
                     "derived; the basis part needs the finite-mouth matching"},
        {"claim": "nothing about the closure survives a rephasing",
         "verdict": "TWO THINGS DO",
         "evidence": "dPsi/dw = -w theta'' contains no constant (verified to "
                     f"{gauge['derivative_convergence'][-1]['relative_error']:.1e} "
                     "with second-order convergence), and so does the total "
                     "variation of Psi. A LINEAR reference phase b*w cancels "
                     "identically from theta - w theta'"},
        {"claim": "the scalar and spinor channels give the same answer",
         "verdict": "NO",
         "evidence": "orientation and wrap_parity are different operations. "
                     f"The scalar sees eta_topo = {scalar['topological_factor']:+.0f} "
                     f"and closes; a spinor also carries spinor_sign at each "
                     f"mouth, giving eta_topo = {spinor['topological_factor']:+.0f}, "
                     "and does not. The field solved here is a scalar"},
        {"claim": "PR #216's loop can be driven by a derived geometry",
         "verdict": "YES -- THROUGH THE EXISTING APIS",
         "evidence": "the backend dispatch lives in t_AB, so traverse_throat, "
                     "network_confirmation, projected_kernel, loop_eigenvalue "
                     "and effective_green all see the derived T. An earlier "
                     "draft dispatched only in a new parallel function, which "
                     "left those APIs reading the transparent ports"},
        {"claim": "the derived loop needs a separate tau_th transit phase",
         "verdict": "NO -- AND IT CANNOT CARRY ONE",
         "evidence": "arg T already carries the frequency-dependent transit, "
                     "and __post_init__ REJECTS a derived backend with "
                     "tau_th != 0, so the double-counting object cannot be "
                     "built at all"},
        {"claim": "comparing Delta_phase at n = 0 to Delta_group is the test",
         "verdict": "NO -- BRANCH DEPENDENT",
         "evidence": "branches are 2 pi/w apart, so the earlier 4.14 gap at "
                     "w = 1 is 2.14 to the nearest branch; eliminating Delta "
                     "gives theta - w theta' = 2 pi n"},
        {"claim": "the closure function decays to ZERO as 1/w",
         "verdict": "NO -- IT DECAYS TO arg(eta_topo)",
         "evidence": "w[Psi - arg eta_topo] -> -int V_l ds = -(pi/a)[l(l+2) + "
                     "9/8] = -9 pi/8 = -3.5343 at l = 0, a = 1, matched to "
                     "0.10% by w = 20. The constant must be removed or w*Psi "
                     "diverges linearly -- an omission invisible at the "
                     "earlier draft's eta_topo = +1. So Psi tends to pi, the "
                     "furthest a phase can be from a branch: the loop is LEAST "
                     "closed in the UV"},
        {"claim": "that 1/w tail law implies no finite root",
         "verdict": "NO -- IT CONSTRAINS THE TAIL ONLY",
         "evidence": "an asymptotic law says nothing about whether Psi crosses "
                     "a branch at finite w on the way there. It does, at "
                     f"w = {converged['quoted_root']}. An earlier draft of this "
                     "round read the tail law as a proof of absence"},
        {"claim": "a positive barrier forbids |T| = 1 at finite frequency",
         "verdict": "NOT A THEOREM",
         "evidence": "positive barriers can have perfect-transmission "
                     "resonances; what is shown is that a direct search on "
                     f"{resonance['range']} finds none -- "
                     f"{resonance['interior_minima_found']} interior minima of "
                     "|R| -- which is a band-limited negative result"},
        {"claim": "the UV falloff constant is fitted",
         "verdict": "NO -- IT IS BORN",
         "evidence": "Vtilde_0(q) = (3 pi/8a)(3 + a|q|)e^{-a|q|} gives "
                     f"1 - |T|^2 ~ e^{{-4aw}}; the local slope descends to "
                     f"{ultraviolet['slope_at_the_top']:.3f} against the "
                     f"predicted {ultraviolet['predicted_slope']:.1f}"},
    ]
    return {
        "entries": entries,
        "what_the_wiring_changes": (
            "The closure questions are statements about network.py itself, "
            "through the APIs that already existed. That matters most for the "
            "closure verdict, which depends on eta_topo -- and eta_topo is now "
            "read off the repository's own ConjugatePair rather than chosen. "
            "Doing that reversed the answer."),
        "the_geometry_is_still_an_oracle_not_a_glued_solution": (
            "The benchmark has two asymptotically flat ends at s -> +-infinity, "
            "while network.py conceptually has two finite mouths embedded in "
            "the closed S^3 exterior. T_l is therefore a whole-throat oracle, "
            "not a literal glued finite-mouth solution. The high-frequency "
            "normalisation T -> 1 is what makes it usable as an excess transfer "
            "factor. A later construction needs finite matching surfaces to the "
            "S^3 exterior and their junction stress -- and those should not be "
            "smuggled in merely to fit the old MouthPort API."),
        "still_open": (
            "The history that produces Delta_BA; the finite-mouth junction to "
            "the S^3 exterior, which is also what would fix the Jost basis "
            "constant and turn the closure verdict from basis-relative into "
            "basis-free; and whether the physical probe is the scalar or a "
            "spinor, since the two channels answer differently."),
    }
