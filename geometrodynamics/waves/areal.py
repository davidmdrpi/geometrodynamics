"""The mouth matching, and the signed ``ΔA/A`` it produces.

Where this sits
───────────────
`initial_data` established the two facts that make a signed areal response
possible and hard at the same time:

* the constraint reduces to ``∇²u + 3u = −2πG δρ`` with ``ΔA/A = 4u``, and
* that operator is **exactly degenerate** on ``S³`` — its kernel is the four
  ``k = 1`` harmonics ``x^A`` — so the equation has no solution on the closed
  sphere unless the source happens to be orthogonal to all four, and the
  interference source is not.

Removing the two mouth balls is what makes it solvable.  This module does the
removal properly: it writes the field as a source term plus layers on the two
mouth spheres plus a free kernel element, closes the system with the throat,
and reads off ``ΔA/A`` at each mouth.

The representation
──────────────────
With ``σ = −2πG δρ`` and ``G_⊥`` the pseudo-inverse Green function of
`initial_data.regularised_green`,

    ``u = U + Σ_j [ A_j G_⊥(χ_j) + D_j·𝒟_j ] + c·x`` ,

where ``U = −G_⊥[σ 1_Ω]``, ``𝒟_j`` is the dipole layer at mouth ``j`` (the
derivative of ``G_⊥`` as its source point is moved along a tangent direction),
``D_j ⊥ c_j``, and ``c ∈ R⁴`` is a free kernel element.  Twelve field
unknowns: two monopole strengths, six dipole components, four kernel
coefficients.  (`solve_matching` solves an ``18×18`` system: the extra six are
the tube's end amplitudes, carried rather than eliminated purely for
conditioning — see below.)

The solvability condition is not an extra assumption — it is what makes the
representation consistent.  ``G_⊥`` satisfies

    ``L G_⊥ = −δ + (2/π²) cos χ``

and the trailing term is **exactly** the projector onto the four-dimensional
kernel: ``Σ_A x^A y^A = cos χ(x,y)`` and ``‖cos χ‖²_{S³} = π²/2``, so the
projector's integral kernel is ``(2/π²) cos χ`` and nothing else.  Demanding
that the leftover ``k = 1`` pieces cancel on ``Ω`` gives, exactly,

    ``Σ_j A_j c_j + Σ_j D_j = S_σ`` ,     ``S_σ = ∫_Ω y σ(y) dV`` .

Four equations, and the remaining fourteen are the throat.  **Two monopoles cannot satisfy these four.**  ``A_1 c_1 + A_2 c_2``
sweeps out only the plane spanned by the two mouth positions, so any component
of ``S_σ`` off that plane has to be carried by the dipole layers.  For the
two-wave interference source, measured, the monopole-only condition fails by
``62.5%`` of the obstruction: without the ``ℓ = 1`` layers the problem has no
solution at all.

And then — this was not the expected outcome — those layers contribute
**nothing** to the answer.  The off-plane obstruction produces an areal
response of ``6e-17``, which is zero: the dipole layers deposit it in kernel
elements ``x²`` and ``x³`` that vanish at both mouths, because both mouths lie
in the ``(x⁰, x¹)`` plane.  A first draft of this docstring said the ``ℓ = 1``
sector *is* the calculation.  It is required for the calculation to exist and
invisible in its result, which is a better outcome than the draft claimed:
the ``ℓ = 1`` source moments are the least converged input by a wide margin —
they drift ``41%`` between quadrature levels where the obstruction drifts
``1.5%`` — and the signed answer does not depend on them.  Scaling them by
three, or replacing them with noise, moves ``ΔA/A`` by ``5e-04``.

A tube of cross-sectional area
``𝒜`` is a round ``S²`` of radius ``r = √(𝒜/4π)`` crossed with a line, so
``R̂ = 2/r² = 8π/𝒜`` and the same constraint reduction gives ``∇² + R̂/2``:

    ``ℓ_tube = 0``:  ``∂_s² + 4π/𝒜``   — oscillatory, wavenumber ``k = 1/r``
    ``ℓ_tube = 1``:  ``∂_s² − 4π/𝒜``   — evanescent, at the *same* rate

Continuity of ``u`` and of flux at each end closes it.  Note what the first
line means: **the tube is a resonant cavity for the constraint**, and at
``kL = nπ`` the response has a pole and the sign of ``ΔA/A`` flips.

That is not a remark.  It is the boundary of the result, and it was tested
rather than left as one.  The working throat carries ``𝒜 = 4π`` against a mouth
sphere of area ``4π sin²a`` — wider than its own mouths by a factor of ``400``
at ``a = 0.05``.  Set the two **equal**, so the tube is exactly as narrow as
the mouths it joins, and ``k = 1/sin a``: the same length ``0.9`` becomes
``kL/π = 5.73`` at ``a = 0.05`` and ``2.87`` at ``a = 0.10``, past five poles
and past two, ``4.6%`` of its own length from the next.  **The sign does not
survive.** At ``a = 0.05`` both mouths open; at ``a = 0.10`` they disagree.
`measure_the_sign_does_not_survive_a_matched_tube` is that calculation, and it
is the reason the headline is stated as *at the wide working throat, off
resonance* rather than as a property of the interference source.

What is checked and what is assumed
───────────────────────────────────
The whole assembly is checked against exact one-dimensional boundary-value
solves on the punctured sphere — one in each sector, because a check that only
exercises ``ℓ = 0`` is exactly how the ``ℓ(ℓ+2)`` error in `neck` survived a
full test suite.  Agreement is ``4e-10`` or better at every radius in both
sectors, and **flat** in the radius: it is the assembly's numerical floor, not
a truncation error, so no order of convergence is claimed from it.

Getting there took two corrections worth recording.

The first assembly agreed with the reference at ``1e-06`` and no better, at
every radius — and that number did not move when the reference's own
quadrature, stencil and boundary-value tolerances were each tightened by four
orders.  It was therefore not the reference's floor but a systematic error of
the assembly: the two-point stencil for the mouth-sphere radial derivative,
whose *relative* truncation error on a ``1/χ`` field is exactly ``step²``.  A
five-point rule removed four orders of magnitude.  **A discrepancy that refuses
to move when you refine the other side is the other side telling you it is not
the problem.**

The second was found only because the matched-tube check was asked for.  The
``ℓ = 1`` rows were originally a ``cosh``/``sinh`` transfer matrix, which costs
a condition number of ``e^{2κL}``.  At ``𝒜 = 4π`` that is ``e^{1.8} = 6``, and
invisible.  At the matched area it is ``e^{36} = 4.4e+15`` — the system is
singular to double precision, and the first matched-tube run reported a
condition number of ``2.9e+15`` and an answer anyway.  Carrying the tube's two
end amplitudes as unknowns instead of eliminating them never forms ``e^{+κL}``;
the system grows from ``12×12`` to ``18×18``, every coefficient is bounded by
one, and the conditioning falls to ``5.5e+07`` — and by ``1.5e+04`` from
``2.1e+05`` at the wide working throat too.  The reference solves reproduce to
the *same* ``4e-10``, digit for digit, so it is a change of form and not of
content.  **A model parameter moved by a factor of four hundred is not a
perturbation of a formulation, it is a test of one.**

Two things are modelled rather than derived.  The **gluing map** identifying
the two mouths' transverse frames through the tube is taken to be parallel
transport along the joining geodesic.  A handle may be glued with a twist, and
that freedom turns out to be harmless: a full ``2π`` sweep of the twist moves
the individual dipole strengths by ``1%`` and ``ΔA/A`` by less than ``1e-12``,
because the areal response sees the dipoles only through ``Σ_j D_j``, which
the solvability rows pin.  And the **source's channel data** are taken from
the local expansion of ``U`` about each mouth centre, exact through ``O(a²)``
in the ``ℓ = 0`` channel (the correction is free: ``∇²U = −3U − (2/π²) S_σ·x``
inside the ball, where the source does not reach) and leading-order in
``ℓ = 1``.  Its cost is priced by running the whole solve at two radii — and
bounded above by the fact that the answer barely depends on those data at all.

What the answer turns out to be made of
───────────────────────────────────────
Not what was expected, and worth stating plainly because it is what makes the
result robust.  ``ΔA/A`` is, to ``0.09%``, a **linear functional of the
obstruction ``S_σ`` alone** — of its in-plane part alone, exactly — and the
local source values ``U(c_j)`` and ``∇U(c_j)`` are a tenth-of-a-percent
correction.  Deleting the obstruction and keeping everything else leaves a
response a thousand times smaller.  So the signed answer rests on the single
best-converged number the source quadrature produces, which is not how these
rounds usually go.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Dict, Sequence, Tuple

import numpy as np

from .initial_data import KERNEL_PROJECTOR, regularised_green

__all__ = [
    "TubeModel",
    "WORKING_TUBE",
    "regularised_green_derivative",
    "SourceMoments",
    "INTERFERENCE_MOMENTS",
    "MOUTHS",
    "tangent_frame",
    "direction_rule",
    "mouth_sphere",
    "channel_split",
    "parallel_transport",
    "gluing_matrix",
    "basis_channels",
    "source_channels",
    "solve_matching",
    "measure_the_kernel_projector_is_the_green_functions_own_tail",
    "measure_the_matching_reproduces_an_exact_radial_solve",
    "measure_the_dipole_layers_are_required_not_optional",
    "measure_the_obstruction_carries_the_answer",
    "measure_the_signed_areal_response",
    "measure_the_sign_does_not_survive_a_matched_tube",
    "measure_the_throat_is_a_resonant_cavity",
]

FOUR_PI = 4.0 * math.pi
_Q = 4.0 * math.pi ** 2


# ════════════════════════════════════════════════════════════════════════════
# THE GREEN FUNCTION'S DERIVATIVE
# ════════════════════════════════════════════════════════════════════════════
def regularised_green_derivative(chi) -> np.ndarray:
    """``dG_⊥/dχ`` in closed form.

    Kept exact rather than differenced: the dipole layer *is* this function,
    and the dipole layers are what make the matching system solvable
    at all.
    """
    x = np.asarray(chi, dtype=float)
    s, c = np.sin(x), np.cos(x)
    return (-np.cos(2.0 * x) / s
            - 2.0 * (math.pi - x) * np.sin(2.0 * x) / s
            - (math.pi - x) * np.cos(2.0 * x) * c / s ** 2) / _Q


# ════════════════════════════════════════════════════════════════════════════
# THE TUBE
# ════════════════════════════════════════════════════════════════════════════
@dataclass(frozen=True)
class TubeModel:
    """A round tube of cross-section ``𝒜`` and length ``L``.

    The constraint operator on a product ``S²(r) × R`` with ``𝒜 = 4πr²`` is
    ``∇² + R̂/2 = ∂_s² + (1/r²)Δ_{S²} + 1/r²``, so the two channels the mouths
    can drive are ``∂_s² + k²`` and ``∂_s² − k²`` with the *same*
    ``k = 1/r = √(4π/𝒜)``.  One oscillates and one decays; the mass term shifts
    both.
    """

    area: float = FOUR_PI
    length: float = 0.9
    interior_mass: float = 0.0

    def wavenumber(self) -> float:
        return math.sqrt(FOUR_PI / float(self.area))

    def monopole_rate(self) -> complex:
        k2 = self.wavenumber() ** 2 - float(self.interior_mass) ** 2
        return math.sqrt(k2) if k2 >= 0.0 else 1j * math.sqrt(-k2)

    def dipole_rate(self) -> float:
        return math.sqrt(self.wavenumber() ** 2
                         + float(self.interior_mass) ** 2)

    def monopole_transfer(self) -> Tuple[float, float, float]:
        """``(k, cos kL, sin kL)`` — the oscillatory channel."""
        k = self.monopole_rate()
        if isinstance(k, complex):
            kk = float(k.imag)
            return kk, math.cosh(kk * self.length), -math.sinh(kk * self.length)
        return float(k), math.cos(k * self.length), math.sin(k * self.length)

    def dipole_transfer(self) -> Tuple[float, float, float]:
        """``(κ, cosh κL, sinh κL)`` — the evanescent channel.

        Kept for reference and for the one-dimensional checks.  `solve_matching`
        does **not** use it: see `dipole_attenuation`.
        """
        k = self.dipole_rate()
        return k, math.cosh(k * self.length), math.sinh(k * self.length)

    def dipole_attenuation(self) -> Tuple[float, float]:
        """``(κ, e^{−κL})`` — the evanescent channel written stably.

        A transfer matrix in ``cosh``/``sinh`` costs a condition number of
        ``e^{2κL}``, which is not a numerical detail: with the tube's area
        matched to the mouth's, ``κ = 1/sin a`` and a length of ``0.9`` gives
        ``e^{2κL} = 4.4e+15`` at ``a = 0.05`` — an eliminated system goes
        singular to double precision and no answer can be read out of it at all.

        Writing the tube mode as ``P e^{−κs} + Q e^{κ(s−L)}`` instead never
        forms ``e^{+κL}``: every coefficient is bounded by one, the two end
        amplitudes enter as their own unknowns, and the physical limit — a tube
        many e-foldings long, whose two mouths are decoupled in this channel —
        is the well-conditioned case rather than the singular one.
        """
        k = self.dipole_rate()
        return k, math.exp(-k * self.length)

    def monopole_resonances(self, count: int = 4) -> np.ndarray:
        """Lengths at which the ``ℓ = 0`` channel has a standing wave.

        These are poles of the areal response.  They are a property of the
        throat, not of the source, and the sign of ``ΔA/A`` flips across each.
        """
        k = self.monopole_rate()
        if isinstance(k, complex):
            return np.empty(0)
        return math.pi * np.arange(1, int(count) + 1) / float(k)


WORKING_TUBE = TubeModel(area=FOUR_PI, length=0.9, interior_mass=0.0)


# ════════════════════════════════════════════════════════════════════════════
# MOUTH SPHERES AND CHANNELS
# ════════════════════════════════════════════════════════════════════════════
def tangent_frame(centre: Sequence[float]) -> np.ndarray:
    """A deterministic orthonormal basis of ``c^⊥`` — Householder, not SVD.

    The same construction as `backreaction`, and for the same reason: an SVD
    null space is not canonical across LAPACK builds, and PR #263's CI failure
    was exactly that.
    """
    c = np.asarray(centre, dtype=float)
    c = c / np.linalg.norm(c)
    e = np.zeros(4)
    e[0] = 1.0
    v = c + (1.0 if c[0] >= 0.0 else -1.0) * e
    nv2 = float(v @ v)
    h = np.eye(4) if nv2 < 1e-24 else np.eye(4) - 2.0 * np.outer(v, v) / nv2
    return h[:, 1:].T


def direction_rule(n_theta: int = 16, n_phi: int = 32
                   ) -> Tuple[np.ndarray, np.ndarray]:
    """Gauss × uniform on the two-sphere of directions; weights sum to one.

    Exact for the ``ℓ ≤ 1`` moments this module extracts, by a wide margin —
    which matters, because the whole point is to separate two channels cleanly
    rather than to leak one into the other.
    """
    t, w = np.polynomial.legendre.leggauss(int(n_theta))
    ph = 2.0 * math.pi * np.arange(int(n_phi)) / int(n_phi)
    st = np.sqrt(np.maximum(0.0, 1.0 - t ** 2))
    dirs = np.stack([np.repeat(t, int(n_phi)),
                     np.outer(st, np.cos(ph)).ravel(),
                     np.outer(st, np.sin(ph)).ravel()], axis=1)
    return dirs, np.repeat(w, int(n_phi)) / int(n_phi) / 2.0


def mouth_sphere(centre: Sequence[float], radius: float,
                 dirs: np.ndarray) -> np.ndarray:
    """Points of ``S³`` at geodesic radius ``a`` from ``centre``."""
    c = np.asarray(centre, dtype=float)
    c = c / np.linalg.norm(c)
    return (math.cos(radius) * c[None, :]
            + math.sin(radius) * (np.asarray(dirs, float) @ tangent_frame(c)))


def channel_split(values: np.ndarray, dirs: np.ndarray,
                  weights: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """``(ℓ = 0, ℓ = 1)`` parts of a field sampled on a mouth sphere.

    Normalised so that a field ``α (m̂·n̂)`` returns ``ℓ = 1`` part ``α m̂``:
    the ``3`` is ``1/⟨n̂ n̂⟩``.
    """
    v = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    return (np.einsum("p,p...->...", w, v),
            3.0 * np.einsum("p,pi,p...->i...", w, np.asarray(dirs, float), v))


# ════════════════════════════════════════════════════════════════════════════
# THE GLUING
# ════════════════════════════════════════════════════════════════════════════
def parallel_transport(c1: Sequence[float], c2: Sequence[float]) -> np.ndarray:
    """The ``4×4`` rotation carrying ``c₁`` to ``c₂`` along their geodesic.

    It fixes the two directions orthogonal to both mouths and rotates the
    joining plane, which is what "carry the transverse frame through the
    throat without twisting it" means.
    """
    a = np.asarray(c1, float) / np.linalg.norm(c1)
    b = np.asarray(c2, float) / np.linalg.norm(c2)
    ct = float(np.clip(a @ b, -1.0, 1.0))
    d = math.acos(ct)
    if math.sin(d) < 1e-12:
        return np.eye(4)
    u = (b - ct * a) / math.sin(d)
    return (np.eye(4) + (math.cos(d) - 1.0) * (np.outer(a, a) + np.outer(u, u))
            + math.sin(d) * (np.outer(u, a) - np.outer(a, u)))


def gluing_matrix(c1: Sequence[float], c2: Sequence[float],
                  reflect: bool = False) -> np.ndarray:
    """``3×3`` map from mouth 2's transverse frame to mouth 1's.

    ``reflect`` flips one transverse axis — the discrete freedom a handle has
    that no amount of care about the exterior can remove.  It is a modelling
    choice and is reported as one.
    """
    t1, t2 = tangent_frame(c1), tangent_frame(c2)
    r = t1 @ parallel_transport(c1, c2).T @ t2.T
    if reflect:
        r = r @ np.diag([1.0, 1.0, -1.0])
    return r


# ════════════════════════════════════════════════════════════════════════════
# ASSEMBLY
# ════════════════════════════════════════════════════════════════════════════
def _field_values(kind: str, source: np.ndarray, direction: np.ndarray,
                  points: np.ndarray) -> np.ndarray:
    chi = np.arccos(np.clip(points @ source, -1.0, 1.0))
    if kind == "monopole":
        return regularised_green(chi)
    return -regularised_green_derivative(chi) * (points @ direction) / np.sin(chi)


def basis_channels(mouths: Sequence[Sequence[float]], radius: float,
                   n_theta: int = 16, n_phi: int = 32,
                   step: float = 1e-3) -> Dict[str, np.ndarray]:
    """Channel data of all twelve basis fields at both mouth spheres.

    Assembled by quadrature on the spheres rather than from a hand-written
    asymptotic expansion.  The basis fields are analytic and known in closed
    form, so this costs nothing and removes a whole class of algebra error —
    the one that produced ``ℓ(ℓ+2)``.
    """
    cs = [np.asarray(c, float) / np.linalg.norm(c) for c in mouths]
    frames = [tangent_frame(c) for c in cs]
    dirs, w = direction_rule(n_theta, n_phi)
    h = float(step) * float(radius)

    def columns(points: np.ndarray) -> np.ndarray:
        out = []
        for c, fr in zip(cs, frames):
            out.append(_field_values("monopole", c, c, points))
            for m in fr:
                out.append(_field_values("dipole", c, m, points))
        for a in range(4):
            out.append(points[:, a].copy())
        return np.array(out)                                  # 12 × P

    v0 = np.zeros((2, 12))
    d0 = np.zeros((2, 12))
    v1 = np.zeros((2, 3, 12))
    d1 = np.zeros((2, 3, 12))
    for j, c in enumerate(cs):
        here = columns(mouth_sphere(c, radius, dirs))
        # Fourth order on purpose.  The self-monopole behaves as 1/χ, for which
        # a two-point rule's *relative* truncation error is exactly ``step²`` —
        # 1e-06 at the obvious step, which is a systematic discrepancy against
        # the reference solve rather than the reference's own floor.  The
        # five-point rule puts it at step⁴ and out of the way.
        off = [columns(mouth_sphere(c, radius + s * h, dirs))
               for s in (-2, -1, 1, 2)]
        slope = (off[0] - 8.0 * off[1] + 8.0 * off[2] - off[3]) / (12.0 * h)
        v0[j], v1[j] = channel_split(here.T, dirs, w)
        d0[j], d1[j] = channel_split(slope.T, dirs, w)
    return {"value_0": v0, "slope_0": d0, "value_1": v1, "slope_1": d1}


def source_channels(mouths: Sequence[Sequence[float]], radius: float,
                    value: Sequence[float], gradient: Sequence[Sequence[float]],
                    obstruction: Sequence[float]) -> Dict[str, np.ndarray]:
    """Channel data of the source field ``U`` from its local expansion.

    ``value[j] = U(c_j)`` and ``gradient[j] = ∇U(c_j)`` as a four-vector
    orthogonal to ``c_j``.  Inside the ball the source does not reach, so
    ``∇²U = −3U − (2/π²) S_σ·x`` there exactly, and the ``ℓ = 0`` channel keeps
    its ``O(a²)`` term for free.  The ``ℓ = 1`` channel is leading order.
    """
    cs = [np.asarray(c, float) / np.linalg.norm(c) for c in mouths]
    s = np.asarray(obstruction, float)
    a = float(radius)
    v0 = np.zeros(2)
    d0 = np.zeros(2)
    v1 = np.zeros((2, 3))
    d1 = np.zeros((2, 3))
    for j, c in enumerate(cs):
        u = float(value[j])
        lap = -3.0 * u - KERNEL_PROJECTOR * float(s @ c)
        v0[j] = u + a ** 2 * lap / 6.0
        d0[j] = a * lap / 3.0
        g = tangent_frame(c) @ np.asarray(gradient[j], float)
        v1[j] = a * g
        d1[j] = g
    return {"value_0": v0, "slope_0": d0, "value_1": v1, "slope_1": d1}


def solve_matching(mouths: Sequence[Sequence[float]], radius: float,
                   tube: TubeModel, source: Dict[str, np.ndarray],
                   obstruction: Sequence[float], reflect: bool = False,
                   basis: Dict[str, np.ndarray] | None = None
                   ) -> Dict[str, object]:
    """Close the matching system and read off ``ΔA/A`` at each mouth.

    Rows 0–1 are the tube's oscillatory channel, rows 2–13 its three evanescent
    ones, rows 14–17 the kernel solvability condition.  Columns are
    ``[A₁, D₁, A₂, D₂, c]`` followed by the tube's two end amplitudes in each
    transverse direction.

    The twelve field unknowns are the physics; the six tube amplitudes are
    carried rather than eliminated purely for conditioning.  Eliminating them
    gives a twelve-by-twelve system with ``cosh κL`` and ``sinh κL`` in it, and
    a condition number of ``e^{2κL}`` — fine for a wide tube, fatal for a narrow
    one.  See `TubeModel.dipole_attenuation`.
    """
    cs = [np.asarray(c, float) / np.linalg.norm(c) for c in mouths]
    a = float(radius)
    area_mouth = FOUR_PI * math.sin(a) ** 2
    at = float(tube.area)
    b = basis_channels(cs, a) if basis is None else basis
    r = gluing_matrix(cs[0], cs[1], reflect=reflect)

    def rows(key: str) -> Tuple[np.ndarray, np.ndarray]:
        return b[key], source[key]

    v0, sv0 = rows("value_0")
    d0, sd0 = rows("slope_0")
    v1, sv1 = rows("value_1")
    d1, sd1 = rows("slope_1")

    m = np.zeros((18, 18))
    rhs = np.zeros(18)

    k, cc, ss = tube.monopole_transfer()
    scale = area_mouth / (at * k)
    m[0, :12] = v0[1] - cc * v0[0] + scale * ss * d0[0]
    rhs[0] = -(sv0[1] - cc * sv0[0] + scale * ss * sd0[0])
    m[1, :12] = at * k * ss * v0[0] + area_mouth * (cc * d0[0] + d0[1])
    rhs[1] = -(at * k * ss * sv0[0] + area_mouth * (cc * sd0[0] + sd0[1]))

    kd, x = tube.dipole_attenuation()
    scd = area_mouth / (at * kd)
    g1, sg1 = r @ v1[1], r @ sv1[1]
    h1, sh1 = r @ d1[1], r @ sd1[1]
    # w(s) = P e^{-k s} + Q e^{k(s-L)}: four conditions per transverse
    # direction, with P and Q carried as unknowns rather than eliminated.
    # Rows 2..13, columns 12..17.  The flux rows are divided through by
    # (area * k) so every entry is bounded by one.
    for i in range(3):
        m[2 + 4 * i, :12] = -v1[0][i]
        m[2 + 4 * i, 12 + i] = 1.0
        m[2 + 4 * i, 15 + i] = x
        rhs[2 + 4 * i] = sv1[0][i]

        m[3 + 4 * i, :12] = scd * d1[0][i]
        m[3 + 4 * i, 12 + i] = -1.0
        m[3 + 4 * i, 15 + i] = x
        rhs[3 + 4 * i] = -scd * sd1[0][i]

        m[4 + 4 * i, :12] = -g1[i]
        m[4 + 4 * i, 12 + i] = x
        m[4 + 4 * i, 15 + i] = 1.0
        rhs[4 + 4 * i] = sg1[i]

        m[5 + 4 * i, :12] = -scd * h1[i]
        m[5 + 4 * i, 12 + i] = -x
        m[5 + 4 * i, 15 + i] = 1.0
        rhs[5 + 4 * i] = scd * sh1[i]

    frames = [tangent_frame(c) for c in cs]
    for i, c in enumerate(cs):
        m[14:18, i * 4] = c
        m[14:18, i * 4 + 1: i * 4 + 4] = frames[i].T
    rhs[14:18] = np.asarray(obstruction, float)

    full = np.linalg.solve(m, rhs)
    coef = full[:12]
    u_mouth = v0 @ coef + sv0
    return {
        "coefficients": coef,
        "tube_amplitudes": full[12:].copy(),
        "monopoles": coef[[0, 4]].copy(),
        "dipoles": np.stack([frames[0].T @ coef[1:4],
                             frames[1].T @ coef[5:8]]),
        "kernel": coef[8:12].copy(),
        "conformal_factor": u_mouth,
        "dipole_at_mouth": v1 @ coef + sv1,
        "areal_response": 4.0 * u_mouth,
        "condition_number": float(np.linalg.cond(m)),
        "residual": float(np.linalg.norm(m @ full - rhs)
                          / max(np.linalg.norm(rhs), 1e-300)),
    }


# ════════════════════════════════════════════════════════════════════════════
# THE ONE-DIMENSIONAL REFERENCE SOLVES
# ════════════════════════════════════════════════════════════════════════════
def _radial_reference(radius: float, tube: TubeModel, ell: int,
                      profile: Callable[[float], float],
                      nodes: int = 900) -> object:
    """Exact solve on ``a < χ < π − a`` with the mouths at the two poles.

    A boundary-value problem in one variable, integrated to ``1e-09``, with the
    *same* throat conditions the matching system imposes.  It is the
    only check here that does not share code with the thing it checks.
    """
    from scipy.integrate import solve_bvp

    a = float(radius)
    am = FOUR_PI * math.sin(a) ** 2
    at = float(tube.area)
    if ell == 0:
        k, cc, ss = tube.monopole_transfer()
    else:
        k, cc, ss = tube.dipole_transfer()
    cent = 0.0 if ell == 0 else 2.0

    def rhs(x, y):
        return np.vstack([y[1],
                          np.array([profile(float(t)) for t in x])
                          - 2.0 * np.cos(x) / np.sin(x) * y[1]
                          - (3.0 - cent / np.sin(x) ** 2) * y[0]])

    def bc(ya, yb):
        v1, s1 = ya[0], ya[1]
        v2, s2 = yb[0], -yb[1]
        if ell == 0:
            return np.array([v2 - cc * v1 + (am / (at * k)) * ss * s1,
                             at * k * ss * v1 + am * (cc * s1 + s2)])
        return np.array([v2 - cc * v1 + (am / (at * k)) * ss * s1,
                         at * k * ss * v1 - am * (cc * s1 + s2)])

    lo = a * np.exp(np.linspace(0.0, math.log(0.5 * math.pi / a), int(nodes)))
    x = np.unique(np.concatenate([lo, [0.5 * math.pi], math.pi - lo[::-1]]))
    sol = solve_bvp(rhs, bc, x, np.zeros((2, x.size)), tol=1e-9,
                    max_nodes=900000)
    if sol.status != 0:
        raise RuntimeError(f"the reference solve did not converge: {sol.message}")
    return sol


def _homogeneous_pair(ell: int):
    """The two solutions of ``L f = 0`` in the radial sector, in closed form.

    ``ℓ = 0``: ``cos χ`` and ``−2 cos χ cot 2χ``.
    ``ℓ = 1``: ``sin χ`` and ``−cos χ − cos³χ/(3 sin²χ)``.

    In **both** sectors the regular member is a ``k = 1`` harmonic — that is the
    degeneracy, seen twice.  The three ``ℓ = 1`` partners are as degenerate as
    the ``ℓ = 0`` one and were invisible while the centrifugal term was wrong.
    Normalised so the Wronskian combination ``sin²χ (f₁f₂' − f₂f₁') = 1``.
    """
    if int(ell) == 0:
        return (lambda x: math.cos(x),
                lambda x: -2.0 * math.cos(x) / math.tan(2.0 * x))
    return (lambda x: math.sin(x),
            lambda x: -math.cos(x) - math.cos(x) ** 3 / (3.0 * math.sin(x) ** 2))


def _source_field(ell: int, profile: Callable[[float], float]):
    """``U`` and ``U'`` for a radial source, by variation of parameters.

    The projector term is subtracted first — without it there is no solution
    regular at both poles, which is the degeneracy stated as an obstruction
    rather than as a pole.
    """
    from scipy.integrate import quad

    f1, f2 = _homogeneous_pair(ell)
    if int(ell) == 0:
        weight = math.cos
        moment = FOUR_PI * quad(
            lambda t: math.cos(t) * profile(t) * math.sin(t) ** 2,
            0.0, math.pi, limit=400)[0]
    else:
        weight = math.sin
        moment = (FOUR_PI / 3.0) * quad(
            lambda t: profile(t) * math.sin(t) ** 3,
            0.0, math.pi, limit=400)[0]

    def g(t: float) -> float:
        return math.sin(t) ** 2 * (profile(t)
                                   - KERNEL_PROJECTOR * moment * weight(t))

    def field(x: float) -> float:
        i1 = quad(lambda t: f1(t) * g(t), 0.0, x, limit=400)[0]
        i2 = quad(lambda t: f2(t) * g(t), 0.0, x, limit=400, points=[0.0])[0]
        return f2(x) * i1 - f1(x) * i2

    def slope(x: float, h: float = 1e-6) -> float:
        return (field(x + h) - field(x - h)) / (2.0 * h)

    return field, slope, moment


def measure_the_kernel_projector_is_the_green_functions_own_tail(
) -> Dict[str, object]:
    """``L G_⊥ = −δ + (2/π²) cos χ``, and that tail **is** the projector.

    Two facts, checked separately rather than asserted together: the residual
    of ``G_⊥`` is ``(2/π²)cos χ`` away from the source, and ``(2/π²)cos χ`` is
    the integral kernel of the orthogonal projector onto ``span{x^A}``, because
    ``Σ_A x^A y^A = cos χ(x,y)`` and ``‖cos χ‖²_{S³} = π²/2``.

    This is what makes the solvability condition an identity rather than a
    modelling choice, so it is worth not taking on trust.
    """
    xs = np.array([0.3, 0.9, 1.6, 2.4, 3.0])
    h = 1e-5
    lap = ((regularised_green(xs + h) - 2.0 * regularised_green(xs)
            + regularised_green(xs - h)) / h ** 2
           + 2.0 / np.tan(xs) * regularised_green_derivative(xs))
    residual = lap + 3.0 * regularised_green(xs)
    tail = KERNEL_PROJECTOR * np.cos(xs)
    norm = FOUR_PI * float(np.trapezoid(
        np.cos(np.linspace(0, math.pi, 200001)) ** 2
        * np.sin(np.linspace(0, math.pi, 200001)) ** 2,
        np.linspace(0, math.pi, 200001)))
    near = regularised_green(1e-6) * (4.0 * math.pi * 1e-6)
    return {
        "residual_matches_the_tail": float(np.max(np.abs(residual - tail))),
        "kernel_norm_squared": norm,
        "closed_form_norm": math.pi ** 2 / 2.0,
        "projector_constant": KERNEL_PROJECTOR,
        "one_over_the_norm": 1.0 / norm,
        "unit_source_normalisation": float(near),
        "it_is_the_projector": bool(
            abs(1.0 / norm - KERNEL_PROJECTOR) < 1e-9
            and float(np.max(np.abs(residual - tail))) < 1e-5),
    }


def measure_the_matching_reproduces_an_exact_radial_solve(
        radii: Sequence[float] = (0.2, 0.1, 0.05),
        tube: TubeModel = WORKING_TUBE) -> Dict[str, object]:
    """The matching system against exact one-dimensional solves.

    Antipodal mouths, a radial source, and the same throat conditions — but the
    reference is a boundary-value solve that shares no code with the assembly.
    Run in **both** sectors on purpose.  Every closed-form check in PR #262 was
    ``ℓ = 0``, which is exactly why an ``ℓ = 1`` error survived that round's
    whole suite; a suite that cannot see a sector cannot defend it.
    """
    c1 = np.array([1.0, 0.0, 0.0, 0.0])
    c2 = -c1
    frames = [tangent_frame(c1), tangent_frame(c2)]
    out: Dict[str, object] = {"rows": []}
    worst = {0: 0.0, 1: 0.0}
    for ell, centre in ((0, 1.0), (1, 1.15)):
        def profile(t: float, m: float = centre) -> float:
            return math.exp(-((t - m) / 0.4) ** 2)
        field, slope, moment = _source_field(ell, profile)
        for a in radii:
            ref = _radial_reference(a, tube, ell, profile)
            src = {"value_0": np.zeros(2), "slope_0": np.zeros(2),
                   "value_1": np.zeros((2, 3)), "slope_1": np.zeros((2, 3))}
            if ell == 0:
                src["value_0"] = np.array([field(a), field(math.pi - a)])
                src["slope_0"] = np.array([slope(a), -slope(math.pi - a)])
                obstruction = moment * c1
            else:
                m_hat = frames[0][0]
                src["value_1"] = np.stack([
                    field(a) * (frames[0] @ m_hat),
                    field(math.pi - a) * (frames[1] @ m_hat)])
                src["slope_1"] = np.stack([
                    slope(a) * (frames[0] @ m_hat),
                    -slope(math.pi - a) * (frames[1] @ m_hat)])
                obstruction = moment * m_hat
            got = solve_matching([c1, c2], a, tube, src, obstruction)
            if ell == 0:
                mine = np.asarray(got["conformal_factor"], float)
            else:
                got1 = np.asarray(got["dipole_at_mouth"], float)
                mine = np.array([got1[j] @ (frames[j] @ m_hat)
                                 for j in (0, 1)])
            exact = np.array([float(ref.sol(a)[0]),
                              float(ref.sol(math.pi - a)[0])])
            rel = float(np.max(np.abs(mine / exact - 1.0)))
            worst[ell] = max(worst[ell], rel)
            out["rows"].append({"sector": f"l={ell}", "radius": float(a),
                                "exact": exact.tolist(),
                                "matched": mine.tolist(),
                                "relative_error": rel})
    errors = [float(r["relative_error"]) for r in out["rows"]]
    out.update({
        "worst_l0": worst[0],
        "worst_l1": worst[1],
        "worst_overall": float(max(errors)),
        "radii": [float(a) for a in radii],
        "both_sectors_agree": bool(max(errors) < 1e-8),
        "the_error_is_flat_in_the_radius": bool(
            max(errors) / min(errors) < 1e3),
        "what_flat_means": "the residual no longer falls with a, so it is the "
                           "assembly's numerical floor and not a truncation "
                           "error — an order of convergence would be a claim "
                           "this data cannot support",
        "what_it_defends": "PR #262's suite checked only l=0, so an l=1 error "
                           "in the centrifugal term survived it intact",
    })
    return out


# ════════════════════════════════════════════════════════════════════════════
# THE MEASURED SOURCE
# ════════════════════════════════════════════════════════════════════════════
@dataclass(frozen=True)
class SourceMoments:
    """Everything the matching needs to know about ``ΔT₀₀``, in eight numbers.

    Measured on ``Ω`` by the quadrature `backreaction` built and PR #263's
    review made deterministic, averaged over the ringdown window ``4 < t < 30``:

        ``S^A  = ∫_Ω y^A ΔT₀₀ dV``                       — the obstruction
        ``u0_j = ∫_Ω G_⊥(θ_j) ΔT₀₀ dV``                  — ``U(c_j)/2πG``
        ``u1_j = ∫_Ω G_⊥'(θ_j)(y − cos θ_j c_j)/sin θ_j ΔT₀₀ dV``

    so that ``U(c_j) = 2πG u0_j`` and ``∇U(c_j) = −2πG u1_j``.  Two structural
    checks come for free and both pass to ``1e-17``: ``u1_j ⟂ c_j``, and the
    ``ℓ = 1`` moment about a mouth on a coordinate axis has no component along
    that axis.
    """

    radius: float
    points: int
    volume: float
    energy: float
    obstruction: Tuple[float, ...]
    value: Tuple[float, ...]
    gradient: Tuple[Tuple[float, ...], ...]

    def as_source(self, coupling: float = 1.0) -> Dict[str, np.ndarray]:
        """Channel data at the mouth spheres, in units of ``2πG = coupling``."""
        return source_channels(
            MOUTHS, self.radius,
            [coupling * v for v in self.value],
            [[-coupling * g for g in row] for row in self.gradient],
            [-coupling * s for s in self.obstruction])

    def signed_obstruction(self, coupling: float = 1.0) -> np.ndarray:
        return -coupling * np.asarray(self.obstruction, float)


MOUTHS: Tuple[Tuple[float, ...], ...] = (
    (1.0, 0.0, 0.0, 0.0),
    (0.2674988286245874, 0.9635581854171929, 0.0, 0.0),
)


INTERFERENCE_MOMENTS: Tuple[SourceMoments, ...] = (
    SourceMoments(
        radius=0.05, points=5158, volume=19.682593700508534,
        energy=-0.007854399293613368,
        obstruction=(0.0034132270576741498, 0.0035183413980792045, -0.002979172945367291, -0.0024736064468879205),
        value=(0.003199581947273776, 0.0036432474738761914),
        gradient=((5.933946833581996e-21, -0.00037138341683181595, 0.00048325507527221227, 0.0002275428893235468), (-0.0006359297672851777, 0.000176544053499603, 0.00040381444883674344, 0.0009259667674323997))),
    SourceMoments(
        radius=0.05, points=12630, volume=19.735603393170642,
        energy=-0.010752488233688056,
        obstruction=(0.003365052052516453, 0.003428094097091135, -0.0029247941051175202, -0.0024988541762218553),
        value=(0.0034487468503370242, 0.0036324941611150952),
        gradient=((-4.306234895016468e-20, -0.0003858382056828556, 0.0006827019105937979, 0.00016790213903231614), (-0.0006970477325458289, 0.0001935113569407189, 0.0005170099779381257, 0.000896294884047609))),
    SourceMoments(
        radius=0.1, points=5158, volume=19.667965366268298,
        energy=-0.010593951070866605,
        obstruction=(0.0018862741824174484, 0.002101794586219982, -0.0029790563012842652, -0.0024738685323700666),
        value=(0.0019578105743441668, 0.001987360220029202),
        gradient=((3.986555891235616e-21, -0.0004160206729618793, 0.00045254912709998294, 0.00027318268952510836), (-0.0006934760047303373, 0.00019251978941396615, 0.00040711323785138986, 0.0009447128430025853))),
    SourceMoments(
        radius=0.1, points=12630, volume=19.720975058930406,
        energy=-0.013499300280663082,
        obstruction=(0.0018337143251692743, 0.00200776485850214, -0.0029246774609888566, -0.0024991162617931826),
        value=(0.00219361481516947, 0.001960764869072857),
        gradient=((-3.554957662813278e-20, -0.0004306465008529876, 0.0006519271377257116, 0.00021375791695511215), (-0.0007548228314964676, 0.00020955062838989104, 0.0005203586480752496, 0.00091521410119001))),
)


def measure_the_dipole_layers_are_required_not_optional(
        moments: SourceMoments | None = None) -> Dict[str, object]:
    """Required for the problem to *exist*, and invisible in its answer.

    The first half is exact linear algebra and needs no measurement:
    ``A₁c₁ + A₂c₂`` sweeps the plane ``span{c₁, c₂}`` and nothing else, so a
    two-parameter monopole model can meet at most two of the four solvability
    equations.  The dipole layers reach the rest, because ``D_j`` ranges over
    all of ``c_j^⊥`` and ``c_1^⊥ + c_2^⊥ = R⁴``.

    The second half is the measurement, and it went the other way from the
    expectation.  ``62.5%`` of the interference source's obstruction lies off
    the mouth plane, so without the layers there is no solution at all — but
    the response to that off-plane part is **zero**, not small: the layers
    deposit it in the kernel elements ``x²`` and ``x³``, which vanish at both
    mouths because both mouths lie in the ``(x⁰, x¹)`` plane.  The ``ℓ = 1``
    sector is load-bearing for existence and contributes nothing to the value.
    """
    cs = np.stack([np.asarray(c, float) for c in MOUTHS], axis=1)
    q, _ = np.linalg.qr(cs)
    m = moments or INTERFERENCE_MOMENTS[1]
    s = m.signed_obstruction()
    par = q @ (q.T @ s)
    perp = s - par
    fit, *_ = np.linalg.lstsq(cs, s, rcond=None)
    reach = np.linalg.matrix_rank(np.concatenate(
        [cs.T, tangent_frame(MOUTHS[0]), tangent_frame(MOUTHS[1])]))
    basis = basis_channels(MOUTHS, m.radius)
    quiet = {"value_0": np.zeros(2), "slope_0": np.zeros(2),
             "value_1": np.zeros((2, 3)), "slope_1": np.zeros((2, 3))}

    def response(obstruction: np.ndarray) -> np.ndarray:
        return np.asarray(solve_matching(MOUTHS, m.radius, WORKING_TUBE, quiet,
                                         obstruction, basis=basis
                                         )["areal_response"], float)

    full, in_plane, off_plane = response(s), response(par), response(perp)
    return {
        "monopole_span_dimension": int(np.linalg.matrix_rank(cs)),
        "equations": 4,
        "with_dipoles_dimension": int(reach),
        "obstruction": s.tolist(),
        "norm": float(np.linalg.norm(s)),
        "off_plane_fraction": float(np.linalg.norm(perp)
                                    / max(np.linalg.norm(s), 1e-300)),
        "monopole_only_shortfall": float(np.linalg.norm(cs @ fit - s)
                                         / max(np.linalg.norm(s), 1e-300)),
        "response_to_the_whole": full.tolist(),
        "response_to_the_in_plane_part": in_plane.tolist(),
        "response_to_the_off_plane_part": off_plane.tolist(),
        "monopoles_alone_cannot_close_it": bool(
            int(np.linalg.matrix_rank(cs)) < 4
            and np.linalg.norm(cs @ fit - s) > 0.5 * np.linalg.norm(s)),
        "the_dipoles_close_it": bool(reach == 4),
        "and_then_they_do_not_move_the_answer": bool(
            np.max(np.abs(off_plane)) < 1e-12 * np.max(np.abs(full))),
    }


def measure_the_obstruction_carries_the_answer(
        moments: SourceMoments | None = None) -> Dict[str, object]:
    """``ΔA/A`` is a linear functional of ``S_σ``, and of nothing else that matters.

    Three separations, run rather than argued:

    * drop the local source data ``U(c_j)``, ``∇U(c_j)`` and keep only the
      obstruction — the answer survives to ``0.09%``;
    * drop the obstruction and keep the local data — the answer collapses by
      three orders of magnitude;
    * scale the ``ℓ = 1`` moments by ``3``, by ``0``, or replace them with
      noise — the answer moves by ``5e-04`` at worst.

    That last one is why this is here and not a footnote.  The ``ℓ = 1``
    moments are the worst-converged input the source quadrature produces —
    ``41%`` drift between levels, against ``1.5%`` for the obstruction —
    and the signed result does not rest on them.
    """
    m = moments or INTERFERENCE_MOMENTS[1]
    basis = basis_channels(MOUTHS, m.radius)
    s = m.signed_obstruction()
    quiet = {"value_0": np.zeros(2), "slope_0": np.zeros(2),
             "value_1": np.zeros((2, 3)), "slope_1": np.zeros((2, 3))}

    def run(src, obstruction) -> np.ndarray:
        return np.asarray(solve_matching(MOUTHS, m.radius, WORKING_TUBE, src,
                                         obstruction, basis=basis
                                         )["areal_response"], float)

    full = run(m.as_source(), s)
    rng = np.random.default_rng(20260820)
    noisy = SourceMoments(
        radius=m.radius, points=m.points, volume=m.volume, energy=m.energy,
        obstruction=m.obstruction, value=m.value,
        gradient=tuple(tuple(g * (1.0 + 0.5 * rng.standard_normal())
                             for g in row) for row in m.gradient))
    flat = SourceMoments(
        radius=m.radius, points=m.points, volume=m.volume, energy=m.energy,
        obstruction=m.obstruction, value=m.value,
        gradient=tuple(tuple(0.0 for _ in row) for row in m.gradient))
    variants = {
        "obstruction_only": run(quiet, s),
        "local_data_only": run(m.as_source(), np.zeros(4)),
        "dipole_moments_zeroed": run(flat.as_source(), s),
        "dipole_moments_randomised": run(noisy.as_source(), s),
        "obstruction_doubled": run(quiet, 2.0 * s),
    }
    scale = max(float(np.max(np.abs(full))), 1e-300)
    return {
        "full": full.tolist(),
        "variants": {k: v.tolist() for k, v in variants.items()},
        "obstruction_alone_reproduces_it": float(
            np.max(np.abs(variants["obstruction_only"] / full - 1.0))),
        "without_the_obstruction": float(
            np.max(np.abs(variants["local_data_only"])) / scale),
        "worst_drift_from_the_dipole_moments": float(max(
            np.max(np.abs(variants[k] / full - 1.0))
            for k in ("dipole_moments_zeroed", "dipole_moments_randomised"))),
        "linear_in_the_obstruction": float(np.max(np.abs(
            variants["obstruction_doubled"] / variants["obstruction_only"]
            - 2.0))),
        "the_obstruction_is_the_answer": bool(
            np.max(np.abs(variants["obstruction_only"] / full - 1.0)) < 5e-3
            and np.max(np.abs(variants["local_data_only"])) < 1e-2 * scale),
    }


def measure_the_signed_areal_response(
        moments: Sequence[SourceMoments] | None = None,
        tube: TubeModel = WORKING_TUBE,
        coupling: float = 1.0) -> Dict[str, object]:
    """**The answer** — at the wide working throat, off resonance.

    ``ΔA/A`` at each mouth, with its sign, reported in units of ``2πG``, in
    which the whole problem is linear: the response scales with ``G`` and with
    the square of the wave amplitude, and neither can change a sign.

    The qualifier in the first line is load-bearing.  This throat has
    ``𝒜 = 4π``, four hundred times its own mouth area, which puts it at
    ``kL = 0.9`` — inside the first cavity cell.  Matching the tube's area to
    the mouths' instead flips the sign, and
    `measure_the_sign_does_not_survive_a_matched_tube` is that calculation.
    What follows is a statement about a throat, not about the source.

    Four controls travel with it, because on this problem a number with no
    control attached has twice turned out to be noise:

    * the quadrature level of the source moments — two levels, a factor
      ``2.4`` in points;
    * the mouth radius ``a``, doubled, which prices the local expansion of
      ``U`` and is *not* a regulator: the mouth radius is a parameter of the
      throat, the source goes as ``1/χ⁴`` there, and there is no ``a → 0``
      limit to converge to.  That is the singular point PR #262 removed;
    * the gluing of the two transverse frames, transported or reflected;
    * the conditioning of the matching system and its residual.

    The sign is the same in all eight combinations.  The magnitude at fixed
    ``a`` is stable to ``2.2%`` across quadrature levels, and moves by a
    factor ``1.75`` when the mouth radius is doubled — almost exactly the
    factor by which the source's own mouth-weighted moment moves, because a
    bigger ball excludes more of the ``1/χ⁴`` pile-up.
    """
    ms = list(moments or INTERFERENCE_MOMENTS)
    rows = []
    for m in ms:
        basis = basis_channels(MOUTHS, m.radius)
        for reflect in (False, True):
            got = solve_matching(MOUTHS, m.radius, tube, m.as_source(coupling),
                                 m.signed_obstruction(coupling),
                                 reflect=reflect, basis=basis)
            rows.append({
                "radius": m.radius, "points": m.points,
                "gluing": "reflected" if reflect else "transported",
                "areal_response": np.asarray(got["areal_response"]).tolist(),
                "monopoles": np.asarray(got["monopoles"]).tolist(),
                "dipole_norms": [float(np.linalg.norm(d))
                                 for d in np.asarray(got["dipoles"])],
                "kernel_norm": float(np.linalg.norm(got["kernel"])),
                "condition_number": got["condition_number"],
                "residual": got["residual"],
            })
    vals = np.array([r["areal_response"] for r in rows])
    finest = min(m.radius for m in ms)
    best = max(m.points for m in ms if m.radius == finest)
    head = next(r for r in rows if r["radius"] == finest
                and r["points"] == best and r["gluing"] == "transported")
    ref = np.array(head["areal_response"])
    same_a = np.array([r["areal_response"] for r in rows
                       if r["radius"] == finest])
    return {
        "rows": rows,
        "headline": head,
        "areal_response": ref.tolist(),
        "sign": ["closes" if v < 0 else "opens" for v in ref],
        "every_variant_agrees_in_sign": bool(
            np.all(np.sign(vals) == np.sign(vals[0]))
            and np.all(np.sign(vals[0]) != 0)),
        "quadrature_spread_at_fixed_radius": float(
            np.max(np.abs(same_a / ref - 1.0))),
        "spread_over_the_radius": float(np.max(np.abs(vals / ref - 1.0))),
        "worst_condition_number": float(max(r["condition_number"]
                                            for r in rows)),
        "worst_residual": float(max(r["residual"] for r in rows)),
        "the_answer": "toward a neck AT THIS THROAT — the conformal factor "
                      "falls at both mouths, so both mouth areas contract.  The "
                      "interference energy alone would open them (U(c_j) > 0 at "
                      "both); the throat's monopole layers overshoot that and "
                      "invert it.  Matching the tube's area to the mouths' "
                      "flips the sign, so this is a property of the throat and "
                      "not of the source.",
        "the_qualifier": "wide tube, area 4 pi against a mouth area of "
                         "4 pi sin^2 a; kL = 0.9, inside the first cavity cell",
    }


def measure_the_sign_does_not_survive_a_matched_tube(
        moments: Sequence[SourceMoments] | None = None,
        length: float | None = None) -> Dict[str, object]:
    """**The headline is a statement about a wide tube, and here is the proof.**

    `WORKING_TUBE` carries area ``4π`` while a mouth sphere has area
    ``4π sin²a`` — a factor of ``400`` at ``a = 0.05``.  That is a deliberate
    idealisation, a wide throat entered through small mouths, and it is the one
    `neck` has used since PR #261.  The obvious alternative is to set the two
    equal, which makes the tube exactly as narrow as its own mouths.

    Then ``k = √(4π/𝒜) = 1/sin a``, and the *same* length ``0.9`` becomes
    ``kL/π = 5.73`` at ``a = 0.05`` and ``2.87`` at ``a = 0.10``.  The throat is
    no longer inside the first cavity cell; it is past five poles and past two,
    and it sits ``4.6%`` of its own length from the next one.

    **The sign does not survive.**  At ``a = 0.05`` both mouths *open*; at
    ``a = 0.10`` they disagree — mouth 1 closes and mouth 2 opens.  Neither is a
    numerical artefact: the conditioning is reported alongside, and it is what
    made this check necessary in the first place (see
    `TubeModel.dipole_attenuation`).

    So `measure_the_signed_areal_response` should be read as: *at the wide
    working throat, off resonance, both mouths close.*  Not as a property of the
    interference source.
    """
    ms = [m for m in (moments or INTERFERENCE_MOMENTS)
          if m.points == max(x.points for x in (moments or INTERFERENCE_MOMENTS)
                             if x.radius == m.radius)]
    seen: Dict[float, SourceMoments] = {}
    for m in ms:
        seen.setdefault(m.radius, m)
    lng = float(WORKING_TUBE.length if length is None else length)
    rows = []
    for radius, m in sorted(seen.items()):
        basis = basis_channels(MOUTHS, radius)
        for name, area in (("wide", WORKING_TUBE.area),
                           ("matched", FOUR_PI * math.sin(radius) ** 2)):
            tube = TubeModel(area=area, length=lng)
            k = tube.wavenumber()
            got = solve_matching(MOUTHS, radius, tube, m.as_source(),
                                 m.signed_obstruction(), basis=basis)
            response = np.asarray(got["areal_response"], float)
            phase = k * lng / math.pi
            # n = 0 is not a pole: a tube of zero length is not a cavity.
            near = min((max(1, math.floor(phase)), max(1, math.ceil(phase))),
                       key=lambda n: abs(phase - n))
            rows.append({
                "radius": radius, "model": name, "area": float(area),
                "wavenumber": float(k), "phase_over_pi": float(phase),
                "areal_response": response.tolist(),
                "sign": ["closes" if v < 0 else "opens" for v in response],
                "nearest_pole": int(near),
                "distance_to_the_pole_over_pi": float(abs(phase - near)),
                "distance_as_a_fraction_of_the_length": float(
                    abs(near * math.pi / k - lng) / lng),
                "condition_number": got["condition_number"],
                "residual": got["residual"],
            })
    wide = [r for r in rows if r["model"] == "wide"]
    matched = [r for r in rows if r["model"] == "matched"]
    return {
        "rows": rows,
        "length": lng,
        "wide_signs": [r["sign"] for r in wide],
        "matched_signs": [r["sign"] for r in matched],
        "the_wide_throat_always_closes": bool(
            all(v == "closes" for r in wide for v in r["sign"])),
        "the_matched_throat_does_not": bool(
            any(v == "opens" for r in matched for v in r["sign"])),
        "the_matched_mouths_can_disagree": bool(
            any(len(set(r["sign"])) == 2 for r in matched)),
        "worst_condition_number": float(max(r["condition_number"]
                                            for r in rows)),
        "worst_residual": float(max(r["residual"] for r in rows)),
        "the_sign_is_a_property_of_the_throat": bool(
            all(v == "closes" for r in wide for v in r["sign"])
            and any(v == "opens" for r in matched for v in r["sign"])),
    }


def measure_the_throat_is_a_resonant_cavity(
        lengths: Sequence[float] | None = None,
        moments: SourceMoments | None = None,
        area: float = FOUR_PI) -> Dict[str, object]:
    """The sign is a property of the throat, and it flips at ``kL = nπ``.

    The ``ℓ = 0`` tube channel is ``∂_s² + k²`` with ``k = √(4π/𝒜)``: a cavity.
    At a standing-wave length the exterior sees an infinite impedance, the
    response has a pole, and ``ΔA/A`` changes sign across it.

    This is the honest limit on how far the headline generalises.  It is not a
    numerical instability — the poles sit at closed-form lengths and the scan
    lands on them to the resolution of its own grid.  The working throat is at
    ``kL = 0.9``, well inside the first cell.
    """
    m = moments or INTERFERENCE_MOMENTS[0]
    k = math.sqrt(FOUR_PI / float(area))
    ls = np.asarray(lengths if lengths is not None
                    else np.linspace(0.2, 3.0 * math.pi / k, 240), float)
    vals = []
    for length in ls:
        t = TubeModel(area=area, length=float(length))
        got = solve_matching(MOUTHS, m.radius, t, m.as_source(),
                             m.signed_obstruction())
        vals.append(np.asarray(got["areal_response"], float))
    vals = np.array(vals)
    s = np.sign(vals[:, 0])
    flips = ls[np.where(np.diff(s) != 0)[0]]
    poles = math.pi * np.arange(1, 4) / k
    poles = poles[poles <= float(ls[-1])]
    nearest = [float(np.min(np.abs(flips - p))) if flips.size else float("nan")
               for p in poles]
    return {
        "wavenumber": k,
        "closed_form_poles": poles.tolist(),
        "sign_flips_at": flips.tolist(),
        "distance_to_the_closed_form": nearest,
        "grid_spacing": float(ls[1] - ls[0]),
        "flips_land_on_the_poles": bool(
            all(d <= 2.0 * float(ls[1] - ls[0]) for d in nearest[:2])),
        "working_length": float(WORKING_TUBE.length),
        "working_phase": float(k * WORKING_TUBE.length),
        "the_working_throat_is_off_resonance": bool(
            k * WORKING_TUBE.length < math.pi),
    }
