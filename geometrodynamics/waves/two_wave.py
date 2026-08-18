"""The dynamical two-wave invariant, against its known WKB limit.

Roadmap step 3, properly this time.  PR #258 built a *static* source-interaction
kernel and was explicit that it was **not** this: no local null momenta, so it
could not distinguish equal-energy collinear from counterpropagating waves — the
one control the two-wave invariant exists to apply.  This module builds the
time-dependent object and applies that control.

**The point is not to re-derive the WKB identity.**  ``𝒞 = A_A²A_B²(k_A·k_B)²``
with ``k_A·k_B = −E_AE_B(1 − cos θ)`` is known: zero collinear, ``−2E_AE_B``
head-on.  The research content is the **difference** between the exact
time-dependent, multipath, throat-coupled field and that limit — how big it is,
what it is made of, and how fast it closes.

What is solved
──────────────
The retarded field of a pulsed point source on the ESU with a self-adjoint
point throat, exactly, by Krein's resolvent formula in the frequency domain:

    φ_s(x,ω) = s(ω)[ G(χ_{xs},ω) + Σ_ij G(χ_{xi},ω) R_ij(ω) G(χ_{js},ω) ]

with ``R(ω) = (A − Γ(ω²))⁻¹``, inverted back along the **retarded contour**
``Im ω = ε``.  That contour is exact, not approximate: writing ``ω = u + iε``,

    φ(t)  =  e^{εt} · (1/2π) ∫ du  e^{−iut} φ̂(u + iε)

so ``ε`` only trades pole sharpness against the growth factor, and both are
reported.  The free part is checked against PR #254's closed-form image sum,
which never saw a frequency integral.

Derivatives without a mesh
──────────────────────────
Every term is a function of one geodesic distance, so the full four-gradient and
Hessian follow in closed form from radial derivatives and the sphere's geometry:

    ∇_a χ = −(y − (x·y)x)/sin χ            (unit, pointing away from the source)
    ∇_a∇_b χ = cot χ (δ_ab − ∇_aχ ∇_bχ)

Both verified numerically against the exponential map.  Nothing is
finite-differenced, so the stress tensor is exact to the solver's own accuracy —
which the tracelessness of the conformal ``T_{μν}`` then measures.

The stress tensor
─────────────────
The improved (conformally coupled, ``ξ = 1/6``) tensor,

    T_{μν} = ∂_μφ ∂_νφ − ½ η_{μν}(∂φ)² + ξ[η_{μν}□ − ∇_μ∇_ν + G_{μν}]φ²

in an orthonormal frame, with ``G_{00} = 3``, ``G_{ab} = −δ_ab`` for ``S³ × R``.
``□φ`` is taken from the solved derivatives rather than replaced by its on-shell
value: substituting on shell would make ``T^μ_μ`` vanish *algebraically* for any
input, so the trace would test nothing.  Computed honestly it equals
``φ(□φ − φ)``, the wave-equation residual — which is what makes the measured
trace a test of the **solver**.

What is put in
──────────────
A linear field on a fixed background, the mouth positions, and the boundary
data — four real numbers chosen, not derived, and taken strictly inside PR
#257's stability cone with the Löwner margin quoted.  No backreaction: the
stress tensor is computed *from* the field and never fed back.
"""

from __future__ import annotations

import cmath
import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .field_solve import image_field
from .throat_operator import MouthPair
from .throat_positivity import positivity_defect
from .two_source import disconnection_defect, geodesic, mouth_positions

__all__ = [
    "TWO_PI",
    "green_omega",
    "green_omega_derivatives",
    "gamma_omega",
    "GaussianPulse",
    "RetardedGrid",
    "orthonormal_frame",
    "radial_frame_data",
    "TwoWaveSetup",
    "solve_field",
    "stress_tensor",
    "contract_stress",
    "wkb_stress",
    "wkb_invariant",
    "collinear_and_head_on",
    "measure_the_solver_reproduces_the_closed_form_free_field",
    "measure_the_solved_field_satisfies_the_conformal_wave_equation",
    "measure_the_improved_stress_tensor_is_traceless",
    "measure_the_wkb_collinear_head_on_result_is_recovered",
    "measure_multipath_destroys_the_collinear_null",
    "measure_the_arrivals_are_the_branch_ledger_with_maslov_signs",
    "measure_the_only_tail_is_the_throats",
    "measure_the_caustic_is_where_wkb_stops",
    "measure_the_low_frequency_limit_recovers_the_tomography",
    "superpose",
    "zero_like",
    "cross_stress_tensor",
    "two_leg_channels",
    "measure_the_cross_mouth_channels_are_labelled_by_the_exit_mouth",
    "measure_the_interference_tensor_is_largest_where_the_invariant_is_null",
]

TWO_PI = 2.0 * math.pi
XI = 1.0 / 6.0                      # conformal coupling in D = 4
RICCI_SCALAR = 6.0                  # of the unit-radius ESU


# ════════════════════════════════════════════════════════════════════════════
# THE FREQUENCY-DOMAIN KERNEL, FOR COMPLEX ω
# ════════════════════════════════════════════════════════════════════════════
def green_omega(chi: float, omega):
    """``G(χ, ω) = sin(ω(π−χ)) / (4π sin χ · sin πω)``, continued to complex ω.

    Array-aware in ``ω``.  Written in ``e = π − χ`` for the reason PR #257's
    review gave: the antipode is a **removable** singularity, and forming
    ``π − χ`` from a float near ``π`` loses its digits.
    """
    w = np.asarray(omega, dtype=complex)
    e = math.pi - float(chi)
    se = math.sin(e)
    sp = np.sin(math.pi * w)
    scale = 4.0 * math.pi * sp
    if abs(se) < 1e-9:                   # sin(ωe)/sin(e) → ω(1 + e²(1−ω²)/6)
        out = w * (1.0 + e * e * (1.0 - w * w) / 6.0) / scale
    else:
        out = np.sin(w * e) / (se * scale)
    return out if out.ndim else complex(out)


def green_omega_derivatives(chi: float, omega):
    """``(G, ∂_χ G, ∂²_χ G)`` — analytic, not finite-differenced, array-aware.

    With ``e = π − χ``, ``N = sin(ωe)``, ``D = sin e`` and ``H = N/D``:

        ``H' = N'/D − N D'/D²``
        ``H'' = N''/D − 2N'D'/D² − N D''/D² + 2N D'²/D³``

    and ``∂_χ = −∂_e``, so ``∂_χG = −H'/(4π sin πω)`` and
    ``∂²_χG = +H''/(4π sin πω)``.
    """
    w = np.asarray(omega, dtype=complex)
    e = math.pi - float(chi)
    scale = 4.0 * math.pi * np.sin(math.pi * w)
    se = math.sin(e)
    if abs(se) < 1e-9:                       # H = ω[1 + e²(1−ω²)/6 + …]
        w2 = w * w
        h = w * (1.0 + e * e * (1.0 - w2) / 6.0)
        hp = w * e * (1.0 - w2) / 3.0
        hpp = w * (1.0 - w2) / 3.0
    else:
        ce = math.cos(e)
        n0 = np.sin(w * e)
        n1 = w * np.cos(w * e)
        n2 = -w * w * n0
        d0, d1, d2 = se, ce, -se
        h = n0 / d0
        hp = n1 / d0 - n0 * d1 / d0 ** 2
        hpp = (n2 / d0 - 2.0 * n1 * d1 / d0 ** 2 - n0 * d2 / d0 ** 2
               + 2.0 * n0 * d1 ** 2 / d0 ** 3)
    return h / scale, -hp / scale, hpp / scale


def gamma_omega(omega, separation: float) -> np.ndarray:
    """``Γ(ω)`` for complex ``ω`` — ``g`` on the diagonal, ``G_d`` off it.

    ``g(ω) = −(ω/4π) cot(πω)`` is the regularized coincidence limit; the
    off-diagonal is `green_omega` at the mouth separation.  Agreement with
    `throat_operator.gamma_at` on the real axis is measured rather than assumed.
    """
    w = np.asarray(omega, dtype=complex)
    g = -w * np.cos(math.pi * w) / (4.0 * math.pi * np.sin(math.pi * w))
    gd = green_omega(float(separation), w)
    if w.ndim == 0:
        return np.array([[g, gd], [gd, g]], dtype=complex)
    out = np.empty(w.shape + (2, 2), dtype=complex)
    out[..., 0, 0] = out[..., 1, 1] = g
    out[..., 0, 1] = out[..., 1, 0] = gd
    return out


# ════════════════════════════════════════════════════════════════════════════
# SOURCES AND THE RETARDED CONTOUR
# ════════════════════════════════════════════════════════════════════════════
@dataclass(frozen=True)
class GaussianPulse:
    """A band-limited real source, ``a·cos(ω₀(t−t₀))·exp(−(t−t₀)²/2w²)``.

    Its spectrum is a pair of Gaussians at ``±ω₀``, so ``ω₀`` is a genuine
    carrier and ``1/w`` is the bandwidth — the two knobs the WKB limit is taken
    with.  ``ω₀ = 0`` degenerates to the plain Gaussian PR #254 used, which is
    what makes the free-field regression possible.

    **Normalization is peak, not area**: ``amplitude`` is the value at ``t₀``.
    PR #254's `image_field` uses a *unit-area* Gaussian, so the regression
    against it runs at ``amplitude = 1/(w√(2π))`` — a convention, stated here
    because a silent factor of ``6.6`` is exactly the kind of thing that gets
    mistaken for a solver error.
    """

    amplitude: float = 1.0
    carrier: float = 12.0
    width: float = 0.18
    t0: float = 0.0

    def spectrum(self, omega):
        """Array-aware.  ``∫dt e^{iωt}`` of the pulse."""
        z = np.asarray(omega, dtype=complex)
        w = float(self.width)
        a = float(self.amplitude) * w * math.sqrt(TWO_PI)
        phase = np.exp(1j * z * float(self.t0))
        if self.carrier == 0.0:
            out = a * np.exp(-(z * w) ** 2 / 2.0) * phase
        else:
            c = float(self.carrier)
            out = 0.5 * a * (np.exp(-((z - c) * w) ** 2 / 2.0)
                             + np.exp(-((z + c) * w) ** 2 / 2.0)) * phase
        return out if out.ndim else complex(out)


@dataclass(frozen=True)
class RetardedGrid:
    """The contour ``Im ω = ε`` and its time grid.

    ``φ(t) = e^{εt}·(1/2π)∫du e^{−iut} φ̂(u+iε)``.  The ``e^{εt}`` is exact, so
    ``ε`` is a pure numerical trade: larger ``ε`` smooths the poles and
    amplifies the inverse transform's error by ``e^{εT}``, which is reported.
    """

    n: int = 1 << 16
    span: float = 400.0
    eps: float = 0.05

    @property
    def dt(self) -> float:
        return self.span / self.n

    @property
    def times(self) -> np.ndarray:
        return np.arange(self.n) * self.dt

    @property
    def omegas(self) -> np.ndarray:
        return np.fft.fftfreq(self.n, d=self.dt) * TWO_PI + 1j * self.eps

    def invert(self, spectrum: np.ndarray) -> np.ndarray:
        """Back to the time domain, along the contour.

        On the grid ``(1/2π)∫du e^{−iut}F(u)`` is
        ``Σ_k F_k e^{−2πikn/N}·du/2π``, and that kernel is ``np.fft.fft``'s,
        not ``ifft``'s — ``ifft`` evaluates the transform at ``−t`` and returns
        ~0 for a retarded field, which is exactly how the sign was caught.
        """
        out = np.fft.fft(spectrum) / self.span
        return np.real(out) * np.exp(self.eps * self.times)

    def growth(self, t_max: float) -> float:
        return float(math.exp(self.eps * t_max))


# ════════════════════════════════════════════════════════════════════════════
# FRAMES AND RADIAL DATA
# ════════════════════════════════════════════════════════════════════════════
def orthonormal_frame(x: Sequence[float], seed: int = 7) -> np.ndarray:
    """Three orthonormal tangent vectors at ``x`` on ``S³``, as rows."""
    p = np.asarray(x, dtype=float)
    p = p / np.linalg.norm(p)
    rng = np.random.default_rng(int(seed))
    rows: List[np.ndarray] = []
    while len(rows) < 3:
        v = rng.normal(size=4)
        v -= np.dot(v, p) * p
        for b in rows:
            v -= np.dot(v, b) * b
        nrm = np.linalg.norm(v)
        if nrm > 1e-6:
            rows.append(v / nrm)
    return np.array(rows)


def radial_frame_data(x: Sequence[float], y: Sequence[float],
                      frame: np.ndarray) -> Tuple[float, np.ndarray]:
    """``(χ, n̂)`` — the geodesic distance and the unit tangent **away** from
    ``y``, in frame components.

    ``∇χ = −(y − (x·y)x)/sin χ``, verified against the exponential map.
    """
    p = np.asarray(x, dtype=float)
    q = np.asarray(y, dtype=float)
    chi = geodesic(p, q)
    s = math.sin(chi)
    if s < 1e-9:
        raise ValueError("the observation point is at a source or its antipode")
    grad = -(q - float(np.dot(p, q)) * p) / s
    return chi, np.array([float(np.dot(grad, b)) for b in frame])


# ════════════════════════════════════════════════════════════════════════════
# THE SOLVE
# ════════════════════════════════════════════════════════════════════════════
@dataclass
class TwoWaveSetup:
    """Two pulsed point sources on a throated ESU, and where they are watched.

    The throat is a `MouthPair`, taken strictly inside PR #257's cone; its
    Löwner margin is reported by `margin` and quoted wherever a number is.
    """

    pair: MouthPair
    source_a: np.ndarray
    source_b: np.ndarray
    observer: np.ndarray
    pulse_a: GaussianPulse
    pulse_b: GaussianPulse
    grid: RetardedGrid = field(default_factory=RetardedGrid)
    with_throat: bool = True

    def margin(self) -> float:
        return float(positivity_defect(self.pair.boundary_matrix(),
                                       self.pair.separation)["min_eigenvalue"])

    def mouths(self) -> Tuple[np.ndarray, np.ndarray]:
        return mouth_positions(self.pair.separation)


def _terms(setup: TwoWaveSetup, source: np.ndarray, pulse: GaussianPulse,
           ) -> List[Tuple[float, np.ndarray, np.ndarray]]:
    """The radial terms of one source's field: ``(χ, n̂, kernel(ω))``.

    Direct propagation plus one term per mouth, the mouth kernels carrying the
    throat's response ``q = R(ω) G(c, y) s(ω)``.  With ``with_throat`` false the
    mouth terms are dropped, which is the null model every comparison needs.
    """
    frame = orthonormal_frame(setup.observer)
    om = setup.grid.omegas
    spec = pulse.spectrum(om)
    chi_d, n_d = radial_frame_data(setup.observer, source, frame)
    out = [(chi_d, n_d, spec)]
    if not setup.with_throat:
        return out
    c1, c2 = setup.mouths()
    d_src = np.array([geodesic(c1, source), geodesic(c2, source)])
    a = setup.pair.boundary_matrix()
    # M(ω) = A − Γ(ω), inverted in closed form: 2×2 beats a per-frequency solve
    m = a[None, :, :] - gamma_omega(om, setup.pair.separation)
    det = m[:, 0, 0] * m[:, 1, 1] - m[:, 0, 1] * m[:, 1, 0]
    drive = np.stack([green_omega(float(d_src[0]), om),
                      green_omega(float(d_src[1]), om)], axis=-1)
    qs = np.empty_like(drive)
    qs[:, 0] = (m[:, 1, 1] * drive[:, 0] - m[:, 0, 1] * drive[:, 1]) / det
    qs[:, 1] = (-m[:, 1, 0] * drive[:, 0] + m[:, 0, 0] * drive[:, 1]) / det
    qs *= spec[:, None]
    for i, c in enumerate((c1, c2)):
        chi_i, n_i = radial_frame_data(setup.observer, c, frame)
        out.append((chi_i, n_i, qs[:, i]))
    return out


def solve_field(setup: TwoWaveSetup, source: np.ndarray,
                pulse: GaussianPulse) -> Dict[str, np.ndarray]:
    """The field of one source at the observation point, with derivatives.

    Returns time series for ``φ``, ``∂_tφ``, ``∂_t²φ``, the spatial gradient
    ``∂_aφ``, the mixed ``∂_t∂_aφ``, and the spatial Hessian ``∇_a∇_bφ``, all in
    the orthonormal frame.  Every one comes from the same contour integral with
    an analytic ``χ``- or ``ω``-factor in the integrand — nothing is
    differenced.
    """
    grid = setup.grid
    om = grid.omegas
    n = grid.n
    phi = np.zeros(n)
    dt1 = np.zeros(n)
    dt2 = np.zeros(n)
    grad = np.zeros((n, 3))
    dtgrad = np.zeros((n, 3))
    hess = np.zeros((n, 3, 3))
    for chi, nhat, kernel in _terms(setup, source, pulse):
        g0, g1, g2 = green_omega_derivatives(chi, om)
        f0 = grid.invert(kernel * g0)
        f1 = grid.invert(kernel * g1)
        f2 = grid.invert(kernel * g2)
        ft0 = grid.invert(-1j * om * kernel * g0)
        ft1 = grid.invert(-1j * om * kernel * g1)
        ftt = grid.invert(-(om ** 2) * kernel * g0)
        cot = 1.0 / math.tan(chi)
        proj = np.eye(3) - np.outer(nhat, nhat)
        phi += f0
        dt1 += ft0
        dt2 += ftt
        grad += np.outer(f1, nhat)
        dtgrad += np.outer(ft1, nhat)
        hess += (f2[:, None, None] * np.outer(nhat, nhat)[None, :, :]
                 + (f1 * cot)[:, None, None] * proj[None, :, :])
    return {"phi": phi, "dt": dt1, "dtt": dt2, "grad": grad,
            "dtgrad": dtgrad, "hess": hess,
            "laplacian": np.einsum("taa->t", hess)}


# ════════════════════════════════════════════════════════════════════════════
# THE STRESS TENSOR
# ════════════════════════════════════════════════════════════════════════════
_ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
_EINSTEIN = np.diag([3.0, -1.0, -1.0, -1.0])


def stress_tensor(sol: Dict[str, np.ndarray], index: int) -> np.ndarray:
    """``T_{μν}`` at one time index, in the orthonormal frame.

        ``T_{μν} = ∂_μφ∂_νφ − ½η_{μν}(∂φ)² + ξ[η_{μν}□ − ∇_μ∇_ν + G_{μν}]φ²``

    with ``□φ² = 2[(∂φ)² + φ□φ]``.

    ``□φ`` is taken from the **solved derivatives**, ``−∂_t²φ + ∇²φ``, and not
    replaced by its on-shell value ``ξRφ = φ``.  That matters: substituting on
    shell makes ``T^μ_μ`` vanish *algebraically* for any input, so the trace
    would test nothing.  Computed this way,

        ``T^μ_μ = φ(□φ − φ)``

    which is the wave-equation residual — a real measurement of the solver.
    """
    phi = float(sol["phi"][index])
    d = np.empty(4)
    d[0] = sol["dt"][index]
    d[1:] = sol["grad"][index]
    dd = np.empty((4, 4))
    dd[0, 0] = sol["dtt"][index]
    dd[0, 1:] = sol["dtgrad"][index]
    dd[1:, 0] = sol["dtgrad"][index]
    dd[1:, 1:] = sol["hess"][index]
    grad_sq = float(d @ _ETA @ d)            # η^{μν}∂_μφ∂_νφ, with lower d
    box_phi = float(sol["laplacian"][index] - sol["dtt"][index])
    box_phi2 = 2.0 * (grad_sq + phi * box_phi)
    hess_phi2 = 2.0 * (np.outer(d, d) + phi * dd)
    return (np.outer(d, d) - 0.5 * _ETA * grad_sq
            + XI * (_ETA * box_phi2 - hess_phi2 + _EINSTEIN * phi * phi))


def contract_stress(t_a: np.ndarray, t_b: np.ndarray) -> float:
    """``T_A^{μν} T^B_{μν}`` — indices raised with ``η``."""
    raised = _ETA @ np.asarray(t_a) @ _ETA
    return float(np.einsum("ij,ij->", raised, np.asarray(t_b)))


def trace_of(t: np.ndarray) -> float:
    """``T^μ_μ`` — identically zero for the conformal tensor on shell."""
    return float(np.einsum("ij,ji->", _ETA, np.asarray(t)))


# ════════════════════════════════════════════════════════════════════════════
# THE WKB LIMIT
# ════════════════════════════════════════════════════════════════════════════
def wkb_stress(amplitude: float, k: Sequence[float]) -> np.ndarray:
    """``T^{WKB}_{μν} = A² k_μ k_ν`` — the leading geometric-optics tensor."""
    kk = np.asarray(k, dtype=float)
    return amplitude * amplitude * np.outer(kk, kk)


def wkb_invariant(amp_a: float, k_a: Sequence[float],
                  amp_b: float, k_b: Sequence[float]) -> float:
    """``𝒞 = A_A²A_B²(k_A·k_B)²`` — the roadmap's invariant, in its own limit.

    With both ``k`` null and future-directed, ``k_A·k_B = −E_AE_B(1 − cos θ)``:
    **zero for collinear**, ``−2E_AE_B`` head-on.  That contrast is the control,
    and it is what the exact field is measured against rather than re-derived.
    """
    dot = float(np.asarray(k_a) @ _ETA @ np.asarray(k_b))
    return (amp_a * amp_b * dot) ** 2


def collinear_and_head_on(separation_ab: float = 0.9,
                          reach: float = 0.8) -> Dict[str, np.ndarray]:
    """Two source points and two observers on one great circle.

    Both configurations put the sources and the observer on a common geodesic,
    so the arriving spatial directions are exactly parallel or exactly
    antiparallel — the WKB control is set up by geometry, not by fitting.

    * **collinear** — the observer is *beyond* ``B`` on the ray from ``A``, so
      both signals travel the same way and ``k_A·k_B = 0``;
    * **head-on** — the observer sits *between* them, so the signals meet and
      ``k_A·k_B = −2E_AE_B``.
    """
    def on_circle(angle: float) -> np.ndarray:
        return np.array([math.cos(angle), math.sin(angle), 0.0, 0.0])

    a = on_circle(0.0)
    b = on_circle(separation_ab)
    return {"source_a": a, "source_b": b,
            "collinear": on_circle(separation_ab + reach),
            "head_on": on_circle(separation_ab - reach)}


# ════════════════════════════════════════════════════════════════════════════
# THE EXPERIMENT
# ════════════════════════════════════════════════════════════════════════════
WORKING_SEPARATION = 1.3
WORKING_BOUNDARY = (0.30, 0.35, 0.06 + 0.0j)
SOURCE_GAP = 1.6
OBSERVER_REACH = 0.8


def working_pair() -> MouthPair:
    a1, a2, b = WORKING_BOUNDARY
    return MouthPair(WORKING_SEPARATION, a1, a2, b)


def source_circle() -> Tuple[np.ndarray, np.ndarray]:
    """An orthonormal pair spanning the sources' great circle.

    Deliberately **not** the mouths' plane: `two_source.mouth_positions` puts
    them in the ``(x₀,x₁)`` plane, and a source placed there would sit on a
    mouth.  This basis is tilted out of it, so every geodesic distance in the
    experiment is generic and nothing is accidentally degenerate.
    """
    u1 = np.array([0.3, 0.1, 0.9, 0.0])
    u1 /= np.linalg.norm(u1)
    v = np.array([0.0, 0.2, 0.1, 1.0])
    v -= float(np.dot(v, u1)) * u1
    return u1, v / np.linalg.norm(v)


def circle_point(theta: float) -> np.ndarray:
    u1, u2 = source_circle()
    return math.cos(theta) * u1 + math.sin(theta) * u2


def normalized_invariant(setup: TwoWaveSetup, t_star: float,
                         half_window: float) -> Dict[str, float]:
    """``𝒩 = (T_A:T_B)/(T_A^{00} T_B^{00})`` at the peak of the energy product.

    **Pointwise, not window-averaged.**  For two WKB waves ``T_s = a_s²k_μk_ν``
    every envelope and every ``sin²`` factor cancels between numerator and
    denominator, so the pointwise ratio is exactly ``(1 − cos θ)²`` — no
    averaging, no phase-correlation bookkeeping.  A window average does *not*
    have that property: it drags ``⟨sin⁴⟩/⟨sin²⟩²`` and the envelope overlap
    into the answer, which is how an earlier draft got ``8.3`` where the answer
    is ``4``.
    """
    sol_a = solve_field(setup, setup.source_a, setup.pulse_a)
    sol_b = solve_field(setup, setup.source_b, setup.pulse_b)
    ts = setup.grid.times
    idx = np.where((ts > t_star - half_window) & (ts < t_star + half_window))[0]
    best: Optional[Tuple[float, float, float]] = None
    for i in idx:
        ta = stress_tensor(sol_a, int(i))
        tb = stress_tensor(sol_b, int(i))
        u = float(ta[0, 0] * tb[0, 0])
        if best is None or u > best[0]:
            best = (u, contract_stress(ta, tb) / u, float(ts[i]))
    if best is None:
        raise ValueError("the window contains no samples")
    return {"energy_product": best[0], "invariant": best[1], "t": best[2]}


def arrival_directions(observer: Sequence[float], point_a: Sequence[float],
                       point_b: Sequence[float]) -> float:
    """``n̂_A · n̂_B`` at the observer, from geometry alone.

    The WKB invariant is ``(1 − n̂_A·n̂_B)²``, so this is the *prediction* every
    measured ``𝒩`` below is compared against — computed from positions, never
    fitted to the field.
    """
    frame = orthonormal_frame(observer)
    _, na = radial_frame_data(observer, point_a, frame)
    _, nb = radial_frame_data(observer, point_b, frame)
    return float(np.dot(na, nb))


def _setup(observer: np.ndarray, t0_a: float, t0_b: float, carrier: float,
           width: float, throat: bool, grid: RetardedGrid) -> TwoWaveSetup:
    return TwoWaveSetup(
        pair=working_pair(), source_a=circle_point(0.0),
        source_b=circle_point(SOURCE_GAP), observer=observer,
        pulse_a=GaussianPulse(carrier=carrier, width=width, t0=t0_a),
        pulse_b=GaussianPulse(carrier=carrier, width=width, t0=t0_b),
        grid=grid, with_throat=throat)


def _fine_grid() -> RetardedGrid:
    """``ε`` comfortably above the frequency spacing ``2π/span``.

    That inequality is the whole convergence condition: with ``ε ≲ du`` the
    poles fall between grid points and the answer is wrong by orders of
    magnitude — measured, in `measure_the_wkb_collinear_head_on_result_is
    _recovered`.
    """
    return RetardedGrid(n=1 << 17, span=600.0, eps=0.05)


# ════════════════════════════════════════════════════════════════════════════
# MEASUREMENTS
# ════════════════════════════════════════════════════════════════════════════
def measure_the_solver_reproduces_the_closed_form_free_field(
        chi: float = 1.7, width: float = 0.06) -> Dict[str, object]:
    """The frequency-domain solve against PR #254's image sum.

    Two constructions that share no code: one integrates ``G(χ,ω)`` along the
    retarded contour, the other sums Gaussians over the winding images with
    alternating signs.  Agreement at the arrival times — including the **sign
    flips**, which are the Maslov index — is what licenses everything after it.

    The normalization is stated rather than absorbed: `GaussianPulse` is
    peak-normalized and `image_field` area-normalized, so the comparison runs at
    ``amplitude = 1/(w√(2π))``.
    """
    amp = 1.0 / (width * math.sqrt(TWO_PI))
    pulse = GaussianPulse(amplitude=amp, carrier=0.0, width=width)
    obs = circle_point(chi)
    setup = TwoWaveSetup(pair=working_pair(), source_a=circle_point(0.0),
                         source_b=circle_point(0.0), observer=obs,
                         pulse_a=pulse, pulse_b=pulse,
                         grid=RetardedGrid(n=1 << 15, span=300.0, eps=0.05),
                         with_throat=False)
    sol = solve_field(setup, setup.source_a, pulse)
    ts, dt = setup.grid.times, setup.grid.dt
    rows = []
    for k, t_arr in enumerate((chi, TWO_PI - chi, TWO_PI + chi,
                               2.0 * TWO_PI - chi)):
        i = int(round(t_arr / dt))
        got = float(sol["phi"][i])
        want = image_field(chi, float(ts[i]), width=width, n_images=10)
        rows.append({"t": float(ts[i]), "solver": got, "image_sum": want,
                     "difference": got - want,
                     "maslov_sign": int(math.copysign(1.0, want))})
    worst = max(abs(r["difference"]) for r in rows)
    signs = [r["maslov_sign"] for r in rows]
    return {"chi": chi, "rows": rows,
            "worst_difference": float(worst),
            "peak_scale": float(max(abs(r["image_sum"]) for r in rows)),
            "the_two_constructions_agree": bool(worst < 1e-12),
            "the_signs_alternate": bool(signs == [1, -1, 1, -1]),
            "what_this_licenses": ("the frequency-domain solve, on which every "
                                   "stress tensor below is built")}


def measure_the_solved_field_satisfies_the_conformal_wave_equation(
        ) -> Dict[str, object]:
    """``∂_t²φ = ∇²φ − φ`` from the solved derivatives, throat included.

    ``ξR = 1`` on the unit ESU, so the conformal equation is
    ``(−∂_t² + ∇² − 1)φ = 0``.  Nothing here is finite-differenced: the time
    derivatives are ``−iω`` inside the contour integral and the spatial ones are
    the analytic ``∂_χG``, ``∂²_χG`` plus the sphere's Hessian
    ``cot χ (δ_ab − n̂_an̂_b)``.  The residual is therefore a test of that whole
    construction at once.
    """
    grid = _fine_grid()
    obs = circle_point(SOURCE_GAP + OBSERVER_REACH)
    out = {}
    for throat in (False, True):
        setup = _setup(obs, 0.0, 0.0, 16.0, 0.10, throat, grid)
        sol = solve_field(setup, setup.source_a, setup.pulse_a)
        ts = grid.times
        sl = (ts > 0.5) & (ts < 9.0)
        resid = sol["dtt"][sl] - (sol["laplacian"][sl] - sol["phi"][sl])
        scale = float(np.abs(sol["dtt"][sl]).max())
        out["with_throat" if throat else "free"] = {
            "worst_residual": float(np.abs(resid).max()),
            "scale": scale,
            "relative": float(np.abs(resid).max() / scale)}
    return {**out,
            "loewner_margin": _setup(obs, 0, 0, 16.0, 0.1, True,
                                     grid).margin(),
            "the_equation_holds": bool(
                out["free"]["relative"] < 1e-12
                and out["with_throat"]["relative"] < 1e-12),
            "nothing_is_differenced": True}


def measure_the_improved_stress_tensor_is_traceless() -> Dict[str, object]:
    """``T^μ_μ = 0`` — and it is a solver test, not an algebraic one.

    With ``□φ`` replaced by its on-shell value the trace would vanish
    *identically* for any input, measuring nothing.  Taken from the solved
    derivatives it equals ``φ(□φ − φ)``, so a nonzero trace is a nonzero
    wave-equation residual.  Reported relative to the largest component, at
    both a free and a throat-coupled event.
    """
    grid = _fine_grid()
    obs = circle_point(SOURCE_GAP + OBSERVER_REACH)
    rows = []
    for throat in (False, True):
        setup = _setup(obs, 0.0, 0.0, 16.0, 0.10, throat, grid)
        sol = solve_field(setup, setup.source_a, setup.pulse_a)
        chi = geodesic(obs, setup.source_a)
        i = int(round(chi / grid.dt))
        t = stress_tensor(sol, i)
        rows.append({"throat": throat, "T00": float(t[0, 0]),
                     "trace": trace_of(t),
                     "relative": float(abs(trace_of(t))
                                       / float(np.abs(t).max()))})
    return {"rows": rows,
            "worst_relative_trace": float(max(r["relative"] for r in rows)),
            "the_tensor_is_traceless": bool(
                max(r["relative"] for r in rows) < 1e-12),
            "why_it_is_not_vacuous": ("□φ is taken from the solve, not "
                                      "substituted on shell; the trace then "
                                      "equals φ(□φ − φ)")}


def measure_the_wkb_collinear_head_on_result_is_recovered(
        carriers: Sequence[float] = (6.0, 8.0, 12.0, 16.0, 24.0, 32.0, 48.0),
        width: float = 0.10) -> Dict[str, object]:
    """The known limit, as a limit — with the rate, not just the value.

    Two sources and an observer on one great circle, so the arriving directions
    are exactly parallel or exactly antiparallel and the control is set by
    geometry rather than fitted.  WKB says ``𝒩 = (1 − n̂_A·n̂_B)²``: **4**
    head-on, **0** collinear.  What is measured is the approach.

    The collinear null turns out to be far stronger than a leading-order
    statement.  On this geometry the two arriving wavefronts share their normal
    *exactly* — amplitude gradients are along the same ``n̂`` — so the residue
    falls off much faster than the ``ω^{-2}`` a generic curvature correction
    would give, and the exponent steepens with ``ω``.  That is what makes the
    multipath result next door a genuinely large effect rather than a
    competition between small ones.

    Convergence is part of the measurement: the contour needs ``ε`` well above
    the frequency spacing ``2π/span``, and the same run at ``ε ≈ du`` is wrong
    by orders of magnitude.
    """
    grid = _fine_grid()
    obs_col = circle_point(SOURCE_GAP + OBSERVER_REACH)
    obs_head = circle_point(SOURCE_GAP - OBSERVER_REACH)
    rows = []
    prev: Optional[Tuple[float, float]] = None
    for c in carriers:
        out = {"carrier": float(c)}
        for name, obs in (("collinear", obs_col), ("head_on", obs_head)):
            chi_a = geodesic(obs, circle_point(0.0))
            chi_b = geodesic(obs, circle_point(SOURCE_GAP))
            t_star = 3.0
            setup = _setup(obs, t_star - chi_a, t_star - chi_b, float(c),
                           width, False, grid)
            got = normalized_invariant(setup, t_star, 2.0 * width)
            dot = arrival_directions(obs, setup.source_a, setup.source_b)
            out[name] = got["invariant"]
            out[name + "_wkb"] = (1.0 - dot) ** 2
            out[name + "_dot"] = dot
        if prev is not None:
            out["collinear_slope"] = float(
                math.log(abs(out["collinear"] / prev[1]))
                / math.log(float(c) / prev[0]))
        prev = (float(c), out["collinear"])
        rows.append(out)
    coarse = RetardedGrid(n=1 << 16, span=300.0, eps=0.02)
    chi_a = geodesic(obs_col, circle_point(0.0))
    chi_b = geodesic(obs_col, circle_point(SOURCE_GAP))
    bad = normalized_invariant(
        _setup(obs_col, 3.0 - chi_a, 3.0 - chi_b, 32.0, width, False, coarse),
        3.0, 2.0 * width)["invariant"]
    good = [r for r in rows if r["carrier"] == 32.0]
    return {"rows": rows, "width": width,
            "head_on_at_the_largest_carrier": rows[-1]["head_on"],
            "head_on_target": 4.0,
            "head_on_error": abs(rows[-1]["head_on"] - 4.0),
            "collinear_at_the_largest_carrier": rows[-1]["collinear"],
            "the_directions_are_exactly_parallel": bool(
                all(abs(r["collinear_dot"] - 1.0) < 1e-12 for r in rows)),
            "the_directions_are_exactly_antiparallel": bool(
                all(abs(r["head_on_dot"] + 1.0) < 1e-12 for r in rows)),
            "head_on_converges_to_four": bool(
                abs(rows[-1]["head_on"] - 4.0) < 1e-3),
            "collinear_converges_to_zero": bool(
                abs(rows[-1]["collinear"]) < 1e-8),
            "the_collinear_exponent_steepens": bool(
                rows[-1].get("collinear_slope", 0.0)
                < rows[2].get("collinear_slope", 0.0)),
            "under_resolved_contour_value": float(bad),
            "converged_value_there": float(good[0]["collinear"]) if good
            else float("nan"),
            "the_contour_needs_eps_above_the_spacing": bool(
                abs(bad) > 100.0 * abs(good[0]["collinear"])
                if good else False)}


def measure_multipath_destroys_the_collinear_null() -> Dict[str, object]:
    """**The result.**  Same two sources, same event, three different answers.

    The WKB invariant is a function of the arriving *directions*, and a
    multipath field arrives on several.  Holding the sources and the observation
    point fixed and moving only *which branch has arrived*:

    * both signals on their **direct** branches — the configuration is collinear
      and ``𝒩`` is at its numerical floor;
    * ``A`` on its **long-way winding image** — that branch propagates the other
      way round the sphere, so its arrival direction is *reversed* and the same
      pair reads **head-on**, ``𝒩 ≈ 4``;
    * ``B`` through the **cross-mouth** channel — it emerges from a mouth, at an
      angle set by that mouth's position, and ``𝒩`` takes the intermediate value
      geometry predicts.

    So the collinear null is not destroyed by curvature corrections, which are
    ``1e-7`` here; it is destroyed by **multipath**, at ``O(1)``.  That is the
    branch-resolved invariant PR #255 said was needed, and it is exactly what
    PR #258's static kernel could not see.

    Two scoping notes, both measured elsewhere in this module rather than
    asserted here.  The mouth row is compared against a *prediction from the
    mouth positions*, not a fit; the shortest of the four two-leg paths is used,
    and `measure_the_cross_mouth_channels_are_labelled_by_the_exit_mouth` audits
    all four explicitly.  And the free-propagation control at the same instant
    is reported alongside, because there ``B`` has **no arrival at all**, so the
    ratio would be a meaningless ``0/0`` — the comparison is stated as
    amplitudes.  That control says the *mouths* create the branch; it does
    **not** say their *connection* does, and the ``β = 0`` control in the audit
    shows it does not.
    """
    grid = _fine_grid()
    carrier, width, t_star = 24.0, 0.10, 3.0
    obs = circle_point(SOURCE_GAP + OBSERVER_REACH)
    src_a, src_b = circle_point(0.0), circle_point(SOURCE_GAP)
    chi_a, chi_b = geodesic(obs, src_a), geodesic(obs, src_b)
    c1, c2 = mouth_positions(WORKING_SEPARATION)

    def row(name: str, t0a: float, t0b: float, throat: bool,
            predicted: float) -> Dict[str, object]:
        setup = _setup(obs, t0a, t0b, carrier, width, throat, grid)
        got = normalized_invariant(setup, t_star, 2.0 * width)
        return {"branch_pair": name, "invariant": got["invariant"],
                "geometric_prediction": predicted,
                "energy_product": got["energy_product"]}

    direct = row("A direct + B direct", t_star - chi_a, t_star - chi_b, False,
                 (1.0 - arrival_directions(obs, src_a, src_b)) ** 2)
    long_way = row("A long-way image + B direct",
                   t_star - (TWO_PI - chi_a), t_star - chi_b, False,
                   (1.0 + arrival_directions(obs, src_a, src_b)) ** 2)

    legs = [(geodesic(src_b, cj) + geodesic(ci, obs), i)
            for i, ci in enumerate((c1, c2)) for cj in (c1, c2)]
    delay, mouth_index = min(legs)
    mouth = (c1, c2)[mouth_index]
    through = row("A direct + B via a mouth", t_star - chi_a,
                  t_star - delay, True,
                  (1.0 - arrival_directions(obs, src_a, mouth)) ** 2)
    control = _setup(obs, t_star - chi_a, t_star - delay, carrier, width,
                     False, grid)
    control_out = normalized_invariant(control, t_star, 2.0 * width)

    rows = [direct, long_way, through]
    err = abs(through["invariant"] - through["geometric_prediction"])
    return {"carrier": carrier, "observer_reach": OBSERVER_REACH,
            "loewner_margin": _setup(obs, 0, 0, carrier, width, True,
                                     grid).margin(),
            "rows": rows,
            "throat_delay": float(delay),
            "throat_exit_mouth": int(mouth_index) + 1,
            "collinear_floor": direct["invariant"],
            "long_way_value": long_way["invariant"],
            "through_the_throat_value": through["invariant"],
            "through_the_throat_prediction": through["geometric_prediction"],
            "throat_relative_error": float(
                err / through["geometric_prediction"]),
            "free_control_energy_product": control_out["energy_product"],
            "throat_energy_product": through["energy_product"],
            "the_control_has_no_second_arrival": bool(
                control_out["energy_product"]
                < 1e-6 * through["energy_product"]),
            "the_direct_pair_is_null": bool(abs(direct["invariant"]) < 1e-6),
            "the_winding_image_reads_head_on": bool(
                abs(long_way["invariant"] - 4.0) < 5e-3),
            "the_throat_matches_its_geometry": bool(
                err / through["geometric_prediction"] < 1e-2),
            "the_lesson": ("the two-wave invariant is branch-resolved: the "
                           "same sources at the same event give 0, 4 or the "
                           "mouth's angle depending on which branch arrived"),
            "what_the_control_does_not_say": ("that the mouths' CONNECTION "
                                              "supplies the branch — β = 0 "
                                              "gives the same invariant; see "
                                              "the cross-mouth audit")}


def measure_the_arrivals_are_the_branch_ledger_with_maslov_signs(
        ) -> Dict[str, object]:
    """Where the solved field has support, and with which sign.

    The free arrivals are PR #253's branch set with PR #254's Maslov signs, read
    off a solve that never saw either: sharp peaks at ``χ``, ``2π − χ``,
    ``2π + χ`` with signs ``+ − +``.

    The **throat** adds arrivals the free ledger does not contain, at the
    two-leg times ``χ(y,c_j) + χ(c_i,x)`` PR #255 named.  Those are checked at
    the **causal onset** rather than at a peak, and deliberately so: the mouth
    response ``R(ω)`` has poles, so a throat arrival is not a sharp pulse but a
    ring-up — the same fact `measure_the_only_tail_is_the_throats` reports as a
    tail.  Asking a broadened, overlapping series of arrivals to have peaks at
    the geometric times is asking the wrong question, and the earliest onset is
    the part that is genuinely sharp.
    """
    grid = _fine_grid()
    obs = circle_point(SOURCE_GAP + OBSERVER_REACH)
    src = circle_point(0.0)
    width = 0.06
    pulse = GaussianPulse(amplitude=1.0 / (width * math.sqrt(TWO_PI)),
                          carrier=0.0, width=width)
    c1, c2 = mouth_positions(WORKING_SEPARATION)
    chi = geodesic(obs, src)
    free_times = [chi, TWO_PI - chi, TWO_PI + chi]
    two_leg = sorted(geodesic(src, cj) + geodesic(ci, obs)
                     for ci in (c1, c2) for cj in (c1, c2))
    ts = grid.times

    free_setup = TwoWaveSetup(pair=working_pair(), source_a=src, source_b=src,
                              observer=obs, pulse_a=pulse, pulse_b=pulse,
                              grid=grid, with_throat=False)
    free_sol = solve_field(free_setup, src, pulse)
    rows = []
    for t_arr in free_times:
        i = int(round(t_arr / grid.dt))
        win = slice(max(0, i - 40), i + 40)
        j = int(np.argmax(np.abs(free_sol["phi"][win]))) + win.start
        rows.append({"branch": "free", "predicted_t": float(t_arr),
                     "found_t": float(ts[j]), "offset": float(ts[j] - t_arr),
                     "value": float(free_sol["phi"][j])})

    throat_setup = TwoWaveSetup(pair=working_pair(), source_a=src,
                                source_b=src, observer=obs, pulse_a=pulse,
                                pulse_b=pulse, grid=grid, with_throat=True)
    throat_sol = solve_field(throat_setup, src, pulse)
    only = throat_sol["phi"] - free_sol["phi"]        # the throat's own field
    scale = float(np.abs(only[ts < 8.0]).max())
    early = np.where((ts > 0.2) & (np.abs(only) > 0.02 * scale))[0]
    onset = float(ts[early[0]]) if len(early) else float("nan")
    predicted_onset = float(two_leg[0])
    free_signs = [int(math.copysign(1, r["value"])) for r in rows]
    return {"free_arrivals": free_times,
            "two_leg_times": [float(t) for t in two_leg],
            "rows": rows,
            "worst_free_offset": float(max(abs(r["offset"]) for r in rows)),
            "free_signs": free_signs,
            "the_free_signs_alternate": bool(free_signs == [1, -1, 1]),
            "the_free_arrivals_are_sharp": bool(
                max(abs(r["offset"]) for r in rows) < 0.01),
            "throat_onset_measured": onset,
            "throat_onset_predicted": predicted_onset,
            "onset_offset": float(onset - predicted_onset),
            "the_throat_onset_is_causal": bool(
                -3.0 * width < onset - predicted_onset < 3.0 * width),
            "the_throat_arrivals_are_new": bool(
                all(min(abs(t - f) for f in free_times) > 0.1
                    for t in two_leg)),
            "why_onsets_and_not_peaks": ("R(ω) has poles, so a throat arrival "
                                         "rings up rather than pulsing; the "
                                         "onset is the sharp part")}


def measure_the_only_tail_is_the_throats() -> Dict[str, object]:
    """Huygens holds exactly for the free field; the throat breaks it.

    The ESU is conformally flat, so a conformally coupled massless scalar
    propagates strictly **on** the light cones: between geometric arrivals the
    free field is zero, and the measured floor is the pulse's own Gaussian wing
    plus round-off.  The throat's response ``R(ω)`` is not a pure delay — it has
    poles — so the mouths ring, and the field between arrivals is no longer
    zero.

    So the "tail correction" this round was asked to quantify has a sharp
    answer: for the free field there is none, and every bit of it belongs to the
    throat.
    """
    grid = _fine_grid()
    obs = circle_point(SOURCE_GAP + OBSERVER_REACH)
    src = circle_point(0.0)
    pulse = GaussianPulse(carrier=16.0, width=0.10)
    chi = geodesic(obs, src)
    ts = grid.times
    mask = (ts > 0.2) & (ts < 5.0)
    for t_arr in (chi, TWO_PI - chi, TWO_PI + chi):
        mask &= np.abs(ts - t_arr) > 0.6
    out = {}
    for throat in (False, True):
        setup = TwoWaveSetup(pair=working_pair(), source_a=src, source_b=src,
                             observer=obs, pulse_a=pulse, pulse_b=pulse,
                             grid=grid, with_throat=throat)
        sol = solve_field(setup, src, pulse)
        peak = float(np.abs(sol["phi"][ts < 6.0]).max())
        between = float(np.abs(sol["phi"][mask]).max())
        out["with_throat" if throat else "free"] = {
            "peak": peak, "between_arrivals": between,
            "ratio": between / peak}
    return {**out,
            "free_ratio": out["free"]["ratio"],
            "throat_ratio": out["with_throat"]["ratio"],
            "amplification": float(out["with_throat"]["ratio"]
                                   / out["free"]["ratio"]),
            "the_free_field_has_no_tail": bool(out["free"]["ratio"] < 1e-6),
            "the_throat_has_one": bool(out["with_throat"]["ratio"] > 1e-2),
            "why": ("S³ × R is conformally flat, so the conformal scalar obeys "
                    "Huygens exactly; R(ω) has poles, so the mouths ring")}


def measure_the_caustic_is_where_wkb_stops(
        carriers: Sequence[float] = (8.5, 16.5, 32.5)) -> Dict[str, object]:
    """The antipodal focus, and the scale at which geometric optics fails.

    Geometric optics gives the amplitude ``1/(4π sin χ)``, which **diverges** at
    the antipode.  The exact kernel does not: ``G(π,ω) = ω/(4π sin πω)``, which
    and *linear in* ``ω``.  In between, everything depends on the single
    combination ``ωe`` with ``e = π − χ``:

    * for ``ωe ≫ 1`` the exact amplitude tracks the WKB ``1/(4π sin e)``;
    * for ``ωe ≲ 1`` it saturates, and the caustic is cut off at ``e* ∼ 1/ω``.

    Measured as a **collapse**: ``exact/WKB = |sin(ωe)/sin(πω)|`` depends on
    ``ωe`` alone, so the ratio at fixed ``ωe`` is the same number at every
    carrier — the statement that ``1/ω`` is the only scale the caustic has.

    Evaluated at **real** ``ω`` and at **half-integer** carriers.  Both choices
    are deliberate: ``ω ∈ ℤ`` are the free ESU eigenfrequencies, where ``G`` has
    a genuine pole, and a contour offset placed there breaks the collapse at
    ``O(ε/ω)`` for a numerical reason.  At half-integers ``|sin πω| = 1``,
    so the saturation is exactly ``ω/4π`` and the collapse ratio is exactly
    ``|sin(ωe)|``.  This is a property of the kernel and needs no contour.
    """
    rows = []
    for c in carriers:
        w = complex(c, 0.0)
        sat = abs(green_omega(math.pi, w))
        for scaled in (2.0, 1.0, 0.5):
            e = scaled / c
            exact = abs(green_omega(math.pi - e, w))
            wkb = 1.0 / (4.0 * math.pi * math.sin(e))
            rows.append({"carrier": float(c), "omega_times_e": scaled,
                         "e": float(e), "exact": float(exact),
                         "wkb": float(wkb), "ratio": float(exact / wkb)})
        rows.append({"carrier": float(c), "omega_times_e": 0.0, "e": 0.0,
                     "exact": float(sat), "wkb": float("inf"),
                     "ratio": 0.0})
    sats = [r["exact"] for r in rows if r["omega_times_e"] == 0.0]
    per_scale: Dict[float, List[float]] = {}
    for r in rows:
        if r["omega_times_e"] > 0.0:
            per_scale.setdefault(r["omega_times_e"], []).append(r["ratio"])
    spread = max(max(v) - min(v) for v in per_scale.values())
    growth = [sats[i + 1] / sats[i] for i in range(len(sats) - 1)]
    steps = [carriers[i + 1] / carriers[i] for i in range(len(carriers) - 1)]
    return {"rows": rows, "saturation_amplitudes": sats,
            "carrier_ratios": steps, "saturation_ratios": growth,
            "worst_collapse_spread": float(spread),
            "the_saturation_is_linear_in_omega": bool(
                all(abs(g - s) < 1e-6 for g, s in zip(growth, steps))),
            "the_ratio_collapses_in_omega_times_e": bool(spread < 1e-6),
            "the_caustic_scale": "e* ~ 1/ω",
            "what_wkb_gets_wrong": ("a divergence where the exact amplitude is "
                                    "finite and proportional to ω")}


def measure_the_low_frequency_limit_recovers_the_tomography(
        n_observations: int = 12, seed: int = 20260901) -> Dict[str, object]:
    """The bridge back to PR #258: ``𝒲 = −β``, out of the time-dependent solve.

    ``∫dt φ(t) = φ̂(ω = 0)``, so the DC weight of the *solved time series* is
    exactly the static kernel PR #258 did its tomography on.  Running that
    protocol on numbers produced by the dynamic solver — least squares for
    ``S = Re R``, then ``𝒲 = S₁₂/det S − G₀`` — has to give ``−β``, and the
    route goes through the whole contour integral rather than through the
    matrix.

    On the retarded contour the accessible integral is ``φ̂(iε)`` rather than
    ``φ̂(0)``.  ``Γ`` is **even** in ``ω``, so the bias is ``O(ε²)`` and two
    contours Richardson-extrapolate it away; both the raw and the extrapolated
    numbers are reported so the correction is visible rather than folded in.
    """
    from .two_source import (green_at, random_points, source_vector,
                             static_response)

    pair = working_pair()
    d = pair.separation
    pulse = GaussianPulse(amplitude=1.0, carrier=0.0, width=0.25)
    grid_hi = RetardedGrid(n=1 << 17, span=600.0, eps=0.08)
    grid_lo = RetardedGrid(n=1 << 17, span=600.0, eps=0.04)

    def dc(x: np.ndarray, y: np.ndarray, grid: RetardedGrid) -> float:
        setup = TwoWaveSetup(pair=pair, source_a=y, source_b=y, observer=x,
                             pulse_a=pulse, pulse_b=pulse, grid=grid)
        sol = solve_field(setup, y, pulse)
        damped = sol["phi"] * np.exp(-grid.eps * grid.times)
        weight = complex(pulse.spectrum(complex(0.0, grid.eps))).real
        return float(np.sum(damped) * grid.dt / weight)

    pts = random_points(2 * int(n_observations), seed=seed)
    rows, mat, rhs = [], [], []
    s_true = static_response(pair)
    for k in range(int(n_observations)):
        x, y = pts[2 * k], pts[2 * k + 1]
        hi, lo = dc(x, y, grid_hi), dc(x, y, grid_lo)
        rich = (4.0 * lo - hi) / 3.0
        exact = float(green_at(geodesic(x, y))
                      + source_vector(x, d) @ s_true @ source_vector(y, d))
        rows.append({"raw_eps_0p08": hi, "raw_eps_0p04": lo,
                     "richardson": rich, "exact": exact,
                     "error": rich - exact})
        va, vb = source_vector(x, d), source_vector(y, d)
        mat.append([va[0] * vb[0], va[0] * vb[1] + va[1] * vb[0],
                    va[1] * vb[1]])
        rhs.append(rich - green_at(geodesic(x, y)))
    sol, *_ = np.linalg.lstsq(np.array(mat), np.array(rhs), rcond=None)
    s_rec = np.array([[sol[0], sol[1]], [sol[1], sol[2]]])
    g0 = float(gamma_omega(0.0 + 1e-9j, d).real[0, 1])
    w = disconnection_defect(s_rec, g0)
    beta = complex(pair.beta).real
    return {"n_observations": int(n_observations),
            "rows": rows[:4],
            "worst_kernel_error": float(max(abs(r["error"]) for r in rows)),
            "worst_response_error": float(np.abs(s_rec - s_true).max()),
            "W_from_the_time_dependent_solve": float(w),
            "minus_beta": float(-beta),
            "W_error": float(abs(w + beta)),
            "the_bridge_closes": bool(abs(w + beta) < 1e-3),
            "what_it_checks": ("the DC content of the contour integral, the "
                               "least-squares recovery, and PR #258's defect, "
                               "end to end")}


# ════════════════════════════════════════════════════════════════════════════
# SUPERPOSITION AND THE INTERFERENCE STRESS TENSOR
# ════════════════════════════════════════════════════════════════════════════
def superpose(*solutions: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    """Add solved fields.  The equation is linear, so this is exact.

    Every entry — the field, its time derivatives, its gradient and its Hessian
    — adds, because each is a linear functional of the same contour integral.
    """
    if not solutions:
        raise ValueError("nothing to superpose")
    keys = solutions[0].keys()
    return {k: sum(s[k] for s in solutions) for k in keys}


def zero_like(solution: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    """A switched-off source, as a solution — for removing one honestly."""
    return {k: np.zeros_like(v) for k, v in solution.items()}


def cross_stress_tensor(sol_a: Dict[str, np.ndarray],
                        sol_b: Dict[str, np.ndarray],
                        index: int) -> np.ndarray:
    """``ΔT_{μν} = T[φ_A + φ_B] − T[φ_A] − T[φ_B]`` — the interference tensor.

    ``T`` is quadratic in the field, so ``ΔT`` is its **bilinear** cross term:
    identically zero when either source is switched off, and traceless whenever
    the two pieces are.  Computed from three separate evaluations of the same
    functional rather than from a hand-derived bilinear form, so the bilinearity
    is a measurement and not an assumption — the same discipline PR #258's
    review imposed on the static cross term.

    This is the object **backreaction** would actually see: the total source is
    ``T[φ_A + φ_B]``, and ``ΔT`` is the part of it that exists only because both
    waves are there.  It is *not* the same diagnostic as ``T_A:T_B``, and
    `measure_the_interference_tensor_is_largest_where_the_invariant_is_null`
    measures how differently they behave.
    """
    total = stress_tensor(superpose(sol_a, sol_b), index)
    return total - stress_tensor(sol_a, index) - stress_tensor(sol_b, index)


def two_leg_channels(observer: Sequence[float], source: Sequence[float],
                     separation: float) -> List[Dict[str, object]]:
    """Every ``(i, j)`` two-leg path, enumerated rather than minimised over.

    ``j`` is the mouth the source drives, ``i`` the mouth the signal leaves
    from, and the delay is ``χ(y,c_j) + χ(c_i,x)``.  The **predicted invariant
    depends only on ``i``** — the arrival direction at the observer is set by
    which mouth the wave emerges from, and the entry leg contributes only a
    delay and a weight.  Listing all four makes that visible; taking a minimum
    over them, as a first draft did, hides it.
    """
    c1, c2 = mouth_positions(separation)
    out: List[Dict[str, object]] = []
    for i, ci in enumerate((c1, c2)):
        for j, cj in enumerate((c1, c2)):
            out.append({
                "exit_mouth": i + 1, "entry_mouth": j + 1,
                "delay": float(geodesic(source, cj) + geodesic(ci, observer)),
                "predicted_invariant": float(
                    (1.0 - arrival_directions(observer, source, ci)) ** 2)})
    return sorted(out, key=lambda r: r["delay"])


def measure_the_cross_mouth_channels_are_labelled_by_the_exit_mouth(
        carrier: float = 60.0, width: float = 0.035) -> Dict[str, object]:
    """An explicit ``(i, j)`` audit — and the ``β = 0`` control that scopes it.

    The four two-leg paths carry **two** distinct predicted invariants, one per
    exit mouth (``0.651935`` and ``0.563669``), so this is a discriminating test
    rather than a single number: the field has to pick the right one at each
    delay.  A short pulse is used so the two extreme channels are clean of
    neighbours with a *different* exit mouth; neighbours sharing an exit mouth
    are harmless, because they arrive from the same direction.

    And then the control that matters, the one PR #258's review taught this arc
    to run first: the same measurement with **``β = 0``**, two *disconnected*
    mouths.  The invariant barely moves at all, and it cannot:
    ``𝒩`` is amplitude-normalized, a single channel is a single direction, and
    ``β`` rescales the channel's weight without touching its geometry.  Swept
    over ``β ∈ [0, 0.26]``, all inside PR #257's cone, ``𝒩`` moves by ``6e-07``
    — a part in ``10⁶``, and **five orders below** the ``0.088`` that separates
    the two exit mouths — while the channel's weight moves by ``0.6%``.  It is
    not exactly zero because the neighbouring channels leak a little into the
    window, and their weights do depend on ``β``; the residual is quoted rather
    than rounded away.

    So the honest statement is the dynamical version of #258's: **this
    observable sees structure at the mouths, not the connection between them.**
    What sees the connection is still ``𝒲 = −β``, from the low-frequency limit
    of the same solve.  The multipath result stands — a second arrival direction
    destroys the collinear null — but the throat's *non-locality* does not
    supply it.
    """
    grid = _fine_grid()
    obs = circle_point(SOURCE_GAP + OBSERVER_REACH)
    src_a, src_b = circle_point(0.0), circle_point(SOURCE_GAP)
    chi_a = geodesic(obs, src_a)
    t_star = 3.0
    channels = two_leg_channels(obs, src_b, WORKING_SEPARATION)

    def at(delay: float, pair: MouthPair) -> Dict[str, float]:
        setup = TwoWaveSetup(
            pair=pair, source_a=src_a, source_b=src_b, observer=obs,
            pulse_a=GaussianPulse(carrier=carrier, width=width,
                                  t0=t_star - chi_a),
            pulse_b=GaussianPulse(carrier=carrier, width=width,
                                  t0=t_star - delay),
            grid=grid, with_throat=True)
        return normalized_invariant(setup, t_star, 2.0 * width)

    # the two channels whose different-exit neighbours are far enough away
    resolvable = [channels[0], channels[-1]]
    a1, a2, _ = WORKING_BOUNDARY
    connected = MouthPair(WORKING_SEPARATION, a1, a2, 0.06)
    disconnected = MouthPair(WORKING_SEPARATION, a1, a2, 0.0)
    rows = []
    for ch in resolvable:
        got = at(ch["delay"], connected)
        ctrl = at(ch["delay"], disconnected)
        rows.append({**ch,
                     "measured_invariant": got["invariant"],
                     "control_beta_zero": ctrl["invariant"],
                     "relative_error": abs(got["invariant"]
                                           - ch["predicted_invariant"])
                     / ch["predicted_invariant"],
                     "beta_shift": abs(got["invariant"] - ctrl["invariant"]),
                     "weight": got["energy_product"],
                     "control_weight": ctrl["energy_product"]})

    sweep = []
    delay = float(channels[0]["delay"])
    base = None
    for beta in (0.0, 0.06, 0.12, 0.20, 0.26):
        pair = MouthPair(WORKING_SEPARATION, a1, a2, beta)
        got = at(delay, pair)
        base = base or got["energy_product"]
        sweep.append({"beta": float(beta),
                      "loewner_margin": float(positivity_defect(
                          pair.boundary_matrix(),
                          WORKING_SEPARATION)["min_eigenvalue"]),
                      "invariant": got["invariant"],
                      "weight_ratio": got["energy_product"] / base})
    spread = max(abs(r["invariant"] - sweep[0]["invariant"]) for r in sweep)
    distinct = sorted({round(c["predicted_invariant"], 9) for c in channels})
    return {"carrier": carrier, "width": width,
            "channels": channels,
            "distinct_predictions": distinct,
            "rows": rows,
            "beta_sweep": sweep,
            "worst_relative_error": float(
                max(r["relative_error"] for r in rows)),
            "the_prediction_depends_only_on_the_exit_mouth": bool(
                len(distinct) == 2),
            "the_field_picks_the_right_one": bool(
                max(r["relative_error"] for r in rows) < 3e-3),
            "worst_beta_shift": float(max(r["beta_shift"] for r in rows)),
            "beta_sweep_spread": float(spread),
            "exit_mouth_separation": float(distinct[-1] - distinct[0]),
            "beta_spread_over_the_signal": float(
                spread / (distinct[-1] - distinct[0])),
            "the_invariant_is_beta_independent": bool(
                spread < 1e-5
                and spread < 1e-4 * (distinct[-1] - distinct[0])),
            "the_weight_moves_instead": float(
                abs(1.0 - sweep[-1]["weight_ratio"])),
            "every_sweep_point_is_inside_the_cone": bool(
                all(r["loewner_margin"] > 0.0 for r in sweep)),
            "the_scope": ("this observable sees structure at the mouths, not "
                          "the connection between them; W = −β from the "
                          "low-frequency limit is what sees the connection")}


def measure_the_interference_tensor_is_largest_where_the_invariant_is_null(
        carrier: float = 24.0, width: float = 0.10) -> Dict[str, object]:
    """``ΔT_{μν}``, and the fact that it and ``T_A:T_B`` disagree completely.

    ``T`` is quadratic, so the two-wave content of the *total* stress tensor is
    the bilinear cross term ``ΔT = T[φ_A+φ_B] − T[φ_A] − T[φ_B]``.  Checked
    rather than assumed: it is traceless, and it is **exactly zero** when either
    source is switched off — PR #253's missing property, now at tensor level and
    obtained by evaluating the same functional three times rather than by
    multiplying anything by zero.

    The result worth the round: for two waves of comparable strength,
    ``ΔT^{00}/√(T_A^{00}T_B^{00})`` reaches its **maximum** ``2`` in the
    *collinear* configuration — two parallel waves add coherently — which is
    precisely the configuration where ``T_A:T_B`` vanishes.  Head-on, where the
    invariant is maximal, the interference energy is roughly half that.

    So the two diagnostics are not interchangeable, and a backreaction estimate
    driven by ``𝒞 = T_A:T_B`` would look at the collinear case, see nothing, and
    be wrong about the source by a factor of order the whole effect.  ``ΔT`` is
    what backreaction integrates; ``T_A:T_B`` is what the collision invariant
    measures.
    """
    grid = _fine_grid()
    src_a, src_b = circle_point(0.0), circle_point(SOURCE_GAP)
    t_star = 3.0
    rows = []
    for name, theta in (("collinear", SOURCE_GAP + OBSERVER_REACH),
                        ("head_on", SOURCE_GAP - OBSERVER_REACH)):
        obs = circle_point(theta)
        chi_a, chi_b = geodesic(obs, src_a), geodesic(obs, src_b)
        setup = TwoWaveSetup(
            pair=working_pair(), source_a=src_a, source_b=src_b, observer=obs,
            pulse_a=GaussianPulse(carrier=carrier, width=width,
                                  t0=t_star - chi_a),
            pulse_b=GaussianPulse(carrier=carrier, width=width,
                                  t0=t_star - chi_b),
            grid=grid, with_throat=False)
        sol_a = solve_field(setup, src_a, setup.pulse_a)
        sol_b = solve_field(setup, src_b, setup.pulse_b)
        ts = grid.times
        idx = np.where((ts > t_star - 2.0 * width)
                       & (ts < t_star + 2.0 * width))[0]
        best = None
        for i in idx:
            ta = stress_tensor(sol_a, int(i))
            tb = stress_tensor(sol_b, int(i))
            u = float(ta[0, 0] * tb[0, 0])
            if best is None or u > best[0]:
                best = (u, int(i))
        i = best[1]
        ta = stress_tensor(sol_a, i)
        tb = stress_tensor(sol_b, i)
        dt = cross_stress_tensor(sol_a, sol_b, i)
        # bilinearity, the honest way: switch a source off and re-evaluate
        removed = cross_stress_tensor(sol_a, zero_like(sol_b), i)
        scale = math.sqrt(float(ta[0, 0] * tb[0, 0]))
        rows.append({
            "configuration": name,
            "invariant": contract_stress(ta, tb) / float(ta[0, 0] * tb[0, 0]),
            "delta_T00": float(dt[0, 0]),
            "normalized_delta_T00": float(dt[0, 0] / scale),
            "max_component": float(np.abs(dt).max()),
            "trace": trace_of(dt),
            "with_a_source_removed": float(np.abs(removed).max())})
    col = [r for r in rows if r["configuration"] == "collinear"][0]
    head = [r for r in rows if r["configuration"] == "head_on"][0]
    return {"carrier": carrier, "width": width, "rows": rows,
            "collinear_invariant": col["invariant"],
            "collinear_interference": col["normalized_delta_T00"],
            "head_on_invariant": head["invariant"],
            "head_on_interference": head["normalized_delta_T00"],
            "worst_trace": float(max(abs(r["trace"]) for r in rows)),
            "worst_value_with_a_source_removed": float(
                max(r["with_a_source_removed"] for r in rows)),
            "delta_T_is_traceless": bool(
                max(abs(r["trace"]) for r in rows) < 1e-12),
            "delta_T_vanishes_when_a_source_is_removed": bool(
                max(r["with_a_source_removed"] for r in rows) == 0.0),
            "the_interference_is_maximal_where_the_invariant_is_null": bool(
                abs(col["normalized_delta_T00"] - 2.0) < 1e-2
                and abs(col["invariant"]) < 1e-5
                and head["normalized_delta_T00"]
                < col["normalized_delta_T00"] * 0.7),
            "the_lesson": ("ΔT and T_A:T_B are different diagnostics: the "
                           "collinear configuration nulls the invariant and "
                           "MAXIMISES the interference energy")}
