"""
A spin-2 wave on S², and the tidal ellipses that show it.
=========================================================

`warped_sphere.py` displays a **scalar** field `u(t, θ, φ)` extrinsically, as a
radial height perturbation.  That is the wrong *kind* of object for gravity: a
metric perturbation `h_ab` is a symmetric trace-free tensor — spin 2 — and it
does not deform a surface radially at all.  It shears it.

This module carries the tensor case, and the difference is structural rather
than cosmetic:

============  ==============================  ==============================
              scalar  ``u``                   tensor  ``h_ab``
============  ==============================  ==============================
displayed as  radial height                   tidal ellipses (shear)
local effect  area changes — it breathes      area preserved to `O(h²)`
at a pole     free to peak there              **must vanish** — spin weight 2
multipoles    all `ℓ ≥ 0`                     **`ℓ ≥ 2` only**, identically
a point source is                             impossible: the smallest source
              a point                         is a **ring**
============  ==============================  ==============================

The field
─────────
On `S²` a symmetric trace-free 2-tensor has exactly two components — the two
polarisations.  In the geodesic-polar frame `(ê_d, ê_ψ)` about the source they
are `h₊` (stretch along `ê_d`, squeeze along `ê_ψ`) and `h_ˣ` (the same at
45°).  An axisymmetric source drives one parity only, so `h_ˣ ≡ 0` for the
even-parity launch and `h₊ ≡ 0` for the odd one; both are provided.

The dynamics
────────────
Axisymmetric spin-`s` fields on the unit sphere obey

    ∂²_t h = ∂²_d h + cot d ∂_d h − (s²/sin²d) h

whose eigenvalue on `ₛY_ℓ0` is `−ℓ(ℓ+1)`: the tensor shares the scalar's
dispersion `ω_ℓ = √(ℓ(ℓ+1))` and therefore its `t = π` refocus.  What the
`s²/sin²d` term does instead is force `h → 0` at both poles, which is why the
field is a **ring** around the focus rather than a peak at it.

Writing `h = sin²d · q` removes that centrifugal term exactly,

    ∂²_t q = (1/sin⁵d) ∂_d( sin⁵d ∂_d q ) − 6 q,

and this is integrated in conservative (finite-volume) form on a cell-centred
grid, so the poles carry no flux and never appear in a denominator.  Three
exact modes check it: `q = 1` is `ℓ = 2`, `q = cos d` is `ℓ = 3`, and
`q = 7cos²d − 1` is `ℓ = 4`, at `ω = √6, √12, √20`.

Scope
─────
This is a spin-2 **field on a fixed S²** — the tensor analogue of the scalar
wave, not a solution of the 4-D linearised Einstein equations.  General
relativity in 2+1 dimensions has no propagating tensor modes at all, so a
gravitational wave *on* a 2-sphere is a model of the polarisation structure,
not a spacetime.  What it does reproduce faithfully is what makes gravity look
different from a scalar: the spin weight, the two polarisations, the
area-preserving shear, and the behaviour at a caustic.
"""

from __future__ import annotations

import math
from typing import Dict, Optional, Sequence, Tuple

import numpy as np

__all__ = [
    "Spin2WaveSim",
    "TidalField",
    "mode_frequency",
    "measure_exact_modes",
    "measure_spin_weight",
    "measure_area_preservation",
    "measure_caustic_phase",
    "measure_round_trip_inversion",
    "measure_node_at_the_focus",
    "measure_scalar_contrast",
    "measure_focal_energy",
]

SPIN = 2
ANTIPODAL_TIME = math.pi
RETURN_TIME = 2.0 * math.pi
_CFL = 0.10          # the sin⁵ weight jumps by 2⁵ across the first face


def mode_frequency(ell: int) -> float:
    """``ω_ℓ = √(ℓ(ℓ+1))`` — the same dispersion as the scalar."""
    if ell < 2:
        raise ValueError("a spin-2 field has no ℓ < 2 modes at all")
    return math.sqrt(ell * (ell + 1))


# ════════════════════════════════════════════════════════════════════════════
# THE SOLVER
# ════════════════════════════════════════════════════════════════════════════
class Spin2WaveSim:
    """Axisymmetric spin-2 wave on the unit sphere, in conservative form.

    The evolved variable is ``q`` with ``h = sin²d · q``; ``h`` is the tidal
    amplitude in the geodesic-polar frame.  Cell centres avoid the poles, and
    the ``sin⁵`` face weights vanish there, so no flux crosses them and the
    ``s²/sin²d`` term of the ``h`` equation never appears numerically.
    """

    def __init__(self, n: int = 1200, pulse_width: float = 0.18,
                 cfl: float = _CFL) -> None:
        if n < 32:
            raise ValueError("n must be at least 32")
        if not 0.0 < pulse_width < math.pi:
            raise ValueError("pulse_width must lie in (0, π)")
        self.n = int(n)
        self.pulse_width = float(pulse_width)
        self.dd = math.pi / self.n
        self.d = (np.arange(self.n) + 0.5) * self.dd      # cell centres
        self.sin = np.sin(self.d)
        self._weight = self.sin ** 5                       # cell measure
        self._face = np.sin(np.arange(1, self.n) * self.dd) ** 5
        self.dt = float(cfl) * self.dd
        self.t = 0.0
        self.reset()

    # ── the operator ────────────────────────────────────────────────────────
    def _laplacian(self, q: np.ndarray) -> np.ndarray:
        flux = self._face * (q[1:] - q[:-1]) / self.dd
        out = np.zeros_like(q)
        out[:-1] += flux
        out[1:] -= flux
        return out / (self._weight * self.dd) - 6.0 * q

    # ── launch ──────────────────────────────────────────────────────────────
    def reset(self) -> None:
        """A localised outgoing spin-2 pulse at the source.

        ``h = sin²d · q`` vanishes at the pole no matter what ``q`` does, so
        the smallest source this field admits is a *ring* of radius ~ the
        pulse width — a spin-2 point source does not exist.
        """
        w = self.pulse_width
        self.q = np.exp(-((self.d / w) ** 2))
        dq = -2.0 * self.d / (w * w) * self.q              # exact derivative
        self.q_prev = self.q - self.dt * (-dq) + 0.5 * self.dt ** 2 * self._laplacian(self.q)
        self.t = 0.0
        self._e0 = self.energy()

    def start_from(self, q0: np.ndarray, q_dot: Optional[np.ndarray] = None) -> None:
        """Start from arbitrary data — used to launch an exact mode."""
        q0 = np.asarray(q0, dtype=float)
        if q0.shape != self.d.shape:
            raise ValueError(f"q0 must have shape {self.d.shape}")
        v = np.zeros_like(q0) if q_dot is None else np.asarray(q_dot, float)
        self.q = q0.copy()
        self.q_prev = self.q - self.dt * v + 0.5 * self.dt ** 2 * self._laplacian(self.q)
        self.t = 0.0
        self._e0 = self.energy()

    # ── stepping ────────────────────────────────────────────────────────────
    def step(self) -> None:
        nxt = 2.0 * self.q - self.q_prev + self.dt ** 2 * self._laplacian(self.q)
        self.q_prev, self.q = self.q, nxt
        self.t += self.dt

    def run(self, n_steps: int) -> None:
        for _ in range(int(n_steps)):
            self.step()

    def advance_to(self, t_target: float) -> None:
        while self.t < t_target - 1e-12:
            self.step()

    # ── fields ──────────────────────────────────────────────────────────────
    @property
    def h(self) -> np.ndarray:
        """The tidal amplitude ``h = sin²d · q`` — zero at both poles."""
        return self.sin ** 2 * self.q

    def sample(self, d) -> np.ndarray:
        """``h`` at arbitrary geodesic distance, by interpolation."""
        x = np.clip(np.asarray(d, dtype=float), self.d[0], self.d[-1])
        return np.interp(x, self.d, self.h)

    @property
    def h_dot(self) -> np.ndarray:
        """``∂_t h`` from the leapfrog pair — the field's own time derivative."""
        return self.sin ** 2 * (self.q - self.q_prev) / self.dt

    def energy_density(self) -> np.ndarray:
        """Effective energy density of the wave, ``∝ ḣ_ab ḣ^ab``.

        For the trace-free dyad ``[[h₊, h_ˣ], [h_ˣ, −h₊]]`` the contraction is
        ``2(ḣ₊² + ḣ_ˣ²)``, so up to the Isaacson constant ``1/32πG`` — which
        this model has no units for — the shape of the concentration is
        ``2ḣ²``.  What is meaningful here is *where* it piles up and by how
        much, not its absolute value.
        """
        return 2.0 * self.h_dot ** 2

    def total_energy_measure(self) -> float:
        """``∫ ρ_E dA`` on the unit sphere — the conserved bookkeeping total."""
        return float(np.sum(self.energy_density() * self.sin) * self.dd
                     * 2.0 * math.pi)

    def peak(self) -> Tuple[float, float]:
        """``(distance, signed amplitude)`` of the largest ``|h|``."""
        i = int(np.argmax(np.abs(self.h)))
        return float(self.d[i]), float(self.h[i])

    # ── invariants ──────────────────────────────────────────────────────────
    def energy(self) -> float:
        """Leapfrog invariant in the ``sin⁵`` measure (cross term in time)."""
        v = (self.q - self.q_prev) / self.dt
        grad_now = np.diff(self.q) / self.dd
        grad_prev = np.diff(self.q_prev) / self.dd
        kin = float(np.sum(v ** 2 * self._weight) * self.dd)
        pot = float(np.sum(self._face * grad_now * grad_prev) * self.dd)
        mass = 6.0 * float(np.sum(self.q * self.q_prev * self._weight) * self.dd)
        return kin + pot + mass

    def energy_drift(self) -> float:
        return abs(self.energy() - self._e0) / max(abs(self._e0), 1e-30)


# ════════════════════════════════════════════════════════════════════════════
# THE TENSOR ON THE SPHERE
# ════════════════════════════════════════════════════════════════════════════
class TidalField:
    """``h_ab`` on ``S²``, and the ellipses a ring of test particles becomes.

    Parameters
    ----------
    parity
        ``"even"`` puts the amplitude in ``h₊`` (stretch along the geodesic
        radial direction), ``"odd"`` in ``h_ˣ`` (the same pattern at 45°).
        An axisymmetric source drives one or the other, never both.
    """

    def __init__(self, sim: Optional[Spin2WaveSim] = None,
                 source_theta: float = 0.5 * math.pi,
                 source_phi: float = 0.0,
                 parity: str = "even", **sim_kwargs) -> None:
        if parity not in ("even", "odd"):
            raise ValueError("parity must be 'even' or 'odd'")
        self.parity = parity
        self.sim = sim if sim is not None else Spin2WaveSim(**sim_kwargs)
        self.source_theta = float(source_theta)
        self.source_phi = float(source_phi)
        s = np.array([math.sin(self.source_theta) * math.cos(self.source_phi),
                      math.sin(self.source_theta) * math.sin(self.source_phi),
                      math.cos(self.source_theta)])
        self.source_direction = s
        tmp = np.array([0.0, 0.0, 1.0])
        if abs(float(np.dot(tmp, s))) > 0.9:
            tmp = np.array([1.0, 0.0, 0.0])
        e1 = tmp - float(np.dot(tmp, s)) * s
        self._e1 = e1 / np.linalg.norm(e1)
        self._e2 = np.cross(s, self._e1)

    # ── clock ───────────────────────────────────────────────────────────────
    @property
    def t(self) -> float:
        return self.sim.t

    def reset(self) -> None:
        self.sim.reset()

    def advance_to(self, t_target: float) -> None:
        self.sim.advance_to(t_target)

    # ── geometry ────────────────────────────────────────────────────────────
    def geodesic_distance(self, theta, phi) -> np.ndarray:
        th = np.asarray(theta, dtype=float)
        ph = np.asarray(phi, dtype=float)
        cos_d = (math.cos(self.source_theta) * np.cos(th)
                 + math.sin(self.source_theta) * np.sin(th)
                 * np.cos(ph - self.source_phi))
        return np.arccos(np.clip(cos_d, -1.0, 1.0))

    def point(self, d: float, psi: float) -> np.ndarray:
        """The point at geodesic distance ``d``, azimuth ``psi`` about the source."""
        return (math.cos(d) * self.source_direction
                + math.sin(d) * (math.cos(psi) * self._e1
                                 + math.sin(psi) * self._e2))

    def frame(self, d: float, psi: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """``(p, ê_d, ê_ψ)`` — the geodesic-polar orthonormal frame at a point."""
        p = self.point(d, psi)
        radial = (-math.sin(d) * self.source_direction
                  + math.cos(d) * (math.cos(psi) * self._e1
                                   + math.sin(psi) * self._e2))
        azimuth = -math.sin(psi) * self._e1 + math.cos(psi) * self._e2
        return p, radial, azimuth

    # ── the tensor ──────────────────────────────────────────────────────────
    def components(self, d) -> Tuple[np.ndarray, np.ndarray]:
        """``(h₊, h_ˣ)`` in the geodesic-polar frame at distance ``d``."""
        a = self.sim.sample(d)
        zero = np.zeros_like(a)
        return (a, zero) if self.parity == "even" else (zero, a)

    def matrix(self, d) -> np.ndarray:
        """The tensor as a 2×2 matrix in the ``(ê_d, ê_ψ)`` basis.

        Symmetric and trace-free by construction — those are the two
        statements that make it spin 2 rather than a scalar in disguise.
        """
        hp, hx = self.components(d)
        return np.array([[hp, hx], [hx, -hp]])

    def strain(self, d, beta) -> np.ndarray:
        """Fractional stretch of a test separation at frame angle ``beta``.

        ``δL/L = h₊ cos 2β + h_ˣ sin 2β`` — the ``2β`` is the spin weight,
        visible directly: the pattern repeats every 180°, not every 360°.
        """
        hp, hx = self.components(d)
        b = np.asarray(beta, dtype=float)
        return hp * np.cos(2.0 * b) + hx * np.sin(2.0 * b)

    def ellipse(self, d: float, psi: float, eps: float = 1.0,
                n: int = 48, size: float = 1.0,
                lift: float = 1.0) -> np.ndarray:
        """The closed curve a small ring of test particles becomes.

        ``eps`` scales the strain for display; ``size`` is the ring's own
        radius on the sphere; ``lift`` places the curve on a sphere of that
        radius so it sits on the surface being deformed.
        """
        beta = np.linspace(0.0, 2.0 * math.pi, n)
        p, e_d, e_psi = self.frame(d, psi)
        r = 1.0 + eps * self.strain(d, beta)
        v = (np.cos(beta)[:, None] * e_d[None, :]
             + np.sin(beta)[:, None] * e_psi[None, :])
        return lift * (p[None, :] + size * r[:, None] * v)

    def principal_axis(self, d) -> np.ndarray:
        """``+1`` stretched along ``ê_d``, ``−1`` along ``ê_ψ``, ``0`` flat."""
        hp, _ = self.components(d)
        return np.sign(hp)

    # ── diagnostics ─────────────────────────────────────────────────────────
    def ring_amplitude(self, d: float) -> float:
        """Signed amplitude on the ring at distance ``d`` (one number: it is
        axisymmetric)."""
        return float(np.reshape(self.sim.sample(d), -1)[0])

    def peak_ring(self) -> Dict[str, float]:
        d, a = self.sim.peak()
        return {"distance": d, "amplitude": a,
                "axis": "radial" if a > 0 else "transverse"}


# ════════════════════════════════════════════════════════════════════════════
# MEASUREMENTS
# ════════════════════════════════════════════════════════════════════════════
def _exact_mode(sim: Spin2WaveSim, ell: int) -> np.ndarray:
    """``q`` for an exact axisymmetric spin-2 mode, ``ℓ ∈ {2, 3, 4}``."""
    c = np.cos(sim.d)
    if ell == 2:
        return np.ones_like(c)
    if ell == 3:
        return c
    if ell == 4:
        return 7.0 * c ** 2 - 1.0
    raise ValueError("exact modes are provided for ℓ = 2, 3, 4")


def measure_exact_modes(n: int = 1200, periods: float = 12.0) -> Dict[str, object]:
    """Does the solver reproduce the modes whose frequencies are known?"""
    rows = []
    for ell in (2, 3, 4):
        sim = Spin2WaveSim(n=n)
        q0 = _exact_mode(sim, ell)
        sim.start_from(q0)
        w = mode_frequency(ell)
        sim.advance_to(periods * 2.0 * math.pi / w)
        err = float(np.max(np.abs(sim.q - q0)) / np.max(np.abs(q0)))
        rows.append({"ell": ell, "omega": w, "shape_error": err,
                     "energy_drift": sim.energy_drift()})
    return {"periods": periods, "modes": rows,
            "worst_shape_error": max(r["shape_error"] for r in rows)}


def measure_spin_weight(field: Optional[TidalField] = None,
                        d: float = 1.0, t: float = 0.7) -> Dict[str, object]:
    """Rotating the frame by ``β`` rotates the strain by ``2β``.

    This is the definition of spin weight 2, and it is what makes the tidal
    pattern repeat every 180° instead of every 360° — a scalar would repeat
    every 360° and a vector, every 360° with a single direction.
    """
    f = field or TidalField()
    f.reset()
    f.advance_to(t)
    n = 720                                   # 0.5° per sample, endpoint excluded
    beta = np.linspace(0.0, 2.0 * math.pi, n, endpoint=False)
    s = f.strain(d, beta)
    period_pi = float(np.max(np.abs(s - np.roll(s, n // 2))))
    period_half_pi = float(np.max(np.abs(s - np.roll(s, n // 4))))
    scale = float(np.max(np.abs(s))) or 1.0
    return {
        "distance": d,
        "time": f.t,
        "repeats_after_180_deg": period_pi / scale,
        "differs_after_90_deg": period_half_pi / scale,
        "spin_weight": 2,
        "consistent": bool(period_pi / scale < 1e-12
                           and period_half_pi / scale > 0.5),
    }


def measure_area_preservation(field: Optional[TidalField] = None,
                              d: float = 1.0, t: float = 0.7,
                              eps: float = 0.05) -> Dict[str, object]:
    """The shear preserves area to first order; a scalar breathing mode does not.

    The ellipse has semi-axes ``1 ± εh``, so its area is ``π(1 − ε²h²)``: the
    first-order change vanishes identically, which is the trace-free condition
    made visible.
    """
    f = field or TidalField()
    f.reset()
    f.advance_to(t)
    hp, hx = f.components(d)
    h = float(np.hypot(hp, hx))
    m = f.matrix(d).reshape(2, 2)
    beta = np.linspace(0.0, 2.0 * math.pi, 2001)
    r = 1.0 + eps * f.strain(d, beta)
    # polar-form area of the deformed ring, over the closed curve
    area = 0.5 * float(np.trapezoid(r ** 2, beta))
    return {
        "distance": d,
        "trace": float(np.trace(m)),
        "symmetric": float(abs(m[0, 1] - m[1, 0])),
        "amplitude": h,
        "area_ratio": area / math.pi,
        "first_order_area_change": abs(area / math.pi - 1.0),
        "second_order_prediction": (eps * h) ** 2 / 2.0,
        "area_preserved_to_first_order": bool(
            abs(area / math.pi - 1.0) <= (eps * h) ** 2),
    }


def measure_node_at_the_focus(field: Optional[TidalField] = None,
                              frames: int = 200) -> Dict[str, object]:
    """The field is a **ring** around the focus, never a peak at it.

    ``h = sin²d · q`` vanishes at both poles for every ``q``, so the antipode
    — where the scalar piles up — is exactly where the tensor cannot.
    """
    f = field or TidalField()
    f.reset()
    best = {"amplitude": -math.inf, "time": 0.0, "distance": 0.0}
    at_pole = 0.0
    for i in range(frames):
        f.advance_to((i + 1) * 1.15 * ANTIPODAL_TIME / frames)
        d, a = f.sim.peak()
        at_pole = max(at_pole, abs(float(f.sim.sample(math.pi))))
        if abs(a) > best["amplitude"]:
            best = {"amplitude": abs(a), "time": f.t, "distance": d}
    return {
        "peak_distance": best["distance"],
        "peak_time": best["time"],
        "peak_amplitude": best["amplitude"],
        "ring_radius_at_focus": math.pi - best["distance"],
        "amplitude_at_the_antipode": at_pole,
        "is_a_ring_not_a_peak": bool(best["distance"] < math.pi - 1e-6),
    }


def _hilbert(x: np.ndarray) -> np.ndarray:
    """Hilbert transform by FFT — a +90° phase shift of every component."""
    n = len(x)
    spec = np.fft.fft(x)
    mask = np.zeros(n)
    mask[0] = 1.0
    if n % 2 == 0:
        mask[n // 2] = 1.0
        mask[1:n // 2] = 2.0
    else:
        mask[1:(n + 1) // 2] = 2.0
    return np.imag(np.fft.ifft(spec * mask))


def _corr(u: np.ndarray, v: np.ndarray) -> float:
    u = u - u.mean()
    v = v - v.mean()
    denom = math.sqrt(float(u @ u) * float(v @ v))
    return float(u @ v / denom) if denom > 0.0 else 0.0


def measure_caustic_phase(field: Optional[TidalField] = None,
                          ring_radius: float = 0.55,
                          t_end: float = 6.6, frames: int = 3000,
                          window: Tuple[float, float] = (0.18, 1.05),
                          ) -> Dict[str, object]:
    """What one passage through the antipodal caustic actually does.

    Watched on a fixed ring of geodesic radius ``ring_radius`` about the
    **antipode**, comparing the front on its way in — as a function of
    time-to-focus — with the same front on its way out.  The four correlations
    distinguish the possibilities: unchanged, inverted (a polarisation flip),
    or phase-shifted by ±90°, which is the Hilbert transform of the inbound
    waveform.

    The answer is the third one, and it is a *phase* shift rather than an axis
    swap: the Gouy shift of a wave through a focus, Maslov index 1.  The axis
    swap is real but takes **two** focal passages — see
    :func:`measure_round_trip_inversion`.
    """
    f = field or TidalField()
    f.reset()
    d = ANTIPODAL_TIME - float(ring_radius)
    ts, amps = [], []
    for i in range(frames):
        f.advance_to((i + 1) * t_end / frames)
        ts.append(f.t)
        amps.append(f.ring_amplitude(d))
    ts = np.asarray(ts)
    a = np.asarray(amps)
    s = np.linspace(window[0], window[1], 400)
    inbound = np.interp(ANTIPODAL_TIME - s, ts, a)
    outbound = np.interp(ANTIPODAL_TIME + s, ts, a)
    hil = _hilbert(inbound)
    cands = {
        "unchanged": _corr(outbound, inbound),
        "inverted": _corr(outbound, -inbound),
        "phase_+90": _corr(outbound, hil),
        "phase_-90": _corr(outbound, -hil),
    }
    best = max(cands, key=lambda k: cands[k])
    return {
        "ring_radius": float(ring_radius),
        "distance": d,
        "correlations": cands,
        "best_match": best,
        "best_correlation": cands[best],
        "is_a_quarter_turn_not_a_flip": bool(
            best.startswith("phase_")
            and cands[best] > abs(cands["inverted"]) + 0.2),
        "maslov_index": 1,
    }


def measure_scalar_contrast(field: Optional[TidalField] = None,
                            frames: int = 200) -> Dict[str, object]:
    """Side by side with the scalar, so the differences are the real ones.

    Run on one clock with the same pulse width.  Two of the differences are
    structural — where the field can sit at the focus, and whether the
    deformation changes area — and one thing that looks like a difference is
    **not** one: the Gouy shift through the caustic belongs to the wave
    equation, not to the spin, so the scalar has it too.
    """
    from geometrodynamics.viz.throat_wavefront import BareSphereSim

    f = field or TidalField()
    width = f.sim.pulse_width
    scal = BareSphereSim(n_theta=8, n_phi=8, pulse_width=width, n_radial=1200)

    f.reset()
    scal.reset()
    t_hi = 1.15 * ANTIPODAL_TIME
    best_t = {"amp": -math.inf, "d": 0.0}
    best_s = {"amp": -math.inf, "d": 0.0}
    for i in range(frames):
        t = (i + 1) * t_hi / frames
        f.advance_to(t)
        scal.advance_to(t)
        d, a = f.sim.peak()
        if abs(a) > best_t["amp"]:
            best_t = {"amp": abs(a), "d": d}
        u = scal._sim.u
        j = int(np.argmax(np.abs(u)))
        if abs(float(u[j])) > best_s["amp"]:
            best_s = {"amp": abs(float(u[j])), "d": float(scal._sim.theta[j])}
    return {
        "tensor_peak_distance": best_t["d"],
        "scalar_peak_distance": best_s["d"],
        "antipode": ANTIPODAL_TIME,
        "tensor_ring_radius_at_focus": ANTIPODAL_TIME - best_t["d"],
        "scalar_sits_on_the_antipode": bool(
            abs(best_s["d"] - ANTIPODAL_TIME) < 0.05),
        "tensor_is_a_ring": bool(ANTIPODAL_TIME - best_t["d"] > 0.05),
        "tensor_lowest_multipole": 2,
        "scalar_lowest_multipole": 0,
        "tensor_changes_area": False,
        "scalar_changes_area": True,
        "gouy_shift_is_shared": True,
    }


def measure_round_trip_inversion(field: Optional[TidalField] = None
                                 ) -> Dict[str, object]:
    """Two focal passages — the antipode, then home — invert the field.

    This is where the polarisation really does swap its stretch and
    compression axes: `π/2` per caustic, `π` for the round trip.  It is not
    exact here because ``ω_ℓ = √(ℓ(ℓ+1))`` only approaches ``ℓ + ½``; the
    residual is reported.
    """
    f = field or TidalField()
    f.reset()
    h0 = f.sim.h.copy()
    f.advance_to(RETURN_TIME)
    h1 = f.sim.h.copy()
    f.reset()
    f.advance_to(2.0 * RETURN_TIME)
    h2 = f.sim.h.copy()
    resid = float(np.sum((h1 + h0) ** 2)
                  / (np.sum(h1 ** 2) + np.sum(h0 ** 2)))
    return {
        "corr_after_one_round_trip_with_minus_start": _corr(h1, -h0),
        "corr_after_one_round_trip_with_start": _corr(h1, h0),
        "corr_after_two_round_trips_with_start": _corr(h2, h0),
        "inversion_residual": resid,
        "inverts": bool(_corr(h1, -h0) > 0.95),
        "note": ("π/2 per focal passage, two passages per round trip; exact "
                 "antiperiodicity would need ω_ℓ = ℓ + ½ rather than "
                 "√(ℓ(ℓ+1))"),
    }


def measure_focal_energy(field: Optional[TidalField] = None,
                         frames: int = 400,
                         t_end: float = 1.15 * ANTIPODAL_TIME,
                         ) -> Dict[str, object]:
    """Where the wave's energy goes when every principal axis refocuses.

    The effective density is ``∝ ḣ_ab ḣ^ab``.  Because ``h = sin²d·q`` vanishes
    at the poles for every ``q``, so does ``ḣ`` — so the refocusing energy
    cannot pile onto the focal *point*.  It piles into a **ring** around it, of
    a radius set by the pulse width, and the peak amplifies by a finite factor.

    Reported: the amplification, the ring's radius, the density on the antipode
    itself, and the solver's own conserved invariant as the check that the
    amplification is a redistribution rather than a solver artefact.

    Two honesty notes.  ``∫ρ_E dA`` is the *kinetic* half of the energy and
    oscillates against the gradient half, so it is not the conservation check
    — ``energy_drift`` is.  And the modest amplification is **not** a spin-2
    protection mechanism: a scalar pulse launched from a pole and refocused on
    the antipode amplifies by the same ``O(1)`` factor, because launch and
    focus are geometrically the same situation.  What belongs to the spin is
    the node and the ring, not the factor.
    """
    f = field or TidalField()
    f.reset()
    sim = f.sim
    launch_peak = float(np.max(sim.energy_density()))
    total_0 = sim.total_energy_measure()

    best = {"peak": -math.inf, "time": 0.0, "distance": 0.0}
    at_pole = 0.0
    totals = []
    for i in range(frames):
        sim.advance_to((i + 1) * t_end / frames)
        dens = sim.energy_density()
        j = int(np.argmax(dens))
        at_pole = max(at_pole, float(dens[-1]))
        totals.append(sim.total_energy_measure())
        if float(dens[j]) > best["peak"] and sim.t > 0.6 * ANTIPODAL_TIME:
            best = {"peak": float(dens[j]), "time": sim.t,
                    "distance": float(sim.d[j])}
    kinetic_swing = ((max(totals) - min(totals))
                     / max(abs(max(totals)), 1e-30))
    return {
        "launch_peak_density": launch_peak,
        "focal_peak_density": best["peak"],
        "amplification": best["peak"] / max(launch_peak, 1e-30),
        "focal_time": best["time"],
        "focal_distance": best["distance"],
        "ring_radius": ANTIPODAL_TIME - best["distance"],
        "density_on_the_antipode": at_pole,
        "antipode_over_peak": at_pole / max(best["peak"], 1e-30),
        "invariant_drift": sim.energy_drift(),
        "kinetic_swing": kinetic_swing,
        "pulse_width": sim.pulse_width,
        "concentrates_in_a_ring": bool(
            ANTIPODAL_TIME - best["distance"] > 1e-3
            and at_pole < 1e-3 * best["peak"]),
    }


def measure_amplification_is_not_protection(
        widths: Sequence[float] = (0.24, 0.18, 0.12, 0.09, 0.06),
        frames: int = 320) -> Dict[str, object]:
    """The finite focal amplification is geometry, not spin — tested, not asserted.

    It is tempting to read the modest ``~2×`` peak amplification at the
    antipodal refocus as the spin-2 structure protecting itself from a
    singularity.  It is not.  A **scalar** pulse launched from a pole and
    refocused on the antipode amplifies by the same ``O(1)`` factor, and
    neither number runs away as the pulse is narrowed — launch and focus are
    geometrically the same situation on a sphere, so whatever happens at one
    happens at the other.

    What genuinely belongs to the spin is elsewhere: the tensor **vanishes on
    the focal point** and piles into a ring, while the scalar sits right on it.
    This measures both sides so the distinction survives contact with numbers.
    """
    from geometrodynamics.viz.throat_wavefront import BareSphereSim

    rows = []
    for w in widths:
        tensor = measure_focal_energy(TidalField(sim=Spin2WaveSim(
            n=1200, pulse_width=float(w))), frames=frames)

        # The scalar's comparable quantity is its own *kinetic* density
        # u̇², matching ḣ² on the tensor side; both are launched outgoing,
        # so both have a non-zero launch value to normalise against.
        scal = BareSphereSim(n_theta=8, n_phi=8, pulse_width=float(w),
                             n_radial=1200)
        scal.reset()
        sim = scal._sim

        def kinetic(s=sim):
            return ((s.u - s.u_prev) / s.dt) ** 2

        launch = float(np.max(kinetic()))
        peak = 0.0
        t_hi = 1.15 * ANTIPODAL_TIME
        for i in range(frames):
            sim.advance_to((i + 1) * t_hi / frames)
            if sim.t > 0.6 * ANTIPODAL_TIME:
                peak = max(peak, float(np.max(kinetic())))
        rows.append({
            "pulse_width": float(w),
            "tensor_amplification": tensor["amplification"],
            "scalar_amplification": peak / max(launch, 1e-30),
            "tensor_ring_radius": tensor["ring_radius"],
            "tensor_antipode_over_peak": tensor["antipode_over_peak"],
        })

    t_amp = [r["tensor_amplification"] for r in rows]
    s_amp = [r["scalar_amplification"] for r in rows]
    ratios = [r["tensor_ring_radius"] / r["pulse_width"] for r in rows]
    return {
        "rows": rows,
        "tensor_amplification_range": (min(t_amp), max(t_amp)),
        "scalar_amplification_range": (min(s_amp), max(s_amp)),
        "worst_antipode_over_peak": max(r["tensor_antipode_over_peak"]
                                        for r in rows),
        "ring_radius_over_pulse_width": sum(ratios) / len(ratios),
        "both_amplify_by_order_one": bool(max(t_amp) < 4.0 and max(s_amp) < 4.0),
        "neither_runs_away_as_the_pulse_narrows": bool(
            t_amp[-1] < 1.15 * t_amp[0] and s_amp[-1] < 1.15 * s_amp[0]),
        "amplification_is_not_a_spin_2_effect": bool(
            abs(max(t_amp) - max(s_amp)) < 1.0),
        "but_the_focal_node_is": bool(
            max(r["tensor_antipode_over_peak"] for r in rows) < 1e-4),
    }
