"""Angular kicks a scalar wave hands to the two mouths of a glued handled sphere.

The surface is ``viz.klein_handles.build_surface`` (#338): a unit sphere with one
antipodal mouth pair, A at ``+AXIS`` and B at ``-AXIS``, joined by the finite
catenoid neck ``r(s) = sqrt(s^2 + b^2)``. The gluing is orientation-reversing
(Klein, rim angle psi -> pi - psi) or orientation-preserving (torus, psi -> psi + pi).
The field obeys ``u_tt = Delta u - m2 u`` on the intrinsic cotangent mesh with
lumped mass. Nothing here is fitted, and no quantization rule is used.

Kick bookkeeping
----------------
Rotations are Killing fields on the sphere band (the sphere minus both mouth
caps). For each generator ``xi_k = e_k x p`` the density ``l_k = -u_t (xi_k . grad u)``
obeys ``d_t l_k + div S_k = 0`` with

    S_k = grad u (xi_k . grad u) - xi_k (|grad u|^2 + m2 u^2 - u_t^2) / 2 .

The kick on a mouth is the angular momentum the band hands to its rim. For a
weight ``chi`` that is 1 on that rim and 0 near the other mouth,

    J(T) = int_0^T dt int grad chi . S_k dA  -  [ int chi l_k dA ]_0^T ,

which needs only area integrals. The geometry is held fixed, so this is the
leading-order kick. Each mouth is treated as its own body: neck stiffness is not
modelled. An angular impulse J moves A (at +AXIS) along ``J x AXIS`` and B (at
-AXIS) along ``J x (-AXIS)``.

From antipodes the separation is ``pi - |a + b|`` to first order in the two
displacements, so any kick that is not a rigid rotation brings the mouths
closer. Whether the mouths close in together is the sign of ``a . b``:
``follow = a . b / |a|^2 > 0`` means B moves the same way as A.

Flux through the throat
-----------------------
A massless field (``m2 = 0``) admits a static harmonic field whose potential has
a period around the loop through the handle and back over the sphere: Wheeler's
"charge without charge" for a scalar potential. ``FluxField`` gives it in closed
form. Field lines leave one mouth into the sphere and enter the other, so the
mouths act as a source and a sink. The kick then gains a term bilinear in that
field and the wave, which flips sign with either one; ``RimKicks.cross_rates``
evaluates it.

Clock offsets
-------------
``offset_history`` runs #338's ``HandleWaveSim`` with a mouth clock offset (the
Morris-Thorne-Yurtsever ingredient): a signal crossing the neck from A to B
emerges ``offset`` later in external time, earlier when the offset is negative.
The band is static in external time, so the same kick bookkeeping applies to
the external-time frames.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
import scipy.sparse as sp

from geometrodynamics.viz.klein_handles import HandledSphere, HandleSpec, HandleWaveSim, build_surface

AXIS = np.array([1.0, 0.0, 0.0])
MOUTH_ANGLE = 0.3


def smootherstep(q):
    q = np.clip(q, 0.0, 1.0)
    return q * q * q * (10 + q * (-15 + 6 * q))


def glued_surface(handle: str = "klein", n_sphere: int = 24000, n_ring: int = 48, ds: float = 0.012,
                  a: float = MOUTH_ANGLE) -> HandledSphere:
    """#338's sphere with one antipodal handle at ±x; ``handle`` is 'klein', 'torus' or 'none'."""
    if handle not in ("klein", "torus", "none"):
        raise ValueError("handle must be 'klein', 'torus' or 'none'")
    handles = [] if handle == "none" else [HandleSpec(tuple(AXIS), (0.0, 0.0, 1.0), handle == "klein")]
    return build_surface(a=a, handles=handles, n_sphere=n_sphere, n_ring=n_ring, ds_target=ds)


def mouth_frame(surf: HandledSphere):
    """(e1, e2): the rim-angle frame of the handle; e1 is the Klein gluing's mirror axis."""
    if surf.drawn:
        return surf.drawn[0].e1, surf.drawn[0].e2
    e1 = np.array([0.0, 0.0, 1.0])
    return e1, np.cross(AXIS, e1)


class RimKicks:
    """Band operators and rim weights for the angular kick bookkeeping."""

    def __init__(self, surf: HandledSphere, m2: float = 0.25, weight=(0.35, 0.9)):
        self.surf, self.m2 = surf, float(m2)
        P = surf.points
        T = surf.tris[surf.sphere_tris]
        p0, p1, p2 = P[T[:, 0]], P[T[:, 1]], P[T[:, 2]]
        nrm = np.cross(p1 - p0, p2 - p0)
        twice = np.linalg.norm(nrm, axis=1)
        nh = nrm / twice[:, None]
        self.area = twice / 2
        grads = [np.cross(nh, P[T[:, (i + 2) % 3]] - P[T[:, (i + 1) % 3]]) / twice[:, None] for i in range(3)]
        n, nt = len(P), len(T)
        rows = np.repeat(np.arange(nt), 3)
        self.G = [sp.csr_matrix((np.stack([grads[i][:, c] for i in range(3)], 1).ravel(), (rows, T.ravel())),
                                shape=(nt, n)) for c in range(3)]
        self.avg = sp.csr_matrix((np.full(3 * nt, 1 / 3), (rows, T.ravel())), shape=(nt, n))
        cen = (p0 + p1 + p2) / 3
        cen /= np.linalg.norm(cen, axis=1)[:, None]
        self.centroid = cen
        self.xi = np.stack([np.cross(np.eye(3)[k], cen) for k in range(3)])      # (3, nt, 3)
        self.weights = [self.weight(c, *weight) for c in (AXIS, -AXIS)]
        self.gchi = [self.grad(w) for w in self.weights]
        self.cbar = [self.avg @ w for w in self.weights]

    def weight(self, centre, th1: float, th2: float) -> np.ndarray:
        x = self.surf.points / np.linalg.norm(self.surf.points, axis=1)[:, None]
        th = np.arccos(np.clip(x @ centre, -1, 1))
        return np.where(self.surf.on_sphere, 1 - smootherstep((th - th1) / (th2 - th1)), 0.0)

    def grad(self, f: np.ndarray) -> np.ndarray:
        """P1 gradient of vertex values on each band triangle, as ambient vectors (nt, 3)."""
        return np.stack([G @ f for G in self.G], 1)

    def rates(self, u: np.ndarray, v: np.ndarray):
        """Quadratic part: weighted fluxes (2, 3), weighted contents (2, 3) and band total (3,)."""
        g, ub, vb = self.grad(u), self.avg @ u, self.avg @ v
        xg = np.einsum("kti,ti->kt", self.xi, g)
        lag = 0.5 * (np.einsum("ti,ti->t", g, g) + self.m2 * ub * ub - vb * vb)
        ell = -vb[None, :] * xg
        flux = np.array([(np.einsum("ti,ti->t", g, gc)[None, :] * xg
                          - np.einsum("kti,ti->kt", self.xi, gc) * lag[None, :]) @ self.area for gc in self.gchi])
        content = np.array([(ell * cb[None, :]) @ self.area for cb in self.cbar])
        return flux, content, ell @ self.area

    def cross_rates(self, u: np.ndarray, v: np.ndarray, gs: np.ndarray):
        """Terms bilinear in the wave (u, v) and a static field with band gradient ``gs`` (m2 = 0)."""
        if self.m2 != 0:
            raise ValueError("a static field with a period needs m2 = 0")
        g, vb = self.grad(u), self.avg @ v
        xg = np.einsum("kti,ti->kt", self.xi, g)
        xs = np.einsum("kti,ti->kt", self.xi, gs)
        dot = np.einsum("ti,ti->t", gs, g)
        ell = -vb[None, :] * xs
        flux = np.array([(np.einsum("ti,ti->t", gs, gc)[None, :] * xg + np.einsum("ti,ti->t", g, gc)[None, :] * xs
                          - np.einsum("kti,ti->kt", self.xi, gc) * dot[None, :]) @ self.area for gc in self.gchi])
        content = np.array([(ell * cb[None, :]) @ self.area for cb in self.cbar])
        return flux, content, ell @ self.area


class KickLedger:
    """Trapezoid accumulation of rim kicks along a sampled (u, u_t) history."""

    def __init__(self, kicks: RimKicks, static: Optional[np.ndarray] = None):
        self.k, self.static = kicks, static
        self.t = None
        self.out = dict(t=[], J=[], Jx=[], L=[], Lx=[])

    def push(self, t: float, u: np.ndarray, v: np.ndarray, record: bool = True):
        f, c, L = self.k.rates(u, v)
        fx = cx = Lx = None
        if self.static is not None:
            fx, cx, Lx = self.k.cross_rates(u, v, self.static)
        if self.t is None:
            self.c0, self.int = c, np.zeros_like(c)
            if fx is not None:
                self.cx0, self.intx = cx, np.zeros_like(cx)
        else:
            dt = t - self.t
            self.int = self.int + 0.5 * dt * (self.f + f)
            if fx is not None:
                self.intx = self.intx + 0.5 * dt * (self.fx + fx)
        self.t, self.f, self.fx = t, f, fx
        if record:
            self.out["t"].append(t)
            self.out["J"].append(self.int - (c - self.c0))
            self.out["L"].append(L)
            if fx is not None:
                self.out["Jx"].append(self.intx - (cx - self.cx0))
                self.out["Lx"].append(Lx)

    def history(self) -> "KickHistory":
        o = {k: np.array(v) for k, v in self.out.items()}
        return KickHistory(t=o["t"], J=o["J"], Jx=o["Jx"] if len(o["Jx"]) else None, L=o["L"],
                           Lx=o["Lx"] if len(o["Lx"]) else None)


@dataclass
class KickHistory:
    """Cumulative rim kicks: ``J[i, m]`` is mouth m's (0 = A, 1 = B) angular impulse at ``t[i]``."""

    t: np.ndarray
    J: np.ndarray
    Jx: Optional[np.ndarray]
    L: np.ndarray
    Lx: Optional[np.ndarray] = None
    energy: Optional[np.ndarray] = None

    def at(self, T: float) -> int:
        return int(np.argmin(np.abs(self.t - T)))


class GluedWave:
    """Direct leapfrog evolution on the glued surface, with the kick ledger run every step."""

    def __init__(self, surf: HandledSphere, m2: float = 0.25, weight=(0.35, 0.9)):
        self.surf, self.m2 = surf, float(m2)
        self.K, self.M = surf.operators()
        self.kicks = RimKicks(surf, m2, weight)
        self.e1, self.e2 = mouth_frame(surf)

    def energy(self, u: np.ndarray, v: np.ndarray) -> float:
        return 0.5 * float(v @ (self.M * v) + u @ (self.K @ u) + self.m2 * u @ (self.M * u))

    def packet(self, d0: float, psi0: float, mode: str, width: float = 0.25):
        """Unit-energy single-signed Gaussian bump ``d0`` from A on side ``u(psi0)``.

        ``mode``: 'release' starts at rest; 'in' travels along the great circle toward A;
        'out' travels away from A, toward B. Returns (u, v, u_dir).
        """
        u_dir = np.cos(psi0) * self.e1 + np.sin(psi0) * self.e2
        s0 = np.cos(d0) * AXIS + np.sin(d0) * u_dir
        away = -np.sin(d0) * AXIS + np.cos(d0) * u_dir
        x = self.surf.points / np.linalg.norm(self.surf.points, axis=1)[:, None]
        c = np.clip(x @ s0, -1, 1)
        d = np.arccos(c)
        u = np.where(self.surf.on_sphere, np.exp(-(d / width) ** 2), 0.0)
        v = np.zeros_like(u)
        if mode in ("in", "out"):
            e = -away if mode == "in" else away
            xi = np.cross(np.cross(s0, e), x)
            gd = -(s0[None, :] - c[:, None] * x) / np.sqrt(np.maximum(1 - c * c, 1e-30))[:, None]
            gu = (-2 * d / width ** 2 * u)[:, None] * gd
            v = np.where(self.surf.on_sphere, -np.einsum("ij,ij->i", xi, gu), 0.0)
        elif mode != "release":
            raise ValueError("mode must be 'release', 'in' or 'out'")
        scale = 1 / np.sqrt(self.energy(u, v))
        return u * scale, v * scale, u_dir

    def ring(self, d0: float, width: float = 0.25):
        """Unit-energy bump on the circle at angle ``d0`` from A: symmetric about the A-B axis."""
        x = self.surf.points / np.linalg.norm(self.surf.points, axis=1)[:, None]
        th = np.arccos(np.clip(x @ AXIS, -1, 1))
        u = np.where(self.surf.on_sphere, np.exp(-((th - d0) / width) ** 2), 0.0)
        v = 0 * u
        return u / np.sqrt(self.energy(u, v)), v

    def run(self, u, v, t_end: float, static: Optional[np.ndarray] = None, record_dt: float = 0.02,
            cfl: float = 0.6) -> KickHistory:
        acc = lambda w: -(self.K @ w) / self.M - self.m2 * w
        x = np.random.default_rng(0).standard_normal(len(u))
        lam = 0.0
        for _ in range(60):
            y = -acc(x)
            lam = np.linalg.norm(y) / np.linalg.norm(x)
            x = y / np.linalg.norm(y)
        n_rec = int(round(t_end / record_dt))
        sub = int(np.ceil(record_dt / (cfl / np.sqrt(1.05 * lam))))
        dt = record_dt / sub
        u, v = np.array(u, float), np.array(v, float)
        a = acc(u)
        ledger = KickLedger(self.kicks, static)
        ledger.push(0.0, u, v)
        energy = [self.energy(u, v)]
        for k in range(n_rec):
            for j in range(sub):
                v += 0.5 * dt * a
                u += dt * v
                a = acc(u)
                v += 0.5 * dt * a
                ledger.push((k * sub + j + 1) * dt, u, v, record=j == sub - 1)
            energy.append(self.energy(u, v))
        h = ledger.history()
        h.energy = np.array(energy)
        return h


@dataclass
class FluxField:
    """Static harmonic field with flux ``2 pi c`` through the throat (massless field only).

    On the band ``phi = c log tan(theta_A / 2)``, with theta_A the angle from mouth A;
    in the neck ``phi = c asinh(s / b) + const``. For c > 0, ``grad phi`` points from A
    toward B over the sphere and back from B to A through the neck: A is the source
    mouth and B the sink. The potential gains ``period`` once around the handle loop.
    Because ``asinh(cot a) = log cot(a/2)``, the neck carries exactly half of it.
    The field is axisymmetric, so it is the same for Klein and torus gluings.
    """

    c: float = 1.0
    a: float = MOUTH_ANGLE

    @property
    def period(self) -> float:
        return 4 * self.c * np.log(1 / np.tan(self.a / 2))

    @property
    def energy(self) -> float:
        """(1/2) int |grad phi|^2 over band and neck: pi c^2 times the period per unit c."""
        return np.pi * self.c * self.period

    def band_potential(self, surf: HandledSphere) -> np.ndarray:
        x = surf.points / np.linalg.norm(surf.points, axis=1)[:, None]
        th = np.arccos(np.clip(x @ AXIS, -1, 1))
        with np.errstate(divide="ignore", invalid="ignore"):
            phi = self.c * np.log(np.tan(th / 2))
        return np.where(surf.on_sphere, phi, 0.0)

    def band_gradient(self, kicks: RimKicks) -> np.ndarray:
        return kicks.grad(self.band_potential(kicks.surf))

    @classmethod
    def unit_energy(cls, a: float = MOUTH_ANGLE) -> "FluxField":
        return cls(c=1 / np.sqrt(4 * np.pi * np.log(1 / np.tan(a / 2))), a=a)


def release_source(psi0: float, d0: float = 1.0, surf: Optional[HandledSphere] = None):
    """Sphere point ``d0`` from A on side u(psi0), in the handle's rim frame."""
    e1, e2 = mouth_frame(surf) if surf is not None else (np.array([0.0, 0.0, 1.0]), np.array([0.0, -1.0, 0.0]))
    return np.cos(d0) * AXIS + np.sin(d0) * (np.cos(psi0) * e1 + np.sin(psi0) * e2)


def offset_history(surf: HandledSphere, source, offsets: Sequence[float] = (), width: float = 0.25,
                   t_max: float = 2 * np.pi + 0.5, frame_dt: float = 0.0025, flux: Optional[FluxField] = None,
                   weight=(0.35, 0.9)):
    """Release a bump at ``source`` with #338's HandleWaveSim and mouth clock offsets.

    The field is massless (HandleWaveSim solves u_tt = Delta u). Kicks are taken from
    the external-time frames on the band, with u_t by central differences, and are
    normalised to unit wave energy (the cross term to unit wave and unit flux energy).
    Frames ``frame_dt`` apart reproduce the direct leapfrog kicks to 0.5%. ``info`` also
    gives the first external time the field at each rim exceeds HandleWaveSim's arrival
    threshold (0.5% of the launch amplitude). Returns (KickHistory, info).
    """
    offsets = tuple(float(o) for o in offsets)
    span = t_max + max([0.0] + [-o for o in offsets]) + 0.1
    sim = HandleWaveSim(surf, source=source, width=width, t_max=t_max, n_frames=int(np.ceil(span / frame_dt)) + 1,
                        offsets=offsets if any(offsets) else None)
    sim.run()
    K, _ = surf.operators()
    u0 = sim.initial()
    e0 = 0.5 * float(u0 @ (K @ u0))
    n = surf.n_vertices
    U = sim.frames[:, :n].astype(float)
    t = np.asarray(sim.times, float)
    V = np.gradient(U, t, axis=0, edge_order=2)
    kicks = RimKicks(surf, m2=0.0, weight=weight)
    static = flux.band_gradient(kicks) / np.sqrt(flux.energy) if flux is not None else None
    ledger = KickLedger(kicks, static)
    for i in range(len(t)):
        ledger.push(float(t[i]), U[i], V[i])
    h = ledger.history()
    h.J = h.J / e0
    h.L = h.L / e0
    if h.Jx is not None:
        h.Jx = h.Jx / np.sqrt(e0)
        h.Lx = h.Lx / np.sqrt(e0)
    ra, rb = surf.rims[0] if surf.rims else (np.array([], int), np.array([], int))
    first = lambda ids: float(np.nanmin(sim.arrival[ids])) if len(ids) and np.isfinite(sim.arrival[ids]).any() else None
    info = dict(offsets=list(offsets), tilt=float(sim.clock.tilt), energy_drift=float(sim.energy_drift),
                t_min=float(sim.t_min), dt=float(sim.dt), wave_energy=e0,
                arrival_rim_A=first(ra), arrival_rim_B=first(rb), arrival_threshold=float(sim.threshold))
    return h, info


def displacements(J: np.ndarray):
    """Ambient displacement directions of A and B from angular impulses ``J[..., m, :]``."""
    a = np.cross(J[..., 0, :], AXIS)
    b = np.cross(J[..., 1, :], -AXIS)
    return a, b


def follow(a: np.ndarray, b: np.ndarray) -> float:
    """b . a / |a|^2: positive when B moves the same way as A, so the mouths close in together."""
    return float(b @ a / (a @ a))
