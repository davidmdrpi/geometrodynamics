"""Scalar waves on a sphere carrying antipodal finite-mouth handles.

THE OBJECT
──────────
A unit sphere with one or two antipodal mouth pairs of angular radius ``a``.
Each pair is joined by the finite-mouth neck of ``bulk.finite_mouth``::

    r(s) = sqrt(s² + b²),    b = R sin²a,    |s| <= S = R sin a cos a

which meets the sphere with matching radius and slope at both rims (a C¹
join; the Gaussian curvature jumps from +1/R² to −1/R² there). A handle is
glued either orientation-preservingly (**torus-type**) or orientation-reversingly
(**Klein-type**). One Klein-type handle makes a Klein bottle (χ = 0); two make
the connected sum of two Klein bottles (χ = −2). The antipodal map of a mouth
``S^k`` has degree ``(−1)^(k+1)``; for the S² mouths of a 3D space it reverses
orientation, which is what the Klein-type gluing stands in for in this 2D
picture.

WHAT IS COMPUTED ON WHAT
───────────────────────
The scalar wave equation ``u_tt = Δu`` runs on the **intrinsic** glued surface.
Sphere triangles use their embedded (true) edge lengths. Neck triangles use the
neck metric ``ds² + r(s)² dψ²``, so a wave crosses a neck in its true length
``2S``. Rim vertices are shared, and the far rim is attached through the
gluing permutation, so the field is continuous across both mouths.

The **drawn** handle is only where values are displayed. No surface in 3D can
show the neck at true length (``2S < 2R cos a``, the straight chord), so its
thin middle is stretched along an arc. Each Klein-type handle is immersed: it
passes through the sphere wall and joins the far mouth from inside. That
crossing is a drawing artifact; the abstract surface does not intersect itself.
A scalar field needs no normal, so it is well defined on the non-orientable
surface; a "height along the normal" would not be, which is why the field is
shown as colour and isolines.

MOUTH CLOCK OFFSETS (Morris–Thorne–Yurtsever)
─────────────────────────────────────────────
Each handle may carry a clock offset ``Δ``, the ``Δ_BA = o_B − o_A`` of
``transaction.network``: a signal crossing the neck from mouth A to mouth B
emerges ``Δ`` later in external (Killing) time, and one crossing B → A emerges
``Δ`` earlier. The jump sits at the neck waist. The spacetime is then stationary,
not static, and has no global ``t = const`` Cauchy surface. Waves are evolved in
a tilted time ``τ = t − θ`` (``ClockField``), which exists exactly when no closed
causal curve does. The two-handle case shows the Roman-ring effect: offsets
that are each safe alone can close a causal loop through both handles
(``loop_lengths``, ``chronology_ratio``).

Units: R = 1 and wave speed 1. Nothing here is fitted, and no quantum rule is
inserted: this is the classical wave equation on a classical spacetime.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import scipy.sparse as sp
from scipy.spatial import ConvexHull

__all__ = [
    "HandleSpec",
    "Neck",
    "HandledSphere",
    "HandleWaveSim",
    "ClockField",
    "ChronologyError",
    "CutMesh",
    "clock_field",
    "source_clock",
    "LoopLength",
    "loop_lengths",
    "chronology_ratio",
    "neck_crossing_length",
    "default_handles",
    "build_surface",
    "is_orientable",
    "euler_characteristic",
    "orientation_loop",
    "symlog_curvature",
]


# ─────────────────────────────────────────────────────────────────────────────
# Specifications
# ─────────────────────────────────────────────────────────────────────────────

@dataclass(frozen=True)
class HandleSpec:
    """One antipodal mouth pair: mouth A at ``+axis``, mouth B at ``−axis``.

    ``side`` is the direction the drawn outside arc bulges toward; it only
    affects the picture. ``klein`` selects the orientation-reversing gluing.
    """

    axis: Tuple[float, float, float]
    side: Tuple[float, float, float]
    klein: bool = True


def default_handles(n: int = 2, klein: Sequence[bool] | bool = True) -> Tuple[HandleSpec, ...]:
    """Up to two non-colliding handles: mouths at ±z (arc on +y) and ±x (arc on −y)."""
    if n not in (0, 1, 2):
        raise ValueError("this layout supports 0, 1 or 2 handles")
    kinds = [klein] * n if isinstance(klein, bool) else list(klein)
    if len(kinds) != n:
        raise ValueError("one klein flag per handle")
    base = [((0.0, 0.0, 1.0), (0.0, 1.0, 0.0)), ((1.0, 0.0, 0.0), (0.0, -1.0, 0.0))]
    return tuple(HandleSpec(base[i][0], base[i][1], bool(kinds[i])) for i in range(n))


@dataclass(frozen=True)
class Neck:
    """The finite-mouth neck for mouth angular radius ``a`` (R = 1)."""

    a: float

    @property
    def sa(self) -> float:
        return math.sin(self.a)

    @property
    def ca(self) -> float:
        return math.cos(self.a)

    @property
    def b(self) -> float:
        return self.sa ** 2

    @property
    def S(self) -> float:
        return self.sa * self.ca

    @property
    def s1(self) -> float:
        """Where the exact flare hands over to the stretched middle (drawing only)."""
        return min(math.sqrt(3.0) * self.b, 0.6 * self.S)

    def r(self, s):
        return np.sqrt(np.asarray(s, float) ** 2 + self.b ** 2)

    def gaussian_curvature(self, s):
        """Gaussian curvature of the 2D slice, −b²/r⁴."""
        s = np.asarray(s, float)
        return -(self.b ** 2) / (s * s + self.b ** 2) ** 2

    def flare_height(self, s):
        """Exact catenoid rise above the rim plane for |s| in [s1, S]."""
        s = np.abs(np.asarray(s, float))
        return self.b * (np.arcsinh(self.S / self.b) - np.arcsinh(s / self.b))

    @property
    def flare_rise(self) -> float:
        return float(self.flare_height(self.s1))

    @property
    def area(self) -> float:
        """Exact neck area ∫ 2π r ds over [−S, S]."""
        S, b = self.S, self.b
        return 2 * math.pi * (S * math.sqrt(S * S + b * b) + b * b * math.asinh(S / b))


# ─────────────────────────────────────────────────────────────────────────────
# Drawn handle (display only)
# ─────────────────────────────────────────────────────────────────────────────

def _catmull_rom(points: np.ndarray, n: int = 4000) -> np.ndarray:
    """Centripetal Catmull–Rom through 2D control points, densely sampled."""
    P = np.asarray(points, float)
    ext = np.vstack([2 * P[0] - P[1], P, 2 * P[-1] - P[-2]])
    out = []
    per = max(8, n // (len(P) - 1))
    for i in range(1, len(ext) - 2):
        p0, p1, p2, p3 = ext[i - 1], ext[i], ext[i + 1], ext[i + 2]
        t0 = 0.0
        t1 = t0 + np.linalg.norm(p1 - p0) ** 0.5
        t2 = t1 + np.linalg.norm(p2 - p1) ** 0.5
        t3 = t2 + np.linalg.norm(p3 - p2) ** 0.5
        ts = np.linspace(t1, t2, per, endpoint=(i == len(ext) - 3))[:, None]
        a1 = (t1 - ts) / (t1 - t0) * p0 + (ts - t0) / (t1 - t0) * p1
        a2 = (t2 - ts) / (t2 - t1) * p1 + (ts - t1) / (t2 - t1) * p2
        a3 = (t3 - ts) / (t3 - t2) * p2 + (ts - t2) / (t3 - t2) * p3
        b1 = (t2 - ts) / (t2 - t0) * a1 + (ts - t0) / (t2 - t0) * a2
        b2 = (t3 - ts) / (t3 - t1) * a2 + (ts - t1) / (t3 - t1) * a3
        out.append((t2 - ts) / (t2 - t1) * b1 + (ts - t1) / (t2 - t1) * b2)
    return np.vstack(out)


def _rmf(P: np.ndarray, T: np.ndarray, N0: np.ndarray) -> np.ndarray:
    """Rotation-minimizing frames by double reflection."""
    N = np.empty_like(P)
    N[0] = N0 - T[0] * (T[0] @ N0)
    N[0] /= np.linalg.norm(N[0])
    for i in range(len(P) - 1):
        v1 = P[i + 1] - P[i]
        c1 = v1 @ v1
        if c1 < 1e-16:
            N[i + 1] = N[i]
            continue
        rL = N[i] - (2 / c1) * (v1 @ N[i]) * v1
        tL = T[i] - (2 / c1) * (v1 @ T[i]) * v1
        v2 = T[i + 1] - tL
        c2 = v2 @ v2
        n = rL if c2 < 1e-16 else rL - (2 / c2) * (v2 @ rL) * v2
        n = n - T[i + 1] * (T[i + 1] @ n)
        N[i + 1] = n / np.linalg.norm(n)
    return N


@dataclass
class _DrawnHandle:
    """Display embedding of one handle, in the plane spanned by its axis and side."""

    spec: HandleSpec
    neck: Neck
    c: np.ndarray = field(init=False)
    e1: np.ndarray = field(init=False)
    e2: np.ndarray = field(init=False)

    def __post_init__(self):
        c = np.asarray(self.spec.axis, float)
        c /= np.linalg.norm(c)
        e1 = np.asarray(self.spec.side, float)
        e1 = e1 - c * (c @ e1)
        e1 /= np.linalg.norm(e1)
        self.c, self.e1, self.e2 = c, e1, np.cross(c, e1)
        nk = self.neck
        zt = nk.ca + nk.flare_rise
        if self.spec.klein:
            self.rho_t = (nk.sa - float(nk.r(nk.s1))) / (1 + nk.sa)
            self.z_end = -nk.ca * (1 - self.rho_t)
            self.punt_len = self.rho_t * (math.pi / 2 + nk.a)
            ctrl = [(0, zt), (0, zt + .25), (.38, zt + .62), (1.22, 1.02), (1.42, .38), (.72, 0.0),
                    (.28, self.z_end + .6), (0, self.z_end + .28), (0, self.z_end)]
        else:
            ctrl = [(0, zt), (0, zt + .25), (.38, zt + .62), (1.25, 1.08), (1.68, 0.0), (1.25, -1.08),
                    (.38, -(zt + .62)), (0, -(zt + .25)), (0, -zt)]
        vu = _catmull_rom(np.array(ctrl, float))
        P = vu[:, 1:2] * c + vu[:, 0:1] * e1          # u along c, v along e1
        seg = np.linalg.norm(np.diff(P, axis=0), axis=1)
        arc = np.concatenate([[0.0], np.cumsum(seg)])
        T = np.gradient(P, axis=0)
        T /= np.linalg.norm(T, axis=1)[:, None]
        T[0] = c
        T[-1] = -c if self.spec.klein else c
        N = _rmf(P, T, e1)
        self.length = float(arc[-1])
        self._arc, self._P, self._T, self._N = arc / arc[-1], P, T, N
        Ne = N[-1] - c * (c @ N[-1])
        Ne /= np.linalg.norm(Ne)
        Te = T[-1]
        Be = np.cross(Te, Ne)
        self._Nend, self._Bend = Ne, Be
        q0 = Ne
        q1 = Be
        alpha0 = math.atan2(q0 @ self.e2, q0 @ self.e1)
        sign = 1.0 if (math.atan2(q1 @ self.e2, q1 @ self.e1) - alpha0) % (2 * math.pi) < math.pi else -1.0
        if self.spec.klein:
            if sign > 0:
                raise RuntimeError("drawn Klein handle does not reverse the rim orientation")
            self.twist = alpha0 - math.pi            # lands column psi on rim-B angle pi − psi
        else:
            if sign < 0:
                raise RuntimeError("drawn torus handle reverses the rim orientation")
            self.twist = math.pi - alpha0            # lands column psi on rim-B angle psi + pi

    def target_angle(self, psi):
        """Rim-B angle glued to tube column psi (in the (e1, e2) frame)."""
        return (math.pi - psi) if self.spec.klein else (psi + math.pi)

    def _middle(self, f):
        i = np.interp(f, self._arc, np.arange(len(self._arc)))
        k = np.clip(np.floor(i).astype(int), 0, len(self._arc) - 2)
        w = (i - k)[..., None]
        P = (1 - w) * self._P[k] + w * self._P[k + 1]
        N = (1 - w) * self._N[k] + w * self._N[k + 1]
        T = (1 - w) * self._T[k] + w * self._T[k + 1]
        N = N - T * np.sum(N * T, axis=-1, keepdims=True)
        N /= np.linalg.norm(N, axis=-1, keepdims=True)
        B = np.cross(T / np.linalg.norm(T, axis=-1, keepdims=True), N)
        return P, N, B

    def point(self, s, psi):
        """Drawn position of neck point (s, psi); psi is the intrinsic column angle."""
        nk = self.neck
        s = np.asarray(s, float)
        psi = np.asarray(psi, float)
        s, psi = np.broadcast_arrays(s, psi)
        out = np.empty(s.shape + (3,))
        r = nk.r(s)
        A = s >= nk.s1
        if np.any(A):
            h = nk.ca + nk.flare_height(s[A])
            ang = psi[A]
            out[A] = h[:, None] * self.c + r[A][:, None] * (np.cos(ang)[:, None] * self.e1 + np.sin(ang)[:, None] * self.e2)
        M = (s < nk.s1) & (s > -nk.s1)
        if np.any(M):
            f = (nk.s1 - s[M]) / (2 * nk.s1)
            P, N, B = self._middle(f)
            ang = psi[M] + self.twist * f
            out[M] = P + r[M][:, None] * (np.cos(ang)[:, None] * N + np.sin(ang)[:, None] * B)
        Bm = s <= -nk.s1
        if np.any(Bm):
            if self.spec.klein:
                sig = self.punt_len * (nk.S + s[Bm]) / (nk.S - nk.s1)
                beta = math.pi + nk.a - sig / self.rho_t
                rho = nk.sa - self.rho_t * (np.sin(beta) + nk.sa)
                z = -nk.ca + self.rho_t * (np.cos(beta) + nk.ca)
                ang = self.target_angle(psi[Bm])
                out[Bm] = z[:, None] * self.c + rho[:, None] * (np.cos(ang)[:, None] * self.e1 + np.sin(ang)[:, None] * self.e2)
            else:
                h = -(nk.ca + nk.flare_height(s[Bm]))
                ang = psi[Bm] + self.twist
                out[Bm] = h[:, None] * self.c + r[Bm][:, None] * (np.cos(ang)[:, None] * self._Nend + np.sin(ang)[:, None] * self._Bend)
        return out


# ─────────────────────────────────────────────────────────────────────────────
# The glued surface
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class CutMesh:
    """Display mesh cut at the neck waists.

    ``orig``    surface vertex behind each display vertex (originals come first)
    ``handle``  −1 for originals; for a waist copy, the handle it belongs to
    ``tris``    triangles in display indices
    ``jump``    per triangle corner: the handle whose clock jump applies there, else −1
    """

    orig: np.ndarray
    handle: np.ndarray
    tris: np.ndarray
    jump: np.ndarray


@dataclass
class HandledSphere:
    """Abstract glued surface plus its drawing.

    ``points``     drawn positions (sphere vertices are exact)
    ``tris``       triangles (vertex indices; orientation not meaningful)
    ``lengths``    intrinsic edge lengths per triangle, opposite each vertex
    ``on_sphere``  vertex lies on the sphere band (rims included)
    ``tube_of``    handle index for neck vertices, −1 elsewhere
    ``s``          neck coordinate for neck vertices (nan on the sphere interior)
    """

    a: float
    handles: Tuple[HandleSpec, ...]
    neck: Neck
    points: np.ndarray
    tris: np.ndarray
    lengths: np.ndarray
    on_sphere: np.ndarray
    tube_of: np.ndarray
    s: np.ndarray
    sphere_tris: np.ndarray
    rims: List[Tuple[np.ndarray, np.ndarray]]
    drawn: List[_DrawnHandle]
    n_ring: int
    tube_ids: List[np.ndarray] = field(default_factory=list)
    tri_tube: Optional[np.ndarray] = None
    tri_layer: Optional[np.ndarray] = None
    n_s: int = 0

    @property
    def waist(self) -> int:
        """Ring index of the neck waist (s = 0), where each handle's clock jump sits."""
        return self.n_s // 2

    @property
    def n_vertices(self) -> int:
        return len(self.points)

    def curvature(self) -> np.ndarray:
        """Gaussian curvature per vertex: +1 on the sphere, −b²/r⁴ in the necks."""
        K = np.ones(self.n_vertices)
        m = self.tube_of >= 0
        K[m] = self.neck.gaussian_curvature(self.s[m])
        return K

    def local_stiffness(self) -> Tuple[np.ndarray, np.ndarray]:
        """Per-triangle area and 3×3 cotangent matrices, θᵀ K_T θ = ∫_T |∇θ|²."""
        l = self.lengths
        l0, l1, l2 = l[:, 0], l[:, 1], l[:, 2]
        sp_ = 0.5 * (l0 + l1 + l2)
        area = np.sqrt(np.maximum(sp_ * (sp_ - l0) * (sp_ - l1) * (sp_ - l2), 1e-300))
        cot = np.stack([(l1 ** 2 + l2 ** 2 - l0 ** 2), (l0 ** 2 + l2 ** 2 - l1 ** 2),
                        (l0 ** 2 + l1 ** 2 - l2 ** 2)], axis=1) / (4 * area[:, None])
        KT = np.zeros((len(l), 3, 3))
        for k, (i, j) in enumerate([(1, 2), (0, 2), (0, 1)]):   # edge opposite vertex k
            w = 0.5 * cot[:, k]
            KT[:, i, j] -= w
            KT[:, j, i] -= w
            KT[:, i, i] += w
            KT[:, j, j] += w
        return area, KT

    def assemble(self, KT: np.ndarray) -> sp.csr_matrix:
        t = self.tris
        rows = np.repeat(t, 3, axis=1).ravel()
        cols = np.tile(t, (1, 3)).ravel()
        n = self.n_vertices
        return sp.csr_matrix((KT.ravel(), (rows, cols)), shape=(n, n))

    def operators(self) -> Tuple[sp.csr_matrix, np.ndarray]:
        """Intrinsic cotangent stiffness K (positive semidefinite) and lumped mass."""
        area, KT = self.local_stiffness()
        mass = np.zeros(self.n_vertices)
        for k in range(3):
            np.add.at(mass, self.tris[:, k], area / 3)
        return self.assemble(KT), mass

    def cut(self) -> "CutMesh":
        """The surface cut open at every neck waist, for display in external time.

        Each waist ring is duplicated; triangles on the mouth-B side of the waist use
        the copies. A field shown in external (Killing) time jumps there by the
        handle's clock offset, and the cut lets the picture show that jump.
        """
        n = self.n_vertices
        orig = [np.arange(n)]
        handle = [np.full(n, -1)]
        tris = self.tris.copy()
        jump = np.full(self.tris.shape, -1)
        nxt = n
        for k, ids in enumerate(self.tube_ids):
            ring = ids[self.waist]
            dup = np.full(n, -1)
            dup[ring] = nxt + np.arange(len(ring))
            orig.append(ring)
            handle.append(np.full(len(ring), k))
            nxt += len(ring)
            bside = np.where((self.tri_tube == k) & (self.tri_layer >= self.waist))[0]
            sub = tris[bside]
            on = dup[sub] >= 0
            sub[on] = dup[sub][on]
            tris[bside] = sub
            jump[bside] = np.where(on, k, -1)
        return CutMesh(orig=np.concatenate(orig), handle=np.concatenate(handle), tris=tris, jump=jump)

    def total_area(self) -> float:
        return float(self.operators()[1].sum())


def _fibonacci_sphere(n: int) -> np.ndarray:
    i = np.arange(n) + 0.5
    z = 1 - 2 * i / n
    phi = math.pi * (3 - math.sqrt(5)) * i
    rr = np.sqrt(1 - z * z)
    return np.stack([rr * np.cos(phi), rr * np.sin(phi), z], axis=1)


def build_surface(a: float = 0.3, handles: Optional[Sequence[HandleSpec]] = None,
                  n_ring: int = 48, n_sphere: int = 6000, ds_target: Optional[float] = None) -> HandledSphere:
    """Build the glued surface with exact rims, intrinsic neck metric and drawing."""
    handles = tuple(default_handles(2) if handles is None else handles)
    if n_ring % 2:
        raise ValueError("n_ring must be even so both gluings map rim vertices to rim vertices")
    if not 0.05 < a < math.pi / 4:
        raise ValueError("mouth angle must lie in (0.05, pi/4)")
    neck = Neck(a)
    centers = []
    for h in handles:
        c = np.asarray(h.axis, float)
        c /= np.linalg.norm(c)
        centers += [c, -c]
    for i in range(len(centers)):
        for j in range(i + 1, len(centers)):
            if math.acos(np.clip(centers[i] @ centers[j], -1, 1)) < 2 * a + 0.05:
                raise ValueError("mouths overlap; reduce a or move the handles apart")

    drawn = [_DrawnHandle(h, neck) for h in handles]
    psi = 2 * math.pi * np.arange(n_ring) / n_ring
    rim_pts, rims = [], []
    idx = 0
    for d in drawn:
        ring = []
        for sign in (+1, -1):
            pts = sign * neck.ca * d.c + neck.sa * (np.cos(psi)[:, None] * d.e1 + np.sin(psi)[:, None] * d.e2)
            rim_pts.append(pts)
            ring.append(np.arange(idx, idx + n_ring))
            idx += n_ring
        rims.append((ring[0], ring[1]))
    n_rim = idx
    spacing = math.sqrt(4 * math.pi / n_sphere)
    fib = _fibonacci_sphere(n_sphere)
    keep = np.ones(len(fib), bool)
    for c in centers:
        keep &= np.arccos(np.clip(fib @ c, -1, 1)) > a + 0.7 * spacing
    sphere_pts = np.vstack(rim_pts + [fib[keep]]) if rim_pts else fib
    hull = ConvexHull(sphere_pts)
    if n_rim and not np.all(np.isin(np.arange(n_rim), hull.vertices)):
        raise RuntimeError("a rim vertex was dropped by the hull")
    stri = hull.simplices.copy()
    rim_id = np.full(len(sphere_pts), -1)
    for k, (ra, rb) in enumerate(rims):
        rim_id[ra] = 2 * k
        rim_id[rb] = 2 * k + 1
    if n_rim:
        same = (rim_id[stri[:, 0]] >= 0) & (rim_id[stri[:, 0]] == rim_id[stri[:, 1]]) & (rim_id[stri[:, 1]] == rim_id[stri[:, 2]])
        stri = stri[~same]
    pts = [sphere_pts]
    on_sphere = [np.ones(len(sphere_pts), bool)]
    tube_of = [np.full(len(sphere_pts), -1)]
    svals = [np.full(len(sphere_pts), np.nan)]
    for k, (ra, rb) in enumerate(rims):
        svals[0][ra] = neck.S
        svals[0][rb] = -neck.S
    sphere_len = np.linalg.norm(sphere_pts[stri[:, [1, 2, 0]]] - sphere_pts[stri[:, [2, 0, 1]]], axis=2)
    tris = [stri]
    lengths = [sphere_len]
    nv = len(sphere_pts)
    ds = ds_target if ds_target is not None else min(0.02, neck.S / 14)
    n_s = max(16, int(math.ceil(2 * neck.S / ds)))
    n_s += n_s % 2                       # even, so a ring sits exactly on the waist s = 0
    tube_ids, tri_tube, tri_layer = [], [np.full(len(stri), -1)], [np.full(len(stri), -1)]
    s_grid = neck.S - 2 * neck.S * np.arange(n_s + 1) / n_s
    half = n_ring // 2
    for k, (d, (ra, rb)) in enumerate(zip(drawn, rims)):
        perm = (half - np.arange(n_ring)) % n_ring if d.spec.klein else (np.arange(n_ring) + half) % n_ring
        inner = s_grid[1:-1]
        S2, P2 = np.meshgrid(inner, psi, indexing="ij")
        drawn_pts = d.point(S2.ravel(), P2.ravel())
        ids = np.empty((n_s + 1, n_ring), int)
        ids[0] = ra
        ids[1:-1] = nv + np.arange(len(inner) * n_ring).reshape(len(inner), n_ring)
        ids[-1] = rb[perm]
        tube_ids.append(ids)
        pts.append(drawn_pts)
        on_sphere.append(np.zeros(len(drawn_pts), bool))
        tube_of.append(np.full(len(drawn_pts), k))
        svals.append(S2.ravel())
        nv += len(drawn_pts)
        for i in range(n_s):
            j = np.arange(n_ring)
            jn = (j + 1) % n_ring
            v00, v01, v10, v11 = ids[i, j], ids[i, jn], ids[i + 1, j], ids[i + 1, jn]
            sa_, sb_ = s_grid[i], s_grid[i + 1]
            dpsi = 2 * math.pi / n_ring
            # each quad -> (v00, v10, v11) and (v00, v11, v01), coordinates (s, psi_local)
            for tri, coords in (((v00, v10, v11), ((sa_, 0.0), (sb_, 0.0), (sb_, dpsi))),
                                ((v00, v11, v01), ((sa_, 0.0), (sb_, dpsi), (sa_, dpsi)))):
                tris.append(np.stack(tri, axis=1))
                tri_tube.append(np.full(n_ring, k))
                tri_layer.append(np.full(n_ring, i))
                L = np.empty((n_ring, 3))
                for m, (p, q) in enumerate([(1, 2), (2, 0), (0, 1)]):
                    (s_p, y_p), (s_q, y_q) = coords[p], coords[q]
                    rbar = float(neck.r(0.5 * (s_p + s_q)))
                    L[:, m] = math.sqrt((s_p - s_q) ** 2 + (2 * rbar * math.sin(0.5 * abs(y_p - y_q))) ** 2)
                lengths.append(L)
    surf = HandledSphere(
        a=a, handles=handles, neck=neck,
        points=np.vstack(pts), tris=np.vstack(tris), lengths=np.vstack(lengths),
        on_sphere=np.concatenate(on_sphere), tube_of=np.concatenate(tube_of), s=np.concatenate(svals),
        sphere_tris=np.arange(len(stri)), rims=rims, drawn=drawn, n_ring=n_ring,
        tube_ids=tube_ids, tri_tube=np.concatenate(tri_tube), tri_layer=np.concatenate(tri_layer), n_s=n_s)
    return surf


# ─────────────────────────────────────────────────────────────────────────────
# Topology
# ─────────────────────────────────────────────────────────────────────────────

def _edge_map(tris: np.ndarray) -> Dict[Tuple[int, int], List[Tuple[int, int]]]:
    em: Dict[Tuple[int, int], List[Tuple[int, int]]] = {}
    for t, (i, j, k) in enumerate(tris):
        for a, b in ((i, j), (j, k), (k, i)):
            key = (a, b) if a < b else (b, a)
            em.setdefault(key, []).append((t, 1 if a < b else -1))
    return em


def euler_characteristic(tris: np.ndarray) -> int:
    nv = len(np.unique(tris))
    return nv - len(_edge_map(tris)) + len(tris)


def is_orientable(tris: np.ndarray) -> bool:
    """Combinatorial test: can all triangles be oriented consistently?

    Raises if the surface is not closed (an edge without exactly two faces).
    """
    em = _edge_map(tris)
    adj: List[List[Tuple[int, int]]] = [[] for _ in range(len(tris))]
    for key, uses in em.items():
        if len(uses) != 2:
            raise ValueError(f"edge {key} has {len(uses)} faces; surface is not closed")
        (t1, d1), (t2, d2) = uses
        # consistent orientation needs the shared edge traversed in opposite directions
        rel = -d1 * d2            # +1: same orientation sign required, −1: opposite
        adj[t1].append((t2, rel))
        adj[t2].append((t1, rel))
    sign = np.zeros(len(tris), int)
    for start in range(len(tris)):
        if sign[start]:
            continue
        sign[start] = 1
        stack = [start]
        while stack:
            t = stack.pop()
            for u, rel in adj[t]:
                want = sign[t] * rel
                if sign[u] == 0:
                    sign[u] = want
                    stack.append(u)
                elif sign[u] != want:
                    return False
    return True


# ─────────────────────────────────────────────────────────────────────────────
# Mouth clock offsets
# ─────────────────────────────────────────────────────────────────────────────

class ChronologyError(ValueError):
    """The clock offsets admit (or come too close to) closed causal curves."""


@dataclass
class ClockField:
    """A time function for a surface whose handles carry clock offsets.

    Handle ``k`` identifies its two waist sides with a time shift: a signal that
    crosses from the mouth-A side to the mouth-B side arrives ``offsets[k]`` later
    in external (Killing) time, the ``Δ_BA = o_B − o_A`` of ``transaction.network``.
    Negative offsets are advances.

    ``theta`` is single-valued on the surface cut at the waists and jumps by
    ``offsets[k]`` across waist ``k``. External time is ``t = τ + theta``, and
    the spacetime metric is the stationary

        ds² = −(dτ + dθ)² + g ,

    with ``dθ`` the smooth part. The slices ``τ = const`` are spacelike exactly
    when ``|dθ| < 1`` everywhere, so they are Cauchy surfaces then. ``theta`` is
    chosen to make the largest slope ``tilt = max |dθ|`` as small as possible.
    By duality that minimum is the largest ratio of clock gain to length over
    all closed loops, so ``tilt < 1`` exactly when no closed causal curve exists.
    The computed tilt overestimates the exact minimum by a few percent (finite
    p and mesh), so it errs on the safe side.
    """

    offsets: np.ndarray
    theta: np.ndarray
    corner: np.ndarray
    slope: np.ndarray
    display: np.ndarray
    tube_tilt: np.ndarray

    @property
    def tilt(self) -> float:
        return float(self.slope.max()) if self.slope.size else 0.0

    @property
    def chronology_safe(self) -> bool:
        return self.tilt < 1.0

    def binding(self, frac: float = 0.97) -> str:
        """Which closed loops set the tilt: through one handle or through both."""
        if not np.any(self.offsets):
            return "none"
        hot = [k for k, t in enumerate(self.tube_tilt) if t >= frac * self.tilt]
        if len(hot) >= 2:
            return "both handles (Roman ring)"
        if hot:
            return f"handle {hot[0] + 1}"
        return "sphere"


def _slopes(surf: HandledSphere, corner: np.ndarray, area: np.ndarray, KT: np.ndarray) -> np.ndarray:
    return np.sqrt(np.maximum(np.einsum("ti,tij,tj->t", corner, KT, corner), 0.0) / area)


def clock_field(surf: HandledSphere, offsets: Sequence[float], anchor: int = 0,
                fixed: Optional[np.ndarray] = None, schedule: Sequence[int] = (4, 8, 16, 32, 64, 128),
                maxiter: int = 800, gauge: Optional[np.ndarray] = None,
                cut: Optional[CutMesh] = None) -> ClockField:
    """Minimal-tilt clock function for the given waist offsets (see ``ClockField``).

    Starts from the harmonic solution and minimises the p-mean of ``|dθ|`` with
    L-BFGS for increasing p, which tends to the minimax. ``theta = 0`` on the
    vertex mask ``fixed`` (default: ``anchor`` alone). Holding it at zero over a
    source disk makes the slice flat there, so the source emits at one external
    time. ``gauge`` adds a single-valued function (an exact form). That changes
    the slicing but not the physics, and is used to test exactly that.
    """
    from scipy.optimize import minimize
    import scipy.sparse.linalg as spl

    cut = surf.cut() if cut is None else cut
    off = np.zeros(len(surf.tube_ids))
    off[:len(offsets)] = np.asarray(offsets, float)[:len(off)]
    area, KT = surf.local_stiffness()
    tris, n = surf.tris, surf.n_vertices
    off_ = np.append(off, 0.0)                     # index −1 (no jump) reads the trailing 0
    J = off_[cut.jump]
    theta = np.zeros(n)
    if fixed is None:
        fixed = np.arange(n) == anchor
    if np.any(off):
        keep = ~np.asarray(fixed, bool)
        rhs = np.zeros(n)
        KJ = np.einsum("tij,tj->ti", KT, J)
        for k in range(3):
            np.add.at(rhs, tris[:, k], -KJ[:, k])
        K = surf.assemble(KT).tocsr()
        theta[keep] = spl.spsolve(K[keep][:, keep].tocsc(), rhs[keep])
        A = area.sum()
        flat = tris.ravel()
        kt = [[np.ascontiguousarray(KT[:, i, j]) for j in range(3)] for i in range(3)]
        for p in schedule:
            def f(x, p=p):
                theta[keep] = x
                c = theta[tris] + J
                Kc = np.stack([kt[i][0] * c[:, 0] + kt[i][1] * c[:, 1] + kt[i][2] * c[:, 2] for i in range(3)], axis=1)
                g = np.sqrt(np.maximum((c * Kc).sum(axis=1), 0.0) / area)
                gm = g.max()
                r = g / gm
                F = gm * ((area * r ** p).sum() / A) ** (1.0 / p)
                coef = (F / gm) ** (1 - p) * r ** (p - 2) / (gm * A)
                gc = coef[:, None] * Kc
                grad = np.bincount(flat, weights=gc.ravel(), minlength=n)
                return F, grad[keep]
            res = minimize(f, theta[keep], jac=True, method="L-BFGS-B",
                           options=dict(maxiter=maxiter, gtol=1e-12, ftol=1e-14))
            theta[keep] = res.x
    if gauge is not None:
        theta = theta + np.asarray(gauge, float)
    theta = theta - theta[anchor]
    corner = theta[tris] + J
    slope = _slopes(surf, corner, area, KT)
    display = theta[cut.orig] + off_[cut.handle]
    tube_tilt = np.array([slope[surf.tri_tube == k].max() for k in range(len(surf.tube_ids))])
    return ClockField(offsets=off, theta=theta, corner=corner, slope=slope, display=display,
                      tube_tilt=tube_tilt)


def neck_crossing_length(neck: Neck, n: int = 300):
    """Shortest neck crossing from rim A to rim B with angular displacement ``Δψ``.

    Geodesics of ``ds² + r(s)² dψ²`` obey Clairaut's ``r² ψ' = c``. Those with
    ``c < b`` cross the waist after turning by

        Δψ(c) = 2c ∫₀^S ds / (r √(r² − c²)),   length(c) = 2 ∫₀^S r ds / √(r² − c²),

    and ``Δψ → ∞`` as ``c → b``, so every displacement is reachable. Returns a
    vectorised ``N(Δψ)`` (``Δψ`` taken mod 2π) with ``N(0) = 2S``.
    """
    from scipy.integrate import quad

    b, S = neck.b, neck.S
    cs = np.concatenate([[0.0], b * (1 - np.geomspace(1e-9, 1, n))[::-1]])
    cs = np.unique(cs[cs < b])
    turn, length = [], []
    for c in cs:
        k2 = b * b - c * c
        turn.append(2 * quad(lambda s: c / (math.sqrt(s * s + b * b) * math.sqrt(s * s + k2)), 0, S, limit=200)[0])
        length.append(2 * quad(lambda s: math.sqrt(s * s + b * b) / math.sqrt(s * s + k2), 0, S, limit=200)[0])
    turn, length = np.array(turn), np.array(length)

    def N(dpsi):
        d = np.abs((np.asarray(dpsi, float) + math.pi) % (2 * math.pi) - math.pi)
        return np.interp(d, turn, length)
    return N


@dataclass(frozen=True)
class LoopLength:
    """Shortest closed loop in one class: through handle ``i`` (A→B), then optionally ``j``.

    ``sign`` is +1 when handle ``j`` is also crossed A→B (the loop gains
    ``Δ_i + Δ_j``), −1 when it is crossed B→A (gains ``Δ_i − Δ_j``), and 0 for a
    loop through handle ``i`` alone. ``obstructed`` flags that a sphere leg of
    the minimiser runs through another mouth. The length is then only a lower
    bound for the class.
    """

    i: int
    j: Optional[int]
    sign: int
    length: float
    obstructed: bool

    def gain(self, offsets: Sequence[float]) -> float:
        d = list(offsets) + [0.0] * 4
        return d[self.i] + (self.sign * d[self.j] if self.j is not None else 0.0)

    @property
    def label(self) -> str:
        if self.j is None:
            return f"handle {self.i + 1}"
        return f"ring {self.i + 1}{'+' if self.sign > 0 else '−'}{self.j + 1}"


def loop_lengths(surf: HandledSphere) -> List[LoopLength]:
    """Exact shortest loops through one handle and through each pair (both senses).

    A loop is built from neck crossings (``neck_crossing_length``) and
    great-circle legs between rim points, glued as in ``build_surface``. Its
    length is minimised over the rim angles. A closed causal curve exists in a
    class exactly when its clock gain is at least its length, so with
    ``ratio = max |gain| / length`` the offsets are chronology-safe for these
    classes when ``ratio < 1``. The minimal-tilt clock field gives the same
    number from the other side, by duality, and covers every class.
    """
    from scipy.optimize import minimize

    nk = surf.neck
    N = neck_crossing_length(nk)
    ca, sa = nk.ca, nk.sa
    fr = []
    for d in surf.drawn:
        beta = (lambda p: math.pi - p) if d.spec.klein else (lambda p: p + math.pi)
        fr.append((d.c, d.e1, d.e2, beta))
    centers = [s * f[0] for f in fr for s in (1, -1)]

    def P(f, psi, side):
        psi = np.asarray(psi, float)[..., None]
        return side * ca * f[0] + sa * (np.cos(psi) * f[1] + np.sin(psi) * f[2])

    def arc(p, q):
        return np.arccos(np.clip(np.sum(p * q, axis=-1), -1, 1))

    def blocked(p, q, own):
        t = np.linspace(0, 1, 60)[:, None]
        w = arc(p, q)
        if w < 1e-9:
            return False
        path = (np.sin((1 - t) * w) * p + np.sin(t * w) * q) / math.sin(w)
        return any(arc(path, m).min() < nk.a - 1e-6 for k, m in enumerate(centers) if k not in own)

    def solve(fn, legs, dims, own):
        g = np.linspace(0, 2 * math.pi, 25 if dims == 4 else 73)[:-1]
        grid = np.meshgrid(*([g] * dims), indexing="ij")
        vals = fn(*grid)
        best = None
        for k in np.argsort(vals.ravel())[:12 if dims == 4 else 4]:
            x0 = [G.ravel()[k] for G in grid]
            r = minimize(lambda v: float(fn(*v)), x0, method="Nelder-Mead",
                         options=dict(xatol=1e-10, fatol=1e-12, maxiter=4000))
            if best is None or r.fun < best.fun:
                best = r
        return best.fun, any(blocked(p, q, own) for p, q in legs(*best.x))

    out = []
    for i, f in enumerate(fr):
        fn = lambda a, b: N(b - a) + arc(P(f, f[3](b), -1), P(f, a, 1))
        legs = lambda a, b: [(P(f, f[3](b), -1), P(f, a, 1))]
        L, obs = solve(fn, legs, 2, {2 * i, 2 * i + 1})
        out.append(LoopLength(i, None, 0, L, obs))
    for i in range(len(fr)):
        for j in range(i + 1, len(fr)):
            fi, fj = fr[i], fr[j]
            own = {2 * i, 2 * i + 1, 2 * j, 2 * j + 1}
            plus = lambda a, b, p, q: (N(b - a) + arc(P(fi, fi[3](b), -1), P(fj, p, 1))
                                       + N(q - p) + arc(P(fj, fj[3](q), -1), P(fi, a, 1)))
            plus_legs = lambda a, b, p, q: [(P(fi, fi[3](b), -1), P(fj, p, 1)), (P(fj, fj[3](q), -1), P(fi, a, 1))]
            minus = lambda a, b, p, q: (N(b - a) + arc(P(fi, fi[3](b), -1), P(fj, fj[3](q), -1))
                                        + N(q - p) + arc(P(fj, p, 1), P(fi, a, 1)))
            minus_legs = lambda a, b, p, q: [(P(fi, fi[3](b), -1), P(fj, fj[3](q), -1)), (P(fj, p, 1), P(fi, a, 1))]
            for sign, fn, legs in ((1, plus, plus_legs), (-1, minus, minus_legs)):
                L, obs = solve(fn, legs, 4, own)
                out.append(LoopLength(i, j, sign, L, obs))
    return out


def chronology_ratio(loops: Sequence[LoopLength], offsets: Sequence[float]) -> Tuple[float, Optional[LoopLength]]:
    """``max |clock gain| / length`` over the loops, and the loop that attains it."""
    best, arg = 0.0, None
    for lp in loops:
        r = abs(lp.gain(offsets)) / lp.length
        if r > best:
            best, arg = r, lp
    return best, arg


def source_clock(surf: HandledSphere, source, offsets: Sequence[float], width: float = 0.25,
                 margin: float = 0.05, gauge: Optional[np.ndarray] = None,
                 cut: Optional[CutMesh] = None) -> ClockField:
    """Clock field held at ``θ = 0`` on the disk of radius ``width + margin`` around ``source``.

    With the slice flat over the source, a bump laid on it is an emission at a
    single external time. That is what makes runs with different offsets
    comparable.
    """
    p = np.asarray(source, float)
    p = p / np.linalg.norm(p)
    m = np.where(surf.on_sphere)[0]
    dist = np.full(surf.n_vertices, np.inf)
    dist[m] = np.arccos(np.clip(surf.points[m] @ p, -1, 1))
    return clock_field(surf, offsets, anchor=int(m[np.argmin(dist[m])]), fixed=dist < width + margin,
                       gauge=gauge, cut=cut)


# ─────────────────────────────────────────────────────────────────────────────
# Waves
# ─────────────────────────────────────────────────────────────────────────────

class HandleWaveSim:
    """Scalar waves ``□u = 0`` on the glued surface, with optional mouth clock offsets.

    Without offsets this is leapfrog for ``u_tt = Δu``. With offsets the field
    is evolved in the clock function's time ``τ`` (see ``ClockField``), where the
    wave equation of ``−(dτ + dθ)² + g`` reads

        (1 − |dθ|²) u_ττ + 2 dθ·∇u_τ + (Δθ) u_τ − Δu = 0 .

    Its finite-element form ``M ü + C u̇ + K u = 0`` has a skew gyroscopic ``C``
    that does no work. The centred scheme conserves ``½v·Mv + ½uₙ·K uₙ₊₁``
    exactly, and ``energy_drift`` reports it.

    The source is a compact bump placed on the slice ``τ = 0`` around a sphere
    point, with ``θ = 0`` there. It is corrected (see ``initial``) so that it
    excites neither the constant mode nor its linearly growing partner. Before
    ``τ = 0`` the field is zero, which makes the source a retarded emission.

    ``frames`` and ``arrival`` are given on the ``CutMesh`` display vertices
    (surface vertices first) in **external time** ``t = τ + θ``. Each neck
    therefore shows its clock jump at the waist. ``arrival`` holds the first
    external time ``|u|`` exceeds ``threshold`` times the launch amplitude. The
    default is low (0.5%) because a narrow neck reflects most of an incoming
    wave: what crosses it is a few percent of the launch amplitude.

    At that threshold the recorded arrival respects the causal bound (geodesic
    distance minus the source width) on handled surfaces. The one exception is
    the exact antipodal caustic of a bare sphere, where the coherently focused
    discretisation tail arrives up to about 0.09 early at the default mesh; the
    error halves when the vertex count is quadrupled.
    """

    max_tilt = 0.95

    def __init__(self, surf: HandledSphere, source=(-0.55, 0.62, 0.56), width: float = 0.25,
                 t_max: float = 2 * math.pi + 0.6, n_frames: int = 300, threshold: float = 0.005,
                 offsets: Optional[Sequence[float]] = None, gauge: Optional[np.ndarray] = None,
                 clock: Optional[ClockField] = None):
        self.surf = surf
        p = np.asarray(source, float)
        self.source = p / np.linalg.norm(p)
        self.width = float(width)
        self.t_max = float(t_max)
        self.n_frames = int(n_frames)
        self.threshold = float(threshold)
        self.cut = surf.cut()
        m = np.where(surf.on_sphere)[0]
        self.source_vertex = int(m[np.argmax(surf.points[m] @ self.source)])
        if clock is None:
            clock = source_clock(surf, self.source, offsets if offsets is not None else (),
                                 width=self.width, gauge=gauge, cut=self.cut)
        self.clock = clock
        if clock.tilt >= 1.0:
            raise ChronologyError(
                f"closed causal curves: max |dθ| = {clock.tilt:.3f} >= 1 (loops through {clock.binding()})")
        if clock.tilt > self.max_tilt:
            raise ChronologyError(
                f"within {100 * (1 - clock.tilt):.0f}% of the chronology horizon (max |dθ| = "
                f"{clock.tilt:.3f}); the evolution would need a much finer mesh and step")
        self.K, self.mass0 = surf.operators()
        area = self.mass0.sum()
        self._area_T, KT = surf.local_stiffness()
        w = self._area_T / 3 * (1 - clock.slope ** 2)
        self.mass = np.zeros(surf.n_vertices)
        for k in range(3):
            np.add.at(self.mass, surf.tris[:, k], w)
        self.C = None
        if np.any(clock.corner):
            g = np.einsum("tij,tj->ti", KT, clock.corner) / 3          # B_ij = Σ_T (K_T θ_T)_j / 3, i ∈ T
            rows = np.repeat(surf.tris, 3, axis=1).ravel()
            cols = np.tile(surf.tris, (1, 3)).ravel()
            vals = np.repeat(g[:, None, :], 3, axis=1).ravel()
            B = sp.csr_matrix((vals, (rows, cols)), shape=(surf.n_vertices,) * 2)
            self.C = (B - B.T).tocsr()
        self._inv_mass = 1.0 / self.mass
        self._area_total = area
        self.frames: Optional[np.ndarray] = None
        self.arrival: Optional[np.ndarray] = None
        lo = min(0.0, float(clock.display.min()))
        self.t_min = math.floor(lo / 0.05) * 0.05 if lo < 0 else 0.0
        self.times = np.linspace(self.t_min, self.t_max, self.n_frames)

    def initial(self) -> np.ndarray:
        """Compact launch data ``(1 − x²)⁴ (1 + c₁x² + c₂x⁴)``, ``x = d/w < 1``.

        Without offsets only ``c₁`` is used, to make the mass-weighted mean zero.
        That removes the constant pedestal a closed surface would otherwise carry,
        and with it an acausal offset everywhere. With offsets the zero frequency
        has a Jordan pair (the constant ``1`` and ``τ·1 + φ₁`` with ``Kφ₁ = −C1``).
        Both conserved charges, ``Q = 1·C u`` and ``W = 1·M u − φ₁·C u``, are then
        set to zero.
        """
        x = self.surf.points
        d = np.full(len(x), np.inf)
        m = self.surf.on_sphere
        d[m] = np.arccos(np.clip(x[m] @ self.source, -1, 1))
        q = np.where(d < self.width, (d / self.width) ** 2, 1.0)
        p = np.where(d < self.width, (1 - q) ** 4, 0.0)
        basis = [p, p * q, p * q * q]
        if self.C is None:
            return p - (self.mass @ p) / (self.mass @ (p * q)) * p * q
        import scipy.sparse.linalg as spl
        n = self.surf.n_vertices
        keep = np.arange(n) != self.source_vertex
        phi = np.zeros(n)
        rhs = -(self.C @ np.ones(n))
        phi[keep] = spl.spsolve(self.K[keep][:, keep].tocsc(), rhs[keep])
        ones = np.ones(n)
        Q = lambda u: ones @ (self.C @ u)
        W = lambda u: self.mass @ u - phi @ (self.C @ u)
        rows = [[Q(b) for b in basis], [W(b) for b in basis]]
        scale = max(abs(v) for v in rows[1])
        rows = [r for r in rows if max(abs(v) for v in r) > 1e-12 * scale]   # Q ≡ 0 if dθ vanishes at the source
        A = np.array([r[1:] for r in rows])
        c = np.linalg.lstsq(A, -np.array([r[0] for r in rows]), rcond=None)[0]
        return basis[0] + c[0] * basis[1] + c[1] * basis[2]

    def max_frequency(self, iters: int = 60) -> float:
        rng = np.random.default_rng(0)
        v = rng.standard_normal(self.surf.n_vertices)
        lam = 0.0
        for _ in range(iters):
            w = self._inv_mass * (self.K @ v)
            lam = float(np.linalg.norm(w) / np.linalg.norm(v))
            v = w / np.linalg.norm(w)
        return math.sqrt(1.05 * lam)

    def energy(self, u_prev, u_next, dt) -> float:
        """The centred scheme's exactly conserved energy for the step u_prev -> u_next."""
        v = (u_next - u_prev) / dt
        return 0.5 * float(v @ (self.mass * v) + u_prev @ (self.K @ u_next))

    def run(self) -> "HandleWaveSim":
        frame_dt = (self.t_max - self.t_min) / (self.n_frames - 1)
        theta_min = float(self.clock.display.min())
        tau_end = self.t_max - min(0.0, theta_min) + 2 * frame_dt
        n_tau = int(math.ceil(tau_end / frame_dt)) + 1
        dt0 = 1.8 / self.max_frequency()
        sub = max(1, int(math.ceil(frame_dt / dt0)))
        dt = frame_dt / sub
        self.dt = dt
        u0 = self.initial()
        amp = float(np.max(np.abs(u0)))
        M = self.mass
        if self.C is None:
            step = lambda cur, prev: 2 * cur - prev - dt * dt * self._inv_mass * (self.K @ cur)
        else:
            import scipy.sparse.linalg as spl
            Ap = (sp.diags(M / dt ** 2) + self.C / (2 * dt)).tocsc()
            Am = (sp.diags(M / dt ** 2) - self.C / (2 * dt)).tocsr()
            lu = spl.splu(Ap)
            step = lambda cur, prev: lu.solve(2 * M / dt ** 2 * cur - self.K @ cur - Am @ prev)
        prev = u0
        cur = u0 - 0.5 * dt * dt * self._inv_mass * (self.K @ u0)
        tau_frames = np.empty((n_tau, len(u0)), np.float32)
        tau_frames[0] = u0
        arrival = np.full(len(u0), np.nan)
        arrival[np.abs(u0) > self.threshold * amp] = 0.0
        e0 = self.energy(prev, cur, dt)
        t = dt
        first = True
        for f in range(1, n_tau):
            for _ in range(sub):
                if not first:
                    prev, cur = cur, step(cur, prev)
                first = False
                hit = np.isnan(arrival) & (np.abs(cur) > self.threshold * amp)
                arrival[hit] = t
                t += dt
            tau_frames[f] = cur
        self.energy_drift = abs(self.energy(prev, cur, dt) - e0) / e0
        self.tau_frames, self.tau_arrival, self.amplitude = tau_frames, arrival, amp
        self.frame_dt = frame_dt
        self.frames = self._external_frames()
        self.arrival = arrival[self.cut.orig] + self.clock.display
        return self

    def _external_frames(self) -> np.ndarray:
        """Resample τ-frames at external times ``t = τ + θ`` on the display vertices."""
        orig, th = self.cut.orig, self.clock.display
        out = np.zeros((self.n_frames, len(orig)), np.float32)
        last = len(self.tau_frames) - 1
        for i, t in enumerate(self.times):
            x = (t - th) / self.frame_dt
            ok = x >= 0
            i0 = np.minimum(np.floor(x).astype(int), last - 1)
            fr = np.clip(x - i0, 0, 1)
            i0 = np.clip(i0, 0, last - 1)
            val = (1 - fr) * self.tau_frames[i0, orig] + fr * self.tau_frames[i0 + 1, orig]
            out[i] = np.where(ok, val, 0.0)
        return out

    def frame_at(self, t: float) -> int:
        return int(np.clip(round((t - self.t_min) / (self.t_max - self.t_min) * (self.n_frames - 1)),
                           0, self.n_frames - 1))

    def compressed(self, frame: int, eps: float = 0.005) -> np.ndarray:
        """Sign-preserving log compression of u, in [−1, 1], so weak transmitted fronts stay visible."""
        u = self.frames[frame] / self.amplitude
        return np.sign(u) * np.log1p(np.abs(u) / eps) / math.log1p(1 / eps)

    def value_at(self, point, frame: int) -> float:
        x = self.surf.points
        p = np.asarray(point, float)
        i = int(np.argmin(np.linalg.norm(x - p, axis=1)))
        return float(self.frames[frame, i])


# ─────────────────────────────────────────────────────────────────────────────
# Orientation loop on the drawing
# ─────────────────────────────────────────────────────────────────────────────

def orientation_loop(surf: HandledSphere, handle: int = 0, n: int = 120):
    """Carry the surface normal once around a loop through one handle.

    Returns ``(points, normals, flipped)``. The loop leaves a latitude 0.6 rad
    above the handle's equator, enters mouth A, runs down the neck column
    psi = 0, leaves mouth B, and returns over the sphere. ``flipped`` is True
    when the normal comes back reversed, which happens exactly for a
    Klein-type handle.
    """
    d = surf.drawn[handle]
    nk = surf.neck
    pts, nrm = [], []

    def push(x, nvec):
        nvec = nvec / np.linalg.norm(nvec)
        if nrm and nvec @ nrm[-1] < 0:
            nvec = -nvec
        pts.append(x)
        nrm.append(nvec)

    th_r = math.pi / 2 - 0.6

    def sph(th, ph):
        return math.cos(th) * d.c + math.sin(th) * (math.cos(ph) * d.e1 + math.sin(ph) * d.e2)

    for th in np.linspace(th_r, nk.a, n):
        x = sph(th, 0.0)
        push(x, x)
    eps = 1e-5
    for s in np.linspace(nk.S, -nk.S, 4 * n):
        x = d.point(s, 0.0)
        s_lo, s_hi = max(-nk.S, s - eps), min(nk.S, s + eps)
        dsv = d.point(s_hi, 0.0) - d.point(s_lo, 0.0)
        dpv = d.point(s, eps) - d.point(s, -eps)
        push(x, np.cross(dsv, dpv))
    phB = d.target_angle(0.0)
    for th in np.linspace(math.pi - nk.a, th_r, n):
        x = sph(th, phB)
        push(x, x)
    back = ((-phB + math.pi) % (2 * math.pi)) - math.pi
    for ph in np.linspace(phB, phB + back, n):
        x = sph(th_r, ph)
        push(x, x)
    pts, nrm = np.array(pts), np.array(nrm)
    return pts, nrm, bool(nrm[-1] @ nrm[0] < 0)


def symlog_curvature(surf: HandledSphere) -> np.ndarray:
    """Curvature compressed to [−1, 1]: +1/R² on the sphere maps to a fixed blue, the waist to −1."""
    K = surf.curvature()
    kmax = 1.0 / surf.neck.b ** 2
    return np.sign(K) * np.sqrt(np.log1p(np.abs(K)) / math.log1p(kmax))
