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

Units: R = 1 and wave speed 1. Nothing here is fitted, and no quantum rule is
inserted. Mouth clock offsets (the Morris–Thorne–Yurtsever ingredient of
``transaction.network``) are not modelled: the surface is static.
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

    @property
    def n_vertices(self) -> int:
        return len(self.points)

    def curvature(self) -> np.ndarray:
        """Gaussian curvature per vertex: +1 on the sphere, −b²/r⁴ in the necks."""
        K = np.ones(self.n_vertices)
        m = self.tube_of >= 0
        K[m] = self.neck.gaussian_curvature(self.s[m])
        return K

    def operators(self) -> Tuple[sp.csr_matrix, np.ndarray]:
        """Intrinsic cotangent stiffness K (positive semidefinite) and lumped mass."""
        t = self.tris
        l = self.lengths
        l0, l1, l2 = l[:, 0], l[:, 1], l[:, 2]
        sp_ = 0.5 * (l0 + l1 + l2)
        area = np.sqrt(np.maximum(sp_ * (sp_ - l0) * (sp_ - l1) * (sp_ - l2), 1e-300))
        cot = np.stack([(l1 ** 2 + l2 ** 2 - l0 ** 2), (l0 ** 2 + l2 ** 2 - l1 ** 2),
                        (l0 ** 2 + l1 ** 2 - l2 ** 2)], axis=1) / (4 * area[:, None])
        rows, cols, vals = [], [], []
        for k, (i, j) in enumerate([(1, 2), (0, 2), (0, 1)]):   # edge opposite vertex k
            w = 0.5 * cot[:, k]
            a, b = t[:, i], t[:, j]
            rows += [a, b, a, b]
            cols += [b, a, a, b]
            vals += [-w, -w, w, w]
        n = self.n_vertices
        K = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n))
        mass = np.zeros(n)
        for k in range(3):
            np.add.at(mass, t[:, k], area / 3)
        return K, mass

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
        sphere_tris=np.arange(len(stri)), rims=rims, drawn=drawn, n_ring=n_ring)
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
# Waves
# ─────────────────────────────────────────────────────────────────────────────

class HandleWaveSim:
    """Leapfrog ``u_tt = Δu`` on the intrinsic glued surface, frames precomputed.

    The source is a compact, zero-mean bump around a sphere point (see ``initial``).
    ``arrival`` holds the first time ``|u|`` exceeds ``threshold`` times the launch
    amplitude. The default is low (0.5%) because a narrow neck reflects most of
    an incoming wave: what crosses it is a few percent of the launch amplitude.

    At that threshold the recorded arrival respects the causal bound (geodesic
    distance minus the source width) on handled surfaces. The one exception is
    the exact antipodal caustic of a bare sphere, where the coherently focused
    discretisation tail arrives up to about 0.09 early at the default mesh; the
    error halves when the vertex count is quadrupled.
    """

    def __init__(self, surf: HandledSphere, source=(-0.55, 0.62, 0.56), width: float = 0.25,
                 t_max: float = 2 * math.pi + 0.6, n_frames: int = 300, threshold: float = 0.005):
        self.surf = surf
        p = np.asarray(source, float)
        self.source = p / np.linalg.norm(p)
        self.width = float(width)
        self.t_max = float(t_max)
        self.n_frames = int(n_frames)
        self.threshold = float(threshold)
        self.K, self.mass = surf.operators()
        self._inv_mass = 1.0 / self.mass
        self.frames: Optional[np.ndarray] = None
        self.times = np.linspace(0.0, self.t_max, self.n_frames)
        self.arrival: Optional[np.ndarray] = None

    def initial(self) -> np.ndarray:
        """Compact, monopole-free launch data: (1 − x²)⁴ (1 − c x²), x = d/w < 1.

        ``c`` is set on the mesh so the mass-weighted mean is exactly zero. The
        field is then zero outside the light cone, with no constant pedestal
        (a closed surface can never shed a monopole, and subtracting it would
        put an acausal offset everywhere).
        """
        x = self.surf.points
        d = np.full(len(x), np.inf)
        m = self.surf.on_sphere
        d[m] = np.arccos(np.clip(x[m] @ self.source, -1, 1))
        q = np.where(d < self.width, (d / self.width) ** 2, 1.0)
        p = np.where(d < self.width, (1 - q) ** 4, 0.0)
        c = float(self.mass @ p) / float(self.mass @ (p * q))
        return p * (1 - c * q)

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
        """Leapfrog's exactly conserved energy for the step u_prev -> u_next."""
        v = (u_next - u_prev) / dt
        return 0.5 * float(v @ (self.mass * v) + u_prev @ (self.K @ u_next))

    def run(self) -> "HandleWaveSim":
        frame_dt = self.t_max / (self.n_frames - 1)
        dt0 = 1.8 / self.max_frequency()
        sub = max(1, int(math.ceil(frame_dt / dt0)))
        dt = frame_dt / sub
        self.dt = dt
        u0 = self.initial()
        amp = float(np.max(np.abs(u0)))
        acc = lambda u: -self._inv_mass * (self.K @ u)
        prev = u0
        cur = u0 + 0.5 * dt * dt * acc(u0)
        frames = np.empty((self.n_frames, len(u0)), np.float32)
        frames[0] = u0
        arrival = np.full(len(u0), np.nan)
        arrival[np.abs(u0) > self.threshold * amp] = 0.0
        e0 = self.energy(prev, cur, dt)
        t = dt
        step = 1
        for f in range(1, self.n_frames):
            for _ in range(sub):
                if step > 1:
                    nxt = 2 * cur - prev + dt * dt * acc(cur)
                    prev, cur = cur, nxt
                hit = np.isnan(arrival) & (np.abs(cur) > self.threshold * amp)
                arrival[hit] = t
                t += dt
                step += 1
            frames[f] = cur
        self.energy_drift = abs(self.energy(prev, cur, dt) - e0) / e0
        self.frames, self.arrival, self.amplitude = frames, arrival, amp
        return self

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
