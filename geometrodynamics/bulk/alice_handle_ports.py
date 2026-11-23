"""What crosses an Alice handle, booked at its two mouths.

Freeze: ``docs/alice_handle_ports_prereg.md`` (commit ``36ba9c2``).

One handle joins antipodal mouths of the unit sphere on the z-axis. A real
scalar obeys ``u_tt = Δu`` on the glued surface. Separating the azimuth,

    u = Σ_m [ A_m(ℓ, t) cos mφ + B_m(ℓ, t) sin mφ ],

each ``(A_m, B_m)`` lives on a closed loop ``ℓ``: the sphere band from rim A
(``θ = a``) to rim B (``θ = π − a``), then the neck from rim B back to rim A.
The circle radius ``ρ`` is ``sin θ`` on the band and ``r(s) = sqrt(s² + b²)``
on the neck. The gluing enters only at rim B, as a 2×2 map ``T`` from the
neck's frame to the band's:

* torus-type, twist ``c``: ``T`` is the rotation by ``m c``;
* Klein-type: ``T = diag(1, −1)``.

Discretization: linear finite elements on the loop, lumped (diagonal) mass,
consistent stiffness and ``m²/ρ²`` potential by Gauss quadrature. Every node
stores ``(A, B)`` in its own frame: band nodes in ``φ``, neck nodes and rim B
in ``ψ``, rim A in both (they agree there). With lumped mass the angular
momentum and energy of any set of nodes change only through the edges that
leave it, so the semi-discrete continuity equations are exact; see
``edge_flux_L`` and ``edge_flux_E``.

Time evolution is a sum over eigenmodes, so there is no time-stepping error.
Eigenvalues are evaluated as Rayleigh quotients in energy form (sums of
nonnegative terms) of the computed eigenvectors; the raw eigensolver values
are kept as well. That choice was made before any run.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property

import numpy as np
from scipy import linalg

GAUSS = np.polynomial.legendre.leggauss(4)
J2 = np.array([[0.0, -1.0], [1.0, 0.0]])


def rotation(angle: float) -> np.ndarray:
    """Band-frame (A, B) from neck-frame (A, B) when ψ = φ + c, angle = m c."""
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, s], [-s, c]])


@dataclass(frozen=True)
class Gluing:
    kind: str                 # "torus" or "klein"
    twist: float = 0.0        # c, torus only

    def matrix(self, m: int) -> np.ndarray:
        if m == 0:
            return np.ones((1, 1))
        if self.kind == "torus":
            return rotation(m * self.twist)
        if self.kind == "klein":
            return np.diag([1.0, -1.0])
        raise ValueError(self.kind)

    @property
    def label(self) -> str:
        return "klein" if self.kind == "klein" else f"torus c={self.twist:.6g}"


TORUS = Gluing("torus", 0.0)
KLEIN = Gluing("klein")


class Loop:
    """The band-plus-neck loop for mouth angle ``a`` with ``n`` nodes."""

    def __init__(self, a: float, n: int):
        self.a, self.n = float(a), int(n)
        self.b = np.sin(a) ** 2
        self.S = np.sin(a) * np.cos(a)
        self.band_length = np.pi - 2 * a
        self.length = self.band_length + 2 * self.S
        self.n_band = int(round(n * self.band_length / self.length))
        self.n_neck = n - self.n_band
        self.h_band = self.band_length / self.n_band
        self.h_neck = 2 * self.S / self.n_neck
        self.rim_a, self.rim_b = 0, self.n_band
        # node coordinate: ℓ on the band (0 at rim A), s on the neck.
        self.ell = np.full(n, np.nan)
        self.ell[: self.n_band + 1] = self.h_band * np.arange(self.n_band + 1)
        self.s = np.full(n, np.nan)
        self.s[self.n_band:] = self.S - self.h_neck * np.arange(self.n_neck)
        self.s[0] = -self.S
        self.on_band = np.zeros(n, bool)
        self.on_band[: self.n_band] = True          # rim B node belongs to the neck
        # elements (i, j, length, rho at the Gauss points)
        x, w = GAUSS
        u = 0.5 * (x + 1)                            # local coordinate in [0, 1]
        elems = []
        for i in range(self.n_band):
            l0 = i * self.h_band
            rho = np.sin(a + l0 + u * self.h_band)
            elems.append((i, i + 1, self.h_band, rho, True))
        for k in range(self.n_neck):
            i = self.n_band + k
            j = (i + 1) % n
            s0 = self.S - k * self.h_neck
            rho = np.sqrt((s0 - u * self.h_neck) ** 2 + self.b ** 2)
            elems.append((i, j, self.h_neck, rho, False))
        self.elements = elems
        self.weights = 0.5 * w
        self.shape = np.stack([1 - u, u])           # φ_i, φ_j at the Gauss points
        mass = np.zeros(n)
        for i, j, h, rho, _ in elems:
            mass[i] += h * np.sum(self.weights * rho * self.shape[0])
            mass[j] += h * np.sum(self.weights * rho * self.shape[1])
        self.mass = mass

    def rho_exact(self):
        """ρ at every node from the formula, for checking the C¹ join."""
        r = np.empty(self.n)
        band = np.arange(self.n_band + 1)
        r[band] = np.sin(self.a + self.ell[band])
        neck = np.arange(self.n_band + 1, self.n)
        r[neck] = np.sqrt(self.s[neck] ** 2 + self.b ** 2)
        return r


class Sector:
    """One azimuthal number ``m`` with one gluing: operators, modes, fluxes."""

    def __init__(self, loop: Loop, m: int, gluing: Gluing):
        self.loop, self.m, self.gluing = loop, int(m), gluing
        self.d = 1 if m == 0 else 2
        self.T = gluing.matrix(m)
        n, d = loop.n, self.d
        # couplings: (i, j, w, Tij) with (K x)_i += w Tij x_j, frames included
        stiff = np.zeros((n * d, n * d))
        self.edges = {}
        self.potential = []
        for idx, (i, j, h, rho, band) in enumerate(loop.elements):
            kint = h * np.sum(loop.weights * rho) / h ** 2
            pot = (m ** 2) * h * np.einsum("q,aq,bq->ab", loop.weights / rho,
                                          loop.shape, loop.shape)
            local = kint * np.array([[1.0, -1.0], [-1.0, 1.0]]) + pot
            Ti = np.eye(d)
            Tj = self.T if (band and j == loop.rim_b) else np.eye(d)
            Ts = (Ti, Tj)
            for p, gp in enumerate((i, j)):
                for q, gq in enumerate((i, j)):
                    block = local[p, q] * (Ts[p].T @ Ts[q])
                    stiff[gp * d:(gp + 1) * d, gq * d:(gq + 1) * d] += block
            self.edges[(i, j)] = (local[0, 1], Ti.T @ Tj)    # maps j's frame to i's
            self.edges[(j, i)] = (local[1, 0], Tj.T @ Ti)
            self.potential.append((i, j, h, rho, Ti, Tj, kint))
        self.K = stiff
        self.M = np.repeat(loop.mass, d)

    # ---- spectrum --------------------------------------------------------

    @cached_property
    def modes(self):
        scale = 1 / np.sqrt(self.M)
        C = self.K * scale[:, None] * scale[None, :]
        raw, Y = linalg.eigh(C)
        V = Y * scale[:, None]                       # M-orthonormal
        refined = self.energy_form(V) / np.einsum("ik,i,ik->k", V, self.M, V)
        order = np.argsort(refined)
        return refined[order], raw[order], V[:, order]

    def energy_form(self, V, chunk=256):
        """vᵀKv for each column, as a sum of nonnegative terms.

        ∫ρ |w'|² + m² ∫ |w|²/ρ element by element, with the rim-B frame map
        applied. Accurate to rounding relative to the value itself, which the
        eigensolver's eigenvalues are not for the lowest modes.
        """
        loop, d = self.loop, self.d
        I = np.array([e[0] for e in self.potential])
        Jn = np.array([e[1] for e in self.potential])
        h = np.array([e[2] for e in self.potential])
        rho = np.array([e[3] for e in self.potential])           # (E, q)
        kint = np.array([e[6] for e in self.potential])          # ∫ρ / h²
        special = [k for k, e in enumerate(self.potential)
                   if not np.allclose(e[5], np.eye(d))]
        V = np.atleast_2d(V.T).T
        out = np.empty(V.shape[1])
        for lo in range(0, V.shape[1], chunk):
            X = V[:, lo:lo + chunk].reshape(loop.n, d, -1)
            xi, xj = X[I], X[Jn].copy()                          # (E, d, K)
            for k in special:
                xj[k] = self.potential[k][5] @ xj[k]
            total = np.einsum("e,edk->k", kint, (xj - xi) ** 2)
            if self.m:
                vals = (loop.shape[0][None, :, None, None] * xi[:, None]
                        + loop.shape[1][None, :, None, None] * xj[:, None])
                total += self.m ** 2 * np.einsum(
                    "e,q,eq,eqdk->k", h, loop.weights, 1 / rho, vals ** 2)
            out[lo:lo + chunk] = total
        return out

    def sector_of(self, V):
        """For d = 2: fraction of each mode's norm in the A (cosine) component."""
        x = V.reshape(self.loop.n, 2, -1)
        a = np.einsum("i,ik->k", self.loop.mass, x[:, 0, :] ** 2)
        b = np.einsum("i,ik->k", self.loop.mass, x[:, 1, :] ** 2)
        return a / (a + b)

    # ---- fluxes (semi-discrete, exact) ----------------------------------

    def angular_momentum(self, x, xdot, nodes=None):
        """π m Σ m_i (A Ḃ − Ȧ B) over ``nodes`` (all by default), own frames."""
        n = self.loop.n
        x = x.reshape(*x.shape[:-1], n, 2)
        xd = xdot.reshape(*xdot.shape[:-1], n, 2)
        weights = self.loop.mass if nodes is None else np.where(
            np.isin(np.arange(n), nodes), self.loop.mass, 0.0)
        return np.pi * self.m * np.einsum(
            "i,...i->...", weights, x[..., 0] * xd[..., 1] - xd[..., 0] * x[..., 1])

    def energy_flux(self, i, j, x, xdot):
        """Energy leaving the node set containing i across edge (i, j)."""
        w, Tij = self.edges[(i, j)]
        d = self.d
        xi, xj = x[..., i * d:(i + 1) * d], x[..., j * d:(j + 1) * d]
        vi, vj = xdot[..., i * d:(i + 1) * d], xdot[..., j * d:(j + 1) * d]
        a = np.einsum("...a,ab,...b->...", xi, Tij, vj)
        b = np.einsum("...a,ab,...b->...", vi, Tij, xj)
        return -0.5 * w * (a - b) * (np.pi if self.m else 2 * np.pi)


def lflux(sector: Sector, i: int, j: int, xi, xj):
    """Angular momentum leaving the node set that contains i, across (i, j).

    dL_R/dt = −π m Σ_{i∈R} (K x)_iᵀ J x_i, and the pairs inside R cancel, so the
    outflow across an edge to j ∉ R is  π m w_ij (T_ij x_j)ᵀ J x_i.
    """
    w, Tij = sector.edges[(i, j)]
    return np.pi * sector.m * w * np.einsum("...a,ab,...b->...",
                                            xj @ Tij.T, J2, xi)


TWO_PI = 2 * np.longdouble(np.pi)


def phase(omega, t):
    """ω t reduced modulo 2π in extended precision, returned as float64.

    At t ~ 1e8 a float64 product loses ~1e-8 of a radian; reduction in long
    double keeps the evaluated trigonometric functions at rounding level.
    """
    prod = np.multiply.outer(np.asarray(omega, np.longdouble),
                             np.asarray(t, np.longdouble))
    return np.asarray(np.mod(prod, TWO_PI), float)


def _segment_integral(Omega, tau):
    """∫₀^τ e^{iΩt} dt for an array Ω and a scalar τ, stable for all Ωτ."""
    x = 0.5 * Omega * tau
    xr = phase(0.5 * Omega, tau)
    small = np.abs(x) < 1e-3
    ratio = np.where(small, 1 - x ** 2 / 6 + x ** 4 / 120,
                     np.sin(xr) / np.where(small, 1.0, x))
    return tau * ratio * np.exp(1j * xr)


class ModeState:
    """x(t) = Re Σ_k V_k z_k e^{iω_k t}: an exact semi-discrete solution."""

    def __init__(self, sector: Sector, x0, v0, keep=None, omega=None, V=None):
        self.sector = sector
        lam, _, Vall = sector.modes
        if omega is None:
            omega = np.sqrt(np.maximum(lam, 0.0))
            V = Vall
        if np.any(omega <= 0):
            raise ValueError("ModeState needs ω > 0; m = 0 is handled separately")
        alpha = V.T @ (sector.M * x0)
        beta = V.T @ (sector.M * v0) / omega
        z = alpha - 1j * beta
        if keep is not None:
            sel = np.abs(z) * np.maximum(omega, 1) > keep * np.max(
                np.abs(z) * np.maximum(omega, 1))
            V, omega, z = V[:, sel], omega[sel], z[sel]
        self.V, self.omega, self.z = V, omega, z

    def rows(self, nodes):
        d = self.sector.d
        idx = np.concatenate([np.arange(k * d, (k + 1) * d) for k in nodes])
        return idx

    def state(self, t, rows=None):
        t = np.atleast_1d(t)
        V = self.V if rows is None else self.V[rows]
        ph = self.z[:, None] * np.exp(1j * phase(self.omega, t))
        x = (V @ ph).real.T
        xd = (V @ (1j * self.omega[:, None] * ph)).real.T
        return x, xd

    def angular_momentum(self, t, nodes=None):
        loop, m = self.sector.loop, self.sector.m
        nodes = np.arange(loop.n) if nodes is None else np.asarray(nodes)
        x, xd = self.state(t, self.rows(nodes))
        x = x.reshape(len(np.atleast_1d(t)), len(nodes), 2)
        xd = xd.reshape(len(np.atleast_1d(t)), len(nodes), 2)
        return np.pi * m * np.einsum("i,ti->t", loop.mass[nodes],
                                     x[..., 0] * xd[..., 1] - xd[..., 0] * x[..., 1])

    def lflux_integral(self, i, j, t):
        """∫₀ᵗ (L leaving i's side across (i, j)) dτ, in closed form."""
        w, Tij = self.sector.edges[(i, j)]
        G = np.pi * self.sector.m * w * (Tij.T @ J2)    # outflow = x_jᵀ G x_i
        a = self.V[self.rows([j])] * self.z[None, :]     # 2 × K for node j
        b = self.V[self.rows([i])] * self.z[None, :]     # 2 × K for node i
        P1 = a.T @ G @ b                                 # pairs: e^{i(ω+ω')}
        P2 = a.T @ G @ b.conj()                          # pairs: e^{i(ω−ω')}
        S1 = self.omega[:, None] + self.omega[None, :]
        S2 = self.omega[:, None] - self.omega[None, :]
        out = []
        for tau in np.atleast_1d(t):
            E1 = _segment_integral(S1, tau)
            E2 = _segment_integral(S2, tau)
            out.append(0.5 * np.real(np.sum(P1 * E1) + np.sum(P2 * E2)))
        return np.array(out)

    def angular_momentum_terms(self, nodes=None):
        """L(t) = Re Σ c_k e^{iΩ_k t}: the frequencies and coefficients.

        From L = ẋᵀ D x with D = π m (mass ⊗ J) and x = Re V z e^{iωt}.
        """
        sec = self.sector
        n = sec.loop.n
        weight = sec.loop.mass.copy()
        if nodes is not None:
            keep = np.zeros(n, bool)
            keep[np.asarray(nodes)] = True
            weight[~keep] = 0.0
        from scipy import sparse
        D = np.pi * sec.m * sparse.kron(sparse.diags(weight), J2, format="csr")
        Q = self.V.T @ (D @ self.V)
        wz = 1j * self.omega * self.z
        c1 = 0.5 * wz[:, None] * Q * self.z[None, :]
        c2 = 0.5 * wz[:, None] * Q * self.z.conj()[None, :]
        O1 = self.omega[:, None] + self.omega[None, :]
        O2 = self.omega[:, None] - self.omega[None, :]
        return (np.concatenate([O1.ravel(), O2.ravel()]),
                np.concatenate([c1.ravel(), c2.ravel()]))


# ---- initial data -----------------------------------------------------------

def bump(x, centre, half_width):
    """cos² bump, C¹, support |x − centre| ≤ half_width; returns (f, f')."""
    y = (x - centre) / half_width
    inside = np.abs(y) <= 1
    f = np.where(inside, np.cos(0.5 * np.pi * y) ** 2, 0.0)
    df = np.where(inside, -0.5 * np.pi / half_width * np.sin(np.pi * y), 0.0)
    return f, df


def packet(loop: Loop, d: int, lo=0.05, hi=0.45, k=14.0):
    """Band packet on ℓ ∈ [lo, hi] moving toward rim A, with L(0) > 0.

    A = e cos kℓ, B = e sin kℓ, ∂_t = ∂_ℓ (a profile of ℓ + t).
    """
    n = loop.n
    ell = np.where(loop.on_band, loop.ell, -10.0)     # off the band: outside
    e, de = bump(ell, 0.5 * (lo + hi), 0.5 * (hi - lo))
    x = np.zeros((n, d))
    v = np.zeros((n, d))
    x[:, 0] = e * np.cos(k * ell)
    v[:, 0] = de * np.cos(k * ell) - k * e * np.sin(k * ell)
    if d == 2:
        x[:, 1] = e * np.sin(k * ell)
        v[:, 1] = de * np.sin(k * ell) + k * e * np.cos(k * ell)
    return x.ravel(), v.ravel()


def generic(loop: Loop, d: int):
    """Deterministic smooth data touching band and neck, both components."""
    n = loop.n
    pos = np.where(loop.on_band, loop.ell, np.nan)
    pos[loop.rim_b:] = loop.band_length + (loop.S - loop.s[loop.rim_b:])
    pos[0] = 0.0
    f1, _ = bump(pos, 0.6, 0.35)
    f2, _ = bump(pos, loop.band_length + loop.S, 0.6 * loop.S + 0.05)
    f3, _ = bump(pos, 1.6, 0.5)
    f4, _ = bump(pos, loop.band_length - 0.2, 0.4)
    x = np.zeros((n, d))
    v = np.zeros((n, d))
    x[:, 0] = f1 + 0.4 * f2 * np.cos(9 * pos)
    v[:, 0] = 0.7 * f3 * np.sin(5 * pos)
    if d == 2:
        x[:, 1] = 0.8 * f4 + 0.3 * f2
        v[:, 1] = 1.1 * f1 * np.cos(4 * pos) - 0.5 * f2
    return x.ravel(), v.ravel()


# ---- quadratic forms in time --------------------------------------------------

def _w(state: ModeState, t):
    return state.z[:, None] * np.exp(1j * state.omega[:, None] * np.atleast_1d(t)[None, :])


class TrigSum:
    """f(t) = Re Σ c_k e^{iΩ_k t}, evaluated with reduced phases."""

    def __init__(self, Omega, coef):
        # Re(c e^{-iΩt}) = Re(c̄ e^{iΩt}): fold to Ω ≥ 0, then merge equal Ω so
        # that terms which cancel are bounded after cancelling.
        Omega = np.asarray(Omega, float)
        coef = np.asarray(coef, complex)
        neg = Omega < 0
        Omega = np.where(neg, -Omega, Omega)
        coef = np.where(neg, coef.conj(), coef)
        uniq, inv = np.unique(Omega, return_inverse=True)
        merged = np.zeros(len(uniq), complex)
        np.add.at(merged, inv, coef)
        keep = np.abs(merged) > 0
        self.Omega, self.coef = uniq[keep], merged[keep]

    def __call__(self, t, chunk=20_000):
        t = np.atleast_1d(np.asarray(t, float))
        out = np.empty(len(t))
        for lo in range(0, len(t), chunk):
            ph = phase(self.Omega, t[lo:lo + chunk])
            out[lo:lo + chunk] = np.real(self.coef @ np.exp(1j * ph))
        return out

    def split(self, cut):
        slow = np.abs(self.Omega) <= cut
        return (TrigSum(self.Omega[slow], self.coef[slow]),
                float(np.sum(np.abs(self.coef[~slow]))))


def angular_momentum_series(state: ModeState, t, nodes=None):
    return TrigSum(*state.angular_momentum_terms(nodes))(t)


class UncertifiedSearch(RuntimeError):
    """The search met an interval it could neither clear nor resolve."""


def first_sign_change(f: TrigSum, t_end, cut, fine, width=None):
    """First t in (0, t_end] where f leaves the sign of f(0); None if none.

    Every interval [t0, t1] is either cleared by a bound or split, so no
    crossing is skipped. Two bounds clear an interval on which f keeps the
    sign s0 of f(0):

    * slow/fast: f = slow + fast with |fast| ≤ B (frequencies above ``cut``)
      and |slow'| ≤ D_slow = Σ|cΩ| over the slow terms. If slow(t0) has sign
      s0 and |slow(t0)| − (t1 − t0) D_slow > B, f cannot vanish there.
    * Lipschitz: |f'| ≤ D = Σ|cΩ| over all terms. If f(t0) and f(t1) both have
      sign s0 and |f(t0)| + |f(t1)| > (t1 − t0) D, f cannot vanish there.
    * curvature: |f''| ≤ D2 = Σ|cΩ²|, so from either end s0 f(t0 ± τ) ≥
      s0 f(t0) ± s0 f'(t0) τ − D2 τ²/2, a concave bound; if it is positive at
      both ends of the interval, f cannot vanish there. This one clears every
      interval lying wholly to one side of a simple root, however shallow.

    Intervals longer than ``fine`` that the first bound does not clear are
    halved; shorter ones are sampled at 51 points and each sub-interval is
    tested with the Lipschitz bound, uncleared ones pushed back in order.
    An interval whose end has the other sign is bisected to the crossing;
    being the leftmost uncleared interval, that crossing is the first. An
    interval narrower than ``width`` (default 1e-12 (1 + t1)) that is neither
    cleared nor shows the other sign (f touches zero without a visible sign
    change) raises :class:`UncertifiedSearch` instead of being passed over.
    """
    s0 = np.sign(f(0.0)[0])
    slow, bound = f.split(cut)
    D_slow = float(np.sum(np.abs(slow.coef * slow.Omega)))
    D = float(np.sum(np.abs(f.coef * f.Omega)))
    D2 = float(np.sum(np.abs(f.coef) * f.Omega ** 2))
    df = TrigSum(f.Omega, 1j * f.Omega * f.coef)
    top = float(np.max(np.abs(slow.Omega))) if len(slow.Omega) else 0.0
    coarse = 0.01 / top if top > 0 else t_end
    stack = []
    edges = np.arange(0.0, t_end, coarse)
    for t0 in edges[::-1]:
        stack.append((t0, min(t0 + coarse, t_end)))

    def crossing(left, right):
        for _ in range(200):
            mid = 0.5 * (left + right)
            if mid in (left, right):
                break
            if np.sign(f(mid)[0]) == s0:
                left = mid
            else:
                right = mid
        return 0.5 * (left + right)

    while stack:
        t0, t1 = stack.pop()
        sv = slow(t0)[0]
        if np.sign(sv) == s0 and abs(sv) - (t1 - t0) * D_slow > bound:
            continue
        if t1 - t0 > fine:
            mid = 0.5 * (t0 + t1)
            stack.append((mid, t1))
            stack.append((t0, mid))
            continue
        tiny = width if width is not None else 1e-12 * (1.0 + abs(t1))
        if t1 - t0 <= tiny:
            v1 = f(t1)[0]
            if np.sign(v1) != s0:
                return crossing(t0, t1)
            raise UncertifiedSearch(f"f touches zero near t = {t0!r} without a "
                                    "certified sign change")
        grid = np.linspace(t0, t1, 51)
        vals, der = s0 * f(grid), s0 * df(grid)
        h = np.diff(grid)
        signed = vals > 0
        lip = vals[:-1] + vals[1:] > h * D
        left = vals[:-1] + der[:-1] * h - 0.5 * D2 * h ** 2 > 0
        right = vals[1:] - der[1:] * h - 0.5 * D2 * h ** 2 > 0
        clear = signed[:-1] & signed[1:] & (lip | left | right)
        open_ = np.flatnonzero(~clear)
        if not open_.size:
            continue
        k = open_[0]
        if not signed[k + 1] and signed[k]:
            return crossing(grid[k], grid[k + 1])
        for q in open_[::-1]:
            stack.append((grid[q], grid[q + 1]))
    return None


# ---- the predictions ------------------------------------------------------------

def neck_nodes(loop: Loop):
    return np.arange(loop.rim_b, loop.n)


def band_nodes(loop: Loop):
    return np.arange(0, loop.rim_b)


def p1_m0(a, n, times=np.linspace(0.0, 4.0, 2001)):
    """m = 0: spectra and energy crossing each rim for three gluings."""
    loop = Loop(a, n)
    x0, v0 = generic(loop, 1)
    out = {}
    for g in (TORUS, Gluing("torus", 1.0), KLEIN):
        sec = Sector(loop, 0, g)
        lam, raw, V = sec.modes
        om = np.sqrt(np.maximum(lam, 0.0))
        alpha, gamma = V.T @ (sec.M * x0), V.T @ (sec.M * v0)
        tt = times[None, :]
        safe = np.where(om > 1e-9, om, 1.0)[:, None]
        sinc = np.where(om[:, None] > 1e-9, np.sin(om[:, None] * tt) / safe, tt)
        X = V @ (alpha[:, None] * np.cos(om[:, None] * tt) + gamma[:, None] * sinc)
        Xd = V @ (-alpha[:, None] * om[:, None] * np.sin(om[:, None] * tt)
                  + gamma[:, None] * np.cos(om[:, None] * tt))
        X, Xd = X.T, Xd.T
        fa = sec.energy_flux(loop.n - 1, 0, X, Xd)            # leaving neck at A
        fb = sec.energy_flux(loop.rim_b, loop.rim_b - 1, X, Xd)  # leaving neck at B
        from scipy.integrate import simpson
        out[g.label] = dict(spectrum=lam, crossed_a=float(simpson(fa, x=times)),
                            crossed_b=float(simpson(fb, x=times)))
        del sec
    ref = out[TORUS.label]
    spec, cross = 0.0, 0.0
    for v in out.values():
        spec = max(spec, float(np.max(np.abs(v["spectrum"] - ref["spectrum"])
                                      / np.maximum(np.abs(ref["spectrum"]), 1.0))))
        for key in ("crossed_a", "crossed_b"):
            cross = max(cross, abs(v[key] - ref[key]) / max(abs(ref[key]), 1e-300))
    return dict(a=a, n=n, spectrum_relative=spec, crossing_relative=cross,
                crossed={k: [v["crossed_a"], v["crossed_b"]] for k, v in out.items()})


def p2_torus(sec: Sector, times=np.linspace(0.0, 10.0, 41)):
    lam = sec.modes[0]
    pairs = lam.reshape(-1, 2)
    degeneracy = float(np.max(np.abs(pairs[:, 0] - pairs[:, 1]) / pairs.mean(axis=1)))
    x0, v0 = generic(sec.loop, 2)
    st = ModeState(sec, x0, v0)
    L = np.array([st.angular_momentum(t)[0] for t in times])
    return dict(degeneracy=degeneracy, L0=float(L[0]),
                L_drift=float(np.max(np.abs(L - L[0])) / abs(L[0])))


def klein_sectors(sec: Sector):
    lam, raw, V = sec.modes
    frac = sec.sector_of(V)
    cos = frac > 0.5
    purity = float(max(np.max(1 - frac[cos]), np.max(frac[~cos])))
    return lam, V, cos, purity


def p3_union(klein: Sector, torus0: Sector, torus_pi: Sector):
    lam, V, cos, purity = klein_sectors(klein)
    t0 = torus0.modes[0].reshape(-1, 2).mean(axis=1)
    tp = torus_pi.modes[0].reshape(-1, 2).mean(axis=1)
    lc, ls = np.sort(lam[cos]), np.sort(lam[~cos])
    ok_counts = len(lc) == len(t0) and len(ls) == len(tp)
    rel_c = float(np.max(np.abs(lc - t0) / t0)) if ok_counts else np.inf
    rel_s = float(np.max(np.abs(ls - tp) / tp)) if ok_counts else np.inf
    # the whole Klein spectrum against the union, sorted together
    union = np.sort(np.concatenate([t0, tp]))
    rel_u = float(np.max(np.abs(np.sort(lam) - union) / union))
    return dict(purity=purity, counts=[int(cos.sum()), int((~cos).sum())],
                cos_vs_c0=rel_c, sin_vs_cpi=rel_s, union=rel_u)


def p4_p5(klein: Sector, pairs=6, samples=41):
    lam, V, cos, _ = klein_sectors(klein)
    loop = klein.loop
    ic, is_ = np.flatnonzero(cos), np.flatnonzero(~cos)
    rows = []
    for n in range(pairs):
        kp, km = ic[n], is_[n]
        wp, wm = np.sqrt(lam[kp]), np.sqrt(lam[km])
        delta, mean = abs(wp - wm), 0.5 * (wp + wm)
        row = dict(n=n, omega_cos=float(wp), omega_sin=float(wm),
                   split=float(delta / mean))
        if delta <= 0.2 * mean and delta > 0:
            Vp = V[:, [kp, km]]
            x0 = V[:, kp]
            v0 = wm * V[:, km]
            st = ModeState(klein, x0, v0, omega=np.array([wp, wm]), V=Vp)
            L = TrigSum(*st.angular_momentum_terms())
            L0 = float(L(0.0)[0])
            t_pred = np.pi / (2 * delta)
            t_found = first_sign_change(L, 1.5 * t_pred, cut=0.5 * mean,
                                        fine=0.5 / mean)
            ts = np.linspace(0.0, np.pi / delta, samples)
            Lt = L(ts)
            flux = st.lflux_integral(loop.rim_b, loop.rim_b - 1, ts)
            identity = np.abs(Lt - L0 + 2 * flux)
            # diagnostics, not gates: how exact the two eigenvectors are, and
            # how the identity's error grows with time
            resid = []
            for k_ in (kp, km):
                v_ = V[:, k_]
                r_ = klein.K @ v_ - lam[k_] * klein.M * v_
                resid.append(float(np.linalg.norm(r_) / np.linalg.norm(klein.K @ v_)))
            row.update(eigen_residual=max(resid),
                       identity_series=(identity / abs(L0)).tolist(),
                       identity_times=ts.tolist())
            row.update(L0=L0, t_predicted=float(t_pred),
                       t_found=None if t_found is None else float(t_found),
                       timing_relative=None if t_found is None
                       else float(abs(t_found - t_pred) / t_pred),
                       timing_bound=float(delta / mean),
                       identity=float(identity.max() / abs(L0)),
                       port_fraction=float(np.max(np.abs(flux)) / abs(L0)),
                       L_end=float(Lt[-1]))
        rows.append(row)
    return rows


def p6_packet(a, n, m, samples=401):
    loop = Loop(a, n)
    x0, v0 = packet(loop, 2)
    R = np.flatnonzero(loop.on_band & (loop.ell >= loop.band_length - 0.25))
    t0 = 0.05 + 2 * loop.S
    t1 = np.pi - 2 * a - 0.7
    times = np.linspace(t0, t1, samples)
    res = {}
    for g in (TORUS, KLEIN):
        sec = Sector(loop, m, g)
        full = ModeState(sec, x0, v0)
        L0 = float(full.angular_momentum(0.0)[0])
        LR = full.angular_momentum(times, R)
        entry = dict(L0=L0, LR=LR)
        if g.kind == "klein":
            st = ModeState(sec, x0, v0, keep=1e-14)
            ts = np.linspace(0.0, t1, 21)
            Lt = st.angular_momentum(ts)
            flux = st.lflux_integral(loop.rim_b, loop.rim_b - 1, ts)
            entry.update(identity=float(np.max(np.abs(Lt - Lt[0] + 2 * flux)) / abs(L0)),
                         kept_modes=int(len(st.omega)),
                         truncation=float(abs(Lt[0] - L0) / abs(L0)),
                         port_fraction=float(np.max(np.abs(flux)) / abs(L0)))
        res[g.kind] = entry
        del sec, full
    LT, LK = res["torus"]["LR"], res["klein"]["LR"]
    peak = float(np.max(np.abs(LT)))
    return dict(a=a, n=n, m=m, window=[float(t0), float(t1)], region_nodes=int(len(R)),
                L0=res["torus"]["L0"], L0_klein=res["klein"]["L0"],
                peak_fraction=peak / abs(res["torus"]["L0"]),
                antisymmetry=float(np.max(np.abs(LK + LT)) / peak),
                symmetry=float(np.max(np.abs(LK - LT)) / peak),
                identity=res["klein"]["identity"],
                identity_truncation=res["klein"]["truncation"],
                kept_modes=res["klein"]["kept_modes"],
                times=times.tolist(), L_region_torus=LT.tolist(),
                L_region_klein=LK.tolist())


def expected_split(a, m):
    return float(np.exp(-2 * m * np.arcsinh(1 / np.tan(a))))


# ---- P7: the #338 mesh ------------------------------------------------------------

def mesh_gluing(surf):
    """Measured far-rim/near-rim azimuth relation of the single handle."""
    ids = surf.tube_ids[0]
    pts = surf.points
    near = np.arctan2(pts[ids[0], 1], pts[ids[0], 0])
    far = np.arctan2(pts[ids[-1], 1], pts[ids[-1], 0])

    def spread(angles):
        z = np.exp(1j * angles)
        centre = np.angle(z.mean())
        dev = np.angle(z * np.exp(-1j * centre))
        return float(centre % (2 * np.pi)), float(np.max(np.abs(dev)))

    diff, sdiff = spread(far - near)
    tot, stot = spread(far + near)
    return dict(difference=diff, difference_spread=sdiff,
                sum=tot, sum_spread=stot)


def mesh_spectrum(klein, n_sphere, n_ring, count=20, a=0.3):
    from scipy.sparse import diags
    from scipy.sparse.linalg import eigsh
    from geometrodynamics.viz.klein_handles import HandleSpec, build_surface
    surf = build_surface(a=a, handles=(HandleSpec((0.0, 0.0, 1.0), (0.0, 1.0, 0.0),
                                                  klein),),
                         n_ring=n_ring, n_sphere=n_sphere)
    K, mass = surf.operators()
    vals = eigsh(K.tocsc(), k=count + 6, M=diags(mass).tocsc(), sigma=-0.05,
                 which="LM", return_eigenvectors=False)
    vals = np.sort(vals)
    return surf, vals[1:count + 1], float(vals[0])


def separable_spectrum(gluing: Gluing, a, n, m_max, count=20):
    loop = Loop(a, n)
    vals, label = [], []
    for m in range(m_max + 1):
        lam = Sector(loop, m, gluing).modes[0]
        lam = lam[lam < 200]
        vals.append(lam)
        label.append(np.full(len(lam), m))
    vals, label = np.concatenate(vals), np.concatenate(label)
    order = np.argsort(vals)
    vals, label = vals[order], label[order]
    return vals[1:count + 1], label[1:count + 1], float(vals[0])


def p7_mesh(a=0.3, n=2000, m_max=8, meshes=((6000, 48), (12000, 64)), count=20):
    out = {}
    for klein in (False, True):
        kind = "klein" if klein else "torus"
        rows = []
        geom = None
        reference = None
        for n_sphere, n_ring in meshes:
            surf, mesh_vals, zero = mesh_spectrum(klein, n_sphere, n_ring, count, a)
            g = mesh_gluing(surf)
            if geom is None:
                geom = g
                glu = KLEIN if klein else Gluing("torus", g["difference"])
                reference = separable_spectrum(glu, a, n, m_max, count)
            sep, labels, sep_zero = reference
            rel = np.abs(mesh_vals - sep) / sep
            rows.append(dict(n_sphere=n_sphere, n_ring=n_ring, gluing=g,
                             mesh=mesh_vals.tolist(), mesh_zero=zero,
                             max_relative=float(rel.max())))
        sep, labels, sep_zero = reference
        out[kind] = dict(geometry=geom, separable=sep.tolist(),
                         separable_m=labels.tolist(), separable_zero=sep_zero,
                         largest_m=int(labels.max()), runs=rows)
    return out
