"""A localized, responding Klein handle. Freeze ``5d15963``.

A lattice ring with two ports joined by a tube whose seam may reverse time
and may reverse the fibre of a charged field. The time-orientation cover is
two sheets with two cross-connected tubes; in cover time it is an ordinary
static network, and the sheet swap ``sigma`` is a spatial symmetry of it.

One waist coordinate per tube is solved from the same Hamiltonian as the
field. Nothing is prescribed except a static potential in the charged model.

    H = 1/2 |p|^2 + 1/2 X^T K(q, q') X + sum_i A_i rho_i
        + sum_tubes (p_q^2 / 2M + K_s q^2 / 2)
    K = graph Laplacian + m^2 + V,   V_i = U_0 exp(-2 q h_i) on tube sites
    rho_i = x_i p_{y,i} - y_i p_{x,i}

State layout: ``[X (c N), P (c N), q (2), p_q (2)]`` with ``c`` field
components (1 neutral, 2 charged) and ``N`` cover sites.

See ``docs/klein_handle_response_prereg.md`` for the frozen predictions and
``docs/klein_handle_response.md`` for the result. This is not Stage 1 of the
boundary-pairing protocol.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property

import numpy as np

SEED = 2026111002
MASS_SQ = 0.25
BARRIER = 2.0
MOUTH_MASS = 1.0
MOUTH_STIFFNESS = 4.0
STEP = 0.02
STEPS = 200                      # histories on [-4, 4]
FIXED_POINT_TOL = 1e-13
FIXED_POINT_MAX = 200
AMPLITUDES = (1 / 64, 1 / 32, 1 / 16, 1 / 8)
STRESS_AMPLITUDE = 0.5
ALPHAS = (1e-3, 2e-3, 4e-3)
SAMPLE_EVERY = 10

PRIMARY = dict(n_ext=12, port_a=0, port_b=5, n_tube=4)
REPLICATE = dict(n_ext=20, port_a=0, port_b=7, n_tube=6)

CELLS = ((1, 1), (1, -1), (-1, 1), (-1, -1))      # (s_t, s_f)


# --------------------------------------------------------------------------
# Network
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Network:
    """Two-sheeted cover of a ring with one handle."""

    n_ext: int
    port_a: int
    port_b: int
    n_tube: int

    @property
    def n(self) -> int:
        """Sites per sheet: the exterior and the tube that sheet owns."""
        return self.n_ext + self.n_tube

    @property
    def size(self) -> int:
        return 2 * self.n

    def ext(self, sheet: int, e: int) -> int:
        return sheet * self.n + e % self.n_ext

    def tube(self, sheet: int, k: int) -> int:
        return sheet * self.n + self.n_ext + k

    @cached_property
    def tube_sites(self) -> np.ndarray:
        """Shape ``(2, n_tube)``: sites of ``theta`` and of ``theta'``."""
        return np.array([[self.tube(s, k) for k in range(self.n_tube)]
                         for s in (0, 1)])

    @cached_property
    def profile(self) -> np.ndarray:
        """``h = (1 - s^2)^3`` at stations equally spaced in [-0.6, 0.6]."""
        s = np.linspace(-0.6, 0.6, self.n_tube)
        return (1 - s ** 2) ** 3

    @cached_property
    def edges(self) -> list:
        out = []
        for s in (0, 1):
            for e in range(self.n_ext):
                out.append((self.ext(s, e), self.ext(s, e + 1)))
            out.append((self.ext(s, self.port_a), self.tube(s, 0)))
            for k in range(self.n_tube - 1):
                out.append((self.tube(s, k), self.tube(s, k + 1)))
            # The far end of each tube lands on the *other* sheet.
            out.append((self.tube(s, self.n_tube - 1),
                        self.ext(1 - s, self.port_b)))
        return out

    @cached_property
    def laplacian(self) -> np.ndarray:
        lap = np.zeros((self.size, self.size))
        for i, j in self.edges:
            lap[i, i] += 1
            lap[j, j] += 1
            lap[i, j] -= 1
            lap[j, i] -= 1
        return lap

    @cached_property
    def swap(self) -> np.ndarray:
        """Permutation matrix of the sheet swap ``sigma``."""
        p = np.zeros((self.size, self.size))
        for i in range(self.size):
            p[i, (i + self.n) % self.size] = 1.0
        return p

    def quotient_operator(self, twist: int) -> np.ndarray:
        """Independent single-sheet handle network.

        ``twist = +1`` is the ordinary handle. ``twist = -1`` gives the last
        tube edge a sign: a field that changes sign once around the handle.
        Built without reference to the cover.
        """
        n = self.n
        op = np.zeros((n, n))

        def link(i, j, sign=1.0):
            op[i, i] += 1
            op[j, j] += 1
            op[i, j] -= sign
            op[j, i] -= sign

        for e in range(self.n_ext):
            link(e, (e + 1) % self.n_ext)
        link(self.port_a, self.n_ext)
        for k in range(self.n_tube - 1):
            link(self.n_ext + k, self.n_ext + k + 1)
        link(self.n_ext + self.n_tube - 1, self.port_b, float(twist))
        op += MASS_SQ * np.eye(n)
        op[np.arange(self.n_ext, n), np.arange(self.n_ext, n)] += BARRIER
        return op

    def potential_profile(self, parity: int) -> np.ndarray:
        """``w = 1`` on tube sites; ``parity = -1`` flips it on ``theta'``."""
        w = np.zeros(self.size)
        w[self.tube_sites[0]] = 1.0
        w[self.tube_sites[1]] = float(parity)
        return w


# --------------------------------------------------------------------------
# Coupled model
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Model:
    net: Network
    components: int = 1
    potential: np.ndarray | None = None
    respond: bool = True
    _a: np.ndarray = field(init=False, repr=False, compare=False,
                           default=None)

    def __post_init__(self):
        a = np.zeros(self.net.size) if self.potential is None \
            else np.asarray(self.potential, dtype=float)
        if self.components == 1 and np.any(a):
            raise ValueError("a neutral field does not couple to a potential")
        object.__setattr__(self, "_a", a)

    @property
    def field_dim(self) -> int:
        return 2 * self.components * self.net.size

    @property
    def dim(self) -> int:
        return self.field_dim + 4

    def split(self, z: np.ndarray):
        c, n = self.components, self.net.size
        z = np.atleast_2d(z)
        b = z.shape[0]
        x = z[:, :c * n].reshape(b, c, n)
        p = z[:, c * n:2 * c * n].reshape(b, c, n)
        return x, p, z[:, 2 * c * n:2 * c * n + 2], z[:, 2 * c * n + 2:]

    def barrier(self, q: np.ndarray) -> np.ndarray:
        """``V_i`` for each state of a batch, shape ``(B, N)``."""
        v = np.zeros((q.shape[0], self.net.size))
        for t in (0, 1):
            v[:, self.net.tube_sites[t]] = BARRIER * np.exp(
                -2.0 * q[:, t, None] * self.net.profile)
        return v

    def vector_field(self, z: np.ndarray) -> np.ndarray:
        x, p, q, pq = self.split(z)
        v = self.barrier(q)
        kx = x @ self.net.laplacian.T + (MASS_SQ + v[:, None, :]) * x
        dx, dp = p.copy(), -kx
        if self.components == 2:
            a = self._a
            dx[:, 0] -= a * x[:, 1]
            dx[:, 1] += a * x[:, 0]
            dp[:, 0] -= a * p[:, 1]
            dp[:, 1] += a * p[:, 0]
        dq = np.zeros_like(q)
        dpq = np.zeros_like(pq)
        if self.respond:
            dq = pq / MOUTH_MASS
            density = np.sum(x ** 2, axis=1)
            for t in (0, 1):
                sites = self.net.tube_sites[t]
                dpq[:, t] = -MOUTH_STIFFNESS * q[:, t] + np.sum(
                    self.net.profile * v[:, sites] * density[:, sites], axis=1)
        b = x.shape[0]
        return np.concatenate([dx.reshape(b, -1), dp.reshape(b, -1), dq, dpq],
                              axis=1)

    def energy(self, z: np.ndarray) -> np.ndarray:
        x, p, q, pq = self.split(z)
        v = self.barrier(q)
        kx = x @ self.net.laplacian.T + (MASS_SQ + v[:, None, :]) * x
        e = 0.5 * np.sum(p ** 2, axis=(1, 2)) + 0.5 * np.sum(x * kx,
                                                             axis=(1, 2))
        if self.components == 2:
            e += np.sum(self._a * (x[:, 0] * p[:, 1] - x[:, 1] * p[:, 0]),
                        axis=1)
        e += np.sum(pq ** 2, axis=1) / (2 * MOUTH_MASS)
        e += 0.5 * MOUTH_STIFFNESS * np.sum(q ** 2, axis=1)
        return e

    def jacobian(self, z: np.ndarray) -> np.ndarray:
        """Derivative of ``vector_field`` at one state."""
        c, n = self.components, self.net.size
        x, _, q, _ = self.split(z)
        x, q = x[0], q[0]
        v = self.barrier(q[None])[0]
        k = self.net.laplacian + np.diag(MASS_SQ + v)
        jac = np.zeros((self.dim, self.dim))
        xs = lambda comp: slice(comp * n, (comp + 1) * n)
        ps = lambda comp: slice(c * n + comp * n, c * n + (comp + 1) * n)
        qi, pi = 2 * c * n, 2 * c * n + 2
        for comp in range(c):
            jac[xs(comp), ps(comp)] = np.eye(n)
            jac[ps(comp), xs(comp)] = -k
        if c == 2:
            a = np.diag(self._a)
            jac[xs(0), xs(1)] = -a
            jac[xs(1), xs(0)] = a
            jac[ps(0), ps(1)] = -a
            jac[ps(1), ps(0)] = a
        if self.respond:
            jac[qi:qi + 2, pi:pi + 2] = np.eye(2) / MOUTH_MASS
            h = self.net.profile
            for t in (0, 1):
                sites = self.net.tube_sites[t]
                density = np.sum(x[:, sites] ** 2, axis=0)
                jac[pi + t, qi + t] = -MOUTH_STIFFNESS - 2.0 * np.sum(
                    h ** 2 * v[sites] * density)
                for comp in range(c):
                    coupling = 2.0 * h * v[sites] * x[comp, sites]
                    rows = np.arange(n)[sites] + c * n + comp * n
                    cols = np.arange(n)[sites] + comp * n
                    jac[rows, qi + t] = coupling
                    jac[pi + t, cols] = coupling
        return jac

    # ---- bilinear forms --------------------------------------------------

    @cached_property
    def symplectic(self) -> np.ndarray:
        """``Omega(u, v) = u^T W v``."""
        c, n = self.components, self.net.size
        w = np.zeros((self.dim, self.dim))
        f = c * n
        w[:f, f:2 * f] = np.eye(f)
        w[f:2 * f, :f] = -np.eye(f)
        w[2 * f:2 * f + 2, 2 * f + 2:] = np.eye(2)
        w[2 * f + 2:, 2 * f:2 * f + 2] = -np.eye(2)
        return w

    def charge_form(self, sites=None) -> np.ndarray:
        """Symmetric matrix of ``Q = sum_i rho_i`` over ``sites``."""
        if self.components != 2:
            raise ValueError("charge needs two components")
        n = self.net.size
        mask = np.zeros(n)
        mask[np.arange(n) if sites is None else sites] = 1.0
        m = np.zeros((self.dim, self.dim))
        d = np.diag(mask)
        x0, x1, p0, p1 = (slice(0, n), slice(n, 2 * n),
                          slice(2 * n, 3 * n), slice(3 * n, 4 * n))
        m[x0, p1] += d
        m[x1, p0] -= d
        return (m + m.T) / 2

    def charge(self, z: np.ndarray, sites=None) -> np.ndarray:
        z = np.atleast_2d(z)
        return np.einsum("bi,ij,bj->b", z, self.charge_form(sites), z)

    @cached_property
    def rotate_quarter(self) -> np.ndarray:
        """Multiplication of ``Psi`` by ``i``."""
        n = self.net.size
        r = np.eye(self.dim)
        eye, zero = np.eye(n), np.zeros((n, n))
        block = np.block([[zero, -eye], [eye, zero]])
        r[:2 * n, :2 * n] = block
        r[2 * n:4 * n, 2 * n:4 * n] = block
        return r


def deck(model: Model, s_t: int, s_f: int = 1, eta: int = 1) -> np.ndarray:
    """Deck map of the cover as a matrix on states.

    ``s_t = -1``: the seam reverses time. ``s_f = -1``: it reverses the
    fibre, an antilinear action on ``Psi = x + i y``.
    """
    c, n = model.components, model.net.size
    p = model.net.swap
    d = np.zeros((model.dim, model.dim))
    for comp in range(c):
        fibre = s_f if comp == 1 else 1
        d[comp * n:(comp + 1) * n, comp * n:(comp + 1) * n] = eta * fibre * p
        lo = c * n + comp * n
        d[lo:lo + n, lo:lo + n] = eta * s_t * fibre * p
    f = 2 * c * n
    flip = np.array([[0.0, 1.0], [1.0, 0.0]])
    d[f:f + 2, f:f + 2] = flip
    d[f + 2:, f + 2:] = s_t * flip
    return d


# --------------------------------------------------------------------------
# Integrator: implicit midpoint
# --------------------------------------------------------------------------

def midpoint_step(model: Model, z: np.ndarray, h: float) -> np.ndarray:
    z = np.atleast_2d(z)
    new = z + h * model.vector_field(z)
    for _ in range(FIXED_POINT_MAX):
        nxt = z + h * model.vector_field(0.5 * (z + new))
        change = float(np.abs(nxt - new).max())
        new = nxt
        if change < FIXED_POINT_TOL:
            return new
    raise RuntimeError("implicit midpoint did not converge")


def flow(model: Model, z0: np.ndarray, steps: int, h: float = STEP
         ) -> np.ndarray:
    """States at ``0, h, ..., steps h``; shape ``(steps + 1, B, dim)``."""
    z = np.atleast_2d(np.asarray(z0, dtype=float))
    out = np.empty((steps + 1,) + z.shape)
    out[0] = z
    for j in range(steps):
        z = midpoint_step(model, z, h)
        out[j + 1] = z
    return out


def tangent_flow(model: Model, z0: np.ndarray, basis: np.ndarray, steps: int,
                 h: float = STEP, sample: int = SAMPLE_EVERY) -> dict:
    """Linearized scheme along one history. Returns sampled tangent frames."""
    z = np.atleast_2d(np.asarray(z0, dtype=float))
    frame = basis.copy()
    eye = np.eye(model.dim)
    frames = {0: frame.copy()}
    for j in range(steps):
        new = midpoint_step(model, z, h)
        jac = model.jacobian(0.5 * (z + new))
        frame = np.linalg.solve(eye - 0.5 * h * jac,
                                (eye + 0.5 * h * jac) @ frame)
        z = new
        if (j + 1) % sample == 0:
            frames[j + 1] = frame.copy()
    return frames


# --------------------------------------------------------------------------
# Linear algebra helpers
# --------------------------------------------------------------------------

def null_space(a: np.ndarray, rtol: float = 1e-9) -> np.ndarray:
    _, s, vt = np.linalg.svd(a)
    scale = s[0] if s.size and s[0] > 0 else 1.0
    return vt[int(np.sum(s > rtol * scale)):].T.copy()


def rank(a: np.ndarray, rtol: float = 1e-9) -> int:
    if a.size == 0:
        return 0
    s = np.linalg.svd(a, compute_uv=False)
    return int(np.sum(s > rtol * s[0])) if s[0] > 0 else 0


def fixed_space(model: Model, s_t: int, s_f: int = 1, eta: int = 1
                ) -> np.ndarray:
    d = deck(model, s_t, s_f, eta)
    return null_space(d - np.eye(model.dim))


def field_fixed_space(model: Model, s_t: int, s_f: int = 1, eta: int = 1
                      ) -> np.ndarray:
    """Fixed space of the field block alone, embedded with zero mouths."""
    f = model.field_dim
    d = deck(model, s_t, s_f, eta)[:f, :f]
    basis = null_space(d - np.eye(f))
    return np.vstack([basis, np.zeros((4, basis.shape[1]))])


def admissible_state(model: Model, s_t: int, s_f: int, amplitude: float,
                     rng: np.random.Generator, eta: int = 1,
                     real_only: bool = False,
                     mouth: tuple = (0.0, 0.0)) -> np.ndarray:
    """Deck-fixed mirror data with field norm ``amplitude``.

    ``mouth = (q, p_q)`` sets ``q = q'`` and ``p_q = -s_t ... `` consistently
    with the deck map.
    """
    basis = field_fixed_space(model, s_t, s_f, eta)
    z = basis @ rng.standard_normal(basis.shape[1])
    if real_only:
        n = model.net.size
        z[n:2 * n] = 0.0
        z[3 * n:4 * n] = 0.0
    z *= amplitude / np.linalg.norm(z)
    f = model.field_dim
    z[f:f + 2] = mouth[0]
    z[f + 2], z[f + 3] = mouth[1], s_t * mouth[1]
    return z


def seam_defect(model: Model, z0: np.ndarray, d: np.ndarray,
                steps: int = STEPS) -> float:
    """``max_tau |z(tau) - D z(-tau)|`` relative to the field scale."""
    forward = flow(model, z0, steps, STEP)[:, 0]
    backward = flow(model, z0, steps, -STEP)[:, 0]
    scale = float(np.linalg.norm(z0))
    return float(np.abs(forward - backward @ d.T).max()) / scale


# --------------------------------------------------------------------------
# H1  Localized linear handle
# --------------------------------------------------------------------------

def localized_report(net: Network, rng: np.random.Generator) -> dict:
    model = Model(net, components=1, respond=False)
    w = model.symplectic
    out = {"size": net.size, "twists": []}
    f = model.field_dim
    for eta in (1, -1):
        d = deck(model, -1, 1, eta)
        fixed = null_space(d[:f, :f] - np.eye(f))
        wf = w[:f, :f]
        z0 = admissible_state(model, -1, 1, 1.0, rng, eta)
        out["twists"].append({
            "eta": eta,
            "involution_error": float(np.abs(d @ d - np.eye(model.dim)).max()),
            "antisymplectic_error": float(np.abs(d.T @ w @ d + w).max()),
            "field_dim": f, "fixed_dim": fixed.shape[1],
            "omega_on_fixed": float(np.abs(fixed.T @ wf @ fixed).max()),
            "frozen_seam_defect": seam_defect(model, z0, d),
        })
    # Spectra: sigma-even against the ordinary handle, odd against the twisted.
    cover = net.laplacian + MASS_SQ * np.eye(net.size)
    cover[net.tube_sites.ravel(), net.tube_sites.ravel()] += BARRIER
    p = net.swap
    out["swap_commutator"] = float(np.abs(p @ cover - cover @ p).max())
    spectra = {}
    for name, sign in (("even", 1), ("odd", -1)):
        basis = null_space(p - sign * np.eye(net.size))
        got = np.sort(np.linalg.eigvalsh(basis.T @ cover @ basis))
        want = np.sort(np.linalg.eigvalsh(net.quotient_operator(sign)))
        spectra[name] = float(np.abs(got - want).max())
        spectra[name + "_dim"] = basis.shape[1]
    out["spectrum_error"] = spectra
    out["even_odd_spectra_differ"] = float(np.abs(
        np.sort(np.linalg.eigvalsh(net.quotient_operator(1)))
        - np.sort(np.linalg.eigvalsh(net.quotient_operator(-1)))).max())
    # Control: a time-preserving seam.
    d = deck(model, 1, 1, 1)
    fixed = null_space(d[:f, :f] - np.eye(f))
    restricted = fixed.T @ w[:f, :f] @ fixed
    out["control"] = {
        "symplectic_error": float(np.abs(d.T @ w @ d - w).max()),
        "fixed_dim": fixed.shape[1], "restricted_rank": rank(restricted)}
    return out


# --------------------------------------------------------------------------
# H2  Local data and the sheet Wronskian
# --------------------------------------------------------------------------

def local_bracket_report(net: Network, rng: np.random.Generator) -> dict:
    model = Model(net, components=1, respond=False)
    n, size = net.n, net.size
    f = model.field_dim
    fixed = null_space(deck(model, -1)[:f, :f] - np.eye(f))
    sheet = np.concatenate([np.arange(n), size + np.arange(n)])
    a = admissible_state(model, -1, 1, 1.0, rng)
    b = admissible_state(model, -1, 1, 1.0, rng)
    start = np.vstack([a, b])
    forward = flow(model, start, STEPS, STEP)
    backward = flow(model, start, STEPS, -STEP)

    def omega_1(states):
        xa, pa = states[:, 0, :n], states[:, 0, size:size + n]
        xb, pb = states[:, 1, :n], states[:, 1, size:size + n]
        return np.sum(xa * pb - xb * pa, axis=1)

    def omega_total(states):
        xa, pa = states[:, 0, :size], states[:, 0, size:2 * size]
        xb, pb = states[:, 1, :size], states[:, 1, size:2 * size]
        return np.sum(xa * pb - xb * pa, axis=1)

    fwd, bwd = omega_1(forward), omega_1(backward)
    return {
        "sheet_data_rank": rank(fixed[sheet]), "sheet_data_dim": 2 * n,
        "mirror_symmetry_error": float(np.abs(fwd - bwd).max()),
        "relative_variation": float(fwd.max() - fwd.min()),
        "value_at_mirror": float(fwd[0]),
        "cover_wronskian": float(np.abs(omega_total(forward)).max()),
    }


# --------------------------------------------------------------------------
# H3  Response preserves the structure
# --------------------------------------------------------------------------

def response_structure_report(net: Network, rng: np.random.Generator) -> dict:
    model = Model(net, components=1, respond=True)
    d = deck(model, -1)
    w = model.symplectic

    generic = rng.standard_normal(model.dim)
    generic[:model.field_dim] *= 0.125 / np.linalg.norm(
        generic[:model.field_dim])
    generic[model.field_dim:] *= 0.05
    once = midpoint_step(model, generic, STEP)[0]
    back = midpoint_step(model, d @ once, STEP)[0]
    reversibility = float(np.abs(d @ back - generic).max()
                          / np.abs(generic).max())

    direction = admissible_state(model, -1, 1, 1.0, rng)
    rows = []
    for amp in AMPLITUDES + (STRESS_AMPLITUDE,):
        z0 = amp * direction
        forward = flow(model, z0, STEPS, STEP)[:, 0]
        energy = model.energy(forward)
        rows.append({
            "amplitude": amp,
            "seam_defect": seam_defect(model, z0, d),
            "energy_drift": float(np.abs(energy - energy[0]).max()
                                  / energy[0]),
            "max_q": float(np.abs(forward[:, model.field_dim]).max()),
        })

    # Tangent space along a history with moving mouths.
    z0 = admissible_state(model, -1, 1, 0.125, rng, mouth=(0.05, 0.03))
    basis = fixed_space(model, -1)
    frames = tangent_flow(model, z0, basis, STEPS)
    lagrangian = max(float(np.abs(fr.T @ w @ fr).max())
                     for fr in frames.values())
    tangent_rank = min(rank(fr) for fr in frames.values())
    eps = 1e-5
    plus = flow(model, z0 + eps * basis.T, STEPS)[-1]
    minus = flow(model, z0 - eps * basis.T, STEPS)[-1]
    finite = ((plus - minus) / (2 * eps)).T
    exact = frames[STEPS]
    return {
        "reversibility_error": reversibility,
        "amplitudes": rows,
        "moving_mouth_seam_defect": seam_defect(model, z0, d),
        "phase_space_dim": model.dim, "tangent_dim": basis.shape[1],
        "tangent_rank": tangent_rank,
        "tangent_lagrangian_error": lagrangian,
        "tangent_fd_error": float(np.abs(finite - exact).max()
                                  / np.abs(exact).max()),
    }


# --------------------------------------------------------------------------
# H4  Order counting
# --------------------------------------------------------------------------

def order_counting_report(net: Network, rng: np.random.Generator) -> dict:
    model = Model(net, components=1, respond=True)
    frozen = Model(net, components=1, respond=False)
    direction = admissible_state(model, -1, 1, 1.0, rng)
    f, size = model.field_dim, net.size
    rows = []
    for amp in AMPLITUDES:
        z0 = amp * direction
        fwd = flow(model, z0, STEPS, STEP)[:, 0]
        bwd = flow(model, z0, STEPS, -STEP)[:, 0]
        ref = flow(frozen, z0, STEPS, STEP)[:, 0]
        rows.append({
            "amplitude": amp,
            "max_q": float(np.abs(fwd[:, f:f + 2]).max()),
            "field_deviation": float(np.linalg.norm(
                fwd[:, :size] - ref[:, :size], axis=1).max()),
            "mouth_asymmetry": float(np.abs(fwd[:, f] - fwd[:, f + 1]).max()),
            "mouth_mirror_error": float(np.abs(fwd[:, f + 1]
                                               - bwd[:, f]).max()),
        })

    def slopes(key):
        return [float(np.log2(rows[i + 1][key] / rows[i][key]))
                for i in range(len(rows) - 1)]

    return {"rows": rows, "q_slopes": slopes("max_q"),
            "deviation_slopes": slopes("field_deviation")}


# --------------------------------------------------------------------------
# H5  Four seam types for a charged field
# --------------------------------------------------------------------------

def seam_cell_report(net: Network) -> list:
    model = Model(net, components=2)
    w = model.symplectic
    charge = model.charge_form()
    quarter = model.rotate_quarter
    f = model.field_dim
    rows = []
    for s_t, s_f in CELLS:
        d = deck(model, s_t, s_f)
        fixed = fixed_space(model, s_t, s_f)
        ffix = field_fixed_space(model, s_t, s_f)
        restricted = fixed.T @ w @ fixed
        form = ffix.T @ charge @ ffix
        eig = np.linalg.eigvalsh(form)
        minus = deck(model, s_t, s_f, -1)
        rows.append({
            "s_t": s_t, "s_f": s_f,
            "involution_error": float(np.abs(d @ d - np.eye(model.dim)).max()),
            "symplectic_sign_error": float(np.abs(d.T @ w @ d - s_t * w).max()),
            "charge_parity_error": float(np.abs(
                d.T @ charge @ d - s_t * s_f * charge).max()),
            "phase_space_dim": model.dim, "fixed_dim": fixed.shape[1],
            "omega_on_fixed": float(np.abs(restricted).max()),
            "omega_rank_on_fixed": rank(restricted) if np.abs(
                restricted).max() > 1e-9 else 0,
            "charge_on_fixed": float(np.abs(form).max()),
            "charge_positive": int(np.sum(eig > 1e-6)),
            "charge_negative": int(np.sum(eig < -1e-6)),
            "field_fixed_dim": ffix.shape[1],
            # eta = -1 against the quarter rotation of the eta = +1 map.
            "twist_is_a_rotation_error": float(np.abs(
                minus[:f, :f] - (quarter @ d @ quarter.T)[:f, :f]).max()),
            "fix_meet_i_fix": 2 * ffix.shape[1] - rank(
                np.hstack([ffix, quarter @ ffix])),
        })
    return rows


# --------------------------------------------------------------------------
# H6  Potential parity
# --------------------------------------------------------------------------

def potential_parity_report(net: Network, rng: np.random.Generator,
                            alpha: float = 0.5) -> dict:
    out = {"alpha": alpha}
    seed_state = None
    for name, parity in (("matching", 1), ("mismatched", -1)):
        model = Model(net, components=2,
                      potential=alpha * net.potential_profile(parity))
        if seed_state is None:
            seed_state = admissible_state(model, -1, -1, 0.125, rng)
        d = deck(model, -1, -1)
        out[name + "_seam_defect"] = seam_defect(model, seed_state, d)
    return out


# --------------------------------------------------------------------------
# H7  Where the pairing is
# --------------------------------------------------------------------------

def pairing_report(net: Network, rng: np.random.Generator) -> dict:
    n, size = net.n, net.size
    sheet_sites = np.arange(n)
    even = net.potential_profile(1)

    # (a) The charge is the Wronskian pairing of the two neutral sectors.
    charged = Model(net, components=2)
    neutral = Model(net, components=1)
    fn = neutral.field_dim
    wn = neutral.symplectic[:fn, :fn]
    plus = null_space(deck(neutral, -1, 1, 1)[:fn, :fn] - np.eye(fn))
    minus = null_space(deck(neutral, -1, 1, -1)[:fn, :fn] - np.eye(fn))
    z = admissible_state(charged, -1, -1, 1.0, rng)
    re = np.concatenate([z[:size], z[2 * size:3 * size]])
    im = np.concatenate([z[size:2 * size], z[3 * size:4 * size]])
    in_plus = float(np.linalg.norm(re - plus @ (plus.T @ re)))
    in_minus = float(np.linalg.norm(im - minus @ (minus.T @ im)))
    identity = float(abs(charged.charge(z)[0] - re @ wn @ im))

    # (b) Conservation and the sheet relation, potential on, mouths moving.
    model = Model(net, components=2, potential=0.5 * even)
    z0 = admissible_state(model, -1, -1, 0.125, rng, mouth=(0.05, 0.03))
    fwd = flow(model, z0, STEPS, STEP)[:, 0]
    bwd = flow(model, z0, STEPS, -STEP)[:, 0]
    total = model.charge(fwd)
    scale = float(np.linalg.norm(z0[:model.field_dim]) ** 2)
    sheet_f = model.charge(fwd, sheet_sites)
    sheet_b = model.charge(bwd, sheet_sites)
    energy = model.energy(fwd)

    # (c) Mechanical response alone: no imaginary part, ever.
    bare = Model(net, components=2)
    direction = admissible_state(bare, -1, -1, 1.0, rng, real_only=True)
    leak, moved = [], []
    for amp in AMPLITUDES + (STRESS_AMPLITUDE,):
        hist = flow(bare, amp * direction, STEPS, STEP)[:, 0]
        leak.append(float(max(np.abs(hist[:, size:2 * size]).max(),
                              np.abs(hist[:, 3 * size:4 * size]).max())))
        moved.append(float(np.abs(hist[:, bare.field_dim:
                                       bare.field_dim + 2]).max()))

    # (d) The potential couples the sectors at first order.
    ratios = []
    for alpha in ALPHAS:
        pot = Model(net, components=2, potential=alpha * even)
        hist = flow(pot, 0.125 * direction, STEPS, STEP)[:, 0]
        imag = np.sqrt(np.sum(hist[:, size:2 * size] ** 2, axis=1)
                       + np.sum(hist[:, 3 * size:4 * size] ** 2, axis=1))
        ratios.append(float(imag.max() / 0.125))
    return {
        "re_outside_untwisted": in_plus, "im_outside_twisted": in_minus,
        "charge_minus_wronskian": identity,
        "charge_value": float(total[0]), "charge_scale": scale,
        "charge_drift": float(np.abs(total - total[0]).max() / scale),
        "sheet_relation_error": float(np.abs(
            sheet_f + sheet_b - total[0]).max() / scale),
        "sheet_charge_variation": float((sheet_f.max() - sheet_f.min())
                                        / scale),
        "energy_drift": float(np.abs(energy - energy[0]).max()
                              / abs(energy[0])),
        "mechanical_leak": leak, "mechanical_mouth_motion": moved,
        "alphas": list(ALPHAS), "potential_ratios": ratios,
        "potential_slopes": [float(np.log2(ratios[i + 1] / ratios[i]))
                             for i in range(len(ratios) - 1)],
    }


# --------------------------------------------------------------------------
# Supplementary measurements (not frozen predictions)
# --------------------------------------------------------------------------

def supplementary_report(net: Network, rng: np.random.Generator) -> dict:
    """Three measurements added after the frozen checks were first run.

    They are reported separately and gate nothing.

    * a strongly nonlinear history, because the frozen stress amplitude
      moves the mouths by less than 0.01;
    * the potential's Re-Im coupling on the ordinary handle, to see whether
      that coupling is specific to the time-reversing seam;
    * conditioning of the sheet-1 Cauchy data away from the mirror time.
    """
    n, size = net.n, net.size
    model = Model(net, components=1, respond=True)
    d = deck(model, -1)
    f = model.field_dim
    strong = []
    direction = admissible_state(model, -1, 1, 1.0, rng)
    for amp in (2.0, 4.0):
        z0 = amp * direction
        fwd = flow(model, z0, STEPS, STEP)[:, 0]
        energy = model.energy(fwd)
        barrier = np.exp(-2.0 * fwd[:, f:f + 2].max() * net.profile.max())
        strong.append({
            "amplitude": amp, "seam_defect": seam_defect(model, z0, d),
            "max_q": float(np.abs(fwd[:, f:f + 2]).max()),
            "min_barrier_factor": float(barrier),
            "energy_drift": float(np.abs(energy - energy[0]).max()
                                  / energy[0])})

    # Potential coupling on the ordinary handle (s_t = s_f = +1).
    even = net.potential_profile(1)
    bare = Model(net, components=2)
    plain = admissible_state(bare, 1, 1, 1.0, rng, real_only=True)
    ratios = []
    for alpha in ALPHAS:
        pot = Model(net, components=2, potential=alpha * even)
        hist = flow(pot, 0.125 * plain, STEPS, STEP)[:, 0]
        imag = np.sqrt(np.sum(hist[:, size:2 * size] ** 2, axis=1)
                       + np.sum(hist[:, 3 * size:4 * size] ** 2, axis=1))
        ratios.append(float(imag.max() / 0.125))

    # Sheet-1 data of admissible fields at times away from the mirror.
    frozen = Model(net, components=1, respond=False)
    basis = field_fixed_space(frozen, -1)
    hist = flow(frozen, basis.T, STEPS, STEP)
    sheet = np.concatenate([np.arange(n), size + np.arange(n)])
    conditioning = []
    for step in (0, 50, 100, 200):
        s = np.linalg.svd(hist[step][:, sheet].T, compute_uv=False)
        conditioning.append({"time": step * STEP, "rank": rank(
            hist[step][:, sheet].T), "smallest_singular": float(s[-1]),
            "largest_singular": float(s[0])})
    return {
        "strong_response": strong,
        "ordinary_handle_potential_ratios": ratios,
        "ordinary_handle_potential_slopes": [
            float(np.log2(ratios[i + 1] / ratios[i]))
            for i in range(len(ratios) - 1)],
        "sheet_data_conditioning": conditioning,
    }


# --------------------------------------------------------------------------
# Verdict
# --------------------------------------------------------------------------

INCOMPLETE = "INCOMPLETE"
UNRESOLVED_NUMERICALLY = "UNRESOLVED_NUMERICALLY"
COUNTEREXAMPLE = "COUNTEREXAMPLE"

VERDICT_FIELDS = ("numerical_validity", "support_dynamics",
                  "electric_response", "localized_structure", "local_bracket",
                  "neutral_sector_pairing", "charged_sector_pairing",
                  "positive_response", "spinor_matching", "event_law",
                  "bam_cell")
CONFIRMED = {
    "localized_structure":
        "LAGRANGIAN_HALF_SURVIVES_A_LOCALIZED_RESPONDING_HANDLE",
    "local_bracket": "SHEET_DATA_FREE_AT_MIRROR_WRONSKIAN_NOT_CONSERVED",
    "neutral_sector_pairing": "NOT_SUPPLIED_BY_MECHANICAL_RESPONSE",
    "charged_sector_pairing":
        "CHARGE_PAIRS_THE_SECTORS_ONLY_WITH_BOTH_REVERSALS",
    # The H8 checks establish that the pairing is indefinite and that no
    # complex structure acts on the admissible set. That is what
    # ``NOT_DERIVED`` records here.
    "positive_response": "NOT_DERIVED",
}
FIXED = {
    "support_dynamics": "EFFECTIVE_SUPPORT_ASSUMED",
    "electric_response": "EXTERNAL_POTENTIAL_PRESCRIBED",
    "spinor_matching": "NOT_DERIVED",
    "event_law": "NOT_DERIVED",
    "bam_cell": "NOT_ESTABLISHED",
}


def verdict(checks: dict, groups: dict, counterexamples: dict | None = None
            ) -> dict:
    """Status per field; a missed tolerance is not a refutation."""
    counterexamples = counterexamples or {}
    out, failed, missing_all = {}, [], []
    for name, names in groups.items():
        missing = [c for c in names if c not in checks]
        bad = [c for c in names if c in checks and not checks[c]]
        if name in counterexamples:
            out[name] = COUNTEREXAMPLE
        elif missing or not names:
            out[name] = INCOMPLETE
        elif bad:
            out[name] = UNRESOLVED_NUMERICALLY
        else:
            out[name] = CONFIRMED[name]
        failed += bad
        missing_all += missing
    absent = [g for g in CONFIRMED if g not in groups]
    for name in absent:
        out[name] = INCOMPLETE
    out.update(FIXED)
    out["numerical_validity"] = (
        INCOMPLETE if missing_all or absent
        else UNRESOLVED_NUMERICALLY if failed else "PASS")
    out["failed_checks"] = failed
    out["missing_checks"] = missing_all
    out["counterexamples"] = dict(counterexamples)
    return out
