"""A responding electric field on the Klein handle. Freeze ``4ab8bcd``.

Corrected after review of ``562af71``: trapped flux is sufficient but not
necessary for sheet pairing, the frozen tangent test left the Gauss surface,
and the imaginary part of the field is gauge dependent.

Lattice scalar electrodynamics in temporal gauge on the cover network of
``klein_handle_response``, with the two waist coordinates of that round:

    H = 1/2 sum_i |Pi_i|^2 + 1/2 sum_{e: i->j} |Psi_j - exp(i g a_e) Psi_i|^2
        + 1/2 sum_i (m^2 + V_i(q)) |Psi_i|^2 + 1/2 sum_e E_e^2
        + sum_tubes (p_q^2 / 2M + K_s q^2 / 2)

``a_e`` is the link potential, ``E_e`` its conjugate electric field. The
Gauss function ``G_i = (div E)_i - g rho_i`` is conserved. Nothing is
prescribed: field, gauge field and mouths are solved from this one ``H``.

On a graph the gauge field has one free mode per independent cycle. Flux
around the handle loop is source-free flux trapped by the topology.

State layout: ``[x, y, p_x, p_y (N each), a, E (Ne each), q (2), p_q (2)]``.

See ``docs/klein_handle_electric_prereg.md`` and
``docs/klein_handle_electric.md``. Not Stage 1 of the boundary-pairing
protocol.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property

import numpy as np

from geometrodynamics.bulk import klein_handle_response as kh

SEED = 2026111003
COUPLING = 1.0
STEP = kh.STEP
STEPS = kh.STEPS
AMPLITUDES = kh.AMPLITUDES
STRESS_AMPLITUDE = kh.STRESS_AMPLITUDE
FLUXES = (1e-3, 2e-3, 4e-3)
CELLS = kh.CELLS
PRIMARY, REPLICATE = kh.PRIMARY, kh.REPLICATE

INCOMPLETE = kh.INCOMPLETE
UNRESOLVED_NUMERICALLY = kh.UNRESOLVED_NUMERICALLY
COUNTEREXAMPLE = kh.COUNTEREXAMPLE


@dataclass(frozen=True)
class Electric:
    """Charged field, link gauge field and two mouths on the cover."""

    net: kh.Network
    coupling: float = COUPLING

    # ---- layout ----------------------------------------------------------

    @property
    def n_sites(self) -> int:
        return self.net.size

    @cached_property
    def n_edges(self) -> int:
        return len(self.net.edges)

    @property
    def dim(self) -> int:
        return 4 * self.n_sites + 2 * self.n_edges + 4

    @cached_property
    def slices(self) -> dict:
        n, m = self.n_sites, self.n_edges
        marks = np.cumsum([0, n, n, n, n, m, m, 2, 2])
        names = ("x", "y", "px", "py", "a", "E", "q", "pq")
        return {name: slice(int(marks[i]), int(marks[i + 1]))
                for i, name in enumerate(names)}

    def part(self, z: np.ndarray, name: str) -> np.ndarray:
        return np.atleast_2d(z)[:, self.slices[name]]

    # ---- graph operators -------------------------------------------------

    @cached_property
    def tail(self) -> np.ndarray:
        return np.array([i for i, _ in self.net.edges])

    @cached_property
    def head(self) -> np.ndarray:
        return np.array([j for _, j in self.net.edges])

    @cached_property
    def tail_matrix(self) -> np.ndarray:
        t = np.zeros((self.n_edges, self.n_sites))
        t[np.arange(self.n_edges), self.tail] = 1.0
        return t

    @cached_property
    def head_matrix(self) -> np.ndarray:
        h = np.zeros((self.n_edges, self.n_sites))
        h[np.arange(self.n_edges), self.head] = 1.0
        return h

    @cached_property
    def divergence(self) -> np.ndarray:
        """``(div E)_i = sum_out E - sum_in E``; shape ``(N, Ne)``."""
        return (self.tail_matrix - self.head_matrix).T

    @cached_property
    def degree(self) -> np.ndarray:
        return (self.tail_matrix + self.head_matrix).sum(axis=0)

    @cached_property
    def edge_swap(self) -> np.ndarray:
        """Permutation of oriented edges induced by the sheet swap."""
        half = self.n_edges // 2
        p = np.zeros((self.n_edges, self.n_edges))
        for k in range(self.n_edges):
            p[k, (k + half) % self.n_edges] = 1.0
        return p

    @cached_property
    def seam_edges(self) -> tuple:
        """Indices of the last edge of ``theta`` and of ``theta'``."""
        half = self.n_edges // 2
        return half - 1, self.n_edges - 1

    @cached_property
    def tube_edges(self) -> np.ndarray:
        """Edges with at least one tube site."""
        tube = set(self.net.tube_sites.ravel().tolist())
        return np.array([k for k, (i, j) in enumerate(self.net.edges)
                         if i in tube or j in tube])

    @cached_property
    def handle_cycle(self) -> np.ndarray:
        """Unit flux around the lift of the handle loop.

        ``theta``, then the forward ring arc ``b -> a`` on sheet 2, then
        ``theta'``, then the same arc on sheet 1. Every edge is traversed
        along its orientation, so the vector is an indicator.
        """
        net, half = self.net, self.n_edges // 2
        c = np.zeros(self.n_edges)
        for s in (0, 1):
            base = s * half
            c[base + net.n_ext:base + half] = 1.0          # tube of sheet s
            e = net.port_b
            while e % net.n_ext != net.port_a:
                c[base + e % net.n_ext] = 1.0              # ring arc b -> a
                e += 1
        return c

    # ---- dynamics --------------------------------------------------------

    def barrier(self, q: np.ndarray) -> np.ndarray:
        v = np.zeros((q.shape[0], self.n_sites))
        for t in (0, 1):
            v[:, self.net.tube_sites[t]] = kh.BARRIER * np.exp(
                -2.0 * q[:, t, None] * self.net.profile)
        return v

    def _fields(self, z):
        z = np.atleast_2d(z)
        s = self.slices
        psi = z[:, s["x"]] + 1j * z[:, s["y"]]
        mom = z[:, s["px"]] + 1j * z[:, s["py"]]
        return z, psi, mom, z[:, s["a"]], z[:, s["E"]], z[:, s["q"]], z[:, s["pq"]]

    def vector_field(self, z: np.ndarray) -> np.ndarray:
        z, psi, mom, a, e, q, pq = self._fields(z)
        g = self.coupling
        link = np.exp(1j * g * a)
        carried = link * psi[:, self.tail]            # U Psi_i, on each edge
        at_head = psi[:, self.head]
        hop = carried @ self.head_matrix + (np.conj(link) * at_head
                                            ) @ self.tail_matrix
        v = self.barrier(q)
        force = -(self.degree * psi - hop) - (kh.MASS_SQ + v) * psi
        current = np.imag(np.conj(carried) * at_head)
        density = np.abs(psi) ** 2
        dpq = np.empty_like(pq)
        for t in (0, 1):
            sites = self.net.tube_sites[t]
            dpq[:, t] = -kh.MOUTH_STIFFNESS * q[:, t] + np.sum(
                self.net.profile * v[:, sites] * density[:, sites], axis=1)
        return np.concatenate([mom.real, mom.imag, force.real, force.imag,
                               e, g * current, pq / kh.MOUTH_MASS, dpq],
                              axis=1)

    def energy(self, z: np.ndarray) -> np.ndarray:
        z, psi, mom, a, e, q, pq = self._fields(z)
        link = np.exp(1j * self.coupling * a)
        hop = psi[:, self.head] - link * psi[:, self.tail]
        v = self.barrier(q)
        return (0.5 * np.sum(np.abs(mom) ** 2, axis=1)
                + 0.5 * np.sum(np.abs(hop) ** 2, axis=1)
                + 0.5 * np.sum((kh.MASS_SQ + v) * np.abs(psi) ** 2, axis=1)
                + 0.5 * np.sum(e ** 2, axis=1)
                + np.sum(pq ** 2, axis=1) / (2 * kh.MOUTH_MASS)
                + 0.5 * kh.MOUTH_STIFFNESS * np.sum(q ** 2, axis=1))

    def charge_density(self, z: np.ndarray) -> np.ndarray:
        z, psi, mom, *_ = self._fields(z)
        return np.imag(np.conj(psi) * mom)

    def gauss(self, z: np.ndarray) -> np.ndarray:
        """``G_i = (div E)_i - g rho_i``."""
        e = self.part(z, "E")
        return e @ self.divergence.T - self.coupling * self.charge_density(z)

    def sheet_charge(self, z: np.ndarray) -> np.ndarray:
        return self.charge_density(z)[:, :self.net.n].sum(axis=1)

    def seam_flux(self, z: np.ndarray) -> np.ndarray:
        """Electric flux out of sheet 1 through the two seam edges."""
        e = self.part(z, "E")
        first, second = self.seam_edges
        return e[:, first] - e[:, second]

    def handle_flux(self, z: np.ndarray) -> np.ndarray:
        """Component of ``E`` along the handle-loop cycle."""
        c = self.handle_cycle
        return self.part(z, "E") @ c / (c @ c)

    def link_current(self, z: np.ndarray) -> np.ndarray:
        """Gauge-invariant ``J_e = Im(conj(U_e Psi_i) Psi_j)`` on each edge."""
        z, psi, _, a, *_ = self._fields(z)
        carried = np.exp(1j * self.coupling * a) * psi[:, self.tail]
        return np.imag(np.conj(carried) * psi[:, self.head])

    def loop_holonomy(self, z: np.ndarray) -> np.ndarray:
        """``sum of a_e`` around the handle loop; gauge invariant."""
        return self.part(z, "a") @ self.handle_cycle

    def gauss_differential(self, z: np.ndarray) -> np.ndarray:
        """``dG`` at ``z``, shape ``(N, dim)``.

        ``G`` is quadratic in the state, so a central difference is exact.
        """
        z = np.asarray(z, dtype=float).ravel()
        eye = np.eye(self.dim)
        return (self.gauss(z + eye) - self.gauss(z - eye)).T / 2.0

    # ---- bilinear forms and symmetries -----------------------------------

    @cached_property
    def symplectic(self) -> np.ndarray:
        w = np.zeros((self.dim, self.dim))
        s = self.slices
        for pos, mom in (("x", "px"), ("y", "py"), ("a", "E"), ("q", "pq")):
            size = s[pos].stop - s[pos].start
            w[s[pos], s[mom]] = np.eye(size)
            w[s[mom], s[pos]] = -np.eye(size)
        return w

    def deck(self, s_t: int, s_f: int) -> np.ndarray:
        """Deck map: field and mouths as in round 13, and
        ``a -> s_f a``, ``E -> s_t s_f E`` on the swapped edge."""
        p, pe, s = self.net.swap, self.edge_swap, self.slices
        d = np.zeros((self.dim, self.dim))
        d[s["x"], s["x"]] = p
        d[s["y"], s["y"]] = s_f * p
        d[s["px"], s["px"]] = s_t * p
        d[s["py"], s["py"]] = s_t * s_f * p
        d[s["a"], s["a"]] = s_f * pe
        d[s["E"], s["E"]] = s_t * s_f * pe
        flip = np.array([[0.0, 1.0], [1.0, 0.0]])
        d[s["q"], s["q"]] = flip
        d[s["pq"], s["pq"]] = s_t * flip
        return d

    @cached_property
    def conjugation(self) -> np.ndarray:
        """Charge conjugation ``C``: flips ``y, p_y, a, E``."""
        sign = np.ones(self.dim)
        for name in ("y", "py", "a", "E"):
            sign[self.slices[name]] = -1.0
        return np.diag(sign)

    def gauge_transform(self, z: np.ndarray, chi: np.ndarray) -> np.ndarray:
        """Time-independent gauge transformation with site phases ``chi``."""
        z = np.atleast_2d(z).copy()
        s, g = self.slices, self.coupling
        phase = np.exp(1j * g * chi)
        psi = (z[:, s["x"]] + 1j * z[:, s["y"]]) * phase
        mom = (z[:, s["px"]] + 1j * z[:, s["py"]]) * phase
        z[:, s["x"]], z[:, s["y"]] = psi.real, psi.imag
        z[:, s["px"]], z[:, s["py"]] = mom.real, mom.imag
        z[:, s["a"]] += chi[self.head] - chi[self.tail]
        return z


# --------------------------------------------------------------------------
# Integrator and admissible data
# --------------------------------------------------------------------------

def flow(model: Electric, z0: np.ndarray, steps: int, h: float = STEP
         ) -> np.ndarray:
    """Implicit midpoint; shape ``(steps + 1, B, dim)``."""
    return kh.flow(model, z0, steps, h)


def admissible_state(model: Electric, s_t: int, s_f: int, amplitude: float,
                     rng: np.random.Generator, real_only: bool = False,
                     flux: float = 0.0, link_scale: float = 0.0,
                     mouth: tuple = (0.0, 0.0)) -> np.ndarray:
    """Deck-fixed mirror data on the Gauss surface.

    The field is drawn from the deck-fixed field space. Unless ``real_only``
    its imaginary part is adjusted so the total charge vanishes, which a
    closed network requires. ``E`` is the minimum-norm solution of Gauss's
    law plus ``flux`` around the handle loop.
    """
    s, n = model.slices, model.n_sites
    d = model.deck(s_t, s_f)
    block = slice(0, 4 * n)
    basis = kh.null_space(d[block, block] - np.eye(4 * n))

    def draw():
        v = basis @ rng.standard_normal(basis.shape[1])
        out = np.zeros(model.dim)
        out[block] = v
        return out

    z = draw()
    if real_only:
        z[s["y"]] = 0.0
        z[s["py"]] = 0.0
    else:
        other = draw()

        def total(re, im):
            return re[s["x"]] @ im[s["py"]] - im[s["y"]] @ re[s["px"]]

        q0, q1 = total(z, z), total(z, other)
        if abs(q0) > 1e-14:
            ratio = q0 / q1
            z[s["y"]] -= ratio * other[s["y"]]
            z[s["py"]] -= ratio * other[s["py"]]
    z[block] *= amplitude / np.linalg.norm(z[block])
    rho = model.charge_density(z)[0]
    z[s["E"]] = np.linalg.lstsq(model.divergence, model.coupling * rho,
                                rcond=None)[0] + flux * model.handle_cycle
    if link_scale:
        raw = link_scale * rng.standard_normal(model.n_edges)
        z[s["a"]] = 0.5 * (raw + s_f * model.edge_swap @ raw)
    z[s["q"]] = mouth[0]
    z[s["pq"]] = (mouth[1], s_t * mouth[1])
    return z


def two_sided(model: Electric, z0: np.ndarray, steps: int = STEPS):
    """Forward and backward histories of a batch of mirror states."""
    return flow(model, z0, steps, STEP), flow(model, z0, steps, -STEP)


# --------------------------------------------------------------------------
# E1  Structure
# --------------------------------------------------------------------------

def deck_report(model: Electric, rng: np.random.Generator) -> list:
    w = model.symplectic
    probe = 0.3 * rng.standard_normal((4, model.dim))
    rows = []
    for s_t, s_f in CELLS:
        d = model.deck(s_t, s_f)
        mapped = probe @ d.T
        e0, e1 = model.energy(probe), model.energy(mapped)
        gauss_error = float(np.abs(
            model.gauss(mapped)
            - s_t * s_f * model.gauss(probe) @ model.net.swap.T).max())
        fixed = kh.null_space(d - np.eye(model.dim))
        restricted = fixed.T @ w @ fixed
        rows.append({
            "s_t": s_t, "s_f": s_f,
            "involution_error": float(np.abs(d @ d - np.eye(model.dim)).max()),
            "symplectic_sign_error": float(np.abs(d.T @ w @ d - s_t * w).max()),
            "energy_invariance_error": float(np.abs(e1 - e0).max()
                                             / np.abs(e0).max()),
            "gauss_covariance_error": gauss_error,
            "phase_space_dim": model.dim, "fixed_dim": fixed.shape[1],
            "omega_on_fixed": float(np.abs(restricted).max()),
            "omega_rank_on_fixed": kh.rank(restricted) if np.abs(
                restricted).max() > 1e-9 else 0,
        })
    return rows


def structure_report(model: Electric, rng: np.random.Generator) -> dict:
    out = {"cells": []}
    for s_t, s_f in CELLS:
        d = model.deck(s_t, s_f)
        rows = []
        for amp in AMPLITUDES:
            z0 = admissible_state(model, s_t, s_f, amp, rng, link_scale=0.05,
                                  mouth=(0.02, 0.01))
            fwd, bwd = two_sided(model, z0)
            fwd, bwd = fwd[:, 0], bwd[:, 0]
            scale = float(np.linalg.norm(z0))
            if s_t == -1:
                defect = float(np.abs(fwd - bwd @ d.T).max()) / scale
            else:
                defect = float(max(np.abs(fwd - fwd @ d.T).max(),
                                   np.abs(bwd - bwd @ d.T).max())) / scale
            energy = model.energy(fwd)
            rows.append({
                "amplitude": amp,
                "initial_fixed_error": float(np.abs(d @ z0 - z0).max()),
                "symmetry_defect": defect,
                "gauss_max": float(max(np.abs(model.gauss(fwd)).max(),
                                       np.abs(model.gauss(bwd)).max())),
                "total_charge": float(np.abs(
                    model.charge_density(fwd).sum(axis=1)).max()),
                "local_charge": float(np.abs(model.charge_density(fwd)).max()),
                "energy_drift": float(np.abs(energy - energy[0]).max()
                                      / energy[0]),
                "max_mouth": float(np.abs(model.part(fwd, "q")).max()),
                "max_field_strength": float(np.abs(model.part(fwd, "E")).max()),
            })
        out["cells"].append({"s_t": s_t, "s_f": s_f, "amplitudes": rows})

    # Reversibility and the tangent space in the doubly reversed cell.
    d = model.deck(-1, -1)
    generic = 0.1 * rng.standard_normal(model.dim)
    once = kh.midpoint_step(model, generic, STEP)[0]
    back = kh.midpoint_step(model, d @ once, STEP)[0]
    out["reversibility_error"] = float(np.abs(d @ back - generic).max()
                                       / np.abs(generic).max())
    z0 = admissible_state(model, -1, -1, 0.125, rng, link_scale=0.05,
                          mouth=(0.02, 0.01))
    basis = kh.null_space(d - np.eye(model.dim))
    eps = 1e-5
    plus = flow(model, z0 + eps * basis.T, STEPS)[-1]
    minus = flow(model, z0 - eps * basis.T, STEPS)[-1]
    frame = ((plus - minus) / (2 * eps)).T
    w = model.symplectic
    out["tangent_dim"] = basis.shape[1]
    out["tangent_rank"] = kh.rank(frame)
    out["tangent_lagrangian_error"] = float(np.abs(frame.T @ w @ frame).max())
    # Sanity scale: the same frame against a symplectic partner is O(1).
    out["tangent_frame_scale"] = float(np.abs(frame).max())
    return out


# --------------------------------------------------------------------------
# E2  Admissible trapped flux
# --------------------------------------------------------------------------

def flux_admissibility_report(model: Electric) -> dict:
    cycles = kh.null_space(model.divergence)
    pe = model.edge_swap
    out = {"cycle_dim": cycles.shape[1],
           "handle_cycle_divergence": float(np.abs(
               model.divergence @ model.handle_cycle).max()),
           "handle_cycle_swap_error": float(np.abs(
               pe @ model.handle_cycle - model.handle_cycle).max()),
           "cells": []}
    for s_t, s_f in CELLS:
        sign = s_t * s_f
        # Cycles with E o sigma = sign E.
        block = cycles.T @ pe @ cycles
        vals, vecs = np.linalg.eigh((block + block.T) / 2)
        keep = cycles @ vecs[:, np.abs(vals - sign) < 1e-9]
        on_tubes = float(np.abs(keep[model.tube_edges]).max()) \
            if keep.shape[1] else 0.0
        out["cells"].append({"s_t": s_t, "s_f": s_f,
                             "admissible_flux_dim": keep.shape[1],
                             "max_on_tube_edges": on_tubes})
    return out


# --------------------------------------------------------------------------
# E3  No coupling from response alone
# --------------------------------------------------------------------------

def spontaneous_report(model: Electric, rng: np.random.Generator,
                       s_t: int, s_f: int) -> dict:
    direction = admissible_state(model, s_t, s_f, 1.0, rng, real_only=True)
    rows = []
    for amp in AMPLITUDES + (STRESS_AMPLITUDE,):
        hist = flow(model, amp * direction, STEPS)[:, 0]
        rows.append({
            "amplitude": amp,
            "imaginary": float(max(np.abs(model.part(hist, "y")).max(),
                                   np.abs(model.part(hist, "py")).max())),
            "gauge": float(max(np.abs(model.part(hist, "a")).max(),
                               np.abs(model.part(hist, "E")).max())),
            "max_mouth": float(np.abs(model.part(hist, "q")).max()),
        })
    return {"s_t": s_t, "s_f": s_f, "rows": rows}


# --------------------------------------------------------------------------
# E4 - E7  Trapped flux
# --------------------------------------------------------------------------

def _log2_slopes(values):
    return [float(np.log2(values[i + 1] / values[i]))
            for i in range(len(values) - 1)]


def flux_report(model: Electric, rng: np.random.Generator,
                s_t: int, s_f: int) -> dict:
    """Real field plus trapped handle flux, in one seam cell."""
    s = model.slices
    direction = admissible_state(model, s_t, s_f, 1.0, rng, real_only=True)
    cycle = model.handle_cycle
    d = model.deck(s_t, s_f)

    def start(amp, phi):
        z = amp * direction
        z[s["E"]] = phi * cycle
        return z

    def imag_norm(hist):
        return np.sqrt(np.sum(model.part(hist, "y") ** 2, axis=1)
                       + np.sum(model.part(hist, "py") ** 2, axis=1))

    # E4: first order in the flux, exactly odd.
    ratios, flux_change_by_phi = [], []
    for phi in FLUXES:
        hist = flow(model, start(0.125, phi), STEPS)[:, 0]
        ratios.append(float(imag_norm(hist).max() / 0.125))
        flux_change_by_phi.append(float(np.abs(
            model.handle_flux(hist) - phi).max()))
    phi = FLUXES[-1]
    z_plus = start(0.125, phi)
    fwd, bwd = two_sided(model, np.vstack([z_plus, start(0.125, -phi)]))
    pos, neg = fwd[:, 0], fwd[:, 1]
    field_scale = 0.125
    odd_error = float(max(
        np.abs(model.part(pos, "y") + model.part(neg, "y")).max(),
        np.abs(model.part(pos, "py") + model.part(neg, "py")).max()))
    even_error = float(max(
        np.abs(model.part(pos, "x") - model.part(neg, "x")).max(),
        np.abs(model.part(pos, "px") - model.part(neg, "px")).max()))
    imag_scale = float(imag_norm(pos).max())

    # E5: the flux responds, second order in amplitude.
    change_by_amp, signed = [], []
    for amp in AMPLITUDES:
        hist = flow(model, start(amp, phi), STEPS)[:, 0]
        delta = model.handle_flux(hist) - phi
        change_by_amp.append(float(np.abs(delta).max()))
        signed.append(float(delta[np.argmax(np.abs(delta))]))

    # E6: sheet charge against the flux through the seam.
    back = bwd[:, 0]
    q1_f, q1_b = model.sheet_charge(pos), model.sheet_charge(back)
    g = model.coupling
    return {
        "s_t": s_t, "s_f": s_f, "fluxes": list(FLUXES),
        "initial_fixed_error": float(np.abs(d @ z_plus - z_plus).max()),
        "imaginary_ratios": ratios,
        "imaginary_slopes": _log2_slopes(ratios),
        "odd_error": odd_error / field_scale,
        "even_error": even_error / field_scale,
        "imaginary_scale": imag_scale,
        "flux_change_by_flux": flux_change_by_phi,
        "flux_slopes_in_flux": _log2_slopes(flux_change_by_phi),
        "flux_change_by_amplitude": change_by_amp,
        "flux_slopes_in_amplitude": _log2_slopes(change_by_amp),
        "flux_change_signed": signed,
        "screening": bool(all(v * phi < 0 for v in signed)),
        "gauss_identity_error": float(np.abs(
            g * q1_f - model.seam_flux(pos)).max()),
        "gauss_max": float(np.abs(model.gauss(pos)).max()),
        "total_charge": float(np.abs(
            model.charge_density(pos).sum(axis=1)).max()),
        "sheet_charge_max": float(np.abs(q1_f).max()),
        "sheet_charge_antisymmetry_error": float(np.abs(q1_f + q1_b).max()),
        "seam_relation_defect": float(np.abs(pos - back @ d.T).max()
                                      / np.linalg.norm(z_plus))
        if s_t == -1 else None,
    }


# --------------------------------------------------------------------------
# Review of 562af71
# --------------------------------------------------------------------------

def matter_phase_state(model: Electric, amplitude: float,
                       rng: np.random.Generator) -> np.ndarray:
    """R1. Admissible Klein data with no charge, no field and no flux.

    ``x`` swap even, ``y`` swap odd, ``Pi = a = E = 0``. The field has a
    phase structure, so the gauge-invariant link currents are not zero.
    """
    s, p = model.slices, model.net.swap
    r1, r2 = (rng.standard_normal(model.n_sites) for _ in range(2))
    z = np.zeros(model.dim)
    z[s["x"]] = 0.5 * (r1 + p @ r1)
    z[s["y"]] = 0.5 * (r2 - p @ r2)
    z *= amplitude / np.linalg.norm(z)
    return z


def matter_phase_report(model: Electric, rng: np.random.Generator) -> dict:
    """R1. Sheet pairing generated with no trapped flux at any time."""
    from scipy.integrate import solve_ivp

    d = model.deck(-1, -1)
    direction = matter_phase_state(model, 1.0, rng)
    g = model.coupling
    start = 0.125 * direction
    rows = []
    for amp in AMPLITUDES:
        hist = flow(model, amp * direction, STEPS)[:, 0]
        back = flow(model, amp * direction, STEPS, -STEP)[:, 0]
        e = model.part(hist, "E")
        rows.append({
            "amplitude": amp,
            "max_sheet_charge": float(np.abs(model.sheet_charge(hist)).max()),
            "max_handle_flux": float(np.abs(model.handle_flux(hist)).max()),
            "max_loop_holonomy": float(np.abs(model.loop_holonomy(hist)).max()),
            "max_field": float(np.abs(e).max()),
            "field_swap_odd_error": float(np.abs(e @ model.edge_swap.T
                                                 + e).max()),
            "gauss_max": float(np.abs(model.gauss(hist)).max()),
            "gauss_identity_error": float(np.abs(
                g * model.sheet_charge(hist) - model.seam_flux(hist)).max()),
            "seam_defect": float(np.abs(hist - back @ d.T).max() / amp),
        })
    charges = [r["max_sheet_charge"] for r in rows]

    # Charge conjugation reverses the sheet charge: the matter phase orients.
    hist = flow(model, start, STEPS)[:, 0]
    conj = flow(model, model.conjugation @ start, STEPS)[:, 0]
    q1 = model.sheet_charge(hist)

    # With the coupling off there is still charge transport, and no field.
    free = Electric(model.net, coupling=0.0)
    free_hist = flow(free, start, STEPS)[:, 0]

    # Step refinement and an independent integrator.
    fine = flow(model, start, 2 * STEPS, STEP / 2)[:, 0]
    reference = solve_ivp(
        lambda t, y: model.vector_field(y)[0], (0.0, STEPS * STEP), start,
        method="DOP853", rtol=1e-11, atol=1e-13,
        t_eval=np.linspace(0.0, STEPS * STEP, 2 * STEPS + 1))
    ref_q1 = model.sheet_charge(reference.y.T)
    return {
        "initial": {
            "deck_fixed_error": float(np.abs(d @ start - start).max()),
            "charge_density": float(np.abs(model.charge_density(start)).max()),
            "field": float(np.abs(model.part(start, "E")).max()),
            "momentum": float(np.abs(start[model.slices["px"].start:
                                           model.slices["py"].stop]).max()),
            "link_potential": float(np.abs(model.part(start, "a")).max()),
            "max_link_current": float(np.abs(model.link_current(start)).max()),
        },
        "rows": rows,
        "sheet_charge_slopes": _log2_slopes(charges),
        "conjugate_sum": float(np.abs(model.sheet_charge(conj) + q1).max()),
        "uncoupled_max_sheet_charge": float(np.abs(
            free.sheet_charge(free_hist)).max()),
        "uncoupled_max_field": float(np.abs(free.part(free_hist, "E")).max()),
        "step_refinement_change": float(
            abs(np.abs(model.sheet_charge(fine)).max() - np.abs(q1).max())
            / np.abs(q1).max()),
        "dop853_max_sheet_charge": float(np.abs(ref_q1).max()),
        "midpoint_max_sheet_charge": float(np.abs(
            model.sheet_charge(fine)).max()),
        "dop853_disagreement": float(
            np.abs(ref_q1 - model.sheet_charge(fine)).max()
            / np.abs(ref_q1).max()),
        "dop853_gauss_max": float(np.abs(model.gauss(reference.y.T)).max()),
    }


def constrained_tangent_report(model: Electric, rng: np.random.Generator
                               ) -> dict:
    """R2. The tangent test on the Gauss surface and in the reduced space.

    The frozen test perturbed along every deck-fixed direction. Those leave
    ``G = 0`` at first order. Here the directions satisfy ``dG delta = 0``,
    the gauge directions are identified, and the count is done in the
    reduced space ``{G = 0} / gauge``.
    """
    d = model.deck(-1, -1)
    w = model.symplectic
    z0 = admissible_state(model, -1, -1, 0.125, rng, link_scale=0.05,
                          mouth=(0.02, 0.01))
    dg = model.gauss_differential(z0)
    ambient = kh.null_space(d - np.eye(model.dim))
    surface = kh.null_space(dg)                     # tangent to {G = 0}
    gauge = w @ dg.T                                # Hamiltonian fields of G_i
    fixed = kh.null_space(np.vstack([d - np.eye(model.dim), dg]))
    n_gauss = kh.rank(dg)
    restricted = surface.T @ w @ surface
    fixed_gauge = 2 * 0 + fixed.shape[1] + kh.rank(gauge) - kh.rank(
        np.hstack([fixed, gauge]))                  # dim(fixed meet gauge)
    reduced_dim = model.dim - 2 * n_gauss

    eps = 1e-5
    plus = flow(model, z0 + eps * fixed.T, STEPS)[-1]
    minus = flow(model, z0 - eps * fixed.T, STEPS)[-1]
    frame = ((plus - minus) / (2 * eps)).T
    end = flow(model, z0, STEPS)[-1, 0]
    dg_end = model.gauss_differential(end)
    return {
        "phase_space_dim": model.dim, "gauss_rank": n_gauss,
        "ambient_fixed_dim": ambient.shape[1],
        "ambient_linearized_gauss_violation": float(np.abs(dg @ ambient).max()),
        "surface_dim": surface.shape[1],
        "presymplectic_nullity": surface.shape[1] - kh.rank(restricted),
        "gauge_in_surface_error": float(np.abs(dg @ gauge).max()),
        "gauge_in_presymplectic_kernel": float(np.abs(
            surface.T @ w @ gauge).max()),
        "gauge_rank": kh.rank(gauge),
        "constrained_fixed_dim": fixed.shape[1],
        "constrained_omega": float(np.abs(fixed.T @ w @ fixed).max()),
        "gauge_directions_in_fixed": int(fixed_gauge),
        "reduced_dim": reduced_dim,
        "reduced_fixed_dim": fixed.shape[1] - int(fixed_gauge),
        "propagated_gauss_violation": float(np.abs(dg_end @ frame).max()),
        "propagated_omega": float(np.abs(frame.T @ w @ frame).max()),
        "propagated_rank": kh.rank(frame),
    }


def gauge_report(model: Electric, rng: np.random.Generator) -> dict:
    """R3. What is and is not gauge invariant in the flux experiment.

    A time-independent gauge transformation with swap-odd phases preserves
    the strict Klein deck condition. It changes the imaginary part of the
    field and leaves every gauge-invariant observable alone.
    """
    s, p = model.slices, model.net.swap
    d = model.deck(-1, -1)
    cycle = model.handle_cycle
    direction = admissible_state(model, -1, -1, 1.0, rng, real_only=True)
    raw = rng.standard_normal(model.n_sites)
    chi = 0.3 * 0.5 * (raw - p @ raw)

    def start(phi):
        z = 0.125 * direction
        z[s["E"]] = phi * cycle
        return z

    def imag_norm(h):
        return np.sqrt(np.sum(model.part(h, "y") ** 2, axis=1)
                       + np.sum(model.part(h, "py") ** 2, axis=1))

    def observables(h):
        return {"field": model.part(h, "E"), "current": model.link_current(h),
                "density": model.charge_density(h),
                "sheet": model.sheet_charge(h)[:, None],
                "holonomy": model.loop_holonomy(h)[:, None],
                "energy": model.energy(h)[:, None]}

    phi = FLUXES[-1]
    z0 = start(phi)
    moved = model.gauge_transform(z0, chi)[0]
    plain, gauged = flow(model, z0, STEPS)[:, 0], flow(model, moved, STEPS)[:, 0]
    a, b = observables(plain), observables(gauged)
    invariant = {k: float(np.abs(a[k] - b[k]).max()) for k in a}

    # Gauge-invariant version of E4: currents, sheet charge and holonomy.
    currents, sheets, holonomies = [], [], []
    for f in FLUXES:
        h = flow(model, start(f), STEPS)[:, 0]
        currents.append(float(np.abs(model.link_current(h)).max()))
        sheets.append(float(np.abs(model.sheet_charge(h)).max()))
        holonomies.append(float(np.abs(model.loop_holonomy(h)).max()))
    neg = flow(model, start(-phi), STEPS)[:, 0]
    na = observables(neg)

    # Gauge-invariant version of E3: real data, no flux.
    quiet = flow(model, start(0.0), STEPS)[:, 0]
    return {
        "gauged_deck_fixed_error": float(np.abs(d @ moved - moved).max()),
        "imaginary_before": float(imag_norm(z0[None])[0]),
        "imaginary_after": float(imag_norm(moved[None])[0]),
        "imaginary_history_difference": float(np.abs(
            imag_norm(plain) - imag_norm(gauged)).max()),
        "invariant_differences": invariant,
        "current_scale": float(np.abs(a["current"]).max()),
        "current_slopes": _log2_slopes(currents),
        "sheet_charge_slopes": _log2_slopes(sheets),
        "holonomy_slopes": _log2_slopes(holonomies),
        "current_odd_error": float(np.abs(a["current"] + na["current"]).max()),
        "sheet_charge_odd_error": float(np.abs(a["sheet"] + na["sheet"]).max()),
        "holonomy_odd_error": float(np.abs(
            a["holonomy"] + na["holonomy"]).max()),
        "quiet_current": float(np.abs(model.link_current(quiet)).max()),
        "quiet_density": float(np.abs(model.charge_density(quiet)).max()),
        "quiet_field": float(np.abs(model.part(quiet, "E")).max()),
    }


# --------------------------------------------------------------------------
# Verdict
# --------------------------------------------------------------------------

VERDICT_FIELDS = ("numerical_validity", "support_dynamics",
                  "electric_dynamics", "structure", "constrained_structure",
                  "flux_admissibility", "spontaneous_pairing",
                  "matter_phase_route", "flux_pairing", "gauge_dependence",
                  "flux_dynamics", "pairing_carrier", "seam_specificity",
                  "orienting_datum", "positive_response", "spinor_matching",
                  "event_law", "bam_cell")
# Three labels were narrowed and three fields added after review of
# ``562af71``; see the write-up.
CONFIRMED = {
    "structure": "AMBIENT_LAGRANGIAN_HALF_SURVIVES_ELECTRIC_RESPONSE",
    "constrained_structure": "LAGRANGIAN_IN_THE_GAUSS_REDUCED_SPACE",
    "flux_admissibility": "HANDLE_FLUX_ADMISSIBLE_IFF_CHARGE_UNTWISTED",
    "spontaneous_pairing": "NOT_GENERATED_FROM_REAL_DATA_WITHOUT_FLUX",
    "matter_phase_route": "MATTER_CURRENT_GENERATES_SHEET_PAIRING_WITHOUT_FLUX",
    "flux_pairing": "TRAPPED_FLUX_DRIVES_GAUGE_INVARIANT_CURRENT_ODD_IN_FLUX",
    "gauge_dependence": "IMAGINARY_PART_IS_GAUGE_DEPENDENT",
    "flux_dynamics": "FLUX_RESPONDS_AT_SECOND_ORDER",
    "pairing_carrier": "SHEET_PAIRING_IS_THE_SEAM_FLUX_BY_GAUSS",
    "seam_specificity": "FLUX_COUPLING_ALSO_ON_THE_ORDINARY_HANDLE",
}
FROZEN_LABELS = {
    "structure": "LAGRANGIAN_HALF_SURVIVES_ELECTRIC_RESPONSE",
    "spontaneous_pairing": "NOT_GENERATED_BY_RESPONSE_WITHOUT_FLUX",
    "flux_pairing": "TRAPPED_FLUX_COUPLES_SECTORS_AND_SETS_ORIENTATION",
    "flux_origin": "INITIAL_DATA_NOT_DERIVED",
}
FIXED = {
    "support_dynamics": "EFFECTIVE_SUPPORT_ASSUMED",
    "electric_dynamics": "LATTICE_GAUGE_FIELD_SOLVED_FROM_THE_ACTION",
    "orienting_datum":
        "CONJUGATION_ODD_INITIAL_DATA_FLUX_OR_MATTER_CURRENT_NOT_DERIVED",
    "positive_response": "NOT_DERIVED",
    "spinor_matching": "NOT_DERIVED",
    "event_law": "NOT_DERIVED",
    "bam_cell": "NOT_ESTABLISHED",
}


def verdict(checks: dict, groups: dict, counterexamples: dict | None = None
            ) -> dict:
    """Status per field; a missed tolerance is not a refutation."""
    counterexamples = counterexamples or {}
    out, failed, missing_all = {}, [], []
    for name in CONFIRMED:
        names = groups.get(name, [])
        missing = [c for c in names if c not in checks]
        bad = [c for c in names if c in checks and not checks[c]]
        if name in counterexamples:
            out[name] = COUNTEREXAMPLE
        elif missing or not names:
            out[name] = INCOMPLETE
            missing_all += missing or [f"<no checks for {name}>"]
        elif bad:
            out[name] = UNRESOLVED_NUMERICALLY
        else:
            out[name] = CONFIRMED[name]
        failed += bad
    out.update(FIXED)
    out["numerical_validity"] = (
        INCOMPLETE if missing_all
        else UNRESOLVED_NUMERICALLY if failed else "PASS")
    out["failed_checks"] = failed
    out["missing_checks"] = missing_all
    out["counterexamples"] = dict(counterexamples)
    return out


# --------------------------------------------------------------------------
# Supplementary measurements (not frozen predictions)
# --------------------------------------------------------------------------

def supplementary_report(model: Electric, rng: np.random.Generator,
                         draws: int = 12) -> dict:
    """Added after the frozen checks were first run; gates nothing.

    * The sign of the flux response was left unfrozen and came out opposite
      at the two network sizes. This surveys it over ``draws`` real fields.
    * The sheet charge on the ordinary handle, for comparison with the
      time-reversing seam, over the same fields.
    """
    s = model.slices
    phi, amp = FLUXES[-1], 0.125
    out = {"draws": draws, "flux": phi, "amplitude": amp, "cells": []}
    for s_t, s_f in ((-1, -1), (1, 1)):
        starts = []
        for _ in range(draws):
            z = admissible_state(model, s_t, s_f, amp, rng, real_only=True)
            z[s["E"]] = phi * model.handle_cycle
            starts.append(z)
        hist = flow(model, np.array(starts), STEPS)
        steps, batch = hist.shape[:2]
        flat = hist.reshape(steps * batch, -1)
        delta = model.handle_flux(flat).reshape(steps, batch) - phi
        sheet = model.sheet_charge(flat).reshape(steps, batch)
        peak = delta[np.argmax(np.abs(delta), axis=0), np.arange(batch)]
        final = delta[-1]
        out["cells"].append({
            "s_t": s_t, "s_f": s_f,
            "screening_at_peak": int(np.sum(peak * phi < 0)),
            "screening_at_final_time": int(np.sum(final * phi < 0)),
            "peak_changes": [float(v) for v in peak],
            "max_sheet_charge": float(np.abs(sheet).max()),
            "min_over_draws_of_max_sheet_charge": float(
                np.abs(sheet).max(axis=0).min()),
        })
    return out
