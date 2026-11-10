"""The matching map from the recorded BAM geometry to the seam cells.

Freeze ``3d3ee38``. Three parts.

I.   On the bulk mouth ``S^3``, with the maps the repository records: the
     antipode ``-I_4`` (angular part of PR #129's ``iota``) and the transport
     ``sigma = L_{-j}`` (PR #279), acting on fibre-charged functions.
II.  A lattice multiplet on the handle cover with a sign ``W`` around the
     doubled handle loop, to test when the transport seam admits a field.
III. The Pin sign against the trapped-flux sign, in round 14's model.

See ``docs/klein_handle_matching_prereg.md`` and
``docs/klein_handle_matching.md``.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property

import numpy as np

from geometrodynamics.bulk import klein_handle_electric as ke
from geometrodynamics.bulk import klein_handle_response as kh
from geometrodynamics.bulk import mouth_spin_frame, mouth_topology
from geometrodynamics.embedding import transport

SEED = 2026111004
MAX_DEGREE = 4
SAMPLE_POINTS = 400
STEP, STEPS = kh.STEP, kh.STEPS
AMPLITUDES = kh.AMPLITUDES
PRIMARY, REPLICATE = kh.PRIMARY, kh.REPLICATE

INCOMPLETE = kh.INCOMPLETE
UNRESOLVED_NUMERICALLY = kh.UNRESOLVED_NUMERICALLY
COUNTEREXAMPLE = kh.COUNTEREXAMPLE


# ==========================================================================
# Part I  The recorded maps on S^3
# ==========================================================================

def random_sphere(n: int, rng: np.random.Generator) -> np.ndarray:
    """``n`` points of ``S^3`` as complex pairs ``(z1, z2)``."""
    z = rng.standard_normal((n, 2)) + 1j * rng.standard_normal((n, 2))
    return z / np.linalg.norm(z, axis=1, keepdims=True)


def random_equator(n: int, rng: np.random.Generator) -> np.ndarray:
    """Points of the spatial brane mouth ``{Re q = 0}``."""
    v = rng.standard_normal((n, 3))
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    return np.stack([1j * v[:, 0], v[:, 1] + 1j * v[:, 2]], axis=1)


def as_quaternion(z: np.ndarray) -> np.ndarray:
    """``q = z1 + z2 j`` in the basis ``(1, i, j, k)``."""
    return np.stack([z[:, 0].real, z[:, 0].imag, z[:, 1].real, z[:, 1].imag],
                    axis=1)


def antipode(z: np.ndarray) -> np.ndarray:
    return -z


def sigma(z: np.ndarray) -> np.ndarray:
    """The repository's transport, point by point through its own function."""
    out = np.empty_like(z)
    for n, (z1, z2) in enumerate(z):
        out[n] = transport.orientation_reversal_on_s3(complex(z1), complex(z2))
    return out


def fibre_rotate(z: np.ndarray, phi: float) -> np.ndarray:
    return np.exp(1j * phi) * z


def hopf_base(z: np.ndarray) -> np.ndarray:
    """``h(q) = q^{-1} i q`` through the repository's spin-frame map."""
    return np.array([mouth_spin_frame.spin_frame(q)[0]
                     for q in as_quaternion(z)])


def bidegree_exponents(p: int, r: int) -> list:
    """Monomials ``z1^a z2^(p-a) conj(z1)^b conj(z2)^(r-b)``."""
    return [(a, p - a, b, r - b) for a in range(p + 1) for b in range(r + 1)]


def evaluate(p: int, r: int, z: np.ndarray) -> np.ndarray:
    cols = [z[:, 0] ** a1 * z[:, 1] ** a2
            * np.conj(z[:, 0]) ** b1 * np.conj(z[:, 1]) ** b2
            for a1, a2, b1, b2 in bidegree_exponents(p, r)]
    return np.stack(cols, axis=1)


def _fit(target: np.ndarray, values: np.ndarray):
    """Least-squares coefficients and relative residual."""
    coeff, *_ = np.linalg.lstsq(target, values, rcond=None)
    resid = float(np.linalg.norm(target @ coeff - values)
                  / max(np.linalg.norm(values), 1e-300))
    return coeff, resid


def antipode_report(rng: np.random.Generator) -> dict:
    """M1. Pullback by ``-I_4`` on every bidegree up to ``MAX_DEGREE``."""
    z = random_sphere(SAMPLE_POINTS, rng)
    worst_scalar, worst_fit, rows = 0.0, 0.0, []
    for p in range(MAX_DEGREE + 1):
        for r in range(MAX_DEGREE + 1 - p):
            coeff, resid = _fit(evaluate(p, r, z), evaluate(p, r, antipode(z)))
            k = p - r
            err = float(np.abs(coeff - (-1) ** k * np.eye(coeff.shape[0])).max())
            worst_scalar, worst_fit = max(worst_scalar, err), max(worst_fit, resid)
            rows.append({"p": p, "r": r, "charge": k, "scalar_error": err})
    phi = 0.7
    return {
        "bidegrees": len(rows), "rows": rows,
        "scalar_error": worst_scalar, "fit_residual": worst_fit,
        "is_fibre_translation_by_pi": float(np.abs(
            antipode(z) - fibre_rotate(z, np.pi)).max()),
        "fibre_commutator": float(np.abs(
            antipode(fibre_rotate(z, phi)) - fibre_rotate(antipode(z), phi)
        ).max()),
        "hopf_base_moved": float(np.abs(hopf_base(antipode(z))
                                        - hopf_base(z)).max()),
    }


def two_spheres_report(rng: np.random.Generator) -> dict:
    """M2. The spatial brane mouth against the Hopf base."""
    z = random_sphere(SAMPLE_POINTS, rng)
    e = random_equator(SAMPLE_POINTS, rng)
    p = as_quaternion(e)[:, 1:]                 # the imaginary 3-vector
    formula = 2 * p[:, [0]] * p - np.array([1.0, 0.0, 0.0])
    phi = 0.7
    matrix = mouth_topology.hopf_transport_matrix()
    return {
        "equator_real_part": float(np.abs(e[:, 0].real).max()),
        "antipode_leaves_equator": float(np.abs(antipode(e)[:, 0].real).max()),
        "hopf_two_to_one_error": float(np.abs(hopf_base(antipode(e))
                                              - hopf_base(e)).max()),
        "hopf_formula_error": float(np.abs(hopf_base(e) - formula).max()),
        "sigma_leaves_equator": float(np.abs(sigma(e)[:, 0].real).max()),
        "sigma_covers_base_antipode": float(np.abs(hopf_base(sigma(z))
                                                   + hopf_base(z)).max()),
        "sigma_reverses_fibre": float(np.abs(
            sigma(fibre_rotate(z, phi)) - fibre_rotate(sigma(z), -phi)).max()),
        "sigma_squared_is_antipode": float(np.abs(sigma(sigma(z)) + z).max()),
        "sigma_matches_recorded_matrix": float(np.abs(
            as_quaternion(sigma(z)) - as_quaternion(z) @ matrix.T).max()),
        "sigma_determinant": float(np.linalg.det(matrix)),
    }


def transport_report(rng: np.random.Generator) -> dict:
    """M3. ``sigma`` on fibre-charged functions; ``A = conj o sigma^*``."""
    z = random_sphere(SAMPLE_POINTS, rng)
    sz = sigma(z)
    rows, worst = [], {"square": 0.0, "fit": 0.0, "leak": np.inf}
    for p in range(MAX_DEGREE + 1):
        for r in range(MAX_DEGREE + 1 - p):
            k = p - r
            pulled = evaluate(p, r, sz)
            _, to_opposite = _fit(evaluate(r, p, z), pulled)
            m, fit = _fit(evaluate(p, r, z), np.conj(pulled))
            square = float(np.abs(m @ np.conj(m)
                                  - (-1) ** k * np.eye(m.shape[0])).max())
            row = {"p": p, "r": r, "charge": k, "square_error": square,
                   "lands_in_opposite_charge": to_opposite}
            if k != 0:
                _, same = _fit(evaluate(p, r, z), pulled)
                row["stays_in_same_charge"] = same
                worst["leak"] = min(worst["leak"], same)
            rows.append(row)
            worst["square"] = max(worst["square"], square)
            worst["fit"] = max(worst["fit"], fit, to_opposite)
    # The doublet in the explicit basis (z1, z2): column j holds A(basis_j).
    doublet, _ = _fit(z, np.conj(sz))
    i_sigma_y = np.array([[0.0, 1.0], [-1.0, 0.0]])
    return {
        "bidegrees": len(rows), "rows": rows,
        "square_error": worst["square"], "fit_residual": worst["fit"],
        "smallest_same_charge_residual": float(worst["leak"]),
        "doublet_basis": ["z1", "z2"],
        "doublet_matrix": doublet.real.tolist(),
        "doublet_imaginary_part": float(np.abs(doublet.imag).max()),
        "doublet_is_minus_i_sigma_y": float(np.abs(doublet + i_sigma_y).max()),
        "doublet_is_plus_i_sigma_y": float(np.abs(doublet - i_sigma_y).max()),
    }


def pr129_cell_report(net: kh.Network) -> dict:
    """M4. Time reversed, fibre kept, both twists: the charge vanishes."""
    model = kh.Model(net, components=2)
    charge, f = model.charge_form(), model.field_dim
    out = []
    for eta in (1, -1):
        d = kh.deck(model, -1, 1, eta)
        fixed = kh.null_space(d[:f, :f] - np.eye(f))
        fixed = np.vstack([fixed, np.zeros((4, fixed.shape[1]))])
        out.append({
            "eta": eta, "fixed_dim": fixed.shape[1], "field_dim": f,
            "charge_on_fixed": float(np.abs(fixed.T @ charge @ fixed).max()),
            "charge_parity_error": float(np.abs(d.T @ charge @ d
                                                + charge).max()),
            "omega_on_fixed": float(np.abs(
                fixed.T @ model.symplectic @ fixed).max()),
        })
    return {"twists": out}


# ==========================================================================
# Part II  A multiplet on the handle cover with loop sign W
# ==========================================================================

I2 = np.eye(2, dtype=complex)
PAULI = (np.array([[0, 1], [1, 0]], dtype=complex),
         np.array([[0, -1j], [1j, 0]], dtype=complex),
         np.array([[1, 0], [0, -1]], dtype=complex))


@dataclass(frozen=True)
class Multiplet:
    """``d`` complex components per site, two mouths, loop sign ``W``.

    State layout ``[Re Psi, Im Psi, Re Pi, Im Pi (d N each), q, p_q]``.
    ``W = -1`` puts a sign on the last seam edge, so the field is
    antiperiodic around the doubled handle loop.
    """

    net: kh.Network
    components: int = 2
    loop_sign: int = 1

    @property
    def n_sites(self) -> int:
        return self.net.size

    @property
    def field_dim(self) -> int:
        return 4 * self.components * self.n_sites

    @property
    def dim(self) -> int:
        return self.field_dim + 4

    @cached_property
    def operator(self) -> np.ndarray:
        """Signed graph Laplacian plus mass."""
        n = self.n_sites
        k = np.zeros((n, n))
        last = len(self.net.edges) - 1
        for index, (i, j) in enumerate(self.net.edges):
            sign = self.loop_sign if index == last else 1
            k[i, i] += 1
            k[j, j] += 1
            k[i, j] -= sign
            k[j, i] -= sign
        return k + kh.MASS_SQ * np.eye(n)

    @cached_property
    def sheet_sign(self) -> np.ndarray:
        """``+1`` on sheet 1; ``W`` on sheet 2."""
        eps = np.ones(self.n_sites)
        eps[self.net.n:] = self.loop_sign
        return eps

    # ---- packing ---------------------------------------------------------

    def unpack(self, z: np.ndarray):
        z = np.atleast_2d(z)
        b, d, n = z.shape[0], self.components, self.n_sites
        block = d * n
        psi = (z[:, :block] + 1j * z[:, block:2 * block]).reshape(b, d, n)
        mom = (z[:, 2 * block:3 * block]
               + 1j * z[:, 3 * block:4 * block]).reshape(b, d, n)
        return psi, mom, z[:, 4 * block:4 * block + 2], z[:, 4 * block + 2:]

    def pack(self, psi, mom, q, pq) -> np.ndarray:
        b = psi.shape[0]
        return np.concatenate([psi.real.reshape(b, -1), psi.imag.reshape(b, -1),
                               mom.real.reshape(b, -1), mom.imag.reshape(b, -1),
                               q, pq], axis=1)

    # ---- dynamics --------------------------------------------------------

    def barrier(self, q: np.ndarray) -> np.ndarray:
        v = np.zeros((q.shape[0], self.n_sites))
        for t in (0, 1):
            v[:, self.net.tube_sites[t]] = kh.BARRIER * np.exp(
                -2.0 * q[:, t, None] * self.net.profile)
        return v

    def vector_field(self, z: np.ndarray) -> np.ndarray:
        psi, mom, q, pq = self.unpack(z)
        v = self.barrier(q)
        force = -(psi @ self.operator.T) - v[:, None, :] * psi
        density = np.sum(np.abs(psi) ** 2, axis=1)
        dpq = np.empty_like(pq)
        for t in (0, 1):
            sites = self.net.tube_sites[t]
            dpq[:, t] = -kh.MOUTH_STIFFNESS * q[:, t] + np.sum(
                self.net.profile * v[:, sites] * density[:, sites], axis=1)
        return self.pack(mom, force, pq / kh.MOUTH_MASS, dpq)

    def energy(self, z: np.ndarray) -> np.ndarray:
        psi, mom, q, pq = self.unpack(z)
        v = self.barrier(q)
        kpsi = psi @ self.operator.T + v[:, None, :] * psi
        return (0.5 * np.sum(np.abs(mom) ** 2, axis=(1, 2))
                + 0.5 * np.sum(np.real(np.conj(psi) * kpsi), axis=(1, 2))
                + np.sum(pq ** 2, axis=1) / (2 * kh.MOUTH_MASS)
                + 0.5 * kh.MOUTH_STIFFNESS * np.sum(q ** 2, axis=1))

    # ---- linear maps on states -------------------------------------------

    def _matrix(self, fn) -> np.ndarray:
        """Matrix of a real-linear map given on ``(psi, mom, q, pq)``."""
        return np.asarray(self.pack(*fn(*self.unpack(np.eye(self.dim))))).T

    def internal(self, name: str):
        """A map on the ``d`` components, as a function of ``(B, d, N)``."""
        d = self.components
        if name == "identity":
            return lambda v: v
        if name == "conjugate":
            return np.conj
        if name == "i":
            return lambda v: 1j * v
        if d != 2:
            raise ValueError(f"{name} needs a doublet")
        if name == "J":                       # i sigma_y K
            return lambda v: np.stack([np.conj(v[:, 1]), -np.conj(v[:, 0])],
                                      axis=1)
        if name in ("sigma_x", "sigma_y", "sigma_z"):
            pauli = PAULI["xyz".index(name[-1])]
            return lambda v: 1j * np.einsum("ab,nbs->nas", pauli, v)
        raise ValueError(name)

    def spatial_deck(self, field: np.ndarray) -> np.ndarray:
        """``P diag(eps)``: sheet swap with the gauge sign of the twist."""
        return (field * self.sheet_sign) @ self.net.swap.T

    def deck(self, internal: str, s_t: int, sign: int = 1) -> np.ndarray:
        """``D = T o (B (x) P_eps)`` with ``B = sign * internal``."""
        b = self.internal(internal)

        def fn(psi, mom, q, pq):
            return (sign * b(self.spatial_deck(psi)),
                    s_t * sign * b(self.spatial_deck(mom)),
                    q[:, ::-1], s_t * pq[:, ::-1])

        return self._matrix(fn)

    def generator(self, name: str) -> np.ndarray:
        b = self.internal(name)
        return self._matrix(lambda psi, mom, q, pq: (
            b(psi), b(mom), np.zeros_like(q), np.zeros_like(pq)))

    @cached_property
    def symplectic(self) -> np.ndarray:
        w = np.zeros((self.dim, self.dim))
        half = self.field_dim // 2
        w[:half, half:2 * half] = np.eye(half)
        w[half:2 * half, :half] = -np.eye(half)
        f = self.field_dim
        w[f:f + 2, f + 2:] = np.eye(2)
        w[f + 2:, f:f + 2] = -np.eye(2)
        return w

    def noether_form(self, name: str) -> np.ndarray:
        """Symmetric matrix of ``N_X = 1/2 Omega(z, X z)``."""
        m = 0.5 * self.symplectic @ self.generator(name)
        return 0.5 * (m + m.T)

    def noether(self, name: str, z: np.ndarray) -> np.ndarray:
        z = np.atleast_2d(z)
        return np.einsum("bi,ij,bj->b", z, self.noether_form(name), z)

    def sheet_charge(self, z: np.ndarray) -> np.ndarray:
        psi, mom, *_ = self.unpack(z)
        n = self.net.n
        return np.sum(np.imag(np.conj(psi[:, :, :n]) * mom[:, :, :n]),
                      axis=(1, 2))


def flow(model, z0, steps, h=STEP):
    return kh.flow(model, z0, steps, h)


def consistency_report(net: kh.Network, rng: np.random.Generator) -> list:
    """M5. ``D^2 = B^2 W`` for the singlet and the doublet."""
    rows = []
    for d, internal, b_squared in ((1, "conjugate", 1), (2, "J", -1)):
        for w in (1, -1):
            model = Multiplet(net, components=d, loop_sign=w)
            deck = model.deck(internal, -1)
            eye = np.eye(model.dim)
            square = deck @ deck
            f = model.field_dim
            fixed = kh.null_space(deck[:f, :f] - np.eye(f))
            probe = 0.3 * rng.standard_normal((4, model.dim))
            e0, e1 = model.energy(probe), model.energy(probe @ deck.T)
            # D^2 on the mouths is always +1; the sign lives on the field.
            field_square = square[:f, :f]
            rows.append({
                "components": d, "internal": internal, "loop_sign": w,
                "predicted_square": b_squared * w,
                "square_error": float(np.abs(
                    field_square - b_squared * w * np.eye(f)).max()),
                "mouth_square_error": float(np.abs(
                    square[f:, f:] - eye[f:, f:]).max()),
                "energy_invariance_error": float(np.abs(e1 - e0).max()
                                                 / np.abs(e0).max()),
                "field_dim": f, "fixed_dim": fixed.shape[1],
            })
    return rows


def seam_charge_report(net: kh.Network) -> list:
    """M6 table. Which Noether forms survive on the fixed space."""
    rows = []
    for internal, w in (("J", -1), ("identity", 1)):
        model = Multiplet(net, components=2, loop_sign=w)
        f = model.field_dim
        forms = {name: model.noether_form(name)
                 for name in ("i", "sigma_x", "sigma_y", "sigma_z")}
        for s_t in (-1, 1):
            deck = model.deck(internal, s_t)
            fixed = kh.null_space(deck[:f, :f] - np.eye(f))
            fixed = np.vstack([fixed, np.zeros((4, fixed.shape[1]))])
            restricted = fixed.T @ model.symplectic @ fixed
            row = {"internal": internal, "loop_sign": w, "s_t": s_t,
                   "involution_error": float(np.abs(
                       deck @ deck - np.eye(model.dim)).max()),
                   "symplectic_sign_error": float(np.abs(
                       deck.T @ model.symplectic @ deck
                       - s_t * model.symplectic).max()),
                   "field_dim": f, "fixed_dim": fixed.shape[1],
                   "omega_on_fixed": float(np.abs(restricted).max())}
            for name, form in forms.items():
                eig = np.linalg.eigvalsh(fixed.T @ form @ fixed)
                row[name] = {
                    "max": float(np.abs(eig).max()),
                    "positive": int(np.sum(eig > 1e-6)),
                    "negative": int(np.sum(eig < -1e-6)),
                    "parity_error": float(min(
                        np.abs(deck.T @ form @ deck - form).max(),
                        np.abs(deck.T @ form @ deck + form).max())),
                }
            rows.append(row)
    return rows


def transport_seam_dynamics(net: kh.Network, rng: np.random.Generator) -> dict:
    """M6 first row: doublet, ``J``, ``W = -1``, time reversed, responding."""
    model = Multiplet(net, components=2, loop_sign=-1)
    deck = model.deck("J", -1)
    f = model.field_dim
    basis = kh.null_space(deck[:f, :f] - np.eye(f))
    names = ("i", "sigma_x", "sigma_y", "sigma_z")
    direction = np.zeros(model.dim)
    direction[:f] = basis @ rng.standard_normal(basis.shape[1])
    direction[:f] /= np.linalg.norm(direction[:f])
    rows = []
    for amp in AMPLITUDES:
        z0 = amp * direction
        z0[f:f + 2] = 0.02
        z0[f + 2], z0[f + 3] = 0.01, -0.01
        fwd = flow(model, z0, STEPS)[:, 0]
        bwd = flow(model, z0, STEPS, -STEP)[:, 0]
        scale = amp ** 2
        energy = model.energy(fwd)
        charge = model.noether("i", fwd)
        q1f, q1b = model.sheet_charge(fwd), model.sheet_charge(bwd)
        rows.append({
            "amplitude": amp,
            "initial_fixed_error": float(np.abs(deck @ z0 - z0).max()),
            "seam_defect": float(np.abs(fwd - bwd @ deck.T).max()
                                 / np.linalg.norm(z0)),
            "energy_drift": float(np.abs(energy - energy[0]).max()
                                  / energy[0]),
            "drift": {name: float(np.abs(model.noether(name, fwd)
                                         - model.noether(name, fwd[:1])).max()
                                  / scale) for name in names},
            "total_spin": float(max(np.abs(model.noether(name, fwd)).max()
                                    for name in names[1:]) / scale),
            "charge": float(charge[0] / scale),
            "sheet_relation_error": float(np.abs(
                q1f + q1b - model.sheet_charge(fwd[:1])[0] * 2).max() / scale),
            "sheet_charge_variation": float((q1f.max() - q1f.min()) / scale),
            "max_mouth": float(np.abs(fwd[:, f:f + 2]).max()),
        })
    w = model.symplectic
    full = kh.null_space(deck - np.eye(model.dim))
    # The opposite Pin sign: -J. Conjugate to J by the quarter rotation.
    minus = model.deck("J", -1, sign=-1)
    quarter = model.generator("i")
    quarter[f:, f:] = np.eye(4)
    return {
        "rows": rows, "phase_space_dim": model.dim,
        "fixed_dim": full.shape[1],
        "omega_on_fixed": float(np.abs(full.T @ w @ full).max()),
        "minus_j_is_quarter_rotated_j": float(np.abs(
            quarter @ deck @ quarter.T - minus).max()),
    }


# ==========================================================================
# Part III  Pin sign against trapped-flux sign (round 14's model)
# ==========================================================================

def label_independence_report(net: kh.Network, rng: np.random.Generator
                              ) -> dict:
    """M7. Which symmetry flips which sign."""
    model = ke.Electric(net)
    s = model.slices
    deck = model.deck(-1, -1)
    field_sign = np.ones(model.dim)
    for name in ("x", "y", "px", "py"):
        field_sign[s[name]] = -1.0
    deck_minus = deck * field_sign[:, None]          # the opposite Pin sign

    quarter = np.eye(model.dim)                      # Psi -> i Psi
    n = model.n_sites
    eye, zero = np.eye(n), np.zeros((n, n))
    rot = np.block([[zero, -eye], [eye, zero]])
    quarter[:2 * n, :2 * n] = rot
    quarter[2 * n:4 * n, 2 * n:4 * n] = rot
    conj = model.conjugation

    probe = 0.3 * rng.standard_normal((4, model.dim))
    energy = model.energy(probe)

    def rel(a):
        return float(np.abs(a - energy).max() / np.abs(energy).max())

    flux = model.handle_cycle
    cycle = np.zeros(model.dim)
    cycle[s["E"]] = flux / (flux @ flux)             # handle_flux as a covector

    z0 = ke.admissible_state(model, -1, -1, 0.125, rng, flux=4e-3,
                             link_scale=0.05, mouth=(0.02, 0.01))

    def seam_defect(start, d):
        fwd, bwd = ke.two_sided(model, start)
        return (float(np.abs(fwd[:, 0] - bwd[:, 0] @ d.T).max()
                      / np.linalg.norm(start)),
                model.handle_flux(fwd[:, 0]))

    base_defect, base_flux = seam_defect(z0, deck)
    rot_defect, rot_flux = seam_defect(quarter @ z0, deck_minus)
    wrong_defect, _ = seam_defect(quarter @ z0, deck)
    conj_defect, conj_flux = seam_defect(conj @ z0, deck)
    return {
        "quarter_preserves_energy": rel(model.energy(probe @ quarter.T)),
        "quarter_preserves_flux": float(np.abs(cycle @ quarter - cycle).max()),
        "quarter_flips_pin_sign": float(np.abs(
            quarter @ deck @ quarter.T - deck_minus).max()),
        "quarter_commutes_with_deck": float(np.abs(
            quarter @ deck @ quarter.T - deck).max()),
        "conjugation_preserves_energy": rel(model.energy(probe @ conj.T)),
        "conjugation_flips_flux": float(np.abs(cycle @ conj + cycle).max()),
        "conjugation_commutes_with_deck": float(np.abs(
            conj @ deck @ conj.T - deck).max()),
        "initial_flux": float(model.handle_flux(z0)[0]),
        "base_seam_defect": base_defect,
        "rotated_seam_defect_opposite_pin": rot_defect,
        "rotated_seam_defect_same_pin": wrong_defect,
        "rotated_flux_change": float(np.abs(rot_flux - base_flux).max()),
        "conjugated_seam_defect": conj_defect,
        "conjugated_flux_sum": float(np.abs(conj_flux + base_flux).max()),
        "flux_scale": float(np.abs(base_flux).max()),
    }


# ==========================================================================
# Verdict
# ==========================================================================

VERDICT_FIELDS = ("numerical_validity", "antipodal_map", "two_spheres",
                  "transport_on_fibre_modes", "pr129_cell",
                  "spinor_seam_consistency", "transport_seam_charges",
                  "flux_sign_vs_pin_label", "bam_cell", "support_dynamics",
                  "positive_response", "event_law")
CONFIRMED = {
    "antipodal_map": "FIBRE_TRANSLATION_FIXING_THE_HOPF_BASE",
    "two_spheres": "SPATIAL_BRANE_MOUTH_IS_NOT_THE_HOPF_BASE",
    "transport_on_fibre_modes":
        "ANTILINEAR_SQUARE_MINUS_ONE_ON_ODD_CHARGE_J_ON_THE_DOUBLET",
    "pr129_cell": "FIBRE_PRESERVING_NO_PAIRING",
    "spinor_seam_consistency":
        "TRANSPORT_SEAM_NEEDS_ANTIPERIODIC_DOUBLED_LOOP",
    "transport_seam_charges":
        "CHARGE_SURVIVES_SPIN_VANISHES_WITH_TIME_REVERSAL",
    "flux_sign_vs_pin_label": "INDEPENDENT_LABELS",
}
FIXED = {
    "bam_cell": "DEPENDS_ON_FIELD_TYPE_AND_TIME_EXTENSION",
    "support_dynamics": "EFFECTIVE_SUPPORT_ASSUMED",
    "positive_response": "NOT_DERIVED",
    "event_law": "NOT_DERIVED",
}


def verdict(checks: dict, groups: dict, counterexamples: dict | None = None
            ) -> dict:
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
