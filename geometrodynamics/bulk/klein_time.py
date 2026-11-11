"""Classical field on a time-nonorientable (Klein) quotient. Freeze ``bb6ff2c``.

Corrected after review of ``5b79975``: the frozen "no free field" reading
holds only in the frozen truncation, and two verdict labels were narrowed.

Background: Einstein static universe ``R x S^3`` (radius 1), conformally
coupled real scalar, zonal sector ``phi = u(t, chi) / sin(chi)``, for which

    u_tt = u_chichi,   u(t, 0) = u(t, pi) = 0

holds exactly. The antipode is ``chi -> pi - chi``.

The lattice ``chi_i = i pi / N`` with ``dt = dchi`` and the leapfrog rule
``u^{j+1} = L u^j - u^{j-1}`` is exact d'Alembert propagation; its discrete
frequencies are exactly ``omega_k = k``. A state is ``s_j = (u^j, u^{j+1})``.

Deck maps, with twist ``eta = +-1``:

    orientable control      (P u)^j   = F u^j
    Klein cap at m dt / 2   (T_m u)^j = eta F u^{m - j}

``F`` is the spatial flip. Everything here is classical and linear. See
``docs/klein_time_attachment_prereg.md`` for the frozen predictions and
``docs/klein_time_attachment.md`` for the result and its limits.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

N_FROZEN = 24
SEED = 2026111001
TOL_ZERO = 1e-9
TOL_NONDEGENERATE = 1e-6
RANK_RTOL = 1e-9

CAPS_P1 = (0, 1, 7, 24)
TABLE_P3 = {
    # m: (expected untwisted modes, expected twisted modes)
    48: (tuple(range(1, 24)), ()),
    24: (tuple(range(2, 24, 2)), tuple(range(1, 24, 2))),
    16: (tuple(range(3, 24, 3)), ()),
    12: (tuple(range(4, 24, 4)), tuple(range(2, 24, 4))),
    7: ((), ()),
}

VERDICT_FIELDS = ("phase_space", "complex_structure", "phase_closure",
                  "field_freedom", "continuum_free_field",
                  "irrational_spacing", "attachment_contribution",
                  "sector_pairing", "hbar",
                  "born_rule", "counting_function", "lepton_index")
# Two labels were narrowed after review of ``5b79975``; see the write-up.
CONFIRMED = {
    "phase_space": "ADMISSIBLE_FIELDS_ARE_A_LAGRANGIAN_HALF",
    "complex_structure": "TIME_DECK_IS_A_CONJUGATION_OF_THE_COVER_MODE_SPACE",
    "phase_closure": "TWO_CAPS_FORCE_PHASE_0_OR_PI_FIXED_BY_TWIST",
    "field_freedom": "UNIQUE_SOURCED_SOLUTION_IN_THE_FROZEN_TRUNCATION",
    "attachment_contribution":
        "QUADRATURE_SELECTION_NOT_THE_TIME_SYMMETRIC_RESPONSE",
}
FROZEN_LABELS = {
    "complex_structure": "TIME_DECK_IS_A_CONJUGATION_OF_THE_MODE_SPACE",
    "field_freedom": "GENERIC_CLOSED_HISTORY_HAS_NO_FREE_FIELD",
}
NOT_CLAIMED = {
    "irrational_spacing": "OPEN_SMALL_DIVISORS",
    "sector_pairing": "NOT_SUPPLIED",
    "hbar": "NOT_DERIVED",
    "born_rule": "NOT_DERIVED",
    "counting_function": "NOT_SELECTED",
    "lepton_index": "NOT_IDENTIFIED",
}
# Status vocabulary. A failed tolerance is not a mathematical refutation.
INCOMPLETE = "INCOMPLETE"
UNRESOLVED_NUMERICALLY = "UNRESOLVED_NUMERICALLY"
COUNTEREXAMPLE = "COUNTEREXAMPLE"

REFINEMENT_SIZES = (24, 48, 72, 96)
REFINEMENT_TAU = (7, 48)          # tau / pi = 7 / 48, held fixed
REFINEMENT_EXPECTED = {
    24: ((), ()),
    48: ((), (24,)),
    72: ((48,), (24,)),
    96: ((48,), (24, 72)),
}


# --------------------------------------------------------------------------
# Lattice operators
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Lattice:
    """Operators of the unit-Courant string with ``n - 1`` interior sites."""

    n: int

    @property
    def dim(self) -> int:
        return self.n - 1

    @property
    def neighbour_sum(self) -> np.ndarray:
        d = self.dim
        return np.eye(d, k=1) + np.eye(d, k=-1)

    @property
    def flip(self) -> np.ndarray:
        return np.eye(self.dim)[::-1].copy()

    @property
    def step(self) -> np.ndarray:
        """One leapfrog step on ``(u^j, u^{j+1})``."""
        d = self.dim
        eye, zero = np.eye(d), np.zeros((d, d))
        return np.block([[zero, eye], [-eye, self.neighbour_sum]])

    @property
    def wronskian(self) -> np.ndarray:
        """``Omega(s, s') = s^T W s'``, conserved by ``step``."""
        d = self.dim
        eye, zero = np.eye(d), np.zeros((d, d))
        return np.block([[zero, eye], [-eye, zero]])

    @property
    def swap(self) -> np.ndarray:
        d = self.dim
        eye, zero = np.eye(d), np.zeros((d, d))
        return np.block([[zero, eye], [eye, zero]])

    @property
    def flip2(self) -> np.ndarray:
        d = self.dim
        zero = np.zeros((d, d))
        return np.block([[self.flip, zero], [zero, self.flip]])

    def mode(self, k: int) -> np.ndarray:
        """Unit-norm spatial profile ``sin(k chi_i)`` on the interior sites."""
        chi = np.arange(1, self.n) * np.pi / self.n
        v = np.sin(k * chi)
        return v / np.linalg.norm(v)


def klein_deck(lat: Lattice, m: int, eta: int = 1) -> np.ndarray:
    """State-space matrix of ``(T u)^j = eta F u^{m-j}``.

    ``(T s)_0 = eta (F u^m, F u^{m-1}) = eta F S s_{m-1}``.
    """
    power = np.linalg.matrix_power(lat.step, m - 1) if m >= 1 else \
        np.linalg.inv(lat.step)
    return eta * lat.flip2 @ lat.swap @ power


def orientable_deck(lat: Lattice) -> np.ndarray:
    """State-space matrix of the time-preserving antipode."""
    return lat.flip2


def march(lat: Lattice, state: np.ndarray, steps: int) -> np.ndarray:
    """Leapfrog by explicit loop. Returns levels ``u^0 .. u^steps``.

    Deliberately independent of ``Lattice.step`` and of matrix powers.
    """
    d = lat.dim
    levels = np.zeros((max(steps, 1) + 1, d))
    levels[0], levels[1] = state[:d], state[d:]
    for j in range(1, steps):
        padded = np.concatenate(([0.0], levels[j], [0.0]))
        levels[j + 1] = padded[2:] + padded[:-2] - levels[j - 1]
    return levels[:steps + 1]


# --------------------------------------------------------------------------
# Linear algebra helpers
# --------------------------------------------------------------------------

def null_space(a: np.ndarray, rtol: float = RANK_RTOL) -> np.ndarray:
    """Orthonormal basis of ``ker a`` as columns."""
    _, s, vt = np.linalg.svd(a)
    scale = s[0] if s.size and s[0] > 0 else 1.0
    rank = int(np.sum(s > rtol * scale))
    return vt[rank:].T.copy()


def rank(a: np.ndarray, rtol: float = RANK_RTOL) -> int:
    if a.size == 0:
        return 0
    s = np.linalg.svd(a, compute_uv=False)
    return int(np.sum(s > rtol * max(s[0], 1e-300))) if s[0] > 0 else 0


def _smallest_singular(a: np.ndarray) -> float:
    return float(np.linalg.svd(a, compute_uv=False)[-1]) if a.size else 0.0


# --------------------------------------------------------------------------
# P1  Lagrangian halving
# --------------------------------------------------------------------------

def halving_report(lat: Lattice, m: int, eta: int) -> dict:
    t = klein_deck(lat, m, eta)
    w = lat.wronskian
    size = 2 * lat.dim
    plus = null_space(t - np.eye(size))
    minus = null_space(t + np.eye(size))
    return {
        "cap": m, "eta": eta,
        "involution_error": float(np.abs(t @ t - np.eye(size)).max()),
        "antisymplectic_error": float(np.abs(t.T @ w @ t + w).max()),
        "dim_plus": plus.shape[1], "dim_minus": minus.shape[1],
        "omega_on_plus": float(np.abs(plus.T @ w @ plus).max()),
        "omega_on_minus": float(np.abs(minus.T @ w @ minus).max()),
        "pairing_min_singular": _smallest_singular(plus.T @ w @ minus),
    }


def orientable_control(lat: Lattice) -> dict:
    p = orientable_deck(lat)
    w = lat.wronskian
    size = 2 * lat.dim
    fixed = null_space(p - np.eye(size))
    restricted = fixed.T @ w @ fixed
    return {
        "symplectic_error": float(np.abs(p.T @ w @ p - w).max()),
        "dim_fixed": fixed.shape[1],
        "restricted_rank": rank(restricted),
        "restricted_min_singular": _smallest_singular(restricted),
    }


def direct_admissibility(lat: Lattice, m: int, eta: int,
                         rng: np.random.Generator) -> dict:
    """March by explicit loop and test ``u^j = eta F u^{m-j}`` level by level."""
    size = 2 * lat.dim
    t = klein_deck(lat, m, eta)
    plus = null_space(t - np.eye(size))
    inside = plus @ rng.standard_normal(plus.shape[1])
    generic = rng.standard_normal(size)

    def defect(state):
        span = max(m, 2)
        levels = march(lat, state, span)
        worst = 0.0
        for j in range(0, m + 1):
            worst = max(worst, float(np.abs(
                levels[j] - eta * levels[m - j][::-1]).max()))
        return worst / float(np.abs(levels).max())

    return {"cap": m, "eta": eta, "admissible_defect": defect(inside),
            "generic_defect": defect(generic)} if m >= 2 else \
        {"cap": m, "eta": eta, "admissible_defect": 0.0,
         "generic_defect": float("nan"), "skipped": "cap inside first step"}


# --------------------------------------------------------------------------
# P2  Complex structure
# --------------------------------------------------------------------------

def complex_structure(lat: Lattice) -> np.ndarray:
    """Quarter turn of every mode plane in the sense of time evolution."""
    values, vectors = np.linalg.eig(lat.step)
    if np.min(np.abs(values.imag)) < 1e-9:
        raise ValueError("step matrix has a real eigenvalue; j is undefined")
    j = vectors @ np.diag(1j * np.sign(values.imag)) @ np.linalg.inv(vectors)
    if np.abs(j.imag).max() > 1e-8:
        raise ValueError("complex structure is not real")
    return j.real


def conjugation_report(lat: Lattice, m: int, eta: int) -> dict:
    j = complex_structure(lat)
    t = klein_deck(lat, m, eta)
    w = lat.wronskian
    size = 2 * lat.dim
    plus = null_space(t - np.eye(size))
    g = w @ j
    eig = np.linalg.eigvalsh((g + g.T) / 2)
    sign = 1.0 if eig[0] > 0 else -1.0
    image = j @ plus
    return {
        "cap": m, "eta": eta,
        "j_squared_error": float(np.abs(j @ j + np.eye(size)).max()),
        "commutes_with_evolution": float(np.abs(j @ lat.step
                                                - lat.step @ j).max()),
        "antilinear_error": float(np.abs(t @ j + j @ t).max()),
        "j_plus_in_minus": float(np.abs((t + np.eye(size)) @ image).max()),
        "metric_symmetry_error": float(np.abs(g - g.T).max()),
        "metric_definite_min": float(np.min(sign * eig)),
        "hermitian_imaginary_on_plus": float(np.abs(plus.T @ w @ plus).max()),
    }


# --------------------------------------------------------------------------
# P3  Two caps
# --------------------------------------------------------------------------

def two_cap_space(lat: Lattice, m: int, eta0: int, eta1: int) -> np.ndarray:
    """Fields admissible at a cap at 0 and a cap at ``m dt / 2``."""
    size = 2 * lat.dim
    eye = np.eye(size)
    stacked = np.vstack([klein_deck(lat, 0, eta0) - eye,
                         klein_deck(lat, m, eta1) - eye])
    return null_space(stacked)


def surviving_modes(lat: Lattice, basis: np.ndarray) -> tuple:
    """Harmonics ``k`` whose mode plane meets the span of ``basis``."""
    d = lat.dim
    found = []
    for k in range(1, lat.n):
        v = lat.mode(k)
        plane = np.zeros((2 * d, 2))
        plane[:d, 0], plane[d:, 1] = v, v
        weight = float(np.linalg.norm(plane.T @ basis)) if basis.size else 0.0
        if weight > 0.5:
            found.append(k)
    return tuple(found)


def predicted_modes(n: int, m: int, twisted: bool) -> tuple:
    """``2 k tau = k m pi / n`` must be ``0`` or ``pi`` mod ``2 pi``."""
    target = 1 if twisted else 0
    return tuple(k for k in range(1, n) if (k * m) % n == 0
                 and ((k * m) // n) % 2 == target)


def two_cap_report(lat: Lattice, m: int) -> dict:
    out = {"m": m, "tau_over_pi": m / (2 * lat.n)}
    for name, (e0, e1) in {"pp": (1, 1), "pm": (1, -1),
                           "mp": (-1, 1), "mm": (-1, -1)}.items():
        basis = two_cap_space(lat, m, e0, e1)
        out[name] = {"dim": basis.shape[1],
                     "modes": surviving_modes(lat, basis)}
    out["predicted_untwisted"] = predicted_modes(lat.n, m, False)
    out["predicted_twisted"] = predicted_modes(lat.n, m, True)
    return out


def loop_phase(k: int, m: int, n: int) -> float:
    """Phase ``2 k tau`` accumulated around the time loop, mod ``2 pi``."""
    return float((k * m * np.pi / n) % (2 * np.pi))


# --------------------------------------------------------------------------
# P4  Closed history as a boundary-value problem
# --------------------------------------------------------------------------

def periodic_operator(lat: Lattice, m: int) -> np.ndarray:
    """``D = (shift_+ + shift_-) (x) 1 - 1 (x) L`` on ``m`` time levels."""
    cyc = np.roll(np.eye(m), 1, axis=0)
    return np.kron(cyc + cyc.T, np.eye(lat.dim)) - np.kron(
        np.eye(m), lat.neighbour_sum)


def periodic_deck(lat: Lattice, m: int) -> np.ndarray:
    """``(T_0 u)^j = F u^{-j mod m}`` on the periodic time lattice."""
    reverse = np.zeros((m, m))
    for j in range(m):
        reverse[j, (-j) % m] = 1.0
    return np.kron(reverse, lat.flip)


def kernel_report(lat: Lattice, m: int) -> dict:
    d_op = periodic_operator(lat, m)
    t0 = periodic_deck(lat, m)
    kernel = null_space(d_op)
    invariant = 0
    if kernel.shape[1]:
        eig = np.linalg.eigvalsh((kernel.T @ t0 @ kernel
                                  + kernel.T @ t0.T @ kernel) / 2)
        invariant = int(np.sum(eig > 0.5))
    return {"m": m, "nullity": kernel.shape[1], "invariant_nullity": invariant,
            "deck_commutator": float(np.abs(t0 @ d_op - d_op @ t0).max()),
            "operator_symmetry": float(np.abs(d_op - d_op.T).max())}


def slaved_field_report(lat: Lattice, m: int = 7) -> dict:
    """Generic closed history: the source fixes the field."""
    d_op = periodic_operator(lat, m)
    t0 = periodic_deck(lat, m)
    green = np.linalg.inv(d_op)
    d = lat.dim

    def levels_of(source):
        return (green @ source).reshape(m, d)

    def ratio(field, level):
        before = np.linalg.norm([field[(level - s) % m] for s in (1, 2, 3)])
        after = np.linalg.norm([field[(level + s) % m] for s in (1, 2, 3)])
        return float(before / after)

    # Invariant source: one time level, on the cap, flip-even in space.
    profile = np.zeros(d)
    profile[4], profile[d - 1 - 4] = 1.0, 1.0
    on_cap = np.zeros((m, d))
    on_cap[0] = profile
    field_cap = levels_of(on_cap.ravel())

    # Supplementary: one level off the cap, one site (not deck invariant).
    off_cap = np.zeros((m, d))
    off_cap[2, 6] = 1.0
    field_off = levels_of(off_cap.ravel())

    return {
        "m": m,
        "smallest_singular_value": _smallest_singular(d_op),
        "green_symmetry_error": float(np.abs(green - green.T).max()),
        "green_deck_commutator": float(np.abs(t0 @ green - green @ t0).max()),
        "field_invariance_error": float(np.abs(
            t0 @ field_cap.ravel() - field_cap.ravel()).max()),
        "before_after_ratio_on_cap": ratio(field_cap, 0),
        "before_after_ratio_off_cap": ratio(field_off, 2),
        "field_before_source_norm": float(np.linalg.norm(
            [field_off[(2 - s) % m] for s in (1, 2, 3)])),
    }


def fredholm_report(lat: Lattice, m: int, rng: np.random.Generator) -> dict:
    """Resonant closed history: solvable iff the source avoids ``ker D``."""
    d_op = periodic_operator(lat, m)
    kernel = null_space(d_op)
    source = rng.standard_normal(d_op.shape[0])
    in_kernel = kernel @ (kernel.T @ source)
    clean = source - in_kernel

    def residual(s):
        x = np.linalg.lstsq(d_op, s, rcond=None)[0]
        return float(np.linalg.norm(d_op @ x - s))

    return {"m": m, "nullity": kernel.shape[1],
            "residual_orthogonal_source": residual(clean),
            "residual_generic_source": residual(source),
            "kernel_component_norm": float(np.linalg.norm(in_kernel))}


# --------------------------------------------------------------------------
# Review round: fixed-geometry refinement and the continuum
# --------------------------------------------------------------------------

def continuum_modes(p: int, q: int, k_max: int, twisted: bool) -> tuple:
    """Continuum survivors below ``k_max`` for ``tau / pi = p / q``.

    ``2 k tau = 2 pi k p / q`` must be ``0`` (untwisted) or ``pi`` (twisted)
    mod ``2 pi``. In lowest terms: untwisted iff ``q | k``; twisted iff ``q``
    is even and ``k`` is an odd multiple of ``q / 2``.
    """
    g = int(np.gcd(p, q))
    p, q = p // g, q // g
    if not twisted:
        return tuple(k for k in range(1, k_max) if k % q == 0)
    if q % 2:
        return ()
    half = q // 2
    return tuple(k for k in range(1, k_max)
                 if k % half == 0 and (k // half) % 2 == 1)


def refinement_report(sizes=REFINEMENT_SIZES, tau=REFINEMENT_TAU) -> list:
    """Hold the physical cap separation fixed and refine the lattice."""
    p, q = tau
    rows = []
    for n in sizes:
        if (2 * n * p) % q:
            raise ValueError(f"tau = {p} pi / {q} is not on the N = {n} lattice")
        m = 2 * n * p // q
        lat = Lattice(n)
        rows.append({
            "n": n, "m": m,
            "untwisted": surviving_modes(lat, two_cap_space(lat, m, 1, 1)),
            "twisted": surviving_modes(lat, two_cap_space(lat, m, 1, -1)),
            "continuum_untwisted": continuum_modes(p, q, n, False),
            "continuum_twisted": continuum_modes(p, q, n, True),
        })
    return rows


def continuum_counterexample(p: int, q: int, twisted: bool = False) -> dict:
    """An explicit homogeneous field at ``tau / pi = p / q``.

    The continuum mode ``sin(k chi) cos(k t - alpha)`` with the smallest
    surviving ``k``, sampled on the smallest lattice that resolves it, where
    the leapfrog is exact. Both cap conditions are then tested by marching.
    """
    g = int(np.gcd(p, q))
    p, q = p // g, q // g
    if twisted and q % 2:
        return {"p": p, "q": q, "twisted": True, "exists": False}
    k = q // 2 if twisted else q
    n = 2 * q                      # k < n, and m = 4 p is an integer
    m = 2 * n * p // q
    lat = Lattice(n)
    chi = np.arange(1, n) * np.pi / n
    dt = np.pi / n
    even = (k + 1) % 2 == 0        # spatial parity (-1)^(k+1)
    time = (lambda t: np.cos(k * t)) if even else (lambda t: np.sin(k * t))
    levels = march(lat, np.concatenate([np.sin(k * chi) * time(0.0),
                                        np.sin(k * chi) * time(dt)]),
                   max(m, 2))
    scale = float(np.abs(levels).max())
    eta1 = -1 if twisted else 1

    def defect(cap, eta):
        # (T u)^j = eta F u^{cap - j}; for cap = 0 use the exact period 2n.
        full = march(lat, np.concatenate([levels[0], levels[1]]), 2 * n)
        return max(float(np.abs(full[j] - eta * full[(cap - j) % (2 * n)][::-1]
                                ).max()) for j in range(0, 2 * n)) / scale

    return {"p": p, "q": q, "twisted": twisted, "exists": True, "k": k,
            "n": n, "m": m, "amplitude": scale,
            "quadrature": "cos" if even else "sin",
            "cap_0_defect": defect(0, 1), "cap_tau_defect": defect(m, eta1)}


def continuum_gap(alpha: float, k_max: int = 100_000) -> dict:
    """Diagnostic only. Smallest ``|k^2 - (q / alpha)^2|`` for ``k <= k_max``.

    ``alpha = tau / pi``. The untwisted continuum operator on a closed
    history of period ``2 tau`` is diagonal with symbol ``nu_q^2 - k^2``,
    ``nu_q = q / alpha``. A gap bounded away from zero means a bounded
    inverse on ``L^2``; it is not an existence theory for the lattice limit.
    """
    k = np.arange(1, k_max + 1, dtype=float)
    # The objective is a difference of squares, so the nearest integer to
    # k alpha is not always the minimiser. Try floor and ceiling, q = 0
    # included (corrected after review of c433f91).
    lower = np.floor(k * alpha)
    candidates = np.stack([lower, lower + 1.0])
    gaps = np.abs(k ** 2 - (candidates / alpha) ** 2)
    pick = np.argmin(gaps, axis=0)
    gap = gaps[pick, np.arange(k.size)]
    q = candidates[pick, np.arange(k.size)]
    worst = int(np.argmin(gap))
    return {"alpha": float(alpha), "k_max": k_max,
            "min_gap": float(gap[worst]), "at_k": int(k[worst]),
            "at_q": int(q[worst]),
            "min_k_times_distance": float(np.min(
                k * np.abs(k * alpha - np.rint(k * alpha))))}


# --------------------------------------------------------------------------
# Review round: what the attachment adds to temporal closure
# --------------------------------------------------------------------------

def antiperiodic_operator(lat: Lattice, m: int) -> np.ndarray:
    """Time-orientable control with ``u^{j+m} = -u^j``."""
    cyc = np.roll(np.eye(m), 1, axis=0)
    cyc[0, m - 1] = -1.0           # the wrap-around link carries the sign
    return np.kron(cyc + cyc.T, np.eye(lat.dim)) - np.kron(
        np.eye(m), lat.neighbour_sum)


def antiperiodic_deck(lat: Lattice, m: int) -> np.ndarray:
    """``(T_0 u)^j = F u^{-j}`` with ``u^{-j} = -u^{m-j}`` for ``j != 0``."""
    reverse = np.zeros((m, m))
    reverse[0, 0] = 1.0
    for j in range(1, m):
        reverse[j, m - j] = -1.0
    return np.kron(reverse, lat.flip)


def _invariant_count(kernel: np.ndarray, deck: np.ndarray) -> int:
    if not kernel.shape[1]:
        return 0
    block = kernel.T @ deck @ kernel
    return int(np.sum(np.linalg.eigvalsh((block + block.T) / 2) > 0.5))


def closure_control_report(lat: Lattice, m: int) -> dict:
    """Same period, same operator: orientable closed time against the caps.

    The orientable control imposes only ``u^{j+m} = +-u^j``. The Klein
    history adds invariance under the cap at 0.
    """
    periodic = null_space(periodic_operator(lat, m))
    anti_op = antiperiodic_operator(lat, m)
    anti = null_space(anti_op)
    anti_deck = antiperiodic_deck(lat, m)
    return {
        "m": m,
        "orientable_periodic_modes": surviving_time_modes(lat, periodic, m),
        "orientable_periodic_dim": periodic.shape[1],
        "klein_untwisted_dim": _invariant_count(periodic,
                                                periodic_deck(lat, m)),
        "orientable_antiperiodic_modes": surviving_time_modes(lat, anti, m),
        "orientable_antiperiodic_dim": anti.shape[1],
        "klein_twisted_dim": _invariant_count(anti, anti_deck),
        "antiperiodic_deck_commutator": float(np.abs(
            anti_deck @ anti_op - anti_op @ anti_deck).max()),
        "antiperiodic_deck_involution": float(np.abs(
            anti_deck @ anti_deck - np.eye(anti_deck.shape[0])).max()),
    }


def surviving_time_modes(lat: Lattice, basis: np.ndarray, m: int) -> tuple:
    """Spatial harmonics present in a space of fields on ``m`` time levels."""
    if not basis.shape[1]:
        return ()
    stacked = basis.reshape(m, lat.dim, -1)
    return tuple(k for k in range(1, lat.n)
                 if float(np.linalg.norm(np.einsum("i,jic->jc", lat.mode(k),
                                                   stacked))) > 0.5)


def response_control_report(lat: Lattice, m_generic: int = 7,
                            m_resonant: int = 24) -> dict:
    """Which features of the sourced problem need the attachment at all."""
    d = lat.dim
    # Generic period: the unrestricted orientable problem is already unique.
    d_op = periodic_operator(lat, m_generic)
    green = np.linalg.inv(d_op)
    t0 = periodic_deck(lat, m_generic)
    source = np.zeros((m_generic, d))
    source[2, 6] = 1.0             # one site, off the cap, not deck invariant
    field = (green @ source.ravel()).reshape(m_generic, d)
    mirror = max(float(np.abs(field[(2 - s) % m_generic]
                              - field[(2 + s) % m_generic]).max())
                 for s in (1, 2, 3))
    size = d_op.shape[0]
    invariant_sources = rank(t0 + np.eye(size))

    # Resonant period: count the solvability conditions.
    r_op = periodic_operator(lat, m_resonant)
    r_deck = periodic_deck(lat, m_resonant)
    kernel = null_space(r_op)
    r_size = r_op.shape[0]
    projector = (np.eye(r_size) + r_deck) / 2
    return {
        "m_generic": m_generic, "m_resonant": m_resonant,
        "source_is_deck_invariant": bool(np.abs(
            t0 @ source.ravel() - source.ravel()).max() < 1e-12),
        "unrestricted_smallest_singular": _smallest_singular(d_op),
        "unrestricted_green_symmetry": float(np.abs(green - green.T).max()),
        "time_mirror_defect_without_deck": mirror,
        "field_scale": float(np.abs(field).max()),
        "source_space_dim": size,
        "deck_invariant_source_dim": invariant_sources,
        "conditions_unrestricted": rank(kernel.T),
        "conditions_on_invariant_sources": rank(kernel.T @ projector),
    }


def sector_action_report(lat: Lattice, m: int, eta: int = 1) -> dict:
    """Does the cover's complex structure act inside one twist sector?

    ``j L_+ = L_-``, so the answer is no: ``L_+`` meets ``j L_+`` only in 0.
    """
    size = 2 * lat.dim
    plus = null_space(klein_deck(lat, m, eta) - np.eye(size))
    j = complex_structure(lat)
    dim = plus.shape[1]
    return {
        "cap": m, "eta": eta, "dim": dim,
        "dim_plus_meet_j_plus": 2 * dim - rank(np.hstack([plus, j @ plus])),
    }


# --------------------------------------------------------------------------
# Verdict
# --------------------------------------------------------------------------

def verdict(checks: dict, groups: dict, counterexamples: dict | None = None,
            extra: dict | None = None) -> dict:
    """Status per field.

    ``COUNTEREXAMPLE``           an explicit counterexample is recorded
    ``INCOMPLETE``               a required check was not run
    ``UNRESOLVED_NUMERICALLY``   a check ran and missed its tolerance
    confirmed label              every check of the group is present and true

    A missed tolerance is a numerical outcome, not a refutation.
    """
    counterexamples = counterexamples or {}
    out, failed, missing_all = {}, [], []
    for field, names in groups.items():
        missing = [n for n in names if n not in checks]
        bad = [n for n in names if n in checks and not checks[n]]
        if field in counterexamples:
            out[field] = COUNTEREXAMPLE
        elif missing or not names:
            out[field] = INCOMPLETE
        elif bad:
            out[field] = UNRESOLVED_NUMERICALLY
        else:
            out[field] = CONFIRMED.get(field, "CONFIRMED")
        failed += bad
        missing_all += missing
    out.update(NOT_CLAIMED)
    out.update(extra or {})
    out["failed_checks"] = failed
    out["missing_checks"] = missing_all
    out["counterexamples"] = dict(counterexamples)
    return out
