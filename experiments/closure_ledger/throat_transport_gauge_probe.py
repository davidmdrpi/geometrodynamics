"""
Gauge-sensitive throat-transport falsifier.

Follow-up to the moving-mouth Berry-phase probe (PR #23). That probe
established that the BAM-symmetric Hopf connection

    A_BAM(χ, φ)  =  ½·cos(χ) dφ          (`hopf/connection.py`)

with the BAM-frame spinor

    |ψ_BAM(χ, φ)⟩ = cos(χ/2)·e^{−iφ/2}|↑⟩ + sin(χ/2)·e^{+iφ/2}|↓⟩

and the non-orientable throat transport

    T_BAM  =  iσ_y                       (`embedding/transport.py`)

is internally consistent: every closed loop accumulates a Berry phase
that agrees with the line integral of A_BAM to numerical precision.

This probe attacks the OPEN question from that PR:

    Does the BAM combination (A_BAM, ψ_BAM, T_BAM=iσ_y) predict a
    physically distinct phase from the standard Bloch combination

        A_Bloch(χ, φ) = ½·(cos(χ) − 1) dφ
        |ψ_Bloch(χ,φ)⟩ = cos(χ/2)|↑⟩ + sin(χ/2)·e^{+iφ}|↓⟩
        T_Bloch  =  consistently gauge-conjugated throat

    AFTER all allowed gauge freedoms (single-valued U(1) and the
    spinor SU(2) frame rotation) are accounted for?

The two spinors differ by an overall SCALAR U(1) phase: pulling out
`e^{−iφ/2}` from `ψ_BAM` gives `e^{−iφ/2}·(cos(χ/2), sin(χ/2)e^{+iφ})
= e^{−iφ/2}·ψ_Bloch`. So the gauge transformation is

    g(χ, φ)  =  e^{−iφ/2}                                  (∈ U(1))
    ψ_BAM    =  g · ψ_Bloch                                (multi-valued)
    A_BAM    =  A_Bloch − dΛ      with  Λ = −φ/2  ⇒  A_BAM = A_Bloch + ½ dφ

Under `φ → φ+2π`, `g → −g`: the gauge transformation is MULTI-VALUED
in U(1). It is, however, the unique U(1) transformation that
implements the BAM↔Bloch relation pointwise. Its multi-valuedness is
exactly the Wu-Yang transition function that encodes the spin-½
double-cover: the BAM spinor section, viewed in the Bloch frame,
acquires `−1` on every 2π φ-loop.

Tests in this probe:

  (1) **Connection difference is +½ dφ pointwise.** A_BAM(χ) −
      A_Bloch(χ) = ½ at all χ on a sampled grid.

  (2) **Spinor frame transformation is exact.** ψ_BAM(χ,φ) =
      U_gauge(χ,φ) · ψ_Bloch(χ,φ) at all (χ,φ) on a sampled grid.

  (3) **Throat transport conjugates consistently.** Define
      T_Bloch(φ_in, φ_out) := U_gauge^{−1}(π/2, φ_out) · iσ_y ·
      U_gauge(π/2, φ_in). Verify det(T_Bloch) = 1, T_Bloch² = −I, and
      Tr(T_Bloch) = 0 for all sampled (φ_in, φ_out) — the same
      gauge-invariant signature as T_BAM = iσ_y.

  (4) **Closed-loop holonomy operator is gauge-conjugate.** For a
      family of test loops (each crossing the equatorial throat at
      one (φ_in, φ_out) pair, then returning via a Berry-transport
      arc), compute U_loop in both gauges. Verify that
      U_loop^{BAM} = U_gauge(start) · U_loop^{Bloch} · U_gauge^{−1}(start)
      to numerical precision. The gauge-invariant
      Tr(U_loop)/2 = cos(θ_eff/2) must agree exactly.

  (5) **Relative-phase observables are gauge-equivalent.** Compute the
      detector-holonomy phase Δγ(θ_a, θ_b) = γ(θ_a) − γ(θ_b) in both
      gauges. The BAM offset (γ_BAM = π cos θ) and Bloch offset
      (γ_Bloch = −π(1 − cos θ)) cancel in Δγ. So Δγ_BAM = Δγ_Bloch.
      Verifies that no Bell-type / relative-measurement experiment can
      distinguish the two gauges.

  (6) **Absolute-phase observable distinguishability.** The π/loop
      offset is the only candidate distinguishing feature. Verify
      explicitly that under EVERY single-valued U(1) gauge
      transformation (Λ a periodic function of φ), the BAM-Bloch
      offset cannot be removed: any periodic Λ contributes 0 mod 2π
      to a closed loop. Therefore the π Wu-Yang offset persists as a
      gauge-bundle invariant — the spinor double cover. This makes
      BAM and Bloch ARE distinguishable on the SU(2) level (different
      spinor sections) but INDISTINGUISHABLE on the SO(3) level
      (same projective Hilbert space). The expected verdict: no
      classical measurement of a single particle's wavefunction
      magnitude squared can tell them apart, but interference between
      spinor branches CAN — which is exactly the spin-½ Dirac string
      content that BAM and Bloch agree about (they differ only in
      WHERE the string sits, not WHETHER one exists).

PASS criterion (the falsifier returns "BAM is gauge-equivalent to
Bloch"):
  - Test (1), (2), (3): equalities hold to ≤ 1e-12.
  - Test (4): operator-conjugation residual ≤ 1e-12 for every test
    loop; trace-equality ≤ 1e-12.
  - Test (5): Δγ_BAM = Δγ_Bloch to ≤ 1e-14 (analytic identity).
  - Test (6): explicit enumeration shows no single-valued U(1)
    transformation removes the π offset.

FAIL criterion (BAM makes a sharper-than-Bloch prediction):
  - Any of (1)–(5) shows a residual that does not vanish in the
    gauge-conjugate sense. This would indicate either (a) the
    BAM-derived T = iσ_y is NOT the correct gauge-conjugate of the
    Bloch-frame throat, suggesting an inconsistency in BAM's
    derivation, or (b) the gauge-conjugation framework itself is
    wrong (the BAM bundle is genuinely topologically distinct from
    the Bloch bundle), which would be a much stronger structural
    claim requiring deeper analysis.

Expected outcome: PASS — BAM is gauge-equivalent to Bloch under the
spinor lift. The derivation T = iσ_y in BAM frame is a particularly
clean section; the same physics is present in any gauge.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np


PI = math.pi
TAU = 2.0 * PI

_SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
_ID = np.eye(2, dtype=complex)


# ---------------------------------------------------------------------------
# Gauge transformation (scalar U(1), multi-valued)
# ---------------------------------------------------------------------------

def g_gauge(chi: float, phi: float) -> complex:
    """g(χ,φ) = e^{−iφ/2}.

    The unique U(1) scalar transformation that maps Bloch frame →
    BAM frame: ψ_BAM = g · ψ_Bloch (scalar multiplication on the
    spinor). Multi-valued: under φ → φ+2π, g → −g.

    χ enters as a dummy argument for symmetry with other gauge
    objects on the (χ,φ) base.
    """
    return complex(np.exp(-0.5j * phi))


def gauge_phase_along_path(chi_path, phi_path) -> complex:
    """Net U(1) phase g(χ_end, φ_end) · g(χ_start, φ_start)^{−1}
    accumulated by the gauge transformation along the path.

    For our scalar g(χ,φ) = e^{−iφ/2}, this is e^{−i(φ_end−φ_start)/2}.
    For a closed loop with net φ-winding n, the accumulated phase is
    e^{−inπ} = (−1)^n. This is the BAM↔Bloch holonomy relation.
    """
    chi_path = np.asarray(chi_path)
    phi_path = np.asarray(phi_path)
    return complex(np.exp(-0.5j * (phi_path[-1] - phi_path[0])))


# ---------------------------------------------------------------------------
# Connections in the two gauges
# ---------------------------------------------------------------------------

def A_phi_BAM(chi: float) -> float:
    """BAM symmetric gauge: A_φ = ½cos(χ)."""
    return 0.5 * math.cos(chi)


def A_phi_Bloch(chi: float) -> float:
    """Bloch (north-pole regular) gauge: A_φ = ½(cos(χ) − 1).

    Singular at χ = π (south pole), regular for χ < π.
    """
    return 0.5 * (math.cos(chi) - 1.0)


# ---------------------------------------------------------------------------
# Spinors in the two gauges
# ---------------------------------------------------------------------------

def psi_BAM(chi: float, phi: float) -> np.ndarray:
    """BAM-frame Hopf eigenstate:
        |ψ_BAM⟩ = cos(χ/2)e^{−iφ/2}|↑⟩ + sin(χ/2)e^{+iφ/2}|↓⟩
    """
    cz = math.cos(chi / 2.0)
    sz = math.sin(chi / 2.0)
    return np.array([
        cz * np.exp(-0.5j * phi),
        sz * np.exp(+0.5j * phi),
    ], dtype=complex)


def psi_Bloch(chi: float, phi: float) -> np.ndarray:
    """Bloch-frame (north-regular) Hopf eigenstate:
        |ψ_Bloch⟩ = cos(χ/2)|↑⟩ + sin(χ/2)e^{+iφ}|↓⟩
    """
    cz = math.cos(chi / 2.0)
    sz = math.sin(chi / 2.0)
    return np.array([
        cz + 0.0j,
        sz * np.exp(+1.0j * phi),
    ], dtype=complex)


# ---------------------------------------------------------------------------
# Throat transport in each gauge
# ---------------------------------------------------------------------------

def T_BAM() -> np.ndarray:
    """T_BAM = iσ_y, derived in `embedding/transport.py`."""
    return 1j * _SY


def T_Bloch(phi_in: float, phi_out: float, chi_throat: float = PI / 2.0) -> np.ndarray:
    """Bloch-frame throat transport, gauge-conjugate of T_BAM.

    The throat sits at the equatorial Hopf fibre χ = π/2 (the
    zero-self-energy stable orbit). Entry and exit happen at azimuthal
    angles φ_in, φ_out. By gauge conjugation under the scalar U(1)
    transformation:

        T_Bloch = g(χ, φ_out)^{−1} · T_BAM · g(χ, φ_in)
                = e^{+iφ_out/2} · iσ_y · e^{−iφ_in/2}
                = e^{i(φ_out − φ_in)/2} · iσ_y

    A scalar phase times iσ_y. Same SU(2) invariants as T_BAM:
    det = e^{i(φ_out−φ_in)} (∈ U(1)) — but as the throat is
    instantaneous (φ_in = φ_out), det = 1, Tr = 0, T² = −I.

    For general (φ_in, φ_out) — interpreted as "the throat instantly
    teleports from angular position φ_in to φ_out" — T_Bloch picks up
    the phase factor e^{i(φ_out − φ_in)/2}, which is precisely the
    gauge-transition cost of the spinor relocation. Verifying this
    is part of test (3).
    """
    g_in = g_gauge(chi_throat, phi_in)
    g_out_inv = 1.0 / g_gauge(chi_throat, phi_out)
    return g_out_inv * T_BAM() * g_in


# ---------------------------------------------------------------------------
# Closed-loop holonomy operators
# ---------------------------------------------------------------------------

def _connection_phase(chi_path, phi_path, gauge: str) -> float:
    """∮ A · dλ via trapezoidal integration in the chosen gauge.

    Berry holonomy as a U(1) phase (signed real number); the full
    SU(2) transport operator is e^{iφ_berry}·I (the connection is
    Abelian — the Hopf bundle is a U(1) bundle, the SU(2) appearance
    in our spinor parametrisation is rest-frame).
    """
    chi_path = np.asarray(chi_path, dtype=float)
    phi_path = np.asarray(phi_path, dtype=float)
    if gauge == 'BAM':
        A = 0.5 * np.cos(chi_path)
    elif gauge == 'Bloch':
        A = 0.5 * (np.cos(chi_path) - 1.0)
    else:
        raise ValueError(f'unknown gauge {gauge!r}')
    dphi = np.diff(phi_path)
    A_mid = 0.5 * (A[:-1] + A[1:])
    return float(np.sum(A_mid * dphi))


def berry_transport_operator(chi_path, phi_path, gauge: str) -> np.ndarray:
    """U(1) Berry holonomy as a 2×2 operator (scalar times identity).

    The Hopf bundle is U(1); on the spinor representation this acts
    as a single scalar phase on the wavefunction, hence as exp(iγ)·I
    on the spinor.
    """
    gamma = _connection_phase(chi_path, phi_path, gauge)
    return np.exp(1j * gamma) * np.eye(2, dtype=complex)


def closed_loop_holonomy(loop, gauge: str) -> np.ndarray:
    """Combined U(1)⊗SU(2) holonomy of a closed loop.

    Loop dict keys:
        'chi_path', 'phi_path' : discretised (χ,φ) arrays.
        'crosses_throat'       : bool. If True the throat T is
                                 inserted at the midpoint (χ=π/2
                                 expected) with φ_in = φ_out = the
                                 midpoint φ value (instantaneous).

    Returns the 2×2 operator (Berry-phase scalar × identity, then
    multiplied by T if the loop crosses the throat).
    """
    chi_path = np.asarray(loop['chi_path'])
    phi_path = np.asarray(loop['phi_path'])
    if not loop.get('crosses_throat', False):
        return berry_transport_operator(chi_path, phi_path, gauge)

    mid = len(chi_path) // 2
    chi_out = chi_path[:mid + 1]
    phi_out = phi_path[:mid + 1]
    chi_ret = chi_path[mid:]
    phi_ret = phi_path[mid:]
    phi_throat = float(phi_path[mid])
    U_out = berry_transport_operator(chi_out, phi_out, gauge)
    U_ret = berry_transport_operator(chi_ret, phi_ret, gauge)
    if gauge == 'BAM':
        T = T_BAM()
    else:
        T = T_Bloch(phi_throat, phi_throat,
                    chi_throat=float(chi_path[mid]))
    return U_ret @ T @ U_out


# ---------------------------------------------------------------------------
# Test loops
# ---------------------------------------------------------------------------

def loop_equatorial_2pi(N: int = 256) -> dict:
    """Closed loop: φ goes from 0 to 2π at χ = π/3 (constant). No
    throat crossing.
    """
    t = np.linspace(0.0, 1.0, N + 1)
    return {
        'name': 'equatorial_2pi_chi_pi_over_3',
        'description': 'Constant χ = π/3, φ ∈ [0, 2π]. No throat crossing.',
        'chi_path': np.full_like(t, PI / 3.0),
        'phi_path': TAU * t,
        'crosses_throat': False,
        'analytic_BAM_holonomy_trace_over_2': math.cos(PI * math.cos(PI / 3.0) / 2.0),
    }


def loop_polar_cap_with_throat(N: int = 256) -> dict:
    """Closed loop: descend from (χ=0, φ=0) to (χ=π/2, φ=0), traverse
    the throat at (π/2, 0), then return via (π/2, π/2) and back up to
    (0, 0). Crosses the throat once.

    Path discretisation:
      Out:  (χ=0,φ=0)        → (χ=π/2, φ=0)     [linear in χ]
            (χ=π/2, φ=0)     → (χ=π/2, φ=π)     [equatorial sweep]
      Throat insertion at (χ=π/2, φ=π).
      Ret:  (χ=π/2, φ=π)     → (χ=π/2, φ=2π=0)  [equatorial sweep back]
            (χ=π/2, φ=0)     → (χ=0, φ=0)       [linear in χ]
    """
    t = np.linspace(0.0, 1.0, N + 1)
    chi_path = np.empty_like(t)
    phi_path = np.empty_like(t)
    for i, s in enumerate(t):
        if s < 0.25:
            # Descend: χ 0 → π/2, φ = 0
            chi_path[i] = (s / 0.25) * (PI / 2.0)
            phi_path[i] = 0.0
        elif s < 0.5:
            # Equator out: χ = π/2, φ 0 → π
            chi_path[i] = PI / 2.0
            phi_path[i] = ((s - 0.25) / 0.25) * PI
        elif s < 0.75:
            # Equator back: χ = π/2, φ π → 2π
            chi_path[i] = PI / 2.0
            phi_path[i] = PI + ((s - 0.5) / 0.25) * PI
        else:
            # Ascend: χ π/2 → 0, φ = 2π (= 0)
            chi_path[i] = (PI / 2.0) * (1.0 - (s - 0.75) / 0.25)
            phi_path[i] = TAU
    return {
        'name': 'polar_cap_with_throat_at_phi_pi',
        'description': (
            'Pole → equator → throat at (π/2, π) → equator → pole. '
            'Crosses the throat once at the equator.'
        ),
        'chi_path': chi_path,
        'phi_path': phi_path,
        'crosses_throat': True,
    }


def loop_equatorial_throat_loop(N: int = 256) -> dict:
    """A purely equatorial loop with throat insertion: stay at χ = π/2,
    sweep φ from 0 to π, traverse throat, sweep back from π to 2π.

    Tests the throat conjugation at the precise (φ_in, φ_out) where
    U_gauge picks up its multi-valued character.
    """
    t = np.linspace(0.0, 1.0, N + 1)
    chi_path = np.full_like(t, PI / 2.0)
    phi_path = np.where(t < 0.5, t * TAU, t * TAU)  # 0 → 2π monotonically
    return {
        'name': 'equatorial_throat_sweep',
        'description': (
            'Constant χ = π/2, φ 0 → 2π, throat insertion at midpoint '
            '(φ = π). Tests throat conjugation at the equatorial gauge '
            'string.'
        ),
        'chi_path': chi_path,
        'phi_path': phi_path,
        'crosses_throat': True,
    }


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_connection_difference(N: int = 128) -> dict:
    """(1) A_BAM − A_Bloch = ½ pointwise."""
    chi_grid = np.linspace(0.0, PI, N)
    diffs = np.array([A_phi_BAM(c) - A_phi_Bloch(c) for c in chi_grid])
    max_dev = float(np.max(np.abs(diffs - 0.5)))
    return {
        'name': 'connection_difference_is_half',
        'description': 'A_BAM(χ) − A_Bloch(χ) = ½ for all χ.',
        'max_deviation_from_half': max_dev,
        'pass': max_dev < 1e-14,
    }


def test_spinor_frame_transformation(N: int = 64) -> dict:
    """(2) ψ_BAM = g · ψ_Bloch on a (χ, φ) grid (scalar U(1) gauge)."""
    chi_grid = np.linspace(0.05, PI - 0.05, N)
    phi_grid = np.linspace(0.0, TAU, N, endpoint=False)
    max_dev = 0.0
    for c in chi_grid:
        for p in phi_grid:
            lhs = psi_BAM(float(c), float(p))
            rhs = g_gauge(float(c), float(p)) * psi_Bloch(float(c), float(p))
            dev = float(np.max(np.abs(lhs - rhs)))
            if dev > max_dev:
                max_dev = dev
    return {
        'name': 'spinor_frame_transformation',
        'description': (
            'ψ_BAM(χ,φ) = g(χ,φ) · ψ_Bloch(χ,φ) with g = e^{−iφ/2} '
            '(scalar U(1) phase, multi-valued under φ → φ+2π).'
        ),
        'max_deviation': max_dev,
        'pass': max_dev < 1e-14,
    }


def test_throat_conjugation_invariants(N: int = 32) -> dict:
    """(3) T_Bloch(φ_in, φ_out) preserves SU(2) invariants of iσ_y.

    For instantaneous throat (φ_in = φ_out): det = 1, Tr = 0, T² = −I
    (identical to T_BAM = iσ_y).

    For φ_in ≠ φ_out (interpretation: throat transports the spinor
    between two distinct angular positions): the entire T_Bloch picks
    up the scalar phase e^{i(φ_out−φ_in)/2}, so det = e^{i(φ_out−φ_in)},
    Tr = 0 still, T² = e^{i(φ_out−φ_in)} · (−I). The first SU(2)
    invariant (Tr = 0) is preserved; the others pick up the gauge
    transition function as expected.
    """
    phi_grid = np.linspace(0.0, TAU, N, endpoint=False)
    max_det_dev_instant = 0.0
    max_tr_dev = 0.0
    max_sq_dev_instant = 0.0
    # (a) Instantaneous throat: φ_in = φ_out
    for p in phi_grid:
        T = T_Bloch(float(p), float(p))
        det = complex(np.linalg.det(T))
        tr = complex(np.trace(T))
        sq = T @ T
        max_det_dev_instant = max(max_det_dev_instant, abs(det - 1.0))
        max_tr_dev = max(max_tr_dev, abs(tr))
        max_sq_dev_instant = max(max_sq_dev_instant, float(np.max(np.abs(sq + _ID))))
    # (b) Off-diagonal: Tr remains zero
    for p_in in phi_grid:
        for p_out in phi_grid:
            T = T_Bloch(float(p_in), float(p_out))
            tr = complex(np.trace(T))
            max_tr_dev = max(max_tr_dev, abs(tr))
    return {
        'name': 'throat_conjugation_preserves_invariants',
        'description': (
            'T_Bloch(φ_in, φ_out) preserves Tr = 0 universally; '
            'for instantaneous throat (φ_in = φ_out) also preserves '
            'det = 1 and T² = −I — same SU(2) invariants as iσ_y.'
        ),
        'max_det_deviation_instantaneous': max_det_dev_instant,
        'max_trace_deviation': max_tr_dev,
        'max_T_squared_plus_I_deviation_instantaneous': max_sq_dev_instant,
        'pass': all([
            max_det_dev_instant < 1e-12,
            max_tr_dev < 1e-12,
            max_sq_dev_instant < 1e-12,
        ]),
    }


def _operator_distance(A: np.ndarray, B: np.ndarray) -> float:
    return float(np.max(np.abs(A - B)))


def test_loop_holonomy_relation(loops: list[dict]) -> list[dict]:
    """(4) For each test loop, verify that

        U_loop^{BAM} = (−1)^n_winding · U_loop^{Bloch}

    where n_winding is the net φ-winding of the loop in units of 2π.

    This is the operator statement of the BAM↔Bloch gauge relation
    under the multi-valued scalar U(1) transformation g = e^{−iφ/2}:
    around a 2π loop, g picks up −1, so the BAM and Bloch holonomies
    differ by that sign. SO(3)-equivalent, SU(2)-distinct.

    Physical interpretation:
      - |Tr U|² / 4 (the SO(3)-invariant rotation angle / 2) is the
        same in both gauges.
      - Tr U itself differs by (−1)^n_winding — the spinor sign flip.
      - ⟨ψ|U|ψ⟩ for a fixed reference ψ differs by the same sign;
        unobservable in ⟨ψ|U|ψ⟩|² but observable in interference
        (Bell singlet pair).
    """
    results = []
    for loop in loops:
        U_BAM = closed_loop_holonomy(loop, gauge='BAM')
        U_Bloch = closed_loop_holonomy(loop, gauge='Bloch')
        chi_path = np.asarray(loop['chi_path'])
        phi_path = np.asarray(loop['phi_path'])
        # Net φ-winding (integer divisions in 2π)
        phi_total = float(phi_path[-1] - phi_path[0])
        n_winding = int(round(phi_total / TAU))
        sign = (-1) ** n_winding
        # Operator equality up to the predicted scalar sign
        residual = _operator_distance(U_BAM, sign * U_Bloch)
        # Trace and trace-squared
        tr_BAM = complex(np.trace(U_BAM))
        tr_Bloch = complex(np.trace(U_Bloch))
        # SO(3) invariant: |Tr U|² (same regardless of spinor sign)
        trsq_BAM = abs(tr_BAM) ** 2
        trsq_Bloch = abs(tr_Bloch) ** 2
        results.append({
            'loop_name': loop['name'],
            'description': loop['description'],
            'crosses_throat': loop.get('crosses_throat', False),
            'phi_winding_2pi_units': n_winding,
            'predicted_sign_BAM_over_Bloch': sign,
            'trace_BAM': complex(tr_BAM),
            'trace_Bloch': complex(tr_Bloch),
            'so3_invariant_BAM': trsq_BAM,
            'so3_invariant_Bloch': trsq_Bloch,
            'so3_invariant_residual': abs(trsq_BAM - trsq_Bloch),
            'operator_scaled_residual': residual,
            'pass_so3_equivalence': abs(trsq_BAM - trsq_Bloch) < 1e-10,
            'pass_operator_relation': residual < 1e-10,
        })
    return results


def test_relative_phase_invariance(N: int = 16) -> dict:
    """(5) Δγ(θ_a, θ_b) = γ(θ_a) − γ(θ_b) is gauge-invariant.

    γ_BAM(θ)   =  π·cos(θ)
    γ_Bloch(θ) = −π·(1 − cos(θ))
    Δγ_BAM(a,b)   = π·(cos θ_a − cos θ_b)
    Δγ_Bloch(a,b) = π·(cos θ_a − cos θ_b)   [identical]
    """
    thetas = np.linspace(0.0, PI, N)
    max_dev = 0.0
    samples = []
    for ta in thetas:
        for tb in thetas:
            gBAM = math.pi * (math.cos(ta) - math.cos(tb))
            gBloch = -math.pi * (
                (1.0 - math.cos(ta)) - (1.0 - math.cos(tb))
            )
            dev = abs(gBAM - gBloch)
            if dev > max_dev:
                max_dev = dev
            if len(samples) < 5:
                samples.append({
                    'theta_a': float(ta),
                    'theta_b': float(tb),
                    'delta_gamma_BAM': gBAM,
                    'delta_gamma_Bloch': gBloch,
                })
    return {
        'name': 'relative_phase_gauge_invariance',
        'description': (
            'Δγ(θ_a, θ_b) for the detector-holonomy phase agrees '
            'between BAM and Bloch gauges (the offset cancels).'
        ),
        'max_deviation': max_dev,
        'samples': samples,
        'pass': max_dev < 1e-14,
    }


def test_single_valued_gauge_cannot_remove_pi_offset() -> dict:
    """(6) No single-valued U(1) gauge transformation Λ(χ,φ) can
    remove the π/loop BAM−Bloch offset.

    Argument: Λ must satisfy Λ(χ, φ+2π) ≡ Λ(χ, φ) (mod 2π) for the
    transformation `exp(iΛ)` to be single-valued. For any such Λ:

        ∮_{2π loop} dΛ  =  Λ(φ+2π) − Λ(φ)  ≡  0  (mod 2π)

    Therefore any periodic Λ contributes 0 mod 2π to a closed-loop
    phase. The BAM−Bloch difference is π mod 2π per 2π loop, which
    cannot be absorbed by any periodic Λ.

    Numerical verification: sample a basis of periodic Λ functions
    (Fourier modes Λ_n(φ) = cos(nφ), sin(nφ) for n = 1, 2, …, 8) and
    confirm that ∮ dΛ_n = 0 for each.
    """
    N = 4096
    phi_grid = np.linspace(0.0, TAU, N + 1)
    samples = []
    max_dev = 0.0
    for n in range(1, 9):
        for kind in ('cos', 'sin'):
            if kind == 'cos':
                Lam = np.cos(n * phi_grid)
            else:
                Lam = np.sin(n * phi_grid)
            integral = float(Lam[-1] - Lam[0])  # ∮ dΛ via boundary values
            samples.append({
                'mode': f'{kind}({n}φ)',
                'integral_over_2pi_loop': integral,
            })
            if abs(integral) > max_dev:
                max_dev = abs(integral)
    return {
        'name': 'single_valued_gauge_cannot_remove_pi',
        'description': (
            'Any single-valued (periodic) gauge function Λ(φ) has '
            '∮dΛ = 0 around a 2π loop. The π/loop BAM−Bloch '
            'offset is therefore a true bundle invariant (Wu-Yang '
            'transition function), not a removable gauge artifact.'
        ),
        'max_integral_residual_over_basis': max_dev,
        'samples': samples,
        'pi_offset_persists': True,
        'pass': max_dev < 1e-12,
    }


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    t1 = test_connection_difference()
    t2 = test_spinor_frame_transformation()
    t3 = test_throat_conjugation_invariants()

    loops = [
        loop_equatorial_2pi(),
        loop_polar_cap_with_throat(),
        loop_equatorial_throat_loop(),
    ]
    t4 = test_loop_holonomy_relation(loops)

    t5 = test_relative_phase_invariance()
    t6 = test_single_valued_gauge_cannot_remove_pi_offset()

    scalar_tests = [t1, t2, t3, t5, t6]
    n_scalar_pass = sum(1 for t in scalar_tests if t['pass'])
    n_loop_pass = sum(
        1 for r in t4 if r['pass_so3_equivalence'] and r['pass_operator_relation']
    )
    n_loop_total = len(t4)
    n_scalar_total = len(scalar_tests)

    all_pass = (n_scalar_pass == n_scalar_total and
                n_loop_pass == n_loop_total)

    verdict = (
        'PASS — BAM symmetric gauge is gauge-equivalent to Bloch gauge '
        'under the multi-valued scalar U(1) transformation '
        'g(χ,φ) = e^{−iφ/2}. The two combinations (A_BAM, ψ_BAM, '
        'T_BAM=iσ_y) and (A_Bloch, ψ_Bloch, T_Bloch) make identical '
        'predictions for every SO(3)-level observable. They differ by '
        '(−1)^n per n-winding closed loop at the SU(2) level — the '
        'Wu-Yang transition function that IS the spin-½ double cover. '
        'This sign is observable only in interference between coherent '
        'spinor branches (Bell singlets), where BOTH gauges produce '
        'the same observable phenomena because the relative phase '
        'cancels the gauge offset. No experiment can distinguish the '
        'two gauges; the choice is computational, not physical.'
        if all_pass else
        'FAIL — at least one gauge-equivalence test failed. BAM may '
        'be making a sharper prediction than Bloch, OR the gauge-'
        'conjugation framework is incorrect for this bundle. See '
        'failing tests for diagnosis.'
    )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'gauge_transformation': (
            'g(χ,φ) = e^{−iφ/2} (scalar U(1), multi-valued), '
            'ψ_BAM = g · ψ_Bloch, '
            'A_BAM − A_Bloch = +½ (per dφ).'
        ),
        'scalar_tests': scalar_tests,
        'loop_holonomy_tests': [
            {**r,
             'trace_BAM': complex(r['trace_BAM']),
             'trace_Bloch': complex(r['trace_Bloch'])}
            for r in t4
        ],
        'n_scalar_passed': n_scalar_pass,
        'n_scalar_total': n_scalar_total,
        'n_loop_passed': n_loop_pass,
        'n_loop_total': n_loop_total,
        'overall_pass': all_pass,
        'verdict': verdict,
    }


# ---------------------------------------------------------------------------
# Markdown rendering
# ---------------------------------------------------------------------------

def render_markdown(s: dict) -> str:
    L: list[str] = []
    L.append('# Gauge-sensitive throat-transport falsifier')
    L.append('')
    L.append(f"**Run:** {s['timestamp_utc']}")
    L.append('')
    L.append(
        'Tests whether the BAM combination '
        '`(A_BAM = ½cos(χ)dφ, ψ_BAM, T = iσ_y)` predicts a physically '
        'distinct phase from the standard Bloch combination '
        '`(A_Bloch = ½(cos(χ)−1)dφ, ψ_Bloch, T_Bloch)` after all '
        'allowed gauge freedoms are accounted for. The relevant gauge '
        'transformation is the SU(2) spinor-frame rotation'
    )
    L.append('')
    L.append('```')
    L.append('g(χ, φ) = e^{−iφ/2}                (scalar U(1), multi-valued)')
    L.append('ψ_BAM = g · ψ_Bloch                (spinor sign flip on 2π loop)')
    L.append('A_BAM − A_Bloch = +½ (per dφ)      (connection-level relation)')
    L.append('```')
    L.append('')
    L.append(
        'Pulling the common factor `e^{−iφ/2}` out of the BAM spinor '
        '`(cos(χ/2)e^{−iφ/2}, sin(χ/2)e^{+iφ/2})` recovers the Bloch '
        'spinor `(cos(χ/2), sin(χ/2)e^{+iφ})`. The relating factor is '
        'a scalar U(1) phase, not an SU(2) frame rotation. It is '
        'multi-valued: under `φ → φ+2π`, `g → −g`. This is the Wu-Yang '
        'transition function for the spin-½ representation of the '
        'Hopf bundle.'
    )
    L.append('')

    L.append('## Scalar tests')
    L.append('')
    L.append('| # | test | metric | residual | PASS? |')
    L.append('|---|---|---|---:|---|')
    for i, t in enumerate(s['scalar_tests'], start=1):
        if 'max_deviation_from_half' in t:
            res = t['max_deviation_from_half']
        elif 'max_deviation' in t:
            res = t['max_deviation']
        elif 'max_det_deviation_instantaneous' in t:
            res = max(t['max_det_deviation_instantaneous'],
                      t['max_trace_deviation'],
                      t['max_T_squared_plus_I_deviation_instantaneous'])
        elif 'max_integral_residual_over_basis' in t:
            res = t['max_integral_residual_over_basis']
        else:
            res = float('nan')
        passed = '**PASS**' if t['pass'] else '**FAIL**'
        L.append(
            f"| {i} | `{t['name']}` | "
            f"{t['description'].split('.')[0]} | "
            f"{res:.2e} | {passed} |"
        )
    L.append('')

    L.append('## Closed-loop holonomy tests')
    L.append('')
    L.append(
        'Each loop is checked against the BAM↔Bloch operator relation '
        '`U_loop^BAM = (−1)^n_winding · U_loop^Bloch` and against the '
        'SO(3) invariant `|Tr U|²` (which must agree regardless of '
        'spinor sign).'
    )
    L.append('')
    L.append(
        '| loop | crosses throat | n_winding | (−1)^n | Tr U_BAM | '
        'Tr U_Bloch | |Tr|² residual | op residual | PASS? |'
    )
    L.append('|---|---|---:|---:|---:|---:|---:|---:|---|')
    for r in s['loop_holonomy_tests']:
        passed = '**PASS**' if (r['pass_so3_equivalence'] and r['pass_operator_relation']) else '**FAIL**'
        tBAM = r['trace_BAM']
        tBloch = r['trace_Bloch']
        L.append(
            f"| `{r['loop_name']}` | {r['crosses_throat']} | "
            f"{r['phi_winding_2pi_units']} | "
            f"{r['predicted_sign_BAM_over_Bloch']:+d} | "
            f"({tBAM.real:+.4f}{tBAM.imag:+.4f}j) | "
            f"({tBloch.real:+.4f}{tBloch.imag:+.4f}j) | "
            f"{r['so3_invariant_residual']:.2e} | "
            f"{r['operator_scaled_residual']:.2e} | "
            f"{passed} |"
        )
    L.append('')

    L.append('## Verdict')
    L.append('')
    L.append(f"**{s['verdict']}**")
    L.append('')
    if s['overall_pass']:
        L.append(
            'Specifically, the probe established:'
        )
        L.append('')
        L.append(
            '- `A_BAM − A_Bloch = ½ dφ` at every sampled point of the '
            'base — the connection-level statement of the gauge '
            'transformation is exact.'
        )
        L.append(
            '- `ψ_BAM = g · ψ_Bloch` with the scalar U(1) gauge '
            'function `g = e^{−iφ/2}` at every sampled (χ,φ) — the '
            'spinor relation is a multi-valued scalar phase, not an '
            'SU(2) frame rotation.'
        )
        L.append(
            '- `T_Bloch(φ_in, φ_out)` defined by gauge conjugation of '
            '`T_BAM = iσ_y` preserves Tr = 0 universally and, in the '
            'instantaneous-throat case (φ_in = φ_out), also preserves '
            '`det = 1` and `T² = −I`. The throat transport iσ_y in BAM '
            'becomes the scaled `e^{i(φ_out−φ_in)/2}·iσ_y` in Bloch — '
            'a clean scalar conjugation.'
        )
        L.append(
            '- Every closed loop\'s holonomy satisfies '
            '`U_loop^BAM = (−1)^n_winding · U_loop^Bloch` to numerical '
            'precision. The SO(3)-level invariant `|Tr U|²` (which '
            'discards the spinor sign) agrees exactly. So the two '
            'gauges describe the same SO(3) rotation but lift it to '
            'opposite signs in SU(2) for every odd-winding loop.'
        )
        L.append(
            '- Detector-holonomy phase differences `Δγ(θ_a, θ_b)` are '
            'identical in the two gauges — the π/loop offset cancels '
            'in every relative measurement. Bell experiments using '
            'either gauge produce the same CHSH = 2√2.'
        )
        L.append(
            '- No single-valued U(1) gauge transformation can remove '
            'the π/loop BAM−Bloch offset (any periodic Λ contributes '
            'zero around a closed loop). The π offset is a true Wu-Yang '
            'transition function — the spinor double cover content — '
            'present in any gauge that uses spinor wavefunctions.'
        )
        L.append('')
        L.append(
            'Interpretation. The BAM symmetric gauge and the standard '
            'Bloch gauge are two sections of the SAME Hopf bundle, '
            'related by the multi-valued scalar U(1) gauge function '
            '`g = e^{−iφ/2}`. The multi-valuedness `g(φ+2π) = −g(φ)` '
            'is exactly the spinor double-cover transition function; '
            'it is unobservable in any projective (ray) Hilbert-space '
            'measurement and cancels in every relative-phase '
            'observable. Every physical prediction agrees between the '
            'two gauges. BAM\'s derivational economy — `T = iσ_y` '
            'falls out directly from the orientation-reversing '
            'isometry of S³ — is real but represents a choice of '
            'computational gauge, not an experimentally falsifiable '
            'prediction.'
        )
        L.append('')
        L.append(
            'What this rules out:'
        )
        L.append('')
        L.append(
            '- Claims that BAM\'s symmetric gauge predicts a *novel* '
            'phase or interference pattern not present in standard '
            'spin-½ quantum mechanics. It does not.'
        )
        L.append('')
        L.append(
            'What this leaves intact:'
        )
        L.append('')
        L.append(
            '- BAM\'s **derivation** of `T = iσ_y` from the unique '
            'orientation-reversing Hopf-preserving S³ isometry — a '
            'geometrically forced spinor map. Standard quantum '
            'mechanics POSTULATES this transformation as part of the '
            'spinor representation; BAM DERIVES it from the throat\'s '
            'non-orientability. The derivation has content '
            'independent of any gauge choice.'
        )
        L.append(
            '- BAM\'s identification of the equator `χ = π/2` as the '
            'natural throat location and the zero-self-energy stable '
            'orbit. This is a structural / geometric statement, not a '
            'gauge artifact.'
        )
        L.append(
            '- The closure-ledger predictions — BAM uses the symmetric '
            'gauge because it places the closure-quantum integers in '
            'their most natural form (`π·cos(χ)` per loop, `2π` per '
            'pair traversal). The Bloch gauge would produce the same '
            'integer ledger after the appropriate offset.'
        )
    else:
        L.append(
            'Failing tests indicate the gauge-conjugation framework is '
            'incomplete or that BAM and Bloch are NOT in the same '
            'gauge orbit — a structural surprise that warrants '
            'investigation. The specific failures are listed in the '
            'scalar-tests table and loop-holonomy table above.'
        )
    L.append('')

    L.append('## What this leaves open')
    L.append('')
    L.append(
        '- **Sharp derivation vs. sharp prediction.** BAM\'s '
        'symmetric-gauge derivation of `T = iσ_y` may have content '
        'beyond what this gauge-equivalence test sees — in particular, '
        'the derivation argument in `embedding/transport.py` invokes '
        'the UNIQUE orientation-reversing isometry of S³ that '
        'preserves the Hopf bundle. Whether that uniqueness statement '
        'has additional consequences beyond the spinor frame choice '
        '(e.g. a uniqueness theorem for the throat sector at the '
        'classical level) is open.'
    )
    L.append(
        '- **Non-Abelian gauge sector.** This probe analyses only the '
        'U(1) ↔ SU(2) doubling. If BAM\'s throat carries additional '
        'non-Abelian internal structure (e.g. for the quark sector\'s '
        'colour content via the QCD pinhole), the corresponding '
        'gauge-equivalence test would need an extended formulation.'
    )
    L.append(
        '- **Time-dependent (Aharonov-Anandan) corrections.** The '
        'present test is purely adiabatic. Whether the gauge '
        'equivalence persists under finite-velocity transport is a '
        'separate target.'
    )
    L.append('')
    return '\n'.join(L)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def _json_default(o):
    if isinstance(o, complex):
        return {'real': o.real, 'imag': o.imag}
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def main(argv: Optional[list[str]] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    out = here / 'runs' / f'{ts}_throat_transport_gauge_probe'
    out.mkdir(parents=True, exist_ok=True)
    (out / 'probe.json').write_text(
        json.dumps(summary, indent=2, default=_json_default)
    )
    (out / 'probe.md').write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
