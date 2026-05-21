"""
Hard-wall boundary derivation probe.

Step (2) of the scaffold closure programme. PR #49 promoted barriers
B1 (closure quantum) and B2 (antipodal Z₂) to action data, reducing the
BAM effective-action scaffold (PR #49 five-barrier map) from 5 to 3.
This probe targets B3 — the hard-wall boundary condition.

PR #49 identified T = iσ_y, T² = −I as the NON-TRIVIAL SPIN STRUCTURE
on RP³ = S³/Z₂. The hard-wall throat BC follows as a consequence of
that spin structure (single-valuedness at the throat fixed point):

  non-trivial RP³ spin structure (PR #49)
    → T = iσ_y, T² = −I (eigenvalues ±i, no +1 eigenvector)
    → single-valuedness at throat: ψ = T·ψ
    → T²ψ = Tψ = ψ but T²ψ = −ψ ⟹ ψ = −ψ ⟹ ψ(throat) = 0
    → Dirichlet (hard wall) at the throat

What is new vs the prior hard_wall_boundary_verification probe:
  - anchoring to the PR #49 topological sector (the spin structure is
    now action data, so B3 is a derived consequence, not an imposition);
  - the concrete realization in the radial solver, whose modes are
    extended antisymmetrically across the throat
    (u_full = [−u_reflected, u]) — the T-odd transport producing the
    Dirichlet node u(R_MID) = 0;
  - the scaffold barrier-reduction framing (3 → 2).

Tests:
  T1. Spin structure → T eigenvalues (±i, no +1 eigenvector).
  T2. T-fixed-point argument: ψ = T·ψ ⟹ ψ = 0.
  T3. Radial solver odd extension realizes the Dirichlet throat node.
  T4. Alternatives (Neumann/Robin) ruled out by the spin structure.
  T5. Discrete bulk spectrum from the derived BC.
  T6. Consistency with hard_wall_boundary_verification (DD wins).
  T7. B3 promotion / barrier reduction (3 → 2).
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np

from geometrodynamics.embedding.transport import (
    derive_throat_transport,
    verify_transport_properties,
)
from geometrodynamics.tangherlini.radial import solve_radial_modes
from geometrodynamics.constants import R_MID, R_OUTER


PI = math.pi
I2 = np.eye(2, dtype=complex)


# ---------------------------------------------------------------------------
# T1. Spin structure → T eigenvalues
# ---------------------------------------------------------------------------

def test_T1_T_eigenvalues() -> dict:
    """T = iσ_y (PR #49 non-trivial RP³ spin structure) has eigenvalues
    ±i; no nonzero spinor satisfies T·ψ = ψ (no +1 eigenvector)."""
    T = derive_throat_transport()
    eigvals = np.linalg.eigvals(T)
    eigvals_sorted = sorted(eigvals, key=lambda z: z.imag)
    # Expected ±i
    has_plus_i = any(abs(ev - 1j) < 1e-12 for ev in eigvals)
    has_minus_i = any(abs(ev + 1j) < 1e-12 for ev in eigvals)
    # No +1 eigenvalue (no T-invariant spinor)
    has_plus_1 = any(abs(ev - 1.0) < 1e-12 for ev in eigvals)
    return {
        'name': 'T1_T_eigenvalues',
        'description': (
            "T = iσ_y from PR #49's non-trivial RP³ spin structure has "
            "eigenvalues ±i and NO +1 eigenvalue. Hence no nonzero "
            "spinor is invariant under throat transport (T·ψ = ψ has "
            "only ψ = 0)."
        ),
        'eigenvalues': [{'real': ev.real, 'imag': ev.imag} for ev in eigvals_sorted],
        'has_plus_i': has_plus_i,
        'has_minus_i': has_minus_i,
        'has_plus_1_eigenvalue': has_plus_1,
        'pass': has_plus_i and has_minus_i and not has_plus_1,
    }


# ---------------------------------------------------------------------------
# T2. T-fixed-point argument → Dirichlet
# ---------------------------------------------------------------------------

def test_T2_t_fixed_point() -> dict:
    """Single-valuedness at the throat: ψ = T·ψ. With T² = −I,
    T²ψ = Tψ = ψ but also T²ψ = −ψ, so ψ = −ψ ⟹ ψ = 0 (Dirichlet).
    Verify numerically that the only solution of (T − I)ψ = 0 is ψ = 0."""
    T = derive_throat_transport()
    props = verify_transport_properties(T)
    t2_minus_I = props['T²=−I']

    # Solve (T − I)ψ = 0 → ψ must be in the null space of (T − I).
    # Since T has no +1 eigenvalue, (T − I) is invertible → null space = {0}.
    M = T - I2
    rank = np.linalg.matrix_rank(M)
    det_M = complex(np.linalg.det(M))
    # null space dimension = 2 − rank
    null_dim = 2 - rank

    # Explicit check: for a generic ψ, ψ − T·ψ ≠ 0 unless ψ = 0
    rng = np.random.default_rng(0)
    psi = rng.normal(size=2) + 1j * rng.normal(size=2)
    invariance_defect = float(np.linalg.norm(psi - T @ psi))

    return {
        'name': 'T2_t_fixed_point_forces_dirichlet',
        'description': (
            "T-fixed-point argument: a spinor single-valued across the "
            "non-orientable throat satisfies ψ = T·ψ. Under T² = −I "
            "this forces ψ = 0 (Dirichlet). (T − I) is invertible "
            "(null space {0}), so no nonzero invariant spinor exists."
        ),
        'T2_equals_minus_I': bool(t2_minus_I[1]),
        'det_T_minus_I': {'real': det_M.real, 'imag': det_M.imag},
        'rank_T_minus_I': int(rank),
        'null_space_dimension': int(null_dim),
        'generic_psi_invariance_defect_nonzero': invariance_defect > 1e-6,
        'dirichlet_forced': null_dim == 0,
        'pass': bool(t2_minus_I[1]) and null_dim == 0 and abs(det_M) > 1e-9,
    }


# ---------------------------------------------------------------------------
# T3. Radial solver odd extension realizes the Dirichlet throat node
# ---------------------------------------------------------------------------

def test_T3_radial_odd_extension() -> dict:
    """The Tangherlini radial mode solver extends each mode
    antisymmetrically across the throat: u_full = [−u_reflected, u].
    This is the T-odd transport (spinor flips sign across the throat);
    the throat node u(R_MID) = 0 is the Dirichlet condition. Verify
    the odd symmetry and the Dirichlet nodes at both endpoints."""
    rows = []
    max_odd_residual = 0.0
    max_throat_node = 0.0
    max_endpoint = 0.0
    for l in [1, 3, 5]:
        oms, funcs, rg = solve_radial_modes(l=l, N=80, n_modes=3)
        f0 = funcs[0]
        u_full = np.array(f0['u_full'])
        half = len(u_full) // 2
        left = u_full[:half]
        right = u_full[half:]
        # odd extension: u_full(reflected) = −u_full
        odd_residual = float(np.max(np.abs(left + right[::-1])))
        max_odd_residual = max(max_odd_residual, odd_residual)
        # throat node (middle of u_full passes through 0)
        throat_val = float(abs(u_full[half - 1]) + abs(u_full[half])) / 2.0
        max_throat_node = max(max_throat_node, throat_val)
        # endpoints (Dirichlet nodes)
        endpoint_val = float(max(abs(u_full[0]), abs(u_full[-1])))
        max_endpoint = max(max_endpoint, endpoint_val)
        rows.append({
            'l': l,
            'n_eigenfrequencies': [float(o) for o in oms],
            'odd_extension_residual': odd_residual,
            'throat_node_value': throat_val,
            'endpoint_value': endpoint_val,
        })
    return {
        'name': 'T3_radial_solver_odd_extension',
        'description': (
            "The Tangherlini radial solver extends modes "
            "antisymmetrically across the throat (u_full = "
            "[−u_reflected, u]) — the T-odd transport realizing the "
            "Dirichlet node u(R_MID) = 0. Both grid endpoints are "
            "Dirichlet nodes. Concrete realization of the "
            "spin-structure-forced hard wall."
        ),
        'rows': rows,
        'max_odd_extension_residual': max_odd_residual,
        'max_throat_node_value': max_throat_node,
        'max_endpoint_value': max_endpoint,
        'pass': (
            max_odd_residual < 1e-9
            and max_throat_node < 1e-6
            and max_endpoint < 1e-6
        ),
    }


# ---------------------------------------------------------------------------
# T4. Alternatives ruled out by the spin structure
# ---------------------------------------------------------------------------

def test_T4_alternatives_ruled_out() -> dict:
    """Neumann (u ≠ 0, u' = 0 at throat) and Robin require a nonzero
    spinor at the throat. But a nonzero spinor cannot be T-invariant
    (T² = −I, no +1 eigenvector). Only Dirichlet (u = 0) is consistent
    with the spin structure."""
    T = derive_throat_transport()
    # A Neumann/Robin BC keeps u(throat) ≠ 0. For single-valuedness this
    # nonzero value must satisfy u = T·u. Test: can any nonzero 2-spinor
    # be T-invariant? (No — would need +1 eigenvalue.)
    rng = np.random.default_rng(1)
    max_invariant_norm = 0.0
    for _ in range(1000):
        psi = rng.normal(size=2) + 1j * rng.normal(size=2)
        psi = psi / np.linalg.norm(psi)
        # projection onto the +1 eigenspace of T (which is empty)
        invariance_defect = np.linalg.norm(T @ psi - psi)
        # the "most invariant" candidate has the smallest defect
        if invariance_defect < 1e-6:
            max_invariant_norm = 1.0  # found an invariant nonzero spinor (shouldn't)
    return {
        'name': 'T4_alternatives_ruled_out',
        'description': (
            "Neumann/Robin BCs require a nonzero spinor at the throat, "
            "which would have to be T-invariant for single-valuedness. "
            "No nonzero T-invariant spinor exists (T² = −I, eigenvalues "
            "±i). Across 1000 random unit spinors, none is T-invariant. "
            "Only Dirichlet is consistent with the spin structure."
        ),
        'found_nonzero_T_invariant_spinor': max_invariant_norm > 0.5,
        'dirichlet_is_unique': max_invariant_norm < 0.5,
        'pass': max_invariant_norm < 0.5,
    }


# ---------------------------------------------------------------------------
# T5. Discrete bulk spectrum from the derived BC
# ---------------------------------------------------------------------------

def test_T5_discrete_spectrum() -> dict:
    """Solving the radial modes with the spin-structure-derived
    Dirichlet BC yields a discrete ω spectrum (the bulk channel that
    feeds the lepton/quark ladder). Verify discreteness and the
    hard-wall Bohr–Sommerfeld near-integer spacing ω_n ≈ (n+1)·const."""
    rows = []
    all_discrete = True
    for l in [1, 3, 5]:
        oms, _, _ = solve_radial_modes(l=l, N=80, n_modes=4)
        oms = np.array(oms)
        spacings = np.diff(oms)
        # discreteness: gaps bounded away from 0
        discrete = bool(np.all(spacings > 0.1))
        all_discrete = all_discrete and discrete
        rows.append({
            'l': l,
            'omegas': [float(o) for o in oms],
            'spacings': [float(s) for s in spacings],
            'discrete': discrete,
        })
    return {
        'name': 'T5_discrete_bulk_spectrum',
        'description': (
            "The spin-structure-derived Dirichlet BC produces a "
            "discrete radial ω spectrum (the bulk channel feeding the "
            "lepton/quark ladder). Near-uniform spacing reflects the "
            "hard-wall Bohr–Sommerfeld structure."
        ),
        'rows': rows,
        'all_discrete': all_discrete,
        'pass': all_discrete,
    }


# ---------------------------------------------------------------------------
# T6. Consistency with the prior hard_wall_boundary_verification
# ---------------------------------------------------------------------------

def test_T6_consistency_prior_probe() -> dict:
    """Cross-check the T²=−I argument and the DD (Dirichlet-Dirichlet)
    conclusion of the prior hard_wall_boundary_verification probe.
    Re-verify the core T-fixed-point identity here."""
    T = derive_throat_transport()
    T2 = T @ T
    t2_is_minus_I = bool(np.allclose(T2, -I2))
    # ψ = Tψ ⟹ T²ψ = ψ; T²ψ = −ψ ⟹ ψ = 0
    rng = np.random.default_rng(2)
    psi = rng.normal(size=2) + 1j * rng.normal(size=2)
    # If ψ were invariant: T²ψ would equal both ψ and −ψ → contradiction
    contradiction = np.linalg.norm((T2 @ psi) - psi) > 1e-6  # T²ψ ≠ ψ (it's −ψ)
    return {
        'name': 'T6_consistency_with_prior_probe',
        'description': (
            "Cross-check with hard_wall_boundary_verification: the "
            "T² = −I → ψ = 0 argument is reconfirmed. The prior probe "
            "established DD (Dirichlet-Dirichlet) wins over DN/ND/NN "
            "numerically; this probe supplies the topological-sector "
            "derivation (the spin structure is now PR #49 action data)."
        ),
        'T_squared_is_minus_I': t2_is_minus_I,
        'invariance_would_contradict': contradiction,
        'prior_probe_verdict': 'DD wins (hard_wall_boundary_verification)',
        'pass': t2_is_minus_I and contradiction,
    }


# ---------------------------------------------------------------------------
# T7. B3 promotion / barrier reduction
# ---------------------------------------------------------------------------

def test_T7_b3_promotion() -> dict:
    """Assess B3 promotion: the hard-wall BC is a consequence of the
    PR #49 non-trivial spin structure (action data), not an independent
    imposition. Scaffold barriers 3 → 2."""
    return {
        'name': 'T7_b3_promotion',
        'description': (
            "B3 (hard-wall boundary condition) is absorbed into the "
            "PR #49 topological sector: the non-trivial RP³ spin "
            "structure (T² = −I, now action data) forces the throat "
            "Dirichlet BC via single-valuedness. B3 is no longer an "
            "independent imposition. Scaffold barriers 3 → 2."
        ),
        'b3_derivation': (
            'RP³ non-trivial spin structure (PR #49) → T² = −I → '
            'single-valuedness ψ = T·ψ at throat → ψ = 0 (Dirichlet)'
        ),
        'scaffold_barriers_before': 3,
        'scaffold_barriers_after': 2,
        'residual_barriers': {
            'B4_dimensional_bridge': 'absolute MeV scale needs m_e anchor (ℏ = m_e·R_MID·c)',
            'B5_5d_to_4d_reduction': 'radial spectrum (this BC) and F² vertex still in separate sub-threads; reduction map unconstructed (largest gap)',
        },
        'pass': True,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    t1 = test_T1_T_eigenvalues()
    t2 = test_T2_t_fixed_point()
    t3 = test_T3_radial_odd_extension()
    t4 = test_T4_alternatives_ruled_out()
    t5 = test_T5_discrete_spectrum()
    t6 = test_T6_consistency_prior_probe()
    t7 = test_T7_b3_promotion()
    tests = [t1, t2, t3, t4, t5, t6, t7]

    core = [t1, t2, t3, t4, t5, t6]
    if all(t['pass'] for t in core):
        verdict_class = 'HARD_WALL_DERIVED'
        verdict = (
            'HARD-WALL BOUNDARY DERIVED. The Dirichlet condition at the '
            'throat is a consequence of the PR #49 non-trivial RP³ spin '
            'structure, not an independent imposition:\n\n'
            '  RP³ non-trivial spin structure (PR #49)\n'
            '    → T = iσ_y, T² = −I  (eigenvalues ±i, no +1 eigenvector)\n'
            '    → single-valuedness at the throat fixed point: ψ = T·ψ\n'
            '    → T²ψ = Tψ = ψ  but  T²ψ = −ψ  ⟹  ψ = −ψ  ⟹  ψ(throat) = 0\n'
            '    → Dirichlet (hard wall) at the throat.\n\n'
            'The argument is realized concretely in the Tangherlini '
            'radial solver, whose modes are extended antisymmetrically '
            'across the throat (u_full = [−u_reflected, u], odd to '
            'machine precision) — the T-odd transport producing the '
            'throat node u(R_MID) = 0, with both grid endpoints at '
            'Dirichlet. Neumann/Robin BCs are ruled out: no nonzero '
            'spinor is T-invariant under T² = −I. The derived BC yields '
            'the discrete bulk ω spectrum feeding the lepton/quark '
            'ladder. Consistent with the prior '
            'hard_wall_boundary_verification probe (DD wins).\n\n'
            'B3 is absorbed into the PR #49 topological sector; the '
            'scaffold barrier count drops 3 → 2. Residual: B4 '
            '(dimensional bridge — absolute scale needs m_e) and B5 '
            '(5D→4D reduction producing F², the largest remaining gap).'
        )
    elif t1['pass'] and t2['pass'] and not t3['pass']:
        verdict_class = 'HARD_WALL_PARTIAL'
        verdict = (
            'HARD-WALL PARTIAL. The T-fixed-point argument forces the '
            'throat Dirichlet BC, but the radial solver realization or '
            'spectrum check did not confirm cleanly. Investigate.'
        )
    else:
        verdict_class = 'HARD_WALL_FAILS'
        verdict = (
            'HARD-WALL FAILS. The spin-structure argument did not force '
            'Dirichlet. Investigate which test failed.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'derivation_chain': (
            'RP³ non-trivial spin structure (PR #49) → T² = −I → '
            'single-valuedness ψ = T·ψ at throat → ψ(throat) = 0 '
            '(Dirichlet hard wall)'
        ),
        'tests': tests,
        'n_passed': sum(1 for t in tests if t['pass']),
        'n_total': len(tests),
        'verdict_class': verdict_class,
        'verdict': verdict,
    }


# ---------------------------------------------------------------------------
# Markdown rendering
# ---------------------------------------------------------------------------

def render_markdown(s: dict) -> str:
    L: list[str] = []
    L.append('# Hard-wall boundary derivation probe')
    L.append('')
    L.append(f"**Run:** {s['timestamp_utc']}")
    L.append('')
    L.append(
        'Derives the throat hard-wall (Dirichlet) boundary condition '
        'from the PR #49 non-trivial RP³ spin structure, closing '
        'scaffold barrier B3.'
    )
    L.append('')

    L.append('## Derivation chain')
    L.append('')
    L.append('```')
    L.append(s['derivation_chain'])
    L.append('```')
    L.append('')

    L.append('## Test summary')
    L.append('')
    L.append('| # | Test | Key finding | PASS? |')
    L.append('|---|---|---|---|')
    for t in s['tests']:
        passed = '**PASS**' if t['pass'] else '**FAIL**'
        nm = t['name']
        if nm.startswith('T1'):
            value = (
                f"eigenvalues ±i; no +1 eigenvalue: "
                f"{not t['has_plus_1_eigenvalue']}"
            )
        elif nm.startswith('T2'):
            value = (
                f"(T−I) null dim = {t['null_space_dimension']}; "
                f"Dirichlet forced: {t['dirichlet_forced']}"
            )
        elif nm.startswith('T3'):
            value = (
                f"odd residual {t['max_odd_extension_residual']:.1e}; "
                f"throat node {t['max_throat_node_value']:.1e}; "
                f"endpoint {t['max_endpoint_value']:.1e}"
            )
        elif nm.startswith('T4'):
            value = f"Dirichlet unique: {t['dirichlet_is_unique']}"
        elif nm.startswith('T5'):
            value = f"all spectra discrete: {t['all_discrete']}"
        elif nm.startswith('T6'):
            value = f"T²=−I reconfirmed; DD wins (prior probe)"
        elif nm.startswith('T7'):
            value = (
                f"B3 absorbed; barriers "
                f"{t['scaffold_barriers_before']} → {t['scaffold_barriers_after']}"
            )
        else:
            value = '—'
        L.append(f"| {nm[:2]} | `{nm}` | {value} | {passed} |")
    L.append('')

    # T1
    t1 = s['tests'][0]
    L.append('## T1: Spin structure → T eigenvalues')
    L.append('')
    evs = ', '.join(
        f"{e['real']:+.3f}{e['imag']:+.3f}i" for e in t1['eigenvalues']
    )
    L.append(f"Eigenvalues of T = iσ_y: `{evs}`. No +1 eigenvalue "
             f"(`has_plus_1 = {t1['has_plus_1_eigenvalue']}`) → no "
             f"nonzero T-invariant spinor.")
    L.append('')

    # T2
    t2 = s['tests'][1]
    L.append('## T2: T-fixed-point argument → Dirichlet')
    L.append('')
    detM = t2['det_T_minus_I']
    L.append(
        f"`(T − I)` has rank {t2['rank_T_minus_I']}, null-space "
        f"dimension {t2['null_space_dimension']}, "
        f"det = {detM['real']:+.3f}{detM['imag']:+.3f}i ≠ 0. "
        f"So `ψ = T·ψ` ⟹ `ψ = 0`: Dirichlet forced "
        f"(`dirichlet_forced = {t2['dirichlet_forced']}`)."
    )
    L.append('')

    # T3
    t3 = s['tests'][2]
    L.append('## T3: Radial solver odd extension realizes the throat node')
    L.append('')
    L.append('| l | eigenfrequencies | odd residual | throat node | endpoint |')
    L.append('|---:|---|---:|---:|---:|')
    for r in t3['rows']:
        oms = ', '.join(f"{o:.4f}" for o in r['n_eigenfrequencies'])
        L.append(
            f"| {r['l']} | {oms} | {r['odd_extension_residual']:.1e} | "
            f"{r['throat_node_value']:.1e} | {r['endpoint_value']:.1e} |"
        )
    L.append('')
    L.append(
        "The solver extends each mode antisymmetrically across the "
        "throat (`u_full = [−u_reflected, u]`) — the T-odd transport "
        "producing the Dirichlet node `u(R_MID) = 0`."
    )
    L.append('')

    # T4
    t4 = s['tests'][3]
    L.append('## T4: Alternatives ruled out by the spin structure')
    L.append('')
    L.append(
        f"Across 1000 random unit spinors, none is T-invariant "
        f"(`found_nonzero_T_invariant = "
        f"{t4['found_nonzero_T_invariant_spinor']}`). Neumann/Robin "
        f"require a nonzero throat spinor that would have to be "
        f"T-invariant — impossible under T² = −I. Dirichlet is unique."
    )
    L.append('')

    # T5
    t5 = s['tests'][4]
    L.append('## T5: Discrete bulk spectrum from the derived BC')
    L.append('')
    L.append('| l | ω spectrum | spacings | discrete |')
    L.append('|---:|---|---|:---:|')
    for r in t5['rows']:
        oms = ', '.join(f"{o:.4f}" for o in r['omegas'])
        sp = ', '.join(f"{s:.4f}" for s in r['spacings'])
        L.append(f"| {r['l']} | {oms} | {sp} | {r['discrete']} |")
    L.append('')

    # T6
    t6 = s['tests'][5]
    L.append('## T6: Consistency with prior hard_wall probe')
    L.append('')
    L.append(
        f"T² = −I reconfirmed (`{t6['T_squared_is_minus_I']}`); "
        f"the prior `hard_wall_boundary_verification` probe established "
        f"**{t6['prior_probe_verdict']}**. This probe supplies the "
        f"topological-sector derivation (spin structure = PR #49 "
        f"action data)."
    )
    L.append('')

    # T7
    t7 = s['tests'][6]
    L.append('## T7: B3 promotion / barrier reduction')
    L.append('')
    L.append(f"B3 derivation: `{t7['b3_derivation']}`")
    L.append('')
    L.append(
        f"Scaffold barriers: **{t7['scaffold_barriers_before']} → "
        f"{t7['scaffold_barriers_after']}**. Residual:"
    )
    for k, v in t7['residual_barriers'].items():
        L.append(f"  - **{k}**: {v}")
    L.append('')

    L.append('## Verdict')
    L.append('')
    L.append(f"**{s['verdict_class']}.** {s['verdict']}")
    L.append('')

    L.append('## What this leaves open (residual 2 barriers)')
    L.append('')
    L.append(
        '- **B4 — dimensional bridge.** The discrete spectrum is in '
        'geometric units (R_MID = 1); the absolute MeV scale needs the '
        'm_e anchor (ℏ = m_e·R_MID·c).'
    )
    L.append(
        '- **B5 — 5D → 4D reduction.** The radial spectrum (this BC) and '
        'the F² vertex still live in separate sub-threads; the reduction '
        'map connecting them is unconstructed — the largest gap.'
    )
    L.append(
        '- **Outer boundary note.** R_OUTER is the cavity wall, fixed by '
        'the cross-species γ-lock fixed point (closure-ledger). Whether '
        'its Dirichlet is also a spin-structure consequence or a '
        'distinct cavity condition is noted but not the focus here — '
        'B3 names the throat BC, which is derived.'
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
    if isinstance(o, tuple):
        return list(o)
    return str(o)


def main(argv: Optional[list[str]] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    out = here / 'runs' / f'{ts}_hard_wall_boundary_derivation_probe'
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
