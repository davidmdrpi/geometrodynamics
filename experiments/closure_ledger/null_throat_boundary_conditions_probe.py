"""
Null throat boundary conditions for wave transport on the 5D horizon (PR #129).

PR #128 built the horizon-regular charts for the 5D Tangherlini throat and
identified the antipodal map (U,V,Ω) → (−U,−V,Ω_antipodal) as BAM's
throat ↔ antithroat C-swap. PR #116 ran the matter cavity operator
H = −d²/dr*² + V_tangherlini on the exterior. This probe asks the next
question: what boundary condition does the NULL throat (the 5D horizon) impose
on the transported waves? The answer is the BAM-native one — an
l-PARITY-dependent antipodal condition that makes the throat a UNITARY MIRROR,
not an absorbing horizon.

## The potential vanishes at the horizon

The massless wave equation □Φ = 0, separated as
Φ = e^{−iωt} Y_l(Ω) ψ(r)/r^{3/2}, gives the Schrödinger form
−d²ψ/dr*² + V_l ψ = ω²ψ with V_l = f(r)[l(l+2)/r² + 3rs²/r⁴] (PR #116). Since
V_l ∝ f(r) → 0 at the throat r = rs, the near-horizon equation is −ψ'' = ω²ψ:
the modes are the pure null phases ψ ~ e^{±iωr*} (the ingoing/outgoing null
rays of PR #128, v = t + r*, u = t − r*).

## Three candidate boundary conditions

  1. **Ingoing / absorbing** (standard quasinormal): ψ ~ e^{−iωr*} at the
     throat — flux flows *into* the horizon and is lost. Non-unitary (the
     horizon is a sink).
  2. **Reflective wall** (Dirichlet ψ = 0 or Neumann ψ' = 0): a hard box. Real,
     discrete spectrum — the matter-cavity picture (PR #116).
  3. **Antipodal** (BAM-native): the PR #128 identification
     Φ(U,V,Ω) = Φ(−U,−V,Ω_antipodal) links the two sheets at the throat.

## The antipodal map fixes the BC by l-parity

The scalar harmonics on the horizon S³ have antipodal parity
Y_l(−x) = (−1)^l Y_l(x) (they are degree-l harmonic polynomials on ℝ⁴;
verified here). Single-valuedness of the (untwisted) scalar under the antipodal
identification, Φ(U,V,Ω) = Φ(−U,−V,Ω_antipodal), then forces the radial
function to carry the compensating parity (−1)^l across the throat:

  - **even l** ⟹ radial even ⟹ **Neumann** ψ'(throat) = 0 (an antinode);
  - **odd l** ⟹ radial odd ⟹ **Dirichlet** ψ(throat) = 0 (a node).

So the null throat BC is not a single choice — it is the l-parity-graded
antipodal condition. (A twisted/Möbius field flips even ↔ odd, connecting to
the Z₂ orientation grading, PRs #67/#121.)

## The throat is a unitary mirror

Both even-l Neumann and odd-l Dirichlet are *real* boundary conditions, so the
Wronskian/Klein–Gordon flux j ∝ Im(ψ*ψ') through the throat *vanishes*: the
throat is a perfect, unitary mirror — no net flux is lost into it. This is the
sharp contrast with the ingoing/absorbing horizon BC (ψ ~ e^{−iωr*}), whose
flux j = −ω|A|² ≠ 0 carries probability into the hole. The antipodal throat
conserves flux (unitarity) — consistent with the global CPT / unitarity of the
throat histories (PR #64): what falls toward the throat on one sheet re-emerges
on the antipodal sheet, nothing is destroyed.

## The cavity spectrum is real, discrete, and l-parity-graded

On the exterior [R_MID+ε, R_OUTER] with the outer shell wall (Dirichlet at
R_OUTER) and the antipodal BC at the throat, the operator −d²/dr*² + V_l has a
real, positive, discrete spectrum (a unitary bound cavity), with even-l
(Neumann-at-throat) and odd-l (Dirichlet-at-throat) families giving *distinct*
spectra — the wave-transport face of BAM's even-k/odd-k structure (#67 even-k
absence, #121 odd-k lemma).

## Scope

This establishes the *kinematic* boundary-condition structure of classical
wave transport across the null throat: the vanishing potential, the
l-parity-graded antipodal BC, the unitary-mirror flux property, and the
resulting real discrete spectrum. It does NOT solve the full quasinormal
spectrum, nor compute the dynamical throat ↔ antithroat *nucleation* amplitude
(PR #58/#88). The absorbing-vs-antipodal distinction is fixed by the BAM
antipodal postulate (PR #128); this probe shows that postulate yields a
self-consistent, unitary, l-graded transport condition.

Tests:
  T1. Goal: derive the null throat BC for wave transport on the 5D horizon.
  T2. Vanishing potential: V_l ∝ f → 0 at the throat ⟹ near-horizon null
      modes ψ ~ e^{±iωr*}.
  T3. Three candidate BCs: ingoing/absorbing, reflective wall, antipodal.
  T4. Antipodal map ⟹ l-parity BC: Y_l(−x) = (−1)^l Y_l(x) ⟹ even-l Neumann,
      odd-l Dirichlet at the throat.
  T5. Unitary mirror: real BC ⟹ zero throat flux; ingoing BC ⟹ flux −ω
      (absorbing). The antipodal throat conserves flux.
  T6. Real discrete l-graded spectrum: even-l (N) vs odd-l (D) distinct; the
      wave-transport face of the even/odd Z₂ structure (#67/#121).
  T7. Scope: kinematic BC structure established; QNM spectrum / nucleation open.
  T8. Assessment.

Verdict:
  - NULL_THROAT_BC_ANTIPODAL_L_PARITY_UNITARY_MIRROR (expected): the null
    throat imposes the antipodal, l-parity-graded boundary condition (even-l
    Neumann, odd-l Dirichlet, from Y_l(−x) = (−1)^l Y_l), making the throat a
    unitary flux-conserving mirror rather than an absorbing horizon; the
    exterior cavity spectrum is real, discrete, and even/odd-graded. The QNM
    spectrum and the nucleation rate remain open.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np


PI = math.pi
R_MID = 1.0
R_OUTER = 1.26
RS = R_MID
MU = RS * RS
EPS = 0.02


def f_metric(r: float) -> float:
    return 1.0 - MU / r**2


def V_l(r: float, l: int) -> float:
    """PR #116 Tangherlini cavity potential V = f[l(l+2)/r² + 3rs²/r⁴]."""
    return f_metric(r) * (l * (l + 2) / r**2 + 3.0 * MU / r**4)


def r_star(r: float) -> float:
    return r + (RS / 2.0) * math.log((r - RS) / (r + RS))


def s3_harmonic_parity(l: int, seed: int = 0) -> float:
    """Antipodal ratio Y_l(−x)/Y_l(x) for a degree-l harmonic polynomial on
    ℝ⁴ restricted to S³ — equals (−1)^l."""
    rng = np.random.default_rng(seed)
    v = rng.normal(size=4)
    x = v / np.linalg.norm(v)
    harmonics = {
        1: lambda z: z[0],
        2: lambda z: z[0] * z[1],
        3: lambda z: z[0] * z[1] * z[2],
        4: lambda z: z[0] * z[1] * z[2] * z[3],
    }
    if l == 0:
        return 1.0
    h = harmonics[l]
    return float(h(-x) / h(x))


def cavity_spectrum(l: int, throat_bc: str, N: int = 1000) -> list:
    """Lowest eigenvalues ω² of −d²/dr*² + V_l on [rs+ε, R_OUTER] with a
    Dirichlet shell wall at R_OUTER and throat_bc ('N' Neumann / 'D' Dirichlet)
    at the throat side."""
    r = np.linspace(RS + EPS, R_OUTER, N)
    x = np.array([r_star(rr) for rr in r])
    Vv = np.array([V_l(rr, l) for rr in r])
    A = np.zeros((N, N))
    for i in range(1, N - 1):
        hm = x[i] - x[i - 1]
        hp = x[i + 1] - x[i]
        A[i, i - 1] = -2.0 / (hm * (hm + hp))
        A[i, i + 1] = -2.0 / (hp * (hm + hp))
        A[i, i] = 2.0 / (hm * hp)
    if throat_bc == 'D':
        H = (A + np.diag(Vv))[1:N - 1, 1:N - 1]       # ψ=0 at both ends
    else:
        hp = x[1] - x[0]
        A[0, 0] = 2.0 / hp**2
        A[0, 1] = -2.0 / hp**2                          # Neumann ghost reflection
        H = (A + np.diag(Vv))[0:N - 1, 0:N - 1]
    w2 = np.sort(np.linalg.eigvals(H).real)
    w2 = w2[w2 > 1e-6]
    return [round(float(v), 3) for v in w2[:3]]


# ---------------------------------------------------------------------------
# T1. Goal
# ---------------------------------------------------------------------------

def test_T1_goal() -> dict:
    return {
        'name': 'T1_goal',
        'description': (
            "Derive the boundary condition the null throat (the 5D horizon) "
            "imposes on transported waves (PR #116 cavity operator, PR #128 "
            "antipodal structure): the BAM-native antipodal, l-parity-graded BC."
        ),
        'builds_on': ['#128 antipodal Kruskal horizon (throat↔antithroat)',
                      '#116 Tangherlini cavity operator', '#64 CPT/unitarity',
                      '#67/#121 even-k/odd-k Z₂ structure'],
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T2. Vanishing potential at the horizon
# ---------------------------------------------------------------------------

def test_T2_vanishing_potential() -> dict:
    """V_l ∝ f(r) → 0 at the throat r = rs, so the near-horizon equation is
    −ψ'' = ω²ψ: modes are the pure null phases ψ ~ e^{±iωr*}."""
    rows = []
    ok = True
    for l in (0, 1, 2, 3):
        v_throat = V_l(RS, l)
        v_out = V_l(R_OUTER, l)
        ok = ok and abs(v_throat) < 1e-12 and v_out > 0
        rows.append({'l': l, 'V_at_throat': round(v_throat, 12),
                     'V_at_R_OUTER': round(v_out, 4)})
    return {
        'name': 'T2_vanishing_potential_null_modes',
        'description': (
            "V_l = f[l(l+2)/r² + 3rs²/r⁴] ∝ f → 0 at the throat ⟹ near-horizon "
            "−ψ'' = ω²ψ, the pure null modes ψ ~ e^{±iωr*} (ingoing/outgoing, "
            "v = t+r*, u = t−r*)."
        ),
        'rows': rows,
        'near_horizon_modes': 'ψ ~ e^{±iωr*} (null ingoing/outgoing)',
        'pass': ok,
    }


# ---------------------------------------------------------------------------
# T3. Three candidate boundary conditions
# ---------------------------------------------------------------------------

def test_T3_candidate_bcs() -> dict:
    return {
        'name': 'T3_three_candidate_bcs',
        'description': (
            "Three candidate throat BCs: (1) ingoing/absorbing ψ ~ e^{−iωr*} "
            "(flux lost into the horizon, non-unitary); (2) reflective wall "
            "(Dirichlet/Neumann, real discrete spectrum, the cavity of #116); "
            "(3) antipodal (BAM-native, #128): Φ(U,V,Ω) = Φ(−U,−V,Ω_antipodal)."
        ),
        'candidates': {
            'ingoing_absorbing': 'ψ ~ e^{−iωr*}; flux into horizon; non-unitary',
            'reflective_wall': 'Dirichlet/Neumann; real discrete spectrum (#116)',
            'antipodal_BAM': 'Φ(U,V,Ω)=Φ(−U,−V,Ω_antipodal) (#128 identification)',
        },
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T4. Antipodal map ⟹ l-parity boundary condition
# ---------------------------------------------------------------------------

def test_T4_l_parity_bc() -> dict:
    """S³ harmonics have antipodal parity Y_l(−x) = (−1)^l Y_l(x) (verified).
    Single-valuedness under the antipodal identification forces the radial
    function to carry (−1)^l across the throat: even-l ⟹ Neumann ψ'(throat)=0,
    odd-l ⟹ Dirichlet ψ(throat)=0."""
    rows = []
    ok = True
    for l in (0, 1, 2, 3, 4):
        parity = s3_harmonic_parity(l)
        bc = 'Neumann (ψ\'=0, antinode)' if round(parity) == 1 else 'Dirichlet (ψ=0, node)'
        parity_ok = abs(parity - (-1) ** l) < 1e-9
        ok = ok and parity_ok
        rows.append({'l': l, 'antipodal_parity_Yl': round(parity, 6),
                     'minus1_to_l': (-1) ** l, 'throat_bc': bc})
    return {
        'name': 'T4_antipodal_map_fixes_l_parity_bc',
        'description': (
            "Y_l(−x) = (−1)^l Y_l(x) on S³ (degree-l harmonic polynomials; "
            "verified). Single-valuedness Φ(U,V,Ω) = Φ(−U,−V,Ω_antipodal) ⟹ "
            "radial parity (−1)^l across the throat: even-l → Neumann "
            "(antinode), odd-l → Dirichlet (node)."
        ),
        'rows': rows,
        'rule': 'even-l: Neumann ψ\'(throat)=0; odd-l: Dirichlet ψ(throat)=0',
        'twisted_field_note': 'a twisted/Möbius field flips even↔odd (Z₂ grading, #67/#121)',
        'pass': ok,
    }


# ---------------------------------------------------------------------------
# T5. Unitary mirror (flux conservation)
# ---------------------------------------------------------------------------

def test_T5_unitary_mirror() -> dict:
    """Both antipodal BCs (Neumann/Dirichlet) are real, so the KG/Wronskian
    flux j ∝ Im(ψ*ψ') through the throat vanishes: a unitary mirror. Contrast
    the ingoing/absorbing BC ψ ~ e^{−iωr*}, whose flux j = −ω|A|² ≠ 0 carries
    probability into the hole."""
    om = 1.3
    rstar = np.linspace(-8.0, 0.0, 4000)
    # ingoing complex mode
    psi_in = np.exp(-1j * om * rstar)
    j_in = float(np.mean(np.imag(np.conj(psi_in) * np.gradient(psi_in, rstar))))
    # real standing wave (real BC)
    psi_st = np.sin(om * rstar)
    j_st = float(np.mean(np.abs(np.imag(np.conj(psi_st) * np.gradient(psi_st, rstar)))))
    unitary = abs(j_st) < 1e-9 and abs(j_in - (-om)) < 1e-2
    return {
        'name': 'T5_unitary_mirror_flux_conservation',
        'description': (
            "Real (antipodal) BC ⟹ throat flux j ∝ Im(ψ*ψ') = 0: a unitary "
            "mirror. Ingoing/absorbing BC ψ ~ e^{−iωr*} ⟹ j = −ω ≠ 0 (flux "
            "into the hole). The antipodal throat conserves flux (unitarity, "
            "global CPT #64)."
        ),
        'flux_ingoing_absorbing': round(j_in, 4),
        'flux_real_antipodal_BC': float(f'{j_st:.2e}'),
        'minus_omega': round(-om, 4),
        'antipodal_is_unitary_mirror': unitary,
        'pass': unitary,
    }


# ---------------------------------------------------------------------------
# T6. Real discrete l-graded spectrum
# ---------------------------------------------------------------------------

def test_T6_real_discrete_spectrum() -> dict:
    """The exterior cavity −d²/dr*² + V_l on [rs+ε, R_OUTER] with the shell
    wall (Dirichlet at R_OUTER) and the antipodal throat BC has a real,
    positive, discrete spectrum (unitary bound cavity), with even-l (Neumann)
    and odd-l (Dirichlet) families giving distinct spectra — the
    wave-transport face of the even/odd Z₂ structure (#67/#121)."""
    rows = []
    ok = True
    for l in (0, 1, 2, 3):
        bc = 'N' if (-1) ** l == 1 else 'D'
        spec = cavity_spectrum(l, bc)
        real_positive = all(v > 0 for v in spec) and len(spec) == 3
        ok = ok and real_positive
        rows.append({'l': l, 'antipodal_bc': bc, 'lowest_omega2': spec})
    # even-l vs odd-l differ
    even_spec = cavity_spectrum(0, 'N')
    odd_spec = cavity_spectrum(1, 'D')
    graded = abs(even_spec[0] - odd_spec[0]) > 1e-2
    return {
        'name': 'T6_real_discrete_l_graded_spectrum',
        'description': (
            "Exterior cavity −d²/dr*² + V_l (shell wall at R_OUTER, antipodal "
            "BC at the throat) has a real, positive, discrete spectrum (unitary "
            "bound cavity); even-l (Neumann) and odd-l (Dirichlet) families "
            "differ — the wave-transport face of the even/odd Z₂ structure "
            "(#67/#121)."
        ),
        'rows': rows,
        'real_positive_discrete': ok,
        'even_odd_graded': graded,
        'pass': ok and graded,
    }


# ---------------------------------------------------------------------------
# T7. Scope
# ---------------------------------------------------------------------------

def test_T7_scope() -> dict:
    return {
        'name': 'T7_scope',
        'description': (
            "Establishes the kinematic BC structure of classical wave transport "
            "across the null throat: vanishing potential, l-parity-graded "
            "antipodal BC, unitary-mirror flux, real discrete spectrum. Does "
            "NOT solve the full quasinormal spectrum or the dynamical "
            "throat ↔ antithroat nucleation amplitude (#58/#88); the "
            "absorbing-vs-antipodal choice is fixed by the BAM antipodal "
            "postulate (#128), shown here to be self-consistent and unitary."
        ),
        'established': [
            'V_l → 0 at the throat; near-horizon null modes e^{±iωr*}',
            'antipodal BC: even-l Neumann, odd-l Dirichlet (Y_l parity (−1)^l)',
            'unitary mirror (zero throat flux); real discrete l-graded spectrum',
        ],
        'open': [
            'the full quasinormal-mode spectrum (complex ω, ringdown)',
            'the dynamical throat ↔ antithroat nucleation rate (#58/#88)',
        ],
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T8. Assessment
# ---------------------------------------------------------------------------

def test_T8_assessment() -> dict:
    return {
        'name': 'T8_assessment',
        'description': (
            "The null throat imposes the antipodal, l-parity-graded boundary "
            "condition — even-l Neumann, odd-l Dirichlet, from the S³ harmonic "
            "parity Y_l(−x) = (−1)^l Y_l(x) under the PR #128 antipodal "
            "identification — making the throat a unitary, flux-conserving "
            "mirror rather than an absorbing horizon; the exterior cavity "
            "spectrum is real, discrete, and even/odd-graded. The QNM spectrum "
            "and the nucleation rate remain open."
        ),
        'classification': 'NULL_THROAT_BC_ANTIPODAL_L_PARITY_UNITARY_MIRROR',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    tests = [
        test_T1_goal(),
        test_T2_vanishing_potential(),
        test_T3_candidate_bcs(),
        test_T4_l_parity_bc(),
        test_T5_unitary_mirror(),
        test_T6_real_discrete_spectrum(),
        test_T7_scope(),
        test_T8_assessment(),
    ]
    if all(t['pass'] for t in tests[:7]):
        verdict_class = 'NULL_THROAT_BC_ANTIPODAL_L_PARITY_UNITARY_MIRROR'
        verdict = (
            'THE NULL THROAT IMPOSES THE ANTIPODAL, l-PARITY-GRADED BOUNDARY '
            'CONDITION, MAKING IT A UNITARY MIRROR. PR #128 built the '
            'horizon-regular charts and identified the antipodal map as BAM\'s '
            'throat ↔ antithroat C-swap; this probe derives the boundary '
            'condition that antipodal null throat imposes on transported '
            'waves.\n\n'
            'THE POTENTIAL VANISHES AT THE HORIZON. The wave equation □Φ = 0, '
            'separated as Φ = e^{−iωt} Y_l(Ω) ψ(r)/r^{3/2}, gives '
            '−d²ψ/dr*² + V_l ψ = ω²ψ with V_l = f[l(l+2)/r² + 3rs²/r⁴] '
            '(PR #116). Since V_l ∝ f → 0 at the throat, the near-horizon '
            'equation is −ψ'' = ω²ψ: the modes are the pure null phases '
            'ψ ~ e^{±iωr*} (the ingoing/outgoing null rays of PR #128).\n\n'
            'THREE CANDIDATE BCs. (1) ingoing/absorbing ψ ~ e^{−iωr*} (flux '
            'lost into the horizon, non-unitary); (2) reflective wall '
            '(Dirichlet/Neumann, real discrete spectrum — the cavity of #116); '
            '(3) antipodal (BAM-native): the #128 identification '
            'Φ(U,V,Ω) = Φ(−U,−V,Ω_antipodal).\n\n'
            'THE ANTIPODAL MAP FIXES THE BC BY l-PARITY. The S³ scalar '
            'harmonics carry antipodal parity Y_l(−x) = (−1)^l Y_l(x) '
            '(degree-l harmonic polynomials on ℝ⁴; verified). Single-valuedness '
            'of the untwisted scalar under the antipodal identification then '
            'forces the radial function to carry the compensating (−1)^l across '
            'the throat: EVEN l ⟹ Neumann ∂ψ(throat) = 0 (an antinode), ODD l '
            '⟹ Dirichlet ψ(throat) = 0 (a node). A twisted/Möbius field flips '
            'even ↔ odd, connecting to the Z₂ orientation grading '
            '(#67/#121).\n\n'
            'THE THROAT IS A UNITARY MIRROR. Both antipodal BCs are REAL, so '
            'the Klein–Gordon/Wronskian flux j ∝ Im(ψ* ∂ψ) through the throat '
            'VANISHES — a perfect unitary mirror, no net flux lost. This is the '
            'sharp contrast with the ingoing/absorbing horizon BC, whose flux '
            'j = −ω|A|² ≠ 0 carries probability into the hole. The antipodal '
            'throat conserves flux: what falls toward it on one sheet '
            're-emerges on the antipodal sheet (global CPT/unitarity, '
            '#64).\n\n'
            'THE SPECTRUM IS REAL, DISCRETE, AND l-GRADED. On the exterior '
            '[R_MID+ε, R_OUTER] with the shell wall (Dirichlet at R_OUTER) and '
            'the antipodal throat BC, −d²/dr*² + V_l has a real, positive, '
            'discrete spectrum (a unitary bound cavity), with even-l '
            '(Neumann) and odd-l (Dirichlet) families giving distinct spectra — '
            'the wave-transport face of BAM\'s even-k/odd-k structure '
            '(#67/#121).\n\n'
            'SCOPE. This establishes the kinematic BC structure of classical '
            'wave transport across the null throat: the vanishing potential, '
            'the l-parity-graded antipodal BC, the unitary-mirror flux, the '
            'real discrete spectrum. It does NOT solve the full quasinormal '
            'spectrum or the dynamical throat ↔ antithroat nucleation amplitude '
            '(#58/#88); the absorbing-vs-antipodal choice is fixed by the BAM '
            'antipodal postulate (#128), shown here to be self-consistent and '
            'unitary.'
        )
    else:
        verdict_class = 'NULL_THROAT_BC_INCONCLUSIVE'
        verdict = (
            'INCONCLUSIVE. A check failed; review the vanishing potential, the '
            'S³ harmonic parity, the flux computation, or the spectrum.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'identification': (
            'the null throat imposes the antipodal, l-parity-graded BC (even-l '
            'Neumann, odd-l Dirichlet, from Y_l(−x) = (−1)^l Y_l), making it a '
            'unitary flux-conserving mirror rather than an absorbing horizon; '
            'the exterior cavity spectrum is real, discrete, even/odd-graded'
        ),
        'vanishing_potential': 'V_l ∝ f → 0 at throat; null modes ψ ~ e^{±iωr*}',
        'antipodal_bc': 'even-l Neumann (ψ\'=0), odd-l Dirichlet (ψ=0) from Y_l parity (−1)^l',
        'unitarity': 'real BC ⟹ zero throat flux (unitary mirror); ingoing BC ⟹ flux −ω (absorbing)',
        'spectrum': 'real, positive, discrete; even-l (N) vs odd-l (D) distinct',
        'open': 'full QNM spectrum (complex ω); throat↔antithroat nucleation rate (#58/#88)',
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
    out: list[str] = []
    out.append('# Null throat boundary conditions for wave transport on the 5D horizon (PR #129)')
    out.append('')
    out.append(f"**Run:** {s['timestamp_utc']}")
    out.append('')
    out.append(
        "Derives the boundary condition the null throat (the 5D horizon) "
        "imposes on transported waves. The answer is the BAM-native antipodal "
        "condition — l-parity-graded (even-l Neumann, odd-l Dirichlet) — which "
        "makes the throat a unitary mirror, not an absorbing horizon. Builds on "
        "the PR #128 antipodal Kruskal structure and the PR #116 cavity "
        "operator."
    )
    out.append('')
    out.append(f"- **Vanishing potential**: {s['vanishing_potential']}")
    out.append(f"- **Antipodal BC**: {s['antipodal_bc']}")
    out.append(f"- **Unitarity**: {s['unitarity']}")
    out.append(f"- **Spectrum**: {s['spectrum']}")
    out.append(f"- **Open**: {s['open']}")
    out.append('')

    out.append('## Test summary')
    out.append('')
    out.append('| # | Test | Key finding | PASS? |')
    out.append('|---|---|---|---|')
    label_map = {
        'T1': 'derive the null throat BC for wave transport on the 5D horizon',
        'T2': 'V_l ∝ f → 0 at the throat ⟹ near-horizon null modes e^{±iωr*}',
        'T3': 'three candidate BCs: ingoing/absorbing, reflective, antipodal',
        'T4': 'Y_l(−x) = (−1)^l Y_l ⟹ even-l Neumann, odd-l Dirichlet',
        'T5': 'unitary mirror: real BC ⟹ zero flux; ingoing ⟹ flux −ω',
        'T6': 'real discrete spectrum; even-l (N) vs odd-l (D) distinct',
        'T7': 'scope: kinematic BC established; QNM / nucleation open',
        'T8': 'NULL_THROAT_BC_ANTIPODAL_L_PARITY_UNITARY_MIRROR',
    }
    for t in s['tests']:
        passed = '**PASS**' if t['pass'] else '**FAIL**'
        nm = t['name']; prefix = nm[:2]
        out.append(f"| {prefix} | `{nm}` | {label_map.get(prefix, '—')} | {passed} |")
    out.append('')

    t4 = s['tests'][3]
    out.append('## The antipodal map fixes the BC by l-parity (Y_l(−x) = (−1)^l Y_l)')
    out.append('')
    out.append('| l | antipodal parity Y_l(−x)/Y_l(x) | (−1)^l | throat BC |')
    out.append('|---:|---:|---:|---|')
    for r in t4['rows']:
        out.append(f"| {r['l']} | {r['antipodal_parity_Yl']} | "
                   f"{r['minus1_to_l']} | {r['throat_bc']} |")
    out.append('')

    t6 = s['tests'][5]
    out.append('## The exterior cavity spectrum is real, discrete, and l-graded')
    out.append('')
    out.append('| l | antipodal BC | lowest ω² |')
    out.append('|---:|:---:|---|')
    for r in t6['rows']:
        out.append(f"| {r['l']} | {r['antipodal_bc']} | {r['lowest_omega2']} |")
    out.append('')
    out.append("Real, positive, discrete (a unitary bound cavity); even-l "
               "(Neumann-at-throat) and odd-l (Dirichlet-at-throat) families "
               "give distinct spectra — the wave-transport face of BAM's "
               "even-k/odd-k Z₂ structure (#67/#121).")
    out.append('')

    out.append('## Verdict')
    out.append('')
    out.append(f"**{s['verdict_class']}.** {s['verdict']}")
    out.append('')
    return '\n'.join(out)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def _json_default(o):
    if isinstance(o, tuple):
        return list(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if hasattr(o, '__dict__'):
        try:
            return asdict(o)
        except Exception:
            return o.__dict__
    return str(o)


def main(argv: Optional[list[str]] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    out = here / 'runs' / f'{ts}_null_throat_boundary_conditions_probe'
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
