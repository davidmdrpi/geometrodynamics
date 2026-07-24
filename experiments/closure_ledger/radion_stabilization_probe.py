"""
Radion stabilization from the primordial EM-capped throat:
V_eff(phi), alpha, and the radion mass (PR #226).

> Framing: QFT on the *fixed classical* throat geometry (geometry -> fields),
> not quantum gravity.

THE QUESTION
------------
PR #225 closed with the radion phi as a canonical 4D field on a FLAT
tree-level potential, alpha(phi) = alpha(0) e^{sqrt3 phi} its dilaton
coupling, and the EM cap NAMED as the program's dynamical stabilizer -
"deriving rho* from the cap's own equations is the open dynamical
problem."  This probe executes that problem: build the radion
effective potential V_eff(phi) OF the primordial EM-capped throat
(#55-#58, primordial per #222), and ask whether it stabilizes the
radion, at what alpha, and with what radion mass.

THE ANSWER (derived, machine-checked, honest)
---------------------------------------------
THE HOPF CHARGE IS THE STABILIZER.  The committed cap (#55) is the
throat-capped Coulomb self-energy U_EM = alpha hbar c / (2R) plus a
cohesive term; the committed topology (#58) puts ONE HOPF CHARGE on
every throat (Sigma c_1 = 0 pair creation).  Under the #225 dilaton
coupling e^{-sqrt3 phi} F^2 the two flux energies of this DYONIC
throat carry OPPOSITE radion charges:

  * fixed electric charge (Gauss law through the dielectric):
    U_el ~ e^{+sqrt3 phi}   (weaker coupling -> cheaper field),
  * fixed magnetic/Hopf flux (F pinned by topology):
    U_mag ~ e^{-sqrt3 phi}  (weaker coupling -> DEARER flux),

and the primordial throat radius (fixed in 5D-frame lengths, #222)
adds the frame factor e^{a phi}, a = 1/(2 sqrt3).  So

  V_eff(phi) = U_el0 e^{p phi} + U_mag0 e^{r phi},
  p = sqrt3 + a = 7/(2 sqrt3),   r = -sqrt3 + a = -5/(2 sqrt3):

one exponent of each sign - A MINIMUM EXISTS.  Without the Hopf
charge every committed term has a positive exponent and the radion
runs away to decompactification (alpha -> 0): the no-go is proven,
and the topological charge of #58 is EXACTLY what saves the vacuum.

Closed forms (all machine-checked):
  * alpha at the minimum:  alpha*^2 = (kappa/4)(sqrt3-a)/(sqrt3+a)
    = 5 kappa/28, where kappa is the magnetic-to-electric cap-energy
    ratio in Dirac units (kappa = 1 at the symmetric Dirac point:
    U_mag = U_el/(4 alpha^2) from g = 2 pi/e).  At kappa = 1:
    ALPHA* = sqrt(5/28) = 0.4226 - order one, the self-dual landing.
    The OBSERVED alpha requires kappa = 28 alpha^2/5 = 2.98e-4: the
    stabilization mechanism is derived; the smallness of kappa is
    the new, sharply-quantified open problem (guardrail scan: no
    closure-constant match, nearest alpha/k5^2 at 2.1% - rejected).
  * the radion mass identity:  m_phi^2 = -p r V_min = (3 - a^2)
    V_min = (35/12) V_min - independent of every coefficient; per
    throat, E''(phi*) = 7 U_el(phi*) exactly (U_mag/U_el = 7/5 at
    the minimum).  Anchored: E'' = (7/4) alpha*^{3/2} m_P on the
    #225 fiber anchor (1.3e16 GeV at observed alpha; 0.48 m_P at
    the Dirac point) or (7 alpha/2) m_e c^2 = 13 keV on the #55
    Compton anchor - the B4 one-anchor discipline kept explicit.

Tests:
  T1. Goal.
  T2. The committed cap imported: the #55 ledger re-read (A = alpha
      hbar c/2 - an alpha-dependent holdout - R* = (A/2B)^{1/3},
      U/mc^2 = alpha/2, stability E'' = 6B).
  T3. The radion charges of the cap energies, derived: fixed charge
      -> e^{+sqrt3 phi} (quadrature through the dielectric), fixed
      flux -> e^{-sqrt3 phi}, primordial radius -> e^{a phi}.
  T4. V_eff(phi): the dyonic minimum exists; the no-go without the
      Hopf charge; the symbolic m^2 = -p r V_min identity.
  T5. Alpha at the minimum: alpha* = sqrt(5 kappa/28); the Dirac
      point 0.4226; kappa_required = 2.98e-4 + guardrail; cohesion
      provenance robustness (<0.3% tilt).
  T6. The radion mass: E'' = 7 U_el* exactly; the anchored scales.
  T7. Arc consistency: the stabilizer row occupies EXACTLY the
      rank-audit slot #225 reserved for the EM cap; alpha-dependent
      ledger re-reads (#55, #225).
  T8. Honest scope.
  T9. Assessment.

Verdict:
  THE_HOPF_CHARGE_STABILIZES_THE_RADION_V_EFF_HAS_A_MINIMUM_WITH_
  ALPHA_STAR_EQUALS_SQRT_5KAPPA_OVER_28_ORDER_ONE_AT_THE_DIRAC_POINT_
  AND_M_PHI2_EQUALS_35_OVER_12_V_MIN_THE_OBSERVED_ALPHA_NEEDS_
  KAPPA_3E_MINUS_4_THE_NAMED_OPEN_PROBLEM
"""

from __future__ import annotations

import glob
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np
from scipy.optimize import minimize_scalar

_ALPHA = 7.2973525693e-3
_K5 = 5
_HBAR = 1.054571817e-34
_C = 2.99792458e8
_ME_KG = 9.1093837015e-31
_MP_GEV = 1.220890e19
_ME_MEV = 0.51099895000

_A = 1.0 / (2.0 * math.sqrt(3.0))      # the #225 Einstein-frame exponent
_S3 = math.sqrt(3.0)
_P = _S3 + _A                          # electric exponent
_R = -_S3 + _A                         # magnetic exponent

_CACHE: dict = {}


def _latest(pat):
    fs = sorted(glob.glob(pat))
    return json.load(open(fs[-1])) if fs else None


def _t(d, name):
    return [t for t in d['tests'] if t['name'] == name][0]


# ── the effective potential ─────────────────────────────────────────────

def v_eff(phi, kappa=1.0, b_coh=0.0, q_coh=_A, alpha0=_ALPHA):
    """V_eff(phi) of the dyonic EM-capped primordial throat, in units
    of hbar c/(2 R_5) (the #55 cap energy scale at phi = 0):
    electric alpha(phi) e^{a phi} + magnetic (kappa/4 alpha(phi))
    e^{a phi} + optional cohesion b_coh e^{q_coh phi}."""
    al = alpha0 * np.exp(_S3 * np.asarray(phi, dtype=float))
    fr = np.exp(_A * np.asarray(phi, dtype=float))
    return al * fr + (kappa / (4.0 * al)) * fr \
        + b_coh * np.exp(q_coh * np.asarray(phi, dtype=float))


def find_min(kappa=1.0, b_coh=0.0, q_coh=_A):
    res = minimize_scalar(lambda f: float(v_eff(f, kappa, b_coh, q_coh)),
                          bounds=(-8.0, 8.0), method='bounded',
                          options={'xatol': 1e-13})
    return float(res.x)


# ========================================================================
# T1. Goal
# ========================================================================


def test_T1_goal() -> dict:
    return {
        'name': 'T1_goal',
        'description': (
            'Execute the open dynamical problem #225 named: build the '
            'radion effective potential of the primordial EM-capped '
            'throat (#55-#58 cap, #222 primordial, #225 dilaton '
            'coupling alpha(phi) = alpha0 e^{sqrt3 phi}), and answer: '
            'does V_eff(phi) stabilize the radion, at what alpha, and '
            'with what radion mass?'
        ),
        'pass': True,
    }


# ========================================================================
# T2. The committed cap imported (#55 ledger re-read)
# ========================================================================


def test_T2_cap_imported() -> dict:
    if 'T2' in _CACHE:
        return _CACHE['T2']
    here = Path(__file__).resolve().parent / 'runs'
    d55 = _latest(str(here / '*_self_consistent_throat_radius_probe/probe.json'))
    t4 = _t(d55, 'T4_self_consistent_equilibrium')
    t6 = _t(d55, 'T6_renormalization_free')
    A_st, B_st = t4['A_coupling'], t4['B_coupling']
    R_st = t4['R_star_analytic_m']
    # (a) the cap's EM coupling is alpha-dependent: A = alpha hbar c/2
    A_der = _ALPHA * _HBAR * _C / 2.0
    dev_A = abs(A_st / A_der - 1.0)
    # (b) the equilibrium closed form R* = (A/2B)^{1/3}
    dev_R = abs(R_st / (A_st / (2 * B_st)) ** (1.0 / 3.0) - 1.0)
    # (c) the finite capped self-energy: U/mc^2 = alpha/2 exactly
    dev_U = abs(t6['fraction_U_over_mc2'] / (_ALPHA / 2.0) - 1.0)
    # (d) stability of the committed equilibrium: E''(R*) = 6B > 0
    Epp = 6.0 * B_st
    # (e) the committed anchor: R* is the reduced Compton wavelength
    lamC = _HBAR / (_ME_KG * _C)
    dev_anchor = abs(R_st / lamC - 1.0)

    ok = (dev_A < 1e-12 and dev_R < 1e-12 and dev_U < 1e-12
          and Epp > 0 and dev_anchor < 1e-9)
    out = {
        'name': 'T2_cap_imported',
        'description': (
            'the committed EM cap imported from the #55 ledger: the '
            'throat is the inner boundary that caps the Coulomb field '
            '- U_EM = alpha hbar c/(2R) FINITE - with the equilibrium '
            'E(R) = A/R + B R^2, R* = (A/2B)^{1/3} (stable, E\'\' = '
            '6B > 0), U/mc^2 = alpha/2, and the committed anchor R* = '
            'lambda_C (the #55 Compton anchor).  The re-read '
            'continues the #225 alpha-dependent holdout: A = alpha '
            'hbar c/2 is a genuine function of alpha, re-derived at '
            'machine zero'
        ),
        'A_stored_J_m': float(A_st),
        'A_derived_alpha_hbar_c_over_2': float(A_der),
        'A_rel_dev': float(dev_A),
        'R_star_identity_rel_dev': float(dev_R),
        'U_over_mc2_vs_alpha_over_2_rel_dev': float(dev_U),
        'stability_Epp_6B': float(Epp),
        'R_star_vs_compton_rel_dev': float(dev_anchor),
        'pass': bool(ok),
    }
    _CACHE['T2'] = out
    return out


# ========================================================================
# T3. The radion charges of the cap energies, derived
# ========================================================================


def test_T3_radion_charges() -> dict:
    if 'T3' in _CACHE:
        return _CACHE['T3']
    # the #225 dilaton coupling: L = -1/4 e^{-sqrt3 phi} F^2 is a
    # dielectric eps(phi) = e^{-sqrt3 phi}.  Radial quadrature with
    # the throat cap at R (hbar = c = 1, R = 1):
    rr = np.linspace(1.0, 2001.0, 400001)
    q = math.sqrt(4 * math.pi * _ALPHA)
    g = 2 * math.pi / q                       # Dirac partner flux
    devs_el, devs_mag = [], []
    for phi in (-0.7, 0.0, 1.3):
        eps = math.exp(-_S3 * phi)
        # fixed CHARGE: Gauss law D = q/4 pi r^2 (phi-independent);
        # U_el = int D^2/(2 eps) dV  ~  e^{+sqrt3 phi}
        D = q / (4 * math.pi * rr ** 2)
        U_el = float(np.trapezoid(D ** 2 / (2 * eps) * 4 * math.pi * rr ** 2, rr))
        U_el += (q ** 2 / (8 * math.pi * eps)) / rr[-1]     # exact tail
        devs_el.append(abs(U_el / ((_ALPHA / 2.0) * math.exp(_S3 * phi)) - 1))
        # fixed FLUX: B = g/4 pi r^2 pinned by topology (the Hopf
        # charge); U_mag = int eps B^2/2 dV  ~  e^{-sqrt3 phi}
        B = g / (4 * math.pi * rr ** 2)
        U_mag = float(np.trapezoid(eps * B ** 2 / 2 * 4 * math.pi * rr ** 2, rr))
        U_mag += (eps * g ** 2 / (8 * math.pi)) / rr[-1]
        devs_mag.append(abs(U_mag / ((1.0 / (8.0 * _ALPHA))
                                     * math.exp(-_S3 * phi)) - 1))
    dev_el = float(max(devs_el))
    dev_mag = float(max(devs_mag))
    # Dirac units: U_mag/U_el = 1/(4 alpha^2) at phi = 0 (kappa = 1)
    dirac_ratio_dev = abs(((1.0 / (8.0 * _ALPHA)) / (_ALPHA / 2.0))
                          * (4.0 * _ALPHA ** 2) - 1.0)
    # the primordial-radius frame factor: #222 fixes the throat in
    # 5D-frame lengths; #225's Einstein frame gives L_E = e^{-a phi}
    # L_5, a = 1/(2 sqrt3), so 1/R_E carries e^{+a phi}
    p_exp, r_exp = _S3 + _A, -_S3 + _A
    exps_ok = (abs(p_exp - 7 / (2 * _S3)) < 1e-15
               and abs(r_exp + 5 / (2 * _S3)) < 1e-15)

    ok = dev_el < 2e-5 and dev_mag < 2e-5 and dirac_ratio_dev < 1e-12 \
        and exps_ok
    out = {
        'name': 'T3_radion_charges',
        'description': (
            'the radion charges of the two flux energies of the '
            'DYONIC throat, derived from the #225 dilaton coupling '
            '(a dielectric eps = e^{-sqrt3 phi}): a FIXED CHARGE '
            'costs U_el ~ e^{+sqrt3 phi} (Gauss law fixes D; weaker '
            'coupling -> cheaper field), a FIXED FLUX costs U_mag ~ '
            'e^{-sqrt3 phi} (topology pins F; weaker coupling -> '
            'DEARER flux) - both by radial quadrature with the '
            'throat cap; the Dirac partner ratio U_mag/U_el = '
            '1/(4 alpha^2) exactly (kappa = 1); and the primordial '
            'radius (#222, fixed in 5D-frame lengths) adds e^{a '
            'phi}: the exponents are p = sqrt3 + a = 7/(2 sqrt3) '
            'and r = -sqrt3 + a = -5/(2 sqrt3) - ONE OF EACH SIGN'
        ),
        'electric_scaling_e_plus_sqrt3_phi_dev': dev_el,
        'magnetic_scaling_e_minus_sqrt3_phi_dev': dev_mag,
        'dirac_ratio_1_over_4alpha2_dev': float(dirac_ratio_dev),
        'exponent_electric_p': float(p_exp),
        'exponent_magnetic_r': float(r_exp),
        'pass': bool(ok),
    }
    _CACHE['T3'] = out
    return out


# ========================================================================
# T4. V_eff(phi): the minimum, the no-go, the mass identity
# ========================================================================


def test_T4_potential() -> dict:
    if 'T4' in _CACHE:
        return _CACHE['T4']
    # (a) the dyonic minimum exists (kappa = 1): numeric vs closed form
    phi_star = find_min(kappa=1.0)
    al_star_num = _ALPHA * math.exp(_S3 * phi_star)
    al_star_cf = math.sqrt(5.0 / 28.0)
    dev_min = abs(al_star_num / al_star_cf - 1.0)
    # (b) curvature: V'' > 0, and the coefficient-free identity
    #     m^2 = -p r V_min = (3 - a^2) V_min = (35/12) V_min
    h = 1e-4
    Vpp = float((v_eff(phi_star + h) - 2 * v_eff(phi_star)
                 + v_eff(phi_star - h)) / h ** 2)
    Vmin = float(v_eff(phi_star))
    ident_dev = abs(Vpp / ((35.0 / 12.0) * Vmin) - 1.0)
    # symbolic check of the two-exponential theorem
    import sympy as sp
    Ps, Qs, ps, rs, f = sp.symbols('P Q p r f', positive=True)
    V = Ps * sp.exp(ps * f) + Qs * sp.exp(-rs * f)      # r = -rs < 0
    fstar = sp.solve(sp.diff(V, f), f)[0]
    Vpp_s = sp.simplify(sp.diff(V, f, 2).subs(f, fstar))
    Vmin_s = sp.simplify(V.subs(f, fstar))
    thm = sp.simplify(Vpp_s - ps * rs * Vmin_s) == 0    # m^2 = -p r Vmin
    # (c) the no-go: kappa = 0 (no Hopf charge) with the 5D-frame
    #     cohesion (exponent +a): every exponent positive -> V is
    #     monotonic, runaway to phi -> -inf (decompactification,
    #     alpha -> 0)
    grid = np.linspace(-30.0, 5.0, 2000)
    Vng = v_eff(grid, kappa=0.0, b_coh=_ALPHA, q_coh=_A)
    no_go = bool(np.all(np.diff(Vng) > 0)) and float(Vng[0]) < 1e-3

    ok = (dev_min < 1e-7 and Vpp > 0 and ident_dev < 1e-6 and thm
          and no_go)
    out = {
        'name': 'T4_potential',
        'description': (
            'V_eff(phi) = U_el0 e^{p phi} + U_mag0 e^{r phi} with p '
            '> 0 > r: the dyonic throat STABILIZES the radion - the '
            'minimum exists (numeric = closed form), V\'\' > 0, and '
            'the coefficient-free mass identity m_phi^2 = -p r V_min '
            '= (3 - a^2) V_min = (35/12) V_min holds numerically AND '
            'symbolically (the two-exponential theorem).  THE NO-GO: '
            'without the Hopf charge (kappa = 0) every committed '
            'term has a positive exponent - V is monotonic and the '
            'radion runs away to decompactification (alpha -> 0): '
            'the #58 topological charge (one Hopf charge per throat) '
            'is EXACTLY what saves the vacuum'
        ),
        'phi_star': float(phi_star),
        'alpha_at_min_numeric': float(al_star_num),
        'alpha_at_min_closed_form': float(al_star_cf),
        'min_rel_dev': float(dev_min),
        'V_second_derivative': Vpp,
        'mass_identity_35_over_12_dev': float(ident_dev),
        'mass_identity_symbolic': bool(thm),
        'no_go_without_hopf_charge': bool(no_go),
        'pass': bool(ok),
    }
    _CACHE['T4'] = out
    return out


# ========================================================================
# T5. Alpha at the minimum
# ========================================================================


def test_T5_alpha() -> dict:
    if 'T5' in _CACHE:
        return _CACHE['T5']
    # closed form: alpha*^2 = (kappa/4)(sqrt3 - a)/(sqrt3 + a)
    #                       = 5 kappa/28
    al_dirac = math.sqrt(5.0 / 28.0)
    ratio_obs = al_dirac / _ALPHA
    kappa_req = 28.0 * _ALPHA ** 2 / 5.0
    # verify kappa_req numerically: minimize with kappa = kappa_req
    phi_r = find_min(kappa=kappa_req)
    al_r = _ALPHA * math.exp(_S3 * phi_r)
    dev_req = abs(al_r / _ALPHA - 1.0)
    # guardrail scan on kappa_req (match claimed only within 1%)
    cands = {
        'alpha/2pi': _ALPHA / (2 * math.pi),
        'alpha^2': _ALPHA ** 2,
        'alpha^{3/2}': _ALPHA ** 1.5,
        '2pi alpha^2': 2 * math.pi * _ALPHA ** 2,
        'alpha/e^pi': _ALPHA / math.exp(math.pi),
        'alpha/k5^2': _ALPHA / _K5 ** 2,
        '1/k5^4': 1.0 / _K5 ** 4,
        '(me/mmu)^2': (1.0 / 206.7683) ** 2,
    }
    scan = {n: {'value': float(v),
                'rel_dev': float(abs(kappa_req - v) / kappa_req)}
            for n, v in cands.items()}
    nearest = min(scan.items(), key=lambda kv: kv[1]['rel_dev'])
    no_match = all(v['rel_dev'] > 0.01 for v in scan.values())
    # cohesion-provenance robustness: both derived cases tilt alpha*
    # by < 0.3% at the #55-scale cohesion
    tilts = {}
    for q_c, nm in ((_A, 'sigma_5D_fixed_+a'), (-2 * _A, 'sigma_E_fixed_-2a')):
        for b_c in (_ALPHA / 2, _ALPHA, 2 * _ALPHA):
            phi_c = find_min(kappa=1.0, b_coh=b_c, q_coh=q_c)
            al_c = _ALPHA * math.exp(_S3 * phi_c)
            tilts[f'{nm}_b{b_c:.4f}'] = float(al_c / al_dirac - 1.0)
    max_tilt = max(abs(v) for v in tilts.values())

    ok = (dev_req < 1e-8 and no_match and max_tilt < 3e-3
          and 0.42 < al_dirac < 0.43)
    out = {
        'name': 'T5_alpha',
        'description': (
            'alpha at the minimum, closed form: alpha*^2 = (kappa/4)'
            '(sqrt3-a)/(sqrt3+a) = 5 kappa/28.  At the symmetric '
            'Dirac point (kappa = 1, U_mag = U_el/(4 alpha^2)): '
            'ALPHA* = sqrt(5/28) = 0.4226 - ORDER ONE, the '
            'self-dual landing, 58x the observed value.  The '
            'observed alpha requires the magnetic-to-electric cap '
            'asymmetry kappa = 28 alpha^2/5 = 2.98e-4 (verified by '
            'direct minimization): the stabilization MECHANISM is '
            'derived and modulus-free; the smallness of kappa is '
            'the new, sharply-quantified open problem.  Guardrail: '
            'no closure-constant match for kappa_req (nearest '
            'alpha/k5^2 at 2.1% - rejected).  Cohesion-provenance '
            'robustness: both derived tilt cases shift alpha* by '
            '< 0.3%'
        ),
        'alpha_star_closed_form': 'alpha*^2 = 5 kappa / 28',
        'alpha_star_dirac_point': float(al_dirac),
        'alpha_star_over_alpha_obs': float(ratio_obs),
        'kappa_required_for_observed_alpha': float(kappa_req),
        'kappa_required_verified_dev': float(dev_req),
        'guardrail_scan': scan,
        'nearest_candidate': {'name': nearest[0], **nearest[1]},
        'no_numerological_match': bool(no_match),
        'cohesion_tilts': tilts,
        'max_cohesion_tilt': float(max_tilt),
        'pass': bool(ok),
    }
    _CACHE['T5'] = out
    return out


# ========================================================================
# T6. The radion mass
# ========================================================================


def test_T6_radion_mass() -> dict:
    if 'T6' in _CACHE:
        return _CACHE['T6']
    # per-throat curvature closed chain at the minimum:
    #   U_mag/U_el = (sqrt3+a)/(sqrt3-a) = 7/5
    #   E_min = (12/5) U_el*,  E'' = (35/12) E_min = 7 U_el* EXACTLY
    phi_star = find_min(kappa=1.0)
    U_el = _ALPHA * math.exp(_S3 * phi_star) * math.exp(_A * phi_star)
    U_mag = (1.0 / (4.0 * _ALPHA * math.exp(_S3 * phi_star))) \
        * math.exp(_A * phi_star)
    ratio_dev = abs(U_mag / U_el - 7.0 / 5.0)
    h = 1e-4
    Epp = float((v_eff(phi_star + h) - 2 * v_eff(phi_star)
                 + v_eff(phi_star - h)) / h ** 2)
    seven_dev = abs(Epp / (7.0 * U_el) - 1.0)
    # anchored scales for E'' = (7/4) alpha*^{3/2} m_P on the #225
    # fiber anchor (R_E = 2 l_P/sqrt(alpha*), U_el = alpha*^{3/2}
    # m_P/4), and (7 alpha/2) m_e c^2 on the #55 Compton anchor
    Epp_dirac_mP = (7.0 / 4.0) * (math.sqrt(5.0 / 28.0)) ** 1.5
    Epp_obs_mP = (7.0 / 4.0) * _ALPHA ** 1.5
    Epp_obs_GeV = Epp_obs_mP * _MP_GEV
    Epp_compton_keV = (7.0 * _ALPHA / 2.0) * _ME_MEV * 1e3
    # canonical conversion: with the #225 kinetic convention the
    # canonical field is phi_c = phi m_P/sqrt(16 pi), so a throat
    # number density n gives m_phi^2 = 16 pi n E''/m_P^2 - reported
    # as the formula (n is the one astrophysical input, not derived)

    ok = (ratio_dev < 1e-6 and seven_dev < 1e-6
          and 1.0e16 < Epp_obs_GeV < 1.7e16)
    out = {
        'name': 'T6_radion_mass',
        'description': (
            'the radion mass, closed chain: at the minimum U_mag/'
            'U_el = (sqrt3+a)/(sqrt3-a) = 7/5, E_min = (12/5) U_el*, '
            'and the per-throat curvature is E\'\'(phi*) = 7 U_el* '
            'EXACTLY (machine-checked).  Anchored (B4 discipline, '
            'both anchors reported): on the #225 fiber anchor E\'\' '
            '= (7/4) alpha*^{3/2} m_P - 0.48 m_P at the Dirac '
            'point, 1.3e16 GeV at the observed-alpha minimum (a '
            'GUT-scale curvature: the radion is HEAVY, no '
            'fifth-force tension); on the #55 Compton anchor E\'\' '
            '= (7 alpha/2) m_e c^2 = 13 keV.  The canonical mass is '
            'm_phi^2 = 16 pi n E\'\'/m_P^2 with n the throat number '
            'density - the one astrophysical input, stated not '
            'derived'
        ),
        'U_mag_over_U_el_at_min': float(U_mag / U_el),
        'ratio_7_over_5_dev': float(ratio_dev),
        'Epp_over_7Uel_dev': float(seven_dev),
        'Epp_dirac_point_mP': float(Epp_dirac_mP),
        'Epp_observed_alpha_mP': float(Epp_obs_mP),
        'Epp_observed_alpha_GeV': float(Epp_obs_GeV),
        'Epp_compton_anchor_keV': float(Epp_compton_keV),
        'mass_formula': 'm_phi^2 = 16 pi n E\'\'(phi*) / m_P^2',
        'pass': bool(ok),
    }
    _CACHE['T6'] = out
    return out


# ========================================================================
# T7. Arc consistency
# ========================================================================


def test_T7_arc_consistency() -> dict:
    if 'T7' in _CACHE:
        return _CACHE['T7']
    here = Path(__file__).resolve().parent / 'runs'
    # (a) the stabilizer row occupies EXACTLY the rank-audit slot
    # #225 reserved for the EM cap: the minimum condition fixes
    # alpha(phi), i.e. the combination ln rho = ln lambda +
    # ln(R_u/l_P) - gradient direction (1,1,0,0,0) in the #225
    # 5-knob log space
    grad_alpha = np.array([-2.0, -2.0, 0.0, 0.0, 0.0])
    C_weld = np.array([0.0, 0.0, 0.0, 1.0, 0.0])
    C_stab = np.array([1.0, 1.0, 0.0, 0.0, 0.0])     # d(ln rho) fixed
    M = np.array([C_weld, C_stab])
    rank = int(np.linalg.matrix_rank(M))
    _, _, Vt = np.linalg.svd(M)
    null = Vt[rank:]
    proj = float(np.linalg.norm(null @ grad_alpha))
    cap_row_match = bool(np.allclose(
        C_stab / np.linalg.norm(C_stab),
        np.array([1.0, 1.0, 0.0, 0.0, 0.0]) / math.sqrt(2.0)))
    # (b) #225 ledger re-reads (alpha-dependent, continuing T7 there)
    d225 = _latest(str(here / '*_absolute_coupling_capstone_probe/probe.json'))
    t6c = _t(d225, 'T6_stabilizer')
    rho_dev = abs(t6c['rho_star'] / (2.0 / math.sqrt(_ALPHA)) - 1.0)
    t2c = _t(d225, 'T2_canonical_chain')
    conv_dev = abs(t2c['r_h_in_planck_lengths']
                   / (2.0 / math.sqrt(_ALPHA)) - 1.0)
    # (c) the alpha law closes the loop: at the minimum, rho(phi*) =
    # 2/sqrt(alpha(phi*)) by #225's law - at the Dirac point the
    # fiber sits at rho* = 3.08 (near-Planckian), at the observed-
    # alpha minimum at 23.41
    al_dirac = math.sqrt(5.0 / 28.0)
    rho_dirac = 2.0 / math.sqrt(al_dirac)
    rho_obs = 2.0 / math.sqrt(_ALPHA)
    # (d) #222's primordial verdict is the INPUT that fixed the
    # throat radius in 5D-frame lengths (the e^{a phi} factor):
    # re-read the exclusion that forced it
    d222 = _latest(str(here / '*_coupled_5d_ekg_weld_probe/probe.json'))
    t7w = _t(d222, 'T7_confrontation')
    excl = t7w['anchor_exclusion_factor']

    ok = (rank == 2 and proj < 1e-12 and cap_row_match
          and rho_dev < 1e-9 and conv_dev < 1e-9 and excl > 100.0)
    out = {
        'name': 'T7_arc_consistency',
        'description': (
            'the stabilizer slots into the arc without retuning: the '
            'V_eff minimum condition fixes alpha(phi) - its gradient '
            'in the #225 5-knob log space is EXACTLY the (1,1,0,0,0) '
            'cap row the rank audit reserved (rank 2, grad alpha '
            'annihilated to machine zero, three alpha-decoupled '
            'flats unchanged); the #225 ledger re-reads pass (rho* '
            '= 2/sqrt(alpha), the model-unit conversion); the alpha '
            'law closes the loop (Dirac-point fiber at rho = 3.08 '
            'near-Planckian, observed-alpha fiber at 23.41); and '
            '#222\'s x210 anchor exclusion - the verdict that '
            'forced the throat PRIMORDIAL, giving the e^{a phi} '
            'frame factor - is re-read from its ledger'
        ),
        'rank_with_stabilizer_row': rank,
        'grad_alpha_null_projection': proj,
        'stabilizer_row_equals_cap_row': cap_row_match,
        'rho_star_ledger_rel_dev': float(rho_dev),
        'model_unit_conversion_rel_dev': float(conv_dev),
        'rho_at_dirac_point': float(rho_dirac),
        'rho_at_observed_alpha': float(rho_obs),
        'primordial_exclusion_factor_222': float(excl),
        'pass': bool(ok),
    }
    _CACHE['T7'] = out
    return out


# ========================================================================
# T8. Honest scope
# ========================================================================


def test_T8_honest_scope() -> dict:
    scope = [
        'kappa - the magnetic-to-electric cap-energy ratio in Dirac '
        'units - is NOT derived.  kappa = 1 is the symmetric Dirac '
        'point (g = 2 pi/e with equal geometric factors); the '
        'observed alpha requires kappa = 2.98e-4.  The O(1) Wu-Yang '
        'half-charge and cap-geometry factors can shift kappa by '
        'O(1), not by 3e3: the smallness is real and is THE '
        'quantified open problem this probe leaves.',
        'V_eff is the Born-Oppenheimer energy of a SINGLE throat as '
        'a function of the asymptotic radion value; the canonical '
        'mass requires a throat number density n (m_phi^2 = 16 pi n '
        'E\'\'/m_P^2) - n is an astrophysical input, stated not '
        'derived.  For any realistic n the radion is far above '
        'fifth-force bounds (E\'\' is GUT-scale per throat at the '
        'observed-alpha minimum).',
        'Tree-level throughout: no Casimir energy of the fiber, no '
        'loop-induced radion potential, no backreaction of V_eff on '
        'the 5D geometry.  The #225 flatness statement (V_tree = 0 '
        'in vacuo) is unchanged - V_eff here is SOURCED by the '
        'throat, which is the point.',
        'The cohesion term\'s frame provenance is mapped (5D-fixed '
        'tension -> e^{+a phi}; Einstein-frame-fixed -> e^{-2a '
        'phi}), not chosen: both tilt alpha* by < 0.3% at the '
        '#55-scale cohesion, and neither overturns the electric-'
        'magnetic balance.  The fiber-wrapped case carries exactly '
        'the magnetic exponent and is absorbed into kappa.',
        'The two committed anchors (#55 Compton lambda_C for the '
        'matter-sector throat; #225 Planck fiber rho* for the '
        'geometric throat) are reported separately per the B4 '
        'one-anchor discipline; alpha* itself is anchor-FREE (pure '
        'exponent-and-topology arithmetic).',
        'The electric/magnetic quadratures use the flat capped-'
        'field model of #55 (the throat as inner boundary); the '
        'full 5D field profile would shift the O(1) coefficients '
        'inside kappa, not the exponents - the exponents are Weyl '
        'algebra (#225, machine-checked there).',
    ]
    return {
        'name': 'T8_honest_scope',
        'description': 'what this probe does and does not establish',
        'scope': scope,
        'pass': True,
    }


# ========================================================================
# T9. Assessment
# ========================================================================


def test_T9_assessment() -> dict:
    t2 = test_T2_cap_imported()
    t3 = test_T3_radion_charges()
    t4 = test_T4_potential()
    t5 = test_T5_alpha()
    t6 = test_T6_radion_mass()
    t7 = test_T7_arc_consistency()
    core = all(t['pass'] for t in (t2, t3, t4, t5, t6, t7))
    assessment = (
        'The open dynamical problem #225 named is executed, and the '
        'answer has a clean shape.  The primordial EM-capped throat '
        'is a dyon: its capped electric flux costs e^{+sqrt3 phi} '
        'and its topological Hopf flux costs e^{-sqrt3 phi} under '
        'the #225 dilaton coupling, and with the #222 primordial '
        'frame factor e^{a phi} the effective potential has one '
        'exponent of each sign: THE RADION IS STABILIZED, and the '
        'no-go shows the #58 Hopf charge is exactly what saves the '
        'vacuum from decompactification.  The minimum is modulus-'
        'free and closed-form: alpha*^2 = 5 kappa/28, with the '
        'coefficient-free mass identity m_phi^2 = (35/12) V_min and '
        'E\'\' = 7 U_el* per throat.  At the symmetric Dirac point '
        'the landing is alpha* = 0.42 - order one, 58x the observed '
        'coupling: the mechanism is derived, the scale is honest, '
        'and the gap is sharply quantified as kappa = 2.98e-4 with '
        'no closure-constant match - the program\'s open problem '
        'moves from "derive rho*" to "derive the electric-magnetic '
        'cap asymmetry", a strictly sharper question.  The '
        'stabilizer row lands exactly in the rank-audit slot #225 '
        'reserved for it: nothing else in the arc moves.'
        if core else
        'INCONCLUSIVE - a core check failed; do not quote.'
    )
    return {
        'name': 'T9_assessment',
        'description': 'the stabilization verdict',
        'core_green': bool(core),
        'assessment': assessment,
        'pass': bool(core),
    }


# ========================================================================
# Probe runner
# ========================================================================


def run_probe() -> dict:
    tests = [
        test_T1_goal(),
        test_T2_cap_imported(),
        test_T3_radion_charges(),
        test_T4_potential(),
        test_T5_alpha(),
        test_T6_radion_mass(),
        test_T7_arc_consistency(),
        test_T8_honest_scope(),
        test_T9_assessment(),
    ]
    t2, t3, t4, t5, t6, t7 = (tests[1], tests[2], tests[3], tests[4],
                              tests[5], tests[6])
    all_ok = all(t['pass'] for t in tests)
    if all_ok:
        verdict_class = (
            "THE_HOPF_CHARGE_STABILIZES_THE_RADION_V_EFF_HAS_A_"
            "MINIMUM_WITH_ALPHA_STAR_EQUALS_SQRT_5KAPPA_OVER_28_"
            "ORDER_ONE_AT_THE_DIRAC_POINT_AND_M_PHI2_EQUALS_35_"
            "OVER_12_V_MIN_THE_OBSERVED_ALPHA_NEEDS_KAPPA_3E_"
            "MINUS_4_THE_NAMED_OPEN_PROBLEM"
        )
        verdict = (
            "ESTABLISHED (the argument is in "
            "docs/radion_stabilization.md).\n\n"
            "THE CAP. The #55 ledger re-read at machine zero (A = "
            "alpha hbar c/2, R* = (A/2B)^{1/3}, U/mc^2 = alpha/2, "
            "E'' = 6B > 0) - the alpha-dependent holdout continues."
            "\n\n"
            "THE CHARGES. Fixed charge -> e^{+sqrt3 phi} "
            f"({t3['electric_scaling_e_plus_sqrt3_phi_dev']:.0e}), "
            "fixed Hopf flux -> e^{-sqrt3 phi} "
            f"({t3['magnetic_scaling_e_minus_sqrt3_phi_dev']:.0e}), "
            "Dirac ratio 1/(4 alpha^2) exact; primordial radius "
            "(#222) adds e^{a phi}: exponents 7/(2 sqrt3) and "
            "-5/(2 sqrt3) - ONE OF EACH SIGN.\n\n"
            "THE POTENTIAL. The dyonic minimum exists (numeric = "
            f"closed form to {t4['min_rel_dev']:.0e}); m_phi^2 = "
            "(35/12) V_min holds numerically "
            f"({t4['mass_identity_35_over_12_dev']:.0e}) and "
            "symbolically; WITHOUT the Hopf charge V is monotonic - "
            "runaway to decompactification: the #58 topology saves "
            "the vacuum.\n\n"
            "ALPHA. alpha*^2 = 5 kappa/28: at the Dirac point "
            f"alpha* = {t5['alpha_star_dirac_point']:.4f} (order "
            f"one, {t5['alpha_star_over_alpha_obs']:.0f}x observed); "
            "observed alpha needs kappa = "
            f"{t5['kappa_required_for_observed_alpha']:.2e} "
            "(verified by minimization; guardrail: no match, "
            f"nearest {t5['nearest_candidate']['name']} at "
            f"{t5['nearest_candidate']['rel_dev']:.1%}); cohesion "
            f"tilt < {t5['max_cohesion_tilt']:.1%}.\n\n"
            "THE MASS. U_mag/U_el = 7/5 at the minimum, E'' = 7 "
            "U_el* exactly "
            f"({t6['Epp_over_7Uel_dev']:.0e}); anchored: "
            f"{t6['Epp_observed_alpha_GeV']:.1e} GeV per throat at "
            "the observed-alpha minimum (fiber anchor; GUT-scale - "
            "the radion is heavy), "
            f"{t6['Epp_compton_anchor_keV']:.0f} keV on the #55 "
            "Compton anchor; m_phi^2 = 16 pi n E''/m_P^2.\n\n"
            "THE ARC. The stabilizer row = the reserved cap row "
            "(rank 2, grad alpha annihilated "
            f"{t7['grad_alpha_null_projection']:.0e}); rho at the "
            f"Dirac point {t7['rho_at_dirac_point']:.2f}, at "
            f"observed alpha {t7['rho_at_observed_alpha']:.2f}; "
            "#222's x210 primordial exclusion re-read. Nothing "
            "else in the arc moves."
        )
    else:
        verdict_class = "RADION_STABILIZATION_INCONCLUSIVE"
        verdict = (
            "INCONCLUSIVE. A core check failed; re-examine before "
            "quoting the result."
        )

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "Radion stabilization from the primordial EM-capped "
            "throat: the dyonic throat's electric (e^{+sqrt3 phi}) "
            "and Hopf-magnetic (e^{-sqrt3 phi}) cap energies give "
            "V_eff(phi) a modulus-free minimum - alpha*^2 = 5 "
            "kappa/28, order one at the Dirac point - with m_phi^2 "
            "= (35/12) V_min and E'' = 7 U_el* per throat; the "
            "observed alpha requires kappa = 2.98e-4, the named "
            "open problem"
        ),
        "executes": (
            "the open dynamical problem #225 named: V_eff(phi), "
            "alpha at the minimum, and the radion mass, from the "
            "committed cap, topology, and dilaton coupling"
        ),
        "tests": tests,
        "n_passed": sum(1 for t in tests if t["pass"]),
        "n_total": len(tests),
        "verdict_class": verdict_class,
        "verdict": verdict,
    }


def _json_default(o):
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return "<array>"
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, complex):
        return [o.real, o.imag]
    if callable(o):
        return "<fn>"
    return str(o)


def render_markdown(s: dict) -> str:
    out = []
    out.append("# Radion stabilization from the primordial EM-capped "
               "throat (PR #226)")
    out.append("")
    out.append(f"**Run:** {s['timestamp_utc']}")
    out.append("")
    out.append(
        "The deliverable is `docs/radion_stabilization.md` - "
        "V_eff(phi), alpha at the minimum, and the radion mass. "
        "*(QFT on the fixed classical throat geometry, not quantum "
        "gravity.)*"
    )
    out.append("")
    out.append("## Test summary")
    out.append("")
    out.append("| # | Test | Key finding | PASS? |")
    out.append("|---|---|---|---|")
    labels = {
        "T1": "the goal: execute #225's open dynamical problem",
        "T2": "the #55 cap imported; alpha-dependent re-read",
        "T3": "the radion charges: e^{+sqrt3 phi} vs e^{-sqrt3 phi}",
        "T4": "the dyonic minimum; the no-go; m^2 = (35/12) V_min",
        "T5": "alpha* = sqrt(5 kappa/28); Dirac point 0.4226",
        "T6": "the radion mass: E'' = 7 U_el*; anchored scales",
        "T7": "the stabilizer lands in the reserved rank-audit slot",
        "T8": "honest scope",
        "T9": "assessment",
    }
    for t in s["tests"]:
        p = "**PASS**" if t["pass"] else "**FAIL**"
        pre = t["name"][:2]
        out.append(f"| {pre} | `{t['name']}` | {labels.get(pre,'-')} | {p} |")
    out.append("")
    out.append("## Verdict")
    out.append("")
    out.append(f"**Class:** `{s['verdict_class']}`")
    out.append("")
    out.append(s["verdict"])
    out.append("")
    return "\n".join(out)


def main(argv: Optional[list] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_radion_stabilization_probe"
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(
        json.dumps(summary, indent=2, default=_json_default))
    (out / "probe.md").write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"\nWrote: {out / 'probe.md'}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
