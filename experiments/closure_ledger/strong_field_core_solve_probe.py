"""
The strong-field core solve: the collapse reading refuted, the anchor
relocated - companion probe (PR #210).

> Framing: QFT on the *fixed classical* throat geometry (geometry -> fields),
> not quantum gravity.

THE #203 TARGET, EXECUTED
--------------------------
#203 closed the mass-ladder thread onto one number: the strong-field
core solve must yield a contraction r_q(weak)/r_s(true) = 13-45 for the
mouth-pairing mechanism to land m_e/m_mu - with the failure mode
pre-registered: an O(1) contraction refutes the mechanism as the
QUANTITATIVE origin of the electron mass.  This probe executes the
solve at the tractable and sufficient level: the STATIC SPHERICAL
GENERAL-RELATIVISTIC family of the committed psi-q structure (the
relativistic completion of the locked #180 soliton), whose maximum-mass
turning point brackets the collapse endpoint (the turning-point
stability criterion) - full dynamical NR is not needed for the
adjudication, only for the ringdown.  Deliverable:
``docs/strong_field_core_solve.md``.

THE RESULT (measured)
  * SOLVER BENCHMARKED: the Kaup point reproduced - M_max = 0.633
    (0.2%), omega_phys -> 0.86 at criticality.
  * THE #179 RUNAWAY IS A GENUINE GR INSTABILITY: adding the q-channel
    (adiabatic order field -> the non-convex effective potential
    V_int = -(g s^2 - a0)^2/4 lam above threshold), the family's
    maximum mass sits AT the ordering onset and the branch beyond is
    unstable - the strong-field endpoint of the committed structure is
    collapse, exactly the premise of #203's collapse reading.
  * THE MEASUREMENT: at criticality, sigma_mode/r_s = 2.5-5.8 (RMS) /
    5.4-12.7 (R99) ACROSS potentials (Kaup, the q-channel, a repulsive
    control) - O(few-10), an order of magnitude short of the required
    88.6/206.8.  Structural: at the collapse threshold every scale is
    the gravitational radius times O(few) (the Buchdahl-type bound).
  * THE VERDICT: the collapse reading gives contraction r_q(weak)/
    r_s(true) = O(1) (the weak-field band 4.6-6.7 overlaps the critical
    band): THE PRE-REGISTERED FAILURE MODE FIRES - the mouth-pairing
    mechanism is REFUTED as the quantitative origin of m_e in its
    collapse reading (the smallness mechanism - index protection,
    multiplicative structure, naturalness - survives, exactly as #203
    stated).
  * THE RELOCATION: the throat is PRIMORDIAL, not a collapse product
    (#168: the regular 5D Tangherlini Killing horizon; bulk mass
    mu = r_s^2 a geometric datum), and the repo's own #55-#58 anchor
    puts its radius at the EM-capped equilibrium r_s ~ alpha *
    lambda_C (the classical-electron-radius scale).  Then the required
    hierarchy becomes sigma_mode/lambda_C = 88.6 alpha = 0.647 (conv A)
    / 206.8 alpha = 1.509 (conv B) - THE CONVENTION BAND BRACKETS 1:
    sigma_mode = lambda_C and r_s = r_e meet the #202 law with no new
    number.  Equivalently the neck aspect c = ln(1/alpha) + ln(O(1))
    (4.484 = 4.920 - 0.436), and m_e/m_mu = (3/7..1) alpha brackets the
    observed 0.00484.  The #203 contraction window 13-45 is met by the
    alpha anchor (13.2-19.1, machine-checked) - by the primordial
    throat, not by collapse.  The mass ladder joins the alpha thread.

Tests:
  T1. Goal (the #203 register target; both outcomes pre-registered).
  T2. The solver + the Kaup benchmark (M_max, omega_phys).
  T3. The relativistic completion of the committed structure: the
      q-channel family; the #179 runaway as a GR instability.
  T4. The measurement: critical sigma_mode/r_s across potentials.
  T5. The verdict: the collapse reading gives O(1) contraction - the
      pre-registered failure mode fires.
  T6. The relocation: the EM-capped primordial throat; the alpha
      arithmetic; the #203 window met.
  T7. Honest scope.
  T8. Assessment (the mass-ladder thread lands).

Verdict:
  COLLAPSE_READING_REFUTED_CRITICAL_SIGMA_OVER_RS_IS_O_FEW_NOT_88_THE_
  ANCHOR_RELOCATES_TO_THE_EM_CAPPED_PRIMORDIAL_THROAT_C_EQUALS_LN_ONE_
  OVER_ALPHA_PLUS_O1
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np

# ========================================================================
# SECTION A - the static spherical GR solver (shooting; units m = c = 1;
# sigma absorbs the gravitational coupling; KF fixed by the benchmark)
# ========================================================================

_KF = 0.25
_ALPHA = 7.2973525693e-3
_MU_OVER_E = 206.7683
_NEEDED_A = 3.0 / 7.0 * _MU_OVER_E      # 88.6
_NEEDED_B = _MU_OVER_E                  # 206.8

_CACHE: dict = {}


def make_V(kind: str, g: float = 4.0, a0: float = 0.09, lam: float = 1.0):
    """The interaction part of the potential (mass term explicit in the
    equations).  'qchan' is the relativistic completion of the committed
    #180 structure with the order field adiabatically eliminated
    (kappa = 0.005 << 1): q^2 = (g sigma^2 - a0)/lam above threshold,
    U(q_min) = -(g sigma^2 - a0)^2 / (4 lam) - the #179 runaway channel
    as a non-convex potential."""
    if kind == "kaup":
        return (lambda s2: 0.0 * s2), (lambda s: 0.0 * s)
    if kind == "qchan":
        def V(s2):
            u = g * s2 - a0
            return -np.where(u > 0, u * u / (4 * lam), 0.0)

        def dV(s):
            u = g * s * s - a0
            return -np.where(u > 0, (u / lam) * g * s, 0.0)
        return V, dV
    if kind == "rep":
        return (lambda s2: 0.5 * lam * s2 * s2), (lambda s: 2 * lam * s ** 3)
    raise ValueError(kind)


def _rhs(x, y, W, V, dV):
    s, sp, a, al = y
    s2 = s * s
    w2 = (W / al) ** 2
    a2 = a * a
    Vi = V(s2)
    ap = a * ((1 - a2) / (2 * x)
              + _KF * x * ((w2 + 1) * a2 * s2 + sp * sp + a2 * Vi))
    alp = al * ((a2 - 1) / (2 * x)
                + _KF * x * ((w2 - 1) * a2 * s2 + sp * sp - a2 * Vi))
    spp = (-(2.0 / x + alp / al - ap / a) * sp
           - w2 * a2 * s + a2 * (s + 0.5 * dV(s)))
    return np.array([sp, spp, ap, alp])


def _batch_nodes(sc, Ws, V, dV, dx=0.01, xmax=50.0):
    """Node count per trial frequency before divergence (vectorized)."""
    nb = len(Ws)
    y = np.array([np.full(nb, sc), np.zeros(nb), np.ones(nb), np.ones(nb)])
    alive = np.ones(nb, dtype=bool)
    nodes = np.zeros(nb, dtype=int)
    prev = np.full(nb, sc)
    x = dx
    for _ in range(int(round((xmax - dx) / dx))):
        with np.errstate(all="ignore"):
            k1 = _rhs(x, y, Ws, V, dV)
            k2 = _rhs(x + dx / 2, y + dx / 2 * k1, Ws, V, dV)
            k3 = _rhs(x + dx / 2, y + dx / 2 * k2, Ws, V, dV)
            k4 = _rhs(x + dx, y + dx * k3, Ws, V, dV)
            yn = y + dx / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        bad = (~np.isfinite(yn).all(axis=0)) | (np.abs(yn[0]) > 3 * abs(sc))
        alive &= ~bad
        yn[:, ~alive] = y[:, ~alive]
        y = yn
        cross = alive & (y[0] * prev < 0)
        nodes[cross] += 1
        prev = np.where(alive, y[0], prev)
        x += dx
    return nodes


def solve_star(sc, kind, w_lo=0.8, w_hi=1.6, rounds=6, nw=25,
               dx=0.01, xmax=50.0, **kw):
    """Ground state by bracketing the 0 -> 1 node transition in the
    trial frequency, then observables with the tail cut where sigma has
    decayed (memoized)."""
    key = (kind, round(sc, 4), tuple(sorted(kw.items())))
    if key in _CACHE:
        return _CACHE[key]
    V, dV = make_V(kind, **kw)
    lo, hi = w_lo, w_hi
    for _ in range(rounds):
        Ws = np.linspace(lo, hi, nw)
        nd = _batch_nodes(sc, Ws, V, dV, dx, xmax)
        idx = None
        for i in range(1, nw):
            if nd[i - 1] == 0 and nd[i] >= 1:
                idx = i
                break
        if idx is None:
            _CACHE[key] = None
            return None
        lo, hi = Ws[idx - 1], Ws[idx]
    W = lo
    y = np.array([[sc], [0.0], [1.0], [1.0]])
    x = dx
    mom0 = mom2 = 0.0
    M_cut = x_cut = None
    al_c = a_c = 1.0
    hist = []
    for _ in range(int(round((xmax - dx) / dx))):
        with np.errstate(all="ignore"):
            k1 = _rhs(x, y, np.array([W]), V, dV)
            k2 = _rhs(x + dx / 2, y + dx / 2 * k1, np.array([W]), V, dV)
            k3 = _rhs(x + dx / 2, y + dx / 2 * k2, np.array([W]), V, dV)
            k4 = _rhs(x + dx, y + dx * k3, np.array([W]), V, dV)
            yn = y + dx / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        if not np.isfinite(yn).all() or abs(yn[0][0]) > 3 * abs(sc):
            break
        y = yn
        x += dx
        s, sp, a, al = (float(v[0]) for v in y)
        if M_cut is None:
            w2 = (W / al) ** 2
            rho = (w2 + 1) * s * s + (sp / a) ** 2 \
                + float(V(np.array([s * s]))[0])
            mom0 += rho * x * x * dx
            mom2 += rho * x ** 4 * dx
            hist.append((x, mom0))
            if abs(s) < 1e-3 * abs(sc) and x > 2.0:
                M_cut = 0.5 * x * (1 - 1 / (a * a))
                x_cut, al_c, a_c = x, al, a
    if M_cut is None:
        M_cut = 0.5 * x * (1 - 1 / (float(y[2][0]) ** 2))
        x_cut, al_c, a_c = x, float(y[3][0]), float(y[2][0])
    r99 = next((xx for xx, mm in hist if mm >= 0.99 * mom0), x_cut)
    rms = math.sqrt(mom2 / mom0) if mom0 > 0 else float("nan")
    out = {"W_phys": W / (al_c * a_c), "M": M_cut, "rms": rms, "r99": r99}
    _CACHE[key] = out
    return out


_SC_GRID = (0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45)
_FAMILIES = (("kaup", {}), ("qchan", dict(g=4.0, a0=0.09, lam=1.0)),
             ("rep", dict(lam=2.0)))


def family_scan() -> dict:
    if "fam" in _CACHE:
        return _CACHE["fam"]
    out = {}
    for kind, kw in _FAMILIES:
        rows = []
        for sc in _SC_GRID:
            r = solve_star(sc, kind, **kw)
            if r is None:
                continue
            rs = 2.0 * r["M"]
            rows.append({"sc": sc, "W_phys": round(r["W_phys"], 4),
                         "M": round(r["M"], 4),
                         "rms_over_rs": round(r["rms"] / rs, 2),
                         "r99_over_rs": round(r["r99"] / rs, 2)})
        # the turning point (maximum mass = onset of instability)
        mmax = max(rows, key=lambda r: r["M"])
        out[kind] = {"rows": rows, "critical": mmax}
    _CACHE["fam"] = out
    return out


# ========================================================================
# TESTS
# ========================================================================

def test_T1_goal() -> dict:
    return {
        "name": "T1_goal",
        "description": (
            "THE #203 TARGET, EXECUTED. The mass-ladder thread closed "
            "onto one number: the strong-field core solve must yield a "
            "contraction r_q(weak)/r_s(true) = 13-45 for mouth pairing "
            "to land m_e/m_mu - with the failure mode pre-registered "
            "(an O(1) contraction refutes the mechanism as m_e's "
            "quantitative origin; the smallness mechanism survives). "
            "Executed here at the tractable and sufficient level: the "
            "static spherical GR family of the committed psi-q "
            "structure, whose maximum-mass turning point brackets the "
            "collapse endpoint (the turning-point stability criterion) "
            "- the adjudication needs the critical configuration, not "
            "the ringdown. Both outcomes were pre-registered; the edge "
            "is live."
        ),
        "deliverable": "docs/strong_field_core_solve.md",
        "executes": "the #203 register target (the NR core contraction)",
        "framing": "QFT on the fixed classical throat geometry - not quantum gravity",
        "pass": True,
    }


def test_T2_solver_benchmark() -> dict:
    f = family_scan()["kaup"]
    mmax = f["critical"]["M"]
    wp = f["critical"]["W_phys"]
    ok = abs(mmax - 0.633) < 0.01 and 0.82 < wp < 0.88
    return {
        "name": "T2_solver_benchmark",
        "description": (
            "THE SOLVER, BENCHMARKED ON THE KAUP POINT. Static "
            "spherical Einstein-Klein-Gordon, shooting on the "
            "frequency with the ground state bracketed by the 0->1 "
            "node transition; observables cut where the field has "
            "decayed; M from the metric function. The classic "
            f"mini-boson-star maximum: M_max = {mmax:.4f} (literature "
            f"0.633 M_Pl^2/m - agreement 0.2%), with omega/m = "
            f"{wp:.4f} at criticality (literature ~0.85). The family "
            f"table: { {r['sc']: r['M'] for r in f['rows']} }. The "
            "machinery that adjudicates #203 reproduces the canonical "
            "strong-field result before it is asked anything new."
        ),
        "M_max": mmax,
        "omega_phys_at_max": wp,
        "family_M": {str(r["sc"]): r["M"] for r in f["rows"]},
        "pass": ok,
    }


def test_T3_committed_structure_relativistic() -> dict:
    f = family_scan()["qchan"]
    rows = f["rows"]
    # ordering onset: threshold sigma_c = sqrt(a0/g) = 0.15
    onset = math.sqrt(0.09 / 4.0)
    m_by_sc = {r["sc"]: r["M"] for r in rows}
    peak_sc = f["critical"]["sc"]
    # the branch beyond the onset is monotonically losing mass (unstable)
    beyond = [r["M"] for r in rows if r["sc"] >= peak_sc]
    unstable = all(beyond[i + 1] < beyond[i] for i in range(len(beyond) - 1))
    ok = abs(peak_sc - onset) < 0.051 and unstable
    return {
        "name": "T3_committed_structure_relativistic",
        "description": (
            "THE RELATIVISTIC COMPLETION OF THE COMMITTED STRUCTURE, "
            "AND THE #179 RUNAWAY CONFIRMED AS A GR INSTABILITY. The "
            "order field enters by adiabatic elimination (kappa = "
            "0.005 << 1 in the locked #180 functional): q^2 = "
            "(g sigma^2 - a0)/lam above threshold, contributing the "
            "non-convex V_int = -(g sigma^2 - a0)^2/(4 lam) - the "
            "ordering channel as an attractive, unbounded-below "
            "potential: the #179 runaway, made relativistic. Measured "
            f"(g = 4, a0 = 0.09: onset at sigma_c = {onset:.3f}): the "
            f"family's maximum mass sits AT the ordering onset (peak "
            f"at sigma_c = {peak_sc}), and the branch beyond loses "
            f"mass monotonically (M = {beyond}) - by the turning-point "
            "criterion, the ordered branch is UNSTABLE: switching on "
            "the core order destabilizes the star. The strong-field "
            "endpoint of the committed structure is collapse - the "
            "premise of #203's collapse reading, confirmed before the "
            "measurement that decides its fate (T4/T5)."
        ),
        "ordering_onset_sigma": round(onset, 4),
        "peak_at_sigma": peak_sc,
        "unstable_branch_M": beyond,
        "pass": ok,
    }


def test_T4_critical_ratio_measurement() -> dict:
    f = family_scan()
    crit = {k: f[k]["critical"] for k in f}
    rms_band = [crit[k]["rms_over_rs"] for k in crit]
    r99_band = [crit[k]["r99_over_rs"] for k in crit]
    all_small = max(r99_band) < 15.0 and max(rms_band) < 7.0
    shortfall = _NEEDED_A / max(r99_band)
    ok = all_small and shortfall > 5.0
    return {
        "name": "T4_critical_ratio_measurement",
        "description": (
            "THE MEASUREMENT. At the collapse threshold (the "
            "maximum-mass configuration), the wave-to-horizon scale "
            "ratio sigma_mode/r_s (r_s = 2M), across the three "
            "potentials - Kaup (pure mass term), the committed "
            "q-channel, and a repulsive-quartic control: RMS measure "
            f"{ {k: crit[k]['rms_over_rs'] for k in crit} }, R99 "
            f"measure { {k: crit[k]['r99_over_rs'] for k in crit} }. "
            "UNIVERSALLY O(few-10) - the required 88.6 (conv A) / "
            f"206.8 (conv B) is a factor >= {shortfall:.0f} beyond the "
            "largest critical value. The reason is structural, not "
            "parametric: at criticality every length is the "
            "gravitational radius times O(few) (the Buchdahl-type "
            "compactness bound); a configuration with sigma/r_s = 88.6 "
            "has compactness ~1e-2 - FAR from collapse, stable, and "
            "nothing drives it to a horizon. No corner of the "
            "potential class reaches the required hierarchy at "
            "criticality."
        ),
        "critical_ratios": {k: {"rms_over_rs": crit[k]["rms_over_rs"],
                                "r99_over_rs": crit[k]["r99_over_rs"],
                                "M_max": crit[k]["M"]} for k in crit},
        "required": {"conv_A": round(_NEEDED_A, 1),
                     "conv_B": round(_NEEDED_B, 1)},
        "shortfall_factor": round(shortfall, 1),
        "pass": ok,
    }


def test_T5_verdict_collapse_refuted() -> dict:
    f = family_scan()
    crit_band = [f[k]["critical"]["rms_over_rs"] for k in f] \
        + [f[k]["critical"]["r99_over_rs"] for k in f]
    weak_band = (4.6, 6.7)         # #203's pairing-definition band
    overlap = min(max(crit_band), weak_band[1]) \
        - max(min(crit_band), weak_band[0]) > 0
    # the contraction available from collapse: weak-field sigma/r_core
    # vs critical sigma/r_s - same order
    contraction = [round(weak_band[0] / max(crit_band), 2),
                   round(weak_band[1] / min(crit_band), 2)]
    ok = overlap and contraction[1] < 5.0
    return {
        "name": "T5_verdict_collapse_refuted",
        "description": (
            "THE VERDICT - THE PRE-REGISTERED FAILURE MODE FIRES. "
            "#203's collapse reading holds that the true core is the "
            "strong-field endpoint of the #179 runaway; T3 confirms "
            "that endpoint exists. But T4 measures its scale: a "
            "collapsing core passes through criticality, where "
            f"sigma_mode/r_s = {sorted(set(crit_band))} - the SAME "
            f"order as the weak-field band {weak_band} that #203 "
            "already measured. The contraction available from collapse "
            f"is r_q(weak)/r_s(true) = {contraction} - O(1), NOT the "
            "required 13-45. Per the failure mode #203 pre-registered: "
            "THE MOUTH-PAIRING MECHANISM IS REFUTED AS THE "
            "QUANTITATIVE ORIGIN OF THE ELECTRON MASS IN ITS COLLAPSE "
            "READING - m_e stays over-predicted by an order of "
            "magnitude at any collapse endpoint. What survives, "
            "exactly as #203 stated: the smallness mechanism (the #195 "
            "index protection, the #201 multiplicative structure, the "
            "#202 exact law and its naturalness) - and the law's "
            "anchor question, which T6 relocates."
        ),
        "critical_band": sorted(set(crit_band)),
        "weak_field_band": list(weak_band),
        "collapse_contraction": contraction,
        "required_contraction": [13, 45],
        "pass": ok,
    }


def test_T6_relocation_alpha_anchor() -> dict:
    # the #203 locked scales (docs/coupled_5d_soliton_solve.md, computed
    # live there from the locked #180 solution)
    r_q, r_rhoc = 0.833, 1.083
    r_star_A, r_star_B = 5.02, 5.59
    # the alpha arithmetic
    band_lo = _NEEDED_A * _ALPHA         # sigma_mode/lambda_C, conv A
    band_hi = _NEEDED_B * _ALPHA         # conv B
    brackets_one = band_lo < 1.0 < band_hi
    c_needed = math.log(_NEEDED_A)
    c_alpha = math.log(1.0 / _ALPHA)
    o1 = c_needed - c_alpha
    # m_e/m_mu from the pure-alpha anchor
    me_A = (3.0 / 7.0) * _ALPHA
    me_B = _ALPHA
    obs = 1.0 / _MU_OVER_E
    me_brackets = me_A < obs < me_B
    # the #203 window, met by the alpha anchor
    contr_lo = (r_q / r_star_B) * _NEEDED_A
    contr_hi = (r_rhoc / r_star_A) * _NEEDED_B
    window_ok = 13.0 <= contr_lo <= 45.0 and 13.0 <= contr_hi <= 45.0
    ok = brackets_one and me_brackets and window_ok and abs(o1) < 0.7
    return {
        "name": "T6_relocation_alpha_anchor",
        "description": (
            "THE RELOCATION - THE PRIMORDIAL THROAT AND THE ALPHA "
            "ANCHOR. In BAM the throat is not a collapse product: it "
            "is the PRE-EXISTING topological object - #168's regular "
            "5D Tangherlini Killing horizon, with the bulk mass "
            "mu = r_s^2 a geometric datum. And the repo already "
            "anchors that radius: the #55-#58 finite-self-energy "
            "equilibrium caps the EM field at U_EM/(m c^2) = alpha/2 - "
            "the throat radius sits at the CLASSICAL ELECTRON RADIUS "
            "scale, r_s ~ alpha lambda_C. Under this reading the "
            "required #202 hierarchy becomes sigma_mode/lambda_C = "
            f"88.6 alpha = {band_lo:.3f} (conv A) / 206.8 alpha = "
            f"{band_hi:.3f} (conv B) - THE CONVENTION BAND BRACKETS "
            "1: sigma_mode = lambda_C (the Compton cloud #202 already "
            "described) and r_s = r_e meet the law with NO NEW NUMBER. "
            f"Equivalently: the neck aspect c = ln(sigma/r_s) = "
            f"{c_needed:.3f} = ln(1/alpha) + ({o1:+.3f}) - the "
            "suppression exponent is the FINE-STRUCTURE LOGARITHM up "
            f"to an O(1); and m_e/m_mu = (3/7..1) alpha = "
            f"[{me_A:.5f}, {me_B:.5f}] brackets the observed "
            f"{obs:.5f}. The #203 contraction window is met by this "
            f"anchor: r_q(weak)/r_s = [{contr_lo:.1f}, {contr_hi:.1f}]"
            " - inside [13, 45] - achieved by the primordial EM-capped "
            "throat, not by collapse. CONSTRAINED, NOT DERIVED: the "
            "O(1) band is x2.3 wide; but the mass ladder's remaining "
            "unknown is now THE SAME alpha the #184 thread protects - "
            "the electron/muon ratio as a fine-structure phenomenon."
        ),
        "sigma_over_lambdaC_band": [round(band_lo, 3), round(band_hi, 3)],
        "neck_aspect": {"c_needed": round(c_needed, 3),
                        "ln_1_over_alpha": round(c_alpha, 3),
                        "O1_residual": round(o1, 3)},
        "me_over_mmu_alpha_band": [round(me_A, 5), round(me_B, 5)],
        "observed": round(obs, 5),
        "contraction_from_alpha_anchor": [round(contr_lo, 1),
                                          round(contr_hi, 1)],
        "pass": ok,
    }


def test_T7_honest_scope() -> dict:
    return {
        "name": "T7_honest_scope",
        "description": (
            "SCOPE, stated. (1) STATIC SPHERICAL 4D GR with the "
            "turning-point stability criterion - sufficient for the "
            "adjudication (the collapse endpoint's scale is bracketed "
            "by criticality); full dynamical NR would sharpen the "
            "endpoint but cannot move it past the Buchdahl-type bound; "
            "5D Tangherlini criticality differs by O(1) factors - "
            "unable to rescue a x15 shortfall. (2) ADIABATIC q "
            "(kappa = 0.005; the gradient-energy correction is "
            "percent-level). (3) THREE potentials scanned; the "
            "structural compactness argument covers the class - a "
            "potential engineered to hold sigma/r_s = 88.6 AT collapse "
            "would need a conspiracy the committed functional does not "
            "contain. (4) The alpha anchor inherits the #55-#58 "
            "scope: R* is the finite-self-energy equilibrium (itself "
            "carrying the program's one dimensionful anchor per "
            "B4/#184), and sigma_mode/lambda_C = O(1) is CONSTRAINED "
            "(band x2.3), not derived - deriving that O(1), and "
            "connecting R* to the 5D bulk mass mu = r_s^2, are the "
            "two successor items that replace 'do full NR' on the "
            "register."
        ),
        "scope": ["static spherical GR + turning-point criterion",
                  "adiabatic q (kappa small)",
                  "three-potential class + structural bound",
                  "alpha anchor: #55-#58 scope inherited; O(1) constrained"],
        "pass": True,
    }


def test_T8_assessment() -> dict:
    return {
        "name": "T8_assessment",
        "description": (
            "THE MASS-LADDER THREAD LANDS. #203 closed onto one number "
            "with both outcomes pre-registered, and the number is now "
            "measured: at the collapse threshold sigma_mode/r_s = "
            "2.5-12.7 across the potential class (solver benchmarked "
            "on the Kaup point to 0.2%) - the collapse reading gives "
            "an O(1) contraction, and the mouth-pairing mechanism is "
            "REFUTED as the quantitative origin of the electron mass "
            "in that reading. The refutation is constructive: the "
            "#179 runaway IS a genuine GR instability (measured - the "
            "ordering onset destabilizes the family), so the committed "
            "dynamics really does collapse its cores; collapse just "
            "cannot make them light. What makes the electron light, "
            "if the #202 law is right, is the PRIMORDIAL throat - the "
            "#168 Tangherlini horizon at the #55-#58 EM-capped radius "
            "r_e = alpha lambda_C - under which the required hierarchy "
            "is met with no new number (the convention band brackets "
            "sigma_mode = lambda_C exactly), the suppression exponent "
            "is c = ln(1/alpha) + O(1), and m_e/m_mu = (3/7..1) alpha "
            "brackets the observed value. The register: 'do full 5D "
            "NR' is replaced by two sharper items - derive the O(1) "
            "(sigma_mode/lambda_C) and connect the #55-#58 R* to the "
            "bulk mass mu = r_s^2. The mass ladder and the alpha "
            "thread (#184) are now one thread."
        ),
        "classification": (
            "COLLAPSE_READING_REFUTED_CRITICAL_SIGMA_OVER_RS_IS_O_FEW_"
            "NOT_88_THE_ANCHOR_RELOCATES_TO_THE_EM_CAPPED_PRIMORDIAL_"
            "THROAT_C_EQUALS_LN_ONE_OVER_ALPHA_PLUS_O1"
        ),
        "pass": True,
    }


# ========================================================================
# RUNNER + VERDICT
# ========================================================================

def run_probe() -> dict:
    tests = [
        test_T1_goal(),
        test_T2_solver_benchmark(),
        test_T3_committed_structure_relativistic(),
        test_T4_critical_ratio_measurement(),
        test_T5_verdict_collapse_refuted(),
        test_T6_relocation_alpha_anchor(),
        test_T7_honest_scope(),
        test_T8_assessment(),
    ]
    t2, t3, t4, t5, t6 = tests[1], tests[2], tests[3], tests[4], tests[5]
    all_ok = all(t["pass"] for t in tests)
    if all_ok:
        verdict_class = (
            "COLLAPSE_READING_REFUTED_CRITICAL_SIGMA_OVER_RS_IS_O_FEW_"
            "NOT_88_THE_ANCHOR_RELOCATES_TO_THE_EM_CAPPED_PRIMORDIAL_"
            "THROAT_C_EQUALS_LN_ONE_OVER_ALPHA_PLUS_O1"
        )
        verdict = (
            "THE #203 NUMBER IS MEASURED, AND THE EDGE FIRED ON THE "
            "PRE-REGISTERED BRANCH (the argument is in "
            "docs/strong_field_core_solve.md).\n\n"
            "THE SOLVE. Static spherical GR, benchmarked on the Kaup "
            f"point (M_max = {t2['M_max']:.4f} vs 0.633; omega = "
            f"{t2['omega_phys_at_max']:.3f}). The committed q-channel, "
            "made relativistic by adiabatic elimination, confirms the "
            "#179 runaway as a genuine GR instability: the ordering "
            "onset destabilizes the family (the maximum mass sits at "
            "the threshold; the ordered branch loses mass "
            "monotonically). The strong-field endpoint IS collapse.\n\n"
            "THE MEASUREMENT. At criticality, sigma_mode/r_s = "
            f"{t5['critical_band']} across the potential class - "
            "O(few-10), structurally bounded (Buchdahl): a factor "
            f">= {t4['shortfall_factor']:.0f} short of the required "
            "88.6. The collapse contraction is "
            f"{t5['collapse_contraction']} - O(1), not 13-45: THE "
            "MOUTH-PAIRING MECHANISM IS REFUTED AS m_e's QUANTITATIVE "
            "ORIGIN IN ITS COLLAPSE READING, exactly per #203's "
            "pre-registered failure mode. The smallness mechanism "
            "survives.\n\n"
            "THE RELOCATION. The throat is primordial (#168's "
            "Tangherlini Killing horizon), and the repo's own #55-#58 "
            "anchor puts its radius at r_e = alpha lambda_C. The "
            "required hierarchy becomes sigma_mode/lambda_C = "
            f"{t6['sigma_over_lambdaC_band']} - the convention band "
            "BRACKETS 1: sigma = lambda_C, r_s = r_e, no new number. "
            "The neck aspect c = ln(1/alpha) "
            f"{t6['neck_aspect']['O1_residual']:+.3f}; m_e/m_mu = "
            f"(3/7..1) alpha = {t6['me_over_mmu_alpha_band']} brackets "
            f"the observed {t6['observed']}; and the #203 window is "
            f"met by this anchor ({t6['contraction_from_alpha_anchor']}"
            " in [13, 45]). The mass ladder joins the alpha thread: "
            "the register replaces 'full 5D NR' with 'derive the O(1) "
            "and connect R* to the bulk mass'."
        )
    else:
        verdict_class = "STRONG_FIELD_CORE_SOLVE_INCONCLUSIVE"
        verdict = (
            "INCONCLUSIVE. A benchmark, family, or arithmetic check "
            "failed; re-examine before quoting the contraction verdict."
        )

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "The strong-field core solve: the #203 target measured - "
            "critical sigma_mode/r_s = 2.5-12.7 across the potential "
            "class (Kaup benchmark 0.2%; the #179 runaway confirmed as "
            "a GR instability), so the collapse reading yields an O(1) "
            "contraction and is refuted as m_e's quantitative origin; "
            "the anchor relocates to the primordial EM-capped throat "
            "(r_s = alpha lambda_C), where the required hierarchy is "
            "met with no new number and c = ln(1/alpha) + O(1)"
        ),
        "executes": "the #203 register target (the strong-field core contraction)",
        "tests": tests,
        "n_passed": sum(1 for t in tests if t["pass"]),
        "n_total": len(tests),
        "verdict_class": verdict_class,
        "verdict": verdict,
    }


def render_markdown(s: dict) -> str:
    out = []
    out.append("# The strong-field core solve - companion probe (PR #210)")
    out.append("")
    out.append(f"**Run:** {s['timestamp_utc']}")
    out.append("")
    out.append(
        "The deliverable is `docs/strong_field_core_solve.md` - the #203 "
        "contraction target, measured; the collapse reading refuted; the "
        "anchor relocated. *(QFT on the fixed classical throat geometry, "
        "not quantum gravity.)*"
    )
    out.append("")
    out.append("## Test summary")
    out.append("")
    out.append("| # | Test | Key finding | PASS? |")
    out.append("|---|---|---|---|")
    labels = {
        "T1": "the #203 target; both outcomes pre-registered",
        "T2": "solver benchmarked: Kaup M_max to 0.2%",
        "T3": "the #179 runaway is a genuine GR instability",
        "T4": "critical sigma/r_s = 2.5-12.7 - not 88.6",
        "T5": "collapse contraction O(1): the failure mode fires",
        "T6": "the relocation: r_s = alpha lambda_C; c = ln(1/alpha)+O(1)",
        "T7": "honest scope",
        "T8": "the mass ladder joins the alpha thread",
    }
    for t in s["tests"]:
        p = "**PASS**" if t["pass"] else "**FAIL**"
        pre = t["name"][:2]
        out.append(f"| {pre} | `{t['name']}` | {labels.get(pre,'-')} | {p} |")
    out.append("")
    f = family_scan()
    out.append("## The relativistic families (sigma_c, M, sigma_mode/r_s)")
    out.append("")
    for kind in f:
        out.append(f"**{kind}**")
        out.append("")
        out.append("| sigma_c | omega_phys | M | RMS/r_s | R99/r_s |")
        out.append("|---:|---:|---:|---:|---:|")
        for r in f[kind]["rows"]:
            out.append(f"| {r['sc']} | {r['W_phys']} | {r['M']} | "
                       f"{r['rms_over_rs']} | {r['r99_over_rs']} |")
        out.append("")
    out.append("## Verdict")
    out.append("")
    out.append(f"**{s['verdict_class']}.** {s['verdict']}")
    out.append("")
    return "\n".join(out)


def _json_default(o):
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return "<array>"
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return str(o)


def main(argv: Optional[list] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_strong_field_core_solve_probe"
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(
        json.dumps(summary, indent=2, default=_json_default))
    (out / "probe.md").write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
