"""Pre-flight checks for the Z2 link-twist plan (PR #227 research plan).

NOT a probe.  These are minimal reduced models that establish the five
claims of ``docs/nonorientable_er_link_z2_gauge_research_plan.md`` are
non-vacuous and fix the expected signs and magnitudes BEFORE the real
probe is written.  Every number printed here is quoted in that document.

The probe proper must redo all of it on the repo's real machinery (the
#224 two-throat network, the #166 conformal zonal solver, the #137
cubic-vertex ledger); a pass here is a licence to build that, not a
result about the BAM geometry.

Runs in a few seconds, no repo imports, numpy only.

    python -m experiments.closure_ledger.nonorientable_er_link_z2_preflight
"""

from __future__ import annotations

import numpy as np


# ── C1: a twisted link freezes the mouth-exchange doublet ────────────────────

def check_frustrated_doublet() -> list:
    """Two mouths joined by TWO exterior arcs (the #223/#224 ring).

    The doublet coupling is the sum over paths, J = J1 + s*J2 with
    s = W(gamma) the loop holonomy, so a twisted arc makes the two paths
    subtract.  On the symmetric network that is exact cancellation:
    Delta_omega = 0, exchange period pi/Delta_omega -> infinity.
    """
    rows = []
    for j1, j2 in [(1.0, 1.0), (1.0, 0.7)]:
        for s in (+1, -1):
            j = j1 + s * j2
            ev = np.linalg.eigvalsh(np.array([[0.0, -j], [-j, 0.0]]))
            rows.append({'J1': j1, 'J2': j2, 'holonomy': s,
                         'splitting': float(ev[1] - ev[0])})
    return rows


# ── C2/C3: Z2 gauge invariance, and the twisted loop is antiperiodic ─────────

def _ring(signs: np.ndarray, bonds: np.ndarray) -> np.ndarray:
    n = len(bonds)
    h = np.zeros((n, n))
    for i in range(n):
        j = (i + 1) % n
        h[i, j] = h[j, i] = -signs[i] * bonds[i]
    return np.linalg.eigvalsh(h)


def check_gauge_invariance(n: int = 6, trials: int = 4, seed: int = 0) -> dict:
    """Only the Wilson loop is physical.

    A mouth reframing eps_b -> s_a eps_b s_a' is a similarity transform of
    the ring operator, so it cannot move the spectrum at all; a genuine
    holonomy flip does.
    """
    rng = np.random.default_rng(seed)
    bonds = rng.uniform(0.5, 1.5, n)
    base = _ring(np.ones(n), bonds)

    drifts = []
    for _ in range(trials):
        g = rng.choice([-1.0, 1.0], n)
        signs = np.array([g[i] * g[(i + 1) % n] for i in range(n)])
        drifts.append(float(np.max(np.abs(_ring(signs, bonds) - base))))

    tw = np.ones(n)
    tw[0] = -1.0
    return {
        'gauge_drifts': drifts,
        'max_gauge_drift': max(drifts),
        'twisted_shift': float(np.max(np.abs(_ring(tw, bonds) - base))),
    }


def check_antiperiodic_shift(m: int = 240) -> dict:
    """The twisted loop carries half-integer modes.

    The same integer -> half-integer shift that lifts the Moebius glueball
    tower by pi*sigma in M^2 (#100) and gives the Moebius baryons their
    2*sqrt(sigma) gap (#103).
    """
    out = {}
    for name, signs in (('untwisted', np.ones(m)),
                        ('twisted', np.concatenate([[-1.0], np.ones(m - 1)]))):
        e = _ring(signs, np.ones(m))
        k = np.arccos(np.clip(-e / 2.0, -1.0, 1.0)) * m / (2 * np.pi)
        out[name] = [round(float(v), 3) for v in np.sort(k)[:6]]
    return out


# ── C4: the quotient revival — period pi*R, sign set by the twist ────────────

def check_quotient_revival(lmax: int = 160, npts: int = 4001,
                           l0: int = 40, width: float = 12.0) -> list:
    """Conformal zonal wave on S^3, restricted to one twist sector.

    Reduced field f = sin(chi)*psi has modes sin((l+1)chi) with
    omega_l = (l+1)/R (R = 1).  Sections of the trivial flat bundle over
    the quotient are the even-l modes, of the twisted bundle the odd-l
    modes.  P is the mirror chi -> pi - chi.

    Expected: the full tower reproduces #166 (f(pi R) = -P f0); each twist
    sector instead revives EXACTLY at pi R, with sign -1 (untwisted) or
    +1 (twisted).
    """
    chi = np.linspace(1e-6, np.pi - 1e-6, npts)
    ls = np.arange(lmax + 1)
    coef = np.exp(-((ls - l0) ** 2) / (2.0 * width ** 2))
    basis = np.sin(np.outer(ls + 1, chi))

    def field(mask, t):
        return (coef[mask, None] * basis[mask]
                * np.cos((ls[mask] + 1)[:, None] * t)).sum(axis=0)

    masks = {'full': np.ones_like(ls, dtype=bool),
             'untwisted_even_l': ls % 2 == 0,
             'twisted_odd_l': ls % 2 == 1}

    rows = []
    for name, mask in masks.items():
        f0 = field(mask, 0.0)
        n0 = np.linalg.norm(f0)
        for t, lbl in ((np.pi, 'pi R'), (2 * np.pi, '2 pi R')):
            ft = field(mask, t)
            rows.append({
                'sector': name, 't': lbl,
                'rel_f_minus_f0': float(np.linalg.norm(ft - f0) / n0),
                'rel_f_plus_f0': float(np.linalg.norm(ft + f0) / n0),
                'rel_f_plus_Pf0': float(np.linalg.norm(ft + f0[::-1]) / n0),
            })
    return rows


# ── C5: the vertex selection rule generalises ───────────────────────────────

def check_vertex_rule(lmax: int = 3) -> list:
    """Deck invariance of the integrand needs (-1)^{Sum l} = prod eps_i,
    i.e. Sum l == n_twist (mod 2).  #137's 'Sum l even' is n_twist = 0.
    """
    rows = []
    for n_tw in range(4):
        keep = [(a, b, c)
                for a in range(lmax + 1) for b in range(lmax + 1)
                for c in range(lmax + 1)
                if (-1) ** (a + b + c) == (-1) ** n_tw]
        rows.append({'n_twist': n_tw,
                     'surviving': len(keep),
                     'total': (lmax + 1) ** 3,
                     'sum_l_parity': sorted({sum(t) % 2 for t in keep})})
    return rows


def main() -> int:
    print("Z2 link-twist pre-flight (research plan PR #227)\n")

    print("C1  a twisted link freezes the mouth-exchange doublet")
    for r in check_frustrated_doublet():
        print(f"     J1={r['J1']} J2={r['J2']}  W={r['holonomy']:+d}"
              f"   splitting = {r['splitting']:.6f}")

    g = check_gauge_invariance()
    print("\nC2  Z2 gauge invariance (only the Wilson loop is physical)")
    print(f"     max spectral drift under gauge transforms = "
          f"{g['max_gauge_drift']:.2e}")
    print(f"     shift under a genuine W = -1 twist        = "
          f"{g['twisted_shift']:.4f}")

    print("\nC3  the twisted loop is antiperiodic (mode index k*M/2pi)")
    for name, ks in check_antiperiodic_shift().items():
        print(f"     {name:10s} {ks}")

    print("\nC4  quotient revival: period pi R, sign set by the twist class")
    for r in check_quotient_revival():
        print(f"     {r['sector']:18s} t={r['t']:7s}"
              f" |f-f0|={r['rel_f_minus_f0']:.2e}"
              f" |f+f0|={r['rel_f_plus_f0']:.2e}"
              f" |f+Pf0|={r['rel_f_plus_Pf0']:.2e}")

    print("\nC5  vertex selection rule: Sum l == n_twist (mod 2)")
    for r in check_vertex_rule():
        print(f"     n_twist={r['n_twist']}  surviving {r['surviving']}"
              f"/{r['total']}  Sum-l parity {r['sum_l_parity']}")

    print("\nAll five checks are reduced models, not BAM-geometry results.")
    print("See docs/nonorientable_er_link_z2_gauge_research_plan.md sec. 4.")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
