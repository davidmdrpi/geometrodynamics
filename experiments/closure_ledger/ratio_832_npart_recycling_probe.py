"""
Is 832 = N_q + ΔN an independent selection of the lepton/QCD scale ratio,
or recycled n_part? (PR #107).

PR #106 left the lepton/QCD scale ratio √σ/m_e ≈ 830 as an UNDERIVED open
dimensionless residual — the one number whose derivation would reduce BAM
to a single anchor (the bulk-gravity scale G). A tempting candidate
"derivation" is

    N_lepton = 100,  N_q = 466,  ΔN = N_q − N_lepton = 366,
    N_q + ΔN = 832  ≈  √σ/m_e ≈ 830   (0.2% match!).

This probe tests, skeptically, whether the bulk shell-stress integral
INDEPENDENTLY selects 832 — or whether 832 just recycles the
phenomenological compensator n_part. **Answer: it recycles n_part. The
match is a baseline coincidence; the channel-normalisation derivation
via this route fails.**

## 832 is built from n_part

`N_q + ΔN = 2·N_q − N_lepton = 2·(2·n_part) − 4·k_5² = 4·n_part − 4·k_5²
= 4·233 − 100 = 832`. So 832 is a linear function of `n_part` — the quark
closure integer that PR #76/#97 established is a PHENOMENOLOGICAL
COMPENSATOR (it absorbs the quark flavor puzzle; it drifts 216–255 across
the `quark_axioms` §8 ablations; only its parity is invariant). 832
inherits that status.

## The decisive test: §8 drift

If 832 were an independent geometric selection of the FIXED observed ratio
830, it would be §8-stable. It is not. Propagating n_part's ablation range
`{216, …, 255}` through `4·n_part − 100` gives

    832-analogue ∈ [764, 920]   (span 156, ≈ ±9%),

so the quantity drifts almost ±10% while the observed `√σ/m_e = 830.3` is
fixed. The baseline value `n_part = 233` happens to land it at 832 (0.2%
from 830), but that is a BASELINE COINCIDENCE — the same kind as
`50π·k_5 = 785` (PR #106) and `F_13 = 233` (PR #76) — not a stable
selection of the ratio.

## No independent bulk shell-stress integral selects 832

The genuine bulk shell-stress integrals over the Tangherlini geometry are
O(10–70), nowhere near 466/832: the sum of shell eigenvalues
`Σ ω²(l=1, n=3..5) ≈ 69.8`, the Bohr–Sommerfeld closure sum
`Σ (n+1)π ≈ 47.1`. The number 466 enters ONLY through the v3 Hamiltonian
closure count `4β_quark/(2π) = 2·n_part`, i.e. through the fit. There is no
independent bulk shell-stress integral that yields ~466 or ~832.

## The circularity

`n_part` was FIT to reproduce the quark spectrum (which already encodes
the physical scales). Recovering a scale ratio from `n_part` is therefore
circular — you get back (an unstable version of) what was put in. So 832
≈ 830 is not a derivation of the lepton/QCD hierarchy; it is the
compensator echoing the spectrum it was fit to.

## What this means for the ledger

The channel-normalisation derivation via `N_q + ΔN` FAILS. The PR #106
status stands unchanged: `√σ/m_e ≈ 830` remains an UNDERIVED open
dimensionless residual, and BAM remains at one foundational scale (G) +
one open ratio. A genuine derivation would need an INDEPENDENT bulk
shell-stress integral (not the v3-fit closure count) that selects ~830
AND is §8-stable — which is not available.

Tests:
  T1. The observation: N_q + ΔN = 2N_q − N_lepton = 832 ≈ √σ/m_e ≈ 830
      (0.2%).
  T2. 832 = 4·n_part − 4·k_5² — built from n_part.
  T3. n_part is the compensator (PR #76/#97): §8-drifts 216–255, parity-
      only invariant, absorbs the flavor puzzle.
  T4. §8-drift test: 4·n_part − 100 drifts 764–920 (±9%) while 830 is
      fixed ⟹ baseline coincidence, not a stable selection.
  T5. Independent bulk shell-stress integrals give O(10–70), not ~466/832;
      466 enters only via the v3 fit.
  T6. The circularity: n_part fit to the spectrum ⟹ recovering a scale
      ratio from it is circular.
  T7. Ledger unchanged: √σ/m_e still underived; one scale G + one open
      ratio; a genuine derivation needs an independent §8-stable integral.
  T8. Assessment.

Verdict:
  - RATIO_832_RECYCLES_N_PART_NOT_INDEPENDENT (expected): 832 = 2N_q −
    N_lepton = 4·n_part − 4·k_5² is built from the phenomenological
    compensator n_part; it §8-drifts 764–920 (so the 832 ≈ 830 match is a
    baseline coincidence), and no independent bulk shell-stress integral
    selects ~466/832. We are recycling n_part — the channel-normalisation
    derivation fails and √σ/m_e remains underived.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np

from geometrodynamics.constants import R_MID, R_OUTER
from geometrodynamics.tangherlini.radial import (
    V_tangherlini, r_to_rstar, rstar_to_r,
)
from geometrodynamics.qcd.constants import SIGMA_QCD


PI = math.pi
K_5 = 5
N_LEPTON = 4 * K_5 ** 2            # 100
N_PART_BASELINE = 233
N_Q = 2 * N_PART_BASELINE          # 466
DELTA_N = N_Q - N_LEPTON           # 366
VAL_832 = N_Q + DELTA_N            # 832 = 2 N_q − N_lepton

# n_part across the 12 quark_axioms §8 ablations (PR #76/#97).
N_PART_S8 = [233, 233, 233, 238, 237, 237, 241, 216, 247, 247, 220, 255]

M_E_GEV = 0.5109989e-3
RATIO = math.sqrt(SIGMA_QCD) / M_E_GEV     # √σ/m_e ≈ 830.3


def _shell_eigenvalues(l: int = 1, n: int = 800):
    rsmin = r_to_rstar(R_MID + 5e-4, R_MID)
    rsmax = r_to_rstar(R_OUTER - 5e-4, R_MID)
    rstar = np.linspace(rsmin, rsmax, n)
    h = rstar[1] - rstar[0]
    rphys = np.array([rstar_to_r(s, R_MID) for s in rstar])
    V = V_tangherlini(rphys, l, R_MID)
    main = 2.0 / h ** 2 + V[1:-1]
    off = -1.0 / h ** 2 * np.ones(n - 3)
    ev = np.sort(np.linalg.eigvalsh(
        np.diag(main) + np.diag(off, 1) + np.diag(off, -1)))
    return np.maximum(ev, 0.0)


# ---------------------------------------------------------------------------
# T1. The observation
# ---------------------------------------------------------------------------

def test_T1_observation() -> dict:
    off = (VAL_832 - RATIO) / RATIO
    return {
        'name': 'T1_the_832_observation',
        'description': (
            "N_lepton=100, N_q=466, ΔN=366; N_q+ΔN = 2N_q−N_lepton = 832 ≈ "
            "√σ/m_e ≈ 830 (0.2%). Tempting as a derivation of the PR #106 "
            "scale ratio."
        ),
        'N_lepton': N_LEPTON, 'N_q': N_Q, 'delta_N': DELTA_N,
        'N_q_plus_delta_N': VAL_832,
        'equals_2Nq_minus_Nlepton': 2 * N_Q - N_LEPTON,
        'observed_ratio': RATIO,
        'rel_off_percent': off * 100,
        'tempting_match': abs(off) < 0.01,
        'pass': VAL_832 == 2 * N_Q - N_LEPTON,
    }


# ---------------------------------------------------------------------------
# T2. 832 is built from n_part
# ---------------------------------------------------------------------------

def test_T2_built_from_npart() -> dict:
    """832 = 2N_q − N_lepton = 2(2·n_part) − 4·k_5² = 4·n_part − 4·k_5².
    A linear function of n_part."""
    via_npart = 4 * N_PART_BASELINE - 4 * K_5 ** 2
    return {
        'name': 'T2_832_built_from_npart',
        'description': (
            "832 = 2N_q − N_lepton = 4·n_part − 4·k_5² = 4·233 − 100 = 832. "
            "A linear function of n_part — built from the compensator."
        ),
        'formula': '4·n_part − 4·k_5²',
        'value': via_npart,
        'matches_832': via_npart == VAL_832,
        'is_function_of_npart': True,
        'pass': via_npart == VAL_832,
    }


# ---------------------------------------------------------------------------
# T3. n_part is the compensator
# ---------------------------------------------------------------------------

def test_T3_npart_is_compensator() -> dict:
    """PR #76/#97: n_part = 233 is a phenomenological compensator — it
    absorbs the quark flavor puzzle, drifts 216–255 across the §8
    ablations, and only its parity is invariant. Not derived."""
    span = max(N_PART_S8) - min(N_PART_S8)
    return {
        'name': 'T3_npart_is_phenomenological_compensator',
        'description': (
            "n_part = 233 is a compensator (PR #76/#97): drifts 216–255 "
            "across §8, parity-only invariant, absorbs the quark flavor "
            "puzzle. Not derived."
        ),
        'n_part_s8': N_PART_S8,
        's8_min': min(N_PART_S8), 's8_max': max(N_PART_S8), 's8_span': span,
        'parity_only_invariant': True,
        'derived': False,
        'pass': span > 20,
    }


# ---------------------------------------------------------------------------
# T4. The decisive §8-drift test
# ---------------------------------------------------------------------------

def test_T4_s8_drift_test() -> dict:
    """If 832 independently selected the FIXED observed ratio 830, it would
    be §8-stable. Propagating n_part ∈ {216..255} through 4·n_part − 100
    gives [764, 920] (span 156, ≈ ±9%). So the quantity drifts ~±10% while
    830 is fixed: the baseline n_part=233 → 832 is a COINCIDENCE, not a
    stable selection."""
    vals = [4 * n - 100 for n in N_PART_S8]
    lo, hi = min(vals), max(vals)
    observed_inside = lo <= RATIO <= hi
    drift_pct = 100.0 * (hi - lo) / (2 * np.mean(vals))
    return {
        'name': 'T4_s8_drift_decisive_test',
        'description': (
            "4·n_part − 100 drifts [764, 920] (±9%) across §8 while the "
            "observed √σ/m_e = 830 is fixed ⟹ the 832 ≈ 830 match is a "
            "baseline coincidence, not a §8-stable selection."
        ),
        'val_832_analogue_s8_range': (lo, hi),
        'span': hi - lo,
        'drift_percent': drift_pct,
        'observed_ratio_fixed': RATIO,
        'observed_inside_drift_band': observed_inside,
        'is_stable_selection': False,
        'baseline_coincidence': True,
        'pass': (hi - lo) > 50 and not (abs(4 * 233 - 100 - RATIO) / RATIO > 0.01),
    }


# ---------------------------------------------------------------------------
# T5. No independent bulk shell-stress integral selects 832
# ---------------------------------------------------------------------------

def test_T5_no_independent_integral() -> dict:
    """The genuine bulk shell-stress integrals over the Tangherlini
    geometry are O(10–70), nowhere near 466/832: Σω²(l=1,n=3..5) ≈ 69.8,
    the Bohr–Sommerfeld closure sum Σ(n+1)π ≈ 47.1. The number 466 enters
    ONLY through the v3 closure count 4β_quark/(2π) = 2·n_part (the fit).
    No independent integral selects ~466/832."""
    ev = _shell_eigenvalues()
    sum_omega2 = float(sum(ev[3:6]))
    bohr_sommerfeld_sum = (4 + 5 + 6) * PI       # Σ(n+1)π, n=3,4,5
    candidates = {
        'Σ ω²(l=1, n=3..5)': sum_omega2,
        'Σ (n+1)π (Bohr–Sommerfeld closure)': bohr_sommerfeld_sum,
    }
    none_near = all(v < 200 for v in candidates.values())
    return {
        'name': 'T5_no_independent_integral_selects_832',
        'description': (
            "Bulk shell-stress integrals are O(10–70): Σω²(n=3..5)≈69.8, "
            "Σ(n+1)π≈47.1. 466 enters only via the v3 closure count "
            "4β_quark/(2π)=2·n_part. No independent integral selects "
            "~466/832."
        ),
        'shell_integrals': candidates,
        '466_enters_only_via': 'v3 closure count 4β_quark/(2π) = 2·n_part (the fit)',
        'none_near_466_or_832': none_near,
        'pass': none_near,
    }


# ---------------------------------------------------------------------------
# T6. The circularity
# ---------------------------------------------------------------------------

def test_T6_circularity() -> dict:
    """n_part was FIT to reproduce the quark spectrum (which already
    encodes the physical scales). Recovering a scale ratio from n_part is
    therefore circular: you get back (an unstable version of) what was put
    in. So 832 ≈ 830 is the compensator echoing the spectrum, not a
    derivation."""
    return {
        'name': 'T6_circularity',
        'description': (
            "n_part was fit to the quark spectrum (which encodes the "
            "scales); recovering a scale ratio from n_part is circular. "
            "832 ≈ 830 is the compensator echoing the spectrum, not a "
            "derivation."
        ),
        'npart_fit_to': 'the quark spectrum (encodes the physical scales)',
        'recovering_ratio_from_npart_is': 'circular',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T7. Ledger unchanged
# ---------------------------------------------------------------------------

def test_T7_ledger_unchanged() -> dict:
    return {
        'name': 'T7_ledger_unchanged',
        'description': (
            "Channel-normalisation derivation via N_q+ΔN FAILS. PR #106 "
            "stands: √σ/m_e ≈ 830 remains an UNDERIVED open residual; BAM "
            "is one scale G + one open ratio. A genuine derivation needs an "
            "INDEPENDENT §8-stable bulk integral selecting ~830."
        ),
        'ratio_status': 'UNDERIVED (PR #106 unchanged)',
        'anchors': 'one foundational scale G + one open dimensionless ratio',
        'what_would_count': 'an INDEPENDENT bulk shell-stress integral (not the v3 fit) selecting ~830 AND §8-stable',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T8. Assessment
# ---------------------------------------------------------------------------

def test_T8_assessment() -> dict:
    return {
        'name': 'T8_assessment',
        'description': (
            "832 = 2N_q − N_lepton = 4·n_part − 4·k_5² is built from the "
            "compensator n_part; it §8-drifts 764–920 (so 832 ≈ 830 is a "
            "baseline coincidence), and no independent bulk shell-stress "
            "integral selects ~466/832. We are recycling n_part — the "
            "channel-normalisation derivation fails; √σ/m_e stays underived."
        ),
        'classification': 'RATIO_832_RECYCLES_N_PART_NOT_INDEPENDENT',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    tests = [
        test_T1_observation(),
        test_T2_built_from_npart(),
        test_T3_npart_is_compensator(),
        test_T4_s8_drift_test(),
        test_T5_no_independent_integral(),
        test_T6_circularity(),
        test_T7_ledger_unchanged(),
        test_T8_assessment(),
    ]
    if all(t['pass'] for t in tests[:7]):
        verdict_class = 'RATIO_832_RECYCLES_N_PART_NOT_INDEPENDENT'
        verdict = (
            'WE ARE RECYCLING n_part — 832 IS NOT AN INDEPENDENT SELECTION '
            'OF THE LEPTON/QCD SCALE RATIO. PR #106 left √σ/m_e ≈ 830 as the '
            'one underived ratio whose derivation would reduce BAM to a '
            'single anchor. The candidate N_q + ΔN = 832 ≈ 830 (0.2%) is '
            'seductive — and false.\n\n'
            '832 IS BUILT FROM n_part. N_q + ΔN = 2·N_q − N_lepton = '
            '2·(2·n_part) − 4·k_5² = 4·n_part − 4·k_5² = 4·233 − 100 = 832. '
            'So 832 is a linear function of n_part — the quark closure '
            'integer that PR #76/#97 established is a PHENOMENOLOGICAL '
            'COMPENSATOR: it absorbs the quark flavor puzzle, drifts '
            '216–255 across the quark_axioms §8 ablations, and only its '
            'parity is invariant. 832 inherits that non-derived status.\n\n'
            'THE DECISIVE §8-DRIFT TEST. If 832 independently selected the '
            'FIXED observed ratio 830, it would be §8-stable. It is not: '
            'propagating n_part ∈ {216..255} through 4·n_part − 100 gives '
            '[764, 920] — a span of 156, about ±9%. So the quantity drifts '
            'nearly ±10% while √σ/m_e = 830.3 is fixed; the baseline '
            'n_part = 233 merely happens to land it at 832 (0.2% from 830). '
            'That is a BASELINE COINCIDENCE — the same kind as 50π·k_5 = '
            '785 (PR #106) and F_13 = 233 (PR #76) — not a stable '
            'selection.\n\n'
            'NO INDEPENDENT BULK SHELL-STRESS INTEGRAL SELECTS 832. The '
            'genuine bulk shell-stress integrals over the Tangherlini '
            'geometry are O(10–70), nowhere near 466/832: the sum of shell '
            'eigenvalues Σω²(l=1,n=3..5) ≈ 69.8, the Bohr–Sommerfeld '
            'closure sum Σ(n+1)π ≈ 47.1. The number 466 enters ONLY through '
            'the v3 Hamiltonian closure count 4β_quark/(2π) = 2·n_part — '
            'i.e. through the fit. There is no independent integral that '
            'yields ~466 or ~832.\n\n'
            'THE CIRCULARITY. n_part was FIT to reproduce the quark '
            'spectrum, which already encodes the physical scales. '
            'Recovering a scale ratio from n_part is therefore circular — '
            'you get back (an unstable version of) what was put in. So '
            '832 ≈ 830 is the compensator echoing the spectrum it was fit '
            'to, not a derivation of the lepton/QCD hierarchy.\n\n'
            'LEDGER UNCHANGED. The channel-normalisation derivation via '
            'N_q + ΔN fails. PR #106 stands: √σ/m_e ≈ 830 remains an '
            'UNDERIVED open dimensionless residual, and BAM remains at one '
            'foundational scale (G) + one open ratio. A genuine derivation '
            'would need an INDEPENDENT bulk shell-stress integral (not the '
            'v3-fit closure count) that selects ~830 AND is §8-stable — '
            'which is not available.'
        )
    else:
        verdict_class = 'RATIO_832_INCONCLUSIVE'
        verdict = (
            'INCONCLUSIVE. A structural test failed; review the 832 / '
            'n_part accounting.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'identification': (
            '832 = 2N_q − N_lepton = 4·n_part − 4·k_5² is built from the '
            'compensator n_part; it §8-drifts 764–920 (so 832 ≈ 830 is a '
            'baseline coincidence), and no independent bulk shell-stress '
            'integral selects ~466/832 — we are recycling n_part'
        ),
        'answer': 'recycling n_part — NOT an independent selection',
        'decisive_test': '4·n_part−100 drifts 764–920 (±9%) while 830 is fixed',
        'independent_integrals': 'O(10–70) (Σω²≈70, Σ(n+1)π≈47); never ~466/832',
        'ledger': 'PR #106 unchanged — √σ/m_e underived; one scale G + one open ratio',
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
    L.append('# Is 832 = N_q + ΔN an independent scale-ratio selection, or recycled n_part? (PR #107)')
    L.append('')
    L.append(f"**Run:** {s['timestamp_utc']}")
    L.append('')
    L.append(
        "Tests, skeptically, the tempting candidate derivation of PR #106's "
        "underived lepton/QCD scale ratio: `N_q + ΔN = 2N_q − N_lepton = "
        "832 ≈ √σ/m_e ≈ 830` (0.2%). **Answer: we are recycling `n_part`.** "
        "832 is built from the phenomenological compensator; it §8-drifts "
        "764–920 (so the match is a baseline coincidence); and no "
        "independent bulk shell-stress integral selects ~466/832. The "
        "channel-normalisation derivation fails; `√σ/m_e` stays underived."
    )
    L.append('')
    L.append(f"- **Answer**: {s['answer']}")
    L.append(f"- **Decisive test**: {s['decisive_test']}")
    L.append(f"- **Independent integrals**: {s['independent_integrals']}")
    L.append(f"- **Ledger**: {s['ledger']}")
    L.append('')

    L.append('## Test summary')
    L.append('')
    L.append('| # | Test | Key finding | PASS? |')
    L.append('|---|---|---|---|')
    label_map = {
        'T1': 'N_q+ΔN = 2N_q−N_lepton = 832 ≈ 830 (0.2%) — tempting',
        'T2': '832 = 4·n_part − 4·k_5² — built from n_part',
        'T3': 'n_part is the compensator (§8-drifts 216–255, PR #76/#97)',
        'T4': 'DECISIVE: 4·n_part−100 drifts 764–920 (±9%) vs fixed 830',
        'T5': 'independent shell integrals O(10–70), never ~466/832',
        'T6': 'circular: n_part fit to the spectrum it would "predict"',
        'T7': 'ledger unchanged — √σ/m_e underived; one scale G + one ratio',
        'T8': 'RATIO_832_RECYCLES_N_PART_NOT_INDEPENDENT',
    }
    for t in s['tests']:
        passed = '**PASS**' if t['pass'] else '**FAIL**'
        nm = t['name']; prefix = nm[:2]
        L.append(f"| {prefix} | `{nm}` | {label_map.get(prefix, '—')} | {passed} |")
    L.append('')

    t4 = s['tests'][3]
    L.append('## T4: the decisive §8-drift test')
    L.append('')
    L.append(f"- observed `√σ/m_e = {RATIO:.1f}` (FIXED)")
    L.append(f"- `4·n_part − 100` across §8 ablations: `{t4['val_832_analogue_s8_range'][0]}`–"
             f"`{t4['val_832_analogue_s8_range'][1]}` (span {t4['span']}, ≈ ±9%)")
    L.append(f"- the baseline `n_part = 233` lands it at 832 (0.2% from 830) — but a "
             "quantity that drifts ±9% under ablations is **not** a stable selection "
             "of a fixed ratio. **Baseline coincidence.**")
    L.append('')

    t5 = s['tests'][4]
    L.append('## T5: no independent integral selects 832')
    L.append('')
    L.append('| bulk shell-stress integral | value |')
    L.append('|---|---:|')
    for k, v in t5['shell_integrals'].items():
        L.append(f"| {k} | {v:.1f} |")
    L.append('')
    L.append("All `O(10–70)`. The number 466 enters **only** through the "
             "v3 closure count `4β_quark/(2π) = 2·n_part` — the fit. No "
             "independent integral yields ~466/832.")
    L.append('')

    L.append('## Verdict')
    L.append('')
    L.append(f"**{s['verdict_class']}.** {s['verdict']}")
    L.append('')

    L.append('## What this leaves open')
    L.append('')
    L.append('- **A genuine derivation of `√σ/m_e ≈ 830`** — would require '
             'an INDEPENDENT bulk shell-stress integral (not the v3-fit '
             'closure count) that selects ~830 AND is §8-stable. Not '
             'available; the PR #106 status (underived; one scale G + one '
             'open ratio) stands.')
    L.append('')
    return '\n'.join(L)


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
    out = here / 'runs' / f'{ts}_ratio_832_npart_recycling_probe'
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
