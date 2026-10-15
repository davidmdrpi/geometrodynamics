# Resonance-breaking ladder of the LRS circle family: results

Date: 2026-10-10.

- **Specification and code:** [`4c924cc`](r3_ladder_prereg.md), pushed 18:27:19 UTC, before any scan or λ evaluation at any rung.
- **Archives:** `experiments/closure_ledger/runs/20261010_r3_ladder/`. Each binds the SHA-256 of all ten sources. They were committed one by one as written, without being read.
- **Replay:** `python -m experiments.closure_ledger.r3_ladder_replay --full` returns VERIFIED.
  - All registered labels re-derive.
  - Every 10th high-precision solution of every rung re-closes with the unchanged 160-bit map to ≤ 5.6e-36.
  - Tamper tests are in `tests/test_r3_ladder.py`.
- **Post-hoc diagnostics:** `experiments/closure_ledger/r3_ladder_posthoc.py` writes `posthoc.json`. It was written after the labels were fixed and is descriptive only.

## Verdict: informative negative

The weak resonance breaking in the LRS sector is **ordinary**.

Counted with the half-clock order Q_h (P = h∘h; the P-resonance p/q is h's (p+q)/(2q)), every rung of the family breaks the way an analytic near-integrable twist map does:
- The lowest-order rung, 3/7 (Q_h = 7), breaks at 5.6e-9. The uniform-radius extrapolation from 2/5, registered in advance, predicted 6.5e-9.
- Every rung shows exactly the harmonic structure that the half-clock order implies.
- The breaking strengths across Q_h = 7…24 (Λ from 5.6e-9 down to 5.6e-27) follow a smooth per-harmonic decay law.

**The "extraordinarily weak breaking" was a counting error.** Mine: I treated the 2/5 resonance as order 5 when it is order 10. There is no evidence for a hidden integral, and I recommend dropping near-integrability of the LRS sector as a lead.

## 1. Registered labels

| label | value |
|---|---|
| **primary (rung 3/7)** | **ORDINARY_BREAKING**: Λ = 5.61e-9, inside the window [2.4e-12, 6.5e-7] |
| S1, the 2/5 signal | **SIGNAL_CONFIRMED**: Λ = 1.54e-11, double-precision error 1.7e-12, dominant harmonic 10 |
| S2, harmonic selection | **HALF_MAP_SELECTION**: dominant harmonic = Q_h at all 7 rungs |
| S3, exponent | **NEITHER_FITS**: RMS 1.12 decades with exponent Q_h, 4.02 with exponent q (threshold 1 decade) |
| S4 | **BREAKING_RESOLVED**: all 7 rungs resolved, at r = 1e-30 |

**Numerical gates.**
- All 420 phases (60 per rung) converged in every stage.
- High-precision noise ν ≤ 1.1e-35 per rung, so the 1e-30 floor sets the resolution.
- Chord Newton took 5–6 iterations everywhere, with no high-precision Jacobian fallback.
- Scan condition numbers ≤ 1.6e4; seed distances ≤ 5.1e-3.

**S3 failed its registered criterion.** The pre-specified uniform-radius law misses by 1.12 decades RMS against a 1-decade threshold. It does fit 3.6× better with exponent Q_h than with q.

## 2. Per rung

The M1 and M2 predictions were registered in advance (§3 of the specification).

| P-rung | a | Q_h | Λ (48-digit) | M1 | M2 | dominant harmonic | sign changes | double-precision Λ | double-precision error |
|---|---|---|---|---|---|---|---|---|---|
| 5/11 | .1484 | 11 | 6.83e-16 | 1.3e-14 | 4.5e-15 | 11 | 22 | 2.3e-12 | 2.3e-12 |
| 4/9 | .1712 | 18 | 1.57e-22 | 2.8e-19 | 4.3e-23 | 18 | 36 | 2.3e-12 | 2.3e-12 |
| **3/7** | .2022 | **7** | **5.61e-9** | 2.4e-10 | **6.5e-9** | 7 | 14 | 5.61e-9 | 1.5e-12 |
| 5/12 | .2233 | 24 | 5.64e-27 | 4.1e-21 | 8.9e-28 | 24 | 48 | 1.4e-12 | 1.4e-12 |
| 2/5 | .2486 | 10 | 1.54e-11 | anchor | anchor | 10 | 20 | 1.57e-11 | 1.7e-12 |
| 3/8 | .2839 | 16 | 2.07e-15 | 3.1e-14 | 4.3e-17 | 16 | 32 | 2.2e-12 | 2.2e-12 |
| 4/11 | .2982 | 22 | 1.07e-18 | 4.8e-17 | 9.3e-23 | 22 | 44 | 3.4e-12 | 3.4e-12 |

**Harmonic structure.** At every rung:
- the sign-change count is exactly 2Q_h;
- harmonic Q_h is 4.2–5.9× the largest other harmonic.

These are complete chains of the half-clock map h at its own resonance, Q_h elliptic plus Q_h hyperbolic points. No rung shows a further symmetry (a dominant harmonic at a multiple of Q_h).

**Double precision resolves only two rungs.** The double-precision scan (the r3_breaking machinery) resolves only 3/7 and, marginally, 2/5. Everywhere else its λ is pure map error of 1.4–3.4e-12. That error is the size r3_breaking estimated (ν = 2.8e-12), so its noise model was right.

## 3. What the ladder shows (post hoc, descriptive)

Write the breaking as Λ = exp(−Q_h σ). For an analytic family, σ is the width of the analyticity strip of the dynamics in the circle angle. The order-Q resonant coefficient decays like e^{−Qσ}.

| rung | 5/11 | 4/9 | 3/7 | 5/12 | 2/5 | 3/8 | 4/11 |
|---|---|---|---|---|---|---|---|
| a | .148 | .171 | .202 | .223 | .249 | .284 | .298 |
| σ = −ln Λ / Q_h | 3.18 | 2.79 | 2.71 | 2.52 | 2.49 | 2.11 | 1.88 |

**σ decreases monotonically with amplitude.** The linear fit is σ ≈ 4.22 − 7.53a, with RMS 0.089 per harmonic. That is ≤ 0.85 decades in Λ even at Q_h = 22. This is the textbook picture for a near-integrable analytic twist map: the strip narrows as the invariant circles grow toward breakdown.

A linear extrapolation puts σ = 0 near a ≈ 0.56. That is a rough indication of where the family stops being regular, not a measurement.

**Why S3 failed.** The registered uniform-radius model is the special case σ = ln(R/a), with R constant.
- It holds for a ≤ .25: R = 2.8–3.6.
- It fails at the two largest circles: R = 2.35 at 3/8 and 1.96 at 4/11. These break more strongly, as the strip narrows faster near breakdown.

**Implications for the earlier 2/5 result.**
- The 2/5 rung itself is unremarkable. Its σ lies on the same smooth curve as every other rung.
- Its high-precision λ reproduces the archived double-precision scan to 1.7e-12, with correlation 0.998.
- The post-hoc harmonic-10 signal reported in `docs/r3_breaking.md` was real.

## 4. What this does not establish

- **LRS only.** Nothing about the diagonal 6D family, inhomogeneous modes or n ≥ 3.
- **Not a proof of non-integrability.** The analytic-family picture fits but is not derived. It does, however, remove the reason to suspect a hidden integral.
- **The #319 control is not covered.** It is the 6D two-return loop, unbroken to about 1e-13 in r3_breaking. Under the same identity it lifts to an order-4 h-resonance, so its near-perfect closure is not explained by this study. Its amplitude and any D3 constraint need their own check. That is the one R3 near-integrability question still open.

## 5. Relation to PR #325

- **Priority.** #325 (frozen at `82c832f`, 18:16:46 UTC) derived the same half-clock identity and lift ten minutes before this freeze. I had not seen it, and priority for the identity is theirs.
- **Its prediction at 3/7.** #325 registered a held-out prediction at 3/7: dominant harmonic 7, 14 sign changes, and harmonic-7 amplitude ≥ 10r at double precision. This study measures the identical scan system there:
  - harmonic 7 at 1.9e-9;
  - 14 sign changes;
  - double-precision λ equal to the exact value to 1.5e-12.

  #325's prediction should therefore come out supported if its numerical gates pass.
- **Pending cross-check.** The 12 shared phases will be compared once #325's archive exists.

## 6. Disclosures

**Production was resumed after a stop.**
- The first production process started at 18:27 UTC. At 20:27 it was killed by the session's 2-hour limit on background commands, during the 4/11 scan.
- The six scans already written were kept as written.
- `experiments/closure_ledger/r3_ladder_resume.py` was committed (`59f4219`) before it was used. It calls only the frozen probe functions, one rung at a time, skipping existing archives, and changes no setting.
- Two resumption records are archived and hash-pinned in the replay: `resume_r1.json` (4/11 scan plus the high-precision stage) and `resume_r2.json` (noise stage plus scoring).
- No result had been read at any point before `result.json` was written.

**Registered vs post hoc.**
- The specification's M1/M2 anchors and the 3/7 window were fixed before measurement.
- The σ(a) analysis in §3 was not registered.
- The verdict rests on the registered primary, S1 and S2 labels. S3 failed its criterion and is reported as such.

**Draft history.** This is disclosed in the specification (§1).
- My first, uncommitted draft scored a double-precision ladder against a^q.
- The half-clock identity was found before any measurement, and the draft was replaced before the freeze.
- Under the corrected orders, double precision would have resolved only 3/7 and 2/5. The high-precision stage is what made the ladder informative.
