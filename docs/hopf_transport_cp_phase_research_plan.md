# Hopf-connection derivation of the quark CP phase (PR #158)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity. The CP phase
> is the Hopf-fiber transport phase of the same-partition shell couplings.

PR #156 calibrated a CP phase on the partition-mixing element and left the
db-triangle shape (β = 22°) as the acceptance test for a Hopf-derived phase.
This PR runs that test — and finds the mechanism itself must move,
**correcting #156** along the way.

## The correction to #156 (stated first)

The partition-mixing route's apparent J was an **artifact**:

| φ_q(k) form | J raw | J_db raw | non-unitarity | J unitarized |
|---|---|---|---|---|
| linear k | 7.7e-6 | −3.2e-9 | 0.158 | ≈ 0 |
| winding k² | ~7e-6 | ~−3e-9 | 0.146 | ≈ 0 |
| Casimir k(k+2) | ~5e-6 | ~−2e-9 | ~0.15 | ≈ 0 |

The charged-current CKM is non-unitary at the 16% level (the u–d
near-degeneracy amplifies the cross-partition leakage), the quartet Jarlskog
invariants disagree by ×1000 (a unitary CKM has one J), and the unitarized
core carries J ≈ 0 for **every** phase form. Independently, the #156 mixing
strength implies a first-row unitarity deficit ~×40 over experimental
bounds. **Partition mixing is doubly excluded as the CP origin.** (#156's
ceiling identity and J = 0 baseline are |V| arithmetic and stand; its
calibration interpretation is superseded.)

## The relocation

The locked same-partition coupling `−t·e^{−α·dk}·cos(phase·dk)` is the
**real part** of the Hopf transport factor `e^{iφ·dk}` — the handoff's own
placeholder note pointed here. The geometry forces the complexification with
**opposite orientation per partition class** (the #63 C-swap flips c₁):

    (H±)_{kk'} = −t·e^{−α·dk}·e^{±i·φ_h·dk},   dk = max(k, k')

Both blocks stay Hermitian ⟹ V = U₊†U₋ **exactly unitary** (6.7e-16) with
**quartet-consistent J** (|J₁₂ + J_db| = 2.4e-18): genuine CP. The locked
baseline is recovered exactly at φ_h = 0; π itself is the χ = 0 fiber
holonomy (`hopf/connection.py`: ∮A = π cos χ).

## One parameter, the full triangle

| quantity | calibrated (φ_h\* = 0.611, J-fit only) | pure π/k₅ = 0.6283 (no fit) | observed |
|---|---|---|---|
| J/target | 1.000 | 0.969 | 1.0 |
| β | 22.89° | 22.78° | 22.2° |
| γ | 65.79° | 63.48° | 65.9° |
| α | 91.32° | 93.75° | 91.9° |
| sin δ | 0.905 | **0.888** | 0.887 |
| max mass shift | 8.7e-4 | 9.0e-4 | < 0.016 |
| V_cb | 0.0377 | 0.0377 | (stiff pred) |

Calibrating to J **alone** predicts the entire triangle to ~1°; masses shift
0.09%, V_cb is untouched, and V_us moves toward the data (0.112 → 0.123).
**The #156 acceptance test (β = 22°) is passed.**

## The π/k₅ candidate (flagged per #107/#108)

The calibration sits 2.7% from **π/k₅** — the χ = 0 fiber holonomy π divided
by the k₅ = 5 winding quanta. The pure, uncalibrated value reproduces
**five CP observables from zero free parameters** (table above; sin δ to
0.1%). Per the anti-numerology discipline this is a **candidate** closure
identification — a clean one-step combination of two derived quantities, no
ad-hoc factor, tested on observables it was not fitted to — pending an
independent transport derivation (the #152 modelled→derived path).

## Budget impact

No input consumed; if π/k₅ is confirmed, the #156-consumed input (the CP
phase content) is **returned** — quark CP fully derived. Today: the input
downgrades to a one-parameter family with a principled candidate at 2.7%. A
new independent bound emerges: partition mixing ε ≲ 0.004 from first-row
CKM unitarity.

## Reproduce

```bash
python -m experiments.closure_ledger.hopf_transport_cp_phase_probe
# Verdict: QUARK_CP_FROM_HOPF_TRANSPORT_TRIANGLE_PREDICTED_PI_OVER_K5_CANDIDATE
```
