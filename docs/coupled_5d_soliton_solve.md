# The coupled 5D+soliton solve: the confrontation, the bound, and the NR target (PR #203)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. This PR executes the final register
> item of the mass-ladder thread: couple the exact 5D suppression law
> (#202) to the #180 soliton energetics and confront the prediction.
> The refutation edge was live, and it **fired at the weak-field
> level** — the result is a quantified over-prediction, a clean
> negative on closability within the weak-field model, and a single
> falsifiable numerical-relativity target. No new fitted numbers appear
> anywhere in this PR.

## 0. The setup: a prediction with no knobs left

#202 established, exactly: `ε₁ = r_s/σ_mode` (the matching constant is
zero), hence

```
m_e/m_μ  =  (3/7) · (r_s / σ_mode)        (convention A).
```

The observed ratio requires `σ_mode/r_s = 88.6` (206.8 in convention B).
The coupled solve's job: extract **both scales from the one locked
#180 ψ–Φ–q solution** — σ_mode from the wave, r_core from the ordered
throat core — with every definition carried as a band and nothing fit.

## 1. The scales, extracted from the locked solution

From `relax(3.5, 0.05)` (the #180 fixed point, untouched):

| quantity | definition | value |
|---|---|---:|
| σ_mode | RMS of ψ | 1.273 |
| σ_mode | R* (the #201 kernel inversion, conv A/B) | 5.02 / 5.59 |
| r_core | q half-central radius | 0.833 |
| r_core | ψ² = ρ_c (the #178/#179 ordering threshold) | 1.083 |
| — | Φ(0) | −4.22 |
| — | M(ρ > ρ_c) | 1.65 |

## 2. The confrontation: the edge fired

| σ_mode | r_core | ratio | needed (A / B) |
|---|---|---:|---:|
| RMS | ρ_c radius | 1.18 | 88.6 / 206.8 |
| RMS | q-half radius | 1.53 | |
| R* (A) | ρ_c radius | 4.63 | |
| R* (A) | q-half radius | 6.02 | |
| R* (B) | q-half radius | 6.71 | |

**The weak-field coupled solve does not land the electron mass ratio.**
With the pairing-relevant definitions (R* over the q-core band) the
ratio is 4.6–6.7 versus the required 88.6 — a gap factor of **13–19**
(conv A; up to ~45 across the definition band). Equivalently: the
weak-field solve predicts `m_e/m_μ ≈ 0.071`, an **over-prediction of
m_e by ≈ ×15**. Stated plainly, as the result.

## 3. The bound, and the clean negative on weak-field closability

Two structural facts organize the failure:

1. **The direction is right: the weak-field value is an upper bound.**
   The physical claim of the whole arc is that the true 5D core is the
   *strong-field* endpoint of the #179 runaway (the q self-gravity
   channel that, above threshold, drives the core toward collapse) —
   necessarily *smaller* than the weak-field q-core. Since
   m_e ∝ r_core, the weak-field solve bounds m_e **from above**, and
   the observed value lies on the allowed side. The bound holds; the
   value does not land.
2. **The gap is not closable inside the weak-field model — measured.**
   Sweeping the binding strength (the soliton mass knob at locked
   couplings): the ratio moves the **wrong way** —

   | M | RMS | r_q(half) | RMS/r_q | Φ(0) |
   |---:|---:|---:|---:|---:|
   | 2.75 | 1.584 | 0.500 | 3.17 | −2.58 |
   | 3.50 | 1.273 | 0.833 | 1.53 | −4.22 |
   | 4.50 | 1.038 | 0.833 | 1.25 | −6.87 |

   Stronger binding compacts the wave *faster* than the core: the
   ratio **decreases** toward strong binding. No tuning of the
   weak-field family reaches 88.6. This is a genuine, cleanly
   established negative result — the missing factor is *physics absent
   from the model*, not a corner of its parameter space.

(Honest aside: Φ(0) = −4.2 at the locked point — the "weak-field" label
is already strained; the Newtonian model is being used at the edge of
its validity, which independently argues that the strong-field solve is
not optional.)

## 4. The named resolution and the falsifiable NR target

The missing physics has a name and a location in the repo: **the
strong-field core contraction** — the #179 runaway branch, where the
order field's self-gravity overwhelms the quartic saturation and the
core collapses toward the true 5D Tangherlini horizon r_s ≪ r_q. The
weak-field q-core is the *seed* of that collapse, not its endpoint.

The register item therefore resolves into **one number with a pass/fail
window**:

> **The NR target.** The full 5D numerical-relativity core solve must
> yield a core contraction `r_q(weak)/r_s(true) ≈ 13–45` (the
> convention/definition band) for the pairing mechanism to land
> m_e/m_μ. If NR instead gives an O(1) contraction, **the mouth-pairing
> mechanism is refuted as the quantitative origin of the electron
> mass** — the smallness mechanism (#195/#201/#202: index protection,
> multiplicative structure, naturalness) would survive, but its
> numerical anchor would not.

Everything upstream of this number is now exact or measured: the index
protection (#195), the multiplicative chain (#201), the suppression law
with its constants (#202), and the weak-field bound with its
non-closability (this PR). The chain from geometry to m_e/m_μ is
complete up to a single dimensionless output of a well-posed (if hard)
GR computation.

## 5. Honest scope

- All numbers from the locked #180 parameters and the locked #201/#202
  conventions; definition dependence carried as bands; **no new fits**.
- r_q(half) saturates at 0.833 for M ≥ 3.5 — partly a grid effect (the
  q-core is ~100 grid points, but its half-radius moves in coarse
  steps); the band [0.83, 1.08] absorbs this.
- The M-sweep probes binding strength at fixed couplings; a full
  weak-field parameter sweep would widen the family but cannot change
  the sign of the trend (compacting the wave compacts it faster than
  the threshold core — a structural feature of the ψ²>ρ_c core
  definition).
- The strong-field statement (the #179 runaway as the contraction
  mechanism) is a *location*, not a computation — precisely what the NR
  target formalizes.

## Reproduce

```bash
python -m experiments.closure_ledger.coupled_5d_soliton_solve_probe
# Verdict: WEAK_FIELD_COUPLED_SOLVE_OVERPREDICTS_M_E_BY_15X_GAP_NOT_CLOSABLE
#          _IN_WEAK_FIELD_THE_NR_CORE_CONTRACTION_13_TO_45X_IS_THE_TARGET
```
