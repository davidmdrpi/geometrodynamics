# Resonance breaking at the LRS 2/5 crossing: results

Date: 2026-10-05.
- Specification and code: [`4aa4cd7`](r3_breaking_prereg.md), pushed 05:08:38 UTC, before any map evaluation at or near the 2/5 crossing and before any phase scan of the #319 loop.
- Archives: `experiments/closure_ledger/runs/20261005_r3_breaking/` (`scan_control.json`, `scan_main.json`, `noise.json`, `orbits.json`, `result.json`). Each binds the SHA-256 of all seven sources.
- Replay: `python -m experiments.closure_ledger.r3_breaking_replay --full` returns VERIFIED. Labels re-derive, and every converged scan point re-closes with the unchanged maps to ≤ 1.4e-12.
- Post-hoc diagnostics: `experiments/closure_ledger/r3_breaking_posthoc.py`, written to `posthoc.json` after the labels were fixed. Descriptive only.

## 1. Registered labels

| case | label | Λ = max\|λ\| | noise ν | resolution r | sign changes |
|---|---|---|---|---|---|
| control: #319 two-return loop (q = 2, 6D) | **UNBROKEN_LOOP** | 8.1e-14 | 9.7e-13 | 1.0e-11 | 6 |
| main: LRS circle at rotation 2/5 (q = 5, 4D) | **UNBROKEN_LOOP** | 1.57e-11 | 2.8e-12 | 2.8e-11 | 20 |

**Gates.**
- All 60 phases converge in both scans:
  - closure residuals ≤ 2.0e-12;
  - constraint ≤ 1.9e-12;
  - condition number ≤ 7.7e3.
- Main scan: seed distance ≤ 5.1e-3, from the interpolated circle, s = .3122.

**Main case.** UNBROKEN_LOOP is assigned because Λ ≤ r. The circle crosses 2/5 at a ≈ .249 with action **I = 0.087172**. That action is post hoc: it comes from a periodic spline through the converged nodes, and using half the points changes it by 1.6e-7.

**This was not my prior.** I registered BROKEN_CHAIN as the expected outcome for the main case.

## 2. What the main λ(φ) shows (post hoc, descriptive)

The registered rule applies a factor-10 margin above the noise. Under that margin, the main scan is classified UNBROKEN_LOOP. But λ(φ) is not structureless:
- it is a clean oscillation, dominated by harmonic 10, with amplitude 5.75e-12 and Λ = 1.57e-11;
- the next harmonics (8 and 12) are about 1.1e-12;
- harmonic 5 is 3.7e-13;
- there are 20 sign changes, a multiple of 2q.

A smooth systematic error of the map would produce the same signature. Summed over five nodes spaced 2π/5 apart, it keeps only harmonics that are multiples of 5. So λ was re-solved with three map evaluations:

| re-evaluation | phases | max\|λ\| | max\|Δλ\| vs archive | correlation |
|---|---|---|---|---|
| Radau, 1e-12/1e-14 (registered noise stage) | 6 | 1.4e-11 | 8.9e-13 (plus residual) | n/a |
| DOP853, tightened to 1e-13/1e-15 | 60 | 1.54e-11 | 1.5e-12 | 0.998 |
| fixed-step RK4, Richardson (2048, 4096), Newton return time | 12 | 1.35e-11 | 8.8e-13 | 0.999 |

The RK4 map shares no integration or event code with `esu_map`. It uses:
- `r3_extension.rhs`, not `conformal_rhs`;
- `r3_family`'s eigenbasis coordinate maps, not `esu_map`'s;
- no `solve_ivp` and no event location.

All three reproduce the pattern at the 1e-12 level. **The structure is therefore a property of the flow**, as far as these checks can tell. It is not an artefact of the integrator.

**Isolated orbits.**
- The orbit stage solved the unconstrained closure from all 20 brackets:
  - 14 converged to ≤ 1e-11;
  - 6 stalled between 1.5e-11 and 4e-8.
- The smallest singular values of the closure Jacobian are 6e-10 to 1e-6. Compare 3e-13 for #319's loop, which is a continuum.
- Residues from exact jets are about ±(1–8)e-10.
  - In 12 of 14 orbits the sign follows the crossing direction: negative at upward crossings, positive at downward.
  - That is the pattern of alternating hyperbolic and elliptic chain orbits.
  - Coarse and fine jets differ by up to about 1e-10.
- De-duplication at 1e-8 kept 14 "distinct" orbits. In a near-continuum, least-squares Newton lands at slightly different points, so that count is not meaningful.

**Reading.**
- **Main case.** The LRS 2/5 resonant circle is most likely broken, but only at the **~1e-11 level**: Λ/a ≈ 6e-11. The leading breaking harmonic is 10, so harmonic 5 is suppressed, which is consistent with two pairs of chain orbits. This lies below the registered resolution, so it does not change the label. It sits about 5× above the registered noise and is reproduced by independent maps.
- **Control.** The #319 loop shows nothing above about 1e-13. That tightens #319's continuum from its 1e-10 multiplier spread to a closure-obstruction bound of about 1e-13.

## 3. What this establishes

Within the exact homogeneous system:

1. **Resonant circles survive almost perfectly in a family with no known protecting symmetry.**
   - Large LRS amplitude: a ≈ .25, I ≈ .087.
   - Low order: q = 5.
   - Any breaking is ≲ 1.6e-11.

   A generic nonintegrable map would break this resonance at a level set by the low-order resonant harmonics, many orders of magnitude larger. The homogeneous dynamics is therefore **integrable or extremely close to integrable**, at least in the LRS sector, though the post-hoc analysis suggests it is not exactly integrable. The extra integral has not been identified.
2. **#319's constant multipliers are not evidence of a symmetry specific to the circular family.** Near-integrability of the homogeneous system explains them generically. §7c of the specification guessed at a Z3 selection rule for the circular family; it is not needed.
3. **Closure is not special to the 3/2 condition.** Closed period-5 histories exist at the LRS 2/5 crossing, at action I ≈ 0.0872. They form either a loop to ≤ 1.6e-11 or, post hoc, a chain of four orbits. Part D of `r3_extension` listed 2/5 only as a candidate; it is now verified. A closure condition fixes an action only once its rational is chosen. Closure as such does not select a unique action.

## 4. What it does not establish

- **The integral is not identified, and exact integrability is not shown.** The post-hoc signal points the other way, at the 1e-11 level.
- **Only two resonances were tested:** LRS 2/5, and the circular 1/2 as the control.
- **Nothing is claimed about inhomogeneous perturbations or the n ≥ 3 sectors.**
- **No action selection or quantisation.** Item 3 above points against selection by closure alone.

## 5. Consequence for the experiment proposed in #320

A nonlinear "symmetry-reduced neighbourhood" study was proposed. The family direction is not a geometric symmetry. To within about 1e-11 to 1e-13, it is the flow of an (approximate) extra integral. The natural nonlinear framework is therefore action–angle (normal-form) coordinates of a near-integrable system. Over any feasible horizon, a direct nonlinear evolution would show quasi-periodic motion, as integrability predicts, and could not discriminate anything.

The informative next step is to **identify the integral**:
- Fit the Birkhoff normal form of the diagonal section map to high order, building on #317's multimode normal-form machinery.
- Or test explicit candidate conserved quantities along orbits.

Either would show:
- whether the integral is exact or only formal;
- whether its breaking is the harmonic-10 term seen here;
- whether it carries any physical meaning for the program.

## 6. Disclosures

- **Two specification items worked as intended:**
  - the registered noise includes the final Radau residual, which dominated ν (2.8e-12) in the main case;
  - the factor-10 margin was fixed before the freeze.

  The harmonic-10 signal was found by inspecting λ after the labels were fixed. It is reported here only as a post-hoc observation.
- **Commit sequence.** The archives were committed stage by stage while the run progressed, without being read. The run used the frozen code unchanged.
- **Replay scaffold.** The replay scaffold and tamper tests were committed during the run, before any result was read. The archive hashes were pinned afterwards.
- **What the run used.** The run used only the frozen probe. The post-hoc script reuses the archived Jacobians and nodes (chord Newton) and changes only the map evaluation.

## Note (2026-10-10): the section map is a square

This note was added before the ladder study (`docs/r3_ladder_prereg.md`) measured anything.

**The structure.** The flow commutes with (q, q') → (−q, −q'). That symmetry maps downward clock crossings to upward ones and fixes z. So P = h∘h, where h is the half-return map. This is verified to 3e-40 with the 48-digit Taylor map `lrs_taylor`. On the archived circle a = .2348, h rotates by ρ/2 + 1/2.

**What it means for the 2/5 result.** The LRS 2/5 crossing is h's 7/10 resonance. Its leading resonant order is 10, not 5. The post-hoc observations in §2 are exactly the signature of that order-10 h-chain:
- harmonic 10 dominant;
- harmonic 5 near noise;
- 20 sign changes.

**What changes in §3 item 1.** The comparison with "low-order resonant harmonics" should use a^10 ≈ 9e-7, not a^5 ≈ 1e-3. The breaking (≲ 1.6e-11) is then about 5 orders below that yardstick, not about 8. That is equivalent to an effective analyticity radius of about 3 in units of a. Whether this is anomalous is the question the ladder study tests. The labels in §1 are unaffected.
