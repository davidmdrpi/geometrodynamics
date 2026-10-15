# Specification: the resonance-breaking ladder of the LRS circle family, at 48-digit precision

Date: 2026-10-10. Parent: `claude/geometrodynamics-qft-audit-vpktax` at `9015ad8`.

This document and its code are frozen in one commit, published before any
phase scan or λ evaluation at any rung below. The 2/5 rung's λ is the
archived double-precision one from `docs/r3_breaking.md`. Corrections will be
dated notes only.

## 1. Why this study, and what changed before the freeze

`docs/r3_breaking.md` found the LRS resonant circle at rotation 2/5
(a ≈ .249) intact to Λ ≤ 1.6e-11. A post-hoc harmonic-10 signal of about
6e-12 was reproduced by three map evaluations. My summary compared Λ with
a^5 ≈ 1e-3, which suggested about eight orders of unexplained suppression.

**Finding before measurement: the section map is a square.** The flow
commutes with S: (q, q') → (−q, −q'):
- q enters only through q², and linearly in q'' = −G q;
- in the full model S is −1 ∈ O(4) acting on the quartet.

S maps upward clock crossings to downward ones and fixes
z = (A, p_A, x, p_x). So P = h∘h, where h is the half-return map from the
downward section to the next upward crossing, read in z.

**Verification** (new module `geometrodynamics/waves/lrs_taylor.py`, test
points off the circle family):
- |P − h∘h| ≤ 3e-40;
- return times add to 1e-41.

On the archived circle a = .2348 (rotation ρ = .40966), h advances the
circle by .70485 turns. The two candidates are ρ/2 = .2048 and
ρ/2 + 1/2 = .70483, so ρ_h = ρ/2 + 1/2.

**Consequence.** The P-resonance p/q is the h-resonance (p+q)/(2q), whose
order Q_h is the reduced denominator: Q_h = q if p and q are both odd, and
2q otherwise.
- At 2/5 this gives h-resonance 7/10, of order 10. An h-chain of order 10
  gives λ(φ) exactly the 20 sign changes and dominant harmonic 10 that were
  observed. The post-hoc signal is the expected chain signature, not an
  unexplained selection.
- The generic yardstick at 2/5 is therefore a^10 ≈ 9e-7, not a^5. The
  apparent suppression shrinks from about 8 to about 5 orders, and is
  equivalent to an effective analyticity radius R ≈ 3 in units of a.
- Whether the weak breaking is anomalous at all is now open.

**Draft history (disclosed).** My first, uncommitted draft scored a
double-precision ladder against a^q.

**Why that draft was replaced.** It had two flaws:
- **Wrong orders:** it counted the resonance order as q instead of Q_h.
- **Predetermined null:** under a uniform-suppression model anchored at 2/5,
  every other rung was predicted to break below the double-precision floor,
  so its null outcome was fixed in advance.

This specification replaces it.

## 2. Question

With orders counted correctly, is the breaking along the family:
- ordinary for an analytic map (consistent with the 2/5 level), or
- anomalously suppressed?

**The decisive rung is 3/7.** It is the lowest-order resonance in the
family, h-order 7 at a = .2022. Anchored at 2/5, both models predict breaking
at 2e-10 to 7e-9 there. That is above the double-precision floor and fifteen
orders above the high-precision floor.

## 3. Rungs and anchored predictions

The rungs are every P-rational with q ≤ 12 between the archived rotation
numbers .4847 and .3587. The amplitude a and action I are interpolated
linearly in ω between the bracketing archived circles.

The models are anchored at Λ* = 1.57e-11 (archived, 2/5) with exponent Q_h:
- **M1, uniform coefficient:** Λ = ε a^{Q_h}, with ε = 1.74e-5.
- **M2, uniform radius:** Λ = (a/R)^{Q_h}, with R = 2.99.

| P-rung | a | I | ρ_h | Q_h | a^q | a^{Q_h} | M1 | M2 |
|---|---|---|---|---|---|---|---|---|
| 5/11 | .1484 | .0315 | 8/11 | 11 | 7.7e-10 | 7.7e-10 | 1.3e-14 | 4.5e-15 |
| 4/9 | .1712 | .0419 | 13/18 | 18 | 1.3e-7 | 1.6e-14 | 2.8e-19 | 4.3e-23 |
| **3/7** | .2022 | .0582 | 5/7 | **7** | 1.4e-5 | 1.4e-5 | **2.4e-10** | **6.5e-9** |
| 5/12 | .2233 | .0703 | 17/24 | 24 | 1.5e-8 | 2.4e-16 | 4.1e-21 | 8.9e-28 |
| 2/5 | .2486 | .0871 | 7/10 | 10 | 9.5e-4 | 9.0e-7 | anchor | anchor |
| 3/8 | .2839 | .1121 | 11/16 | 16 | 4.2e-5 | 1.8e-9 | 3.1e-14 | 4.3e-17 |
| 4/11 | .2982 | .1233 | 15/22 | 22 | 1.7e-6 | 2.8e-12 | 4.8e-17 | 9.3e-23 |

**Counting the order as q gives different magnitudes and a different
ordering.** Anchored the same way, the q-ordering predicts 3/8 ≥ 3/7. The
Q_h-ordering predicts 3/7 ≫ 3/8 by 4 to 8 orders.

## 4. Method

**Stage `scan` (double precision, frozen machinery).**
- Uses `r3_breaking.scan_point` on the unchanged `r3_breaking_probe.P4`
  (`esu_map`, DOP853, 1e-12/1e-14).
- Has q nodes and 60 node-0 phases φ_j = 2πj/60. The seeds are
  K*(φ + 2πi·p/q) on the interpolated circle K*, with unfolding along
  g = Ωᵀc'(φ) and a phase row.
- Is accepted at residual ≤ 1e-11.
- For 2/5 this recomputes the archived scan.

**Stage `hp` (config C1).** Every converged scan point is re-solved by chord
Newton with the Taylor-series map `lrs_taylor.hp_map`:
- 160-bit mpfr, order 40, step ρ·(1e-40)^{1/40};
- the residual is computed in mpfr and the correction with the archived
  double Jacobian;
- it switches once to a 160-bit finite-difference Jacobian (h = 1e-20) if
  an iteration contracts by less than 1e3;
- the target is residual < 1e-35, in at most 14 iterations.

The scan system (g, t0, c0 and the slice) is the archived double one, read
as exact binary numbers. So the hp stage solves exactly the same equations
to roughly 23 more digits.

**Stage `hpnoise` (config C2).** Each C1 solution is continued with the
200-bit, order-50, 1e-50 map to residual < 1e-45. The noise is
ν = max_j |λ_C2 − λ_C1| + max_j (both residuals). The resolution is
r = max(10ν, 1e-30).

**Validation already done** (`tests/test_r3_ladder.py`; no ladder rung is
evaluated):
- **Isotropic sector:** return time π to 1e-36. A matches mpmath's
  independent Taylor solution of A'' = −A + A³ to 1e-30.
- **Agreement:** with `esu_map` to ≤ 1e-11; between C1 and C2 to ≤ 1e-37.
- **Constraint:** ≤ 1e-38.
- **Structure:** P = h∘h to 1e-37; symplectic to 1e-12 (finite-difference
  Jacobian).
- **Toy map (gmpy2 port of `r3_breaking.toy_map`):**
  - the unbroken case gives λ < 1e-38;
  - the broken case matches the double λ and agrees across C1 and C2 to
    1e-37;
  - the Jacobian fallback converges;
  - the full pipeline gives the analytic Λ and harmonic 5.

## 5. Registered labels

**Per rung.**
- **INDETERMINATE:** more than 6 of the 60 phases fail in any stage.
- **RESOLVED:** Λ = max_j |λ_C1| ≥ 10r.
- **UNRESOLVED:** otherwise. Its upper bound is 10r.

For each rung we record:
- Λ, ν and r;
- the double-precision error max_j |λ_double − λ_C1|;
- sign changes;
- the harmonic spectrum of λ_C1 (all 60 phases);
- log10 deviations from M1 and M2.

**Primary label (rung 3/7).** The window is [10⁻² · min(M1, M2),
10² · max(M1, M2)] = [2.4e-12, 6.5e-7].

| label | condition |
|---|---|
| **ORDINARY_BREAKING** | resolved and Λ within the window |
| **ANOMALOUS_SUPPRESSION** | the upper bound (Λ if resolved, 10r if not) is below 2.4e-12 |
| **ENHANCED_BREAKING** | resolved and Λ > 6.5e-7 |
| **INCONCLUSIVE** | otherwise, or INDETERMINATE |

**Secondary labels.**

**S1, the 2/5 signal at high precision.**

| label | condition |
|---|---|
| SIGNAL_CONFIRMED | resolved, with 5e-12 ≤ Λ ≤ 5e-11, double-precision error ≤ 5e-12 and dominant harmonic 10 |
| SIGNAL_ARTEFACT | upper bound ≤ 1.6e-13 |
| OTHER | otherwise |

**S2, harmonic selection** (resolved rungs with a full spectrum; at least
two are needed, else UNTESTED). Here Q_h is the predicted dominant harmonic.

| label | condition |
|---|---|
| HALF_MAP_SELECTION | every dominant harmonic equals Q_h |
| EXTRA_SELECTION | all are multiples of Q_h, at least one larger (a further hidden symmetry) |
| VIOLATED | otherwise |

Multiples above 30 alias, so EXTRA is detectable only for Q_h ≤ 15.

**S3, exponent.** Over the resolved rungs (at least four, else UNTESTED),
fit log10 Λ − n·log10 a = α − n·log10 R by least squares, for n = Q_h and
for n = q. Compare the RMS residuals in decades, RMS_h and RMS_q.

| label | condition |
|---|---|
| EXPONENT_QH | RMS_h ≤ 1 and RMS_h ≤ RMS_q/2 |
| EXPONENT_Q | the mirror condition |
| NEITHER_FITS | both RMS > 1 |
| UNDISCRIMINATED | otherwise |

**S4.** INTEGRABLE_TO_HP_RESOLUTION if no rung is resolved, otherwise
BREAKING_RESOLVED.

## 6. Interpretation map (stated before measurement)

| outcome | reading and next step |
|---|---|
| ORDINARY_BREAKING (with S3 EXPONENT_QH) | The "extraordinarily weak breaking" was mostly a miscount of the resonance order. Once the half-map is accounted for, the family breaks like an ordinary analytic map with R ≈ 3. Near-integrability should be dropped as a lead, and I will say so plainly. |
| ANOMALOUS_SUPPRESSION | Suppression beyond the half-map order at the lowest-order rung. This is evidence for a hidden approximate integral. Next step: construct it from the hp data (for example, a fitted invariant on the 48-digit orbits). |
| ENHANCED_BREAKING | 2/5 is the anomalously quiet rung, so look for a rung-specific mechanism. |
| S2 EXTRA_SELECTION at any rung | A further symmetry exists beyond P = h∘h. |
| S4 INTEGRABLE_TO_HP_RESOLUTION | Would contradict S1 and the three-integrator agreement in r3_breaking. It would be checked before any claim. |

**Priors.**
- Primary: ORDINARY 0.6, ANOMALOUS 0.3, ENHANCED 0.1.
- S1 CONFIRMED: 0.85.
- S2 HALF_MAP_SELECTION: 0.85.
- S3: EXPONENT_QH 0.55, EXPONENT_Q 0.1, rest 0.35.

## 7. What will not follow
- **LRS only.** Nothing about the diagonal 6D family, inhomogeneous modes or
  n ≥ 3.
- **No integral is identified.**
- **M1 and M2 are not normal-form computations.** They are one-parameter
  extrapolations from one anchor, so the window is two decades either side.
- **No action selection.** Closure at several rationals is expected either
  way.

## 8. Disclosures
- **Pre-freeze map evaluations:**
  - `lrs_taylor` at five off-family test points (timing: 0.14 s per
    evaluation at C1, 0.2 s at C2);
  - P and h at one point, c(0), of the archived circle a = .2348 (for ρ_h);
  - the toy maps.
- **No scan, Newton solve or λ** has been computed at any rung.
- **Brackets** come from the archive only.
- **gmpy2** (2.3.2, pip) was installed for this study and is listed as the
  optional extra `precision`. Results do not depend on mpmath's backend.
- **The 2/5 anchor is the archived value**, not the one re-measured here.
  The predictions in §3 are therefore fixed now.
- **Runtime estimate:** scan about 70 min, hp about 30 min, hpnoise about
  15 min, on 4 processes.
- **Archives** are append-only (the producer refuses to overwrite) and bind
  the SHA-256 of all sources. They will be committed as written. A replay
  with pinned hashes follows the run.
- **Correction note.** A dated note on `docs/r3_breaking.md` (half-map
  structure) is committed with this freeze.

## 9. Reproduction

    pip install -e ".[precision]"
    OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.r3_ladder_probe
    python -m pytest -q tests/test_r3_ladder.py
