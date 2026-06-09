# One-loop photon vacuum polarisation and the running of α (PR #144)

> **Framing (to avoid a category error).** QFT on the *fixed classical* throat
> geometry — geometry → fields, **not** quantum gravity. The photon self-energy
> is computed on the classical antipodal cavity; the metric stays a classical
> input.

PRs #141–#143 completed the gauge–matter **structure**: the minimal coupling and
the Σl-even vertex (#141), the Ward identity / current conservation / photon
masslessness — all structurally (#142), and the α ledger separating the derived
EM structure from the one input, the value α (#143). Both #142 and #143 left the
same flagged open item: the **running of α** was classified as "derived" but
never **computed**. Meanwhile the matter sector already has its one-loop
two-point function (the self-energy Σ, #136); the photon self-energy Π — the
vacuum polarisation — is the one missing one-loop object. This PR computes it.

## Why this probe, now

The one-loop ledger has three corners: the matter propagator correction
(Σ, #136), the vertex (#137/#141), and the gauge propagator correction (Π).
Two of three existed; the Ward identity (#142) ties all three together but was
verified only structurally. Computing Π closes the one-loop two-point sector
and upgrades the #142 Ward statement from structural to quantitative.

## The polarisation bubble

The charged-pair loop over the antipodal cavity modes (#135): pair (n, m) opens
at the threshold `s_nm = (ω_n + ω_m)²`, with photon–pair vertex
`v_nm = ∫ φ_γ ψ_n ψ_m dr*` (the #137/#141 triple overlap with one photon leg)
and spectral density `ρ_nm = c_nm |v_nm|² ≥ 0`. The angular part carries the
#141 antipodal Z₂ selection rule — the photon couples only to **even-Σl** pair
channels (re-verified exactly via the S³ monomial integral).

## The cavity Ward identity, computed

Gauge invariance under minimal substitution `p → p − c₁A` forces the O(A²)
energy shift to vanish: the diamagnetic (seagull) `+1` cancels the paramagnetic
(current–current) sum,

    1 − S = 0,    S = 4 Σ_{m≠n} |⟨m|∂|n⟩|² / (E_m − E_n),

which is the Thomas–Reiche–Kuhn sum rule `Σ (E_m − E_n)|x_mn|² = 1` in disguise
(`p_mn = ±(E_m − E_n) x_mn / 2`). On the cavity's Dirichlet matter tower the
cancellation is verified **numerically to ~3e-5** — the quantitative face of the
#142 Ward identity. Consequence: `Π(0) = 0` — no photon mass, the photon pole
stays exactly at `q² = 0`, and the `1/q²` kernel (#42–#44) is protected through
one loop.

## The absorbing counterfactual

With an absorbing throat the matter modes are complex (#130) and charge leaks
(#142): the pair thresholds move off the real axis, `Im Π ≠ 0` at all real `s`
(the photon acquires an absorption width below every pair threshold — computed:
`−0.042` vs the antipodal `−6e-09`), and the real-mode orthonormality enforcing
the Ward cancellation is gone. Gauge protection of the massless photon
**requires** the unitary antipodal throat — the one-loop face of #129/#142.

## Screening and the running

The Ward-protected (once-subtracted) polarisation is the dispersion sum
`Δ(Q²) = Σ ρ_nm Q²/(s_nm(s_nm + Q²))`: manifestly ≥ 0 and monotone increasing in
spacelike `Q²`, so `α_eff = α/(1 − Δ)` **increases** with `Q²` — the QED
screening direction, with the discrete pair thresholds the cavity analogue of
the lepton thresholds in the running. Feeding the **same** dispersion machinery
the flat-space 4D pair density `ρ_QED(s) = (α/3π)√(1−4m²/s)(1+2m²/s)` reproduces
the textbook log running with slope `dΔ/d ln Q² = α/3π` (verified numerically to
0.97% over three decades). The running's form, sign, and coefficient follow from
the Ward-protected spectral representation; the boundary value `α(μ₀)` stays the
one EM input (#143) — this probe deliberately does **not** hunt for 137 (the
#107/#108 anti-numerology discipline).

## Scope and epistemic ledger

- **Derived:** the Ward cancellation (computed), photon masslessness / `1/q²`
  protection, spectral positivity, no width below the lowest pair threshold,
  monotone screening, and the flat-limit log coefficient `α/3π`.
- **Modelled:** the photon-leg radial profile (the soft cavity photon mode —
  the same posture as the #136 cubic vertex).
- **Input:** the boundary value `α(μ₀) ≈ 1/137` (#143, the 137 problem).
- **Open:** higher loops; the full 4D tensor `Π^μν` beyond the partial-wave
  scalar; the absolute normalisation (#133); the flavor residuals (#134).

## Reproduce

```bash
python -m experiments.closure_ledger.vacuum_polarization_running_probe
# Verdict: VACUUM_POLARIZATION_WARD_MASSLESS_SCREENING_LOG_RUNNING_ALPHA_INPUT
```
