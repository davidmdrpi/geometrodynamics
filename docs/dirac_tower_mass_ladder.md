# The Dirac-tower mass ladder: un-dialing the electron (PR #201)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. This PR executes the first item of
> the #200 open-items register: rebuild the lepton mass ladder with the
> electron level on the index-protected Dirac zero mode (#195), with
> the mouth coupling computed from the throat overlap machinery
> (#185/#190) — and re-run the three #192/#194 fine-tuning diagnostics
> on the rebuilt model to demonstrate that the dial is gone.

## 0. The result, up front

The surrogate's electron level was a **dialed cancellation**
(#194: `6.8754 − 6.6758 = 0.1996`, Barbieri–Giudice Δ = 74.7,
sign-flipping under ±2%). The rebuilt electron level is a
**multiplicative chain**

```
m_e-level  =  ε₁ · o₁ · S₁ ,
```

- `o₁ = 1` — the two mouths' zero-mode pairing overlap (#195; exact
  for k = 1, both zero modes are the base constants);
- `S₁ = (3/7)·m_μ-level` — the k = 1 bare scale tied to the k = 3 scale
  by the **Dirac tower ratio** m₁/m₃ = (3/2)/(7/2) (#197, closed form;
  a stated convention — the band across conventions is reported);
- `ε₁` — the **mouth coupling**: the tunneling amplitude of the
  k = 1 winding mode between the throat's two mouths.

Fitting μ/e = 206.77 fixes `ε₁ = (7/3)/206.77 = 0.01129`
(convention A; 0.00484 in convention B) — and the point is what this
number *is*:

- **the WKB neck aspect:** ε₁ = e^{−ℓ/a} gives ℓ/a = 4.48 (A) – 5.33
  (B): a throat neck a few radii long — an **O(1) geometric aspect
  ratio**, not a tuned cancellation;
- **the soliton overlap:** inverting the actual #185 kernel
  K(R*) = ε₁ on the #180 self-gravitating throat-soliton gives
  R* = 5.0 (A) – 5.6 (B) code units = **3.9 – 4.4 × the soliton RMS**
  (1.274) — the mouth separation is a few soliton radii. The required
  coupling sits exactly where the repo's own GR overlap machinery
  puts a few-radii separation.

Neither inversion is claimed as a derivation of 206.77 (the
anti-rigging rule): the claim is **structural** — the electron's
smallness is now a *product of geometric O(1-to-few) factors*, with no
subtraction anywhere in the chain.

## 1. The rebuilt model

- **The electron sector is superselected.** The k = 1 level sits on the
  index-protected zero mode: its bare energy is *exactly* 0
  (Atiyah–Singer, re-verified at 10⁻¹⁰), the one-mouth additive lift is
  *forbidden by angular momentum* (#195: no opposite-chirality
  j = q−½ partner; the D²₋ gap = 2 re-verified), and first-order mixing
  with the k = 3, 5 sectors is forbidden by winding-charge
  superselection. **The surrogate's dangerous t₁₃ = 14.85 — the term
  that produced the cancellation — has no counterpart.** The only mass
  channel is the two-mouth pairing: multiplicative by structure.
- **The heavy sector stays the surrogate's natural dynamics.** With the
  k = 1 row decoupled, the heavy block is the surrogate's 2×2 {k=3,5}
  block (base action + resistance + pinhole + uplift + t₃₅ transport).
  Removing the k = 1 repulsion shifts the μ level, so the uplift is
  honestly refit to τ/μ = 16.82 (β: 50π → ≈ 0.79·50π — the uplift
  remains the fitted dynamical knob it always was; #194 showed its
  sensitivities are natural, Δ < 1).
- **Calibration:** overall scale to m_e, c = −ln ε₁ to μ/e, β to τ/μ —
  three fitted numbers for three masses, same count as before. The
  *content* is not the fit; it is what the three diagnostics below say
  about the fitted point.

## 2. The un-dialing, measured (the #192/#194 diagnostics re-run)

| diagnostic | surrogate (#192/#194) | rebuilt (this PR) |
|---|---|---|
| Barbieri–Giudice Δ(m_e), worst | **74.7** (transport −57, base +31, slope +31) | **≈ 4.5** (= c; all other entries ≤ 1) |
| sign stability of m_e | flips under **±2%** transport | **positive identically** (a product of positive factors); ±25% moves it multiplicatively ×0.33–×3.1 |
| Berger λ-sensitivity, d ln m_e/dλ at round | **−70.9** = −1/(1−λ_break) | **+4.15** (= c − 1/3; the fiber-riding neck: ε₁(λ) = e^{−c/λ}) |
| λ_break (zero crossing) | **0.98598** (a 1.4% squash) | **none on (0, ∞)** — m_e(λ) > 0 identically |
| Monte Carlo P(m_e ≤ observed), ±25% priors | 7.7% — a linear-measure **sliver** | **O(50%)** — the observed value is *typical*: ln m_e is smoothly distributed (width c/4 ≈ 1.1), no cliff, no sign flips anywhere in the prior box |

The #194 fine-tuning — an unprotected cancellation, exactly as rare as
generic linear measure — is replaced by a technically natural product:
small because a tunneling exponent is a few, stable because nothing can
cancel.

## 3. The clean impossibility: the hierarchy is NOT mouth pairing

Could the inter-generation ratios also come from the pairing? **No —
provably, at the structural level.** Pairing-generated masses scale as
ε_k = e^{−k·c} (the winding-k mode decays k times faster through the
neck): *decreasing* in k. Reproducing μ/e = 206.77 from pairing alone
would need ε₃/ε₁ ≈ 88.6 — a tunneling amplitude that *grows* with the
barrier. Tunneling amplitudes cannot grow with the barrier:

> **The electron's smallness and the inter-generation hierarchy have
> different origins.** Smallness = index protection + mouth pairing
> (geometric, natural — this PR). Ratios = the dynamical uplift (the
> #122/#136 phase budget — fitted, and *naturally* fitted: its
> sensitivities were < 1 already in #194).

This is the same division of labor #193/#197 found spectrally
(structure kinematic, hierarchy dynamical), now realized in the mass
model itself: after the rebuild, **every Barbieri–Giudice sensitivity
in the ladder is O(few) or below** — the ladder is fully natural, with
the hierarchy parametrized but not tuned.

## 4. Honest scope

- The mouth coupling is *constrained*, not derived: ε₁ is fit to μ/e
  and then shown to correspond to O(1) geometry (neck aspect 4.5–5.3;
  mouth separation 3.9–4.4 soliton radii on the actual #180 profile).
  Deriving ℓ/a from the 5D throat solution is the remaining step of
  this item (it requires the #199-register 5D core dynamics).
- The S₁ convention (tying the k=1 bare scale to k=3 by the Dirac-tower
  ratio) is stated and the convention band is carried through every
  number; no conclusion depends on the choice.
- The heavy sector remains the surrogate's fitted dynamics (uplift
  refit β ≈ 0.79·50π after the k=1 decoupling); its naturalness is
  #194's result, re-verified here on the rebuilt model.
- ε₁(λ) = e^{−c/λ} (fiber-riding neck) is the declared Berger map for
  the electron factor, in the spirit of #192's ingredient maps; the
  conclusion (no zero crossing, O(1) log-sensitivity) is robust to the
  map's details because the factor is a positive exponential.

## Reproduce

```bash
python -m experiments.closure_ledger.dirac_tower_mass_ladder_probe
# Verdict: ELECTRON_LEVEL_REBUILT_ON_THE_INDEX_PROTECTED_ZERO_MODE
#          _FINE_TUNING_REMOVED_HIERARCHY_REMAINS_DYNAMICAL_AND_NATURAL
```
