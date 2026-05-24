# Throat-to-shell transition: leptons vs the QCD shell channel

Tests a physical hypothesis: the odd-`k` fermionic ladder (the charged
leptons, #67) does not stay a localized **pointlike throat** forever — as
energy/excitation rises, the mode delocalizes from the throat into a
**shell/ring** standing wave, leaving the localized charged-lepton
channel and entering a shell-coupled (QCD-like) channel. In the user's
framing: a low-energy **focused pulse** converges to the pointlike throat
(lepton); a high-energy **wavefront** spreads into a ring/shell (QCD).
This would also bear on the three-generation cutoff.

## What is tested

The radial overtone ladder at `l = 1` (the closure-ledger lepton reading
`e, μ, τ = n = 0, 1, 2`), extended to higher `n`, on the cavity
`r ∈ [R_MID, R_OUTER]` (shell thickness `ΔR = R_OUTER − R_MID`). For each
mode we measure where its probability `|u|²` lives:

  - `⟨r⟩ − R_MID` — mean displacement from the throat;
  - **throat fraction** — `|u|²` in the inner third (the throat zone);
  - **participation ratio** — `1 / (N Σ p²)` (`→ 2/3` for a uniform shell
    standing wave, smaller for a throat-focused mode);
  - the radial **wavelength** `λ = 2π/ω` relative to `ΔR` (focused pulse
    `λ ≫ ΔR` vs shell wavefront `λ → ΔR`).

## The finding (honest)

  - **The three leptons (n = 0, 1, 2) are the throat-localized end.** The
    localization metrics rise *monotonically* through them — the electron
    (`n=0`) is the most focused (`⟨r⟩−R_MID ≈ 0.021`, throat-frac `0.96`,
    `λ/ΔR ≈ 23` — a long-wavelength focused pulse converging to the
    pointlike throat); the muon and tau progressively delocalize
    (throat-frac `0.85`, `0.75`).

  - **From n ≈ 3 the modes saturate into shell standing waves.** `⟨r⟩`
    plateaus (`≈ 0.046`), the participation ratio reaches the
    uniform-standing-wave value `2/3`, and the wavelength approaches the
    shell scale `ΔR` — the delocalized **shell/ring** channel (the
    QCD-side candidate).

So the throat→shell transition is **confirmed as a real trend, and it
saturates right after the third generation** — the localized
charged-lepton throat ladder gives way to shell-coupled modes. Honest
caveat: it is a **saturating crossover, not a razor-sharp cutoff**; the
metrics plateau rather than hard-stop at `n = 2`. The exact
three-generation boundary therefore involves *both* this delocalization
crossover *and* the closure-quantum / β-uplift cutoff (the #67 follow-on)
— the two mechanisms are complementary, not redundant.

## Energy reading (the user's framing)

The wavelength `λ = 2π/ω` makes the focused-pulse vs wavefront picture
quantitative: the electron (`n=0`, lowest energy) has `λ ≈ 23 ΔR` — a
long-wavelength pulse focused onto the point throat; rising `n` (energy)
shrinks `λ` toward the shell scale `ΔR`, where the mode becomes a
multi-node wavefront resolving and filling the shell. Low energy →
pointlike throat (lepton); high energy → shell/ring wavefront (QCD).

## What it is and is not

  - **Is:** a confirmed delocalization trend (throat → shell) along the
    fermionic ladder, saturating after the three leptons; the localized
    charged-lepton modes and the shell-coupled modes are geometrically
    distinct channels of the same radial structure.
  - **Is not:** a hard three-generation cutoff (the crossover is smooth),
    nor a full identification of the shell modes with the QCD/quark
    spectrum (that needs the quark sector — `docs/quark_beta_status.md` —
    and is future work).

## B4 accounting

The localization metrics are **dimensionless ratios** (`⟨r⟩/ΔR`,
participation, `λ/ΔR`); the transition is a geometric/structural feature,
independent of the single anchor `m_e`. The mass *values* carry the
scale.

## Tests

  T1. **The overtone ladder + metrics.** `n = 0…8` at `l=1`: `⟨r⟩`,
      throat fraction, participation ratio, `λ/ΔR`.
  T2. **Leptons are the throat-localized end.** `n=0,1,2` (e,μ,τ):
      localization metrics rise monotonically (progressive
      delocalization); the electron is the most focused.
  T3. **Shell saturation from n≈3.** participation `→ 2/3`, `⟨r⟩`
      plateaus — shell-filling standing waves (the shell/ring channel).
  T4. **Focused-pulse vs wavefront.** electron `λ/ΔR ≈ 23` (focused →
      pointlike throat); rising `n` → `λ → ΔR` (wavefront → shell).
  T5. **Transition after the third generation.** the localization
      saturates right after `n=2` (τ).
  T6. **Honest assessment.** confirmed trend + saturating crossover (not
      a hard cutoff); complements the closure cutoff; shell↔QCD
      identification future work.
  T7. **Falsification / B4.** no delocalization would falsify; BAM shows
      a clear monotone trend. Metrics dimensionless/scale-independent.
  T8. **Assessment.**

## Verdict structure

  - **THROAT_TO_SHELL_TRANSITION_CONFIRMED** (expected): higher odd-`k`
    fermionic excitations *do* leave the localized charged-lepton throat
    channel and delocalize into shell/ring standing waves — the
    localization metrics rise monotonically through the three leptons
    (`n=0,1,2`, the focused throat end) and saturate from `n≈3` into
    uniform shell standing waves (participation `→ 2/3`), with the
    wavelength shrinking from `λ ≈ 23 ΔR` (electron, focused pulse) toward
    the shell scale. The transition is a real trend that saturates right
    after the third generation. Honest: a saturating crossover (not a
    razor-sharp cutoff); it complements the closure-quantum cutoff, and
    the shell↔QCD identification is future work.

  - **NO_TRANSITION**: the modes stay equally localized at all `n` — the
    hypothesis would be falsified.

## What this leaves open

  - **Shell ↔ QCD identification.** Whether the shell-saturated modes
    *are* the quark/QCD spectrum (the repo's quark sector,
    `docs/quark_beta_status.md`) — a full match is not done here.
  - **The sharp three-generation cutoff.** The crossover locates the
    transition but does not hard-cut at `k=5`; the closure-quantum /
    β-uplift cutoff (#67 follow-on) is the complementary piece.

## Cross-references

  - `docs/even_k_absence_research_plan.md` — the odd-`k` fermionic ladder
    (#67).
  - `docs/quark_beta_status.md` — the quark/QCD sector (the shell-channel
    candidate).
  - `geometrodynamics/tangherlini/radial.py` — `V_tangherlini`, the
    radial cavity.
  - `experiments/closure_ledger/throat_to_shell_transition_probe.py` —
    this probe.
