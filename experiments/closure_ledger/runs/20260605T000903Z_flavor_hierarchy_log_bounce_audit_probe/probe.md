# Flavor hierarchy audit from logarithmic throat bounce lengths (PR #134)

**Run:** 2026-06-05T00:09:03+00:00

Audits whether the three-generation flavor hierarchy follows from the logarithmic throat bounce lengths `L*(ε) = (rs/2) ln(1/ε)` (#88/#132). The log-bounce gives a power-law mass `m ∝ ε^p`, which governs the chargeless neutrino sector in form and ordering but overshoots its value; the charged-lepton (winding) and quark (cavity) hierarchies are different mechanisms.

- **Mechanism**: log-bounce S = c(rs/2)ln(1/ε) ⟹ m ∝ ε^p (power law)
- **Neutrino**: the log-bounce sector: m_ν ∝ ε^{4.8} (#112), right ordering, value overshoots ×28 (#113)
- **Charged leptons**: winding ladder β·k² (#71), irregular ln(m) — NOT log-bounce
- **Quarks**: cavity overtones / n_part (#77–#80), irregular ln(m) — NOT log-bounce
- **Why residual**: m ∝ ε^p hypersensitivity (∂ln m/∂ln ε = p) ⟹ flavor values residual (#108)
- **Open**: neutrino value overshoot (#113); charged/quark irregular magnitudes (flavor puzzle #108)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | audit the flavor hierarchy from logarithmic throat bounce lengths | **PASS** |
| T2 | `T2_log_bounce_power_law_mass` | log-bounce S = c(rs/2)ln(1/ε) ⟹ m ∝ ε^p (power law, verified) | **PASS** |
| T3 | `T3_neutrino_sector_is_log_bounce` | neutrino = log-bounce: m_ν ∝ ε^{4.8}, ordering right, value ×28 over | **PASS** |
| T4 | `T4_charged_leptons_winding_not_bounce` | charged leptons = winding (#71), NOT log-bounce (ln m irregular) | **PASS** |
| T5 | `T5_quarks_cavity_not_bounce` | quarks = cavity / n_part (#77), NOT log-bounce (ln m irregular) | **PASS** |
| T6 | `T6_log_bounce_hypersensitivity` | hypersensitivity ∂ln m/∂ln ε = p ⟹ flavor values residual (#108) | **PASS** |
| T7 | `T7_scope` | scope: three-mechanism structure; log-bounce governs ν form only | **PASS** |
| T8 | `T8_assessment` | FLAVOR_HIERARCHY_LOG_BOUNCE_GOVERNS_NEUTRINO_SECTOR_FORM_NOT_VALUE | **PASS** |

## The neutrino sector: log-bounce amplification (#113)

`m_ν ∝ ε^{4.8}` (#112), `ε_n ∝ 1/χ_n` (#79). `χ_n = [0.304, 0.097, 0.039]` ⟹ `ε_n` ratios `[1.0, 3.13, 7.79]` (normal ordering: True). The steep power amplifies a ×2 ε spread into `2^4.8 ≈ 27.9×` in mass — the ×28 (= 2^4.8) overshoot of #113.

## Why the flavor values are residual: hypersensitivity

| p | ∂ln m/∂ln ε | mass spread for ×2 ε |
|---:|---:|---:|
| 1.0 | 1.0 | 2.0× |
| 4.8 | 4.8 | 27.9× |

`m ∝ ε^p` makes masses hypersensitive to the throat penetration depth, so the flavor values' irreducibility (#108) is a consequence of the exponential mass-action relation — not a separate mystery.

## Verdict

**FLAVOR_HIERARCHY_LOG_BOUNCE_GOVERNS_NEUTRINO_SECTOR_FORM_NOT_VALUE.** THE LOGARITHMIC THROAT BOUNCE GOVERNS THE NEUTRINO FLAVOR SECTOR — ITS FORM AND ORDERING, NOT ITS VALUE — AND IS NOT THE FLAVOR HIERARCHY'S UNIFYING MECHANISM. PRs #88/#132 made the bounce action logarithmic, S = c(rs/2)ln(1/ε); this probe audits whether the three-generation mass hierarchy follows from those logarithmic bounce lengths.

THE MECHANISM: LOG-BOUNCE ⟹ POWER-LAW MASS. With S = c·L*(ε) = c(rs/2)ln(1/ε), a tunnelling mass is m = m_0 e^{−S} = m_0 ε^{c rs/2} = m_0 ε^p — the logarithm turns the exponential into a PURE POWER LAW in the throat penetration depth ε (verified e^{−cL*} = ε^p). Masses are powers of ε, not exponentials of a linear quantity.

THE NEUTRINO SECTOR IS THE LOG-BOUNCE SECTOR. The neutrino is the only genuine tunnelling sector — k = 0, chargeless, the neck not EM-propped (#86/#88) — so m_ν = m_D e^{−S} ∝ ε^p with p ≈ 4.8 (#112). The generation healing lengths ε_n ∝ 1/χ_n (#79/#113) give the right ORDERING (normal), but the steep power amplifies the modest χ_n spread: a ~×2 spread in ε becomes 2^4.8 ≈ 28× in mass — the ×28 overshoot of #113. So the log-bounce governs the neutrino hierarchy's form and ordering, with the value residual.

CHARGED LEPTONS AND QUARKS ARE NOT LOG-BOUNCE. Charged leptons are Dirac (c₁ = ±1, neck EM-propped, no tunnelling, #86/#88); their masses come from the winding ladder β·k² (#71), and their ln(m) is irregular (gen-differences 5.33, 2.82 — ratio 0.53). Quarks are shell-resolving cavity overtones (#77–#80, the n_part sector); their ln(m) is irregular too (up-type 6.37, 4.91 — ratio 0.77). So the flavor hierarchy is a THREE-mechanism structure — bounce (ν), winding (charged leptons), cavity (quarks) — not a single log-bounce phenomenon.

WHY THE FLAVOR VALUES ARE RESIDUAL. The power law m ∝ ε^p has ∂ln m/∂ln ε = p, so with p a few the masses are HYPERSENSITIVE to the throat penetration depth: a ×2 ambiguity in ε becomes 2^p ≈ 28× in mass. The irreducibility of the flavor values (#108) is therefore a CONSEQUENCE of the exponential mass-action relation — the log-bounce amplification — not a separate mystery: even a fully-determined geometry would leave the masses hypersensitive to the ε residual (#112).

SCOPE. An audit/classification: it establishes the log-bounce ⟹ power-law mechanism, locates it in the neutrino sector (form & ordering), separates the charged-lepton (winding) and quark (cavity) sectors, and explains the flavor values' residual-ness as log-bounce hypersensitivity. It does NOT solve the flavor hierarchy, predict any mass value, or remove the universal flavor puzzle (#104/#108).
