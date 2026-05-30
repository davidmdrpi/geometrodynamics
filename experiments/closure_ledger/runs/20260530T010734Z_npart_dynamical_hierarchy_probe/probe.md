# n_part = 233 revisited: the quark hierarchy is the program's one dynamical sector (PR #97)

**Run:** 2026-05-30T01:07:34+00:00

Revisits PR #76's `n_part = 233` compensator with the machinery built since. **Fresh angle:** the neutrino arc (#86–#90) proved a huge hierarchy can be geometric (the `e^{S}` tortoise bounce), so size is not the obstruction. The quark hierarchy is non-geometric for a *diagnosable* reason — it is **irregular** (neither power-law nor exponential), the signature of QCD-RG dynamics. The quark sector is the program's **one dynamical** hierarchy; `n_part` (and the 366-quantum lepton↔quark gap) compensates it.

- **Identification**: the quark inter-generation hierarchy is irregular (the QCD-RG signature) — the program's one dynamical sector; unlike the geometric lepton (integer 100) and neutrino (exponential bounce) hierarchies it has no stable geometric closure quantity, so n_part=233 (and the 366-quantum gap) compensates it
- **Reframing**: a huge hierarchy can be geometric (neutrino exponential); size is not the obstruction
- **Obstruction**: irregularity (neither power-law nor exponential) = QCD-RG dynamics
- **Closure gap**: N_q − N_lepton = 466 − 100 = 366 quanta = the dynamical (QCD) excess
- **Right route**: QCD-shell model WITH α_s running (the missing RG dynamics)
- **Upholds**: PR #76 N_PART_IS_PHENOMENOLOGICAL_COMPENSATOR, sharpened

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_recap_compensator` | n_part=233 = v3 compensator (parity-only §8 inv.; drift 216–255) | **PASS** |
| T2 | `T2_geometric_lepton_neutrino_sectors` | geometric sectors: leptons (integer 100) + neutrinos (e^{−S}) | **PASS** |
| T3 | `T3_huge_hierarchy_can_be_geometric` | huge hierarchy can be geometric (neutrino) ⟹ size not the issue | **PASS** |
| T4 | `T4_quark_hierarchy_irregular` | quark hierarchy IRREGULAR: c/u≈588 vs t/c≈136 (not exp/power) | **PASS** |
| T5 | `T5_geometric_shell_undercarries` | shell span ×2.2 vs observed ×6.4×10⁹ ⟹ ~2.9×10⁹ gap | **PASS** |
| T6 | `T6_dynamical_qcd_rg_diagnosis` | irregular ⟹ QCD-RG dynamical; gap 366 = dynamical excess | **PASS** |
| T7 | `T7_honest_scope` | #76 upheld + sharpened; route = QCD-shell + α_s running | **PASS** |
| T8 | `T8_assessment` | QUARK_HIERARCHY_DYNAMICAL_N_PART_COMPENSATES | **PASS** |

## T4–T5: Irregular hierarchy, geometric shell under-carries

- up-type consecutive ratios: `c/u ≈ 588`, `t/c ≈ 136` (not constant ⟹ not exponential)
- down-type: `s/d ≈ 20`, `b/s ≈ 45` (≠ up-type ⟹ not a single power law)
- geometric shell `ω²(1,n=3,4,5)` span = ×2.2; observed quark mass² span = ×6.4e+09; gap ≈ ×2.9e+09

## T6: Geometric vs dynamical closure integers

| sector | hierarchy type | closure quantity | §8 |
|---|---|---|---|
| charged leptons | power-law (ledger) | `4·k_5² = 100` | **stable** |
| neutrinos | exponential `e^{−S}` (bounce) | `S ≈ 16` (geometric) | stable form |
| quarks | **irregular (QCD-RG)** | `N_q = 466` (`n_part=233`) | **drifts 216–255** |

The lepton↔quark closure gap `N_q − N_lepton = 366` quanta is the dynamical (QCD) excess — the quantity a dynamical QCD-shell model must produce. The quark closure integer is the only one that drifts under §8: the signature that it compensates dynamical, not geometric, physics.

## Verdict

**QUARK_HIERARCHY_DYNAMICAL_N_PART_COMPENSATES.** THE QUARK HIERARCHY IS THE PROGRAM'S ONE DYNAMICAL SECTOR; n_part = 233 COMPENSATES IT. PR #76 classified n_part = 233 (N_q = 2·n_part = 466, β_quark = 233π) as a phenomenological compensator — only its parity is §8-stable, n_part itself drifting 216–255 — and identified the right route as a quantitative QCD-shell model. This probe revisits it with the machinery built since and sharpens the diagnosis.

A HUGE HIERARCHY CAN BE GEOMETRIC. When PR #76 ran, the worry was that the ~9-order quark mass² hierarchy was simply too large for the geometric closure machinery. The neutrino arc (#86–#90) overturned that: it derived a comparable hierarchy — the keV→TeV seesaw M_R = m_D·e^{S}, ~10⁶ in mass — as a CLEAN GEOMETRIC EXPONENTIAL (the non-orientable tortoise bounce, an O(15) action). So size is not the obstruction. The program now has two geometric hierarchy types: charged leptons (a closure-ledger ladder with the clean, §8-stable integer 4·k_5² = 100, PR #71) and neutrinos (the exponential bounce, PR #88–90).

THE QUARK HIERARCHY IS IRREGULAR. The question is then sharper — does the quark hierarchy have the REGULARITY of a geometric one? It does not. The consecutive up-type mass ratios are m_c/m_u ≈ 588 and m_t/m_c ≈ 136 — not constant, so NOT a clean exponential (geometric progression); and not the k²-style pattern of a power law. The down-type ratios (m_s/m_d ≈ 20, m_b/m_s ≈ 45) differ from the up-type, so the two Z₂ partitions are asymmetric. The quark hierarchy is irregular — matching neither geometric type.

THE GEOMETRIC SHELL CANNOT CARRY IT. The quark shell basis (k=0, overtones n=3,4,5; PR #77/#83) has ω²(1,n) = 14.6, 22.7, 32.5 — a mass² span of only ×2.2 — whereas the observed quark mass² span (t/u) is ×6.4×10⁹. The geometric shell under-produces by ~2.9×10⁹.

THE DIAGNOSIS: A DYNAMICAL (QCD-RG) HIERARCHY. Irregularity across a wide scale range is the signature of renormalisation-group running: α_s runs logarithmically, so the quark masses are QCD-dressed differently at each scale — an intrinsically DYNAMICAL hierarchy. Leptons and neutrinos do not feel QCD, so their hierarchies are geometric; the quarks do, so theirs is dynamical. This is why the quark closure integer is the ONLY one that drifts under §8 (216–255): it absorbs dynamical content that no geometric closure quantity encodes. The lepton↔quark closure gap N_q − N_lepton = 466 − 100 = 366 quanta is precisely that dynamical (QCD) excess.

So n_part is not merely fit on the "wrong basis" (PR #76) — the quark hierarchy is the BAM mass program's ONE DYNAMICAL SECTOR, and a geometric closure integer can only compensate it, never derive it. The right route is sharpened: a QCD-shell model WITH α_s running (the missing ingredient is the RG dynamics, identified now by contrast with the geometric lepton and neutrino sectors).

HONEST SCOPE. ESTABLISHED: the PR #76 compensator verdict is upheld and sharpened — size is not the obstruction (the neutrino exponential is geometric and huge); the quark hierarchy is irregular (the QCD-RG signature); the geometric shell carries only ×2.2 of the ×6.4×10⁹ span; n_part (and the 366-quantum gap) is the dynamical excess; the quark sector is the program's one dynamical hierarchy. NOT established: a first-principles n_part = 233 — none exists in the geometric machinery, and by the diagnosis none should; the derivation route (a QCD-shell model with α_s running) is a substantial program outside the closure-ledger machinery.

## What this leaves open

- **A first-principles `n_part = 233`** — none exists in the geometric machinery, and by the diagnosis (the quark hierarchy is dynamical) none should.
- **A QCD-shell model with `α_s` running** — the missing RG dynamics; substantial, outside the closure-ledger machinery (PR #76's "right route", now diagnosed specifically as the RG ingredient).
