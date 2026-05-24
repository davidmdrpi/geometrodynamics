# Geometric gyromagnetic-ratio probe

**Run:** 2026-05-24T03:37:55+00:00

Derives the electron throat's gyromagnetic ratio g=2 from the BAM geometry — the throat's Pauli/SU(2) spinor structure (T=iσ_y) minimally coupled to the Hopf monopole (A_φ=½ cos χ) — and checks the Schwinger anomaly a=α/2π. Extends the spin thread (PR #60) from kinematics to the magnetic moment.

- **g factor**: 2.0
- **Origin**: `(σ·D)²=D²−eσ·B (Pauli/SU(2), T=iσ_y) + Hopf monopole (½)`
- **Thomas link**: g=2 ⟺ ω_a=0 (spin tracks momentum, #60)
- **Anomaly**: a=(g−2)/2=α/2π (one loop; tree g=2 geometric)
- **B4 caveat**: g, a dimensionless; μ_B=eℏ/2m carries the single anchor

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_pauli_algebra_g2` | (σ·D)²=D²−eσ·B; {σ,σ}=2δ → g=2 | **PASS** |
| T2 | `T2_hopf_monopole_g2` | Hopf monopole ½ → μ=(e/m)S, g=2 (g·s=1) | **PASS** |
| T3 | `T3_bmt_spin_tracks_momentum` | ω_a=(g/2−1)(eB/m); g=2 → spin tracks momentum | **PASS** |
| T4 | `T4_magnetic_moment_bohr_magneton` | μ = g·μ_B·s = μ_B (g=2, s=½) | **PASS** |
| T5 | `T5_schwinger_anomaly` | a=α/2π=0.0011614 (one loop) | **PASS** |
| T6 | `T6_falsification_criterion` | spinor → g=2 (not classical g=1); BAM passes | **PASS** |
| T7 | `T7_b4_accounting` | g, a dimensionless; μ_B carries the anchor | **PASS** |
| T8 | `T8_assessment` | g=2 geometric; spin tracks momentum; μ=μ_B | **PASS** |

## T1: g = 2 from the Pauli/SU(2) algebra

- Pauli identity (σ·a)(σ·b)=a·b+iσ·(a×b) holds: True
- anticommutator {σ_i,σ_j}=2δ_ij: True (the factor 2)
- (σ·D)²=D²−eσ·B; σ·B carries σ=2S → g = 2

## T2: g = 2 from the Hopf monopole

- Hopf monopole charge (A_φ/cos χ) = 0.5000 (= spin ½)
- minimal coupling → μ=(e/m)S → g = 2
- g·s = 1.0000 → magnetic moment = 1 Bohr magneton

## T3: BMT — g = 2 ⟺ spin tracks momentum

| g | ω_a = (g/2−1)(eB/m) | spin tracks momentum |
|---:|---:|:---:|
| 0.0000 | -1.0000 | False |
| 1.0000 | -0.5000 | False |
| 2.0000 | +0.0000 | True |
| 2.0023 | +0.0011 | False |

At g=2: ω_a = 0.0 → spin locked to momentum (the #60 Thomas/Wigner result).

## T4: Magnetic moment = Bohr magneton

- μ_B = eℏ/2m = 9.2740e-24 J/T
- μ = g·μ_B·s = 9.2740e-24 J/T (g=2, s=½) = μ_B

## T5: Schwinger anomaly

- a = (g−2)/2 = α/2π = 0.00116141 (one loop)
- g (one loop) = 2.00232282
- measured a_e = 0.00115965 (rel diff 0.15%)
- tree g=2 geometric: True; α/2π needs the loop: True

## T6: Falsification criterion

- classical/scalar g = 1; spinor (BAM) g = 2; distinct: True
- σ·B Pauli term nonzero: True
- **BAM passes the falsifier: True**

## T7: B4 accounting

- g = 2, a = 0.0011614 (dimensionless)
- μ_B = eℏ/2m carries the single anchor (m): True

## T8: Assessment

- g factor: 2
- origin: Pauli/SU(2) σ·B term (T=iσ_y) + Hopf monopole (½)
- Thomas link: g=2 ⟺ ω_a=0 (spin tracks momentum, #60)
- magnetic moment: μ = μ_B (g·s = 1)
- anomaly: a = α/2π (one loop; tree g=2 geometric)
- remaining: α/2π from the throat loop; higher-order a_e; throat spinor from S_BAM

## Verdict

**G_FACTOR_DERIVED.** g = 2 DERIVED. The electron throat's gyromagnetic ratio g=2 follows from the BAM geometry — extending the spin thread (PR #60, the Wigner rotation) from kinematics to the magnetic moment.

PAULI/SU(2) ORIGIN. Minimally coupling the throat spinor (D=p−eA) and squaring the Dirac/Weyl operator gives, via the Pauli identity (σ·a)(σ·b)=(a·b)I+iσ·(a×b) and [D_i,D_j]=−ieε_ijk B_k, (σ·D)²=D²−eσ·B, so H=(p−eA)²/2m−(e/2m)σ·B. The magnetic-moment term carries the full σ (=2S), so μ=(e/2m)σ=g(e/2m)S with g=2. The factor 2 is the SU(2) anticommutator {σ_i,σ_j}=2δ_ij — the throat's non-orientable transport T=iσ_y=ε (B2). A classical (scalar) moment would give g=1; the spinor structure forces g=2.

HOPF MONOPOLE. The Hopf connection A_φ=½ cos χ is the spin-½ monopole (charge ½); minimal coupling gives μ=(e/m)S, g=2, and g·s=2·½=1 → magnetic moment = 1 Bohr magneton (the electron).

THOMAS/WIGNER LINK. The BMT anomalous precession is ω_a=(g/2−1)(eB/m); g=2 ⟹ ω_a=0, the spin stays locked to the momentum — exactly the Thomas/Wigner result of #60 (the kinematic Thomas precession conspires with Larmor so that at g=2 spin and momentum rotate together).

SCHWINGER ANOMALY. The leading anomalous moment is a=(g−2)/2=α/2π≈0.0011614 (vs measured a_e=0.00115965). This is a ONE-LOOP QED vertex correction: BAM's TREE geometry gives g=2 exactly; the α/2π requires the throat vertex/self-energy loop (beyond the tree-level geometric structure). B4: g and a are dimensionless; the moment scale μ_B=eℏ/2m carries the single anchor (m). g=2 is geometric/topological, independent of the anchor's value. Remaining: α/2π from the throat loop, the higher-order a_e series, and the explicit throat spinor from S_BAM.

## What this leaves open

- **α/2π from BAM.** The one-loop anomaly requires the explicit throat vertex/self-energy loop; the tree probe gives g=2 only.
- **Higher-order anomaly.** The full a_e series (α², α³, …).
- **The throat spinor from S_BAM.** The explicit minimally-coupled boosted throat spinor (shared with #59/#60).
