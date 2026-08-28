[![DOI](https://zenodo.org/badge/1181274003.svg)](https://doi.org/10.5281/zenodo.20225786)
# Geometrodynamics

**A research framework implementing and testing Wheeler's geometrodynamic program.**

This package computationally explores the hypothesis that structures
physicists call electromagnetism, charge, spin, confinement, **black
holes**, and **Bell correlations** may emerge from the geometry of
spacetime itself — specifically the Hopf fibration on S³, 5D Tangherlini
wormholes, topological flux-tube networks, coherent wormhole-throat
condensates, and non-orientable throat topology.

## Direction of the program: GR → QFT, *not* quantum gravity

**BAM derives quantum field theory *from* continuous (classical) general
relativity — it is the opposite of a quantum-gravity program, and does not
attempt to quantise gravity.** The foundation is a *classical*, continuous GR
geometry: the wormhole throat, the 5D Tangherlini bulk, the metric `f(r)`. The
quantum field theory — the matter spectrum, the propagator/exchange kernel, the
self-energy, the interaction vertices (PRs #116, #129–#140) — is the *derived*
output, obtained as standard field theory **on that fixed classical background**.

So the arrow runs **geometry → fields**, never **fields → geometry**:

  - the metric is a classical input, never a quantised dynamical field;
  - "throat", "horizon", "5D Tangherlini" name a *classical GR background*, and
    the probes that build propagators, vertices, and self-energies are deriving
    *QFT on that background*, in the precise sense of QFT-on-curved-spacetime;
  - asking BAM to "tackle quantum gravity" is therefore a **category error** —
    it would invert the program. Gravity here is the foundational *classical*
    layer from which quantum matter is reconstructed, not a thing to be
    quantised.

When the probes below speak of the path-integral measure `S_BAM`, the
one-loop determinant, or the bounded interacting vacuum, these are statements
about the **matter QFT** read off the classical throat geometry — not about a
quantum theory of the metric.

## Where ℏ enters: scale-free closure ledger + one geometric anchor

The closure-ledger arc (`experiments/closure_ledger/`, PRs #11–#74)
reduces most dimensionless parameters in the locked lepton surrogate
to closure-quantum invariants (`action_base = 2π`, `transport = 8π`,
`resistance = 7π/100`, `β_lepton = k_5²·(2π) = 50π`,
`ε = 7π/(100·k_5⁴)`). **The `pinhole γ = Σ V_max[1..5]` entry is reopened by
PR #271** — correcting the radial scalar operator removed the numerical support
for deriving `γ = 22.5` from the canonical geometry, and the chain's sensitivity
to `γ` (`d ln m_μ/d ln γ = −17.5`) turns the corrected −0.75 % residual into a
**measured 15.2 %** muon error. `γ` remains *required*; it is no longer
*derived*, and an audit (`maslov_dimensional_bridge_probe`,
PR #52) then established that the machinery is **scale-free**:
rescaling `R_MID → λ·R_MID` leaves every dimensionless output
invariant. By dimensional analysis, **exactly one external dimensionful
anchor is mathematically required** (B4 is irreducible). The Compton
bridge then collapses to

```
ℏ  =  m_e · R_MID · c              (equivalently  m_e = f_closure · ℏ / (ΔR·c))
```

That anchor need not be a particle mass: it is **relocatable to the
invariant bulk separation** `ΔR = R_OUTER − R_INNER` (PR #53,
`delta_r_scale_modulus_probe`), a cosmologically fixed length (the
throat is a static bound vacuole, decoupled from Hubble flow), with
`f_closure = ΔR/R_MID = 0.52`. The chain
**imposed `R_MID` → invariant geometric length `ΔR` → finite-self-energy
equilibrium** has each step more physical (PRs #55–#58):
`self_consistent_throat_radius_probe` recasts `R*` as a stable
equilibrium `E(R) = A/R + B·R²` of EM repulsion vs cohesion (`U_EM/(mc²) = α/2`,
no UV divergence); `cohesive_tension_derivation_probe` derives
`B = 4πσ` as the throat brane tension (the unique `R²` power by
power-counting); `brane_tension_tuning_probe` sharpens the bulk-gravity
relation to the **exact** RS fine-tuning `λ_crit = √(6|Λ₅|)/κ₅²`
(dimensionless factor √6, the flat / static-throat condition); and
`pair_production_threshold_probe` makes `2 m_e c²` the lowest stable
configuration (one Hopf charge per throat → C-conjugate
throat–antithroat pair).

**Scaffold status:** four of five mismatch terms (B1, B2, B3, B5) closed;
B4 audited as irreducible-by-dimensional-necessity. Full ledger:
`docs/bam_scaffold_status.md`. Release note:
`docs/scaffold_closure_release_note.md`.

**Reproduce in seconds:**

```bash
python -m experiments.closure_ledger.maslov_dimensional_bridge_probe
# Verdict: B4_IRREDUCIBLE — scale-free invariance verified.
```

## Why progress is possible beyond Wheeler's geometrodynamics

Wheeler's original geometrodynamic programme had the right *instinct*
— that what we call "matter" should ultimately be a property of
spacetime itself — but it stalled in the 1960s and 70s for a concrete
reason: it lacked the **global / topological machinery** needed to
turn that instinct into a quantitative spectrum.  The continuum
Einstein equations alone do not pick out discrete spectra; they
admit far too many solutions.  Wheeler's "charge without charge" and
"mass without mass" remained slogans precisely because there was no
mechanism to make them *count* anything.

The line continued here is concrete: discreteness arises from three
independent topological/geometric channels, all of which can be
written down explicitly and integrated numerically.

1. **Antipodal S³ closure.**  Compactifying the spatial slice as
   S³ replaces the open continuum with a closed cavity, so any
   field that closes on itself does so over a great circle of fixed
   length 2π.  Resonance on a closed cavity is intrinsically
   discrete; the closure constants (`action_base = 2π`, the
   integer-winding lock `4β = 100·(2π)` for the τ lepton) are
   *exact* topological invariants of this antipodal closure.  The
   closure constants are not fitted; they are read off from the
   global structure.
2. **Non-orientable throat/shell spectra.**  A wormhole throat
   that is non-orientable carries a Z₂ partition class (`p = ±`)
   which is a real topological label, not a continuous parameter.
   The unique orientation-reversing isometry of S³ that preserves
   the Hopf bundle is `T = iσ_y` (derived in `embedding/transport.py`
   without ansatz).  T² = −I is the 4π periodicity of spinors; the
   partition splitting drives every mass-ordering inversion in the
   shelled sector (the m_u < m_d but m_c > m_s pattern).  The
   throat orientation is what makes spin-½ unavoidable rather
   than imposed.
3. **Uniform bulk distance from outer to inner.**  The throat
   confines a radial coordinate to the finite shell `[R_INNER,
   R_OUTER]` (geometric units; throat at `R_MID = 1`).  In tortoise
   coordinates this becomes a finite interval with regular
   boundary conditions, which produces a discrete eigenmode
   spectrum (`tangherlini.radial.solve_radial_modes`) — bound
   modes `u_{l,n}(r*)` with frequencies `ω(l,n)`.  This is the
   bulk geometry's own quantization channel, independent of the
   S³ closure but composing with it.

What was missing in Wheeler's day — and what this package now
demonstrates operationally — is that these three channels **compose**.
The lepton ladder is a "minimal closure" spectrum where channel 1
(S³ closure) dominates: each lepton mass scales with its global
pass-count winding `β·k²` on a nearly bare closure skeleton, locked
by `4β_lepton = 100·(2π)`.  The quark ladder (added in this work)
is a "shell-coupled closure" spectrum where channel 1 picks up the
heaviest shell only and channels 2 and 3 — partition asymmetry on
the throat and bulk-mode coupling — determine the lighter shells.
Three of the four quark-sector residuals derive from
`tangherlini.radial.solve_radial_modes` and
`tangherlini.alpha_q.derive_alpha_q` to within 1%, on the same
tortoise grid that defines the radial bound modes (see
`docs/quark_axioms.md` §8 for the full derivation log and the
quantitative match per residual).

This is what allows progress: the right machinery for *quantitative*
geometrodynamics exists, it is just not the differential-geometric
machinery Wheeler had at hand.  Antipodal closure on a compact 3-space,
non-orientable throat topology, and bulk-mode confinement are each
old and individually well understood; what is new here is putting
them together and showing that they reproduce charged-lepton masses
to sub-percent and the six-quark mass ladder to ~1.6%.

## What the Code Validates

| Claim | Status | Evidence |
|-------|--------|----------|
| Charge quantisation from topology | **Verified** | c₁ = 1 to < 1e-9 error |
| Spin-½ from Hopf holonomy | **Demonstrated** | SU(2) sign-flip at 2π, illustrative |
| Coulomb law from throat eigenmode | **Verified** | BVP matches Q/r to rel_err < 3e-9 |
| Two-throat Coulomb force on S³ (finite separation) | **Demonstrated** | S³ Green response → V ∝ 1/r, F ∝ 1/r² (flat limit); F ∝ 1/sin²ψ on S³; Gauss law exact (`two_throat_coulomb_probe`) |
| α_q coupling ratios (no free parameters) | **Computed** | Forced-origin slope extraction |
| Möbius half-integer spectrum | **Verified** | Numerical vs analytic < 5% |
| Meson energy conservation | **Verified** | Drift < 1% over test window |
| Bridge nucleation / string breaking | **Verified** | Correct daughter topology |
| Hayward metric from throat density | **Derived** | n(r) → ρ(r) → m(r) → f(r) matches Hayward to < 1% |
| de Sitter EOS from Einstein eqs | **Derived** | p_r/ρ = −1 exact at all radii |
| SEC violation for regularity | **Derived** | Penrose-required SEC violation confirmed (~85% of interior) |
| Singularity avoidance (Hayward core) | **Derived** | K(0) = 24/l⁴ finite; metric now derived from throat density |
| Geodesic completeness | **Modeled** | Hayward infaller decelerates; heuristic completeness criterion |
| BH entropy from throat counting | **Consistent** | S_throat matches S_BH by construction (N set from area law) |
| Charge without charge (BH) | **Modeled** | Net Q from orientation sum, Q/N → 0 for large M |
| First law dM = T dS | **Checked** | Residual < 5%, Schwarzschild limit only |
| T from collective modes | **Derived** | T_mode matches T_surface_gravity to < 1% for M ≥ 3 |
| Core scale l ≈ Planck | **Derived** | l = 2M/√N ≈ 0.47 l_P, independent of mass |
| Schwarzschild recovery | **Verified** | Hayward → Schwarzschild as l → 0 |
| Two-horizon structure | **Verified** | Inner + outer horizons for 0 < l < l_crit |
| Singlet from throat transport | **Constructed** | T=iσ_y → |Ψ⟩ built from transport; E(a,b) = −cos(a−b) |
| T = iσ_y from Hopf fibration | **Derived** | Unique orientation-reversing Hopf-preserving map; 7 properties verified |
| Bell phases from Hopf holonomy | **Derived** | π/2 baseline + π[cos(θ_a)−cos(θ_b)]/2 from connection A = ½cos(χ)dφ |
| History closure → E = −cos(a−b) | **Derived** | SU(2) amplitudes × closure weights reproduce singlet; CHSH = 2√2 |
| History no-signaling | **Derived** | Marginals = ½ from branch enumeration; independent of remote setting |
| History conservation | **Verified** | Charge balance exact for Bell and transaction histories |
| Bulk identity Bell (kinematic) | **Verified** | Same E(a,b) from pure topology, no time stepping; both paths match |
| CHSH S = 2√2 (topological) | **Verified** | Exact Tsirelson value; topology determines correlations, cavity determines dynamics |
| No-signaling | **Verified** | Marginals = ½ from singlet; cavity dynamics don't alter spin correlations |
| Cavity detector-conditioned dynamics | **Dynamical** | Derived Hopf phases drive cavity ODE; packets fire on 0/π branches |
| Cavity persistent memory | **Verified** | Energy persists between steps; slow ring-down |
| Green kernel derivative | **Fixed** | Now matches analytic dG/dψ to < 10⁻⁴ |
| Lepton mass ladder (e, μ, τ) | **Closed** | Sub-percent all three generations from locked S³ axioms (see below) |
| S³ action base `action_base = 2π` | **Locked** | Hard topological invariant; default in all lepton scans |
| k=5 uplift `4β = 200π` (100 × 2π) | **Locked** | τ uplift equals exactly 100 S³ winding quanta |
| Closure cycle integer-quantised in 2π | **Verified** | `(N_e, N_μ, N_τ) = (3, 6, 109)` from antipodal + Hopf-throat + radial BS + τ-uplift |
| R_OUTER selected by cross-species fixed point | **REOPENED (PR #271)** | The bisection agreement was computed with the pre-#271 radial operator. The lepton observables do **not** currently select a unique `R_OUTER`: the locked Hamiltonian discards `r_outer` and sees only `γ`, so the two corrected `γ = 22.5` roots (`1.24614`, `1.26788`) give **bit-identical** masses. See `docs/scalar_operator_audit.md` |
| Pinhole γ ≈ Σ V_max[1..5] on Chebyshev grid | **REOPENED (PR #271)** | Under the corrected scalar operator the residual improves to −0.75 %, but `d ln m_μ/d ln γ = −17.5` at the lock makes that a **measured 15.2 %** muon error (linearised 14.0 %), and no channel set at `R = 1.26` lands near 22.5. `γ = 22.5` is required by the locked surrogate; its **geometric derivation is reopened**. The `9 %` first published here was a generic half-percent illustration misapplied to the actual residual — corrected in `docs/quark_residual_reaudit.md`. See `docs/scalar_operator_audit.md` |
| Transport = 8π = 4·(2π) | **Verified** | +0.13 % off the locked transport = 25.1; 4th closure quantum |
| Resistance = 7π / 100 | **Verified; selector superseded (PR #272)** | +0.94 % off the locked resistance = 0.2179. The stated selector — `R_OUTER` bisection — is reopened by PR #271, and under the legacy operator it had chosen the *worse*-fitting candidate (`4·(ω−1)` = +0.48 %). Under the corrected operator the competitor degrades to **+2.50 %**, so `7π/100` now wins on proximity as well: **conclusion survives, stated reason does not**. See `docs/quark_residual_reaudit.md` |
| Inner cutoff `ε = resistance / k_5⁴` | **Verified** | Closes the Compton bridge `ℏ = m_e R_MID c` to 0.04 % |
| Closure-quantum ledger closes modulo m_e | **PARTIALLY REOPENED (PR #271)** | The `2π`-quantised entries stand. The `pinhole γ` entry does not: its geometric derivation is reopened, so the ledger no longer reduces **every** dimensionless parameter to closure quanta. See `docs/scalar_operator_audit.md` |
| `D = 5` scalar quasinormal frequency (`ℓ = 1`) | **Settled against published values (PR #274)** | `ω = 1.01601691 − 0.36232802i` from Matyjasek 2021 (continued fractions + Hill determinants, arXiv:2107.04815). An independent characteristic evolution here reproduces it to `0.018%` in damping. This **confirms PR #270's Kerr–Schild code (`0.005%`) and excludes its tortoise damping (`27.1%` off)** — #270's own prime suspect was the wrong code. No autopsy: neither #270 code was landed. See `docs/ringdown_cross_validation.md` |
| Self-convergence as an error bar | **REFUTED BY MEASUREMENT (PR #274)** | This round's own step-size study gave a last successive difference of `4.0e-5` while the true error against the published value is `1.1e-4` — **2.7× larger**, with `h = 0.05` closer than `h = 0.025`. A convergence study measures only the error it refines away. Nothing internal to a solver can expose this |
| Frequency-domain QNM shooting in real `r` | **Cannot settle it (PR #274)** | Reproduced #270's non-convergence rather than fixing it: the root moves with every knob because for `Im ω < 0` the outgoing piece grows like `e^{\|Im ω\|R}` and swamps the coefficient being zeroed. Sixth-order WKB by finite differences also **diverges** under refinement (`9.0 → 18.6 → 623`) |
| Quark mass ladder (u, d, s, c, b, t) | **Fitted** | 1.6% max rel err on s, c, b, t with d-anchor, four shell-index axioms, and one phenomenological β |
| Quark CKM / flavor-CP realization (v4) | **Locally flexible realization, NOT a prediction (PR #273)** | `rank K = 4` where `K = J_F ker(J_M)`: the mass-preserving parameter freedom spans the *entire* physical flavor space, so the CKM agreement is not evidence for the Hamiltonian. Zero left-null vectors ⟹ no predicted relation. Holding the derived `φ_h = π/k₅` fixed leaves rank 4. See `docs/flavor_identifiability.md` |
| Quark v4 counting: "+3 parameters, +5 independent observables, net +2" | **REFUTED (PR #273)** | A unitary 3×3 CKM has exactly **four** physical parameters, so "+5 independent" exceeds the ceiling. Measured calibration dimension of the v4 additions (`φ_h` fixed) = **4**; net surplus **≤ 0**. The six diagonal-shift symbols measure flavor rank 2+2, not 3+3 — the trace of each triple is an exact CKM gauge |
| Quark shell-index axioms (ε, η, χ, phase) | **Geometric** | All four expressible in `k_5 = 5` only: `(1−1/k_5², k_5, (k_5−1)·k_5, 0)` |
| Quark residual sector (transport, pinhole, resistance) | **Derived — and IMPROVED by PR #271's correction** | All three move *toward* their locked values under the corrected operator, opposite to the lepton sector. The difference is elasticity: `d ln m_s/d ln pinhole = +4.8` against the lepton chain's `−17.5`. See `docs/quark_residual_reaudit.md` |
| Pinhole = `Σ V_max(l=1..5)` (tortoise grid) | **Verified (re-derived)** | `−1.09%` (legacy) → **`+0.36%`** off the fitted 22.25; the residual changes sign |
| Transport = `mean ⟨u_l\|V_{l+2}−V_l\|u_{l+2}⟩` | **Verified (re-derived)** | `+0.88%` (legacy) → **`+0.70%`** off the fitted 0.54; the cross-ℓ operator itself is exactly invariant, only the eigenvectors drift |
| Resistance = `transport · ln(α_q(k_5)/α_q(k_1))` | **Verified (re-derived)** | `+0.49%` (legacy) → **`−0.02%`** off the fitted 0.14. The previously published `−0.43%` used the *locked* transport 0.54 in place of the derived one |
| The derived residuals **reproduce** the quark ladder | **NOT ESTABLISHED (PR #272)** | Never established under *either* operator: substituted into the locked Hamiltonian with nothing retuned, the derived triples give `3.44%` (legacy) and `3.78%` (corrected) against the fitted lock's `1.61%`. Per-knob agreement is not ladder agreement; the legacy triple's errors merely cancel better. See `docs/quark_residual_reaudit.md` |
| `R_OUTER = 1.26` bracketed by the lepton and quark sectors | **New, and weak (PR #272)** | Legacy puts 1.26 *outside* the bracket the two sectors define (both demand more: 1.2723, 1.2874); corrected **straddles** it (1.2564, 1.2679), with 1.26 at 31% across a 0.91% window. Suggestive — a 0.91% bracket admits anything inside it — and **not** the single-sector fixed point PR #271 reopened |
| Quark winding β = N·π/2 with N=466 | **Phenomenological (scope sharpened)** | `N = 2·n_part`, parity (Z₂) topological; `n_part = 233` is fit compensator absorbing the inter-generation hierarchy — diagnosed as **dynamical** (irregular, neither power-law like leptons nor exponential like neutrinos), and specifically the **flavor puzzle**: quark mass ratios are RG-invariant ⟹ not `α_s` running but the (irregular) Yukawa couplings, which overflow the geometric shell capacity (`quark_beta_*` probes, PRs #76, #80, #97, #98) |
| Compton antipodal kinematics | **Verified** | Closure-compatible: front + back-mouth 4-momentum conservation under (E, **p**) → (E, −**p**); inter-mouth γ skew vanishes identically; throat-pinch skew is recoil-induced `O(ω²/m²)` |
| Compton S³-propagator pole `1/(s−m²)` | **Verified** | S³ Green function `G(ψ) ∼ 1/ψ` with `ψ ∝ s−m²` reproduces QED propagator pole; fitted exponent 1.0002 across five ω-decades |
| Thomson `(1+cos²θ)` angular factor | **Derived** | Polarization-summed BAM amplitude reproduces Klein-Nishina at ω → 0 from transverse photon polarisations on the tangent bundle |
| Compton vertex coupling `γ = −3/2` at O(ω/m) | **Derived** | Exact analytic solution to the 4-equation linear system in {1, c, c², c³} basis; clean rational coefficient |
| `γ = −3/2` is d-independent | **Verified** | Numerical γ(d) = −3/2 in d ∈ {3, 4, 5, 6, 8} to 7-digit precision; falsifies the embedding-dim/polarization-count origin |
| Compton vertex closed-form resummation | **Derived** | `F²(x, c) = 4·x³·(x²+1−x·sin²θ) / [(1+c²)·(1+x)²]` with `x = ω'/ω` reproduces Klein-Nishina to all orders in ε up to ε ~ 2 (machine precision); the perturbative PRs #31–34 are Taylor expansions of this closed form |
| F² and masses from one master integral | **Derived** | Single `C × S³` master functional `ℳ = G_C ⊗ 𝒢_{S³}`: ω-poles → mass spectrum, throat boundary → `K(x)`, S³ Hopf → `Q(x,c)`; vertex residue = `F²=K²·Q` to `2e-14`. Closes scaffold barrier B5′ (`master_integral_probe`, `docs/bam_scaffold_status.md`) |
| Dimensional anchor (B4) is structural, not a gap | **Audited** | Closure-ledger/Maslov machinery is scale-free (rescale `R_MID → λ·R_MID` → all dimensionless outputs invariant), so exactly one external dimensionful anchor is required; relocatable to the cosmologically-invariant bulk separation `ΔR`, giving `m_e = 0.52·ℏ/(ΔR·c)` (`maslov_dimensional_bridge_probe`, `delta_r_scale_modulus_probe`) |
| Finite-self-energy throat equilibrium | **Derived / Modeled** | `R* = (A/2B)^{1/3}` stable minimum of `E(R) = A/R + B·R²`; throat caps the EM field so `U_EM/(mc²) = α/2` (finite, no UV divergence) (`self_consistent_throat_radius_probe`, PR #55) |
| Cohesive brane tension `B·R²` | **Derived** | `E = σ·Area = 4πσR²` (`B = 4πσ`); `R²` power uniquely selected by power-counting (Tangherlini junction is `R¹`, EH is `R¹`, bag is `R³`) (`cohesive_tension_derivation_probe`, PR #56) |
| RS-like √6 brane tuning | **Derived** | `λ_crit = √(6\|Λ₅\|)/κ₅² = 6k/κ₅²` from `Z₂` Israel junction `K_μν = −κ₅²λ/6 h_μν` + bulk `AdS₅` (`Λ₅ = −6k²`); flat / static-throat condition `Λ₄ = 0` (`brane_tension_tuning_probe`, PR #57) |
| Pair-production threshold `2 m_e c²` | **Derived** | One Hopf charge per throat (`\|c₁\| = 1`) ⟹ `Σc₁ = 0` forces C-conjugate throat–antithroat pair; bubble-nucleation barrier `R_c = 2σ/ρ`; Schwinger critical field `eE_S R_MID = m_e c²` (`pair_production_threshold_probe`, PR #58) |
| Moving throat = relativistic particle | **Verified** | Dispersion `ω(k)=√(ω₀²+c²k²)` ⟹ `E²−(pc)²=(mc²)²` with `mc²` = static eigenvalue `ω(1,0)` to machine precision; closed `S³` breaks global Lorentz, suppressed by `(R_MID/R_cosmo)² ~ 10⁻⁷⁸` (`stable_moving_throat_probe`, PR #59) |
| Spin-½ Wigner rotation (relativistic) | **Verified** | Hopf-holonomy `∮A = π cos χ` reproduces Wigner `SU(2)` rotation from two non-collinear boosts (`SL(2,C)` composition); the same `½` factor / spinor double cover / `½ × solid angle` (`spin_wigner_rotation_probe`, PR #60) |
| Throat `g = 2` | **Derived** | Pauli/SU(2) `T = iσ_y` + Hopf monopole `A_φ = ½ cos χ`; `(σ·D)² = D² − eσ·B` with `σ = 2S` (the `SU(2)` anticommutator factor of 2); BMT anomalous precession vanishes ⟺ `g = 2` (`gyromagnetic_ratio_probe`, PR #61) |
| Schwinger anomaly `a = α/2π` | **Reconstructed** | One-loop dressing: virtual photon = `S³` Green-function exchange (flat `1/q²`), vertex = throat pinch, Feynman-parameter `∫₀¹ 2z dz = 1` ⟹ `F₂(0) = α/2π = 0.0011614`; vs `a_e = 0.00115965` to ~0.15% (`throat_vertex_loop_probe`, PR #62) |
| `S_BAM` loop measure `1/(2π)` | **Structurally identified** | The `1/(2π)` in `a = α/(2π)` = BAM closure-quantum loop measure factor — same `2π` as `action_base`, `Φ_avail(k) = 2π(k+1)+…`, `β_lepton = k_5²·(2π)`, Hopf, throat dwell, `ε`'s `4β/(2π) = 100`; closed cycle of length `2π` → measure `dk/(2π)` per loop dim. Full covariant `(2π)^d` path-integral derivation open (`s_bam_loop_measure_probe`, PR #74) |
| The hard `S_BAM` path-integral measure: full loop-measure construction | **Structurally defined; analytic core open** | Takes up PR #74's flagged open work — builds the full measure `Z = Σ_sectors ∫ Dμ[X] e^{−S_BAM[X]}` around the `1/(2π)` factor. **Arena:** loop space `LS³ / (Diff S¹ ⋉ U(1)_Hopf ⋉ Z₂)`, `Dμ ~ Π dk/(2π)`. **Fixed (computable):** closure quantum `2π` = loop holonomy; superselection sectors = the closure ledger (homotopy `k`, `c₁∈π₃(S²)=ℤ`, `n_part`); **odd-k lemma UPGRADED to the Z₂ orientation-anomaly condition** `e^{ikπ}=−1 ⟹ k odd` (even `k` = torus cover only); the PRs #87–#90 bounces = the leading saddle. **Hard part:** `Diff(S¹)` gauge-fixing ⟹ FP(`bc`-ghost) × fluctuation-det; the fluctuation operator (= 2nd variation of `S_BAM` = Tangherlini cavity operator) is stable (min `ω²≈1.11>0`). **Open:** the bare determinant `Π ω_n` diverges (log-det → ∞) ⟹ needs zeta/heat-kernel regularization; `Z` not yet rigorously constructed. Prior saddle results (leading `e^{−S}`) unaffected (`s_bam_path_integral_measure_probe`, PR #115) |
| Regularize the Tangherlini fluctuation determinant | **Analytic core CLOSED — finite, two ways** | Resolves PR #115's open one-loop factor: the divergent bare `Π_n ω_n` is regularized to a **finite, scheme-independent** value by two independent standard methods that agree. **Gel'fand–Yaglom** (no mode sum — one IVP solve): `det(H)/det(H_free) = y(L)/L = 1.57437` (log `0.45386`), converged to 6 digits `N = 2000 → 32000`, zero interior nodes (no negative modes). **Zeta/heat-kernel:** `ζ(0) = a₀ = −1/2` (the universal Dirichlet-interval value — finite, no zero mode, no anomaly), Weyl leading coeff `a_{−1/2} ≈ L/√(4π)` to 0.9%, counting `N(λ) ≈ (L/π)√λ` confirmed. The `S_BAM` one-loop measure factor is finite and computable. **Still open:** closed-form expression (the value is numerical) and the absolute `Z` normalization (the `κ₅²/Λ₅` anchor, PR #112) (`tangherlini_fluctuation_determinant_probe`, PR #116) |
| Diff(S¹) Faddeev–Popov / ghost determinant | **Gauge sector complete — finite, anomaly-free** | Supplies the measure's gauge sector (PR #115 flagged it; PR #116 did matter). Worldline reparametrization: gauge-fixing the loop einbein leaves **1 Teichmüller modulus** `L` (circumference = Schwinger proper time) + **1 CKV** (rigid `U(1)` rotation). FP operator `P = d/dτ` (vector ghost ↦ einbein variation), `P†P = −d²/dτ²`, kernel = constants = the 1 CKV. **The FP ghost determinant is the `bc`-ghost integral `Δ_FP = det'(P) = det'(P†P)^{1/2} = L`** — the **square root** of the intermediate `det'(P†P) = L²` (`ζ(0) = −1`; both verified to machine precision). **Corrected measure** `Z = Σ ∫ (dL/L)·det^{−1/2}_matter·e^{−S}`: `Δ_FP = L` is the einbein→proper-length Jacobian (⟹ modulus measure `dL`), ghost L-power **`L¹`** (not the `L²` of the first draft); the `1/L` is the CKV factor. **PR #74 unchanged:** `1/L = 1/(2π)` at the closure loop `L = 2π` is the CKV (c-ghost zero-mode) factor, independent of the determinant power. **Anomaly-free:** 1D worldline has no conformal anomaly (vs 2D string `c = −26`); the only nontrivial anomaly is the discrete `Z₂` (odd-k, PR #115). Open: abs `Z` (`κ₅²/Λ₅`), multi-loop (`diff_s1_ghost_determinant_probe`, PR #117) |
| First-order Diff(S¹) FP ghost audit | **L-power fixed: ghost is L¹ (first order)** | Rigorous audit distinguishing the 4 objects: `P = ∂_τ` (first order, eigenvalues `2πin/L`, 1 zero mode = CKV), `P†P = −∂_τ²` (second order), `det'(P)`, `det'(P†P)`. **`det'(P†P) = L²`; `det'(P) = det'(P†P)^{1/2} = L`** (verified). **η-invariant:** `η(−i∂_τ) = 0` (spectrum symmetric `n↔−n`) ⟹ `det'(∂_τ) = +L`, no anomalous phase (antiperiodic/Möbius sector: `η = 0` too but **no zero mode ⟹ no CKV**). **Convention:** the physical FP is the first-order `bc` system, `Δ_FP = det'(P) = L`; `det'(P†P) = L²` arises **only** under an explicit second-order ghost convention (over-counts by one `L`). **No double-counting (proof):** the ghost space splits `ker(P)`[CKV] ⊕ `ker(P†)`[modulus] ⊕ nonzero; `det'(P)` is the **primed** det over **nonzero modes only** (SVD: exactly **1 zero singular value**, right-null = CKV), so the CKV norm enters **only** `Vol(CKG)` and the modulus norm **only** `dL` — each divided **once**. (The first draft's extra `√L·√L` division alongside `1/Vol(CKG)` double-counted the CKV; removed.) **Measure table:** `Z = Σ ∫ (dL/L)·det^{−1/2}_matter·e^{−S}`, single `1/L = 1/Vol(CKG)` (= PR #74's `1/(2π)` at `L=2π`); `det'(P)=L` folds into the matter heat kernel; net L-power `dL·L^{−1−d/2}`. Open: abs `Z`, multi-loop (`diff_s1_first_order_ghost_audit_probe`, PR #118) |
| Phase / η-invariant framework for `det′(∂_τ)` | **Phase = local ζ(0) + topological η; both BAM sectors η=0** | Builds the full framework PR #118 only asserted. `P=∂_τ` anti-self-adjoint (eigenvalues `2πin/L`), `A=−i∂_τ` self-adjoint; modulus `|det′(∂_τ)| = det′(P†P)^{1/2} = L` unambiguous. **Singer/APS phase formula:** `det′(A) = |det′|·exp[±i(π/2)(ζ_{|A|}(0) − η_A(0))]` — phase splits into a **local** (heat-kernel/scaling) `ζ(0)` piece and a **topological** (spectral-asymmetry) `η(0)` piece. **η with flux** (Hopf holonomy `a = kχ/2π`): `η_A(0) = 1 − 2a` (Hurwitz `ζ_H(0,a)=½−a`); reduced `η ≡ 0` for periodic (zero mode = CKV removed) and antiperiodic. **Concrete:** `det(∂_τ+m)_periodic = 2sinh(mL/2) → det′(∂_τ) = L` (residue); `det(∂_τ+m)_AP = 2cosh(mL/2) → det = 2`. **BAM:** orientable `a=0` and Möbius `a=1/2` both `η=0` ⟹ `det′(∂_τ)` real (rigorously justifies PR #118's `+L`); generic holonomy gives an η-phase `exp[−i(π/2)(1−2a)]` (open) (`detprime_dtau_eta_invariant_phase_probe`, PR #119) |
| High-resolution **lattice validation** (discrete ≡ continuum) | **Software reproduces the analytic derivation** | Verifies the discrete finite-difference operators reproduce the continuum analytic results of PRs #116–#119. **Eigenvalues** `−∂_τ²` → `(2πk/L)²`, relative error **`O(1/N²)`** (ratio exactly 16 per `N×4`). **Ghost det** (periodic): lattice log-det `Σ log[2−2cos(2πk/N)+(mh)²]` → continuum `(2sinh(mL/2))²`, `O(1/N²)`; transfer-matrix `2(cosh Nα−1)` cross-check at `N=10⁶`. **Antiperiodic** → `(2cosh(mL/2))²`; `m→0` ⟹ `det′(−∂_τ²)=L²`, `det_AP=2`. **Generic holonomy `a∈{1/4,1/3,2/3,3/4}`** (twisted BC `e^{2πia}`): twisted eigenvalues → `2π(n+a)/L` `O(1/N²)`; **`|det P_a|=2sin(πa)` EXACT on the lattice** (identity `Π 2(1−cos(2π(k+a)/N))=|1−e^{2πia}|²=4sin²(πa)` → `√2,√3,√3,√2`); `η(a)=1−2a`; **branch convention** `ζ(0)=0` ⟹ phase `(π/2)(1−2a)`, `det P_a=2sin(πa)e^{i(π/2)(1−2a)}` (= `1+i, 1.5+0.866i, …`; `a=1/2` ⟹ real `2`). **η = 0 EXACT at finite N** (centered `∂_τ`, odd `N`, 1 zero mode). **Tangherlini GY** `det(H)/det(H_free) → 1.574370` (PR #116). Structural/symmetry quantities (incl. `|det P_a|`) exact at finite `N`; finite-difference `O(1/N²)` (`lattice_validation_probe`, PR #120) |
| BAM **sector-phase ledger** (continuous η vs discrete Z₂) | **Factorizes; no double-counting** | Converts the validated `det'(∂_τ)` η-machinery into a ledger of the loop-measure phase. **Two independent structures:** U(1) holonomy `a` (connection, continuous) and orientation `w₁`/odd-k parity (discrete). **Continuous η-phase** `e^{i(π/2)(1−2a)}` from the holonomy — `θ(a)=(π/2)(1−2a) ∈ (−π/2,+π/2)` for `a∈(0,1)`, confined to the **open right half-circle** (`Re>0`), `=+1` at `a=1/2`, **never `−1`**. **Discrete Z₂ sign** `(−1)^k` from the Möbius/odd-k orientation (`+1` torus, `−1` Möbius). **No double-counting (proof):** (a) different groups U(1) vs Z₂; (b) different geometry connection vs orientation; (c) the η-phase never reaches `−1` (closest ≈ `√2`), so the Möbius `−1` is purely Z₂ — and at `a=1/2` the η-phase is `+1`, so the antiperiodic det's sign is entirely `(−1)^k`. Factorized: `det_full = |det P_a|·e^{i(π/2)(1−2a)}·(−1)^k`, each factor once (`bam_sector_phase_ledger_probe`, PR #121) |
| **Factorized sector sum Z** (full one-loop measure assembled) | **Discrete Z₂ × continuous η; graded UV cancels** | Assembles PRs #74,#115–#121 into `Z = Σ_{k odd, c₁, n_part} (−1)^k ∫(dL/L) det^{−1/2}_matter · e^{i(π/2)(1−2a)} · e^{−S_BAM}`. **Factorizes:** the Z₂ orientation sign `(−1)^k` is a sector-constant (winding parity, not `L`/`a`), so it pulls **out** of the continuous integral ⟹ `Z = Σ_{discrete} (−1)^k × [continuous moduli integral]` — discrete Z₂-signed (topological) sum ⊗ continuous η-phased (analytic) integral, no double-counting (PR #121). **Z₂-graded UV cancellation:** the Weyl term `a_{−1/2}=L/√(4π)` is BC-independent, so it cancels between orientable (+) and Möbius (−) sectors — each heat trace `θ ~ L/√(4πt) → ∞` but `θ_per − θ_anti ~ e^{−π²/t} → 0` (UV-finite). Every factor finite/validated; open: absolute normalization (`κ₅²/Λ₅`), non-perturbative convergence, multi-loop (`bam_factorized_sector_sum_probe`, PR #122) |
| **APS quark partition index** (from the factorized sum) | **Fixes the topological doubling, not n_part's value** | Reads the Witten/APS index off the factorized sector sum (PR #122). The Z₂ grading `(−1)^k` makes `I = Tr(−1)^k` topological; the **APS ξ-invariant** `ξ(a) = (η+h)/2 = 1/2 − a = ζ_H(0,a)` is the η-boundary term. **Integer index = spectral flow:** as `a:0→1` one mode crosses zero ⟹ `ξ(0⁺) − ξ(1⁻) = 1` (integer). **Applied to quarks:** `N_q = 2·n_part = 466` — the **even doubling** IS the Z₂-graded structure (the orientation index pairs/doubles the modes). **Topological vs residual:** the doubling `N_q = 2·n_part` (even across all 12 §8 ablations) + the integer index are §8-**stable** (the mod-2 / APS topological content); the bare value `n_part` (drifts 216–255) is the non-topological **residual** — formalising the PR #97/#107 compensator split. The index derives the structure, not the value (`aps_quark_partition_index_probe`, PR #123) |
| **APS lepton partition index** (the clean contrast) | **Fully determined — structure AND value, no residual** | The same APS audit (PR #123) on the **lepton** sector. `N_lepton = 4·k₅² = 100` with `k₅ = 5` the **derived** bulk dimension `dim(S³)+2` (PR #73), `β_lepton = k₅²·2π = 50π` (PR #71); 3 generations = `(k₅+1)/2` (odd-k `k∈{1,3,5}`). Same machinery: `ξ(a) = 1/2 − a`, spectral flow `= 1` (universal). **But the outcome flips:** because `k₅` is a fixed derived integer (not a compensator), `N_lepton = 4·k₅²` is fixed in **both structure AND value** — **no residual**. Contrast: quark `N_q = 2·n_part` fixes structure only (`n_part` drifts 216–255); lepton `N_lepton = 4·k₅²` fixes everything. **Leptons are the clean APS case; the quark `n_part` is the program's lone compensator residual** (`aps_lepton_partition_index_probe`, PR #124) |
| **Combined matter-sector APS ledger** (the capstone) | **Leptons derived; quarks one residual; budget assembled** | Combines #123/#124 and ties to the input budget (#104–#108, #112). Every matter partition is **(derived topological factor) × (feeding integer)**, with the topology (structural factor + integer spectral flow `=1`, `ξ(a)=1/2−a`) derived **everywhere**; only the feeding integer can be residual. Ledger: **lepton** `4·k₅²=100` (k₅ derived, **no residual**), **quark** `2·n_part=466` (`n_part` residual), **neutrino** `ε` (order-of-mag derived, value residual). ⟹ exactly **one matter-partition residual** (`n_part`). **Full input budget:** 1 dimensionful anchor `G` + 4 dimensionless residuals {`n_part`, `√σ/m_e≈830`, `ε`, `α`} + the universal flavor puzzle. APS isolates `n_part` as the unique matter-**partition** residual (the others are a ratio, a compliance, a coupling); it organizes the residuals, does not remove them (`combined_matter_sector_aps_ledger_probe`, PR #125) |
| **Non-perturbative convergence audit** of the Z₂-graded sector sum | **Converges — finite in all three pieces** | Audits the PR #122 open item: does `Z = Σ_{k odd, c₁, n_part} (−1)^k ∫(dL/L) det^{−1/2}_matter · e^{i(π/2)(1−2a)} · e^{−S_BAM}` converge non-perturbatively? It factorizes over three independent labels, each finite. **Winding sum FINITE:** the odd-k lemma + available phase `Φ_avail(k) = 2π(k+1) + 50π·max(0,k−3)²` cap `k ∈ {1,3,5}` (3 generations, `k₅=5` the bound) — a 3-term sum, not a tower; `k=7` costs 2563.5, far over budget. **Hopf-charge sum CONVERGENT:** `Σ_{c₁∈ℤ} e^{−A c₁²} = √(π/A)·θ₃ → √(π/A)` (verified `A=0.5,1,2`), Gaussian `c₁²` cost ⟹ absolutely convergent; `Σc₁=0` (PR #58) constrains further. **Moduli integral FINITE at both ends:** `∫(dt/t)[θ_per−θ_anti]e^{−m²t}` — UV (`t→0`) killed by the Z₂ cancellation `θ_per−θ_anti ~ e^{−π²/t}→0` (the grading removes the Weyl divergence the individual BCs carry; integrand `~9·10⁻¹⁴` at `t=0.02`); IR (`t→∞`) killed by the mass gap `e^{−m²t}` (`0.61, 0.17, 0.0075` at `m=0.3,0.5,1.0`). ⟹ `(finite winding)×(convergent Hopf theta)×(finite moduli) ⟹ converges`. Open: absolute normalization (`κ₅²/Λ₅`), multi-loop measure (`z2_graded_sector_sum_convergence_probe`, PR #126) |
| **5D Tangherlini bulk lift** (the throat's parent geometry) | **Genuine D=5 vacuum; cavity curvature-regular; AdS₅ reconciled** | Lifts the PR #116 Tangherlini cavity operator (`V = f[l(l+2)/r² + 3rs²/r⁴]`, `f = 1 − (rs/r)²`) to its explicit 5D parent metric `ds² = −f dt² + f⁻¹dr² + r²dΩ₃²` and verifies the throat is the boundary trace of a real D=5 geometry (curvature computed by a self-contained numerical GR routine). **Ricci-flat vacuum:** `R_μν = 0`, `Λ = 0` (verified across the cavity) — asymptotically flat, distinct from the AdS₅ RS bulk. **Cavity curvature-regular:** Kretschmann `K = 72 rs⁴/r⁸` (numeric ≈ analytic to 1e-3), finite on the whole cavity (72 at the throat → 11.3 at `R_OUTER`); the only true singularity is at `r=0`, behind the throat (`r=rs` is a coordinate/horizon singularity). **Throat = S³ horizon** at `r=rs=R_MID` = BAM's Hopf base `S¹→S³→S²`. **Potential descends from D=5:** `l(l+2)` = S³ Casimir `l(l+D−3)` (`D−3=2`), `3rs²/r⁴ = (D−2)/(2r)·f'` (coeff `D−2=3`) ⟹ `k₅ = D_bulk = 5` (PR #73) realised as the genuine bulk dimension. **Hawking period carries 2π:** `κ = 1/rs`, `T_H = 1/(2π rs)`. **AdS₅/RS reconciliation:** the Schwarzschild–Tangherlini–AdS₅ metric `f = 1 − rs²/r² + k²r²` is Einstein (`R_μν = −4k²g_μν`, `Λ₅ = −6k²`, verified), interpolating the Tangherlini neck (`k²r²→0`, #116) to the AdS₅/RS asymptote (#57, √6); cavity correction `O(10⁻²)` for `k·rs ≲ 0.1`. Open: exact AdS scale `k = κ₅²/Λ₅` (#112), global brane-localised solution (`five_d_tangherlini_bulk_lift_probe`, PR #127) |
| **Horizon-regular coordinate lift** for the throat | **Coord. singularity removable; antipodal bifurcation = C-swap** | Builds the horizon-regular charts that remove the throat's *coordinate* singularity (flagged in PR #127), make the crossing smooth, and exhibit the antipodal structure. **Removable:** `g_rr = 1/f → ∞` at `r=rs` while `K = 72 rs⁴/r⁸` finite ⟹ coordinate artifact. **Eddington–Finkelstein regular:** `ds² = −f dv² + 2 dv dr + r²dΩ₃²`; at the throat `g_vv=0` but `g_vr=1` ⟹ `det g = −r⁶ sin⁴χ sin²θ` finite/nonzero, and `K = 72 rs⁴/r⁸` computed in EF coords (same regular geometry, nondegenerate metric — verified by the numerical GR routine). **Tortoise vs proper:** `r* → −∞` (infinite optical distance) but proper `∫dr/√f ≈ √(2 rs Δr)` finite = the ε healing length `√(2 rs ε)` (#112). **Surface gravity & Kruskal:** `κ = f'(rs)/2 = 1/rs` (so `κ·rs = 1`); the Kruskal factor `F = (r+rs)²/r²·e^{−2r/rs}` is finite/nonzero at the throat (`F(rs) = 4 e⁻²`) because `κ·rs = 1` cancels `f`'s simple zero; `T_H = 1/(2π rs)`. **Maximal extension:** `UV = −(1/κ²)e^{2κr*} → 0` at the throat — the bifurcate Killing horizon `U=V=0`; four regions (I exterior, II interior, III antipodal exterior, IV white hole). **Antipodal = C-swap:** the isometry `(U,V,Ω) → (−U,−V,Ω̄)` preserves `UV` (region I ↔ III) — the geometric home of BAM's throat ↔ antithroat identification (`C` inner/outer swap #63, `c₁→−c₁` #58); **"Bulk Antipodal Mechanics" is the antipodal identification of the throat's Kruskal horizon**. Open: nucleation rate (#58/#88), exact AdS scale `k` (#112), global brane solution (#127) (`five_d_tangherlini_throat_horizon_lift_probe`, PR #128) |
| **Null throat boundary conditions** for wave transport | **Antipodal l-parity BC; unitary mirror, not absorbing horizon** | Derives the BC the null throat (5D horizon) imposes on the transported matter waves (PR #116 cavity, PR #128 antipodal structure). **Vanishing potential:** `V_l = f[l(l+2)/r² + 3rs²/r⁴] ∝ f → 0` at the throat ⟹ near-horizon `−ψ''=ω²ψ`, pure null modes `ψ ~ e^{±iωr*}`. **Three candidate BCs:** ingoing/absorbing (`e^{−iωr*}`, flux-losing, non-unitary), reflective wall (Dirichlet/Neumann box, #116), antipodal (BAM-native, #128). **Antipodal map fixes the BC by l-parity:** S³ harmonics carry `Y_l(−x) = (−1)^l Y_l(x)` (degree-l harmonic polynomials; verified), so single-valuedness `Φ(U,V,Ω) = Φ(−U,−V,Ω̄)` forces radial parity `(−1)^l` across the throat — **even-l ⟹ Neumann `ψ'(throat)=0` (antinode), odd-l ⟹ Dirichlet `ψ(throat)=0` (node)** (twisted/Möbius field flips even↔odd, #67/#121). **Unitary mirror:** both antipodal BCs are real ⟹ throat flux `j ∝ Im(ψ*ψ') = 0` (verified) — a perfect mirror, no flux lost; vs the ingoing BC's `j = −ω ≠ 0` (absorbing sink). The antipodal throat conserves flux (global CPT/unitarity, #64): what falls in on one sheet re-emerges on the antipodal sheet. **Spectrum:** real, positive, discrete (unitary cavity); even-l (N) vs odd-l (D) families distinct (lowest `ω²`: l=0→1.37, l=1→5.27, l=2→2.03, l=3→6.73) — the wave-transport face of the even-k/odd-k Z₂ structure (#67/#121). Open: full QNM spectrum (complex ω), throat↔antithroat nucleation rate (#58/#88) (`null_throat_boundary_conditions_probe`, PR #129) |
| **Antipodal vs absorbing throat QNM spectrum** | **Antipodal → real undamped (stable matter); absorbing → complex ringdown** | Computes the full frequency spectrum of the BAM cavity `−d²/dr*² + V_l = ω²` on `[R_MID+ε, R_OUTER]` (shell wall at `R_OUTER`) under the two throat BCs — the spectral fingerprint distinguishing BAM's antipodal throat (#129) from an ordinary absorbing horizon. The absorbing case (ingoing `ψ'(throat)=−iωψ`) is a quadratic eigenvalue problem solved by companion linearisation. **Antipodal ⟹ real ω:** the real l-parity BC (Neumann even-l / Dirichlet odd-l) is self-adjoint ⟹ `Im(ω)=0` (verified `max|Im ω|≈0`) — undamped normal modes, quality factor `Q=∞`, sharp zero-width lines, l-parity graded. **Absorbing ⟹ complex ω:** the ingoing BC is non-self-adjoint ⟹ `ω = ω_R − i|ω_I|`, `Im(ω)<0` — damped quasinormal ringdown (fundamental `≈1.89−1.24i`), lifetime `τ=1/|ω_I|`, `Q=ω_R/(2|ω_I|)~O(1)` (thin cavity leaks fast into the throat). **Physical consequence:** a matter state is a sharp mass (stable particle) only if its mode frequency is real — the absorbing throat gives every state a width/complex mass (decaying resonance), so **stable matter (the lepton/quark bound states) requires the unitary antipodal throat** — the spectral face of global CPT/unitarity (#64). Open: idealised `r*→−∞` horizon QNMs, GW coupling, absolute normalisation (`antipodal_vs_absorbing_qnm_probe`, PR #130) |
| **Geometric throat arc synthesis** (capstone of #116, #127–#130) | **One primitive, five faces: the antipodal 5D-horizon identification** | Capstone re-verifying a keystone from each arc member together (a cross-arc consistency check) and consolidating the unified picture. **The keystones (re-run together):** `K=72` at the throat (regular, #116/#127), `T_H=1/(2π rs)=0.159`, EF `det g=−0.299` (nondegenerate, #128), Kruskal `F(rs)=4e⁻²=0.541`, proper distance `√(2 rs ε)=0.2` = ε healing length, antipodal `l=0` mode real `ω=1.186` (#129), absorbing `l=0` mode complex `ω=1.893−1.159i` (#130) — all mutually consistent. **One primitive, five faces:** the antipodal identification of the 5D Tangherlini horizon appears as (1) the `C` inner/outer swap (#63, `c₁→−c₁`), (2) the throat↔antithroat nucleation channel (#58), (3) the antipodal Kruskal map `(U,V,Ω)→(−U,−V,Ω̄)` (#128), (4) the l-parity unitary-mirror BC (#129), (5) the selector of the real stable-matter spectrum (#130) — **"Bulk Antipodal Mechanics" is the mechanics of this one identification.** **Epistemic ledger:** DERIVED — the throat's parent is a genuine curvature-regular D=5 Tangherlini vacuum (Ricci-flat, `S³` horizon=Hopf base, `k₅=D_bulk`), the coordinate singularity is removable, the antipodal identification fixes the l-parity BC/unitary mirror, the antipodal spectrum is real (stable matter) vs absorbing complex; POSTULATED — the antipodal identification itself (BAM's defining axiom), shown self-consistent not forced; OPEN — exact AdS scale `k=κ₅²/Λ₅` (#112), nucleation rate (#58/#88), global brane solution, idealised horizon QNM tower (`geometric_throat_arc_synthesis_probe`, PR #131) |
| **Throat↔antithroat nucleation rate** on the regular 5D background | **Antipodal instanton on a smooth cigar; S ∝ ln(1/ε) is the horizon tortoise divergence** | Closes the #131 lead open item: the dynamical nucleation rate, placed on the horizon-regular background (#128) and tied to the Majorana bounce arc (#87–#90). The transition (ΔL=2 Majorana/pair-production, #58) is the Kruskal region I↔III crossing (#128) via the odd `c₁→−c₁` instanton (#63); rate `Γ ~ [det(H)/det(H_free)]^{−1/2} e^{−S}`. **Smooth Euclidean cigar (Gibbons–Hawking):** the near-horizon `ds²_E ≈ dρ²+κ²ρ²dτ²` (`ρ=√(2rs(r−rs))`, `κ=1/rs`) is smooth — deficit `2π−κβ=0` — iff the imaginary-time period `β=2π/κ=2π rs`; so the nucleation temperature `T_nuc=1/β=1/(2π rs)=T_H` carries the closure quantum 2π (#127). **Action = horizon tortoise divergence:** the bounce tortoise length `L*(ε)=(rs/2)ln(1/ε)+const` (slope `rs/2=0.5`, verified to 4 digits), so `S∝ln(1/ε)` — the exact-horizon limit `ε→0` costs infinite length ⟹ `S→∞`, `m_ν→0` (the "rigid throat ⟹ massless ν" of #88, read off geometrically; regulated by the finite ε healing length #112). **Rate:** with `t∈[2π, k₅√(2π)]` (#89) and `ε~R_c³` (#112), `S≈15–18`, `m_ν=m_D e^{−S}~few meV` (#87/#90). **The prefactor closes the arc:** the one-loop prefactor is the #116 Tangherlini fluctuation determinant `1.574370` — #116 prefactor, #127/#128 stage, #58/#87–#90 bounce. Open (inherited): exact ε, absolute scale `κ₅²/Λ₅`, precise S/m_ν (#88–#90, #112) (`throat_antithroat_nucleation_rate_probe`, PR #132) |
| **Bulk scale ledger** for κ₅²/Λ₅ and ΔR | **The recurring absolute-scale residual = one bounded number; ΔR = the unit** | Consolidating ledger for the absolute bulk scale that surfaced open at every step (#57/#112/#127/#132). Counts the D=5 dimensionful content (`κ₅²[L³]` the 5D Newton constant, `Λ₅[L⁻²]` ⟺ `k=√(|Λ₅|/6)[L⁻¹]`, `L_AdS=1/k`; `λ_crit=6k/κ₅²[L⁻⁴]`; `R_MID, ΔR[L]`) into **three categories**: **(1) ΔR = the scale modulus** (`ΔR=R_OUTER−R_INNER=0.52 R_MID`) — the one dimensionful anchor the B4 theorem requires (#52), a proper invariant length (#53); it sets the length unit, geometry ratios `ΔR/R_MID=0.52`, `R_OUTER/R_MID=1.26` fixed — **units, not a residual**; **(2) √6 = the one fixed tuning** `λ_crit κ₅²/√|Λ₅|=√6≈2.449` (RS flatness, #57); **(3) the open bulk number** = the AdS scale `k·R_MID=R_MID/L_AdS` (`=κ₅²/Λ₅` in throat units) — **bounded ≲0.1** by the cavity correction `(k r)²~O(10⁻²)` (#127), so `R_MID≲L_AdS/10` (throat deep in the near-flat AdS region, why pure-Tangherlini #116/#127 is a good approx). **Ledger:** `{κ₅²,Λ₅} → {G=κ₅²/ΔR³ anchor} + {√6 fixed} + {k·rs open bounded ≲0.1}` with ΔR the unit ⟹ the recurring `κ₅²/Λ₅` residual is **one bounded dimensionless number**, not a multi-parameter mystery. Bounds and isolates it; does NOT pin `k·rs` (still the #112 residual) or add a free parameter (`bulk_scale_ledger_probe`, PR #133) |
| **Flavor hierarchy audit** from logarithmic throat bounce lengths | **Log-bounce governs the neutrino sector only (form/ordering, not value)** | Audits whether the three-generation flavor hierarchy follows from the logarithmic throat bounce lengths `L*(ε)=(rs/2)ln(1/ε)` (#88/#132), via tunnelling masses `m~e^{−S}`. **Mechanism:** `S=c·L*(ε)=c(rs/2)ln(1/ε)` ⟹ `m=m_0 e^{−S}=m_0 ε^{c rs/2}=m_0 ε^p` — the logarithm turns the exponential into a **power law** in the throat penetration depth ε (identity `e^{−cL*}=ε^p` verified). **Neutrino = the log-bounce sector:** the only genuine tunnelling sector (`k=0` chargeless, neck not EM-propped #86/#88), `m_ν∝ε^{4.8}` (#112); `ε_n∝1/χ_n` (#79) gives the right **ordering** (normal), but the steep power amplifies the modest χ_n spread — a ×2 ε spread → `2^4.8≈28×` in mass = the **×28 overshoot** (#113). Form/ordering governed, value residual. **Charged leptons NOT log-bounce:** Dirac (`c₁=±1` EM-propped, no tunnelling #86/#88), masses from the winding ladder `β·k²` (#71); `ln m` irregular (gen-diffs 5.33, 2.82 → ratio 0.53). **Quarks NOT log-bounce:** shell-resolving cavity overtones / `n_part` (#77–#80); `ln m` irregular (up-type ratio 0.77). ⟹ the flavor hierarchy is a **three-mechanism structure** (bounce ν, winding charged-lep, cavity quark), NOT a single log-bounce phenomenon. **Why residual:** `m∝ε^p` ⟹ `∂ln m/∂ln ε=p` ⟹ masses **hypersensitive** to the throat depth (×2 ε → 2^p mass), so the flavor values' irreducibility (#108) is a consequence of the exponential mass-action relation, not a separate mystery. Open: the ν value overshoot (#113), the charged/quark irregular magnitudes (the flavor puzzle #97/#107/#108) (`flavor_hierarchy_log_bounce_audit_probe`, PR #134) |
| **Antipodal-horizon exchange kernel** (matter-sector propagator) | **The antipodal-BC cavity resolvent: reciprocal, unitary, parity-graded** | Builds the matter-sector exchange kernel — the two-point Green's function / resolvent of the matter cavity operator (#116) with the antipodal horizon boundary data (#129). (The gauge-sector photon kernel `1/q²` is the separate PR #42–#44 `bam_exchange_kernel_probe`.) **Kernel:** `K_l(r,r';ω) = ⟨r|(H_l − ω²)⁻¹|r'⟩`, `H_l = −d²/dr*² + V_l` with the #129 antipodal BC (even-l Neumann / odd-l Dirichlet, Dirichlet shell wall); `H_l` exactly self-adjoint. **Spectral representation:** `K_l = Σ_n ψ_n(r)ψ_n(r')/(ω_n² − ω²)` — a sum over the **stable modes**, poles = the real #130 spectrum (mode sum = matrix resolvent to ~1e-14): the propagator is built as an exchange of stable modes. **Reciprocity:** `K_l(r,r') = K_l(r',r)` (self-adjoint ⟹ symmetric kernel, ~1e-14). **Unitary vs lossy — the boundary data decides:** antipodal (real BC) ⟹ Hermitian ⟹ real poles ⟹ **unitary** undamped kernel; absorbing (ingoing BC) ⟹ non-Hermitian ⟹ complex poles ⟹ **lossy** kernel (#130) — the propagator-level face of the unitary mirror (#129) and global CPT/unitarity (#64). **Angular parity grading:** `K(x,x') = Σ_l K_l(r,r';ω) C_l(Ω·Ω')`; under the throat↔antithroat exchange `Ω'→AΩ'`, `C_l(−Ω·Ω') = (−1)^l C_l(Ω·Ω')` ⟹ each l-channel carries the antipodal sign `(−1)^l` (even-l symmetric, odd-l antisymmetric under the C-swap #63) — the same `(−1)^l` that fixed the BC (#129/#134). Open: the interacting/multi-loop kernel (vertices, self-energy), absolute normalisation (#133), flavor residuals (#134) (`antipodal_horizon_exchange_kernel_probe`, PR #135) |
| **One-loop self-energy audit** for the antipodal matter kernel | **Finite real mass shift; lightest mode stays exactly stable; unitarity survives** | Audits the leading interacting correction to the #135 free kernel — the one-loop self-energy `Σ`. **Dyson dressing:** `G(s) = 1/(s − ω_k² − Σ(s))`, `s=ω²`; `Re Σ` = mass renormalisation, `Im Σ` = width. **One-loop Σ = the two-particle bubble:** `Σ_k(s) = Σ_{n≤m} c_{nm}|g_{knm}|²/(s − (ω_n+ω_m)² + i0⁺)`, with the cubic vertex `g_{knm} = ∫ ψ_k ψ_n ψ_m dr*` the triple overlap of the antipodal modes. **Lightest mode exactly stable:** by the optical theorem `Im Σ` is the two-particle phase space — nonzero only above a threshold `(ω_n+ω_m)²`; the lowest is `2ω_0`, and the lightest mode at `ω_0 < 2ω_0` has its pole `s=ω_0²=1.36` below `s_thr=(2ω_0)²=5.45` ⟹ `Im Σ_0(ω_0²)=0` ⟹ cannot decay (energy conservation), stays a sharp real-pole stable particle through one loop. **Finite mass shift:** `Re Σ_0(ω_0²)` converges with the mode cutoff (−0.277→−0.280 for cutoff 10→40), the residual UV piece being the #116 zeta/heat-kernel regularisation — a finite mass renormalisation (×coupling²), no UV catastrophe. **Unitarity survives + no horizon-absorption width:** `Im Σ ≤ 0` above threshold, `=0` below (optical theorem); because the throat is a unitary mirror (#129) there is **no** horizon-absorption contribution — the only width is genuine multi-particle decay (above `2ω_0`), vs the absorbing horizon's tree-level width on every mode (#130). One loop extends the tree-level stable spectrum (#130/#135). Open: the interaction vertex/coupling (modelled, not derived from S_BAM), higher loops, absolute normalisation (#133), flavor residuals (#134) (`antipodal_kernel_one_loop_self_energy_probe`, PR #136) |
| **Cubic vertex ledger** for the antipodal matter kernel | **Antipodal Z₂ selection rule + geometric shape DERIVED; coupling INPUT** | Ledger for the cubic vertex `g_{knm} = ∫ψ_kψ_nψ_m` the #136 self-energy modelled — separating its derived structure from its input magnitude. **Factorises:** `V = λ · [∫_{S³} Y_{l1}Y_{l2}Y_{l3} dΩ] · [∫ψ_kψ_nψ_m dr*]` (angular × radial × coupling). **Angular selection rule (DERIVED):** the S³ harmonic triple integral is nonzero only if **(a) `l1+l2+l3` even** — the antipodal parity: under `x→−x` (the throat↔antithroat C-swap #63) `Y_l→(−1)^l Y_l`, so `(−1)^{Σl}=+1` over the inversion-symmetric S³ — **AND (b) the triangle** `|l1−l2|≤l3≤l1+l2` (SO(4)). Verified exactly via the S³ monomial integral (odd-Σl → 0; triangle-violating → 0; allowed → nonzero). **The parity rule IS the arc's `(−1)^l`** that fixed the BC (#129), graded the kernel (#135), and sorted the flavor sectors (#134); the #136 bubble connects only even-Σl triples. **Radial overlap (DERIVED shape):** `∫ψ_kψ_nψ_m dr*` is a definite geometric number (the #116 cavity modes), totally symmetric in (k,n,m) (Bose, ~1e-14) and real. **INPUT/residual:** the overall coupling `λ` (dimensionless, #136 set it to 1), and whether S_BAM (#115–#122) generates the cubic term at all. So the vertex STRUCTURE (selection rule + geometric shape + symmetry) is BAM-native; only its MAGNITUDE is input. Open: `λ` not derived, quartic/higher vertices, S_BAM cubic generation (`cubic_vertex_ledger_probe`, PR #137) |
| **Quartic vertex ledger + bounded interaction audit** | **Same antipodal Z₂ rule; positive overlap ⟹ bounded-below stable vacuum** | Extends the #137 cubic ledger to the quartic vertex and audits whether the matter interaction is bounded below. **Factorises:** `V_4 = λ_4·[∫_{S³}Y_{l1}Y_{l2}Y_{l3}Y_{l4}dΩ]·[∫ψ_kψ_lψ_mψ_n dr*]`. **Quartic angular rule (DERIVED):** nonzero only if **(a) `l1+l2+l3+l4` even** — the **same antipodal Z₂** as the cubic (#137, the `x→−x` C-swap #63) — **AND (b) a common SO(4) channel** `∃L∈[|l1−l2|,l1+l2]∩[|l3−l4|,l3+l4]`. Verified exactly (odd-Σl→0); the Z₂ parity persists cubic→quartic. **Positive self-overlap (DERIVED):** `g_4 = ∫ψ_k⁴ dr* > 0` manifestly (integral of a 4th power; `1.03, 1.02,…`). **Bounded interaction ⟹ stable vacuum:** the single-mode potential `V(a) = ½ω²a² + (λ_3 g_3/6)a³ + (λ_4 g_4/24)a⁴` has a⁴ coefficient `λ_4 g_4/24 > 0` (g_4>0, λ_4>0) ⟹ `V→+∞` as `|a|→∞` for **any** cubic ⟹ bounded below, a **stable vacuum** (the cubic only tilts, never unbounds — verified `V(±10⁴)>0` up to λ_3=200). **Boundedness = measure convergence (#122):** a bounded-below action is exactly the condition for `∫Dμ e^{−S}` to converge (established non-perturbatively, #122), so the positive quartic is required by, not added to, the measure's existence — extending the stability thread (#130 stable spectrum, #136 unitary self-energy, #138 bounded vacuum). Open: coupling magnitudes `λ_3, λ_4` input (sign `λ_4>0` from #122), quintic/higher vertices, S_BAM generation (`quartic_vertex_bounded_interaction_probe`, PR #138) |
| **Antipodal matter interaction synthesis** (capstone of #129–#138) | **Two threads, one postulate: the antipodal Z₂ + the unitary mirror** | Capstone re-verifying a keystone from each arc member together and organising the whole arc into two threads from the single antipodal postulate. **Keystones (re-run together):** `Y_l` parity `[1,−1,1,−1]=(−1)^l` (#129); antipodal fundamental real `≈1.17` vs absorbing complex `≈1.89−1.16i` (#130); kernel reciprocity `~1e-14`, real poles (#135); lightest-mode `Im Σ≈0` (stable, #136); `∫ψ⁴=1.03>0` (bounded vacuum, #138) — all mutually consistent. **Thread A — the antipodal Z₂ `(−1)^l`:** the C-swap inversion `x→−x` (#63) carrying `Y_l→(−1)^l Y_l` fixes the BC (#129), grades the exchange kernel (#135), and selects the cubic (#137) and quartic (#138) vertices (`Σl` even). **Thread B — unitarity/stability:** the unitary mirror (#129) ⟹ real stable spectrum (#130) ⟹ unitary reciprocal propagator (#135) ⟹ unitarity-preserving self-energy + stable lightest mode (#136) ⟹ bounded-below vacuum (#138 = #122 measure convergence) — stable at every order. **One postulate, two faces:** the real l-parity BC (#129) IS both the Z₂ grading and the unitary mirror. **Epistemic ledger:** DERIVED (given the antipodal BC) — the Z₂ selection structure + the unitary stable propagator/self-energy/vacuum; POSTULATED — the antipodal identification (#128, self-consistent not forced); INPUT — the coupling magnitudes `λ_3, λ_4` (sign `λ_4>0` from #122); OPEN — S_BAM vertex generation, higher loops/vertices, scale (#133), flavor (#134) (`antipodal_matter_interaction_synthesis_probe`, PR #139) |
| **S_BAM vertex generation** (vertices derived, not modelled) | **Vertices = action Taylor coefficients; selection rule = antipodal Ward identity** | Closes the #137–#139 open item (the vertices were modelled). **Vertices = Taylor coefficients of S_BAM:** expanding `S_BAM[φ_cl+φ] = S_cl + S_2 + S_3 + S_4 + …`, `S_n = (1/n!)∫(δⁿS/δφⁿ)φⁿ` — `S_2` the #116 determinant / #135 propagator, `S_3=(g/3!)∫φ³√g` ⟹ `∫ψ_kψ_nψ_m` (#137), `S_4=(λ/4!)∫φ⁴√g` ⟹ `∫ψ⁴` (#138); the geometric (non-quadratic) S_BAM generates the tower, a free action none. **Selection rule = the antipodal Ward identity:** the S_BAM measure carries the `Diff(S¹)⋉U(1)⋉Z₂` quotient (#74), whose Z₂ is the C-swap `A: x→−x` (#63/#128); under `A` a mode amplitude transforms `a_l→(−1)^l a_l`, so a vertex picks up `(−1)^{Σl}` and is `A`-invariant ⟺ `Σl` even. Since S_BAM is `A`-invariant, every vertex has `Σl` even — the #137/#138 rule as a **Ward identity**, not a modelling choice (verified: `A`-invariance matches the explicit S³ odd-Σl vanishing). **Quartic sign = measure consistency:** `∫Dμ e^{−S}` exists (reflection-positive ⟹ unitary kernel #135; convergent #122) ⟺ `S` bounded below ⟺ `λ_4>0` (#138) — the positive sign is required, not chosen, realised by `∫ψ⁴>0`. **Ledger:** DERIVED — vertex existence (action expansion), Σl-even selection (antipodal Ward identity), positive quartic sign (measure consistency #122); INHERITED — the coupling magnitudes `g, λ` (the action's higher derivatives), carrying the `κ₅²/Λ₅` scale (#133). So the vertex STRUCTURE is generated; only the MAGNITUDES inherit #133. Open: exact S_BAM form, coupling magnitudes, scale (#133), higher vertices, flavor (#134) (`s_bam_vertex_generation_probe`, PR #140) |
| **Gauge–matter coupling** from the antipodal throat | **Minimal coupling, Z₂-selected vertex, charge conserved; only α input** | Joins the gauge sector (the U(1)_Hopf photon `1/q²`, #42–#44) to the matter sector (the antipodal cavity modes, #129–#140) at the throat. **Minimal coupling:** matter of Hopf charge `c₁` (`|c₁|=1`, #58/#74) couples through `D_μ=∂_μ−ic₁A_μ`, vertex `c₁∫A_μ j^μ`. **The C-swap = inversion × charge conjugation:** the antipodal map `A:x→−x` (#63) acts at once on the matter harmonics (`Y_l→(−1)^l Y_l`, #129/#140) and the Hopf charge (`c₁→−c₁`, #63) — one operation, two effects, so **the throat is the particle↔antiparticle (C) surface** (#63/#64). **The gauge vertex inherits the antipodal Z₂ selection rule:** the photon-matter-matter angular part is the triple overlap `∫Y_{l_γ}Y_{l₁}Y_{l₂}` (the cubic-vertex structure #137 with a gauge leg); antipodal invariance (the #140 Ward identity) ⟹ `Σl=l_γ+l₁+l₂` even (verified: odd-Σl forbidden). **U(1) charge conserved at the throat:** the unitary mirror (#129, zero net matter flux) conserves the charge flux; with the C-swap flip outgoing charge re-emerges as the conjugate ⟹ `Σc₁=0` (#58) — the gauge face of the mirror. **Coupling strength = α:** the structure (covariant derivative, Σl-even vertex, charge conservation, C-surface) is derived; the strength is the EM coupling `α` (the "137 problem", #105), the universal **input** residual. (QFT on the classical throat, not quantum gravity.) Open: `α` (#105/#108), EM normalisation, higher gauge vertices, scale (#133), flavor (#134) (`gauge_matter_coupling_probe`, PR #141) |
| **Gauge Ward identity + current conservation audit** | **Gauge invariance is the gauge face of the unitary mirror; only α input** | Audits the consistency of the #141 coupling — current conservation, the Ward–Takahashi identity, photon masslessness — and ties them to the unitary mirror (#129). **Conserved Noether current:** global U(1)_Hopf ⟹ `j^μ=i(ψ*∂^μψ−…)`, `∂_μ j^μ=0`; stationary mode ⟹ `ρ=2ω_n|ψ_n|²` static. **Current conserved at the throat:** the antipodal modes are **real** (#135) ⟹ radial charge current `j^r∝Im(ψ_n*∂_rψ_n)=0` **exactly** (verified) — no charge flux through the throat, charge static and conserved; this IS the zero-flux unitary-mirror property (#129). **An absorbing throat would break it:** complex quasinormal modes (#130) carry `j^r≠0` (verified `≈0.014` at the throat) — charge leaks into the horizon, current conservation fails ⟹ gauge invariance broken, so **gauge invariance REQUIRES the antipodal throat**. **Ward–Takahashi:** `q_μ Γ^μ=S⁻¹(p_out)−S⁻¹(p_in)` ties the gauge vertex (#141) to the matter inverse propagator (#135) — the gauge coupling fixed by the matter dynamics. **Photon masslessness:** transversality `q_μ Π^μν=0` ⟹ no photon mass ⟹ the `1/q²` photon (#42–#44) is protected. **One postulate, both:** current conservation, the Ward identity, and masslessness all follow from the unitary antipodal throat (#129) — the same real/self-adjoint/zero-flux structure that gave the stable spectrum (#130), unitary propagator (#135), stable self-energy (#136), and bounded vacuum (#138). Gauge invariance is not an extra assumption; only `α` (#105) is input. (QFT on the classical throat, not quantum gravity.) Open: `α` (#105/#108), higher-order Ward identities, running of `α`, scale (#133), flavor (#134) (`gauge_ward_identity_probe`, PR #142) |
| **Alpha normalization ledger** for the gauge–matter coupling | **Charge quantum + 1/2π measure + running derived; the value α is one residual** | Consolidating ledger for the EM coupling normalisation `α` (the strength left input by #141/#142) — parallel to the bulk-scale ledger (#133). **How α enters:** `A_EM=α·ℏc/2` (#105), the vertex strength `∝c₁²α` (#141), the Schwinger anomaly `a=α/2π` (#74) — a single dimensionless number. **Charge quantum DERIVED:** `|c₁|=1`, the integer Hopf number (#58/#74) — charge quantisation topological, the charge unit geometric. **1/2π measure DERIVED:** in `a=α/2π` the `2π` is the closure-quantum loop measure (#74); BAM fixes the `1/2π`, only `α` is the input prefactor. **Running DERIVED, value not:** the RG flow (vacuum polarisation, transverse by the #142 Ward identity) is derived — BAM derives **how** `α` runs; the boundary value `α(μ_0)≈1/137` is input — not **where** it starts (#105). **Value = one EM residual (the 137 problem):** a fit-independent scan of `α⁻¹=137.036` vs the closure numbers (`2π`, `k₅`, `50π`) finds **no clean match** — the near-misses (`50π−20=137.08`, `4·k₅²+37=137.0`) each need an ad-hoc additive `O(20–37)` term (fits, not derivations, the #107/#108 failure mode); `α` plausibly irreducible like `√σ/m_e` (#108). So the EM sector contributes exactly **one** residual, the value `α` ∈ `{n_part, √σ/m_e, ε, α}` (#104). (QFT on the classical throat, not quantum gravity.) Open: the value `α` (137 problem), EM normalisation, scale (#133), flavor (#134) (`alpha_normalization_ledger_probe`, PR #143) |
| **One-loop photon vacuum polarisation** and the running of α | **Ward-protected; massless; screening; flat-limit log slope α/3π — α(μ₀) input** | Computes the running that #142/#143 classified as derived but never computed — the one missing one-loop two-point function (matter Σ was #136). **Bubble:** Π = charged-pair loop over the antipodal modes (#135), vertex `v_nm = ∫φ_γψ_nψ_m dr*` (the #137/#141 triple overlap, one photon leg), density `ρ_nm = c\|v_nm\|² ≥ 0`; the photon couples only to **even-Σl** pair channels (the #141 antipodal Z₂ rule, re-verified exactly). **Cavity Ward identity COMPUTED:** the diamagnetic +1 cancels the paramagnetic sum, `1 − S = 3.1e-05` (the TRK sum rule in disguise; the #142 structural identity made quantitative) ⟹ `Π(0) = 0`: the photon stays **exactly massless**, the `1/q²` kernel (#42–#44) protected through one loop. **Absorbing counterfactual breaks it:** complex pair thresholds (#130) ⟹ `Im Π ≠ 0` below threshold (photon absorption width, `−0.042` vs antipodal `−6e-09`) and the real-mode cancellation is lost — gauge protection REQUIRES the unitary antipodal throat (#129/#142). **Unitarity:** `Im Π = 0` below the lowest pair threshold `(2ω_0)²` (the #136 pattern, now on the photon). **Screening/running COMPUTED:** the Ward-protected `Δ(Q²) = Σρ Q²/(s(s+Q²))` is monotone ⟹ `α_eff` increases with `Q²` (the QED direction, discrete pair thresholds = the cavity analogue of lepton thresholds); the same dispersion machinery fed the flat 4D pair density reproduces the textbook log running with slope `α/3π` to 0.97%. **No 137-hunting** (the #107/#108 discipline): the boundary value α(μ₀) stays the one EM input (#143). Open: higher loops, the 4D tensor `Π^μν`, normalisation (#133), flavor (#134) (`vacuum_polarization_running_probe`, PR #144) |
| **Charge non-renormalization** `Z₁ = Z₂` (the renormalization triangle closed) | **Ward identity computed: dressed charge = c₁ exactly, universal; e = √Z₃·e₀** | With Π in hand (#144) all three renormalization constants exist: Z₂ (#136 Σ), Z₃ (#144 Π), Z₁ (#141/#142 vertex); `e = (Z₂/Z₁)·√Z₃·e₀`, so Z₁ = Z₂ ⟺ the matter sector does not renormalize charge. **Computed on the cavity:** charged χ (odd-l Dirichlet tower, `c₁ = 1`) × neutral φ (even-l Neumann) with the charge-conserving #136/#137 cubic vertex; `Σ′(s₀)` two independent ways (spectral sum vs finite difference, agree ~1e-12), `Z₂ = 1/(1−Σ′) ≈ 0.986` — genuinely renormalized. **The Ward identity term by term:** the q=0 photon insertion doubles the charged propagator, `Λ(0) = Σ c₁\|g\|²/(s₀−s_nm)² = −c₁Σ′(s₀)` to machine precision (neutral line contributes 0) — the #142 Ward–Takahashi identity at one loop ⟹ `Z₁ = Z₂`. **Dressed charge exact + universal:** `F(0) = Z₂(c₁+Λ) = c₁` to machine precision across species (l_χ ∈ {1,3}, l_φ ∈ {0,2}, g ∈ {0.5,0.7,1}): Z₂ varies 0.9855–0.9963, `F(0) − c₁ = 0` identically — each sector’s self-interaction cancels out of its own charge (why all generations k ∈ {1,3,5} carry the same `\|c₁\| = 1`, #71); charge renormalization collapses to `e = √Z₃·e₀` — **the running of α is purely the #144 vacuum polarisation**. **Counterfactual:** a charge-violating vertex (`c′ ≠ c₁`) shifts F(0) (−0.003 to −0.014) and makes it species-dependent — the protection is exact charge conservation at the unitary throat (`Σc₁ = 0`, #58/#141/#142); the absorbing throat leaks charge and loses it. Open: q ≠ 0 form factors, higher loops, normalisation (#133), flavor (#134) (`charge_non_renormalization_probe`, PR #145) |
| **Finite-momentum charge form factor** F(q) on the antipodal cavity | **Bethe sum rule = q² at every q; F(0) = c₁ anchored; charge radius GEOMETRIC** | Supplies the q ≠ 0 structure #145 left open. **Form factor:** the dressed charge density `ρ(x)` integrates to `c₁` exactly; `F(q) = ∫ρ e^{iq(x−x̄)}dx` falls monotonically — the throat charge has spatial structure (not pointlike). **Finite-q Ward identity:** the Bethe sum rule `Σ(E_m−E_n)\|⟨m\|e^{iqx}\|n⟩\|² = q²` verified to ~1e-4 across `q ∈ [0.5, 10]` (the double commutator `[e^{−iqx},[H,e^{iqx}]] = 2q²` is V-independent) — the finite-q generalization of the #144 TRK sum rule (its q²→0 limit), current conservation at EVERY momentum transfer. **Dressed density:** the one-loop dressed state reproduces the #145 Dyson Z₂ exactly (`1/(1+Σa²) = 1/(1−Σ′)`, machine precision — the two one-loop pictures agree) and its total charge is `c₁` exactly at every coupling (the #145 anchor in real space). **Charge radius GEOMETRIC:** `r_c = 0.2649` (tortoise units) from the density variance = from the small-q fall-off of F; the one-loop cloud shifts it only ~9e-5 (× coupling²) ⟹ the radius is the bare cavity mode profile — finite, no UV divergence (the form-factor face of the #55 finite self-energy), set by the classical geometry with the QFT dressing a small correction (geometry → fields). **Counterfactual:** a charge-violating cloud (`c′ ≠ c₁`) shifts the total charge — the protection is exact charge conservation at the unitary throat (#58/#141/#142/#145). Open: recoil, F₁/F₂ (g−2, #62), higher loops, normalisation (#133), flavor (#134) (`charge_form_factor_probe`, PR #146) |
| **Electric/magnetic form-factor decomposition** `F₁/F₂` — the **EM gauge-arc capstone** | **Gordon split exact; Ward pins F₁ only; g = 2 + α/2π keystones re-verified; radii geometric** | Assembles `Γ^μ = γ^μF₁ + iσ^{μν}q_νF₂/2m` on the cavity and capstones #141–#146 (+ keystones #61/#62, the #131 convention). **Gordon decomposition:** `ū′γ^μu = ū′[(p+p′)^μ + iσ^{μν}q_ν]u/2m` verified with explicit Dirac spinors to ~1e-15 — the E/M split is the Dirac algebra of the #141 minimal coupling, not an ansatz. **Why charge is exact and the moment is not — ONE identity:** the Ward contraction kills the F₂ term twice (`q_μσ^{μν}q_ν = 0` exactly; on-shell `ū′q̸u = 0` ~1e-16) ⟹ #142/#145 pins F₁ only: `F₁(0) = c₁` exact + coupling-independent, while F₂ is gauge-FREE and dresses at every loop. **Keystones re-verified together:** tree `(σ·D)² = D² − σ·B` (~1e-6) ⟹ `g_s = 2`, `F₂(0) = 0` (#61); loop Schwinger simplex `∫₀¹2z dz = 1` (0.9999998) ⟹ `a = α/2π = 0.00116141` vs measured `a_e = 0.00115965` (+0.15%, the α² Sommerfield term and beyond), `g = 2.0023228` vs `2.0023193` (#62). **Sachs assembly:** `G_E(0) = c₁ = 1` EXACT (Ward-pinned), `G_M(0) = 1 + α/2π = 1.0011614` DRESSED, `g/2 = G_M(0)/c₁`; the magnetization rides the same charged-mode profile ⟹ `r_M = r_E = 0.2649` (geometric, #146) and `G_M/G_M(0) = G_E/G_E(0)` (scaling, minimal model). **The arc, one primitive:** every face #141→#147 derives from the unitary antipodal throat with integer Hopf charge; the single EM input is the value α(μ₀) (#143). Open: α² term, r_E−r_M splitting, recoil, normalisation (#133), flavor (#134) (`em_form_factor_decomposition_probe`, PR #147) |
| **Bulk-scale residual audit** for `k·r_s` (the #133 open number) | **Two-sided bracket: 0 < k·r_s ≲ 0.0064–0.070 — the #133 estimate tightened ~16×; value residual** | Makes the #133 bound quantitative. **Background under the operator:** the #127 interpolating `f_k = 1 − r_s²/r² + k²r²` (Einstein, `Λ₅ = −6k²`) with `V_l = f[l(l+2)/r² + (3/2r)f′]` — reduces to the #116 Tangherlini potential at machine precision at k = 0, and the k = 0 pinhole operator reproduces the documented −2.2% γ-lock residual (`Σ V_max = 22.02` vs locked 22.5) — the audit stands on the locked machinery itself. **Quadratic scaling DERIVED:** ω(1,0), ω(0,0), and the pinhole sum all shift as `(k·r_s)²` (fitted exponents 1.98–2.00 across a decade). **Spectrum bound:** sensitivities `c ≈ 9.9` (ω) / `4.5` (pinhole) convert the locked precisions into `k·r_s ≤ √(tol/c)`: the 0.04% Compton-bridge closure ⟹ **≲ 0.0064**; the 2.2% pinhole lock ⟹ ≲ 0.070 — the throat sits deep in the near-flat AdS region (why pure Tangherlini #116/#127 worked so well, now quantified). **Lower bound:** `k → 0` ⟹ `λ_crit → 0` (#57) ⟹ `B = 4πσ → 0` (#56) ⟹ `R* = (A/2B)^{1/3} → ∞`, `E(R) = A/R` monotone (no minimum, computed) — a static throat requires `k > 0`. **Bracket:** `0 < k·r_s ≲ 0.006–0.07`, the #89 two-sided ε pattern: structure and bracket derived, value residual (`κ₅²/Λ₅`, #112); input budget unchanged. Open: the value (absolute normalisation), global brane solution (#127) (`bulk_scale_k_rs_audit_probe`, PR #148) |
| **Neutrino log-bounce sensitivity audit** and ε_n overshoot bracket | **Required profile data-pinned to ~0.3%; spread ∈ [1.32, 1.44]/step; χ-driven law and any single power law EXCLUDED** | The #148-pattern audit applied to the #113 ×28 overshoot. **Keystones re-verified:** `L*(ε) = (r_s/2)ln(1/ε)` (slope re-fit 0.500, #88/#132) and the power-law identity `e^{−cL*} ∝ ε^{c·r_s/2}` constant over three ε decades (`p = 4.8`, #112). **Hypersensitivity INVERTED:** forward, the steepness produced the overshoot; inverted, `δln ε = δln m/p` compresses the oscillation-data errors (±2.8%/±1.1% on Δm²₂₁/Δm²₃₁) into **~0.3–0.4%** on the required ε ratios — the residual is sharply localized, not fuzzy (the #148 inversion). **Attribution bracket:** the data fix only `m_D,n·ε_n^p`; between the pure-bounce (uniform m_D) and #113-implied (1.88, 1.48) endpoints the required spread is `ε₃/ε₂ ∈ [1.32, 1.44]` — the principled χ-driven `χ₂/χ₃ = 2.49` sits far outside both (×1.7–1.9 in ε ⟹ ×14–21 in mass at this data vintage; #113: ×28). **Power-law exclusion (scanned):** per-pair exponents in `ε_n ∝ χ_n^{−q}` disagree under BOTH attributions (q ratio 1.5–2.1); best single `q = 0.32` still misses ×1.38 in mass (~25× the data bracket); `q = 1` overshoots ×21 — the spread is NOT a power law in the boundary stress (#113 sharpened). **Consistency:** normal ordering (derived); `Σm_ν ≈ 61 meV` inside the program’s own 59–65 meV window. Residual: the gentle profile’s origin (plausibly mixing/anarchy, #92); budget unchanged. Open: deriving the profile, the #112 anchor, the m_D attribution (`neutrino_eps_n_overshoot_bracket_probe`, PR #149) |
| **Residual-bracket synthesis** and the categorized **input-budget ledger** | **One anchor + 4 residuals + 2 brackets + flavor puzzle — constant since #104/#125; #144–#149 added ZERO inputs** | The synthesis capstone for the program’s input accounting (#104/#105/#106/#107/#108/#123–#125/#133/#143/#148/#149 consolidated; one keystone re-verified per category, the #131 convention). **The categorized budget:** 1 dimensionful ANCHOR (`G` → `ΔR = 0.52·R_MID` unit, B4-mandatory; `√6` a FIXED derived tuning #57); 2 UNIVERSAL residuals (`α`, `√σ/m_e` — structure/running derived, values scan-excluded: best principled candidates `k₅³+2π` −4.2% and `2π·k₅³` −5.4%, every sub-% match ad-hoc, scans re-run); 2 PROGRAM residuals (`n_part` — APS doubling topological, value compensator, the #107 circularity re-derived `4n_part−100 ∈ [764, 920]` vs fixed 830; `ε` — order-of-mag derived, window `[2π, k₅√(2π)]`); 2 BRACKETED sub-residuals (`k·r_s ∈ (0, 0.0064–0.070]` #148, sensitivity re-checked `c = 9.86`; `ε_n` spread `[1.32, 1.44]`/step ~0.3% #149, inversion re-derived); 1 UNIVERSAL open problem (the flavor puzzle); and the lepton sector as the NO-residual contrast (`N = 4k₅² = 100` fully derived, #124). **The no-loose-knobs claim, checkable:** #144–#149 — six probes, zero inputs added; the budget is the SAME as #104/#125 while the derived ledger grew by the full one-loop EM sector and two bracket audits. Every residual row carries derived structure: a residual here is a number boxed by geometry, not a free knob. The consolidated table is added to `docs/THESIS.md`. Scope: organizes, does not remove (#125’s honesty) (`residual_bracket_synthesis_probe`, PR #150) |
| **Mixing/anarchy origin of the ε_n profile** (the #149 hypothesis, tested) | **POSITIVE: channel dominance makes the measured r₃₂ anarchy-natural (~77th pct), grows large mixing, and predicts m₁ ≈ 0.04 meV** | Builds the test #149 pointed at. **Model (hierarchical inputs all derived):** seesaw in the overtone basis `M_ij = m_D,i m_D,j·c_ij·G_ij(β)` — Dirac growth = the #91 cavity floors (= the #149 m_D-attribution endpoint to <1%, a verified identity), channel suppressions = the χ-driven compliances through the p = 4.8 bounce, `c_ij` = anarchic O(1) cross-channel overlaps (#91), and β interpolating the pair-tunneling saddle from FACTORIZED to CHANNEL-DOMINANT (each element tunnels through the widest neck available to the pair). **The measured ratio selects channel dominance:** β = 0 re-derives the #113/#149 overshoot in matrix form (ensemble median r₃₂ ≈ 113; observed 5.66 at the **0.1th percentile — excluded**); β = 1 collapses the steep hierarchy out of the heavy pair (median ≈ 2.8; observed at the **~77th percentile — natural**). **One rule, two observables:** the same β grows the mixing indicator 0.085 → 0.44 (aligned → large-mixing anarchic, the #91 cross-channel consistency). **The unmeasured ratio becomes a falsifiable prediction:** r₂₁ stays ≈ 200 ⟹ `m₁ ≈ 0.04 meV`, `Σm_ν ≈ 58.8 meV` (mixing solution) vs `m₁ = 2.08`, `Σ = 61.1` (#112 uniform anchor) — a ~2 meV cosmology discriminator, m_ββ shifted; both normal ordering. **Residual relocation:** the three-number profile (#149) → derived compliances + derived floors + one discrete saddle rule + an anarchic O(1) draw — ratios become percentile-natural statistics (the flavor puzzle’s BAM face, localized); no new continuous knob, the #150 budget unchanged. Open: derive the saddle rule from the bounce path integral; explicit PMNS angles/CP; the anarchic draw (`eps_n_mixing_anarchy_origin_probe`, PR #151) |
| **Bounce path-integral derivation** of the channel-dominant saddle | **DERIVED — the #151 β knob retired: mouth conversion ⟹ A_nm ≍ O(1)·e^{−min(S_n,S_m)}; counterfactual flips the rule** | Supplies #151’s lead open item. **Two-path decomposition:** the conversion vertex (the anarchic O(1) overlap #91) has support only where the overtones COEXIST — the cavity mouths; the neck interior is the single-channel tunneling region (#88/#132) ⟹ `A_nm = c_near·e^{−S_m} + c_far·e^{−S_n}` (convert-then-tunnel ⊕ tunnel-then-convert), dominated by the cheaper segment — channel dominance; mid-neck conversion would give the factorized `e^{−(S_n+S_m)/2}` but the vertex has no support there. **Exact computation:** a controlled 3-channel double well (WKB actions 15.4/11.4/8.4, splittings spanning ×2000) with mouth-localized coupling, solved exactly + Löwdin-projected (extraction faithful to <10%): `t/(Δ_max/2)` CONSTANT across pairs (×1.22 spread — the O(1) conversion factor) while `t/(Δ_geo/2)` varies ×8.65 (= `e^{\|ΔS\|/2}`, the two-path prediction). **The counterfactual decides:** vertex moved INSIDE the barrier ⟹ the rule FLIPS (t/geo constant ×1.59; t/max varies ×6.79) — the vertex location selects the saddle, and #151’s data exclusion of factorized corroborates mouth conversion (the BAM cavity/neck structure itself). **One vertex:** `t ∝ W₀` exactly (<1%). **Closure:** the #151 chain stands derived-footed — r₃₂ natural (~77th pct), mixing 0.44, `m₁ ≈ 0.04 meV / Σm_ν ≈ 58.8 meV` falsifiable; the β knob is RETIRED (one modelling assumption removed, zero inputs added; #150 budget unchanged). Open: the full 5D bounce path integral, the anarchic prefactor distribution, PMNS angles/CP (`channel_dominant_saddle_derivation_probe`, PR #152) |
| **PMNS angle extraction** from mouth-localized cross-channel overlaps | **All three angles anarchy-natural (62/56/27th pct); CP generic; the e-row hierarchy-protected — why θ₁₃ is small while θ₂₃ is large** | Assembles `U = R₂₃(φ_ℓ)·O_geom·U_ν` from the derived #151/#152 structure: U_ν = the channel-dominant complex anarchic ensemble; `O_geom` = the computed winding↔overtone mouth-overlap rotation (near-diagonal, ~8–13°); φ_ℓ = ONE charged-side μ–τ rotation. **The two failure modes bracket the structure:** near-diagonal O alone ⟹ θ₁₂/θ₁₃ natural but θ₂₃ far too small (98th pct); fully anarchic O ⟹ θ₂₃ natural but θ₁₃ far too large (≤7th pct) — the data select a specific intermediate. **The exact resolution:** a left μ–τ rotation leaves `sin²θ₁₂`/`sin²θ₁₃` EXACTLY invariant (machine zero — never touches the e-row) and moves only `sin²θ₂₃` ⟹ the data demand exactly ONE charged rotation — and it is the one the winding hierarchy PERMITS (`m_μ/m_τ = 0.060`, ×12 less hierarchical than `m_e/m_μ`; an O(m_τ) off-diagonal ⟹ O(1) left μ–τ rotation while the e-row stays protected). The natural window is broad (φ_ℓ ∈ ~[25°, 65°], ~45% of a uniform draw — no fine-tuning). **Assembled (φ_ℓ = 45°):** `sin²θ₁₂` observed at the 62nd percentile, `sin²θ₂₃` at the 56th, `sin²θ₁₃` at the 27th — the full observed point anarchy-typical; **CP generic** (median `\|J\| = 0.015` vs data max 0.033; `P(\|J\| > 0.01) = 61%`) — the "generic CP" claim quantified at the PMNS level. **Predictions:** θ₁₃ not-too-small preserved (e-row protection); no θ₂₃ octant preference; the #151/#152 `m₁ ≈ 0.04 meV / Σm_ν ≈ 58.8 meV` prediction unchanged. Structure derived, values statistical (the anarchic draw — the localized flavor residual); no new input (#150 budget unchanged). Open: derive the charged-side matrix; Majorana phases/m_ββ; the CKM intra-channel analogue (`pmns_angle_extraction_probe`, PR #153) |
| **Majorana phase and m_ββ prediction** from the PMNS flavor ensemble | **m_ββ ≈ 3 meV (68% [1.5, 5.9]); EXACTLY φ_ℓ-invariant; generic Majorana phases; detection > ~10 meV falsifies** | Completes the neutrino-sector card (#153’s lead open item). **Exact shortcut:** `m_ββ = \|(W M W^T)_ee\|` — the (e,e) element of the flavor-basis Majorana matrix (no mixing-matrix approximation, all phases included); the Takagi decomposition cross-checks it term by term (`Σt_i = M_fl,ee` to ~1e-12); each draw rescaled to the measured `m₃ = 50.14 meV`. **Exact invariance:** the charged-side μ–τ rotation never touches the e-row and `M_fl,ee` depends only on the e-row of W ⟹ m_ββ is EXACTLY φ_ℓ-independent (machine zero across 3000 draws) — the one modelled O(1) angle of #153 drops out entirely; m_ββ is MORE robust than the angles. **Prediction:** self-consistent median 3.2 meV (68% [1.5, 5.9], 95% [0.5, 8.7]); data-anchored 2.9; conditioned on data-compatible spreads 3.1 — robust. The few-meV scale is structural: the light m₁ (#151/#152, median 0.074 meV — negligible) makes m_ββ a TWO-TERM interference `\|t₂ + t₃\|`. **Falsification card:** P(< 1 meV) = 7.9% (cancellation uncommon); P(> 10 meV) = 0.5%; **P(> 20 meV) = 0** — a detection above ~10 meV falsifies the ensemble; ton-scale experiments predicted to see nothing or a floor-level signal; the earlier `m_ββ ≲ 8 meV` claim sharpened into a distribution. **Generic Majorana CP:** `P(\|Φ₂₃\| > π/2) = 69%` — the Majorana-sector face of #153’s generic Dirac CP. **The card complete:** NO + m₁ ≈ 0.05 meV + Σm_ν ≈ 58.8 meV + anarchy-natural angles + generic CP (Dirac & Majorana) + m_ββ ≈ 1.5–6 meV. No new input (#150 budget unchanged). Open: sharpen the O_geom e-row; CKM analogue; joint neutrino-sector test (`mbb_majorana_phase_prediction_probe`, PR #154) |
| **CKM intra-channel analogue** from mouth-overlap alignment | **OUT-OF-SAMPLE PREDICTION, zero new inputs: every element ≤ ×2.0; V_cb & V_ts within 10% (stiff); the PMNS/CKM dichotomy quantified** | The quark mirror of #153, computed from the LOCKED quark Hamiltonian (calibrated on the six masses ALONE). **Construction:** `partition_mixing = 0` ⟹ the 6×6 is exactly block-diagonal in the Z₂ partition: (+) = (u,c,t), (−) = (d,s,b) over the shared shells k = 1,3,5; `V_CKM = U₊†U₋`, unitary to machine precision. **Prediction:** `V_us = 0.112` (obs 0.225, ×0.50), `V_cb = 0.0377` (obs 0.0418, **×0.90**), `V_ub = 0.0020` (obs 0.0037, ×0.55), `V_td = 0.0063` (×0.73), `V_ts = 0.0372` (**×0.91**) — every element within ≤ ×2.0, the heavy pair at 10%, hierarchy `\|V_us\| > \|V_cb\| > \|V_ub\|` exact. **Mechanism:** the up sector aligned to 0.008, the down sector carries the mixing (0.12) — the minus-partition couplings that order the masses order the mixing; the anarchic cross-channel counterfactual gives `\|V_us\| ~ 0.46` ⟹ large PMNS (#153) and small CKM from ONE framework (#91 quantified at matrix level). **Stiffness audit:** `V_cb` moves < 1% under ±10% coupling shifts (STIFF — the 10% agreement is a sharp falsifiable prediction); `V_us` swings ×0.55–×3.3 under pinhole (SOFT — the small d–s splitting amplifies sensitivity ~×8; its ×2.0 deficit sits in the soft direction). **CP:** the locked baseline has phase = 0 ⟹ `J = 0` exactly vs observed 3×10⁻⁵ — the quark phase is the open item, constrained to reproduce J without disturbing the fixed `\|V\|`. Budget untouched (zero inputs consumed). Open: the quark CP phase; partition couplings ↔ #152 mouth machinery (`ckm_intra_channel_probe`, PR #155) |
| **Quark CP phase calibration** on the locked Hamiltonian | **Calibrated without disturbing \|V\| or masses; J-ceiling deficit = the #155 soft direction EXACTLY; triangle shape (β = 22°) = the Hopf-phase acceptance test** | Solves the constrained problem #155 posed. **Extension:** the global phase knob is unusable (it enters transport as `cos(phase·dk)` — phase 0.5 collapses `\|V_us\|` ×20; why the mass calibration locked it at 0 and the #155 baseline had J = 0 structurally); CP rides the v3 §4 partition-mixing element with the Hopf-placeholder `φ_q(k) = φ·k`: `H(ε,φ) = H_locked − ε Σ_k e^{iφk}\|k,+⟩⟨k,−\| + h.c.` (locked blocks exactly intact). **Scaling derived:** `J ∝ ε^1.9` (quadratic — one insertion per sector side), sinusoidal sign-changing φ-dependence; shifts O(ε²). **The ceiling identity:** `J ≤ \|V_us·V_cb·V_ub\|` ⟹ predicted ceiling 8.64e-6 vs observed 3.47e-5 — ratio **0.249 = 0.498 × 0.902 × 0.555 exactly** (the per-element #155 soft-direction ratios): the J shortfall IS the V_us/V_ub soft direction, not an independent CP failure; observed CP near-maximal (0.887); consistency lock — when the soft directions land, the ceiling rises to 3.5e-5. **The calibration:** targeting the observed phase content on the predicted `\|V\|` (J_target 7.67e-6): `ε* = 0.0528, φ* = 0.80` — `\|V_cb\|` shifts −0.0% (stiff prediction untouched), `\|V_us\|/\|V_ub\|` −4% (inside the soft direction), masses ≤ 0.5% — **the locked structure survives**; calibrated sin δ = 0.967, near-maximal like the data. **The sharp edge:** the placeholder phase reproduces the AREA (J) but squashes the db-triangle — (β, γ) ≈ (0°, 180°) vs observed (22.2°, 65.9°): **β = 22° is the quantitative acceptance test for the true Hopf-connection φ_q(k)** (the v3 §4 TODO). Consumed: ONE input (the CP phase content — the flavor puzzle’s CP entry made explicit). Open: the Hopf phase (β target); the soft `\|V\|` directions (J-ceiling target) (`quark_cp_phase_calibration_probe`, PR #156) |
| **Flavor-sector synthesis capstone** (#149–#156 assembled) | **Masses + both mixing matrices + CP in both sectors; net +1 input, −1 modelling knob; FIVE falsifiable targets** | The #150/#131-convention capstone for the flavor arc: one keystone re-verified from every member — #149 inversion `ε₃/ε₂ = 1.435`; #151 observed r₃₂ at the ~77th pct; #152 saddle (max rule holds, geo fails ×9, light re-check); #153 angles natural; #154 `m_ββ ≈ 3.2 meV` with the EXACT φ_ℓ invariance; #155 `V_cb = 0.0377` from the mass-locked blocks; #156 ceiling identity exact + calibration point reproduced. **THE CARD (10 rows):** normal ordering (derived) · m₁ ≈ 0.04–0.07 meV (predicted) · `Σm_ν ≈ 58.8 meV` (falsifiable) · ε_n spread (derived, β retired) · three PMNS angles (anarchy-natural 62/56/27th pct) · lepton Dirac CP + Majorana phases (generic) · `m_ββ` 3.2 meV 68% [1.5, 5.9] (falsifiable) · CKM (out-of-sample, zero inputs, V_cb/V_ts stiff 10%) · quark CP (calibrated; β = 22° Hopf test). **Bookkeeping:** eight probes, ONE input consumed (the quark CP phase content) and ONE knob RETIRED (the #151 β) — the modelled-assumption count DECREASED while the sector was assembled. **Mechanism map (the #134 structure, matrixed):** bounce (ν: channel-dominant anarchy) · winding (charged: hierarchy-protected e-row, one permitted μ–τ rotation) · shell (quarks: Z₂ partition alignment) — one geometry, every asymmetry traced. **Targets:** Σm_ν discriminator; m_ββ > 10 meV falsifies; β = 22°; the J ceiling → 3.5e-5; V_cb stiff. The card added to `docs/THESIS.md`. Residual locus: the anarchic draw, one CP content, soft V_us/V_ub, the O_geom e-row (`flavor_sector_synthesis_probe`, PR #157) |
| **Hopf-connection derivation of the quark CP phase** | **One parameter (J-fit only) predicts the FULL triangle to ~1°; pure π/k₅ reproduces all five CP observables uncalibrated (candidate); #156 corrected** | Runs the #156 acceptance test — and relocates the mechanism. **CORRECTION to #156:** the partition-mixing CP was an artifact — the charged-current CKM is non-unitary at ~16% (u–d near-degeneracy leakage), the quartet invariants disagree ×1000 (`J₁₂ ~ 7.7e-6` vs `J_db ~ −3e-9`), the UNITARIZED core has `J ≈ 0` for EVERY φ_q(k) form (linear/k²/Casimir), and the #156 ε* implies a first-row unitarity deficit ~×40 over bounds — doubly excluded (the ceiling identity and J = 0 baseline stand; the calibration interpretation is superseded; new bound: partition mixing `ε ≲ 0.004`). **The relocation:** the locked coupling `−t·e^{−α·dk}·cos(phase·dk)` is the REAL PART of the Hopf transport factor `e^{iφ·dk}`; the two Z₂ partition classes traverse the fiber with OPPOSITE orientation (the #63 C-swap flips c₁) ⟹ `(H±) ∝ e^{±iφ_h·dk}` — exactly unitary V (7e-16), quartet-consistent J (2e-18), baseline = φ_h = 0. **One parameter, the full triangle:** φ_h* = 0.611 calibrated to J ALONE ⟹ `(β, γ, α) = (22.9°, 65.8°, 91.3°)` vs observed `(22.2°, 65.9°, 91.9°)`, sin δ = 0.905; masses shift 0.09%, V_cb untouched, V_us moves TOWARD data (0.112 → 0.123) — the #156 β = 22° acceptance test PASSED. **The π/k₅ candidate (flagged per #107/#108):** the calibration sits 2.7% from `π/k₅ = 0.6283` (the χ = 0 fiber holonomy π over the k₅ winding quanta); the PURE uncalibrated value gives J at 0.969 of target, `(β, γ, α) = (22.8°, 63.5°, 93.8°)`, and **sin δ = 0.888 vs observed 0.887** — five CP observables, zero parameters; candidate pending an independent transport derivation (the #152 path). Budget: no input consumed; the #156 input conditionally RETURNED if π/k₅ confirms (`hopf_transport_cp_phase_probe`, PR #158) |
| **Explicit Hopf-transport derivation of φ_h = π/k₅** | **DERIVED — transport exact to 1e-15; all alternative sector counts excluded; the #156 input RETURNED: quark CP is a calibration-free five-observable prediction** | Supplies the derivation #158 flagged (the #152 modelled→derived path). **The chain:** the rate ½ = `A_φ(0)` (the spin-½ factor of the connection — DERIVED); the sign ± = the #63 C-swap (opposite fiber orientation per Z₂ partition class; explicit opposite transport conjugates — VERIFIED); the winding content `dk = max(k,k′)` = the SAME `winding_mode=max` rule the MASS calibration locked (phase and magnitude rules share their dk — independent corroboration); the arc = one winding-sector `2π/k₅` (capacity k₅ = 5, #73/#126) ⟹ `phase = ±dk·(½)·(2π/k₅) = ±dk·π/k₅`. **Explicit path-ordered transport:** full circuit k = 1 → exactly π (the spinor sign flip — the module’s own consistency anchor); sector arc → `k·π/k₅` exact to 1e-15 for k = 1, 3, 5. **The exclusion scan:** π/3 (generation count) kills J outright (−0.02); π/4 misses γ by 24°; π/6 misses γ by 12° and J by 14%; 2π/5 flips the CP sign; π/10 misses γ by 49° — **π/k₅ is the unique survivor of six principled candidates** (the anti-numerology discipline in its positive mode). **The derived prediction (no calibration):** J at 0.969 of target, `(β, γ, α) = (22.8°, 63.5°, 93.8°)` vs `(22.2°, 65.9°, 91.9°)`, **sin δ = 0.888 vs 0.887**, masses 0.09%, V_cb untouched. **Budget:** the #156 input RETURNED — the flavor card’s last open row closes (quark CP: derived); net flavor-arc bookkeeping #149–#159: inputs +0, modelling knobs −1; THESIS postscript added. Open: the hop arc from explicit shell wavefunctions (the final mile); the soft V_us direction (`phi_h_transport_derivation_probe`, PR #159) |
| **The final geometric mile**: the Hopf sector arc (Weyl algebra) + the pinhole refinement | **Arc = the Weyl commutator quantum (machine-exact); the soft V_us direction reduces to ONE angle (γ); the J-ceiling lock VERIFIED** | Closes the two items flagged at #159. **Part A:** the capacity-k₅ winding space carries the canonical clock–shift pair with `UVU†V† = e^{2πi/k₅}·1` EXACT — the fiber discretizes into k₅ sites `θ_n = 2πn/k₅` and the shell hop is the shift whose minimal step is one site: **the sector arc 2π/k₅ is the Weyl quantum** (algebra, not identification; no radial profile enters) ⟹ `φ_h = (½)·(2π/k₅) = π/k₅` algebraic end-to-end (the #159 caveat removed). **Part B exclusions:** the pinhole single-knob lands V_us but breaks m_s −22.5%; the exact transport rescale (dk = 3 element ×2, dk = 5 fixed) self-defeats via level repulsion (V_us only 0.133, m_s +50%; invariant `sin2θ = 2\|H_ds\|/Δλ` verified). **The mass-preserving joint solution** (eigenvector rotations at fixed eigenvalues, masses to 1e-15; `(δθ_u, δθ_d) = (−5.2°, +9.9°)`) lands SEVEN of eight observables: V_us ×1.00, β = 22.2° exact-fit, **J ×1.05**, V_ub ×1.10, V_td ×1.19, V_cb/V_ts ×0.90 — with **γ = 104° vs 65.9° the single remaining misfit**. **The J-ceiling lock VERIFIED:** #156/#158 predicted the ceiling rises to 3.5e-5 when the soft elements land — at the refined point: ceiling 0.99 of observed, J ×1.05 (a prediction about a then-nonexistent state, checked in that state, passed). **Re-lock targets tabulated** (`H_ds ×1.84` + %-level diagonal compensation) for the next v3+CP joint lock. Residual after the final mile: the γ angle + the knob-level re-lock; no new inputs (#150 budget unchanged) (`final_mile_sector_arc_pinhole_probe`, PR #160) |
| **The γ misfit resolved** — the full flavor-CP dataset realized | **ALL NINE observables land ≤ 1% (four of them PREDICTED) at exactly preserved masses and the derived φ_h = π/k₅; the quark flavor-CP sector closes** | Attacks the #160 residual. **Diagnosis:** γ lives in the ub corner, coupled to the (1,3)/(2,3) rotation planes the #160 two-plane family never used; the full mass-preserving family is **SO(3)×SO(3)** (6 Euler angles at exactly fixed eigenvalues) vs 5 data constraints ⟹ a one-parameter solution manifold exists — the misfit was a restricted-family artifact. **The solution** (residual 0.005): the five constrained observables land sub-percent (`V_us = 0.2256, V_cb = 0.0419, V_ub = 0.00368, β = 22.3°, γ = 65.9°`) and the four UNCONSTRAINED observables are predicted and land — `V_td ×1.01, V_ts ×1.00, J ×1.00, α = 91.8° (obs 91.9°), sin δ = 0.889 (obs 0.887)`. Masses preserved to 1e-14. **The physical branch:** the manifold’s up-dominant end is excluded (the 5768 eigenvalue amplifies sub-degree rotations into ×−181 elements); the down-dominant branch reaches the same residual with O(1–2) targets — up block `H₊[12] ×1.287` (others exactly unchanged), down block `×1.832/×1.996/×1.111`, angles ≤ 6.1° — **the complete targets for the knob-level v3+CP re-lock**. **What closes:** the #157 card’s last quantitative misfit — a single target state realizes masses + the full CKM + the full triangle + J + sin δ at the derived phase, zero new inputs. Remaining: the knob-level re-lock (engineering, targets complete) + the lepton anarchic draw. Honesty: V_td/V_ts/α land partly via unitarity — the nontrivial content is EXISTENCE at fixed masses and fixed derived phase; branch choice documented (`gamma_misfit_resolution_probe`, PR #161) |
| **Flavor phase addendum** — the Hopf CP derivation and the full CKM realization, consolidated | **Every #156→#161 keystone re-verified in one run; the thesis carries the complete account; the #157 card’s quark-CP row → DERIVED** | The consolidating addendum (the #131/#150 convention) for the CP arc: the #156 correction (#158: 16% non-unitarity, quartet-inconsistent J, unitarized core `J ≈ 0`, first-row unitarity ×40 — partition mixing excluded, `ε ≲ 0.004`) → the relocation to the Hopf-fiber transport (`(H±) ∝ e^{±iφ_h·dk}`, the #63 orientation) → the derivation `φ_h = π/k₅` (#159: connection ½ × C-swap × mass-locked dk × sector arc; #160: the arc = the Weyl commutator quantum) → the realization (#160/#161: the mass-preserving family lands all NINE flavor-CP observables at ≤ 1%). **Keystones re-run together:** #158 unitarity 7e-16 / quartet 2e-18 / excluded-route `J ≈ 0`; #159 sector transport `k·π/k₅` exact; #160 Weyl commutator exact; #161 re-solved — residual 0.005, masses 1e-14, predicted `V_td ×1.01, V_ts ×1.00, J ×1.00, α = 91.8°, sin δ = 0.889`. **Bookkeeping:** thirteen probes #149–#161, net ZERO inputs (the #156 input consumed then returned) and ONE knob retired. **Thesis edits:** the interim postscript replaced by the full addendum subsection (chain table + realization + bookkeeping); the #157 card’s quark-CP row updated to derived. Remaining: the knob-level v3+CP re-lock (targets complete); the lepton anarchic draw (`flavor_phase_addendum_probe`, PR #162) |
| **The v3+CP joint re-lock** — the v4 candidate lock | **Targets realized: minimal-law NO-GO (exact, partition-asymmetric); three new targeted couplings + one retune; nine observables ≤ 1%; net surplus +2; CP at ZERO parameters** | The flagged #161 successor. **The no-go (exact):** the v3 off-diagonal law enforces partition-symmetric transport (`H₊[12] = H₋[12]`, machine-exact at the lock) and the `dk = max` degeneracy — the targets break both (partition split ratio **1.424**; minus d–b enhancement **×1.996**) while the up block, where data permits, keeps the law EXACTLY (5e-6). The breaking pattern is the partition asymmetry on the minus block’s d-row — where #155 located the mixing, and where the v3 lock already carries targeted couplings (χ_k3, η_35^−). **The v4 lock:** element level — twelve tabulated numbers + the derived phases reproduce the v3 masses EXACTLY (1e-15) and all nine flavor-CP observables at ≤ 1% (the first complete flavor state in one parameter set); structured level — the v3 law + `η_12^+ = −0.102, η_12^− = −0.295, η_13^− = −0.267` (new) + `η_35^−: 5.0 → 5.586` (retune, +11.7%) + diagonal retunes inside the existing diagonal law (the extension CONTINUES the lock’s own targeted-coupling pattern). **Counting:** +3 parameters for +5 independent observables (the other four follow from unitarity + the derived phase) — net predictive surplus **+2**; the CP sector costs ZERO parameters (`φ_h = π/k₅` derived); masses inherited; the #150 budget unchanged. **Migration staged** (four steps: QuarkParams fields, complexified transport with the derived default, the lock update, regression re-baseline) — the library is untouched here (`v4_relock_realization_probe`, PR #163) |
| **The σ_z-readout capstone** — the geometric object behind σ_z-sensitive readout, identified and classified | **THE MEASUREMENT ARC CLOSES WITH AN EXACT THREE-WAY FACTORIZATION. The meaning is K's (#238): the derived fiber U(1) dial is GENERATED by σ_z and the derived σ_z detector commutes with it (the abelian trap, CHSH = 2); sensitivity needs the unique Δk = 2 mixer (L: pure σ_x at 1/√6; M: Z₈-conserving with a winding-2 carrier) — leaving open only the carrier's geometry. THE OBJECT IS AN ELLIPTIC MOUTH: a quadrupole deformation of the mouth fiber circle — ⟨+1|H|−1⟩ = (2−√2)·ε·e^{iφ₀} EXACTLY on the committed N_χ = 8 fiber (linear in ellipticity, phase = orientation; multipole selection rule exact — only m = 2 connects the qubit). THE PARITY HOME: odd modes have the EXACT #202 neck node (2e-10, #221 re-read), even modes finite neck amplitude (contrast 1e7) — the carrier lives in the even sector the lepton ladder is parity-EXCLUDED from: apparatus and matter in complementary sectors of ONE field. THE HONEST NEGATIVE: no committed source populates the even sector coherently — M's n̄ ≈ 16 condensate (re-read: 1→2.29, 16→2.76) is APPARATUS, the substrate's Stern–Gerlach magnet. THE PHASE IS THE DIAL: qubit restriction = 0.408248(cosφ σ_x − sinφ σ_y), |c| constant 3e-12, angle = −φ exact, Tsirelson at EVERY orientation (4e-15); [H, K_total] = 0 exact with total charge frozen 3e-16 and field coherence nonzero — the #208 no-go untouched. THE CLASSIFICATION: operator = DERIVED GEOMETRY; occupation = ADDITIONAL DETECTOR STRUCTURE; orientation = the IMPORTED MEASUREMENT FREEDOM reduced to ONE U(1), physically located at the mouth** | The capstone of the H–M measurement arc (deliverable `docs/sigma_z_readout_capstone.md`): the '#238 one non-derived ingredient' factors exactly as operator ⊗ occupation ⊗ orientation — what geometry gives, what the apparatus adds, and what freedom remains are three labeled, machine-anchored lines instead of one imported axiom. Scope: K/L/M lattice level; first order in ε; populating the condensate = named successor. Suite 278 passed (`sigma_z_readout_capstone_probe`, PR #241) |
| **The physical gate set** — the qubit clock, the gravitational Kerr, and the twist switch on the committed network | **#232'S STRUCTURAL PLACEHOLDER EXECUTED AND BOTH ARCS WELDED: every gate rate measured on the committed #224 network. THE CLOCK: the doublet beat IS the one-qubit X rotation (Δω = 1.714e-3, period 1833; ledger re-derived at machine agreement; basins 0.977). THE SWITCH: one twisted link (#227's Z₂) freezes the beat 4.9e-11 (the committed row re-derived EXACTLY) — a single topological bit toggles GATE ↔ MEMORY mode (leakage time 3.8e13) — and the twisted basins shift the dressing integral < 0.2%: the twist kills the linear (leaky) channel and SPARES the Kerr channel — an error switch, not a gate hazard. THE GRAVITATIONAL KERR: the committed nonlinearity is the EKG backreaction (the #223 δμ = cA² law, re-read: exactly quadratic); per one quantum δμ_q = 5.80e-2, measured dω/dμ = −0.947 ⟹ χ = 5.49e-2, t_CZ = 57.2 — THE CZ IS 32× FASTER THAN THE CLOCK: GRAVITY IS THE ENTANGLING RESOURCE. HONESTY: δμ_q/μ = 0.64 — one quantum dresses the throat at O(1); χ is leading-order in a STRONG coupling (the self-consistent dressed gate = the successor). THE SELECTION RULE: same-qubit Kerr terms vanish IDENTICALLY on the dual-rail sector (Fock, machine zero) — the encoding auto-silences self-phases; the calibrated CZ under frozen linear coupling is exact (1e-16, zero leakage); residual = 0.13 rad coherent ZZ per CZ, compensable. BUDGET: GHZ₃ = 4700 model units, twist-frozen leakage error 1.6e-20. THE WELD: seven constants from FOUR ledgers spanning BOTH arcs (#224, #227, #223, #232) re-derived at machine agreement** | The #232 successor (deliverable `docs/physical_gate_set.md`): the gate set is physical — clock = the network's own tunneling, switch = one topological bit, entangler = gravitational cross-phase modulation, calibrated per quantum on committed geometry; the two development arcs meet in one machine without retuning. Scope: model units; leading-order χ (strong coupling is the finding); shared-throat two-qubit geometry = construction work. Suite 278 passed (`physical_gate_set_probe`, PR #234) |
| **Tensor-product emergence** — the 2^N state space from one field on the mouth-doublet network | **THE DECISIVE ONTOLOGY QUESTION ANSWERED CONSTRUCTIVELY, WITH TWO PROVEN NO-GOS: the qubit is the MOUTH DOUBLET (which-mouth dual rail; the #224 pair re-read from its ledger), one scalar field supplies 2N modes, and the N-excitation sector with one quantum per doublet has dimension EXACTLY 2^N (N = 1–5) — the exponential space is the N-QUANTUM COHERENCE OF ONE FIELD ON PHYSICAL SPACE. NO-GO 1 (counting): a classical field's 4N real dof lose to the 2·2^N−2 pure manifold at N = 3, exponentially after (#206's LHV theorem re-read as the N = 2 face). THE TENSOR PRODUCT DERIVED: decoupled doublets evolve as the EXACT Kronecker product of their one-body unitaries (1e-15, zero leakage) — LOCALITY IS THE TENSOR PRODUCT. NO-GO 2 (the measured KLM boundary): across the bridge family + 200 random linear-optics evolutions, every sector-unitary map has ZERO entangling power; every entangling map rides O(1) Hong–Ou–Mandel leakage — a linear field on ANY topology entangles only probabilistically. THE COMMITTED QUARTIC COMPLETES THE GATE SET: the #210/#222 |φ|⁴ term as a bridge-overlap cross-Kerr gives CZ EXACTLY at χt = π (1e-16, zero leakage, Bell entropy 1.000000) — the soliton nonlinearity IS the entangling resource. UNIVERSALITY: GHZ_N circuits at fidelity 1−1e-12 for N = 2–5; Mermin = 4 LIVE at N = 3 (ledger #208: 4.0). THE ONTOLOGY: the pilot-wave configuration-space object is DERIVED — one field's N-quantum sector read through N doublet frames; the quantum computer's workspace is physically carried coherence, never parallel worlds** | The general-N successor to #206/#207/#208 (deliverable `docs/tensor_product_emergence.md`): the 2^N structure emerges from one quantized field plus topology — tensor products from locality, entanglement from the committed quartic, universality to N = 5 — and the no-gos prove no classical (counting/Bell) and no linear-deterministic (KLM) shortcut exists. Scope: graph-level network; χ structural; bosonic scalar; quantization imported (structure emerges, QM not derived from classical fields); measurement = #209/#219. Suite 278 passed (`tensor_product_emergence_probe`, PR #232) |
| **Radion stabilization from the primordial EM-capped throat** — V_eff(φ), α at the minimum, and the radion mass | **THE ONE DYNAMICAL PROBLEM #225 LEFT, EXECUTED: THE HOPF CHARGE STABILIZES THE RADION. The committed throat is a DYON: the #55 cap threads it with fixed electric flux (A = αℏc/2 re-read from the committed ledger at machine zero — the α-dependent holdout continues) and the #58 topology puts one Hopf charge on every throat (Σc₁ = 0). Under the #225 dilaton coupling e^{−√3φ}F² the two flux energies carry OPPOSITE radion charges (quadrature through ε = e^{−√3φ}): fixed CHARGE → e^{+√3φ}, fixed FLUX → e^{−√3φ} (Dirac ratio 1/(4α²) exact), plus the #222 primordial frame factor e^{aφ}: V_eff = U_el·e^{7φ/2√3} + U_mag·e^{−5φ/2√3} — ONE EXPONENT OF EACH SIGN: THE MINIMUM EXISTS (numeric = closed form 1e-8), with the coefficient-free identity m_φ² = (35/12)V_min (numeric + symbolic) and E″ = 7U_el* exactly. THE NO-GO: without the Hopf charge every exponent is positive — runaway to decompactification (α → 0): the #58 topology saves the vacuum. THE LANDING: α*² = 5κ/28 — Dirac point (κ=1) gives α* = 0.4226, ORDER ONE, 58× observed (mechanism derived, modulus- and anchor-free); observed α needs κ = 2.98e-4 (verified; guardrail: no match, nearest α/k₅² 2.1% — rejected; cohesion tilt < 0.3%): the open problem sharpens from ρ* to κ. THE MASS: E″ = (7/4)α^{3/2}m_P = 1.33e16 GeV per throat at observed α (fiber anchor, GUT-scale — heavy radion, no fifth-force tension); 13 keV on the #55 Compton anchor; m_φ² = 16πn·E″/m_P². THE ARC: the stabilizer row is EXACTLY the (1,1,0,0,0) rank-audit slot #225 reserved (∇α annihilated 2e-16); ρ = 3.08 Dirac / 23.41 observed; #222's ×210 exclusion re-read — NOTHING ELSE MOVES** | The #225 successor (deliverable `docs/radion_stabilization.md`): the radion effective potential of the dyonic EM-capped primordial throat stabilizes the modulus with closed-form α and a coefficient-free mass identity; the observed coupling is quantified into a single small asymmetry κ with no numerological shortcut. Scope: κ not derived; single-throat Born–Oppenheimer; tree-level; anchors reported separately per B4. Suite 278 passed (`radion_stabilization_probe`, PR #226) |
| **The absolute-coupling capstone (FINAL; revised)** — canonical Hopf–KK normalization, geometric α, the Einstein-frame radion, the α-dependent holdout | **THE CLOSING QUESTION ANSWERED HONESTLY: A MODULUS REMAINS. The canonical chain determines α_k = 4k²(l_P/R_f)² = 4k²/ρ² EXACTLY — ONE continuous modulus ρ = R_f/l_P (the radion). MACHINE-CHECKED: the fiber KK tower (FD = closed form 1e-11) WELDS to the #193 Berger spectra with the Wu–Yang half-charge q = k/2; THE EXPLICIT GEOMETRY MAP: the Hopf fiber of the committed Berger metric has proper length 2πλR_u by quadrature — R_f = λR_u, round point R_f = R_u — the KK-monopole base is the HALF-radius sphere, and with the committed R_u = r_h = 1 the EM cap converts the unit: 1 model unit = 23.41 l_P; charge = fiber winding (flow slope 2(k+η)/R_f² to 3e-4; EXACT flux periodicity 4e-11; adiabatic-ramp force check 3e-4). THE FULL EINSTEIN-FRAME REDUCTION WITH THE RADION (symbolic 5×5 Ricci): the 4D curvature terms carry e^{(2a+b)φ} — Einstein frame ⟺ b = −2a; kinetic coefficient −6a² (canonical −½ at a = 1/(2√3)); gauge kinetic −(R₀³/2)e^{3bφ} = the KK dilaton e^{−√3φ}F²; V_tree = 0 EXACTLY; and the dilaton exponent EQUALS the geometric-law exponent 6a = √3: α(φ) = α(0)e^{√3φ} — THE α LAW IS THE DILATON COUPLING. THE MODULUS-RANK AUDIT: five log-knobs, committed constraints (#222 weld; EM cap) — before the cap rank 1 and ∇α SURVIVES in the flat space (α undetermined); after the cap rank 2 and ∇α is annihilated to 2e-16: the cap fixes EXACTLY the one α-coupled direction, all three remaining flats α-DECOUPLED (rank-proven no-retuning). THE STABILIZER AND THE LEPTONS: ρ* = 2/√α = 23.4125; #165 guardrail: no closure-constant match (nearest e^π 1.2% — rejected); the implied KK scale (√α/2)m_P = 5.2e17 GeV sits 5e18 above m_μ — leptons are NOT fiber KK modes and need not be: COMPATIBLE BY DECOUPLING, with the #197 winding-ladder slope 1/λ = 1/R_f (4e-5) and β_lepton/L_fiber = k₅² exact. THE α-DEPENDENT HOLDOUT (replacing the ρ-independent check): twelve committed ledger constants that GENUINELY carry α — linear, quartic (α⁴ transit protection), inverse-quartic (α⁻⁴ frozen period) — ALL re-derive at machine zero; inverted, ONE COMMON α to 4e-16 relative. THE ARC CLOSES WITH ITS BOOKS BALANCED** | The final capstone, revised pre-merge (deliverable `docs/absolute_coupling_capstone.md`): the modulus is now a canonical 4D field on a flat tree-level potential whose dilaton coupling IS the α law; the fiber is explicitly the committed geometry; the rank audit makes no-retuning a theorem about constraint geometry; the lepton sector is compatible by decoupling; and the holdout is genuinely α-dependent. The Bessel-root digits misprint (2.29948 → 2.29991) is corrected. **What is derived is derived, what is selected is named, and nothing was retuned along the way.** Scope: tree-level (α = the #184-protected boundary invariant); term-isolating test metrics for the reduction; l_P via G₄ (B4 discipline); ℏ still not derived. Runs in CI. Suite 278 passed (`absolute_coupling_capstone_probe`, PR #225) |
| **Mouth-exchange dynamics** — P_other-mouth(t) and its laws on the two-throat network | **THE FIVE REQUESTED EXCHANGE OBSERVABLES, ON THE NETWORK'S GENUINE TWO-LEVEL SYSTEM. AN HONEST TOPOLOGICAL FINDING FIRST: a single two-mouth bridge ring has ONE connected exterior cavity — no mouth-localized doublet exists on it (#223's even/odd splitting is the neck-BC sensitivity, the correct transit-coupling measure, not a beatable splitting); the genuine mouth-to-mouth system is the TWO-THROAT NETWORK (two #221 interior throats on the shared S³ exterior; segment-commensurate grid for the EXACT A↔B swap symmetry — without it the tiny coupling loses to grid detuning). Working doublet (r_s = 0.3, ω ≈ 2.19): basins localized 0.977. P_OTHER-MOUTH(t): the EXACT identity P_B(t) = P_B(0) + [P_max − P_B(0)]sin²(Δωt/2) (machine zero) + direct leapfrog evolution: A depletes to 1e-4, max at t/T = 0.998, half-period at P_max/2. MAXIMUM TRANSFERRED PROBABILITY = the localization weight 0.977 (= 4a²b², limited only by the dressing tail). EXCHANGE PERIOD π/Δω = 1833 (evolution 1829). THE COUPLING LAW: Δω(r_s) exponents 3.3→3.7 → the #223 (ωr_s)⁴ limit — at the primordial anchor the period scales as α⁻⁴: MOUTH-TO-MOUTH TRANSFER IS DYNAMICALLY FROZEN for the physical electron. SURVIVAL WITH REALISTIC ASYMMETRY (clock-rate bias = the MTY aging proxy): the Rabi identity P_max·Δω² = Δω₀² to 1%, the Lorentzian P_max = Δω₀²/(Δω₀²+δ²) pointwise to 1%, survival > 98.6% at a 1.5% bias — LOCALIZATION BY ASYMMETRY closes the exact-degeneracy loophole. THE COMPLEX WIDTH: compact network Γ = 0 EXACTLY (real symmetric operator — the #218 persistence); with an open channel (lead beyond throat B) the QNM doublet has Γ = 4.0e-4/4.2e-4 (outgoing-BC complex Newton, lead-independent < 1%), STRONG COUPLING J > Γ/4 (ratio ~8: the exchange survives under the decay envelope), Γ_pair ≈ Γ_direct/2 (width shared by hybridization)** | Executes the requested observables (deliverable `docs/mouth_exchange_dynamics.md`): P_other-mouth(t), maximum transferred probability, exchange period π/Δω, survival with realistic asymmetry, complex decay width with open channels. **What lands:** the two-level reduction is exact for the doublet and confirmed by direct wave evolution; asymmetry completes the #223 transit protection dynamically (no two physical throats are identical → the dressed eigenhistory stays home); the width statement is sharp (zero on the compact network, computed and lead-independent when a channel opens, in the regime where exchange survives). **Scope:** classical norm fraction as probability (ℏ nowhere); interior-channel reading supplies the bound basins; bias asymmetry (length is grid-quantized); model lead for the open channel; frozen background. Suite 274 passed (`mouth_exchange_dynamics_probe`, PR #224) |
| **The bridge dressing and the two-mouth port** — a closed-form universal, the π/2 step, transit protection | **BOTH REQUESTED STEPS ON THE #222-FORCED PRIMORDIAL BACKGROUND. THE GEOMETRY: the ultrastatic MTY bridge (the #202 t=0 slice ρ = √(r_s²+σ²), unit lapse, traversable) with the exact ρ³-measure reduction u = ρ^{3/2}φ — neck height (ℓ(ℓ+2)+3/2)/r_s² closed-form, far tail = the #215 form to 1e-3, the #202 φ ∝ σ near-neck law to 0.4% — glued into the S³ exterior ring (u = R sinχψ exact). THE PORT: the two-mouth spectrum solved; same-parity spacing = the #217 closure comb π/L to 2%. A NEW UNIVERSAL: the #202 matching radius on the true bridge gives X = 2.2995 — invariant across modes (4 digits), r_s (anchor limit), curvature, seam — WITH THE CLOSED FORM z* = root of z·J₁(z) = 3·J₂(z) (the 5D two-mouth analog of the quarter wave). THE INTERIOR-DEPTH FAMILY UNIFIES #221 AND #223: any regulated tortoise channel past the quarter wave steps X to π/2 EXACTLY (1e-4) — and X ≥ π/2 ACROSS THE WHOLE FROZEN CLASS: the required conv-B 1.5089 sits 3.9% below the class infimum — THE RESIDUAL IS STRUCTURAL (m_e/m_μ ≤ 2α/π = 0.004646 vs observed 0.004836; the anchor or S₁ carries a ~4% correction — the named successor derives it; #165 forbids matching without derivation). TRANSIT PROTECTION: the even/odd splitting through the neck obeys Δω ∝ (ωr_s)^4.1/4.0 (predicted 2ν = 4) — at the primordial anchor the electron mode's non-local network exchange is O(α⁴) ~ 1e-12: THE DRESSED SOLITON SITS STILL while the carrier waves (#213–#221) transit. THE PERTURBATIVE DRESSING: δμ (ρ³-measure, first-order) EXACTLY quadratic in amplitude; throat-local share 0.2% of the cloud — the neck undisturbed, μ the primordial datum; the cloud energy = the particle's exterior-mass contribution; 10%-window A ≈ 0.014** | Executes the two requested steps (deliverable `docs/bridge_dressing_network.md`): the soliton as light perturbative dressing of the fixed-μ bridge, and the multi-mouth port mapping the scale-welded soliton under global non-local transits. **What lands:** the audit's ρ³-measure item executed on the TRUE geometry; #221's 1D-cavity π/2 and the bridge's Bessel z* are one interior-depth family; the class-wide bound X ≥ π/2 makes the 4% residual structural and localizes it to the anchor/S₁ — a sharp successor target; the network answer to "how does the welded soliton behave under transits" is quantitative: protected in identity as (ωr_s)⁴, coupled only through the carrier comb. **Scope:** ultrastatic vs horizon readings bracket X (adjudication = successor); transparent mouth seam (robustness checked); first-order back-reaction; classical, frozen background; α imported only for the anchor. Suite 274 passed (`bridge_dressing_network_probe`, PR #223) |
| **The coupled 5D EKG weld** — the scales lock, and the lock refutes self-sourcing | **THE #221 AUDIT'S SUCCESSOR, EXECUTED — THE ONE PR THAT UNFREEZES THE GEOMETRY: the soliton core's energy density directly sources the local 5D Tangherlini metric (static spherical Einstein–Klein–Gordon, shooting; the #210 solver generalized to n = D−2 — n = 2 IS #210: Kaup M_max = 0.6327 reproduced; tt-Einstein identity pointwise 3e-6; dx/xmax-converged 2.5e-5/1e-5). THE WELD EXISTS: R*, R_RMS, and r_s = √μ∞ emerge from ONE solve in ONE unit — convention invariance EXACT (the (KF, φ)→(KF/4, 2φ) rescale shifts every geometric observable by machine ZERO; the audit's free radial rescale is forbidden by the field equations). THE ISRAEL JUNCTION VERIFIED (option B): [K^θ] ~ 1e-16, [K^t] ~ 5e-7 at the core boundary, each tracking the scalar tail (no shell); the exterior reads the interior mass to 1e-6. THE 5D CRITICAL-MASS MARGINALITY MEASURED: dilute μ → μ_crit = 7.695 (difference ratio 2.0), ω → 1, size FREE at fixed mass (X ∝ φ_c^{−1/2}, slope −0.503) — the 5D Newtonian 1/r² marginality, the same marginality as #221's tails. THE CAVITY FROM THE SAME METRIC (vacuum-validated vs the #215 closed form to 2e-6): R* tracks the THROAT (R*/r_s = 0.78–1.00) while R*/R_RMS falls 0.39 → 0.17 — R* and R_RMS emerge together and DECOUPLE: the cavity belongs to the gravitational core, not the cloud. THE CONFRONTATION: r_s·ω PINNED INTO [1.53, 2.774] over the entire family (spiral included; q-channel in-band) vs the EM-cap anchor r_s·ω = α — EXCLUDED ×210; jointly σ/r_s = 206.8 forces X = 574 vs required 1.51 (×380): A 5D SCALAR CANNOT BE MUCH LIGHTER THAN THE THROAT ITS OWN FIELD CREATES — self-sourcing refuted at the coupled level, #210's PRIMORDIAL relocation FORCED by the weld, the #221 cavity is the primordial throat's** | Executes the requested action modification (option A) with the Israel junction as the verification layer (option B) — both options in one system (deliverable `docs/coupled_5d_ekg_weld.md`). **What the weld decides:** the audit asked whether any equation ties the soliton length unit to the bulk unit — the coupled solve ties them (exactly, with the convention dropping out) and the tied values FORCE the primordial reading: the throat mass is a geometric datum, not the electron's self-field. **The program's framing survives its own stress test:** this PR unfreezes the geometry precisely to test the frozen-background posit, and the coupled system confirms it as the only consistent reading of the anchor. **Successor (sharply defined):** the soliton as perturbative dressing of the fixed-μ primordial bridge — this exact solver on the two-mouth topology. **Scope:** boson-star ansatz for #180's real ψ-φ-q (potential-class robustness: q-channel spot checks land in the same r_s·ω band, per #210's compactness argument); one-mouth topology (the refutation is of SELF-SOURCING); ground states (spiral entered, not traced — tracing shrinks the band from inside, never to α); test-wave cavity; classical; α imported for the confrontation only. Suite 274 passed (`coupled_5d_ekg_weld_probe`, PR #222) |
| **The derived O(1) lepton mass coefficient** — the quarter-wave invariant, with the landing conditional on the soliton↔cavity weld | **THE QUARTER-WAVE INVARIANT IS DERIVED (the m_e/m_μ landing CONDITIONAL on the soliton↔cavity weld — see the independent audit at the end of this cell): X = σ_mode·ω* IS A RATIO OF TWO PROPERTIES OF ONE GEOMETRIC OBJECT — the eigenhistory whose frequency is the mass (λ̄_C = 1/ω*) and whose #202 matching radius is σ_mode. THE MODE IS FORCED: the arc (#217/#218) puts eigenhistories on the throat's interior resonance comb (ring interior = #202's bridge coordinate, neck = cross-cap) and #202's Pin-twisted boundary condition forces the electron (k = 1) transit mode ODD — the node at the neck reproduced to 1e-12 with the φ ∝ σ near-neck law (drift 2.0%); the even mode touches the neck (the k = 0 uncharged channel — #202's dichotomy realized). THE HARD-WALL THEOREM: the odd fundamental's antinode is the quarter wave — X = ω·σ_match = π/2 EXACTLY, cavity-length independent (the analytic reason the coefficient is O(1); the definitional family in closed form, verified 5e-9). THE MEASUREMENT (real Tangherlini barriers): X_match = 1.5783 (reference D = 12); regulator scan depth 8–48 band [1.579, 1.638] (4%, all within 6% of π/2); exterior spread 0.035 (the matching radius robust exactly where RMS is hybridization-contaminated); dx-converged 7e-4. THE EIGENHISTORY CARRIES IT: the #220 Gauss–Newton orbit on the odd mode (residual 3e-13) with the source EXACTLY decoupled by parity (q* = 3e-17 — the charged transit mode is source-transparent at the crossing), complete monodromy unit-circle to 2e-14, orbit-averaged X unchanged to 6e-7, exactly linear (amplitude-independent, 1e-13). THE INVARIANCE SUITE: X unchanged under energy budget ×16 (independent GN solves — X and T identical to 6e-15/5e-15), source coupling g = 0→4× (the SAME orbit periodic to 3e-13, monodromy at 4g unit-circle 3e-14), mode amplitude, cavity depth (4% regulator band), and BRANCH CHOICE — three odd branches give X_match = 1.620/1.575/1.571 → π/2 (spread 0.049) while the RMS definition grows ×3.3: branch invariance singles out the #202 matching radius as the physical definition. THE ALTERNATIVE BRANCHES RUN: the second odd branch's complete GN orbit (residual 6e-13, source still exactly decoupled, monodromy unit-circle 2e-13) carries X_match = 1.57070 — π/2 TO 1e-4; the even (k = 0 channel) branch's full source-COUPLED orbit (q* = −6.2e-3 nonzero — the coupling contrast; residual 4e-13, unit-circle 3e-12) lands the parity-excluded conv-A alternate at 0.004891 = 101.1% of observed — numerically closer than the primary, excluded PURELY by #202's parity theorem (the sharpest falsification target for the 5D successor). THE CONFRONTATION (inputs: geometry + α via #210's r_s = α·λ̄_C; m_e used ONLY for comparison): convention A (required 0.6467) EXCLUDED ×2.4; convention B (required 1.5089) SELECTED — **m_e/m_μ = α/X = (2α/π)(1 + throat shift) = 0.004455–0.004624 vs observed 0.004836: 92–96% OF THE OBSERVED RATIO WITH ZERO FITTED NUMBERS**; the #210 neck aspect c = ln(1/α) + O(1) now has the O(1) computed (ln X ≈ ln π/2 = 0.45). THE INDEPENDENT IDENTIFIABILITY AUDIT (`lepton_o1_identifiability_audit_probe`, 8/8) NARROWS THE CLAIM, ADOPTED: the quarter-wave invariant is UNCONDITIONAL, but no existing equation welds the soliton length unit (#202/#203's σ_mode is the SOLITON IR scale) to the cavity unit — the soliton length family spans ×4 under the cavity frequency, a radial-unit rescale moves X linearly — so the m_e/m_μ landing is CONDITIONAL on the eigenhistory-particle identification (σ_mode = the cavity antinode); the audit's successor contract (coupled Pin-Dirac/soliton bridge state, ρ³ measure, action-derived 3D→5D map, antinode-free overlap functional, coefficient locked before comparison) derives or refutes the weld** | Executes the #210 register item "derive the O(1) (σ_mode/λ̄_C)" on the requested #220 background (deliverable `docs/lepton_o1_coefficient.md`). **What changed in kind:** #210 left σ_mode/λ̄_C as a ×2.3 constrained band; it is now a derived value with a 4% regulator band, the #201 convention ambiguity is resolved BY the measurement (B selected, A excluded — not assumed), and m_e/m_μ becomes an outright prediction landing at 92–96%. **The definitional discipline (the #165 guardrail):** parity from #202's theorem (not chosen), length from #202's matching radius (not chosen), the hard-wall π/2 is a theorem (not a matched constant), the even/conv-A alternate lands in-band numerically but is parity-EXCLUDED and recorded as such. **The residual is one-sided and owned:** the throat correction pushes X up from π/2 while the observation sits 4% below it — the 1D transit measure is not the 5D ρ³ bridge measure, and redoing X on the true 5D radial operator is the named successor. **Scope:** winding k not represented (zonal scalar — parity imported); interior depth a scanned regulator (the α-capped tortoise depth a successor refinement); m_μ enters only as the #201 anchor scale; classical, frozen geometry. Suite 274 passed (`lepton_o1_coefficient_probe`, PR #221) |
| **The PDE-ring eigenhistory** — every reduction retired, the field confirms the source | **THE FULL HAMILTONIAN SUCCESSOR TO #219: THE FIELD IS THE PDE — the two-direction wave on the periodic ring (C = 2π + 8, N = 192) with the glued finite-width Tangherlini barriers — AND THE DUFFING (q, p) EVOLVE IN EXPLICIT TIME DOMAIN, H = Σdx[π²/2 + (Du)²/2 + Vu²/2] + p²/2 + ω₀²q²/2 + μq⁴/4 + g·q·u(0) WITH THE INTERACTION TERM IN THE CONSERVED ENERGY. THE TERM IS LOAD-BEARING: under symplectic leapfrog the full H is conserved to a bounded non-secular 2e-5 over ONE HUNDRED periods (O(dt²) ratio 4.0) while the ledger without it fails 99× worse in a SINGLE period. THE PERIODIC ORBIT, LITERALLY AS REQUESTED: unknowns (X(0), T), X ∈ R^386 the complete source-field state; conditions X(T) − X(0) = 0, H[X] − E₀ = 0, p(0) = 0 (the single phase condition removing time-translation degeneracy); Gauss–Newton with the full FD Jacobian converges QUADRATICALLY to 2e-13 — E₀ = 20.834 hit exactly, ω_orbit = 2.669296 (nonlinear pull −3e-3), dt-converged 3e-7. THE ONE-PERIOD MONODROMY OF THE COMPLETE STATE: the tangent leapfrog (exact linearization of the symplectic map) gives the 386×386 monodromy — ALL eigenvalues on the unit circle to 4e-14 (NO parametric instability of ANY field mode against the eigenhistory), det = 1 to 2e-13, the dx-weighted symplectic form preserved to 4e-15, the trivial pair to 9e-11. THE SOURCE PAIR CONFIRMS #219 FROM THE FIELD ITSELF: identified by (q,p)-eigenvector weight (0.36 vs 0.008 next), rotating at the dressed frequency 3.2124 — WITHIN 0.07% OF #219'S REDUCED-MODEL 3.2102 (bare ω₀ = 3.2); the low field modes at their folded angles. THE LEDGER, TIME-RESOLVED: E_int oscillates about −0.0467 (negative — q and u(0) antiphased), equal to the #219 harmonic-balance value g·a·U₀/2 to 0.04%; the full H constant along the orbit to 5e-5. SPATIAL CONVERGENCE: the complete solve repeated at N = 128/192/256/384 (ring, barriers, coupling, E₀ fixed) — ω_orbit, q, u(0), the partition, the interpolated profile, the source pair, and the low Floquet angles ALL O(dx²) with the exact predicted diff ratios (6.400/2.400; measured 6.06–6.66/2.36–2.42), every monodromy unit-circle; Richardson continuum ω_orbit = 2.673464, source frequency 3.2132 — 0.09% from #219. Harmonic balance, single-mode reduction, point-barrier transfer matrices — retired, and the #219 answer survives** | The requested full Hamiltonian successor (deliverable `docs/pde_ring_eigenhistory.md`): explicit time-domain evolution of the Duffing variables coupled to the two-direction field ring, the interaction term in the conserved energy, the one-period monodromy of the complete source–field state, and the periodic-orbit conditions solved with one phase condition. **What the full system buys:** #219's Floquet analysis lived on a reduced 2-dof model with the caveat that higher ring modes were dropped — the complete monodromy now carries EVERY discretized field mode's Floquet pair, and all of them are marginal: the eigenhistory is stable against the whole field, not just its own mode. **The confirmation is nontrivial:** the ring here uses the physical finite-width #215 barrier potentials, not #219's point-scatterer transfer matrices — the source pair matching to 0.07% shows the local source physics is robust to the fuller geometry. **Scope:** the discrete symplectic map at n_s = 2048 (continuum approached O(dt²), verified by refinement); E₀ a specified budget, ℏ not derived; classical, zonal scalar, frozen geometry — the ring is the time-closed loop's covering model, the CTC/network aspects live in #216–#218. Suite 274 passed (`pde_ring_eigenhistory_probe`, PR #220) |
| **The Hamiltonian source eigenhistory** — the imposed law derived, the joint solve | **NOTHING ABOUT THE SOURCE IS IMPOSED ANYMORE: THE MINIMAL CONSERVATIVE SOURCE H = p²/2 + ω₀²q²/2 + (μ/4)q⁴ + g·q·u(0) (side-coupled Duffing, harmonic balance) YIELDS ITS SCATTERING AS CONSEQUENCES — unitary BECAUSE conservative (1e-13), reactive (zero net power, 1e-15), the amplitude-dependent phase carried by D(a) = ω₀²−ω²+(3/4)μa² — AND #218'S IMPOSED LAW IS ITS WEAK-COUPLING LIMIT (ratio 1.000000). THE RING: the reflecting source makes the loop a genuine two-direction ring (real trace, 1e-10); the #217 gap tangency RESOLVED (tr > 2 on [2.7325, 2.7390] — the modes are the split pair at the gap edges); the homogeneous condition defines a NONLINEAR MODE BRANCH ω*(A) — homogeneity alone cannot fix the amplitude. THE JOINT SOLVE U(X)X = X + TOTAL-ENERGY CLOSURE ON THE CORRECTED LEDGER (⟨g·q·u(0)⟩ = g·a·U₀/2 included): 2D Newton → (ω*, U₀*) = (2.732375, 0.914367) — shifted 0.2% by the interaction term — a* = −0.1963; FIXED-POINT RESIDUAL ‖U(X*)X* − X*‖ = 3e-14 (Newton 9e-16/1e-13, slaving 1e-13); energy partition E_field = 11.060 + E_source = 0.171 + E_int = −0.05397 = E₀ to 1e-13; flux constant around the ring (2e-15); 10⁴-pass drift 5e-11. THE FULL HAMILTONIAN STABILITY SPECTRUM — the Duffing (q, p) EVOLVED AS INDEPENDENT VARIABLES in the 4×4 variational monodromy of the reduced 2-dof Hamiltonian (ω_r = 2.738858 and g_eff = g·ψ(0) = 0.3203 both DERIVED from the bare ring) about its shooting-refined periodic orbit (residual 4e-14, NNM frequency within 8.8e-5 of the full-ring ω*): the Floquet-trivial pair at 1 (8e-9) + THE SOURCE PAIR 0.4541 ± 0.8910i AT ITS DRESSED FREQUENCY 3.2102 (bare ω₀ = 3.2) — the dof the slaved harmonic-balance map froze — ALL \|λ\| = 1, det M = 1 − 6e-14, symplectic 5e-14, orbit energy drift 6e-11: the conservative eigenhistory is MARGINAL in the full Hamiltonian sense — Novikov-passive, no runaway anywhere in its phase space** | The requested replacement and joint solve (deliverable `docs/hamiltonian_source_eigenhistory.md`): source state and energy in the loop, the homogeneous condition solved together with total-energy closure, fixed-point residual and stability eigenvalues reported. **What the Hamiltonian buys:** #218's theorem survives with its one modeling element removed — the eigenhistory now rests on (ω₀, μ, g, E₀) and geometry alone; the amplitude is fixed by the ENERGY BUDGET (the mode branch parameterized by amplitude, the budget picking the point), sharpening #218's phase-only fixing into the requested energy-closure form. **The marginal spectrum is the right answer:** a conservative periodic orbit cannot be attracting — an all-unit-circle symplectic Floquet spectrum IS Hamiltonian marginality, and the route to attraction is the #209 dissipative registration open. **The winding map vs Floquet distinction:** the slaved harmonic-balance map iterates windings (its spectrum, retained for comparison, answers winding convergence); the Floquet monodromy evolves time over one carrier period — only it sees the source pair. **Numerics honesty:** mouth ports are spline interpolants of the unitarized Tangherlini greybody, re-unitarized pointwise; raw-port systematic 2e-5 at the fixed point; third-harmonic neglect 2e-4. **Scope:** E₀ specified (the #58 nucleation quantum is program-level input); ℏ not derived; monochromatic skeleton (packet eigenhistory = named successor); classical, frozen geometry (`hamiltonian_source_eigenhistory_probe`, PR #219) |
| **The eigenhistory transaction** — the existence theorem | **A FINITE, SOURCE-INCLUSIVE, ENERGY-CONSERVING, SELF-CONSISTENT WORMHOLE TRANSACTION EXISTS. The completed transaction is not a response to a source — it is a HOMOGENEOUS, GLOBALLY CONSTRAINED EIGENHISTORY F = Λ_tot(ω, I)F, existing where the null space of the full transfer system opens, its amplitude fixed by the two closures. THE SOURCE JOINS THE LOOP: a reactive scatterer (\|s\| = 1 to 1e-16, zero net absorbed power) whose phase is pulled by the field — φ_s(I) = βI/(1+I/I_sat), the classical internal-state shift with saturation. ENERGY CLOSURE (\|Λ\| = 1): the throat is lossless EXACTLY on its interior resonance by a unitarity identity (loop factor real positive ⟹ \|t_net\| = T/T = 1: exact tier 1e-16; Tangherlini deficit 2.0e-4 = 10× the solver's port flux error) — the eigenhistory can live ONLY on the resonance comb, forced by energy closure alone. STATE CLOSURE (arg Λ = 0) FIXES THE AMPLITUDE: φ_s(I*) = −χ mod 2π with χ = 2 arg t = 2.8631; the IVT (source range 12 > 2π) guarantees a root — I* = 1.5944 (residual 8e-14), Λ_tot − 1 = 3e-16; σ_min(I − M_tot) = 1e-16 AT the eigenhistory vs 0.14/0.9 detuned — the homogeneous solution exists at exactly one point per branch, and THE BRANCHES ARE DISCRETE: I*₀ = 1.594, I*₁ = 16.90 — a discrete transaction-amplitude spectrum from classical closure. DYNAMICS: energy conserved per pass to 1e-16; 10⁴-pass persistence (drift 1e-16/1e-11); the null eigenvector carries the source's state shift AND the mouths' resonantly enhanced interior — one closed history; neighbors dephase at exactly dφ_s/dI; the physical tier decays at exactly its solver deficit. THE CORRESPONDENCE: #217's driven pole IS the eigenhistory (response 1e15 there vs O(10) detuned) — the marginal Novikov point is populated; weak driving shadows it** | The requested formulation and existence theorem (deliverable `docs/eigenhistory_transaction.md`). **The arc closes conceptually:** #213 derived the propagator from complete histories; #216 gave the advanced half a globally causal mechanism; #217 resummed the loop self-consistently; #218 exhibits the completed transaction itself — not a limit, not a pole of a response function, but a finite self-consistent object with its own fixed amplitude, including the source's internal-state shift and the mouths' excited interiors as parts of one closed history. **Two-tier verification:** the unitary model class (exact, 1e-16) where the losslessness identity is algebraic, and the Tangherlini instance (solver precision, deficit = finesse-amplified port flux error). **Honesty:** ℏ NOT derived — the amplitude scale is the source's nonlinear saturation (β, I_sat) and the discreteness is branch structure, not quantization (the successor question is connecting the eigenhistory scale to the quantum of action); monochromatic skeleton (a localized packet eigenhistory needs the #217 group closure too — named successor); marginal (Novikov-passive) stability — selection by dephasing, with a dissipative registration mechanism (#209 opens) the route to attraction; single-channel reactive source; frozen geometry; classical zonal scalar (`eigenhistory_transaction_probe`, PR #218) |
| **The self-consistent network loop** — G_eff derived before I± | **THE LOOP IS SOLVED, NOT SAMPLED: F = F₀ + ΛF ⟹ G_eff = g/(1 − Λ), WITH EVERYTHING IN Λ DERIVED — AND THREE CORRECTIONS OF #216 FORCED IN THE PROCESS. (1) CLOCK-RATE-CORRECT TRAVERSAL: the throat runs at ω_τ = ω/rate_A; the emergent frequency is ω·rate_B/rate_A — a slow-clock exit mouth REDSHIFTS (correcting an inverted ratio in #216); elastic confirmation REQUIRES rate-matched mouths — state closure derives the equal-rate analysis. (2) THE VALUE-TRANSPORT EIGENVALUE Λ = t_net(ω_τ)·e^{iωD_loop}, validated FROM FIRST PRINCIPLES against a ring spectrum (transfer-matrix modes 2.7429/3.6000 vs the Λ = 1 comb 2.7390/3.6071 — 0.2%; the opposite convention misses by 0.16, EXCLUDED) — correcting #216's closure bookkeeping (offset from the engine convention exactly e^{−iω(dA+dB+τ)}, 1e-16; deform knob and magnitudes unchanged): at time closure the carrier closes on the throat's own scattering phase. (3) WIGNER-DELAY-CORRECT CLOSURE: solving group (D_loop = −τ_W) and carrier (arg Λ = 0) closure TOGETHER — the packet lands ON the crossing (0.006 vs 0.175 uncorrected, phase-aligned; Λ = 0.999998 + 0i), and below the barrier both equations are solved by the throat alone: THE COMPLETED TRANSACTION SITS ON THE INTERIOR FABRY–PÉROT RESONANCE (τ* = 5.1105, residual 2e-12, \|t_net\|² = 0.969) — completed transactions live on the throat's resonance comb, a solved fixed point. THE TRANSFER SYSTEM: matrix resolvent = 1/(1−Λ) to 9e-15; nested interior×winding resummation 1e-15. STATE EVOLUTION: input–output rates match the derived linewidth and storage time (30.7 vs 31.1); the source fixed point converges at rate \|Λ\| EXACTLY (0.23239907 vs 0.23239900) — marginal at completion: a persistent self-consistent history. G_EFF, THEN I±: the completed-transaction line is QUARTIC (group closure makes the loop phase stationary; FWHM 0.0116 vs 0.0118 predicted — anomalously flat: completed transactions are ROBUST to detuning), Lorentzian at carrier-only closure (0.0439 vs 0.0456); PASSIVITY max\|Λ\| = 0.9999993 ≤ 1 — THE NOVIKOV FIXED POINT CANNOT RUN AWAY; O(Λ) = #216's K₁ = I₊ + ΛI₋; the resummed confirmation weight is Λ/(1−Λ) (at completion ~5e5 — the divergence IS the transaction pole); the deficit 1−\|Λ\| → 0 at completion (2e-6 vs 0.11): THE COMPLETED TRANSACTION IS THE #213 COHERENT ε → 0 LIMIT AND THE DEFICIT IS #214'S ABSORBER DAMPING — THREE ε'S, ONE NUMBER** | The decisive #216 successor (deliverables `docs/self_consistent_network_loop.md` + module upgrades in `transaction/network.py`): both mouth barriers, clock-rate-correct traversal, Wigner-delay-correct packet closure, repeated winding resummation, source and mouth state evolution, and the effective Green function derived before comparison with I±. **The module fixes:** `emergent_frequency` inverted ratio corrected (redshift out of the well); `t_AB`/`r_AA`/`loop_expansion`/`traverse_throat` run at the throat proper frequency with the exact clock composition of the exit time; `loop_eigenvalue` (value transport, ring-validated) and `effective_green` added; 16 unit tests. **The honesty ledger:** the #216 convention correction affects closure-phase bookkeeping only — magnitudes, the advanced projection, and the deform knob survive; the engine weld remains true of the engine's own convention; migration staged. **Scope:** continuous ω interpolates the physical tower; the above-barrier carrier tuning uses the mouths' geometric transfer phase (free network parameter) while the below-barrier point is throat-solved; CMT rates a moderate-finesse consistency check; frozen linear network per winding (no back-reaction of stored cavity energy); classical zonal scalar; MTY aging history posited. Suite 274 passed (`self_consistent_network_loop_probe`, PR #217) |
| **The wormhole-network confirmation** — the advanced half gets a mechanism (two-port throat) | **THE ADVANCED HALF OF THE WHEELER TRANSACTION IS AN EVERYWHERE-FUTURE-DIRECTED NETWORK TRAVERSAL — THROUGH A TWO-PORT THROAT WITH REPEATED-LOOP INTERIOR PHYSICS. Waves only propagate forward — but the network reaches the past: one retarded C-wave → future antipode (t = π) → mouth A's interface (#215 greybody t_A) → FORWARD traversal, looping k times between the two interior barrier faces (each loop +2τ_th local, future-directed; echo train at (2k+1)τ_th + Δ_BA damped by \|r_in\|² = 0.667) → mouth B's OWN interface (t_B) → emergence at t < 0 (differential aging, Morris–Thorne–Yurtsever) → forward on S³ → THE ORIGINAL CROSSING (t_return = t_emit to 1e-12). THE TWO-PORT COMPOSITION IS DERIVED AND VALIDATED DIRECTLY: t_net = t_At_B/(1 − r_inAr_inB e^{2iωτ}) matches a from-scratch solve of the glued two-barrier potential to 3e-4 in \|T\|² across a sweep through resonances; the mouth S-matrix measured from BOTH sides (reciprocity 2e-4; unitarity relation r_in = −r̄t/t̄ at 6e-5); the Airy average = the incoherent echo sum (1.5e-4). RESONANT CONFIRMATION (new physics): identical mouths transmit PERFECTLY on the interior resonance comb even deep below the barrier — single-port T = 0.33 → resonant T_net = 0.9996 vs 0.040 off (= T²/(1+R)² exactly); live packet confirms 0.88 ON resonance (stored in the cavity, released over the storage time ≈ 58) vs 0.040 off — 22× contrast: the #215 IR transparency is punctured by the throat's interior modes. THE CENTRAL IDENTITY: clock continuity FORCES U_BA = e^{iωΔ_BA} (1e-12) — the clock offset IS the transfer phase — so the absorption→return kernel is K = Λ(ω)e^{−iωD′} with D′ < 0 and \|Λ\| = \|t_net\|: THE RETARDED RULE ANALYTICALLY CONTINUED TO A NEGATIVE EXTERIOR INTERVAL = THE ADVANCED KERNEL. PHASE CLOSURE = #213'S COHERENCE CONDITION: a discrete frequency comb (spacing 2π/\|Δ_BA + τ_W\| to 0.06%); the engine's retro_phase_match accepts the loop at the comb (1.000) and rejects it detuned (0.03). THE POLE: at closure I₊ + ΛI₋ = the covariant pole; DETUNING THE NETWORK IS THE #213 DEFORM TEST (φ = ωδ to 1e-15) — the refutation edge has a geometric knob; off the interior comb the exact split (1+\|t_net\|)/2 coherent + (1−\|t_net\|)/2 deficit** | The requested model plus the review items — the second interface and the repeated-loop physics (deliverables `geometrodynamics/transaction/network.py` + `docs/wormhole_network_confirmation.md`): `advanced_confirm_amplitude` replaced by an explicit network traversal. **The data structures:** `NetworkMouth` (psi, link_id, clock_rate, clock_offset, orientation, transfer_phase), `MouthPort` (t, r_out, r_in — the full one-sided S-matrix of each interface), `NetworkThroat` (mouth_A, mouth_B, τ_th, port_A, port_B) with Δ_BA, U_BA(ω) AND the composite t_AB(ω) DERIVED, not free. **Time-machine honesty:** global hyperbolicity broken by construction — the closure condition is exactly the Novikov self-consistency fixed point, and transactions are its solutions; no chronology-protection dynamics addressed (frozen background per program framing). **The network is posited, not solved:** mouth pair from the #58/#200 nucleation channel, Δ_BA from MTY differential aging; the 1D interior channel is model-level — the direct glued-barrier solve validates the composition on exactly that model; the true 5D interior is #168/#200 territory. **Mode selectivity:** the closure comb ∩ the interior resonance comb — demonstrated, not spectrum-matched. **Migration staged:** advanced_confirm_amplitude retained; the weld shown. 14 unit tests; suite 272 passed (`wormhole_network_confirm_probe`, PR #216) |
| **The throat's spectral density** — the flat bank retired | **THE ABSORBER'S SPECTRUM IS THE THROAT'S OWN GREYBODY TRANSMISSION, WITH BOTH WINGS PINNED BY THEOREMS OF THE GEOMETRY. The scattering problem on the 5D Tangherlini mouth (ψ_xx + [ω² − V]ψ = 0, V = f[(ℓ(ℓ+2)+3/4)/r² + (9/4)r_h²/r⁴], ingoing at the horizon, Hankel-matched; flux conservation 2e-5) yields T_ℓ(ω) and: (1) THE IR WING IS THE AREA THEOREM — the universal σ_s(ω→0) = A_horizon = 2π²r_h³ MEASURED (σ/A_h = 1.012 at ωr_h = 0.04, monotone → 1, extrapolated 0.983; T₀ ∝ ω³, slope 3.02): J(ω→0) ∝ A_h·ω³/4π is the horizon AREA, not a parameter — the universe is IR-transparent; (2) THE UV WING IS THE PHOTON SPHERE — r_c = √2·r_h, b_c = 2r_h (machine-checked), T = ½ at the eikonal (ℓ+1)/b_c (ratios 1.04/1.02/1.01), T → 1 above the barrier: UV-BLACK, a perfect absorber for every mode that resolves it; (3) THE HORIZON IS THE CONTINUUM — the tortoise depth diverges as −(r_h/2)ln δ (fitted −0.5005 vs −½ exact): level spacing → 0, recurrence → ∞ — #214's order of limits (N → ∞ first) is a property of the horizon, and classically the throat NEVER REVIVES; (4) THE WELD IS PARAMETER-FREE — cavity quasimodes (complex secant: ingoing at horizon, node at wall) obey γ = T(ω_q)/(2L_cav) at 0.994/1.013/1.080 down the three-mode ladder, spacing π/L_cav to 4%: ε = ω·T(ω)/(2L_cav), FULLY GEOMETRIC. The tower density (r_h/R = 0.1, one antipodal delivery per half-period): 0.0017 (n = 1, rising as n³) → 0.85 (n = 8) → 1.000 (n = 24). AFTER #213/#214/#215 NOTHING ABOUT THE ABSORBER IS CHOSEN — spectrum, address, and ε all read off the frozen bulk** | The named #214 successor (deliverable `docs/throat_spectral_density.md`): the last hand-chosen element of the absorber sector removed. **The solver:** tortoise-coordinate integration (x(r) = r + (r_h/2)ln((r−r_h)/(r+r_h)) analytic) with least-squares Hankel matching over spread points (a 2-point solve is phase-degenerate at ωΔx ≈ nπ) and an ω-scaled far radius (the metric's long-range 1/r² phase otherwise contaminates the UV at 7%); regulator independence 4e-6. **The cavity length subtlety:** L_cav runs wall → barrier peak, not wall → horizon — the sub-barrier region is the drain, not the cavity (with L to the horizon the transit ratio is 1.25; to the barrier peak, 0.994). **The WKB bracket:** the transit ratio drifts 0.994 → 1.08 as ω_q climbs toward the barrier top — the transit picture degrading exactly where WKB says it should, bracketing the law's validity domain. **Physical reading:** with the #210 anchor r_s ~ αλ̄_C, laboratory modes sit far into the IR-transparent wing (everyday fields propagate freely); the throat blackens at its own Compton-edge scale — the same scale the #211 capstone located as the natural domain's edge. **Falsifiers checked:** σ/A_h ↛ 1 (no), crossings off the eikonal (no — 1–4%), widths violating the transit law (no — 0.6–8%), finite tortoise depth (no — log-divergent at the exact −½ coefficient). **Scope:** classical scalar greybody (no Hawking, no spin/tensor); zonal 1D transit reduction (full S³-with-throat matched asymptotics named successor); r_h/R a parameter; horizon interior #168/#200 territory; frozen geometry (`throat_spectral_density_probe`, PR #215) |
| **The dynamical absorber** — S = S_field + S_absorber + S_coupling | **ABSORPTION IS AN OUTCOME OF THE DYNAMICS, THE iε IS A COMPUTED DAMPING RATE, AND THE GEOMETRY ITSELF PICKS THE ABSORBER'S ADDRESS. The absorber promoted from #213's boundary condition to a degree of freedom: a bank of oscillators (the throat's internal continuum — the engine's `MouthState.modes`, given an action for the first time) coupled to the tower through Φ(ψₐ). EXACT ELIMINATION (Gaussian): G_eff(Ω)⁻¹ = ω₀² − Ω² − Σ(Ω) verified against direct inversion of the coupled system (1e-10); stiffness positive-definite (stable absorber, no ghost). THE ε DERIVED: the pole moves below the axis by a computed amount Ω* = ω̃ − iγ/2 — the complex pole, the LIVE time-domain energy decay, and the golden rule γ = πg²ρ/(2ω₀²) agree (0.0499/0.0499/0.0500), scaling EXACTLY as g² (slope 1.9996); the orderings' iε filled in as the physical rate iγ/2 (quadrature 3e-8); the Δ²-pole form agrees to O(γ²); g → 0 recovers the #213 kernel monotonically — THE BOUNDARY CONDITION WAS THE ZERO-COUPLING SHADOW OF THE DYNAMICS. THE ADDRESS ADJUDICATED (antipodal vs distributed): γₙ = π(gYₙ(ψₐ))²ρ/(2ωₙ²) — (a) the ANTIPODE is impedance-matched: Yₙ(π) ∝ n never vanishes and exactly cancels the kinematic suppression ⟹ γₙ = g²/(4π), THE SAME RATE FOR EVERY TOWER MODE (the #166 caustic as impedance matching); (b) a GENERIC point leaks (near-zero \|sin 22\| = 0.009, rates ∝ 1/n² — long-lived remnants); (c) a UNIFORM distribution is FORBIDDEN by the exact selection rule ∫sin(nψ)sinψ dψ = (π/2)δₙ₁ — it couples only to the ground mode. Live three-way adjudication (24 kicked modes + 3000-oscillator bank, exact normal-mode evolution): residual field energy 0.014 (antipode, per-mode uniform, = e^{−γt} predicted 0.0136) vs 0.95 (generic) vs 1.00 (uniform). WHAT ε MEANS: a finite bank revives near T_rec = 2π/dν (linear in N, ratio 1.90) — on the closed bulk finite ε IS Poincaré recurrence, and ε → 0⁺ is the continuum limit of the throat's spectrum taken first** | The decisive #213 successor (deliverable `docs/dynamical_absorber_propagator.md`): everything the transactional derivation imposed as a boundary condition now follows from an equation of motion. **The model:** S_field = the conformal tower ωₙ = n/R; S_absorber = the oscillator bank; S_coupling = −Σgⱼqⱼ·Φ(ψₐ) with Yₙ(ψ) = sin(nψ)/(sinψ√(2π²)) — all classical, all Gaussian (that is what makes the elimination exact rather than perturbative). **The falsifiers checked:** placement-independent decay (no — factors of ~70× between placements), n-dependent antipodal rates (no — per-mode residuals within [0.7, 1.03]× of the flat prediction), γ failing g² scaling or the three-way agreement (no — slope 1.9996, 1% agreement), revivals not scaling with N (no — linear). **The successor named:** derive the absorber's spectral density from the throat geometry itself — the quasinormal spectrum of the Tangherlini mouth replacing the flat bank. **Scope:** Gaussian/bilinear by design; anharmonic registration/irreversibility remains the #209 open; the source point shares the antipode's matching property (the two distinguished zonal points, recorded); classical, per-mode, frozen geometry (`dynamical_absorber_propagator_probe`, PR #214) |
| **The transactional Compton propagator** — the time structure derived from complete bulk histories | **THE ONE IMPORTED ELEMENT OF THE COMPTON AMPLITUDE — THE FEYNMAN iε TIME STRUCTURE — IS NOW A THEOREM ABOUT COMPLETE HISTORIES ON THE FROZEN BULK. The construction: (1) the bulk is STATIC ⟹ t → −t is an isometry ⟹ the elementary field is the time-symmetric Wheeler–Feynman Ḡ = (G_ret+G_adv)/2 (tower identity G_ret(t) = G_adv(−t) at machine zero on ω_ℓ = (ℓ+1)/R); (2) the bulk is CLOSED ⟹ the perfect-absorber hypothesis is a THEOREM of the geometry — every retarded front refocuses at the antipode at t = πR (50% of the \|G\|² mass within 0.3 rad vs 1e-28 at mid-flight) and returns at 2πR (the #166 focusing as absorber condition); (3) complete histories leave no free remnants ⟹ a well-posed 2×2 system (condition number 1.0) with the UNIQUE solution G_F = Ḡ + (i/2ω)cos(ωt) = (i/2ω)e^{−iω\|t\|} — THE FEYNMAN PROPAGATOR IS THE TIME-SYMMETRIC FIELD PLUS THE ABSORBER RESPONSE OF THE CLOSED UNIVERSE. THE TWO COMPLETION ORDERINGS are the offer and confirmation segments of one history: each half-line transform is an OFPT energy denominator — individually regulator-dependent (truncation spread 0.50, never converges) and non-covariant (Δ→−Δ violation ≥ 0.21) — and the coherent sum is the EXACT covariant pole −(ω−iε)/[ω(Δ²−(ω−iε)²)] (7e-15). THE RELATIVE PHASE IS FORCED: the isometry maps offer into confirmation ⟹ φ = 0; the deform test e^{iφ} breaks evenness and the pole form at O(1) for EVERY φ ≠ 0 (0.81/0.59 at 0.3, 2.0/3.9 at π) and fails the Compton denominators by ≥ 7%. The tie-in: the two-ordering sums equal 1/(s−m²) and 1/(u−m²) EXACTLY across lab kinematics (1e-14; the u-channel entirely off-shell — the "virtual" regime is the coherent pair, no on-shell intermediary): the denominators #35/#211 assumed are now derived** | The stronger successor the capstone pointed to: from RECONSTRUCTING the QED amplitude on the fixed geometry to DERIVING its transactional structure from complete bulk histories (deliverable `docs/compton_transactional_propagator.md`). **The factorization after this PR:** amplitude = denominators (this PR: complete histories) × numerators/vertices (#37–#44, #46) × pole/magnitude/tensor geometry (#35/#45/#46). **The engine connection:** the repo's `retro_phase_match` (offer + confirm phase closure) peaks EXACTLY at the coherent point — the transactional engine has been enforcing the derived phase all along; its second π branch is the #48 Möbius exchange sign, a distinct discrete sector, recorded honestly. **The regulator is physical:** ε is the inverse absorption time of the closed universe — the finite-T truncation of a single ordering never converges (spread stays at the full oscillation amplitude) while the ε-damped ordering converges to the closed form (7e-12); the frequency-splitting identity G̃_F = θ(Ω)G̃_ret + θ(−Ω)G̃_adv holds with linear ε-convergence (ratio 10.0/decade). **Falsifiers checked:** a second complete-absorption solution (none — unique), a static bulk failing time reversal (machine zero), a closed bulk failing return (full return measured), any φ ≠ 0 surviving covariance (none). **Scope:** per-mode, tree-level, free propagator; spinor numerators/vertices not rederived; ε → 0 still an idealization; frozen geometry, no backreaction (`compton_transactional_propagator_probe`, PR #213) |
| **The Compton-edge capstone** — Release II | **THE #204–#210 ARC CLOSES GREEN, WITH ONE NEW PARAMETER-FREE LAW. The Compton-edge law: adding the mass term to #202's bridge equation `(ρ³φ′)′ = [k(k+2)ρ + m²ρ³]φ` deforms the exact suppression law UNIVERSALLY in x = σ_mode/λ̄_C alone (m·r_s = 0.02 vs 0.05 coincide to 2.6e-4; massless limit re-derives φ = σ to machine zero — #202's c₀(1) = 0). The law is ε₁ = (r_s/σ_mode)·D(x) with D = 0.966 (conv A) / 0.831 (conv B); the sensitivity S(x) = \|d ln ε₁/d ln σ\| = 1 EXACTLY below the Compton scale (the #202 naturalness, its validity domain now located) and grows beyond — so NATURALNESS CAPS THE #210 O(1): S ≤ 2 confines x ≤ 2.58, the ladder's own worst-tolerated 4.48 confines x ≤ 5.58, and BOTH convention anchors sit inside (S = 1.07, 1.36). The deformation-corrected self-consistent band tightens the #210 bracket to σ_mode/λ̄_C ∈ [0.626, 1.307] — STILL BRACKETING 1: the electron/muon ratio is a fine-structure phenomenon with its one O(1) confined to a compact, derived window. Plus the whole arc re-verified green in one run** | The closing PR of the arc (deliverable `docs/compton_edge_capstone.md`); does what #200 did for the previous arc — the chain green in one run, the register updated — and carries the new law rather than a dialed number (closing a research arc with a dialed O(1) would undo what twenty PRs removed). **Ledger 1 — commitment chain (#204/#205):** a kicked RETARDED run with real potentials conserves the norm to 6e-14 and closes continuity at 2e-3; the classical mean-field channel keeps a product state at entropy 0 vs the quantized pairwise operator's 0.026; SN phase at the interferometry record 4.8e-17 rad, BMV witness 0.79 rad where BAM predicts zero (the two near-term nulls armed). **Ledger 2 — topology chain (#206–#208):** ψ_eff = (I⊗T)\|Φ⁺⟩ is the singlet (fidelity 1.0, E = −cos(a−b) to 2e-16); the swapping law leaves (1,4) in Ψ_{a+b+c} (fidelity 1.0, separable holonomy mixture 2e-17); enumerated LHV bounds CHSH = 2 / Mermin = 2 vs the bridge's 2√2 and the Y-junction's Mermin = 4 with exactly-empty pairwise marginals (negativity 0). **Ledger 3 — measurement + mass (#209/#210):** the k-odd dispersion identity to 2e-16, a live Stern–Gerlach Born check within 0.0024 of cos²β, a live Kaup point M = 0.632, the α bands σ_mode/λ̄_C ∈ [0.647, 1.509] and m_e/m_μ = (3/7…1)·α ∈ [0.00313, 0.00730] bracketing their targets. **The register after the arc:** derive σ_mode/λ̄_C within the Compton-edge window; connect the #55–#58 R* to the bulk mass μ = r_s²; the 5D pants nucleation; W-class reachability; registration/irreversibility; standing negative kept — the cosmological constant (#165). **Falsification card:** SN-null + BMV-null (near-term), m_e/m_μ = (3/7…1)·α at ×1.5 with the O(1) confined, the neutrino cards unchanged. **Scope:** full-GR no-signaling a completion argument; CM-wave emergence derived only in sourcing + entangled structure; pinch/pants as topological content; the α anchor constrained not derived; equilibrium a hypothesis with a mechanism. Release II: the deepest machinery is theorems, measurements, and bounded windows — not imports (`compton_edge_capstone_probe`, PR #211) |
| **The strong-field core solve** — the collapse reading refuted, the anchor relocated | **THE #203 NUMBER IS MEASURED, AND THE EDGE FIRED ON THE PRE-REGISTERED BRANCH. The static spherical GR family of the committed ψ–q structure (Kaup benchmark: M_max = 0.6327 vs 0.633 — 0.2%; ω/m = 0.848 at criticality) shows: (1) THE #179 RUNAWAY IS A GENUINE GR INSTABILITY — the adiabatic order channel (V_int = −(gσ²−a₀)²/4λ) puts the maximum mass AT the ordering onset with the ordered branch losing mass monotonically (0.542 → 0.295): the strong-field endpoint IS collapse; (2) THE MEASUREMENT — at criticality σ_mode/r_s = 2.45–5.83 (RMS) / 5.42–12.72 (R99) across Kaup / the committed q-channel / a repulsive control: O(few–10), a factor ≥ 7 short of the required 88.6 (structural — the Buchdahl-type bound); (3) THE VERDICT — the collapse contraction is O(1) (the weak-field band 4.6–6.7 overlaps the critical band), NOT 13–45: **THE MOUTH-PAIRING MECHANISM IS REFUTED AS THE QUANTITATIVE ORIGIN OF m_e IN ITS COLLAPSE READING** (the smallness mechanism survives, per #203's own pre-registration); (4) THE RELOCATION — the throat is PRIMORDIAL (#168's Tangherlini Killing horizon) at the #55–#58 EM-capped radius r_s ~ α·λ̄_C: the required hierarchy becomes σ_mode/λ̄_C = 88.6α = 0.647 (A) / 206.8α = 1.509 (B) — **THE CONVENTION BAND BRACKETS 1**: σ_mode = λ̄_C, r_s = r_e, NO NEW NUMBER; the neck aspect c = ln(1/α) + O(1) (4.484 = 4.920 − 0.436); m_e/m_μ = (3/7…1)·α = [0.00313, 0.00730] brackets the observed 0.00484; the #203 window met by this anchor ([13.2, 44.6] ∈ [13, 45])** | Executes the #203 register target (deliverable `docs/strong_field_core_solve.md`). **Why static + turning point suffices:** a collapsing core passes through the maximum-mass configuration of its static family; the turning-point criterion locates the instability; the critical scales bracket the endpoint's — full dynamical NR sharpens the ringdown, not the adjudication, and cannot move the answer past the Buchdahl-type bound (5D criticality differs by O(1) — unable to rescue ×15). **The solver:** shooting with the ground state bracketed by the 0→1 node transition in the trial frequency; observables cut at field decay; vectorized frequency batches. **The families** (σ_c grid 0.1–0.45): Kaup M: 0.475 → 0.633 (max at 0.4) → 0.629; q-channel M: 0.475 → 0.542 (max AT the onset 0.15) → 0.295 (unstable branch); repulsive M_max = 0.757. **The refutation is constructive:** the committed dynamics really does collapse its cores — collapse just cannot make them light: at any collapse endpoint m_e stays over-predicted by an order of magnitude. **The α-anchor arithmetic** (machine-checked): 88.6α = 0.6466, 206.8α = 1.5091, band ∋ 1; c_needed = ln 88.6 = 4.484 vs ln(1/α) = 4.920 (O(1) residual −0.436 = ln 0.647); contraction from the anchor = (r_q/R*)·needed = [13.2, 44.6] using #203's locked scales — inside the pre-registered [13, 45]. **Constrained, not derived** (the O(1) band ×2.3 wide); the mass ladder's unknown is now the SAME α the #184 thread protects. **The register:** 'do full 5D NR' → two sharper items: derive σ_mode/λ̄_C (the O(1)); connect the #55–#58 R* to the bulk mass μ = r_s². The mass ladder and the α thread are one thread (`strong_field_core_solve_probe`, PR #210) |
| **The measurement sector** — pointer outcomes for the entangled sector | **THE ENTANGLED-SECTOR THREAD CLOSES OPERATIONALLY: CHSH = 2.824 READ OUT OF THE POSITIONS OF CLASSICAL BEABLES, with no imported quantum rule anywhere in the chain. The three links, derived/measured: (1) THE POINTER COUPLING EXISTS IN THE COMMITTED STRUCTURE — winding couples minimally to the fiber connection (winding = charge, #42–#44: the KK gauge coupling); the per-channel potential is exact on the lattice and its k-ODD part is an identity (1e-16): a connection-gradient region exerts opposite forces on opposite windings — THE WINDING STERN–GERLACH IS CHARGE MEASUREMENT BY DEFLECTION (live from the raw dispersion: k = +1 crosses to +28, k = −1 turned back to −23); (2) FIBER-INTEGRATED EQUIVARIANCE — #198's theorem extends verbatim to the multichannel wave: continuity residual 1e-4 live, Born ensemble at noise through branch separation; (3) BORN FOR INTERNAL STATES, MEASURED — P(+) = cos²β to ≤ 0.003 across the sweep (branches 7σ apart): the internal Born rule is spatial equivariance routed through the committed coupling, not a postulate; POINTER PERMANENCE: the empty branch's influence dies with the branch overlap (1e-4 → 7e-9) — effective collapse, from geometry, without collapse. THE CLOSING: #206-singlet + local rotations + SG at both wings + dBB beables: E(0,0) = −1.0000 exact, four correlators within 0.012 of −cos(a−b), CHSH = 2.824, marginals setting-independent (operational no-signaling). THE SPATIAL SECTOR: positional EPR from conservation at nucleation (Duan–Simon 0.53 < 2); the #205 guiding-without-gravitating split realized exactly in measurement** | Closes the last open of #206–#208 (deliverable `docs/measurement_sector.md`). **The chain:** internal (fiber/spin) state → spatial pointer branches → position beables → Born. **The device, derived:** V_k(x) = −2t_χcos(2πk/N − θ(x)) exactly (χ-translation invariance at each x ⟹ channels decouple); k-odd part −2t_χ sin(2πk/N) sinθ (identity); the spin/frame analog is the fiber-geometry (Berger-squash, #192/#197) gradient. **The measurement theorem:** for orthogonal internal channels with real per-channel potentials, the fiber-integrated ρ, J close continuity exactly — a throat transported by the fiber-integrated v = J/ρ is equivariant THROUGH the measurement interaction: P(k) = \|a_k\|² follows from #198, measured to 0.003 at 20 000 throats. **The operational Bell:** exact accelerated-Gaussian branches, Heun transport of 20 000 pairs on (x₁,x₂), outcomes = sign(x): sanity E(0,0) = −1.0000; CHSH = 2.824 vs Tsirelson 2.8284; Alice's marginal shift under Bob's setting change ≤ noise — the #204/#206 marginal theorems at the pointer level. **The spatial sector:** the #58/#200 C-conjugate pair is born at a local event with conserved momentum → anticorrelated momenta → the EPR spatial state: Var(x₁−x₂)+Var(p₁+p₂) = 0.53 < 2 (any separable state ≥ 2) — structure from symmetry, widths as inputs; the empty pointer branch guides until separation and never gravitates (#205 conditional sourcing; back-action ~1e-17). **Scope:** SG = co-moving connection-gradient transit; registration/irreversibility beyond branch separation (amplification, radiative decoherence — the #204 dissipative machinery) not modeled; equilibrium hypothesis; nucleation widths underived. **Program-wide remaining opens:** the 5D pants nucleation, W-class reachability, the strong-field NR core contraction (`measurement_sector_probe`, PR #209) |
| **The GHZ sector** — multipartite entanglement is bridge valence | **A THREE-MOUTH BRIDGE MAKES GHZ, AND CHARGED GHZ IS SUPERSELECTION-FORBIDDEN — BOTH DERIVED. The no-go (enumerated): a flux-label (winding = charge) Y-junction has ZERO conserving channels in the ±1 doublet (all sums ±1, ±3), and GHZ components straddle distinct total-charge sectors — QM's charge superselection, derived; the pairwise sector never felt it ((k,−k) is zero-sum): THE MULTIPARTITE SECTOR IS WHERE CHARGE AND SPIN PART WAYS — GHZ lives only in the transported-frame (spin/Pin) label. The junction (live, on the #206/#207 lattice): one bulk fiber read by three mouths — exactly symmetric distribution (1e-16), per-leg deck phases ±i as transport demands, one shared variable, leg-cut collapsing to the #206 pair with the composed law (valence 2 = the junction's special case). THE STATE: GHZ at fidelity 1.00000 with the MULTIPARTITE HOLONOMY LAW φ = −π(s_B+s_C)/2 (phases 0/π/π/2 exact at (2,2)/(0,2)/(2,1)) — one composition rule φ = −πΣs/2 at every valence; 3-tangle τ = 1.0000. MERMIN = 4.0000 (the algebraic maximum) against the exactly-enumerated local bound 2 — from PAIRWISE-EMPTY marginals (negativities 0, pairwise CHSH exactly 2): all correlation in the triple, none in any pair — the exact opposite of #207's matching. BRIDGE VALENCE IS THE ENTANGLEMENT CLASS** | Closes the first named open of #207 (deliverable `docs/multi_mouth_bridge_ghz.md`). **The object:** the trousers (pair-of-pants) nucleation — three boundary mouths sharing one bulk fiber; lattice: 192-ring × 8-fiber + a junction fiber coupled to mouths A (identity leg) and B, C (conjugating legs with holonomies s_B, s_C). **The embedding generalizes:** W\|k⟩ = \|k⟩⊗T\|k⟩⊗…⊗T\|k⟩ is an isometry for n = 2 and 3 alike (1e-16) — the N-party tensor structure is the N-frame description of one bulk object (#206's lemma, at every valence). **The witness chain:** GHZ-fidelity witness > ½ (genuine tripartite); τ = C²(A\|BC) − C²_AB − C²_AC = 1.0000 with pairwise concurrences 0; Mermin optimization over x–y settings reaches 4.0000 while every reduced pair has negativity 0 and CHSH = 2.000000 exactly. **The valence ledger:** valence 2 → Bell pairs, monogamy = matching, swapping (the a+b+c law); valence 3 → GHZ, pairwise-empty, Mermin 4; holonomies select the state within each class; charge superselection prunes the hypergraph (flux labels admit only zero-sum channels — the pairwise sector; frame labels populate every valence). **Scope:** the junction fiber stands in for the 5D trousers cobordism (Sorkin class + nucleation rate vs the #58/#200 pair channel unsolved — whether tri-mouth nucleation is dynamically REALIZED is the physical successor question); W-class reachability via networks + surgery named open; the spatial/measurement sector is now the ONLY standing open of the entangled-sector thread (`multi_mouth_bridge_ghz_probe`, PR #208) |
| **Entanglement swapping is bridge surgery** — linking never-co-nucleated throats | **LOCAL PAIR ANNIHILATION RELINKS BRIDGES, AND THE COMPOSITION LAW IS THE QUANTUM SWAPPING LAW. QM side machine-checked (16-dim): projecting mouths (2,3) of Ψ_a ⊗ Ψ_b onto Ψ_c leaves (1,4) in Ψ_{a+b+c} (fidelity 1−1e-12, probabilities exactly ¼). Bulk side measured (4-mouth lattice): φ₁₄ = φ_a + φ_b + φ_c with the BELL OUTCOME SUPPLIED BY THE JUNCTION HOLONOMY OF THE PINCH — the four-outcome orbit reproduced to ≤0.024 rad (fidelities ≥ 0.9998), end-to-end conjugation purity 0.999832 (winding superselection: outcomes confined to the charge-zero sector). THE EVENT, LIVE: two populated disjoint bridges; the proximity pair (2,3) annihilates LOCALLY (fibers glue, wells vanish); the distant (1,4) response to Alice's phase rises ×4650 (3e-6 → 1.3e-2), lands in exactly the swapped channel (power ratio ~290), carries her phase linearly (D(α) ∝ e^{iα}−1: 1.99 vs 1.99, 0.100 vs 0.100), and the middle drains 0.50 → 0.29 → 0.20 to propagating waves (control ~0.54) — THE PAIR RETURNS TO WAVE FRONTS, THE LINKAGE PERSISTS. The swapped pair saturates Tsirelson (CHSH ≥ 2.8280) between mouths that NEVER SHARED A GLUING; the unconditioned mixture over the four holonomies is SEPARABLE (negativity 0 — no outcome record, no usable correlation: classical communication + no-signaling in one identity); MONOGAMY IS THE MATCHING (partner negativity 1/2 exact, non-partner 0 exact); repeater chains compose associatively — pairwise entanglement distribution over the throat network is MECHANISM-COMPLETE** | Opens the dynamical half of #206's register consequence (deliverable `docs/bridge_surgery_entanglement_swapping.md`). **The mechanism** (user-level statement made quantitative): throats interact — a proximity pair pinches off in pair annihilation (the #58/#200 topology-change channel run in reverse), relinking the bridges of two previously unrelated pairs just as if they shared a nucleation. **The lattice**: 192-ring × 8-fiber, four mouths, bridges (1↔2), (3↔4) with holonomies s₁₂ = s₃₄ = 2; junction gluing (m₂,χ) ↔ (m₃,(s_j−χ) mod 8) — a LOCAL operation (8 sites); two-stage evolution (junction off/wells on → junction on/wells off = the pinch). **The orbit**: s_j = 0/1/2/3 → φ₁₄ = 0.0000/1.5874/3.1650/4.7288 vs predicted 0/π/2/π/3π/2 — the four Bell outcomes of swapping = the four fiber holonomies at which the pair can pinch. **The before-control**: response 3e-6 = the ballistic Lieb–Robinson tail (mouths 1,4 are 64+ sites apart; no path exists pre-surgery). **Monogamy** (three-bridge six-mouth state, machine): a mouth is maximally entangled with exactly its bridge partner and in a product state with everything else — the matching IS the monogamy of maximal entanglement; surgery is the only rewiring move (two pairs in → one pair + radiation out). **Scope:** the pinch modeled by its topological content (fiber gluing + well removal; the full 5D pinch-off with its Sorkin causal-continuity structure not solved — the #200 cobordism machinery is where it lives); the holonomy outcome distribution rests on equilibrium (¼-each at dBB grade, #198/#204); internal sector; **named opens: multipartite (≥3-mouth) junctions (the GHZ sector) + the spatial/measurement sector** (`bridge_surgery_swapping_probe`, PR #207) |
| **Configuration-space emergence** — entanglement is bridge topology | **CLASSICAL ER=EPR, MADE QUANTITATIVE. The sharp target stated and machine-checked: a single classical field on 3-space with local dynamics is an LHV model — all 16 deterministic local strategies enumerated, max CHSH = 2 EXACTLY (Bell); the quantum sector needs 2√2 (Tsirelson). So the entangled sector must come from the one nonlocal element BAM owns: THE BRIDGE (two mouths, one object through the bulk — #168/#169). Derived: the tensor-product structure is the TWO-FRAME DESCRIPTION OF ONE BULK OBJECT (the shared-fiber embedding W\|k⟩ = \|k⟩_A ⊗ T\|k⟩_B is an isometry — this is where configuration space comes from); with the repo's derived transport T = iσ_y (T² = −1, the SAME Pin⁻ sign that forces Fermi statistics) the symmetric bridge gives ψ_eff = (I⊗T)\|Φ⁺⟩ = THE SINGLET — `bell.bulk_identity`'s postulated state, now DERIVED (fidelity 1.000000000000; E(a,b) = −cos(a−b) matching to 2e-16). THE LAW: Schmidt weights = bridge-mode amplitudes (concurrence tracks sin 2β to 2e-3); CHSH = 2√(1+C²) EXACTLY — 2 with the bridge cut, 2√2 with the symmetric π-holonomy bridge. Entanglement is bridge topology, measured** | Attacks #198 condition 2 at full strength (deliverable `docs/configuration_space_emergence.md`). **The identification, live on a lattice** (universal field on a 128-ring × 8-fiber, local everywhere except the mouths glued through the bulk as χ → s−χ): (a) CHARGE CONJUGATION — k=+1 into A arrives at B as k=−1, purity 1.000000 (the C-conjugate pair, Σc₁ = 0, #58/#200); (b) PHASE LOCKING — B's inter-channel phase tracks A's with slope 1 exactly (one phase, not two); (c) THE HOLONOMY IS THE BELL-STATE SELECTOR — the handle's π fiber holonomy (s = 2, a topological datum) offsets the locked phase by exactly π: the singlet sign; (d) ONE OBJECT — both mouths read the same shared mode (ratio spread 1e-8); (e) CUT CONTROL — remove the handle: B receives 1e-30; everything is carried by the bridge. **The ER=EPR curve** (β sweep of the bridge preparation): C_ext = {0, 0.500, 0.707, 0.866, 1.000} vs sin 2β exact; S_ent = bridge participation entropy; CHSH = {2.000, 2.236, 2.449, 2.646, 2.828} = 2√(1+C²) (Horodecki, cross-checked by direct settings optimization ≤ 5e-3); singlet fidelity at β = π/4: 1.000000; bridge-cut product state: entropy 0, CHSH = 2.000000000. **The budget:** marginals invariant under all Alice-side unitaries (1e-16) — Bell correlation without a telegraph; operational statistics at dBB grade on equilibrium (#198), whose signal-locality #204 measured; non-traversability honest — the identification is imprinted at nucleation (the mouths bound ONE 2-handle, #200) and conserved topologically; the lattice handle is a stand-in. **The register consequence:** #198 condition 2 SPLITS — the Bell-sector half DISCHARGED (non-factorizable structure has a topological origin: the bridge), the dynamical half (spatial sector, N-body bridge networks/swapping, measurement transport) remains, sharply bounded. **Scope:** internal (fiber/winding = spin/charge) sector — where Bell lives; k = ±1 doublet (scalar reduction of #195/#197); N_χ = 8; equilibrium hypothesis (`configuration_space_emergence_probe`, PR #206) |
| **The SN-phenomenology audit** — the committed Φ[ρ] meets the laboratory | **THE DYNAMICS MAKES THE PICK, AND THE REGISTER GAINS ITS NEAREST-TERM FALSIFIERS. #204 made Φ[ρ] LOAD-BEARING, so SN-type self-gravity is now a COMMITMENT — and the lab predictions hinge on the sourcing subtlety #198 flagged: UNIVERSAL (empty branches gravitate → full Schrödinger–Newton phenomenology) vs CONDITIONAL (mass rides the actual branch → SN nulls). Adjudicated BY THE COMMITTED DYNAMICS: the beamsplitter experiment — the bound throat-soliton transmits/reflects AS A WHOLE (max(R,T) ≥ 0.95 at 7/8 velocities; co-occupation confined to one grid point; deep-sub-threshold leakage ~1e-4) while the linear control co-occupies branches across the ENTIRE sweep (max ≤ 0.82); the lab regime (kinetic/binding ~ 5e-13) sits eleven orders below the sandbox threshold: THE MASS NEVER CO-OCCUPIES INTERFEROMETER ARMS — effective-level sourcing is CONDITIONAL (#198's aside, measured). The branch-attraction channel is REAL and exactly Newtonian (measured/predicted = 1.000–1.005, zero fit) with lab weight f² ~ 0: SN SIGNATURES PREDICTED NULL. The classical Φ CANNOT ENTANGLE (S = 0 machine vs 0.148 for the quantized pairwise comparator at the same coupling): BMV WITNESS PREDICTED STRICTLY NULL where quantized gravity gives 0.79 rad. Existing data exclude NOTHING (SN phase at the Fein record 4.8e-17 rad; mass margin ~1e5 below m* ≈ 3–6e9 amu; ω_SN(Si) ≈ 0.05/s beyond current reach)** | The #204 follow-through (deliverable `docs/sn_phenomenology_audit.md`). **Why decidable:** at the FUNDAMENTAL level sourcing is unambiguous (the committed Poisson/wave equation sources ψ² + μq²); the ambiguity lives in the EFFECTIVE CM description, where it reduces to a dynamical question the sandbox answers: does the mass-carrying field co-occupy branches when a bound soliton meets a beamsplitter in the lab regime (kinetic ≪ binding)? **E1 (the pick):** relaxed 1D throat-soliton (M = 4, μ_c = −10.5, ordered core q_max = 0.88) vs barrier, velocity sweep 0.3–1.2: R = 1.0000/0.9999/0.91 below, T = 0.9998–1.0000 above, transition bracket (0.5, 0.6) — whole-body transport; linear control (binding off): R,T = (0.82, 0.13) → (0.18, 0.64), never clean — the whole-body behavior IS the self-binding; splitting the mass belongs to the relativistic/QFT regime (#200's domain), not interferometry. **E2 (the rate):** hand-prepared 50/50 split contracts at the exact periodic-pair Newtonian rate d̈ = −2πGM(1−2d/L) (ratio 1.000 at t ≤ 3, 1.005 at t = 5) — the channel is real; signatures scale as f²; BAM's f exponentially zero. **E3 (the discriminator):** product two-throat state under mean-field Φ[ρ₁+ρ₂]: entanglement entropy 0 to machine precision (structural: c-number field = sum of local potentials = zero entangling power; conditional variant = stochastic local channel, LOCC); pairwise operator 2πG\|x₁−x₂\| at the SAME coupling: S = 0.148 by t = 4. **E4 (SI):** f=1 SN phase at Fein-2019 (2.7e4 amu, 266 nm, 10 ms) = 4.8e-17 rad; m*(σ) = (ℏ²/2Gσ)^{1/3} = 5.7/3.3/2.6e9 amu at σ = 100 nm/500 nm/1 μm; ω_SN(Si) = 0.0505/s (Yang formula, Δx_zp = 4.9 pm); BMV witness phase 0.79 rad (1e-14 kg, 200 μm, 2.5 s). **The discrimination matrix:** existing interferometry — everyone says fringes (no exclusion); SN-scale tests (1e9–1e10 amu, ω_SN searches) — QM null / SN SIGNAL / **BAM NULL**; BMV witness — quantized gravity ENTANGLES / **BAM STRICTLY NULL**. **The register gains two near-term nulls, both actively hunted:** SN-null (an SN detection refutes the committed sourcing) and BMV-null (an observed witness refutes classical Φ outright) — the first BAM predictions addressable by living experimentalists. **Consistency:** #204 unaffected (it audited the MAXIMAL channel — a fortiori causal); fundamental equivariance untouched (Born measure and source coincide at field level, part ways only effectively); effective-level corrections ~1e-17 rad. **Scope:** 1D structural/regime adjudication (3D repeat a follow-up); the CM pilot-wave emergence (guiding without gravitating) remains the standing open item, narrowed; the committed signature is the COMBINATION SN-null + BMV-null + fringes at every mass; if the throat sector ever requires f > 0, these bounds propagate back (`sn_phenomenology_audit_probe`, PR #205) |
| **The nonlinear no-signaling audit** — the Gisin edge, faced by construction | **THE EDGE FIRED WHERE IT SHOULD — AND THE THEORY SURVIVES IT. BAM's pilot equation is GENUINELY nonlinear (superposition defect `O(1)` vs `2e-14` for the stripped linear control), so the Gisin/Polchinski threat applies — and indeed a local unitary kick at A reaches B (separation 30) at `5e6×` the kinematic floor, INSTANTANEOUSLY, in the Newtonian model. But the channel is EXACTLY the gravitational field — O(G) (ratio 2.14 under G→G/2), and clamping Φ to the no-kick history collapses it BACK to the floor (suppression `6e6`; q is local and carries nothing) — and RETARDATION CONFINES IT TO THE CONE: with the causal `□Φ = 4πGρ`, the response is MACHINE-FLOOR quiet (`~1e-15`) outside the light cone, the front scales as d/c (`c·t_front = 25.2/26.4/27.2` at c = 8/12/16 vs geometric 25), and c→∞ recovers Newton. The completion is FREE: the retarded Φ is still REAL, so #198 equivariance survives untouched (continuity residual `1e-4`; 20 000-throat Born ensemble at noise through the kicked retarded evolution). The entangled sector obeys the same structure: exact marginal invariance in the linear part (`3e-15`), equilibrium dBB ensembles SIGNAL-LOCAL (KS 0.015 ≈ noise), NON-equilibrium ensembles SIGNAL (KS 0.045 — the Valentini boundary), and the mean-field gravitational dressing violates marginal invariance at O(G) instantaneously (the Gisin channel, EXHIBITED) but confined BEHIND THE FRONT when retarded (arrival 1.9 vs geometric 2.0)** | Executes the no-signaling edge of the #200 register item "nonlinear measurement theory" (#198 condition 2) — the second named frontier item after #203 (deliverable `docs/nonlinear_no_signaling_audit.md`). **Why by construction:** the linear no-signaling theorem cannot protect BAM (it is nonlinear), and the Gisin proof cannot automatically convict it (it uses the projection postulate; BAM measurement is dynamical transport, #198) — so the audit is run channel by channel on the live ψ–Φ–q dynamics. **The design:** two throat-packets separated by 30; Alice applies a local REAL potential pulse (a local unitary); Bob measures `D_B(t) = ‖ρ_kick − ρ_no-kick‖_{L²(B)}`; every comparison kick-vs-no-kick at the same couplings. **The channel (t = 3):** floor (G=0) `3.3e-13`; Newton `1.5e-6`; half-G `6.9e-7` (ratio 2.14 — linear response); Φ-clamp `2.4e-13` — at the floor: one channel, gravitational strength, carried by Φ alone. **The cone:** the wave-equation solver recovers the Poisson runs exactly at c→∞ (same mean-subtracted spectral operator); pre-front response `~5e-15` (e.g. c=8 at t=2.5 vs Newton's `9e-7`); fronts at `c·t = 25.2/26.4/27.2` vs geometric 25, monotone across a factor-2 speed sweep — the superluminality belongs to the NEWTONIAN APPROXIMATION, not the theory (the 5D theory is causal: #199 verified the Bianchi structure whose weak-field gauge-fixed form IS the wave equation). **The coexistence:** #198's theorem needed only the REALITY of the potentials — retardation preserves it: norm exact, continuity `1e-4`, Born ensemble at noise on the kicked retarded flow — no-signaling and the Born rule hold SIMULTANEOUSLY. **The entangled sector** (the #198 effective two-throat ψ(x₁,x₂), where Gisin's theorem actually lives): (a) linear — a local unitary on throat 1 leaves the x₂-marginal invariant to `3e-15` (the partial-trace theorem on the discrete flow); (b) equilibrium dBB trajectories through a branch-crossing: KS 0.015 vs noise 0.018 — signal-local; (c) branch-only (non-Born) preparation: KS 0.045 = 2.5× noise — Alice's choice IS visible at Bob out of equilibrium: Valentini's signal-locality boundary reproduced on the BAM transport (equilibrium is load-bearing, with the #198 H-theorem relaxation as the mechanism — exactly dBB's position); (d) nonlinear — the BAM mean field `Φ[ρ₁+ρ₂]` violates marginal invariance at O(G) from t = 0.6 when instantaneous (the Gisin channel, exhibited against an EXACTLY-zero linear background) and stays at machine floor until t = 1.9 ≈ d/c = 2.0 when retarded. **Scope:** minimal causal completion (full GR constraint analysis not run); the configuration-space pair is the #198 EFFECTIVE description — the register item is NARROWED (no-signaling edge audited), not closed (emergence of configuration space remains); equilibrium is a hypothesis with a mechanism; non-equilibrium signaling is a prediction shared with the whole dBB program; 1D/2D reductions, c a model parameter (`nonlinear_no_signaling_audit_probe`, PR #204) |
| **The coupled 5D+soliton solve** — the confrontation, the bound, and the NR target | **THE REFUTATION EDGE FIRED AT THE WEAK-FIELD LEVEL — AND THE RESULT IS CLEAN. With no knobs left (the #202 law is EXACT: `m_e/m_μ = (3/7)·(r_s/σ_mode)`, matching constant zero; the #180 solution LOCKED), the coupled weak-field solve gives `σ_mode/r_core = 4.6–6.7` (pairing definitions; full band 1.2–6.7) versus the REQUIRED `88.6` (conv A): **m_e OVER-PREDICTED by ~×13–19**. The direction is right — the true 5D core is the strong-field endpoint of the #179 runaway, SMALLER than the weak-field q-core, so the weak-field value is an UPPER BOUND on m_e (which holds) — and the gap is PROVABLY NOT CLOSABLE inside the weak-field model: the binding sweep gives `RMS/r_q = 3.17 → 1.53 → 1.25` as `Φ(0) = −2.6 → −6.9` — the ratio moves the WRONG way. The mass-ladder thread closes onto ONE falsifiable number: the NR core contraction `r_q/r_s = 13–45×`, with the failure mode stated (an O(1) contraction refutes the pairing mechanism as m_e's quantitative origin)** | Executes the #202 register item (deliverable `docs/coupled_5d_soliton_solve.md`; everything computed live from the locked solution — NO new fits anywhere). **The setup:** #202 made the law exact, so the observed μ/e REQUIRES `σ_mode/r_s = 88.6` (conv A; 206.8 in B); the coupled solve extracts BOTH scales from the one locked `relax(3.5, 0.05)` solution: wave extent RMS = 1.273 and the live #201 kernel inversion `R* = 5.02/5.59`; core radius = the repo's own throat definitions — q half-central `r_q = 0.833`, ordering threshold `r(ψ²=ρ_c) = 1.083`. **The confrontation:** the ratio table spans 1.2–6.7; with the pairing-relevant definitions 4.6–6.7 vs 88.6 — the weak-field solve predicts `m_e/m_μ ≈ 0.065–0.093` vs observed `0.00484`: over-prediction ×13–19. Stated plainly, as the result. **The structure of the failure:** (a) the BOUND direction is right — m_e ∝ r_core and the true core is the strong-field collapse endpoint (#179's runaway branch: q self-gravity overwhelming the quartic saturation), so weak-field bounds m_e from above ✓; (b) NON-CLOSABILITY, measured — sweeping binding strength at locked couplings, the wave compacts FASTER than the threshold core (`ψ² > ρ_c` structural), so the ratio decreases toward strong binding: no corner of the weak-field family reaches 88.6; the missing factor is physics ABSENT from the model. **The NR target:** the full 5D numerical-relativity core solve must yield `r_q(weak)/r_s(true) = 13–45` (convention/definition band) for the pairing mechanism to land m_e/m_μ; if NR gives O(1) contraction, the mechanism is REFUTED as the quantitative origin of the electron mass (the smallness mechanism — index protection, multiplicative structure, naturalness — survives; its numerical anchor does not). **The thread, closed onto one number:** #192 (fine-tuning found) → #193 (not kinematic) → #194 (dialed) → #195 (the index mechanism) → #201 (the rebuild; diagnostics collapse 74.7→4.5) → #202 (the exact law; sensitivity 1) → #203 (the confrontation): every step exact, measured, or bounded; no knobs added anywhere; either outcome of the one remaining computation is progress. **Scope:** locked #180 parameters (`Φ(0) = −4.2`: the weak-field label is itself strained — NR is not optional); `r_q` grid saturation absorbed in the band; the #179 runaway is a LOCATION, not a computation — exactly what the target formalizes (`coupled_5d_soliton_solve_probe`, PR #203) |
| **The 5D throat-core solve** — the exact suppression law on the Tangherlini bridge | **#201'S FITTED EXPONENT IS NOW A LAW. On the t=0 slice of the J-quotiented Tangherlini throat the k-winding problem is 1D on the bridge `ρ(σ)=√(r_s²+σ²)`, and three exact results follow: (1) THE PIN-TWISTED BOUNDARY CONDITION — the deck map = (σ→−σ)×((−1)^k on the fiber): ODD-k MODES HAVE A NODE AT THE CROSS-CAP (quotient-spectrum equality verified to 1e-5; the geometric realization of #195's forbidden one-mouth lift); (2) THE EXACT E≈0 RADIAL LAW — for k=1 the regular mode is EXACTLY `φ(σ)=σ` (closed form, verified 3e-8), far exponents = k to 1e-3, the shallow-well tail theorem at 2%; (3) THE SUPPRESSION LAW `ε_k = (r_s/σ_mode)^k·e^{c₀(k)}`, `c₀ = {0 (exact), −0.405, −0.783}` — the 5D DERIVATION of #201's `e^{−kc}` with `c = ln(σ_mode/r_s)`. Consequences: the electron mass is a LINEAR readout of the throat-to-mode hierarchy (physical sensitivity = −1 EXACTLY: the naturalness chain closes `74.7 → 4.48 → 1`); the inversion is exact (`c₀(1)=0`): `σ_mode/r_s = 88.6` (conv A) — ONE dimensionless number now separates the program from an outright m_e/m_μ prediction; the impossibility bound hardens (`ε₃/ε₁ ≈ 2e-4` vs 88.6 needed)** | Executes the #200/#201 register item (deliverable `docs/throat_core_5d_solve.md`). **The bridge:** the proper-distance identity `d/dρ√(ρ²−r_s²) = 1/√f` verified pointwise to 1e-16; the throat interior is the 1D problem `(ρ³φ′)′ = k(k+2)ρφ` (+E) with the winding barrier `k(k+2)/ρ²`. **The boundary condition, derived:** the quotient projects modes onto `φ(−σ)·(−1)^k = φ(σ)` — odd k ⟹ Dirichlet at the neck, even k ⟹ Neumann; machine-checked by QUOTIENT-SPECTRUM EQUALITY (half-bridge Dirichlet spectrum ≡ full-bridge odd-parity spectrum to 1e-5, and Neumann ≡ even-parity): the fermionic modes vanish at the identification locus — #195's angular-momentum prohibition of the one-mouth mass lift, geometrized. **The exact law:** `φ = σ` solves the k=1 equation identically (`(ρ³)′ = 3ρσ`); general k grows as `σ^k` (the decaying partner `σ^{−(k+2)}`; `p(p+2)=k(k+2)`); a genuinely shallow bound state (`E = −5.6e-4`, `\|E\| ≪ k(k+2)/σ_c²`) has its interior tail EQUAL to the regular solution pointwise (ratio constant to 1.9% over `σ ∈ [2,20]`) — the E≈0 power-law regime is the physical one (the deep-binding `√(k(k+2))` WKB form does not apply). **The suppression law:** slopes `dS/d ln σ = {1.0000, 2.9996, 4.9992}` and matching constants `c₀ = {0.0000, −0.4054, −0.7827}` — computed, not fitted; hence `ε_k = (r_s/σ_mode)^k e^{c₀(k)}`: #201's ansatz derived with `c = ln(σ_mode/r_s)`. **Naturalness in the physical variable:** `d ln ε_k/d ln σ_mode = −k` (measured `{1.0000, 2.9999, 4.9997}` at the physical point) — for the electron `\|Δ\| = 1` exactly: the #201 value 4.48 was a log-parametrization artifact; the chain from #194 closes `74.7 (dialed) → 4.48 (log) → 1 (physical)`. **The inversion, exact:** `c₀(1)=0` leaves no O(1) fudge: `σ_mode/r_s = 1/ε₁ = 88.6` (conv A; 206.8 in B) — the throat's GR core radius is ~0.5–1% of the particle's wave extent, and m_e/m_μ measures that hierarchy LINEARLY; physically the right shape (point-like core inside the Compton-scale cloud); the microphysical ratio awaits the coupled 5D+soliton solve — the remaining unknown is ONE dimensionless number governed by an exact law. **The impossibility bound, exact:** the full computed `ε₃/ε₁ = 1.9e-4` (the law: `(r_s/σ)²e^{−0.405}`) vs the 88.6 a pairing hierarchy would need — hardened. **The zero-winding contrast:** the k=0 (vacuum/boson) Neumann solution is exactly FLAT (`φ = const`): O(1) neck amplitude, ×40 above the k=1 sector — only the uncharged channel touches the cross-cap; the mass protection is specific to the winding sectors. **Scope:** rigid vacuum Tangherlini; E≈0 regime (verified); `σ_mode` from the soliton sector (the vacuum barrier is repulsive, as it must be); scalar reduction of the winding sector (`throat_core_5d_solve_probe`, PR #202) |
| **The Dirac-tower mass ladder** — un-dialing the electron (#200 register item 1) | **THE #194 DIAL IS REMOVED. The electron level is rebuilt as a MULTIPLICATIVE chain on the #195 index-protected zero mode — `m_e = ε₁·o₁·S₁` with the bare level EXACTLY zero (Atiyah–Singer, re-verified `1e-10`), the one-mouth additive lift FORBIDDEN (angular momentum), and inter-sector mixing SUPERSELECTED (the surrogate's `t₁₃ = 14.85` — the term whose near-cancellation was the dial — has NO counterpart). The three fine-tuning diagnostics COLLAPSE: Barbieri–Giudice `74.7 → 4.48` (= the neck aspect; every other entry ≤ 1; ZERO sign flips in 2000 ±25% draws vs flipping at ±2%); Berger sensitivity `−70.9 → +4.15` with NO λ_break on `(0,∞)` (vs 0.986); Monte Carlo `P(m_e ≤ obs) = 7.7% sliver → 49.6% TYPICAL`. The fitted coupling IS O(1) geometry: neck aspect `ℓ/a = 4.5–5.3`; soliton-kernel inversion `K(R*) = ε₁` on the ACTUAL #180 profile at `R* = 3.9–4.4 × RMS`. And the IMPOSSIBILITY BOUND: pairing-only μ/e needs `ε₃/ε₁ = 88.6`, tunneling gives `e^{−2c} ≈ 1e-4` — the hierarchy CANNOT be pairing: smallness = geometry, ratios = dynamics (natural per #194)** | Executes the first item of the #200 open-items register (deliverable `docs/dirac_tower_mass_ladder.md`). **The rebuilt model:** heavy sector = the surrogate's natural 2×2 {k=3,5} block with the uplift honestly REFIT after the k=1 decoupling (removing the electron's level repulsion shifts μ; `β → 0.79·50π` — the uplift remains the fitted dynamical knob, whose sensitivities #194 already showed natural); electron = `ε₁·o₁·S₁` with `o₁ = 1` (the #195 pairing overlap, exact for k=1) and `S₁ = (3/7)·m_μ` (the k=1 bare scale tied to k=3 by the #197 Dirac-tower ratio; convention stated, band carried: conv B `S₁ = m_μ` shifts `ℓ/a` 4.48→5.33 with no conclusion changing); fitting μ/e = 206.77 fixes `ε₁ = (7/3)/206.77 = 0.01129` — INDEPENDENT of the heavy refit. Three fitted numbers for three masses, same count as the surrogate — the content is the structure of the fitted point. **What ε₁ is** (constrained, NOT derived — anti-rigging): `e^{−ℓ/a}` with `ℓ/a ≈ 4.5` (a throat neck a few radii long); and the #185 GR overlap kernel on the actual #180 self-gravitating soliton (RMS 1.274) gives `K(R*) = ε₁` at `R* ≈ 5.0` = 3.9 soliton radii — the coupling the electron mass requires sits exactly where the repo's own machinery puts a few-radii mouth separation; deriving `ℓ/a` from the 5D core solve is the remaining register step. **The diagnostics** (all re-run): BG table — worst `Δ(m_e) = 4.484` (= c, the exponent of a positive exponential), heavy entries ≤ 0.89, and m_e is a product of positive factors: 0 sign flips in 2000 ±25% draws; Berger — `m_e(λ) = e^{−c/λ}(1/λ + λ/2)` (fiber-riding neck × the #197 tower factor): `d ln m_e/dλ\|₁ = c − 1/3 = 4.15`, positive on all of `(0,∞)`: the #192 signature (sensitivity = inverse distance to a zero crossing) is eliminated because THERE IS NO ZERO CROSSING; MC — `ln m_e` smoothly distributed (width `c/4 ≈ 1.1`), `P(≤ obs) = 0.496`, no cliff. **The impossibility bound:** pairing masses scale `ε_k = e^{−kc}` (the winding-k mode decays k× faster through the neck) — DECREASING in k; μ/e from pairing alone needs `ε₃/ε₁ = 88.6` but tunneling amplitudes cannot grow with the barrier (`e^{−2c} = 1.3e-4`): the inter-generation hierarchy provably CANNOT be mouth pairing — it stays with the dynamical uplift (the #122/#136 phase budget), fitted and natural. The #193/#197 division of labor (structure kinematic, hierarchy dynamical) is now realized in the mass model itself: after the rebuild EVERY sensitivity in the ladder is O(few) or below — fully natural, hierarchy parametrized but not tuned. **Scope:** ε₁ constrained not derived (the 5D core solve is the register item); S₁ convention stated with band; heavy sector remains the surrogate's fitted-and-natural dynamics (`dirac_tower_mass_ladder_probe`, PR #201) |
| **PR #200 — the pair-creation cobordism, constructed** — completing "Pauli from GR"; the release capstone | **THE MILESTONE: the ONE open construction of the #196 adjudication is CLOSED — the BAM pair-creation cobordism is EXHIBITED: `M = (S³×I) + H²₍₊₂₎ + H²₍₋₂₎` (two 2-handles on a two-component UNLINK, framings ±2), boundary `S³ → L(2,1) # L(−2,1) = RP³ # RP³` — pair creation from the ACTUAL closed background. Machine-checked: linking matrix `diag(2,−2)`, `\|det\| = 4 = \|H₁(RP³#RP³)\|`, Smith form `diag(2,2)` (coker = Z₂⊕Z₂ EXACTLY), EVEN framings ⟹ SPIN cobordism, signature 0, ONLY index-2 Morse points ⟹ causally continuous (avoids the PROVEN-discontinuous classes {1,3}) — the Sorkin selection rule SELECTS pair creation, the channel #58's threshold always assumed; the ±2 mirror structure IS the C-conjugate throat–antithroat pair. With it, "PAULI FROM GR + the forced Pin⁻ framing" is a CONSTRUCTED THEOREM, and the RELEASE LEDGER is re-verified green in one run (#183/#196 algebra; #193 scalar ladder; #195 index zero mode; #197 Dirac ladder; #198/#199 guidance current)** | The capstone (deliverable `docs/pair_creation_cobordism_capstone.md`; the probe machine-checks the construction and re-verifies one core invariant per milestone live). **The construction:** integer n-surgery on an unknot = `L(n,1)`; a split link = the connected sum; so attaching 2-handles along the ±2-framed unlink to the final slice of `S³×I` gives `∂M = S³ ⊔ RP³#RP³` — no ℝ³ transplant: built over the S³ background directly. `L(−2,1)` = the mirror `RP³` ≅ `RP³` (amphichirality, the #196 Lemma 2b criterion) — the pair are IDENTICAL geons AND the #58 C-conjugate pair (Hopf charges ±1, Σc₁=0): the cobordism's mirror structure is the charge conservation of the pair threshold. **The Dowker–Sorkin ledger, all four rows closed:** mirror-pair identity ✓ (the framings), bordism existence ✓ (explicit), spin extension ✓ (even framings — the 2-handle obstruction vanishes; the Pin⁻-framed SSC theorem gets its structure), the explicit BAM 4-manifold ✓ (THIS PR). **The selection rule:** Morse index 1 and n−1=3 points are PROVEN causally discontinuous (Dowker–Garcia–Surya) and suppressed in the SOH; our cobordism has only index-2 points — causally continuous under the Borde–Sorkin criterion (our channel needs only the proven part) — while single-throat creation would require the discontinuous classes: topology-change kinematics SELECTS pair creation, agreeing with the #58 energetics for a reason. **The completed chain:** 5D Einstein equations → RP³ geon prime (#169/#196) → prime/non-chiral/abelian ⟹ SSC (#196) → the 4-manifold explicit, spin, causally continuous (#200) → SOH abelian weights, exchange = rotation (#196) → the forced Pin⁻ framing lifts the rotation to `−I` (#196/#188, re-verified) → exchange `−1`: FERMI, Pauli (#185/#187) — end to end, no unconstructed step. **The release ledger (green in this run):** `½trT² = −1` + deck dets; `E_k = 2k+(k/λ)²`; the Atiyah–Singer zero mode (monopole ground = q to 1e-7); the Dirac ladder (`(n+1)(n+2)` multiplicities, `λ×(1)=√6` to 1e-9); the guidance-current identity `T_{μχ} = k·Im(ψ*∂_μψ)` to 1e-10. **The open-items register (stated, 6 items):** the Dirac-tower mass ladder (the hierarchy is still a dialed fit, #194, with the mechanism found, #195); the nonlinear measurement theory (#198 cond 2); the full 5D core dynamics (#199); the cosmological-constant failure (#165) and the γ misfit (#160) kept on the books; the Borde–Sorkin remainder (cited conjecture). **The defining claim of the release:** not that the program is finished — that its deepest steps (geon statistics, the generation spectrum, the Born rule, the guidance law) are now THEOREMS WITH HYPOTHESES instead of imports (`pair_creation_cobordism_capstone_probe`, PR #200) |
| **The guidance law from the 5D bulk** — the throat rides its own conserved charge (discharging #198's condition 1) | **THE GUIDANCE LAW IS DERIVED, NOT ASSUMED: (i) the de Broglie current IS bulk stress-energy — for the fiber-winding KK mode `Ψ = ψe^{ikχ}`, `T_{μχ} = k·Im(ψ*∂_μψ)` identically (verified to `1e-12`, k=1,2,3); (ii) its conservation is the χ-component of the contracted BIANCHI identity of the 5D Einstein equations — verified SYMBOLICALLY AND EXACTLY (sympy) on the weak-field KK metric `diag(−(1+2Φ),(1−2Φ)I₃,λ²)` with arbitrary Φ(t,x): all five components of `∇·G ≡ 0`; (iii) the throat is the quantized, topologically conserved unit of winding (#178/#181/#182) and MUST ride its own conserved current — demonstrated POINTWISE: quantized 2D cores transported through the live nonlinear dynamics move with the ambient `J/ρ` at the core (background `0.3142` and mutual `1/d = 0.1` parts both reproduced; winding exactly ±1 throughout). GR taken whole SELECTS the Bohmian flow: the transport differs from geodesic motion exactly by the quantum potential Q — which sits inside the SAME stress tensor (Madelung balance verified to `1e-3`; `\|∇Q\|` dominates `\|∇V\|` ×27 on the live dynamics) — and the law is species-universal (k cancels in `T^i_χ/T^0_χ`). With #198: 5D Einstein → Bianchi → conserved fiber current → topological transport → `v = ∇S` → Born rule** | Target B of the foundational-realism program (deliverable `docs/guidance_law_from_5d.md`; the probe verifies every step). **The chain:** #198 proved `\|ψ\|²` is the unique equivariant measure IF throats move with `v = ∇S`; that identification is here derived from the bulk. (1) The throat mode is the fiber-winding KK mode (the standing #83/#193 reduction). (2) THE STRESS-TENSOR IDENTITY: the mixed fiber components of the 5D stress tensor are identically the de Broglie current with the winding number as charge unit (winding = charge, #42–#44) — the current is not an extra structure, it is a COLUMN of the T that sources the 5D Einstein equations; verified on explicit `(x,χ)` grids, worst deviation `1e-12`, χ-independent (the Killing symmetry). (3) THE BIANCHI STEP: `∂_χ` is Killing, so `𝒥^μ = T^μ_ν ξ^ν` is conserved — and since `∇·G ≡ 0` identically (computed exactly in sympy for the 5D weak-field KK metric: Christoffels, Ricci, G, divergence — all five components identically zero, no expansion), the 5D Einstein equations FORCE `∇_μT^μ_χ = 0`: the continuity equation behind the Born rule (#198 Theorem 1, live residual `3e-4` re-verified) is a BIANCHI IDENTITY of the bulk — GR's 'automatic conservation machine', not a model assumption. (4) THE TOPOLOGICAL TRANSPORT: the throat is a localized INTEGER charge (`Q = (1/2π)∮∇arg ψ`, discrete and unfractionable); it cannot sit still while its winding flows away, so `v^i = 𝒥^i/𝒥^0 = ∇S`; demonstrated pointwise — a vortex–antivortex pair (charges ±1) on a uniform condensate with a quantized background flow, evolved 8 time units: winding EXACTLY ±1 at every sample, measured core velocities `(0.3125, 0.098)/(0.3125, 0.117)` vs ambient-`J/ρ` ring predictions `(0.307, 0.091)/(0.322, 0.092)` — both the background part and the partner-induced `1/d` part reproduced; nothing about the core's motion put in by hand. **The geodesic contrast:** the WKB route gives only geodesics; the derived flow obeys the Madelung–Euler law `∂ₜv + v∂ₓv = −∂ₓ(V_eff + Q)` (verified `1e-3` on the live #198 mid-collision dynamics) with the quantum force DOMINANT (×27) — and Q is the gradient part of the same wave stress-energy: geodesic transport (dropping Q) is exactly the wrong flow #198 showed destroys equivariance; GR taken whole selects the Bohmian flow. **Universality:** `v = T^i_χ/T^0_χ` cancels k (spread `1e-16`, k=1,2,3) — one guidance law for all species, no per-particle constant, exactly as dBB requires. **The upshot:** the chain now runs 5D Einstein equations → Bianchi → conserved fiber current → topological charge transport → `v = ∇S` → (#198) `\|ψ\|²` equivariant and unique → Born rule: the dBB-grade interpretation is a CONSEQUENCE of the bulk field equations. **Scope:** the full 5D throat-core dynamics is not solved (the core is represented by its exact conserved winding + soliton-scale localization; the topological argument constrains the geometric core to its charge, but constrained ≠ solved); #198's equilibrium/measurement conditions unchanged (this discharges condition 1, not 2); 2D transverse + 1D longitudinal demonstrations (a combined 3D+fiber run is computational, not conceptual); requires sympy (`guidance_law_from_5d_probe`, PR #199) |
| **The Born-rule equivariance test** — what measure does the BAM transport preserve? | **THE DEEPEST REFUTATION VECTOR DID NOT FIRE: `\|ψ\|²` is EXACTLY equivariant under the real-time BAM ψ–Φ–q transport flow, and it is the UNIQUE equivariant density among functionals of ρ. Theorem 1: the polar decomposition of the pilot equation gives `∂ₜρ + ∇·(ρ∇S) = 0` exactly — BECAUSE every BAM coupling (the self-consistent Φ, the order field q) enters as a REAL potential; verified on the live dynamics (residual `3e-4`; a 20 000-throat Born ensemble stays at sampling noise, KS ≤ 0.009 vs noise 0.007, through a two-soliton collision). Theorem 2: `∂ₜh + ∇·(hv) = (h−ρh′)∇·v` ⟹ only `h ∝ ρ` (compressible flow; `√ρ`- and `ρ²`-ensembles FAIL at KS 0.10/0.15; wrong flows FAIL at 0.21/0.40). FALSIFIABLE: a dissipative deformation produces EXACTLY the predicted continuity source and breaks the transport. RELAXATION: a uniform ensemble falls 3.48 nats toward Born — fixed point AND attractor. The Born rule at dBB GRADE — the flagged deepest import replaced by a theorem with stated hypotheses** | Item 3 of the theorem-shaped program (deliverable `docs/born_rule_equivariance.md`; the probe verifies every claim on the running dynamics). **The question:** for a classical-field-foundational theory the one honest route to the Born rule is EQUIVARIANCE (Dürr–Goldstein–Zanghì): is `\|ψ\|²` the measure preserved by the BAM transport, or does the ψ–Φ–q dynamics equivariantly transport some OTHER functional? Anything else ⟹ BAM refuted at its foundation. **The flow:** the real-time Hamiltonian flow of the #179/#180 energy functional — `i∂ₜψ = −½∇²ψ + [Φ + ½gq²]ψ`, self-consistent Poisson gravity, live real order field — with the throat transported by the local momentum field `v = ∇S = J/ρ` (the guidance identification: Galilean structure, the exact Ehrenfest relation `d⟨x⟩/dt = ⟨∇S⟩_ρ` verified on the full nonlinear dynamics for a moving soliton, and uniqueness of the current that closes continuity); norm conserved exactly (real potentials, unitary split-step). **Theorem 1 (equivariance):** `∂ₜρ + ∇·(ρ∇S) = 0` EXACTLY — the nonlinearity of Φ[ρ] and the q-coupling are irrelevant because real potentials move the phase, never the modulus; grid residual `3.3e-4` (integrator error) on the live evolution; the Born ensemble's KS series `[0.0065, 0.0068, 0.0091, 0.0087]` vs sampling noise `0.0071` through the collision. **Theorem 2 (uniqueness):** for any `h(ρ)` the transported residual is `(h−ρh′)∇·v` — zero for all configurations iff `h ∝ ρ`, given compressibility (verified: the identity pointwise on the live flow; `max\|∇·v\|` huge near the collision); the `√ρ` and `ρ²` ensembles depart immediately (final KS `0.095/0.148` vs Born's `0.0087`), the half-speed and frozen flows fail (`0.205/0.396`) — equivariance is a property of the PAIR (`\|ψ\|²`, `∇S`). **The teeth:** deforming by `iγW(x)ψ` (γ=0.15, the structure of an imaginary-time/absorbing coupling) produces EXACTLY the predicted source `−2γWρ` (match to 0.1%), 17% norm loss and a drifting ensemble (KS 0.11) — the repo's own imaginary-time RELAXATION flow is of precisely this non-equivariant type: the solver tool and the physical Hamiltonian flow are distinct, and only the latter carries the Born rule; had the BAM functional contained any complex/dissipative/derivative coupling, this PR would have refuted the Born rule. **Relaxation (Valentini):** a UNIFORM ensemble transported through the collision falls `4.13 → 0.65` in the coarse-grained H-function (3.48 nats) — the Born rule is an attractor of non-equilibrium, not just a fixed point. **The label:** the Born rule enters BAM at dBB grade — equivariance + uniqueness + relaxation, the same epistemic status as in Bohmian mechanics — CONDITIONAL on (1) the guidance identification (its 5D derivation is its own program) and (2) the linear measurement regime (test throat in an external pilot wave; the general nonlinear measurement theory is open); 1D demonstration, dimension-blind theorems (verbatim for the radial 3D #180 system) (`born_rule_equivariance_probe`, PR #198) |
| **The analytic Berger–Dirac ladder** — closed-form spectral geometry for the {1,3,5} sectors | **The Dirac spectrum on the whole Berger family, in CLOSED FORM (Peter–Weyl; every fiber-momentum sector an exact 2×2 block; the spectrum classical — Hitchin/Bär — re-derived and checkpointed): family A `a_j = (2j+1)/λ + λ/2`, family B `b± = λ/2 ± 2√((j+½)² + m′²(λ⁻²−1))`. The odd-k ladder is the WINDING TOWER `m_k(λ) = k/λ + λ/2` with UNIFORM gaps `2/λ` — ordered and gapped at EVERY λ, O(1) round-point sensitivities (`−1/3, −5/7, −9/11`: the #192 fine-tuning has NO spectral counterpart), and every crossing LOCATED: character changes at `λ×(k) = √(2k+4)` (√6 first), harmonic-spinor masslessness at `λ*(k) = 5.668, 8.035, 9.851` — the ELECTRON SECTOR COLLAPSES FIRST; nothing below λ=2 (Lichnerowicz). The k₅ = 5 cutoff has NO spectral counterpart at any λ — the generation count stays dynamical** | Item 2 of the theorem-shaped program: upgrade the ladder protection from algebra (#183), surrogate numerics (#192), and the scalar operator (#193) to the ACTUAL spinor field content, analytically, on the whole family. **The derivation** (deliverable `docs/berger_dirac_analytic_ladder.md`): left trivialization on SU(2), Koszul coefficients `Γ₁₂₃=Γ₂₃₁=λ, Γ₃₁₂=2/λ−λ`, `𝒟 = (σ₊J₋+σ₋J₊) + (2/λ)σ₃J₃ + (λ/2+1/λ)`; total fiber momentum `m′` is a good quantum number, each sector an exact 2×2 block. **Machine-validated on four independent checkpoints:** (i) closed forms = the assembled operator matrices to `1e-15` (j ≤ 3, four λ); (ii) at λ=1 the families assemble EXACTLY to the round S³ Dirac spectrum `±(3/2+n)` with multiplicities `(n+1)(n+2)` (counted through n=7); (iii) the λ→0 Gromov–Hausdorff collapse limit reproduces the Dirac spectrum of the base S²(½) (`±2(l+1)`, the Ammann–Bär circle-bundle limit) with the m′≠0 modes diverging as the KK momenta `2\|m′\|/λ`; (iv) LICHNEROWICZ: `𝒟² ≥ scal/4` with `scal = 8−2λ²` forbids zeros for λ<2 — every zero of the closed forms lies at λ>2 (first harmonic spinor at λ=4, the Hitchin phenomenon, located), a theorem-level cross-check; spectral asymmetry off round = the nonzero Berger eta invariant. **The ladder:** the odd-k sector grounds are the winding tower `m_k = k/λ + λ/2` (the KK fiber momentum LINEAR for the first-order operator — the #83 throat winding term again), round values `3/2, 7/2, 11/2`, gaps `2/λ` uniform in k and positive at every λ; ordering `m₁<m₃<m₅` verified on 800 grid points over the whole window `(0, λ*₁)`; squash side absolutely protected (`m_k → k/λ`); NO infinitesimal-squash breakdown — analytically (smooth closed form; sensitivities `(½−k)/(k+½)`, metric-SOFT, closing the #194/#195 loop: the fine-tuning was dynamics, never spectral geometry). **Every crossing located:** the sector minimum changes character at `λ×(k)=√(2k+4)` (exact; below it — the whole squash side + a 145% stretch margin — every sector ground is the pure winding state); masslessness at `λ*(k)² = 8(k+1)+2√(16(k+1)²+k²)` ≈ 5.668, 8.035, 9.851 (verified vs roots to 1e-9) — extreme stretch makes the ELECTRON massless first, μ and τ in order. **The k₅ question answered:** the gap `2/λ` is INDEPENDENT of k at every λ (k=5→7 spacing = k=1→3 spacing exactly) and the collapse boundaries GROW with k — stretch removes sectors from the bottom, never truncates the top: NO spectral counterpart of the 3-generation cutoff anywhere on the Berger family; the cutoff is dynamical (the #122/#136 phase budget), now as a closed-form statement off the round point. **Scope:** round ratios `7/3, 11/7` O(1) (hierarchy stays dynamical, #193–#195); S³ cover with its unique spin structure (on RP³ the odd-k sectors are the Pin-twisted modes; the deck map is a fiber translation — an isometry of every Berger metric — so everything descends); the Berger axis is the fiber/base-separating deformation; the spectrum is classical (Hitchin 1974; Bär GAFA 1996; Ammann–Bär), the checkpointed derivation and the ladder application are the contribution (`berger_dirac_analytic_ladder_probe`, PR #197) |
| **The geon-statistics adjudication** — is the exchange −1 a theorem? (mathematics against the Sorkin-school literature) | **The spin-statistics CORRELATION is a THEOREM for BAM throats — the throat prime is the RP³ GEON (the repo's own #169 quotient) and it passes ALL THREE Dowker–Sorkin hypotheses: PRIME, NON-CHIRAL (`q²≡−1 mod 2`), ABELIAN (`π₁=Z₂`). But the SIGN is framing-dependent, and the honest label changes: RP³ is NON-SPINORIAL in bare Diff (Hendriks: cyclic-π₁ primes — the 2π rotation is isotopic to the identity; a CORRECTION to #170/#171), so bare metric GR + SSC selects BOSE; in PIN⁻-FRAMED GR — the framing BAM's own quotient forces — the rotation lifts to −1 on the pin bundle and SSC selects FERMI. "Pauli from GR" → "Pauli from GR + the forced Pin⁻ framing"** | The theorem-shaped adjudication (deliverable: `docs/geon_statistics_adjudication.md`; the probe machine-checks the arithmetic). **Lemma 1 (the topology):** the #169 involution on the Einstein–Rosen neighborhood S²×ℝ is free and orientation-PRESERVING (`(−1)·(−1)=+1`, verified) — the quotient is the twisted ℝ-bundle over RP² = RP³∖{pt}: the throat prime is the **RP³ geon** (one-sided RP² cross-cap slice), THE canonical example of the geon-statistics literature; n throats ⟹ `Σ = #ⁿRP³`. **Lemma 2 (the hypotheses):** prime ✓ (elliptic, irreducible); non-chiral ✓ (lens criterion `q²≡−1 mod p` at L(2,1); the reflection `diag(−1,1,1,1)` commutes with the deck map and descends — the throat is its own mirror, so the pair-created partner is an IDENTICAL geon); abelian ✓ (`π₁=Z₂`) — the Dowker–Sorkin spin-statistics theorem APPLIES. **Lemma 3 (the correction):** RP³=L(2,1) is NON-SPINORIAL (Hendriks; Giulini: non-spinorial ⟺ lens spaces + handles; cyclic π₁ suffices) — the 2π rotation is isotopic to the identity in Diff, so the #170/#171 sentence "the single geon's 2π rotation acts as −I (spinorial)" is FALSE as a bare-Diff statement; the Friedman–Sorkin spin-½ from pure topology is not available for RP³. **Lemma 4 (the sectors):** the two-throat mapping class group is `G = ⟨s₁,s₂,E \| s₁²=s₂²=E²=1, Es₁E=s₂⟩ = (Z₂∗Z₂)⋊S₂` (Sorkin–Surya; slides s²=e from π₁=Z₂; internal diffeos trivial; the slide subgroup is the FREE product — G infinite): exactly 4 scalar sectors (Bose/Fermi × slide-parity — the Fermi sector EXISTS, #171's consistency survives) plus a CONTINUUM of 2-dim indefinite-statistics sectors (exchange and slides non-commuting, verified in a faithful model, irreducible) — bare frozen-topology GR does NOT select the −1. **Theorem (the selection):** in the SOH with topology change the weights are ABELIAN (kills the indefinite sectors) and pair-created geons obey SSC (exchange = 2π-rotation phase; all hypotheses hold): bare GR — rotation +1 ⟹ **BOSE** (if BAM were bare geometrodynamics the arc would be refuted here); Pin⁻-framed GR — the framing FORCED by #169 (non-orientable RP² slice) + #170 (Pin⁻ unique, `w₂+w₁²=0`) + #195 (the modes ARE pin spinors) — the trivializing isotopy traces the `π₁(SO(3))=Z₂` generator (RP³≅SO(3)) whose SU(2) lift ends at **−I** (computed): spinorial WITH FRAMING ⟹ **FERMI**, the #188 holonomy correctly reinterpreted as the pin lift. **The #58 cobordism ledger:** mirror-pair identity ✓ (amphichirality), bordism existence ✓ (`Ω₃^{SO}=0`), framing extension ✓ (RP³ spin — SO(3) parallelizable; `Ω₃^{Spin}=0`); OPEN: the explicit BAM 4-manifold (the Dowker–Sorkin ℝ³ construction transplanted to the S³ background with their causal/Morse conditions) — a construction, not an obstruction. **The adjudication:** BAM owns the spin-statistics CORRELATION as a theorem (the strongest available form of the claim) and owns the −1 CONDITIONALLY — on the Pin⁻ framing its own construction forces and on the SOH transplant; the bare-GR reading is not survivable (bosons). References: Friedman–Sorkin PRL 44 (1980); Sorkin–Surya IJMPA 13 (1998) 3749; Dowker–Sorkin CQG 15 (1998) 1153; Hendriks (1977); Giulini math-ph/0606066, 0910.2574; Louko–Marolf PRD 58 (1998) (`geon_statistics_adjudication_probe`, PR #196) |
| **The index mechanism** — a Pin/Dirac zero mode for the k=1 sector (answering #194) | **MECHANISM FOUND, with NO new ingredients: BAM throats are Pin⁻ (#183/#188), so the throat mode is a SPINOR — and on the #193 monopole reduction (winding k = charge `q=k/2`) the Dirac operator's Atiyah–Singer index pins EXACTLY `2q = k` chiral zero modes in sector k: sectors {1,3,5} carry {1,3,5} zero modes (verified, residuals ≤ 1e-7), the k=1 electron level is ZERO BY TOPOLOGY. Protection CERTIFIED (two-sided, since `D²⪰0`): gauge wobble pins the zero in `[0, 6e-10]`; metric deformation leaves ker D conformally rigid (`9e-16`, Ω-independent) while the SCALAR ground on the same metric moves `−0.099` — 8 ORDERS of contrast. The natural mass: one-mouth lift FORBIDDEN by angular momentum; the two-mouth pairing gives `\|E_e\| = ε·o` (o = 1.000) — linear, MULTIPLICATIVE, sign-stable, 't Hooft-natural. The #194 'geometric identity' IS the SUSY factorization `D² = A†A` — available only to the spinor** | Answers #194's sharply-posed question: is there structure in the throat geometry that pins a k=1 zero mode? **The mechanism:** the Pin⁻ grading BAM already requires makes the throat mode a spinor; on the base S² it is a charged spinor (monopole `q=k/2`; half-integer q for odd k = the Pin-twisted bundle, #193), and Atiyah–Singer gives index `2q = k`. Computable via the SUSY decomposition `D²₊ = L_{q−½} − (q−½)` (kernel: 2q modes at `l=q−½`), `D²₋ = L_{q+½} + (q+½)` (gapped at `2q+1`) — verified through the validated #193 monopole solver: {1,3,5} zero modes in sectors {1,3,5}, opposite chirality gapped, towers matching the exact Dirac spectrum `λ²=(j+½)²−q²`. **The protection (the discriminator vs #194's cancellation):** flux-preserving gauge wobble (ε=0.05) — the deformed zero mode is CONSTRUCTED in closed form (`v₀e^{−ε∫f}`) and certified variationally: energy pinned in `[0, 6e-10]` while the wavefunction deforms; metric deformation (in 2D all metric deformations are conformal) — ker D conformally rigid, certificate `9e-16`, Ω-INDEPENDENT; the scalar contrast on the same deformed metric: ground `1.500 → 1.401` (O(ε)) — 8 orders between energy-PINNED and energy-TUNED (the #192 squash moved the surrogate's electron because the surrogate is a scalar-type cancellation; the spinor zero would not have moved); flux-change control: the count jumps ONLY with a flux quantum (1→3→5) — the #183 smooth-vs-topology dichotomy realized at the spectrum level. **The natural mass:** within one mouth a first-order lift is FORBIDDEN (zero modes at `j=q−½`; the opposite-chirality tower starts at `j=q+½` — no partner, a Weyl-like protection); the lift requires the throat's TWO MOUTHS (±k winding ⟹ opposite charge ⟹ opposite chirality — the BAM wormhole supplies the Dirac pair): for k=1 both zero modes are the l=0 constants, antipodal pairing overlap `o=1.000`, `\|E_e\| = ε·o` — linear in the mouth coupling, multiplicative and sign-stable ('t Hooft: ε→0 restores the chiral pair, radiatively protected). CONTRAST #194: a sign-flipping DIFFERENCE of two O(7) numbers (Δ=74.7) vs a PRODUCT — small because the mouths couple weakly, not because two big numbers cancel; the surrogate's dialing is an artifact of treating a spinor problem with scalar dynamics. **Noted, not built on:** {1,3,5} zero modes coincide with the generation depths. **Scope:** mechanism established; the mass LADDER is NOT re-derived — identifying the surrogate's k=1 state with the spinor zero mode, computing ε from the throat overlap machinery (#185/#190), and rebuilding the ladder on the Dirac tower is the follow-up; base = round/conformally-deformed S² (the full Berger 3-geometry Dirac spectrum is a further check) (`k1_zero_mode_index_mechanism_probe`, PR #195) |
| **Attacking the fine-tuning** — the electron near-zero: stabilized or dialed? | **DIALED, NOT STABILIZED — and the dial is identified. Every protection candidate is EXCLUDED: no chiral/sublattice conjugation (`tr H = 736 ≠ 0`), no spectral reflection, no index-like structured zero-mode (best overlap `0.917`), and E_e FLIPS SIGN under ±2% transport (`+0.43 → +0.20 → −0.03`) — a cancellation, not a sign-stable seesaw. The tuning is quantified (global `Δ = \|∇ln E_e\| = 74.7`; a ~7% linear-measure accident, Monte Carlo `P = 0.077` vs estimate `0.060`, FLAT through zero — no attractor) and LOCALIZED to ONE codimension-1 combination (76% transport vs 42% base action + 41% slope), whose fiber-map contraction `+71.1` reproduces the #192 Berger λ-sensitivity — the two probes see the same dial. Origin: the CALIBRATION imports the hierarchy (`μ/e = 207` with an O(10) matrix forces `\|E_e\| = E_μ/207`)** | The follow-up to #192/#193: is the near-cancellation protected? **Anatomy:** `E_e = h₁₁ − repulsion = 6.8754 − 6.6758 = 0.1996` — a 2.9% residue between the k=1 diagonal (2π base action + resistance) and the transport repulsion from the k=3,5 rungs (`t₁₃=14.85`, `t₁₅=10.46`); the zero locus `det H=0` is codimension 1, needing no symmetry, only one tuned combination. **Exclusions:** all 8 sign conjugations `SHS=−H` fail (chiral needs zero diagonal); the spectrum is not reflection-symmetric; the near-zero eigenvector `(0.917, 0.398, 0.020)` is pinned to no topological direction; the seesaw-vs-cancellation discriminator (multiplicative suppression is sign-stable) lands on cancellation. **Only the near-zero is tuned:** Barbieri–Giudice `Δ(E_e)` up to 57 (transport), 31/31 (base action/slope), 18 (pinhole), while EVERY `Δ(E_μ)` and `Δ(E_τ)` is < 1 — the heavy rungs are natural; the whole fine-tuning is concentrated in the electron's cancellation. **The Monte Carlo null** (20 000 samples, log-uniform ±25%, fixed seed): `P(\|E_e\| ≤ 0.1996) = 0.077` vs the naive codimension-1 linear-measure estimate `0.060` — same order — with a FLAT E_e histogram through zero (no accumulation = no attractor; no gap = no repulsion): exactly as rare as generic; `P(μ/e ≥ 206.7) = 0.040`. **The cross-check:** `∇ln E_e` contracted with the #192 fiber-map direction = `+71.1`, reproducing the measured Berger λ-sensitivity (+71 on E_e, −70.9 on μ/e) — #192's squash probed an 87%-overlap with the dialed direction. **Origin:** the calibration transfers the tuning from data to parameters — the surrogate CARRIES the hierarchy problem, it does not solve it; and the #193 operator has no near-zero at all (`E₁ ≥ 2`), so the tuning lives entirely in the instanton/transport dynamics, exactly where the dialed direction points. **Numerology guardrail held** (#165 rule): det-zero roots measured vs round constants and NOT matched (transport root `25.538` vs `8π = 25.133`: 1.6%, rejected as derivation). **The sharpened question:** a real solution needs new structure pinning a k=1 zero mode — an index, a chiral grading of winding sectors, or a geometric identity tying the 2π base action to the transport repulsion — finding or refuting it in the throat geometry is the follow-up. **Scope:** the attack is on the locked surrogate (where #192 measured the tuning); MC prior ±25% declared; exclusions cover 3×3-expressible mechanisms — a sub-resolution mechanism in the throat field theory is what the follow-up must probe (`electron_near_zero_naturalness_probe`, PR #194) |
| **The field-theoretic odd-k ladder** — the actual wave operator on the Berger sphere (the follow-up #192 promised) | **The genuine SU(2) Berger Laplacian, sectored by Hopf-fiber winding `k=2m`, has closed-form sector grounds `E_k(λ) = 2k + (k/λ)²` — a Kaluza–Klein split DERIVED from the spectrum (the `(k/λ)²` fiber term IS the #83 mass-operator throat winding term; the `2k` base part is the charge-`q=k/2` MONOPOLE zero-point on the base S², verified by an independent Wu–Yang finite-volume solve to `~2e-7`). The {1,3,5} ladder is ABSOLUTELY protected — positivity `E_k ≥ 2k` and gaps `> 4` in closed form for EVERY `λ ∈ (0,∞)`, NO λ_break — while the mass ratios are pinned O(1) at every λ (μ/e-analog ≤ 3 vs observed 207): STRUCTURE IS KINEMATIC, THE HIERARCHY IS DYNAMICAL — the #192 fine-tuning has no operator counterpart** | Replaces the #192 instanton SURROGATE with the actual deformed wave operator — no ingredient map. **The sector reduction:** a mode's Hopf-fiber winding is `k=2m` (fiber phase `e^{ikθ}`), so the imported #165 spectrum `4j(j+1)+4m²(λ⁻²−1)` is exactly `E(j,k;λ)=4j(j+1)−k²+k²/λ²`, sector ground (`j=k/2`) `E_k(λ)=2k+(k/λ)²`; for fixed k the operator reduces to the Wu–Yang charged Laplacian on the base S² with monopole charge `q=k/2` (the winding IS the charge — the #42–#44 Hopf⟷charge geometry; half-integer q for odd k = the Pin-twisted monopole bundle), solved independently (cell-centered finite volume, natural no-flux pole boundaries): grounds `l(l+1)−q²` at `l=q` give `4q=2k` to `~2e-7`, excited levels `3q+2` match. **Absolute protection:** in closed form for every `λ ∈ (0,∞)`: `E_k ≥ 2k ≥ 2`, `E_{k+2}−E_k = 4+(4k+4)/λ² > 4` — squash stiffens the fiber momenta, stretch floors at the base part, sectors can never cross or vanish (sweep λ∈[0.05,20]: min E₁ `2.0025`, min gaps `4.02/4.04`); the deck grading `(−1)^k` is λ-independent because the ANTIPODE LIES ON THE HOPF FIBER (θ→θ+π), a fiber translation and hence an isometry of every Berger metric — odd k = the antipodally-odd Pin-twisted sector of RP³, the #183 algebra realized spectrally at every λ. Contrast #192: the surrogate's electron level broke at a 1.4% squash; the operator's k=1 sector is bounded below by 2 at ANY squash — the fine-tuning lives in the dynamics, not the kinematics. **The hierarchy is not kinematic:** conformal ratios over ALL λ: `ω₃/ω₁ ∈ [1.53, 3.0]` (round: 2.0) vs observed μ/e `206.8` — factor ≥ 69 away at closest approach; `ω₅/ω₃ ∈ [1.25, 1.67]` vs `16.8`; no Berger deformation of the bare operator reaches the lepton hierarchy. **The two-sided bracket (with #192):** structure from kinematics/topology (absolutely protected), hierarchy from the instanton dynamics (metric-fine-tuned near-cancellation) — the division of labor is now measured. **Controls:** minimal vs conformal coupling (conclusions unchanged); even-k untwisted sectors obey the same closed form; the #192 λ_break reproduced by import. **Scope:** scalar operator, throat = winding sector (not a solved throat geometry); one deformation axis (the fiber/base-separating one); code units (`field_theoretic_odd_k_ladder_probe`, PR #193) |
| **The spectral deformation test** — the {1,3,5} ladder on the Berger-squashed S³ (upgrading #183 from algebra to spectrum) | **The locked lepton Hamiltonian's geometric ingredients, rebuilt on the Berger sphere S³_λ (the #165 SU(2) machinery), keep the {1,3,5} structure over a FINITE window `λ ∈ (0.986, ≥3]` with SMOOTH, LINEAR-at-round-point mass-ratio deformation — the ladder does NOT break at infinitesimal squash, so the #183 protection upgrades from algebra to spectrum. THE DISCOVERY: the electron level is a metric-fine-tuned NEAR-ZERO (`0.1996` vs μ `41.26`, τ `694.98`) that crosses zero at a 1.4% fiber squash, and the μ/e log-sensitivity `−70.9` EQUALS `−1/(1−λ_break)=−71.3` — the steepness IS the proximity to the spectral boundary — while τ/μ is gentle (`+0.82`). Topology guarantees three generations; the round metric tunes the e–μ hierarchy** | Converts the #183 protection claim into the claim the thesis needs. #183 proved the odd-k sector is protected by metric-independent ALGEBRA (deck det, `½trT²=−1`, odd parity) — but algebra is not spectrum: it cannot distinguish "topology protects the spectrum" from "the round metric was doing spectral work the topology can't protect". **The test:** each locked ingredient of the generation block is classified by what it rides on — fiber-riding (× λ: the 2π base action = fiber circumference, the tunnel slope, the winding resistance, the τ uplift in fiber 2π-quanta) vs metric-blind (the Hopf HOLONOMY phase — the very #183 invariant — the base-S² pinhole, the transport prefactor), the map declared BEFORE the sweep; at λ=1 the deformed Hamiltonian reproduces `solved_lepton_masses_mev` to machine precision. **The pass:** three positive ordered levels persist over `λ ∈ (0.986, ≥3]` (dense 81-point stretch-side check); the response at the round point is LINEAR (central-difference slope converges to `−1.47e4`, stable <0.1% from h=1e-4 to 1e-5) — no jump, no divergence, the fail-mode ruled out. **The discovery:** the squash-side boundary is the ELECTRON level's zero crossing at `λ_break=0.98598` (1.4% squash — finite but close); the round-point electron eigenvalue is a near-zero (`0.1996` in action units), so μ/e ≈ 207 is the ratio of a fine-tuned near-cancellation to a normal level; the identity `\|d ln(μ/e)/dλ\| = 1/(1−λ_break)` holds in EVERY ingredient map (`71/71`, `64/64`, `70/70`, `31/31` — default, flip-resistance, flip-uplift, minimal-2π-only), while `d ln(τ/μ)/dλ` stays O(1) or below in all four. **Upshot:** the topology guarantees THREE GENERATIONS (structure metric-robust); the round metric does real spectral work for the e–μ HIERARCHY (metric-fine-tuned, 1.4% from a spectral catastrophe) — the two claims are now separated by measurement. **Scope:** the Hamiltonian is the locked instanton-transition surrogate (not a first-principles wave operator on S³_λ); the fiber/base map is a declared modeling choice (conclusion invariant under flipping the ambiguous assignments); the field-theoretic ladder on S³_λ is the follow-up; action units, electron-calibrated ratios only (`odd_k_ladder_spectral_deformation_probe`, PR #192) |
| **Dynamic two-throat exchange path with back-reaction** | **The exchange `−1` of #188 is the ADIABATIC limit of a real dynamical swap: the throat's Pin spinor driven around the swap loop at FINITE speed `1/T`, under a field that BACK-REACTS on the moving throat. As the swap slows the dynamical geometric phase → `−π` (the exchange `−1`) and the non-adiabatic excitation → 0; at finite speed both deviate by an explicit `O(1/T)` cost that VANISHES adiabatically. The back-reacting field is sourced by the moving throat and acts back, but does NOT alter the adiabatic `−1`** | Asks the dynamical question #188 (adiabatic holonomy) and #189 (static relaxation) leave open: is the exchange `−1` robust at finite swap speed, with a back-reacting field? **The model (effective):** the internal Pin spinor `ψ` driven around the swap loop `n̂(s)=(cos2πs, sin2πs, 0)`, `s=t/T`, over finite duration `T` (speed `1/T`) under `H(s,x)=(Δ/2)n̂·σ + g·x·σ_z`, coupled to a back-reacting field `x(t)` with `ẍ+γẋ+ω²x=−κ⟨σ_z⟩` — the moving throat SOURCES the field (`⟨σ_z⟩` drives `x`) and the field ACTS BACK (`g·x·σ_z` in `H`), a two-way coupling; spinor by RK4, field by Verlet. The geometric phase is extracted by continuously factoring the dynamical phase `φ_geo=arg(⟨ψ(0)|ψ(T)⟩·e^{+iφ_dyn})` (robust, `e^{iφ_dyn}` periodic) and `P_exc=1−|⟨ψ_lower(H_f)|ψ(T)⟩|²` the non-adiabatic leakage. **The invariant:** the loop's exact adiabatic Berry phase is `−½·2π=−π` — the exchange `−1` (the #188 holonomy `e^{iπ}=−1`), independent of speed and back-reaction. **Adiabatic recovery:** as `T→∞`, `φ_geo→−π` (deviation `≈0.02` at `T=1000`) and `P_exc→0` — the full dynamics reproduce the `−1`. **The non-adiabatic cost:** at finite speed `φ_geo` deviates from `−π` by `O(1/T)` (`deviation×T≈29.7` constant across `T=10…300`) and `P_exc` grows with speed — the quantified price of a fast swap. **Back-reaction:** the moving throat sources the field (peak energy `≈0.028` at `T=20`, → `1.5e-4` slow), energy is exchanged and it acts back, but the adiabatic limit (the `−1`) is UNCHANGED — the dynamics add an adiabatically-vanishing cost, not a change of sign. **Scope:** an effective model (the internal spinor + one field mode); the field-resolved real-time two-throat solve (orbitals translating, the #190 Coulomb field re-solved each step, the swap geometry in 3D) is the follow-up; weak-field, code units (`dynamic_two_throat_exchange_probe`, PR #191) |
| **The BAM Coulomb-photon kernel for the two-throat HF** — replacing the Yukawa stand-in | **Replaces the screened-photon (Yukawa) stand-in of #187/#189 with the genuine BAM Coulomb-photon kernel — the UNSCREENED `1/(4πd) ⟷ 1/q²` (the photon propagator, #42–#44, the flat-space limit of the S³ Green function `G(ψ)=((π−ψ)cotψ−½)/(4π²R)`), regulated PROPERLY for an isolated system (the Hockney zero-padded Coulomb, validated to `~0.08%`). The direct energy is now correctly LONG-RANGED `J(R)→1/(4πR)` (the Coulomb tail) while the exchange stays short-ranged; and the #187 overlap-normalized physics is ROBUST — fermion-lower for the repulsive photon, the zero-vector Pauli state at coincidence** | Upgrades the two-throat interaction from a stand-in to the real photon. **The kernel:** the unscreened Coulomb `V(d)=1/(4πd)` ⟷ `1/q²` — the photon propagator BAM derives from the throat-fibre exchange geometry (#42–#44) — verified: the isolated-system Coulomb of a unit point source reproduces `1/(4πd)` to machine precision (`Φ={1.25:0.0637, 2.5:0.0318, …}=1/(4πd)`), and it is the flat-space limit of the BAM S³ scalar Green function `G(ψ)=((π−ψ)cotψ−½)/(4π²R)` (the repo's `s3_green_potential`): near the source `G·4πs=0.957→1`, the Coulomb coefficient (`s=Rψ` the geodesic distance), so on the local weak-field patch the throats see the unscreened photon with S³ curvature corrections `O(1/R²)` carried by `G`. **The regulator:** the unscreened Coulomb is solved by the Hockney zero-padded convolution (density padded to a 2× box, convolved with the free-space `1/(4πr)` Green function) — a proper OPEN-BOUNDARY regulator, NOT physical screening; validated against the analytic Gaussian Coulomb self-energy (`U=0.0167` vs exact `0.0167`, ratio `0.9992`, ~0.08%). **The recomputed energies:** on the #180 orbitals, `J(R)={0:0.063,1:0.054,2:0.038,3:0.026,4:0.020,6:0.013}` is now correctly LONG-RANGED — `J(6)=0.0133≈1/(4π·6)=0.01326` (ratio `1.003`), the point-charge Coulomb tail — unlike the Yukawa stand-in's exponential decay; the exchange `K_ex(R)={0:0.063,1:0.038,2:0.0096,3:0.0012,4:0.0001,6:0}` stays SHORT-ranged (overlap-set), so far-apart throats feel the Coulomb direct field but not the exchange (the physically correct long-range structure the screened stand-in lacked). **The #187 physics survives:** with the overlap-normalized `E±(R)=(J±K_ex)/(1±S²)`, for the repulsive photon the antisymmetric (Pin⁻) branch `E₋` lies BELOW the symmetric `E₊` at every finite separation (`R=1: 0.041<0.056`; `R=2: 0.034<0.040`; `R=3: 0.026<0.027`) — the fermion-lower result of #187, established there with the Yukawa, holds with the real unscreened photon; and at coincidence (`R=0`) `J=K_ex=0.0627` (numerator `0`) and `S→1`, so the antisymmetric state is the ZERO VECTOR (Pauli-forbidden). The statistics are a property of the geometry (the Pin⁻ sign + the overlap structure), not of the interaction's screening. **Scope:** the kernel is the genuine BAM photon (the flat Coulomb limit; the `O(1/R²)` S³ curvature corrections are carried by `G` but not applied — the weak-field local patch); the Hockney is a numerical open-boundary regulator (validated), not screening; the orbitals are the rigid #180 throat-solitons (the self-consistent #189 SCF with the Coulomb kernel is the follow-up); the Yukawa was a faithful short-range stand-in, now upgraded to the real long-ranged photon; weak-field, code units (`bam_coulomb_two_throat_hf_probe`, PR #190) |
| **Self-consistent two-throat Hartree–Fock relaxation** — relax orbitals in each other's direct + exchange field | **Does the #187 follow-up: a genuine HF SCF lets the two same-spin throats DEFORM in each other's direct (Hartree) + exchange (Fock) field, with an orbital-specific SELF-INTERACTION-FREE Fock operator `F_i=h+J_{≠i}−K_{≠i}` consistent with the reported energy. The imaginary-time relaxation lowers the energy MONOTONICALLY to a self-consistent variational fixed point (machine-converged, robust across seeded restarts); the relaxed energy lies `2.54%` below the rigid 1D #187-style reference (`−3.808→−3.905`); the two-throat density polarizes (RMS `2.65→3.06`, fidelity `0.978`); and turning off the non-local Fock exchange (consistent control) raises the energy by `0.567` — the same-spin exchange hole doing real work in the mean field** | Relaxes #187's rigid two-throat orbitals self-consistently. **The mean field:** each orbital is relaxed in an orbital-specific, self-interaction-free Fock operator `F_i = h + J_{≠i} − K_{≠i}` — kinetic + confinement `h`, the DIRECT (Hartree) field `J_{≠i}(x)=∫|φ_{≠i}(x')|²V(x−x')dx'` of the OTHER throat (self-interaction excluded), and the non-local Fock EXCHANGE with the other orbital — which is the EXACT variational derivative of the reported energy `E=Σ⟨i|h|i⟩+(J₀₁−K₀₁)` (the same functional for the operator and the energy, in both the full-HF and Hartree-only branches); two same-spin fermions (the two throats, whose spatial state is the Pin⁻ antisymmetric sector of #185/#188) occupy the two orthonormal orbitals, relaxed by imaginary-time gradient descent. **Convergence + robustness:** the SCF lowers the energy monotonically (`−3.808→−3.905`, final `ΔE~10⁻⁸→0`) to a self-consistent fixed point — the orbitals come to rest in the field they themselves produce — and the fixed point is ROBUST across seeded restarts (five random localized initial orbitals converge to `E≈−3.905`, spread `~8×10⁻³`); so it is a self-consistent VARIATIONAL FIXED POINT (robustly reached, not certified the global ground state). The monotone descent confirms a genuine variational relaxation, not the eigenstate-swapping oscillation a naive diagonalization-SCF shows. **Relaxation (energy lowering):** the rigid 1D #187-STYLE reference (the unrelaxed orbitals evaluated with the full HF energy — the 1D analogue of #187's rigid-orbital evaluation, not the 3D number) gives `E_rigid=−3.808`; relaxing self-consistently gives `E=−3.905`, LOWER by `2.54%` — the variational gain from optimizing the orbital shapes (in #187 the orbitals were held rigid; here they deform). **The orbitals deform:** the two-throat density polarizes/spreads in the mean field — RMS width `2.65→3.06`, and its overlap (fidelity) with the rigid density drops to `0.978<1`; the throats respond to each other's field. **The exchange field matters:** with the CONSISTENT self-interaction-free control (the exchange dropped from BOTH the Fock operator `F_i=h+J_{≠i}` and the energy `E=Σ⟨i|h|i⟩+J₀₁`, so its operator and energy are the same functional), running the SCF Hartree-ONLY gives `E=−3.338`, HIGHER by `0.567` than the full `E_HF=−3.905` — for the two same-spin throats the Fock exchange substantially LOWERS the energy (the exchange hole of #186/#187 keeps the like throats apart, reducing the repulsive direct energy; the `−1` of #185/#188 doing real work in the self-consistent mean field, not only in the rigid kernel). **Scope:** a 1D sandbox SCF — external double-well confinement (a stand-in for the throats' self-binding, the #180 self-gravity), a screened-photon (Yukawa) interaction stand-in, spatial-orbital same-spin (unrestricted) HF, in 1D for tractability; the SCF itself is genuine (imaginary-time, monotone, machine-converged, self-interaction-free operator consistent with the energy) and the qualitative physics robust; the relaxed state is a variational fixed point (robust across seeded restarts, not certified the global ground state) and the 'rigid' reference is the 1D #187-style evaluation (not the 3D number); the full 3D self-gravitating two-throat SCF (relaxing actual #180 ψ–Φ–q throat-solitons) is the follow-up; weak-field, code units (`self_consistent_two_throat_hf_probe`, PR #189) |
| **Adiabatic two-throat exchange holonomy** — measure the Pin⁻ sign along a swap path | **Operationalizes the #185 exchange sign: adiabatically transports the throat spin-½ state along an explicit two-throat SWAP PATH and MEASURES the holonomy = `−I` (the spinor returns to MINUS itself; `⟨ψ\|Hol\|ψ⟩=−1`, Berry phase `π`) — to machine precision (`‖Hol+I‖~10⁻⁶`); robust to a wandering rotation axis (the ℤ₂ homotopy class), with a double-swap (4π) → `+1` and a contractible loop → `+1`; the `−1` is the Pin⁻ monodromy `T²=−I` realized along the path** | Makes the #185 Pin⁻ exchange sign operational — measured, not asserted. **The swap path & spin-statistics:** the relative-coordinate config space of two identical throats is `(ℝ³∖0)/ℤ₂ ≃ RP²×ℝ₊`, whose angular factor `S²/antipodal=RP²` is the BAM antipodal closure itself (#169/#170); the exchange `r→−r` is the generator of `π₁(RP²)=ℤ₂`, and by the Finkelstein–Rubinstein / Friedman–Sorkin spin-statistics theorem for solitons/geons it is HOMOTOPIC to a 2π rotation of one throat (the belt-trick / tether twist). **The holonomy measured:** path-ordering the spin connection along the swap (2π) loop — `dU/ds=−i(ω·σ/2)U` with `ω` the loop's angular velocity — gives the adiabatic holonomy `Hol=−I` to machine precision (`‖Hol+I‖~10⁻⁶`); the throat's spin-½ state returns to minus itself, the measured exchange sign `⟨ψ\|Hol\|ψ⟩=−1`, and the Berry phase (arg of the holonomy eigenvalue) is `π`. **Topological (path-independent):** a swap path with a WANDERING rotation axis (a non-commuting, genuinely path-ordered transport) gives the SAME `−I`, and it converges to `−I` as the transport is refined (`N=250→4000` steps) — any way of carrying out the exchange gives the same `−1`; it is the ℤ₂ homotopy class that matters. **Controls:** a DOUBLE-swap (two exchanges, a 4π rotation) gives `+I` (two fermion exchanges = a boson, `−1×−1=+1`) and a CONTRACTIBLE loop (no net rotation, the throats wiggle but don't exchange) gives `+I` — so the `−1` is specific to the single-swap (odd) class. **The Pin⁻ identification:** the `−1` is the throat monodromy `T=iσ_y`, `T²=−I` (`½ tr T²=−1`; #170/#174/#183) — the Pin⁻ structure on the non-orientable RP² closure (the unique spin structure RP² admits), which makes the throat a SPINOR; the 2π/swap holonomy of a spin-`j` object is `(−1)^{2j}`, so the spin-½ throat gives `−1` (a scalar/bosonic throat, `j=0`, would give `+1` along the same path). The adiabatic holonomy IS this `T²=−I`, now transported along the swap rather than read off the algebra. **Scope:** operationalizes the FR/geon-statistics result — the holonomy is computed exactly (path-ordered SU(2)) and is TOPOLOGICAL (the ℤ₂ class), so the measured `−1` is exact; the swap path is the reduced relative-coordinate/frame model, the spin-statistics connection (exchange ≃ 2π rotation) is the FR theorem cited not re-derived, the throat's spinor (Pin⁻) nature is the #170 result; the adiabatic limit is assumed; complements the #185–#187 spatial exchange kernel / Hartree–Fock energies (`adiabatic_exchange_holonomy_probe`, PR #188) |
| **Two-throat Hartree–Fock sandbox** — direct plus exchange terms | **Convolves the #186 hardened overlap kernels with an interaction to build the OVERLAP-NORMALIZED two-throat HF energy `E±(R)=(J(R)±K_ex(R))/(1±S²)` — the orbital overlap `S=⟨φ_a\|φ_b⟩`, the DIRECT (Hartree) numerator `J=∫∫ρ_aρ_b V` and the EXCHANGE numerator `K_ex=∫∫ττV` (`τ=φ_aφ_b`) all from the actual #180 throat-soliton orbitals; for the REPULSIVE screened `V` the GR-selected Pin⁻ antisymmetric branch `E₋` sits below the boson `E₊` at every finite separation, and at coincidence the antisymmetric state is the ZERO VECTOR (Pauli-FORBIDDEN), not a zero-energy state** | Completes the multi-throat mechanics: dress the #185/#186 exchange kernel with an interaction `V` (a screened-photon Yukawa stand-in for the BAM throat-fibre exchange) and assemble the two-throat Hartree–Fock energy with both terms. **The structure (overlap-normalized):** two displaced throats are NON-orthogonal (overlap `S(R)=⟨φ_a\|φ_b⟩≠0`), so the physical two-body energy is `E±(R)=(J(R)±K_ex(R))/(1±S²)` — the direct (Hartree) numerator `J=∫∫ρ_a(r₁)ρ_b(r₂)V(r₁−r₂)` (sign-independent density–density, the #186 direct channel) and the exchange numerator `K_ex=∫∫τ(r₁)τ(r₂)V(r₁−r₂)` (the overlap-density self-energy, `τ=φ_aφ_b`, the #186 exchange channel); `J, K_ex` are the unnormalized HF numerators and `(1±S²)` the overlap normalization. **The integrals** (3D-FFT Coulomb on the #180 orbitals): `S(R)={0:1.0,1:0.79,2:0.41,3:0.15,4:0.05}`, `J(R)={0:0.039,1:0.031,2:0.017,…}`, `K_ex(R)={0:0.039,1:0.024,2:0.006,…}` — all positive (repulsive `V`), decaying, the direct dominating (`K_ex/J` from 1.0 at contact to ~0.01 at R=4). **The normalized energies & fermion ordering (repulsive V):** `E±=(J±K_ex)/(1±S²)`; for the repulsive screened interaction the antisymmetric (fermion) branch `E₋` sits BELOW the symmetric (boson) `E₊` at every finite separation (`R=1: E₋=0.019<E₊=0.034`; `R=2: 0.014<0.020`; `R=3: 0.008<0.009`) — the exchange hole lowers the GR-selected antisymmetric (Pin⁻) state, the gap closing as the overlap `S→0` (distinguishable, nearly degenerate); the ordering is SCOPED to a repulsive `V` (an attractive interaction reverses it; and the normalization itself works against it — dividing `E₋` by `1−S²<1` raises it — yet the fermion branch stays lower across the tested range). **Pauli at coincidence (the zero vector):** as `R→0`, `S→1` and BOTH the numerator `J−K_ex→0` AND the normalization `1−S²→0` — the antisymmetric combination `Ψ₋=(φ_aφ_b−φ_bφ_a)/√(2(1−S²))` is the ZERO VECTOR, i.e. FORBIDDEN (two identical fermions cannot occupy the same orbital), **not** a state with zero interaction energy; the boson `E₊=(J+K_ex)/(1+S²)` survives (bunching). For a CONTACT `V` the numerator `J−K_ex=0` at all R (`J=K_ex=g·D(R)`, the hardened #186 direct overlap), so `E₋=0` at every finite separation (the exchange exactly cancels the direct — the Pauli hole removes the contact interaction), the state forbidden only at exact coincidence. **Controls + convergence:** at far separation `S,J,K_ex→0` (orthogonal, distinguishable throats, `E₊≈E₋`); under 3D-grid refinement `N=64→80→96` the energies converge to `<0.1%` (the precise values carry the #186 soliton-profile ~3% uncertainty). **Scope:** a SANDBOX — rigid #180 orbitals (the self-consistent two-throat HF solve is a follow-up), a screened-photon (Yukawa) regulated stand-in for the BAM Coulomb/photon `1/q²` interaction (#42–#44), spatial (orbital) exchange only (the spin/statistics factor is the separate Pin⁻ `−1`), energies in code units; the overlap-normalized structure (S(R), the direct/exchange numerators, the `(1±S²)` normalization, fermion-lower for a repulsive V, and the forbidden zero-vector antisymmetric state at coincidence) is robust; weak-field (`two_throat_hartree_fock_probe`, PR #187) |
| **Rigid soliton exchange-kernel hardening** — convergence, normalization, direct-term controls | **Hardens #185's rigid-soliton exchange kernel into a trustworthy benchmark: NORMALIZED (orbital `∫\|φ\|²d³r=1.000000`, self-overlap `K(0)=1.001` reproducing the norm to 0.1%, parity `K(R)=K(−R)` exact, Cauchy–Schwarz bound `K(R)≤1`); CONVERGENT (the overlap integral converges to `<0.01%` under quadrature refinement, with the `~3%` soliton-profile core grid-sensitivity — the #180 caveat — honestly identified as the dominant uncertainty); DIRECT-TERM CONTROLLED (the sign-independent direct density-overlap `D(R)` separated from the ±-carrying exchange `K(R)`, both vanishing at far separation, the Pin⁻ `−1` living purely in the exchange channel)** | The hardening probe for #185 (as #177 was for #176), with the three things a credible overlap kernel needs. **Normalization:** the single-throat orbital is normalized `∫\|φ\|²d³r=1.000000`; the self-overlap `K(0)=1.001` reproduces the norm to 0.1% (the residual is the overlap-quadrature discretization); the kernel is parity-symmetric `K(2)=K(−2)=0.40963` (φ radial); and obeys the Cauchy–Schwarz bound `0≤K(R)≤K(0)` (`K(2)=0.41<1`). **Overlap-grid convergence:** refining the angular/radial quadrature `(Nr,Nu)`, `K(2)=0.40977→0.40963→0.40960→0.40959` — the production grid `(320×160)` differs from the finest `(1000×500)` by `0.01%`, converged to `<0.1%`. **Soliton-grid convergence:** rebuilding the #180 soliton at `N=240→320`, `K(2)=0.4096→0.4215`, a `~2.9%` shift — NOT the overlap integral but the soliton profile's core grid-sensitivity (the documented #180 caveat; the sharp core carries ~10% per-point grid uncertainty, integral overlaps a few %); the kernel's shape/scale are trustworthy to a few %, the soliton profile the limiting factor (an inherited, known limitation, not a flaw in the kernel). **Direct vs exchange:** the DIRECT density-overlap `D(R)=∫ρ_aρ_b d³x` (the sign-independent Hartree channel) and the EXCHANGE amplitude-overlap `K(R)` (the ±-carrying channel) are distinct GR-geometric kernels — `D={0:0.065, 1:0.038, 2:0.008, 3:0.001,…}` vs the exchange weight `X=K²={0:1.0, 1:0.63, 2:0.17, 3:0.023,…}` — both decaying monotonically to zero at large R, the direct FASTER (product of densities vs amplitude overlap). **Direct-term controls:** at far separation `R=6` both `D` and `X→0` (distinguishable throats, no exchange); the direct term is a positive density overlap with NO sign — identical for the boson (+) and fermion (−) sectors — while the exchange carries the Pin⁻ `−1`, so the −1 lives PURELY in the exchange channel (the direct is the boson=fermion control isolating it); consequently the two-body energy splits as `E=E_direct ∓ E_exchange`, the structure the #187 Hartree–Fock sandbox evaluates against an interaction `V`. **Scope:** hardens the RIGID soliton-overlap kernel (two rigid #180 orbitals at separation R); convolving the overlaps with an interaction `V` to get the Hartree (direct) and exchange ENERGIES is PR #187; relaxing the orbitals self-consistently is a further follow-up; weak-field (`rigid_soliton_exchange_kernel_hardening_probe`, PR #186) |
| **Multi-throat mechanics & the exchange kernel from the GR ψ–Φ–q soliton** | **The two-throat exchange kernel derived from GR: `K_exchange(R) = (−1)·K(R)` — a GR-geometric SPATIAL overlap `K(R)` computed from the ACTUAL #180 self-gravitating throat-soliton (decays over the soliton size RMS≈1.27, setting the exchange range) times a TOPOLOGICAL SIGN `−1` from the non-orientable Pin⁻ geon statistics (the swap large-diffeomorphism ≃ a 2π rotation = `T²=−I`); the geometry SELECTS the antisymmetric Fermi sector — the two-throat state vanishes at coincidence (Pauli, exact) and N throats give `P∝n^{5/3}` (Γ=5/3, the measured #172 EoS)** | Takes two #180 self-gravitating ψ–Φ–q throat-solitons — multi-throat mechanics — and derives the exchange kernel from GR (no postulated statistics), factorized into a geometric overlap and a topological sign. **The exchange operator:** P swaps two identical throats; `P²=1` ⟹ eigenvalues `±1` (boson/fermion); which sector is NOT a free choice — the GR geometry of the swap fixes it. **The spatial kernel:** `K(R)`, the overlap of two actual #180 throat-solitons separated by R, decays smoothly and monotonically `K̂ = {0:1.00, 1:0.79, 2:0.41, 3:0.15, 4:0.045, 6:0.003}` over the throat-soliton size (RMS≈1.27) — a GR-geometric exchange RANGE, not a postulated form factor (throats exchange strongly only when they overlap). **The sign:** `−1` (fermionic), derived from GR — the large diffeomorphism that swaps two throats is homotopic to a 2π ROTATION of one throat (the Friedman–Sorkin / Dowker–Sorkin spin-statistics theorem for geons), and a 2π rotation on the non-orientable Pin⁻ throat is the monodromy `T²=−I` (`½ tr T²=−1`), so the exchange phase is `−1`; a boson would need the orientable `T²=+I` closure the throat does not have (#170/#174/#183). **Pauli exclusion:** with the GR-selected `−1`, the antisymmetric two-throat state `Ψ₋(z₁,z₂)=φ_a(z₁)φ_b(z₂)−φ_a(z₂)φ_b(z₁)` vanishes IDENTICALLY at coincidence `z₁=z₂` (max `\|Ψ₋(z,z)\|=0` to machine precision — the determinant of two equal rows) — two identical throats cannot occupy the same state — while the boson `Ψ₊` does not (it bunches). **The exchange hole + Fermi pressure:** the exchange term `∝K(R)²` carves an exchange hole (coincidence suppressed by `{0.5:0.89, 1:0.63, 2:0.17, 4:0.002}`) of GR range = the soliton size; macroscopically the exclusion fills a degenerate Fermi tower — with the 3D DOS `g(E)∝√E`, `N∝E_F^{3/2}` and `E∝E_F^{5/2}` ⟹ `E∝N^{5/3}` ⟹ `P=(2/3)(E/V)∝n^{5/3}`, polytropic `Γ=1.6667=5/3` — exactly the Fermi EoS measured in #172. The GR-derived exchange kernel is the microscopic origin of the Fermi pressure of throat matter. **Scope:** the exchange SIGN is exact/topological (the Pin⁻ geon statistics, a GR large-diffeomorphism / mapping-class-group representation); the SPATIAL kernel is the rigid #180 soliton-overlap model (the orbitals) — the full two-body GR problem (the actual two-throat metric, the gravitational direct/Hartree term alongside the exchange, and the dynamical swap with back-reaction) is a follow-up; the Fermi index 5/3 is the standard degenerate-gas result, here attributed to the GR-derived exchange kernel (and matching the independently measured #172 EoS); weak-field/semi-dynamical soliton (`multi_throat_exchange_kernel_probe`, PR #185) |
| **α as a protected boundary invariant** — not a continuous tuning parameter | **Reframes the #105/#143 "137 problem": α's DERIVED structure — the charge quantum (the boundary Chern number `c₁=±1`, `\|c₁\|=1`, the Gauss-law charge / integer Hopf number) and the `1/2π` closure loop measure — is a PROTECTED BOUNDARY INVARIANT, topologically robust to smooth deformation (30 smooth boundary diffeomorphisms keep `c₁` the same integer to 5×10⁻⁷, while a generic continuous coupling functional drifts ~7% under the same deformations) and changing ONLY at a topology change (the Berry gap closing); the VALUE `α⁻¹≈137` remains the single EM residual** | Applies the #181/#182/#183 protected-invariant test to the EM coupling. **The decomposition (#105/#143):** the geometry DERIVES α's structure — the charge quantum `\|c₁\|=1` (the integer Hopf number; charge quantization is topological), the `1/2π` closure loop measure (the `2π` in the Schwinger anomaly `a=α/2π`), and the running — but NOT the value `α⁻¹≈137` (the residual; the closure-number scans failed). The sharp question: is that derived structure a protected boundary invariant or a tunable continuum? **The charge quantum is a boundary invariant:** the first Chern number of the BAM Hopf/spin-½ monopole over the boundary S² (the Gauss-law charge `(1/2π)∮F`), computed by the exactly-quantized Fukui–Hatsugai–Suzuki method, `c₁=+1` — an exact integer. **Protected:** across 30 smooth diffeomorphisms of the boundary geometry (`θ→θ+a sinθ+b sin2θ`, monotone) `c₁` stays the same integer (all round to `+1`) to `5×10⁻⁷` — it does not drift. **Not a tuning parameter:** under the SAME deformations a generic continuous coupling functional (the mean monopole potential `⟨A_φ⟩`) drifts `6.8%` on average (up to `15.8%`) while `c₁` moves `5×10⁻⁷` — the discriminator (quantized + deformation-invariant = protected; continuous + drifts = tuning) puts α's charge quantum cleanly on the protected side. **The loop measure + topology change:** the total boundary flux `∮F=2π·c₁` is quantized in units of the closure quantum `2π` (fixing the `2π` of `a=α/2π`); and the charge integer changes ONLY when the Berry gap closes — sweeping the gap parameter `m`, `C(m)=1` while the gap is open (`m<1`) and jumps to `0` exactly at `m=1` where `min\|d\|→0` (the degeneracy crosses the boundary), the EM-boundary analog of `\|q\|=0` (#182) and `½ tr T²=0` (#183). **Unity:** the EM charge quantum is to the boundary what the order-field winding is to the soliton (#181/#182) and the generation sector is to the bulk (#183) — a protected topological charge robust to smooth deformation, changing only at a topology-change event. **Honest scope:** this does NOT derive the value `α⁻¹≈137` — that residual input (the 137 problem) stands; what is established is that the structure around it is specifically PROTECTED, so α should be tested as protected-boundary-structure × one residual scale, not fit as a continuous tuning family. Refines the #105/#143 ledger; completes the protected-invariant picture at the EM sector (`alpha_protected_boundary_invariant_probe`, PR #184) |
| **Odd-k / generation-sector survival under a deformed bulk geometry** — topologically protected | **The odd-k {1,3,5} charged-lepton ladder (#174) is set by metric-INDEPENDENT invariants — the orientability of the antipodal quotient (deck det `−1` for the brane RP², `+1` for the bulk RP³) and the spin-closure class (`½ tr T²=−1`, Pin⁻) — so it survives ANY smooth bulk deformation to machine precision (1000 random orientation-preserving deformations preserve `½ tr T²` and the deck dets to ~10⁻¹⁵), the 3-generation count rides through (`D_bulk=5` + odd parity, both topological), and it changes ONLY at a genuine topology change (`½ tr T²` crosses 0 at θ=π/4) — the bulk-level analog of #181/#182** | The roadmap capstone: does the odd-k generation sector survive when the bulk geometry is deformed? **The grading is topological:** the antipodal deck map is `−I` in any linear frame, so `det=(−1)^dim` — the brane angular slice `S²/antipodal=RP²` is non-orientable (`det=−1`), the bulk `S³/antipodal=RP³` orientable (`det=+1`), for any metric; the throat closure `T=iσ_y` has `T²=−I`, `½ tr T²=−1` (the Pin⁻ structure forced by `w₁²=w₂`); the grading `tr(T^k)=2cos(kπ/2)=0` for odd k (off-diagonal, fermionic — the realized lepton rungs) and `±2` for even k (diagonal, bosonic). None reference the metric. **Survival:** a smooth metric/frame deformation acts on the holonomy by orientation-preserving conjugation and on the deck map by a GL⁺ frame change — neither can flip a determinant sign or a trace; across 1000 random deformations `½ tr T²` stays `−1` (to 10⁻¹⁵), the brane deck det stays `−1` and the bulk `+1` (to 10⁻¹⁵); named deformations (a Berger squash of the S³ Hopf fiber, a tidal-charge radial stretch) too. **The generation count:** odd k ≤ `k₅=D_bulk=5` ⟹ `{1,3,5}=(k₅+1)/2=3` generations (matching the repo's `LEPTON_BASELINE_DEPTHS`); `D_bulk` (an integer) and the odd-parity selection are topological, so 3 generations survive every smooth deformation — not an artifact of the round metric. **Changes only at a topology change:** the only sector-flipping path (`T²:−I→+I`, non-orientable→orientable, odd→even, fermion→boson) is a NON-metric path `T(θ)=exp(iθσ_y)` whose orientability invariant crosses 0 at θ=π/4 — a degenerate spin structure, the topology-change event; smooth deformations act by conjugation and never move θ, so they can never reach it (exactly as continuous evolution keeps `\|q\|>0` in #182). **Robustness:** 500/500 random deformations preserved the odd-k grading and the orientability class. **Unity:** the generation sector is to the bulk what the order-field winding is to the soliton (#181/#182) — a topological charge robust to smooth deformation, changing only at a singular/topology-change event; the #174 round-metric derivation is not special, the 3-generation odd-k structure is topologically protected. **Scope:** the invariance is EXACT (topological: deck determinant + spin-closure / Stiefel–Whitney class, metric-independent); the deformations are within the orientability/spin-preserving class (a genuine topology change is outside the smooth family by construction); this establishes robustness, not a re-derivation of #174; purely topological — weak-field is not invoked (`odd_k_generation_survival_probe`, PR #183) |
| **The phase-slip / topology-change event** — exactly how the invariant changes when q hits zero | **Dissects the one event that changes the discrete charge: the topological OBSTRUCTION (changing Q forces an exact zero — a straight homotopy winding-1→winding-0 is driven through `min\|q\|=2.5×10⁻¹⁷` at `s*=0.5, φ*=π`), the QUANTUM (across a simple zero `ΔQ=±1`, the integrated `∮∇φ` changing by exactly −2π), and the DYNAMICS (Q(t) is a quantized staircase whose every step coincides with a `min\|q\|→0` event); away from the zeros Q is an exact integer (to 10⁻¹⁵) — the invariant is ambiguous ONLY at amplitude zeros** | Resolves the #181 handle (the unsustained winding driven to `\|q\|=0` slips) into the full anatomy of the topology-change event. **The obstruction:** Q is a homotopy invariant of `q:S¹→ℂ∖{0}`, so to change it the field MUST leave `ℂ∖{0}` — the straight homotopy `(1−s)·[winding 1]+s·[winding 0]` is FORCED through an exact zero (`min\|q\|=2.5×10⁻¹⁷` at `s*=0.5`, located at `φ*=π`), where Q jumps 1→0; there is no nowhere-zero path between sectors (the dynamical content of the #175 gate). **The quantum:** across that simple zero `ΔQ=−1` exactly — the integrated winding density `∮∇φ` changes by `−2π` (one full turn of phase removed at the zero point); the winding is an integer before and after, ambiguous only AT the zero. A generic simple zero carries unit topological charge in (space×parameter), so each elementary slip is `±1`. **The single event:** in a genuine ψ–Φ–q evolution the unsustained k=5 holds FLAT at 5, then the first slip steps `ΔQ=−1` to 4 EXACTLY when `min\|q\|(t)→0` (`min\|q\|=5×10⁻⁴` at the step). **The staircase:** a strongly over-wound state (k=8) cascades down a quantized staircase `[8,7,5,4,3,2]`, every step coinciding with a `min\|q\|→0` event — the over-wound throat sheds winding one quantum at a time through amplitude-zero nodes (a recorded `−2` step is two elementary slips unresolved in sampling time, not a `ΔQ=2` event). **Localization:** away from the slips (where `min\|q\|>0.1`) the unrounded winding equals an integer to `10⁻¹⁵` — Q is a rigid integer between events (#181), jumping by `±1` at them (#182); the non-integer-ness is confined to the measure-zero set of amplitude zeros. **Physical meaning:** the phase slip IS the throat changing its winding / generation sector (`k→k∓1`) through the amplitude-zero node — the #175 antipodal node, the #178 defect core; the #175 'gate' is here the sector-CHANGING event itself, resolved. With #181 (survival between events), the throat's winding is a conserved topological charge that transitions ONLY at nodes; the odd-k ladder is the #174 orientability grading, its survival under a deformed bulk geometry the subject of #183. **Scope:** the obstruction and the `±1` quantum are EXACT (topological); the dynamical staircase is on the reduced vortex-on-soliton loop; the full 2D/3D vortex-line reconnection (the zero as a moving vortex line threading the loop) is a follow-up; weak-field (`phase_slip_topology_change_probe`, PR #182) |
| **Discrete invariant survival on the ψ–Φ–q soliton** — the continuous geometry carries the discrete charge | **The throat's discrete winding charge `Q=(1/2π)∮∇φ` (the #178 vortex charge / #174 odd-k ladder) survives on the #180 continuous self-consistent soliton: dress its ordered core with winding-k phase and Q is conserved to MACHINE PRECISION (ΔW~10⁻¹⁶) under continuous ψ–Φ–q evolution while `|q|>0` — rigid under every `|q|>0`-preserving homotopy, changing ONLY through `|q|=0` (the unsustained k=5 phase-slips — the #182 event)** | Bridges the continuous soliton (#179/#180) and the discrete ladder (#174/#178): does the continuous geometry CARRY the discrete invariant? **The setup:** take the #180 self-consistent soliton (M=3.5, actually built) and an equatorial loop of radius `R=0.75` in its ordered core (`ρ_loop=0.36>ρ_c=0.20`, so `|q|>0`); dress it with winding `q=|q|e^{ikφ}`. **Quantized:** `Q=(1/2π)∮∇φ=k` exactly (to ~10⁻¹⁵) for k∈{1,3,5}. A winding-k vortex is SUSTAINED when the well beats the centrifugal cost, `A²=(gρ−a₀)−(κ/R²)k²>0`: the soliton sustains k=1 (A²=0.153), k=3 (A²=0.082); k=5 exceeds it (A²=−0.061). **Survival — wave:** under continuous norm-conserving (unitary) wave evolution Q is conserved to machine precision (ΔW~10⁻¹⁶) for ALL k∈{1,3,5}, `min|q|>0` throughout — the charge rides the continuous dynamics untouched. **Survival — dissipative:** under the order field's OWN gradient-flow dynamics (the same relaxation that built the soliton), the sustained windings k=1,3 survive — a perturbed vortex relaxes back to the clean one, Q conserved to ~10⁻¹⁵, `min|q|>0`. **The criterion (→#182):** Q changes ONLY through `|q|=0` — the unsustained k=5 is driven to a zero (`min|q|→10⁻⁴`) and the charge SLIPS (5→2); survival ⟺ `|q|>0`, exactly. That slip is the #182 topology-change event. **Rigidity:** under 40 random `|q|>0`-preserving homotopies per sector, Q is unchanged in every case (40/40 for k=1,3,5) — a superselection charge outside the continuous moduli (the #173/#174 rigidity, now on the dynamical soliton). **Honest scope:** homotopy-invariance is exact; the geometry is the reduced vortex-on-soliton (amplitude from the #180 radial soliton, winding azimuthal — the full 2D/3D vortex-LINE soliton with the winding back-reacting on ψ,Φ is a follow-up); which rungs survive is set by the soliton's capacity (k=1,3 dissipatively, all three under wave evolution); the realized PHYSICAL ladder is odd-k {1,3,5} by the #174 orientability grading (even-k conserved too but excluded by orientability), with its survival under a DEFORMED bulk geometry the subject of #183; weak-field (`discrete_invariant_survival_probe`, PR #181) |
| **ψ–Φ–q soliton hardening** — stationarity, branch scan, and basin map | **Hardens #179's two-way throat-soliton into a trustworthy object AND corrects it: with a spectral ψ kinetic the soliton is a genuine STATIONARY eigenstate (residual ~10⁻⁴; ψ stays stationary in the frozen self-consistent (Φ,q) background, mass-conserving), a SMOOTH everywhere-convergent BRANCH in mass and self-gravity over the tested range μ∈[0.05,2], and a robust ATTRACTOR (the full width×seed initial grid → the same state, ~1%); #179's reported high-μ "runaway collapse" is shown to be a finite-difference DISCRETIZATION ARTIFACT — the spectral ψ kinetic finds no collapse up to μ=2, only a branch deepening out of weak-field validity** | The hardening probe for #179 (as #177 was for #176), with the three things a self-consistent PDE state needs plus a re-examination of the collapse claim. **Stationarity:** put ψ's kinetic on a spectral basis (`u=rψ`, DST; the order field `q` keeps its finite-difference Laplacian — this is **not** a fully spectral ψ–q solver) so the imaginary-time relaxation and the real-time step share the SAME ψ Laplacian — the relaxed state is then a genuine eigenstate (`‖Hψ−μψ‖/‖ψ‖≈10⁻⁴`, chemical potential `μ≈−1.45`), and evolving **ψ alone in the frozen self-consistent (Φ,q) background** by a unitary real-time split-step leaves it stationary (profile drift `~4×10⁻⁵`, mass conserved to machine precision): ψ is a stationary eigenstate of its self-consistent potential, a real bound soliton (the fully coupled real-time ψ–Φ–q dynamics is a follow-up). **Branch scan (mass):** the soliton is a smooth monotone family — the order field switches on where `ρ_peak` crosses `ρ_c=0.20` (near M≈2.7, the spatial onset just above `ρ_c` by the GL droplet barrier), `max\|q\|` and the well deepening monotonically. **Branch scan (self-gravity μ):** smooth, monotone, and EVERYWHERE-CONVERGENT over the tested range — `max\|q\|` 0.42→2.62 and `Φ(0)` −3.09→−24.6 across `μ∈[0.05,2]` with residuals `≤10⁻³`, no blow-up. **The #179 correction:** #179 reported a runaway collapse at super-critical μ (`\|q\|→31`, `Φ(0)→−252`) — but that used a finite-difference (`np.gradient`) Laplacian; the spectral ψ kinetic finds NO collapse up to μ=2, so the runaway was a DISCRETIZATION ARTIFACT. The genuine large-μ limit is not a runaway but the soliton deepening out of weak-field validity (`Φ(0)`: −3.09→−24.6 across the tested μ) — the strong-field domain for full NR. **Basin map:** the soliton is a robust attractor — the FULL initial-condition grid (Gaussian widths `w∈{1.2,1.8,2.6}` crossed with order seeds `∈{10⁻²,10⁻¹}`, all six combinations) flows to the SAME state (`max\|q\|` spread **~1%**, `Φ(0)` spread **~0.1%**); a tiny seed `10⁻³` reaches the same attractor, only more slowly (convergence-time, not a different basin). **Robustness:** the well depth `Φ(0)` converges under grid refinement (`N=160→240→320`: −3.34→−3.09→−2.98, successive change ~3.6% at the finest step), while the pointwise core `max\|q\|` is more grid-sensitive (the sharp core, ~10% per refinement — converging but slowly), so the core amplitude carries grid uncertainty while the soliton's existence and structure are robust. **What survives #179:** the soliton's existence, two-way back-reaction, and threshold continuity — confirmed and hardened; the specific "runaway" claim does not survive as stated. **Scope:** weak-field/semi-dynamical, spherically reduced, effective constants; the deep-μ branch and the strong-field endpoint are for full NR (`psi_phi_q_soliton_hardening_probe`, PR #180) |
| **Two-way ψ–Φ–q evolution** — the self-consistent matter–metric–order system | **Closes #178's one-way coupling into the full two-way system of three co-evolving fields (matter wave ψ, gravitational potential Φ, throat-order field q) from ONE energy functional: the coupled flow CONVERGES to a self-consistent throat-soliton whose well is ~5% deeper and core ~13% denser than the pure self-gravitating state (the throat traps the wave AND gravitates); the binding feedback saturates into a stable bound object while q's self-gravity drives runaway collapse above a critical coupling; sub-threshold it reduces exactly to Schrödinger–Newton** | Closes the dynamical-coupling follow-up #178 flagged (where `ρ→q` was one-way: q neither gravitated nor acted on the wave). **One functional:** the whole system descends from `E[ψ,q]=∫[½\|∇ψ\|²+½κ\|∇q\|²+½a₀q²+¼λq⁴−½g\|ψ\|²q²]+W_grav[\|ψ\|²+μq²]`, whose fixed-mass gradient flow is `∂_τψ=½∇²ψ−Φψ+½gq²ψ`, `∂_τq=κ∇²q−(a₀−g\|ψ\|²)q−λq³`, `∇²Φ=4πG(\|ψ\|²+μq²)` — so the FOUR channels (ψ↔Φ Schrödinger–Newton; ψ→q ordering; q→ψ the throat binds the wave; q→Φ q gravitates) are consistently coupled (the ordering and binding terms share the same `g`). **Self-consistent:** the coupled flow converges — energy monotone→plateau, q stationarity residual `~10⁻⁴` — a self-consistent throat-soliton EXISTS. **Two-way back-reaction:** at M=3 (super-threshold) the order field nucleates (`max\|q\|≈0.49`) and vs the pure Schrödinger–Newton soliton (q=0) the self-consistent state has a deeper well (`Φ(0)=−3.18` vs `−3.03`, **5.0% deeper**) and denser core (**13.4% denser**) — the throat traps the wave, concentrating it, strengthening the order (a self-reinforcing loop). **Saturation vs collapse:** with sub-critical self-gravity (`μ=0.05`) the quartic `λq⁴` saturates the feedback into a STABLE bound soliton (`\|q\|` plateaus, residual `10⁻⁴`; intermediate μ just gives a denser soliton), but super-critical self-gravity (`μ=2.0`) has NO weak-field fixed point — the flow diverges (`max\|q\|→31`, `Φ(0)→−252`, residual `→6500`): the onset of strong-field gravitational collapse. **Continuity:** below the ordering threshold (M=1, `ρ_peak≈0.005<ρ_c=0.20`) the order field vanishes and the system reduces EXACTLY to the Schrödinger–Newton soliton of #176/#177; the #176→#178→#179 arc is one continuous system, switched by the matter concentration. **Honest scope:** weak-field/semi-dynamical and spherically reduced (self-gravity sphericalizes, #176); the constants `a₀,g,λ,κ,μ` are EFFECTIVE (the result is the STRUCTURE — a self-consistent two-way throat-soliton with a stable binding branch and a runaway self-gravity branch — not the numbers, whose microscopic values await `V(q)` and the q–metric coupling from the 5D bulk); the stable soliton is sub-critical and the strong-field runaway endpoint is for full NR (`two_way_psi_phi_q_probe`, PR #179) |
| **Self-gravity-driven throat-order instability** — bind, or drive a new order parameter? | **Weak-field concentration does NOT merely bind: coupling the #176/#177 self-gravity solver to the #178 order field through a density-dependent Landau potential, the matter density is the control parameter of a geometric phase transition — bind-only below a critical concentration ρ_c, order-NUCLEATING above it (reached only by the gravitational collapse), and the ordering vanishes with gravity off** | Answers the dynamical-coupling follow-up the #178 throat-order field flagged: does the self-gravitating concentration of #176/#177 merely BIND the wave, or DRIVE the order parameter? **The coupling:** the matter density `ρ=\|ψ\|²` (from the #176/#177 solver, actually run) is the control field of a density-dependent Landau potential `V(q;ρ)=½(a₀−gρ)\|q\|²+(λ/4)\|q\|⁴`, so the order field's effective mass² `a(ρ)=a₀−gρ` changes sign at a critical concentration `ρ_c=a₀/g` (here 0.30): below it `q=0` is the only minimum (disordered, merely bound), above it `q=0` destabilizes and the order parameter rolls to `\|q\|=√((gρ−a₀)/λ)` (ordered). **Merely bound:** a sub-threshold packet (M=1, G=1) concentrates only to `ρ_peak≈0.06<ρ_c`, so the order field relaxed under the Ginzburg–Landau gradient flow stays at zero (`max\|q\|~10⁻¹⁶`) — bound, but NO geometric order. **Drives order:** above the mass threshold (M=3, G=1) the collapse drives `ρ_peak≈0.90>ρ_c` and the order field NUCLEATES a localized symmetry-broken domain (`max\|q\|≈0.68`) sitting exactly at the density peak (the throat core of #178). **It is gravity:** with gravity off (`G=0`) the same M=3 packet never concentrates past `ρ_c` (`ρ_peak≈0.18`) and no order nucleates; restoring `G` it does — the ordering inherits the `M_c∝1/G` gravity of #176/#177. **Dynamical:** driving `q` by the time-dependent `ρ_peak(t)` of the collapse, the order parameter switches on (causally) only after the density crosses `ρ_c` — a moving order front following the gravitational concentration. **Honest scope:** ONE-WAY coupling (`ρ→q`; `q`'s back-reaction on the metric — the self-consistent q–metric system — is the next step); the constants `a₀,g,λ` and hence the numerical `ρ_c` are EFFECTIVE (the result is the EXISTENCE of a gravitationally-crossed concentration threshold separating bind-only from drive-order, not the specific `ρ_c`, whose microscopic value awaits `V(q)` from the 5D bulk); the spatial nucleation carries the usual GL droplet-size barrier (the ordered core must exceed the coherence length `ξ=√(κ/\|a\|)`); still weak-field/semi-dynamical. The collapse of #176/#177 carries the matter density across the order transition and nucleates the throat-order field of #178 (`self_gravity_driven_order_probe`, PR #178) |
| **The throat-order field q(t,r,θ)** — the throat as a topological defect of a Ginzburg–Landau order parameter | **Introduces ONE complex order field whose defects ARE the throats, unifying the arc: the defect core (`\|q\|→0`) is the antipodal node (#175), the phase winding (`∮∇φ=2πk`) is the discrete odd-k charge (#174), and the disorder→defect nucleation is the focusing threshold (#176/#177)** | Introduces the throat-order field `q(t,r,θ)=\|q\|e^{iφ}` — a single complex Ginzburg–Landau order parameter — to unify three discrete facts the arc established in three separate languages. **The two phases:** the Mexican-hat Landau potential `V(q)=(λ/4)(\|q\|²−q₀²)²` has the disordered point `q=0` as an UNSTABLE local maximum (`V″=−1<0`, the symmetric phase) and the ordered amplitude `\|q\|=q₀` as a STABLE degenerate minimum (`V″=+2>0`, the broken-symmetry vacuum that fills the orientable bulk; the free phase φ is the U(1) a defect winds). **The throat is a vortex:** the radial GL profile `f(r)=\|q\|(r)` solving `f″+f′/r−k²f/r²=λf(f²−q₀²)` with `f(0)=0, f(∞)=q₀` exists for each winding — `\|q\|=0` at the core, healing to `q₀` in the bulk, the core widening with k (core size r(`\|q\|=q₀/2`) = 0.8/1.4/2.0 for k=1,3,5): a localized defect, the throat itself. **Winding = the discrete k:** the topological charge `∮∇φ/2π` is the integer winding (measured 1,3,5), conserved while `\|q\|>0` (`π₁(S¹)=ℤ`); the realized sector is ODD-k — the #174 orientability grading. **Core = the antipodal node:** the field must vanish where the phase winds, so the defect core `\|q\|=0` IS the forced amplitude-zero node of #175 — acquiring the discrete charge from the continuous (winding-0) sector REQUIRES passing through a zero (the gate is topological, not dynamical). **Nucleation = the threshold:** the disordered `q=0` is unstable, so under the GL gradient flow `∂_t\|q\|=−λ\|q\|(\|q\|²−q₀²)` any perturbed region rolls off zero to `q₀` and a fixed-winding defect nucleates; the trigger that drives a region off zero is precisely the self-gravitating focusing of #176/#177 (`M_c∝1/G`). **Honest scope:** the EFFECTIVE Ginzburg–Landau level — q is introduced as the coarse-grained order field whose defects are the throats; the microscopic `V(q)` (λ, q₀ from the 5D bulk action) and the dynamical q–metric coupling are follow-ups. The throat's three discrete facts are now one object: a vortex of `q(t,r,θ)` (`throat_order_field_probe`, PR #178) |
| **Weak-field self-gravity threshold hardening** — controls, scaling, robustness | **Turns #176's promising self-gravity proxy into a trustworthy PDE benchmark: G=0/repulsive controls give no collapse (it is gravity); the threshold coincides with the energy binding criterion E=T+W=0; M_bind·G=const to 0.69% (the 1/G law, sharp); grid-converged, mass-conserving** | Hardens #176 with the three things a credible PDE result needs. **Controls:** with `G=0` (gravity off) the packet never concentrates at any mass (M=1,3,5 all disperse), and with `G<0` (repulsive) even M=5 does not collapse — the threshold requires attractive gravity, not an artifact of the packet/grid. **Energy anchor:** the total energy `E=T+W` (kinetic + gravitational self-energy) defines the rigorous binding threshold `M_bind` (E=0); dynamically, below it (E>0, unbound) the core mass drains (disperse), above it (E<0, bound) the core holds — the dynamical disperse/bound transition tracks the energy sign (an independent physics check on the integrator). **Scaling, sharp:** `M_bind·G` is constant across `G∈{0.5,1,2}` — `1.134, 1.133, 1.127`, a spread of only **0.69%** (versus #176's coarse 0.48) — so `M_bind ∝ 1/G` to <1%; the full Schrödinger–Newton invariant `M_bind·G·w ≈ const` holds across widths (~10%, the residual from the Gaussian not being the exact eigenstate). **Robustness:** `M_bind` converges to ~1–2% under radial-grid refinement (N=120→160→220: 1.119→1.133→1.144) and the split-step integrator conserves mass to ~10⁻³ (the energy drift is dt-independent — a diagnostic, not an instability). **Now trustworthy:** gravitational (controls), energy-validated, 1/G to <1%, SN-scaling, grid-converged. **Standing scope:** still weak-field/semi-dynamical (not full NR; the strong-field horizon/throat endpoint is for full numerical relativity); the fixed-width Gaussian makes the `M·G·w` invariant approximate (`self_gravity_threshold_hardening_probe`, PR #177) |
| **Real GR back-reaction** — a semi-dynamical axisymmetric self-gravitating wave packet | **Past the 1D proxy: under a metric that responds to the energy density, the disperse/collapse threshold of the focusing arc survives — and the critical mass scales as 1/G, proving it is GRAVITY, not a toy nonlinearity** | Moves past #175's 1D ring proxy (whose focusing was an ad-hoc `g\|ψ\|^p`) to ask whether real general relativity backs it up. **The model:** a semi-dynamical, axisymmetric self-gravitating scalar `i∂_tψ=−½∇²ψ+Φψ`, `∇²Φ=4πG\|ψ\|²` (the weak-field Einstein–Klein–Gordon / Schrödinger–Newton system; the metric `g_tt=−(1+2Φ)` responds to the energy density `ρ=\|ψ\|²`) — `ψ(r,θ,t)` in the `(r,ℓ)` Legendre basis, radial Laplacian by Dirichlet sine transform, centrifugal `ℓ(ℓ+1)/2r²` diagonal in ℓ, and `Φ(r,θ)` from the axisymmetric multipole Poisson each step; split-step, mass-conserving. **Stable:** mass conserved to ~10⁻³. **The gravitational threshold:** below a critical mass the packet disperses (the metric stays shallow), above it the self-gravity concentrates it (the metric deepens, runaway) — the disperse/persist threshold of #58/#166/#175, now under actual gravitational back-reaction rather than a phenomenological nonlinearity. **It is gravity (the decisive test):** the critical mass (interpolated from the peak-growth crossing) scales as `1/G` — `G=0.5→>3.2`, `G=1→2.29`, `G=2→1.10`, halving from G=1 to G=2 (ratio 0.48≈0.5). The threshold is set by the gravitational coupling, not a toy knob. **The metric responds:** the central potential well deepens as the field concentrates (the back-reaction in action). **Honest scope:** semi-dynamical weak-field GR — the field evolves while the metric responds quasi-statically (not full NR); the collapse confirms the threshold and the concentration (the throat-formation analog), but the strong-field endpoint (a horizon / resolved throat) is for full numerical relativity; self-gravity sphericalizes (the monopole dominates), so the collapse is predominantly radial — the axisymmetric machinery is exercised, not a directional jet claimed (`self_gravitating_axisymmetric_probe`, PR #176) |
| **The nonlinear antipodal focusing PDE sandbox** — can a continuous geometry evolve into the discrete sector? | **Yes, but only through the caustic: smooth evolution conserves the winding (the discrete sector is locked out, #174 dynamically); the gate is an amplitude-zero node forced exactly at the antipode; the nonlinear focusing reaches it only above a critical mass (#58/#166); the jump is quantized ±1** | The dynamical capstone — #166 simulated the *linear* focusing and deferred the nonlinear throat formation; #174 showed the discrete sector sits outside the continuous deformation manifold. **The sandbox:** a complex field `ψ(χ,t)` on the antipodal ring evolved by a focusing NLS `i∂_tψ=−∂_χχψ−g\|ψ\|^pψ` (split-step Fourier, mass-conserving); the discrete sector = the winding `Q=(1/2π)∮d(arg ψ)`, an integer proxy for the throat's `k`. **Smooth conserves:** a smooth `Q=1` field (`\|ψ\|>0`) keeps its winding EXACTLY (1.00000→1.00000, mass conserved) — a continuous geometry cannot smoothly deform into a different discrete sector (the dynamical confirmation of #174). **The gate, at the antipode:** the winding is a homotopy invariant of maps to `ℂ∖{0}`, so Q changes only across an amplitude-zero node; interpolating between `Q=0` and `Q=1` forces that node EXACTLY at `χ=π` (the focus) — so the discrete sector is gated by a singular core at the antipodal caustic, and the focusing is what drives the field there. **The threshold:** with the critical (quintic) nonlinearity, mass <~1.6 disperses (stays continuous, Q frozen) while mass >~1.9 concentrates toward the core (the nucleation onset) — the disperse/persist threshold of #58/#166, now actually simulated nonlinearly (which #166 deferred). **The jump is quantized:** crossing the antipodal node changes Q by exactly ±1 — a discrete response to a smooth focusing drive. **Answer:** a continuous time-dependent geometry enters the discrete sector ONLY by developing a focusing singularity at the antipode, never by smooth deformation. **Honest scope:** a reduced 1D ring model (Q proxies the discrete k, the collapse core proxies throat nucleation, the critical-NLS collapse is marginal); the conceptual answer is robust, the numbers model-dependent (`nonlinear_antipodal_focusing_pde_probe`, PR #175) |
| **The odd-k ladder** — forced, rigid, unique to the non-orientable 5D geometry | **The cleanest discrete feature audited: the odd-k lepton ladder {1,3,5} is forced (fermion ⟹ odd via the T²=−I orientability grading; k≤5 ⟹ 3 generations), rigid against every active/null/mixed continuous deformation (the rank story survives nonlinearity), and a unique signature of the non-orientable antipodal 5D geometry** | Continues the #173 inverse problem on one discrete feature, asking whether the geometry forces it and whether anything but this geometry could produce it. **The origin** (recap #67/#169/#170): the throat monodromy `T=iσ_y` (`T²=−I`) makes `T^k` off-diagonal for odd k (the orientation-reversing, non-orientable RP²/Pin⁻ closure of a spin-½ fermion) and diagonal for even k (orientable, bosonic) — so `k mod 2` is the orientability grading; charged leptons are fermions ⟹ odd k, and `k≤k_5=D_bulk=5` ⟹ {1,3,5}=(k_5+1)/2=3 generations. **The three direction sets** (from the #173 Jacobian): *active* (the rank-10 subspace) moves the masses/CKM **linearly** (exponent 1.03); *null* (the 10-dim kernel) is **flat to first order** (exponent 2.0, ~10⁴× smaller response than active at ε=10⁻²); *mixed* (active+null) is **active-dominated** (1.02) — so **nonlinear effects do not break the local rank story** (the null leakage stays quadratic, not linear). **Forced & rigid:** the odd-k labels and the generation count are integer winding + the ℤ₂ grading (`T²=−I`), discrete topological data *outside* the entire continuous deformation manifold (no generation-number knob) — not an emergent near-integer that could drift, but structurally forced; the continuous geometry deforms only the masses and the CKM, the ladder is rigid against all of it. **Unique:** an orientable geometry (`T²=+I`) gives the even/bosonic sector, not an odd-only fermion ladder; the specific {1,3,5} needs `k≤5=D_bulk` — so odd-{1,3,5} is the joint signature of the non-orientable antipodal spin structure and the 5D bulk (an exclusion/signature argument within BAM, not a no-go against all alternatives) (`odd_k_ladder_rigidity_probe`, PR #174) |
| **The sensitivity audit** — Jacobian rank, the forced core, the isolation dimension | **The dynamical inverse problem, measured: of 14 live observables the free geometry dials rank 10; the FORCED CORE is 4 (entirely CKM unitarity, structural V=U₊†U₋); the masses are FITTED; a 10-dim compensator redundancy (the diagonal shifts); deriving φ_h saves 0 effective inputs** | Lets continuous geometry search the spectrum/signature landscape instead of adding more static proofs: compute the Jacobian `J_ij=∂O_i/∂I_j` of the live observables at the lock and read its SVD. **Observables** (14): 4 quark mass ratios (s,c,b,t/d), the 5 CKM magnitudes, J, β, γ (`LOCKED_QUARK_PARAMS_V4`), and the 2 lepton ratios (μ/e, τ/e). **Inputs:** the free *fitted* knobs — NOT the k₅-derived locks (φ_h=π/k₅, χ=k₅(k₅−1), uplift=1−1/k₅², action=π, winding=max), which are zero-cost. **Isolation dimension** rank(J)=10 (quark 8 + lepton 2), with a clean singular-value gap (quark sv falls 22.6→1.2e-2 then drops ~5e-6 to zero). **Forced core** = n_obs−rank = **4**, entirely CKM combinations (forced weight on CKM >0.99, on masses ~1e-16) — the **CKM unitarity relations**: `V=U₊†U₋` is exactly unitary (‖V†V−I‖~1e-16), so the 8 CKM observables on the 4-parameter unitary manifold force 8−4=4 relations. The largest set forced at zero input cost — a genuine structural prediction, but the *standard* unitarity, not a BAM-specific relation. **The masses are FITTED** (zero weight in the forced core, quark and lepton): the ladder sets their values, the knobs span them; no forced mass relation. **Compensator redundancy** = n_inputs−rank = **10**, dominated by the mass-preserving diagonal shifts (kernel share ~0.8) — the *loose knob / compensator* structure the program flagged qualitatively, now measured; the v4 quark parametrization is substantially over-complete. **CP-at-zero-cost test:** adding φ_h as an input leaves the rank unchanged (8→8) — the CP-phase direction is already spanned, so deriving φ_h saves no effective input; *"CP at zero parameters"* is a counting economy, not a Jacobian reduction. **Robustness:** the rank is stable across a 4×3 grid of finite-difference step (1e-3…1e-6) and SVD tolerance (1e-4…1e-8) — total rank 10 in 11/12 cells (the lone outlier is the largest step at the tightest tolerance, truncation noise crossing the cut), with the rank-8 separation gap ≥1e3 everywhere, so the forced core of 4 is not a tuning of h or tol. An audit — a measurement of the predictive content, honest where it is not flattering (`sensitivity_jacobian_audit_probe`, PR #173) |
| **The measured Fermi equation of state** — a many-throat ensemble (companion to #171, same branch) | **The EoS *measured*, not assumed: free fermions in a box with the −1 exchange sign as Pauli single-occupancy give Γ=5/3 (NR, measured 1.6665) and 4/3 (UR, 1.3332) and P/u=2/3, 1/3; the Bose control gives Γ=1 with vanishing degeneracy pressure — three routes (assumed #170 / topological #171 / measured #172) agree** | The second of the two closing options, kept on the same branch as #171 to compare. #170 *assumed* antisymmetry and read 5/3 off the analytic Fermi integral; #171 derived the −1 sign topologically; here the −1 is imposed as **Pauli single-occupancy** in a many-throat box ensemble and the EoS is **measured** from the simulated ground-state energies. **The ensemble:** N identical throats as free fermions in a cubic box; each spatial mode `(n_x,n_y,n_z)` holds `g=2`; the many-body ground state is the filled Fermi sea built by filling the N lowest of ~1.3M enumerated modes — no occupation distribution assumed. **Measured P/u:** from `P=−dE/dV`, `0.6667=2/3` (NR) and `0.3333=1/3` (UR) — the virial relation emerges. **Measured Γ:** the local log-slope of the filled-mode energy sum `K(N)`, finite-size-extrapolated via the Weyl correction `Γ(N)=Γ∞−a·N^{−1/3}`, gives **1.6665≈5/3** (NR) and **1.3332≈4/3** (UR), 0.01% from target — outputs of the simulation, not a formula. **The control:** Bose statistics (all N in the ground mode) give `Γ=1` and a `T=0` degeneracy pressure that vanishes (mode energy `∝1/L²→0`) — so the 5/3 stiffening is a *measured consequence of the −1 exchange sign*, not a property of the box. **Comparison:** the measured indices reproduce #170's assumed values and confirm the EoS implied by #171's topological sign — three routes, one answer. **Honest scope:** the exchange sign is the input (Pin⁻/geon, #170/#171, as Pauli occupancy), not re-derived here; idealizations are the standard degenerate-gas ones (free, T=0, box) (`measured_fermi_eos_ensemble_probe`, PR #172) |
| **π₁ of the two-mouth configuration space** — FR-homotopy survival for the Pin⁻ throat (geon statistics) | **Replaces #170's *orientable* FR citation with the geon-statistics framework: π₁ exchange σ²=e (no 3D braiding), spinorial 2π=−1, Pin⁻ reflection²=−1 (the ingredient orientable FR lacks), achiral → the −1 (Fermi) survives, conditional on the cited Dowker–Sorkin exchangeability hypothesis** | Closes the cited gap in #170: the Finkelstein–Rubinstein homotopy it invoked is the **orientable** result, but the throat mouth is the non-orientable RP² (Pin⁻), so the correct framework is **geon statistics** (Friedman–Sorkin; Aneziris–Balachandran et al.; Dowker–Sorkin), where spin–statistics is a theorem *with hypotheses, known to fail for some geons*. **π₁ of the two-mouth config space** (model): exchange `σ` with `σ²=e` (≥3D ⇒ symmetric group, **no braiding**, only ±1 statistics), per-geon 2π rotation `R_i`, and — because the mouth is non-orientable — an orientation-reversing loop `τ_i`. **Spinorial:** the single geon's 2π rotation = −I (4π=+I), the Pin⁻ holonomy of #170 / Friedman–Sorkin's spin-½ from gravity. **The new ingredient orientable FR lacks:** the non-orientable exchange carries an orientation reversal — a *reflection* — and RP² admits **Pin⁻ only** (`w₂≠0` kills Pin⁺), in which a reflection **squares to −1** (computed; Pin⁺ gives +1) — exactly what makes the non-orientable exchange sign well-defined and fermionic. **Achirality:** a non-orientable geon is its own mirror image, meeting the geon spin–statistics theorem's handedness hypothesis automatically. **Result:** spinorial + Pin⁻ reflection²=−1 + achiral ⇒ exchange = −1 (**Fermi survives**), now on the correct non-orientable footing. **Honest scope:** conditional on the Dowker–Sorkin exchangeability ("slide") hypothesis — holds for identical asymptotically-flat throats, *cited not derived* from the BAM field theory; the literature has spin–statistics violation examples, so this is a genuine check BAM passes (Pin⁻ + achiral + exchangeable), not automatic. The remaining gap is that hypothesis (and the field-theory mapping class group, modeled here), not the spinor sign or the reflection algebra (`geon_statistics_pi1_probe`, PR #171) |
| **The Z₂ link twist** — non-orientable ER links as a *gauge field* over the network | **The orientation Z₂ is a gauge field, not a parity: it freezes mouth exchange EXACTLY, halves the quotient revival, and IS the repo's own Möbius sector — but leaves the nucleation threshold untouched** | Joins two arcs that never met: the antipodal wave interaction (#128/#129→#135→#137/#138→#166/#175), which treats the antipodal Z₂ `(−1)^l` as a **single global parity of one throat** (#139's Thread A), and the bulk ER-link network (#206–#224), where orientation is an **inert decoration** (`NetworkMouth.orientation`, multiplied into `U_BA` and never varied). A Z₂ on a *network* is not a parity but a **gauge field**: link variables `ε_b`, gauge orbits `ε_b → s_a ε_b s_a'`, Wilson loops `W(γ)=∏ε_b`, only holonomies physical. **THE HEADLINE:** a twisted link freezes the #224 mouth-exchange doublet — `Δω` ratio `2e-11…5e-11` over `r_s=0.25–0.5`, `T_exchange` 1833 → 3.8e+13 — and it is a **theorem**, not a two-level effect: the half-ring translation obeys `S²=W`, so `W=−1` gives eigenvalues `±i` and a **real** ring operator forces *every* level into an exactly degenerate pair (7.1e-12 vs 2.1e-03). This is the complement of #224's own two freezing mechanisms, both approximate and both degrading toward symmetry — the twist is exact precisely there, closing the "exact-degeneracy loophole" #224 names. **Gauge invariance exact** (0.00e+00 drift across gauge copies on the real 1114-site ring). **The Möbius sector was already in the repo**: `make_glueball_ring` (periodic,+1) and `make_mobius_tube` (mobius,−1) are the two holonomy sectors of one cycle, built independently in the QCD arc. **On the quotient the antipodal revival period HALVES to `πR`** with the twist class setting its sign (`−f₀` untwisted, `+f₀` twisted, both 0.00e+00). **Corrected by measurement, not defended:** the conjectured factor-of-4 energy advantage at the focus is **wrong** (energy conserved exactly, ratio 1.000000; sector caustic gains differ by exactly `2/l_max` — so the `2 m_e c²` threshold of #58/#166 is **untouched**), and the advertised "new vertex rule" `Σl ≡ n_twist (mod 2)` is a **tautology** (`ε=(−1)^l` is forced), leaving the reinterpretation of #137 as bundle descent plus a genuine *network* junction rule `∏ε_b=+1`. **Precision note:** RP³ is **orientable** (deck det +1, #183) — the twist is the flat Z₂ bundle in `H¹(RP³;Z₂)=Z₂` and the Pin⁻ RP² mouth (#169/#170), *not* ambient non-orientability. The `b₁` audit is the declared falsifier and would fire on a bridges-only reading (#207's perfect matching is a forest, `b₁=0`); it does not, because the physical graph includes the shared S³ exterior arcs (`nonorientable_er_link_z2_probe`, PR #227) |
| **The surgery Z₂ term** — does bridge surgery feel the twist? | **NO — the composition law is unchanged, and the reason is algebraic** | #227 §8.1 conjectured `φ₁₄ = φ_a+φ_b+φ_c + π·[n_twist mod 2]` — surgery across an odd number of twisted bridges swapping the Bell outcome. Computed on #207's own 4-mouth lattice over all 8 sign assignments at every junction holonomy: outer-bridge twists move the pair phase by at most **2.3e-13** (identically zero, not `π`), and every odd-`n_twist` row misses `π` by ≥ **3.09 rad**. **The corrected law is the original law:** `φ₁₄ = φ_a+φ_b+φ_c`. Not a null result awaiting a better lattice — #207 extracts the phase as `arg(c_mp/c_pm)`, a **ratio of two branches traversing the same links**, so a link-uniform sign cancels identically. Not gaugeability either: the twist moves the spectrum by ~2.5e-02. **Where it does act:** the junction's 4.68e-02 residual is bridge-vs-exterior interference around the M2–M3 **cycle** — it dies with cycle length (2.9e-01 → 1.8e-13 as the gap goes 4→48) and scales continuously with the exterior hopping (2.5e-05 → 1.8e-01). A topological phase does neither. So the twist acts **only through cycles**, confirming #227's `b₁` principle *dynamically* (`er_link_twist_surgery_probe`, PR #228) |
| **Is the Möbius identification quantitative?** — and a correction to #100 | **Quantitative in its topological half only; #100's Möbius intercept is wrong** | #227 showed the link twist and the QCD Möbius label are the same Z₂ *structurally*. The `+πσ` tower shift factorizes as `(1/2) × (2πσ)`. **The `½` is derived, shared and topological:** #100's own T4 only *asserts* it (returns `pass: True` with no computation); computed here three independent ways — the #227 twisted ring, the repo's own `MobiusSpectrum`, a direct antiperiodic solve — all exactly ½, with deviation falling 1.5e-05 → 2.2e-07 over grid 120→960 in successive ratios **4.00, 4.00, 4.00** (exactly `O(1/M²)`, finite-difference dispersion and nothing else) and moving by only **1.4e-11** across a factor-of-5 change in circumference. **The `2πσ` is not supplied by the network** — it is the closed-string Regge slope #100 imports, with `σ` the lattice anchor. **And #100's tower carries an unjustified constant:** it gives both towers the *same* intercept `M₀²`, but antiperiodic moding also shifts the zero point — `E₀` periodic `= −π/(6C)`, antiperiodic `= +π/(12C)`, so `ΔE₀ = +π/(4C)` exactly (ζ-regularization and an exponential-cutoff extraction agreeing to 4.1e-07), scaling as `1/C` and therefore not absorbable into a constant. **The Möbius ground state is not `M₀² + πσ`.** The corrected intercept needs the closed-loop dynamics #100 defers, so the gap is reported, not patched. **Scope:** glueballs only (unobserved); the heavy Möbius *baryon* search table (#103/#109/#114 — Λ_c 3135, Λ_b 6469, the 849 MeV dipion endpoint) is built from `Δ=2√σ` and **stands as published** (`mobius_identification_quantitative_probe`, PR #229) |
| **What selects the twist sector?** — energetics select, topology freezes | **The two sectors differ by `π/(4C)` per d.o.f. and the sector is then topologically frozen at nucleation — but energetics favour UNTWISTED in *both* sectors, and matter is selected by the zero-mode INDEX, not by energy** | The last open item of the #227 arc: nothing said *which* networks carry `W=−1`. **Selection is energetic but not relaxational.** The holonomy sectors differ in zero-point energy by exactly `π/(4C)` per degree of freedom (verified to 1e-12), **and the sign was originally read off the statistics — which is wrong.** Computing BAM's actual Pin⁻ spinor (rather than asserting its statistics) shows the holonomy composes with the **intrinsic B2 spin structure** (`T²=−I`, `η=−1`), which **swaps the spinor's moding** relative to the scalar; the statistics flip and the spin-structure flip then **cancel**, so `ΔE=+π/(4C)` for the spinor too and **energetics favour untwisted universally**. **But it cannot relax:** changing `W` continuously drives the link amplitude through **zero**, severing the cycle (`b₁`: 1→0) and leaving the holonomy undefined — the same amplitude-zero gate #175 found for the winding number. And `W` is exactly deformation-invariant: across 200 random re-weightings the untwisted ring keeps an exact zero mode (≤1.6e-15) and the twisted ring never has one (≥4.8e-02), a gap that never closes. The twisted sector does carry a zero eigenvalue (`dim ker D = 1`, `det D = 0`, vs a `π/C` gap and `det D = 2` untwisted), and a first draft read that as #195's Atiyah–Singer mechanism seen from the network side. **That reading is withdrawn on audit:** an index-protected mode survives deformation, and this one does not — a mass lifts it exactly (`lowest = m²`), any potential lifts it (`0.000 → 1.949`), the mode **is** the constant function (overlap `1.000000000000`), and the **index is 0, not 1** (the 1D spectrum is symmetric about zero, so `n₊ = n₋`). `dim ker` is not an index; #195's modes are protected because they live on S² with a monopole `q = k/2` giving `2q = k` modes *of one chirality*, and the 1D cycle has no monopole, no chirality grading, and no index theorem. **Net:** energy explains the bosonic QCD Möbius sector (untwisted cheaper ⟹ excitations above the ground, as #100/#103 present them); **why matter sits in the twisted sector is reopened** — two candidate explanations have now been tried and both failed. Also honest: this is the *free* cycle mode (superpotential `W=0`), not BAM's throat Dirac operator with its Tangherlini superpotential, and the SUSY square root discards the spectrum's sign so it cannot express an index even in principle (`bam_spinor_spectrum_effective_action_probe`, same PR). **Does not explain:** absolute sector populations (set by the nucleation ensemble, not computed), or why a given network nucleated as it did — selection is a preference, not a mechanism (`twist_sector_selection_probe`, PR #230) |
| **The RP³ spin-structure arc (Probes A–D)** — the genuine first-order operator | **`h = 0`, identical determinant modulus, and `η = ±1/4` as the ONLY distinction — so energetics, anomaly inflow and back-reaction ALL fail to select, and the arc's own moding swap is left conditional** | Replaces the 1D SUSY-square-root construction (which delivers `|spec(D)|` only) with the real objects. **A (spin lifts):** `J = −I₄` has `det = +1` (SO(4) — why RP³ is orientable); the only lifts are `J̃₊ = (−1,1)` and `J̃₋ = (1,−1)`, and **both square to +1**, differing by the covering kernel — so the two spin structures are distinguished by *which SU(2) factor carries the −1*, not by `J̃² = ±1` (refuting the shortcut the earlier `η = ±1` bookkeeping used). BAM's `σ` is **orientation-preserving** with `σ² = J`, and `T = iσ_y` gives `T² = J̃₊`, so **B2 selects a specific lift**. But B2 lives in `H¹(RP³;Z₂)` and `W` in `H¹(Γ;Z₂)` for the ER graph — **independent data**, so the rule `θ = 0 iff W·η = +1` is **unjustified** unless the network cycle generates `π₁(RP³)` (unproven). **B (spectrum):** from the analytic S³ eigenspinors `±(3/2+k)`, mult `(k+1)(k+2)`, with the two branches carrying **opposite** antipodal parity at the same `|λ|` — the origin of the asymmetry, and exactly what `|spec|` cannot see. Projection keeps exactly half the modes at every level (`vol(RP³) = vol(S³)/2` ✓). Results: **`η(0) = ±1/4` exact** (cross-checked against the convergent mode sum to 2e-19), **`h = 0` for both** (`|λ| ≥ 3/2` — the earlier 1D 'zero mode' was an artefact, confirmed independently), `ζ'_{D²}(0) = 0.2189594807`, and **identical `|det D|`** — the sectors differ *only* by a determinant phase `π/4` apart. **This retires the energetic route exactly:** the `π/(4C)` difference is identically zero in 3D. **C (inflow):** `Y³ = RP³ = L(2,1)`, `X⁴` = the D²-bundle over S² of even Euler number; APS boundary defect `∓1/8`; `Ω₃^Spin = 0` so **both extend** — inflow does not select. But `Ω₂^{Pin⁻} = Z₈` with `ABK(RP²) = ±1` **forces the two mouths of a throat to carry opposite Pin⁻ structures**, *deriving* the assertion `ThroatDefect` had been making — though only as a *relative* rule. **D (back-reaction):** precondition met; `ΔE_vac = 0` exactly and `Δ⟨T_AB⟩ = 0` **pointwise** by homogeneity, so identical radial system, regularity and stability — only the action *phase* differs, and a phase does not source `⟨T_AB⟩`. **Net: three independent selection mechanisms tested, none selects.** The lone survivor is `η = ±1/4`, which selects nothing energetically. Honest limits: round S³ exterior (not the Tangherlini throat), massless free fermion, and the bordism facts are cited not derived **E (does any of this apply?):** a loop in `RP³=S³/{±1}` is classified by lifting to `S³` — closes ⟹ class 0, antipodal ⟹ generator. Read from the constructions: #224's `Ltot = 2*_D_CH+4*_SIG_M+2*_LARC` with `_LARC=3.14` gives **two** π-arcs, total 2π, so the lift **closes**: **deck class 0, NOT the generator**. #223's single π-arc **is** the generator. And the ambient is not RP³ — a wormhole is a 1-handle, so `π₁ = π₁(exterior) * Zⁿ` and the #224 cycle lives in the **free handle part**; robust to reading the exterior as `S³` or `RP³`. **So Probes B–D, which target the deck Z₂, do NOT bear on the #224 network twist** — their numbers stand, their target narrows. Two consequences: the `W·η` composition rule is justified on the **#223** ring but not the **#224** ring (where the exchange-freeze headline lives), and the "three mechanisms all fail" result is about the **spin-structure** sectors — **the network twist sector remains untested by any of them** (`rp3_spin_lift_probe` / `rp3_dirac_eta_probe` / `bulk_aps_spin_structure_probe` / `twist_sector_einstein_dirac_probe` / `network_cycle_pi1_class_probe`, PR #230) |
| **The #223 ring — where the two Z₂ labels actually coincide** | **The moding shift is REAL there (measured), the composition rule is confirmed on an actual construction — and with the labels identified, B–D apply to the network twist and STILL nothing selects. Plus: the phenomenon and the coupling never co-occur** | Probe E found exactly one network with nontrivial deck class: the **#223** single-bridge ring (one π-arc), as against #224 (two π-arcs, class 0). So #223 is the one place the `W·η` rule is justified. **What `W` is there:** `build_ring` shoots from the neck to the source point at `χ = π`, and the `end` parity (Neumann `deriv` vs Dirichlet `val`) is exactly the two ways to close the deck-class-1 loop — so `end` **is** `W`. **The measurement:** flipping `W` must shift the comb by half a spacing, and it does — the two sectors **interleave** with fractional offset 0.529, 0.504, 0.499, 0.506, 0.498 → **0.5** (converging as modes climb out of the neck region), at both `r_s = 0.3` and `0.4`. *The first time in this investigation that `W` and the spin structure are shown to be the same label on a real construction.* **But nothing selects:** with the labels identified, Probe B–D results now bear on `W` itself — identical `\|det D\| = 0.8963003229`, `h = 0`, `η = ±1/4`, so `ΔE_vac = 0` exactly, both extend under inflow, back-reaction degenerate. The selection question is **answered for the one network where it is well posed: nothing selects the twist** — no longer a scoping caveat. **The structural by-product, sharper than the result:** #223 couples the labels but has **no mouth doublet** (the exterior is one connected cavity; the pure ultrastatic bridge has no interior state to exchange — #224's own words), while #224 has the doublet and the exchange freeze but deck class 0. **The phenomenon and the coupling do not co-occur on either of these two rings** — which is why the question kept slipping: every mechanism was computed on one ring and applied to the other. **Scope correction (PR #233):** that is a fact about *these two rings*, not about throat networks in general — a ring of two throats whose mouth pairs are a quarter circle apart carries the doublet, the exchange freeze **and** deck class 1 simultaneously, so a selection argument for the freeze sector **can** draw on RP³ spin-structure data. It still selects nothing there (`rp3_spin_structure_on_223_ring_probe`, PR #231) |
| **The freeze/deck "parity obstruction" — claimed, then RETRACTED by its own probe** | **The freeze and deck class 1 are INDEPENDENT, not mutually exclusive: two throats with quarter-circle arcs carry both — so the construction #231 recommended does exist** | Probe F observed that the exchange freeze (#224, 2 throats) and the spin-structure coupling (#223, 1 throat) never co-occur, and named a geometry carrying both as its next construction. PR #233 first answered that this was a **parity obstruction** — deck class = `N mod 2` versus freeze ⟺ `N` even — and withdrew #231's construction as impossible. **That answer was wrong, and this entry records the retraction.** **What stands:** the cyclic translation `S` advancing the ring by one throat obeys `S^N = W`; the freeze is *every* level forced into a degenerate multiplet, needing `S` to have **no real eigenvalue**, and for `W = −1` the eigenvalues are the `N`-th roots of `−1`, one of which is real iff `N` is odd — so **full freeze ⟺ `N` even**. Argued, measured on real rings (`N=2` twisted **2.2e-12 frozen**, `N=3` **4.0e-01 not frozen**, `N=4` **1.7e-12 frozen**), and this time verified at operator level: `‖[H,S]‖ < 1e-8` and `S^N = W` **exactly**. **What was wrong:** "an `N`-throat ring carries `N` exterior π-arcs, so deck class = `N mod 2`" silently assumed **every arc has length π** — and an arc of length π joins a point to its *antipode*, so assuming it for every arc forces all `N` throats to share one antipodal mouth pair. That is #224's configuration, not the general one. Probe E's actual criterion is **deck class = (Σ exterior arc lengths)/π mod 2**, which never mentions `N`: it is set by **where the mouths sit**. **The counterexample:** take `N` even for the freeze, then place the mouths so the arcs sum to an odd multiple of π — **two throats rotated a quarter circle apart**, each arc π/2, total π, **deck class 1**, two interior channels hence a **mouth doublet**, and twisted max intra-pair gap **8.482e-13 — FROZEN** (untwisted 1.9e-02, not frozen), by the same `S² = W` mechanism verified operator-wise. Likewise `N=2` arcs 3π/2 (3.3e-12) and `N=4` arcs 3π/4 (2.1e-12). **The other evasion classes, computed rather than asserted:** one throat with a finite interior channel (`D = 0.5, 4, 8, 16` — #231's own phrasing) gives a **nondegenerate interior ladder** (0.3785, 0.7568, 1.1350, 1.5131, 1.8910 at `D=8`), twisted and untwisted agreeing to 5 digits, **no doublet and no freeze** — the earlier probe's assertion here was right but was an assertion, and is moot anyway; **inhomogeneous or internally structured** rings (unequal arcs π/2π, unequal channels 4/5, an internal bump in one throat) **destroy** the freeze even at even `N` (gaps 2.5e-01, 3.8e-01, 4.3e-01) — it needs the cyclic symmetry *exactly*; a **non-cyclic theta network** of three throats between two tri-mouth junctions (`b₁ = 2`, four twist sectors) shows **no full freeze in any sector**, its degeneracies being `S₃` edge-permutation degeneracies present in the **untwisted** sector too. **Consequences:** Probe F's non-co-occurrence is **demoted** back to an observation about the two rings this repo contains; #231's construction is **restored** in corrected form (move the mouths, don't lengthen the channel); and the claim that the freeze sector is cut off from RP³ spin-structure data *permanently* is **retracted** — on the quarter-circle ring the cycle **is** the π₁(RP³) generator, so Probes B–D finally apply to a network that actually freezes. They still find nothing selecting the twist, but that null result is now a statement *about the freeze sector*. **Honest limit, running the other way from before:** the freeze is **fragile**, so any physical reading of it must explain why the exact cyclic symmetry would hold; and one theta network is a data point, not a classification of non-cyclic networks (`freeze_deck_parity_obstruction_probe`, PR #233) |
| **One throat, two mouth states — and why the twist still cannot touch them** | **A single throat DOES give two mouth-localized states at deck class 1; what it cannot give is a translation for the holonomy to act on** | #233 answered the one-throat case by looking at the interior channel *ladder* and finding it nondegenerate. But #224's doublet is a **mouth** doublet and a single throat already has **two mouths**. Put a barrier in the middle of the interior channel and each sub-well abuts one mouth: the pair rotates to states carrying **0.978 → 0.999** of their weight in a single well, with the splitting tunable **1.281e-01 → 5.849e-03** by the bump height, and one exterior π-arc keeps **deck class 1** throughout. *So #233's "no doublet on one throat" was too narrow — it generalized from the ladder.* **But the twist has no purchase on them.** The two wells are joined by **two paths** — through the central bump, or the long way round through the exterior arc — and only the second carries the holonomy, so `split = 2|t_bump + W·t_arc|` and the two twist sectors separate them. Extracted: **`|t_arc| = 2.0…2.9e-6`, independent of the bump height** (the two-path decomposition *confirmed*, not assumed) while `|t_bump|` falls **47×** over the same scan. The exterior path crosses **two mouth barriers** and is 3–5 orders of magnitude weaker, so the interference never balances: ratio twisted/untwisted stays **> 0.999** throughout. The resonant alternative is no better — tuning the channel so a channel level crosses an arc level gives a genuine avoided crossing (`D ≈ 8.25`, gap 3.2e-3, basins swapping across it) that the twist moves by only ~10%. **The mechanism, which is the real result:** the freeze is the `S^k = W` mechanism — a cyclic **translation** whose `k`-th power is the holonomy, `k` even so no eigenvalue is real. A one-throat ring has no such translation: its two mouth walls **both face the same interior channel**, so they are **mirror images**, and the ring carries a **reflection** `R` — with **`R² = +1` in both twist sectors, always**. A reflection has no holonomy-dependent square, so it can never do what `S² = W` does. Shown by building the same ring both ways, identical in lengths, basin profiles, deck class 1 and two mouths, differing **only** in wall orientation: **mirror** (physical) `‖[H,S]‖ = 4.90e+01`, `‖[H,R]‖ = 0.00e+00`, **neither sector freezes** (gap 5.299e-01); **translate** (not realizable by one throat) `‖[H,S]‖ = 0.00e+00`, `S² = W`, twisted sector **FROZEN at 2.954e-12**. The symmetry swaps exactly and the freeze follows the *translation*. **Consequences:** the governing integer is the **order of the cyclic translation**, which equals the throat count *because* each throat contributes one interior basin bounded by a mirrored pair of walls — so two throats are the minimum for a translation to exist at all, and #231's construction (#233's quarter-circle ring) stands for a structural rather than numerological reason. And **#234's gate set cannot move onto a one-throat network** — not because the doublet is missing (it is not) but because the Z₂ twist has nothing to act on, so the memory/gate switch has nothing to switch. The arc's through-line: the freeze was attributed first to **arc counting** (#233 v1, wrong), then to **channel structure** (#233 v2, too narrow), and now to whether a cyclic translation commutes with `H` — the first of the three that is checkable on *any* network (`single_throat_mouth_doublet_probe`, PR #235) |
| **What caps BAM's nonlocality at Tsirelson — and why the answer matters** | **The bridge saturates 2√2 and can never exceed it; but the ceiling is NOT geometric — it follows from the readout being two-valued, so the residue of "quantization is imported" is exactly one identity** | #206 escapes Bell's theorem by spending the bridge as its one nonlocal resource, reaching CHSH = 2√2. But nonlocality is not a dial that stops where quantum mechanics stops: **Popescu–Rohrlich boxes are perfectly no-signaling and reach the algebraic maximum 4**. Tsirelson's number appears in eight documents of this repo, always as the value being *matched*, never as something the geometry *explains*. So: can the bridge be driven past it, and if not, **which** structure forbids it? **The ceiling is real.** Maximizing over every bridge the model has — all 8 gluing holonomies (the π-holonomy handle is only one), all preparations (β, α) — the largest CHSH is **2.828427125** against Tsirelson 2.828427125: excess **+8.9e-16**. BAM does not predict super-quantum correlations; on this axis it agrees with experiment *for a reason*, and is falsifiable in the other direction. **The settings are exactly sufficient, neither deficient nor excessive.** BAM's local setting is a *fiber-frame rotation*, which on the `k = ±1` qubit is `diag(e^{iθ}, e^{−iθ})` — a **one-parameter U(1)**, not SU(2) — so conjugating the fixed transverse readout sweeps the x–y plane and nothing else. That is exactly enough, for a reason worth stating: **CHSH's optimum lies in a plane**. The lattice-derived state is the singlet at fidelity **1.000000000000**, its T-matrix x–y block has singular values **(1, 1)**, and the plane-restricted Horodecki value is 2√2 — identical at grid resolutions n = 25, 49, 97, 193, so a saturation and not a discretization artifact. **But locality is NOT what caps it — this is the result.** The obvious guess is that the mouths' operations act on commuting algebras, and they do overwhelmingly: `‖[P_A,P_B]‖ = 4.7e-113`, because the committed bridge is a **link in `H` between distinct lattice sites, not an identification of degrees of freedom** (worth stating plainly, since #206's prose — "one object", "two frames of one bulk defect" — reads like an identification and the implementation is not one). Yet **dropping commutativity does not help**: for *any* dichotomic Hermitian observables on **one** Hilbert space, with no tensor split at all, `(B₀+B₁)² + (B₀−B₁)² = 4I` identically (residual **5.3e-15** over 800 random draws), and Cauchy–Schwarz then caps CHSH at 2√2 **without invoking locality anywhere** — measured, taking the optimal `A` for each `B`, the best attainable value is **2.828427095**, Tsirelson again. **So the ceiling is supplied by the readout being a ±1-valued observable on a Hilbert space, not by the bridge, the mouths' separation, or anything else geometric.** **What that does to "quantization is imported" (#232):** it prices it. Three of the four ingredients of a Bell test *are* genuinely geometric — the state (field + bulk gluing), the nonlocality (raw B-side amplitude **3.747e-01** with the handle vs **4.140e-16** without, with local strategies capped at 2), and the settings. The fourth, the ceiling, is not, and the **entire** missing ingredient is one algebraic identity: `B² = I`. Not Hilbert space in general, not the Born rule, not linearity — just the readout being two-valued. *Naming it is the progress; deriving it is not attempted and is not claimed* — an argument that a classical throat readout must be two-valued would close the gap, and the Pin⁻/`T² = −1` machinery (#195/#196, #227–#231) is where it would have to come from, but none of it currently addresses the readout algebra. **A second finding, on the nonlocality budget:** #206's T6 checks that Bob's reduced state is invariant under Alice-side unitaries — a **partial-trace identity true of every bipartite state**, so it tests nothing about the geometry. Tested on the lattice instead (Alice's setting applied as an operator supported at mouth A, then evolution, then Bob's channel populations), no-signaling is an **equal-time** statement: exactly 0 at `dt = 0` (disjoint supports), but the committed handle is *always on*, so shifts of **2.273e-02 / 1.019e-01 / 2.315e-01** appear at `dt = 0.5 / 2 / 8`, scaling with the coupling (3.6e-02 at gb = 0.2) and vanishing to ~1e-15 when cut. This substantiates with numbers the caveat #206 states only in prose ("the physical bridge is non-traversable; the lattice handle is a modeling stand-in") and converts it into a **requirement**: the bridge must be dynamically inert at measurement time, which the committed model is not — a constraint on the model, not a refutation of the program, but one to discharge rather than assert. **Honest scope:** one 2-outcome, 2-setting scenario on one lattice; nothing here bounds multipartite or higher-input scenarios, where the quantum set is not characterized by a single inequality (`tsirelson_ceiling_probe`, PR #236) |
| **Which correlation tables BAM's encoding actually admits — and why the first answer was wrong** | **The representation admits exactly the Gram/TLM body; BAM's *single-shot* family is only ~55% of it, but the **convex hull** of its `c=1` tables is the **whole body**, so the only surviving restriction is identically zero marginals** | #236 located the Tsirelson ceiling at `B² = I` by probing the CHSH *maximum*. A maximum is one number; the encoding represents a whole table. **Level A — the representation in the abstract** admits exactly the **unit-vector Gram matrices** `E_xy = ⟨u_x, v_y⟩`, equivalently the Tsirelson–Landau–Masanes body `max over sign variants of |Σ ± arcsin E_xy| ≤ π`. Verified in *both* directions rather than cited: an independent fit agrees with the TLM predicate on **60/60** random tables, and 4000 tables built *from* random unit vectors never violate TLM (min slack **+2.6e-04**). The three bodies by Monte Carlo: local **66.83%** of the correlation cube, quantum **92.64%**, no-signaling **100%**; the PR box at slack **−π**. **Level B — what BAM implements single-shot is strictly smaller.** Read off the committed lattice rather than assumed: the pair state's T-matrix x–y block is **exactly `c·I₂`** (deviation ≤ **3.1e-15**) with **`c = −sin 2β`** to 1e-6, so with #236's U(1) settings the achievable tables are **exactly `E_xy = c·cos(θ_x − φ_y)`**. At maximal preparation that is a 3-parameter surface in a 4-dimensional body — **measure zero**, `d=2` Gram fit succeeding on **0.0%** of random quantum tables against **100%** at `d=3` — rising to **55.2% ± 2.0%** with `β` free, so `β` is the knob buying the missing dimension. Of the `c=1` surface, **33.7%** lies exactly *on* the Tsirelson boundary, which is why #236 saw saturation there: **saturating the boundary and filling the body are different things.** **But the single-shot gap is not an operational restriction — and this corrects the PR's own first framing.** Shared randomness over preparations is an ordinary resource, and it makes the reachable set the **convex hull** of the `c=1` cosine tables. Tested by **linear programming** (`min ‖Gᵀλ − E‖₁` over the simplex): coverage **99.2 / 100.0 / 100.0%** at 2k / 8k / 32k generators, **98.7%** on **300** fresh random TLM tables, with support size at most **5** — exactly the **Carathéodory bound** for a 4-dimensional body. The residual failures are boundary resolution rather than a gap: refining to 128k generators *on those tables only* takes 4 failures → **1** and the worst residual **7.2e-03 → 4.2e-03**. The PR box stays outside (**1.18**), as it must. **And the committed 'unreachable' witness decomposes into exactly 5 `c=1` cosine tables** — residual **0.0e+00**, reconstruction error **4.4e-16**, weights `[0.2124, 0.3447, 0.0409, 0.0102, 0.3918]`. So *"~45% of quantum correlation tables are unreachable" was a **single-shot** statement dressed as an operational one*; with mixing there is no correlation-table restriction to report. **The one restriction that survives:** under BAM's own U(1) settings `max |⟨A_θ⟩| = 0.000e+00`, so **no behavior with biased marginals is representable at all** — not for want of structure in the state (it carries a z-marginal `cos 2β`, 0.866 at `β = π/12`) but because the settings never leave the x–y plane. **Unlike the correlation gap this one is immune to mixing** — a mixture of zero-marginal behaviours has zero marginals — so it is the whole of the probe's falsifiable content, and the restriction sits in the **readout**, the same place #236 found the ceiling. **A methodological trap, recorded:** single-shot membership was first tested with generic optimizers, which got it wrong **in both directions** (L-BFGS-B 48.3%, Nelder-Mead 63.3%) — *a failed fit is not proof of non-representability*. Replaced by an **exact** test: fixing the global rotation with `θ₀ = 0` leaves four unknowns for four equations, three solving in closed form, leaving `E₁₁` as one consistency condition in `c` alone — a 1-D root find on 8 sign branches, bracketed and bisected, validated by recovering **100%** of family-built tables (max residual 3.8e-14). Scopes #236's "exactly sufficient" to the CHSH maximum. Open: whether the throat geometry supplies any mouth operation rotating out of the x–y plane; whether BAM has access to shared randomness at all (assumed, not derived); the multipartite and higher-input cases (`admissible_correlation_tables_probe`, PR #237) |
| **The detector algebra and the marginals — and why the fully derived pairing does not violate Bell** | **The measurement algebra is full; what is restricted is which effects are *dialable* — and with the setting the docs derive and the detector #209 derives, CHSH = 2.000000 exactly** | #237 established that correlator coverage is **not** the remaining problem: with convex mixing the hull of BAM's tables is the whole quantum body. What remains is the **detector marginals** and the **measurement algebra** — and this probe finds they are *one* question. **The generated algebra, computed separately for each pairing.** For the **fully derived** pairing — the fiber U(1) setting with the `σ_z` winding Stern–Gerlach — the accessible observable set is a **single point** (a fiber rotation cannot move a winding-diagonal detector) and the `*`-algebra it generates is exactly **`span_C{I, σ_z}`**: complex dimension **2**, and **abelian**. That is the real explanation of the `CHSH = 2.000000` below, and it is a **theorem rather than an observation** — a commutative algebra of observables admits a joint distribution over all of them at once, hence a local hidden-variable model, hence `CHSH ≤ 2` necessarily. Both violating pairings generate all of `M₂(C)` (complex dim 4, non-abelian), which is what makes a violation possible at all. *(This corrects the probe's own first answer, which computed the algebra for one pairing only and concluded "the algebra is not the restriction".)* **The docs and the committed code disagree about what a setting is.** The documents describe *"the fiber-frame rotation before the device"*; a rotation of the `χ` fiber multiplies the winding-`k` channel by `e^{2πikδ/N}`, so on the `k = ±1` qubit it is `diag(e^{iθ}, e^{−iθ})` — a **z-rotation**. But committed `measurement_sector_probe._rot` is `[[c,−s],[s,c]] = exp(−iθσ_y/2)` — a **y-rotation**, imported **directly** and verified identical to **0.0e+00**. (Method check: the code rotates the *state*, this probe conjugates the *observable*; identical correlators to **3.3e-16**.) The difference is exactly the decisive one: the fiber rotation **commutes exactly** with the `σ_z` winding Stern–Gerlach (`‖[fiber, σ_z]‖ = 0.00e+00` at every angle) while the y-rotation does not (up to 1.95). **A fiber rotation cannot move a winding-diagonal detector.** **THE TRILEMMA, measured on the committed #206 singlet:** fiber U(1) × `σ_z` winding SG — *both derived* — gives **2.000000** with dialable span **1**; fiber U(1) × `σ_x` transverse *(assumed by #236/#237)* gives **2.828427**; y-rotation *(the code)* × `σ_z` gives **2.828427**. *The only pairing in which both the setting and the detector are physically derived by this repository gives exactly the local bound — no Bell violation at all*, because the setting commutes with the detector and the dial does nothing. The violation requires **exactly one** non-derived ingredient: a transverse detector, or a winding-mixing setting. **The marginals follow the same fork, and this corrects #237.** #237 reported *"identically zero marginals, immune to convex mixing"* as its one surviving falsifiable content; that was measured with the `fiber × σ_x` pairing and is correct **there**, but under the other two the marginals are **|cos 2β|** exactly — 0.8660 at `β = π/12`, matching to 1e-6. So the zero-marginal restriction is an artifact of one of three modelling choices and is **retracted as a general claim about BAM**; #237's falsifiable content was contingent on the same undetermined choice the Bell violation itself hangs on. **Both missing ingredients are the same physical requirement.** `σ_x` in the winding basis is purely off-diagonal between `k = +1` and `k = −1` — a coherence between distinct **topological charges** — and the y-rotation setting mixes those same sectors. So the one thing BAM must supply is *an operation or observable at a mouth that coherently connects distinct winding sectors*. And if winding is superselected that is unavailable: restricting both parties to `k`-diagonal dichotomic observables gives max CHSH **2.000000** over 4000 random draws, exactly the local bound. The repo already carries a superselection structure here — **#208's charged-GHZ zero-sum no-go** — and nothing in it has been reconciled with the Bell chain. **What this does and does not say:** it does *not* show BAM's Bell violation is wrong. It shows the violation rests on one ingredient the program has not derived, that the two candidate ingredients are the *same* requirement, and that the requirement is in tension with a superselection structure committed elsewhere. Naming which single thing must be supplied is the deliverable. Open: reconciling #208 with the Bell chain; whether relocating onto the transported-frame / spin doublet (#192/#197 Berger-squash, named in #209 as the device for *that* carrier) avoids the problem — the natural next probe (`detector_algebra_probe`, PR #238) |
| **The minimal missing interaction, tested directly — it exists, and both costs first reported were overstated** | **#238's gap is fillable by exactly one term (the second χ-harmonic at a mouth), which restores CHSH to 2.828427 with the *derived* detector — and review corrected both costs: total charge IS conserved with a winding-2 carrier, and the operational Bell value is 2.3309 broadly, not a narrow 2.13** | #238 showed the Bell chain needs one undelivered thing: an operation coherently connecting distinct winding sectors. The obvious move was to relocate the chain to the spin doublet — but **assuming the interaction cannot be had would repeat the error #233 made about the one-throat ring**, so it was tested directly first. **It is unique at first order.** Among mouth-localized fiber-angle potentials `λ cos(2πmχ/N_χ + φ)`, only **m = 2** connects the qubit: `|⟨+1|V|−1⟩| = 0.4082`, against **exactly 0** for m = 1, 3 and 4. The reason is elementary — the qubit is `k = ±1`, so the coupling must carry `Δk = 2` and the m-th harmonic carries `Δk = m`. So the minimal missing interaction is the **second χ-harmonic at a mouth**: physically a **fiber-angle anisotropy**, a mouth that is not round in the fiber direction. **It is exactly the missing generator.** Restricted to the qubit it is **pure `σ_x`** with coefficient **0.408248** — the `I`, `σ_y` and `σ_z` coefficients are all zero to **1e-16** — and it fails to commute with the derived detector (`‖[V, σ_z]‖ = 0.4082`), so conjugating `σ_z` by its settings sweeps the x–z circle instead of standing still. **It works, with the derived detector and no transverse readout assumed anywhere:** max CHSH goes **2.000000 → 2.828427**, i.e. Tsirelson. *So #238's gap is not a no-go, and the question becomes what filling it costs.* **CHARGE — RETRACTED.** The first version reported `‖[H, K]‖ = 6.6e-16` for the committed dynamics against `‖[V₂, K]‖ = **1.414**` and concluded the interaction is *"a charge-non-conserving mouth term"*. **That overstates what the number shows.** A *prescribed external potential* breaks the winding of the **throat subsystem treated as closed**; it does not show *total* charge must be abandoned, any more than an external field means angular momentum is not conserved — the compensating charge sits in whatever sources the potential. Tested directly: extend the model with a **winding-2 carrier** and `K_total = k_field + 2·n_carrier` is conserved **exactly** — `‖[H_ext, K_total]‖ = **0.0e+00**` at every truncation — with the large-amplitude (mean-field) limit reproducing the prescribed second harmonic at the **same coefficient 0.408248**. *So the prescribed term is the mean-field limit of a charge-conserving interaction*; the real cost is a **carrier requirement**, not a broken conservation law. **THE BELL WINDOW — RETRACTED AND REPLACED.** The first version combined an **ideal projected-qubit CHSH** (with no leakage in it at all) with a **separate leakage proxy**, fed both into `S > 4/η − 2`, and reported *"a narrow loophole-free window peaking at S = 2.13"*. **Those are two different experiments and the combination is not a Bell test.** Computed properly — each mouth having **three outcomes** (`k=+1`, `k=−1`, or **no click** when the excitation has left the encoded sector), the full winding space evolved so leakage sits in the dynamics rather than bolted on, no-click assigned to a fixed outcome with **both** assignments tested and agreeing exactly, and probabilities verified to normalize to **1.000000000000** — the operational CHSH is **2.032289** at `|t| = 0.2`, **2.117169** at 0.4, **2.306045** at 0.8, peaking at **2.330905** at `|t| = 1.0` (`P(both click) = 0.709620`) and still **2.326141** at 1.9. *The detection-loophole-free violation is real and **broad** — above the local bound at every span tested — and **stronger** than first reported: 2.3309, not a marginal 2.13.* **What survives is only a ceiling:** leakage caps the operational violation at **2.3309** against Tsirelson **2.8284**, conceding roughly half the available margin — the one cost of the two that review left standing. **Net.** Testing beat assuming — the interaction is real and does the whole job. But **neither cost is what this probe first claimed**: charge conservation is *relocated into a carrier* rather than lost, and the Bell violation is *broad* rather than marginal. The case for relocating to the spin doublet is therefore **weakened** — it now rests on a carrier requirement and a 2.33 ceiling, not on a broken conservation law. Open: whether a **winding-2 carrier** at a mouth is realizable in any BAM geometry (that, not charge conservation, is what the construction requires); whether a `k = ±3`-resolving detector raises the 2.33 ceiling; whether the spin doublet supplies a rotation generator without breaking any conservation law — the relocation this probe was run to inform, still untested (`minimal_mixing_interaction_probe`, PR #239) |
| **The charge-conserving apparatus and the complete pointer statistics — the winding carrier stands** | **#239's open item executed: the apparatus conserves the Z₈ charge exactly, every 'leaked' channel is actually detected, and the complete statistics reach Tsirelson — so the 2.33 ceiling was an artifact and no relocation is needed** | #239 closed by asking for a **charge-conserving throat–apparatus interaction** and the **complete multi-outcome pointer statistics** before deciding whether the winding carrier must give way to the spin doublet. **The apparatus.** The *setting* interaction moves the field's winding by `Δk = +2` while absorbing one quantum from a **carrier** holding 2 units of winding each; the *measurement* interaction is the committed winding Stern–Gerlach `g·K_field ⊗ p̂`, winding-**diagonal** and therefore charge-conserving identically. **A correction to #239 falls out:** winding on a discrete `N_χ = 8` fiber is a **Z₈** charge, not an integer one, so the conserved operator is `exp(2πiK_total/8)` — `‖[H, Z₈]‖ = **1.1e-15**` across truncations, while the integer-`K` test gives **9.24, 13.06, 18.48**. #239's T5 obtained `0.0e+00` only because it *excluded the wrap-around transitions*; with them included the integer test fails. The conclusion (total charge is conserved) survives — the operator expressing it does not. **Every channel is detected.** Under `Δk = ±2` the orbit of `k = +1` closes on the **4-cycle {+1, +3, −3, −1}**, so the population #239 called *leakage* sits in `k = ±3` and nowhere else — and a winding SG deflects **∝ k**, so those land at their own pointer positions rather than failing to register (at `g_p = 2`: centres 2, 6, −6, −2, minimum separation **4.0σ**, branch overlap **4.6e-02**). *The apparatus has four outcomes per mouth and no no-click outcome at all; #239's three-outcome model was **wrong about the detector**, not conservative about it.* **The ceiling disappears.** Binning the complete statistics with the natural setting-independent rule `sign(k)`: **2.134677** at `|t| = 0.4`, **2.457397** at 0.8, **2.758947** at 1.2, **2.824508** at 1.5, **2.828028** at 1.9 — Tsirelson (2.828427) to grid accuracy, and itself the best of **all 16** setting-independent binnings. #239's rule (send `k = ±3` to a fixed outcome as though undetected) is one of those 16 and yields only **2.33** on the same data. *So the '2.33 leakage ceiling' was not a property of the winding carrier; it was the cost of throwing away channels the apparatus resolves.* **Finite-carrier back-action, exact** (carrier kept and traced out, no mean-field step): monotone **2.290317 / 2.601231 / 2.760053** at `n̄ = 1 / 4 / 16`, converging to the mean-field **2.8245**, with coherent-state truncation ≤ 2.5e-13. *A bug worth recording: an earlier run showed a spurious **non-monotonic drop** at large `n̄`, caused by a coherent state truncated below its own mean.* **No-signaling** on the complete four-outcome statistics holds to **2.2e-16**. **THE DECISION: the winding carrier need not be replaced.** Across #239 and this probe, **three separate costs were reported and all three dissolved on closer calculation** — a prescribed potential mistaken for broken charge conservation, an efficiency proxy mistaken for a Bell test, and detected channels mistaken for lost ones. What survives is a single requirement *on the geometry*: a **winding-2 carrier at the mouth, populated to `n̄ ≈ 16`**. The spin doublet may still be worth exploring on its own merits, but **not as a rescue** (`throat_apparatus_pointer_probe`, PR #240) |
| **Pin⁻ on the throat's RP² mouth** — the exchange sign and the Fermi equation of state | **The deferred calculation, done: RP² admits Pin⁻ only; its spin-½ spinor (2π=−1) + Finkelstein–Rubinstein gives the −1 exchange sign; antisymmetry → Pauli → the Fermi EoS (P=⅔u Γ=5/3 NR; P=⅓u Γ=4/3 UR; T=0 degeneracy pressure)** | Takes the Pin structure on the non-orientable throat mouth (#169) and shows it *delivers* the physics, not just the topology. **The Pin structure:** the Stiefel–Whitney classes `w(RP²)=(1+a)³=1+a+a²` give `w₁=a` (non-orientable, no Spin) and `w₂=a²`, so RP² admits **Pin⁻ only** (`w₂+w₁²=0`), not Pin⁺ (`w₂≠0`) — a unique, definite spinor structure on the mouth. **The exchange sign:** the Pin⁻ spinor is spin-½ — `R(2π)=exp(−iπσ_z)=−I`, `R(4π)=+I` (explicit) — and by Finkelstein–Rubinstein the two-throat exchange is homotopic to a 2π rotation of one (the orientation-entanglement / belt trick; two-particle config-space `π₁=ℤ₂`), so the exchange sign is **−1** and the wavefunction antisymmetric; the spin-statistics connection is *realised* by the same holonomy, not assumed. **The Fermi EoS:** antisymmetry → Pauli exclusion (`n_p∈{0,1}` vs Bose `{0,1,2,…}`) → filling the Fermi sphere gives the degenerate equation of state `P=⅔u, Γ=5/3` (non-relativistic) and `P=⅓u, Γ=4/3` (ultra-relativistic), with a **strictly positive T=0 degeneracy pressure** (the support of white dwarfs / neutron stars) that a Bose gas (collapsing to `p=0`, `P=0`) lacks. **Honest scope:** computed — the Pin⁻ classification, the spinor 2π sign, and the Fermi EoS integrals/indices; cited (not re-derived) — the Finkelstein–Rubinstein exchange↔rotation homotopy, the one configuration-space theorem linking the throat's internal Pin holonomy to the physical exchange. The same orientability grading as #63 (`C=iσ_y`, `T²=−1`) and #67 (even-`k` absence), carried through to statistics and the equation of state (`pin_rp2_fermi_statistics_probe`, PR #170) |
| **Tangherlini J-quotient consistency** — the topological root of the non-orientable throat | **One free isometric antipodal involution: bulk S³/J = RP³ orientable (det +1), brane mouth S²/J = RP² non-orientable (det −1); the #167 non-orientable throat is the RP² cross-cap forced by the bulk→mouth dimension drop** | Explains WHY the throat is non-orientable (#167) while the bulk is not, as a dimension-parity statement about the antipodal (J) quotient. **The split:** the antipodal involution `J: x↦−x` on `Sⁿ` has orientation determinant `(−1)^{n+1}` — orientation-preserving for ODD `n` (orientable `RPⁿ`), reversing for EVEN `n` (non-orientable). The bulk angular sphere is `S³` (odd → `RP³` orientable, det +1, computed explicitly via the pushed tangent frame); the throat mouth is the brane's angular `S²` (even → `RP²` non-orientable, det −1). **The consistency:** the *same* `J` acts oppositely because the two spheres sit one dimension apart, on opposite sides of the parity; `J` is FREE (`−x=x⟹x=0∉Sⁿ` → smooth manifold quotient) and an ISOMETRY (`JᵀJ=I` → the round-angular Tangherlini metric descends). So the non-orientable throat is *forced* by the single-dimension drop from bulk to mouth, not assumed. **The #168 realization:** in those coordinates `J=(χ,θ,φ)↦(π−χ,π−θ,φ+π)` fixes the equatorial `χ=π/2` brane and restricts to the `S²` antipodal map `(θ,φ)↦(π−θ,φ+π)`; the metric descends, so the #167 non-orientable throat is exactly the `RP²` cross-cap inside the orientable `RP³` bulk of #168. **Consistency (remark, not a new derivation):** `RP³≅SO(3)` is orientable/spin; `RP²` is non-orientable and admits only a Pin structure — the half-twist carrying the spin-½ character, the same orientability grading as the C-swap (`C=iσ_y`, `T²=−1`; #63) and the even-`k` absence (#67) (`tangherlini_j_quotient_probe`, PR #169) |
| **The global regular 5D embedding** — the BAM throat as the equatorial Tangherlini slice | **The 5D derivation #167 flagged, supplied: an explicit GLOBAL REGULAR exact bulk — the BAM throat is the totally-geodesic equatorial slice of the 5D Tangherlini vacuum; three checks + the regularity gate pass; the gap closes** | Closes PR #167's gap (necessary conditions → sufficient) by building the **global regular** embedding the program flagged — not Campbell–Magaard local existence. **The construction:** the BAM brane is the equatorial `χ=π/2` slice of the 5D Schwarzschild–Tangherlini bulk `ds²₅=−F dt²+dρ²/F+ρ²dΩ₃²`, `F=1−μ/ρ²`, with `μ=r_s²`. The equator is a Z₂ fixed-point set, hence **totally geodesic** (`K_μν=0`): a tension-free, matter-free brane; the construction works only for the pure-tidal `M=0` form (a Schwarzschild `1/r` term has no 5D-Tangherlini counterpart), so the gate has teeth and ties the bulk mass to the throat scale. **The three printed checks:** (1) the induced 4D metric is exactly `f=1−(r_s/r)²` with `K_μν=0`; (2) the projected bulk Weyl `E_μν=−G⁴_μν` to ~1e-8 — the brane's effective exotic stress (`ρ_eff<0`, the tidal fluid) **is** the projected Weyl of the ordinary 5D vacuum (the bulk-Weyl mechanism, made explicit by an actual solution); (3) the bulk is Ricci-flat to ~3e-7 — an ordinary 5D vacuum, no 5D matter or exotic source. **The regularity gate:** the coordinate-invariant Kretschmann `K₅=72μ²/ρ⁸` (closed form validated numerically to ~1e-6 at `ρ≥1.5r_s`, away from the `1/F` coordinate breakdown) is **finite throughout** the exterior `ρ≥r_s`, max `72/r_s⁴` at the throat — the only singularity `ρ=0` is behind the regular 5D Killing horizon, and the extra dimension `χ∈[0,π]` is compact and regular. **GLOBAL and REGULAR — the gate passes.** **What closes:** #167's bulk-Weyl reading is now *realised*, not consistent-with — no exotic brane matter, no brane gauge field, and the `f=0` throat is identified as the **regular 5D Killing horizon** (improving #167's caveat: regular, not singular). Honest residue: the throat sits at a (regular) horizon; the brane is the tension-free totally-geodesic slice (`μ=r_s²` fixed); it is the exterior embedding `ρ≥r_s` (`global_regular_5d_embedding_probe`, PR #168) |
| **Israel junction audit for the non-orientable throat** — the braneworld Weyl split | **Throat stress is the tidal-charge / bulk-Weyl form; on-brane exotic matter avoidable *if* the 5D embedding sources E_μν and BAM has no fundamental brane gauge field — necessary conditions met, 5D derivation pending** (closed positively by PR #168) | The naive Israel result is settled (thin-shell σ<0, WEC-violated; the non-orientable Z₂/C-swap gluing does NOT rescue the sign), so the real, *non-predetermined* question is the **braneworld split**: is the throat's negative effective-4D σ supplied by the projected bulk Weyl term `E_μν` of an ordinary 5D bulk, or is exotic brane matter irreducible? **The eight deliverables:** the Lanczos `S^a_b=diag(−σ,p_t,p_t)`; `σ=−√f(a)/(2πa)<0`; `p_t=(1/4π)[f'/(2√f)+√f/a]`; `S_ab kᵃkᵇ∝σ+p_t` (null); WEC violated on the shell (exotic); σ has the wormhole sign (wrong for ordinary matter) and the inverse-throat scale; and the discrete Lanczos σ is recovered from a tanh wall as thickness→0. **The decisive fact (computed):** `f=1−(r_s/r)²` is **Ricci-flat** (R≤1e-17), and its effective 4D stress is **traceless** with the `r⁻⁴` form (`ρ_eff=−r_s²/(8πG r⁴)<0`, `p_r=−ρ`, `p_t=+ρ`; `r⁴ρ_eff=−r_s²`) — the tidal-charge / bulk-Weyl form. **The split (honest):** the stress is 100% of the bulk-Weyl FORM, but the no-brane-exotic ATTRIBUTION is *consistent-with, not proven*. Its **necessary** conditions are met — `R=0` (required for a vacuum brane `G_μν=−E_μν`, Shiromizu–Maeda–Sasaki) and a **negative** tidal charge (`ρ_eff<0`) that *excludes a real on-brane Maxwell source* (which gives ρ>0; the same r⁻⁴ form is Reissner–Nordström, only the sign distinguishes them), so the reading also requires BAM carry no fundamental brane gauge field. The **sufficient** step — the explicit 5D embedding whose Weyl projection sources exactly this `E_μν` — is **pending** (Dadhich/Bronnikov–Kim cited, not re-solved for BAM's bulk). **And the f=0 horizon is not evaded:** the surgical surface term vanishes there, but f=0 is a null/degenerate locus, so this relocates σ rather than removing it. This is the strongest "consistent-with" the audits have reached — narrow, specific, closable by a 5D embedding calculation (`israel_junction_weyl_split_probe`, PR #167) |
| **Antipodal wave-packet focusing threshold** — the focusing computed dynamically | **A conformal wave packet on S³ refocuses EXACTLY at the antipode at t=πR (machine precision); the 1/sin²χ caustic is the geometric trigger for the 2 m_e c² pair-nucleation threshold** | Closes the THESIS "antipodal focusing" gap: the reconvergence was asserted (and the 2 m_e c² threshold derived *statically* in PR #58), but the focusing itself was never simulated. **The reduction:** the zonal sector of S³ reduces EXACTLY to a 1D string (modes `sin((ℓ+1)χ)`, Dirichlet ends), with the physical field `ψ=f/sinχ` carrying the geometric focusing factor `1/sinχ`. **Exact refocus:** a conformal packet (`ω_ℓ=(ℓ+1)/R`) launched near χ₀ refocuses at the antipode `π−χ₀` at `t=πR` (half the great-circle period) — the identity `ψ(χ,πR)=−ψ(π−χ,0)` holds to **3e-15**, amplitude recovery ×1.0000; at `t=2πR` it **fully revives** (`0e+00`) — the sub-threshold focus passes through and re-disperses. **Conformal required:** the sharp focus needs the equally-spaced conformal tower (recovery ×1.000) — the minimally-coupled tower `√(ℓ(ℓ+2))` dephases (×0.877); the same coupling that makes the S³ vacuum tower equally spaced (PR #165) makes the caustic sharp. **The caustic:** density `∝1/sin²χ` diverges as the wavefront converges, regularized by `ℓ_max~R/R_MID` — it lets a delocalized, S³-wide wave reconcentrate onto the throat scale (concentration ~`R/R_MID`), the dynamical bridge from a diffuse wave to a local nucleation density. **The threshold (inherited, PR #58):** focused energy ≥ `E(R*)=m_e c²`; the C-conjugate pair (`Σc₁=0`) → `2 m_e c²=1.022 MeV`; disperse-below / persist-above. **Honest scope:** linear conformal focusing is computed exactly; the *nonlinear* throat formation is named, not simulated (zero fit constants; the focus is the trigger, the particle the persistent response) (`antipodal_focusing_threshold_probe`, PR #166) |
| **Berger deformation audit of R-unification** — *clean failure (negative result)* | **ρ(1) lands ~35 orders off the measured cosmic/particle ratio: the global cosmic-cavity Casimir and the local throat self-energy do NOT ride on one R** | An AUDIT (not a quantum / throat-formation / wave-propagation test): does BAM's unified mass operator `m²=(k·2π/L_throat)²+((n+1)π/L_cavity)²` really ride the throat (Hopf-fiber winding, `L_throat=√(2π)/k₅`) and the cavity (radial/base) on ONE S³ radius? The Berger sphere `S³_λ` squashes the fiber alone — the one move that separates the two scales. **Guardrails (anti-rigging):** no derived inversions (no fitted constant relabelled as a π-multiple; enforced by source scan), no hidden Born/singlet imports (`_forbid_quantum_inputs` raises instead), no false victories from the `A/R+B·R²` well's stability (computed, then DISCOUNTED). **Global Casimir** `E_cav(λ)`: zeta-regularized conformal scalar on the genuine SU(2) Berger spectrum `Δ=4j(j+1)+4m²(λ⁻²−1)`, **validated at λ=1 against the exact closed form `1/240R`** (anomaly-free; a residual `1/n` log-ambiguity grows with the squash). **Local self-energy** `λ_min(λ)=√((2π/(λL_throat))²+ω₀²)` **moves** ×1.99 across λ∈[0.7,1.4]. **The clean failure:** R-unification forces `ρ=E_cav/E_self` to be a pure number — computed `ρ(1)=3.3e-4` — but the measured global/local ratio `λ_C/R_Hubble=3.0e-39` is **~35 orders of magnitude off** (the cosmological-constant problem, geometrically). **Survives only as** scale-free bookkeeping: `ρ(λ)` is parameter-free but NOT flat (cavity and throat respond differently to the same deformation), so they are not one dynamical object even in shape; consistent with the B4 single-anchor audit. All three guardrails held (`berger_r_unification_audit_probe`, PR #165) |
| **The v4 library migration** — the flavor-CP lock lands in `geometrodynamics/qcd` | **Migrated additively over the FROZEN v3 lock; v3 bit-reproducible; v4 inherits the masses and realizes the nine observables ≤ 1% from a library call** | The #163 successor: the v4 candidate lock moves from probe-local code into the calibrated library, in the four staged steps. **The surface:** six new `QuarkParams` fields (the Hopf phase `phi_h`; the three targeted couplings `eta_k1k3_plus / eta_k1k3_minus / eta_k1k5_minus`, subtract convention; the per-shell `diag_shift_plus / diag_shift_minus`), the `extract_ckm_matrix()` reader, and `LOCKED_QUARK_PARAMS_V4` at the derived `φ_h = π/k₅` — all DEFAULT-OFF. **Two views (the #158 relocation, in code):** the holonomy is a pure phase, so `extract_physical_spectrum` STRIPS `φ_h` (the v4 lock inherits the v3 masses to ~3e-9) while `extract_ckm_matrix` KEEPS it (the CKM with the physical Jarlskog). **Default-off:** at `φ_h = 0` the v3 Hamiltonian is exactly real and its CKM is a real rotation with `J = 0` — every PR #155–#162 probe pins to the frozen v3 lock and is bit-for-bit untouched; the migration is ADDITIVE, not a re-baseline. **Verified from the library:** `\|V_us\| ×1.00, \|V_cb\| ×1.00, \|V_ub\| ×1.00, \|V_td\| ×1.01, \|V_ts\| ×1.00, J ×1.00`, `(β, γ, α) = (22.3, 65.9, 91.8)°` vs `(22.2, 65.9, 91.9)°`, `sin δ = 0.889` vs 0.887, unitary to 7e-16. **Counting unchanged** (#163): +3 parameters for +5 independent observables (net +2); CP at ZERO parameters; the #150 budget unchanged. **The unmixed reference:** the three targeted couplings are *mixing* couplings (they generate the CKM rotation), so `_unmixed_params` zeroes them and the adiabatic ramp turns them on alongside `transport` — the species-labeling reference stays cleanly block-diagonal; the structural `eta_k3k5_minus` (not a mixing knob) stays on, so v3's mass path is bit-untouched. New regression test `tests/test_quark_v4_lock.py` (19 tests, incl. the unmixed-reference regression); suite 245 passed, 1 xfailed (`v4_library_migration_probe`, PR #164) |
| `C` = inner/outer swap | **Derived** | `C = S: r ↦ 2R_MID − r` involution fixing the throat; reverses mouth normal `n̂ = ±r̂` ⟹ flips Hopf curvature `c₁ → −c₁` (throat → antithroat); `C² = id`, consistent with `T = iσ_y` (B2) and pair-production antithroat (`charge_conjugation_swap_probe`, PR #63) |
| CPT on throat histories | **Assembled** | `q→−, p→+, x→−, s→−, t→−, E→+` with `C²=P²=+I`, `T²=−I`; throat → antithroat run backwards (Feynman–Stückelberg); guaranteed by local Lorentz, global violation `~ 10⁻⁷⁸` (`cpt_assembly_probe`, PR #64) |
| Explicit CPT operator `Θ = −iγ⁵` | **Constructed** | Total spacetime inversion `Θ = γ⁰γ¹γ²γ³ = −iγ⁵`; built from `C = iγ²γ⁰`, `P = γ⁰`, `T = γ¹γ³K`; anticommutes with every `γ^μ` (`j^μ → −j^μ`); matrix `Θ_m² = −I` but antiunitary `Θ² = +I` ((CPT)²=+1) (`cpt_dirac_operator_probe`, PR #65) |
| Throat Dirac 4-spinor from `S_BAM` | **Derived** | Radial `H = −d²/dr*² + V` is a perfect square `A†A + E₀` (SUSY factorization, `W² − W′ = V − E₀`); two SUSY-partner sectors = two wormhole mouths (joined by B3 odd extension); `4 = 2 (mouths) × 2 (SU(2) spin, B2)` = `Ψ_inner ⊕ Ψ_outer` (`throat_dirac_spinor_probe`, PR #66) |
| Even-`k` absence (spin-statistics) | **Classified** | `k mod 2` is the orientability/spin-statistics grading: `T^k` off-diagonal for odd `k` (spin-½ fermion, orientation-reversing) vs diagonal for even `k` (bosonic, orientable double cover); charged leptons = odd class. Not arithmetic — `Φ_avail(k) ≡ 0 mod 2π` for every `k` (`even_k_absence_probe`, PR #67) |
| Throat-to-shell transition | **Demonstrated** | Higher excitations delocalize from the focused lepton-throat pulse into the QCD shell channel (extended-character wavefront); same `S³` closure skeleton, different mode geometry (`throat_to_shell_transition_probe`, PR #68) |
| Shell ↔ QCD structural match | **Partial / Structural** | Shell modes reproduce the documented quark-sector invariants: `Z₂` partition (B2), `3 × 2 = 6` flavors, heavier scale, extended character (`shell_to_qcd_match_probe`, PR #69) |
| Three-generation boundary (sharp `k ≤ 5`) | **Derived / Pinned** | β-uplift quadratic growth `(k−3)²` + throat–shell availability combine to forbid `k ≥ 7`; the sharp `k ≤ 5` cap is the structural three-generation boundary (`three_generation_boundary_probe`, PR #70) |
| `β_lepton = k_5²·(2π) = 50π` | **Derived structurally** | The closure-quantum face of the topological charge: one closure quantum (`2π`) per pair of throat passes (`k_5²`); closes the PR #70 follow-on (`beta_lepton_derivation_probe`, PR #71) |
| `#generations = (k_5+1)/2 = 3` | **Derived structurally** | The linear face of the same `k_5`: number of allowed odd-`k` modes in `{1, 3, …, k_5}` (same primitive as `β_lepton`'s quadratic face) (`three_throat_modes_probe`, PR #72) |
| `k_5 = dim(S³) + 2 = 5` | **Derived structurally** | `k_5 = D_bulk = time + radial + dim(S³) = 1 + 1 + 3 = 5`; `D = 5` is the minimal bulk above 4D giving `f(r) = 1 − (rs/r)²` (squared, matches spin-½ double cover `T² = −I`); reduces "why `k_5 = 5`" to "why the Hopf bundle / S³" (`k5_origin_probe`, PR #73) |
| Quark `n_part = 233` is phenomenological | **Classified** | Extended candidate catalog (Fibonacci, Lucas, Padovan, Perrin, tribonacci, color × flavor × generation, QCD β₀, Tangherlini QCD-shell modes); only baseline coincidences (`F_13 = 233`, `9·k_5²+k_5+3 = 233`), no enumeration survives §8 drift; v3 Hamiltonian is lepton-shaped — wrong machinery for the quark sector (`quark_npart_origin_probe`, PR #76) |
| `n_part` compensates a *dynamical* hierarchy | **Diagnosed (PR #76 sharpened)** | The neutrino arc proved a huge hierarchy can be geometric (the `e^{S}` bounce, ~10⁶), so *size* isn't the obstruction. The quark hierarchy is non-geometric because it is **irregular** (up-type `c/u≈588` vs `t/c≈136` ⟹ not exponential; up/down asymmetric ⟹ not power-law). Geometric shell `ω²(1,n=3,4,5)` carries only ×2.2 of the ×6.4×10⁹ observed mass² span. Quarks are the program's **one dynamical sector**; the lepton↔quark gap `N_q−N_lepton=366` is the dynamical excess `n_part` absorbs (`npart_dynamical_hierarchy_probe`, PR #97) |
| The quark hierarchy is the *flavor puzzle* | **Refined (PR #97 sharpened)** | First step on #97's "right route", testing the mechanism. Quark mass *ratios* are **RG-invariant** (QCD's `γ_m` is flavor-universal ⟹ the common running factor cancels), so the hierarchy is **not** `α_s` running — it is the **flavor puzzle** (the irregular Yukawa couplings, free SM inputs, open across all physics). The quark Yukawas overflow the compressed shell-overtone capacity (mass range ×1.49) by ~×5×10⁴ ⟹ `n_part` compensates; the charged leptons (also a flavor puzzle) instead fit the winding ladder `k∈{1,3,5}` that has the range. BAM captures the quark **structure** (counting), not the Yukawa **magnitudes**. #97 core (dynamical/non-geometric) stands (`quark_hierarchy_flavor_puzzle_probe`, PR #98) |
| QCD confinement: Cornell / flux-tube audit | **Geometric (one scale anchored)** | Cornell `V(L)=σL − A·ℏc/L`: linear `σL` = flux-tube **wormhole bridge** of constant tension; Coulomb = short-distance throat/gluon exchange. **String breaking = Schwinger pair nucleation `exp(−πm_q²/(σL))` = the PR #58 throat-pair mechanism with `eE→σ`** (the string snaps when `σL ≈ 2m_q`). The BAM `σ` reproduces the Regge slope `α'=1/(2πσ)=0.884 GeV⁻²` (obs ~0.88–0.93) and the string-breaking length (~1.4 fm vs lattice 1.35). `√σ ≈ 0.42 GeV` = the single QCD scale anchor (B4 analogue: lepton `m_e` ↔ QCD `√σ`); form geometric, scale calibrated (`qcd_confinement_cornell_audit_probe`, PR #99) |
| Glueballs: pure-confinement benchmark + Möbius tower | **Benchmark + topological prediction** | Closed flux loops (no valence quarks ⟹ no flavor puzzle) are the cleanest confinement probe. BAM orientable ground `√(4πσ)≈1.50 GeV` (3.5√σ) benchmarks lattice 0++ (4.1√σ) to ~13%; closed-string glueball Regge slope = half the meson. **BAM-specific:** the non-orientable **Möbius** sector (`make_mobius_tube`, antiperiodic) gives an *extra* glueball tower (half-integer modes, shifted `+πσ` in `M²`) interleaving the orientable one — ≈2× the states. Glueballs are **not experimentally observed**, so this topological divergence is testable against lattice, not contradicted by experiment (`glueball_closed_flux_loop_probe`, PR #100). **CORRECTED by PR #229:** the `½` moding shift is confirmed (computed three independent ways, shown topological — #100 had only asserted it), so the `+πσ` *level* shift stands; but #100 gave both towers a common intercept `M₀²`, and antiperiodic moding also shifts the zero point by `+π/(4C)` per polarization (two independent regularizations, agreeing to 4.1e-07), which scales as `1/C` and cannot be absorbed into a constant — so the Möbius **ground state is not `M₀² + πσ`** and its absolute masses are open. The heavy Möbius *baryon* search table (#103/#109/#114) is built from `Δ = 2√σ`, not this intercept, and **is unaffected** (`mobius_identification_quantitative_probe`, PR #229) |
| Möbius flux tube ⟹ exotic `J^PC`; observed hybrids match | **Matches data** | Flux-network topology = hadron taxonomy (meson/baryon/tetraquark/pentaquark/hybrid/glueball + Möbius Z₂). A **non-orientable (Möbius) flux tube** carries the antiperiodic phonon that opens the **exotic `1-+`** (forbidden to ordinary qq̄: `P=(−1)^{L+1}`, `C=(−1)^{L+S}`). The observed exotic hybrids `π₁(1600)`, `η₁(1855)` (both `1-+`) match at the right `J^PC` and at `ρ/ω + 2√σ ≈ 1.62, 1.85 GeV`; the tetraquarks (`X, Z_c, T_cc`) / pentaquarks (`P_c`) fit multi-junction networks. **Unlike glueballs, exotics are observed** — so this is where BAM's non-orientable topology meets data, and matches (`mobius_exotic_sector_probe`, PR #101) |
| BAM baryonic exotics: classification + constraints | **Most-constrained corner** | Unlike mesons (smoking-gun `1-+`), **baryons have no forbidden `J^P`** (`P=(−1)^L`, `S∈{½,3/2}`, no `C`) — so BAM's Möbius/hybrid baryons are **supernumerary ordinary-`J^P`** states, identifiable only by counting. They sit in the light N*/Δ* region (`nucleon/Δ + 2√σ ≈ 1.79, 2.08 GeV`), the densest, best-measured spectrum — the **most experimentally constrained** corner of BAM's non-orientable predictions (opposite extreme from glueballs). The Möbius doubling must coincide with observed resonances or decouple (`πN`), else be excluded. Constraint ranking: light N*/Δ* > strange hyperons > charm/bottom baryons (freest) (`baryonic_exotics_classification_probe`, PR #102) |
| Heavy-quark Möbius baryon: prediction in the freest channel | **Findable / unconstrained** | By heavy-quark symmetry (heavy quark = spectator) the Möbius/flux gap `Δ=2√σ≈0.85 GeV` is **flavor-independent** (same for c and b) — the cross-flavor signature replacing the absent exotic-`J^P`. Predictions: Λ_c **~3.14**, Ω_c ~3.54, Λ_b **~6.47**, Ω_b ~6.89, Ξ_cc ~4.47 GeV — all just **above** current excitation ceilings (findable at LHCb/Belle II, not excluded) and above the orbital tower. Doubly-heavy `Ξ_cc` and `Ω_b` have no measured excitations → entirely unconstrained. Exact mass (lattice hybrid gap 0.8–1.3 GeV) / `J^P` open (`heavy_mobius_baryon_probe`, PR #103) |
| Heavy Möbius baryon: decay channels + search strategy | **Twist-unwinding → hybrid selection rule (falsifiable)** | Completes #103: how the state decays and how to find it. Decay = **twist-unwinding** (non-orientable `−1` → orientable `+1` ground state sheds `2√σ` as light isoscalar hadrons), so it inherits the flux-tube **hybrid selection rule**: single-S-wave-π-to-ground **SUPPRESSED**; `Σ_Q π` / isoscalar dipion `Λ_Q(ππ)` / P-wave+π **PREFERRED** — the branching **pattern** that distinguishes it from a radial excitation (which does the opposite). Cross-flavor clincher: all-light Q-values **identical** for c and b (`Λ_Q ππ` **569**, `Λ_Q η` **301** MeV; `Σ_Q π` offset only by hyperfine 167/194). Broad (~tens–150 MeV, open channels) → best in LHCb/Belle II amplitude analyses of `Λ_Q ππ`, `Σ_Q π`, `DN`/`BN` (`Ξ_cc`/`Ω_b` wide open). Branching fractions / width / `J^P` open (`heavy_mobius_baryon_decay_probe`, PR #109) |
| Non-orientable sector: compact **experimental note** | **Compiled (reference card)** | Consolidates the whole Möbius / closed-flux-loop sector (PRs #100–#109) into one LHCb/Belle II/BESIII-style note — predicted masses, Q-values, preferred/suppressed modes, analysis handles — every number a pushforward of the single input `√σ`. **Masses:** mesonic `1⁻⁺` π₁ **~1.62**, η₁ **~1.85** GeV (matched to π₁(1600)/η₁(1855)); glueball `0⁺⁺` `√(4πσ)` ~1.50 GeV (unobserved, freest); heavy Möbius baryons Λ_c 3135 … Ω_b 6894 MeV. **Decays:** twist-unwinding → hybrid selection rule (single-π-to-ground suppressed), cross-flavor Q-match (`Λ_Q ππ` 569, `Λ_Q η` 301 MeV identical c=b). **Handles:** branching pattern vs radial, isoscalar high-`m(ππ)` dipion, broad→amplitude fits, `1⁻⁺` smoking gun (mesons). Standalone at `docs/bam_nonorientable_experimental_note.md` (`nonorientable_experimental_note_probe`, PR #110) |
| Heavy Möbius baryon: sharper **LHCb / Belle II search table** | **Tiered, actionable** | Converts #109/#110 into a ranked search table. **New handle:** the `Λ_Q(ππ)` **dipion endpoint** `m(ππ)_max = M_Möbius − M_ground = 2√σ ≈ 849 MeV` is **flavor-independent** (same edge above charm and bottom, peaking high) — a fixed edge in a directly-plotted observable, one overlay tests the framework. **Tier 1** (discovery pair): Λ_c (3135, `Λ_c⁺π⁺π⁻`, `Λ_c⁺→pK⁻π⁺`, LHCb+Belle II) + Λ_b (6469, `Λ_b⁰π⁺π⁻`, LHCb b-decays) — the cross-flavor clincher. **Tier 2** (unexplored, rare): Ξ_cc (4471, `Ξ_cc⁺⁺→Λ_c⁺K⁻π⁺π⁺`), Ω_b (6894). **Tier 3** (calibratable): Ω_c (3544, above 2017 excitations). Discriminators: suppressed single-π-to-ground, 849 MeV endpoint, cross-flavor Q-match. Standalone at `docs/heavy_mobius_baryon_search_table.md`; masses ±band / broad / BFs / `J^P` open (`heavy_mobius_baryon_search_table_probe`, PR #114) |
| **Program-wide synthesis: the input budget** | **Capstone** | Classifies every result into 5 epistemic tiers. **The whole dimensionful content reduces to 2 B4 anchors** — `m_e = ℏc/R_MID` (QED/lepton) and `√σ ≈ Λ_QCD` (confinement) — the irreducible minimum (one scale/sector, PR #52). Open dimensionless inputs are localized to 2 (neutrino compliance `ε`, quark `n_part`); the only other open input is the **universal flavor puzzle** (Yukawa hierarchy — not BAM-specific). The APS partition audit (PRs #123–#125) sharpens the status of `n_part`: it is **not** an unexplained compensator but the **unique matter-partition residual after APS reduction** — the one feeding integer the index machinery cannot derive (leptons `N_lepton=4·k₅²=100` are fully derived from the bulk dimension `k₅`, quarks `N_q=2·n_part` keep `n_part`). The rest is ~22 derived-geometry results + 6 non-orientable topological predictions (matched → falsifiable → findable → free) (`program_synthesis_probe`, PR #104) |
| **α and G in the ledger** | **G = anchor, α = universal residual** | **G** is the dimensionful **anchor** — the GR-foundational scale (the throat's size, the one B4 length, set by bulk gravity `λ_crit=√(6\|Λ₅\|)/κ₅²`, PR #57) and the root the #104 sector anchors `m_e`/`√σ` descend from. **α** is a **universal residual** — used as input everywhere (`A_EM=α·ℏc/2`, `a=α/2π`); BAM derives the charge unit `\|c₁\|=1`, the `1/2π` measure, and α's *running*, but the *value* 1/137 is a free input (the "137 problem"), sitting with the flavor puzzle. **ℏ** is geometric (the closure quantum, `ℏ=m_e·R_MID·c`); **c** is units (`alpha_G_ledger_classification_probe`, PR #105) |
| **How many scales? `m_e` vs `√σ`** | **Not independent — one G + an underived ratio** | `m_e` and `√σ` both descend from the single bulk-gravity scale `G` (PR #57: `R_MID` and `σ` from `λ_crit=√(6\|Λ₅\|)/κ₅²`), so the **dimensionful-anchor count reduces 2→1**. But their ratio `√σ/m_e ≈ 830` (the lepton-throat / QCD-confinement hierarchy) is **underived** — no clean closure number (nearest `50π·k_5=785`, 5.4% off, a near-coincidence like `F_13=233`). So it's a **repackaging, not a free reduction**: a dimensionful anchor becomes a dimensionless residual (joining `ε`, `n_part`, `α`), total inputs unchanged. The gain: the sole fundamental *scale* is now `G` (`scale_count_anchors_probe`, PR #106) |
| Is `832 = N_q+ΔN` an independent ratio, or recycled `n_part`? | **Recycled n_part (negative result)** | A tempting candidate derivation of the #106 ratio: `N_q+ΔN = 2N_q−N_lepton = 832 ≈ √σ/m_e ≈ 830` (0.2%). **Rejected.** `832 = 4·n_part − 4·k_5²` is built from the `n_part` compensator. Decisive §8-drift test: propagating `n_part∈{216..255}` makes "832" drift **764–920 (±9%)** while 830 is fixed → a baseline coincidence (like `50π·k_5=785`, `F_13=233`). No independent bulk shell-stress integral selects ~466/832 (natural ones are `O(10–70)`); 466 enters only via the v3 fit. Circular. `√σ/m_e` stays underived; the #106 ledger is unchanged (`ratio_832_npart_recycling_probe`, PR #107) |
| The legitimate search: does any fit-independent, §8-stable bulk quantity select `√σ/m_e ≈ 830`? | **No — search fails; ratio plausibly irreducible** | Ran the fit-independent route #107 called for: quantities built **only** from fixed geometry (`k_5=5`, `β_lepton=50π`, `2π`), scored against 830.3 under 4 criteria (C1 fit-independent, C2 §8-stable, C3 <1%, C4 principled). **C2 is automatic** for geometric candidates (they never touch the quark ablations). But C3∧C4 fail: best **principled** candidate `2π·k_5³ = β_lepton·k_5 = 785.4` (**−5.4%**); every sub-% match needs an ad-hoc factor (`π·265`, `(4/3)·k_5⁴`, `k_5⁵/3.77` — 265, 4/3, 3.77 reverse-engineered). Exponential route: `ln(830)=6.72` vs clean action `2π=6.28` (7% off). Cavity integrals `O(10–350)`, select nothing near 830. **`√σ/m_e` stays UNDERIVED — now plausibly IRREDUCIBLE, like `α`.** BAM does **not** collapse to a single anchor: one scale `G` + this ratio + `α` + the flavor puzzle (`lepton_qcd_ratio_legitimate_search_probe`, PR #108) |
| Shell waveguide basis + operator scaffold | **Constructed** | Quarks reframed as cavity wavefronts that resolve the shell (NOT throat traversals). 6-state `(l, n, p)` basis with `H = H_kin + H_Z2 + H_couple`; `H_kin = ω²(l, n)` cavity-eigenfrequency-squared, not the lepton `β·k²·(2π)` winding cost (`qcd_shell_waveguide_scaffold_probe`, PR #77) |
| Shell mass-ordering / `n_part` audit | **Sharpened** | Shell basis structurally better than v3 in 4 ways (cavity wavefronts; ω² kinetic; Z₂ partition slot; 6 flavors). Uniform `χ·σ_z` cannot reproduce within-generation inversion (best 2/3 blocks); sign-flipping χ_n can (existence proof). Coverage gap: shell kinetic ×2.2 vs observed ×6.4·10⁹ — `n_part` NOT resolved at #78 alone (`shell_mass_ordering_audit_probe`, PR #78) |
| Boundary-stress `χ_n` + singlet placeholder | **Derived structurally** | `χ_n = T_odd(n) = (T_inner − T_outer)/2` from Z₂-antisymmetric piece of cavity-mouth boundary stress (PR #63's inner/outer swap). NO free parameter once cavity geometry fixed. Uniform-positive sign (no flip), shell-suppressed magnitude — 30–100× too small for observed splittings; PR #78 sign-flipping ansatz overruled (`boundary_stress_chi_n_probe`, PR #79) |
| BAM-native color algebra = `SU(2) × Z₂` | **Identified** | SU(2) from B2 / Hopf holonomy (PRs #59–#66; `T = iσ_y`, `T² = −I`) + Z₂ from PR #63 inner/outer swap. SU(2) acts on partition index; Z₂ swaps n=3 ↔ n=5. SU(3) NOT BAM-derivable from current scaffold (all natural triplets give SO(3)/SU(2)); Pati-Salam SU(4) requires throat↔shell algebra map (open). v3 species map revised: `+ = heavier` uniformly. Inter-generation hierarchy outside BAM color scope; `n_part = 233` residual with sharply identified scope (`color_algebra_shell_probe`, PR #80) |
| Throat ↔ shell `n + 3` Pati-Salam bridge | **Built (partial)** | Each generation has a lepton at `n = g−1` (throat) and a quark-pair at `n = g+2` (shell); shift `+3` = PR #68 shell threshold (no free parameter). Unified 12-state `(l, n, p)` basis + throat-shell Z₂. Full SU(4) PS needs 3 open extensions: BAM-native neutrinos, 3-fold quark color, lepton-quark mass-operator unification (`pati_salam_throat_shell_bridge_probe`, PR #82) |
| **Lepton + quark masses = ONE Bohr-Sommerfeld operator** | **Unified** | `m²(k, n) = (k·2π/L_throat)² + ((n+1)·π/L_cavity)²`, `L_throat = √(2π)/k_5`. Lepton `β·k²` (PR #71) and quark `ω²(l,n)` (PR #77) are the same operator `m² = (S/L_eff)²`. Cavity Bohr-Sommerfeld `∮√(ω²−V)dr* = (n+1)·π` verified to machine precision; `(2π/L_throat)² = k_5²·(2π) = 50π = β_lepton` recovered. `k = 0` for quarks = "don't pass through the throat"; closure quanta `2π` (throat) vs `π` (cavity) = BAM full/half-cycle (`throat_shell_mass_operator_unification_probe`, PR #83) |
| `(k≠0, n≥3)` quadrant = leptoquark sector | **Mapped** | The unified `(k, n)` operator's fourth quadrant (winding **and** shell-saturated) is the leptoquark sector, completing the four-quadrant reading: lepton `(k≠0, n<3)`, quark `(k=0, n≥3)`, neutrino `(k=0, n<3)`, leptoquark `(k≠0, n≥3)` (`winding_shell_quadrant_probe`, PR #85) |
| Neutrino = Majorana (seesaw) | **Derived structurally** | `k=0 ⟹ c₁=0 ⟹ C-invariant` (PR #63) ⟹ neutrino is its own antiparticle ⟹ **Majorana**; suppression = seesaw `m_ν = m_D²/M_R`, available **only** to the chargeless sector (charged leptons `c₁=±1` are Dirac and keep `β·k²`) — explains why only ν is light; required `M_R ≈ 0.3–1.8 TeV` open (`neutrino_quadrant_suppression_probe`, PR #86) |
| Seesaw scale `M_R` from throat-nucleation tunnelling | **Grounded / scale recast** | `ΔL=2` Majorana = PR #58 throat↔antithroat (antipodal `Z₂`) transition; PR #58's `Σc₁=0` on a single state **is** PR #86's only-neutrino rule. `M_R` ≠ barrier height (`E_c ≈ 2.8 keV`, ~10⁸ too small); suppression = tunnelling through the barrier `m_ν = m_D·e^{−S}` ⟹ `M_R = m_D·e^{S}`, recasting the open ~TeV scale as a modest, generation-stable bounce action `S ≈ 15–18` (the PR #58 instanton follow-on) (`seesaw_scale_nucleation_compliance_probe`, PR #87) |
| Majorana bounce `S` = non-orientable tortoise log | **Sharpened / open** | Reduced Euclidean bounce `S = √(2 μ E_c)·L*(ε)` on the odd (`c₁→−c₁`) tortoise path: the tortoise coord diverges logarithmically at the throat ⟹ **rigid throat = massless ν** (compliance `ε` is the mass-generating parameter), and `S ∝ ln(1/ε)` is naturally `O(10)`/gen-stable — the form PR #87 required. But the EM-throat tension **under-produces** by ~40× (`S ≲ 1`); `S ≈ 15–18` needs a `ΔL=2` (B−L) tension `~6–12×` stiffer. Open input localised: ~TeV mass (#86) → `O(15)` action (#87) → `O(10)` tension ratio (#88) (`majorana_bounce_action_probe`, PR #88) |
| ΔL=2 / B−L tension ratio `t` bracketed | **Constrained** | The `ΔL=2` flip reverses orientation (`c₁→−c₁`) ⟹ a **global** operation, so `t` is a global-closure enhancement of the **local** EM surface tension. Bracketed parameter-free by the **closure quantum `2π`** (minimal orientation reversal, lower) and the **winding action `k_5√(2π) = √β_lepton`** (full winding, upper): `t ∈ [2π, k_5√(2π)] ≈ [6.28, 12.53]` — exactly PR #88's required `6–12` (computed `[6.41, 12.05]` sits inside). Residual = where in the window = compliance `ε`; `m_charged/m_D ≈ 11.9 ≈ √β` cross-check (`b_minus_l_tension_ratio_probe`, PR #89) |
| Boundary compliance `ε` from bulk geometry → `m_ν` scale | **Chain closed (order-of-mag)** | `ε` is the chargeless throat's sub-throat **healing length** (`ε = ℓ²/2rs` from the neck warp `f≈2(r−rs)/rs`); sub-throat *for the neutrino* because the `c₁=0` neck is not EM-propped (the charged `c₁=±1` neck is, and stays Dirac). Natural BAM scales (`R_c³, Δ³, (m_D/m_ch)²`) land `ε` in the PR #89 window; with the winding-edge tension `t≈√β` (cross-check-favoured) the chain gives `S ≈ 15–19`, **`m_ν ~ few meV`** — the observed scale, untuned (`2π` edge gives `S≈4`, too small). The full chain `~TeV → S → t → window → ε → meV` is closed; precise `m_ν` / generation spread residual (`boundary_compliance_bulk_geometry_probe`, PR #90) |
| Is `ε` computed from bulk compliance, or inferred from meV? | **Smallness derived; precise value residual** | Sharpens PR #90's question. **Computed (meV-free):** the neck healing length `ℓ ~ R_c = 2σ/ρ` (with `σ,ρ` from the **electron** calibration PR #58, `R_c = 2/9`) gives `ε ~ R_c³ ≈ 0.011` — sub-throat, `O(10⁻²)`, no neutrino input. With `t = k_5√(2π) = √β_lepton` (PR #89), `S ≈ 16.85` ⟹ `m_ν ≈ 2.1 meV` — the meV **scale** *output* (retrodiction), structurally deriving the lightness (`ε≪1 ⟹ S large ⟹ m_ν = m_D e^{−S}` tiny). **Residual:** the *precise* `ε`. Since `m_ν ∝ ε^{4.8}`, the `O(1)` ambiguity (`R_c³`→2, `Δ³`→20, `R_c²/2`→108 meV) spans ×50; the absolute compliance normalization is the unpinned `κ₅²/Λ₅` (only `√6` fixed, PR #57). **So the smallness is derived from bulk compliance; the exact value is not** (`epsilon_bulk_compliance_probe`, PR #112) |
| Generation-dependent `ε_n` and the hierarchy spread | **Direction derived; magnitude overshoots → residual** | Tests PR #91's fix for the spread PR #112 left open. Generations = cavity overtones `n`; the overtone boundary stress `χ_n` (PR #79) decreases (0.304, 0.097, 0.039), so `ε_n ∝ 1/χ_n` (compliance = 1/stiffness). **Direction right:** `ε_n` increases with `n` ⟹ less suppression ⟹ heavier ⟹ **normal ordering**, untuned. **Magnitude overshoots:** the observed spread needs gentle `ε_n` ratios `(1, 1.18, 1.57)`, but `1/χ_n` gives `(1, 3.13, 7.79)` ⟹ `m_ν3/m_ν2 ≈ 162` vs observed 5.85 (**×28**). Cause: the steep bounce (`m_ν ∝ ε^{4.8}`, PR #112) amplifies the ×8 `χ_n` variation into ~10⁴ in mass; the required power `p ≈ 0.15–0.31` (≠ principled 1). So `ε_n` **accommodates** the spread (fit) but does not **predict** it — the spread stays a residual, plausibly the mixing/anarchy sector (PR #92) (`generation_dependent_eps_n_probe`, PR #113) |
| Generation spread + `PMNS ≫ CKM` from channels | **Structural** | Generations = cavity overtones ⟹ bare `m_ν ∝ m_D` (normal ordering `1:1.87:2.74`); the spread is widened in the right direction by the overtone-dependent neck coupling (PR #79 `χ_n` ↓ with `n` ⟹ higher-`n` less suppressed ⟹ heavier). **Headline:** large PMNS vs small CKM is the **cross-channel** (leptons: charged throat-winding `k≠0` × neutrino cavity-resolving `k=0`) vs **intra-channel** (quarks: up & down both cavity-shell `k=0`) distinction — the BAM reason `PMNS ≫ CKM`. Exact angles/spectrum open (`generation_spread_pmns_mixing_probe`, PR #91) |
| PMNS anarchic, CKM aligned (quantitative) | **Tested** | A naive radial mode overlap gives near-permutation (small) mixing — so large PMNS is **not** a literal overlap. The lepton generation labels live in **different coordinates** (charged: closure-winding `k`; neutrino: radial-overtone `n`) ⟹ no alignment ⟹ **anarchic** map. Observed PMNS (33.4°, 49°, 8.6°) is **typical** of a Haar-random `U(3)` (30th/57th/4th percentile); CKM (13°, 2.4°, 0.2°) is **extremely atypical** (joint `p ≈ 0`) = aligned (up & down share the radial-overtone coordinate). PMNS ∈ anarchy class, CKM ∈ aligned class; specific angles not pinned (θ13 mild tension) (`cross_channel_pmns_overlap_probe`, PR #92) |
| θ13 suppressed by residual alignment | **Tension resolved** | θ13 = `U_e3` is the corner / most coordinate-distant (**two-hop**) element (lowest winding `k=1` × highest overtone `n=2`, gap 2); θ12, θ23 are adjacent (gap 1). The throat↔shell coupling is **local** in the `(k,n)` lattice (PR #82 `+3` shift, PR #83 operator), so the corner `U_e3` is a suppressed two-hop amplitude — a residual **nearest-neighbour** alignment. A structured-anarchy model with `μ≈3` makes θ13 robustly the smallest angle (frac 0.50→0.72) and moves observed θ13=8.6° from the 4th to ~21st percentile (PR #92 tension resolved), θ12/θ23 staying typical. Exact θ13 (μ; median saturates ~14–16°) open (`theta13_residual_alignment_probe`, PR #93) |
| CP violation generic; two Majorana phases exist | **Structural** | CP violation is **generic**: the winding amplitudes carry the complex Hopf holonomy `e^{ikχ}` (PR #60), so the PMNS is generically complex (`δ_CP ≠ 0, π`; CP conservation is measure-zero). The **Jarlskog dichotomy** mirrors the angles: `|J_PMNS| ≈ 0.026` is typical of anarchy (51st/81st percentile, large CP violation), `|J_CKM| ≈ 3×10⁻⁵` is extremely atypical (~0.1th, aligned/suppressed). **Two Majorana phases exist** because the neutrino is Majorana (`c₁=0`, PR #86) — CP phases of the ΔL=2 throat↔antithroat sector, observable in 0νββ; Dirac would have none. Specific values anarchic/not pinned (`cp_majorana_phase_probe`, PR #94) |
| 0νββ effective mass `m_ββ ≲ 8 meV` (falsifiable) | **Predicted** | Combines the arc: 0νββ **occurs** (neutrino Majorana ⟸ `c₁=0`, PR #86); **normal ordering** (PR #91) selects the NO band; **anarchic Majorana phases** (PR #94) populate it incl. cancellation to ~0; the **light scale** (PR #90, ~few meV) gives `m_ββ ≲ 8 meV`. Below current bound (KamLAND-Zen 28–122 meV, null result expected) and largely below next-gen reach (~9–20 meV), and below the inverted-ordering floor (~19 meV). **Falsifier:** a discovery at `m_ββ ≳ 19 meV` ⟹ inverted/degenerate, contradicting BAM (`zeronubb_effective_mass_probe`, PR #95) |
| Cosmological `Σm_ν ≈ 59–65 meV` (falsifiable) | **Predicted** | The same light, normal-ordered spectrum fixes `Σm_ν = m1+m2+m3`: the NO floor is `√Δm²_21 + √Δm²_31 ≈ 58.7 meV` (IO floor ≈ 99 meV), and the light scale (PR #90) keeps `Σm_ν ≈ 59–65 meV`, pinned near the floor. Consistent with Planck (<120 meV), just inside DESI DR1+CMB (<72 meV), right at the DESI DR2+CMB frontier (~60–64 meV). **Falsifier:** robust `Σm_ν < 58.7 meV` ⟹ NO excluded; `Σm_ν ≳ 100 meV` ⟹ not light. Cross-checks the 0νββ prediction (one spectrum) (`cosmological_sigma_mnu_probe`, PR #96) |
| meV-scale spectrum **sharpened** (NuFIT 6.0 + DESI DR2) | **Pinned; only Σm_ν testable** | Sharpens the #96 band into a full pinned spectrum. NuFIT 6.0 fixes `m₂ = 8.65`, `m₃ = 50.34 meV` (NO floor `Σm_ν = 59.0`); DESI DR2 + CMB (≲60–64 meV) corners `m₁ ≲ 3 meV` ⟹ **`Σm_ν ∈ [59.0, 62.6] meV`** (tightened from 59–65, toward the floor). Lab effective masses: `m_β ≈ 8.8–9.3 meV`; **`m_ββ` has a nonzero floor `[1.5, 3.7] meV`** — NO contributions can't fully cancel (`s12²c13²m₂ = 2.60 > s13²m₃ = 1.10 meV`). **Honest reachability:** only `Σm_ν` is near-term testable (DESI, at the floor now); `m_β` ~4–5× below Project 8, `m_ββ` ~3–10× below LEGEND-1000/nEXO. Flag: some 2025 DESI+CMB fits prefer Σm_ν at/below the floor → tension for all NO models (`neutrino_mev_scale_sharpening_probe`, PR #111) |

### Research goals (not yet fully derived)

| Physics | Proposed geometry |
|---------|-------------------|
| Electromagnetism | Curvature of the Hopf connection on S³ |
| Charged-lepton ladder (e, μ, τ) | Eigenvalues of a k-pass instanton-transition matrix with S³ action base `2π` and k=5 uplift `200π` — **sub-percent fit achieved** |
| Particle mass (general) | One Bohr-Sommerfeld closure operator `m² = (S/L_eff)²` over both fermion sectors: leptons = throat-winding (`k ≠ 0`), quarks = cavity-resolving (`k = 0`); inter-generation hierarchy still open (PR #83) |
| QCD confinement | 1D flux-tube network with bridge nucleation — Cornell `σL−A/L` audited (PR #99): flux tube = wormhole bridge, string breaking = PR #58 Schwinger throat-pair (`eE→σ`); `√σ` the one QCD anchor |
| Retrocausal photon exchange | Wheeler–Feynman absorber theory on S³ |
| Black-hole interior | Coherent condensate of non-orientable wormhole throats |
| Bell correlations | Non-orientable throat transport + Hopf SU(2) projection |
| Entanglement = wormholes | Bell correlations from throat connectivity |
| Quantisation from resonance | S³ antipodal cavity selecting discrete spectrum |
| Topological censorship | Non-orientable throats evading standard no-go theorems |
| QFT event reinterpretation (Compton) | Antipodal `S³` Green function as propagator + Hopf-fibre photon polarisation + closed-form vertex resummation reproducing Klein-Nishina exactly — see [QFT-event-reinterpretation thread](#qft-event-reinterpretation-thread-compton-scattering) below |

## Package Structure

```
geometrodynamics/
├── geometrodynamics/
│   ├── constants.py          # Shared physical & simulation constants
│   ├── hopf/                 # Hopf fibration on S³
│   │   ├── connection.py     # A = ½cos(χ)dφ, curvature, holonomy
│   │   ├── chern.py          # First Chern number c₁ = 1
│   │   └── spinor.py         # SU(2) spinor transport (spin-½)
│   ├── tangherlini/          # 5D wormhole eigenmodes
│   │   ├── radial.py         # Chebyshev spectral solver
│   │   ├── maxwell.py        # Sourced Maxwell BVP (Coulomb validation)
│   │   ├── alpha_q.py        # Throat flux ratios (no free parameters)
│   │   └── lepton_spectrum.py # Locked e/μ/τ instanton-transition matrix
│   ├── transaction/          # Wheeler–Feynman absorber theory on S³
│   │   ├── particles.py      # ThroatMode, MouthState, Particle4, GravWave
│   │   ├── s3_geometry.py    # Geodesics, Green function, antipodal map
│   │   ├── handshake.py      # Offer/confirm/transaction protocol
│   │   └── cavity.py         # CavityMode, CavityPacket, AntipodalCavity
│   ├── embedding/            # Non-orientable throat topology
│   │   ├── topology.py       # ThroatDefect, ConjugatePair, transport ops
│   │   └── transport.py      # T = iσ_y derived from Hopf fibration
│   ├── bell/                 # Bell correlations from geometry
│   │   ├── pair_state.py     # BellPair with cavity history evolution
│   │   ├── analyzers.py      # Detector settings as SU(2) operators
│   │   ├── correlations.py   # E(a,b), CHSH, no-signaling
│   │   ├── hopf_phases.py    # Bell closure phases from Hopf holonomy
│   │   └── bulk_identity.py  # Kinematic Bell from shared bulk topology
│   ├── history/              # Closed-history framework (unifying backend)
│   │   └── closure.py        # Events, Worldlines, History, branch enumeration
│   ├── qcd/                  # Geometrodynamic QCD
│   │   ├── constants.py      # σ, α_s, ℏc, SAT parameters
│   │   ├── color.py          # SU(3) color algebra, generators
│   │   ├── bridge.py         # BridgeField, Cornell potential
│   │   ├── network.py        # Node, Branch, Junction, HadronicNetwork
│   │   ├── topology.py       # Meson, baryon, glueball, hybrid, …
│   │   ├── solver.py         # Störmer–Verlet + SAT boundaries
│   │   ├── spectrum.py       # Möbius modes, throat–branch crosswalk
│   │   └── diagnostics.py    # String tension, mode shifts, calibration
│   ├── blackhole/            # Black holes as wormhole-throat condensates
│   │   ├── condensate.py     # CoherentCondensate, ThroatState, constructors
│   │   ├── interior.py       # Hayward regular metric, geodesics, horizons
│   │   ├── entropy.py        # Bekenstein-Hawking from throat counting
│   │   └── derivation.py     # Condensate → metric via Einstein equations
│   └── viz/                  # Watchable classical wave studies
│       ├── antipodal_focusing.py  # open disperses vs closed refocuses (#166)
│       ├── antipodal_crossing.py  # antipodal crossing / absorption events
│       ├── throat_wavefront.py    # bare S2 / sealed mouths / open catenoid
│       │                          #   throat, on one clock (#242)
│       ├── radial_caustic.py      # focal geometry: which fronts fold (#243)
│       ├── warped_sphere.py       # one continuous S2 warped by the solved
│       │                          #   wave, nested between two shells
│       ├── spin2_tidal.py         # the tensor case: spin-2 wave as tidal
│       │                          #   shear, ring at the focus, no l<2
│       ├── embedded_wave.py       # h_ab as a continuous R^3 deformation:
│       │                          #   shape = -Lap(W)/2, shear = TF Hess W
│       └── geometry_panels.py     # Hopf / throat / Green / handshake panels
├── tests/                    # pytest validation suite
├── notebooks/                # Jupyter notebooks (per-topic)
├── scripts/                  # Lepton-ladder calibration CLIs
├── docs/                     # Lepton axioms + scan archaeology
└── pyproject.toml
```

## Installation

```bash
git clone https://github.com/davidmdrpi/geometrodynamics.git
cd geometrodynamics
pip install -e ".[dev]"
```

## Running the Validation Suite

```bash
# All tests (fast)
pytest

# Include slow tests (bridge nucleation, string tension scans)
pytest -m ""

# Skip slow tests
pytest -m "not slow"
```

## Lepton mass ladder (e, μ, τ) from a locked S³ action

The lepton surrogate now ships with a **fully locked topological baseline**
that reproduces all three charged-lepton masses to sub-percent accuracy with
**zero free parameters at scan time** — only the electron mass is used to set
the overall MeV scale.

### Locked axioms

- `action_base = 2π`  — the S³ great-circle action (circumference invariant).
- `k_uplift_beta = 50π`  — k-selective uplift coefficient.
  For `k=5`, the uplift is `4·β = 200π`, i.e. **exactly 100 × (2π)** S³
  winding quanta.
- `winding_mode = "max"`  — off-diagonal tunneling cost scales with the deeper
  branch, `Δk = max(kᵢ, kⱼ)`.
- `depth_cost_mode = "tunnel_only"`  — the S³ base action enters only through
  the tunneling suppression, not as an additional diagonal offset.
- `resistance_model = "exponential"`  — re-entry cost `κ·(eᵏ − 1)` captures
  exponential geometric writhe/curvature build-up with generation depth.
- Baseline anchor `(phase, transport, pinhole, resistance) ≈
  (0.001, 25.1, 22.5, 0.217869)`. As of the closure-ledger sequence
  (`docs/hbar_origin_note.md`), all four are now identified with
  closure-quantum / Tangherlini-grid invariants:
  `transport = 8π`, `pinhole γ = Σ V_max[1..5]`,
  `resistance = 7π/100`, with the phase channel decoupled.

The generation-block diagonal takes the form

```
H_kk = action_base + resistance_scale · k²  +  res_diag(k)
                  +  pinhole(k ∈ {3, 5})   +  β · max(0, k−3)²
```

and off-diagonals are `−transport · exp(−α_eff · Δk) · cos(phase · Δk)`.
See `docs/lepton_axioms.md` for the full matrix construction.

### Validated predictions (locked baseline, no tuning)

| Lepton | k | Predicted (MeV) | Observed (MeV) | Relative error |
|--------|---|-----------------|----------------|----------------|
| e      | 1 | 0.510999        | 0.510999       | 0.0000% (anchor) |
| μ      | 3 | 105.61260       | 105.65838      | **0.0433%** |
| τ      | 5 | 1778.93809      | 1776.86        | **0.1170%** |

Muon/electron ratio: predicted **206.6787**, observed **206.7683**
(relative error **4.33 × 10⁻⁴**).

Reproduce directly from Python:

```python
from geometrodynamics.tangherlini import solved_lepton_masses_mev
masses = solved_lepton_masses_mev()           # read-only np.ndarray
print(masses)   # [0.51099895, 105.6126..., 1778.9381...]
```

Or by CLI (no `PYTHONPATH` needed):

```bash
python scripts/lock_beta_50pi_probe.py --n-points 32
```

which additionally pins `β = 50π` exactly and optimizes only the four
sub-leading knobs; it reports `mu/e` error ≈ 1 × 10⁻⁶% and
`τ` relative error ≈ 0.161%.

### Geometric implications

1. **Three generations correspond to odd pass depths `k = 1, 3, 5`.** The
   ladder is labelled by the number of S³ passes before closure; the locked
   baseline scans exactly these three depths. **Even-`k` absence is now
   classified** as a spin-statistics selection rule (`even_k_absence_probe`,
   PR #67): `k mod 2` is the orientability/spin-statistics grading
   (`T^k` off-diagonal for odd `k` = orientation-reversing closure across
   the non-orientable throat = spin-½ fermion; diagonal for even `k` =
   orientable double cover = bosonic). Charged leptons are spin-½, hence
   the odd class. The sharp upper bound `k ≤ 5` is the
   **three-generation boundary** (`three_generation_boundary_probe`, PR
   #70), and `k_5 = 5 = D_bulk = dim(S³) + 2` is the BAM bulk dimension
   (`k5_origin_probe`, PR #73), with `β_lepton = k_5²·(2π) = 50π`
   (`beta_lepton_derivation_probe`, PR #71) and `#generations = (k_5+1)/2 = 3`
   (`three_throat_modes_probe`, PR #72) both derived from the same `k_5`
   primitive.
2. **τ uplift is exactly 100 quanta of the S³ action.** The k=5 uplift is
   `4β = 200π = 100·(2π)`, a pure integer multiple of the great-circle action
   `2π`. No tuning is required; removing the integer lock degrades `τ` by an
   order of magnitude (see `docs/lepton_tau_target.md`).
3. **The μ/e ratio is a structural eigenvalue ratio, not a coupling.** With
   `action_base = 2π` locked and the exponential resistance profile, the
   calibration scan finds exact μ/e roots on a broad resistance basin
   (±1% resistance keeps `mu_err` < 8%), replacing the earlier
   "attractor needle" regime (see `docs/lepton_tau_target.md`, "Hard S³ lock
   experiment").
4. **Quadratic diagonal `∝ k²` plus quadratic uplift `∝ (k−3)²`** together
   reproduce the observed `m_e : m_μ : m_τ ≈ 1 : 207 : 3477` hierarchy: the
   `k²` term sets the `μ/e` split and the `(k−3)²` term independently lifts
   the τ sector without disturbing the `μ/e` root.
5. **Tunneling-side depth cost dominates diagonal depth cost.** The ablation
   scan showed `tunnel_only` outperforms `diag_only` by nearly 2× on best
   μ/e (see `docs/lepton_ablation_results.md`) — consistent with a picture in
   which the inter-generation transition amplitude, not the on-generation
   mass term, sets the ratio.
6. **A `max` winding rule beats a `delta` winding rule.** Setting
   `Δk = max(kᵢ, kⱼ)` (rather than `|kᵢ − kⱼ|`) in the tunneling action was
   the change that first pushed `μ/e` from ~10 toward the experimental
   ~206.77, because it penalises transitions into deeper branches by the full
   target winding — a topological-cost interpretation consistent with the S³
   action base.

### Script map

| Script | Purpose |
|--------|---------|
| `scripts/calibrate_muon_ratio.py` | Coarse grid; solves resistance for exact μ/e root at each (phase, transport, pinhole). |
| `scripts/sweep_k_uplift_beta.py`  | Sweeps `β` with exact μ/e enforced; locates best τ fit. |
| `scripts/map_basin_k_uplift.py`   | Local gradient probe around an exact-μ/e point; reports basin width. |
| `scripts/refine_locked_tau.py`    | Dense locked scan with action_base fixed to 2π; reports integer-winding β family. |
| `scripts/lock_beta_50pi_probe.py` | Hard `β = 50π` lock; optimizes only (phase, transport, pinhole, resistance). |

See `docs/lepton_ablation_results.md`, `docs/lepton_tau_target.md`, and
`docs/lepton_next_steps.md` for the full scan archaeology, and
`docs/hbar_origin_note.md` for the closure-ledger reduction of the
locked surrogate's parameters to closure-quantum invariants.

## Quark mass ladder (u, d, s, c, b, t) from a shell-coupled S³ closure

Parallel to the lepton sector, the six observed quark masses are
fit by a 6×6 Hermitian Hamiltonian on the closure basis
`{(k=1,±), (k=3,±), (k=5,±)}`.  The minimal v3 ansatz did not
suffice; three opt-in structural extensions (`uplift_mode =
"partition_asymmetric"`, `spectrum_zero_mode = "min_eigenvalue"`,
`chi_q_k3`, `eta_k3k5_minus`), all with defaults that recover
the minimal lepton-style ansatz, give the locked spectrum.

### Locked spectrum (d-anchor, max rel err 1.6%)

Anchored on `d = 4.67 MeV`; `u` is at zero by construction under
min-eigenvalue spectrum zero.

| species | predicted (MeV) | observed (MeV) | rel err |
|---------|----------------:|---------------:|--------:|
| u | 0           | 2.16    | 1.00 (by construction) |
| d | 4.67        | 4.67    | 0 (anchor)             |
| s | 94.82       | 93.4    | **1.5%**               |
| c | 1290.92     | 1270    | **1.7%**               |
| b | 4219.92     | 4180    | **0.95%**              |
| t | 170342.41   | 172690  | **1.4%**               |

### Locked parameters (constraint-reduced)

The full residual sector is *derivable from existing geometry*
on the eigensolver's tortoise grid:

| sector | reading |
|--------|---------|
| `action_base = π` | structural |
| `uplift_asymmetry ε = 1 − 1/k_5² = 24/25` | partition asymmetry from inverse-square shell scaling |
| `eta_k3k5_minus η = k_5 = 5` | (3,−)–(5,−) targeted off-diagonal coupling |
| `chi_q_k3 χ = (k_5 − 1)·k_5 = 20` | k = 3 partition splitter |
| `phase = 0` | partition-mixing channel inactive at the lock |
| `gamma_q = 1/10` | empirical clean rational |
| `transport ≈ 0.54` | mean `⟨u_l\|V_{l+2}−V_l\|u_{l+2}⟩` on tortoise grid (+0.87% off) |
| `pinhole ≈ 22.25` | `Σ_{l=1..5} V_max(l)` on tortoise grid (−1.09% off) |
| `resistance ≈ 0.14` | `transport · ln(α_q(k_5)/α_q(k_1))` (−0.43% off) |
| `β = N · π/2 with N=466` | **remaining phenomenological parameter** |

### Shell-coupled vs minimal closure

The diagonal-Hamiltonian decomposition shows what makes the
quark ladder structurally distinct from the lepton ladder:

| species | β contribution |
|---------|---------------:|
| u, d (k=1) | 0% |
| s         | +11% (level mixing only) |
| c         | **−27%** (pushed *down* by level repulsion) |
| b         | +76% via β·4·(1−ε) = β·4/k_5² |
| t         | **+99%** via β·4·(1+ε) ≈ β·4·(49/25) |

`β` only enters at the heaviest shell (k=5), via the
partition-asymmetric `(1±ε)` factor.  The lighter shells (u, d,
s, c) are determined entirely by the chamber-coupling sector
(pinhole, χ, γ_q).  This is the operational signature of the
"shell-coupled closure" picture: the same S³ closure skeleton
that drives the lepton ladder is, in the quark sector, primarily
expressed through how the closure interacts with an interior
chamber rather than through global pass-count winding.

### Calibration archaeology

| Script | Purpose |
|--------|---------|
| `scripts/calibrate_quark_ratios.py` | Coarse grid over the residual sector; identifies γ_q regime where positivity holds. |
| `scripts/sweep_quark_beta.py` | Integer-winding β sweep (now known to be a fit knob, not a topological lock). |
| `scripts/map_basin_quark_uplift.py` | Basin probe around the best β. |
| `scripts/lock_quark_beta_probe.py` | Final lock with β hard-fixed (legacy from the integer-N attempt). |
| `scripts/experiment_partition_asymmetric_uplift.py` | Tests the k=5 b/t splitter. |
| `scripts/experiment_min_eigenvalue_zero.py` | Tests d-anchor with min-eigenvalue spectrum zero. |
| `scripts/experiment_k3_splitter.py` | Tests χ for the c/s splitter. |
| `scripts/experiment_refined_k3k5.py` | Pass-2 refinement crossing the user-named "serious candidate" threshold (max rel err < 0.3 → 0.13). |
| `scripts/basin_probe_topological_locks.py` | Verifies N, χ, η are basin features, not grid coincidences. |
| `scripts/refine_pass3_coord_descent.py` | Coordinate-descent refinement to 1.6%. |
| `scripts/experiment_constraint_search.py` | Constraint-reduction pass: 9 free knobs → 4 + 1. |
| `scripts/experiment_n_ablation.py` | First N-stability check (residuals free); N drifts. |
| `scripts/experiment_residuals_from_geometry.py` | Substitutes residuals with broad geometric scalars. |
| `scripts/experiment_transport_pinhole_search.py` | 1D refinement of transport and pinhole derivations. |
| `scripts/experiment_transport_overlap.py` | Derives transport from QM perturbation overlap to within 0.87%. |
| `scripts/experiment_resistance_wkb.py` | WKB tunneling-derived resistance (negative result), then discovers `resistance = transport · ln(α_q ratio)` to within 0.43%. |
| `scripts/experiment_n_ablation_geometric.py` | Decisive N-stability check with all residuals derived; N still drifts → β is phenomenological. |

See `docs/quark_axioms.md` (full v3 spec, calibration log §8,
phenomenological interpretation §9) and the JSON archive in
`docs/calibration_runs/` for the raw outputs of every scan.

## QFT-event-reinterpretation thread (Compton scattering)

An 11-PR thread (PRs #25 – this PR) testing whether BAM's three
composable dynamical elements — **throat worldlines + time dilation
at mouth + antipodal closure** — reproduce QFT event structure for a
canonical local interaction, Compton scattering `γ + e → γ + e`. The
thread progressively identified the BAM-native ingredients needed
to reproduce Klein-Nishina, then resummed the perturbative result
into a closed-form vertex factor.

### Result chain

  - **Kinematics** (PR #25): closure-compatible. The antipodal map
    `(E, **p**) → (E, −**p**)` automatically conserves the
    back-vertex when the front does. Inter-mouth proper-time skew
    vanishes; throat-pinch skew is a recoil-induced `O(ω²/m²)`
    quantity, not a topological closure quantum.

  - **Propagator** (PR #26): the `S³` Green function
    `G(ψ) ∼ 1/(4πψ)` with `ψ = (s − m²)/(2m²)` reproduces the QED
    propagator pole `1/(s − m²)` exactly (fitted exponent 1.0002).

  - **Photon structure** (PR #28): giving the photon two transverse
    polarisations on the `S³` tangent bundle and treating the
    electron as a scalar charge in the Thomson limit reproduces
    `(1 + cos²θ)/2` exactly — the full Klein-Nishina angular factor.

  - **Finite-energy gap** (PR #29): the natural BAM construction
    fails at `O(ω/m)`. The recoil sign is qualitatively wrong
    (BAM enhances backscatter, KN suppresses it), localised to the
    missing per-channel kinematic weighting.

  - **Vertex coupling** (PRs #30, #31): an extended Family B vertex
    modification `V = (ε·ε'*)·(1 + ε·μ₁ + ...)` with
    `μ₁ = γ·(ω/m)·(1 − cos θ)` closes the `O(ε)` gap exactly at
    `γ = −3/2` — derived analytically from a 4-equation linear
    system over `{1, c, c², c³}` basis.

  - **Coefficient origin** (PRs #32, #33): 8 natural BAM ingredients
    evaluate to `−3/2`; the dimensional-scaling test in `d ∈ {3, 4,
    5, 6}` falsifies the embedding-dim / polarisation-count origin
    (candidate C), leaving 7 surviving candidates rooted in
    group-theoretic invariants of SU(2).

  - **`O(ε²)` extension** (PR #34): polynomial leading-order
    closure with `(ν₀, ν₁, ν₂, ξ) = (9/4, −4, 7/4, −1/2)`, with
    structural patterns `ν₀ = γ² = (−3/2)²` (recursive) and
    `ξ = −A_φ(0)` (Hopf-charge link).

  - **Resummation** (PR #35): the closed form

      F²(x, c) = 4·x³·(x² + 1 − x·sin²θ) / [(1 + c²)·(1 + x)²]
              = (2x/(1+x))² · [x·(x²+1−x·sin²θ) / (1+c²)]

    with `x = ω'/ω = 1/(1 + ε(1 − cos θ))` reproduces Klein-Nishina
    **exactly at all orders in ε up to ε ~ 2** (machine precision).
    The perturbative results of PRs #31–34 are Taylor expansions
    of this closed form.

  - **Cross-process validation via Breit–Wheeler** (this PR): the
    same closed-form F, expressed in Lorentz invariants and
    analytically continued via standard Mandelstam crossing
    (`s_C → u_BW`, `t_C → s_BW`, `u_C → t_BW`), exactly reproduces
    the Breit–Wheeler pair-production amplitude `γγ → e⁺e⁻`.
    Crossed variables `x_⊗ = −(1−β·cosθ)/(1+β·cosθ) < 0` and
    `c_⊗ = (2β² − β²cos²θ − 1)/(1−β²cos²θ)` carry the construction
    from Compton lab kinematics to BW CM kinematics; the
    BAM-predicted `|M̄|²_BW = −2·(f_baseline · F²)/x_⊗²` agrees
    with the textbook formula to machine precision at all sampled
    `(β, cosθ)`, and the integrated differential reproduces the
    textbook BW total at threshold (`β → 0` linear) and in the
    ultra-relativistic regime (`β → 1` logarithmic). The vertex F
    is therefore **not a Compton-specific algebraic fit** — it is
    the closed form of the invariant QED amplitude carried by
    crossing to its tree-level partners.

### Structural reading

The `(1 + c²)` denominator in the angular factor IS the
polarisation-sum factor. The closed-form F must be derived AS a
modification of the polarisation-sum projector, not as a separate
amplitude factor. The two-factor decomposition

  - kinematic Padé `(2x/(1+x))²` — pure x-function
  - angular polarisation modification `[x·(x²+1−x·sin²θ) / (1+c²)]`

suggests two BAM-native ingredients combine to produce the full
vertex coupling. The clean half-integer/integer rationals appearing
at every order (γ = −3/2, ν₀ = 9/4, ν₁ = −4, ν₂ = 7/4, ξ = −1/2)
indicate a deeper geometric origin awaiting first-principles
derivation from the Hopf-bundle / throat-transport algebra.

### What survives and what is still open

  - Survives: BAM's antipodal-`S³` propagator + Hopf-fibre photon
    polarisation + closed-form vertex `F²` together reproduce
    Klein-Nishina exactly. The same closed form, crossed via
    Mandelstam permutation, reproduces Breit–Wheeler `γγ → e⁺e⁻`
    (PR #36) and pair annihilation `e⁺e⁻ → γγ` (this PR); the full
    Compton/BW/annihilation crossing triangle closes (loop is
    identity at both the Mandelstam-label and amplitude level).
  - Open: first-principles BAM derivation of `F²` from a BAM
    Lagrangian / action. Two-channel tree processes (Bhabha, Møller)
    with interfering s+t diagrams; loop corrections requiring the
    bulk radial channel.

### Probe sequence

| # | Probe | Outcome |
|---|---|---|
| PR #25 | `compton_antipodal_kinematics_probe.py` | closure-compatible |
| PR #26 | `compton_amplitude_structure_probe.py` | propagator ✓, polarization ✗ |
| PR #28 | `compton_photon_structure_probe.py` | Thomson KN ✓ |
| PR #29 | `compton_finite_energy_kn_probe.py` | recoil ✗ at `O(ω/m)` |
| PR #30 | `compton_vertex_structure_probe.py` | empirical finite-ε fit |
| PR #31 | `compton_vertex_derivation_probe.py` | exact γ = −3/2 |
| PR #32 | `compton_coefficient_origin_probe.py` | 8 plausible derivations |
| PR #33 | `compton_dimensional_scaling_probe.py` | C falsified, 7 survive |
| PR #34 | `compton_eps2_extension_probe.py` | `O(ε²)` polynomial fit |
| PR #35 | `compton_vertex_resummation_probe.py` | exact closed-form F² |
| PR #36 | `breit_wheeler_cross_process_probe.py` | F process-general under Compton → BW crossing |
| PR #37 | `pair_annihilation_crossing_probe.py` | full Compton/BW/annihilation crossing triangle closes |
| PR #38 | `throat_nucleation_caustic_derivation_probe.py` | F² = K(x)²·Q(x, c) BAM-geometric decomposition |
| PR #39 | `two_mouth_flux_action_probe.py` | K(x) = 2x/(1+x) from equal-action throat-rate splitting |
| PR #40 | `hopf_helicity_transport_probe.py` | Q(x, c) from Hopf-fibre helicity spinor (A_pres, A_flip) |
| PR #41 | `throat_action_derivation_probe.py` | BAM throat action: both equal-action postulates derived from S³ antipodal symmetry + closure quantum + stationary action |
| PR #42 | `bhabha_moller_interference_probe.py` | 4-fermion gap identified: scalar Compton kernel insufficient for Bhabha/Møller |
| PR #43 | `dirac_trace_geometry_probe.py` | 4-fermion diagonal numerators (s²+u²), (u²+t²), (s²+t²) from SU(2) Hopf-bundle Pauli traces |
| PR #44 | `mobius_exchange_sign_probe.py` | Bhabha/Møller interference signs from T = iσ_y = ε non-orientable throat transport |
| PR #45 | `bam_exchange_kernel_probe.py` | photon propagator magnitude 1/q² from S³ Green function (flat limit) |
| PR #46 | `hopf_vector_exchange_kernel_probe.py` | **photon propagator Lorentz tensor −η^{μν}/q² from Hopf-bundle U(1) connection** |
| PR #48 | `two_throat_coulomb_probe.py` | inverse-square Coulomb force from the S³ Green response; Gauss law exact |
| PR #49 | `topological_discrete_sector_probe.py` | scaffold B1+B2 promoted to action data (RP³ + spin structure + winding θ-term) |
| PR #50 | `radial_reduction_bridge_probe.py` | scaffold B5 factorized: 5D→4D into three channels; F² not a radial overlap |
| PR #51 | `bulk_boundary_interaction_probe.py` | scaffold B5′: radial (masses) + throat (K) unified by one bulk-boundary cavity |
| PR #51 | `master_integral_probe.py` | **scaffold B5 closed: masses and F²=K²·Q from one C×S³ master functional** |
| PR #52 | `maslov_dimensional_bridge_probe.py` | scaffold B4 audit: irreducible by scale-freeness; Maslov closure-ledger (radial +1 = μ=4) |
| PR #53 | `delta_r_scale_modulus_probe.py` | scaffold B4 anchor: ΔR is a cosmologically-invariant bulk separation |

**Synthesis / release note:** `docs/tree_qed_status.md` collects the
PR #35 → #46 result — all tree-level `2 → 2` QED scalar intensities
(Compton, Breit–Wheeler, pair annihilation, Bhabha, Møller)
reproduced from BAM-geometric primitives.

The Compton derivation rests on the algebraic identity

  x² + 1 − x·sin²θ ≡ (1 − x)² + x · (1 + c²)

which yields two equivalent decompositions:

  F²(x, c) = [2x/(1+x)]² · [x² + x·(1−x)²/(1+c²)]
  |M̄|²_KN/(8e⁴) = (1+c²) + (1−x)²/x

with BAM-geometric interpretation:

  - **P(x) = 2x/(1+x)** = harmonic mean of in/out photon frequencies
    = standard classical bottleneck-flux average through the throat;
    squared because both throat-pair mouths pinch. Uniquely
    polynomial — alternative throat-rates (arithmetic, geometric mean,
    linear x) leave Q non-polynomial at x → −1.
  - **(1+c²)/2 = cos⁴(θ/2) + sin⁴(θ/2)** = sum of squared Wigner-d¹₁,±₁
    matrix elements = Hopf-fibre spin-1 helicity transport through θ.
  - **Q = |a|² + |b|²** = orthogonal sum of helicity-preserving
    (a = x) and helicity-flipping (b = √x(1−x)/√(1+c²)) channels,
    each non-negative across the physical region.
  - The Hopf connection at the BAM lock `A_φ(0) = 1/2` (from
    `geometrodynamics.hopf.connection`) matches the PR #34 perturbative
    coefficient `ξ = −1/2` exactly.
  - Decomposition survives analytic continuation under crossing
    (full Compton ↔ BW ↔ annihilation triangle, PR #37).

The full F² closed form is derived from three foundational
principles via a single BAM throat action functional (PR #41):

  (P1) closure quantum `S = 2π` (BAM `action_base`)
  (P2) S³ antipodal symmetry `σ(p) = −p` (involution swapping mouths)
  (P3) stationary action under the antipodally-symmetric ansatz

Both equal-action postulates (PR #39 energy → K, PR #40 spin/Hopf → Q)
follow as consequences. Alternative principles (broken antipodal
symmetry; wrong closure quantum; wrong action functional) all fail
to reproduce K(x), confirming the principles are necessary.

The thread then extends to 4-fermion tree QED (Bhabha, Møller,
PRs #42–#46): SU(2) Hopf-bundle Pauli traces give the Dirac-trace
diagonal numerators (#43), the non-orientable throat transport
`T = iσ_y = ε` gives the Fermi-statistics interference signs (#44),
and the `S³` Green function (scalar #45, Hopf-bundle vector #46)
gives the photon propagator `1/q²` with full Lorentz tensor
structure. End-to-end Bhabha and Møller `|M̄|²` match QED to machine
precision from BAM-geometric ingredients alone.

See `docs/tree_qed_status.md` for the full synthesis. Per-PR
research plans: `docs/compton_vertex_resummation_research_plan.md`
(#35), `docs/breit_wheeler_cross_process_research_plan.md` (#36),
`docs/pair_annihilation_crossing_research_plan.md` (#37),
`docs/throat_nucleation_caustic_derivation_research_plan.md` (#38),
`docs/two_mouth_flux_action_research_plan.md` (#39),
`docs/hopf_helicity_transport_research_plan.md` (#40),
`docs/throat_action_derivation_research_plan.md` (#41),
`docs/bhabha_moller_interference_research_plan.md` (#42),
`docs/dirac_trace_geometry_research_plan.md` (#43),
`docs/mobius_exchange_sign_research_plan.md` (#44),
`docs/bam_exchange_kernel_research_plan.md` (#45), and
`docs/hopf_vector_exchange_kernel_research_plan.md` (#46).

### BAM effective-action scaffold — barrier closure (PRs #49–#53)

The tree-QED ingredients above were assembled into a single covariant
5D effective-action scaffold and its five mismatch terms (B1–B5) were
worked off one by one. Four are now **closed**:

| barrier | what it was | now |
|---|---|---|
| **B1** closure quantum `∮A = 2πn` | imposed constraint | winding θ-term `S_top = 2π·n` |
| **B2** antipodal `Z₂` (`T = iσ_y`) | imposed identification | `RP³ = S³/Z₂` + non-trivial spin structure |
| **B3** hard-wall throat BC | imposed by hand | single-valuedness under `T² = −I` ⟹ `ψ(throat) = 0` |
| **B5** 5D→4D reduction producing F² | unconstructed | one master functional yields masses **and** `F²=K²·Q` |

B5 is closed by the **master integral**: a single separable functional
on the warped-product internal geometry `M_int = C × S³`
(`C` = radial cavity `[R_MID, R_OUTER]`),

```
ℳ(ω; x, c) = G_C(r, r′; ω) ⊗ 𝒢_{S³}(Ω, Ω′)
```

read three ways from one object —

  - **poles in ω** → the mass spectrum `ω(l,n)` (radial ladder `n` ×
    S³ Casimir `l`, the centrifugal term of the warp);
  - **throat boundary of `G_C`** → `K(x) = 2x/(1+x)` (dwell-time
    impedance `Z(ω)=π/ω` in series);
  - **S³ Hopf reduction of `𝒢_{S³}`** → `Q(x,c) = x²+x(1−x)²/(1+c²)`
    (Hopf-fibre helicity spinor).

The vertex residue reproduces `F²(x,c) = K²·Q` to machine precision
(`2e-14`) while the poles give the masses — **masses and the full
vertex from one functional**. The `F²=K²·Q` factorization is the direct
consequence of the product internal geometry (separation of variables),
not a failure to unify.

The fifth barrier **B4** (the dimensional bridge `ℏ = m_e·R_MID·c`) is
not a gap but a **structural necessity**: the closure-ledger/Maslov
machinery is *scale-free* (rescaling `R_MID → λ·R_MID` leaves every
dimensionless output invariant), so exactly one external dimensionful
anchor is mathematically required — **B4 is irreducible** (#52). That
anchor need not be a particle mass: it can be the **invariant bulk
separation** `ΔR = R_OUTER − R_INNER`, a proper (cosmologically fixed)
length, giving `m_e = f_closure·ℏ/(ΔR·c)` with `f_closure = 0.52` (#53).
The scaffold is therefore complete: four barriers derived, the fifth the
single mandatory dimensionful unit. Full ledger:
`docs/bam_scaffold_status.md`; closure release note:
`docs/scaffold_closure_release_note.md`; per-probe plans:
`docs/topological_discrete_sector_research_plan.md` (#49),
`docs/radial_reduction_bridge_research_plan.md` (#50),
`docs/bulk_boundary_interaction_research_plan.md` and
`docs/master_integral_research_plan.md` (#51),
`docs/maslov_dimensional_bridge_research_plan.md` (#52),
`docs/delta_r_scale_modulus_research_plan.md` (#53).

### Throat-as-particle arc (PRs #55–#74)

With the scaffold closed, the same primitives extend through the
lepton/QCD sector arc:

| arc | PRs | summary |
|---|---|---|
| **Throat as anchor** | #55–#58 | `R_MID` recast as finite-self-energy equilibrium (#55), cohesive `B·R²` = brane tension (#56), bulk-gravity tuning factor √6 (#57), pair-threshold `2 m_e c²` (#58). |
| **Throat = relativistic spin-½ particle** | #59–#62 | Moving throat dispersion `E²−(pc)²=(mc²)²` (#59), Hopf-holonomy Wigner rotation (#60), `g = 2` from Pauli/SU(2) + Hopf monopole (#61), one-loop `a = α/2π` reconstructed (#62). |
| **C, CPT, throat Dirac spinor** | #63–#66 | `C` = inner/outer swap `c₁ → −c₁` (#63), CPT on throat histories (#64), explicit `Θ = γ⁰γ¹γ²γ³ = −iγ⁵` on throat spinor (#65), throat 4-spinor from `S_BAM` SUSY factorization (#66). |
| **Even-k absence → QCD shell** | #67–#69 | Even-`k` absence = spin-statistics selection rule (#67), higher excitations transition into QCD shell channel (#68), shell ↔ QCD structural match (#69). |
| **Three generations / `k_5 = 5`** | #70–#74 | Sharp `k ≤ 5` boundary (#70), `β_lepton = k_5²·(2π) = 50π` (#71), `#generations = (k_5+1)/2 = 3` (#72), `k_5 = D_bulk = dim(S³)+2 = 5` (#73), `1/(2π)` in Schwinger anomaly = BAM closure-quantum loop measure (#74). |

### QCD-shell arc (PRs #76–#80) — quarks as cavity wavefronts

The quark sector is reframed via the user's physical insight:
**"Quarks do not pass through the throat; they are the wavefronts
that resolve the cavity itself."** This is the quantitative
development of PRs #68–#69 (throat-to-shell transition + shell ↔ QCD
structural match) that PR #76 identified as the right derivation
route.

| arc | PRs | summary |
|---|---|---|
| **`n_part = 233` diagnosis** | #76 | Extended candidate catalog (Fibonacci, color × flavor × generation, QCD β₀, Tangherlini QCD-shell modes); no enumeration survives §8 drift. v3 Hamiltonian is **lepton-shaped** — wrong machinery for the quark sector. Right derivation route is the QCD shell waveguide. |
| **Shell waveguide scaffold** | #77 | 6-state `(l, n, p)` basis: `l` = S³ Casimir, `n` = shell-saturated radial overtone (≥ 3 for l=1), `p ∈ {+, −}` = Z₂ partition. Operator scaffold `H = H_kin + H_Z2 + H_couple` with `H_kin = ω²(l, n)` cavity-eigenfrequency-squared (NOT lepton `β·k²·(2π)`). 3 × 2 = 6 flavors matches PR #69. |
| **Mass-ordering audit** | #78 | Shell basis structurally better than v3 in 4 ways. Uniform `χ·σ_z` cannot reproduce within-generation inversion (best 2/3 blocks); sign-flipping `χ_n` can (existence proof). Coverage gap: shell kinetic ×2.2 vs observed ×6.4·10⁹; `n_part` not resolved at #78 alone. |
| **Boundary-stress `χ_n`** | #79 | `χ_n = T_odd(n) = (T_inner − T_outer)/2` from Z₂-antisymmetric piece of cavity-mouth stress (PR #63's inner/outer swap). NO free parameter. Uniform-positive sign (no flip), shell-suppressed — 30–100× too small for observed splittings. PR #78 sign-flipping ansatz **overruled** by the structural derivation. |
| **Color algebra** | #80 | **BAM-native color algebra = SU(2) × Z₂** (SU(2) from B2 / Hopf, Z₂ from PR #63). SU(3) NOT derivable from current scaffold (all natural triplets give SO(3)/SU(2)). Pati-Salam SU(4) requires throat↔shell algebra map (open extension). v3 species map revised: `+ = heavier` uniformly. Inter-generation mass hierarchy (~9 orders in mass²) is **outside the scope** of any BAM color algebra on the shell basis. |

**Arc closure summary.** The four-PR arc (#77 → #80) closes
structurally — the shell basis is the right machinery, `χ_n` is
derived without a free parameter, the BAM-native color algebra is
identified, and the v3 species map is settled. What remained open at
#80: the inter-generation mass hierarchy and the Pati-Salam SU(4)
extension.

### Pati-Salam bridge + mass-operator unification (PRs #82–#83)

| arc | PRs | summary |
|---|---|---|
| **Throat ↔ shell `n+3` bridge** | #82 | Each generation has a lepton at `n = g−1` (throat) and a quark-pair at `n = g+2` (shell); shift `+3` = PR #68 shell threshold (no free parameter). Unified 12-state `(l, n, p)` basis + throat-shell Z₂. Full SU(4) PS needs 3 open extensions: BAM-native neutrinos, 3-fold quark color, **lepton-quark mass-operator unification**. |
| **Bohr-Sommerfeld mass-operator unification** | #83 | The third extension is **closed at the structural-form level**: the lepton `β·k²` (PR #71) and quark `ω²(l,n)` (PR #77) mass operators are the SAME Bohr-Sommerfeld operator `m²(k,n) = (k·2π/L_throat)² + ((n+1)·π/L_cavity)²`, `L_throat = √(2π)/k_5`. Cavity `∮√(ω²−V)dr* = (n+1)·π` verified to machine precision; `(2π/L_throat)² = k_5²·(2π) = 50π = β_lepton` recovered. |

**The unification, in one line.** Leptons and quarks are not two kinds
of object with two mass formulas. They are **one Bohr-Sommerfeld
closure operator** `m² = (S/L_eff)²` read in two channels of the
closure ledger (PR #52's `N_total = N_layer1 + N_radial`):

  - **Leptons wind through the throat** — winding number `k ∈ {1,3,5}`,
    closure quantum `2π` (full S³ great circle) → `m² ≈ β·k²`.
  - **Quarks resolve the cavity** — `k = 0` (no throat traversal),
    radial overtone `n ∈ {3,4,5}`, closure quantum `π` (half-cycle
    Bohr-Sommerfeld node) → `m² ≈ ω²(l, n)`.

The user's physical insight — *"quarks do not pass through the throat;
they are the wavefronts that resolve the cavity itself"* — is exactly
`k = 0` in this single operator. The `2π`-vs-`π` distinction between the
two channels is BAM's pervasive full/half-cycle structure (throat dwell
`τ = π/ω`, Hopf holonomy `∮A = π cos χ`, B3 reflection phase `π`).
What remains open: an independent derivation of the two `L_eff` from one
principle, and the inter-generation hierarchy (the cross-channel /
mixed-mode question).

### Neutrino & full-quadrant sector (PRs #85–#87)

With the lepton/quark mass operator unified (PR #83), the `(k, n)`
plane splits into four quadrants, and the chargeless `k = 0` corner
turns out to be the neutrino — the long-open "BAM-native neutrino"
extension of the Pati-Salam bridge (PR #82).

| arc | PRs | summary |
|---|---|---|
| **Four-quadrant map / leptoquark** | #85 | The unified `(k, n)` operator's fourth quadrant (winding **and** shell-saturated, `k≠0, n≥3`) is the **leptoquark** sector, completing the reading: lepton `(k≠0, n<3)`, quark `(k=0, n≥3)`, neutrino `(k=0, n<3)`, leptoquark `(k≠0, n≥3)`. |
| **Neutrino = Majorana seesaw** | #86 | The `(k=0, n<3)` quadrant gives the lightest states but ~10⁵–10⁶ too heavy. The fix is BAM-native: `k=0 ⟹ c₁=0 ⟹ C-invariant` (PR #63) ⟹ the neutrino is its own antiparticle ⟹ **Majorana**, so it admits the seesaw `m_ν = m_D²/M_R`. The seesaw is available **only** to the chargeless sector — charged leptons (`c₁=±1`) are Dirac and keep their full winding mass — which is precisely why only neutrinos are anomalously light. Required `M_R ≈ 0.3–1.8 TeV` was left open (no BAM scale at ~TeV). |
| **`M_R` from throat-nucleation tunnelling** | #87 | The `ΔL=2` Majorana coupling **is** the PR #58 throat↔antithroat (antipodal `Z₂`) transition, and PR #58's `Σc₁=0` applied to a single state **is** PR #86's only-neutrino selection rule (`k=0` flips `0→0`, allowed; `k≠0` gives `Σc₁=∓2`, forbidden). The literal `M_R = `barrier-height hypothesis is **falsified** — with the electron-throat `σ, ρ` the barrier is `E_c ≈ 2.8 keV`, ~10⁸ too small. Instead the suppression is **tunnelling through** the barrier, `m_ν = m_D·e^{−S}`, so `M_R = m_D²/m_ν = m_D·e^{S}`: the ~TeV scale is the keV Dirac floor exponentially lifted, and the open input is recast from a mysterious ~TeV mass to a modest, generation-stable bounce action `S ≈ 15–18` — exactly the instanton-rate follow-on PR #58 flagged. |
| **Bounce action `S` = non-orientable tortoise log** | #88 | A reduced Euclidean bounce `S = √(2 μ E_c)·L*(ε)` for the flip, run along the odd (`c₁→−c₁`) tortoise path. The 5D tortoise coordinate diverges logarithmically at the throat, giving two structural results: a **rigid throat ⟹ exactly massless neutrino** (the boundary compliance `ε` is the mass-generating parameter, and the smallness of `m_ν` is the near-rigidity of the throat), and `S ∝ ln(1/ε)` — naturally `O(10)` and generation-stable, the form PR #87 required. **Honest magnitude:** the EM-throat tension under-produces `S` by ~40× (`S ≲ 1` even near-rigid); matching `S ≈ 15–18` needs a `ΔL=2` (B−L) throat tension `~6–12×` stiffer. The open input is localised once more: ~TeV mass (#86) → `O(15)` action `S` (#87) → `O(10)` B−L/EM tension ratio (#88). |
| **B−L tension ratio bracketed by closure & winding** | #89 | The `ΔL=2` flip reverses orientation (`c₁→−c₁`) — a **global** operation — so `t` is a global-closure enhancement of the **local** EM surface tension. It is bracketed, parameter-free, by the two basic BAM action scales: the **closure quantum `2π`** (a single great-circle orientation reversal, lower) and the **winding action `k_5√(2π) = √β_lepton`** (a full throat winding, upper), giving `t ∈ [2π, k_5√(2π)] ≈ [6.28, 12.53]` — **exactly** PR #88's required `6–12` (the computed `[6.41, 12.05]` sits inside). So the `6–12` band was not a fit but the BAM closure-to-winding window. The residual is "where in the window" = the compliance `ε` (`t=2π ⟹ ε≈6e-7`, `t=√β ⟹ ε≈1.3e-2`); the winding/cavity mass ratio `m_charged/m_D ≈ 11.9 ≈ √β` corroborates the winding edge. |
| **Compliance `ε` from bulk geometry → `m_ν` scale** | #90 | The capstone. `ε` is the chargeless throat's sub-throat **healing length** (`ε = ℓ²/2rs` from the neck warp `f≈2(r−rs)/rs`), sub-throat *for the neutrino* because the `c₁=0` neck has no EM term to prop it open (the charged `c₁=±1` neck is propped open and stays Dirac) — the same chargelessness that makes the neutrino Majorana makes its `ε` tiny, hence its mass tiny. Natural BAM sub-throat scales (`R_c³, Δ³, (m_D/m_ch)²`) land `ε` in the PR #89 window; with the winding-edge tension `t≈√β` (cross-check-favoured) the chain gives `S ≈ 15–19` and **`m_ν ~ few meV`** — the observed scale, with no input outside the throat geometry. At the `2π` edge `S≈4` (too small): the chain closes only at the winding edge. |
| **Generation spread + `PMNS ≫ CKM`** | #91 | Generations are the cavity radial overtones `n`, so the bare prediction is **normal ordering** with `m_ν ∝ m_D` (cavity-floor ratios `1:1.87:2.74`). The spread is widened in the right direction by the overtone-dependent neck coupling — PR #79's boundary stress `χ_n` decreases with `n` (0.304, 0.097, 0.039), so higher-`n` neutrinos are less throat-coupled ⟹ more compliant ⟹ less suppressed ⟹ relatively heavier (lifting `m₃` toward the observed spread). **Headline:** large PMNS vs small CKM is the BAM **cross-channel** (leptons: charged throat-winding `k≠0` × neutrino cavity-resolving `k=0`) vs **intra-channel** (quarks: up & down both cavity-shell `k=0`) distinction — the structural reason `PMNS ≫ CKM`. Precise spectrum (`ε_n(χ_n)` `O(1)`, absolute scale unmeasured) and explicit angles open. |
| **PMNS anarchic vs CKM aligned (quantitative)** | #92 | Computes the cross-channel overlap. A naive radial overlap gives near-permutation (small) mixing — large PMNS is *not* a literal mode overlap. The lepton generation labels live in **different coordinates** (charged: closure-winding `k`; neutrino: radial-overtone `n`), so the map has no preferred alignment ⟹ **anarchic** (Haar-random). Observed PMNS (33.4°, 49°, 8.6°) is **typical** of a Haar `U(3)` (30th/57th/4th percentile); CKM (13°, 2.4°, 0.2°) is **extremely atypical** (joint `p ≈ 0`) = aligned — up & down share the radial-overtone (shell) coordinate. So PMNS ∈ anarchy class (cross-coordinate), CKM ∈ aligned class (intra-coordinate) — a falsifiable separation matching observation. Specific angles not pinned (anarchy is statistical; θ13 at the 4th percentile is the mild tension). |
| **θ13 suppression / residual alignment** | #93 | Resolves the PR #92 θ13 tension. θ13 = `U_e3` is the corner element — it links the lowest winding (`k=1`) to the highest overtone (`n=2`), the most coordinate-distant (**two-hop**) pair — while θ12, θ23 are adjacent (one hop). Since the throat↔shell coupling is **local** in the `(k,n)` lattice (PR #82 `+3` shift, PR #83 operator), the corner is a suppressed two-hop amplitude — a residual **nearest-neighbour** alignment. A structured-anarchy model (corner variance `exp(−μ)`, `μ=0` = pure anarchy) with `μ≈3` shifts the θ13 distribution down (median 33°→~16°), makes θ13 robustly the *smallest* angle (frac 0.50→0.72), and moves observed θ13=8.6° from the 4th to the ~21st percentile — **tension resolved** — while θ12 (~44th) and θ23 (~70th) stay typical. The exact value (μ; θ13 median saturates ~14–16°) is open. |
| **CP / Majorana phases** | #94 | The phase sector. **CP violation is generic**: the winding amplitudes carry the Hopf holonomy `e^{ikχ}` (PR #60), so the cross-channel overlaps are intrinsically complex and `δ_CP ≠ 0, π` with probability 1 (CP conservation is measure-zero). The **Jarlskog invariant** mirrors the angle dichotomy: `|J_PMNS| ≈ 0.026` is typical of anarchy (51st/81st percentile → large CP violation), `|J_CKM| ≈ 3×10⁻⁵` is extremely atypical (~0.1th → aligned/suppressed). And the **two Majorana phases exist** because the neutrino is Majorana (`c₁=0`, PR #86) — CP phases of the ΔL=2 throat↔antithroat sector (PRs #87–#90), observable in 0νββ; a Dirac neutrino would have none. The specific phase values are anarchic (uniform) — not pinned (`δ_CP` is itself poorly measured, consistent with uniform). |
| **0νββ effective mass** | #95 | Turns the whole arc into one falsifiable number-range. The effective Majorana mass `m_ββ = |Σ U_ei² m_i|` combines: 0νββ **occurs** (neutrino Majorana ⟸ `c₁=0`, PR #86; Dirac would forbid it); **normal ordering** (PR #91) selects the NO band (`m_ββ ≈ 1.5–3.7 meV` at zero lightest mass); **anarchic Majorana phases** (PR #94) populate the full band incl. a cancellation trough (`m_ββ → ~0` around `m_lightest ~ 3–5 meV`); and the **light scale** (PR #90, ~few meV) gives `m_ββ ≲ 8 meV`. This sits below the current bound (KamLAND-Zen 28–122 meV — null result expected), largely below next-gen reach (LEGEND-1000 / nEXO ~9–20 meV), and below the inverted-ordering floor (~19 meV). **Sharp falsifier:** a 0νββ discovery with `m_ββ ≳ 19 meV` would imply inverted ordering or a quasi-degenerate scale, contradicting the BAM normal-ordering + light-scale prediction. |
| **Cosmological Σm_ν** | #96 | The cosmological companion to #95: the same light, normal-ordered spectrum fixes `Σm_ν = m1+m2+m3`. The NO floor is `√Δm²_21 + √Δm²_31 ≈ 58.7 meV` (the IO floor ≈ 99 meV), and the light scale (PR #90, ~few meV) keeps the sum pinned near it: **`Σm_ν ≈ 59–65 meV`**. This is consistent with Planck 2018 + BAO (<120 meV), just inside DESI DR1 + CMB (<72 meV), and **right at the DESI DR2 + CMB frontier (~60–64 meV)** — exactly where current cosmology is probing. **Falsifiers:** a robust `Σm_ν < 58.7 meV` excludes normal ordering (and is in tension with the oscillation `Δm²` themselves); a quasi-degenerate `Σm_ν ≳ 100 meV` contradicts the light scale. `Σm_ν` and `m_ββ` (#95) are one spectrum's two observables — a joint, cross-checkable prediction. |

**Where it lands.** The neutrino sector is now structurally complete:
the only-neutrino-Majorana selection rule, the seesaw mechanism, and a
BAM-native home for the seesaw scale (the throat↔antithroat nucleation
tunnelling) are all in place. The headline is the reframing: **`M_R` is
no longer a free ~TeV mass but an instanton action**. Because
`M_R = m_D·e^{S}`, the entire 6-order gap between the keV Dirac floor
and the TeV seesaw scale is carried by a single dimensionless number
`S ≈ 15–18` — the Euclidean bounce/instanton action for the `ΔL=2`
throat↔antithroat tunnelling. PR #88 then builds that bounce explicitly
and shows it is the **non-orientable tortoise logarithm**: a rigid
throat gives an exactly massless neutrino, and `S ∝ ln(1/ε)` is
naturally `O(10)` and generation-stable. PR #89 then constrains the
remaining tension ratio: because the flip reverses orientation it is a
*global* operation, so `t` is bracketed parameter-free by the **closure
quantum `2π`** and the **winding action `k_5√(2π) = √β_lepton`** —
`t ∈ [6.28, 12.53]`, exactly PR #88's required `6–12`. The open input
has now been localised four times — ~TeV mass (#86) → `O(15)` action
(#87) → `O(10)` tension ratio (#88) → the BAM closure-to-winding window
(#89) — leaving a single residual number: *where in that window*, i.e.
the boundary compliance `ε`. PR #90 closes the chain: `ε` is the
chargeless throat's sub-throat **healing length** (`ε = ℓ²/2rs`), tiny
*for the neutrino* precisely because its `c₁=0` neck is not propped open
by charge — the same chargelessness that makes it Majorana. With the
winding-edge tension the natural bulk scales give `S ≈ 15–19` and
**`m_ν ~ few meV`**, the observed scale, with no input outside the
throat geometry. So the whole chain — `~TeV` mass → `O(15)` action →
`O(10)` tension ratio → closure-to-winding window → sub-throat healing
length → `meV` — is closed at order-of-magnitude: **the neutrino mass
scale is geometric, not tuned.** PR #91 then takes up the spread and the
mixing: generations are the cavity overtones, so the bare prediction is
normal ordering with `m_ν ∝ m_D`, widened in the right direction by the
overtone-dependent neck coupling (PR #79's `χ_n` falls with `n`, so
higher-`n` neutrinos are less suppressed, hence heavier). And the
long-standing `PMNS ≫ CKM` puzzle is the **cross-channel vs
intra-channel** distinction: leptons mix across the throat-winding
(`k≠0`) / cavity-resolving (`k=0`) divide — large; quarks mix within the
single cavity-shell channel — small. What remains open is the precise
neutrino spectrum (an `O(1)` coefficient; the absolute scale is
unmeasured — only `Δm²`) and the explicit mixing angles. PR #92 takes up
the angles and finds the cross-channel structure is **anarchic**: because
the charged-lepton generation lives in the closure-winding coordinate and
the neutrino generation in the radial-overtone coordinate — different,
unaligned coordinates — the PMNS matrix is effectively Haar-random, and
the observed angles (33.4°, 49°, 8.6°) are *typical* of that anarchic
distribution, while CKM is *extremely atypical* (aligned, joint
`p ≈ 0`), as expected for up/down quarks sharing the one shell
coordinate. The class-level separation (PMNS anarchic, CKM aligned) is a
firm BAM prediction; the specific angles, being statistical, are not
pinned (θ13 sitting on the small side is the one mild tension). PR #93
resolves that last tension: θ13 = `U_e3` is the corner element — the
lowest winding (`k=1`) × highest overtone (`n=2`), the most
coordinate-distant pair — so it is reached by *two* channel-hops, and a
residual nearest-neighbour alignment (the throat↔shell coupling is local
in the `(k,n)` lattice) suppresses that two-hop amplitude. This makes θ13
robustly the smallest angle and moves the observed value from the 4th to
the ~21st percentile, with θ12, θ23 staying typical — leaving only the
exact value (one residual-alignment parameter) and the CP/Majorana phases.
PR #94 closes that last item structurally: CP violation is **generic**
(the winding amplitudes carry the complex Hopf holonomy `e^{ikχ}`, so the
PMNS is generically complex and CP conservation is measure-zero), the
Jarlskog invariant mirrors the angle dichotomy (`|J_PMNS|` typical of
anarchy, `|J_CKM|` extremely atypical/aligned), and the **two Majorana
phases exist** because the neutrino is Majorana (`c₁=0`, PR #86) —
observable in 0νββ, with none for a Dirac neutrino. The phase *values*,
like the angles beyond the dichotomy, are anarchic and not pinned. With
this, the neutrino arc (#85–#94) closes: the sector's *structure* —
Majorana nature, mass scale, ordering, mixing class, θ13 hierarchy, CP
genericity, Majorana-phase existence — is BAM-native, while the precise
spectrum and the specific phases/angles remain the (statistical /
one-parameter) residuals. PR #95 then collapses that structure into a
single falsifiable observable, the 0νββ effective Majorana mass: 0νββ
*occurs* (Majorana), in *normal ordering*, with *anarchic phases* and a
*light scale*, giving `m_ββ ≲ 8 meV` (with cancellation to ~0) — below
current bounds and the inverted-ordering floor (~19 meV), so a discovery
at `m_ββ ≳ 19 meV` would falsify the prediction. The neutrino sector thus
ends not just structurally complete but with a concrete experimental
target for next-generation tonne-scale 0νββ searches. PR #96 adds the
cosmological companion from the *same* spectrum: `Σm_ν ≈ 59–65 meV`,
pinned at the normal-ordering floor and sitting right at the DESI DR2 +
CMB frontier (~60–64 meV) — so the two flagship neutrino observables,
`m_ββ` (≲ 8 meV) and `Σm_ν` (~60 meV), are a joint, cross-checkable pair
that current and near-term experiments are now testing.

## Watch the geometry route a wave (PR #242)

Most of the recent arc reached the throat through algebra. `viz/throat_wavefront.py`
goes back the other way and asks what a classical wave does when the geometry
alone is allowed to act on it — and every stage of it is watchable.

**Three surfaces on one clock**, because the comparison needs all three terms:
the **bare S²** (the canonical picture: point → expanding ring → great circle →
contracting ring → antipodal focus), the same sphere with both caps cut out and
**sealed by a mirror** (what merely cutting holes does), and the same cut sphere
with the mouths **joined by a neck** (what opening a second route does).

The neck is a genuine **catenoid** — the minimal surface of revolution, `H = 0`,
`r = b cosh(z/b)`. Matching circumference *and* arclength slope at the mouth
fixes it from the mouth angle alone, `b = sin²a`, `L = sin 2a`, and its
curvature genuinely **varies**: exactly `−1` at each mouth (the sphere's `+1`
with its sign flipped) deepening to `−1/sin⁴a` at the waist. Each mouth is one
finite-volume face shared by a sphere cell and a neck cell, so the coupled solve
conserves its discrete energy to round-off (`~1e-15`).

With no fitted parameter anywhere, the wave reports the handle:

| finding | measured |
|---|---|
| the bare front sweeps each point once — a pulse cannot meet itself | max arrival count **exactly 1**, no second front anywhere |
| a sealed mouth sends a front back home; an open one does not | source-side second fronts: **12.7%** sealed vs **0%** open |
| the open/sealed echo delay **is** the neck length | `1.0024` vs `L = 0.9975` (0.49%) |
| of the energy that reaches the mouth, most crosses | transmission **91.4%** at `a = 0.75`, by integrated flux |
| a gluing twist aims where the bulk energy lands | antipodal precursor `0.477` ahead of the geodesic focus, `6.3×` the untwisted throat |
| torus vs Klein bottle is hidden at `τ ∈ {0, π}` | difference `0.0000` there, `~0.3` elsewhere |

Two things worth stating plainly. `∫K dA = −2π[r']` for any surface of
revolution, so **`χ = 0` tests the `C¹` join and not the profile** — it is not
evidence for the catenoid. And **mirror suppression is an amplitude ratio at one
point in time while transmission is an energy ratio at the mouth**; they are two
views of one fact and are not interchangeable.

The last row is the sharpened inner/outer statement: intrinsic spherical
symmetry and extrinsic orientation asymmetry coexist, and it takes a twist
that breaks the source's meridian mirror to make the orientation observable.

```bash
python -m experiments.closure_ledger.geometric_wave_refocusing_probe
# Verdict: GEOMETRY_ALONE_ROUTES_THE_WAVE  (8/8)
```

```python
from geometrodynamics.viz import (
    BareSphereSim, ThroatWaveSim, plot_wavefront_panel, run_throat_animation)

anim = run_throat_animation(ThroatWaveSim(mode="throat", mouth_angle=0.75,
                                          twist_steps=96))

sim = ThroatWaveSim(mode="throat", mouth_angle=0.75, twist_steps=96)
sim.advance_to(1.4)
plot_wavefront_panel(sim)   # both hemispheres, the neck interior, the map
```

Full write-up: `docs/geometric_wave_refocusing_research_plan.md`.

## What a wavefront has to be like to fold (PR #243)

The companion to the wave study above asks the *prior* question — not what a
wave does on a surface that already has a throat, but **what kind of wavefront
can fold at all**. It is answered with focal geometry; the closed forms are
exact, and the topology is measured *independently of them*, from the front's
own area element.

Two concentric spheres with a **flat** bulk between them. Along a straight ray
`b = r sin α` is conserved, so a ray launched inward reaches radius `b` and no
deeper.

- **A point source does not fold *here*.** Its front is the metric sphere, whose
  signed area element `t² sin θ` never changes sign — it crosses at `t = ΔR` and
  fills its own void. This is a fact about the *flat bulk*, not about point
  sources: on `S²`/`S³` the same source focuses on the antipode at `t = πR`,
  which the wave study measures directly.
- **A circle folds *coherently*.** Any curved extended source has a focal set,
  so a ring is not special for folding — it is special for folding all at once.
  Its tube's area element `t(ρ + t cos v)` first changes sign at `t = ρ`, where
  the *whole ring* arrives at the centre together (equidistant to `4.4e-16`);
  for `t > ρ` it stays singular at two axis points separating as `√(t² − ρ²)`.

Detected from the area element alone — no radius of curvature consulted — the
fold lands at `1.019806` against `ρ = 1.019804`, and reproduces `ρ` on four
unrelated rings to ~`1e-6`.

**The core result.** The ring whose caustic lands on the inner sphere
(`cos θ₀ = R_in/R_out`) launches at exactly the grazing angle
`sin α = R_in/R_out`, tangent to that sphere — both errors identically `0`. The
ring that focuses on the throat and the ray that grazes it are the same ring,
and it forms at `t = √(R_out² − R_in²)`.

**Acceptance is asymmetric; propagation is not.** Outer → inner accepts **19.1%**
of the hemisphere (Monte-Carlo 19.3%), inner → outer **100%** — a **5.2×**
difference in the *measure of launch directions that connect*. Every individual
ray is reversible, and the probe asserts it. No symmetry is broken; the ordering
of the two radii is the whole of it.

```bash
python -m experiments.closure_ledger.ring_caustic_defect_probe
# Verdict: WHOLE_RING_FOCUSES_AT_ONE_POINT  (8/8)

python scripts/geometrodynamics_v40_ring_defect.py --still ring_defect.png
```

Scope: a flat bulk (a curved one is accepted only where `r/√f(r)` is monotone,
and refused otherwise — a photon sphere breaks the argument). Nothing here is
dynamical; no throat is shown forming, which would need backreaction. A
wavefront caustic is not yet a topological defect of the geometry. Full
write-up: `docs/ring_caustic_defect_research_plan.md`.

## The geometry itself, restored (one continuous S²)

The archive rendered BAM as **one continuous surface** whose radius carried the
field, nested between two fixed shells like Russian dolls
(`archive/geometrodynamics_v39.py`). #242 and #243 replaced that with maps,
strips and meridional sections — correct, measurable, and no longer the thing
you could look at. `viz/warped_sphere.py` puts the object back, and this time
the warp is **solved**:

```
r(θ, φ, t) = R_mid + Δ · tanh( g · u(θ, φ, t) / u_ref )
R_inner = 0.74   R_mid = 1.00   R_outer = 1.26
```

The radii are deliberately the vacuole of `radial_caustic.py`, so the shell that
module's ring caustic lands on is the same inner doll drawn here.

- **One closed surface.** It carries its own poles and its `φ` seam matches to
  `2.4e-16` — nothing cut out of it. That is what makes *"a pulse sweeps every
  point once and fills its own void"* a statement about a closed manifold, and
  therefore why a **ring** is what a throat needs.
- **Nested, never touching.** Over a full return the radius stays in
  `[0.7594, 1.2226]`, clearing the inner doll by `0.0194` and the outer by
  `0.0374`. `tanh` cannot leave the gap.
- **The focus is measured.** The deepest deformation sits at geodesic distance
  `3.141593` from the source — the antipode, to `0.0e+00` — at `t = 2.9814`
  against `π`. The shortfall is the pulse's own width and shrinks monotonically
  as the pulse narrows. v39 grew its mound on the frame number instead.
- **And it rings.** The arrival drives the surface `85.6%` of the way to the
  outer shell, then inverts and pulls it `79.6%` of the way to the **inner**
  one — toward the shell the ring caustic lands on.

```bash
python -m experiments.closure_ledger.warped_sphere_geometry_probe
# Verdict: THE_SURFACE_ITSELF_CARRIES_THE_WAVE  (6/6)

python scripts/geometrodynamics_v41_warped_sphere.py --still sheet.png
```

Scope: a *display* of a solved field as a radial displacement of a **fixed**
surface — not backreaction, so no throat forms here. Sign and amplitude
ordering survive the display; ratios do not, and the probe says so. Full
write-up: `docs/warped_sphere_restoration.md`.

## Spin 0 against spin 2 on the same S²

The wave in the picture above is a **scalar**, displayed extrinsically as a
radial height. A metric perturbation is not that kind of object: `h_ab` is
symmetric and trace-free — spin 2 — and it does not push a surface outward at
all. It *shears* it. `viz/spin2_tidal.py` carries the tensor case and
`scripts/geometrodynamics_v43_tidal_sphere.py` runs both on one clock.

| | scalar `u` | tensor `h_ab` |
|---|---|---|
| displayed as | radial height | tidal ellipses |
| local effect | area changes — it breathes | area preserved to `O(h²)` |
| at the focus | peaks **on** the antipode | a **ring** around it |
| multipoles | all `ℓ ≥ 0` | **`ℓ ≥ 2` only** |
| smallest source | a point | a **ring** |

An axisymmetric spin-`s` field on the unit sphere obeys
`∂²_t h = ∂²_d h + cot d ∂_d h − (s²/sin²d) h`, eigenvalue `−ℓ(ℓ+1)` on
`ₛY_ℓ0` — the same dispersion as the scalar. Substituting `h = sin²d·q` removes
the centrifugal term exactly, and the result is integrated in conservative form
so the poles carry no flux. Three exact modes check it: `q = 1` is `ℓ = 2`,
`q = cos d` is `ℓ = 3`, `q = 7cos²d − 1` is `ℓ = 4`, at `ω = √6, √12, √20`.

- **A spin-2 field cannot sit on a pole.** `h = sin²d·q` vanishes there for
  every `q`, so at the focus it is a ring of radius `0.198` about the antipode
  with `2.2e-06` *on* it — exactly where the scalar peaks. The same fact at the
  other end: **there is no point source of tidal shear.**
- **Pure shear.** Trace `0.0e+00`, and the first-order area change vanishes
  identically — measured `1.17e-08` against the `ε²h²/2` prediction `1.17e-08`.
- **Spin weight 2, directly.** The strain pattern is identical after a 180°
  rotation of the frame (`1.1e-15`) and inverted after 90° (2.00×).
- **The caustic is a quarter turn, not a flip.** One passage does *not* swap
  the stretch and compression axes: the outbound front is the Hilbert transform
  of the inbound one (correlation `+0.82` against `−0.35` for an inversion) —
  the Gouy shift, Maslov index 1, which the scalar shares. The axes do swap on
  the **round trip**: at `t = 2π` the field is minus its start (`+0.9974`).

```bash
python -m experiments.closure_ledger.spin2_tidal_probe
# Verdict: A_TENSOR_WAVE_CANNOT_SIT_ON_ITS_OWN_FOCUS  (7/7)
```

Scope: a spin-2 field on a **fixed** `S²`, not linearised GR on a spacetime —
2+1 gravity has no propagating tensor modes at all. Full write-up:
`docs/spin2_tidal_field.md`.

## Drawing the tensor wave in the embedding (continuous, not sampled)

Tidal ellipses are faithful but flat — they never touch the embedding, so they
say nothing about how shear meets a bulk. `viz/embedded_wave.py` projects
`h_ab` into a **continuous** surface deformation instead, and the projection is
forced rather than chosen.

A height field alone cannot do it: for `X = r n̂` the induced metric is
`g_ab = r²ĝ_ab + ∂_a r ∂_b r`, whose gradient term is *second* order, so at
first order a radial deformation is purely conformal (measured trace-free part
`1.8e-04` against a trace of `3.16`). Shape carries the trace and nothing else.
Adding the tangential part and demanding tracelessness fixes the rest:

```
X = (R + ερ) n̂ + ε∇W       ρ = −½ΔW       h_ab = [2∇₍ₐ∇_b₎W]^TF
```

**One potential carries both** — the shear is its trace-free Hessian, the shape
is minus half its Laplacian.

- **The theorem.** The induced metric perturbation of the drawn surface *is*
  the solved `h_ab`: component by component to `4.7e-04` of the peak, trace
  `1.1e-03`, off-diagonal `6.7e-08`.
- **The quadrupole is the textbook shape.** `ℓ = 2` returns `ρ = P₂(cos d)` to
  `1.0e-07` — the prolate–oblate picture, derived rather than drawn.
- **The free constant is a rigid translation** (`1.1e-16`), and `ℓ = 0` cannot
  appear at all: a spin-2 wave can never breathe the sphere's area.
- **Area holds** to `4.8e-06` at gain `1e-03` — second order, which is what
  trace-free means.
- **It reaches the bulk**: `74.8%` of the way to `R_outer`, `44.3%` toward
  `R_inner`, touching neither.
- **Principal-axis bars** at sparse points show the eigenvectors of `h_ab` —
  red for the positive eigenvalue, blue for the negative, length `∝ |λ|`. For a
  trace-free 2×2 the eigenvalues are `±√(h₊² + h_ˣ²)`, so the two bars are
  always equal; the stretch axis sits at `β = ½ atan2(h_ˣ, h₊)` and *swaps*
  between `ê_d` and `ê_ψ` wherever `h₊` changes sign.

```bash
python -m experiments.closure_ledger.embedded_wave_probe
# Verdict: ONE_POTENTIAL_CARRIES_SHAPE_AND_SHEAR  (7/7)

python scripts/geometrodynamics_v44_embedded_wave.py --still sheet.png
```

Scope: a representation of `h_ab` in the embedding, not backreaction — the wave
gains an extrinsic amplitude, it still does not act on the sphere. Full
write-up: `docs/embedded_tidal_wave.md`.

## What the refocus does to a trace-free deformation

Every point of the sphere runs its own principal-strain history, and on a
compact `S²` they all reconverge at the antipode at `t = π`. Because
`h = sin²d · q` vanishes at both poles for **every** `q`, so does `ḣ` — so the
effective density `ρ_E ∝ ḣ_ab ḣ^ab` cannot pile onto the focal *point*.

**The strains refocus on a ring, never on the point.** The density measured on
the antipode itself is `2.9e-08` of the peak, and the ring around it has a
radius that tracks the pulse width at `0.952 ×` across a factor of four in
width. Same fact that forbids a spin-2 point *source*, seen at the far end of
the trip.

- **The amplification is *not* a spin-2 effect.** The peak grows by a finite
  `2.18–2.32×`, which is tempting to read as the spin protecting itself from a
  singularity. A **scalar** refocused the same way amplifies by `1.86–2.07×`,
  and neither runs away as the pulse narrows — launch and focus are the same
  situation on a sphere. What belongs to the spin is the node and the ring.
- **A material patch changes shape without changing size**: area held to
  `1.9e-07` at gain `1e-2`, while its long axis matches the local stretch
  eigenvector to `1.000`. Near the focus that alignment is a question of patch
  *size* — a patch straddling the focal ring averages a sign change, and
  shrinking it recovers the eigenvector (`0.959 → 1.000`).
- **The area law fails first, and hardest, at the focus.** At the display gain
  the same patch moves its area by `1.93%` at mid-latitude and `25.95%` on the
  focal ring, while distorting `8.5×` harder. The fitted residual exponent is
  `ε^2.00`: the second-order term carries the *local gradient* of the field,
  which is the wavelength away from the focus and the pulse width on it.

```bash
python -m experiments.closure_ledger.focal_refocus_probe
# Verdict: THE_STRAINS_REFOCUS_ON_A_RING  (7/7)

python scripts/geometrodynamics_v45_focal_patches.py --still sheet.png
```

The v45 renderer shows three things and nothing else: the continuous deformed
surface, sparse principal-strain vectors, and two advected constant-area
material patches. The camera follows the wave, so the refocus faces the viewer.

Scope: no singularity forms here and none can — a linear field on a fixed round
background, with no backreaction and no bulk crossing rule. This establishes
the geometry such a rule would act on, and the amplitude at which the linear
description stops being trustworthy. Full write-up: `docs/focal_refocus.md`.

## A circle slice, and a bulk the wave can wrap through

Cut the sphere with the great circle through the source and its antipode. Along
it the geodesic distance is just `d = |σ|`, so **one circle carries both halves
of the wave** — two lobes running opposite ways and meeting head-on at `σ = ±π`.
Nothing is re-solved: the field is the 2-D solve sampled at `d(σ)`, verified
against the full `(θ, φ)` route to `1.4e-14`, so the slice inherits the sphere's
real caustic rather than the `2×` superposition a 1-D wave on a circle gives.

The slice is drawn inside the vacuole, with the obvious crossing rule: **a
radius past `R_outer` re-enters at `R_inner`**. So the wave that reaches up into
the bulk comes back *inside* the circle. That glues the boundaries, makes the
radial direction periodic, and turns the curve's home into a torus `S¹_σ × S¹_ρ`.
The wrap threshold is exactly the half-gap over the run's peak (`0.220420`
predicted, same to `3.8e-16` by bisection).

**And then nothing accumulates.** Driving the gain up buys crossings and never
buys charge:

| gain / threshold | unsigned | signed | winding |
|---:|---:|---:|---:|
| 1.6 | 2 | `+0` | `+0` |
| 3.6 | 4 | `+0` | `+0` |
| 5.0 | 6 | `+0` | `+0` |

The drawn curve is a **graph** `r = f(σ)` with `f` single-valued, so `f(π) =
f(−π)` and its degree as a map `S¹ → S¹` is identically zero — every outward
crossing of the seam is paid for by an inward one. **A height field cannot
wind.** Checked two independent ways (a signed crossing ledger and a degree from
unwrapped increments), agreeing at every gain and time.

Different pulses driven at one common gain cross the seam the *same* number of
times but put `0.155 → 0.033` of the circle on the far sheet — a `4.7×` spread
that is `2.61 ×` the pulse width for all of them. The far-sheet arc *is* the
pulse. None of them winds.

```bash
python -m experiments.closure_ledger.circle_slice_probe
# Verdict: THE_SEAM_IS_CROSSED_IN_PAIRS  (7/7)

python scripts/geometrodynamics_v46_circle_slice.py --still sheet.png
python scripts/geometrodynamics_v46_circle_slice.py --waves waves.png
```

Scope: the crossing rule is a representation choice, not a derived boundary
condition — nothing here makes the wave dynamically aware of the seam. What it
rules out is general: any representation drawing the bulk excursion as a height
over the slice inherits the same zero, so a stable topological object has to
come from a curve free to stop being a graph. Full write-up:
`docs/circle_slice_bulk.md`.

## The scaling at the seam is a choice

That fold — `r → r − gap` — carries a radial *offset* across unchanged, and the
two boundary circles do not have the same circumference. So a feature emerging
at `R_inner` keeps its full height on an arc shorter by `R_outer/R_inner`:
**the emerging wave was not the same wave**, and the scaling that made it so was
never chosen deliberately. The alternative is to translate in `ln r` instead:

| rule | map | aspect distortion | inward sheets |
|---|---|---:|---|
| `translate` | `r → r − gap` | **1.7027** | `0.22 → −0.30 → −0.82` |
| `conformal` | `r → r·(R_inner/R_outer)` | **1.0000** | `0.435 → 0.255 → 0.150` |

The conformal rule scales the offset with the boundary it crosses, so height and
arc shrink together and the feature returns a **faithful scaled copy**. Its
sheets are geometric, accumulating at the origin without reaching it — while
arithmetic sheets march straight through `r = 0` into negative radius. It pairs
with a multiplicative radial law `r = R_mid·exp(εu)`, positive by construction.

**And the choice decides what winding would look like.** A curve that genuinely
winds is a logarithmic spiral on the conformal seam: it returns to the same point
of the quotient magnified by `1.7027`, by the same factor from every starting
radius (spread `2.2e-16`). On the translate seam it returns *displaced*, with a
ratio that drifts by `0.217` with where it started — not a scale at all.

**What the choice does not change is the winding number.** Rebuilt conformally
it is still identically zero, at gains driving `274` unsigned crossings —
`ρ(σ)` comes from a single-valued function on the circle whichever coordinate
the seam translates in. The earlier negative result was not an artefact of an
arbitrary scaling; what the conformal rule adds is that the winding you cannot
have would have been *visible*, as `1.703` per turn.

```bash
python -m experiments.closure_ledger.seam_scale_probe
# Verdict: THE_GLUING_SETS_THE_SCALE_BUT_NOT_THE_TOPOLOGY  (6/6)

python scripts/geometrodynamics_v47_seam_scale.py --still sheet.png
```

Scope: the conformal rule is preferred on grounds of consistency, not because
any dynamics picks it out; `translate` stays the default so v46 is reproducible.
Full write-up: `docs/seam_scale.md`.

## The ring reaches across, but only a fold crosses

The wave is not only its focal pulse. A ring leaves the source, thins to `0.156`
at the equator, then **grows** to `0.933` at the focus — a factor of `5.97`,
following `1/√(sin d)` to `1.034 ± 0.071`. So: can *that* reach across the
vacuole, and can shrinking the gap or raising the energy make it intersect? Two
questions, two different knobs.

**Reaching across — yes, and the ring gets there first.** The threshold is
exactly `δ/max|u|`, and shrinking the gap buys *lead* on the focal pulse:

| δ | ε | dwell | spans from `d` | lead |
|---:|---:|---:|---:|---:|
| 0.26 | 0.40 | 0.062 | 3.142 | 0.209 |
| 0.16 | 0.80 | 0.371 | 2.437 | 0.811 |
| 0.09 | 0.80 | 1.000 | 1.832 | 1.413 |

At the bottom the ring spans from just past the equator for the whole converging
leg — a sustained state, not an instant at the focus.

**Intersecting — no, at any gap and any energy.** Swept with a real segment
test (validated against a limaçon, a lemniscate and a folded loop first):
`0` self-intersections at up to `346` seam crossings. A graph `r = f(σ)` is
*embedded* by construction — the winding obstruction seen from the side.

**What crosses is a fold.** Give the material tangential freedom,
`σ = σ₀ + λε ∂_σu`, and the map folds where `1 + λε ∂²_σu < 0`. Threshold
`λε = 0.012692` predicted, matched by bisection to `1.8e-12`, converged under
joint grid refinement — and it folds first **on the converging ring**, `0.0157`
from the antipode. Folding is necessary but not sufficient: crossing-without-fold
is `0` at every drive, fold-without-crossing happens.

**And the two knobs are orthogonal.** The fold threshold is gap-independent
(spread `0.0` across a fivefold range of `δ`) while the span threshold scales
with `δ` directly. What the fold threshold scales with is the pulse:
`λε ≈ 0.385 w²`, spread `3.7%`. **Reducing the shell separation changes when the
wave arrives; it changes nothing about whether it can cross itself.**

```bash
python -m experiments.closure_ledger.ring_and_fold_probe
# Verdict: THE_RING_REACHES_BUT_ONLY_A_FOLD_CROSSES  (8/8)

python scripts/geometrodynamics_v48_ring_and_fold.py --still sheet.png
```

Scope: `λ` is a modelling choice, not derived from the scalar equation, and
`λ = 0` is exactly the earlier height field. Full write-up:
`docs/ring_and_fold.md`.

## What a drawn wave has to obey

An audit of the height-field representation against the physics it stands in
for. Three objections, three different fates.

**"The antipode moves before the front arrives" — not in the solve.** The
amplitude there runs `5e-133 → 6e-40 → 7e-17 → 0.93`, climbing 130 orders of
magnitude in step with `t → π`; signal ahead of the front sits at the scheme's
`1e-07` noise floor.

**But a constant never leaves — and that is a real defect.** A Gaussian carries
a monopole `w²/4` (`0.008056` measured against `0.008100`), and `ℓ = 0` has
`ω = 0`, so a closed surface can never shed it. Not an early response — ahead of
the front the higher modes cancel it exactly — a **permanent** one: every point
carries a time-averaged displacement of `w²/4`, and the quietest instant of a
whole run still leaves `max|u| = 0.094`. *Nothing moves early, and nothing ever
stops moving.* EM forbids monopole radiation and gravity forbids `ℓ = 0, 1`, so
this mode is outside the analogy, not a blemish on it — the spin-2 field's DC is
`2e-06` against the scalar's `8.3e-03`.

**The fix has to stay inside the pulse.**

| source | monopole | far side before arrival |
|---|---:|---:|
| gaussian | `+8.1e-03` | `1.5e-20` |
| gaussian, monopole-free | `-4.3e-07` | `1.4e-05` |
| compact, monopole-free | `-3.1e-07` | `3.8e-22` |

Subtracting the mean leaves `−w²/4` at the antipode — exactly what it removed.
A wider Gaussian corrector costs four orders of far-side quiet. Compactly
supported bumps give both at once, with a far side that is *exactly* zero.

**"Height should grow as the ring compresses" — right physics, hidden.**
`A²·(circumference)` is conserved to `8.4%`, but `1/√(sin d)` is flat across the
middle (`72%` of the trip shows nothing) and, because the launch is itself a
focus on a compact surface, the focal height comes in at `0.935×` the *launch*
height. Unbounded growth belongs to an open geometry.

**"A deforming surface is the wrong representation" — right, and measurably
so.** The energy density is `u̇² + |∇u|²`, so the constant offset displaces every
point and carries exactly `0.000` energy against the pulse's `0.497`. The most
global feature of the drawn shape is invisible to the physics.

```bash
python -m experiments.closure_ledger.wave_constraints_probe
# Verdict: NOTHING_MOVES_EARLY_AND_NOTHING_EVER_SETTLES  (8/8)
```

Full write-up: `docs/wave_constraints.md`.

## Draw the wave as vectors, and they intersect

Four rounds established that a height field cannot wind and cannot cross itself.
All four were about the same object — **the graph of the displacement's tips**,
`r = f(σ)`, which is embedded by construction. Draw the **vectors themselves**
and the obstruction is gone, for a classical reason: *neighbouring normals to a
curve meet at its centre of curvature*. A normal of length `L` crosses its
neighbours as soon as `L` exceeds the local radius of curvature `ρ = 1/κ`.
Nothing is added by hand.

The same wave, the same instant, two objects: the graph gives **0**
self-intersections, the normal field gives **520**.

**And the threshold is the ring concentration, finally visible.**

| t | `ρ_min` | first drawn crossing |
|---:|---:|---:|
| 1.2 | 0.1408 | 0.492 |
| 2.6 | 0.1087 | 0.340 |
| 3.0 | **0.0540** | **0.189** |

The converging ring sharpens its own surface by `2.61×`. The focusing shows up
not as height — which barely moves, and never beats the launch — but as
*curvature*, which is what decides whether the normals meet.

**The reset is a second, separable mechanism.** A normal leaving through
`R_outer` re-enters at `R_inner` at the angle where it left, shooting outward
from deep inside the annulus: at the focus, `306 → 398` crossings at `L = 0.35`.

**And the gap matters again** — what `ring_and_fold.md` had severed. The vector
length *is* what spans the gap, so they are one knob. At `L = δ`:

| δ | normals alone | with reset |
|---:|---:|---:|
| 0.40 | 382 | 382 |
| 0.16 | 0 | **180** |
| 0.09 | 0 | **472** |

At the tightest gap the normals are too short to reach each other, but almost all
of them wrap and the stubs cross everything. **Reducing the shell separation now
produces intersections rather than being unable to.**

```bash
python -m experiments.closure_ledger.normal_field_probe
# Verdict: THE_NORMALS_INTERSECT  (6/6)
```

Scope: the vector length is a display choice; the directions and curvature are
the surface's own. Full write-up: `docs/normal_field.md`.

## Focus the congruence to the pinch, and let the equations decide

A singularity is a failure of evolution, not a bright spot — so drawing one as a
glowing dot assumes the answer. The object that does not is a **congruence with
a deforming cross-section**: integrate `F̈ = ½ḧF` from `F(0) = I` and watch
`J = det F`, the cross-sectional area of the bundle. `J → 0` is a caustic **of
the map**, and says nothing about the metric.

Every line the renderer draws is a principal axis of `F`, so the spin-2 tensor
builds its own bulk structure — no normals and no invented vector length.

**Three cases that were being conflated, now separated:** an ordinary focus
where `J` dips and recovers; a caustic where `J` reaches zero and the geometry
stays regular; and a curvature singularity, which **cannot occur in this
program** — the background is a fixed round `S²` with curvature `1` at every
time, with no Einstein equation and no backreaction.

**Of passage, singular termination, and finite-radius reconnection, the
equations give passage.** At the source ring `J` crosses zero with slope
`−17.877`, converged to `−17.836` under a halved timestep — a tangency would
have driven it to zero — plunges to `−471`, and the evolution continues with the
solver's invariant unmoved at `2.5e-14`. The other two were never available, and
for different reasons: termination needs the geometry to fail, and reconnection
needs the congruence to act back on something, but each point's `F` is driven
only by the external `h`, never by its neighbours. Not "we did not find them" —
"this program could not have produced them".

**Two rings, a factor of ten apart** (window `1.2π`):

| ring | threshold peak strain |
|---|---:|
| source | **0.026** |
| converging (`d > π/2`) | **0.247** |

Even at `0.247` the antipodal crossing only *grazes* zero and the depth of that
excursion does not converge. So "does focusing drive the neck radius to zero?"
— barely, and only at a strain nothing physical would reach.

**The neck is a ring, and spin weight is why.** `h = sin²d·q` vanishes at both
poles, so the tidal field is identically zero *at* the antipode and the
congruence there is never driven. The neck sits at radius `0.44 w`, the same
ratio across a `3.3×` range of pulse width to within one grid cell. A spin-2
focus has no centre.

The transverse Raychaudhuri equation turns out to be *the same statement* as the
deviation equation here — the residual is an algebraic identity, so its
`6.7e-15` is round-off rather than a convergence result. The content is that the
**Ricci term vanishes identically** for trace-free `h`, so none of the focusing
is matter and all of it is shear-squared: second order in the amplitude, which
is why a weak wave barely focuses and a strong one does so abruptly.

**Two errors this round had to catch first.** Forming `ḧ` from a time difference
seeded with `h(−dt) = h(0)` injects a spurious `½ḣ(0)` impulse whose `1/dt²`
cancels against the update's `dt` — so **refinement reproduces it instead of
exposing it**, and the two forms converge to `0.628` (no caustic at all) against
`−12.84`. And a Gaussian launch deforms the far side at `t = 2.00` against a
causal bound of `2.77`, *grid-converged*: an analytic tail, fixed only by the
compact bump `(1−x²)⁴`.

```bash
python -m experiments.closure_ledger.congruence_probe
# Verdict: THE_CAUSTIC_IS_A_PASSAGE  (10/10)
```

Scope: `gain` is a strength dial, reported as a peak strain throughout; the
deviation equation is exact in `ξ` and linear in `h`. Full write-up:
`docs/congruence.md`.

## Can a detached shell do the throat's exotic work for it?

**Scope: Einstein gravity, Darmois–Israel thin shells, spherical symmetry,
vacuum bulk, `G = 1`.** The dimension is a parameter — `D = 4` is the regression
case, `D = 5` (Tangherlini) the one this program cares about.

A throat needs negative surface energy. Could a *detached* closed shell, glued
with the opposite orientation, supply the exotic-looking restoring stress while
itself being ordinary matter? That is three questions wearing one coat, and
**they do not agree**, so they are reported separately.

**The orientation is derived from the gluing, not set by hand.** Each side
retains either the INNER (`r ≤ R`) or OUTER (`r ≥ R`) branch, and `ε` follows.
With `σ = −(D−2)(ε₊β₊ − ε₋β₋)/(8πG R)` that gives **four** gluings:

| `−` | `+` | `η` | what it is | `σ` |
|---|---|---:|---|---|
| INNER | OUTER | `+1` | ordinary bubble | either sign |
| OUTER | OUTER | `−1` | **minimal surface** | `< 0` always |
| INNER | INNER | `−1` | **maximal surface** | `> 0` always |
| OUTER | INNER | `+1` | anti-bubble | either sign |

So `η = −1` alone decides nothing — it covers two gluings whose forced signs are
*opposite*. What is forced is sharper: a minimal surface has
`σ = −(D−2)(β₊+β₋)/8πGR < 0` and a maximal one the same with the other sign,
both identities, neither violated once in 40,000 random Tangherlini / de Sitter
/ charged pairs across `D = 4, 5, 6`.

**The dichotomy that follows is the answer.** A detached surface that
*connects* to the throat's region does so through a minimal surface and is
necessarily exotic. One that is non-exotic by its gluing is a **maximal**
surface — it caps off on both sides, shares no bulk with the throat, and is
non-exotic precisely *because* it is disconnected, so it cannot support
anything. Within Einstein–Israel spherical thin shells the exotic matter is
relocated, never removed.

**The three observables still disagree on one system.** An ordinary bubble has
`σ = +6.2e-05`, its screening does push the throat outward, and `σ_throat =
−0.027` regardless. Screening raises the critical `β²` monotonically
(`−1.083 → −0.652`), enlarging the stability window but never reaching
`β² ≥ 0`.

Two things are carried in the output rather than buried: the shell's effect is a
potential **gradient**, not an equilibrium-consistent force (fixed throat rest
mass, no equation-of-state response); and Birkhoff's vanishing `∂²V/∂a∂b` is
**structural**, imported the moment the intervening region is written with a
constant mass parameter. What *is* measured is that a family of shells differing
**701×** in surface density leaves the throat bit-for-bit unchanged — which
establishes no separation-dependent coupling *in this model*, not that every
spherical trapped resonator is impossible. `ℓ ≥ 2` is where such a coupling
would have to live.

```bash
python -m experiments.closure_ledger.shell_junction_probe
# Verdict: CONNECTED_MEANS_EXOTIC  (10/10)
```

Imported rather than derived: Birkhoff, the Darmois–Israel formalism, and `β²`
as a free parameter at the equilibrium. Stiffnesses are stiffnesses, not
normal-mode frequencies — no kinetic metric is derived. Full write-up:
`docs/shell_junction.md`.

## Where the two-shell coupling starts

**Scope: the static Newtonian (Laplace) two-shell model** — the weak-field
analogue of the junction problem. It establishes the shell-theorem/multipole
structure. **Birkhoff is a GR result and remains what PR #249 relies on**; the
`ℓ = 0` statement here is its Newtonian analogue, not a replacement.

In this model the **monopole mutual stiffness vanishes** while higher angular
multipoles couple, with the coupling **suppressed geometrically by separation**:

```
∂²U/∂α∂γ  =  G m_b m_a · ℓ(ℓ + D − 3) · b^ℓ / a^{ℓ+D−3} · κ_ℓ(D)
```

The prefactor `ℓ(ℓ + D − 3)` is **the eigenvalue of the Laplacian on `S^{D−2}`**,
so the `ℓ = 0` decoupling is that zero eigenvalue — in every dimension, not as a
four-dimensional accident. **`D` is derived, not assumed**: each case is checked
by brute force *in its own dimension*.

| | `D = 4` | `D = 5` — the BAM case |
|---|---|---|
| sphere / kernel | `S²`, `1/r` | `S³`, `1/r²` |
| `κ_ℓ(D)` | `1/(2ℓ+1)²` | `1/(ℓ+1)` |
| closed form | `ℓ(ℓ+1)(b/a)^ℓ/(a(2ℓ+1)²)` | `ℓ(ℓ+2)/(ℓ+1) · b^ℓ/a^{ℓ+2}` |
| brute force agrees to | `9e-06` | `3.3e-04` |

with `ℓ = 0` vanishing to `1.7e-12` in the `D = 5` control, whose undeformed
value reproduces the four-dimensional shell theorem `−G m_b m_a / a²`.

**Coupling starts at `ℓ = 2`, not `ℓ = 1`.** A first draft claimed "everything
`ℓ ≥ 1` couples", having checked translation invariance of the **area** — the
wrong quantity. Run on the mutual **energy**, a rigidly displaced inner sphere
leaves it at exactly `−G m_b m_a / a` (Newton's shell theorem, held to `1e-15`),
so the cross-derivative for rigid translations is `8.3e-13`: **the translation
mode does not couple.** The pure `P₁` *shape* mode is a different object and
does, at `1.78e-02`. So `ℓ = 0` decouples by the vanishing eigenvalue, the
`ℓ = 1` translation mode by the shell theorem, and genuine coupling begins at
`ℓ = 2` — which is what #249 guessed.

**And it is screened.** `(b/a)^ℓ` gives a factor of **544** from `ℓ = 1` to
`ℓ = 8` at `b/a = 0.4`. The multipoles carrying a spin-2 signal are the ones
separation suppresses hardest.

**The shear response is an input.** A perfect-fluid shell resists area change
and nothing else; resisting *shape* change at fixed area needs an elastic
modulus no equation of state supplies. Carried explicitly, never fitted.

```bash
python -m experiments.closure_ledger.multipole_coupling_probe
# Verdict: COUPLING_STARTS_AT_ELL_TWO  (10/10)
```

Not Regge–Wheeler, no quasinormal frequencies, no claim that a trapped `ℓ ≥ 2`
resonance exists. Full write-up: `docs/multipole_coupling.md`.

## One conserved wave, seen in pieces

**Scope: kinematics and bookkeeping on a fixed round `S³`.** The wormhole is an
*identification map* — a pair of mouths and a time offset — not a solution; `κ`
and the delay `Δ` are parameters. PR #249 priced such a throat (a minimal
surface has `σ < 0` identically), and this round **inherits that bill without
paying it**. No Einstein equation, no backreaction, no cross-section.

An emitter fires a shell; while it expands, a *collapsing* shell sweeps past it
and lands on a receiver, which recoils. Locally, two unrelated events. Globally,
one object:

```
E ──expand──▶ M_future ──throat, Δt < 0──▶ M_past ──collapse──▶ R
```

**The staging is geometry, and the same fact places both ends.** A geodesic
sphere at distance `χ` has area `4π sin²χ`, so a shell refocuses exactly at
`χ = π` — checked against uniform sampling of `S³` through the enclosed volume
and scored against the estimator's own binomial standard error (worst
`z = 2.02`). That puts the future mouth at the emitter's antipode *and* the
receiver at the past mouth's, which is the only place the arriving shell is
genuinely **collapsing**: `dA/dχ = 4π sin 2χ < 0` all the way in. Against a
receiver displaced to `χ = 1.2` the same wave is still *expanding* when it
lands — so "collapse" is a measured sign, not a caption.

**The content is the self-consistency.** A closed timelike loop carrying a
**linear** wave has exactly one amplitude, `A = A_src/(1 − κ)`, matching brute
iteration of 4000 round trips to `7.1e-13` even at `κ = 0.99` (amplification
`×100`). It is unique for **every** `κ ≠ 1`, and exists even where `|κ| > 1` and
the *iteration* diverges — divergence of a summation method is not absence of a
solution. Nothing is tuned and no paradox is available.

**That is a fact about linear equations, not about time travel**, and it is
demonstrated rather than asserted: the same loop with a quadratic return
`A = S + κA²` has two solutions below a source threshold of `1/4κ` and none
above it.

**The picture measures rather than illustrates.** Every drawn point sits at
geodesic distance `χ` from its centre to `6.7e-15`, and the shadow's screen
extent is proportional to `sin χ` with **one constant to `3.6e-16`** — that is
`√(A/4π)`, so the apparent size in the figure *is* the area law.

Getting there meant replacing the projection, which is this round's real
correction. A stereographic chart is unbounded at its own pole, and a shell
launched from a point sweeps **all** of `S³`, so whatever pole is chosen the
shell crosses it once: the radius grows as `2/ε` across four decades and never
converges. The first renderer projected from `q₃ = +1` — the emitter's own
position — so the emitter was a division by zero and never got drawn. The
orthographic **shadow** replaces it: bounded by the unit ball everywhere, with
the discarded coordinate kept as brightness. It is **2-to-1**, so a crossing on
screen is not a crossing on `S³`, and nothing here rests on one.

**And what closing the ledger is not evidence for.** Flux conservation through
the throat is *put in* when the mouths are identified, so the `1.1e-16` residual
checks the arithmetic. What it does establish is the identification: neither
local event conserves anything alone — the emitter loses flux, the receiver
gains it — and the pair closes exactly. The delay decides the entire story and
changes no conserved quantity: amplitude spread across `Δ` from `−2` to `−10⁴`
is `0.0`.

```bash
python -m experiments.closure_ledger.wormhole_ledger_probe
# Verdict: ONE_WAVE_AND_LINEARITY_IS_WHY  (10/10)

python scripts/geometrodynamics_v51_wormhole_ledger.py --still v51.png
```

Full write-up: `docs/wormhole_ledger.md`.

## Pair creation is a collision, not a focus

**Scope: kinematics on a fixed round sphere.** The Breit–Wheeler threshold and
cross-section are **imported QED**; rays-as-photons is a correspondence; **no
rate is computed** (a rate needs a photon number density a classical amplitude
does not supply); and calling the crossing chord a throat is this program's
reading, with PR #249's exotic-matter bill inherited, not paid.

Every earlier wave round drew **one** wavefront converging on its antipode and
treated the caustic as the interesting event. That was the wrong *quantity*. A
caustic is where the amplitude gets large; pair creation is a threshold on an
invariant:

```
γ γ → e⁺ e⁻      s = 2 E₁E₂ (1 − cos θ) ≥ (2m)²
```

`E` is what focusing raises; `θ` is what a collision supplies. So focusing is
**neither sufficient** — collinear momenta have `s = 0` identically, still
exactly zero after amplifying by `10¹²` — **nor necessary**, since two beams
crossed head-on with no focusing anywhere clear threshold as soon as `E ≥ m`.
The honest complication is recorded rather than buried: a spherically converging
front *does* contain diametrically opposed rays, so the real distinction is
**independence of the sources**, not brightness.

**The invariant is an identity of geodesic triangles.** Two sources a geodesic
distance `δ` apart, both firing:

```
1 − cos θ = (1 − cos δ)/sin²t     ⟹     s(t) = 4 E₁E₂ sin²(δ/2)/sin²t
```

verified to `2.0e-14` against a control that never uses the law of cosines — the
crossing point solved as a linear system, the momenta built as great-circle
tangents in the embedding — and checked on `S²` **and** `S³`, because a geodesic
triangle lies in a great 2-sphere whatever it is embedded in.

**Which makes the collision head-on twice.** `s(t)` is **U-shaped**: maximal at
*both* ends of the crossing window and minimal at the equator, where the opening
angle is exactly `δ`. The moment the wavefronts are largest is the moment the
invariant is smallest. So the threshold cuts **two disjoint windows, never one**:

| `E/m` | windows | where |
| ---: | ---: | --- |
| 0.6 | **0** | even head-on, `s_max = 4E² < 4m²` |
| 1.0 | **0** | zero width — touched only at the two head-on instants |
| 1.4 | **2** | `[0.210, 0.296]` and `[2.845, 2.932]` |
| 6.0 | **1** | merged, above `E = m/sin(δ/2) = 4.797 m` |

**And only the far window is a collision of independent waves.** The near one
sits on top of the sources — the fronts have travelled `0.296` against a
separation of `0.42`, so nothing there has propagated independently. The far one
is reached only after each front has crossed a half-circumference: a factor of
**9.6** in path length. *That* is why the second interaction has to be
antipodal, and it is measured rather than staged.

**One further trap, caught.** The momenta are exact — perpendicular to their own
front to `1e-15`, matching the closed form to `2e-13` — but a figure shows their
**projection**, and projection does not preserve angles. Measured off the
picture the opening angle is wrong by up to **67.4°**, and differs by up to
**56.4°** between the two crossing points, whose true opening angle is
*identical*. The renderer therefore draws the angle in the plane the two momenta
span, and labels the arrows on the sphere as projections.

```bash
python -m experiments.closure_ledger.pair_creation_probe
# Verdict: THE_SECOND_INTERACTION_MUST_BE_ANTIPODAL  (9/9)

python scripts/geometrodynamics_v52_pair_creation.py --still v52.png
```

Full write-up: `docs/pair_creation.md`.

## Two closed histories, sewn at one interaction

**Scope: a counting result on a kinematic skeleton.** Throats are identification
maps with **given** mouths and delays; PR #249 priced a connected throat as
necessarily exotic and this round puts **two** on the books, paying for neither.
**No action principle, no field equations, no topology change, no dynamics, no
rate, and no worldline** — whether a particle trajectory is the locus where
expanding and collapsing components stay consistent is untouched. Conjugacy
`Q_A + Q_B = 0` is a label, carried and checked, never derived.

PR #251 built one closed history; PR #252 showed pair creation needs two
independently propagated waves. Sewing two closed histories at one shared
interaction asks the question that comes *before* attempting topology change:
**is that event constrained by the closure conditions, or still put in by hand?**

Every leg is null, so a history closes in coordinate time — **on the principal
branch** exactly on a geodesic ellipsoid, the locus whose summed distance to its
two mouths is `|Δ|`, feasible on `[d, 2π − d]` and checked against 40,000
uniform samples of `S³`. That branch scope is load-bearing: `d` is the
*principal* geodesic distance, and off the principal branch a mixed leg
assignment fixes the **difference** of distances — a hyperboloid, not an
ellipsoid. Inside the band the principal branch is the only feasible one, so the
rest is principal-branch **by construction of its prior**. Discreteness survives
per branch. The global system is then

```
|c|² = 1                                normalisation
d(S_A, c) = t_C − τ_A                   C lies on front A
d(S_B, c) = t_C − τ_B                   C lies on front B
d(c, M_A⁺) + d(M_A⁻, c) + Δ_A = 0       history A closes
d(c, M_B⁺) + d(M_B⁻, c) + Δ_B = 0       history B closes
```

**Five equations, five unknowns.** Solved blind from random starts, every root
found sits at **full Jacobian rank 5** (12 of 12). Precisely: that shows each
root *found* is locally isolated — not that all roots were found, nor that the
event is unique. Existence is restrictive too: only about half of configurations
*drawn from this module's prior* admit a closed pair-history.

**Removing a wave is a dimensionality control, not physics.** Deleting one
scalar equation from a square nondegenerate system drops the rank by one, for
*any* equation — deleting a **closure** instead gives the identical result:

| | equations | rank | solutions |
| --- | ---: | ---: | --- |
| both waves | 5 | **5** | isolated events |
| wave B removed | 4 | **4** | a **one-parameter family** (~159 sampled) |

So this establishes nondegeneracy, and is **not** evidence that pair creation
needs two photons — that content lives in the invariant `s`, which needs two
independent momenta. What survives as interesting is only the direction: the
solutions do not vanish, they stop being isolated, and dropping the constraint
can even **create** solutions where two waves admitted none.

**And in this model the conjugate pair needs two distinct throats.** One shared
throat fails both ways: traversed oppositely, history B's closure demands a
*negative* sum of leg lengths — infeasible on **every** branch, the one
conclusion here independent of branch scope; traversed the same way on the same
branch, the two closure equations coincide and the rank drops to 4. That second
half is scoped to the minimal single-pass model and **scanned** over branches
rather than argued: no counterexample at winding ≤ 1.

**The non-circularity check is the one that matters.** With the delays free the
nullity is **measured on the actual 5×7 Jacobian** (rank 5, nullity 2), and
**100% of 345 sampled** events on both fronts close by choosing `Δ` afterwards.
The whole result rests on the throat being data — a version that solved for the
delays would constrain nothing while looking identical from outside.

**Closure constrains where; the invariant decides whether** — with two warnings.
That no event clears `s ≥ 4m²` at `E = m` is **forced, not measured**
(`s ≤ 4E²`, equality only at exactly head-on, measure zero); and every fraction
reported is conditioned on an arbitrary prior, so they are **regression
diagnostics, not predictions**.

**And it is now branch-complete.** A leg length is `ℓ ≥ 2πk` and the legs must
sum to `|Δ|`, so `k₁ + k₂ ≤ ⌊|Δ|/2π⌋`: the feasible branch set is **finite and
bounded by the delay**, with no cutoff to choose (brute-enumerated to winding 11,
zero violations). Taking the union over *every* feasible branch pair, all 18
roots found are still at full rank 5; the shared-throat obstruction survives **51
distinct branch pairs** with none restoring full rank, so that check is now
exhaustive rather than scanned; and the delay dependence is untouched.

**The surviving claim:** *with given throat data, intersecting two null fronts
with two independent closure hypersurfaces produces locally isolated candidate
events, branch-completely; removing one front constraint restores a continuous
degree of freedom.*

**And that is where rank counting ends.** It cannot supply a quantity that
*vanishes* when a source is removed rather than merely becoming underdetermined
— that needs a field.

```bash
python -m experiments.closure_ledger.pair_history_probe
# Verdict: BRANCH_COMPLETE_AND_STILL_DISCRETE  (14/14)

python scripts/geometrodynamics_v53_pair_history.py --still v53.png
```

Full write-up: `docs/pair_history.md`.

## A solved field reproduces the branch ledger — and signs it

**Scope: a linear scalar field on a fixed background.** The throat is still an
**identification map**, not a solution, with PR #249's exotic-matter bill
inherited and unpaid. **No backreaction, no topology change, no rate, and no
two-source invariant yet.** Repeated throat traversals are PR #251's fixed point
and are not redone here.

PR #253 ended by conceding that rank counting had reached its limit. Before
building the quantity that replaces it, the ray ledger has to survive contact
with a solved field. It does — and the branches turn out to be **exact support**
rather than stationary-phase contributions.

On the Einstein static universe the scalar Laplacian has eigenvalues `−n(n+2)`
and `R = 6`, so the **conformally** coupled massless field has

```
ω² = n(n+2) + 1 = (n+1)²        ⟹        ω_n = n + 1
```

**Integer frequencies**, so the retarded Green function is exactly periodic and
is a sum of images:

```
G(χ,t) = 1/(4π sin χ) [ Σ_k δ(t − χ − 2πk) − Σ_k δ(t + χ − 2πk) ]
```

A truncated **mode** sum (which never sees an image) against the closed-form
**image** sum (which never sees a mode): **`8.3e-13`**.

**The field's support is the ray ledger.** Peaks read off the solved field land
on PR #253's branch times to `3.0e-04` — half a grid cell, so grid-limited —
with every branch matched and no peak spurious.

**The amplitude is PR #251's shell law, derived rather than imposed.** That round
set `A ∝ 1/sin χ` by conserving energy across a shell of area `4π sin²χ`; here
`peak × sin χ` is the same constant at every `χ` to **`7.0e-16`**.

**And the field supplies phases the ray ledger could not.** Every arrival carries
a sign, and it is `(−1)^m` with `m` the number of focal crossings — the antipode
at `t = π`, the source point again at `t = 2π`, and so on. That is the **Maslov
index**, and **12 of 12** signs agree:

| `t` | branch | crossings | field sign | Maslov |
| ---: | --- | ---: | :---: | :---: |
| `0.700` | short, `k=0` | 0 | `+` | `+` |
| `5.583` | long, `k=0` | 1 | `−` | `−` |
| `6.983` | short, `k=1` | 2 | `+` | `+` |
| `11.866` | long, `k=1` | 3 | `−` | `−` |

A path-length ledger gives arrival times and has no way to produce a sign. This
is the first quantity in the arc the ray picture could not in principle have
carried.

**And the ledger belongs to the *conformal* field specifically.** The minimally
coupled field has `ω = √(n(n+2))`, irrational, so no images and no sharp
branches — **63%** of the peak amplitude sits *between* the arrivals, against
`4.0e-08` for the conformal field. PR #253 never said which field its ledger
described, because rays cannot tell the two apart.

**The throat reproduces the closure condition.** A through-throat contribution
lands at `ℓ₁ + Δ + ℓ₂`; setting `Δ` to minus a branch-pair sum — exactly PR
#253's closure condition — puts an arrival back on the emission event, **9 of
9**, with the field adding the sign as `η` times the two Maslov factors.

```bash
python -m experiments.closure_ledger.field_solve_probe
# Verdict: THE_FIELD_REPRODUCES_THE_LEDGER_AND_ADDS_ITS_PHASES  (7/7)

python scripts/geometrodynamics_v54_field_solve.py --still v54.png
```

Full write-up: `docs/field_solve.md`.

## The mouth transfer solved for, not applied

PR #254 solved the field but kept the mouth relation on the **outside**:
`φ(M⁺,t) = η φ(M⁻,t+Δ)` was applied *to the free branches after they were
computed*. One traversal, by construction — a post-processing step cannot notice
that what it re-emits will come back. Here the relation enters the equation that
is solved:

```
a(ω) = ηκ e^{−iωΔ}[ S(ω) + T_d(ω) a(ω) ]   ⟹   a = ηκ e^{−iωΔ} S / (1 − L)
```

**Scope first**, because the resolvent being exact for a model says nothing about
which model: this is a **self-consistent rank-one mouth-transfer model**, *not* a
throat boundary operator and *not* a quotient of the manifold. It relates field
values through the free Green function — no normal-derivative (flux) matching, no
reflected channel (the mouth scattering object is `1×1` where a flux-conserving
two-mouth junction needs at least `2×2` unitary), and power out over power in is
`κ²` exactly, so for `κ < 1` it is lossy and not an identification of anything.

What *is* claimed is narrow and holds.

**The resolvent is the sum over every traversal** — against an explicit walk over
400 of them, `3.5e-18` — and PR #254's answer is its `n = 0` term, whose relative
error is *exactly* the round-trip gain `|L|`, as an identity rather than a fit.

**The branch series sums in closed form, and its poles are the spectrum.** The
short-way images all carry Maslov factor `+1` and the long-way images all carry
`−1`, so the winding sum is two geometric series:

```
Σ_b s_b e^{−u ℓ_b} = (e^{−uχ} − e^{−u(2π−χ)}) / (1 − e^{−2πu}),   u = γ + iω
```

verified term-by-term to `2.7e-15`. As the regulator goes, its poles sit at
`ω = 1, 2, 3, …` — the conformal ESU eigenfrequencies — with residues equal to
the mode functions over `2ω`. **The image representation and the mode
representation are one function**, which is the strongest statement in the arc
that the branch labels are a representation rather than an approximation.

**The solve adds events, not amplitudes.** The solved waveform equals the sum
over history words `(a, c₁…c_n, b)` to `5.4e-06`; at echo times
`ℓ_a + Δ + n(ℓ_c + Δ) + ℓ_b` that no one-traversal word reaches, the solved field
stands `3.3e+12` above the control, on a `κⁿ` ladder, each echo signed by every
Maslov factor in its word. Those are arrivals at times PR #254's ledger does not
contain.

**The primitive is indexed by a pair of branches.**

```
K_ab(ω) = ηκ · s_a A₁ e^{−u ℓ_a} · e^{−iωΔ} · s_b A₂ e^{−u ℓ_b}
```

`K_ab` carries the phase `e^{−iω(ℓ_a + Δ + ℓ_b)}`, so PR #253's closure condition
is *exactly* the statement that it does not depend on `ω`: closed pairs have band
coherence `1.000`, every other pair below `0.091`. And the reason the **pair** is
the primitive is that the amplitude factorizes over that index while the
condition does not: at `Δ = −(χ₁+χ₂+4π)` **three** pairs close inside the **nine**
any single-index rule would have to admit. An anti-diagonal, not a rectangle.

**Rank counts transfer channels, not histories.** One throat already carries
`144` distinct `(a,b)` histories and `K` is still rank one — an outer product
counts independent *separable channels*, and one value-feedback throat supplies
one. A second throat adds a second, in a shared topological branch-label basis
(checked, since the two throats have different `χ` on every leg), and the
interference between the channels is a full fringe, bilinear, identically zero
without either.

| configuration | singular values of `K` (normalised) |
| --- | --- |
| one throat | `1`, `1.3e-16`, `4.6e-17` |
| two throats | `1`, `0.542`, `1.0e-16` |

**And existence, convergence and stability are three different conditions.**
`1/(1−L)` exists for any `L ≠ 1`, `|L| > 1` included; `|L| < 1` is only the radius
of `Σ Lⁿ`, and that radius does not depend on the delay at all. Stability is
`Im ω > 0` for every root of `D(ω) = 1 − L(ω)` in complex `ω`, and the coupling
*displaces* the bare poles `ω = m + iγ` by
`δ_m = −ηκ e^{−imΔ} sin(md)/(4π² sin d)` — matched to `2.2e-04` — whose imaginary
part goes like `sin(md) sin(mΔ)` and **changes sign with the mode**. Stability is
phase-sensitive; no bound on `|L|` can see it:

| `Δ` | `κ_series` | `κ_stability` | ratio |
| ---: | ---: | ---: | ---: |
| `1.0` | `0.7619` | `0.7710` | `1.012` |
| `π` | `0.7619` | **`3.0336`** | **`3.98`** |

At `κ = 1.520` in that gap the traversal series diverges to `1.3e+119` while the
solve is finite and the least-damped pole is still at `Im ω = +0.0145`. Solving
and summing are not the same operation.

Still put in: a linear field on a fixed background and `κ` by hand. When
`Δ + ℓ_c < 0` the loop is closed in time and `1/(1−L)` is a self-consistency
condition rather than a history sum. The two-throat fringe is a *throat–throat*
interference, **not** roadmap step 3's two-source invariant, and the
flux-conserving boundary operator is the next step, not this one.

```bash
python -m experiments.closure_ledger.branch_coupling_probe
# Verdict: THE_TRANSFER_IS_SOLVED_FOR_AND_THE_PRIMITIVE_IS_A_PAIR  (10/10)

python scripts/geometrodynamics_v55_branch_coupling.py --still v55.png
```

Full write-up: `docs/branch_coupling.md`.

## Conservation is not stability

PR #255 owed a boundary operator. A point-supported throat is a **self-adjoint
extension** of the Laplacian on `S³ ∖ {M⁺, M⁻}`, parametrized by `U(2)`; writing
the boundary condition as the *pair* `B φ^reg = C q` — general enough to hold
#255's relation, which is not of the form `φ^reg = A q` — the mouth-active
spectrum is `det(C − BΓ) = 0`.

**It is definable at all.** `G(χ,ω) = sin(ω(π−χ))/(4π sin χ sin(πω))` in closed
form: real on the axis, poles exactly at `ω = n+1`, and *finite* at the antipode
where the numerator's zero cancels `sin χ`. It matches #255's branch series to
`6.3e-12`, and its short-distance split `1/(4πχ) + g(ω) + O(χ)` has remainder
first order in `χ`, so the subtraction a point interaction needs is forced.

**Hermiticity is exactly flux conservation.** The current through a small sphere
at mouth `j` is `Im(q_j* φ_j^reg)`, so the total absorbed is `Im(q† A q)` — zero
for *every* `q` when `A = A†`, at `1.8e-16` over 200 draws, against a median
`0.54` for a non-Hermitian control.

**And that is all it buys.** `Γ` is real symmetric for real `λ = ω²` of *either*
sign, so `λ` is real — but nothing forces it positive, and `λ < 0` means
`ω = ±i√|λ|` with one member of the pair growing. A first version of this round
claimed "real spectrum for every coupling, so a conserving throat cannot ring
up"; that is **false**, and two of its own three examples are the counterexample:

| `(α₁, α₂, β)` | `σ` | `λ` | growing? |
| --- | ---: | ---: | :---: |
| `(0.05, 0.05, 0.03)` | — | — | no |
| `(0.2, −0.13, 0.15+0.07i)` | **`2.470532`** | `−6.104` | **yes** |
| `(−0.4, 0.07, −0.09+0.31i)` | **`7.090982`** | `−50.28` | **yes** |

They were missed because the earlier root search seeded only `Re ω ∈ [1.1, 6.9]`
and discarded roots leaving that window — a search that by construction could not
find a root on the imaginary axis.

**What replaces it is a stability region with a closed form.** Both channel
functions fall monotonically along the imaginary axis from their `λ = 0` values,
so

```
stable  ⟺  α + β ≥ g₀ + G₀ = +0.02308202   and   α − β ≥ g₀ − G₀ = −0.07374262
g₀ = −1/(4π²),   G₀ = (π−d)/(4π² sin d)
```

verified against a negative-`λ` scan at all **221** grid points with **0**
mismatches — and only **56** of them stable, so positivity is a real restriction.

**Scope: `det(C − BΓ) = 0` is the rank-two mouth-active sector**, not the
spectrum. Level `n` has degeneracy `(n+1)²` and only two combinations can move,
so **23 of 25** modes at level 4 never leave the free eigenvalue. Inside the
sector there is also a mode *below* the free ground state (`λ = 0.311`) that an
`ω`-scan starting above 1 cannot see, then two per interlacing gap. And the
convenient claim that both channels run `−∞ → +∞` across every gap is false: the
`n = 0` constant mode is equal at both mouths, so the antisymmetric channel's
pole at `ω = 1` cancels and a first-gap root is conditional on `α − β`.

**#255's relation embeds exactly** as `B = [[0,0],[gain,0]]`, `C = I`, giving
`det(C − BΓ) = 1 − gain·G_d` — its own `1 − L`, matched to `3.5e-18`. Maximal,
but `BC† = B` is not Hermitian, and no finite Hermitian `A` reproduces it. This
is a *classification* of that boundary condition, **not** a diagnosis of its
off-axis poles: a self-adjoint throat can be unstable too.

Still put in: the boundary data — four real numbers chosen, not derived. The
throat is **point-supported**, so no interior, no proper length and **no delay**:
the `Δ` that carried #251–#255 does not survive into a point extension.

```bash
python -m experiments.closure_ledger.throat_operator_probe
# Verdict: SELF_ADJOINTNESS_IS_CONSERVATION_NOT_STABILITY  (8/8)

python scripts/geometrodynamics_v56_throat_operator.py --still v56.png
```

Full write-up: `docs/throat_operator.md`.

## The positive sector is a light cone

PR #256 showed that flux conservation does not imply stability, and mapped the
stable region on a two-parameter slice by scanning. The full four-parameter
answer is one inequality:

```
non-negative   ⟺   A ⪰ Γ(0)          (Löwner order)

Γ(0) = [[g₀, G₀], [G₀, g₀]] ,  g₀ = −1/(4π²) ,  G₀ = (π−d)/(4π² sin d)
```

— for **distinct non-antipodal mouths in the finite-`A` self-adjoint chart**.

**Why**, in one line: `dΓ/dλ ≻ 0` below threshold, so every eigenvalue of
`A − Γ(λ)` is strictly decreasing in `λ` while both run to `+∞` as `λ → −∞`; one
crosses zero below threshold **iff** it is already negative at it. Checked
against an actual negative-`λ` root scan on 200 random Hermitian `A` — all with
complex `β` and unequal mouths — **0 mismatches**, 19 stable and 181 not.

**And the monotonicity is a theorem, not a sample.** `dΓ_ij/dλ =
⟨δ_i, (H₀−λ)⁻² δ_j⟩` is a **Gram matrix** — PSD for free, positive definite
whenever the two mouths are distinct. Rebuilt mode by mode from the `S³`
addition theorem and agreeing with the closed form to **`8.1e-12`**, antipode
included.

**And the same argument counts.** For any `λ*` below the free ground state,
`#{eigenvalues < λ*} = #{negative eigenvalues of A − Γ(λ*)}` — a Krein-type
inertia theorem, **0 mismatches in 160 tests** at `λ* = −2, 0, 0.5, 0.9`.

**The geometry is a forward light cone.** Hermitian `2×2` matrices are `ℝ⁴` under
`A − Γ(0) = x₀I + x·σ`, and PSD is `x₀ ≥ |x|`:

| | |
| --- | --- |
| apex `A = Γ(0)` | a *doubly* degenerate zero mode |
| null boundary `x₀ = \|x\| > 0` | exactly one zero mode, `λ = 0` in the spectrum |
| interior | strictly positive |

with `x₀ = (α₁+α₂)/2 − g₀`, `x₁ = Re β − G₀`, `x₂ = −Im β`, `x₃ = (α₁−α₂)/2`.
Tested as a cone: convex, and closed under positive scaling *from the apex*.

**The boundary is detectable.** On it the secular function vanishes to `1.8e-17`
and the marginal mode is found by independent root-finding at `1.4e-14`. Outside,
the instability turns on continuously — `λ` linear in the distance `ε`, so
`σ = √|λ|` rises with exponent **`0.50001`**, coefficient predicted from the
eigenvalue slope (`−7.37443`) rather than fitted (`−7.37448`).

**#256's wedge is the `x₂ = x₃ = 0` slice** — exact on all 143 sampled slice
points, and **wrong on 65 of 400** general draws when reused by averaging the
mouths and dropping `Im β`. Those are precisely the two dimensions it cannot see.

**Where the apex sits.** `tr Γ(0) = 2g₀ = −1/(2π²)` at *every* mouth separation;
its eigenvalues are exactly #256's two channel thresholds; and `det Γ(0) < 0`
for `0 < d < π`, so `Γ(0)` is indefinite there — **`A = 0` is unstable wherever
the mouths are actually apart**, which no placement short of the antipode fixes.

**The exact antipode is a different statement.** `G_d` has a *removable*
singularity at `d = π` — not a pole — with `G_π(0) = +1/(4π²) = −g₀`. So
`Γ(0) = g₀[[1,−1],[−1,1]]` has eigenvalues `(2g₀, 0)`: **negative
semidefinite**, and `A = 0` is **marginally non-negative** there, sitting on the
cone's boundary with a zero mode in the symmetric channel. For a through-throat
geodesic on `S³` that is the natural configuration, so it gets its own test.

**And `A ⪰ Γ(0)` is the criterion in a chart.** `φ^reg = A q` needs `B`
invertible; the strata it misses are Dirichlet directions, reached only as
`‖A‖ → ∞`. The general criterion is `A_eff ⪰ P†Γ(0)P` on the allowed-charge
subspace — agreeing with the cone on 60 chart draws and with a root scan on 60
`k = 1` stratum draws, **0 mismatches** either way.

Within the stated box `|α_j|, |Re β|, |Im β| ≤ 0.2`, the stable fraction is
**`0.083`** — a genuine restriction, not a formality.

Still put in: the boundary data itself. `A` is four real numbers chosen, not
derived; which point *inside* the cone a physical throat corresponds to is
exactly as open as before.

```bash
python -m experiments.closure_ledger.throat_positivity_probe
# Verdict: THE_POSITIVE_SECTOR_IS_A_LIGHT_CONE_AT_GAMMA_ZERO  (10/10)

python scripts/geometrodynamics_v57_throat_positivity.py --still v57.png
```

Full write-up: `docs/throat_positivity.md`.

## Static two-source throat tomography

**Not the roadmap's two-wave invariant** — that step stays open. The object here
is a *static* source-interaction kernel at a fixed spectral parameter: it carries
no local null momenta, so it cannot distinguish equal-energy collinear from
counterpropagating waves, which was the whole control behind
`𝒞 = I_A I_B (k_A·k_B)²`. The index `(i,j)` labels **mouth channels**, not the
geodesic/winding branches of #253–#255. What it *is* is an exact static inverse
result.

PR #253 ended rank counting by naming what it could not supply: a quantity that
**vanishes** when a source is removed rather than merely becoming
underdetermined. Superposition makes every linear functional additive, so the
object has to be quadratic, and its cross term is

```
𝒞(y_A, y_B) = G(y_A,y_B) + Re Σ_ij G(y_A,c_i) R_ij G(c_j,y_B) ,  R = (A − Γ(λ))⁻¹
```

Computed from a functional that carries its own self-energy terms:
`Q[a,b] − Q[a,0] − Q[0,b]` matches `ab·𝒞` to `2.8e-17`, and removing a source
means evaluating the same functional at `b = 0` — not multiplying by zero.

**The throat channel is rank two at any source count** — `Vᵀ S V` with `V` of
shape `2 × N`, rank **2** against a direct table of rank **12** for 12 sources.
Off the chart `rank R = rank B`, but static sources see only `Re R`, whose rank
is two even for a *complex* rank-one boundary condition and one for a real one.

**Three things that look like the signature and are not.** The cross term being
nonzero is interference. **Anisotropy** — the interaction depending on more than
the geodesic separation, which no free field here can do (`8e-17`) — is real at
`66%` of the mean, and two **disconnected** scatterers give `69%`. And the
**off-diagonal response block** is nonzero for `β = 0` too, because `Γ` couples
the mouths through the ambient field: it is a *cross-mouth* channel, not
"through the throat".

**What discriminates is a parameter count.** The static invariant determines
three numbers, the entries of `S = Re R`; two independent scatterers have two
knobs, so their image is a surface with the exact equation `S₁₂ = G₀ det S`
(`1.4e-16` on 200 draws). The **disconnection defect** `𝒲 = S₁₂/det S − G₀` is
zero on it, and on real `β` equals **`−β`** to `5.0e-16` — independent of the
self-energies, the separation, and the **Löwner margin**: driven from margin
`0.4` to `0.004` the invariant grows `3.8×` and `𝒲` drifts `2e-17`, which
answers #255's caution that a resummed field measures the pole rather than the
source. Stated exactly, `𝒲` detects **off-diagonal mouth-boundary mixing
relative to the diagonal two-scatterer null model** — not topology, not an
interior.

**And it is a protocol.** An observer who measures interaction energies and knows
the background and the mouth positions, but is not told the boundary data,
recovers `S` by least squares and gets `𝒲` to `1.1e-16` from 24 placements.

**Which field is being solved decides the blind spot.** `𝒲 = 0` off `β = 0`
needs complex `β` — and a **real** scalar (what #254 solves) requires the
self-adjoint domain to be conjugation-invariant, `A = A*`, hence `β` real.
Measured: with complex `β` a real unit static source produces a field with
`Im = −2.4e-3`. So for the arc's field content **there is no blind family**;
it belongs to a deliberately time-reversal-breaking complex extension. Inside
that extension #257's gate removes the `Re β > G_d` branch, the surviving
couplings are *smaller* than their own self-energies (`0.215–0.254` against
`0.25–0.40`), and even there the limit is the **protocol**: phase-sensitive
complex sources give the full complex `R` and hence `A = Γ + R⁻¹` at **one**
spectral parameter (`3.9e-15`). The real-static-source reconstruction needs two
spectral parameters — both positive and below `λ = 1`, since `λ = ω²` makes a
negative `λ` an imaginary frequency — and returns `A` to `1.1e-16`.

**The antipodal endpoint, on its own.** At `d = π`, `Γ(0)` is negative
semidefinite, so the static response is singular as `A → 0` and the invariant
**diverges** like `1/ε` — while `𝒲` stays exactly zero through four decades.
**Size is not evidence.**

Everything is evaluated at `A = (0.30, 0.35, β = 0.06)` with `d = 1.3` —
strictly inside #257's cone, Löwner margin `0.323`.

```bash
python -m experiments.closure_ledger.two_source_probe
# Verdict: STATIC_THROAT_TOMOGRAPHY_MEASURES_THE_MOUTH_MIXING  (12/12)

python scripts/geometrodynamics_v58_two_source.py --still v58.png
```

Full write-up: `docs/static_throat_tomography.md`.

## The two-wave invariant is branch-resolved

**Roadmap step 3, properly.** #258 built a static kernel and said plainly it was
not this: no local null momenta, so it could not tell equal-energy collinear from
counterpropagating waves. This round solves the time-dependent field on the
throated ESU, builds the improved conformal stress tensors, and applies exactly
that control.

The known WKB result is the **control, not the result**. The research content is
the difference between the exact multipath throat-coupled field and that limit.

**What is solved.** The retarded field of a pulsed point source, exactly, by
Krein's resolvent formula in the frequency domain, inverted along the retarded
contour `Im ω = ε` — which is exact, since `φ(t) = e^{εt}(1/2π)∫du e^{−iut}
φ̂(u+iε)`. Derivatives are analytic, not differenced: the four-gradient and
Hessian close in form from `∇χ = −(y−(x·y)x)/sin χ` and
`∇∇χ = cot χ (δ_ab − n̂_an̂_b)`.

| | |
| --- | ---: |
| free field vs #254's closed-form image sum | **`3.3e-16`** |
| conformal wave equation, with and without the throat | `4e-16` relative |
| trace of the improved stress tensor | **`1.9e-15`** relative |

The trace is a real test: `□φ` comes from the solve rather than its on-shell
value, so `T^μ_μ = φ(□φ − φ)` instead of vanishing algebraically.

**The known limit, as a limit.** With the arriving directions exactly parallel
or antiparallel by construction (`1e-12`), the pointwise
`𝒩 = (T_A:T_B)/(T_A⁰⁰T_B⁰⁰)` converges to WKB's `(1 − n̂_A·n̂_B)²`: head-on
`3.99995`, collinear `1.8e-10` at `ω₀ = 48`. The collinear null is *stronger*
than leading order — the two wavefronts share their normal exactly, so amplitude
gradients cannot tilt either `k`. Convergence is measured too: with `ε` at the
frequency spacing the answer is four orders wrong and looks plausible.

**The result — multipath destroys the collinear null.** Sources fixed,
observation point fixed, only the branch changes:

| branch pair | `𝒩` exact | geometry |
| --- | ---: | ---: |
| `A` direct + `B` direct | **`1.905e-07`** | `0` |
| `A` **long-way winding image** + `B` direct | **`3.99806`** | `4` |
| `A` direct + `B` **via a mouth** | **`0.56501`** | `0.56367` |

The winding image runs the other way round the sphere, so its arrival direction
reverses and a collinear pair reads head-on. The free control at the same instant
has *no* second arrival (energy product `4e-29` vs `1.2e-02`), so the mouths
**create** the branch rather than bending one. **The collinear null is not
spoiled by curvature corrections — those are `1e-7` here — but by multipath, at
`O(1)`.**

**The `(i,j)` audit, and the control that scopes it.** Writing `j` for the mouth
the source drives and `i` for the mouth the signal leaves from, all four two-leg
paths are enumerated rather than minimised over. The prediction depends on `i`
alone, so the four carry **two** values (`0.563669`, `0.651935`) and the field
has to pick, not merely match a number — it does, to `8.1e-04`.

Then the control PR #258's review taught this arc to run first: the same
channels with **`β = 0`**, the mouths *disconnected*. Swept over `β ∈ [0, 0.26]`,
all inside #257's cone, `𝒩` moves by `6.2e-07` — `7e-06` of the `0.0883` that
separates the two exit mouths — while the channel's weight moves `0.6%`. It has
to: `𝒩` is amplitude-normalized, and a single channel is a single arrival
direction. **This observable sees structure at the mouths, not the connection
between them.** The multipath result stands; the throat's *non-locality* is not
what supplies it. What sees the connection is `𝒲 = −β`, below.

**`ΔT_{μν}` disagrees with `T_A:T_B` completely.** The bilinear cross term
`ΔT = T[φ_A+φ_B] − T[φ_A] − T[φ_B]`, built from three evaluations of the same
functional, is traceless to `1.8e-15` and vanishes **exactly** when either
source is switched off:

| configuration | `T_A:T_B/(T_A⁰⁰T_B⁰⁰)` | `ΔT⁰⁰/√(T_A⁰⁰T_B⁰⁰)` |
| --- | ---: | ---: |
| collinear | `1.9e-07` | **`2.000`** |
| head-on | `3.998` | `1.044` |

The interference energy hits its **maximum possible value, 2** — two parallel
waves adding coherently — exactly where the invariant is null. **A backreaction
estimate driven by `𝒞 = T_A:T_B` would look at the collinear case, see nothing,
and be wrong about its own source by the size of the whole effect.**

**The other corrections, quantified.** Free arrivals land on #253's ledger to
`1.3e-03` with Maslov signs `+ − +`; the throat adds two-leg arrivals, checked at
the causal onset because `R(ω)` has poles and a throat arrival rings up. **Tail:**
`S³ × R` is conformally flat, so Huygens is exact — between arrivals the free
field is `1.4e-08` of its peak against `8.1e-02` with the throat, a factor of
`5.7e+06`. Every tail here is the throat's. **Caustic:** WKB's `1/(4π sin χ)`
diverges at the antipode where the exact kernel is finite and *linear in* `ω`;
the exact/WKB ratio is `|sin(ωe)|`, a function of `ωe` alone, collapsing across
three carriers to `6.6e-15`, so the caustic is cut off at `e* ∼ 1/ω`.

**And it closes back on #258.** `∫dt φ = φ̂(0)`, so the DC content of the solved
time series *is* the static kernel that round did its tomography on. Running the
protocol on numbers from the dynamic solver returns `𝒲 = −0.060010` against
`−β = −0.06`, with the `O(ε²)` contour bias Richardson-extrapolated and both
numbers reported.

**No backreaction:** the stress tensor is computed from the field and never fed
back. That is the next step, and it now has a concrete object to feed.

```bash
python -m experiments.closure_ledger.two_wave_probe
# Verdict: THE_TWO_WAVE_INVARIANT_IS_BRANCH_RESOLVED  (12/12)

python scripts/geometrodynamics_v59_two_wave.py --still v59.png
```

Full write-up: `docs/two_wave_invariant.md`.

## The throat has an interior, and the interior is the delay

**The flux-conserving throat operator, finally.** Every round from #253 to #259
carried the same disclaimer — *point-supported, no interior, no proper length,
no delay* — and what stood in for one was a rank-one **mouth-transfer** model:
field values only, no normal-derivative matching, no reflected channel, `1×1`
where a conserving junction needs `2×2`, and lossy for `κ < 1`.

**Two lines are put in.** A tube of length `L`, cross-section `𝒜` and interior
mass `m` joins the mouths; its Dirichlet-to-Neumann map is exact,

```
N(λ) = 𝒜k [[cot kL, −csc kL], [−csc kL, cot kL]],    k² = λ − m²
```

and the matching is value and flux continuity, `q = −NΦ`. Everything follows.
Since `det N = −(𝒜k)²`, the chart is closed-form: `A(λ) = −N(λ)⁻¹`, so the
transmission amplitude is `β(λ) = csc(kL)/(𝒜k)` and the self-energy is
`α(λ) = cot(kL)/(𝒜k)`. **The boundary condition is now frequency-dependent, and
that dependence *is* the interior.**

**Where the self-adjointness lives.** The conservative object is the **enlarged
system**, ambient `⊕` tube, with the `λ`-*independent* matching above — one
self-adjoint operator on `L²(S³) ⊕ L²([0,L])`. Eliminating the tube leaves a
`λ`-*dependent* boundary condition: `A(λ)` is the **Weyl (`M`-) function** of
that elimination, not itself a self-adjoint operator on the ambient space, which
an energy-dependent boundary condition never is. What it is, is a matrix
**Nevanlinna** function — monotone in `λ` between its poles — and that
monotonicity is the enlarged system's self-adjointness showing through. What is
checked pointwise is that the elimination is *faithful*:

| | |
| --- | ---: |
| `‖BC† − CB†‖`, seven `λ` on both sides of zero | **`0.0`** |
| DtN vs its own interior, **sesquilinear** Green's identity | `1.5e-07` |
| **control** — #255's rank-one transfer model | **`0.30`** |

**The result — the throat transmits at the traversal time.** The measured object
is the **two-mouth block's** impulse response: the source and observer legs are
gone, but `Γ` — the ambient's own mouth-to-mouth propagator — stays in, so this
is the coupled ambient+tube response and not the throat alone. `r₁₁` (same mouth
in and out) starts at **`t = 0`** — a wave that reaches a mouth is partly
reflected instantaneously — and `r₁₂` (opposite mouths) starts at
**`min(L, d)`**, with `d(onset)/dL = **1.0071**` against a predicted `1` below
the ambient path and a spread of `0.0` above it. The ambient *also* connects the
two mouths, along a geodesic of length `d`, whether or not they are joined:
#258's cross-mouth channel and #259's `β = 0` control, now **separated in time**
instead of by rank counting — and the reason the onset saturates. A frozen `A`
transmits at `0.0`, which is what a point throat is.

The ledger is a derivation, not a fit: on the contour `cot x = −i − 2iΣe^{2ikx}`
and `csc x = −2iΣe^{i(2k+1)x}` to `4.5e-16` and `1.7e-15`, so the same-mouth
entry carries `0, 2L, 4L…` and the cross-mouth entry `L, 3L, 5L…`. **The
parities are the physics**, and the reflected channel is the one the rank-one
model does not have at all.

**There *is* a point limit — and it is not a finite `A`.** Freezing `A` at
`A(λ₀)` is exact at `λ₀` and `4.3%`, `17%`, `73%`, `121%` wrong at
`1.05, 1.2, 2, 3 λ₀`. Everything in `A` varies through `kL`, so the range over
which freezing it is defensible is an `O(1/L)` **frequency** scale.
As `L → 0` the antisymmetric channel converges to `−L/(2𝒜)` while the symmetric
one **diverges** like `2/(𝒜λL)`. A first draft concluded from that that the
limit does not exist. It does: a boundary pair is defined up to
`(B, C) → (MB, MC)`, so a diverging *chart matrix* means the limit has **left
the chart**. Row-scaled, the pair converges linearly in `L` (rate `𝒜λ/2 = 2π`,
measured) to

```
(B, C)  ⟶  (P_anti, −P_sym)  ,   i.e.   Φ_anti = 0   and   q_sym = 0
```

a **mixed Dirichlet–Neumann** stratum: maximal (`rank[B|C] = 2` throughout),
self-adjoint, and reached by **no finite Hermitian `A`** since both blocks are
singular — exactly the stratum #257's review said the chart does not cover. So
the correct statement is *"no finite-`A` point limit"*. It is also what a very
short pipe should do: short the two mouths together and store nothing.

**And that zero mode breaks #258's tomography.** At `λ = 0` the static response
collapses onto `[[1,−1],[−1,1]]` to `4e-05`, `det S → 0` **linearly** in `λ`
(coefficient `149.08`, constant to `1e-3` over four decades), so
`𝒲 = S₁₂/det S − G₀` diverges like `1/λ`. **What that falsifies is the generic
finite-`A` family** — every member of which has `rank S = 2` — and *not*
point-ness: the tube's own short-tube stratum gives `R → diag(0, −1/Γ_anti)`,
rank one as well, and the tube converges to it. (The first draft claimed the
stronger and wrong version.) An interior mass restores the rank
(`det S ∝ −148.7 m²`), and off the collapse **`𝒲 = −β(λ)` exactly**, to
`3.1e-13`: #258's theorem survives the generalization and returns the interior's
own amplitude.

**And the candidate fails the stability gate.** `A(λ)` decreases and `Γ(λ)`
increases (#257's Gram identity), so `A − Γ` is strictly monotone between poles
and each channel has **at most one root** — a count, not a scan. The symmetric
channel always has exactly one, and it is at `λ < 0`: **an exponentially growing
mode, for every choice of parameters.** In the `σL, σd ≫ 1` limit its rate
matches **`σ* = 2√(π/𝒜)`** to `1.5e-03` with **no `L` in it**, two mouth
separations agree to `3.9e-09`, and the channel splitting is `1.04·e^{−σ*d}` —
so there the mode **localizes to a single mouth**. It is generated at the
**point-mouth/tube interface**.

The working throat is *not* in that limit, and the qualification matters: at
`𝒜 = 4π` the asymptotic form gives `σ* = 1` while `L = 0.9` gives `1.417`, and
`σ*` runs `1.769 → 1.152` across `L = 0.4 → 3`, a spread of `54%`. An earlier
draft said the mode "belongs to the mouth and not the interior"; the interface
statement plus asymptotic localization is what the data support, and it is what
makes a finite-radius mouth the right discriminator.

**That is the round's closure result, and it gates the roadmap.** An action or a
backreaction computed on a background with a growing mode inherits the mode, so
the next construction is not #261 but a **finite-radius mouth or neck** — the
ambient solved outside two small balls rather than a point interaction with a
radius parameter — with one question to answer: *does the negative mode survive?*

**And the contour rule has two parts, not one.** `ε > σ*` is the analytic
Bromwich condition; at `0.03` *below* `σ*` the inversion returns a field with
support before its own light cone, a pedestal at 99% of the peak. But at a
clearance of `+0.02` the contour is *above* the mode and the pedestal is still
`2.6e-03` — that clearance is `0.95` of the frequency spacing `2π/span`, so the
pole is cleared but **unresolved by the grid**, #259's lesson a second time. Both
`ε > σ*` and `ε − σ* ≫ 2π/span` are needed; at `14–72` spacings the pedestal is
`1.0e-16` and the recovered onset converges to `0.0092`, four time steps. Neither
condition stabilizes anything: above `σ*` the inversion returns the correct
retarded solution *of an unstable system*. Whether a finite-radius mouth or neck
geometry removes the mode is open, and should be settled **before**
stationary-action or backreaction work.

**Which frequencies.** The delay and the bounce ledger are statements about the
model's **analytic structure at all frequencies** — a causal onset is a UV
object, and the probe pulse carries content to `ω ∼ 30`, far above `σ* ∼ 1.4` —
so they are exact results *about this model*, not predictions about a resolved
physical mouth. The static and low-frequency results sit inside the band. And
`𝒜` is a **one-dimensional coupling**, not an area with a radius attached.

```bash
python -m experiments.closure_ledger.finite_throat_probe
# Verdict: THE_INTERIOR_GIVES_A_DELAY_AND_THE_POINT_MOUTH_IS_UNSTABLE  (10/10)
#   (the mode is interface-generated; mouth-localized only for σL, σd ≫ 1)

python scripts/geometrodynamics_v60_finite_throat.py --still v60.png
```

Full write-up: `docs/finite_conservative_throat.md`.

## The negative mode does not survive a finite mouth

**PR #261 — the gate #260 set, answered.** That round found its conservative
throat carried an exponentially growing mode for *every* choice of parameters and
stopped the roadmap on one question: **does it survive a finite-radius mouth?**

**No — and the statement is structural, not parametric.**

**What a resolved mouth changes.** A point interaction must subtract the
`1/(4πχ)` divergence and keeps the **renormalized** self-energy `g(λ)`, which is
*negative* — a leftover of an infinite subtraction, and what the mode fed on. A
sphere needs no subtraction. Smearing the coupling over `∂B_a` — the same
operator on both sides, so the composite stays manifestly self-adjoint — gives

```
𝒢_self(λ) = f(a,λ)·G(a,λ)        𝒢_cross(λ) = f(a,λ)²·G(d,λ)
```

with `G` the **unsubtracted** Green function and `f(χ,λ) = sin(ωχ)/(ω sin χ)` the
regular radial solution. Both are **mean-value identities**, checked against
direct quadrature on `S³`: the cross one to **`1.0e-10`**, the self one to
`4.1e-04` (grid-limited by the singularity at coincidence, and reported as such).

**The signs decide it.** At `λ = −σ²` the tube gives `−coth(κL/2)/(𝒜κ)` and
`−tanh(κL/2)/(𝒜κ)` — strictly **negative**, a passive interior has no restoring
force — while the ambient gives `f·G(a) ± f²·G(d)`, strictly **positive**, every
bracket in `(0,1]` once `a < d/2`, which disjoint mouths require anyway. A
difference of a negative and a positive number has no zero. **3078 samples** over
`(a, d, L, 𝒜, m, σ)`: **0 roots**, worst approach `−5.1e-04`.

**And #260's mode was the linearization.** That round froze the mouth at the
*constant* `1/(4πa)` — the leading term of `G(a,λ) = 1/(4πa) + g(λ) + O(a)`. The
exact `G(a,−κ²)` is **screened**, `≈ e^{−κa}/(4πa)`; the constant is not, so it
eventually beats the tube's `−1/(𝒜κ)` and crosses:

| `a` | linearized root `κ*` | `κ*·a` | exact |
| ---: | ---: | ---: | :---: |
| `0.02` | `50.02` | **`1.0004`** | none |
| `0.05` | `20.05` | **`1.0025`** | none |
| `0.15` | `6.814` | **`1.0221`** | none |
| `0.35` | `3.220` | **`1.1269`** | none |

**The root sits at `κa ≈ 1` — the edge of its own approximation.** The two models
agree to `0.8%` for `κa ≤ 0.1` and differ by `1000%` at `κa = 3`, disagreeing not
in magnitude but in **sign**. #260 suspected exactly this and could only record
the suspicion; this is the demonstration.

**Where the mode went: soft, and positive.** One state below the free gap
`λ = 1`, in the symmetric channel, with

```
λ₀  ⟶  8πa/(𝒜L)
```

two mouth capacitances `4πa` restoring a tube of volume `𝒜L` — ratio `0.998` at
`a = 0.005`. **The point limit drives the mode to zero from *above*.** #260 did
not get a rate slightly wrong; it took a mode approaching `0⁺` like `a` and put
it on the other side of zero, at `λ ≈ −1/a²`.

**The good results survive.** The traversal delay keeps slope **`1.0010`** in `L`
and saturates at the ambient path to `0.0`, the mouth adding only a sub-leading
`O(a)` shift (a first draft predicted `−2a` from an ambient block missing the
shell form factor; the measured slope is quoted and the prediction recorded as
wrong). The static response is still rank one and #258's `𝒲 = −β(λ)` still holds
to `3.6e-12` — all of which came from the *tube's* zero mode, which the mouth
does not touch. The contour is easier too: `ε = 0.4` where #260 needed `ε > 2`.

**What is still put in.** One channel per mouth, so only `ℓ = 0` couples; the
dropped multipoles obey #250's screening law, dipole/monopole `= 0.934·(a/d)`
across a decade in `a`, dropped power `6.9e-05` at the working radius. The mouths
are **spheres in a fixed ambient, not a solved neck**. No backreaction.

**This ungates the roadmap.** #260 blocked stationary action and backreaction
because an integral over a field on a growing background measures the mode. That
reason is gone.

*(#262 removes the balls and finds the same answer as a theorem. It also
corrects one claim here: "one state below the gap" holds for `L < π`, not
structurally — above that the tube's own harmonics enter the gap and each brings
another.)*

```bash
python -m experiments.closure_ledger.finite_mouth_probe
# Verdict: THE_NEGATIVE_MODE_DOES_NOT_SURVIVE_A_FINITE_MOUTH  (9/9)

python scripts/geometrodynamics_v61_finite_mouth.py --still v61.png
```

Full write-up: `docs/finite_radius_mouth.md`.

## The balls removed, and the answer made a theorem

**PR #262 — the one limitation #261 named about itself, closed.** That round
answered the gate, and said plainly where it was weakest: its mouths were
**spheres in a fixed ambient**. The balls were never removed — it smeared the
coupling over `∂B_a` while still using the *whole sphere's* Green function — and
only `ℓ = 0` coupled. Two things could have hidden there: a self-energy wrong
because the ball is still in, and a multipole that goes soft where the monopole
does not.

This round removes them: `Ω = S³ ∖ (B_a(c₁) ∪ B_a(c₂))`, tube glued along the
boundary spheres. **The answer is still no — and it is now a theorem.**

**The theorem.** With the balls removed there is **no subtraction anywhere**, so

```
E[φ,u] = ∫_Ω (|∇φ|² + φ²) dV  +  𝒜 ∫₀^L (|u'|² + m²|u|²) ds
```

is a sum of non-negative terms. `E = 0` forces `φ ≡ 0` on `Ω`; matching then
gives `u(0) = u(L) = 0`, and Poincaré gives `𝒜∫|u'|² ≥ (π/L)²𝒜∫|u|²`. Hence
**`λ > 0` for every configuration — all multipoles, no truncation, no sweep.**
That is a change of footing rather than a refinement: #261 established a *sign*
on a reduced `2×2`; this is positivity of the form itself, and the renormalized
`g(λ) < 0` that #260's mode fed on has nowhere to enter.

**The object it is about, checked.** The exterior DtN `N_ℓ(λ) = −4π sin²a·ψ'/ψ`
comes from shooting `v'' + [λ − ℓ(ℓ+2)/sin²χ]v = 0` from the far pole, and
matches an independent `ℓ=0` closed form to **`1.7e-14`**. It is positive in
every channel and **increasing in `ℓ`** — the monopole is the softest, so the
higher channels cannot be the first to go soft — and `N₀ → 4πa`, the capacitance
of a sphere. Explicit trial configurations give Rayleigh quotients from `0.359`
up, all above the computed lowest mode `0.1075`, and the degenerate purely
interior case lands on the Poincaré floor `π²/L²` to `2e-07`. A **1197-sample**
sweep agrees: `0` roots, worst approach `−1.6e-03`, reaching #261's conclusion
from a different construction (`1/N₀` on the diagonal, not `f(a)G(a)`).

**And the monopole truncation was never a stability limitation.** A one-channel
tube presents a single number at each mouth, so it drives `ℓ = 0` and nothing
else; the `ℓ ≥ 1` sectors are the exterior's own modes, `1.45×` stiffer or more.
#250's `(a/d)^ℓ` screening bounds missed **amplitude**, not the answer.

**What the fixed ambient cost, priced.** `f(a)G(a)` against `1/N₀`: `1.3e-04` at
`a = 0.02, λ = 0`, `3.8e-03` at the working radius, `11%` at `a = 0.35, λ = −4`
— the fraction of the sphere wrongly left in. A limitation with a measured size
is a different object from one with a caveat. The same measurement prices the
approximation this round has *not* removed: the reduced `2×2`'s cross term is
**single-scattering**, its neglected series expands in `cross/self = 0.8·(a/d)`
= `9.5e-04` at the working point and at worst `0.16` anywhere sampled — too
small to flip a sign, and irrelevant to the theorem, which does not go through
the reduced model.

**The soft mode is forced, not found.** The same style of argument that kills the
growing mode produces this one, from the two ends of the gap: `F_sym → +∞` as
`λ → 0⁺` (the tube's `2/(𝒜λL)`), and `→ −∞` as `λ → 1⁻`, because the exterior
stiffness **vanishes** at the free ESU threshold —

```
N₀(λ)  ⟶  2π (π − a + sin a cos a) · (1 − λ)
```

reproduced to `3.1e-05`, first order exactly (error `×0.1` per decade). A
continuous function running `+∞ → −∞` has a zero.

**One correction to #261.** Its "exactly one state below the gap" is a statement
about `L < π`, not a structural one. The channel functions have **poles** at
`λ = (2πj/L)²` and `(π(2j−1)/L)²`; above `L = π` these enter the gap and each
brings a genuine extra state just above it — three at `L = 8`. A pole is a *sign
change with no zero*, so crossing-counting alone reports states that are not
there; the residual separates roots from poles by fifteen orders of magnitude
(`1e-15` against `1e+15`), so the discrimination is not a tuned threshold. None
of it touches the stability answer — every one of those states is positive, as
the form requires.

**What is still put in.** The tube has **one transverse channel**, so `𝒜` is a
coupling and the neck is a quantum-graph edge, not a solved cross-section. The
ambient metric is **fixed**. **No backreaction.** Those are the next round's
subject, not gates on it.

```bash
python -m experiments.closure_ledger.neck_probe
# Verdict: THE_NEGATIVE_MODE_DOES_NOT_SURVIVE_THE_NECK  (9/9)

python scripts/geometrodynamics_v62_neck.py --still v62.png
```

Full write-up: `docs/finite_radius_neck.md`.

## A + B does what rescaling A or B cannot

**PR #263 — the gate #260 set is closed, so this is the roadmap's first GR
question.** And deliberately *not* "does spacetime pinch off?":

> does `A + B` produce a metric response that rescaling `A` or `B` alone cannot?

**Yes.** At the working point **`0.921`** of the interference response lies
outside everything rescaling can reach, and the interference term is
**comparable in size** to the single-wave responses (`‖β_ΔT‖/‖β_A‖ = 1.02`)
rather than a correction to them.

**Why it is a linear-algebra question.** The field equation is linear so the
fields add; `T` is quadratic so `T[A+B] = T[A] + T[B] + ΔT` with `ΔT` bilinear;
linearized Einstein is linear so the responses add. Rescaling `A → cA` sends
`β_A → c²β_A`, so **everything reachable is the two-parameter family
`{c²β_A + d²β_B}`** and the question is a projection residual. Measured off the
full linear *span*, which strictly contains that cone — so the figure is
conservative.

**The channel is forced, not chosen.** The ESU is held static by a fluid this
arc never specifies. A perfect fluid carries no anisotropic stress — in an
orthonormal frame `T_ab = diag(ρ,p,p,p)` whatever the anisotropy — so neither it
nor `Λ` touches the traceless spatial part. The **transverse-traceless** sector
is the only channel whose answer does not depend on what was never put in.

**The response, derived.** Cartan about the ESU in the homogeneous anisotropy:

```
δG^TT_ij = β̈_ij + (8/a²) β_ij        ⟹    ω₃ = 2√2,   ω₃² > 0
```

so the tensor sector is **stable** — #260's gate, applied to this round's own
background. The connection comes from *solving* the torsion-free condition after
a first draft's remembered formula produced a `G₀₀` containing `ä`. Three
validations pass (round `S³`; ESU, independently reproducing `two_wave`'s
`_EINSTEIN`; closed FRW), the first-order piece is automatically traceless, and
`δG_{0i} = 0` identically.

**It is not a universal constant.** `0.88–1.00` across time windows, `0.56–0.99`
across carriers, `0.929` with the throat and `0.986` without. Large everywhere,
constant nowhere — so the range is the headline.

**The control is the round's real content.** A first attempt reported `0.982`
unreachable and it was **pure quadrature noise**: independent rules for the same
quantity correlated at `−0.04`, and nothing about the run looked wrong. Two
uncorrelated noise vectors give a residual of `≈1`, so the failure mode and the
desired result are *the same number*. The cause was the singular set missing the
two **mouths** — and `two_source` puts the first at `(1,0,0,0)`, exactly the
natural quadrature axis, so the rule's own pole sat on a `1/χ⁴` divergence and
refining made it worse. With all eight singular points under a smooth partition
of unity, two refinement levels now agree to correlation `0.970`–`0.999` and
magnitude drift `0.027`, residual moving `0.0029` — under a **deterministic** Householder basis, after an `svd`-derived one made the whole rule platform-dependent.

**And the resonance test, done properly, reverses itself.** A first version
argued the channel was off resonance *by construction*: the conformal scalar on
the ESU has integer spectrum `ω_n = n+1`, so a source built from the free field
rings on integers, and `ω₃ = 2√2` is irrational — `0.172` from the nearest. All
true, **of the uncoupled ambient**. The throat rings where `det(A − Γ(ω))`
vanishes, and those zeros are *not* integers: `0.875, 1.854, 1.872, 2.706,
2.878, …`, the nearest sitting **`0.050`** from `ω₃`, with a local spacing of
only `0.17`. Across sixteen throat configurations the nearest is always within
`0.086`, and at `d = 0.9` it is **`0.001`**. So the corrected statement is a
working-point one pointing the other way — off resonance with the free ambient,
**generically near-resonant with the coupled system** — and that is the
mechanism the first version lacked for the carrier sensitivity.

**And what this channel cannot say.** A traceless shear preserves volume
*exactly* (`det e^β = 1`) and mouth area *to first order* — the measured areal
change is `0.403 ε²`, and **positive**, so what second-order effect exists opens
rather than pinches. **The `n = 3` channel therefore cannot answer whether the
interaction metric moves toward a neck or away from one**; it distorts the mouth
into an equal-area ellipse. That question needs the areal sector on a resolved
neck.

**What is still put in.** The `n = 3` harmonic only; a **fixed** ESU; **point**
sources and #257's **point** throat rather than #261/#262's resolved mouths; and
a **linear** response, not a solved coupled system.

```bash
python -m experiments.closure_ledger.backreaction_probe

python scripts/geometrodynamics_v63_backreaction.py --still v63.png
```

Full write-up: `docs/metric_backreaction.md`.

## The interference metric deforms toward a neck

**PR #264 — the geometric verdict #263 could not give.** #263 answered *does
`A + B` do something rescaling cannot* with a large yes, and then proved its own
channel could not answer the question actually asked: `δA/A = −½⟨h_nn⟩` vanishes
**identically** for any transverse-traceless field. So the question —

> does the interaction metric deform **toward** a neck, **away** from one, or
> merely **oscillate**?

— moves to the scalar sector, posed as an **initial-data constraint solve** and
not an evolution. On a maximal slice the `K` terms in the Hamiltonian constraint
are quadratic, so `δR⁽³⁾ = 16πG δρ` with no time derivatives in it: a constraint
has no sound speed and no Eddington mode, which is exactly what made #263 avoid
the scalar sector.

**Toward a neck — at the wide working throat, off resonance. Both mouths
close.**

```
ΔA/A = ( −2.06e−03 ,  −1.88e−03 )      in units of 2πG
```

Negative in all eight controls — two quadrature levels, two mouth radii, two
gluings.

**And the qualifier is load-bearing.** The working tube carries `𝒜 = 4π` against
a mouth sphere of area `4π sin²a` — wider than its own mouths by `400×` at
`a = 0.05`. Set them equal, so the tube is as narrow as the mouths it joins, and
`k = 1/sin a` puts the *same* length `0.9` at `kL/π = 5.73` and `2.87` — past
five poles and past two. **The sign does not survive:**

| `a` | `kL/π` | `𝒜 = 4π` | matched `𝒜 = 4π sin²a` |
|-----|--------|----------|------------------------|
| 0.05 | 5.732 | closes / closes | **opens / opens** |
| 0.10 | 2.870 | closes / closes | **closes / opens** |

So this is a statement about a throat, not about the interference source. **And the mechanism is not the obvious one:** the interference energy
*alone* would **open** the mouths (`U(c_j) > 0` at both). The throat's own
monopole layers overshoot that and invert it.

**The problem, and why it is hard.** With `g = ψ⁴ĝ`, `ψ = 1+u`, the constraint
is `∇²u + 3u = −2πG δρ` and `ΔA/A = 4u`. That operator is **exactly degenerate**
on `S³` — `4 = (n+1)²` at `n = 1`, the ESU's dipole level — so it has no
solution on the closed sphere. Removing the two mouth balls is what makes it
solvable, and the field is then a source term plus monopole and dipole layers on
the two mouth spheres plus a free kernel element: **twelve unknowns, twelve
equations.** The solvability condition `Σ_j A_j c_j + Σ_j D_j = S_σ` is an
*identity*, not a modelling choice, because `(2/π²) cos χ` is exactly the
projector onto the kernel.

**Two results that came out the other way from expectation.** The dipole layers
are **required** — two monopoles sweep only the plane of the two mouths, and the
monopole-only condition fails by **62.5%** of the obstruction, so without them
there is no solution at all. And then they **do not move the answer**: the
off-plane response is `6e−17`. `ΔA/A` is, to `0.09%`, a linear functional of the
obstruction alone — which is lucky, because the `ℓ = 1` source moments drift
`41%` between quadrature levels where the obstruction drifts `1.5%`.

**Why.** The tube's `ℓ = 0` constraint channel is `∂_s² + 4π/𝒜` — a **cavity**.
At `kL = nπ` the response has a pole and the sign flips; the scan finds flips at
`3.133` and `6.260` against `π` and `2π`. The working throat sits at `kL = 0.9`,
inside the first cell; the matched throat does not.

**One bug that check found.** The `ℓ = 1` rows were a `cosh`/`sinh` transfer
matrix, which costs a condition number of `e^{2κL}` — invisible at `𝒜 = 4π`
(`e^{1.8} = 6`), fatal at the matched area (`e^{36} = 4.4e+15`). The first
matched-tube run reported `cond = 2.9e+15` **and an answer anyway**. Carrying
the tube's end amplitudes as unknowns never forms `e^{+κL}`: `18×18`, every
coefficient bounded by one, `5.5e+07` at the matched area and `1.5e+04` (from
`2.1e+05`) at the wide one. The reference solves reproduce to the *same*
`4e−10`, so it is a change of form and not of content.

**What is still put in.** The ESU's fluid is held **rigid**; the source is #263's,
on a **fixed** background with **point** sources; the response is **linear** in
`G` and quadratic in the wave amplitude.

```bash
python -m experiments.closure_ledger.areal_probe

python scripts/geometrodynamics_v64_areal.py --still v64.png
```

Full write-up: `docs/signed_areal_response.md`.

## Which throat is physical — and the sign reverses

**PR #265 — the geometry stops being decoration.** #264 found that matching the
tube's area to its own mouths *flips* the sign of `ΔA/A`. So the question could
not be deferred: `𝒜` and `L` were free parameters, and which values are physical?

**They were never free.** On a maximal slice the background constraint is
`R̂ = 16πGρ̄`, so a profile does not choose its matter — the matter is whatever
the profile implies. A product tube of area `𝒜` has `R̂ = 8π/𝒜`:

| throat | `ρ_tube/ρ̄` |
|--------|-------------|
| #261–#264, `𝒜 = 4π` | **`1/3`** |
| matched, `𝒜 = 4π sin²a` | **`133`** |

Neither is the ambient's own fluid.

**The throat that is forced.** Ask for one needing *no* matter (`R̂ = 0`) and
gluing on with *no* surface layer. `R̂ = 0` integrates to `f'² = 1 − f₀/f`, and
smooth gluing at mouth radius `a` forces

```
f₀ = sin³a          L = 2[sin³a·arccosh(1/sin a) + sin a cos a] ≈ 2a
                    I = ∫ds/f² = 4 cos a / sin³a
```

**No free parameter is left.** The closed forms check against quadrature to
`1e-12`, and the conductance is exactly `N₀(a,4)/4` at every radius.

**And it is an Einstein–Rosen bridge, which derives its mass.** `R̂ = 0`, `K = 0`
and a spherical neck don't merely *permit* one — they are one. The
time-symmetric Schwarzschild slice `ds² = dr²/(1−2M/r) + r²dΩ²` in proper radial
distance is exactly `f'² = 1 − 2M/f`, so **`f₀ = 2M`** and

```
M = sin³a / 2          the throat's mass, from the size of the excised mouth
```

Three quasi-local masses agree exactly — the Schwarzschild parameter, the
irreducible mass `√(A/16π)` (the neck area is `16πM²`), and the Hawking mass,
which is *constant* along the vacuum piece. And the gluing condition **is** a
mass statement: `(f/2)(1−f'²)` is `f₀/2` on the throat and `sin³χ/2` on the
ambient, so *the seam is smooth exactly when the Hawking mass doesn't jump.*

What it does not say, each asserted in the tests: no asymptotic region so **no
ADM mass** (the derived mass is quasi-local, unambiguous only because the
Hawking mass is constant); it is **dimensionless**, `M/R`, which is all PR #52's
scale-modulus theorem allows; and both ends sew into the *same* `S³`, so it is a
handle of Misner's kind.

**And the neck — narrowly.** What is proved is an identity:

```
H = 0 at the neck ,  K_ij = 0 on the slice   ⟹   θ₊ = θ₋ = 0
```

so the neck is a minimal surface and a **marginal sphere (MOTS) of this slice**.
That is a statement about one surface in one slice, and it is all of it. It is
**not** shown to be an *apparent horizon* — that is the **outermost** MOTS, a
global condition on the slice nothing here evaluates — and it is **not** shown
to be *non-traversable*, which is a property of the Lorentzian development,
while this is spatial initial data with no lapse chosen. Non-traversability
*does* follow if the development is additionally taken to be the standard vacuum
Schwarzschild/Einstein–Rosen one, and it is under **that** added assumption that
this becomes the vacuum face of the arc's earlier "connected implies exotic".

**There is no cavity.** `∇² + R̂/2` with `R̂ = 0` is the plain Laplacian:
`(f²u')' = 0`, solutions monotone. #264's resonances at `kL = nπ` and its sign
flips were properties of **matter in the tube**, not of a throat.

**And the sign reverses.**

```
ΔA/A = ( +6.64 , +8.58 )      in units of 2πG      — both mouths OPEN
```

**The mechanism is one number.** Split any symmetric two-port as
`Y = G·[[−1,1],[1,−1]] + shunt·I`. The shunt is the flux a *uniform* potential
drives into the throat, and `(f²u')' = 0` makes it vanish **identically** for a
vacuum tube. Scanned over eight orders, the conductance never changes the sign;
the shunt passes through a pole near `2e−03` and flips it. #264's tube sat at
`6.07`.

**What it costs.** The response is `3000×` larger and grows as `a⁻³`, because a
throat with zero shunt barely lifts the constraint's exact degeneracy. The sign
is robust; **the amplitude at which linearising is legitimate is now the binding
question**, and this round does not answer it.

### Release hardening — the two-port in closed form

*A correctness repair to what #265 shipped, not a new result. The geometry and
the answer above are unchanged.*

`f'² = 1 − f₀/f` is `f = f₀cosh²x` with `ds = 2f dx`, which turns
`(f²u')' = ℓ(ℓ+1)u` into `y'' = (2ℓ+1)²y` under `R = y/cosh x` — constant
coefficients. The half-length has `e^{−X} = tan(a/2)` **exactly**, so with
`k = 2ℓ+1` and `q = tan^{2k}(a/2)` the whole two-port is rational in `q`:

```
D_ℓ = −2π sin a [ k(1+q²)/(1−q²) − cos a ]      C_ℓ = +4π sin a · kq/(1−q²)
```

That is now the production `admittance`. The Riccati solve is **kept as a
validator** (`admittance_riccati`), not deleted — because the closed form needs
something independent to be checked against.

**What was wrong with it.** It formed the cross term as `½(s − t)` from two
eigenchannel values that agree to more digits than the solver carries. At
`a = 0.05` it is right to `1e-13` at `ℓ = 0` and to four figures at `ℓ = 1`, but
at `ℓ = 2` it returns **`−1.17e-14` for a true `+3.00e-16`** — wrong sign, and
larger than the answer. The diagonal was never affected (`1e-14` in every
channel); the floor is `~1e-12` in `|Y₁₂|`. The tests pin **the boundary, not
the `39×` factor** — that factor is one solver's step sequence in one build.

**#265's answer is unchanged**: `solve_matching` uses only `ℓ = 0` and `ℓ = 1`,
both above the floor, and `ΔA/A` moves in the *thirteenth* digit.

**And the closed form gives the hierarchy for free.**

```
C_ℓ → 4π(2ℓ+1) sin a · tan^{4ℓ+2}(a/2)  ~  a^{4ℓ+3}
```

— fitted exponents `3.000000 / 7.000004 / 11.000009 / 15.000013`. Each unit of
angular momentum costs **four powers of the mouth radius**, and `C₀/C₁ = 8.5e5`
at `a = 0.05`. The four `n = 1` ESU harmonics split locally as `1 ⊕ 3`
(`X⁰ = cos χ` is `ℓ = 0` at the mouth, the three `Xⁱ = sin χ n̂ⁱ` are `ℓ = 1`),
so the kernel's two pieces cross four powers of `a` apart.

> The statement this supports, and the only one: **the static scalar Laplacian
> on this scalar-flat spatial throat suppresses the local `ℓ = 1`
> mouth-to-mouth channel by `~10⁻⁹` at `a = 0.05`, while preserving a much
> stronger monopole channel.**

It is **not** a claim that orientation information cannot cross the throat —
one operator, one slice, **no lapse chosen**, and the `ℓ = 1` channel is small
rather than zero. Two scope corrections go with it: `f_min > 0` is forced
within #265's class (spherical, scalar-flat, `K_ij = 0`, `C¹`-matched), not by
Einstein's equations generally; and `physical_throat` supplies **spatial initial
data only**.

```bash
python -m experiments.closure_ledger.physical_throat_probe

python scripts/geometrodynamics_v65_throat.py --still v65.png
```

Full write-up: `docs/which_throat_is_physical.md`.

## Two waves, and where they connect inner to outer

**PR #266 — revisiting v46 before the nonlinear arc.** v46 put one scalar wave
on the great circle through a source and its antipode and found a negative
result: the curve is a **graph** `r = f(σ)`, so its radial winding is identically
zero — **a height field cannot wind**, and one wave never meets itself.

Everything since has needed *two* waves. So: one pulsing **outward**, one
pulsing **inward**, both refocusing at the antipode — do they connect?

**Yes. At the antipode, on the seam, at the refocus.** One curve reaches exactly
`R_inner` and the other exactly `R_outer`, and glued those are one point.

**And the threshold is not a new number.** A single wave crosses the seam when
`εu = gap/2`; the pair spans `2ε|u|` of the radial circle and touches through it
when that reaches `gap` — the same inequality:

| | |
|--|--|
| single-wave wrap gain | `0.220059` |
| pair-contact gain | `0.220059` |
| difference | `0.0` |

v46's *"the wave comes back inside the circle"* and this round's *"the two
pulses connect inner to outer"* are one event described twice.

**What two waves can do that one cannot.** Above threshold the tangency opens
into an arc, bounded by two crossings, on which the band between the two curves
covers the *whole* radial circle — no radius is left outside the pair. A single
wave past its own wrap threshold still leaves every radius outside itself,
because it is a graph. **Two graphs bound a band, and a band can be radially
surjective.**

Two details worth having: the antipodal refocus is a **rarefaction**, so it is
the *inward*-driven wave that reaches `R_outer`; and meeting mid-flight — the
two travelling pulses crossing at the quarter points — is the *worst* place, not
the best, costing `7–9×` more because they partially cancel.

**What is still put in.** The crossing rule is a **representation** choice, the
field is **linear** on a **fixed** background so the two waves do **not
interact**, and the gain is a **display** amplitude. This is not a claim that two
physical waves reconnect a throat — it is that v46's obstruction does not apply
to two of them.


### Off the degenerate axis: the offset, and the signs

The co-located pair is the **most degenerate** configuration the construction
has: both wave histories hang off one antipodal axis, so bringing them together
at one pole invites either exact overlap or exact cancellation and tests
neither. Two knobs move off it — the source separation `α`, and the radial sense
each wave is driven in:

    δ = r_A − r_B = ε (s_A u_A − s_B u_B)

**Opposed** signs give the *sum* of the two fields; **like** signs give the
*difference*. That one line is the whole asymmetry.

**A correction to the framing first.** Inner–inner and outer–outer are not two
cases, they are **one** — `|δ|` agrees exactly, as a difference of fields, to the
bit. Flipping both signs is a reflection about `R_mid`, an isometry of the glued
radial circle. So the picture cannot distinguish them, and the reason is worth
being blunt about: *the radial direction here carries the field's amplitude, not
its direction of propagation.* That is a limitation of the representation, stated
rather than worked around.

**The bisector.** `σ = α/2` is equidistant from both sources, so `u_A = u_B` on
it identically. A **like-signed** pair therefore has `δ ≡ 0` there — the two
curves are the *same curve*, never separated by a hair, so no gain however large
carries them through the seam. An **opposed** pair has `δ = 2εu(α/2)`, as large
as it can be. There are two such axes, `α/2` and `α/2 − π`, and the far one is
the cheaper because it sits nearer the antipodal caustic.

**So yes — and the offset is what produces it.** Above threshold the opposed
pair's contact opens into an arc **centred on the bisector to machine zero**,
off both the sources and their antipodes, on which the like-signed pair's contact
set is **empty at every offset tested**. At `α = 0` the bisector collapses onto
the source axis and there is nothing off-axis to find: the degeneracy, recovered
as a coordinate fact.

| | |
|--|--|
| threshold at `α = 0` → `α = π` | `0.2201` → `1.6639` (`7.56×`) |
| timing, against the pulse crossing `t = α/2` | `0.0031π` |
| price of exclusivity | `1.68–3.74×` |
| cheapest point pins to one of the four axes from | `α = 0.125π` |

**Exclusive is not cheap.** The globally cheapest connection stays on a source
axis or an antipode and is available to *both* pairs. Both numbers are reported.

```bash
python -m experiments.closure_ledger.two_wave_slice_probe

python scripts/geometrodynamics_v66_two_wave_slice.py --still v66.png
```

Full write-up: `docs/two_wave_slice.md`.

## One field on one surface — and the antipode is parity-dependent

**The objection lands.** v66 drew **two** curves and asked whether their images
meet through the glued seam. Two curves in one frame are two surfaces, and
reading their overlap as a connection is a statement about a picture, not about
a field. If the two contributions are two pieces of one scalar deformation of
one surface there is only ever **one** curve, `r = R_mid + ε(s_A u_A + s_B u_B)`,
and the question is whether *it* reaches `R_outer` at one `θ` and `R_inner` at
another.

**The repair costs nothing, which is the surprise.** `δ = r_A − r_B` — the
quantity v66 plotted as a separation between two curves — **is** the one-surface
deformation with the second sign flipped, the same array to `2` ulp of `R_mid`.
Every v66 number survives, with the two configurations **swapping names**:

    v66 "like-signed"  ==  one surface OPPOSED
    v66 "opposed"      ==  one surface LIKE-signed

which inverts v66's headline: its cheapest-when-co-located result belongs to the
*like* pair, and its identically-zero bisector is the **node of the opposed
field**, which is where it belongs.

**Coincidence cancels exactly.** `α = 0` with opposite orientation gives `u ≡ 0`
at every time — zero, not small — so no amplitude connects it.

**The monochromatic law.** `u = −2A sin(mα/2) sin(mθ − ωt)`, verified
symbolically and on a grid, so `B = 2A|sin(mα/2)|` and the optimum is
`α* = π/m` — **half a wavelength**, with the antipode simply the `m = 1` member.
And `sin(mπ/2)` is `±1` for odd `m` and `0` for even `m`:

> **the antipode is parity-dependent** — maximal for odd modes, exactly
> cancelling for even ones. Not a visualization effect: on `S³` the same parity
> is `Z_n(π) = (−1)ⁿ`, checked at `2.0000` and `0.0000`.

**Where the two models part company.** A zonal harmonic is *centred*
(`Z_n(0) = 1`, `|Z_n| ≤ 1`), so `|Z_A − Z_B| ≤ 2` with equality only where one
focus sees `+1` and the other `−1` — exactly the antipode with odd `n`. For
zonal modes `α* = π` for **every** odd `n` and it *saturates* the bound; half a
wavelength reaches only `1.10–1.41`. For even `n` the antipode cancels and
nothing reaches the bound at all. **The parity carries across the two models
exactly; the location of the optimum does not.** (The kernel here is `n = 1`,
odd — so for it the antipode is optimal *and* saturating.)

**A pulse is not a mode.** v46's field carries a power-weighted mean `n ≈ 10`
with fifteen modes holding 90% of the power. For it the `1/|sin|` divergence is
confined to about one pulse width; past that the pulses stop overlapping and the
threshold **saturates** at `0.2163` instead of falling to `0.13`.

**The chord, and what it costs.** At the optimum the two extrema sit `π/m` apart,
so `L = √(D² + 4 R_out R_in sin²(π/2m))` falls from `2.000` to the purely radial
gap `0.520`. At fixed *display amplitude* the span is flat at `2.0000` across the
whole family — same deformation, shorter connection. But `E ∝ ω²A²`, so at fixed
**energy** `A ∝ 1/ω` and the span falls as fast as the chord: the highest mode
that still spans the gap is `m = 7`. No favourable frequency is claimed — that
needs an energy normalisation and a packet focusing law this model lacks. What is
claimed is narrower: **a frequency slider cannot hold displacement fixed and then
be read as constant-energy physics.**

```bash
python -m experiments.closure_ledger.one_surface_probe
```

Full write-up: `docs/one_surface.md`.

### Where each front sits on that surface

Asked inside the one-surface object: as the two axes move apart, where do `A`
and `B` individually sit, what are their signs, and how does the surface answer?
Only the **surface** is ever a closed curve in the annulus; the two
contributions appear on graphs of field against `σ`, and the annulus panels
colour the single curve by which front owns each arc. An inward dent is a
negative contribution and an outward one positive, so each front's sign reads
straight off the surface.

| offset `α/π` | peak `c_A` | peak `u` | amplification | overlap arc |
|--|--|--|--|--|
| `0.00` | `0.6710` | `0.0000` | `0.0000` | `0.794` |
| `0.15` | `1.1796` | `1.1994` | `1.0168` | `0.148` |
| `0.25` | `1.1796` | `1.1978` | `1.0155` | `0.000` |
| `1.00` | `1.1796` | `1.1935` | `1.0118` | `0.000` |

At `α = 0` the surface is a **perfect circle** — the contributions cancel
identically. Past one pulse width they stop overlapping at all, and the total is
`1.012–1.017×` **one** contribution, peaking exactly where a single one does.

> The offset does not turn interference on. It turns the **cancellation** off,
> and what is left is two nearly independent dents in one surface.

```bash
python scripts/geometrodynamics_v68_two_fronts.py --still v68.png
```

## The centre as a finite bearing, not a point

**The point was doing work, and it cost `f = 0`.** Every picture in this arc put
a point in the middle, because a point is where the clock-hands story works: two
radial arms `P_A → O → P_B` change direction at `O` for free, since at a point
there is no angular direction left to change. That is the property the
connection is wanted for. It is bought where the geometry stops existing.

**Blow it up.** Keep `dℓ² = ds² + f(s)²dΩ²` and set `f_min = f₀ > 0`. The middle
becomes a small circle in the 2-D cross-section, or the space of radial
directions `S^{d−1}` — `RP^{d−1}` if the clock hand is an unoriented axis.

**The arms are the repo's own geometry with the symmetry dropped.**

```
L(F) = √(F(F − f₀)) + f₀ arcosh√(F/f₀)        I(F) = (2/f₀)√(1 − f₀/F)
```

Set `f_o = f_i = sin a`, `f₀ = sin³a` and these reproduce
`VacuumThroat.length()` and `.resistance()` **bit for bit**. So `L_o ≠ L_i` is
ordinary — scanned to `L_o/L_i = 437` with nothing breaking — and the old
vacuole's one shared arbitrary radial gap is gone.

**Two caveats that the first draft of this section missed.** Unequal arms are
an *intrinsic* statement: they are truncations of one scalar-flat profile. They
are **not** a fully matched #265 throat with the symmetry removed, because `C¹`
matching to a unit round `S³` forces `f₀ = sin³a_j = F_j³` at each end, and one
common `f₀` fixes `F = f₀^{1/3}` uniquely — so *both matched mouths have the
same scale*, and asymmetric arms need asymmetric ambient attachments. And
`R_inner/R_outer` does not lose all significance: the *endpoint scale ratio* is
load-bearing, since `w_i/w_o = f_i/f_o`. What loses significance is the
vacuole's **arbitrary drawing ratio**.

A feature of angular width `Δθ` has physical width `f(s)Δθ`, so it is
**squeezed into the bearing and let out again** rather than teleporting at
fixed size; `w_i/w_o = f_i/f_o`, with `f₀` nowhere in it.

**And the result is a correction to the proposal.** Turning through `α` around
a bearing of radius `f₀` looks like it should cost the arc, `f₀α`. That is exact
for the down-turn-up route and an honest upper bound — but the geodesic spreads
the turn instead. **Why it wins is Pythagoras, not leverage:** the corner pays
its angle as *pure transverse motion* (`f₀α`, first order), while the geodesic
*tilts* motion that is already radial, and a tilt costs `≈ ½f²θ′²ds` — second
order.

```
T(α) = α²/(2I₂) + O(α⁴)  ,   I₂ = ∫ds/f²      quadratic, not linear
```

Leading order, not an exact law: the exact object is the integrated
`turn_cost`, which is `0.9248` of the quadratic value at `π`. And `I₂` is **not
a new quantity** — it is `physical_throat`'s own resistance, so
`T(α) = α²(4π/I₂)/(8π)` at angular dimension `q = 2`. Long arms give
`I₂ → 4/f₀` and `T → f₀α²/8`. Checked against the integrated geodesic, not the
expansion: exact to `8e-05` at `α = 0.1`. The geodesic spends `1.25%` of the
arc at `α = 0.1` and `36%` at `π`.

*Where* the turn happens also came out backwards in the first draft: `θ′ = h/f²`
puts the angular rate **highest where `f` is smallest**, so the geodesic hugs
the neck — `76%` of the turn inside `f < 2.4 f₀`.

**So the property the point was wanted for survives, more strongly than
proposed.** A half-turn costs `8.4e-04` of the arms — and `π` is the largest
separation there is, since `bearing_distance` reduces any pair of directions to
`[0, π]`, so that is the worst case over the whole configuration space. **No
reachable orientation makes the hinge cost as much as the journey.** (An
earlier draft quoted a `104` rad "break-even angle"; withdrawn — it extrapolates
outside both the law's domain and the configuration space.) And the point model
is the **limit**, not a rival: `T` is linear in `f₀`.

**Intersection becomes overlap.** Two fronts land on the bearing at angular
positions. *Whether* they meet is a question about angles with `f₀` nowhere in
it; *how big* the meeting is, is `f₀ × (overlap angle)`. So `f₀ → 0` does **not**
make every route meet. **As `f₀ → 0` the angular incidence survives and the
physical interaction region collapses:** direction of arrival and angular
overlap never involved `f₀`, while the region where the fronts actually share
space is `f₀ ×` (overlap angle) and goes to zero — as does the separation of two
that miss. (The *size* law is dimension-dependent: a length on the drawn `S¹`
cross-section, an area `~f₀²` on `S²`, a volume `~f₀³` on `S³`. The yes/no
criterion is not.) The distinction survives as a yes/no and disappears as a length.

**The deeper identity.** `I₂` appearing in both the monopole conductance and
the hinge is not two calculations sharing a number. They are **one Dirichlet
form** on the tube:

```
minimise  E[φ] = ∫ f² φ'² ds   at fixed increment   ⟹   (f²φ')' = 0
```

Read it at `φ = u` and the conserved current *is* the monopole flux `4πf²u′`;
read it at `φ = θ` and it *is* Clairaut's `h = f²θ′`. **Static monopole flux
and infinitesimal throat rotation follow the same spatial weighting** — but
only at one dimension, and the first draft did not say so. The azimuth's weight
is the *metric* coefficient `f²` for every bearing dimension; the monopole's is
the *volume* element `f^q` for an `S^q` cross-section, so its resistance is
`∫ds/f^q`. They coincide **exactly at `q = 2`** — the physical case here, a
3-D spatial throat, and what `physical_throat` carries. The great-circle
reduction of the hinge *is* dimension-free; this identity is not.

The sharpest form is about the *profiles*, not the numbers: normalised to their
own total, the monopole potential and the geodesic azimuth are the **same
function of position** along the tube, with the deviation falling as `α²`
(`4.9e-03` at `α = 1` → `4.9e-07` at `α = 0.01`, a clean 100× per decade).
That is why *infinitesimal* is the operative word.

**And the moment hierarchy says what is universal.** Expanding exactly,

```
T(α) = α²/(2I₂) − α⁴I₄/(8I₂⁴) + O(α⁶)      I_n = ∫ds/fⁿ
shape = T/(α²/2I₂) = 1 − α²I₄/(4I₂³)
```

- **The universal thing is the leading functional form** `α²/(2I₂)` — every
  neck obeys it with its own `I₂`, which is *not* itself universal (`4/f₀` here
  against `π/f₀` there). What `I₂` has is a second job: at `q = 2` it is also
  the monopole resistance.
- **`I₄` is the first *additional independent* moment**, and where the neck's
  shape first shows. Not the entire profile dependence: `I₆` and beyond enter
  at `O(α⁶)` and matter by `α = π`.

| profile | `I₂` | `I₄` | shape |
|--|--|--|--|
| scalar-flat | `4/f₀` | `32/(15f₀³)` | `1 − α²/120` |
| hyperbolic `√(f₀²+s²)` | `π/f₀` | `π/(2f₀³)` | `1 − α²/(8π²)` |

`1/120` against `1/79` sets how fast they separate — the shape difference is
`α²(1/79 − 1/120)`, so `4.3e-05` at `α = 0.1` growing to `3.5e-02` at `π`. (An
earlier draft said they "agree to eight digits" at `α = 0.1`; **withdrawn** —
eight digits is how well each profile matches *its own* quartic law, not how
well the two match *each other*.) Every moment has a closed form: substitute
`t = f′` and `I_n = (2/f₀^{n−1})∫₀^T(1−t²)^{n−2}dt`.

**How far it generalises, measured.** `T = α²/(2I)` never used the profile, and
holds to `1.3e-04` on an unrelated one (`f = √(f₀²+s²)`, `I → π/f₀`). The
`O(α⁴)` correction does not: the shape at `α = π` is `0.9250` scalar-flat
against `0.8886` hyperbolic. **The quadratic law is about necks; the ~8–11%
deficit at a half turn is about a particular neck.**

**Scope.** Geometry only — a metric, its geodesics, and an angular width carried
along them. No field equation is solved, nothing evolves, and this does not
choose between the three candidate bulks; it works the finite bearing out far
enough to be compared with the finite caustic ring and the finite neck with
moving attachment maps.

```bash
python -m experiments.closure_ledger.regularized_center_probe

python scripts/geometrodynamics_v69_regularized_center.py --still v69.png
```

Full write-up: `docs/regularized_center.md`.

## What higher dimension does to the bulk picture

**3-D intuition models an extra dimension as "the same sphere, with another
direction available."** It is not, and the ways it is not bear directly on the
finite bearing.

**A correction first.** The familiar "spheres peak at 5D" is the **unit ball**,
`V_d = π^{d/2}/Γ(d/2+1)`, peaking at `d = 5`. The unit *sphere*'s surface,
`A_{d−1} = 2π^{d/2}/Γ(d/2)`, peaks two later at `d = 7` — the sphere being
`S⁶ ⊂ ℝ⁷`. And neither is a fact about dimension alone: the two carry different
units at different `d`, so comparing across `d` picks a length scale and the
peak follows it — `d = 1` at `R = 0.5`, `d = 5` at `R = 1`, `d = 100` at
`R = 4`.

**Almost all of a sphere is at the equator of any chosen point.** The shell at
angle `χ` has measure `∝ sin^{n−1}χ`, and `χ = π/2 + δ` gives
`e^{−(n−1)δ²/2}` — a band of width `~1/√n`. Measured, `std(χ)·√n` runs
`0.9669 → 1.000000` over `n = 2 … 1000`.

**So the antipode is an extremely non-generic relation.** Random directions have
`x·y` of width `1/√n`, so `α → π/2`; the fraction with `α > 0.99π` is `3.2e-04`
on `S²` and *zero* in `2e+05` samples by `n = 10`. Selecting `x ↔ −x` is not
"pairing a point with the far one" — it picks a vanishing-measure relation out
of an overwhelming nearly-orthogonal majority, and gets **more** non-generic as
dimension rises. (That does not make it correct; it removes the bland reading.)

**And the collapse is `fⁿ`, not `f`.** For `dℓ² = ds² + f²dΩ_n²` the transverse
measure of an angular patch is `fⁿ dΩ_n`:

| squeeze `f₀/F` | `S¹` | `S²` | `S³` | `S⁴` |
|--|--|--|--|--|
| `1e-03` | `1.0e-03` | `1.0e-06` | `1.0e-09` | `1.0e-12` |

> **The angular overlap can stay finite while the physical overlap collapses as
> `f₀ⁿ`.**

**Which `n` is physical depends on which object is drawn**, and that has to stay
explicit or the exponent will migrate between objects that do not share it. `n`
is the dimension of the object's own transverse sphere: the drawn cross-section
is `S¹`, **PR #265's spatial throat is `S²`** — `n = 2`, understating the drawing
by a **thousand**, and its neck area is the measured `4π f₀²` — while the `S³`
that gives the millionfold figure is a *bearing in a four-spatial-dimensional
embedding*, a different object. Nothing here licenses carrying an exponent from
one to the other.

That changes what PR #268's picture is a picture *of*: not two thick ribbons
squeezing until they touch, but **large angular structure → tiny proper measure →
large angular structure**, angular labels intact throughout — two angular sectors
becoming coincident on a greatly compressed direction space, then re-expanding
into different radial sectors. The yes/no overlap criterion is untouched; it was
always angular, so the drawing need not force a dramatic macroscopic crossing:
`ΔΩ` stays constant while `V_overlap ∝ f(s)ⁿ`, and the intersection is real in
the topology of the coordinate mapping and extraordinarily small in measure.

**And the finite centre is a routing manifold, not a hub.** Directional capacity
— how many directions a bearing resolves at a given angular resolution — is
**dimensionless**: at `20°` an `S³` bearing distinguishes `113.5` directions and
an `S²⁰` one `2.2e+10`, and the number is *bit-identical* across `f₀ = 1e-01 …
1e-07` while the bearing's proper measure runs `1.974e-02 → 1.974e-20`. So the
singular centre is not obtained because every direction becomes equivalent
there. It is obtained because an entire finite direction space is compressed to
zero proper measure with its angular structure intact — which is a different
statement, and the `f₀ → 0` limit separates three things that were being run
together: angular incidence survives, directional capacity survives, proper
interaction measure collapses as `f₀ⁿ`.

**Orientability flips with parity — and this repo uses two quotients that are
always opposite.** `ℝP^n` is orientable iff `n` is odd (`det(−I_{n+1}) =
(−1)^{n+1}`). The repo carries the **spatial** quotient `S^d/± = ℝP^d` *and* the
**two-body exchange** space `(ℝ^d∖0)/± ≃ ℝP^{d−1}`, one apart:

| spatial `d` | spatial `ℝP^d` | exchange `ℝP^{d−1}` | `π₁`(exchange) |
|--|--|--|--|
| 2 | non-orientable | orientable | `ℤ` (braid) |
| **3** | **orientable** | **non-orientable** | `ℤ₂` |
| 4 | non-orientable | orientable | `ℤ₂` |

At `d = 3` the Pin⁻ structure that makes the throat a spinor lives on the
*exchange* `ℝP²`. Raising the spatial dimension by one **swaps which is which**,
so that mechanism would have to be re-derived rather than carried across. Said
as a consequence *if* the dimension moved — not an argument that it should.

**And `S³` is not a generic sphere.** `S³ ≅ SU(2)`; it is parallelizable (the
frame `q·i, q·j, q·k`, orthonormal to `6.7e-16` everywhere), which only `S¹`,
`S³` and `S⁷` are; and it carries the Hopf fibration `S¹ → S³ → S²`. `S²` admits
no nowhere-zero tangent field at all. So none of the above is a smooth trend to
extrapolate.

**Scope.** Measure, orientation and frames on round spheres. No field equation
is solved, nothing evolves, and the embedding-centre reading of the bearing —
that `f₀ > 0` restores at finite size the direction space the origin crushed —
is marked as an interpretation, not a result.

```bash
python -m experiments.closure_ledger.hyperspherical_probe

python scripts/geometrodynamics_v70_hyperspherical.py --still v70.png
```

Full write-up: `docs/hyperspherical.md`.

## The first evolved Einstein equations (PR #270)

Every gravity result above this line is stationary, weak-field or linearized, and
`THESIS.md` has said in the same words across five rounds that *"the strong-field
endpoint (a horizon / a resolved throat) is left for full numerical relativity."*
Nothing in the tree evolved the Einstein equations in time. This round is the
first that does, at the highest-symmetry 4+1 problem — `D = 5`, spherical
symmetry, one minimally coupled massless scalar — in **horizon-penetrating**
ingoing Eddington–Finkelstein coordinates,

    ds² = −A(v,r) e^{2δ(v,r)} dv² + 2 e^{δ(v,r)} dv dr + r² dΩ_n² ,   n = 3

Vacuum is not an option: Birkhoff in `D` dimensions makes Tangherlini the unique
spherically symmetric vacuum, so the scalar is the dynamical content.

**Derived, not recalled.** The metric, connection, Ricci and Einstein tensors are
built in sympy for **general `n`** and the system self-checks at both `n = 3` and
`n = 2` — the latter being the known `D = 4` system, which is what validates the
general-`n` derivation. What comes out:

    rr    →   ∂_r δ = (κ/n) · r · (∂_r φ)²
    vr    →   (r^{n−1} e^δ A)' = (n−1) r^{n−2} e^δ
    wave  →   2 r^n ∂_r(∂_v φ) + n r^{n−1} ∂_v φ + ∂_r(e^δ A r^n ∂_r φ) = 0
    vv    →   an independent equation containing ∂_v A — never used

The `vr` result is the surprise: it is not an ODE to integrate alongside `A`, it
is an exact **quadrature**.

**The Einstein equation the code never solves converges at second order.** The
hierarchy *solves* `rr` and `vr` on every slice, so their residuals are
identically zero and testing them would be circular. `vv` is the one independent
component left over, it carries `∂_v A`, and the code never forms `∂_v A` for any
other purpose:

| points | spacing | max abs `vv` residual | rate |
|--|--|--|--|
| 400 | `0.0501` | `1.5511e-04` | — |
| 800 | `0.0250` | `3.9070e-05` | **`1.989`** |
| 1600 | `0.0125` | `9.7862e-06` | **`1.997`** |
| 3200 | `0.0063` | `2.4478e-06` | **`1.999`** |

Stated as the characteristic-scheme *analogue* of a Hamiltonian/momentum
constraint test, not as one. Two exact solutions pin the scheme first:
Tangherlini comes back at machine precision (`1.6e-15`, `δ ≡ 0`, `ψ ≡ 0`), and
the closed-form flat mode `φ = cos(ω(v−r))J₁(ωr)/r` is reproduced at rate
`2.003`.

**And a regular centre forbids a trapped surface.** The `vr` quadrature reads
`r^{n−1}e^{δ(r)}A(r) = (n−1)∫₀^r s^{n−2}e^{δ(s)}ds` — a positive integrand over a
positive interval — so `A > 0` strictly for `r > 0`, **identically**. Four
profile families driven to `min A = 5.63e-03` confirm it never crosses. So
**horizon formation is not observable in this gauge**, and the criterion has to
be posed as the loss of central regularity rather than as `A` changing sign. That
is a statement about the chart, not the physics: collapse still happens, and it
is why production characteristic codes use *outgoing* null cones or excise the
centre.

**A discrepancy found in passing, reported and not acted on.** The
Schrödinger-form potential for a minimally coupled massless scalar with
`ψ = r^{n/2}φ` is `A[(ℓ(ℓ+2) + 3/4)/r² + (9/4)r_h²/r⁴]`, while
`tangherlini.radial.V_tangherlini` carries `A[ℓ(ℓ+2)/r² + 3r_h²/r⁴]` — a
difference of exactly `3A²/(4r²)`, reproduced to `5.4e-16`. The flat limit
settles it: `ψ = r^{1/2}J_{ℓ+1}(ωr)` gives `V → ((ℓ+1)² − ¼)/r²`, matched to
`4.3e-16`. **Nothing was changed** — `V_tangherlini` is consumed by roughly fifty
probes and by several derived constants, so acting on it is a decision about the
repository's published numbers, not a side effect of a dynamics round.

**What this round did not earn.** The perturbation spectrum and the retarded
outer→inner transfer function are **not delivered**. Two horizon-penetrating
time-domain constructions — a Kerr–Schild slicing of the same chart and a
tortoise `(t,r*)` evolution with the derived potential — are both stable and both
*converged*, and they disagree: real parts within `0.3%` at `ℓ = 1`, damping
rates apart by `37%` (`1.01622 − 0.36231i` against `1.01876 − 0.26404i`). So no
quasinormal frequency is reported, and the transfer function is not built because
it is a ratio of the same two signals. **A converged number is not a correct
number.** *(Settled by PR #274 against published values: Matyjasek 2021 gives
`1.01601691 − 0.36232802i` at `ℓ = 1`, confirming the Kerr–Schild code to
`0.005%` and excluding the tortoise damping, which is `27.1%` off. #270's own
prime suspect — the Kerr–Schild inner cut — was the wrong code. See
`docs/ringdown_cross_validation.md`.)*

**Scope.** Classical, spherically symmetric, one massless scalar, second-order
accurate and stated as such. Horizon *persistence* is shown only on a seeded
background, where it is exact; a dynamically formed horizon is not evolved.

```bash
python -m experiments.closure_ledger.tangherlini_dynamics_probe

python scripts/geometrodynamics_v71_tangherlini_dynamics.py --still v71.png
```

Full write-up: `docs/tangherlini_dynamics.md`.

## The radial scalar operator, corrected (PR #271)

PR #270 found — while doing something else — that `tangherlini.radial.V_tangherlini`
was not the master potential of a minimally coupled massless scalar, and reported
it without changing anything. This round makes the correction and **prices** it.

**The result, first.** Correcting the operator leaves the radial spectrum and the
cross-`ℓ` structure nearly intact, but **removes the numerical support for the
claim that the canonical `R_OUTER = 1.26` geometry generates the locked lepton
pinhole `γ = 22.5`.** Eigenvalues move at `10⁻³`; the cross-`ℓ` operator is
unchanged to `3.6e-15`. What moves is the barrier sum — and the lepton chain is
far more sensitive to it than had been measured: `d ln m_μ / d ln γ = −17.5` at the lock (secant `−16.6` over `22.331…22.836`), so
a half-percent error in `γ` is a nine-percent error in the muon. Under the
legacy operator the canonical geometry nearly produced the lock (`22.453` at
`R = 1.26`, `0.21%` away, masses within `3.8%`); under the corrected operator
**no channel set at `R = 1.26` lands near `22.5`** — `22.331` or `22.836` — and
both damage the ladder at the `15–21%` level. The locked Hamiltonian never sees
`R_OUTER` or the channel set (`del l, n_points, rs, r_outer`), only `γ`, so
enforcing `γ = 22.5` makes alternative geometric roots **observationally
indistinguishable**. *`γ = 22.5` remains required by the locked model; its
claimed derivation from the canonical radial barrier geometry is reopened.*

For `ds² = −A dt² + A⁻¹dr² + r²dΩ_n²` with `ψ = r^{n/2}R`, the unique
first-derivative-free Schrödinger form carries

    V_scalar = A[ ℓ(ℓ+n−1)/r² + n(n−2)A/(4r²) + n A'/(2r) ]

verified symbolically at `n = 2 … 6`. The repository carried
`A[ℓ(ℓ+2)/r² + 3r_h²/r⁴]`, short by exactly `3A²/(4r²)`. The flat limit settles
which is which with no appeal to authority: `ψ = r^{1/2}J_{ℓ+1}(ωr)` gives
`V → ((ℓ+1)² − ¼)/r²`, matched to `2.2e-16` — and the legacy operator **fails**
that test. **A bug, not a convention.**

`V_tangherlini_legacy` is frozen for archived runs; `V_scalar_tangherlini` is the
corrected general-`n` operator; `V_tangherlini` now delegates to it.

**The eigenvalues barely move, and move less as `ℓ` rises** — `+0.1320%` at
`ℓ = 0` down to `+0.0192%` at `ℓ = 5`, overlaps above `0.999998`. An eigenvalue
averages the potential against a bound state, so an `ℓ`-independent shift matters
least where the centrifugal term already dominates.

**The barrier sums are not protected, and the `γ` story swaps:**

| channels | legacy | corrected | `R_OUTER` legacy → corrected |
|--|--|--|--|
| `ℓ = 1..5` | `22.00824` (`−2.19%`) | `22.33119` (**`−0.75%`**) | `1.28737 → 1.26788` |
| `ℓ = 0..5` | `22.45268` (`−0.21%`) | `22.83642` (**`+1.50%`**) | `1.26227 → 1.24614` |

The canonical README claim improves threefold with nothing tuned; the claim that
adding the `ℓ = 0` 5D channel closes the pinhole gap **breaks**, and the sum
closest to `22.5` swaps channel sets. **Withdrawn, not replaced.**

**Exactly invariant:** `ΔV` carries no `ℓ`, so the cross-`ℓ` operator
`V_{ℓ+2} − V_ℓ` is unchanged to `3.6e-15`. Its matrix elements still drift,
because the eigenfunctions do — structure invariant, numbers shifted. Hopf, Pin⁻,
the odd-`k` ladder and antipodal parity have no dependence on this operator and
are **not** re-run: proximity is not dependence.

**One narrow downstream re-derivation, run before merging.** Three geometries
through the *locked* lepton Hamiltonian with nothing retuned: **A** fix
`R_OUTER = 1.26`; **B** enforce `Σ[1..5] = 22.5`; **C** enforce `Σ[0..5] = 22.5`.
**B and C come out bit-identical** — `compute_knotted_lepton_spectrum` discards
`r_outer` outright, so the locked block sees the geometry *only* through the
scalar `γ`, and the channel-set choice leaves no trace in any observable.

> **`γ = 22.5` is the selector; `R_OUTER` is downstream of it.**

Fixing `R_OUTER` and letting `γ` float is what breaks the ladder: `+15.16%` and
`−20.52%` on `m_μ`, against the legacy geometry's `+3.78%`. So **the correction
weakens the "geometry supplies `γ`" story** even while improving the `1..5`
residual in isolation — because `d ln m_μ / d ln γ = −17.5` at the lock (secant `−16.6` over `22.331…22.836`), and a sub-percent
geometric residual is not a small residual in this chain. The channel-set
question is **not decidable by the lepton observables**.

**And the suite never protected any of it.** Flipping the operator broke exactly
**2 tests out of 1582**, both PR #270's own bookkeeping. The `γ` sums, the
`R_OUTER` fixed point and the `1.054` factor are not regression-locked anywhere;
they live in prose. A silent replacement would have sailed through CI.

```bash
python -m experiments.closure_ledger.scalar_operator_audit_probe
```

Full write-up: `docs/scalar_operator_audit.md`.

### The deferred re-derivation, carried out (PR #272)

#271 classified the downstream chains and re-derived none of them. Doing it for
the **quark residual sector** reverses the expectation: all three residuals
derived from the same eigensolver move **toward** their locked values —
`pinhole −1.09% → +0.36%` (sign change), `transport +0.88% → +0.70%`,
`resistance +0.49% → −0.02%`.

One barrier feeds both sectors and the correction moves it once, so the split
verdict is not about geometry — it is **elasticity**: `d ln m_μ/d ln γ = −17.5`
against `d ln m_s/d ln pinhole = +4.8`. The lepton's `−0.75%` residual is a
measured `15.2%` muon error; the quark's `+0.36%` is a measured `1.79%` strange
error.

> **A percentage agreement between a geometric quantity and a fitted knob means
> nothing until multiplied by the elasticity of what it feeds.**

Two things get *weaker* on inspection, under both operators. The derived triple
**never** reproduced the ladder (`3.44%` legacy, `3.78%` corrected, against the
fitted lock's `1.61%`) — and the ordering reverses, because the legacy triple's
larger individual errors partially cancel. And `N` still drifts by O(50) with
the residuals pinned, so `β` remains the model's last fit knob.

One thing gets stronger: read as demands on `R_OUTER`, the corrected operator has
the two sectors **straddle** `1.26` (`1.25645`, `1.26788`) where the legacy
operator put it outside their bracket entirely. Different evidence from the
single-sector fixed point #271 reopened, and weak — a `0.91%` window admits
anything inside it.

```bash
python -m experiments.closure_ledger.quark_residual_reaudit_probe
```

Full write-up: `docs/quark_residual_reaudit.md`.

### The fit manifold, not the residuals (PR #272, same round)

The residual round measured one knob at a time and guessed the gap between
individually-right and jointly-wrong residuals was a scalar relation between
`transport` and `resistance`. **The guess was wrong and so was the object.**
The right one is the response map `J_ia = ∂ln(m_i/m_d)/∂ln p_a` and its SVD.

**Three statements survive, none metric-dependent:**

1. **Four scored masses cannot identify the current parameterisation.**
   `rank J = 4` — capped by the observable count — against **8** first-order
   knobs; `JᵀJ` carries four exact zeros. `pinhole`, `transport`, `resistance`
   and `N` are **not separately constrained**; at most four combinations are.
   Every compensation seen since PR #76 is this.
2. **The positive-knob Jacobian is numerically full row rank**, converged to
   five digits across three decades of step size — a real derivative, not a
   difference artefact.
3. **Per-knob proximity to a fitted value does not determine the effect on the
   spectrum.** The three radial residuals do not compose linearly into the
   ladder.

**Does correcting the operator move toward the data?** Two different questions,
and the first draft of this round ran them together. As *candidate
displacements from the lock*, `cos(J·g, r)` is `+0.464` (legacy) and `−0.616`
(corrected). But the move the correction actually makes is
`Δg = g_corrected − g_legacy`, measured against the residual the legacy triple
leaves — and that is **`+0.873`**, toward it. The two objectives also disagree:
the L2 log-residual **improves** (`0.0548 → 0.0433`) while the repository's
max-relative-error score **worsens** (`3.19% → 3.80%`). Both are true; the
direction of "improvement" is objective-dependent.

**The `1.61%` floor is not structural, and removing it proves nothing.** A
displacement whose largest single knob moves `1.78%` reaches `0.0179%` on an
exact nonlinear re-run — but the rank deficiency guarantees some such
displacement exists.

Three nulls were separated into three different objects: `action_base` is an
**exact gauge** (`H(a) = H(0) + a·I`, killed by the zero-point subtraction
upstream of the anchor); `phase` is an **antiunitary** `Z₂` of the spectrum
(`H(−φ,p) = H(φ,p)*` for *arbitrary* `p`, so the masses are even in `φ` at
every mixing); `partition_mixing` is a **unitary-conjugation** `Z₂`
(`H(−p) = D H(p) D†`). The latter two are quadratic, visible in `q = x²`.

The four **nonzero** singular values span only `22.6×`, so there is no long
hierarchy among identifiable directions — but the rank deficiency is exact, so
the pathology is **structural non-identifiability**, not ill-conditioning. And
which knobs drift is itself structure: identifiable share runs from `1.0000`
(`uplift_asymmetry`) to `0.0003` (`eta_k3k5_minus`).

```bash
python -m experiments.closure_ledger.quark_response_jacobian_probe
```

Full write-up: `docs/quark_response_jacobian.md`.

### Is the CKM a prediction? (PR #273)

The flavor half of the same question, and the audit arc's terminal round. Build
two response maps over **one** parameter chart `x` — `J_M` for the mass ratios
and `J_F` for four *genuinely independent* flavor coordinates
`y_F = (θ₁₂, θ₂₃, θ₁₃, δ)` — take the mass-preserving tangent space
`N_M = ker J_M`, and form **`K = J_F N_M`**. Its rank counts the physically
independent CKM directions reachable without disturbing the masses. Rank is
chart-independent, so unlike a pseudoinverse it smuggles in no metric.

**`rank K = 4`** — the full dimension of the physical flavor space of a unitary
3×3 matrix. Stable over a 3×3 grid of step and cutoff, singular values spread
only `379×`.

**And on the decisive restricted question too.** v4 was built additively over a
*frozen* v3 lock, so the surplus test is whether the v4 calibration freedoms
**alone** span the CKM at fixed masses: `K_v4 = J_F[:,G]·ker(J_M[:,G])` with the
v3 knobs held fixed. It gives **`rank K_v4 = 4`** with `φ_h` also fixed, over
the 10 quantities selected in the re-lock — the three targeted `η`s, the
`η₃₅` **retune**, and the six diagonal shifts. Confirmed nonlinearly: scaling a
mass-preserving direction by `ε` and re-running the Hamiltonian gives clean
`O(ε²)` (ratios `4.00`).

> The mass-preserving parameter freedom spans everything the CKM can be.
> **Fitting it is a successful realisation, not a prediction** — and there are
> zero left-null vectors, so no first-order relation `wᵀδy_F = 0` exists to
> compare against experiment.

**The `φ_h` A/B test.** The library treats `φ_h = π/k₅` as derived and as the
source of CP structure. Holding it fixed leaves **`rank K = 4`** — the other
fitted matrix elements absorb arbitrary CKM data on their own. Its response **is** genuinely a `δ` direction — `0.99978` of its `J_F` column
lies along `δ`, and that share is chart-independent — and releasing it
multiplies the leading singular value by `4.8`, though **that ratio is
chart-dependent** and is kept only as a scoped diagnostic. Efficient, not
identifying. This sharpens PR #173, which found that
*adding* `φ_h` left the observable rank unchanged.

**The census, measured rather than counted.** `rank J_F` restricted to each
knob group: v4 targeted `η`s **3**, `eta_k3k5_minus` retune **1**,
`diag_shift_plus` **2**, `diag_shift_minus` **2**, `φ_h` **1** — and the nine v4
additions with `φ_h` fixed measure **4**, saturating the flavor space. The trace
direction of each diagonal triple is an **exact CKM gauge** (`|J_F·1| ≈ 2e-10`
against `|J_M·1| = 12.5`), which is why three symbols measure rank two; both
realised triples are traceless to `~1e-10`.

**The "+3 parameters for +5 independent observables, net +2" claim is refuted
on the ceiling alone:** a unitary 3×3 CKM has exactly four physical parameters,
so at most four of the nine quoted flavor-CP observables are independent.
Measured net surplus **≤ 0**.

This does not say the v4 numbers are wrong — the `≤ 1%` agreement across nine
observables is real. It says the agreement is not evidence *for* the
Hamiltonian, because the same Hamiltonian could have reproduced any *locally
neighbouring* CKM equally well at the same masses. **Scope:** local and first-order at the v4
lock; the excursions required are not small (`|δx| ≈ 50–80`).

```bash
python -m experiments.closure_ledger.flavor_identifiability_probe
```

Full write-up: `docs/flavor_identifiability.md`.

## Settling PR #270's ringdown cross-validation (PR #274)

PR #270 left one thing explicitly unearned: it had two horizon-penetrating
time-domain codes for a test scalar on a fixed `D = 5` Tangherlini background,
both stable, both converged, and they disagreed — real parts within `0.3%` at
`ℓ = 1`, damping rates apart by `37%` of the smaller value. It refused to quote a
frequency, which was right, and it named a prime suspect: the Kerr–Schild
operator's inner cut.

**The Kerr–Schild code was right. The tortoise code's damping was wrong. The
prime suspect was the wrong code.**

The decisive evidence is **external to this repository**: Matyjasek, *Phys. Rev.
D* **104**, 084066 (2021), [arXiv:2107.04815](https://arxiv.org/abs/2107.04815),
which computes these modes by continued fractions cross-checked against Hill
determinants, agreeing to 11 digits. Converting its scaled frequency
(`ω̃ = ω/T_H`, `T_H = 1/(2π)`) gives `ω = 1.01601691149 − 0.36232802385i` at
`ℓ = 1`, `r_h = 1`. All errors below are relative to it.

| source | `ω` at `ℓ = 1` | damping error |
|--|--|--|
| **published** (CF / Hill) | `1.01601691 − 0.36232802i` | — |
| #270 Kerr–Schild | `1.01622 − 0.36231i` | **`0.005%`** |
| **this round** (characteristic, `h → 0`) | `1.01612 − 0.36244i` | **`0.031%`** |
| #270 tortoise | `1.01876 − 0.26404i` | **`27.1%`** ← excluded |

**So this round is not merely a tie-breaker.** An independent
Gundlach–Price–Pullin characteristic evolution, written from scratch and sharing
no code with either #270 implementation, reproduces a published high-precision
spectrum: `0.018%` and `0.028%` in damping at `ℓ = 1` and `ℓ = 2`. That is a
considerably stronger check on PR #271's corrected radial operator and on the GPP
machinery than internal arbitration could be.

**Note the ordering.** #270's Kerr–Schild code is about **six times more
accurate** than this round's solver. The characteristic scheme was chosen not for
accuracy but because **it applies no spatial boundary condition at all** — the
domain of dependence is the null diamond, so the horizon and infinity are limits,
never boundaries, which is exactly the excision question #270 suspected. An
arbitrator that could fail the same way would have settled nothing. It arbitrated
correctly; it did not out-resolve.

**Naming the denominator.** "`X%` off" is ambiguous, and earlier drafts of this
round used both conventions in different places. The tortoise damping is
**`27.1%`** off in the conventional sense (divided by the published value) — the
number to quote — and equivalently the correct damping is **`37.3%` larger** than
the tortoise value. #270 measured the latter because it had two codes and no
reference. Both are now reported together, neither bare.

**What the external reference exposed about this solver.** The step-size study's
last successive difference in damping is `4.0e-5`, while the finest value's
actual distance to the published one is `1.07e-4` — **2.7× larger** — and
`h = 0.05` lands *closer* to the truth than `h = 0.025`. Discretization is not
what limits this solver; extraction systematics are. **Self-convergence measures
only the error it is refining away: it is a consistency check, not an error bar.**
Nothing internal to a solver can reveal that.

Two `D = 5` identities make the problem checkable, and one of them needed
correcting in this round's own prose. The tortoise correction has the **exact**
closed form `r* − r = −artanh(1/r) = −1/r − 1/(3r³) − ⋯`, so `−1/r` is the
*leading asymptotic behaviour*, not an exact equality as first written. The
physics is untouched — every term decays, unlike 4D's growing `2M ln r`, so the
far field is a pure Hankel wave — and the test is now sharper than the prose was,
checking `r(r*−r) + 1` against the predicted `−1/(3r²)` rather than merely
requiring it to shrink. The second identity, `V → [(ℓ+1)² − ¼]/r²`, is the same
flat-limit identity PR #271 used to settle which radial operator was correct,
reused here as a boundary condition.

**Two honest negatives, reported as negatives.** Frequency-domain shooting
**reproduced** #270's non-convergence rather than fixing it — the root moves with
every knob because the QNM boundary-value problem is exponentially ill-conditioned
in real `r`. Sixth-order WKB by finite differences **diverges** under refinement
(`9.01 → 18.63 → 623.09`). The published continued-fraction calculation is the
frequency-domain reference both were reaching for, so neither is worth forcing.

**What this round cannot do.** There is **no autopsy**. Neither #270 code was
landed in the tree, only their reported numbers, so this establishes *which*
number is right — not which line of the unlanded tortoise code did it. **Scope:**
test scalar on a fixed background, no backreaction, fundamental mode only. `ℓ = 0`
is `0.21%` off the published value, an order of magnitude looser than `ℓ = 1, 2`,
which independently vindicates having quoted it with a wider uncertainty.

**Next.** The retarded transfer function `G_ℓ^ret(t; r_obs, r_src)` from the same
evolution with a compact purely ingoing excitation, gated on three checks before
any physical reading: causal support `G(t < t_null) = 0`; flux conservation
`|R_ℓ|² + |T_ℓ|² = 1`; and late-time ringdown consistent with the **external**
`1.01601691149 − 0.36232802385i` rather than this solver's own fitted number.

```bash
python -m experiments.closure_ledger.ringdown_cross_validation_probe
```

Full write-up: `docs/ringdown_cross_validation.md`.

## The geometric-visualization arc, end to end

Nine rounds (PRs #242–#250) asked one question repeatedly: *given a geometry and
a field on it, which object should be drawn, and what does the choice smuggle
in?* **Four negative results and three positive ones, and the negatives are the
load-bearing half** — almost every apparent obstruction turned out to be a
property of the object being drawn, until the last two, where it turned out to
be an identity no drawing choice can move.

| round | question | answer |
|---|---|---|
| circle slice | can a height field wind through a glued bulk? | **no** — a graph has degree 0 identically |
| seam scale | does the gluing scale rescue it? | **no** — it sets shape, not winding |
| ring & fold | can the ring reach across and intersect? | **no** — it reaches, never crosses |
| normal field | does drawing the *vectors* change that? | **yes** — 520 crossings, threshold `ρ = 1/κ` |
| congruence | does focusing reach a singularity? | **no** — a caustic, and it passes through |
| shell junction | can a detached shell replace exotic matter? | **no** — connected implies exotic, by identity |
| multipole | where does the two-shell coupling start? | **`ℓ = 2`**, screened as `(b/a)^ℓ` |
| wormhole ledger | are the four apparent objects one wave? | **yes** — one balance, and linearity makes it free |
| pair creation | is the antipodal caustic a creation event? | **no** — it is a venue; the threshold needs two waves, and then *forces* a second antipodal interaction |
| pair history | do two closed histories constrain their shared event? | **yes, discretely** — 5 equations, 5 unknowns, rank 5, branch-completely |
| field solve | does a solved field reproduce that ledger? | **yes, exactly** — and it adds the Maslov phase rays could not carry |

The recurring methodological lesson: **a converged number is not a correct
number.** Three of the nine errors the arc caught survived grid refinement and
fell only to an independent construction — a closed form against brute force, an
operator form against a difference, an exact surface against a truncated family.
The tenth round sharpens it: both of its corrections were *representational* —
a receiver placed where the caption said "collapse" but the geometry said
"arrive", and a projection that sent the very point the scene is about to
infinity. **An object drawn correctly can still be the wrong object.**

Full write-up, including what is imported rather than derived and what each
closing result names as its own missing ingredient:
`docs/geometric_visualization_capstone.md`.

## Quick Start

### Verify charge quantisation from pure geometry

```python
from geometrodynamics.hopf import compute_c1

result = compute_c1()
print(f"|c₁| = {result['c1_abs']:.10f}  (error: {result['err_abs']:.2e})")
# |c₁| = 1.0000000000  (error: 9.99e-14)
```

### Verify spin-½ from Hopf holonomy

```python
from geometrodynamics.hopf import compute_spinor_monodromy

result = compute_spinor_monodromy()
print(f"⟨ψ₀|U(2π)|ψ₀⟩ = {result['overlap_2pi']:.6f}  (should be −1)")
print(f"⟨ψ₀|U(4π)|ψ₀⟩ = {result['overlap_4pi']:.6f}  (should be +1)")
```

### Validate Coulomb law from eigenmode throat flux

```python
from geometrodynamics.tangherlini import solve_radial_modes, solve_maxwell_from_eigenmode

modes = {}
for l in [1, 3, 5]:
    oms, fns, rg = solve_radial_modes(l)
    modes[l] = {"omega": oms, "funcs": fns}

result = solve_maxwell_from_eigenmode(modes)
print(f"Q = {result['Q']:.6f}")
print(f"Relative error vs exact Coulomb: {result['rel_err']:.2e}")
```

### Reproduce the full charged-lepton ladder

```python
from geometrodynamics.tangherlini import (
    solved_lepton_masses_mev, S3_ACTION_BASE, TAU_BETA_50PI, tau_uplift_2pi_quanta,
)

masses = solved_lepton_masses_mev()   # locked baseline, no tuning
print(f"m_e  = {masses[0]:.6f} MeV")
print(f"m_mu = {masses[1]:.6f} MeV   (obs 105.658376)")
print(f"m_tau= {masses[2]:.6f} MeV   (obs 1776.860000)")

print(f"action_base = 2π         = {S3_ACTION_BASE:.6f}")
print(f"k_uplift β  = 50π        = {TAU_BETA_50PI:.6f}")
print(f"τ uplift    = 4β = 200π  = {tau_uplift_2pi_quanta(TAU_BETA_50PI):.0f} × (2π)")
```

### Run a QCD meson simulation

```python
import numpy as np
from geometrodynamics.qcd import make_meson_tube, HadronicNetworkSolver

net = make_meson_tube(L=1.0, v=1.0, N=100, dt=0.004)
s = np.linspace(0, 1.0, 100)
net.initialize_fields(psi0={0: 0.5 * np.sin(np.pi * s)})

solver = HadronicNetworkSolver(net, antipodal_coupling=0.05)
history = solver.run(n_steps=1000, record_every=50)
print(f"Energy drift: {np.std(history['energy']) / history['energy'][0]:.4f}")
```

### Build a black-hole condensate and verify entropy

```python
from geometrodynamics.blackhole import (
    build_schwarzschild_condensate, compute_entropy_balance,
    find_horizons, integrate_radial_geodesic,
)

# Schwarzschild BH as a coherent wormhole-throat condensate
bh = build_schwarzschild_condensate(mass=5.0)
bal = compute_entropy_balance(bh)
print(f"N throats: {bh.N}")
print(f"S_BH  = {bal.S_BH:.2f}")
print(f"S_thr = {bal.S_throat:.2f}  (relative error: {bal.relative_error:.2e})")
print(f"Net charge Q = {bh.net_charge}  (neutral)")

# Nonsingular interior: Hayward metric with core scale from throat network
l = bh.core_scale
horizons = find_horizons(bh.mass, l)
print(f"\nCore scale l = {l:.4f}")
print(f"Horizons: {['%.4f' % h for h in horizons]}")

# Geodesic completeness: infalling worldline decelerates, never hits r=0
geo = integrate_radial_geodesic(M=bh.mass, l=l, r_start=3*bh.mass, tau_max=100)
print(f"Geodesic complete: {geo.is_complete}  (r_min = {geo.r_min:.2e})")
```

## Lineage

This package refactors and unifies three monolithic scripts:

| Original file | Package modules |
|---|---|
| `geometrodynamics_v39.py` | `hopf/`, `tangherlini/`, `transaction/`, `constants.py` |
| `s3_spin2_closure_toy_solver_v22.py` | `tangherlini/` (shared spectral solver) |
| `qcd_topology_solver_v30.py` | `qcd/` (entire subpackage) |
| New in v0.41.0 | `blackhole/` (condensate, interior, entropy, derivation) |
| New in v0.42.0 | `embedding/`, `bell/`, `transaction/cavity.py` |
| New in v0.43.0 | `embedding/transport.py`, `bell/hopf_phases.py`, `history/` |
| New in v0.44.0 | `tangherlini/lepton_spectrum.py` (locked e/μ/τ ladder) + `scripts/` (calibration CLIs) |
| New in v0.45.0 | `qcd/quark_spectrum.py` + `qcd/hadron_spectrum.py` (shell-coupled six-quark ladder; residual sector geometrized to ~1% via Tangherlini eigenmode) |
| New in v0.46.0 | `experiments/closure_ledger/` (closure-ledger sequence; reduces the locked lepton surrogate's residual external input from six phenomenological parameters to one anchor m_e). Paper draft in `docs/hbar_origin_note.md`. |
| New in v0.47.0 | BAM effective-action scaffold (PRs #49–#53): five mismatch terms B1–B5; four closed (B1+B2 topological/discrete sector, B3 hard-wall BC, B5 master integral); B4 audited as irreducible-by-dimensional-necessity. Closure release note in `docs/scaffold_closure_release_note.md`. |
| New in v0.48.0 | Throat-as-anchor arc (PRs #55–#58): self-consistent equilibrium `R*`, cohesive brane tension `B·R²`, RS-like √6 brane tuning, pair threshold `2 m_e c²`. |
| New in v0.49.0 | Throat-as-relativistic-spin-½-particle arc (PRs #59–#62): moving-throat covariance, Hopf-holonomy Wigner rotation, `g = 2`, one-loop Schwinger `a = α/2π` reconstructed. |
| New in v0.50.0 | C / CPT / throat Dirac arc (PRs #63–#66): `C` = inner/outer swap, CPT on throat histories, explicit `Θ = −iγ⁵`, throat 4-spinor from `S_BAM` SUSY factorization. |
| New in v0.51.0 | Even-k absence + QCD shell arc (PRs #67–#69): spin-statistics classification of even-`k` absence, throat → QCD-shell transition, shell ↔ QCD structural match. |
| New in v0.52.0 | Three-generation / `k_5 = 5` arc (PRs #70–#74): sharp `k ≤ 5` boundary, `β_lepton = k_5²·(2π) = 50π`, `#gen = (k_5+1)/2 = 3`, `k_5 = D_bulk = dim(S³)+2 = 5`, `1/(2π)` in Schwinger anomaly = BAM closure-quantum loop measure (PR #74). |
| New in v0.53.0 | QCD-shell arc (PRs #76–#80): quark `n_part = 233` diagnosed as phenomenological compensator (PR #76, v3 lepton-shaped Hamiltonian is wrong machinery); quarks reframed as cavity wavefronts that resolve the shell with 6-state `(l, n, p)` basis + 6×6 operator scaffold (PR #77); shell mass-ordering / `n_part` audit identifies structural slots but not closure (PR #78); `χ_n` derived from cavity-mouth boundary stress (Z₂-antisymmetric piece, no free parameter; PR #79); BAM-native color algebra identified as `SU(2) × Z₂` from B2 + Hopf + PR #63 inner/outer swap (PR #80); inter-generation hierarchy outside BAM color scope, `n_part` remains residual compensator with sharply identified scope. |
| New in v0.54.0 | Pati-Salam bridge + mass-operator unification (PRs #82–#83): throat ↔ shell `n+3` Z₂ bridge unifying the lepton (throat) and quark (shell) sectors on a 12-state basis, with 3 open extensions identified for full SU(4) PS (PR #82); **the lepton `β·k²` and quark `ω²(l,n)` mass operators unified as one Bohr-Sommerfeld operator** `m² = (S/L_eff)²` with `L_throat = √(2π)/k_5` recovering `β_lepton = k_5²·(2π) = 50π`, `k = 0` for quarks = "don't pass through the throat", closure quanta `2π` (throat) vs `π` (cavity half-cycle) (PR #83). |

## License

MIT
