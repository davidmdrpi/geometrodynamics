[![DOI](https://zenodo.org/badge/1181274003.svg)](https://doi.org/10.5281/zenodo.20225786)
# Geometrodynamics

**A research framework implementing and testing Wheeler's geometrodynamic program.**

This package computationally explores the hypothesis that structures
physicists call electromagnetism, charge, spin, confinement, **black
holes**, and **Bell correlations** may emerge from the geometry of
spacetime itself — specifically the Hopf fibration on S³, 5D Tangherlini
wormholes, topological flux-tube networks, coherent wormhole-throat
condensates, and non-orientable throat topology.

## Where ℏ enters: scale-free closure ledger + one geometric anchor

The closure-ledger arc (`experiments/closure_ledger/`, PRs #11–#74)
reduces every dimensionless parameter in the locked lepton surrogate
to closure-quantum invariants (`action_base = 2π`, `transport = 8π`,
`resistance = 7π/100`, `pinhole γ = Σ V_max[1..5]`, `β_lepton = k_5²·(2π) = 50π`,
`ε = 7π/(100·k_5⁴)`), and an audit (`maslov_dimensional_bridge_probe`,
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
| R_OUTER selected by cross-species fixed point | **Verified** | Bisection on each lepton gives same R* ≈ 1.262 to 0.008 % |
| Pinhole γ ≈ Σ V_max[1..5] on Chebyshev grid | **Verified** | −2.2 % off the locked γ = 22.5; same operator as the QCD-sector γ_q |
| Transport = 8π = 4·(2π) | **Verified** | +0.13 % off the locked transport = 25.1; 4th closure quantum |
| Resistance = 7π / 100 | **Verified** | +0.94 % off the locked resistance = 0.2179; selected over `4·(ω−1)` by R_OUTER bisection |
| Inner cutoff `ε = resistance / k_5⁴` | **Verified** | Closes the Compton bridge `ℏ = m_e R_MID c` to 0.04 % |
| Closure-quantum ledger closes modulo m_e | **Established** | Every locked parameter is a closure-quantum invariant; m_e is the unique remaining external input |
| Quark mass ladder (u, d, s, c, b, t) | **Fitted** | 1.6% max rel err on s, c, b, t with d-anchor, four shell-index axioms, and one phenomenological β |
| Quark shell-index axioms (ε, η, χ, phase) | **Geometric** | All four expressible in `k_5 = 5` only: `(1−1/k_5², k_5, (k_5−1)·k_5, 0)` |
| Quark residual sector (transport, pinhole, resistance) | **Derived** | Each matches Tangherlini eigenmode quantity within ~1% on the tortoise grid |
| Pinhole = `Σ V_max(l=1..5)` (tortoise grid) | **Verified** | −1.09% off the fitted lock |
| Transport = `mean ⟨u_l\|V_{l+2}−V_l\|u_{l+2}⟩` | **Verified** | +0.87% off the fitted lock |
| Resistance = `transport · ln(α_q(k_5)/α_q(k_1))` | **Verified** | −0.43% off the fitted lock |
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
| Glueballs: pure-confinement benchmark + Möbius tower | **Benchmark + topological prediction** | Closed flux loops (no valence quarks ⟹ no flavor puzzle) are the cleanest confinement probe. BAM orientable ground `√(4πσ)≈1.50 GeV` (3.5√σ) benchmarks lattice 0++ (4.1√σ) to ~13%; closed-string glueball Regge slope = half the meson. **BAM-specific:** the non-orientable **Möbius** sector (`make_mobius_tube`, antiperiodic) gives an *extra* glueball tower (half-integer modes, shifted `+πσ` in `M²`) interleaving the orientable one — ≈2× the states. Glueballs are **not experimentally observed**, so this topological divergence is testable against lattice, not contradicted by experiment (`glueball_closed_flux_loop_probe`, PR #100) |
| Möbius flux tube ⟹ exotic `J^PC`; observed hybrids match | **Matches data** | Flux-network topology = hadron taxonomy (meson/baryon/tetraquark/pentaquark/hybrid/glueball + Möbius Z₂). A **non-orientable (Möbius) flux tube** carries the antiperiodic phonon that opens the **exotic `1-+`** (forbidden to ordinary qq̄: `P=(−1)^{L+1}`, `C=(−1)^{L+S}`). The observed exotic hybrids `π₁(1600)`, `η₁(1855)` (both `1-+`) match at the right `J^PC` and at `ρ/ω + 2√σ ≈ 1.62, 1.85 GeV`; the tetraquarks (`X, Z_c, T_cc`) / pentaquarks (`P_c`) fit multi-junction networks. **Unlike glueballs, exotics are observed** — so this is where BAM's non-orientable topology meets data, and matches (`mobius_exotic_sector_probe`, PR #101) |
| BAM baryonic exotics: classification + constraints | **Most-constrained corner** | Unlike mesons (smoking-gun `1-+`), **baryons have no forbidden `J^P`** (`P=(−1)^L`, `S∈{½,3/2}`, no `C`) — so BAM's Möbius/hybrid baryons are **supernumerary ordinary-`J^P`** states, identifiable only by counting. They sit in the light N*/Δ* region (`nucleon/Δ + 2√σ ≈ 1.79, 2.08 GeV`), the densest, best-measured spectrum — the **most experimentally constrained** corner of BAM's non-orientable predictions (opposite extreme from glueballs). The Möbius doubling must coincide with observed resonances or decouple (`πN`), else be excluded. Constraint ranking: light N*/Δ* > strange hyperons > charm/bottom baryons (freest) (`baryonic_exotics_classification_probe`, PR #102) |
| Heavy-quark Möbius baryon: prediction in the freest channel | **Findable / unconstrained** | By heavy-quark symmetry (heavy quark = spectator) the Möbius/flux gap `Δ=2√σ≈0.85 GeV` is **flavor-independent** (same for c and b) — the cross-flavor signature replacing the absent exotic-`J^P`. Predictions: Λ_c **~3.14**, Ω_c ~3.54, Λ_b **~6.47**, Ω_b ~6.89, Ξ_cc ~4.47 GeV — all just **above** current excitation ceilings (findable at LHCb/Belle II, not excluded) and above the orbital tower. Doubly-heavy `Ξ_cc` and `Ω_b` have no measured excitations → entirely unconstrained. Exact mass (lattice hybrid gap 0.8–1.3 GeV) / `J^P` open (`heavy_mobius_baryon_probe`, PR #103) |
| Heavy Möbius baryon: decay channels + search strategy | **Twist-unwinding → hybrid selection rule (falsifiable)** | Completes #103: how the state decays and how to find it. Decay = **twist-unwinding** (non-orientable `−1` → orientable `+1` ground state sheds `2√σ` as light isoscalar hadrons), so it inherits the flux-tube **hybrid selection rule**: single-S-wave-π-to-ground **SUPPRESSED**; `Σ_Q π` / isoscalar dipion `Λ_Q(ππ)` / P-wave+π **PREFERRED** — the branching **pattern** that distinguishes it from a radial excitation (which does the opposite). Cross-flavor clincher: all-light Q-values **identical** for c and b (`Λ_Q ππ` **569**, `Λ_Q η` **301** MeV; `Σ_Q π` offset only by hyperfine 167/194). Broad (~tens–150 MeV, open channels) → best in LHCb/Belle II amplitude analyses of `Λ_Q ππ`, `Σ_Q π`, `DN`/`BN` (`Ξ_cc`/`Ω_b` wide open). Branching fractions / width / `J^P` open (`heavy_mobius_baryon_decay_probe`, PR #109) |
| Non-orientable sector: compact **experimental note** | **Compiled (reference card)** | Consolidates the whole Möbius / closed-flux-loop sector (PRs #100–#109) into one LHCb/Belle II/BESIII-style note — predicted masses, Q-values, preferred/suppressed modes, analysis handles — every number a pushforward of the single input `√σ`. **Masses:** mesonic `1⁻⁺` π₁ **~1.62**, η₁ **~1.85** GeV (matched to π₁(1600)/η₁(1855)); glueball `0⁺⁺` `√(4πσ)` ~1.50 GeV (unobserved, freest); heavy Möbius baryons Λ_c 3135 … Ω_b 6894 MeV. **Decays:** twist-unwinding → hybrid selection rule (single-π-to-ground suppressed), cross-flavor Q-match (`Λ_Q ππ` 569, `Λ_Q η` 301 MeV identical c=b). **Handles:** branching pattern vs radial, isoscalar high-`m(ππ)` dipion, broad→amplitude fits, `1⁻⁺` smoking gun (mesons). Standalone at `docs/bam_nonorientable_experimental_note.md` (`nonorientable_experimental_note_probe`, PR #110) |
| Heavy Möbius baryon: sharper **LHCb / Belle II search table** | **Tiered, actionable** | Converts #109/#110 into a ranked search table. **New handle:** the `Λ_Q(ππ)` **dipion endpoint** `m(ππ)_max = M_Möbius − M_ground = 2√σ ≈ 849 MeV` is **flavor-independent** (same edge above charm and bottom, peaking high) — a fixed edge in a directly-plotted observable, one overlay tests the framework. **Tier 1** (discovery pair): Λ_c (3135, `Λ_c⁺π⁺π⁻`, `Λ_c⁺→pK⁻π⁺`, LHCb+Belle II) + Λ_b (6469, `Λ_b⁰π⁺π⁻`, LHCb b-decays) — the cross-flavor clincher. **Tier 2** (unexplored, rare): Ξ_cc (4471, `Ξ_cc⁺⁺→Λ_c⁺K⁻π⁺π⁺`), Ω_b (6894). **Tier 3** (calibratable): Ω_c (3544, above 2017 excitations). Discriminators: suppressed single-π-to-ground, 849 MeV endpoint, cross-flavor Q-match. Standalone at `docs/heavy_mobius_baryon_search_table.md`; masses ±band / broad / BFs / `J^P` open (`heavy_mobius_baryon_search_table_probe`, PR #114) |
| **Program-wide synthesis: the input budget** | **Capstone** | Classifies every result into 5 epistemic tiers. **The whole dimensionful content reduces to 2 B4 anchors** — `m_e = ℏc/R_MID` (QED/lepton) and `√σ ≈ Λ_QCD` (confinement) — the irreducible minimum (one scale/sector, PR #52). Open dimensionless inputs are localized to 2 (neutrino compliance `ε`, quark `n_part`); the only other open input is the **universal flavor puzzle** (Yukawa hierarchy — not BAM-specific). The rest is ~22 derived-geometry results + 6 non-orientable topological predictions (matched → falsifiable → findable → free) (`program_synthesis_probe`, PR #104) |
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
│   └── viz/                  # Visualisation (placeholder)
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
