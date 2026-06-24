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
| **The odd-k ladder** — forced, rigid, unique to the non-orientable 5D geometry | **The cleanest discrete feature audited: the odd-k lepton ladder {1,3,5} is forced (fermion ⟹ odd via the T²=−I orientability grading; k≤5 ⟹ 3 generations), rigid against every active/null/mixed continuous deformation (the rank story survives nonlinearity), and a unique signature of the non-orientable antipodal 5D geometry** | Continues the #173 inverse problem on one discrete feature, asking whether the geometry forces it and whether anything but this geometry could produce it. **The origin** (recap #67/#169/#170): the throat monodromy `T=iσ_y` (`T²=−I`) makes `T^k` off-diagonal for odd k (the orientation-reversing, non-orientable RP²/Pin⁻ closure of a spin-½ fermion) and diagonal for even k (orientable, bosonic) — so `k mod 2` is the orientability grading; charged leptons are fermions ⟹ odd k, and `k≤k_5=D_bulk=5` ⟹ {1,3,5}=(k_5+1)/2=3 generations. **The three direction sets** (from the #173 Jacobian): *active* (the rank-10 subspace) moves the masses/CKM **linearly** (exponent 1.03); *null* (the 10-dim kernel) is **flat to first order** (exponent 2.0, ~10⁴× smaller response than active at ε=10⁻²); *mixed* (active+null) is **active-dominated** (1.02) — so **nonlinear effects do not break the local rank story** (the null leakage stays quadratic, not linear). **Forced & rigid:** the odd-k labels and the generation count are integer winding + the ℤ₂ grading (`T²=−I`), discrete topological data *outside* the entire continuous deformation manifold (no generation-number knob) — not an emergent near-integer that could drift, but structurally forced; the continuous geometry deforms only the masses and the CKM, the ladder is rigid against all of it. **Unique:** an orientable geometry (`T²=+I`) gives the even/bosonic sector, not an odd-only fermion ladder; the specific {1,3,5} needs `k≤5=D_bulk` — so odd-{1,3,5} is the joint signature of the non-orientable antipodal spin structure and the 5D bulk (an exclusion/signature argument within BAM, not a no-go against all alternatives) (`odd_k_ladder_rigidity_probe`, PR #174) |
| **The sensitivity audit** — Jacobian rank, the forced core, the isolation dimension | **The dynamical inverse problem, measured: of 14 live observables the free geometry dials rank 10; the FORCED CORE is 4 (entirely CKM unitarity, structural V=U₊†U₋); the masses are FITTED; a 10-dim compensator redundancy (the diagonal shifts); deriving φ_h saves 0 effective inputs** | Lets continuous geometry search the spectrum/signature landscape instead of adding more static proofs: compute the Jacobian `J_ij=∂O_i/∂I_j` of the live observables at the lock and read its SVD. **Observables** (14): 4 quark mass ratios (s,c,b,t/d), the 5 CKM magnitudes, J, β, γ (`LOCKED_QUARK_PARAMS_V4`), and the 2 lepton ratios (μ/e, τ/e). **Inputs:** the free *fitted* knobs — NOT the k₅-derived locks (φ_h=π/k₅, χ=k₅(k₅−1), uplift=1−1/k₅², action=π, winding=max), which are zero-cost. **Isolation dimension** rank(J)=10 (quark 8 + lepton 2), with a clean singular-value gap (quark sv falls 22.6→1.2e-2 then drops ~5e-6 to zero). **Forced core** = n_obs−rank = **4**, entirely CKM combinations (forced weight on CKM >0.99, on masses ~1e-16) — the **CKM unitarity relations**: `V=U₊†U₋` is exactly unitary (‖V†V−I‖~1e-16), so the 8 CKM observables on the 4-parameter unitary manifold force 8−4=4 relations. The largest set forced at zero input cost — a genuine structural prediction, but the *standard* unitarity, not a BAM-specific relation. **The masses are FITTED** (zero weight in the forced core, quark and lepton): the ladder sets their values, the knobs span them; no forced mass relation. **Compensator redundancy** = n_inputs−rank = **10**, dominated by the mass-preserving diagonal shifts (kernel share ~0.8) — the *loose knob / compensator* structure the program flagged qualitatively, now measured; the v4 quark parametrization is substantially over-complete. **CP-at-zero-cost test:** adding φ_h as an input leaves the rank unchanged (8→8) — the CP-phase direction is already spanned, so deriving φ_h saves no effective input; *"CP at zero parameters"* is a counting economy, not a Jacobian reduction. **Robustness:** the rank is stable across a 4×3 grid of finite-difference step (1e-3…1e-6) and SVD tolerance (1e-4…1e-8) — total rank 10 in 11/12 cells (the lone outlier is the largest step at the tightest tolerance, truncation noise crossing the cut), with the rank-8 separation gap ≥1e3 everywhere, so the forced core of 4 is not a tuning of h or tol. An audit — a measurement of the predictive content, honest where it is not flattering (`sensitivity_jacobian_audit_probe`, PR #173) |
| **The measured Fermi equation of state** — a many-throat ensemble (companion to #171, same branch) | **The EoS *measured*, not assumed: free fermions in a box with the −1 exchange sign as Pauli single-occupancy give Γ=5/3 (NR, measured 1.6665) and 4/3 (UR, 1.3332) and P/u=2/3, 1/3; the Bose control gives Γ=1 with vanishing degeneracy pressure — three routes (assumed #170 / topological #171 / measured #172) agree** | The second of the two closing options, kept on the same branch as #171 to compare. #170 *assumed* antisymmetry and read 5/3 off the analytic Fermi integral; #171 derived the −1 sign topologically; here the −1 is imposed as **Pauli single-occupancy** in a many-throat box ensemble and the EoS is **measured** from the simulated ground-state energies. **The ensemble:** N identical throats as free fermions in a cubic box; each spatial mode `(n_x,n_y,n_z)` holds `g=2`; the many-body ground state is the filled Fermi sea built by filling the N lowest of ~1.3M enumerated modes — no occupation distribution assumed. **Measured P/u:** from `P=−dE/dV`, `0.6667=2/3` (NR) and `0.3333=1/3` (UR) — the virial relation emerges. **Measured Γ:** the local log-slope of the filled-mode energy sum `K(N)`, finite-size-extrapolated via the Weyl correction `Γ(N)=Γ∞−a·N^{−1/3}`, gives **1.6665≈5/3** (NR) and **1.3332≈4/3** (UR), 0.01% from target — outputs of the simulation, not a formula. **The control:** Bose statistics (all N in the ground mode) give `Γ=1` and a `T=0` degeneracy pressure that vanishes (mode energy `∝1/L²→0`) — so the 5/3 stiffening is a *measured consequence of the −1 exchange sign*, not a property of the box. **Comparison:** the measured indices reproduce #170's assumed values and confirm the EoS implied by #171's topological sign — three routes, one answer. **Honest scope:** the exchange sign is the input (Pin⁻/geon, #170/#171, as Pauli occupancy), not re-derived here; idealizations are the standard degenerate-gas ones (free, T=0, box) (`measured_fermi_eos_ensemble_probe`, PR #172) |
| **π₁ of the two-mouth configuration space** — FR-homotopy survival for the Pin⁻ throat (geon statistics) | **Replaces #170's *orientable* FR citation with the geon-statistics framework: π₁ exchange σ²=e (no 3D braiding), spinorial 2π=−1, Pin⁻ reflection²=−1 (the ingredient orientable FR lacks), achiral → the −1 (Fermi) survives, conditional on the cited Dowker–Sorkin exchangeability hypothesis** | Closes the cited gap in #170: the Finkelstein–Rubinstein homotopy it invoked is the **orientable** result, but the throat mouth is the non-orientable RP² (Pin⁻), so the correct framework is **geon statistics** (Friedman–Sorkin; Aneziris–Balachandran et al.; Dowker–Sorkin), where spin–statistics is a theorem *with hypotheses, known to fail for some geons*. **π₁ of the two-mouth config space** (model): exchange `σ` with `σ²=e` (≥3D ⇒ symmetric group, **no braiding**, only ±1 statistics), per-geon 2π rotation `R_i`, and — because the mouth is non-orientable — an orientation-reversing loop `τ_i`. **Spinorial:** the single geon's 2π rotation = −I (4π=+I), the Pin⁻ holonomy of #170 / Friedman–Sorkin's spin-½ from gravity. **The new ingredient orientable FR lacks:** the non-orientable exchange carries an orientation reversal — a *reflection* — and RP² admits **Pin⁻ only** (`w₂≠0` kills Pin⁺), in which a reflection **squares to −1** (computed; Pin⁺ gives +1) — exactly what makes the non-orientable exchange sign well-defined and fermionic. **Achirality:** a non-orientable geon is its own mirror image, meeting the geon spin–statistics theorem's handedness hypothesis automatically. **Result:** spinorial + Pin⁻ reflection²=−1 + achiral ⇒ exchange = −1 (**Fermi survives**), now on the correct non-orientable footing. **Honest scope:** conditional on the Dowker–Sorkin exchangeability ("slide") hypothesis — holds for identical asymptotically-flat throats, *cited not derived* from the BAM field theory; the literature has spin–statistics violation examples, so this is a genuine check BAM passes (Pin⁻ + achiral + exchangeable), not automatic. The remaining gap is that hypothesis (and the field-theory mapping class group, modeled here), not the spinor sign or the reflection algebra (`geon_statistics_pi1_probe`, PR #171) |
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
| Glueballs: pure-confinement benchmark + Möbius tower | **Benchmark + topological prediction** | Closed flux loops (no valence quarks ⟹ no flavor puzzle) are the cleanest confinement probe. BAM orientable ground `√(4πσ)≈1.50 GeV` (3.5√σ) benchmarks lattice 0++ (4.1√σ) to ~13%; closed-string glueball Regge slope = half the meson. **BAM-specific:** the non-orientable **Möbius** sector (`make_mobius_tube`, antiperiodic) gives an *extra* glueball tower (half-integer modes, shifted `+πσ` in `M²`) interleaving the orientable one — ≈2× the states. Glueballs are **not experimentally observed**, so this topological divergence is testable against lattice, not contradicted by experiment (`glueball_closed_flux_loop_probe`, PR #100) |
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
