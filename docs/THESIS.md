# Bulk Antipodal Mechanics — Thesis

This document describes the foundational thesis of **Bulk Antipodal
Mechanics (BAM)**, a research program in classical-geometric physics on
a closed `S³` universe.

The name foregrounds the channels by which discrete quantum-mechanical
spectra arise from continuous, GR-consistent geometry. *Bulk* refers to
the 5D Tangherlini extension and the non-orientable wormhole throats
that live in it; *antipodal* refers to `S³` closure with antipodal wave
focusing; *mechanics* signals a quantitative framework for computing
observables, in the same family of usage as quantum mechanics or
statistical mechanics. BAM is a descendant and extension of Wheeler's
geometrodynamics rather than a simple revival of it: it is what you get
when the global and topological machinery Wheeler did not have is added
to the geometrodynamic instinct that *matter is a property of
spacetime, not an independent field attached to it*.

One framing point should be stated at the outset, because it is easy to invert.
**BAM derives quantum field theory from continuous, classical general
relativity; it is the opposite of a quantum-gravity program, and it does not
quantise gravity.** The foundational layer is a *classical* GR geometry — the
`S³` slice, the wormhole throat, the 5D Tangherlini bulk and its metric `f(r)` —
and the quantum structure (the discrete matter spectrum, the propagator and
exchange kernel, the self-energy, the interaction vertices) is *reconstructed
on that fixed classical background*, in the precise sense of quantum field
theory on a curved spacetime. The arrow runs geometry → fields and never the
reverse: the metric is a classical input, not a quantised dynamical variable.
Consequently, asking BAM to address quantum gravity is a category error that
would turn the program upside down — gravity is the foundational classical layer
*from which* quantum matter is derived, not an object the program seeks to
quantise. Where this thesis later speaks of a path-integral measure `S_BAM`, a
one-loop fluctuation determinant, or a bounded interacting vacuum, those are
statements about the emergent *matter* field theory, read off the classical
throat geometry — not about a quantum theory of the metric.

## The BAM conjecture

> Quantum particles are self-consistent topological boundary conditions on
> a closed spatial geometry — stable when at rest, and persisting under
> motion via time-dependent throat transport.

This package is a computational investigation of that conjecture in the
specific setting of an `S³` spatial slice carrying the Hopf fibration,
with non-orientable wormhole throats and a 5D Tangherlini bulk. The new
ingredient — and the reason BAM can make quantitative progress where
Wheeler's program stalled — is a set of global and topological tools
that were not available in the 1960s.

## What BAM is trying to demonstrate

The aim is to test, computationally, whether quantum-mechanical
observables can be computed from purely geometric inputs on a closed
`S³`, without canonical quantization as an independent postulate. BAM
does **not** claim to have shown this in full. It claims to have
demonstrated that several of the most distinctive observables — charge
quantization, spin-½, the Coulomb radial response, the charged-lepton
mass ladder, the six-quark mass ladder, Bell's inequality saturation,
regular black-hole interiors — can each be derived from, or made
consistent with, classical geometry on `S³` alone. (The
finite-separation Coulomb claim — that two throat mouths on `S³`
produce `F(ψ) ∝ 1/sin²(ψ)` — has now been **demonstrated**: the S³
Green response reproduces the Coulomb potential and inverse-square
force in the flat-space limit, with the `1/sin²ψ` form confirmed as
the leading behaviour and refined by a compact-`S³` antipodal-image
modulation; see `docs/two_throat_coulomb_research_plan.md` and the
`two_throat_coulomb_probe`.) Whether this set is closeable into a
complete derivation of quantum field theory is the open question BAM
exists to test.

A reader looking for "BAM replaces QFT" should read the program as a
falsification campaign for that claim. A reader looking for "geometric
structures that happen to reproduce quantum observables" should read it
as evidence that the geometric channels BAM identifies are real,
whatever the final relationship to QFT turns out to be.

## Why this direction — GR-foundational rather than QM-foundational

Most quantum-gravity research takes quantum mechanics as foundational
and asks how to quantize gravity. String theory, loop quantum gravity,
asymptotic safety, causal dynamical triangulations, and the various
holographic programs all share this orientation: the quantum is given,
gravity must be made compatible with it. BAM goes in the opposite
direction. It treats classical general relativity as foundational and
asks whether quantum-mechanical observables emerge from it.

There is a macro-scale intuition that motivates the inversion. The only
robust classical GR structures we see at astrophysical scales are a
short list: gravitational waves, horizons, black holes, topology, and
global causal structure. That vocabulary is already enough to describe
nearly everything we observe in the strong-field regime. If a
microscopic "particle" is *not* a new category of object — if it is a
small, resonant, topologically constrained version of the same
gravitational vocabulary — then quantum mechanics is what the resonance
spectrum of that vocabulary looks like at small scales. The microscopic
world is then not a different physics added on top of GR; it is the
same physics, run on a closed compact slice, where antipodal closure
and bulk confinement force the spectrum to be discrete.

This is a bet, not a theorem. BAM bets that:

- *Particle* ≈ a topologically constrained, resonant configuration of
  the same wave / horizon / throat machinery that produces black holes
  and gravitational waves at large scales.
- *Charge* ≈ a topological winding number of the geometry, not a
  separate field attached to it.
- *Spin* ≈ holonomy of the geometry's natural connection.
- *Quantization* ≈ closure on a compact spatial slice.

Each identification is testable. The validation table records which
checks have passed; the falsification tests below record which ones are
next. If the bet is right, quantum mechanics is not foundational; it is
the small-scale resonance theory of GR on a closed universe. If the bet
is wrong, it will be wrong in a specific, diagnosable way — which
channel was overcredited, which mechanism failed to compose — and a
corrected program may still be available.

What makes the bet worth taking seriously now, beyond Wheeler, is that
the topological machinery this requires has already proven productive
elsewhere. Chern numbers, holonomy of non-orientable bundles, and
spectral methods on finite intervals are the standard toolkit of
topological condensed matter, where they reproduce robust quantization
phenomena from continuous geometric inputs. BAM's proposal is that the
same toolkit, applied to `S³` rather than to a band structure, has an
unexploited use in fundamental physics: not as a calculation device
layered onto an independently quantized substrate, but as the
geometric vocabulary that *is* the quantum description, when run on a
closed universe.

## Why BAM is a research program, not a theory

Three things distinguish BAM from a theory-of-everything announcement.

First, the program is **explicit about its open problems** and lists
them in the validation table with status flags: *Verified*, *Derived*,
*Modeled*, *Constructed*, *Fitted*, *Phenomenological*. A claim labelled
`Modeled` is not the same as a claim labelled `Derived`, and the README
is uniform about the distinction.

Second, BAM is **falsifiable in concrete tests**. Moving-throat Berry
phase, even-`k` closure absence, and quark `β` derivation are tests
where a clean negative result would invalidate parts of the program;
these are listed below and have not yet been run. The **two-throat
Coulomb force on `S³`** has now been run (`two_throat_coulomb_probe`):
it reproduces the Coulomb potential and inverse-square force in the
flat-space limit and confirms the `F(ψ) ∝ 1/sin²(ψ)` leading form, so
this particular falsification test is **passed** rather than open.

Third, BAM **rests on classical geometric machinery only**. No canonical
commutators are imposed, no Hilbert space is assumed at the outset, and
no path integral is performed. Any recovery of quantum-like observables
therefore constrains the geometry, not a hidden quantization step.
Where Planck's constant enters the quantitative comparisons (lepton and
quark masses in MeV, for example) is itself an open question called out
in the validation table — the absolute MeV scale is set by anchoring
the electron mass; the *ratios* are derived.

A compact status map, expressed at the level of *claim classes* rather
than per-test rows (the README's validation table covers the per-test
detail):

| Claim class | Current status |
|---|---|
| Hopf charge / Chern structure | exact geometric identities |
| Spinor transport from `T = iσ_y` | verified throat-orientation structure |
| Bell correlations / CHSH | derived from throat transport; `2√2` verified |
| Lepton mass ladder | locked spectral model; `β_lepton = k_5²·(2π) = 50π` now derived from the topological charge (PR #71) |
| Quark mass ladder | 1.6% fitted ladder; residual sector geometrized; quark `β = 233π` (with `n_part = 233`) diagnosed as phenomenological compensator absorbing the inter-generation mass hierarchy (PR #76); shell-waveguide arc PRs #77–#80 reframes the sector structurally |
| Quark sector reframed as cavity wavefronts | PR #76 diagnosis + four-PR shell arc (PRs #77–#80): quarks are the shell-saturated wavefronts that resolve the cavity (not throat traversals like leptons); 6-state `(l, n, p)` basis; `χ_n` derived from cavity-mouth boundary stress (no free parameter, PR #79); BAM-native color algebra `SU(2) × Z₂` from B2 + Hopf + PR #63 (PR #80) |
| QCD color algebra | BAM-native = `SU(2) × Z₂` from established primitives (PR #80); standard SU(3) NOT derivable from current scaffold — natural triplet candidates all give SO(3)/SU(2); Pati-Salam SU(4) extension (with throat↔shell algebra map) is the most plausible route to SU(3), genuine open work |
| **Lepton + quark mass operators unified** | The lepton `β·k²` (PR #71) and quark `ω²(l, n)` (PR #77) mass operators are **one Bohr-Sommerfeld operator** `m² = (S/L_eff)²` (PR #83): `m²(k,n) = (k·2π/L_throat)² + ((n+1)·π/L_cavity)²`, `L_throat = √(2π)/k_5`. Leptons wind (`k ∈ {1,3,5}`); quarks resolve the cavity (`k = 0`). Cavity Bohr-Sommerfeld verified to machine precision; `β_lepton = k_5²·(2π)` recovered. Throat↔shell `n+3` Pati-Salam bridge built (PR #82); inter-generation hierarchy still open |
| Coulomb radial response | verified by Tangherlini/Maxwell BVP |
| Coulomb force at finite separation | verified (`two_throat_coulomb_probe`) |
| Black-hole interior / entropy | regular metric derived; entropy currently a consistency check |
| Compton tree amplitude (Klein-Nishina) | reproduced exactly via closed-form F² (PRs #25–#35); F² = K²·Q from one C×S³ master functional (B5 closed, PR #51) |
| Tree QED (BW, annihilation, Bhabha, Møller) | reproduced from BAM-geometric primitives (PRs #36–#46); see `docs/tree_qed_status.md` |
| BAM effective-action scaffold (B1–B5) | four closed (B1+B2, B3, B5); B4 audited as irreducible-by-dimensional-necessity (PRs #49–#53); `docs/bam_scaffold_status.md` |
| Throat as finite-self-energy equilibrium | derived `R* = (A/2B)^{1/3}` (PR #55); cohesive `B·R²` = brane tension (PR #56); bulk-gravity tuning √6 (PR #57); pair threshold `2m_e c²` (PR #58) |
| Throat = relativistic spin-½ particle | dispersion `E²−(pc)²=(mc²)²` (PR #59), Hopf-holonomy Wigner rotation (PR #60), `g = 2` from Pauli + Hopf monopole (PR #61), Schwinger `a = α/2π` reconstructed (PR #62), `1/(2π)` = BAM closure quantum (PR #74) |
| C / CPT / throat Dirac spinor | `C` = inner/outer swap `c₁ → −c₁` (PR #63), CPT on throat histories (PR #64), explicit `Θ = −iγ⁵` on throat spinor (PR #65), throat 4-spinor from `S_BAM` SUSY (PR #66) |
| Even-`k` absence → QCD shell | classified as spin-statistics selection rule (PR #67); higher-`k` excitations transition into QCD shell channel (PR #68); shell ↔ QCD structural match (PR #69) |
| Three generations / `k_5 = 5` | sharp `k ≤ 5` boundary (PR #70), `β_lepton = k_5²·(2π)` (PR #71), `#gen = (k_5+1)/2 = 3` (PR #72), `k_5 = D_bulk = dim(S³)+2 = 5` (PR #73) |
| `ℏ` origin | B4 audited (#52): closure-ledger machinery scale-free, so exactly one external dimensionful anchor required; relocatable to invariant bulk separation `ΔR` (#53), giving `m_e = 0.52·ℏ/(ΔR·c)`. Predicting ℏ in SI is gated solely by the value of that single geometric anchor |
| Full QFT / loop measure | `1/(2π)` in Schwinger anomaly identified as BAM closure quantum (PR #74); the full `S_BAM` path-integral measure is then **structurally constructed** as a loop-measure — a closure-ledger sector sum over a `Diff(S¹)`-gauge-fixed loop-space integral, with odd-k upgraded to the `Z₂` orientation-anomaly condition and the bounces as leading saddle (PR #115) — but its **analytic core remains open**: the bare fluctuation determinant diverges, so the normalisation needs regularization and is not yet rigorously constructed |

## Three mechanisms that compose

Wheeler's original geometrodynamics had the right instinct but lacked
the global machinery to make discrete spectra count anything. The
continuum Einstein equations admit far too many solutions; "charge
without charge" and "mass without mass" remained slogans because there
was no mechanism by which a continuous theory could pick out a discrete
spectrum. BAM proposes that three independent topological and geometric
channels each contribute discreteness, and that they **compose** rather
than competing. The name foregrounds two geometric arenas — the
antipodal `S³` cavity and the higher-dimensional bulk — but the
discreteness mechanism has three parts: antipodal closure, non-orientable
throat/shell transport, and finite bulk radial confinement.

**1. Antipodal `S³` closure.** Compactifying the spatial slice as `S³`
replaces the open continuum with a closed cavity. Any field that closes
on itself does so over a great circle of fixed length `2π`. Resonance
on a closed cavity is intrinsically discrete. Some closure constants
are direct geometric invariants — `action_base = 2π` is the
great-circle circumference of the cavity. Others, such as the
lepton-sector integer-winding lock `4β = 100·(2π)`, are sharp spectral
regularities identified by the solver; deriving why the multiplier is
exactly 100 remains an open analytic problem.

**2. Non-orientable throat / shell spectra.** A wormhole throat that is
non-orientable carries a `Z₂` partition class `p = ±` — a real
topological label, not a continuous parameter. The unique
orientation-reversing isometry of `S³` that preserves the Hopf bundle
is `T = iσ_y`, derived without ansatz in
`embedding/transport.py`. The identity `T² = −I` is the 4π
periodicity of spinors. The partition splitting drives every
mass-ordering inversion in the shelled sector (the
`m_u < m_d` but `m_c > m_s` pattern). Throat orientation is the
mechanism by which spin-½ behavior may become unavoidable rather than
imposed; the static holonomy result is verified, while the dynamic
moving-mouth version is the Berry-phase falsification test below.

**3. Uniform bulk distance from outer to inner.** The Tangherlini
throat confines a radial coordinate to the finite shell
`[R_INNER, R_OUTER]`. In tortoise coordinates this is a finite interval
with regular boundary conditions, producing a discrete eigenmode
spectrum independent of the `S³` closure but composing with it.

The new claim is that **these three channels compose into a quantitative
spectrum**. The lepton ladder is a "minimal closure" spectrum where
channel 1 dominates: each lepton mass scales with its global pass-count
winding `β·k²` on a nearly bare closure skeleton, locked by
`4β_lepton = 100·(2π)`. The quark ladder is a "shell-coupled closure"
spectrum where channel 1 picks up the heaviest shell only and channels 2
and 3 — partition asymmetry on the throat and bulk-mode coupling —
determine the lighter shells. Three of the four quark-sector residuals
are derivable from the existing eigensolver on the same tortoise grid
to within ~1%.

As of PR #83 these two ladders are recognized as **one Bohr-Sommerfeld
mass operator** read in two channels, not two separate spectra. The
unified operator is

```
m²(k, n)  =  (k·2π / L_throat)²  +  ((n+1)·π / L_cavity)²
```

with `L_throat = √(2π)/k_5`. The first term is channel 1 (throat
winding, closure quantum `2π`); the second is channel 3 (radial cavity,
Bohr-Sommerfeld half-cycle `π`). **Leptons wind through the throat**
(`k ∈ {1,3,5}`, lowest radial mode `n=0`) so the winding term dominates
and `m² ≈ β·k²`, recovering `β_lepton = (2π/L_throat)² = k_5²·(2π) = 50π`
(PR #71). **Quarks resolve the cavity** (`k = 0`, radial overtones
`n ∈ {3,4,5}`) so the winding term vanishes and `m² ≈ ω²(l, n)`, with
the cavity eigenvalues verified to be Bohr-Sommerfeld
(`∮√(ω²−V) dr* = (n+1)·π` to machine precision). The throat-traversal /
cavity-resolution dichotomy is the single quantum number `k`, and the
`2π`-vs-`π` channel quanta are the program's pervasive full/half-cycle
distinction. The two channels are exactly the closure ledger's
`N_total = N_layer1 + N_radial` (the B4 Maslov audit).

**A boundary on the unified operator (PR #165).** The two channels share
one S³ radius `R` — "everything rides on one R". A Berger-sphere
deformation audit (squashing the Hopf fiber by `λ` while keeping the base
round — the one move that separates the throat scale from the cavity
scale) maps where that shorthand breaks. The global cosmic-cavity Casimir
energy `E_cav(λ)` (zeta-regularized conformal scalar on the genuine SU(2)
Berger spectrum, validated at `λ=1` against the exact `1/240R`) and the
local throat self-energy `λ_min(λ)` both vary with `λ`, but **differently**:
the parameter-free ratio `ρ(λ) = E_cav/E_self` is not flat, so the two are
not one dynamical object even in shape. Decisively, `ρ(1) ≈ 3.3·10⁻⁴` (the
geometric one-R prediction) sits **~35 orders of magnitude** off the
measured global/local ratio `λ_C/R_Hubble ≈ 3·10⁻³⁹` — the cosmic cavity
and the local throat cannot ride on one `R` (the cosmological-constant
problem, in geometric form). R-unification is therefore a valid
**scale-free bookkeeping** device — fully consistent with the B4 audit
that the machinery is dimensionless and needs exactly one external anchor —
but **not** a physical single-`R` identification. *(Honest negative result;
the audit enforced anti-rigging guardrails: no relabelled constants, no
imported Born rule/singlet, and the `A/R+B·R²` well's stability explicitly
discounted as non-evidence.)*

## Why antipodal focusing matters

On a closed surface, a wavefront does not dissipate into infinity; it
reconverges at the antipode. In an embedded `S³`, that reconvergence
can have inner and outer bulk components — wavefronts arriving at the
antipode from opposite normal directions in the embedding. If those
focused wavefronts deform the embedding shell strongly enough, the
caustic may nucleate a non-orientable throat rather than reflect.

This picture promotes "particle creation" from a postulated quantum
event to a geometric instability. The threshold energy for nucleation
should correspond to the lowest stable Tangherlini eigenvalue; below
threshold, the antipodal focus disperses and the geometry relaxes; above
threshold, a self-consistent throat persists. The lepton-sector
threshold `2 m_e c²` is now **derived**
(`pair_production_threshold_probe`, PR #58) as twice the lowest stable
throat configuration, with the pair forced by Hopf-charge / antipodal-`Z₂`
conservation (one Hopf charge per throat ⟹ `Σ c₁ = 0` ⟹ C-conjugate
throat–antithroat pair), a bubble-nucleation barrier `R_c = 2σ/ρ` giving
the disperse-below / persist-above dichotomy, and the Schwinger critical
field `e E_S R_MID = m_e c²` tying the throat scale to the threshold.

**The exotic-matter question, and the narrowest gap (PR #167).** A throat
is a wormhole junction, and a thin-shell wormhole classically needs
NEC/WEC-violating *exotic* surface matter — the Israel/Lanczos surface
density `σ = −√f(a)/(2πa) < 0` — which the non-orientable (antipodal `Z₂` /
C-swap) gluing does **not** rescue. The honest resolution is *consistent
with* braneworld. The BAM throat metric `f(r) = 1 − (r_s/r)²` is
**Ricci-flat** (`R = 0`), and its effective 4D stress is **traceless** with
the `r⁻⁴` form `ρ_eff = −r_s²/(8πG r⁴) < 0`, `p_r = −ρ_eff`, `p_t = +ρ_eff`
— the **tidal-charge / bulk-Weyl form**: exactly what a projected bulk Weyl
tensor `E_μν` (traceless by construction) takes. Crucially it is *also*
what a real on-brane Maxwell field (Reissner–Nordström) takes, and **only
the sign distinguishes them** — a real brane gauge field gives `ρ > 0`,
whereas here `ρ_eff < 0`. So on-brane exotic matter is **avoidable if**
(i) the 5D embedding sources `E_μν` — by Shiromizu–Maeda–Sasaki a vacuum
brane obeys `G_μν = −E_μν`, forcing `R = 0`, which is **met** — and
(ii) BAM carries **no fundamental brane gauge field** that would force the
Reissner–Nordström reading. The **necessary conditions are met**; the
**sufficient** step — the explicit 5D embedding (the
Dadhich–Maartens–Papadopoulos–Rezania / Bronnikov–Kim tidal-charge
construction) — is **cited, not re-solved**, so this is a *consistent-with*,
not a proof. Nor does it evade the `f = 0` **horizon**: that locus is null
and degenerate, and the surgical surface term merely **vanishes** there,
relocating σ rather than removing it. This is the strongest "consistent-
with" the audits have reached — narrow, specific, and closable by a 5D
embedding calculation (`israel_junction_weyl_split_probe`, PR #167).

**The embedding, supplied — the gap closes (PR #168).** That 5D
calculation is now done, as an explicit **global regular** embedding (not
a Campbell–Magaard local-existence series). The BAM throat is the
**equatorial (`χ = π/2`) totally-geodesic slice** of the 5D
Schwarzschild–Tangherlini bulk `ds²₅ = −F dt² + dρ²/F + ρ²dΩ₃²`,
`F = 1 − μ/ρ²`, with `μ = r_s²`. The equator is a Z₂ fixed-point set, so
`K_μν = 0` (a tension-free, matter-free brane); the induced 4D metric is
exactly `f = 1 − (r_s/r)²`; the projected bulk Weyl equals the brane tidal
stress, `E_μν = −G⁴_μν` (verified to `~10⁻⁸`); the bulk is Ricci-flat (an
ordinary 5D vacuum); and the 5D Kretschmann `K₅ = 72 μ²/ρ⁸` is **finite
throughout** the exterior `ρ ≥ r_s` (max `72/r_s⁴` at the throat), the only
singularity `ρ = 0` lying behind the regular 5D Killing horizon `ρ = r_s`,
with the extra dimension `χ` compact and regular. Three checks and the
regularity gate pass: the bulk-Weyl reading is **realised**, not merely
consistent-with — no exotic brane matter, no brane gauge field — and the
`f = 0` throat is identified as the **regular 5D Killing horizon** (an
improvement on the #167 caveat: regular, not singular). The honest residue
is that the throat sits at that (regular) horizon, the brane is the
tension-free totally-geodesic slice (`μ = r_s²` fixing the bulk mass), and
it is the exterior embedding `ρ ≥ r_s`
(`global_regular_5d_embedding_probe`, PR #168).

**Why the throat is non-orientable but the bulk is not (PR #169).** The
non-orientability of the throat is not an extra assumption — it is forced
by a dimension-parity property of the antipodal (J) quotient. The antipodal
involution `J: x ↦ −x` on `Sⁿ` has orientation determinant `(−1)^{n+1}`:
orientation-**preserving** for odd `n` (orientable `RPⁿ`),
orientation-**reversing** for even `n` (non-orientable `RPⁿ`). The bulk's
angular sphere is `S³` (odd), so `S³ / J = RP³` is **orientable** (det +1);
the throat mouth is the brane's angular `S²` (even), so `S² / J = RP²` is
**non-orientable** (det −1). The *same* free isometric involution thus acts
oppositely on the two — they sit one dimension apart — and since `J` is
free (`−x = x ⟹ x = 0 ∉ Sⁿ`) and an isometry of the round angular metric,
the 5D Tangherlini geometry descends to the quotient cleanly. In the #168
coordinates `J = (χ,θ,φ) ↦ (π−χ, π−θ, φ+π)` fixes the equatorial `χ = π/2`
brane and restricts there to the `S²` antipodal map, so the #167
non-orientable throat is exactly the `RP²` cross-cap inside the orientable
`RP³` bulk. The split also lands where BAM's spinor structure lives:
`RP³ ≅ SO(3)` is orientable and spin, while `RP²` admits only a **Pin**
structure — the half-twist of the spin-½ / fermionic character, the same
orientability grading as the C-swap (`C = iσ_y`, `T² = −1`; PR #63) and the
even-`k` absence (PR #67) (`tangherlini_j_quotient_probe`, PR #169).

**The Pin⁻ mouth delivers Fermi statistics (PR #170).** The Pin structure
is not just a topological remark — it carries the physics. The throat mouth
`RP²` has Stiefel–Whitney classes `w₁ = a`, `w₂ = a²`, so it admits **no
Spin and no Pin⁺ structure, only Pin⁻** (`w₂ + w₁² = 0`): a unique, definite
spinor structure. That Pin⁻ spinor is spin-½ — a 2π rotation acts as
`R(2π) = exp(−iπσ_z) = −I`, with only `R(4π) = +I` — and by the
Finkelstein–Rubinstein construction the exchange of two identical throats is
homotopic to a 2π rotation of one (the two-particle configuration space in
≥3D has `π₁ = ℤ₂`, the exchange generator mapping to the 2π-rotation
generator). The exchange sign is therefore **−1**: the two-throat
wavefunction is antisymmetric, the spin-statistics connection realised by the
*same* holonomy that gives `2π = −1`. Antisymmetry forces Pauli exclusion
(occupation `n_p ∈ {0,1}`, against the Bose `{0,1,2,…}`), and filling the
Fermi sphere yields the degenerate **Fermi equation of state**: `P = ⅔u`,
`P ∝ n^{5/3}` (`Γ = 5/3`, non-relativistic) and `P = ⅓u`, `P ∝ n^{4/3}`
(`Γ = 4/3`, ultra-relativistic), with a strictly positive `T = 0` degeneracy
pressure — the support of white dwarfs and neutron stars — that a Bose gas
lacks. Computed here: the Pin⁻ classification, the spinor 2π sign, and the
Fermi-gas EoS; cited (not re-derived): the Finkelstein–Rubinstein
exchange↔rotation homotopy, the one configuration-space theorem linking the
throat's internal Pin holonomy to the physical exchange
(`pin_rp2_fermi_statistics_probe`, PR #170).

**The exchange, on the correct non-orientable footing (PR #171).** The
Finkelstein–Rubinstein homotopy cited above is the *orientable* result, and
the throat mouth is the non-orientable `RP²` (Pin⁻), so it does not transfer
for free. The correct framework is **geon statistics** (Friedman–Sorkin;
Aneziris–Balachandran–Bourdeau–Jo–Ramadas–Sorkin; Dowker–Sorkin), where a
geon's statistics is a representation of `π₁` of the configuration space and
the spin–statistics correlation is a theorem *with hypotheses, known to fail
for some geons*. Computing `π₁` of the two-mouth configuration space: the
exchange `σ` has `σ² = e` (in ≥3 spatial dimensions the symmetric group, no
braiding — only the ±1 statistics), the single geon's 2π rotation acts as
`−I` (spinorial; the Pin⁻ holonomy and Friedman–Sorkin's spin-½), and —
because the mouth is non-orientable — there is an orientation-reversing loop
`τ_i` the orientable argument never sees. That reversal carries a
**reflection**, and `RP²` admits **Pin⁻ only**, in which a reflection
**squares to −1** (Pin⁺, which `RP²` does not admit, would give `+1`) — the
ingredient that makes the non-orientable exchange sign well-defined and
fermionic. Non-orientability also makes the geon **achiral** (its own mirror
image), meeting the theorem's handedness hypothesis automatically. So the
−1 (Fermi) **survives** the Pin⁻ mouth, now on the right footing —
**conditional** on the Dowker–Sorkin exchangeability ("slide") hypothesis,
which holds for identical asymptotically-flat throats and is cited, not
derived from the full BAM field theory. The remaining honest gap is that
hypothesis (and the field-theory mapping class group), not the spinor sign
or the reflection algebra (`geon_statistics_pi1_probe`, PR #171).

**The equation of state, measured (PR #172).** Where #170 *assumed*
antisymmetry and read the index `5/3` off the analytic Fermi integral, and
#171 derived the `−1` exchange sign topologically, a companion simulation
(kept on the same branch for comparison) **measures** the equation of state
from a many-throat ensemble. `N` identical throats are free fermions in a
cubic box, the `−1` sign realised as Pauli single-occupancy of the box
modes, and the filled Fermi sea built by level-filling. From the volume
derivative `P = −dE/dV` the virial ratios `P/u = 2/3` (non-relativistic)
and `1/3` (ultra-relativistic) emerge; from the filled-mode energy sum,
finite-size-extrapolated by the Weyl correction `Γ(N) = Γ∞ − a·N^{−1/3}`,
the polytropic index is measured as `Γ = 1.6665 ≈ 5/3` and `1.3332 ≈ 4/3`
(0.01% from target) — outputs of the simulation, not a formula. A Bose
control (all `N` in the ground mode) gives `Γ = 1` with a vanishing `T = 0`
degeneracy pressure, so the stiffening is a measured consequence of the
exchange sign. The three routes — assumed-analytic (#170), topological-sign
(#171), and measured (#172) — agree
(`measured_fermi_eos_ensemble_probe`, PR #172).

The same throat↔antithroat nucleation channel later supplies the
neutrino's Majorana suppression (`seesaw_scale_nucleation_compliance_probe`,
PR #87): a `ΔL=2` Majorana mass *is* a throat↔antithroat flip, and the
single-state version of `Σ c₁ = 0` selects exactly the chargeless
`k = 0` (neutrino) sector — `0 → −0 = 0` is allowed, `±1 → ∓1` is not.
The seesaw scale is *not* the static barrier height (`E_c ≈ 2.8 keV`,
~10⁸ too small for the required ~TeV) but the **tunnelling amplitude
through** the barrier, `m_ν = m_D·e^{−S}`, so `M_R = m_D·e^{S}` with a
modest, generation-stable bounce action `S ≈ 15–18` — recasting PR #86's
open ~TeV scale as the instanton number this nucleation picture already
owed.

The refinement that distinguishes the present program from "particles
as static defects" is that the throats produced this way are not
required to remain at rest. A particle is a **moving topological
boundary condition** — two mouth positions `X₁(t), X₂(t)` on `S³`, a
bulk throat length `L_throat(t)`, and a time-dependent transport map
`T(t): T_{X₁}S³ → T_{X₂}S³` between tangent frames. The covariance
(`stable_moving_throat_probe`, PR #59) and Hopf-holonomy Wigner rotation
(`spin_wigner_rotation_probe`, PR #60) now verify that the boosted
throat is a genuine relativistic spin-½ particle. The focus is the
*trigger*; the particle is the persistent topological response.

**The focusing, computed (PR #166).** The antipodal reconvergence is no
longer only asserted. The zonal sector of `S³` reduces exactly to a 1D
wave on the string `[0,π]` (modes `sin((ℓ+1)χ)`), and the physical field
`ψ = f/sin χ` carries the geometric focusing factor `1/sin χ`. A
**conformal** wave packet (`ω_ℓ = (ℓ+1)/R`) launched near `χ₀` refocuses
**exactly** at the antipode `π−χ₀` at `t = πR` (half the great-circle
period) — the identity `ψ(χ,πR) = −ψ(π−χ,0)` holds to machine precision
and the amplitude fully recovers — then **revives** to its initial state
at `t = 2πR` (the sub-threshold focus passes through and re-disperses, the
geometry relaxing). The sharp focus **requires** conformal coupling: the
minimally-coupled tower `√(ℓ(ℓ+2))` dephases and blurs the caustic, so the
same conformal coupling that makes the `S³` vacuum tower equally spaced
(`berger_r_unification_audit_probe`, PR #165) is what makes the antipodal
caustic sharp. The caustic energy density `∝ 1/sin²χ` is regularized by
the spectral cutoff `ℓ_max ∼ R/R_MID`, so a delocalized, `S³`-wide wave
reconcentrates onto the throat scale — the dynamical bridge that lets a
diffuse wave reach the local nucleation density of the inherited `2 m_e c²`
threshold (`antipodal_focusing_threshold_probe`, PR #166). The probe maps
the *trigger* and applies the threshold; the nonlinear throat formation is
named, not simulated.

**The focused pulse / extended-wavefront bridge to the QCD shell**
(`throat_to_shell_transition_probe`, PR #68; `shell_to_qcd_match_probe`,
PR #69) extends the same antipodal-focusing story to the quark sector:
higher excitations of the focused lepton-throat pulse delocalize into a
QCD shell channel (extended-character wavefront), reproducing the
documented quark-sector structural invariants (`Z₂` partition,
`3 × 2 = 6` flavors, heavier scale, extended character). The lepton
throat and the QCD shell are two mode geometries of the same `S³`
closure skeleton.

**Quantitative QCD-shell development (PRs #76–#80).** The shell
channel is built out into a quantitative basis: 6-state `(l, n, p)`
shell waveguide (PR #77), `χ_n` derived from cavity-mouth boundary
stress (PR #79), BAM-native color algebra `SU(2) × Z₂` (PR #80). The
user's reframe makes the picture sharp: **"Quarks do not pass through
the throat; they are the wavefronts that resolve the cavity itself."**
What closes: the structural machinery. What remains open: the
inter-generation mass hierarchy (~9 orders in mass²) is outside the
scope of any BAM color algebra acting on the 6-state shell basis;
the phenomenological compensator `n_part = 233` (from the v3 lepton-
shaped fit) survives but with sharply identified scope (PR #76's
diagnosis). Most plausible derivation route: **Pati-Salam SU(4)**
unifying throat-leptons and shell-quarks via a quantitative throat↔
shell algebra map (beyond PR #68's structural transition; genuine
open work).

## What success looks like — falsification tests

The next phase of BAM was organized around demonstrations, not parameter
fitting. Each is a test the existing framework can be put to that admits
a clean pass-or-fail. Through PR #74, the program's most exposed tests
have largely **passed** — the entries below now record the outcome of
each falsification, not its prospectus. The remaining genuinely open
items are concentrated in the loop-measure sector and the quark `β`
lock.

**Odd-`k` classification (spin-statistics).** _Closed_ —
`docs/odd_k_closure_lemma.md` (the original closure lemma); upgraded to
a **classification** by `even_k_absence_probe` (PR #67). The closure
lemma states (i) even `k` and odd `k` both admit valid closure boundary
conditions on the throat — even `k` is orientation-preserving closure
on the doubled cover, odd `k` is orientation-reversing closure across
the non-orientable throat; (ii) under the locked baseline the Layer-1
ledger sum is identically zero mod 2π for every integer `k`. The
upgrade is that **`k mod 2` is the orientability / spin-statistics
grading**: each throat pass applies `T = iσ_y` (`T² = −I`, B2); the
spinor monodromy `T^k` is off-diagonal for odd `k` (opposite `Z₂`
class — orientation-reversing across the non-orientable throat = a
spin-½ fermion) and diagonal for even `k` (same class —
orientation-preserving on the orientable double cover `S³` = bosonic).
Charged leptons are spin-½ Dirac fermions (PRs #59–#66), hence the odd
class; even-`k` (bosonic) is excluded. So the original "choice of
sector" is upgraded to a **selection rule**: which `k` enter the
charged-lepton spectrum is forced by spin-statistics on `T² = −I`, not
chosen. The upper bound `k ≤ 5` is the three-generation boundary
(PR #70), with `k_5 = D_bulk = dim(S³)+2 = 5` (PR #73).

**Rigidity and uniqueness of the ladder (PR #174).** The selection rule
can be stressed: is the discreteness rigid against the continuous geometry,
and is it unique to this geometry? Running the #173 inverse problem on the
odd-`k` ladder as the discrete feature answers both. The continuous
deformation space splits (from the #173 Jacobian) into 10 *active*
directions that move the masses and CKM **linearly** (scaling exponent
≈ 1.0), 10 *null/compensator* directions that are **flat to first order**
(exponent ≈ 2.0, ~10⁴× smaller response), and *mixed* directions that are
active-dominated (≈ 1.0) — so nonlinear effects do **not** break the local
rank story; the null leakage stays quadratic. The odd-`k` labels and the
generation count are invariant under *every* active, null, and mixed
deformation, because they are integer winding plus the ℤ₂ orientability
grading (`T² = −I`) — discrete topological data that lives **outside** the
entire continuous deformation manifold (there is no generation-number
knob). The discreteness is therefore structurally **forced**, not an
emergent near-integer that could drift. And it is **unique**: an orientable
geometry (`T² = +I`) gives the orientation-preserving even/bosonic sector,
not an odd-only fermion ladder, while the specific `{1,3,5}` needs
`k ≤ k_5 = 5 = D_bulk` — so odd-`{1,3,5}` is the joint signature of the
non-orientable antipodal spin structure and the 5D bulk (an
exclusion/signature argument within BAM, not a no-go against every
conceivable alternative) (`odd_k_ladder_rigidity_probe`, PR #174).

**Can a continuous geometry evolve into the discrete sector? (PR #175).**
The static rigidity of #174 has a dynamical counterpart: a nonlinear
antipodal-focusing PDE sandbox (`nonlinear_antipodal_focusing_pde_probe`)
— a focusing nonlinear Schrödinger field `i∂_t ψ = −∂_χχ ψ − g|ψ|^p ψ` on
the antipodal ring, with the discrete sector taken as the winding number
`Q = (1/2π)∮ d(arg ψ)`. The answer is **yes, but only through the caustic**.
Smooth evolution **conserves** `Q` exactly while `|ψ| > 0` (the discrete
sector is locked out of continuous evolution — the dynamical confirmation of
#174). The only gate into it is an amplitude-zero **node**: because the
winding is a homotopy invariant of maps to `ℂ∖{0}`, interpolating between
`Q = 0` and `Q = 1` forces a zero of `|ψ|` located **exactly at the antipode**
— the focus — so the antipodal focusing of #166 is precisely what drives the
field toward the gate. Whether the nonlinear focusing reaches that core
depends on a **critical mass** (below it the field disperses and stays
continuous; above it concentrates toward the core — the disperse/persist
threshold of #58/#166, now simulated nonlinearly where #166 had deferred
it), and the winding jump at the core is **quantized ±1**: a discrete
response to a smooth focusing drive. The honest scope is a reduced 1D ring
model (`Q` proxies the discrete `k`, the collapse core proxies throat
nucleation, the critical-NLS collapse is marginal) — the conceptual answer
is robust, the numbers model-dependent
(`nonlinear_antipodal_focusing_pde_probe`, PR #175).

**Does real GR back the focusing threshold? (PR #176).** The #175 sandbox
used a 1D ring with an *ad-hoc* focusing nonlinearity `g|ψ|^p`; the next
step replaces the proxy with real (weak-field) general relativity. A
semi-dynamical, **axisymmetric** self-gravitating scalar
(`self_gravitating_axisymmetric_probe`) evolves `ψ(r,θ,t)` under
`i∂_t ψ = −½∇²ψ + Φ ψ` with the metric potential `∇²Φ = 4πG|ψ|²`
(`g_tt = −(1+2Φ)`, the weak-field Einstein–Klein–Gordon / Schrödinger–Newton
system) — the field in the `(r,ℓ)` Legendre basis, the radial Laplacian by a
Dirichlet sine transform, and `Φ(r,θ)` from the axisymmetric multipole
Poisson each step (split-step, mass-conserving to ~10⁻³). The disperse/
collapse **threshold survives** under actual gravitational back-reaction:
below a critical mass the packet disperses (the metric stays shallow), above
it the self-gravity concentrates it (the metric well deepens, runaway). The
decisive check that this is GRAVITY and not a tuned nonlinearity is the
**`1/G` scaling** of the critical mass — `G=0.5: >3.2`, `G=1: 2.29`,
`G=2: 1.10`, halving from `G=1` to `G=2` (ratio `0.48 ≈ 0.5`). So real
weak-field GR backs the antipodal-focusing threshold of #166/#175 and the
nucleation of #58. Honest scope: semi-dynamical weak-field GR — the field
evolves while the metric responds quasi-statically (not full numerical
relativity); the threshold and concentration (the throat-formation analog)
are confirmed, but the strong-field endpoint (a horizon / a resolved throat)
is for full NR, and self-gravity sphericalizes (the monopole dominates), so
the collapse is predominantly radial — the axisymmetric machinery is
exercised, not a directional jet claimed
(`self_gravitating_axisymmetric_probe`, PR #176).

**The self-gravity threshold, hardened into a benchmark (PR #177).** PR #176
was a promising proxy; `self_gravity_threshold_hardening_probe` turns it
into a trustworthy PDE benchmark with controls, scaling, and robustness.
**Controls:** with `G = 0` (gravity off) the packet never concentrates at
any mass, and with repulsive gravity (`G < 0`) it never collapses — the
threshold requires attractive gravity, not the packet or the grid.
**Energy anchor:** the total energy `E = T + W` (kinetic + gravitational
self-energy) defines the rigorous binding threshold `M_bind` (where
`E = 0`); below it the core mass drains (disperse, `E > 0`), above it the
core holds (bound, `E < 0`) — the dynamical disperse/bound transition
tracks the energy sign, an independent physics check on the integrator.
**Scaling:** the product `M_bind·G` is constant to `0.69%` across
`G ∈ {0.5, 1, 2}` (`1.134, 1.133, 1.127`) — the `1/G` law sharpened from
#176's coarse 0.48 to <1% — and the Schrödinger–Newton invariant
`M_bind·G·w ≈ const` holds across widths. **Robustness:** `M_bind`
converges to ~1–2% under radial-grid refinement and the split-step
integrator conserves mass to ~10⁻³. The weak-field self-gravity collapse
threshold is therefore gravitational, energy-validated, `1/G`-scaling to
<1%, and grid-converged — a benchmark, not a proxy. Standing scope
unchanged: weak-field / semi-dynamical, with the strong-field endpoint
(a horizon / a resolved throat) left for full numerical relativity
(`self_gravity_threshold_hardening_probe`, PR #177).

**The throat-order field q(t,r,θ) (PR #178).** The arc had established three
discrete facts about the throat in three separate languages — the odd-k
winding ladder (#174), the forced antipodal amplitude-zero node (#175), and
the self-gravitating focusing threshold (#176/#177). `throat_order_field_probe`
introduces a single field that unifies them: a complex Ginzburg–Landau order
parameter `q(t,r,θ) = |q| e^{iφ}` with the Mexican-hat potential
`V(q) = (λ/4)(|q|² − q₀²)²`, whose ordered vacuum `|q| = q₀` fills the
orientable bulk and whose **topological defects ARE the throats**. **Two
phases:** `q = 0` is an unstable maximum (`V″ = −1 < 0`, the disordered
symmetric phase) and `|q| = q₀` a stable degenerate minimum (`V″ = +2 > 0`,
the broken-symmetry vacuum); the free phase φ is the U(1) a defect winds.
**The throat is a vortex:** the radial GL profile `f(r)` solving
`f″ + f′/r − k²f/r² = λ f(f² − q₀²)` with `f(0) = 0, f(∞) = q₀` exists for
each winding — `|q| = 0` at the core, healing to `q₀` in the bulk, the core
widening with k (core size 0.8/1.4/2.0 for k = 1, 3, 5). **Winding = the
discrete k:** the charge `∮∇φ/2π` is the integer winding (`π₁(S¹) = ℤ`),
conserved while `|q| > 0`; the realized sector is odd-k — the #174
orientability grading. **Core = the antipodal node:** the field must vanish
where the phase winds, so the defect core `|q| = 0` is precisely the forced
amplitude-zero node of #175 — reaching the discrete sector from the
continuous (winding-0) sector requires passing through a zero. **Nucleation
= the threshold:** the disordered `q = 0` is unstable, so under the GL
gradient flow any perturbed region rolls off zero to `q₀` and a fixed-winding
defect nucleates; the trigger that drives a region off zero is the
self-gravitating focusing of #176/#177 (`M_c ∝ 1/G`). Scope: this is the
**effective** Ginzburg–Landau level — q is introduced as the coarse-grained
order field whose defects are the throats; the microscopic `V(q)` (λ, q₀
from the 5D bulk action) and the dynamical q–metric coupling are the
follow-ups. The throat's three discrete facts become one object: a vortex of
`q(t,r,θ)` (`throat_order_field_probe`, PR #178).

**Self-gravity-driven throat-order instability (PR #178).** The #178
throat-order field introduced `q` but coupled it to the geometry only by
hand. `self_gravity_driven_order_probe` closes the loop: does the
self-gravitating concentration of #176/#177 merely BIND the wave, or DRIVE
the order parameter? The matter density `ρ = |ψ|²` (from the #176/#177
solver, actually run) becomes the control field of a density-dependent
Landau potential `V(q; ρ) = ½(a₀ − gρ)|q|² + (λ/4)|q|⁴`, so the order
field's effective mass² `a(ρ) = a₀ − gρ` changes sign at a critical
concentration `ρ_c = a₀/g`: below it `q = 0` is the only minimum
(disordered, merely bound), above it `q = 0` destabilizes and the order
parameter rolls to `|q| = √((gρ − a₀)/λ)` (ordered). **Merely bound:** a
sub-threshold packet (M = 1) reaches only `ρ_peak ≈ 0.06 < ρ_c`, and the
order field relaxed under the Ginzburg–Landau gradient flow stays at zero —
bound, no geometric order. **Drives order:** above the mass threshold
(M = 3) the collapse drives `ρ_peak ≈ 0.90 > ρ_c` and the order field
NUCLEATES a localized symmetry-broken domain (`max|q| ≈ 0.68`) at the
density peak (the throat core of #178). **Gravitational:** with gravity off
(`G = 0`) the same mass never crosses `ρ_c` (`ρ_peak ≈ 0.18`) and no order
nucleates; restoring `G` it does — the ordering inherits the `M_c ∝ 1/G`
gravity of #176/#177. **Dynamical:** driving `q` by the time-dependent
`ρ_peak(t)` of the collapse, the order parameter switches on only after the
density crosses `ρ_c` — a moving order front following the gravitational
concentration. So weak-field concentration does NOT merely bind: above a
critical concentration, reached only by the gravitational collapse, it
drives the throat-order parameter off zero and nucleates geometric order.
Scope: ONE-WAY coupling (`ρ → q`; the self-consistent q–metric back-reaction
is the next step); the constants `a₀, g, λ` — and so `ρ_c` — are effective
(the existence of a gravitationally-crossed concentration threshold is the
result, not its microscopic value); the spatial nucleation carries the usual
GL droplet-size barrier; still weak-field / semi-dynamical
(`self_gravity_driven_order_probe`, PR #178).

**Two-way ψ–Φ–q evolution: the self-consistent throat-soliton (PR #179).**
PR #178 coupled the geometry to the order field ONE way (`ρ → q`; `q`
neither gravitated nor acted on the wave). `two_way_psi_phi_q_probe` closes
the loop into the full two-way system of three co-evolving fields — the
matter wave `ψ`, the gravitational potential `Φ`, and the throat-order field
`q` — all descending from one energy functional `E[ψ,q] = ∫[½|∇ψ|² +
½κ|∇q|² + ½a₀q² + ¼λq⁴ − ½g|ψ|²q²] + W_grav[|ψ|² + μq²]`. Its fixed-mass
gradient flow is `∂_τψ = ½∇²ψ − Φψ + ½g q²ψ`, `∂_τq = κ∇²q − (a₀−g|ψ|²)q −
λq³`, `∇²Φ = 4πG(|ψ|² + μq²)`, so the four back-reaction channels are all
live: ψ↔Φ (Schrödinger–Newton, #176/#177), ψ→q (the density orders q,
#178), q→ψ (the ordered throat core binds the wave, NEW), q→Φ (the order
field gravitates, NEW); the ordering and binding terms share the same `g`,
so the coupling is consistent, not hand-wired. **Self-consistent:** the
coupled flow converges — energy monotone → plateau, q stationarity residual
`~10⁻⁴` — a self-consistent throat-soliton exists. **Two-way back-reaction:**
at super-threshold mass the order field nucleates and, versus the pure
Schrödinger–Newton soliton, the self-consistent state has a deeper well
(`Φ(0) = −3.18` vs `−3.03`, ~5% deeper) and a denser core (~13% denser) —
the throat traps the wave, which concentrates it, which strengthens the
order. **Saturation vs collapse:** with sub-critical self-gravity the
quartic `λq⁴` saturates the feedback into a stable bound soliton (`|q|`
plateaus; intermediate μ gives a denser soliton), but super-critical
self-gravity has no weak-field fixed point and the flow diverges (`max|q| →
31`, `Φ(0) → −252`) — the onset of strong-field gravitational collapse.
**Continuity:** below the ordering threshold the order field vanishes and the
system reduces exactly to the Schrödinger–Newton soliton of #176/#177 — the
#176 → #178 → #179 arc is one continuous system, switched by the matter
concentration. So the throat-order field is not a passive readout of the
geometry; closing the loop it back-reacts both ways, forming a
self-consistent throat-soliton. Scope: weak-field, semi-dynamical,
spherically reduced (the self-gravity sphericalizes, #176); the constants
are effective (the structure is the result); the stable soliton is
sub-critical and the strong-field runaway endpoint is for full numerical
relativity (`two_way_psi_phi_q_probe`, PR #179).

**ψ–Φ–q soliton hardening: stationarity, branch scan, basin map (PR #180).**
`psi_phi_q_soliton_hardening_probe` hardens the #179 two-way throat-soliton
(as #177 hardened #176) and re-examines its collapse claim with a
better-conditioned solver. **Stationarity:** putting `ψ`'s kinetic on a
spectral basis (`u = rψ`, DST; the order field `q` keeps its finite-difference
Laplacian — this is not a fully spectral ψ–q solver) so the relaxation and the
real-time step share the same `ψ` Laplacian, the relaxed state is a genuine
eigenstate (`‖Hψ − μψ‖/‖ψ‖ ≈ 10⁻⁴`, chemical potential `μ ≈ −1.45`), and
evolving `ψ` alone in the frozen self-consistent `(Φ, q)` background by a
unitary real-time split-step leaves it stationary (profile drift `~4×10⁻⁵`,
mass conserved to machine precision) — `ψ` is a stationary eigenstate of its
self-consistent potential, a real bound soliton (the fully coupled real-time
ψ–Φ–q dynamics is a follow-up). **Branch scan:** the soliton is a smooth
monotone family in mass (the order field switches on where `ρ_peak` crosses
`ρ_c`, near M ≈ 2.7) and in `q`'s self-gravity `μ` (`max|q|` 0.42 → 2.62,
`Φ(0)` −3.09 → −24.6 across the tested range `μ ∈ [0.05, 2]`, residuals
`≤ 10⁻³`, everywhere convergent). **A correction to #179:** #179 reported a
runaway collapse at super-critical `μ` (`|q| → 31`, `Φ(0) → −252`), but that
used a finite-difference (`np.gradient`) Laplacian; the spectral `ψ` kinetic
finds no collapse up to `μ = 2`, so the runaway was a discretization artifact
— the genuine large-`μ` limit is the soliton deepening out of weak-field
validity (`Φ(0)` −3.09 → −24.6 across the tested μ; the strong-field domain
for full NR), not a numerical runaway. **Basin:** the soliton is a robust
attractor — the full initial-condition grid (widths `w ∈ {1.2, 1.8, 2.6}`
crossed with seeds `∈ {10⁻², 10⁻¹}`, all six) flows to the same state
(`max|q|` spread ~1%, `Φ(0)` spread ~0.1%; a tiny seed `10⁻³` reaches the
same attractor more slowly). **Robustness:** the well depth
`Φ(0)` grid-converges (`N = 160 → 240 → 320`: −3.34 → −3.09 → −2.98) while
the pointwise core `max|q|` is more grid-sensitive (~10% per refinement, the
sharp core) — an honest caveat. What survives #179 — the soliton's
existence, two-way back-reaction, and threshold continuity — is confirmed and
hardened; the specific "runaway" claim does not survive as stated. Scope
unchanged: weak-field, semi-dynamical, spherically reduced, effective
constants; the deep-`μ` branch and the strong-field endpoint are for full
numerical relativity (`psi_phi_q_soliton_hardening_probe`, PR #180).

**Discrete invariant survival on the ψ–Φ–q soliton (PR #181).** The arc has a
continuous object (the #179/#180 self-consistent throat-soliton) and a
discrete one (the #174/#178 winding ladder). `discrete_invariant_survival_probe`
shows the continuous geometry CARRIES the discrete charge: dress the #180
soliton's ordered core (an equatorial loop of radius `R = 0.75`, where
`ρ = |ψ|² = 0.36 > ρ_c`, so `|q| > 0`) with a winding-k phase
`q = |q| e^{ikφ}`, and the topological charge `Q = (1/2π)∮∇φ = k` (exact to
~10⁻¹⁵). A winding-k vortex is sustained when the well beats the centrifugal
cost, `A² = (gρ − a₀) − (κ/R²)k² > 0`: the soliton sustains k = 1, 3; k = 5
exceeds it. **Survival:** under continuous norm-conserving (wave) evolution Q
is conserved to MACHINE PRECISION (`ΔW ~ 10⁻¹⁶`) for all k ∈ {1, 3, 5} with
`min|q| > 0`; under the order field's own dissipative gradient flow the
sustained windings k = 1, 3 survive (a perturbed vortex relaxes back, Q
conserved to ~10⁻¹⁵). **The criterion:** Q changes ONLY through `|q| = 0` —
the unsustained k = 5 is driven to a zero (`min|q| → 10⁻⁴`) and the charge
slips (5 → 2); survival ⟺ `|q| > 0`, exactly. That slip is the phase-slip /
topology-change event of PR #182. **Rigidity:** under 40 random
`|q| > 0`-preserving homotopies per sector the charge is unchanged in every
case (40/40 for k = 1, 3, 5) — a superselection charge outside the continuous
moduli (the #173/#174 rigidity, now on the dynamical soliton). So the
#174/#178 winding ladder rides the #179/#180 soliton untouched, except at the
amplitude zeros where topology changes. Scope: homotopy-invariance is exact;
the geometry is the reduced vortex-on-soliton (amplitude from the radial #180
soliton, winding azimuthal — the full 2D/3D self-consistent vortex-line
soliton is a follow-up); which rungs survive is set by the soliton's capacity;
the realized PHYSICAL ladder is odd-k {1, 3, 5} by the #174 orientability
grading, with its survival under a deformed bulk geometry the subject of #183;
weak-field (`discrete_invariant_survival_probe`, PR #181).

**The phase-slip / topology-change event (PR #182).** PR #181 showed the
winding charge survives the continuous evolution while `|q| > 0` and can
change only where `|q| = 0`. `phase_slip_topology_change_probe` dissects that
event — exactly how the invariant changes when `q` hits zero. **The
obstruction:** `Q` is a homotopy invariant of `q: S¹ → ℂ∖{0}`, so to change it
the field must leave `ℂ∖{0}`; the straight homotopy `(1−s)·[winding 1] +
s·[winding 0]` is FORCED through an exact zero (`min|q| = 2.5×10⁻¹⁷` at
`s* = 0.5`, located at `φ* = π`), where `Q` jumps 1 → 0 — there is no
nowhere-zero path between the sectors (the dynamical content of the #175
gate). **The quantum:** across that simple zero `ΔQ = −1` exactly — the
integrated winding density `∮∇φ` changes by `−2π` (one full turn removed at
the zero point); a generic simple zero carries unit topological charge, so
each elementary slip is `±1`. **The dynamics:** in a genuine ψ–Φ–q evolution
`Q(t)` is piecewise-constant and steps by `±1` EXACTLY at the instants
`min|q|(t) → 0` — a single slip (the unsustained k = 5 holds flat then steps
−1 to 4 at a zero) or a quantized STAIRCASE (k = 8 cascading
`[8,7,5,4,3,2]`, every step at an amplitude-zero event, shedding winding one
quantum at a time; a recorded `−2` step is two elementary slips unresolved in
sampling time). **Localization:** away from the slips (`min|q| > 0.1`) the
unrounded winding equals an integer to `10⁻¹⁵` — `Q` is a rigid integer
between events, ambiguous only at the measure-zero set of amplitude zeros. The
phase slip is the throat changing its winding / generation sector
(`k → k∓1`) through the amplitude-zero node — the #175 antipodal node, the
#178 defect core: the #175 gate made into the sector-CHANGING event itself.
With #181 (survival between events), the throat's winding is a conserved
topological charge that transitions ONLY at nodes; the realized ladder is
odd-k by the #174 orientability grading, its survival under a deformed bulk
geometry the subject of #183. Scope: the obstruction and the `±1` quantum are
exact (topological); the dynamical staircase is on the reduced
vortex-on-soliton loop (the full 2D/3D vortex-line reconnection is a
follow-up); weak-field (`phase_slip_topology_change_probe`, PR #182).

**Odd-k / generation-sector survival under a deformed bulk geometry
(PR #183).** PR #174 derived the odd-k charged-lepton ladder {1, 3, 5} (= 3
generations) from the non-orientable bulk (the throat closure `T = iσ_y`,
`T² = −I`). `odd_k_generation_survival_probe` shows that derivation is
topologically PROTECTED — it survives any smooth deformation of the bulk
geometry — closing the #181/#182 structure one level up, at the bulk. **The
grading is metric-independent:** the antipodal deck map is `−I` in any linear
frame, so `det = (−1)^dim` — the brane angular slice `S²/antipodal = RP²` is
non-orientable (`det = −1`), the bulk `S³/antipodal = RP³` orientable
(`det = +1`); the closure `T = iσ_y` has `T² = −I`, `½ tr T² = −1` (the Pin⁻
structure forced by `w₁² = w₂`); the grading `tr(T^k) = 2cos(kπ/2) = 0` for
odd k (off-diagonal, fermionic) and `±2` for even k (diagonal, bosonic).
**Survival:** a smooth metric/frame deformation acts on the holonomy by
orientation-preserving conjugation and on the deck map by a GL⁺ frame
change, neither of which can flip a determinant sign or a trace; across 1000
random deformations `½ tr T²` stays `−1` and the deck dets stay `∓1` to
machine precision (`~10⁻¹⁵`), with named squash/tidal deformations likewise.
**The generation count:** odd k ≤ `k₅ = D_bulk = 5 ⟹ {1, 3, 5} = (k₅+1)/2 =
3` (matching `LEPTON_BASELINE_DEPTHS`); `D_bulk` and the odd-parity selection
are topological, so the count survives every smooth deformation — not an
artifact of the round metric. **Changes only at a topology change:** the only
sector-flipping path, the non-metric `T(θ) = exp(iθσ_y)` driving `T² : −I →
+I`, has its orientability invariant cross zero at `θ = π/4` — a degenerate
spin structure, the topology-change event; smooth deformations act by
conjugation and never move `θ`, so they can never reach it (the exact
bulk-level analog of the #182 amplitude zero). **Unity:** the generation
sector is to the bulk geometry what the order-field winding is to the soliton
(#181/#182) — a topological charge robust to smooth deformation, changing only
at a singular / topology-change event. So the #174 round-metric derivation is
not special; the odd-k, three-generation structure is topologically protected
against any smooth deformation of the bulk. Scope: the invariance is exact
(topological: the deck determinant and the spin-closure / Stiefel–Whitney
class are metric-independent); the deformations are within the
orientability/spin-preserving class; this establishes robustness, not a
re-derivation; the result is purely topological — weak-field is not invoked
(`odd_k_generation_survival_probe`, PR #183).

**α as a protected boundary invariant, not a continuous tuning parameter
(PR #184).** PR #105/#143 classified the EM coupling α: the geometry derives
its STRUCTURE — the charge quantum `|c₁| = 1` (the integer Hopf number), the
`1/2π` closure loop measure (the `2π` in the Schwinger anomaly `a = α/2π`),
and the running — but not the VALUE `α⁻¹ ≈ 137` (the residual "137 problem").
`alpha_protected_boundary_invariant_probe` applies the #181/#182/#183
protected-invariant test to that derived structure: it is a PROTECTED BOUNDARY
INVARIANT, not a tunable continuum. **The charge quantum is a boundary
invariant:** the first Chern number of the BAM Hopf / spin-½ monopole over the
boundary S² (the Gauss-law charge `(1/2π)∮F`), by the exactly-quantized
Fukui–Hatsugai–Suzuki method, `c₁ = +1` (`|c₁| = 1`) — an exact integer.
**Protected:** across 30 smooth diffeomorphisms of the boundary geometry `c₁`
stays the same integer to `5×10⁻⁷` — it does not drift. **Not a tuning
parameter:** under the same deformations a generic continuous coupling
functional (the mean monopole potential `⟨A_φ⟩`) drifts `6.8%` on average
(up to `15.8%`) while `c₁` moves `5×10⁻⁷` — the discriminator (quantized +
deformation-invariant = protected; continuous + drifts = tuning) puts α's
charge quantum on the protected side. **The loop measure + topology change:**
the boundary flux `∮F = 2π·c₁` is quantized in units of the closure quantum
`2π` (fixing the `2π` of `a = α/2π`); and the charge integer changes ONLY when
the Berry gap closes — sweeping the gap parameter `m`, `C(m) = 1` while the gap
is open and jumps to `0` exactly at `m = 1`, where `min|d| → 0` (the
degeneracy crosses the boundary), the EM-boundary analog of `|q| = 0` (#182)
and `½ tr T² = 0` (#183). **Unity:** the EM charge quantum is to the boundary
what the order-field winding is to the soliton (#181/#182) and the generation
sector is to the bulk (#183) — a protected topological charge robust to smooth
deformation, changing only at a topology-change event. This does NOT derive
the value `α⁻¹ ≈ 137` — that residual stands (the 137 problem is unchanged) —
but it refines the #105/#143 ledger: the structure around α is specifically
PROTECTED, so α should be tested as protected-boundary-structure × one residual
scale, not fit as a continuous tuning family
(`alpha_protected_boundary_invariant_probe`, PR #184).

**Multi-throat mechanics & the exchange kernel from the GR soliton
(PR #185).** The arc built a single self-gravitating ψ–Φ–q throat-soliton
(#176–#180); `multi_throat_exchange_kernel_probe` takes TWO of them and
derives the EXCHANGE KERNEL from GR — no postulated statistics. It factorizes
as `K_exchange(R) = (−1)·K(R)`: a GR-geometric SPATIAL overlap times a
TOPOLOGICAL sign. **The spatial kernel** `K(R)` is the overlap of two actual
#180 throat-solitons separated by `R`, decaying smoothly from `K̂(0) = 1` over
the soliton size (RMS ≈ 1.27) — a GR exchange RANGE, not a postulated form
factor (`K̂`: `1.00, 0.79, 0.41, 0.15, 0.045, 0.003` at `R = 0,1,2,3,4,6`).
**The sign** is `−1` (fermionic), derived from GR: the large diffeomorphism
that swaps two throats is homotopic to a 2π rotation of one throat (the
Friedman–Sorkin / Dowker–Sorkin spin-statistics theorem for geons), and a 2π
rotation on the non-orientable Pin⁻ throat is `T² = −I` (`½ tr T² = −1`), so
the exchange phase is `−1`; a boson would need the orientable `T² = +I`
closure the throat does not have (#170/#174/#183). So the geometry SELECTS the
antisymmetric (Fermi) eigenvalue of the exchange operator `P` (`P² = 1`,
eigenvalues `±1`). **Pauli exclusion:** the antisymmetric two-throat state
`Ψ₋(z₁,z₂) = φ_a(z₁)φ_b(z₂) − φ_a(z₂)φ_b(z₁)` vanishes identically at
coincidence (`max|Ψ₋(z,z)| = 0` to machine precision — the determinant of two
equal rows), so two identical throats cannot occupy the same state; the boson
`Ψ₊` does not (it bunches). **The exchange hole + Fermi pressure:** the
exchange term `∝ K(R)²` carves an exchange hole of GR range = the soliton
size; macroscopically the exclusion fills a degenerate Fermi tower — with the
3D DOS `g(E) ∝ √E`, `N ∝ E_F^{3/2}`, `E ∝ E_F^{5/2}` ⟹ `E ∝ N^{5/3}` ⟹
`P = (2/3)(E/V) ∝ n^{5/3}`, polytropic `Γ = 5/3` — exactly the Fermi EoS
measured in #172. The GR-derived exchange kernel is the microscopic origin of
the Fermi pressure of throat matter. Scope: the exchange sign is exact /
topological (the Pin⁻ geon statistics, a GR large-diffeomorphism / mapping-
class-group representation); the spatial kernel is the rigid #180
soliton-overlap model (the single-particle orbitals) — the full two-body GR
problem (the two-throat metric, the gravitational direct/Hartree term, the
dynamical swap with back-reaction) is a follow-up; the Fermi index 5/3 is the
standard degenerate-gas result, here attributed to the GR-derived exchange
kernel; weak-field / semi-dynamical soliton
(`multi_throat_exchange_kernel_probe`, PR #185).

**Rigid soliton exchange-kernel hardening (PR #186).**
`rigid_soliton_exchange_kernel_hardening_probe` hardens the #185
rigid-soliton exchange kernel (as #177 hardened #176). **Normalization:** the
single-throat orbital is normalized (`∫|φ|² d³r = 1.000000`); the self-overlap
`K(0) = 1.001` reproduces the norm to 0.1% (the overlap-quadrature residual);
the kernel is parity-symmetric (`K(2) = K(−2) = 0.40963`, φ radial) and obeys
the Cauchy–Schwarz bound `K(R) ≤ K(0) = 1`. **Convergence:** refining the
overlap quadrature, `K(2)` converges to `< 0.01%` (the overlap integral is
well-resolved); the dominant uncertainty is the soliton PROFILE — rebuilding
the #180 soliton at `N = 240 → 320` shifts `K(2)` by `~2.9%`, the documented
#180 core grid-sensitivity, honestly identified as inherited (not a flaw in
the kernel). **Direct-term controls:** the DIRECT density-overlap
`D(R) = ∫ ρ_a ρ_b d³x` (the sign-independent Hartree channel) and the
EXCHANGE amplitude-overlap `K(R)` (the ±-carrying channel) are distinct
GR-geometric kernels, both decaying to zero at large R (the direct faster);
at far separation both vanish (distinguishable throats); and the direct is a
positive density overlap with no sign — identical for the boson (+) and
fermion (−) sectors — so the Pin⁻ `−1` lives purely in the exchange channel,
the direct being the control that isolates it. Consequently the two-body
energy splits as `E = E_direct ∓ E_exchange`, the structure the #187
Hartree–Fock sandbox evaluates against an interaction. Scope: the rigid
soliton-overlap kernel is now a trustworthy benchmark; convolving the
overlaps with `V` to get the Hartree and exchange energies is PR #187;
weak-field (`rigid_soliton_exchange_kernel_hardening_probe`, PR #186).

**Two-throat Hartree–Fock sandbox: direct plus exchange terms (PR #187).**
`two_throat_hartree_fock_probe` convolves the #186 hardened overlap kernels
with an interaction `V` (a screened-photon Yukawa stand-in for the BAM
throat-fibre exchange) to build the two-throat Hartree–Fock energy. Because
two displaced throats are NON-orthogonal (orbital overlap
`S(R) = ⟨φ_a|φ_b⟩ ≠ 0`), the physical, OVERLAP-NORMALIZED energy is
`E±(R) = (J(R) ± K_ex(R))/(1 ± S²)`, with the DIRECT (Hartree) numerator
`J = ∫∫ ρ_a(r₁) ρ_b(r₂) V(r₁−r₂)`, the EXCHANGE numerator
`K_ex = ∫∫ τ(r₁) τ(r₂) V(r₁−r₂)` (`τ = φ_a φ_b`), and `(1 ± S²)` the overlap
normalization — `J, K_ex` are the unnormalized HF numerators, computed from
the actual #180 throat-soliton orbitals by a 3D-FFT Coulomb solve.
**The integrals** (positive, decaying): `S = 1, 0.79, 0.41, …`,
`J = 0.039, 0.031, 0.017, …`, `K_ex = 0.039, 0.024, 0.006, …` at `R = 0,1,2`
— the direct dominating (the exchange has the shorter, overlap-set range).
**The normalized energies & ordering (repulsive V):** for the repulsive
screened interaction the antisymmetric (fermion) branch
`E₋ = (J − K_ex)/(1 − S²)` sits BELOW the symmetric (boson)
`E₊ = (J + K_ex)/(1 + S²)` at every finite separation (the exchange hole
lowers the GR-selected antisymmetric Pin⁻ state; the gap closes as `S → 0`) —
an ordering SCOPED to a repulsive `V` (an attractive interaction reverses it).
**Pauli at coincidence — the zero vector:** as `R → 0`, `S → 1` and BOTH the
numerator (`J − K_ex`) → 0 AND the normalization (`1 − S²`) → 0, so the
antisymmetric combination `Ψ₋ = (φ_a φ_b − φ_b φ_a)/√(2(1 − S²))` is the
ZERO VECTOR — two identical fermions cannot occupy the same orbital, the
antisymmetric state is Pauli-FORBIDDEN (not a state with zero interaction
energy); the boson `E₊` survives (bunching). For a contact `V` the numerator
`J − K_ex = 0` at all `R` (`J = K_ex = g·D(R)`, the hardened #186 direct
overlap), so `E₋ = 0` at every finite separation, the state forbidden only at
exact coincidence. The overlap and numerators vanish at far separation
(orthogonal, distinguishable), and the energies are grid-convergent to
`< 0.1%`. So the multi-throat mechanics close: the GR-derived exchange kernel,
dressed by an interaction and properly overlap-normalized, gives a
Hartree–Fock energy whose antisymmetric (Pin⁻ fermion) branch sits below the
symmetric (boson) branch for a repulsive interaction and is the
Pauli-forbidden zero vector at coincidence — the exchange interaction and the
Pauli physics of throat matter, from GR. Scope: a sandbox — rigid #180
orbitals (the self-consistent two-throat solve is a follow-up), a
screened-photon (Yukawa) regulated stand-in for the BAM Coulomb/photon
interaction, spatial exchange only (the spin factor is the Pin⁻ `−1`),
energies in code units; the structure is robust, the precise numbers carry
the #186 soliton-profile ~3% uncertainty; weak-field
(`two_throat_hartree_fock_probe`, PR #187).

**Adiabatic two-throat exchange holonomy (PR #188).**
`adiabatic_exchange_holonomy_probe` makes the #185 Pin⁻ exchange sign
operational — measured, not asserted. The relative-coordinate configuration
space of two identical throats is `(ℝ³∖0)/ℤ₂ ≃ RP² × ℝ₊`, whose angular factor
`S²/antipodal = RP²` is the BAM antipodal closure itself (#169/#170); the
exchange `r → −r` is the generator of `π₁(RP²) = ℤ₂`, and by the
Finkelstein–Rubinstein / Friedman–Sorkin spin-statistics theorem it is
homotopic to a 2π rotation of one throat. **Measured:** path-ordering the
spin connection along the swap (2π) loop — `dU/ds = −i(ω·σ/2)U` — gives the
adiabatic holonomy `Hol = −I` to machine precision (`‖Hol + I‖ ~ 10⁻⁶`): the
throat's spin-½ state returns to MINUS itself, the measured exchange sign
`⟨ψ|Hol|ψ⟩ = −1`, the Berry phase `π`. **Topological:** a wandering-axis swap
gives the same `−I` (the ℤ₂ homotopy class), converging as the transport is
refined — any way of doing the exchange gives the same `−1`. **Controls:** a
double-swap (4π, two exchanges) gives `+I` (two fermion exchanges = a boson)
and a contractible loop gives `+I`, so the `−1` is the single-swap (odd)
class. **Pin⁻:** the `−1` is the monodromy `T = iσ_y`, `T² = −I`
(`½ tr T² = −1`; #170/#174/#183) — the throat is a spin-½ spinor via the
non-orientable RP² closure, so its 2π/swap holonomy is `(−1)^{2j} = −1`
(a scalar throat would give `+1`). The adiabatic holonomy IS this `T² = −I`,
now transported along an explicit swap path rather than read off the algebra.
Scope: operationalizes the FR/geon-statistics result (the holonomy is exact
and topological — the ℤ₂ class — so the `−1` is exact); the swap path is the
reduced relative-coordinate / frame model, the spin-statistics connection
(exchange ≃ 2π rotation) the FR theorem cited, the throat's Pin⁻ spinor nature
the #170 result; the adiabatic limit assumed; it complements the #185–#187
spatial exchange kernel and Hartree–Fock energies — the statistics/holonomy
layer of the multi-throat mechanics
(`adiabatic_exchange_holonomy_probe`, PR #188).

**Self-consistent two-throat Hartree–Fock relaxation (PR #189).** PR #187
assembled the two-throat HF energy with RIGID orbitals and flagged the
self-consistent solve as a follow-up; `self_consistent_two_throat_hf_probe`
does it. A genuine self-consistent-field iteration relaxes the two same-spin
throats (the Pin⁻ antisymmetric sector of #185/#188) in each other's DIRECT
(Hartree) + EXCHANGE (Fock) field. Each orbital is relaxed in an
orbital-specific, self-interaction-free Fock operator `F_i = h + J_{≠i} −
K_{≠i}` — kinetic + confinement `h`, the direct (Hartree) field
`J_{≠i}(x) = ∫ |φ_{≠i}(x')|² V(x−x') dx'` of the OTHER throat, and the
non-local exchange with the other orbital — which is the EXACT variational
derivative of the reported energy `E = Σ⟨i|h|i⟩ + (J₀₁ − K₀₁)` (the same
functional for the operator and the energy), by imaginary-time gradient
descent with the two orbitals kept orthonormal. **Convergence + robustness:**
the energy descends MONOTONICALLY (`−3.808 → −3.905`, final `ΔE ~ 10⁻⁸ → 0`)
to a self-consistent fixed point — the orbitals come to rest in the field they
produce — and the fixed point is ROBUST across seeded restarts (five random
localized initial orbitals converge to `E ≈ −3.905`, spread `~8×10⁻³`); so it
is a self-consistent VARIATIONAL FIXED POINT (robustly reached, not certified
the global ground state). The imaginary-time relaxation avoids the
eigenstate-swapping oscillation a naive diagonalization-SCF shows for these
near-degenerate orbitals. **Relaxation:** the self-consistent energy lies
`2.54%` below the rigid 1D #187-STYLE reference (the unrelaxed orbitals with
the full HF energy — the 1D analogue of #187's rigid-orbital evaluation, not
the 3D number) — the variational gain from optimizing the orbital shapes that
#187 omitted. **Deformation:** the two-throat density polarizes in the mean
field (RMS width `2.65 → 3.06`, overlap with the rigid density `0.978 < 1`).
**Exchange:** with the CONSISTENT self-interaction-free control (the exchange
dropped from BOTH the Fock operator and the energy, so its operator and energy
are the same functional), turning off the non-local Fock exchange raises the
energy by `0.567` (`E_HF = −3.905` vs `E_Hartree = −3.338`) — for the two
same-spin throats the exchange substantially lowers the energy (the exchange
hole of #186/#187 keeps the like throats apart, reducing the repulsive direct
energy), so the `−1` of #185/#188 does real work in the self-consistent mean
field. So the two throats are no longer rigid: they settle into a
self-consistent variational fixed point of their mutual direct + exchange
field. Scope: a 1D sandbox SCF (external double-well confinement standing in
for the throats' self-gravity, a screened-photon interaction stand-in,
spatial-orbital same-spin HF, in 1D for tractability); the SCF itself is
genuine (self-interaction-free operator consistent with the energy) and the
qualitative physics robust; the relaxed state is a variational fixed point
(not certified the global ground state); the full 3D self-gravitating
two-throat SCF (relaxing actual #180 solitons) is the follow-up; weak-field,
code units (`self_consistent_two_throat_hf_probe`, PR #189).

**The BAM Coulomb-photon kernel for the two-throat HF (PR #190).** PR #187 and
#189 used a screened-photon (Yukawa) interaction as a regulated stand-in;
`bam_coulomb_two_throat_hf_probe` replaces it with the genuine BAM
Coulomb-photon kernel — the UNSCREENED Coulomb `V(d) = 1/(4πd) ⟷ 1/q²`, the
photon propagator BAM derives from the throat-fibre exchange geometry
(#42–#44), the flat-space limit of the S³ scalar Green function
`G(ψ) = ((π−ψ)cotψ − ½)/(4π²R)` (the repo's `s3_green_potential`). **The
kernel** is verified: the isolated-system Coulomb of a unit point source
reproduces `1/(4πd)` to machine precision, and near the source `G·4πs = 0.957
→ 1`, the Coulomb coefficient (`s = Rψ` the geodesic distance), so on the
local weak-field patch the throats see the unscreened photon with S³ curvature
corrections `O(1/R²)` carried by `G`. **The regulator** is the Hockney
zero-padded open-boundary Coulomb (density padded to a 2× box, convolved with
the free-space `1/(4πr)` Green function) — a proper isolated-system scheme,
NOT physical screening; validated against the analytic Gaussian Coulomb
self-energy (ratio `0.9992`, ~0.08%). **The recomputed energies** (on the #180
orbitals): the direct `J(R)` is now correctly LONG-RANGED — `J(6) = 0.0133 ≈
1/(4π·6)` (ratio `1.003`), the point-charge Coulomb tail — unlike the Yukawa
stand-in's exponential decay; the exchange `K_ex(R)` stays short-ranged
(overlap-set), so far-apart throats feel the Coulomb direct field but not the
exchange. **The #187 physics survives:** with the overlap-normalized
`E±(R) = (J ± K_ex)/(1 ± S²)`, for the repulsive photon the antisymmetric
(Pin⁻) branch `E₋` lies below the symmetric `E₊` at every finite separation
(fermion-lower), and at coincidence `J = K_ex` and `S → 1`, so the
antisymmetric state is the zero vector (Pauli-forbidden). The statistics are a
property of the geometry (the Pin⁻ sign + the overlap structure), not of the
interaction's screening — the Yukawa was a faithful short-range stand-in, now
upgraded to the real long-ranged photon. Scope: the flat Coulomb limit (the
`O(1/R²)` S³ curvature corrections are carried by `G` but not applied — the
weak-field local patch); the Hockney is a numerical open-boundary regulator,
not screening; the orbitals are the rigid #180 throat-solitons (the
self-consistent #189 SCF with the Coulomb kernel is the follow-up);
weak-field, code units (`bam_coulomb_two_throat_hf_probe`, PR #190).

**The transactional Compton propagator (PR #213).** The Compton arc
(#35/#45/#46, #211) derived the propagator's *spatial* structure from
geometry but **imported its time structure** — the Feynman iε
prescription, the split into two time-ordering completions, and their
coherent relative phase. `compton_transactional_propagator_probe`
derives all three from the **frozen classical bulk** via complete
(Wheeler–Feynman/transactional) histories
(`docs/compton_transactional_propagator.md`). **The construction**: (1)
the bulk is *static*, so t → −t is an isometry and the elementary field
is the time-symmetric Ḡ = (G_ret + G_adv)/2 — the geometry cannot
prefer retardation (tower identity G_ret(t) = G_adv(−t) at machine
zero on ω_ℓ = (ℓ+1)/R); (2) the bulk is *closed*, so Wheeler–Feynman's
perfect-absorber hypothesis is a **theorem of the geometry** — every
retarded front refocuses at the antipode at t = πR (measured: 50% of
the |G|² mass within 0.3 rad, vs 1e-28 at mid-flight) and returns at
2πR (the #166 focusing as absorber condition); (3) complete histories
leave no free remnants (pure positive-frequency far future = every
emission confirmed), a **well-posed 2×2 system** (condition number 1.0)
whose *unique* solution is `G_F = Ḡ + (i/2ω)cos(ωt) =
(i/2ω)e^{−iω|t|}` — **the Feynman propagator is the time-symmetric
field plus the absorber response of the closed universe**. **The two
orderings**: the θ(±t) split of G_F is offer and confirmation; each
half-line transform is an OFPT energy denominator — individually
regulator-dependent (finite-T truncation never converges, spread 0.50)
and non-covariant (Δ→−Δ violation ≥ 0.21) — and the coherent sum is
the *exact* covariant pole `−(ω−iε)/[ω(Δ²−(ω−iε)²)]` (7e-15 on a
random grid), exactly even in Δ. **The phase is forced**: the isometry
maps offer into confirmation (the same history read in the two time
directions), so relative phase 0; the deform test S_φ = I₊ + e^{iφ}I₋
breaks evenness and the pole form at O(1) for every φ ≠ 0 (0.81/0.59
at φ = 0.3 up to 2.0/3.9 at π), and the engine's `retro_phase_match`
already peaks exactly at the coherent point (its π branch = the #48
exchange sign, a distinct discrete sector, recorded). **The Compton
tie-in**: the two-ordering sums equal `1/(s−m²)` and `1/(u−m²)` exactly
across lab kinematics (1e-14; the u-channel entirely off-shell — the
"virtual" regime is the coherent pair, no on-shell intermediary), and
the deformed phase fails them by ≥ 7%: the denominators the arc
assumed are now derived, and the amplitude factors as denominators
(this PR: complete histories) × numerators/vertices (#37–#44, #46) ×
pole/magnitude/tensor geometry (#35/#45/#46). What QED postulates as
the iε prescription is, on this geometry, a theorem about complete
histories. Scope: per-mode, tree-level, free propagator; spinor
numerators/vertices not rederived; ε → 0 is the complete-absorption
idealization (on the closed bulk a physical inverse absorption time);
frozen geometry, no backreaction
(`compton_transactional_propagator_probe`, PR #213).

**The Compton-edge capstone: Release II (PR #211).** The closing PR of
the #204–#210 arc — the whole chain re-verified green in one run (as
#200 did for the previous arc), the register updated, plus **one new
parameter-free law** aimed at the arc's final unknown
(`docs/compton_edge_capstone.md`). **The Compton-edge law**: #210 left
the mass ladder with one O(1) — σ_mode/λ̄_C. Adding the mass term to
#202's bridge equation, `(ρ³φ′)′ = [k(k+2)ρ + m²ρ³]φ`, deforms the
exact suppression law **universally in x = σ_mode/λ̄_C alone** (the
m·r_s = 0.02 and 0.05 curves coincide to 2.6e-4; the massless limit
re-derives φ = σ to machine zero, #202's c₀(1) = 0): the law is
`ε₁ = (r_s/σ_mode)·D(x)` with D computed (0.966 at conv A, 0.831 at
conv B), and the **sensitivity S(x) = |d ln ε₁/d ln σ| = 1 exactly
below the Compton scale** (the #202 naturalness, its validity domain
now located) growing beyond — so **naturalness caps the O(1)**: S ≤ 2
confines x ≤ 2.58, the ladder's own worst-tolerated 4.48 confines
x ≤ 5.58, and **both convention anchors sit inside** (S = 1.07, 1.36).
The deformation-corrected self-consistent band tightens the #210
bracket to **σ_mode/λ̄_C ∈ [0.626, 1.307] — still bracketing 1**: the
successor derivation has a compact, derived target at the Compton edge
of the natural domain, not an unbounded O(1); the electron/muon ratio
is a fine-structure phenomenon with its correction controlled and small
exactly where the anchor sits. **The ledger, green in one run**: (1)
the commitment chain — a kicked *retarded* run with real potentials
conserves the norm (6e-14) and closes continuity (2e-3), the classical
mean-field channel keeps a product state at entropy 0 vs the quantized
pairwise operator's 0.026, the SN/BMV nulls armed (SN phase 4.8e-17
rad, BMV witness 0.79 rad where BAM predicts zero); (2) the topology
chain — the singlet from (I⊗T)|Φ⁺⟩ (fidelity 1.0, E = −cos(a−b) to
2e-16), the swapping law (Ψ_{a+b+c}, separable mixture 2e-17), the
enumerated LHV bounds CHSH = 2 / Mermin = 2 vs the Y-junction's
Mermin = 4 with exactly-empty pairwise marginals (negativity 0); (3)
measurement + mass — the k-odd dispersion identity (2e-16), a live
Stern–Gerlach Born check (within 0.0024 of cos²β), a live Kaup point
(M = 0.632), the α bands bracketing their targets. **The register
after the arc** (all named and bounded): derive σ_mode/λ̄_C within the
Compton-edge window; connect the #55–#58 R* to the bulk mass μ = r_s²;
the 5D pants nucleation (Sorkin class, rate); W-class reachability;
registration/irreversibility (radiative decoherence); standing
negative kept — the cosmological constant (#165). **The falsification
card**: two near-term nulls (SN-scale signatures; a BMV witness), the
m_e/m_μ = (3/7…1)·α prediction at the ×1.5 level with its O(1)
confined, the neutrino cards unchanged. Scope: full-GR no-signaling is
a completion argument; the CM-wave emergence is derived only in its
gravitational sourcing and entangled structure; the pinch/pants are
topological content; the α anchor is constrained, not derived;
equilibrium is a hypothesis with a mechanism throughout. **Release II**:
the program's deepest interpretive machinery is now theorems,
measurements, and bounded windows — not imports
(`compton_edge_capstone_probe`, PR #211).

**The strong-field core solve: the collapse reading refuted, the anchor
relocated (PR #210).** Executes #203's register target — the
strong-field core contraction — at the tractable and sufficient level:
the **static spherical GR family of the committed ψ–q structure**
(shooting; ground state by the 0→1 node transition; the turning-point
criterion brackets the collapse endpoint), benchmarked on the **Kaup
point to 0.2%** (M_max = 0.6327 vs 0.633; ω/m = 0.848 at criticality)
(`docs/strong_field_core_solve.md`). (1) **The #179 runaway is a
genuine GR instability**: the order field enters by adiabatic
elimination (κ = 0.005), contributing the non-convex
V_int = −(gσ²−a₀)²/4λ above threshold; measured — the family's maximum
mass sits AT the ordering onset (σ_c = 0.15) and the ordered branch
loses mass monotonically (0.542 → 0.295): switching on the core order
destabilizes the star — the strong-field endpoint of the committed
structure IS collapse, confirming the premise of #203's collapse
reading before the measurement that decides it. (2) **The
measurement**: at criticality, σ_mode/r_s = **2.45–5.83 (RMS) /
5.42–12.72 (R99)** across the potential class (Kaup / the committed
q-channel / a repulsive control) — universally O(few–10), a factor ≥ 7
short of the required 88.6 (conv A) and ≥ 16 of 206.8 (conv B);
structural, not parametric (the Buchdahl-type bound: at criticality
every length is the gravitational radius × O(few); σ/r_s = 88.6 means
compactness ~1e-2 — far from collapse). (3) **The verdict — the
pre-registered failure mode fires**: the collapse contraction is
r_q(weak)/r_s(true) = **O(1)** (the #203 weak-field band 4.6–6.7
overlaps the critical band 2.5–12.7), not 13–45: **the mouth-pairing
mechanism is refuted as the quantitative origin of the electron mass in
its collapse reading** — constructively (the cores really collapse;
collapse just cannot make them light); the smallness mechanism (#195
index protection, #201 multiplicative structure, #202 exact law)
survives, exactly as #203 stated. (4) **The relocation — the
primordial EM-capped throat**: in BAM the throat is the PRE-EXISTING
object (#168's regular 5D Tangherlini Killing horizon; bulk mass
μ = r_s² a geometric datum), and the repo's own #55–#58 anchor puts its
radius at the EM-capped equilibrium (U_EM/mc² = α/2): r_s ~ α·λ̄_C —
the classical electron radius. The required hierarchy becomes
σ_mode/λ̄_C = 88.6α = **0.647** (A) / 206.8α = **1.509** (B) — **the
convention band brackets 1**: σ_mode = λ̄_C and r_s = r_e meet the #202
law with **no new number**; the neck aspect **c = ln(1/α) + O(1)**
(4.484 = 4.920 − 0.436); **m_e/m_μ = (3/7…1)·α = [0.00313, 0.00730]
brackets the observed 0.00484**; and the #203 window is met by this
anchor (r_q(weak)/r_s = [13.2, 44.6] ∈ [13, 45]) — by the primordial
throat, not by collapse. Constrained, not derived (the O(1) band is
×2.3); but the mass ladder's remaining unknown is now the same α the
#184 thread protects — **the electron/muon ratio as a fine-structure
phenomenon; the mass ladder and the α thread are one thread**. Scope:
static spherical 4D GR + turning point (sufficient for the
adjudication; 5D criticality differs by O(1) — cannot rescue ×15);
adiabatic q; three-potential class + the structural bound; the α anchor
inherits #55–#58's scope. **The register**: "do full 5D NR" is replaced
by two sharper items — derive σ_mode/λ̄_C (the O(1)) and connect the
#55–#58 R* to the bulk mass μ = r_s²
(`strong_field_core_solve_probe`, PR #210).

**The measurement sector: pointer outcomes for the entangled sector
(PR #209).** Closes the last standing open of the entangled-sector
thread. #206–#208 derived the entangled STATES from bridge topology,
but their operational content rested on "Born statistics for internal
states at dBB grade" — and #198's equivariance covers spatial transport
only. The missing chain — internal state → spatial pointer branches →
position beables → Born — is delivered in three links
(`docs/measurement_sector.md`): (1) **the pointer coupling exists in
the committed structure**: winding couples minimally to the fiber
connection (winding = charge, #42–#44 — the KK gauge coupling); the
per-channel potential V_k = −2t_χcos(2πk/N − θ(x)) is exact on the
lattice, its **k-odd part** −2t_χ sin(2πk/N) sinθ an identity (1e-16) —
a connection-gradient region exerts opposite forces on opposite
windings: **the winding Stern–Gerlach is charge measurement by
deflection**, derived not postulated (live, from the raw dispersion:
k = +1 crosses to +28, k = −1 turned back to −23); (2)
**fiber-integrated equivariance** — the #198 theorem extends verbatim
to the multichannel wave (orthogonal internal channels, real
per-channel potentials): continuity residual 1e-4 on the live
Stern–Gerlach evolution, a 20 000-throat Born ensemble at sampling
noise through branch separation; (3) **Born statistics for internal
states, measured**: P(+) = cos²β to ≤ 0.003 across the sweep (branches
7σ apart) — the internal-sector Born rule is #198's spatial
equivariance routed through the committed coupling, not a new
postulate; **pointer permanence**: the empty branch's influence on
occupied-branch guidance dies with the Gaussian branch overlap
(1e-4 → 7e-9): **effective collapse, from geometry, without
collapse**. **The operational closing**: the #206-derived singlet +
local setting rotations + Stern–Gerlach branches at both wings + dBB
position beables: E(0,0) = −1.0000 exact; the four CHSH correlators
within 0.012 of −cos(a−b); **CHSH = 2.824 read out of the POSITIONS OF
CLASSICAL BEABLES** with setting-independent marginals (operational
no-signaling) — Bell violation as pointer statistics, with no imported
quantum rule anywhere: state (#206 topology), coupling (KK), statistics
(#198 extended), equilibrium (#198/#204). **The spatial sector**: the
pointer IS the spatial part of the pair wave in a measurement;
**positional EPR follows from conservation at nucleation** (the
#58/#200 C-conjugate pair: co-located birth, anticorrelated momenta):
Duan–Simon sum 0.53 < 2 — spatially entangled, structure from symmetry
(widths as inputs); and the **#205 guiding-without-gravitating split is
realized exactly here**: the empty pointer branch guides until
separation and never gravitates (back-action ~1e-17, nil). Scope: the
SG window = co-moving connection-gradient transit;
registration/irreversibility beyond branch separation (amplification,
radiative decoherence) not modeled; equilibrium hypothesis; nucleation
widths underived. **The entangled-sector thread closes operationally**;
program-wide remaining opens: the 5D pants nucleation, W-class
reachability, the strong-field NR target
(`measurement_sector_probe`, PR #209).

**The GHZ sector: multipartite entanglement is bridge valence
(PR #208).** Closes the first named open of #207. A matching cannot
make genuinely multipartite entanglement; GHZ requires a **junction
joining three mouths into one bridge** — the trousers (pair-of-pants)
nucleation (`docs/multi_mouth_bridge_ghz.md`). Delivered with a derived
no-go attached: (1) **charged GHZ is superselection-forbidden** — if
the doublet label is a conserved flux (winding = charge), the
Y-junction has **zero** conserving channels in the ±1 doublet
(enumerated: all sums ±1, ±3) and GHZ components straddle distinct
total-charge sectors — exactly QM's charge superselection; the pairwise
sector never felt this ((k,−k) is zero-sum), so **the multipartite
sector is where charge and spin part ways**: GHZ lives only in the
transported-frame (spin/Pin, #195/#197) label — a derived distinction.
(2) **The Y-junction, live** (one bulk junction fiber read by three
mouths on the #206/#207 lattice): exactly symmetric three-way
distribution (spread 1e-16); per-leg deck phases η_k = ±i as transport
demands; one shared variable (ratio spreads 1e-16); **leg-cut
consistency** — removing a leg collapses the Y to the #206 two-mouth
bridge (cut mouth receives 1e-33; the (A,B) pair phase equals the
composed two-leg law): valence 2 recovered as the junction's special
case. (3) **GHZ emerges with a holonomy law**: the N-frame embedding
W|k⟩ = |k⟩⊗T|k⟩⊗…⊗T|k⟩ is an isometry for n = 2, 3 alike; the
extracted state is GHZ at **fidelity 1.00000** with relative phase
obeying the **multipartite holonomy law φ = −π(s_B+s_C)/2** (verified
at (2,2)/(0,2)/(2,1): phases 0/π/π/2 exact) — the #206 pair law and the
#207 swapping composition extended: **one composition rule,
φ = −πΣs/2, at every valence**; genuinely tripartite: **3-tangle
τ = 1.0000** (pairwise concurrences 0). (4) **Mermin = 4.0000 from
pairwise-empty marginals**: the local bound enumerated exactly (all 64
deterministic strategies cap at 2); the Y-state reaches the algebraic
maximum 4 (beyond any pairwise mechanism — Tsirelson caps two-party at
2√2); every two-mouth marginal is **unentangled** (negativities 0,
pairwise CHSH exactly 2) — all correlation in the triple, none in any
pair: the exact opposite of #207's matching. **The valence ledger:**
valence 2 → Bell pairs + strict monogamy + swapping; valence 3 → GHZ +
empty pairs; holonomies select the state within the class; charge
superselection prunes the hypergraph (flux labels: zero-sum channels
only; frame labels: every valence). Scope: the junction fiber stands in
for the 5D trousers cobordism (its Sorkin class and nucleation rate vs
the #58/#200 pair channel unsolved — whether tri-mouth nucleation is
dynamically *realized* is the physical successor question); W-class
reachability via networks + surgery named open; the
spatial/measurement sector is now the **only** standing open of the
entangled-sector thread (`multi_mouth_bridge_ghz_probe`, PR #208).

**Entanglement swapping is bridge surgery: linking never-co-nucleated
throats (PR #207).** Opens the dynamical half of #206's register
consequence. By #206's own theorem, throats that never shared a
nucleation cannot become Bell-correlated through any 3-space-local
dynamics — the mechanism that closes the gap is that **throats
interact**: a proximity pair pinches off in pair annihilation (the
#58/#200 topology-change channel run in reverse, returning the pair to
wave fronts), and the pinch **relinks the bridges** of two previously
unrelated pairs, leaving the distant survivors connected exactly as if
they shared a nucleation — entanglement swapping, as bridge surgery
(`docs/bridge_surgery_entanglement_swapping.md`). Quantitative, both
sides: (1) **the composition law is the quantum swapping law** — QM
side machine-checked in 16-dim (projecting mouths (2,3) of Ψ_a⊗Ψ_b onto
Ψ_c leaves (1,4) in Ψ_{a+b+c}, fidelity 1−1e-12, probabilities exactly
¼); bulk side measured on the 4-mouth lattice: φ₁₄ = φ_a + φ_b + φ_c
with the Bell-outcome phase **supplied by the junction holonomy of the
pinch** (the four-outcome orbit reproduced to ≤0.024 rad, state
fidelities ≥0.9998); (2) **winding superselection**: annihilation
conserves winding, so outcomes are confined to the charge-zero Ψ
sector, and the **unconditioned mixture over the four holonomies is
separable** (negativity 0, machine) — swapping's
classical-communication requirement and no-signaling (#204) in one
identity; (3) **the event, live**: two populated disjoint bridges; the
proximity pair (2,3) annihilates locally (fibers glue, binding wells
vanish); the distant (1,4) response to Alice's phase rises from 3e-6
(the ballistic tail) to 1.3e-2 — **×4650** — landing in exactly the
swapped winding channel (power ratio ~290) and carrying her phase
linearly (D(α) ∝ e^{iα}−1: magnitude ratio 1.99 vs 1.99, phase step
0.100 vs 0.100), while the middle bound population drains 0.50 → 0.29 →
0.20 (no-annihilation control ~0.54) into propagating waves — **the
annihilated pair returns to wave fronts, the linkage persists**; (4)
**the swapped pair saturates Tsirelson** (CHSH ≥ 2.8280 across the
orbit) between mouths that never shared a gluing; **repeater chains
compose associatively** (the a+b+c law applied twice = the total sum,
machine) — pairwise entanglement distribution over the whole throat
network is **mechanism-complete**; (5) **monogamy is the matching**:
bridges pair mouths (a perfect matching); partner negativity = 1/2
exactly, non-partner = 0 exactly (machine, three-bridge state) — the
monogamy of maximal entanglement is the statement that a mouth has one
bridge; surgery is the only rewiring move. Scope: the pinch modeled by
its topological content (full 5D pinch-off dynamics — the #58/#200
channel in reverse — not solved); the holonomy outcome distribution
rests on the equilibrium hypothesis (¼-each recovered at dBB grade);
internal sector; **multipartite (≥3-mouth) junctions — the GHZ
sector — are the named open construction**, alongside the
spatial/measurement sector (`bridge_surgery_swapping_probe`, PR #207).

**Configuration-space emergence: entanglement is bridge topology
(PR #206).** Attacks #198's condition 2 at its sharpest point, stated
as the target: **a single classical field on 3-space with local
dynamics and local readout is an LHV model and cannot violate CHSH —
Bell's theorem itself** (machine-checked: all 16 deterministic local
strategies enumerated, max CHSH = 2 exactly; the quantum sector needs
2√2). So BAM's entangled sector cannot emerge from field correlations
alone — it must come from the one nonlocal element the theory
legitimately owns: **the bridge**, whose two mouths are one object
through the bulk (#168's regular 5D Killing horizon; #169's J-quotient;
`bell.bulk_identity`'s kinematic reading). The concrete problem —
derive the effective ψ(x₁,x₂) of a throat pair from the universal
3-space field plus the bulk identification — is solved in three
measured/machine-checked layers
(`docs/configuration_space_emergence.md`). (1) **The identification,
live on a lattice** (a universal field on a 128-site ring × 8-site
fiber, local everywhere except the mouth cells glued through the bulk
as χ → s−χ): winding k=+1 into mouth A arrives at B as k=−1 with purity
**1.000000** (charge conjugation — the C-conjugate pair, Σc₁ = 0); the
inter-channel phase at B tracks A with slope **1 exactly** (one phase,
not two); the handle's **π fiber holonomy (s = 2) is the Bell-state
selector** (offsets the locked phase by exactly π — the singlet sign);
both mouths read the same shared mode (ratio spread 1e-8); cutting the
bridge kills everything (1e-30). (2) **The emergence lemmas:** the map
`W|k⟩ = |k⟩_A ⊗ T|k⟩_B` — one shared fiber read at two mouths — is an
**isometry**: the tensor-product structure is the two-frame description
of one bulk object (this is where configuration space comes from); with
the repo's derived non-orientable transport **T = iσ_y (T² = −1 — the
same Pin⁻ sign that forces Fermi statistics)**, the symmetric bridge
gives ψ_eff = (I⊗T)|Φ⁺⟩ = **the singlet**, equal to
`bell.bulk_identity`'s *postulated* pair state to fidelity
1.000000000000 — what that module assumed from topology is now
**derived** from field + identification (two independent paths, one
state; E(a,b) = −cos(a−b) matching to 2e-16). (3) **The quantitative
ER=EPR law:** sweeping the bridge preparation, the extracted
concurrence tracks the bridge-mode amplitudes to 2e-3 (**Schmidt
weights are bridge-mode amplitudes**; entanglement entropy = bridge
participation entropy) and **CHSH(ψ_eff) = 2√(1+C²) exactly**
(Horodecki, verified by direct settings optimization) — from the Bell
bound 2 at zero bridge coherence (and exactly 2 for the bridge-cut
product state) to Tsirelson **2√2** at the symmetric π-holonomy bridge
(extracted singlet fidelity 1.000000). **The nonlocality budget
closes:** marginals invariant under all Alice-side unitaries (1e-16) —
Bell correlation without a telegraph; statistics at dBB grade on
equilibrium (#198), whose signal-locality #204 measured; the physical
bridge is non-traversable — the identification is **imprinted at
nucleation** (the mouths bound one 2-handle, #200) and conserved
topologically (the lattice handle is a stand-in). **The register
consequence: #198's condition 2 splits** — its Bell-sector half (where
does non-factorizable structure come from?) is **discharged**:
entanglement is bridge topology, quantitatively; its dynamical half
(the spatial sector, N-body bridge networks/swapping, measurement
transport) remains the standing item, now sharply bounded. Scope:
internal (fiber/spin) sector; k = ±1 doublet (scalar reduction of
#195/#197); N_χ = 8; equilibrium hypothesis
(`configuration_space_emergence_probe`, PR #206).

**The SN-phenomenology audit: the committed Φ[ρ] meets the laboratory
(PR #205).** The #204 follow-through. Until #204, Φ[ρ] was an internal
modeling device; #204's no-signaling verdict is a statement *about* it,
so the SN-type self-coupling is now a **commitment of the theory** —
and with it BAM enters the reach of laboratory tests (macromolecule
interferometry, levitated optomechanics, gravitational-entanglement
witnesses): the program's **nearest-term falsification channel**, far
closer than the neutrino-sector cards, now on the register. The lab
predictions hinge on the sourcing subtlety #198 flagged and never
resolved ("a throat's own self-field is not a pilot wave for itself"):
**universal-wave sourcing** (empty interferometer branches gravitate at
~M/2 → full Schrödinger–Newton phenomenology) versus
**conditional/actual-configuration sourcing** (the mass rides the lump
through one branch → SN signatures absent). The two give different lab
predictions; #204 forces the pick, in writing — and this audit makes
the pick the only honest way: **by running the experiment on the
committed dynamics** (`docs/sn_phenomenology_audit.md`). (1) **The
beamsplitter adjudication:** the relaxed 1D throat-soliton (μ_c =
−10.5, ordered core) sent at a barrier transmits or reflects **as a
whole** — max(R,T) ≥ 0.95 at 7/8 velocities, co-occupation confined to
one grid point (bracket (0.5, 0.6)), deep-sub-threshold leakage into
the untaken branch ~1e-4 — while the **linear control** (same profile,
binding off) co-occupies both branches across the entire sweep
(max(R,T) ≤ 0.82): whole-body transport *is* the self-binding. The
regime mapping: the sandbox threshold sits at kinetic/binding ≈ 0.01; a
300 m/s molecular beam sits at ~5e-13 — eleven orders below: **the
mass-carrying field never co-occupies interferometer arms**; splitting
the mass is the relativistic/QFT regime (#200's pair-creation domain).
The pick, made by the equations: effective-level sourcing is
**conditional** — #198's aside, now measured. (2) **The channel is real
and exactly calculable:** a hand-prepared 50/50 branch split contracts
at exactly the field-equation Newtonian rate (measured/predicted =
1.000–1.005 to merger, zero fit) — the SN signatures scale as f², and
BAM's co-occupied fraction f is exponentially zero in the lab regime:
**SN signatures predicted null**. (3) **The classical channel cannot
entangle:** the BAM mean-field Φ[ρ₁+ρ₂] keeps a product two-throat
state at entanglement entropy **0 (machine)** while the
quantized-gravity pairwise operator at the same coupling entangles it
to S = 0.148 — BAM predicts **strictly null** in
Bose–Marletto–Vedral-type witnesses (phase 0.79 rad for the standard
proposal scale, where quantized gravity predicts a witness). (4) **The
SI confrontation:** existing data exclude nothing — the f=1 SN phase at
the Fein-2019 record (2.7e4 amu) is 4.8e-17 rad, the record sits ~1e5
below the SN inhibition scale m*(σ) = (ℏ²/2Gσ)^{1/3} ≈ 3–6e9 amu, and
ω_SN(Si) ≈ 0.05 s⁻¹ is beyond current optomechanics. **The register
gains two near-term nulls, both actively hunted:** the **SN-null**
(detection of SN signatures refutes the committed sourcing) and the
**BMV-null** (an observed gravitational entanglement witness refutes
classical Φ outright) — the first BAM predictions addressable by living
experimentalists. Scope: 1D structural/regime adjudication (3D repeat a
follow-up); the emergence of the CM pilot wave (how the empty branch
guides without gravitating) remains the standing open item, narrowed
(its gravitational sourcing now fixed); the committed signature is the
combination SN-null + BMV-null + fringes intact at every mass; if the
throat sector ever requires f > 0, these bounds propagate back as
constraints (`sn_phenomenology_audit_probe`, PR #205).

**The nonlinear no-signaling audit: the Gisin edge, faced by
construction (PR #204).** The second named frontier item after #203 —
the no-signaling edge of the #200 register's "nonlinear measurement
theory" (#198 condition 2). The stakes: BAM's pilot equation is
genuinely nonlinear (superposition defect `O(1)` measured vs `2e-14`
for the stripped linear control), and by Gisin/Polchinski nonlinear
quantum evolutions generically signal superluminally; the linear
no-signaling theorem cannot protect BAM, and the Gisin proof (which
uses the projection postulate BAM doesn't have) cannot automatically
convict it — so the audit is run **by construction** on the live
dynamics (`docs/nonlinear_no_signaling_audit.md`). The edge **fired
where it should, and the theory survives it**: (1) **the channel is
exactly the gravitational field** — a local unitary kick at A produces
a density response at B (separation 30) at `5e6×` the kinematic floor
in the Newtonian model (the weak-field model *does* signal, stated
plainly), but the response is O(G) (ratio 2.14 under G→G/2) and
**clamping Φ to the no-kick history collapses it back to the floor**
(suppression `6e6`; the order field q is local and carries nothing);
(2) **retardation confines it to the cone** — replacing Poisson by the
causal `□Φ = 4πGρ` (the gauge-fixed weak-field form of the 5D Einstein
equations whose Bianchi structure #199 verified), the response is
*machine-floor* quiet (`~1e-15`) outside the light cone, the front
scales as d/c (`c·t_front = 25.2/26.4/27.2` at c = 8/12/16 vs geometric
distance 25), and c→∞ recovers Newton — the superluminality belongs to
the **approximation, not the theory**; (3) **the completion is free** —
the retarded Φ is still a *real* potential, so #198 equivariance
survives untouched (norm exact, continuity residual `1e-4`, a
20 000-throat Born ensemble at sampling noise through the kicked
retarded evolution): no-signaling and the Born rule hold
*simultaneously*; (4) **the entangled sector obeys the same
structure** — on the two-throat state ψ(x₁,x₂) the x₂-marginal is
invariant under Alice's local unitary to machine precision (`3e-15`,
the partial-trace theorem on the discrete flow), equilibrium dBB
trajectory ensembles are signal-local (KS 0.015 ≈ noise 0.018) while
**non-equilibrium ensembles signal** (KS 0.045 — the Valentini
boundary, reproduced on the BAM transport: equilibrium is load-bearing,
exactly as in dBB), and the BAM mean-field gravitational dressing
violates marginal invariance at O(G) instantaneously (the Gisin
channel, exhibited) but **confined behind the front** when retarded
(arrival 1.9 vs geometric 2.0; machine floor before). Scope: minimal
causal completion (full GR constraint analysis not run); the
configuration-space pair is the #198 *effective* description — the
register item is **narrowed** (no-signaling edge audited), not closed
(the emergence of configuration space from the single 3-space wave
remains); equilibrium is a hypothesis with the #198 relaxation as
mechanism — non-equilibrium signaling is a prediction shared with the
whole dBB program (`nonlinear_no_signaling_audit_probe`, PR #204).

**The coupled 5D+soliton solve: the confrontation (PR #203).** The
final register item of the mass-ladder thread — and the refutation edge
**fired at the weak-field level**, cleanly. With no knobs left (the
#202 law exact, the #180 solution locked), the coupled solve gives
`σ_mode/r_core = 4.6–6.7` (pairing definitions) versus the required
`88.6`: **m_e over-predicted by ×13–19**. The structure of the failure
is itself informative: the direction is right — the true 5D core is the
strong-field endpoint of the #179 runaway, smaller than the weak-field
q-core, so the weak-field value is an **upper bound** on m_e (which
holds) — and the gap is **provably not closable inside the weak-field
model**: the binding sweep moves the ratio the wrong way
(`RMS/r_q = 3.17 → 1.53 → 1.25` as `Φ(0)` deepens `−2.6 → −6.9`; the
`ψ² > ρ_c` core compacts slower than the wave, structurally). The
thread — #192 → #193 → #194 → #195 → #201 → #202 → #203 — therefore
closes onto **one falsifiable number**: the NR core contraction
`r_q(weak)/r_s(true) = 13–45×` (convention/definition band), with the
failure mode stated: an O(1) contraction refutes the mouth-pairing
mechanism as the quantitative origin of the electron mass (the
smallness mechanism — index protection, multiplicative structure,
naturalness — survives; its numerical anchor does not). No new fits
anywhere; `Φ(0) = −4.2` at the locked point independently strains the
weak-field label — the strong-field solve is not optional
(`coupled_5d_soliton_solve_probe`, PR #203).

**The 5D throat-core solve: the exact suppression law (PR #202).**
#201's fitted mouth-coupling exponent becomes a law. On the t=0 slice
of the J-quotiented Tangherlini throat the k-winding problem is 1D on
the bridge `ρ(σ) = √(r_s²+σ²)` (proper-distance identity verified to
1e-16), and: (1) the deck parity forces the **Pin-twisted boundary
condition** — odd-k modes have a **node at the cross-cap**
(quotient-spectrum equality to 1e-5; #195's forbidden one-mouth lift,
geometrized), even-k are Neumann, and the k=0 (vacuum/boson) channel is
exactly flat — only the uncharged mode touches the identification
locus; (2) the **exact E≈0 radial law**: for k=1 the regular mode is
exactly `φ(σ) = σ` (closed form), general k grows as `σ^k`, and a
genuinely shallow bound state's tail *is* the regular solution
(pointwise, 2%) — the power-law regime is the physical one; (3) the
**suppression law** `ε_k = (r_s/σ_mode)^k·e^{c₀(k)}` with computed
constants `c₀ = {0 (exact), −0.405, −0.783}` — #201's `e^{−kc}`
derived, `c = ln(σ_mode/r_s)`. Consequences: the electron mass is a
**linear readout of the throat-to-mode scale hierarchy** — the
sensitivity to the physical parameter is `−1` exactly, closing the
naturalness chain `74.7 (dialed, #194) → 4.48 (log-parametrized, #201)
→ 1 (physical, #202)`; the inversion is exact (`c₀(1) = 0`):
`σ_mode/r_s = 88.6` (conv A) — the throat core is ~1% of the wave
extent, and **one dimensionless number, governed by an exact law, now
separates the program from an outright m_e/m_μ prediction** (the
coupled 5D+soliton solve); the pairing-impossibility bound hardens
(`ε₃/ε₁ ≈ 2e-4` vs the 88.6 a pairing hierarchy would need). Scope:
rigid vacuum Tangherlini, E≈0 regime, `σ_mode` from the soliton sector,
scalar reduction (`throat_core_5d_solve_probe`, PR #202).

**The Dirac-tower mass ladder: un-dialing the electron (PR #201).**
Executes the first item of the #200 register. The electron level is
rebuilt as a **multiplicative chain** on the #195 index-protected zero
mode — `m_e = ε₁·o₁·S₁` with the bare level *exactly* zero
(Atiyah–Singer), the one-mouth additive lift forbidden (angular
momentum), and inter-sector mixing superselected: the surrogate's
`t₁₃ = 14.85`, whose near-cancellation was the #194 dial, has **no
counterpart** — nothing in the chain can cancel. Heavy sector: the
surrogate's natural 2×2 {k=3,5} block, uplift honestly refit after the
k=1 decoupling (`β → 0.79·50π`). Fitting μ/e fixes
`ε₁ = (7/3)/206.77 = 0.0113` (convention stated, band carried), and the
fitted number **is O(1) geometry**: a WKB neck aspect `ℓ/a ≈ 4.5–5.3`,
and the #185 overlap kernel on the actual #180 soliton gives
`K(R*) = ε₁` at `R* ≈ 3.9–4.4 soliton radii` — constrained, not derived
(the 5D core solve remains the register item). **The three fine-tuning
diagnostics collapse:** Barbieri–Giudice `74.7 → 4.48` (= the neck
aspect; all heavy entries ≤ 0.89; zero sign flips in 2000 ±25% draws vs
flipping at ±2%); Berger sensitivity `−70.9 → +4.15` with **no
λ_break** on `(0,∞)` (the #192 signature — sensitivity = inverse
distance to a zero crossing — eliminated because there is no crossing);
Monte Carlo `P(m_e ≤ observed)`: `7.7% sliver → 49.6% typical`, smooth,
no cliff. **The impossibility bound:** pairing masses scale
`ε_k = e^{−kc}` (decreasing in k); μ/e from pairing alone would need
`ε₃/ε₁ = 88.6` — a tunneling amplitude growing with the barrier —
impossible: the inter-generation hierarchy provably cannot be mouth
pairing and stays with the dynamical uplift (fitted, natural per #194).
After the rebuild every sensitivity in the ladder is O(few) or below:
the ladder is fully natural — smallness from geometry (index +
pairing), ratios from dynamics (the phase budget), the #193/#197
division of labor realized in the mass model itself
(`dirac_tower_mass_ladder_probe`, PR #201).

**PR #200 — the pair-creation cobordism, constructed: the capstone.**
The one open construction of the #196 adjudication is closed: the BAM
pair-creation cobordism is exhibited explicitly —
`M = (S³×I) + H²₍₊₂₎ + H²₍₋₂₎`, two 2-handles attached along a
two-component unlink with framings ±2, boundary
`S³ → L(2,1) # L(−2,1) = RP³ # RP³` — **pair creation from the actual
closed S³ background** (integer surgery on an unknot = L(n,1); split
link = connected sum). Machine-checked: linking matrix `diag(2,−2)`
with `|det| = 4 = |H₁(RP³#RP³)|` and Smith form `diag(2,2)` (coker
= Z₂⊕Z₂ exactly); **even framings ⟹ spin cobordism** (the structure the
Pin⁻-framed SSC theorem requires); signature 0; **only index-2 Morse
points** — causally continuous, avoiding the proven-discontinuous
classes {1, 3} (Dowker–Garcia–Surya): the Sorkin selection rule
*selects* pair creation, the channel the #58 threshold always assumed,
and the ±2 mirror structure *is* the C-conjugate throat–antithroat
pair (Σc₁ = 0). With this, the exchange-statistics chain runs end to
end with no unconstructed step — 5D Einstein equations → RP³ geon →
SSC hypotheses → the explicit spin, causally-continuous 4-manifold →
SOH abelian weights → the forced Pin⁻ lift `−I` → **Fermi, Pauli** —
and **"Pauli from GR + the forced Pin⁻ framing" is a constructed
theorem**. The capstone probe also re-verifies the release ledger green
in one run (the #183/#196 algebra, the #193 scalar ladder, the #195
index zero mode, the #197 Dirac ladder, the #198/#199 guidance-current
identity) and encodes the open-items register (the Dirac-tower mass
ladder; the nonlinear measurement theory; the 5D core; the standing
negative results #165 and #160; the Borde–Sorkin remainder). The
defining claim at #200: the program's deepest steps are theorems with
hypotheses, not imports
(`pair_creation_cobordism_capstone_probe`, PR #200).

**The guidance law from the 5D bulk (PR #199 — discharging #198's
condition 1).** The #198 Born rule was conditional on the guidance
identification `v = ∇S`; it is here **derived** from the bulk field
equations (`docs/guidance_law_from_5d.md`). The chain: (1) the throat
mode is the fiber-winding KK mode `Ψ = ψe^{ikχ}` (#83/#193); (2) the de
Broglie current **is** bulk stress-energy — `T_{μχ} = k·Im(ψ*∂_μψ)`
identically (verified to `1e-12`; winding = charge, #42–#44); (3) its
conservation is the χ-component of the **contracted Bianchi identity**
of the 5D Einstein equations — `∇·G ≡ 0` verified symbolically and
*exactly* on the weak-field KK metric with arbitrary Φ(t,x), all five
components — so the continuity equation behind the Born rule is a
Bianchi identity of the bulk, not a model assumption; (4) the throat is
the quantized, topologically conserved unit of winding (#178/#181/#182)
and must ride its own conserved current: `v = 𝒥^i/𝒥^0 = ∇S` —
demonstrated **pointwise**: quantized 2D cores transported through the
live nonlinear dynamics move with the ambient `J/ρ` at the core
(background and partner-induced parts both reproduced; winding exactly
±1 throughout). GR taken whole **selects the Bohmian flow**: the
transport differs from geodesic motion exactly by the quantum potential
Q, which sits inside the same wave stress tensor (Madelung balance
verified at `1e-3`; the quantum force dominates ×27 on the live
dynamics) — geodesic transport is precisely the wrong flow #198 showed
destroys equivariance. The law is species-universal (k cancels in the
ratio). The chain now runs: **5D Einstein equations → Bianchi →
conserved fiber current → topological transport → `v = ∇S` → (#198)
`|ψ|²` equivariant and unique → Born rule** — the dBB-grade
interpretation as a consequence of the bulk field equations,
conditional only on the #198 equilibrium/measurement conditions and the
throat-core scale identification (the full 5D core dynamics is
constrained by its charge, not solved)
(`guidance_law_from_5d_probe`, PR #199).

**The Born-rule equivariance test (PR #198 — the deepest import
replaced by a theorem).** Item 3 of the theorem-shaped program. For a
classical-field-foundational theory the one honest route to the Born
rule is **equivariance** (Dürr–Goldstein–Zanghì): is `|ψ|²` the measure
preserved by the BAM transport flow, or does the ψ–Φ–q dynamics
equivariantly transport some other functional? Anything else would
refute BAM at its foundation. The answer
(`docs/born_rule_equivariance.md`): **Theorem 1** — the real-time
Hamiltonian flow of the #179/#180 functional (`i∂ₜψ = −½∇²ψ +
[Φ + ½gq²]ψ`, self-consistent Poisson gravity, live real order field)
satisfies `∂ₜρ + ∇·(ρ∇S) = 0` **exactly**, because every BAM coupling
enters as a *real* potential — the nonlinearity of Φ[ρ] and q is
irrelevant; verified on the live dynamics (residual `3e-4`; a
20 000-throat Born ensemble stays at sampling noise through a
two-soliton collision, KS ≤ 0.009 vs noise 0.007). The throat transport
is the local momentum field `v = ∇S` (Galilean structure; the exact
Ehrenfest relation verified on the nonlinear dynamics; the unique
current closing continuity). **Theorem 2** — `∂ₜh + ∇·(hv) =
(h−ρh′)∇·v`: only `h ∝ ρ` is equivariant (compressible flow, verified;
`√ρ`/`ρ²` ensembles and wrong flows all fail demonstrably).
**Falsifiable:** a dissipative deformation `iγWψ` produces exactly the
predicted continuity source `−2γWρ` (0.1% match) and breaks the
transport — the repo's own imaginary-time *relaxation* flow is of
precisely this non-equivariant type; only the Hamiltonian flow carries
the Born rule. **Relaxation:** a uniform ensemble falls 3.48 nats in
the coarse-grained H-function toward Born (the subquantum H-theorem on
the BAM dynamics) — fixed point *and* attractor. The label: the Born
rule enters BAM at **dBB grade** — equivariance + uniqueness +
relaxation, the same epistemic status as in Bohmian mechanics —
conditional on the guidance identification (throat velocity = phase
gradient; its 5D derivation is its own program) and the linear
measurement regime (test throat in an external pilot wave). The
`bell`-module import is now labeled: derived at dBB grade with stated
hypotheses (`born_rule_equivariance_probe`, PR #198).

**The analytic Berger–Dirac ladder (closed-form spectral geometry,
PR #197).** Item 2 of the theorem-shaped program: the Dirac spectrum on
the Berger-squashed S³ derived in closed form (Peter–Weyl; the spectrum
classical — Hitchin 1974, Bär — re-derived self-contained in
`docs/berger_dirac_analytic_ladder.md`) and applied to the {1,3,5}
ladder for the **actual spinor field content**. The spectrum: family A
`a_j = (2j+1)/λ + λ/2` (the winding tower), family B `b± = λ/2 ±
2√((j+½)² + m′²(λ⁻²−1))` — machine-validated against the assembled
operator (`1e-15`), the round spectrum `±(3/2+n)` with exact
multiplicities `(n+1)(n+2)`, the λ→0 collapse limit (the S²(½) Dirac
spectrum), and **Lichnerowicz** (`scal = 8−2λ²`: every zero lies at
λ > 2; the first harmonic spinor — the Hitchin phenomenon — located at
λ = 4). The odd-k ladder is the winding tower `m_k(λ) = k/λ + λ/2` with
**uniform gaps 2/λ**: ordered and gapped at every λ, absolutely
protected on the squash side, with O(1) round-point sensitivities
`(½−k)/(k+½) = −1/3, −5/7, −9/11` — the infinitesimal-squash refutation
is **analytically excluded**, and the #192 surrogate's fine-tuning has
no counterpart in the true spinor spectrum (closing the #194/#195
loop). Every crossing is located in closed form: character changes at
`λ×(k) = √(2k+4)` (√6, √10, √14), harmonic-spinor masslessness at
`λ*(k) ≈ 5.668, 8.035, 9.851` — under extreme stretch **the electron
sector collapses first**. The k₅ = 5 cutoff has **no spectral
counterpart at any λ** (the gap is independent of k; the collapse
boundaries grow with k — stretch removes sectors from the bottom, never
truncates the top): the three-generation count is dynamical (the
#122/#136 phase budget), now as a closed-form statement on the whole
family. Round ratios 7/3, 11/7 — O(1); the hierarchy stays dynamical.
Scope: S³ cover, unique spin structure; on RP³ the odd-k sectors are
the Pin-twisted modes and everything descends (the deck map is a fiber
translation, an isometry of every Berger metric)
(`berger_dirac_analytic_ladder_probe`, PR #197).

**The geon-statistics adjudication: is the exchange −1 a theorem?
(PR #196 — a correction and a sharpening.)** Careful mathematics against
the Sorkin-school literature (Sorkin–Surya gr-qc/9605050; Dowker–Sorkin
gr-qc/9609064; Hendriks; Giulini), applied to the actual BAM topology:
the #169 quotient makes the throat prime the **RP³ geon** (the
antipodally identified Einstein–Rosen neighborhood; one-sided RP²
cross-cap slice) — the canonical example of that literature. **What is a
theorem:** the throat passes all three Dowker–Sorkin hypotheses — prime
(elliptic, irreducible), non-chiral (`q² ≡ −1 mod 2`; an explicit
reflection descends), abelian (`π₁ = Z₂`) — so the spin-statistics
**correlation** is a theorem for pair-created BAM throats. **The
correction (to #170/#171):** RP³ is **non-spinorial** in bare Diff
(Hendriks: cyclic-π₁ primes; lens spaces) — the 2π rotation is *isotopic
to the identity*, so the sentence "the single geon's 2π rotation acts as
−I (spinorial)" is false as a bare-diffeomorphism statement; in bare
metric GR the SSC selects **Bose**, and frozen-topology canonical GR
leaves statistics a sector choice (the two-throat mapping class group is
`(Z₂∗Z₂)⋊S₂`: four scalar sectors plus a continuum of
indefinite-statistics sectors). **What survives, sharpened:** in
**Pin⁻-framed GR** — the framing BAM's own construction forces (#169
non-orientable RP² slice; #170 Pin⁻ unique; #195 pin-spinor modes) — the
trivialized rotation lifts to **−1** on the pin bundle (the SU(2) lift of
the `π₁(SO(3)) = Z₂` generator; RP³ ≅ SO(3)), the throat is spinorial
*with framing*, and the same SSC theorem selects **Fermi** — the #188
holonomy correctly reinterpreted as the pin lift. The honest label:
**"Pauli from GR + the forced Pin⁻ framing"**, not "Pauli from GR". The
#58 nucleation channel passes every bordism-level check (mirror-pair
identity by amphichirality; `Ω₃^{SO} = 0`; RP³ spin and `Ω₃^{Spin} = 0`);
the explicit BAM pair-creation 4-manifold (the Dowker–Sorkin construction
transplanted from ℝ³ to the S³ background) is the honest open
construction. Full argument: `docs/geon_statistics_adjudication.md`;
machine-checked arithmetic: (`geon_statistics_adjudication_probe`,
PR #196).

**The index mechanism: a Pin/Dirac zero mode for the k=1 sector.**
Answers #194's mechanism question YES, with no new ingredients. BAM
throats are Pin⁻ (#183/#188), so the throat mode is a **spinor**; on the
#193 sector reduction it is a charged spinor on the base S² (monopole
`q = k/2`), and the Atiyah–Singer index pins exactly `2q = k` chiral zero
modes in winding sector k — verified via the SUSY decomposition
`D²₊ = L_{q−½} − (q−½)`, `D²₋ = L_{q+½} + (q+½)` through the validated
#193 monopole solver: sectors {1,3,5} carry {1,3,5} zero modes
(residuals ≤ 1e-7), the opposite chirality gapped at `2q+1`, towers
matching the exact Dirac spectrum. The k=1 electron level is **zero by
topology**: under a flux-preserving gauge wobble the zero energy is
certified pinned in `[0, 6e-10]` (the wavefunction deforms, the energy
does not); under metric deformation of the base, ker D is conformally
rigid (`9e-16`, Ω-independent) while the scalar ground on the same
metric moves by O(ε) — 8 orders of contrast, energy-pinned vs
energy-tuned; the count changes only with a flux quantum (1→3→5). The
natural mass: a one-mouth first-order lift is forbidden by angular
momentum (no opposite-chirality `j = q−½` partner); the lift requires
the throat's **two mouths** (±k winding, opposite chirality — the BAM
wormhole supplies the Dirac pair): `|E_e| = ε·o` with `o = 1.000` —
linear, multiplicative, sign-stable, 't Hooft-natural. The #194
"geometric identity tying the diagonal to the repulsion" is the SUSY
factorization `D² = A†A`, available only to the spinor — the surrogate's
dialing is an artifact of treating a spinor problem with scalar
dynamics. Scope: mechanism established; the mass ladder is *not*
re-derived — rebuilding it on the Dirac tower with ε computed from the
throat overlap machinery (#185/#190) is the follow-up
(`k1_zero_mode_index_mechanism_probe`, PR #195).

**Attacking the fine-tuning: the electron near-zero — stabilized or
dialed?** The follow-up to #192/#193: is the surrogate's electron
near-cancellation (`E_e = 6.8754 − 6.6758 = 0.1996`, a 2.9% residue
between the k=1 diagonal and the transport repulsion) protected by
anything? Every candidate is tested and **excluded**: no chiral/sublattice
conjugation (`tr H = 736 ≠ 0`), no spectral reflection, no index-like
structured zero-mode (best overlap `0.917`), and E_e flips sign under a
±2% transport change — a cancellation, not a sign-stable seesaw. The
tuning is quantified and localized: Barbieri–Giudice `Δ(E_e)` up to 57
while every `Δ(E_μ), Δ(E_τ) < 1` (only the near-zero is tuned); the one
codimension-1 dialed combination is 76% transport vs 42% base action +
41% slope (global `Δ = 74.7`), and its fiber-map contraction `+71.1`
reproduces the #192 Berger λ-sensitivity — the two probes see the same
dial. The Monte Carlo null (20 000 samples, ±25% log-uniform) gives
`P(|E_e| ≤ observed) = 0.077` vs the linear-measure estimate `0.060`,
flat through zero — no attractor, no repulsion: exactly as rare as
generic. Origin: the **calibration imports the hierarchy** — fitting
`μ/e = 206.77` with an O(10)-scale matrix forces `|E_e| = E_μ/206.77`;
the surrogate carries the hierarchy problem, it does not solve it (and
the #193 operator has no near-zero at all — the tuning lives in the
instanton dynamics, where the dialed direction points). Numerology
guardrail held (det-zero roots vs round constants NOT matched; transport
root vs 8π: 1.6%, rejected). The sharpened question: a real solution
needs new structure pinning a k=1 zero mode — an index, a chiral grading
of winding sectors, or a geometric identity tying the 2π base action to
the transport repulsion (`electron_near_zero_naturalness_probe`,
PR #194).

**The field-theoretic odd-k ladder (the actual wave operator on the
Berger sphere).** The follow-up #192 promised: the genuine SU(2) Berger
Laplacian sectored by Hopf-fiber winding `k = 2m` replaces the instanton
surrogate — no ingredient map. The sector grounds are the closed-form
Kaluza–Klein split `E_k(λ) = 2k + (k/λ)²`, **derived** from the spectrum:
the `(k/λ)²` fiber term is the #83 unified-mass-operator throat winding
term, and the `2k` base part is the charge-`q = k/2` **monopole**
zero-point on the base S² (the winding is the charge — the #42–#44
Hopf⟷charge geometry; half-integer q for odd k = the Pin-twisted monopole
bundle), verified by an independent Wu–Yang finite-volume solve to
`~2e-7`. The {1,3,5} ladder is **absolutely protected**: `E_k ≥ 2k ≥ 2`
and gaps `> 4` in closed form for every `λ ∈ (0,∞)` — no λ_break — and
the deck grading `(−1)^k` is λ-independent because the antipode lies *on*
the Hopf fiber (a fiber translation, an isometry of every Berger metric):
the #183 algebra realized spectrally at every λ. But the same operator's
mass ratios are pinned O(1) at every λ (`ω₃/ω₁ ∈ [1.53, 3.0]` vs observed
μ/e `206.8`; `ω₅/ω₃ ∈ [1.25, 1.67]` vs `16.8`), so the hierarchy is
**not kinematic**. With #192 the claim is bracketed from both sides:
structure from kinematics/topology (absolutely protected — the surrogate's
1.4%-squash electron break has no operator counterpart), hierarchy from
the instanton dynamics (the metric-fine-tuned near-cancellation). Scope:
scalar operator, throat = winding sector, one deformation axis
(`field_theoretic_odd_k_ladder_probe`, PR #193).

**The spectral deformation test (upgrading #183 from algebra to
spectrum).** #183 proved the odd-k {1,3,5} sector is protected by
metric-independent algebra (deck determinant, `½trT²=−1`, odd parity) —
but algebra is not spectrum. Rebuilding the locked lepton Hamiltonian's
geometric ingredients on the Berger-squashed S³_λ (the #165 SU(2)
machinery; fiber-riding ingredients × λ, connection/base ingredients
fixed, the map declared before the sweep): the {1,3,5} structure — three
positive ordered levels — survives a **finite** window `λ ∈ (0.986, ≥3]`
with smooth, linear-at-round-point mass-ratio deformation, so the ladder
does **not** break at infinitesimal squash and the protection claim
upgrades from algebra to spectrum. The discovery the algebra could not
see: the electron level is a metric-fine-tuned **near-zero** (`0.1996` in
action units vs μ `41.26`, τ `694.98`) that crosses zero at a 1.4% fiber
squash (`λ_break = 0.98598`), and the μ/e log-sensitivity `−70.9` equals
`−1/(1−λ_break) = −71.3` — the steepness *is* the proximity to the
spectral boundary, in every ingredient map (default, flipped, minimal) —
while τ/μ is gentle (`+0.82`). So the topology guarantees three
generations; the round metric tunes the e–μ hierarchy — the protected
structure and the metrically-tuned numbers are now separated by
measurement. Scope: the locked instanton-transition surrogate under a
declared fiber/base map, not a first-principles wave operator on S³_λ
(that is the follow-up); action units, electron-calibrated ratios
(`odd_k_ladder_spectral_deformation_probe`, PR #192).

**Dynamic two-throat exchange path with back-reaction.** The exchange `−1`
of #188 (the adiabatic Pin⁻ holonomy) is the **adiabatic limit of a real
dynamical swap**: the throat's internal Pin spinor driven around the swap loop
`n̂(s)=(cos2πs, sin2πs, 0)` at finite speed `1/T`, coupled to a back-reacting
field `ẍ+γẋ+ω²x=−κ⟨σ_z⟩` that the moving throat sources and that acts back
through `g·x·σ_z`. The loop's exact Berry phase is `−½·2π=−π` (the exchange
`−1`); as the swap slows the *dynamical* geometric phase `φ_geo→−π` and the
non-adiabatic excitation `P_exc→0`, while at finite speed both deviate by an
explicit `O(1/T)` cost (`deviation×T≈29.7`, constant) that vanishes
adiabatically. The back-reacting field is sourced at finite speed (peak energy
`≈0.028` at `T=20`, → `1.5e-4` slow) and acts back on the spinor, but leaves
the adiabatic `−1` **unchanged** — the dynamics add a quantified,
adiabatically-vanishing cost, not a change of sign. Honest scope: an effective
model (the internal spinor + one field mode); the field-resolved real-time
two-throat solve (orbitals translating, the #190 Coulomb field re-solved each
step, the swap geometry in 3D) is the follow-up; weak-field, code units
(`dynamic_two_throat_exchange_probe`, PR #191).

**Moving-mouth Berry phase.** _Closed_ (`spin_wigner_rotation_probe`,
PR #60). The Hopf-holonomy result `A_φ = ½ cos χ` (`∮A = π cos χ`)
reproduces the relativistic **Wigner rotation** from two non-collinear
boosts (composed in `SL(2,C)` to an `SU(2)` rotation), with the same
`½` factor / spinor double cover (`2π → −1`) / "rotation = ½ × solid
angle" law. The boosted throat is a genuine relativistic spin-½
particle. This closes the static-spin-½ → dynamic-particle bridge in
the strongest possible form: not just `2π → −1` under a closed loop,
but the full Wigner-rotation law under composed boosts.

**Two-throat Coulomb limit.** _Run and passed (`two_throat_coulomb_probe`)._
Two wormhole mouths at angular separation `ψ` on `S³` interact through
the `S³` Green response. The probe confirms `F(ψ) ∝ 1/sin²(ψ)` as the
leading form, reducing to `F ∝ 1/r²` (inverse square) and `V ∝ 1/r`
(Coulomb) in the small-separation flat limit, both to machine
precision. The exact force carries a compact-`S³` modulation
`N(ψ) = (π−ψ) + sin ψ cos ψ` — the antipodal-image correction — which
makes the force vanish at the antipode; Gauss's law `Φ(ψ) = Q_enclosed(ψ)`
holds exactly on the closed manifold. This was the cleanest single test
of "EM from `S³` geometry" and the place the program was most exposed to
falsification; it is now a demonstrated result. See
`docs/two_throat_coulomb_research_plan.md`.

**Quark `β = N · π/2` with `N = 466` and `n_part = 233`.**
_Scope sharpened across a 6-PR arc; absolute derivation outside the
current closure-ledger scope._

The original five-probe sequence (origin → boundary → decomposition
→ audit → sub-block stability; `docs/quark_beta_status.md`)
localized the irreducible structural piece: across all 12 logged §8
ablations the only preserved invariant is `N_q ≡ 0 (mod 2)`. The
structural reading is `N_q = 2 · n_part`, the factor of 2 topological
(the Z₂ partition multiplicity), `n_part` phenomenological.

PR #76 (`quark_npart_origin_probe`) extended the candidate catalog
to Fibonacci, Lucas, Padovan, Perrin, tribonacci, color × flavor ×
generation, QCD β₀, and Tangherlini QCD-shell mode counts — finding
only baseline coincidences (`F_13 = 233`, `9·k_5²+k_5+3 = 233`) and
no enumeration that survives §8 drift. The structural diagnosis:
**the v3 quark Hamiltonian is lepton-shaped machinery** (basis
`{(k=1,±), (k=3,±), (k=5,±)}` — the same odd-`k` throat-traversal
modes that give the lepton ladder), but quarks live in the QCD shell
channel (#68–#69). `n_part = 233` is the empirical price of fitting
QCD-confined quarks on closure-quantum throat basis vectors.

PR #97 revisits this with the now-complete lepton and neutrino sectors
and sharpens the diagnosis. The implicit worry behind PR #76 — that the
~9-order quark mass² hierarchy is simply too large for the geometric
closure machinery — is overturned by the neutrino arc, which derived a
comparable hierarchy (the keV→TeV seesaw `M_R = m_D·e^{S}`, ~10⁶ in mass)
as a *clean geometric exponential* (the tortoise bounce, an O(15)
action). So size is not the obstruction. The program now has two
geometric hierarchy types — charged leptons (the closure-ledger ladder
with the §8-stable integer `4·k_5² = 100`) and neutrinos (the bounce
exponential). The quark hierarchy, by contrast, is **irregular**: the
consecutive up-type ratios `m_c/m_u ≈ 588` and `m_t/m_c ≈ 136` are not
constant (so not a clean exponential), and the down-type ratios differ
(so not a single power law). Irregularity across a wide scale range is
the signature of renormalisation-group running (`α_s` logarithmic), so
the quark masses are QCD-dressed dynamically — the quark sector is the
mass program's **one dynamical (non-geometric) hierarchy**. This is why
the quark closure integer is the only one that §8-drifts (216–255): it
absorbs dynamical content no geometric closure quantity encodes, and the
lepton↔quark gap `N_q − N_lepton = 466 − 100 = 366` quanta is precisely
that dynamical (QCD) excess. The compensator verdict stands, sharpened: a
geometric closure integer can only *compensate* a dynamical hierarchy,
never derive it — and the right route is a QCD-shell model *with* `α_s`
running.

PR #98 takes the first step on that route and, in testing it, sharpens
the mechanism. QCD's mass anomalous dimension `γ_m` is flavor-universal,
so under renormalisation-group running every quark mass is multiplied by
the *same* factor and quark mass *ratios* are RG-invariant — `α_s`
running sets the overall scale, not the hierarchy. So the irregular
hierarchy is not a running effect at all; it is the **flavor puzzle** —
the Yukawa couplings, free inputs in the Standard Model, hierarchical and
irregular for reasons no current theory derives. The quark Yukawas
overflow the compressed shell-overtone capacity (`ω(1,n=3,4,5)` spans
only ×1.49 in mass) by ~×5×10⁴, which is precisely why `n_part`
compensates; the charged leptons (also a flavor puzzle, `m_τ/m_e ≈ 3477`)
instead fit the winding ladder (`k ∈ {1,3,5}`, PR #71) that *has* the
dynamic range. So BAM captures the quark sector's structure — six quarks
(`3×2`), the Z₂ up/down partition, the `k=0` shell channel, three
generations from `k_5` — geometrically, while the Yukawa magnitudes are
the flavor puzzle, open across all of physics. The #97 core (the quark
hierarchy is dynamical / non-geometric; `n_part` can only compensate it)
stands; the mechanism is refined from "QCD-RG running" to "the flavor
puzzle".

**QCD confinement: Cornell potential and the flux-tube string tension
(PR #99).** Where the quark *masses* are the (non-geometric) flavor
puzzle, the quark *confinement* is geometric in BAM. The Cornell static
energy `V(L) = σ·L − A·ℏc/L` reads, term by term, as BAM geometry: the
linear `σ·L` is a flux tube — a 1D wormhole-bridge connecting the
quark–antiquark with constant energy per unit length (the defining
property of a confining string), and the Coulomb `−A·ℏc/L` is
short-distance one-gluon exchange (the QCD analogue of the lepton Coulomb
law derived from eigenmode throat flux). The flux tube breaks by Schwinger
pair nucleation `Γ ∝ exp(−π m_q²/(σL))` — the QED Schwinger form
`exp(−π m_e²/(eE))` with the electric field replaced by the string
tension, `eE → σ`. This is precisely the pair-production / throat-pair
nucleation of PR #58 (`e E_S · R_MID = m_e c²`) transported to QCD: the
string snaps when its work `σ·L` reaches the pair threshold `≈ 2 m_q`. So
QCD string breaking and lepton pair production are the *same* BAM
nucleation physics with `eE ↔ σ`. The BAM string tension reproduces the
Nambu–Goto Regge slope `α' = 1/(2πσ) = 0.884 GeV⁻²` (observed ~0.88–0.93)
and the lattice string-breaking length (`L ≈ 1.4 fm` vs ~1.35), with
`√σ ≈ 0.42 GeV` — the confinement (Λ_QCD) scale — as the single
dimensionful anchor, the B4 analogue (lepton `m_e = ℏc/R_MID` ↔ QCD
`√σ`). The confinement *form* is geometric; only the one scale is
anchored.

**Glueballs: a pure-confinement benchmark, and where BAM's topology
diverges (PR #100).** Glueballs — closed flux loops with no valence
quarks — are the cleanest confinement probe: no quark masses, untouched
by the flavor puzzle, mass set entirely by `σ`. The BAM orientable
closed-loop ground state `√(4πσ) ≈ 1.50 GeV` (3.5√σ) benchmarks the
lattice 0++ `√σ` scale (4.1√σ) to ~13%, and the closed-string glueball
Regge slope is half the meson — both parameter-free given `σ`. The
BAM-specific content is topological: the framework carries *two*
closed-loop sectors, orientable (the glueball ring, periodic) and
non-orientable (the Möbius tube, antiperiodic). The Möbius half-twist
makes the modes antiperiodic — half-integer rather than integer — so the
non-orientable glueball tower is shifted by `πσ` in `M²` and interleaves
the orientable one, effectively doubling the glueball spectrum.
Orientable-string lattice QCD has no such sector. Crucially, glueballs
are *not experimentally observed* — they mix with ordinary qq̄ mesons of
the same `J^PC` and have never been cleanly isolated — so the Möbius
tower is a legitimate BAM-vs-lattice difference for a non-observable:
testable against lattice (which can isolate pure glue), but contradicted
by no experiment. This is the first place where BAM's non-orientable
topology is *expected* to predict differently from QCD, precisely where
nature has not yet ruled — a feature of the program, not a defect.

**The Möbius / exotic sector, where the topology meets data (PR #101).**
Pursuing the non-orientable topology into the *open* flux networks turns
the glueball logic on its head: hybrids and multiquark exotics, unlike
glueballs, *are* experimentally observed, so here BAM's topology must
match. It does. BAM's flux-network topology is the hadron taxonomy —
meson (open tube), baryon (Y-junction), tetraquark and pentaquark
(multi-junction networks), hybrid (tube + excited/twisted flux), glueball
(closed loop). The key BAM-native statement is that the *exotic* J^PC
quantum numbers are the signature of a non-orientable flux tube: an
ordinary orientable qq̄ meson is restricted to `P = (−1)^{L+1}`,
`C = (−1)^{L+S}` and so cannot have `1-+` (or `0--, 0+-, 2+-`), whereas a
Möbius (antiperiodic) flux tube carries the phonon that opens exactly
those channels. The lightest observed exotic hybrids — `π₁(1600)` and
`η₁(1855)`, both `1-+` — match this both in quantum numbers and in mass:
the Möbius/hybrid excitation gap is one flux-tube quantum, `≈ 2√σ ≈
0.85 GeV`, placing them at `ρ + 2√σ ≈ 1.62 GeV` and `~1.0 + 2√σ ≈
1.85 GeV`. The observed tetraquarks (`X(3872)`, `Z_c`, `T_cc`) and
pentaquarks (`P_c`) fill the multi-junction network types. The same
non-orientable Z₂ that gives the throat its spin-½ (PRs #63–#67) is the
half-twist that marks these exotics; the Möbius baryon remains a
BAM-specific prediction. So the exotic sector is the first place BAM's
non-orientable topology confronts data — and passes.

**Baryonic exotics: the most-constrained corner, and why it lacks a
smoking gun (PR #102).** The baryon analogue is subtler and more
exposed. For mesons the Möbius twist gave a manifestly exotic `J^PC`
(`1-+`, C-forbidden); baryons have no good `C` and `P = (−1)^L`,
`S ∈ {½, 3/2}`, so *every* half-integer `J^P` is reachable by an ordinary
qqq — there is no forbidden, exotic baryon `J^P`. A BAM Möbius / hybrid
baryon therefore carries ordinary quantum numbers and shows up only as a
*supernumerary* state, an extra resonance beyond the quark-model count.
Its natural mass is one flux-tube quantum above the ground baryon,
`nucleon/Δ + 2√σ ≈ 1.79, 2.08 GeV` — squarely in the light N*/Δ* region,
the densest and best-measured part of the entire hadron spectrum. So the
baryonic exotics are the *opposite* extreme from the glueballs: not BAM's
freest topological prediction but its most experimentally constrained
one. The Möbius doubling (a Z₂-twisted partner per state) must either
coincide with observed-but-unexplained resonances — filling the "missing
resonances" the quark model under-predicts — or decouple from `πN`
production (the standard missing-resonance mechanism); unmitigated
over-prediction would be excluded. Ranked by data density the channels
run light N*/Δ* (strongest constraint) → strange hyperons → charm/bottom
baryons (the freest, where a clean new state is most likely findable).
The honest scope: with no smoking-gun `J^P`, this is a *counting*
prediction, testable only against the dense spectrum or via decoupling —
the sharpest, but also the least clean, of BAM's non-orientable tests.

**The heavy-quark Möbius baryon: a concrete, findable prediction
(PR #103).** The freest channel of PR #102 — the heavy-quark baryons — is
where BAM can make a clean prediction, and heavy-quark symmetry supplies
the missing handle. With the heavy quark a near-static spectator, the
Möbius / flux-tube excitation lives entirely in the light/flux sector, so
its gap is the flavor-INDEPENDENT flux-tube quantum `Δ = 2√σ ≈ 0.85 GeV`
— the same above the charm and bottom ground baryons. That
flavor-independence is the heavy-sector signature that replaces the
absent exotic `J^P`: a supernumerary state at the *same* ~0.85 GeV above
both the charm and the bottom ground baryon. The concrete predictions —
`Λ_c ~3.14`, `Ω_c ~3.54`, `Λ_b ~6.47`, `Ω_b ~6.89`, `Ξ_cc ~4.47 GeV` —
sit above the ordinary orbital tower (`Λ_c` P-wave at ~+0.31 GeV, the
`Λ_c(2940)` at ~+0.65 GeV) and just above the currently-measured
excitation ceilings, so they are unexplored, within LHCb / Belle II
reach, and not excluded. The doubly-heavy `Ξ_cc` and the `Ω_b` have no
measured excitation spectrum at all — entirely unconstrained. So across
the non-orientable hadron sector BAM runs the full gamut of testability:
the unobserved glueballs (freest), the observed mesonic `1-+` hybrids
(matched), the densely-constrained light baryons (a counting test), and
now the heavy Möbius baryons (a clean, correlated, findable prediction in
the sparse heavy spectrum) — the exact mass (lattice hybrid gaps span
~0.8–1.3 GeV) and `J^P` remaining open.

**How the heavy Möbius baryon decays — and how to find it (PR #109).** A
mass prediction is only half a discovery program; the other half is the
decay pattern that tells the state apart from an ordinary excitation. Here
the topology does the work. The Möbius excitation is the non-orientable
(orientation `−1`) flux sector, while the ground heavy baryon is orientable
(`+1`), so to decay the half-twist must **unwind**, shedding the stored
`2√σ ≈ 0.85 GeV` as light isoscalar hadrons — a hybrid de-excitation with
the heavy quark a spectator. That mechanism inherits the flux-tube model's
**hybrid selection rule**: a hybrid cannot decay into two ground-state
(both-S-wave) hadrons, so the Möbius baryon's naive single-S-wave-pion
transition to the ground state is *suppressed*, while `Σ_Q π` (the spin-1
light diquark), the coherent isoscalar S-wave dipion `Λ_Q(ππ)` (like
`ψ(2S) → J/ψ ππ`), and P-wave-baryon + π are *preferred*. An ordinary
radial excitation does the opposite — single pion to the ground state — so
the branching **pattern**, not the mass, is the discriminator. The
cross-flavor handle sharpens further in the decays: because both the gap
`2√σ` and the light-meson thresholds are flavor-independent, the all-light
release energies are *identical* for charm and bottom — `Λ_Q ππ` at
`Q = 569 MeV` and `Λ_Q η` at `Q = 301 MeV` — the same dipion spectrum above
both ground baryons, with the `Σ_Q π` channels offset only by the small
`Σ_Q − Λ_Q` hyperfine splitting (167 MeV for c, 194 for b). Honesty about
the cost: with several open channels at `Q ≈ 0.5 GeV` the state is *broad*
(lattice hybrid widths run ~tens–150 MeV), so it is best resolved in LHCb
and Belle II amplitude (Dalitz) analyses of `Λ_Q ππ`, `Σ_Q π`, and the
open-flavor `D N` / `B N` channels (`Ξ_cc` and `Ω_b` wide open), not as a
sharp peak. The absolute branching fractions, the total width, and the
`J^P` need the flux-tube decay amplitudes and remain open; what BAM
delivers is the branching pattern and the cross-flavor Q-structure — a
falsifiable search program, not just a bump.

**One page for the experimentalist (PR #110).** With the non-orientable
sector now spanning mesonic `1⁻⁺` hybrids, glueballs, and heavy Möbius
baryons with their decays (PRs #100–#109), the useful consolidation is a
compact *experimental note* an LHCb / Belle II / BESIII analyst can read
off — predicted masses, Q-values, preferred/suppressed modes, and analysis
handles, with every number a pushforward of the single confinement scale
`√σ`. It collects the matched mesonic hybrids (π₁ ~1.62, η₁ ~1.85 GeV, the
exotic `1⁻⁺` smoking gun), the unobserved `0⁺⁺` glueball at `√(4πσ) ≈ 1.50
GeV`, the heavy Möbius baryon masses (Λ_c ~3135 … Ω_b ~6894 MeV), and their
twist-unwinding decays (single-pion-to-ground suppressed; the cross-flavor
Q-match at 569 / 301 MeV), with the open items — exact masses in the
~0.8–1.3 GeV band, branching fractions and widths, and the baryon `J^P` —
carried forward unchanged. It adds no physics; it makes the sector's
existing content searchable on one page.

**Sharper still: a tiered search table (PR #114).** Where #110 is a
reference card, #114 turns the heavy-baryon channels into a ranked,
actionable LHCb / Belle II search table, with one genuinely sharper handle.
The preferred twist-unwinding channel `Λ_Q(ππ)` has a dipion invariant-mass
endpoint `m(ππ)_max = M_Möbius − M_ground = 2√σ ≈ 849 MeV` that is
*flavor-independent* — the same edge above the charm and the bottom ground
baryon, with the spectrum peaking high (the coherent isoscalar S-wave, like
`ψ(2S) → J/ψ ππ`). That is sharper than a Q-value: a fixed edge in a
directly-plotted observable, so a single overlay of the charm and bottom
dipion spectra tests the whole picture at once. The table is tiered by
discovery feasibility — Tier 1 the `Λ_c` + `Λ_b` cross-flavor pair (highest
yield, golden `pK⁻π⁺` reconstruction, and together the Q-match clincher);
Tier 2 the entirely-unexplored doubly-heavy `Ξ_cc` and `Ω_b` (a clean bump
would be discovery, but production is rare); Tier 3 the calibratable `Ω_c`,
sitting above the well-mapped 2017 excitations — with the suppressed
single-pion-to-ground branch and the cross-flavor Q-match as the
discriminators. It adds no physics beyond #109/#110; it makes the search
concrete and ranked.

The four-PR QCD-shell arc (PRs #77–#80) develops the right machinery
quantitatively. The user's physical reframe: **quarks do not pass
through the throat; they are the wavefronts that resolve the cavity
itself.**

  - **PR #77 (`qcd_shell_waveguide_scaffold_probe`).** 6-state
    `(l, n, p)` basis where `l` = S³ Casimir, `n` = shell-saturated
    radial overtone (n ≥ 3 for l=1), `p ∈ {+, −}` = Z₂ partition.
    Operator scaffold `H = H_kin + H_Z2 + H_couple` with `H_kin =
    ω²(l, n)` cavity-eigenfrequency-squared — NOT the lepton
    `β·k²·(2π)` winding cost. Matches PR #69's 3 × 2 = 6 flavor
    count.
  - **PR #78 (`shell_mass_ordering_audit_probe`).** Shell basis is
    structurally better than v3 in four ways (cavity wavefronts; ω²
    kinetic; Z₂ partition slot for within-generation inversion; 6
    flavors). Uniform `χ·σ_z` cannot reproduce the inversion (best
    2/3 blocks); sign-flipping `χ_n` can (existence proof). Coverage
    gap: shell kinetic ×2.2 vs observed ×6.4·10⁹; `n_part` not
    resolved at PR #78 alone.
  - **PR #79 (`boundary_stress_chi_n_probe`).** `χ_n = T_odd(n) =
    (T_inner − T_outer)/2`, the Z₂-antisymmetric piece of cavity-
    mouth boundary stress under PR #63's inner/outer swap. NO free
    parameter once cavity geometry is fixed. Uniform-positive sign
    (no flip), shell-suppressed magnitude — 30–100× too small for
    observed splittings. PR #78's sign-flipping ansatz is
    structurally overruled. Singlet projector added as placeholder.
  - **PR #80 (`color_algebra_shell_probe`).** **BAM-native color
    algebra = `SU(2) × Z₂`** (SU(2) from B2 / Hopf holonomy in
    PRs #59–#66; Z₂ from PR #63's inner/outer swap). 4 generators vs
    SU(3)'s 8 — substantive structural difference. SU(3) **not**
    derivable from current scaffold: all natural triplet candidates
    (3 generations from `(k_5+1)/2`, three Hopf fibrations of S³,
    S³'s SO(4) isometries, Hopf U(1), bulk 5D structure) give
    SU(2)/SO(3) algebras. v3 species map revised: `+ = heavier`
    uniformly. `n_part` re-audit: even with large illustrative
    couplings, the full Hamiltonian's eigenvalue range factor stays
    single- to two-digit, while observed inter-generation mass²
    range is ~6.4·10⁹. **The inter-generation hierarchy is outside
    the scope of any BAM color algebra acting on the shell basis.**

**Arc closure.** What closed: shell basis is the right machinery;
`χ_n` has a no-free-parameter structural origin; BAM-native color
algebra is identified; v3 species map is settled. What remains open:
the inter-generation mass hierarchy and `n_part = 233` as a residual
phenomenological compensator with sharply identified scope (it
absorbs the hierarchy). The most plausible extension route is
**Pati-Salam SU(4)** with a BAM-native throat↔shell algebra map —
extending PR #68's structural transition into a quantitative
unification of the throat (lepton) and shell (quark) sectors.

**Pati-Salam bridge + mass-operator unification (PRs #82–#83).** The
throat↔shell bridge was built (PR #82): each generation has a lepton
at radial overtone `n = g−1` (throat) and a quark-pair at `n = g+2`
(shell), the shift `+3` being PR #68's shell-saturation threshold;
the unified 12-state `(l, n, p)` basis carries a throat-shell Z₂.
PR #82 identified three open extensions for a full SU(4): BAM-native
neutrinos, 3-fold quark color, and **lepton-quark mass-operator
unification**. The third — the deepest — is now **closed at the
structural-form level** (PR #83): the lepton `β·k²` (PR #71) and quark
`ω²(l, n)` (PR #77) mass operators are the SAME Bohr-Sommerfeld
operator

```
m²(k, n)  =  (k·2π / L_throat)²  +  ((n+1)·π / L_cavity)²,   L_throat = √(2π)/k_5.
```

The cavity sector is verified Bohr-Sommerfeld
(`∮√(ω²−V) dr* = (n+1)·π` to machine precision); the lepton sector's
`β·k² = (k·2π/L_throat)²` is exact and recovers `β_lepton =
(2π/L_throat)² = k_5²·(2π) = 50π`. **Leptons wind through the throat**
(`k ≠ 0`, closure quantum `2π`); **quarks resolve the cavity** (`k =
0`, closure quantum `π` per Bohr-Sommerfeld node). The whole
throat-vs-shell distinction collapses to the single winding quantum
number `k`, and the `2π`-vs-`π` channel quanta are the program's
pervasive full/half-cycle structure. This is genuine open work only
in two remaining respects: an independent derivation of the two
`L_eff` scales from one principle (the lepton `L_throat` re-expresses
PR #71's already-derived `β_lepton`), and the inter-generation
hierarchy (the cross-channel / mixed-mode question, still comparable
in scope to deriving the QCD hadron spectrum from geometry).

**Neutrino sector — the first of PR #82's three extensions (PRs
#85–#87).** The unified `(k, n)` operator splits the plane into four
quadrants; the chargeless `k = 0, n < 3` corner is the neutrino, and
the winding-and-saturated `k ≠ 0, n ≥ 3` corner the leptoquark (PR
#85). The neutrino quadrant gives the lightest states, but ~10⁵–10⁶
too heavy — until the BAM-native fix: `k = 0 ⟹ c₁ = 0`, so under `C`
(`c₁ → −c₁`, PR #63) the neutrino is invariant — **its own
antiparticle, necessarily Majorana** — and admits the seesaw
`m_ν = m_D²/M_R`, available *only* to the chargeless sector (charged
leptons, `c₁ = ±1`, are Dirac), which is exactly why only neutrinos
are anomalously light (PR #86). The seesaw scale `M_R ≈ 0.3–1.8 TeV`
is then grounded in the PR #58 nucleation channel (PR #87): a `ΔL = 2`
Majorana mass *is* a throat↔antithroat flip; `Σ c₁ = 0` on a single
state reproduces the only-neutrino rule; the scale is *not* the static
barrier height (`E_c ≈ 2.8 keV`, ~10⁸ too small) but the **tunnelling
amplitude through** the barrier, `M_R = m_D·e^{S}`, with a modest,
generation-stable bounce action `S ≈ 15–18`. The conceptual upshot is
that the seesaw scale is **reframed from a free ~TeV mass into an
instanton action**: the whole keV→TeV gap is carried by the single
dimensionless tunnelling exponent `S`, not by a new heavy particle. This
closes the first of PR #82's three extensions structurally; what remains
is `S` from first principles (the Euclidean throat-action / instanton
normalisation), which would promote the absolute `m_ν` to a prediction.

PR #88 builds that bounce explicitly and identifies it as the
**non-orientable tortoise logarithm**. The `ΔL = 2` flip runs along the
odd extension across the throat (`c₁ → −c₁`), and the 5D tortoise
coordinate `r* = r + (rs/2)ln|(r−rs)/(r+rs)|` diverges logarithmically
there. Two structural results follow: a perfectly **rigid throat gives
an exactly massless neutrino** (`ε → 0 ⟹ L* → ∞ ⟹ S → ∞`), so the
boundary compliance `ε` is the mass-generating parameter and the
smallness of `m_ν` is the near-rigidity of the throat; and the reduced
bounce `S = √(2 μ E_c)·L*(ε) ∝ ln(1/ε)` is naturally `O(10)` and
generation-stable — exactly the form `S` required. Honestly, though, the
electron-throat tension under-produces: `S ≲ 1`, some `~40×` short of
`~16`. Matching needs the `ΔL = 2` (B−L) throat tension to be `~6–12×`
stiffer than the EM-throat tension. So the open input is localised once
more — a mysterious `~TeV` mass (PR #86) → an `O(15)` instanton action
(PR #87) → an `O(10)` dimensionless tension ratio (PR #88).

PR #89 then constrains that tension ratio. Because the `ΔL = 2` flip
reverses the throat's orientation (`c₁ → −c₁`), it is a **global**
operation on S³, so `t` is a global-closure enhancement of the **local**
EM surface tension (PR #56) — not a free coupling. BAM has exactly two
fundamental action scales for such a closure, and they bracket `t`: the
**closure quantum `2π`** (a single great-circle orientation reversal —
the cheapest global flip; lower bound) and the **winding action
`√β_lepton = k_5√(2π)`** (a full throat winding to the antipode — the
costliest lepton-sector route; upper bound). Hence
`t ∈ [2π, k_5√(2π)] ≈ [6.28, 12.53]`, parameter-free — *exactly* PR #88's
required `6–12` (the computed `[6.41, 12.05]` sits inside). The `6–12`
band was not a fit but the BAM closure-to-winding window. The residual
freedom is reduced to a single number — *where in the window* — which is
the boundary compliance `ε` (the window edges map to `ε ≈ 6×10⁻⁷` at the
closure-quantum end and `ε ≈ 1.3×10⁻²` at the winding end); the
winding/cavity mass ratio `m_charged/m_D ≈ 11.9 ≈ √β_lepton` corroborates
the winding edge. So the open input has been localised four times —
`~TeV` mass → `O(15)` action → `O(10)` ratio → the closure-to-winding
window — leaving the compliance `ε` as the last undetermined number; an
`(t, ε)` degeneracy and the bounce normalisation are the honest caveats.

PR #90 closes the chain by deriving `ε` from the bulk throat geometry.
Near the neck the warp is `f(r) = 1 − (rs/r)² ≈ 2(r − rs)/rs`, so the
proper distance from the neck to `rs + ε` is `ℓ = √(2 rs ε)`, i.e.
`ε = ℓ²/(2 rs)`: the compliance is the throat's (neck-warped) **healing
length**. Crucially it is sub-throat *for the neutrino* and only the
neutrino: the charged-lepton throat (`c₁ = ±1`) is propped open by its
EM self-repulsion `A/R` at `R* ≈ R_MID` (and so cannot flip — it stays
Dirac), whereas the neutrino throat (`c₁ = 0`) has `A = 0`, so nothing
props its neck open and the bounce approaches it down to the bulk
healing length. The chargelessness that makes the neutrino Majorana is
the *same* property that makes its compliance sub-throat — and hence its
mass tiny; the smallness of `m_ν` is the unobstructed near-rigidity of
the chargeless neck. The natural BAM sub-throat scales (`R_c³`, `Δ³`,
`(m_D/m_charged)²`, `E_c`) all land `ε` inside the PR #89 window, and at
the winding-edge tension `t ≈ √β` — the edge the PR #89 mass-ratio
cross-check already favoured — the chain yields `S ≈ 15–19` and
`m_ν ~ few meV`, squarely the observed scale, with no input outside the
throat geometry (the `2π` edge gives `S ≈ 4`, too small, so the chain
closes only at the winding edge — the same one the cross-check picked).
So the entire chain — `~TeV` seesaw scale → `O(15)` instanton action →
`O(10)` tension ratio → closure-to-winding window → sub-throat healing
length → `meV` — is closed within BAM throat geometry at the
order-of-magnitude level: **the neutrino mass scale is geometric, not
tuned.** What remains is the precise `m_ν` and its generation spread (a
geometry-only `(t, ε)` gives a uniform `S`, hence `m_ν ∝ m_D` — a ×2.7
spread — whereas the observed `m_ν/m_D` spans ×18, calling for a
generation-dependent healing length or the mixing sector).

PR #112 presses on the one number this leaves implicit — can `ε` be
*computed* from the bulk compliance, or is it still being *read back* from
the meV scale it is meant to predict? The honest answer is a genuine
partial. The compliance is `ε = ℓ²/(2 rs)` with the neck healing length
`ℓ ~ R_c = 2σ/ρ`, and the surface tension `σ` and bag density `ρ` are fixed
by the *electron* rest-energy calibration (PR #58: `σ = 1/(12π)`,
`ρ = 3/(4π)`, so `R_c = 2/9`) — with no neutrino mass anywhere. The
candidate compliances are all sub-throat, `O(10⁻²)` (`R_c³ ≈ 0.011`,
`Δ³ ≈ 0.018`, `R_c²/2 ≈ 0.025`), so the *order of magnitude and sub-throat
character of `ε` are computable from bulk geometry alone* — and with the
winding-edge tension `t = k_5√(2π) = √β_lepton` the chain returns
`S ≈ 16.85` and `m_ν ≈ 2.1 meV`, the observed scale as an *output*. That
genuinely derives the neutrino's exponential lightness (sub-throat `ε ≪ 1`
⟹ large `S` ⟹ `m_ν = m_D e^{−S}` tiny) without using the meV scale. What it
does *not* do is fix the precise value: the bounce action is steep in `ε`
(`m_ν ∝ ε^{4.8}` at this tension), so the `O(1)` ambiguity among the
healing-length candidates (`R_c³`, `Δ³`, `R_c²/2`) already spreads `m_ν`
from ~2 to ~108 meV, and the absolute compliance normalisation is precisely
the bulk stiffness `κ₅²/Λ₅` that BAM never pins (only the dimensionless
`√6` RS combination is fixed, PR #57). So `ε` is rightly reclassified from
"inferred from the meV scale" to "bulk-geometric to order of magnitude":
**the smallness is derived from bulk compliance; the exact value remains a
residual**, of the same kind as every other dimensionless number tied to
the single anchor.

PR #91 takes up that spread and the mixing. The neutrino generations are
the cavity radial overtones `n`, so the bare prediction is **normal
ordering** with `m_ν ∝ m_D` (the cavity-floor ratios `1 : 1.87 : 2.74`).
The spread is widened in the right direction by the overtone-dependent
neck coupling: the bounce suppression grows with the throat-neck
coupling, which is precisely PR #79's boundary stress `χ_n`, and `χ_n`
*decreases* with `n` (0.304, 0.097, 0.039) — so higher-overtone
neutrinos are less throat-coupled, more compliant, less suppressed, hence
relatively heavier, lifting `m₃` toward the `Δm²`-implied value. The
headline result is the mixing dichotomy: the PMNS matrix is the overlap
of the charged-lepton mass basis (throat-winding, `k≠0`) with the
neutrino mass basis (cavity-resolving, `k=0`) — *different* channels of
the unified operator, generically strongly misaligned ⟹ **large PMNS** —
whereas up- and down-type quarks are *both* cavity-shell modes (`k=0`,
same channel) ⟹ nearly aligned ⟹ **small CKM**. So the long-standing
`PMNS ≫ CKM` puzzle is the BAM **cross-channel** (leptons: throat-winding
× cavity-resolving) vs **intra-channel** (quarks: shell × shell)
distinction. The spread direction and the mixing dichotomy are
structural; the precise spectrum (an `O(1)` coefficient; the absolute
scale unmeasured) and the explicit angles (the cross-channel overlap
integrals), and the CP/Majorana phases, are open.

PR #113 presses PR #91's spread mechanism for a quantitative prediction —
can a generation-dependent healing length `ε_n` driven by `χ_n` actually
reproduce the hierarchy, not merely point in the right direction? The
honest answer is no, and the reason is instructive. Taking the natural
law `ε_n ∝ 1/χ_n` (compliance is inverse stiffness), the direction is
indeed correct — `ε_n` rises with the overtone, so higher-`n` neutrinos
are less suppressed and heavier, giving normal ordering untuned — but the
magnitude overshoots wildly. The observed splittings require only a gentle
`ε_n` profile, ratios `(1, 1.18, 1.57)` across the three generations,
whereas `1/χ_n` supplies `(1, 3.13, 7.79)`; pushed through the bounce this
yields `m_ν3/m_ν2 ≈ 162` against the measured `5.85`, a factor of ~28 too
much spread (and orders of magnitude in the absolute masses). The culprit
is the very steepness identified in PR #112 — `m_ν ∝ ε^{4.8}` at the
winding-edge tension — which amplifies the factor-eight variation in `χ_n`
into four orders of magnitude in mass; the power that would fit the data,
`p ≈ 0.15` to `0.31` in `ε_n ∝ χ_n^{−p}`, is an inconsistent fraction, not
the principled unity. So a generation-dependent `ε_n` can *accommodate*
the spread by fitting a gentle profile but cannot *derive* it from the
overtone stress: the same bounce steepness that made `ε`'s absolute value
a residual now blocks the natural overtone variation from setting the
spread, which therefore stays a residual — plausibly the business of the
mixing/anarchy sector rather than the healing length.

PR #92 takes up the angles. A literal same-coordinate mode overlap turns
out to give *small* mixing — the cavity overtones are near-orthonormal
sinusoids, so a winding-imprint overlap is a near-permutation matrix
(mixing ≲ 5°). The largeness of PMNS is therefore not a literal radial
overlap; it is that the two lepton generation labels live in *different*
coordinates of the S³ × radial space — charged leptons in the
closure-winding `k = 1, 3, 5` (the Hopf fibre), neutrinos in the
radial-overtone `n = 0, 1, 2` (the cavity) — so the map between them has
no preferred alignment and is effectively **anarchic** (Haar-random in
generation space). This is the BAM realisation of neutrino anarchy.
Quantitatively, a Haar-random `U(3)` has angle medians `θ12 ≈ θ23 ≈ 45°`,
`θ13 ≈ 33°`, and the observed PMNS angles (33.4°, 49°, 8.6°) sit at the
~30th / 57th / 4th percentiles — broadly typical of anarchy — whereas the
CKM angles (13°, 2.4°, 0.2°) sit at the ~5th / 0.2th / 0.0th percentiles,
with a joint probability ≈ 0 of a Haar matrix being so aligned. So the
quark mixing is extremely atypical of anarchy — aligned — exactly as
expected when up- and down-type generations share the single
radial-overtone (shell) coordinate (intra-channel). The class-level
separation — PMNS anarchic (cross-coordinate), CKM aligned
(intra-coordinate) — is a firm BAM prediction matching observation; the
specific angles, being statistical, are not pinned, and `θ13` sitting on
the small side of anarchy (4th percentile) is the one mild tension.

PR #93 resolves that last tension. `θ13 = |U_e3|` is the *corner*
element of the generation/channel lattice — it links the lowest winding
(`k = 1`, the electron flavour) to the highest overtone (`n = 2`, the
heaviest neutrino), the most coordinate-distant pair (generation gap 2),
whereas `θ12` and `θ23` are adjacent (gap 1). Because the throat↔shell
coupling (the PR #82 `+3` shift, the PR #83 unified operator) is *local*
in the `(k, n)` lattice, the `g = 1 ↔ g = 3` corner is reached only by
*two* channel-hops, so `U_e3` is a suppressed two-hop amplitude — a
residual nearest-neighbour alignment of the otherwise-anarchic map. A
structured-anarchy model (corner variance `exp(−μ)`, `μ = 0` being pure
anarchy) with a modest `μ ≈ 3` shifts the `θ13` distribution down (median
33° → ~16°), makes `θ13` robustly the smallest angle (the observed
hierarchy `θ13 < θ12, θ23`), and moves the observed `θ13 = 8.6°` from the
4th to the ~21st percentile — relieving the tension — while `θ12` (~44th)
and `θ23` (~70th) stay typical. The suppression mechanism is robust; the
exact value (the residual-alignment strength `μ`; the `θ13` median
saturates near 14–16°) remains open.

PR #94 closes the phase sector. CP violation is *generic*: the winding
(charged-lepton) amplitudes carry the Hopf holonomy `e^{ikχ}` (PR #60,
the throat Berry phase `∮A = π cos χ`), so the cross-channel overlaps
that build the PMNS matrix are intrinsically complex, and `δ_CP ≠ 0, π`
with probability 1 — CP conservation (a real PMNS) is measure-zero, with
no BAM symmetry forcing it. The Jarlskog invariant
`J = Im(U_e1 U_μ2 U*_e2 U*_μ1)` mirrors the angle dichotomy: the observed
`|J_PMNS| ≈ 0.026` is typical of an anarchic `U(3)` (51st–81st
percentile, large CP violation), whereas `|J_CKM| ≈ 3×10⁻⁵` is extremely
atypical (~0.1th percentile) — aligned, CP-suppressed. And the *two
Majorana phases exist* precisely because the neutrino is Majorana
(`c₁ = 0`, PR #86): they are CP-violating phases of the ΔL=2
throat↔antithroat sector (the bounce of PRs #87–#90), observable in 0νββ,
where a Dirac neutrino would have none. The specific phase values, like
the angles beyond the dichotomy, are anarchic and not pinned. With this
the neutrino arc closes at the structural level — Majorana nature, mass
scale, ordering, mixing class, the `θ13` hierarchy, CP genericity, and
Majorana-phase existence are all BAM-native — leaving the precise
spectrum and the specific phases/angles as the (statistical /
one-parameter) residuals.

PR #95 collapses that structure into a single falsifiable observable —
the effective Majorana mass `m_ββ = |Σ U_ei² m_i|` measured in
neutrinoless double-beta decay. It combines the whole arc: 0νββ *occurs*
because the neutrino is Majorana (`c₁ = 0`, PR #86; a Dirac neutrino
would forbid it); the *normal ordering* (PR #91) selects the NO band of
`m_ββ` (≈ 1.5–3.7 meV at zero lightest mass); the *anarchic Majorana
phases* (PR #94) populate the whole band, including the cancellation
trough where the three terms partially cancel and `m_ββ → 0` (around a
lightest mass of a few meV); and the *light absolute scale* (PR #90)
places us there, giving `m_ββ ≲ 8 meV`. That lies below the current bound
(KamLAND-Zen, `m_ββ ≲ 28–122 meV` — so the null result is expected),
largely below next-generation reach (LEGEND-1000 / nEXO, ~9–20 meV), and
below the inverted-ordering floor (~19 meV). It is a sharp falsifier: a
0νββ discovery with `m_ββ ≳ 19 meV` would imply inverted ordering or a
quasi-degenerate spectrum, contradicting the BAM normal-ordering +
light-scale prediction. So the neutrino sector ends not merely
structurally complete but with a concrete experimental target for the
coming tonne-scale 0νββ searches — the exact `m_ββ` remaining a band
because the lightest mass is unmeasured and the Majorana phases are
anarchic.

PR #96 adds the cosmological companion from the *same* spectrum. The sum
of neutrino masses `Σm_ν = m1 + m2 + m3`, probed by CMB lensing and
large-scale structure, has a normal-ordering floor `√Δm²_21 + √Δm²_31 ≈
58.7 meV` (the inverted-ordering floor being ≈ 99 meV); the BAM light
scale (PR #90) keeps the sum pinned just above it, `Σm_ν ≈ 59–65 meV`,
out of the quasi-degenerate regime. This is consistent with Planck 2018 +
BAO (< 120 meV), just inside DESI DR1 + CMB (< 72 meV), and right at the
DESI DR2 + CMB frontier (~60–64 meV) — exactly where cosmology is now
probing. A robust cosmological `Σm_ν` below the normal-ordering floor
would exclude normal ordering (and sit in tension with the oscillation
`Δm²` themselves), while a quasi-degenerate `Σm_ν ≳ 100 meV` would
contradict the light scale. So `m_ββ` (≲ 8 meV) and `Σm_ν` (~60 meV) are
the two observables of one light, normal-ordered, Majorana spectrum — a
joint, cross-checkable pair that current and near-term experiments are
already testing.

PR #111 sharpens these meV-scale numbers from a band into a pinned
spectrum. Updating the oscillation inputs to the latest global fit (NuFIT
6.0, 2024) fixes two of the three masses outright — `m₂ = 8.65`, `m₃ =
50.34 meV`, so the normal-ordering floor is `Σm_ν = 59.0 meV` — and the
2025 DESI DR2 + CMB bound (`≲ 60–64 meV`) corners the lightest mass against
that floor, `m₁ ≲ 3 meV`, tightening the sum to `Σm_ν ∈ [59.0, 62.6] meV`.
The sharpest BAM statement is the hierarchical limit `m₁ → 0`, where the
laboratory effective masses follow directly: the β-decay mass `m_β ≈ 8.8
meV` and the 0νββ mass `m_ββ`, which in normal ordering has a *nonzero*
floor — the solar contribution `s12²c13² m₂ = 2.60 meV` exceeds the
reactor one `s13² m₃ = 1.10 meV`, so the terms cannot fully cancel and
`m_ββ ∈ [1.5, 3.7] meV` over the Majorana phases. The honest other half is
reachability: only `Σm_ν` is near-term testable — DESI is cornering it at
the floor now — while `m_β` sits ~4–5× below the best foreseeable β-decay
reach (Project 8) and `m_ββ` ~3–10× below next-gen 0νββ (LEGEND-1000,
nEXO). The spectrum is pinned, but the program's testable handle is
cosmological, not laboratory; and a 2025 caveat is worth flagging — some
DESI + CMB fits already prefer central `Σm_ν` at or below the floor, which
if it hardens is tension for every normal-ordered model, BAM included.

**QFT event reinterpretation: Compton scattering from BAM.** _Closed
at the analytic level._ An 11-PR thread (PRs #25–#35) constructed a
BAM amplitude for Compton scattering by progressively identifying
the BAM-native ingredients needed to reproduce Klein-Nishina. The
thread reaches a closed form at the resummation stage: the vertex
modification factor is

```
F²(x, c) = 4·x³·(x² + 1 − x·sin²θ) / [(1 + c²)·(1 + x)²]
```

with `x = ω'/ω = 1/(1 + ε(1 − cos θ))` and `c = cos θ`. The BAM
amplitude `f_BAM_baseline · F²` reproduces the full Klein-Nishina
differential cross section **exactly at all orders in ε up to
ε ~ 2** (highly relativistic Compton), to machine precision. The
perturbative results — γ = −3/2 at O(ε), the patterns `ν₀ = γ²` and
`ξ = −A_φ(0)` at O(ε²) — are Taylor expansions of this closed form.

What this establishes:

  - The antipodal `S³` Green function `G(ψ) ∼ 1/(4πψ)` reproduces
    the QED propagator pole exactly.
  - Two transverse photon polarisations on the `S³` tangent bundle
    reproduce the Thomson angular factor `(1 + cos²θ)/2`.
  - The closed-form F² resums all finite-energy corrections; no
    `ε·k` vertex contraction is needed (α = 0 in the parametric
    scans, at both O(ε) and O(ε²)).

What this does *not* yet establish:

  - Loop corrections. Tree level only; vertex/self-energy/vacuum
    polarisation would couple to the bulk radial channel.

Both of the open items flagged at the time of the resummation result
have since been addressed:

  - **First-principles BAM derivation of F²** — _now established._ The
    two-factor decomposition `F² = K(x)²·Q(x, c)` is derived from a
    single master functional on the warped-product internal geometry
    `M_int = C × S³` (`C` = radial cavity `[R_MID, R_OUTER]`):

    ```
    ℳ(ω; x, c) = G_C(r, r′; ω) ⊗ 𝒢_{S³}(Ω, Ω′)
    ```

    read three ways from one object — the mass spectrum as its ω-poles
    (radial ladder × S³ Casimir), `K(x) = 2x/(1+x)` as the throat-
    boundary dwell-time impedance series, and `Q(x, c)` as the S³
    Hopf-fibre helicity reduction. The vertex residue reproduces F² to
    machine precision while the poles give the masses — masses and the
    vertex from one functional. The `F² = K²·Q` factorization is the
    consequence of the product internal geometry (separation of
    variables), not a failure to unify. This closes the **B5′ residual**
    of the BAM effective-action scaffold; see
    `docs/bam_scaffold_status.md` and `docs/master_integral_research_plan.md`.
  - **Cross-process consistency** — _established._ The same closed-form
    F is process-general under crossing (Breit–Wheeler `γγ → e⁺e⁻`,
    pair annihilation `e⁺e⁻ → γγ`) and the thread extends to the
    4-fermion tree processes (Bhabha, Møller); see
    `docs/tree_qed_status.md` (PRs #36–#46).

### BAM effective-action scaffold — barrier closure

The first-principles programme was organised as a covariant 5D
effective-action scaffold with five mismatch terms (B1–B5), and is now
**complete**. Four are closed: **B1** (closure quantum `∮A = 2πn`) and
**B2** (antipodal `Z₂`, `T = iσ_y`) promoted to a topological/discrete
action sector (`RP³ + spin structure + winding θ-term`); **B3**
(hard-wall throat BC) forced by single-valuedness under `T² = −I`; and
**B5** (the 5D→4D reduction producing F²) closed by the master integral
above. The fifth, **B4** (the dimensional bridge `ℏ = m_e·R_MID·c`), is
not a gap but a **structural necessity**: the closure-ledger/Maslov
machinery is scale-free, so exactly one external dimensionful anchor is
mathematically required (B4 irreducible). That anchor is relocatable
from the particle mass to the **invariant bulk separation**
`ΔR = R_OUTER − R_INNER` — a proper, cosmologically fixed length (the
throat is a static bound vacuole, decoupled from Hubble flow; comoving
co-expansion is observationally excluded) — giving
`m_e = f_closure·ℏ/(ΔR·c)`, `f_closure = 0.52`. The full ledger is in
`docs/bam_scaffold_status.md`; the closure release note (through PR #53)
is in `docs/scaffold_closure_release_note.md`.

This is the strongest amplitude-level result in BAM so far. It
demonstrates that the Compton amplitude — historically the cleanest
non-trivial QED prediction — is reproducible from a BAM-native
construction (antipodal Green function + transverse Hopf-fibre
polarisation + closed-form vertex). With the master integral (above)
deriving `F² = K²·Q` from one `C × S³` functional and the tree-QED
thread (PRs #36–#46) extending the same primitives to BW, annihilation,
Bhabha and Møller, BAM's amplitude-level reach has extended from
"reproduces Compton" to "derives QED tree amplitudes from geometry".

## Open problems

The README's validation table is the authoritative status snapshot.
Several items deserve to be called out here as research-level open
problems rather than implementation TODOs:

- **Where does `ℏ` enter?** _Closure-cycle structurally complete;
  R_OUTER physically selected; SI conversion reduced to a single
  1.054 factor._ See `docs/hbar_origin_status.md` for the closing
  summary. An eight-probe sequence in `experiments/closure_ledger/`
  established three connected results: (1) **the closure cycle is
  integer-quantized in units of 2π for every species** —
  `N_total = N_layer_1 + N_radial` with all constituent channels
  (antipodal closure, Hopf-throat partnership at χ = 0, β-uplift
  closure quantum, hard-wall radial Bohr-Sommerfeld) integer-quantized
  individually; the Hopf-throat partnership and the hard-wall
  Dirichlet condition at the throat are both forced by `T² = −I`.
  (2) **The Compton-bridge geometry (ω(1, 0) = 1 exactly at
  R_OUTER ≈ 1.449) is physically vetoed** — re-running the locked
  lepton surrogate at this geometry breaks the muon and tau masses
  by ~46 %, and no β re-tuning recovers both species. The γ-lock
  geometry (R_OUTER ≈ 1.262, Σ V_max = 22.5, ω(1, 0) = 1.054) is the
  unique physical selection. (3) **The γ-lock R_OUTER is selected
  by a cross-species self-consistency loop**: bisecting μ and τ
  independently against `γ = Σ V_max(R)` gives the same R* to within
  0.008 %, confirming that the radial barrier-sum geometry
  reproduces both lepton mass ratios at a single R_OUTER. The
  remaining open piece is **the SI conversion factor 1.054** —
  ω(1, 0) at the cross-species fixed point. Whether 1.054 has a
  closed form in `(k_5, π, barrier invariants)`, and whether the
  R_OUTER fixed point can be lifted from "phenomenological-
  parameter-dependent" to "fully geometric" (current sensitivity to
  transport / resistance is ~1–7 %), are concrete next-pass targets
  identified in `docs/hbar_origin_status.md`. The B4 audit (PRs
  #52–#53) sharpens this: the 1.054 factor is *dimensionless* and
  orthogonal to the anchor (even a closed form would not supply the MeV
  scale); the only dimensionful input is a single length, irreducible by
  scale-freeness and relocatable to the invariant bulk separation `ΔR`.
  Predicting ℏ in SI is therefore gated solely by the value of that one
  geometric anchor.
- **Self-consistent throat radius.** _Addressed (PRs #55–#58)._ The
  imposed `R_MID` is recast as a **finite-self-energy stable equilibrium**
  `R* = (A/2B)^{1/3}` of `E(R) = A/R + B·R²` (EM repulsion vs cohesion);
  the throat caps the EM field so `U_EM/(m c²) = α/2` (finite, no UV
  divergence). The cohesive `B·R²` is **derived** as the throat brane
  tension `4πσR²` (the `R²` power uniquely selected by power-counting),
  with `σ` set by the bulk gravity sector at the **exact** RS fine-tuning
  `λ_crit = √(6|Λ₅|)/κ₅²` (dimensionless factor `√6` from the `Z₂`
  Israel junction + bulk `AdS₅`; the flat / static-throat condition
  `Λ₄ = 0`). The **pair-production threshold** `2 m_e c²` falls out as
  twice the lowest stable throat, forced into a C-conjugate
  throat–antithroat pair by Hopf-charge / antipodal-`Z₂` conservation,
  with the Schwinger critical field `e E_S R_MID = m_e c²` tying the
  scale to the threshold. Consistent with B4, the absolute `R*` still
  rides on one dimensionful coupling — the chain
  *imposed `R_MID` → invariant `ΔR` → finite-self-energy equilibrium*
  recasts and relates the anchor, it does not derive the value.
  Remaining: matching the canonical RS brane to the exact BAM throat
  junction from `S_BAM`, the full instanton/tunneling nucleation rate,
  and the heavier-lepton thresholds (`2 m_μ c²`, `2 m_τ c²`). The same
  instanton/tunneling rate is now doubly motivated: PR #87 shows the
  neutrino's Majorana scale `M_R = m_D·e^{S}` is set by the
  throat↔antithroat bounce action `S ≈ 15–18`, so deriving the
  nucleation rate would simultaneously fix the absolute neutrino mass.
  See `docs/self_consistent_throat_radius_research_plan.md`,
  `docs/cohesive_tension_derivation_research_plan.md`,
  `docs/brane_tension_tuning_research_plan.md`,
  `docs/pair_production_threshold_research_plan.md`, and
  `docs/seesaw_scale_nucleation_compliance_research_plan.md`.
- **Stable moving throats.** _Addressed (`stable_moving_throat_probe`)._
  A boosted throat obeys the relativistic dispersion
  `ω(k)=√(ω₀²+c²k²)`, so `E²−(pc)²=(mc²)²` with the invariant mass equal
  to the static rest eigenvalue `ω(1,0)` to machine precision — `m c²`
  for a moving throat agrees with the static eigenvalue (the throat is a
  particle). It contracts as `R*/γ` with a boost-invariant proper frame
  and stays stable (`d²E/dR²>0` is a rest-frame condition). The closed
  `S³` breaks *global* Lorentz invariance (a preferred frame), but the
  finite-size violation is suppressed by `(R_MID/R_cosmo)² ~ 10⁻⁷⁸` —
  local Lorentz covariance holds. The companion **spin** test
  (`spin_wigner_rotation_probe`) confirms the throat's Hopf-holonomy spin
  (`A_φ = ½ cos χ`, `∮A = π cos χ`) reproduces the relativistic **Wigner
  rotation**: both are the spin-½ `SU(2)` holonomy — the same `½` factor,
  the spinor double cover (`2π → −1`, the Hopf/RP³ structure), and the
  geometric-phase law "rotation = ½ × solid angle"; two non-collinear
  boosts compose (in `SL(2,C)`) to the Wigner `SU(2)` rotation matching
  the closed form. So the boosted throat is a genuine relativistic
  spin-½ particle. The **magnetic moment** completes the spin sector
  (`gyromagnetic_ratio_probe`): `g = 2` follows from the throat's
  Pauli/SU(2) spinor structure (`T = iσ_y`) minimally coupled to the
  Hopf monopole (`A_φ = ½ cos χ`) — `(σ·D)² = D² − eσ·B`, the σ·B term
  carrying the full `σ = 2S` (the factor 2 = the `SU(2)` anticommutator),
  giving `μ = μ_B`; and `g = 2 ⟺` the BMT anomalous precession vanishes
  (spin tracks momentum, the Thomas/Wigner link). The Schwinger anomaly
  `a = (g−2)/2 = α/2π`: the **one-loop** correction is reconstructed
  (`throat_vertex_loop_probe`, PR #62) as the throat dressing its moment
  by one virtual-photon self-exchange — the virtual photon an S³
  Green-function exchange (flat limit `1/q²`), the vertex the throat
  pinch — with the Feynman-parameter integral `∫₀¹ 2z dz = 1` giving
  `F₂(0) = α/2π` (`g = 2.00232…`, vs `a_e = 0.00115965` to ~0.15%).
  PR #62 inherited the `1/(2π)` silently from the tree normalization;
  **`s_bam_loop_measure_probe` (PR #74)** identifies that `1/(2π)`
  explicitly as the **BAM closure-quantum loop measure factor** — the
  same `2π` that underlies `action_base`, `Φ_avail(k) = 2π(k+1) + …`,
  `β_lepton = k_5²·(2π) = 50π`, the Hopf holonomy, the throat dwell,
  and `ε`'s `4β/(2π) = 100`. The structural correspondence: a closed
  cycle of length `L` has density of states `L/(2π)`; for `L = 2π`
  (BAM's S³ great-circle quantum) the loop integration measure is
  `dk/(2π)` — **one closure quantum per loop momentum dimension**, the
  same primitive in QFT and BAM. **Honest scope:** this is the
  structural identification of `1/(2π)` = BAM closure quantum; a fully
  rigorous covariant `(2π)^d` Fourier measure derivation from a
  written-out `S_BAM` path integral on the throat configuration space
  (explicit path integral, gauge fixing, Jacobians) is taken up by PR #115
  as a loop-measure formalism — a throat *is* its closure loop, so the
  configuration space is loop space `LS³` and the measure is
  `Z = Σ_{k odd, c₁∈ℤ, n_part} ∫_{LS³/(Diff S¹ ⋉ U(1)_Hopf ⋉ Z₂)} Dμ[X]
  e^{−S_BAM[X]}` (the sector sum = the closure ledger). It fixes the
  structure — the closure quantum `2π` is the loop holonomy, the odd-k
  lemma is upgraded to the `Z₂` orientation-anomaly condition
  (`e^{ikπ} = −1 ⟹ k odd`), and the PRs #87–#90 bounces are the leading
  saddle — and sets up the `Diff(S¹)` Faddeev–Popov (`bc`-ghost) gauge
  fixing, with the fluctuation operator (the second variation of `S_BAM`,
  the Tangherlini cavity operator) stable (`min ω² ≈ 1.11`). PR #116 then
  closes that analytic core: the divergent bare determinant is regularized
  to a finite, scheme-independent value by two independent standard methods
  that agree — the Gel'fand–Yaglom boundary-value construction gives
  `det(H)/det(H_free) = y(L)/L = 1.574` (log `0.454`) from a single
  initial-value solve with no mode sum (converged to six digits, zero
  interior nodes), and the zeta/heat-kernel method gives `ζ(0) = −1/2`, the
  universal Dirichlet-interval value (finite, no zero mode, no anomaly),
  with the Weyl law `a_{−1/2} ≈ L/√(4π)` confirmed to under a percent. So
  the `S_BAM` one-loop measure factor is *finite and computable*, not merely
  formal. What remains is a closed-form analytic expression (the determinant
  is a definite number, computed numerically) and the absolute normalisation
  of `Z` (which still carries the bulk `κ₅²/Λ₅` anchor); the saddle results,
  normalisation-independent, are unaffected throughout. PR #117 then closes
  the *gauge* half: the reparametrization group `Diff(S¹)` is gauge-fixed
  (worldline/Polyakov) to one Teichmüller modulus — the loop circumference
  `L`, i.e. the Schwinger proper time — and one conformal Killing vector, the
  rigid `U(1)` rotation. The Faddeev–Popov operator is `P = d/dτ` (the vector
  ghost `c` mapped to the einbein variation), with `P†P = −d²/dτ²` on
  periodic fields and kernel `= ` the constants `= ` the one CKV. The ghost
  determinant is the `bc`-ghost path integral `Δ_FP = det'(P) =
  det'(P†P)^{1/2}` (the two coincide in 1D by the `±n` mode pairing): the
  intermediate Laplacian determinant is `det'(P†P) = det'(−d²/dτ²) = L²`
  (`ζ(0) = −1`), so the ghost determinant is its *square root*, `Δ_FP = L` —
  not `L²` (a correction made on review). `Δ_FP = L` is the Jacobian of the
  einbein → proper-length gauge fixing, so the modulus measure is the proper
  circumference `dL`; dividing out the CKV (`Vol U(1) = L`) gives the
  symmetry factor `1/L`, and the worldline measure is `dL/L`, whose `1/L =
  1/(2π)` for the closure loop is precisely PR #74's per-loop factor — so
  PR #74's `1/(2π)` is the CKV (ghost zero-mode) factor of the `Diff(S¹)`
  quotient, independent of the determinant power. And, unlike the 2D string
  (`bc`-ghost `c = −26 ⟹ D = 26`), the 1D worldline carries no
  Weyl/conformal anomaly, so this gauge-fixing is anomaly-free — the only
  nontrivial anomaly being the discrete `Z₂` orientation (odd-k) condition.
  With both halves in hand the one-loop measure `Z = Σ_sectors ∫ (dL/L) ·
  det^{−1/2}_matter · e^{−S}` (the ghost determinant `Δ_FP = L` being the
  proper-length Jacobian, an `L¹` not `L²` power) is finite and computable
  factor by factor; the absolute normalisation (the `κ₅²/Λ₅` anchor) and the
  multi-loop measure remain the standing open pieces. PR #118 audits this
  ghost sector in full to fix the `L`-power unambiguously: it separates the
  four objects `P = ∂_τ`, `P†P = −∂_τ²`, `det'(P)`, `det'(P†P)`; confirms
  `det'(P†P) = L²` and `det'(P) = det'(P†P)^{1/2} = L`; checks the phase via
  the η-invariant (`η(−i∂_τ) = 0` by the symmetric spectrum, so `det'(∂_τ) =
  +L` with no anomalous phase, and in the antiperiodic Möbius sector `η = 0`
  but there is no zero mode and hence no CKV); and proves that the
  conformal Killing vector is divided exactly once. The last point is the
  subtle one: the ghost field space splits orthogonally as `ker(P)` (the
  CKV) ⊕ `ker(P†)` (the Teichmüller modulus) ⊕ the nonzero modes, and the
  Faddeev–Popov determinant `det'(P)` is the *primed* determinant over the
  nonzero modes alone (the SVD of `∂_τ` shows exactly one zero singular
  value, whose right-null vector is the CKV and left-null the modulus). The
  CKV norm therefore enters *only* the gauge-orbit volume `Vol(CKG)`, the
  modulus norm *only* the `dL` measure, and `det'(P)` excludes both — so each
  zero mode is divided once. (An earlier draft divided additionally by the
  zero-mode norms `√L·√L` alongside `1/Vol(CKG)`, double-counting the single
  CKV whose norm is already inside `Vol(CKG)`; corrected.) The conclusion is
  that the FP ghost is first-order, contributing `L¹` through `det'(P)`; the
  `L²` value is reached only by adopting an explicit second-order ghost
  convention, which over-counts by one power of `L`. The measure is
  `Z = Σ_sectors ∫ (dL/L) det^{−1/2}_matter e^{−S}` — the single `1/L` the
  CKV factor `= 1/Vol(CKG)` (PR #74's `1/(2π)` at `L = 2π`), with `det'(P) =
  L` folding into the matter heat kernel — net `dL · L^{−1−d/2}` (`d` the
  matter dimension). PR #119 then supplies the mathematical framework behind
  the one ingredient PR #118 had only asserted — the *phase* of the
  first-order determinant `det'(∂_τ)`. Writing the determinant of the
  self-adjoint `A = −i∂_τ` with a branch choice for the negative
  eigenvalues, the Singer/Atiyah–Patodi–Singer formula gives `det'(A) =
  |det'(A)| · exp[±i(π/2)(ζ_{|A|}(0) − η_A(0))]`, so the phase splits into a
  *local* (heat-kernel / scaling) piece set by `ζ_{|A|}(0)` and a
  *topological* piece set by the η-invariant `η_A(0)`, the intrinsic
  spectral asymmetry. Threading a `U(1)` holonomy `a` (the Hopf/Wilson
  holonomy `∮A = e^{ikχ}`, `a = kχ/2π`) gives, via the Hurwitz zeta `ζ_H(0,a)
  = ½ − a`, the linear law `η_A(0) = 1 − 2a`. The two physical BAM sectors
  sit exactly at its zeros: the orientable (periodic) loop at `a = 0` — where
  the reduced η vanishes once the zero mode (the CKV) is removed — and the
  Möbius (antiperiodic) loop at `a = 1/2`, where there is no zero mode and
  the half-integer spectrum is symmetric. So both are `η = 0`, the
  determinant is real, and the closed forms make it concrete: `det(∂_τ +
  m)_periodic = 2 sinh(mL/2)` whose residue is `det'(∂_τ) = L`, and
  `det(∂_τ + m)_antiperiodic = 2 cosh(mL/2)` giving `det(∂_τ) = 2`. This
  derives PR #118's real `+L` rather than asserting it; a genuine η-phase
  `exp[−i(π/2)(1−2a)]` survives only for an intermediate Hopf holonomy, which
  is the open handle. Finally, PR #120 closes the loop on the *numerics*:
  a high-resolution lattice validation confirms that the discrete
  finite-difference operators in the code reproduce these continuum analytic
  results — the `−∂_τ²` eigenvalues converge to `(2πk/L)²` and the lattice
  log-determinants to `(2 sinh(mL/2))²` (periodic) and `(2 cosh(mL/2))²`
  (antiperiodic) at the second-order rate `O(1/N²)` (error ratio 16 per
  `N×4`), the Tangherlini Gel'fand–Yaglom ratio sits at `1.574370`, and the
  structural quantities (the vanishing η-invariant, the single zero mode, the
  spectrum's `k ↔ N−k` symmetry) hold *exactly* at finite `N` — so the
  software behaves exactly as the analytic derivation. The validation also
  covers the generic twisted-holonomy sectors `a ∈ {1/4, 1/3, 2/3, 3/4}`
  (boundary condition `ψ(τ+L) = e^{2πia}ψ(τ)`, eigenvalues `2πi(n+a)/L`):
  there the magnitude `|det P_a| = 2 sin(πa)` is reproduced *exactly* on the
  lattice at any `N` through the product identity `Π_k 2(1−cos(2π(k+a)/N)) =
  |1−e^{2πia}|² = 4 sin²(πa)`, the twisted eigenvalues and massive
  determinant converge `O(1/N²)`, and `η(a) = 1 − 2a` with — because the
  twisted operator has no zero mode and `ζ(0) = 0` — a phase that is purely
  the η piece, `arg det P_a = (π/2)(1 − 2a)`, so `det P_a = 2 sin(πa)
  e^{i(π/2)(1−2a)}` (reducing to the real `2` at the antiperiodic point
  `a = 1/2`). PR #121 then organises all of this into a **sector-phase
  ledger** with a clean separation theorem: the BAM loop-measure phase
  factorises as a *continuous* η-phase `e^{i(π/2)(1−2a)}` — from the U(1)
  holonomy `a` (the connection / Hopf–Wilson line) — times a *discrete* `Z₂`
  sign `(−1)^k` — from the Möbius / odd-k orientation (the first
  Stiefel–Whitney class `w₁`). These never double-count, for three
  independent reasons: they live in different groups (`U(1)` vs `Z₂`), arise
  from different geometry (connection vs orientability), and — most sharply —
  the continuous phase `θ(a) = (π/2)(1−2a)` is confined to the open right
  half-circle `θ ∈ (−π/2, +π/2)` for `a ∈ (0,1)`, so it *never* reaches `−1`;
  the Möbius `−1` is therefore inaccessible to it and is purely topological.
  At `a = 1/2` the η-phase is exactly `+1`, so the antiperiodic determinant's
  Möbius character is carried entirely by `(−1)^k` — the cleanest witness
  that the continuous and discrete sectors are independent. The measure phase
  is the product `det_full = |det P_a| · e^{i(π/2)(1−2a)} · (−1)^k`, each
  factor counted once. PR #122 closes the arc by *assembling* all of these
  validated pieces into the full factorized sector sum
  `Z = Σ_{k odd, c₁∈ℤ, n_part} (−1)^k ∫₀^∞ (dL/L) det^{−1/2}_matter
  e^{i(π/2)(1−2a)} e^{−S_BAM}` — the closure-ledger sum, the discrete `Z₂`
  orientation sign, the gauge-fixed `dL/L` moduli measure (whose `1/L` is the
  closure-quantum CKV factor), the finite matter and ghost determinants, the
  continuous η-phase, and the leading bounce. Because the `Z₂` sign is a
  sector-constant (the winding parity, independent of the moduli `L` and
  holonomy `a`), it pulls out of the continuous integral, so `Z` factorises
  cleanly into a *discrete* `Z₂`-signed (topological) sum of *continuous*
  η-phased (analytic) moduli integrals — the two never double-counting
  (PR #121). The grading even has teeth: the leading heat-kernel (Weyl)
  coefficient `a_{−1/2} = L/√(4π)` is a bulk quantity independent of the
  boundary condition, so it is identical in the orientable (periodic) and
  Möbius (antiperiodic) sectors and *cancels* in their `Z₂`-graded
  difference — each heat trace diverges as `L/√(4πt)` as `t → 0`, but
  `θ_per − θ_anti ~ e^{−π²/t} → 0` is UV-finite, so the orientation grading
  renders the bulk UV of the sector sum finite. What the assembly does not
  fix is the overall scale: the absolute normalisation (the `κ₅²/Λ₅` anchor),
  the full non-perturbative convergence of the sum, and the multi-loop
  measure remain the standing open pieces. PR #123 then puts this grammar to
  work in the quark sector. A `Z₂`-graded partition sum carries a Witten /
  Atiyah–Patodi–Singer *index* — the graded trace `Tr(−1)^k`, a topological
  invariant whose boundary correction is the η-invariant `ξ(a) = (η+h)/2 =
  1/2 − a` of PRs #119–#121. As the holonomy winds once, one closure
  eigenvalue crosses zero, so the index (the spectral flow `ξ(0⁺) − ξ(1⁻) =
  1`) is an integer. Reading this off for the quark sector is illuminating
  precisely because of *what* it fixes: the quark closure count is `N_q =
  2·n_part = 466`, and the factor of two — the even doubling — is exactly the
  `Z₂`-graded structure, the orientation index pairing and doubling the
  modes. So the APS index pins down the §8-*stable* topological content of
  the quark partition (the doubling, even across all twelve `quark_axioms`
  §8 ablations, and the integer spectral flow) while leaving the bare value
  `n_part` — the continuous, ξ-type residual that drifts `216–255` — exactly
  where the compensator audits (PRs #97/#107) put it: undetermined. The index
  formalises the empirical split: the topology is protected and derived, the
  phenomenological value is not. PR #124 runs the identical audit on the
  lepton sector, and the outcome flips in the most informative way. The
  lepton partition is `N_lepton = 4·k₅² = 100`, and here the feeding integer
  is `k₅ = 5` — the bulk dimension `dim(S³)+2`, a *derived* structural number
  (PR #73), not a fit. The APS machinery is the same (the universal
  spectral-flow integer `1`, the boundary term `ξ(a) = 1/2 − a`), but because
  `k₅` is fixed by geometry there is no §8 ablation that can move it: the
  lepton partition is determined in *both* its structure (the `4k₅²` closure
  form) *and* its value, with no residual. The contrast is the point — quark
  `N_q = 2·n_part` has its structure (the doubling) protected but its value
  `n_part` free and drifting, whereas lepton `N_lepton = 4·k₅²` is fixed
  outright. So the same index, applied to both sectors, isolates exactly
  where the program's single undetermined dimensionless integer lives: the
  leptons are the clean, fully-derived case, and the quark `n_part` is — after
  the APS reduction — the *unique matter-partition residual*, no longer an
  unexplained compensator but the one feeding integer the index cannot fix.
  PR #125 collects the two sector audits into a
  single matter-sector APS ledger and reads off the input budget. The
  pattern is uniform: every matter partition factorises as a derived
  topological piece (a structural factor times the integer spectral flow `1`,
  with `ξ(a) = 1/2 − a` the boundary term) times a single feeding integer,
  and only that feeding integer can be a residual. So the lepton partition
  `4·k₅²` is fully derived (`k₅` the bulk dimension), the quark partition
  `2·n_part` carries the one residual `n_part`, and the neutrino sector
  carries `ε` (derived to order of magnitude, value residual). The combined
  picture, tied to the earlier input-budget audits, is then sharp: one
  dimensionful anchor `G` (the bulk-gravity scale, with `m_e` and `√σ`
  descending from it), four dimensionless residuals — `n_part` (the lone
  matter-*partition* residual), `√σ/m_e ≈ 830` (the irreducible lepton/QCD
  ratio), `ε` (the neutrino compliance value), and `α` (the universal
  coupling) — and the universal flavour puzzle. The APS audit does not remove
  any of these; what it adds is the clean statement that, among the matter
  sectors' closure-partition counts, exactly one is undetermined, and it is
  `n_part`.
  PR #126 then audits whether the factorised Z₂-graded sector sum that carries
  all of this — the sum of PR #122, `Z = Σ_{k odd, c₁, n_part} (−1)^k ∫ (dL/L)
  det^{−1/2}_matter · e^{i(π/2)(1−2a)} · e^{−S_BAM}` — actually *converges*
  non-perturbatively, rather than being a formal expression that diverges. It
  factorises over three independent labels, and each piece is finite. The
  winding sum is not an infinite tower: the odd-`k` lemma and the available
  closure phase `Φ_avail(k) = 2π(k+1) + 50π·max(0, k−3)²` cap it at
  `k ∈ {1, 3, 5}` — the three generations, with `k₅ = 5` the bound — so it is a
  finite three-term sum. The Hopf-charge sum is a convergent Jacobi theta,
  `Σ_{c₁∈ℤ} e^{−A c₁²} = √(π/A) · θ₃ → √(π/A)`, the Gaussian `c₁²` action cost
  making it absolutely convergent. The moduli integral
  `∫₀^∞ (dt/t)[θ_per − θ_anti] e^{−m²t}` is finite at *both* ends: at the UV
  (`t → 0`) the Z₂ grading cancels the boundary-condition-independent Weyl
  divergence, leaving `θ_per − θ_anti ~ e^{−π²/t} → 0`, and at the IR
  (`t → ∞`) the mass gap `e^{−m²t}` kills the tail. The grading is doing real
  work — it is exactly the orientation signs `(−1)^k` that remove the UV
  divergence the individual boundary conditions would carry. So the Z₂-graded
  sector sum is `(finite winding) × (convergent Hopf theta) × (finite moduli
  integral)`: it converges non-perturbatively. What stays open is the absolute
  normalisation (the bulk `κ₅²/Λ₅` anchor) and the multi-loop measure — the
  finiteness is established, the overall scale is not.
  PR #127 then lifts the matter background itself to its parent geometry. The
  Tangherlini fluctuation determinant (PR #116) runs the radial cavity operator
  `V = f(r)[l(l+2)/r² + 3 rs²/r⁴]`, `f = 1 − (rs/r)²` — a reduced, radial
  object — and the question is whether that throat is the boundary of a genuine
  five-dimensional geometry or a 4D ansatz dressed up. The lift makes it
  explicit: the parent is the D=5 Schwarzschild–Tangherlini metric
  `ds² = −f dt² + f⁻¹dr² + r² dΩ₃²` with `f = 1 − (rs/r)^{D−3} = 1 − (rs/r)²`,
  and a self-contained numerical curvature computation (metric → Christoffel →
  Riemann → Ricci/Kretschmann) confirms it is a real vacuum: `R_μν = 0`,
  `Λ = 0`, asymptotically flat, with Kretschmann `K = 72 rs⁴/r⁸` finite on the
  whole cavity — the only true curvature singularity is at `r = 0`, behind the
  throat, while `r = rs = R_MID` is a coordinate (horizon) singularity. The
  throat is that 5D horizon, and its spatial section is the round `S³` — exactly
  the Hopf base `S¹ → S³ → S²` the spin/CPT arc was built on. The two
  coefficients of the PR #116 potential are precisely the D=5 reductions of this
  metric — the centrifugal `l(l+2)` is the S³ Casimir `l(l+D−3)` (`D−3 = 2`),
  the curvature term `3 rs²/r⁴` is `(D−2)/(2r)·f'` (`D−2 = 3`) — so `k₅ = D_bulk
  = 5` (PR #73) is realised as the genuine bulk dimension of the metric, not a
  fitted label. The Hawking period carries the closure quantum: surface gravity
  `κ = f'(rs)/2 = 1/rs`, so `T_H = κ/2π = 1/(2π rs)`. Finally the lift reconciles
  this Ricci-flat bulk with the AdS₅ Randall–Sundrum bulk of PR #57 (the `√6`
  tuning): the Schwarzschild–Tangherlini–AdS₅ metric `f = 1 − rs²/r² + k²r²` is
  Einstein with `R_μν = −4k² g_μν`, `Λ₅ = −6k²` (verified), interpolating the
  Tangherlini neck (`k²r² → 0` near the throat) to the AdS₅/RS asymptote
  (`f → k²r²` far away); on the cavity the AdS correction is `O(10⁻²)` for
  `k·rs ≲ 0.1`, so the pure-Tangherlini cavity is the near-throat limit, good to
  ~1%. What stays open is the exact AdS scale `k` — the unpinned bulk ratio
  `κ₅²/Λ₅` (PR #112) — and the full global brane-localised solution; the
  classical bulk geometry of the throat is established, its absolute scale is
  not.
  PR #128 then makes the throat crossing manifestly smooth. The bulk lift left
  the throat `r = rs` as a *coordinate* (horizon) singularity — the Kretschmann
  scalar is finite there, but in Schwarzschild-type coordinates the metric still
  degenerates (`g_rr = 1/f → ∞`). The horizon-regular charts remove it. In
  Eddington–Finkelstein coordinates, with the tortoise `r* = r + (rs/2)
  ln|(r−rs)/(r+rs)|` and `v = t + r*`, the metric is `ds² = −f dv² + 2 dv dr +
  r² dΩ₃²`: at the throat `g_vv = 0` but `g_vr = 1`, so the determinant
  `det g = −r⁶ sin⁴χ sin²θ` is finite and nonzero, and the Kretschmann scalar
  computed in these coordinates is still `72 rs⁴/r⁸` — the same regular
  geometry, now with a nondegenerate metric. The throat is infinitely far in the
  tortoise (optical) coordinate (`r* → −∞`) but only a finite *proper* distance
  away, `∫dr/√f ≈ √(2 rs (r−rs))` — exactly the ε healing length `√(2 rs ε)`
  (PR #112). The Kruskal–Szekeres extension completes the picture: the surface
  gravity `κ = f'(rs)/2 = 1/rs` gives `κ·rs = 1`, and the Kruskal conformal
  factor `F = −f·e^{−2κr*} = (r+rs)²/r²·e^{−2r/rs}` is finite and nonzero at the
  throat (`F(rs) = 4 e⁻²`) precisely because `κ·rs = 1` makes the
  `(r−rs)^{−κrs}` factor cancel the simple zero of `f`. The product
  `UV = −(1/κ²) e^{2κr*}` vanishes at the throat — the bifurcate Killing horizon
  `U = V = 0` — and the maximal extension has the four regions (exterior I,
  interior II, antipodal exterior III, white hole IV). The deepest point is the
  last: the antipodal map `(U, V, Ω) → (−U, −V, Ω_antipodal)` is an isometry
  that preserves `UV` (hence `r`) and exchanges region I with region III, and
  this is exactly BAM's throat ↔ antithroat identification — the `C = inner/outer
  swap` (PR #63) with `c₁ → −c₁` (PR #58). The maximally-extended 5D Tangherlini
  horizon with its antipodal gluing is the geometric stage of *Bulk Antipodal
  Mechanics* itself: the program's defining antipodal structure is the antipodal
  identification of the throat's Kruskal horizon. What this lift establishes is
  the kinematic stage — the smooth crossing, the finite proper distance, the
  antipodal bifurcation; it does not compute the dynamical throat ↔ antithroat
  nucleation rate (the bounce action, PRs #58/#88), which lives on that stage.
  PR #129 then asks what that null throat does to the waves crossing it — the
  boundary condition the 5D horizon imposes on the matter modes of PR #116. The
  separated wave equation `−d²ψ/dr*² + V_l ψ = ω²ψ` has `V_l = f[l(l+2)/r² +
  3rs²/r⁴] ∝ f → 0` at the throat, so near the horizon the modes are the pure
  null phases `ψ ~ e^{±iωr*}` — the ingoing and outgoing null rays. Three
  boundary conditions compete: the ingoing/absorbing one of a standard
  quasinormal horizon (`ψ ~ e^{−iωr*}`, flux lost into the hole), the reflective
  wall of a hard box (the matter cavity of #116), and the antipodal one dictated
  by the #128 identification `Φ(U,V,Ω) = Φ(−U,−V,Ω_antipodal)`. The antipodal
  postulate settles it, and in a way that is graded by angular parity: the
  scalar harmonics on the horizon `S³` are degree-`l` harmonic polynomials, so
  they carry `Y_l(−x) = (−1)^l Y_l(x)`, and single-valuedness of the field under
  the antipodal map forces the radial function to compensate with the same
  `(−1)^l` across the throat — even-`l` modes meet the throat as a Neumann
  antinode (`ψ'(throat) = 0`), odd-`l` modes as a Dirichlet node (`ψ(throat) =
  0`). Both are *real* conditions, so the Klein–Gordon flux `j ∝ Im(ψ*ψ')`
  through the throat vanishes: the throat is a perfect *unitary mirror*, not a
  sink, in sharp contrast with the ingoing horizon whose flux `j = −ω` carries
  probability into the hole. This is the wave-transport face of the program's
  global CPT and unitarity (PR #64) — what falls toward the throat on one sheet
  re-emerges on the antipodal sheet, nothing destroyed — and the resulting
  exterior cavity has a real, discrete spectrum split by parity into even-`l`
  (Neumann) and odd-`l` (Dirichlet) families, the wave-transport echo of the
  even-`k`/odd-`k` Z₂ structure (PRs #67/#121). What stays open is the full
  quasinormal spectrum (complex `ω`, ringdown) and, again, the dynamical
  nucleation rate; the kinematic transport law across the throat is fixed.
  PR #130 computes that quasinormal spectrum and turns the antipodal-vs-absorbing
  distinction into a sharp spectral fingerprint. On the same cavity
  `−d²ψ/dr*² + V_l ψ = ω²ψ` (shell wall at `R_OUTER`), the throat is given either
  the antipodal real l-parity BC of PR #129 or the absorbing ingoing condition
  `ψ'(throat) = −iω ψ(throat)` of an ordinary horizon; the latter puts `ω` in the
  boundary condition, making it a quadratic eigenvalue problem solved by
  companion linearisation. The two spectra could not be more different. The
  antipodal BC is self-adjoint, so its spectrum is exactly real — `Im(ω) = 0` to
  numerical precision — a tower of *undamped* normal modes, sharp zero-width
  lines of infinite quality factor `Q`, split by parity into even-`l` (Neumann)
  and odd-`l` (Dirichlet) families. The absorbing BC is non-self-adjoint, and its
  frequencies are *complex*, `ω = ω_R − i|ω_I|` with `Im(ω) < 0` — damped
  quasinormal ringdown, the fundamental sitting near `1.89 − 1.24i`, with a finite
  lifetime `τ = 1/|ω_I|` and a quality factor `Q = ω_R/(2|ω_I|) ∼ O(1)` because
  the thin cavity leaks fast into the horizon. The physical reading is the
  payoff: a matter state is a sharp mass — a stable or long-lived particle — only
  if its cavity mode frequency is real, and the absorbing throat gives every mode
  a width, a complex mass, a decaying resonance. Only the antipodal, unitary
  throat yields the real, stable spectrum the BAM matter sectors — the
  lepton/quark bound states — actually have. The undamped-versus-ringdown
  contrast is therefore the spectral face of the program's global CPT and
  unitarity (PR #64): BAM matter is stable precisely because the throat reflects
  antipodally rather than absorbing. What stays open is the idealised
  `r* → −∞` horizon quasinormal tower, the coupling to gravitational radiation,
  and the absolute mode normalisation; the absorbing case is the counterfactual
  that shows what the antipodal postulate buys.
  PR #131 is the capstone of this geometric throat arc, and it is worth stating
  plainly what the arc, taken together, amounts to. The five steps — the cavity
  operator (#116), the 5D bulk lift (#127), the horizon-regular charts (#128),
  the null-throat boundary condition (#129), and the quasinormal spectrum
  (#130) — re-verify together as a mutually consistent set: `f(rs) = 0` with a
  finite Kretschmann `K = 72` at the throat, `T_H = 1/2πrs`, an
  Eddington–Finkelstein determinant `det g = −0.299` that is nondegenerate
  across the throat, a Kruskal factor `F(rs) = 4 e⁻²`, a proper distance
  `√(2 rs ε)` to the throat equal to the ε healing length, an antipodal
  fundamental mode that is real (`ω ≈ 1.19`) and an absorbing one that is
  complex (`ω ≈ 1.89 − 1.16i`). The unifying recognition is that all of this is
  one geometric object seen from several sides: the antipodal identification of
  the 5D Tangherlini horizon. That single primitive is the charge conjugation
  `C` (the inner/outer swap, #63), the throat ↔ antithroat nucleation channel
  (#58), the antipodal map `(U,V,Ω) → (−U,−V,Ω̄)` on the maximal Kruskal
  extension (#128), the l-parity unitary-mirror boundary condition (#129), and
  the selector of the real, stable matter spectrum (#130) — five faces of the
  same gluing. *Bulk Antipodal Mechanics* is, quite literally, the mechanics of
  this one identification on the bulk Tangherlini horizon. The honest ledger is
  equally plain: what the arc *derives* is that the throat's parent is a genuine
  curvature-regular D=5 Tangherlini vacuum (Ricci-flat, `S³` horizon = the Hopf
  base, `k₅ = D_bulk`), that its coordinate singularity is removable, and that —
  *given* the antipodal gluing — the boundary condition is the l-parity unitary
  mirror and the matter spectrum is real and stable rather than a decaying
  ringdown. What the arc *postulates* is the antipodal identification itself,
  BAM's defining axiom; the arc shows that axiom is self-consistent (unitary,
  stable-matter-supporting), not that it is forced by anything more primitive.
  And what stays *open* is unchanged by the synthesis: the exact AdS scale
  `k = κ₅²/Λ₅` (PR #112), the dynamical nucleation rate (PRs #58/#88), the
  global brane-localised solution, and the idealised horizon quasinormal tower.
  PR #132 takes up the first of those open items — the dynamical nucleation
  rate — and connects this geometric arc back to the Majorana bounce arc
  (#87–#90), which had computed the bounce action `S` controlling
  `m_ν = m_D e^{−S}` on the EM/tortoise picture without a regular background to
  stand on. Placed on the horizon-regular geometry, the throat ↔ antithroat
  transition — the `ΔL = 2` Majorana / pair-production channel (#58) — is the
  region I ↔ III crossing of the maximal Kruskal extension (#128), mediated by
  the odd `c₁ → −c₁` instanton (the C-swap #63), with the standard bounce rate
  `Γ ∼ [det(H)/det(H_free)]^{−1/2} e^{−S}`. The geometry supplies three things
  the earlier arc could only posit. First, the Euclidean section is a *smooth
  cigar*: Wick-rotating, the near-horizon metric in the proper radius
  `ρ = √(2 rs(r−rs))` is `ds²_E ≈ dρ² + κ²ρ² dτ²`, the flat plane in polar
  coordinates `(ρ, κτ)`, regular with no conical defect precisely when the
  imaginary-time period is `β = 2π/κ = 2π rs` — the Gibbons–Hawking condition —
  so the nucleation temperature is the Hawking temperature `T_nuc = 1/β =
  1/(2π rs) = T_H` and the period is the closure quantum `2π`. Second, the
  bounce action's logarithm is the horizon's own tortoise divergence: the
  tortoise length of the odd path in to the `ε` healing length is
  `L*(ε) = (rs/2) ln(1/ε) + const` (asymptotic slope `rs/2`, verified to four
  digits), so `S ∝ ln(1/ε)` and the exact-horizon limit `ε → 0` costs infinite
  tortoise length, sending `S → ∞`, `Γ → 0`, `m_ν → 0` — the "rigid throat ⟹
  massless neutrino" of #88 now read off directly from the metric, regulated by
  the finite healing length (#112). Third, the one-loop prefactor is the
  Tangherlini fluctuation determinant of #116, `1.574370` — so the geometric arc
  closes on itself: #116 is the prefactor, #127/#128 the regular stage, and
  #58/#87–#90 the bounce. With the `ΔL = 2` tension window `t ∈ [2π, k₅√(2π)]`
  (#89) and `ε ~ R_c³` (#112), the chain still gives `S ≈ 15–18` and
  `m_ν ~ few meV` to order of magnitude; what this PR adds is the regular stage,
  the smoothness condition, the geometric origin of the `ln(1/ε)`, and the
  prefactor, while the inherited residuals — the exact `ε`, the absolute scale
  `κ₅²/Λ₅`, and hence the precise `S` and `m_ν` — are unchanged and stay open.
  PR #133 takes that recurring `κ₅²/Λ₅` residual head-on, not by pinning it but
  by drawing up its ledger. The absolute bulk scale has surfaced as an open knob
  at every step — the RS tuning (#57) fixed only the dimensionless `√6`, the ε
  healing length (#112) left its absolute normalisation to `κ₅²/Λ₅`, the bulk
  lift (#127) and the nucleation rate (#132) both left the absolute scale open —
  and the ledger asks what, exactly, is open. The 5D content is two dimensionful
  parameters, `κ₅² [L³]` (the 5D Newton constant) and `Λ₅ [L⁻²]` (equivalently
  the AdS inverse radius `k = √(|Λ₅|/6)`), against the geometric lengths
  `R_MID` and `ΔR`. Sorting these honestly gives three categories rather than
  one mystery. First, `ΔR = R_OUTER − R_INNER = 0.52 R_MID` is the *scale
  modulus* — the single dimensionful anchor the B4 theorem (#52) proved is
  required, a proper cosmologically-invariant length (#53) — and it sets the
  unit, so it is units, not a residual; the geometry ratios `ΔR/R_MID = 0.52`,
  `R_OUTER/R_MID = 1.26` are fixed. Second, `√6 = λ_crit κ₅²/√|Λ₅|` is the one
  fixed dimensionless tuning, the Randall–Sundrum flatness condition (#57).
  Third — and this is the whole of the recurring residual — the only remaining
  dimensionless freedom is the AdS scale in throat units, `k·R_MID = R_MID/L_AdS
  = κ₅²/Λ₅` expressed in the unit. It is not pinned, but it is *bounded*: the
  cavity correction to the pure-Tangherlini background is `(k r)²` (#127), so
  `k·R_MID ≲ 0.1` keeps it below about `1.6%` across the cavity, which means
  `R_MID ≲ L_AdS/10` — the throat sits deep in the near-flat region of the AdS
  bulk, and that is exactly why the pure-Tangherlini cavity (#116/#127) was a
  good approximation all along. So the bookkeeping is `{κ₅², Λ₅} → {G, the
  gravity-strength anchor `κ₅²/ΔR³`} + {√6, fixed} + {k·R_MID, open but bounded
  ≲ 0.1}`, with `ΔR` the unit: the "`κ₅²/Λ₅` mystery" is one bounded
  dimensionless number, not a multi-parameter freedom. The ledger bounds and
  isolates the residual; it does not pin it, and it adds no new free parameter —
  it is the same #112 residual, now singular and constrained.
  PR #134 turns the same logarithmic bounce length on a different question — the
  flavor hierarchy — and the result is a clean classification rather than a
  solution. If the bounce action is `S = c·L*(ε) = c·(rs/2) ln(1/ε)` (#88/#132)
  and a tunnelling mass is `m = m_0 e^{−S}`, then the logarithm collapses the
  exponential into a *power law* in the throat penetration depth,
  `m = m_0 ε^{c·rs/2} = m_0 ε^p`. Masses are powers of `ε`, not exponentials of
  a linear quantity, and that single observation sorts the three generations'
  three sectors. The neutrino is the only genuine tunnelling sector — chargeless,
  `k = 0`, the neck not propped open — so `m_ν ∝ ε^p` with `p ≈ 4.8` (#112), and
  the generation healing lengths `ε_n ∝ 1/χ_n` (#79) give the correct normal
  ordering; but the steep power amplifies the modest `χ_n` spread, turning a
  roughly twofold spread in `ε` into `2^{4.8} ≈ 28×` in mass — precisely the
  overshoot #113 had found. So the log-bounce governs the neutrino hierarchy's
  *form and ordering*, with the value residual. The other two sectors are not
  log-bounce at all: the charged leptons are Dirac, their masses set by the
  winding ladder `β·k²` (#71), and the quarks are shell-resolving cavity
  overtones (#77–#80, the `n_part` sector) — and both have irregular `ln m`
  spacings, the signature of the flavor puzzle (#97/#107). The flavor hierarchy
  is therefore a *three-mechanism* structure — bounce, winding, cavity — not a
  single log-bounce phenomenon. What the audit does add is an explanation of
  *why* the flavor values are residual: because `m ∝ ε^p` has
  `∂ln m/∂ln ε = p`, a few-fold ambiguity in the throat depth becomes an
  order-of-magnitude ambiguity in mass, so the irreducibility of the flavor
  values (#108) is a consequence of the exponential mass–action relation rather
  than a separate mystery. The audit does not predict any mass; the neutrino
  overshoot and the charged/quark irregular magnitudes stand.
  PR #135 returns to the antipodal horizon and builds the object its boundary
  data defines — the matter-sector exchange kernel, the program's propagator.
  (The gauge sector already had its exchange kernel, the photon `1/q²` derived
  from the S³ Green function in PRs #42–#44; this is the complementary matter
  kernel.) For each angular channel the kernel is the resolvent of the matter
  cavity operator with the antipodal boundary data of #129,
  `K_l(r,r';ω) = ⟨r|(H_l − ω²)^{−1}|r'⟩`, and three properties follow directly
  from that boundary data. First, because the antipodal operator is
  self-adjoint, the kernel is the mode sum
  `K_l = Σ_n ψ_n(r)ψ_n(r')/(ω_n² − ω²)`, with poles at the real normal-mode
  spectrum of #130 — the propagator is literally a sum over the stable
  exchanged modes, with no decaying contribution. Second, self-adjointness makes
  the kernel symmetric, `K_l(r,r') = K_l(r',r)`: the exchange is reciprocal.
  Third, the boundary data decides unitarity: the antipodal (real) condition
  gives real poles and an undamped, unitary kernel, whereas the absorbing
  horizon would give complex poles and a lossy one — so the antipodal horizon is
  exactly what makes the matter propagator unitary, the two-point face of the
  unitary mirror (#129) and the global CPT and unitarity (#64). And the kernel
  carries the same parity grading that fixed the boundary condition: writing the
  full kernel as `Σ_l K_l(r,r';ω) C_l(Ω·Ω')`, the throat ↔ antithroat exchange
  sends `C_l(Ω·Ω') → (−1)^l C_l(Ω·Ω')`, so each angular channel is graded by the
  antipodal sign `(−1)^l` — even channels symmetric, odd channels antisymmetric
  under the C-swap (#63). What stays open is the same as everywhere else at this
  layer: this is the free, one-loop kernel on the fixed antipodal background —
  the propagator of the S_BAM fluctuation measure — not the interacting
  multi-loop kernel, and it does not fix the absolute normalisation; the
  bulk-scale (#133) and flavor (#134) residuals are untouched.
  PR #136 takes the propagator one order further, to its leading interacting
  correction — the one-loop self-energy `Σ` — and asks the natural question:
  does dressing the free antipodal kernel spoil the stability and unitarity it
  had at tree level? The self-energy enters through the Dyson form
  `G(s) = 1/(s − ω_k² − Σ(s))`, with `s = ω²`, so that `Re Σ` shifts the mass
  and `Im Σ` gives a width; a mode remains a sharp, stable particle exactly when
  `Im Σ` vanishes at its pole. For a cubic self-interaction on the cavity the
  one-loop `Σ` is the two-particle bubble
  `Σ_k(s) = Σ_{n≤m} c_{nm}|g_{knm}|²/(s − (ω_n+ω_m)² + i0⁺)`, with the vertex the
  triple overlap `g_{knm} = ∫ ψ_k ψ_n ψ_m dr*` of the antipodal modes, and the
  optical theorem makes `Im Σ` the two-particle phase space: it is nonzero only
  once `s` reaches a threshold `(ω_n+ω_m)²`. The lowest such threshold is
  `2ω_0`, and the lightest mode sits at `ω_0 < 2ω_0`, so its pole `s = ω_0²`
  lies below `(2ω_0)²` and `Im Σ_0(ω_0²) = 0`: the lightest matter state cannot
  decay — energy conservation forbids it — and stays a sharp, real-pole, stable
  particle through one loop. The real part `Re Σ_0` is a finite mass
  renormalisation: the vertex overlaps fall off with mode index, the mode sum
  converges, and the residual UV piece is the same zeta/heat-kernel
  regularisation that gave the #116 fluctuation determinant. The decisive point
  is unitarity. Above threshold `Im Σ ≤ 0` is a genuine decay width, below it is
  zero — the optical theorem holds — and because the throat is a unitary mirror
  (#129) there is *no* horizon-absorption contribution to `Σ` at all: the only
  width is real multi-particle decay, which the lightest mode is kinematically
  forbidden from. That is the sharp contrast with an absorbing horizon, which
  would hand every mode a width already at tree level (#130). So the one-loop
  self-energy extends the tree-level stable spectrum (#130/#135) intact: BAM
  matter is stable not only as a free spectrum but through its leading
  interaction. What stays open is honest and familiar — the interaction vertex
  is modelled rather than derived from the S_BAM measure, the coupling is an
  input, and higher loops, the absolute normalisation (#133), and the flavor
  residuals (#134) are untouched.
  PR #137 takes up exactly that flagged input — the modelled cubic vertex — and
  draws its ledger, asking how much of it the antipodal structure actually
  fixes. The vertex of three matter modes factorises into an angular integral,
  a radial overlap, and an overall coupling, `V = λ · [∫_{S³} Y_{l1}Y_{l2}Y_{l3}
  dΩ] · [∫ ψ_k ψ_n ψ_m dr*]`, and the first two factors turn out to be derived.
  The angular integral obeys a selection rule with two parts: it vanishes unless
  `l1 + l2 + l3` is even, and unless the SO(4) triangle inequality
  `|l1−l2| ≤ l3 ≤ l1+l2` holds. The even-sum condition is the decisive one — it
  is the antipodal parity itself: under the inversion `x → −x`, which is the
  throat ↔ antithroat C-swap (#63), each harmonic carries `Y_l → (−1)^l Y_l`, so
  the integrand over the inversion-symmetric three-sphere survives only when
  `(−1)^{l1+l2+l3} = +1`. This is the *same* `(−1)^l` Z₂ that fixed the antipodal
  boundary condition (#129), graded the exchange kernel (#135), and sorted the
  flavor sectors (#134); the cubic vertex respects it too, so the one-loop
  self-energy bubble of #136 connects only even-sum mode triples. The radial
  factor is geometric: a definite overlap of the antipodal cavity modes (#116),
  totally symmetric in its three indices and real, so the vertex *shape* is fixed
  by the geometry. What is *not* fixed is the overall coupling `λ` — the
  dimensionless strength #136 set to one — and whether the S_BAM measure
  (#115–#122) generates a cubic term at all. So the ledger reads cleanly: the
  vertex's structure — its selection rule, its geometric shape, its symmetry and
  reality — is BAM-native, while its magnitude is input. The quartic and higher
  vertices, and the bulk-scale (#133) and flavor (#134) residuals, stand.
  PR #138 takes the next vertex — the quartic — and with it answers a question
  the cubic alone could not: whether the interacting vacuum is stable. The
  quartic factorises the same way, into an angular integral, a four-mode radial
  overlap, and a coupling, and its angular selection rule carries the same
  antipodal Z₂: the integral `∫_{S³} Y_{l1}Y_{l2}Y_{l3}Y_{l4} dΩ` vanishes unless
  `l1+l2+l3+l4` is even (and a common SO(4) channel exists), the even-sum
  condition being once more the inversion parity `(−1)^{Σl} = +1` of the C-swap
  (#63). So the `(−1)^l` Z₂ that fixed the boundary condition, graded the
  propagator, and selected the cubic vertices governs the quartic too. The new
  content is the stability audit. A purely cubic potential is unbounded below, so
  the cubic ledger of #137 left the vacuum's stability open; the quartic settles
  it, because the diagonal quartic overlap `g_4 = ∫ ψ_k⁴ dr*` is manifestly
  positive — an integral of a fourth power. The single-mode effective potential
  `V(a) = ½ ω_k² a² + (λ_3 g_3/6) a³ + (λ_4 g_4/24) a⁴` is then bounded below
  whenever its `a⁴` coefficient `λ_4 g_4/24` is positive, which it is for any
  positive coupling, so `V → +∞` at large field for *any* cubic strength: the
  cubic can tilt the minimum but never unbound it, and the vacuum is stable. This
  is not an extra assumption bolted on. A bounded-below action is exactly the
  condition for the path-integral measure `∫ Dμ e^{−S}` to converge, and that
  convergence was already established non-perturbatively for the Z₂-graded sector
  sum (#122); the positive geometric quartic is what realises it. So the
  stability thread of the program closes a loop: the free modes are stable (#130),
  the one-loop self-energy preserves that stability and unitarity (#136), and the
  full interacting vacuum is bounded below (#138) — the same condition the measure
  itself demanded. What remains input is, as always, the magnitudes: the
  couplings `λ_3, λ_4` are not derived from S_BAM (only the sign `λ_4 > 0` is
  forced, by convergence), the quintic and higher vertices are untouched, and the
  bulk-scale (#133) and flavor (#134) residuals stand.
  PR #139 is the capstone of this matter-interaction arc, and it is worth saying
  plainly what the arc, taken whole, amounts to. The six steps — the boundary
  condition (#129), the spectrum (#130), the free propagator (#135), the
  one-loop self-energy (#136), and the cubic and quartic vertices (#137/#138) —
  re-verify together as a mutually consistent set: the harmonic parity is
  `(−1)^l`, the exchange kernel is reciprocal with real poles, the lightest
  self-energy correction has vanishing imaginary part, the quartic overlap is
  positive, and the antipodal fundamental frequency is real where the absorbing
  one is complex. The recognition is that the whole arc is two threads from a
  single postulate. The first thread is the antipodal `Z₂`: the inversion
  `x → −x` of the C-swap (#63) carries `Y_l → (−1)^l Y_l`, and that one parity
  fixes the boundary condition, grades the propagator, and selects which cubic
  and quartic vertices can exist. The second thread is unitarity and stability:
  the antipodal boundary condition is a unitary mirror, and from it follow a real
  stable spectrum, a unitary reciprocal propagator, a self-energy that preserves
  stability with an exactly-stable lightest mode, and a bounded-below interacting
  vacuum — the same boundedness the measure itself required for convergence
  (#122). The two threads are one object seen twice: the real l-parity boundary
  condition is at once the `Z₂` grading and the unitary mirror, both faces of the
  single antipodal identification (#128). The honest ledger is then clean. Given
  the antipodal boundary condition, the selection structure and the unitary,
  stable, bounded interacting theory are *derived*; the antipodal identification
  itself is *postulated* — BAM's axiom, shown here to yield a self-consistent
  interacting theory rather than forced by anything more primitive; the coupling
  magnitudes are *input* (only the sign of the quartic is fixed, by convergence);
  and the open pieces are unchanged by the synthesis — the S_BAM generation of
  the vertices, the higher loops and higher vertices, and the standing bulk-scale
  (#133) and flavor (#134) residuals.
  PR #140 takes up the first of those open pieces — the S_BAM generation of the
  vertices, which every step from #137 to #139 had to flag as modelled rather
  than derived — and closes it structurally. The point is simply that the
  vertices are not separate objects bolted onto the theory: they are the Taylor
  coefficients of the S_BAM action expanded about the throat background. Writing
  `S_BAM[φ_cl + φ] = S_cl + S_2 + S_3 + S_4 + …`, the quadratic piece `S_2` is the
  fluctuation action already met as the #116 determinant and the #135 propagator,
  while `S_3` and `S_4` are the cubic and quartic vertices of #137 and #138 — the
  higher functional derivatives of one geometric action. A free, purely quadratic
  action would have no vertices at all; the geometric, non-quadratic S_BAM
  generates the whole tower. Two of the vertices' properties, moreover, are not
  free but forced. The first is the selection rule. The S_BAM measure carries the
  loop quotient `Diff(S¹) ⋉ U(1) ⋉ Z₂`, whose `Z₂` is the antipodal map — the
  C-swap `x → −x` of #63 — and under it a mode of angular momentum `l` carries the
  harmonic parity, so its amplitude transforms `a_l → (−1)^l a_l`. A vertex of
  several modes therefore picks up `(−1)^{Σl}` and survives the symmetry only when
  `Σl` is even. Because S_BAM is invariant under that antipodal `Z₂`, every vertex
  it generates must have `Σl` even — which is exactly the selection rule #137 and
  #138 found, now read as a Ward identity of the antipodal symmetry rather than an
  assumption. The second forced property is the sign of the quartic. The measure
  exists — it is reflection-positive, it yields the unitary kernel of #135, and it
  converges non-perturbatively (#122) — only if the action is bounded below, and
  that fixes the quartic coupling positive; the geometric overlap `∫ψ⁴ > 0`
  realises it. So the structure of the interaction is generated and constrained by
  the action's symmetry and the measure's consistency. What remains genuinely
  input is narrower than before: not the existence of the vertices, not their
  selection rule, not the quartic sign, but only the coupling magnitudes — the
  numerical higher derivatives of S_BAM — which carry the overall normalisation
  and so inherit the `κ₅²/Λ₅` bulk scale (#133). The exact functional form of
  S_BAM, the higher vertices, and the standing scale and flavor residuals stand.
  PR #141 turns from the matter self-interaction to the gauge sector and joins
  the two. The program had built the photon long ago — the exchange kernel
  `1/q²` read off the S³ Green function (PRs #42–#44) — and the matter sector on
  the antipodal throat across #129–#140; this step asks how the U(1) Hopf gauge
  field couples to that matter at the throat. The coupling is the minimal one,
  `D_μ = ∂_μ − i c₁ A_μ` with `c₁` the Hopf charge, and the interesting content is
  that the antipodal throat is exactly the right place for it to live. The C-swap
  `x → −x` of #63 is a single operation with two effects: on the matter harmonics
  it is the spatial inversion that gives `Y_l → (−1)^l Y_l`, and on the Hopf
  charge it is charge conjugation, `c₁ → −c₁`. So the throat is the
  particle ↔ antiparticle surface, the locus where `C` acts — which is why a
  gauge field that carries charge can couple the matter there at all. Two
  consequences follow. First, the gauge–matter vertex couples a photon to two
  matter legs through the angular triple overlap `∫ Y_{l_γ} Y_{l₁} Y_{l₂}`, which
  is structurally the cubic matter vertex of #137 with one leg now the gauge
  boson, so it inherits the *same* antipodal Ward identity: `l_γ + l₁ + l₂` must
  be even. The one `(−1)^l` Z₂ that threaded the boundary condition, the
  propagator, and the self-interactions now also selects the gauge coupling.
  Second, U(1) charge is conserved at the throat: the antipodal mirror that lets
  no net matter flux through (#129) likewise conserves the charge flux, and the
  C-swap sends outgoing charge back as its conjugate on the antipodal sheet, so
  `Σ c₁ = 0` (#58) and the throat balances particle against antiparticle —
  charge conservation is the gauge face of the unitary mirror. What is *not*
  derived is the one thing that was never going to be: the coupling strength is
  the fine-structure constant `α`, the 137 problem (#105), which the geometry
  organises but does not fix. So the gauge–matter coupling's structure — minimal
  form, the Σl-even vertex, charge conservation, the throat as the C-surface — is
  BAM-native, while its magnitude `α` is the standing universal residual, beside
  the bulk-scale (#133) and flavor (#134) ones.
  PR #142 audits whether that gauge coupling is consistent — and the answer joins
  the gauge sector to the stability thread that has run through the whole matter
  arc. A gauge coupling is consistent only if the matter current it couples to is
  conserved, and the test of that at the throat is direct: the conserved Noether
  current of the global U(1) phase symmetry has, for a stationary cavity mode, a
  time-independent charge density, so conservation comes down to the radial charge
  current `j^r ∝ Im(ψ* ∂_r ψ)`. Because the antipodal cavity modes are real — the
  same self-adjointness that made the throat a unitary mirror (#129) and the
  spectrum stable (#130) — that current is *exactly* zero: no charge flows through
  the throat, the charge sits static, and the particle is a stable charged state.
  The counterfactual is just as sharp. An absorbing horizon would give complex,
  ringing modes whose radial current does not vanish, so charge would leak into
  the hole and current conservation would fail; a charged black-hole-style throat
  is simply not gauge-consistent. Gauge invariance therefore *requires* the
  antipodal throat, exactly as stable matter did. From current conservation the
  rest follows in the textbook way: the Ward–Takahashi identity `q_μ Γ^μ = S⁻¹(p')
  − S⁻¹(p)` ties the gauge vertex of #141 to the matter inverse propagator of
  #135, so the coupling is fixed by the matter dynamics rather than chosen, and
  the vacuum polarisation is transverse, `q_μ Π^μν = 0`, which forbids a photon
  mass and protects the `1/q²` photon (#42–#44). The synthesis is that one
  postulate carries both sectors: the unitary antipodal throat that gave the
  stable spectrum, the unitary propagator, the stable self-energy, and the bounded
  vacuum also gives current conservation, the Ward identity, and the massless
  photon. Gauge invariance is the gauge face of the unitary mirror, not an extra
  assumption — and once again the only thing left input is the coupling strength,
  the fine-structure constant `α` (#105), beside the standing bulk-scale (#133)
  and flavor (#134) residuals.
  PR #143 draws the ledger for that one remaining input, the fine-structure
  constant, exactly as #133 did for the bulk scale — separating what the geometry
  fixes about the electromagnetic coupling from what it does not. A great deal is
  fixed. The charge itself is quantised geometrically: the Hopf number is an
  integer, `|c₁| = 1`, so the unit of charge is topological, not chosen. The `1/2π`
  that famously sits in the one-loop anomaly `a = α/2π` is the closure-quantum loop
  measure of #74, so of that celebrated number the geometry supplies the measure
  and leaves only the prefactor. And the running of `α` — the way the coupling
  flows with scale, through the transverse vacuum polarisation of #142 — is derived
  too. What is *not* derived is the one number the program has always been honest
  about: the value `α ≈ 1/137`, the boundary condition of that running. A
  fit-independent scan against the closure numbers — `2π`, `k₅`, `β_lepton = 50π` —
  finds no clean landing near `137`; the tempting near-misses, `50π − 20` and
  `4·k₅² + 37`, each smuggle in an ad-hoc additive integer of order twenty to
  thirty, which is precisely the reverse-engineering the program rejected for
  `√σ/m_e` in #107 and #108. So `α` is plausibly irreducible in the same sense, and
  the electromagnetic sector contributes exactly one dimensionless residual, the
  value of the coupling, taking its place beside `n_part`, `√σ/m_e`, and `ε` in the
  input budget of #104. The ledger derives the charge quantum, the loop measure,
  the coupling structure, and the running; the value `α` — the 137 problem — stays
  the single open input, beside the bulk-scale (#133) and flavor (#134) residuals.
  The higher-order `a_e` series and the full bulk spinor are the
  related open pieces. See
  `docs/stable_moving_throat_research_plan.md`,
  `docs/spin_wigner_rotation_research_plan.md`,
  `docs/gyromagnetic_ratio_research_plan.md`,
  `docs/throat_vertex_loop_research_plan.md`, and
  `docs/s_bam_loop_measure_research_plan.md`.
- **Charge conjugation from inner/outer swap.** _Addressed
  (`charge_conjugation_swap_probe`)._ C is the inner/outer reflection
  `S: r ↦ 2R_MID − r` — an involution fixing the throat that exchanges
  `R_INNER ↔ R_OUTER`, under which the throat modes are odd (the B3
  antisymmetric extension). The wormhole mouth's induced orientation is
  set by its outward normal `n̂ = ±r̂` (opposite for inner/outer), so the
  swap reverses the mouth orientation and flips the integrated Hopf
  curvature `c₁ → −c₁` (the two orientations `c1_chiphi = −1`,
  `c1_phichi = +1` of `compute_c1`), taking a throat to its antithroat.
  So `C = S`, `C² = id` — charge conjugation as geometry, consistent with
  the antipodal `Z₂` / `T = iσ_y` (B2) and the pair-production antithroat
  (#58). C, P (spatial `S³` reflection), and T (`iσ_y`, B2) then assemble
  (`cpt_assembly_probe`) into the antiunitary **CPT** symmetry on throat
  histories — `q→−, p→+, x→−, s→−, t→−, E→+` with `C²=P²=+1`, `T²=−I` —
  mapping a throat to the antithroat run backwards (the Feynman–Stückelberg
  antiparticle = the pair-production "V" in time, #58). CPT is guaranteed
  by the throat's local Lorentz invariance (#59–#60); the closed `S³`
  breaks *global* Lorentz invariance, suppressing CPT violation by
  `(R_MID/R_cosmo)² ~ 10⁻⁷⁸`. The **explicit CPT operator** on the throat
  Dirac spinor (`cpt_dirac_operator_probe`) is the total-spacetime-
  inversion product `Θ = γ⁰γ¹γ²γ³ = −iγ⁵` (∝ the chiral matrix), built
  from `C = iγ²γ⁰` (the #63 swap), `P = γ⁰`, `T = γ¹γ³K` (the B2 `iσ_y`,
  `T²=−I`): it anticommutes with every `γ^μ` (`j^μ → −j^μ`, the sign
  table above), with matrix `Θ_m² = −I` but antiunitary `Θ² = +I`
  ((CPT)²=+1; the fermionic `−1` is `T²=−I`). The throat 4-spinor itself
  is in turn **derived** from `S_BAM` (`throat_dirac_spinor_probe`): the
  radial operator `H = −d²/dr*² + V` is a perfect square `A†A + E₀`
  (`A = d/dr* + W` the first-order radial Dirac operator,
  `V − E₀ = W² − W′`), its two SUSY-partner sectors (`A†A`, `AA†`,
  isospectral on the nonzero spectrum) are the two wormhole mouths
  (joined by the B3 odd extension, #63), and `4 = 2 (mouths) × 2 (SU(2)
  spin, B2) = Ψ_inner ⊕ Ψ_outer`; parity (`γ⁰`, radial) and the antipodal
  `Z₂` (angular) are disentangled. Remaining: the full closed-form bulk
  spinor with the S³ angular coupling. See
  `docs/charge_conjugation_swap_research_plan.md`,
  `docs/cpt_assembly_research_plan.md`,
  `docs/cpt_dirac_operator_research_plan.md`, and
  `docs/throat_dirac_spinor_research_plan.md`.
- **Even-`k` absence.** _Classified (`even_k_absence_probe`)._ Even-`k`
  modes are absent from the charged-lepton sector by a **spin-statistics
  selection rule**, upgrading the odd-k closure lemma from a
  "choice of sector" to a genuine rule. Each throat pass applies
  `T = iσ_y` (`T² = −I`, B2); the spinor monodromy `T^k` is off-diagonal
  for odd `k` (opposite `Z₂` class — the orientation-reversing closure
  across the non-orientable throat = a spin-½ fermion) and diagonal for
  even `k` (same class — orientation-preserving on the orientable double
  cover `S³` = bosonic). So `k mod 2` is the orientability/spin-statistics
  grading. Charged leptons are spin-½ Dirac fermions (#59–#66), hence the
  odd class; even `k` (bosonic) is excluded — and *not* arithmetically,
  since `Φ_avail(k) ≡ 0 mod 2π` for every integer `k`. The even-`k`
  absence is the spin-statistics face of the same `T² = −I` fermionic
  throat as #60/#61/#65/#66. Remaining: the even-`k` (bosonic) spectrum,
  and why exactly three generations (`k ≤ 5`). See
  `docs/even_k_absence_research_plan.md`.
- **Quark `β` lock.** Listed above. The README correctly flags this as
  a phenomenological compensator under all current ablations.

## Why this matters

If BAM closes — if `ℏ` is geometric, if moving throats are stable, if
the Coulomb law comes out from the connection at finite separation —
then quantum mechanics is a consequence of closed-universe classical
geometry, and Wheeler's geometrodynamic instinct was correct in detail.
That is a strong claim and BAM does not assert it yet.

If BAM partially closes — if the spectra come out cleanly but `ℏ` does
not, or if the static results hold but moving throats are unstable —
then geometry is doing more work than the standard QFT picture credits
it with, and BAM has identified specific geometric channels through
which it does that work. That is itself a result.

If BAM fails on its remaining falsification tests — if `β = 466·π/2`
resists every principled enumeration — then the proposal is wrong in a
way that points to which of the three channels was overcredited, and a
sharper version of the program may still be available. (Two of the
program's most exposed predictions have now passed. The two-throat
Coulomb force test: the force goes as `1/sin²(ψ)` and reduces to the
inverse-square law. And the relativistic-particle tests for a moving
throat: the energy–momentum obeys `E²−(pc)²=(mc²)²` with the invariant
mass equal to the static eigenvalue (`stable_moving_throat_probe`), and
the Hopf-holonomy Berry phase reproduces spin-½ under motion — the
relativistic Wigner rotation (`spin_wigner_rotation_probe`).)

The package is a tool for distinguishing these three outcomes.

## Synthesis: the input budget

Taking the program as a whole, a five-tier epistemic accounting (PR #104)
makes its structure concrete. BAM's *entire dimensionful content* reduces
to **two anchors** — `m_e = ℏc/R_MID` (the QED/lepton scale) and
`√σ ≈ Λ_QCD` (the confinement scale). The B4 scale-modulus theorem
(PR #52) holds that one dimensionful input per sector is mandatory, so
two is the irreducible minimum, and the program sits at it. The
genuinely-open *dimensionless* inputs are localized to two — the neutrino
boundary compliance `ε` (the seesaw/bounce residual, itself bracketed to
`[2π, k_5√(2π)]`) and the quark `n_part = 233` (a compensator for the
flavor puzzle). Beyond these there is a single *universal* open problem,
the flavor puzzle — the quark Yukawa hierarchy, which has RG-invariant
ratios (so it is not a running effect) and is irregular, derivable by no
current theory and so not a BAM-specific failing. Everything else falls
into two productive piles: roughly two dozen *derived* geometric or
topological results (charge quantization, spin-½, `g=2`, the one-loop
Schwinger term, `k_5=5`, `β_lepton=50π`, three generations, the
Bohr-Sommerfeld mass operator, the Cornell confinement form, the
Regge slope, the neutrino Majorana selection rule, PMNS anarchy, generic
CP, …), and about half a dozen *topological predictions* from the
non-orientable structure that span the full testability gamut — matched
(the mesonic `1-+` hybrids), falsifiable (neutrino normal ordering,
`m_ββ ≲ 8 meV`, `Σm_ν ≈ 59–65 meV`), accommodated (the multiquark exotic
zoo), constrained (the light baryonic exotics), findable (the heavy
Möbius baryon), and free (the Möbius glueball tower). In one line: **two
mandatory B4 anchors, a couple of localized open dimensionless residuals,
and the universal flavor puzzle — with the rest derived geometry and a
set of falsifiable non-orientable predictions.** Whether the program
ultimately closes, partially closes, or fails, that is the honest ledger
it has reached.

The same accounting places the fundamental constants (PR #105). Of the
four, `c` is a unit convention and `ℏ` is the closure quantum — the
closure ledger reduces every dimensionless parameter to `2π`-invariants,
and the Compton bridge `ℏ = m_e·R_MID·c` makes `ℏ` geometric once the
single dimensionful anchor is fixed. That anchor is gravitational:
because BAM is GR-foundational, the throat is a gravitational wormhole
whose size — the one mandatory B4 length — is set by the bulk gravity via
the Randall–Sundrum tuning `λ_crit = √(6|Λ₅|)/κ₅²`, so **`G` is the
dimensionful anchor**, the GR scale that the sector anchors `m_e` and
`√σ` descend from and that no gravity-foundational theory can derive from
within. The fine-structure constant is the opposite kind of object: `α`
appears throughout only as a numerical *input* (`A_EM = α·ℏc/2`,
`a = α/2π`), with BAM deriving the charge *unit* (`|c₁| = 1`), the `1/2π`
loop measure, and `α`'s *running* — but never the *value* `1/137`. As in
the Standard Model and every current framework, that value is a free
input, the "137 problem"; so **`α` is a universal dimensionless
residual**, sitting beside the flavor puzzle rather than among the
program's own residuals (`ε`, `n_part`). In short: `ℏ` geometric, `c`
units, `G` the anchor, `α` a universal residual — the two genuinely
irreducible inputs being one gravitational scale and one electromagnetic
coupling, both of which every physical theory must currently take as
given.

Pressing the gravitational anchor one level further (PR #106) asks
whether the two sector scales `m_e` and `√σ` are genuinely independent or
both readouts of the single `G`. They are not independent: both are
brane scales of the one bulk geometry, descending from the same
gravity-tuned tension (`R_MID` and `σ` alike trace to `λ_crit =
√(6|Λ₅|)/κ₅²`). So the *dimensionful*-anchor count collapses from two to
one — the sole fundamental scale is `G`. The catch is that their
dimensionless ratio `√σ/m_e ≈ 830` — the lepton-throat to
QCD-confinement hierarchy — is not derived: it is no clean closure
number (the nearest, `50π·k_5 = 785`, is a 5.4% near-coincidence in the
spirit of `F_13 = 233`). So the reduction is a *repackaging*, not a free
reduction: a dimensionful anchor has been traded for a dimensionless
residual, and the total count of irreducible inputs is unchanged. What it
buys is the cleaner statement of the GR-foundational posture — one
gravitational scale `G` sets the units, and everything else, including
the `m_e/√σ` hierarchy, is a dimensionless number the program either
derives or, for now, carries as a residual alongside `ε`, `n_part`, and
`α`. Deriving that one ratio — fixing the relative normalisation of the
throat-winding and cavity-confinement channels — would reduce BAM to a
single irreducible input.

The first attempt at exactly that derivation (PR #107) is a cautionary
negative result worth recording, because the trap is seductive. The
closure integers offer a near-perfect candidate: with `N_lepton = 100`,
`N_q = 466`, and the gap `ΔN = N_q − N_lepton = 366`, the combination
`N_q + ΔN = 2N_q − N_lepton = 832` sits 0.2% from the observed
`√σ/m_e ≈ 830`. But `2N_q − N_lepton = 4·n_part − 4·k_5²` is built
directly from `n_part`, the quark closure integer that is a
phenomenological compensator — fit to the quark spectrum, drifting
216–255 across the `quark_axioms` §8 ablations. The decisive test is to
propagate that drift: `4·n_part − 100` ranges over `[764, 920]` (±9%)
while the observed ratio is fixed, so `832` is a baseline coincidence
(of the same family as `50π·k_5 = 785` and `F_13 = 233`), not a stable
geometric selection. And no *independent* bulk shell-stress integral
lands near 466 or 832 — the natural ones (`Σω²(n=3..5) ≈ 70`, the
Bohr–Sommerfeld closure sum `Σ(n+1)π ≈ 47`) are `O(10–70)`; the 466
enters only through the v3-fit closure count itself. Recovering the
scale ratio from `n_part` is therefore circular, since `n_part` was fit
to the spectrum that encodes the scales. The lesson is methodological:
a genuine reduction to one input must come from a bulk integral that is
*independent of the spectral fit* and *§8-stable*, not from
recombining the compensator. Until then `√σ/m_e` stays an honest open
residual, and the program's irreducible content remains one
gravitational scale plus a short list of undetermined dimensionless
numbers.

So we ran the search PR #107 called for (PR #108): a quantity built
*only* from BAM's fixed geometry — the structural integers `k_5 = 5` and
`β_lepton = 50π`, and the closure constant `2π` — that selects `~830`
while being §8-stable and free of any tuned factor. The encouraging part
is that §8-stability, the bar that killed the `832` candidate, is
automatic here: a genuinely geometric combination never touches the quark
ablations, so it cannot drift. The discouraging part is everything else.
The best *principled* candidate is `2π·k_5³ = β_lepton·k_5 = 785.4`, which
is `−5.4%` from the target; the next-best principled forms are far worse
(`k_5⁴ = 625`, `−25%`; `e^(2π) = 535`, `−36%`). The only combinations that
reach sub-percent agreement do so through a factor reverse-engineered from
the answer — `π·265 = 832.5`, `(4/3)·k_5⁴ = 833.3`, `k_5⁵/3.77 = 828.9` —
and none of `265`, `4/3`, `3.77` is a fixed BAM quantity, so each is a fit
in disguise, not a derivation. The dimensional-transmutation route fares
no better: `ln(830.3) = 6.72`, while the only clean geometric action
nearby, the closure quantum `2π = 6.28`, is `7%` off, so `830` is not
`e^(action)` for any principled action either. And the genuine Tangherlini
cavity eigenvalue sums are `O(10–350)` — they select nothing near `830`.
The honest reading of this exhausted search is that `√σ/m_e ≈ 830` is not
merely *undetermined* but plausibly *irreducible*, belonging to the same
class as `α` and the electron's anomalous lightness — the universal flavor
puzzle, a pure number the geometry does not fix. BAM therefore does **not**
collapse to a single anchor; it rests at one gravitational scale `G` plus
this one open ratio, `α`, and the flavor puzzle. PR #107's caution is
vindicated: with the `n_part` recycling rejected, the fit-independent
route that would have closed the gap comes up empty, and that empty result
is itself the finding.

### The categorized input budget (PR #150)

The residual-bracket synthesis (`residual_bracket_synthesis_probe`,
PR #150) consolidates the accounting above — #104's five tiers, #105/#106's
constants placement, #107/#108's negative results, #123–#125's APS partition
ledger, #143's α ledger, #133/#148's bulk scale, #113/#149's flavor audits —
into one categorized table, re-verifying a keystone from every category:

| category | item | status | source PRs |
|---|---|---|---|
| **Anchor** (dimensionful) | `G` (→ `ΔR = 0.52·R_MID` unit) | mandatory (B4), relocatable | #52/#53/#57/#106/#133 |
| **Fixed tuning** | `√6` (RS flatness) | derived constant, not a knob | #57 |
| **Universal residual** | `α ≈ 1/137` | structure/measure/running derived; value scan-excluded | #74/#141–#147; #143 |
| **Universal residual** | `√σ/m_e ≈ 830` | one-`G` repackaging derived; value scan-excluded | #106; #107/#108 |
| **Program residual** | `n_part = 233` | doubling topological (APS); value compensator | #97/#123/#125 |
| **Program residual** | `ε` (ν compliance) | order-of-magnitude derived; window `[2π, k₅√(2π)]` | #89/#112 |
| **Bracketed sub-residual** | `k·r_s` | `(0, 0.0064–0.070]` two-sided | #133/#148 |
| **Bracketed sub-residual** | `ε_n` spread | `[1.32, 1.44]`/step, ~0.3%; power laws excluded | #113/#149 |
| **Universal open problem** | flavor puzzle | RG-invariant ⟹ not running; no theory derives it | #97/#107/#108/#134 |
| **No residual** (contrast) | lepton `N = 4k₅² = 100` | structure AND value derived | #124 |

Two features of this table carry the program's epistemic weight. First,
every residual row has **derived structure attached**: the charge quantum,
the `1/2π` measure, and the full one-loop EM sector for `α`; the APS
doubling for `n_part`; the bounce mechanism and ordering for `ε`; two-sided
brackets — derived from the program's own locked spectrum and the
oscillation data — for `k·r_s` and the `ε_n` spread. A residual here is not
a free knob; it is a number boxed by structure the geometry fixes.

Second, the budget is **constant**. The recent arc — #144 (vacuum
polarisation and the running), #145 (`Z₁ = Z₂`), #146 (the charge form
factor), #147 (the `F₁/F₂` capstone), #148 and #149 (the two bracket
audits) — added six probes of derived structure and **zero new inputs**.
The budget today is the #104/#125 budget: one gravitational anchor, two
universal dimensionless residuals shared with every current theory, two
program residuals with derived structure, two bounded sub-residuals, and
the universal flavor puzzle. Whether the program ultimately closes,
partially closes, or fails, it is not failing by knob accumulation — the
ledger is short, categorized, and audited.

**The sensitivity audit (PR #173).** The input budget above is counted by
hand; it can also be *measured*. The dynamical inverse problem — vary the
continuous geometry and read the Jacobian `J_ij = ∂O_i/∂I_j` of the live
observables at the lock — turns the predictive accounting into a
singular-value decomposition. On the 14 currently-reproduced observables (4
quark mass ratios, the 5 CKM magnitudes, `J`, `β`, `γ`, and the 2 lepton
mass ratios) against the free *fitted* knobs (the k₅-derived locks φ_h, χ,
uplift, action, winding excluded as zero-cost), the result is honest and
mixed. The **isolation dimension** rank(J) = 10, with a clean
singular-value gap. The **forced core** — n_obs − rank = **4** — is entirely
CKM combinations: the **CKM unitarity relations** (`V = U₊†U₋` is exactly
unitary, so the 8 CKM observables lie on the 4-parameter unitary manifold,
forcing 8 − 4 = 4 relations). This is the largest observable set the rigid
core forces at zero input cost — a genuine structural prediction, but the
*standard* unitarity, not a BAM-specific numerical relation. The **masses
are fitted** (quark and lepton): no forced mass relation appears. The
**compensator redundancy** — n_inputs − rank = **10** — is dominated by the
mass-preserving diagonal shifts, which is the `n_part`/loose-knob
compensator structure flagged above, now measured: the v4 quark
parametrization is substantially over-complete. And a direct test of the
*"CP at zero parameters"* claim — adding φ_h as an input — leaves the rank
unchanged, so deriving φ_h saves no effective input; the CP economy is a
counting statement, not a Jacobian reduction. The audit confirms the ledger
is not failing by knob accumulation *and* quantifies precisely where the
predictive content is thin — the forced core is real but modest, the masses
are calibrated, and the flavor parametrization carries genuine redundancy
(`sensitivity_jacobian_audit_probe`, PR #173).

### The flavor sector, assembled (PRs #149–#157)

The flavor arc (`flavor_sector_synthesis_probe`, PR #157) converted the
flavor residuals the input budget above carries into an assembled,
falsifiable sector: bracket the residual (#149) → test the mixing/anarchy
hypothesis (#151) → derive the channel-dominant saddle from the bounce
(#152, retiring the one modelling knob) → extract both mixing matrices
(#153 PMNS, #155 CKM) → complete CP in both sectors (#154 Majorana, #156
quark). The card:

| observable | prediction | status | source |
|---|---|---|---|
| mass ordering | normal | derived | #113/#151 |
| `m₁` | ≈ 0.04–0.07 meV | predicted | #151/#152 |
| `Σm_ν` | ≈ 58.8 meV (vs 61.1 uniform-anchor) | falsifiable (~1–2 meV cosmology) | #151 |
| `ε_n` spread | channel dominance (β knob retired) | derived | #149→#152 |
| `sin²θ₁₂/θ₂₃/θ₁₃` | anarchy-natural (62/56/27th pct) | statistical | #153 |
| lepton Dirac CP | generic (`P(\|J\|>0.01) = 61%`) | derived | #153 |
| Majorana phases | generic (`P(\|Φ₂₃\|>π/2) = 69%`) | derived | #154 |
| `m_ββ` | 3.2 meV, 68% [1.5, 5.9]; > 10 meV falsifies | falsifiable | #154 |
| CKM `\|V\|` | all ≤ ×2.0; `V_cb/V_ts` 10% (stiff) | out-of-sample, zero inputs | #155 |
| quark CP | **derived**: `φ_h = π/k₅`; full dataset realized | derived (see addendum) | #156→#161 |

Three features carry the weight. First, the **mechanism map**: the #134
three-mechanism flavor structure is realized at matrix level — the bounce
sector (neutrinos: channel-dominant anarchy through the most compliant
neck), the winding sector (charged leptons: a hierarchy-protected e-row
with exactly one permitted μ–τ rotation), and the shell sector (quarks: Z₂
partition alignment). Large PMNS and small CKM, small θ₁₃ with large θ₂₃ —
each asymmetry traces to derived structure within one geometry.

Second, the **bookkeeping**: eight probes consumed net ONE new input (the
quark CP phase content — the flavor puzzle's CP entry made explicit) and
RETIRED one modelling assumption (the #151 β interpolation, derived in
#152). The modelled-assumption count went *down* while the sector was
assembled.

Third, the **falsifiable targets**: (1) `Σm_ν` 58.8 vs 61.1 meV at ~1–2 meV
cosmology precision; (2) an `m_ββ` detection above ~10 meV falsifies the
ensemble; (3) β = 22° is the acceptance test for the Hopf-connection
`φ_q(k)`; (4) the Jarlskog ceiling must rise to 3.5×10⁻⁵ when the soft
`V_us/V_ub` directions land; (5) `V_cb = 0.038` is stiff at 10%. The
residual locus after the arc: the anarchic draw (statistical), one CP phase
content (input), the soft `V_us/V_ub` direction, and the `O_geom` e-row.

### Flavor phase addendum: the Hopf CP derivation and the full CKM realization (PRs #158–#162)

The card's quark-CP row closed through a correction, a derivation, and a
realization (`flavor_phase_addendum_probe`, PR #162, re-verifies every
keystone in one run).

**The correction (PR #158).** The #156 partition-mixing calibration was an
artifact: with partition mixing on, the charged-current CKM is non-unitary
at ~16% (the u–d near-degeneracy amplifies cross-partition leakage), the
quartet Jarlskog invariants disagree ×1000, the *unitarized* core carries
J ≈ 0 for every `φ_q(k)` form, and the required mixing violates first-row
CKM unitarity ×40. Partition mixing is excluded as the CP origin (and newly
bounded: `ε ≲ 0.004`).

**The relocation and derivation (PRs #158–#160).** The locked
same-partition coupling `−t·e^{−α·dk}·cos(phase·dk)` is the real part of
the Hopf transport factor `e^{iφ·dk}`; the two Z₂ partition classes
traverse the fiber with opposite orientation (the #63 C-swap), giving
`(H±) ∝ e^{±iφ_h·dk}` — exactly unitary CKM, quartet-consistent J. The
scale is derived end-to-end:

| ingredient | value | status |
|---|---|---|
| rate ½ | `A_φ(χ=0)` — the spin-½ factor | derived (connection) |
| sign ± | Z₂ partition orientation | derived (#63 C-swap) |
| winding dk | `max(k, k′)` | locked (the v3 mass calibration) |
| arc 2π/k₅ | the Weyl commutator quantum of the capacity-k₅ fiber | derived (#160 algebra) |
| **φ_h = π/k₅** | 0.6283 | **derived** (#159; all alternative sector counts excluded by data) |

One parameter — and it is not free — yields the full unitarity triangle:
uncalibrated `π/k₅` gives J at 0.97 of target, `(β, γ, α)` within ~2°, and
`sin δ = 0.888` vs the observed 0.887.

**The realization (PRs #160–#161).** The soft `V_us` direction resolved in
stages: single-knob routes excluded exactly (the pinhole breaks `m_s`
−22.5%; the transport rescale self-defeats via level repulsion); the
mass-preserving `SO(3)×SO(3)` family (eigenvector rotations at exactly
fixed eigenvalues) then realizes the **complete nine-observable flavor-CP
dataset at ≤ 1%** — five constrained (`V_us, V_cb, V_ub, β, γ`), four
*predicted and landing* (`V_td ×1.01, V_ts ×1.00, J ×1.00, α = 91.8°,
sin δ = 0.889`) — at the derived phase, with physical down-dominant
re-lock targets tabulated (down-block elements ×1.83/×2.00/×1.11; up-block
×1.29) and the #156/#158 J-ceiling consistency lock verified along the way.

**The bookkeeping.** Across #149–#161 (thirteen probes): **net zero new
inputs** — the one input consumed (#156, the CP phase content) was returned
by the #159 derivation — and **one modelling knob retired** (the #151 β
interpolation, derived away in #152). The quark flavor-CP sector stands as
a consistency statement: locked masses + derived CP phase + the realized
dataset + complete re-lock targets.

**The re-lock, realized and migrated (#163–#164).** The knob-level v3+CP
re-lock is now done. #163 realized the tabulated targets as the **v4
candidate lock** — and found an exact *minimal-law no-go*: the v3
off-diagonal law enforces partition-symmetric transport (`H₊[12] = H₋[12]`)
and the `dk = max` degeneracy (`H[13] = H[23]`), but the targets break both
(a partition split of ratio 1.424 and a minus-block d–b enhancement ×1.996),
while the up block keeps the law *exactly* where data permits (5e-6). The
breaking pattern is the partition asymmetry on the minus block's d-row —
exactly where #155 located the mixing. The realization costs the v3 law +
three new targeted couplings (`η_12^+, η_12^−, η_13^−`) + one retune
(`η_35^−: 5.0 → 5.586`) + diagonal retunes inside the existing law: **+3
parameters buying +5 independent observables (net surplus +2), with the CP
sector at zero parameters** (`φ_h = π/k₅` derived). #164 then migrated this
lock into `geometrodynamics/qcd` *additively* — the v3 lock stays frozen and
bit-reproducible (`φ_h = 0` ⇒ real Hamiltonian, real CKM, no CP), while
`LOCKED_QUARK_PARAMS_V4` and `extract_ckm_matrix()` deliver the masses
(inherited to ~3e-9, the holonomy stripped) and the nine observables (≤ 1%,
unitary) from the library directly. The quark flavor-CP sector is now closed
in code. Remaining: the lepton sector's anarchic draw.
