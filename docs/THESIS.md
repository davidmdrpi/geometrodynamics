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
better-conditioned solver. **Stationarity:** putting `ψ`'s kinetic on an
operator-consistent spectral basis (`u = rψ`, DST) so relaxation and
real-time evolution share the same Laplacian, the relaxed state is a genuine
eigenstate (`‖Hψ − μψ‖/‖ψ‖ ≈ 10⁻⁴`, chemical potential `μ ≈ −1.45`) that
persists under a unitary real-time split-step (profile drift `~4×10⁻⁵`, mass
conserved to machine precision) — a real bound soliton. **Branch scan:** the
soliton is a smooth monotone family in mass (the order field switches on
where `ρ_peak` crosses `ρ_c`, near M ≈ 2.7) and in `q`'s self-gravity `μ`
(`max|q|` 0.42 → 2.62, `Φ(0)` −3.09 → −24.6 across `μ ∈ [0.05, 2]`, residuals
`≤ 10⁻³`, everywhere convergent). **A correction to #179:** #179 reported a
runaway collapse at super-critical `μ` (`|q| → 31`, `Φ(0) → −252`), but that
used a finite-difference (`np.gradient`) Laplacian; the operator-consistent
spectral solver finds no collapse up to `μ = 2`, so the runaway was a
discretization artifact — the genuine large-`μ` limit is the soliton
deepening out of weak-field validity (the strong-field domain for full NR),
not a numerical runaway. **Basin:** the soliton is a robust attractor —
initial widths `w ∈ {1.2, 1.8, 2.6}` and seeds `∈ {10⁻², 10⁻¹}` all flow to
the same state (`max|q|` spread 1.1%, `Φ(0)` spread 0.16%; a tiny seed `10⁻³`
reaches the same attractor more slowly). **Robustness:** the well depth
`Φ(0)` grid-converges (`N = 160 → 240 → 320`: −3.34 → −3.09 → −2.98) while
the pointwise core `max|q|` is more grid-sensitive (~10% per refinement, the
sharp core) — an honest caveat. What survives #179 — the soliton's
existence, two-way back-reaction, and threshold continuity — is confirmed and
hardened; the specific "runaway" claim does not survive as stated. Scope
unchanged: weak-field, semi-dynamical, spherically reduced, effective
constants; the deep-`μ` branch and the strong-field endpoint are for full
numerical relativity (`psi_phi_q_soliton_hardening_probe`, PR #180).

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
