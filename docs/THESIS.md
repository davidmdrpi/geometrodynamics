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
| Lepton mass ladder | locked spectral model; `β` multiplier still needs analytic origin |
| Quark mass ladder | 1.6% fitted ladder; residual sector geometrized; `β` phenomenological |
| Coulomb radial response | verified by Tangherlini/Maxwell BVP |
| Coulomb force at finite separation | falsification test, not yet closed |
| Black-hole interior / entropy | regular metric derived; entropy currently a consistency check |
| Compton tree amplitude (Klein-Nishina) | reproduced exactly via closed-form F² (PRs #25–#35); F² = K²·Q derived from a single C×S³ master functional (scaffold closure, B5′ closed) |
| Full QFT / Born rule / `ℏ` origin | tree level closed (Compton, BW, annihilation, Bhabha, Møller, PRs #36–#46); loops open; B4 audited — `ℏ` needs exactly one dimensionful anchor (scale-free machinery), relocatable to the invariant bulk separation `ΔR` |

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
threshold, a self-consistent throat persists. A natural candidate for
the lepton-sector bifurcation threshold is `2m_e c²`.

The refinement that distinguishes the present program from "particles
as static defects" is that the throats produced this way are not
required to remain at rest. A particle is a **moving topological
boundary condition** — two mouth positions `X₁(t), X₂(t)` on `S³`, a
bulk throat length `L_throat(t)`, and a time-dependent transport map
`T(t): T_{X₁}S³ → T_{X₂}S³` between tangent frames. This is a stronger
claim than "antipodal focus = particle." It says the focus is the
*trigger*; the particle is the persistent topological response to that
trigger.

## What success looks like — falsification tests

The next phase of BAM is organized around demonstrations, not parameter
fitting. Each is a test the existing framework can be put to that admits
a clean pass-or-fail.

**Odd-`k` closure from topology.** _Closed_ — see `docs/odd_k_closure_lemma.md`.
The lemma states (i) even `k` and odd `k` both admit valid closure
boundary conditions on the throat — even `k` is orientation-preserving
closure on the doubled cover, odd `k` is orientation-reversing closure
across the non-orientable throat. BAM's lepton sector chooses the
orientation-reversing branch, so the physical particle states sit at
`k ∈ {1, 3, 5, …}`. (ii) Under the locked baseline (`action_base = 2π`,
`χ = 0`, `T²` convention, `4β/(2π) = 100 ∈ ℤ`) the Layer-1 ledger sum
is identically zero mod 2π for every integer `k`; the choice of sector
selects which `k` enter the spectrum. The empirical "all three lepton
ledgers close to 0 mod 2π" result of `experiments/closure_ledger` is
the direct verification.

**Moving-mouth Berry phase.** Drag a single throat mouth around a
closed loop in `S³` in the existing solver and verify that the
accumulated geometric phase reproduces `2π → −1`, `4π → +1`. Bridges
the static-spin-½ result (already verified) to the dynamic-particle
picture. Operationalizes the moving-boundary refinement and gives the
project a clean, single-experiment statement of "spin from geometry."

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

**Quark `β = N · π/2` with `N = 466` derived or principled-bounded.**
_Principled-bounded; full derivation open._ See
`docs/quark_beta_status.md` for the closing summary. A five-probe
sequence in `experiments/closure_ledger/` (origin → boundary →
decomposition → audit → sub-block stability) localized the
irreducible structural piece: across all 12 logged §8 ablations the
only preserved invariant is `N_q ≡ 0 (mod 2)`. The structural reading
is `N_q = 2 · n_part` with the factor of 2 topological (the Z₂
partition multiplicity from the v3 Hamiltonian basis `{(k, ±)}`) and
`n_part` the phenomenological compensator. The clean negative
result for the listed enumerations — torus-knot crossings, SU(3)
representation dimensions, S³/S² harmonic counts, Tangherlini barrier
sums — is recorded by the origin probe. Deriving `n_part = 233` from
first principles is open and is unlikely to be the next-most-tractable
work in this framework.

**Three structural quark axioms reduced to one partition principle.**
The quark Hamiltonian currently rests on four shell-index axioms
expressible in `k_5 = 5` only: `(1 − 1/k_5², k_5, (k_5 − 1)·k_5, 0)`.
The first, `ε = 1 − 1/k_5²`, is recognizably an inverse-square shell
asymmetry. If `χ` and `η` can be derived from the same `Z₂` partition
structure rather than just expressed in `k_5`, the residual sector
goes from "geometric expressions" to "consequences of one partition
principle" — a much stronger claim with no new physics required.

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
- **Self-consistent throat radius.** `R_MID` is currently imposed. The
  deeper version determines `R_MID` dynamically as the equilibrium
  throat radius for a given excitation amplitude, with the
  pair-production threshold falling out as the lowest stable
  configuration.
- **Stable moving throats.** A boosted throat solution must remain
  self-consistent. The "is the throat actually a particle" test is
  whether `m c²` for a moving throat agrees with the static eigenvalue.
- **Charge conjugation from inner/outer swap.** Promote the C-symmetry
  from a postulate to a geometric statement — that swapping
  `r < R_MID ↔ r > R_MID` in the throat eigenmodes flips the sign of
  the integrated Hopf curvature.
- **Even-`k` absence.** Listed above as the highest-leverage near-term
  result.
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

If BAM fails on its remaining falsification tests — if Berry phases do
not reproduce spin-½ under motion, if `β = 466·π/2` resists every
principled enumeration — then the proposal is wrong in a way that
points to which of the three channels was overcredited, and a sharper
version of the program may still be available. (The two-throat
Coulomb force test, once the program's most exposed prediction, has
now passed: the force goes as `1/sin²(ψ)` and reduces to the
inverse-square law.)

The package is a tool for distinguishing these three outcomes.
