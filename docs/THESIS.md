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
  (explicit path integral, gauge fixing, Jacobians) remains future work.
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
