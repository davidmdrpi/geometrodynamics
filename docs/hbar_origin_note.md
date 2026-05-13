# From closed geometry to a single scale gap

*Where ℏ enters in BAM, after the closure-ledger sequence.*

**Status:** Working paper draft. Companion: `docs/hbar_origin_status.md`
(probe-by-probe ledger), `docs/transport_resistance_research_plan.md`
(derivation thread for the off-diagonal channels), and
`experiments/closure_ledger/runs/` (timestamped JSON+markdown
archives for every numerical result quoted below).

---

## Abstract

The BAM framework (compact S³, non-orientable throat, 5D Tangherlini
bulk) reproduces the three charged-lepton masses to sub-percent by
the eigenvalues of a small, geometry-derived Hamiltonian. Until the
closure-ledger sequence began, six parameters of that Hamiltonian
were phenomenological. We show in five steps that they are not free:
(i) the closure cycle is exactly integer-quantised in units of 2π
for every species; (ii) the cross-species lepton spectrum selects a
single outer radius R\* ≈ 1.262 by a self-consistency loop to
0.008 %; (iii) the pinhole γ, off-diagonal transport, and diagonal
resistance of the locked surrogate are all closure-quantum
invariants (`Σ V_max[1..5]`, `8π`, `7π/100`); (iv) the inner-
boundary regularization ε admits a closed form `ε = resistance/k_5⁴
= 7π/(100·5⁴)` that closes the Compton bridge `ℏ = m_e R_MID c` to
0.04 %; (v) the m_e anchor is the unique remaining external input.
The closure-ledger framework cannot derive m_e itself; the open
question is whether throat dynamics (THESIS.md "self-consistent
throat radius") can lift it.

---

## §1. Introduction

The BAM closure ledger and the locked lepton surrogate together
predict m_μ/m_e ≈ 207 and m_τ/m_e ≈ 3477 to sub-percent. Until
recently, the prediction came at a cost: the diagonal pinhole γ,
the off-diagonal transport, the diagonal resistance, the outer
radius R_OUTER, the closure-quantum integer 100 in the τ uplift,
and the inner-boundary regularization ε were all fit-or-fixed by
hand. The dimensional bridge to ℏ ran through a residual numerical
factor 1.054 = ω(l=1, n=0) of the lowest Tangherlini eigenmode at
R\*. The framework was "dimensional-ratio-complete and dimensional-
scale-incomplete", with several external constants.

This note shows that under a sequence of closure-ledger probes,
the residual external input collapses to a single anchor: the
electron mass m_e. The chain of structural derivations is

```
2π ledger  →  R*  →  (γ, transport, resistance)  →  ε  →  Compton bridge  →  m_e
```

Each arrow is a probe with a quantitative test; each parameter is
identified with a closure-quantum invariant or a Tangherlini-grid
quantity computed on the same geometry. The narrative below tracks
this chain through §§2–5; §6 discusses what remains.

The three structural inputs of BAM are summarised once:

1. **Closed spatial S³** with great-circle action `2π` — the
   antipodal cavity that quantises closure.
2. **Non-orientable throat** with transport map `T = iσ_y`,
   satisfying `T² = −I` — the geometric source of the Z₂ partition
   class and the spinor double cover.
3. **5D Tangherlini bulk** in the radial direction with metric
   factor `f(r) = 1 − r_s² / r²`, discrete radial bound states
   at angular indices `l` and radial labels `n`.

The closure-ledger experiment tracks the total phase accumulated
along a closed worldline through four wired channels (antipodal
closure, Hopf holonomy, throat phase, β-uplift). The earlier
closure-ledger experiment showed the Layer-1 sum closes mod 2π
universally for odd k (the odd-k lemma; `docs/odd_k_closure_lemma.md`).
The new contribution is treating the closure cycle as an *action*
and asking whether it is integer-valued in units of 2π for each
species individually.

---

## §2. Integer 2π closure ledger

**Layer 1 — closed-form channels.** Summing the four wired channels
at the locked baseline `(action_base = 2π, χ = 0, T², β = 50π)`:

```
Φ_total(k)  =  2π·k  +  π·cos(0)  +  π  +  50π · max(0, k − 3)²
            =  2π·(k + 1)  +  50π · max(0, k − 3)²
```

Per species: electron Φ = 4π, muon Φ = 8π, tau Φ = 212π. Dividing
by 2π gives integer counts `(N_e, N_μ, N_τ) = (2, 4, 106)`. The
Hopf holonomy π and the throat phase π *partner* at χ = 0 to form
a single closure quantum 2π; the τ-uplift contributes 100·(2π) from
the integer-winding lock `4β = 100·(2π)`. Layer 1 is closure-cycle
integer-quantised by construction.

**Layer 2 — the radial bulk.** The closed-orbit Bohr-Sommerfeld
action for the (l, n) bound mode satisfies

```
S_full(l, n) / (2π)  →  (n + 1)   as n → ∞
```

— hard-wall Bohr-Sommerfeld with the eigensolver's 0-indexed `n`
mapping to the standard BS quantum number `N = n + 1`. The hard-
wall identification is forced by the throat condition `ψ = T ψ`
combined with `T² = −I`: at any T-fixed point ψ_throat = T ψ_throat
= −ψ_throat, so ψ = 0. Numerically, the eigenfunctions vanish at
both grid endpoints to machine precision; the alternative boundary
conditions (DN, ND, NN, soft+soft) all fail by more than 25 % at
high n, while DD matches to better than 0.001 % at n ≥ 2.

Earlier closure-ledger Layer-2 failures (residual closure spread
of 0.326 rad and 0.577 rad in two candidate couplings) are isolated
to *WKB-at-ground-state* error: the exact action is `(n + 1)·2π` per
bound mode; the WKB approximation at n = 0 underestimates it by
10–25 %. Probes that used WKB ground-state actions saw the
underestimate as a closure-spread residue.

**Angular channel.** The Hopf connection `A = ½ cos(χ) dφ` produces
a holonomy `π·cos(χ)` per fibre loop, verified numerically to
machine precision. Under the spinor double cover (two traversals),
this becomes `2π·cos(χ)` — integer-quantised at exactly the three
polar fibres χ ∈ {0, π/2, π}. At χ = 0 (the canonical fibre, where
the BAM lock sits), the Hopf-throat partnership contributes a
single closure quantum.

**Result.** The closure cycle decomposes per species as

```
N_total(species)  =  k                          [antipodal: k·2π]
                  +  1                          [Hopf-throat partnership at χ = 0]
                  +  100 · [k = 5]              [τ-uplift closure quantum]
                  +  Σ_{(l, n)} (n + 1)         [radial Layer-2]
```

With the cleanest coupling (B2_radial_ladder):

| species | k | coupled mode | N_total |
|---|---:|---|---:|
| electron | 1 | (l=1, n=0) | **3** |
| muon | 3 | (l=1, n=1) | **6** |
| tau | 5 | (l=1, n=2) | **109** |

The τ-uplift quantum 100 reappears as 109 − 3 − 6 = 100 — the same
closure-quantum integer that anchored the lepton sector all along.

---

## §3. Self-consistent R_OUTER

The closure cycle is integer-quantised in *geometric* units. The
canonical identification

```
R_MID · m_e c² / ℏ  =  ω(1, 0)
```

ties the throat radius to the electron Compton scale through the
lowest Tangherlini eigenfrequency. If `ω(1, 0) = 1`, the throat
radius equals the reduced Compton wavelength and the dimensional
bridge closes with no fit parameters. We call this the *Compton
bridge condition*. The question for §3 is whether the geometry
already selects R_OUTER, irrespective of whether the bridge closes.

**Two natural conditions on R_OUTER.** On the canonical R_MID = 1
grid:

| condition | R_OUTER | ω(1, 0) | Σ V_max[0..5] |
|---|---:|---:|---:|
| Compton bridge `ω = 1` | 1.449005 | 1.000000 | 23.6308 |
| γ-lock `Σ V_max = 22.5` | 1.262266 | 1.053694 | 22.5000 |

These differ by 14.79 % and cannot both hold — the Tangherlini
spectrum and the barrier-sum γ are both monotonic in R_OUTER but
pass through their natural targets at different geometries.

**The Compton bridge is physically vetoed.** Re-running the locked
lepton surrogate at the Compton-bridge γ = 23.63 breaks both μ and
τ by ~46 %. Sweeping β over `{30π, ..., 200π}` finds no value
recovering both species at sub-percent (best joint fit: 42 % error).
The Compton geometry is *mathematically definable but incompatible
with the locked lepton spectrum*.

**The γ-lock is selected by cross-species self-consistency.** We
bisect independently on each species' mass error:

| anchor | R* | residual on the other species |
|---|---:|---:|
| muon (zero of err_μ) | **1.262239** | err_τ = +0.16 % |
| tau (zero of err_τ) | **1.262338** | err_μ = −0.16 % |

The two species select the *same* R_OUTER to within
`|ΔR\*|/R\* = 0.008 %`. At the fixed point, `Σ V_max[0..5] = 22.4994`
— the locked γ value to better than 0.01. A single geometric
R_OUTER simultaneously fits both lepton mass ratios; γ is not fit
but is the radial barrier-sum on the same grid.

This is non-trivial: γ is one parameter, and m_μ/m_e and m_τ/m_e
are two independent anchors. The cross-species agreement at 0.008 %
is the geometry telling us *which* R_OUTER it wants.

The fixed point is structurally selected: phase-per-pass perturbations
shift R\* by less than 10⁻⁴ % (decoupled); resistance perturbations
shift R\* by up to 3 % at ±5 %; transport perturbations by up to
7 % at ±1 %. The closure-phase channel is decoupled; the remaining
sensitivities live in the cross-shell mixing parameters
`(transport, resistance)` of the locked surrogate. §4 closes those
channels.

---

## §4. γ, transport, and resistance as closure-quantum invariants

The R_OUTER fixed point of §3 still uses the surrogate with γ, T, ρ
as phenomenological inputs. The pinhole-origin probe established
the structural reading for γ:

```
γ_lepton  ≈  Σ_{l = 1..5} V_max(l)  =  22.0    (canonical Chebyshev grid)
```

— the same radial barrier-sum operator that the QCD residual sector
locks at γ_quark = 22.25. The pinhole is structurally the *barrier-
sum on a Tangherlini grid*, not a free parameter.

The transport / resistance derivation thread
(`docs/transport_resistance_research_plan.md`) closes the remaining
two channels. The origin probe scans six categories of candidates
for each parameter; the two leading readings are:

| parameter | locked value | best closure-quantum reading | %Δ |
|---|---:|---|---:|
| `transport_strength` | 25.1 | `8π = 4·(2π)` | +0.13 % |
| `resistance_scale`   | 0.2179 | `7π / 100` *or* `4·(ω(1,0) − 1)` | +0.94 % / +0.48 % |

The transport identification is decisive — `8π` is the 4th closure
quantum, structurally identical in form to the antipodal `k·2π`,
the Hopf-throat `1·2π`, and the τ-uplift `100·2π`. Resistance has
two within-1 % candidates that the origin probe alone cannot
distinguish on the locked geometry.

**Disambiguation by R_OUTER self-consistency.** Re-bisecting R_OUTER
under each resistance reading, paired with `transport = 8π`:

| reading                          | R*_μ      | cross-species | R*-match to locked | γ at R* |
|----------------------------------|----------:|--------------:|-------------------:|--------:|
| locked baseline (control)        | 1.262239  | 0.0078 %      | (ref)              | 22.499  |
| closure-quantum (`8π`, `7π/100`) | 1.262636  | **0.0021 %**  | **0.031 %**        | 22.508  |
| eigenfrequency (`8π`, `4·(ω−1)`) | 1.258316  | 0.0025 %      | 0.311 %            | 22.417  |

Both readings produce tighter cross-species agreement than the
locked baseline — they are cleaner mathematical objects than the
fitted values. The closure-quantum reading wins decisively on the
discriminating criteria: R\*-match to the locked baseline (0.031 %
vs 0.311 %) and γ-match to canonical 22.5 (0.034 % vs 0.370 %). The
eigenfrequency reading is ruled out — the numerical coincidence
`0.218 ≈ 4·(1.054 − 1)` is not preserved when R_OUTER varies through
the bisection.

The closure-quantum ledger of the locked surrogate (sans the inner
boundary, treated in §5) is then complete:

| parameter             | locked value | structural identification |
|-----------------------|-------------:|---------------------------|
| `action_base`         | 2π           | S³ great-circle action     |
| `transport_strength`  | 25.1         | 8π = 4·(2π)               |
| `resistance_scale`    | 0.2179       | 7π / 100                  |
| `pinhole γ`           | 22.5         | Σ V_max[1..5] ≈ 22.0       |
| `β` (τ-uplift)        | 50π          | locked closure quantum     |
| `4β` (τ-uplift integer) | 100·2π     | τ-uplift quantum          |

With this reading, R_OUTER is structurally selected by the BAM
closure-quantum scaffolding alone — no fitted transport or
resistance constants enter the self-consistency loop.

---

## §5. Inner cutoff and Compton bridge restoration

After §§2–4, the BAM framework predicts mass ratios at sub-percent
and integer closure quanta `(N_e, N_μ, N_τ) = (3, 6, 109)`. The
absolute MeV scale runs through the lowest Tangherlini eigenmode:

```
ℏ · ω(1, 0)  =  1.054 · m_e c²    at the self-selected R_OUTER ≈ 1.262
```

The residual numerical factor 1.054 was the apparent remaining
external input. A direct closed-form search over small
`(k_5, π, integer)` combinations returned a clean negative result
(`factor_1054_search_probe`).

### 1.054 is not a Tangherlini eigenvalue

The follow-up `scale_bridge_regularization_probe` asked the deeper
question: *is 1.054 even a converged eigenvalue?* The 5D Tangherlini
radial problem has a singular inner boundary at the throat
(`r = r_s`, where `f(r) → 0`). In tortoise coordinates the throat
maps to `r* → −∞`; the closure-ledger code regularizes by truncating
the grid at `r = r_s + ε` with `ε = 5×10⁻⁴`.

Sweeping ε from 5×10⁻³ to 10⁻⁴ at the closure-quantum R\*:

| ε       | R*(ε)    | γ at R* | ω(1, 0)  |
|---------|---------:|--------:|---------:|
| 5×10⁻³  | 1.267136 | 22.508  | 1.591606 |
| 2×10⁻³  | 1.264136 | 22.508  | 1.325684 |
| 1×10⁻³  | 1.263137 | 22.508  | 1.175060 |
| 5×10⁻⁴  | 1.262636 | 22.508  | 1.053527 |
| 2×10⁻⁴  | 1.262338 | 22.508  | 0.924840 |
| 1×10⁻⁴  | 1.262236 | 22.508  | 0.845462 |

Two facts emerge. **(R\*, γ) are ε-invariant** to better than 0.4 %
over 1.5 orders of magnitude in ε — the mass-ratio prediction
decouples completely from the regularization. **ω at R\* is NOT
ε-invariant** — it drifts from 1.59 down to 0.85, ruling out any
reading of 1.054 as a Sturm-Liouville eigenvalue of the bare
Tangherlini operator.

The 1.054 factor of the prior framing is therefore a regularization
artefact, not a structural eigenvalue. The proper structural object
is the regularization ε itself.

### The Compton bridge is restorable

At `ε* ≈ 3.51×10⁻⁴` (between 5×10⁻⁴ and 2×10⁻⁴),
`ω(1, 0; R*, ε*) = 1` exactly. The dimensional bridge
`ℏ = m_e R_MID c` would close with *no 1.054 factor* at this
regularization — and the closure-quantum machinery still predicts
the lepton mass ratios at the same precision (self-consistency
holds across the entire ε sweep). The Compton bridge that was
physically vetoed under the canonical ε is therefore *recovered*
at the Compton-bridge ε.

### ε is a closure-quantum invariant

The `inner_boundary_derivation_probe` enumerated 319 candidates over
the closure-quantum scaffolding established in §4 and verified each
by computing ω at the candidate ε. The result:

```
ε  =  resistance / k_5⁴  =  7π / (100 · 5⁴)  =  3.5186 × 10⁻⁴
```

closes the Compton bridge to **0.04 %** — within the tight 0.1 %
tolerance. Every coefficient is a closure-quantum invariant from
§4: `7π/100` is the resistance reading, `k_5 = 5` is the τ closure-
quantum integer, and the exponent 4 is the same '4' that appears
in `transport = 8π = 4·(2π)`.

Both (prefactor 7, exponent 4) are uniquely selected by the bridge.
Within the family `Nπ / (100·k_5^M)`:

| family                     | M=3   | **M=4**        | M=5   |
|----------------------------|------:|---------------:|------:|
| ω at ε = `7π/(100·k_5^M)`  | 1.296 | **1.000403**   | 0.810 |

| family                     | N=6        | **N=7**        | N=8        |
|----------------------------|-----------:|---------------:|-----------:|
| ω at ε = `Nπ/(100·k_5^4)`  | 0.979      | **1.000403**   | 1.020      |

Neighbouring (N, M) miss the bridge by 2 % or more. The closure-
quantum form (N=7, M=4) is uniquely picked out at the same
precision level as the prior closure-quantum identifications
(transport 0.13 %, γ 0.034 %, resistance 0.94 %).

### The full closure-quantum ledger

With ε identified, every geometric parameter of the locked
surrogate is structurally determined:

| parameter         | locked value | structural identification |
|-------------------|-------------:|---------------------------|
| `action_base`     | 2π           | S³ great-circle action     |
| `transport`       | 25.1         | 8π = 4·(2π)               |
| `resistance`      | 0.2179       | 7π / 100                  |
| `pinhole γ`       | 22.5         | Σ V_max[1..5] ≈ 22.0       |
| `β` (τ-uplift)    | 50π          | locked closure quantum     |
| `4β` (τ-uplift integer) | 100·2π | τ-uplift quantum          |
| `R*` (outer radius) | 1.2626     | cross-species fixed point  |
| `ε` (inner cutoff) | 3.51×10⁻⁴   | resistance / k_5⁴          |

The dimensional bridge collapses to the clean Compton form:

```
ℏ  =  m_e · R_MID · c                    (Compton bridge at ε = 7π/62500)
```

BAM is *dimensional-scale-incomplete only modulo m_e*.

---

## §6. Discussion

### What the closure-ledger chain achieved

The closure-ledger sequence reduced the residual external input to
the locked lepton surrogate from six parameters to one (m_e):

- six phenomenological constants at start (γ, transport, resistance,
  phase-per-pass, R_OUTER, ε),
- → two (transport, resistance) after the cross-species
  R_OUTER fixed point (§3),
- → zero closure-quantum constants after the transport/resistance
  derivation (§4),
- → one factor (1.054) at the locked-ε reading (start of §5),
- → m_e alone after the inner-cutoff identification (end of §5).

Every closure-quantum identification in this chain is at the
0.04 %–2 % precision level — comparable to the precision of the
mass-ratio predictions themselves. The chain is internally
consistent: the same closure-quantum integers (`2π`, `8π`, `7π/100`,
`50π`, `100`, `5`, `4`) appear in every step.

### Caveats

The closure-quantum identifications are not exact closed forms in
the strict sense. The pinhole γ vs `Σ V_max[1..5]` gap is ~2 %; the
resistance vs `7π/100` gap is 0.94 %; the inner-cutoff vs
`resistance/k_5⁴` gap is 0.28 % on ε (0.04 % on the Compton bridge).
These gaps are at the precision of the locked surrogate's
ability to predict mass ratios. Whether they represent missing
physics or are irreducible (analogous to numerical-precision limits
in a closed-form identity) is open per parameter.

Additionally, the closure-quantum derivation of ε is a *numerical*
identification — it does not derive the hard-wall regularization
scheme. The physical inner boundary of the Tangherlini radial
problem remains a quantum-throat question.

### What this leaves open

Three structural questions remain, in order of depth:

1. **The 0.28 % residual gap on ε.** Whether the gap admits a small
   correction from an un-read channel of the closure ledger, or is
   irreducible.

2. **The physical inner boundary.** A boundary condition derived
   from throat dynamics — either a quasi-regular Frobenius expansion
   around the singular point, or a finite-thickness throat from
   quantum fluctuations (THESIS.md "self-consistent throat radius")
   — would remove the regularization dependence entirely and
   derive ε from independent physics.

3. **The deeper R_MID self-consistency.** R_MID = 1 by convention.
   Determining R_MID dynamically from a self-consistency condition
   (equilibrium throat radius for the locked mass spectrum) would
   lift the m_e anchor itself.

The closure-ledger framework cannot reach any of these — they
require throat dynamics outside the closure-ledger scope. The
closure-ledger has done what it can: every geometric parameter is
now a closure-quantum invariant, and the residual external input
is reduced to a single mass scale.

---

## Methods and reproducibility

Every numerical result quoted above is generated by a single probe
in `experiments/closure_ledger/`; results are archived as paired
JSON (machine-readable) and markdown (human-readable) files under
`experiments/closure_ledger/runs/<timestamp>_<probe>/`. To
reproduce the full chain:

```bash
# §2 — Integer 2π closure ledger
python -m experiments.closure_ledger.closure_cycle_action_probe
python -m experiments.closure_ledger.closed_orbit_radial_action_probe
python -m experiments.closure_ledger.hard_wall_boundary_verification
python -m experiments.closure_ledger.aharonov_bohm_hopf_fibre_probe

# §3 — Self-consistent R_OUTER
python -m experiments.closure_ledger.tangherlini_omega_m_e_probe
python -m experiments.closure_ledger.tangherlini_electron_scale_bridge
python -m experiments.closure_ledger.compton_bridge_feasibility_probe
python -m experiments.closure_ledger.R_outer_self_consistency_probe

# §4 — γ, transport, resistance as closure-quantum invariants
python -m experiments.closure_ledger.pinhole_origin_probe
python -m experiments.closure_ledger.transport_resistance_origin_probe
python -m experiments.closure_ledger.resistance_disambiguation_probe

# §5 — Inner cutoff and Compton bridge
python -m experiments.closure_ledger.factor_1054_search_probe
python -m experiments.closure_ledger.scale_bridge_regularization_probe
python -m experiments.closure_ledger.inner_boundary_derivation_probe
```

Each probe is self-contained and re-runs end-to-end in seconds to
minutes on a single CPU. The eigensolver is the Chebyshev spectral
method in `geometrodynamics/tangherlini/radial.py`; the locked
lepton surrogate is `geometrodynamics/tangherlini/lepton_spectrum.py`.

## Cross-references

- `docs/hbar_origin_status.md` — eight-probe ℏ-origin summary table
  with archive pointers, updated after every closure step.
- `docs/hbar_origin_research_plan.md` — the research-plan document
  that opened the thread.
- `docs/transport_resistance_research_plan.md` — derivation thread
  for the off-diagonal channels of the locked surrogate (§4).
- `docs/odd_k_closure_lemma.md` — the Layer-1 universality lemma
  that the integer counts refine.
- `docs/quark_axioms.md`, `docs/quark_beta_status.md` — companion
  results for the quark sector (parallel closure-quantum analysis).
- `experiments/closure_ledger/runs/` — per-probe JSON + markdown
  archives, one timestamped directory per run.
