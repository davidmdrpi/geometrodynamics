# Shelled-closure axioms (the hadronic-constituent mass ladder, v3)

This document is the design-spec companion to
`geometrodynamics/qcd/quark_spectrum.py` and
`geometrodynamics/qcd/hadron_spectrum.py`. It mirrors
`docs/lepton_axioms.md` in structure and tone, and it inherits the
methodological rule that separates geometric vocabulary from Standard
Model vocabulary.

---

## §0 Governing principle — emergent quantization from fixed bulk connection distances

Spacetime is continuous and classical. There is no postulated Planck-
scale discreteness. What produces discrete particle spectra is that any
given stable defect topology admits only one family of via-bulk closure
paths, and those paths have geometrically fixed lengths.

Concretely: a closed particle state is a wave configuration that
traverses the outer S³ surface, reaches through the bulk to the inner
S³ surface, traverses the inner surface, and reconnects to the outer
surface — potentially multiple times — before returning to its starting
configuration. For a given topology, the total geodesic length of this
closure is not free. It is determined by the embedding geometry, and
it is the same every time the same topology forms.

That fixed path length makes the resonance condition repeatable:

> same topology → same bulk closure length → same resonance condition
> → same particle mass.

Electrons have the same mass everywhere not because spacetime is
discretized but because the narrow-throat closure topology has exactly
one possible via-bulk path length, and that length is a property of
S³ itself, not of the electron.

Everything else in this document — the matrix, the locked action base,
the integer β-windings, the six closure species — is a consequence of
this principle, not an independent axiom. In particular:

- `action_base = 2π` (lepton sector) is the geodesic length of an S³
  great circle, the minimal via-bulk closure path.
- The integer-winding hard-lock `4β = 200π = 100·(2π)` for the τ is the
  statement that the τ closure traverses that minimal via-bulk path
  exactly 100 times before closing.
- Three lepton generations exist because the narrow-throat topology
  admits exactly three stable via-bulk pass counts (1, 3, 5). Not three
  as a parameter, three as a topological census.

---

## §0.5 Methodological rule (binding)

> The first design spec describes topology and its consequences. Labels
> that denote geometric objects — pass count, partition class, shell
> circumference, coupling path, mouth, neck, chamber — are permitted.
> Labels that denote Standard Model constructs — isospin, weak charge,
> flavor-as-SU(2)-doublet, CKM, W, gluon, color-as-gauge-field — are
> **not** permitted inside the axioms, the Hamiltonian, or any public
> identifier introduced in this document or the accompanying modules.

Analogies to Standard Model structures may be noted in a clearly
separated phenomenological-interpretation section (§9 below), but never
inside the axioms or the Hamiltonian. The lepton module already follows
this rule; the hadronic-constituent module inherits it verbatim.

Module and subpackage filenames (`qcd/`, `quark_spectrum.py`) are
pragmatic conventions the repo has already established and are not
constrained by this rule; the rule applies to spec text and module
docstrings.

---

## §1 The S² visual analogy — why the shelled topology is not "lepton with a bigger parameter"

Picture an S² embedded in R³, with a wave propagating on its outer
surface. Antipodal focusing of that wave can happen two qualitatively
different ways:

**Low-energy focus — the Alice handle.** The wavefront converges to a
point at the antipode. If the focal energy is modest, only a narrow
point-spike reaches through the embedding bulk, pierces the outer
embedding sphere, and reconnects to the inner surface of the embedded
S². What forms is a narrow tube — a non-orientable wormhole handle. No
interior cavity. All of the particle's energy lives on the boundary and
in the throat. **This is the lepton.**

**Higher-energy ring focus — the embedded shell.** If the focal energy
is higher, the wavefront does not converge to a point. The caustic is a
*ring* around the antipode, and that ring has enough amplitude that the
entire ring reaches through the bulk, not just a central spike. Mapped
through the bulk, the ring produces a detached closed surface inside
the throat region — an **embedded shell**. The mouth on the outer S² is
still a point-like puncture. But inside the throat, there is now a full
2-surface with its own geometry, supporting its own harmonics, which
contain energy that is invisible from outside.

Observed from the outside, both defects look pointlike at the mouth.
The difference is entirely interior. **Point mouth outside, embedded
shell inside — that is the hadronic constituent.**

The generalization to S³ is direct: replace "point on outer S²" with
"point on outer S³", replace the 2-shell with a 3-shell (the bulk
chamber), and keep the non-orientable identification that the lepton
case already has. The resulting defect has a pointlike mouth consistent
with deep-inelastic scattering observations and a bulk chamber
consistent with most of the composite mass sitting somewhere that is
not the throat.

**The shell is the single most important piece of geometric intuition
in this sector.** It is what makes this defect distinct from a lepton
defect, what supports the two partition classes of §2, what makes the
interior harmonics identical to the Möbius modes already validated in
`qcd/spectrum.py`, and what releases its energy back into the universe
when the defect is forced to collapse.

---

## §2 Applying the governing principle to each sector

### Lepton sector — the minimal closure census

Two geometric length scales: throat traversal length and outer-surface
arc between successive throat encounters. The full via-bulk closure is
therefore a pass count. Non-orientable closure selects odd pass counts,
`k ∈ {1, 3, 5, …}`; only the first three are stable against self-
intersection cost. Three stable pass counts, three fixed closure
lengths, three generations. Already locked in the repo.

### Shelled sector — the shelled closure census

Inherits the pass-count structure `k ∈ {1, 3, 5}` and introduces two
new fixed geometric lengths from the embedded shell:

1. **Shell circumference `L_shell`** — a closed geodesic on the
   embedded shell's surface. Sets the integer-wavelength condition for
   the trapped interior harmonics.
2. **Neck-to-shell coupling path `L_coupling`** — the geometric distance
   from the mouth through the narrow neck to the shell surface. Gates
   energy exchange between exterior and interior.

Both fixed by topology. Neither free.

The Hopf Chern flux of a shelled defect can distribute itself between
the exterior point mouth and the interior shell in two topologically
distinct stable ways, `p = +` and `p = −`. Which of `L_shell` /
`L_coupling` dominates the total closure phase differs between the two
partitions, and the relative contribution is `k`-dependent. That
`k`-dependence is the geometric origin of the observed sign-pattern
inversion across pass counts: the `p`-ordering at the lowest pass count
differs from the `p`-ordering at higher pass counts.

**Six stable closures = 3 pass counts × 2 partition classes.** Color
(already in `qcd/color.py`) enters as a 3-fold degeneracy on each
closure and does not contribute to the mass eigenvalues.

---

## §3 Hamiltonian — structure

Basis: `|k, p⟩` with `k ∈ {1, 3, 5}`, `p ∈ {+, −}`. Canonical ordering:

    (1,+), (1,−), (3,+), (3,−), (5,+), (5,−)

Diagonal:

    H_{(k,p),(k,p)} = action_base
                    + r_q · k²
                    + pinhole · 𝟙[k ∈ {3, 5}]
                    + β · max(0, k − 3)²
                    + γ_q · σ(p) · u_q(k)

with `σ(+) = +1`, `σ(−) = −1`, and minimal ansatz `u_q(k) = k − 2`.

Off-diagonal cases:

1. Same partition, different pass count:

        H_{(k,p),(k',p)} = −transport · exp(−α_eff · Δk) · cos(phase · Δk),
        Δk = max(k, k')  (inherited `winding_mode = max` from lepton sector)

2. Same pass count, different partition — **partition-mixing amplitude**:

        H_{(k,+),(k,−)} = −w_q · exp(i · φ_q(k))

   No lepton analog. The geometric phase `φ_q(k)` comes from the Hopf
   connection (`geometrodynamics.hopf.connection`, `.spinor`). It is
   not a free complex phase.

3. Different pass count AND different partition: zero in the minimal
   v3 model. Turned on only if demonstrated necessary by fit residuals.

---

## §3.5 Spectrum-zero semantics and physical-spectrum extraction (r2)

Raw eigenvalues of the 6×6 Hamiltonian are closure-cost path lengths
relative to an arbitrary zero, not physical masses. To extract masses
we need two pieces of machinery that are absent in the bare
`np.linalg.eigvalsh` result.

### Spectrum zero

The physically meaningful zero is `action_base` itself. This is the
topological minimum closure cost for any shelled defect (§3): every
closed wave configuration accumulates at least this much path length
before it can close at all. The mass of a physical species is the
*excess* closure cost over this topological minimum:

    mass = eigenvalue − action_base

The MeV anchor then converts path-length excess into physical mass,
mirroring how the lepton sector anchors against `m_e`.

### Positivity constraint and the physical parameter regime

A parameter point is only physically meaningful if the anchor species
(and ideally all six species) sit strictly above `action_base` after
the shift. Two mechanisms can push an eigenvalue below `action_base`:

1. **Direct depression from the γ_q term.** At `k=1, p=+`,
   `σ(p)·u_q(k) = +1·(−1) = −1`, so the diagonal entry at that basis
   state is `action_base + r_q·1 − γ_q`. If `γ_q > r_q`, this diagonal
   entry lies *below* `action_base` before any mixing. Therefore the
   scan must enforce `γ_q < r_q · k²` at `k=1`, i.e. `γ_q < r_q`, to
   keep the unmixed (1,+) diagonal above the topological minimum.

2. **Level repulsion from transport and partition mixing.** Large
   off-diagonals can pull the lowest eigenvalue below the lowest
   diagonal entry via level repulsion. This is a physical effect: two
   nearby levels repel, and the lower one descends. The calibration
   must verify that after full mixing the anchor species remains above
   `action_base`.

The calibration scripts handle this by **rejecting parameter points**
where the anchor species has non-positive mass after the shift. They
treat such points as outside the physical regime rather than as fit
failures. Published predictions in §8 must come from a parameter point
strictly inside the physical regime; a "best fit" that relies on a
negative anchor is meaningless.

### Adiabatic species identification

Raw eigenvalue order is not species order under mixing. The module
identifies each species by adiabatic continuation from the fully-
unmixed reference state (all off-diagonals zero), where each
eigenvector is exactly a basis vector `|k,p⟩` and species
identification via `BASIS_TO_SPECIES` is unambiguous. All mixing knobs
(transport, γ_q, partition_mixing) are then ramped together from zero
to target, and each eigenvector is tagged at every step by maximum
overlap with its previous tagged ancestor.

This means the reported `m_u` is the mass of whichever final
eigenvector adiabatically descends from the initial `|1,+⟩` basis
vector — *not* the lowest eigenvalue of the mixed spectrum. Under
strong mixing these can differ (level crossings), and the tracking-
based tag is the physically meaningful one.

### Open question: the γ_q / r_q tension

The tension between the `γ_q · σ(p) · u_q(k)` depression term and the
positivity constraint at `k=1, p=+` is a real structural question the
calibration must resolve. Three possible outcomes:

1. **The constraint `γ_q < r_q` is natural.** The scan finds its best
   fit in the small-γ regime and the constraint is never binding.
   Best outcome; confirms the ansatz is compatible with positivity.

2. **The fit wants `γ_q > r_q`.** In that regime the (1,+) diagonal
   sits below `action_base`, and the species reported as `u` has
   negative path-length excess. Two ways to interpret:
   - The minimal ansatz `u_q(k) = k − 2` is wrong; try
     `u_q(k) = k · (k−2)` which at k=1 gives `−1 · 1 = −1` — same sign
     depression, same problem. Or `|k − 2|` which removes the
     depression entirely but breaks the generation-1 inversion.
   - Or: reinterpret `action_base` not as the absolute floor but as
     a *typical* closure cost, and accept that some physical species
     sit below it. In that case the spectrum zero should shift to
     `min(eigenvalues)` and the MeV-anchor formula becomes
     `mass = scale · (eigenvalue − min_eigenvalue)` with no claim of
     topological significance for the zero.

3. **The tension is a real physical result.** The `γ_q` term captures
   a sub-topological contribution that the minimal v3 spec didn't
   anticipate, and the calibration's best-fit value of `γ_q` is
   itself an integer multiple of some smaller topological quantity
   (e.g., `γ_q = π/something`). In that case the spec gets a new
   locked axiom and the apparent tension resolves.

Which of these is right is an empirical question that cannot be
settled in the spec. It must come out of running the calibration and
examining how the best fit behaves. The code drop is structured so
that all three outcomes can be distinguished by the scan output;
interpretation is left to the human reviewer of the scan log.

---

## §4 Parameter locking — everything is a path length

Five-step locking discipline, mirroring the lepton playbook:

1. **Anchor** one lightest closure mass to set the MeV scale (default
   anchor: the `(1,+)` closure, identified with the lightest species
   in `OBSERVED_MASSES_MEV`).
2. **Lock** `action_base_q` to a topological value (`π`, `2π`, or a
   small integer multiple). The calibration scan reports which.
3. **Lock** `β_q` to an integer-winding value: `β = N · π / 2` such
   that `4β = N · (2π)`, exactly as the lepton sector locked `N = 100`
   for the τ. The calibration scan reports `N`.
4. **Lock** `γ_q` and the shape of `u_q(k)` to whatever minimal
   topological form reproduces the observed sign-pattern across pass
   counts. If `u_q(k) = k − 2` suffices, it is a new discovered
   invariant.
5. **Optimize only** the residual continuous knobs (phase, transport,
   pinhole, resistance, partition_mixing, gamma_q) against the
   remaining five masses with one anchor consumed.

The success metric is not "did the fit come out good." The success
metric is **how few remaining continuous parameters were needed after
the topological locks**. If, like the lepton sector, four or fewer
residual continuous knobs fit five residual masses at sub-percent
accuracy, the model is working as the framework predicts.

---

## §5 Reuse map — what comes from where

| Need | Module | Role |
|---|---|---|
| SU(3) color algebra | `qcd/color.py` | degeneracy on each closure |
| Cornell flux tubes | `qcd/bridge.py` | composite-binding |
| Hadronic network | `qcd/network.py` | Y-junction and meson-tube backbone |
| Composite topologies | `qcd/topology.py` | meson / baryon / glueball / hybrid |
| Möbius shell harmonics | `qcd/spectrum.py` | interior shell modes of §1 |
| String-breaking | `qcd/solver.py` + `qcd/bridge.py` | decay channels |
| Hopf connection / transport | `hopf/connection.py`, `hopf/spinor.py` | partition-mixing phase |

New in v3:

| Module | Role |
|---|---|
| `qcd/quark_spectrum.py` | shelled-closure census (this document) |
| `qcd/hadron_spectrum.py` | composition: sum closure energies + bridge + shell |
| `scripts/calibrate_quark_ratios.py` | coarse grid on residual knobs |
| `scripts/sweep_quark_beta.py` | integer-winding β-lock hunt |
| `scripts/map_basin_quark_uplift.py` | basin-width sanity check |
| `scripts/lock_quark_beta_probe.py` | final hard-lock optimizer |

---

## §6 Success criteria

For `qcd/quark_spectrum.py`:

1. **Emergent quantization is visible**, not imposed. Every Hamiltonian
   entry traces to a concrete closure path length. No entry is a free
   parameter at scan time that cannot be rewritten as an integer
   multiple of `2π`, a topological degeneracy factor, or a Hopf
   transport phase.
2. **Reproduces the six observed masses** in the hadronic-constituent
   sector to within 1 % for 4 of 6 species after one mass anchor and
   all topological locks.
3. **Exhibits at least one new integer-winding discovery** parallel to
   the `100·(2π)` τ lock — a genuinely predicted integer, not a fit.
4. **Lepton limit.** With `γ_q = 0`, `w_q = 0`, and `action_base = 2π`,
   the 6×6 Hamiltonian collapses to two decoupled 3×3 blocks, each
   structurally identical to the lepton Hamiltonian. Eigenvalues pair
   up, and the distinct-eigenvalue ratios match the lepton mass ratios
   to `1 × 10⁻⁶`.
5. **Color-independence.** Eigenvalues are independent of color label
   (structurally: `build_quark_hamiltonian` takes no color argument).

For `qcd/hadron_spectrum.py`:

6. **Lightest-pseudoscalar universality.** The ~130 MeV lightest
   pseudoscalar mass is reproduced independently of the specific pair
   drawn from the first-pass-count sector. Species independence of the
   bridge-plus-shell contribution is the sharpest available test that
   the shell-as-chamber picture of §1 is quantitatively correct.
7. **Lightest three-defect composite at ≈ 938 MeV** with defect-matrix
   contribution of only ≈ 10 MeV, with the remaining ≈ 928 MeV coming
   from the Y-junction bridge configuration already in `qcd/network.py`.

---

## §7 Non-goals for v3

This spec is **not** attempting any of the following, and attempts to
add them should be rejected until the six criteria of §6 are met:

- Deriving the QCD Lagrangian from the topology.
- Reproducing the full hadron zoo (excited states, exotic states).
- Predicting mixing structures between the three pass-count sectors
  beyond what the minimal same-partition off-diagonal already provides.
- Addressing the strong-CP problem.
- Explaining weak-sector phenomenology.

These are all reasonable future work, not preconditions for v3 shipping.

---

## §8 Calibration log — to be populated

Once the calibration pipeline runs, record here:

- The locked `action_base` (π, 2π, or numeric).
- The locked integer winding `N`.
- The basin width around `N` within 2× of the best residual.
- The locked `QuarkParams` as a code snippet.
- The six predicted masses and their relative errors.

This section is the direct analog of the calibration-results section in
`docs/lepton_axioms.md`; it is what the community can cite when
discussing the discovered invariants.

---

## §9 Phenomenological interpretation (post-topology, separated by rule)

**This section is separated from the axioms by the methodological rule
of §0.5.** Content here is commentary on what the topology corresponds
to in Standard Model language, explicitly *after* the topology is
developed and *never* fed back into the Hamiltonian or the axioms.

- The two partition classes `p = +` and `p = −` map onto the
  up-type / down-type distinction in each generation.
- The `k`-dependence of `u_q(k)` that produces the partition-ordering
  inversion between the first pass count and the higher ones matches
  the observed mass-ordering inversion: `m_u < m_d` but `m_c > m_s`
  and `m_t > m_b`.
- The partition-mixing amplitude `w_q` at fixed `k` is the natural
  candidate for the weak-sector charged-current transition amplitude
  at tree level.
- Off-diagonal entries of type (3) — different pass count AND different
  partition — would correspond to cross-generation partition-flipping
  transitions, i.e. the off-diagonal entries of what the Standard
  Model calls the CKM matrix. Keeping them zero in the minimal v3
  model is the topological statement that such transitions are
  suppressed relative to same-pass-count mixings.

None of the above is an axiom. It is post-hoc interpretation. If the
calibration reveals that the natural mapping differs from what this
section suggests, the correct response is to update this section, not
the axioms.
