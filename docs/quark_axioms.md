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

## §8 Calibration log — first pipeline run (2026-04-24)

The four-step pipeline from `HANDOFF.md` §7 was executed end-to-end.
**Outcome: negative for the minimal v3 ansatz.** The pipeline ran to
completion and produced a nominal `LOCKED_QUARK_PARAMS`, but the
residual against the observed masses is O(1), not O(10⁻³) as in the
lepton sector. This section records the numbers so the next iteration
has a starting point; §3.5 outcome #2/#3 is the structural reason.

### Step 1 — coarse grid scan
Command: `python scripts/calibrate_quark_ratios.py --n-points 8 --verbose`
Grid: 5 residual axes × 4 γ_q values = 131072 points.

- max rel err:          **9.99 × 10⁻¹** (anchor u=2.16 MeV locked)
- scan points:          131072
- rejected (unphysical):115160 (87.9% — positivity rejector firing hard)
- best γ_q:             0.050 (smallest in {0.05, 0.1, 0.2, 0.5})

Predicted vs observed (MeV):

| species | predicted | observed |
|---------|----------:|---------:|
| u | 2.16  | 2.16     |
| d | 3.95  | 4.67     |
| s | 181.3 | 93.4     |
| c | 183.1 | 1270     |
| b | 229.9 | 4180     |
| t | 233.3 | 172690   |

Interpretation: as `HANDOFF.md` warned, the coarse grid cannot discover
the heavy-sector uplift; `β` is still zero here, so the heavy sector
collapses to a tight cluster near the strange-region diagonal.

### Step 2 — integer-winding β hunt
Command: `python scripts/sweep_quark_beta.py --input-json /tmp/quark_step1.json --integer-min 10 --integer-max 1000000 --n-samples 600 --verbose`

- best integer winding: **N = 122** (log-spaced sample)
- β = N · π/2 =         191.64
- max rel err:          **9.53 × 10⁻¹**
- rejected:             0 / 556
- observed split:       max_rel_err decreases monotonically ≈ 4 × 10⁻⁴
  per unit N over N ∈ [10, 122]; larger log-spaced samples past N=122
  gave worse error (never triggered a new-best log line), so the basin
  is bounded.

Predicted vs observed (MeV) at N=122:

| species | predicted | observed | rel err |
|---------|----------:|---------:|--------:|
| u | 2.16    | 2.16     | 0       |
| d | 3.95    | 4.67     | 0.154   |
| s | 180.8   | 93.4     | 0.936   |
| c | 182.6   | 1270     | 0.856   |
| b | 8097.7  | 4180     | 0.937   |
| t | 8101.2  | 172690   | 0.953   |

The heavy uplift brings both k=5 states up to ~8100 MeV — bracketing b
and massively undershooting t. **The t/b ratio is essentially 1.0005 in
the prediction versus 41.3 in the observed spectrum.** This is the
structural signature of outcome #2/#3 from §3.5: the β · (k−3)² term
gives the SAME uplift to both partitions at k=5, so the b–t split must
come from `γ_q · σ · u_q(k=5) = γ_q · σ · 3`, and γ_q is held small by
the (1,+) positivity constraint.

### Step 3 — basin-width probe around N=122
Command: `python scripts/map_basin_quark_uplift.py --integer-winding 122 --half-width 50 --n-points 101 --input-json /tmp/quark_step1.json`

- best N in window:     **N = 123**
- best max rel err:     9.53 × 10⁻¹
- basin within 2× best: [72, 172]  (width 100)

The basin is wide and shallow — not a needle, but also not a well: the
floor is O(1), so the "attractor" is topologically real in the sense
that error changes smoothly across a 100-wide integer neighborhood,
but the attractor bottom is nowhere near a good fit.

### Step 4 — final lock (β and action_base hard-locked)
Command: `python scripts/lock_quark_beta_probe.py --integer-winding 123 --action-base-label pi --n-points 6 --verbose`

(Run with `--n-points 6` for tractability — the default n=12 is a
12⁶ = 2.99M-point grid that timed out at 20 min; n=6 is 46656 points.
The comparison against step 2 is still apples-to-apples because step 4
optimizes a log-RMS residual, not max rel err, and the tighter n=12
grid would not have rescued the underlying structural limitation.)

- integer_winding:      123
- action_base:          π = 3.14159…
- residual (log-RMS):   1.33
- **max rel err:        4.25**  (worse than step 2's max-rel-err
  minimum because step 4 targets log-RMS, not the L∞ error)
- rejected points:      42282 / 46656 (90.6% — positivity fires even
  harder on this grid because γ_q is extended up to 0.5)

Locked `QuarkParams` written to `quark_spectrum.py`:

```python
LOCKED_QUARK_PARAMS = QuarkParams(
    action_base=3.141592653589793,
    beta=193.20794819577227,       # 123·π/2
    gamma_q=0.108,
    u_q_form='k_minus_2',
    phase=0.0005,
    transport=1.24,
    pinhole=10.0,
    resistance=0.22,
    partition_mixing=0.0,
    winding_mode='max',
    resistance_model='exponential',
    depth_cost_mode='tunnel_only',
    spectrum_zero=None,
)
```

Predicted vs observed (MeV):

| species | predicted | observed | rel err |
|---------|----------:|---------:|--------:|
| u |     2.16 |      2.16 | 0       |
| d |     8.14 |      4.67 | 0.74    |
| s |    331.6 |     93.4  | 2.55    |
| c |    337.6 |   1270    | 0.73    |
| b |  21945   |   4180    | **4.25**|
| t |  21963   | 172690    | 0.87    |

### Verdict

Per `HANDOFF.md` §7 step 8: **the integer winding N≈122–123 is not a
clean topological invariant and the minimal v3 ansatz needs revision.**
The diagnostic signatures are:

1. **Heavy-sector rigidity.** Both b and t predictions sit within 0.1%
   of each other regardless of β, because `max(0, k−3)² = 4` identically
   at k=5 for both partitions. No β can open a t/b gap.

2. **Positivity squeeze.** Raising γ_q (which would split b and t)
   pushes the (1,+) diagonal below `action_base`, and the
   `extract_physical_spectrum` rejector fires. Rejection rate climbs
   from 88% (step 1) to 91% (step 4, wider γ_q range).

3. **Shallow basin, bad floor.** The integer-winding locus is a
   100-wide smooth attractor around N=123, but the attractor's floor
   is max-rel-err ≈ 0.95, not ≈ 10⁻³.

These three together point to §3.5 outcome #2 or #3 — the minimal
ansatz `u_q(k) = k − 2` combined with a k-dependent uplift that is
partition-degenerate at k=5 cannot reproduce the observed heavy-sector
hierarchy. Next-session candidates:

- Try `u_q(k) = k · (k − 2)` or `|k − 2|` variants (the former has the
  same sign problem at k=1; the latter breaks the generation-1
  inversion — see §3.5 outcome #2).
- Introduce a partition-dependent uplift term, e.g. `β · σ(p) · (k−3)²`
  or `β · (k−3)² · (1 + α σ)` — both would break the k=5 degeneracy.
- Reinterpret `action_base` as a typical rather than absolute floor,
  and take the spectrum zero to be `min(eigenvalues)`. This removes the
  positivity constraint and frees γ_q to grow.
- Treat the integer-winding hunt as a false positive and revert to
  searching over continuous β.

The nominal lock is recorded in `quark_spectrum.py` with a prominent
"NOT a successful fit" comment so downstream code that calls
`solved_quark_masses_mev()` does not silently claim a calibrated
result.

This section is the direct analog of the calibration-results section in
`docs/lepton_axioms.md`; it is what the community can cite when
discussing the discovered invariants — or, as in this first pass, the
discovered obstruction.

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
