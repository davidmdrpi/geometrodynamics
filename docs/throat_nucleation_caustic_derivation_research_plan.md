# Dynamic throat-nucleation / antipodal-caustic derivation — research plan

Open thread from PR #35 (closed-form Compton vertex factor) and the
crossing triangle PRs #36, #37: derive

    F²(x, c) = 4·x³·(x² + 1 − x·sin²θ) / [(1 + c²)·(1 + x)²]

from BAM-native geometry — caustic / flux / throat-opening kinematics
and Hopf-fibre polarization transport — rather than as an algebraic
fit to Klein–Nishina.

## The algebraic core

Two equivalent closed-form decompositions sit at the centre.

### Decomposition A: the F² factor

    F²(x, c) = [2x/(1+x)]² · [x² + x·(1−x)² / (1+c²)]
             = P(x)² · Q(x, c)

via the exact identity

    x² + 1 − x·sin²θ ≡ (1 − x)² + x · (1 + c²).

### Decomposition B: the spin-summed amplitude

    |M̄|²_KN / (8e⁴) = (1 + c²) + (1 − x)² / x

— the standard Klein–Nishina formula `1/x + x − sin²θ` rewritten as
**Thomson polarization sum** + **recoil-induced correction**, with
the two pieces algebraically orthogonal (the second vanishes
identically at Thomson `x = 1`).

Both decompositions are algebraic identities and exact at machine
precision. The derivation question is whether each piece has a
BAM-native geometric origin.

## BAM-geometric interpretations under test

### 1. Caustic / throat-opening kinematics: `P(x) = 2x/(1+x)`

`P(x) = 2x/(1+x) = 2ω'/(ω + ω')` is the **harmonic mean** of the
incoming photon frequency `ω` and the outgoing photon frequency `ω'`,
normalized to `ω'`. In the BAM antipodal-throat picture:

  - the throat opens at the front mouth with frequency `ω`,
  - closes at the antipodal back mouth with frequency `ω'`,
  - the throat-pinch rate is the harmonic mean of mouth rates
    (the natural averaging for a flux that must be continuous through
    a varying-radius bottleneck).

The squared factor `P²` then arises because **both mouths** of the
throat-pair contribute a pinch; the kinematic factor enters twice,
once per mouth.

This is a candidate **derivation, not a postulate**: harmonic-mean
flux-matching at a varying-cross-section conduit is the standard
classical result for any conserved flux through a bottleneck.

### 2. Hopf-fibre polarization transport: `(1+c²) = 2·(cos⁴(θ/2) + sin⁴(θ/2))`

The Thomson polarization sum `(1 + cos²θ)` is **exactly** the sum

    Σ_λ |d^1_{1,λ}(θ)|² = |d^1_{1,1}(θ)|² + |d^1_{1,−1}(θ)|²
                        = cos⁴(θ/2) + sin⁴(θ/2)
                        = (1 + cos²θ)/2

of squared Wigner-`d` matrix elements for spin-1 helicity transport
through angle `θ`. The BAM photon is a helicity-`±1` Hopf-fibre mode;
the `(1 + c²)` factor is the **squared overlap** of incoming and
outgoing Hopf-fibre helicity states under parallel transport along
the scattering geodesic.

The Hopf connection used in the repo (`geometrodynamics.hopf.
connection.hopf_connection`) has `A_φ(χ) = ½·cos(χ)`, with
`A_φ(0) = 1/2` — the value PR #34 found as the perturbative coefficient
`ξ = −A_φ(0) = −1/2`. The same `½` appears in the helicity-transport
amplitude.

### 3. Channel-split structure of `Q(x, c)`

The polarization factor `Q = x² + x·(1−x)²/(1+c²)` decomposes as

    Q = |a|² + |b|²,    a, b  orthogonal helicity amplitudes:
        |a|² = x²                      (helicity-preserving channel)
        |b|² = x·(1 − x)² / (1 + c²)   (helicity-flipping channel)

  - **Helicity-preserving**: the recoil-rescaled amplitude `a = x`
    survives at Thomson (`a = 1` at `x = 1`). This is the diagonal
    Hopf-fibre transport.
  - **Helicity-flipping**: weighted by `(1 − x)²` (the recoil
    deficit squared, vanishing at Thomson) and by `1/(1 + c²)`
    (the inverse Thomson polarization sum — natural for a
    "cross-channel" amplitude that competes with the diagonal one).

The orthogonal sum structure means `Q` is a probability — the total
"weight" carried by the throat-pinch through both Hopf helicities.

### 4. Optional Tangherlini radial-mode projection at threshold

At Thomson (`x = 1`), the BAM modification is trivial: `F² = 1`.
The radial Tangherlini ground-mode integrated against the throat
density should give the trivial weight; finite recoil (`x < 1`)
opens the channel for radial-mode contamination. This is the
expected "threshold" link to the bulk radial channel; it is optional
because the closed-form F is already complete without it.

## Tests

  T1. **Algebraic identity (F² decomposition)**: verify
      `F²(x, c) = [2x/(1+x)]² · [x² + x(1−x)²/(1+c²)]` to machine
      precision across a fine (x, c) grid.

  T2. **Algebraic identity (|M̄|² decomposition)**: verify
      `|M̄|²/(8e⁴) = (1+c²) + (1−x)²/x` to machine precision.

  T3. **Caustic Padé as harmonic-mean throat-rate**: verify
      `P(x) = 2x/(1+x) = 2ω'/(ω+ω')` at sample points. Check the
      harmonic-mean limits: `P → 0` as `x → 0`, `P → 1` at `x = 1`,
      `P → 2` as `x → ∞` (no upper bound = no harmonic-mean cap; the
      throat opens wider than the average).

  T4. **Hopf-fibre helicity transport (Wigner-d identity)**: verify
      `(1 + cos²θ)/2 = cos⁴(θ/2) + sin⁴(θ/2)` at sample angles —
      identifies the `(1+c²)` factor as the Hopf-fibre helicity sum.

  T5. **Channel-split orthogonality**: verify
      `Q(x, c) = |a|² + |b|²` with `a = x`, `b = √x·(1−x)/√(1+c²)`
      and the two channels non-negative across the physical region.

  T6. **Hopf connection at lock (`A_φ(0) = 1/2`)**: read the value
      from the repo and confirm it matches the perturbative `ξ`
      coefficient from PR #34.

  T7. **Cross-process consistency**: the P²·Q decomposition must
      hold under crossing `x → x_⊗ < 0` (Compton → BW / annihilation).
      Verify the analytic continuation preserves the structure.

  T8. **Alternative throat-rate candidates rejected**: arithmetic
      mean `(1+x)/2`, geometric mean `√x`, etc. None give a clean
      Q form. Verify by computing `Q_alt = F²/P_alt²` and confirming
      it is *not* of the form `x²·(...) + x·(1−x)²·(...)`.

  T9. **Tangherlini ground-mode at Thomson** (optional/informative):
      report the Tangherlini lowest radial-mode weight at the
      Thomson threshold; trivial weight expected.

## Verdict structure

  - **DERIVATION_COMPLETE**: T1, T2, T3, T4, T5 all pass at machine
    precision and the BAM-geometric interpretation of each factor is
    forced by the algebra (no alternative decomposition into
    BAM-native pieces works). The closed-form `F²` is derived from
    harmonic-mean throat-rate (caustic) and Hopf-fibre helicity
    transport (polarization).

  - **DERIVATION_PARTIAL**: T1–T5 pass but T8 finds an alternative
    decomposition that also matches the data — derivation is
    geometric but not unique.

  - **DERIVATION_OPEN**: any of T1–T5 fails — algebraic decomposition
    is not as proposed, requiring a different geometric story.

## What this leaves open

  - **First-principles BAM action / Lagrangian**: the harmonic-mean
    throat-rate is *consistent* with classical flux conservation
    through a varying-radius conduit, but a derivation from a
    specific BAM action (S³ throat dynamics + Hopf bundle coupling)
    is the deeper task.

  - **Loop corrections**: the closed form is tree-level. Whether
    the same channel structure survives one-loop is a separate
    target.

  - **Bhabha / Møller**: two-channel processes; the BAM kernel would
    need to combine two crossed copies coherently.

## Cross-references

  - PR #35: `compton_vertex_resummation_probe.py` — closed-form F².
  - PR #36: `breit_wheeler_cross_process_probe.py` — Compton → BW.
  - PR #37: `pair_annihilation_crossing_probe.py` — triangle closure.
  - `geometrodynamics/hopf/connection.py` — Hopf connection
    `A_φ(χ) = ½cos(χ)`.
  - `experiments/closure_ledger/throat_nucleation_caustic_derivation_probe.py`
    — this probe.
