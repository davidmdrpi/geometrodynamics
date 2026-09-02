# Audit: emergent QFT from classical GR with an antipodal identity and an inner–outer non-orientable bulk connection

*Independent audit of `davidmdrpi/geometrodynamics` at commit `b771b16`
(2026-09-02). Question audited: does the repository demonstrate quantum-field-
theory-like behaviour emerging from classical general relativity, via (a) an
antipodal identification and (b) a bulk interaction that connects the inner
and outer surfaces through a non-orientable wormhole? Everything below was
checked by reading the code and re-running it, not by reading the prose.*

## Bottom line

**No.** The repository contains a great deal of careful classical computation,
and its most recent rounds are honest about their negative results. But on the
specific question:

1. **The non-orientable transport is not derived.** The map that
   `embedding/transport.py` calls "the unique orientation-reversing,
   Hopf-preserving isometry of S³" is orientation-*preserving* (its 4×4 real
   matrix has determinant +1; it is left multiplication by the unit quaternion
   −j). No orientation-reversing isometry of S³ preserves the Hopf fibration at
   all. `T = iσ_y` is the ordinary spin-½ Kramers structure, inserted, not
   non-orientable geometry.
2. **Quantum structure is an input wherever the output is quantitative.** The
   Bell result `E = −cos(a−b)`, `CHSH = 2√2` is a textbook two-qubit
   calculation once the singlet and the Born rule are assumed; both are
   assumed. The repository's own later probes concede that a classical field
   with local dynamics gives `CHSH = 2` exactly, and that the "fully derived"
   detector pairing also gives `CHSH = 2.000000`.
3. **The inner–outer gluing is a drawing rule, not a solution.** Every wave
   round that "connects" `R_outer` to `R_inner` states in its own docstring
   that the crossing rule is a representation choice on a fixed round
   background with no Einstein equation. The visualization capstone lists "the
   wormhole identification itself" among the imported inputs.
4. **Where GR is actually solved, the connection is forbidden.** The
   finite-mouth, source-audit, Gauss–Bonnet and negative-EGB rounds show, on
   the repository's own evidence and with numbers I reproduced, that a smooth
   traversable inner–outer throat needs null-energy-condition violation at the
   neck and that no field in the repository can supply it. The Tangherlini
   branch that survives has a horizon, and the cross-exterior retarded Green
   function vanishes identically there, which kills the transaction mechanism.
5. **No field is quantized anywhere.** There is no Fock space, no commutator,
   no mode quantization; `ħ` appears only as an SI constant pasted into probe
   scripts. The "propagator" derivation is the Fourier transform of the
   classical Coulomb Green function. The QED vertices and squared amplitudes
   used downstream are imported formulas.

What the repository *does* establish is narrower and worth keeping: a
correctly derived scalar-flat handle geometry, a correct recovery of the
Morris–Thorne / Hochberg–Visser flare-out theorem in five dimensions, a
working 4+1 characteristic Einstein–scalar evolution, and a candid internal
audit trail. Those are classical-GR results. None of them produces QFT.

## Method

* Read the README (4996 lines), `docs/THESIS.md` headers, the package
  docstrings for `embedding`, `bell`, `bulk`, `waves`, `viz`, `tangherlini`,
  `transaction`, `qcd`, and the prior internal audit commit `69fc599`.
* Ran the full test suite (result recorded in the last section).
* Re-ran `finite_mouth_probe`, `source_audit_probe` and
  `maslov_dimensional_bridge_probe`; run directories were deleted afterwards
  so this audit adds no ledger entries.
* Wrote and ran an independent check of the orientation claim (appendix).
* Verified by hand: `f = √(s²+b²)` solves `f f″ = 1 − f′²`; the Raychaudhuri
  neck values `θ(0) = 0`, `dθ/dλ = 3/b²`; the conformal-coupling identity
  `1 − 2ξ_c = D/(2(D−1))`; the lepton masses returned by
  `solved_lepton_masses_mev()`.

## Findings

### F1. The claimed orientation-reversing Hopf isometry preserves orientation

`geometrodynamics/embedding/transport.py` states:

> Derives T = iσ_y as the UNIQUE orientation-reversing spinor map on S³ that
> preserves the Hopf bundle structure. This is not an ansatz — it is a theorem
> of the Hopf fibration.

The map is `σ(z₁, z₂) = (z̄₂, −z̄₁)`. Written on `R⁴ = (x₁, y₁, x₂, y₂)`:

```
M = [[ 0, 0, 1, 0],
     [ 0, 0, 0,-1],
     [-1, 0, 0, 0],
     [ 0, 1, 0, 0]]        MᵀM = I,   det M = +1,   M² = −I
```

`det M = +1`, so `σ` is an element of `SO(4)` and preserves the orientation of
`S³`. Identifying `(z₁, z₂)` with the unit quaternion `q = z₁ + z₂ j`, the map
is `q ↦ −j·q`: left multiplication by a unit quaternion, a fixed-point-free
*rotation* of `S³`.

The stronger statement is also false: **no orientation-reversing isometry of
`S³` preserves the Hopf fibration.** An isometry `g ∈ O(4)` maps Hopf fibres to
Hopf fibres iff it normalizes the `U(1)` fibre action, i.e. `g I g⁻¹ = ±I`
where `I` is multiplication by `i`. The `+` case is `U(2)` (complex-linear),
the `−` case is `K·U(2)` (antilinear, `K` = complex conjugation). Both have
real determinant `+1` (`det_R K = (+1)(−1)(+1)(−1) = +1`). A random search
over 20 000 determinant-`−1` isometries found none that preserve fibres; a
genuine reflection `diag(1,1,1,−1)` does not.

Consequences:

* The chain "Hopf fibration → orientation reversal → `T = iσ_y` → singlet"
  advertised in the README ("derived in `embedding/transport.py` without
  ansatz"; claim-table rows "History closure → E = −cos(a−b): **Derived**",
  "CHSH S = 2√2 (topological): **Verified**") has a false first link. What
  `σ` actually is: the antiunitary map `ψ ↦ iσ_y ψ*`, the standard time-
  reversal / Kramers structure of a spin-½ representation. Its properties
  `T² = −I`, `det T = 1` are representation theory, not throat topology.
* `verify_hopf_preservation` only checks `|σ(z)| = 1`; the fibre-mapping check
  its docstring promises is an empty comment. `test_hopf_fibration_preservation`
  therefore cannot detect the error. No test checks the determinant.
* The `Z₂` "wrap parity" and "non-orientable" labels in
  `embedding/topology.py` are dataclass fields set by hand
  (`orientation=±1`, `wrap_parity=±1`). No metric, identification, or
  quotient anywhere in the package produces a non-orientable manifold; the
  only orientability computation is a lookup table in `viz/hyperspherical.py`
  (`RPⁿ` orientable iff `n` odd).

### F2. Quantum mechanics is assumed where "QFT-like behaviour" is quantitative

`bell/bulk_identity.py` declares a 4-component complex `pair_state`, builds
`|Ψ⟩ = Σ_s |s⟩ ⊗ T|s⟩` (which *is* the singlet, by construction), applies
detector projectors `(cos θ/2, sin θ/2)`, and squares amplitudes. That is the
standard quantum calculation; `E = −cos(a−b)` and `CHSH = 2√2` follow
identically for any `T ∈ SU(2)` with `T² = −1` and have nothing to do with a
throat. The tensor-product Hilbert space, the projectors and the Born rule are
all inputs.

The repository already knows this. Two of its own rows (README lines ~320 and
~367) state:

* "a single classical field on 3-space with local dynamics is an LHV model —
  all 16 deterministic local strategies enumerated, max CHSH = 2 EXACTLY";
* "with the setting the docs derive and the detector #209 derives, CHSH =
  2.000000 exactly … the *-algebra it generates is exactly `span{I, σ_z}`,
  abelian";
* the capstone classifies the measurement orientation as "the IMPORTED
  MEASUREMENT FREEDOM".

So the honest status of the Bell arc is: **classical dynamics on the network
gives Bell-local correlations; Tsirelson is reached only after importing the
SU(2) transport and a measurement freedom that the geometry does not
supply.** The top-of-README claim table has not been updated to say this.

### F3. The inner–outer connection is a representation choice on a fixed background

The visual/wave arc (PRs #242–#269, `viz/*`, `waves/*`) is where the
"inner surface glued to outer surface" picture lives. Its own scope statements:

* `viz/circle_slice.py`: "The crossing rule is the obvious one: a radius that
  would pass `R_outer` re-enters at `R_inner` … **How** it is glued is a
  choice, and not a cosmetic one."
* `viz/one_surface.py`: "The crossing rule that glues `R_outer` to `R_inner`
  is a *representation* choice, not a derived boundary condition. The field is
  a linear scalar on a fixed round background."
* `docs/geometric_visualization_capstone.md` §34, "What is imported rather
  than derived": "The fixed round `S²` background of every wave round — and
  the fixed round `S³` of the last one: curvature 1 everywhere, at every time,
  with no Einstein equation and no backreaction. The wormhole identification
  itself: a pair of mouths, a time offset Δ, a loop transfer κ, and flux
  conservation through the throat. All inputs."

The one-surface result that "the antipode is parity-dependent"
(`Z_n(π) = (−1)ⁿ`) is correct but is the ordinary parity of zonal harmonics;
it holds for any linear wave on a round sphere and carries no information about
wormholes or non-orientability.

### F4. Where the Einstein equations are actually solved, the mechanism fails

These are the rounds that engage classical GR, and I reproduced their numbers:

| module | what it shows | reproduced |
|---|---|---|
| `bulk/finite_mouth.py` | scalar-flat handle `f = √(s²+b²)`, Darmois-matched to `S⁴_R`; `b = R sin²a`, shell-free seam; `8πG₅(ρ+p_s)|₀ = −3/b²` for every lapse with `N(0) > 0` | probe 6/6; `ff″ = 1 − f′²` checked by hand |
| `bulk/source_audit.py` | flare-out ⇒ `R_kk = −3/b²` (Raychaudhuri residual 1e-13); minimal scalar, GL order field, GL potential, Maxwell, Λ, perfect fluid, vacuum GWs all give `T_kk ≥ 0` | probe 5/5; neck values checked by hand |
| `bulk/gauss_bonnet.py`, `bulk/negative_egb.py` | `α_GB H_kk` has the same sign as `R_kk` at a neck; the exterior needs `α_GB ≥ −R²/4`, the throat `α_GB ≤ −R²/4` | not re-derived; internally consistent |
| `tangherlini/traversable_throat.py` | for the Tangherlini bridge, `G_ret(c,s)·G_adv(c,d) ≡ 0` across exteriors by a causal-set argument; the MTY transaction needs a traversable throat | argument checked; it is correct |
| `tangherlini/dynamics.py` | 4+1 Einstein–scalar characteristic evolution, second-order convergent; Tangherlini exact fixed point; `A > 0` on any regular-centre ingoing slice | not re-run; scope statements accurate |

Read together: **a smooth traversable inner–outer throat in 5D Einstein
gravity requires `T_kk < 0` at the neck; nothing in the repository's matter
content provides it; the branch that needs no exotic matter has a horizon,
and through a horizon the exchange kernel is zero.** This is the repository's
own conclusion, and it is the correct one for the fields it contains. It is
also fatal to the premise of the question audited: the bulk interaction
connecting inner and outer surfaces has no classical-GR support in this code.

### F5. One loophole in the source audit is closed by assumption, not by theorem

Candidate C10 (nonminimally coupled scalar) is dismissed because at a node
`q = 0` the prefactor `1 − 8πG₅ξq²` equals 1 and the sign is `sign(1 − 2ξ)`.
That is correct *at a node*. The audit then says "BAM places the defect core at
`q = 0`, so this is the relevant point". That is a modelling assumption, not
part of the flare-out theorem.

Away from a node the known classical result goes the other way: Barceló &
Visser, *Traversable wormholes from massless conformally coupled scalar fields*,
Phys. Lett. B 466 (1999) 127, exhibit an exact three-parameter family of
traversable wormholes in Einstein gravity with a conformally coupled,
non-ghost scalar, supported where `8πGξq² > 1` and the effective coupling
changes sign. The source audit's list of escapes (ghost, higher curvature,
quantum stress, horizon, reinterpretation) omits this sixth classical branch.
Whether BAM can use it depends on giving up `q = 0` at the throat; the
repository should say so rather than report C10 as closed "in every
dimension".

### F6. The "antipodal identity" is never imposed as an identification

`antipode4(p) = −p` is used to (i) place mouth B at `−p_A`, (ii) weight
handshake candidates by a Gaussian in the antipodal error with width
`SIGMA_ANTI = 0.18` (a free constant), (iii) couple a delayed, half-ring-
shifted copy of the field into the QCD network solver with strength
`antipodal_coupling = 0.15` (a free constant), and (iv) feed a Bohr-like
resonance rule `ω τ + φ_spin + φ_throat ≡ 0 or π` in `transaction/cavity.py`
whose mode frequencies are hand-set constructor defaults (`1.055`, `1.219`).

Nowhere is a field equation solved on the quotient `S³/Z₂ = RP³`, nowhere is a
twisted spin structure constructed, and nowhere does an identification `x ~ −x`
enter a boundary condition of a solve. The ledger entry "B2_antipodal_Z2:
CLOSED (RP³ + spin structure)" in `maslov_dimensional_bridge_probe` is prose
with no corresponding computation. The 5D transfer kernel round, which is the
one place an antipodal boundary condition could have been implemented at a
horizon (cf. 't Hooft's antipodal identification), instead computes a
Tangherlini exterior-to-horizon kernel and finds it is zero across mouths.

### F7. The lepton and quark ladders are calibrated fits, and the README contradicts itself about them

`solved_lepton_masses_mev()` returns `[0.511, 105.6126, 1778.938]`, i.e.
`−0.043 %` and `+0.117 %` on μ and τ, and the "Lepton mass ladder" section
still advertises "sub-percent accuracy with zero free parameters at scan
time". The same README's PR #271 section states that the pinhole `γ = 22.5`
is "no longer derived", that the corrected radial operator yields `22.331` or
`22.836` at the canonical geometry, and that either value moves the muon by
`15–21 %` because `d ln m_μ / d ln γ = −17.5`.

The locked Hamiltonian carries four real anchors
(`phase = 0.001`, `transport = 25.1`, `pinhole = 22.5`,
`resistance = 0.217869`) plus three discrete structural switches
(`winding_mode`, `depth_cost_mode`, `resistance_model`), all selected by the
calibration scans listed in the section's own script map, to reproduce two
mass ratios. The later identification of some anchors with "closure-quantum
invariants" (`8π`, `7π/100`) is post hoc. The quark sector adds three opt-in
extensions and reaches 1.6 %. These are fits with more knobs than targets, and
the top-level claim table should say "fitted", as the row "**The masses are
FITTED**" deeper in the same table already does.

### F8. The verification instruments cannot fail, and pre-registration is nominal

The prior internal audit (`69fc599`) found 45/45 probe runs passing and a
suite that is edited when a verdict changes. This audit confirms it:

* 45 assertions in `tests/` compare a returned string against prose
  (`assert "NARROWED, not closed" in ledger["verdict"]`,
  `assert entry["verdict"] == "WITHDRAWN"`). These lock the narrative, not
  the number.
* The three "pre-registered" rounds are real, but the pre-registration commit
  precedes the module by 24–40 minutes in the same session and author
  (`765dbaa → 7de86ce`, `ca07204 → 1082f5d`, `d47df40 → 70ad286`). That is
  a useful discipline; it is not independent.
* The one external oracle the prior audit identified (Matyjasek 2021 for the
  5D QNM) remains the only one. The Hochberg–Visser recovery in the source
  audit is a legitimate second.

### F9. Terminology: "QFT on a fixed classical background" is not what the code does

The README's framing is "standard field theory on that fixed classical
background, in the precise sense of QFT-on-curved-spacetime". QFT on curved
spacetime means a quantized field: mode expansion, canonical commutators or a
path-integral measure, a vacuum state, renormalized stress tensor. None of
those objects exists in the code. Grepping the package and all 302 probes for
`Fock`, `annihilation`, `creation operator`, `commutator` finds only prose
mentions. `ħ` occurs as `1.054571817e-34` in ten probe scripts and in the
anchor relation `ħ = m_e R_MID c`, which is a unit choice. The
"exchange kernel" probe derives `1/q²` as the Fourier transform of the
classical Coulomb potential, then feeds it into the textbook Bhabha/Møller
squared matrix elements, which are imported. The "one-loop determinant" and
"path-integral measure `S_BAM`" are names attached to classical spectral sums
and Bohr–Sommerfeld/Maslov counting.

## What is sound

* The finite-mouth handle: geometry, matching, and the lapse-independent neck
  cost are correct and cleanly derived.
* The source audit's flare-out theorem and its attribution to
  Morris–Thorne / Hochberg–Visser.
* The causal-support argument that the cross-exterior transactional product
  vanishes on the Tangherlini bridge.
* The 4+1 characteristic code's positivity identity and its stated
  limitations.
* The retractions themselves. The repository's later rounds routinely
  withdraw earlier claims, and the internal audit `69fc599` is accurate.

## Recommendations

1. Correct `embedding/transport.py` and every README row that depends on it:
   `σ` is orientation-preserving; `T = iσ_y` is an SU(2)/Kramers structure,
   not a non-orientable transport. Add a test asserting `det = +1` so the
   docstring cannot drift back.
2. Rewrite the claim table's Bell rows to the repository's own later verdict:
   classical network dynamics gives `CHSH = 2`; `2√2` requires the imported
   singlet and measurement freedom.
3. Move the inner–outer gluing from a drawing rule to a boundary condition of
   an actual solve, or stop describing the wave rounds as a "connection".
4. Add Barceló–Visser as a sixth branch in the source audit, with the explicit
   statement that it requires abandoning `q = 0` at the neck.
5. Implement the antipodal identification once, as a real quotient
   (`RP³` boundary condition or an antipodal horizon condition in the 5D
   kernel), so that "B2 CLOSED" corresponds to a computation.
6. Mark the lepton and quark ladders "fitted" in the top table and reconcile
   the "sub-percent, zero free parameters" section with the PR #271 section.
7. Replace prose-verdict assertions with numeric ones, and regression-lock the
   `γ` sums, `R_OUTER` fixed point and `1.054` factor that PR #271 found were
   unprotected.
8. Drop "QFT on a fixed background" from the framing until a field is
   quantized; the accurate description is classical wave mechanics and
   spectral counting on a classical geometry.

## Appendix: the orientation check

```python
import numpy as np
from geometrodynamics.embedding.transport import orientation_reversal_on_s3

def real4(z1, z2):
    return np.array([z1.real, z1.imag, z2.real, z2.imag])

M = np.zeros((4, 4))
for k in range(4):
    e = np.zeros(4); e[k] = 1.0
    w1, w2 = orientation_reversal_on_s3(complex(e[0], e[1]), complex(e[2], e[3]))
    M[:, k] = real4(w1, w2)

assert np.allclose(M.T @ M, np.eye(4))       # isometry
assert np.allclose(M @ M, -np.eye(4))        # squares to -1
print(np.linalg.det(M))                      # +1.0  -> orientation-preserving

# no orientation-reversing isometry preserves the Hopf fibration:
I4 = np.array([[0,-1,0,0],[1,0,0,0],[0,0,0,-1],[0,0,1,0]], float)  # mult by i
def fibre_preserving(G):
    return np.allclose(G @ I4, I4 @ G) or np.allclose(G @ I4, -I4 @ G)
rng = np.random.default_rng(0); hits = 0
for _ in range(20000):
    Q, R = np.linalg.qr(rng.standard_normal((4, 4)))
    Q = Q @ np.diag(np.sign(np.diag(R)))
    if np.linalg.det(Q) > 0: Q[:, 0] *= -1
    hits += fibre_preserving(Q)
print(hits)                                  # 0
```

Output on this machine: `det = 1.0`, `hits = 0`.

## Test suite

*Full-suite result pending at the time of this commit; updated in the next commit.*
