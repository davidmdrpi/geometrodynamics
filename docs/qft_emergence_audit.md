# Audit: emergent QFT from classical GR with an antipodal identity and an inner–outer non-orientable bulk connection

*Independent audit of `davidmdrpi/geometrodynamics` at commit `b771b16`
(2026-09-02), revised after review. Question audited: does the repository
demonstrate quantum-field-theory-like behaviour emerging from classical
general relativity, via (a) an antipodal identification and (b) a bulk
interaction that connects the inner and outer surfaces through a
non-orientable wormhole? Everything below was checked by reading the code and
re-running it, not by reading the prose.*

*The **Bottom line** below is the **current** status, updated through PR #280.
The audit's original findings F1–F9 and its original bottom line are kept
verbatim as the historical record, and the **Revision note** at the end
records what the first draft got wrong, what has since been resolved, and how
each correction was verified.*

## Bottom line — current status (through PR #280)

**Not yet, but the gap is now localised to one place.** Two of the three
originally-missing arrows have been substantially closed, conditionally; the
third has been narrowed from "the quantization map is not derived" to three
named inputs plus composition.

1. **Topology → transport: substantially resolved, conditionally.** The
   `RP²` mouth's induced Pin⁻ holonomy squares to `−1` — the spin holonomy of
   a `2π` rotation on the round neck, forced by the ambient Pin⁺ structure and
   the two twisted normal lines — and lies in the same fibre-reversing
   component of `Pin(2) ⊂ SU(2)` that contains `J = iσ_y K`, up to a `Spin(2)`
   direction that is gauge and a sign that is the Pin⁻ sector. The Hopf fibre
   **is** the mouth's `Spin(2)` (`q ↦ (q⁻¹iq, q⁻¹jq, q⁻¹kq)` is
   `Spin(3) → SO(3) = F_SO(S²)`; fibre angle `φ` = frame angle `2φ`;
   `c₁ = 1` against `e = 2`), so the antilinear `K` is canonical rather than
   inserted. `J² = −1` is geometric. Conditional on the antipodal quotient
   construction. (`docs/finite_mouth_bundle_transport.md`,
   `docs/mouth_spin_frame_prereg.md`; PR #279.)
2. **Finite mouth → antipodal boundary condition: substantially resolved,
   conditionally.** PR #129's `Φ(U,V,Ω) = Φ(−U,−V,Ω̄)` is the `η = +1` scalar
   sector of the unique free involution `ι(s,Ω) = (−s,−Ω)` of the finite
   handle, obtained at the finite ultrastatic neck with no horizon and no
   limit, and reproduced against the PR #277 admittance oracle to `5e-7`.
   Conditional on `P_B = −P_A`, on taking the quotient rather than the double
   cover, and on `η`. (Same documents; PR #279.)
3. **The transport bug is fixed.** `embedding/transport.py` no longer claims
   `σ` is orientation-reversing on `S³`; it is `L_{−j}`, the base antipode
   composed with fibre reversal, `det = +1`, with Spin(4) lift `(−j, 1)`.
4. **Traversable causal exchange remains unsupported**, and that is unchanged:
   5D Einstein gravity with the repository's matter cannot hold a smooth
   traversable neck open, the Tangherlini bridge does not transmit across
   exteriors, and constant-coupling Gauss–Bonnet reaches its critical value
   only where the linearised principal symbol degenerates. PR #206's
   non-traversable global-constraint reading is untouched by this and is now
   the branch being developed.
5. **The current substantive gap is the probability and composition problem,
   in four named parts.** Sharp phase closure of the closed history *does*
   make the source ensemble depend on both analyzer settings, with exact
   detector no-signalling and `S = 2.1423 > 2` — a Bell violation computed
   from a classical global constraint rather than an inserted singlet — and
   the holonomy-weighted current of the derived loop gives the quantum joint
   law `P(s_A,s_B) = (1 − s_A s_B cos γ)/4`, `S = 2√2`, with no projectors.
   What remains underived:
   - **branch aggregation** — positive count `|D|` (`S_max = 2.1423`) or the
     oriented, holonomy-weighted sum `D` (`2√2`);
   - **the relative sector coefficients** — the equal outcome-sector prior is
     a chosen counting measure and moves the correlation on *both* branches;
     the quantum law appears only at `r = 1`;
   - **the readout** — why observed frequencies would be *linear* in the
     integrated current rather than quadratic in it;
   - **composition** — everything above is one pair; `Γ₁ × Γ₂` versus
     `H₁ ⊗ H₂` is untouched.
   (`docs/closure_measurement_dependence.md`, `docs/closure_current.md`,
   `docs/classical_born_rule.md`; PR #280.)
6. **The action route to branch aggregation is closed, and the causality gate
   fails.** The holonomy trace `S_H = −½ Tr G = −cos(Ω/2)` has the derived
   closure set as its **Morse–Bott critical manifold**, and its transverse
   stationary-phase magnitude is exactly the positive coarea density — with
   nothing tuned. It is nevertheless not an action: closure holonomies at a
   common basepoint compose additively in `θ`, so `e^{iθ}` is a genuine
   holonomy carrying the branch sign but has **no critical points**, while
   `S_H` has the critical set and is **not additive**, so `e^{iκS_H}` is not
   the holonomy of any connection and `κ` is free. `A_π/A_0 = e^{2iκ}e^{∓iπ/2}`
   has magnitude `1` for every real `κ`, so the fork is purely a phase
   question; `κ = 1` selects neither branch. Separately: the fixed-setting
   symmetry group has two sector orbits except at `γ = π/2`, so `r` is free at
   every CHSH angle; every classical coupling in the repository is degree-2
   homogeneous, and a quadratic readout is *superquantum*
   (`S_max = 8√2/3 = 3.7712`); and **operational no-signalling to the past
   fails** — the conditioned source density is exactly antipodally even, so odd
   observables are blind, but every even one signals, and non-coplanar settings
   give mutually singular supports (total variation `1`). Since every BAM
   coupling is even, the readable class is exactly the one BAM has.
   (`docs/history_action.md`; round 7.)

So the organising question is no longer whether the geometry supplies spin
structure or a Bell-violating global constraint — it does, conditionally on
named choices — but whether anything in the classical theory selects how
closed histories aggregate into observed event frequencies. Round 7 sharpens
the negative: the branch holonomy and the coarea magnitude come from two
provably different functionals, and no single object supplies both.

### Original bottom line, as written at `b771b16`

*Kept for the record; items 1 and 2 are superseded by the current status
above, and the findings F1–F9 below are the audit as originally written.*

**Not yet.** The repository contains a large body of careful classical
computation and an unusually candid retraction record, but the specific
chain "classical GR + antipodal identity + non-orientable inner–outer
connection ⟹ QFT" is not closed. Four statements summarise the state:

1. **A real bug in the transport derivation.** `embedding/transport.py`
   identifies `J = iσ_y K` with "the unique orientation-reversing
   Hopf-preserving isometry of S³". The map it writes down has real
   determinant `+1`: it is the antipode on the Hopf base composed with a
   reversal of the Hopf fibre, and the product preserves the orientation of
   `S³`. No orientation-reversing isometry of `S³` preserves the Hopf
   fibration. `J` remains a natural candidate for the lift of a non-orientable
   mouth transition, but the derivation the README cites for it is false, and
   the theorem that the physical `RP²` mouth transition lifts to `J` has not
   been proved.
2. **The antipodal identification exists as a boundary condition, not as a
   derived geometry.** PR #129 imposes `Φ(U,V,Ω) = Φ(−U,−V,Ω̄)` at the
   Tangherlini horizon and derives even-ℓ Neumann / odd-ℓ Dirichlet throat
   conditions; PR #135 builds the matter resolvent on them. What is missing is
   any derivation of that identification from the finite-mouth 5D geometry of
   PR #277, or a solve of the Einstein system on a global antipodal quotient.
3. **Traversable causal exchange through the throat is in genuine trouble.**
   With the fields the repository contains, 5D Einstein gravity cannot hold a
   smooth traversable neck open; the Tangherlini bridge does not transmit a
   retarded signal across exteriors; constant-coupling Gauss–Bonnet reaches
   its critical value only where the linearised principal symbol degenerates.
   This kills the Morris–Thorne–Yurtsever transaction mechanism. It does not
   by itself kill the alternative already in the repository (PR #206), where
   the bridge is non-traversable and the pair identity is a global constraint
   imprinted at nucleation; that alternative is simply underived.
4. **The quantization map is not derived.** The repository has both explicit
   Hilbert-space calculations (`bell/bulk_identity.py`) and a formal
   path-integral arc (`Z = Σ ∫ (dL/L) det^{−1/2} e^{−S_BAM}`, PRs #74,
   #115–#122). Nothing derives from the classical Einstein equations why
   histories are weighted by amplitudes, why probabilities are `|A|²`, or why
   the one-object classical configuration space becomes a complex tensor-
   product Hilbert space with the measurement structure Bell needs. The
   repository's own probes record that the fully derived detector pairing
   gives `CHSH = 2.000000`. `THESIS.md` contradicts itself on whether a path
   integral is performed.

## Method

* Read the README (4996 lines), `docs/THESIS.md`, the package docstrings for
  `embedding`, `bell`, `bulk`, `waves`, `viz`, `tangherlini`, `transaction`,
  `qcd`, the PR #129/#135 probes in `experiments/closure_ledger/`, and the
  prior internal audit commit `69fc599`.
* Ran the full test suite (result in the last section).
* Re-ran `finite_mouth_probe`, `source_audit_probe` and
  `maslov_dimensional_bridge_probe`; run directories were deleted afterwards
  so this audit adds no ledger entries.
* Wrote and ran an independent check of the orientation claim (appendix).
* Verified by hand: `f = √(s²+b²)` solves `f f″ = 1 − f′²`; the Raychaudhuri
  neck values `θ(0) = 0`, `dθ/dλ = 3/b²`; the conformal-coupling identity
  `1 − 2ξ_c = D/(2(D−1))`; `h(σz) = −h(z)` and `σ(e^{iφ}z) = e^{−iφ}σ(z)`
  for the Hopf map `h`; the lepton masses returned by
  `solved_lepton_masses_mev()`.
* Checked GitHub Actions for the CI history claim (run 32330296533).

## Findings

*These are the findings as written at `b771b16`. Three of them have since
changed status — see "What has been resolved since the audit was written" in
the Revision note, and the current bottom line at the top.*

### F1. The claimed orientation-reversing Hopf isometry preserves orientation; the transport derivation is reopened, not deleted

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

`det M = +1`, so `σ ∈ SO(4)` and it preserves the orientation of `S³`. With
`(z₁, z₂) ↔ q = z₁ + z₂ j`, the map is `q ↦ −j·q`: left multiplication by a
unit quaternion, a fixed-point-free rotation.

What `σ` actually does, with `h: S³ → S²` the Hopf map:

```
h(σz) = −h(z)                    antipode on the base S²   (orientation-reversing on S²)
σ(e^{iφ} z) = e^{−iφ} σ(z)       reverses the fibre S¹     (orientation-reversing on S¹)
```

Two orientation reversals whose product preserves the orientation of the total
space. That is exactly why `det = +1`, and it is a more useful description of
`σ` than the one in the file: `J = iσ_y K` is the canonical quaternionic
structure on `C²`, covering the base antipode, reversing the fibre, with
`J² = −1`.

The stronger statement also holds: **no orientation-reversing isometry of
`S³` preserves the Hopf fibration.** `g ∈ O(4)` maps fibres to fibres iff
`g I g⁻¹ = ±I` with `I` = multiplication by `i`; the `+` case is `U(2)`, the
`−` case is `K·U(2)`, and both have real determinant `+1`. A random search
over 20 000 determinant-`−1` isometries found none that preserve fibres.

Consequences, stated at the right strength:

* The README chain "Hopf fibration → orientation reversal of S³ → `T = iσ_y`
  → singlet" has a false link. The claim-table rows that cite it as *derived*
  ("History closure → E = −cos(a−b): **Derived**", "CHSH S = 2√2
  (topological): **Verified**") are not supported by this derivation.
* The repository already separates the orientability questions: the spatial
  `RP³` quotient is orientable (`viz/hyperspherical.py`), and the
  non-orientable / Pin structure is assigned to the `RP²` mouth/exchange
  sector and to a flat-bundle twist class on the `Z₂` network, not to `RP³`.
  `J` is the natural candidate lift of base-antipode-plus-fibre-reversal.
  **The missing theorem is that the physical `RP²` mouth transition lifts to
  `J = iσ_y K`.** F1 reopens that derivation; it does not discard `J`.
* `verify_hopf_preservation` only checks `|σ(z)| = 1`; the fibre-mapping
  check its docstring promises is an empty comment, so
  `test_hopf_fibration_preservation` cannot detect the error. No test checks
  the determinant.
* In `embedding/topology.py` the `orientation` and `wrap_parity` labels are
  dataclass fields set by hand; no metric or quotient in the package produces
  them.

### F2. `bulk_identity.py` assumes the quantum rules it claims to derive

`bell/bulk_identity.py` declares a 4-component complex `pair_state`, builds
`|Ψ⟩ = Σ_s |s⟩ ⊗ T|s⟩` (the singlet, by construction), applies the spin-½
projectors `(cos θ/2, sin θ/2)`, and sets `P = |A|²`. The tensor product, the
projectors and the Born rule are inputs; for any `T ∈ SU(2)` with `T² = −1`
the result `E = −cos(a−b)`, `CHSH = 2√2` follows identically and says nothing
about a throat. The module's statements that the law follows "from pure
topology alone … with nothing else" are untenable.

The repository's own later rounds already contain the sharper physical fact
(README rows at lines ~320 and ~367):

* a single classical field on 3-space with local dynamics is an LHV model,
  `max CHSH = 2` exactly over all 16 deterministic strategies;
* with the setting the docs derive and the detector PR #209 derives, the
  generated algebra is `span{I, σ_z}`, abelian, and `CHSH = 2.000000`;
* a winding-2 carrier (Probes L/M) can supply the missing non-commuting
  operation and operational models then violate Bell; the *geometric
  production* of that carrier/apparatus is classified as "the IMPORTED
  MEASUREMENT FREEDOM".

Accurate status: **BAM has identified geometric structures compatible with
the required measurement algebra, but has not derived the full Bell
experiment from GR.** The top-of-README claim table has not been updated to
say this.

### F3. Three different inner–outer objects, and only one of them is a drawing

The audit question conflates three constructions the repository keeps
separate, and they should be judged separately:

| object | where | status |
|---|---|---|
| visual crossing `R_outer → R_inner` | `viz/circle_slice.py`, `viz/one_surface.py`, PRs #242–#269 | representation convention on a fixed round background, by the modules' own docstrings and `docs/geometric_visualization_capstone.md` §34 ("the wormhole identification itself … All inputs") |
| horizon antipode `Φ(U,V,Ω) = Φ(−U,−V,Ω̄)` | `null_throat_boundary_conditions_probe.py` (PR #129), `antipodal_horizon_exchange_kernel_probe.py` (PR #135) | imposed physical boundary condition on a linear matter-wave problem, on the Tangherlini background; yields even-ℓ Neumann / odd-ℓ Dirichlet, a unitary mirror, a real ℓ-graded spectrum, and a self-adjoint resolvent |
| finite mouth `S⁴_R` handle | `bulk/finite_mouth.py` (PR #277) | classical geometric junction construction, Darmois-matched, with the discrete identification deliberately left out |

The one-surface result that "the antipode is parity-dependent"
(`Z_n(π) = (−1)ⁿ`) is correct but is ordinary zonal-harmonic parity; it is
the same `(−1)^ℓ` that PR #129 uses, and it carries no information about
non-orientability on its own.

The real gap is the composition
**finite-mouth geometry → orientation/topology lift → field boundary
condition.** The horizon identification of PR #129 is a postulate (its own
docstring: "fixed by the BAM antipodal postulate (PR #128)"), and PR #277
explicitly declines to fold `Φ_spatial`, `(−1)^ℓ`, `η_orientation`,
`η_wrap`, `U_spin` into the metric. Nothing connects the two.

### F4. Traversable causal exchange is unsupported; a global constraint is not ruled out

Rounds that engage the Einstein equations, with numbers reproduced here:

| module | what it shows | reproduced |
|---|---|---|
| `bulk/finite_mouth.py` | scalar-flat handle `f = √(s²+b²)`, Darmois-matched to `S⁴_R`; `b = R sin²a`; shell-free seam; `8πG₅(ρ+p_s)|₀ = −3/b²` for every lapse with `N(0) > 0` | probe 6/6; `ff″ = 1 − f′²` by hand |
| `bulk/source_audit.py` | flare-out ⇒ `R_kk = −3/b²` (Raychaudhuri residual 1e-13); minimal scalar, GL order field, GL potential, Maxwell, Λ, perfect fluid, vacuum GWs all give `T_kk ≥ 0` | probe 5/5; neck values by hand |
| `bulk/gauss_bonnet.py`, `bulk/negative_egb.py` | `α_GB H_kk` has the sign of `R_kk` at a neck; exterior needs `α_GB ≥ −R²/4`, throat `α_GB ≤ −R²/4`; at the critical value the linearised principal symbol `P = C_t ω² + C_s κ²` has `C_t = 0`, so the system stops being an evolution equation, and for `−R²/4 < α < 0` the tensor cone lies outside the metric null cone | not re-derived; internally consistent |
| `tangherlini/traversable_throat.py` | for the Tangherlini bridge `G_ret(c,s)·G_adv(c,d) ≡ 0` across exteriors by support composition | argument checked; correct |
| `tangherlini/dynamics.py` | 4+1 Einstein–scalar characteristic evolution, second-order convergent; `A > 0` on any regular-centre ingoing slice | scope statements accurate |

Read together: **a smooth traversable inner–outer throat in 5D Einstein
gravity with the present BAM matter needs `T_kk < 0` at the neck and nothing
supplies it; the constant-EGB escape ends in a classical ill-posedness; and
the horizon branch that needs no exotic matter does not transmit a retarded
signal across exteriors.** The Morris–Thorne–Yurtsever transaction mechanism
of PR #216/#276 is therefore unsupported. This is a real negative result, not
a visualization artefact.

It is *not* fatal to every reading of the premise. PR #206 and
`docs/configuration_space_emergence.md` already take the physical bridge to
be non-traversable, with the pair identification "a property of the
geometry" imprinted at nucleation rather than transmitted afterwards, and the
later mouth-freeze probe converts that into a requirement (the handle must be
dynamically inert at measurement). Whether such a global topological
constraint can reproduce Bell/QFT without importing quantum state rules is
unresolved (that is F2 and F9), but the no-traversability theorems do not
decide it:

```
no traversable causal exchange  ≠  no possible global topological state constraint
```

### F5. The nonminimal-scalar loophole is closed by a modelling assumption, not by the theorem

Candidate C10 is dismissed because at a node `q = 0` the prefactor
`1 − 8πG₅ξq²` equals 1 and the sign is `sign(1 − 2ξ)`. Correct at a node. The
audit then says "BAM places the defect core at `q = 0`, so this is the
relevant point"; that is a modelling choice, not part of the Hochberg–Visser
flare-out theorem.

Away from a node the classical literature goes the other way: Barceló &
Visser, *Traversable wormholes from massless conformally coupled scalar
fields*, Phys. Lett. B 466 (1999) 127, give an exact three-parameter family of
traversable wormholes in Einstein gravity with a non-ghost conformally
coupled scalar, supported where `8πGξq² > 1` and the effective Newton
constant changes sign. Two qualifications: those solutions are 4D, not BAM's
finite 5D `S⁴`-matched handle; and using the branch means abandoning the
identification of the throat core with the `q = 0` order-parameter node and
accepting the effective-gravity-sign region. So it is a legitimate classical
escape branch, not presently a BAM solution. The source audit should list it
as a sixth branch and reword "closes in every dimension" to "closes at a
node".

### F6. The antipodal identification is implemented as an imposed boundary condition; it is not derived

**Correction of the first draft of this audit**, which said the identification
is "never imposed in a solve". That was wrong; the search had been confined to
the `geometrodynamics/` package and missed `experiments/closure_ledger/`.

What exists:

* `null_throat_boundary_conditions_probe.py` (PR #129, 615 lines) takes
  the PR #128 identification `Φ(U,V,Ω) = Φ(−U,−V,Ω̄)` at the 5D Tangherlini
  horizon, uses `Y_ℓ(−Ω) = (−1)^ℓ Y_ℓ(Ω)`, and derives even-ℓ Neumann /
  odd-ℓ Dirichlet throat conditions, zero throat flux (unitary mirror), and a
  real discrete ℓ-graded spectrum on `[R_MID+ε, R_OUTER]`.
* `antipodal_horizon_exchange_kernel_probe.py` (PR #135) builds
  `K_ℓ(r,r′;ω) = ⟨r|(H_ℓ − ω²)⁻¹|r′⟩` with those conditions and checks the
  spectral representation, reciprocity and `(−1)^ℓ` parity grading to ~1e-14.
* The QNM comparison (`antipodal_vs_absorbing_qnm_probe.py`) and one-loop
  self-energy probe build on the same boundary data.

What does not exist:

* a derivation of that identification from the finite-mouth 5D geometry;
* a solve of the nonlinear Einstein system on a global antipodal quotient;
* any twisted spin structure on the quotient (the probe itself notes that a
  Möbius field would flip even ↔ odd and leaves it there).

The package-level uses of `antipode4` remain what the first draft described:
mouth placement at `−p`, a Gaussian match weight with free width
`SIGMA_ANTI = 0.18`, a delayed half-ring coupling in the QCD solver with free
strength `antipodal_coupling = 0.15`, and a cavity resonance rule with
hand-set default frequencies. The ledger entry "B2_antipodal_Z2: CLOSED
(RP³ + spin structure)" still has no computation behind the words "spin
structure".

Accurate status: **BAM has implemented the antipodal identification as a
boundary condition on a linear matter-wave problem on the imposed Tangherlini
horizon geometry. It has not derived that identification from the finite-mouth
solution, nor solved the Einstein system on a global antipodal quotient.**

### F7. The mass ladders are calibrated surrogate fits, and the README contradicts itself about them

`solved_lepton_masses_mev()` returns `[0.511, 105.6126, 1778.938]`
(`−0.043 %`, `+0.117 %` on μ, τ), and the "Lepton mass ladder" section still
advertises "sub-percent accuracy with zero free parameters at scan time". The
same README's PR #271 section states that `γ = 22.5` is "no longer derived",
that the corrected radial operator gives `22.331` or `22.836` at the canonical
geometry, and that either value moves the muon by `15–21 %` because
`d ln m_μ / d ln γ = −17.5`.

The locked Hamiltonian carries four real anchors (`phase = 0.001`,
`transport = 25.1`, `pinhole = 22.5`, `resistance = 0.217869`) and three
discrete structural switches, selected by the calibration scans in the
section's own script map. Raw parameter counting is not the decisive test,
since some anchors may be independently sourced; the decisive test is the
repository's own dependency ledger (derived | chosen | calibrated | holdout),
and on that ledger `γ` moved from *derived* to *chosen* in PR #271 while the
headline stayed *predicted*. The quark sector adds three opt-in extensions and
reaches 1.6 %. Accurate status: **excellent calibrated surrogate fits, not
parameter-free mass predictions.** The row "**The masses are FITTED**"
deeper in the same table already says so.

### F8. The suite catches numerical failures; it has not caught entitlement errors

**Correction of the first draft**, whose heading repeated the prior internal
audit's "the instruments cannot fail". GitHub Actions run 32330296533 on the
PR #263 branch (2026-08-20, head `9cbfa25`) concluded `failure`
(1 failed, 1185 passed, 1 xfailed; a quadrature convergence test). The suite
can say no.

What survives, and is more interesting:

* The failures the suite has caught are implementation/numerical
  (convergence, bookkeeping). The failures that recent reviews caught — the
  radial operator short by `3A²/(4r²)`, an unfixed `η_topo = +1`, a chosen
  tube area carrying an answer, an underidentified fit, and now F1 — passed
  through a green suite. 45 assertions in `tests/` compare a returned string
  against prose (`assert "NARROWED, not closed" in ledger["verdict"]`); they
  lock the narrative, not the number. The `γ` sums, `R_OUTER` fixed point and
  `1.054` factor are unprotected (PR #271's own finding).
* The three pre-registered rounds (`765dbaa → 7de86ce`, `ca07204 → 1082f5d`,
  `d47df40 → 70ad286`) commit the criterion 24–40 minutes before the build,
  by the same session. That is internal prospective testing: it prevents
  "see answer → choose criterion", which is its purpose, and it is not
  independent evidence. The PR #277 sequence did that job correctly.
* The external oracles remain few: Matyjasek 2021 for the 5D QNM, and the
  Hochberg–Visser recovery in the source audit.

### F9. The central gap: the quantization and probability map is not derived

**Correction of the first draft**, which said "no field is quantized
anywhere". The absence of a Fock space or canonical commutators does not
settle that; the repository has a substantial formal path-integral arc:
`Z = Σ_{k odd, c₁, n_part} (−1)^k ∫ (dL/L) det^{−1/2}_matter ·
e^{i(π/2)(1−2a)} · e^{−S_BAM}` with Faddeev–Popov gauge fixing, a
`Diff(S¹)` ghost determinant, ζ / Gel'fand–Yaglom determinants, η phases and
a sector sum (README rows for PRs #74, #115–#122;
`s_bam_path_integral_measure_probe.py`, `diff_s1_ghost_determinant_probe.py`,
`tangherlini_fluctuation_determinant_probe.py`).

So there are two opposite problems, not one:

* some modules use quantum states, tensor products, amplitudes and Born
  probabilities explicitly (`bell/bulk_identity.py`, `history/closure.py`'s
  closure weights, the Bhabha/Møller squared amplitudes fed by the
  exchange-kernel probe);
* other modules build path-integral / one-loop machinery on the classical
  background and describe the result as emergent.

`THESIS.md` holds both positions. Its opening frames `S_BAM`, the one-loop
determinant and the bounded interacting vacuum as matter QFT read off the
classical geometry; lines 149–151 say "No canonical commutators are imposed,
no Hilbert space is assumed at the outset, and no path integral is
performed." Those need reconciling.

The question neither side answers is:

> **Why does classical GR imply the quantum composition and probability
> rules?**

Nothing in the repository derives from the Einstein equations the weight
`e^{iS/ħ}` on histories, the rule `P = |A|²`, or the passage from the
one-object classical configuration space to a complex tensor-product Hilbert
space with the measurement structure Bell requires. Geometry has produced
discrete spectra, holonomies, Green functions, self-adjoint kernels and
algebraic correspondences; the bridge

```
classical geometric phase space  ⟶  quantum amplitudes / probabilities / multiparticle structure
```

is not closed. Accurate status: **the repository contains quantum-formal
structures, but has not derived the quantization/probability map from
classical GR; in the quantitative sectors that map is imported.**

## Dependency structure: the boxes exist, the arrows do not

The first draft made several central ingredients look absent. The corrected
picture is that they exist at different evidentiary levels and have not been
composed. Read top to bottom, each box is real; the arrows between them are
the missing derivations.

*The diagram below is the picture as of `b771b16`; the arrows marked `?` are
annotated with their current status, established on this branch.*

```
  finite classical mouth geometry          (PR #277, Darmois-matched, derived)
            │  ✓ the antipodal BC is the η=+1 sector of the unique free involution
            │    (conditional on P_B = −P_A, quotient over cover, η)          [PR #279]
            ▼
  antipodal boundary condition             (PR #129/#135, imposed postulate; consequences derived)
            │  ✓ the RP² mouth's Pin⁻ holonomy squares to −1 in J's component;
            │    the Hopf fibre IS the mouth's Spin(2), so the K is canonical;
            │    the direction is gauge, the sign is the Pin⁻ sector          [PR #279]
            ▼
  Hopf / quaternionic transport candidate  (J = iσ_y K = L_{−j}; base antipode × fibre reversal; Spin(4) lift (−j, 1))
            │  ? branch aggregation, sector coefficients, readout, composition
            │    all still underived                                          [PR #280]
            ▼
  discrete spectral / Green-function machinery   (compact spectra, Z₂ sectors, holonomies, self-adjoint kernels, η phases: derived)
            │  ? "represented using QM" versus "equivalent to QM" not distinguished
            ▼
  formal quantum / path-integral machinery       (Z = Σ ∫ (dL/L) det^{−1/2} e^{−S}, Hilbert-space Bell: present, assumed)
```

What the repository has demonstrated is that classical geometry naturally
produces several ingredients that quantum theories use: compact spectra, `Z₂`
sectors, holonomies, self-adjoint kernels, `η` phases, topological charge.
Those are not trivial. What it has not demonstrated is that those structures
organise themselves *uniquely* into quantum mechanics rather than merely
permitting a rewrite in quantum-mechanical language. That is a clean
falsification target, and it should be the organising question of the
project.

*Follow-up (same branch):* links 1 and 2 were tested in
`docs/finite_mouth_bundle_transport.md`. Result:
`FINITE_MOUTH_ADMITS_BUT_DOES_NOT_SELECT_THE_BAM_LIFT`. The antipodal
boundary condition is the `η = +1` sector of the unique free involution of
the finite handle, conditional on `P_B = −P_A` and on taking the quotient;
`J = iσ_y K` is the spin lift of the rotation `L_{−j}`, not of `ι`; and the
induced Pin⁻ holonomy of the `RP²` mouth lies in the same component of
`Pin(2) ⊂ SU(2)` with `H̃² = −1` (outcome B of
`docs/mouth_pin_holonomy_prereg.md`); the `U(1)` direction is gauge, the
Hopf fibre is the mouth's `Spin(2)` (so the `K` is canonical), and the sign
is a Pin⁻ sector that pair creation must produce in opposite pairs
(`docs/mouth_spin_frame_prereg.md`). Link 1 has therefore moved from
"assumed" to a derivation conditional on the antipodal quotient
construction. Link 3 (the quantum probability law) was then tested in
`docs/classical_born_rule.md`: with the mouth's classical measure and
deterministic detectors, Born is reached only by importing a Haar-on-`S²`
variable and a tuned weight
(`CURRENT_BAM_PREPARATION_AND_DETECTOR_DYNAMICS_DO_NOT_DERIVE_BORN`, the
pre-registered label narrowed after review); the classical route that
works is a local hidden-variable model and an open derivation route, so
the pair problem falls under Bell's theorem unless the BAM global boundary
problem, solved with both settings as boundary data, makes the ensemble
setting-dependent. That was then computed in
`docs/closure_measurement_dependence.md`: sharp phase closure of the closed
history does make the source ensemble depend on both settings with exact
no-signalling and `S = 2.14 > 2`, but yields a closed-form correlation that
is not `cos γ` under the positive count of closed histories; the quantum
correlation is the same closure set with the two branches summed with
their holonomy (`e^{iΩ/2} = sgn D`), which gives the quantum joint law with
no projectors (`docs/closure_current.md`). Nothing classical in the
repository selects count over oriented sum, and two further inputs — the
relative outcome-sector coefficients and the current-to-frequency readout —
are equally underived: that is where the quantization gap now sits.

### Ranked missing links — current

*Original links 1 and 2 (topology → transport, finite mouth → antipodal
boundary condition) are substantially resolved, conditionally, and are
recorded in the historical list below rather than here.*

1. **Branch aggregation.** Derive whether observed event frequencies are the
   positive count of closed histories, `Σ|D|`, or their oriented sum with the
   closure holonomy, `Σ e^{iΩ/2}|D| = ΣD`. Nothing classical in the
   repository selects between them: the Pin structure supplies the branch
   holonomy as a *label*, the extremal-action condition is unimplemented, and
   its phase-stationarity proxy is analytically disjoint from sharp closure.
   The sharper form: do the Pin/Hopf data make the closure locus an oriented
   current with local coefficients (sign mandatory), or is the physical
   object a measure on histories (positivity forces `|D|`)?
   *Round 7: the stationary-action route is closed. The functional whose
   saddle gives the coarea magnitude (`−cos θ`) is not additive, and the
   additive one (`θ`) has no critical points, so `κ` in `e^{iκS_H}` is a free
   normalisation — a fourth underived input, not a route to the first.*
2. **The relative sector coefficients.** Derive `π_like = π_unlike` from a
   symmetry or chain argument rather than adopting the counting measure. It
   moves the correlation on both branches and the quantum law appears only at
   `r = 1`; no-signalling does not constrain it.
   *Round 7: the full fixed-setting isometry group has order 4 and two sector
   orbits at every angle except `γ = π/2`, where it is order 8 and transitive.
   `r = 1` is forced there and nowhere else — and `π/2` is not a CHSH angle.
   Fibre symmetries cannot help; the sector weights are fibre-blind.*
3. **The readout.** Derive why observed frequencies would be linear in the
   integrated current rather than quadratic in it — the usual classical
   detector response to an amplitude — or given by another functional.
   *Round 7: measured under `φ ↦ cφ`, all five null stresses of
   `source_audit.py` and the conserved mouth current `Im(q*Aq)` have degree
   exactly `2`; none is linear. A quadratic readout keeps marginals at `1/2`
   but gives `E = −2cos γ/(1+cos²γ)` and `S_max = 8√2/3 = 3.7712` —
   superquantum. Linearity is what holds the model at Tsirelson.*
4. **Composition.** Everything above concerns one pair. Derive `H₁ ⊗ H₂` from
   `Γ₁ × Γ₂` for the opposite-Pin-sector pair, rather than assuming it.
5. **Operational no-signalling to the past — now a live failure, not an open
   question.** *Round 7 ran the gate.* The conditioned density is exactly
   antipodally even, so every **odd** source observable is blind
   (`E[x·m|a,b] = 0`, `P(x·m>0) = 1/2` exactly). Every **even** one signals
   (`E[(x·m)²]` spread `0.0103`), and for non-coplanar settings the supports
   are mutually singular — a one-shot signal, total variation `1`. Every
   classical coupling BAM has is degree 2, hence even, so the readable class
   is exactly the one the repository supplies. Either a dynamical reason for
   source inaccessibility must be derived, or `x` must be declared gauge — which
   collides with round 5's use of `x` as a physical source variable. That
   collision is now the sharpest open problem on this strand.
6. **Only then revisit Bell in full**, and **keep the mass and QCD fits
   frozen** until the trunk is resolved.

### Ranked missing links, as written at `b771b16`

1. **Topology → transport.** Derive the transition function of the actual
   finite `RP²`-type mouth and determine whether its natural lift is
   `J = iσ_y K`, instead of inserting that lift afterwards.
   *(Substantially resolved, conditionally — see the current bottom line.)*
2. **Finite mouth → antipodal boundary condition.** Derive PR #129's
   `Φ(U,V,Ω) = Φ(−U,−V,Ω̄)` (or its twisted variant) from the same physical
   finite-mouth geometry and bundle gluing.
   *(Substantially resolved, conditionally — see the current bottom line.)*
3. **Classical global dynamics → quantum composition and statistics.**
   Derive, rather than assume, why alternatives combine as complex amplitudes,
   why `P = |A|²`, why composition is a tensor product, and why the relevant
   measure is the proposed path-integral one rather than a classical
   statistical ensemble. This is the genuine QFT-emergence hurdle.
   *(Narrowed to the four current items above.)*
4. **Only then revisit Bell.** If 1–3 succeed, the singlet and the Born
   calculation stop being imported quantum mechanics and become consequences;
   if they fail, the Bell arc is downstream of an unproved premise.
5. **Keep the mass and QCD fits frozen** until the trunk is resolved; further
   numerical matches do not bear on links 1–3.

## What is sound

* The finite-mouth handle: geometry, matching, and the lapse-independent neck
  cost.
* The source audit's flare-out theorem and its attribution to
  Morris–Thorne / Hochberg–Visser.
* The support-composition argument that the cross-exterior transactional
  product vanishes on the Tangherlini bridge.
* The negative-EGB closure by degeneration of the classical principal symbol.
* The PR #129/#135 boundary-value construction as a self-consistent, unitary,
  parity-graded transport condition, given the postulate.
* The 4+1 characteristic code's positivity identity and stated limits.
* The retraction record, including the internal audit `69fc599`.

## Recommendations

1. Correct `embedding/transport.py` and the README rows that cite it: `σ` is
   the base antipode composed with fibre reversal, orientation-preserving on
   `S³`; `J = iσ_y K` is the quaternionic structure on `C²`. Add a test on
   `det = +1` and on `h(σz) = −h(z)`, `σ(e^{iφ}z) = e^{−iφ}σ(z)` so the
   docstring cannot drift back. Then attempt the actual theorem: that the
   `RP²` mouth transition lifts to `J`.
2. Rewrite the claim table's Bell rows to the repository's own later verdict:
   derived network dynamics gives `CHSH = 2`; `2√2` requires the imported
   singlet and measurement freedom; the carrier's geometry is open.
3. Keep the three inner–outer objects distinct in the README (visual
   convention / imposed horizon BC / finite-mouth junction), and build the
   missing composition: derive the PR #129 boundary condition, or its
   twisted variant, from the PR #277 finite-mouth geometry instead of leaving
   it as an independent postulate.
4. Add Barceló–Visser as a sixth branch in the source audit, with the
   explicit statement that it requires abandoning `q = 0` at the neck and is
   a 4D result.
5. State the traversability no-go at its correct scope: it closes MTY causal
   exchange, not the PR #206 global-constraint reading, which then needs its
   own derivation.
6. Mark the lepton and quark ladders "calibrated" in the top table and
   reconcile the "sub-percent, zero free parameters" section with the PR #271
   section using the dependency ledger.
7. Replace prose-verdict assertions with numeric ones, and regression-lock
   the `γ` sums, `R_OUTER` fixed point and `1.054` factor.
8. Reconcile `THESIS.md` with itself on whether a path integral is performed,
   and state in one place what is derived versus assumed about amplitudes,
   `|A|²`, and the tensor product.

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

Output on this machine: `det = 1.0`, `hits = 0`. The base-antipode and
fibre-reversal identities follow from `z̄₁′z₂′ = z₂(−z̄₁) = −z̄₁z₂` and
`σ(e^{iφ}z) = (e^{−iφ}z̄₂, −e^{−iφ}z̄₁)`.

## Revision note

The first draft of this audit was reviewed against the repository history.
Four of its statements were wrong or overstated, and each correction was
verified before being adopted:

| draft statement | correction | verified by |
|---|---|---|
| "`T = iσ_y` is inserted, not non-orientable geometry" | `σ` is base antipode × fibre reversal; `J` is a natural lift candidate whose derivation from the `RP²` mouth is missing | `h(σz) = −h(z)`, `σ(e^{iφ}z) = e^{−iφ}σ(z)` by hand |
| "nowhere does `x ~ −x` enter a boundary condition of a solve" | PR #129 imposes it at the horizon; PR #135 builds the resolvent on it | `null_throat_boundary_conditions_probe.py` (615 lines), `antipodal_horizon_exchange_kernel_probe.py` |
| "the instruments cannot fail" | CI run 32330296533 failed on 2026-08-20 | GitHub Actions API, `conclusion: failure` |
| "no field is quantized anywhere" | a formal path-integral / determinant arc exists; the gap is the underived quantization map, and `THESIS.md` contradicts itself | README rows for PRs #74, #115–#122; `docs/THESIS.md` lines 149–151 |

Two further scope corrections: the traversability no-go was narrowed from
"fatal to the premise" to "fatal to MTY causal exchange", since PR #206's
non-traversable global-constraint reading is not touched by it; and the
Barceló–Visser branch was qualified as 4D and as incompatible with the
`q = 0` core.

### What has been resolved since the audit was written

The findings F1–F9 above are the audit as written at `b771b16`. Subsequent
rounds on this branch changed the status of three of them; the current
bottom line at the top of this document supersedes the original one.

| original finding | current status | where |
|---|---|---|
| F1 — the `RP²` mouth → `J` theorem is unproved | substantially resolved, conditionally: the mouth's induced Pin⁻ holonomy squares to `−1` and lies in `J`'s component, the Hopf fibre **is** the mouth's `Spin(2)` so the `K` is canonical, and the direction is gauge | `docs/finite_mouth_bundle_transport.md`, `docs/mouth_spin_frame_prereg.md` (PR #279) |
| F1 — the transport docstring is false | fixed in the module; `σ = L_{−j}`, `det = +1`, Spin(4) lift `(−j, 1)` | `geometrodynamics/embedding/transport.py` (PR #279) |
| F6 — the antipodal BC is imposed, not derived from the finite mouth | substantially resolved, conditionally: it is the `η = +1` scalar sector of the unique free involution, at the finite neck, no horizon | same (PR #279) |
| F9 — the quantization map is not derived | narrowed, not closed: sharp closure gives setting-dependent source measure with detector no-signalling and `S = 2.14`, and the holonomy-weighted current gives the quantum joint law with no projectors; what remains is branch aggregation, sector coefficients, readout, and composition | `docs/closure_measurement_dependence.md`, `docs/closure_current.md`, `docs/classical_born_rule.md` (PR #280) |

F2–F5, F7 and F8 stand as written.

## Test suite

Run at commit `b771b16` on Python 3.11 (numpy 2.4.6, scipy 1.17.1, sympy
1.14.0), one `pytest -q --timeout=900` invocation per test file, 64 files:

| passed | failed | errors | xfailed | skipped | wall time |
|---|---|---|---|---|---|
| 2007 | 0 | 0 | 1 | 0 | ~26 min |

The one expected failure is the pre-existing `xfail` in
`tests/test_quark_spectrum.py`. A single monolithic `pytest -q -x` run was
also started and was killed by its own 590 s timeout before finishing; the
per-file run above is the complete result. As stated in F8, a green suite
validates the implementation, not the derivations this audit is about.
