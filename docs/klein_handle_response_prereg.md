# Pre-registration: a localized, responding Klein handle and the pairing of its two sectors

Publish this file before implementing or running anything. Baseline: branch
head `096a78d` (round 12 corrected, `main` at `840183f` merged). Seed
`2026111002`.

**Prior state.** Nothing below has been computed. The predictions are
analytic, derived in this file. As in round 12, freeze and implementation
are written in one session by one author; the only ordering evidence is that
this commit precedes the code.

## Question

Round 12 ([klein_time_attachment.md](klein_time_attachment.md)) left, after
review, one constructive result: a time-reversing identification constrains
complementary quadratures. It was a **global** quotient of a **linear**
field with nothing responding. The review asked for the next calculation:

> whether a localized, responding mouth pair preserves that structure — and
> whether its physical response supplies a pairing between the
> complementary sectors.

## Relation to the boundary-pairing protocol

This is **not** Stage 1 of
[conjugate_boundary_pairing_protocol.md](conjugate_boundary_pairing_protocol.md).
It does not use the supported Maxwell bridge, the protocol's preparation
grid or incoming characteristic data, and it closes none of that protocol's
gates. It is a structural calculation on a lattice network, meant to hand
Stage 1 a symmetry oracle and a located candidate for the pairing variable.
It adopts the protocol's rules where they apply: the mouth coordinate is
solved from the same action as the field and is never prescribed; the
support is an explicit positive oscillator and is labelled
`EFFECTIVE_SUPPORT_ASSUMED`; a symplectic or charge form is not called a
positive norm; nothing is identified with a spinor.

## Model

**Network.** An exterior ring of `n_E` sites with two ports `a ≠ b`, and a
tube of `n_T` sites joining `a` to `b` through a seam. The seam may reverse
time. The time-orientation cover is two copies of the exterior (sheets 1 and
2) with two tubes cross-connected, `θ: a_1 → b_2` and `θ': a_2 → b_1`. In
cover time this is an ordinary static network. The sheet swap `σ` is a
spatial symmetry of it, `σ(θ) = θ'`. Sheet 1 is exterior 1 with tube `θ`.

Frozen sizes: primary `n_E = 12`, `a = 0`, `b = 5`, `n_T = 4` (32 cover
sites); replicate `n_E = 20`, `a = 0`, `b = 7`, `n_T = 6` (52 cover sites).
Every structural check below must hold at both.

**Action.** Field `Ψ = x + iy` on the cover sites, one waist coordinate per
tube (`q` on `θ`, `q'` on `θ'`), static potential `A_i`:

    H = ½|p_x|² + ½|p_y|² + ½ xᵀK x + ½ yᵀK y + Σ_i A_i ρ_i
        + Σ_tubes ( p_q²/2M + K_s q²/2 ),
    K(q,q') = graph Laplacian + m² + V,   V_i = U_0 exp(−2 q h_i) on tube sites,
    ρ_i = x_i p_{y,i} − y_i p_{x,i}.

`V` is the lattice form of the protocol's `λN ℓ(ℓ+1)/r_q²` with
`r_q = r_0 exp(q h)`, `h = (1 − s²)³`. The mouth equation is
`M q̈ = −K_s q + Σ h_i V_i |Ψ_i|²`, a waist-opening load, and its leading
coupling is `q|Ψ|²`. The **neutral** model is `y = p_y = 0`, `A = 0`.

Frozen parameters: `m² = 0.25`, `U_0 = 2`, tube stations `s` equally spaced
in `[−0.6, 0.6]`, `M = 1`, `K_s = 4`. Potential profile `w = 1` on tube
sites and `0` elsewhere, in a `σ`-even version (`+1` on both tubes) and a
`σ`-odd version (`+1` on `θ`, `−1` on `θ'`); `A = α w`.

**Integrator.** Implicit midpoint, step `0.02`, fixed-point tolerance
`1e-13`, mirror time `0`, histories on `[−4, 4]`. It is symplectic,
symmetric, equivariant under linear symmetries and reversing symmetries, and
preserves quadratic invariants.

**Deck maps.** With `P` the permutation matrix of `σ`, signs
`s_t, s_f ∈ {+1, −1}` (`−1` = the seam reverses time, resp. the fibre) and
twist `η`:

    x → η P x,          y → η s_f P y,
    p_x → η s_t P p_x,  p_y → η s_t s_f P p_y,
    (q, q') → (q', q),  (p_q, p_q') → s_t (p_q', p_q).

`s_f = −1` is an antilinear action on `Ψ`: the bundle's transition function
through the handle reverses the fibre.

## Analytic predictions

### H1. The localized handle keeps the Lagrangian half

The deck map scales the symplectic form by `s_t`. For `s_t = −1` it is an
antisymplectic involution of the cover's phase space and, when it preserves
`H`, conjugates the flow to its inverse. Admissible histories are those
whose state at the mirror time is deck-fixed.

Predicted, neutral model, mouths frozen at `q = 0`, both twists, both sizes:

* deck map involutive and antisymplectic; fixed space of exactly half the
  dimension (primary: field `32` of `64`); symplectic form zero on it;
* the cover operator's `σ`-even spectrum equals the spectrum of an
  independently built **ordinary** handle network, and its `σ`-odd spectrum
  that of the same network with one **sign-twisted** tube edge. So
  admissible fields are time-even solutions on the ordinary handle plus
  time-odd solutions on the twisted one (untwisted; reversed for `η = −1`);
* control, `s_t = +1`: the fixed space is a symplectic subspace.

### H2. Locally nothing is halved, and the sheet Wronskian is not conserved

At the mirror time the sheet-1 Cauchy data of admissible fields are
**free**: the restriction map has full rank (primary `32`). For two
admissible histories the sheet-1 Wronskian satisfies

    Ω_1(τ) = Ω_1(−τ)

exactly, and is predicted **not constant**: symplectic flux crosses the
seam. Gate for non-constancy: relative variation `> 1e-6` on `[−4, 4]`. Its
size is a measurement, not a prediction.

### H3. Response preserves the structure

`H` is even in `Ψ` and in the momenta, so the deck map is an exact
reversing symmetry of the coupled nonlinear system. Predicted for the
neutral responding model at amplitudes `A = 1/64, 1/32, 1/16, 1/8` and at a
stress amplitude `1/2`:

* `D Φ_h D = Φ_h^{-1}` on generic states;
* admissible histories satisfy `z(τ) = D z(−τ)` to `1e-9`, independent of
  amplitude;
* the tangent space to the admissible histories is Lagrangian at every
  sampled time (`≤ 1e-9`) and has half the dimension (primary `34` of `68`
  with the two mouths);
* tangent by the linearized scheme agrees with central finite differences
  to `1e-6`;
* total energy is conserved to `1e-5` of the prepared energy for
  `A ≤ 1/8` (the protocol's gate).

### H4. The response is real and has the protocol's order counting

From mirror data with `q = q' = 0`, `p_q = p_q' = 0` and field amplitude
`A`, at fixed time: `q = A² q_2 + O(A⁴)` and
`ψ − ψ_frozen = O(A³)`. Predicted adjacent log-2 slopes over the four
amplitudes: `2 ± 0.15` for `max|q|`, `3 ± 0.3` for the field deviation.
The two mouths must respond differently, `max|q − q'| > 0`, with
`q'(τ) = q(−τ)`.

### H5. Four seam types for a charged field

| `s_t` | `s_f` | deck map | charge density `ρ` | where `Re Ψ`, `Im Ψ` sit |
|---:|---:|---|---|---|
| `+` | `+` | symplectic, linear | even | both in the symplectic fixed space |
| `+` | `−` | symplectic, antilinear | **odd** | spatial Alice handle |
| `−` | `+` | antisymplectic, linear | **odd** | both in `L_+` |
| `−` | `−` | antisymplectic, antilinear | even | `Re` in `L_+`, `Im` in `L_−` |

`ρ` picks up `s_t s_f`. Predicted at both sizes: fixed space of half the
dimension in all four cells (primary `66` of `132`); Lagrangian for
`s_t = −1`, symplectic for `s_t = +1`; total charge identically zero on the
fixed space of the two mixed cells (`≤ 1e-12`).

In the antilinear cells a twist `η` is a rotation of `Ψ` by `π` and is not
a separate label: the two twist sectors of round 12 are the real and
imaginary parts of one field.

### H6. The potential must have the matching parity

The deck map preserves `H` if and only if `A_{σ i} = s_t s_f A_i`.
Predicted in the `(−,−)` cell at `α = 0.5`: with the `σ`-even profile
admissible histories satisfy the seam relation to `1e-9`; with the
`σ`-odd profile the same test fails by at least `1e-4` of the field scale.

### H7. Where the pairing is, and what does not supply it

In the `(−,−)` cell the fixed space has `x, p_y` `σ`-even and `y, p_x`
`σ`-odd, so

    Q = Σ ρ_i = ⟨x, p_y⟩ − ⟨y, p_x⟩

is the round-12 Wronskian pairing between the untwisted sector (the data of
`Re Ψ`) and the twisted sector (the data of `Im Ψ`). Predicted:

* `Q` restricted to the field part of the fixed space is nondegenerate with
  signature `(32, 32)` at the primary size, and equals `Ω(Re, Im)` to
  `1e-12`;
* `Q` is conserved to `1e-9` along responding histories with `A ≠ 0`, and
  the sheet-1 charge satisfies `Q_1(τ) + Q_1(−τ) = Q`;
* **the mechanical response supplies no coupling between the sectors**:
  from admissible data with `Im Ψ = 0` and `A = 0`, `Im Ψ` stays zero
  (`≤ 1e-12`) at every amplitude including the stress amplitude, with the
  mouths responding;
* **the potential does**: with the `σ`-even profile, the generated
  `‖Im Ψ‖ / ‖Re Ψ‖` is first order in `α` — adjacent log-2 slopes
  `1 ± 0.1` over `α = 1e-3, 2e-3, 4e-3`.

### H8. It is a pairing, not a norm and not a complex structure

`Q` is indefinite (H7). Multiplication by `i` does not preserve the
admissible set: predicted `dim(Fix ∩ i·Fix) = 0` on the field part in the
`(−,−)` cell.

## Verdict fields

    numerical_validity         PASS | UNRESOLVED_NUMERICALLY | INCOMPLETE
    support_dynamics           EFFECTIVE_SUPPORT_ASSUMED
    electric_response          EXTERNAL_POTENTIAL_PRESCRIBED
    localized_structure        LAGRANGIAN_HALF_SURVIVES_A_LOCALIZED_RESPONDING_HANDLE
    local_bracket              SHEET_DATA_FREE_AT_MIRROR_WRONSKIAN_NOT_CONSERVED
    neutral_sector_pairing     NOT_SUPPLIED_BY_MECHANICAL_RESPONSE
    charged_sector_pairing     CHARGE_PAIRS_THE_SECTORS_ONLY_WITH_BOTH_REVERSALS
    positive_response          NOT_DERIVED
    spinor_matching            NOT_DERIVED
    event_law                  NOT_DERIVED
    bam_cell                   NOT_ESTABLISHED

Status vocabulary as corrected in round 12: a required check that was not
run is `INCOMPLETE`; a check that ran and missed its tolerance is
`UNRESOLVED_NUMERICALLY`; only an explicit counterexample is
`COUNTEREXAMPLE`. A prediction of H1–H8 that is contradicted by a resolved
measurement is reported as wrong in the first paragraph of the write-up.

## What may not be concluded

* **The cell BAM occupies is not established.** The repository records a
  time-reversing antipodal map (PR #129, PT-type) and a fibre-reversing
  mouth transport (PR #279). If both carried over to a fibre-charged field
  on a handle, that would be the `(−,−)` cell. No matching map from the
  Hopf/spin-frame data to this lattice field is constructed here, and the
  protocol forbids identifying two structures because their descriptions
  agree. `bam_cell` stays `NOT_ESTABLISHED` whatever the numbers are.
* A charged field is **additional field content**. The protocol's Maxwell
  master channel is neutral.
* The potential is prescribed and static. No electric degree of freedom
  responds, and no Wheeler-flux solution is constructed.
* The support oscillator is not GR support dynamics.
* An indefinite conserved pairing is not a probability, a detector law or a
  Hilbert-space inner product.
* A lattice network is not a solved handle geometry. The symmetry
  statements are exact at any size; the scaling slopes are finite-time
  statements about this model.
* Round 12's cutoff lesson applies: nothing here is to be read as a
  continuum uniqueness or existence statement.
