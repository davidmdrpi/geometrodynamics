# Pre-registration: a responding electric field on the Klein handle

Publish this file before implementing or running anything. Baseline: branch
head `fe35525` (round 13). Seed `2026111003`.

**Prior state.** Nothing below has been computed. The predictions are
analytic, derived in this file. Freeze and implementation are written in one
session by one author; the only ordering evidence is that this commit
precedes the code.

## Question

Round 13 ([klein_handle_response.md](klein_handle_response.md)) found that
on a seam reversing both time and fibre, the conserved charge is the pairing
between the two quadrature sectors, and that the variable coupled to it is
the electric potential. That potential was **prescribed**, so no response
was shown to supply the pairing. This round makes the electric field a
dynamical variable solved from the same action and asks:

1. Does the structure survive electric response?
2. Does electric response to the field, by itself, couple the sectors?
3. If not, what does, and which electric variable carries the pairing?

As before this is **not** Stage 1 of
[conjugate_boundary_pairing_protocol.md](conjugate_boundary_pairing_protocol.md).

## Model

The cover network of round 13 (`Network` in
`bulk/klein_handle_response.py`), same sizes: primary 32 sites / 34 edges,
replicate 52 sites / 54 edges. Edges are oriented: ring edges forward, tube
edges from the `a` end to the `b` end. The sheet swap `σ` maps oriented
edges to oriented edges.

**One Hamiltonian**, lattice scalar electrodynamics in temporal gauge, with
the waist coordinates of round 13:

    H = ½ Σ_i |Π_i|² + ½ Σ_{e: i→j} |Ψ_j − e^{i g a_e} Ψ_i|²
        + ½ Σ_i (m² + V_i(q)) |Ψ_i|² + ½ Σ_e E_e²
        + Σ_tubes (p_q²/2M + K_s q²/2).

`a_e` is the link potential and `E_e` its conjugate electric field:
`ȧ_e = E_e`, `Ė_e = g J_e`, `J_e = Im( conj(e^{i g a_e} Ψ_i) Ψ_j )`. The
Gauss function `G_i = (div E)_i − g ρ_i` is conserved, and physical states
have `G = 0`. Nothing is prescribed. Parameters as in round 13, `g = 1`.
Integration by implicit midpoint, step `0.02`, histories on `[−4, 4]`.

On a graph the gauge field has no radiation. Its free content is one mode
per independent cycle. The cover has three; the handle network itself has
two, the ring and the **handle loop**. Electric flux around the handle loop
is source-free flux trapped by the topology: Wheeler's charge without
charge.

**Deck map.** As in round 13 for the field and the mouths, and

    a_{σe} → s_f a_e,      E_{σe} → s_t s_f E_e.

Fibre reversal is charge conjugation and flips `a`; time reversal flips `a`
again relative to `E`. With a dynamical gauge field there is no parity
condition left to impose: all four seam types are symmetries of the one `H`.

## Analytic predictions

### E1. The structure survives electric response

The deck map is an involution, scales the symplectic form by `s_t`,
preserves `H`, and sends `G` to `s_t s_f G∘σ`. For `s_t = −1` it reverses
the flow. Predicted, both time-reversing cells, both sizes, field
amplitudes `1/64, 1/32, 1/16, 1/8`, generic complex admissible data on the
Gauss surface with field, gauge field and mouths all responding:

* `H(Dz) = H(z)` to `1e-12` relative on generic states, all four cells;
* admissible histories satisfy `z(τ) = D z(−τ)` to `1e-9`;
* the fixed space is half the phase space (primary `100` of `200`) and the
  tangent space along a responding history is Lagrangian (finite-difference
  tangents, tolerance `1e-7`);
* `G` stays zero to `1e-9`, total charge is zero to `1e-12`, energy closes
  to `1e-5`;
* control, `s_t = +1`: the fixed set is invariant under the flow.

### E2. Trapped handle flux is admissible exactly when the charge is untwisted

Source-free electric configurations form the cycle space, dimension `3` on
the cover: `σ`-even part of dimension `2` (ring and handle loop), `σ`-odd
part of dimension `1` (the two rings opposed, zero on the tubes).
Admissibility at the mirror time requires `E∘σ = s_t s_f E`. Predicted:

| `s_t` | `s_f` | admissible flux dimension | flux through the tubes |
|---:|---:|---:|---|
| `+` | `+` | 2 | allowed |
| `+` | `−` | 1 | **forbidden** |
| `−` | `+` | 1 | **forbidden** |
| `−` | `−` | 2 | allowed |

The cells that admit trapped handle flux are exactly the cells where round
13 found a nonzero charge.

### E3. Electric response to the field alone does not couple the sectors

Charge conjugation `C: (x, y, p_x, p_y, a, E) → (x, −y, p_x, −p_y, −a, −E)`
is a symplectic symmetry of `H`. Its fixed set — real field, no gauge
field — is invariant. Predicted, in the `(−,−)` cell and on the ordinary
handle, at every amplitude including `1/2`, with gauge field and mouths
responding: from real admissible data with `a = 0`, `E = 0`, the imaginary
part, the link potentials and the electric field all stay zero to `1e-12`.

### E4. Trapped flux couples the sectors, and its sign sets the orientation

With handle-loop flux `Φ` and real field data, `ȧ = E` winds the holonomy
around the handle loop and the hopping phases rotate `Re` into `Im`.
Predicted in the `(−,−)` cell at amplitude `1/8`:

* generated `‖Im Ψ‖ / ‖Re Ψ‖` is first order in `Φ`: adjacent log-2 slopes
  `1 ± 0.1` over `Φ = 1e-3, 2e-3, 4e-3`;
* exactly odd: the history with `−Φ` has the opposite imaginary part and
  the same real part, to `1e-9`. The sign of the trapped flux fixes the
  sense of the `Re → Im` rotation.

### E5. The flux responds

`Ė = gJ` with `J` bilinear in the field, so the handle-loop flux is not
conserved once charged matter is present. With `ΔΦ(t)` the change in the
handle-loop component of `E`, predicted for `max_t |ΔΦ|`:

* second order in field amplitude at fixed `Φ = 4e-3`: slopes `2 ± 0.15`
  over the four amplitudes;
* first order in `Φ` at amplitude `1/8`: slopes `1 ± 0.1`.

The sign is **not** frozen. Screening (`ΔΦ·Φ < 0`) is expected for a
gapped field; whatever sign is measured is reported.

### E6. The sheet value of the pairing is the flux through the seam

Summing Gauss's law over sheet 1, whose boundary is the two seam edges
(the last edge of each tube):

    g Q_1 = E_{θ, last} − E_{θ', last}.

Predicted along the E4 histories: this identity holds to `1e-9`; total
charge is zero; `Q_1(τ) = −Q_1(−τ)`; and `Q_1` is **not** identically zero
(charge separates across the seam), gate `max|Q_1| > 1e-9`. Round 13's
sheet pairing `Ω_1(Re, Im)` is therefore the difference of the electric
fields on the two seam edges, divided by `g`.

### E7. Control: the ordinary handle

Frozen as a control because round 13's corresponding control was added only
afterwards. On the ordinary handle `(+,+)` trapped flux also generates an
imaginary part at first order (`1 ± 0.1`) and with exact oddness. The
coupling of `Re` to `Im` by flux is therefore **not** specific to the
time-reversing seam.

## Verdict fields

    numerical_validity         PASS | UNRESOLVED_NUMERICALLY | INCOMPLETE
    support_dynamics           EFFECTIVE_SUPPORT_ASSUMED
    electric_dynamics          LATTICE_GAUGE_FIELD_SOLVED_FROM_THE_ACTION
    structure                  LAGRANGIAN_HALF_SURVIVES_ELECTRIC_RESPONSE
    flux_admissibility         HANDLE_FLUX_ADMISSIBLE_IFF_CHARGE_UNTWISTED
    spontaneous_pairing        NOT_GENERATED_BY_RESPONSE_WITHOUT_FLUX
    flux_pairing               TRAPPED_FLUX_COUPLES_SECTORS_AND_SETS_ORIENTATION
    flux_dynamics              FLUX_RESPONDS_AT_SECOND_ORDER
    pairing_carrier            SHEET_PAIRING_IS_THE_SEAM_FLUX_BY_GAUSS
    seam_specificity           FLUX_COUPLING_ALSO_ON_THE_ORDINARY_HANDLE
    flux_origin                INITIAL_DATA_NOT_DERIVED
    positive_response          NOT_DERIVED
    spinor_matching            NOT_DERIVED
    event_law                  NOT_DERIVED
    bam_cell                   NOT_ESTABLISHED

Status vocabulary as in rounds 12–13: `INCOMPLETE` for a check not run,
`UNRESOLVED_NUMERICALLY` for a missed tolerance, `COUNTEREXAMPLE` only for
an explicit one. A prediction contradicted by a resolved measurement is
reported as wrong in the first paragraph of the write-up.

## What may not be concluded

* **The trapped flux is initial data.** Nothing here derives it. E3 says
  response does not create it. A GR origin for the flux, or for its sign,
  is not addressed.
* **The reduced phase space is not analysed.** Statements are made on the
  unreduced space with `G = 0` imposed. States fixed by the deck map only
  up to a gauge transformation are not classified.
* A graph gauge field has no radiation; this is not Maxwell theory on a
  handle geometry, and not the protocol's radiative master channel.
* An orientation for the `Re → Im` rotation is not a complex structure on
  the admissible set, and an indefinite pairing is not a norm.
* A charged field remains additional field content; `bam_cell` stays
  `NOT_ESTABLISHED` whatever the numbers are; the matching map is not
  attempted.
* Scaling slopes are finite-time statements about this lattice model.
