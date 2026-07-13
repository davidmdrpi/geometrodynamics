# The transactional Compton propagator (PR #213)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. The Compton arc (#35/#45/#46, capstone
> #211) reconstructed the standard QED amplitude on the frozen geometry;
> this PR replaces the one element that was still *imported* — the
> propagator's time structure — with a derivation from complete classical
> bulk histories. The companion probe machine-checks every identity
> (~10 s).

## 0. What was imported, until now

The arc derived the propagator's **spatial** structure from geometry:

- the **pole** from the antipodal S³ Green function (#35),
- the **magnitude** 1/q² from the flat limit of the S³ scalar Green
  function (#45),
- the **Lorentz tensor** −η^{μν}/q² from the Hopf-bundle connection (#46).

But its **time** structure — the Feynman *iε* prescription, i.e. the
split of G_F into two time-ordering completions and the assertion that
they add coherently with relative phase +1 — was taken from QED. In the
standard theory that structure is the "virtual particle running forward
and backward in time". A program whose only primitive is classical
geometry cannot leave it as an import: it must come from the one place
the program allows — **complete classical histories on the frozen bulk**
(the Wheeler–Feynman / transactional structure the repo's engine,
`geometrodynamics/transaction/`, already implements as offer →
confirmation → phase closure).

## 1. The construction

Work per mode of frequency ω (conventions: (∂_t² + ω²)G(t) = δ(t);
G_ret(t) = θ(t) sin(ωt)/ω; G_adv(t) = −θ(−t) sin(ωt)/ω).

**Step 1 — time symmetry is not a choice.** The bulk is *frozen*:
static, so t → −t is an isometry. The elementary field of a source is
therefore the time-symmetric Wheeler–Feynman propagator
Ḡ = (G_ret + G_adv)/2 — the geometry cannot prefer retardation over
advancement. (Machine check: on the S³ conformal tower ω_ℓ = (ℓ+1)/R,
the mode-sum kernels satisfy G_ret(t, ψ) = G_adv(−t, ψ) to machine
zero.)

**Step 2 — complete absorption is a geometric fact, not an
assumption.** Wheeler–Feynman theory needs a *perfect absorber* — in
flat space a contingent hypothesis about distant matter. On the closed
S³ it is a **theorem of the geometry**: every retarded wavefront
refocuses at the antipode at t = πR (#166's 1/sin ψ focusing) and
returns to the source at t = 2πR. Nothing escapes; every offer meets an
absorber. Measured on the smeared tower kernel: **50%** of the |G|²
mass concentrates within 0.3 rad of the antipode at t = πR, versus
**10⁻²⁸** at mid-flight t = πR/2, with the peak back at the launch
radius on the return leg.

**Step 3 — the absorber response selects G_F uniquely.** Complete
histories leave no free remnants: the total field must be purely
positive-frequency (e^{−iωt}) in the far future — every emitted quantum
is a confirmed transaction — and purely negative-frequency in the far
past. Adding the general homogeneous response a·cos(ωt) + b·sin(ωt) to
Ḡ, that boundary condition is a 2×2 linear system (condition number
1.0 — well-posed) with the **unique** solution a = i/(2ω), b = 0:

```
G_F(t) = Ḡ(t) + (i/2ω) cos(ωt) = (i/2ω) e^{−iω|t|} .
```

**The Feynman propagator is the time-symmetric field plus the absorber
response of the closed universe.** Equivalently (machine-checked with
linear ε-convergence, mismatch 1.9×10⁻² → 1.9×10⁻⁴ over ε = 10⁻¹ →
10⁻³): G̃_F(Ω) = G̃_ret(Ω) for Ω > 0 and G̃_adv(Ω) for Ω < 0 —
"positive frequencies propagate forward, negative backward" is the
frequency-domain face of the same boundary condition.

## 2. The two completion orderings

G_F(t) splits into its θ(t) and θ(−t) pieces — in transactional terms,
the **offer segment** (the exchanged excitation lives *after* the first
vertex) and the **confirmation segment** (it lives *before* it). Their
half-line Fourier transforms are the two **old-fashioned
perturbation-theory energy denominators**:

```
I₊(Δ) = −(1/2ω) · 1/(Δ − ω + iε)        (offer ordering)
I₋(Δ) = +(1/2ω) · 1/(Δ + ω − iε)        (confirmation ordering)
```

Three measured facts about them:

1. **Individually regulator-dependent.** The finite-T truncation of
   either ordering never converges — its late-time spread stays at the
   full oscillation amplitude (0.50 measured vs 0.49 predicted) — while
   the ε-damped ordering converges to the closed form to 7×10⁻¹².
   ε is not decoration: it is the complete-absorption condition itself,
   the inverse absorption time of the closed universe.
2. **Individually non-covariant.** Neither ordering is a function of
   Δ²: the Δ → −Δ evenness violation of a single ordering is ≥ 0.21
   (relative) across the test grid.
3. **Coherently covariant.** The sum is *exactly* (finite-ε closed
   form, verified to 7×10⁻¹⁵ over a random (Δ, ω) grid):

```
I₊ + I₋ = −(ω − iε) / [ ω (Δ² − (ω − iε)²) ]  →  −1/(Δ² − ω² + iε) ,
```

   and exactly even in Δ → −Δ (machine zero). The covariant pole exists
   *only* as the coherent pair.

## 3. The relative phase is forced

Why do the two orderings enter with relative phase +1? Because the
t → −t isometry of the frozen bulk maps the offer segment into the
confirmation segment (G_ret(t) = G_adv(−t), §1 step 1): they are the
*same history read in the two time directions*, and enter every
complete history with equal weight.

**The deform test (the refutation edge).** Replace the coherent sum by
S_φ = I₊ + e^{iφ} I₋:

| φ | evenness violation | pole-form deviation |
|---|---:|---:|
| 0 | 3×10⁻¹⁵ | 7×10⁻¹⁶ |
| 0.3 | 0.81 | 0.59 |
| π/2 | 1.89 | 2.79 |
| π | 2.00 | 3.94 |

Any φ ≠ 0 destroys the Δ²-dependence (Lorentz/crossing structure) and
the pole form at O(1). Only the coherent point survives.

The repo's transactional engine already enforces this discretely:
`retro_phase_match` (the phase-closure score of every confirmed
handshake) is maximal exactly at total phase 0 and decreases
monotonically away from it. Its second maximum at the π branch is the
**#48 Möbius exchange sign** — the discrete topological sector
(fermionic exchange), distinct from the continuous propagator phase
derived here, and recorded honestly in the probe.

## 4. The Compton tie-in: the denominators, derived

For Compton scattering off an electron at rest (photon ω_γ, scattered
ω′ at angle θ), the two-ordering coherent sums reproduce the covariant
denominators the arc assumed, *exactly* across the kinematics grid
(4 energies × 7 angles, errors ≤ 10⁻¹⁴):

```
s-channel:  (1/2Eₙ)[ 1/(m+ω_γ−Eₙ) − 1/(m+ω_γ+Eₙ) ] = 1/(s − m²) ,
            Eₙ = √(ω_γ² + m²)
u-channel:  (1/2Eₙ′)[ 1/(m−ω′−Eₙ′) − 1/(m−ω′+Eₙ′) ] = 1/(u − m²) ,
            Eₙ′ = √(ω′² + m²)
```

— the u-channel identity holding entirely off-shell (u < m² always),
which is precisely the "virtual" regime: no on-shell intermediate ever
exists, only the coherent pair of ordering segments of one complete
history. The deformed phase (φ = π/2) fails the s-channel denominator
by ≥ 7% everywhere on the grid.

So the Compton amplitude's structure now factors as: **denominators**
(this PR: complete histories) × **numerators/vertices** (#37–#44, #46)
× **pole/magnitude/tensor geometry** (#35/#45/#46). What QED postulates
as the iε prescription is, on this geometry, a theorem about complete
histories.

## 5. Honest scope

- **Derived**: the time structure — the Feynman iε, the two ordering
  completions, their coherent relative phase — from *static* (time
  symmetry) + *closed* (complete absorption) + the complete-history
  boundary condition.
- **Not rederived**: spinor numerators and vertex factors (#37–#44,
  #46) — this PR supplies the denominators those results multiply.
- **The ε → 0 limit** is the complete-absorption idealization. On the
  closed bulk ε is a physical inverse absorption time, not a formal
  regulator — but the limit is still taken.
- **Per-mode, tree-level**: the free propagator only; no loop or
  interacting-history statement.
- **Frozen geometry**: no backreaction of the exchanged field on the
  bulk (the program's standing framing).
- **The π branch**: the engine's phase closure also admits total phase
  π — the #48 exchange sign, a discrete sector this derivation does not
  fix, and does not need to.

## 6. What would falsify this

- A static bulk whose tower kernels fail G_ret(t) = G_adv(−t) — the
  time-symmetry step would collapse. (Checked: machine zero.)
- A closed geometry that does *not* return all radiation — complete
  absorption would revert to an assumption. (Checked: 50% antipodal
  concentration at t = πR, full return at 2πR.)
- A second solution of the complete-absorption system — G_F would not
  be selected uniquely. (Checked: condition number 1.0, unique.)
- Any φ ≠ 0 that still reproduces the covariant pole and the Compton
  denominators — the phase would not be forced. (Checked: O(1)
  violations at every deformed φ.)

## 7. Companion probe

`experiments/closure_ledger/compton_transactional_propagator_probe.py`
(T1–T8; runs in ~10 s) machine-checks: the jump conditions and the
absorber-response decomposition G_F − Ḡ = (i/2ω)cos(ωt) (5.6×10⁻¹⁷,
homogeneous to 5×10⁻⁸); the 2×2 uniqueness (exact); the
frequency-splitting ε-scaling (linear, ratio 10.0/decade); the
truncation instability and the exact covariant-pole identity
(7×10⁻¹⁵); the tower time-reversal identity and the antipodal
refocus/return; the deform table; the engine's phase-closure maximum;
and the s/u-channel OFPT identities (10⁻¹⁴) with the deformed-phase
failure.

**Verdict:**
`FEYNMAN_PROPAGATOR_FROM_COMPLETE_HISTORIES_THE_TWO_ORDERINGS_ARE_OFFER_AND_CONFIRMATION_THE_COHERENT_PHASE_FORCED_BY_THE_FROZEN_BULK`
