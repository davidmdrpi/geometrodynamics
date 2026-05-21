# Bulk-boundary interaction probe — research plan

Targets the **B5′ residual** left by the radial reduction bridge
(PR #50): a single master integral unifying the bulk radial modes
(masses) and the boundary throat-pinch (the F² vertex). PR #50
factorized the 5D→4D reduction into three channels (radial → masses,
S³ → gauge+propagator, throat → F²) and found, honestly, that F² is
the throat-channel form factor, *not* a radial overlap. The residual:
treat the bulk radial modes and the boundary throat-pinch on the same
footing.

This probe formulates the **bulk-boundary interaction** that does so
for the radial+throat channels: the same throat cavity produces both
the bulk mass spectrum and the throat K factor, via one bulk Green
function and its boundary behavior.

## The bulk-boundary structure

The radial channel is a cavity `r ∈ [R_MID, R_OUTER]` with Dirichlet
hard walls (the throat BC from B3). Define the bulk Green function

```
G(r, r'; ω) = Σ_n u_n(r) u_n(r') / (ω² − ω_n²)
```

where `u_n` are the Dirichlet radial modes and `ω_n = ω(l,n)` the
eigenfrequencies. This single object carries **two** kinds of data:

  - **Bulk spectrum (poles).** `G` has poles at `ω = ω_n` — the
    radial mass ladder. Residues = mode products `u_n(r)u_n(r')`.

  - **Boundary coupling (normal derivatives).** At the Dirichlet
    throat `u_n(R_MID) = 0`, but the normal derivative
    `u_n'(R_MID) ≠ 0` is the throat coupling of each bulk mode. The
    throat-to-throat response is

    ```
    Π(ω) = Σ_n [u_n'(R_MID)]² / (ω² − ω_n²)
    ```

    — the boundary-to-boundary propagator, poles at the masses.

## The throat impedance and the K factor

At the Dirichlet throat, a wave of frequency `ω` dwells for the
equal-action time `τ(ω) = π/ω` (the closure-quantum half-split, PR #41).
This dwell time is the **throat impedance** `Z(ω) = π/ω`. The Compton
in/out photons (frequencies `ω`, `ω'`) traverse the throat with
impedances `Z(ω)`, `Z(ω')` **in series**; the effective rate is the
harmonic mean (PR #39's series-impedance result, now grounded in the
bulk-boundary cavity):

```
K(x) = 2·(1/(Z(ω)+Z(ω'))) · π = 2ω'/(ω+ω') = 2x/(1+x)
```

So the **same throat cavity** yields:

  - the **mass spectrum** as the poles of `G` (bulk), and
  - the **K factor** as the series of throat impedances (boundary),

both from one bulk-boundary structure. This is the master integral for
the radial+throat channels — the unification the B5′ residual asked for,
restricted to those two channels.

## What this unifies vs leaves open

**Unifies (radial + throat):** masses (bulk poles) and the K factor
(boundary impedance) are two faces of one throat-cavity Green function.

**Leaves open:** the **Q factor** (Hopf-fibre helicity) is the S³
angular channel, not the radial cavity — it is not part of the
bulk-boundary throat structure. The full vertex `F² = K²·Q` therefore
still combines the bulk-boundary (radial+throat → K) with the S³ channel
(→ Q). The B5′ residual is **narrowed** from "three separate channels"
to "radial+throat unified by the bulk-boundary cavity; S³ (Q) combined
separately."

## Tests

  T1. **Bulk Green function poles = masses** (P1): `G(ω)` poles at the
      radial eigenfrequencies `ω(l,n)`; residues = mode products.
  T2. **Boundary normal derivatives** (P2): `u_n'(R_MID) ≠ 0` are the
      throat couplings of the bulk modes (Dirichlet ⟹ `u_n = 0` but
      `u_n' ≠ 0`).
  T3. **Throat-to-throat response Π(ω)** (P2): `Π(ω) = Σ [u_n'(R_MID)]²
      /(ω²−ω_n²)`; poles at the masses.
  T4. **Throat impedance → K** (P3): `Z(ω) = π/ω` (dwell time); series
      of in/out impedances → `K(x) = 2x/(1+x)`.
  T5. **One cavity, two outputs**: the same throat cavity produces the
      mass spectrum (bulk poles) and the K factor (boundary impedance).
  T6. **Shared substrate**: `R_MID`, the hard-wall BC (B3), and the
      closure quantum (B1) are shared by the bulk spectrum and the
      boundary impedance.
  T7. **B5′ assessment**: the bulk-boundary interaction unifies the
      radial (mass) and throat (K) channels in one structure; the
      residual is narrowed to the S³ (Q) combination.

## Verdict structure

  - **BULK_BOUNDARY_FORMULATED**: the bulk-boundary interaction is
    formulated; the same throat-cavity Green function yields the mass
    spectrum (poles) and the throat impedance (boundary) → K. The
    radial and throat channels are unified; the residual is narrowed to
    combining with the S³ (Q) channel. B5′ advanced.

  - **FORMULATION_INCOMPLETE**: a piece (poles, boundary response, or
    the impedance → K link) does not hold.

## What this leaves open

  - **S³ (Q) combination.** The Hopf-helicity Q factor is the angular
    channel; combining it with the bulk-boundary (K + masses) to write
    the *complete* `F²` and the mass spectrum in one integral is the
    remaining piece of B5′. The three-channel factorization (PR #50)
    already shows they are reductions of one action; a single covariant
    integral over (radial, throat, S³) together is the final step.

  - **B4 — dimensional bridge.** Unaffected; the single `m_e` anchor
    remains.

## Cross-references

  - PR #50: `radial_reduction_bridge_probe` — the three-channel
    factorization and the B5′ residual.
  - PR #39: `two_mouth_flux_action_probe` — the series-impedance →
    harmonic mean (K), now grounded in the bulk-boundary cavity.
  - PR #41: `throat_action_derivation_probe` — the dwell time
    `τ = π/ω` (the throat impedance).
  - B3 hard-wall derivation — the Dirichlet throat BC.
  - `geometrodynamics/tangherlini/radial.py` — the radial cavity modes.
  - `experiments/closure_ledger/bulk_boundary_interaction_probe.py` —
    this probe.
