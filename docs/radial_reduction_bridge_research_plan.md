# 5D → 4D radial reduction bridge probe — research plan

Step (3) of the scaffold closure programme — the largest remaining gap.
PR #49 (topological/discrete sector) promoted B1 + B2 to action data;
the hard-wall derivation closed B3. The BAM effective-action scaffold
is down to two barriers:

  - **B4 — dimensional bridge** (the single `m_e` anchor;
    `ℏ = m_e·R_MID·c`); long-standing, documented in
    `docs/hbar_origin_status.md`.
  - **B5 — 5D → 4D reduction** producing F². "The radial spectrum (the
    mass sub-thread) and the F² vertex (the amplitude sub-thread) live
    in separate sub-threads; the reduction map connecting them is
    unconstructed."

This probe targets B5: build the Kaluza–Klein-like reduction of the 5D
Tangherlini bulk action and determine **what the reduction produces and
what it does not** — connecting the mass and amplitude sub-threads, and
honestly identifying any residual that prevents a single master
integral.

## The reduction

Write the 5D field on `M₅ = (4D spacetime) × (radial channel) × (S³)`:

```
Ψ(x^μ, r, Ω) = Σ_{l,n} ψ_{l,n}(x^μ) · u_{l,n}(r) · Y_l(Ω)
```

and integrate the 5D action over the internal geometry (`r` and `Ω`).
The reduction **factorizes into three channels**:

| channel | integrate over | produces | thread |
|---|---|---|---|
| **radial** | `r ∈ [R_MID, R_OUTER]` | KK masses `ω(l,n)` | closure-ledger (mass) |
| **S³ angular** | `Ω ∈ S³` | gauge coupling `c₁ = 1` + propagator `1/q²` | tree-QED (gauge/prop) |
| **throat** | `r → R_MID` pinch | form factor `F²(x, c)` | tree-QED (vertex) |

All three are reductions of **one** 5D action on the **same** internal
geometry, sharing `R_MID` (throat radius), the closure quantum `2π`,
and the spin structure `T² = −I`.

## The key structural question

Does the reduction produce `F²(x, c)` as a **radial overlap integral**
`∫ u_m u_n u_p u_q dr` (the naive KK expectation)? 

**No** — and this is the central honest finding. Radial overlap
integrals are *constants* (numbers), independent of the external
scattering kinematics. But `F²(x, c)` is a non-trivial function of
`x = ω'/ω` and `c = cos θ`. So `F²` **cannot** be a radial overlap; it
is the **throat-channel form factor** (the throat-pinch dynamics from
PRs #39–#41 + the S³ propagator from PRs #45–#46), a 4D-effective
vertex shape, not a KK constant.

The reduction therefore connects the two sub-threads **structurally**
(shared substrate, one 5D action, three factorized channels) but does
**not** unify masses and `F²` into a single master integral: they are
produced by different channels of the same reduction.

## Predictions

### P1. Radial modes form an orthonormal KK basis

The Sturm–Liouville radial modes `u_{l,n}(r*)` (Dirichlet BCs, from B3)
are orthonormal; they are the KK reduction basis.

### P2. Radial channel → mass spectrum

The KK masses are the radial eigenfrequencies `ω(l,n)` — the
closure-ledger lepton/quark ladder.

### P3. S³ channel → gauge coupling + propagator

The S³ reduction gives the Hopf charge `c₁ = 1` (gauge coupling) and the
S³ Green function (propagator `1/q²` in the flat limit, PRs #45–#46).

### P4. Throat channel → F² form factor

The throat-pinch (at `R_MID`) gives `F²(x, c) = K(x)²·Q(x, c)` — the
throat-channel form factor, distinct from the radial modes.

### P5. F² is NOT a radial overlap (the central finding)

Radial overlaps are kinematics-independent constants; `F²(x, c)` varies
with `(x, c)`. Hence `F²` is the throat form factor, not a KK overlap.
The naive "F² from radial integration" is falsified.

### P6. Shared-substrate consistency

`R_MID`, the closure quantum `2π`, and `T² = −I` appear identically in
all three channels — the bridge that connects the sub-threads.

## Tests

  T1. **KK basis** (P1): build a clean symmetric-FD radial operator;
      verify orthonormal modes; cross-check the Chebyshev solver
      spectrum.
  T2. **Radial → masses** (P2): the KK masses are `ω(l,n)`; spectrum
      for `l = 1, 3, 5`.
  T3. **S³ → gauge + propagator** (P3): `c₁ = 1` (`compute_c1`); S³
      Green function flat limit → `1/q²`.
  T4. **Throat → F²** (P4): `F²(x, c) = K²·Q`.
  T5. **F² ≠ radial overlap** (P5): radial overlaps are constants;
      `F²(x, c)` varies with kinematics. The honest falsification of
      the naive reduction expectation.
  T6. **Shared substrate** (P6): `R_MID`, `2π`, `T² = −I` consistent
      across all three channels.
  T7. **B5 assessment**: the reduction factorizes into three consistent
      channels (mass / gauge+propagator / F²) sharing the substrate.
      Residual: a single master integral producing masses AND F²
      together is not written. B5 reduced from "unconstructed" to
      "factorized framework + one residual."

## Verdict structure

  - **BRIDGE_FACTORIZED** (expected): the 5D → 4D reduction factorizes
    into three consistent channels — radial → masses, S³ → gauge +
    propagator, throat → F² — all reductions of one 5D action on the
    shared internal geometry. The mass and amplitude sub-threads are
    structurally connected. The honest residual: `F²` is the
    throat-channel form factor, *not* a radial overlap, so a single
    master integral unifying masses and F² is not achieved. B5 is
    substantially reduced (framework built; residual narrowed to the
    master integral), not fully closed.

  - **BRIDGE_FULL** (not expected): a single overlap integral produces
    both masses and `F²`. Would require `F²` to be a radial overlap,
    contradicting its kinematic dependence.

  - **BRIDGE_FAILS**: the channels do not share the substrate or a
    channel does not reduce as claimed.

## What this leaves open

  - **The master integral (B5 residual).** A single covariant
    reduction integral producing the mass spectrum *and* the F² vertex
    shape together. The three-channel factorization shows they are all
    reductions of one action, but writing them as one integral requires
    treating the throat-pinch (boundary) dynamics and the bulk radial
    modes on the same footing — open.

  - **B4 — dimensional bridge.** Unaffected; the single `m_e` anchor
    remains.

## Cross-references

  - PR #49: `topological_discrete_sector_probe` (B1+B2), and the
    hard-wall derivation (B3) — the spin structure setting the radial
    Dirichlet BC.
  - PRs #39–#41: throat action → `K`, `Q`, `F²`.
  - PRs #45–#46: S³ Green function → `1/q²` propagator + Lorentz tensor.
  - `geometrodynamics/tangherlini/radial.py` — radial modes.
  - `geometrodynamics/hopf/chern.py` — `c₁ = 1`.
  - `experiments/closure_ledger/radial_reduction_bridge_probe.py` —
    this probe.
