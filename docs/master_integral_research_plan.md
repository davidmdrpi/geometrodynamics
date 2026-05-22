# Master integral — completing B5′ with the S³ Hopf Q-channel

Closes the **B5′ residual** left open by the bulk-boundary interaction
probe (PR #51): integrate the **S³ Hopf Q-channel** into the
bulk-boundary master functional, so that a *single* functional on the
internal geometry produces all three channels together —

  - the **mass spectrum** (radial × S³-Casimir poles),
  - the **K factor** (throat dwell-time impedance), and
  - the **Q factor** (S³ Hopf-fibre helicity),

i.e. the mass ladder **and** the full vertex `F²(x,c) = K(x)²·Q(x,c)`
from one object.

## Where B5′ stands

| step | what it did | residual |
|---|---|---|
| PR #50 radial reduction bridge | factorized 5D→4D into 3 channels; F² is the throat form factor, not a radial overlap | a single master integral producing masses **and** F² |
| PR #51 bulk-boundary interaction | unified **radial (masses) + throat (K)** via one throat-cavity Green function | the **S³ (Q)** channel still combined separately |
| **this probe** | integrate the **S³ Hopf Q-channel** into the master functional | — (B5′ closed; B4 the only survivor) |

## The internal geometry is a (warped) product

The 5D bulk internal space is

```
M_int  =  C × S³ ,    C = radial cavity [R_MID, R_OUTER]
```

a **warped product**: the S³ radius warps with `r`, which is exactly why
the S³ Casimir `l(l+…)` enters the radial potential `V_tangherlini(r,l)`
as the centrifugal barrier. Fields separate,

```
Ψ(x^μ, r, Ω) = Σ_{l,n} ψ_{l,n}(x^μ) · u_{l,n}(r) · 𝒴_l(Ω) ,
```

so the internal Helmholtz Green function **separates** into a sum of
factor products. **This is the mechanism**: the F² = K²·Q factorization
is not a failure to unify — it is the direct consequence of the product
internal geometry. One separable kernel, three reductions.

## The master functional

Define the single bulk-boundary-angular master functional

```
ℳ(ω; x, c) = G_C(r, r′; ω) ⊗ 𝒢_{S³}(Ω, Ω′)
```

read off three ways:

  1. **Poles in ω** → the **mass spectrum** `ω(l,n)`. The radial factor
     `G_C(r,r′;ω) = Σ_n u_{l,n}(r)u_{l,n}(r′)/(ω²−ω²(l,n))` has poles at
     the eigenfrequencies; the S³ Casimir `l` enters via the warp. The
     masses are the **product spectrum** (radial ladder × S³ Casimir).

  2. **Throat boundary of `G_C`** → the **K factor**. At the Dirichlet
     throat `u_{l,n}(R_MID)=0` but `u′_{l,n}(R_MID)≠0`; the dwell-time
     impedance `Z(ω)=π/ω` of the in/out photons in series gives the
     harmonic mean `K(x)=2x/(1+x)` (PR #51).

  3. **S³ Hopf reduction of `𝒢_{S³}`** → the **Q factor**. The Hopf
     fibre carries the helicity spinor `(A_pres, A_flip)` with
     `A_pres = x`, `A_flip = √x(1−x)/√(1+c²)`; the `(1+c²)` is the
     Hopf-fibre helicity sum `cos⁴(θ/2)+sin⁴(θ/2)`. Then
     `Q = A_pres²+A_flip² = x²+x(1−x)²/(1+c²)` (PR #40).

The **vertex residue** of `ℳ` — the throat boundary (→ K) dressed by the
S³ Hopf reduction (→ Q) — is

```
F²(x,c) = K(x)² · Q(x,c) ,
```

and the **poles** of the *same* `ℳ` are the mass spectrum. Masses **and**
F² from one functional — the master integral B5′ asked for.

## Shared substrate across all three channels

| substrate | radial (mass) | throat (K) | S³ (Q) |
|---|---|---|---|
| `R_MID` | cavity inner wall | throat radius (dwell) | S³ radius scale |
| closure quantum `2π` | mode normalization | dwell time `τ=π/ω` | Hopf holonomy `π cos χ` |
| `T²=−I` | Dirichlet hard wall (B3) | throat node `u(R_MID)=0` | helicity-flip ε (A_flip) |

The same three pieces appear in every channel — the bridge that makes
them one functional, not three.

## Tests

  T1. **Separable master kernel.** The internal modes separate as
      `u_{l,n}(r)·𝒴_l(Ω)`; the master Green function is the sum of
      factor products. Its ω-poles are the radial eigenfrequencies
      independent of the angular factor.
  T2. **Radial × S³ Casimir → mass spectrum.** `ω(l,n)` is monotone in
      both `l` (S³ Casimir) and `n` (radial ladder) — the product
      spectrum of the warped product.
  T3. **Throat boundary → K.** `Z(ω)=π/ω` in series → `K(x)=2x/(1+x)`
      (the bulk-boundary throat channel, PR #51).
  T4. **S³ Hopf reduction → Q.** `(A_pres,A_flip)` → `Q(x,c)` to
      machine precision; the `(1+c²)` is the Hopf-fibre helicity sum.
  T5. **Master vertex residue = F² = K²·Q.** The single functional's
      vertex (throat × S³ Hopf) reproduces the closed-form F²(x,c) to
      machine precision over an (x,c) grid.
  T6. **One functional, masses AND F².** From the *same* `ℳ`: read the
      mass spectrum (poles) and the full F² vertex (residue). Both, one
      object — the B5′ master integral.
  T7. **Product geometry is the mechanism.** The warped-product
      separation of variables is *why* F²=K²·Q factorizes; verify the
      product spectrum identity (joint mode = radial × angular).
  T8. **Shared substrate.** `R_MID`, `2π`, `T²=−I` appear in all three
      channels (table above).
  T9. **B5′ assessment.** Radial + throat + S³ unified in one master
      functional; B5′ closed; B4 (m_e anchor) the only survivor.

## Verdict structure

  - **MASTER_INTEGRAL_COMPLETE**: the S³ Hopf Q-channel is integrated
    into the bulk-boundary master functional. One separable functional
    `ℳ` on the warped-product internal geometry `C × S³` yields the mass
    spectrum (poles), the K factor (throat boundary impedance), and the
    Q factor (S³ Hopf helicity); its vertex residue reproduces
    `F²=K²·Q` and its poles give the masses — masses **and** F² from one
    object. The factorization is the product-geometry consequence. B5′
    closed.

  - **INTEGRATION_INCOMPLETE**: a channel does not reduce from the
    master functional, or the vertex residue does not reproduce F².

## What this leaves open

  - **B4 — dimensional bridge.** Unchanged: the single `m_e` anchor
    (`ℏ = m_e·R_MID·c`) sets the absolute MeV scale. The master integral
    is dimensionless (`F²`, masses-in-units-of-R_MID); B4 is orthogonal.
  - **First-principles internal action.** The master functional uses the
    established channel reductions (radial Sturm–Liouville, throat
    dwell-time, Hopf helicity); deriving all three from one written-out
    5D Lagrangian density with the throat boundary term is the natural
    follow-on. The factorization shows they are reductions of one
    structure; the explicit covariant Lagrangian is not re-derived here.

## Cross-references

  - PR #50: `radial_reduction_bridge_probe` — three-channel factorization.
  - PR #51: `bulk_boundary_interaction_probe` — radial+throat unified.
  - PR #40: `hopf_helicity_transport_probe` — Q from Hopf helicity.
  - PR #39: `two_mouth_flux_action_probe` — K from series impedance.
  - `docs/bam_scaffold_status.md` — barrier ledger (B1–B5).
  - `geometrodynamics/tangherlini/radial.py` — radial cavity modes.
  - `geometrodynamics/hopf/connection.py` — Hopf connection/holonomy.
  - `experiments/closure_ledger/master_integral_probe.py` — this probe.
