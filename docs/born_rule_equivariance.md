# The Born-rule equivariance test: what measure does the BAM transport preserve? (PR #198)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. This document is item 3 of the
> theorem-shaped program: the deepest refutation vector. It derives the
> right question, proves the two theorems that answer it for the actual
> BAM ψ–Φ–q dynamics, and delimits honestly what is and is not thereby
> established. The companion probe verifies every claim on the running
> dynamics.

## 0. The question, and the stakes

For a classical-field-foundational theory there is essentially one honest
route to the Born rule: **equivariance** (Dürr–Goldstein–Zanghì). If
throats are transported by the BAM dynamics, there is a definite question
with a definite answer:

> Is |ψ|² the measure preserved by the BAM transport flow — as in Bohmian
> mechanics, where |ψ|² is the unique equivariant density — or does the
> ψ–Φ–q dynamics equivariantly transport some **other** functional of the
> field?

The stakes are asymmetric and real. If the equivariant measure is |ψ|²,
BAM acquires its first genuine foothold on the Born rule — the item every
audit has flagged as the deepest import (the `geometrodynamics.bell`
module imports it; the #165 guardrail *forbids* importing it into derived
paths precisely because it was underived). If the equivariant measure is
anything else, BAM is refuted at its foundation: no amount of spectral
success survives predicting wrong measurement statistics.

**The answer, up front:** |ψ|² is equivariant, and it is the *unique*
equivariant density among local functionals of the field — but this is
**not automatic**: it holds because of a specific structural property of
the BAM energy functional (all couplings enter the ψ-equation as *real*
potentials), and the probe demonstrates concretely that deformations
violating that property — including the repo's own *imaginary-time
relaxation flow* — destroy equivariance. The test had teeth; the BAM
functional passes.

## 1. The BAM transport flow

The #179/#180 throat-soliton system descends from one energy functional

```
E[ψ, q] = ∫ [ ½|∇ψ|² + ½κ|∇q|² + ½a₀q² + ¼λq⁴ + V_int ] + W_grav ,
V_int = ±½ g |ψ|² q² ,     ∇²Φ = 4πG (|ψ|² + μ q²) .
```

The repo's solvers run its **imaginary-time gradient flow** (relaxation
to the soliton). The Born-rule question lives in the **real-time
Hamiltonian flow of the same functional**:

```
i ∂ₜψ = −½∇²ψ + [ Φ ± ½ g q² ] ψ ,          (the pilot equation)
∂ₜq   = κ∇²q − (a₀ ∓ g|ψ|²) q − λq³ ,        (the order field, real)
∇²Φ   = 4πG (|ψ|² + μ q²) .                  (the potential, real)
```

(The ± is the repo's coupling-sign convention — #179 vs the #180
hardening; nothing below depends on it.)

**The guidance identification.** A throat is a localized excitation
riding the wave; its transport velocity is the local momentum field of ψ:

```
v = J/ρ = ∇S ,      ψ = √ρ e^{iS} ,   J = Im(ψ*∇ψ) .
```

This is not an extra postulate grafted on: (i) a soliton with phase ramp
e^{iv·x} translates at velocity v (Galilean structure of the kinetic
term); (ii) the exact Ehrenfest relation d⟨x⟩/dt = ⟨∇S⟩_ρ holds for the
full nonlinear equation (verified on the running dynamics in the probe);
(iii) it is the unique velocity field whose current closes the continuity
equation (§2). The BAM transport flow is ẋ = v(x, t).

## 2. Theorem 1 (equivariance): the BAM flow preserves |ψ|² exactly

Write ψ = √ρ e^{iS}. The imaginary part of the pilot equation is

```
∂ₜρ + ∇·(ρ∇S) = 0        — exactly, with no correction terms,
```

**because and only because** every coupling in the pilot equation is a
*real multiplicative* potential:

- the self-consistent gravity Φ — real (the Poisson equation with real
  source |ψ|² + μq²); its nonlinearity (Φ depends on ρ) is irrelevant:
  real potentials move the phase, never the modulus;
- the throat-order coupling ±½gq² — real (q is a real field);
- the kinetic term is the standard −½∇² (the #180 hardening's spectral
  kinetic is exactly this operator).

Hence the density ρ = |ψ|² satisfies the same continuity equation as an
ensemble of throats transported by ẋ = ∇S: **if the throat ensemble is
|ψ|²-distributed at any time, it is |ψ|²-distributed at all times.**
Equivariance — the dBB property, surviving the full nonlinear,
self-gravitating, order-coupled BAM dynamics untouched.

What would have broken it (the refutation edge, made explicit):

1. **Any imaginary/dissipative coupling** iγW(x)ψ adds a source
   R = 2γWρ to the continuity equation: no functional of ρ is then
   equivariant. In particular the repo's own **imaginary-time relaxation
   flow** ∂_τψ = ½∇²ψ − V_eff ψ is *maximally* non-equivariant (it does
   not even conserve the norm — it is renormalized by hand each step).
   The distinction between the relaxation flow (a solver tool) and the
   Hamiltonian flow (the physics) is load-bearing.
2. **Derivative couplings** (q·∇ψ, ∇Φ·∇ψ) would deform the current away
   from ρ∇S.
3. **Non-potential gravity** (Φ entering with a complex or velocity-
   dependent coupling) likewise.

The BAM energy functional contains none of these. The probe verifies the
continuity residual vanishes to integrator precision on the running
ψ–Φ–q evolution, and demonstrates control (1) numerically: switching on
iγWψ produces exactly the predicted residual 2γWρ and a drifting
ensemble.

## 3. Theorem 2 (uniqueness): only |ψ|²

Could BAM have equivariantly transported some *other* functional — |ψ|,
|ψ|⁴, the energy density? For any smooth h(ρ) transported by the same
flow:

```
∂ₜ h(ρ) + ∇·(h(ρ) v) = h′(ρ)[∂ₜρ + v·∇ρ] + h∇·v
                     = [h(ρ) − ρ h′(ρ)] ∇·v .
```

This vanishes identically iff h = ρh′, i.e. **h ∝ ρ** — provided the flow
is compressible (∇·v ≢ 0), which the probe verifies on the actual
dynamics (|∇·v| = O(1); the flow is strongly compressible during soliton
breathing and collision). So among densities that are functions of ρ,
|ψ|² is the **unique** equivariant one; Goldstein–Struyve extend
uniqueness to all local functionals of ψ in the linear theory, cited not
re-proved. The probe demonstrates the failure concretely: ensembles
prepared as √ρ- and ρ²-distributed and transported by the same flow
depart from those functionals immediately (KS distance grows), while the
|ψ|² ensemble stays at sampling noise through the same evolution —
including a two-soliton collision.

## 4. Relaxation: the Born rule as attractor

Equivariance says |ψ|² is preserved *if attained*. The stronger dBB-style
statement (Valentini): non-equilibrium ensembles **relax toward** |ψ|² in
the coarse-grained sense. The probe prepares a uniform (maximally
non-Born) ensemble, transports it with the BAM flow through the
soliton-collision dynamics, and measures the coarse-grained H-function

```
H̄(t) = ∫ f̄ ln( f̄ / ρ̄ ) ,
```

which decreases substantially from its initial value — the subquantum
H-theorem operating on the BAM dynamics. The Born rule is not only a
fixed point of the transport; it is an attractor of coarse-grained
non-equilibrium.

## 5. What is and is not established (honest scope)

**Established (the foothold):**
- |ψ|² is *exactly* equivariant under the real-time BAM ψ–Φ–q transport
  (Theorem 1, structural; verified on the running dynamics);
- it is the *unique* equivariant density among functions of ρ (Theorem
  2, given compressibility — verified); alternatives fail demonstrably;
- non-equilibrium ensembles relax toward it (coarse-grained H-theorem,
  demonstrated);
- the property is *not* automatic: it is equivalent to the reality of
  the BAM potentials, and the probe exhibits the failure under
  deformation. The refutation vector was live and did not fire.

This is precisely the epistemic status of the Born rule in Bohmian
mechanics — no weaker, and no stronger. The `bell`-module import can now
be labeled: *derived at dBB grade, conditional on the following.*

**Not established (the conditions):**
1. **The guidance identification** — throat velocity = ∇S — is motivated
   (Galilean boost, Ehrenfest, uniqueness of the closing current) and is
   the only choice under which any measure is equivariant; but a
   first-principles derivation of throat motion from the 5D geometry
   (the throat as geodesic-following geon vs wave-following soliton) is
   its own program.
2. **Operational statistics.** Translating an equivariant ensemble
   measure into measurement probabilities requires the measurement
   regime to be effectively linear (test throat in an external pilot
   wave; subsystem effective wavefunctions). BAM's pilot equation is
   nonlinear through Φ[ρ] and q[ρ]; in the weak-field test-throat regime
   the pilot equation linearizes and the standard dBB measurement
   analysis applies, but the general nonlinear measurement theory is
   open. (Relatedly: a throat's *own* self-field is not a pilot wave for
   itself in superposition; the operative regime is throat-in-external-ψ.)
3. The demonstration dynamics is the 1D reduction with the same
   structure (standard kinetic + 1D Poisson gravity + real q-coupling);
   the theorem is dimension-blind and applies verbatim to the radial 3D
   #180 system, but the 3D ensemble demonstration is not run here.

## References

- D. Dürr, S. Goldstein, N. Zanghì, *Quantum equilibrium and the origin
  of absolute uncertainty*, J. Stat. Phys. 67 (1992) 843. [Equivariance;
  quantum equilibrium.]
- S. Goldstein, W. Struyve, *On the uniqueness of quantum equilibrium in
  Bohmian mechanics*, J. Stat. Phys. 128 (2007) 1197. [Uniqueness of
  |ψ|² among local functionals.]
- A. Valentini, *Signal-locality, uncertainty, and the subquantum
  H-theorem*, Phys. Lett. A 156 (1991) 5; A. Valentini, H. Westman,
  *Dynamical origin of quantum probabilities*, Proc. R. Soc. A 461
  (2005) 253. [Relaxation to quantum equilibrium.]
- L. de Broglie (1927); D. Bohm, Phys. Rev. 85 (1952) 166, 180. [The
  guidance equation.]
- The BAM dynamics: PR #176–#180 (the ψ–Φ–q functional and soliton);
  PR #191 (real-time methods).

## Reproduce (the machine-checkable demonstrations)

```bash
python -m experiments.closure_ledger.born_rule_equivariance_probe
# Verdict: PSI_SQUARED_IS_THE_UNIQUE_EQUIVARIANT_DENSITY_OF_THE_BAM_TRANSPORT
#          _FLOW_BORN_RULE_AT_DBB_GRADE_CONDITIONAL_ON_GUIDANCE_AND_EQUILIBRIUM
```
