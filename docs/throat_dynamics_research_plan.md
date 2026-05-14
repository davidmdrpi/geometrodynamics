# Throat-dynamics research plan

Opens the follow-up to the closure-ledger sequence (PRs #11–18,
summary in `docs/hbar_origin_note.md`). The closure-ledger reduced
the locked lepton surrogate's residual external input from six
phenomenological parameters to one anchor (m_e), and identified the
inner-boundary regularization ε as a closure-quantum invariant
`ε = resistance / k_5⁴ = 7π / (100·5⁴)`. The Compton bridge
`ℏ = m_e R_MID c` closes to 0.04 % at this ε.

The identification is **numerical**, not physical. The hard wall at
`r = r_s + ε` is a mathematical convenience: the 5D Tangherlini
radial operator has a singular point at the throat (`r = r_s`, where
`f(r) = 0`), and truncating the grid at `r_s + ε` is the simplest
way to get a discrete spectrum. The closure-quantum match
`ε ≈ resistance / k_5⁴` is structurally elegant but does not derive
the hard-wall regularization scheme itself.

This thread asks the deeper physics question: **what is the right
inner boundary condition for the radial eigenproblem, derived from
throat dynamics rather than imposed by hand?**

## Statement of the problem

In tortoise coordinates `r* = r + (r_s/2)·ln|(r − r_s)/(r + r_s)|`,
the throat maps to `r* → −∞`. The radial equation

```
−u''(r*) + V(r*) u(r*) = ω² u(r*),     with V(r* → −∞) → 0
```

has *plane-wave asymptotics* at the throat. Two linearly independent
solutions `u ~ exp(±iω r*)` exist at `r* → −∞`, and neither decays.
The throat is an *asymptotically free* boundary, not a confining one.

The current scheme imposes a hard wall at `r = r_s + ε` (Dirichlet
`u = 0` at the corresponding `r*`), which discretizes the spectrum
by confining the wavefunction. As `ε → 0`, the box width
`L = r*(R_OUTER − ε) − r*(r_s + ε)` diverges logarithmically; the
spectrum becomes continuous and the lowest eigenvalue ω(1, 0) → 0.

The closure-ledger result that the surrogate's mass-ratio prediction
is ε-invariant (probe `scale_bridge_regularization_probe`) shows
that the closure-quantum machinery does not depend on the
regularization — but the ω(1, 0) factor that sets the absolute
MeV scale does. Therefore: **the physical content of the BAM
dimensional bridge lives in the inner boundary condition**, and the
closure-ledger framework cannot reach it.

## Candidate routes

Three routes to a physically-motivated inner boundary, in roughly
increasing order of physical commitment:

### (A) Boundary-condition substitution

Replace Dirichlet `u(r_s + ε) = 0` with a different self-adjoint
boundary condition:

  - **Neumann:** `u'(r_s + ε) = 0`  (free end).
  - **Robin:** `u'(r_s + ε) + κ·u(r_s + ε) = 0`  for some κ.
  - **Periodic:** identify the wavefunction across the throat with a
    T-action: `u(r_s + ε) = ±u_antipode(r_s + ε)` for some
    antipodal map. This is the BAM-natural "non-orientable throat"
    BC.

A discrete spectrum emerges for any of these at finite ε. The
question is whether any candidate is **ε-converged** (eigenvalue
stable as ε → 0) and whether it agrees with the closure-quantum
spectrum.

This is the natural first probe — it tests the simplest physical
hypothesis (the throat is a wall with non-trivial reflection
phase, not absolute Dirichlet) without introducing new free
parameters.

### (B) Throat-thickness regularization

Replace the hard wall with a *smooth* confining potential
`V_throat(r)` peaked near the throat. The simplest ansatz: add a
Gaussian or step function

```
V_full(r, l)  =  V_Tangherlini(r, l)  +  V_0 · θ(r_throat − r)
```

where `r_throat = r_s + δ_throat` is a finite throat thickness and
`V_0` is the throat-confinement strength.

In the limit `V_0 → ∞`, this recovers the hard wall. For finite
`V_0`, the spectrum is well-defined without an ε regularization;
the wavefunction has a tunneling tail into the throat region.

The question: do any natural `(δ_throat, V_0)` derive from BAM
ingredients, and do they give the closure-quantum spectrum?

### (C) Quasi-regular Frobenius asymptotics

The throat `r = r_s` is an *irregular* singular point of the radial
equation in r-coordinate (the coefficient `ω²/f(r)` diverges as
`1/(r − r_s)`). Standard Frobenius doesn't give a regular series
solution. However, a *generalized* asymptotic expansion may exist
that selects one of the two oscillating solutions at the throat.

If the right asymptotic behaviour can be characterised (e.g., a
specific reflection phase φ(ω) that the throat geometry imposes),
the discrete spectrum emerges from matching that phase to the
outer Dirichlet BC. This is the most physically natural route but
the hardest to make concrete.

## Sub-targets, in order

### (1) BC substitution probe (concrete)

**Status.** Open. Implement Neumann and Robin BCs at the inner
endpoint and scan κ. Compute:

  - ω(1, 0; R*, ε, κ) for several κ values and ε from 1e-2 down
    to 1e-5.
  - Convergence in ε for each κ.
  - Whether any (κ, ε) combination matches the closure-quantum
    spectrum or closes the Compton bridge.

**Concrete probe.** `experiments/closure_ledger/throat_boundary_condition_probe.py`.

Expected outcome: no simple BC change removes the ε-dependence,
because the throat is asymptotically free and any local BC at
finite ε just shifts the reflection phase. This is a *useful*
negative result — it identifies what kind of physics is required
(non-local matching, or finite throat thickness).

### (2) Throat-thickness probe

**Status (2026-05-14): closed — negative result.**
`experiments/closure_ledger/throat_thickness_probe.py` replaces the
hard wall with a smooth confining sigmoid `V_throat(r) = V_0 / (1 +
exp((r − r_t)/σ))` centered at r_t = r_s + δ_throat, then scans
(V_0, δ_throat) and tests both ε-convergence and Compton-bridge
closure under closure-quantum natural choices.

Outcome:

- **σ is a third arbitrary parameter.** The smoothing scale shifts
  ω by O(σ/δ): default σ = δ/30 puts the V_0 → ∞ limit at ω ≈ 1.04,
  while σ = δ/100 (sharp step) gives ω ≈ 1.004 (matching the
  hard-wall reference). Apparent "Compton-clean hits" at the
  default σ (`V_0 = 100·γ` at +0.13 % from ω = 1) are σ-fitting
  artifacts; at σ → 0 they shift by 1–3 %.
- **No ε-convergence at finite V_0.** Even for strong barriers
  (V_0 = 100, 10⁶) the spectrum still drifts with ε at the
  few-% level. The throat remains asymptotically free; tunneling
  into the throat region depends on the cutoff.
- **No closure-quantum natural (V_0, δ) on the Compton-bridge
  surface at σ → 0.** Best within-5 % at the sharp limit:
  V_0 = ∞ at +1.46 %, V_0 = 100·γ at −1.84 %. Comparable to
  precisions of other closure-quantum identifications, but
  insufficient to claim structural derivation.

The thickness model is a valid mathematical alternative to the
hard wall but not parameter-reducing (3 parameters vs 1) and not
derivable from closure-quantum scaffolding. Sub-target (B) is
closed.

### (3) Quasi-regular reflection-phase analysis

**Status.** Open. Hard analytical work. Compute the *reflection
phase* φ(ω) at the throat for various BCs and see if a natural
phase derives from BAM physics (e.g., from the T = iσ_y action,
or from the Hopf holonomy at the throat).

### (4) R_MID self-consistency (deepest)

**Status.** Open. Beyond this thread. R_MID = 1 by convention. If
the throat thickness is structurally determined by quantum
fluctuations of the throat geometry, then so is R_MID. Closing
this would lift the m_e anchor — fully outside the closure-ledger
scope (THESIS.md "self-consistent throat radius").

## Stopping condition

The thread closes when:

  (a) A physical inner BC is identified that (i) gives a discrete
      spectrum without ε regularization, (ii) yields the
      closure-quantum mass ratios with the same precision as the
      hard-wall scheme, and (iii) closes the Compton bridge
      cleanly. With this, BAM is dimensional-scale-incomplete only
      modulo m_e, with no remaining unphysical regularization.

  (b) A clean negative result is established: no simple BC or
      throat-thickness model reproduces the closure-quantum
      spectrum. The hard-wall scheme is then either (i) the
      correct effective description and the closure-quantum
      ε = `resistance/k_5⁴` is its physical content, or (ii) the
      problem requires the deeper R_MID self-consistency route.

Either outcome sharpens the framework.

### Status (2026-05-14): partial outcome (b)

Sub-targets (1) BC substitution and (2) thickness regularization
both returned negative results. No purely local prescription at
the inner endpoint reproduces the closure-quantum spectrum without
external input. The remaining route within the closure-ledger
scope is sub-target (3) reflection-phase analysis, which is non-
local in the radial coordinate. If (3) also fails, the thread
closes on outcome (b)(ii): the hard-wall scheme with
`ε = resistance/k_5⁴` is the correct effective description, and
the deeper R_MID self-consistency (sub-target 4, THESIS.md) is
the next layer of physics.

## Cross-references

- `docs/hbar_origin_note.md` — closure-ledger paper draft; §5 and
  §6 motivate this thread.
- `docs/hbar_origin_status.md` — probe-by-probe ledger of the
  closure-ledger sequence.
- `experiments/closure_ledger/scale_bridge_regularization_probe.py`
  — establishes the ε-dependence of ω(1, 0) and the Compton-bridge
  recovery at ε* ≈ 3.51×10⁻⁴.
- `experiments/closure_ledger/inner_boundary_derivation_probe.py`
  — identifies ε = resistance/k_5⁴ as the closure-quantum form.
- `experiments/closure_ledger/throat_boundary_condition_probe.py`
  — first probe in this thread (sub-target 1, negative result).
- `experiments/closure_ledger/throat_thickness_probe.py`
  — second probe (sub-target 2, negative result).
