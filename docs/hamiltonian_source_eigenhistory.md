# The Hamiltonian source eigenhistory (PR #219)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. #218 established the eigenhistory
> existence theorem with an *imposed* source phase law. This PR replaces
> it with the minimal conservative source that has an explicit
> Hamiltonian, includes the source's state and energy in the loop —
> **with the time-averaged interaction energy ⟨g·q·u(0)⟩ in the
> ledger** — solves the homogeneous condition U(X)X = X **together
> with** total-energy closure, and reports the fixed-point residual and
> the **full Hamiltonian stability spectrum** (the Duffing (q,p)
> evolved as independent variables, not slaved). The companion probe
> machine-checks every claim (~25 s).

## 0. The source, derived

The minimal conservative source is a side-coupled Duffing oscillator:

```
H = p²/2 + ω₀²q²/2 + (μ/4)q⁴ + g·q·u(0)
```

(the quartic is the minimal anharmonicity — the cubic is absent by
q → −q symmetry). Harmonic balance at frequency ω, with q = Re[a·e^{−iωt}]
and the gauge u(0) = U₀ real:

```
D(a)·a = −g·U₀ ,          D(a) = ω₀² − ω² + (3/4)μa²
```

and the field jump [u′] = g·a makes the source a point scatterer of
**derived** strength κ_eff = g·a/U₀:

```
t_s = 2iω/(2iω − κ) ,     r_s = κ/(2iω − κ)
```

Three things #218 imposed are now consequences:

- **unitarity** |t_s|² + |r_s|² = 1 (to 10⁻¹³ over the working grid) —
  because the Hamiltonian is conservative, the scattering must be;
- **reactivity** — zero net power (g·ω/2)·Im(a·Ū₀) = 0 on the real
  branch (10⁻¹⁵);
- **the amplitude-dependent phase** — arg t_s = −atan(κ/2ω) with the
  amplitude dependence carried by D(a); the **#218 imposed law is its
  weak-coupling limit** (φ → −κ/2ω = g²/(2ωD(I)); measured ratio
  1.000000 at g → 0).

## 1. The ring eigenproblem with the source inside

Because the Hamiltonian source also *reflects*, the loop is a genuine
two-direction ring — [source] —π— [mouth A] —τ— [mouth B] —π— [source] —
with state X = (right-mover, left-mover; source amplitude a). The
homogeneous condition U(X)X = X is the ring monodromy eigenproblem with
the source's S-matrix depending on the state itself.

Machine facts (unitarized Tangherlini ports, spline-interpolated and
re-unitarized pointwise):

- the ring trace is **real** (time-reversal; 10⁻¹⁰);
- **the #217 gap tangency is resolved**: a fine scan finds the tr > 2
  gap segment ([2.7325, 2.7390], max 2.0006) — the ring modes are the
  split pair at the gap edges, exactly where tr crosses 2;
- **the homogeneous condition defines a branch, not a point**: the
  nonlinear mode curve ω*(A) traced from U₀ = 0.2 to 1.6 (eigen-residual
  ≤ 10⁻¹⁴ everywhere, the frequency pulled by the source's
  nonlinearity). Homogeneity alone cannot fix the amplitude.

## 2. The joint solve: U(X)X = X + total-energy closure

The amplitude is fixed by the energy budget. Two equations, two
unknowns:

```
tr T_ring(ω, U₀) = 2                       (the homogeneous condition)
E_field + E_source + ⟨g·q·u(0)⟩ = E₀       (total-energy closure,
                                            corrected ledger)
```

with the time-averaged interaction energy ⟨g·q·u(0)⟩ = g·a·U₀/2
(negative: a and u(0) are antiphased) included in the total
Hamiltonian energy. Solved by 2D Newton at E₀ = 11.177:

| quantity | value |
|---|---|
| ω* | 2.732375 |
| U₀* | 0.914367 (shifted from 0.912168 by the corrected ledger) |
| a* (source amplitude) | −0.1963 |
| Newton residuals | (9×10⁻¹⁶, 1.0×10⁻¹³) |
| **‖U(X*)X* − X*‖/‖X*‖** | **3×10⁻¹⁴** |
| source slaving residual | 10⁻¹³ |
| E_field / E_source / **E_int** | 11.060 / 0.171 / **−0.05397** |
| energy closure residual | 10⁻¹³ |
| raw-port systematic | 2×10⁻⁵ (interpolation, anchored) |

The state X* *includes the source*: its amplitude a*, its internal-state
(frequency-pull) shift, and its energy share are components of the one
closed history.

**Energy closure, element-wise and dynamically**: the conserved flux
form |R|² − |L|² is constant around the ring (2×10⁻¹⁵); the source
absorbs zero net power; the full nonlinear state iterated 10⁴ passes
drifts by 5×10⁻¹¹.

## 3. The full Hamiltonian stability spectrum

Harmonic balance slaves the source algebraically to the instantaneous
field; its winding-map spectrum (retained for comparison) cannot see
the source's own dynamics. The full analysis evolves the Duffing
**(q, p) as independent phase-space variables**: the reduced 2-dof
Hamiltonian

```
H_red = (P² + ω_r²Q²)/2 + (p² + ω₀²q²)/2 + μq⁴/4 + g_eff·q·Q
```

with **ω_r = 2.738858** (the bare ring mode, source removed) and
**g_eff = g·ψ(0) = 0.3203** (ψ the normalized bare mode function) both
*derived* from the ring. Its periodic orbit — the eigenhistory in
reduced form — is found by shooting (residual 4×10⁻¹⁴), with the NNM
frequency landing within **8.8×10⁻⁵** of the full-ring ω* (the
single-mode reduction consistency metric). The **4×4 variational
monodromy** over one period gives the full Hamiltonian Floquet
spectrum:

```
λ = { 1 (double, to 8×10⁻⁹ — the Floquet-trivial pair:
       along-flow + energy) ,
      0.454058 ± 0.890972i  (|λ| = 1 exactly) }
```

- **the nontrivial pair is the source**: its unwrapped rotation
  frequency (θ + 2π)/T = **3.2102** ≈ the Duffing's dressed frequency
  (bare ω₀ = 3.2 plus anharmonic + coupling dressing) — exactly the
  degree of freedom the slaved map froze;
- **symplectic to machine precision**: det M = 1 − 6×10⁻¹⁴,
  MᵀJM − J at 5×10⁻¹⁴; energy drift along the orbit 6×10⁻¹¹;
- **all four on the unit circle** — the conservative eigenhistory is
  marginal in the *full* Hamiltonian sense: Novikov-passive, no runaway
  direction anywhere in its phase space.

The winding-map (slaved) spectrum answers a different question — the
convergence of the winding iteration — and is reported alongside; only
the Floquet monodromy exposes the source pair.

## 4. Honest scope

- Harmonic balance keeps the fundamental; the neglected third-harmonic
  weight at the fixed point is O(μa²/8ω²) ≈ 2×10⁻⁴.
- The corrected ledger's interaction term is ~0.5% of E₀ here; its
  inclusion shifts U₀* by 0.2% — small, but the ledger is now the full
  Hamiltonian.
- The Floquet analysis lives on the reduced 2-dof Hamiltonian (one ring
  mode + source); higher ring modes are dropped, with the NNM-vs-ring
  frequency match (8.8×10⁻⁵) as the consistency metric.
- The Duffing quartic is minimal; any conservative source with a
  nonlinear frequency pull plays the same role.
- E₀ is a specified budget (its physical value — e.g. the #58 nucleation
  quantum — is program-level input). ℏ is still not derived; the scales
  are (ω₀, μ, g, E₀).
- The mouth ports are spline interpolants of the unitarized Tangherlini
  greybody; the raw-port systematic at the fixed point is 2×10⁻⁵.
- Monochromatic (carrier-closed) skeleton; packet eigenhistories need
  the group condition (named successor since #218).
- Classical, zonal scalar, frozen geometry; MTY network history posited.

## 5. What would falsify this

- Source scattering violating unitarity or absorbing net power — the
  Hamiltonian derivation would be wrong. (Checked: 10⁻¹³/10⁻¹⁵.)
- No tr = 2 crossings — no ring modes. (Checked: the gap segment
  resolved, crossings bracketed.)
- A Newton fixed point failing U(X)X = X or the energy budget.
  (Checked: 2.5×10⁻¹⁴ and 5×10⁻¹⁴.)
- A Floquet eigenvalue off the unit circle — growth or decay of the
  conservative history in its full phase space. (Checked: all within
  10⁻⁹; det M = 1 and symplectic to 5×10⁻¹⁴.)
- A missing source pair — the (q,p) dynamics would be spurious.
  (Checked: present, at the dressed frequency 3.2102 vs bare 3.2.)

## 6. Companion probe

`experiments/closure_ledger/hamiltonian_source_eigenhistory_probe.py`
(T1–T8, ~25 s): the derived source checks; the ring modes and the
resolved gap; the joint Newton on the corrected ledger with all
residuals; the element-wise and dynamic energy closure; the full
Hamiltonian Floquet spectrum with the slaved comparison.

**Verdict:**
`THE_SOURCE_PHASE_LAW_IS_DERIVED_FROM_A_HAMILTONIAN_THE_EIGENHISTORY_SOLVES_UXX_EQUALS_X_WITH_TOTAL_ENERGY_CLOSURE_AND_ITS_STABILITY_SPECTRUM_IS_MARGINAL`
