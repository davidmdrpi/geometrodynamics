# The Hamiltonian source eigenhistory (PR #219)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. #218 established the eigenhistory
> existence theorem with an *imposed* source phase law. This PR replaces
> it with the minimal conservative source that has an explicit
> Hamiltonian, includes the source's state and energy in the loop,
> solves the homogeneous condition U(X)X = X **together with**
> total-energy closure, and reports the fixed-point residual and the
> stability eigenvalues. The companion probe machine-checks every claim
> (~20 s).

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
tr T_ring(ω, U₀) = 2          (the homogeneous condition)
E_field + E_source = E₀       (total-energy closure)
```

solved by 2D Newton. At E₀ = 11.177:

| quantity | value |
|---|---|
| ω* | 2.732375 |
| U₀* | 0.912168 |
| a* (source amplitude) | −0.1963 |
| Newton residuals | (1.3×10⁻¹⁵, 4.6×10⁻¹⁴) |
| **‖U(X*)X* − X*‖/‖X*‖** | **2.5×10⁻¹⁴** |
| source slaving residual | 10⁻¹³ |
| E_field / E_source | 11.006 / 0.171 |
| energy closure residual | 5×10⁻¹⁴ |
| raw-port systematic | 2×10⁻⁵ (interpolation, anchored) |

The state X* *includes the source*: its amplitude a*, its internal-state
(frequency-pull) shift, and its energy share are components of the one
closed history.

**Energy closure, element-wise and dynamically**: the conserved flux
form |R|² − |L|² is constant around the ring (2×10⁻¹⁵); the source
absorbs zero net power; the full nonlinear state iterated 10⁴ passes
drifts by 5×10⁻¹¹.

## 3. The stability eigenvalues

The Jacobian of the loop map X → T_ring(ω*, U₀(X))·X at the fixed point
(4-dim real):

```
λ = { 1 (double — a Jordan block, numerically split by O(√ε_FD),
       reciprocal product 1 + 3×10⁻¹⁰) ,
      0.9999882 ± 0.0048669i  (|λ| = 1 to 10⁻¹⁰) }
```

**All four on the unit circle** — the conservative eigenhistory is
marginal, as a Hamiltonian fixed point must be: Novikov-passive, no
runaway. The Jordan pair at 1 is the parabolic direction along the mode
branch: perturbations **shear at most linearly** (measured log-log
growth slope 0.999 — exactly linear; late-time ratio 1.82 ≈ 2 for a
doubling of pass count), never exponentially. The rotating pair
(angle 0.00487 rad/pass) is the amplitude–phase libration about the
center. Selection of the eigenhistory against its neighbors remains by
dephasing/shear — a weakly dissipative registration mechanism (the #209
opens) would convert marginal into attracting.

## 4. Honest scope

- Harmonic balance keeps the fundamental; the neglected third-harmonic
  weight at the fixed point is O(μa²/8ω²) ≈ 2×10⁻⁴.
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
- A stability eigenvalue off the unit circle — growth or decay of the
  conservative history. (Checked: all within 5×10⁻⁷, the deviation
  being finite-difference splitting of the Jordan pair with reciprocal
  product 1 + 3×10⁻¹⁰.)
- Exponential perturbation growth — a runaway time loop. (Checked:
  linear shear, slope 0.999.)

## 6. Companion probe

`experiments/closure_ledger/hamiltonian_source_eigenhistory_probe.py`
(T1–T8, ~20 s): the derived source checks; the ring modes and the
resolved gap; the joint Newton with all residuals; the element-wise and
dynamic energy closure; the stability spectrum and shear analysis.

**Verdict:**
`THE_SOURCE_PHASE_LAW_IS_DERIVED_FROM_A_HAMILTONIAN_THE_EIGENHISTORY_SOLVES_UXX_EQUALS_X_WITH_TOTAL_ENERGY_CLOSURE_AND_ITS_STABILITY_SPECTRUM_IS_MARGINAL`
