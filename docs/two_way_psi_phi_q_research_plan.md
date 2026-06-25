# Two-way ψ–Φ–q evolution: the self-consistent matter–metric–order system (PR #179)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## Closing the loop

PR #178 coupled the self-gravity solver to the throat-order field **one
way**: the matter density `ρ = |ψ|²` drove `q`, but `q` did not feed back —
it neither gravitated nor acted on the wave. This probe closes the loop into
the full **two-way** self-consistent system of three co-evolving fields: the
matter wave `ψ`, the gravitational potential `Φ`, and the throat-order field
`q`.

The whole coupled system descends from **one energy functional** (so the
coupling is consistent, not hand-wired):

```
E[ψ, q] = ∫ d³r [ ½|∇ψ|²  +  ½κ|∇q|² + ½a₀q² + ¼λq⁴  −  ½g|ψ|²q² ]
        + W_grav[ρ_tot] ,     ρ_tot = |ψ|² + μ q² .
```

Its fixed-mass gradient (imaginary-time) flow is the coupled evolution

```
∂_τ ψ = ½∇²ψ − Φψ + ½g q² ψ      (renormalize ψ to mass M)
∂_τ q = κ∇²q − (a₀ − g|ψ|²) q − λ q³
∇²Φ   = 4πG (|ψ|² + μ q²) .
```

### The four back-reaction channels (all live)

| channel | effect | origin |
|---|---|---|
| `ψ ↔ Φ` | Schrödinger–Newton self-gravity | #176/#177 |
| `ψ → q` | the matter density orders `q` where `ρ > ρ_c = a₀/g` | #178 |
| `q → ψ` | the ordered throat core `+½g q²` **binds the wave** | **NEW** |
| `q → Φ` | the order field's energy density `μ q²` **gravitates** | **NEW** |

The `ψ→q` ordering term and the `q→ψ` binding term share the **same** `g`
(both come from the single `−½g|ψ|²q²` in `E`) — the two-way coupling is
consistent, not two independent knobs. Constants: `a₀=0.20, g=1, λ=1,
κ=0.005` (so `ρ_c = 0.20`), `G=1`.

## What each test shows (the coupled flow actually run, N=300 radial grid)

### T3 — self-consistency

The coupled gradient flow **converges**: the energy decreases monotonically
and plateaus (final-segment `ΔE ≈ −6×10⁻⁴`), the order amplitude `max|q|`
plateaus (last-segment growth `0.3%`, `|q| ≈ 0.49`), and the order-field
stationarity residual drops to `~10⁻⁴`. A self-consistent throat-soliton —
`ψ`, `Φ`, `q` mutually consistent — **exists**. The loop closes.

### T4 — two-way back-reaction

At `M = 3` (super-threshold) the order field nucleates (`max|q| ≈ 0.49`).
Versus the pure Schrödinger–Newton soliton (`q = 0`):

| quantity | pure SN (q=0) | two-way ψ–Φ–q | change |
|---|---:|---:|---:|
| well depth `Φ(0)` | −3.03 | −3.18 | **5.0% deeper** |
| core density `ρ_peak` | — | — | **13.4% denser** |

The ordered throat core traps the matter wave (`+½g q²`), concentrating it,
which strengthens the order — a self-reinforcing two-way loop. Below the
threshold (`M = 1`) the order field relaxes to zero and the two-way system
reduces **exactly** to the pure SN soliton.

### T5 — saturation vs collapse

Two branches, set by `q`'s self-gravity `μ`:

- **sub-critical** (`μ = 0.05`): the quartic `λq⁴` **saturates** the binding
  feedback — `max|q|` plateaus at `0.49`, residual `10⁻⁴`: a **stable**
  self-consistent bound throat-soliton. (Intermediate `μ` simply gives a
  denser, more-bound soliton — the well deepens with `q`'s gravity.)
- **super-critical** (`μ = 2.0`): the positive feedback (more `q` ⟹ deeper
  well ⟹ denser `ψ` ⟹ more `q`) overwhelms the saturation and there is **no
  weak-field fixed point** — the flow diverges (`max|q| → 31`, `Φ(0) → −252`,
  residual `→ 6500`): the onset of strong-field gravitational collapse the
  weak-field scheme cannot resolve.

There is a critical self-gravity separating the stable throat-soliton from
collapse.

### T6 — threshold continuity

Below the ordering threshold (`M = 1`: `ρ_peak ≈ 0.005 < ρ_c = 0.20`) the
order field relaxes to zero and the system is **exactly** the pure
Schrödinger–Newton soliton of #176/#177. Above it (`M = 3`: `ρ_peak ≈ 0.48 >
ρ_c`) the order field nucleates and the self-consistent state is a
throat-soliton. The #176 (self-gravity) → #178 (one-way ordering) → #179
(two-way soliton) arc is one continuous system, switched by the matter
concentration.

## The answer

The throat-order field is **not a passive readout** of the geometry. Closing
the `ψ–Φ–q` loop, it back-reacts both ways: the ordered core traps the matter
wave and gravitates, forming a self-consistent **throat-soliton** whose well
is deeper and core denser than the pure self-gravitating state. The binding
feedback saturates into a stable bound object; `q`'s own gravity, pushed past
a critical coupling, drives runaway collapse.

## Honest scope

- **Weak-field, semi-dynamical, spherically reduced.** The self-gravity
  sphericalizes (#176), so the radial system is the relevant one; the full
  axisymmetric/relativistic two-way problem is not solved here.
- **Effective constants.** `a₀, g, λ, κ, μ` are effective: the result is the
  **structure** — a self-consistent two-way throat-soliton with a stable
  binding branch and a runaway self-gravity branch — not the specific
  numbers (whose microscopic values await `V(q)` and the q–metric coupling
  derived from the 5D bulk).
- The **stable soliton is sub-critical**; the strong-field throat (the
  runaway/collapse endpoint) is for full numerical relativity.

## Reproduce

```bash
python -m experiments.closure_ledger.two_way_psi_phi_q_probe
# Verdict: TWO_WAY_PSI_PHI_Q_CONVERGES_TO_A_SELF_CONSISTENT_THROAT_SOLITON_BACK_REACTION_REAL
```
