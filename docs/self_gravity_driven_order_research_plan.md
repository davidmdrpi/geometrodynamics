# Self-gravity-driven throat-order instability: bind, or drive a new order parameter? (PR #178)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## The question

PR #176/#177 established that weak-field self-gravity concentrates a wave
packet above a critical mass (`M_c ∝ 1/G`) — it **binds**. PR #178-prev
introduced the throat-order field `q` (a Ginzburg–Landau order parameter
whose defects are the throats) but coupled it to the geometry only by hand —
the focusing was *named* as the trigger, not solved together. The standing
follow-up:

> does the self-gravitating **concentration** merely produce a bound lump
> (the field stays in the disordered `q = 0` phase), or does it actually
> **drive** the order parameter — nucleating geometric order?

This probe couples the two and answers it.

## The coupling

The local matter concentration `ρ = |ψ|²` (from the #176/#177 self-gravity
solver, **actually run** here) becomes the **control field** of a
density-dependent Landau potential

```
V(q; ρ) = ½ (a₀ − g ρ) |q|² + (λ/4) |q|⁴ ,     a(ρ) = a₀ − g ρ .
```

The order field's effective mass² `a(ρ)` changes sign at the **critical
concentration** `ρ_c = a₀/g`:

| regime | `a(ρ)` | order field | meaning |
|---|---:|---|---|
| `ρ < ρ_c` (dilute) | `> 0` | `q = 0` only minimum | disordered — **merely bound** |
| `ρ > ρ_c` (concentrated) | `< 0` | `\|q\| = √((gρ−a₀)/λ)` | ordered — **order nucleates** |

The constants used: `a₀ = 0.30, g = 1.0, λ = 1.0`, so `ρ_c = 0.30`. The
order field's amplitude relaxes under the Ginzburg–Landau gradient flow
`∂_t f = −[(a₀ − gρ(r)) f + λ f³] + κ ∇²f` (`κ = 0.02`).

## What each test shows (the self-gravity solver actually run)

Measured peak (spherically-averaged) densities of the #176/#177 collapse:

| regime | M | G | `ρ_peak` | vs `ρ_c = 0.30` | `max |q|` | outcome |
|---|---:|---:|---:|---|---:|---|
| sub-threshold | 1 | 1 | 0.06 | below | ~10⁻¹⁶ | **merely bound** (no order) |
| super-threshold | 3 | 1 | 0.90 | above | 0.68 | **order nucleates** |
| gravity off | 3 | 0 | 0.18 | below | ~10⁻¹⁶ | **no order** (it is gravity) |

- **T3 — merely bound.** A sub-threshold packet (M = 1) concentrates only to
  `ρ_peak ≈ 0.06 < ρ_c`; with `ρ` below the critical concentration
  everywhere, the order field stays at zero. The wave is bound (or
  dispersing), but the geometry carries **no** order.
- **T4 — drives order.** Above the mass threshold (M = 3) the collapse drives
  `ρ_peak ≈ 0.90 > ρ_c`, and the order field **nucleates** a localized
  symmetry-broken domain (`max |q| ≈ 0.68`) sitting **at the density peak**
  (the throat core of #178). The gravitational collapse drives the
  throat-order parameter off zero.
- **T5 — gravitational.** With gravity off (`G = 0`) the same high-mass
  packet never concentrates past `ρ_c` (`ρ_peak ≈ 0.18`), so no order
  nucleates; restoring gravity it does. The ordering is driven by the
  gravitational concentration — not by the bare mass or amplitude — and
  inherits the `M_c ∝ 1/G` gravity of #176/#177.
- **T6 — dynamical.** Driving `q` by the *time-dependent* peak density
  `ρ_peak(t)` of the collapse, the order parameter switches on (causally)
  only **after** `ρ_peak(t)` crosses `ρ_c` — a moving order front following
  the gravitational concentration. The geometric order is not pre-existing;
  it is switched on by the collapse.

## The answer

Weak-field concentration does **not** merely bind. Above a critical
concentration — reached only by the gravitational collapse — it **drives**
the throat-order parameter off zero and nucleates geometric order. The
binding threshold of #176/#177 and the ordering threshold are linked: the
collapse is what carries the matter density across the order transition.

## Honest scope

- **One-way coupling.** The self-gravitating `ρ(t)` drives `q`, but `q`'s
  back-reaction on the metric is **not** yet included — the fully
  self-consistent q–metric system (`q` sourcing the geometry it nucleates in)
  is the next step.
- **Effective constants.** `a₀, g, λ` — and hence the numerical `ρ_c` — are
  effective, chosen to place the order transition inside the collapse's
  reachable density range. The physics result is the **existence** of a
  gravitationally-crossed concentration threshold separating bind-only from
  drive-order, **not** the specific `ρ_c` (its microscopic value awaits
  `V(q)` derived from the 5D bulk).
- **Droplet-size barrier.** The spatial nucleation carries the usual
  Ginzburg–Landau droplet barrier — the ordered core must exceed the
  coherence length `ξ = √(κ/|a|)` — so the spatial ordering threshold sits
  somewhat *above* the local `a(ρ) = 0` crossing; the sub-/super-threshold
  cases (M = 1 vs M = 3) are chosen well-separated so the bind-vs-drive split
  is robust.
- Still **weak-field / semi-dynamical** (not full numerical relativity).

## Reproduce

```bash
python -m experiments.closure_ledger.self_gravity_driven_order_probe
# Verdict: SELF_GRAVITY_DRIVES_THROAT_ORDER_ABOVE_A_CRITICAL_CONCENTRATION_NOT_MERELY_BINDING
```
