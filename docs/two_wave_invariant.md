# The two-wave invariant is branch-resolved

`geometrodynamics/waves/two_wave.py` · renderer
`scripts/geometrodynamics_v59_two_wave.py` · probe
`experiments/closure_ledger/two_wave_probe.py`

## 0. What this round is for

PR #258 built a *static* source-interaction kernel and said plainly that it was
not the two-wave invariant: it carried no local null momenta, so it could not
tell equal-energy collinear from counterpropagating waves — the one control the
two-wave invariant exists to apply. This round builds the time-dependent object
and applies that control.

**The point is not to re-derive the WKB identity.** That

```
𝒞 = A_A²A_B²(k_A·k_B)² ,     k_A·k_B = −E_AE_B(1 − cos θ)
```

is zero collinear and `−2E_AE_B` head-on is known. It is used here as the
*control*. The research content is the **difference** between the exact
time-dependent, multipath, throat-coupled field and that limit — how big it is,
what it is made of, and where it is not small.

The short answer: the departure that matters is not a curvature correction. It
is **multipath**, and it is `O(1)`.

## 1. What is solved

The retarded field of a pulsed point source on the throated ESU, exactly, by
Krein's resolvent formula in the frequency domain:

```
φ_s(x,ω) = s(ω)[ G(χ_xs,ω) + Σ_ij G(χ_xi,ω) R_ij(ω) G(χ_js,ω) ] ,
R(ω) = (A − Γ(ω²))⁻¹
```

inverted along the **retarded contour** `Im ω = ε`. That contour is exact rather
than approximate: with `ω = u + iε`,

```
φ(t) = e^{εt} · (1/2π) ∫ du e^{−iut} φ̂(u + iε)
```

so `ε` only trades pole sharpness against a growth factor, and both are
reported.

**Derivatives without a mesh.** Every term is a function of one geodesic
distance, so the four-gradient and Hessian close in form:

```
∇_a χ = −(y − (x·y)x)/sin χ ,        ∇_a∇_b χ = cot χ (δ_ab − ∇_aχ ∇_bχ)
```

both checked against the exponential map. Nothing is finite-differenced —
time derivatives are `−iω` inside the contour integral, radial ones are the
analytic `∂_χG`, `∂²_χG`.

**The stress tensor** is the improved conformal one,

```
T_{μν} = ∂_μφ∂_νφ − ½η_{μν}(∂φ)² + ξ[η_{μν}□ − ∇_μ∇_ν + G_{μν}]φ² ,  ξ = 1/6
```

with `G_{00} = 3`, `G_{ab} = −δ_ab` for `S³ × R`.

## 2. The solver earns its results

| | |
| --- | ---: |
| free field vs PR #254's closed-form image sum, 4 arrivals | **`3.3e-16`** |
| Maslov signs `+ − + −` | reproduced |
| conformal wave equation `∂_t²φ = ∇²φ − φ`, free | `4e-16` relative |
| … with the throat | `4e-16` relative |
| trace of the improved stress tensor | **`1.9e-15`** relative |

The trace is a real test and not an algebraic one. `□φ` is taken **from the
solve** rather than replaced by its on-shell value `ξRφ = φ`; substituting on
shell makes `T^μ_μ` vanish identically for any input. Computed honestly it
equals `φ(□φ − φ)`, so a nonzero trace *is* a wave-equation residual.

## 3. The known limit, as a limit

Two sources and an observer on one great circle, so the arriving spatial
directions are exactly parallel or exactly antiparallel — set by geometry, to
`1e-12`, not fitted. The diagnostic is the **pointwise** ratio

```
𝒩 = (T_A:T_B)/(T_A^{00} T_B^{00})       WKB:  (1 − n̂_A·n̂_B)²
```

pointwise because for two WKB waves every envelope and every `sin²` factor
cancels between numerator and denominator. A window average does *not* have that
property — it drags `⟨sin⁴⟩/⟨sin²⟩²` and the envelope overlap into the answer,
which is how a first attempt got `8.3` where the answer is `4`.

| `ω₀` | `𝒩` head-on (WKB 4) | `𝒩` collinear (WKB 0) |
| ---: | ---: | ---: |
| `6` | `3.95291` | `1.994e-04` |
| `12` | `3.98171` | `2.860e-05` |
| `24` | `3.99817` | `1.905e-07` |
| `48` | **`3.99995`** | **`1.789e-10`** |

**The collinear null is much stronger than a leading-order statement.** On this
geometry the two arriving wavefronts share their normal *exactly*, so amplitude
gradients are along the same `n̂` and cannot tilt either `k`. The residue falls
faster than any fixed power and the measured exponent steepens with `ω₀`. That
matters for the next section: it means the multipath effect is not competing
with a comparable curvature correction.

**Convergence is part of the measurement.** The contour needs `ε` well above the
frequency spacing `2π/span`. The same run at `ε ≈ du` returns `−6.7e-04` where
the converged value is `1.4e-08` — four orders wrong, and entirely plausible
looking.

## 4. The result: multipath destroys the collinear null

Sources fixed, observation point fixed, carrier fixed. **Only which branch has
arrived changes.**

| branch pair | `𝒩` exact | geometry `(1 − n̂_A·n̂_B)²` |
| --- | ---: | ---: |
| `A` direct + `B` direct | **`1.905e-07`** | `0` |
| `A` **long-way winding image** + `B` direct | **`3.99806`** | `4` |
| `A` direct + `B` **through the throat** | **`0.56501`** | `0.56367` |

Three different values of the two-wave invariant for the same pair of sources at
the same event.

* the **winding image** propagates the other way round the sphere — its phase is
  `t + χ`, so the front moves toward the source and its arrival direction is
  *reversed*. A collinear pair reads head-on;
* the **throat** leg emerges from a mouth, so its direction is set by the
  mouth's position. The exact field agrees with that prediction to **0.24%**,
  and the prediction comes from positions alone — it was never fitted.

The free-propagation control at the same instant is reported with it, because
there `B` has **no arrival at all**: energy product `4.1e-29` against `1.2e-02`.
The throat is *creating* the second branch, not bending an existing one, so the
comparison is stated as amplitudes and not as a ratio there.

**So the collinear null is not spoiled by curvature corrections, which are
`1e-7` here. It is spoiled by multipath, at `O(1)`.** The invariant has to carry
the branch index PR #255 named; a single-branch WKB formula is not merely
approximate on this background, it is answering a different question.

## 5. The other corrections, quantified

### Arrivals

The free arrivals land at `χ`, `2π−χ`, `2π+χ` to `1.3e-03` with signs `+ − +`,
out of a solve that never saw PR #253's ledger. The throat adds arrivals the
free ledger does not contain, at the two-leg times `χ(y,c_j) + χ(c_i,x)` PR #255
named — checked at the **causal onset** rather than at a peak, because `R(ω)`
has poles and a throat arrival rings up instead of pulsing. Onset measured
`2.8839` against `2.8445` predicted, inside one pulse width.

### Tail — there isn't one, except the throat's

`S³ × R` is conformally flat, so the conformally coupled massless scalar obeys
**Huygens' principle exactly**: it propagates strictly *on* the light cones.

| | between geometric arrivals, relative to peak |
| --- | ---: |
| free field | **`1.4e-08`** (the Gaussian source's own wing) |
| with the throat | **`8.1e-02`** |
| amplification | `5.7e+06` |

Every bit of tail in this model belongs to the throat.

### Caustic — where geometric optics stops, with a scale

Geometric optics gives amplitude `1/(4π sin χ)`, divergent at the antipode. The
exact kernel is finite there, `G(π,ω) = ω/(4π sin πω)`, and **linear in `ω`**.
In between, everything depends on the single combination `ωe` with `e = π − χ`:

```
exact / WKB  =  |sin(ωe) / sin(πω)|
```

Measured as a **collapse** — at half-integer carriers the ratio is exactly
`|sin(ωe)|`, identical across `ω = 8.5, 16.5, 32.5` to **`6.6e-15`**. The
caustic is cut off at `e* ∼ 1/ω`, and the saturation amplitude is `ω/4π`.

## 6. And the round closes back on PR #258

`∫dt φ(t) = φ̂(ω = 0)`, so the DC content of the *solved time series* is exactly
the static kernel PR #258 did its tomography on. Running that protocol on
numbers produced by the dynamic solver — least squares for `S = Re R`, then
`𝒲 = S₁₂/det S − G₀`:

| | |
| --- | ---: |
| worst kernel error vs the exact static value | `1.5e-06` |
| worst error in the recovered `S` | `2.8e-04` |
| `𝒲` from the time-dependent solve | **`−0.060010`** |
| `−β` | `−0.06` |

The route goes through the whole contour integral. On the retarded contour the
accessible integral is `φ̂(iε)` rather than `φ̂(0)`; `Γ` is **even** in `ω`, so
the bias is `O(ε²)` and two contours Richardson-extrapolate it — both the raw
and corrected numbers are reported so the correction is visible.

## 7. What this closes, and what it does not

**Closes:** roadmap step 3, properly. There is now a time-dependent two-wave
invariant built from improved conformal stress tensors on the stable
self-adjoint throat background; it reproduces the known WKB collinear/head-on
result as a limit, with rates; and its departure from that limit is measured,
attributed, and dominated by a single effect — multipath — which is `O(1)` where
every other correction is `1e-7` or smaller. The branch decomposition of PRs
#253–#255 is what indexes it, and PR #258's static defect `𝒲 = −β` is recovered
from the same solve.

**Does not close:**

* **backreaction.** The stress tensor is computed *from* the field and never fed
  back. That is the next step, and it now has a concrete object to feed;
* the boundary data is still four real numbers **chosen**, not derived —
  `shells/junction.py` (PR #249) is what would fix them from matter;
* the throat is still **point-supported**: no interior, no proper length, no
  delay of its own beyond the two geodesic legs;
* everything is a **linear scalar** on a **fixed** background. No topology
  change, no rate, no gravitational radiation.

## 8. What this changes for the next step

The stationary-action and backreaction rounds inherit three things:

* **an object to vary.** `T_{μν}` from a solved field, traceless to `1.9e-15`,
  is exactly what an on-shell action and a backreaction source are made of;
* **a branch index that is not optional.** Any quantity built by integrating
  over the field must state which branches were present. A number quoted without
  it can be off by the full range `0 → 4` — measured here, not argued;
* **a scale where geometric optics stops.** `e* ∼ 1/ω` at the caustic, and
  `𝒩` collapsing in `ωe`. Backreaction near the antipodal focus is precisely
  where a WKB estimate would be worst, and now there is a number for how much.
