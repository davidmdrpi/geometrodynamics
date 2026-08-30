# PR #216's loop, driven by a derived geometry

**6/6 checks pass.**

`traversable_throat` computes `T_ℓ(ω)` from a supported traversable 5D metric. `network.py` carries PR #216's loop eigenvalue. This wires them, so the closure questions are statements about **the BAM module itself** rather than a reconstruction beside it.

```
Λ_ℓ(ω, Δ) = η_topo · T_ℓ(ω) · e^{iω(d_A + d_B + Δ)}
```

`η_topo` is `NetworkThroat.topological_factor` — the module's own deck orientations and mouth phases. **No separate `tau_th` phase**: a whole-throat `T` already carries the transit in `arg T`, and adding one would double-count the Wigner delay.

> ## The result

> **No finite frequency lets one clock offset serve both a monochromatic carrier and a localised packet.** The branch-free residual `C = Arg exp[i(θ − ωθ′)]` has **no root** over `[0.2, 12.0]` at 900 points, and its decay is analytic: `ωC → −∫V_ℓ ds = -3.534292` `= −9π/8`. So `C` vanishes as `1/ω` — simultaneous closure is a UV limit, the same limit in which `|T| → 1`.

---

## N0 — the loop is the module's own

| `ω` | `\|T\|` | `\|Λ\|` | difference |
|--|--|--|--|
| `0.5` | `0.1265140646` | `0.1265140646` | `0.0e+00` |
| `1` | `0.6414248997` | `0.6414248997` | `1.1e-16` |
| `2` | `0.9965714953` | `0.9965714953` | `0.0e+00` |
| `4` | `0.9999995410` | `0.9999995410` | `0.0e+00` |

`|η_topo| = 1` ⟹ `|Λ| = |T|`, confirmed to `1e-12`. The batched scan path is asserted equal to the scalar `network.py` path, so the continuous searches below the same object the module exposes.

> derived_loop_eigenvalue applies eta_topo * T * e^{i w (d_A + d_B + Delta)} and no tau_th. A whole-throat T already carries the transit in arg T; adding tau_glob on top, as the MouthPort backend must, would double-count the Wigner delay.

## N1 — the branch-free closure residual has no root

Eliminating `Δ` between phase closure `ω(D+Δ) + θ = 2πn` and group closure `D + Δ + θ′ = 0` gives `θ − ωθ′ = 2πn`, so

```
C_ℓ(ω) = Arg exp[i(θ_ℓ − ω θ_ℓ′)]
```

vanishes exactly when one offset serves both. It searches over `n` automatically.

| `ω` | `C_ℓ(ω)` |
|--|--|
| `0.20` | `+3.10855` |
| `2.17` | `-1.82070` |
| `4.14` | `-0.87519` |
| `6.11` | `-0.58503` |
| `8.08` | `-0.44033` |
| `10.04` | `-0.35325` |

Roots found: **none**. Smallest `|C|` = `0.29533`.

> Eliminating Delta between phase closure and group closure gives theta - w theta' = 2 pi n, so C = Arg exp[i(theta - w theta')] searches over n automatically. Comparing Delta_phase at n = 0 against Delta_group -- as an earlier draft did -- is branch dependent: at w = 1 the branches are 2 pi/w = 6.28 apart, so a raw gap of 4.14 is 2.14 to the nearest branch.

> A constant rephasing theta -> theta + c shifts C by c, so this is only meaningful evaluated end-to-end in the network's own convention with its actual eta_topo -- which is why it is computed here rather than from a bare T.

## N1b — and its decay is analytic, not fitted

`int V_l ds = (pi/a)[l(l+2) + 9/8]`, so `ωC → −∫V_ℓ ds`:

| `ω` | `C` | `ωC` |
|--|--|--|
| `4` | `-0.906983` | `-3.627931` |
| `6` | `-0.595661` | `-3.573965` |
| `8` | `-0.444534` | `-3.556271` |
| `10` | `-0.354826` | `-3.548263` |
| `12` | `-0.295330` | `-3.543958` |
| `16` | `-0.221232` | `-3.539709` |
| `20` | `-0.176888` | `-3.537753` |

Predicted `-3.534292`, reached to `0.10%` by `ω = 20`, monotonically.

> **C vanishes as 1/w and not faster, so there is no finite frequency at which one clock offset serves both a monochromatic carrier and a localised packet. Simultaneous closure is a UV limit -- the SAME limit in which |T| -> 1, so the loop's magnitude and its two phase conditions all close only in the ultraviolet.**

## N2 — searching for perfect transmission, not assuming it away

> A positive barrier CAN have perfect-transmission resonances, so '|T| < 1 at every finite frequency' is not a general theorem. What is established here is narrower: over the scanned range |R| falls monotonically with no interior zero, so no finite-frequency perfect-transmission point was found, and |T| -> 1 only in the ultraviolet.

Interior minima of `|R|` found: **0**. Smallest `|R|` = `7.157e-11` at `ω = 12.00`, which is the top of the scanned range (True) — i.e. the minimum is the UV limit, not an interior resonance.

> |Lambda| = |eta_topo| |T| and |eta_topo| = 1, so |Lambda| < 1 wherever |T| < 1: 1 - Lambda does not vanish and G_eff has no pole on the scanned range. PR #216's completed transaction is a high-Q limit approached in the UV, not an attainable resonance.

## N3 — the UV falloff against a no-fit oracle

`Vtilde_0(q) = (3 pi / 8a)(3 + a|q|) e^{-a|q|}`, so first-Born reflection gives `1 − |T|² ~ e^{−4aω}`: slope `-4.0` with nothing fitted.

| `ω` | local slope |
|--|--|
| `2.00` | `-4.7188` |
| `2.75` | `-4.4924` |
| `3.50` | `-4.3337` |
| `4.25` | `-4.2421` |
| `5.00` | `-4.1844` |
| `5.75` | `-4.1456` |
| `6.50` | `-4.1169` |

> An earlier draft quoted a FITTED slope of -4.25 over 1.5 < w < 8. The Born calculation predicts -4a with nothing fitted, and the LOCAL slope descends monotonically toward it (-4.72, -4.43, -4.27, -4.18, -4.13 ...), so -4.25 was the finite-frequency approach to the analytic asymptote rather than a new constant.

## N4 — the ledger

| claim | verdict | evidence |
|--|--|--|
| PR #216's loop can be driven by a derived geometry | **YES, END TO END** | NetworkThroat gained a whole_throat_transfer backend; Lambda comes from network.derived_loop_eigenvalue with the module's own eta_topo, and the MouthPort backend is retained untouched |
| the derived loop needs a separate tau_th transit phase | **NO -- IT WOULD DOUBLE-COUNT** | arg T already carries the frequency-dependent transit; the Fabry-Perot backend needs tau_glob only because its t_AB is an excess factor over free interior propagation |
| |Lambda| = |T| when the topological factor is unit modulus | **CONFIRMED** | agreement to <1e-12 at every probe; |eta_topo| = 1 |
| one clock offset serves both carrier and packet | **NO ROOT FOUND** | C = Arg exp[i(theta - w theta')] scanned over [0.2, 12.0] at 900 points; smallest |C| = 0.2953 |
| comparing Delta_phase at n = 0 to Delta_group is the test | **NO -- BRANCH DEPENDENT** | branches are 2 pi/w apart, so the earlier 4.14 gap at w = 1 is 2.14 to the nearest branch; eliminating Delta gives the invariant theta - w theta' = 2 pi n |
| the closure residual's decay is a fitted observation | **NO -- IT IS ANALYTIC** | w C -> -int V_l ds = -(pi/a)[l(l+2) + 9/8]; at l = 0, a = 1 that is -9 pi/8 = -3.5343, matched to 0.1% by w = 20. C vanishes as 1/w, so simultaneous closure is a UV limit and never finite |
| a positive barrier forbids |T| = 1 at finite frequency | **NOT A THEOREM** | positive barriers can have perfect-transmission resonances; what is shown is that a direct search over [0.05, 12.0] finds none -- 0 interior minima of |R| |
| the UV falloff constant is fitted | **NO -- IT IS BORN** | Vtilde_0(q) = (3 pi/8a)(3 + a|q|)e^{-a|q|} gives 1 - |T|^2 ~ e^{-4aw}; the local slope descends to -4.145 against the predicted -4.0 |

**What the wiring changes.** The closure questions are now statements about network.py itself rather than about a reconstruction beside it. That matters most for the closure residual, whose value shifts under a constant rephasing of the transfer -- so it is only well posed inside the network's own convention, with its own eta_topo.

**Scope.** The benchmark has two asymptotically flat ends at s -> +-infinity, while network.py conceptually has two finite mouths embedded in the closed S^3 exterior. T_l is therefore a whole-throat oracle, not a literal glued finite-mouth solution. The high-frequency normalisation T -> 1 is what makes it usable as an excess transfer factor. A later construction needs finite matching surfaces to the S^3 exterior and their junction stress -- and those should not be smuggled in merely to fit the old MouthPort API.

**Still open.** The history that produces Delta_BA, and the finite-mouth junction to the S^3 exterior.
