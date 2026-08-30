# PR #216's loop, driven by a derived geometry

**8/8 checks pass.**

`traversable_throat` computes `T_ℓ(ω)` from a supported traversable 5D metric. `network.py` carries PR #216's loop eigenvalue. This wires them — through the APIs that already existed — so the closure questions are statements about **the BAM module itself**.

```
Λ_ℓ(ω, Δ) = η_topo · T_ℓ(ω) · e^{iω(d_A + d_B + Δ)}
```

`η_topo` is `NetworkThroat.topological_factor`, and its orientations are **read off `embedding.topology.make_singlet_pair()`** rather than chosen. That is not cosmetic: `ConjugatePair` asserts the two mouths of one throat carry *opposite* orientations, so the scalar channel has `η_topo = −1`.

> ## The result, which reverses an earlier draft of this round

> **One clock offset does serve both the carrier and the packet — at `ω = 1.4617`.** An earlier draft searched with `η_topo = +1`, a *chosen* sign, and reported no root at all; the declared topology shifts `Ψ = θ − ωθ′` by `π` and a root appears.

> **But the verdict is gauge dependent, and this says so.** `Ψ` sweeps only `3.9676` across the band, less than `2π = 6.2832`, so a constant rephasing of the Jost basis can create or remove the root. Neither *closes* nor *never closes* is a property of the geometry alone. What is invariant is `dΨ/dω = −ωθ″` and the total variation.

---

## N0 — every existing API is the one being driven

| `ω` | `\|T\|` | `\|Λ\|` | difference |
|--|--|--|--|
| `0.5` | `0.1265140646` | `0.1265140646` | `0.0e+00` |
| `1` | `0.6414248997` | `0.6414248997` | `1.1e-16` |
| `2` | `0.9965714953` | `0.9965714953` | `0.0e+00` |
| `4` | `0.9999995410` | `0.9999995410` | `0.0e+00` |

`|η_topo| = 1` ⟹ `|Λ| = |T|`. And the dispatch is in `t_AB`, so the pre-existing entry points agree with the derived transfer exactly:

| API | disagreement with the derived `T` |
|--|--|
| `t_AB` | `0.0e+00` |
| `traverse_throat` | `0.0e+00` |
| `loop_eigenvalue_vs_derived` | `0.0e+00` |
| `effective_green_uses_it` | `5.6e-17` |

> An earlier draft of this round dispatched only inside a new derived_loop_eigenvalue, leaving traverse_throat, network_confirmation, projected_kernel, loop_eigenvalue and effective_green reading the TRANSPARENT ports -- so effective_green(derived_throat, ...) was not the G0/(1-Lambda) whose behaviour the round discussed. The dispatch now lives in t_AB, the one primitive all of them already call, so no caller can pick an API that sees a different throat.

> NetworkThroat.__post_init__ REJECTS a derived backend with tau_th != 0, so the traversal leg's free transit phase e^{-i w tau_th} is exactly 1 and arg T is counted once. loop_expansion and r_AA raise instead of answering from the transparent ports, since a smooth barrier has no echo train.

## N1 — where the loop closes

Eliminating `Δ` between phase closure `ω(D+Δ) + θ = 2πn` and group-delay closure `D + Δ + θ′ = 0` gives

```
Ψ_ℓ(ω) = θ_ℓ − ω θ_ℓ′ = 2πn ,   θ = arg(η_topo T_ℓ)
```

which searches over `n` automatically. (`dΦ/dω = 0` is group-delay closure *at the carrier* — the necessary first-order condition for a finite-band packet, not exact packet closure, which would also constrain the amplitude and every higher derivative.)

| channel | `η_topo` | roots | closest approach |
|--|--|--|--|
| scalar | `-1` | `1.4617` | `1.97e-03` |
| spinor | `+1` | none | `1.77e-01` |

> **An earlier draft of this round scanned with eta_topo = +1 and reported NO root, concluding that simultaneous closure was a UV limit never reached at finite frequency. That eta was chosen, not derived. ConjugatePair asserts opposite mouth orientations and make_singlet_pair builds (+1, -1), so the scalar channel carries eta_topo = -1, which shifts Psi by pi. On the declared topology the loop DOES close, at a finite frequency.**

> orientation and wrap_parity are different operations. A scalar sees the deck orientation only, product -1. A spinor also picks up ThroatDefect.spinor_sign() at each mouth, product -1 again, which returns eta_topo to +1 and removes the root. The field solved here is a scalar, so the scalar row is the applicable one -- but the difference is physical, not a convention, and is reported rather than collapsed.

> Sign changes of C are bracketed with Brent against every branch 2 pi n the band reaches, AND every interior local minimum of the smooth objective |e^{i Psi} - 1| is refined with a bounded minimiser. The second search is what would catch a tangential zero; a sign-change scan alone would step past one.

## N1c — what a rephasing can and cannot move

`Ψ` runs over `[-1.0029, 2.9647]`, a total variation of `3.9676` against `2π = 6.2832`.

> Whether Psi reaches a level 2 pi n. Because the swing (3.9676) is less than 2 pi (6.2832), the band does not cover a full branch, so a constant can create or remove the root. NEITHER 'it closes' NOR 'it never closes' is a property of the geometry by itself.

> dPsi/dw = -w theta'' contains no constant, and the total variation of Psi is likewise rephasing independent. A LINEAR reference-plane phase b*w cancels identically from theta - w theta', so the entire residual freedom is one constant.

`dΨ/dω = −ωθ″` verified against its own refinement:

| step | relative error | ratio |
|--|--|--|
| `1e-02` | `5.61e-04` | — |
| `5e-03` | `1.40e-04` | `4.00` |
| `2e-03` | `2.25e-05` | `6.25` |

> The topological part is: eta_topo comes from ConjugatePair's opposite mouth orientations, not from a choice. The Jost basis constant is NOT yet physically fixed -- it needs the finite-mouth matching surfaces to the S^3 exterior. Until then the closure verdict is stated relative to that basis, and is not promoted to a basis-free statement about BAM.

## N1b — the UV tail law is analytic, not fitted

`int V_l ds = (pi/a)[l(l+2) + 9/8]`, so `ω[Ψ − arg η_topo] → −∫V_ℓ ds`. The topological constant must come out first — otherwise `ωΨ` diverges linearly instead of tending to a limit:

| `ω` | `Ψ` | `Ψ − arg η_topo` | `ω(Ψ − arg η_topo)` |
|--|--|--|--|
| `4` | `+2.234610` | `-0.906983` | `-3.627931` |
| `6` | `+2.545932` | `-0.595661` | `-3.573965` |
| `8` | `+2.697059` | `-0.444534` | `-3.556271` |
| `10` | `+2.786766` | `-0.354826` | `-3.548263` |
| `12` | `+2.846263` | `-0.295330` | `-3.543958` |
| `16` | `+2.920361` | `-0.221232` | `-3.539709` |
| `20` | `+2.964705` | `-0.176888` | `-3.537753` |

Predicted `-3.534292`, reached to `0.10%` by `ω = 20`, monotonically.

> **theta = arg(eta_topo T) carries arg(eta_topo) at EVERY frequency, so w*theta diverges linearly instead of tending to a limit. Psi therefore decays to arg(eta_topo), not to zero -- for the scalar channel that is pi, which is the furthest a phase can be from a branch 2 pi n. The loop is LEAST closed in the ultraviolet.**

> The TAIL: Psi - arg(eta_topo) decays as 1/w and not faster, so far in the ultraviolet Psi sits near arg(eta_topo) and, for the scalar channel, cannot reach a branch 2 pi n.

> **Whether Psi crosses a branch at FINITE w on the way to that tail. An earlier draft of this round read the 1/w law as proving no finite root existed. It does not: an asymptotic law constrains the tail, not the interior. On the declared topology there IS a finite root -- see measure_where_the_loop_closes.**

## N2 — searching for perfect transmission, not assuming it away

> A positive barrier CAN have perfect-transmission resonances, so '|T| < 1 at every finite frequency' is not a general theorem. What is established is narrower AND band limited: on [0.05, 12] |R| falls monotonically with no interior zero, so NO PERFECT-TRANSMISSION POINT WAS FOUND ON THE TESTED BAND. Nothing here rules one out below 0.05 or above 12.

Interior minima of `|R|` found: **0**. Smallest `|R|` = `7.157e-11` at `ω = 12.00`, which is the top of the scanned range (True) — i.e. the minimum is the UV limit, not an interior resonance.

> |Lambda| = |eta_topo| |T| and |eta_topo| = 1, so |Lambda| < 1 wherever |T| < 1: on the tested band 1 - Lambda does not vanish and G_eff has no pole there. That is a statement about this band, not about all finite frequencies.

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

## N5 — the phase-sensitive root, refined

> |R|^2 + |T|^2 = 1 to 1e-13 constrains MODULI only, and says nothing about the error in arg T. Psi differentiates arg T, so the phase-sensitive verdict needs its own study -- against the matching edge, the spatial step, and the finite-difference step together.

| variant | `edge` | `steps` | fd step | root | shift |
|--|--|--|--|--|--|
| baseline | `200` | `60000` | `1e-04` | `1.461703899` | `0.0e+00` |
| edge 150 | `150` | `60000` | `1e-04` | `1.461703753` | `1.5e-07` |
| edge 300 | `300` | `60000` | `1e-04` | `1.461704316` | `4.2e-07` |
| steps 30000 | `200` | `30000` | `1e-04` | `1.461704898` | `1.0e-06` |
| steps 120000 | `200` | `120000` | `1e-04` | `1.461703649` | `2.5e-07` |
| fd 1e-3 | `200` | `60000` | `1e-03` | `1.461704151` | `2.5e-07` |
| fd 1e-5 | `200` | `60000` | `1e-05` | `1.461703897` | `2.2e-09` |

> The root moves by at most 9.99e-07 across all three knobs, so it is quoted as 1.4617 -- four decimals, not the twelve the bracketing solver happens to return.

## N4 — the ledger

| claim | verdict | evidence |
|--|--|--|
| eta_topo may be chosen as +1 | **NO -- IT IS DERIVED, AND IT IS -1** | ConjugatePair asserts the two mouths of one throat carry OPPOSITE orientations and make_singlet_pair builds (+1, -1), so the scalar channel has eta_topo = -1. network_mouth_from_defect now maps ThroatDefect onto NetworkMouth instead of taking the signs as arguments |
| one clock offset serves both carrier and packet | **YES, AT A FINITE FREQUENCY** | on the declared topology (eta_topo = -1) Psi = theta - w theta' reaches 2 pi n at w = 1.4617, stable to 1.0e-06 across edge, step and finite-difference refinement. THIS REVERSES the earlier draft of this round, which searched at eta_topo = +1 |
| that verdict is a property of the geometry alone | **NO -- IT IS GAUGE DEPENDENT** | Psi sweeps 3.9676 over the band, less than 2 pi = 6.2832, so a constant rephasing of the Jost basis can create or remove the root. The topological part of that constant is now derived; the basis part needs the finite-mouth matching |
| nothing about the closure survives a rephasing | **TWO THINGS DO** | dPsi/dw = -w theta'' contains no constant (verified to 2.2e-05 with second-order convergence), and so does the total variation of Psi. A LINEAR reference phase b*w cancels identically from theta - w theta' |
| the scalar and spinor channels give the same answer | **NO** | orientation and wrap_parity are different operations. The scalar sees eta_topo = -1 and closes; a spinor also carries spinor_sign at each mouth, giving eta_topo = +1, and does not. The field solved here is a scalar |
| PR #216's loop can be driven by a derived geometry | **YES -- THROUGH THE EXISTING APIS** | the backend dispatch lives in t_AB, so traverse_throat, network_confirmation, projected_kernel, loop_eigenvalue and effective_green all see the derived T. An earlier draft dispatched only in a new parallel function, which left those APIs reading the transparent ports |
| the derived loop needs a separate tau_th transit phase | **NO -- AND IT CANNOT CARRY ONE** | arg T already carries the frequency-dependent transit, and __post_init__ REJECTS a derived backend with tau_th != 0, so the double-counting object cannot be built at all |
| comparing Delta_phase at n = 0 to Delta_group is the test | **NO -- BRANCH DEPENDENT** | branches are 2 pi/w apart, so the earlier 4.14 gap at w = 1 is 2.14 to the nearest branch; eliminating Delta gives theta - w theta' = 2 pi n |
| the closure function decays to ZERO as 1/w | **NO -- IT DECAYS TO arg(eta_topo)** | w[Psi - arg eta_topo] -> -int V_l ds = -(pi/a)[l(l+2) + 9/8] = -9 pi/8 = -3.5343 at l = 0, a = 1, matched to 0.10% by w = 20. The constant must be removed or w*Psi diverges linearly -- an omission invisible at the earlier draft's eta_topo = +1. So Psi tends to pi, the furthest a phase can be from a branch: the loop is LEAST closed in the UV |
| that 1/w tail law implies no finite root | **NO -- IT CONSTRAINS THE TAIL ONLY** | an asymptotic law says nothing about whether Psi crosses a branch at finite w on the way there. It does, at w = 1.4617. An earlier draft of this round read the tail law as a proof of absence |
| a positive barrier forbids |T| = 1 at finite frequency | **NOT A THEOREM** | positive barriers can have perfect-transmission resonances; what is shown is that a direct search on [0.05, 12.0] finds none -- 0 interior minima of |R| -- which is a band-limited negative result |
| the UV falloff constant is fitted | **NO -- IT IS BORN** | Vtilde_0(q) = (3 pi/8a)(3 + a|q|)e^{-a|q|} gives 1 - |T|^2 ~ e^{-4aw}; the local slope descends to -4.145 against the predicted -4.0 |

**What the wiring changes.** The closure questions are statements about network.py itself, through the APIs that already existed. That matters most for the closure verdict, which depends on eta_topo -- and eta_topo is now read off the repository's own ConjugatePair rather than chosen. Doing that reversed the answer.

**Scope.** The benchmark has two asymptotically flat ends at s -> +-infinity, while network.py conceptually has two finite mouths embedded in the closed S^3 exterior. T_l is therefore a whole-throat oracle, not a literal glued finite-mouth solution. The high-frequency normalisation T -> 1 is what makes it usable as an excess transfer factor. A later construction needs finite matching surfaces to the S^3 exterior and their junction stress -- and those should not be smuggled in merely to fit the old MouthPort API.

**Still open.** The history that produces Delta_BA; the finite-mouth junction to the S^3 exterior, which is also what would fix the Jost basis constant and turn the closure verdict from basis-relative into basis-free; and whether the physical probe is the scalar or a spinor, since the two channels answer differently.
