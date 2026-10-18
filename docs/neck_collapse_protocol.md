# Focused-wave neck collapse: effective protocol (draft)

**Status: draft for discussion. Not frozen.** Freeze it in a local commit before the first coupled run, and publish its hash with the results. As in #340 and #341, that is not external preregistration.

**Revised after the review of `d864d0b`:**
- one energy zero before and after surgery;
- the surgery constraint on the complete surviving Hamiltonian;
- a conservative Hamiltonian retraction variant;
- the waist domain continued past r = a, with an outward stop;
- a threshold search that does not assume monotonicity, with brackets and censoring;
- a post-collapse horizon separate from the collapse search, with counterfactual attribution for H4.

## Question
Can a focused scalar wave push a metastable neck past its collapse threshold? If it can, where do the neck's stored energy and angular momentum go?

The picture under test is a handle whose neck carries surface energy. The neck is metastable: an energy barrier separates it from a lower-energy state with the neck gone. A wave focused onto its mouths delivers energy to the neck's breathing mode. Past the barrier, the waist contracts until its walls meet (self-contact). The handle then splits, the two halves retract into the sphere, and the released energy leaves as a wave.

This protocol tests that chain in the round-S² effective model of #339, #340 and #341. It is classical and conservative throughout. Topology change enters only through an explicit, stated rule (§4), never through the dynamics.

## 1. What is already known, and what this does not claim
- **Soap-film physics supports the mechanism.**
  - A catenoid neck between two coaxial rings of radius a, a distance 2h apart, is metastable for 0.5277 < h/a < 0.6627. Over that range two flat discs have less area than the stable catenoid, and a barrier (the unstable thin catenoid) separates them. This is Goldschmidt's result, checked numerically in §2.
  - A long tube at fixed pressure is unstable to pinching (Plateau–Rayleigh).
  - So a metastable neck is physically realisable only as a short neck.
- **In relativity a throat needs support, and topology change is singular.**
  - Ellis-type throats are dynamically unstable: a pulse collapses them or makes them expand (Shinkai & Hayward 2002; González, Guzmán & Sarbach 2009).
  - Classical spatial topology change requires closed timelike curves (Geroch 1967), or a singularity under the energy conditions (Tipler 1977).
  - In #200's terms, pinching a 3D throat's neck is an index-3 topology change, the class #200 treated as causally discontinuous.
- **Therefore this protocol does not claim:**
  - a GR neck pinch;
  - a derivation of topology change;
  - neck stability itself, which is an assumption parameterised in §2;
  - any quantization.
  
  The "mass to energy" reading means only that the neck's stored surface energy is released as field energy.
- **Exact refocusing is a property of the operator.** The conformal scalar (−Δ + ¼) refocuses exactly at the antipode. A tension-membrane wave (ω_l² ∝ (l−1)(l+2)) would not. Every result here inherits that assumption.

## 2. Model (Stage A: effective neck mode)
Use the unchanged #340 Hamiltonian, with M = 0.1 per pair, cutoff 32, κ = 1, μ = 0.05 and Ω = 5. Make one change: each pair's symmetric shape mode becomes a metastable neck waist.

**Coordinates.** For pair h, let Q_h = (q_h+ + q_h−)/2 and P_h = (q_h+ − q_h−)/2. In #340 the shape energy μΩ²(q_h+² + q_h−²)/2 equals μΩ²(Q_h² + P_h²), and the shape kinetic energy is μ(Q̇_h² + Ṗ_h²).

**Waist.** Let r_h = r₀ + ζQ_h. Replace μΩ²Q_h² by

  V(r_h) = γ [A(r_h) − A(r₀)],

where A(r) is the area of a reduced shape family:
- the family is the cosh profile ρ(z) = r cosh(βz), with β = arccosh(a/r)/h, so that ρ(±h) = a;
- it contains both exact catenoids;
- its stationary points are exactly the thin and thick catenoid waists r_b < r₀.

The P_h mode, the mouth–field coupling κ(q − ū)²/2 and everything else stay as in #340.

**Domain.** The cosh family is real only for 0 < r ≤ a. For r > a, use its analytic continuation ρ(z) = r cos(βz), with β = arccos(a/r)/h: barrel-shaped necks, since cosh(ix) = cos x.
- The area depends on β², and arccosh(a/r)² continues analytically to −arccos(a/r)², so A(r) is smooth through r = a. The Stage-0 script checks this: A′(a±) = 5.026548 and A″(a±) = 27.5902 from both sides.
- Outward excursions are reachable. At the baseline, r = a is reached at Q = (a − r₀)/ζ = 0.1118, after a potential rise of only 0.0259.
- Beyond r = a the potential climbs steeply: V = 0.101 at r = 1.2a and V = 0.157 at r = 1.3a (Q = 0.243).
- **Outward limit.** Stop the run at r = 1.3a, bracketed like the collision stop; that run is INCOMPLETE (OUTWARD_LIMIT).
- No substep may evaluate V at r ≤ 0 or r > 1.3a. A step that would do so is the crossing step: the run either stops at the previous state (outward) or takes the surgery state of §3 (inward).

**Nesting.** Choose γ and ζ so that:
- the small-oscillation limit is #340 exactly: γA''(r₀)ζ²/2 = μΩ², with waist inertia 2μ/ζ²;
- the barrier ΔV = γ[A(r_b) − A(r₀)] takes the preset value.

One parameter, ΔV, then fixes both γ and ζ. The barrier in field units is Q_b = (r_b − r₀)/ζ < 0.

**Baseline geometry: h/a = 0.60** (a = 1 in neck units).

| quantity | value |
|---|---|
| thick (stable) catenoid r₀ | 0.745071 |
| thin (barrier) catenoid r_b | 0.351888 |
| A(r₀) | 6.999643 |
| A(r_b) | 7.207912 (barrier = 2.98% of A(r₀)) |
| two discs, 2πa² | 6.283185 (release = 10.2% of A(r₀)) |
| family area at r = 1e-4 | 6.375373 |

The family reaches the two-disc value only logarithmically, which §4 accounts for.

**Stage-0 gates.** `experiments/closure_ledger/neck_stage0_check.py` reproduces every number in this section and asserts these gates. Before any coupled run, all of the following must hold:
- the code reproduces the critical h/a = 0.66274 to 1e-4;
- the area crossing h/a = 0.52770 is reproduced to 1e-4;
- the family's stationary points equal the catenoid waists to 1e-6 at h/a = 0.55, 0.60 and 0.65.

**Barrier sweep.** ΔV ∈ {0.0025, 0.005, 0.01, 0.02} in model energy units. #340's prepared packet carries 0.1. The baseline is ΔV = 0.01. With A''(r₀) = 10.012 and A(r_b) − A(r₀) = 0.208270 at h/a = 0.60, the nesting conditions give:

| ΔV | γ | ζ | barrier Q_b (field units) | harmonic estimate −√(ΔV/μΩ²) | waist inertia 2μ/ζ² |
|---|---|---|---|---|---|
| 0.0025 | 0.01200 | 4.561 | −0.086 | −0.045 | 0.0048 |
| 0.005 | 0.02401 | 3.225 | −0.122 | −0.063 | 0.0096 |
| 0.01 | 0.04801 | 2.280 | −0.172 | −0.089 | 0.0192 |
| 0.02 | 0.09603 | 1.613 | −0.244 | −0.126 | 0.0385 |

The well is softer than its harmonic approximation, so the barrier sits about twice as far out. For scale, #340's mouths reach shape amplitudes of order 0.07 at λ = 1.

**Integrator.** Keep #340's exact harmonic drift, and put V(r₀ + ζQ) − μΩ²Q² into the kick. The splitting stays symplectic, and with V quadratic it reduces to #340 exactly.

## 3. Collapse event
Collapse of pair h occurs when the waist reaches the surgery radius r_s = 0.05a, moving inward. Past r_b, V decreases monotonically to r_s, so crossing the barrier with inward velocity leads to collapse.

Detect the crossing after every integration step.
- **The surgery state is the first integration state with r ≤ r_s.** Record the bracket [t_prev, t_c]. Every quantity in §4 is evaluated at that state.
- If the crossing step would evaluate V at r ≤ 0, halve that step until the state lands in (0, r_s].

A waist that crosses r_b and returns is a **sub-threshold excursion**, not a collapse. Record each excursion's depth.

## 4. Surgery rule (imposed, not derived)
At the collapse time t_c:

1. **Remove pair h:**
   - its axis and L_h;
   - its waist mode Q_h and antisymmetric mode P_h, with their momenta;
   - its two mouth–field couplings.

   The field domain is unchanged: it is the full sphere, as in #340.
2. **One energy zero throughout.** Keep the shifted convention of §2, V(r) = γ[A(r) − A(r₀)].
   - The two closed caps left after surgery then have V_caps = γ[2πa² − A(r₀)], which is negative: −0.034400 at the baseline.
   - The surface energy released is V(r_s) − V_caps = γ[A(r_s) − 2πa²]. At the baseline, with r_s = 0.05a, it is 0.023287.
   - Never subtract the absolute cap energy γ·2πa² from the shifted potential. At the baseline that would give −0.31280 in place of +0.02329.
3. **Energy balance on the complete Hamiltonian.** Let H_before be the full Hamiltonian at the surgery state: field, the surviving pair (kinetic, shape and interaction terms), and the removed pair (kinetic, shape, V(r_s) and its interaction terms). Let H_after(δa, δp) be the full Hamiltonian after removal: field (a + δa, p + δp), the surviving pair with its interaction energy κΣ(q − ū)²/2 evaluated at a + δa, and V_caps. Require H_after(δa, δp) = H_before. A nonzero δa changes the surviving mouths' ū, so a field-energy constraint alone does not conserve the total.
4. **Angular momentum.** J_field(a + δa, p + δp) − J_field(a, p) = L_h, with J = −Re⟨p, D a⟩. The surviving pair's L is unchanged at the impulse, so total J is conserved.
5. **Deposition form.** Apply an impulsive field displacement δa and velocity δp at the two former mouth sites s = ±. Both lie in the span of {w_hs, D_i w_hs}: the axisymmetric footprints and their rotational dipoles.
   - Including δa lets the deposit carry its own angular momentum, through the ⟨δp, D δa⟩ term, instead of relying on the field already present.
   - Take the solution of minimum norm in the field-energy metric that satisfies the two constraints of steps 3 and 4 exactly.
   - If no real solution exists, report SURGERY_INFEASIBLE. Do not substitute another rule.
6. **Variant: Hamiltonian retraction (conservative by construction).** Do not spread a prescribed force over time: a fixed impulse delivered over an interval τ does not deliver fixed work. A resting oscillator given total impulse J that way receives (J²/2)·sinc²(ωτ/2), not J²/2. Instead:
   - At t_c, keep pair h's rotor and its two mouth oscillators, now cap oscillators, and drop the waist potential.
   - Reset Q_h to the cap rest value 0, with harmonic potential μΩ²Q_h² from then on. The cap oscillators' energy is measured from V_caps, the closed-cap energy in the common zero, so they hold V_caps + μΩ²(Q_h² + P_h²) + μ(Q̇_h² + Ṗ_h²).
   - Choose Q̇_h so that the complete Hamiltonian is unchanged, including the change in the pair's own interaction terms. This is an internal conversion of the released surface energy into cap motion. The field and J are unchanged. If no real Q̇_h exists, report RETRACTION_INFEASIBLE.
   - Evolve the unchanged #340 Hamiltonian for a retraction time τ_r ∈ {0.05, 0.1, 0.2}. The caps ring down into the field through the existing coupling, and the full energy and angular-momentum ledger holds at every step.
   - At t_c + τ_r, apply steps 1–5 to what remains in the pair (the reservoir): rotor kinetic energy, cap-oscillator energy and interaction terms. The state after removal again holds V_caps, as in the baseline.
   - Report the reservoir's energy and L_h at the hand-over. τ_r = 0 recovers the baseline rule.

   For scale, the soap-film retraction speed is √(2γ/σ) (Taylor–Culick). Here τ_r is a free parameter, and results must be insensitive to it (§7).

## 5. Runs
**Packets.** Use #340's preparations scaled by λ, so the energy is 0.1λ².

**Threshold search.** Bisection would assume that collapse is monotone in λ. That has not been established: amplitude-dependent phase, recoil and repeated foci can open and close collapse windows.
1. Scan 24 log-spaced values of λ in [0.1, √10], plus λ = √10. Classify each run as COLLAPSE, NO_COLLAPSE or INCOMPLETE (collision stop or outward limit).
2. Refine every boundary between adjacent COLLAPSE and NO_COLLAPSE points by bisection to 1% in λ.
3. λ* is the lower edge of the lowest collapse interval, reported as a bracket [λ*_lo, λ*_hi]. Report every interval: a non-monotone window is a result, not noise.
4. A boundary next to an INCOMPLETE point is UNRESOLVED. If nothing collapses up to √10, report λ* > √10 (censored).

Run:

| case | preparation | role |
|---|---|---|
| C1 | symmetric x-packet alone, normalised to energy 0.1λ² by itself; pair 2 at the focus, pair 1 at 90° as spectator | focused, no first-order recoil |
| C2 | directional x-packet alone, same layout | focused, with recoil carrying the pair off-focus |
| C3 | C1 packet with pair 2 rotated to ±y, off the source–target axis | off-focus control |
| C4 | #340's two-packet preparation (directional and symmetric) | both necks at risk |
| C5 | mirror of C2 (η → −η) | λ* must agree to 1e-8 relative |
| N1 | empty field, waist at r₀ | no collapse over the horizon |
| N2 | κ = 0 at λ = √10 | no collapse |
| N3 | V replaced by its quadratic Taylor term, at δ = 0 | must reproduce #340's directional CSV to 1e-12 (nesting gate) |

**Start convention.** Backpropagate the prepared packets by π + δ, with δ = 1.2, so that both the source focus (at t = δ) and the target focus (at t = π + δ) are complete events with resting, unexcited mouths beforehand. #341's review found that starting in mid-focus (δ = 0) halves the leading exchange coefficient. Also report δ = 0 for comparability with #340 and #341, labelled as the half-focus convention.

**Horizons.**
- **Collapse search:** to 2π + δ + 1.2, so it includes the source-side return focus at 2π + δ. Only collapses before this horizon count.
- **Post-collapse observation:** after a collapse, continue the run to t_c + π + 0.6, subject to the validity limits. If a limit stops the run earlier, H4 is UNAVAILABLE for that run.
- **Sampling:** exactly at δ ± 0.6, at kπ + δ ± 0.6 and kπ + δ for k = 1, 2, and at t_c + π/2 and t_c + π.
- Stop at #340's disjoint-mouth limit (0.65 rad), bracketed as in #341.
- A collision stop makes that run INCOMPLETE.

**Recorded per run:**
- every energy term, plus all-step maxima of the energy and angular-momentum errors;
- per-mouth impulses in every window;
- the waist trajectory and the deepest excursion;
- t_c and its bracket;
- the deposited E and J, and δa, δp;
- at t_c + π, the refocused field amplitude and energy in 0.3-rad caps around both former mouth sites;
- after t_c, the other pair's impulse.

## 6. Hypotheses, set before running
- **H1 (focus selectivity).** The off-focus threshold energy is at least 5× larger than the focused one, tested conservatively on brackets:
  - holds if λ*_hi(C1) < λ*_lo(C3)/√5;
  - if C3 is censored (no collapse up to √10), holds if λ*_hi(C1) < √2;
  - UNRESOLVED if either bracket touches an UNRESOLVED boundary.

  #339 found that the focus deposits 40× to 440× more shape energy than an off-focus pass.
- **H2 (symmetric collapses more easily than directional).** λ*(C1) < λ*(C2). A symmetric single packet gives no first-order recoil, so the pair stays on the focus. #339's focus shape-energy maxima were 0.0164 (symmetric) against 0.0133 (directional).
- **H3 (where collapse happens).** At threshold, t_c falls in a focus window: target [π + δ ± 0.6], or source return [2π + δ ± 0.6]. Report which, and whether the source focus at t = δ pre-loads the neck.
  - Both mouths of a pair drive the same waist mode Q. With Ω = 5 the source and target foci are 2.5 ringing periods apart, so whether the second focus adds to or cancels the first depends on the waist's phase.
  - Record the waist phase at each focus.
  - This is distinct from #341's rotational exchange, which accrues at mid-transit crossings rather than at foci.
- **H4 (re-convergence).** The released pulse re-converges on the two former mouth sites at t_c + π. Each site's antipode is the other site, and that site is where the pulse's other half started.
  - Cap energy alone cannot test this, because the original packets also refocus.
  - **Attribution:** evolve a paired counterfactual, the identical post-surgery state with the deposit omitted, under the same Hamiltonian. This is for attribution only; the counterfactual does not conserve energy. Use the difference field.
  - **Test:** the difference field's energy in the two 0.3-rad caps at t_c + π must exceed its value at t_c + π/2 by a preset factor of 5.
  - **Check:** also evolve the deposit alone in the free field, with no mouths.
- **H5 (barrier scaling).** For single-pass transfer, λ*² ∝ ΔV across the sweep. A significantly weaker growth means multi-focus pumping. Report the log-log slope with its uncertainty.

## 7. Gates
**Numerical (otherwise UNRESOLVED_NUMERICALLY):**
- relative energy error < 1e-4 over each trajectory, including across the surgery;
- normalised angular-momentum error < 1e-8;
- surgery ledger: the complete Hamiltonian (field, surviving pair with its interaction terms, caps) and total angular momentum unchanged across the surgery to 1e-12;
- halving dt and raising the cutoff to 40 each change λ* by < 1% and t_c by < 0.01;
- the τ_r variants keep the complete ledger within the numerical tolerances above at every step of the retraction stage, and change the re-convergence factor by < 10%.

**Controls:**
- the N1–N3 gates hold;
- the C5 mirror holds;
- the Stage-0 gates hold.

**Verdict:**
- **SUPPORTED_FOR_THIS_EFFECTIVE_FAMILY:** H1 holds at the baseline ΔV with converged thresholds, and every ledger closes.
- **NOT_SUPPORTED:** C1 and C2 both censored (no collapse for λ ≤ √10), or H1 fails while converged.
- **UNRESOLVED:** H1 UNRESOLVED, or the C1 threshold bracket touches an UNRESOLVED boundary.
- **Reported, not gated:** H2–H5, each with its uncertainty and stated as confirmed or not.

## 8. Practices carried over from #340 and #341
- Stop the integration (not the display) at every validity limit, bracketed at step resolution.
- Write archives with sorted keys, so a rerun reproduces the manifest bytes.
- Any statement about a cross term between packets needs amplitude scaling (λ, λ/2, λ/4) and a leading-order reference.
- Draw large motion on a fixed sphere.

## 9. Drawing
- Use the fixed-sphere renderer from the #340 review, with colour-only field and moving cutouts.
- Draw each neck with its waist scaled by r_h(t)/r₀, so that contraction and pinch are visible.
- At t_c, replace the neck by two retracting cap bumps whose height follows the deposited footprint, then show the outgoing ring and its re-convergence at t_c + π.
- Label the surgery frame "imposed topology-change rule".
- Klein-handle crossings are drawing artifacts and are never treated as contact.

## 10. Stage B (separate protocol, before physical interpretation)
Before Stage A thresholds are read physically, solve the neck as an axisymmetric membrane: a surface of revolution with tension and surface density between two rims driven by an imposed rim pulse. Compare its collapse threshold and released energy with the reduced potential V(r) and waist inertia; they must agree within 10%.

Stage B removes the one-mode reduction. It does not remove the surgery rule, which any continuum calculation also needs at the singularity.

## 11. Scope
This is a classical, conservative test in an assumed effective Hamiltonian. The topology change is imposed (§4). Its 3+1 analogue is a singular, index-3 event in #200's classification. A positive result would establish, for this family, a focused-wave route from stored neck energy to outgoing field energy, with exact bookkeeping. It would not establish neck dynamics in GR, a particle-to-radiation law, or quantization.
