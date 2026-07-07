# The nonlinear no-signaling audit: the Gisin edge, faced by construction (PR #204)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. This PR executes the no-signaling
> edge of the #200 register item "nonlinear measurement theory beyond
> the linear test-throat regime" (#198's condition 2) — the second of
> the two named frontier items left standing after #203. The companion
> probe measures every claim on the live dynamics.

## 0. The question, and the stakes

The BAM pilot equation is **nonlinear**: the gravitational potential
Φ[ρ] and the order field q[ρ] both feed the density back into the
evolution. By the Gisin/Polchinski theorems, deterministic nonlinear
modifications of quantum mechanics *generically* allow superluminal
signaling through entangled states. That puts a live refutation edge
under the whole #198/#199 interpretive chain:

> If BAM permits superluminal signaling in quantum equilibrium in its
> causal formulation, it is refuted — by special relativity and by
> experiment — regardless of every spectral success upstream.

Neither standard theorem settles this. The **no-signaling theorem**
protects only linear quantum mechanics; BAM cannot appeal to it. The
**Gisin proof** convicts nonlinear theories *via the projection
postulate* — which BAM does not have (measurement is dynamical
transport, #198). So the question must be settled **by construction**,
on the running ψ–Φ–q dynamics, channel by channel. That is what this
audit does.

**The answer, up front:** the edge fired exactly where it should — the
Newtonian model *does* signal instantaneously — and the theory survives
it: the one and only channel the nonlinearity opens is the
**gravitational field itself**, its instantaneity is a property of the
*Poisson approximation* (removed, at zero cost to the Born rule, by the
causal wave-equation completion), and signal-locality holds in quantum
equilibrium at exactly dBB grade, with the out-of-equilibrium Valentini
signal exhibited as the boundary of the guarantee.

## 1. The audit design

Two far-separated throat-packets (separation 30, widths 2) in the 1D
#198 ψ–Φ–q dynamics. **Alice** applies a local unitary — a real
potential pulse ε·e^{−(x−x_A)²/2w²} for t ≤ 0.5 — at her packet.
**Bob** measures the density response
`D_B(t) = ‖ρ_kick − ρ_no-kick‖_{L²(B)}` on his region. Every comparison
is kick-vs-no-kick *at the same couplings*, so common evolution cancels
and D_B isolates what Alice's choice transmits.

First, the edge is confirmed live: the superposition defect of the flow
map is O(1) — `1.9` — for the full dynamics versus `2×10⁻¹⁴` for the
stripped linear control (G = 0, q = 0): BAM is genuinely (not
perturbatively) nonlinear; Gisin's threat genuinely applies.

## 2. The channel, identified (the 1D field sector)

At t = 3 (before any kinematic ψ-tail reaches Bob):

| configuration | D_B(3) |
|---|---:|
| kinematic floor (G = 0) | 3.3×10⁻¹³ |
| Newtonian (instantaneous Poisson) gravity | **1.5×10⁻⁶** |
| half G | 6.9×10⁻⁷ |
| **Φ-clamp control** (same kick, same live q, gravity blinded) | **2.4×10⁻¹³** |

Three measured facts:

1. **The response is real.** With instantaneous gravity, Alice's local
   choice reaches Bob at 5×10⁶ times the kinematic floor —
   instantaneously. The weak-field model signals. Stated plainly.
2. **It is O(G).** Halving G halves it (ratio 2.14) — linear response
   in the coupling, the signature of a physical interaction.
3. **Nothing but Φ carries it.** Clamping Φ to the no-kick history
   while keeping the kick and the live order field collapses the
   response back *to the kinematic floor* (suppression 6×10⁶). The
   order field q is local (diffusive, κ = 0.5) and transmits nothing.

One channel, at gravitational strength, carried by the gravitational
field. This is not a measurement-theoretic pathology (a Gisin-type
collapse nonlinearity acting on configuration-space states); it is
Newtonian action-at-a-distance — which raises the only question that
matters: is the instantaneity the theory's, or the approximation's?

## 3. The cone (the causal completion)

The Poisson equation is the c → ∞ limit of the causal wave equation
`□Φ = 4πG ρ` — and the 5D theory is causal: the #199 chain derived the
guidance current from the contracted Bianchi identity of the 5D
Einstein equations, whose weak-field, gauge-fixed form *is* a wave
equation, not a Poisson equation. Replacing the solver accordingly
(same mean-subtracted spectral operator, so c → ∞ recovers the Newton
runs exactly):

- **Outside the cone: machine floor.** The response at Bob is
  ~10⁻¹⁵ — not small, *machine zero* — until the front arrives
  (e.g. c = 8 at t = 2.5: 5.5×10⁻¹⁵ versus Newton's 9.0×10⁻⁷ at the
  same instant).
- **The front scales as d/c.** First-response times for
  c ∈ {8, 12, 16}: c·t_front = {25.2, 26.4, 27.2} versus the geometric
  kick-to-Bob distance 25 — monotone in c across a factor-2 sweep.
- **c → ∞ recovers Newton.** At c = 20 the t = 3 amplitude is 1.2× the
  instantaneous one (wave-transient ringing accounts for the modulation).

The superluminal signaling of §2 is therefore the action-at-a-distance
**of the Newtonian approximation, not of the theory**: in the causal
completion, the one nonlinear channel is an ordinary retarded
interaction — Alice can "signal" Bob through it exactly as she can with
a flashlight, and no faster.

## 4. The coexistence theorem: retardation costs the Born rule nothing

The #198 equivariance theorem needed exactly one structural property:
every coupling enters the pilot equation as a **real potential**. The
retarded Φ is precisely as real as the instantaneous one. Verified on
the kicked, retarded (c = 12) live dynamics: norm drift 10⁻¹²
(unitary), continuity residual `∂ₜρ + ∇·(ρ∇S)` at integrator error
(1.5×10⁻⁴), and a 20 000-throat Born ensemble transported by v = ∇S
stays at sampling noise throughout (KS ≤ 0.0073 vs noise 0.0071). **No-signaling and Born equivariance hold
simultaneously** — the fix for the one does not touch the other,
structurally.

## 5. The entangled sector (where Gisin's theorem actually lives)

The two-throat entangled state ψ(x₁,x₂) = φ_L(x₁)φ_R(x₂) +
φ_R(x₁)φ_L(x₂) — the #198 *effective* description of the linear
measurement regime (its emergence from the single 3-space wave is the
part of the register item that stays open, §6).

**Linear part.**
- *The theorem, on the discrete flow:* a local unitary on throat 1
  (potential pulse on x₁ only) leaves the x₂-marginal invariant to
  max|Δρ₂| = 3×10⁻¹⁵ — the exact factorization
  Tr₁[(K⊗1)ρ(K⊗1)†] = ρ̂₂, at machine precision.
- *Equilibrium trajectories are signal-local:* 12 000 throat pairs
  sampled from |ψ|² and transported by the dBB flow through a
  branch-crossing show no kick dependence in the x₂-marginal
  (two-sample KS 0.015 vs noise 0.018).
- *Non-equilibrium trajectories signal:* the same flow with a
  branch-only (non-Born) preparation gives KS 0.045 — 2.5× noise.
  Alice's local choice **is** visible at Bob out of equilibrium. This
  is Valentini's signal-locality boundary, reproduced on the BAM
  transport: no-signaling is an *equilibrium* property, with the #198
  H-theorem relaxation as the mechanism that attains equilibrium —
  exactly dBB's epistemic position, no weaker and no stronger.

**Nonlinear part — the Gisin channel, exhibited and confined.**
Dressing the pair with the BAM mean-field gravity (Φ(x) sourced by the
physical-space density ρ₁(x) + ρ₂(x), V = Φ(x₁) + Φ(x₂)) and kicking
throat 1 only, so the linear background is *exactly zero*:

| mean field | x₂-marginal response in Bob's region |
|---|---|
| none (linear control) | ≤ 6×10⁻¹⁵ at all times (the theorem) |
| instantaneous (Newtonian) | nonzero from t = 0.6, rising to 1.9×10⁻⁶ |
| retarded (c = 6) | machine floor until t = 1.9 ≈ d/c = 2.0, then rising |

The marginal-invariance theorem **is** violated by the nonlinearity —
at O(G), immediately, exactly as Gisin says it generically must be —
and the violation is **confined behind the light front** in the causal
completion. At the entangled level too, BAM's nonlinear "signaling" is
ordinary retarded gravity.

## 6. What is and is not established (honest scope)

**Established:**
- the BAM nonlinearity opens exactly one channel from a local operation
  to a distant region — the gravitational field (O(G); Φ-clamp kills
  it; q is local) — in both the plain and the entangled sector;
- its instantaneity belongs to the Poisson approximation; the minimal
  causal completion confines it to the cone (front = d/c at three
  speeds; machine-floor quiet outside; c → ∞ recovers Newton);
- the completion is free: the retarded potential is real, so #198
  equivariance, unitarity, and the Born ensemble survive untouched;
- equilibrium signal-locality and the non-equilibrium Valentini signal,
  both demonstrated on the BAM transport.

**Not established (the conditions):**
1. The wave-equation Φ is the *minimal* causal completion — the
   gauge-fixed weak-field form of the 5D Einstein equations whose
   Bianchi structure #199 verified symbolically. The full GR constraint
   analysis (lapse/shift; gauge potentials vs gauge-invariant
   observables) is not run here; a no-signaling *theorem* at full-GR
   strength would require it.
2. The configuration-space pair is the #198 **effective** description
   of the linear measurement regime — now shown to admit a causal
   gravitational dressing. The derivation of configuration-space
   structure from the single 3-space wave remains the standing open
   item: the register item is **narrowed** (its no-signaling edge is
   audited), not closed (its emergence question is not).
3. Equilibrium is a hypothesis with a mechanism (#198's relaxation),
   not a theorem. Non-equilibrium ensembles signal — a prediction BAM
   *shares with the entire dBB program* (Valentini), not a defect
   specific to it.
4. 1D/2D reductions; c is a model parameter (physically, the light
   speed); the theorems invoked (partial-trace invariance, continuity
   from real potentials) are dimension-blind.

## References

- N. Gisin, *Weinberg's non-linear quantum mechanics and supraluminal
  communications*, Phys. Lett. A 143 (1990) 1.
- J. Polchinski, *Weinberg's nonlinear quantum mechanics and the
  Einstein–Podolsky–Rosen paradox*, Phys. Rev. Lett. 66 (1991) 397.
- A. Valentini, *Signal-locality, uncertainty, and the subquantum
  H-theorem*, Phys. Lett. A 156 (1991) 5; 158 (1991) 1.
- D. Dürr, S. Goldstein, N. Zanghì, J. Stat. Phys. 67 (1992) 843.
  [Quantum equilibrium.]
- The BAM dynamics and Born rule: PR #176–#180 (the ψ–Φ–q functional),
  PR #198 (equivariance), PR #199 (the guidance law from the 5D bulk).

## Reproduce

```bash
python -m experiments.closure_ledger.nonlinear_no_signaling_audit_probe
# Verdict: NO_SIGNALING_SURVIVES_AUDIT_THE_ONLY_NONLINEAR_CHANNEL_IS_THE_RETARDED
#          _GRAVITATIONAL_FIELD_EQUIVARIANCE_INTACT_EQUILIBRIUM_SIGNAL_LOCALITY_AT_DBB_GRADE
```
