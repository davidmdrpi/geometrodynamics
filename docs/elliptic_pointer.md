# A classical record at an elliptic mouth

At equal channel intensities the field exerts zero force. The switch supplies work and the bistable detector amplifies its own initial offset. At unequal intensities the field biases the detector by gz. For the stationary two-mode preparations the exact outcome law is P₊(z)=$F(w_*z)$. Choosing the uniform width 0.12 fixes the maximum discrepancy of 0.125 on the chosen input grid before integration. Reproducing it checks the integrator and basin assignment. The independent question is what physics prepares F. Restoring the eight-site hopping tests sensitivity to the two-mode approximation.

The [protocol commit](https://github.com/davidmdrpi/geometrodynamics/commit/9be2cc4ae401bc3fba637416c234988dabf7877e) precedes the results commit in Git history. That ordering alone does not date the run or demonstrate independent public preregistration. The eight-site calculation was added after the primary run. This follows [#356](https://github.com/davidmdrpi/geometrodynamics/pull/356) and uses the local ellipse coupling from [#241](sigma_z_readout_capstone.md). It addresses the distinction in [#236](https://github.com/davidmdrpi/geometrodynamics/pull/236) between an operator on a mode pair and the outcome of one trial. The [earlier Born probe](classical_born_rule.md) already showed how a selected hidden-variable measure can give Born probabilities. Here a mechanical trajectory, energy barrier and explicit work budget determine each record.

## Field, detector and action

Write the normalized classical winding amplitudes as a column u. The ellipse at orientation φ selects the Hermitian matrix with upper-right entry exp(iφ), conjugate lower-left entry and zero diagonal. In its eigenbasis denote the amplitudes by a₊ and a₋. Then z = |a₊|² − |a₋|². Equal winding intensities with relative phase δ give z = cos(δ+φ). Equal **measurement-channel** intensities give z=0 regardless of their relative phase. The preparation fixes the field norm to unity. The underlying classical field amplitude is continuously scalable.

In dimensionless units the action during the write and latch stages is

$$
 S=\int dt\,\left[\frac{i}{2}(a^\dagger\dot a-\dot a^\dagger a)
       +\frac12\dot q^2-V(q,t)+g(t)qz\right]+S_{\rm load}.
$$

Here q is the signed detector displacement along a selected ellipticity axis. During writing V=q²/2 and g=0.1(2−√2), using #241's eight-site slope with the chosen ellipticity ε=−0.1q. At the switch, g becomes zero and V becomes (q²−1)²/4. The latch potential, gain scale, external switch and initial preparation are supplied by the apparatus. The field uses classical canonical quadratures. No projection or random field collapse is applied.

The passive load is a unit-wave-speed string with action (γ/2)∫dt dx [(∂ₜy)²−(∂ₓy)²], attached by y(0,t)=q(t). Initially the string has uniform displacement q₀, zero velocity and no incident wave. Its outgoing solution gives reaction −γq̇ and power γq̇², with γ=0.4. A string longer than half the observation time has the same response before an end reflection returns. Asymptotic relaxation assumes an outgoing load. A permanent record in an isolated compact universe is a separate physical question.

The write equations give q̈+γq̇+q=gz and ȧ±=±igq a±. The field phases respond while the channel populations are exactly conserved in this reduction. At the switch the stored-energy change at fixed state is external work. During both stages the source-free budget is mechanical plus field energy plus radiated energy minus external work. This apparatus can make a record even for a balanced field because the switch supplies energy.

## Individual outcomes and the basin boundary

Let Ω=√(1−γ²/4), T=π/Ω and ρ=exp(−γT/2). For an initially stationary pointer the exact write solution reaches

$$
 q(T)=(1+\rho)gz-\rho q_0,\qquad \dot q(T)=0.
$$

All sampled tuned-time launches satisfy |q(T)|<√2. Every nonzero launch then has energy below the central barrier of the unbiased latch. Dissipation keeps it on its side and brings it to that minimum. The trajectory at q=q̇=0 stays at the unstable stationary point. The numerical readout tests distance, speed and energy at T+80 and counts unresolved trials in the denominator.

| Same balanced field, z=0 | q at switch | q after 80 | Record |
|---:|---:|---:|---|
| q₀=−0.06 | +0.0315972 | +0.9999999 | positive |
| q₀=0 | 0 | 0 | none |
| q₀=+0.06 | −0.0315972 | −0.9999999 | negative |

The switch delivers 0.25 energy units to the centered preparation. It radiates zero and gives no record. Each displaced example radiates 0.250802, including the initially stored spring energy, and switches with work 0.249002. The initial pointer microstate selects the individual balanced outcome. The field retains both channels after the record forms.

## Outcome law and controls

For a continuous initial-displacement distribution with CDF F, the measured positive fraction at the specified write time is exactly

$$
 P_+(z)=F(w_*z),\qquad w_*=(1+\rho)g/\rho=0.16981365
$$

on preparations satisfying the barrier condition. An atom exactly at the threshold instead gives no record. Therefore Born frequencies over the whole interval require $F(w_*z)=(1+z)/2$ for every z in [−1,1], which fixes a uniform displacement law on $[-w_*,w_*]$ for this protocol. This conclusion is scoped to the frozen action and stationary preparations. Symmetry alone fixes the half split at balance for a symmetric continuous preparation.

| Channel difference z | Target | Primary uniform, width 0.12 | Matched uniform, width $w_*$ | Projected circle, width $w_*$ |
|---:|---:|---:|---:|---:|
| 0 | 0.500 | 256/512 | 256/512 | 256/512 |
| 0.25 | 0.625 | 347/512 | 320/512 | 297/512 |
| 0.50 | 0.750 | 437/512 | 384/512 | 341/512 |
| 0.75 | 0.875 | 512/512 | 448/512 | 394/512 |
| 1 | 1.000 | 512/512 | 512/512 | 512/512 |

The negative-input counts follow by reflection and are archived. Every primary ensemble trajectory latched. Deterministic midpoint quantiles integrate the declared preparation law, so these counts have quadrature error. The exact primary continuum curve saturates at |z|≈0.707, and already gives unit probability at z=0.75. The discrepancy of 0.125 persisted at 256 and 512 quantiles. All tuned-time basin assignments agreed with the evolved trajectories. Moving the switch to 0.9T and 1.1T gave maximum discrepancies of 0.09375 and 0.19922, respectively.

Choosing $w_*$ after deriving the basin gives the Born curve analytically for the uniform preparation. Keeping that width and projecting a uniform circle angle gives P₊=1/2+arcsin(z)/π. At z=0.5 it gives 2/3. Exact uniform-law agreement is a calibration of the apparatus distribution. The action supplies no preparation principle selecting that calibration.

## Eight-site control and reproduction

The additional control evolves the eight complex site amplitudes with link energy Σⱼ |uⱼ₊₁−uⱼ|²/Rⱼ² and Rⱼ=1−0.1q cos(2χⱼ,mid). The pointer feels minus the actual q derivative of that energy. At disconnection the field returns to the free ring and the mechanical latch retains its switch position and speed. At 256 quantiles the full ring changed one positive count at each of z=±0.5. All 2,304 trials latched. The largest population outside k=±1 at disconnection was 5.35×10⁻⁵ or less. Maximum norm and energy residuals were below 1.4×10⁻¹³ and 7×10⁻¹². The ellipse orientation is fixed and the preparations use the real eigenchannel amplitudes specified in the probe.

Run from the repository root with NumPy and SciPy. The two modules print the [archived measurements](../experiments/closure_ledger/runs/20261123_elliptic_pointer/).

```bash
export OPENBLAS_NUM_THREADS=1
python -m experiments.closure_ledger.antipodal_network.elliptic_pointer
python -m experiments.closure_ledger.antipodal_network.elliptic_pointer_ring
python -m unittest discover -s tests -p test_elliptic_pointer.py
```

This experiment gives an explicit mechanical readout for a local mode of a prescribed mouth. The bistability, load and switch remain apparatus assumptions. A throat derived from a gravitational action, an intrinsic preparation law and a repeatable field-state measurement remain open.
