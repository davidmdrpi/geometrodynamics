# Does the timelike bulk stress current guide a mixed winding?

The fiber-averaged energy density stays positive at the balanced midpoint where the signed winding density vanishes. The energy flux divided by that density approaches the positive norm-current velocity in a nonrelativistic regime. At finite carrier energy they differ, both on the [#358 field](pointer_priors_and_transport.md) and in an independently evolved five-dimensional scalar. A small difference between the channel potentials can produce a finite energy velocity at a point with exactly zero norm current.

## Specifying the comparison

[#209](measurement_sector.md) supplies a co-moving Schrödinger force window. It does not fix a unique relativistic apparatus, carrier energy, or timelike symmetry. Here the 5D action is the minimally coupled free complex scalar, with two spectator directions, periodic fiber, and externally specified metric

$$
 ds^2=-c^2dt^2+dx^2+dy^2+dz^2+R^2(d\chi-A\,dt)^2,\qquad R=1/c .
$$

The radius sets the effective mass to one for |k|=1 in units with action scale one. The time vector is the fixed laboratory $\partial_t$. The one-form shift gives $D_t=\partial_t+A\partial_\chi$ and the reduced force potential $V_k=kA$. In the driven region $A=-F(t)x$ with #209's parameters. Space is periodic on [−56,56), with a smooth first-derivative taper from |x|=36 to 54 before the seam. For c=8,12,16 the time vector is everywhere timelike. A constant offset of A allows an explicit check of the choice of laboratory time generator. Other five-dimensional embeddings require their own stress tensors.

The #241 elliptic mixer at ε=0.02 is applied as in #358. For the θ=π/8 comparison the ±1 pair retains all but 0.00005864 of the population. Renormalizing this pair gives plus weight 0.853513. We also test an exactly balanced pair. A continuous KK fiber gives a different rest energy to |k|=3, so its leakage is excluded explicitly. The experiment starts at the end of mixing.

## Which density and flux?

For the canonical Hilbert tensor $T_{AB}=\mathrm{Re}(\partial_A\Psi^*\partial_B\Psi)-g_{AB}|\nabla\Psi|^2/2$, average over the fiber **before** taking either ratio. For each mode write $b_k=D_t\Psi_k$, $d_k=\partial_t\Psi_k=b_k-ikA\Psi_k$, and $s_k=\partial_x\Psi_k$. The time column gives

$$
 E=-\langle T^t{}_t\rangle
 =\sum_k\left[\frac{\mathrm{Re}(b_k^*d_k)}{c^2}
 +\frac12\left(-\frac{|b_k|^2}{c^2}+|s_k|^2+\frac{k^2}{R^2}|\Psi_k|^2\right)\right],
 \qquad S=-\langle T^x{}_t\rangle=-\sum_k\mathrm{Re}(s_k^*d_k).
$$

Thus $v_E=S/E$, whereas #209 uses $v_N=\sum_k\mathrm{Im}(\Psi_k^*s_k)/\sum_k|\Psi_k|^2$. The stress trace in the time component and the distinction between $D_t$ and $\partial_t$ both matter. Direct contraction of the full five-dimensional metric and tensor on a 16-point fiber independently checks the modal implementation. With a dominant common rest carrier $b_k\simeq-ic^2\Psi_k$, the leading terms are $E\simeq c^2\rho$ and $S\simeq c^2J$. This gives the positive norm law to leading order.

**Exact split-potential counterexample.** Consider stationary local modes $u_k(x)e^{ik\chi-i\omega_k t}$. At the sample point give both modes equal intensity, opposite spatial phase gradients ±p, and intrinsic frequency $W=\sqrt{c^4+c^2p^2}$. Let $A=A_0$ there and $\omega_k=W+kA_0$. The local KG equation fixes $u_k''=[c^2-(\omega_k-kA)^2/c^2]u_k$. These data therefore apply even in a smooth stationary connection gradient. They are exact global plane waves for constant A. With unit total norm the gradient trace vanishes at the point and

$$
 E=W^2/c^2,\qquad S=pA_0,\qquad v_N=0,\qquad v_E=c^2pA_0/W^2 .
$$

At c=8, p=1 and $A_0=0.1$, the channel-energy gap is 0.2 compared to the common rest energy 64. The energy density is 65 and $v_E=1/650$. Both velocities vanish when the offset is zero in this symmetric control. A spatially constant offset can be removed by a rotating fiber coordinate, which also changes the time generator; the energy flow is relative to the specified laboratory generator. The driven test below already gives unequal velocities with zero constant offset.

## Test on the archived superposition and on a KG solution

First evaluate the stress on the exact Gaussian envelopes of #358, restoring the common rest carrier and using the Schrödinger equation for their time derivatives. These fields are approximate relativistic solutions. Second, start the 5D wave equation from the same spatial fields with positive free-frequency data for $D_t\Psi$. In each winding, the exact KG evolution is $D_t^2\Psi=c^2\partial_x^2\Psi-c^4\Psi$. Split the exact free Fourier evolution and the midpoint shift evolution symmetrically. Remove only the common fast phase from the stored amplitudes. The stress uses the restored derivatives.

Density-weighted rms velocity differences at t=1, near the maximum imposed gradient, for θ=π/8:

| Evolution | c | Mesh / time step | Energy velocity minus norm velocity on the same field |
|---|---:|---:|---:|
| #358 Gaussian envelope lift | 8 | 1024 / analytic | 0.00794456 |
| Full KG | 8 | 1024 / 0.01 | 0.00794260 |
| Full KG | 8 | 2048 / 0.005 | 0.00794336 |
| Full KG | 12 | 1024 / 0.01 | 0.00354619 |
| Full KG | 16 | 1024 / 0.01 | 0.00199710 |

The rms includes points with total density exceeding one millionth of its peak. Doubling c reduces the difference by about four. Mesh refinement changes the c=8 result by less than one millionth. For the exactly balanced packet at t=1, x=0, the two velocities vanish by reflection and E=12.6707. At the same time and x=0.984375 the energy and norm velocities are 0.025234 and 0.032873, respectively. The instantaneous potentials coincide after the drive turns off; at t=2 the mixed packet still has rms difference 0.005856. Common instantaneous potential is insufficient for equality on a packet with gradients and frequency spread.

During the prescribed drive $\partial_t$ is not Killing. On shell, the coordinate energy balance has a metric work source $-\sqrt{-g}\,T^{AB}\partial_tg_{AB}/2$. In this reduction the integrated source is $\int Q\,\partial_tA\,dx$, with $Q=\sum_k\mathrm{Re}(b_k^*(-ik\Psi_k))/c^2$. At c=8 the drive changes the integrated energy by 0.497588; the relative work-budget residual falls from $6.37\times10^{-7}$ to $1.60\times10^{-7}$ on refinement. After the drive the maximum relative energy drift in the archived runs is below $10^{-12}$. Conservation of a stationary energy current does not on its own identify a geometric core's worldline.

## Scope and replay

A timelike generator and positive-energy scalar supply a regular midpoint denominator for these populated states. The mixed-channel guidance postulate of #209 remains an additional input at finite c. This calculation tests a particular 5D lift and an exact local counterexample; it does not model the localized core, its gravitational backreaction, or a closed compact solution of the Einstein equations.

From the repository root:

~~~bash
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.antipodal_network.timelike_current
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -p test_timelike_current.py
~~~

[Archived probe](../experiments/closure_ledger/runs/20261123_timelike_current/probe.json) · [Merged gravity closure record](classical_waist_status.md)
