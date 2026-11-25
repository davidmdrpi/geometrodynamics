# Does a finite classical pulse leave exactly one durable record?

Two initially ready latches receive the output of the elliptic mixer. Each latch can store a positive energy and remain fired after the field is gone. With the registered pulse energy, energy conservation forbids two fired latches. The apparatus readily makes **zero** records: at the balanced setting both detectors returned to ready in all 512 primary trials. Both registered physical criteria fail. The calibration, energy accounting and refinement checks pass.

This is a test of an operational principle for an ideal detector. Every prepared trial should leave one durable record, and the unconditional right-record probability should follow the angle law. Previous experiments computed a branch intensity, a transported position, or the sign of a single latch. Here each of two output ports has a separate state that can be read later. Both, neither, and unresolved states enter the denominator. [Pre-run protocol](two_detector_records_protocol.md) · [Mixer and transport](pointer_priors_and_transport.md) · [Single latch and outgoing load](elliptic_pointer.md)

## Action and incident preparation

For endpoint displacements $q=(q_R,q_L)$, add two independent string actions to
$$
 L_{\rm mech}=\tfrac12(\dot q_R^2+\dot q_L^2)+\tfrac{\mu}{2}(\dot q_R-\dot q_L)^2
 -U(q_R)-U(q_L),\qquad U(q)=6q^4-14q^3+9q^2.
$$
Each lossless, unit-speed, semi-infinite string has Lagrangian $(\gamma/2)\int_0^\infty(y_t^2-y_x^2)\,dx$, endpoint $y_i(0,t)=q_i(t)$ and impedance $\gamma=2$. The inertial link is a lumped model of one apparatus with two terminations. It allows energy exchange while leaving both static well locations intact. For incoming characteristic velocity $f_i$, the outgoing velocity is $\dot q_i-f_i$. Varying the endpoint gives
$$
 M\ddot q+\gamma\dot q+\nabla U_{\rm pair}=2\gamma f,\qquad
 \dot E_{\rm mech}=\gamma\sum_i[f_i^2-(\dot q_i-f_i)^2].
$$
This uses the same time-independent action on every trial; the pulse is initial data on the leads. Incoming modal intensities and separation into two ports are supplied by #358's mixer and ideal analyzer. The nonlinear memory, impedance and link are effective detector parameters. No geometric waist moves in this calculation.

The ready minimum is $q=0$ with energy zero, the fired minimum is $q=1$ with energy one, and the intervening saddle is $q=3/4$ at $B=1.0546875$. A fixed smooth reference trajectory prepares the pulse that takes one isolated cold port to the fired minimum in $T=2.5$. Substitution in the endpoint equation fixes the waveform and gives incident energy $E_0=1.54284583$. At each angle, the same waveform goes to the two ports with amplitudes $\sqrt{w}$ and $\sqrt{1-w}$ from the full eight-site mixer at ellipticity 0.02. No gain is adjusted by angle. The incoming pulses cease at $T$.

For each port, prepare canonical $(q_i,p_i)$ uniformly inside $18q_i^2+p_i^2\le2e$ in dimensionless units. This is a specified bounded ensemble in the ready well. The primary uses $e=0.02$, $\mu=0.25$ and 512 deterministic quadrature points per angle. Scrambled Sobol points and their port-swapped partners enforce exchange symmetry of the sample; runs at different angles share the same launches. A second scramble uses 1,024 points. A broader preparation uses $e=0.06$.

## A certificate for a record

Read both sides at $t=24$. At a first escape from the present basin pair, one latch must reach the saddle while the other still has at least its own basin's minimum energy. Thus a sufficient locking condition is
$$
 E_{\rm mech}(24)<B+\min\{U_{\min,R},U_{\min,L}\}.
$$
Energy after the pulse is nonincreasing, so this certifies the record for every later time under the stated load. A failure of the strict bound is reported as uncertified. In a finite closed apparatus, the same argument applies only until stored radiation returns. All primary outcomes, including all the misses, were certified.

There is an independent **at-most-one** theorem for the registered preparations. The mass eigenvalues are $1$ and $1+2\mu$, so the kinetic energy is at most $|p|^2/2$. Put $d=\sqrt{2e/18}$. The total initial mechanical energy is at most $2[U(-d)+e]$, even using the looser rectangular bounds of the ellipse. At the broader preparation this is 0.255775. Adding the entire incident pulse still leaves less than the energy 2 needed to occupy both fired basins. This bound holds throughout the trajectory and on the entire preparation support.

For uncoupled ports the same argument works separately. Even at $e=0.06$, a balanced port receives at most $E_0/2+U(-d)+e<0.900$, below the saddle. An open interval of settings around balance therefore gives exactly zero records for every allowed launch. For a coupled cold pair with exactly balanced input, uniqueness and exchange symmetry give $q_R(t)=q_L(t)$. Together with the energy bound this also excludes a single record. The finite warm coupled ensemble is tested numerically.

## What the detectors retained

All nine primary angles are in the archive. The other four angles mirror the first four with left and right exchanged. These are counts out of 512 prepared trials per angle.

| Angle | Right only | Left only | Neither | Both | Uncertified | Ideal right probability |
|---|---:|---:|---:|---:|---:|---:|
| 0 | 490 | 0 | 22 | 0 | 0 | 1.0000 |
| π/16 | 359 | 0 | 153 | 0 | 0 | 0.9619 |
| π/8 | 4 | 0 | 508 | 0 | 0 | 0.8536 |
| 3π/16 | 0 | 0 | 512 | 0 | 0 | 0.6913 |
| π/4 | 0 | 0 | 512 | 0 | 0 | 0.5000 |

The maximum unconditional right-frequency error is 0.845741, against the registered tolerance 0.01. Comparing with the actual mixer output gives the same verdict. At balance the mean outgoing energy is 1.56119, including the mechanical energy initially in the warm detector. The largest final mechanical energy there is below $7\times10^{-15}$. The leads carry the energy away and the detectors return to ready.

Halving the time step changed zero individual outcomes across 4,608 angle–launch pairs. The primary energy-balance residual fell from $2.06\times10^{-6}$ to $1.39\times10^{-7}$ in pulse-energy units. The second scramble gave zero exclusive records at balance and 17 right-only records out of 1,024 at π/8. Broadening the preparation increased the primary π/8 count from 4 to 89 out of 512; all balanced trials still missed. The dark, one-full-port and two-full-port calibration runs certified neither, right only and both, respectively.

## Could a brighter balanced pulse repair the misses?

A separate post-hoc sweep kept the two inputs equal and varied the total incident energy from $E_0$ to $2.25E_0$ in increments of $0.125E_0$. At $1.75E_0$ the coupled apparatus gave **118 right only, 118 left only, 98 both and 178 neither**. Discarding both and neither would give a perfect 50:50 split. In the experiment as prepared, only 236 of 512 trials wrote exactly one record. At $2E_0$ every trial wrote both. The largest exclusive fraction on this grid was 46.1% with the link and 49.6% without it. These are observations on the chosen grid, not an optimization bound for coupled receivers.

For identical **independent** detectors at balance the obstruction is exact for any gain or local response curve. If either latches with probability $p$, the exclusive probability is $2p(1-p)\le1/2$. For a fair ideal pair with exactly one record, the binary detector outputs instead require covariance $-1/4$. Reproducible local response and a shared limit on available energy alone do not supply that joint law. A global BAM history could in principle correlate the two records. To explain the principle it must produce both the exclusion of misses and doubles and a measure of histories giving the branch weights. This calculation establishes an obstruction for independent ports and a measured failure for one specified coupled apparatus.

[Primary experiment and controls](../experiments/closure_ledger/runs/20261125_two_detector_records/probe.json) · [Exploratory balanced-fluence sweep](../experiments/closure_ledger/runs/20261125_two_detector_records/fluence.json)

~~~sh
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.antipodal_network.two_detector_records
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.antipodal_network.balanced_record_fluence
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -p test_two_detector_records.py
~~~
