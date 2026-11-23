# Does the laboratory energy flow shift the branch frequencies?

Yes. On the two-winding #359 packet with the original zero-offset Stern–Gerlach drive, the maximum absolute shift of the right-record probability from its same-packet nonrelativistic limit is approximately $0.065/c^2$. At $c=8$, it exceeds #358's 0.00049 window in both the continuous probability and the 2,048-position count. The conserved scalar number current gives a different answer on exactly the same fields.

## From a velocity to a recorded branch

Use #359's positive-free-frequency $k=\pm1$ initial data, #241's mixer at ellipticity 0.02, all nine angles from #358, and its sign-of-position readout at $t=20$. Set the initial position density proportional to the density of each candidate current. With $\rho=\sum|u_k|^2$ and $J=\sum\mathrm{Im}(u_k^*\partial_xu_k)$, the three comparisons are $J/\rho$, $J/n$ and $S/E$. Here $n=-\sum\mathrm{Im}(u_k^*D_tu_k)/c^2$ and $(E,S)$ are the fiber-integrated densities and fluxes of the laboratory energy defined in the [stress derivation](timelike_current.md). All angles use the same two evolved component fields with their mixer weights. No apparatus distribution is fitted.

In one spatial dimension, trajectories preserve their ordering wherever the velocity is regular. Integrate the trajectory ending at $x=0$ backwards with RK4 on the evolved KG fields, then integrate the initial density above its starting point. The domain ends in the packet tails at $|x|=56$. Each density is integrated separately, using cubic interpolation of fields before forming velocities. For a continuous right probability $P$, ordering gives the exact count on #358's deterministic midpoint ranks by rounding $2048P$ to the nearest integer. These are deterministic numerical counts, with one-count resolution $1/2048$.

The reference $P_{\rm NR}=q+(1-2q)w_+$ uses the actual two-channel mixer weight and the single-branch wrong-side Gaussian tail $q=0.0002092$. We also archive the original #358 full-ring limit; retaining only the pair changes the continuous reference by at most 0.0000676 across angles. The #358 value 0.00048828125 was its maximum counted discrepancy from $\cos^2\theta$. Its continuous discrepancy was 0.0002768.

## Measured shifts

Entries are maxima over the nine angles. The series uses 1,024 spatial sites and time step 0.02. $P_E$ denotes the transported energy-weighted **position** ensemble.

| $c$ | $\max\lvert P_E-P_{\rm NR}\rvert$ | $c^2\max\lvert P_E-P_{\rm NR}\rvert$ | Maximum continuous error from $\cos^2\theta$ | Maximum counted error from $\cos^2\theta$ |
|---:|---:|---:|---:|---:|
| 8 | 0.00100759 | 0.06449 | 0.00089902 | 0.00104110 |
| 12 | 0.00045096 | 0.06494 | 0.00034238 | 0.00055282 |
| 16 | 0.00025339 | 0.06487 | 0.00020900 | 0.00006453 |
| 24 | 0.00011266 | 0.06489 | 0.00020913 | 0.00006453 |

At $c=12$, quantile counting puts the result just outside the window although the underlying probability is inside it. Doubling the spatial resolution and halving the time step changes $P_E$ by at most $5.38\times10^{-7}$ at $c=8$ and $1.99\times10^{-7}$ at $c=24$. At the refined $c=8$, $\theta=3\pi/16$, $\cos^2\theta=0.691342$, $P_{\rm NR}=0.691233$, and $P_E=0.692241$. This gives 1,418 right records. Both the scalar number flow and #358 give 1,416. The shift changes sign at the reflected angle; the balanced packet stays at one half for zero offset.

## Independent checks on the interpretation

The scalar phase current obeys $\partial_t n+\partial_xJ=0$ on the time-dependent metric. Its transported position fraction and the directly integrated final number fraction agree within $4.51\times10^{-7}$ on the refined $c=8$ run. Its position fraction differs from the same-pair nonrelativistic reference by at most $5.70\times10^{-7}$ on that run. The initialized complex packet supplies positive number density on the measured support. A Noether current for global $U(1)$ phase symmetry requires a complex field (or an equivalent real doublet).

During the drive the laboratory energy obeys $\partial_t E+\partial_x S=Q\partial_t A$. The energy of trajectories that end on the right therefore changes. Integrate the work over the moving region to the right of the energy characteristic. The largest residual in final right energy minus transported initial right energy minus this work, normalized to total initial energy, is $2.64\times10^{-6}$ on the primary $c=8$ grid and $8.07\times10^{-7}$ after refinement. At the example angle the final right **energy fraction** is 0.691234. Simply integrating final energy would miss almost all of the shift in the recorded positions.

Changing the constant fiber shift from zero to 0.1 changes the choice of laboratory time generator while preserving the force. At $c=8$ the maximum energy-guided shift grows to 0.0016711. The number-guided probabilities are invariant to numerical precision. Taken together, the fixed-drive scaling and the generator control support a finite relativistic correction to this particular energy weighting. The coefficient 0.065 also depends on packet shape, force history and readout time. The experiment does not derive a universal coefficient linear in the local potential. Determining whether a physical throat follows either current requires a model of its motion.

~~~sh
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.antipodal_network.branch_current_shift
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -p test_timelike_current.py
~~~

[Per-angle probabilities, counts, energy budget and refinements](../experiments/closure_ledger/runs/20261123_equivariance/probe.json)
