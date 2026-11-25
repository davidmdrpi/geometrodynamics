# What happens to a detector record when the wave comes back?

The shared loop returned energy at the predicted times and erased writes that survived on open leads. At the strongest right-biased setting, 489 of 512 primary runs had one detector above the saddle at the first transit. On open leads those 489 remained written; on the loop the single-write count fell to zero by the third transit. Across all nine primary angles, **zero of 4,608 launches kept one detector written and the other ready throughout the registered readout window**. The independent sampling repeat also produced zero persistent single records. [Pre-run protocol](closed_loop_records_protocol.md) · [Open detector experiment](two_detector_records.md)

## One apparatus, two return paths

Place identical #363 receivers at opposite nodes of a two-arm ring. Both arms have characteristic admittance $\gamma/2$ and one-way travel time $D$. In transmission-line units the node voltage is $v_i=\dot q_i$. A shunt capacitance of one at each node, a mutual capacitance $\mu=0.25$, and nonlinear restoring energy $U(q)=6q^4-14q^3+9q^2$ implement the same passive endpoint equations as #363. Here $q_i$ is an integrated node voltage (flux coordinate). Initial receptor flux and traveling voltage/current fields specify the circuit state. The lumped mutual inertia is present in both matched geometries.

Prepare the two arms symmetrically, so each carries half the incoming energy at a port. On $0\le t<D$, the incoming wave at each port is exactly the #363 incident waveform, scaled by the square root of that port's mixer weight, and zero after the pulse of duration 2.5. These first-transit histories specify the initial wave packets already on the ring. No source replenishes them. With $f_i$ and $b_i$ the incident and outgoing voltage waves, respectively,

$$
 M\dot v=2\gamma f-\gamma v-\nabla U_{\rm pair},\quad b_i=v_i-f_i,
 \qquad f_R(t+D)=b_L(t),\quad f_L(t+D)=b_R(t).
$$

The two equal arms give two identical incident signals at each port. This reduction retains both propagation directions and the phase of the reflected signal. The sum of the two arm admittances is $\gamma$. The open control starts with exactly the same packets and sends outgoing waves away. Before time $D$ its detector solution therefore coincides with the ring solution. This setup tests causal recirculation in a closed classical wave circuit. A two-time boundary condition on complete histories would be a different experiment.

The node work rate is $\gamma\sum_i(f_i^2-b_i^2)$. The line energy at $t$ equals $\gamma\sum_i\int_t^{t+D}f_i(s)^2ds$, since these future incident waves are precisely the packets currently in flight. Its derivative cancels the node work rate. Thus the primary #363 total-energy cap against simultaneous occupation of the two fired wells still applies. Closure allows radiation to return and push a previously written detector back over its saddle.

## How a persistent record is counted

The protocol fixes the common primary delay $D=4$, #363's bounded readiness ensemble ($e=.02$), all nine mixer angles, and every prepared trial in the denominator. From $t=12$ through $36$, a written port must stay above 0.85 and a ready port below 0.65. The interval spans three full circuits at the primary length and at least two at every control length. A trajectory that fails the margins enters the changing/unsettled category. Final-time basin counts use the single saddle at 0.75. A final-time crossing need not be a lasting record.

| Angle | Matched open right only | Loop right only | Loop left only | Loop neither | Loop changing/unsettled |
|---|---:|---:|---:|---:|---:|
| 0 | 489 | 0 | 0 | 512 | 0 |
| π/16 | 358 | 0 | 0 | 487 | 25 |
| π/8 | 4 | 0 | 0 | 510 | 2 |
| 3π/16 | 0 | 0 | 0 | 512 | 0 |
| π/4 | 0 | 0 | 0 | 512 | 0 |

There were zero persistent doubles. The other four angles exchange the ports. At balanced illumination the ring exclusive fraction and the open exclusive fraction were both zero. The largest unconditional right-only discrepancy from $\cos^2\theta$ was **1.0**, at $\theta=0$. All three registered physical gates failed on both meshes and both preparation samples. These counts describe this pulse, initial wave phase, apparatus and finite readout interval. They give no universal half ceiling for coupled systems.

## The return can be seen directly

The balanced warm sample's tracked trajectory agreed exactly with its matched open trajectory up through $D$. The first sampled position difference above $10^{-8}$ occurred at 3.025, 4.025 and 6.025 for $D=3,4,6$, respectively, one time step after each physical arrival. At every length all 512 balanced preparations stayed in the ready margins over the readout interval. The exactly symmetric cold control remained symmetric; the dark pair stayed at rest. Removing the mutual inertia still gave 512 persistent misses at balance.

The following **instantaneous** exclusive counts show the fate of the strong right-biased write at $D=4$.

| Time | Closed loop | Matched open |
|---:|---:|---:|
| 2.5, pulse ends | 501 | 501 |
| 4, first transit | 489 | 489 |
| 8, first full circuit | 492 | 489 |
| 12, readout starts | 0 | 489 |
| 24 | 0 | 489 |
| 36 | 0 | 489 |

All counts above use the primary grid. On the finer open grid the endpoint setting gives 490 writes, reproducing #363's count. The archived closed-loop persistent counts are identical on the two grids. The second sample (1,024 launches per angle) even has a few one-detector basin occupancies at $t=36$, while the persistent single count is still zero.

**Symmetry and finite readout.** With equal inputs and exchange-symmetric initial data, uniqueness gives $q_R(t)=q_L(t)$ at every time. The two positive record margins require a nonzero separation, so by continuous dependence an open neighborhood of that history also fails to give one record at the specified time. A preparation with positive measure on that neighborhood cannot give perfect single-record efficiency. This argument leaves the size of the missed fraction open; the measured zero persistent success is a stronger, model-specific outcome.

## Exploratory brighter-pulse control

After reading the primary result, we repeated #363's balanced fluence grid from $1$ to $2.25$ times its nominal energy in steps of $0.125$, at all three loop lengths. At $1.75$ times the nominal energy, open leads reproduced #363's **118 right only, 118 left only, 98 both and 178 neither**. The loop gave zero persistent singles at all three lengths. At $D=4$, 18 trials had exactly one detector above the saddle at $t=36$, while over the full window 186 were consistently ready and 326 were changing/unsettled. No loop setting on this exploratory grid exceeded the half ceiling. This control compares finite sampled fluences and does not optimize the circuit.

## Numerical checks and replay

The step stores one transit of midpoint wave samples and delivers the outgoing sample to the opposite port after exactly $D/h$ steps. The node update uses midpoint velocity $z=(v_1+v_0)/2$ and the divided difference of $U$ between $q_0$ and $q_1=q_0+hz$. Multiplying the update by $hz$ gives the mechanical energy increment $h\gamma\sum(f_i^2-b_i^2)$ up to the implicit solve residual. At each transit and the terminal time the implementation separately sums the wave energy still in the buffer. It adds escaped energy in the open control. The worst registered residual was $7.23\times10^{-13}$ of the nominal pulse energy. Midpoint quadrature gives initial incident energy 1.542845854 at $h=.025$ and 1.542845830 at $h=.0125$ (continuum 1.542845829).

Persistence is evaluated on the discrete trajectory. The archive also reports a more conservative sampled-margin diagnostic enlarged by $h\sqrt{2E_{\rm initial}}/2$, using the continuum speed bound; it gives the same zero persistent singles. This is a sensitivity check, not a rigorous enclosure of the numerical solution. A separate adaptive method-of-steps test integrates a warm asymmetric launch through the entire first echo and checks convergence of the buffered trajectory. No inference of infinite-time memory is made.

[All registered runs, counts, traces and gates](../experiments/closure_ledger/runs/20261125_closed_loop_records/probe.json) · [Brighter-pulse grid](../experiments/closure_ledger/runs/20261125_closed_loop_records/fluence.json)

~~~sh
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.antipodal_network.closed_loop_records
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.antipodal_network.closed_loop_fluence
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -p test_closed_loop_records.py
~~~
