# Does an expanding loop preserve an early detector write?

At the strongest scanned expansion, 71 of 490 early single writes survived every sampled return, with the same detector retaining its record throughout the three-circuit readout window. An energy-preserving echo with the same travel times erased all 490. The registered primary rate erased all early writes. Nominal balanced illumination produced zero lasting singles at every tested rate. Expansion can protect an existing record in this apparatus. Producing one record from a balanced split remains open.

[Before-run protocol](expanding_detector_protocol.md) · [Static loop and latch definitions](closed_loop_records.md) · [Registered results](../experiments/closure_ledger/runs/20261125_expanding_loop/probe.json) · [Strong-rate follow-up](../experiments/closure_ledger/runs/20261125_expanding_loop/strong_rate.json)

## Propagation and the local energy account

The two comoving receivers sit a conformal travel time D=4 apart along each arm. Detector proper time obeys dt=a dη. Set a=1 through the pulse end at η=2.5 and then a²=1+2g(η−2.5)/D. The dimensionless parameter g controls the imposed expansion. This piecewise prescription gives every rate the same initial packet, detector masses, local potential, admittance and mutual inertia. At g=3 the scale factor reaches 1.80 by the leading edge of the first echo and 7.16 by the final readout. The effect therefore pertains to a guide that expands appreciably during the measurement.

A guided massless packet emitted at η returns after D in conformal time. Its proper arrival duration stretches by the scale-factor ratio a_rec/a_emit. Assuming that the guide captures the entire packet and conserves its wave action, its local energy falls by a_emit/a_rec. Since detector power is proportional to squared voltage, the returning voltage equals the emitted voltage times a_emit/a_rec. The timing-only control uses the square root of that voltage factor. It retains the pulse stretch and return time while preserving the integrated energy. [FLRW redshift and arrival intervals, §1.1.3](https://www.damtp.cam.ac.uk/user/tong/cosmo/cosmohtml/S1.html)

In this guide, local node work balances the difference of incoming and outgoing powers, and freely traveling energy falls as 1/a. On the closed ring these give d(E_receivers+E_wave)/dt=−H E_wave, where H is the proper-time expansion rate. The simulation sums the in-flight inventory, receiver energy and accumulated redshift term separately. The largest sampled balance residual over all reported runs was 7.28×10⁻¹² of nominal pulse energy. Monotone expansion preserves the upper bound from the initial energy at nominal illumination. The model prescribes the background and does not solve the Einstein equations or a dynamical throat-support problem.

## Record preservation

Each cell has 512 matched warm launches. The first-echo count uses the saddle at 0.75. A surviving write keeps the same exclusive basin at every sampled instant from the first echo through the end and meets the 0.85/0.65 margins over η=12 through 36. The durable count includes writes established at any time before the start of that window. It spans three full circuits in each geometry, including the slowest returns.

| g, right-biased input | First-echo singles | Same write survives | Durable right only | Mean initial energy fraction lost by end |
|---|---:|---:|---:|---:|
| 0 | 489 | 0 | 0 | 0% |
| 0.000001 | 489 | 0 | 0 | 0.00056% |
| 0.1 | 490 | 0 | 0 | 27.45% |
| 0.5 | 490 | 0 | 0 | 49.04% |
| 1 (primary) | 490 | 0 | 0 | 55.59% |
| 3 | 490 | 71 | 72 | 57.52% |

At g=1, all 4,608 nominal launches over nine angles were consistently ready in the readout window. The registered preservation criterion failed at this primary rate. All three selection gates also failed, including an unconditional right-only Born discrepancy of 1. Both meshes gave the same counts. The static limit reproduced the archived endpoint and persistent counts from #364 at every angle. The expanding and expanding-open reference traces agreed up to their first return.

## What survived the strong-rate controls?

The g=3 sentinel was on the registered grid. After observing the 71 survivors, we ran the full angular sweep and finer mesh at g=3, plus timing-only and open controls. Of the 490 early singles in the right-biased cell, the redshifted echo preserved 71 (14.49%), the timing-only echo preserved zero and open leads preserved all 490. Both closed controls gave identical counts when the maximum proper step was halved from .025 to .0125. The additional durable right-only trial brings the readout-window count to 72/512.

| Angle | Unconditional durable right only, g=3 | Same early write survives, either port |
|---|---:|---:|
| 0 | 72/512 | 71/512 |
| π/16 | 143/512 | 143/512 |
| π/8 | 0/512 | 0/512 |
| 3π/16 | 0/512 | 0/512 |
| π/4 | 0/512 | 0/512 |

The other four angles interchange the ports. The measured right-only frequency rises between 0 and π/16 while the Born target falls. Its maximum absolute error is 0.859375. At nominal balance all 512 pairs stayed ready. At 1.75× balanced illumination the strongly expanding ring retained zero singles and two doubles, with 498 pairs ready and 12 unsettled. The matched open leads retained 236 singles, 98 doubles and 178 misses. Thus some writes survive attenuation, but the unclicked branch has acquired no general mechanism that makes it dark whenever its partner fires.

## Replay

The PR workflow replays all categorical results exactly and budgets within a numerical tolerance. The tests check the null travel distance, pulse energy, and an entire asymmetric first echo computed by an independent proper-time ODE solver. All readouts describe a finite observation window. A spacetime boundary condition selecting complete MTY histories and a self-consistent cosmological solution both remain untested.

~~~sh
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.antipodal_network.expanding_detector_loop
OPENBLAS_NUM_THREADS=1 python -c 'from experiments.closure_ledger.antipodal_network.expanding_detector_loop import strong_rate_followup; print(strong_rate_followup())'
python -m unittest discover -s tests -p test_expanding_detector_loop.py
~~~
