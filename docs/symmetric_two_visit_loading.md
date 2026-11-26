# Can two identical detector loads protect complementary ranks?

Equal intensity loads at the two visits break the rank sum even when **both visits are free**. A different reciprocal coupling to a shared pointer has an exact protected regime, but its pointer then has zero displacement on every selected history. Above its response threshold, the pointer moves and complementarity fails on its invariant ensemble. Two identical pointers with specified ready states provide a third comparison in which both ranks and both readouts vary together.

The positive result from [#367](antipodal_two_time.md) is the variational selection of complementary ranks for the unloaded minimum. An invariant ensemble supplies a uniform rank and comparison with the mixer's incident energy supplies the conditional Born curve. The distinction matters because neither an unconditional probability law nor detector backreaction follows from the two-pass sign. Single-qubit classical probability representations have a long history, including the Kochen–Specker model reviewed by [Leifer](https://arxiv.org/abs/1409.1570).

[Before-run protocol](symmetric_two_visit_protocol.md) · [Code](../experiments/closure_ledger/antipodal_network/symmetric_two_visit.py) · [All path fits and predictions](../experiments/closure_ledger/runs/20261125_symmetric_two_visit/probe.json)

## The equal local loads

Use #367's unit-time Berger kinetic action S with λ=2. Write x=X(0), y=X(1/2), X(1)=−x, u=|x₁|² and v=|y₁|². Every path obeys S≥A=π²/2, and every unloaded minimizer has u+v=1. The loaded action is S−κ(u+v). The power setting p does not appear. Under the two barrier probes, u+v>1 opens an interval of misses, and u+v<1 opens an interval of doubles. Its length in p is |u+v−1|, with any separatrix endpoints handled separately.

**Pinned first visit.** The extra term −κu is constant. The entire set of minimizing paths is therefore identical to #367's one-sided load at the same x. For any interior initial rank, that earlier trial proof gives κ(u+v−1)>0 on every minimum.

**Floating first visit.** Allow x to vary along with every interior node, preserving antipodal closure. Every complementary path still costs at least A−κ. For a trial path in the real horizontal plane, take its visit angles to be α=π/4−δ/2 and α+π/2+δ, and finish at α+π. Use constant angular speed on each half of the path. Its cost improvement over A−κ is exactly

$$
 G(\delta)=\kappa\sin\delta-2\delta^2 .
$$

Indeed both visited ranks are (1+sin δ)/2. For any κ≠0 some sufficiently small δ of the same sign gives G>0. A global minimizer exists by coercivity and compactness. The bare bound S≥A implies κ(u+v−1)≥G for **every** global minimizer, even outside the real trial family. This proves the sign of the defect without asserting a distribution for the loaded paths.

| Positive load | Pinned guaranteed defect, u=1/2 | Pinned fit | Floating guaranteed defect | Floating fit |
|---:|---:|---:|---:|---:|
| 0.1 | ≥0.0124895 | 0.02496 | ≥0.0124973 | 0.02499 |
| 0.3 | ≥0.0372190 | 0.07390 | ≥0.0374297 | 0.07472 |

The displayed bounds use δ=κ/4. Reversing the sign reverses the measured defect; exchanging the component labels relates the two actions up to a constant. Equal treatment of visits removes a port preference while allowing both ranks to move to the same side of the complementary line.

## A shared pointer with an exact response threshold

Couple an auxiliary real coordinate Q to the overlap c=Re〈x,y〉 by adding kQ²/2−gQc with k>0. Both Q and the path are varied. This observable is invariant under common U(2) rotations and interchanging x and y. Relabeling the origin by half a circuit on the antipodal lift sends (x,y) to (y,−x), flipping c. In that convention the signed pointer flips with c, and the reduced coupling is invariant. The pointer reads a two-visit correlation. An equilibrated displacement is not by itself a which-port record.

Eliminating Q gives Q=gc/k and the load −χc²/2 with χ=g²/k. For fixed real overlap let d=arccos c. The minimum round lengths of the half paths are d and π−d. Comparing the Berger and round metrics on each half gives the bare lower bound d²+(π−d)². For each c in the interior, a real horizontal plane realizes the bound with constant speed on each half. Equality requires horizontal shortest arcs. In particular 〈ix,y〉=0 and the complex overlap has modulus |c|. With δ=arcsin |c| the **exact** reduced global problem is therefore

$$
 \min_{0\le\delta\le\pi/2}\left[A+2\delta^2-\tfrac12\chi\sin^2\delta\right].
$$

For χ≤4, the unique minimizing overlap is zero, including at the quartically flat critical point. The original complementary ranks survive and Q=0 on every history. For χ>4, the unique minimizing magnitude solves sin(2δ)/(2δ)=4/χ. Two signs of c occur, with nonzero pointer magnitude. They need no chosen relative weights for the following statistics, which agree on either orbit.

Under the common U(2)-invariant orbit measure, the eigenvalues of xx*+yy* are 1±|c| and u+v is its first diagonal entry. Its distribution is uniform between those eigenvalues. Thus the mean defect is |c|/2. The balanced double and miss probabilities are each δ/π by the hemisphere calculation of #367. Each marginal rank remains uniform, illustrating why looking only at the right detector's total click frequency would conceal the doubles. These are eventual barrier-probe outcomes.

| χ | Equilibrated pointer magnitude (k=1) | Mean rank defect | Balanced doubles | Balanced misses |
|---:|---:|---:|---:|---:|
| 2 | 0 | 0 | 0 | 0 |
| 4 | 0 | 0 | 0 | 0 |
| 4.1 | 0.38635 | 0.09540 | 6.111% | 6.111% |
| 6 | 1.66588 | 0.34005 | 23.806% | 23.806% |
| 10 | 2.76275 | 0.43683 | 33.826% | 33.826% |

## Follow-up with two independently readable, identical pointers

After the registered scan, we gave each visit its own unit-mass pointer. Each starts at zero, has one local unit of time, and pays kinetic action plus a terminal spring cost (Qᵢ−bᵢ)²/2, with bᵢ=u or v. These auxiliary pointer paths and geometry minimize a joint positive functional. This is an endpoint interaction in a variational model. It does not model an irreversible laboratory amplifier.

For fixed final Qᵢ, the minimizing pointer path is straight, with kinetic cost Qᵢ²/2. Varying Qᵢ then gives Qᵢ=bᵢ/2. Back on the geometry the two pointer loads contribute (u²+v²)/4. Every complementary history now costs at least A+1/8, attained within the complementary class at u=v=1/2. The same floating real trial at δ=−1/16 has cost A+2δ²+(1+sin δ)²/8, which improves that bound by G=0.00731468. For any loaded minimizer, (u+v)²≤2(u²+v²)≤1−8G. **The guaranteed defect 1−u−v is at least 0.02969.** The numerical fits give u≈v≈0.47065 and both final readings ≈0.23533. Their implied double interval includes balanced illumination in those fits. No invariant probability is assigned to this loaded minimum set.

The ready-state constraint is material: if the two initial pointer positions can instead float, the kinetic and spring costs can both vanish on constant paths Qᵢ=bᵢ. That construction simply starts the pointers with the answers.

## Numerical checks and interpretation

The solver varies eight segments in four real coordinates, imposes the antipodal endpoint, and uses #367's exact action of each round-arc interpolant. It varies the shared pointer independently in all seven susceptibility cells and varies both terminal readings independently in the follow-up. Three seeds are archived for each case. At χ=4 the optimizer returns small residual overlaps (up to 0.00042 in the archive) along the flat direction. The analytic lower envelope establishes zero there. Tests differentiate the arc action independently, reconstruct the trial certificates from paths, compare numerical actions and pointer readings against the predictions, and sample 150,000 invariant histories.

The tested local loads break complementarity. The tested shared pointer protects it at weak coupling while providing no displacement to read. This establishes a limitation of these particular couplings. The unloaded theorem remains intact. A detector that retains which-port information and preserves the correlation under load is still missing.

~~~sh
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.antipodal_network.symmetric_two_visit
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -p test_symmetric_two_visit.py
~~~
