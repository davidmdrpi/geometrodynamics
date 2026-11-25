# Which two-detector histories can produce one record at every setting?

In this classical boundary model, complementary detector thresholds are necessary and sufficient for exact single records across the full range of power splits. A norm-preserving exchange between two quadrature pairs realizes this correlation. With a fixed, uniform disk preparation, the areas of its two basins are exactly cos²θ and sin²θ. The same preparation is used for every angle. Every trial counts.

[Lean proofs](ClosedHistory.lean) · [Two open receivers (#363)](https://github.com/davidmdrpi/geometrodynamics/pull/363) · [Recirculating receivers (#364)](https://github.com/davidmdrpi/geometrodynamics/pull/364)

## The operational obstruction

Let a complete trial carry threshold ranks u and v in [0,1]. Supply power fractions p and 1−p at the two ports. The right port fires if u < p; the left port fires if v ≤ 1−p. The opposite strictness assigns a threshold tie to the left on a complementary pair. It affects pointwise statements and has zero probability under the disk preparation.

For this response law the switches occur at p = u and p = 1−v. Their ordering gives the full classification.

| Threshold sum | What happens between the switches | Length in power-split space |
|---|---|---:|
| u + v < 1 | Both detectors fire | 1−u−v |
| u + v = 1 | Exactly one fires at every setting | 0 |
| u + v > 1 | Both detectors miss | u+v−1 |

Lean identifies the entire failure set as (min(u,1−v), max(u,1−v)] and its Lebesgue length as |u+v−1|. This gives a quantitative diagnostic for a proposed physical return map. The length uses uniform weight in power split p. An ensemble uniform in angle has a different weighting.

The classification also handles the operational order of quantifiers. Assume the same distribution of histories for every setting and zero probability of failure at each fixed setting. Then u+v=1 almost surely. For every history off the anti-diagonal the failure interval contains a rational setting. Only countably many null sets need be discarded. The proof applies to any joint preparation law on bounded ranks. In particular, broadening two independent continuous readiness distributions cannot produce this coupling. The one-dimensional support in the joint threshold plane is a physical constraint on any realization in this class.

## A one-pass exchange with antipodal two-pass closure

Take a normalized pair of planar modes (z,w) with |z|²+|w|²=1. At each port the detector samples the squared radius of the **first** mode. Between the two visits apply the same exchange

    J(z,w) = (−w,z),      J²(z,w) = (−z,−w).

A history closes after two passes on identifying antipodal states. The right visit sees u=|z|² and the left visit sees v=|w|². The normalized budget therefore supplies the required rank reversal v=1−u. In the two-rank description the complete-boundary map is (u,v) ↦ (1−v,1−u), whose fixed points are exactly these pairs.

The control K(z,w)=(iz,iw) also squares to the antipode. It preserves the sampled rank at both visits. At the normalized state (0,1), balanced illumination gives two clicks. This distinguishes an exchange of modes from the two-pass sign alone. These maps act on classical real quadratures written in complex notation. Connecting either map to a particular throat holonomy requires a physical transport law.

## Why the basin has the Born area

The prepared coordinate z ranges over the open unit disk with normalized planar Lebesgue area. The proof lifts it to the normalized pair (z,√(1−|z|²)); the phase of the second mode drops out of the readouts. With the supplied mixer angle, the right test is equivalent to |z| < |cos θ|. Thus its basin is a concentric disk and its relative area is cos²θ. The left basin is the complementary annulus and has relative area sin²θ. The proof checks the set equalities and both normalized measures, including the endpoint settings. No probability of a detector click appears among the hypotheses of the witness.

There is also a geometric interpretation of the chosen measure. Write a point of the unit round three-sphere as (r exp(iα),√(1−r²) exp(iβ)). In the interior the induced metric in (r,α,β) is diagonal with entries 1/(1−r²), r² and 1−r². The volume element is r dr dα dβ. Integrating out β gives uniform area in the first disk. This coordinate calculation explains the relation to round-shell volume. The Lean witness starts with planar area and checks its lift and basin measures. The identification of a physical BAM preparation with round-shell volume is an additional hypothesis.

## What a throat calculation would have to show

The toy specifies the conserved two-mode budget, the one-pass exchange, the mixer response, the threshold readout and the preparation measure. The discriminator for a candidate throat action is to compute the paired return ranks before assigning click frequencies. Any nonzero |u+v−1| opens a calculable interval of misses or doubles. Even exact complementarity leaves the marginal rank distribution to be determined. On complementary histories the right frequency is precisely Pr(u<p); recovering p for every split requires uniform rank.

This is a consistency result for records assigned to complete histories. Energy transfer to a detector, latching, persistence over repeated circuits, and the boundary-value selection of histories for the apparatus in #364 are separate physical problems. There is no inference here of multi-setting Bell correlations or of a general failure of classical dynamics.

## Reproduce

The dedicated PR check uses Lean 4.19.0 and mathlib commit c44e0c8ee63ca166450922a373c7409c5d26b00b. From a checkout of mathlib at that commit, with this repository in a sibling directory named geometrodynamics, run

~~~sh
lake exe cache get Mathlib/MeasureTheory/Measure/Lebesgue/VolumeOfBalls.lean Mathlib/Tactic.lean
lake env lean -DwarningAsError=true ../geometrodynamics/formal/closed_history/ClosedHistory.lean
~~~

The file prints the axiom dependencies of the main results. The build treats warnings as errors and rejects a proof log containing sorryAx.
