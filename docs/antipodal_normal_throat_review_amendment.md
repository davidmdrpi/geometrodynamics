# Review amendment: separate angular reversal from normal reversal

Applies to PR #353 at 5bb26f2. Original freeze 533f925 and its first archive
remain unchanged. This amendment is published before its implementation.
Seed: 2026111202.

## Information already known

The reviewer has already run all four cells and strong-response amplitudes.
Their reported results are not novel blind predictions: shifts allow
nonzero total L_theta with either normal action; reflections force zero
total. Net normal transfer occurs for shift/normal-kept and
reflect/normal-reversed, and vanishes for the other two. Generic compatible
A data with nonzero mouth momenta also preserve the seam. Amplitudes
0.8,1.5,2.5 produce appreciable openings; the reported 2.5 opening is 1.06.
The existing tests and original 11-gate archive were independently reproduced.

We withdraw the causal reading "normal reversal permits circulation".
The replacement claim is compatibility, with the transfer obstruction as
the result attributable to normal reversal at fixed shift.

## Four controls fixed for the amended experiment

All maps advance phi by pi. Let c=+1 denote theta+pi, c=-1 denote -theta;
let r=+1 keep u and r=-1 reverse it. The spatial determinant is c*r.
Time remains preserved. The Hamiltonian and parameters are unchanged.
The previous "antipodal" and "reflective" APIs retain their meanings.

| c | r | total angular momentum on fixed data | face values |
|---|---|---|---|
| +1 | +1 | allowed | independent; transfer allowed |
| +1 | -1 | allowed | equal and constant; transfer forbidden |
| -1 | +1 | zero | each zero; transfer forbidden |
| -1 | -1 | zero | opposite; transfer allowed |

The derivation, not just the four samples, is required. If faces are swapped,
L_minus=c L_plus. If they are individually preserved, L_s=c L_s.
Together with total conservation these give the table. In this family the
two non-orientable cells forbid net normal transfer. This is not a theorem
about arbitrary non-orientable manifolds, boundary conditions or generators.

## Revised gates

1. Actual coordinate map: finite-difference Jacobians (wrapped angle
   differences), involution, and reconstructed real-field pullback against
   the state deck, for all four cells. Error limits: Jacobian 1e-8,
   coordinate involution/pullback 1e-10. Test H invariance and canonicality.
2. Keep spectral checks. Replace the susceptibility eigenchannel-only probe
   by a one-face source, retaining both measured face responses against the
   independent full-grid solve, error <1e-10.
3. Rename continuum_convergence to phi_symbol_convergence. Test m=1,2,3,
   n=1,2 at N=16,24,32,48; order 1.8--2.2. This makes no claim about normal
   continuum refinement or nonlinear continuum convergence.
4. Retain the 32 original responding runs, but tighten their sampled relative
   energy threshold from 2e-3 to 2e-5. No frozen threshold is rewritten.
5. Gated generic runs: (N,n)=(16,1),(24,2), four cells, eta=+/-1,
   four draws per eta/cell/configuration. Project independently drawn field
   positions and momenta; normalize their joint norm to 0.4. Set nonzero
   random opening positions and momenta (scales .03 and .07), then project.
   Integrate [0,40] with h=.005, record every 20 steps. Require >10% of field
   power outside the strongest phi Fourier bin, nonzero initial mouth
   momenta, seam/H-symmetry errors <1e-9, normalized global L drift <1e-9,
   and sampled relative energy excursion <5e-4. This energy bound accounts
   for the resolved high-frequency content at this timestep.
   Normalize quadratic momentum/flux diagnostics by the initial squared
   field-phase-space norm, including zero-total cases.
6. On generic runs, require forbidden total/face variation/net flux <1e-9;
   allowed total has a resolved witness >1e-4; allowed transfer has face
   variation >1e-5 in every draw. Test the face Noether derivative against
   the explicit normal flux on arbitrary non-compatible states.
   Inject broken-deck and frozen-flow data into validation to demonstrate
   that seam and non-vacuity gates reject them.
7. Strong response: A, N=16,n=1,m=0, amplitude=2.5, q=pq=0, [0,40].
   Require opening >0.5, seam and relative L error <1e-9, sampled energy
   excursion <5e-5. Compare h=.01,.005 with independent DOP853
   (rtol=1e-10,atol=1e-12); both final-state errors <1e-2 and halving improves.
   The A-compatible opening is ONE independent scalar coordinate.
   The normal-kept cells can retain two independent opening coordinates.

The general coupled symmetry still implies preservation in exact arithmetic;
these are implementation and discretization validations, not independent
evidence for the symmetry theorem. Generic data and negative controls make
the coverage and ways those gates can fail explicit.

## Reporting

Archive amended results separately, link this public amendment, update the
title, README, changelog and write-up around the four-cell causal comparison.
Report original and amended thresholds distinctly. Keep unknown charge,
spin, external-mouth geometry and Einstein dynamics underived. The original
25 tests, first data archive and preregistration remain part of the history.
