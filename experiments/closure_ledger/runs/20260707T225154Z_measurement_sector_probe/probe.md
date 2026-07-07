# The measurement sector: pointer outcomes for the entangled sector - companion probe (PR #209)

**Run:** 2026-07-07T22:51:54+00:00

The deliverable is `docs/measurement_sector.md` - the measurement chain for the entangled sector: coupling, equivariance, Born, and the operational Bell test. *(QFT on the fixed classical throat geometry, not quantum gravity.)*

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the last open: internal state -> pointer -> Born | **PASS** |
| T2 | `T2_coupling_derived` | the KK connection gradient is the pointer coupling | **PASS** |
| T3 | `T3_fiber_integrated_equivariance` | fiber-integrated equivariance (the #198 extension) | **PASS** |
| T4 | `T4_born_for_internal_states` | Born for internal states: P = cos^2 to 0.005; permanence | **PASS** |
| T5 | `T5_operational_bell` | CHSH 2.82 from beable positions; marginals flat | **PASS** |
| T6 | `T6_spatial_sector` | the spatial sector: pointer + EPR from nucleation | **PASS** |
| T7 | `T7_honest_scope` | honest scope | **PASS** |
| T8 | `T8_assessment` | the thread closes operationally | **PASS** |

## Born statistics for internal states (the SG sweep)

| beta | P(+) measured | cos^2(beta) | diff |
|---|---:|---:|---:|
| pi/8 | 0.8516 | 0.8536 | -0.002 |
| pi/6 | 0.7477 | 0.75 | -0.0024 |
| pi/4 | 0.4974 | 0.5 | -0.0027 |
| pi/3 | 0.25 | 0.25 | -0.0 |

## The operational Bell test

- sanity E(0,0) = -1.0
- correlators: {'(0.00,0.79)': -0.7019, '(0.00,-0.79)': -0.7019, '(1.57,0.79)': -0.7019, '(1.57,-0.79)': 0.7187}
- **CHSH = 2.8244** (Tsirelson 2.8284); marginal shift 0.0 vs noise 0.01

## Verdict

**THE_MEASUREMENT_CHAIN_CLOSES_THE_KK_COUPLING_MAKES_THE_POINTER_FIBER_INTEGRATED_EQUIVARIANCE_MAKES_BORN_OPERATIONAL_CHSH_2SQRT2_FROM_BEABLE_POSITIONS.** THE LAST OPEN OF THE ENTANGLED-SECTOR THREAD CLOSES (the argument is in docs/measurement_sector.md; this probe measures each link).

THE COUPLING. The committed structure already contains the pointer device: winding couples to the fiber connection (the KK gauge coupling), whose gradient exerts a k-odd force (dispersion identity to 2e-16; live: opposite windings deflected to opposite sides). The Stern-Gerlach is charge measurement by deflection - derived, not postulated.

THE STATISTICS. Fiber-integrated equivariance (the #198 theorem extended; residual 1e-04, Born ensemble at noise through branch separation) delivers Born statistics for INTERNAL states: P(+) = cos^2(beta) to 0.0027 across the sweep, with the empty branch's influence dying with the branch overlap (1e-04 -> 6e-09: effective collapse from geometry, without collapse).

THE CLOSING. The #206 bridge singlet + local rotations + Stern-Gerlach branches + dBB beables: E(0,0) = -1.0, CHSH = 2.8244 from POINTER POSITIONS, marginals setting-independent (0.0 vs noise 0.01). Bell violation as classical-beable statistics, no imports anywhere in the chain.

THE SPATIAL SECTOR. The pointer IS the spatial sector; positional EPR follows from conservation at nucleation (Duan-Simon 0.5313 < 2); and the #205 guiding-without-gravitating split is realized exactly here. Remaining program-wide: the 5D pants nucleation, W-class reachability, the strong-field NR target.
