# Response to review of PR #347 at 03bd7ca

[Review](https://github.com/davidmdrpi/geometrodynamics/pull/347#issuecomment-6331643198)

| Item | Resolution |
| --- | --- |
| P1: caustic phase | Complex coefficients carry −i at every empty forward antipodal focus, and +i backward. Open-mouth capture before/after this phase is explicit boundary data. Two forward empty foci give −1. Checked against the integrated conformal field and `pair_birth.bare_timeline`. Closed boundaries independently select Dirichlet −1 or Neumann +1. |
| P2-1: chart-dependent topology | Explicit `handle='klein'|'torus'` replaces the ambiguous parity argument. Chart matrices are derived from topology and both chart conventions. Actual rim-column and orientability tests cover both handles and all four chart choices. |
| P2-2: loop, gain, time direction | With post-caustic capture, the old route is explicitly a closed null route imposed by the clock data plus channel-forward resets. Its new gain is −1, its unit-seed mismatch is 2, and only zero solves the homogeneous front condition. Default transport now carries time orientation through both attachments; a paired control shows that the same clocks then do not close the four-leg loop. Neither case is a solved driven field/recoil problem. |
| P2-3: aperture interior | Disc-to-disc transport is explicitly identified as a portal extension in the API documentation, exported data, replay, and report. The mesh test validates the rim only. |
| P3-4: absolute manifest key | Generator paths resolve from the repository root; every manifest key is relative to it. A test copies the complete manifest file set to another root and verifies all hashes there. |
| P3-5: rounded recurrence | Rounded keys are replaced by a stated tolerance on time and spatial return, with channel, endpoint and time direction matched. Closest-return residuals are retained even when the leg budget is reached. Closure is checked after the last allowed leg too. |
| P3-6: archive size | Deliberately retained. The user requested publication and preservation of the negative result. The optional size reduction is deferred; the original commit pin and historical labels remain. |

Validation: **44 tests pass** (32 transport/reproduction and 12 historical-action
regression tests). The eight controls reproduce the updated table. Canvas QA
exercises 638 sampled frames, selection, stepping, rotation, cutaway and playback
termination, with no nonfinite draw coordinates. The reviewer independently
verified headless Chromium playback without network access at `03bd7ca`; that
result is attributed to the reviewed revision.

The [updated report](antipodal_throat_transport.md) states the phase convention,
portal assumption, time policies, residuals and remaining dynamics limitations.


## Later follow-ups: ad746d6 and 0380308

Inspected the pinned #345 source and its §10 report. The latest comment corrects
the earlier all-angular-momentum-reversal claim. Under that specific tube rule:

| Component in the stated antipodal birth frame | Klein | Torus |
| --- | --- | --- |
| Orbital contribution along the birth circle | Reversed | Reversed |
| Orbital contribution from sideways motion | Kept | Reversed |
| Spin about the mouth axis | Reversed | Kept |
| Assumed compression sign | Flipped | Kept |

Its reported attractive Klein exchange is too weak to bind (about 3–8% of the
needed symmetric stopping change). Those runs were not repeated here, and the
stated capture/recoil rules are not treated as a field derivation.

The source comparison found a further distinction: #345 applies caustic phase
only on empty-focus continuation, before testing subsequent captures. #347 now
makes pre-/post-caustic capture explicit. The paired imposed-loop controls have
gains +1 and −1 respectively, demonstrating why that boundary choice matters.
The default carries time orientation; the channel-forward reset remains an
explicit control. The next static-flux exchange proposal is documented as
unbuilt here, without silently adding continuous attraction or an ad hoc impulse.
