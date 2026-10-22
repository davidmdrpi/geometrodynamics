# PR #347: antipodal-only transport, phase and topology

The earlier local scalar-port action coupled fields to mouths away from
antipodal arrivals. Its trajectory changes are **not evidence for the corrected
BAM interaction rule**. The current `geometrodynamics.waves.history_throat`
entry point exports `AntipodalTransport`; the old action is isolated in
`legacy_history_throat.py` for reproduction.

The [replay](images/history_throat/replay.html) follows one front through the
existing connected Klein mesh, showing outer/inner attachments, surface focusing,
and throat passage. The curved handle is an immersion drawing of an
extradimensional connection. Its visible arc length and wall crossing are not
physical propagation distances.

![Connected-mouth transport](images/history_throat/antipodal_transport.png)

## Interaction support and complex phase

For a front launched at (n,t), the next surface interaction is at
**(−n,t+dπ)**, where d=±1 is its exterior time direction on the unit-sphere
reduction. No mouth is tested or driven between these endpoints. At arrival:

- An open mouth containing the focus admits it into the throat, using the
  mouth's position at **arrival**, not its earlier position.
- An empty antipode refocuses the same front onward. No new source is injected.
- A closed mouth at the focus reflects it there. A closed throat outlet
  reflects it inside the throat.

The previous real −1 per focus was wrong. The positive-frequency conformal
sphere obeys

    exp[−i(l+1/2)π] = (−i)(−1)^l.

The spatial antipodal map supplies (−1)^l; the remaining caustic factor is
**−i forward**, **+i backward**. Every empty focus carries this phase. Open-mouth capture uses an explicit
`capture_phase`: `after_caustic` (default) includes it before capture, whereas
`before_caustic` admits the incoming front without that phase. The free-sphere
identity does not determine an open-mouth boundary condition. Thus two forward empty
foci give −1, not +1. For a pulse released at rest, the real field vanishes
at π while its conjugate quadrature is nonzero; at 2π the real field is
sign-flipped. Tests compare against both the integrated free conformal field
and `pair_birth.bare_timeline`. This phase convention is specific to the
conformal S² reduction; it is not an exact phase law for the massless spectrum
or a derivation for the full BAM S³ geometry.

Amplitudes are arbitrary complex numbers stored as `[real, imaginary]`. They
carry phase bookkeeping, not quanta or a computed full field profile. Closed
boundary reflection is separate: `end_condition='dirichlet'` gives −1 and
`'neumann'` gives +1. A closed mouth at a surface focus applies both the
caustic factor and its boundary factor. Open ports still transmit perfectly
by assumption.

## Topology is explicit; aperture interior is an extension

`handle='klein'` (default) or `'torus'` specifies topology independently of
`chart_side`. In sphere-oriented charts, the existing mesh's Klein map is −I;
the torus map is diag(1,−1). For source and destination chart changes C_s and
C_d, the API derives **M_chart = C_d M_geometry C_s**. Changing coordinates
therefore cannot silently change the handle's orientability. The old
`spatial_parity` parameter has been removed: its sign alone did not specify
topology.

With the displayed outer/inner chart choices (+1,−1), the Klein map becomes
(u,v)→(u,−v). Tests follow actual shared rim columns for **both topologies and
all four chart-side combinations**, and check mesh orientability independently.

**The interior aperture map is an assumed disc-to-disc portal extension.**
The glued two-dimensional mesh contains the rim circle; its removed disc has
no throat image. Only |local|=a is validated by that mesh test. Applying the
same matrix at |local|<a is an explicit transport assumption, not validation
of an interior throat field, metric, or volume.

## Time direction and the imposed loop

At mouth i, σ_i(t)=ε_i s_i(t)+b_i. Each throat leg obeys
**σ_out=σ_in+d_σ L**; the destination's actual clock history is inverted without
extrapolation. Two policies are now distinguished:

| Policy | Entry into throat | Exit onto surface |
| --- | --- | --- |
| `transport_orientation` (default) | d_σ=ε_in d_t | d_t,out=ε_out d_σ |
| `channel_forward` (explicit control) | Reset d_σ to +1 | Reset d_t to +1 |

The default transports the characteristic's time orientation through both
attachments. The other convention resets the arrow in each channel and is
additional boundary data. Playback follows the route parameter; exterior or
throat time may decrease. It never silently changes an advanced edge into a
forward-time flight.

With `channel_forward` and `capture_phase=after_caustic`, ε_B=−1 and anchor 16, the supplied data impose the
closed null route

    A:1 → B:4.141593 → A:12.423407 → B:15.565000 → A:1.

This is a **closed null route built into this boundary model**, not a dynamical
result or a derived GR causal geometry. Its two caustics give loop gain G=−1.
For the homogeneous scalar front condition a=G a, **only a=0 is consistent**;
the displayed unit seed returns as −1, with mismatch 2. A driven/broadband
boundary problem has not been solved. The replay reports that mismatch and
stops at the recurrence instead of treating the nonzero seed as consistent.

The paired `reversed_precaustic` control retains channel-forward resets but
captures before the caustic. Its gain is +1 and any homogeneous front amplitude
is consistent. This is a boundary-convention comparison, not a selected nonzero
field solution.

The paired `reversed_oriented` control uses the same clocks with the default
orientation transport. It contains past-directed surface legs and does not
close that four-leg loop; it reaches the finite history boundary.

Recurrence detection uses explicit absolute tolerances (default 1e−9 in unit
sphere/time units), matching channel, time direction and throat endpoint. It
reports time and position residuals, gain, and seed mismatch. Near-returns are
reported with their residuals even if the leg budget is exhausted. A return
within tolerance is never presented as an exact geometric theorem.

## Eight prescribed-motion controls

| Control | Surface captures | Throat exits | Antipodal reflections | Throat reflections |
| --- | ---: | ---: | ---: | ---: |
| Fixed antipodal mouths | 5 | 5 | 0 | 0 |
| Mouth displaced from first focus | 3 | 2 | 3 | 0 |
| Moving mouth reaches first focus | 1 | 1 | 5 | 0 |
| Moving mouth leaves first focus | 0 | 0 | 6 | 0 |
| Reversed clocks, channel arrows reset | 2 | 2 | 0 | 0 |
| Reversed clocks, pre-caustic capture and reset arrows | 2 | 2 | 0 | 0 |
| Reversed clocks, transported arrow | 5 | 4 | 0 | 0 |
| Closed Dirichlet throat outlet | 2 | 3 | 2 | 3 |

Counts cover finite supplied histories. They are route bookkeeping, not a
spectrum or discrete momentum exchanges. Moving mouths follow prescribed
inertial histories with proper clocks, chosen to exercise a hit and a miss.
The closed-outlet control starts with a front in the throat as separate initial
data. Unsupported continuations remain in JSON and are never called absorption.

## Scope and archive decision

This remains characteristic transport on supplied worldtubes, **not a
self-consistent field/recoil solver**. No old local-port force, Coulomb
acceleration, impulse estimate, or energy ledger is imported as a result of
the corrected support. Self-acceleration, radius shrinkage, collapse, particle
cycles and QED-like momentum exchange remain unestablished.

The original fourteen-case experiment and its negative results remain both
in the tree and pinned at `48d9dc3ede4edaa6220dd86f3d162aa58df2944c`. Retaining
the archives is deliberate: the user requested their publication and preservation
of the failed sustained-approach/collapse result. The review's optional removal
would reduce checkout size, but is not required for correctness. Its historical
manifest describes the original commit, not this correction.

## Reproduction and review validation

```bash
python -m experiments.closure_ledger.antipodal_throat_replay
python -m pytest -q tests/test_antipodal_throat_transport.py tests/test_history_throat.py
```

The replay generator resolves paths from the repository root. All manifest
keys are repository-relative and are tested by copying and verifying the files
under a different checkout root. The eight controls and hashes are in
[`runs/20261021_antipodal_transport`](../experiments/closure_ledger/runs/20261021_antipodal_transport).

The 32 transport/reproduction tests and 12 historical-action regression tests
pass (44 total). New checks cover conformal quadrature and return phase,
Dirichlet/Neumann outlets, chart-invariant topology, transported time direction,
loop gain, near-recurrence tolerances, and portable manifests.

The review of `03bd7ca` independently checked that the replay loaded and played
in headless Chromium without page errors or network fetches. That browser result
applies to that revision. The revised eight-case replay is checked with the real
canvas renderer over 638 sampled frames across all supported legs, finite-coordinate
checks, interaction controls and playback termination. It is self-contained.

Optional canvas QA (requires the Node package `@napi-rs/canvas`):

```bash
node experiments/closure_ledger/check_antipodal_throat_viewer.cjs
python -m experiments.closure_ledger.antipodal_throat_replay
```

The final command refreshes the manifest after the preview image is rendered.


## Relation to the newer #345 signed-exchange experiment

The follow-ups at commits `ad746d6` and `0380308` add moving-front energy and
angular momentum, a signed capture rule, tube-mapped exit recoil, and contact/
rebirth rules to #345. The pinned [source](https://github.com/davidmdrpi/geometrodynamics/blob/0380308e273042155c3c31aa755f48ed5f67d0a9/geometrodynamics/waves/antipodal_pair_birth.py)
and [report §10](https://github.com/davidmdrpi/geometrodynamics/blob/0380308e273042155c3c31aa755f48ed5f67d0a9/docs/pair_birth_timeline.md#10-the-klein-throat-as-a-reversal)
were inspected; their 22 trajectories were not rerun for this routing repair.

Their latest correction supersedes the earlier claim that Klein transport
simply reverses all angular momentum. Under their stated tube rule and
antipodal birth-frame convention, both handles turn J by π. Along the birth
circle, both reverse the orbital contribution; their distinction in the
planar runs is the separately assumed compression sign flip for Klein.
The reported signed exchange is attractive in the symmetric Klein controls
but provides only about 3–8% of the required stopping change; no pair binds.
Those are results of that specified exchange rule, not a derivation of a force
from orientability alone.

There is a material boundary convention difference: #345's `_focus` applies
−i only on the **no-hit** branch, so its signed capture uses the pre-caustic
phase. #347 now exposes that choice rather than silently treating its default
post-caustic capture as the same model. #347's complex scalar coefficient also
does not include #345's separate orientation-sensitive compression sign.
Its throat-coordinate map is not a momentum or normal-displacement transport law.

A static-flux-mediated exchange is the proposed next model. It still needs a
specified field/boundary coupling and recoil ledger; routing the old continuous
Coulomb force by name would not derive it. This PR repair does not add that
force law or import the reported binding fractions as its own predictions.
