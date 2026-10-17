# PR #339 review follow-up

The [independent review at 44cd287](https://github.com/davidmdrpi/geometrodynamics/pull/339#issuecomment-6206173894)
reports eight passing tests, exact reproduction of all 240 numeric summary
entries, and byte-identical CSVs and spectral render states across different
Python/NumPy/SciPy versions. Its regenerated figures have the same content
with rasterization differences. This strengthens confidence in the numerical
result without changing the scalar comparator's physical scope.

We replicated the review's additional experiments in a separate
[supplementary archive](../experiments/closure_ledger/runs/20261016_antipodal_kick_review/).
These are **post-review replications of known measurements**, not additional
pre-registered successes. The original protocol and numerical archive remain
unchanged.

## The field supplies the momentum

The directional field initially carries J_z=0.0560747662. The attachment
receives 0.0074792496 by the first-transit horizon: **13.34%** of that amount.
The dephased case receives 16.07%; the half-amplitude case receives 12.61%
of its smaller initial field momentum. "Absorbed" means transferred by that
horizon, not irreversibly dissipated. No direction is selected from a
directionless state. Symmetric nulls and mirrored recoil are symmetry
requirements; the nontrivial measurements are how much transfers and when.

## Same packet, attachment off the focus

Initial a,p are identical within each family. Only attachment position changes;
q and all attachment momenta start at zero. Both off-focus positions are π/2
from source and target. The original horizon and parameters apply.

| Packet | Attachment | Final norm of ΔL | Activity in W | Maximum sampled shape energy |
|---|---|---:|---:|---:|
| Symmetric | Focus | 1.59e-32 | Numerical zero | 0.01638 |
| Symmetric | Pole (0,0,1) | 5.43631e-4 | 0.1964% | 3.90033e-4 |
| Symmetric | Equator (0,1,0) | 5.43631e-4 | 0.1964% | 3.90033e-4 |
| Directional | Focus | 7.47925e-3 | 92.8800% | 0.01330 |
| Directional | Pole (0,0,1) | 4.39812e-4 | 0.1898% | 3.13404e-4 |
| Directional | Equator (0,1,0) | 2.04000e-5 | 0.5896% | 3.02242e-5 |

Directional focus/off-focus recoil ratios are **17.01** and **366.63**;
shape-energy ratios are about 42.5 and 440.2. This is stronger evidence for
focus-selective transfer than the dephased comparison because the launch is
unchanged. These are two placements, not an angular scan or a universal lower
bound. Local fluence and incidence differ; directional cancellation may also
contribute to the large equatorial ratio. The symmetric focus still has zero
recoil, so focusing does not always maximize momentum transfer.

Symmetric off-focus kicks initially point away from the source. This does
**not** establish an antipodal attractor: a local restoring response near the
antipode and a long-time budget remain missing. Damping and repeated driving
would change the conservative model.

## Self-emission and return-time braking

Start the field at zero, q=0.2, p_q=0, and either L=0 or L_z=0.03 (angular
speed 0.1). The unexcited moving control has q=0 and the same L_z. Evolve to
2π+1.2, integrating both π and 2π window endpoints exactly.

| Case | ΔL_z by t=2 | ΔL_z in π±0.6 | ΔL_z in 2π±0.6 | Final ΔL_z |
|---|---:|---:|---:|---:|
| Rest emitter | 0 | 0 | 0 | 0 |
| Moving emitter | −7.11331e-4 | −1.69993e-4 | −1.95279e-4 | −1.76531e-3 |
| Moving, unexcited | 0 | 0 | 0 | 0 |

The rest emitter's maximum |L| is below 1.2e-32. The moving emitter loses
5.884% of its initial angular momentum; the field gains it. This is braking
for this emitter and initial state, not a theorem that conservative systems
cannot recoil from directional emission.

**Initial interaction energy matters.** The moving emitter starts with
E_shape=0.025, E_interaction=0.020, E_kinetic=0.0015, E_field=0, E_total=0.0465.
At t=2 its field energy is 0.0245515, shape energy 0.0179644, interaction
energy 0.00255368, and kinetic energy 0.00142971. The shape energy's 28.14%
decrease is **not the amount radiated**: field energy also comes from the
coupling energy. Shape energy oscillates and temporarily exceeds its initial
value; it does not empty monotonically.

The review's return-window value −1.99e-4 is recovered at sampled times 5.68
and 6.88: −1.9931e-4. Exact endpoints 2π±0.6 give −1.9528e-4. This roughly
2% endpoint effect exceeds the time-step error, motivating exact event times.

Moving-emitter relative final-impulse changes are 8.26e-9 for cutoff 32→40
and 7.15e-5 for dt→dt/2. The return-window impulse changes by 5.82e-4
relatively under time refinement. Maximum relative energy drift is 2.45e-5
(6.12e-6 at dt/2); normalized angular-momentum drift is below 2.4e-14.

## What still needs isolation

Nonzero π-window braking is continuing local interaction: the emitter is not
at its own antipode. Around 2π, a free disturbance can return to its emission
region (with a sign reversal for this S² operator). Late braking is consistent
with that return, but comparing windows does not isolate its causal share.
Emission continues, the emitter moves, and the field samples a history of
source positions with a globally tailed footprint. "Nothing happens at π"
and "the memory channel is exactly one 2π delay" are too strong literally.
Before calling the entire late impulse an **echo-only kick**, decompose
retarded torque by emission time or design an explicitly budgeted suppression
control. Removing/absorbing waves changes the conservation accounting.

The next unknown is a broader same-launch placement/angle response, followed
by a source–receiver calculation: one attachment emits, another sits near the
antipodal arrival region, and both react to the field. The receiver's π event
is distinct from the emitter's 2π self-return. Isotropic emission at rest and
an exactly centered symmetric receiver must remain a recoil null; motion or
a directional internal mode can break that symmetry.

Freeze ranges and failure gates before those new scans. Driven emission,
dipolar modes, and damped repeated-pulse attractors are further extensions,
not established consequences. Include drive work/torque, a dipole's internal
angular momentum, and damping sinks in their budgets. Uniform scalar driving
cannot select motion from exact symmetry without a perturbation/instability
test. No driven-emitter or attractor scan was run here.

## Changes and reproduction

The movie now includes a separate **magnified 0–3 degree surface-track inset**
using actual X(t). Sphere positions remain unmagnified; simulation states
were not changed to make motion visible.

The off-focus extension exposed a diagnostic assumption: displacement had
been measured from (1,0,0). It now uses the actual initial X. An analytic
geodesic test checks that and exact diagnostic endpoints. Nine targeted tests
pass; the equations and original protocol are unchanged.

```bash
PYTHONPATH=. python experiments/closure_ledger/antipodal_kick_review_probe.py
python -m pytest -q tests/test_antipodal_kick.py
PYTHONPATH=. python scripts/render_antipodal_kick.py
```

The supplementary archive contains nine CSVs, a summary with full energy
budgets, and SHA256 hashes, separate from the frozen first-transit archive.
