# Isolating the exchange mechanism

This follows merged #340. The calculation isolates the leading mixed rotational
impulse of two wave packets in the existing effective Hamiltonian. It compares
an independent fixed-axis response calculation against the weak-amplitude limit
of two independently moving antipodal mouth pairs.

[Numerical results, frozen gates and figure](mouth_exchange_numerics.md).

## What is computed

Write the smoothed scalar at mouth hs as U_hs=<w_hs,a>_real, its shape coordinate
as q_hs, and e_hs=q_hs-U_hs. The interaction is kappa sum(e_hs^2)/2.
At fixed mouth positions the field/shape equations are linear:

\[
\dot a=p,\quad \dot p=-\omega^2a+\kappa\sum_{hs}e_{hs}w_{hs},\qquad
\dot q_{hs}=\pi_{hs}/\mu,\quad
\dot\pi_{hs}=-\mu\Omega^2q_{hs}-\kappa e_{hs}.
\]

Propagate the x- and z-packet initial data separately through these **coupled**
equations, with the same normalization as their sum. Holding axes fixed removes
the moving-geometry correction, not the reciprocal scalar response. Their
explicit mixed torque is

\[
\boldsymbol\tau^{xz}_{hs}=\kappa\left[
 e^x_{hs}\operatorname{Re}\langle w_{hs},\mathbf D a^z\rangle+
 e^z_{hs}\operatorname{Re}\langle w_{hs},\mathbf D a^x\rangle\right].
\]

The new adaptive ODE solver integrates this expression directly. It does not
subtract two finite-amplitude trajectories or call the moving-mouth time
stepper. All four shape oscillators and reciprocal forcing of the scalar are
present. The mixed field angular momentum is

\[
\mathbf J^{xz}_{f}=-\operatorname{Re}\left[
 \langle p^x,\mathbf D a^z\rangle+\langle p^z,\mathbf D a^x\rangle\right].
\]

Its change plus the integrated mouth torques must vanish. Fixed-axis torques
are response coefficients; finite-amplitude trajectories still use freely
moving independent pair axes, with their own kinetic energy and full ledger.

For fixed finite time, smooth dependence gives
L_h^mixed(lambda)=lambda^2 C_h+lambda^4 D_h+O(lambda^6).
The proper rotation (x,y,z)->(z,-y,x) exchanges the identical pair locations
and packet profiles and reverses y angular momentum. The mixed quadratic term
is invariant under interchange of packet weights, even though those weights
are unequal. Therefore C_1y=-C_2y for this preparation. This symmetry predicts
a specific coefficient relation; the computed nonzero coefficient and its
agreement with freely moving trajectories are the substantive numerical tests.
There is no analogous assertion that all components of the full recoil vectors
are equal and opposite. The finite-amplitude field can retain mixed y momentum. In the directional
case, mixed x and z field momentum already changes at quadratic order; the
quartic result is specific to mixed y, not the whole field vector.

## Reproduction and scope

Run from the repository root with NumPy, SciPy and Matplotlib installed:

```bash
python experiments/closure_ledger/mouth_exchange_probe.py
python experiments/closure_ledger/mouth_exchange_report.py
python -m pytest -q tests/test_mouth_exchange.py tests/test_moving_mouth_pairs.py tests/test_antipodal_kick.py
```

The protocol was committed locally before runs at `0206381` (published protocol commit `7426483`). PR340's reviewer
had already reported the coarser scaling trend. This is a replication and a
new independent reference calculation, not a claim of blind preregistration.
The archive preserves six exact event states per run, per-mouth and pair
impulses, field angular momentum, every energy component, and every-step
maximum conservation errors. It is an event archive, not a full field replay.

A first crossing of .65 radians stops integration and retains the last valid
state and bracket; the horizon is always bounded by pi+1.2. No post-collision
results are admitted. Long-time behavior and displaced-pair attraction are
outside this first-transit test.

This is classical angular-momentum exchange under an assumed round-S2
scalar/shape coupling. Neither the drawn necks nor a changing handled metric
carry the computed scalar. No through-handle transport, GR completion,
quantized action, or localized measurement law is established here.

Each moving object is an antipodally constrained mouth pair. The two pair axes
are independent; the four individual mouth positions are not unconstrained.
