# Pair birth to annihilation: antipodal returns, attraction and release

A classical experiment on #345, answering the acceptance test set in its review (end of
[mouth_kick_experiments.md](mouth_kick_experiments.md)).
- **The story:** two mouths are born touching at a point P and fly apart with the momentum
  left over from pair creation. Then:
  - the wave emitted at birth refocuses at the antipode and comes back to meet them;
  - the pair slows, turns round, and is pulled back together;
  - the rims touch, the handle is removed, and its energy leaves as a wave.
- **The mouths move independently.** Every force on them comes from one action, and nothing
  about the motion is scheduled.
- **Classical throughout.** Nothing is quantized; "exchange" here means the retarded field of
  the partner.

![Timeline](images/pair_birth/timeline.png)

**Answers in brief.**
1. **Without mouths,** a pulse released at P refocuses at −P at t = π. It comes back to P at
   t = 2π with its sign flipped, and repeats (panel a).
   - This is exact for the conformal operator.
   - The massless operator, which carries the mouths' field, refocuses to within a few per cent
     over two loops.
2. **The mouths attract because the flux is conserved.**
   - Two mouths sharing one handle carry equal and opposite flux (Wheeler's charge without
     charge).
   - The static energy of their field grows with separation, so they attract. The same pair
     treated as scalar sources would repel.
   - On the sphere the pull vanishes at antipodes, where the field lines are equally long both
     ways round.
3. **Born touching, apart, back: two outcomes.**
   - Without a field, the mouths always meet again on the far side, at −P, at t = (π − 2a)/v.
   - With the field, pairs slower than v ≈ 1/3 turn round and annihilate back at the birth point
     P.
   - Faster pairs pass the antipodal configuration and meet at −P, later than free motion.
   - At v = 1/3 the pair reaches 163°, stays beyond 143° for about 13 time units, and is
     pulled back to P.
4. **Waves decelerate.**
   - The sudden birth radiates a flash, which carries 13% of the birth kinetic energy.
   - The mouths pay for it, so the speed needed to cross over rises from 0.30–0.32 (no flash) to
     0.333–0.34.
   - The partner's flash reaches each mouth at t ≈ 0.8 and pulls it towards the partner. That is
     the first exchange.
5. **The flash comes back at the predicted times.**
   - After one or two antipodal foci, each returning ray pushes the mouth it crosses along the
     ray's direction of travel.
   - The returning ray from −P meets the outgoing mouth head-on and pushes it back towards P.
   - This held in 40 of 40 event windows across ten runs. Windows placed at the wrong
     (unshifted) times agree in 21 of 32, with half the impulse.
   - Later returns keep their timing but not a fixed sign.
6. **Near the threshold, timing decides where annihilation happens.**
   - Take the same birth at v = 1/3, with the flash re-timed after it has left the mouths.
   - Returns 0.5 or 1.0 later send the pair over to −P.
   - Returns earlier, or the flash absorbed, bring it back to P sooner.
   - The energy change of the intervention does not order these outcomes.
7. **Annihilation follows a stated rule.**
   - The glued surface is regular while the rims are disjoint. The dynamics brings them into
     contact at a computed time.
   - At contact the handle is removed. Its energy (rest, kinetic and Coulomb; 2.19 at v = 0.30)
     and its angular momentum leave as one deformation pulse.
   - The pulse refocuses at the antipode of the contact point and returns with its sign flipped,
     as in #344.
8. **The complete history closes the cycle, and the throat clocks follow it.**
   - The released pulse refocuses at the antipode after π, and the pair is re-born there. That
     is the archived v39 cycle, now with exact energy and angular momentum: five births and
     four re-births in 60 time units.
   - Each mouth carries its proper time since its own birth. The throat offset is the
     integrated difference, so it depends on the whole velocity history, never on position.
   - Symmetric pairs have no offset at all.
   - Even a lopsided birth (0.45 / 0.15) peaks at 0.215 mid-flight and is back below 0.0011 at
     contact. The offset resets at every birth.
   - A closed loop through the throat would need |offset| ≥ d + L_neck. The margin never falls
     below 1.155.
9. **With a neck channel, the offset acts back on the motion.**
   - The neck carries the field through the throat. Signals cross it in L + A one way and
     L − A the other, and the offset's reaction is a clock momentum μ that shifts the two rest
     energies by ±μ.
   - In the lopsided birth, the history clock brings the first contact 0.012 earlier than the
     same neck with its clock off, and the reversed clock 0.007 later. That is robust to the
     port stiffness and the neck mesh.
   - The reaction is small: μ ≤ 3 × 10⁻³. The effect compounds over later cycles to contact
     shifts of 0.05–0.4.
   - The neck itself, with or without its clock, changes the history more. Symmetric pairs feel
     no clock at all.

**What this is not.**
- **The field and mouths are effective.** The sphere-side field is 2+1 Maxwell on the round
  unit sphere, and the mouths are point particles with Gaussian footprints.
  - The neck is not part of the field domain, unlike #338/#345's glued mesh. The optional neck
    channel of §8 adds its interior as a 1-D wave channel attached by ports.
  - Its shape is rigid, and its static energy is a fixed rest energy M per mouth.
- **The removal is a rule, not dynamics.** The surgery at contact is a stated rule. #342's
  waist-barrier collapse is not coupled here.
- **The static part looks instantaneous only in this gauge.** The calculation uses Coulomb
  gauge. The gauge-invariant total field of each mouth is retarded.
- **None of this is QED.** It is a classical analogue, with no quantization step anywhere.

## Model
`geometrodynamics/waves/pair_birth.py`.

**Mouths.**
- There are two mouths, with independent positions x_A, x_B on the unit sphere and independent
  momenta.
- They are relativistic particles of rest energy M. That energy stands for the handle: neck
  surface energy plus the flux energy inside the neck. It is held fixed.
- Each mouth has a Gaussian footprint of width f, which smooths its charge and current as #340's
  footprint does.

**Charges.**
- A carries flux +q through the handle and B carries −q.
- On the sphere side, Gauss's law forces the two charges to be opposite whatever the gluing, so
  Klein and torus handles do not differ here.

**Field.** 2+1 Maxwell on the sphere, in Coulomb gauge.
- The longitudinal part has energy U(x_A, x_B) = ½∫|∇φ|².
  - It is fixed by the mouth positions.
  - It is a Legendre series in cos d, where d is the mouths' separation.
- The transverse part is A = n × ∇χ, with modes ω_l = √(l(l+1)) for l ≤ 32.
- Dual description: in 2+1 dimensions this is one scalar ψ, with E = n × ∇ψ and B = ψ_t.
  - ψ winds by ±q around the two rims, so the mouths are a vortex and an antivortex of one
    scalar wave field.
  - That is the 2-D form of Wheeler's charge without charge, and the dual of #345's static flux
    field.

**Action.**

  L = Σ_k [−M√(1 − v_k²) + q_k v_k·A(x_k)] − U(x_A, x_B) + ½ Σ_lm l(l+1)(|ċ_lm|² − l(l+1)|c_lm|²)

- **Conservation.** Total energy and total canonical angular momentum are conserved: the
  mouths' x × (p + qA), plus the field's.
- **Two-sided forces.** Every change of mouth momentum is matched by one in the field.
- **The field is split by source.** It is kept as one copy per mouth, each driven only by that
  mouth's current. Linearity then splits the force on A exactly into:
  - the partner's field: Coulomb plus transverse, the retarded exchange;
  - A's own transverse field: local reaction, plus A's own emission returning.

**Integrator.**
- It alternates exact free-field flow (mouths frozen, transverse electric kick −qΔA) with mouth
  motion in the frozen field.
- That motion is a Coulomb kick, a magnetic rotation and an exact geodesic drift. The mouths'
  currents feed the field on a trapezoid rule.
- The step is dt = 0.004.

**Birth.**
- The rims touch at P, so the separation is 2a with a = 0.3 (#338).
- The field starts as the static Coulomb field of the touching pair, with no transverse field.
- Each mouth is kicked away from P along one great circle with speed v.
- The sudden start of the currents radiates the **birth flash**.
  - By linearity the flash is exactly the free wave launched by −Q_qs.
  - Q_qs is the quasi-static (magnetostatic) field of the initial currents.
  - So the flash's force on the mouths along any trajectory can be computed after the fact.

**Regularity law and release.**
- The glued surface is regular while the two rim discs are disjoint, d > 2a.
- **The contact time is dynamical.** When the dynamics brings the separation down to 2a, the run
  stops on the step that reaches 2a exactly (secant refinement).
- **Then the surgery rule:**
  - the mouths, their neck and their longitudinal field are removed;
  - their energy (rest + kinetic + Coulomb) goes into one pulse of the conformal deformation field
    (−Δ + ¼, #339/#344) at the contact point;
  - the pulse is shaped by least squares to carry the pair's canonical angular momentum, and its
    amplitude fixes the energy.
- **What continues afterwards:** the transverse field continues freely.

**Parameters.**

| quantity | value |
|---|---|
| rest energy per mouth M | 1 |
| flux q | 0.75 |
| rim radius a | 0.3 |
| footprint f | 0.12 |
| mode cutoff | l ≤ 32 |
| release pulse width | 0.25 |
| horizon | t = 34 |
| static barrier U(π) − U(2a) | 0.109 |

**Controls.**
- **no field:** q = 0.
- **instantaneous Coulomb:** no transverse field, so no waves.
- **no flash (dressed birth):** the transverse field starts as Q_qs.
- **flash absorbed at t = 3:** the known free flash is removed when it is near its foci, far from
  the mouths, as if it escaped to infinity.
- **flash re-timed at t = 3:** the flash is removed and put back advanced or delayed by up to 1.
  The birth is identical, so only the returns move.

## 1. The timeline without mouths
`bare_timeline` evolves an axisymmetric pulse released at rest at P. Panel (a) of the figure above
shows its energy density along the birth circle.
- **Conformal operator** (ω = l + ½):
  - the pulse focuses at −P at t = π;
  - it is back at P at t = 2π with its sign flipped;
  - it is back at P unchanged at t = 4π.
- **Massless operator** (ω = √(l(l+1))), which the mouths' field uses:
  - the peak at −P reaches 97% of the conformal one at t = π, and 91% at t = 3π.
- **Free mouths** (dashed, v = 1/3) meet the returning front at t = (2π − a)/(1 + v) = 4.49, at
  103° from P.
- **The flash has two sources.** With mouths, it leaves from the two start points ±a rather than
  from P.
  - Each mouth's own flash focuses at the antipode of its own start.
  - Its partner's flash focuses on the mouth's own side of the circle.

## 2. Why opposite flux mouths attract
![Outcomes](images/pair_birth/outcomes.png)

Panel (a) shows the static interaction energy against separation.
- **Flux mouths attract.** For the flux mouths it rises monotonically, by 0.109, from contact
  (d = 2a) to antipodes.
- **The pull vanishes at antipodes.** The slope there is zero, so the antipodal configuration is
  force-free.
- **Scalar sources repel.** The same charges treated as sources of the conformal scalar
  (#339/#340-style exchange) give the opposite slope: a source–sink pair is pushed apart.
- **So the attraction comes from conservation.** Opposite-orientation mouths attract because the
  flux through the handle is conserved (Maxwell), not because a scalar is exchanged.

## 3. Birth, separation and return
Contact time and place for each birth speed v. The free column is (π − 2a)/v.

| v | free (q = 0) | full field | no flash | instantaneous Coulomb |
|---|---|---|---|---|
| 0.200 | 12.71 at −P | 3.85 at P | 4.83 at P | 4.07 at P |
| 0.250 | 10.17 at −P | 6.12 at P | 8.36 at P | 6.69 at P |
| 0.280 | 9.08 at −P | 8.36 at P | 13.21 at P | 9.64 at P |
| 0.300 | 8.47 at −P | 10.78 at P | 33.90 at P | 13.52 at P |
| 0.320 | 7.94 at −P | 15.30 at P | 14.71 at −P | 23.36 at −P |
| 0.333 | 7.62 at −P | **26.92 at P** | 12.43 at −P | 15.07 at −P |
| 0.340 | 7.48 at −P | 21.38 at −P | 11.66 at −P | 13.57 at −P |
| 0.360 | 7.06 at −P | 13.94 at −P | 9.99 at −P | 10.97 at −P |
| 0.400 | 6.35 at −P | 9.89 at −P | 8.02 at −P | 8.44 at −P |
| 0.450 | 5.65 at −P | 7.76 at −P | 6.60 at −P | 6.79 at −P |
| 0.500 | 5.08 at −P | 6.44 at −P | 5.67 at −P | 5.77 at −P |

- **Free mouths never come back.** On a sphere they meet again anyway, but always on the far
  side.
- **With the field there are two outcomes.**
  - **Bound** (back to P): the pair turns round before antipodes.
  - **Crossing over** (to −P): the pair passes antipodes and meets at −P, later than free motion.
- **Where the threshold sits:**
  - between 1/3 and 0.34 with the full field;
  - between 0.30 and 0.32 without the flash, and also without waves at all.
- **The figure.** Panel (b) shows separations, and panel (c) contact times and places.
- **The user's sequence happens near threshold.** It is born touching, moves into near-antipodal
  positions, and is attracted back.
  - **Example:** at v = 1/3 the pair reaches d = 2.85 (163°) at t = 13.7.
  - **Lingering:** it stays beyond d = 2.5 from t ≈ 7.1 to t ≈ 19.8.
  - **Return:** it annihilates back at P at t = 26.92.

## 4. The birth flash and its returns
![Flash returns](images/pair_birth/flash_returns.png)

**Waves decelerate.**
- **The flash's energy.** It is a free wave, so its energy equals that of the dressing it removes:
  0.0128 at v = 0.30 and 0.0158 at v = 1/3, about 13% of the birth kinetic energy at every speed.
- **Where the energy comes from.** The mouths pay for it.
  - In panels (a) and (b), the reaction of each mouth's own flash in the first 0.3 time units
    (blue dip) decelerates it.
  - This is why the full-field threshold is the highest.
- **The first exchange.** The partner's flash reaches mouth A from behind at t = 0.79
  (v = 0.30) and t = 0.83 (v = 1/3), before any focus.
  - It pulls A towards its partner, with impulse −0.014 and −0.016.

**First returns.**
- After one antipodal focus, a ray heading back to P meets the outgoing mouth head-on and pushes
  it towards P.
- After two foci, a ray heading away from P overtakes the mouth and pushes it outwards.
- At v = 0.30 the four first-return crossings of mouth A are:

| time | flash from | ray heading | foci | impulse on A (window ±0.35) |
|---|---|---|---|---|
| 5.09 | partner | to P | 1 | −0.0081 |
| 5.69 | own | to P | 1 | −0.0187 |
| 6.82 | own | away from P | 2 | +0.0132 |
| 7.39 | partner | away from P | 2 | +0.0109 |

- **The statistics.** Across the two natural runs and eight re-timed runs:
  - **at the flash's own timing:** 40 of 40 first-return windows point along the ray, with mean
    impulse 0.0116;
  - **at the unshifted (wrong) timing:** in the re-timed runs, 21 of 32 windows agree, with mean
    impulse 0.0061.

**Later returns.**
- They keep their timing: panel (b) shows four loops at v = 1/3.
- The sign of their push is not fixed. The transverse modes are not exactly periodic
  (ω = √(l(l+1))), and the returns overlap the mouths' own continuing radiation (grey).
- **Size near antipodes.** The static pull falls to 0.006–0.01 there, while the returning flash
  pushes by ±0.02–0.05. Near antipodes, the force is mostly the returning wave.

## 5. Same birth, re-timed flash
At t = 3 the flash is near its foci and at least 1.4 rad from both mouths. It is removed and,
except in the absorbed case, put back advanced or delayed. The birth and its energy cost are
identical. The energy change at t = 3 comes from cross terms with the rest of the field.

| v | flash after t = 3 | returns shifted by | energy change at t = 3 | contact |
|---|---|---|---|---|
| 0.300 | natural | 0 | 0 | 10.78 at P |
| 0.300 | absorbed | — | −0.0113 | 10.31 at P |
| 0.300 | re-timed | +1.0 (later) | +0.0004 | 11.14 at P |
| 0.300 | re-timed | +0.5 (later) | −0.0028 | 11.05 at P |
| 0.300 | re-timed | −0.5 (earlier) | +0.0048 | 10.40 at P |
| 0.300 | re-timed | −1.0 (earlier) | +0.0042 | 9.98 at P |
| 0.333 | natural | 0 | 0 | **26.92 at P** |
| 0.333 | absorbed | — | −0.0145 | 20.81 at P |
| 0.333 | re-timed | +1.0 (later) | +0.0020 | **27.22 at −P** |
| 0.333 | re-timed | +0.5 (later) | −0.0024 | **30.78 at −P** |
| 0.333 | re-timed | −0.5 (earlier) | +0.0045 | 21.09 at P |
| 0.333 | re-timed | −1.0 (earlier) | +0.0028 | 18.12 at P |

- **Below threshold** (v = 0.30), the returning flash only moves the contact time, by −0.8 to
  +0.4, with the same place.
- **Near threshold** (v = 1/3), it decides the place.
  - Returns delayed by 0.5 or 1.0 push the lingering pair over to −P.
  - Earlier returns, or no return, bring it back to P sooner.
  - With the natural timing the pair still returns to P, but 6.1 later than without the
    return.
- **The energy change does not explain this.** A change of −0.0024 still ends at −P, and
  +0.0045 ends at P.
- **Compact versus open.** On an open space the flash would escape: that is the absorbed run.
  Annihilation would then come at 20.81 rather than 26.92.

## 6. Ledger, contact and release
![Ledger](images/pair_birth/ledger.png)

**Ledger rows:**
- rest energy 2M (neck surface and flux energy, fixed);
- mouth kinetic energy;
- Coulomb (longitudinal) field energy;
- transverse field: the waves;
- the released deformation pulse.

**Clock support:** no offset is imposed. The throat offset that the mouths' own histories
generate is tracked in §7. #345 §3 treats imposed offsets separately.

| v | kinetic at birth | birth flash | kinetic at contact | transverse field at contact | released |
|---|---|---|---|---|---|
| 0.200 | 0.0412 | 0.0057 | 0.0296 | 0.0116 | 2.1360 |
| 0.250 | 0.0656 | 0.0089 | 0.0642 | 0.0014 | 2.1705 |
| 0.280 | 0.0833 | 0.0111 | 0.0615 | 0.0218 | 2.1679 |
| 0.300 | 0.0966 | 0.0128 | 0.0826 | 0.0139 | 2.1889 |
| 0.320 | 0.1110 | 0.0145 | 0.0917 | 0.0193 | 2.1980 |
| 0.333 | 0.1213 | 0.0158 | 0.1012 | 0.0202 | 2.2075 |
| 0.340 | 0.1267 | 0.0164 | 0.1171 | 0.0096 | 2.2234 |
| 0.360 | 0.1437 | 0.0184 | 0.1168 | 0.0269 | 2.2231 |
| 0.400 | 0.1822 | 0.0227 | 0.1249 | 0.0573 | 2.2312 |
| 0.450 | 0.2396 | 0.0287 | 0.1901 | 0.0495 | 2.2964 |
| 0.500 | 0.3094 | 0.0354 | 0.2343 | 0.0751 | 2.3406 |

- **The kinetic deficit is in the waves.** Contact happens at the birth separation, so the Coulomb
  energy is the same at birth and contact. The mouths' kinetic deficit at contact equals the
  transverse field energy then.
- **Exchange during the flight.** Panel (a) shows that during the flight the kinetic energy goes
  into the Coulomb field and comes back. The waves take a share at birth, and exchange it at each
  return.
- **Release.** At contact the rest, kinetic and Coulomb energy go into one deformation pulse.
  - It focuses at the antipode of the contact point at t_c + π and returns with its sign flipped
    at t_c + 2π (panels b and c of the timeline figure).
  - With independent mouths the pulse also carries the pair's angular momentum. The residual is
    below 10⁻⁹.

**Independent mouths.**
- In this run the speeds are unequal and one mouth is kicked off the birth circle (panel c).
- The mouths move apart on separate tracks, come back together and touch off the circle at
  t = 10.26.
- Energy drift is 1.3 × 10⁻⁷. The angular-momentum drift is 2.2 × 10⁻⁹, through the release.

## 7. Complete history: throat clocks and the closed cycle
![History](images/pair_birth/history.png)

The archived v39 cycle (in PR #1's `archive/geometrodynamics_v39.py`) runs as follows:
- separation from the birth point;
- exchange across the sphere;
- attraction back together;
- annihilation;
- the annihilation wave travelling to the antipode;
- re-creation there, and the cycle repeats.

`run_history` now does the same with the action-derived model.

**Re-birth.**
- At t_c + π the conformal release pulse has refocused at −c.
- The rule (`recreate`) removes it and places a touching pair at the antipodal images of the
  contact positions. It starts from the contact momenta, which already point apart there.
- Least squares in the tangent planes then correct those momenta, so the new pair carries
  exactly the pulse's energy and canonical angular momentum.
- For a symmetric pair the mirror needs no correction.
- The transverse field is left untouched.
- Like the release, this is a stated rule.

**Throat clocks.**
- Each mouth carries its proper time τ_k since its own birth, integrated with the midpoint
  1/γ of every step. The throat identifies equal proper times on the two mouths (MTY).
- The throat's external-time offset is therefore t_B − t_A with τ_B(t_B) = τ_A(t_A). To first
  order it is γ_B(τ_A − τ_B).
- So the offset is a functional of the two velocity histories since birth. It is zero at
  every birth, and clock memory does not survive the surgery.
- Chronology check: a closed null loop through the throat needs |offset| ≥ d + L_neck, with
  L_neck = 2S = 0.565 (#338).

| run | births | contacts | max throat offset in each cycle | offset at contact | min chronology margin |
|---|---|---|---|---|---|
| symmetric, v = 0.30 | 5 | P, −P, P, −P, P | 0 | 0 | 1.165 |
| symmetric, v = 1/3 | 3 | P, −P, P | 0 | 0 | 1.165 |
| unequal speeds + off-circle kick | 5 | off the circle, wandering | 0.045, 0.048, 0.048, 0.044, 0.041 | ≤ 0.0093 | 1.155 |
| lopsided, 0.45 and 0.15 | 5 | on the circle at 113°, 24°, −92°, 146°, 14° | 0.215, 0.170, 0.115, 0.094, 0.072 | ≤ 0.0011 | 1.164 |

**What the history shows.**
- **Every cycle closes.**
  - Contacts alternate P, −P, P, … for symmetric pairs, as in v39. With net angular momentum
    the contact point wanders.
  - Total energy holds to 5 × 10⁻⁷ over 60 time units and nine surgeries.
  - Angular momentum holds to 3 × 10⁻⁸ (10⁻¹⁷ for symmetric pairs).
- **Waves accumulate.**
  - At v = 0.30 the pair is re-born at speeds 0.300, 0.279, 0.262, 0.246 and 0.221.
  - The transverse field grows from 0.014 after the first flight to 0.041 after the fifth.
    Each birth inherits the previous cycles' waves (panel d).
- **No clock offset for symmetric pairs.** Both mouths age identically, so the offset is
  exactly zero throughout. A symmetric cycle has no MTY ingredient at all.
- **Asymmetric births give small, transient offsets.**
  - The offset grows while one mouth is faster. It cancels once the attraction has reversed
    the roles: the mouth born fast arrives slow, and the slow one arrives fast.
  - Even the lopsided birth peaks at 0.215 and is back within 0.0011 at contact. The peak
    falls cycle by cycle as energy goes into the waves.
- **No chronology horizon.** Separation plus neck length always exceeds the offset by at least
  1.155.
- **An imposed offset is not reachable here.** An offset of order 1, such as the −1.2 imposed
  at t = 0 in #346's port model, cannot arise from these histories. It would need a sustained
  velocity-squared difference of about 0.1 over T ≈ 24 without the attraction reversing it.

## 8. Neck channel: the offset acts back
![Neck channel](images/pair_birth/neck.png)

In §7 the offset was bookkeeping. Nothing travelled through the throat, so nothing could feel
it. `NeckChannel` makes the throat interior a wave channel, from the same action, with its
timing set by the history offset.

**The neck.**
- **Geometry.** #338's catenoid, fixed by the mouth angle alone:
  - the circumference is 2πf(s), with f = √(s² + b²), b = sin²a = 0.087 and s ∈ [−S, S];
  - S = sin a cos a, so the length is L = sin 2a = 0.565;
  - the circumference equals the rim's at both ends.
- **Field.** The transverse field's dual scalar ψ (B = ψ_t, |E_T| = |∇ψ|) continues into the
  tube. The neck keeps two parts of it:
  - **m = 0**, axisymmetric, with weight 2πf;
  - **m = 1**, the cos and sin of the rim angle, with weight πf and the centrifugal term 1/f².
    The narrow waist (radius 0.087) makes this a barrier, ∫ds/f = 3.78. On the glued mesh it
    reflected most of the sideways flux ([mouth_kick_experiments.md](mouth_kick_experiments.md)).
- **Ports.** Each neck end is tied by a spring of stiffness g to the matching rim moment of the
  sphere-side ψ around its mouth:
  - for m = 0, the rim average;
  - for m = 1, the cos and sin moments in the mouth's rim frame.

  Large g is continuity of ψ across the junction.
- **Rim frames and spin.** Each mouth carries a rim frame. The frame is transported with the
  mouth and turned by a spin s_k, with moment of inertia M sin²a.
  - So the ports' twist is booked, and angular momentum stays exact.
  - A spinning mouth feels the gyroscopic force s (x × v) of the sphere's curvature.
- **Handle.** The neck joins rim angle ψ at A to ψ + π at B (torus) or to π − ψ (Klein). The
  charges cannot tell the two apart (Model). The neck can, through its m = 1 part.
- **Mirror symmetry.** A pair moving in one plane through the centre makes B odd under
  reflection in that plane. So for coplanar flights (symmetric and lopsided) the m = 0 part
  vanishes identically, and the channel is the out-of-plane m = 1 row.

**Clocks in the neck.**
- **The slice.** The neck is ultrastatic in its own time T, and the throat joins equal proper
  times. So at external time t the neck slice runs from T = τ_A at A to T = τ_B at B. In that
  chart it has:
  - a lapse that interpolates the two mouths' proper-time rates;
  - a tilt dT/dσ = −A/L, with A = τ_A − τ_B the history offset of §7.
- **Crossing times.** Signals cross in L + A from A to B and in L − A from B to A. This is
  checked to 10⁻³.
- **The constraint.** Ȧ = 1/γ_A − 1/γ_B is imposed with a conjugate clock momentum μ. The
  momentum shifts the rest energies to M + μ and M − μ.
- **The Hamiltonian.** The lapse folds the neck and port energies into each mouth's rest
  energy:

      H = Σ_k √(m_k² + |p_k − q_k A_T|²) + U + H_T,   m_k = M ± μ + N_k + V_k + s_k²/(2M sin²a)

  - N_k is mouth k's lapse share of the neck energy and V_k its port energy.
  - H is autonomous, so energy and canonical angular momentum are conserved.
- **The clock force is the throat's energy flux.** dμ/dt = −∂H/∂A = (1/L)∫ N S dσ, where S is
  the energy flux along the neck.
  - A parcel of energy that crosses the whole neck from A to B adds its energy to μ.
  - μ is the offset's only reaction.

**Integrator.**
- The neck and ports evolve exactly, by matrix exponential, with frozen inputs: the mouths'
  rates, the offset, the frames and the rim moments.
- Simpson's rule integrates the reactions on the field, the mouths, their spins and μ.
- These flows are Strang-split with the field and mouth flows at dt = 0.004.
- A newborn neck starts at rest in the static continuation of the field it is born into: the
  minimum of neck plus port energy, which is empty for a field-free birth. It starts with μ = 0
  and no spin.

**Controls.** Each birth is run four ways:
- with no neck;
- with the neck and the clock off: the neck sees no offset, and μ ≡ 0;
- with the history clock;
- with the clock reversed: the neck sees −A.

The handle is a torus; the lopsided birth is also run with a Klein handle. The comparisons are:
- **history − off** isolates the offset's feedback;
- **off − no neck** is the neck's own effect;
- **reversed − off** should flip the first-order feedback.

**Calibration** (`neck_checks`: lopsided birth to its first contact, history clock against clock off).

| port g | nodes | dt | mean neck + port energy | μ at contact | contact shift, history − off | energy drift |
|---|---|---|---|---|---|---|
| 25 | 65 | 0.004 | 6.34e-04 | +2.88e-04 | −1.14e-02 | 6.6e-07 |
| 25 | 129 | 0.004 | 6.34e-04 | +2.67e-04 | −1.13e-02 | 6.6e-07 |
| **100** | **65** | **0.004** | **1.06e-03** | **+5.50e-04** | **−1.23e-02** | **1.4e-06** |
| 100 | 65 | 0.002 | 1.06e-03 | +5.47e-04 | −1.23e-02 | 4.6e-07 |
| 100 | 129 | 0.004 | 1.06e-03 | +4.02e-04 | −1.13e-02 | 1.4e-06 |
| 400 | 65 | 0.004 | 1.24e-03 | +1.02e-03 | −1.47e-02 | 4.3e-06 |
| 400 | 129 | 0.004 | 1.24e-03 | +6.04e-04 | −1.20e-02 | 4.3e-06 |

- **Port stiffness.** The neck's energy is set by g and converged in nodes. At g = 100 it holds
  85% of the energy at g = 400.
- **μ converges only at first order in h for stiff ports.** The port acts as a discrete Robin
  boundary.
  - At g = 100, μ at contact is 5.5 × 10⁻⁴ at 65 nodes and 4.0 × 10⁻⁴ at 129.
- **The feedback on the motion is robust.** The first contact comes earlier by 0.011–0.015 for
  every g and mesh, and by 0.011–0.012 at 129 nodes. Halving dt changes it by 3 × 10⁻⁵.
- **Production settings** (bold row): g = 100 with 65 nodes. The shift there is 9% larger
  than at 129 nodes.

**Results** (five births each, 60 time units, torus handle unless stated).

| birth | neck and clock | contacts at t | max \|A\| | max \|μ\| | neck + port, max | neck / waves (mean) | energy drift |
|---|---|---|---|---|---|---|---|
| lopsided (0.45 / 0.15) | no neck | 12.261, 25.317, 35.273, 44.553, 52.817 | 0.2121 | — | — | — | 5.0e-07 |
| | clock off | 12.235, 24.351, 36.665, 47.781, 57.928 | 0.2108 | 0 | 3.2e-03 | 0.050 | 3.3e-06 |
| | history | 12.223, 24.407, 36.598, 47.918, 57.896 | 0.2106 | 1.0e-03 | 3.3e-03 | 0.058 | 3.5e-06 |
| | reversed | 12.243, 24.306, 36.676, 47.651, 57.608 | 0.2110 | 1.1e-03 | 3.3e-03 | 0.054 | 3.6e-06 |
| | Klein, clock off | 12.234, 24.344, 36.197, 47.520, 57.400 | 0.2107 | 0 | 3.6e-03 | 0.057 | 3.8e-06 |
| | Klein, history | 12.223, 24.347, 36.270, 47.591, 57.661 | 0.2105 | 7.9e-04 | 3.5e-03 | 0.056 | 3.4e-06 |
| unequal + off-circle kick | no neck | 10.257, 22.048, 32.874, 43.048, 53.214 | 0.0481 | — | — | — | 4.0e-07 |
| | clock off | 10.224, 21.816, 32.235, 41.546, 50.002 | 0.0464 | 0 | 5.2e-03 | 0.075 | 4.5e-06 |
| | history | 10.224, 21.814, 32.186, 41.598, 50.378 | 0.0464 | 2.7e-03 | 4.3e-03 | 0.074 | 4.4e-06 |
| | reversed | 10.224, 21.818, 32.288, 41.489, 49.512 | 0.0463 | 2.8e-03 | 5.9e-03 | 0.074 | 4.4e-06 |
| symmetric, v = 0.30 | no neck | 10.781, 23.067, 34.219, 44.680, 54.084 | 0 | — | — | — | 4.1e-07 |
| | clock off, history or reversed | 10.755, 22.731, 34.466, 44.908, 53.156 | 0 | 0 | 6.1e-03 | 0.064 | 2.8e-06 |

Contact-time shifts by cycle:

| birth | contact | history − off | reversed − off | off − no neck |
|---|---|---|---|---|
| lopsided | 1 | **−0.0123** | **+0.0074** | −0.026 |
| | 2 | +0.056 | −0.045 | −0.97 |
| | 3 | −0.067 | +0.010 | +1.39 |
| | 4 | +0.136 | −0.130 | +3.23 |
| | 5 | −0.032 | −0.320 | +5.11 |
| unequal + off-circle kick | 1 | −0.00003 | −0.0002 | −0.033 |
| | 2 | −0.0014 | +0.0023 | −0.23 |
| | 3 | −0.049 | +0.053 | −0.64 |
| | 4 | +0.052 | −0.057 | −1.50 |
| | 5 | +0.376 | −0.490 | −3.21 |
| symmetric | 1–5 | < 10⁻¹¹ | < 10⁻¹¹ | −0.026, −0.34, +0.25, +0.23, −0.93 |

**What the neck shows.**
- **The offset acts back, at first order.**
  - With the history clock, the lopsided pair's first contact comes 0.0123 earlier than with
    the clock off. With the clock reversed it comes 0.0074 later.
  - So the part odd in the offset is −0.0099, and the even part is −0.0024.
  - The separation difference grows steadily through the first flight (panel b).
- **The reaction itself is small.**
  - μ stays at or below 1.1 × 10⁻³ (lopsided) and 2.8 × 10⁻³ (unequal + kick). The rest
    energies move by at most 0.3%.
  - The kicked birth has the smaller offset (0.046) but the larger μ. Its off-plane motion
    drives the m = 0 part, which crosses the neck without a barrier. The lopsided pair's
    sideways part has to tunnel through the waist.
  - Even so, its first contact moves by only 3 × 10⁻⁵. Its offset is a fifth of the lopsided
    one, and its feedback shows only in later cycles.
- **The feedback compounds along the history.**
  - By the second to fifth contacts the clock shifts contact times by 0.05–0.14 (lopsided) and
    up to 0.38 (kicked), mostly with the opposite sign for the reversed clock.
  - These are amplified differences of a sensitive multi-cycle history, not a per-cycle effect
    that grows: each changed contact moves the next birth.
- **The neck itself matters more than its clock.**
  - Clock off against no neck, every birth's first contact comes 0.026–0.033 earlier.
  - By the third to fifth contacts the histories differ at order one, with contacts at
    different places.
  - The neck holds 5–7.5% of the wave energy, and is emptied into the release at each contact.
- **Symmetric pairs feel no clock at all.** A ≡ 0, so the three clock modes agree to 10⁻¹¹.
  The neck still acts on them.
- **No runaway, no horizon.**
  - The feedback changes the lopsided pair's peak offset by 2 × 10⁻⁴ (0.2106 against 0.2108).
  - The chronology margin stays above 1.154.
- **The handle matters only in the details.**
  - Torus against Klein moves the first lopsided contact by 0.001, and the clock feedback is
    −0.0123 against −0.0109.
  - Later contacts differ by up to 0.5.

## 9. Checks

| check | result |
|---|---|
| energy, all 51 single-flight runs | maximum relative drift 2.4 × 10⁻⁷ |
| energy and angular momentum over repeated cycles | 5 × 10⁻⁷ and 3 × 10⁻⁸ through up to nine surgeries per run; each re-birth matches the pulse it replaces to rounding |
| free clocks (q = 0) | τ = t √(1 − v²) to 10⁻⁹; zero offset for equal speeds |
| angular momentum | symmetric runs below 10⁻¹⁸; independent mouths 2.2 × 10⁻⁹ |
| mirror symmetry of the forces on A and B | exact (difference 0) |
| step halved (dt = 0.002) | contact time changes by 4 × 10⁻⁶ (v = 0.30) and 2.5 × 10⁻⁴ (v = 1/3) |
| cutoff 40 instead of 32 | contact time changes by 4 × 10⁻⁸ and 5 × 10⁻⁷ |
| no field (q = 0) | contact at −P at 7.62478 = (π − 2a)/v |
| vector potential and magnetic force | match finite differences of the action |
| flash decomposition | for weak charge, default − dressed force equals the post hoc flash force (2%) |
| instantaneous Coulomb turning point | matches the energy budget |
| history offset | equals the integrated dilation difference; the external offset agrees with γ_B(τ_A − τ_B) |
| neck channel: energy and angular momentum | ≤ 4.5 × 10⁻⁶ and ≤ 3.2 × 10⁻⁸ over 60 time units in all 14 runs; every re-birth matches its pulse to 10⁻¹⁴; each neck-and-port flow conserves angular momentum exactly |
| neck crossing times | L + A from A to B and L − A from B to A, to 10⁻³ |
| neck clock and port forces | ∂H/∂A, the port forces on the mouths, the spin torque and the pull on the field all match finite differences of H (relative 10⁻⁵) |
| neck controls | a decoupled neck (g = 0) reproduces the bare pair to 10⁻¹³; the clock-off neck keeps μ ≡ 0; symmetric births keep A ≡ 0 and μ < 10⁻¹⁵ with either handle |
| coplanar flights | the m = 0 row stays below 10⁻¹⁵, as the mirror symmetry requires |
| tests | `tests/test_pair_birth.py`: 20 checks |

## Against the acceptance test

| requirement (review of 3428e48) | status here |
|---|---|
| independently moving, non-antipodal mouths | yes: two rotors with their own positions and momenta; an asymmetric run |
| a two-sided force, or a recoil from an action | yes: one Lagrangian; energy and canonical angular momentum conserved |
| event-aligned returns of the same field, with mistiming controls | first returns: 40/40 along the ray at the predicted times, against 21/32 at wrong times; re-timing changes the outcome. Later returns: timing only |
| ledger of wave, mouth, neck and clock-support energy | waves (transverse, deformation), mouths (kinetic), Coulomb; neck rest energy fixed; no imposed clock support, and the history offset is tracked (§7). With the neck channel: the neck field, ports, spins and the clock momentum μ, which carries the offset's reaction (§8) |
| collapse from a stated regularity law, not a schedule | rim contact sets the time dynamically; the surgery itself is a stated rule |

## Limits and next steps
- **Effective geometry.** The round sphere, point mouths and a rigid neck are used, not the
  glued mesh of #338/#345. The mouths move on the sphere only.
- **No collapse by amplitude yet.** The neck has no waist dynamics. The next step is to couple
  #342's metastable waist to the field at the mouths. Annihilation would then be a wave
  amplitude crossing the barrier, with contact as the backstop.
- **Approximate refocusing.** The transverse modes are massless, so they refocus only
  approximately. The released deformation pulse refocuses exactly.
- **Parameters not scanned.** The thresholds depend on q, M and the footprint. Only the birth
  speed and the flash timing were scanned.
- **The neck channel is a hybrid.**
  - The sphere stays whole, including the disc inside each rim. The neck is attached by port
    springs rather than glued, and only its m = 0 and m = 1 parts are kept.
  - The port stiffness is calibrated by convergence, not against the glued mesh. For stiff
    ports, μ converges only at first order in the neck spacing.
  - The neck's shape stays rigid.

## Reproduction
```bash
python experiments/closure_ledger/pair_birth_probe.py featured   # v = 0.30, 1/3 with field maps, ~1 min on 4 cores
python experiments/closure_ledger/pair_birth_probe.py scan       # 11 speeds x 3 field models, ~3 min
python experiments/closure_ledger/pair_birth_probe.py retimed    # flash absorbed or re-timed at t = 3, ~2 min
python experiments/closure_ledger/pair_birth_probe.py checks     # step, cutoff, free motion, independent mouths
python experiments/closure_ledger/pair_birth_probe.py history    # repeated cycles with throat clocks, ~1 min
python experiments/closure_ledger/pair_birth_probe.py neck_checks    # port stiffness, neck nodes, step, ~9 min
python experiments/closure_ledger/pair_birth_probe.py neck       # neck channel, four clock modes, ~13 min
python experiments/closure_ledger/pair_birth_report.py           # tables, figures, summary.json, sha256.json
python -m pytest -q tests/test_pair_birth.py
```
Results are in `experiments/closure_ledger/runs/20261019_pair_birth/`.
