# Deformation contact in a narrow neck

[Interactive study](visualizations/finite_neck.html) — download and open locally
when viewing on GitHub. No network or external packages are needed by the HTML.

![Wave-contact sequence](images/finite_neck/finite_neck.gif)

This revises the visualization from #342 in one essential way: **the resting
neck no longer shrinks independently of the wave**. A signed inward crest closes
the gap between opposing finite-thickness surfaces. Reconnection is triggered
by contact with that band, rather than by a preset pinch time. The wave-off and
weak-wave controls leave the same resting neck intact throughout the time window.

## Geometry and event

Rim radius is 1, half-separation h=.6, and the fixed resting waist is b=.16:

    R0(x) = b cosh[acosh(1/b) x/h].

This narrow resting shape is assumed; it is not asserted to be a stable
catenoid equilibrium. The prescribed signed displacement is

    u(x,t) = -A exp[-((x-c(t))/.18)^2] cos[14(x-c(t))] [1-(x/h)^2]^2,
    c(t) = -.85 + .20 t,   R(x,t)=R0(x)+u(x,t),   A=.18.

The rims remain fixed. Pink denotes inward deformation, teal outward deformation.
A radial display half-thickness epsilon=.025 defines opposing surface bands.
Their facing gap is 2[min_x R(x,t)-epsilon]. The first zero of this gap triggers
the contact rule, located with a time bracket and bisection on an 801-point
axial grid. This is **contact of prescribed geometry**, not a dynamical instability
calculation. `contact.json` records the parameters and event.

The default contact occurs off the center of the waist as the crest arrives.
The pre-contact connected radius remains positive; no zero-radius throat is
rendered. The contact ring is highlighted, and the axial close-up represents
the local contact region with a finite-width blend.

## Rounded reconnection and rim motion

The renderer blends from the connected surface to two rounded closed caps over
.5 display-time units. The close-up zoom follows the separating caps. The two cap tips start with a positive axial separation;
they are never drawn intersecting each other at an axis point. Near each cap tip,
the radius scales as the square root of distance from the tip, giving a rounded
surface of revolution rather than a conical endpoint. Rounded shoulders and
curved gold motion guides represent rolling/retracting rims. They are illustrative
kinematics, not material trajectories, vortices, or a computed velocity field.
No detached torus or energy reservoir is silently introduced.

The visible topology transition is an **explicit regularization**, not a claim
that a smooth embedded surface can change topology while remaining everywhere
regular. Finite film thickness, coalescence, reconnection and post-contact flow
require a material/continuum model. The local blend is not a phase-field PDE.
No surface mass, energy or momentum conservation is established by this visual.
It also does not claim azimuthal filtering, shoaling, or inevitable pinch-off.

## Controls and reproduction

The HTML has time scrubbing, wave amplitude and resting-neck width controls,
wave-off/weak-wave/contact presets, and a jump-to-contact button. Changing the
wave or neck recomputes the contact event; no contact means no reconnection.
Reduced-motion preference starts playback paused.

```bash
python scripts/render_finite_neck.py
python -m pytest -q tests/test_neck_wave_contact.py
```

Python rendering requires NumPy and Pillow. Model equations are in
`scripts/neck_wave_contact.py`. Checks cover no-wave/weak-wave nulls, first contact,
fixed rims, positive pre-contact radii, and separated closed post-contact tips.
The static sequence and animation are visually inspected; the HTML's matching
math and controls are checked with JavaScript. A local browser was not available
for a browser-specific layout check.

The #342 protocol and sphere storyboard are unchanged. This study replaces only
the earlier local visualization and its independent shrinking-waist schedule.
The next physical step would solve the wave and surface response together, then
check whether they actually reach the specified finite contact layer.
