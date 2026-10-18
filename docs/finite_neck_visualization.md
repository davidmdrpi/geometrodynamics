# Finite neck: a geometric visualization

[Open the interactive study](visualizations/finite_neck.html) (download and open
locally if viewing on GitHub). It is self-contained and makes no network requests.
Play/pause, scrub time, turn the packet off, rotate the cutaway, or directly
adjust the pre-surgery waist. Reduced-motion preference starts playback paused.

![Finite-neck sequence](images/finite_neck/finite_neck.gif)

## What this shows

A finite bridge between two fixed rims, a prescribed signed deformation packet,
a narrowing waist, and an explicitly imposed replacement by two closed caps.
The cutaway removes a front sector only for visibility. Both the 3D surface and
axial section use the same radius profile. The geometry is independent of the
packet: switching off the packet does not prevent the scripted narrowing.

Rim radius a=1, half-separation h=.6, initial waist b0=.745071089. The fixed-rim
cosh family is

    R_c(z,b) = b cosh[acosh(a/b) z/h].

Only its initial member is identified as the equilibrium catenoid. As b falls
from .5 to .2, a smooth interpolation changes the profile to

    R_d(z,b) = sqrt[b² + (a²-b²)(z/h)²].

Both profiles have radius a at z=+/-h and radius b at z=0. The latter approaches
a double cone as b tends to zero. This geometric limiting statement does not
establish a fluid pinch-off solution. The visualization never reaches b=0:
at t=6.8, b=.045, an explicit topology rule replaces the connected surface by
two closed caps with a positive axial gap. Each cap uses R/a=sqrt[(u²+.03u)/1.03],
where u runs from zero at the closed tip to one at the fixed rim. The square-root
behavior rounds the tip. Their axial lengths decrease from .98h to .02h,
approaching flat discs according to a prescribed schedule.
This replacement need not be a continuous material trajectory.

The incoming packet multiplies the displayed radius by 1+u(z,t), with |u|<=.08
and u=0 at both rims. It is a prescribed moving envelope, **not a solved wave
on a narrowing membrane**. Color encodes the sign of u; it is not an energy
scale. No azimuthal filtering, inverse-taper shoaling, unavoidable pinch-off,
released pulse, or conservation ledger is claimed.

This supplements the #342 sphere storyboard with an inspectable local neck.
It does not replace that renderer or modify the neck-collapse protocol.
The next physical calculation is to determine R(z,t) from a specified membrane
and surrounding-fluid dynamics, and then test whether an incoming wave drives
it to the admission boundary. The regularized topology rule remains separate.

For experimental context, see Goldstein, Pesci, Raufaste and Shemilt (2021),
[Geometry of catenoidal soap film collapse induced by boundary deformation](https://doi.org/10.1103/PhysRevE.104.035105).
Real collapse can develop a central thread, multiple pinch sites and satellite
bubbles; the present rounded-double-cone-to-caps sequence is deliberately simpler.

## Reproduce

```bash
python scripts/render_finite_neck.py
```

Requires NumPy and Pillow. The script verifies fixed rims, the prescribed nominal
waist, positivity before surgery, and closed tips after surgery, then writes
the GIF, a close-up, a sequence still and their SHA-256 manifest. The HTML uses
the same equations with user-adjustable display controls. No coupled simulation
is run and no existing experimental archive is changed.

Validation: the Python geometry checks pass. The HTML JavaScript was executed
with a mock canvas to check time/waist/wave controls and geometry; a local browser
was unavailable, so browser-specific layout was not visually tested. The rendered
GIF and PNG stills were visually inspected.
