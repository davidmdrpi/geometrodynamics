# Review assets for PR #343: wave contact in a narrow neck

All scripts import PR #343's own model, `scripts/neck_wave_contact.py`, with its default
parameters (b = 0.16, A = 0.18, ε = 0.025). Run them from a checkout of #343:

    PYTHONPATH=scripts python docs/reviews/pr343/contact_checks.py
    PYTHONPATH=scripts python docs/reviews/pr343/contact_figure.py docs/reviews/pr343/contact.png
    PYTHONPATH=scripts python docs/reviews/pr343/diffuse_shell.py docs/reviews/pr343/diffuse.png
    node docs/reviews/pr343/browser_check.cjs docs/visualizations/finite_neck.html /tmp

## The scripts

- **`contact_checks.py`** prints every number in the review:
  - the contact time from continuous minimisation;
  - radius and curvature of the outer, mid and facing surfaces at contact;
  - the approach rate of the facing surface;
  - the cap separation and radius jump at reconnection, and the cap-tip curvature;
  - the reopening time under the prescribed wave;
  - the amplitude and waist thresholds.
- **`contact_figure.py`** draws `contact.png`. It shows:
  - the ε-band through contact;
  - the PR's caps against the band;
  - curvature against time to contact;
  - a level-set (Morse) alternative;
  - the contact window.
- **`diffuse_shell.py`** draws `diffuse.png`. It shows the same prescribed walls as a smooth density of width ε, which needs no reconnection rule.
- **`browser_check.cjs`** loads the HTML in headless Chromium with a globally installed Playwright. It records:
  - console errors;
  - the contact readout;
  - the wave-off control;
  - horizontal overflow at 390 px;
  - screenshots.

The alternatives are drawings of the same prescribed R(x, t). They are not solved surface
dynamics, and under that wave all of them reopen at t = 4.971. This is a classical
geometric comparison; nothing here is a GR calculation or makes a quantization claim.
