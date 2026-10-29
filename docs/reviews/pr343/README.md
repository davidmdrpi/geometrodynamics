# Review assets for PR #343: wave contact in a narrow neck

All scripts import PR #343's own model, `scripts/neck_wave_contact.py`, with its default
parameters (b = 0.16, A = 0.18, ε = 0.025).

## First review (`df9a66e`)
Run these from a checkout of `df9a66e`. They use that version's `branches()` cap rule, which
the next revision removed.

    PYTHONPATH=scripts python docs/reviews/pr343/contact_checks.py
    PYTHONPATH=scripts python docs/reviews/pr343/contact_figure.py docs/reviews/pr343/contact.png
    PYTHONPATH=scripts python docs/reviews/pr343/diffuse_shell.py docs/reviews/pr343/diffuse.png

## Follow-up review (`cb4b5ec`, option C)
Run these from a checkout of `cb4b5ec`:

    PYTHONPATH=scripts python docs/reviews/pr343/field_checks.py
    PYTHONPATH=scripts python docs/reviews/pr343/crossing_figure.py docs/reviews/pr343/crossing.png
    node docs/reviews/pr343/browser_check.cjs docs/visualizations/finite_neck.html /tmp \
        "t=4.25;playing=false;draw()" \
        "document.getElementById('amp').value='.30';document.getElementById('amp').oninput();t=4.25;playing=false;draw()"

- **`field_checks.py`** checks the following:
  - recomputes both event pairs and their time symmetry about t = 4.25;
  - checks that ψ is even in ρ and in R;
  - locates the sampled Hessian maximum: it sits at a rim, and the value near the contact point is also reported;
  - counts the axis bridges at A from 0.18 to 0.30;
  - tests the smooth positive-part clamp R₊ = (R + √(R² + (ε/4)²))/2.
- **`crossing_figure.py`** draws `crossing.png`: the wall crossing at A = 0.30, and the clamped alternative.

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
- **`browser_check.cjs`** loads an HTML page in headless Chromium, using a globally installed Playwright. It reports console errors, horizontal overflow at 390 px and the status line. It also takes a screenshot after each optional JavaScript snippet. The first review used an earlier, page-specific version, which is kept at commit `6207f4d`.

The alternatives are drawings of the same prescribed R(x, t). They are not solved surface
dynamics, and under that wave all of them reopen at t = 4.971. This is a classical
geometric comparison; nothing here is a GR calculation or makes a quantization claim.
