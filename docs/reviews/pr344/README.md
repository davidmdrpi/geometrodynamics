# Review assets for PR #344: Klein throat release on the sphere

Run these from a checkout of #344 (`158113f`). They import the PR's `scripts/klein_release_model.py`
and the repository's `geometrodynamics.viz.klein_handles`:

    PYTHONPATH=.:scripts python docs/reviews/pr344/ledger_checks.py
    PYTHONPATH=.:scripts python docs/reviews/pr344/glued_drive.py 12000     # writes glued_drive.npz in the cwd
    PYTHONPATH=.:scripts python docs/reviews/pr344/glued_figure.py docs/reviews/pr344/glued.png

- **`ledger_checks.py`** independently checks the sphere-side numbers:
  - the release work, from an ODE integration of the forced modes: ∫F·p dt against the change in field energy;
  - the emitted cap energy at midflight and at the antipode, with a separate quadrature;
  - the 2π sign return;
  - the contact position, and the share of the release pulse emitted before a unit-speed signal from it could reach each mouth;
  - the #338 and #343 neck geometries.
- **`glued_drive.py`** evolves the PR's incident packet with −Δ + ¼ on two surfaces:
  - the bare sphere, as validation against the modal solution;
  - #338's glued surface, with one Klein handle at ±x, using the repository's intrinsic cotangent operators and lumped mass.

  It records the field at the neck waist, the energy inside the neck, and the energy on the sphere within 1 rad of mouth B. The argument is the number of sphere vertices; 8,000 and 16,000 agree to 0.3%.
- **`glued_figure.py`** draws `glued.png` from `glued_drive.npz`.

The glued run uses only the classical scalar wave equation on a fixed surface. It
tests what the PR's incident data does when the throat is part of the space. It is not
a model of neck dynamics, and makes no GR or quantization claim.
