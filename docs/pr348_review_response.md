# Response to PR #348 review at 11bf1d3

[Review](https://github.com/davidmdrpi/geometrodynamics/pull/348#issuecomment-6334818665)

The two P1 findings change the interpretation, not just the presentation. The
centred neutral single-mode model does not establish orbital attraction.

| Item | Resolution |
| --- | --- |
| P1-1: spectral source degeneracy | Replace the main source with `(1,0,0,0)`. Retain the antisymmetric source as an explicit analytic degeneracy control; add a biased pair. Report impulse/source energy and exposure alongside impulse/exposure. Add complete-throat-phase controls showing the old handle fractions were a sampling artifact. |
| P1-1: physical length and spectrum | Add `L=sin(2a)` controls weighted by the actual `PairBirth.flash` modal energies, with source copies summed before squaring. Preserve its Maxwell frequencies `sqrt(l(l+1))`. The mapping from that field to incident port amplitudes is unknown, so one-sided, biased and antisymmetric source choices remain labelled projections, not a reconstructed birth trajectory. |
| P1-2: zero tangential force | State the centred m=0 tangential resultant is exactly zero by azimuthal symmetry. Test the angular integral and its change under an m=±1 perturbation. Reclassify the finite-mass example as neck-axis motion, not orbital braking. No arbitrary dipole coupling is added. |
| P1-2: collapse interpretation | Normal/neck stress is the supported interpretation. A radius force still requires the shape-work projection of traction; negative axial stress alone does not prove shrinkage or cross a waist barrier. |
| P2-1: charged mouths | Restrict the sign-insensitivity statement to an uncharged passive junction. Implement the electric static–wave stress cross term and verify its integral gives `q E` in a local tangent patch. The global charged-field/throat coupling and impulse normalization remain unresolved. |
| P2-2: dimension | Compare both `z=b/a` (2D rim) and `z=(b/a)^2` (3D area). Balanced optimal ratios are 0.1716 and 0.4142 respectively; regenerate geometry controls and figure. |
| P2-3: caustic bins | Explicitly identify `closed_paths` as a positive-frequency carrier-envelope model. It does not propagate a real broadband birth flash or its Hilbert-transform reshaping. The weighted frequency control is separate, with the exact-focus phase assumption labelled. |
| P3-4: grid identity | Remove doubled-bin rescaling from validation gates and summary. It was exact by construction, not convergence evidence. |
| P3-5: surface seed | Describe the long-duration normalized stress as oscillating near zero, at about 10^-4 in the recorded later samples; retain all durations and avoid claiming one negative case as a mechanism. |

The main report and PR description supersede the initial candidate-attraction
interpretation. The original commit and its data remain available at `11bf1d3`.
No prior negative trajectory result is removed or rerun.
