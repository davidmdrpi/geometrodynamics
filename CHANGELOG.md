# Changelog

## Unreleased

- Two classical throats and the tensor product (note, not pre-registered;
  `docs/two_throat_nogo.md`; corrected after review). On S³_A × S³_B with
  intensity-only setting loads:
  - fibre coupling 2λ_c A_A A_B is **not** inert. A counter-rotating
    non-horizontal pair beats every horizontal one at λ_s = 2, λ_c = 1.99
    (S − loads 9.619 against π²). The claimed factorization is withdrawn, and
    the coupled Bell behaviour is open;
  - with a fixed preparation, port-symmetric loads break complementarity, with
    σ ≈ −2κu(1−u)(2u−1); u = ½ is stable only for κ < 1;
  - with λ_c = 0, shared-holonomy conditions are trivial, lopsided or
    unstable, because c(θ) is maximal at π. With λ_c ≠ 0 the costs do not
    add;
  - Bloch-sphere coupling reaches |S| = 3.34 with flat binary marginals, but
    Bob's setting changes Alice's local intensity statistics (0.667 against
    0.499).

- Alice-handle exchange, residual at the spectral nulls (freeze `f03b034`):
  an in-sample truncation refuted the radial-doublet candidate. The doublet
  alone reproduces the residual, which splits into an intrinsic term,
  second order in the drain rate and falling as (ω̄τ)⁻³, and a strip-profile
  mismatch term K·Δ/ω̄, with the mismatch measured at 3.9 Δ/ω̄. Verdict
  `CONFIRMED` on held-out (0.45, 2) and (0.55, 2): the four-term predictor
  gives all 16 null values to within 0.85–1.13 (frozen ×2). The doublet
  alone matches the full field to 1.02%. The mismatch and intrinsic terms
  match within 0.84–1.10 and 0.99–1.06. The depletion boundary is
  independent of the aperture to 0.0009, and the frozen validity rule (peak
  drain rate at most 0.5ω̄) held at every point. See `docs/alice_nulls.md`.

- Alice-handle exchange, V-curve and aperture scaling (freeze `291ab5f`):
  the defect floor of #370 is the detector window's spectral weight at 2ω̄,
  the counter-rotating coupling. A two-part predictor, found in-sample,
  predicts held-out points on (0.5, 2) and (0.35, 2) to within 0.79–1.76.
  The V minima fall at the predicted grid points, and the smaller aperture
  lowers the best defect 15× and 12×, to 1.1e-5 at η = 0.5 and 5.1e-5 at
  η = 0.9. The envelope ranking holds: sin⁴ beats sin², and the C^∞ bump is
  worse. Verdict `COUNTEREXAMPLE` through G4: the sin² spectral nulls fall
  where predicted, but a residual of up to 2.8e-4 makes them shallower
  (depth 0.50 against the frozen 0.4). That residual is unexplained. Two
  η = 0.9 points at the shortest window cannot reach that depletion. See
  `docs/alice_vcurve.md`.

- Alice-handle exchange under detector depletion (freeze `0fb3b6e`): in the
  circulating basis the lowest Klein doublet swaps its ±m channels over
  `π/Δ`, so a port reading the +m intensity at two visits sees u and 1 − u.
  A counter-rotating dashpot on the band, which absorbs only the channel it
  reads, keeps `|u + v − 1|` between 5e-5 and 2e-2 at depletions of up to
  90%. A static control breaks it by 0.500. Verdict `COUNTEREXAMPLE`
  through G2: the frozen claim that the defect vanishes as Δτ → 0 fails at
  the frozen window, where a floor of 2.8e-4 remains at (0.3, 2), η = 0.5.
  Post hoc, the floor falls about 400-fold as the window lengthens from 4
  to 32, pointing to non-adiabatic switching. That is a candidate for a new
  pre-registration, not a result. See `docs/alice_exchange.md`.

- Exact-neck Klein splitting law (freeze `834ac70`): Herring's law with the
  neck equation solved exactly, `E = 4 sin^{2m}(a) |D(0)D′(0)| / (ω N_m)`,
  where `D` is the neck solution vanishing at the far rim; nothing fitted.
  All 17 measured cells were in-sample. Verdict `CONFIRMED` on 14 held-out
  cells (`a = 0.25–0.55`, `m = 1–5`), after reproducing the 17 to 5.1e-8:
  the splitting, spanning 5.5e-2 to 5.8e-8, is predicted to at most 0.73%,
  within the frozen `0.2 a⁴ + 1e-3` (worst use 0.38). It removes at least
  93% of `H`'s discrepancy, against the frozen 80%, and N = 2000/4000 agree
  to 3e-4. At the diagnostic `a = 0.65–0.7`, the worst error is 1.7%. The
  residual of about `0.11 a⁴` is measured, not derived. On review, #360's
  continuum solver re-measured every cell: it agrees with the finite
  element to 1e-4, and the gates also pass on its values. See
  `docs/klein_exact_neck.md`.

- Klein splitting law (freeze `56055ab`): a closed-form tunnelling law for
  the lowest Klein doublet on #355's handle,
  `H = 4 κ_rim sin^{2m}(a) e^{−W_neck} / (ω N_m)`, derived with Herring's
  formula in the conformal coordinate and nothing fitted. Verdict
  `CONFIRMED` on six held-out cells (`a = 0.2, 0.4, 0.5`, `m = 1–4`): the
  measured splitting, spanning 1.5e-2 to 4e-7, exceeds `H` by
  `m(r − 1) = 0.02–0.15`, within the frozen `0.25/m`, with the error ordered
  in `a` as predicted and N = 2000/4000 agreeing to 7e-4. `a = 0.6` is
  diagnostic and exceeds the bound for `m ≤ 3`. Explains, without
  changing, #355's P4 counterexample. See `docs/klein_splitting_law.md`.

- Alice handle ports (PR #355, freeze `36ba9c2`): a real scalar on the unit
  sphere with one antipodal handle, glued torus-type or Klein-type, by a
  separable finite-element solver with #338's mesh as a cross-check.
  Verdict `COUNTEREXAMPLE`: the frozen split threshold `1e-6` fails at
  `a = 0.3`, `m = 3`, where the converged split is `1.01e-8` at both
  `N = 2000` and `4000`. The P5 identity also misses `1e-8` over runs of
  `ω̄t ~ 1e8`, and the miss does not shrink with `N`. Through the handle,
  `m = 0` cannot tell the gluings apart; torus gluing conserves `L`; the
  Klein spectrum is the union of two torus spectra; a packet arrives with
  its angular momentum reversed (to `3e-8`), and the reversal time is
  `π/(2Δ)` in all 36 pair runs. The mesh agrees to 0.18%. Classical, 2D,
  not GR. See `docs/alice_handle_ports.md`.

- Four-cell throat comparison (PR #353; original freeze `533f925`,
  review amendment `31cfc67`): the angular action selects circulation;
  at fixed shift, normal reversal blocks net face transfer. The two
  orientable controls admit transfer and the two non-orientable controls
  forbid it in this Hamiltonian family. 20 amended checks and 45 tests
  pass, including 64 generic trajectories and a strong A response reaching
  opening 1.065. A has one independent opening coordinate, not two
  independent exterior mouths. Actual coordinate-map, one-face-source
  and tighter energy checks supplement the unchanged original archive.
  Charge and Einstein dynamics remain underived.
  See `docs/antipodal_normal_throat.md`.

- Round 12, Klein time attachment (`bulk/klein_time.py`, freeze `bb6ff2c`,
  corrected after review of `5b79975`): a classical scalar on the
  time-nonorientable quotient of the Einstein static universe keeps one
  quadrature of every mode, a Lagrangian half of the cover's phase space,
  and two caps compose to a twist-dependent resonance. The frozen "no free
  field" reading is withdrawn as a cutoff artefact, with an explicit
  continuum counterexample; the time-symmetric response is shown to exist
  without the attachment; the complex structure is the cover's. Verdicts now
  separate `INCOMPLETE`, `UNRESOLVED_NUMERICALLY` and `COUNTEREXAMPLE`.
  28 checks, 80 tests. No quantum rule is derived.

- Round 13, localized responding Klein handle
  (`bulk/klein_handle_response.py`, freeze `5d15963`): on a lattice handle
  with a time-reversing seam and a solved waist coordinate, the
  complementary-quadrature structure survives nonlinear response exactly;
  the waist response supplies no sector pairing; for a charged field the
  conserved charge pairs the sectors only when the seam reverses both time
  and fibre. Indefinite pairing, prescribed potential, no probability rule.
  32 checks at two sizes, 26 tests.

- Round 14, responding electric field on the Klein handle
  (`bulk/klein_handle_electric.py`, freeze `4ab8bcd`, corrected after
  review of `562af71`): lattice scalar electrodynamics solved from one
  Hamiltonian with the charged field and the mouths. The sheet value of the
  sector pairing is the electric flux through the seam, by Gauss's law. It
  is driven by conjugation-odd initial data of two kinds: flux trapped around
  the handle loop, or a matter current with no flux. Review withdrew the
  claim that trapped flux is necessary, moved the tangent test onto the
  Gauss surface and into the reduced space, and marked the imaginary part of
  the field as gauge dependent. 45 checks at two sizes, 34 tests.
- Round 15, the matching map (`bulk/klein_handle_matching.py`, freeze
  `3d3ee38`): PR #129's antipodal map is a Hopf fibre translation, so on a
  bulk scalar it gives the fibre-preserving cell with no pairing; PR #279's
  transport acts on fibre-charged functions antilinearly with square
  `(-1)^k` and is `J` on the doublet; the two are attached to different
  spheres. The spinor doublet is admissible on the transport seam only on an
  antiperiodic doubled loop; with time reversal charge survives and total
  spin vanishes. The trapped-flux sign is independent of the Pin label,
  retracting a round-14 guess. 30 checks, 23 tests.

- Round 16, the two routes on equal terms (`bulk/klein_handle_routes.py`,
  freeze `fd78693`, corrected after review of `76c2067`): at matched total
  energy a matter current drives a larger sheet pairing than trapped flux,
  by a factor shrinking as the square root of energy at short times and
  about three over a long window, while the flux keeps its sign in every
  history. A matter current separates charge on every seam except the
  ordinary handle, with the parity in time fixed by the seam, so the pairing
  and the Lagrangian halving are distinct. Review withdrew the claim that
  fibre reversal is necessary, narrowed additivity to leading order in the
  matter amplitude, and replaced a median bound on flux retention by
  per-draw extrema. One frozen tolerance is contradicted and recorded.
  32 tests.

- Round 17, temporal symmetry of exchange (`bulk/klein_handle_parity.py`,
  freeze `95d7a80`): the parity in time of a conserved quantity's sheet
  value is fixed by whether the seam reverses time and whether its generator
  commutes with the deck map. On the transport seam with time reversal
  charge exchange is odd and spin exchange even. The field type and the time
  extension select it, not the geometry alone. 33 checks, 6 tests.

- Review of PR #352 at `c433f91`: round 15 no longer equates zero total
  charge with no sheet pairing (the fibre-preserving cell exchanges charge
  evenly in time); the Hopf map on the spatial brane mouth is generically
  two-to-one with a circle fibre over `-i`; round 17's exchange gate now
  requires variation; the continuum-gap diagnostic minimises over both
  neighbouring integers.

- Visualization (`viz/wave_mouth_network.py`): a figure integrated from one
  Hamiltonian in which a real wave moves relativistic mouths, triggers the
  transition of a loaded throat that releases its stored energy into the
  scalar field, and an energy ledger closes to `2.8e-05`. Comparator model; the metastable support is assumed.
  Follows PR #340, whose drawing has no propagating wave, collapse or
  ledger. Review of `89ceeed`: the collapse detector no longer fires on a
  stable throat, the 6.7x figure is named total field energy rather than
  emitted radiation, and the full-window DOP853 check (`7.8e-05`) is
  archived. Outward energy flux through fixed boundaries around the
  throat is measured from the exact lattice flux: net `3.315` by `t = 70`,
  balancing the energy lost inside to `1.2e-04`. The interaction is also
  drawn on the sphere with two Klein handles of #338/#340, from its own
  solved S² model (`waves/wave_mouth_sphere.py`, `viz/wave_mouth_sphere.py`),
  with the net flux out of a fixed cap around the loaded mouth measured by
  exact quadrature (`+1.167` at `t = 15`, non-monotonic on a closed sphere).
  Review of `324c2e4`: gross flows are split at each point of the circle
  (`3.318` out, `2.152` in), footprint-tail work closes the balance at
  second order, and an unclosed README code fence is fixed. 44 tests.

## 1.3.5 — release candidate (prepared 2026-10-29)

Capstone of PR #301–#350 since v1.3.0. See the
[data-generated release notes](docs/releases/bam-foundations-1.3.5.md).

- Collected verified compact-handle constraints, localized four-scalar initial
  data, Floquet and packet transport, R3 return-family and symmetry results,
  reciprocal scalar recoil, Maxwell throat response and Einstein radius evolution.
- Imported final PR #345 pair-birth, strict-router, conformal S³ and two-port code
  and tests at `cd70f0f`; retained the original import snapshot for historical audits.
- Added pinned evidence snapshots for research lines outside the `main` baseline.
- Bound headline measurements to JSON paths and formatting rules; added automatic
  release-table verification and a figure spanning six research lines.
- Aligned package metadata and planned release tag at 1.3.5.
- Added citation metadata, deterministic evidence packaging, wheel/checksum build
  instructions and an artifact-preparation workflow.
- Preserved all original archive verdicts, historical negative results and model
  limitations in the linked research record.

This candidate has not itself created a release tag, published a GitHub release
or assigned a new DOI.
