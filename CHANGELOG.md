# Changelog

## Unreleased

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
