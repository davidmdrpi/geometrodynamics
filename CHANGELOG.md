# Changelog

## Unreleased

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
  (`bulk/klein_handle_electric.py`, freeze `4ab8bcd`): lattice scalar
  electrodynamics solved from one Hamiltonian with the charged field and the
  mouths. Response alone never couples the two quadrature sectors; flux
  trapped around the handle loop does, with orientation set by its sign, and
  is admissible only where the charge exists; the sheet pairing is the seam
  flux by Gauss's law. The flux is initial data. Expected screening was not
  found. 27 checks at two sizes, 24 tests.

- Round 15, the matching map (`bulk/klein_handle_matching.py`, freeze
  `3d3ee38`): PR #129's antipodal map is a Hopf fibre translation, so on a
  bulk scalar it gives the fibre-preserving cell with no pairing; PR #279's
  transport acts on fibre-charged functions antilinearly with square
  `(-1)^k` and is `J` on the doublet; the two are attached to different
  spheres. The spinor doublet is admissible on the transport seam only on an
  antiperiodic doubled loop; with time reversal charge survives and total
  spin vanishes. The trapped-flux sign is independent of the Pin label,
  retracting a round-14 guess. 30 checks, 23 tests.

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
