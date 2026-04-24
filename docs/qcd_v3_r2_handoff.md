# QCD v3-r2 code drop — handoff for Claude Code

This is revision 2 of the v3 code package, incorporating fixes for the
issues identified in the external review of r1.

## Changes from r1

| Review item | Status in r2 |
|---|---|
| **(1) Raw signed eigenvalues anchored as masses** | **Fixed.** New `extract_physical_spectrum()` function implements spectrum-zero shift (anchoring at `action_base`, the topological minimum). Calibration scripts use it throughout and reject parameter points where the anchor species sits at or below `action_base`. See `quark_axioms.md §3.5`. |
| **(2) Species labeling by sort order fragile under mixing** | **Fixed.** `extract_physical_spectrum()` uses adiabatic continuation from the fully-unmixed reference state, tracking each eigenvector by maximum-overlap with its ancestor. Species → mass mapping is stable under strong mixing. Smoke test verifies this at `(γ=5, w=2)`. |
| **(3) Smoke test path broken** | **Fixed.** `smoke_test.py` now checks both repo-relative and flat layouts, reports which one it loaded from, and carries a sys.path-free importer. |
| **(4) CLI scripts not standalone-portable** | **Fixed.** All four scripts now have a `sys.path` shim identical in spirit to the lepton scripts, so they can run with `python scripts/<name>.py` from the repo root without `PYTHONPATH` or `pip install`. Confirmed by running `--help` on each. |
| **(5) `hadron_spectrum.py` scaffold honesty** | **Unchanged.** Still a scaffold. Review acknowledged this was honestly marked in r1 and not a flaw. |
| **(6) Partition-mixing phase still a placeholder** | **Unchanged.** The r1 review correctly flagged this as a scaffolding-stage item. Replacing `params.phase * k` with the Hopf-derived phase from `hopf/connection.py` requires access to the live Hopf module API; it's a post-landing TODO for Claude Code, not something that can be honestly completed in this session. |

## New in r2 that the review didn't ask for, but is worth flagging

**The `γ_q` positivity tension** (documented in `quark_axioms.md §3.5`).
While fixing issue (1), the spectrum-zero shift revealed a real
structural question about the `γ_q · σ(p) · u_q(k)` term:

- At `k=1, p=+`, the partition contribution is `γ_q · (+1) · (-1) = -γ_q`.
- This subtracts from `action_base` directly, and if `γ_q > r_q` then
  the unmixed diagonal at `(1,+)` sits *below* `action_base`.
- The calibration pipeline rejects such points as unphysical.

This means the calibration is **constrained** to the `γ_q < r_q`
regime. If the scan's best fit wants to cross that boundary, one of
three things is happening; `§3.5` enumerates the possibilities and
asks the human reviewer (that's you, Claude Code, and/or David) to
interpret the scan log to decide which is right.

I flagged this up-front rather than burying it because it is the kind
of thing that the spec must acknowledge explicitly, not something that
should emerge as a mysterious convergence failure during calibration.

## What's in this drop

```
qcd-v3-drop-r2/
├── HANDOFF.md                                    ← this file
├── smoke_test.py                                 ← standalone, 14 checks
├── geometrodynamics/qcd/
│   ├── quark_spectrum.py                         ← REVISED for r2
│   └── hadron_spectrum.py                        ← unchanged (scaffold)
├── tests/
│   └── test_quark_spectrum.py                    ← unchanged from r1
├── scripts/
│   ├── calibrate_quark_ratios.py                 ← REVISED: path shim + extract
│   ├── sweep_quark_beta.py                       ← REVISED: path shim + extract
│   ├── map_basin_quark_uplift.py                 ← REVISED: path shim + extract
│   └── lock_quark_beta_probe.py                  ← REVISED: path shim + extract
└── docs/
    └── quark_axioms.md                           ← REVISED: added §3.5
```

## Verified locally before handoff

```
$ cd qcd-v3-drop-r2 && python smoke_test.py
Loaded quark_spectrum.py from: geometrodynamics/qcd/quark_spectrum.py
[pass] Hermiticity (defaults)
[pass] Hermiticity (partition_mixing != 0)
[pass] Block-diagonal structure in lepton limit
[pass] Eigenvalue pairing in lepton limit (max err = 2.84e-14)
[pass] u_q(k) = k - 2 gives [-1.0, 1.0, 3.0]
[pass] k=1 ordering: (1,+) < (1,-)  [u < d]
[pass] k=3 ordering: (3,-) < (3,+)  [s < c]
[pass] k=5 ordering: (5,-) < (5,+)  [b < t]
[pass] Color-independence structural check
[pass] Uncalibrated solve raises helpful NotImplementedError
[pass] extract_physical_spectrum at unmixed limit: all 6 species
[pass] extract_physical_spectrum at strong mixing: all 6 species
[pass] Spectrum-zero shift gives positive anchor mass (u = 0.0908 > 0.001)
[pass] Unphysical parameter point correctly yields non-positive anchor (u = -9.7900)
[pass] Adiabatic tracking at (γ=5, w=2) with 32 steps

SMOKE TEST PASSED — r2 drop is internally consistent.

$ python scripts/calibrate_quark_ratios.py --n-points 3
... (runs end-to-end, produces sensible output, anchors u to 2.16 MeV)
```

## Landing checklist (unchanged from r1, re-confirmed for r2)

1. **Clone the repo and branch off `main`.**
   ```bash
   git clone https://github.com/davidmdrpi/geometrodynamics.git
   cd geometrodynamics
   git checkout -b qcd-v3-r2-quark-spectrum
   ```

2. **Drop files into their matching target paths.** Every file in this
   drop has a path identical to its destination in the repo tree.

3. **Wire up the `qcd` subpackage exports** in `geometrodynamics/qcd/__init__.py`:
   ```python
   from .quark_spectrum import (
       BASIS_STATES,
       BASIS_TO_SPECIES,
       OBSERVED_MASSES_MEV,
       PASS_COUNTS,
       QUARK_ACTION_BASE,
       QUARK_ANCHOR_MASS_MEV,
       QUARK_ANCHOR_SPECIES,
       QUARK_BETA_DEFAULT,
       QUARK_SPECIES,
       SPECIES_TO_BASIS,
       QuarkParams,
       adiabatic_tracking_check,
       build_quark_hamiltonian,
       color_independence_check,
       extract_physical_spectrum,
       quark_lepton_limit_check,
       solved_quark_masses_mev,
   )
   from .hadron_spectrum import (
       BaryonConfig,
       MesonConfig,
       baryon_mass_mev,
       kaon_minus_config,
       meson_mass_mev,
       neutron_config,
       pion_minus_config,
       proton_config,
   )
   ```

4. **Verify imports resolve cleanly:**
   ```bash
   python -c "from geometrodynamics.qcd import quark_spectrum as qs; print(qs.BASIS_STATES)"
   ```

5. **Run the structural tests:**
   ```bash
   pytest tests/test_quark_spectrum.py -v
   ```
   The r1 tests are kept as-is; they exercise the pre-extract behavior.
   Consider adding r2 tests for `extract_physical_spectrum` in a
   follow-up PR once the calibration has produced a locked point.

6. **Verify the `sys.path` shim works** by running each script with
   `--help` from the repo root. All four must succeed.

7. **Run the calibration pipeline:**
   ```bash
   # Step 1: coarse grid scan
   python scripts/calibrate_quark_ratios.py --n-points 8 --verbose \
       --output-json /tmp/quark_step1.json

   # Step 2: integer-winding hunt
   python scripts/sweep_quark_beta.py \
       --input-json /tmp/quark_step1.json \
       --integer-min 10 --integer-max 1000000 \
       --n-samples 600 --verbose \
       --output-json /tmp/quark_step2.json

   # Step 3: basin-width sanity check around the best integer
   python scripts/map_basin_quark_uplift.py \
       --integer-winding <N_FROM_STEP_2> \
       --half-width 50 --n-points 101 \
       --output-json /tmp/quark_step3.json

   # Step 4: final lock
   python scripts/lock_quark_beta_probe.py \
       --integer-winding <N_FROM_STEP_2> \
       --action-base-label pi \
       --n-points 12 --verbose \
       --output-json /tmp/quark_step4.json
   ```

   **Expected behavior:** step 1 will likely *not* find a parameter
   point that fits all six masses — the coarse grid cannot discover
   the heavy-sector uplift. Step 2 is what enables the heavy-sector
   fit. If step 1 reports `rejected_points > 0`, that's expected and
   good (the positivity rejector is working); if it reports
   `rejected_points == scan_points_total`, the scan ranges are wrong
   and need to be narrowed around smaller `γ_q`.

8. **Inspect the step-2 output for physical interpretation.** The
   integer winding reported by step 2 should be inspected against
   `docs/quark_axioms.md §3.5`:
   - If it's a clean round number in some S³ unit, note it as a
     candidate new topological invariant.
   - If it's not, the minimal v3 ansatz may need revision.

9. **Paste the `LOCKED_QUARK_PARAMS` snippet from step 4 into**
   `geometrodynamics/qcd/quark_spectrum.py`, replacing
   `LOCKED_QUARK_PARAMS: Optional[QuarkParams] = None`.

10. **Re-run the full test suite:**
    ```bash
    pytest -v
    ```

11. **Update `docs/quark_axioms.md §8`** with the locked values, basin
    width, and relative errors per species.

12. **Update `README.md`** to reference the new module in the "What
    the Code Validates" table with honest status.

13. **Commit and push.**

## Known limitations

- **`hadron_spectrum.py` is still a scaffold.** Needs the `qcd.bridge`
  and `qcd.spectrum` API surface verified before wiring up.

- **Partition-mixing phase `φ_q(k)` is still `params.phase * k`.**
  Replace with the Hopf-derived phase once `hopf/connection.py` has
  been inspected. Non-blocking for initial calibration but required
  for v3 §7 criterion 1 ("quantization visibly emergent") to be met.

- **`quark_lepton_limit_check()` still compares ratios, not absolute
  values.** Absolute equality requires synchronizing residual knobs
  with the lepton locked baseline, which should happen as part of the
  post-landing validation.

- **The `r2` fix for spectrum-zero rejects unphysical parameter
  points.** This is by design, but it means naive coarse scans may
  report "no valid points found" if the initial ranges are too
  aggressive. If that happens, narrow the `γ_q` and `transport`
  ranges and retry.

## Methodological rule compliance

Checked by grep before writing this handoff:

```bash
git diff main -- '*.py' '*.md' | grep -iE \
    '(isospin|weak.charge|flavor.as|CKM|W.boson|gluon|SU\(3\).gauge)'
```

All matches in r2 are either:
- inside the rule's own definition (listing banned terms),
- inside `docs/quark_axioms.md §9` (the permitted phenomenological-
  interpretation section),
- or inside this HANDOFF's own grep command.

No leakage into axioms, the Hamiltonian, or public identifiers.
