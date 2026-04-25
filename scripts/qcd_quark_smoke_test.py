"""
Standalone smoke test for the v3-r2 code drop.

Tests:
  1. Hamiltonian Hermiticity (defaults)
  2. Hamiltonian Hermiticity (with partition_mixing)
  3. Block structure in lepton limit
  4. Eigenvalue pairing in lepton limit
  5. u_q(k) = k - 2 gives expected values
  6. Species ordering at k=1,3,5 matches BASIS_TO_SPECIES
  7. Color-independence structural check
  8. Uncalibrated solve raises helpful error
  9. (r2) extract_physical_spectrum produces all six species at unmixed limit
 10. (r2) extract_physical_spectrum produces all six species at strong mixing
 11. (r2) spectrum_zero shift gives positive masses at anchor
 12. (r2) adiabatic_tracking_check passes at moderate and strong mixing

Run from anywhere:
    python smoke_test.py

The test locates quark_spectrum.py relative to its own file location,
so it works regardless of where the drop is unpacked.
"""

from __future__ import annotations

import importlib.util
import math
import os
import sys


def _load_module():
    """
    Robustly locate quark_spectrum.py relative to this script.

    Looks in (in order):
      - ./geometrodynamics/qcd/quark_spectrum.py       (repo-relative layout)
      - ../geometrodynamics/qcd/quark_spectrum.py      (invoked from scripts/)
      - ./quark_spectrum.py                            (flat layout)
    """
    here = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.path.join(here, "geometrodynamics", "qcd", "quark_spectrum.py"),
        os.path.join(here, "..", "geometrodynamics", "qcd", "quark_spectrum.py"),
        os.path.join(here, "quark_spectrum.py"),
    ]
    for path in candidates:
        if os.path.isfile(path):
            spec = importlib.util.spec_from_file_location("quark_spectrum", path)
            mod = importlib.util.module_from_spec(spec)
            # Register before exec so dataclass introspection works on 3.12+
            sys.modules["quark_spectrum"] = mod
            spec.loader.exec_module(mod)
            print(f"Loaded quark_spectrum.py from: {path}")
            return mod
    raise FileNotFoundError(
        f"Could not locate quark_spectrum.py in any of: {candidates}"
    )


def main() -> int:
    import numpy as np
    qs = _load_module()

    failures = []

    # ── 1. Hermiticity (defaults) ────────────────────────────────────
    H = qs.build_quark_hamiltonian()
    if not np.allclose(H, H.conj().T, atol=1e-14):
        failures.append("Hamiltonian not Hermitian at defaults")
    else:
        print("[pass] Hermiticity (defaults)")

    # ── 2. Hermiticity (with partition mixing) ───────────────────────
    params_mixed = qs.QuarkParams(partition_mixing=0.5, phase=0.01)
    Hm = qs.build_quark_hamiltonian(params_mixed)
    if not np.allclose(Hm, Hm.conj().T, atol=1e-14):
        failures.append("Hamiltonian not Hermitian with partition_mixing=0.5")
    else:
        print("[pass] Hermiticity (partition_mixing != 0)")

    # ── 3. Block structure in lepton limit ──────────────────────────
    limit_params = qs.QuarkParams(
        action_base=2.0 * math.pi,
        gamma_q=0.0,
        partition_mixing=0.0,
    )
    Hl = qs.build_quark_hamiltonian(limit_params)
    plus_idx = [0, 2, 4]
    minus_idx = [1, 3, 5]
    max_cross = 0.0
    for i in plus_idx:
        for j in minus_idx:
            max_cross = max(max_cross, abs(Hl[i, j]))
    if max_cross > 1e-14:
        failures.append(
            f"Cross-partition entries nonzero in lepton limit: max = {max_cross:.2e}"
        )
    else:
        print("[pass] Block-diagonal structure in lepton limit")

    # ── 4. Eigenvalue pairing in lepton limit ───────────────────────
    eigs = np.linalg.eigvalsh(Hl)
    max_pair_err = 0.0
    for i in range(3):
        max_pair_err = max(max_pair_err, abs(eigs[2 * i + 1] - eigs[2 * i]))
    if max_pair_err > 1e-12:
        failures.append(
            f"Eigenvalues don't pair up in lepton limit: max pair err = {max_pair_err:.2e}"
        )
    else:
        print(f"[pass] Eigenvalue pairing in lepton limit (max err = {max_pair_err:.2e})")

    # ── 5. u_q(k) = k - 2 ordering ───────────────────────────────────
    u_vals = [qs._u_q(k, "k_minus_2") for k in qs.PASS_COUNTS]
    if u_vals != [-1.0, 1.0, 3.0]:
        failures.append(f"u_q values wrong: {u_vals}")
    else:
        print(f"[pass] u_q(k) = k - 2 gives {u_vals}")

    # ── 6. Ordering consistent with BASIS_TO_SPECIES ────────────────
    params_split = qs.QuarkParams(gamma_q=1.0, partition_mixing=0.0)
    Hs = qs.build_quark_hamiltonian(params_split)
    diag = np.real(np.diag(Hs))
    idx = {state: i for i, state in enumerate(qs.BASIS_STATES)}

    if diag[idx[(1, "+")]] >= diag[idx[(1, "-")]]:
        failures.append("At k=1, (1,+) is not lighter than (1,-)")
    else:
        print("[pass] k=1 ordering: (1,+) < (1,-)  [u < d]")
    if diag[idx[(3, "-")]] >= diag[idx[(3, "+")]]:
        failures.append("At k=3, (3,-) is not lighter than (3,+)")
    else:
        print("[pass] k=3 ordering: (3,-) < (3,+)  [s < c]")
    if diag[idx[(5, "-")]] >= diag[idx[(5, "+")]]:
        failures.append("At k=5, (5,-) is not lighter than (5,+)")
    else:
        print("[pass] k=5 ordering: (5,-) < (5,+)  [b < t]")

    # ── 7. Color check is structural ────────────────────────────────
    cc = qs.color_independence_check()
    if cc["passed"] and cc["color_multiplicity"] == 3:
        print("[pass] Color-independence structural check")
    else:
        failures.append(f"Color check failed: {cc}")

    # ── 8. Calibration-state gate ────────────────────────────────────
    # Behaves differently depending on whether the four-step pipeline
    # has populated LOCKED_QUARK_PARAMS yet:
    #   - unlocked (None)  → solved_quark_masses_mev must raise
    #                        NotImplementedError with a helpful message.
    #   - locked           → solved_quark_masses_mev must return a
    #                        length-6 ndarray in QUARK_SPECIES order.
    if qs.LOCKED_QUARK_PARAMS is None:
        try:
            qs.solved_quark_masses_mev()
            failures.append("solved_quark_masses_mev should have raised")
        except NotImplementedError as exc:
            if "not yet calibrated" in str(exc):
                print("[pass] Uncalibrated solve raises helpful "
                      "NotImplementedError")
            else:
                failures.append(f"Wrong error message: {exc}")
    else:
        try:
            masses = qs.solved_quark_masses_mev()
            # Two valid regimes: u-anchor (action_base zero, masses[0]
            # is the u observed mass) and d-anchor (min_eigenvalue
            # zero, masses[0] = 0 and masses[1] is the d observed mass).
            mode = qs.LOCKED_QUARK_PARAMS.spectrum_zero_mode
            ok_u_anchor = (
                mode == "action_base"
                and masses.shape == (6,)
                and masses[0] > 0
            )
            ok_d_anchor = (
                mode == "min_eigenvalue"
                and masses.shape == (6,)
                and masses[0] == 0
                and masses[1] > 0
            )
            if ok_u_anchor:
                print(f"[pass] Locked solve returns 6 masses "
                      f"(u-anchor, u = {masses[0]:.4f} MeV)")
            elif ok_d_anchor:
                print(f"[pass] Locked solve returns 6 masses "
                      f"(d-anchor, d = {masses[1]:.4f} MeV)")
            else:
                failures.append(f"solved_quark_masses_mev returned "
                                f"unexpected shape/anchor: {masses}")
        except Exception as exc:
            failures.append(f"Locked solve raised unexpectedly: {exc}")

    # ── 9. (r2) Extract spectrum at unmixed limit ───────────────────
    try:
        params_unmixed = qs.QuarkParams(gamma_q=0.0, partition_mixing=0.0)
        spec_unmixed = qs.extract_physical_spectrum(params_unmixed)
        if set(spec_unmixed.keys()) == set(qs.QUARK_SPECIES):
            print(f"[pass] extract_physical_spectrum at unmixed limit: all 6 species")
        else:
            failures.append(
                f"Unmixed-limit species set wrong: {set(spec_unmixed.keys())}"
            )
    except Exception as exc:
        failures.append(f"extract at unmixed limit raised: {exc}")

    # ── 10. (r2) Extract spectrum at strong mixing ──────────────────
    try:
        params_strong = qs.QuarkParams(gamma_q=5.0, partition_mixing=2.0)
        spec_strong = qs.extract_physical_spectrum(
            params_strong, n_adiabatic_steps=32,
        )
        if set(spec_strong.keys()) == set(qs.QUARK_SPECIES):
            print(f"[pass] extract_physical_spectrum at strong mixing: all 6 species")
        else:
            failures.append(
                f"Strong-mixing species set wrong: {set(spec_strong.keys())}"
            )
    except Exception as exc:
        failures.append(f"extract at strong mixing raised: {exc}")

    # ── 11. (r2) Spectrum zero produces strictly positive anchor mass ──
    # Physics: mass = eigenvalue - action_base.  At small γ_q and small
    # transport, every diagonal entry of the Hamiltonian exceeds
    # action_base, and level repulsion stays bounded, so every shifted
    # mass is strictly positive.  At large γ_q or transport, level
    # depression can push the lowest state below action_base — the
    # calibration pipeline rejects such points as unphysical.
    #
    # This check verifies that the positive regime exists.
    try:
        light_params = qs.QuarkParams(
            gamma_q=0.1,            # small, keeps (1,+) above action_base
            partition_mixing=0.05,  # small
            transport=1.0,          # small, keeps level repulsion mild
            phase=0.001,
            pinhole=22.5,
            resistance=0.217869,
        )
        spec = qs.extract_physical_spectrum(light_params)
        anchor_val = spec[qs.QUARK_ANCHOR_SPECIES]
        POSITIVITY_FLOOR = 1e-3
        if anchor_val > POSITIVITY_FLOOR:
            print(f"[pass] Spectrum-zero shift gives positive anchor mass "
                  f"({qs.QUARK_ANCHOR_SPECIES} = {anchor_val:.4f} > {POSITIVITY_FLOOR})")
            # Report the full spectrum as diagnostic
            print(f"       Full spectrum (path-length excess over action_base):")
            for s in qs.QUARK_SPECIES:
                print(f"         {s}: {spec[s]:.4f}")
        else:
            failures.append(
                f"Anchor mass not strictly positive even at small-γ point: "
                f"{qs.QUARK_ANCHOR_SPECIES} = {anchor_val:.4e}.  "
                f"This suggests a structural issue with the spectrum-zero choice."
            )
    except Exception as exc:
        failures.append(f"Spectrum-zero check raised: {exc}")

    # ── 11b. (r2) Unphysical parameters trigger anchor rejection ────
    # The pipeline must reject parameter points where the anchor sits
    # at or below action_base.  Test that large γ_q triggers the check.
    try:
        unphysical = qs.QuarkParams(
            gamma_q=10.0,  # large, depresses (1,+) well below action_base
            partition_mixing=0.0,
            transport=1.0,
        )
        spec_un = qs.extract_physical_spectrum(unphysical)
        anchor_un = spec_un[qs.QUARK_ANCHOR_SPECIES]
        if anchor_un <= 0:
            print(f"[pass] Unphysical parameter point correctly yields "
                  f"non-positive anchor ({qs.QUARK_ANCHOR_SPECIES} = {anchor_un:.4f})")
        else:
            failures.append(
                f"Expected non-positive anchor at γ_q=10, got {anchor_un:.4f}"
            )
    except Exception as exc:
        failures.append(f"Unphysical-point check raised: {exc}")

    # ── 12. (r2) Adiabatic tracking check ───────────────────────────
    try:
        at = qs.adiabatic_tracking_check(
            gamma_q=5.0, partition_mixing=2.0, n_steps=32,
        )
        if at["passed"]:
            print(f"[pass] Adiabatic tracking at (γ=5, w=2) with 32 steps")
        else:
            failures.append(f"Adiabatic tracking failed: {at.get('error')}")
    except Exception as exc:
        failures.append(f"Adiabatic check raised: {exc}")

    # ── Summary ──────────────────────────────────────────────────────
    print()
    if failures:
        print(f"SMOKE TEST FAILED ({len(failures)} failures):")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("SMOKE TEST PASSED — r2 drop is internally consistent.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
