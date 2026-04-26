"""
geometrodynamics.qcd.hadron_spectrum
=====================================

Composite-mass layer for the shelled-closure framework.

Composes defect closure energies from ``quark_spectrum`` with the existing
flux-tube bridge from ``qcd.bridge``, the Y-junction / linear-tube network
from ``qcd.network``, and the Möbius shell harmonics from ``qcd.spectrum``
to produce composite masses for mesons, baryons, glueballs, and hybrids.

METHODOLOGICAL RULE (binding, inherited from v3 spec §0.5)
───────────────────────────────────────────────────────────
Same neutrality rule as ``quark_spectrum``: only geometric vocabulary in
the axioms and the composition recipes. "Meson", "baryon", "glueball",
"hybrid" are used because the existing ``qcd.topology`` module already
uses these names as structural labels (not as gauge-theoretic claims),
and this module inherits that convention.

GOVERNING PRINCIPLE (v3 spec §6, §7 criteria 6 and 7)
──────────────────────────────────────────────────────
Composite masses decompose as:

    M_composite = Σ_defects E_closure(defect)      # from quark_spectrum
                + E_bridge(network_configuration)  # from qcd.bridge
                + E_shell(Möbius modes)            # from qcd.spectrum

The defect-level contribution is typically a small fraction of the total
for ground-state composites built from the lightest species; the bulk
comes from bridge static-tension and shell binding. This is the quark-
sector analog of the current-vs-constituent mass decomposition, now
distributed across three modules.

STATUS
───────
Scaffold only. The concrete composition functions (``meson_mass``,
``baryon_mass``, etc.) are stubs pending:

    1. Calibration of ``quark_spectrum.LOCKED_QUARK_PARAMS`` (blocking).
    2. Verification of the existing ``qcd.bridge`` / ``qcd.network`` API
       surface against the call sites sketched below.

Once those are in place, the stubs become straightforward delegation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

# Imports deferred to function bodies so that the module remains importable
# even before LOCKED_QUARK_PARAMS is populated.


# ════════════════════════════════════════════════════════════════════════
# COMPOSITE DEFINITIONS
# ════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class MesonConfig:
    """
    Minimal meson configuration: two defects connected by a single bridge.

    Parameters
    ----------
    species_a, species_b : str
        Keys into ``quark_spectrum.QUARK_SPECIES``. Partition classes of
        the two defects are required to be opposite for ground-state
        mesons (closure requires p_a + p_b = 0 on the flux-partition
        sum rule).
    bridge_length : float
        Static flux-tube length in the natural units of ``qcd.bridge``.
    """
    species_a: str
    species_b: str
    bridge_length: float = 1.0


@dataclass(frozen=True)
class BaryonConfig:
    """
    Minimal baryon configuration: three defects at a Y-junction.

    Parameters
    ----------
    species_trio : tuple[str, str, str]
        Three keys into ``quark_spectrum.QUARK_SPECIES``. The color-triplet
        closure constraint (degeneracy sum = 0 mod 3) is enforced by
        ``qcd.color`` at composition time.
    junction_arm_length : float
        Length of each of the three Y-junction arms.
    """
    species_trio: Tuple[str, str, str]
    junction_arm_length: float = 1.0


# ════════════════════════════════════════════════════════════════════════
# MASS COMPOSITION (§6)
# ════════════════════════════════════════════════════════════════════════

def _defect_contribution_mev(species: str) -> float:
    """
    Total closure energy of a single defect in MeV.

    Reads from ``quark_spectrum.solved_quark_masses_mev`` once calibrated.
    Raises ``NotImplementedError`` if the spectrum is not yet locked.
    """
    from .quark_spectrum import (
        QUARK_SPECIES,
        solved_quark_masses_mev,
    )
    masses = solved_quark_masses_mev()
    if species not in QUARK_SPECIES:
        raise ValueError(f"Unknown species: {species!r}. Expected one of {QUARK_SPECIES}.")
    return float(masses[QUARK_SPECIES.index(species)])


def _bridge_contribution_mev(config: object) -> float:
    """
    Static bridge energy for a composite configuration.

    Delegates to ``qcd.bridge``. Scaffold: the actual call depends on the
    bridge API, which should be verified before wiring up. The signature
    sketched here assumes a BridgeField-like object with a static-energy
    method.
    """
    # TODO: wire up to the actual qcd.bridge API. Placeholder:
    #
    #   from .bridge import BridgeField
    #   bridge = BridgeField(...)
    #   return bridge.static_energy_mev(config)
    #
    raise NotImplementedError(
        "Bridge-contribution wiring pending verification against qcd.bridge API. "
        "See hadron_spectrum.py TODO in _bridge_contribution_mev."
    )


def _shell_binding_contribution_mev(config: object) -> float:
    """
    Shell-harmonic binding energy from ``qcd.spectrum``.

    The Möbius half-integer modes of the existing ``qcd.spectrum`` module
    ARE the interior shell harmonics of the v3 §1 picture. This function
    sums the lowest stable Möbius-mode occupations consistent with the
    composite configuration.
    """
    # TODO: wire up to qcd.spectrum. Placeholder:
    #
    #   from .spectrum import mobius_modes, stable_occupation
    #   modes = mobius_modes(config)
    #   return float(sum(stable_occupation(modes)))
    #
    raise NotImplementedError(
        "Shell-binding wiring pending verification against qcd.spectrum API. "
        "See hadron_spectrum.py TODO in _shell_binding_contribution_mev."
    )


def meson_mass_mev(config: MesonConfig) -> float:
    """
    Total mass of a meson-configured composite, in MeV.

    Decomposition (v3 §6):
        M = E_defect(a) + E_defect(b) + E_bridge + E_shell

    For the lightest-pseudoscalar test of §7 criterion 6, the sum of
    defect contributions is ~10 MeV while the total is ~130 MeV; the
    remaining ~120 MeV comes from the bridge + shell binding and must
    be species-independent within the u/d sector. That species-
    independence is the primary test of the shell-as-chamber picture.
    """
    e_a = _defect_contribution_mev(config.species_a)
    e_b = _defect_contribution_mev(config.species_b)
    e_bridge = _bridge_contribution_mev(config)
    e_shell = _shell_binding_contribution_mev(config)
    return e_a + e_b + e_bridge + e_shell


def baryon_mass_mev(config: BaryonConfig) -> float:
    """
    Total mass of a baryon-configured (three-defect Y-junction) composite.

    Decomposition (v3 §6):
        M = Σ_{i=1..3} E_defect(species_i) + E_Y-bridge + E_shell

    For the proton/neutron test of §7 criterion 7, the defect sum is
    ~10 MeV while the total is ~940 MeV; ~930 MeV from the Y-junction
    bridge and shell binding.
    """
    e_defects = sum(_defect_contribution_mev(s) for s in config.species_trio)
    e_bridge = _bridge_contribution_mev(config)
    e_shell = _shell_binding_contribution_mev(config)
    return e_defects + e_bridge + e_shell


# ════════════════════════════════════════════════════════════════════════
# STANDARD-COMPOSITE CONSTRUCTORS
# ════════════════════════════════════════════════════════════════════════

def pion_minus_config(bridge_length: float = 1.0) -> MesonConfig:
    """Pseudoscalar composite from first-pass-count sector, opposite partitions."""
    return MesonConfig(species_a="d", species_b="u", bridge_length=bridge_length)


def kaon_minus_config(bridge_length: float = 1.0) -> MesonConfig:
    """Strange-light pseudoscalar composite."""
    return MesonConfig(species_a="s", species_b="u", bridge_length=bridge_length)


def proton_config(arm_length: float = 1.0) -> BaryonConfig:
    """Three-defect composite: two (1,+) and one (1,-) at Y-junction."""
    return BaryonConfig(species_trio=("u", "u", "d"), junction_arm_length=arm_length)


def neutron_config(arm_length: float = 1.0) -> BaryonConfig:
    """Three-defect composite: one (1,+) and two (1,-) at Y-junction."""
    return BaryonConfig(species_trio=("u", "d", "d"), junction_arm_length=arm_length)
