"""
Geometrodynamic QCD — flux-tube network solver and shelled-closure census.

Models QCD confinement as a network of 1D flux tubes (branches) joining
at junctions, with topology-dependent boundary conditions, bridge fields
for string breaking, and loop currents for hybrid excitations.  Also
hosts the v3 shelled-closure defect census for the hadronic-constituent
sector (``quark_spectrum``) and the composite-mass scaffold built on
top of it (``hadron_spectrum``).

Submodules
----------
constants         QCD-specific physical constants (σ, α_s, ℏc, etc.)
color             SU(3) color algebra: singlet checks, generator projection
bridge            BridgeField double-well, Cornell potential
network           Node, Branch, Junction, HadronicNetwork dataclasses
topology          Topology constructors: meson, baryon, glueball, hybrid, …
solver            HadronicNetworkSolver (Störmer–Verlet + SAT boundaries)
spectrum          MöbiusSpectrum, AlphaSStringAnsatz, ThroatBranchCrosswalk
diagnostics       LatticeStringTension, HybridModeShift, BridgeCalibrator
quark_spectrum    6×6 Hermitian Hamiltonian for the six hadronic-constituent
                  closure states (v3 spec, see docs/quark_axioms.md)
hadron_spectrum   composite-mass scaffold composing quark_spectrum with
                  the existing bridge/network/spectrum modules
"""

from geometrodynamics.qcd.network import (
    Node,
    Branch,
    Junction,
    HadronicNetwork,
)
from geometrodynamics.qcd.solver import HadronicNetworkSolver
from geometrodynamics.qcd.topology import (
    make_meson_tube,
    make_baryon_y_network,
    make_glueball_ring,
    make_mobius_tube,
    make_hybrid_excitation,
    make_tetraquark_double_y,
    make_mobius_baryon_y_network,
    make_mobius_baryon_v12,
)
from geometrodynamics.qcd.quark_spectrum import (
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
from geometrodynamics.qcd.hadron_spectrum import (
    BaryonConfig,
    MesonConfig,
    baryon_mass_mev,
    kaon_minus_config,
    meson_mass_mev,
    neutron_config,
    pion_minus_config,
    proton_config,
)
