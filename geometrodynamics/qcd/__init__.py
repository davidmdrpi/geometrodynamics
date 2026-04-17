"""
Geometrodynamic QCD — flux-tube network solver.

Models QCD confinement as a network of 1D flux tubes (branches) joining
at junctions, with topology-dependent boundary conditions, bridge fields
for string breaking, and loop currents for hybrid excitations.

Submodules
----------
constants    QCD-specific physical constants (σ, α_s, ℏc, etc.)
color        SU(3) color algebra: singlet checks, generator projection
bridge       BridgeField double-well, Cornell potential
network      Node, Branch, Junction, HadronicNetwork dataclasses
topology     Topology constructors: meson, baryon, glueball, hybrid, …
solver       HadronicNetworkSolver (Störmer–Verlet + SAT boundaries)
spectrum     MöbiusSpectrum, AlphaSStringAnsatz, ThroatBranchCrosswalk
diagnostics  LatticeStringTension, HybridModeShift, BridgeCalibrator
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
