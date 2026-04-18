"""
Geometrodynamics
================

A research framework implementing and testing Wheeler's geometrodynamic
program: exploring whether QED structures (charge, spin-½, Coulomb law),
QCD structures (confinement, string breaking), black-hole structure
(nonsingular interior, entropy, charge), and Bell correlations can
emerge from pure spacetime geometry — specifically from non-orientable
wormhole-throat topology and the Hopf fibration on S³.

Subpackages
-----------
hopf          Hopf fibration on S³: connection, curvature, Chern number, spinors
tangherlini   5D wormhole eigenmodes, Maxwell solver, throat flux ratios
transaction   Retrocausal Wheeler–Feynman handshake protocol on S³ + cavity
embedding     Non-orientable throat topology, transport derivation from Hopf geometry
bell          Bell correlations from throat transport + Hopf SU(2) projection
history       Closed-history framework: unifying Bell, QED, and BH through closure
qcd           QCD flux-tube network: topology, Störmer–Verlet solver, diagnostics
blackhole     Black holes as coherent wormhole-throat condensates
viz           Visualisation helpers (animation, plotting)
"""

__version__ = "0.43.0"
