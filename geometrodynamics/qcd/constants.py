"""
QCD-specific physical constants.

All dimensionful quantities are in GeV (energies/masses) and fm (lengths),
with ℏc = 0.197327 GeV·fm for conversion.
"""

from math import pi, sqrt
from typing import Literal

HBAR_C = 0.197327          # GeV·fm
SIGMA_QCD = 0.18           # string tension (GeV²)
SIGMA_FM = SIGMA_QCD / HBAR_C  # string tension (GeV/fm)
A_COULOMB = 0.30           # Coulomb coefficient in Cornell potential
L_BREAK_LAT = 1.35         # lattice string-breaking length (fm)

# ── Cornell amplitude calibration ────────────────────────────────────────────
R_HAD_FM = 0.84            # hadronic size anchor (fm)
L_STRING_ONSET_FM = 1.10   # flux-tube onset scale (fm)
AMP_EQ_MIN = 0.26          # hadronic anchor amplitude
AMP_EQ_POWER = 2.2

# ── Schwinger tunnelling ────────────────────────────────────────────────────
M_Q_SCHWINGER_GEV = 0.300  # constituent quark mass for Schwinger rate

# ── Bridge field ─────────────────────────────────────────────────────────────
BRIDGE_ALPHA = -0.50
BRIDGE_BETA = 0.25
BRIDGE_GAMMA = 0.20
BRIDGE_THRESHOLD = 1.0

# ── Confinement masses ──────────────────────────────────────────────────────
MU0_CONF_SQ = 0.50
MUJ_CONF_SQ = 0.30
KAPPA_J = 0.10

# ── SAT boundary treatment ──────────────────────────────────────────────────
SAT_PENALTY = 0.35
SAT_DAMPING = 5.0
SAT_HYBRID_SCALE = 0.29
SAT_RING_SCALE = 0.25

# ── Wormhole embedding ──────────────────────────────────────────────────────
R_MID = 1.00
R_OUTER = 1.26

# ── Type aliases ─────────────────────────────────────────────────────────────
TopologyKindBranch = Literal["open", "periodic", "mobius", "attached_loop"]
TopologyKind = Literal[
    "meson", "baryon", "glueball", "hybrid", "mobius",
    "mobius_baryon", "tetraquark",
]
BranchRole = Literal[
    "quark_arm", "ring", "hybrid_link", "mobius_arm", "junction_bridge",
]
