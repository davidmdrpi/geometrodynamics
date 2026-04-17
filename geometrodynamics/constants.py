"""
Physical and simulation constants for the geometrodynamic QED model.

These parameterise the S³ spatial topology, the 5D Tangherlini wormhole,
and the retrocausal transaction protocol.
"""

import numpy as np

# ── Wormhole geometry ────────────────────────────────────────────────────────
R_MID = 1.00        # throat radius (geometric units)
DELTA = 0.26        # half-width of the radial domain
R_OUTER = R_MID + DELTA
R_INNER = R_MID - DELTA

# ── Geon / wave parameters ──────────────────────────────────────────────────
M_GEON = 0.050
SIGMA_G = 0.40
SIGMA_R = 0.15
THETA_MAX = np.pi / 2.0
C_GW = 0.52          # gravitational-wave propagation speed on S³
WAVE_W = 0.38
WAVE_SEP = 0.10
A0_RING = 0.44 * DELTA

# ── Animation timing ────────────────────────────────────────────────────────
REFORM_OFFSET = 0.20
GROW_START = 0.82 * np.pi
RECONNECT_THRESHOLD = 0.82
SPIN_TWIST = np.pi
SPIN_ELLIPSE = 2.0
OMEGA_SPIN = 0.42
CONE_ALPHA = 0.28
TUBE_R = 0.040
TUBE_LEN = 0.22
TUBE_N = 18
TUBE_L = 8
RIBBON_N = 60
FPS = 30
T_HALF = 5.5
T_EXCHANGE = np.pi / C_GW
T_WAVE = np.pi / C_GW
T_DETACH = 0.70
T_RECONNECT = 0.40
T_CYCLE = T_EXCHANGE + T_HALF + T_DETACH + T_WAVE + T_RECONNECT + T_HALF
FRAMES = int(round(2.0 * T_CYCLE * FPS))

# ── Mode coupling & solver ──────────────────────────────────────────────────
GAMMA_MODE = 0.08
DT_DEFAULT = 0.015

# ── S³ cavity & transaction protocol ────────────────────────────────────────
S3_RADIUS = 1.0
S3_GREEN_EPS = 0.08
EPS_HIT = 0.06
SIGMA_ANTI = 0.18
SIGMA_SRC = 0.025
PHASE_MATCH_SIGMA = 0.60
PHASE_MATCH_MAX = 1.20
OFFER_TTL = np.pi / C_GW + 2.0
CONFIRM_TTL = 1.35

# ── Cavity mode oscillator ──────────────────────────────────────────────────
CAVITY_GAMMA = 0.018
CAVITY_BMIN = 0.015
CAVITY_LAMBDA = 0.012
CAVITY_SOFT_COUP = 0.002
CAVITY_ALPHA_EMIT = 0.60
CAVITY_ALPHA_ADV = 1.10
CAVITY_PACKET_FRAC = 0.50

# ── Visualisation grid ──────────────────────────────────────────────────────
GRID_RES = 160
SHELL_RES = 28
