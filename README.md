# Geometrodynamics

**A research framework implementing and testing Wheeler's geometrodynamic program.**

This package computationally explores the hypothesis that structures
physicists call electromagnetism, charge, spin, confinement, **black
holes**, and **Bell correlations** may emerge from the geometry of
spacetime itself — specifically the Hopf fibration on S³, 5D Tangherlini
wormholes, topological flux-tube networks, coherent wormhole-throat
condensates, and non-orientable throat topology.

## What the Code Validates

| Claim | Status | Evidence |
|-------|--------|----------|
| Charge quantisation from topology | **Verified** | c₁ = 1 to < 1e-9 error |
| Spin-½ from Hopf holonomy | **Demonstrated** | SU(2) sign-flip at 2π, illustrative |
| Coulomb law from throat eigenmode | **Verified** | BVP matches Q/r to rel_err < 3e-9 |
| α_q coupling ratios (no free parameters) | **Computed** | Forced-origin slope extraction |
| Möbius half-integer spectrum | **Verified** | Numerical vs analytic < 5% |
| Meson energy conservation | **Verified** | Drift < 1% over test window |
| Bridge nucleation / string breaking | **Verified** | Correct daughter topology |
| Hayward metric from throat density | **Derived** | n(r) → ρ(r) → m(r) → f(r) matches Hayward to < 1% |
| de Sitter EOS from Einstein eqs | **Derived** | p_r/ρ = −1 exact at all radii |
| SEC violation for regularity | **Derived** | Penrose-required SEC violation confirmed (~85% of interior) |
| Singularity avoidance (Hayward core) | **Derived** | K(0) = 24/l⁴ finite; metric now derived from throat density |
| Geodesic completeness | **Modeled** | Hayward infaller decelerates; heuristic completeness criterion |
| BH entropy from throat counting | **Consistent** | S_throat matches S_BH by construction (N set from area law) |
| Charge without charge (BH) | **Modeled** | Net Q from orientation sum, Q/N → 0 for large M |
| First law dM = T dS | **Checked** | Residual < 5%, Schwarzschild limit only |
| T from collective modes | **Derived** | T_mode matches T_surface_gravity to < 1% for M ≥ 3 |
| Core scale l ≈ Planck | **Derived** | l = 2M/√N ≈ 0.47 l_P, independent of mass |
| Schwarzschild recovery | **Verified** | Hayward → Schwarzschild as l → 0 |
| Two-horizon structure | **Verified** | Inner + outer horizons for 0 < l < l_crit |
| Singlet from throat transport | **Constructed** | T=iσ_y → |Ψ⟩ built from transport; E(a,b) = −cos(a−b) |
| T = iσ_y from Hopf fibration | **Derived** | Unique orientation-reversing Hopf-preserving map; 7 properties verified |
| Bell phases from Hopf holonomy | **Derived** | π/2 baseline + π[cos(θ_a)−cos(θ_b)]/2 from connection A = ½cos(χ)dφ |
| History closure → E = −cos(a−b) | **Derived** | SU(2) amplitudes × closure weights reproduce singlet; CHSH = 2√2 |
| History no-signaling | **Derived** | Marginals = ½ from branch enumeration; independent of remote setting |
| History conservation | **Verified** | Charge balance exact for Bell and transaction histories |
| Bulk identity Bell (kinematic) | **Verified** | Same E(a,b) from pure topology, no time stepping; both paths match |
| CHSH S = 2√2 (topological) | **Verified** | Exact Tsirelson value; topology determines correlations, cavity determines dynamics |
| No-signaling | **Verified** | Marginals = ½ from singlet; cavity dynamics don't alter spin correlations |
| Cavity detector-conditioned dynamics | **Dynamical** | Derived Hopf phases drive cavity ODE; packets fire on 0/π branches |
| Cavity persistent memory | **Verified** | Energy persists between steps; slow ring-down |
| Green kernel derivative | **Fixed** | Now matches analytic dG/dψ to < 10⁻⁴ |

### Research goals (not yet fully derived)

| Physics | Proposed geometry |
|---------|-------------------|
| Electromagnetism | Curvature of the Hopf connection on S³ |
| Particle mass | Eigenvalue of the 5D Tangherlini operator |
| QCD confinement | 1D flux-tube network with bridge nucleation |
| Retrocausal photon exchange | Wheeler–Feynman absorber theory on S³ |
| Black-hole interior | Coherent condensate of non-orientable wormhole throats |
| Bell correlations | Non-orientable throat transport + Hopf SU(2) projection |
| Entanglement = wormholes | Bell correlations from throat connectivity |
| Quantisation from resonance | S³ antipodal cavity selecting discrete spectrum |
| Topological censorship | Non-orientable throats evading standard no-go theorems |

## Package Structure

```
geometrodynamics/
├── geometrodynamics/
│   ├── constants.py          # Shared physical & simulation constants
│   ├── hopf/                 # Hopf fibration on S³
│   │   ├── connection.py     # A = ½cos(χ)dφ, curvature, holonomy
│   │   ├── chern.py          # First Chern number c₁ = 1
│   │   └── spinor.py         # SU(2) spinor transport (spin-½)
│   ├── tangherlini/          # 5D wormhole eigenmodes
│   │   ├── radial.py         # Chebyshev spectral solver
│   │   ├── maxwell.py        # Sourced Maxwell BVP (Coulomb validation)
│   │   └── alpha_q.py        # Throat flux ratios (no free parameters)
│   ├── transaction/          # Wheeler–Feynman absorber theory on S³
│   │   ├── particles.py      # ThroatMode, MouthState, Particle4, GravWave
│   │   ├── s3_geometry.py    # Geodesics, Green function, antipodal map
│   │   ├── handshake.py      # Offer/confirm/transaction protocol
│   │   └── cavity.py         # CavityMode, CavityPacket, AntipodalCavity
│   ├── embedding/            # Non-orientable throat topology
│   │   ├── topology.py       # ThroatDefect, ConjugatePair, transport ops
│   │   └── transport.py      # T = iσ_y derived from Hopf fibration
│   ├── bell/                 # Bell correlations from geometry
│   │   ├── pair_state.py     # BellPair with cavity history evolution
│   │   ├── analyzers.py      # Detector settings as SU(2) operators
│   │   ├── correlations.py   # E(a,b), CHSH, no-signaling
│   │   ├── hopf_phases.py    # Bell closure phases from Hopf holonomy
│   │   └── bulk_identity.py  # Kinematic Bell from shared bulk topology
│   ├── history/              # Closed-history framework (unifying backend)
│   │   └── closure.py        # Events, Worldlines, History, branch enumeration
│   ├── qcd/                  # Geometrodynamic QCD
│   │   ├── constants.py      # σ, α_s, ℏc, SAT parameters
│   │   ├── color.py          # SU(3) color algebra, generators
│   │   ├── bridge.py         # BridgeField, Cornell potential
│   │   ├── network.py        # Node, Branch, Junction, HadronicNetwork
│   │   ├── topology.py       # Meson, baryon, glueball, hybrid, …
│   │   ├── solver.py         # Störmer–Verlet + SAT boundaries
│   │   ├── spectrum.py       # Möbius modes, throat–branch crosswalk
│   │   └── diagnostics.py    # String tension, mode shifts, calibration
│   ├── blackhole/            # Black holes as wormhole-throat condensates
│   │   ├── condensate.py     # CoherentCondensate, ThroatState, constructors
│   │   ├── interior.py       # Hayward regular metric, geodesics, horizons
│   │   ├── entropy.py        # Bekenstein-Hawking from throat counting
│   │   └── derivation.py     # Condensate → metric via Einstein equations
│   └── viz/                  # Visualisation (placeholder)
├── tests/                    # pytest validation suite
├── notebooks/                # Jupyter notebooks (per-topic)
└── pyproject.toml
```

## Installation

```bash
git clone https://github.com/davidmdrpi/geometrodynamics.git
cd geometrodynamics
pip install -e ".[dev]"
```

## Running the Validation Suite

```bash
# All tests (fast)
pytest

# Include slow tests (bridge nucleation, string tension scans)
pytest -m ""

# Skip slow tests
pytest -m "not slow"
```

## Quick Start

### Verify charge quantisation from pure geometry

```python
from geometrodynamics.hopf import compute_c1

result = compute_c1()
print(f"|c₁| = {result['c1_abs']:.10f}  (error: {result['err_abs']:.2e})")
# |c₁| = 1.0000000000  (error: 9.99e-14)
```

### Verify spin-½ from Hopf holonomy

```python
from geometrodynamics.hopf import compute_spinor_monodromy

result = compute_spinor_monodromy()
print(f"⟨ψ₀|U(2π)|ψ₀⟩ = {result['overlap_2pi']:.6f}  (should be −1)")
print(f"⟨ψ₀|U(4π)|ψ₀⟩ = {result['overlap_4pi']:.6f}  (should be +1)")
```

### Validate Coulomb law from eigenmode throat flux

```python
from geometrodynamics.tangherlini import solve_radial_modes, solve_maxwell_from_eigenmode

modes = {}
for l in [1, 3, 5]:
    oms, fns, rg = solve_radial_modes(l)
    modes[l] = {"omega": oms, "funcs": fns}

result = solve_maxwell_from_eigenmode(modes)
print(f"Q = {result['Q']:.6f}")
print(f"Relative error vs exact Coulomb: {result['rel_err']:.2e}")
```

### Run a QCD meson simulation

```python
import numpy as np
from geometrodynamics.qcd import make_meson_tube, HadronicNetworkSolver

net = make_meson_tube(L=1.0, v=1.0, N=100, dt=0.004)
s = np.linspace(0, 1.0, 100)
net.initialize_fields(psi0={0: 0.5 * np.sin(np.pi * s)})

solver = HadronicNetworkSolver(net, antipodal_coupling=0.05)
history = solver.run(n_steps=1000, record_every=50)
print(f"Energy drift: {np.std(history['energy']) / history['energy'][0]:.4f}")
```

### Build a black-hole condensate and verify entropy

```python
from geometrodynamics.blackhole import (
    build_schwarzschild_condensate, compute_entropy_balance,
    find_horizons, integrate_radial_geodesic,
)

# Schwarzschild BH as a coherent wormhole-throat condensate
bh = build_schwarzschild_condensate(mass=5.0)
bal = compute_entropy_balance(bh)
print(f"N throats: {bh.N}")
print(f"S_BH  = {bal.S_BH:.2f}")
print(f"S_thr = {bal.S_throat:.2f}  (relative error: {bal.relative_error:.2e})")
print(f"Net charge Q = {bh.net_charge}  (neutral)")

# Nonsingular interior: Hayward metric with core scale from throat network
l = bh.core_scale
horizons = find_horizons(bh.mass, l)
print(f"\nCore scale l = {l:.4f}")
print(f"Horizons: {['%.4f' % h for h in horizons]}")

# Geodesic completeness: infalling worldline decelerates, never hits r=0
geo = integrate_radial_geodesic(M=bh.mass, l=l, r_start=3*bh.mass, tau_max=100)
print(f"Geodesic complete: {geo.is_complete}  (r_min = {geo.r_min:.2e})")
```

## Lineage

This package refactors and unifies three monolithic scripts:

| Original file | Package modules |
|---|---|
| `geometrodynamics_v39.py` | `hopf/`, `tangherlini/`, `transaction/`, `constants.py` |
| `s3_spin2_closure_toy_solver_v22.py` | `tangherlini/` (shared spectral solver) |
| `qcd_topology_solver_v30.py` | `qcd/` (entire subpackage) |
| New in v0.41.0 | `blackhole/` (condensate, interior, entropy, derivation) |
| New in v0.42.0 | `embedding/`, `bell/`, `transaction/cavity.py` |
| New in v0.43.0 | `embedding/transport.py`, `bell/hopf_phases.py`, `history/` |

## License

MIT
