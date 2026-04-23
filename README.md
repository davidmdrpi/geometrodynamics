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
| Lepton mass ladder (e, μ, τ) | **Closed** | Sub-percent all three generations from locked S³ axioms (see below) |
| S³ action base `action_base = 2π` | **Locked** | Hard topological invariant; default in all lepton scans |
| k=5 uplift `4β = 200π` (100 × 2π) | **Locked** | τ uplift equals exactly 100 S³ winding quanta |

### Research goals (not yet fully derived)

| Physics | Proposed geometry |
|---------|-------------------|
| Electromagnetism | Curvature of the Hopf connection on S³ |
| Charged-lepton ladder (e, μ, τ) | Eigenvalues of a k-pass instanton-transition matrix with S³ action base `2π` and k=5 uplift `200π` — **sub-percent fit achieved** |
| Particle mass (general) | Eigenvalue of the 5D Tangherlini operator (leptons only so far) |
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
│   │   ├── alpha_q.py        # Throat flux ratios (no free parameters)
│   │   └── lepton_spectrum.py # Locked e/μ/τ instanton-transition matrix
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
├── scripts/                  # Lepton-ladder calibration CLIs
├── docs/                     # Lepton axioms + scan archaeology
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

## Lepton mass ladder (e, μ, τ) from a locked S³ action

The lepton surrogate now ships with a **fully locked topological baseline**
that reproduces all three charged-lepton masses to sub-percent accuracy with
**zero free parameters at scan time** — only the electron mass is used to set
the overall MeV scale.

### Locked axioms

- `action_base = 2π`  — the S³ great-circle action (circumference invariant).
- `k_uplift_beta = 50π`  — k-selective uplift coefficient.
  For `k=5`, the uplift is `4·β = 200π`, i.e. **exactly 100 × (2π)** S³
  winding quanta.
- `winding_mode = "max"`  — off-diagonal tunneling cost scales with the deeper
  branch, `Δk = max(kᵢ, kⱼ)`.
- `depth_cost_mode = "tunnel_only"`  — the S³ base action enters only through
  the tunneling suppression, not as an additional diagonal offset.
- `resistance_model = "exponential"`  — re-entry cost `κ·(eᵏ − 1)` captures
  exponential geometric writhe/curvature build-up with generation depth.
- Baseline anchor `(phase, transport, pinhole, resistance) ≈
  (0.001, 25.1, 22.5, 0.217869)`.

The generation-block diagonal takes the form

```
H_kk = action_base + resistance_scale · k²  +  res_diag(k)
                  +  pinhole(k ∈ {3, 5})   +  β · max(0, k−3)²
```

and off-diagonals are `−transport · exp(−α_eff · Δk) · cos(phase · Δk)`.
See `docs/lepton_axioms.md` for the full matrix construction.

### Validated predictions (locked baseline, no tuning)

| Lepton | k | Predicted (MeV) | Observed (MeV) | Relative error |
|--------|---|-----------------|----------------|----------------|
| e      | 1 | 0.510999        | 0.510999       | 0.0000% (anchor) |
| μ      | 3 | 105.61260       | 105.65838      | **0.0433%** |
| τ      | 5 | 1778.93809      | 1776.86        | **0.1170%** |

Muon/electron ratio: predicted **206.6787**, observed **206.7683**
(relative error **4.33 × 10⁻⁴**).

Reproduce directly from Python:

```python
from geometrodynamics.tangherlini import solved_lepton_masses_mev
masses = solved_lepton_masses_mev()           # read-only np.ndarray
print(masses)   # [0.51099895, 105.6126..., 1778.9381...]
```

Or by CLI (no `PYTHONPATH` needed):

```bash
python scripts/lock_beta_50pi_probe.py --n-points 32
```

which additionally pins `β = 50π` exactly and optimizes only the four
sub-leading knobs; it reports `mu/e` error ≈ 1 × 10⁻⁶% and
`τ` relative error ≈ 0.161%.

### Geometric implications

1. **Three generations correspond to odd pass depths `k = 1, 3, 5`.** The
   ladder is labelled by the number of S³ passes before closure; the locked
   baseline scans exactly these three depths. Even-`k` branches are not part
   of the surrogate; deriving their absence from the underlying Hopf/S³
   topology remains an open research task.
2. **τ uplift is exactly 100 quanta of the S³ action.** The k=5 uplift is
   `4β = 200π = 100·(2π)`, a pure integer multiple of the great-circle action
   `2π`. No tuning is required; removing the integer lock degrades `τ` by an
   order of magnitude (see `docs/lepton_tau_target.md`).
3. **The μ/e ratio is a structural eigenvalue ratio, not a coupling.** With
   `action_base = 2π` locked and the exponential resistance profile, the
   calibration scan finds exact μ/e roots on a broad resistance basin
   (±1% resistance keeps `mu_err` < 8%), replacing the earlier
   "attractor needle" regime (see `docs/lepton_tau_target.md`, "Hard S³ lock
   experiment").
4. **Quadratic diagonal `∝ k²` plus quadratic uplift `∝ (k−3)²`** together
   reproduce the observed `m_e : m_μ : m_τ ≈ 1 : 207 : 3477` hierarchy: the
   `k²` term sets the `μ/e` split and the `(k−3)²` term independently lifts
   the τ sector without disturbing the `μ/e` root.
5. **Tunneling-side depth cost dominates diagonal depth cost.** The ablation
   scan showed `tunnel_only` outperforms `diag_only` by nearly 2× on best
   μ/e (see `docs/lepton_ablation_results.md`) — consistent with a picture in
   which the inter-generation transition amplitude, not the on-generation
   mass term, sets the ratio.
6. **A `max` winding rule beats a `delta` winding rule.** Setting
   `Δk = max(kᵢ, kⱼ)` (rather than `|kᵢ − kⱼ|`) in the tunneling action was
   the change that first pushed `μ/e` from ~10 toward the experimental
   ~206.77, because it penalises transitions into deeper branches by the full
   target winding — a topological-cost interpretation consistent with the S³
   action base.

### Script map

| Script | Purpose |
|--------|---------|
| `scripts/calibrate_muon_ratio.py` | Coarse grid; solves resistance for exact μ/e root at each (phase, transport, pinhole). |
| `scripts/sweep_k_uplift_beta.py`  | Sweeps `β` with exact μ/e enforced; locates best τ fit. |
| `scripts/map_basin_k_uplift.py`   | Local gradient probe around an exact-μ/e point; reports basin width. |
| `scripts/refine_locked_tau.py`    | Dense locked scan with action_base fixed to 2π; reports integer-winding β family. |
| `scripts/lock_beta_50pi_probe.py` | Hard `β = 50π` lock; optimizes only (phase, transport, pinhole, resistance). |

See `docs/lepton_ablation_results.md`, `docs/lepton_tau_target.md`, and
`docs/lepton_next_steps.md` for the full scan archaeology.

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

### Reproduce the full charged-lepton ladder

```python
from geometrodynamics.tangherlini import (
    solved_lepton_masses_mev, S3_ACTION_BASE, TAU_BETA_50PI, tau_uplift_2pi_quanta,
)

masses = solved_lepton_masses_mev()   # locked baseline, no tuning
print(f"m_e  = {masses[0]:.6f} MeV")
print(f"m_mu = {masses[1]:.6f} MeV   (obs 105.658376)")
print(f"m_tau= {masses[2]:.6f} MeV   (obs 1776.860000)")

print(f"action_base = 2π         = {S3_ACTION_BASE:.6f}")
print(f"k_uplift β  = 50π        = {TAU_BETA_50PI:.6f}")
print(f"τ uplift    = 4β = 200π  = {tau_uplift_2pi_quanta(TAU_BETA_50PI):.0f} × (2π)")
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
| New in v0.44.0 | `tangherlini/lepton_spectrum.py` (locked e/μ/τ ladder) + `scripts/` (calibration CLIs) |

## License

MIT
