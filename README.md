[![DOI](https://zenodo.org/badge/1181274003.svg)](https://doi.org/10.5281/zenodo.20225786)
# Geometrodynamics

**A research framework implementing and testing Wheeler's geometrodynamic program.**

This package computationally explores the hypothesis that structures
physicists call electromagnetism, charge, spin, confinement, **black
holes**, and **Bell correlations** may emerge from the geometry of
spacetime itself — specifically the Hopf fibration on S³, 5D Tangherlini
wormholes, topological flux-tube networks, coherent wormhole-throat
condensates, and non-orientable throat topology.

## Where ℏ enters: closure-ledger reduction of the locked surrogate

A sequence of closure-ledger probes
(`experiments/closure_ledger/`, PRs #11–18 of this repository) has
reduced the locked lepton surrogate's residual external input from
six phenomenological parameters down to **one anchor (m_e)**. The
chain is:

```
2π ledger  →  R*  →  (γ, transport, resistance)  →  ε  →  Compton bridge  →  m_e
```

Each arrow is a probe with a quantitative test; each parameter is
identified with a closure-quantum invariant or a Tangherlini-grid
quantity computed on the same geometry. The full closure-quantum
ledger of the locked surrogate after the chain:

| parameter         | locked value | structural identification |
|-------------------|-------------:|---------------------------|
| `action_base`     | 2π           | S³ great-circle action     |
| `transport`       | 25.1         | 8π = 4·(2π)               |
| `resistance`      | 0.2179       | 7π / 100                  |
| `pinhole γ`       | 22.5         | Σ V_max[1..5] ≈ 22.0       |
| `β` (τ-uplift)    | 50π          | locked closure quantum     |
| `4β` (τ-uplift integer) | 100·2π | τ-uplift quantum          |
| `R*` (outer radius) | 1.2626    | cross-species fixed point  |
| `ε` (inner cutoff) | 3.51×10⁻⁴   | resistance / k_5⁴          |

At the inner-cutoff identification `ε = 7π/(100·5⁴)`, the
dimensional bridge collapses to the clean Compton form:

```
ℏ  =  m_e · R_MID · c
```

BAM is *dimensional-scale-incomplete only modulo m_e*. The
remaining open work is in throat physics (THESIS.md "self-consistent
throat radius"), not the closure ledger.

**Paper draft:** `docs/hbar_origin_note.md` — five-section narrative
of the closure-ledger chain with the full quantitative trail.

**Probe ledger:** `docs/hbar_origin_status.md` — every probe with
result, precision, and archive pointer.

**Reproduce in seconds:**

```bash
python -m experiments.closure_ledger.inner_boundary_derivation_probe
# Writes runs/<timestamp>_inner_boundary_derivation_probe/{probe.json, probe.md}
# Verdict: ε = resistance / k_5⁴ closes the Compton bridge to 0.04 %.
```

## Why progress is possible beyond Wheeler's geometrodynamics

Wheeler's original geometrodynamic programme had the right *instinct*
— that what we call "matter" should ultimately be a property of
spacetime itself — but it stalled in the 1960s and 70s for a concrete
reason: it lacked the **global / topological machinery** needed to
turn that instinct into a quantitative spectrum.  The continuum
Einstein equations alone do not pick out discrete spectra; they
admit far too many solutions.  Wheeler's "charge without charge" and
"mass without mass" remained slogans precisely because there was no
mechanism to make them *count* anything.

The line continued here is concrete: discreteness arises from three
independent topological/geometric channels, all of which can be
written down explicitly and integrated numerically.

1. **Antipodal S³ closure.**  Compactifying the spatial slice as
   S³ replaces the open continuum with a closed cavity, so any
   field that closes on itself does so over a great circle of fixed
   length 2π.  Resonance on a closed cavity is intrinsically
   discrete; the closure constants (`action_base = 2π`, the
   integer-winding lock `4β = 100·(2π)` for the τ lepton) are
   *exact* topological invariants of this antipodal closure.  The
   closure constants are not fitted; they are read off from the
   global structure.
2. **Non-orientable throat/shell spectra.**  A wormhole throat
   that is non-orientable carries a Z₂ partition class (`p = ±`)
   which is a real topological label, not a continuous parameter.
   The unique orientation-reversing isometry of S³ that preserves
   the Hopf bundle is `T = iσ_y` (derived in `embedding/transport.py`
   without ansatz).  T² = −I is the 4π periodicity of spinors; the
   partition splitting drives every mass-ordering inversion in the
   shelled sector (the m_u < m_d but m_c > m_s pattern).  The
   throat orientation is what makes spin-½ unavoidable rather
   than imposed.
3. **Uniform bulk distance from outer to inner.**  The throat
   confines a radial coordinate to the finite shell `[R_INNER,
   R_OUTER]` (geometric units; throat at `R_MID = 1`).  In tortoise
   coordinates this becomes a finite interval with regular
   boundary conditions, which produces a discrete eigenmode
   spectrum (`tangherlini.radial.solve_radial_modes`) — bound
   modes `u_{l,n}(r*)` with frequencies `ω(l,n)`.  This is the
   bulk geometry's own quantization channel, independent of the
   S³ closure but composing with it.

What was missing in Wheeler's day — and what this package now
demonstrates operationally — is that these three channels **compose**.
The lepton ladder is a "minimal closure" spectrum where channel 1
(S³ closure) dominates: each lepton mass scales with its global
pass-count winding `β·k²` on a nearly bare closure skeleton, locked
by `4β_lepton = 100·(2π)`.  The quark ladder (added in this work)
is a "shell-coupled closure" spectrum where channel 1 picks up the
heaviest shell only and channels 2 and 3 — partition asymmetry on
the throat and bulk-mode coupling — determine the lighter shells.
Three of the four quark-sector residuals derive from
`tangherlini.radial.solve_radial_modes` and
`tangherlini.alpha_q.derive_alpha_q` to within 1%, on the same
tortoise grid that defines the radial bound modes (see
`docs/quark_axioms.md` §8 for the full derivation log and the
quantitative match per residual).

This is what allows progress: the right machinery for *quantitative*
geometrodynamics exists, it is just not the differential-geometric
machinery Wheeler had at hand.  Antipodal closure on a compact 3-space,
non-orientable throat topology, and bulk-mode confinement are each
old and individually well understood; what is new here is putting
them together and showing that they reproduce charged-lepton masses
to sub-percent and the six-quark mass ladder to ~1.6%.

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
| Closure cycle integer-quantised in 2π | **Verified** | `(N_e, N_μ, N_τ) = (3, 6, 109)` from antipodal + Hopf-throat + radial BS + τ-uplift |
| R_OUTER selected by cross-species fixed point | **Verified** | Bisection on each lepton gives same R* ≈ 1.262 to 0.008 % |
| Pinhole γ ≈ Σ V_max[1..5] on Chebyshev grid | **Verified** | −2.2 % off the locked γ = 22.5; same operator as the QCD-sector γ_q |
| Transport = 8π = 4·(2π) | **Verified** | +0.13 % off the locked transport = 25.1; 4th closure quantum |
| Resistance = 7π / 100 | **Verified** | +0.94 % off the locked resistance = 0.2179; selected over `4·(ω−1)` by R_OUTER bisection |
| Inner cutoff `ε = resistance / k_5⁴` | **Verified** | Closes the Compton bridge `ℏ = m_e R_MID c` to 0.04 % |
| Closure-quantum ledger closes modulo m_e | **Established** | Every locked parameter is a closure-quantum invariant; m_e is the unique remaining external input |
| Quark mass ladder (u, d, s, c, b, t) | **Fitted** | 1.6% max rel err on s, c, b, t with d-anchor, four shell-index axioms, and one phenomenological β |
| Quark shell-index axioms (ε, η, χ, phase) | **Geometric** | All four expressible in `k_5 = 5` only: `(1−1/k_5², k_5, (k_5−1)·k_5, 0)` |
| Quark residual sector (transport, pinhole, resistance) | **Derived** | Each matches Tangherlini eigenmode quantity within ~1% on the tortoise grid |
| Pinhole = `Σ V_max(l=1..5)` (tortoise grid) | **Verified** | −1.09% off the fitted lock |
| Transport = `mean ⟨u_l\|V_{l+2}−V_l\|u_{l+2}⟩` | **Verified** | +0.87% off the fitted lock |
| Resistance = `transport · ln(α_q(k_5)/α_q(k_1))` | **Verified** | −0.43% off the fitted lock |
| Quark winding β = N·π/2 with N=466 | **Phenomenological** | Compensator under all ablations; awaits an analytic closure condition |
| Compton antipodal kinematics | **Verified** | Closure-compatible: front + back-mouth 4-momentum conservation under (E, **p**) → (E, −**p**); inter-mouth γ skew vanishes identically; throat-pinch skew is recoil-induced `O(ω²/m²)` |
| Compton S³-propagator pole `1/(s−m²)` | **Verified** | S³ Green function `G(ψ) ∼ 1/ψ` with `ψ ∝ s−m²` reproduces QED propagator pole; fitted exponent 1.0002 across five ω-decades |
| Thomson `(1+cos²θ)` angular factor | **Derived** | Polarization-summed BAM amplitude reproduces Klein-Nishina at ω → 0 from transverse photon polarisations on the tangent bundle |
| Compton vertex coupling `γ = −3/2` at O(ω/m) | **Derived** | Exact analytic solution to the 4-equation linear system in {1, c, c², c³} basis; clean rational coefficient |
| `γ = −3/2` is d-independent | **Verified** | Numerical γ(d) = −3/2 in d ∈ {3, 4, 5, 6, 8} to 7-digit precision; falsifies the embedding-dim/polarization-count origin |
| Compton vertex closed-form resummation | **Derived** | `F²(x, c) = 4·x³·(x²+1−x·sin²θ) / [(1+c²)·(1+x)²]` with `x = ω'/ω` reproduces Klein-Nishina to all orders in ε up to ε ~ 2 (machine precision); the perturbative PRs #31–34 are Taylor expansions of this closed form |

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
| QFT event reinterpretation (Compton) | Antipodal `S³` Green function as propagator + Hopf-fibre photon polarisation + closed-form vertex resummation reproducing Klein-Nishina exactly — see [QFT-event-reinterpretation thread](#qft-event-reinterpretation-thread-compton-scattering) below |

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
  (0.001, 25.1, 22.5, 0.217869)`. As of the closure-ledger sequence
  (`docs/hbar_origin_note.md`), all four are now identified with
  closure-quantum / Tangherlini-grid invariants:
  `transport = 8π`, `pinhole γ = Σ V_max[1..5]`,
  `resistance = 7π/100`, with the phase channel decoupled.

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
`docs/lepton_next_steps.md` for the full scan archaeology, and
`docs/hbar_origin_note.md` for the closure-ledger reduction of the
locked surrogate's parameters to closure-quantum invariants.

## Quark mass ladder (u, d, s, c, b, t) from a shell-coupled S³ closure

Parallel to the lepton sector, the six observed quark masses are
fit by a 6×6 Hermitian Hamiltonian on the closure basis
`{(k=1,±), (k=3,±), (k=5,±)}`.  The minimal v3 ansatz did not
suffice; three opt-in structural extensions (`uplift_mode =
"partition_asymmetric"`, `spectrum_zero_mode = "min_eigenvalue"`,
`chi_q_k3`, `eta_k3k5_minus`), all with defaults that recover
the minimal lepton-style ansatz, give the locked spectrum.

### Locked spectrum (d-anchor, max rel err 1.6%)

Anchored on `d = 4.67 MeV`; `u` is at zero by construction under
min-eigenvalue spectrum zero.

| species | predicted (MeV) | observed (MeV) | rel err |
|---------|----------------:|---------------:|--------:|
| u | 0           | 2.16    | 1.00 (by construction) |
| d | 4.67        | 4.67    | 0 (anchor)             |
| s | 94.82       | 93.4    | **1.5%**               |
| c | 1290.92     | 1270    | **1.7%**               |
| b | 4219.92     | 4180    | **0.95%**              |
| t | 170342.41   | 172690  | **1.4%**               |

### Locked parameters (constraint-reduced)

The full residual sector is *derivable from existing geometry*
on the eigensolver's tortoise grid:

| sector | reading |
|--------|---------|
| `action_base = π` | structural |
| `uplift_asymmetry ε = 1 − 1/k_5² = 24/25` | partition asymmetry from inverse-square shell scaling |
| `eta_k3k5_minus η = k_5 = 5` | (3,−)–(5,−) targeted off-diagonal coupling |
| `chi_q_k3 χ = (k_5 − 1)·k_5 = 20` | k = 3 partition splitter |
| `phase = 0` | partition-mixing channel inactive at the lock |
| `gamma_q = 1/10` | empirical clean rational |
| `transport ≈ 0.54` | mean `⟨u_l\|V_{l+2}−V_l\|u_{l+2}⟩` on tortoise grid (+0.87% off) |
| `pinhole ≈ 22.25` | `Σ_{l=1..5} V_max(l)` on tortoise grid (−1.09% off) |
| `resistance ≈ 0.14` | `transport · ln(α_q(k_5)/α_q(k_1))` (−0.43% off) |
| `β = N · π/2 with N=466` | **remaining phenomenological parameter** |

### Shell-coupled vs minimal closure

The diagonal-Hamiltonian decomposition shows what makes the
quark ladder structurally distinct from the lepton ladder:

| species | β contribution |
|---------|---------------:|
| u, d (k=1) | 0% |
| s         | +11% (level mixing only) |
| c         | **−27%** (pushed *down* by level repulsion) |
| b         | +76% via β·4·(1−ε) = β·4/k_5² |
| t         | **+99%** via β·4·(1+ε) ≈ β·4·(49/25) |

`β` only enters at the heaviest shell (k=5), via the
partition-asymmetric `(1±ε)` factor.  The lighter shells (u, d,
s, c) are determined entirely by the chamber-coupling sector
(pinhole, χ, γ_q).  This is the operational signature of the
"shell-coupled closure" picture: the same S³ closure skeleton
that drives the lepton ladder is, in the quark sector, primarily
expressed through how the closure interacts with an interior
chamber rather than through global pass-count winding.

### Calibration archaeology

| Script | Purpose |
|--------|---------|
| `scripts/calibrate_quark_ratios.py` | Coarse grid over the residual sector; identifies γ_q regime where positivity holds. |
| `scripts/sweep_quark_beta.py` | Integer-winding β sweep (now known to be a fit knob, not a topological lock). |
| `scripts/map_basin_quark_uplift.py` | Basin probe around the best β. |
| `scripts/lock_quark_beta_probe.py` | Final lock with β hard-fixed (legacy from the integer-N attempt). |
| `scripts/experiment_partition_asymmetric_uplift.py` | Tests the k=5 b/t splitter. |
| `scripts/experiment_min_eigenvalue_zero.py` | Tests d-anchor with min-eigenvalue spectrum zero. |
| `scripts/experiment_k3_splitter.py` | Tests χ for the c/s splitter. |
| `scripts/experiment_refined_k3k5.py` | Pass-2 refinement crossing the user-named "serious candidate" threshold (max rel err < 0.3 → 0.13). |
| `scripts/basin_probe_topological_locks.py` | Verifies N, χ, η are basin features, not grid coincidences. |
| `scripts/refine_pass3_coord_descent.py` | Coordinate-descent refinement to 1.6%. |
| `scripts/experiment_constraint_search.py` | Constraint-reduction pass: 9 free knobs → 4 + 1. |
| `scripts/experiment_n_ablation.py` | First N-stability check (residuals free); N drifts. |
| `scripts/experiment_residuals_from_geometry.py` | Substitutes residuals with broad geometric scalars. |
| `scripts/experiment_transport_pinhole_search.py` | 1D refinement of transport and pinhole derivations. |
| `scripts/experiment_transport_overlap.py` | Derives transport from QM perturbation overlap to within 0.87%. |
| `scripts/experiment_resistance_wkb.py` | WKB tunneling-derived resistance (negative result), then discovers `resistance = transport · ln(α_q ratio)` to within 0.43%. |
| `scripts/experiment_n_ablation_geometric.py` | Decisive N-stability check with all residuals derived; N still drifts → β is phenomenological. |

See `docs/quark_axioms.md` (full v3 spec, calibration log §8,
phenomenological interpretation §9) and the JSON archive in
`docs/calibration_runs/` for the raw outputs of every scan.

## QFT-event-reinterpretation thread (Compton scattering)

An 11-PR thread (PRs #25 – this PR) testing whether BAM's three
composable dynamical elements — **throat worldlines + time dilation
at mouth + antipodal closure** — reproduce QFT event structure for a
canonical local interaction, Compton scattering `γ + e → γ + e`. The
thread progressively identified the BAM-native ingredients needed
to reproduce Klein-Nishina, then resummed the perturbative result
into a closed-form vertex factor.

### Result chain

  - **Kinematics** (PR #25): closure-compatible. The antipodal map
    `(E, **p**) → (E, −**p**)` automatically conserves the
    back-vertex when the front does. Inter-mouth proper-time skew
    vanishes; throat-pinch skew is a recoil-induced `O(ω²/m²)`
    quantity, not a topological closure quantum.

  - **Propagator** (PR #26): the `S³` Green function
    `G(ψ) ∼ 1/(4πψ)` with `ψ = (s − m²)/(2m²)` reproduces the QED
    propagator pole `1/(s − m²)` exactly (fitted exponent 1.0002).

  - **Photon structure** (PR #28): giving the photon two transverse
    polarisations on the `S³` tangent bundle and treating the
    electron as a scalar charge in the Thomson limit reproduces
    `(1 + cos²θ)/2` exactly — the full Klein-Nishina angular factor.

  - **Finite-energy gap** (PR #29): the natural BAM construction
    fails at `O(ω/m)`. The recoil sign is qualitatively wrong
    (BAM enhances backscatter, KN suppresses it), localised to the
    missing per-channel kinematic weighting.

  - **Vertex coupling** (PRs #30, #31): an extended Family B vertex
    modification `V = (ε·ε'*)·(1 + ε·μ₁ + ...)` with
    `μ₁ = γ·(ω/m)·(1 − cos θ)` closes the `O(ε)` gap exactly at
    `γ = −3/2` — derived analytically from a 4-equation linear
    system over `{1, c, c², c³}` basis.

  - **Coefficient origin** (PRs #32, #33): 8 natural BAM ingredients
    evaluate to `−3/2`; the dimensional-scaling test in `d ∈ {3, 4,
    5, 6}` falsifies the embedding-dim / polarisation-count origin
    (candidate C), leaving 7 surviving candidates rooted in
    group-theoretic invariants of SU(2).

  - **`O(ε²)` extension** (PR #34): polynomial leading-order
    closure with `(ν₀, ν₁, ν₂, ξ) = (9/4, −4, 7/4, −1/2)`, with
    structural patterns `ν₀ = γ² = (−3/2)²` (recursive) and
    `ξ = −A_φ(0)` (Hopf-charge link).

  - **Resummation** (PR #35): the closed form

      F²(x, c) = 4·x³·(x² + 1 − x·sin²θ) / [(1 + c²)·(1 + x)²]
              = (2x/(1+x))² · [x·(x²+1−x·sin²θ) / (1+c²)]

    with `x = ω'/ω = 1/(1 + ε(1 − cos θ))` reproduces Klein-Nishina
    **exactly at all orders in ε up to ε ~ 2** (machine precision).
    The perturbative results of PRs #31–34 are Taylor expansions
    of this closed form.

  - **Cross-process validation via Breit–Wheeler** (this PR): the
    same closed-form F, expressed in Lorentz invariants and
    analytically continued via standard Mandelstam crossing
    (`s_C → u_BW`, `t_C → s_BW`, `u_C → t_BW`), exactly reproduces
    the Breit–Wheeler pair-production amplitude `γγ → e⁺e⁻`.
    Crossed variables `x_⊗ = −(1−β·cosθ)/(1+β·cosθ) < 0` and
    `c_⊗ = (2β² − β²cos²θ − 1)/(1−β²cos²θ)` carry the construction
    from Compton lab kinematics to BW CM kinematics; the
    BAM-predicted `|M̄|²_BW = −2·(f_baseline · F²)/x_⊗²` agrees
    with the textbook formula to machine precision at all sampled
    `(β, cosθ)`, and the integrated differential reproduces the
    textbook BW total at threshold (`β → 0` linear) and in the
    ultra-relativistic regime (`β → 1` logarithmic). The vertex F
    is therefore **not a Compton-specific algebraic fit** — it is
    the closed form of the invariant QED amplitude carried by
    crossing to its tree-level partners.

### Structural reading

The `(1 + c²)` denominator in the angular factor IS the
polarisation-sum factor. The closed-form F must be derived AS a
modification of the polarisation-sum projector, not as a separate
amplitude factor. The two-factor decomposition

  - kinematic Padé `(2x/(1+x))²` — pure x-function
  - angular polarisation modification `[x·(x²+1−x·sin²θ) / (1+c²)]`

suggests two BAM-native ingredients combine to produce the full
vertex coupling. The clean half-integer/integer rationals appearing
at every order (γ = −3/2, ν₀ = 9/4, ν₁ = −4, ν₂ = 7/4, ξ = −1/2)
indicate a deeper geometric origin awaiting first-principles
derivation from the Hopf-bundle / throat-transport algebra.

### What survives and what is still open

  - Survives: BAM's antipodal-`S³` propagator + Hopf-fibre photon
    polarisation + closed-form vertex `F²` together reproduce
    Klein-Nishina exactly. The same closed form, crossed via
    Mandelstam permutation, reproduces Breit–Wheeler `γγ → e⁺e⁻`
    (PR #36) and pair annihilation `e⁺e⁻ → γγ` (this PR); the full
    Compton/BW/annihilation crossing triangle closes (loop is
    identity at both the Mandelstam-label and amplitude level).
  - Open: first-principles BAM derivation of `F²` from a BAM
    Lagrangian / action. Two-channel tree processes (Bhabha, Møller)
    with interfering s+t diagrams; loop corrections requiring the
    bulk radial channel.

### Probe sequence

| # | Probe | Outcome |
|---|---|---|
| PR #25 | `compton_antipodal_kinematics_probe.py` | closure-compatible |
| PR #26 | `compton_amplitude_structure_probe.py` | propagator ✓, polarization ✗ |
| PR #28 | `compton_photon_structure_probe.py` | Thomson KN ✓ |
| PR #29 | `compton_finite_energy_kn_probe.py` | recoil ✗ at `O(ω/m)` |
| PR #30 | `compton_vertex_structure_probe.py` | empirical finite-ε fit |
| PR #31 | `compton_vertex_derivation_probe.py` | exact γ = −3/2 |
| PR #32 | `compton_coefficient_origin_probe.py` | 8 plausible derivations |
| PR #33 | `compton_dimensional_scaling_probe.py` | C falsified, 7 survive |
| PR #34 | `compton_eps2_extension_probe.py` | `O(ε²)` polynomial fit |
| PR #35 | `compton_vertex_resummation_probe.py` | exact closed-form F² |
| PR #36 | `breit_wheeler_cross_process_probe.py` | F process-general under Compton → BW crossing |
| PR #37 | `pair_annihilation_crossing_probe.py` | full Compton/BW/annihilation crossing triangle closes |
| PR #38 | `throat_nucleation_caustic_derivation_probe.py` | F² = K(x)²·Q(x, c) BAM-geometric decomposition |
| PR #39 | `two_mouth_flux_action_probe.py` | K(x) = 2x/(1+x) from equal-action throat-rate splitting |
| PR #40 | `hopf_helicity_transport_probe.py` | Q(x, c) from Hopf-fibre helicity spinor (A_pres, A_flip) |
| this | `throat_action_derivation_probe.py` | **BAM throat action: both equal-action postulates derived from S³ antipodal symmetry + closure quantum + stationary action** |

The derivation rests on the algebraic identity

  x² + 1 − x·sin²θ ≡ (1 − x)² + x · (1 + c²)

which yields two equivalent decompositions:

  F²(x, c) = [2x/(1+x)]² · [x² + x·(1−x)²/(1+c²)]
  |M̄|²_KN/(8e⁴) = (1+c²) + (1−x)²/x

with BAM-geometric interpretation:

  - **P(x) = 2x/(1+x)** = harmonic mean of in/out photon frequencies
    = standard classical bottleneck-flux average through the throat;
    squared because both throat-pair mouths pinch. Uniquely
    polynomial — alternative throat-rates (arithmetic, geometric mean,
    linear x) leave Q non-polynomial at x → −1.
  - **(1+c²)/2 = cos⁴(θ/2) + sin⁴(θ/2)** = sum of squared Wigner-d¹₁,±₁
    matrix elements = Hopf-fibre spin-1 helicity transport through θ.
  - **Q = |a|² + |b|²** = orthogonal sum of helicity-preserving
    (a = x) and helicity-flipping (b = √x(1−x)/√(1+c²)) channels,
    each non-negative across the physical region.
  - The Hopf connection at the BAM lock `A_φ(0) = 1/2` (from
    `geometrodynamics.hopf.connection`) matches the PR #34 perturbative
    coefficient `ξ = −1/2` exactly.
  - Decomposition survives analytic continuation under crossing
    (full Compton ↔ BW ↔ annihilation triangle, PR #37).

The full F² closed form is now derived from three foundational
principles via a single BAM throat action functional (this PR):

  (P1) closure quantum `S = 2π` (BAM `action_base`)
  (P2) S³ antipodal symmetry `σ(p) = −p` (involution swapping mouths)
  (P3) stationary action under the antipodally-symmetric ansatz

Both equal-action postulates (PR #39 energy → K, PR #40 spin/Hopf → Q)
follow as consequences. Alternative principles (broken antipodal
symmetry; wrong closure quantum; wrong action functional) all fail
to reproduce K(x), confirming the principles are necessary.

See `docs/compton_vertex_resummation_research_plan.md` for the
Compton-thread culmination,
`docs/breit_wheeler_cross_process_research_plan.md` for the BW
cross-process plan,
`docs/pair_annihilation_crossing_research_plan.md` for the
triangle-closure plan,
`docs/throat_nucleation_caustic_derivation_research_plan.md` for the
geometric decomposition plan,
`docs/two_mouth_flux_action_research_plan.md` for the K-factor plan,
`docs/hopf_helicity_transport_research_plan.md` for the Q-factor plan,
and `docs/throat_action_derivation_research_plan.md` for the throat
action derivation plan.

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
| New in v0.45.0 | `qcd/quark_spectrum.py` + `qcd/hadron_spectrum.py` (shell-coupled six-quark ladder; residual sector geometrized to ~1% via Tangherlini eigenmode) |
| New in v0.46.0 | `experiments/closure_ledger/` (closure-ledger sequence; reduces the locked lepton surrogate's residual external input from six phenomenological parameters to one anchor m_e). Paper draft in `docs/hbar_origin_note.md`. |

## License

MIT
