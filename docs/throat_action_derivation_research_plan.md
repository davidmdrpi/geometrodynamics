# Equal-action derivation from BAM throat action — research plan

Closes the deepest open question from PRs #38/#39/#40: the equal-action
splitting at the two throat mouths was *postulated* (flux continuity)
in PR #39 (for energy → K factor) and PR #40 (for spin → Q channels).
This probe derives **both** postulates from a single BAM throat action
functional via stationary action under S³ antipodal symmetry.

## The BAM throat action

For a photon traversing a closed orbit through a throat connecting
two antipodal mouths on S³, the action functional is

```
S[γ; ω(s), A(s)]  =  ∮_γ [ω(s) · dt/ds + A_Hopf(s) · dχ/ds] ds
```

where:

  - `γ` is the closed orbit (two segments, one per mouth)
  - `ω(s)` is the local photon angular frequency on the orbit
  - `A_Hopf(s) = ½·cos(χ(s))` is the Hopf connection
  - `dt/ds` is the proper-time rate along the orbit
  - `dχ/ds` is the Hopf-fibre rotation rate

Parameterising by per-segment time:

```
S = S_energy + S_Hopf
S_energy  =  ω_1·τ_1 + ω_2·τ_2
S_Hopf    =  A_φ · Δχ_1 + A_φ · Δχ_2  =  (1/2)·(Δχ_1 + Δχ_2)   at χ=0 lock
```

## The three principles

### P1. Closure quantum

The closed orbit must accumulate exactly one BAM closure quantum:

```
S_energy = 2π·ℏ     (one S³ great-circle action, action_base from repo)
S_Hopf   = 2π·ℏ     (one full Hopf-fibre revolution)
```

### P2. S³ antipodal symmetry

The two mouths are at antipodal points: `p_2 = antipode4(p_1) = −p_1`.
The antipodal map `σ: p → −p` is an involution (`σ² = id`) that
exchanges mouth 1 and mouth 2. The throat action `S` is invariant
under `σ`:

```
S[σ(γ)] = S[γ]   for all closed orbits γ
```

### P3. Stationary action under symmetric ansatz

By P2, the extremum of `S` lies in the antipodally-symmetric
configuration space. Stationary action subject to the closure
constraint (P1) and symmetry (P2) gives:

```
∂S/∂τ_i = λ·ω_i  (Lagrange multiplier λ from closure)
ω_1·τ_1 = ω_2·τ_2  (antipodal symmetry: equipartition of action)
```

Combining with `ω_1·τ_1 + ω_2·τ_2 = 2π`:

```
ω_i·τ_i = π     (equal-action splitting, FORCED by symmetry)
```

## Consequences

### For energy (reproduces PR #39)

  - `τ_i = π/ω_i`
  - Total period: `T = π·(ω_1 + ω_2)/(ω_1·ω_2)`
  - Effective angular frequency: `ω_eff = 2π/T = 2ω_1ω_2/(ω_1+ω_2)`
    (harmonic mean of mouth frequencies)
  - `K(x) = ω_eff/ω_1 = 2x/(1+x)` (Padé factor)

### For Hopf rotation (reproduces PR #40)

Same symmetry argument applied to `S_Hopf`:

  - Per-mouth Hopf rotation `Δχ_i = π` (half-revolution)
  - Per-mouth helicity transport amplitude follows Wigner-d¹:
    `cos²(θ_i/2)` for preserving, `sin²(θ_i/2)` for flipping
  - Recoil-deficit coupling: when ω_1 ≠ ω_2 (Compton recoil), the
    helicity-flip amplitude picks up a factor (1−x) (energy
    transferred to electron drives angular momentum kick)
  - Hopf normalisation: the (1+c²)/2 = Wigner-d² sum normalises
    the flip channel amplitude

Result: `A_pres = x`, `A_flip = √x·(1−x)/√(1+c²)` as in PR #40.

## Tests

  T1. **Closure quantum**: `action_base = 2π` from repo (recap PR #39).

  T2. **S³ antipodal involution**: verify `antipode4(antipode4(p)) = p`
      and `antipode4(p) = −p` for sample S³ points.

  T3. **Antipodal-symmetric action extremum (energy)**: numerically
      extremise `S = ω_1·τ_1 + ω_2·τ_2` subject to closure and
      symmetry; verify `ω_i·τ_i = π`.

  T4. **K(x) from extremum**: derive `K = 2x/(1+x)` from the
      stationary-action solution.

  T5. **Antipodal-symmetric action extremum (Hopf)**: same principle
      applied to `S_Hopf = A_φ·(Δχ_1 + Δχ_2)` with `Δχ_1 + Δχ_2 = 2π`;
      derive `Δχ_i = π`.

  T6. **A_pres = x derived**: per-mouth amplitude √x from PR #39 +
      two-mouth product = x. (Recap PR #40 T3.)

  T7. **A_flip = √x(1−x)/√(1+c²) derived from coupled energy–Hopf
      action**: the recoil deficit (1−x) provides the angular
      momentum kick; the (1+c²)/2 Hopf sum normalises the flip
      amplitude. Per-mouth derivation chain.

  T8. **Alternative principles rejected**: test variations that
      break (a) antipodal symmetry, (b) closure quantum value,
      (c) stationary-action condition. Verify none reproduce
      K and Q simultaneously.

  T9. **End-to-end F² reconstruction**: combine derived K and Q;
      verify `F²(x, c) = K²·Q` to machine precision (PR #38 +
      PR #39 + PR #40 + this).

## Verdict structure

  - **ACTIONS_DERIVED**: T1–T9 all pass. Both equal-action postulates
    (energy splitting → K, Hopf splitting → Q) follow from a single
    BAM throat action functional via stationary action under S³
    antipodal symmetry.

  - **PARTIAL_DERIVATION**: energy splitting derives cleanly but
    the Hopf-rotation splitting (or the recoil-coupling for A_flip)
    requires additional ansatz beyond the action functional.

  - **DERIVATION_INCOMPLETE**: alternative principles also produce
    K and Q; the proposed derivation is not unique.

## What this leaves open

  - **The Hopf-connection coupling form**: the BAM throat action uses
    `A_Hopf = ½·cos(χ)` from `geometrodynamics.hopf.connection`. The
    *coefficient* ½ (= `A_φ(0)`) is what links the Hopf rotation
    closure to the PR #34 `ξ = −½` coefficient, but its first-principles
    origin (why ½ and not another value?) is open.

  - **Loop corrections**: still tree-level.

  - **Higher closure modes (n > 1)**: this probe targets the lowest
    closure mode `S = 2π`; whether higher modes give consistent
    spectra is a separate (longer) target.

## Cross-references

  - PR #38: `throat_nucleation_caustic_derivation_probe.py` — F² = K²·Q.
  - PR #39: `two_mouth_flux_action_probe.py` — K from equal-energy-action.
  - PR #40: `hopf_helicity_transport_probe.py` — Q from helicity spinor.
  - `geometrodynamics/transaction/s3_geometry.py` — `antipode4`.
  - `geometrodynamics/hopf/connection.py` — `A_φ(χ) = ½cos(χ)`.
  - `experiments/closure_ledger/ledger.py` — `action_base = 2π`.
  - `experiments/closure_ledger/throat_action_derivation_probe.py`
    — this probe.
