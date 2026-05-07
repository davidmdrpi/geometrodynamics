# Aharonov-Bohm Hopf-fibre action probe

**Run:** 2026-05-07T22:39:37+00:00

**Connection 1-form:** `A = (1/2) cos(χ) dφ`
**Closed-form holonomy:** `∮A = π · cos(χ)`

Sub-target #2 of the ℏ-origin research plan: verify the Aharonov-Bohm Hopf-fibre holonomy numerically and document its role in the closure-cycle integer count.

## (a) Numerical holonomy vs closed form

Integrating `A_φ dφ` along a φ-loop at fixed χ on a 1000-point trapezoidal grid:

| χ | χ in π | closed form | numerical | deviation | matches? |
|---:|---:|---:|---:|---:|:---:|
| 0.0000 | 0.000π | 3.141593 | 3.141593 | 0.00e+00 | **✓** |
| 0.7854 | 0.250π | 2.221441 | 2.221441 | 0.00e+00 | **✓** |
| 1.0472 | 0.333π | 1.570796 | 1.570796 | 2.22e-16 | **✓** |
| 1.5708 | 0.500π | 0.000000 | 0.000000 | 2.47e-32 | **✓** |
| 2.0944 | 0.667π | -1.570796 | -1.570796 | 4.44e-16 | **✓** |
| 2.3562 | 0.750π | -2.221441 | -2.221441 | 0.00e+00 | **✓** |
| 3.1416 | 1.000π | -3.141593 | -3.141593 | 0.00e+00 | **✓** |

**Match:** numerical holonomy ≡ `π·cos(χ)` to machine precision at every canonical χ. The closed-form expression in `hopf/connection.py` is verified.

## (b) Spinor double-cover quantization

A spinor traversing a closed Hopf-fibre loop ONCE picks up phase π·cos(χ) (= half a closure quantum at χ=0). Under the natural fermionic double-cover (4π closure), the loop is traversed TWICE, giving accumulated phase 2π·cos(χ). The double-cover phase is an integer multiple of 2π only at the polar fibres χ ∈ {0, π/2, π}.

| χ | label | single-closure (2π units) | double-closure (2π units) | single integer? | double integer? |
|---:|---|---:|---:|:---:|:---:|
| 0.0000 | 0 (canonical north pole, χ=0) | +0.5000 | +1.0000 | — | ✓ |
| 0.7854 | π/4 | +0.3536 | +0.7071 | — | — |
| 1.0472 | π/3 | +0.2500 | +0.5000 | — | — |
| 1.5708 | π/2 (equator) | +0.0000 | +0.0000 | ✓ | ✓ |
| 2.0944 | 2π/3 | -0.2500 | -0.5000 | — | — |
| 2.3562 | 3π/4 | -0.3536 | -0.7071 | — | — |
| 3.1416 | π (south pole) | -0.5000 | -1.0000 | — | ✓ |

**Per-fibre interpretation:**

- χ = 0 (canonical north pole, χ=0): Hopf flux at maximum; single closure picks up π = 1/2 quantum, spinor double cover gives 2π = 1 full quantum.
- χ = π/4: Intermediate; cos(π/4) = √2/2 ≈ 0.707.
- χ = π/3: Intermediate; cos(π/3) = 1/2 = exact half.
- χ = π/2 (equator): Hopf flux vanishes; both single and double closure give 0 — trivial AB phase. Equatorial fibre has zero electromagnetic self-interaction.
- χ = 2π/3: Intermediate; cos(2π/3) = −1/2.
- χ = 3π/4: Intermediate; cos(3π/4) = −√2/2.
- χ = π (south pole): Antipode of canonical; single closure picks up −π = −1/2 quantum, double cover gives −2π = −1 full quantum (mod 2π = 0).

**Integer-quantization at:** χ ∈ {0.0000, 1.5708, 3.1416} — exactly the polar fibres of the Hopf base S².

## (c) Full closure-cycle integer count at χ = 0

Per-species breakdown of the Layer-1 closure cycle in units of 2π. The Hopf and throat channels each contribute 1/2 quantum at χ = 0, partnering to give a full quantum:

| species | k | antipodal (k·2π) | Hopf (π) | throat (π) | uplift (β·…) | sum (2π units) | N_layer_1 | B2_radial (n+1) | N_total |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| electron | 1 | 1 | 0.5000 | 0.5000 | 0.0000 | 2.0000 | **2** | 1 | **3** |
| muon | 3 | 3 | 0.5000 | 0.5000 | 0.0000 | 4.0000 | **4** | 2 | **6** |
| tau | 5 | 5 | 0.5000 | 0.5000 | 100.0000 | 106.0000 | **106** | 3 | **109** |

**Layer-1 integer count verified:** the sum of (antipodal + Hopf + throat + uplift) channels equals the integer N_layer_1 from the closure_cycle_action_probe. The Hopf-throat partnership at χ = 0 (each contributes 1/2 quantum, summing to one full quantum) is the angular structural piece that completes the closure-cycle integer reading.

## Verdict

**Hopf-fibre Aharonov-Bohm holonomy verified.** The closed-form `π·cos(χ)` matches numerical integration of the connection 1-form to machine precision. The spinor double-cover doubles this to `2π·cos(χ)`, integer-quantized at the three polar fibres χ ∈ {0, π/2, π}.

**Closure-cycle integer reading complete.** At the locked χ = 0:

- antipodal closure: `k · 2π` — k closures of the great circle on S³ (integer per species).
- Hopf AB: `π·cos(0) = π = 1/2 · 2π` — half quantum from the canonical-fibre holonomy.
- throat T²: `π = 1/2 · 2π` — half quantum from the spinor double-cover sign flip (T² = −I has eigenvalue arg π per closure pass).
- β-uplift: `4β·max(0, k − 3)² / (2π)` — closure-quantum integer 100 for the τ row only.

The Hopf and throat channels **partner to give a single full quantum** at χ = 0: each contributes 1/2 quantum, summing to 1·2π. This is the angular structural piece complementing the radial channel's hard-wall BS reading (closed-orbit probe) and the topological antipodal closure.

**Conceptually:** the closure cycle is built from THREE structural integer-quantum sources — antipodal closure (`k`), angular Hopf-throat partnership (`+1` at χ = 0), and the τ-uplift closure quantum (`+100` at k = 5). Plus the Layer-2 radial contribution `(n + 1)` per coupled mode. **All four pieces are integer-valued in units of 2π.** This is the cleanest formulation of the closure-cycle action quantum picture.

## What's next

With both radial (closed-orbit) and angular (Hopf AB) channels verified at integer quantization, the closure-cycle action is fully integer-quantized in units of 2π. The remaining ℏ-origin sub-targets are:

- **Sub-target #3: Tangherlini eigenvalue ↔ m_e relation.** Test whether `ℏ ω(l, n) = m_e c² · (something natural)` is consistent across species. If so, the dimensional bridge from geometric units to SI is closed.
- **Identify the physical species → (l, n) coupling.** The B2_radial_ladder gives (3, 6, 109); cross-checking against an independent observable (e.g. lepton-anomalous-magnetic-moment or some QCD identity) would pin which coupling is physical.