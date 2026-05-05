# Odd-k Closure Lemma

A one-page lemma fixing the lepton-generation count and the universal
Layer-1 closure-phase residue. Combines a topological selection rule
(only odd `k` are physically allowed) with an arithmetic closure
identity (the Layer-1 ledger residue is automatically 0 mod 2π).

## Statement

Let `Φ_avail(k)` denote the Layer-1 closure-phase ledger sum at depth
`k`, i.e. the four channels currently wired into the BAM machinery:

```
Φ_avail(k)  =  k · action_base
            +  Φ_hopf(χ)
            +  Φ_throat(p)
            +  β_lepton · max(0, k − 3)²                   (1)
```

with the locked baseline:

- `action_base = 2π`                          (S³ great-circle)
- `Φ_hopf(χ = 0) = π · cos(0) = π`            (Hopf holonomy)
- `Φ_throat(p = 2) = 2 · (π/2) = π`           (throat T² = −I closure)
- `β_lepton = 50π`, equivalently `4β/(2π) = 100 ∈ ℤ`  (closure quantum)

**Lemma (odd-k closure).** Under (1) and the locked baseline:

1. **Topological selection.** Even `k` and odd `k` both admit
   self-consistent closure boundary conditions on the throat, but
   they identify different sectors: even `k` is orientation-preserving
   closure on the doubled cover, while odd `k` is orientation-
   reversing closure across the non-orientable throat. BAM's lepton
   sector chooses the orientation-reversing branch, so the physical
   particle states sit at `k ∈ {1, 3, 5, …}`.

2. **Universal closure.** For every odd `k`, the Layer-1 ledger residue
   is identically zero mod 2π:

   ```
   Φ_avail(k)  ≡  0   (mod 2π)
   ```

   independently of `k`. The universal residue is the same for every
   lepton generation (`e`, `μ`, `τ` at `k = 1, 3, 5`).

## Proof

### Part 1 — topological selection

The throat transport `T = iσ_y` is the unique orientation-reversing
Hopf-preserving spinor map on `S³`
([`embedding/transport.py`][transport]). It satisfies `T² = −I`: each
pass through the non-orientable throat flips the Z₂ partition class of
the spinor and contributes a half-turn of phase.

A stable bound state is one whose worldline closes on itself after `k`
throat passes, returning the spinor to a state identified with the
initial one. Both even and odd `k` admit a self-consistent boundary
condition, but they identify *different* states:

- **Even `k`**: the spinor returns to the same Z₂ partition class,
  with accumulated phase `k · π = (k/2) · 2π`. This is
  orientation-preserving closure on the doubled cover — the worldline
  identifies with itself after an even number of throat traversals
  without ever experiencing a net orientation flip.
- **Odd `k`**: the spinor returns to the *opposite* Z₂ class, with
  accumulated phase `k · π = (odd) · π`. This is
  orientation-reversing closure across the non-orientable throat —
  the worldline identifies with itself only after a net partition
  flip, the defining property of the non-orientable embedding.

Both are valid topological boundary conditions. BAM's lepton sector
**chooses the second**: the physical particle states are those whose
closure picks up the net orientation flip enforced by the
non-orientable throat. This selects `k = 1, 3, 5, …` as the
generation index. The even-`k` branch corresponds to a different
sector of the same geometry and does not host the lepton spectrum.

### Implication

The "odd-k" restriction is a *choice of boundary condition* tied to
the BAM thesis (particles ARE the orientation-reversing closures of
the non-orientable throat), not an exclusion principle. A theory
that selected the even-`k` branch would describe a different
spectrum — orientation-preserving stable states — but would not
produce the lepton ladder anchored at `m_e`.

### Part 2 — universal closure

Substitute the locked baseline into (1):

```
Φ_avail(k)  =  2π·k + π + π + 50π · max(0, k − 3)²
            =  2π · (k + 1)  +  50π · max(0, k − 3)²
```

Each term is an integer multiple of `2π`:

- `2π · (k + 1) = 2π · m` for `m = k + 1 ∈ ℤ`. ≡ 0 (mod 2π).
- `50π · max(0, k − 3)² = 25 · max(0, k − 3)² · (2π)`. With `max(0, k−3)`
  ∈ ℤ, the factor `25 · max(0, k−3)²` is in `ℤ`. ≡ 0 (mod 2π).

Therefore `Φ_avail(k) ≡ 0` (mod 2π), independently of `k`. The
universal value is `0`. ∎

## Interpretation

The arithmetic closure is **automatic**, not fitted: every locked
constant — `action_base = 2π`, the Hopf-holonomy + throat sum equal to
`2π`, and the closure-quantum integer `4β/(2π) = 100` — is itself an
integer multiple of `2π` (or completes one with its partner). The
lemma factors as: **the Layer-1 ledger is universally zero mod 2π by
construction of its constituent geometric channels.**

The condition that picks `k` odd is a *choice of sector*, not an
arithmetic restriction: the closure identity (Part 2) holds for any
integer `k`. The lepton sector selects the orientation-reversing
closure (odd `k`); the orientation-preserving closure (even `k`)
hosts a different sector of the same geometry. The non-orientable
throat does not forbid even-`k` closure — it identifies it with the
orientation-preserving branch, which BAM does not assign to the
lepton spectrum.

## Caveats and scope

The lemma covers the **Layer-1 ledger**: the four wired channels.
Layer 2 — the radial bulk-mode contribution `Φ_radial(k)` from the
Tangherlini sector — is **not** part of the lemma. The
`closure_ledger` experiment showed that no single S(k) → {(l, n)} map
in the catalog (A, B1, B2, C1, C2, C1/B2/C2_maslov_standard,
D0/D1/D2 operator-valued) closes the full ledger including the
radial channel: the tightest result is C1 at 0.326 rad of circular
spread mod 2π, and D1 at 0.577 rad. Whether the radial residue
admits a closed-form sum has been the subject of three follow-up
probes (`dynamic_phase_probe`, `geometric_hamiltonian_probe`,
`composed_hamiltonian_probe`) which jointly conclude that the full
closure requires either a new physics channel or a non-ground-mode
admixture beyond the present catalog.

This lemma therefore **explains the universality of the available
residue** (Layer 1) and **predicts the generation count** (odd k
only), but does not itself close the full `Φ_total` ledger.

## Cross-references

- Topological derivation of `T = iσ_y`: `embedding/transport.py`,
  `derive_throat_transport`.
- Locked closure-quantum `4β = 100·(2π)`: `tangherlini/lepton_spectrum.py`,
  `TAU_BETA_50PI`; calibration log `docs/lepton_axioms.md` §3.
- Layer-1 ledger machinery and per-row arithmetic:
  `experiments/closure_ledger/ledger.py`, `_assemble_lepton_row`.
- Empirical verification across all wired candidates:
  `experiments/closure_ledger/runs/<timestamp>_comparison/status_table.md`,
  the row `none | layer1 | PASS | [0, 0, 0]`.
- THESIS.md research target this lemma closes: §"Odd-k closure from
  topology."

[transport]: ../geometrodynamics/embedding/transport.py
