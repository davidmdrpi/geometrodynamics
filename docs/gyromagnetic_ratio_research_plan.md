# Geometric gyromagnetic-ratio probe

Extends the spin thread (PR #60, the Wigner rotation) from *kinematics*
to the *magnetic moment*: derive the electron throat's gyromagnetic ratio
`g = 2` from the BAM geometry — the throat's Pauli/SU(2) spinor structure
(`T = iσ_y`, B2) minimally coupled to the Hopf monopole connection
(`A_φ = ½ cos χ`) — and check the leading **Schwinger anomaly**
`a = (g−2)/2 = α/2π`. A falsifier: a classical magnetic moment would give
`g = 1`; only the spinor structure gives `g = 2`.

## Why g = 2 (the Pauli/SU(2) origin)

Minimally coupling the throat spinor to the connection (`D = p − eA`) and
squaring the Dirac/Weyl operator gives the Pauli identity

```
(σ·D)² = D² + i σ·(D×D) = D² − e σ·B ,
```

using `(σ·a)(σ·b) = (a·b) I + i σ·(a×b)` and `[D_i, D_j] = −ie ε_{ijk}
B_k`. The non-relativistic Hamiltonian is then

```
H = (σ·D)² / 2m = (p−eA)²/2m − (e/2m) σ·B .
```

The magnetic-moment term is `−(e/2m) σ·B = −μ·B` with
`μ = (e/2m)·σ = (e/2m)·g·S` and `S = ½σ` — so **`g = 2`**: the Pauli term
carries the full `σ` (= `2S`), not `S`. The factor 2 is the `2` in the
`SU(2)` anticommutator `{σ_i, σ_j} = 2δ_{ij}` — the throat's
non-orientable transport `T = iσ_y = ε` (B2). g = 2 is geometric, from
the spinor algebra; a *classical* (scalar) magnetic moment would give
`g = 1`.

## g = 2 ⟺ spin tracks momentum (the BMT / Thomas link to #60)

In a magnetic field the spin precesses (Larmor) and the momentum
precesses (cyclotron); the BMT *anomalous* precession is

```
ω_a = (g/2 − 1)(eB/m) ,
```

so `g = 2 ⟺ ω_a = 0`: the spin stays locked to the momentum. This is
exactly the Thomas/Wigner result of PR #60 — the kinematic Thomas
precession (`γ²/(γ+1)`) conspires with the Larmor precession so that, at
`g = 2`, spin and momentum rotate together. The `(g/2 − 1)` is the
anomaly that `g − 2` experiments isolate.

## The Schwinger anomaly (one loop)

The leading anomalous moment is `a = (g−2)/2 = α/(2π) ≈ 0.0011614`
(Schwinger), close to the measured `a_e = 0.00115965`. This is a *one-
loop* QED vertex correction: BAM's **tree-level** geometry gives `g = 2`
exactly; the `α/2π` requires the throat vertex/self-energy loop, beyond
the tree-level geometric structure. The probe records `g = 2` as derived
and `α/2π` as the known leading correction (with the honest note that its
BAM derivation needs the loop).

## B4 accounting

`g` and the anomaly `a` are **dimensionless** (pure numbers); the
magnetic moment scale is the Bohr magneton `μ_B = eℏ/2m`, which carries
the single dimensionful anchor (`m`, via `m_e c² = ℏc/R_MID`). g = 2 is a
geometric/topological property, independent of the anchor's value —
consistent with the whole arc.

## Tests

  T1. **g = 2 from the Pauli/SU(2) algebra.** `(σ·a)(σ·b) = a·b + iσ·(a×b)`;
      `(σ·D)² = D² − eσ·B`; the σ·B term gives `g = 2` (σ = 2S);
      `{σ_i,σ_j} = 2δ_{ij}` is the factor 2.
  T2. **g = 2 from the Hopf monopole (minimal coupling).** The spin-½
      monopole (`A_φ = ½ cos χ`) minimally coupled → `μ = (e/m)S`, `g = 2`.
  T3. **BMT: g = 2 ⟺ spin tracks momentum.** `ω_a = (g/2−1)(eB/m)`;
      `g = 2 → ω_a = 0` (the #60 Thomas/Wigner link).
  T4. **Magnetic moment = Bohr magneton.** `μ = g·μ_B·s = μ_B` for
      `g = 2`, `s = ½`.
  T5. **Schwinger anomaly.** `a = (g−2)/2 = α/2π ≈ 0.0011614` (one loop);
      tree `g = 2` geometric, `α/2π` needs the loop; vs `a_e`.
  T6. **Falsification criterion.** spinor → `g = 2` (not classical
      `g = 1`); the σ-algebra forces 2; BAM passes.
  T7. **B4 accounting.** `g`, `a` dimensionless; `μ_B` carries the anchor.
  T8. **Assessment.**

## Verdict structure

  - **G_FACTOR_DERIVED** (expected): `g = 2` follows from the throat's
    Pauli/SU(2) spinor structure (`T = iσ_y`, B2) minimally coupled to
    the Hopf monopole (`A_φ = ½ cos χ`) — `(σ·D)² = D² − eσ·B`, the
    σ·B term carrying the full σ = 2S, the factor 2 being the SU(2)
    anticommutator. `g = 2 ⟺ ω_a = 0` (spin tracks momentum, the #60
    Thomas/Wigner link), and `μ = μ_B`. The Schwinger anomaly
    `a = α/2π` is the known leading loop correction (tree `g = 2`
    geometric; `α/2π` needs the loop). `g`, `a` dimensionless; scale is
    the single anchor (B4).

  - **G_FACTOR_FAILS**: the geometry gives `g ≠ 2` (e.g. the classical
    `g = 1`), or the σ·B term does not emerge — the throat is not a
    Dirac spin-½ magnetic moment.

## What this leaves open

  - **`α/2π` from BAM.** The one-loop anomaly requires the explicit
    throat vertex/self-energy loop; the tree probe gives `g = 2` only.
  - **Higher-order anomaly.** The full `a_e` series (`α²`, `α³`, …).
  - **The throat spinor from `S_BAM`.** The explicit minimally-coupled
    boosted throat spinor (shared with the #59/#60 follow-on).

## Cross-references

  - `docs/spin_wigner_rotation_research_plan.md` — the spin Wigner
    rotation / Thomas precession (#60).
  - `geometrodynamics/hopf/connection.py` — `hopf_connection`
    (`A_φ = ½ cos χ`, the spin-½ monopole).
  - `geometrodynamics/embedding/transport.py` — `T = iσ_y` (B2, the
    SU(2)/ε structure).
  - `docs/maslov_dimensional_bridge_research_plan.md` — the B4 theorem.
  - `experiments/closure_ledger/gyromagnetic_ratio_probe.py` — this probe.
