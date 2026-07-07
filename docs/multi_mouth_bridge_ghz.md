# The GHZ sector: multipartite entanglement is bridge valence (PR #208)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. This PR closes the first named open
> of #207 — multipartite (≥3-mouth) junctions — and attaches the no-go
> that comes with it. The companion probe measures every claim (~22 s).

## 0. The open, and the object

#207 made pairwise entanglement distribution mechanism-complete and
proved its limit: bridges pair mouths into a perfect matching, and a
matching cannot produce genuinely multipartite entanglement. GHZ
requires a **junction joining three or more mouths into one bridge** —
topologically, the trousers (pair-of-pants) nucleation: three boundary
mouths sharing one bulk fiber. This PR builds that junction on the
#206/#207 lattice and derives what it makes.

## 1. The charged no-go, derived (where spin and charge part ways)

If the doublet label is a conserved **flux** (winding = charge,
#42–#44), a Y-junction must conserve it: k_A + k_B + k_C = 0.
Enumerated: the ±1 doublet has **zero** conserving triples (all sums
are ±1, ±3). Independently, GHZ-type states superpose **distinct total
charges** (our components carry totals {−1, +1}; the textbook
|+++⟩+|−−−⟩ carries {±3}) — forbidden by charge superselection, exactly
as in QM.

The pairwise sector never felt this: (k, −k) is zero-sum, so #206/#207's
Bell pairs are consistent with *both* readings of the doublet. **The
multipartite sector is where charge and spin part ways**: GHZ can only
live in the **transported-frame label** — the spin/Pin structure a
bridge carries without a flux sum (#195/#197's spinor sector). BAM thus
predicts: multipartite entanglement of charge states is impossible; of
spin states, a topology — matching the superselection structure of
quantum mechanics. A derived distinction, not an input.

## 2. The Y-junction, live

One bulk junction fiber (the tri-mouth bridge core) read by three
mouths — leg A by the identity map, legs B and C by the conjugating
transport with holonomies (s_B, s_C). A nucleation superposition
injected into the junction:

- distributes **exactly symmetrically** (three-way population spread
  10⁻¹⁶);
- imprints the per-leg **deck phases** precisely as transport theory
  demands (the conjugating legs carry η_k = e^{iπk/2} = ±i on the
  k = ±1 channels, measured);
- is read as **one shared variable** by all three mouths (readout-ratio
  spreads 10⁻¹⁶ across channels and legs);
- **leg-cut consistency**: removing leg C collapses the Y to the #206
  two-mouth bridge — mouth C receives 10⁻³³, and the extracted (A,B)
  pair phase equals the composed two-leg holonomy law
  φ = −π(s_A + s_B)/2. The valence-2 sector of #206/#207 is recovered
  as the special case of the junction.

## 3. GHZ emerges, with a holonomy law

- **The embedding generalizes.** W|k⟩ = |k⟩ ⊗ T|k⟩ ⊗ … ⊗ T|k⟩ — one
  shared junction mode read through n legs — is an **isometry** for
  n = 2 and n = 3 alike (10⁻¹⁶): the N-party tensor structure is the
  N-frame description of one bulk object.
- **The state is GHZ**, at fidelity 1.00000, with the relative phase
  obeying the **multipartite holonomy law** φ = −π(s_B + s_C)/2:

| legs (s_B, s_C) | φ measured | predicted | fidelity |
|---|---:|---:|---:|
| (2,2) | 0.0000 | 0.0000 | 1.00000 |
| (0,2) | 3.1416 | 3.1416 | 1.00000 |
| (2,1) | 1.5708 | 1.5708 | 1.00000 |

  — the #206 pair law and the #207 swapping composition, extended to
  valence 3: **one composition rule, φ = −π Σ s_legs / 2, at every
  valence**.
- **Genuinely tripartite:** 3-tangle τ = 1.0000 (the GHZ maximum;
  C²(A|BC) = 1 with pairwise concurrences 0), and the GHZ-fidelity
  witness (> ½) certifies genuine multipartite entanglement. The
  junction does not make three Bell pairs; it makes one irreducible
  triple.

## 4. Mermin = 4, from pairwise-empty marginals

- **The local bound, enumerated exactly:** all 64 deterministic local
  strategies cap the Mermin functional
  M = E(a′bc) + E(ab′c) + E(abc′) − E(a′b′c′) at **2** — the tripartite
  analog of #206's CHSH enumeration.
- **The Y-state reaches M = 4.0000** — the algebraic maximum, double
  the local bound and beyond any pairwise mechanism (two-party
  correlations cap at Tsirelson's 2√2 ≈ 2.83).
- **The irreducibility:** every two-mouth marginal of the same state is
  **unentangled** — negativities 0, pairwise CHSH exactly 2. The
  tri-mouth bridge stores all of its correlation in the triple and none
  in any pair — the exact opposite of #207's matching, where pairs
  carry everything.

## 5. The valence ledger

| bridge valence | entanglement class | measured |
|---|---|---|
| 2 (the tube) | Bell pairs; monogamy = matching; swapping | CHSH 2√2; partner ½ / non-partner 0; the a+b+c law |
| 3 (the trousers) | GHZ; pairwise-empty | Mermin 4; negativities 0; φ = −πΣs/2 |

**Bridge valence is the entanglement class**; the leg holonomies select
the state within the class (one law at every valence); charge
superselection prunes the hypergraph (flux labels admit only zero-sum
channels — the pairwise sector — while frame labels populate every
valence). W-class states (pairwise-robust multipartite) are *not*
reachable by a single junction — whether bridge networks plus surgery
reach them is a sharp successor question, stated not answered.

## 6. What is and is not established (honest scope)

**Established:** the charged no-go (enumerated); the Y-junction's
identification structure (symmetry, deck phases, one object, leg-cut
consistency with #206/#207); the GHZ state with the multipartite
holonomy law (fidelity 1.0000); 3-tangle 1; Mermin 4 against the
enumerated local bound 2; pairwise marginals exactly empty.

**Not established (the conditions):**
1. The junction fiber is the stand-in for the **trousers cobordism**;
   the full 5D pants geometry — its Sorkin causal-continuity class, its
   nucleation channel and rate relative to the #58/#200 pair channel —
   is not solved. Whether tri-mouth nucleation is dynamically
   *realized* (or only kinematically consistent) is the physical
   successor question.
2. The scalar lattice fiber models the transported-**frame** label; the
   physical GHZ carrier is the Pin spinor sector (#195/#197), per §1's
   own derived split.
3. Equilibrium hypothesis for operational statistics (#198/#204).
4. The **spatial/measurement sector** (#206 scope) remains — now the
   only standing open of the entangled-sector thread.

## References

- D. M. Greenberger, M. A. Horne, A. Zeilinger (1989); N. D. Mermin,
  PRL 65 (1990) 1838. [GHZ; the Mermin inequality.]
- V. Coffman, J. Kundu, W. K. Wootters, PRA 61 (2000) 052306. [The
  3-tangle; monogamy.]
- W. Dür, G. Vidal, J. I. Cirac, PRA 62 (2000) 062314. [GHZ vs W
  classes.]
- G. C. Wick, A. S. Wightman, E. P. Wigner, Phys. Rev. 88 (1952) 101.
  [Charge superselection.]
- The BAM chain: PR #206 (entanglement is bridge topology), #207
  (swapping is bridge surgery; monogamy is the matching), #195/#197
  (the Pin spinor sector), #200 (the cobordism machinery).

## Reproduce

```bash
python -m experiments.closure_ledger.multi_mouth_bridge_ghz_probe
# Verdict: GHZ_IS_A_THREE_MOUTH_BRIDGE_MERMIN_4_MEASURED_CHARGED_GHZ
#          _SUPERSELECTION_FORBIDDEN_PAIRWISE_MARGINALS_EMPTY_BRIDGE
#          _VALENCE_IS_THE_ENTANGLEMENT_CLASS
```
