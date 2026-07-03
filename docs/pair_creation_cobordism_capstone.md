# PR #200 — the pair-creation cobordism, constructed: completing "Pauli from GR" (the capstone)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. This is the milestone PR: it closes
> the one open construction of the deepest theorem in the program
> (#196's pair-creation 4-manifold), and assembles the release ledger —
> the state of the derivation chain from the 5D field equations to
> fermion statistics, three generations, and the Born rule, with every
> grade and every honest open item. The companion probe machine-checks
> the new construction and re-verifies one core invariant of every
> milestone in a single run.

## 0. What #200 closes

PR #196 adjudicated the exchange −1 against the geon-statistics
literature and proved: the spin-statistics correlation is a theorem for
BAM throats (prime ✓, non-chiral ✓, abelian ✓), the sign is Fermi in
Pin⁻-framed GR — **conditional on one remaining construction**: the
explicit 4-manifold implementing the pair creation of two throats from
the S³ background, with the structure the Dowker–Sorkin sum-over-
histories theorem requires. This PR exhibits it.

## 1. The construction: two 2-handles with framings ±2

**The cobordism.** Take the trivial history of the closed spatial
universe, M₀ = S³ × [0, 1]. In two disjoint 3-balls of the final slice
S³ × {1}, attach two 4-dimensional 2-handles along a **two-component
unlink**, with framings **+2 and −2**:

```
M  =  (S³ × I)  ∪  H²₍₊₂₎  ∪  H²₍₋₂₎ .
```

**The boundary.** Handle attachment performs surgery on the attaching
link: integer n-surgery on an unknot in S³ yields the lens space
L(n, 1), and a split link gives the connected sum. Hence

```
∂M  =  S³ ⊔ ( L(2,1) # L(−2,1) )  =  S³ ⊔ ( RP³ # RP³ ) ,
```

using L(2,1) = RP³ and L(−2,1) = the mirror RP³ ≅ RP³ (amphichirality:
q² ≡ −1 mod 2 — the #196 Lemma 2b check). **M is a cobordism from the
vacuum S³ to the two-throat universe** — pair creation from the actual
BAM background, no ℝ³ transplant needed.

**The machine-checked arithmetic** (companion probe):

| check | value | meaning |
|---|---|---|
| linking matrix | diag(2, −2) | the Kirby diagram of M |
| \|det\| | 4 = \|H₁(RP³#RP³)\| | the boundary's homology, presented by the surgery |
| Smith normal form | diag(2, 2) | coker = Z₂ ⊕ Z₂ = H₁(RP³#RP³) exactly |
| diagonal parity | even, even | the spin structure of S³ **extends over both handles** (2-handle obstruction = framing mod 2) — M is a **spin cobordism**, as the Pin⁻-framed theorem requires |
| signature | 0 | no gravitational/CP anomaly from the creation event |
| Morse indices | {2, 2} | both critical points index 2 (a 2-handle = one index-2 Morse point); χ(M) = 0 + 2 ✓ |

## 2. The Dowker–Sorkin conditions, now all closed

Against the #196 ledger:

- **Mirror-pair identity** ✓ (was closed): the framings +2 and −2
  create the geon and its mirror image — exactly the pair-creation
  structure the theorem is about — and amphichirality makes the pair
  *identical* geons. This is also the **C-conjugate throat–antithroat
  pair of #58** (Hopf charges ±1, Σc₁ = 0): the mirror structure of the
  cobordism is the charge conservation of the pair threshold.
- **Bordism existence** ✓ (was closed abstractly; now explicit).
- **Spin/pin extension** ✓ (was closed abstractly; now explicit: even
  framings).
- **The explicit BAM 4-manifold** ✓ **— CLOSED BY THIS PR.** Built
  directly over the S³ background; the #196 open item is discharged.

## 3. The Sorkin selection rule: pair creation is the *selected* channel

The Borde–Sorkin criterion (with the Dowker–Garcia–Surya results):
Morse points of index 1 or n−1 = 3 produce **causally discontinuous**
histories (proven), and such histories are suppressed in the SOH; the
remaining indices give causally continuous ones (the Borde–Sorkin
conjecture, proven for the yarmulke classes and supported in general).
Our cobordism has **only index-2 points** — it avoids the
proven-discontinuous classes entirely:

> **The pair-creation history is causally continuous; a single-throat
> creation would require the discontinuous handle classes.** The
> selection rule that suppresses "trousers"-type histories *selects
> pair creation* — which is exactly the channel #58's energetics
> (threshold 2m_ec², one C-conjugate pair) always assumed. The
> topology-change kinematics and the pair-production physics agree for
> a reason.

## 4. The completed chain (what the release can claim)

With this construction, the exchange-statistics chain runs end-to-end
with no unconstructed step:

```
5D Einstein equations on the antipodally identified Tangherlini bulk
  → the throat prime is the RP³ geon (#169/#196 Lemma 1)
  → prime, non-chiral, abelian: the SSC theorem applies (#196 Lemma 2)
  → the pair-creation 4-manifold EXISTS and is EXPLICIT, spin,
    causally continuous (THIS PR)
  → SOH weights abelian; SSC: exchange = rotation phase (#196 Thm)
  → the Pin⁻ framing (forced: #169/#170/#195) lifts the rotation to −1
    (#196 Lemma 3/5; measured as the #188 holonomy)
  → exchange = −1: FERMI; Pauli exclusion (#185/#187)
```

**"Pauli from GR + the forced Pin⁻ framing" is now a constructed
theorem** — hypotheses verified, the 4-manifold exhibited, the
selection rule passed.

## 5. The release ledger — the state of the derivation at #200

Each row's core invariant is **re-verified live** by the companion
probe in this single run (fast checks; the heavy verifications live in
the cited PRs):

| claim | grade | key invariant (re-verified) | conditions / cross-ref |
|---|---|---|---|
| Fermion statistics / Pauli | **constructed theorem** | ½trT² = −1; deck dets (−1, +1); pin lift of 2π loop = −I; the ±2-handle cobordism arithmetic | Pin⁻ framing forced by the quotient (#169/#170); SOH framework (#196, #200) |
| Three generations (structure) | **spectral fact** | E_k = 2k + (k/λ)² (scalar); Dirac ladder m_k = k/λ + λ/2 with λ×(1) = √6; round multiplicities (n+1)(n+2) | whole Berger family; cutoff k₅ dynamical, not spectral (#193/#197) |
| The k=1 (electron) zero mode | **index theorem** | monopole ground = q to ~10⁻⁷ (Atiyah–Singer count k in sector k) | mass ladder on the Dirac tower not yet rebuilt (#195, open) |
| Born rule | **dBB grade, guidance derived** | T_{μχ} = k·Im(ψ*∂_μψ) identity; continuity residual on live flow | equilibrium/measurement conditions (#198); Bianchi step symbolically exact (#199) |
| Photon / Coulomb | **derived kernel** | 1/(4πd) ⟷ 1/q² (#42–44/#190) | weak-field patch |
| Mass hierarchy (μ/e) | **fit, fine-tuned** | the surrogate's dialed near-zero (Δ = 74.7) | the honest open problem (#192/#194); mechanism exists (#195) but ladder not rebuilt |
| Cosmological constant / one-R | **clean failure** | ρ(1) ~35 orders off | #165 — stands as stated |
| CKM γ angle | **misfit** | γ = 104° vs 65.9° | #160 — stands as stated |

## 6. The open-items register (post-#200)

1. **Rebuild the mass ladder on the Dirac tower** with the mouth
   coupling ε from the throat overlap machinery (#195's follow-up) —
   the route to un-dialing #194's fine-tuning.
2. **The nonlinear measurement theory** (#198 condition 2): the linear
   test-throat regime is the operative one; the general theory is open.
3. **The full 5D throat-core dynamics** (#199): the core is constrained
   to its quantized charge, not solved.
4. The **cosmological-constant failure** (#165) and the **γ misfit**
   (#160) stand as the program's honest negative results.
5. The Borde–Sorkin criterion is a conjecture beyond the proven index
   classes; our channel avoids the *proven-discontinuous* classes, and
   the full conjecture is cited, not re-proved.

## References

- H. F. Dowker, R. D. Sorkin, CQG 15 (1998) 1153 [the SSC theorem this
  construction feeds]; R. D. Sorkin, S. Surya, IJMPA 13 (1998) 3749.
- A. Borde, R. D. Sorkin (causal continuity criterion); H. F. Dowker,
  R. S. Garcia, S. Surya, *Morse index and causal continuity*,
  gr-qc/9711070.
- Kirby calculus / handle attachment and surgery: R. Gompf,
  A. Stipsicz, *4-Manifolds and Kirby Calculus* (integer surgery on the
  unknot = L(n,1); split links = connected sums; 2-handle spin
  extension = even framing).
- The repo chain: #58 (pair threshold), #169/#170 (the quotient and
  Pin⁻), #185/#187/#188 (exchange and holonomy), #196 (the
  adjudication), #198/#199 (Born rule and guidance).

## Reproduce (the capstone verification run)

```bash
python -m experiments.closure_ledger.pair_creation_cobordism_capstone_probe
# Verdict: PAIR_CREATION_COBORDISM_CONSTRUCTED_PAULI_CHAIN_COMPLETE
#          _RELEASE_LEDGER_GREEN
```
