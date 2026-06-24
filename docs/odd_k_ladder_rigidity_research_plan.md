# The odd-k ladder: forced, rigid, unique to the non-orientable 5D geometry (PR #174)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## One discrete feature, audited

The sensitivity audit (PR #173) measured the *continuous* predictive
content. This probe takes the cleanest **discrete** feature — the odd-k
charged-lepton ladder `k ∈ {1, 3, 5}` (e, μ, τ; even k absent) — and asks
the two questions that decide whether a discrete feature is physics or
bookkeeping:

- **(a) does the geometry force it?** — is the discreteness rigid against
  every continuous deformation of the geometry, or could it drift?
- **(b) could anything but *this* geometry produce it?** — is odd-{1,3,5} a
  signature of the non-orientable antipodal 5D geometry, or generic?

It answers them with the #173 direction analysis, split into the three sets
the inverse problem provides.

## The geometric origin (recap: PR #67 / #169 / #170)

The throat monodromy is `T = iσ_y` (`T² = −I`, the B2 spin structure).
`T^k` has period 4:

| k | T^k | sector |
|---|---|---|
| 1 | `+iσ_y` (off-diag) | odd → fermion (non-orientable, RP², Pin⁻) |
| 2 | `−I` (diag) | even → boson (orientable) |
| 3 | `−iσ_y` (off-diag) | odd → fermion |
| 4 | `+I` (diag) | even → boson |
| 5 | `+iσ_y` (off-diag) | odd → fermion |
| 6 | `−I` (diag) | even → boson |

So `k mod 2` is the orientability (= spin-statistics) grading. Charged
leptons are spin-½ **fermions** ⟹ odd k; the bulk boundary `k ≤ k_5 =
D_bulk = 5` ⟹ exactly `{1, 3, 5}` = `(k_5+1)/2` = **3 generations**.

## The three direction sets (measured)

| set | dimension | scaling exponent | reading |
|---|---:|---:|---|
| **active** | 10 | **1.03** | linear — these are the directions that move the masses and CKM |
| **null** | 10 | **2.0** | quadratic — flat to first order (~10⁴× smaller response than active at ε=10⁻²) |
| **mixed** | — | **1.02** | active-dominated — the null component hides no observable motion |

- **Active** singular directions (the rank-10 subspace of the #173
  Jacobian): perturbing along them moves the observables **linearly**.
- **Null / compensator** directions (the 10-dim kernel): **flat to first
  order** — the only motion is the second-order curvature.
- **Mixed** directions test whether nonlinear effects break the local rank
  story. **They do not**: active stays linear, null stays quadratic, and an
  (active+null) combination is dominated by its active content — the
  rank-10 picture is robust beyond the local linear approximation.

## Forced & rigid

The odd-k labels `{1,3,5}` and the generation count (3) are **invariant
under every active, null, and mixed continuous deformation** — because they
are integer winding plus the ℤ₂ orientability grading (`T² = −I`, fixed by
the B2 spin structure), discrete topological data that lives **outside** the
entire continuous deformation manifold. There is no generation-number or
k-label among the continuous inputs. The discreteness is therefore *not* an
emergent near-integer that could drift; it is structurally forced — the
right way for a discrete feature to be geometric. The continuous geometry
(rank-10 active + 10 null) deforms only the masses and the CKM.

## Unique

The odd-**only** grading is the ℤ₂ orientability of the non-orientable
antipodal quotient (`T² = −I`): a fermion's orientation-reversing closure
exists only for odd k. An **orientable** geometry (`T² = +I`) gives the
orientation-preserving even/bosonic sector — no odd-only fermion ladder.
And the specific `{1,3,5}` requires the bulk boundary `k ≤ k_5 = 5 =
D_bulk`. So odd-{1,3,5} is the **joint signature** of the non-orientable
antipodal spin structure and the 5D bulk — not generic, not tunable.

*Honest scope:* this is an exclusion/signature argument **within BAM** (the
grading requires non-orientability), not a no-go theorem against every
conceivable alternative mechanism.

## Reproduce

```bash
python -m experiments.closure_ledger.odd_k_ladder_rigidity_probe
# Verdict: ODD_K_LADDER_FORCED_RIGID_UNIQUE_TO_NON_ORIENTABLE_5D_GEOMETRY
```
