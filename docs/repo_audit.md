# Repository audit: is this program converging?

*Prompted by the observation that results are repeatedly announced, reviewed,
and reversed. This audit tests that claim against the history rather than
accepting it, and reports what it finds — including where the framing is wrong
and where the underlying problem is worse than the framing suggests.*

**Scope.** 82 merged PRs, 173 non-merge commits, 277 documents, 1405 test
functions, 324 claim rows in `README.md`, 45 recorded probe runs.

---

## The verdict in one paragraph

**It is not going in circles.** No claim in this repository has gone
`A → ¬A → A`; every reversal traced below ratchets to a strictly narrower or
better-posed statement, and several dissolve the original question rather than
answering it the other way. But that is not reassurance, because the real
finding is worse than circling: **the repository has no internal mechanism that
can say "wrong."** 45 probe runs have produced 45 passes and 0 failures. 1405
tests are green and have always been green. Every substantive correction in the
recent history was initiated from outside — by review — not by the
instruments. The rate at which bad claims are caught is therefore set by review
bandwidth, not by the science, and claim production (324 rows) outruns external
anchoring (essentially one published oracle) by more than two orders of
magnitude.

---

## Finding 1 — the reversals ratchet, they do not cycle

Three chains, traced through the commits:

**The areal-sign chain.**

| round | claim | what the next round did to it |
|--|--|--|
| #263 | TT sector gives `dA/A = 0` **identically** (algebraic; `4.5e-15`, control `3.3e13` larger) | **never touched — still stands** |
| #263 | scalar sector: `dA/A < 0`, "toward a neck", both mouths, 8 controls | #264: with the tube area *matched* to its mouths, the sign **reverses** |
| #264 | so the sign is a property of the throat | #265: the area was **never free** — `R̂ = 16πGρ` fixes `ρ` from the profile; `A = 4π` implies `ρ/3`, matched implies `133ρ`, neither is the ambient's fluid |
| #265 | the forced throat **is** an Einstein–Rosen bridge, `M = sin³(a)/2`, three quasi-local masses agreeing to `1.3e-13` | stands |

That is not a circle. It is: *contingent number → the contingency exposed → the
contingency removed → an exact result*. The endpoint is stronger than the start.

**The resonance chain.** T6 claimed the tensor channel was off resonance *by
construction* — every step true, conclusion false, because all of it described
the **uncoupled** ambient. The coupled system rings where
`det(A − Γ(ω))` vanishes, which is near `ω₃` generically. One reversal, no
return.

**The closure chain (this session).** #276 claimed no finite frequency closes
the loop. `η_topo` had been *chosen* as `+1`; `ConjugatePair` **asserts**
opposite mouth orientations, giving `−1`, and a root appears at `ω = 1.4617`.
But the same round then showed `Ψ` sweeps `3.9676 < 2π`, so **neither answer is
a property of the geometry** — a constant decides. The question was dissolved,
not flipped back.

> **So the "circles" framing is inaccurate, and I want to be precise about that
> rather than agreeable.** What is actually happening is a one-way ratchet
> running at the speed of external review.

---

## Finding 2 — the instruments cannot fail

This is the core defect.

| instrument | count | times it has ever reported a problem |
|--|--|--|
| probe runs (`experiments/*/runs/*/probe.md`) | 45 recorded | **0** |
| test suite | 1405 tests | 0 red at any commit |

A probe is written *after* the result is known, and its checks are constructed
around what the calculation already produced. `run_probe()` returns
`passed == total` by construction in every recorded run in the repository's
history. **An instrument that has never once fired is a recording device, not a
verification device**, and the `8/8 checks pass` banner at the top of each
write-up reads as evidence when it carries none.

The test suite has the same problem in a subtler form. Most of the 1405 tests
pin *machinery* — exact closed forms, unitarity, symmetry, convergence order —
and those are genuinely load-bearing. But the handful of **verdict** tests are
edited whenever the verdict changes. In `5750069` I renamed
`test_no_finite_frequency_closes_carrier_and_packet_together` to
`test_the_loop_closes_at_a_finite_frequency_on_the_declared_topology` — the test
followed the conclusion. A suite that is regenerated to match the current answer
cannot ratchet against a wrong one, and its size (1890 passing) actively
misleads about stability.

---

## Finding 3 — error detection is externalized

Of the corrections in the visible history, at least 16 commits name an external
prompt (`Address review…`, `Asked to…`, `Review corrections…`); 12 are titled
"address review" outright. In this session specifically, **all three PRs
(#274, #275, #276) had their conclusions materially changed by review**, and in
#276 the headline was reversed by review point 1.

There are genuine self-caught negatives — the TT `dA/A = 0` theorem was found
by building the tower and discovering the observable was identically zero, and
the T6 reversal was self-initiated. So the capability exists. But it is not
systematic, and it is not what the probes and tests are doing.

---

## Finding 4 — one failure mode, two flavours

Every reversal traced above has the same shape: **a quantity that was chosen,
approximated, or mis-derived entered a reported number, and nothing in the
pipeline asked whether it was entitled to.**

*Type A — an unfixed choice treated as fixed:*

| instance | the choice | consequence |
|--|--|--|
| #276 | `η_topo = +1` | verdict reversed once derived (`−1`) |
| #276 | Jost basis constant | verdict **still** only basis-relative |
| #263→#264 | tube area `A = 4π` | sign reversed under matched area |
| T6 | uncoupled vs coupled spectrum | "off resonance" → near-resonant |

*Type B — an outright derivation or estimator error:*

| instance | the error | consequence |
|--|--|--|
| #271 | radial scalar operator short by `3A²/4r²` | `γ = 22.5`, `R_OUTER` selection reopened |
| #275 | `Re[iω(T−1)]` biased by `−c³/6ω²` | the entire "0.047% deficit" attribution withdrawn |
| #275 | left-endpoint time quadrature | sum rule `−0.978`, never converging |
| #274 | diamond centre off by `h/4` | real bug — but it did **not** explain the floor |

**In every one of these, the full suite was green and the probe read 100%.**

---

## Finding 5 — what is actually stable

The repository does know things, and they share a property: they are exact,
structural, or externally anchored — never a numerical verdict on an
under-determined setup.

- `dA/A = 0` identically in the TT sector (algebraic; control `3.3e13` larger)
- support composition: `d ∉ J⁺(s)` ⟹ the transactional product vanishes for
  every intermediate event — proved, not computed
- `∫V_ℓ dr* = ℓ(ℓ+2) + 3/2` (Tangherlini) and
  `∫V_ℓ ds = (π/a)[ℓ(ℓ+2) + 9/8]` (throat), both exact
- `∫T_ab k^a k^b dλ = −3/(16G₅a)` exactly
- `M = sin³(a)/2`, three quasi-local masses agreeing to `1.3e-13`
- `V_max = 100/81` at `r² = 9/5`, the only rational `ℓ`
- `I₃ = 2/a²`, `g = π²a²`
- constraint operator exactly degenerate on `S³`: `4 = (n+1)²` at `n = 1`
- the `D = 5` fundamental QNM against Matyjasek 2021 (arXiv:2107.04815) —
  **the one genuine external oracle in the repository**
- unitarity to `6.6e-14` and structural reciprocity, imposed nowhere

Note the asymmetry: **21 "Derived" and 23 "Verified" rows against 324 total**,
and one external anchor. The exact results are also not immune — #271's operator
error was in a closed form — but they are the only category with a track record.

---

## Finding 6 — the same missing construction is behind multiple rounds

This is the most useful synthesis, and it was not visible from inside any single
round. Two independent arcs stalled on **the same** unbuilt object: the
**junction between a finite mouth and the closed `S³` exterior**, with its
surface layer.

- The areal-sign arc ended at *"the one area that is the ambient's own fluid
  (`A = 4π/3`) cannot be glued without a surface layer."*
- The closure arc ends at *"the Jost basis constant is not physically fixed;
  fixing it needs the finite-mouth matching surfaces."* `T_ℓ` is explicitly a
  whole-throat **oracle** with two asymptotically flat ends, not a glued
  finite-mouth solution.

So the provisionality of the current headline is not bad luck — it is the
**same** gap, reached from two directions. Another round of scattering
calculations on the same benchmark cannot remove it.

---

## Recommendations

Ordered by how much they would change the failure rate.

**1. Make probes capable of failing, or stop calling them verification.**
Pre-register the discriminating number *before* computing it: write the
predicted value and the falsifier into the probe, commit that, then run. A probe
whose checks are authored after the answer is a write-up with a scoreboard.

**2. Require a dependency ledger for every reported number.** Before any
headline: enumerate every constant, parameter, convention and basis it depends
on, and mark each **derived** or **chosen**. A headline that depends on a
"chosen" entry is not a result, it is a result *conditional on that entry* — and
must be stated that way. This single rule would have caught `η_topo`, the tube
area, and the uncoupled-spectrum error before publication, not after review.

**3. Split the README's claim table by evidential class.** Exact/structural,
externally anchored, and numerical-verdict claims currently sit in one table
with identical visual weight. They have wildly different survival rates.

**4. Drop the one-headline-per-round expectation.** The requirement to produce a
striking result each round is the proximate generator of over-claim; the strong
sentence is usually available, and it is usually the one that gets reversed.
"This round built the machinery and found nothing decisive" is a legitimate and
under-used outcome.

**5. Build the finite-mouth junction next, not another scattering round.**
It is the common blocker of Finding 6, and it is what would convert the current
basis-relative closure verdict into a statement about BAM.

**6. When a verdict test is edited, the commit must say what the old assertion
got wrong.** Otherwise the suite silently tracks the conclusion.

---

## What this audit does not claim

It does not assess whether the physical program — GR reproducing QED, particles
as non-orientable wormholes — is correct or promising. That is a physics
judgement outside what the commit history can settle. It assesses only whether
the *process* is accumulating reliable results, and the answer is: the
machinery and the exact results are accumulating; the numerical verdicts are
not, because nothing in the pipeline is trying to falsify them.
