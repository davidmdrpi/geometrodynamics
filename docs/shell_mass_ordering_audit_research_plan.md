# Shell Hamiltonian mass-ordering and `n_part` audit (PR #78)

Uses the PR #77 shell waveguide basis to test whether shell
eigenvalues and Z₂ partition splitting reproduce the qualitative
quark ordering structure better than the v3 lepton-shaped basis,
and whether the `n_part = 233` compensator (PR #76) shrinks or
disappears.

## The mass-ordering target

Per the v3 basis-to-species map:

| state | species | observed mass (MeV) |
|---|---|---:|
| (k=1, +) | u | 2.16 |
| (k=1, −) | d | 4.67 |
| (k=3, +) | c | 1270.0 |
| (k=3, −) | s | 93.4 |
| (k=5, +) | t | 172690.0 |
| (k=5, −) | b | 4180.0 |

The qualitative pattern:

  - **Within-generation inversion.** The `+` partition is LIGHTER in
    generation 1 (u < d at k=1) but HEAVIER in generations 2, 3
    (c > s at k=3, t > b at k=5).
  - **Inter-generation hierarchy.** Mass² spans ~6.4·10⁹ from u to t
    — about 12 orders of magnitude.

The v3 lepton-shaped basis fits this via `uplift_asymmetry ε = 1 −
1/k_5²` (asymmetric β-uplift per partition), absorbing the unmodeled
QCD physics into `n_part = 233` (PR #76).

## What this audit tests

The PR #77 shell basis (n_varied: `(l=1, n=3,4,5, ±)`) with the
operator scaffold `H = H_kin + H_Z2 + H_couple`:

  - T1. **Shell kinetic spectrum range** — ω²(l=1, n=3,4,5).
  - T2. **Coverage gap** — shell ×2.2 vs observed ×6.4·10⁹.
  - T3. **Uniform χ is structurally insufficient** — at most 2 of 3
    blocks correct.
  - T4. **Sign-flipping χ_n is sufficient** (existence proof) — χ_3 < 0
    (+ lighter), χ_4, χ_5 > 0 (+ heavier).
  - T5–T6. **Compensator measure** — log₁₀ residual after best
    single-scale fit; kinetic-only ≈ 3.3, +χ_n ≈ 3.3 (modest
    reduction), observed log₁₀ range ≈ 9.8.
  - T7. **Comparison with v3** — shell structurally cleaner; n_part
    not reduced at PR #78 alone.
  - T8. **Honest scope** — PR #78 sharpens scope, does not close.

## Honest finding

**The shell basis is structurally better than v3** in four ways:

  1. Basis = cavity wavefronts (PR #77), not throat-traversal modes.
  2. Kinetic = `ω²(l, n)` cavity eigenfrequency squared (the natural
     "wavefront mass operator"), not phenomenological `β·k²·(2π)`
     winding cost.
  3. Z₂ partition slot is the right structural home for the
     within-generation inversion (T4 existence proof).
  4. 3 × 2 = 6 flavor count matches PR #69's structural QCD-shell
     reading.

**But three things are NOT resolved at PR #78:**

  - **Coverage gap.** Shell kinetic spans ×2.2 in mass²; observed
    spans ×6.4·10⁹. ~9 orders of magnitude unaccounted-for.
  - **Inversion mechanism.** Uniform `χ·σ_z` cannot reproduce the
    u<d, c>s pattern. Sign-flipping `χ_n` can (T4), but the values
    are illustrative — derivation is PR #79's job (boundary stress
    tensor on the cavity wall).
  - **`n_part` compensator.** PR #78 alone does NOT reduce or remove
    the v3 `n_part = 233` compensator. The shell basis identifies
    the structural slots that must populate the operator (`χ_n`
    from PR #79, `H_couple` from PR #79 + PR #80), but the values
    that span the observed range are not derivable at PR #78.

The verdict is therefore:

  - **`SHELL_BASIS_STRUCTURALLY_BETTER_N_PART_NOT_YET_RESOLVED`**.

PR #78 SHARPENS the scope of what PRs #79–#80 must do; it does not
close `n_part`.

## What PR #78 sharpens for PR #79–#80

  - **PR #79** must derive a generation-dependent `χ_n` from a
    boundary stress tensor on the cavity wall, with signs
    `χ_3 < 0, χ_4 > 0, χ_5 > 0` (or whichever sign convention
    matches the actual partition assignment) and magnitudes large
    enough to push the n=5 block up by a factor of ~10⁹ in mass²
    — or to provide the structural input that, combined with PR
    #80's `H_couple`, achieves this.
  - **PR #80** must identify the color algebra acting on `(l, n, p)`
    (SU(3) color triplet, SU(2) × Z₂ partition-flavored, Pati-Salam
    SU(4), or other) and populate `H_couple` with inter-mode mixing
    terms transforming under it. This is the channel for the
    inter-generation hierarchy.

After PR #79 + PR #80 populate the operator, the `n_part` audit can
be re-run honestly. If the spread is then structurally accounted for,
`n_part = 233` is replaced by a derived value. Otherwise it persists
as a residual phenomenological compensator with reduced scope.

## Honest scope

  - **Is:** structural comparison of the PR #77 shell basis against
    v3; explicit demonstration that uniform `χ` cannot reproduce the
    within-generation inversion; existence proof that sign-flipping
    `χ_n` can; quantitative compensator measure (log₁₀ residual);
    sharpened scope statement for PR #79 and PR #80.
  - **Is not:** a derivation of quark masses, a derivation of `χ_n`
    from physics, or a removal/reduction of the `n_part`
    compensator. The 6×6 operator's `χ` and `H_couple` slots are
    explored but not populated from a physical principle. Mass
    values shown in T4 are illustrative, not predictions.

## B4 accounting

Shell eigenfrequencies `ω(l, n)` are dimensionful (`1/length`); the
mass-ratio audit is scale-free. The absolute MeV scale rides on the
single B4 anchor `m_e = f_closure · ℏ / (ΔR·c)` (PR #53). The shell
basis does not change B4's status.

## Verdict structure

  - **`SHELL_BASIS_STRUCTURALLY_BETTER_N_PART_NOT_YET_RESOLVED`**
    (expected): T1–T7 pass; the shell basis is the right machinery,
    the structural slots are identified, but `n_part` cannot be
    closed at PR #78 alone.
  - **`AUDIT_INCONCLUSIVE`**: a structural test fails; investigate
    before proceeding to PR #79.

## What this leaves open

  - **PR #79** — `χ_n` from boundary stress tensor; singlet
    constraint.
  - **PR #80** — color algebra acting on `(l, n, p)`; `H_couple`
    inter-mode mixing.
  - **Re-audit `n_part`** after PR #79–#80 populate the operator.

## Cross-references

  - `docs/qcd_shell_waveguide_scaffold_research_plan.md` — PR #77,
    the basis + operator scaffold.
  - `docs/quark_npart_origin_research_plan.md` — PR #76, the
    diagnosis that the v3 Hamiltonian is the wrong machinery.
  - `docs/quark_axioms.md` — v3 axioms (k → species map, uplift
    asymmetry, partition mixing).
  - `geometrodynamics/qcd/quark_spectrum.py` — v3 Hamiltonian
    reference implementation.
  - `experiments/closure_ledger/qcd_shell_waveguide_scaffold_probe.py`
    — PR #77 basis + operator that this audit reuses.
  - `experiments/closure_ledger/shell_mass_ordering_audit_probe.py`
    — this probe.
