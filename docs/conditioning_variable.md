# Round 8 — is the `|D|` density forced by the closure set?

Pre-registered in `docs/conditioning_variable_prereg.md` (`39be3e3`) before
`geometrodynamics/bulk/conditioning_variable.py` existed. Module, tests
(`tests/test_conditioning_variable.py`, 11), probe
(`experiments/closure_ledger/conditioning_variable_probe.py`, 10/10).

**Verdict: `CONDITIONING_VARIABLE_IS_A_CHOSEN_INPUT_JUSTIFIED_BY_THE_PHASE_AXIOM`.**

A correction round. It changes the recorded *status* of the input every result
of rounds 5–7 rests on, and by pre-registered rule it changes **no number**
those rounds produced (worst movement `3.2e-11`, Monte-Carlo noise).

## What was wrong

The repository said, in four places, that

> Haar conditioned on `N = 0` is the coarea measure, density `|D|/(2|u×v|)`

and ledgered it **`[derived; window limit]`**. But `window_monte_carlo`
conditions on `|Ω mod 2π| < ε` — the **phase** — and the two prescriptions have
the same support and different limits. The claim as written is false, and the
ledger entry attributed to derivation something that is a choice.

Found by an external audit of `92a915b` dated 5 September 2026.

## Why the two windows disagree

Coarea density is `1/|∇(conditioning variable)|`. On the closure circle `Γ`,
where `x ⊥ q`:

| variable | gradient on `Γ` | coarea density |
|---|---|---|
| `N = x·(u×w)` | `∇N = P_x(q) = q`, so `\|∇N\| = \|a×b\| = sin γ`, **constant** | uniform in arclength |
| `θ = atan2(N,D)` | `∇θ = q/D`, so `\|∇θ\| = \|q\|/\|D\|` | proportional to **`\|D\|`** |

So the `|D|` density is the coarea measure *with respect to the phase*. It is
not the coarea measure of the closure set, because a set does not have one.
Conditioning a measure on a measure-zero set is not determined by the set —
Borel–Kolmogorov — and the limiting family is part of the specification.

There is a second, sharper reason the `N` window cannot be the intended one:
with `u = s_A a` and `w = −s_B b`,

```
u × w = −s_A s_B (a × b),
```

so **`|N|` is identical in all four outcome sectors**. An `|N| < ε` window
selects literally the same set of `x` in every sector. Measured at `γ = 1`,
`ε = 0.01`, `n = 5×10⁵`: all four counts `5956` — equal, not approximately
equal — giving `E = 0` exactly, at every angle and every `ε`.

| `ε` | `E`, `|N| < ε` | `E`, phase window |
|---|---:|---:|
| `0.05` | `0.000000` | `0.397021` |
| `0.02` | `0.000000` | `0.400039` |
| `0.01` | `0.000000` | `0.398457` |
| `0.005` | `0.000000` | `0.394911` |
| `0.002` | `0.000000` | `0.398601` |

against the closed form `0.3984966504`. Same support, two limits.

## Why the repository's choice is nevertheless the right one

Pre-registered rule 1 forbade preferring the phase window because it gives the
interesting answer; a justification had to be a structure already present,
named by file and line. There is one, and the module reads it out of the source
rather than quoting it:

```
geometrodynamics/history/closure.py:11
  3. Phase closure: total phase around every loop ≡ 0 or π (mod 2π)
```

The repository's third closure condition is stated **on phase**. A phase
tolerance is the natural regularisation of that axiom; `N = 0` is a derived
description of the same locus, not the axiom itself. So the choice is
justified — and justified is not derived. The ledger now carries it as an
input:

```
conditioning variable = phase, not N   [chosen; justified by the axiom
                                        being stated on phase]
coarea density given that variable     [derived: 1/|∇θ|]
```

This matters because `|D|` is what feeds round 5's positive branch, round 6's
oriented current, and round 7's Morse–Bott component masses. The single most
load-bearing input in the chain was recorded as derived; it is chosen.

## What did not change

Pre-registered rule 3: no round 5–7 number may move, since this corrects a
status and not a computation.

| quantity | expected | got | Δ |
|---|---:|---:|---:|
| round-5 `correlation(abs)` at `γ=1` | `+0.3984966504` | `+0.3984966504` | `2.2e-11` |
| round-5 `correlation(signed)` | `+0.5403023059` | `+0.5403023059` | `3.2e-11` |
| round-6 oriented singlet `E` | `−0.5403023059` | `−0.5403023059` | `0.0` |
| round-7 `M_0` | `+11.6708058364` | `+11.6708058364` | `0.0` |
| round-7 `M_π` | `+0.1695122784` | `+0.1695122784` | `0.0` |

## Two further narrowings from the same audit

**`κ` is route-local.** Round 7's prose said the round "adds a fourth (`κ`)",
which reads as a fourth *universal* underived input beside branch aggregation,
sector coefficients and readout. `κ` is the normalisation of `e^{iκS_H}` and
exists only inside the holonomy-trace route that round 7 closed. Narrowed in
the README, the round-7 write-up and the audit.

**The guidance velocity is not unique.** `docs/born_rule_equivariance.md` gave
as reason (iii) that `v = J/ρ` "is the unique velocity field whose current
closes the continuity equation". It is not: for any `K` with `∇·K = 0`,
`∇·(ρv') = ∇·(J+K) = ∇·J`, so `v' = (J+K)/ρ` closes the same equation wherever
`ρ > 0`. Shown with an explicit compactly supported `K = ∇×(f ẑ)` on a Gaussian
wavepacket — `|∇·K|/|∇·J| = 5e-16`, `|∇·(J+K) − ∇·J|/|∇·J| = 9e-16`. An
Ehrenfest or mean-velocity check does not exclude it either, since
`∫K d³x = 0` for compactly supported `∇×A` (`1.7e-16` relative). Reasons (i)
and (ii) stand, and Theorem 2's uniqueness over *densities* `h(ρ)` is a
different statement and unaffected; but the cited Goldstein–Struyve result
assumes Bohmian dynamics and locality conditions rather than deriving them, so
it cannot stand in for (iii).

## Where this leaves the chain

Rounds 5–7 are arithmetically untouched and epistemically one notch weaker: the
density they all rest on follows from the closure axiom **plus** a choice of
conditioning variable, where the axiom's own phrasing supplies the reason for
that choice. The open list is unchanged — branch aggregation, sector
coefficients, readout, composition, and the source-readability hazard — with
the conditioning variable now recorded beside them rather than above them.
