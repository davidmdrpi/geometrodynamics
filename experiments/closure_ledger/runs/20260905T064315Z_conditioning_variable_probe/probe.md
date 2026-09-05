# Round 8 — is the `|D|` density forced by the closure set?

**10/10 checks pass.**

## Verdict

* **A_conditioning_variable** — `CONDITIONING_VARIABLE_IS_A_CHOSEN_INPUT_JUSTIFIED_BY_THE_PHASE_AXIOM`
* **ledger_change** — `coarea conditioning: derived -> chosen (the variable)`
* **downstream_numbers_unchanged** — `True`
* **B_kappa** — `route-local parameter, not a fourth universal input`
* **C_velocity_uniqueness** — `born_rule_equivariance.md:73 reason (iii) is withdrawn`

## Checks

* PASS  O1   |N| is identical in all four outcome sectors
* PASS  O2   |grad N| is constant on Gamma -> uniform arclength
* PASS  O3   |grad theta| = |q|/|D| -> the |D| density
* PASS  O4   the N-window gives E = 0 exactly
* PASS  O5   the phase window tracks the repository law
* PASS  BK   same support, two limits (Borel-Kolmogorov)
* PASS  R1   the justification is in the repository, cited from source
* PASS  R3   NO round 5-7 number moved
* PASS  C1   div K = 0 leaves continuity unchanged
* PASS  C2   int K = 0, so a mean-velocity check cannot exclude it

## The two conditioning families

| eps | E, `\|N\| < eps` window | E, phase window |
|---|---:|---:|
| 0.05 | +0.000000 | +0.397021 |
| 0.02 | +0.000000 | +0.400039 |
| 0.01 | +0.000000 | +0.398457 |
| 0.005 | +0.000000 | +0.394911 |
| 0.002 | +0.000000 | +0.398601 |

Limit of the phase window (closed form): `0.3984966504`. The `N`-window limit is `0` at every `eps`, because `|N|` does not depend on the outcome signs.

Counts at `eps = 0.01`: `{'(+1,+1)': 5956, '(+1,-1)': 5956, '(-1,+1)': 5956, '(-1,-1)': 5956}` — equal, not approximately equal.

## The justification, read out of the source

* `geometrodynamics/history/closure.py` line 11: `3. Phase closure: total phase around every loop ≡ 0 or π (mod 2π)`
* the closure axiom is a condition on total phase, so a phase tolerance is the natural regularisation of it; N = 0 is a derived description of the same locus

## Rule 3 — nothing downstream moved

| quantity | expected | got | delta |
|---|---:|---:|---:|
| round-5 correlation(abs) at gamma=1 | +0.3984966504 | +0.3984966504 | 2.2e-11 |
| round-5 correlation(signed) at gamma=1 | +0.5403023059 | +0.5403023059 | 3.2e-11 |
| round-6 oriented singlet E at gamma=1 | -0.5403023059 | -0.5403023059 | 0.0e+00 |
| round-7 M_0 | +11.6708058364 | +11.6708058364 | 0.0e+00 |
| round-7 M_pi | +0.1695122784 | +0.1695122784 | 0.0e+00 |

## C — the velocity field is not unique

* `max|div K| / max|div J|` = `5.0e-16`
* `max|div(J+K) - div J| / max|div J|` = `8.7e-16`
* `int K d^3x` relative to `int|J|` = `1.7e-16` — a mean-velocity check cannot exclude it either

## Corrected dependency ledger

* `Haar prior on S^2` — **chosen** (invariant prior; physicality of x is itself a choice)
* `closure axiom: total phase = 0 or pi (mod 2pi)` — **repository axiom** (history/closure.py:12)
* `geodesic-realignment detection model` — **chosen** (round 5)
* `the conditioning VARIABLE (phase, not N)` — **chosen** (was ledgered 'derived'. Conditioning on a measure-zero set is not fixed by the set: an N-window gives the uniform measure and E = 0, a phase window gives |D| and the repository's law. Justified by the axiom being stated on phase, but justified is not derived)
* `coarea density |D|/(2|u x v|) GIVEN the phase variable` — **derived** (1/|grad theta| with grad theta = q/D)
* `outcome signs as history boundary data` — **chosen** (D-type)
* `equal prior on the four outcome sectors` — **chosen** (counting measure; only r = 1 gives the quantum law)
