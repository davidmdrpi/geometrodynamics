# Correction note to the R3 preregistration (freeze `2e984ac`)

Date: 2026-09-29. Published before any implementation or measurement. The
freeze document is unchanged. This note amends one procedural step. No
gate, threshold, ladder or decision rule changes.

## 1. Bracket widening in the centre-manifold bisection (section 4, step 3)

The freeze fixes bisection brackets of ±1e-2 for the first window and
±1e-6 after that. On the centre manifold, A at the section differs from 1 by
O(eps^2). At eps = .16, that offset can exceed .01. The two ends of the
first bracket would then classify alike, and the run would fail for a
reason unrelated to the question.

Amended rule: if both bracket ends classify alike, widen the bracket
symmetrically by a factor of 4 and retry.
- The first window may widen to at most ±.32.
- Later windows may widen to at most ±1e-3.
- If the bracket still fails to straddle, the run fails, and N3 treats it as
  failing, as before.
- The widening count of every window is recorded in the archive.

## 2. Event filtering

The section event q0 = 0 is also satisfied at the start of every
integration. Any section event within conformal time 1 of the integration
start is ignored. This is an implementation detail, needed so that the
registered definition of a clock period holds.
