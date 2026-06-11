# CKM intra-channel analogue from mouth-overlap alignment (PR #155)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity. The CKM
> matrix is read off the partition blocks of the mass-locked shelled-closure
> Hamiltonian.

The quark mirror of the #153 PMNS extraction. #91 made the structural claim:
quarks mix **within** one channel (up-type and down-type are the two Z₂
partition classes of the same cavity shells), so their mass bases are nearly
aligned and the CKM is small — in contrast to the cross-channel (large)
PMNS. This PR makes it quantitative, and the result is a genuine
**out-of-sample prediction**: the CKM matrix extracted from the LOCKED quark
Hamiltonian, whose parameters were calibrated on the six quark **masses**
alone. **Zero new inputs.**

## The construction

`LOCKED_QUARK_PARAMS` has `partition_mixing = 0`, so the 6×6 Hamiltonian is
exactly block-diagonal in the partition label: (+) = (u, c, t), (−) =
(d, s, b), over the shared shells k = 1, 3, 5. The CKM matrix is the
shell-space misalignment of the two partition eigenbases, `V = U₊†U₋` —
unitary to machine precision.

## The predicted matrix

| element | predicted | observed (PDG) | ratio |
|---|---|---|---|
| V_us | 0.112 | 0.225 | 0.50 |
| V_cd | 0.112 | 0.225 | 0.50 |
| V_cb | **0.0377** | 0.0418 | **0.90** |
| V_ts | **0.0372** | 0.0411 | **0.91** |
| V_td | 0.0063 | 0.0086 | 0.73 |
| V_ub | 0.0020 | 0.0037 | 0.55 |

Every off-diagonal element within a factor ≤ 2.0 (V_us/V_cd sit exactly at
×2.0 — the soft direction, see below); **V_cb and V_ts within 10%**; the
hierarchy |V_us| > |V_cb| > |V_ub| exactly reproduced; unitarity at machine
precision.

## The mechanism and the dichotomy

The up sector is aligned to 0.008 while the down sector carries the mixing
(0.12): the minus-partition's asymmetric couplings — the same terms the mass
calibration fixed to order the spectrum — order the mixing. The anarchic
cross-channel counterfactual (the PMNS situation) gives |V_us| ~ 0.46:
**large PMNS (#153) and small CKM (here) emerge from one framework** —
cross-channel anarchy vs intra-channel partition alignment.

## The stiffness audit

Under ±10% shifts of the locked couplings, |V_cb| moves **< 1%** (stiff —
the 10% agreement is a sharp falsifiable prediction) while |V_us| swings
×0.55–×3.3 under pinhole shifts (soft — the small d–s splitting, 3.34 vs
6.40, amplifies the sensitivity ~×8). V_us's factor-2 deficit sits inside
the soft direction's calibration slop; **V_cb/V_ts are the falsifiable
core**.

## CP

The locked baseline has `phase = 0` ⟹ J = 0 exactly, vs the observed
J ≈ 3×10⁻⁵: the quark Dirac phase is the open item — and a constrained one,
since a phase calibration must reproduce J **without disturbing the
already-fixed |V|**. Contrast the leptonic generic CP, which came free from
the anarchic phases (#153/#154).

## Ledger and scope

- **Derived (zero new inputs):** CKM smallness, hierarchy ordering, all
  elements ≤ ×2.0, V_cb/V_ts at 10% (stiff), up-aligned/down-mixing
  structure, the PMNS/CKM dichotomy.
- **Soft:** V_us precision (the d–s splitting direction).
- **Open:** the quark CP phase (J); the link from the partition-asymmetric
  couplings to the #152 mouth-overlap machinery. The #150 budget is
  untouched — this probe consumed no inputs at all.

## Reproduce

```bash
python -m experiments.closure_ledger.ckm_intra_channel_probe
# Verdict: CKM_SMALL_HIERARCHICAL_INTRA_CHANNEL_VCB_10PCT_ZERO_NEW_INPUTS
```
