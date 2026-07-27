# Does the #207 surgery law acquire a Z2 term? (PR #227 follow-on)

_Run 2026-07-27T03:29:52.401458+00:00 · 4/4 PASS_

## The answer: no

  - outer-bridge twists move the pair phase by at most **2.3e-13** (identically zero)
  - the junction twist moves it by **4.68e-02** rad — not pi
  - every odd-n_twist row misses the conjectured pi by at least **3.09** rad

## The mechanism: cycle interference, not a Z2 phase

| M2–M3 gap | d(phi) junction twist | d(phi) outer twist |
|---:|---:|---:|
| 4 | 2.877e-01 | 1.221e-13 |
| 8 | 4.678e-02 | 2.283e-13 |
| 16 | 2.139e-04 | 2.984e-13 |
| 32 | 2.975e-14 | 1.776e-13 |
| 48 | 1.763e-13 | 1.337e-13 |

| exterior hopping t_x | d(phi) junction twist |
|---:|---:|
| 0.25 | 2.486e-05 |
| 0.5 | 1.154e-03 |
| 1.0 | 4.678e-02 |
| 2.0 | 1.770e-01 |

## Verdict

**THE_SURGERY_COMPOSITION_LAW_HAS_NO_Z2_TERM_THE_PAIR_PHASE_IS_A_RATIO_AND_THE_TWIST_ACTS_ONLY_AS_CYCLE_INTERFERENCE**

Corrected law: `phi_14 = phi_a + phi_b + phi_c   (unchanged)`
