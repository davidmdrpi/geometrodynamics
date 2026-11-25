# Prediction for the lowest Alice-handle doublet

This protocol concerns the ground pair on the unit handled surface in #355. The historical cells at a=0.3,0.5 and m=1,2,3 are development data. They will not score the prediction. The aperture is a in radians. Use the relative frequency split |ωD−ωN|/[(ωD+ωN)/2].

Conformal distance from the equator to a rim and from the rim to the waist is L=log cot(a/2) in each case. The decaying tail of the equator-normalized sphere mode sech(x)^m has amplitude 2^m exp(−mx). Reflection at conformal distance 2L doubles the value of the even tail and doubles the derivative of the odd tail. Their boundary Wronskian divided by I_m=∫₀^∞ sech(x)^(2m+2) dx gives the eigenvalue splitting. Dividing by 2m(m+1) predicts the relative frequency split

    D_pred = C(m) exp(−4mL)
    C(m) = 2 · 4^m / [(m+1) I_m].

This is the fixed-m, small-aperture leading law. No fitted coefficient, power, or aperture correction will enter the scored comparison. On all 16 cells of a={0.08,0.12,0.18,0.24} crossed with m={1,2,3,4}, require |D_pred/D_meas−1| ≤ 0.03 individually. Additionally the absolute relative error at a=0.08 must be smaller than at a=0.24 for each m. Report all cells and the worst discrepancy. Report a=0.4,0.6 crossed with m=2,4 separately as finite-aperture stress cases; these do not decide the primary gate.

Compute the continuum splitting by a normalized waist Wronskian for independent Neumann and Dirichlet solutions, without subtracting nearly equal double-precision frequencies. Cross-check by directly subtracting arbitrary-precision roots using a sphere-special-function/neck series formulation at (0.08,4), (0.24,1) and (0.4,2); require relative disagreement below 1e−7. Any failed gate is reported with its original threshold. #355's frozen COUNTEREXAMPLE remains the verdict for that PR.

## Numerical predictions

| a | m=1 | m=2 | m=3 | m=4 |
|---:|---:|---:|---:|---:|
| 0.08 | 1.539281524e-05 | 1.31632645e-10 | 1.181948241e-15 | 1.091610654e-20 |
| 0.12 | 7.813446059e-05 | 3.391663295e-09 | 1.545867062e-13 | 7.247129342e-18 |
| 0.18 | 0.0003979427104 | 8.797688931e-08 | 2.042236104e-11 | 4.876157823e-15 |
| 0.24 | 0.001268360825 | 8.937439901e-07 | 6.612607543e-10 | 5.032303415e-13 |
| 0.4 | — | 5.702062764e-05 | — | 2.048351746e-09 |
| 0.6 | — | 0.001676781688 | — | 1.771306002e-06 |
