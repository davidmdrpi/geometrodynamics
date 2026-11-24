# The lowest Alice-handle frequency split

All 16 held-out cells passed the 3% relative-error bound registered at [73cca3f](https://github.com/davidmdrpi/geometrodynamics/commit/73cca3f014ab2d0d7c507d30b0e6fd63d5a9d185). The largest error was 1.0783%, at aperture 0.24 and angular winding 4. The two parity sectors have different waist conditions. The spherical band and the neck each supply a conformal distance over which their common ground mode decays. The result applies to the linear scalar on the fixed unit surface in [#355](alice_handle_ports.md). That experiment's six-pair threshold failed. This follow-up predicts the magnitude for the lowest pair in each winding.

## Deriving the coefficient

Write the meridional metric as $d\ell^2+r^2d\varphi^2$. With $dx=d\ell/r$ the mode equation becomes $-u''+m^2u=\lambda r^2u$, with $\lambda=\omega^2$. Let $x=0$ at the equator and increase $x$ toward the waist. Integration along the sphere gives $L=\log\cot(a/2)$ and $r=\operatorname{sech}x$. On the neck, $r^2=s^2+b^2$, $b=\sin^2 a$, so $r=b\cosh(2L-x)$ and the waist is at $x=2L$. Both the mode and its conformal derivative match at the rim $x=L$.

With Klein gluing, the cosine winding is periodic and the sine winding antiperiodic. The ground modes are even at the equator, with Neumann (N) and Dirichlet (D) waist conditions respectively. In the small-aperture limit the equator-normalized ground mode is $u_0=\operatorname{sech}^m x$ with $\lambda_0=m(m+1)$. Substitution gives both the eigenvalue and the asymptotic tail $A e^{-mx}$, $A=2^m$. Its half-sphere weighted norm is

$$
 I_m=\int_0^\infty\operatorname{sech}^{2m+2}x\,dx
 =\frac{\sqrt\pi\,\Gamma(m+1)}{2\Gamma(m+3/2)} .
$$

For fixed $m\ge1$ the neck equation has evanescent solutions. The integral of $r^2$ across the half-neck is $b^2[L/2+\sinh(2L)/4]=O(a^2)$. Variation of constants about $u''=m^2u$ therefore gives a vanishing relative correction to the growing waist-to-rim transfer coefficient. In matching from the equator, one can first fix a large $X$ on the spherical tail, then send $a$ to zero. On $[0,X]$ the ground solution approaches the sphere mode. Beyond $X$, the spherical weight has integral $O(e^{-2X})$ and the neck adds $O(a^2)$. Matching the incoming decaying tail and its waist reflection gives $u_N(2L)\sim2A e^{-2mL}$ and $u_D'(2L)\sim-2mA e^{-2mL}$. The factors of two follow by adding or subtracting the reflected exponential. The weighted overlap tends to $I_m$.

Subtract the two differential equations after cross-multiplication. The equator term vanishes and the exact boundary term is $-u_N(2L)u_D'(2L)$. Hence $\lambda_D-\lambda_N\sim4m A^2 e^{-4mL}/I_m$. Converting from squared frequencies to the registered relative frequency split gives

$$
 \frac{\omega_D-\omega_N}{(\omega_D+\omega_N)/2}
 \sim C(m)\tan^{4m}(a/2),\qquad
 C(m)=\frac{2\,4^m}{(m+1)I_m}.
$$

For example, $I_1=2/3$ and $C(1)=6$. The coefficient multiplying $a^{4m}$ is $C(m)/16^m$. The factor $e^{-4mL}$ includes the decay from the equator to the rim as well as the decay through the neck. No measured frequencies enter this prediction.

## Held-out results

Each entry is the signed percentage error $100(D_{\rm pred}/D_{\rm meas}-1)$. The full precision predictions, measured splits and parity eigenvalues are archived. The magnitude of the error decreases with aperture for each tested winding.

| Aperture | m=1 | m=2 | m=3 | m=4 |
|---:|---:|---:|---:|---:|
| 0.08 | −0.0101 | −0.0151 | −0.0201 | −0.0252 |
| 0.12 | −0.0423 | −0.0635 | −0.0846 | −0.1057 |
| 0.18 | −0.1697 | −0.2544 | −0.3388 | −0.4231 |
| 0.24 | −0.4330 | −0.6500 | −0.8642 | −1.0783 |

In the four wider-aperture stress cases, the prediction falls below the measured split by 2.84% and 4.64% at $a=0.4$ and by 5.72% and 9.14% at $a=0.6$ (respectively $m=2,4$). The accuracy degrades outside the registered small-aperture range.

## Two continuum calculations

The primary solver starts at the waist with N data $(1,0)$ or D data $(0,1)$, in a coordinate increasing toward the equator. It locates each equator-Neumann root. An unwrapped Prüfer phase verifies the ground-state index and that the competing equator-odd antiperiodic sector lies higher. The unit initial Wronskian yields $(\lambda_D-\lambda_N)\int r^2u_Nu_D\,dx=1$ exactly. Integrating the overlap resolves a split even when the two stored double-precision eigenvalues round to the same value.

A separate solver uses 45-decimal-digit arithmetic. It Taylor-propagates through the neck and matches to the even associated-Legendre solution on the band, expressed as a Gauss hypergeometric function. It subtracts the frequencies of the two independently located roots directly. The relative differences from the overlap solver were $8.72\times10^{-13}$ or smaller in all three registered checks. At $(a,m)=(0.08,4)$ the prediction was $1.091610654\times10^{-20}$ and the directly measured split was $1.091885331\times10^{-20}$.

Replay from the repository root with Python, SciPy, NumPy and mpmath.

~~~sh
python -m experiments.closure_ledger.alice_ground_split
python -m unittest discover -s tests -p test_alice_ground_split.py
~~~

[Registered prediction table](alice_ground_split_protocol.md) · [All 20 cells and three precision checks](../experiments/closure_ledger/runs/20261123_herring_split/probe.json)
