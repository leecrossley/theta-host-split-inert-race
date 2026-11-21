# A short prime-power sum from a prime-race extremum

## Statement and conventions

Let $`K=\mathbf Q(\zeta_{16})^+`$, $`F=K(\sqrt{397669})`$ and $`E=K(\sqrt{-547},\sqrt{-727})`$. For a prime ideal $`v`$ of $`F`$ put $`\epsilon_v=1`$ for split and $`-1`$ for inert in $`E/F`$. Omit places above $`2,547,727,343057911569,343070636977`$. The normalized first-step theta traces give $`\epsilon_v`$ on all the other places. The eight primitive factors $`\chi`$ of $`L=\zeta_E/\zeta_F`$ are the quartic even twists of the quadratic characters of conductors 547 and 727. They are all odd and are closed under conjugation. Write $`u_j(p)=\sum_\chi\chi(p)^j`$ with primitive local values, including zero when a character is ramified. [Field factorization](https://github.com/leecrossley/theta-race-short-shifts/blob/483c31f2df2132ed2165620ea55e20e969228bfb/PROOF.md) · [Trace lemma](https://github.com/leecrossley/theta-race-short-shifts/blob/483c31f2df2132ed2165620ea55e20e969228bfb/RaceTrace.lean#L83)

Fix $`b\gt 1`$, $`a=\log b`$, $`r=\sqrt b`$, $`c=2(r-1)^2`$. Define

```math
W_b(t)=\max(0,\min(b-b/t,b^2/t-1)),\qquad
S_b(y)=-\frac{1}{rc\sqrt y}\sum_v\epsilon_v\log(Nv)W_b(Nv/y).
```

For $`0\lt h\lt \log2`$ put $`d=\cosh(h/2)`$ and define two limiting observables of the prime-only score,

```math
J_b(h)=\limsup_{y\to\infty}\{2dS_b(y)-S_b(e^hy)-S_b(e^{-h}y)\},\qquad
C_b(h)=\limsup_{y\to\infty}\{2S_b(y)-S_b(e^hy)-S_b(e^{-h}y)\}.
```

Let $`x_+=\max(x,0)`$. For $`n=p^j`$ write $`\beta_n=(\log p)u_j(p)/\sqrt n`$. Define the known background $`B_b(h)=2(d-1)+2[d\Gamma_b(0)-\Gamma_b(h)]`$ using the explicit gamma term below. Assume the generalised Riemann hypothesis (GRH) for the primitive Dirichlet functions in the 547/727 race and its base field, including the Riemann zeta function. Then for every allowed width and displacement,

```math
J_b(h)-B_b(h)=\frac{4r}{c}\sum_{n=p^j}\beta_n\,
 \sinh\!\left(\frac{(h-|\log(b/n)|)_+}{2}\right).
```

The band at $`b=n=p^j`$, $`h=\log(1+1/(2n))`$ contains exactly one integer, since its upper endpoint is $`n+1/2`$ and its lower endpoint exceeds $`n-1/2`$. Hence one measurement gives

```math
u_j(p)=\frac{c\,[J_{p^j}(h)-B_{p^j}(h)]}{4(\log p)\sinh(h/2)}.
```

This includes negative traces, invisible prime powers with zero trace and all conductor primes. Absolute error strictly below $`2(\log p)\sinh(h/2)/c`$ in the corrected extremum permits rounding to the unique integer. At odd primes the local degree $`w`$ is eight except at 547 and 727, where it is four. At the prime two it is two. From $`u=u_1(p)`$ and $`v=u_2(p)`$ the multiplicities of $`1,-1,i,-i`$ are $`(w+v+2u)/4,(w+v-2u)/4,(w-v)/4,(w-v)/4`$. Thus two limiting extrema at widths $`p,p^2`$ determine the Euler factor, or three at widths $`p,p^2,p^4`$ determine it without the degree. This counts extremal values of prime sums. [Fixed-width inverse theorem and conductors](https://github.com/leecrossley/theta-race-curvature-euler-factors/blob/af7093b56a7f6b4155dba27e735b5b1d471ff707/PROOF.md) · [Certified examples](aperture.py)

Alpöge's announced Dirichlet GRH theorem supplies this input. [Announcement](https://alpo.ge/remark.pdf) [Independent formal checks](https://math.mit.edu/~drew/grh/) The cusp law for $`C_b`$ below uses the same assumption. For each nontrivial zero $`\rho`$ set $`t_\rho=(\rho-1/2)/i`$ and write $`V_b^{\mathrm{EF}}(s)=\sum_\rho G_b(t_\rho)e^{it_\rho s}`$, counting multiplicities. The prime-side identity gives the exact band formula and the $`h\log(1/h)`$ expansion for the corresponding expressions in $`V_b^{\mathrm{EF}}`$ unconditionally. Under GRH each $`t_\rho`$ is real, $`V_b^{\mathrm{EF}}=V_b`$, and these expressions give the limiting extrema $`J_b`$ and $`C_b`$.

## The prime-only limit at a fixed width

The log centre is $`s=\log(by)`$. Set

```math
g_b(x)=\frac{2r}{c}\sinh\!\left(\frac{(a-|x|)_+}{2}\right),\qquad
G_b(t)=\frac{b+1-2r\cos(at)}{c(t^2+1/4)},\qquad
V_b(s)=\sum_\nu G_b(\nu)e^{i\nu s}.
```

Here $`\nu`$ runs over the nontrivial zero ordinates of the eight primitive factors, with multiplicity. Direct Mellin integration gives the displayed transform $`G_b`$ of $`g_b`$. The numerator is at least $`(r-1)^2`$, so every weight is positive. Under GRH, which Alpöge's announced Dirichlet GRH theorem supplies, the zeros lie on the critical line. The zero count and the quadratic denominator give absolute uniform convergence. Complex conjugation makes the series real and even. [Mellin transform at general width](https://github.com/leecrossley/theta-race-short-shifts/blob/483c31f2df2132ed2165620ea55e20e969228bfb/PROOF.md#euler-factors-and-the-kernel) [Independent formal checks](https://math.mit.edu/~drew/grh/)

For each fixed width, $`S_b(e^s/b)=1+V_b(s)+o_b(1)`$. Here are the details extending the dyadic inversion. The transform of the centred score is the continuation of the transform of $`g_b`$ times the signed prime series at $`1/2+z`$. Split the prime series into $`L'/L(1/2+z)-\zeta_F'/\zeta_F(1+2z)`$ and a Dirichlet series whose absolute coefficient mass at $`7/20`$ is finite. The latter starts at third prime-ideal powers, apart from finitely many local corrections which may equivalently be handled before sieving. Bounds for each of those finitely many Euler-factor series are geometric for positive real part. On the fixed compact support its prime contribution after normalization is $`O_b(e^{-3s/20})`$. The original omitted primes eventually leave the support.

The pole of $`\zeta_F`$ gives $`G_b(0)/2=1`$. A base zero has normalized real part $`-1/4`$, and the transformed kernel on that vertical line is $`O_b((1+t^2)^{-1})`$. Its total contribution is $`O_b(e^{-s/4})`$. Trivial zeros give a geometrically convergent sum once $`s\gt a+1`$. Equivalently first convolve with a compact four-box approximate identity. Dominated convergence removes it from the zero sums using quadratic decay on the two fixed vertical lines, and from the original continuous compact prime sum. This proves the assertion at fixed width. On every compact interval of widths inside $`(1,\infty)`$ the supports, transform bounds and error majorants can be chosen uniformly. [Prime-square cancellation and remainder](https://github.com/leecrossley/theta-race-fixed-factor/blob/12b7c94cfb9586c25fb29e7c7aaf0a20adaf5b2f/PROOF.md#global-input-and-transform) · [Unshifted inversion](https://github.com/leecrossley/theta-race-isolated-peaks/blob/ee35ecfd698ab0ad193669a8109bbbfb1aac6049/PROOF.md#the-unshifted-prime-profile)

For $`J_b`$ the profile has constant $`2(d-1)`$ and coefficients $`2G_b(\nu)(d-\cos(\nu h))\ge0`$. For $`C_b`$ replace $`d`$ by one. The sum of the coefficients bounds the profile above. Arbitrarily late simultaneous returns of finitely many phases to one, followed by an absolute tail bound, show that the upper bound is exact. This argument uses no linear independence of zeros. In particular

```math
J_b(h)=2(d-1)+2[dV_b(0)-V_b(h)],\qquad
C_b(h)=2[V_b(0)-V_b(h)].
```

## The local prime-power sum

Put $`A_*=\sum_\chi\log(q_\chi/\pi)-8(\gamma_E+\log 2)`$. The eight conductors are $`547,4376,8752,8752,727,5816,11632,11632`$. With $`z=e^{-|x|}`$ define

```math
T(x)=\frac{(1+z)\log(1+z)+(1-z)\log(1-z)}{2\sqrt z},\quad T(0)=\log2,\qquad
\Gamma_b(x)=A_*g_b(x)-\frac8c\{r[T(x+a)+T(x-a)]-(b+1)T(x)\}.
```

The even explicit formula gives

```math
V_b(x)=\Gamma_b(x)-\sum_{n=p^j}\beta_n\,[g_b(x-\log n)+g_b(x+\log n)].
```

For clarity, the conductor term before integration is $`[\sum_\chi\log(q_\chi/\pi)+8\psi(3/4)]g_b(x)`$. The gamma integral subtracts eight times the integral of $`g_b(x-t)+g_b(x+t)-2g_b(x)`$ against $`e^{-t/2}/(2\sinh t)`$ on $`t\gt 0`$. Expand the latter into exponentials with exponents $`v_k=2k+3/2`$. The second differences are Lipschitz at zero, so their absolute integrals against the individual exponentials are summable. Off its three knots $`g_b''=g_b/4`$. Its derivative jumps at $`-a,0,a`$ are $`r/c,-(b+1)/c,r/c`$. Convolution against $`e^{-v_k|x|}`$ therefore contributes the three displayed translations divided by $`v_k^2-1/4`$, apart from the multiple of $`g_b(x)`$. The constant series $`\sum_k[v_k(v_k^2-1/4)]^{-1}=\pi-4\log2`$ combines with the digamma value to give $`A_*`$. The translation series is $`T(x)=\sum_k e^{-v_k|x|}/(v_k^2-1/4)`$, whose elementary sum is displayed above. These are the same odd gamma factors at every width. [Convolution calculation](https://github.com/leecrossley/theta-race-distant-interval/blob/527fffe7b7808457b47a292144fea96b556158a8/PROOF.md) · [Directed implementation](width_kernel.py)

Because $`0\lt h\lt \log2\le\log n`$, both translates of $`\log n`$ stay positive. Set $`v=\log(b/n)`$. Inserting the explicit formula in the exact expression for $`J_b`$ leaves, for each prime power, $`2\beta_n`$ times

```math
g_b(\log n-h)+g_b(\log n+h)-2d\,g_b(\log n)
=\frac{2r}{c}\sinh\!\left(\frac{(h-|v|)_+}{2}\right).
```

For $`v\le-h`$ all three terms vanish. For $`v\ge h`$ the hyperbolic addition formula cancels them. For $`-h\le v\le0`$ the only positive argument is $`v+h`$. For $`0\le v\le h`$ subtraction leaves $`\sinh((h-v)/2)`$. This proves the exact band formula, including equality at the two endpoints. At $`b=n`$ its central term simplifies because $`r=\sqrt n`$. [Supporting kernel identities in Lean](WidthCusp.lean)

## The cusp at every width

Let $`h\downarrow0`$ and remove the universal terms by putting

```math
R_b(h)=\frac{c C_b(h)}{2h}-4(b+1)\log(1/h)-\frac{(8+A_*)(b+1)}2 .
```

For each $`b_0\gt 1`$ define $`m(b_0)=(\log p)u_j(p)`$ if $`b_0=p^j`$ and zero otherwise. The limit for fixed width and the simultaneous limit for a moving width are

```math
\lim_{h\downarrow0}R_{b_0}(h)=m(b_0),\qquad
\lim_{h\downarrow0}R_{b_0 e^{\lambda h}}(h)
       =m(b_0)(1-|\lambda|)_+ .
```

The second convergence is uniform for $`\lambda`$ in each bounded interval. In particular $`C_b(h)\sim [8(b+1)/c]h\log(1/h)`$ for every fixed $`b\gt 1`$. The ordinary second difference is even in $`h`$, giving the same asymptotic from either direction with $`|h|`$.

To prove this, the formula for $`T`$ gives $`T(h)=\log2+\tfrac h2(\log h-1)+O(h^2|\log h|)`$. Indeed $`1-e^{-h}=h+O(h^2)`$ and the term from differentiating the denominator cancels the linear multiple of $`\log2`$. The pair $`T(a+h)+T(a-h)-2T(a)`$ is $`O(h^2)`$ since $`a\gt 0`$. Direct differentiation gives $`g_b(0)-g_b(h)=(b+1)h/(2c)+O_b(h^2)`$. All of these error bounds are uniform when $`b`$ ranges in a compact neighbourhood of $`b_0`$ inside $`(1,\infty)`$. Also $`V_b(0)`$ is uniformly bounded there by the summable spectral majorant. Consequently $`J_b(h)-C_b(h)=2(d-1)(1+V_b(0))=O_{b_0}(h^2)`$. The same estimates in the known background leave

```math
R_b(h)=\frac{2\sqrt b}{h}\sum_{n=p^j}\beta_n
 \sinh\!\left(\frac{(h-|\log(b/n)|)_+}{2}\right)+O_{b_0}(h|\log h|).
```

For bounded $`\lambda`$ and $`b=b_0e^{\lambda h}`$, discreteness of the integers leaves only the possible term $`n=b_0`$ when $`h`$ is sufficiently small. Its normalized contribution is $`(\log p)u_j(p)(1-|\lambda|)_++O_{b_0,\lambda}(h)`$ with the constant uniform on each bounded set of detunings. If $`b_0`$ is not a prime power the band eventually contains none. This also handles vanishing traces. [Complete-profile checks at 243 width and displacement pairs](cusp_check.py)

## Verification

The first implementation sums the full prime profile up to each compact endpoint using modular powers and prefix sums. It gives 32-decimal outward enclosures of 20 extremal values from ten primes. The second implementation computes the known background using a convergent gamma series away from zero and separate hyperbolic logarithms near zero. It recovers all 20 traces from the stored extrema using the eight primitive Conrey characters. Seven other intervals, including an empty interval and one containing 112 prime powers, check the whole-band formula to 83 decimal places. The moving-width test includes both edges $`\lambda=\pm1`$, zero, positive and negative traces, and displacements down to $`2\log(1+10^{-18})`$. These finite checks accompany the analytic proof for all widths. [Independent calculation](independent.py) · [Exact sample intervals](evidence/weighted-extrema.tsv) · [Replay transcript](evidence/check.txt)
