# Recovering the local Euler polynomial from a limiting second difference

## The observable and the inverse problem

Take $`K=\mathbf Q(\zeta_{16})^+`$, $`F=K(\sqrt{397669})`$ and $`E=K(\sqrt{-547},\sqrt{-727})`$. With $`a=\log2`$, $`r=\sqrt2`$, $`c=2(r-1)^2`$ and $`B=rc`$, use the inert-minus-split score of width two,

```math
S(y)=\frac1{B\sqrt y}\sum_v-\epsilon_v\log(Nv)\,W(Nv/y),\qquad
W(u)=\max(0,\min(2-2/u,4/u-1)).
```

Here $`\epsilon_v=1`$ for split and $`-1`$ for inert in $`E/F`$. Remove places above $`2,547,727,343057911569,343070636977`$. The theta trace supplies these signs at the remaining places. For every real $`h`$ define

```math
C(h)=\limsup_{y\to\infty}\left(2S(y)-S(e^hy)-S(e^{-h}y)\right).
```

**Theorem.** Assume the generalised Riemann hypothesis (GRH) for the primitive Dirichlet functions in the 547/727 race and its base field, including the Riemann zeta function. For every rational prime $`p`$, six values of $`C`$ at specified log dilations determine the local polynomial $`D_p(T)=\prod_\chi(1-\chi(p)T)`$ of the eight primitive Dirichlet factors of $`\zeta_E/\zeta_F`$. The Euler factor is $`D_p(p^{-s})^{-1}`$. The reconstruction uses known Archimedean data, the constant $`V(0)`$ below and, for odd $`p`$, the degree determined by the conductors. Three further values at an odd prime recover that degree without using the local conductor count. This is a statement about the primitive product, including its local factors at primes omitted from the original prime-only score. [Field and trace](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-short-shifts/PROOF.md) [Trace identity](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-short-shifts/RaceTrace.lean#L83)

The explicit-formula reconstruction is unconditional. For each nontrivial zero $`\rho`$, evaluate the Fourier transform of the translated test function at $`(\rho-1/2)/i`$. Write $`V_{\mathrm{EF}}(h)`$ for the resulting zero sum. With the same conductor and Archimedean data, six values of $`2(V_{\mathrm{EF}}(0)-V_{\mathrm{EF}}(h))`$ determine $`D_p(T)`$ by the prime-side identity below. Under GRH the spectral parameters are real, $`V_{\mathrm{EF}}=V`$, and the six inputs are values of the limiting curvature $`C`$. [Explicit formula, Theorem 7.8](https://www.math.u-bordeaux.fr/~hecohen/rennes2.pdf#page=54) [GRH announcement](https://alpo.ge/remark.pdf) [Independent formal checks](https://math.mit.edu/~drew/grh/)

## The exact extremum at all log dilations

The score at log centre $`s=\log(2y)`$ equals $`1+V(s)+o(1)`$, where

```math
G(t)=\frac{3-2r\cos(at)}{c(t^2+1/4)},\qquad
V(s)=\sum_\nu G(\nu)e^{i\nu s}.
```

The real, even series uses the zero ordinates of the eight primitive factors with multiplicity and converges absolutely and uniformly. Under GRH, which Alpöge's announced Dirichlet GRH theorem supplies, their zeros lie on the critical line. The unshifted prime-only error tends to zero uniformly on each higher tail. [Proof and exact upper endpoint](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-isolated-peaks/PROOF.md) [Independent formal checks](https://math.mit.edu/~drew/grh/)

For a fixed $`h`$, the profile of the centred second difference has Fourier coefficients $`2G(\nu)(1-\cos(\nu h))\ge0`$. Their sum is $`2(V(0)-V(h))`$. This supplies an upper bound for the absolute value at every centre. Arbitrarily late simultaneous returns of the phases to one attain the upper bound in the limit. More precisely, if $`S(y_n)\to1+V(0)`$ and $`y_n\to\infty`$, the local-shape theorem gives convergence to that bound uniformly for $`h`$ in each compact interval. In particular,

```math
C(h)=\limsup_{y\to\infty}|2S(y)-S(e^hy)-S(e^{-h}y)|
     =2(V(0)-V(h)).
```

The same sequence of centres works at all finitely many dilations prescribed by the inverse theorem. This conclusion uses positivity and compactness of the phase orbit. It permits arbitrary rational relations among the ordinates. [Local-shape theorem](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-isolated-peaks/PROOF.md) · [Finite extremal inequality](CurvatureEuler.lean)

## Remove the Archimedean term

Put $`g(x)=4\sinh(\max(0,a-|x|)/2)/B`$ and $`A_*=\sum_\chi\log(q_\chi/\pi)-8(\gamma_E+a)`$. For real $`s`$ put $`z=e^{-|s|}`$ and

```math
T(s)=\frac{(1+z)\log(1+z)+(1-z)\log(1-z)}{2\sqrt z},\quad T(0)=a,\qquad
\Gamma_*(h)=A_*g(h)-\frac8c\{r[T(h+a)+T(h-a)]-3T(h)\}.
```

The eight conductors are $`547,4376,8752,8752,727,5816,11632,11632`$. The compact-test explicit formula gives $`V(0)=\Gamma_*(0)`$. For every $`h\gt 0`$ the corrected function is

```math
P(h)=\tfrac12 C(h)-V(0)+\Gamma_*(h)
    =\sum_{n=p^j}\beta_n\,g(h-\log n),\qquad
\beta_{p^j}=\frac{\log p}{\sqrt{p^j}}\,u_j(p),\quad
u_j(p)=\sum_\chi\chi(p)^j .
```

To see the signs, the even explicit formula expresses $`V(h)`$ as $`\Gamma_*(h)`$ minus the forward and reflected prime sums. The reflected sum vanishes for $`h\gt 0`$ since $`-\log n-h\lt -a`$ for $`n\ge2`$. At $`h=0`$ the first possible prime is at the kernel boundary where $`g(a)=0`$. The gamma correction follows by expanding $`e^{-t/2}/(2\sinh t)`$ into exponentials and convolving each with $`g`$. Partial fractions give the displayed three translates of $`T`$. These convolution identities hold at every real argument. Only finitely many prime powers contribute locally. [Even explicit formula and closed gamma kernel](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-distant-interval/PROOF.md) · [Archimedean correction](curvature.py)

## Three values give one derivative jump exactly

For $`n=p^j`$ choose $`x=\log n`$ and $`d=\log(1+1/(4n))`$. All possible derivative discontinuities of the prime sum occur at scales $`m/2,m,2m`$ for an integer prime power $`m`$. These scales are half-integers. The open scale window, including the three sampling points, lies inside $`(n-1/2,n+1/2)`$ since its endpoints are $`n/(1+1/(4n))`$ and $`n+1/4`$. Its only possible knot is exactly $`n`$. Both flanking segments solve $`P''=P/4`$. The kernel is continuous at its three knots.

Write $`J_n=P'_+(x)-P'_-(x)`$. Evaluating the two local exponential solutions gives $`P(x\pm d)=P(x)\cosh(d/2)\pm2P'_\pm(x)\sinh(d/2)`$. Consequently,

```math
J_n=\frac{P(x+d)-2\cosh(d/2)P(x)+P(x-d)}{2\sinh(d/2)} .
```

This identity removes every smooth prime-power contribution on the two flanks exactly. It includes knots of zero jump. The sampling points in the scale variable are rational, which permits direct comparisons with kernel boundaries. [Three-point identity in Lean](CurvatureEuler.lean) · [Sampling and decoding](recover.py)

## Odd primes and the ramified prime two

The derivative jumps of $`g`$ at $`-a,0,a`$ are $`r/c,-3/c,r/c`$. At an odd prime power $`n=p^j`$, the other possible contributing centres $`n/2`$ and $`2n`$ are respectively a noninteger and a number divisible by two distinct primes. Thus

```math
u_j(p)=-\frac{c\sqrt{p^j}}{3\log p}J_{p^j}.
```

All nonzero character values are fourth roots of unity and conjugation preserves the multiset. Write $`A,B_1,D`$ for the multiplicities of $`1,-1,i`$ respectively, with $`-i`$ also having multiplicity $`D`$. If $`w`$ is the number of nonzero characters at $`p`$, power sums of orders 1 and 2 give $`u=A-B_1`$, $`v=A+B_1-2D`$, and $`w=A+B_1+2D`$. Therefore

```math
A=(w+v+2u)/4,\qquad B_1=(w+v-2u)/4,\qquad D=(w-v)/4,\qquad
D_p(T)=(1-T)^A(1+T)^{B_1}(1+T^2)^D.
```

For odd $`p`$ the conductors give $`w=4`$ at 547 or 727 and $`w=8`$ otherwise. This proves the six-sample formula. If the local degree is unknown, use $`w=u_4(p)`$ and sample near $`4\log p`$ as well. [Multiplicity identities](CurvatureEuler.lean) · [Original primitive characters](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-factor-400/PROOF.md)

At $`p=2`$ set $`u_0=0`$ by convention for the absent integer-one prime term. The kernel and the neighboring dyadic powers give

```math
z_j=\frac{c\,2^{j/2}}{\log2}J_{2^j}
    =u_{j+1}(2)-3u_j(2)+2u_{j-1}(2).
```

Only the primitive odd-conductor quadratic characters can be nonzero at two, and their values are real. Write $`U=u_1(2)=u_3(2)`$ and $`W=u_2(2)`$. Then $`z_1=W-3U`$ and $`z_2=3U-3W`$, giving $`W=-(z_1+z_2)/2`$ and $`U=(W-z_1)/3`$. The two multiplicities are $`(W+U)/2`$ and $`(W-U)/2`$. The actual data give $`z_1=2,z_2=-6`$, recovering $`D_2(T)=1-T^2`$. The correction at $`\log2`$ removes the Archimedean cusp before the jump is evaluated. [Two-adic inversion](CurvatureEuler.lean) · [Sample enclosures](evidence/curvatures.tsv)

## Precision and independent checks

If each input value of $`C`$ has absolute error at most $`\epsilon`$ and the explicit constants are known exactly, each corrected value of $`P`$ has error at most $`\epsilon/2`$. The jump error is at most $`\tfrac12\coth(d/4)\epsilon`$. An odd-prime power sum can therefore be rounded unambiguously if

```math
\epsilon\lt \frac{3\log p}{c\sqrt n}\tanh(d/4),\qquad n=p^j .
```

For the nine-sample method the common sufficient accuracy $`\epsilon\lt (\log p)/(20p^6)`$ follows from $`d\ge1/(4n+1)\ge1/(5n)`$, $`\tanh(d/4)\ge d/8`$ and $`c\lt 1`$. For instance $`\tanh x\ge x/2`$ on $`[0,1/4]`$ by integrating $`1/\cosh^2x\ge1/2`$ there. The same reasoning with $`j\le2`$ gives $`\epsilon\lt (\log p)/(20p^3)`$ for the six-sample odd-prime method. At two, require $`\epsilon\lt (\log2)\tanh(d/4)/(c\sqrt n)`$ at $`n=2,4`$ to round $`z_j`$. The included implementation instead carries intervals and extracts a unique enclosed integer. [Rounding lemma](CurvatureEuler.lean) · [Interval decoder](recover.py)

The sample table contains 81 values at ten rational primes, each in an outward interval of width $`10^{-36}`$. The generator enumerates every prime power in the compact supports, evaluates the quadratic characters by modular exponentiation and the quartic sum by reduction modulo 16, and adds the gamma terms. The independent check evaluates $`T`$ by its gamma series $`\sum_{k\ge0}z^{2k+3/2}/[(2k+1)(2k+2)]`$ with a rigorous tail, recovers the two one-sided slopes, and compares the jumps with the eight primitive Conrey characters. The rank inversion is checked for all 95 conjugation-stable fourth-root multisets of degree at most eight. The sample values are evaluations of the proved explicit formula. The time or height required to estimate a limiting extremum directly from prime counts is a separate problem. [Sample generator](sample_values.py) · [Independent arithmetic](independent.py) · [Replay transcript](evidence/check.txt)
