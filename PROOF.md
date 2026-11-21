# The factor 18 obstruction

## Two statements about the same race

Let $`K=\mathbf Q(\zeta_{16})^+`$, $`F=K(\sqrt{397669})`$ and $`E=K(\sqrt{-547},\sqrt{-727})`$. Let $`\epsilon_v=1`$ for split and $`-1`$ for inert base prime ideals. Omit the places above $`2,547,727,343057911569,343070636977`$. Put $`a=\log2`$, $`B=2\sqrt2(\sqrt2-1)^2`$, and

```math
W(u)=\max(0,\min(2-2/u,4/u-1)),\quad
S(y)=\frac{\sum_v-\epsilon_v\log(Nv)W(Nv/y)}{B\sqrt y}.
```

A gain or loss means $`S(2y)-S(y)\gt 1.9735a`$ or $`S(2y)-S(y)\lt -1.9735a`$. For a starting height $`H`$ and factor $`C`$, the permitted starts satisfy $`H\lt y\le CH/2`$.

1. For every real $`1\le H\le10^9`$, both moves finish by $`18H`$. For $`10^6\le H\le10^9`$, both finish by $`(13603772/2048491)H\lt 7H`$.
2. Assume the generalised Riemann hypothesis (GRH) for the primitive Dirichlet $L$-functions in the 547/727 race and its base field, including the Riemann zeta function. For arbitrarily large $`H`$ there is no loss finishing by $`2e^{11/5}H`$. For arbitrarily large $`H`$ there is no gain finishing by that factor. These may be different heights. Since $`2e^{11/5}\gt 18.05`$, neither $`18`$ nor $`18.05`$ works universally.

The finite count is independent of any global zero theorem. The second assertion invokes Alpöge's announced Dirichlet GRH theorem for primitive Dirichlet $`L`$-functions. Its strip statement uses Mathlib's analytic continuation of the Dirichlet series. Sutherland reports that all modules replayed, the solution export passed independent kernel checks and con-leche matched the exact primitive theorem against the challenge. Its axiom dependencies are the three standard Lean axioms. [Primitive theorem and its statement](https://github.com/leanprover-community/riemann-hypothesis/blob/371a1a26d5e17f604428e3a8ce6ce341cf517434/src/Statements.lean#L14) [Independent statement comparison](https://math.mit.edu/~drew/grh/con-leche-compare.log) [Verification record](https://math.mit.edu/~drew/grh/) [Concrete theta trace](https://github.com/leecrossley/theta-race-short-shifts/blob/483c31f2df2132ed2165620ea55e20e969228bfb/RaceTrace.lean#L83)

## The limiting oscillation

Use $`x=\log(2\sqrt2\,y)`$ for the log centre of a doubling and $`f(x)=[S(2y)-S(y)]/a`$. Write $`L=\zeta_E/\zeta_F`$. Its eight odd primitive factors are $`\chi_p\xi^j`$, $`p=547,727`$, $`j=0,1,2,3`$, where $`\chi_p=(\cdot/p)`$ and $`\xi`$ is the even quartic character of conductor 16 with $`\xi(5)=i`$. Write $`\nu`$ for the ordinates of their nontrivial zeros with multiplicity. The real kernel for $`f`$ and its bilateral transform are

```math
g(v)=\frac4B\sinh\!\left(\frac{\max(0,a-|v|)}2\right),\quad
j(v)=\frac{g(v+a/2)-g(v-a/2)}a,\quad
J(z)=\frac{2\sinh(az/2)}a
 \frac{2\sqrt2\cosh(az)-3}{2(\sqrt2-1)^2(z^2-1/4)}.
```

In particular $`j`$ is odd, with support $`[-j_0,j_0]`$, $`j_0=3a/2`$. Define

```math
U(x)=\sum_\nu J(i\nu)e^{i\nu x}.
```

On the critical line the coefficients are $`O((1+\nu^2)^{-1})`$ with summable majorant by the Dirichlet zero count. The series is uniformly and absolutely convergent. The set of factors is conjugation-stable. Pairing opposite ordinates shows that $`U`$ is real and odd.

We have $`f(x)-U(x)\to0`$ as $`x\to+\infty`$. To see this, use the prime-only decomposition $`D(s)=(L'/L)(s)-(\zeta_F'/\zeta_F)(2s)+D_R(s)`$ in the bilateral transform $`J(z)D(1/2+z)`$. The odd nontrivial zeros give $`U`$. Each base-field factor is also a primitive Dirichlet function, including ordinary zeta. Its nontrivial zeros in the square term have shifted real part $`-1/4`$, so that term is $`O(e^{-x/4})`$ by absolute summability. The base pole cancels at $`J(0)=0`$. The trivial-zero sums tend to zero for $`x\gt j_0`$. The remaining prime powers cost $`O(e^{-x/10})`$ by absolute convergence at $`2/5`$ and the compact window. The finitely many omitted places eventually leave the windows. One may first convolve with four narrow centred boxes and then let their widths tend to zero. The prime sums are continuous, the kernel multipliers have summable bounds, and the trivial sums converge uniformly on compact subsets of $`x\gt j_0`$. [Prime-only decomposition and explicit bounds](https://github.com/leecrossley/theta-race-waiting-bound/blob/9193ef9459374da9f68156ad6a1aa8db088882b0/PROOF.md)

## Evaluate the zero sum near the origin

The obstruction follows from a strict inequality on a closed real interval.

```math
U(x)\gt -9/5\quad (1.014\le x\le3.214).
\tag{1}
```

Near its left endpoint the support of the translated kernel crosses the origin. The conductor term and gamma integral must be retained there. We give two evaluations.

For each factor set $`\Lambda(s)=(q/\pi)^{(s+1)/2}\Gamma((s+1)/2)L(s,\chi)`$. Let $`P(x)`$ be the contribution of the Dirichlet logarithmic derivative:

```math
P(x)=\sum_{n=p^k}\frac{c_n\log p}{\sqrt n}\,j(x-\log n),\qquad
 c_{p^k}=-\sum_{\chi}\chi(p)^k .
```

For every real $`x`$ the symmetrized explicit formula gives

```math
U(x)=P(x)-P(-x)+A j(x)
 -8\int_0^\infty \frac{e^{-t/2}}{2\sinh t}
    [j(x-t)+j(x+t)-2j(x)]\,dt,
\quad
A=\sum_\chi\log(q/\pi)+8\psi(3/4).
\tag{2}
```

Here $`\psi=\Gamma'/\Gamma`$. For a derivation, integrate the completed logarithmic derivatives against $`J(z)e^{xz}`$ on the two sides of the critical strip in the shifted variable $`z=s-1/2`$. The functional equation for the conjugation-stable family and $`J(-z)=-J(z)`$ leave $`J(z)(e^{xz}-e^{-xz})`$ on the right line. The Euler series there gives $`P(x)-P(-x)`$. The constant conductor terms give $`j(x)\sum\log(q/\pi)`$. With $`u_k=3/2+2k`$ the gamma logarithmic derivative has the expansion

```math
\tfrac12\psi((z+3/2)/2)=-\gamma_E/2+
 \sum_{k\ge0}\left[\frac1{2(k+1)}-\frac1{z+u_k}\right].
```

Inverting the reciprocal kernels and subtracting the constant part at $`z=0`$ gives (2), since $`\sum e^{-u_k t}=e^{-t/2}/(2\sinh t)`$. Lipschitz continuity gives absolute convergence near zero, and exponential decay gives it at infinity. To justify the contour step, smooth $`j`$ by compact symmetric approximate identities. For a smooth test the usual logarithmic derivative contour shift applies, on heights avoiding the zeros using the Dirichlet zero count and partial fractions. In the limit the Fourier transforms are dominated by the original summable majorant, prime contributions stay finite, and the uniformly bounded Lipschitz constants give dominated convergence for (2). [Guinand–Weil formula for Dirichlet functions](https://www.cs.uleth.ca/~kadiri/articles/explicit_zerofree_regions_for_dirichlet_l_functions.pdf)

On the entire interval, $`P(-x)=0`$ since $`e^{-x+j_0}\lt 2`$, while $`P(x)`$ uses $`n\le e^{3.214+j_0}\lt 71`$. The nonzero pairs $`(n,c_n\log p)`$ are

```math
(4,-2\log2),\ (16,-2\log2),\ (31,8\log31),\quad
(47,-8\log47),\ (49,-8\log7),\ (64,-2\log2).
```

At 2 only the two odd-conductor characters contribute. At an odd prime the quartic powers sum to 4 if $`p^k\equiv\pm1\bmod16`$ and to zero otherwise. This accounts for prime powers and rational ramification in the completed logarithmic derivative used to recover $`U`$.

**First computation.** The kernel $`j`$ has Lipschitz constant below 14 and supremum below 9, by differentiation of the two branches of $`g`$. For $`0\lt t\lt .001`$ the absolute value of the integrand in (2) is at most 14 because $`e^{-t/2}/(2\sinh t)\le1/(2t)`$. The tail from 8 costs at most $`36e^{-12}/[(3/2)(1-e^{-16})]`$ inside the integral. On the middle range, directed rectangle sums enclose the integral with $`x`$ ranging over each complete cell of width $`.0005`$ from 1.014 to 1.040. For $`x\ge1.040\gt j_0`$, the odd-factor explicit formula simplifies to $`U=P-8\sum_{k\ge0}J(-u_k)e^{-u_kx}`$. All the added terms are positive since the kernel for $`G`$ is nonnegative and $`2\sinh(-u_ka/2)/a\lt 0`$. Keeping just $`k=0`$ gives (1) on the remaining 4,348 cells. [Integral calculation](profile.py)

**Second computation.** Put $`c=2(\sqrt2-1)^2`$ and $`D_u(v)=\sqrt2(e^{-u|v+a|}+e^{-u|v-a|})-3e^{-u|v|}`$. The transform of $`e^{-u|t|}`$ is $`2u/(u^2-z^2)`$. Multiplying it by the transform of $`g`$ and splitting into partial fractions gives

```math
(e^{-u|\cdot|}*g)(v)=\frac{2u\,g(v)+D_u(v)/c}{u^2-1/4}.
```

The $`k`$th integral in the exponential expansion of (2) therefore equals

```math
\frac{j(x)/(2u_k)+[D_{u_k}(x+a/2)-D_{u_k}(x-a/2)]/(ac)}
     {u_k^2-1/4}.
```

Its absolute value is at most $`28/u_k^2`$ directly from the original integral and the Lipschitz bound. After 256 terms the omitted contribution to $`U`$ is at most $`224[u_{256}^{-2}+(2u_{256})^{-1}]\lt .22`$ by integral comparison. A second verifier uses this series, the original tent expression for $`g`$, and a grid of 2,200 closed cells of width $`.001`$. Every lower bound after subtracting the full tail exceeds $`-9/5`$. [Exponential-series calculation](profile_series.py)

## Transfer the entire interval

For each $`\eta\gt 0`$ there exist arbitrarily large positive $`\tau`$ with $`\sup_x |U(x+\tau)-U(x)|\lt \eta`$. First make the absolute tail of the series less than $`\eta/4`$. For the finite remaining frequencies, the sequence of phase vectors at integer times lies on a compact torus. Take a convergent subsequence. Differences between sufficiently widely separated members tend to the identity along arbitrarily large time differences. The finite head then changes by less than $`\eta/2`$, uniformly in $`x`$. This argument permits every rational relation among the frequencies.

Choose $`\tau`$ arbitrarily large so that this difference is less than $`.05`$ and $`|f-U|\lt .05`$ on the translated interval. By (1), $`f(\tau+t)\gt -1.9`$ for **every** $`t\in[1.014,3.214]`$. Taking $`H=e^{\tau+1.014}/(2\sqrt2)`$, every permitted doubling with $`2y\le2e^{11/5}H`$ has its log centre in that interval translated by $`\tau`$, whose length is $`3.214-1.014=11/5`$. Thus no permitted move loses more than $`1.9735\log2`$. Oddness gives $`U\lt 9/5`$ on $`[-3.214,-1.014]`$. Translating the reflected interval yields arbitrarily large starting heights with no gain. The theorem locates arbitrarily late failures without giving their first numerical height. [Lean transfer and endpoint deductions](FactorObstruction.lean)

## Direct finite certificate

For unramified rational $`p`$ put $`u=(p/547)`$ and $`v=(p/727)`$. The common base residue degree $`d`$ is the least of $`1,2,4`$ with $`p^d\equiv\pm1\bmod16`$ and $`(uv)^d=1`$. There are $`8/d`$ prime ideals of norm $`p^d`$, each of splitting sign $`u^d`$. Their total signed log weight is $`-8u^d\log p`$. The segmented sieve marks multiples in the two classes $`16k+1`$ and $`16k+15`$ for degree one, and separately includes higher-degree ideals from primes through the square root of the largest norm. The first multiple of a sieving prime in each class is chosen at or beyond its square. [Sieve](sieve.cpp)

The logarithms use rigorous integer enclosures. For $`k=\lfloor\log_2p\rfloor`$, $`D=2^{15}`$ and $`j=\lfloor D(p/2^k-1)\rfloor`$, set $`z=Dp/[(D+j)2^k]-1`$. Then $`0\le z\lt 1/D`$ and $`\log p=k\log2+\log(1+j/D)+\log(1+z)`$. Arb certifies each lookup floor in units of $`2^{-40}`$. The rational bounds $`z-z^2/2\le\log(1+z)\le z`$ give the integer correction with outward rounding. The sums of signed $`\log p/(Nv)`$ use units of $`2^{-72}`$. Integer bin differences enclose each of the two branches of $`W`$ before conversion to Arb. Independent primality tests, Euler's criterion and direct logarithms check the bin bounds through norm 250,000.

The certificate lists outward enclosures of $`(S(2y)-S(y))/\log2`$ for 625 integer starts. The largest norm sieved is 32,467,011,360. Given $`H`$, take the **next strictly later** certified witness of each sign. Between consecutive witness starts this pair stays fixed, so the largest ratio of their endpoints to $`H`$ is attained at the left endpoint. Checking these breakpoints and both range endpoints proves the two finite bounds for every real $`H`$ in the stated ranges. At $`H=10^9`$ a gain runs from 1,016,187,975 to 2,032,375,950 and a loss from 1,714,594,397 to 3,429,188,794. [Integer table](evidence/windows.tsv) · [Recomputation and coverage](finite.py) · [Expected output](evidence/check.txt)

The explicit formula and almost-period calculation are written proofs. The two profile verifiers use distinct representations and the finite check uses actual prime ideals. The Lean kernel checks the transfer and log-window implications. The global theorem is supported by the linked external kernel runs.
