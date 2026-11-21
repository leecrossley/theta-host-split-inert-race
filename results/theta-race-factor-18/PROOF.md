# Waiting factors in the 547/727 prime race

## Two statements about the same race

Let $`K=\mathbf Q(\zeta_{16})^+`$, $`F=K(\sqrt{397669})`$ and $`E=K(\sqrt{-547},\sqrt{-727})`$. Let $`\epsilon_v=1`$ for split and $`-1`$ for inert base prime ideals. Omit the places above $`2,547,727,343057911569,343070636977`$. Put $`a=\log2`$, $`B=2\sqrt2(\sqrt2-1)^2`$, and

```math
W(u)=\max(0,\min(2-2/u,4/u-1)),\quad
S(y)=\frac{\sum_v-\epsilon_v\log(Nv)W(Nv/y)}{B\sqrt y}.
```

A gain or loss means $`S(2y)-S(y)\gt 1.9735a`$ or $`S(2y)-S(y)\lt -1.9735a`$. For a starting height $`H`$ and factor $`C`$, the permitted starts satisfy $`H\lt y\le CH/2`$.

1. For every real $`1\le H\le10^9`$, both moves finish by $`18H`$. For $`10^6\le H\le10^9`$, both finish by $`(13603772/2048491)H\lt 7H`$.
2. Assume the generalised Riemann hypothesis (GRH) for the primitive Dirichlet $`L`$-functions in the 547/727 race and its base field, including the Riemann zeta function. For arbitrarily large $`H`$ there is no loss finishing by $`2e^{11/5}H`$. For arbitrarily large $`H`$ there is no gain finishing by that factor. These may be different heights. Since $`2e^{11/5}\gt 18.05`$, neither $`18`$ nor $`18.05`$ works universally.

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

## The sharp finite deadline

Fix the field, the excluded primes and the weight in the first section. Set $`r=\sqrt2`$, $`a=\log 2`$, $`B_0=2r(r-1)^2`$ and $`D(y)=(S(2y)-S(y))/a`$. For $`t\ge0`$ call $`C`$ admissible if for every real $`1\le H\le10^9`$ there are starts $`y_+,y_-\gt H`$ with $`2y_\pm\le CH`$, $`D(y_+)\gt t`$ and $`D(y_-)\lt-t`$. Write $`C_*(t)`$ for the infimum of these factors.

Put $`m=201007/2`$ and $`\tau=-D(m)`$. The following formulas give the sharp answer for $`0\le t\le5`$. In both branches the admissible factors form the open interval $`(C_*(t),\infty)`$.

For $`0\le t\lt\tau`$ put $`g=8/(B_0a)`$ and define

```math
A=g\big[(1+2/r)\log31+(\log47+\log7)/r\big],
```

```math
B=g\big[4(1+1/r)\log31/31+(8/r)(\log47/47+\log7/49)\big].
```

Then $`C_*(t)=2((t+\sqrt{t^2+4AB})/(2B))^2`$. The controlling height is $`H=1`$.

For $`\tau\le t\le5`$ put $`p=20068`$ and $`q=214710`$. Define $`u(t)`$ as the supremum of the strict-loss starts in $`[p,m)`$ and $`v(t)`$ as the infimum of the strict-loss starts in $`(m,q]`$. Then

```math
C_*(t)=\frac{2v(t)}{u(t)}.
```

Only prime norms up to $`8q=1717680`$ enter these two endpoints. Each is a root of a quadratic in $`\sqrt y`$ on an explicitly enumerable cell. The controlling height is now $`H=u(t)`$. At $`t=\tau`$ one has $`21033.5229566\lt u\lt21033.5229567`$ and $`207428.0927825\lt v\lt207428.0927826`$. The left limit $`17.848747144\ldots`$ and the value $`19.723571102\ldots`$ enclose 18 strictly between them. Since admissibility decreases with the threshold, **18 covers every real height through one billion exactly for $`0\le t\lt\tau`$**, including all thresholds beyond the range of the displayed formulas in the exclusion.

### Between two prime norms

For a good rational prime $`\ell`$, let $`h=(\ell/547)`$, $`j=(\ell/727)`$ denote Legendre symbols. Its residue degree is the least $`f\in\{1,2,4\}`$ with $`\ell^f\equiv\pm1\bmod16`$ and $`(hj)^f=1`$. Its prime ideals together contribute the signed logarithm $`-8h^f\log\ell`$ at norm $`n=\ell^f`$. This is the fibre trace used in the original finite sieve. In $`\sqrt{y}D(y)`$ that mass is multiplied by $`(W(n/(2y))/r-W(n/y))/(B_0a)`$. The tent changes branch only at $`y=n/8,n/4,n/2,n`$. Between consecutive knots the normalized move has the form $`(c+dy)/\sqrt y`$.

In $`x=\log y`$ each piece therefore solves $`F''=F/4`$. Its value between endpoints $`s\lt x\lt w`$ is

```math
\frac{\sinh((w-x)/2)}{\sinh((w-s)/2)}F(s)
+\frac{\sinh((x-s)/2)}{\sinh((w-s)/2)}F(w).
```

The two weights are positive and sum to at most one, by the hyperbolic addition formula. Thus if both endpoints are at least $`-t`$ with $`t\ge0`$, the whole cell is at least $`-t`$. This also holds on pieces clipped at an endpoint of the search range. A strictly negative minimum can be certified by testing the knots. At an equality knot with negative value, a strict bound at each adjoining knot gives a strict bound in both open cells.

### The first loss

Before $`y=31/4`$ the knots, including the left boundary, are $`1,31/8,47/8,49/8,31/4`$. All five values are nonnegative. At the first support boundary every weight vanishes. The knot bound covers the intervening starts. On $`[31/4,79/8]`$ the contributing norms in the open cell are $`31,47,49`$, with signed masses $`8\log31,-8\log47,-8\log7`$. Inserting the two tent branches gives exactly $`D(y)=(A-By)/\sqrt y`$. It decreases strictly, from a positive value to a value below $`-5`$. The first boundary of a strict loss at height one is $`y_0(t)=((t+\sqrt{t^2+4AB})/(2B))^2`$. A gain is already available at $`y=5`$ with $`D(5)\gt5`$. The bounds $`2y_0(0)\gt17`$ and $`2y_0(5)\lt18.1`$ hold directly from the coefficients.

### The loss that disappears

The finite enumeration uses integer coordinates $`k=8y`$. Among all 26,960 knots from $`k=168286`$ to $`k=1659407`$, the unique minimum is at $`k=804028`$, corresponding to $`y=m`$. Every other knot lies more than $`0.000006`$ above that minimum. The cell bound proves $`D(y)\ge-\tau`$ throughout this closed central interval, with equality only at $`m`$. Both endpoint values at $`p,q`$ are strictly below $`-5`$. In particular, $`p\lt u(t)\lt21100`$ and $`207000\lt v(t)\lt q`$ for $`\tau\le t\le5`$. The bounds follow from the certified central interval and continuity at $`p,q`$. At and above $`\tau`$ the whole closed interval $`[u(t),v(t)]`$ is free of strict losses. Strict losses approach its endpoints from the outside. The endpoints can be found by inspecting consecutive knots outwards from the centre and solving $`c+dy=-t\sqrt y`$ at the boundary of the first strict-loss component on each side.

To prove the global upper bounds, the 625 original enclosures supply witnesses with $`D\gt5`$ or $`D\lt-5`$. For real $`2\le H\le10^9`$ the next strictly later gain always finishes by $`41H/4`$. The next loss finishes by $`10H`$ outside $`[p,q)`$. For a sorted witness list, it suffices to check the left boundary, the upper endpoint and each witness height in the range. On each gap the next available witness is constant and its ratio to $`H`$ decreases. Including the witness height itself enforces the strict start condition. [Original grid](evidence/windows.tsv) · [Coverage calculation](sharp_deadline.py)

If $`t\lt\tau`$, the additional strict loss at $`m`$ splits the exceptional gap. Both signs then finish by $`41H/4`$ for all $`2\le H\le10^9`$. For $`1\le H\lt2`$ use the gain at 5 and losses just after $`y_0(t)`$. Every factor strictly greater than $`2y_0(t)`$ works. At $`H=1`$ the endpoint is only an equality and every earlier start misses the loss. This proves the first branch, including nonattainment.

Suppose now that $`\tau\le t\le5`$. The factor $`2v/u`$ is greater than $`2(207000)/21100\gt19.6`$, which covers the early heights and both grid bounds outside the gap. For $`p\le H\lt u`$ take a strict loss approaching $`u`$ from below. Its endpoint ratio is at most $`2(21100)/p`$. For $`u\le H\le v`$ take a strict loss just above $`v`$. Every $`C\gt2v/u`$ gives $`CH\gt2v`$, so this start fits. For $`v\lt H\lt q`$ use the loss at $`q`$, with ratio at most $`2q/207000`$. Both outside ratios are below three. The gain still comes from the grid. At $`H=u`$ a factor of $`2v/u`$ permits starts only in $`(u,v]`$, throughout which the loss threshold is missed. This proves the second branch.

### Numerical and formal checks

One implementation sieves rational primes by Eratosthenes, finds their signs from residue tables and accumulates the changes of $`c,d`$ at each knot. An independent implementation tests primality with FLINT, applies Euler's criterion and forms the two original tent sums from prefixes of signed $`\log\ell`$ and $`(\log\ell)/n`$. Both run with 176-bit Arb arithmetic. They agree in overlapping intervals of radius below $`10^{-30}`$ at every one of 31,776 knots and endpoints through $`q`$. Comparison against the independently regenerated integer sieve checks a further 305 common grid values. The verifier brackets the boundary roots at $`\tau,4,4.5,5`$, confirms the signs of the derivatives on the four pairs of crossing cells, and compares all witness ratios as exact fractions. Four local Lean lemmas check the interpolation inequality, the failed endpoint, the strict upper bound immediately beyond a gap and the touching bridge. The prime enumeration and the analytic cell representation are covered by the programs and the foregoing argument. [Independent evaluations](prime_knots.py) · [Root and coverage checks](sharp_deadline.py) · [Lean deductions](FactorObstruction.lean) · [Replay transcript](evidence/check.txt)

## All failing starts from one million to one billion

Retain the field, omitted places, tent and score from the first section. Put $`t=1.9735`$, $`D(y)=(S(2y)-S(y))/\log2`$ and fix $`10^6\le H\le10^9`$. For a sign $`\sigma\in\{+1,-1\}`$ an eligible move requires $`H\lt y\le CH/2`$ and $`\sigma D(y)\gt t`$. The following six real numbers give the complete answer for $`C\ge6`$. Each is the unique root of $`\sigma D(y)=t`$ in the specified prime-knot cell. The integers in square brackets are bounds on $`8y`$. The decimals below identify the roots to seven places by truncation.

| $`i`$ | $`\sigma_i`$ | Cell for $`8a_i`$ | $`a_i`$ | Cell for $`8b_i`$ | $`b_i`$ |
|---|---:|---|---:|---|---:|
| 1 | $`+1`$ | $`[50308994,50309071]`$ | $`6288626.8228400\ldots`$ | $`[162452254,162452386]`$ | $`20306533.5953377\ldots`$ |
| 2 | $`-1`$ | $`[160619068,160619119]`$ | $`20077385.9489463\ldots`$ | $`[488563106,488563336]`$ | $`61070400.3947486\ldots`$ |
| 3 | $`-1`$ | $`[4499851921,4499852033]`$ | $`562481498.8459634\ldots`$ | $`[13610271236,13610271343]`$ | $`1701283915.0779512\ldots`$ |

**Theorem.** For each real $`C\ge6`$ the starting heights missing the required sign are precisely the intervals

```math
\left[a_i,\frac{2b_i}{C}\right]\qquad (i=1,2,3),
```

taken only for the sign in that row and interpreted as empty when the upper endpoint is smaller. These intervals are disjoint whenever nonempty. At $`C=c_i=2b_i/a_i`$ the corresponding interval is the single height $`a_i`$. The critical factors in decreasing order are $`6.458177331046\ldots,6.083501163950\ldots,6.049208439987\ldots`$. For real starts both signs are covered exactly for $`C\gt c_1`$. The strictness of the move matters at equality.

At $`C=6`$ the failing integers in row $`i`$ run from $`\lceil a_i\rceil`$ through $`\lfloor b_i/3\rfloor`$. This gives the three ranges in the README and exactly 5,372,773 failures among the 999,000,001 permitted integer starts. With only integer starts the critical factors are $`2b_i/\lceil a_i\rceil`$. For general $`C\ge6`$ the total logarithmic length of the failing real starts is $`\sum_{i=1}^3\max(0,\log(c_i/C))`$. This is a finite integral against $`dH/H`$ over the stated range.

### Reducing the full range to three intervals

Sort the original certified starts of each sign and look up the first strictly later witness after each breakpoint. Outside the following three gaps the next same-sign endpoint is strictly below $`6H`$. The test includes both domain endpoints and every witness height in the domain. On a gap the endpoint is constant and the ratio to $`H`$ falls. [Original witnesses](evidence/windows.tsv) · [Exact comparison](late_roots.py)

| Sign | Left witness $`p`$ | Right witness $`q`$ |
|---|---:|---:|
| $`+1`$ | 6,145,473 | 20,405,658 |
| $`-1`$ | 19,787,304 | 61,780,789 |
| $`-1`$ | 549,154,839 | 1,714,594,397 |

For each row the finite calculation proves $`p\lt a_i\lt b_i\lt q`$, $`\sigma_iD(y)\le t`$ on $`[a_i,b_i]`$ and strict qualifying values approaching $`a_i`$ from the left and $`b_i`$ from the right. It also proves $`2a_i/p\lt3`$ and $`2q/b_i\lt3`$. Here and below intersect each range of $`H`$ with $`[10^6,10^9]`$. If $`p\le H\lt a_i`$, a strict witness near $`a_i`$ fits before $`3H`$. If $`b_i\lt H\lt q`$ use $`q`$. If $`a_i\le H\le b_i`$, a strict witness fits exactly when $`CH\gt2b_i`$. At equality every permitted start is still in the closed gap. Thus the failures are exactly the displayed closed intervals. The verifier checks $`2b_i/6\lt a_{i+1}`$ for adjacent rows and the inclusion of all three intervals at 6 in the permitted range. This proves disjointness for every higher factor and justifies adding their integer counts and logarithmic lengths.

### Certifying the interiors without sampling assumptions

Write $`G(y)=\sqrt yD(y)`$. Between the prime knots $`G`$ is affine. Let $`[\alpha,\beta]`$ be a trial interval with length $`h`$. The secant interpolation error of a unit hinge $`(y-k)_+`$ inside the interval has magnitude at most $`(k-\alpha)(\beta-k)/h\le h/4`$. On each side of $`k`$ the interpolation error is linear and attains its maximum magnitude at $`k`$. A piecewise affine function is its initial affine function plus the sum of its slope jumps times unit hinges, so its absolute secant error is at most $`h/4`$ times the sum of the absolute jumps inside the interval.

For a good norm $`n=\ell^f`$, its coefficient in $`G`$ is the signed $`8\log\ell/(B_0\log2)`$ from the fibre trace. At a knot $`y=n/j`$ the absolute change of slope has the additional factor $`M_j/n`$. In the order $`j=8,4,2,1`$ these constants are

```math
M_8=8/\sqrt2,\qquad M_4=12/\sqrt2+4,\qquad
M_2=4/\sqrt2+6,\qquad M_1=2.
```

These follow by differentiating the two linear tent branches on either side of each knot. There are at most $`jh+1`$ integer norms in $`[j\alpha,j\beta]`$. Each good norm comes from a unique rational prime after the fibre has been summed. Since $`n\ge j\alpha`$ and $`\ell\le8\beta`$, the following bounds the absolute secant error of either signed $`G`$.

```math
E(\alpha,\beta)=
\frac{h}{4}\frac{8\log(8\beta)}{B_0\log2}
\sum_{j=1,2,4,8}M_j\frac{jh+1}{j\alpha}.
```

Both endpoint margins $`(t-\sigma D(v))\sqrt v`$ exceed $`E`$ on every accepted cell. Linear interpolation therefore keeps the entire signed graph below the secant through $`t\sqrt\alpha,t\sqrt\beta`$. The concavity of $`t\sqrt y`$ puts that secant below the threshold curve. This proves the strict exclusion $`\sigma D(y)\lt t`$ throughout each interior cell. The integer endpoint partition closes every interval and is reconstructed from three regular meshes with 1,213 further abscissae. The printed totals count cells separately for each sign. [Certified mesh](evidence/late-subdivisions.txt) · [Error and coverage checks](late_slabs.py)

### Resolving the six crossings

The knot search in the six boundary neighborhoods uses $`k=8y`$ and tests every possible norm in each of the four relevant prime-power ranges. For the three interiors the integer endpoints in ordinary $`y`$ units are $`[6288755,20306405]`$, $`[20077514,61070272]`$ and $`[562485616,1701279810]`$. Each boundary search reaches the adjacent interior endpoint. Two implementations enclose the normalized move at all 1,530 resulting boundary knots and endpoints.

The first forms 160-bit Arb logarithms of the primes in two independently sieved residue progressions. It uses signed prefix sums of $`\log\ell`$ and $`(\log\ell)/n`$ for the exact tent branches. Higher degrees are determined separately using primality tests and Euler’s criterion. The second implementation uses the integer segmented sieve with a finer fixed-point log table. It encloses $`z,z^2,z^3,z^4`$ by directed integer arithmetic and bounds $`\log(1+z)`$ between its fourth and third Taylor sums for $`0\le z\le2^{-15}`$. The log units are $`2^{-58}`$ and the reciprocal units $`2^{-96}`$. The earlier sieve precision and its stored billion-height table are also replayed. [Independent prime sums](near_roots.py) · [Integer logarithms](sieve.cpp) · [Rational-height scores](band_values.py)

On each candidate crossing cell, $`\sqrt y(\sigma D(y)-t)`$ is an affine function minus $`t\sqrt y`$. The derivative has one strict sign throughout the reported cell. The two methods enclose its root independently, starting from the same prime-knot cells and crossing code. Their root intervals overlap and their factor intervals agree within radius $`10^{-11}`$ for the integer implementation. The high-precision method supplies the decimal bounds in the table. All intervening boundary cells are below threshold by the endpoint interpolation inequality from the preceding finite theorem. On the opposite sides of the roots continuity gives strictly qualifying values arbitrarily close. These checks together with the interior certificate establish all the hypotheses in the three-gap reduction. [Endpoint, integer-count and coverage verifier](late_roots.py) · [Verification output](evidence/check.txt)

Three further local Lean lemmas check the one-hinge ceiling, interpolation with the certified error margin and the exact condition $`CH\le2b_i`$ for failure inside a gap. The prime enumeration, knot locations and numerical inequalities have the separate reproducible checks above. [Lean deductions](FactorObstruction.lean)
