# Exceptional powers of two in the 547/727 prime race

## Statement

For $`K=\mathbf Q(\zeta_{16})^+`$, $`F=K(\sqrt{397669})`$ and $`E=K(\sqrt{-547},\sqrt{-727})`$, put $`\epsilon_v=1`$ for a split prime ideal of $`F`$ and $`-1`$ for an inert one in $`E/F`$. Omit the primes above $`2,547,727,343057911569,343070636977`$. At every good prime the normalized first-step theta trace jumps recover this sign. Set $`a=\log2`$, $`r_0=\sqrt2`$, $`c=2(r_0-1)^2`$ and $`B=r_0c`$. Fix the score width $`b=2`$. [Field and trace identity](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-short-shifts/RaceTrace.lean#L83)

```math
W(u)=\max(0,\min(2-2/u,4/u-1)),\qquad
S(y)=\frac1{B\sqrt y}\sum_v-\epsilon_v\log(Nv)\,W(Nv/y),\qquad
D_n(y)=S(2^ny)-S(y).
```

Let $`\tau=1.9735\log2`$, $`N=10^9`$, $`M=2^{90000}`$ and $`Y_*=e^{1000}/(2\sqrt2)`$. Define the exceptional exponents and a more restrictive subset by

```math
\mathcal Q=\{n\ge1:\limsup_{y\to\infty}|D_n(y)|\le\tau\},\qquad
\mathcal P=\{n\ge1:\ |D_n(y)|\lt4/3\ \text{for every }y\ge Y_*\}.
```

Here and below $`n`$ is an integer. Since $`4/3\lt\tau`$, one has $`\mathcal P\subseteq\mathcal Q`$.

**Theorem.** Assume the generalised Riemann hypothesis (GRH) for the primitive Dirichlet functions in the 547/727 race and its base field, including the Riemann zeta function.

1. For each $`1\le n\le N`$ there are arbitrarily high starts with $`D_n(y)\gt14/5`$ and arbitrarily high starts with $`D_n(y)\lt-14/5`$.
2. Any two elements of $`\mathcal Q`$ differ by more than $`N`$ if they are distinct. Every block of $`N+1`$ consecutive positive exponents therefore contains at most one exceptional exponent. At every exponent outside $`\mathcal Q`$ both signs exceed $`\tau`$ in magnitude arbitrarily far out.
3. For every positive integer $`m`$ there is $`1\le j\le M`$ with $`mj\in\mathcal P`$. The set of such multipliers also satisfies

```math
\liminf_{T\to\infty}\ \inf_{A\in\mathbf Z_{\ge0}}
  \frac{\#\{j:A\lt j\le A+T,\ mj\in\mathcal P\}}{T}\ge 2^{-90000},
```

Here $`T`$ runs through the positive integers. In particular the first exponent in $`\mathcal P`$ lies strictly above a billion and at most $`2^{90000}`$.

Thus the lower uniform density of $`\mathcal P`$ is positive and the upper uniform density of $`\mathcal Q`$ is at most $`1/1000000001`$. All heights in part 3 use the same $`Y_*`$, however large the exponent. These bounds leave the first exceptional exponent unidentified. All three parts use the hypothesis, part 1 included, because the supremum $`L_n`$ below is defined through the limiting profile. Alpöge's announced Dirichlet GRH theorem supplies this input. [Announcement](https://alpo.ge/remark.pdf) [Independent formal checks](https://math.mit.edu/~drew/grh/) The arithmetic on the certified zero brackets is unconditional.

## Telescope the dyadic profile

Write $`s=\log(2y)`$, $`x=s+a/2`$ and $`f(x)=D_1(y)/a`$. For the multiset of all ordinates of the eight odd primitive factors of $`\zeta_E/\zeta_F`$, put

```math
G(t)=\frac{3-2r_0\cos(at)}{c(t^2+1/4)},\qquad
V(s)=\sum_\nu G(\nu)e^{i\nu s},\qquad
U(x)=\sum_\nu \frac{2i}a G(\nu)\sin(a\nu/2)e^{i\nu x}.
```

Under GRH the ordinates are real. With multiplicities the family is closed under conjugation. The weights $`G`$ are strictly positive, and both series converge absolutely and uniformly. In particular $`V`$ is real and even. The coefficient identity gives $`aU(x)=V(s+a)-V(s)`$, including any central zeros. The pointwise prime-only estimate is $`|f(x)-U(x)|\le E(x)`$ for $`x\ge1000`$, with

```math
E(x)=2\cdot10^6e^{-x/4}+
 128e^{-(x-2)/2}/(1-e^{-(x-2)})+22000e^{-3x/20}.
```

The constants include the contributions of the even base factors, trivial zeros and higher prime powers. The base-factor contribution uses GRH for the base field, including zeta. The finite excluded primes have left the support. [Dyadic explicit formula, character family and prime-only error](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-factor-400/PROOF.md) [Independent upstream verification](https://math.mit.edu/~drew/grh/)

Let $`P_n(s)=V(s+na)-V(s)`$. Telescoping the actual scores and then the Fourier series gives

```math
|D_n(y)-P_n(s)|\le a\sum_{k=0}^{n-1} E(x+ka)
 \le a\sum_{k\ge0}E(x+ka)\lt10^{-19}\qquad(x\ge1000).
```

For the last inequality, retain the denominator $`1-e^{-(x-2)}`$ at the initial point. The three geometric ratios are $`e^{-a/4},e^{-a/2},e^{-3a/20}`$. Substitution at 1,000 gives the stated bound. Its right-hand majorant tends to zero as $`x\to\infty`$, uniformly in $`n`$. [Arithmetic check](power_phases.py) · [Finite telescoping identity in Lean](DyadicPowers.lean)

For fixed $`n`$, write $`L_n=\sup_{s\in\mathbf R}|P_n(s)|`$. Centred at $`s+na/2`$ the difference is odd. Either sign therefore attains values arbitrarily close to $`L_n`$ in magnitude. These values recur arbitrarily far to the right. Indeed truncate the absolutely convergent series uniformly, then use compactness of its finite phase torus to find unbounded positive translations changing every retained coordinate by arbitrarily little. Infinite pigeonholing in a finite cover of the torus yields pairs of integers arbitrarily far apart in the same small box. The discarded mass bounds the remaining translation error. Together with the vanishing prime error this proves

```math
\limsup_{y\to\infty}D_n(y)=L_n,\qquad
\liminf_{y\to\infty}D_n(y)=-L_n.
```

## Why the first billion exponents oscillate

At positive ordinate $`\nu`$ the Fourier coefficient of $`P_n`$ is $`m_\nu G(\nu)(e^{ina\nu}-1)`$, where $`m_\nu\ge1`$ is its multiplicity. The negative ordinate gives its complex conjugate. If $`\|\nu na/(2\pi)\|_{\mathbf R/\mathbf Z}\ge1/6`$ and $`G(\nu)\gt7/5`$, the coefficient has magnitude greater than $`7/5`$. The full coefficient has a common phase for every occurrence of the same ordinate.

Take two **distinct** such positive ordinates. Averaging $`|P_n|^2`$ in $`s`$ kills cross terms between distinct real frequencies. Multiplication and integration of the series are justified by absolute uniform convergence, or by taking finite sums and passing to the limit. The two positive ordinates and their negatives give four distinct terms. Thus $`L_n^2`$ is at least their sum of squared coefficient magnitudes, which exceeds $`4(7/5)^2`$. Hence $`L_n\gt14/5`$.

The 6,000 sign-bracket table comes from the one-way window certificate. Within it there are 25 low, mutually disjoint positive brackets with a strict lower bound $`G(\nu)\gt7/5`$. For each, a directed calculation produces a binary phase numerator $`p`$ such that

```math
10^9\left|\frac{a\nu}{2\pi}-\frac{p}{2^{64}}\right|\lt1/1000 .
```

For each of the first billion exponents the binary scan finds two brackets with $`np\bmod2^{64}`$ between $`11\cdot2^{58}`$ and $`2^{64}-11\cdot2^{58}`$, inclusively. Their actual distances from the nearest integer exceed $`11/64-1/1000\gt1/6`$. The multiplications are unsigned 64-bit arithmetic, so machine wraparound is exactly the required modulus. [Bracket table](evidence/zero-brackets.tsv) · [Certified phase numerators](evidence/low-phases.tsv) · [Exhaustive scan](scan_powers.c)

An independent check uses the half-angle expression for $`G`$, calculates numerators to the nearest unit of $`10^{-16}`$ in turns and encloses the error afresh. It filters blocks of 8,192 candidate exponents with exact 128-bit products modulo $`10^{16}`$, using the inner arc $`[0.17,0.83]`$, in reverse bracket order. It verifies two witnesses for every exponent. Here $`0.17-0.001\gt1/6`$. Both scans check disjointness of the contributing intervals through their associated coefficient verifiers. [Independent coefficient calculation](independent.py) · [Decimal sieve](decimal_check.c)

For positive integers $`n\gt m`$, translation invariance of the supremum norm and telescoping imply $`L_{n-m}\le L_n+L_m`$. Since $`\tau\lt7/5`$, two exceptional exponents at a distance at most $`N`$ would yield $`L_{n-m}\le2\tau\lt14/5`$, contradicting the finite certificate. This proves the gap and the upper density bound, including translated blocks. The same argument allows zero as an additional exceptional index for the purpose of excluding the first billion.

## Why some arbitrarily large powers suppress every move

Choose one positive ordinate from each of the 6,000 brackets with denominator $`10^{14}`$ and its negative partner. For their total even weight $`A`$ and the remaining weight $`R`$ the compact-test moment gives $`A\lt175`$ and $`0\le R\lt0.647`$. In detail, put $`g_0=4\sinh(a/2)/(r_0c)`$, $`C_*=\sum_\chi\log(q_\chi/\pi)-8(\gamma_E+a)`$ and

```math
T(a)=\frac{(1+1/2)\log(1+1/2)+(1-1/2)\log(1-1/2)}{2\sqrt{1/2}}.
```

The full mass is $`V(0)=C_*g_0-8(2r_0 T(a)-3a)/c`$. Subtract the directed enclosures of $`2G(\nu)`$ over the chosen brackets. Under GRH any repeated or unlisted zero stays in the positive remainder. As a second calculation, the original digamma integral at the origin is integrated in 65,536 directed rectangles per dyadic piece from $`a/2^{20}`$ to $`32a`$. The pieces below and above this range are bounded by the published Lipschitz and exponential estimates, with endpoint $`a/2^{20}`$ in the former. This gives $`R\lt0.649`$, independently sufficient for the claimed bound. [Compact moment and integral estimates](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-factor-400/PROOF.md) [Numerical transcript](evidence/check.txt)

Partition each true coordinate $`a\nu/(2\pi)`$ modulo one into $`q=32768`$ equal arcs of the circle when testing its integer multiples. If every coordinate of the multiple $`n`$ is within $`1/q`$ of an integer, the selected part contributes at most $`2\pi A/q`$ to $`|P_n(s)|`$, for every $`s`$. The remaining part contributes at most $`2R`$. Even using the weaker integral bound,

```math
|D_n(y)|\le 2\pi A/32768+2R+10^{-19}\lt4/3\qquad(y\ge Y_*).
```

For each $`m\ge1`$ the phase vectors of $`0,m,2m,\ldots,Mm`$ occupy $`q^{6000}=M`$ boxes. Two share a box and the positive difference $`mj`$ has every coordinate within $`1/q`$ of zero. This proves the effective bound for a multiple of every prescribed exponent.

For density use $`k(v)=\max(1-q\|v\|,0)`$ on the circle, with distance to the nearest integer. Its Fourier coefficients are nonnegative, absolutely summable and total one. Its constant coefficient is $`1/q`$, as follows by convolution of a length-$`1/q`$ interval with its reflection, divided by its length. The product of these tents at $`jma\nu/(2\pi)`$ is bounded by the indicator of $`mj\in\mathcal P`$ for positive integers $`j`$. The constant multi-index contributes $`q^{-6000}`$. In a discrete average all characters whose total phase is an integer contribute their nonnegative coefficients. All others average to zero uniformly in the first index by the finite geometric-series bound $`2/(T|1-e^{2\pi i\theta}|)`$ for a nonintegral phase $`\theta`$. Absolute summability permits passage to the limit. No independence assumption on the zero ordinates is needed. [Continuous version of the tent argument](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-paired-density/PROOF.md)
