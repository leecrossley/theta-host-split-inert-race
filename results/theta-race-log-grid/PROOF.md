# Two whole intervals of winning starts on every logarithmic grid

## The theorem

On $`K=\mathbf Q(\zeta_{16})^+`$, $`F=K(\sqrt{397669})`$ and $`E=K(\sqrt{-547},\sqrt{-727})`$, take the weighted inert-minus-split prime-ideal score of width two. At good $`v`$ write $`\epsilon_v=1`$ for split and $`-1`$ for inert. These are the signs from the first-step theta traces. Omit primes above $`2,547,727,343057911569,343070636977`$. Put $`a=\log2`$, $`r_0=\sqrt2`$, $`c=2(r_0-1)^2`$, $`B=r_0c`$ and

```math
W(u)=\max(0,\min(2-2/u,4/u-1)),\quad
S(y)=\frac1{B\sqrt y}\sum_v-\epsilon_v\log(Nv)\,W(Nv/y),\quad
f(x)=\frac{S(2y)-S(y)}a,\quad x=\log(2r_0y).
```

Assume the generalised Riemann hypothesis (GRH) for the primitive Dirichlet functions in the 547/727 race and its base field, including the Riemann zeta function. Let $`C=10^{18}`$, $`\eta=10^{-6}`$ and $`d=\eta/2`$. For every real $`H\ge1`$ there are two closed intervals of **log starting heights** $`\log y`$, each of length $`\eta`$, inside $`(\log H,\log(CH/2))`$. At every point of the first the doubling gains more than $`1.9735a`$. At every point of the second it loses more than $`1.9735a`$. Alpöge's announced Dirichlet GRH theorem supplies this input. The range $`1\le H\le10^9`$ does not use it. [Announcement](https://alpo.ge/remark.pdf) [Independent formal checks](https://math.mit.edu/~drew/grh/) [Score and trace identity](https://github.com/leecrossley/theta-race-short-shifts/blob/483c31f2df2132ed2165620ea55e20e969228bfb/RaceTrace.lean#L83) [Pointwise deadline theorem](https://github.com/leecrossley/theta-race-filtered-deadline/blob/2400dc5ac2a7969d41e1a4bde1eb44dd6c63b7a2/PROOF.md)

In particular, on **every** geometric grid $`y=\exp(v+nj)`$, with arbitrary real offset $`v`$, integer $`n`$ and $`0\lt j\le\eta`$, both moves have eligible grid starts after every $`H`$. For $`j=\eta`$ there are at most 40,753,385 eligible grid points. Indeed any closed interval $`[t,t+j]`$ contains the grid point indexed by $`\lceil(t-v)/j\rceil`$. The log length of the eligible range is $`\log(C/2)`$, so its point count is at most $`1+\lfloor\log(C/2)/j\rfloor`$. [Grid lemma in Lean](LogGrid.lean)

## A quantitative modulus from the even zero moment

For the multiset of ordinates of the eight odd primitive factors of $`\zeta_E/\zeta_F`$, put

```math
G(t)=\frac{3-2r_0\cos(at)}{c(t^2+1/4)},\quad
\alpha(t)=\frac{2G(t)\sin(at/2)}a,\quad
V(x)=\sum_\nu G(\nu)e^{i\nu x},\quad
U(x)=\sum_\nu i\alpha(\nu)e^{i\nu x}.
```

All $`G`$ are positive, the ordinates pair under conjugation, and both series converge absolutely and uniformly. GRH makes every ordinate real. The predecessor bounds the difference between $`f`$ and $`U`$ at each $`x\ge55`$. In its notation $`|f(x)-U(x)|\le E(x)`$, with $`E`$ decreasing, $`E(58)\lt6`$ and $`E(100)\lt0.007`$. The earliest linked derivation states the bound for $`x\ge1000`$. It already holds for $`x\gt53+3a/2`$, because every omitted norm is below $`e^{53}`$ and the other terms do not depend on the window. This release needs only $`x\gt58`$. The base-field term of $`E`$ uses GRH for the base field, including zeta. [Explicit remainder and arithmetic inputs](https://github.com/leecrossley/theta-race-filtered-deadline/blob/2400dc5ac2a7969d41e1a4bde1eb44dd6c63b7a2/PROOF.md) [Independent primitive-theorem verification](https://math.mit.edu/~drew/grh/)

The Fourier coefficients of $`|\sin(s/2)|`$ give an identity with nonnegative terms

```math
|e^{is}-1|=\frac8\pi\sum_{k\ge1}\frac{1-\cos(ks)}{4k^2-1}.
```

Integrate the even cosine coefficients on $`[0,\pi]`$ by the product-to-sum identity, or compute the constant term $`4/\pi`$ for the chord and the cosine coefficients $`-8/[\pi(4k^2-1)]`$. The series is absolutely convergent and $`\sum_{k\ge1}(4k^2-1)^{-1}=1/2`$ by telescoping. Since $`|\alpha|\le2G/a`$, this gives the uniform majorant

```math
|U(x+h)-U(x)|\le\frac{16}{a\pi}
   \sum_{k\ge1}\frac{V(0)-V(kh)}{4k^2-1}.
```

Positivity permits interchanging the two sums. It remains to bound the increments of the even companion, including when $`h`$ tends to zero.

Write $`g(t)=4\sinh(\max(0,a-|t|)/2)/B`$, $`K_g=2\cosh(a/2)/B`$ and $`D=\sum_\chi\log(q_\chi/\pi)+8\psi(3/4)`$. Here $`K_g`$ bounds the Lipschitz constant of $`g`$, and the conductors are $`547,4376,8752,8752,727,5816,11632,11632`$. For $`0\le s\le1/4`$ the compact explicit formula reduces to

```math
V(s)=Dg(s)-8\int_0^\infty\frac{e^{-t/2}}{2\sinh t}
   [g(s-t)+g(s+t)-2g(s)]\,dt .
```

Only the integer 2 can occur in the prime sum since $`a+1/4\lt\log3`$. Its two odd-conductor characters have opposite values and the even-conductor factors vanish at 2. The same cancellation applies to the reflected test. [Even explicit formula and prime terms](https://github.com/leecrossley/theta-race-distant-interval/blob/527fffe7b7808457b47a292144fea96b556158a8/PROOF.md)

Subtract the integral at zero. The absolute difference between the two bracketed second differences is at most $`4K_g\min(t,s)`$. Split the integral at $`s`$ and 1. On $`0\lt t\le1`$ its kernel is at most $`1/(2t)`$. Above 1 the kernel is at most $`e^{-3t/2}/(1-e^{-2})`$. Including the conductor term gives

```math
0\le V(0)-V(s)\le s(309+70\log(1/s))\quad(0\lt s\le1/4).
```

Indeed the constant part is bounded by $`K_g[|D|+16+32e^{-3/2}/((3/2)(1-e^{-2}))]\lt309`$ and the logarithmic coefficient by $`16K_g\lt70`$. At zero the increment is zero. For all arguments the bound $`V(0)-V(s)\le2V(0)\lt350`$ holds.

For $`|h|\le d`$ take $`J=500000`$ so that $`Jd=1/4`$. The local bound increases in $`s`$ on this range. With $`L=\log(1/d)`$ and $`K=\log J`$, the terms through $`J`$ in the chord sum before its prefactor are at most

```math
d\left[\frac{309+70L}{3}
 +\frac4{15}\big((309+70L)K-35K^2\big)\right].
```

The first term treats $`k=1`$. For $`k\ge2`$, $`k/(4k^2-1)\le4/(15k)`$. The decreasing function $`(309+70L-70\log t)/t`$ permits bounding its sum from 2 through $`J`$ by its integral from 1 through $`J`$. The tail is at most $`350/[2(2J+1)]`$. After multiplication by $`16/(a\pi)`$ the total is below 0.0141. [Directed integral and harmonic calculation](modulus.py)

An independent bound uses the closed gamma series

```math
T(s)=\sum_{j\ge0}\frac{e^{-(2j+3/2)s}}{(2j+3/2)^2-1/4}
 =\frac{(1+z)\log(1+z)+(1-z)\log(1-z)}{2\sqrt z},\quad
 z=e^{-s},\quad T(0)=a .
```

Set $`D_*=\sum_\chi\log(q_\chi/\pi)-8(\gamma_E+a)\gt0`$. For $`0\le t\le s\le1/4`$, the closed moment formula and monotonicity of $`T`$ give

```math
V(0)-V(t)\le D_*K_gs+\frac8c\left[
 r_0\{T(a-s)-T(a)\}+3\{a-T(s)\}\right].
```

Indeed $`T(a+t)-T(a)\le0`$ in the exact difference. Group the integers through $`J`$ into 344 successive blocks $`[l,r]`$, with $`r=\min(J,\max(l,\lfloor1.03l\rfloor))`$, starting at 1. In each use $`s=rd`$ and the exact reciprocal mass $`[(2l-1)^{-1}-(2r+1)^{-1}]/2`$. Bound the tail using the closed value of $`V(0)`$. This proves the stronger uniform profile increment of 0.0085. [Closed-moment derivation](https://github.com/leecrossley/theta-race-distant-interval/blob/527fffe7b7808457b47a292144fea96b556158a8/PROOF.md) · [Independent block calculation](closed_moment.py)

## Every height and both signs

For $`H\ge e^{100}`$ use the positive filter in the deadline theorem. Its support has radius $`R=155/8`$. The carrier ordinate $`\gamma\gt3.287`$ takes both opposite extremal phases over a period $`2\pi/\gamma`$. Buffer the support by $`d`$ on each side. There is still an open interval of centres longer than one period because $`2R+2d+2\pi/3.287\lt\log(C/2)`$. For each sign the average of the signed **profile** $`U`$ exceeds $`4.917-2.85-1.44/(48\cdot3/32)^2`$. Some point in the buffered support therefore exceeds that number. Everywhere within $`d`$ of that point the signed prime score remains above the threshold, since

```math
4.917-2.85-\frac{1.44}{(48\cdot3/32)^2}-0.0141-0.007\gt1.9735 .
```

The second proof uses 1.52 for the omitted spectral mass and 0.0085 for the profile increment, still clearing the same threshold. The slightly enlarged supports lie above 100 and inside the eligible range. Log centres and log starting heights differ by the constant $`3a/2`$, so interval lengths agree. [Filter margins and certificates](https://github.com/leecrossley/theta-race-filtered-deadline/blob/2400dc5ac2a7969d41e1a4bde1eb44dd6c63b7a2/PROOF.md) · [Exact inequalities in Lean](LogGrid.lean)

For $`10^9\le H\le e^{100}`$ use the three pairs of log abscissae $`(60.50,58.40)`$, $`(83.66,88.00)`$ and $`(109.10,111.44)`$ from the deadline theorem. The normalized signed changes exceed 20 there, every enlarged interval stays above 58 and all of them fit their published ranges after adding the buffer. One pointwise error at the centre, one at a nearby point, and the uniform modulus cost less than $`12+0.0141`$. [Height comparisons](heights.py)

For $`1\le H\le10^9`$ use the finite sieve at $`y_+=1016187975`$ and $`y_-=2893002095`$, with a normalized gain of at least 33.1637 and a loss of at least 42.7045 in magnitude. On a smooth branch the derivative with respect to $`\log z`$ of $`z^{-1/2} W(Nv/z)`$ is at most $`(5/2)z^{-1/2}`$ in absolute value, since $`|uW'(u)|\le2`$ and $`0\le W\le1`$. The weight is continuous at its corners. Over a bounded window there are finitely many terms, so this derivative bound also gives an interval Lipschitz bound. The sum of $`\log Nv`$ above $`p`$ is at most $`8\log p`$.

If $`4ze^d\lt2^k`$, the log-derivative bound for the contribution of $`S`$ at that scale on the whole perturbation interval is at most $`20\cdot2^{k+1}\log2/(B\sqrt{ze^{-d}})`$. Indeed the sum of $`\log p`$ up to $`2^k`$ is at most $`2^{k+1}\log2`$. At each dyadic level the primes in its upper half divide the central binomial coefficient, whose logarithm is at most $`2^j\log2`$ at level $`j`$. Sum over the levels. Use the bound for $`z=y_\pm`$ and $`2y_\pm`$, add, divide by $`a`$ and multiply by $`d`$. The costs are below 14 and 32 respectively, leaving both signs strictly above the target magnitude. Their perturbed starts exceed a billion and their finishes are below $`10^{18}`$ even at $`H=1`$. [Pinned finite sieve](https://github.com/leecrossley/theta-race-factor-18/blob/31ea5d495ca29f7f17869d1a3f321000266e4867/PROOF.md) · [Finite rows](evidence/finite-windows.tsv) · [Verification transcript](evidence/check.txt)
