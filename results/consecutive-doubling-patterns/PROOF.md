# Prescribing consecutive doublings on two theta fields

## The scores and the arithmetic input

Use $`A=(547,727)`$ and $`B=(571,599)`$. For each pair $`i=(u,v)`$ let $`K=\mathbf Q(\zeta_{16})^+`$, $`F_i=K(\sqrt{uv})`$ and $`E_i=K(\sqrt{-u},\sqrt{-v})`$. A good prime ideal $`w`$ of $`F_i`$ has sign $`\epsilon_w=1`$ for splitting and $`-1`$ for inertia in $`E_i/F_i`$. This is the product over its fibre of the normalized first-step theta trace jumps. The two admissible host data and the common-height convention have been established separately. Omit the places above 2 and the two conductors in the finite prime examples. Any additional fixed finite omissions, including the exceptional primes of the trace data, have left the compact support on a sufficiently late tail. [547/727 trace](https://github.com/leecrossley/theta-race-short-shifts/blob/483c31f2df2132ed2165620ea55e20e969228bfb/RaceTrace.lean#L83) · [571/599 trace](https://github.com/leecrossley/theta-host-split-inert-race/blob/8fca7cad37185e3d7951f0dffff5e2442d1dfa24/HostAdmissible.lean) · [Common-height finite-zero theorem](https://github.com/leecrossley/theta-host-split-inert-race/blob/joint-levels-and-moves-v1.0.0/results/joint-levels-and-moves/PROOF.md)

Set $`a=\log 2`$, $`r=\sqrt2`$, $`c=2(r-1)^2`$, $`B_0=rc`$ and $`W(t)=\max(0,\min(2-2/t,4/t-1))`$. Use

```math
S_i(y)=-\frac{1}{B_0\sqrt y}\sum_w\epsilon_w\log(Nw)W(Nw/y),
\qquad D_i(y)=\frac{S_i(2y)-S_i(y)}a .
```

The odd primitive constituents are the quadratic characters of conductors $`u,v`$, each twisted by the four powers of the primitive even quartic character modulo 16. For each odd prime conductor the primitive conductors are $`p,8p,16p,16p`$. The bound $`|V_i|\le M_i`$ below, the dyadic certificate and the infinite-height conclusions invoke Alpöge's announced theorem that each primitive Dirichlet L-function has its zeros in the critical strip on the critical line. This is GRH for these functions, and the announcement has not yet been refereed. Its statement, the printed axiom audit and the independent kernel records are pinned by digest. The finite prime-ideal catalogue uses no zero-location theorem. The density argument permits every rational relation among the ordinates. [Announcement](https://alpo.ge/remark.pdf) · [Upstream theorem](https://github.com/leanprover-community/riemann-hypothesis/blob/371a1a26d5e17f604428e3a8ce6ce341cf517434/src/Main.lean) · [Independent verification](https://math.mit.edu/~drew/grh/) · [Pinned record check](check_source.py)

## Profiles for the sixteen characters

The compact test and its Fourier transform are

```math
g(h)=\frac4{B_0}\sinh\!\left(\frac{(a-|h|)_+}{2}\right),\qquad
G(t)=\frac{3-2r\cos(at)}{c(t^2+1/4)}.
```

Let $`V_i(x)=\sum_{\gamma} G(\gamma)e^{i\gamma x}`$ over all ordinates of the eight primitive factors for $`i`$, with multiplicity. The family is closed under conjugation, so $`V_i`$ is real and even. The coefficients are positive. The standard zero count and the quadratic decay give absolute uniform convergence. In particular $`|V_i(x)|\le M_i=V_i(0)`$. For $`z\ge1`$ set $`x=\log z`$, $`s=e^{-|h|}`$ in $`T(h)`$ below, and define

```math
T(h)=\frac{(1+s)\log(1+s)+(1-s)\log(1-s)}{2\sqrt s},\quad T(0)=a,\qquad
C_i=\sum_\chi\log(q_\chi/\pi)-8(\gamma_{\rm E}+a).
```

For a rational prime power $`n=p^j`$, let $`t_i(n)=\sum_\chi\chi(p)^j`$ using the primitive local values, zero when the character is ramified. The compact-test explicit formula gives a finite arithmetic evaluation,

```math
V_i(x)=C_i g(x)-\frac8c\{r[T(x+a)+T(x-a)]-3T(x)\}
 -\sum_{n=p^j}t_i(n)\frac{\log p}{\sqrt n}\,g(x-\log n).
```

Only $`z/2\lt n\lt 2z`$ contributes. The reflected prime sum vanishes on $`x\ge0`$, because its argument is at or beyond the kernel boundary. Evenness handles negative $`x`$. To obtain the gamma terms, expand the odd gamma-factor kernel as a sum of exponentials of rates $`2k+3/2`$. Between its three knots, $`g''=g/4`$. Convolution against each exponential leaves three translates, with weights proportional to $`r,-3,r`$. Their scalar sum is $`T(h)=\sum_{k\ge0}s^{2k+3/2}/[(2k+1)(2k+2)]`$. The remaining multiple of $`g`$ gives $`C_i`$. Only the conductors and prime traces change between the two hosts. [Compact-test explicit formula and normalization](https://github.com/leecrossley/theta-race-curvature-euler-factors/blob/af7093b56a7f6b4155dba27e735b5b1d471ff707/PROOF.md#remove-the-archimedean-term) · [Odd gamma convolution at fixed width](https://github.com/leecrossley/theta-race-width-cusp/blob/52eee965886958193d2ff38595f92389b90f4b93/PROOF.md)

For clarity, the passage back to ordinary primes only needs the classical prime ideal theorem for the unsigned square term. At $`x=\log(2y)`$ one has $`g(x-\log n)/\sqrt n=W(n/y)/(B_0\sqrt y)`$. Expanding the primitive Euler factors over the base field gives the signed prime powers with coefficient $`\epsilon_w^j\log Nw`$ away from the finitely many exceptional primes. At each exceptional place the powers in the compact window have total normalized size $`O(1/\sqrt y)`$ for large $`y`$. The archimedean terms tend to zero as $`x\to\infty`$. To pass from the negative prime-power sum to the negative prime sum, restore the positive unsigned squares. The prime ideal theorem and partial summation give their limit as $`B_0^{-1}\int_1^2 W(t^2)\,dt=1`$. The integral follows by integrating the two branches, or from the weight's Mellin value at one half. The total absolute contribution from powers of order at least three is $`O(y^{-1/6}\log y)`$ by summing the prime-ideal Chebyshev bound up to $`(4y)^{1/j}`$ over $`3\le j\le\log(4y)/\log 2`$. Hence, with $`P_i(x)=(V_i(x+a)-V_i(x))/a`$,

```math
S_i(y)=1+V_i(\log(2y))+o(1),\qquad
D_i(y)=P_i(\log(2y))+o(1).
```

Both errors tend to zero on the entire positive tail. The prime-ideal convention for the finite examples permits all unramified rational primes. [Fixed-width prime-only limit](https://github.com/leecrossley/theta-race-width-cusp/blob/52eee965886958193d2ff38595f92389b90f4b93/PROOF.md#the-prime-only-limit-at-a-fixed-width)

## Sixteen pairs of consecutive moves

The [witness file](evidence/witnesses.tsv) has five integer columns. They are a code $`w=0,\ldots,15`$, the numerator and denominator of a positive ratio $`z_w`$, then the numerator and denominator of a finite start $`y_w`$. Write the code in four binary digits, with leading zeros. In order the digits prescribe field A then B on the first doubling, followed by A then B on the next. A one prescribes gain and a zero prescribes loss. Write the corresponding signs as $`s_{w,i,j}`$ for $`j=0,1`$.

Two independent interval evaluations prove that every specified sign satisfies

```math
s_{w,i,j}D_i(2^j y_w)\gt 1.6,\qquad
s_{w,i,j}P_i(\log z_w+j a)\gt 9.1 .
```

All endpoints are exact rationals before transcendental evaluation. The largest finite norm is 747,539 and the explicit formula uses prime powers only through 54,638. The direct calculation builds quadratic residue tables. For odd $`p`$, the sum over the quartic twists vanishes unless $`p^j\equiv\pm1\bmod16`$, when it is four times the sum of the two quadratic values to exponent $`j`$. At 2 only the two odd-conductor characters survive. For finite prime ideals the base residue degree $`f`$ is the least exponent in $`1,2,4`$ that satisfies both $`p^f\equiv\pm1\bmod16`$ and $`(\chi_u(p)\chi_v(p))^f=1`$. There are $`8/f`$ primes of this norm and their common sign is $`\chi_u(p)^f`$. Prefix sums of the coefficients and of the coefficients divided by the norm evaluate the two branches exactly up to the directed logarithm enclosures. [Primary verifier](patterns.py)

The independent verifier factors the polynomials $`X^4-4X^2+2`$, $`Y^2-uv`$ and $`Z^2+u`$ over each good residue field. Least common multiples of the local degrees give the residue degree in the base and in the extension. It sums the weights directly at the rational starts. For the spectral calculation it evaluates the sixteen primitive Conrey characters and the gamma series. The remainder after 128 terms of that series is at most $`s^{257+1/2}/[(257)(258)(1-s^2)]`$ for $`0\lt s\lt 1`$. The value at one is handled separately. Both calculations also enclose $`0\lt M_i\lt 180`$. [Independent verifier](independent.py) · [Expected transcript](evidence/check.txt)

## Every six-step word on the dyadic grid

For a word of length $`m`$ write each joint direction as a base-four digit, with its two binary bits giving the signs for fields A and B in that order. The first step is the most significant digit. The [dyadic table](evidence/dyadic-starts.tsv) contains a length, code and signed integer $`k`$ for **every** code of lengths four, five and six. The log centre for the first move is $`k a`$. The associated prime start after an integer translation $`h a`$ is $`y=2^{h+k-1}`$. A negative $`k`$ is a legitimate profile sample by evenness of $`V_i`$. All translated prime starts are positive integer powers of two for sufficiently large $`h`$.

For each code the certificate proves a strict signed lower bound for every $`P_i((k+j)a)`$, $`0\le j\lt m`$. The respective lower bounds are 14.1, 9.1 and 4.1 for lengths four, five and six. There are 5,376 tabulated witnesses and 61,440 signed component inequalities. The zeros and initial exponents can be chosen by floating-point search. The two replays below operate on the exact published integers. [Starting exponents](evidence/dyadic-starts.tsv) · [Rational thresholds](evidence/dyadic-margins.tsv)

The [bracket file](evidence/dyadic-zeros.tsv) gives the host index, primitive conductor, Conrey character label and the two endpoints in units of $`10^{-14}`$. The checker verifies all sixteen primitive odd characters from the conductor frame, all 28,000 Hardy signs, and disjointness within each factor. Select one actual positive ordinate from each of 6,000 brackets in field A and 8,000 in field B and pair it with its negative in the aggregate family. Multiplicity and any additional positive zeros remain in the tail. The companion of a nonreal character belongs to its conjugate. The total selected term at step $`j`$ from an ordinate $`t`$ is

```math
w(t)\sin((k+j+1/2)at),\qquad
w(t)=-4G(t)\sin(at/2)/a .
```

The all-ordinate mass is $`M_i=V_i(0)`$. Its sine-square moment is $`(V_i(0)-V_i(a))/2`$. From each subtract the selected pairs' contributions, respectively $`2G(t)`$ and $`2G(t)\sin^2(at/2)`$. Write $`R_0,R_2`$ for the two residuals. Cauchy–Schwarz gives omitted move amplitude at most $`2\sqrt{R_0R_2}/a`$. This subtraction needs one sign-changing zero per bracket. The proof permits additional zeros anywhere, with arbitrary multiplicity. [Endpoint verifier](check_grid_zeros.py) · [Full moment and primary grid calculation](grid_certificate.py)

For the finite head, set $`\theta=t a/(4\pi)`$ in turns. Approximate $`\theta`$ by an integer tick divided by $`2^{64}`$, and $`w(t)`$ by a signed integer divided by $`2^{24}`$. The script computes enclosing balls over the full root bracket. If the quantized values are $`\hat\theta,\hat w`$ and the odd multiplier is $`d=2k+2j+1`$, the error of one term is bounded by

```math
|w-\hat w|+|\hat w|\left(2\pi |d|\,|\theta-\hat\theta|+\frac{\pi}{32768}+2^{-24}\right).
```

The last two terms bound nearest-node phase sampling and the rounding of its sine. The code checks every sine-table rounding. Unsigned 64-bit multiplication computes turns modulo one and therefore performs exact phase reduction even at the large exponents. The signed dot product fits in 64 bits by the checked sum of absolute weights. For negative starts the primary verifier uses $`P_i(-x-a)=-P_i(x)`$, reversing the order and both signs. The transcript checks every required component after the phase errors and the omitted spectral mass are subtracted. [Fixed-point evaluator](grid_sum.cpp) · [Reflection identity](DoublingPatterns.lean)

The independent evaluator uses cosine at the centres of 65,536 cells, with signed weights and ordinates quantized by taking floors. It evaluates negative integer exponents directly. Its tail bound uses the cubic $`J(u)=(11+75u-48u^2+16u^3)/54`$. For $`z=\sqrt u`$ the difference $`J(z^2)-z`$ equals $`(8/27)(z-1/2)^2(z-1)^2((z+3/2)^2+1/2)`$. Thus the unlisted amplitude is bounded by the unlisted sum of $`(2/a)G(t)J(\sin^2(at/2))`$. The full mass for this polynomial is determined by $`V_i(0),V_i(a),V_i(2a),V_i(3a)`$ via the elementary power-reduction identities. The independent source evaluates these four moments from the gamma series and the primitive terms at 2. Through norm 16 the odd-prime trace vanishes and the boundary terms have zero weight. [Independent shifted-grid and cubic calculation](grid_independent.py) · [Original cubic moment proof](https://github.com/leecrossley/theta-race-factor-400/blob/f1773b07a55e76947d4e0c6ab19a112fd91fe6cf/PROOF.md)

## A common return with an explicit density

Write $`Q=11632`$. Let $`N_\chi(t)`$ count zeros at both signs of the ordinate, with multiplicity. Corollary 1.2 of Bennett, Martin, O'Bryant and Rechnitzer bounds the absolute error in $`(t/\pi)\log(qt/(2\pi e))`$ by $`0.247\log(qt)+6.894`$. For each of $`T=600,250000`$ this gives $`N_\chi(t)\le(t/\pi)\log(Qt)`$ for every $`t\ge T`$. At $`T`$ the discarded term $`(T/\pi)\log(2\pi e)`$ dominates the error, and the difference increases thereafter by differentiation. [Counting zeros of Dirichlet L-functions](https://doi.org/10.1090/mcom/3599) · [Directed constant checks](return_bounds.py)

Since $`G(t)\lt 17/t^2`$ for $`t\gt 0`$, partial summation, dropping the nonpositive boundary term at $`T`$, bounds the omitted mass in either host by

```math
R(T)\le 8\cdot34\int_T^\infty\frac{\log(Qt)}{\pi t^2}\,dt
 =\frac{272(\log(QT)+1)}{\pi T}.
```

Across the two hosts there are at most $`d=8T\log(QT)/\pi`$ positive phase coordinates with ordinate at most $`T`$. Indeed there are sixteen factors, and the aggregate is conjugation stable. Central terms have constant phase. Repetitions can be counted once. If every retained positive phase changes by at most $`2\pi/H`$ in angle modulo $`2\pi`$, the change in each $`V_i`$ is bounded uniformly by $`180(2\pi/H)+2R(T)`$. The change in the normalized profile difference $`P_i`$ is consequently at most

```math
E(T,H)=\frac2a\{180(2\pi/H)+2R(T)\}.
```

For $`T=250000`$, $`H=2^{18}`$, the certificate gives $`E(T,H)\lt 0.06`$ and $`18d\lt 252000000`$. Consider the translations $`t`$ for which the distance from $`\gamma t/(2\pi)`$ to the nearest integer is below $`1/H`$ for all retained coordinates. The product of circle tents $`(1-H\|\gamma t/(2\pi)\|)_+`$ is a nonnegative minorant of the indicator of this event. Each tent has constant Fourier coefficient $`1/H`$, nonnegative absolutely summable Fourier coefficients and value one at the origin. This follows by convolving two interval indicators of width $`1/H`$ and scaling to peak one. The product's zero multi-index has coefficient at least $`H^{-d}`$. Coincident zero frequencies add nonnegative terms.

On any long real interval the mean of a nonconstant character tends uniformly in the initial point to zero. On the lattice $`t=k h`$ with any fixed real $`h\gt 0`$ the finite geometric-series formula gives the corresponding uniform limit for every frequency whose sampled character is nonconstant. Sampled constants again contribute positively. Absolute summability permits taking limits in both cases. Thus the returning real shifts have lower uniform density at least $`2^{-252000000}`$, and the returning integer exponents have the same lower uniform density on every fixed lattice. Here lower uniform density means the limit inferior, as interval length tends to infinity, of the infimum over starting positions on the positive half-line of the fraction of that interval occupied by the set.

At any sufficiently large returning shift use the sixteen actual starts $`y=z_w e^t/2`$. The error in each prime-only move can be taken below $`0.01`$, simultaneously over both fields and both steps. Their signed values exceed $`9.1-0.06-0.01\gt 9`$. So all sixteen two-step patterns occur at arbitrarily large real starts, with every normalized move exceeding 9. For the dyadic words take the lattice of translations with spacing equal to the logarithm of two. The same shifts preserve all 5,376 finite-head witnesses simultaneously. They exceed the table thresholds minus 0.07 in the actual prime moves. Consequently every four-, five- and six-step code occurs at arbitrarily large integer power-of-two starts with normalized component margins greater than 14, 9 and 4, respectively. Each code has its own fixed integer offset, so its successful starting exponents themselves have lower uniform natural density at least the same bound. A common real translation similarly gives positive lower uniform logarithmic density among real starting heights. [Earlier one-host common-return argument](https://github.com/leecrossley/theta-race-paired-density/blob/5181e3a8545c67f1bc93f344b39683cab4a10383/PROOF.md)

There is an explicit height-dependent upper bound on the common return. Let $`K=2143656793`$ be the largest absolute exponent in the published table. For a sufficiently large integer $`N`$ put $`L=N+K+2`$. Partition the retained phase torus into $`H^n`$ cells as above, with $`n`$ the number of distinct retained positive ordinates. Among the $`H^n+1`$ points belonging to the log shifts $`uLa`$ for $`0\le u\le H^n`$, two share a cell. Their positive index difference $`m`$ satisfies $`1\le m\le 2^{252000000}`$. The shift $`mLa`$ therefore preserves all the certified margins. Every prescribed word then starts at a power of two with exponent $`mL+k-1\gt N`$. Its last doubling finishes at an exponent bounded above by $`2^{252000000}(N+K+2)+K+5`$. This inequality bounds exponents. The lower cutoff for $`N`$ comes from the prime-only convergence.

## What longer patterns can occur

Fix a tail on which the errors in the normalized moves are below $`0.01`$. For $`T=600`$ and $`H=2^{10}`$ the directed check gives $`E(T,H)\lt 17.9`$ and $`10d\lt 241000`$. Divide the circle into $`H`$ equal half-open intervals at each retained frequency. If two starting log centres fall in the same phase cell, their angular differences are below $`2\pi/H`$. Adding $`ja`$ to both centres preserves each angular difference. In every position of the two actual words their normalized moves therefore differ by less than $`17.9+0.02\lt 18`$ in each field. Opposite signs both exceeding 9 in magnitude would differ by more than 18. Thus two distinct words of any common length with all moves of magnitude greater than 9 require different cells. There are at most $`2^{241000}`$ such words, independently of their length.

There is also a common integer repetition lag for every sufficiently late word at the smaller margin $`1/25`$. With the finer partition $`T=250000`$, $`H=2^{18}`$, place the phase vectors of $`ka`$ for integers $`0\le k\le H^n`$ in the $`H^n`$ cells, where $`n`$ is the actual number of retained distinct positive ordinates. Two share a cell. Their index difference $`m`$ satisfies $`1\le m\le2^{252000000}`$. It moves every phase by less than $`2\pi/H`$, at every starting height. Thus whenever positions $`j,j+m`$ in a late word both have normalized magnitude greater than $`1/25`$ in a field, their signs agree, because the difference is less than $`0.06+0.02=2/25`$. The same $`m`$ works for both fields, every length and every sufficiently late start.

Finally $`2M_A\lt 57\cdot9a`$ and $`2M_B\lt 56\cdot9a`$ by the directed mass calculation. Adding that many consecutive normalized moves of the same sign and magnitude greater than 9 would give a score difference larger than the corresponding right side. The actual endpoint scores each lie between $`1-M_i-o(1)`$ and $`1+M_i+o(1)`$, a contradiction sufficiently far out. These run lengths are explicit upper obstructions. The exact maximum run length at margin 9 is open. [Budget checks](return_bounds.py) · [Sign, telescoping and counting deductions in Lean](DoublingPatterns.lean)
