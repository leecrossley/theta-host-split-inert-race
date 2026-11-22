# Which primes survive dyadic sampling?

## The two readings

Set $`a=\log 2`$, $`B_0=2\sqrt2(\sqrt2-1)^2`$ and $`W(t)=\max(0,\min(2-2/t,4/t-1))`$ for $`t\gt0`$. For $`i=A,B`$ take $`(u,v)=(547,727),(571,599)`$ respectively. Put $`K_0=\mathbf Q(\zeta_{16})^+`$, $`F_i=K_0(\sqrt{uv})`$ and $`E_i=K_0(\sqrt{-u},\sqrt{-v})`$. Omit $`2,u,v`$ on each field. Also omit 343057911569 and 343070636977 on A, and 413718278401 and 603339156001 on B. These are fixed exclusions from the trace data. All logarithms without subscripts are natural.

For a remaining rational prime $`p`$, put $`\chi_u(p)=(p/u)`$ and define $`f_i(p)`$ as the least $`d\in\{1,2,4\}`$ with $`p^d\equiv\pm1\bmod16`$ and $`(\chi_u(p)\chi_v(p))^d=1`$. Write $`m_i(p)=p^{f_i(p)}`$, $`s_i(p)=-\chi_u(p)^{f_i(p)}`$ and $`L_p=\log(p)/a`$. The residue degree in the real cyclotomic subfield is the order modulo $`\{\pm1\}`$. The quadratic factor adds the second condition. Over $`F_i`$ the Frobenius of $`E_i/F_i`$ has sign $`\chi_u(p)^{f_i(p)}`$ by quadratic reciprocity. There are $`8/f_i(p)`$ primes of norm $`m_i(p)`$ above $`p`$. The family’s local trace identity identifies their split/inert signs with the products of the normalized first-step theta trace jumps. Thus the race with $`b=r=2`$ and its normalized change have the exact finite expansions

```math
S_i(y)=\frac{8}{B_0\sqrt y}\sum_p s_i(p)\log(p)W(m_i(p)/y),
\qquad D_i(y)=\frac{S_i(2y)-S_i(y)}{a}.
```

[547/727 fibre-product lemma](https://github.com/leecrossley/theta-host-split-inert-race/blob/66ba5ad188650e38cbf2e70e4ef6c4960d96d8c7/results/theta-race-short-shifts/RaceTrace.lean#L83) · [571/599 host and sign](https://github.com/leecrossley/theta-host-split-inert-race/blob/8fca7cad37185e3d7951f0dffff5e2442d1dfa24/PROOF.md) · [Prime-degree conventions](https://github.com/leecrossley/theta-host-split-inert-race/blob/consecutive-doubling-patterns-v1.4.1/results/consecutive-doubling-patterns/PROOF.md)

Use the field $`K=\overline{\mathbf Q}\cap\mathbf R`$ and only **finite** $`K`$-linear combinations of $`D_i(2^j)`$, $`j\ge0`$. Readings with negative exponents would be zero. At a dyadic argument all coefficients of $`L_p`$ are algebraic. Unique factorization makes the ordinary logarithms of distinct primes linearly independent over the rationals. Baker’s theorem makes them linearly independent over the algebraic numbers. Dividing by $`a`$ preserves independence, including $`L_2=1`$ when an independent constant term is allowed. The **support** of a value is consequently its uniquely determined set of odd prime coordinates with nonzero coefficient. Isolation means that the resulting value equals $`L_p`$, after algebraic rescaling if necessary. [Baker, *Linear forms in the logarithms of algebraic numbers*](https://doi.org/10.1112/S0025579300003971)

## How many grids isolate a prescribed prime?

A grid has a real algebraic phase $`1\le\alpha\lt2`$ and supplies every $`D_i(\alpha2^j)`$ for integers $`j\ge0`$. A finite collection supplies the union of those readings. The phases and the coefficients of a finite linear combination may depend on the target $`p`$. The [table](README.md) gives the exact minimum on each field. Here a good prime is any prime retained in that field's finite expansion above.

### Three rational phases always suffice

Put $`T_i(z)=B_0\sqrt z\,S_i(z)/(8a)`$. This is the sum of $`s_i(q)L_q W(m_i(q)/z)`$. The initial level at every phase vanishes because $`4\alpha\lt8`$ and the smallest good norms are 31 and 49. Telescoping the changes therefore obtains every level on the grid. The final scaling is algebraic at an algebraic argument.

Fix a good prime $`p`$ and write $`m=m_i(p)`$. In the level argument $`z`$ each summand is a piecewise linear function with possible knots $`m_i(q)/4,m_i(q)/2,m_i(q)`$. Multiplication by four takes every knot to an integer. The target norm $`m`$ is odd, so the only knot strictly between $`(m-1)/4`$ and $`(m+1)/4`$ is the first knot of the target. Unique factorization identifies its rational prime. Its slope increases there by $`4/m`$. All other terms have zero centered second difference, including a term with a knot at an endpoint. Hence

```math
T_i((m+1)/4)-2T_i(m/4)+T_i((m-1)/4)=\frac{s_i(p)}{m}L_p .
```

Each of the three arguments is a positive rational larger than two and belongs to a rational dyadic grid. Multiplication by $`m/s_i(p)`$ recovers the target. The three-score inversion and its interpretation for local Euler factors were proved in the earlier paper. We use its odd-norm case at the first norm of $`p`$. [Finite Euler-factor readout](https://doi.org/10.5281/zenodo.24413208) · [Kernel-checked second differences](DyadicPrimeSupport.lean)

### Why two grids have a finite cutoff

On one field any finite combination of levels has coefficient $`s_i(p)F(m_i(p))`$ at $`L_p`$, where $`F`$ is continuous, compactly supported and affine in $`1/m`$ between consecutive knots of the two grids. The knot sequences are $`\alpha2^j`$ and $`\beta2^j`$. A repeated grid can be augmented by a distinct grid with zero coefficients. Extend the sequences to negative exponents for convenience, with zero values before the first contributing level. If the combination isolates a prime, Baker's theorem requires every other prime coefficient to vanish.

Let $`\kappa`$ be the largest knot at which $`F`$ is nonzero, and let $`u,v`$ be its neighbours. The union of two distinct dyadic grids alternates, so $`v=2u`$. On $`(u,\kappa)`$ the function has at most one zero. On $`(\kappa,v)`$ it interpolates between a nonzero value and zero and has no zero. At $`\kappa`$ its value is nonzero. Thus three different prime norms strictly inside $`(u,2u)`$ force at least two prime logarithms in the support. Beyond $`v`$ the spline is identically zero. [Affine zero lemmas](DyadicPrimeSupport.lean)

On each field, every real $`t\gt289/2`$ has at least three good norms in $`(t,2t)`$. The strict lower cutoff is sharp. At $`t=289/2`$ the only such norms are 223 and 271 on A, and 239 and 241 on B. The norm 289 of 17 is at the excluded right endpoint on both. Three checks cover all larger real starts.

- Up to 6554 an exact sieve lists every good norm through 13108. At $`x=2t`$ a norm $`m`$ is present precisely when $`m\lt x\lt2m`$. At each event the sweep removes closing norms, counts at the point, adds opening norms, and counts to the next event. A separate implementation checks the distance to the third successor after each good norm. [Event sweep](minimum_grids.py) · [Consecutive-norm check](minimum_grids_recheck.py)
- For $`6554\le t\le2^{128}`$ there are 359 consecutive closed intervals $`[l,r]`$. Each carries three distinct rational primes of specified common residue degrees whose norms satisfy $`r\lt m\lt2l`$. Trial division proves primality. Euler's criterion and the orders modulo 16 check the norms on both fields. The largest rational prime in the table is 4,294,967,387. [1,077 witnesses](evidence/common-octaves.tsv) · [Deterministic reconstruction](grid_octaves.py)
- For $`t\ge2^{128}`$ put $`X=t^{1/4}`$. The [modulo-eight estimate](#minimum-support-on-one-field) supplies more than eight primes congruent to 3 modulo 8 in $`(X,19X/16]`$. Each retained prime has degree four on both fields, and $`(19/16)^4\lt2`$. Omitting the five exceptions of either field still leaves at least three for that field. [Explicit prime progression bounds](https://arxiv.org/abs/1802.00085)

Consequently a two-grid singleton must have $`v\le289`$ and its norm must be strictly less than 289. This argument allows arbitrarily large indices and any finite number of readings from the two grids. It applies also to the case of a single grid augmented by an unused phase.

### The small norms and every possible grid placement

The good norms below the cutoff are as follows. A superscript gives the residue degree when it exceeds one.

| Field | Norms in ascending order |
| --- | --- |
| A | $`31,47,7^2,79,3^4,97,113,127,223,271`$ |
| B | $`7^2,79,3^4,113,127,239,241`$ |

Take all union knots up to 289. Relax the problem by allowing arbitrary nodal values except that the final knot is zero and the right tail is zero. Every true singleton obeys those constraints because its last nonzero knot has a successor at most 289. For each proposed target impose a zero at every other listed norm. A zero at a knot fixes its value. Two different zeros in an open cell force both endpoints to zero. A single interior zero propagates a forced zero at either endpoint to the other endpoint. Iterate along the finite path of knots. A target at a forced-zero knot or between two forced-zero endpoints is impossible. A target beyond the final knot is impossible as well. This is a necessary test and explicit formulas below establish sufficiency for every surviving prime.

The test uses only the relative order of knots, norms and the ceiling. For one phase the only changes of order occur at the normalized phases $`m/2^{\lfloor\log_2 m\rfloor}`$ of the listed norms and 289. For two phases also record their order, including equality. Within any one open phase interval, two rational representatives cover both a coincidence and a strict ordering, because the grids can be exchanged. Take each boundary separately, including 1 and excluding 2. This gives 36 representatives and 666 unordered pairs with repetition on A, and 27 representatives and 378 pairs on B. The propagated constraints leave exactly the seven primes per field in the [overview](README.md). In particular 3, 79 and 97 are eliminated on A despite their small norms. A second computation changes the interior representatives and computes exact rational row spaces in reciprocal norm. [Order computation](minimum_grids.py) · [Independent collocation matrices](minimum_grids_recheck.py)

For the remaining primes there are fourteen explicit lists of rational coefficients. A row with phase $`\alpha`$, exponent $`j`$ and coefficient $`c`$ contributes $`cT_i(\alpha2^j)`$. Direct summation gives coefficient 1 at the prescribed $`L_p`$ and 0 at every other prime. The largest support endpoint of every formula is below 512. The independent evaluation uses the positive part $`x_+=\max(x,0)`$ and the identity

```math
m W(m/z)=(4z-m)_+-6(z-m/2)_++2(z-m)_+ .
```

For example on B, $`L_{79}=-(2449/15)T_B(16)+79T_B(20)`$. The phases are 1 and $`5/4`$. The full one-grid classification below gives the unique one-grid target on either field. The other six successes on each field therefore need exactly two. The three-score identity achieves the lower bound for every remaining good prime. [Exact coefficients](evidence/minimum-grid-formulas.tsv) · [Single-phase classification](#the-thirteen-phase-cells)

### Two grids using both fields

The pooled table has three one-grid primes, nine further two-grid primes and every other odd prime has minimum three. In particular 97 needs three grids on either field individually and only two jointly. The two phases are $`3/2,7/4`$. Writing $`T_i`$ for the same normalized levels as above gives the explicit identity

```math
\frac{90}{97}L_{97}
=-2T_A(12)+5\bigl(T_A(24)-T_B(24)\bigr)-6\bigl(T_A(28)-T_B(28)\bigr).
```

The largest contributing norm is strictly less than 112. Both signed expansions have been checked at every contributing prime. This supplies the one target outside the union of the two single-field lists. [Five coefficients](evidence/pooled-grid-formula.tsv) · [Independent evaluation](pooled_grid_recheck.py)

For two grids take the vector of splines $`F=(F_A,F_B)`$ and its last nonzero knot $`\kappa`$ with neighbours $`u,v=2u`$. A prime $`q\in(u,v)`$ of degree one on at least one host has ordinary-prime trace vector $`e(q)`$. Use only primes good on both hosts here. When $`u\gt2`$ every higher-degree occurrence of the same prime has norm at least $`q^2\gt v`$, beyond both splines. There are four projective directions of a nonzero trace vector. For each direction $`d`$ with $`F(\kappa)\cdot d\ne0`$ the scalar spline $`F\cdot d`$ can lose at most one prime in this entire interval, by the two-cell argument. At most one direction is perpendicular to the knot value. If $`n_d`$ counts the ordinary primes of each direction in $`(u,2u)`$, a lower bound for the number of surviving prime coordinates is therefore

```math
B(u)=\sum_d(n_d-1)_+-\max_d(n_d-1)_+ .
```

Every $`u\gt431/2`$ has at least two forced survivors. Up to $`u=13000`$ an endpoint sweep computes $`B(u)`$ exactly. Its last value below two occurs precisely at $`431/2`$, where the four counts are 0, 1, 2 and 4. An independent count bisects four sorted prime lists at every possible endpoint and between endpoints. For $`26211/4\le u\lt2^{1585}`$ the unconditional sixteen-pattern certificate supplies at least 4, 4, 2 and 2 primes of the four directions in $`(u,4u/3]`$, so $`B(u)\ge5`$. For $`\log u\ge1095`$ the estimate in the next section supplies more than $`u/(500\log(2u))`$ primes with nonzero projection at the knot. At most four are lost. This last range overlaps the certified chain since $`1585\log2\gt1095`$. The finite chain is the published list of 66,016 witnesses. Its replay checks every connecting ratio. The complete argument is unconditional. It also proves the sharp ordinary-prime statement that *every* direction sees at least two primes in $`(u,2u)`$ for every $`u\gt271/2`$, using the earlier [initial sweep](#the-last-scale-that-needs-symbolic-arithmetic) through $`431/2`$. The endpoint $`271/2`$ misses. [Ordinary-prime sweep](pooled_octaves.py) · [Independent bisection](pooled_grid_recheck.py) · [Finite chain](https://github.com/leecrossley/theta-host-split-inert-race/blob/joint-frobenius-primes-v1.0.2/results/joint-frobenius-primes/PROOF.md)

### An unconditional overlap using paired zeros

For $`\log X\ge1095`$ every real unit direction sees more than $`X/(500\log(2X))`$ good ordinary primes in $`(X,2X)`$. Use the same base classes, 24 primitive pair twists, exceptional-zero checks and expansion of the squared trace as in the [fourfold smoothing proof](#ordinary-primes-in-every-direction). Now put $`h=(\log2)/2`$ and use the convolution of **two** uniform densities on $`[0,h]`$. Its density $`g`$ has supremum at most $`1/h\lt3`$ and total variation at most $`2/h\lt6`$. Partial summation on the two base classes gives their prime mass greater than $`1/4-14/(160\log X)\gt249/1000`$. The transform is $`K(z)=((e^{hz}-1)/(hz))^2`$. The same Mellin shift is absolutely summable with decay of order two. For a zero $`\rho=\beta+i\gamma`$ its modulus is bounded by $`\min(1,(b/|\gamma|)^2)`$ with $`b=29/5`$. [Progression bounds](https://arxiv.org/abs/1802.00085) · [Zero count, Corollary 1.2](https://arxiv.org/abs/2005.02989)

Take $`Q=7000000`$, $`R=9.65`$ and $`L=1095`$. McCurley's denominator $`9.645908801\log\max(q,q|\gamma|,10)`$ gives the bound $`X^{\beta-1}\le\exp(-L/(R\log(Q\max(1,|\gamma|))))`$ after the [twelve real-character exclusions](#a-zero-free-strip-and-its-real-exception). All conductors exceed $`300000`$, so the applicable real strip boundary lies above 0.99. The functional equation and conjugation pair every zero with $`1-\bar\rho`$ of the **same** primitive character at the same ordinate, preserving multiplicity. The zeros with $`\beta\gt1/2`$ thus contribute at most half the sum of any common nonnegative ordinate majorant over all zeros. Zeros with $`\beta\le1/2`$ receive the factor $`e^{-L/2}`$. [McCurley's region](https://www.sciencedirect.com/science/article/pii/0022314X84900891) · [Functional equation](https://dlmf.nist.gov/25.15.i)

Write $`\nu=8/5`$ and $`t=|\gamma|`$. For $`t\le b`$ the exponential bound is less than $`1/625`$, since $`\log(Qb)\lt17.6`$. For $`t\ge b`$, differentiate $`-L/(R\log(Qt))`$ with respect to $`\log t`$. Its derivative is decreasing and is at most $`L/(17^2 R)\lt2-\nu`$. Multiplying by the transform bound gives for every zero the common envelope

```math
X^{\beta-1}|K(\rho-1)|\le \frac1{625}\min\bigl(1,(b/t)^\nu\bigr).
```

Let $`N(t)`$ count both signs of ordinate. Integrating this last minimum without its factor $`1/625`$ against $`dN(t)`$ gives $`\int_b^\infty \nu b^\nu N(t)t^{-\nu-1}\,dt`$. The BMOR bound with $`q\le Q`$, $`\pi\gt3.14`$ and $`\log(2\pi e)\gt2.8`$ bounds this integral by

```math
\frac{\nu b}{3.14(\nu-1)}
  \left(\log(Qb)-2.8+\frac1{\nu-1}\right)
 +0.247\left(\log(Qb)+\frac1\nu\right)+6.894
 <93 .
```

Thus the right-half zeros cost less than $`93/1250`$, and the others less than $`93e^{-1095/2}`$. The latter uses the same integrable majorant without the factor $`1/625`$. The trivial-zero series and the removal of proper prime powers each cost less than $`10^{-6}`$, by the bounds in the [explicit formula](#the-smoothed-character-sum) with this $`h`$. The prime-power bound is proportional to $`e^{-(\log X)/2}(\log X+\log 2)^2`$, decreasing throughout this range. Each twist has absolute prime mass below $`3/40`$. The squared-trace expansion has total absolute twist coefficient at most $`3/2`$. Its base coefficient is $`1/2`$. The resulting weighted mass exceeds $`249/2000-(3/2)(3/40)=3/250`$, and a single prime contributes at most $`6\log(2X)/X`$. This gives the claimed lower count. A separate majorant with $`b=231/40`$, $`R=9.646`$ and $`\nu=13/8`$ gives twist error below $`1/14`$. Both computations check the constants by rational inequalities. [First envelope](pooled_octaves.py) · [Second envelope](pooled_grid_recheck.py)

### Four candidates and every pair of phase cells

A pooled singleton has terminal endpoint at most 431 by the preceding octave obstruction. In fact **each level used can be taken to have right endpoint at most 431**. Combine repeated levels into one coefficient vector. For two different phases the right endpoints $`4\alpha2^j`$ are distinct across phases and indices. The level with largest right endpoint is the only contributing level immediately to the left of that endpoint, hence cannot cancel there. For coincident phases first combine their coefficients. An unused phase may be added when applying the two-cell bound.

The union of both fields contains sixteen candidate primes of norm strictly below 431. The eleven single-field successes and the formula for 97 leave only 17, 337, 353 and 367 to exclude. Whenever a prime occurs on both fields in this prefix its norm is the same there. The full phase partition has eighteen boundary points in $`[1,2)`$ and eighteen open intervals. They come from the normalized candidate norms, the cutoff 431, and the endpoints 1 and 2. Thus there are 666 unordered pairs of cells, with coincidences included. For each open cell each allowed row of level coefficients is affine in its phase. Eligibility of a level with $`4z\le431`$ is constant in the open cell. For a boundary point use precisely those levels eligible there. Interpolation between the two endpoint rows of an open cell uses the open cell's own fixed set of levels.

For a given pair of cells let the columns correspond to the sixteen primes. First stack both endpoint rows from every level available to either grid, allowing their coefficients to vary independently. This enlarged row space contains every actual row space on the pair of cells. If it excludes the desired coordinate vector, that prime is impossible throughout. This proves 2,407 of the 2,664 exclusions, including every same-cell pair. The other 257 have explicit rational column-dependence certificates of rank at most six. [Certificate table](evidence/pooled-grid-minors.tsv)

Each certificate chooses columns consisting of the target and $`r`$ others. It proves that the full chosen group has rank at most $`r`$ and the other columns have rank $`r`$. A relation expressing the target column in terms of the others prevents any row combination from isolating it. For a variable phase write each row up to positive scaling as $`U+sV`$ with $`s\gt0`$ and the two endpoint rows $`U,V`$. Distinct open cells use independent variables, and a point cell uses a constant row. Each minor of size $`r+1`$ has degree in each variable at most the smaller of $`r+1`$ and the number of its available variable rows. Ranks at a product of more points than those two degree bounds prove the upper bound identically by polynomial interpolation. An insufficient number of rows gives the upper bound immediately. For the lower bound the certificate gives $`r`$ row indices. In the corresponding determinant all nonzero monomial coefficients have the same sign, so it is nonzero throughout the open cell pair. The primary replay obtains these coefficients by multilinearity in the rows. The second uses exact collocation at equally spaced phase coordinates and interpolation in the tensor-product Bernstein basis. Both replays check boundary cells explicitly. [Homogeneous determinants](pooled_grids.py) · [Independent interpolation](pooled_grid_recheck.py)

The two-grid pooled list is therefore complete. The [one-grid classification](#the-thirteen-phase-cells) identifies the three primes needing only one phase. The hosts' excluded odd primes are disjoint, so every odd prime is good on at least one field and the three-level identity applies to it.

### Fixing the grids in advance

Finally fix any finite nonempty collection of algebraic phases, removing duplicates. Their union of knots is a fixed multiplicative pattern repeated every octave. The minimum ratio $`g`$ of successive knots is strictly larger than one, with $`g=2`$ for one phase. For every nonzero combination take its last nonzero vector knot $`\kappa`$. On $`(\kappa,g\kappa)`$ the vector is a strictly positive scalar multiple of its value at $`\kappa`$. The two axial degree-one trace directions are each realized by fixed reduced residue classes of positive density, by the Chinese remainder theorem. For every sufficiently large $`\kappa`$ the prime number theorem in these fixed progressions supplies at least two primes of each direction in this interval. At least one of the directions has nonzero projection at the knot. Higher-degree contributions of each such prime lie above the next knot, which is at most $`2\kappa`$. A singleton consequently has bounded last knot, hence bounded prime. This argument covers any fixed finite bank even if the coefficients and reading indices can depend on the target. [Residue classes and densities](#how-many-ordinary-primes-must-remain)

## Moving the sampling grid

Fix $`1\le\alpha\lt2`$ and now read $`D_i(\alpha2^j)`$ for all nonnegative integers $`j`$. First suppose $`\alpha`$ is real algebraic. Algebraic row scalings and telescoping replace all changes through exponent $`n`$ by the levels $`T_i(\alpha2^j)`$ for $`1\le j\le n+1`$, where $`T_i(z)=B_0\sqrt z\,S_i(z)/(8a)`$. The initial level is zero because the smallest good norms on A and B are 31 and 49. Baker’s theorem identifies numerical isolation with isolation of a prime column. The table in the [overview](README.md) gives every algebraic phase and all its boundary conventions.

### The last scale that needs symbolic arithmetic

For an arbitrary real last start $`y=\alpha2^n\gt271/8`$, every nonzero terminal coefficient pair forces at least two distinct rational-prime coordinates, regardless of the earlier coefficients. This statement concerns the formal prime coefficients and holds for **all** real phases. The coefficient of any degree-one prime $`4y\lt p\lt8y`$ comes entirely from the terminal pair, with a common strictly positive weight. For $`271/8\lt y\le2000`$ list the projective directions of the nonzero degree-one trace vectors. Under the parameter $`t=8y`$ an eligible prime is active exactly for $`p\lt t\lt2p`$. At each event remove closing primes, test the point, add opening primes, and test the following interval. The number visible in the least favourable direction is the total minus the largest count in one projective direction. The final point where this count is below two is exactly $`t=271`$. An exact sieve through 16,000 proves this over the entire real range. [Endpoint sweep](phase_bridge.py)

For $`2000\lt y\lt2^{1583}`$ put $`L=4y`$. The prior unconditional prime chains give all sixteen four-sign patterns in $`(L,4L/3]`$ whenever $`26211/4\le L\lt2^{1585}`$. This forces three projectively distinct trace directions in $`(4y,8y)`$. For $`y\ge2^{2900}`$ the [ordinary-prime estimate](#ordinary-primes-in-every-direction) applies directly to $`X=4y`$. Both parts cover arbitrary real $`y`$. [Finite Frobenius-chain proof](https://github.com/leecrossley/theta-host-split-inert-race/blob/joint-frobenius-primes-v1.0.2/results/joint-frobenius-primes/PROOF.md)

### A moving knot in the intermediate range

It remains to handle every $`2^{1583}\le y\le2^{2900}`$. If $`\lambda_A+\lambda_B\ne0`$, use common degree-four primes $`p\equiv3\bmod8`$ in the outer norm octave $`(4y,8y)`$. The [modulo-eight argument](#minimum-support-on-one-field) uses the real variable $`X=(4y)^{1/4}`$ and supplies at least three. Their coefficients have the same sign on both fields and only the final pair can contribute.

If $`\lambda_A+\lambda_B=0`$ and the pair is nonzero, set $`x_0=2\sqrt y`$ and take the first geometric bin with lower endpoint $`l\gt x_0`$, starting from $`2^{24}`$ and replacing $`l`$ by $`\lfloor13l/10\rfloor`$ at each step. Then $`x_0\lt l\le13x_0/10`$. The twelve labelled primes used in this bin satisfy $`l\le p\le65l/64`$. Consequently $`4y\lt p^2\lt8y`$, since $`(169/128)^2\lt2`$. The bin lies inside the tabulated bit lengths 785 through 1460. Its twelve primes have nonzero degree-one vectors $`(a,b)`$ with entries in $`\{-1,0,1\}`$. A zero coordinate here means that this field contributes at the square. The four omitted labels have degree two on both fields. [21,432 prime certificates](evidence/phase-proth-witnesses.tsv)

The ratio of the largest to the smallest selected prime is at most $`65/64\lt2`$. Thus at most one grid knot can separate the ordinary primes, even as $`\alpha`$ varies. Writing $`x=1/p`$, the ground-level contribution after dividing by $`L_p`$ is $`a(A_0+A_1x)+b(B_0+B_1x)`$ on either side. Permit four unrelated coefficients on each side, a relaxation of the actual continuous kernel. The square-norm contribution on A is the outer branch of the final level alone. On B the final coefficient is its opposite. After one common nonzero rescaling the combined square term is $`z(c+x^2)`$, where $`z=\mathbf1_{a=0}-\mathbf1_{b=0}`$ and $`c`$ is common to the cluster. The genuine $`c`$ is fixed by $`y`$, but allowing arbitrary $`c`$ strengthens the obstruction.

Order the twelve primes and take any of the thirteen splits, including an empty side. A prime at the knot can be assigned to either side by continuity. For each split delete any one row. If at most one prime survived, one of these systems of eleven equations would have a solution, with the four ground coefficients on either side and the constant $`c`$ as its nine unknowns. A row contains $`(a,ax,b,bx)`$ in the four columns for its side, zeros in the other four, and $`z`$ in the ninth column. Its right side is $`-zx^2`$.

On either side let $`m_d`$ count rows with projective direction $`d`$. The four ground columns have rank at most $`\min(4,\sum_d\min(2,m_d))`$, since all rows in a single direction span at most two dimensions. Add the two side bounds and one for the constant, calling the result $`U`$. For **every** deletion at every split, the augmented matrix has rank strictly greater than $`U`$ modulo 65537. A nonzero minor modulo this prime is nonzero over the rationals, since the denominators are powers of the selected primes and are invertible. The rational coefficient matrix has rank at most $`U`$, so the system is inconsistent even over the reals. The second checker uses modulus 1,000,000,007, eliminates the four columns on each side separately and obtains two independent residual constraints on the constant and the normalized quadratic coefficient. There are 1,786 certified clusters and 278,616 systems. This proves at least two prime coordinates survive throughout the intermediate range. [Primary rank and prime checks](phase_bridge.py) · [Independent elimination](phase_bridge_recheck.py) · [Verification output](evidence/check.txt)

### The thirteen phase cells

Only last starts $`y\le271/8`$ now need to be considered for isolation. Thus the needed prefix stops at exponent 5 for $`\alpha\le271/256`$ and at 4 otherwise. Its norm ceiling is 271. List all good norms through this ceiling, the endpoint $`271/256`$, and the phases $`m/2^j`$ between 1 and 2 at which a kernel branch changes. This partitions the interval at

```math
1,\ 271/256,\ 79/64,\ 81/64,\ 47/32,\ 97/64,\ 49/32,\ 223/128,\ 113/64,\ 239/128,\ 241/128,\ 31/16,\ 127/64,\ 2 .
```

In a fixed open cell every entry of the level-by-prime matrix is affine in $`\alpha`$ over the rationals. The standard vector at $`p`$ lies in the row space exactly when deleting column $`p`$ lowers the rank. This follows because the kernel dimension of the projection of the row space along that coordinate is either zero or one. At algebraic phases the rank test and all normalizations take place over $`K`$.

The primary replay checks the full matrix and each single-column deletion. At $`r+2`$ distinct interior rational points it certifies rank $`r`$. All minors of size $`r+1`$ have degree at most $`r+1`$, so they vanish identically. For the lower bound choose an $`r`$-row minor, compute its determinant at $`r+1`$ equally spaced points including both endpoints, and recover its Bernstein coefficients by exact interpolation. They have a common weak sign and at least one is nonzero. The Bernstein basis elements are strictly positive inside the cell, so the rank is at least $`r`$ everywhere inside it. The case $`r=0`$ needs only the upper bound. All 344 polynomial rank tests pass. Each included boundary point has a separate rational matrix calculation. The optional independent replay builds the character tables with Euler’s criterion, computes Smith invariant factors over $`\mathbf Q[t]`$ and counts their real roots in each cell. It also checks the endpoints directly. [Rational minor proof](phase_rank.py) · [Independent Smith calculation](phase_smith_recheck.py) · [Smith transcript](evidence/phase-smith.txt)

The single-field and pooled isolation lists follow. The ranks also show that every nonzero raw reading at a fixed algebraic phase is independent of all the others, including the independent constant 1. Through any exponent $`n\ge2`$ the pooled readings have rank

```math
2n-3+\mathbf1_{\alpha\gt49/32}+\mathbf1_{\alpha\gt31/16}.
```

Indeed both fields have a zero at exponent zero, B has another at exponent one, A has a zero at exponent one through phase $`31/16`$, and B has a zero at exponent two through phase $`49/32`$. The first norms and the strict support endpoints force exactly these initial zero candidates. The computed prefix rank equals the number of remaining readings. Every later linear relation would have a last nonzero pair at a start beyond $`271/8`$, where at least two prime coordinates are forced. Baker’s theorem then rules out the relation. Adding $`1=L_2`$ preserves independence.

### A pair of fields that can isolate a prime only together

Take $`\alpha=3/2`$, $`U=T_A(12)`$ and $`V=T_A(24)-T_B(24)`$. In the second expression the coordinates at 3, 7 and 79 cancel. The two remaining equations are

```math
U=\frac{17}{31}L_{31}-\frac1{47}L_{47},\qquad
V=\frac{14}{31}L_{31}-\frac{46}{47}L_{47}.
```

They give $`L_{31}=31(46U-V)/768`$ and $`L_{47}=47(14U-17V)/768`$. The one-field classification shows that neither field alone isolates even one prime at this phase. At $`\sqrt{7/3}`$ the pooled list is empty because $`(97/64)^2\lt7/3\lt(49/32)^2`$. [Exact rows](phase_rank.py) · [Kernel-checked inversion](DyadicPrimeSupport.lean)

### What changes at transcendental phases

Here we use the physical readings $`D_i(\alpha2^j)`$ with their factor $`1/\sqrt\alpha`$, and the combination coefficients are still restricted to real algebraic numbers. For each odd prime $`p`$, the set of all real phases where both fields’ combined span contains $`L_p`$ is countable and dense. The subset realized by a rational multiple of just **one** reading, $`D_A(256\alpha)`$, is itself countable, dense and entirely transcendental. Different primes give disjoint such subsets for this one reading, since rational proportionality of distinct prime logarithms is impossible. Consequently almost every real phase recovers no prime logarithm at all.

For countability, fix a finite list of exponents and algebraic coefficients. There are finitely many kernel branches as $`\alpha`$ runs from 1 to 2. On each open branch the physical combination is $`U/\sqrt\alpha+V\sqrt\alpha`$ with fixed real constants $`U,V`$. If it equals $`L_p\ne0`$, then $`Vt^2-L_pt+U=0`$ at $`t=\sqrt\alpha`$. This polynomial is nonzero and has at most two roots. Include the finitely many endpoints and take the countable union over finite lists, algebraic coefficients and primes. [Quadratic-root deduction](DyadicPrimeSupport.lean)

For density, on any nonempty open interval choose a smaller interval inside a single branch of $`D_A(256\alpha)`$. The ordinary-prime octave calculation supplies at least two active terminal prime columns at every height here. At any algebraic point in that smaller interval Baker’s theorem gives a nonzero reading. Thus the coefficients $`U,V`$ are not both zero, and the function $`U/\sqrt\alpha+V\sqrt\alpha`$ is nonconstant on every such interval. Its continuous image contains a nonempty open interval. For any $`p`$ this interval contains $`qL_p`$ with nonzero rational $`q`$. The intermediate value theorem supplies a phase, and division by $`q`$ recovers $`L_p`$. An algebraic phase cannot do this using the specified single reading, since it has at least two prime coordinates. This constructs the claimed dense set of transcendental offsets. The phase is part of the choice of observation.

## Minimum support on one field

Fix the field, let $`N_k`$ count its good rational primes with $`2^k\lt m_i(p)\lt2^{k+1}`$ and set $`y=2^n`$. Normalize the last coefficient, on $`D_i(y)`$, to one. As all earlier coefficients vary in $`K`$, the minimum support is $`N_{n+1}+N_{n+2}`$, with one subtracted exactly at $`n=3,4`$ on A and $`n=4,5`$ on B. For every $`n\ge8`$ the unique minimizing **value** is $`S_i(2y)/a`$, the sum of all changes through exponent $`n`$. Coefficients on initial zero readings can be chosen freely.

To see this, change coordinates from differences to levels. Since $`S_i(1)=0`$, the coefficient of $`S_i(2y)/a`$ is fixed to one and those of $`S_i(2^j)/a`$ for $`1\le j\le n`$ are arbitrary. Each level is a triangular function of $`1/m`$ with support $`[2^j,2^{j+2}]`$. It is positive on the interior, vanishes at the endpoints and takes the value one at the central knot, $`m=2^{j+1}`$, before algebraic scaling. The scaled coefficient at each such knot is therefore freely prescribable apart from the specified terminal one. The smallest knots carry fixed zero values where there are no observations. No good norm is a power of two. Distinct rational primes have distinct norms on a fixed field.

On $`(4y,8y)`$ only the fixed terminal level contributes. On $`(2y,4y)`$ the kernel is affine in $`1/m`$ and has fixed nonzero value at $`4y`$. It can lose at most one prime coefficient. Erasing that coefficient prescribes a nonzero value at $`2y`$. To save one prime overall all smaller sampled norms must vanish. Starting at $`2y`$ and proceeding left, a cell with no samples can absorb any incoming knot value and leave all earlier knots zero. One sample uniquely determines a nonzero outgoing value from a nonzero incoming one. Two distinct zeros force both endpoints zero. Thus a saving is possible exactly when $`N_{n+1}\gt0`$ and, moving back from cell $`n`$ past cells of size one, the first other size is zero. If the scan runs out of cells, the condition is vacuous. The initial empty cells below the first good norm make the fixed low-end boundary harmless. Every constructed knot is algebraic, since the reciprocal norms and interpolation coordinates are rational. [Kernel-checked propagation rule](DyadicPrimeSupport.lean) · [Direct linear-system check](support_recheck.py)

For $`0\le n\le8`$ the minima in order are

| Field | Minimum support at exponents 0 through 8 |
| --- | --- |
| A | 0, 0, 1, 2, 6, 6, 8, 17, 21 |
| B | 0, 0, 0, 1, 4, 5, 7, 16, 27 |

Every octave $`k\ge8`$ contains at least three good norms per field. The certificate gives three primes per octave for $`8\le k\le127`$. For the rest take $`X=2^{k/4}\ge2^{32}`$ and $`t=19/16`$, so $`t^4\lt2`$. Bennett, Martin, O’Bryant and Rechnitzer bound the error in $`\theta(x;8,3)`$ from $`x/4`$ by $`0.0005840x/\log x`$ for $`x\ge2261078657`$. As $`\log X\gt20`$, subtracting the endpoint bounds gives prime weight greater than $`X/24`$ in $`(X,tX]`$. Since $`\log(tX)\lt\sqrt X`$ here, this supplies more than $`\sqrt X/24\gt8`$ primes. At most five are excluded on a field and every remaining prime $`3\bmod8`$ has degree four. This also shows uniqueness of the minimizing value for $`n\ge8`$. A nonzero knot at $`2y`$ introduces at least two earlier coefficients while it can erase at most one later one. With this knot zero, all terminal coefficients remain and any new support strictly raises the count. [720 witnesses](evidence/prime-norm-witnesses.tsv) · [BMOR, Theorem 1.2 and Table 1](https://personal.math.ubc.ca/~bennett/BeMaObRe-IJM-2018.pdf)

## Explicit recovery and early rank

Set $`T_i(z)=B_0\sqrt z\,S_i(z)/(8a)`$. Telescoping expresses each such dyadic level as an algebraic combination of readings. There are three exact identities

```math
T_A(8)=L_{31}/31,\qquad T_B(16)=-15L_7/49,\qquad
T_A(16)=30L_{31}/31-17L_{47}/47-15L_7/49 .
```

Thus $`L_{47}=(47/17)(30T_A(8)+T_B(16)-T_A(16))`$. The one-field minima prove neither single field can isolate it. To check all earlier possibilities for the joint space, list the good norms below $`8\cdot2^n`$ on each field. In a column for $`p`$ the row for $`T_i(2^j)`$ is $`s_i(p)W(m_i(p)/2^j)`$, or zero at an omitted prime. Test whether adjoining its standard coordinate vector leaves the rational rank unchanged. A separate elimination modulo 65537 certifies matching nonzero minors, including every non-isolated coordinate. Rational matrices have the same rank over $`K`$. The nonempty rows through exponents 0 to 7 have ranks $`0,0,1,3,5,7,9,11`$ and see $`0,0,1,3,8,11,20,36`$ prime columns. The first isolated coordinate is 31 at exponent 2. At exponent 3 the three listed coordinates are isolated, and the list stays fixed through 7. [Exact row-space tests](pooled_rank.py) · [Expected ranks](evidence/check.txt)

## The finite bridge

Consider any linear combination of the changes whose latest nonzero pair of coefficients is $`\lambda=(\lambda_A,\lambda_B)`$ at exponent $`n`$, and write $`y=2^n`$. At an ordinary prime $`4y\lt p\lt8y`$, only latest changes can contribute. Its vector, up to the same positive scalar on both fields at this prime, is the signed pair $`e(p)`$ of degree-one coefficients. On field $`i`$, set $`e_i(p)=0`$ if that field omits $`p`$ or if $`p\not\equiv\pm1\bmod16`$. In all other cases $`e_i(p)=-(\chi_u(p)+\chi_v(p))/2`$, with $`u,v`$ the conductors of that field. In increasing order, the four listed trace exclusions have degree-one pairs $`(0,-1),(0,0),(1,0),(1,0)`$. For example, 343057911569 is omitted on A and splits on B. [Calculation for the four primes](support_recheck.py)

A nonzero direction can annihilate at most one of three distinct projective directions. Direct search through exponent 21 and labelled prime certificates through exponent 1580 give at least three directions in every such octave from exponent 8. At least two primes survive.

For $`1581\le n\le2899`$ the witnesses use the two different norms of the **same** rational prime. Four distinct primes satisfy $`4y\lt p^2\lt8y`$ with degree two on A and degree one on B, all in one ordinary-prime octave. Two have each sign on B. Four more witnesses exchange A and B. For $`\lambda_A\ne0`$ the degree-two contribution is the outer branch of the latest A reading, with a common expression $`C+D/p^2`$, $`D\ne0`$, after factoring out $`L_p`$. No earlier A level reaches this norm. Over the common ordinary-prime octave all B kernels are affine in $`1/p`$ and have the same sign factor, so their sum contributes $`\pm(u+v/p)`$ for some real $`u,v`$. If the two primes of each sign all vanished, each of $`Dx^2+vx+(C+u)`$ and $`Dx^2-vx+(C-u)`$ would have two distinct positive roots at $`x=1/p`$. Subtracting the equations at the first pair gives $`v=-D(x_1+x_2)`$. The second pair gives $`v=D(x_3+x_4)`$, a contradiction. Hence at least one prime from this group survives. The reverse group gives a different prime if $`\lambda_B\ne0`$. If only one last coefficient is nonzero, the three common fourth-power norms in the terminal octave give three survivors. Primes $`3\bmod8`$ have sign $`-1`$ on both fields. Here these primes exceed every omission. [Opposite-sign quadratic lemma](DyadicPrimeSupport.lean) · [Finite octave checks](dyadic_witnesses.py)

The 13,911 compressed entries come from the earlier 66,016-prime computation. Each reconstructs a labelled prime from its geometric bin, a coefficient offset and a base. Trial division handles the small entries. Each larger entry has $`p=1+k2^m`$, $`k`$ odd, $`0\lt k\lt2^m`$, and an exact certificate $`b^{(p-1)/2}\equiv-1\bmod p`$. For every prime divisor of $`p`$ the order of $`b`$ is a multiple of $`2^m`$, so that divisor is at least $`2^m+1\gt\sqrt p`$. This proves primality. The replay also checks the four quadratic residues, all distinctness requirements, the open endpoints on the relevant powers, and the common ordinary-prime octave. [Witness data](evidence/dyadic-proth-witnesses.tsv) · [Source edition](https://github.com/leecrossley/theta-host-split-inert-race/releases/tag/joint-frobenius-primes-v1.0.2)

## Ordinary primes in every direction

For every real $`X`$ with $`\log X\ge2000`$ and every real unit vector $`\lambda`$, more than $`X/(144\log(2X))`$ primes $`p\in(X,2X)`$ have $`\lambda\cdot e(p)\ne0`$. All exclusions are below $`X`$. The argument uses only pairwise products of conductor characters.

Put $`h=(\log2)/4`$ and let $`g`$ be the fourfold convolution of the probability density $`h^{-1}\mathbf1_{[0,h]}`$. It is nonnegative, supported on $`[0,\log2]`$, vanishes at the two endpoints and integrates to one. Young’s inequality gives $`\|g\|_\infty\le1/h\lt6`$. Differentiating one box gives two translates of the other three, hence $`\|g'\|_1\le2/h`$. Use the prime weight $`w(p)=g(\log(p/X))\log(p)/p`$. On each of the two classes $`\pm1\bmod16`$ partial integration against $`\theta(t;16,c)-t/8`$ gives error at most $`(1+2/h)/(160\log X)`$ from its mass $`1/8`$. Indeed the test derivative is $`(g'-g)/t^2`$, and the BMOR bound is $`t/(160\log t)`$ for $`t\ge8\cdot10^9`$. The base prime weight on their union consequently exceeds $`249/1000`$. [BMOR, Theorem 1.2](https://arxiv.org/abs/1802.00085)

For the four conductor characters in the order A then B, let $`v=(\lambda_A,\lambda_A,\lambda_B,\lambda_B)`$. Expanding $`(\lambda\cdot e(p))^2`$ gives the base indicator with coefficient $`1/2`$ and six pair products with coefficients $`v_jv_k/2`$ on that base. Its indicator is $`(1/4)\sum_{\ell=0}^3\psi^\ell(p)`$, where $`\psi`$ is the primitive even quartic character modulo 16. The sum of absolute coefficients of the resulting 24 twists is $`(1+4|\lambda_A\lambda_B|)/2\le3/2`$. The primitive inducing conductors for each pair $`u,v`$ are $`uv,8uv,16uv,16uv`$ and are all below $`Q=7000000`$. Each of the 24 characters is even and nonprincipal. At primes in the window the induced and primitive values agree.

### A zero-free strip and its real exception

McCurley’s zero-free theorem gives at most one exceptional zero, real and on a real character, in the region to the right of $`1-1/(9.645908801\log\max(q,q|t|,10))`$. We use the weaker denominator $`10\log(q(|t|+2))`$. Each real character among the 24 has the form $`\chi_u\chi_v`$ or its quadratic twist of conductor 8. Its $`L`$-function exceeds 0.46 on the real interval $`[0.99,1]`$. Here is the finite proof of that assertion. [McCurley, *Explicit zero-free regions for Dirichlet L-functions*](https://www.sciencedirect.com/science/article/pii/0022314X84900891) · [Twelve computations](zero_strip.py)

Write $`A(t)=\sum_{n\le t}\chi(n)`$, $`M=2^{16}`$, and compute $`B=\max|A(t)|`$ over a full period, which ends at zero. With $`A_+=\max(A,0)`$, $`A_-=\max(-A,0)`$, bound $`P,N`$ from the sums of $`A_+(n)/(n(n+1))`$ and $`A_-(n)/(n(n+1))`$ for $`1\le n\lt M`$, rounding each first summand down and each second one up in units of $`2^{-36}`$. For $`0.99\le s\le1`$, partial summation gives $`L(s,\chi)=s\int_1^\infty A(t)t^{-s-1}\,dt`$. Since $`M^{0.01}\lt9/8`$, this implies

```math
L(s,\chi)\ \ge\ 0.99P-(9/8)N-9B/(8M).
```

The smallest of the twelve directed lower bounds exceeds 0.468. A second computation evaluates the series at 1, bounds its derivative through the differentiated integral using $`1+\log t`$ and obtains $`L(s,\chi)\gt1/4`$ throughout all twelve intervals. The real exceptional zero therefore cannot enter the common strip used below. [Independent computation and bounds](zero_strip_recheck.py)

### The smoothed character sum

For a primitive character just described set $`K(z)=((e^{hz}-1)/(hz))^4`$, continued at zero. Mellin inversion and the pole expansion of the logarithmic derivative give

```math
\sum_{m\ge1}\frac{\Lambda(m)\chi(m)}m g(\log(m/X))
=-\sum_\rho X^{\rho-1}K(\rho-1)
-\sum_{k\ge0}X^{-1-2k}K(-1-2k).
```

Here $`\rho`$ ranges over nontrivial zeros with multiplicity. To obtain the identity, integrate $`-(L'/L)(s,\chi)X^{s-1}K(s-1)`$ to the right of 1 and use the absolutely convergent Dirichlet series there. Shift the line to successively more negative odd integers. The residues at zeros are negative and the trivial zeros are $`0,-2,-4,\ldots`$. On bounded vertical strips one may take horizontal ordinates tending to infinity away from zeros, with logarithmic derivative $`O_q(\log^2(2+|t|))`$, by the local partial fraction expansion and the unit-height zero count. The transform decays to fourth order in the ordinate. On each negative odd line the functional equation and Stirling’s formula bound the logarithmic derivative by $`O(\log(q(2+|s|)))`$ uniformly, since the line stays separated from the gamma poles. Its integral vanishes in the limit by the factor $`X^{\Re s-1}`$. Absolute convergence of the displayed zero sum follows also from the bound below. [Functional equation and zero count](https://arxiv.org/abs/2005.02989) · [Mellin inversion and explicit formulas](https://doi.org/10.1017/CBO9780511618314)

For a zero $`\rho=\beta+i\gamma`$ the integral transform of the density gives $`|K(\rho-1)|\le1`$ and its explicit formula gives $`|K(\rho-1)|\le(12/|\gamma|)^4`$ when $`\gamma\ne0`$. Write $`N(t)`$ for the number of nontrivial zeros with $`|\gamma|\le t`$. BMOR’s Corollary 1.2 bounds the error about $`(t/\pi)\log(qt/(2\pi e))`$ by $`0.247\log(qt)+6.894`$ for $`t\ge5/7`$. In particular $`N(t)\le(t/2)\log(Qt)`$ for $`t\ge12`$, using $`\pi\gt3`$ and $`\log(Qt)\gt16`$. Integration of the envelope, including zeros below 12, gives

```math
\sum_\rho |K(\rho-1)|\le\int_{12}^\infty 4\cdot12^4N(t)t^{-5}\,dt
 \le 8(\log(12Q)+1/3)\lt155 .
```

At cutoff $`H=512`$ the tail satisfies

```math
\sum_{|\gamma|\gt H}|K(\rho-1)|\le\int_H^\infty 4\cdot12^4N(t)t^{-5}\,dt
 \le\frac{2\cdot12^4}{3H^3}(\log(QH)+1/3)\lt1/400 .
```

For $`|\gamma|\le H`$ the strip and $`\log(Q(H+2))\lt22.5`$ give $`X^{\beta-1}\lt1/7000`$. The common strip line here lies to the right of 0.99, where the real exception has been excluded. The trivial-zero series costs at most $`1/(X(1-X^{-2}))`$. To pass to ordinary primes, the proper prime powers in the support cost at most $`\sqrt{2/X}(\log(2X))^2/(h\log2)`$. There are at most $`\sqrt{2X}`$ possible bases and at most $`\log(2X)/\log2`$ exponents per base. Both errors are less than $`10^{-6}`$ when $`\log X\ge2000`$. Substituting $`t=\log X`$ gives decreasing bounds, and $`e^{1000}\gt2^{1000}`$ suffices at the endpoint. Thus each of the 24 ordinary-prime sums $`\sum_p\chi(p)w(p)`$ has absolute value below $`155/7000+1/400+2\cdot10^{-6}\lt1/40`$. The second replay uses cutoff 640 and obtains $`1/50`$. [Zero count, Corollary 1.2](https://arxiv.org/abs/2005.02989) · [Directed rational constants](zero_strip.py) · [Second cutoff](zero_strip_recheck.py)

Combining the expansion with the base prime mass gives $`\sum_pw(p)(\lambda\cdot e(p))^2\gt249/2000-3/80\gt1/12`$. Each contributing prime supplies at most $`12\log(2X)/X`$, since $`|\lambda\cdot e(p)|^2\le2`$. This proves the asserted strict prime count.

With $`X=4\cdot2^n`$ the estimate applies to every $`n\ge2900`$ and gives at least two surviving ordinary primes in the open terminal octave. Together with the finite bridge it proves that any combination with a nonzero last pair at exponent $`n\ge8`$ has at least two nonzero prime coordinates. The early ranks now give the complete list of isolated primes. They also show that through exponent 7 every nonzero reading is independent. A linear relation involving a later reading would have a last pair and hence a nonzero prime coordinate. Baker’s theorem proves independence for the numerical readings, including the independent constant 1. In particular the two one-field spans intersect only at zero. There are exactly five zero readings, at exponents 0 and 1 on A and 0, 1 and 2 on B. Thus the rank of the readings through $`n\ge3`$ is exactly $`2n-3`$.

## How many ordinary primes must remain?

Let $`\lambda\ne0`$ be the prescribed terminal coefficient pair, $`y=2^n`$, and let $`M_\lambda(y)`$ be the minimum prime support over all earlier algebraic coefficients, with $`\lambda\in K^2`$. Write $`P_\lambda(y)`$ for the number of good ordinary primes $`2y\lt p\lt8y`$ with $`\lambda\cdot e(p)\ne0`$, using primes good for both fields. For sufficiently large $`y`$ all omitted rational primes are below $`2y`$. We have

```math
P_\lambda(y)-4\ \le\ M_\lambda(y)\ \le\ P_\lambda(y)+\pi(\sqrt{8y}).
```

The outer octave forces every active ordinary prime. In the inner octave the terminal and preceding levels combine as an affine function of $`1/p`$ for each projective direction of $`e(p)`$. For any active direction its value at the right knot $`4y`$ is nonzero. There are four projective directions, so at most four primes in the inner octave can be lost. A higher norm of an ordinary prime in these octaves exceeds $`8y`$ for $`y\gt2`$. For the upper bound, use constant coefficients $`\lambda`$ in every step. The levels telescope, leaving only the two final levels. Any additional supported rational prime has residue degree at least two and therefore lies below $`\sqrt{8y}`$.

The prime number theorem in fixed progressions and the Chinese remainder theorem make the four quadratic signs independent and equiprobable on primes $`\pm1\bmod16`$, a set of density $`1/4`$. Each race coordinate there is $`-1,0,1`$ with multiplicities $`1,2,1`$. Conditional on the base, the fractions with nonzero projection are $`1/2`$ for an axial direction, $`5/8`$ for equal absolute components, and $`3/4`$ for any other direction with both components nonzero. Multiplying by the base density and the interval length $`6y`$ gives the sharp equivalents

```math
M_\lambda(y)\sim c_\lambda\,y/\log y .
```

| Terminal direction | $`c_\lambda`$ |
| --- | --- |
| One nonzero component | $`3/4`$ |
| Two nonzero components of equal absolute value | $`15/16`$ |
| Two nonzero components of unequal absolute value | $`9/8`$ |

There are finitely many active sets of directions, so the asymptotic is uniform even when $`\lambda`$ varies with $`n`$. This counts distinct rational primes, including a single shared coordinate when the two fields give different powers of the same prime. [Independent multiplicity checks](support_recheck.py) · [Prime number theorem in progressions](https://kskedlaya.org/ant/chapter-11.html)
