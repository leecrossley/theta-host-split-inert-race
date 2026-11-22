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

Thus $`L_{47}=(47/17)(30T_A(8)+T_B(16)-T_A(16))`$. The one-field minima prove neither single field can isolate it. To check all earlier possibilities for the joint space, list the good norms below $`8\cdot2^n`$ on each field. In a column for $`p`$ the row for $`T_i(2^j)`$ is $`s_i(p)W(m_i(p)/2^j)`$, or zero at an omitted prime. Test whether adjoining its standard coordinate vector leaves the rational rank unchanged. A separate elimination modulo 65537 certifies matching nonzero minors, including all excluded coordinates. Rational matrices have the same rank over $`K`$. The nonempty rows through exponents 0 to 7 have ranks $`0,0,1,3,5,7,9,11`$ and see $`0,0,1,3,8,11,20,36`$ prime columns. The first isolated coordinate is 31 at exponent 2. At exponent 3 the three listed coordinates are isolated, and the list stays fixed through 7. [Exact row-space tests](pooled_rank.py) · [Expected ranks](evidence/check.txt)

## The finite bridge

Consider any linear combination of the changes whose latest nonzero pair of coefficients is $`\lambda=(\lambda_A,\lambda_B)`$ at exponent $`n`$, and write $`y=2^n`$. At an ordinary prime $`4y\lt p\lt8y`$, only latest changes can contribute. Its vector, up to the same positive scalar on both fields at this prime, is the signed pair $`e(p)`$ of degree-one coefficients. On field $`i`$, set $`e_i(p)=0`$ if that field omits $`p`$ or if $`p\not\equiv\pm1\bmod16`$. In all other cases $`e_i(p)=-(\chi_u(p)+\chi_v(p))/2`$, with $`u,v`$ the conductors of that field. In increasing order, the four listed trace exclusions have degree-one pairs $`(0,-1),(0,0),(1,0),(1,0)`$. For example, 343057911569 is omitted on A and splits on B. [Calculation for the four primes](support_recheck.py)

A nonzero direction can annihilate at most one of three distinct projective directions. Direct search through exponent 21 and labelled prime certificates through exponent 1580 give at least three directions in every such octave from exponent 8. At least two primes survive.

For $`1581\le n\le2899`$ the witnesses use the two different norms of the **same** rational prime. Four distinct primes satisfy $`4y\lt p^2\lt8y`$ with degree two on A and degree one on B, all in one ordinary-prime octave. Two have each sign on B. Four more witnesses exchange A and B. For $`\lambda_A\ne0`$ the degree-two contribution is the outer branch of the latest A reading, with a common expression $`C+D/p^2`$, $`D\ne0`$, after factoring out $`L_p`$. No earlier A level reaches this norm. Over the common ordinary-prime octave all B kernels are affine in $`1/p`$ and have the same sign factor, so their sum contributes $`\pm(u+v/p)`$ for some real $`u,v`$. If the two primes of each sign all vanished, each of $`Dx^2+vx+(C+u)`$ and $`Dx^2-vx+(C-u)`$ would have two distinct positive roots at $`x=1/p`$. Subtracting the equations at the first pair gives $`v=-D(x_1+x_2)`$. The second pair gives $`v=D(x_3+x_4)`$, a contradiction. Hence at least one prime from this group survives. The reverse group gives a different prime if $`\lambda_B\ne0`$. If only one last coefficient is nonzero, the three common fourth-power norms in the terminal octave give three survivors. Primes $`3\bmod8`$ have sign $`-1`$ on both fields. Here these primes exceed every omission. [Opposite-sign quadratic lemma](DyadicPrimeSupport.lean) · [Finite octave checks](dyadic_witnesses.py)

The 13,911 compressed entries come from the earlier 66,016-prime computation. Each reconstructs a labelled prime from its geometric bin, a coefficient offset and a base. Trial division handles the small entries. Each larger entry has $`p=1+k2^m`$, $`k`$ odd, $`0\lt k\lt2^m`$, and an exact certificate $`b^{(p-1)/2}\equiv-1\bmod p`$. For every prime divisor of $`p`$ the order of $`b`$ is a multiple of $`2^m`$, so that divisor is at least $`2^m+1\gt\sqrt p`$. This proves primality. The replay also checks the four quadratic residues, all distinctness requirements, the open endpoints on the relevant powers, and the common ordinary-prime octave. [Witness data](evidence/dyadic-proth-witnesses.tsv) · [Source edition](https://github.com/leecrossley/theta-host-split-inert-race/releases/tag/joint-frobenius-primes-v1.0.2)

## Ordinary primes in every direction

For every real $`X`$ with $`\log X\ge2000`$ and every real unit vector $`\lambda`$, more than $`X/(144\log(2X))`$ primes $`p\in(X,2X)`$ have $`\lambda\cdot e(p)\ne0`$. All exclusions are below $`X`$. The argument uses only pairwise products of conductor characters.

Put $`h=(\log2)/4`$ and let $`g`$ be the fourfold convolution of the probability density $`h^{-1}\mathbf1_{[0,h]}`$. It is nonnegative, supported on $`[0,\log2]`$, vanishes at the two endpoints and integrates to one. Young’s inequality gives $`\|g\|_\infty\le1/h\lt6`$. Differentiating one box gives two translates of the other three, hence $`\|g'\|_1\le2/h`$. Use the prime weight $`w(p)=g(\log(p/X))\log(p)/p`$. On each of the two classes $`\pm1\bmod16`$ partial integration against $`\theta(t;16,c)-t/8`$ gives error at most $`(1+2/h)/(160\log X)`$ from its mass $`1/8`$. Indeed the test derivative is $`(g'-g)/t^2`$, and the BMOR bound is $`t/(160\log t)`$. The base prime weight on their union consequently exceeds $`249/1000`$. [BMOR, Theorem 1.2](https://arxiv.org/abs/1802.00085)

For the four conductor characters in the order A then B, let $`v=(\lambda_A,\lambda_A,\lambda_B,\lambda_B)`$. Expanding $`(\lambda\cdot e(p))^2`$ gives the base indicator with coefficient $`1/2`$ and six pair products with coefficients $`v_jv_k/2`$ on that base. Its indicator is $`(1/4)\sum_{\ell=0}^3\psi^\ell(p)`$, where $`\psi`$ is the primitive even quartic character modulo 16. The sum of absolute coefficients of the resulting 24 twists is $`(1+4|\lambda_A\lambda_B|)/2\le3/2`$. The primitive inducing conductors for each pair $`u,v`$ are $`uv,8uv,16uv,16uv`$ and are all below $`Q=7000000`$. Each of the 24 characters is even and nonprincipal. At primes in the window the induced and primitive values agree.

### A zero-free strip and its real exception

McCurley’s zero-free theorem gives at most one exceptional zero, real and on a real character, in the region to the right of $`1-1/(9.645908801\log\max(q,q|t|,10))`$. We use the weaker denominator $`10\log(q(|t|+2))`$. Each real character among the 24 has the form $`\chi_u\chi_v`$ or its quadratic twist of conductor 8. Its value on the real interval $`[0.99,1]`$ exceeds 0.46. Here is the finite proof of that assertion. [McCurley, *Explicit zero-free regions for Dirichlet L-functions*](https://www.sciencedirect.com/science/article/pii/0022314X84900891) · [Twelve computations](zero_strip.py)

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
