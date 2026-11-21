# Joint levels and moves of the two theta races

## The two observables

Let $`A=(547,727)`$ and $`B=(571,599)`$. In each field of the [statement](README.md) put $`\epsilon_v=1`$ for split and $`\epsilon_v=-1`$ for inert. Omit the primes over 2 and the two conductors. Additional fixed finite omissions have no effect on the tail. For the arithmetic trace interpretation omit the exceptional set of the chosen admissible datum. The first-step parameter recursion subtracts the previous twisted trace and divides by the trace of the added segment. The fibre product of these normalized jumps equals $`\epsilon_v`$. Both hosts have explicit admissible data. [First trace extraction](https://github.com/leecrossley/theta-prime-trace/blob/7a7909d2dc7bb9d8d0ec7b2dc8c78e70f93782d7/ThetaPrimeTrace.lean#L68) · [547/727 host](https://github.com/leecrossley/theta-race-short-shifts/blob/483c31f2df2132ed2165620ea55e20e969228bfb/RaceTrace.lean#L83) · [571/599 host](https://github.com/leecrossley/theta-host-split-inert-race/blob/8fca7cad37185e3d7951f0dffff5e2442d1dfa24/HostAdmissible.lean)

For $`b\gt 1`$ put $`W_b(t)=\max(0,\min(b(1-1/t),b^2/t-1))`$ and $`B_b=2\sqrt b(\sqrt b-1)^2`$. The normalized score is

```math
S_{i,b}(y)=-\frac{1}{B_b\sqrt y}\sum_v\epsilon_v\log(Nv)\,W_b(Nv/y).
```

The score is continuous for $`y\gt 0`$. Every compact interval involves finitely many prime ideals. Write $`a=\log 2`$, fix independent widths $`b_A,b_B\gt 1`$ and put $`x=\log y`$. The four coordinates are the midpoint levels $`L_i=S_{i,b_i}(e^{x+a/2})`$ and the normalized endpoint differences $`D_i=(S_{i,b_i}(e^{x+a})-S_{i,b_i}(e^x))/a`$.

## Residues for a real projection

For $`i=(p,q)`$ the eight primitive characters of $`\zeta_{E_i}/\zeta_{F_i}`$ are the two quadratic characters of conductors $`p,q`$, each twisted by all four powers of the even quartic character modulo 16. Take each twist at primitive conductor. Write their product of L-functions as $`\mathcal L_i`$. Expanding the logarithmic derivatives at the prime powers gives the following continuation from $`\Re s\gt 1`$ to $`\Re s\gt 1/3`$,

```math
-\sum_v\epsilon_v\log(Nv)(Nv)^{-s}
=M_i(s)=\mathcal L_i'(s)/\mathcal L_i(s)-\zeta_{F_i}'(2s)/\zeta_{F_i}(2s)+R_i(s).
```

The square correction restores the unsigned prime squares. The signed powers of order at least three and the higher terms of the correction are absolutely convergent in this half-plane. The omitted finite factors are regular there. The real-axis certificate below excludes real zeros from $`[1/2,1]`$ for all sixteen factors. Euler products handle $`s\gt 1`$. Thus the real positive axis of $`M_i(1/2+z)`$ is regular, and the residue at $`z=0`$ is $`1/2`$. At $`z=i\gamma\ne0`$ the prime ideal theorem keeps the base zeta nonzero on $`\Re s=1`$. The residue equals the nonnegative total zero multiplicity $`m_i(\gamma)`$ of the eight factors at $`1/2+i\gamma`$. [Euler products and continuation](https://github.com/leecrossley/theta-race-short-shifts/blob/483c31f2df2132ed2165620ea55e20e969228bfb/PROOF.md#euler-factors-and-the-kernel) · [Real-axis cover](real_axis.py)

In centred log coordinates $`Q_{i,b}(x)=S_{i,b}(e^x/b)`$, direct integration of the weight gives the bilateral Laplace transform $`G_b(z)M_i(1/2+z)`$, with an entire $`G_b`$. On the imaginary axis,

```math
G_b(i\gamma)=
\frac{b+1-2\sqrt b\cos(\gamma\log b)}
 {2(\sqrt b-1)^2(1/4+\gamma^2)}
\ \ge\ \frac{1}{2(1/4+\gamma^2)},\qquad G_b(0)=2.
```

Take a real projection $`f=\sum_i(u_i L_i+v_i D_i)`$. Changing from the bilateral transform to the transform on $`[0,\infty)`$ adds an entire function because the scores vanish sufficiently far to the left. Translation of $`Q_{i,b}`$ gives the residue from host $`i`$ at $`i\gamma`$ as

```math
m_i(\gamma)e^{i\gamma(\log b_i+a/2)}G_{b_i}(i\gamma)
   \bigl(u_i+i\,k(\gamma)v_i\bigr),\qquad
 k(\gamma)=2\sin(a\gamma/2)/a .
```

Conjugate symmetry follows from the full set of twists. At zero the projected residue is $`u_A+u_B`$. The phase includes the independent width of each field. Phase alignment between the two fields is unnecessary.

## Nine frequencies on one height axis

Choose an ordinate inside each of the nine disjoint brackets where the appropriate Hardy function changes sign. The first five rows belong to $`A`$ and the remaining four to $`B`$. Bracket endpoints in the [data](evidence/brackets.tsv) are integers in units of $`10^{-14}`$. For each signed combination with coefficients $`0,\pm1`$ and at least two nonzero terms, the entire resulting interval misses zero. Its absolute value excludes all sixteen spectra on the critical line. Up to simultaneous sign reversal there are $`(3^9-1-18)/2=9832`$ combinations. An additional $`9\cdot8=72`$ tests exclude each singleton from the opposite host. Coincidences or higher multiplicities within its own host only increase its residue amplitude. These assertions cover every choice of a zero in each open bracket.

Fix a starting log height $`T`$. If $`f`$ is unbounded above on this tail, the projected lower bound on its supremum is immediate. Otherwise take a constant $`C`$ above it. The transform of the nonnegative tail $`C-f`$ is regular at every positive real point by the continuation above. The original transform converges at sufficiently large real part, since the prime sums have at most exponential growth in $`x`$. Landau's theorem on the abscissa of convergence gives absolute convergence of the tail transform for every $`\sigma\gt 0`$. Indeed, a finite abscissa would be singular there, as follows by monotone convergence of the exponential Taylor series at a point slightly to its left. This proves the required convergence for $`f`$ and its imaginary translates. No global zero location enters this step. [Diamond's oscillation theorem](https://personal.math.ubc.ca/~gerg/teaching/592-Fall2018/papers/1971.Diamond.pdf)

Let $`P(x)=\prod_j(1+\cos(\gamma_j x+\theta_j))`$. It is nonnegative, and the certified signed combinations show that its constant coefficient is one. In the Abel limit $`\lim_{\sigma\downarrow0}\sigma\int_0^\infty f(x)P(x)e^{-\sigma x}\,dx`$, every mixed frequency contributes zero by regularity at that frequency. Each singleton contributes the real part of its projected residue multiplied by a phase of absolute value one. Choose each $`\theta_j`$ to make this its absolute value, with arbitrary phase for a zero projection. The Abel mass of $`P`$ alone tends to one. A bounded initial interval has limiting mass zero. With $`J_i`$ the rows of host $`i`$, this gives on every tail

```math
\sup_{x\gt T} f(x)\ \ge\ u_A+u_B+
 \sum_{i=A,B}\sum_{j\in J_i}G_{b_i}(i\gamma_j)
              \sqrt{u_i^2+k(\gamma_j)^2 v_i^2}.
```

The inequality uses $`m_i(\gamma_j)\ge1`$. Applying it to the negative projection gives the opposite sign. Every strict consequence for the supremum holds at an actual height by the definition of supremum.

## Common starts and the number of readings

Put $`g_j=G_{b_i}(i\gamma_j)`$ in $`J_i`$, $`k_j=k(\gamma_j)`$, $`a_i=\sum_{J_i}g_j`$ and $`d_i=\sum_{J_i}g_j k_j`$. All are positive since the nine ordinates are between 0 and 0.58. The triangle inequality in the Euclidean plane gives

```math
\sum_{j\in J_i}g_j\sqrt{u_i^2+k_j^2 v_i^2}
\ \ge\ \sqrt{(a_i u_i)^2+(d_i v_i)^2}.
```

The right side is the maximum of $`u_i X+v_i Y`$ on the ellipse centred at the origin with semiaxes $`a_i,d_i`$. Thus every real linear functional has at least as large a supremum on the actual four-coordinate curve as on the product of these two ellipses translated to $`(1,0,1,0)`$. The finite-dimensional separating hyperplane theorem puts the product in the closed convex hull of the tail. The printed axes are strictly smaller in every coordinate, so their whole closed product lies in the interior of that hull. Each such interior point belongs to the ordinary convex hull. For example, enclose it by a small simplex strictly inside the closed hull, then approximate the simplex vertices by points in the ordinary convex hull. Strict interior containment persists.

Carathéodory's theorem now gives at most five actual starts in four dimensions. To reduce the count by one, consider any interior representation by five affinely independent vertices. Move the first continuously along the actual curve from its start to the second start. At the first failure of interior membership in the moving simplex, compactness of the barycentric weights gives a representation there. Either a weight vanishes or the simplex is affinely dependent. In the latter case eliminate a vertex by affine dependence. Interior membership of a nondegenerate simplex is open in its vertices, and at the end two vertices coincide. The whole path stays above $`T`$. This is the connected-set refinement of Carathéodory's theorem. It gives four common starts with one set of weights. Projection onto the levels alone or the moves alone gives a continuous planar curve and the same proof gives two common starts. [Fenchel's connected-set bound](https://www.renyi.hu/~barany/cikkek/131.pdf)

At $`b_A=b_B=2`$ the two unrounded sums $`a_A+a_B`$ and $`d_A+d_B`$ exceed 17.9106 and 5.6552. Apply the projected bounds directly to $`L_A-L_B`$ and $`D_A-D_B`$, and to their negatives. The two levels have the same constant residue, so it cancels. This proves the strict single-height assertions for the differences. The table in the [statement](README.md) rounds the individual axes down separately. For arbitrary widths the common lower bound for $`G_b`$ above gives the last two columns. [Interval bounds for the axes and sums](amplitudes.py) · [Support inequalities in Lean](JointOscillation.lean)

## Reproducing the numerical certificate

The first computation verifies primitive conductors and opposite Hardy signs. For a mixed or opposite-host interval it interpolates the complex L-function at twelve consecutive nodes of the grid of spacing $`1/128`$. With the target in the central cell the product of the twelve node distances is at most $`(6!)^2/128^{12}`$. In the complex ordinate disc of radius $`1/4`$ about each required real point, $`\Re s\ge1/4`$ and $`|s|\lt 6`$. Periodic partial sums bounded by $`q`$ in the integral $`L(s)=s\int_1^\infty C(t)t^{-s-1}\,dt`$ give $`|L(s)|\le24q`$. Cauchy's estimate and real interpolation applied separately to the two components bound the complex remainder by $`48q(6!)^2/32^{12}`$. Forward differences give the interpolation polynomial. A direct interval evaluation handles any unresolved case. [FLINT evaluation](zeros.py) · [FLINT documentation](https://flintlib.org/doc/acb_dirichlet.html)

The second computation constructs each character by modular exponentiation and an explicit quartic table. After the first $`q`$ terms, the remaining series can be expanded around the half-integer in each residue block. The coefficient of order $`k\ge1`$ is

```math
q^{-s}\frac{(-1)^k(s)_k}{k!}\zeta(s+k,3/2)
         \sum_{a=1}^q\chi(a)(a/q-1/2)^k .
```

The constant term vanishes. For $`\Re s=1/2`$ and $`|s|\lt 4`$ the rising factorial divided by $`k!`$ is at most $`(k+3)^3`$ in absolute value. Integral comparison bounds the Hurwitz factor by $`4(3/2)^{-k}`$ and the finite moment by $`q/2^k`$. After $`K=54`$ terms, summing the geometric majorant costs at most $`8\sqrt q(K+4)^3/3^K`$. The ratio of successive terms past the cutoff is at most $`((K+5)/(K+4))^3/3\lt 1/2`$. Normalizing the completed function by the conjugate of its central value makes its values real on the critical line by the functional equation. The signs enclose the same nine crossings. [Residue classes and Hurwitz tail](residues.py)

This second evaluator samples the 26 Chebyshev roots in each cell of length $`1/4`$ and encloses the interpolating polynomial by the Clenshaw recurrence. On that cell the monic node polynomial has absolute value bounded by $`2(1/16)^{26}`$. The same analytic-disc bound gives complex interpolation error at most $`96q\,4^{-26}`$. A separate rational enumeration generates the signed intervals. The real-axis check constructs the characters on every residue and covers $`[1/2,1]`$ by derivative-bounded cells. It uses $`|C(t)|\le\min(t,M)`$ from one full period, whence $`|L'(s)|\le4\sqrt M+4\sqrt M\log M+4`$ for real $`1/2\le s\le1`$, by splitting the integral at $`M`$. [Independent interpolation](independent.py) · [Real segments](real_axis.py)

For the four finite examples the residue degree $`f`$ of an unramified rational prime $`\ell`$ is its Frobenius order modulo complex conjugation. All $`8/f`$ base ideals have norm $`\ell^f`$ and sign $`(\ell/p)^f`$ in Legendre-symbol notation. An independent local computation factors the real quartic, the real quadratic and the imaginary quadratic defining polynomials over the residue field. The common degree in the real compositum is the least common multiple of the first two degrees. Adjoining the third determines the sign. The resulting coefficients agree prime by prime through $`8\cdot658`$. The height weights are exact rational numbers before multiplication by logarithms. These four finite counts certify all four sign patterns at the specified starts. [Prime counts](finite_examples.py) · [Full verification transcript](evidence/check.txt)
