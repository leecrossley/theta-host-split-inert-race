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

## The finite common windows

Fix both widths at two and retain the prime and score conventions of the two explicit theta hosts. In field A omit primes over 2, 547, 727, 343057911569 and 343070636977. In field B omit primes over 2, 571, 599, 413718278401 and 603339156001. Write $`r=\sqrt2`$, $`a=\log2`$, $`B_0=2r(r-1)^2`$, and now write $`x=\log(2r y)`$ for the log centre of a doubling with start $`y`$. Then the midpoint score is $`L_i=S_i(ry)`$ and the normalized move is $`D_i=(S_i(2y)-S_i(y))/a`$. A completion occurs at $`2y=e^x/r`$. [547/727 trace](https://github.com/leecrossley/theta-race-short-shifts/blob/483c31f2df2132ed2165620ea55e20e969228bfb/RaceTrace.lean#L83) · [571/599 trace](https://github.com/leecrossley/theta-host-split-inert-race/blob/8fca7cad37185e3d7951f0dffff5e2442d1dfa24/HostAdmissible.lean)

**Theorem.** At each $`1\le H\le e^{2000}`$, the convex hull of actual common observations with $`H\lt y`$ and $`2y\lt2^{18}H`$ contains the closed square $`[-4,4]^2`$ of normalized moves in its interior. With $`2y\lt2^{80}H`$ the corresponding four-dimensional hull contains the closed box $`|L_i-1|\le2`$, $`|D_i|\le4`$ in its interior.

### Prime-power and zero errors

For each field, replace the prime-only level by the fully signed odd-factor prime-power level $`P_i`$. The local Euler factors above give the eight odd primitive Dirichlet characters $`\mathcal A_i`$. The comparison score is

```math
P_i(y)=-\frac{1}{B_0\sqrt y}
\sum_{p\ {\rm prime}}\sum_{k\ge1}\log p\,
   \sum_{\chi\in\mathcal A_i}\chi(p)^k W_2(p^k/y).
```

Use the midpoint or take the normalized dyadic difference as appropriate. For a nontrivial zero with $`z=\rho-1/2`$, the corresponding term is $`G(z)e^{zx}`$ for the midpoint and $`J(z)e^{zx}`$ for the move. Here $`G=G_2`$ as in the original result and $`J(z)=2G(z)\sinh(az/2)/a`$. These formulas follow on integrating the two tent branches and applying the logarithmic derivatives, including the odd trivial zeros. No principal factor occurs. The strip mass bound for each family is $`\sum_\rho(1+(\Im\rho)^2)^{-1}\lt1000`$, and both strip kernels are bounded by $`1000/(1+(\Im z)^2)`$. Indeed the completed exterior derivative estimate uses at most eight primitive factors of conductor at most 11632 in either field. The level kernel has mass 2 and support within $`[-a,a]`$, and its rational expression gives the same strip bound used for the move. Both kernels have total variation at most 6 and support in $`[-2,2]`$. [Single-field derivation and mass estimate](https://github.com/leecrossley/theta-host-split-inert-race/blob/theta-race-waiting-bound-v1.1.0/results/theta-race-waiting-bound/PROOF.md#total-mass-of-the-omitted-zeros)

Here are uniform physical bounds for replacing prime powers by primes at starts above $`Y\ge10^9`$. Put $`Z=\sqrt Y`$ and $`B=B_0`$ temporarily. A degree-one prime splits in the real quartic field and lies in one of the classes $`\pm1\bmod16`$. For both classes and every $`u\ge Z`$, Ramaré and Rumely give $`|\vartheta(u;16,l)-u/8|\le\epsilon u/8`$. When $`Z\ge10^{10}`$ take $`\epsilon=0.009`$ by Table 1. Otherwise take the larger of this and $`8(2.072)/\sqrt Z`$, using their Theorem 2 for $`u\le10^{10}`$. [Ramaré and Rumely, Theorems 1 and 2](https://ramare-olivier.github.io/Maths/rumely.pdf)

The square correction to a level is nonnegative and bounded by

```math
Q=2+8\epsilon/B+d,\qquad d=8.16\sqrt{2/Z}/B .
```

There are at most eight degree-one ideals per rational prime. Partial summation against $`W_2(v^2)`$ on $`1\le v\le2`$ bounds their contribution by $`2+8\epsilon/B`$. The ordinary theta bound $`\vartheta(t)\le1.02t`$ controls the higher-degree ideals by $`d`$. For the square change use the positive or negative part of $`r^{-1}W_2(v^2/2)-W_2(v^2)`$. Both have area $`A=6r-2(1+r)\sqrt{2r}`$, and weighted variation at most $`2+r`$. The change costs at most $`Q_D=(2A+2\epsilon(2+r))/B+d`$ before division by $`a`$. Bounding higher powers and the five omitted primes at either scale gives

```math
E=\frac{8.16}{B\sqrt Y}\left((4Y)^{1/3}+(4Y)^{1/4}
       +\frac{\log(4Y)}{a}(4Y)^{1/5}\right),\qquad
E_b=\frac{8(3\log2+\sum_{p\in{\cal B}_i,\ p\ne2}\log p)}{B\sqrt Y}.
```

Here $`\mathcal B_i`$ is the five-prime omission set. All these bounds decrease, or stay constant, as the floor increases. The centred level therefore differs from the signed prime-power level by at most $`\max(1,Q-1)+E+E_b`$. The normalized move differs from its prime-power counterpart by at most $`(Q_D+2E+2E_b)/a`$. At $`Y=e^{100}`$ these two costs are below 1.151 and 2.353 in both fields. The finite calculation uses each row's own floor. The finitely many omitted prime powers are charged at every scale. [Derivation of the signed square bound](https://github.com/leecrossley/theta-host-split-inert-race/blob/theta-race-waiting-bound-v1.1.0/results/theta-race-waiting-bound/PROOF.md#cancellation-in-the-prime-squares) · [Arithmetic calculation](deadline_arithmetic.py)

Platt's finite theorem puts every odd-factor zero through height 8590 on the line because $`8590\cdot11632\lt10^8`$. We include 1,584 positive sign brackets for A and 1,634 for B, all below 125. At each $`3/2+iv`$ for $`v=0,2,\ldots,48`$ take the real part of the completed logarithmic derivative. It is the sum of $`(3/2-\beta)/((3/2-\beta)^2+(v-t)^2)`$ over all true zeros. Subtract the known brackets and their negative partners, supplied by conjugation of the complete family. Every residual is below 0.49. An omitted line zero through 48 would contribute at least one half at a nearby centre. At zero each residual is below 0.30. Thus the total unlisted critical-line mass is below 0.30 in each field and all such ordinates exceed 48 in absolute value. Counts include multiplicity. [Zero verifier](deadline_spectrum.py) · [Integer brackets](evidence/common-deadline-zeros.tsv) · [Platt, Theorem 7.1](https://arxiv.org/pdf/1305.3087#page=14)

### Signed common averages

Centre one probability on each listed $`x`$. Convolve 1,150 uniform boxes of half-width $`1/3000`$ and two of half-width $`1/8`$. Its support radius is $`R=1150/3000+2/8`$. On a retained positive ordinate $`t`$ its transfer is $`U(t)=\mathrm{sinc}(t/3000)^{1150}\mathrm{sinc}(t/8)^2`$, using $`\mathrm{sinc}(v)=\sin(v)/v`$. The positive and negative pair contributes

```math
2G(it)U(t)\cos(tx)\quad\hbox{and}\quad
-\frac{4G(it)\sin(at/2)}a\,U(t)\sin(tx)
```

to the level and move. For unlisted line ordinates $`|t|\gt48`$, the rational formula gives bounds $`18/(1+t^2)`$ for $`G(it)`$ and $`50/(1+t^2)`$ for $`J(it)`$. The two wider boxes give a further factor $`1/(48/8)^2`$. Multiply by the omitted mass cap. Zeros above height 8,000 anywhere in the strip cost at most $`10^6 e^{1040+1/8}(\cosh(1/6000)/(8000/3000))^{1150}`$ for the full allowed centre range $`x\le2080`$. At the trivial zeros use the total variation and the support bound. Their total is at most $`48 e^{-3(x-R-2)/2}/(1-e^{-2(x-R-2)})`$. At our smallest allowed $`x-R`$ the latter and the higher-zero estimate together are below 0.002. Counting some high line zeros twice in these upper bounds preserves the inequalities. [Enclosures in both coordinates](finite_common.py)

For each of four move orthants choose a row whose two signed average moves exceed 4. For four coordinates choose rows for all sixteen orthants, with both signed centred levels above 2 and signed moves above 4. All margins include the arithmetic and spectral uncertainties above. One evaluation uses direct sine, cosine and sinc, and another uses split phases, half-angles and cosine differences. Each row $`(m,s,t)`$ gives the centre $`m/100`$, the orthant as a binary integer in lexicographic sign order with negative first, and the log-height floor. Floor zero denotes height $`10^9`$. The two-coordinate order is $`(D_A,D_B)`$ and the four-coordinate order is $`(L_A-1,D_A,L_B-1,D_B)`$. [860 move rows](evidence/common-2.tsv) · [695 four-coordinate rows](evidence/common-4.tsv)

Each row covers the open interval $`(x+R-\log(C/2),x-R)`$ of possible values $`\log(2rH)`$. Its entire support starts strictly after $`H`$ and finishes before $`CH`$. For each orthant the open intervals overlap strictly, through $`\log(2r)+2000`$. With four coordinates they start at $`H=1`$, and for the two moves they start at $`H=e^{15}`$. The remaining two-move range is covered by eight direct prime-ideal counts. The first four starts exceed 2,048 and finish before $`2^{18}`$. The next four exceed $`e^{15}`$ and finish before $`2048\cdot2^{18}`$. A congruence-class sieve handles good degree-one primes modulo 16 and the residue-degree calculation adds all higher-degree ideals. Interval logarithms enclose both scores on every selected doubling, through norm 166,322,920. Each of the four sign pairs exceeds magnitude 4 in both groups. [Direct count](early_common.py) · [Endpoint inequalities](JointOscillation.lean)

## A common test at every height

**Theorem.** For every $`H\ge1`$, with both completions before $`e^{14004}H^{17001/16999}`$, the move hull contains $`[-3.5,3.5]^2`$ in its interior. The level and move hull contains the product of the two closed ellipses centred at $`(1,0)`$ with semiaxes $`1.17`$ and $`2.45`$. All statements concern actual prime-only scores.

For $`h=\log H\ge100`$ select the following two actual ordinates from the bracket table. All endpoints in the last two columns have denominator $`10^{14}`$.

| Field | Conductor | Conrey label | Lower numerator | Upper numerator |
|---|---:|---:|---:|---:|
| A | 11632 | 1453 | 236949851779030 | 236949851779034 |
| B | 9136 | 5709 | 231411042788551 | 231411042788555 |

Call the ordinates $`g_A,g_B`$. For $`M=8`$ let $`w_k=1-|k|/9`$ when $`|k|\le8`$. The Fejér polynomial $`F(\theta)=|\sum_{j=0}^8e^{ij\theta}|^2/9`$ is nonnegative with coefficients $`w_k`$. Use the product $`F(g_Ax+\theta_A)F(g_Bx+\theta_B)`$. Its 289 terms have frequencies $`\nu=kg_A+lg_B`$ and weights $`w_kw_l`$, with arbitrary phases determined by two independently chosen angles. Their absolute weights sum to 81.

Use two uniform boxes of half-width 3,500 convolved with the compact positive sphere probability of the waiting-bound proof, now with $`c=8510`$. Its Laplace transform is

```math
\Psi_L(z)=\frac{c}{\sinh(cL)}
  \frac{\sinh(L\sqrt{c^2+z^2})}{\sqrt{c^2+z^2}} .
```

The law is that of $`Lv_3`$ on the unit sphere with surface measure tilted by $`e^{cLv_1}`$ and normalized. Rotation proves the formula on the real axis and analytic continuation extends it to the plane. Its modulus on the imaginary axis is at most one. Translate the convolution to centre $`X`$ where

```math
L=(h+14084)/16999,\quad R=7000+L,\quad X=h+7002+L.
```

The support has $`X-R=h+2`$ and $`X+R\lt14004+(17001/16999)h`$. Thus every possible common start exceeds $`H`$ and its doubling completes before the stated deadline. [Sphere kernel and the earlier single-field test](https://github.com/leecrossley/theta-host-split-inert-race/blob/theta-race-waiting-bound-v1.1.0/results/theta-race-waiting-bound/PROOF.md#a-smaller-explicit-exponent-at-every-height) · [Lean support deductions](JointOscillation.lean)

Take a real linear functional with coefficients $`(u_i,v_i)`$ on $`(L_i-1,D_i)`$ for each host. On its own selected zero the complex coefficient is $`u_iG(ig_i)+v_iJ(ig_i)`$, its real and imaginary parts in quadrature. The cancelling harmonics and the conjugate zero contribute twice $`w_1`$ times the absolute value of this coefficient when the phase is chosen to align it. Each host has its own angle, so these choices can be made simultaneously. If a coefficient is zero, any angle for that host suffices. All other terms are bounded uniformly in both angles.

In detail set $`V_i=|J(ig_i)|`$ and $`U(v)=\mathrm{sinc}(3500v)^2`$. The mass of the modulated test lies between $`1-E_0`$ and $`1+E_0`$, where $`E_0=\sum_{(k,l)\ne(0,0)}w_kw_l U(kg_A+lg_B)\lt0.01`$. For a retained positive ordinate $`t`$, put

```math
B_i(t)=\sum_{k,l}' w_k w_l U(t-kg_A-lg_B).
```

The prime on the sum omits exactly $`(k,l)=(1,0)`$ at the selected A row in A, or $`(0,1)`$ at the B row in B. Negative ordinates give conjugates. For the level leakage use $`\sum_t2G(it)B_i(t)`$ and for the move use $`\sum_t2|J(it)|B_i(t)`$. All evaluated nonselected frequency differences exclude zero. Missing line terms add $`18m`$ and $`50m`$ respectively, where

```math
m=\frac{0.30\cdot81}{3500^2(48-8(g_A+g_B))^2}.
```

An off-line zero has absolute ordinate above 8590. Since $`8(g_A+g_B)\lt38`$, every shifted ordinate has magnitude at least 8552. For $`z=u+iv`$ at that distance with $`|u|\le1/2`$, a square root $`P+iQ`$ of $`c^2+z^2`$ has $`Q^2\ge v^2-c^2-1/4\gt700^2`$ and $`|P|\le|v|/(2\sqrt{v^2-c^2-1/4})\lt7`$. Hence $`|\Psi_L(z)|\lt37e^{-8503L}`$. Including the two boxes and the translation gives the exponent $`(X+7000)/2-8503L\le-41`$. The off-line cost per coordinate is below $`81\cdot37\cdot10^6e^{-41}\lt0.0001`$. The trivial zeros have real shifted abscissae $`-3/2-2j`$. The pointwise bound 81 on the polynomial product and $`X-R-2=h\ge100`$ give a trivial-zero contribution below 0.0001. The calculation assigns 0.001 to these two sources together. These estimates apply for every centre, regardless of the positions of the higher zeros. [Mixed-frequency checks](all_height_common.py)

Divide the main amplitudes by $`1+E_0`$, obtaining $`b_i=2w_1G(ig_i)/(1+E_0)`$ and $`d_i=2w_1 V_i/(1+E_0)`$. Divide the two spectral error bounds by $`1-E_0`$, then add the corresponding prime-square and higher-power costs. Call the resulting positive errors $`e_i,f_i`$. The normalized average of the chosen functional is at least

```math
\sum_i\left(\sqrt{(b_iu_i)^2+(d_iv_i)^2}
                  -e_i|u_i|-f_i|v_i|\right).
```

By Cauchy–Schwarz the $`i`$-th summand is at least $`p_i\sqrt{(b_iu_i)^2+(d_iv_i)^2}`$ where $`p_i=1-\sqrt{(e_i/b_i)^2+(f_i/d_i)^2}`$. Two implementations check $`p_ib_i\gt1.17`$, $`p_id_i\gt2.45`$ and $`d_i-f_i\gt3.5`$ for both fields. The last inequality applies to a move-only functional directly. For $`H\lt e^{100}`$ use the stronger boxes of the finite theorem. Their completions are before $`e^{100}2^{80}\lt e^{14004}`$. The large-height tests and the finite calculations together prove the theorem. [Sinc and cosine evaluations](all_height_common.py)

## Why these are actual prime observations

At any fixed starting height all witnesses for the relevant comparison live in a compact interval of eligible starts, strictly inside its completion limit. In the finite theorem there are finitely many support intervals or direct witnesses. Enclose their union in its compact interval hull. In the infinite theorem the kernel itself has one compact support interval. The four-coordinate score curve is continuous throughout either interval, so its image and its convex hull are compact.

For a real projection of the finite box choose the orthant agreeing with the signs of its coefficients. The certified positive mean, or the early exact observation, gives a strict support-function margin in that direction. In the all-height theorem choose the two modulating phases as above. Separation puts the specified boxes or ellipse products in the actual hull, in its interior since every nonzero direction has a strict margin. The connected-set refinement of Carathéodory's theorem, applied along the same interval of eligible starts, gives at most two common starts for the two moves and at most four for all coordinates. This reproduces any target in the claimed sets with one common vector of nonnegative weights. [Original connected-curve proof](#common-starts-and-the-number-of-readings) · [Bárány and Karasev](https://www.renyi.hu/~barany/cikkek/131.pdf)

A single linear functional takes a value at least as large as any of its positive averages. With direction $`(1,-1)`$ or $`(-1,1)`$ the strict support inequalities give a single doubling where one normalized move exceeds the other by more than 8 in the finite theorem or 7 at every height. Multiplication by $`\log2`$ restores the original units. More generally the two-sided margins are $`4(|s|+|t|)`$ and $`3.5(|s|+|t|)`$ for any fixed nonzero linear combination $`sD_A+tD_B`$.

## Larger contrasts on power intervals

**Theorem.** Fix $`(s,t)\ne(0,0)`$ and $`\delta\gt0`$. For all sufficiently large $`H`$ the projection $`sD_A+tD_B`$ takes a value above $`5.85(|s|+|t|)`$ and another below its negative, at starts greater than $`H`$ whose doublings finish before $`H^{1+\delta}`$. The height can depend on the direction and on $`\delta`$.

For each field the good-square correction to the level tends to 1 by its unconditional prime ideal theorem and compactly supported partial summation. Explicitly the unsigned sum against $`W_2((Nv/Z)^2)`$ divided by $`BZ`$ tends to $`\int_1^2W_2(v^2)\,dv/B=1`$. Higher powers and finite omissions contribute $`o(1)`$. The move's arithmetic correction consequently tends to zero. Removing the arithmetic error from the mixed-frequency computation gives a move amplitude after spectral errors above 5.85 in each field. This check allows arbitrary modulating phases and both signs.

The Mellin transform of the fixed real projection on the half-plane $`\Re z\gt0`$ is meromorphic. Its possible poles there are shifted odd-factor zeros, with coefficient $`J(z)(s m_A+t m_B)`$ times a nonzero centering exponential, where $`m_i`$ counts multiplicity in the $`i`$-th family. The prime-square logarithmic derivative is holomorphic at $`1+2z`$ and the higher powers are holomorphic there as well. The tent's possible zeros have real part $`\pm1/2`$ and the move multiplier's zeros have real part zero. Every uncancelled off-line residue on the right is therefore nonzero. There are no real poles there by finite verification.

If there is no such pole, the off-line terms in the projected odd prime-power formula cancel with their multiplicities. The functional equations reflect the multiplicities within each family and hence also cancel the left-hand partners. The remaining nontrivial terms give an absolutely convergent trigonometric series on the log axis by the total mass bound. For instance, convolve first with small centred probability boxes, apply the explicit formula, and shrink the boxes at fixed positive log height; the line series and the trivial zeros are dominated there. The latter vanish in the limit at large log height. With the arithmetic correction restored, the actual projection is this uniformly almost periodic series plus $`o(1)`$. Apply the same mixed-frequency test to the series, omitting the sphere factor (its earlier imaginary-axis majorant was one). For each sign the supremum exceeds $`5.85(|s|+|t|)`$ strictly. Almost periodicity gives such values with a fixed margin in each sufficiently long log interval. The actual prime scores inherit both signs.

If there is a pole on the right, let $`\Theta\gt0`$ be the supremum of the real parts of these poles. For any $`0\lt\eta\lt\Theta`$ add and subtract $`\mathbf1_{y\ge1}y^\eta`$ from the projected observable separately. Its Mellin transform converges absolutely far to the right, since prime sums have polynomial growth and vanish sufficiently near zero. In $`\Re z\gt\Theta-c_0`$ with $`0\lt c_0\lt\Theta-\eta`$, the two transforms have precisely the same nonreal poles, with no finite accumulation, conjugate in pairs and with ordinary pole expansions. These are the hypotheses of the power-interval theorem of Kaczorowski and Pintz. Each modified function changes sign on every sufficiently late $`[Y^{1-\varepsilon},Y]`$ for each fixed $`0\lt\varepsilon\lt1`$. The original projection exceeds $`y^\eta`$ and falls below $`-y^\eta`$ there. Put $`Y=H^{1+\delta/2}`$ and $`\varepsilon=(\delta/4)/(1+\delta/2)`$. All starts exceed $`H`$ and $`2Y\lt H^{1+\delta}`$ eventually. The two alternatives prove the claim without a global zero-location hypothesis. [Kaczorowski and Pintz, theorem (1.7)](https://personal.math.ubc.ca/~gerg/teaching/592-Fall2018/papers/1987.Kaczorowski_2.pdf) · [Single-field application](https://github.com/leecrossley/theta-host-split-inert-race/blob/theta-race-waiting-bound-v1.1.0/results/theta-race-waiting-bound/PROOF.md#every-positive-excess-exponent)
