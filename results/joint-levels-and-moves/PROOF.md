# Joint levels and moves of the two theta races

## The two observables

Let $`A=(547,727)`$ and $`B=(571,599)`$. In each field of the [statement](README.md) put $`\epsilon_v=1`$ for split and $`\epsilon_v=-1`$ for inert. Omit the primes over 2 and the two conductors. Additional fixed finite omissions have no effect on the tail. For the arithmetic trace interpretation omit the exceptional set of the chosen admissible datum. The first-step parameter recursion subtracts the previous twisted trace and divides by the trace of the added segment. The fibre product of these normalized jumps equals $`\epsilon_v`$. Both hosts have explicit admissible data. [First trace extraction](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-prime-trace/ThetaPrimeTrace.lean#L68) · [547/727 host](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-short-shifts/RaceTrace.lean#L83) · [571/599 host](https://github.com/leecrossley/theta-host-split-inert-race/blob/8fca7cad37185e3d7951f0dffff5e2442d1dfa24/HostAdmissible.lean)

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

The square correction restores the unsigned prime squares. The signed powers of order at least three and the higher terms of the correction are absolutely convergent in this half-plane. The omitted finite factors are regular there. The real-axis certificate below excludes real zeros from $`[1/2,1]`$ for all sixteen factors. Euler products handle $`s\gt 1`$. Thus the real positive axis of $`M_i(1/2+z)`$ is regular, and the residue at $`z=0`$ is $`1/2`$. At $`z=i\gamma\ne0`$ the prime ideal theorem keeps the base zeta nonzero on $`\Re s=1`$. The residue equals the nonnegative total zero multiplicity $`m_i(\gamma)`$ of the eight factors at $`1/2+i\gamma`$. [Euler products and continuation](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-short-shifts/PROOF.md#euler-factors-and-the-kernel) · [Real-axis cover](real_axis.py)

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

Fix both widths at two and retain the prime and score conventions of the two explicit theta hosts. In field A omit primes over 2, 547, 727, 343057911569 and 343070636977. In field B omit primes over 2, 571, 599, 413718278401 and 603339156001. Write $`r=\sqrt2`$, $`a=\log2`$, $`B_0=2r(r-1)^2`$, and now write $`x=\log(2r y)`$ for the log centre of a doubling with start $`y`$. Then the midpoint score is $`L_i=S_i(ry)`$ and the normalized move is $`D_i=(S_i(2y)-S_i(y))/a`$. A completion occurs at $`2y=e^x/r`$. [547/727 trace](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-short-shifts/RaceTrace.lean#L83) · [571/599 trace](https://github.com/leecrossley/theta-host-split-inert-race/blob/8fca7cad37185e3d7951f0dffff5e2442d1dfa24/HostAdmissible.lean)

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

## One fixed window for all directions

**Theorem.** Fix the two widths $`b_A,b_B\gt1`$. There is a finite $`C`$ such that after every real $`H\ge1`$ the product of the two closed level–move ellipses in the original table is attained by positive averages of at most four **common** starts. The weights for the two fields are identical. Every start satisfies $`y\gt H`$ and every doubling finishes before $`CH`$. At width two for both fields the move hull also contains $`[-5.85,5.85]^2`$ in its interior. Each point of this square takes at most two common starts. On single doublings either field can outgain the other by more than $`11.7\log2`$, and either midpoint level can exceed the other by more than 17.9106. One constant works for all targets and all directions at the specified widths. Its numerical value is open.

There is also a rate independent of any unknown zero positions. For every fixed $`\varepsilon\gt0`$, all these bounded targets fit before $`H(\log H)^{8+\varepsilon}`$ for sufficiently large $`H`$. The proof gives no numerical onset. The earlier finite certificates give numerical small-height factors.

### The vector of residues

Use the four centred real observations $`V=(L_A-1,D_A,L_B-1,D_B)`$ at log starting height $`x=\log y`$. Both completions are at $`2e^x`$. The same argument applies to the projection onto the two moves. Write $`v_{\alpha,g}\in\mathbf C^n`$ for the combined coefficient at a positive ordinate $`g`$ and real exponent $`\alpha=\beta-1/2`$, with $`n=4`$ or $`2`$. For each host this coefficient is the multiplicity times the evaluated level and move kernels, including the translations in the residue formula above. Add contributions at the same complex zero before forming the vector. Allow one of the two hosts to have multiplicity zero. The negative ordinate contributes the conjugate vector.

For fixed widths the kernels at complex exponents have size $`O_b((1+g^2)^{-1})`$ throughout the strip. This follows from the two linear denominator factors, bounded hyperbolic numerators, removable points, and the bounded translation factors at fixed widths. The total zero estimate gives $`\sum_{\alpha,g}\|v_{\alpha,g}\|\lt\infty`$. After centering, the signed vector is the sum of $`2\Re(v_{\alpha,g}e^{(\alpha+ig)x})`$ plus a term tending to zero as $`x\to\infty`$. Indeed the unsigned prime-square level tends to its bias by the unconditional prime ideal theorem in either fixed field, and its dyadic change tends to zero. Higher powers, trivial zeros and excluded primes have vanishing errors. The nontrivial odd zeros give the series by contour inversion with two boxes. Shrinking the boxes on any compact interval of positive log heights recovers the continuous prime scores by absolute convergence there. A central odd zero is excluded by the real-axis certificate. Consequently the critical-line contribution is a bounded uniformly convergent Fourier series. All negative-real-part terms together tend to zero by dominated convergence. [Explicit formula and arithmetic corrections](#prime-power-and-zero-errors) · [Real-axis certificate](real_axis.py)

Let $`W\subseteq\mathbf R^n`$ be the real span of the real and imaginary parts of the vectors $`v_{\alpha,g}`$ for $`\alpha\gt0`$. Write $`P_W u`$ for orthogonal projection onto $`W`$. A unit real functional $`u`$ satisfies $`|u\cdot v_{\alpha,g}|\le\|P_W u\|\|v_{\alpha,g}\|`$ on this part of the spectrum. Thus the same small parameter controls **every** growing coefficient of a direction near $`W^\perp`$. This uniformity will be needed when the direction changes with $`H`$.

If $`W\ne0`$, choose finitely many actual pairs $`(\alpha_j,g_j)`$ whose residue vectors have real and imaginary parts spanning $`W`$. Set $`d_*=\min_j\alpha_j\gt0`$ and fix $`0\lt d_0\lt c\lt d_*`$. There is a constant $`b\gt0`$ such that

```math
\max_j |u\cdot v_{\alpha_j,g_j}|\ \ge\ b\|P_Wu\|.
```

The maximum on the left is a continuous norm on the unit sphere of $`W`$. Compactness gives a positive minimum. Montgomery’s Dirichlet zero-density theorem implies that the positive ordinates of all zeros with $`\alpha\ge d_0`$ have summable reciprocals, including multiplicity, for this finite family of sixteen characters. The exponent $`3(1-\sigma)/(2-\sigma)`$ is below one at $`\sigma=1/2+d_0`$. Summing the dyadic count divided by the lower endpoint proves convergence. [Montgomery’s theorem](https://www.mathnet.ru/eng/mat564) · [Single-field proof and the density argument](https://github.com/leecrossley/theta-host-split-inert-race/blob/theta-race-waiting-bound-v1.2.0/results/theta-race-waiting-bound/PROOF.md#a-fixed-factor-without-a-global-zero-hypothesis)

### Keeping the low-frequency tests

The nine-frequency product above has a finite set of harmonic shifts of absolute value below 6. At the dyadic widths also take every one of the 289 shifts of the two-carrier Fejér product, all of absolute value below 38. Include their negatives and zero. Call this finite set $`S`$. Every positive growing ordinate exceeds 8590 by Platt’s finite theorem for these conductors.

For $`W\ne0`$ enlarge $`S`$ by the finitely many selected $`\pm g_j`$ and call the enlarged set $`S'`$. For each $`\alpha\ge d_0`$, each positive ordinate $`\gamma`$ and every $`s\in S'`$ with $`\gamma+s\ne0`$, convolve a positive symmetric probability of radius $`\pi/|\gamma+s|`$ killing the exponent $`\alpha+i(\gamma+s)`$ and its conjugate. Take the reflected pair of skew two-point probabilities from the single-field proof. Above twice the largest shift each radius is at most $`2\pi/\gamma`$. The sum of the radii is finite, including the finitely many smaller ordinates. Independent draws therefore define one compactly supported probability with every required zero moment. If $`W=0`$, use an atom at zero. Convolve also with two fixed uniform boxes, whose half-width can be as large as required below. Call the resulting probability $`P`$ and its full radius bound $`R`$. It is the same probability for all directions and centres. At every real $`t`$, $`\widehat P(t)=\int e^{tx}\,dP(x)\gt0`$. Its real-frequency responses have modulus at most one.

Translate to $`X=h+1+R`$ for $`h=\log H`$. Multiply pointwise by either low-frequency trigonometric product, choosing phases for the given unit direction. The support is strictly beyond $`h`$. All growing terms with $`\alpha\ge d_0`$ are killed, since $`|\gamma|\gt8590`$ and $`|s|\le38`$ for these tests. The remaining growing terms cost at most $`K e^{d_0 h}\|P_Wu\|`$, where $`K`$ is fixed: the support has a fixed width, the weights are bounded, and all residue norms are summable. If $`W=0`$ this error is zero. The arithmetic and left-half errors are $`o(1)`$ uniformly on the test interval and in unit directions.

For the nine-frequency product the mass tends uniformly to one as the two box widths increase. Every nonconstant harmonic has nonzero frequency. On the critical line only the matched singleton harmonics coincide with zeros, by the 157,384 exclusions above. The sinc-square multiplier tends to zero on all the other shifted frequencies. Their total absolute coefficients are summable, uniformly in the direction and phases by the triangle inequality. Hence the weights can be kept nonnegative while bringing the support lower bound arbitrarily close to the sum of the nine matched absolute coefficients. This remains true with the zero-cancelling probability because its other imaginary-frequency responses have modulus at most one, and its value at zero is one. The strict margins between the unrounded and printed axes allow **one finite box width for all unit directions**. For example, if both semiaxes of each ellipse strictly exceed the printed axes, the support difference is positive on the unit sphere by continuity. Include in the margin the strict unrounded 17.9106 assertion in the two level-difference directions as well.

At width two the previous 289-term calculation similarly gives a strict $`5.85(|u_A|+|u_B|)`$ margin on the critical-line move series. Use two boxes of half-width 3500 for this bound. Additional boxes for the nine-frequency test only reduce the absolute leakage and preserve each matched harmonic and the uniform mass bound. Thus one common measure handles both conclusions. When $`\|P_Wu\|\le e^{-ch}`$, the growing-term error is $`O(e^{-(c-d_0)h})`$. This proves the required strict support bound uniformly in these directions for all sufficiently large $`h`$. It proves every direction when $`W=0`$. [Nine true ordinates](evidence/brackets.tsv) · [Two-harmonic margin calculations](all_height_common.py)

### Resolving collisions at a common ordinate

Suppose now $`\|P_Wu\|\gt e^{-ch}`$. An unmodulated average at any of the finitely many translated centres $`X+\ell`$ below kills every term with $`\alpha\ge d_0`$. Its absolute value is $`O(1+e^{d_0 h}\|P_Wu\|)`$. This includes the bounded critical-line series. For a selected ordinate $`g`$ form the *complex* average of $`u\cdot V(x)`$ times $`e^{-igx}`$ under that translated probability. Let $`d_1,\ldots,d_m`$ be the distinct exponents $`\alpha\ge d_0`$ actually contributing at this ordinate. Then for $`\ell=0,\ldots,m-1`$ these averages are

```math
\sum_{j=1}^m e^{d_j\ell}\,e^{d_jX}
       \widehat P(d_j)(u\cdot v_{d_j,g})
 + O(1+e^{d_0 h}\|P_Wu\|).
```

The negative ordinate moves away from zero in this modulation and is killed along with the other unmatched positive and negative ordinates. For complex averages the displayed leading coefficient comes once from the positive ordinate. The $`m`$ real exponentials $`e^{d_j}`$ are distinct. The matrix with rows $`\ell`$ is an invertible Vandermonde matrix. Its inverse has a fixed finite norm, so at least one of the complex sums before the error has modulus bounded below by a positive constant times the largest $`e^{d_jX}|u\cdot v_{d_j,g}|`$. This remains true if modes at that ordinate cancel at one translation. Use all the finitely many selected ordinates and the spanning inequality. Some centre and ordinate therefore supply a lower bound $`b' e^{d_*h}\|P_Wu\|`$ before errors, for a fixed $`b'\gt0`$. [Vandermonde lemma](JointOscillation.lean) · [Worked example of the collision inversion](fixed_windows.py)

At that centre use the positive weight $`1+\cos(gx+\theta)`$ with phase chosen to extract the full modulus of the complex average. Adding back the unmodulated average costs at most another $`O(1+e^{d_0h}\|P_Wu\|)`$. The mass is at most two. For all sufficiently large $`h`$, uniformly in this part of the unit sphere, the resulting supremum is at least a constant times $`e^{d_*h}\|P_Wu\|`$. Indeed the ratio of the error to the main term is $`O(e^{-(d_*-c)h}+e^{-(d_*-d_0)h})`$. In particular this supremum exceeds every fixed target support. Opposite directions are included.

All the centres have $`0\le\ell\lt M`$ for one fixed integer $`M`$. Their supports lie in a common compact interval strictly above $`h`$, with the doubling endpoint below $`H\exp(2+2R+M+\log2)`$. For all sufficiently large $`H`$ the support in every direction strictly exceeds that of the promised bodies. Separation puts them in the interior of the convex hull of the prime curve on that interval. The connected-curve argument above supplies two common starts for moves and four for the four-coordinate curve. To cover earlier $`1\le H\le H_0`$, take the entire working interval from one sufficiently large $`H_0`$, which lies above every such start, and increase the fixed factor past its largest numerical completion. This single interval already handles *all* targets.

There is a further consequence when $`W\ne0`$. For any $`0\lt\eta\lt d_*`$ choose $`\eta\lt c\lt d_*`$ and $`d_0\lt c`$. In the first part of the sphere the extra target support $`H^\eta\|P_Wu\|`$ tends uniformly to zero. In the other part it is absorbed by the exponential main term. Thus in every sufficiently late fixed-factor window the hull contains the stated bounded body plus the radius-$`H^\eta`$ Euclidean ball **in $`W`$**. This gives growth precisely in the selected spectral directions as a lower bound. The projection onto $`W^\perp`$ stays bounded everywhere, since it annihilates every growing residue. These statements also cover shared spectra and proportional residue vectors.

### A numerical exponent independent of the zeros

For completeness, a different construction handles the *bounded* targets without the spectral case split. Write $`T=10^6+4h`$, $`L=\log(12000T)`$ and take $`U=L^4`$ for sufficiently large $`h`$. Use the separated positive-interpolation lemma for the compact cosine density from the single-field paper. With $`D=16`$ primitive factors and the quartet symmetry, the number of growing positive representatives in a closed interval of length $`d`$ up to $`T`$ is bounded by

```math
\frac{D}{4\pi}\,dL+O_D(L).
```

The error is uniform in $`d`$ by subtracting the two endpoint zero estimates. Greedy grouping makes each group's adjacent positive ordinates at least $`d`$ apart. Below $`U`$ use separate groups for the finitely many shifts of the bounded trigonometric tests, with the fixed cosine scale 1 and spacing 4. They have total radius $`O_S(\log U)=O_S(\log L)`$. Convolve them before modulation. All positive ordinates after shifting stay above 8552. Above $`U`$ through $`T`$ take unshifted groups at scale $`a=L^{-2}`$ and spacing $`d=s/a`$ with $`s=2+1/L`$. The scaled first frequency is at least $`L^2`$. The cosine-density lemma applies, since its column norm is at most $`1-1/(2L)`$ and its perturbation norm is less than one after inversion. Use these groups *after* multiplying by the trigonometric weight. The radius is at most

```math
\frac{\pi a}{2}\left(\frac{D}{4\pi}\,\frac{s}{a}L+O_D(L)\right)
=\frac{D}{8}sL+O_D(aL)=\frac D4 L+O_D(1).
```

At every one of the finitely many matched critical frequencies the product of the outgoing responses has modulus $`1-O(aL)=1-o(1)`$, by the support-based cosine bound. Choose each singleton phase to compensate its response phase. All critical-line leakage estimates remain valid. The terms to the left of the line tend to zero uniformly after translating past $`h`$, by the same summable majorant for each zero and the positivity and bounded total mass of the tests. Add $`\lceil T/4\rceil`$ uniform boxes of half-width $`4/T`$ before modulation. Their radius is $`O(1)`$. For $`|\gamma|\gt T`$ their bound after any of the shifts $`|s|\lt38`$ includes $`\exp(-\lceil T/4\rceil)`$ once $`4(1-38/T)\gt e`$. After the centre and total radius cost, the unverified growing tail tends to zero uniformly since $`h/2+O(\log T)-T/4\to-\infty`$. The earlier absolute strip mass supplies the finite multiplier. The fixed long boxes and bounded weights allow the same explicit-formula inversion. All error bounds are uniform on the relevant unit sphere.

The full log overhead is twice the total radius plus an absolute centering and completion allowance. It is at most $`(D/2)L+O_S(\log L)`$. For every fixed $`\varepsilon\gt0`$ this gives completion before $`H(\log H)^{8+\varepsilon}`$ for all sufficiently large $`H`$. The same strict support and connected-curve argument applies. The quantitative exponent does not assert a bound for the size of the fixed factor. [Positive interpolation and box estimates](https://github.com/leecrossley/theta-host-split-inert-race/blob/theta-race-waiting-bound-v1.2.0/results/theta-race-waiting-bound/PROOF.md#a-positive-density-for-separated-complex-frequencies) · [Bennett–Martin–O’Bryant–Rechnitzer](https://doi.org/10.1090/mcom/3599)

### How many observations are necessary

For any finite interval of positive $`y`$ only finitely many prime norms meet the moving tents. Between breakpoints each coordinate is a linear combination of $`\sqrt y`$, $`1/\sqrt y`$ and a constant. At a breakpoint it is continuous. The actual move-pair curve is therefore locally Lipschitz, and its image in the plane has area zero, even after taking a countable union of bounded starting intervals. Almost every point of the guaranteed move square needs two starts. For two common starts and one common weight, the four-coordinate image is a locally Lipschitz image of three real parameters in four dimensions, so has four-dimensional measure zero. This follows by covering bounded parameter boxes by cubes of side $`1/N`$, whose images have diameter $`O(1/N)`$. Taking a countable union covers all starts and weights. Thus almost every point in the ellipse product needs at least three common starts. The connected-set bound alone leaves the three-start question open. The next section handles the dyadic ellipse product through an explicit larger box under its stated hypotheses.

## Three common starts for the four measurements

Write $`a=\log2`$. Order the centred observations as $`Z(y)=(L_A(y)-1,L_B(y)-1,D_A(y),D_B(y))`$. The four normalized targets form the closed box

```math
\mathcal R=[-10,10]\times[-8,8]\times[-7/2,7/2]\times[-11/5,11/5].
```

The two level intervals in the original coordinates are $`[-9,11]`$ and $`[-7,9]`$. Both previously certified dyadic ellipses, as a product, lie in this box. The result concerns positive real starting heights and the same weights for the two fields. Each start supplies both midpoint levels and both normalized moves.

### Three readings from two arcs and one point

Take a continuous four-dimensional curve $`W`$, a singleton $`v_0=W(a_0)`$, and two closed parameter intervals with endpoints $`b_1,c_1`$ and $`b_2,c_2`$. Enumerate their images as $`v_1,v_2,v_3,v_4`$ in that order. Assume the five vertices are affinely independent. Let $`\lambda_i`$ be their barycentric affine coordinates, so $`\lambda_i(v_j)=\delta_{ij}`$ and $`\sum_i\lambda_i=1`$ everywhere. Denote the two whole arc images by $`A_1,A_2`$. Define five upper bounds

```math
m_0=\max(0,\sup_{A_1\cup A_2}\lambda_0),\qquad
m_i=\max(0,\sup_{A_2}\lambda_i)\quad(i=1,2),\qquad
m_i=\max(0,\sup_{A_1}\lambda_i)\quad(i=3,4).
```

If $`\lambda_i(p)\gt m_i`$ for every $`i`$, there are three points, one on each designated piece, whose common positive average is $`p`$. To prove this, take $`\alpha`$ in the standard four-simplex, put $`s_1=\alpha_1+\alpha_2`$ and $`s_2=\alpha_3+\alpha_4`$ and form

```math
\Phi(\alpha)=\alpha_0v_0+
s_1W\!\left(\frac{\alpha_1b_1+\alpha_2c_1}{s_1}\right)+
s_2W\!\left(\frac{\alpha_3b_2+\alpha_4c_2}{s_2}\right).
```

At zero group mass take its summand to be zero. Continuity there follows from boundedness on each compact arc. On the face $`\alpha_0=0`$, both contributing arcs have zeroth coordinate at most $`m_0`$. On $`\alpha_i=0`$ with $`i=1,2`$, the first arc is pinned to the other endpoint or has zero mass. The singleton and the pinned endpoint have $`i`$-th barycentric coordinate zero. The only remaining contribution comes from the second arc and is at most $`m_i`$ in that coordinate. The identical assertion for $`i=3,4`$ follows by exchanging the arcs. Thus $`\lambda_i(\Phi)\le m_i`$ everywhere on face $`i`$.

Put $`d_i(\alpha)=\lambda_i(p)-\lambda_i(\Phi(\alpha))`$. Their sum vanishes. Normalize the five nonnegative entries $`\max(0,\alpha_i+d_i(\alpha))`$ by their sum. This sum is at least one, so the normalized map is a continuous self-map of the simplex. Brouwer's theorem gives a fixed point. At a zero coordinate the corresponding unnormalized entry would be strictly positive by the face bound. All coordinates of the fixed point are therefore positive. Clipping changes no entry there, the normalizing sum is one, and $`d_i=0`$ for all $`i`$. The affine coordinate map is injective, hence $`\Phi(\alpha)=p`$. The Lean file proves the fixed-point algebra and the three-term face inequality. The use of Brouwer's theorem and continuity of the grouped map are written steps. [Lean deductions](JointOscillation.lean) · [Intermediate value theorem for simplices](https://www.sciencedirect.com/science/article/pii/S0166864119304420)

The same five affine coordinates may be used after a perturbation of the parameterized pieces. If the perturbation changes each of its four entries by at most $`\varepsilon`$, every face coordinate of the grouped average changes by at most $`\varepsilon\max_i\|\nabla\lambda_i\|_1`$. This follows by applying the affine function to a convex average with weights summing to one. Strict face margins then give an exact solution on the perturbed arcs. Compactness alone would only give approximations without this last surjectivity argument.

### Applying the criterion to prime ideals

Use the centred coordinate $`x=\log(2\sqrt2\,y)`$. Put $`B=2\sqrt2(\sqrt2-1)^2`$ and $`g(t)=4\sinh((a-|t|)_+/2)/B`$. The level kernel at scale $`z=e^x`$ applied to a norm $`n`$ is $`g(x-\log n)/\sqrt n`$. For each host $`(p,q)`$ omit the rational primes $`2,p,q`$. The two auxiliary type primes of each displayed admissible host exceed three million and exceed the norm support of the finite trace certificate below. For every remaining rational prime $`\ell`$ put $`u=(\ell/p)`$, $`v=(\ell/q)`$. Its residue degree $`f`$ is the least element of $`\{1,2,4\}`$ satisfying $`\ell^f\equiv\pm1\bmod16`$ and $`(uv)^f=1`$. The $`8/f`$ base primes have common norm $`n=\ell^f`$ and sign $`u^f`$. Their total signed log weight is $`-8u^f\log\ell`$. The primary evaluation computes the characters by listing nonzero squares. The second evaluation factors the real quartic, the real quadratic defining the compositum and the imaginary quadratic modulo $`\ell`$. It uses their degrees to determine the splitting sign. [Host coefficients and exceptional primes](deadline_arithmetic.py) · [Trace product at a good prime](https://github.com/leecrossley/theta-host-split-inert-race/blob/2449f510a992b76307299d527c4e7e1826be2c8e/results/theta-race-short-shifts/RaceTrace.lean#L83)

Write $`P_i(x)`$ for the prime-only midpoint level obtained from these weights. The actual curve in this coordinate has first two entries $`P_i(x)-1`$ and last two entries $`(P_i(x+a/2)-P_i(x-a/2))/a`$. The finite test covers this exact curve. Its certificate uses 320 unit cells in the two level coordinates. In every cell, both move coordinates range independently across their entire target intervals. For each fixed chart and face the minimum of an affine target coordinate is attained at a corner of the four-dimensional cell. Every one of the 320 cells has a chart on which all five of its minima exceed the corresponding whole-face upper bound by more than $`1/250`$. Thus a chart can be chosen from the two target levels alone, before either requested move is specified. [Chart endpoints in thousandths of a centred log height](evidence/three-start-charts.tsv) · [Full finite verifier](three_starts.py)

### One timing table for both calculations

The 23 rows of the table specify three disjoint parameter sets, a singleton and two closed intervals. Use thousandths of $`x=\log(2\sqrt2\,y)`$. The full table runs from 8.092 to 13.808. Before rounding, in every row the three pieces are separated from their neighbours by at least 0.720 and the largest total span is 4.004. In the physical calculation replace the singleton by the nearest integer starting height. Both verifiers evaluate the original tents at this exact integer and bound its log displacement by 0.0005. The other two parameter intervals retain their listed endpoints. Consequently chronological adjacent starts obey $`y_{j+1}\gt(41/20)y_j\gt 2y_j`$ and the last doubling finishes before $`110y_1`$. All three measurements use disjoint doublings. In the physical table each start exceeds 1,155 and every completion is less than 701,816. The rounded singleton stays strictly inside the common interval. For a working translated table all targets together fit inside a factor of 608 measured from the left edge of that table. The wait to the next working table from an arbitrary large height is a separate constant.

Use the very same 23 timing rows for the finite and limiting curves, rounding the singleton in the physical calculation as just described. Each curve determines its own barycentric coordinate functions from its five vertex values. For a fixed target cell take the smaller of the two verified clearances on each row. At least one row has common clearance greater than $`1/250`$. The printed 20-by-16 lookup gives one such row in each cell. A target on a boundary can use either adjacent closed cell. Only its two target levels enter the choice. The actual starts within the indicated arcs and their weights may vary with all four targets. [Timing rows](evidence/three-start-charts.tsv) · [Lookup and checks](evidence/three-starts.txt)

Every affine row in either calculation has gradient of $`\ell^1`$-norm less than $`3/20`$ and its value on every point of either of its own chart arcs is less than $`3/2`$. At the five vertices these values are zero or one as prescribed. For the singleton take coordinate 0, for the first arc coordinate 1 and for the second arc coordinate 3. If the mass of that designated piece is $`w`$, the other two pieces contribute a weighted average bounded above by $`m_i`$, while that piece contributes at most $`3/2`$ before weighting. The whole coordinate is at most $`m_i+3w/2`$. A perturbed group average differs in each affine coordinate by at most $`3/20000`$ if each underlying measurement error is at most $`1/1000`$. Therefore a target with clearance $`1/250`$ has

```math
w>\frac{1/250-3/20000}{3/2}=\frac{77}{30000}>\frac1{400}
```

on each of its three pieces. The perturbation bound applies even when the realization on the perturbed trace uses a different point in each arc. This proves the weight and disjoint-doubling assertions on every good translated table as well as directly on the original physical table.

### The limiting curves and late physical starts

For each host use the eight primitive odd Dirichlet factors defined earlier. Suppose all their nontrivial zeros lie on the critical line. Alpöge’s announced Dirichlet GRH theorem supplies this input. The directly enumerated finite prime theorem above has no zero-location premise. The square correction has limiting level one by the prime ideal theorem for the fixed base field. Higher prime powers, the trivial-zero remainder and fixed finite omissions vanish on every sufficiently high tail. The centred level is therefore $`V_i(x)+o(1)`$ where $`V_i`$ is its critical-line Fourier series. Its coefficients are absolutely summable by the dyadic kernel decay and the zero count. The functions $`V_i`$ are real, even and uniformly almost periodic. The limiting four-vector has first two entries $`V_i(x)`$ and last two entries $`(V_i(x+a/2)-V_i(x-a/2))/a`$. These last entries are odd. The four functions are jointly uniformly almost periodic as uniform limits of finite trigonometric sums on the same real axis. This uses no independence among zero ordinates. [Explicit formula on both hosts](#the-vector-of-residues) [Announcement](https://alpo.ge/remark.pdf) [Upstream primitive theorem](https://github.com/leanprover-community/riemann-hypothesis/blob/371a1a26d5e17f604428e3a8ce6ce341cf517434/src/Main.lean#L12) [Independent statement and kernel checks](https://math.mit.edu/~drew/grh/)

For $`x\ge0`$ the same limiting level can be evaluated without a zero table. Use $`V_i(x)=\Gamma_i(x)-\sum_{n=\ell^k}\log(\ell)u_k(\ell)g(x-\log n)/\sqrt n`$, with $`u_k(\ell)=\sum_\chi\chi(\ell)^k`$. The reflected prime terms disappear by the support of $`g`$. For odd $`\ell`$ the multiplier is $`4((\ell/p)^k+(\ell/q)^k)`$ when $`\ell^k\equiv\pm1\bmod16`$, and zero otherwise. At two only the two untwisted quadratic characters remain. In $`\Gamma_i`$ use the conductors $`p,8p,16p,16p,q,8q,16q,16q`$. If $`r=\sqrt2`$, $`c=2(r-1)^2`$, $`A_i=\sum_\chi\log(q_\chi/\pi)-8(\gamma_E+a)`$ and $`s=e^{-|t|}`$, the compact-test formula is

```math
T(t)=\frac{(1+s)\log(1+s)+(1-s)\log(1-s)}{2\sqrt s},\quad T(0)=a,
\qquad
\Gamma_i(t)=A_i g(t)-\frac8c\{r[T(t+a)+T(t-a)]-3T(t)\}.
```

The local Euler calculation uses primitive factors even at the omitted physical primes. The physical weight uses unramified base ideals and the stated omissions. The two finite formulas enter different steps of the proof. [Archimedean kernel and its derivation](https://github.com/leecrossley/theta-host-split-inert-race/blob/2449f510a992b76307299d527c4e7e1826be2c8e/results/theta-race-curvature-euler-factors/PROOF.md#remove-the-archimedean-term) · [First interval evaluator](three_profiles.py) · [Conrey-character evaluation](three_recheck.py)

Let $`I`$ enclose all limiting chart parameters. There is one uniform face margin valid for all target cells. Choose a set of common almost periods that makes the translate of the entire limiting vector within half this allowable error. Such periods have bounded gaps. For sufficiently late translates, the actual four-vector differs from the limiting one on the whole shifted interval by less than the other half. The face argument with the same affine coordinates gives every target **exactly** at actual prime starts in each such translated interval. There is one finite multiplicative waiting factor beyond every $`H\ge1`$ since the log intervals have fixed length, their good translates have bounded gaps, and one sufficiently late interval also serves all bounded earlier starts. The proof gives no numerical bound for that waiting factor.

The late physical singleton may also be taken to be the nearest integer to its nominal translated starting height. Its log displacement tends uniformly to zero over the finite list of singletons. Uniform continuity of the limiting trace and the uniform error on late physical tails make the change arbitrarily small. Reserve one third of the perturbation budget each for the almost period, the prime-to-limiting error and this rounding. For sufficiently late translates the rounding also costs less than 0.0005 in log height, which preserves the separation and completion bounds above. The good centres even have positive lower uniform density on the log axis on sufficiently high tails. Take the almost periods with half the preceding error. Uniform continuity allows every sufficiently late such period to be moved by $`|h|\le\delta`$ for one fixed $`\delta\gt0`$, while retaining all targets. Every sufficiently late interval of one fixed length contains a whole good subinterval of length $`2\delta`$ by relative density. Disjoint blocks give the asserted positive lower density. This conclusion concerns centres of simultaneous working windows and the whole closed target box.

### Verifying a continuous arc with finitely many operations

At a scale $`z`$, split the prime norms at $`z/2,z,2z`$. On the two halves of the tent the numerator in the original weight is $`2-z/n`$ or $`2z/n-1`$. Two cumulative sums of the coefficients and coefficients divided by $`n`$ evaluate the whole prime sum exactly up to outward rounding. Between consecutive knots its continuation is $`Ue^{-s/2}+Ve^{s/2}`$ in the local log displacement $`s`$. A first-derivative jump of magnitude $`J`$ at $`s_k`$ changes the extrapolation by at most $`2|J|\sinh(|s-s_k|/2)`$. At the three moving cuts the absolute jump multipliers per $`|\log(\ell)u_k(\ell)|/\sqrt n`$ are $`2/B,4\cosh(a/2)/B,2/B`$. The physical coefficients have the same three geometric multipliers. A cumulative sum of absolute coefficients bounds all knots in each closed cell, including those at its endpoints. The Archimedean terms are interval evaluated on the whole cell, with the convergent series $`T=\sum_{j\ge0}s^{2j+3/2}/((2j+1)(2j+2))`$ and a geometric tail for small $`s`$.

The second proof computes the primitive characters in Conrey numbering on every residue, or determines the physical Frobenius by factoring defining polynomials. It forms secants between consecutive nodes of a different grid of step $`h=1/16000`$. A continuous function with smooth second-derivative bound $`M`$ and total absolute slope jumps $`J`$ on that cell differs from its secant by at most $`Mh^2/8+Jh/4`$. Integrating the absolute Dirichlet Green function gives the first term. Its pointwise maximum $`h/4`$ treats the atoms of the second derivative. Between knots the prime contribution satisfies $`P''=P/4`$. On a fixed arc where $`z_0\le e^x\le z_1`$ its absolute size is bounded by $`(2z_1+1)8\log(2z_1)g(0)/\sqrt{z_0/2}`$. At scale multiplier $`k\in\{1/2,1,2\}`$ a node cell can contain at most $`1+kz_1 h`$ integer knots by the mean value theorem. Each has log coefficient of magnitude at most $`8\log(2z_1)/\sqrt{kz_0}`$, before the jump multiplier above. The Archimedean term satisfies $`|\Gamma_i''|\lt100`$ on all nonconstant arcs and their two shifts. There $`e^{-|t|}\lt1/64`$ at the three gamma arguments, the conductor tent vanishes, and the differentiated series obeys $`|T''(t)|\le2s^{3/2}/(1-s^2)`$. The singleton is evaluated directly at its specified scale. The separate secant bounds at the two shifts enclose the move interpolation error. Finally the second program tests all sixteen corners of each target cell against all five functionals. Both methods use outward real balls and compute the barycentric coordinates from their own evaluations. [Closed-cell implementation](three_profiles.py) · [Independent whole-arc recheck](three_recheck.py)

Almost every target in this box requires at least three starts even when all positive real heights are allowed. On each bounded range away from zero, the physical score curve is locally Lipschitz by the explicit finite tent formula. The two-start averages are the image of two heights and one weight. The Lipschitz image of each bounded three-dimensional parameter box has zero four-dimensional volume by subdivision into small cubes. A countable union covers all choices of positive heights. Combined with the three-start existence result this gives the exact minimum for almost every target. [Earlier dimension argument](#how-many-observations-are-necessary)

For almost every target there is a further restriction on the coefficients. With any fixed three weights the three observation heights give only three real parameters, hence the image again has volume zero. A countable union covers all algebraic triples of weights, so at least two of the three weights must be transcendental for almost every target, since they sum to one. If the singleton is restricted to a countable set of anchors, only two heights can vary. Algebraic dependence over the rationals between any two weights confines that pair to a countable union of real plane algebraic curves. Their regular parts have countably many locally Lipschitz one-dimensional charts by the implicit function theorem, and each squarefree curve has only finitely many singular points. The image after including the two heights still has dimension at most three. Two algebraic starting heights in a triple leave only one height and two weights variable, another countable union of three-dimensional images. Hence almost every target needs two transcendental heights whenever the integer singleton is used. Consequently almost every target in the initial physical box requires algebraic independence of every pair of weights in any representation with one of the tabulated anchors. This dimension statement concerns almost every real target, and gives no classification for individual rational targets.

### Recovering the starts by one quartic

All statements in this section concern direct prime counts in the finite window. They use no zero locations. Write $`K_0=(\overline{\mathbf Q}\cap\mathbf R)(\log\ell:\ell\le3{,}000{,}000\text{ prime})`$. Add the two ones to the centred level coordinates and call the resulting target $`P`$. For every target in the box whose four coordinates are algebraically independent over $`K_0`$, the three starts and their common weights can be taken in **one extension of degree at most four** over $`K_0(P)`$. This holds for almost every target since $`K_0`$ is countable. The integer start is prescribed by the lookup. Degree four is attained generically on the particular local inverse described below.

On each open prime-norm band the raw four-vector at a variable start $`y`$ is $`A\sqrt y+B/\sqrt y`$, with $`A,B\in K_0^4`$. At an integer start its four-vector $`P_0`$ lies in $`K_0^4`$. Let $`\lambda`$ be its weight, and $`b,c`$ the weights of the two variable starts $`y,z`$. Set $`u=b\sqrt y`$, $`v=b/\sqrt y`$, $`s=c\sqrt z`$ and $`t=c/\sqrt z`$. For the two pairs of band vectors $`A,B`$ and $`C,D`$ the four measurements say

```math
P=\lambda P_0+uA+vB+sC+tD,\qquad
1-\lambda=\sqrt{uv}+\sqrt{st}.
```

For an algebraically independent target neither variable height can be a prime knot. A knot is algebraic, so together with the fixed anchor it would leave only one height and two freely varying weights to generate four independent observations. The same parameter count excludes any algebraic second height. The four-row matrix of the five vectors $`P_0,A,B,C,D`$ has rank four, since otherwise it would give a nonzero linear relation on $`P`$ over $`K_0`$. Its linear solutions are an affine line over $`K_0(P)`$. Use any free coordinate $`h`$ to parametrize it. All of $`\lambda,u,v,s,t`$ are then affine in $`h`$. The mass equation gives the polynomial

```math
Q(h)=\big((1-\lambda)^2-uv-st\big)^2-4uvst.
```

Its degree is at most four. It cannot vanish identically on this line at a solution with positive weights and both moving heights interior to their bands. Indeed, in a small interval around that solution the four factors of $`Q`$ are $`1-\lambda\pm\sqrt{uv}\pm\sqrt{st}`$. At the solution just one vanishes. If $`Q`$ vanished identically, the positive mass identity would persist on an open subinterval. Pick a parameter there transcendental over the countable field $`K_0(P)`$. Recover the two positive weights and heights by

```math
b=\frac{(1-\lambda)^2+uv-st}{2(1-\lambda)},\qquad
 c=1-\lambda-b,\qquad y=u/v,\qquad z=s/t.
```

The resulting $`\sqrt y,\sqrt z,b,c`$ generate a field of transcendence degree at most four over $`K_0`$. That field contains $`P`$ and the free coordinate, since the original five linear unknowns are rational expressions in those four parameters. This would give transcendence degree at least five. Thus $`Q`$ is a nonzero polynomial and all five unknowns, both heights and the three weights belong to $`K_0(P)(h)`$, of degree at most four. Positivity in the unsquared mass equation selects the appropriate branches. [Elimination identities in Lean](JointOscillation.lean) · [Construction at a specified target](three_inverse.py)

### An exact balanced example and the degree of its local inverse

Take raw target $`P=(1,1,0,0)`$. On timing row one the integer start is **13,471**. The other two starts are approximately **2,457.452761630442** and **61,460.713123523080**. In the order movable, integer, movable their weights are approximately **0.246867572222**, **0.394908054882** and **0.358224372897**. Both weighted midpoint levels equal one and both weighted moves vanish exactly. All three weights exceed $`1/5`$. The 18-place enclosures fix each displayed start and weight more accurately. [Exact brackets](evidence/three-inverse.txt)

Here $`M=[A\ B\ C\ D]`$ is invertible. Use $`\lambda`$ itself as free parameter and put $`(u,v,s,t)^{\mathsf T}=M^{-1}(P-\lambda P_0)`$. The primary verifier obtains the four columns by collecting prime coefficients on each side of the tent. It checks that all active cuts avoid any change of prime support on the two stated bands, and that the coefficient of $`\lambda^4`$ is nonzero. The quartic and its derivative have certified signs on the initial isolating interval, with opposite endpoint signs. Bisection isolates one simple real root, where all required factors have the positive signs for the original equations. [Quartic replay](three_inverse.py)

The second verifier factors the defining polynomials, evaluates the original tents at two values in each band, and interpolates the linear function obtained by multiplying each raw reading by $`\sqrt y`$. On the final four-dimensional box its preconditioned map $`x\mapsto x-N(F(x)-P)`$ is a strict contraction into the interior, where $`N`$ is the verified nonsingular dyadic midpoint of the inverse derivative. The contraction norm is below $`10^{-6}`$. This proves uniqueness within that box for the physical equations without using the elimination polynomial. The same computation checks all targets with $`\|P-(1,1,0,0)\|_\infty\le10^{-4}`$ simultaneously. Take the midpoints of the reported 18-place brackets for the two moving heights and their two weights, and give them respective radii $`1/20,2/5,10^{-5},10^{-5}`$. Every perturbed map takes this larger box strictly into itself, with weighted sup-norm Lipschitz constant below $`1/2`$. Each such target has a unique solution **in this common parameter box** with anchor 13,471. Other timing rows can give other representations. [Independent interval contraction](three_inverse_recheck.py) · [Transcript](evidence/three-inverse.txt)

The quartic on this branch is generically irreducible. Over $`K_0(U,V,S,T)`$ with four indeterminates, $`UV`$ and $`ST`$ have independent square classes by the valuations at $`U`$ and $`S`$. The four signs of $`\sqrt{UV}+\sqrt{ST}`$ give distinct conjugates. The monic polynomial $`((1-h)^2-UV-ST)^2-4UVST`$ is therefore irreducible, also in the polynomial ring $`K_0[U,V,S,T,h]`$ by Gauss’s lemma. Write $`g=M^{-1}P_0`$ and apply the invertible linear substitution $`(U,V,S,T)=H-hg`$ in five variables, where $`H=(H_1,H_2,H_3,H_4)`$. The transformed polynomial stays irreducible. Its fourth-degree coefficient is the constant $`(1-g_1g_2-g_3g_4)^2-4g_1g_2g_3g_4`$, certified nonzero for these prime bands. Gauss’s lemma now gives irreducibility as a quartic over $`K_0(H_1,H_2,H_3,H_4)`$. For any algebraically independent target in the displayed neighbourhood, $`H=M^{-1}P`$ has independent coordinates as well. Its local integer-start weight has exact degree four over the field of its four observations and $`K_0`$. All other solution coordinates lie in that same degree-four field. This sharpness is for the specified local branch. The rational centre has a certified simple quartic root, with no assertion of irreducibility at that specialization.

### A quartic for every target, including the boundary

Let $`P`$ be **any** target in the full closed four-dimensional box of the three-start theorem, and let $`E=K_0(P)`$. There is a representation on its prescribed timing row in which the two moving heights, their positive square roots and all three weights belong to a single real extension of $`E`$ of degree at most four. The two movable starts and the integer from the lookup retain all the earlier timing and weight bounds. This includes every rational target, each corner, and every target with algebraic relations among its four coordinates.

Fix one closed target cell containing $`P`$. The target intervals in this cell all have positive length. For small positive $`t`$, perturb coordinate $`i`$ by adding $`\sigma_i t^{5^{i-1}}`$, where $`\sigma_i=-1`$ at the upper endpoint of that coordinate interval and $`\sigma_i=1`$ elsewhere. The target $`P(t)`$ stays in the same closed cell and has the same prescribed integer. If $`G`$ is any nonzero polynomial over $`E`$ of total degree at most four in the target coordinates, then $`G(P(t))`$ is a nonzero polynomial in $`t`$. Indeed translation by $`P`$ and the four signs give an invertible affine substitution. Distinct remaining multiindices $`\alpha`$ produce distinct exponents $`\alpha_1+5\alpha_2+25\alpha_3+125\alpha_4`$ by uniqueness of base-five digits. There are at most 70 monomials and the exponent never exceeds 500. [Digit identity and limiting root in Lean](JointOscillation.lean)

The two compact moving arcs of the chosen row meet only finitely many prime-norm bands. Include their closures, since the prime sums agree at every common knot. For an ordered pair of such bands let $`N=[P_0\ A\ B\ C\ D]`$ be its constant four-by-five matrix over $`K_0`$. Write $`v=(\lambda,u,v_2,s,v_4)`$ for the five unknowns, with $`v_2=b/\sqrt y`$ and $`v_4=c/\sqrt z`$. The four observations give $`Nv=P(t)`$. If $`N`$ has rank less than four, a nonzero linear form annihilates its image. The perturbed target can meet that image for only finitely many $`t`$.

For each rank-four pair, choose any free coordinate $`h`$ whose four complementary columns are independent. The solution $`v=L(P(t),h)`$ is linear in the five arguments and has coefficients in $`K_0`$. Define its mass polynomial before specialization by

```math
G(P,h)=\bigl((1-\lambda)^2-u v_2-s v_4\bigr)^2-4u v_2 s v_4,
\qquad (\lambda,u,v_2,s,v_4)=L(P,h).
```

It has total degree at most four. The constant-in-$`h`$ coefficient $`G(P,0)`$ is a nonzero polynomial in $`P`$, since $`L(0,0)=0`$ gives $`G(0,0)=1`$. The perturbation makes $`G(P(t),h)`$ nonzero as a polynomial in $`t,h`$ for every rank-four pair. All but finitely many positive perturbation values also avoid its finitely many exceptional zero specializations and all the rank-deficient images.

Take a sequence of such perturbations tending to zero. The three-start theorem for the fixed chart supplies a solution at each one. Some fixed ordered pair of band closures occurs on a subsequence. Their matrix has rank four. By compactness pass to a further subsequence on which all the heights and weights converge. Its limit stays on the two closed arcs, satisfies the four original prime equations by continuity, and retains positive weights at least $`77/30000\gt1/400`$. The strict timing bounds hold throughout the closed arcs. The free coordinate converges to some finite $`h_0`$ because every linear unknown is continuous in the heights and weights and the heights are bounded away from zero.

In the single polynomial $`G(P(t),h)`$ divide by the largest common power of $`t`$ in its coefficients and then set $`t=0`$. The resulting polynomial $`g(h)\in E[h]`$ is nonzero and has degree at most four. Before specializing, the normalized equation vanishes along the chosen sequence. Continuity gives $`g(h_0)=0`$. If the original equation $`G(P,h)`$ was identically zero at this particular target, every surviving coefficient involves a positive power of at least one target increment. Total degree four then gives $`\deg g\le3`$ in this case. All five linear unknowns at the limit belong to $`E(h_0)`$. The positive formulas for $`b,c`$ above have nonzero denominator $`2(1-\lambda)`$. They recover the two square roots as $`u/b,s/c`$ and the heights as $`u/v_2,s/v_4`$. This proves the degree bound without an exceptional target set.

This also describes a finite list of candidate equations for a given target. On each rank-four pair of bands, apply the displayed substitution and read the first nonzero polynomial coefficient in ascending powers of $`t`$. A real root represents the target only if it gives the original positive weights, the assigned prime bands and the prescribed timing. At least one pair passes. The general argument establishes existence and a finite exact procedure given exact sign and equality tests in the coefficient field. The finite prime sums and the 23-row covering certificate supply its hypotheses.

At the 16 vertices the replay isolates an explicit positive witness at each one, with every weight exceeding $`1/100`$. The table gives each movable start and all three weights as an open interval with adjacent 18-decimal endpoints. Its integer entry is exact. In each row the four targets are the two raw midpoint levels and the two moves divided by $`\log2`$. The primary calculation solves a quartic on the indicated open bands, tests the unsquared signs and checks the timing row. A second implementation factors the defining polynomials and contracts all four original equations, with norm below $`10^{-6}`$ on every box. [All 16 corner enclosures](evidence/corner-enclosures.tsv) · [Quartic isolation](target_box_corners.py) · [Independent direct equations](target_box_corners_recheck.py)

## The four conjugates of the integer-start weight

Use the two maximal open prime-norm bands from the balanced example above. They are the bands containing the two specified movable starts, with their order fixed. On them the raw four-vector has the two expressions $`A\sqrt y+B/\sqrt y`$ and $`C\sqrt z+D/\sqrt z`$. Put $`M=[A\ B\ C\ D]`$, take $`P_0`$ at the fixed integer 13,471, and keep the field $`K_0`$ defined above. All entries lie in $`K_0`$. Write $`H=M^{-1}P`$, $`g=M^{-1}P_0`$, and let $`h`$ denote the weight of the integer. The four unknown linear coordinates are $`H_i-g_i h`$. Define

```math
Q_1=(H_1-g_1h)(H_2-g_2h),\quad Q_2=(H_3-g_3h)(H_4-g_4h),
\qquad F(h)=\big((1-h)^2-Q_1-Q_2\big)^2-4Q_1Q_2.
```

Over the entire closed box $`\|P-(1,1,0,0)\|_\infty\le10^{-4}`$ the coefficient of degree four stays nonzero. The verifier tests opposite endpoint signs for $`F`$ on each of the following disjoint intervals. Direct evaluation of the factored expression keeps the same target enclosure in all factors.

| Interval for $`h`$ | Test in the original equations |
| --- | --- |
| $`(-1.81,-1.77)`$ | The integer-start weight is negative. |
| $`(0.38,0.41)`$ | The interval contraction gives a positive solution on the designated bands. |
| $`(0.77,0.81)`$ | Throughout the interval $`(1-h)^2-Q_1-Q_2\lt0`$. |
| $`(1.34,1.38)`$ | The integer-start weight exceeds one. |

A quartic with four distinct real roots has exactly one in each interval. On any positive solution the expression tested in the third row is twice the product of the two positive movable weights. For the second row the existence theorem above gives the solution, since the other three intervals are excluded. The recovery formulas give at most one pair of heights and weights for any admissible $`h`$. This proves uniqueness on the **whole ordered pair of maximal bands** with this anchor, throughout the target box. It also identifies three real extraneous roots produced by squaring. All statements so far use only finite prime arithmetic. [Interval tests](inverse_monodromy.py) · [Existence by contraction](three_inverse_recheck.py) · [Lean sign identity](JointOscillation.lean)

### A simply branched line in target space

Complexify the algebraic coefficients of the fixed prime bands and set $`H=(T,1,1,1)`$. On a dense open subset let $`c`$ denote the signed weight of the second movable start and put $`z=c/(1-g_3h)`$. Then $`c^2=(1-g_3h)(1-g_4h)`$. Introduce the following polynomials in $`z`$.

```math
d=g_4-g_3z^2,\qquad n=1-z^2,\qquad p=d-g_2n,\qquad
m=d-n-(g_4-g_3)z .
```

Solving for the weights gives $`h=n/d`$, $`c=z(g_4-g_3)/d`$ and $`1-h-c=m/d`$. The first squared-weight equation is $`(1-h-c)^2=(T-g_1h)(1-g_2h)`$. It becomes

```math
T=\frac{a(z)}{b(z)}
=\frac{g_1np+m^2}{dp}.
```

Conversely, the two squared weights at a root of $`F`$ recover $`c=((1-h)^2+Q_2-Q_1)/(2(1-h))`$ on the open set where the denominators are nonzero. This recovers $`z`$. Thus the four parameters of this rational map and the four weights give the same function field and the same generic permutations. Exceptional fibres introduced by the rational coordinates are allowed in this argument. The complex parameters here describe algebraic continuation of the two fixed band formulas. Physical starting heights still belong to their specified positive real bands.

The determinant calculation encloses all coefficients from the original prime counts. The quartics $`a,b`$ have coprime homogenizations and four simple poles, all finite. For the binary quartic $`a-Tb`$, its discriminant as a polynomial in $`T`$ has degree at most six. At seven values $`T=1000k`$, $`0\le k\le6`$, the real-ball computation evaluates the Sylvester determinant divided by the leading quartic coefficient. Interpolation gives a sextic with nonzero leading coefficient and nonzero resultant with its derivative. The binary discriminant handles a root passing through infinity in $`z`$. It follows that the map has six simple branch values, with exactly one ramified point of index two above each. The two working precisions are 448 and 576 bits. [Resultants](inverse_monodromy.py)

The second computation factors the defining fields at the rational primes and interpolates the two original tents in each prime band to recover $`g`$. It obtains the numerator of the rational map from five values of the two weight equations. Applied to $`a'b-ab'`$, Rouché’s theorem places one simple critical point in each of six disjoint discs of radius $`10^{-10}`$. The interval images of those discs are pairwise disjoint. Conjugation fixes the critical points in the two discs centred on the real axis. Their images satisfy

```math
-0.025\lt T_-\lt-0.024,\qquad 1375\lt T_+\lt1376 .
```

The other four critical images have imaginary part of absolute value exceeding 10,000. At $`T=0`$ four sign-changing intervals exhibit four real inverse parameters. Above $`T=\infty`$ the two quadratic factors of $`b`$ exhibit precisely two real poles. The number of real points of a regular real fibre is locally constant on the real projective target. The two critical values divide it into two components. The bounded interval $`(T_-,T_+)`$ has four real points in each fibre. Its complement outside the endpoints has two. At either endpoint there are three distinct real points, one of multiplicity two. These counts refer to the normalized rational parameter curve. [Independent complex-disc tests](inverse_branch_recheck.py)

The sphere of source parameters with the inverse images of the branch values removed is connected. Hence the monodromy is transitive. A small loop about each simple branch induces one transposition. Such loops generate the fundamental group of the punctured sphere. The graph whose edges are the generating transpositions is therefore connected. Transpositions along the edges of a connected graph generate every permutation. This proves the group is $`S_4`$ on the line. Removing finitely many extra regular target values to avoid the exceptional birational coordinates preserves this permutation group. Each remaining loop also lies in the complement of the four-variable weight discriminant, so the full group contains $`S_4`$. It acts on four letters, giving equality. The geometric group over complex coefficients is contained in the arithmetic group over $`K_0(H)`$, giving the same equality there. [Pakovich, Theorem 2.2](https://www.math.bgu.ac.il/~pakovich/Publications/imax.pdf) · [Formal transposition deduction](JointOscillation.lean)

Now take any target in the above real box with four coordinates algebraically independent over $`K_0`$. Evaluation identifies the rational function field in the indeterminates with $`K_0(P)`$ since $`M`$ is invertible. The unique admissible weight has four real conjugates. The field containing all four has degree 24 over $`K_0(P)`$. A tower of square-root adjunctions has normal closure of degree a power of two, by adjoining the conjugates of each square root successively. The normal closure of the physical weight has degree 24, so that weight cannot belong to such a tower. This holds for almost every target in the box since $`K_0`$ is countable. The rational balanced target itself has a certified quartic and a unique physical root on these bands. Its specialized Galois group remains undetermined.

### The two elliptic resolvents

On the preceding $`T`$-line choose an unordered pair of the four $`z`$-roots. Forgetting which of the pair and its complement was chosen gives a map to the curve of partitions into two complementary pairs. Use smooth projective normalizations throughout. These connected covers have degrees six and three over the target line since $`S_4`$ acts transitively in each case. In the action on six pairs a transposition fixes two pairs and exchanges the remaining four in two cycles. In the action on three partitions it fixes one and exchanges two. With six simple branch points, Riemann–Hurwitz gives genus one for both curves. The complement map between them has degree two and is unramified. Above $`T=\infty`$ the two quadratic factors of $`b`$ give two distinct $`K_0`$-rational unordered pairs and one partition. Take one of the pairs and its image as origins. The complement quotient is then a degree-two isogeny over $`K_0`$, with the other pole pair as the nontrivial kernel point. For comparison, the ordered distinct-pair curve has degree 12 and genus four, and the Galois closure has degree 24 and genus thirteen. Their ramification contributions at each branch value are five and twelve, respectively. This is the classical trigonal construction applied to the specified prime pencil. [Donagi and Smith, section III.1](https://doi.org/10.1007/BF02392458) · [Ordered pairs](https://www.math.bgu.ac.il/~pakovich/Publications/imax.pdf)

Both curves have explicit plane cubics over $`K_0`$. Write $`a_i,b_i`$ for the ascending coefficients of $`a,b`$. For an unordered pair $`z,w`$ put $`s=z+w`$, $`p=zw`$. The divided difference of $`a(z)b(w)-a(w)b(z)`$ by $`z-w`$ is

```math
\sum_{0\le j\lt i\le4}(a_i b_j-a_j b_i)\,p^j U_{i-j-1}(s,p)=0,\qquad
(U_0,U_1,U_2,U_3)=(1,s,s^2-p,s^3-2sp).
```

Generically its points recover the two roots as the unordered roots of $`X^2-sX+p`$. For a partition use $`S=1/T`$, set $`c_i=Sa_i-b_i`$ and $`v=c_4(z_1z_2+z_3z_4)`$ for the chosen partition. The three generically distinct pair products summed this way yield the second cubic

```math
v^3-c_2v^2+(c_3c_1-4c_4c_0)v
+4c_4c_2c_0-c_4c_1^2-c_3^2c_0=0 .
```

The first cubic contains $`(0,p_0)`$ and $`(0,p_1)`$, where $`p_0=-g_4/g_3`$ and $`p_1=-(g_4-g_2)/(g_3-g_2)`$ are the products of the two pole pairs. Their common image on the second is $`(0,-b_4(p_0+p_1))`$. Project either cubic along lines through its stated rational point. The remaining quadratic in the line coordinate has a branch discriminant that is a polynomial of degree at most four in the slope. The primary calculation encloses its coefficients and checks its binary discriminant is nonzero. For binary quartic coefficients $`(a,b,c,d,e)`$ in descending order use

```math
I=12ae-3bd+c^2,\quad J=72ace+9bcd-27ad^2-27b^2e-2c^3,\qquad
j=\frac{6912I^3}{4I^3-J^2}.
```

The resulting enclosures are

| Curve | Lower bound for $`j`$ | Upper bound for $`j`$ |
| --- | ---: | ---: |
| Unordered pair | 2159.399670638032703450 | 2159.399670638032703451 |
| Complementary partition | 1372.295119267165244835 | 1372.295119267165244836 |

The independent calculation starts from the separately reconstructed prime pencil. It reduces both quartics modulo $`X^2-sX+p`$ on the pair curve, uses the depressed quartic on the partition curve, and computes the cross-ratio of the four certified projection branch points. The two resulting intervals agree with the invariant calculation. They also enclose zero on evaluation of the classical level-two modular equation. The complement quotient proves the exact isogeny. [Both cubics and invariant calculation](inverse_resolvent.py) · [Independent cross-ratios](inverse_resolvent_recheck.py) · [Binary quartic invariants](https://arxiv.org/abs/math/0610318) · [Classical modular polynomial](https://math.mit.edu/~drew/ClassicalModPolys.html)

The first invariant exceeds 1728. The second lies strictly between zero and 1728. In a short real Weierstrass model these ranges force positive and negative elliptic discriminant, respectively. Consequently the unordered-pair curve has two real components and the partition curve has one. The isogeny maps each real component bijectively to the target circle. Indeed the homomorphism restricted to the identity component is open and surjective, as its derivative is invertible. Both components map onto the circle, and the full map has degree two. The nontrivial real kernel point lies on the other component. All these objects and enclosures come from finite sums in the two actual fields.

## One start prescribed before all four targets

Use the physical vector $`W(y)=(L_A(y)-1,L_B(y)-1,D_A(y),D_B(y))`$ with the two common width-two observations defined above. Here $`D_i`$ is the change divided by $`\log2`$. The centred target box is $`B=[-10,10]\times[-8,8]\times[-7/2,7/2]\times[-11/5,11/5]`$. Put $`I=[13450,13490]`$. For every $`k\in I`$ and every $`P\in B`$ there are two further starts $`y,z`$ and three weights $`a,b,c\gt1/400`$ summing to one with

```math
aW(k)+bW(y)+cW(z)=P .
```

The chronological starts are separated by factors greater than 2.017. They all exceed 1,155 and all finish below 583,866 upon doubling. The largest finish divided by the smallest start is less than 148. The prescribed start can be chosen before observing any coordinate of the target. This gives 41 consecutive interchangeable integer choices as well as the intervening real heights. The weights, the other two starts and the timing row may depend on the target. Within each of the certified target cells the same two arcs work for every prescribed height and every target in that cell. The proof counts the physical prime ideals with the omissions and signs in the [trace calculation](#applying-the-criterion-to-prime-ideals).

**Moving the prescribed start.** Each row supplies two closed arcs in the coordinate $`x=\log(2\sqrt2\,y)`$. For its affine coordinates $`\lambda_0,\dots,\lambda_4`$, use the five vertices $`W(13471)`$ and the ordered endpoints of the arcs $`A_1,A_2`$. The inverse is certified separately on each row. Let $`J=W(I)`$, and write $`W(A_j)`$ for the image of the arc in physical starting heights. Define

```math
m_0=\max(0,\sup_{W(A_1)\cup W(A_2)}\lambda_0),
```

```math
m_i=\max(0,\sup_{J\cup W(A_2)}\lambda_i)\quad(i=1,2),
```

```math
m_i=\max(0,\sup_{J\cup W(A_1)}\lambda_i)\quad(i=3,4).
```

The certificate proves $`\lambda_i(P)\gt m_i+1/200`$ for all five coordinates, for every point of the assigned closed target cell. On all three pieces every affine coordinate is less than two. To apply the [simplex face argument](#three-common-starts-for-the-four-measurements), replace the singleton vertex in the grouped map by $`W(k)`$ and retain its original affine coordinate system. On face zero the singleton has zero mass. On every other face the endpoint pinned by that face has coordinate zero and the other arc and the actual singleton are bounded by $`m_i`$. Brouwer's theorem gives an exact solution. For each of the three groups take its designated coordinate $`i=0,1,3`$. A group of mass $`w`$ contributes less than $`2w`$ there. The rest contributes at most $`m_i`$, since $`m_i\ge0`$. The strict target clearance forces $`w\gt1/400`$. This argument works for both endpoints of the prescribed interval and for all shared cell faces. [Mass inequality in Lean](JointOscillation.lean) · [Direct verifier](prescribed_anchor.py)

**The finite cover.** Split each level axis into unit intervals and each move axis into halves. This gives 1,280 closed boxes. A zero in the subdivision table splits all four axes in half. Read children in lexicographic order with the lower half first. Each positive entry selects one of 41 timing rows for a closed leaf. The 1,415 leaves cover the whole target box. The anchor interval is contained strictly between 10.546 and 10.550 in centred log height. The moving-arc endpoints and this outward bound give a smallest gap of at least 0.702, a total log span of at most 4.304 per row, and common endpoints within $`[8.092,13.624]`$. Exponentiating gives the printed timing bounds.

One calculation builds the two curves using quadratic residues and signed norm sums. It evaluates complete intervals using the prime-tent derivative jumps. For each affine gradient it subtracts the exact rectangle cost $`\sum_j |u_j|r_j`$ from the value at its centre. It bisects a parameter interval until all the required strict affine inequalities hold. The other calculation factors defining polynomials to determine Frobenius degrees. On each arc it uses nodal secants with a rigorous interpolation error and visits all 16 corners of every leaf. It checks the prescribed-height endpoints directly at their exact values. Each method constructs its own inverse affine map and checks the whole closed prescribed interval. [Subdivision data](evidence/prescribed-anchor-tree.tsv) · [Independent checker](prescribed_anchor_recheck.py) · [Whole-interval interpolation bound](#verifying-a-continuous-arc-with-finitely-many-operations)

### Connected families of exact solutions

For any closed leaf $`Q`$ use the parameter space $`X=I\times Q`$. Apply the normalized simplex self-map from the [face proof](#three-common-starts-for-the-four-measurements) using $`(k,P)\in X`$. It is continuous in the parameters and the simplex variable, including faces of zero group mass. The parametric fixed-point theorem gives a connected component of its fixed-point set projecting onto all of $`X`$. At each point on this component all group masses exceed $`1/400`$. Reading the two moving starts from the simplex is continuous there. Its image is a connected set of physically exact solutions over the entire five-dimensional parameter rectangle. [Solan and Solan, parametric fixed-point theorem](https://arxiv.org/abs/2210.16369)

The full physical solution set on these compact arcs, with weights at least $`1/400`$, is compact and semialgebraic over the reals. There are only finitely many cuts from prime norms. On each closed prime band the centred vector has the form $`d+u\sqrt y+v/\sqrt y`$ with constant vectors and positive $`y`$. Square roots can be included as positive auxiliary variables. Two formulas agree at their common knots. The physical solution set therefore has a connected semialgebraic component projecting onto the whole parameter rectangle. Connected semialgebraic sets are semialgebraically path connected. For a *fixed* target apply the same argument with parameter space $`I\times\{P\}`$. This supplies a path of exact representations between the two endpoint heights with that target fixed. The path may turn back in the anchor coordinate. [Basu, Pollack and Roy, Theorems 5.22 and 5.23](https://perso.univ-rennes1.fr/marie-francoise.roy/bpr-ed2-posted3.pdf)

The [all-target specialization argument](#a-quartic-for-every-target-including-the-boundary) also applies at each chosen $`k`$. Use its closed leaf in the perturbation and hold $`k`$ fixed. The column $`W(k)`$ is defined over $`K_0(\sqrt k)`$ because its finite prime band has the displayed form. Here $`K_0`$ contains the real algebraic numbers and the prime logarithms through three million. Both closed moving arcs belong to the covering row for every nearby target in this leaf. The rank-four elimination and specialization therefore yield a solution in a real extension of $`K_0(P,\sqrt k)`$ of degree at most four, including both moving square roots and all weights. At an algebraic prescribed height the ground field is just $`K_0(P)`$.

### Intervals of prescribed heights beyond every starting point

For this asymptotic statement assume every nontrivial zero of each of the sixteen primitive Dirichlet factors lies on the critical line. Alpöge's announced primitive Dirichlet theorem supplies that input. In the finite physical theorem every assertion about the original interval and the best approximation with two fixed starts follows directly from the prime counts. [Announcement](https://alpo.ge/remark.pdf) · [Upstream formal checks](https://math.mit.edu/~drew/grh/)

There is a finite $`C`$ with the following property. For every $`H\ge1`$ there is an interval of prescribed starts $`[T,(1349/1345)T]`$ with $`T\gt H`$. For **each real height** in that interval, chosen in advance, three actual doublings realize **every target** in the same four-dimensional box. The two additional starts and all three weights can depend on the target. Every start is strictly after $`H`$, all three doublings finish before $`CH`$, adjacent starts have ratio greater than 2.017, and the last finish is less than 148 times the earliest start. Every weight is greater than $`1/400`$. In particular each such interval contains at least $`\lfloor4T/1345\rfloor`$ consecutive admissible integer starts. The global waiting factor $`C`$ is existential.

Use a second cover on the limiting trace. Its 35 rows cover 1,430 closed leaves, with the same nominal anchor interval and the same reference height 13,471. In scaled template units all starts exceed 1,155 and all finishes are below 586,206. Each leaf has one row covering every anchor in the interval. One verifier evaluates the full prime-power expression and its Archimedean term. The other computes the primitive Conrey characters and verifies independent complete secants for both moves. The same face margins, gradient bound $`3/20`$, coordinate ceiling two and separation and span bounds apply. The two verifiers cover the complete joint target and anchor intervals. [Limiting arc endpoints](evidence/limiting-anchor-charts.tsv) · [Closed-cell table](evidence/limiting-anchor-tree.tsv) · [Two-method receipt](evidence/prescribed-anchor.txt)

Under the stated zero input the four-coordinate limiting vector is uniformly almost periodic, and the physical and limiting vectors differ by $`o(1)`$ uniformly on high log-height tails, by the [prime-power formula and correction](#the-limiting-curves-and-late-physical-starts). The certified inequalities are strict and the lists and parameter sets are finite and compact. Hence one positive uniform tolerance preserves the clearance $`1/200`$ and the ceiling two. Take sufficiently accurate common almost periods. All sufficiently large such translates preserve the entire face test after passing to the physical curve, uniformly in the chosen anchor as well as both moving arcs. At a pinned endpoint the original coordinate is zero and the translated error is controlled by the same uniform tolerance. The fixed-point proof applies with the original affine coordinates and its reserved strict margin. Translating by $`\tau`$ multiplies all heights by $`e^\tau`$. Thus every $`k\in[13450e^\tau,13490e^\tau]`$ works for all targets on that translated table.

The allowable large translates have bounded gaps. For all sufficiently large $`H`$ choose one in an interval of fixed length after $`\log(H/1155)`$. The template's strict lower bound puts every selected start beyond $`H`$. Its strict completion bound supplies one finite $`C`$ independent of $`H`$, the target and the chosen anchor. One fixed high working translate serves the bounded initial range after enlarging $`C`$. Moreover the set of good anchors contains a measurable set of **positive lower uniform density on the log-height axis**. The working intervals there have fixed length $`\delta=\log(1349/1345)\gt0`$. Every sufficiently high interval of one fixed length contains a whole such interval by relative density. Partitioning long intervals into blocks of that length proves the density assertion. The statement gives neither a numerical waiting factor nor a numerical density.

### The obstruction with two prescribed starts

Fix any two observations $`A,B`$ and let only the third start vary. On one of its bands write $`C(y)=d+u t+v/t`$ with $`t=\sqrt y`$. First suppose $`A\ne B`$. Subtract $`B`$ and pass to the quotient of four-space by the line parallel to $`A-B`$. The image of a weighted target there is

```math
w\,\overline{(d-B)}+wt\,\bar u+(w/t)\,\bar v .
```

If the three vectors are independent, these coordinates satisfy $`x_1x_2=x_0^2`$. Their inverse image is a proper quadric containing the line through the two anchors in its singular locus. If the three vectors are dependent the quotient image lies in a proper linear subspace. When $`A=B`$ the unquotiented differences already lie in the at most three-dimensional span of $`d-B,u,v`$. On bounded closed height intervals away from zero there are finitely many bands. The positive half-axis is a countable union of such intervals. Hence for every fixed pair of starts the attainable target set is contained in a countable union of proper quadratic or affine hypersurfaces. It has four-dimensional volume zero and cannot contain the target box. The assertion remains true for any countable list of fixed pairs. [Band identity in Lean](JointOscillation.lean)

### The sharp cost of fixing 2457 and 13471

Prescribe the two exact integer heights 2,457 and 13,471. Set the target to the origin of the centred coordinates, corresponding to two raw midpoint levels equal to one and two zero moves. Allow every triple of nonnegative weights summing to one. Allow the third height $`y\in[1,300000]`$ wherever the three doubled intervals have disjoint interiors. This includes touching endpoints. The unique optimizer has strict gaps, so the conclusion also holds for pairwise disjoint closed doubling intervals. The five quantities below lie in open brackets of width $`10^{-18}`$. Each displayed number is the lower endpoint. [Full enclosures](evidence/fixed-pair-balance.txt)

| Quantity | Lower endpoint |
| --- | --- |
| Third start | 61455.609822534600169779 |
| Minimum worst absolute coordinate error | 0.002274350938622389 |
| Weight at 2,457 | 0.246827861901412719 |
| Weight at 13,471 | 0.394992237138323180 |
| Weight at the third start | 0.358179900960264099 |

Write $`A=W(2457)`$, $`B=W(13471)`$ and $`s=(-1,1,1,-1)`$. On the closed band $`[61455.60,61455.62]`$ put $`C(y)=d+u t+v/t`$ as above. The candidate has weighted error $`e s`$ with all four coordinates active. The determinant condition reduces to

```math
\det(A,B,u,s)t^2+\det(A,B,d,s)t+\det(A,B,v,s)=0 .
```

The first verifier checks unchanged prime support, a nonzero quadratic coefficient, the sign change and a derivative of constant strict sign. Bisection gives the indicated root. The first three equations for the first two weights and $`e`$ have an invertible matrix with rows $`(A_i-C_i,B_i-C_i,-s_i)`$. Its fourth equation follows from the vanishing determinant by subtracting the third column from each of the first two columns of the original determinant and using the invertible minor. The three weights are positive. Independently, direct evaluation of the original tents at two points of the band reconstructs $`u,v`$. An interval contraction of all four original equations in the variables $`(y,a,b,e)`$ encloses the same four quantities, and the third weight follows from mass one. [Quadratic isolation](fixed_pair_balance.py) · [Direct contraction](fixed_pair_recheck.py)

For local optimality construct two linear functionals of dual norm one. Both take the same value on $`A,B,C(y_*)`$ at the candidate. The left functional has its second coordinate zero, the right its first coordinate zero, and all their other coefficients have the corresponding strict signs of $`s`$. Their values at $`s`$ are both one. Positivity and the candidate representation give their common value on the three observations as $`e_*`$. On the whole closed band the left functional is strictly decreasing on $`C(y)`$ and the right one is strictly increasing, by the derivative $`(u-v/y)/(2\sqrt y)`$. Thus for $`y\lt y_*`$ the left functional takes a value greater than $`e_*`$ on every mixture of positive third mass. The right gives the same conclusion for $`y\gt y_*`$. With zero third mass the 547/727 move coordinate is below $`-10`$. At $`y_*`$ equality for both supporting functionals forces every coordinate of the mixture to equal $`e_*s_i`$, since together their supports contain all four coordinates. The invertible minor gives the unique weights.

For global optimality the admissible starts are exactly $`[1,1228.5]`$, $`[4914,6735.5]`$ and $`[26942,300000]`$. Outside the interior of the local band the certificate covers them with 300 closed intervals, with deliberate overlap at the local endpoints. Each interval has an integer covector $`n`$ satisfying

```math
\min\{n\cdot A,\ n\cdot B,\ \inf_y n\cdot C(y)\}
    \gt 0.00227436\,\lVert n\rVert_1 .
```

Every nonnegative mixture then has worst absolute coordinate error greater than 0.00227436, strictly above $`e_*`$. The first implementation checks the binary partition and the whole prime-tent intervals. The second checks exact rational adjacency of the partition, obtains the signs independently by polynomial splitting, and validates every separating inequality with 1,621 secants including all possible prime knots. On a secant of log length $`h`$ the prime part satisfies $`P''=P/4`$ between knots. A bound $`M`$ for its absolute value and total absolute slope jumps $`J`$ give interpolation error at most $`Mh^2/32+Jh/4`$. These are the physical terms of the earlier [secant estimate](#verifying-a-continuous-arc-with-finitely-many-operations), so the estimate applies down to height one here. The 300 covectors and their exact interval addresses suffice for replay. [Global certificate](evidence/fixed-pair-separators.tsv) · [Independent recomputation](fixed_pair_recheck.py)

### One square root suffices for the best fit to any target

Fix two starts and allow the third to vary over any nonempty finite union of compact intervals of positive heights. Use a coefficient field $`E\subset\mathbb R`$ containing the target coordinates, the prime-band coefficients, the two fixed observations and the square roots of the allowed interval endpoints and band cuts. For the two integer starts and the range above one can take $`E=K_0(P)`$. For **every** target $`P`$ there is a minimizer of the worst absolute coordinate error in a single real extension of $`E`$ of degree at most two. This one extension contains the square root of the third height, the height, all three weights and the optimal error. The claim includes targets with zero optimal error and targets on boundaries or singular bands.

A minimizer exists by compactness. On each band, use $`t=\sqrt y\gt0`$ and $`C(t)=d+ut+v/t`$. At a fixed $`t`$ the problem is a linear program in three weights and the error $`e`$. Impose mass one, three nonnegative-weight inequalities and eight signed coordinate inequalities $`\pm(aA_i+bB_i+cC_i(t)-P_i)\le e`$. Its optimal face is nonempty and bounded, since $`e`$ is fixed there and the weights lie in the simplex. At an extreme optimum there are three active inequalities linearly independent together with the mass equality. Indeed a nonzero feasible line through an optimal point would have constant objective there or improve it in one direction. There are at most $`\binom{11}{3}=165`$ choices of basis per band. The same argument applies when a weight is zero or more than three inequalities are active.

Only the column multiplying $`c`$ varies with $`t`$. Its entries have the form $`\alpha+\beta t+\gamma/t`$. The basis determinant and each Cramer numerator have that same form. The numerator of $`c`$ itself is constant because its replaced column contains the fixed right-hand sides. In the residual of any unused inequality, the only variable coefficient multiplies this constant numerator after clearing the determinant. Hence, on each fixed determinant sign, feasibility is determined by polynomial signs of degree at most two after multiplying by $`t`$. Ignore identically zero equations. In particular a nondegenerate basis gives an objective $`N(t)/D(t)`$ with both polynomials of degree at most two. All its weights have the same rational form.

Choose a valid optimal basis at a global minimizer. If it is at a band or allowed endpoint, use that endpoint. If it is a root of a nonzero feasibility or determinant sign polynomial, it has degree at most two. Otherwise the same basis remains feasible on an open interval containing it. A nonconstant objective has derivative numerator $`N'D-ND'`$ of degree at most two and the minimizer is a root. If the objective is constant, choose any rational point in that feasible interval. In each case the chosen basis is invertible at the selected point and rational recovery gives the promised common field. A singleton allowed interval is an endpoint case. [Degree cancellation in Lean](JointOscillation.lean)

**The quadratic degree is attained for almost every target near balance.** For the two fixed integers and the allowed third starts above, the optimum stays unique for all $`\|P\|_\infty\le4\cdot10^{-9}`$ in centred coordinates. Its four error signs stay $`s=(-1,1,1,-1)`$. In the target-dependent determinant use columns $`A-P,B-A,C(y)-A,s`$. The direct check gives a simple sign-changing root in $`[61455.60975,61455.60990]`$ and positive weights throughout the closed target box. Both supporting covectors for that root retain their signs and the required derivatives over the full local band. A separate contraction of the original four equations covers every such perturbed target. Outside the band, distance to each candidate changes by at most $`\|P\|_\infty`$. The balanced candidate itself gives an upper bound $`e_*+\|P\|_\infty`$ and the outside certificate gives a strict lower bound $`0.00227436-\|P\|_\infty`$. The archived bounds imply $`e_*+8\cdot10^{-9}\lt0.00227436`$. This proves global uniqueness throughout the closed neighbourhood.

In the quotient by the linear span of $`B-A`$ and $`s`$, take two independent linear coordinates $`\ell_1,\ell_2`$ and put $`\xi_j=\ell_j(P-A)`$. Write

```math
H_j(t)=\ell_j(u)t^2+\ell_j(d-A)t+\ell_j(v).
```

The determinant condition is $`\xi_1H_2(t)-\xi_2H_1(t)=0`$. The two verifiers compute the coordinates independently. Both quadratics have degree two, with nonzero resultant, and at least one has a nonzero linear coefficient. The polynomial $`H_2(T)-R H_1(T)`$ is linear in $`R`$ and primitive over $`K_0[T]`$ by coprimality. It is irreducible there and is irreducible of degree two in $`T`$ over $`K_0(R)`$ by Gauss's lemma. For four independent target coordinates one can take $`R=\xi_2/\xi_1`$ as one of four rational-function parameters, the others coming from $`\xi_1`$ and two complementary coordinates. Adjoining these independent parameters preserves the degree. The linear term in the generic quadratic is nonzero. Its equation expresses $`t`$ in $`K_0(P,t^2)`$, so the **optimal height itself** has degree exactly two over $`K_0(P)`$. Algebraically independent tuples have full measure here because $`K_0`$ is countable. The numeric certificate isolates the optimizer at the balanced centre without requiring irreducibility at that individual target. [Quadratic pencil and global separation](fixed_pair_balance.py) · [Independent contraction and resultant](fixed_pair_recheck.py)
