# Oscillation for every compact weight

## The measure on prime ideals

For either $`(p,q)=(547,727)`$ or $`(571,599)`$, set $`K=\mathbf Q(\zeta_{16})^+`$, $`F=K(\sqrt{pq})`$ and $`E=K(\sqrt{-p},\sqrt{-q})`$. The sign $`\epsilon_v`$ of a good prime ideal of $`F`$ is one for splitting in $`E/F`$ and minus one for inertia. The first-step theta trace jump has a fibre product equal to this sign. Omit places above 2 and the two conductors. Any additional fixed finite omissions are allowed. For a nonzero real $`w\in C_c((0,\infty))`$, put

```math
S_w(y)=-\frac1{\sqrt y}\sum_v\epsilon_v\log(Nv)\,w(Nv/y),
\quad H_w(z)=\int_0^\infty w(u)u^{z-1/2}\,du,
\quad B_w=\frac12 H_w(0).
```

Define $`F_w(x)=S_w(e^x)-B_w`$. Continuity and finite support give local boundedness. Before the constant is subtracted, the score vanishes for sufficiently negative $`x`$ and is $`O_w(e^{x/2})`$ for large $`x`$ by the prime ideal Chebyshev bound. [547/727 fibre identity](https://github.com/leecrossley/theta-race-short-shifts/blob/483c31f2df2132ed2165620ea55e20e969228bfb/RaceTrace.lean#L83) · [571/599 admissible host and trace](https://github.com/leecrossley/theta-host-split-inert-race/blob/8fca7cad37185e3d7951f0dffff5e2442d1dfa24/HostAdmissible.lean#L80)

Let $`\mathcal L`$ be the product of the eight odd primitive Dirichlet factors of $`\zeta_E/\zeta_F`$. They are the quadratic characters of conductors $`p,q`$, each twisted by the four powers of the even quartic character modulo 16. For $`\Re s\gt1`$ one has the continuation

```math
D(s)=-\sum_v\epsilon_v\log(Nv)(Nv)^{-s}
=\frac{\mathcal L'}{\mathcal L}(s)-\frac{\zeta_F'}{\zeta_F}(2s)+R(s),
```

where $`R`$ is holomorphic for $`\Re s\gt1/3`$. The square correction removes the unsigned squares in the logarithmic derivative; higher powers converge absolutely there, as do the higher powers of the base correction. Finite Euler terms are regular in this half-plane. Mellin integration gives the bilateral transform of the uncentred score as $`H_w(z)D(1/2+z)`$. Restricting to $`x\ge0`$ changes this by an entire function. Subtract $`B_w/z`$ for the centring. [Euler factors, square correction and residues](https://github.com/leecrossley/theta-host-split-inert-race/blob/joint-levels-and-moves-v1.0.0/results/joint-levels-and-moves/PROOF.md#residues-for-a-real-projection)

The finite real-axis enclosure excludes zeros on $`[1/2,1]`$ in all sixteen odd factors across the two fields. It verifies the character partial sums over a complete period and uses an explicit derivative bound to cover the entire segment. The Euler product handles the real axis above one. The base pole at $`s=1/2`$ has residue $`1/2`$ in $`D(s)`$, exactly cancelled in the centred transform. At $`z=i\gamma`$, $`\gamma\ne0`$, the base term is regular by the prime ideal theorem on $`\Re s=1`$. The resulting residue is $`m(\gamma)H_w(i\gamma)`$, where $`m(\gamma)`$ is the total multiplicity of the zeros at $`1/2+i\gamma`$ in the odd factors. Multiplicities add across characters. [Real-axis replay](real_axis.py) · [Derivative bound](https://github.com/leecrossley/theta-host-split-inert-race/blob/joint-levels-and-moves-v1.0.0/results/joint-levels-and-moves/PROOF.md#reproducing-the-numerical-certificate)

## A compact weight cannot cancel all the residues

Write $`p(v)=e^{v/2}w(e^v)`$, a nonzero continuous compactly supported real function. Then $`H_w(z)=\int p(v)e^{zv}\,dv`$ is entire of exponential type. If it vanished identically, integration against every polynomial would give zero by differentiating at zero. Uniform polynomial approximation on the support would then give $`\int p(v)^2\,dv=0`$, a contradiction. Choose a point $`z_0`$ where $`H_w(z_0)\ne0`$. On discs of radius $`2T`$ about that point the log maximum of its absolute value is $`O_w(T+1)`$. Jensen's formula therefore bounds the number of complex zeros in the concentric disc of radius $`T`$ by $`O_w(T+1)`$, including multiplicities. In particular only $`O_w(T+1)`$ distinct ordinates up to height $`T`$ can be cancelled.

Choose the primitive quadratic factor of conductor $`p`$. Bauer's Corollary 1 gives a positive proportion of simple critical-line zeros for this fixed character. The ordinary zero-counting formula gives $`\gg T\log T`$ distinct positive such ordinates. Thus infinitely many survive every nonzero compact weight, even when its integral is zero. No conjecture about zeros off the line or their linear independence enters the argument. [Bauer, Corollary 1](https://matwbn.icm.edu.pl/ksiazki/aa/aa93/aa9313.pdf#page=2) · [Dirichlet zero count](https://doi.org/10.1090/mcom/3599)

Here is the quantitative one-pole argument. If $`\limsup F_w`$ is finite, take any constant $`C`$ strictly above it and discard a bounded initial interval. The Laplace transform of $`C-F_w`$ then comes from a nonnegative function. It is regular on the positive real axis by the preceding calculation. Landau's theorem for a nonnegative Laplace transform forces its abscissa of absolute convergence to be at most zero. For $`\sigma\gt0`$ the absolute value of the transform at $`\sigma+i\gamma`$ is at most its value at $`\sigma`$. Multiply by $`\sigma`$ and take the limit to zero. The left side tends to $`|m(\gamma)H_w(i\gamma)|`$ and the right side to $`C`$. The bounded initial interval contributes zero to either limit. Apply the same argument to the negative, including the case of an infinite extremum. It follows that

```math
\limsup_{x\to\infty}\sigma F_w(x)\ \ge \sup_{\gamma\gt0}m(\gamma)|H_w(i\gamma)|\gt0,\qquad \sigma\in\{-1,1\}.
```

In particular two weights ordered on a tail after centring up to an error tending to zero must be equal. Apply the inequalities to their difference. Likewise a nonzero, compactly supported finite real signed measure of log shifts cannot erase a fixed score. A translate by $`t`$ replaces $`w(u)`$ by $`e^{-t/2}w(ue^{-t})`$ and multiplies its transform by $`e^{tz}`$. Integrating over the measure multiplies $`H_w`$ by its entire Laplace transform. The latter is nonzero by polynomial approximation for measures. Their product is nonzero and is again the transform of a compact continuous annular weight. [Diamond's oscillation theorem](https://personal.math.ubc.ca/~gerg/teaching/592-Fall2018/papers/1971.Diamond.pdf)

## What survives on sparse grids

Assume in this section that the distributional second derivative of $`p(v)=e^{v/2}w(e^v)`$ is a finite signed measure. Continuous, compactly supported, piecewise twice continuously differentiable profiles with bounded one-sided derivatives satisfy this condition when their second derivatives are integrable on the pieces. It includes the dyadic tents and their finite signed mixtures. Let $`d`$ be the supremum of $`\Re\rho-1/2`$ over zeros of the odd factors with $`H_w(\rho-1/2)\ne0`$. The preceding argument ensures $`0\le d\le1/2`$. In repeated-zero sums below $`\rho`$ runs over these factors with multiplicity.

Two integrations by parts in the distributional sense give $`H_w(z)\ll_w(1+|\Im z|)^{-2}`$ uniformly for $`-1\le\Re z\le1`$. To obtain the explicit formula first count all signed prime powers by $`\mathcal L'/\mathcal L`$. Invert on $`\Re z=1`$ and move to $`\Re z=-1`$, so the corresponding left line for the L-function is $`\Re s=-1/2`$. The odd primitive factors have no trivial zero in the region crossed. On the left line the functional equation and absolute convergence of the reflected Euler product give $`O(\log(2+|\Im z|))`$ for the logarithmic derivative. On horizontal lines escaping in both directions choose heights at distance $`\gg1/\log(2+T)`$ from the zeros in the adjacent unit interval. This is possible by the local zero count for the fixed factors. The partial-fraction formula there gives $`O(\log^2(2+T))`$, and the two horizontal integrals tend to zero. The new vertical integral is $`O_w(e^{-x})`$. This argument also follows by smoothing $`p`$, using its second-derivative mass to dominate the zero sum, and passing to the limit. [Primitive explicit formula](https://www.cs.uleth.ca/~kadiri/articles/explicit_zerofree_regions_for_dirichlet_l_functions.pdf) · [Local zero counts](https://doi.org/10.1090/mcom/3599)

Restoring ordinary prime ideals adds the unsigned squares. The prime ideal theorem and partial summation give the limit $`\int_0^\infty w(t^2)\,dt=B_w`$. Higher prime powers are $`O_w(y^{-1/6}\log y)`$ after normalization by summing the prime ideal Chebyshev bound across powers. A fixed exceptional place contributes $`O_w(y^{-1/2})`$ from its powers in a compact log window. Consequently

```math
F_w(x)=\sum_{\rho} H_w(\rho-1/2)e^{(\rho-1/2)x}+o(1).
```

This formula contains zeros throughout the strip. The absolute sum of the coefficients is finite. When $`d=0`$, all surviving terms with negative real shift tend to zero in sum by dominated convergence. The remaining series is uniformly almost periodic and uniformly continuous. Hence every mesh with $`y_{j+1}/y_j\to1`$ has the same two limiting extrema as the full continuous score.

For $`d\gt0`$, let the positive mesh be increasing and unbounded with $`y_{j+1}-y_j=O(y_j^{1-\eta})`$ for some fixed $`\eta\gt0`$. Every large $`x=\log y`$ is within $`O(e^{-\eta x})`$ of a log sample. If the displacement is $`h`$, $`0\lt|h|\le1/2`$, summing the absolute differences gives

```math
|F_w(x+h)-F_w(x)|
\ll_w e^{dx}|h|\log^2(2/|h|)+o(1).
```

Indeed each summand costs at most a constant times $`e^{dx}(1+|\Im\rho|)^{-2}\min(1,|h|(1+|\Im\rho|))`$. Split at $`1/|h|`$ and use the zero count $`O(T\log(2+T))`$ by partial summation. The prime-only remainder is uniformly small on the tail. The resulting rounding error on this mesh is $`O(e^{(d-\eta)x}(1+x)^2)+o(1)`$.

Take $`\max(0,d-\eta)\lt\omega\lt d`$. Some surviving pole has real shift greater than $`\omega`$. For either sign, an eventual bound $`\sigma F_w(x)\le C e^{\omega x}`$ would give a nonnegative tail with transform regular on the real interval $`(\omega,\infty)`$. Landau's theorem would force absolute convergence to the right of $`\omega`$, contradicting the pole there. Thus the continuous signed excursions exceed arbitrarily large multiples of $`e^{\omega x}`$. Rounding to a sample preserves this conclusion because the rounding error is of smaller order. Conversely the explicit formula bounds $`|F_w(x)|\ll_w e^{dx}+o(1)`$. Letting $`\omega`$ increase to $`d`$ gives, for either sign and also in the bounded case,

```math
\limsup_{j\to\infty}
\frac{\log(1+\max(0,\sigma F_w(\log y_j)))}{\log y_j}=d .
```

In particular both signs on any such mesh attain all the strict finite radii proved on continuous tails. This covers $`y_j=j^k`$ for every fixed real $`k\gt0`$ by the mean value theorem. The log grid of powers of two has constant spacing and is outside this mesh condition. [Earlier fixed-power argument for the ordinary-prime tent](https://github.com/leecrossley/prime-window-power-samples/blob/9af696c425c8dfff85416b2a27ff92a2af248d0b/PROOF.md)

For the normalized width-$`b`$ level the Mellin kernel, up to a nonvanishing exponential for the choice of centre, is

```math
G_b(z)=\frac{b+1-2\sqrt b\cosh(z\log b)}
 {2(\sqrt b-1)^2(1/4-z^2)} .
```

The numerator can vanish only at $`z=\pm1/2+2\pi i k/\log b`$. Thus all open-strip zeros with $`\Re z\gt0`$ survive this weight. A normalized dyadic difference adds the factor $`(e^{az}-1)/a`$, $`a=\log2`$, which has no zero in the open right half-plane. For either observable and either sign, an eventual one-sided bound on any of the above meshes is therefore equivalent to all eight odd factors having their strip zeros on the line. The functional equations reflect the strip, and the family contains its conjugates. In particular the cubes suffice. This is a logical equivalence whose proof uses no global zero-location theorem. [Width kernel](https://github.com/leecrossley/theta-host-split-inert-race/blob/joint-levels-and-moves-v1.0.0/results/joint-levels-and-moves/PROOF.md#residues-for-a-real-projection)

## From quadrature to prediction

At width two use $`a=\log2`$, $`r=\sqrt2`$, $`B_0=2r(r-1)^2`$ and $`w(u)=\max(0,\min(2-2/u,4/u-1))/B_0`$. Write $`S(y)=S_w(y)`$ and $`D(y)=(S(2y)-S(y))/a`$. The square bias of the level is one and that of the move is zero. A delay vector of length $`m`$ consists of the levels at $`y,2y,\ldots,2^{m-1}y`$, or the moves at those same starts. For real coefficients $`c_j`$ put $`P(Z)=\sum_{j=0}^{m-1}c_j Z^j`$. At a positive critical-line zero $`t`$, the absolute transform of the linear combination is $`A(t)|P(e^{iat})|`$, where

```math
G(t)=\frac{3-2r\cos(at)}{2(r-1)^2(1/4+t^2)},\qquad
A(t)=G(t)\quad\text{for levels},\qquad
A(t)=G(t)\frac{2|\sin(at/2)|}{a}\quad\text{for moves}.
```

The extra change of centre is a unit complex phase. Each strict change of sign selects a zero of multiplicity at least one. The proof permits coincident ordinates from different characters and higher multiplicities. Every nonzero projection is a nonzero compact weight because the polynomial is nonzero and so is the base Mellin transform. The continuous and sparse-mesh results above apply to it.

For a fixed field, observable and dimension, the certificate supplies positive integer masses $`q_l`$ at certified ordinates $`t_l`$. Let $`Q=\sum_l q_l`$, $`C=\sum_l q_l/A(t_l)^2`$ and $`T_k=\sum_lq_l\cos(kat_l)`$. Expanding the square gives

```math
\sum_l q_l|P(e^{iat_l})|^2
=Q\|c\|_2^2+2\sum_{k=1}^{m-1}T_k\sum_{j=0}^{m-k-1}c_jc_{j+k}
\ \ge\ \left(Q-2\sum_{k=1}^{m-1}|T_k|\right)\|c\|_2^2 .
```

Cauchy–Schwarz bounds each inner correlation in absolute value by $`\|c\|_2^2`$. With $`M=\max_l A(t_l)|P(e^{iat_l})|`$, the same sum is at most $`CM^2`$. The primary calculation checks $`Q-2\sum|T_k|\gt r_m^2C`$ for each of the displayed rational radii $`r_m`$. Therefore $`M\gt r_m\|c\|_2`$ whenever $`c\ne0`$. The one-pole inequality proves both directions with this strict margin. For each field this also places the closed Euclidean ball of radius $`r_m`$ in the closed convex hull of its centred delay vectors on every tail, by separation against real unit vectors. The starting height can depend on the projection. [Quadrature masses](evidence/quadratures.tsv) · [Numerical check](quadrature.py)

The independent calculation uses $`V(\theta)=Q+2\sum_{k=1}^{m-1}T_k\cos(k\theta)`$. Orthogonality gives the same weighted sum as the integral of $`V(\theta)|P(e^{i\theta})|^2`$ against normalized circle measure. At the midpoint of each true bracket it multiplies successive complex unit phases. At order $`k`$ it adds an error bound $`ka`$ times the radius of the bracket. On each of 32 closed subintervals of $`[0,\pi]`$ it evaluates the symbol at the midpoint and subtracts the radius times $`2\sum k|T_k|`$. Evenness covers the rest of the circle. All 896 interval inequalities give $`V(\theta)\gt r_m^2C`$. The table has 710 different positive brackets across both fields. Their endpoints are integers in units of $`10^{-14}`$. Each quadrature row contains the field index, level or move index, dimension, zero-based bracket row and positive mass. [Independent symbol calculation](independent.py) · [Bracket table](evidence/brackets.tsv) · [Endpoint sign check](check_zeros.py)

In an affine prediction of the last coordinate from its predecessors, the coefficient of that last coordinate in the error is one, hence $`\|c\|_2\ge1`$. The centred errors exceed $`r_m`$ in each direction on every tail, including a tail of any fixed power mesh. Adding an arbitrary intercept cannot put both errors inside $`[-r_m,r_m]`$ because two attained values are more than $`2r_m`$ apart. The numerical bounds concern fixed affine predictors. [Radius and affine-error deductions](CompactFilters.lean) · [Full transcript](evidence/check.txt)

## The fixed-window theorem

Fix either of the two fields in the prime-ideal definition above. For a real $`w\in C_c((0,\infty))`$ write $`F_w(x)=S_w(e^x)-B_w`$. Given $`n\ge1`$ weights $`w_j`$, set $`F(x)=(F_{w_1}(x),\ldots,F_{w_n}(x))`$. The weights are linearly independent over the reals if and only if there are $`r\gt0`$ and $`1\lt C\lt\infty`$ for which

```math
\{v\in\mathbb R^n:\|v\|_2\le r\}
\ \subset\ \mathrm{conv}\{F(\log y):H\lt y\lt CH\}
\qquad(H\ge1).
```

Every target on the left is a convex combination of at most $`n`$ actual vectors on the right. A linear relation among the weights places every centred vector in a proper linear hyperplane. It remains to prove sufficiency.

### The coefficients seen by a real direction

In a first step average $`F(x+t)`$ against a single even, nonnegative smooth probability density supported in $`[-\epsilon,\epsilon]`$. Its transform $`L(z)=\int e^{zt}\phi(t)\,dt`$ multiplies each $`H_{w_j}(z)`$. The log profiles become smooth and compact, since logarithmic translation multiplies the original annular weight by the normalizing exponential described above. The bias is preserved because $`L(0)=1`$. One may use any sufficiently small scale of a fixed even bump outside a countable exceptional set. For each prescribed zero, the transform as a function of scale is entire and equals one at scale zero. Thus it can be chosen nonzero at every odd strip zero simultaneously. Denote the smoothed score by $`V`$. Its explicit formula above becomes

```math
V(x)=2\Re\sum_{\gamma\gt0} v_z e^{zx}+o(1),
\qquad z=\alpha+i\gamma,\quad
v_z=m_zL(z)(H_{w_1}(z),\ldots,H_{w_n}(z)).
```

The sum runs over distinct positive-ordinate complex points, including different real parts at a common ordinate. The combined multiplicity $`m_z`$ is positive. The negative-ordinate coefficients are conjugates. The sum of the vector coefficient norms is finite by smoothness and the zero count. The finite real-axis exclusions and functional equations remove real odd zeros. The unsigned square term tends to its bias by the prime ideal theorem, and all remaining prime-power and finite place corrections tend to zero. These limits are uniform on every fixed compact interval translated to infinity. [Explicit formula and arithmetic remainder](#what-survives-on-sparse-grids)

Let $`W`$ be the real span of the real and imaginary parts of all $`v_z`$ with $`\alpha\gt0`$. If it has positive dimension, choose finitely many such coefficient vectors spanning it and write $`a_*`$ for the smallest of their real shifts. Fix $`0\lt a_0\lt a_*`$. For unit real $`u`$, with $`P_W`$ orthogonal projection, compactness of the unit sphere of $`W`$ gives a constant $`b\gt0`$ with $`\max_j|u\cdot v_j|\ge b\|P_Wu\|`$ on this finite list. We may add finitely many growing coefficients from $`\alpha\ge a_0`$. Write

```math
M_X(u)=\max_j e^{\alpha_jX}|u\cdot v_j|,
\qquad \|P_W u\|\le b^{-1}e^{-a_*X}M_X(u)\quad(X\ge0).
```

When $`W=0`$ put $`M_X=0`$ and choose any $`0\lt a_0\lt1/2`$. Montgomery’s fixed-conductor density bound counts the complex points with $`\alpha\ge a_0`$ and positive ordinate at most $`T`$ by $`O(T^\lambda\log^A(2+T))`$ for some $`\lambda\lt1`$. In particular their reciprocal ordinates sum to a finite number. Call a positive line ordinate *clean* here if no point with $`\alpha\ge a_0`$ has that ordinate. [Montgomery's density theorem](https://www.mathnet.ru/eng/mat564) · [Fixed-conductor form](https://www.numdam.org/item/10.5802/jtnb.1162.pdf#page=11)

The clean line coefficients span all real directions. Indeed a nonzero real combination of the independent weights has a nonzero entire Mellin transform of exponential type. Its cancelled ordinates up to $`T`$ number $`O_u(T+1)`$ by Jensen. There are $`\gg T\log T`$ distinct positive simple line ordinates from the real primitive quadratic factor by Bauer. Discarding those at growing ordinates in the fixed half-strip still leaves surviving coefficients. A real direction orthogonal to all clean coefficients must therefore be zero. Select finitely many clean positive ordinates $`g_l`$ so that

```math
A=\min_{\|u\|_2=1}\ \max_l |u\cdot v_{ig_l}|\gt0 .
```

Any given finite clean list with this property will do. We prove the assertion for any prescribed $`0\lt r\lt A`$. This strict form gives the numerical applications. [Bauer's Corollary 1](https://doi.org/10.4064/aa-93-1-37-52) · [Jensen argument](#a-compact-weight-cannot-cancel-all-the-residues)

### Extracting the critical-line coefficients locally

Use the bilateral transform $`\widehat P(z)=\int e^{zt}\,dP(t)`$. For each distinct $`\alpha+i\gamma`$ with $`\alpha\ge a_0`$ and $`\gamma\gt0`$, and each $`s\in\{0,\pm g_l\}`$ over the finite list, form the probability at $`\pm h`$ with weights there proportional to $`e^{\mp\alpha h}`$, where $`h=\pi/(2|\gamma+s|)`$. This kills the moment at $`\alpha+i(\gamma+s)`$ and its conjugate. Cleanliness prevents zero denominators. The finitely many shifts preserve summability of the radii by comparison with $`1/\gamma`$ beyond twice the greatest shift. The infinite convolution exists as the distribution of an absolutely convergent sum of independent bounded random variables. Its compact support is common to all the desired cancellations. Denote it by $`P`$. This uses the classical compact-support construction for trigonometric sums with the two weights tilted to handle a real part. [Kozma and Oravecz](https://arxiv.org/abs/math/0601021) · [Two-point transform](https://github.com/leecrossley/theta-host-split-inert-race/blob/asymptotic-waiting-factors-v1.0.0/results/asymptotic-waiting-factors/PROOF.md#a-positive-average-with-zero-limiting-mean)

Convolve $`P`$ with the uniform probability on $`[-R,R]`$ to get $`\mu`$. For each $`l`$ test $`u\cdot V(X+t)`$ against the positive measure

```math
(1+\Re(\theta e^{ig_l(X+t)}))\,d\mu(t),\qquad |\theta|=1 .
```

Every coefficient with $`\alpha\ge a_0`$ is killed, even after the modulation. On the line the zero-frequency contribution from the matching conjugate pair equals $`\Re(\overline\theta\,u\cdot v_{ig_l})`$. For any nonzero resulting frequency $`k`$ the transfer has modulus at most $`|\sin(Rk)/(Rk)|`$. All unmatched line contributions therefore tend to zero as $`R\to\infty`$, uniformly in $`X,u,l,\theta`$, by absolute summability and the finite list of modulations. The total mass tends uniformly to one by the same bound at $`g_l`$. Fix $`R`$ large enough that the leakage and the mass error use less than one quarter of the strict gap $`A-r`$ after normalization. Negative real shifts give $`o(1)`$ as $`X\to\infty`$ by dominated convergence. For this fixed support, the shifts $`0\lt\alpha\lt a_0`$ cost at most a constant times $`e^{a_0X}\|P_Wu\|`$ by summability. Thus for any fixed $`M_0`$, as soon as $`M_X(u)\le M_0`$ this cost tends uniformly to zero. Choosing $`l`$ to maximize the line projection and choosing either phase gives both oriented support values greater than $`r+\delta`$ on one compact translated interval, for some $`\delta\gt0`$ independent of $`u,X`$. Positive tests transfer these witnesses to the original $`F`$ by enlarging the support by $`\epsilon`$.

### Directions with a large growing coefficient

Here $`W\ne0`$. Construct a positive smooth probability $`K`$ killing the unshifted half-strip $`\alpha\ge a_0`$, using just the unshifted two-point factors and a short smooth probability. There is a closed interval with interior where $`K\ge c\gt0`$. For each of the finitely many chosen growing exponents $`z_j`$ make a complex smooth test $`Q_j`$ supported inside this interval. Its transform vanishes at every point of the half-strip of either sign except $`z_j`$, where it is nonzero. To construct it, take the convolution of the two-point factors for ordinates above $`G`$. Their total radius tends to zero with $`G`$. Their limiting probability has nonzero moment at $`z_j`$ since $`\Re e^{z_j t}=e^{\alpha_jt}\cos(\gamma_jt)\gt0`$ on a sufficiently small support. Use a nonnegative smooth bump near zero with the same property. For the finitely many remaining distinct exponents $`z`$ other than $`z_j`$, apply to the bump the differential operator obtained by substituting $`-d/dt`$ in $`\prod_z (Z-z)/(z_j-z)`$. Integration by parts supplies all the prescribed zeros and keeps the support. Convolve with the tail and translate into the positive interval. This kills the conjugate and every other point at the same ordinate. [Compact interpolation proof](https://github.com/leecrossley/theta-host-split-inert-race/blob/asymptotic-waiting-factors-v1.0.0/results/asymptotic-waiting-factors/PROOF.md#a-signed-test-inside-that-average)

Choose a common $`\lambda\gt0`$ with $`\lambda|Q_j|\le K/2`$. Test with $`K\pm\lambda\Re(\theta Q_j)`$, densities of masses at most $`3/2`$. Select the index attaining $`M_X(u)`$ and orient the phase so its contribution to the signed average is positive. For either sign this contribution is at least $`b_1M_X(u)`$ for a fixed $`b_1\gt0`$. The other half-strip terms cancel. The line and negative shifts and the prime-only remainder cost $`O(1)`$ uniformly. The smaller positive shifts cost $`O(e^{a_0X}\|P_Wu\|)=O(M_X(u)e^{-(a_*-a_0)X})`$. All constants are independent of $`u`$. After division by the mass and increasing the fixed $`M_0`$ if necessary, both raw oriented maxima on the enlarged support exceed $`r+\delta+b_2M_X(u)`$ whenever $`M_X(u)\gt M_0`$ and $`X`$ is sufficiently large, for some $`b_2\gt0`$. We use the preceding line estimate on $`M_X(u)\le M_0`$ with this same cutoff.

Place the two fixed supports, the smoothing translations and an extra unit of slack at each end inside a single translated log interval $`(\log H,\log H+D)`$. One may set $`X=\log H+T_0`$ with fixed $`T_0`$ for all sufficiently large $`H`$. Every real unit direction has support greater than $`r+\delta`$ on one common compact subinterval of this window. Its compact image has compact convex hull in finite dimensions, so separation puts the closed radius-$`r`$ ball in the ordinary convex hull. All smaller $`H\ge1`$ can use the same compact set of starts from the first valid late window, after increasing $`C=e^D`$ by a finite factor.

Carathéodory gives an initial representation by at most $`n+1`$ vectors. If all $`n+1`$ are necessary in this representation, they are affinely independent with strictly positive weights. Slide the first start continuously towards the second. At the first failure of interior membership in the moving simplex, compactness of barycentric weights supplies a representation with a zero weight or an affinely dependent set. Either allows a vertex to be eliminated. If no earlier failure occurs, the final two vertices coincide. This proves the continuous bound $`n`$. For $`n=1`$ this is the intermediate value theorem. [Connected-set version of Carathéodory](https://www.renyi.hu/~barany/cikkek/131.pdf)

### Prime-power starts

Let $`0\lt y_1\lt y_2\lt\cdots\to\infty`$ satisfy $`y_{k+1}-y_k=O(y_k^{1-\kappa})`$ for some fixed $`\kappa\gt0`$. Assume also that the second distributional derivative of each log profile $`e^{v/2}w_j(e^v)`$ is a finite signed measure. In the same radius-$`r`$ conclusion every start can now be restricted to this mesh, with one $`C`$ for all directions and all targets. The exact upper bound becomes $`n+1`$ mesh starts. To prove uniformity even for directions that vary with the height, we keep track of the growing coefficients separately in finitely many bands.

If $`W\ne0`$, choose $`a_0\lt\min(a_*,\kappa/4)`$. Partition $`[a_0,1/2]`$ into finitely many disjoint bands of width at most $`\kappa/4`$. In each band take a finite spanning family among the real and imaginary parts of its smoothed coefficient vectors and add the corresponding complex points to the selected growing list. The span, denoted by $`W_I`$, is also the span for the original coefficients, since multiplication by a nonzero complex scalar preserves the real span of its real and imaginary parts. The projection for a band with lower endpoint $`a`$ satisfies $`\|P_{W_I}u\|\ll e^{-aX}M_X(u)`$. An empty span contributes zero. Keep the original high spanning family in the maximum so the estimate involving $`a_*`$ still holds. If $`W=0`$ take $`a_0\lt\kappa/4`$ and there are no growing contributions.

Any point $`x`$ of the witnessing compact interval has bounded distance from $`X`$, and some mesh point has logarithm $`x+h`$ with $`|h|\ll e^{-\kappa X}`$. For the **original** coefficient at each distinct $`z`$, including its combined multiplicity, the second-derivative bound is $`O(m_z/(1+\gamma^2))`$. A shift costs an extra factor $`\ll\min(1,(1+|\gamma|)|h|)`$. Summed with multiplicity this gives $`O(|h|\log^2(2/|h|))`$ by the zero count. Terms in a band $`[a,b]`$ cost at most this last quantity times $`O(e^{bX}\|P_{W_I}u\|)`$, hence at most $`O(e^{\kappa X/4}M_X(u)|h|\log^2(2/|h|))`$. The line, negative shifts and real shifts below $`a_0`$ are bounded in the same difference calculation by $`O(e^{a_0X}|h|\log^2(2/|h|))`$. The actual prime-only remainder changes by $`o(1)`$, uniformly in unit directions since there are finitely many coordinates. The function $`t\log^2(2/t)`$ is increasing for all sufficiently small positive $`t`$. In particular the total rounding cost is

```math
o(1)\big(1+M_X(u)\big)
```

uniformly in the whole unit sphere as $`X\to\infty`$. For bounded $`M_X(u)`$ the fixed slack $`\delta`$ absorbs the error. Above the fixed cutoff the extra $`b_2M_X(u)`$ absorbs it. The unit of endpoint slack ensures the rounded observations stay inside the original window. Thus the finite set of mesh points in each sufficiently late window has support greater than $`r`$ in every direction. Separation and the finite initial-height bridge prove the claim. There are only countably many choices of at most $`n`$ mesh values, and each of their convex hulls lies in an affine space of dimension at most $`n-1`$. Their union has $`n`$-dimensional measure zero. Almost every target in the ball therefore needs exactly $`n+1`$ mesh starts.

In particular take $`y_k=p_k^s`$ where $`p_k`$ is the $`k`$th ordinary rational prime and $`s\gt0`$ is fixed. Baker, Harman and Pintz give $`p_{k+1}-p_k=O(p_k^{0.525})`$. The mean value theorem gives the required mesh estimate with $`\kappa=0.475/s`$. This includes starts at the primes themselves. Every finite dyadic delay observable above has the required log regularity. For its finite positive quadrature all 710 brackets are below 133, inside Platt’s finite line verification for each of the sixteen conductors. They are therefore clean for every permissible $`a_0`$. Coalescing any shared line zeros only increases the absolute coefficient through the common positive multiplicity. For any even positive smoothing of log radius $`10^{-6}`$ or smaller the real multiplier at these ordinates is at least $`1-(133\cdot10^{-6})^2/2`$ by the cosine bound. The replay verifies that every strict lower bound in the earlier table persists with this loss. Thus $`A`$ exceeds the tabulated $`r`$ separately in all 28 cases by compactness. Choose the common smoothing outside the countable exceptional set as above. [Prime gap theorem](https://www.cs.umd.edu/~gasarch/BLOGPAPERS/BakerHarmanPintz.pdf) · [Platt's Theorem 7.1](https://arxiv.org/abs/1305.3087) · [Finite margin check](fixed_weight_checks.py)

Take the maximum of the finite factors for these 28 cases and enlarge it to include every last dyadic endpoint. It works for all heights, from $`H=1`$ onwards, and for every predictor at the chosen mesh or on the continuum. For the affine error, the last coefficient is one, so the uncentred predictions have error differences exceeding $`2r`$ between some two individual observations. Every intercept then gives at least one absolute error strictly greater than $`r`$. This proves the forecast bounds with fifteen preceding observations, and the corresponding statement for every other length in the table.
