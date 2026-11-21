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
