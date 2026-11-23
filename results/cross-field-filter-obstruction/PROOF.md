# Cross-field prediction of dyadic prime-ideal scores

## The two fields and the error

Take $`(p,q)=(547,727)`$ or $`(571,599)`$. In each case set $`K=\mathbf Q(\zeta_{16})^+`$, $`F=K(\sqrt{pq})`$, $`E=K(\sqrt{-p},\sqrt{-q})`$. At a good prime ideal $`v`$ of $`F`$ write $`\epsilon_v=1`$ for splitting and $`\epsilon_v=-1`$ for inertia in $`E/F`$. The fibre product of the first-step theta trace jumps gives this character. Omit the primes above $`2pq`$. Any further fixed finite omissions give the same limiting statements. [547/727 trace lemma](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-short-shifts/RaceTrace.lean#L83) · [571/599 host and trace](https://github.com/leecrossley/theta-host-split-inert-race/blob/8fca7cad37185e3d7951f0dffff5e2442d1dfa24/HostAdmissible.lean#L80)

Use $`a=\log2`$, $`r=\sqrt2`$ and the same weight on each field,

```math
w(u)=\frac{\max(0,\min(2-2/u,4/u-1))}{2r(r-1)^2},\quad
S(y)=-y^{-1/2}\sum_v\epsilon_v\log(Nv)\,w(Nv/y),\quad
D(y)=\frac{S(2y)-S(y)}a .
```

Let $`X,Y`$ denote the target and comparison fields in either order. Fix $`T_X=S_X`$ or $`T_X=D_X`$. Choose real $`c_0,\ldots,c_{m-1}`$ with $`c_{m-1}=1`$, a real $`b`$ and a real finite signed measure $`\nu`$ of compact support on the log-height line. Write $`P(Z)=\sum_{j=0}^{m-1}c_jZ^j`$. The error, in log coordinates, is

```math
F(x)=\sum_{j=0}^{m-1}c_j T_X(2^je^x)-\int S_Y(e^{x+t})\,d\nu(t)-b .
```

This permits any fixed affine rule on the preceding $`m-1`$ target values and any compactly supported linear filter of the opposite levels. Shifted opposite moves are included by taking differences of translated measures. The logarithmic and power-sampled errors are

```math
R(X_0)^2=\frac1{X_0}\int_0^{X_0}|F(x)|^2\,dx,\qquad
R_k(N)^2=\frac1{\log N}\sum_{n\le N}\frac{|F(k\log n)|^2}{n},\qquad k\gt0 .
```

Every table entry is a strict lower bound for both upper limits, uniformly over the choices of $`c,b,\nu`$. Deleting any fixed initial segment gives the same upper limits. The proof below also permits a general regular compact weight in the comparison field. The exact distance theorem uses shifted copies of its standard level.

## Grouping the two explicit formulas

Each field has eight odd primitive factors in $`\zeta_E/\zeta_F`$. They are the quadratic characters of conductors $`p,q`$, each twisted by the four even characters of $`\mathbf Q(\zeta_{16})^+`$. For a compact annular weight $`v`$ put $`H_v(z)=\int_0^\infty v(u)u^{z-1/2}du`$. When $`e^{t/2}v(e^t)`$ has a finite measure as distributional second derivative, the strip explicit formula for the prime-only score is a constant plus a sum of $`H_v(\rho-1/2)e^{(\rho-1/2)x}`$ over the odd zeros, with remainder $`o(1)`$. The prime ideal theorem gives the square bias. The coefficients satisfy $`H_v(z)\ll_v(1+|\Im z|)^{-2}`$ uniformly on the strip. These facts apply to the tent, its translates, its difference and its convolution with any $`\nu`$ above. Bounded support and finite total variation preserve the measure bound on the second derivative. [Full strip formula and prime-only correction](https://github.com/leecrossley/theta-host-split-inert-race/blob/compact-filter-oscillation-v1.0.0/results/compact-filter-oscillation/PROOF.md#what-survives-on-sparse-grids)

Group contributions at the **same complex point** $`z=\rho-1/2`$ across the two fields. With $`m_X(z),m_Y(z)`$ their total multiplicities, put

```math
b_z=m_X(z)H_{T_X}(z)P(e^{az})-m_Y(z)H_w(z)\int e^{tz}\,d\nu(t).
```

Here $`H_{T_X}=H_w`$ for a level and $`H_w(z)(e^{az}-1)/a`$ for a move. The combined expansion is $`F(x)=b_0^*+\sum_{z\ne0}b_z e^{zx}+o(1)`$ where the real constant includes the bias and intercept. For $`\Re z\gt0`$ a nonzero coefficient gives a pole with that residue in the continued Laplace transform of $`F`$. This also follows directly from the prime-only Dirichlet series times the transforms. The unsigned square term is regular there except for its pole on the boundary. The remainder $`o(1)`$ has a Laplace transform holomorphic on the open right half-plane. Conjugate points have conjugate coefficients. [Prime-only residues](https://github.com/leecrossley/theta-host-split-inert-race/blob/joint-levels-and-moves-v1.0.0/results/joint-levels-and-moves/PROOF.md#residues-for-a-real-projection)

Let $`d=\sup(\{0\}\cup\{\Re z:b_z\ne0\})`$. If $`d=0`$, the left-half-plane terms tend to zero in sum by dominated convergence. The other terms form an absolutely and uniformly convergent trigonometric series. Hence both mean squares have actual limits equal to

```math
|b_0^*|^2+2\sum_{t\gt0}|b_{it}|^2 .
```

For continuous time this follows by averaging a finite sum, then taking uniform limits. For fixed powers, distinct frequencies differ by $`\lambda\ne0`$ and the harmonic mean of $`n^{ik\lambda}`$ tends to zero. Indeed the difference between $`\sum_{n\le N}n^{-1+ik\lambda}`$ and its integral is $`O_{k,\lambda}(1)`$ by integrating the absolute derivative on unit intervals. The integral is bounded and the total harmonic mass is asymptotic to $`\log N`$. This weighting is essential to this calculation. [Akbary, Ng and Shahabi, limiting distributions of prime error terms](https://arxiv.org/abs/1306.1657)

Suppose $`d\gt0`$. Both RMS quantities have upper exponential growth rate $`d`$ on the log-height scale, so each has infinite upper limit. The strip formula gives the upper bound $`|F(x)|\ll e^{dx}+o(1)`$. To prove the reverse inequality for continuous time, suppose the upper rate is smaller. Choose an exponent $`\omega`$ between it and $`d`$, with $`\omega\gt0`$. The squared integral up to $`X_0`$ is then $`O((1+X_0)e^{2\omega X_0})`$. Cauchy–Schwarz on successive unit intervals makes the Laplace integral of $`F`$ absolutely and locally uniformly convergent for $`\Re z\gt\omega`$. This contradicts a surviving pole there.

For power samples use the left-endpoint step function $`J(x)=F(k\log n)`$ on $`[k\log n,k\log(n+1))`$. If the sampled upper rate were less than $`d`$, choose $`\omega`$ above that rate and $`\max(0,d-1/k)`$, still below $`d`$. Because the step lengths are bounded by $`k/n`$, the same mean-square argument gives a holomorphic Laplace transform for $`J`$ on $`\Re z\gt\omega`$. On a large step the log displacement is $`O_k(e^{-x/k})`$. Each grouped zero term changes by at most a fixed constant times $`|b_z|e^{dx}\min(1,(1+|\Im z|)|h|)`$ for displacement $`h`$ with $`|h|\le1/2`$. Counting zeros with multiplicity gives $`O(T\log(2+T))`$ in the two products. Splitting at height $`1/|h|`$ gives

```math
F(x)-J(x)=O_k\!\left(e^{(d-1/k)x}(1+x)^2\right)+o(1).
```

Thus the difference has a holomorphic transform on $`\Re z\gt\omega`$ too, again contradicting a surviving pole. The sampled upper rate is taken as $`\limsup_{N\to\infty}\log(1+R_k(N))/(k\log N)`$ and the continuous rate as $`\limsup_{X_0\to\infty}\log(1+R(X_0))/X_0`$. This proof allows an off-line zero to cancel between the fields before taking the supremum. [Earlier sparse-grid increment estimate](https://github.com/leecrossley/theta-host-split-inert-race/blob/compact-filter-oscillation-v1.0.0/results/compact-filter-oscillation/PROOF.md#what-survives-on-sparse-grids)

## Frequencies that survive every comparison filter

At any positive $`t`$ for which $`m_X(it)\ge1`$ and $`m_Y(it)=0`$, the opposite-field contribution vanishes for every $`\nu`$. The absolute target coefficient is at least $`A(t)|P(e^{iat})|`$, where

```math
G(t)=\frac{3-2r\cos(at)}{2(r-1)^2(1/4+t^2)},\qquad
A(t)=G(t)\ \text{for levels},\qquad
A(t)=G(t)\frac{2|\sin(at/2)|}{a}\ \text{for moves}.
```

The change of centre in the Mellin kernel has absolute value one on the line. For any finite set $`I`$ of distinct such frequencies the squared upper RMS, in continuous time and on every fixed power, is therefore at least $`2\sum_{t\in I}A(t)^2|P(e^{iat})|^2`$. In the bounded case this drops nonnegative terms of the mean square. In the other case the upper RMS is infinite. This is the unconditional inequality.

The 381 positive brackets in 547/727 and 329 in 571/599 come from the compact-filter certificate. They are disjoint within each field across all eight factors. Hardy signs at both endpoints give a critical-line zero in each bracket. One calculation encloses every opposite Dirichlet value on the whole closed interval and excludes zero. A second checks the midpoint magnitude through Hardy's function, whose absolute value on the line equals that of the Dirichlet value. For a primitive nonprincipal character of conductor $`q`$, periodic partial sums have magnitude at most $`q`$. Differentiating the partial-summation integral on $`\Re s=1/2`$ gives $`|L'(s)|\le q(2+4|s|)`$. The upper endpoint bounds $`|s|`$ by $`1/2+t_{\max}`$. The midpoint magnitude exceeds this derivative bound times the radius in all 5,680 cases. No multiplicity, global disjointness or independence of ordinates is assumed. [Closed interval and midpoint calculations](exclusions.py) · [Hardy sign replay](check_zeros.py) · [Bracket integers](evidence/brackets.tsv)

## From 710 brackets to RMS lower bounds

Define $`r_\ell=2\sum_{t\in I}A(t)^2\cos(\ell at)`$ and the real symmetric Toeplitz matrix $`M_{ij}=r_{|i-j|}`$, with indices from zero. Then $`c^TMc=2\sum_{t\in I}A(t)^2|P(e^{iat})|^2`$. For each listed dimension the matrix is positive definite. Positive semidefiniteness follows from this identity and strict positivity also follows from the rational lower matrix below. If $`e`$ is the last coordinate vector, Cauchy–Schwarz in the $`M`$-inner product gives

```math
1=(e^Tc)^2\le(e^TM^{-1}e)(c^TMc).
```

Equality holds for the normalized multiple of $`M^{-1}e`$. The first calculation encloses the true moments directly on all closed intervals and solves $`Mu=e`$ with rigorous preconditioning. It checks $`u_{m-1}\gt0`$ and $`h^2u_{m-1}\lt1`$ for every stated radius $`h`$. Thus the minimum of this finite energy is strictly greater than $`h^2`$. It bounds every comparison filter and intercept at once. [Interval systems](radii.py) · [Algebraic implication](CrossFieldPrediction.lean)

The second calculation evaluates the equivalent positive half-angle form of $`G`$, builds each harmonic by multiplying midpoint phases and adds $`\ell a`$ times the bracket radius to the cosine uncertainty. Round each moment to a specified dyadic rational $`v_\ell`$ and bound its error by $`\delta_\ell`$. Subtract an upward-rounded $`\delta_0+2\sum_{\ell=1}^{127}\delta_\ell`$ from the diagonal of the rational Toeplitz matrix. The resulting rational matrix $`B`$ satisfies $`M\ge B`$ as quadratic forms in every listed dimension, by the row-sum bound for symmetric errors. The choice of rounded centre is arbitrary because every rounding error is enclosed.

For completeness, let $`V_j`$ be the least energy for predicting one entry from the preceding $`j`$ under the rational moments $`s_\ell`$ of $`B`$. Start at $`V_0=s_0`$ with an empty list of predictors. If $`u_1,\ldots,u_j`$ are the optimal coefficients, the forward residual and the reversed backward residual each have squared norm $`V_j`$ by the symmetric Toeplitz identity. Both are orthogonal to the intervening $`j`$ coordinates. Their inner product is $`p=s_{j+1}-\sum_{\ell=1}^ju_\ell s_{j+1-\ell}`$, since projecting out the intervening coordinates leaves the same inner product. Subtracting $`\kappa=p/V_j`$ times the backward residual gives the next residual. Its squared norm is $`V_{j+1}=V_j(1-\kappa^2)`$ and its predictors are $`u_\ell-\kappa u_{j+1-\ell}`$ followed by $`\kappa`$. This is the usual Toeplitz prediction recursion. Equivalently these are successive Schur complements, so verified positive pivots inductively establish positive definiteness, starting with the positive scalar $`s_0`$. All 128 pivots are checked as exact rationals in each of the four systems. The relevant $`V_{m-1}`$ exceed the stated squared thresholds. [Exact rational calculation](prediction_recursion.py) · [Output](evidence/check.txt)

## The exact distance when the two spectra are on the line

For this section only, suppose all sixteen odd factors have their strip zeros on the critical line. This is GRH for these factors. Fix the target polynomial $`P`$ and observable. Write $`\alpha_t`$ for its signed complex coefficient at positive frequency $`t`$ and $`\beta_t`$ for the coefficient of the comparison level there. For continuous time and every harmonically sampled fixed power, the infimum of the limiting **mean square** over $`b,\nu`$ equals

```math
2\sum_{\substack{t\gt0\\\beta_t=0}}|\alpha_t|^2 .
```

The lower bound follows by orthogonality. Choose the intercept to remove the constant. Every shared frequency can be recovered in the closure of compact *past* filters of the other field. Indeed average $`S_Y(e^{x+s})e^{-its}`$ over $`-T\le s\le0`$ with mass $`1/T`$. Each other frequency has a nonzero difference from $`t`$ and its integral tends to zero. Absolute summability permits termwise passage to the limit uniformly in the trigonometric part, leaving $`\beta_t e^{itx}`$. To produce a real conjugate pair use twice the real part after multiplying the kernel by $`\alpha_t/\beta_t`$. For a fixed averaging length the prime-only remainder tends to zero as $`x\to\infty`$. Taking a finite set of shared pairs first and then approximating all shared pairs in mean square proves the reverse bound. Each approximating filter has compact support in the past and finite total variation. Thus past-only and two-sided comparisons have the same infimum. Finite-precision numerical bounds in the table need no critical-line hypothesis.

The horizon of the target history matters. Under the same hypothesis the target has an almost-periodic part with summable coefficients. For any prescribed accuracy, a sufficiently long single dyadic delay approximates its level or move on every late start by simultaneous approximation of finitely many phases and a uniform bound for the remaining tail. This is compatible with positive error bounds at each of the displayed fixed history lengths.
