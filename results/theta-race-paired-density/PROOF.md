# Synchronized one-way windows at positive lower logarithmic density

## Statement

Let $`K=\mathbf Q(\zeta_{16})^+`$, $`F=K(\sqrt{397669})`$ and $`E=K(\sqrt{-547},\sqrt{-727})`$. For an unramified prime ideal $`v`$ of $`F`$, put $`\epsilon_v=1`$ for splitting in $`E/F`$ and $`\epsilon_v=-1`$ for inertia. The normalized first-step theta traces give these signs. Omit the places above $`2,547,727,343057911569,343070636977`$. Throughout the scale and move parameters are $`b=r=2`$. [Field and trace identity](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-short-shifts/RaceTrace.lean#L83)

Put $`a=\log2`$, $`r_0=\sqrt2`$, $`c=2(r_0-1)^2`$ and $`B=r_0c`$. Define the score and normalized change by

```math
W(u)=\max(0,\min(2-2/u,4/u-1)),\quad
S(y)=\frac1{B\sqrt y}\sum_v-\epsilon_v\log(Nv)\,W(Nv/y),\quad
f(x)=\frac{S(2y)-S(y)}a,\quad x=\log(2r_0y).
```

Set $`\ell_+=469417229.340`$, $`\ell_-=484947092.225`$, $`L_+=6`$ and $`L_-=45/8`$. Write $`H_\pm=e^{\ell_\pm}/(2r_0)`$. Let $`\mathcal E`$ be the set of $`t\ge0`$ satisfying **both** conditions

```math
\min_{0\le z\le L_+} f(\ell_++t+z)\gt 1/40,\qquad
\max_{0\le z\le L_-} f(\ell_-+t+z)\lt -1/40.
```

Thus every doubling starting in the first window gains more than $`a/40`$ and every doubling starting in the companion window loses more than $`a/40`$. Their respective completion deadlines are $`2e^6H_+e^t`$ and $`2e^{45/8}H_-e^t`$. The companion starting height always equals the first multiplied by $`\exp(15529862.885)`$. The set is measurable since the score is continuous and extrema over a fixed compact interval vary continuously under translation.

**Theorem.** Assume the generalised Riemann hypothesis (GRH) for the primitive Dirichlet functions in the 547/727 race and its base field, including the Riemann zeta function. With $`|\cdot|`$ denoting Lebesgue measure in log height,

```math
\liminf_{T\to\infty}\ \inf_{A\ge0}
   \frac{|\mathcal E\cap[A,A+T]|}{T}\ \ge\ 2^{-52862},\qquad
\frac{|\mathcal E\cap[0,T]|}{T}\ \ge\ 2^{-52863}\quad(T\gt 0).
```

For the first height $`H=H_+e^t`$, integration against $`dt`$ is integration against $`dH/H`$. The first bound is uniform in the lower endpoint of long log intervals. The second applies to every finite prefix measured from the specified starting height. The prime-score bounds use the critical-line input in the predecessor theorem. Computing the actual proportion and an effective convergence rate are separate questions.

## A uniform phase budget for both windows

Use the same 6,000 chosen positive ordinates $`\nu_j`$ as in the one-way window proof, each lying in a distinct bracket within its primitive factor. Pair each with its negative in the aggregate of eight factors. Extra roots inside or outside the brackets and all possible multiplicities belong to the remaining multiset. For real $`v`$ put

```math
G(v)=\frac{3-2r_0\cos(av)}{c(v^2+1/4)},\qquad
\alpha(v)=2G(v)\sin(av/2)/a,\qquad P(x)=-2\sum_{j=1}^{6000}\alpha(\nu_j)\sin(\nu_j x).
```

Three bounds from the predecessor will be used. With sign $`s=+1`$ on the first base interval and $`s=-1`$ on the second, the selected true-frequency head there satisfies $`sP(x)\gt 1.64-0.043-0.024=1.573`$. Indeed its signed midpoint head exceeds 1.64 at every interpolation node, the interpolation cost is below 0.024 and the bracket error below 0.043. The full absolutely convergent profile $`U`$ satisfies $`|U-P|\lt 1.44`$ **everywhere**. Finally $`|f(x)-U(x)|\lt 10^{-20}`$ for $`x\ge1000`$. Only the first of these bounds is restricted to the base windows. [Signed head, cubic tail and prime-only estimate](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-factor-400/PROOF.md) · [Pinned replay output](evidence/paired-window-check.txt)

Let $`w_j=2|\alpha(\nu_j)|`$. Row $`j`$ of the phase table gives positive integers $`W_j,b_j`$ for which $`w_j\lt W_j2^{-32}`$ and $`q_j=2^{b_j}\ge2`$. Put $`d_j=q_j^{-1}`$. If $`\|\nu_jt/(2\pi)\|_{\mathbf R/\mathbf Z}\le d_j`$ for all $`j`$, periodicity and the unit Lipschitz constant of sine give, for every real $`x`$,

```math
|P(x+t)-P(x)|\le \sum_j2\pi w_jd_j
 \lt \frac{44}{7}\sum_j W_j2^{-32-b_j}
 =\frac{788089118617009}{7881299347898368}\lt \frac1{10}.
```

For a signed base point this leaves more than $`1.573-0.1-1.44-10^{-20}\gt 1/40`$ after translation by $`t\ge0`$. The omitted ordinates are charged just once, at the translated point. Thus every nonnegative phase shift in the indicated box belongs to $`\mathcal E`$. The argument uses the true zero inside each bracket. It imposes no phase constraint on its midpoint in place of that zero. [Partition](evidence/phase-bins.tsv) [Rational checks](check_budget.py) [Lean arithmetic steps](PairedDensity.lean)

Each $`W_j`$ is computed by outward interval evaluation of the entire bracket at 160 bits. A second 192-bit calculation verifies every ceiling using $`3-2r_0\cos(av)=(r_0-1)^2+4r_0\sin^2(av/2)`$. To allocate bits, begin with one per ordinate and halve the largest remaining cost until the bound is strictly below $`1/10`$. The possible marginal reductions of a single coordinate decrease geometrically. Selecting the largest available reduction at each step therefore minimizes the remaining sum for a given bit count. Equivalently, the smallest selected reduction is at least the largest unselected one. Undoing the smallest selected reduction raises the exact rational bound to $`788180326378929/7881299347898368\ge1/10`$. The total $`\sum b_j=52862`$ is minimal for these rounded weights, at least one bit per ordinate, dyadic partitions, linear phase estimate and specified budget. The theorem claims no optimal density. [Partition construction](phase_budget.py)

## A finite bound from collisions

Write $`M=\prod_jq_j=2^{52862}`$. Partition each coordinate circle into $`q_j`$ half-open arcs of equal length. For $`s\in[0,T]`$, take the box of the vector $`(\nu_js/(2\pi))_j`$ modulo one. If two times land in the same box, their difference in either order obeys all phase conditions. Let $`p_k`$ be the measure of the times in box $`k`$. The measure of ordered pairs landing in one box is $`\sum_kp_k^2\ge T^2/M`$ by Cauchy–Schwarz. Let $`\mathcal B`$ denote all real lags satisfying the phase conditions. This set is even, and $`\mathcal B\cap[0,\infty)\subseteq\mathcal E`$. Integrating on the time difference gives

```math
\frac{T^2}{M}\le \sum_kp_k^2
\le \int_{-T}^T (T-|u|)\mathbf 1_{\mathcal B}(u)\,du
\le 2T\,|\mathcal E\cap[0,T]|.
```

This proves the asserted finite bound even if some of the phase coordinates coincide.

## A uniform bound from nonnegative Fourier coefficients

On the circle of circumference one let $`k_d(v)=\max(1-\|v\|/d,0)`$ for $`0\lt d\le1/2`$. The tent is the periodic convolution of the indicator of $`[-d/2,d/2]`$ with its reflection, divided by $`d`$. Its Fourier coefficients are $`d`$ at zero and $`d(\sin(\pi n d)/(\pi n d))^2`$ at nonzero integers $`n`$. They are nonnegative and absolutely summable with sum $`k_d(0)=1`$. The product

```math
K(t)=\prod_{j=1}^{6000} k_{d_j}\big(\nu_jt/(2\pi)\big)
=\sum_{\mathbf n\in\mathbf Z^{6000}}c_{\mathbf n}
 \exp\!\left(it\sum_j n_j\nu_j\right)
```

satisfies $`0\le K(t)\le\mathbf 1_{\mathcal E}(t)`$ for $`t\ge0`$. Its coefficients are products of the one-dimensional coefficients. They sum absolutely to one. The coefficient of the zero multi-index equals $`\prod_j d_j=1/M`$. Other multi-indices whose frequency $`\omega=\sum_jn_j\nu_j`$ vanishes add nonnegative terms to the same constant. For $`\omega\ne0`$ the normalized integral of $`e^{it\omega}`$ on $`[A,A+T]`$ has absolute value at most $`\min(1,2/(T|\omega|))`$, uniformly in $`A`$. Dominated convergence for the summable coefficients shows that the averages of $`K`$ converge uniformly to the sum of its zero-frequency coefficients. This sum is at least $`1/M`$. Integrating the pointwise inequality proves the uniform lower bound. The proof allows arbitrary rational relations among the ordinates and requires no estimate for how close a nonzero frequency can come to zero.
