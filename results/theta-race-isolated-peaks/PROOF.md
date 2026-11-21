# The exact upper limit, local shape and isolation of near-maximum scores

## Statement

Take $`K=\mathbf Q(\zeta_{16})^+`$, $`F=K(\sqrt{397669})`$ and $`E=K(\sqrt{-547},\sqrt{-727})`$. At a good prime ideal $`v`$ of $`F`$, set $`\epsilon_v=1`$ for split and $`-1`$ for inert in $`E/F`$. The first-step theta trace gives these signs. Omit all places above $`2,547,727,343057911569,343070636977`$. Put $`a=\log2`$, $`r_0=\sqrt2`$, $`c=2(r_0-1)^2`$, $`B=r_0c`$ and fix the score width at two. [Field and trace](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-short-shifts/PROOF.md) · [Trace identity in Lean](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-short-shifts/RaceTrace.lean#L83)

```math
W(u)=\max(0,\min(2-2/u,4/u-1)),\qquad
S(y)=\frac1{B\sqrt y}\sum_v-\epsilon_v\log(Nv)\,W(Nv/y),\qquad
s=\log(2y),\quad Q(s)=S(e^s/2).
```

Assume the generalised Riemann hypothesis (GRH) for the primitive Dirichlet functions in the 547/727 race and its base field, including the Riemann zeta function. The upper limit of $`Q(s)`$ as $`s\to\infty`$ is $`M=1+V(0)`$, given explicitly below, with $`175.71160384847\lt M\lt 175.71160384849`$. If $`s,t\ge300`$ and $`1/2\le|s-t|\le10^9`$, the following conclusions hold simultaneously for every such real pair.

1. At most one of $`Q(s),Q(t)`$ can reach 162.
2. If $`Q(s)\gt 170`$ then $`Q(s)-Q(t)\gt 20`$. If $`Q(s)\gt 175`$ then $`Q(s)-Q(t)\gt 40`$.

Thus a large score gives both gains into that centre and losses out of it across an entire continuum of allowed dilations. For centres $`s\ge10^9+300`$ both flanks are available throughout the range. The upper uniform log density of starts with $`S(y)\ge162`$ is at most $`1/(2\cdot10^9)`$. The lower uniform log density of starts with $`S(y)\gt 175`$ is at least $`2^{-48007}`$. Here the lower uniform density is $`\liminf_{T\to\infty}\inf_{A\ge0}T^{-1}\int_A^{A+T}\mathbf1_{\mathcal E}(u)\,du`$ in log height. The upper density uses $`\limsup`$ and $`\sup`$. Shifting all log coordinates by $`a`$ makes no difference to these densities.

Alpöge's announced Dirichlet GRH theorem supplies this input. [Announcement](https://alpo.ge/remark.pdf) [Independent formal checks](https://math.mit.edu/~drew/grh/) Without GRH, the two interval covers still show that the terms $`2G(\nu)(1-\cos(\nu h))`$ from the listed signed-bracket pairs sum to more than 56 for every $`1/2\le h\le10^9`$. The Lean lemmas and the compact-test formula at the origin are unconditional. For each nontrivial zero $`\rho`$ set $`t_\rho=(\rho-1/2)/i`$. The formula evaluates $`\sum_\rho G(t_\rho)`$ with multiplicity. Under GRH each $`t_\rho`$ is real and the sum equals $`V(0)`$ for the profile below.

## The unshifted prime profile

Let $`\nu`$ range, with multiplicities, over the ordinates of nontrivial zeros of the eight odd primitive factors of $`\zeta_E/\zeta_F`$. Their conductors are $`547,4376,8752,8752,727,5816,11632,11632`$. Write

```math
G(t)=\frac{3-2r_0\cos(at)}{c(t^2+1/4)},\qquad
V(s)=\sum_\nu G(\nu)e^{i\nu s},\qquad D(s)=V(0)-V(s).
```

The weights are strictly positive. Under GRH, which Alpöge's announced Dirichlet GRH theorem supplies, the zeros lie on the critical line, including those of the base-field factors in the inversion below. Conjugation pairs the ordinates across the family, so $`V`$ is real and even. The absolute convergence is uniform by the quadratic denominator and the zero count. In particular $`D(s)\ge0`$. [Upstream theorem checks](https://math.mit.edu/~drew/grh/) [Factors and zero weights](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-distant-interval/PROOF.md)

We detail how the original prime-only score recovers this **unshifted** profile. The bilateral transform of $`Q`$ for the full prime sum is $`\mathcal G(z)L_D(1/2+z)`$, where

```math
\mathcal G(z)=\frac{2r_0\cosh(az)-3}{c(z^2-1/4)},\qquad
L_D(w)=\frac{L'}L(w)-\frac{\zeta_F'}{\zeta_F}(2w)+R(w),
\qquad L=\zeta_E/\zeta_F .
```

The singularities at $`z=\pm1/2`$ in $`\mathcal G`$ are removable. The inverse kernel is $`g(v)=4\sinh(\max(0,a-|v|)/2)/B`$, a nonnegative function of mass 2. Here $`\mathcal G(it)=G(t)`$. Euler factors give the displayed comparison because the square of the quadratic sign is one and the extension is unramified at finite places. The correction $`R(w)=\sum c_n n^{-w}`$ starts with powers of order three. The pole of $`\zeta_F`$ at one contributes $`\mathcal G(0)/2=1`$ to the score. These inputs and the contour inversion with extra compact boxes are recorded in the predecessor. The same inversion applies with $`\mathcal G`$ itself. It has the same quadratic decay on the integration strip, and the original score is continuous, allowing the box widths to tend to zero. [Prime Euler factors and transforms](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-waiting-bound/PROOF.md) · [Compact-kernel inversion and remainders](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-fixed-factor/PROOF.md)

For completeness, on $`|\Re z|\le1/2`$ one has $`|\mathcal G(z)|\le100/(1+|\Im z|^2)`$. For $`|\Im z|\le1`$ use $`2e^{a/2}`$ from its positive mass. Otherwise bound the numerator by $`2r_0\cosh(a/2)+3`$ and the product in the denominator by $`|\Im z|^2`$. The sum of $`(1+|\Im\rho|^2)^{-1}`$ over the sixteen factors is below 1,000. Base zeros have real part $`-1/4`$ after the square correction and residue factor $`1/2`$. The sixteen trivial-zero strings begin at shifted real part at most $`-1/2`$, spaced by at least one. Their transforms have modulus at most $`2e^{a|\Re z|}`$. Finally $`\sum |c_n|n^{-7/20}\lt 1600`$. A participating $`n`$ is at most $`e^{s+a}`$, and the single-window normalization is $`r_0 e^{-s/2}/B`$. Every excluded individual prime norm is below $`e^{216}`$, hence has left the support at $`s\ge300`$. These bounds give

```math
|Q(s)-1-V(s)|\le 2\cdot10^6e^{-s/4}
   +128\frac{e^{-(s-2)/2}}{1-e^{-(s-2)}}
   +22000e^{-3s/20}\lt 10^{-12}\quad(s\ge300).
```

The right side decreases to zero. This proves a uniform error on every higher tail. [Numerical constants](peak_profile.py) · [Zero and residual coefficient bounds](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-fixed-factor/PROOF.md)

Put $`A_*=\sum_\chi\log(q_\chi/\pi)-8(\gamma_E+a)`$. The even compact moment gives the exact value

```math
V(0)=A_*g(0)-\frac8c[2r_0T(a)-3a],\qquad
T(a)=\frac{(1+z)\log(1+z)+(1-z)\log(1-z)}{2\sqrt z},\quad z=\tfrac12.
```

It follows that $`Q(s)\le M+o(1)`$. For any finite set of its frequencies there are arbitrarily large positive times making all phases arbitrarily close to one. For example, in the closure of the integer orbit on its finite torus, compactness gives differences of arbitrarily distant orbit indices with arbitrarily small phase error. Uniform absolute convergence handles the omitted frequencies. Thus $`V(s)`$ approaches $`V(0)`$ along a sequence tending to infinity, proving the exact upper limit. [Compact gamma moment](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-distant-interval/PROOF.md) · [Directed enclosure](peak_profile.py)

## A frequency test for every real displacement in a billion-unit band

For $`1/2\le h\le10^9`$ we prove **$`D(h)\gt 56`$**. Paired ordinates contribute $`2G(\nu)(1-\cos(\nu h))\ge0`$. It suffices to lower-bound any subset of the signed-bracket pairs. Both calculations use disjoint true-ordinate intervals from the earlier Hardy-sign certificate. Each interval algorithm prices a whole closed cell and subdivides when its first bound is too small. No set of unexamined gaps lies between sampled points. [Signed brackets and Hardy verification](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-factor-400/PROOF.md)

The first test uses 100 retained positive ordinates chosen by their kernel weights. Its grid unit is $`1/128`$. For each $`t`$ a certified numerator $`p`$ approximates $`t/(128\cdot2\pi)`$ by $`p/2^{64}`$, with accumulated phase error less than $`1/10000`$ turn for every possible centre tick $`m\le128(10^9+1)`$. For a cell centred at $`m/128`$ with radius $`d`$, take the top eight bits $`k`$ of $`mp\bmod2^{64}`$. Throughout the cell the distance of $`th/(2\pi)`$ to an integer is at least

```math
u=\max\left(0,\frac{\min(k,255-k)}{256}-\frac{td}{2\pi}-\frac1{10000}\right),
```

using a lower endpoint for the difference. The price in hundredths is the nonnegative integer lower bound obtained from $`200G(t)(1-\cos(2\pi u))`$. The cosine deficit increases from zero on the distance interval $`[0,1/2]`$. A sum strictly greater than 5,600 pays for the whole cell. Start with $`[1/2,1]`$ and the consecutive unit cells $`[n,n+1]`$ up to $`10^9`$. Bisection keeps every centre integral in ticks through the allowed depth. [Binary phase bounds](cover.py) · [Dyadic interval cover](scan_cells.c)

The second test uses the 128 *lowest* retained positive ordinates, a three-quarter grid starting at $`1/2`$, a tick of $`1/256`$ and modulus $`2^{61}-1`$. It encloses the frequencies anew and builds $`G`$ from its half-angle formula. Its per-centre phase error is below $`1/250`$ turn. With $`r`$ the smaller distance of the residue to zero modulo this modulus, $`k=\lfloor r/2^{53}\rfloor`$ gives lower distance $`k/256`$ since $`2^{61}-1\lt 2^{61}`$. Subtract the cell radius times the frequency in turns and the phase error, then clip at zero. In place of a cosine evaluation, the second test uses

```math
1-\cos u\ \ge\ \frac{u^2}{2}-\frac{u^4}{24}+\frac{u^6}{720}-\frac{u^8}{40320},
\qquad 0\le u\le\pi .
```

The remaining alternating terms decrease in magnitude. The last staggered cell covers the right endpoint, including its small overshoot. The verifier bounds phase errors through that overshoot. Unsigned 128-bit products are reduced by folding at bit 61. They are below $`2^{122}`$, so one fold and at most one subtraction suffice. Each table is generated from whole zero brackets with directed arithmetic at replay time. [Independent phase calculations](second_cover.py) · [Staggered interval cover](scan_staggered.c) · [Both complete transcripts](evidence/check.txt)

## Why any sufficiently high score is isolated

In the Hilbert space with one coordinate per zero, put $`u(s)_\nu=\sqrt{G(\nu)}e^{i\nu s}`$. Translation is unitary and $`\|u(s)-u(t)\|^2=2D(s-t)`$. Thus $`\sqrt{D(s-t)}\le\sqrt{D(s)}+\sqrt{D(t)}`$. If two prime scores reach 162 after log centre 300, then each individual deficit is less than 13.713. Their separation consequently has deficit below $`4\cdot13.713\lt 56`$. This excludes the whole band in either direction. Within any interval of length $`10^9`$ above 300 the high-score set has diameter at most $`1/2`$ and measure at most $`1/2`$. Tiling increasingly long intervals gives its upper uniform density. Initial bounded ranges contribute zero in the limit. [Rational implications in Lean, taking the geometric inequality as a hypothesis](IsolatedPeaks.lean)

For a score exceeding 175 its individual deficit is smaller than 0.713. For a score exceeding 170 the corresponding bound is 5.713. Let $`a_1=\sqrt{D(s-t)}`$ and $`b_1=\sqrt{D(s)}`$. The displacement certificate gives $`a_1^2\gt 56`$. In both cases $`a_1\gt b_1`$, and the triangle inequality gives

```math
V(s)-V(t)=D(t)-D(s)\ \ge\ (a_1-b_1)^2-b_1^2
       =a_1^2-2a_1b_1.
```

In the first case the lower bound exceeds $`40+2\cdot10^{-12}`$ and in the second it exceeds $`20+2\cdot10^{-12}`$. The two prime-only errors preserve the strict inequalities. Further, Cauchy–Schwarz gives

```math
\sup_h|V(s+h)-V(h)|\ \le\ \sqrt{2V(0)D(s)}.
```

Indeed factor $`e^{i\nu s}-1`$ out of the translated sum and use the two weighted square norms. If $`s_n\to\infty`$ and $`Q(s_n)\to M`$, then $`D(s_n)\to0`$. Hence $`Q(s_n+h)\to1+V(h)`$ uniformly on each fixed bounded interval of real $`h`$. The exact peak height fixes the entire local limiting shape.

## Positive density of the highest scores

For the 6,000 selected positive ordinates write $`A`$ for the total weight of both partners and $`R=V(0)-A`$. The closed formula gives $`0\lt R\lt 0.647`$ and the independent digamma integral gives $`R\lt 0.649`$. Set $`q=256`$ and define $`w(s)`$ as the product over the selected positive ordinates of $`\max(1-q\|s\nu/(2\pi)\|,0)`$, where the inner norm denotes distance to the nearest integer. On its support where the product is positive, the selected contribution is at least $`A\cos(2\pi/q)`$. [Independent even-mass calculations](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-silent-powers/PROOF.md) · [Tent Fourier coefficients and uniform means](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-paired-density/PROOF.md)

Each tent is an interval autocorrelation. Its Fourier coefficients are nonnegative and absolutely summable, with constant coefficient $`1/q`$. In products, any extra coincidences of frequencies contribute nonnegatively. Thus the uniform interval mean of $`w`$ is at least $`q^{-6000}`$. The omitted series also has nonnegative Fourier coefficients. The product of this series with $`w`$ therefore has nonnegative uniform mean. These mean assertions follow term by term, since a nonconstant exponential has vanishing long-window average uniformly in the lower endpoint and the coefficient sums converge absolutely.

On $`s\ge300`$, the score is below 176 and $`0\le w\le1`$, so $`\mathbf1_{\{Q\gt 175\}}\ge w(Q-175)`$. After discarding an initial bounded range and taking uniform long-window means, the resulting lower bound is at least

```math
[A\cos(2\pi/256)-174-10^{-12}]\cdot256^{-6000}
       \gt 2^{-48007}.
```

The bracket exceeds $`1/128`$ even with the independent 0.649 omitted-mass estimate. This establishes a positive lower density of peaks to which the forty-unit conclusion applies. Their exact density and the shortest separation between distinct high-score clusters remain open. [Exact budget](peak_profile.py) · [Verification output](evidence/check.txt)
