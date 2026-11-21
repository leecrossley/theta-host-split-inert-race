# The interval and the waiting factor

## Statement

Set $`K=\mathbf Q(\zeta_{16})^+`$, $`F=K(\sqrt{397669})`$ and $`E=K(\sqrt{-547},\sqrt{-727})`$. For a prime ideal $`v`$ of $`F`$ write $`\epsilon_v=1`$ when it splits in $`E/F`$ and $`-1`$ when it is inert. Omit the places above $`2,547,727,343057911569,343070636977`$. The product of the normalized first-step theta trace jumps above each good place recovers $`\epsilon_v`$ for this field. [Concrete fibre theorem](https://github.com/leecrossley/theta-race-short-shifts/blob/483c31f2df2132ed2165620ea55e20e969228bfb/RaceTrace.lean#L83)

Let $`a=\log 2`$, $`r_0=\sqrt2`$, $`c=2(r_0-1)^2`$ and $`B=r_0c`$. For the dyadic race parameters $`b=r=2`$ put

```math
W(u)=\max(0,\min(2-2/u,4/u-1)),\qquad
S(y)=\frac{1}{B\sqrt y}\sum_v-\epsilon_v\log(Nv)\,W(Nv/y).
```

Write $`\theta=1.9735`$, $`\ell=1.008796429610`$, $`h=3.221800332553`$, $`L=h-\ell=2.213003902943`$, and $`C=2e^L`$. All terminating decimals in the proof denote exact rationals.

Assume the generalised Riemann hypothesis (GRH) for the primitive Dirichlet $L$-functions in the 547/727 race and its base field, including the Riemann zeta function. Alpöge's announced Lean theorem supplies this zero-location input. [Announcement](https://alpo.ge/remark.pdf) [Independent formal checks](https://math.mit.edu/~drew/grh/)

For arbitrarily large $`H`$, **every** $`y\gt H`$ with $`2y\le CH`$ satisfies $`S(2y)-S(y)\gt-\theta a`$. For arbitrarily large $`H`$, every such $`y`$ satisfies $`S(2y)-S(y)\lt\theta a`$. In particular no factor at most $`C`$ guarantees both required moves after every height.

## The prime signal and its limiting profile

Use $`x=\log(2r_0 y)`$ and $`f(x)=(S(2y)-S(y))/a`$. In the ratio $`\zeta_E/\zeta_F`$ the eight odd primitive factors have characters $`\chi_p\xi^j`$, with $`p=547,727`$, $`0\le j\lt4`$, $`\chi_p=(\cdot/p)`$, and $`\xi`$ the even quartic character of conductor 16 with $`\xi(5)=i`$. Their conductors for each $`p`$ are $`p,8p,16p,16p`$. The zero theorem for primitive Dirichlet functions puts all their nontrivial zeros on the critical line. It also covers the factors of the base field needed to pass from logarithmic derivatives to prime ideals. Sutherland’s independent checks matched the named upstream Lean statement and reported the three standard axioms. [Primitive theorem](https://github.com/leanprover-community/riemann-hypothesis/blob/371a1a26d5e17f604428e3a8ce6ce341cf517434/src/Main.lean#L12) · [Verification records](https://math.mit.edu/~drew/grh/)

The dyadic kernel and its transform use the following normalization.

```math
g(s)=\frac4B\sinh\!\left(\frac{\max(0,a-|s|)}2\right),\quad
j(s)=\frac{g(s+a/2)-g(s-a/2)}a,\quad
J(z)=\frac{2\sinh(az/2)}a\,\frac{2r_0\cosh(az)-3}{c(z^2-1/4)}.
```

The support of $`j`$ is $`[-j_0,j_0]`$ with $`j_0=3a/2`$. Write $`\nu`$ for the ordinates of the eight odd factors with multiplicity. Then

```math
U(x)=\sum_\nu J(i\nu)e^{i\nu x},\qquad f(x)-U(x)\longrightarrow0\quad(x\longrightarrow+\infty).
```

Indeed $`J(i\nu)=O((1+\nu^2)^{-1})`$, which is summable against the zero count. The base-field square terms have shifted real part $`-1/4`$, the base pole is killed by $`J(0)=0`$, the higher-power remainder is $`O(e^{-x/10})`$, and the finitely many excluded norms eventually leave the support. These estimates and the compact-kernel inversion are proved for this precise field and score in the earlier paper. The series is absolutely and uniformly convergent. The eight characters are stable under conjugation and $`j`$ is real and odd, so $`U`$ is real and odd. [Prime-only limit and explicit formula](https://github.com/leecrossley/theta-race-factor-18/blob/31ea5d495ca29f7f17869d1a3f321000266e4867/PROOF.md#the-limiting-oscillation)

## Sum the gamma correction

For an abscissa $`x`$ put

```math
P(x)=\sum_{n=p^k}\frac{c_n\log p}{\sqrt n}j(x-\log n),
\qquad c_{p^k}=-\sum_\chi\chi(p)^k .
```

The explicit formula for the eight completed odd functions gives, with $`\psi=\Gamma'/\Gamma`$ and conductors $`q_\chi`$,

```math
\begin{aligned}
U(x)&=P(x)-P(-x)+A j(x)-8 I(x),\\[2pt]
A&=\sum_\chi\log(q_\chi/\pi)+8\psi(3/4),\\[2pt]
I(x)&=\int_0^\infty\frac{e^{-t/2}}{2\sinh t}
       \,[j(x-t)+j(x+t)-2j(x)]\,dt .
\end{aligned}
```

The integral converges at zero by the Lipschitz bound on $`j`$ and at infinity by exponential decay. The kernel in it is $`\sum_{k\ge0}e^{-u_k t}`$, where $`u_k=2k+3/2`$. It can be integrated term by term since the integrated absolute second difference is bounded by $`28/u_k^2`$ for each term. [Derivation for this family](https://github.com/leecrossley/theta-race-factor-18/blob/31ea5d495ca29f7f17869d1a3f321000266e4867/PROOF.md#evaluate-the-zero-sum-near-the-origin)

Define $`D_u(v)=r_0(e^{-u|v+a|}+e^{-u|v-a|})-3e^{-u|v|}`$. Convolution of $`g`$ with $`e^{-u|\cdot|}`$ gives $`(2u g+D_u/c)/(u^2-1/4)`$. This follows by multiplying their bilateral transforms and taking partial fractions. Take the same difference of two translates as in $`j`$, then subtract $`2j(x)/u`$. The integral for one exponential becomes

```math
\frac{j(x)/(2u)+[D_u(x+a/2)-D_u(x-a/2)]/(ac)}{u^2-1/4}.
```

Both sums in this expression have elementary evaluations here. First,

```math
\sum_{k\ge0}\frac1{u_k(u_k^2-1/4)}
 =2\psi(3/4)-\psi(1/2)-\psi(1)=\pi-4a .
```

The first equality follows by writing the summand as $`1/(k+1/2)+1/(k+1)-2/(k+3/4)`$. The second follows from the digamma reflection and duplication formulas. Next put $`z=e^{-|s|}`$ and define

```math
T(s)=\sum_{k\ge0}\frac{e^{-u_k|s|}}{u_k^2-1/4}
 =\frac{(1+z)\log(1+z)+(1-z)\log(1-z)}{2\sqrt z},
 \qquad T(0)=\log2 .
```

In fact $`u_k^2-1/4=(2k+1)(2k+2)`$. Subtract the series with denominators $`2k+2`$ from the one with denominators $`2k+1`$, and use continuity at $`z=1`$. The series also shows that $`T`$ decreases with $`|s|`$. The complete formula is

```math
\begin{aligned}
U(x)&=P(x)-P(-x)+A_*j(x)\\[2pt]
 &\quad-\frac8{ac}\Big\{
    r_0[T(x+3a/2)-T(x-3a/2)]
   -(3+r_0)[T(x+a/2)-T(x-a/2)]\Big\},\\[2pt]
A_*&=\sum_\chi\log(q_\chi/\pi)-8(\gamma_E+a).
\end{aligned}
```

This expression includes all of the gamma contribution at the point where the support crosses the origin. It contains four evaluations of the same function $`T`$.

## The interval certificate

On $`[\ell-10^{-12},h+10^{-12}]`$, one has $`e^{-x+j_0}\lt2`$ and $`e^{x+j_0}\lt71`$. Thus $`P(-x)=0`$. For each power $`p^k\le70`$ with nonzero coefficient the pair $`(p^k,c_{p^k}\log p)`$ is one of

```math
(4,-2\log2),\ (16,-2\log2),\ (31,8\log31),\quad
(47,-8\log47),\ (49,-8\log7),\ (64,-2\log2).
```

For odd $`p`$ the quartic character sum vanishes unless $`p^k\equiv\pm1\bmod16`$, when it equals 4. At 2 only the two odd-conductor characters contribute. This gives the six entries by Euler’s criterion.

The first verifier applies the four-term formula at 192-bit precision to a closed interval with rational endpoints. At a kernel breakpoint it uses $`T(0)=\log2`$ and monotonicity in absolute distance. If an interval enclosure fails to clear the required margin, it bisects that exact rational interval and checks both closed children. Its finite successful run proves

```math
U(x)+\theta\gt10^{-11}\qquad(\ell\le x\le h).
```

At the neighbouring points the same calculation gives the following outward bounds. [Closed-form verifier](closed_profile.py) · [Replay output](evidence/check.txt)

| Log coordinate | Enclosure for $`10^{12}(U(x)+\theta)`$ |
|---|---:|
| 1.008796429609 | $`(-313,-312)`$ |
| 1.008796429610 | $`(81,82)`$ |
| 3.221800332553 | $`(40,41)`$ |
| 3.221800332554 | $`(-65,-64)`$ |

### Independent calculation from the real integral

The second verifier reconstructs the prime powers independently and evaluates the original tent at $`y=e^x/(2r_0)`$ and $`2y`$. It certifies the same margin over three adjoining pieces. [Second verifier](integral_check.py)

**Left of the breakpoint.** On $`[\ell,1.035]`$ write each smooth branch of $`j(s)`$ as $`c_i e^{s/2}+d_i e^{-s/2}`$. The breakpoints in $`t\ge0`$ are $`j_0-x,\ x-a/2,\ x+a/2,\ x+j_0`$, in that order. Set $`w=e^{-t/2}`$. If $`j(x-t)+j(x+t)=A_1/w+B_1w`$ on a piece and $`j_x=j(x)`$, its integral runs from the smaller to the larger $`w`$ with integrand $`(2A_1w+2B_1w^3-4j_xw^2)/(1-w^4)`$. An antiderivative is displayed below. On the piece touching $`w=1`$, $`A_1=B_1=j_x`$. There the antiderivative reduces to $`2j_x(\arctan w-\log(1+w))`$, which is regular at 1. The verifier applies these expressions to whole cells by interval arithmetic, including a separate evaluation at the left exterior point.

```math
\frac{A_1-B_1}{2}\log(1+w^2)
 +\left(j_x-\frac{A_1+B_1}{2}\right)\log(1-w)
 -\left(j_x+\frac{A_1+B_1}{2}\right)\log(1+w)
 +2j_x\arctan w .
```

**Across the breakpoint.** On $`[1.035,1.040]`$, ten closed cells cross the support breakpoint. Here $`|j|\le9`$ and $`j`$ has Lipschitz constant at most 14. For the integral above, $`0\lt t\lt0.001`$ costs at most $`14/1000`$ before multiplication by 8. The tail $`t\gt8`$ costs at most $`36e^{-12}/[(3/2)(1-e^{-16})]`$. Directed rectangle sums enclose the intervening integral. The stated constants follow by differentiating the two branches of $`g`$.

**Right of the breakpoint.** For $`x\ge1.040\gt j_0`$ the contour expansion reduces to $`U(x)=P(x)-8\sum_{k\ge0}J(-u_k)e^{-u_k x}`$. Every term added to $`P`$ is positive. The first 16 terms suffice for the lower bound throughout $`[1.040,h]`$. For the upper bound at the right exterior point, the entire omitted tail is smaller than the expression below, with $`d=x-j_0`$ and $`u=33.5`$. Expand the hyperbolic functions in $`J`$, bound each numerator by its largest exponential and sum a geometric progression. The verifier checks that this bound is less than $`10^{-25}`$ at that point.

```math
\frac{8r_0e^{-u d}}{ac(u^2-1/4)(1-e^{-2d})}.
```

## The two consequences

The set $`\{x:U(x)\gt-\theta\}`$ is open. Its component containing $`[\ell,h]`$ is an interval $`(\alpha,\beta)`$. The strict signs at $`\ell-10^{-12}`$ and $`h+10^{-12}`$ put $`\alpha`$ strictly between the two left points and $`\beta`$ strictly between the two right points, by continuity. Hence

```math
2.213003902943\lt\beta-\alpha\lt2.213003902945 .
```

This bounds further extensions of this component. Other components and other limiting profiles can give different bounds.

For every $`\eta\gt0`$ the function $`U`$ has arbitrarily large positive translations $`\tau`$ with $`\sup_t|U(t+\tau)-U(t)|\lt\eta`$. Indeed first truncate the absolutely convergent series with a uniformly small tail. For the finite remaining phase vector at integer times, compactness of the torus gives two arbitrarily widely separated times with phases arbitrarily close. Subtracting these times gives the required uniform translations. The argument accommodates rational relations among the ordinates.

Take $`\eta=10^{-11}/3`$. Choose such translations large enough that $`|f-U|\lt\eta`$ on each translated interval. The certified margin implies $`f(\tau+t)\gt-\theta`$ for every $`t\in[\ell,h]`$. Set $`H=e^{\tau+\ell}/(2r_0)`$. Every $`y\gt H`$ with $`2y\le2e^LH`$ then has its log centre in $`[\tau+\ell,\tau+h]`$. This proves the no-loss assertion, including the right endpoint. Oddness gives $`U(t)\lt\theta-10^{-11}`$ on $`[-h,-\ell]`$ and the same argument with $`H=e^{\tau-h}/(2r_0)`$ proves the no-gain assertion. Lean checks only the final no-loss step, as `recurring_no_loss`, with the certified margin and the close late translates as hypotheses. [Lean transfer step](ProfileInterval.lean)

The resulting factor lies strictly between $`18.286280578160`$ and $`18.286280578161`$. It improves the earlier factor $`2e^{11/5}`$. The proof gives arbitrarily late failures and leaves the first numerical failure undetermined. The analytic deductions are written proofs, including the explicit formula, the limit $`f-U\to0`$, the recurrence of $`U`$ and the no-gain case. The release independently replays the two interval computations and three elementary Lean steps: transfer from a close late translate, admissible starts and component length.
