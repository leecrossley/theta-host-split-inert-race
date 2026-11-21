# A distant interval in the 547/727 prime race

## Theorem

Let $`K=\mathbf Q(\zeta_{16})^+`$, $`F=K(\sqrt{397669})`$ and $`E=K(\sqrt{-547},\sqrt{-727})`$. For a prime ideal $`v`$ of $`F`$, put $`\epsilon_v=1`$ when it splits in $`E/F`$ and $`-1`$ when it is inert. Omit the primes above $`2,547,727,343057911569,343070636977`$. The product of the two normalized first-step theta trace jumps gives this sign at each good place. [Fibre trace identity](https://github.com/leecrossley/theta-race-short-shifts/blob/483c31f2df2132ed2165620ea55e20e969228bfb/RaceTrace.lean#L83)

Write $`a=\log2`$, $`r_0=\sqrt2`$, $`c=2(r_0-1)^2`$ and $`B=r_0c`$. With the race parameters $`b=r=2`$, the weighted prime-ideal sum and its normalized move are

```math
W(u)=\max(0,\min(2-2/u,4/u-1)),\qquad
S(y)=\frac{1}{B\sqrt y}\sum_v-\epsilon_v\log(Nv)\,W(Nv/y),
\qquad f(x)=\frac{S(2y)-S(y)}a,\quad x=\log(2r_0y).
```

Take $`\theta=1.9735`$, $`\ell=3895686.23`$, $`h=3895691.51`$, $`L=h-\ell=132/25`$, $`C=2e^L`$ and $`H_* = e^\ell/(2r_0)`$. All terminating decimals in the proof denote exact rational numbers.

1. At $`H=H_*`$, every real $`y\gt H`$ with $`2y\le CH`$ satisfies $`S(2y)-S(y)\lt\theta a`$.
2. The same conclusion holds for arbitrarily large $`H`$. At other arbitrarily large $`H`$, every such $`y`$ satisfies $`S(2y)-S(y)\gt-\theta a`$.

Thus any factor guaranteeing either strict move after every height must exceed $`C`$. The certificate locates one failure and leaves the first failure and the optimal universal factor open.

## Spectral profile and the eight factors

The eight primitive characters of $`\zeta_E/\zeta_F`$ are $`\chi_p\xi^j`$, $`p\in\{547,727\}`$, $`0\le j\lt4`$. Here $`\chi_p=(\cdot/p)`$ and $`\xi`$ is the even quartic character of conductor 16 with $`\xi(5)=i`$. Their primitive conductors at either prime are $`p,8p,16p,16p`$. The verifier checks the Conrey labels against these definitions on every unit class modulo $`16p`$, as well as primitivity and parity. [Character and zero verifier](check_zeros.py)

The primitive Dirichlet critical-line theorem applies to these characters and to the base-field characters (including the primitive character of modulus one). Sutherland's independent comparison identifies the upstream theorem with the stated primitive Dirichlet problem. The release checks digests of its statement, implementation, axiom report and independent kernel records. This supplies the global zero-location result used here. [Pinned Lean theorem](https://github.com/leanprover-community/riemann-hypothesis/blob/371a1a26d5e17f604428e3a8ce6ce341cf517434/src/Main.lean#L12) · [Independent verification](https://math.mit.edu/~drew/grh/)

For the compact even test $`g`$ and real frequency $`t`$ put

```math
g(s)=\frac4B\sinh\!\left(\frac{\max(0,a-|s|)}2\right),\quad
G(t)=\int_{\mathbb R}g(s)e^{-its}\,ds
=\frac{3-2r_0\cos(at)}{c(t^2+1/4)},\qquad
\alpha(t)=\frac{2G(t)\sin(at/2)}a .
```

The transform of the translated difference $`j(s)=[g(s+a/2)-g(s-a/2)]/a`$ at $`it`$ is $`J(it)=i\alpha(t)`$ with the Fourier convention above. The prime-only limiting profile is $`U(x)=\sum_\nu i\alpha(\nu)e^{i\nu x}`$, where the ordinates of the eight odd factors enter with multiplicity. The coefficient bound $`O((1+t^2)^{-1})`$ and the zero count give absolute uniform convergence. Conjugation pairs $`\nu`$ and $`-\nu`$ across the family, even for a nonreal character. Thus $`U`$ is real and odd. A positive ordinate and its partner contribute $`-2\alpha(\nu)\sin(\nu x)`$. [Prime-only explicit formula](https://github.com/leecrossley/theta-race-fixed-factor/blob/12b7c94cfb9586c25fb29e7c7aaf0a20adaf5b2f/PROOF.md#global-input-and-transform)

## A positive bound for all omitted ordinates

The inequality $`2|\sin v|\le 1+\sin^2v`$ and positivity of $`G(t)`$ give

```math
|\alpha(t)|\le m(t)=\frac{G(t)(3-\cos(at))}{2a},\qquad
M=\sum_\nu m(\nu)=\frac{3V(0)-V(a)}{2a},\qquad
V(x)=\sum_\nu G(\nu)e^{i\nu x}.
```

The identities use the same symmetric multiset of ordinates. The even compact test makes both values of $`V`$ accessible at small prime norms. Write $`A=\sum_\chi\log(q_\chi/\pi)+8\psi(3/4)`$, with $`\psi=\Gamma'/\Gamma`$. Let $`P_g(x)`$ be the sum over prime powers $`n=p^k`$ with summand $`-(\sum_\chi\chi(p)^k)(\log p)g(x-\log n)/\sqrt n`$. The symmetrized explicit formula for these eight completed odd functions is

```math
\begin{aligned}
V(x)&=P_g(x)+P_g(-x)+A g(x)-8I_g(x),\\[2pt]
I_g(x)&=\int_0^\infty \frac{e^{-t/2}}{2\sinh t}
  [g(x-t)+g(x+t)-2g(x)]\,dt .
\end{aligned}
```

Indeed the test on the dual prime term is reflected, giving the plus sign for the even kernel. The completed gamma factors all have parameter $`3/4`$ at the centre. Their symmetric difference gives the displayed integral by the integral formula for the digamma function. For $`x=0,a`$ the support permits only $`2,3,4`$. At 2 only the two odd-conductor characters contribute and their values are opposite. At 3 the quartic sum vanishes. The power 4, when present, lies at $`g(a)=0`$. Hence both prime terms vanish at both points. [Explicit formula and gamma-kernel derivation](https://github.com/leecrossley/theta-race-profile-interval/blob/264a0976afd3a59266400cf8040dbe700692dda8/PROOF.md#sum-the-gamma-correction)

For completeness put $`u_k=2k+3/2`$, $`z=e^{-|s|}`$ and

```math
T(s)=\sum_{k\ge0}\frac{e^{-u_k|s|}}{u_k^2-1/4}
=\frac{(1+z)\log(1+z)+(1-z)\log(1-z)}{2\sqrt z},\qquad T(0)=a.
```

Expand $`e^{-t/2}/(2\sinh t)=\sum_{k\ge0}e^{-u_kt}`$. The integrated absolute second difference for each term is at most $`2\mathrm{Lip}(g)/u_k^2`$, justifying interchange. For $`u=u_k`$, convolution of $`g`$ with $`e^{-u|\cdot|}`$ gives

```math
\frac{2u\,g(x)+[r_0(e^{-u|x+a|}+e^{-u|x-a|})-3e^{-u|x|}]/c}{u^2-1/4}.
```

This follows by multiplying their Fourier transforms and taking partial fractions. Subtract $`2g(x)/u`$ and sum. The remaining constant series is $`\sum_k[u_k(u_k^2-1/4)]^{-1}=\pi-4a`$, by partial fractions and the digamma reflection and duplication identities. The corresponding coefficient of $`g(x)`$ becomes $`A_*=\sum_\chi\log(q_\chi/\pi)-8(\gamma_E+a)`$. At the two required abscissae this gives

```math
\begin{aligned}
V(0)&=A_*g(0)-\frac8c[2r_0T(a)-3a],\\[2pt]
V(a)&=-\frac8c[r_0(T(2a)+a)-3T(a)],\\[2pt]
386.7375649676&\lt M\lt386.7375649677.
\end{aligned}
```

The table gives 5,200 disjoint-within-character positive rational brackets. Its four integer columns are the conductor, Conrey label, and lower and upper endpoints in units of $`10^{-11}`$. For each bracket the Hardy function has strictly opposite signs at its endpoints. This real continuous normalization has the same critical-line zeros as its Dirichlet function. Choose one zero per bracket and its negative partner in the aggregate family. Overlap of ordinates between different factors counts correctly with multiplicity. Any extra zero in a bracket and every unlisted zero continue to contribute to $`M`$. The selected pairs remove exactly $`2\sum m(\nu)`$ from the majorant. The remainder satisfies

```math
\sup_x\left|U(x)+2\sum_{\rm selected}\alpha(\nu)\sin(\nu x)\right|
\le M-2\sum_{\rm selected}m(\nu)\lt1.72.
```

The computation uses an upper bound on $`M`$ and lower interval bounds for the mass of the selected pairs. Completeness of a zero list and simplicity of its members are unnecessary. [Rational brackets](evidence/zero-brackets.tsv) · [Hardy function conventions](https://flintlib.org/doc/acb_dirichlet/l-functions.html) · [Spectral verifier](spectral_bound.py)

## Two whole-interval calculations

For each selected bracket $`I`$, let $`t`$ be its midpoint and $`d`$ its radius. Define the negative midpoint head $`Q(x)=2\sum\alpha(t)\sin(tx)`$. Replacing the actual selected zeros by their midpoints costs less than $`0.03`$ over the full interval. A bound for this replacement is the sum, over brackets, of

```math
2\sup_{v\in I}|\alpha(v)-\alpha(t)|+2|\alpha(t)|\,h\,d.
```

It follows from the Lipschitz constant 1 for the sine. At 4,097 equally spaced nodes, including both endpoints, directed balls give $`Q(x)\gt0.02`$. The second derivative is bounded by $`2\sum|\alpha(t)|t^2\lt320000`$. Between consecutive nodes of spacing $`\Delta=L/4096`$ the linear interpolant differs from $`Q`$ by at most $`320000\Delta^2/8\lt0.067`$. This is the interpolation remainder bound from Rolle's theorem, since the product of the distances to the two endpoints is at most $`\Delta^2/4`$. Consequently $`U(x)\lt1.72+0.03+0.067-0.02\lt1.8`$ throughout $`[\ell,h]`$. [First verifier](spectral_bound.py)

The second calculation starts from the real integral $`I_g`$. It uses 8,192 directed rectangle enclosures on each dyadic piece from $`a/2^{16}`$ to $`32a`$. If $`D=2\cosh(a/2)/B`$, the integral below the first piece costs at most $`Da/2^{16}`$ because $`g`$ is $`D`$-Lipschitz and the integration kernel is at most $`1/(2t)`$. Above $`b_0=32a`$ it costs at most $`4g(0)e^{-3b_0/2}/[(3/2)(1-e^{-2b_0})]`$. Using $`A`$ and the integral directly gives $`386.71\lt M\lt386.76`$ and an omitted mass less than $`1.74`$. For the head it reconstructs $`G(t)`$ through $`3-2r_0\cos(at)=(r_0-1)^2+4r_0\sin^2(at/2)`$. It encloses the true frequency and full coordinate cell simultaneously in every sine, bisecting cells until the upper bound clears $`1.9`$. The 1,179 successful closed cells cover the interval, including their endpoints. Both calculations use Arb ball arithmetic and the same sign-certified table. [Second verifier](interval_check.py) · [Replay output](evidence/check.txt)

## From the profile to prime ideals

For $`x\ge1000`$ the prime-only score has the following pointwise error bound from compact-kernel inversion.

```math
|f(x)-U(x)| \le 2\cdot10^6e^{-x/4}
 +128\,\frac{e^{-(x-2)/2}}{1-e^{-(x-2)}}+22000e^{-3x/20}\lt10^{-20}.
```

For the base zeros, the shifted real part is $`-1/4`$, the residue coefficient is $`1/2`$, the transform bound is $`1000/(1+t^2)`$ and the sum of $`(1+(\Im\rho)^2)^{-1}`$ is below 1,000. Replacing $`t`$ by $`\Im\rho/2`$ costs at most a factor of four. The pole contributes zero since $`J(0)=0`$. The sixteen trivial-zero strings start at shifted real part at most $`-1/2`$ and have spacing at least one. The support of $`j`$ has radius $`3a/2\lt2`$ and its integral of absolute value is at most $`4/a\lt6`$, giving the stated geometric bound. The Euler-product remainder after the prime-square correction starts with third powers and is bounded pointwise by $`22000e^{-3x/20}`$. The omitted prime norms are below $`e^{53}`$ and have left the support. Each term in the displayed upper bound is decreasing. The verifier checks their sum at 1,000. [Bounds for this exact score](https://github.com/leecrossley/theta-race-fixed-factor/blob/12b7c94cfb9586c25fb29e7c7aaf0a20adaf5b2f/PROOF.md#fixed-log-window)

At $`H_*`$ all permitted starts have $`\ell\lt x\le h`$, since $`2y\le2e^LH_*`$ is equivalent to $`y\le e^LH_*`$. The prime-only error therefore proves the first assertion at this explicit height.

For every $`\eta\gt0`$, there are arbitrarily large positive translations $`\tau`$ with $`\sup_x |U(x+\tau)-U(x)|\lt\eta`$. Truncate the absolutely convergent series with uniformly small tail. Compactness of the torus gives two arbitrarily widely separated integer times with all retained phases arbitrarily close, and their difference is a uniform approximate period. This argument allows any rational relations among the ordinates.

Take the translation error below $`1/40`$ and choose translations with the full translated window above 1,000. The total error in transferring the profile to $`f`$ is below $`1/20`$. On translates of $`[\ell,h]`$, this excludes the gain because $`1.9+1/20\lt\theta`$. Oddness gives $`-U(x)\lt1.9`$ on $`[-h,-\ell]`$, so sufficiently late translates of that interval exclude the loss. The respective starting heights are $`e^{\tau+\ell}/(2r_0)`$ and $`e^{\tau-h}/(2r_0)`$. The right endpoint is included in both cases. The five local Lean deductions check the budget arithmetic, explicit gain exclusion, reflection, translation implication and log-centre conversion. [Lean proofs](DistantInterval.lean)
