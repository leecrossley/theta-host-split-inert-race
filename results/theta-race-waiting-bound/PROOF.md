# An explicit waiting interval for the 547/727 race

## Statement and prime convention

Let $`K=\mathbf Q(\zeta_{16})^+`$, $`F=K(\sqrt{397669})`$ and $`E=K(\sqrt{-547},\sqrt{-727})`$. For a prime ideal $`v`$ of $`F`$ let $`\epsilon_v=1`$ for split and $`-1`$ for inert in $`E/F`$. Use the primes outside $`2,547,727,343057911569,343070636977`$. With $`a=\log 2`$ and $`B=2\sqrt2(\sqrt2-1)^2`$ put

```math
W(u)=\begin{cases}2(1-1/u)&1 < u\le2,\\4/u-1&2 < u<4,\\0&\text{otherwise},\end{cases}
\qquad
S(y)=\frac{\sum_v-\epsilon_v\log(Nv)\thinspace W(Nv/y)}{B\sqrt y}.
```

For every real $`H\ge1`$ there are $`y_+,y_-\gt H`$ such that

```math
2\max(y_+,y_-)< e^{162000}H^{3001/2999},\qquad
S(2y_+)-S(y_+)>1.9735a,\qquad S(2y_-)-S(y_-)<-1.9735a.
\tag{1}
```

For $`\log H\ge500000000`$, $`H^{1.001}`$ also bounds both endpoints. For $`1\le H\le10^6`$ a direct certificate gives $`2\max(y_+,y_-)\le18H`$, improving to $`(6851/657)H`$ for $`100\le H\le10^6`$. For $`0\lt H\lt 1`$ the two moves starting at 5 and 9 suffice. The statements concern changes in the lead and allow either ordering of the two moves.

The prior field construction proves the degrees 8 and 16, the unramified quadratic extension at finite places, and the following identity for an admissible datum and every input parameter family. The product of the normalized first-step theta trace jumps in the fibre above a good $`v`$ equals $`\epsilon_v`$. The two large exclusions are the type primes of the level datum. Indeed membership in its exceptional set is obtained by restriction of a place over one of the two selected primes in the quadratic subfield. The other three exclusions fix the convention at rational ramification. Every omitted norm is at most $`343070636977^8\lt e^{216}`$. We can therefore use the full prime sum in the analytic proof wherever the smallest norm in its windows exceeds $`e^{216}`$. [Field and admissibility](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-short-shifts/SecondField.lean) · [Fibre trace theorem](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-short-shifts/RaceTrace.lean#L83) · [Type primes above their rational primes](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Formalism/LocalToGlobal/UnitaryNorms/CharHost.lean#L44)

## The transform and the prime powers

Write $`\chi_p=(\cdot/p)`$ and let $`\xi`$ be the even quartic character of conductor 16 with $`\xi(5)=i`$. At primitive conductors the eight factors of $`L=\zeta_E/\zeta_F`$ are $`\chi_p\xi^j`$, $`p=547,727`$, $`0\le j\le3`$. The eight factors of $`\zeta_F`$ are $`\xi^j`$ and $`\chi_{547}\chi_{727}\xi^j`$, including ordinary zeta for the trivial character. These are the characters odd and even, respectively, on complex conjugation $`(-1,-1,1)`$ in $`C_2\times C_2\times C_4`$. Abelian factorization includes the primitive Euler factors at rational ramification.

In this section use the full score, put $`Q(x)=S(e^x/2)`$ and $`f(x)=[Q(x+a/2)-Q(x-a/2)]/a`$. Then $`f(x)=[S(2y)-S(y)]/a`$ for $`y=e^x/(2\sqrt2)`$. The bilateral Laplace transform of $`f`$ on $`\Re z\gt 1/2`$ is $`J(z)D(1/2+z)`$, where

```math
D(s)=-\sum_v\epsilon_v\log(Nv)(Nv)^{-s},\quad
G(z)=\frac{2\sqrt2\cosh(az)-3}{2(\sqrt2-1)^2(z^2-1/4)},\quad
J(z)=G(z)\frac{2\sinh(az/2)}a.
```

The singularities of $`G`$ at $`\pm1/2`$ are removable. Integrate the two branches of $`W`$ after setting $`u=2Nv e^{-x}`$ to obtain $`G`$. It is the transform of a nonnegative compact kernel on $`[-a,a]`$ with mass $`G(0)=2`$. Thus $`J`$ is the transform of a signed kernel of total variation at most $`4/a\lt 6`$ and support in $`[-3a/2,3a/2]\subset(-2,2)`$. In particular $`J(0)=0`$. For $`\vert \Re z\vert \le1/2`$ we shall use

```math
|J(z)|\le \frac{1000}{1+(\Im z)^2}.
\tag{2}
```

For $`\vert t\vert =\vert \Im z\vert \le1`$ use the total variation bound. For $`\vert t\vert \ge1`$ use $`\vert (z-1/2)(z+1/2)\vert \ge t^2`$, $`\vert \cosh(az)\vert \le\cosh(a/2)`$ and $`\vert \sinh(az/2)\vert \le\cosh(a/4)`$.

Coefficient comparison over prime ideals gives

```math
D(s)=\frac{L'}{L}(s)-\frac{\zeta_F'}{\zeta_F}(2s)+D_R(s).
\tag{3}
```

The prime squares cancel since $`\epsilon_v^2=1`$. The remaining series $`D_R=\sum_n c_n n^{-s}`$ has powers at least three. At $`\sigma=2/5`$ its total absolute mass is at most

```math
\frac{2}{1-2^{-2/5}}\sum_v\frac{\log Nv}{(Nv)^{6/5}}<10000.
```

Indeed $`\sum_{v\mid p}\log Nv\le8\log p`$, and $`\sum_{n\ge2}\log(n)n^{-6/5}\lt 110`$ by $`\log n\le10n^{1/10}`$ and integral comparison. The factor 2 covers the two remaining prime-power series. In physical coordinates the part $`f_R(x)`$ uses $`n\le e^{x+2}`$. Since $`0\le W\le1`$, the two square-root normalizations and $`1/a`$ cost less than $`40e^{-x/2}`$. Consequently

```math
|f_R(x)|<40e^{-x/2}\thinspace10000\thinspace e^{(2/5)(x+2)}<10^6e^{-x/10}.
\tag{4}
```

## An isolated zero and its coefficient

Set $`d=1/1000`$. There is a unique simple zero of $`L`$ at $`1/2+i\gamma`$ with

```math
3287178298086/10^{12}<\gamma<3287178298087/10^{12}.
```

On the critical line the open $`d`$-neighborhood of $`\gamma`$ contains this zero alone, and the $`d`$-neighborhood of zero contains none. The conductor-547 factor has opposite Hardy signs at the endpoints. Throughout the larger interval around $`\gamma`$ the real part of its complex derivative exceeds 1. Thus the imaginary part of $`L(s,\chi_{547})`$ is strictly increasing with the ordinate there, and its zero is simple. Closed interval enclosures exclude zeros of every other odd factor and check the central intervals. Platt's finite-height theorem places every nontrivial zero of these eight factors through height 8500 on the critical line, since $`11632\cdot8500\lt 10^8`$. The family is closed under conjugation, giving the same coverage at negative heights. [Interval verification](verify.py) · [Independent crossing](independent.py) · [Platt, Theorem 7.1](https://arxiv.org/pdf/1305.3087#page=14)

At the selected zero, $`J(i\gamma)=iV`$ with

```math
V=\frac{3-2\sqrt2\cos(a\gamma)}{2(\sqrt2-1)^2(1/4+\gamma^2)}
       \frac{2\sin(a\gamma/2)}a,\qquad 3.342843 < V<3.342844.
```

The square term in (3) requires a **complex disk** test. For $`\omega=0,\gamma`$ and for any nontrivial zero $`\rho`$ of any factor of $`\zeta_F`$ we have

```math
\left|(\rho-1)/2-i\omega\right|\ge d .
\tag{5}
```

Here is the certificate. Use disks of radius $`r=2001/10^6`$ about $`1`$ and $`1+i\thinspace6574356596173/10^{12}`$. They contain both closed disks of radius $`2d`$ at the desired centres. For a nonprincipal character write $`C(u)=\sum_{n\le u}\chi(n)`$. A complete period supplies $`\vert C\vert \le M`$ by taking the maximum of the squared norms of Gaussian-integer partial sums. With $`N=2^{18}`$, partial summation gives errors at each centre $`s`$ of at most $`M(1+\vert s\vert )/N`$ for $`L(s)`$ and $`M[\log N+1+\vert s\vert (\log N+1)]/N`$ for $`L'(s)`$. The truncated sums and their derivatives are evaluated with complex balls.

For the Taylor remainder on the whole disk, $`\Re s\ge.99`$, $`\vert s\vert \lt 8`$ and $`\vert C(u)\vert \le\min(u,q)`$. Put $`I_k=\int_1^\infty\min(u,q) u^{-1-\Re s}(\log u)^k\thinspace du`$. Splitting at $`q`$ and using $`q\le e^{16}`$ gives

```math
I_k\le e^{.16}\left[\frac{16^{k+1}}{k+1}+
             \int_0^\infty e^{-.99v}(16+v)^k\thinspace dv\right].
```

Twice differentiating $`s\int_1^\infty C(u)u^{-s-1}\thinspace du`$ gives $`\vert L''\vert \le2I_1+8I_2\lt 30000`$. Subtract the two truncation errors with their proper weights, $`r\vert L'\vert`$ and the quadratic remainder $`15000r^2`$, as implemented in the verifier. Every resulting lower bound exceeds $`1/100`$. For the zeta factor use $`\zeta(s)=s/(s-1)-s\int_1^\infty\lbrace u\rbrace u^{-s-1}\thinspace du`$. At 1 the regularized function $`(s-1)\zeta(s)`$ stays away from zero by $`1-r-(1+r)r/.99\gt 0.9`$. At the other centre interval evaluation and the integral's derivative give $`\vert \zeta(1+2i\gamma_0)\vert -12r\gt 1/100`$, where $`2\gamma_0=6574356596173/10^{12}`$. This proves (5). The simple pole of zeta gives a term at $`z=0`$ annihilated by $`J(0)`$. [Prime-square disk verifier](square_factors.py) · [Numerical receipt](evidence/check.txt)

## Total mass of the omitted zeros

For all sixteen primitive factors together, including multiplicities, we use

```math
\sum_{\rho}\frac1{1+(\Im\rho)^2}<1000 .
\tag{6}
```

For one factor, $`\Re(2-\rho)^{-1}\ge 1/[4(1+(\Im\rho)^2)]`$. For a primitive nonprincipal Dirichlet character of parity $`k`$ the sum of these real parts is

```math
\tfrac12\log(q/\pi)+\tfrac12\psi((2+k)/2)+\Re(L'/L)(2,\chi).
```

This follows by subtracting the Hadamard logarithmic derivatives of the completed function at $`2`$ and $`-1`$ and using the functional equation. The zeros are paired by $`\rho\mapsto1-\overline\rho`$ for a fixed character. The genus-one subtraction makes the difference absolutely convergent. For zeta, completing with $`s(s-1)`$ adds $`1/2+1`$. Here $`q\le6362704`$, $`\psi(1),\psi(3/2)\lt 1`$, and the Euler product gives $`\vert (L'/L)(2)\vert \le\sum_{n\ge2}\log(n)n^{-2}\lt 2`$. Each sum of real parts is less than 13, proving (6) since $`16\cdot4\cdot13\lt 1000`$.

## A compact average at each starting height

Put $`h=\log H`$, $`X=(h+81000)/(1-1/3000)`$ and $`n=\lceil X\rceil`$. Take $`K_X`$ to be the probability law of $`x=X`$ plus the sum of 40 independent uniforms on $`[-2000,2000]`$ and $`n`$ independent uniforms on $`[-1/3000,1/3000]`$. Its support is $`[X-R,X+R]`$ for $`R=80000+n/3000`$. Its centred exponential moment is

```math
\Phi(w)=\left(\frac{\sinh(2000w)}{2000w}\right)^{40}
        \left(\frac{\sinh(w/3000)}{w/3000}\right)^n .
```

Each factor at zero is interpreted by continuity. Direct algebra gives

```math
X-R>h+999,\qquad X+R<162000+(3001/2999)h.
\tag{7}
```

For $`\omega=0,\gamma`$ form $`M_\omega=\int f(x)e^{-i\omega x}\thinspace dK_X(x)`$ and put $`T_\omega(z)=J(z)e^{(z-i\omega)X}\Phi(z-i\omega)`$. Logarithmic derivatives in (3) contribute the following terms to $`M_\omega`$.

| Source | Term |
|---|---|
| Nontrivial or trivial zero $`\rho`$ of an odd factor | $`T_\omega(\rho-1/2)`$ |
| Nontrivial or trivial zero $`\rho`$ of a base-field factor | $`-T_\omega((\rho-1)/2)/2`$ |
| Pole of ordinary zeta at 1 | $`T_\omega(0)/2=0`$ |

To justify the formula, invert the absolutely convergent prime-power transforms on $`\Re z=1`$ and integrate against the compact measure. The fixed window has $`O((1+\vert \Im z\vert )^{-2})`$ decay there. For each primitive factor expand its logarithmic derivative by the genus-one Hadamard product for the completed function and the partial fractions for its gamma factor. Subtraction at a fixed real $`s_0\gt 1`$ gives terms $`(s-\rho)^{-1}-(s_0-\rho)^{-1}`$ for each zero, and the opposite term at a pole, up to a constant. The gamma terms give the trivial zeros $`-k,-k-2,\ldots`$ for nonprincipal factors, and $`-2,-4,\ldots`$ for zeta. The elementary factors of completed zeta remove the spurious zero at 0 and give its pole at 1.

For clarity, the inverse integral of $`T_\omega(z)/(z-w)`$ on this line equals $`T_\omega(w)`$ for $`\Re w\lt 1`$, and that of $`T_\omega(z)`$ equals zero. This follows by closing to the left. On $`\Re z=-U`$ the exponential decay from $`X\gt R+2`$ dominates $`J`$ and the boxes. The frequency integrals are uniformly controlled there by retaining two long-box denominators, giving a bound by a constant (for fixed $`X`$) times $`e^{-U(X-R-2)}\int_\mathbf R (1+\vert t-\omega\vert )^{-2}\thinspace dt`$ for $`U\ge1`$. Horizontal integrals vanish for each fixed $`U`$ by the same box decay. Alternatively the reciprocal kernel is the transform of $`e^{wu}\mathbf1_{u\gt 0}`$ and the entire averaging window stays on positive log heights. On the original line each subtracted partial fraction, summed absolutely, has at most polynomial growth in $`\vert \Im z\vert`$: split $`\vert \Im\rho\vert`$ at $`2\vert \Im s\vert +4`$, use the distance of $`\Re s`$ from the strip on the near terms, and apply (6) on both parts. The trivial zeros give an absolutely summable subtracted series there as well. This legitimizes integrating term by term. The square factor has the additional $`1/2`$ from its change of scale.

For a box of half-width $`s`$,

```math
\left|\frac{\sinh(sw)}{sw}\right|
\le e^{s|\Re w|}\min\lbrace1,1/(s|w|)\rbrace.
```

This follows from the box integral for the first bound and $`\vert \sinh(u+iv)\vert \le\cosh u`$ for the second. If $`\Re z\le0`$ and $`\vert z-i\omega\vert \ge d`$, (7) and the 40 long boxes give $`\vert e^{(z-i\omega)X}\Phi(z-i\omega)\vert \le2^{-40}`$. This covers all nontrivial base zeros by (5), and all odd-factor zeros on or left of the line apart from the selected term for $`\omega=\gamma`$. If an odd-factor zero gives $`\Re z\gt 0`$, its ordinate has magnitude above 8500 by Platt. Then $`\vert z-i\omega\vert /3000\gt (8500-4)/3000\gt e`$. The short boxes give a bound $`\exp((X+R)/2-n)\le\exp(-(X-R)/2)`$ on its contribution before multiplying by $`J`$. This deals with every higher zero anywhere in the critical strip.

From (2) the kernel weight is at most $`1000/(1+t^2)`$ for an odd-factor zero of ordinate $`t`$ and $`4000/(1+t^2)`$ for a base zero, before its factor $`1/2`$. Thus (6) bounds the total error from all unselected nontrivial zeros by $`4\cdot10^6(2^{-40}+e^{-999/2})\lt .001`$. Each sequence of trivial zeros starts at $`\Re z\le-1/2`$ and has spacing at least 1. Use $`\vert J(z)\vert \le 6e^{2\vert \Re z\vert }`$ and the box growth bound to get a total below $`128e^{-997/2}/(1-e^{-997})\lt .001`$. Finally (4) costs less than $`.001`$ on the entire support. The omitted finite primes cannot enter any of these windows since the smallest possible norm exceeds $`e^{X-R-2}\gt e^{216}`$. We have proved

```math
|M_0|<.01,\qquad |M_\gamma-iV|<.01 .
\tag{8}
```

## Two witnesses

The tests $`P_+(x)=1-\sin(\gamma x)`$ and $`P_-(x)=1+\sin(\gamma x)`$ are nonnegative. Each has mass at most $`1+2^{-40}`$ under $`K_X`$ because $`\vert \Phi(i\gamma)\vert \le 2^{-40}`$. Their $`f`$-integrals are $`M_0+\Im M_\gamma`$ and $`M_0-\Im M_\gamma`$, respectively. By (8) and $`V\gt 3.342843`$, the first exceeds $`1.9735`$ times its mass and the second is below $`-1.9735`$ times its mass. There is consequently a point of each sign and magnitude in the compact interval. No common phase or independence hypothesis on the zeros enters the argument.

For each point take $`y=e^x/(2\sqrt2)`$. The lower inequality of (7) gives $`y\gt H`$. Its upper inequality bounds $`2y=e^x/\sqrt2`$, proving (1). The inequality $`(1.001-3001/2999)\cdot500000000\gt 162000`$ proves the stated large-height simplification. [Lean proofs of the interval and positive-test implications](WaitingBound.lean)

## Finite interval and independent checks

For a rational prime $`p`$ outside the exclusions put $`u=\chi_{547}(p)`$ and $`v=\chi_{727}(p)`$. The residue degree $`f`$ in $`F`$ is the least $`j=1,2,4`$ with $`p^j\equiv\pm1\bmod16`$ and $`(uv)^j=1`$. There are $`8/f`$ ideals of norm $`p^f`$, each of sign $`u^f`$. The finite verifier sieves every relevant rational prime, forms coefficients $`-8u^f\log p`$ and integrates the two branches of the weight using prefix sums of the coefficients and the coefficients divided by the norms. It certifies moves on an integer grid starting at 1 with $`y_{\rm next}=\lceil9y/8\rceil`$. For real $`H`$ the next strictly later witness of each sign is constant between consecutive certified starting points. The endpoint ratio within each half-open interval therefore attains its maximum at the left endpoint. Checking all such endpoints and the upper endpoint covers the entire real range, including equality at a witness. The independent enumerator directly sums $`W`$ separately for the four example moves at the two extreme starting heights. [Continuous coverage](finite_windows.py) · [Direct enumeration](independent.py)

The independent zero check expands the conductor-547 Dirichlet series by residue classes after one period. In the Hurwitz tail centred at $`3/2`$, the $`k`$-th term for $`k\ge1`$ is
$`q^{-s}(-1)^k(s)_ k\zeta(s+k,3/2)\sum_{b=1}^q\chi(b)(b/q-1/2)^k/k!`$.
The constant term cancels. For $`\vert s\vert \lt 4`$, bound $`\vert (s)_k/k!\vert \le {k+3\choose3}`$, the Hurwitz factor by $`4(3/2)^{-k}`$ and the inner sum by $`q/2^k`$. After $`k=70`$ the tail ratio of the bounds is below $`1/2`$. This gives the implemented error and a second pair of opposite signs of the real completed quadratic function. [Independent calculation](independent.py)

The analytic estimates and explicit formula above are written proofs. The release replay checks the complex balls and direct counts, compiles the ten local Lean theorems and independently rechecks the compiled module with the kernel. The pinned field and trace release supplies its own replay. [Prior host receipt](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-short-shifts/evidence/check.txt) · [This release's receipt](evidence/check.txt)

## Both moves through exp(2000)

**Theorem.** For every real $`1\le H\le e^{1000}`$ both qualifying doublings finish before $`256H`$. For $`1\le H\le e^{2000}`$ they finish before $`400H`$. The magnitude in each direction is strictly greater than $`1.9735a`$, with $`a=\log2`$. Both starts must exceed $`H`$. The published direct sieve gives $`18H`$ through $`10^9`$, including the endpoint. Its proof uses actual prime ideals. It therefore suffices here to cover $`10^9\le H`$ in the two new ranges. [Direct sieve and continuous-height coverage](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-factor-18/PROOF.md)

Put $`r=\sqrt2`$, $`B=2r(r-1)^2`$, $`x=\log(2ry)`$ and $`f(x)=(S(2y)-S(y))/a`$. Write $`\mathcal A`$ for the eight odd primitive characters above. First use the signed prime-power score

```math
P(y)=-\frac{1}{B\sqrt y}\sum_{p\ {\rm prime}}\sum_{k\ge1}\log p
 \left(\sum_{\chi\in\mathcal A}\chi(p)^k\right) W(p^k/y).
```

Here the outer prime index ranges over all rational primes. At a good rational prime the local factorization identifies this with the sum over base prime ideals including every power, with sign $`-\epsilon_v^k`$ at the $`k`$-th ideal power. This can also be seen by summing the four even twists and the two quadratic character values. The shift from the signed prime powers to the good prime-only score is bounded below.

### Bounding square contributions without the base-field zero hypothesis

At ideal squares the sign $`\epsilon_v^2`$ is one. The correction to the level is therefore nonnegative. The difference of two such level corrections has absolute value at most the larger of their upper bounds. Let $`Z=\sqrt y`$, and use $`\vartheta(u)\le1.02u`$ for rational primes. The sum of the residue degrees above any rational prime is at most eight. Ideals of residue degree at least two in the square contribution lie above rational primes at most $`\sqrt{2Z}`$. Their total contribution to the normalized level is at most $`8.16\sqrt{2/Z}/B`$. A good degree-one rational prime must split in $`K`$, hence belongs to $`1`$ or $`-1`$ modulo 16. Discarding the remaining splitting condition is a valid upper bound since the weight is nonnegative. [Rosser and Schoenfeld, Theorem 9](https://doi.org/10.1215/ijm/1255631807)

For $`j=0,\ldots,31`$ set $`v=j/32`$ and $`u_j=2/\sqrt{1+v}-\sqrt{2/(2-v)}`$. The strict superlevel set $`W((p/Z)^2)\gt v`$ is an interval of length $`u_j Z`$ starting at $`Z\sqrt{2/(2-v)}`$. Integrate the superlevel counts in $`v`$ and upper-bound their integral by the left-hand sum. Brun–Titchmarsh on each of the two prime classes and $`\log p\le\log(2Z)`$ bound the degree-one part by

```math
Q_{\rm BT}(Z)=\frac4{32B}\sum_{j=0}^{31}u_j\frac{\log(2Z)}{\log(u_jZ/16)}.
```

Each interval has length greater than 32 already at the smallest allowed floor. Each displayed logarithmic ratio decreases with $`Z`$. [Montgomery and Vaughan, The large sieve](https://doi.org/10.1112/S0025579300004708)

For $`Z\ge10^{10}`$ there is another bound. Ramaré and Rumely's Theorem 1 bounds the error of each of the two progression theta sums by $`0.024u/8`$ for every $`u\ge Z`$. Integration by parts against $`W((p/Z)^2)`$ gives a degree-one bound $`2+8(0.024)/B`$. In detail the main term is $`2\int_1^2W(v^2)\,dv/B=2`$, since that integral is $`6r-8=B`$. The error is at most $`(2\cdot0.024/B)\int_1^2v\,|dW(v^2)|\le8(0.024)/B`$. Here the tent rises to one and falls to zero. In each height stratum the verifier uses the applicable upper bound and adds the higher-degree term above. Call the resulting bound $`Q(Y)`$ for both square levels whenever $`y\ge Y`$. [Ramaré and Rumely, Primes in arithmetic progressions](https://ramare-olivier.github.io/Maths/rumely.pdf)

For an ideal power of exponent $`k\ge3`$ use $`Nv\le(4y)^{1/k}`$ and $`W\le1`$. Bound $`k=3,4`$ separately and the remaining exponents by their number times the fifth-root bound. Summing first over rational primes gives, uniformly for $`y\ge Y\ge10^9`$,

```math
E(Y)=\frac{8.16}{B\sqrt Y}
\left((4Y)^{1/3}+(4Y)^{1/4}+\frac{\log(4Y)}a(4Y)^{1/5}\right).
```

Each term after division by $`\sqrt Y`$ decreases over this range, including the logarithmic fifth-root term. At the five omitted rational primes the sum of character coefficients has absolute value at most eight. There are at most three admissible exponents for $`p=2`$ and at most one for each of the other four primes. Their whole contribution is bounded by $`E_b(Y)=8(3\log2+\sum_{p\in\mathcal B,\ p\ne2}\log p)/(B\sqrt Y)`$, where $`\mathcal B`$ is the five-prime omission set. Thus replacing the prime-only move by the signed prime-power move costs at most $`(Q(Y)+2E(Y)+2E_b(Y))/a`$. All comparisons are uniform on the relevant averaging support. [Arithmetic bounds in the verifier](finite_deadlines.py)

### Signed positive averages from a finite zero verification

Use the same $`G(z)`$ and $`J(z)`$ as above. The explicit formula for the normalized move of $`P`$ has a term $`J(\rho-1/2)e^{(\rho-1/2)x}`$ per nontrivial odd-factor zero and the corresponding trivial-zero terms. It follows by the logarithmic-derivative inversion already proved, using only the odd factors because $`P`$ includes their prime powers. Every character is nonprincipal, so there is no pole term. Finite deletions are already charged in $`E_b`$.

For each witness centre $`x`$ use the convolution of $`n`$ uniform probability boxes of half-width $`h=1/3000`$ and two of half-width $`d`$. The support is $`[-R,R]`$ with $`R=nh+2d`$. Its transform at $`z`$ is

```math
\Phi(z)=\left(\frac{\sinh(hz)}{hz}\right)^n
        \left(\frac{\sinh(dz)}{dz}\right)^2 .
```

For the first range take $`(n,d)=(560,1/16)`$, for the second $`(1070,1/20)`$. Platt's unconditional finite computation places all zeros of the eight primitive factors through ordinate 8,000 on the critical line, since $`11632\cdot8000\lt10^8`$. For a higher zero, with $`z=\rho-1/2`$ and $`|\Re z|\le1/2`$, the short-box factor is at most $`\cosh(h/2)/(8000h)`$ and the two long boxes cost at most $`e^d`$. The previously proved strip estimate for $`J`$ and unconditional total zero mass therefore bound the full higher-zero contribution in absolute value by

```math
10^6 e^{x/2+d}\left(\frac{\cosh(h/2)}{8000h}\right)^n .
```

For each range the computation checks this at the largest permitted centre $`x=T+9`$, where $`T=1000`$ or $`2000`$. At the trivial zeros of these odd factors $`z=-3/2-2j`$. The kernel bound $`|J(z)|\le6e^{2|\Re z|}`$ gives a total at most $`48e^{-3(x-R-2)/2}/(1-e^{-2(x-R-2)})`$. The sum of both error bounds is below $`0.002`$. [Platt, Theorem 7.1](https://arxiv.org/pdf/1305.3087#page=14)

The table contains 6,000 positive sign brackets and the aggregate has their negative partners. At $`s=3/2+iv`$, the real part of the completed logarithmic derivative is a sum of positive terms $`(3/2-\beta)/((3/2-\beta)^2+(v-\gamma)^2)`$ over all true zeros, including any that lie off the line at high ordinate. The identity follows by subtracting the genus-one derivatives at reflected points and taking real parts. For $`v=0,3,\ldots,48`$ subtract the selected pairs. Each residual is below $`1/5`$. An omitted critical-line zero of absolute ordinate at most 48 would contribute at least $`1/(1+(3/2)^2)\gt1/5`$ at a covering centre, using conjugation on the negative half. At $`v=0`$ the residual is below $`0.127`$. Thus the unlisted critical-line ordinates satisfy $`\sum (1+t^2)^{-1}\lt0.127`$, with multiplicity. No positivity assertion at a point inside the critical strip is needed. [Exterior mass and sign verifier](finite_spectrum.py) · [Bracket integers](evidence/finite-zeros.tsv)

For $`|t|\gt48`$, the explicit $`J`$ gives $`|J(it)|\le50/(1+t^2)`$. The two long boxes contribute at most $`1/(48d)^2`$. All omitted critical-line terms together therefore cost at most $`50(0.127)/(48d)^2`$. Including high critical-line terms in both error budgets is harmless. For each retained positive ordinate the main signed pair contributes

```math
-\frac{4G(it)\sin(at/2)}a\,
   \mathrm{sinc}(ht)^n\mathrm{sinc}(dt)^2\sin(tx),
\qquad \mathrm{sinc}(v)=\frac{\sin v}{v}.
```

The two calculations evaluate every signed average using the true-frequency intervals. The first uses direct sines and sinc factors. The second uses a split phase, the half-angle form of $`G`$ and two-box cosine differences. Both subtract the arithmetic, omitted-line, high-zero and trivial-zero bounds. Each resulting signed average exceeds $`1.9735`$. Nonnegative unit mass then forces an actual move of the required sign somewhere in the support. [Signed averages](finite_deadlines.py) · [Numerical transcript](evidence/check.txt)

### Covering real starting heights

Write $`b=\log(2r H)`$. A witness centre $`x`$ covers every $`b`$ in the open interval $`(x+R-\log(C/2),x-R)`$. Indeed for any log centre $`u\in[x-R,x+R]`$ its start $`e^u/(2r)`$ exceeds $`H`$ and its end is below $`CH`$. For each sign separately the verifier overlaps successive open intervals, starting at $`\log(2r\cdot10^9)`$ and ending beyond $`T+\log(2r)`$. There are 526 witness rows at $`(T,C)=(1000,256)`$ and 1,017 at $`(T,C)=(2000,400)`$. Each row gives the centre in hundredths, the sign and the lower log-height stratum used in the arithmetic bounds. Stratum zero denotes the exact lower height $`10^9`$. All inequalities at the joins are strict. [First table](evidence/finite-1000.tsv) · [Second table](evidence/finite-2000.tsv) · [Lean endpoint deduction](WaitingBound.lean)

## A smaller explicit exponent at every height

**Theorem.** For every real $`H\ge1`$, both signed dyadic changes of magnitude greater than $`3.32\log2`$ finish strictly before $`e^{161020}H^{17001/16999}`$. The preceding explicit all-height interval used exponent $`3001/2999`$. The new estimate uses the same isolated zero, the same complex disk exclusions for the base-field zeros, and the same 40 long boxes. Here a positive Logan kernel controls the unknown high zeros. [Earlier all-height proof](https://doi.org/10.5281/zenodo.24410807) · [Büthe's use of the Logan kernel in finite-height prime estimates](https://arxiv.org/abs/1410.7015)

Let $`c=8550`$ and $`L\gt0`$. Take a uniform point $`(v_1,v_2,v_3)`$ on the unit two-sphere, tilt its probability law by $`e^{cLv_1}`$ and renormalize. The law of $`Lv_3`$ is a positive symmetric measure supported on $`[-L,L]`$. Its bilateral Laplace transform is

```math
\Psi_L(w)=\frac{c}{\sinh(cL)}
 \frac{\sinh(L\sqrt{c^2+w^2})}{\sqrt{c^2+w^2}} .
```

For real $`w`$, rotate the vector $`(c,0,w)`$ onto the polar axis and integrate the uniform marginal on $`[-1,1]`$. The power series in $`c^2+w^2`$ then gives the identity for all complex $`w`$, independently of the square-root choice. In particular $`\Psi_L(0)=1`$ and $`|\Psi_L(w)|\le e^{L|\Re w|}`$.

For $`w=u+iv`$ with $`|u|\le1/2`$ and $`|v|\ge8586`$, take a square root $`P+iQ`$ of $`c^2+w^2`$. Then $`Q^2=P^2+v^2-c^2-u^2\gt780^2`$ and $`|P|=|uv|/|Q|\lt6`$. For the second inequality use the decreasing function $`v/\sqrt{v^2-c^2-1/4}`$ at positive $`v\ge8586`$. When $`cL\ge1`$ we get

```math
|\Psi_L(w)|\le\frac{3c}{780}e^{-(c-6)L}<34e^{-8544L}.
```

Write $`h=\log H`$, and set

```math
L=\frac{h+162000}{16999},\qquad X=h+81000+L,\qquad R=80000+L .
```

Convolve the new measure with 40 uniform boxes of half-width 2,000, and translate to centre $`X`$. Then $`X-R=h+1000`$ and $`X+R\lt161020+(17001/16999)h`$. In particular its smallest physical prime window lies beyond every omitted norm. For $`\omega=0,\gamma`$ repeat the two moment integrals in the earlier proof.

Every unwanted base zero and every unselected odd zero on or left of the line at distance at least $`1/1000`$ from the test centre is suppressed by the forty long boxes. Each supplies a factor at most $`e^{2000|\Re z|}/2`$, and the total remaining exponential factor on that half-plane is at most one because $`X\gt R`$. The original interval and complex-disk tests establish the required separations, including at zero. For an odd zero to the right of the line, Platt's finite theorem forces $`|\Im z|\gt8590`$, since $`11632\cdot8590\lt10^8`$. As $`|\omega|\lt4`$, the new bound applies at $`w=z-i\omega`$. After including the centre and the long boxes, the absolute multiplier on this term is at most

```math
34\exp((X+80000)/2-8544L)\le34e^{-500}.
```

Indeed $`17000L=X+81000`$. The total mass estimate (6) and the strip bound (2) give an aggregate error below $`0.001`$ for the positive-real-part odd zeros. For the other unwanted zeros the earlier $`2^{-40}`$ estimate applies. The earlier trivial-zero and higher-power estimates apply with $`X-R\ge1000`$. The contour argument still has the integrable decay from the long boxes and the new kernel is entire with its support bound.

Consequently $`|M_0|\lt0.01`$ and $`|M_\gamma-iV|\lt0.01`$ with $`V\gt3.342843`$. The two nonnegative tests $`1\mp\sin(\gamma x)`$ each have mass at most $`1+2^{-40}`$. Their signed integrals exceed $`3.32(1+2^{-40})`$. Thus there is a qualifying move of each sign on the support. Its start is after $`H`$, and its completion is $`e^x/\sqrt2\lt e^{X+R}`$. The exponent and prefactor are effective at every starting height, independently of global GRH. [Kernel constants](power_window.py) · [Lean support deductions](WaitingBound.lean)

## Every positive excess exponent

**Theorem.** Fix any $`\delta\gt0`$. For all sufficiently large $`H`$, a gain and a loss larger than $`3.32\log2`$ finish before $`H^{1+\delta}`$. The threshold height here can depend on $`\delta`$. This is a field-specific application of Kaczorowski and Pintz's power-interval oscillation theorem together with the isolated coefficient already certified above. [Kaczorowski and Pintz, theorem (1.7)](https://personal.math.ubc.ca/~gerg/teaching/592-Fall2018/papers/1987.Kaczorowski_2.pdf)

First suppose that all the odd factors have their strip zeros on the critical line. By the unconditional prime ideal theorem the unsigned square term in the level tends to its constant bias, which disappears in the dyadic move. Higher ideal powers tend to zero. The move therefore differs by $`o(1)`$ at large heights from the absolutely and uniformly convergent, real, odd series in the centred log coordinate $`x`$. This argument for the square correction works for the fixed base field without an assumption on its zeros. At the selected frequency the Fourier coefficient has magnitude $`V\gt3.342843`$, so the uniform norm of the series is at least $`V`$ by averaging against the corresponding character. Oddness gives both signs above $`3.32`$ with a strict margin at some points. Uniform almost periodicity makes each sign with that margin relatively dense on the log axis. Every sufficiently long fixed log interval consequently supplies both signs in the original prime scores. This suffices for the growing log interval in the theorem. [Prime-only strip formula and square bias](https://github.com/leecrossley/theta-host-split-inert-race/blob/compact-filter-oscillation-v1.0.0/results/compact-filter-oscillation/PROOF.md#what-survives-on-sparse-grids)

Otherwise a reflected odd zero has $`\Re\rho\gt1/2`$. Let $`\Theta`$ be the supremum of $`\Re\rho-1/2`$ for such zeros, so $`\Theta\gt0`$. For the normalized prime-only move as a function of its starting height, the Mellin transform converges absolutely far to the right and is meromorphic on $`\Re z\gt0`$. Equation (3), multiplied by the move transform and a nonvanishing centering exponential, shows that its poles there are the shifted odd zeros. The base logarithmic derivative is regular at $`1+2z`$ and the higher-power remainder is holomorphic there. The numerator of the tent transform can vanish only at $`z=\pm1/2+2\pi i k/a`$. The move multiplier has its zeros on the imaginary axis. Hence all odd zeros in the open right half-strip survive, with their multiplicities adding. There are no real poles there by the finite zero verification. The poles form a discrete set, and at each the transform has an ordinary pole with nonzero residue.

Choose $`0\lt\eta\lt\Theta`$. Modify the actual normalized move in turn by adding and subtracting $`\mathbf1_{y\ge1}y^\eta`$. This only adds or subtracts $`1/(z-\eta)`$ to the Mellin transform. In a strip $`\Re z\gt\Theta-c_0`$ with $`0\lt c_0\lt\Theta-\eta`$ the singularities of both transforms are exactly the nonreal poles just described. They have no finite accumulation, their rightmost real-part supremum is $`\Theta`$, and their local pole expansions satisfy the theorem of Kaczorowski and Pintz. For every fixed $`0\lt\varepsilon\lt1`$ both modified real functions change sign in every sufficiently late $`[Y^{1-\varepsilon},Y]`$. Thus the original normalized move exceeds $`y^\eta`$ at one point and is below $`-y^\eta`$ at another in each such interval.

Set $`Y=H^{1+\delta/2}`$ and choose $`\varepsilon=(\delta/4)/(1+\delta/2)`$. The earliest possible start is at least $`H^{1+\delta/4}\gt H`$ and $`2Y\lt H^{1+\delta}`$ eventually. The two cases establish the claimed unconditional result. This proof supplies neither a fixed multiplicative deadline in the off-line case nor a numerical threshold height for arbitrary $`\delta`$.

## A shorter all-height bound

**Theorem.** For every $`H\ge1`$ both a gain and a loss of more than $`3.32\log2`$ finish before

```math
e^{219}H^{17001/16999}.
```

This uses a second zero and the 6,000-bracket table. It retains the exponent of the single-zero argument while reducing its prefactor from $`e^{161020}`$ to $`e^{219}`$. In particular both stronger moves finish before $`H^{1.001}`$ as soon as $`\log H\ge250000`$. The proof uses the published finite Dirichlet verification and treats all higher zeros throughout the strip. [Interval checks for the shorter bound](carrier_window.py)

### Cancellation in the prime squares

Put $`h=\log H`$ and first suppose $`h\ge100`$. At $`Z=\sqrt y`$ the signed square change has weight $`K(v)=W(v^2/2)/\sqrt2-W(v^2)`$ with $`v=N\mathfrak p/Z`$ and normalization $`1/(BZ)`$. Its support lies in $`[1,2\sqrt2]`$. Its sign changes at $`s=\sqrt{2\sqrt2}`$. The negative part peaks at $`\sqrt2`$ with height 1 and the positive part peaks at 2 with height $`1/\sqrt2`$. Substitution in the two integrals of $`W`$ gives $`\int K=0`$. With $`r=\sqrt2`$, integration of the negative part on $`[1,r]`$ and $`[r,s]`$ gives

```math
A=\int K_-=\int K_+=6r-2(1+r)s.
```

A good degree-one rational prime supplies at most eight base ideals and belongs to one of the two classes $`\pm1\bmod16`$. To bound either signed change, drop the other sign and bound the contributing prime ideals by those two entire classes. In Ramaré and Rumely's Table 1 the modulus-16 error is at most $`0.008994\,u/8`$ for $`u\ge10^{10}`$. Use $`0.009`$ and integrate by parts against the positive or negative part. For a unimodal test of height $`m`$ and peak at $`p`$ on $`[a_0,b_0]`$, the weighted variation is at most $`(p+b_0)m`$. For the two parts it is at most $`2+r`$. Thus the absolute degree-one square change costs at most $`(2A+2(0.009)(2+r))/B`$. Higher-degree squares are bounded by the larger of the two positive level caps above. Adding the higher powers and omitted primes and dividing by $`\log2`$ gives a total error below $`2.353`$ at all starts $`y\ge e^{100}`$. This smaller bound uses the signed square difference before the arithmetic estimate. [Ramaré and Rumely, Table 1](https://ramare-olivier.github.io/Maths/rumely.pdf#page=23) · [Arithmetic verification](carrier_window.py)

### Isolating a harmonic with a positive test

Let $`g`$ be the verified positive zero ordinate in the bracket with endpoints $`251203741273248/10^{14}`$ and $`251203741273252/10^{14}`$ of the factor with conductor 8752 and Conrey label 1093. Write $`M=13`$, $`w_k=1-|k|/14`$ for $`|k|\le13`$, and use the nonnegative Fejér polynomial

```math
F(\theta)=\frac1{14}\left|\sum_{j=0}^{13}e^{ij\theta}\right|^2
         =\sum_{k=-13}^{13}w_ke^{ik\theta}.
```

In the positive sphere construction above replace $`c`$ by 8520. Convolve its law with two probability boxes of half-width 54 and translate to centre $`X`$, where

```math
L=\frac{h+300}{16999},\quad R=108+L,\quad X=h+110+L.
```

Integrate the signed prime-power move against this probability times $`F(gx+\theta)`$, once for each choice $`\theta=\pm\pi/2`$. The positive and negative ordinates $`\pm g`$ paired with the cancelling harmonics contribute precisely $`\pm 2w_1|J(ig)|`$. The sign follows from $`J(ig)`$ being positive imaginary. All other contributions have absolute bounds independent of the centre and of the chosen sign.

For a real frequency difference $`v`$ the sphere transfer has modulus at most one. The two boxes supply $`U(v)=\mathrm{sinc}(54v)^2`$. For each of the 6,000 bracketed positive ordinates $`t`$, charge $`2|J(it)|\sum_k w_k U(t-kg)`$ with the selected summand omitted at $`(g,1)`$. The certificate obtains a main contribution above $`6.33`$ and total leakage from these terms below $`0.62`$. The integral of the Fejér weight is at most $`1+2\sum_{k=1}^{13}w_k U(kg)\lt1.00007`$. The reverse triangle inequality also bounds its mass below by $`0.99993`$. The unlisted critical-line terms have absolute ordinate greater than 48, and $`13g\lt33`$. Their combined leakage is at most

```math
\frac{50(0.127)\,14}{54^2(48-13g)^2}<0.00014 .
```

For an off-line zero the ordinate after shifting by any harmonic has absolute value at least $`8590-33=8557`$. With this frequency and $`c=8520`$, the earlier square-root calculation gives $`|Q|\gt700`$ and $`|P|\lt7`$. The sphere factor is therefore at most $`37e^{-8513L}`$. Including the two boxes and translation, the exponential has log at most $`(X+108)/2-8513L\le-41`$. Multiply by the total strip mass bound $`10^6`$, the constant 37 and the weight sum 14. The result is below $`0.0001`$. For trivial zeros use $`X-R-2=h\ge100`$ and the pointwise Fejér bound 14 in the earlier trivial-zero estimate. Together these two costs are below $`0.001`$. The two implementations evaluate $`U`$ by direct sinc and by a cosine difference on all the certified intervals. Subtracting the costs, dividing by the upper mass and subtracting $`2.353`$ leaves strictly more than $`3.32`$ in each direction. [Harmonic verifier and exact comparison](carrier_window.py)

The positive tests consequently locate a move of each sign. Every possible centre in the support has start above $`H`$ since $`X-R=h+2`$ and $`2\gt\log(2\sqrt2)`$. Every completion is below $`e^{X+R}\lt e^{219}H^{17001/16999}`$. For $`0\le h\le100`$, the two finite implementations additionally find a gain and a loss larger than $`3.32\log2`$ after $`e^{100}`$ and before $`e^{219}`$. This covers the remaining starting heights. Finally $`250000(1.001-17001/16999)\gt219`$ gives the stated power cutoff. [Finite witnesses and all-height verification](carrier_window.py) · [Support inequality in Lean](WaitingBound.lean)

## An unconditional log-squared wait

Write $`h=\log H`$, $`T=10^6+4h`$ and $`L(t)=\log(12000t)`$. For every real $`H\ge1`$ each of the two signs of the actual dyadic move occurs with the following strict amplitude and completion deadline.

| Magnitude of each signed move | Both doublings finish before |
|---|---|
| $`3.02\log2`$ | $`H\exp((13/5)L(T)^2)`$ |
| $`3.33\log2`$ | $`H\exp(10000+(13/5)L(T)^2)`$ |

Each start is strictly larger than $`H`$. The bound applies to the 547/727 host defined above, at $`b=r=2`$. The proof uses Platt’s unconditional finite Dirichlet computation and an unconditional upper bound on the total number of zeros. The averaging measures may depend on the actual positions of zeros in a finite strip. Their existence suffices for the stated deadlines. The proof needs no list of the possible off-line zeros. [Platt, Theorem 7.1](https://arxiv.org/abs/1305.3087) · [Bennett–Martin–O’Bryant–Rechnitzer, Corollary 1.2](https://doi.org/10.1090/mcom/3599)

### The positive two-point filter

Use the bilateral moment convention $`\widehat\mu(z)=\int e^{zu}\,d\mu(u)`$. For $`\delta\ge0`$ and $`\eta\gt0`$ take $`d=\pi/(2\eta)`$. Give the atoms at $`-d,d`$ probabilities proportional to $`e^{\delta d},e^{-\delta d}`$, respectively. Call this law $`P`$ and let $`\nu_{\delta,\eta}=P*\check P`$, where the check denotes reflection about zero. This symmetric probability has radius $`2d=\pi/\eta`$ and moment

```math
\widehat\nu_{\delta,\eta}(z)=
\frac{\cosh(d(z-\delta))\cosh(d(z+\delta))}{\cosh^2(d\delta)}.
```

It vanishes at every $`z=\pm\delta\pm i\eta`$. At a real frequency $`g`$ its response is

```math
\widehat\nu_{\delta,\eta}(ig)
=1-\frac{\sin^2(dg)}{\cosh^2(d\delta)}
\ge 1-d^2g^2 .
```

It is real and lies in $`[0,1]`$. A true deficit from one is at most $`d^2g^2`$. The product of finitely many responses in $`[0,1]`$ is at least one minus the sum of their deficits. Convolution adds the radii. Throughout the strip $`|\Re z|\le1/2`$, a probability of radius $`R`$ has moment at most $`e^{R/2}`$ in absolute value. [Formal cancellation and retention lemmas](WaitingBound.lean)

### How many cancellations are needed

Put $`D_0=8590`$. The eight primitive odd factors here have conductors between 547 and 11632. Platt puts every zero through $`D_0`$ on the line. Above it, the functional equation gives reflection about that line and the family is closed under complex conjugation. Count a representative with $`\beta\gt1/2`$ and $`\gamma\gt0`$ from each off-line quartet, including multiplicities. Write $`M(t)`$ for their number up to $`t`$. Every quartet supplies four zeros in the family at the same absolute ordinate. If different factors share a zero, counting each occurrence only increases the allowance.

For $`t\ge D_0`$, Corollary 1.2 of Bennett, Martin, O’Bryant and Rechnitzer bounds the zero count of one factor by $`(t/\pi)L(t)+0.247L(t)+6.894`$. Adding the eight factors gives at most $`3tL(t)`$, since $`8(1/\pi+(0.247+6.894/L(D_0))/D_0)\lt3`$. The corollary counts both signs of the ordinate. In particular $`M(t)\le(3/4)tL(t)`$. Partial summation, dropping the nonpositive lower endpoint term, gives the following for $`D_0\le U\le T`$.

```math
\sum_{U<\gamma\le T}\frac1\gamma
\le \frac34\left(L(T)+\frac{L(T)^2-L(U)^2}{2}\right).
```

The square-reciprocal tail satisfies

```math
\sum_{\gamma>U}\frac1{\gamma^2}
\le 2\int_U^\infty M(t)t^{-3}\,dt
\le \frac{3}{2U}(L(U)+1).
```

Both sums here and below are over the off-line quartet representatives. The latter estimate controls loss at the desired frequency even if zeros beyond $`U`$ turn out to be plentiful. [Total zero count](https://personal.math.ubc.ca/~bennett/BeMaObRe-MathComp-2021.pdf) · [Checks of the bounds](zero_filter.py)

### Filtering before and after modulation

Use the certified ordinate $`g=2.512037412732\ldots`$ and the order-13 nonnegative Fejér polynomial from the shorter all-height bound. In particular $`g\lt2.513`$ and $`13g\lt33`$. Use $`U=D_0`$ for amplitude 3.02 and $`U=200000`$ for amplitude 3.33. In both cases $`T\ge U`$. Write $`\delta=\beta-1/2`$ for a representative. Form two convolutions of probabilities.

1. For $`D_0\lt\gamma\le U`$, convolve one filter $`\nu_{\delta,\gamma+kg}`$ for each integer $`-13\le k\le13`$. Call the resulting probability $`P_0`$ and its radius bound $`r_0`$. All shifted ordinates are positive. When $`U=D_0`$ use a unit atom at zero.
2. For $`U\lt\gamma\le T`$, convolve one unshifted filter $`\nu_{\delta,\gamma}`$ per representative. Call this probability $`P_1`$ and its radius bound $`r_1`$. For an empty product use the same unit atom.

Put $`n=\lceil T/4\rceil`$, $`d_s=4/T`$. Take the law of the sum of two uniform boxes of half-width 54, $`n`$ independent uniform boxes of half-width $`d_s`$ and an independent draw from $`P_0`$. Translate its centre to $`X=h+2+R`$, where $`R=108+nd_s+r_0+r_1`$. Multiply the translated measure pointwise by $`F(gx+\theta)`$, then convolve with $`P_1`$. Use $`\theta=\pi/2`$ and $`-\pi/2`$. Each resulting measure is nonnegative with support in $`[h+2,h+2+2R]`$. Its mass lies between $`0.99993`$ and $`1.00007`$ by the same two long boxes as before, since the other probability factors have real-frequency modulus at most one and the last convolution preserves mass.

In the expansion for a single zero, multiplication by $`F`$ shifts its frequency by every $`kg`$. The factors in $`P_0`$ kill every shifted term of each quartet through $`U`$. For a negative ordinate use the opposite harmonic index. At ordinates from $`U`$ through $`T`$ the final convolution kills the *unshifted* moment of each quartet, after the harmonic terms have been added. For the selected critical-line pair the matched harmonic evaluates all the inner moments at zero. Its contribution is the earlier $`\pm2w_1|J(ig)|`$ multiplied by the real response $`\widehat P_1(ig)`$. All remaining critical-line terms keep the same upper leakage bound. The order of the multiplication and the two convolutions is used in each of these statements.

### The tail throughout the strip

For an off-line zero with $`|\gamma|\gt T`$, each short box at a shifted exponent $`z+ikg`$ contributes at most $`e^{d_s/2}/(d_s(T-33))`$. Indeed the moment is $`\sinh(d_s(z+ikg))/(d_s(z+ikg))`$, and $`|\sinh(u+iv)|\le e^{|u|}`$. Since $`4(1-33/T)\gt e`$, the $`n`$ boxes cost at most $`e^{nd_s/2-n}`$. Accounting for the centre, the long boxes and the two filters, the cost of the whole off-line tail is at most

```math
14\cdot10^6 \exp(h/2+1+R-n).
```

Here 14 is the sum of the Fejér coefficients and $`10^6`$ bounds the strip sum of the absolute odd zero coefficients in the preceding proof. The absolute convergence on each compact support justifies convolution term by term. For the trivial zeros one can integrate the earlier pointwise bound on this same support, with the Fejér weight at most 14. The left endpoint in the centred coordinate remains $`h+2`$. For $`h\ge100`$ their cost and the high off-line cost together are below $`0.001`$. The following radius estimates establish the latter cost uniformly, including the first eligible height.

For the larger cutoff the first reciprocal bound gives

```math
r_0\le27\pi\frac{8590}{8590-33}\frac34
\left(L(U)+\frac{L(U)^2-L(D_0)^2}{2}\right)<5406.
```

for the smaller one $`r_0=0`$. The post-modulation radius satisfies

```math
r_1\le \tfrac{13}{10}L(T)^2-\tfrac{3\pi}{8}L(U)^2 ,
```

since $`(3\pi/4)(1/2+1/L(10^6))\lt1.3`$. Also $`nd_s\le1+4/10^6`$. At $`h=0`$ the crude upper bound $`1+108+1.000004+5406+1.3L(10^6)^2-250000`$ is below $`-100`$. As $`h`$ grows, the derivative of $`1.3L(T)^2`$ is $`10.4L(T)/T\lt1/2`$ and decreases. Since $`n\ge250000+h`$, this proves $`h/2+1+R-n\lt-100`$ throughout. The displayed off-line cost is therefore below $`0.0001`$.

### The two margins and completion factors

At the carrier the product bound gives

```math
\widehat P_1(ig)\ge1-\left(\frac{\pi g}{2}\right)^2
\frac{3(L(U)+1)}{2U}=1-\varepsilon_U .
```

From the two preceding harmonic computations the selected coefficient exceeds $`6.33`$, the bracketed leakage is below $`0.62`$, the unlisted critical-line leakage is below $`0.00014`$ and the signed prime-square and other arithmetic corrections cost less than $`2.353`$ in normalized units. The same bounds apply to the present positive averages. Dividing by the upper mass for either signed inequality, each normalized prime move exceeds in the chosen direction

```math
\frac{6.33(1-\varepsilon_U)-0.62-0.00014-0.001}{1.00007}-2.353 .
```

This is greater than $`3.02`$ at $`U=8590`$ and greater than $`3.33`$ at $`U=200000`$. A positive measure of nonzero mass locates an actual start with that signed move. For $`0\le h\lt100`$, the two finite computations give both stronger witnesses after $`e^{100}`$ and before $`e^{219}`$. [Carrier and early witnesses](carrier_window.py) · [Trivial and square corrections](PROOF.md#cancellation-in-the-prime-squares)

The full support lies strictly after the original starting height because $`2\gt\log(2\sqrt2)`$. The endpoint of any permitted doubling is strictly less than $`e^{X+R}`$. Now $`X+R\lt h+221+2r_0+2r_1`$. At $`U=8590`$ one has $`221-(3\pi/4)L(U)^2\lt0`$. At $`U=200000`$ one has $`221+2r_0-(3\pi/4)L(U)^2\lt10000`$. These give both deadlines, including the finite bridge. At $`h=2\cdot10^9`$ and $`2\cdot10^{10}`$ the stronger log overhead is below $`10^{-5}h`$ and $`10^{-6}h`$ respectively. Its derivative $`(104/5)L(T)/T`$ is smaller than the corresponding excess at each cutoff and decreases thereafter. [Rational and ball checks](zero_filter.py) · [Kernel replay](evidence/check.txt)

### What zero counts buy

For fixed conductors the same zero-count theorem gives the sharper asymptotic $`M(t)\le(2/\pi)t\log t+O(t)`$, by symmetry from the eight factors. Partial summation gives $`r_1\le(\log T)^2+O(\log T)`$. Thus for either amplitude and every fixed $`\varepsilon\gt0`$ both completions eventually lie before $`H\exp((2+\varepsilon)(\log\log H)^2)`$. The leading 2 is an upper bound on the proved logarithmic overhead. Optimality is open.

If the number of off-line quartet representatives through $`t`$ is $`O(t)`$, the same argument gives $`r_1=O(\log T)`$ and a deadline of the form $`CH(1+\log H)^K`$ for fixed constants. If their positive-ordinate reciprocals have finite sum, the radii of all the finite post-modulation convolutions are uniformly bounded by $`\pi\sum1/\gamma`$. The low shifted radius is already bounded. The same construction then gives a *fixed* completion factor for both $`3.33\log2`$ moves at every height. This sufficient hypothesis allows an infinite exceptional set. The unconditional result uses the total zero count instead.

## A polylogarithmic waiting factor

**Theorem.** Put $`h=\log H`$, $`T=10^6+4h`$ and $`L=\log(12000T)`$. For every real $`H\ge1`$ the two strict $`3.33\log2`$ moves finish before

```math
H e^{225}(12000T)^4(3072000\,L^4)^{540}.
```

The two starts exceed $`H`$. In particular for each $`\varepsilon\gt0`$ both ends eventually precede $`H(\log H)^{4+\varepsilon}`$. The exponent 5 is one example. The explicit formula gives the result for all heights, including those preceding any asymptotic comparison. The individual-zero bounds above can be smaller at moderate heights. [Arithmetic and interpolation checks](separated_zeros.py)

### A positive density for separated complex frequencies

Use a probability on $`[-\pi a/2,\pi a/2]`$ with density $`\cos(x/a)/(2a)`$, where $`0\lt a\le1`$. Its bilateral moment at $`z`$ is $`B(az)`$, with $`B(w)=\cosh(\pi w/2)/(1+w^2)`$ by direct integration, extending at the two removable singularities.

Suppose a finite set of real frequencies $`\lambda_i`$ is closed under negation, adjacent frequencies are at least $`s/a`$ apart with $`s\ge2`$, and $`a|\lambda_i|\ge x_0\ge2`$. Assign equal offsets $`0\le\delta_i\le1/2`$ to opposite frequencies. There is a positive probability on the same interval whose moments vanish at every $`\delta_i+i\lambda_i`$ whenever

```math
q=\frac{4\cosh(\pi a/4)}{s^2}<1,\qquad
8\left(x_0^{-2}+(s x_0)^{-1}\right)<1-q.
```

To prove it multiply the cosine density by $`1+\sum_j c_j e^{-i\lambda_j x}`$. The required equations are $`Ac=-b`$, with $`A_{ij}=B(a(\delta_i+i(\lambda_i-\lambda_j)))`$ and $`b_i=B(a(\delta_i+i\lambda_i))`$. The diagonal $`D_{ii}=B(a\delta_i)\ge1`$ because the original law is symmetric. For real $`v,t`$ the squared modulus of $`1+(v+it)^2`$ equals $`(t^2-1)^2+2v^2(1+t^2)+v^4`$. Thus, if $`|t|\gt1`$ and $`|v|\le a/2`$, the absolute moment is at most $`\cosh(\pi a/4)/(t^2-1)`$. In each column the distances to the other frequencies on either side are at least $`ns/a`$ for successive integers $`n\ge1`$. Since $`s\ge2`$,

```math
2\sum_{n\ge1}\frac{1}{n^2s^2-1}
\le\frac{4}{s^2}\,2\sum_{n\ge1}\frac{1}{4n^2-1}=\frac4{s^2}.
```

The last series telescopes. In the induced column-sum norm $`\|D^{-1}(A-D)\|_1\le q`$. At the positive frequencies the scaled ordinates are at least $`x_0+ns`$ for $`n\ge0`$. Use $`\cosh(\pi a/4)\lt2`$, $`1/(t^2-1)\le2/t^2`$ for $`|t|\ge2`$, and the integral bound on the decreasing reciprocal squares. Including both signs gives $`\|b\|_1\le8(x_0^{-2}+(s x_0)^{-1})`$. The Neumann series now gives a unique solution with $`\sum_j|c_j|\lt1`$. Conjugating the system interchanges opposite indices, so the multiplier is real. It is strictly positive throughout. Dividing the resulting measure by its positive mass preserves the zeros. The empty set requires only a unit atom. [Two finite matrix and integral checks](separated_zeros.py)

### Counting groups and controlling the left half-strip

For an upper cutoff $`W\ge8590`$, write $`L_W=\log(12000W)`$. In a closed interval of length $`d`$ inside $`[8590,W]`$ the number of off-line quartet representatives is at most

```math
\frac{2d}{\pi}L_W+0.988 L_W+27.576 .
```

Subtract Corollary 1.2 of Bennett, Martin, O’Bryant and Rechnitzer at the endpoints, using the upper bound at the right and lower bound at the left. For a closed left endpoint take a limit from below. The derivative of the main term is $`\log(qt/(2\pi))/\pi\le L_W/\pi`$. Each of the eight factors costs two error bounds $`0.247L_W+6.894`$. Divide by four for the quartet symmetry. Intersecting an interval with the ordinate range only shortens it. In a finite list of ordinates sorted with multiplicity, assign each in turn to an available group whose last ordinate is at least $`d`$ behind. The largest number of list entries in a closed interval of length $`d`$ bounds the number of groups, including coincident ordinates. Each group has spacing at least $`d`$. This only uses the finite local zero count. [Zero-count theorem](https://doi.org/10.1090/mcom/3599)

Only the growing zeros need interpolation. For a left-hand off-line zero write $`z=\delta+i\gamma`$ with $`\delta\lt0`$. The translated probability averages, including the final convolution, have centre $`X=h+2+R`$ and summed radius $`R`$. The two long boxes give, at a harmonic shift $`v=\gamma+kg`$, an extra factor at most $`(54|v|)^{-2}`$ after extracting their support exponentials. All translation and support exponentials together are at most $`e^{\delta(X-R)}\le1`$. All these ordinates have $`|\gamma|\gt8590`$ and hence $`|v|\gt8557`$. The sum over the entire decaying half, including the Fejér weights and the absolute strip mass, is below $`14\cdot10^6/(54\cdot8557)^2\lt0.00007`$. This estimate holds with the post-modulation convolution because it is a positive probability of its stated radius. It removes any need to impose two nearly coincident real offsets on the interpolation system.

### Two widths for the growing zeros

Take $`U=\min(T,256L^4)`$. Here $`L\gt23`$ and $`U\ge10^6`$. Group the growing representatives through $`U`$ with $`d=4`$. The preceding count gives fewer than $`5L_U`$ groups since $`8/\pi+0.988+27.576/23\lt5`$. For each group and each of the 27 harmonics take the positive density above with $`a=1,s=4`$. Use the shifted positive ordinates $`\lambda=\gamma+kg`$ together with their negatives and the common positive offsets. They are separated by 4 and exceed 8557 in magnitude. The negative-frequency member of a quartet is covered by the opposite harmonic index. The norm condition holds with $`q\lt1/2`$ and $`x_0=8557`$. Convolve these probabilities *before* multiplying by the Fejér test. Every shifted growing zero through $`U`$ is killed, and the total radius satisfies $`r_0\le270L_U`$ since $`\pi/2\lt2`$.

For $`U\lt\gamma\le T`$ use one group per separated list of *unshifted* positive ordinates and their negatives. Set

```math
a=\frac{1}{256L^2},\qquad s=2+\frac1L,\qquad d=s/a.
```

The group count satisfies $`m\le(2d/\pi)L+3L`$. If this list is nonempty then $`aU=L^2`$, the reflected list is still $`d`$-separated, and one can use $`x_0=L^2`$. Since $`\cosh(\pi a/4)\le1+a^2`$ for $`0\le a\le1`$, the column-sum bound gives $`q\le1-1/(2L)`$. One way to check it is to put $`u=1/L`$, bound $`a^2\le u^4`$ and expand $`(1-u/2)(2+u)^2-4(1+u^4)=u(2-u-u^2/2-4u^3)\gt0`$. The cosine-hyperbolic estimate follows from its positive power series and $`\cosh(1)\lt2`$. The right-hand side norm times $`2L`$ is at most $`16/L^3+8/L\lt1`$. The probabilities therefore exist. Convolve them *after* modulation. They annihilate all remaining growing zeros through $`T`$ with radius

```math
r_1\le(\pi a/2)m\le sL+(3\pi/2)aL<2L+2 .
```

For any probability on $`[-\pi a/2,\pi a/2]`$, its response at the selected frequency $`g`$ has real part at least $`1-g^2(\pi a/2)^2/2`$ by the elementary cosine inequality. The modulus of the convolution response is at least one minus the sum of these deficits, hence greater than

```math
1-\frac{g^2\pi a}{4}(2L+2)>0.998 .
```

An asymmetric probability can give a complex response. Choose the two Fejér phases to align the matched conjugate pair with this response, with opposite signs for the two tests. The main contributions then have magnitude above $`6.33(0.998)`$. Neither mass nor any leakage estimate used the phase. In particular all shifted critical-line terms still use the certified two-box bound. [Carrier verification](carrier_window.py)

### Completion, and what would give a fixed factor

Keep the original $`n=\lceil T/4\rceil`$ short boxes and the two long boxes. Now $`R\le108+1.000004+272L+2`$ since $`U\le T`$. The high growing tail has the same bound $`14\cdot10^6\exp(h/2+1+R-n)`$. At $`h=0`$ the crude exponent is below $`-100`$. Its upper envelope decreases thereafter because $`n\ge250000+h`$ and $`4\cdot272/T\lt1/2`$. The high and trivial-zero allowance of 0.001 and the new decaying allowance of 0.00007 together fit inside 0.002. Including both line-leakage bounds and the prime-square correction, the signed normalized margin at $`h\ge100`$ is greater than

```math
\frac{6.33(0.998)-0.62-0.00014-0.002}{1.00007}-2.353>3.33 .
```

The support starts after $`H`$ and every completion precedes $`e^{X+R}`$. The radius bounds give

```math
X+R\lt h+221+2r_0+2r_1
\lt h+225+540\log(3072000L^4)+4L .
```

The independently verified finite pair covers $`0\le h\lt100`$. Since $`T=10^6+4\log H`$, the displayed deadline grows as $`O(H(\log H)^4(\log\log H)^{2160})`$. Every larger fixed logarithmic exponent absorbs the last factor and the constant. This proves the theorem.

A separate sufficient condition for a fixed completion factor is a uniform bound on the number of off-line quartet representatives in each interval of length four. In that case use only the first, shifted construction through $`T`$. It costs at most 27 probabilities of radius $`\pi/2`$ per group, and the group count is uniformly bounded. Choosing the constant part of $`T`$ larger if needed makes the high tail negligible at every $`h\ge100`$. The same finite pair handles the smaller heights. A bounded local count can allow a sequence of exceptional ordinates with divergent reciprocal sum. More generally, a union of such a locally bounded set and a set with summable reciprocals suffices. Apply the separated construction to the former and the 27 shifted two-point factors to the latter before modulation. Their total radius is finite since $`\gamma/(\gamma-33)`$ is bounded for $`\gamma\gt8590`$. None of these sparsity conditions enters the unconditional polylogarithmic bound. The next argument proves unconditional existence by changing the normalization when an off-line zero is present.

## A fixed factor without a global zero hypothesis

**Theorem.** There exists a finite constant $`C`$ such that, after every real $`H\ge1`$, there are $`y_+,y_-\gt H`$ with $`2\max(y_+,y_-)\lt CH`$ and

```math
S(2y_+)-S(y_+)>3.33\log2,\qquad S(2y_-)-S(y_-)<-3.33\log2 .
```

This is unconditional. The proof of existence gives no numerical value of $`C`$ valid in both cases below. The polylogarithmic bound gives a computable deadline for each $`H`$. The finite $`400H`$ certificate and the conditional four-digit obstruction refer to the smaller $`1.9735\log2`$ target.

Let $`\Theta=\sup_\rho(\Re\rho-1/2)`$ over strip zeros of the eight odd factors. If $`0\lt\eta\lt\Theta`$, there is also a finite $`C_\eta`$ such that the two normalized moves exceed $`H^\eta`$ in opposite directions after every sufficiently large $`H`$, both finishing before $`C_\eta H`$. The factor may depend on $`\eta`$ and the field.

### Selecting an ordinate and summing the support radii

Suppose an odd factor has an off-line zero. The reflection and conjugation symmetries give a zero at $`\sigma+ig`$ with $`1/2\lt\sigma\lt1`$ and $`g\gt0`$. For the power assertion first select an ordinate at which some real part exceeds $`1/2+\eta`$. At the selected ordinate take the **largest** real part among all eight odd factors and call it $`\sigma`$. This maximum is over a finite nonempty set. Set $`d=\sigma-1/2\gt0`$. The selection rules out a farther-right zero at exactly the same ordinate. It needs no zero at the global supremum. The coefficient $`A=mJ(d+ig)`$ is nonzero, where $`m`$ is the sum of the multiplicities there. All odd factors enter with the same sign in (3). The zeros of the tent numerator have real parts $`\pm1/2`$ and those of the move multiplier have real part zero. [Transform and arithmetic factors](PROOF.md#the-transform-and-the-prime-powers)

Montgomery’s Dirichlet zero-density theorem with fixed conductor cap $`Q=11632`$ gives

```math
\sum_{\chi} N_\chi(\sigma,t)
\ll_Q t^{\,3(1-\sigma)/(2-\sigma)}(\log(Qt))^{13}\qquad (t\ge2).
```

Here the sum can be restricted to the eight factors by positivity. The exponent is strictly less than one because $`\sigma\gt1/2`$. A dyadic summation therefore gives $`\sum_{\beta\ge\sigma,\,\gamma\gt0}1/\gamma\lt\infty`$, counting multiplicities. There are no real odd zeros off the critical line by the finite verification. [Montgomery, *Zeros of L-functions*](https://www.mathnet.ru/eng/mat564) · [Explicit statement, Lemma 1](https://www.numdam.org/item/10.5802/jtnb.1162.pdf#page=11)

For every positive ordinate of a zero $`\beta+i\gamma`$ with $`\beta\ge\sigma`$, and every $`k\in\{-1,0,1\}`$, put $`v=\gamma+kg`$. If $`v\ne0`$ include the positive two-point convolution $`\nu_{\beta-\sigma,|v|}`$ from above. If $`v=0`$ then $`\gamma=g,k=-1`$ and maximality at this ordinate forces $`\beta=\sigma`$. Leave this matched harmonic and its conjugate in place. For $`\gamma\gt2g`$ each radius is at most $`2\pi/\gamma`$. At smaller ordinates only finitely many positive radii occur. The sum of all radii is finite. Hence the infinite convolution is a compactly supported probability. Indeed, independent draws from its factors form an absolutely convergent series on every outcome, and the sum is bounded by the sum of the radii. At any targeted exponent its moment vanishes by factoring off that one independent summand. In particular no convergence of a conditionally ordered product of complex values is required. Convolve with two additional independent uniforms on $`[-1,1]`$, to provide integrable Fourier decay. Denote the resulting centred probability by $`P`$ and its finite radius bound by $`R`$. [Positive two-point formula](PROOF.md#the-positive-two-point-filter)

### The tilted prime move

Write $`f(x)=(S(2y)-S(y))/\log2`$ with $`y=e^x/(2\sqrt2)`$. Translate $`P`$ to $`X=h+2+R`$ where $`h=\log H`$ and integrate $`e^{-dx}f(x)`$ against $`1\pm\cos(gx+\theta)`$ times the translated probability. These are nonnegative tests of mass at most two. For $`h`$ large the absolute-convergence argument for (2), (3) and (6) applies on this compact support. More explicitly, a nontrivial odd zero contributes at a harmonic $`k`$ its coefficient $`J(\rho-1/2)`$, multiplied by the translation and the centred moment at $`(\beta-\sigma)+i(\gamma+kg)`$ and the harmonic coefficient. The two extra boxes justify the same contour inversion. The tilt changes the test exponent and leaves $`J(\rho-1/2)`$ evaluated at the original zero.

All terms with $`\beta\ge\sigma`$ vanish except the two matched harmonics. For negative ordinates take the opposite harmonic index. The surviving exponent is zero, so the centred moment and the translation both equal one. The surviving pair contributes $`\pm\Re(Ae^{-i\theta})`$ and can be made $`\pm|A|`$. This also handles multiple zeros at the same point.

For $`\beta\lt\sigma`$ the absolute translation-and-moment factor at each harmonic is bounded by $`\exp(-(\sigma-\beta)(X-R))`$. It tends to zero for each such zero and is at most one. Equations (2) and (6) give the summable majorant $`\sum_\rho|J(\rho-1/2)|\lt10^6`$. Dominated convergence for this series gives an $`o(1)`$ total as $`h\to\infty`$, even if real parts accumulate at $`\sigma`$ from below. The unweighted signed square correction, higher prime powers and excluded primes are bounded for $`x\ge h+2`$ by the earlier arithmetic estimates, so multiplying by $`e^{-dx}`$ sends their contribution to zero. The odd trivial-zero terms also tend to zero uniformly there. Equivalently the base terms in (3) have real exponent at most zero before tilting. This gives the two signed integrals $`\pm|A|+o(1)`$ for the actual prime-only move.

Eventually the integral in each chosen direction exceeds $`|A|/2`$. Positivity and mass at most two give a start in the support satisfying $`\pm e^{-dx}f(x)\gt |A|/4`$. All such starts obey $`y\gt H`$ and $`2y\lt e^{2+2R}H`$. Since $`x\ge h+2`$, the two signed magnitudes are at least a positive constant times $`H^d`$. This proves the power assertion whenever $`d\gt\eta`$ and eventually gives both strict $`3.33\log2`$ moves even for an arbitrarily small positive $`d`$.

### The critical-line case and the finite bridge

If all zeros of the odd factors are on the critical line, there is no growing or decaying off-line contribution. Use the two uniform boxes of half-width 54 and the degree-13 Fejér weight at the certified ordinate. A centre $`X=h+110`$ and radius 108 suffice. The matched contribution above 6.33, the two critical-line leakage bounds, the mass bound and the signed arithmetic cost from the harmonic calculation give both normalized moves greater than 3.33 in absolute value, in opposite directions for $`h\ge100`$. The finite pair after $`e^{100}`$ with both completions before $`e^{219}`$ covers $`0\le h\lt100`$. Thus $`e^{219}`$ is one choice for the stronger target in this case. [Carrier margins](carrier_window.py)

In the off-line case take any finite $`H_0`$ beyond which the tilted tests already work. The explicit all-height deadline at $`H_0`$ supplies a single gain and loss starting after $`H_0`$ and finishing before a finite number $`B`$. For every $`1\le H\le H_0`$ the same pair qualifies, with endpoint below $`BH`$. Taking $`C=\max(B,e^{2+2R})`$ covers every height. The possible dependence of $`R`$ on a small gap between two ordinates and of $`H_0`$ on an unknown real-part separation explains why this argument gives no unconditional numerical factor.

Consequently bounded normalized moves and power-sized moves in every sufficiently late window of **some fixed multiplicative width** distinguish the two spectral cases. More precisely, the supremum of the positive $`\eta`$ for which *both* signs exceed $`H^\eta`$ in such windows is $`\Theta`$, with value zero when no positive exponent qualifies. For the upper direction the absolutely convergent strip sum is $`O(e^{\Theta x})`$ and the prime-only corrections are bounded. The same pointwise estimate follows directly by applying the finite-measure identity to smooth approximations to a point: on each compact positive interval (2) and (6) give uniform convergence, and the tent move is continuous. A fixed multiplicative window changes this upper bound by only a constant. Earlier results characterized a limsup on sparse observations. Here each sufficiently late starting height supplies both signs. [Earlier sparse sampling result](https://github.com/leecrossley/theta-host-split-inert-race/blob/compact-filter-oscillation-v1.0.0/results/compact-filter-oscillation/PROOF.md#what-survives-on-sparse-grids)

## Both moves through exp(12500)

**Theorem.** For every real $`1\le H\le e^{12500}`$ both dyadic changes of magnitude greater than $`1.9735\log2`$, one in each direction, occur at starts strictly after $`H`$ and finish strictly before $`400H`$. This uses the field, sign, width, and five omitted rational primes of the statement above. The earlier finite result covers $`H\le e^{2000}`$. It remains to treat $`2000\le\log H\le12500`$. The argument below uses Platt’s finite zero verification and the explicit formula throughout the rest of the strip. [Earlier finite certificate](#both-moves-through-exp2000) · [Arithmetic trace](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-short-shifts/RaceTrace.lean#L83)

### A positive test at each recorded centre

Split the log-height range into 83 closed intervals $`[\ell,T]`$ of length 128, shortening the final interval. Take $`c=8500`$, $`L=(T+82)/(2\cdot8492)`$ and $`d=1/k`$ with $`k\in\{16,24,32,48,64\}`$. Convolve the projected tilted-sphere probability of the earlier proof, of radius $`L`$, with two uniform probability boxes of half-width $`d`$. The convolution is positive, has mass one and is supported on $`[-R,R]`$ for $`R=L+2d`$. Its bilateral transform is

```math
\Psi_L(z)\left(\frac{\sinh(dz)}{dz}\right)^2,\qquad
\Psi_L(z)=\frac{c\,\sinh(L\sqrt{c^2+z^2})}{\sqrt{c^2+z^2}\,\sinh(cL)}.
```

The continuous values at removable singularities are understood. For each recorded centre $`x`$ we have $`x\lt T+7`$ and $`x-R\gt2000`$. Coordinates are $`x=\log(2\sqrt2\,y)`$ for a doubling starting at $`y`$. The positive average of the normalized prime-only move is estimated in each direction. [Probability construction](#a-smaller-explicit-exponent-at-every-height) · [Sphere calculation](finite_sphere.py)

For an odd-factor zero $`z=u+iv`$ measured from the critical line, $`|u|\le1/2`$. All with $`|v|\le8590`$ are on the line by the same finite theorem as before. For $`|v|\ge8590`$ write $`P+iQ=\sqrt{c^2+z^2}`$. Then $`Q^2\ge v^2-c^2-1/4\gt1200^2`$ and $`|P|=|uv|/|Q|\lt8`$, using monotonicity of $`v/\sqrt{v^2-c^2-1/4}`$ for positive $`v`$. Since $`cL\gt1`$, the sphere factor is at most $`30e^{-8492L}`$. Translation costs at most $`e^{x/2}`$ and the two boxes at most $`e^d`$. The total odd strip coefficient mass is less than $`10^6`$. The high-zero bound is therefore $`30\cdot10^6 e^{x/2+d-8492L}`$, whose exponent is below $`-37`$. In the odd trivial-zero bound the earliest support exceeds 2000, so the high-zero and trivial-zero costs together are less than $`0.001`$. This controls every possible off-line zero. [Total coefficient mass](#total-mass-of-the-omitted-zeros) · [Platt, Theorem 7.1](https://arxiv.org/pdf/1305.3087#page=14)

### The omitted line spectrum

The same 6,000 genuine positive brackets and their partners are used. To extend completeness from absolute ordinate 48 through 63, take five further exterior Poisson centres $`v=51,54,57,60,63`$. Evaluate at $`3/2+iv`$ and subtract the two positive terms contributed by each selected bracket and its partner. The resulting enclosures are positive and below $`1/4`$ at all five centres. Any omitted line ordinate within distance $`3/2`$ contributes at least $`4/13`$ to the same positive sum, a contradiction. The existing seventeen centres cover the lower strip and at zero give a total unlisted line mass $`\sum(1+t^2)^{-1}\lt0.127`$. Multiplicity and shared ordinates are counted. [Exterior logarithmic derivatives](finite_spectrum.py) · [Five additional enclosures](finite_sphere.py)

For $`|t|\ge D=63`$, the imaginary-axis sphere response has modulus at most $`\Psi_L(iD)`$. For $`|t|\le c`$ this follows because $`\sinh(Lp)/p`$ is increasing for $`p\ge0`$, by its power series, and $`p=\sqrt{c^2-t^2}`$ decreases. Beyond $`c`$, the sine ratio is bounded by its value at $`p=0`$. Combining this with $`|J(it)|\le50/(1+t^2)`$ and the boxes gives an omitted-line allowance of

```math
50(0.127)\,\Psi_L(iD)\min\{1,(Dd)^{-2}\}.
```

The selected pairs with $`t\gt190`$ are also charged by the sum of their individual absolute transformed amplitudes. The 2,364 remaining pairs are evaluated with their signs. Separately, the signed square difference, higher prime powers and all omitted primes cost less than $`2.353`$ in normalized units everywhere on these supports. The arithmetic estimate uses the two rational residue classes modulo 16 and the signed square kernel. [Arithmetic cancellation](#cancellation-in-the-prime-squares) · [Independent arithmetic verifier](carrier_window.py)

### Cover of every real starting height

Each table row consists of the block endpoint $`T`$, the sign, the numerator $`m`$ of the centre $`x=m/64`$, and $`k=1/d`$. For every row the signed retained sum exceeds $`1.9735`$ plus all four allowances for the omitted line terms, the selected ordinates above 190, high and trivial zeros, and prime arithmetic. Positivity locates an actual signed move above the threshold in that support. For $`h=\log H`$, set $`u=h+3\log2/2`$. The *entire support* gives legal starts strictly after $`H`$ and completions strictly before $`400H`$ whenever

```math
x+R-\log(200)\lt u\lt x-R .
```

For each sign separately the open intervals cover each closed block in $`u`$, including its two endpoints. One verifier advances a frontier through the recorded rows with strict overlap. A second orders the lower endpoints and merges their open intervals. The 2,735 gain rows and 2,743 loss rows thus include every real starting height. [Covering table](evidence/sphere-cover.tsv) · [Full verification transcript](evidence/check.txt)

The first calculation uses 128-bit balls, the direct hyperbolic sphere ratio, sinc boxes and direct phase sines. The second uses 160-bit balls, a rationalized exponential formula for the sphere ratio, cosine differences for the boxes, the half-angle identity for the tent and a phase split at the left block endpoint. Both use the certified frequencies themselves, enclose all discarded terms and check every interval. The Lean module separately checks the high-frequency exponent, the exterior-cell contradiction and the preceding average and support deductions. The explicit formula, estimates and passage to the prime-only score have written proofs. [Numerical implementation](finite_sphere.py) · [Lean deductions](WaitingBound.lean)

## Time to two crossings with transcendental height ratio

Write $`f(x)=[S(2y)-S(y)]/\log2`$ with $`y=e^x/(2\sqrt2)`$, and put $`t=3947/2000`$. For every real $`1\le H\le e^{12500}`$ the following bounds hold independently of a global zero hypothesis. Both strict excursions $`f>t`$ and $`f<-t`$ finish before $`2e^{4.83}H<251H`$. There are two later starts with $`f=t`$ and $`f=-t`$ and transcendental quotient whose doublings both finish before $`2e^{7.38}H<3208H`$. If the order of the signs is prescribed, such a pair finishes before $`2e^{9.53}H<27534H`$. The two orders can use different pairs. At factor 400, the unordered and the two-order guarantees cover at least 99.35816% and 78.18624% of the interval in $`\log H`$, respectively. These percentages are lower bounds.

### An order test on four hits

On field 547/727 at most three hits of the named magnitude can lie in one class for the relation $`y\sim z\Longleftrightarrow y/z\in\overline{\mathbf Q}`$. The single-field proof throughout this range uses Baker, Roy and prime fourth powers and applies with the five primes omitted here. A **complete window** contains both strict excursions. By continuity it contains two different exact hits, one of each sign, strictly between its overshoots. Two chronologically disjoint complete windows supply four different hits. If all four opposite-sign comparisons had algebraic quotient the four hits would lie in one class. Three ordered overshoots of alternating signs also suffice. The two transitions each cross both levels, and the strict middle overshoot separates the two transitions. [Single-field class theorem](https://github.com/leecrossley/theta-host-split-inert-race/blob/algebraic-heights-v1.3.0/results/algebraic-heights/PROOF.md#algebraic-ratios-of-crossing-heights)

For a prescribed order use three successive complete windows, with hits $`P_i,N_i`$ of the indicated signs in window $`i`$. The forward comparisons $`P_1N_2,P_1N_3,P_2N_3`$ connect four different heights. If all three quotients were algebraic, the four would be equivalent. Interchange the signs for the opposite order. Four alternating overshoots also suffice. Their three successive transitions each contain both signed hits. The exact hits can be chosen in the *open* transitions, so they form three ordered groups even though the closed transitions share endpoints. This proves both order tests. [Lean deduction for the three comparisons](WaitingBound.lean)

The same argument gives an all-height consequence. Let $`C`$ be any fixed factor guaranteeing both strict excursions after every height at least one, as in the [fixed-factor theorem](#a-fixed-factor-without-a-global-zero-hypothesis). Starting each successive window at the last permitted start in the previous one gives factors $`C^2/2`$ and $`C^3/4`$ for the unordered and prescribed-order versions. More generally a prescribed word of $`m\ge1`$ signs can be realized in order by pairwise inequivalent exact heights in at most $`3m-2`$ consecutive complete windows. Accept the requested sign whenever its class is new. Every rejection uses one of at most two other distinct members of a previously accepted class. Before the last success there are at most $`2(m-1)`$ such rejections. This yields the factor $`2(C/2)^{3m-2}`$ at all heights, with the same unspecified $`C`$.

### Completing the local zero list

Merge 702 positive disjoint brackets from the [extended table](https://github.com/leecrossley/theta-host-split-inert-race/blob/asymptotic-waiting-factors-v1.0.0/results/asymptotic-waiting-factors/evidence/ordinate-brackets.tsv.xz) with the 6,000 used above. Strict Hardy signs prove existence in every additional bracket. For a primitive odd character $`\chi`$ evaluate its completed logarithmic derivative at $`s=3/2+iv`$ for $`v=0,3,\ldots,294`$. Its real part is the positive sum over strip zeros already used in the [finite zero count](#signed-positive-averages-from-a-finite-zero-verification). Evaluate the Euler sum for $`-L'/L`$ through $`M=65536`$. The entire discarded part is bounded by
```math
\sum_{n>M}\frac{\log n}{n^{3/2}}\le\int_M^\infty\frac{\log x}{x^{3/2}}\,dx
=\frac{2(\log M+2)}{\sqrt M}.
```
Subtract the positive ordinates of $`\chi`$ and the negatives of the positive ordinates of $`\overline\chi`$. For the two quartic characters at each prime this pairs different Conrey labels. In all 792 cases the rigorous **upper** residual is below $`4/13`$. An unaccounted line zero at distance at most $`3/2`$ contributes at least $`4/13`$. Platt’s unconditional theorem places all zeros at these heights on the critical line. Reflection across the family covers negative ordinates. Thus the list accounts for all 4,194 positive zeros through 294, counting multiplicity, and excludes a central zero. As before, the aggregate mass of unlisted line ordinates at the origin is below $`0.127`$, which remains valid on enlarging the list. [Exterior Euler and Hardy checks](deadline_spectrum.py) · [Extra brackets](evidence/extra-zeros.tsv) · [Platt, Theorem 7.1](https://arxiv.org/abs/1305.3087)

### The three error estimates

Use the positive spherical kernel and two boxes from the [previous cover](#cover-of-every-real-starting-height). Set $`c=8500`$, $`L=(T+82)/16984`$, $`d=1/k`$ and $`R=L+2d`$ for a row $`(T,\sigma,m,k)`$, centred at $`x=m/64`$. The earliest true log height used to bound the arithmetic is at least the integer $`\ell=\lfloor x-R-3\log2/2\rfloor\ge21`$. The entire signed average is a positive average of the original prime-only move.

First, all possible odd zeros away from the critical line have ordinate above 8590 in modulus by Platt. The previous strip estimate bounds their total contribution by $`30\cdot10^6\exp((T+7)/2+1/16-8492L)`$. Adding the odd trivial-zero bound at $`x-R>20`$ still gives less than $`0.001`$ on each support. Second, the omitted line terms cost at most $`50(0.127)\Psi_L(294i)\min(1,(294d)^{-2})`$. Charge the absolute transformed contributions of the selected ordinates above 270 separately. Evaluate the remaining 3,820 positive pairs with their signs. [Sphere and strip estimates](finite_sphere.py) · [New enclosures](ordered_averages.py)

Third, use the earlier level and signed-square arithmetic bounds in every stratum where they apply. For $`\ell\ge500`$ there is also the following estimate for the base zeros and prime-power remainder. On the strip $`|\Re z|\le1/2`$ one has $`|J(z)|\le66/(1+(\Im z)^2)`$. For $`|\Im z|\le2`$ use the total variation $`4/\log2`$ and support within $`[-3\log2/2,3\log2/2]`$. Above two use the explicit hyperbolic numerator, bound $`|(z-1/2)(z+1/2)|`$ below by $`(\Im z)^2`$ and use $`(1+t^2)/t^2\le5/4`$. Hence each base coefficient $`|J((\rho-1)/2)|/2`$ is at most $`132/(1+\gamma^2)`$. At $`s=2`$ each base factor has real-part mass less than 13. For $`1\le w=2-\beta\le2`$ one has $`(1+t^2)^{-1}\le2w/(w^2+t^2)`$. The total reciprocal-square mass over the eight base factors is below 208. [Transform](#the-transform-and-the-prime-powers) · [Mass inequality](WaitingBound.lean)

Put $`Q=6362704`$ and $`S(T)=80[\log(2QT)+1]/(T-2)`$. For $`T\ge8590`$ the explicit unit-bin count bounds the reciprocal-square mass above $`T`$ by $`S(T)`$. Take cutoffs $`T_j=e^{u_j}`$ with $`u_j=9.1,12,15,18,21,24,27,30,36`$. All exceed 8590. McCurley’s zero-free region for the *product* modulo $`Q`$ gives $`(1-\beta)/2\ge1/(20\log(QT_j))`$ through $`T_j`$, apart from at most one simple real zero. The primitive conductors divide $`Q`$, and the extra Euler factors in lifting characters only give zeros on the imaginary axis. The exceptional zero costs at most $`20/(e\ell)`$ because $`J(0)=0`$ and $`|J'(-v)|<40`$ for $`0\le v\le1/2`$. With $`E_j=\exp[-\ell/(20\log(QT_j))]`$, the base contribution is at most
```math
132\left(208E_0+\sum_{j=1}^8 S(T_{j-1})E_j+S(T_8)\right)+\frac{20}{e\ell}.
```
Add $`10^6(e^{-\ell/10}+e^{-\ell/3})+128e^{-(\ell-2)/2}/(1-e^{-2(\ell-2)})`$ for the higher powers, five omissions and trivial zeros. The verification uses the smallest of the valid arithmetic allowances in each stratum. Summing the signed modes on each of the 49,868 compact supports gives a strict excursion after every allowance. One implementation evaluates at the full phase. The other splits the phase and uses independent half-angle, cosine and exponential expressions for the weights. [McCurley’s theorem](https://www.sciencedirect.com/science/article/pii/0022314X84900891) · [Unit-bin count and exceptional-zero estimate](https://github.com/leecrossley/theta-host-split-inert-race/blob/theta-race-factor-400-v1.3.0/results/theta-race-factor-400/PROOF.md#comparing-to-the-actual-primes-throughout-the-search) · [Higher prime powers](#the-transform-and-the-prime-powers)

### Both endpoints of every waiting interval

For real $`H`$ put $`u=\log(2\sqrt2 H)`$ and $`i=\lceil100u\rceil`$. Round each certified support outward to integer endpoints $`l=\lfloor100(x-R)\rfloor`$ and $`r=\lceil100(x+R)\rceil`$. A strict averaged inequality supplies an overshoot in the open support by continuity. So a row with $`l\ge i`$ gives an overshoot after $`i/100`$ and before $`r/100`$, including when rounded supports touch. At the early cutoffs $`104\le i\le2176`$ the existing direct prime certificate supplies both strict signs within 220 grid units, since $`H_i=e^{i/100}/(2\sqrt2)\le10^9`$ and $`\log9<2.20`$. [Early prime-count certificate](https://github.com/leecrossley/theta-host-split-inert-race/blob/theta-race-factor-18-v1.2.0/results/theta-race-factor-18/PROOF.md#the-sharp-finite-deadline)

Let $`P(i),N(i)`$ be the earliest available finish for each sign and $`B(i)=\max(P(i),N(i))`$. For the unordered pair take the least of $`B^2(i),PNP(i),NPN(i)`$, where adjacent letters denote composition. For both orders take the least of $`B^3(i),PNPN(i),NPNP(i)`$. The order test just proved applies to every such schedule. A second implementation runs from the finish backwards, maximizing the latest permissible start. It uses a separate integer formula for support rounding, prefix maxima for each sign and the minimum for a complete window. For all 1,250,001 cutoffs it checks that the backward schedule fails one grid step before the forward finish and succeeds at that finish. This optimality is within the stored witnesses and these schedules. [Support table](evidence/ordered-supports.tsv.xz) · [Forward and backward verifiers](ordered_deadlines.py)

For the ordinary, unordered and two-order schedules, respectively, the largest value of $`j-i+1`$ over the grid is 483, 738 and 953. Since $`i-1<100u`$ and the final log centre is at most $`j/100`$, each elapsed log width is strictly smaller than the indicated value divided by 100. Multiplying its exponential by two accounts for the end of the last doubling. Finally, count only full centilog cells contained in the requested real interval. For 1,241,977 and 977,328 such cells the two stronger widths are at most $`5.29<\log200`$. Division by 1,250,000 gives the asserted lower fractions for completion by $`400H`$. [Reproducible counts](evidence/check.txt)
