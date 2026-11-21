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

The prior field construction proves the degrees 8 and 16, the unramified quadratic extension at finite places, and the following identity for an admissible datum and every input parameter family. The product of the normalized first-step theta trace jumps in the fibre above a good $`v`$ equals $`\epsilon_v`$. The two large exclusions are the type primes of the level datum. Indeed membership in its exceptional set is obtained by restriction of a place over one of the two selected primes in the quadratic subfield. The other three exclusions fix the convention at rational ramification. Every omitted norm is at most $`343070636977^8\lt e^{216}`$. We can therefore use the full prime sum in the analytic proof wherever the smallest norm in its windows exceeds $`e^{216}`$. [Field and admissibility](https://github.com/leecrossley/theta-race-short-shifts/blob/483c31f2df2132ed2165620ea55e20e969228bfb/SecondField.lean) · [Fibre trace theorem](https://github.com/leecrossley/theta-race-short-shifts/blob/483c31f2df2132ed2165620ea55e20e969228bfb/RaceTrace.lean#L83) · [Type primes above their rational primes](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Formalism/LocalToGlobal/UnitaryNorms/CharHost.lean#L44)

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

The analytic estimates and explicit formula above are written proofs. The release replay checks the complex balls and direct counts, compiles the ten local Lean theorems and independently rechecks the compiled module with the kernel. The pinned field and trace release supplies its own replay. [Prior host receipt](https://github.com/leecrossley/theta-race-short-shifts/blob/483c31f2df2132ed2165620ea55e20e969228bfb/evidence/check.txt) · [This release's receipt](evidence/check.txt)

## Both moves through exp(2000)

**Theorem.** For every real $`1\le H\le e^{1000}`$ both qualifying doublings finish before $`256H`$. For $`1\le H\le e^{2000}`$ they finish before $`400H`$. The magnitude in each direction is strictly greater than $`1.9735a`$, with $`a=\log2`$. Both starts must exceed $`H`$. The published direct sieve gives $`18H`$ through $`10^9`$, including the endpoint. Its proof uses actual prime ideals. It therefore suffices here to cover $`10^9\le H`$ in the two new ranges. [Direct sieve and continuous-height coverage](https://github.com/leecrossley/theta-race-factor-18/blob/31ea5d495ca29f7f17869d1a3f321000266e4867/PROOF.md)

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

**Theorem.** For every real $`H\ge1`$, both signed dyadic changes of magnitude greater than $`3.32\log2`$ finish strictly before $`e^{161020}H^{17001/16999}`$. The preceding explicit all-height interval used exponent $`3001/2999`$. The new estimate uses the same isolated zero, the same complex disk exclusions for the base-field zeros, and the same 40 long boxes. Here a positive Logan kernel controls the unknown high zeros. [Earlier all-height proof](https://github.com/leecrossley/theta-race-waiting-bound/releases/tag/v1.0.1) · [Büthe's use of the Logan kernel in finite-height prime estimates](https://arxiv.org/abs/1410.7015)

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
