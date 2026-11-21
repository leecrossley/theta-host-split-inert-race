# An explicit waiting interval for the 547/727 race

## Statement and prime convention

Let $K=\mathbf Q(\zeta_{16})^+$, $F=K(\sqrt{397669})$ and $E=K(\sqrt{-547},\sqrt{-727})$. For a prime ideal $v$ of $F$ let $\epsilon_v=1$ for split and $-1$ for inert in $E/F$. Use the primes outside $2,547,727,343057911569,343070636977$. With $a=\log 2$ and $B=2\sqrt2(\sqrt2-1)^2$ put

$$
W(u)=\begin{cases}2(1-1/u)&1<u\le2,\\4/u-1&2<u<4,\\0&\text{otherwise},\end{cases}
\qquad
S(y)=\frac{\sum_v-\epsilon_v\log(Nv)\,W(Nv/y)}{B\sqrt y}.
$$

For every real $H\ge1$ there are $y_+,y_->H$ such that

$$
2\max(y_+,y_-)< e^{162000}H^{3001/2999},\qquad
S(2y_+)-S(y_+)>1.9735a,\qquad S(2y_-)-S(y_-)<-1.9735a.
\tag{1}
$$

For $\log H\ge500000000$, $H^{1.001}$ also bounds both endpoints. For $1\le H\le10^6$ a direct certificate gives $2\max(y_+,y_-)\le18H$, improving to $(6851/657)H$ for $100\le H\le10^6$. For $0<H<1$ the two moves starting at 5 and 9 suffice. The statements concern changes in the lead and allow either ordering of the two moves.

The prior field construction proves the degrees 8 and 16, the unramified quadratic extension at finite places, and the following identity for an admissible datum and every input parameter family. The product of the normalized first-step theta trace jumps in the fibre above a good $v$ equals $\epsilon_v$. The two large exclusions are the type primes of the level datum. Indeed membership in its exceptional set is obtained by restriction of a place over one of the two selected primes in the quadratic subfield. The other three exclusions fix the convention at rational ramification. Every omitted norm is at most $343070636977^8<e^{216}$. We can therefore use the full prime sum in the analytic proof wherever the smallest norm in its windows exceeds $e^{216}$. [Field and admissibility](https://github.com/leecrossley/theta-race-short-shifts/blob/483c31f2df2132ed2165620ea55e20e969228bfb/SecondField.lean) · [Fibre trace theorem](https://github.com/leecrossley/theta-race-short-shifts/blob/483c31f2df2132ed2165620ea55e20e969228bfb/RaceTrace.lean#L83) · [Type primes above their rational primes](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Formalism/LocalToGlobal/UnitaryNorms/CharHost.lean#L44)

## The transform and the prime powers

Write $\chi_p=(\cdot/p)$ and let $\xi$ be the even quartic character of conductor 16 with $\xi(5)=i$. At primitive conductors the eight factors of $L=\zeta_E/\zeta_F$ are $\chi_p\xi^j$, $p=547,727$, $0\le j\le3$. The eight factors of $\zeta_F$ are $\xi^j$ and $\chi_{547}\chi_{727}\xi^j$, including ordinary zeta for the trivial character. These are the characters odd and even, respectively, on complex conjugation $(-1,-1,1)$ in $C_2\times C_2\times C_4$. Abelian factorization includes the primitive Euler factors at rational ramification.

In this section use the full score, put $Q(x)=S(e^x/2)$ and $f(x)=[Q(x+a/2)-Q(x-a/2)]/a$. Then $f(x)=[S(2y)-S(y)]/a$ for $y=e^x/(2\sqrt2)$. The bilateral Laplace transform of $f$ on $\Re z>1/2$ is $J(z)D(1/2+z)$, where

$$
D(s)=-\sum_v\epsilon_v\log(Nv)(Nv)^{-s},\quad
G(z)=\frac{2\sqrt2\cosh(az)-3}{2(\sqrt2-1)^2(z^2-1/4)},\quad
J(z)=G(z)\frac{2\sinh(az/2)}a.
$$

The singularities of $G$ at $\pm1/2$ are removable. Integrate the two branches of $W$ after setting $u=2Nv e^{-x}$ to obtain $G$. It is the transform of a nonnegative compact kernel on $[-a,a]$ with mass $G(0)=2$. Thus $J$ is the transform of a signed kernel of total variation at most $4/a<6$ and support in $[-3a/2,3a/2]\subset(-2,2)$. In particular $J(0)=0$. For $|\Re z|\le1/2$ we shall use

$$
|J(z)|\le \frac{1000}{1+(\Im z)^2}.
\tag{2}
$$

For $|t|=|\Im z|\le1$ use the total variation bound. For $|t|\ge1$ use $|(z-1/2)(z+1/2)|\ge t^2$, $|\cosh(az)|\le\cosh(a/2)$ and $|\sinh(az/2)|\le\cosh(a/4)$.

Coefficient comparison over prime ideals gives

$$
D(s)=\frac{L'}{L}(s)-\frac{\zeta_F'}{\zeta_F}(2s)+D_R(s).
\tag{3}
$$

The prime squares cancel since $\epsilon_v^2=1$. The remaining series $D_R=\sum_n c_n n^{-s}$ has powers at least three. At $\sigma=2/5$ its total absolute mass is at most

$$
\frac{2}{1-2^{-2/5}}\sum_v\frac{\log Nv}{(Nv)^{6/5}}<10000.
$$

Indeed $\sum_{v\mid p}\log Nv\le8\log p$, and $\sum_{n\ge2}\log(n)n^{-6/5}<110$ by $\log n\le10n^{1/10}$ and integral comparison. The factor 2 covers the two remaining prime-power series. In physical coordinates the part $f_R(x)$ uses $n\le e^{x+2}$. Since $0\le W\le1$, the two square-root normalizations and $1/a$ cost less than $40e^{-x/2}$. Consequently

$$
|f_R(x)|<40e^{-x/2}\,10000\, e^{(2/5)(x+2)}<10^6e^{-x/10}.
\tag{4}
$$

## An isolated zero and its coefficient

Set $d=1/1000$. There is a unique simple zero of $L$ at $1/2+i\gamma$ with

$$
3287178298086/10^{12}<\gamma<3287178298087/10^{12}.
$$

On the critical line the open $d$-neighborhood of $\gamma$ contains this zero alone, and the $d$-neighborhood of zero contains none. The conductor-547 factor has opposite Hardy signs at the endpoints. Throughout the larger interval around $\gamma$ the real part of its complex derivative exceeds 1. Thus the imaginary part of $L(s,\chi_{547})$ is strictly increasing with the ordinate there, and its zero is simple. Closed interval enclosures exclude zeros of every other odd factor and check the central intervals. Platt's finite-height theorem places every nontrivial zero of these eight factors through height 8500 on the critical line, since $11632\cdot8500<10^8$. The family is closed under conjugation, giving the same coverage at negative heights. [Interval verification](verify.py) · [Independent crossing](independent.py) · [Platt, Theorem 7.1](https://arxiv.org/pdf/1305.3087#page=14)

At the selected zero, $J(i\gamma)=iV$ with

$$
V=\frac{3-2\sqrt2\cos(a\gamma)}{2(\sqrt2-1)^2(1/4+\gamma^2)}
       \frac{2\sin(a\gamma/2)}a,\qquad 3.342843<V<3.342844.
$$

The square term in (3) requires a **complex disk** test. For $\omega=0,\gamma$ and for any nontrivial zero $\rho$ of any factor of $\zeta_F$ we have

$$
\left|(\rho-1)/2-i\omega\right|\ge d .
\tag{5}
$$

Here is the certificate. Use disks of radius $r=2001/10^6$ about $1$ and $1+i\,6574356596173/10^{12}$. They contain both closed disks of radius $2d$ at the desired centres. For a nonprincipal character write $C(u)=\sum_{n\le u}\chi(n)$. A complete period supplies $|C|\le M$ by taking the maximum of the squared norms of Gaussian-integer partial sums. With $N=2^{18}$, partial summation gives errors at each centre $s$ of at most $M(1+|s|)/N$ for $L(s)$ and $M[\log N+1+|s|(\log N+1)]/N$ for $L'(s)$. The truncated sums and their derivatives are evaluated with complex balls.

For the Taylor remainder on the whole disk, $\Re s\ge.99$, $|s|<8$ and $|C(u)|\le\min(u,q)$. Put $I_k=\int_1^\infty\min(u,q) u^{-1-\Re s}(\log u)^k\,du$. Splitting at $q$ and using $q\le e^{16}$ gives

$$
I_k\le e^{.16}\left[\frac{16^{k+1}}{k+1}+
             \int_0^\infty e^{-.99v}(16+v)^k\,dv\right].
$$

Twice differentiating $s\int_1^\infty C(u)u^{-s-1}\,du$ gives $|L''|\le2I_1+8I_2<30000$. Subtract the two truncation errors with their proper weights, $r|L'|$ and the quadratic remainder $15000r^2$, as implemented in the verifier. Every resulting lower bound exceeds $1/100$. For the zeta factor use $\zeta(s)=s/(s-1)-s\int_1^\infty\{u\}u^{-s-1}\,du$. At 1 the regularized function $(s-1)\zeta(s)$ stays away from zero by $1-r-(1+r)r/.99>0.9$. At the other centre interval evaluation and the integral's derivative give $|\zeta(1+2i\gamma_0)|-12r>1/100$, where $2\gamma_0=6574356596173/10^{12}$. This proves (5). The simple pole of zeta gives a term at $z=0$ annihilated by $J(0)$. [Prime-square disk verifier](square_factors.py) · [Numerical receipt](evidence/check.txt)

## Total mass of the omitted zeros

For all sixteen primitive factors together, including multiplicities, we use

$$
\sum_{\rho}\frac1{1+(\Im\rho)^2}<1000 .
\tag{6}
$$

For one factor, $\Re(2-\rho)^{-1}\ge 1/[4(1+(\Im\rho)^2)]$. For a primitive nonprincipal Dirichlet character of parity $k$ the sum of these real parts is

$$
\tfrac12\log(q/\pi)+\tfrac12\psi((2+k)/2)+\Re(L'/L)(2,\chi).
$$

This follows by subtracting the Hadamard logarithmic derivatives of the completed function at $2$ and $-1$ and using the functional equation. The zeros are paired by $\rho\mapsto1-\overline\rho$ for a fixed character. The genus-one subtraction makes the difference absolutely convergent. For zeta, completing with $s(s-1)$ adds $1/2+1$. Here $q\le6362704$, $\psi(1),\psi(3/2)<1$, and the Euler product gives $|(L'/L)(2)|\le\sum_{n\ge2}\log(n)n^{-2}<2$. Each sum of real parts is less than 13, proving (6) since $16\cdot4\cdot13<1000$.

## A compact average at each starting height

Put $h=\log H$, $X=(h+81000)/(1-1/3000)$ and $n=\lceil X\rceil$. Take $K_X$ to be the probability law of $x=X$ plus the sum of 40 independent uniforms on $[-2000,2000]$ and $n$ independent uniforms on $[-1/3000,1/3000]$. Its support is $[X-R,X+R]$ for $R=80000+n/3000$. Its centred exponential moment is

$$
\Phi(w)=\left(\frac{\sinh(2000w)}{2000w}\right)^{40}
        \left(\frac{\sinh(w/3000)}{w/3000}\right)^n .
$$

Each factor at zero is interpreted by continuity. Direct algebra gives

$$
X-R>h+999,\qquad X+R<162000+(3001/2999)h.
\tag{7}
$$

For $\omega=0,\gamma$ form $M_\omega=\int f(x)e^{-i\omega x}\,dK_X(x)$ and put $T_\omega(z)=J(z)e^{(z-i\omega)X}\Phi(z-i\omega)$. Logarithmic derivatives in (3) contribute the following terms to $M_\omega$.

| Source | Term |
|---|---|
| Nontrivial or trivial zero $\rho$ of an odd factor | $T_\omega(\rho-1/2)$ |
| Nontrivial or trivial zero $\rho$ of a base-field factor | $-T_\omega((\rho-1)/2)/2$ |
| Pole of ordinary zeta at 1 | $T_\omega(0)/2=0$ |

To justify the formula, invert the absolutely convergent prime-power transforms on $\Re z=1$ and integrate against the compact measure. The fixed window has $O((1+|\Im z|)^{-2})$ decay there. For each primitive factor expand its logarithmic derivative by the genus-one Hadamard product for the completed function and the partial fractions for its gamma factor. Subtraction at a fixed real $s_0>1$ gives terms $(s-\rho)^{-1}-(s_0-\rho)^{-1}$ for each zero, and the opposite term at a pole, up to a constant. The gamma terms give the trivial zeros $-k,-k-2,\ldots$ for nonprincipal factors, and $-2,-4,\ldots$ for zeta. The elementary factors of completed zeta remove the spurious zero at 0 and give its pole at 1.

For clarity, the inverse integral of $T_\omega(z)/(z-w)$ on this line equals $T_\omega(w)$ for $\Re w<1$, and that of $T_\omega(z)$ equals zero. This follows by closing to the left. On $\Re z=-U$ the exponential decay from $X>R+2$ dominates $J$ and the boxes. The frequency integrals are uniformly controlled there by retaining two long-box denominators, giving a bound by a constant (for fixed $X$) times $e^{-U(X-R-2)}\int_\mathbf R (1+|t-\omega|)^{-2}\,dt$ for $U\ge1$. Horizontal integrals vanish for each fixed $U$ by the same box decay. Alternatively the reciprocal kernel is the transform of $e^{wu}\mathbf1_{u>0}$ and the entire averaging window stays on positive log heights. On the original line each subtracted partial fraction, summed absolutely, has at most polynomial growth in $|\Im z|$: split $|\Im\rho|$ at $2|\Im s|+4$, use the distance of $\Re s$ from the strip on the near terms, and apply (6) on both parts. The trivial zeros give an absolutely summable subtracted series there as well. This legitimizes integrating term by term. The square factor has the additional $1/2$ from its change of scale.

For a box of half-width $s$,

$$
\left|\frac{\sinh(sw)}{sw}\right|
\le e^{s|\Re w|}\min\{1,1/(s|w|)\}.
$$

This follows from the box integral for the first bound and $|\sinh(u+iv)|\le\cosh u$ for the second. If $\Re z\le0$ and $|z-i\omega|\ge d$, (7) and the 40 long boxes give $|e^{(z-i\omega)X}\Phi(z-i\omega)|\le2^{-40}$. This covers all nontrivial base zeros by (5), and all odd-factor zeros on or left of the line apart from the selected term for $\omega=\gamma$. If an odd-factor zero gives $\Re z>0$, its ordinate has magnitude above 8500 by Platt. Then $|z-i\omega|/3000>(8500-4)/3000>e$. The short boxes give a bound $\exp((X+R)/2-n)\le\exp(-(X-R)/2)$ on its contribution before multiplying by $J$. This deals with every higher zero anywhere in the critical strip.

From (2) the kernel weight is at most $1000/(1+t^2)$ for an odd-factor zero of ordinate $t$ and $4000/(1+t^2)$ for a base zero, before its factor $1/2$. Thus (6) bounds the total error from all unselected nontrivial zeros by $4\cdot10^6(2^{-40}+e^{-999/2})<.001$. Each sequence of trivial zeros starts at $\Re z\le-1/2$ and has spacing at least 1. Use $|J(z)|\le 6e^{2|\Re z|}$ and the box growth bound to get a total below $128e^{-997/2}/(1-e^{-997})<.001$. Finally (4) costs less than $.001$ on the entire support. The omitted finite primes cannot enter any of these windows since the smallest possible norm exceeds $e^{X-R-2}>e^{216}$. We have proved

$$
|M_0|<.01,\qquad |M_\gamma-iV|<.01 .
\tag{8}
$$

## Two witnesses

The tests $P_+(x)=1-\sin(\gamma x)$ and $P_-(x)=1+\sin(\gamma x)$ are nonnegative. Each has mass at most $1+2^{-40}$ under $K_X$ because $|\Phi(i\gamma)|\le 2^{-40}$. Their $f$-integrals are $M_0+\Im M_\gamma$ and $M_0-\Im M_\gamma$, respectively. By (8) and $V>3.342843$, the first exceeds $1.9735$ times its mass and the second is below $-1.9735$ times its mass. There is consequently a point of each sign and magnitude in the compact interval. No common phase or independence hypothesis on the zeros enters the argument.

For each point take $y=e^x/(2\sqrt2)$. The lower inequality of (7) gives $y>H$. Its upper inequality bounds $2y=e^x/\sqrt2$, proving (1). The inequality $(1.001-3001/2999)\cdot500000000>162000$ proves the stated large-height simplification. [Lean proofs of the interval and positive-test implications](WaitingBound.lean)

## Finite interval and independent checks

For a rational prime $p$ outside the exclusions put $u=\chi_{547}(p)$ and $v=\chi_{727}(p)$. The residue degree $f$ in $F$ is the least $j=1,2,4$ with $p^j\equiv\pm1\bmod16$ and $(uv)^j=1$. There are $8/f$ ideals of norm $p^f$, each of sign $u^f$. The finite verifier sieves every relevant rational prime, forms coefficients $-8u^f\log p$ and integrates the two branches of the weight using prefix sums of the coefficients and the coefficients divided by the norms. It certifies moves on an integer grid starting at 1 with $y_{\rm next}=\lceil9y/8\rceil$. For real $H$ the next strictly later witness of each sign is constant between consecutive certified starting points. The endpoint ratio within each half-open interval therefore attains its maximum at the left endpoint. Checking all such endpoints and the upper endpoint covers the entire real range, including equality at a witness. The independent enumerator directly sums $W$ separately for the four example moves at the two extreme starting heights. [Continuous coverage](finite_windows.py) · [Direct enumeration](independent.py)

The independent zero check expands the conductor-547 Dirichlet series by residue classes after one period. In the Hurwitz tail centred at $3/2$, the $k$th term for $k\ge1$ is
$q^{-s}(-1)^k(s)_k\zeta(s+k,3/2)\sum_{b=1}^q\chi(b)(b/q-1/2)^k/k!$.
The constant term cancels. For $|s|<4$, bound $|(s)_k/k!|\le {k+3\choose3}$, the Hurwitz factor by $4(3/2)^{-k}$ and the inner sum by $q/2^k$. After $k=70$ the tail ratio of the bounds is below $1/2$. This gives the implemented error and a second pair of opposite signs of the real completed quadratic function. [Independent calculation](independent.py)

The analytic estimates and explicit formula above are written proofs. The release replay checks the complex balls and direct counts, compiles the five local Lean theorems and independently rechecks the compiled module with the kernel. The pinned field and trace release supplies its own replay. [Prior host receipt](https://github.com/leecrossley/theta-race-short-shifts/blob/483c31f2df2132ed2165620ea55e20e969228bfb/evidence/check.txt) · [This release's receipt](evidence/check.txt)
