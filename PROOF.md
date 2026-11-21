# Proof

## The theorem

Put $K=\mathbf Q(\zeta_{16})^+$, $E=K(\sqrt{-547},\sqrt{-727})$ and $F=K(\sqrt{547\cdot727})$. Let $A_b,B_b,S_b$ and $W_b$ be as in the [statement](README.md). All logarithms are natural. Use log geometric centre $x$ for the window whose lower edge is $e^x/b$, and set $Q_b(x)=S_b(e^x/b)$. For $0< t\le\log 2$ put

$$
D_{b,t}(x)=\frac{Q_b(x+t/2)-Q_b(x-t/2)}{t}.
$$

Set $D_{b,0}(x)=Q'_{b,+}(x)$, the right derivative. The function is piecewise smooth on bounded intervals, and this choice agrees almost everywhere with the ordinary derivative. Let $\mu$ be any compactly supported Borel probability measure on $(1,\infty)\times[0,\log 2]$, fixed independently of $x$. Then

$$
\operatorname*{ess\,liminf}_{x\to\infty}\int D_{b,t}(x)\,d\mu(b,t)<-1.9735,\qquad
\operatorname*{ess\,limsup}_{x\to\infty}\int D_{b,t}(x)\,d\mu(b,t)>1.9735 .
\tag{1}
$$

In particular both strict inequalities are attained on sets of positive Lebesgue measure arbitrarily far out. With a point mass at $(b,\log r)$ this gives the two differences in the README by putting $y=e^{x-t/2}/b$. The point mass at $(b,0)$ gives the derivative statement. Each measure has its own excursions. For each fixed $b>1$ we also have

$$
\liminf_{y\to\infty} S_b(y)<-5.6921,\qquad
\limsup_{y\to\infty} S_b(y)>7.6921.
\tag{2}
$$

## The field and its trace

The quadratic discriminants $-547$ and $-727$ are coprime and odd. The even quartic character of conductor 16 and the two quadratic characters generate independent character groups by the Chinese remainder theorem. Thus $\operatorname{Gal}(E/\mathbf Q)=C_2\times C_2\times C_4$, complex conjugation is $c=(-1,-1,1)$ and its fixed field is $F$. In particular $[E:\mathbf Q]=16$ and $[F:\mathbf Q]=8$. The prime 547 ramifies in $F$ and is unramified in the earlier 571/599 field, so the two base fields are distinct.

The field uses layer $e=2$, with the upstream rank datum `R0`. The two type primes are $343057911569$ and $343070636977$. They are 1 modulo $16\cdot547\cdot727$. The two distinguished split primes are 991 and 1231. For the norm-prime condition write $M=16\cdot547\cdot727\,(343057911569\cdot343070636977)^2$. The element $1+41M\sqrt{-727}$ lies in the prescribed quadratic subfield, is 1 modulo $M$, and its norm from $\mathbf Q(\sqrt{-727})$ to $\mathbf Q$ is the prime $1+727(41M)^2$. Lean verifies Lucas primality certificates, disjointness of the two exceptional sets, both splitting conditions, the norm and congruence, and constructs `exampleRecordOfData` on these inputs. The constant weights at infinity provide `MUQ`. The resulting character frame gives one admissible datum for all the good primes of this field. [Auxiliary primes](HostPrimes.lean) · [Record, degrees and frame](SecondField.lean) · [Upstream constructor](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Formalism/LocalRingCardinality/SpanPowerMembership/ExampleClassGroupCardPos.lean#L82)

For each $w$ of $E$ above $v$ of $F$, the first parameter step gives a normalized jump $J(w)$. Subtract the previous first trace multiplied by the ratio of the source and target characters evaluated at $w$, and divide by the first trace of the added untwisted segment. The singleton trace lemma identifies this jump with the source-character value. The source character at this step is $\mu_1\eta^{-1}$. On restriction to the base, $\mu_1$ has exponent 0 and $\eta$ has exponent 1 for the quadratic character of $E/F$. The extension is unramified at finite places. Taking the product across a good fibre therefore gives $\prod_{w\mid v}J(w)=\epsilon(v)$, for every input parameter family. The product contains two factors for split $v$ and one for inert $v$. The count uses each base ideal once. The input family cancels algebraically. Identifying an individual cuspidal eigensystem with these parameters is a separate question. [Trace extraction](https://github.com/leecrossley/theta-prime-trace/blob/7a7909d2dc7bb9d8d0ec7b2dc8c78e70f93782d7/ThetaPrimeTrace.lean#L68) · [Generic fibre lemma](FiberTrace.lean) · [Concrete theorem](RaceTrace.lean)

## Euler factors and the kernel

Omit the primes above $2,547,727$ and any fixed finite exceptional set. Write $\chi_p(n)=(n/p)$, and let $\xi$ be the primitive even quartic character modulo 16 with $\xi(5)=i$. The eight characters odd on $c$ are $\chi_p\xi^j$ for $p=547,727$ and $0\le j\le3$, taken at primitive conductor. Consequently
$$
L(s)=\zeta_E(s)/\zeta_F(s)=\prod_{p\in\{547,727\}}\prod_{j=0}^3L(s,\chi_p\xi^j).
$$

For a good rational prime of Frobenius $g$, let $f$ be the order of $g$ modulo $\langle c\rangle$. There are $8/f$ base ideals of norm $p^f$ and their signs agree. The sign is $+1$ for $g^f=1$ and $-1$ for $g^f=c$. This also describes the Euler product of $L$ over $F$. The independent enumerator checks the character sums on every Frobenius class. For all 77 good rational primes through 400 it also factors $X^4-4X^2+2$, $Y^2-397669$ and $Z^2+547$ over the residue field. Tensoring fields of degrees $j,k$ gives $\gcd(j,k)$ copies of degree $\operatorname{lcm}(j,k)$. The degree of the last quadratic factor then determines whether $Z^2+547$ splits in each base residue field. This independently reproduces the norms, multiplicities and signs used in the two finite examples. [Euler factors and direct enumeration](independent.py)

Define $D(s)=-\sum_v\epsilon(v)\log Nv\,(Nv)^{-s}$ initially for $\Re s>1$. Comparing prime powers gives
$$
D(s)=\frac{L'}{L}(s)-\frac{\zeta_F'}{\zeta_F}(2s)+R(s),
\qquad \Re s>1/3 ,
\tag{3}
$$
as a meromorphic continuation. The $L'/L$ sum contributes the negative of each signed prime-power term. Since $\epsilon(v)^2=1$, the second term restores the squares. The remaining signed terms have exponent at least three. The extra terms from the zeta factor have exponent at least four. Both remainders converge absolutely on this half-plane. Omitted finite Euler factors are regular there.

The real segment $[1/2,1]$ of every Dirichlet factor is zero-free by the interval cover below. For $\Re s>1$ the Euler products are nonzero. At $s=1/2$, (3) has residue $1/2$ from the pole of $\zeta_F(2s)$. At $s=1/2+i\omega$, $\omega\ne0$, the prime ideal theorem gives $\zeta_F(1+2i\omega)\ne0$. Hence $D$ has residue $m_\omega$, the total zero multiplicity of the eight factors there. A regular point has $m_\omega=0$. The eight factors are stable under complex conjugation and $m_{-\omega}=m_\omega$.

Direct integration on the two branches of the window gives the entire Mellin transform
$$
H_b(s)=\int_1^{b^2}W_b(u)u^{s-1}\,du
       =\frac{(b^s-1)(b^s-b)}{s(s-1)}
$$
with removable values at 0 and 1. The bilateral Laplace transform of $Q_b$ is $G_b(z)D(1/2+z)$, where $G_b(z)=b^{-z}H_b(1/2+z)/B_b$. In particular $G_b(0)=2$. On the imaginary axis the centring cancels the phase, giving
$$
G_b(i\gamma)=
\frac{b+1-2\sqrt b\cos(\gamma\log b)}
 {2(\sqrt b-1)^2(1/4+\gamma^2)}
 \ \ge\ \frac{1}{2(1/4+\gamma^2)}.
\tag{4}
$$

## What a short shift does

Write $U_t(z)=2\sinh(zt/2)/t$ for $t>0$ and $U_0(z)=z$. Translation and, at $t=0$, integration by parts give the following bilateral Laplace transform of $\int D_{b,t}\,d\mu$.
$$
D(1/2+z)\int G_b(z)U_t(z)\,d\mu(b,t).
\tag{5}
$$

Here are the convergence details for arbitrary compactly supported $\mu$. Its widths lie in $[b_0,b_1]$ with $1<b_0\le b_1<\infty$. The scores vanish sufficiently far to the left uniformly in those widths. On bounded intervals of $x$, only finitely many prime ideals can contribute. Each summand is continuous and piecewise smooth with uniformly bounded one-sided log derivatives. The difference quotient is bounded by the supremum of those derivatives on the translated interval, uniformly down to $t=0$. For large $x$, all the norms are $O(e^x)$. Absolute convergence of the Dedekind ideal series gives a bound $O_\varepsilon(e^{(1/2+\varepsilon)x})$ for the scores and their one-sided slopes, absorbing the logarithm into $\varepsilon$. Thus Fubini, translation and integration by parts are justified first for $\Re z>1/2$. The kernel integral in (5) is entire by compactness, since $U_t(z)$ is continuous at $t=0$ locally uniformly in $z$. Taking the transform on $x\ge0$ changes (5) by an entire function.

Every $U_t(0)$ vanishes, so (5) loses the square-prime pole. At $z=i\gamma$ it has residue $i\,m_\gamma V_\gamma$, where
$$
V_\gamma=\int G_b(i\gamma)K_t(\gamma)\,d\mu(b,t),\qquad
K_t(\gamma)=
\begin{cases}
2\sin(\gamma t/2)/t,&t>0,\\
\gamma,&t=0 .
\end{cases}
$$
At $-i\gamma$ the residue is the complex conjugate. All five selected positive ordinates have $0<\gamma<0.58$, so all these multipliers are positive on $0\le t\le\log2$.

## Five ordinates and the test

The table gives open brackets in units of $10^{-9}$. Each row certifies at least one zero. Two rows belong to the same primitive character. Choose any one zero ordinate $\gamma_j$ in each interval.

| Primitive character | Left numerator | Right numerator |
|---|---:|---:|
| 11632.7269 | 48536362 | 48536363 |
| 8752.1093 | 248572181 | 248572182 |
| 4376.1093 | 434216664 | 434216665 |
| 547.546 | 446011529 | 446011530 |
| 11632.7269 | 579615600 | 579615601 |

For each signed sum of the five brackets with coefficients $0,\pm1$ and at least two nonzero coefficients, the resulting interval excludes zero, and its absolute value is zero-free for all eight factors on the critical line. There are $(3^5-1-10)/2=116$ combinations up to sign and $928$ individual factor tests. The brackets are disjoint. Thus in $\prod_j(1\mp\sin(\gamma_j x))$ the constant coefficient is 1 and the only other terms meeting a pole are the positive and negative singleton frequencies. A collision among two nonresonant frequencies contributes zero.

For completeness, let $f$ denote the average slope in (1). To bound its upper extreme it suffices to treat the case that $f$ is essentially bounded above on a tail. A constant minus $f$ is then nonnegative there. Its Laplace transform converges absolutely for every positive real argument. Indeed (3) and (5) give a continuation regular at every positive real point. Landau's theorem forces the abscissa of a nonnegative Laplace transform to be singular if finite. One proof is to expand just to the right of the abscissa and use monotone convergence of the exponential series at a point to its left if continuation were regular there. The added constant has its only pole at zero. The lower extreme follows by using $f$ minus an eventual lower bound.

With that one-sided convergence available, multiply $f(x)$ by the nonnegative products
$$
P_+(x)=\prod_{j=1}^5(1-\sin\gamma_jx),
\qquad
P_-(x)=\prod_{j=1}^5(1+\sin\gamma_jx).
$$
For each frequency $\omega$, the limit of $\sigma\int_0^\infty f(x)e^{-(\sigma+i\omega)x}\,dx$ as $\sigma\downarrow0$ is the residue at $i\omega$, with zero limit at a regular point. The two conjugate terms from $-\sin\gamma_jx$ select $+m_{\gamma_j}V_{\gamma_j}$, and $+\sin$ selects its negative. Each product alone has Abel mean 1. The Abel mean of $fP_\pm$ is therefore $\pm\sum_j m_{\gamma_j}V_{\gamma_j}$. An essential tail bound for $f$ passes to the corresponding mean because the weights are nonnegative and a fixed initial segment has zero limiting weight. This proves the desired two-sided inequalities once the sum is bounded below.

If $l_j,u_j$ are the endpoints after division by $10^9$, put $c_0=69315/100000>\log 2$. The Taylor inequality $\sin v\ge v-v^3/6$ for $v\ge0$ and (4) give, including $t=0$,
$$
G_b(i\gamma_j)K_t(\gamma_j)\ \ge\
\frac{l_j\,[1-(u_jc_0)^2/24]}{2(1/4+u_j^2)} .
$$
The sum of these five rational lower bounds is strictly greater than $19735/10000$. The inequality survives integration against $\mu$, and each $m_{\gamma_j}\ge1$. This proves (1). The interval-to-multiplier inequality for all real $t$ in the stated range and the exact final rational comparison are also checked in Lean. [Uniform bound](ShiftBound.lean)

For (2) use $Q_b$ itself and the nonnegative products $\prod_j(1\pm\cos\gamma_jx)$. Its constant residue is $G_b(0)/2=1$. Each singleton contributes $\pm m_{\gamma_j}G_b(i\gamma_j)$. The same frequency exclusions and (4) apply. The exact sum $\sum_j1/[2(1/4+u_j^2)]$ exceeds $66921/10000$. This gives (2).

## Numerical verification

The first calculation evaluates the primitive Hardy functions with Arb. Opposite endpoint signs give the five crossings by continuity. A ball evaluation on every full mixed interval excludes any zero there. Rational arithmetic also checks the two constants. [Hardy calculation](verify.py) · [FLINT documentation](https://flintlib.org/doc/acb_dirichlet.html)

The second calculation builds the characters using Legendre symbols and the quartic character table. It evaluates the Dirichlet series by residue classes and a Taylor expansion of the Hurwitz zeta tail, without calling the Hardy or Dirichlet L-function evaluator. Split off the first $q$ terms. At $\Re s=1/2$ the tail is
$$
q^{-s}\sum_{k\ge1}\frac{(-1)^k(s)_k}{k!}\zeta(s+k,3/2)
 \sum_{a=1}^q\chi(a)(a/q-1/2)^k.
$$
The constant term vanishes by periodicity. For $|s|<2$ the absolute value of $(s)_k/k!$ is at most $k+1$. Integral comparison bounds the Hurwitz factor by $4(3/2)^{-k}$ and the inner sum by $q2^{-k}$. After $k=54$, the discarded tail is at most $\sqrt q(2k+5)/3^k$ with $k=54$. Every evaluated ordinate satisfies the stated $|s|$ bound. The functional equation makes $\Lambda(1/2+i\gamma,\chi)\overline{\Lambda(1/2,\chi)}$ real for each character, giving an independent sign-crossing test. For a mixed interval, the second verifier subtracts a rigorous variation bound from the midpoint distance to zero. It also checks the sine bound directly at the maximal shift, using that $\sin v/v$ decreases and $\sin v$ increases on the range in question. Both implementations use Arb ball arithmetic. [Residue-class verifier and direct examples](independent.py)

For the real-axis cover, one period of each character supplies $|C(u)|\le\min(u,M)$ for $C(u)=\sum_{n\le u}\chi(n)$. The entries are checked against the Conrey labels. From $L(s)=s\int_1^\infty C(u)u^{-s-1}\,du$ and its derivative, the two bounding integrals at $\Re s\ge1/2$ give $4\sqrt M$ and $4\sqrt M\log M+4$, respectively, after splitting at $M$. Thus for real $1/2\le s\le1$, $|L'(s)|\le4\sqrt M+4\sqrt M\log M+4$. Multiply by the interval radius and compare with the magnitude at the midpoint. For complex mixed-frequency intervals the same calculation uses $M=q$ and $|s|<2$. These bounds cover the full intervals. [Real-axis cover](real_axis.py) · [Published finite-height verification by Platt](https://arxiv.org/pdf/1305.3087#page=14)

## Replay

The release records three arithmetic runs, the source digest, kernel replays and the axiom audits. The host and character witnesses, the trace identity, the width floor, and the bracket-to-uniform-slope implication are Lean theorems. The zero enclosures use ball arithmetic. The Mellin continuation and the positive-weight extremum argument are written mathematics.

For an existing cache put the pinned Lean and Python environment on `PATH`, install `requirements.txt` and run

```sh
bash verify-release.sh /path/to/mathlib4 /path/to/upstream/src /path/to/build-cache /path/to/theta-prime-trace
```

A fresh Linux replay is `bash verify-release.sh`. It downloads the pinned source release, compiler, Mathlib and two trace source files. The manifest covers the scripts, proof, theorem statements and expected output. [Replay script](verify-release.sh) · [Kernel audit](check_lean.py) · [Expected output](evidence/check.txt)
