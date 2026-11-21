# Proof

## Field and local observable

Write $K=\mathbf Q(\zeta_{16})^+$, $a^2=-571$, $d^2=-599$ and $E=K(a,d)$. The real quartic field has conductor 16. The two quadratic fields have coprime odd discriminants, also coprime to 16. The three fields are linearly disjoint. Equivalently their Dirichlet character groups are independent by the Chinese remainder theorem. Thus $G=\mathrm{Gal}(E/\mathbf Q)=C_2\times C_2\times C_4$. Complex conjugation is $c=(-1,-1,1)$. Its fixed field is $F=K(ad)$, of degree 8. These identifications follow directly from the generators of the Lean host. We can discard primes above $2,571,599$ throughout. Any fixed finite set in the prime sum disappears for all sufficiently large $y$. [Field generators](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Formalism/FrobeniusOrbits/Statement.lean#L53) · [Example and ranks](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Formalism/LocalRingCardinality/SpanPowerMembership/ExampleRanksRankComplement.lean#L6)

The record is `exampleRecord` at `H0` with its local-gap witness. Its two infinity exponents are 100 and 118, constant across the CM type. The constant-weight congruence `UpperTheory.Mu.muq_of_const` applies at each index. The source proves that the two required primes, 751 and 991, split completely. The character-frame existence theorem therefore gives a frame on this record. The level and character choices give an admissible `GlobalTheory.ChainOnto.Com`. This step supplies one fixed datum for all primes and widths. [Frame proof](HostFrame.lean) · [Admissibility proof](HostAdmissible.lean) · [Upstream character-frame construction](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AbelianVarietyTheory/Proofs/DivisorProofs/IsTwistPlacesAboveK.lean#L79) · [Upstream admissible-data constructor](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AbelianVarietyTheory/Definitions/LatticeFrames/ExistsStage.lean#L70)

For a prime $w$ of $E$, let $T_0(w),T_1(w)$ be the two first traces of the input family and of its next parameter. At step 1 the source and target characters are $\chi_S=\mu_1\eta^{-1}$ and $\chi_T=\mu_1$. Write $\lambda_w(\chi)$ for evaluation at a local uniformizer and $U(w)>0$ for the first trace of the untwisted segment. The finite trace lemma at the singleton $\lbrace w\rbrace$ gives

$$
J(w)=\frac{T_1(w)-\lambda_w(\chi_S\chi_T^{-1})T_0(w)}{U(w)}
      =\lambda_w(\mu_1)\lambda_w(\eta)^{-1}.
$$

The host extension $E/F$ is unramified at every finite place. The concrete Lean theorem proves this from the CM datum and the ramification index. The admissibility record gives restriction exponents 0 for $\mu_1$ and 1 for $\eta$. At a good $v$ of $F$, embed a base uniformizer into the local factors above $v$. Every valuation is one since $v$ is unramified. Unramified characters ignore the remaining units. The upstream lemma `GlobalTheory.LFun.ResEta.prod_val` identifies the product for a restriction exponent $r$ with $\epsilon(v)^r$, where $\epsilon(v)=1$ at split primes and $-1$ at inert primes. Consequently

$$
\prod_{w\mid v}J(w)=1\cdot\epsilon(v)^{-1}=\epsilon(v).
$$

This is the product that the prime sum uses. The analytic argument can now use its quadratic Euler factors. The kernel-checked statements concern the parameter recursion and its admissible data. The input family is arbitrary because the trace of that family cancels in the one-step difference. [Pinned trace lemma](https://github.com/leecrossley/theta-prime-trace/blob/7a7909d2dc7bb9d8d0ec7b2dc8c78e70f93782d7/ThetaPrimeTrace.lean#L68) · [Restriction product](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AbelianVarietyTheory/Proofs/FoundationProofs/RamificationIdxEqOneOfIsSquare.lean#L545) · [Splitting sign in Lean](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AbelianVarietyTheory/Definitions/Foundations/NormEulerPolySubOneLe.lean#L396) · [Fibre proof](FiberTrace.lean)

## The Mellin transform

Let $A_b(y)=-\sum_v\epsilon(v)\log Nv\thinspace W_b(Nv/y)$ with the weight from the README, omitting the bad primes. The sum is zero near $y=0$. The transform of the weight, continued at 0 and 1, is

$$
H_b(s)=\int_0^\infty W_b(t)t^{s-1}\thinspace dt
=\frac{(b^s-1)(b^s-b)}{s(s-1)} .
$$

This follows by integrating the two branches on $[1,b]$ and $[b,b^2]$. Let $\xi$ be the primitive even quartic character modulo 16 with $\xi(5)=i$. Let $\chi_p$ denote the primitive quadratic character of discriminant $-p$ for $p=571,599$. Artin factorization, or the Euler factors of the abelian group $G$, gives

$$
L(s)=\frac{\zeta_E(s)}{\zeta_F(s)}
 =\prod_{p\in\lbrace571,599\rbrace}\prod_{j=0}^3 L(s,\chi_p\xi^j),
$$

with each character taken primitively. They are exactly the characters taking $c$ to $-1$. At a good rational prime with Frobenius $g$, let $f$ be the order of its image in $G/\langle c\rangle$. There are $8/f$ primes of $F$ of norm $p^f$. Their signs agree and are determined by $g^f\in\lbrace1,c\rbrace$. Summing the eight odd characters on $g^n$ gives zero unless $f\mid n$, and otherwise gives $+8$ or $-8$ according as $g^n=1$ or $c$. This verifies the logarithmic Euler coefficients with the correct base norms and multiplicities. Finite Euler corrections below are holomorphic for $\Re s>0$. [Group check](independent.py)

Write $M_b(s)=\int_1^\infty A_b(y)y^{-s}\thinspace dy/y$. Its difference from the integral starting at zero is entire. Initially for $\Re s>1$, termwise integration and the Euler products give

$$
M_b(s)=H_b(s)\left(\frac{L'}{L}(s)-\frac{\zeta_F'}{\zeta_F}(2s)+R(s)\right)+U(s),
$$

where $R$ is holomorphic in $\Re s>1/3$ and $U$ is entire. Indeed $L'/L$ contributes the negative of every prime-power term. The omitted squares must be restored with positive sign since $\epsilon(v)^2=1$. They are the prime terms of $-\zeta_F'(2s)/\zeta_F(2s)$, whose higher powers converge absolutely for $\Re s>1/4$. The signed powers of exponent at least three converge absolutely for $\Re s>1/3$. These statements also cover the deletion of finitely many local factors.

The transform has no real singularity for $s>1/2$. This follows from nonvanishing on the real segment $[1/2,1]$ for all eight factors, their absolutely convergent Euler products to the right of 1, and the corresponding Euler product of $\zeta_F(2s)$. The computer certificate for the real segment is given below. Platt’s finite-height verification independently covers every conductor here. At $s=1/2$ the simple pole of $\zeta_F$ gives residue $H_b(1/2)/2=B_b$. For $\omega\ne0$, $\zeta_F(1+2i\omega)\ne0$ by the prime ideal theorem. Thus $M_b$ at $1/2+i\omega$ has residue $m_\omega H_b(1/2+i\omega)$, where $m_\omega$ counts zeros of all eight primitive factors there, with multiplicity. Zero residue here means holomorphy. The factors form a conjugation-stable set. [Platt, Theorems 7.1 and 7.2](https://arxiv.org/pdf/1305.3087#page=14)

## Four ordinates and the inequalities

Put $x=\log y$, $f_b(x)=e^{-x/2}A_b(e^x)$, and $z=x+\log b$. Direct simplification gives

$$
H_b(1/2+i\gamma)=e^{i\gamma\log b}h_b(\gamma),
\qquad
h_b(\gamma)=\frac{\sqrt b\thinspace(b+1-2\sqrt b\cos(\gamma\log b))}{1/4+\gamma^2}
\ \ge\ \frac{B_b}{2(1/4+\gamma^2)} .
$$

Take one positive zero ordinate from each row of the table. Endpoints are integers in units of $10^{-8}$. Simplicity or uniqueness in a bracket is unnecessary.

| Primitive character (Conrey label) | Left endpoint | Right endpoint |
|---|---:|---:|
| 9584.5989 | 12294899 | 12294900 |
| 4568.1141 | 27941486 | 27941487 |
| 4792.1197 | 31952267 | 31952268 |
| 9136.1141 | 37800727 | 37800728 |

For each sign, form the nonnegative trigonometric polynomial

$$
P_\pm(z)=\prod_{j=1}^{4}(1\pm\cos(\gamma_j z)).
$$

The constant Fourier coefficient is 1. For each signed sum $\sum_j e_j\gamma_j$, $e_j\in\lbrace-1,0,1\rbrace$, involving at least two nonzero coefficients, the certificate excludes both zero and *every* zero ordinate of $L$ at that frequency. Up to overall sign these are $(3^4-1-8)/2=36$ interval tests for each of eight factors. The four singleton ordinates and their negatives are distinct. Thus the only frequencies of $P_\pm$ meeting poles are zero and the eight singletons. The singleton coefficients are $\pm1/2$.

To bound the liminf, suppose $f_b$ has an eventual lower bound. Otherwise the inequality is immediate. Subtract this lower bound and discard a bounded initial interval. The Laplace transform of the resulting nonnegative function converges absolutely for $\Re z>0$. Indeed it converges initially for $\Re z>1/2$ by the prime sum, and the continuation $M_b(1/2+z)$, adjusted by an entire function and a multiple of $1/z$, is holomorphic at every positive real $z$. Landau’s theorem would give a real singularity at any positive abscissa of convergence. To see this, Taylor-expand just to the right of a purported regular real abscissa. At a point just to its left, all terms under the integral are nonnegative. Monotone convergence sums the exponential series there, extending convergence past the abscissa. The same argument with the sign reversed applies when $f_b$ has an eventual upper bound.

In the respective one-sided cases, evaluating each term of the finite polynomial gives the Abel mean

$$
\lim_{\delta\downarrow0}\delta\int_0^\infty
   f_b(x)P_\pm(x+\log b)e^{-\delta x}\thinspace dx
   = B_b\ \pm\sum_{j=1}^4 m_{\gamma_j}h_b(\gamma_j).
$$

Indeed a term $e^{-i\omega(x+\log b)}$ evaluates $\delta e^{-i\omega\log b}M_b(1/2+\delta+i\omega)$, which tends to the phase-corrected residue. The Abel mean of $P_\pm$ alone is 1. Eventual bounds on $f_b$ can be integrated against the positive weight, since the contribution of any fixed initial interval tends to zero. It follows that the liminf is at most the minus expression and the limsup is at least the plus expression. This argument uses only finitely many boundary residues.

If $u_j$ is the upper endpoint of bracket $j$ as a real number, exact rational arithmetic gives

$$
\sum_{j=1}^4\frac{1}{2(1/4+u_j^2)}>\frac{61027}{10000}.
$$

Each multiplicity is at least one. Dividing by the positive $B_b$ proves the two strict inequalities in the README. Their scope is pointwise in $b$ with universal constants. The locations at which the inequalities are witnessed can depend on the width.

The same test can be used on another quadratic host. For this kernel the data to look for are a finite set of positive zero ordinates with $\sum_j[2(1/4+\gamma_j^2)]^{-1}>1$, a real zero-free segment and exclusion of the mixed signed frequencies in the corresponding product. With those checks, the minus bound is negative for every fixed width and the plus bound is positive. The host and the relevant primitive factors must be identified afresh.

## Reproducing the finite certificates

`verify.py` uses the real Hardy function of each primitive character. Opposite certified signs imply a zero by continuity. Evaluation over an entire mixed-frequency interval excludes a zero anywhere in it. `real_axis.py` independently checks all residue classes against the Conrey labels and covers $[1/2,1]$ by rational intervals. With $C(x)=\sum_{n\le x}\chi(n)$ it computes an integer $M$ with $|C(x)|\le\min(x,M)$ from one full period. Abel summation gives $L(s)=s\int_1^\infty C(x)x^{-s-1}dx$. For $1/2\le s\le1$, the absolute derivative is bounded by $4\sqrt M+4\sqrt M\log M+4$ by splitting the two elementary integrals at $M$. The value at each midpoint exceeds that derivative bound times the interval radius. [Hardy checks](verify.py) · [Real-axis cover](real_axis.py) · [FLINT Hardy function](https://flintlib.org/doc/acb_dirichlet.html)

`independent.py` constructs the characters directly from Legendre symbols and the even character taking 5 to $i$. It evaluates the Dirichlet series by residue classes, using a Hurwitz-zeta expansion, without calling a Dirichlet character, L-function or Hardy evaluator. For $\Re s=1/2$, $|s|<2$, split off the first $q$ terms. The rest is

$$
q^{-s}\sum_{k\ge1}\frac{(-1)^k(s)_k}{k!}\thinspace
\zeta(s+k,3/2)\sum_{a=1}^{q}\chi(a)(a/q-1/2)^k .
$$

The constant term cancels. Indeed in each subsequent period $m\ge1$, expand $(m+1/2+(a/q-1/2))^{-s}$ about $m+1/2$. The series converges absolutely after the cancellation. The $k$-th coefficient has absolute value at most $k+1$, the Hurwitz factor is at most $4(3/2)^{-k}$ by integral comparison, and the inner sum has magnitude at most $q\thinspace2^{-k}$. The tail beyond $k=K$ thus has magnitude at most $\sqrt q(2K+5)/3^K$. The implementation takes $K=54$ and includes this bound in both coordinates.

For the four endpoint tests use $\Lambda(s)=(q/\pi)^{(s+1)/2}\Gamma((s+1)/2)L(s,\chi)$. By the functional equation $\Lambda(1/2+it)\overline{\Lambda(1/2)}$ is real. Its opposite signs certify the crossing independently of any choice of root-number phase. At each mixed interval midpoint the second verifier bounds the variation by Abel summation using $|C(x)|\le\min(x,q)$ and $|s|<2$. A uniform derivative bound there is $4\sqrt q+2(4\sqrt q\log q+4)$. All individual factors are separated from zero over their full intervals. [Second implementation](independent.py)

Each Lean module added here is elaborated afresh and then checked with `leanchecker`. The published axiom receipt covers all five named declarations. A pinned existing source and Mathlib cache can be used as follows, with the pinned Lean and the Python environment on `PATH`. The fourth location holds the two sources from the cited trace revision. [Release script](verify-release.sh) · [Replay output](evidence/check.txt)

```sh
bash verify-release.sh /path/to/mathlib4 /path/to/lean/src /path/to/build-cache /path/to/theta-prime-trace
```
