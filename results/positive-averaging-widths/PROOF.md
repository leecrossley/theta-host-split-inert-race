# Positive averages of the two theta prime races

## The prime races and the question

For $`(p,q)=(547,727)`$ or $`(571,599)`$, put $`K=\mathbf Q(\zeta_{16})^+`$, $`F=K(\sqrt{pq})`$ and $`E=K(\sqrt{-p},\sqrt{-q})`$. At the good prime ideals of $`F`$ take the sign $`\epsilon_v=1`$ for splitting in $`E/F`$ and $`-1`$ for inertia. The local normalized theta trace jumps give this sign. Use the same finite omissions as the respective source races. Put $`a=\log2`$, $`r=\sqrt2`$, and at $`b=r_{\rm move}=2`$ define

```math
w(u)=\frac{\max(0,\min(2-2/u,4/u-1))}{2r(r-1)^2},\qquad
S(y)=-\frac1{\sqrt y}\sum_v\epsilon_v\log(Nv)\,w(Nv/y),\qquad
D(y)=\frac{S(2y)-S(y)}{a}.
```

The odd factors of $`\zeta_E/\zeta_F`$ are the two quadratic characters of conductors $`p,q`$, each twisted by all four even characters of $`K`$. The strip explicit formula for the level has a term at each nontrivial zero $`\rho`$ proportional to $`e^{(\rho-1/2)x}`$ at log height $`x=\log y`$. On the line, the absolute value of a single level coefficient at $`t=\Im\rho`$ is

```math
G(t)=\frac{3-2r\cos(at)}{2(r-1)^2(t^2+1/4)}.
```

A change of centering contributes a phase of modulus one. The move coefficient has the additional factor $`(e^{iat}-1)/a`$. Both coefficient magnitudes are positive when $`0\lt t\le1`$. These facts hold independently of any assumptions on zeros elsewhere. [547/727 trace](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-short-shifts/RaceTrace.lean#L83) · [571/599 trace](https://github.com/leecrossley/theta-host-split-inert-race/blob/8fca7cad37185e3d7951f0dffff5e2442d1dfa24/HostAdmissible.lean#L80) · [Strip formula for prime-only levels and moves](https://github.com/leecrossley/theta-host-split-inert-race/blob/cross-field-filter-obstruction-v1.0.0/results/cross-field-filter-obstruction/PROOF.md)

Let $`\mu`$ be a compactly supported nonnegative probability measure on the log-height line. It averages $`T(e^x)`$ to $`\int T(e^{x+s})\,d\mu(s)`$ for $`T=S`$ or $`D`$. Its transfer at a critical-line zero of ordinate $`t`$ is $`\widehat\mu(t)=\int e^{its}\,d\mu(s)`$. Its width is $`\sup\mathrm{supp}\mu-\inf\mathrm{supp}\mu`$. Translating the measure preserves its width and the zeros of its transfer. We ask for the minimum width that makes the transfer zero at every odd-factor zero with $`|\Im\rho|\le1`$ on either field or on both fields together.

**Theorem.** Write $`m_A,m_B,m_{AB}`$ for these three minimum widths. Then

```math
88.436\lt m_A<88.440,\qquad47.698\lt m_B<47.702,\qquad
103.945\lt m_{AB}<103.949 .
```

For each lower endpoint and every admissible average of at most that width, some indicated positive ordinate has $`|\widehat\mu(t)|\gt1/1000001`$. Each upper endpoint admits a positive symmetric atomic filter and a smooth nonnegative density. The same filter cancels the selected contributions in both the level and the dyadic move. On the joint problem one filter serves both fields. The contribution at zero frequency from the prime-square bias is retained by normalized averaging.

## Counting every zero in the band

For each odd primitive character $`\chi`$ let $`\Lambda_\chi(s)=(q_\chi/\pi)^{(s+1)/2}\Gamma((s+1)/2)L(s,\chi)`$. This entire function of order one has the nontrivial zeros and is nonzero at the strip endpoints. At real $`\sigma\gt1`$,

```math
\Re\frac{\Lambda_\chi'}{\Lambda_\chi}(\sigma)
=\sum_{\rho=\beta+i\gamma}
\frac{\sigma-\beta}{(\sigma-\beta)^2+\gamma^2}.
```

All terms are positive and include multiplicity. To see the identity, subtract the logarithmic derivatives of the genus-one product at $`\sigma`$ and $`1-\sigma`$. The constants cancel. The functional equation makes the real part on the left twice the required real part. On the right the involution $`\rho\mapsto1-\overline\rho`$ in this same character's zero set makes the two real sums equal. They are absolutely convergent by the usual zero count. For negative ordinates of $`\chi`$ the retained sign brackets come from the positive ordinates of $`\overline\chi`$. These are two different uses of symmetry. [Kadiri, explicit logarithmic derivatives](https://www.cs.uleth.ca/~kadiri/articles/zerosL-04-09-07.pdf)

The table supplies 810 positive, disjoint-within-factor Hardy sign brackets. The verifier tests both endpoints, the primitive conductors and the exact quadratic/quartic values on every odd unit class modulo $`16p`$. At $`\sigma=3/2`$ subtract the positive terms contributed by one enclosed critical-line zero at each retained ordinate and the corresponding negative partners. In all sixteen characters the residual mass is below $`1/5`$. This calculation evaluates the completed logarithmic derivative directly using ball series. [Integer brackets](evidence/band-brackets.tsv) · [Hardy signs and mass](band_zeros.py)

A second calculation uses $`\sigma=2`$ and the real Euler sum for $`L'/L`$. It sieves every prime power through 4,096 with the exact real part of its quartic character value. Its absolute remainder is at most $`(\log4096+1)/4096`$ by comparison with the integral of $`(\log u)/u^2`$. The gamma derivative uses $`\psi(3/2)=2-\gamma_{\rm E}-2\log2`$. At this displacement every residual is below $`1/4`$. [Independent Euler-sum count](euler_mass.py)

For either displacement, an unaccounted zero with $`|\gamma|\le1`$ has $`1/2\le u=\sigma-\beta\le2`$. The inequality $`(u-1/2)(2-u)\ge0`$ gives $`u/(u^2+\gamma^2)\ge2/5`$. Either mass check excludes it, including any additional multiplicity at a retained zero. The listed positive brackets in the band are pairwise disjoint across all sixteen factors. There are nine in 547/727 and seven in 571/599. Conjugation gives their negative partners, and the mass check excludes a real zero. The argument covers the full strip and uses no global critical-line hypothesis.

## The continuous extremal problem

For a finite set $`T=\{t_1,\ldots,t_n\}`$ of distinct positive ordinates, translate any measure of width at most $`2R`$ into $`[-R,R]`$. Symmetrizing it preserves vanishing at all selected ordinates. The image of a symmetric measure under absolute value is a probability on $`[0,R]`$. Its constraint is the vector equation

```math
\int v(x)\,d\eta(x)=0,\qquad
v(x)=(\cos(t_1x),\ldots,\cos(t_nx)).
```

Conversely such an $`\eta`$ splits each atom or element of mass equally between the two reflected shifts. Thus feasibility means that zero belongs to the convex hull of $`v([0,R])`$. That hull is compact in finite dimension. A finite feasible radius exists, for instance by convolving the symmetric two-point measures at $`\pm\pi/(2t_j)`$. The minimum is attained by weak compactness of the probability measures after centering.

For a lower bound choose real $`d_j`$ and set $`h(x)=\sum_j d_j\cos(t_j x)`$. At the indicated lower half-width the certificate gives $`h(x)\gt10^{-6}`$ on the entire interval and $`\sum_j|d_j|\lt1.000001`$. For any measure, symmetric or otherwise, center its support first. Then $`\sum_jd_j\Re\widehat\mu(t_j)\gt10^{-6}`$. The triangle inequality gives the claimed strict transfer bound. Centering does not change the magnitudes.

For the upper bound each rational shift $`x_k\ge0`$ labels a pair $`\pm x_k`$. With $`n+1`$ such pairs the system for their total masses is

```math
A_{0k}=1,\qquad A_{jk}=\cos(t_jx_k)\ (1\le j\le n),\qquad
A w=(1,0,\ldots,0)^{\mathsf T}.
```

The true ordinates lie in the certified intervals. The matrix is invertible for all choices in these intervals and each exact weight exceeds $`1/1000`$. Hence the exact solution, rather than rounded decimal masses, defines a probability with exactly zero transfers. By continuity, moving all positive atom locations towards zero by a sufficiently small common factor preserves invertibility and strict positivity after recomputing the weights. This gives the strict upper inequalities. The dual coefficients follow increasing positive ordinate. The data gives integer dual coefficients in units of $`10^{-12}`$ and atom positions in units of $`10^{-9}`$. [Filter data](evidence/design.json)

## Two interval checks of the filters

The first proof evaluates every cosine sum at 16,385 grid vertices. With $`\Delta`$ the mesh length, the possible loss between vertices is at most $`\Delta^2\sum_j|d_j|t_j^2/8`$, by linear interpolation. It solves each upper system with a preconditioned interval linear solve. The second proof adaptively bisects the continuous domain. On a cell with midpoint $`m`$ and radius $`h`$ its lower enclosure is

```math
H(m)-|H'(m)|h-\tfrac12 h^2\sum_j|d_j|t_j^2,\qquad
H(x)=\sum_jd_j\cos(t_jx).
```

It terminates on 233, 163 and 948 closed cells respectively. For an independent upper check take the exact rational inverse $`P`$ of the matrix of dyadic ball midpoints. Direct interval products bound $`\|PA-I\|_\infty\lt10^{-6}`$ and $`\|P e_0\|_\infty\lt1`$. The Neumann series bounds the maximum error between the true weights and $`Pe_0`$ by $`1/999999`$. Each weight still exceeds $`1/1000`$. All matrix and trigonometric inequalities include uncertainty in every zero ordinate. [Mesh and interval solve](positive_widths.py) · [Cell cover and rational inverse](independent.py) · [Verification output](evidence/check.txt)

## Which averages attain the minimum

At a minimal radius the origin is on the boundary of the convex hull. Indeed if a ball about zero were inside, uniform continuity of $`v`$ would keep zero in the hull at all sufficiently close smaller radii. More explicitly, all directional support functions of the hull change uniformly by at most the modulus of continuity of $`v`$ when the outermost interval is shortened. An interior ball gives a positive lower bound on each such support function. The hull has full affine dimension at every positive radius, since $`1,\cos(t_1x),\ldots,\cos(t_nx)`$ are linearly independent there. This follows from their even derivatives at zero and the Vandermonde matrix for the distinct $`t_j^2`$.

A supporting hyperplane at zero gives a nonzero linear combination of the cosines, nonnegative on the whole interval. This real analytic function is not identically zero and has finitely many zeros on the closed interval. Every centered extremizing measure integrates it to zero, even without symmetrizing, so its support is contained in that finite set and its reflections. In particular every minimizer is atomic. Convolving an extremizer with a nonnegative smooth probability bump preserves its Fourier zeros and increases its span by arbitrarily little. This proves that smooth densities have the same infimum and cannot attain it. Taking the strict room in the upper constructions also gives smooth densities at each stated upper width.

Convolving two positive compact measures adds their support widths. For example, points within an arbitrary distance of the two upper endpoints have positive mass in their neighborhoods, so their sums approach the upper endpoint of the convolution support. The lower endpoints behave in the same way. It follows that the best individual filters composed in succession use more than $`32.185`$ log units beyond the best common filter. [Checked width deduction](AveragingWidths.lean)

## What remains visible in the prime scores

For log width at most 20, a translated probability is supported in $`[-10,10]`$. The first positive ordinate is in $`(0.04853,0.04854)`$ on 547/727 and $`(0.12294,0.12296)`$ on 571/599. At either first ordinate $`t`$ the real part of the centered transfer is at least $`\cos(10t)\gt0`$. No comparison-field zero occurs there. For a target level or move averaged with this probability, subtract any fixed signed compact filter of the opposite level and any real intercept. Define its upper logarithmic RMS as

```math
\limsup_{X\to\infty}\left(\frac1X\int_0^X|\mathrm{error}(x)|^2\,dx\right)^{1/2}.
```

A surviving off-line term in the right half-plane gives infinite upper RMS by the pole argument. If none survives, Parseval for the convergent critical-line series bounds the squared RMS below by twice the squared modulus of this one private coefficient. The same dichotomy holds at heights $`n^k`$ with weights $`1/n`$ normalized by $`\log N`$, for each fixed $`k\gt0`$. The previously proved argument allows a compact regular target weight. Translating the tent and integrating against a compact probability preserves the finite-measure bound on its second derivative in log coordinates. The comparison field's coefficient at the selected frequency is zero for every filter. [Strip expansion, pole argument and power sampling](https://github.com/leecrossley/theta-host-split-inert-race/blob/cross-field-filter-obstruction-v1.0.0/results/cross-field-filter-obstruction/PROOF.md)

The lower bound for a level is therefore $`\sqrt2\,G(t)\cos(10t)`$, and the bound for its normalized move has the additional multiplier $`2\sin(at/2)/a`$. The interval check gives respectively $`2.50`$ and $`0.1214`$ on 547/727, and $`0.946`$ and $`0.1162`$ on 571/599, strictly. The bounds concern the actual prime-only scores and permit arbitrary zeros outside the certified low band. [Numerical checks](positive_widths.py)

The five local Lean deductions check the kernel floor, averaging inequality, transfer separation, inverse margin and joint width saving. The averaging and separation lemmas are stated for finite sums. Their extension to compactly supported probability measures, the Neumann-series bound and the interpolation and Taylor enclosures are not in Lean. The analytic arguments above and the ball-arithmetic computations are separate parts of the proof. [Lean source](AveragingWidths.lean) · [Replay instructions](README.md#replay)
