# Prescribing consecutive doublings on two theta fields

## The scores and the arithmetic input

Use $`A=(547,727)`$ and $`B=(571,599)`$. For each pair $`i=(u,v)`$ let $`K=\mathbf Q(\zeta_{16})^+`$, $`F_i=K(\sqrt{uv})`$ and $`E_i=K(\sqrt{-u},\sqrt{-v})`$. A good prime ideal $`w`$ of $`F_i`$ has sign $`\epsilon_w=1`$ for splitting and $`-1`$ for inertia in $`E_i/F_i`$. This is the product over its fibre of the normalized first-step theta trace jumps. The two admissible host data and the common-height convention have been established separately. Omit the places above 2 and the two conductors in the finite prime examples. Any additional fixed finite omissions, including the exceptional primes of the trace data, have left the compact support on a sufficiently late tail. [547/727 trace](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-short-shifts/RaceTrace.lean#L83) · [571/599 trace](https://github.com/leecrossley/theta-host-split-inert-race/blob/8fca7cad37185e3d7951f0dffff5e2442d1dfa24/HostAdmissible.lean) · [Common-height finite-zero theorem](https://github.com/leecrossley/theta-host-split-inert-race/blob/joint-levels-and-moves-v1.0.0/results/joint-levels-and-moves/PROOF.md)

Set $`a=\log 2`$, $`r=\sqrt2`$, $`c=2(r-1)^2`$, $`B_0=rc`$ and $`W(t)=\max(0,\min(2-2/t,4/t-1))`$. Use

```math
S_i(y)=-\frac{1}{B_0\sqrt y}\sum_w\epsilon_w\log(Nw)W(Nw/y),
\qquad D_i(y)=\frac{S_i(2y)-S_i(y)}a .
```

The odd primitive constituents are the quadratic characters of conductors $`u,v`$, each twisted by the four powers of the primitive even quartic character modulo 16. For each odd prime conductor the primitive conductors are $`p,8p,16p,16p`$. The bound $`|V_i|\le M_i`$ below, the six-step spectral certificate and the infinite-height conclusions invoke Alpöge's announced theorem that each primitive Dirichlet L-function has its zeros in the critical strip on the critical line. This is GRH for these functions, and the announcement has not yet been refereed. Its statement, the printed axiom audit and the independent kernel records are pinned by digest. The finite prime-ideal catalogue uses no zero-location theorem. The density argument permits every rational relation among the ordinates. [Announcement](https://alpo.ge/remark.pdf) · [Upstream theorem](https://github.com/leanprover-community/riemann-hypothesis/blob/371a1a26d5e17f604428e3a8ce6ce341cf517434/src/Main.lean) · [Independent verification](https://math.mit.edu/~drew/grh/) · [Pinned record check](check_source.py)

## First occurrence on two consecutive doublings

For a sign word $`s\in\{-1,1\}^4`$ use the coordinate order $`V(y)=(D_A(y),D_B(y),D_A(2y),D_B(2y))`$. Its strict margin exceeds $`q`$ if $`s_jV_j(y)\gt q`$ in all four coordinates. Write $`M(y)=\min_j V_j(y)`$. Among integer starts from 1 through 1,469, fifteen words already occur with margin greater than 2. In binary order, with the high bit corresponding to the first coordinate and a set bit denoting a gain, one witness for each word from 0000 through 1110 is given below. Both sieves check every one of the sixty signed inequalities. The remaining word 1111 consists of four gains. [Two counts and transcripts](evidence/first-catalogue.txt)

| Code | 0000 | 0001 | 0010 | 0011 | 0100 | 0101 | 0110 | 0111 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Start | 1,469 | 9 | 89 | 33 | 17 | 15 | 20 | 28 |

| Code | 1000 | 1001 | 1010 | 1011 | 1100 | 1101 | 1110 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Start | 45 | 8 | 39 | 150 | 56 | 248 | 629 |

**Theorem.** At $`q=3947/2000`$, the infimum of real starts of word 1111 lies strictly between $`20{,}306{,}533.59533`$ and $`20{,}306{,}533.59534`$. The first integer start is 20,306,534. The infimum of heights by which every word can finish lies between 81,226,134.38132 and 81,226,134.38136. It is unattained because the controlling gain equals the threshold at the crossing. Integer starts first complete the catalogue at 81,226,136. Thus for the two-field four-gain problem a deadline from height one must exceed the real completion infimum. The single-field question asking for one gain and one loss has different quantifiers. [Direct integer witness and real exclusion](word_cover.py)

### Exact prime sums and independent counts

For a good odd rational prime $`p`$ let $`\chi_u(p)=(p/u)`$ and similarly for $`v`$. Its residue degree in $`F_i`$ is the least $`f\in\{1,2,4\}`$ with $`p^f\equiv\pm1\pmod {16}`$ and $`(\chi_u(p)\chi_v(p))^f=1`$. The quadratic sign is $`\chi_u(p)^f`$ by reciprocity. The $`8/f`$ base primes above $`p`$ all have norm $`n=p^f`$. Their logarithmic weight, including the inert-minus-split sign, is $`8e_p\log p`$, where $`e_p=-\chi_u(p)^f`$. This is the local trace calculation in the scores defined above. The two extra omitted rational primes for $`A`$ are 343057911569 and 343070636977, and for $`B`$ they are 413718278401 and 603339156001. All exceed 324,944,000, the largest norm in the real-start count. The two conductors and 2 in each field are omitted throughout. [Local host and score definitions](#the-scores-and-the-arithmetic-input) · [Five fixed omissions in each field](https://github.com/leecrossley/theta-host-split-inert-race/blob/joint-levels-and-moves-v1.5.0/results/joint-levels-and-moves/deadline_arithmetic.py)

Write $`P(l,h)=\sum_{l\lt n\le h}e_p\log p`$ and $`Q(l,h)=\sum_{l\lt n\le h}e_p\log(p)/n`$. Substitution of the two tent branches gives

```math
\sqrt y\,S_i(y)=\frac8{B_0}
\bigl(2P(y,2y)-P(2y,4y)+y[4Q(2y,4y)-2Q(y,2y)]\bigr).
```

The primary sieve marks the two progressions $`16m\pm1`$ and enumerates higher residue degrees separately. It builds quadratic-residue tables from squares, uses a logarithm mesh of size $`2^{15}`$ and encloses the remainder by alternating Taylor polynomials of degrees three and four. Its scales for logarithms and reciprocal weights are $`2^{58}`$ and $`2^{96}`$. The second sieve runs through every odd rational prime, computes both Euler criteria and the residue degree independently, and sums into disjoint bins. Its mesh has size $`2^{12}`$ and its scales are $`2^{62}`$ and $`2^{102}`$. For the ratio of a number and the preceding mesh node it uses the positive series for $`2\,\mathrm{arctanh}\,t`$, with $`0\le t\lt 1/8192`$. After the fifth power the tail is below $`2t^7/(7(1-t^2))\lt 2^{-62}`$. The table entries in both runs are directed Arb floors. Cumulative lower endpoints are subtracted from cumulative lower endpoints only when recovering sums of the already enclosed disjoint bins, and likewise for upper endpoints. Individual products and the cumulative signed sums fit in their integer types. For instance $`\sum_{n\le N}\log(N)/n\le\log(N)(1+\log N)`$ bounds the unscaled reciprocal totals. [Primary arithmetic](word_sieve.cpp) · [Independent arithmetic](word_sieve_recheck.cpp) · [Finite-score identity](patterns.py)

### Covering all real heights

The primary cover works in $`x=\log y`$. On every open knot interval a level $`P(x)=S_i(ke^x)`$ satisfies $`P''=P/4`$. On a closed cell of log width at most $`h`$, its departure from its secant is bounded by $`\sup |P|h^2/32+hJ/4`$, where $`J`$ bounds the sum of the absolute jumps of $`P'`$. Indeed a smooth second derivative has Dirichlet Green-kernel integral at most $`h^2/8`$, and a single derivative jump contributes at most its magnitude times $`h/4`$. For log centres $`z=2ky`$ between $`L`$ and $`U`$, knots occur at norms $`n=sz`$, with $`s\in\{1/2,1,2\}`$. There are at most $`1+sUh`$ such integers for each $`s`$ when $`h=(y_{\rm right}-y_{\rm left})/y_{\rm left}`$. The jump coefficients of the log tent $`g`$ are bounded by $`2/B_0,4\cosh(a/2)/B_0,2/B_0`$ in that order. Each carries at most $`8\log(2U)/\sqrt{sL}`$. Bounding the number of active norms by $`2U+1`$ and the tent by its central value gives the smooth term. Adding the two level errors and dividing by $`a=\log 2`$ bounds the error of each move.

The independent cover uses the height itself. For $`H_k(y)=\sqrt y\,D_i(ky)`$ every open knot segment is affine. In $`\sqrt y\,S_i(ky)`$ a term has slope jumps at $`n/(4k),n/(2k),n/k`$, with coefficients $`4,-6,2`$ respectively, multiplied by $`8\sqrt{k}\,e_p\log(p)/(B_0 n)`$. For $`y\in[l,r]`$ and each $`s=4,2,1`$ there are at most $`1+sk(r-l)`$ integers in $`[skl,skr]`$, each with $`\log(p)/n\le\log(skr)/(skl)`$. Add the two levels’ jump bounds and divide by $`a`$ to bound the total jump size $`J_k`$. Since $`-q\sqrt y`$ is convex for $`q\ge0`$, the maximum at the two endpoints of $`H_k(y)-q\sqrt y`$, plus $`(r-l)J_k/4`$, bounds it above on the whole cell. At least one of the four coordinates is below the selected threshold on every exclusion cell. The meshes have different dyadic resolutions, 4,096 and 3,500. [Both closed-interval implementations](word_cover.py)

The slope check also gives inequalities on whole intervals. Let $`u`$ be the right slope of $`H_k`$ at $`l`$. All one-sided slopes on the cell lie in $`[u-J_k,u+J_k]`$, so its range lies in $`H_k(l)+[0,r-l][u-J_k,u+J_k]`$. Away from knots the sign of the move derivative equals that of $`2yH'_k(y)-H_k(y)`$. Strict same-sided bounds on every open segment imply strict monotonicity across the knots by continuity. For $`0\lt y\le6`$ the first $`B`$ move vanishes because its smallest eligible norm is 49.

| Certified interval | Inequality throughout the interval |
| --- | --- |
| $`[6,45700]`$ | At least one move is strictly negative. |
| $`[45700,46701]`$ | The first $`B`$ move increases strictly, the second $`B`$ move decreases strictly, and both $`A`$ moves exceed 2. |
| $`[46701,20288400]`$ | At least one move is strictly below the early integer maximum defined below. |
| $`[20288400,20309000]`$ | The first $`A`$ move increases strictly, and all three other moves exceed 2. |

The primary proof checks 52,935 exclusion cells before the early window and 38,997 between the windows. The second checks 45,267 and 33,890 respectively. An additional 1,966 and 1,970 cells, respectively, give direct exclusions at the benchmark on the part between the two late crossing bands. The two slope tests use respectively 4,004 and 5,005 cells in the early window, and 3,434 and 4,120 in the late window. The finite sums and outward comparisons are executable certificates. Lean checks the abstract deductions about opposite monotone branches, the last missing word and cell intersections. [Replay transcript](evidence/first-catalogue.txt) · [Lean deductions](DoublingPatterns.lean)

### Real and integer first starts

Let $`y_0`$ be the intersection of the two $`B`$ moves in the early window. Near it the moves are $`u\sqrt y+v/\sqrt y`$ and $`s\sqrt y+t/\sqrt y`$. The cell $`[46695.21,46695.22]`$ contains no prime-power knot, as verified at every possible boundary norm. Hence $`y_0=(t-v)/(u-s)`$. Monotonicity gives the unique early maximum $`\mu=M(y_0)`$. Direct evaluation of the two neighbouring integer starts gives $`M(46696)\lt M(46695)=D_B(46695)=\nu`$. The same monotonicity makes 46,695 the unique best early integer. A third calculation forms the two affine lines from ordinary high-precision prime counts at the cell endpoints and encloses

```math
46695.216167402621654621\lt y_0\lt 46695.216167402621654622.
```

```math
1.682449287659191709\lt \mu\lt 1.682449287659191710,\qquad
1.682222888957865575\lt \nu\lt 1.682222888957865576.
```

For $`0\le q\le\mu`$, let $`b(q)`$ be the inverse of the increasing first $`B`$ move on $`[45700,y_0]`$. For $`\nu\le q\le2`$, let $`a(q)`$ be the inverse of the increasing first $`A`$ move on $`[20288400,20309000]`$. The endpoint checks give the advertised domains. Write $`T(q)`$ for the infimum of real starting heights by which every word has appeared, allowing a different start for each word. Write $`N(q)`$ for the first such integer bound when every start must be an integer. The fifteen early witnesses and the four certified intervals give the formulas

```math
T(q)=\begin{cases}b(q)&0\le q\lt \mu,\\ a(q)&\mu\le q\le2,\end{cases}
\qquad
N(q)=\begin{cases}\lfloor b(q)\rfloor+1&0\le q\lt \nu,\\ \lfloor a(q)\rfloor+1&\nu\le q\le2.\end{cases}
```

Indeed for $`q\lt \nu`$ the early inverse is strictly below 46,695, and the next integer lies on the increasing side through that height. At $`q=\nu`$ strict inequality eliminates that last integer. From then on the next integer after the late crossing works, since the late start at margin 2 is below 20,308,191 and the whole late branch extends through 20,309,000. For real starts the early maximum itself is eliminated only at $`\mu`$. With weak positive-margin inequalities the early height $`y_0`$ attains $`\mu`$ and 46,695 attains $`\nu`$ in the integer problem. Completion heights are four times the respective starting heights. [Third peak calculation](word_peak.py) · [Monotonicity and exact boundary checks](word_cover.py)

In particular $`a(\nu)\in(20288500.19344,20288500.19345)`$ and $`a(\mu)\in(20288514.06421,20288514.06422)`$. At the integer threshold the first successful integer jumps from 46,695 to 20,288,501. The real infimum has exactly one jump on $`[0,2]`$, from left limit $`y_0`$ to $`a(\mu)`$. On the whole interval $`\nu\le q\lt \mu`$, of width greater than 0.00022639, one has $`N(q)/T(q)\ge N(\nu)/y_0\gt 434`$. At the benchmark the increasing first $`A`$ move controls the first crossing, and the first $`B`$, second $`A`$ and second $`B`$ moves are respectively greater than 20, 48 and 38 on $`[20306533.58,20306533.61]`$. The fifteen other patterns have already finished long before the all-gain pattern starts. [Certified roots, endpoints and integer tests](evidence/first-catalogue.txt)

For any algebraic $`q\in(0,2]`$, an algebraic crossing height would give an algebraic linear relation between odd prime logarithms and $`\log 2`$, with nonzero coefficient $`-q`$ on $`\log 2`$. Unique factorization gives rational independence of these logarithms. Baker’s theorem rules the relation out. For $`q=0`$, the prime 183,263 has degree one on field $`B`$, has nonzero coefficient in $`S_B(2y)`$ and is absent from $`S_B(y)`$ throughout the bracket $`45700\lt y\lt 45800`$. At the early peak the prime 380,207 has degree one on $`B`$ and appears only in the last term $`S_B(4y)`$ of the second difference $`S_B(4y)-2S_B(2y)+S_B(y)`$. Its weight is nonzero on the full small peak bracket. These two nonzero odd-prime coefficients prove that $`T(0)`$ and $`y_0`$ are transcendental as well. Finally $`\nu`$ is transcendental because it is a nonzero normalized move at the algebraic height 46,695. These deductions determine the arithmetic type of every real onset at an algebraic margin. [Baker’s theorem](https://doi.org/10.1112/S0025579300003971) · [Prime coefficient checks](word_cover.py)

### How many samples detect the early gain?

Suppose the sampling clock has an unknown initial phase. Fix $`\nu\le q\lt \mu`$ and write $`c(q)`$ for the height on the decreasing second $`B`$ move where it equals $`q`$. Before height 20,288,400 the only starts giving four strict gains above $`q`$ form the open interval $`(b(q),c(q))`$. This follows from the four certified intervals, the vanishing move below 6 and the strict maximum. The width $`w=c(q)-b(q)`$ is less than one since $`b(q)\ge46695`$ and $`c(q)\lt 46696`$. Its width in octaves is $`w_2=\log(c(q)/b(q))/\log2`$.

For additive sampling, choose $`m`$ offsets $`t_j`$ modulo one and sample at $`n+\theta+t_j`$ with common unknown $`\theta`$. For octave sampling, use $`2^{n+\theta+t_j}`$. In each case choose any integer $`n`$ that gives a positive start, taking all such observations. On the phase circle the successful offsets for a single clock form an open arc of length $`w`$ or $`w_2`$. If the chosen offsets have cyclic gaps $`g_1,\ldots,g_m`$, every initial phase is caught precisely when the largest gap is **strictly** smaller than the arc length. A gap equal to the length allows the two endpoints to meet at a missed phase. The gaps sum to one, so equal spacing minimizes their maximum and attains the minimum. For uniform initial phase the largest possible success probability is $`\min(1,mw)`$ for additive sampling, and $`\min(1,mw_2)`$ in octaves. The union bound gives the upper bound and equally spaced arcs attain the probability, including when the omitted endpoints have measure zero. Thus the smallest numbers guaranteeing success for **every** phase are

```math
m_1(q)=\lfloor 1/w\rfloor+1,\qquad m_2(q)=\lfloor 1/w_2\rfloor+1.
```

| Strict required margin | Minimum offsets per unit of height | Minimum offsets per octave |
| --- | ---: | ---: |
| Early integer maximum $`\nu`$ | 4 | 97,121 |
| $`\mu-10^{-6}`$ | 680 | 21,987,939 |

At $`q=\nu`$ the window starts at exactly 46,695 and ends between $`46695.333262779289365032`$ and $`46695.333262779289365033`$. Three equally spaced additive samples per unit miss **0.02116621319… percent** of initial phases. No arrangement of three does better. Four equally spaced samples guarantee a hit. The right endpoint is transcendental. If it were algebraic, equality of the second $`B`$ move there and the first $`B`$ move at 46,695 would give an algebraic linear relation in odd prime logarithms. The coefficient of 380,207 occurs only in the former move and is nonzero throughout $`[46695,46696]`$. Baker’s theorem also makes the width and the best miss probability for three samples transcendental. [Direct prime counts](clock_phases.py) · [Independent endpoint signs](clock_phases_recheck.py)

The early passage becomes hard to detect at a rate fixed by two actual prime sums. Put $`p=D_B'(y_0)\gt 0`$ and $`s=(d/dy)D_B(2y)|_{y=y_0}\lt 0`$. The inverse function theorem on the knot-free peak cell gives

```math
w=(p^{-1}-s^{-1})(\mu-q)+O((\mu-q)^2),\qquad
1472.01872576754168\lt p^{-1}-s^{-1}\lt 1472.01872576754169.
```

As the margin increases to the early maximum, $`(\mu-q)m_1(q)\to1/(p^{-1}-s^{-1})`$, while $`(\mu-q)m_2(q)\to y_0\log2/(p^{-1}-s^{-1})\in(21.98793864880736,21.98793864880737)`$. One calculation bisects the original weighted prime sums. Another constructs the two affine coefficients on each local cell by Euler criteria and direct logarithm sums, checks opposite signs at the rational bracket endpoints, and certifies the four integer minima. [Numerical evidence](evidence/clock-phases.txt) · [Lean gap-count deduction](DoublingPatterns.lean)

## Prime and square starting heights

Require all four coordinates $`D_A(y),D_B(y),D_A(2y),D_B(2y)`$ to exceed a common margin $`0\le q\le2`$. Write $`G`$ for the positive rational primes or the positive perfect squares, and $`g_G(t)`$ for the least member of $`G`$ **strictly** greater than $`t`$. Let $`b(q)`$ be the inverse of the first increasing $`B`$ move on the early band, used up to its intersection with the decreasing second $`B`$ move at $`y_0`$. Let $`a(q)`$ be the inverse of the increasing first $`A`$ move on the late band, now defined throughout $`[\sigma,2]`$. Define

```math
\sigma=D_B(216^2),\qquad \rho=D_B(46691),\qquad \nu=D_B(46695).
```

The critical margin $`q_G`$ equals $`\sigma`$ for squares and $`\rho`$ for primes. For **every** $`0\le q\le2`$, the least successful start in $`G`$ is

```math
g_G(b(q))\quad(0\le q\lt q_G),\qquad g_G(a(q))\quad(q_G\le q\le2).
```

Strict successors handle the case where an inverse is itself a grid point. The same least height completes the catalogue of sixteen possible two-step sign patterns when ordered by starting height. For the other fifteen patterns, the signed components all exceed 2 at explicit prime starts through 25,057 and square starts through 25,281. Each of these precedes $`b(0)`$. The order of bits and both witness lists are in the verifier. [Arithmetic witnesses and both covers](arithmetic_start_grid.py) · [Whole-interval transcript](evidence/arithmetic-starts.txt)

To prove the formula, first use the existing exclusion through 45,700. On $`[45700,46740]`$ both $`A`$ moves exceed 2 and the two $`B`$ moves are strictly increasing and strictly decreasing, respectively. The crossing lies strictly between 46,695.21 and 46,695.22. The last prime below it is 46,691. The decreasing move at the next prime is already below $`\rho`$. The last square below the crossing is 46,656 and the next is 47,089, beyond the early band. Thus the last successful grid point just below the critical margin lies on the increasing branch and gives the stated unique early maximum on its grid. For every smaller margin that grid point is eligible, and the first grid point after $`b(q)`$ lies on the increasing branch no later than it.

On the closed interval $`[46740,20284000]`$, at every height at least one coordinate is strictly below $`\sigma`$. The primary secant cover in log height has 38,519 cells. The independent affine cover in height has 33,411. On $`[20284000,20314000]`$ the first $`A`$ move increases from below $`\sigma`$ to above 2 and the other three moves exceed 2. The early and late monotonicity tests use 4,160 and 5,000 cells in the first count, 5,200 and 6,000 in the second. Both sieves evaluate only norms at most 325,024,000. The largest needed successor for primes or squares occurs at $`q=2`$ and is still inside the late band. This proves existence and first occurrence, including equality at each transition. The highest necessary square start 20,313,049 is below the upper endpoint by 951. [Cell bounds](#covering-all-real-heights) · [Prime sieve](word_sieve.cpp) · [Independent sieve](word_sieve_recheck.cpp)

| Margin | Late real crossing, rounded down to five decimals | First start on the indicated grid |
| --- | ---: | ---: |
| Square threshold $`\sigma`$ | 20,286,256.45311 | $`4505^2=20295025`$ |
| Prime threshold $`\rho`$ | 20,288,238.97111 | 20,288,251 |
| Integer threshold $`\nu`$ | 20,288,500.19344 | 20,288,501 |
| 1.9735 | 20,306,533.59533 | Prime 20,306,543 and square 20,313,049 |
| 2 | 20,308,190.27109 | Prime 20,308,199 and square 20,313,049 |

Each of the five-decimal lower endpoints in the table is within $`10^{-5}`$ of the exact crossing. The enclosure at $`\nu`$ is from the integer transition in the earlier calculation. At margin 1.9735 the second doubling finishes at 81,226,172 from the first prime and at 81,252,196 from the first square. [Late crossing enclosures](evidence/arithmetic-starts.txt) · [Real and integer enclosures](evidence/first-catalogue.txt)

For $`\sigma\le q\lt\rho`$ the first square exceeds $`434`$ times the first prime, since $`a(\sigma)\gt434\cdot46691`$ and the first prime is at most 46,691. For $`\rho\le q\lt\nu`$ the first prime exceeds $`434`$ times the first integer, since $`a(\rho)\gt434\cdot46695`$. For $`\nu\le q\lt\mu`$ the first integer exceeds $`434`$ times the real infimum, as proved above. The intervals are nonempty. The two sieves enclose all three grid thresholds, and an independent direct prime calculation resolves them to eighteen decimal places. [Three threshold values](evidence/arithmetic-starts.txt) · [Real infimum](#real-and-integer-first-starts)

## Which early thresholds can later readings recover?

Let $`\mathcal A`$ be the field of real algebraic numbers. For a rational prime $`p`$ good on at least one of the two fields, let $`n_i(p)=p^{f_i(p)}`$ on each field where it is good. Let $`N(p)`$ be the largest of these one or two norms and $`L_p=\log p/\log2`$. The five fixed omissions on each field apply throughout. For a real $`H\ge1`$ define the space of exact linear combinations of later moves by

```math
T_H=\operatorname{span}_{\mathcal A}
\{D_i(y):i=A,B,\ y\in\mathcal A,\ y\ge H\}.
```

**Theorem.** For every $`H\ge1`$ this space has basis exactly the $`L_p`$ with $`N(p)\gt H`$. Rational starting heights alone give the same space. Adding a freely chosen algebraic constant gives the direct sum $`\mathcal A\oplus T_H`$. For any finite vector of earlier algebraic-height readings, one can therefore decide which linear combinations of the vector lie in this enlarged space. Write each reading in the prime-logarithm basis and take the nullspace of the rows with $`N(p)\le H`$. This is a criterion for exact recovery by finite algebraic linear combinations. The coefficient field and the restriction to exact linear combinations are part of the statement. [Prime-log independence](https://doi.org/10.1112/S0025579300003971) · [Earlier finite coefficient test](https://github.com/leecrossley/theta-host-split-inert-race/blob/algebraic-heights-v1.7.1/results/algebraic-heights/PROOF.md#prime-logarithms-at-algebraic-heights)

At algebraic $`y`$, each move is a finite algebraic linear combination of the $`L_p`$. The weight is zero for $`n_i(p)\le y`$, including equality. Thus all terms on the left are in the proposed span. For the converse fix $`p`$ and a field with $`n=n_i(p)\gt H`$. Put $`F(y)=\sqrt y\,D_i(y)`$. The explicit term of $`p`$ is

```math
\frac{8e_{i,p}}{B_0}
\left(\frac{W(n/(2y))}{\sqrt2}-W(n/y)\right)L_p .
```

The function $`F`$ is piecewise affine. Its possible knots on a bounded interval are $`n'/k`$ with $`k=1,2,4,8`$ and $`n'`$ a good norm on this field. All such norms are odd and congruent to $`1`$ or $`-1`$ modulo 16. The only knot inside $`(n-1/2,n+1/2)`$ is $`n`$ itself. Indeed, any other integer norm is at distance at least one. For $`k=2`$ an odd numerator stays at least $`1/2`$ from $`n`$. For $`k=4`$ and $`k=8`$ the nearest possible numerator is at least 3 and 7 integer units respectively from $`kn`$. Distinct rational primes cannot have the same positive prime-power norm.

At $`y=n`$ only the last knot of the subtracted level contributes. Immediately below $`n`$ its tent is $`2-2y/n`$ and immediately above it vanishes. The right slope of $`F`$ minus the left slope is therefore $`-16e_{i,p}L_p/(B_0n)`$. Take rational $`0\lt h\le1/2`$ with $`n-h\ge H`$. Continuity, including possible knots at the sample endpoints, gives the exact reconstruction formula

```math
L_p=-\frac{B_0 n}{16e_{i,p}}
       \frac{F(n-h)-2F(n)+F(n+h)}{h}.
```

All three heights are rational. The factors that convert the moves to $`F`$ are algebraic. Unique factorisation makes the logarithms of 2 and the distinct odd primes rationally independent. Baker’s theorem gives algebraic linear independence, which proves the basis and the assertion about constants. For a fixed finite list of earlier readings, taking the cutoff at least as large as every last norm in their combined support leaves zero as the only recoverable element of their span. [Local signed trace](https://github.com/leecrossley/theta-host-split-inert-race/blob/algebraic-heights-v1.7.0/results/theta-race-short-shifts/RaceTrace.lean#L83) · [Three-point deduction](DoublingPatterns.lean) · [Baker](https://doi.org/10.1112/S0025579300003971)

The same description applies to two practical sampling schedules. On the unshifted half-integer grid, put $`t=\lceil2H\rceil/2`$. Its space is $`T_t`$, since $`n\gt t`$ allows the displayed three samples with $`h=1/2`$. More generally, suppose an unbounded increasing sequence of algebraic sample heights starts at $`t\ge H`$ and its successive gaps are at most $`1/4`$. Its space on the two fields is also $`T_t`$. For each $`n\gt t`$ there is a sampled $`s\in[n-1/2,n)`$ and two distinct samples in $`[n,n+1/2]`$. A sequence increasing to infinity is locally finite. The right affine line $`R`$ is obtained from the two samples on the right. Its continuation to $`s`$ satisfies $`(F(s)-R(s))/(n-s)=F'_+(n)-F'_-(n)`$. This proves the converse without requiring the clock to sample the knot. The cutoff uses the actual first sample, so a prime whose last norm equals that sample has already disappeared. [Exact off-grid deduction](DoublingPatterns.lean)

### Three rank losses

Consider the three early thresholds in the preceding table, in square, prime, integer order. Scale them by the algebraic factors $`\sqrt{46656},\sqrt{46691},\sqrt{46695}`$, respectively. Each threshold is a first $`B`$ move. For a contributing rational prime with $`N(p)\le46993`$, its $`B`$ norm satisfies $`46656\lt n\le46993\lt2\cdot46656`$. The forward level $`S_B(2y)`$ has zero weight. After removing a common nonzero factor depending only on the prime, its coefficient row at $`y=46656,46691,46695`$ is $`(\max(n-y,0))_y`$. The two exact enumerations produce every such row in order of its final norm across **both** fields.

| Final norm $`N(p)`$ | Rational prime | Three scaled coefficients |
| ---: | ---: | --- |
| 46,687 | 46,687 | $`(31,0,0)`$ |
| 46,703 | 46,703 | $`(47,12,8)`$ |
| 46,993 | 46,993 | $`(337,302,298)`$ |

The first calculation sieves rational primes, obtains each residue degree by Euler’s criterion, and takes the larger eligible norm. The second factors each candidate norm in the support range, tests square-residue tables, and discards primes still present later on the other field. It also checks the rational finite differences at all three displayed knots. These rows have determinant 35,960. Therefore $`1,\sigma,\rho,\nu`$ are linearly independent over $`\mathcal A`$. [Sieve and determinant](recover_thresholds.py) · [Factor-by-norm verification](recover_thresholds_recheck.py) · [Transcript](evidence/arithmetic-starts.txt)

For $`1\le H\lt46687`$, all three thresholds belong to $`T_H`$. For $`46687\le H\lt46703`$, the recoverable subspace of their span has basis $`\rho,\nu`$. For $`46703\le H\lt46993`$, it is the single line generated by

```math
3\sqrt{46695}\,\nu-2\sqrt{46691}\,\rho.
```

Indeed the first two rows annihilate the coefficient vector $`(0,-2,3)`$ and the third pairs with it to give 290. From $`H=46993`$ onward the intersection of their span with $`\mathcal A\oplus T_H`$ is zero. This statement permits both fields, arbitrarily many later heights, and an algebraic offset. In particular the three **late real crossings** $`a(\sigma),a(\rho),a(\nu)`$ are all transcendental. Each exceeds 46,993. If one were algebraic its first $`A`$ move would be both its nonzero early threshold and a member of $`T_{46993}`$. The endpoint proof applies even though the three margins are themselves transcendental. [Exact row deductions](DoublingPatterns.lean) · [Crossing bounds](evidence/arithmetic-starts.txt)

### When both fields are needed

Restricting observations to one field changes the list of accessible primes to those with $`n_i(p)\gt H`$ on that field. The same last-knot proof gives the converse on each host. Of the 3,327 rational primes occurring in the three early thresholds, all are good on both fields. The following entries are the first heights at which the dimension of the recoverable threshold space drops to the indicated value. Each drop takes effect at equality.

| Available later moves | Dimension two begins | Dimension one begins | Dimension zero begins |
| --- | ---: | ---: | ---: |
| 547/727 alone | 223 | 271 | 46,687 |
| 571/599 alone | 46,687 | 46,703 | 46,769 |
| Both fields | 46,687 | 46,703 | 46,993 |

For the first row of this table the only active primes with $`A`$ norm at most 46,687 are 223, 271, 353, 367, 463 and 46,687, all of degree one on $`A`$. Their $`B`$ norms are 49,729, 73,441, 124,609, 134,689, 214,369 and 46,687. Put $`w(n,y)=\max(0,\min(2n-2y,4y-n))`$. The scaled coefficient up to a nonzero common factor is $`-2w(n,y)+\sqrt2\,w(n,2y)`$. The first two rows are rational and independent. The coefficient vector $`(4,-39,35)`$ annihilates both rational parts of each of the first five rows and has nonzero pairing with the sixth. For $`B`$ alone there is a further row at 46,769 before the joint loss at 46,993. At $`n=46769`$ the previously scaled integer row is $`(113,78,74)`$. It pairs with $`(0,-2,3)`$ to give 66. This prime still occurs on $`A`$ at norm $`46769^2=2187339361`$. Two full enumerations in opposite orders check all 3,327 primes. [Single-field coefficient tests](single_field_recovery.py) · [Independent norm enumeration](recover_thresholds_recheck.py)

Fix **any** $`46769\le H\lt46993`$ and an upper observation height $`U\ge H`$. On the closed interval $`[H,U]`$ the nonzero surviving threshold combination $`3\sqrt{46695}\,\nu-2\sqrt{46691}\,\rho`$ is recoverable from the two fields precisely when

```math
U\gt \frac{46769^2}{8}=273417420+\frac18 .
```

Indeed its coefficient on $`L_{46769}`$ is nonzero. Field $`B`$ gives this prime zero weight throughout $`y\ge H`$. Its first possible move on field $`A`$ begins at $`y=46769^2/8`$, where the weight is still zero. Baker’s theorem therefore forces the strict lower bound, even with an algebraic offset. This proves failure at the exact upper endpoint.

For the converse, among the participating primes with $`B`$ norms at most $`H`$, the only ones are 46,687, 46,703 and 46,769. The first two have zero coefficient in the surviving combination. Every other required prime can be recovered from $`B`$ near its last knot, using starts at least $`H`$ and below 373,560. To recover $`L_{46769}`$ on $`A`$, take three rational samples around the first knot $`z=46769^2/8`$. All knots on this field lie on the one-eighth grid. At $`z`$ only this prime contributes a slope jump, since its numerator in lowest terms is odd and distinct primes have distinct norms. For $`n=46769^2`$ the outer tent of the forward level is zero to the left of $`z`$ and equals $`8y/n-1`$ immediately to its right. Its jump in the slope of $`\sqrt y\,D_A(y)`$ is $`(8e_{A,46769}/B_0)(4\sqrt2/n)L_{46769}`$. For every $`U\gt z`$ a positive rational step smaller than $`1/8`$, $`z-H`$ and $`U-z`$ gives the exact three-point extraction inside $`[H,U]`$. The single-field ranks show that both fields are necessary even if the upper endpoint is removed. [Exact numerator and entrance slope](recover_thresholds_recheck.py) · [Sharp endpoint certificate](evidence/arithmetic-starts.txt)

Three $`A`$ readings are the **minimum** throughout the following interval of permitted upper endpoints, even if one is allowed arbitrarily late $`B`$ readings:

```math
273417420+\frac18\lt U\le539330952+\frac14 .
```

To see this, use the primes $`46447,46559`$. Both have degree one on $`B`$ and degree two on $`A`$. Each is below 46,656, so it has zero coefficient in all three thresholds, and contributes to no later $`B`$ move. Let $`n_1=46447^2`$ and $`n_2=46559^2`$. Up to $`U\le n_1/4`$, the coefficient of each of these two primes in $`\sqrt yD_A(y)`$ is a fixed nonzero multiple of $`\max(0,y-n_j/8)`$. Any successful reconstruction must include an $`A`$ reading at $`y_2\gt46769^2/8`$, where both these unwanted coefficients are nonzero. To cancel both with one other $`A`$ reading, its height $`y_1`$ must exceed $`n_2/8`$ as well. The determinant of the two background rows, after removing nonzero factors, is $`(y_1-y_2)(n_1/8-n_2/8)`$. Distinct heights give a nonzero determinant. Duplicate readings collapse into a single reading and cancelling both unwanted coefficients then cancels the desired one too. Three $`A`$ readings suffice by the entrance-knot construction. [Two-prime calculation](single_field_recovery.py) · [Independent residue check](recover_thresholds_recheck.py)

## The price of limiting observations on one field

Throughout this section fix any real $`46769\le H\lt46993`$. The target is the nonzero surviving threshold combination from the preceding theorem. Samples must be actual moves at algebraic heights in $`[H,U]`$. A reconstruction is a finite linear combination with real algebraic coefficients and an optional algebraic constant. Let $`m`$ be the maximum number of $`A`$ samples. The number of $`B`$ samples may be any finite number. Write

```math
n=46769^2,\qquad n_1=46447^2,\qquad n_2=46559^2,\qquad
P=17258590369,\quad P_1=17341923439.
```

| Maximum number of $`A`$ samples | Required bound on their latest height if $`B`$ has no upper bound | Exact condition on a common endpoint $`U`$ for both fields |
| ---: | --- | --- |
| 0 | Impossible | Impossible |
| 1 | At least $`n_2=2167740481`$, attained there | $`U\gt P_1^2/8`$ |
| 2 | Strictly above $`n_1=2157323809`$ | $`U\gt P^2/8`$ |
| 3 or more | Strictly above $`n/8=273417420.125`$ | $`U\gt n/8`$ |

For the second column, each entry gives a condition on the permitted upper height for $`A`$ alone. In the third, equality fails in all three finite rows. The exact one- and two-sample common infima are **37,592,788,570,517,198,340.125** and **37,232,367,690,617,444,520.125**. Adding the third $`A`$ reading cuts the joint infimum by a factor between **136,174,087,494** and **136,174,087,495**. The proof uses exact prime coefficients, rational inequalities, and Baker’s linear independence theorem. There is no assumption on zeros. [Primary computation](observation_budget.py) · [Independent computation](observation_budget_recheck.py) · [Transcript](evidence/observation-budget.txt)

### The individual prime weights

The preceding basis theorem allows us to compare the coefficient of each prime logarithm. Use $`F_i(y)=\sqrt yD_i(y)`$, which has the same linear-recovery power as the move at algebraic $`y`$. Its coefficient at a good norm $`m_0`$ is a fixed nonzero algebraic multiple of $`\psi(y/m_0)`$, with

```math
\psi(x)=W(1/(2x))/\sqrt2-W(1/x).
```

Its support with endpoints removed is $`(1/8,1)`$ except for the internal zero $`1/(2\sqrt2)`$. On the three closed pieces from $`1/8`$ to $`1/4`$, from $`1/4`$ to $`1/2`$, and from $`1/2`$ to 1, the formulas are respectively $`(8x-1)/\sqrt2`$, $`1+\sqrt2-(4+2\sqrt2)x`$ and $`2x-2`$. The field-specific sign and the constant of normalization have no effect on the zero locations or proportionality tests.

For a prime good on both hosts, different residue degrees can arise only when the prime is congruent to $`\pm1\bmod16`$. The two degrees then belong to $`\{1,2\}`$. Indeed, in every other odd residue class the order modulo 16 up to sign is 2 or 4, already even, so the extra quadratic condition does not change it. The primes with $`B`$ degree one and $`A`$ degree two are the ones whose weights persist on $`A`$ after their $`B`$ norms expire. Among such primes through 46,992, the largest four are 46,337, 46,447, 46,559 and 46,769. The first three are smaller than all original threshold heights and have zero coefficient in the target. Prime 46,769 has nonzero coefficient in it and must be supplied by $`A`$ at some $`y\in(n/8,n)`$.

The primes omitted on $`B`$ cannot add an unwanted coefficient in the upper constructions near $`n_1,n_2`$. On $`A`$ the two $`B`$ conductors 571 and 599 have norms 106,302,733,681 and 358,801 respectively. The two additional $`B`$ omissions exceed $`8n`$ even as rational primes, and 2 is omitted on both fields. After the other $`B`$-inaccessible coefficients have been cancelled, an unbounded $`B`$ tail can reproduce every remaining finite linear combination by the last-knot formula. [Local splitting rule](#the-scores-and-the-arithmetic-input) · [Prime and omission checks](observation_budget.py)

### One or two observations when the other field may be observed arbitrarily late

A single nonzero $`A`$ reading must include the target prime. At every $`n/8\lt y\lt n_1`$, the unwanted norms $`n_1,n_2`$ are both in the interior of their support. Their internal zeros occur at different heights, so they cannot both vanish. For $`n_1\le y\lt n_2`$ the second norm has nonzero weight on its last affine piece. It follows that one $`A`$ observation must be at least $`n_2`$. At $`y=n_2`$, the target is present and every other $`B`$-inaccessible prime has zero weight. Choosing an algebraic multiplier matches the target prime coefficient. The rest can be supplied by $`B`$. This also works at every algebraic $`y\in[n_2,n)`$.

Now suppose that at most two $`A`$ readings, at distinct heights $`s\lt t`$, recover the target. Equal heights can first be combined. To exclude $`t\le n_1`$, we may assume $`t\gt n/8`$ and both weights are nonzero. At least one of $`n_1,n_2`$ has nonzero coefficient at $`t`$, including when $`t=n_1`$. Cancellation with the observation at $`s`$ forces $`s\gt n_1/8`$. Let $`\mathcal N`$ be the 306 norms $`p^2`$ with $`p\lt46656`$, $`B`$ degree one and $`A`$ degree two. These primes have zero target coefficient and have disappeared from $`B`$ at every allowed height. Let $`b`$ be the second largest norm of $`\mathcal N`$ at most $`t`$. Then $`s\ge b`$. Otherwise at least two different expired norms $`a\le t`$ would lie strictly above $`s`$. Both have $`a/8\le t/8\le n_1/8\lt s`$, so both would have to vanish at their distinct internal zeros at that same $`s`$.

For every $`n/8\lt t\lt n_1`$, two norms $`a_1,a_2`$ of $`\mathcal N`$ satisfy $`t\lt a_1\lt a_2\lt2b`$. This is a finite rational inequality. In ascending order put $`a_j'`$ for the sequence of 306 norms and take $`a_k'\le t\lt a_{k+1}'`$. Use $`b=a_{k-1}'`$ and the two subsequent norms $`a_{k+1}',a_{k+2}'`$. The maximum of $`a_{k+2}'/a_{k-1}'`$ on the 188 required intervals is exactly $`377952481/325766401\lt2`$. Thus at both $`s`$ and $`t`$ the two unwanted prime weights are on the last affine piece, proportional to $`y-a_1`$ and $`y-a_2`$. No two nonzero sample weights can cancel both. Their determinant, up to nonzero factors, is $`(s-t)(a_1-a_2)`$. At $`t=n_1`$, the second largest expired norm is $`46337^2\gt n_1/2`$. If $`s\lt n_1`$, the unwanted norm $`n_1`$ is then nonzero at $`s`$ and zero at $`t`$. This proves the strict lower bound in the two-sample row. [Two interval enumerations](evidence/observation-budget.txt) · [Two-sample determinant](DoublingPatterns.lean)

For every upper bound above $`n_1`$, choose $`s=n_1`$ and a distinct rational $`t`$ sufficiently close to its right, with $`t\lt n_2`$. Only the unwanted norm $`n_2`$ and the target norm $`n`$ occur from the primes inaccessible on $`B`$. Both weights are on their last affine piece. The two affine functions have different zeros, so a linear combination can eliminate the former and match the latter. Equivalently, extrapolate the two affine values to $`n_2`$ using coefficients $`(t-n_2)/(t-s)`$ and $`(n_2-s)/(t-s)`$. Apply the final algebraic scaling for the target prime. All remaining logarithms belong to the $`B`$ tail. [Affine extrapolation](DoublingPatterns.lean)

### What the other field must wait for

A large rational prime that has degree one on $`A`$ and degree two on $`B`$ appears on $`B`$ only after $`p^2/8`$. If its logarithm appears spuriously in the proposed $`A`$ contribution, the $`B`$ observations must reach strictly past this endpoint to remove it. The endpoint itself has zero weight. Two interval sieves give $`P_1`$ and $`P`$ as the largest such primes below $`8n_2`$ and $`8n_1`$, respectively. Their distances below the two even upper bounds are 409 and 103.

For one $`A`$ reading the possible heights are contained in $`[n_2,n)`$. The prime $`P_1`$ satisfies $`4n\lt P_1\lt8n_2`$, so it has nonzero weight in every such reading and zero coefficient in the early target. This forces $`U\gt P_1^2/8`$.

For the two-reading lower bound, take the additional degree-one-on-$`A`$, degree-two-on-$`B`$ primes

```math
Q=17258590801,\qquad R_1=17498715103,\qquad R_2=17498715233.
```

They satisfy $`4n\lt P\lt8n_1\lt Q\lt8n_2\lt8n\lt R_1\lt R_2\lt2P`$. Suppose $`U\le P^2/8`$. On $`B`$, none of $`P,Q,R_1,R_2`$ can contribute. All are absent from the target. If only one $`A`$ term remains, the one-reading bound already precludes recovery since $`P_1\gt P`$. Otherwise use distinct $`s\lt t`$ with nonzero coefficients, and $`t\gt n_1`$ by the preceding argument.

1. Suppose $`n_1\lt t\le Q/8`$. The unwanted norm $`n_2`$ is present on its last affine branch at $`t`$, so $`s\gt n_2/8`$. Both $`46337^2`$ and $`n_1`$ are extinct at $`t`$. If $`s\lt46337^2`$, the two distinct internal zeros could not remove both there. Hence $`s\ge46337^2\gt n_1/2`$, after which the need to remove $`n_1`$ forces $`s\ge n_1`$. The two unwanted prime weights from $`n_2`$ and $`P`$ are now proportional to $`y-n_2`$ and $`y-P/8`$ at both sampled heights. These have different zeros, and two distinct observations cannot cancel both.
2. Suppose $`Q/8\lt t\le P/4`$. Both $`P`$ and $`Q`$ are nonzero on their first affine branches at $`t`$. Both must also be present at $`s`$, which forces $`s\gt Q/8`$. Both then lie on the first branch at both heights. The two different zeros $`P/8,Q/8`$ again give a nonzero determinant.
3. Suppose $`t\gt P/4`$. Then $`t\gt n`$, so $`s`$ alone must supply the target prime. The two unwanted norms $`n_1,n_2`$ vanish at $`t`$. As in the one-reading lower bound, $`n_2\le s\lt n`$. The prime $`P`$ is present at $`s`$, so $`t\lt P`$ is necessary to cancel it. But throughout $`P/4\lt t\lt P`$ each of $`R_1,R_2`$ has $`R_j/8\lt t\lt R_j`$. Their distinct internal zeros cannot coincide, and neither prime occurs at $`s`$ because $`R_j\gt8n`$. At least one unwanted coefficient remains.

This exhausts the heights and forces $`U\gt P^2/8`$. Baker’s theorem justifies each prime-by-prime necessity even when arbitrary algebraic coefficients and a constant are allowed. [Exact prime searches](observation_budget.py) · [Independent segmented sieves](observation_budget_recheck.py)

Both common upper bounds are attained as **infima** by finite reconstructions. For one $`A`$ observation use $`y=n_2`$. For two take $`s=n_1,t=n_1+1/16`$ and the extrapolation above. No integer norm larger than $`8n_1`$ can enter either of these last two readings. After matching the target and eliminating the expired primes, the residual to be supplied by $`B`$ has only finitely many nonzero coefficients, all with $`B`$ norms exceeding $`H`$. In the one-reading case these norms are at most $`P_1^2`$. In the two-reading case they are at most $`P^2`$. Indeed all contributing $`A`$ norms are at most $`8n_2`$ or $`8n_1`$, respectively. An increase of degree from $`A`$ to $`B`$ is possible only from one to two, and the two largest eligible rational primes have just been certified. The original target contributes only $`B`$ norms below 373,560.

For any required $`B`$ norm $`b\le8H`$, its last knot $`b\gt H`$ lies well inside the allowed upper bound. If $`b\gt8H`$, its first knot $`b/8`$ is strictly past $`H`$ and strictly below any permitted $`U`$ in the corresponding row. Three rational observations placed sufficiently close on either side of that knot extract its logarithm, as in the preceding section. Doing this finitely many times supplies the entire residual. The three-reading row and the zero-reading impossibility follow from the previous paired-field theorem. [Three-point recovery](#which-early-thresholds-can-later-readings-recover) · [Exact common endpoints](evidence/observation-budget.txt)

## Opposite runs on the unshifted dyadic grid

An observation starts at exactly $`y=2^n`$. Its four changes, in order, are $`D_A(y),D_B(y),D_A(2y),D_B(2y)`$. A word succeeds at normalized margin $`q\ge0`$ if each change times its prescribed sign exceeds $`q`$. Both doublings finish at $`4y`$. We consider all integer exponents and require strict inequalities even at zero. For $`n\lt0`$ every contributing norm would be below 16. The two races have no such norms. [Scores and trace interpretation](#the-scores-and-the-arithmetic-input)

**Theorem.** Define $`\alpha=D_B(2^{12})`$ and $`\beta=-D_B(2^{35})`$. For every $`0\le q\lt\alpha`$, word `-+-+` first occurs at $`2^{12}`$ and `+-+-` first occurs at $`2^{34}`$. For $`\alpha\le q\lt\beta`$, the first start for `+-+-` is still $`2^{34}`$, and `-+-+` has no start through $`2^{34}`$. For $`q\ge\beta`$ neither orientation has a start through $`2^{34}`$. In particular the starting heights differ by $`2^{22}=4{,}194{,}304`$ at $`q=3947/2000`$. Their respective finishing heights are 16,384 and 68,719,476,736. The enclosures

```math
3.41437510\lt\alpha\lt3.41437512,\qquad
6.75643706\lt\beta\lt6.75643708
```

are obtained by separately enclosing the limiting component. At exponent 12 the four normalized moves, rounded to two decimals, are $`(-17.92,3.41,-4.33,53.23)`$. At exponent 34 they are $`(15.50,-30.07,47.23,-6.76)`$. [Directed prime bins](evidence/dyadic-primary-571.txt) · [Record thresholds](evidence/dyadic-records.tsv)

Here is the first exponent at $`q=3947/2000`$ for each word with a start through exponent 34. The last two words in the table have no strict start in this range at any nonnegative margin. The prime counts make no assertion about their first subsequent occurrences.

| Word | First exponent | Word | First exponent |
| --- | ---: | --- | ---: |
| `+--+` | 3 | `++--` | 6 |
| `--++` | 7 | `++-+` | 8 |
| `---+` | 11 | `-+-+` | 12 |
| `+-++` | 15 | `+++-` | 16 |
| `-+--` | 18 | `-++-` | 21 |
| `+---` | 22 | `-+++` | 24 |
| `----` | 26 | `+-+-` | 34 |
| `--+-` | none through 34 | `++++` | none through 34 |

The table of 23 records gives a complete first-passage rule through exponent 34. For each word read its rows in increasing start exponent. A row reports the first exponent for margins from the preceding record value **inclusive** (zero for the first row) up to the current record value **exclusive**. At or above the final value there is no start in the range. The exact value in a row is the absolute value of the identified field move at the identified limiting exponent. The last column gives its floor after multiplication by ten million. At $`q=0`$ fourteen words occur. All fourteen still occur when $`0\le q\lt\alpha`$. Exactly thirteen remain at $`q=\alpha`$, since all other final record values exceed $`\alpha`$. [Complete table](evidence/dyadic-records.tsv)

### Counting prime ideals on the grid

Set $`A_j=P(2^{j-1},2^j)`$ and $`R_j=Q(2^{j-1},2^j)`$ with the signed logarithmic and reciprocal logarithmic sums defined above. Each rational prime supplies its single eligible ideal norm and the trace determines its sign. At the largest start the last score in the second move uses norms through $`2^{38}`$. The first exceptional auxiliary prime is 343,057,911,569, beyond this cutoff. Both calculations omit the primes above 2 and the respective two conductors. They retain all eligible degree-one, degree-two and degree-four contributions. [Splitting rule](#exact-prime-sums-and-independent-counts) · [Progression sieve](word_sieve.cpp) · [Odd-integer sieve](word_sieve_recheck.cpp)

One calculation forms the three scores and differences them, using directed 58-bit logarithms and 96-bit reciprocal logarithms. The other uses 62 and 102 bits. It reconstructs the move directly from three disjoint octaves. For $`y=2^j`$ multiply the move by $`B_0a\sqrt y/8`$ to obtain

```math
-2A_{j+1}+2yR_{j+1}
 +(1+\sqrt2)A_{j+2}-(4+2\sqrt2)yR_{j+2}
 -A_{j+3}/\sqrt2+4\sqrt2\,yR_{j+3}.
```

The identity follows by substituting both branches of the tent and collecting the coefficients of $`1`$ and $`1/n`$ on each octave. The second calculation encloses $`\sqrt2`$ by an integer square root and $`\log 2`$ by 90 positive terms of $`2\operatorname{arctanh}(1/3)`$ with a geometric tail bound. From this point every endpoint is rational. The signed bounds from the two sieves overlap in every bin, after accounting for their different scales. [Score differences](dyadic_records.py) · [Octave weights and rational intervals](dyadic_records_recheck.py)

All fixed-point products and accumulators fit their integer widths. Below $`2^{38}`$ a rational prime logarithm is below 27. A one-unit allowance in real units bounds the directed rounding. A raw log sum before scaling is consequently below $`28\cdot2^{38}`$ in absolute value. For reciprocal weights, distinct rational primes give distinct norms and $`\sum_{n\le 2^{38}}1/n\le 1+38\log 2`$. After scaling by at most $`2^{102}`$ and allowing one final rounding unit per term the total stays below $`2^{113}`$. Both signed accumulators have 128 bits. The log-normalization shifts in the sieve use inputs below $`2^{38}`$ and at most 15 mantissa bits. [Two sieve implementations](word_sieve.cpp) · [Independent integer bounds](word_sieve_recheck.cpp)

For the first three exponents at least one move is exactly zero. At every other exponent the sign and the unique smallest absolute component are decided by disjoint bounds. The first verifier keeps strict records. The independent verifier checks all 560 word/exponent pairs against the proposed preceding record, including exact ties where adjacent words use the same physical reading. This proves both the successes and every earlier exclusion for each interval of margins. The kernel checks the implication from a strict hit and the earlier weak upper bounds to the least successful exponent. [Replay transcript](evidence/dyadic-starts.txt) · [Kernel deductions](DoublingPatterns.lean)

### Sixteen independent transition margins

There are sixteen different limiting readings among the 23 records. In increasing exponent they occur at $`3,5,6,8,10,12,14,15,17,19,21,23,24,26,31,35`$. **These sixteen positive margins and 1 are linearly independent over the real algebraic numbers.** In particular every transition value is transcendental. Margins shared across different words occur at exactly the repeated readings specified in the table. [Pivot primes and fields](evidence/dyadic-pivots.tsv)

Multiply a proposed algebraic linear relation by $`\log 2`$. Each normalized reading at an algebraic height now contributes a finite linear combination of logarithms of odd rational primes with algebraic coefficients. Its absolute value changes just the overall sign. The coefficient of $`\log 2`$ comes from the constant in the relation. Baker’s theorem and unique factorization give linear independence of all these prime logarithms over the algebraic numbers. [Baker’s theorem](https://doi.org/10.1112/S0025579300003971) · [Prime-logarithm formulation for these races](https://github.com/leecrossley/theta-host-split-inert-race/blob/algebraic-heights-v1.7.1/results/algebraic-heights/PROOF.md#prime-logarithms-at-algebraic-heights)

At the earliest exponent the rational prime 7 has residue degree two on 571/599 and its norm 49 lies strictly between 32 and 64. At every later entry with exponent $`j`$ the listed prime has degree one on the specified field and satisfies $`2^{j+2}\lt p\lt2^{j+3}`$. The earlier exponents use norms at most $`2^{j+2}`$. In the last octave only the second score in the move is active, with strictly positive tent and nonzero signed weight. Thus each pivot occurs with nonzero coefficient in its own reading and zero coefficient in all preceding readings. Taking the largest index with nonzero scalar in a proposed relation gives a contradiction at its pivot prime. The independent program checks primality by trial division, residue degrees by quadratic residues and all octave boundaries by integers. The kernel checks the triangular argument for arbitrary finite coefficient arrays. [Pivot verification](dyadic_records_recheck.py) · [Triangular lemma](DoublingPatterns.lean)

## Profiles for the sixteen characters

The compact test and its Fourier transform are

```math
g(h)=\frac4{B_0}\sinh\!\left(\frac{(a-|h|)_+}{2}\right),\qquad
G(t)=\frac{3-2r\cos(at)}{c(t^2+1/4)}.
```

Let $`V_i(x)=\sum_{\gamma} G(\gamma)e^{i\gamma x}`$ over all ordinates of the eight primitive factors for $`i`$, with multiplicity. The family is closed under conjugation, so $`V_i`$ is real and even. The coefficients are positive. The standard zero count and the quadratic decay give absolute uniform convergence. In particular $`|V_i(x)|\le M_i=V_i(0)`$. For $`z\ge1`$ set $`x=\log z`$, $`s=e^{-|h|}`$ in $`T(h)`$ below, and define

```math
T(h)=\frac{(1+s)\log(1+s)+(1-s)\log(1-s)}{2\sqrt s},\quad T(0)=a,\qquad
C_i=\sum_\chi\log(q_\chi/\pi)-8(\gamma_{\rm E}+a).
```

For a rational prime power $`n=p^j`$, let $`t_i(n)=\sum_\chi\chi(p)^j`$ using the primitive local values, zero when the character is ramified. The compact-test explicit formula gives a finite arithmetic evaluation,

```math
V_i(x)=C_i g(x)-\frac8c\{r[T(x+a)+T(x-a)]-3T(x)\}
 -\sum_{n=p^j}t_i(n)\frac{\log p}{\sqrt n}\,g(x-\log n).
```

Only $`z/2\lt n\lt 2z`$ contributes. The reflected prime sum vanishes on $`x\ge0`$, because its argument is at or beyond the kernel boundary. Evenness handles negative $`x`$. To obtain the gamma terms, expand the odd gamma-factor kernel as a sum of exponentials of rates $`2k+3/2`$. Between its three knots, $`g''=g/4`$. Convolution against each exponential leaves three translates, with weights proportional to $`r,-3,r`$. Their scalar sum is $`T(h)=\sum_{k\ge0}s^{2k+3/2}/[(2k+1)(2k+2)]`$. The remaining multiple of $`g`$ gives $`C_i`$. Only the conductors and prime traces change between the two hosts. [Compact-test explicit formula and normalization](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-curvature-euler-factors/PROOF.md#remove-the-archimedean-term) · [Odd gamma convolution at fixed width](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-width-cusp/PROOF.md)

For clarity, the passage back to ordinary primes only needs the classical prime ideal theorem for the unsigned square term. At $`x=\log(2y)`$ one has $`g(x-\log n)/\sqrt n=W(n/y)/(B_0\sqrt y)`$. Expanding the primitive Euler factors over the base field gives the signed prime powers with coefficient $`\epsilon_w^j\log Nw`$ away from the finitely many exceptional primes. At each exceptional place the powers in the compact window have total normalized size $`O(1/\sqrt y)`$ for large $`y`$. The archimedean terms tend to zero as $`x\to\infty`$. To pass from the negative prime-power sum to the negative prime sum, restore the positive unsigned squares. The prime ideal theorem and partial summation give their limit as $`B_0^{-1}\int_1^2 W(t^2)\,dt=1`$. The integral follows by integrating the two branches, or from the weight's Mellin value at one half. The total absolute contribution from powers of order at least three is $`O(y^{-1/6}\log y)`$ by summing the prime-ideal Chebyshev bound up to $`(4y)^{1/j}`$ over $`3\le j\le\log(4y)/\log 2`$. Hence, with $`P_i(x)=(V_i(x+a)-V_i(x))/a`$,

```math
S_i(y)=1+V_i(\log(2y))+o(1),\qquad
D_i(y)=P_i(\log(2y))+o(1).
```

Both errors tend to zero on the entire positive tail. The prime-ideal convention for the finite examples permits all unramified rational primes. [Fixed-width prime-only limit](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-width-cusp/PROOF.md#the-prime-only-limit-at-a-fixed-width)

## Sixteen pairs of consecutive moves

The [witness file](evidence/witnesses.tsv) has five integer columns. They are a code $`w=0,\ldots,15`$, the numerator and denominator of a positive ratio $`z_w`$, then the numerator and denominator of a finite start $`y_w`$. Write the code in four binary digits, with leading zeros. In order the digits prescribe field A then B on the first doubling, followed by A then B on the next. A one prescribes gain and a zero prescribes loss. Write the corresponding signs as $`s_{w,i,j}`$ for $`j=0,1`$.

Two independent interval evaluations prove that every specified sign satisfies

```math
s_{w,i,j}D_i(2^j y_w)\gt 1.6,\qquad
s_{w,i,j}P_i(\log z_w+j a)\gt 9.1 .
```

All endpoints are exact rationals before transcendental evaluation. The largest finite norm is 747,539 and the explicit formula uses prime powers only through 54,638. The direct calculation builds quadratic residue tables. For odd $`p`$, the sum over the quartic twists vanishes unless $`p^j\equiv\pm1\bmod16`$, when it is four times the sum of the two quadratic values to exponent $`j`$. At 2 only the two odd-conductor characters survive. For finite prime ideals the base residue degree $`f`$ is the least exponent in $`1,2,4`$ that satisfies both $`p^f\equiv\pm1\bmod16`$ and $`(\chi_u(p)\chi_v(p))^f=1`$. There are $`8/f`$ primes of this norm and their common sign is $`\chi_u(p)^f`$. Prefix sums of the coefficients and of the coefficients divided by the norm evaluate the two branches exactly up to the directed logarithm enclosures. [Primary verifier](patterns.py)

The independent verifier factors the polynomials $`X^4-4X^2+2`$, $`Y^2-uv`$ and $`Z^2+u`$ over each good residue field. Least common multiples of the local degrees give the residue degree in the base and in the extension. It sums the weights directly at the rational starts. For the spectral calculation it evaluates the sixteen primitive Conrey characters and the gamma series. The remainder after 128 terms of that series is at most $`s^{257+1/2}/[(257)(258)(1-s^2)]`$ for $`0\lt s\lt 1`$. The value at one is handled separately. Both calculations also enclose $`0\lt M_i\lt 180`$. [Independent verifier](independent.py) · [Expected transcript](evidence/check.txt)

## Every six-step word on the dyadic grid

For a word of length $`m`$ write each joint direction as a base-four digit, with its two binary bits giving the signs for fields A and B in that order. The first step is the most significant digit. The [dyadic table](evidence/dyadic-starts.tsv) contains a length, code and signed integer $`k`$ for **every** code of lengths four, five and six. The log centre for the first move is $`k a`$. The associated prime start after an integer translation $`h a`$ is $`y=2^{h+k-1}`$. A negative $`k`$ is a legitimate profile sample by evenness of $`V_i`$. All translated prime starts are positive integer powers of two for sufficiently large $`h`$.

For each code the certificate proves a strict signed lower bound for every $`P_i((k+j)a)`$, $`0\le j\lt m`$. The respective lower bounds are 14.1, 9.1 and 4.1 for lengths four, five and six. There are 5,376 tabulated witnesses and 61,440 signed component inequalities. The zeros and initial exponents can be chosen by floating-point search. The two replays below operate on the exact published integers. [Starting exponents](evidence/dyadic-starts.tsv) · [Rational thresholds](evidence/dyadic-margins.tsv)

The [bracket file](evidence/dyadic-zeros.tsv) gives the host index, primitive conductor, Conrey character label and the two endpoints in units of $`10^{-14}`$. The checker verifies all sixteen primitive odd characters from the conductor frame, all 28,000 Hardy signs, and disjointness within each factor. Select one actual positive ordinate from each of 6,000 brackets in field A and 8,000 in field B and pair it with its negative in the aggregate family. Multiplicity and any additional positive zeros remain in the tail. The companion of a nonreal character belongs to its conjugate. The total selected term at step $`j`$ from an ordinate $`t`$ is

```math
w(t)\sin((k+j+1/2)at),\qquad
w(t)=-4G(t)\sin(at/2)/a .
```

The all-ordinate mass is $`M_i=V_i(0)`$. Its sine-square moment is $`(V_i(0)-V_i(a))/2`$. From each subtract the selected pairs' contributions, respectively $`2G(t)`$ and $`2G(t)\sin^2(at/2)`$. Write $`R_0,R_2`$ for the two residuals. Cauchy–Schwarz gives omitted move amplitude at most $`2\sqrt{R_0R_2}/a`$. This subtraction needs one sign-changing zero per bracket. The proof permits additional zeros anywhere, with arbitrary multiplicity. [Endpoint verifier](check_grid_zeros.py) · [Full moment and primary grid calculation](grid_certificate.py)

For the finite head, set $`\theta=t a/(4\pi)`$ in turns. Approximate $`\theta`$ by an integer tick divided by $`2^{64}`$, and $`w(t)`$ by a signed integer divided by $`2^{24}`$. The script computes enclosing balls over the full root bracket. If the quantized values are $`\hat\theta,\hat w`$ and the odd multiplier is $`d=2k+2j+1`$, the error of one term is bounded by

```math
|w-\hat w|+|\hat w|\left(2\pi |d|\,|\theta-\hat\theta|+\frac{\pi}{32768}+2^{-24}\right).
```

The last two terms bound nearest-node phase sampling and the rounding of its sine. The code checks every sine-table rounding. Unsigned 64-bit multiplication computes turns modulo one and therefore performs exact phase reduction even at the large exponents. The signed dot product fits in 64 bits by the checked sum of absolute weights. For negative starts the primary verifier uses $`P_i(-x-a)=-P_i(x)`$, reversing the order and both signs. The transcript checks every required component after the phase errors and the omitted spectral mass are subtracted. [Fixed-point evaluator](grid_sum.cpp) · [Reflection identity](DoublingPatterns.lean)

The independent evaluator uses cosine at the centres of 65,536 cells, with signed weights and ordinates quantized by taking floors. It evaluates negative integer exponents directly. Its tail bound uses the cubic $`J(u)=(11+75u-48u^2+16u^3)/54`$. For $`z=\sqrt u`$ the difference $`J(z^2)-z`$ equals $`(8/27)(z-1/2)^2(z-1)^2((z+3/2)^2+1/2)`$. Thus the unlisted amplitude is bounded by the unlisted sum of $`(2/a)G(t)J(\sin^2(at/2))`$. The full mass for this polynomial is determined by $`V_i(0),V_i(a),V_i(2a),V_i(3a)`$ via the elementary power-reduction identities. The independent source evaluates these four moments from the gamma series and the primitive terms at 2. Through norm 16 the odd-prime trace vanishes and the boundary terms have zero weight. [Independent shifted-grid and cubic calculation](grid_independent.py) · [Original cubic moment proof](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-factor-400/PROOF.md)

## A common return with an explicit density

Write $`Q=11632`$. Let $`N_\chi(t)`$ count zeros at both signs of the ordinate, with multiplicity. Corollary 1.2 of Bennett, Martin, O'Bryant and Rechnitzer bounds the absolute error in $`(t/\pi)\log(qt/(2\pi e))`$ by $`0.247\log(qt)+6.894`$. For each of $`T=600,250000`$ this gives $`N_\chi(t)\le(t/\pi)\log(Qt)`$ for every $`t\ge T`$. At $`T`$ the discarded term $`(T/\pi)\log(2\pi e)`$ dominates the error, and the difference increases thereafter by differentiation. [Counting zeros of Dirichlet L-functions](https://doi.org/10.1090/mcom/3599) · [Directed constant checks](return_bounds.py)

Since $`G(t)\lt 17/t^2`$ for $`t\gt 0`$, partial summation, dropping the nonpositive boundary term at $`T`$, bounds the omitted mass in either host by

```math
R(T)\le 8\cdot34\int_T^\infty\frac{\log(Qt)}{\pi t^2}\,dt
 =\frac{272(\log(QT)+1)}{\pi T}.
```

Across the two hosts there are at most $`d=8T\log(QT)/\pi`$ positive phase coordinates with ordinate at most $`T`$. Indeed there are sixteen factors, and the aggregate is conjugation stable. Central terms have constant phase. Repetitions can be counted once. If every retained positive phase changes by at most $`2\pi/H`$ in angle modulo $`2\pi`$, the change in each $`V_i`$ is bounded uniformly by $`180(2\pi/H)+2R(T)`$. The change in the normalized profile difference $`P_i`$ is consequently at most

```math
E(T,H)=\frac2a\{180(2\pi/H)+2R(T)\}.
```

For $`T=250000`$, $`H=2^{18}`$, the certificate gives $`E(T,H)\lt 0.06`$ and $`18d\lt 252000000`$. Consider the translations $`t`$ for which the distance from $`\gamma t/(2\pi)`$ to the nearest integer is below $`1/H`$ for all retained coordinates. The product of circle tents $`(1-H\|\gamma t/(2\pi)\|)_+`$ is a nonnegative minorant of the indicator of this event. Each tent has constant Fourier coefficient $`1/H`$, nonnegative absolutely summable Fourier coefficients and value one at the origin. This follows by convolving two interval indicators of width $`1/H`$ and scaling to peak one. The product's zero multi-index has coefficient at least $`H^{-d}`$. Coincident zero frequencies add nonnegative terms.

On any long real interval the mean of a nonconstant character tends uniformly in the initial point to zero. On the lattice $`t=k h`$ with any fixed real $`h\gt 0`$ the finite geometric-series formula gives the corresponding uniform limit for every frequency whose sampled character is nonconstant. Sampled constants again contribute positively. Absolute summability permits taking limits in both cases. Thus the returning real shifts have lower uniform density at least $`2^{-252000000}`$, and the returning integer exponents have the same lower uniform density on every fixed lattice. Here lower uniform density means the limit inferior, as interval length tends to infinity, of the infimum over starting positions on the positive half-line of the fraction of that interval occupied by the set.

At any sufficiently large returning shift use the sixteen actual starts $`y=z_w e^t/2`$. The error in each prime-only move can be taken below $`0.01`$, simultaneously over both fields and both steps. Their signed values exceed $`9.1-0.06-0.01\gt 9`$. So all sixteen two-step patterns occur at arbitrarily large real starts, with every normalized move exceeding 9. For the dyadic words take the lattice of translations with spacing equal to the logarithm of two. The same shifts preserve all 5,376 finite-head witnesses simultaneously. They exceed the table thresholds minus 0.07 in the actual prime moves. Consequently every four-, five- and six-step code occurs at arbitrarily large integer power-of-two starts with normalized component margins greater than 14, 9 and 4, respectively. Each code has its own fixed integer offset, so its successful starting exponents themselves have lower uniform natural density at least the same bound. A common real translation similarly gives positive lower uniform logarithmic density among real starting heights. [Earlier one-host common-return argument](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-paired-density/PROOF.md)

There is an explicit height-dependent upper bound on the common return. Let $`K=2143656793`$ be the largest absolute exponent in the published table. For a sufficiently large integer $`N`$ put $`L=N+K+2`$. Partition the retained phase torus into $`H^n`$ cells as above, with $`n`$ the number of distinct retained positive ordinates. Among the $`H^n+1`$ points belonging to the log shifts $`uLa`$ for $`0\le u\le H^n`$, two share a cell. Their positive index difference $`m`$ satisfies $`1\le m\le 2^{252000000}`$. The shift $`mLa`$ therefore preserves all the certified margins. Every prescribed word then starts at a power of two with exponent $`mL+k-1\gt N`$. Its last doubling finishes at an exponent bounded above by $`2^{252000000}(N+K+2)+K+5`$. This inequality bounds exponents. The lower cutoff for $`N`$ comes from the prime-only convergence.

## What longer patterns can occur

Fix a tail on which the errors in the normalized moves are below $`0.01`$. For $`T=600`$ and $`H=2^{10}`$ the directed check gives $`E(T,H)\lt 17.9`$ and $`10d\lt 241000`$. Divide the circle into $`H`$ equal half-open intervals at each retained frequency. If two starting log centres fall in the same phase cell, their angular differences are below $`2\pi/H`$. Adding $`ja`$ to both centres preserves each angular difference. In every position of the two actual words their normalized moves therefore differ by less than $`17.9+0.02\lt 18`$ in each field. Opposite signs both exceeding 9 in magnitude would differ by more than 18. Thus two distinct words of any common length with all moves of magnitude greater than 9 require different cells. There are at most $`2^{241000}`$ such words, independently of their length.

There is also a common integer repetition lag for every sufficiently late word at the smaller margin $`1/25`$. With the finer partition $`T=250000`$, $`H=2^{18}`$, place the phase vectors of $`ka`$ for integers $`0\le k\le H^n`$ in the $`H^n`$ cells, where $`n`$ is the actual number of retained distinct positive ordinates. Two share a cell. Their index difference $`m`$ satisfies $`1\le m\le2^{252000000}`$. It moves every phase by less than $`2\pi/H`$, at every starting height. Thus whenever positions $`j,j+m`$ in a late word both have normalized magnitude greater than $`1/25`$ in a field, their signs agree, because the difference is less than $`0.06+0.02=2/25`$. The same $`m`$ works for both fields, every length and every sufficiently late start.

Finally $`2M_A\lt 57\cdot9a`$ and $`2M_B\lt 56\cdot9a`$ by the directed mass calculation. Adding that many consecutive normalized moves of the same sign and magnitude greater than 9 would give a score difference larger than the corresponding right side. The actual endpoint scores each lie between $`1-M_i-o(1)`$ and $`1+M_i+o(1)`$, a contradiction sufficiently far out. These run lengths are explicit upper obstructions. The exact maximum run length at margin 9 is open. [Budget checks](return_bounds.py) · [Sign, telescoping and counting deductions in Lean](DoublingPatterns.lean)
