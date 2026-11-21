# Two-sided proportions on every sufficiently long interval

## The prime score and the conclusion

Use $`K=\mathbf Q(\zeta_{16})^+`$, $`F=K(\sqrt{397669})`$ and $`E=K(\sqrt{-547},\sqrt{-727})`$. For a good prime ideal $`v`$ of $`F`$, put $`\epsilon_v=1`$ for split and $`-1`$ for inert in $`E/F`$. The first-step theta trace gives these signs. Omit places above $`2,547,727,343057911569,343070636977`$. Put $`a=\log2`$, $`r_0=\sqrt2`$, $`c=2(r_0-1)^2`$ and $`B=r_0c`$. Define

```math
W(u)=\max(0,\min(2-2/u,4/u-1)),\qquad
S(y)=\frac1{B\sqrt y}\sum_v-\epsilon_v\log(Nv)\,W(Nv/y),\qquad
f(x)=\frac{S(2y)-S(y)}a,\quad x=\log(2r_0y).
```

For $`\sigma\in\{-1,1\}`$ and $`H\ge1`$ let $`p_\sigma(H,T)`$ be the fraction of log starting heights in $`[\log H,\log H+T]`$ satisfying $`\sigma f(x)\gt 1.9735`$. Lebesgue measure in $`x`$ equals $`dy/y`$. Assume the generalised Riemann hypothesis (GRH) for the primitive Dirichlet functions in the 547/727 race and its base field, including the Riemann zeta function. Then the following inequalities hold for **every** $`H\ge1`$.

| For every $`T`$ at least | $`p_\sigma(H,T)`$ exceeds |
|---:|---:|
| 20,000 | $`1/100`$ |
| 50,000 | $`1/16`$ |
| 100,000 | $`1/12`$ |
| 250,000 | $`1/10`$ |

Further, $`\liminf_{T\to\infty}\inf_{H\ge1}p_\sigma(H,T)\gt 3/20`$ for each sign. The two sets of starts are disjoint. The table and the lower density both use the hypothesis through the comparison with prime scores. Alpöge's announced Dirichlet GRH theorem supplies this input. [Announcement](https://alpo.ge/remark.pdf) [Independent formal checks](https://math.mit.edu/~drew/grh/) [Field and trace in Lean](https://github.com/leecrossley/theta-race-short-shifts/blob/483c31f2df2132ed2165620ea55e20e969228bfb/RaceTrace.lean#L83) [Pinned explicit formula](https://github.com/leecrossley/theta-race-filtered-deadline/blob/2400dc5ac2a7969d41e1a4bde1eb44dd6c63b7a2/PROOF.md)

## A finite set of frequencies determines the first four diagonal moments

The eight odd primitive Dirichlet factors of $`\zeta_E/\zeta_F`$ have conductors $`547,4376,8752,8752,727,5816,11632,11632`$. Under GRH their nontrivial zero ordinates are real. Conjugation pairs the full multiset. Put

```math
G(t)=\frac{3-2r_0\cos(at)}{c(t^2+1/4)},\qquad
A(t)=-4G(t)\frac{\sin(at/2)}a.
```

The limiting profile is the uniformly absolutely convergent sum of $`A(t)\sin(tx)`$ over positive ordinates with multiplicity. The 6,000 signed brackets select frequencies with omitted absolute mass below 1.52. The prime-only error is below 0.007 for $`x\ge100`$. Its base-field term also uses GRH for the base field, including zeta. These bounds allow repeated or unlisted zeros on the critical line. Retain the 800 rows with the largest $`|A|`$, listed in `evidence/head-rows.txt`, in $`P(x)`$ and call the sum of the remaining 5,200 listed sines $`R(x)`$. Thus $`|f-P-R|\lt 1.527`$ for $`x\ge100`$. [Omitted mass and 12,000 Hardy signs](https://github.com/leecrossley/theta-race-factor-400/blob/f1773b07a55e76947d4e0c6ab19a112fd91fe6cf/PROOF.md) [Prime-only bound](https://github.com/leecrossley/theta-race-filtered-deadline/blob/2400dc5ac2a7969d41e1a4bde1eb44dd6c63b7a2/PROOF.md) [Head rows](evidence/head-rows.txt) [Primitive theorem checks](https://math.mit.edu/~drew/grh/)

Write $`X=P/27`$ and $`b_j=|A(t_j)|/54`$ for the 800 rows. Consider the union of the frequencies of $`1,X,X^2`$. Its intended members are zero, $`\pm t_j`$, $`\pm2t_j`$, and both signs of $`t_i+t_j`$ and $`t_i-t_j`$ for $`i\lt j`$. All 1,281,601 listed nodes are pairwise distinct at the true ordinates. To check this, a bracket $`[l,u]/10^{14}`$ is represented by centre $`l+u`$ and radius $`u-l`$ in units $`1/(2\cdot10^{14})`$. Every difference of two node frequencies has uncertainty at most four times the largest original radius. Sort the integer centres and compare their smallest gap with this error. For the 5,200-row residual the analogous signed single-frequency test uses twice its largest radius. [Direct construction](moments.py) · [Signed-pair construction](independent.py)

The zero Fourier coefficients of $`1,X,X^2,X^3,X^4`$ are $`1,0,m_2,0,m_4`$. The separation test excludes every coincidence except the opposite-pair cancellations already counted. Orthogonality gives

```math
m_2=2\sum_j b_j^2,\qquad m_4=3m_2^2-6\sum_j b_j^4,\qquad
m_R=\frac12\sum_{t_j\in R} A(t_j)^2.
```

In the last sum the index runs through the residual rows. One may also obtain $`m_4`$ by squaring the Fourier coefficients of $`X^2`$ to get $`6\sum b_j^4+24\sum_{i\lt j}b_i^2b_j^2`$. Directed calculations on whole brackets give $`1.0065\lt m_2\lt 1.0066`$, $`m_4\lt 3.0034`$ and $`m_R\lt 0.00923`$. Bounds on higher integer relations between zeros are unnecessary. [Two coefficient calculations](evidence/check.txt)

## Uniform finite-window moment errors

Fix a row of the table, with minimum length $`T_0`$. For $`T\ge T_0`$ use a probability measure on the chosen log interval obtained by adding independent uniform variables of lengths $`(1-\beta)T`$ and $`\beta T`$ and translating to the interval. Use $`\beta=1/10,1/20,1/50,1/50`$ in the four rows. The density is at most $`1/((1-\beta)T)`$. At a nonzero frequency difference of magnitude $`d`$ its Fourier transform has modulus at most

```math
K(d)=\min\left(1,\frac{2}{(1-\beta)T_0d}\right)
      \min\left(1,\frac{2}{\beta T_0d}\right).
```

For $`k=1,2,3,4`$ use the pairs of powers $`(i,j)=(0,1),(1,1),(1,2),(2,2)`$, respectively. Then $`X^k=X^i\overline{X^j}`$ and the error between its measured average and its diagonal coefficient is at most

```math
e_{ij}=\sum_{\lambda\ne\mu}w_i(\lambda)w_j(\mu)K(|\lambda-\mu|).
```

Here $`w_i`$ upper-bounds the absolute Fourier coefficients of the indicated power. In $`X^2`$ the coefficient magnitudes are bounded at zero by $`2\sum b_j^2`$, at $`\pm2t_j`$ by $`b_j^2`$, and at distinct-index nodes by $`2b_ib_j`$. Using the signed linear coefficients of $`R`$ in both columns gives an error $`e_R`$ for its second moment. Signs of the columns can be discarded in these error bounds. The diagonal values above use the original coefficients.

The following integer computation encloses these sums. All charges are rounded up in units $`10^{-7}`$. Put $`D=2\cdot10^{14}`$ and let $`u`$ bound the uncertainty of the difference of two nodes in integer ticks. Take the consecutive powers of two $`h`$ from the largest power at most the minimum positive centre gap through the largest power at most the full span. The verifiers check $`h\gt u`$ even at the first power. Put $`k(h)=K((h-u)/D)`$. A pair whose centre distance lies in $`[h,2h)`$ has true weight factor at most $`k(h)`$. If $`C(h)`$ bounds the total ordered charge of distinct pairs with centre distance strictly less than $`2h`$, summation by parts gives the error bound

```math
10^{-14}\left(\sum_h [k(h)-k(2h)]C(h)+k(2h_{\max})C(h_{\max})\right).
```

Indeed a pair first counted at power $`h`$ receives total factor $`k(h)`$ by telescoping. Every summand coefficient is nonnegative, which also permits upper bounds for the cumulative masses. The primary verifier groups centres in bins of width $`h/2`$. Any eligible pair is in bins whose indices differ by at most four. It subtracts self-pairs. The independent verifier takes binary searches at $`p+2h`$ in the sorted centre list, using prefix sums to count the exact cumulative charge for each column pair in both orders. It recomputes the coefficients from the half-angle form of $`G`$ and enumerates all positive sums of signed frequencies. The rational kernel arithmetic and the final inequalities are shared. Integer dot products are checked against total column-mass bounds before use. [Bin verifier](moments.py) · [Direct-distance verifier](independent.py) · [Rational inequalities](certificate.py)

## From moments to signed prime moves

Let $`t=5/27`$, $`b=8/5`$, $`v=39/20`$, $`d=10838/3823`$ and take

```math
q(z)=\frac{16515360000}{325057991399}(z-t)(z+b)^2(d-z)=\sum_{k=0}^4 c_k z^k.
```

Then $`q(z)\le\mathbf1_{\{z\gt t\}}`$. On $`z\le t`$ its product is nonpositive. On the whole line $`1-q(z)`$ factors as $`(z-v)^2`$ times $`(16515360000z^2+67380496000z+91322607600)/325057991399`$. The quadratic is positive by completing the square. In particular the strict threshold causes no endpoint issue. [Quartic bounds in Lean](ExcursionProportions.lean)

For either sign the diagonal value of $`q(\sigma X)`$ is $`D_q=c_0+c_2m_2+c_4m_4`$. The averaged polynomial loses at most $`\sum_{k=1}^4|c_k|e_{ij}`$, with each power matched to its pair. The probability that $`|R|\gt 1.49`$ is bounded by $`(m_R+e_R)/1.49^2`$. On the remaining starts with $`\sigma P\gt 5`$ and $`x\ge100`$, the signed prime score exceeds $`5-1.49-1.52-0.007\gt 1.9735`$. Converting the probability lower bound to Lebesgue fraction costs a factor $`1-\beta`$. At most $`100/T`$ of any allowed interval has $`x\lt 100`$. Consequently the guaranteed fraction is at least

```math
(1-\beta)\left(D_q-\sum_{k=1}^4|c_k|e_{ij}
                     -\frac{m_R+e_R}{1.49^2}\right)-\frac{100}{T_0}.
```

Both verifiers clear every tabulated fraction. For the 250,000 row their rational lower bounds exceed 0.101178 and 0.101716, respectively. [Reproducible transcript](evidence/check.txt)

## The lower uniform density

For the head polynomial use Haar probability on the closure of $`x\mapsto(e^{it_jx})_j`$ in the finite torus. The time averages of every continuous function on that group converge to its Haar mean uniformly in the initial time. Indeed each character is either constant on the group or has interval average of modulus at most $`2/(T|\omega|)`$ for a fixed $`\omega\ne0`$. Finite trigonometric approximation gives the statement for continuous tests. For an open event its lower interval proportion is at least its Haar measure by continuous approximation from below. Inversion preserves Haar measure and sends the sine polynomial $`X`$ to $`-X`$, so the two strict tails outside $`[-t,t]`$ have equal Haar measure.

With $`p=\Pr(|X|\gt t)`$, splitting its second moment at $`t`$ and applying Cauchy–Schwarz gives $`m_2\le t^2+\sqrt{pm_4}`$. For each sign the lower uniform proportion of head excursions is therefore at least $`(m_2-t^2)^2/(2m_4)`$. The interval mean of $`R^2`$ tends to $`m_R`$ uniformly by its finite expansion. Chebyshev's inequality, the same prime-score transfer and the vanishing fraction below $`x=100`$ give

```math
\liminf_{T\to\infty}\inf_{H\ge1}p_\sigma(H,T)
\ \ge\ \frac{(1.0065-25/729)^2}{2\cdot3.0034}
             -\frac{0.00923}{1.49^2}
\ \gt \ \frac3{20}.
```

This argument uses only the four certified diagonal moments, even when other phase relations occur. The exact limiting percentages remain open. [Moment implication and rational margin](ExcursionProportions.lean) · [Verification output](evidence/check.txt)
