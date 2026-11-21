# Atomic width mixtures and smooth synthesis for a split–inert prime race

## Field and observables

Take $`K=\mathbf Q(\zeta_{16})^+`$, $`F=K(\sqrt{397669})`$ and $`E=K(\sqrt{-547},\sqrt{-727})`$. Let $`\epsilon_v`$ be the split sign $`+1`$ or inert sign $`-1`$ in $`E/F`$, supplied by the first-step theta trace on good primes of $`F`$. Remove the places above $`2,547,727,343057911569,343070636977`$. The eight primitive odd Dirichlet factors of $`\zeta_E/\zeta_F`$ have conductors $`547,4376,8752,8752,727,5816,11632,11632`$. Their aggregate multiset of zero ordinates is invariant under negation. Write $`m_t`$ for the total multiplicity at a distinct ordinate $`t`$. [Field and trace lemma](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-short-shifts/RaceTrace.lean#L83) · [Primitive theorem and independent checks](https://math.mit.edu/~drew/grh/)

At width $`b=e^a\gt 1`$, put $`r=\sqrt b`$, $`c_b=2(r-1)^2`$, $`W_b(t)=\max(0,\min(b-b/t,b^2/t-1))`$ and

```math
S_b(y)=-\frac{\sum_v\epsilon_v\log(Nv)W_b(Nv/y)}{r c_b\sqrt y},
\qquad Q_b(x)=S_b(e^x/b).
```

A common geometric centre is essential. Let $`\mu`$ be any finite real signed Borel measure with compact support in $`(0,\infty)`$, in the log-width coordinate $`a`$. Define

```math
P_\mu(x)=\int\frac{c_{e^a}}{2e^{a/2}} Q_{e^a}(x)\,d\mu(a),\qquad
A=\int\cosh(a/2)\,d\mu(a),\qquad
F(t)=A-\int\cos(at)\,d\mu(a).
```

The fixed positive prefactor in the integrand is a normalization. Any finite linear combination of the original aligned scores is included after rescaling the atom masses. The centred score has the limiting profile

```math
P_\mu(x)=P_*(x)+o(1),\qquad
P_*(x)=\int\frac{c_{e^a}}{2e^{a/2}}\,d\mu(a)
   +\sum_{t\ {\rm distinct}}\frac{m_t F(t)}{t^2+1/4}\,e^{itx}.
```

Indeed the known width multiplier is $`[b+1-2r\cos(at)]/[c_b(t^2+1/4)]`$. The prime-only explicit formula has an error tending to zero uniformly on compact sets of widths. Its quadratic spectral majorant is summable against the zero count. Thus both the limit and spectral sum can pass through the finite measure integral. [Prime-only inversion on compact families of widths](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-width-cusp/PROOF.md#the-prime-only-limit-at-a-fixed-width)

## The theorem

Assume the generalised Riemann hypothesis (GRH) for the primitive Dirichlet functions in the 547/727 race and its base field, including the Riemann zeta function. Alpöge's announced Dirichlet GRH theorem supplies this input. [Announcement](https://alpo.ge/remark.pdf) [Independent formal checks](https://math.mit.edu/~drew/grh/) The construction at every finite height below does not use it.

Write $`N(T)=\sum_{0\lt t\le T}m_t`$, $`D_T=N(2T)-N(T)`$, and let $`\mathcal A=A^2+\tfrac12\sum_{a\gt 0}\mu(\{a\})^2`$. Then

```math
\lim_{T\to\infty}\frac1{D_T}\sum_{T\lt t\le2T}m_t|F(t)|^2=\mathcal A .
```

Suppose $`\mathcal A\gt 0`$. This covers every nonzero finite atomic combination, every mixture with at least one atom and every nonzero positive measure. For log displacement $`h\downarrow0`$, put $`\Delta_h^2 f(x)=2f(x)-f(x+h)-f(x-h)`$ and write $`M`$ for the mean in $`x`$, the limit of normalized long-interval integrals. Then

```math
M\!\left(|\Delta_h^2 P_\mu|^2\right)\asymp_\mu h^3\log(1/h).
```

Here $`\asymp_\mu`$ gives positive upper and lower constants for all sufficiently small $`h`$. If $`c_t=m_tF(t)/(t^2+1/4)`$, the spectral mean-square Sobolev norm of order $`s`$ is finite exactly when $`s\lt 3/2`$. By definition its squared nonconstant part is $`\sum_{t\ne0}(1+t^2)^s|c_t|^2`$. At the critical order, the partial sum through $`|t|\le T`$ is of order $`(\log T)^2`$. These statements concern the limiting profile and logarithmic averages of the finite prime score. For balanced atomless measures the first assertion still gives a zero limit. A smooth signed construction appears below.

## Averaging the numerator at the actual zeros

For each fixed $`v\ne0`$, the normalized sum of $`e^{ivt}`$ over the positive ordinates in $`(T,2T]`$, with multiplicity, tends to zero. To see this, take $`z=\exp(|v|)\gt 1`$. The Landau formula for a fixed primitive Dirichlet function bounds $`\sum_{0\lt \Im\rho\le X} z^\rho=O_z(X)`$, including the integer prime-power case. Lemma 2 of Ford, Soundararajan and Zaharescu gives this directly, interpreting the quotient at an integer by continuity. GRH removes the factor $`z^{1/2}`$. Sum over the eight characters, take differences at $`2T,T`$ and conjugate for negative $`v`$. The usual zero count gives $`D_T\sim(4/\pi)T\log T`$. In particular even a frequency equal to the logarithm of a prime power disappears after this normalization. [Landau formula and zero count, Section 3 and Lemma 2](https://arxiv.org/pdf/0805.2745#page=13) · [Primitive Dirichlet theorem](https://github.com/leanprover-community/riemann-hypothesis/blob/371a1a26d5e17f604428e3a8ce6ce341cf517434/src/Main.lean#L12)

Expand $`|F|^2`$ and integrate twice against $`\mu`$. The averages of $`\cos(at)`$ tend to zero for each $`a\gt 0`$. In $`2\cos(at)\cos(bt)=\cos((a-b)t)+\cos((a+b)t)`$, the sum frequency stays positive. The difference frequency survives exactly on the diagonal $`a=b`$. The averages have modulus at most one, so dominated convergence with respect to the product of total variations applies. The diagonal mass is $`\sum_a\mu(\{a\})^2`$, since a finite measure has at most countably many nonzero atoms. This proves the displayed limit, allowing arbitrary rational relations among the widths. It is an application of the diagonal argument in Wiener's lemma. [Wiener's lemma along subsequences](https://arxiv.org/abs/1701.00101)

## Multiplicities and the two-sided energy estimate

Parseval for absolutely convergent almost-periodic series gives

```math
M(|\Delta_h^2P_*|^2)=\sum_{t\ne0}
   \frac{4 m_t^2 |F(t)|^2(1-\cos(ht))^2}{(t^2+1/4)^2}.
```

For a finite trigonometric polynomial integrate its cross terms, then use absolute uniform convergence. The mean is uniform in the starting point. The prime-only $`o(1)`$ error is bounded on finite intervals and tends to zero on the positive half-line, giving the same mean there, uniformly for starting points $`x\ge0`$. Shared and repeated ordinates must be combined before squaring.

For each one of the eight primitive characters, write $`m_{\chi,t}`$ for its multiplicities. The individual pair-correlation theorem of Garunkštis and Paliulionytė gives

```math
\sum_{0\lt t\le T}m_{\chi,t}^2\ll_\chi T\log T,\qquad
\sum_{0\lt t\le T}m_t^2\le8\sum_\chi\sum_{0\lt t\le T}m_{\chi,t}^2\ll T\log T.
```

Here is the deduction, which only uses their proved range $`|\alpha|\le1`$. Integrate their normalized pair sum against $`(1-|\alpha|)_+`$. Its Fourier transform, in angular frequency $`v`$, is $`(\sin(v/2)/(v/2))^2\ge0`$ and equals one at zero. On the critical line their rational weight is $`4/(4+(t-t')^2)\gt 0`$. Every pair at the same ordinate contributes one before normalization and all other integrated pairs are nonnegative. For each fixed conductor $`q`$, Theorem 1.1 applies for all sufficiently large $`T`$ since $`q\le\log T`$. Integrating its even main term $`(qT)^{-2|\alpha|}(\log(qT)+O(1))+|\alpha|`$ gives at most $`1+1/3+o(1)`$ and the uniform error tends to zero. This bounds the diagonal by $`(4/3+o(1))T\log(qT)/(2\pi)`$. The second inequality is finite Cauchy–Schwarz. This use requires neither a pair-correlation conjecture nor simplicity. It does use GRH, since the positivity of the weight and the unit contribution of a coincident pair need the zeros on the critical line. [Garunkštis and Paliulionytė, Theorem 1.1](https://hrcak.srce.hr/en/file/495424)

As $`|F(t)|\le |A|+\|\mu\|_{\rm TV}`$, the last bound prices a positive shell $`T\lt t\le2T`$ by $`O_\mu(T\log T)`$ for $`\sum m_t^2|F(t)|^2`$. When $`\mathcal A\gt 0`$, its matching lower bound follows from $`m_t^2\ge m_t`$ and the proved numerator limit. On $`1/h\lt t\le2/h`$, $`1-\cos(ht)\ge1-\cos1\gt 0`$. This shell alone supplies $`\gg_\mu h^3\log(1/h)`$. For $`0\lt t\le1/h`$, use $`1-\cos(ht)\le (ht)^2/2`$ to bound the sum by $`O_\mu(h^4(1/h)\log(1/h))`$. For larger $`t`$ use $`1-\cos(ht)\le2`$ and sum successive dyadic shells of sizes $`T=2^k/h`$. Their contributions are $`O_\mu(T^{-3}\log T)`$, a convergent geometric sum of the required order. The negative ordinates give an identical contribution.

The shell contribution to the squared Sobolev norm is of order $`T^{2s-3}\log T`$. Sum over $`T=2^k`$. At $`s=3/2`$ the sum through $`T`$ is comparable to $`(\log T)^2`$. Below the threshold it converges and above it diverges. This also proves existence of one derivative in logarithmic mean square.

For the width-two score, let $`C_2(h)`$ be its limiting upper symmetric second difference. It is asymptotic to $`[24/c_2]h\log(1/h)`$. Applying the mean-square estimate to $`Q_2`$ and Chebyshev's inequality shows that for any fixed $`\delta\gt 0`$ the upper uniform log density of $`|\Delta_h^2Q_2(x)|\ge\delta C_2(h)`$ is $`O_\delta(h/\log(1/h))`$. The long-interval mean is taken first, and the bound holds as $`h\downarrow0`$. [Exact small-displacement maximum](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-width-cusp/PROOF.md#the-cusp-at-every-width)

## Constructing a smooth race at every height

Choose any real $`\phi\in C_c^\infty((0,\infty))`$ and put $`d\mu(a)=(\phi''(a)-\phi(a)/4)\,da`$. Integrating by parts twice yields $`A=0`$ and

```math
F(t)=(t^2+1/4)\int_0^\infty\phi(a)\cos(at)\,da,\qquad
c_t=m_t\int_0^\infty\phi(a)\cos(at)\,da .
```

The cosine transform decays faster than every inverse power. The zero count permits termwise differentiation to every order with absolute uniform convergence. There is also an exact identity for the prime sum at each finite geometric centre,

```math
P_\mu(x)=-\frac12\sum_v\frac{\epsilon_v\log(Nv)}{\sqrt{Nv}}\,
                      \phi(|x-\log Nv|).
```

Indeed, writing $`n=Nv`$, the centred width kernel applied to the prime ideal is $`-\epsilon_v(\log n)g_{e^a}(x-\log n)/\sqrt n`$, where $`g_{e^a}(z)=(2e^{a/2}/c_{e^a})\sinh((a-|z|)_+/2)`$. Multiplication by the prefactor in $`P_\mu`$ cancels its normalization. Set $`d_0=|x-\log n|`$. Twice integrating $`(\phi''-\phi/4)\sinh((a-d_0)/2)`$ on $`[d_0,\infty)`$ leaves $`\phi(d_0)/2`$, since the hyperbolic sine vanishes at $`d_0`$ and its right derivative there is $`1/2`$. Only finitely many ideals contribute at any fixed centre. The same identity holds for compactly supported $`C^2`$ choices of $`\phi`$. [Finite-height verification](synthesis.py)

For a nonzero smooth example take $`\phi(a)=\exp(-1/[1-4(a-3/2)^2])`$ for $`1\lt a\lt 2`$ and zero elsewhere. The function and all its derivatives vanish at its endpoints. A certified positive ordinate is below $`1/20`$. At that ordinate the integrand is nonnegative. On $`[5/4,7/4]`$ it is at least $`e^{-4/3}\cos(1/10)`$, so the corresponding cosine transform exceeds $`1/10`$. The constructed profile is genuinely nonconstant. [Zero and coefficient checks](shell_energy.py)

## A balanced two-width example

Use the original-score combination $`c_2 Q_2-(3/5)c_4 Q_4`$. Its numerator is $`F(t)=-2\sqrt2\cos(t\log2)+(12/5)\cos(2t\log2)`$. The constant cancels and its average square at high ordinates tends to $`172/25`$. For $`T=1,2,4,8,16,32`$ respectively, the verified lower bounds on $`T^3M(|\Delta_{1/T}^2P_\mu|^2)/\log(2T)`$ are $`79,51,73,35,31,28`$. These are six finite lower bounds. The asymptotic order comes from the preceding shell argument.

The program independently checks 1,564 Hardy endpoint signs for disjoint brackets in each character, then uses interval arithmetic to lower-bound the selected contributions in each shell. The one additional low bracket verifies the smooth example. Conjugate partners count on the negative side, and $`m_t^2\ge m_t`$ allows coincident brackets from different characters. A second calculation uses half-angle and double-angle formulas. Directed integration tests the finite-height identity on the compact polynomial $`(a-1)^3(2-a)^3`$ on $`[1,2]`$, extended by zero. Exact rational antiderivatives check it over the full interval. [Brackets](evidence/selected-brackets.tsv) · [Independent checks](independent.py) · [Transcript](evidence/check.txt)
