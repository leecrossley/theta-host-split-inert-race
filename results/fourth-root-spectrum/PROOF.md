# A finite spectral formula for the unsmoothed dyadic move

## The fourth-root formula

Use either conductor pair $`(547,727)`$ or $`(571,599)`$. The host $`F`$ has degree eight, with the quadratic extension and signs $`\epsilon_v`$ of the maintained theta trace. Keep the five omitted rational primes on each field as in the joint-level paper. Set $`a=\log2`$, $`r=\sqrt2`$ and $`B=2r(r-1)^2`$. For the usual tent $`W`$, supported on $`[1,4]`$, define

```math
S(y)=-\frac1{B\sqrt y}\sum_v\epsilon_v\log(Nv)W(Nv/y),
\qquad f(x)=\frac{S(2y)-S(y)}a,\qquad x=\log(2r y).
```

The eight primitive odd factors are the two quadratic characters of the chosen conductors, each twisted by the four powers of the even quartic character of conductor 16. For any strip zero write $`z=\rho-1/2`$ and $`t=\Im\rho`$. Include both signs, multiplicities and shared zeros with their separate factor multiplicities. The move transform is entire,

```math
J(z)=\frac{2\sinh(az/2)}{a}\,
      \frac{3-2r\cosh(az)}{2(r-1)^2(1/4-z^2)}.
```

The singularities at $`\pm1/2`$ are removable. The trace identification, explicit formula and finite omissions are recorded for both hosts. [Both fields and their omitted primes](https://github.com/leecrossley/theta-host-split-inert-race/blob/joint-levels-and-moves-v1.2.0/results/joint-levels-and-moves/PROOF.md) · [547/727 trace in Lean](https://github.com/leecrossley/theta-race-short-shifts/blob/483c31f2df2132ed2165620ea55e20e969228bfb/RaceTrace.lean#L83) · [571/599 admissible datum](https://github.com/leecrossley/theta-host-split-inert-race/blob/8fca7cad37185e3d7951f0dffff5e2442d1dfa24/HostAdmissible.lean)

**Theorem.** For every real $`x\ge12494`$ put

```math
n=\lceil x\rceil+1,\quad
\sigma=10^{-6}e^{-x/4},\quad R=\sqrt{3n}\,\sigma,\quad
T=8\cdot10^6\sqrt n\,e^{x/4},\quad
\Phi(z)=\left(\frac{\sinh(zR/n)}{zR/n}\right)^n,\quad m=\Phi(1/2).
```

Use value one for the quotient at zero. On either field the following absolute error bound holds without a zero-location hypothesis.

```math
\left|f(x)-\frac1m\sum_{\substack{\rho\ {\rm odd}\\|\Im\rho|\le T}}
             J(z)e^{zx}\Phi(z)\right|\lt0.0008 .
```

The sum is finite and real by conjugation of the family. It asks for both real and imaginary parts of the zeros. It imposes no independence, simplicity or critical-line assumption. This is an error theorem for a finite sum. The release does not enumerate zeros to this cutoff. [Implementation of the error estimates](spectrum.py)

## What the symmetric average changes

Let $`U`$ be a sum of $`n`$ independent uniform real variables on $`[-R/n,R/n]`$. Its law is even, its variance is $`R^2/(3n)=\sigma^2`$, its absolute first moment is at most $`\sigma`$, and its bilateral transform is $`\Phi`$. In particular $`m=\mathbb E\cosh(U/2)\ge1`$. The support has $`R\lt e^{-3000}\lt1/2`$ throughout the stated range. For this last bound use $`n\le x+2`$ and the decrease of $`\sqrt{x+2}\,e^{-x/4}`$.

On a smooth piece, a term at norm $`N`$ in $`S(y e^u)`$ is a linear combination of $`e^{-u/2}`$ and $`e^{u/2}`$. Indeed $`W(v)=2(1-1/v)`$ on $`1\lt v\le2`$ and $`W(v)=4/v-1`$ on $`2\lt v\lt4`$. Outside the tent both coefficients are zero. At fixed $`0\le s\le R`$ the symmetric two-point average of either exponential is $`\cosh(s/2)`$. Subtracting $`\cosh(s/2)S(y)`$ thus cancels every term whose height interval misses the three breakpoints. If $`\lambda=e^s`$ is rational, clearing the square roots leaves the exact rational expression

```math
W(N/(y\lambda))+\lambda W(N\lambda/y)-(\lambda+1)W(N/y).
```

The release compares this sparse expression with the full prime-ideal sums on eight starts on each host and finds nonzero defects of both signs. [Exact rational and prime checks](prime_check.py) · [Branch cancellation](FourthRoot.lean)

For completeness, bound the symmetric defect of one prime ideal with norm $`N`$ at a breakpoint by $`12s\log N/\sqrt y`$. On $`|u|\le1/2`$, the derivative on each branch is bounded by $`e^{1/4}(1/2+2)\log N/(B\sqrt y)`$, using $`0\le W\le1`$ and $`|vW'(v)|\le2`$. The continuous piecewise function has the same Lipschitz bound. The additional term from $`\cosh(s/2)-1\le s`$ costs at most $`s\log N/(B\sqrt y)`$. Their sum is below the stated factor 12. Only norms in $`[ky e^{-s},ky e^s]`$ can contribute, for $`k=1,2,4`$. Each such interval has length at most $`3kys`$. A fixed norm occurs for at most eight prime ideals, since each lies above the same rational prime and the field has degree eight. All relevant norms for both levels have logarithm at most $`x+3`$.

If $`s\le e^{-x/3}`$, count all integers in each interval by its length plus one. At the two levels $`y`$ and $`2y`$, with the normalization $`1/a`$, this gives a bound at most $`10^5(x+3)e^{-x/6}`$. For example the three intervals at one level give at most $`21ys+3`$ integers. Both $`\sqrt y`$ and $`\sqrt{2y}`$ are at most $`e^{x/2}`$ and their reciprocals are below $`2e^{-x/2}`$. These bounds also cover a norm at a breakpoint.

If $`s\gt e^{-x/3}`$, use Brun–Titchmarsh with modulus one on each norm interval for the residue-degree-one primes. Enlarge its endpoints by one to include both ends. If $`L`$ is the original length, then $`L\ge2ys\gt e^{2x/3}/2`$. Thus $`\log(L+2)\gt x/2`$ and each prime in the original interval has logarithm at most $`x+3`$. The sum of its prime logarithms is at most $`6(L+2)`$. Multiply by the at most eight ideals per prime. At higher residue degrees the underlying rational prime is at most $`\sqrt{4y e^s}`$. The sum of residue degrees over any rational prime is at most eight, so the total ideal log weight from all these primes is bounded by $`8\sqrt{4y e^s}(x+3)`$ at this level. The same argument applies with $`2y`$ in place of $`y`$. [Montgomery and Vaughan, The large sieve](https://doi.org/10.1112/S0025579300004708)

Combining both ranges, the symmetric defect of the normalized move at displacement $`s`$ is at most

```math
10^5\bigl[e^{x/2}s^2+(x+3)s+(x+3)e^{-x/6}\bigr].
```

For the factor $`10^5`$, the degree-one contribution at each level uses at most $`48(21ys+6)`$ in log weight before multiplying by $`12s/\sqrt y`$. The higher-degree bound just given includes every prime power norm. Summing two levels and dividing by $`a`$ fits the displayed constant. Average the inequality over $`s=|U|`$ and divide by $`m\ge1`$. Because $`\mathbb E[U^2]=10^{-12}e^{-x/2}`$, the first term becomes constant. The other two terms decrease beyond the lower endpoint. Consequently $`|f(x)-\mathbb E f(x+U)/m|\lt1.1\cdot10^{-7}`$. Using the variance here determines the square-root logarithm in the spectral cutoff. [Scalar check](spectrum.py)

## The high zeros and arithmetic remainder

The exact odd zero terms of $`f(x+U)`$ give $`J(z)e^{zx}\Phi(z)`$ after averaging. The strip bound $`|J(z)|\le1000/(1+t^2)`$ and the total reciprocal-square mass below 1,000 justify termwise integration on this compact interval. Both constants apply on either host. The completed-function real logarithmic derivative at two gives the mass estimate one factor at a time. All conductors are bounded by $`16\cdot547\cdot727`$ on either field. [Transform and mass proof](https://github.com/leecrossley/theta-host-split-inert-race/blob/theta-race-waiting-bound-v1.3.0/results/theta-race-waiting-bound/PROOF.md#total-mass-of-the-omitted-zeros)

For $`|t|\gt T`$, $`|\sinh(zR/n)|\le\cosh(R/(2n))\lt2`$, whence $`|\Phi(z)|\le(2n/(|t|R))^n\le4^{-n}`$ by the definitions of $`T`$ and $`R`$. The unrestricted high odd spectrum costs less than $`10^6 e^{x/2}4^{-n}\lt10^{-10}`$. This estimate uses the full strip $`|\Re z|\le1/2`$. [Cutoff algebra](FourthRoot.lean)

For the base zeros at $`x'=x+U\ge12493`$ use the primitive factors of the even characters on the selected host. Their conductors divide their own base modulus $`Q_b=16pq`$. Apply McCurley to that modulus, and then use the common upper bound $`Q=6362704`$ in all numerical estimates. With $`S(t)=80[\log(2Qt)+1]/(t-2)`$ for $`t\ge8590`$, $`b=e^{10}`$, $`d=e^{25}`$ and $`x_0=12493`$, the total base-zero allowance is

```math
2\cdot10^6e^{-x_0/[20\log(Qb)]}
+2000 S(b)e^{-x_0/[20\log(Qd)]}
+2000 S(d)+\frac{20}{e x_0}.
```

The final term covers the possible simple real exceptional zero. It follows from $`J(0)=0`$, $`|J'(-u)|\lt40`$ for $`0\le u\le1/2`$, and $`u e^{-ux'}\le1/(e x')`$. There is at most one exceptional zero in each base product. The pole at one also vanishes. The estimates for the other base zeros are the same real-part and reciprocal-square bounds proved in the first-failure certificate. The number of factors, field degree and upper bounds on conductors that enter that argument hold on both fields. [Explicit finite-zero comparison](https://github.com/leecrossley/theta-host-split-inert-race/blob/theta-race-factor-400-v1.2.0/results/theta-race-factor-400/PROOF.md#comparing-to-the-actual-primes-throughout-the-search) · [McCurley](https://www.sciencedirect.com/science/article/pii/0022314X84900891)

The higher good-prime powers cost at most $`10^6 e^{-x_0/10}`$, by the degree-eight bound for their absolutely convergent series at $`2/5`$. Trivial zeros cost at most $`128 e^{-(x_0-2)/2}/(1-e^{-(x_0-2)})`$. The five omitted rational primes in each field are below $`10^{12}`$. At most $`1+4/\log2`$ of the powers of each can meet a log support of width four, and each eight-factor logarithmic coefficient costs at most $`8\log p`$. The omitted-prime terms fit in $`10^6 e^{-x_0/3}`$. The base and all three remaining allowances together are less than $`0.00079`$. Since $`m\ge1`$, averaging and division by $`m`$ preserve their bound. The prime-side and high odd costs now prove the theorem. [Higher prime powers](https://github.com/leecrossley/theta-host-split-inert-race/blob/theta-race-waiting-bound-v1.3.0/results/theta-race-waiting-bound/PROOF.md#the-transform-and-the-prime-powers) · [Numerical checks on both fields](spectrum.py)

## The first missed move at factor 400

For the first-failure theorem on 547/727 put $`X=1620548`$, $`\epsilon=4\cdot10^{-8}`$ and $`p=49/50`$. The localization

```math
1620542.262999\lt\log(2\sqrt2 H_*)\le1620542.276
```

and the identification of the loss as successful and the gain as the first miss hold if there is no odd zero in the smaller region

```math
8590\lt |t|\le e^{405161},\qquad
\Re\rho\gt\tfrac12+\epsilon+\frac{p}{X}\log^+(|t|/10^{12}).
```

Here $`\log^+u=\max(0,\log u)`$. Equivalently, failure of the localization or its boundary conclusions forces such a zero, unconditionally. For the separate interval starting at log centre 3895686.221, replace $`X`$ by 3895692, $`\epsilon`$ by $`2\cdot10^{-8}`$ and the upper endpoint by $`e^{973947}`$. All gains there up to the factor $`2e^{5.299}`$ are then below $`1.962\log2`$. The previous respective endpoints were $`e^{810330}`$ and $`e^{1947900}`$. [Original finite regions and all-height certificate](https://github.com/leecrossley/theta-host-split-inert-race/blob/theta-race-factor-400-v1.2.0/results/theta-race-factor-400/PROOF.md#the-first-missed-move-at-factor-400)

For the transfer, write $`U_0(s)=\sum J(it)e^{its}`$ for the projection of all odd zeros with multiplicity. The earlier proof bounds the difference from this projected profile throughout the whole negative intervals and at every certified positive witness. Keep the same retained zeros, projected moment bounds and witness data. For shifts $`|u|\le R`$ and $`|t|\le e^{100}`$, the phase shift costs at most $`10^6 e^{100}R`$ in the aggregate. The coefficient sum on the complementary projected tail costs at most $`100S_{11632}(e^{100})`$ in the difference. Also $`|1-1/m|\le R`$ and $`|U_0|\le10^6`$. These bounds make the total discrepancy between $`U_0(x)`$ and its normalized symmetric average less than $`10^{-10}`$. The even projected moments at the four small shifts use only the old low-ordinate real-part control to $`10^{12}`$ and an unrestricted high tail. They are still valid for both new regions. [Earlier projected moments](https://github.com/leecrossley/theta-host-split-inert-race/blob/theta-race-factor-400-v1.2.0/results/theta-race-factor-400/PROOF.md#what-a-displacement-of-an-odd-zero-can-change)

For $`|t|\le10^{12}`$ the earlier reflected-pair estimate at the translated coordinate $`x+u`$ applies. Platt gives exact equality of the projected and actual zeros through 8,590. For $`10^{12}\lt|t|\le T_{\max}`$ use the same weighted zero count with exponent $`p'=p(x_{\max}+1)/X\lt1`$, since all shifts are smaller than one. Use $`T_{\max}=e^{405161}`$ or $`e^{973947}`$ for the respective window. Each bounds $`T(x)`$ throughout the associated witness set. Above $`T_{\max}`$ the actual smoothed spectrum costs less than $`10^{-10}`$ by the box estimate, and the projected spectrum costs at most $`50 S_{11632}(T_{\max})`$. Divide by $`m\ge1`$ in the averaged comparison. The prime-side error, the shift of the projected reference and both upper tails together fit inside the former odd budgets at all three numerical cutoffs. [Platt's finite theorem](https://arxiv.org/pdf/1305.3087#page=14) · [Interval check](spectrum.py)

Specifically these budgets are $`0.022`$ up to log centre 1620550, $`0.012`$ at the narrow positive bridge below 1177000, and $`0.03`$ for the other negative interval below 3895692. The base and other arithmetic terms remain under $`0.001`$ on each support. Every previous signed margin therefore survives, including the strict real-height cover and the equality of the boundary doubling. This deduction calls on the earlier witness certificate. The reconstruction theorem for the two fields is independent of that localization. [First-failure replay](https://github.com/leecrossley/theta-host-split-inert-race/releases/tag/theta-race-factor-400-v1.2.0)

## A pure fourth-root cutoff with polynomial weights

For either host take $`x\ge12494`$, put $`\lambda=(5\cdot10^6/6)e^{x/4}`$ and write $`z=\rho-1/2=d+it`$. With $`v_+=\max(v,0)`$, define the quintic spline

```math
B(w)=\frac1{2112}\sum_{j=0}^6(-1)^j\binom6j(w+6-2j)_+^5,
\qquad
P_\lambda(d,t)=B(t/\lambda)-i\frac d\lambda B'(t/\lambda)
               -\frac{d^2}{2\lambda^2}B''(t/\lambda).
```

Put $`m_P=1+5/(88\lambda^2)`$. The following estimate is unconditional for the same unsmoothed prime-only move and true zeros as above.

```math
\left|f(x)-\frac1{m_P}\sum_{\substack{\rho\ {\rm odd}\\|\Im\rho|\le 5\cdot10^6 e^{x/4}}}
       J(\rho-\tfrac12)e^{(\rho-1/2)x}P_\lambda(d,t)\right|\lt0.0008 .
```

The cutoff is a constant times the fourth root of the actual prime height. All three polynomial weights vanish at and beyond the endpoint. The real parts of the zeros enter through the two derivative terms and the original explicit formula. The conjugate terms combine to a real sum. [Polynomial and integral checks](quintic.py)

### The kernel on the strip

Use the even probability density $`K_\lambda(u)=20\lambda\,\operatorname{sinc}^6(\lambda u)/(11\pi)`$, with $`\operatorname{sinc}(v)=\sin(v)/v`$ and value one at zero. Its Fourier transform with exponent $`+itu`$ is $`B(t/\lambda)`$. To see this, start with the indicator density of $`[-1,1]`$, whose characteristic function is sinc. Convolution of six copies gives the central box spline. Inclusion–exclusion for the volume of a sliced cube gives its density proportional to $`\sum_j(-1)^j\binom6j(w+6-2j)_+^5`$, with density at zero $`2112/(2^6 5!)=11/40`$. Fourier inversion therefore gives $`\int_{\mathbb R}\operatorname{sinc}^6(v)\,dv=11\pi/20`$ and the stated normalization. Differentiation under the integral twice is justified by inverse-sixth-power decay. At zero $`B'(0)=0`$, $`B''(0)=-5/11`$, so $`\int u^2K_\lambda(u)\,du=5/(11\lambda^2)`$. The full Fourier integral of $`1+du+d^2u^2/2`$ against this density is exactly $`P_\lambda(d,t)`$. [Schoenberg, Cardinal Spline Interpolation](https://doi.org/10.1137/1.9781611970555) · [Independent real Fourier integrals](quintic.py)

In the actual prime sum average only over $`|u|\le1/2`$. Denote its unnormalized local transform by $`M(d,t)=\int_{|u|\le1/2}K_\lambda(u)e^{(d+it)u}\,du`$, and let $`m_C=M(1/2,0)`$. For every $`|d|\le1/2`$, uniformly in real $`t`$,

```math
|M(d,t)-P_\lambda(d,t)|\le14\lambda^{-3}.
```

Indeed the local Taylor remainder of the real exponential is at most $`e^{1/4}|u|^3/48`$. The absolute third moment is below $`\lambda^{-3}`$ by integrating $`v^3\min(1,v^{-6})`$ over the positive half-line and reflecting. Outside the cutoff, the polynomial costs at most

```math
\frac{40}{11\pi\lambda^5}\int_{1/2}^{\infty}
     (1+u/2+u^2/8)\,u^{-6}\,du
=\frac{40}{11\pi\lambda^5}\left(\frac{32}5+2+\frac13\right)
\lt12\lambda^{-5}.
```

These two bounds prove the uniform estimate, including $`|t|\gt 6\lambda`$ where the polynomial vanishes. Taking $`(d,t)=(1/2,0)`$ gives $`|m_C-m_P|\le14\lambda^{-3}`$. Applying the Taylor expansion to the bounded real part, while keeping the entire oscillatory exponential intact, is what allows a fixed spectral width on the scale $`1/\lambda`$. [Scalar enclosures](quintic.py)

### Recovery at every stated height

The previous breakpoint estimate applies on the full local support. Integrate it against $`K_\lambda`$ and use its mass at most one, full variance $`5/(11\lambda^2)`$ and Cauchy–Schwarz for the first absolute moment. It gives

```math
\left|\int_{|u|\le1/2}K_\lambda(u) f(x+u)\,du-m_C f(x)\right|
\le 10^5\left[\frac{5e^{x/2}}{11\lambda^2}
+\frac{(x+3)\sqrt{5/11}}{\lambda}+(x+3)e^{-x/6}\right]
\lt7\cdot10^{-8}.
```

Averaging the arithmetic remainder over this local measure costs below 0.00079 as $`x+u\ge12493`$. The absolute odd coefficient mass before exponential growth is at most $`10^6`$. The uniform strip estimate changes its averaged sum to the finite polynomial sum at cost at most $`14\cdot10^6 e^{x/2}\lambda^{-3}`$. Also $`|f(x)|\le(10^6+1)e^{x/2}`$, from the same explicit formula and arithmetic bound. Replacing $`m_C`$ by $`m_P`$ costs at most $`14(10^6+1)e^{x/2}\lambda^{-3}`$. The last two errors together are below $`10^{-10}`$ on the whole range. Divide by $`m_P\ge1`$. [Prime-side breakpoint estimate](#what-the-symmetric-average-changes) · [Certified constants](quintic.py)

For the first missed $`400H`$ move in the earlier result, its finite ordinate cutoff now becomes $`e^{405153}`$. For the later separate window it becomes $`e^{973939}`$. Keep the same sloping real-part boundaries. To verify the comparison with the projected odd profile, retain the local positive average divided by $`m_P`$. The low ordinates in the prescribed region obey the previous comparison at $`x+u`$, with a change at most one in its ceiling. Above the new cutoff the actual local transform is at most $`14\lambda^{-3}`$, uniformly in the real part, and the projected tail costs at most the previous reciprocal-square bound. For the difference between a projected profile and its normalized local average, split at $`e^{100}`$. On the lower part the difference costs at most $`10^6[e^{100}\sqrt{5/11}/\lambda+5/(11\lambda^2)+12/\lambda^5]`$. Here the last terms bound the departure of $`m_P`$ from one and the deleted probability mass. The rest costs at most $`100S_{11632}(e^{100})`$. This is below $`10^{-10}`$. The three original signed margins survive with the prime-side and normalization errors included. [Finite inequalities](quintic.py) · [First-failure proof](https://github.com/leecrossley/theta-host-split-inert-race/releases/tag/theta-race-factor-400-v1.2.0)

## A sharp truncation at sufficiently large heights

The original odd sum can also be truncated directly, with coefficient one on each retained term $`J(z)e^{zx}`$. For either fixed host put $`L=(x/\log x)^{1/3}`$ and $`T=\exp(x/4-L/4000)`$. For sufficiently large $`x`$ the sum of the absolute values of the discarded odd terms is $`O(e^{-L/1000})`$. Consequently the unsmoothed prime move equals the directly truncated odd sum up to $`o(1)`$. This improves the fourth-root scale by more than any fixed power of $`\log y`$ eventually. The explicit 0.0008 theorem at every $`x\ge12494`$ uses the quintic weights. The direct truncation has no numerical onset here.

For each fixed conductor the log-free density estimate gives $`N(1-u,V)\ll (2qV)^{5u/2}`$ for $`0\le u\le1/2`$, counting zeros with multiplicity. This follows by taking $`\epsilon=1/10`$ in the Huxley–Jutila estimate with exponent $`12/5+\epsilon`$. Summing over eight fixed factors only changes the implied constant. Khale's Theorem 1.1 gives $`1-\Re\rho\ge[10.5\log q+61.5(\log|t|)^{2/3}(\log\log|t|)^{1/3}]^{-1}`$ for $`|t|\ge10`$. The constants below need only hold eventually for the two fixed sets of conductors. [Density estimate, Harper, Zeros Result 1(iv)](https://arxiv.org/pdf/1208.5992#page=8) · [Jutila](https://doi.org/10.7146/math.scand.a-11701) · [Khale, Theorem 1.1](https://arxiv.org/pdf/2210.06457#page=1)

Split the tail into dyadic blocks $`V\lt |t|\le2V`$ from $`V=T`$. Zeros with $`\beta\le1/2`$ contribute $`O(\log(T)/T)`$ in total by the zero count and $`|J(z)|\ll(1+t^2)^{-1}`$. In a block with $`V\le e^{3x/10}`$ the zero-free region gives $`1-\beta\ge L/(100x)`$ for large $`x`$. Also $`\log(4qV)\le x/3`$. If $`F(u)`$ counts the zeros of the block with $`1-\beta\le u`$, partial summation bounds the unscaled exponential sum by

```math
e^{-x/2}F(1/2)+x\int_{L/(100x)}^{1/2} e^{-ux}F(u)\,du
\ll e^{-x/12}+6e^{-L/600}.
```

Here $`F(u)\ll(4qV)^{5u/2}`$ and $`x-(5/2)\log(4qV)\ge x/6`$. Multiplying by $`e^{x/2}/V^2`$ and summing the geometric series from $`T`$ costs $`O(e^{L/2000})`$. The net saving exceeds $`e^{-L/1000}`$ for sufficiently large $`x`$.

For the remaining blocks put $`V_0=e^{3x/10}`$. In the integral and boundary term, compare $`V^{-2+(5/2)u}`$ with $`V_0^{-2+(5/2)u}(V/V_0)^{-3/4}`$. This is valid throughout $`0\le u\le1/2`$. At $`V_0`$ the integral starting at zero is bounded by a constant, and $`e^{x/2}/V_0^2=e^{-x/10}`$. The geometric series is finite even though the ordinates are unrestricted. Finally the base contribution tends to zero by dominated convergence in its absolutely summable reciprocal-square majorant. Each fixed base zero has real part strictly below one, and the pole term vanishes because $`J(0)=0`$. The other arithmetic terms decay explicitly. This proves the assertion for the prime-only move. [Transform and base terms](#the-high-zeros-and-arithmetic-remainder)

## When a fixed power below one quarter works

Fix one of the two hosts. For $`\theta\gt 0`$ measure the actual absolute error budget of a direct truncation by

```math
R_\theta(x)=\sum_{|t|\gt e^{\theta x}} |J(\rho-1/2)|e^{(\beta-1/2)x}.
```

For $`1/2\lt \sigma\lt 1`$ write $`N_\sigma(T)`$ for the number of its odd zeros with $`\beta\ge\sigma`$ and $`|t|\le T`$. Use only those $`\sigma`$ with infinitely many such zeros, and set $`d(\sigma)=\limsup_{T\to\infty}\log(1+N_\sigma(T))/\log T`$. Define

```math
\theta_*=\sup_\sigma\frac{\sigma-1/2}{2-d(\sigma)},
```

with value zero for the empty supremum. If $`\theta\gt \theta_*`$, then $`R_\theta(x)`$ decays exponentially in $`x`$. If $`0\lt \theta\lt \theta_*`$, then $`R_\theta(x)`$ grows exponentially along a sequence of $`x\to\infty`$. In particular **a fixed exponent below one quarter makes this absolute tail tend to zero if and only if the eight odd factors have a fixed zero-free strip next to real part one at all sufficiently high ordinates**. This statement concerns the sum of absolute errors. Possible cancellation in the signed remainder is a different question. The endpoint $`\theta=\theta_*`$ is unnecessary for the equivalence. [Eight primitive factors](#the-fourth-root-formula)

For the upper bound take a finite partition of $`[1/2,1]`$ of mesh $`h`$. At a left endpoint $`\sigma\gt 1/2`$ with infinitely many zeros, $`N_\sigma(V)\ll_\varepsilon V^{d(\sigma)+\varepsilon}`$. Dyadic summation over $`|t|\gt T`$ costs $`O(T^{-2+d(\sigma)+\varepsilon})`$ before the real exponential. Finite half-strips drop out at large $`T`$. Multiplying by $`e^{(\sigma+h-1/2)x}`$ bounds its bin. Since $`d(\sigma)\le1`$, the definition of $`\theta_*`$ gives a uniform exponent gap at least $`\theta-\theta_*`$ before paying $`h+\theta\varepsilon`$. Choose both small. For $`\beta\le1/2+h`$ use the ordinary zero count to get $`O(e^{hx}\log(T)/T)`$. The finitely many bins give exponential decay.

For the reverse inequality choose a $`\sigma`$ whose quotient exceeds $`\theta`$. On $`\beta\ge\sigma`$ the factor $`|\sinh(az/2)|`$ is bounded below by a positive constant. Also

```math
|3-2\sqrt2\cosh(az)|\ge3-2\sqrt2\cosh(a(\beta-1/2))
\gg_\sigma 1-\beta .
```

The last inequality follows from the positive derivative of $`\cosh`$ on the interval from $`a(\sigma-1/2)`$ to $`a/2`$ and $`2\sqrt2\cosh(a/2)=3`$. The classical high-ordinate zero-free region gives $`1-\beta\gg1/\log(2+|t|)`$. Thus each such zero contributes $`|J(z)|\gg_\sigma |t|^{-2}/\log |t|`$ at large height. If $`d(\sigma)\gt 0`$, infinitely many dyadic shells $`V\lt |t|\le2V`$ contain at least $`V^{d(\sigma)-\varepsilon}`$ zeros for each sufficiently small positive $`\varepsilon`$. Otherwise summing shell counts would contradict the limsup. When $`d(\sigma)=0`$ infinitely many shells contain at least one zero. Take $`x=\log V/\theta`$. The contribution of these shells has exponential rate at least $`\sigma-1/2-\theta(2-d(\sigma)+\varepsilon)`$ up to a vanishing error. Small $`\varepsilon`$ makes this positive. [Explicit high-ordinate zero-free region](https://doi.org/10.1016/0022-314X%2884%2990089-1)

Finally the density estimate for every $`\epsilon\gt 0`$ gives $`d(\sigma)\le\min(1,(12/5)(1-\sigma))`$. If $`b=\limsup_{|t|\to\infty}\beta\lt 1`$, then for all relevant $`\sigma\le b`$ the difference between $`1/4`$ and the quotient is at least $`(1-b)/5`$. Indeed $`2-d\le2`$ and $`4(1-\sigma)-d\ge(8/5)(1-\sigma)`$. If $`b=1`$, there are infinitely many zeros in every fixed half-strip $`\beta\ge\sigma`$ with $`\sigma\lt 1`$. The quotient is at least $`(\sigma-1/2)/2`$, tending to one quarter. Its upper bound is one quarter by the same density estimate. This proves both directions. The preceding subexponential saving is compatible with either value of $`b`$. [Log-free density estimate](https://arxiv.org/pdf/1208.5992#page=8) · [Local cutoff algebra](FourthRoot.lean)
