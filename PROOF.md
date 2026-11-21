# An all-height deadline for both signs in the dyadic 547/727 prime race

## Statement

Put $`K=\mathbf Q(\zeta_{16})^+`$, $`F=K(\sqrt{397669})`$, $`E=K(\sqrt{-547},\sqrt{-727})`$, $`a=\log2`$, $`r_0=\sqrt2`$, $`c=2(r_0-1)^2`$ and $`B=r_0c`$. At each good prime ideal $`v`$ of $`F`$, put $`\epsilon_v=1`$ for splitting in $`E/F`$ and $`\epsilon_v=-1`$ for inertia. These signs equal the products of the normalized first-step theta trace jumps. Omit the places above $`2,547,727,343057911569,343070636977`$. The score width and dilation are both two. [Field and trace identity](https://github.com/leecrossley/theta-race-short-shifts/blob/483c31f2df2132ed2165620ea55e20e969228bfb/RaceTrace.lean#L83)

```math
W(t)=\max(0,\min(2-2/t,4/t-1)),\qquad
S(y)=\frac1{B\sqrt y}\sum_v-\epsilon_v\log(Nv)\,W(Nv/y).
```

**Theorem.** For every real $`H\ge1`$, there exist $`y_+,y_-\gt H`$ with $`2y_\pm\le10^{18}H`$ and

```math
S(2y_+)-S(y_+)\gt1.9735\log2,\qquad
S(2y_-)-S(y_-)\lt-1.9735\log2 .
```

The two starts may differ. The earlier all-height upper factor was $`e^{175}`$. A separate pair of certified one-way windows shows that each universal factor must exceed $`2e^6`$. [Earlier upper bound](https://github.com/leecrossley/theta-race-fixed-factor/blob/12b7c94cfb9586c25fb29e7c7aaf0a20adaf5b2f/PROOF.md) · [Lower obstruction](https://github.com/leecrossley/theta-race-factor-400/blob/f1773b07a55e76947d4e0c6ab19a112fd91fe6cf/PROOF.md)

## The full spectrum and the prime-only error

Write $`x=\log(2r_0y)`$ and $`f(x)=[S(2y)-S(y)]/a`$. The eight odd primitive Dirichlet factors of $`\zeta_E/\zeta_F`$ have conductors $`547,4376,8752,8752,727,5816,11632,11632`$. The family is closed under conjugation. The primitive Dirichlet strip theorem puts every nontrivial zero on the critical line. With multiplicities, define

```math
G(t)=\frac{3-2r_0\cos(at)}{c(t^2+1/4)},\qquad
\alpha(t)=\frac{2G(t)\sin(at/2)}a,\qquad
U(x)=\sum_\nu i\alpha(\nu)e^{i\nu x}.
```

A chosen positive ordinate and its negative partner contribute $`-2\alpha(\nu)\sin(\nu x)`$. The series is absolutely and uniformly convergent. For $`x\ge55`$ the explicit formula gives $`|f(x)-U(x)|\le E(x)`$, where

```math
E(x)=2\cdot10^6e^{-x/4}
  +128e^{-(x-2)/2}/(1-e^{-(x-2)})+22000e^{-3x/20}.
```

The terms cover the nontrivial base-field zeros, trivial zeros and prime powers. Each excluded norm is below $`e^{53}`$ and the log support has half-width $`3a/2\lt2`$, so the excluded primes have left the window. The decreasing bound satisfies $`E(100)\lt0.007`$. The previously certified 6,000 positive brackets select pairs whose omitted absolute $`\alpha`$-mass is below 1.44. A separate gamma-integral calculation bounds the same mass by 1.52. [Pointwise prime-only bound](https://github.com/leecrossley/theta-race-distant-interval/blob/527fffe7b7808457b47a292144fea96b556158a8/PROOF.md#from-the-profile-to-prime-ideals) · [Bracket and tail certificate](https://github.com/leecrossley/theta-race-factor-400/blob/f1773b07a55e76947d4e0c6ab19a112fd91fe6cf/PROOF.md) · [Independent primitive-theorem checks](https://math.mit.edu/~drew/grh/)

## Where the missing zeros can be

For $`d\gt0`$ and real $`z`$, use the completed odd functions $`\Lambda(s,\chi)=(q/\pi)^{(s+1)/2}\Gamma((s+1)/2)L(s,\chi)`$. At $`s=1/2+d+iz`$ their logarithmic derivatives give the positive Poisson sum

```math
\frac1d\Re\sum_\chi\frac{\Lambda'}{\Lambda}(s,\chi)
  =\sum_\nu\frac1{d^2+(\nu-z)^2}.
```

Indeed the difference of logarithmic derivatives at the two points reflected across the critical line is twice the real part by the functional equation. In the genus-one product every zero contributes the difference of its two reciprocal terms. The subtraction cancels the constant and gives the displayed positive kernel on division by $`2d`$. The completed odd factors are entire. Each retained positive ordinate is counted with its negative partner, including for nonreal characters in the aggregate family. [Poisson identity and functional-equation argument](https://github.com/leecrossley/theta-race-fixed-factor/blob/12b7c94cfb9586c25fb29e7c7aaf0a20adaf5b2f/PROOF.md)

At displacement $`d=1/2`$, take $`z=0,4,\ldots,48`$. Subtract the Poisson terms of both partners from each selected bracket using outward arithmetic. At every centre the remaining sum is below $`4/25`$. A missing ordinate in $`[0,48]`$ would be within distance 2 of a centre. It would contribute at least $`1/(1/4+4)=4/17\gt4/25`$. Hence every zero through height 48 is included with its multiplicity. Symmetry supplies the negative half. This also excludes an unlisted central zero. An independent check takes $`d=1`$ and centres $`0,3,\ldots,48`$. Each remaining sum is below $`1/5\lt1/(1+(3/2)^2)`$. [Thirteen half-unit displacement checks](check_spectrum.py) · [Seventeen unit-displacement checks](independent.py)

## A nonnegative compact average with one dominant frequency

Let $`h=3/32`$ and take the convolution of two uniform probability measures on $`[-h,h]`$. Convolve it with the symmetric atomic probability measure that puts total mass $`m_k/10002`$ at the two points $`\pm k/16`$, split equally between them. At zero their masses coincide. The integer file lists the 125 nonzero weights and $`\max k=307`$. The resulting probability $`\mu`$ has support in $`[-R,R]`$ for $`R=307/16+2h=155/8`$. Its even real Fourier transform is

```math
\phi(t)=\left(\frac{\sin(ht)}{ht}\right)^2
       \frac1{10002}\sum_km_k\cos(kt/16),\qquad \phi(0)=1 .
```

All masses are nonnegative. Their values were found by minimizing unwanted amplitudes in a finite linear program, followed by rounding and direct verification. The verifier reads only the integer certificate. [Positive filter](evidence/filter.json)

The carrier $`\gamma`$ is the positive conductor-547 zero in the third bracket, near $`3.28717829808674`$. Set $`A=2|\alpha(\gamma)\phi(\gamma)|`$ and let $`L`$ be the sum of $`2|\alpha(\nu)\phi(\nu)|`$ over the other 5,999 chosen positive ordinates. Rigorous interval calculations give $`A\gt4.917`$ and $`L\lt2.85`$. Every omitted frequency has $`|\nu|\gt48`$, so $`|\phi(\nu)|\le1/(48h)^2`$. Even the independent omitted-mass estimate leaves

```math
A-L-\frac{1.52}{(48h)^2}-E(100)
\gt 4.917-2.85-\frac{1.52}{(48h)^2}-0.007
\gt 1.9735 .
```

The first verifier uses the listed cosine mixture and the sinc directly. The second computes the real part of the polynomial $`\sum m_k z^k`$ by Horner's rule in blocks of 12 at $`z=e^{it/16}`$, computes the triangle factor from $`1-\cos(2ht)`$ and uses the half-angle version of $`G`$. Each encloses the whole frequency bracket before adding amplitudes. [Direct calculation](filter_bound.py) · [Independent calculation](independent.py) · [Replay transcript](evidence/check.txt)

## Both signs before the deadline

For each log centre $`X`$, the $`\mu`$-average of $`U(X+v)`$ equals the carrier sinusoid $`-2\alpha(\gamma)\phi(\gamma)\sin(\gamma X)`$ plus a remainder of absolute value at most $`L+1.52/(48h)^2`$. Interchanging sum and integral is valid by absolute uniform convergence. The carrier takes both values $`+A`$ and $`-A`$ on every closed interval of centres of length $`2\pi/\gamma`$.

For $`H\ge e^{100}`$, the eligible log abscissae form the interval $`(\ell,r]`$ with $`\ell=\log H+3a/2`$ and $`r=\log H+18\log10+a/2`$. Both averaging supports fit strictly inside it because

```math
2R+2\pi/\gamma+\log2\lt18\log10 .
```

Indeed the open interval $`(\ell+R,r-R)`$ has length exceeding one carrier period. Select one centre of each sign there. Every point of its averaging support exceeds 100, so the corresponding signed average of the prime-only score $`f`$ exceeds 1.9735. The measure is nonnegative and has mass one. At some point in its support the signed score must exceed 1.9735 as well. This gives the required strict signs and starting heights. [Local Lean deductions](FilteredDeadline.lean)

The published finite count covers every real $`1\le H\le10^9`$ with factor 18. Three pairs of individual log abscissae cover the remaining heights. The column abscissae are exact terminating decimals. [Finite coverage proof](https://github.com/leecrossley/theta-race-factor-18/blob/31ea5d495ca29f7f17869d1a3f321000266e4867/PROOF.md)

| Range of $`\log H`$ | Gain abscissa | Loss abscissa |
|---|---|---|
| $`[9\log10,55]`$ | 60.50 | 58.40 |
| $`[55,80]`$ | 83.66 | 88.00 |
| $`[80,100]`$ | 109.10 | 111.44 |

At each abscissa the signed true-bracket sum minus 1.44 and $`E(x)`$ exceeds 20, using the actual sine of the full interval product. For each row $`[u,v]`$ and each listed abscissa the verifier checks $`x-3a/2\gt v`$ and $`x-a/2\lt u+18\log10`$. These are precisely the strict start and end comparisons throughout that row. Together the three parts cover every $`H\ge1`$. No independence or simplicity assumption on the zero ordinates enters the proof. [Six-point certificate and height conversions](filter_bound.py)
