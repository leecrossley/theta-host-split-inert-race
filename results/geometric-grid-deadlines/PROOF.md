# Positive filters on geometric grids

## The two prime-ideal races

For $`i=A,B`$ take $`(p,q)=(547,727),(571,599)`$ respectively. Write $`K=\mathbf Q(\zeta_{16})^+`$, $`F_i=K(\sqrt{pq})`$, $`E_i=K(\sqrt{-p},\sqrt{-q})`$. A good prime ideal $`w`$ of $`F_i`$ has sign $`\epsilon_w=1`$ for splitting in $`E_i/F_i`$ and $`-1`$ for inertia. The fibre product of the normalized first-step theta trace jumps equals this sign. In each field omit the places above 2, its two odd conductors and its two type primes. For A the latter are 343057911569 and 343070636977. For B they are 413718278401 and 603339156001. [A trace and auxiliary primes](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-short-shifts/SecondField.lean#L177) · [B trace](https://github.com/leecrossley/theta-host-split-inert-race/blob/8fca7cad37185e3d7951f0dffff5e2442d1dfa24/HostAdmissible.lean#L80) · [B type primes](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Formalism/LocalRingCardinality/SpanPowerMembership/ExampleLayerExistsGenerators.lean#L229)

Let $`a=\log2`$, $`r=\sqrt2`$, $`c=2(r-1)^2`$, $`B_0=rc`$ and $`W(t)=\max(0,\min(2-2/t,4/t-1))`$. Our normalized observation $`f_i`$ is centred at $`x=\log(2r y)`$, so a grid step in $`\log y`$ is the same step in $`x`$.

```math
S_i(y)=-\frac1{B_0\sqrt y}\sum_w \epsilon_w\log(Nw)W(Nw/y),
\qquad f_i(x)=\frac{S_i(2y)-S_i(y)}a.
```

The eight odd factors are the quadratic characters of conductors $`p,q`$ twisted by the four powers of the even quartic character modulo 16. The eight even factors for the base field include ordinary zeta. We use Alpöge's announced primitive Dirichlet critical-line theorem for these factors. This is GRH for them, and the announcement has not yet been refereed. The profiles, the prime-only error and the remainder bound below all use it. Define

```math
G(t)=\frac{3-2r\cos(at)}{c(t^2+1/4)},\qquad
U_i(x)=\sum_{\nu} i_{\mathbf C}\frac{2G(\nu)\sin(a\nu/2)}a\,e^{i_{\mathbf C}\nu x}.
```

Here $`i_{\mathbf C}=\sqrt{-1}`$ and the zeros of all eight odd factors are counted with multiplicity. Conjugate factors pair opposite ordinates. The series converges absolutely and uniformly. One selected positive ordinate $`t`$ and its negative give $`\Re(c_t e^{i_{\mathbf C}tx})`$ with $`|c_t|=4G(t)|\sin(at/2)|/a`$. The sign and phase of the coefficient are incorporated in $`c_t`$. [Compact-test formula](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-factor-400/PROOF.md) · [Primitive theorem and independent checks](https://math.mit.edu/~drew/grh/)

## Error from counting primes only

For $`x\ge100`$ we use

```math
|f_i(x)-U_i(x)|\le
2\cdot10^6 e^{-x/4}+
128\,\frac{e^{-(x-2)/2}}{1-e^{-(x-2)}}+22000e^{-3x/20}
\lt 0.007.
```

The first field's derivation applies with the same constants to the second. Both base fields have degree eight. Every primitive conductor in either completed factorization is at most $`16\cdot547\cdot727`$. At $`s=2`$ the completed logarithmic derivative, including the two pole terms for ordinary zeta, gives a sum of real reciprocal zero distances less than 13 per factor using $`|\Re L'/L(2)|\lt2`$. Consequently $`\sum_\rho (1+(\Im\rho)^2)^{-1}\lt1000`$ across sixteen factors in either field. The dyadic kernel transform on $`|\Re z|\le1/2`$ is bounded by $`1000/(1+(\Im z)^2)`$. The shifted base residues have real part $`-1/4`$ and coefficient $`1/2`$, giving the first term. The base pole has zero dyadic difference. The sixteen trivial-zero strings, kernel radius $`3a/2\lt2`$ and total variation below six give the second term. [Reciprocal zero mass](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-waiting-bound/PROOF.md#total-mass-of-the-omitted-zeros) · [Full inversion and the three bounds](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-filtered-deadline/PROOF.md#the-full-spectrum-and-the-prime-only-error)

After the square correction, powers of order at least three have absolute Dirichlet mass below 1,600 at $`7/20`$, from $`\sum_{w\mid \ell}\log Nw\le8\log\ell`$ and $`-\zeta'(21/20)/\zeta(21/20)\lt21`$. The kernel then bounds their total by the third term. Each type prime is one modulo $`16pq`$ and therefore splits completely. The conductor primes have residue degree at most four, since the base Galois group has exponent four. Their norms and their squared norms are at most $`727^8\lt e^{53}`$. The type primes are below $`e^{28}`$, so even their squares lie below $`e^{56}`$. Thus any omitted first or second power has left the support at $`x\ge100`$. Higher powers are included in the preceding absolute bound. The decreasing right side is checked at 100. [Arithmetic checks](arithmetic.py) · [Higher powers estimate](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-fixed-factor/PROOF.md#fixed-log-window)

## The positive discrete average

Given an integer $`N\ge2`$, take the difference of two independent uniform integers from $`0,\ldots,N-1`$. Its law $`p_n=(N-|n|)/N^2`$ on $`|n|\lt N`$ has real nonnegative Fourier transform

```math
K_N(s)=\sum_n p_ne^{2\pi i_{\mathbf C}ns}
=\left|\frac{\sin(\pi Ns)}{N\sin(\pi s)}\right|^2,\qquad K_N(s)=1\quad(s\in\mathbf Z).
```

For degree $`m`$ the Fejér polynomial $`F_m(\theta)=|\sum_{j=0}^m e^{i_{\mathbf C}j\theta}|^2/(m+1)`$ is nonnegative. Its coefficients are $`u_j=1-|j|/(m+1)`$ for $`|j|\le m`$. At most two distinct selected ordinates $`\gamma,\eta`$ are used per cell. Their degrees are $`m,b`$, with $`b=0`$ indicating the constant second factor. For a log step $`h`$ and arbitrary centre $`X`$ set

```math
P_n=F_m(\gamma hn+\phi)F_b(\eta hn+\psi),\qquad
M=\sum_n p_nP_n .
```

For either desired sign $`\sigma`$ choose $`\phi=\arg(\sigma c_\gamma e^{i_{\mathbf C}\gamma X})`$ and, when present, the analogous $`\psi`$ for $`\eta`$. The $`(-1,0)`$ coefficient for the first carrier contributes exactly $`|c_\gamma|u_1`$ to the signed average. The second contributes $`|c_\eta|v_1`$ from $`(0,-1)`$, where $`v_l=1-|l|/(b+1)`$. Denote their sum by $`A`$. For every other combination bound the real part below by minus its absolute value. Put $`d_{jl}=j\gamma+l\eta`$ and define

```math
M^*=\sum_{j,l}u_jv_lK_N(hd_{jl}/(2\pi)),\qquad
L=\sum_{t\ {\rm selected}} |c_t|
\sum_{j,l}'u_jv_lK_N(h(t+d_{jl})/(2\pi)).
```

The prime on the inner sum omits exactly the resonant term already counted for each carrier. Coincident or nearly coincident sampled frequencies can remain in every other term. The mass satisfies $`0\le M\le M^*`$. If the unlisted profile has uniform absolute size at most $`R`$, nonnegativity gives the full signed average at least $`A-L-RM`$. Thus

```math
A-L\gt(\tau+R)M^*
\quad\Longrightarrow\quad
\max_{|n|\lt N}\sigma U_i(X+nh)\gt\tau .
```

Indeed the signed weighted sum minus $`\tau M`$ is strictly positive, hence some term with positive weight exceeds the threshold. The chosen phases depend on the centre and sign. No relation between different ordinates is assumed. [Positive-average deduction](GeometricGrid.lean) · [Primary evaluator](cover.py)

## Certifying all spacings

The [zero table](evidence/dyadic-zeros.tsv) is byte-identical to the preceding six-step certificate. The 14,000 genuine positive brackets are disjoint within each primitive factor. We select one ordinate from each, with its negative partner in the conjugate family. All other zeros, including additional multiplicities, contribute to the remainder. The earlier full sign transcript is pinned by its release. A standalone endpoint replay is included. [Earlier bracket proof](https://github.com/leecrossley/theta-host-split-inert-race/blob/consecutive-doubling-patterns-v1.0.0/results/consecutive-doubling-patterns/PROOF.md#every-six-step-word-on-the-dyadic-grid) · [Sign replay](check_grid_zeros.py)

Both calculations give $`R\lt1.63`$ in each field. The first subtracts the selected contributions from the full mass $`V_i(0)`$, which is nonnegative because GRH makes every ordinate real, and sine-square moment $`(V_i(0)-V_i(a))/2`$, where $`V_i(x)=\sum_\nu G(\nu)e^{i_{\mathbf C}\nu x}`$. The remaining amplitude is at most $`2\sqrt{R_0R_2}/a`$ by Cauchy–Schwarz. The independent calculation uses the polynomial $`(11+75u-48u^2+16u^3)/54\ge\sqrt u`$ for $`0\le u\le1`$. Its full moment is obtained from $`V_i(0),V_i(a),V_i(2a),V_i(3a)`$. The compact formula through norm 16 has no odd-prime contribution in either field. The two odd-conductor characters at 2 sum to zero at odd powers and two at even powers. One implementation evaluates the gamma convolution in closed form and the other by a series with a geometric bound on its remainder. [Moment formulas](moments.py) · [Derivation and polynomial identity](https://github.com/leecrossley/theta-host-split-inert-race/blob/consecutive-doubling-patterns-v1.0.0/results/consecutive-doubling-patterns/PROOF.md#every-six-step-word-on-the-dyadic-grid)

The cells use $`h=av`$ with integer endpoints for $`2^{21}v`$. Each row gives $`N`$, the two endpoints, a zero-based row index within the host, its degree, and the second index and degree. A second index of $`-1`$ denotes a single carrier. A covers the whole closed interval $`1\le v\le2`$ in 1,063 cells and B in 1,072, with $`N=2048`$. The primary check rounds each turn $`at/(2\pi)`$ to an integer tick within one tick of its enclosure, and rounds the absolute amplitudes outwards. A lower bound $`d`$ for the distance to the nearest integer over a whole cell gives $`K_N\le\min(1,(2Nd)^{-2})`$ by the chord bound for sine, with upper bound one at zero. The verifier uses integer modular distances and checked 64-bit dot products. Every cell satisfies $`A-L\gt(1.985+1.63)M^*`$. [A cells](evidence/cells-0.tsv) · [B cells](evidence/cells-1.tsv) · [Integer bounds](cover.py)

The independent calculation encloses the turn at both ends of every true zero bracket, takes the extreme products with the cell endpoints and tests whether the interval crosses an integer. It evaluates a directed sine table to bound the denominator of $`K_N`$ on each full interval. Separate radial and phase quantizations and the cubic spectral remainder reproduce every strict inequality. Two extra point rows check $`v=1`$ at $`N=216`$ and $`N=315`$, which give blocks of 431 and 629 starts. The sharper point calculation uses $`N=112`$ in A and $`N=160`$ in B, with both Fejér degrees equal to one. It retains the numerator of $`K_N`$. Alongside the sine quotient, an independent division-free evaluation at $`N=d2^j`$ uses the squared normalized odd Dirichlet polynomial for $`d=7,5`$ and the product of $`\cos^2(2^k\pi d s)`$ for $`0\le k\lt j`$. Both evaluations include all selected amplitudes. [Endpoint and sine-table verifier](independent_cover.py) · [Exact dyadic kernels](exact_dyadic.py)

## Grid points, deadlines and densities

Any consecutive block of $`2N-1`$ starts on a log grid of step $`h`$ is the support of the above law translated to its middle point. If its first start $`y`$ is at least $`e^{100}`$, then every centre $`\log(2r y)`$ is above 100. The prime-only error converts $`1.985`$ to a signed normalized move exceeding $`1.9735`$. This proves the two block lengths 223 and 319 at ratio two and the common block length 4,095 at every ratio between two and four.

For a smaller grid step $`0\lt h\lt a`$ take $`s=\lceil a/h\rceil`$. Then $`a\le sh\le2a`$. Every residue class of indices modulo $`s`$ supplies a subgrid with a certified step. Taking $`s=1`$ for $`a\le h\le2a`$ includes all required cases. Partition a long interval of indices in each residue class into whole blocks of 4,095. Each block contains each sign. Summing across the residue classes loses at most a bounded number of endpoint indices for each fixed $`s`$. Thus each sign has lower uniform natural density at least $`1/4095`$ on the full grid. The sharper dyadic densities follow from the same argument in one residue class. Lower uniform density here is the limit inferior of the least fraction in any long interval of nonnegative indices. A bounded initial segment has no effect.

For any grid offset choose the first grid point strictly above $`H'=\max(H,e^{100})`$. It is at most $`e^h H'`$. Take this as the first point of an appropriate subgrid block. The last eligible doubling finishes by $`2e^{h+4094sh}H'\le2e^{100}4^{4095}H`$ for every $`H\ge1`$. At ratio two this calculation also gives $`2e^{100}2^{223}H`$ in A and $`2e^{100}2^{319}H`$ in B. These are existence bounds for the specified sampling rule. [Stride and margin checks](GeometricGrid.lean) · [Constants](arithmetic.py)

## Why unbounded ratios are different

For any prescribed block length $`J`$ choose one profile seed for each of the four possible joint directions from the earlier four-step catalogue. At its first step both signed normalized profiles exceed 14.1. Uniform absolute convergence makes both profiles uniformly almost periodic along the integer dyadic shifts. In detail, truncate both series with uniformly small tails and take integer times whose vectors of retained phases are arbitrarily close. Compactness gives arbitrarily large positive differences $`m`$, with $`\sup_{i,x}|U_i(x+ma)-U_i(x)|\lt1/(2J)`$. One can take an arbitrarily distant pair in the same cell of a finite cover. [Published seeds](https://github.com/leecrossley/theta-host-split-inert-race/blob/consecutive-doubling-patterns-v1.0.0/results/consecutive-doubling-patterns/evidence/dyadic-starts.tsv)

At the ratio $`2^m`$, telescoping over the first $`J`$ samples changes either profile by less than $`1/2`$. A further arbitrarily large common near-period translates the four seeds into the positive prime range with profile error below $`1/10`$. The prime-only error then leaves every prescribed signed move in the four blocks above 13. The same arbitrarily large ratio works for all four choices, with a possibly different offset in each. Consequently no finite block length forces both directions uniformly over **all** grid ratios.
