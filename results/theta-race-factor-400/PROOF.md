# Deadline obstructions in the dyadic 547/727 prime race

The theorem assumes Dirichlet GRH for the eight odd factors of the extension and the factors of the base field, including zeta. Alpöge's announced Dirichlet GRH theorem supplies this input. [Announcement](https://alpo.ge/remark.pdf) [Independent formal checks](https://math.mit.edu/~drew/grh/)

## Statement

Put $`K=\mathbf Q(\zeta_{16})^+`$, $`F=K(\sqrt{397669})`$, $`E=K(\sqrt{-547},\sqrt{-727})`$. For a prime ideal $`v`$ of $`F`$ the sign $`\epsilon_v`$ is $`+1`$ for splitting in $`E/F`$ and $`-1`$ for inertia. Omit the places above $`2,547,727,343057911569,343070636977`$. The product of the normalized first-step theta trace jumps gives this sign at each good place. [Field and trace lemma](https://github.com/leecrossley/theta-race-short-shifts/blob/483c31f2df2132ed2165620ea55e20e969228bfb/RaceTrace.lean#L83)

Write $`a=\log2`$, $`r_0=\sqrt2`$, $`c=2(r_0-1)^2`$ and $`B=r_0c`$. For $`b=r=2`$ use

```math
W(u)=\max(0,\min(2-2/u,4/u-1)),\qquad
S(y)=\frac1{B\sqrt y}\sum_v-\epsilon_v\log(Nv)\,W(Nv/y),\qquad
f(x)=\frac{S(2y)-S(y)}a,\quad x=\log(2r_0y).
```

Let $`I_+=[469417229.340,469417235.340]`$ and $`I_-=[484947092.225,484947097.850]`$. Their lengths are $`6`$ and $`45/8`$. Every terminating decimal here is an exact rational.

**Theorem.** One has $`f(x)\gt1/10`$ throughout $`I_+`$ and $`f(x)\lt-1/10`$ throughout $`I_-`$. The corresponding explicit starting heights are $`H_+=e^{469417229.340}/(2r_0)`$ and $`H_-=e^{484947092.225}/(2r_0)`$. For $`y\gt H_+`$ every doubling finishing by $`2e^6 H_+`$ gains more than $`a/10`$. For $`y\gt H_-`$ every doubling finishing by $`2e^{45/8}H_-`$ loses more than $`a/10`$. In particular the deadline $`400H`$ misses the required loss at the first height and the required gain at the second.

For arbitrarily large $`H`$ every eligible doubling by $`2e^6H`$ gains more than $`a/10`$. For arbitrarily large $`H`$ they all lose more than $`a/10`$. Any factor that guarantees a gain greater than $`1.9735a`$ after every height must exceed $`2e^6`$. The same lower bound holds for a factor that guarantees a loss of that magnitude after every height.

## Prime score and zero ordinates

The eight odd primitive factors of $`\zeta_E/\zeta_F`$ correspond to $`\chi_p\xi^j`$, $`p=547,727`$, $`0\le j\lt4`$, with $`\chi_p=(\cdot/p)`$ and $`\xi`$ the even quartic character of conductor 16 with $`\xi(5)=i`$. At each prime their conductors are $`p,8p,16p,16p`$. The zero verifier matches all eight Conrey labels on **every** unit class modulo $`16p`$, and checks their parity and conductor. The primitive Dirichlet critical-line theorem applies to the eight factors and to the base-field factors, including the primitive character of modulus one. The release checks the named theorem, its statement and its independent verification records. [Zero verifier](check_zeros.py) · [Upstream theorem](https://github.com/leanprover-community/riemann-hypothesis/blob/371a1a26d5e17f604428e3a8ce6ce341cf517434/src/Main.lean#L12) · [Sutherland's verification](https://math.mit.edu/~drew/grh/)

For real $`t`$ set

```math
g(v)=\frac4B\sinh\!\left(\frac{\max(0,a-|v|)}2\right),\quad
G(t)=\frac{3-2r_0\cos(at)}{c(t^2+1/4)},\quad
\alpha(t)=\frac{2G(t)\sin(at/2)}a .
```

Here $`G`$ is the Fourier transform of $`g`$ with exponential $`e^{-itv}`$ and is positive on the real line. The limiting profile and its even companion are

```math
U(x)=\sum_\nu i\alpha(\nu)e^{i\nu x},\qquad
V(x)=\sum_\nu G(\nu)e^{i\nu x}.
```

The sums run over all ordinates of the eight odd factors with multiplicity. Conjugation pairs a positive ordinate with a negative ordinate in the aggregate family. Both series converge absolutely and uniformly. The profile $`U`$ is real and odd, and $`V`$ is real and even. A chosen positive ordinate and its partner contribute $`-2\alpha(\nu)\sin(\nu x)`$ to $`U`$. The previously established pointwise estimate for this exact prime-only score is

```math
|f(x)-U(x)|\le 2\cdot10^6e^{-x/4}
 +\frac{128e^{-(x-2)/2}}{1-e^{-(x-2)}}+22000e^{-3x/20}
 \lt10^{-20}\qquad(x\ge1000).
```

The three terms bound the base-field zeros, trivial zeros and higher prime powers. The omitted norms have left the compact support at these heights. The verifier reevaluates the decreasing bound at 1,000. [Prime-only estimate](https://github.com/leecrossley/theta-race-distant-interval/blob/527fffe7b7808457b47a292144fea96b556158a8/PROOF.md#from-the-profile-to-prime-ideals)

## Four moments for the missing zeros

For $`u=\sin^2(at/2)`$ consider the cubic polynomial

```math
Q(u)=\frac{16u^3-48u^2+75u+11}{54}.
```

It majorizes $`\sqrt{u}`$ on $`[0,1]`$. Indeed for $`z=\sqrt{u}`$ the exact difference factors as

```math
Q(z^2)-z=\frac8{27}(z-1/2)^2(z-1)^2\big((z+3/2)^2+1/2\big)\ge0.
```

This gives zeros of the difference at both $`1/2`$ and $`1`$ in the variable $`z`$. Put $`V_k=V(ka)`$ and $`M_j=\sum_\nu G(\nu)\sin^{2j}(a\nu/2)`$. Elementary power-reduction identities give

```math
M_0=V_0,\quad M_1=(V_0-V_1)/2,\quad
M_2=(3V_0-4V_1+V_2)/8,\quad
M_3=(10V_0-15V_1+6V_2-V_3)/32.
```

The global absolute majorant is $`2(11M_0+75M_1-48M_2+16M_3)/(54a)`$. To calculate its moments, use the even compact-test explicit formula from the previous paper. With $`z=2^{-|j|}`$, write

```math
T_0=a,\qquad
T_j=\frac{(1+z)\log(1+z)+(1-z)\log(1-z)}{2\sqrt z}\quad(j\ne0).
```

The conductor correction is $`A_*=\sum_\chi\log(q_\chi/\pi)-8(\gamma_E+a)`$. Convolution with the odd gamma-factor kernel gives all four values.

```math
V_k=D_k-\frac8c\big[r_0(T_{k+1}+T_{k-1})-3T_k\big],
\qquad (D_0,D_1,D_2,D_3)=(A_*g(0),\,0,\,-a g(0),\,0).
```

To check the prime contributions directly, the support at $`x=ka`$ only permits $`2^{k-1}\le n\le2^{k+1}`$ for $`n=p^j`$. For an odd prime power up to 16 the sum over quartic twists vanishes unless $`n\equiv\pm1\bmod16`$. No such odd prime power contributes. At 2 the two odd-conductor characters have opposite values, so the coefficient is $`-2\log2`$ at even exponents and zero at odd exponents. The endpoints have weight $`g(a)=0`$. Only $`n=4`$ at $`k=2`$ survives, with contribution $`-2a\,g(0)/\sqrt4`$. The term $`A_*g(ka)`$ occurs only at $`k=0`$. Reflection supplies no further prime contribution at these abscissae. [Even explicit formula](https://github.com/leecrossley/theta-race-distant-interval/blob/527fffe7b7808457b47a292144fea96b556158a8/PROOF.md#a-positive-bound-for-all-omitted-ordinates) · [Cubic identity in Lean](PairedWindows.lean)

The table lists 6,000 positive rational brackets. The columns give the conductor, Conrey label and lower and upper numerators with denominator $`10^{14}`$. Strict opposite Hardy signs at the endpoints prove that each contains a zero. Brackets within a factor are pairwise disjoint. Choose one zero per bracket and its negative partner in the aggregate family. Denote these 6,000 positive choices by $`Z`$ and set $`u_\nu=\sin^2(a\nu/2)`$. A uniform upper bound for the absolute omitted part of $`U`$ is

```math
R=\frac2a\left(\sum_\nu G(\nu)Q(u_\nu)-2\sum_{\nu\in Z}G(\nu)Q(u_\nu)\right)\lt1.44 .
```

Positivity accounts for extra roots in a bracket, multiplicities and any missed roots at any height. No assumption on independence or simplicity enters the bound. [Bracket data](evidence/zero-brackets.tsv) · [Hardy function](https://flintlib.org/doc/acb_dirichlet/l-functions.html)

## Whole-interval certificates

For each bracket of centre $`t`$ and radius $`d`$, replace its chosen zero by $`t`$. The uniform cost on both windows is bounded by the sum of

```math
2\sup_{v\ \mathrm{in\ bracket}}|\alpha(v)-\alpha(t)|
   +2|\alpha(t)|\,d\,484947097.850 .
```

The bound follows from the Lipschitz constant 1 for the sine. The sum is less than 0.043. Retaining the signs of the two endpoints required narrowing the largest coefficient brackets to control phase error at these log coordinates. All final endpoint tests are replayed directly.

Take $`s=1`$ on $`I_+`$ and $`s=-1`$ on $`I_-`$. The signed midpoint sum $`-2s\sum\alpha(t)\sin(tx)`$ exceeds 1.64 at 8,193 evenly spaced nodes of **each** closed window. Its second derivative has absolute value at most $`2\sum|\alpha(t)|t^2\lt343000`$. With node spacing $`\Delta\le6/8192`$, linear interpolation costs at most $`343000\Delta^2/8\lt0.024`$ by the two-endpoint remainder formula. Including all errors gives

```math
sU(x)\gt1.64-0.043-0.024-1.44\gt0.13 .
```

[Four-moment verifier](paired_windows.py) · [Reproducible output](evidence/check.txt)

The second computation bounds the omitted part by Cauchy–Schwarz. Over the remaining multiset let $`A=\sum G(\nu)`$ and $`D=\sum G(\nu)\sin^2(a\nu/2)`$. Then $`\sum_{\mathrm{omitted}}|\alpha(\nu)|\le2\sqrt{AD}/a\lt1.52`$. The full moments used before subtracting the selected pairs are $`V_0`$ and $`(V_0-V_1)/2`$. This calculation recovers both from the original real integral with kernel $`e^{-t/2}/(2\sinh t)`$ acting on $`g(x-t)+g(x+t)-2g(x)`$. It uses 16,384 directed rectangles per dyadic piece from $`a/2^{16}`$ to $`32a`$. Below the first piece the integral is bounded in absolute value by $`[2\cosh(a/2)/B]a/2^{16}`$. Above $`b=32a`$ the bound is $`4g(0)e^{-3b/2}/[(3/2)(1-e^{-2b})]`$. The conductor and gamma constant before integration is $`\sum_\chi\log(q_\chi/\pi)+8\psi(3/4)`$. [Integral formula and tail estimates](https://github.com/leecrossley/theta-race-distant-interval/blob/527fffe7b7808457b47a292144fea96b556158a8/PROOF.md#two-whole-interval-calculations)

This second verifier encloses each actual frequency and each full time cell inside the sine, starting with 256 closed cells in each window. Bisection yields 2,194 covering cells for $`I_+`$ and 993 for $`I_-`$. It independently proves $`sU(x)\gt1/20`$ everywhere on the two intervals. The two calculations share the scale, field and sign-certified bracket table. [Direct cell verifier](cell_check.py)

## Heights and recurrence

For a closed interval $`[\ell,\ell+L]`$ put $`H=e^\ell/(2r_0)`$. Then $`y\gt H`$ and $`2y\le2e^LH`$ imply $`\ell\lt\log(2r_0y)\le\ell+L`$. The prime-only error on the two explicit intervals is less than $`10^{-20}`$, so their signed changes exceed $`a/10`$ in the stated direction. Since $`\log200\lt45/8\lt6`$, both explicit witnesses rule out 400. The local Lean deadline proof bounds both exponentials below by 22-term rational Taylor sums. [Lean deductions](PairedWindows.lean)

Absolute uniform convergence makes $`U`$ uniformly almost periodic. To see the required positive translations, truncate its series with uniformly small tail, and take two arbitrarily widely separated integer times whose finite phase vectors are arbitrarily close by compactness of the torus. Their difference translates the entire profile by arbitrarily little, uniformly in $`x`$. On $`I_+`$ the profile is above 0.13. On its reflection it is below $`-0.13`$ by oddness. Choose translations with uniform error below 0.01 and large enough to place either interval above 1,000. The resulting signed prime-only moves still exceed $`a/10`$. This proves both arbitrarily late assertions for the longer factor. The first missed height and the optimal universal factor remain open.

## A four-digit deadline and a shallow loss

Set $`\ell=3855250213.808`$, $`L=1247/200=6.235`$ and $`J=[\ell,\ell+L]`$. The corresponding height is $`H=e^\ell/(2\sqrt2)`$. Every normalized prime-only doubling move with log centre in $`J`$ exceeds $`-1.95`$. One such move exceeds $`12`$ and another is below $`-0.04`$. Thus the actual race changes in both directions inside the interval. Its losses all fall short of $`1.9735\log2`$ in absolute size. The two interior witnesses also finish by $`1000H`$.

Every eligible doubling $`y\gt H`$, $`2y\le 2e^L H`$ has its log centre in $`J`$. Since $`2e^{6.235}=1020.60123689658\ldots`$, even a factor of 1,020 fails. Odd reflection of the limiting profile and its arbitrarily large almost periods give arbitrarily late intervals missing either prescribed direction up to **the same** factor. The single numerical height above misses the large loss. Throughout this argument the Dirichlet critical-line input for the odd and base factors is the same as in the preceding sections. The arithmetic sign check on finite brackets alone cannot control zeros arbitrarily far up the strip. [Input and verification records](https://math.mit.edu/~drew/grh/) · [Previous factor-400 result](https://github.com/leecrossley/theta-race-factor-400/releases/tag/v1.0.0)

### A finite spectral sum with a uniform remainder

The extended table contains 13,946 genuine positive brackets in the eight odd factors. There are no overlaps within a factor. Each contributes one critical-line zero and a reflected partner in the aggregate. Use the same exact four moments and cubic majorant as above, subtracting the new table's lower contribution. Every missing root, including additional multiplicities, is covered by the resulting upper bound $`0.75`$ for the absolute omitted part of $`U`$. The first calculation uses the midpoint $`t`$ and radius $`d_t`$ of each bracket. Write

```math
\beta(t)=-4G(t)\sin(at/2)/a,\qquad
P(x)=\sum_t\beta(t)\sin(tx).
```

The difference between the retained true-frequency sum and $`P`$ is at most the sum of $`\sup|\beta(v)-\beta(t)|+|\beta(t)|(\ell+L)d_t`$ over the brackets. This is less than $`0.14`$ throughout $`J`$. All sine coefficients include both partners. The two estimates are replayed at 144-bit ball precision from the actual endpoints. [Extended brackets](evidence/thousand-brackets.tsv) · [Hardy sign replay](thousand_zeros.py) · [Moment and mesh verifier](thousand_window.py)

To evaluate an entire grid efficiently, put $`N=32768`$, $`x_c=\ell+L/2`$, $`x_k=x_c+kL/N`$ for $`-N/2\le k\le N/2`$. For each midpoint choose an integer $`j`$ nearest to $`tL/(2\pi)`$. Let $`\delta=tL/2-\pi j`$. The verifier checks $`|\delta|\lt8/5`$ on each computed midpoint. Exactly,

```math
e^{itx_k}=e^{itx_c}e^{2\pi ijk/N}\exp\!\left(i\delta\,\frac{2k}{N}\right).
```

Expand the last factor through degree 16. On the imaginary axis the absolute integral Taylor remainder is at most $`(8/5)^{17}/17!`$. After multiplication by the sum of $`|\beta(t)|`$ the error is less than $`10^{-7}`$. Group equal integers $`j`$ modulo $`N`$ at each Taylor order, take the complex ball discrete Fourier transform, and evaluate the resulting polynomial in $`2k/N`$. With the negative-exponent convention the bin at $`-k`$ supplies the positive exponential in the displayed identity. Although the endpoint Fourier bins coincide, their polynomial arguments are respectively minus one and one. Both endpoints are separately checked. [FLINT complex ball DFT](https://flintlib.org/doc/acb_dft/index.html)

The signed midpoint sum exceeds $`-1.06`$ at all 32,769 nodes after Taylor error. Linear interpolation between consecutive values costs at most $`(L/N)^2\sum_t|\beta(t)|t^2/8\lt0.01`$. At each endpoint of every grid cell the verifier subtracts the actual enclosures for all four errors, including the omitted tail and ordinate displacement. These combined inequalities give the stronger uniform conclusion $`U(x)\gt-1.94`$. The prime-only bound above is below $`10^{-20}`$ here, so $`f(x)\gt-1.95`$ on the closed interval. The use of the actual combined enclosure matters because some displayed component bounds have been rounded up independently. [Primary calculation](thousand_window.py) · [Verification transcript](evidence/check.txt)

### Direct cells and interior witnesses

The second proof evaluates the real gamma integrals for $`V(0)`$ and $`V(a)`$ with outward arithmetic. From the remaining nonnegative masses $`M'_0,M'_1`$ after subtraction of the selected pairs, Cauchy–Schwarz bounds the omitted absolute part by $`2\sqrt{M'_0M'_1}/a\lt0.83`$. It evaluates ordinary sines directly with both the full ordinate interval and the full log-time interval as inputs. Starting from 512 equal closed cells, it bisects a cell whenever its true-frequency sum minus the tail has not yet cleared $`-1.97`$. The 2,926 leaves cover $`J`$ and prove $`U(x)\gt-1.97`$ everywhere. This alone, with the prime-only error, rules out a loss at the required threshold. The two calculations share the certified input brackets and the ball arithmetic library. Their tail estimates and their treatment of the full interval are different. [Direct-cell verifier](thousand_cells.py)

The same calculation gives $`U(3855250215.346)\gt12.1`$ and $`U(3855250215.954)\lt-0.05`$. Their distances from $`\ell`$ are $`1.538`$ and $`2.146`$, both smaller than $`\log500`$. After the prime-only correction the corresponding normalized moves exceed $`12`$ and fall below $`-0.04`$. At the latter abscissa the move still exceeds $`-1.95`$ by the whole-interval bound. Continuity gives a zero between the two abscissae. This explains how the threshold can fail while the score does lose ground.

For recurrence, a uniform translation error below $`0.001`$ preserves the gap between $`-1.97`$ and the threshold $`-1.9735`$, even after adding the late prime-only error. Oddness gives the reflected interval on which the large gain is absent. Each orientation occurs arbitrarily far out by the same compact-torus argument as above. The new local Lean deductions check the exponential comparison, preservation of the signed margin and the interior sign implications. [Lean deductions](PairedWindows.lean)

### How often the thousand-fold obstruction returns

The same finite certificate gives a density bound without any assumption that the ordinates are independent. Let $`D=13946`$, $`M=2^{22}`$. For each selected true positive ordinate $`\nu`$ retain the phase of an integer shift $`n`$. Call $`n`$ good when every $`n\nu/(2\pi)`$ is at distance strictly less than $`1/M`$ from an integer. These shifts have lower uniform density at least $`M^{-D}=2^{-306812}`$ in the integers.

Indeed, the closure of the cyclic rotation is a compact subgroup of the $`D`$-torus. Partition the torus into $`M^D`$ half-open cubes of side $`1/M`$. One cube contains at least $`M^{-D}`$ of the subgroup's Haar mass. Subtracting one subgroup point in that cube puts its entire intersection into the open good set. Haar invariance supplies the lower bound. Continuous functions on the subgroup can be approximated uniformly by finite sums of characters. A nonconstant character there takes a value different from one at the generating rotation, so its averages over integer blocks tend to zero uniformly in the initial index by the finite geometric sum. Approximate the indicator of the open good set from below in Haar integral by continuous functions. This proves the stated lower uniform density, allowing all rational relations and torsion. [Earlier paired-density argument](https://github.com/leecrossley/theta-race-paired-density/releases/tag/v1.0.0)

The total absolute move coefficient of the selected pairs is below $`346`$ by the global cubic majorant. Hence a good shift changes the retained sum by at most $`346\cdot2\pi/M\lt0.001`$ uniformly on the entire real line. On the translated intervals use the direct-cell bound on the *retained* sum and the same absolute tail allowance. The omitted terms at the translated point can have any phases. For every sufficiently large good $`n`$, the interval $`n+J`$ still misses the large loss and contains both interior move signs. The paired interval $`n-J`$ still misses the large gain and contains both signs in the reversed orientation. The prime-only remainder is smaller than $`10^{-20}`$ on both by taking $`n`$ large. This argument controls both intervals for the same integer shift.

For the factor 1,000 one may also move the initial log height forward by any $`s\in[0,1/50]`$, since $`1/50+\log500\lt6.235`$. The interior witnesses and their reflections stay eligible. Each good integer therefore supplies a full log-height interval of length $`1/50`$ for each direction. The intervals associated to different integers are disjoint for a given orientation. Consequently the log heights missing a prescribed direction by the thousand-fold deadline have lower uniform continuous density at least $`2^{-306812}/50`$. This lower density is the limit inferior of the smallest relative measure over translated intervals as their lengths tend to infinity, with starting heights on the positive log-height half-line. [Numerical comparisons](thousand_window.py)

## An earlier height for factor 400

Put $`\ell=3895686.221`$, $`W=5.299`$, $`X=3895692`$, $`a=\log2`$, $`H_0=e^\ell/(2\sqrt2)`$, $`B_0=10^{12}`$ and $`\epsilon=2\cdot10^{-8}`$. Write $`\log^+v=\max(0,\log v)`$.

**Theorem.** If a doubling starts at $`y\gt H_0`$, finishes by $`2e^W H_0`$ and gains at least $`1.962\log2`$, one of the eight odd primitive factors has a zero $`\beta+i\gamma`$ satisfying

```math
8590\lt|\gamma|\le T=e^{1947900},\qquad
\beta\gt\frac12+\epsilon+
  \frac{49}{50X}\log^+\!\left(\frac{|\gamma|}{B_0}\right).
```

This implication is unconditional. In particular a gain greater than $`1.9735\log2`$ would force $`\beta\gt0.50000002`$, and $`\beta\gt3/4`$ if $`|\gamma|\ge e^{1000000}`$. With the announced critical-line theorem for primitive characters, every eligible gain falls below $`1.962\log2`$. As $`2e^W=400.27314659968\ldots`$, this starting height rules out 400. The finite interval verification below proves the implication independently of the global theorem. A direct finite exclusion of zeros throughout the indicated region would give the same prime-race conclusion. [Primitive theorem](https://github.com/leanprover-community/riemann-hypothesis/blob/371a1a26d5e17f604428e3a8ce6ce341cf517434/src/Main.lean#L12) · [Independent proof audit](https://math.mit.edu/~drew/grh/)

### Projecting the odd zeros

Use the bilateral transforms $`\mathcal G(z)= [3-2\sqrt2\cosh(az)]/[2(\sqrt2-1)^2(1/4-z^2)]`$ and $`J(z)=2\mathcal G(z)\sinh(az/2)/a`$. Thus $`\mathcal G(it)=G(t)\gt0`$ in the real-ordinate notation above. For each odd-factor zero write $`z=\beta-1/2+it`$. The odd nontrivial-zero contribution to the exact normalized move $`f(x)`$ is $`\sum J(z)e^{zx}`$. Define two absolutely convergent projected sums, retaining every ordinate with its full multiplicity,

```math
U_0(x)=\sum_\rho J(it)e^{itx},\qquad
V_0(s)=\sum_\rho G(t)e^{its}.
```

Conjugation of the family makes $`U_0`$ real and $`V_0`$ real and even. Write $`V_*(s)=\sum_\rho\mathcal G(z)e^{zs}`$. The compact even explicit formula in the four-moment proof evaluates **this exact sum** at $`s=0,a`$ with the same gamma integrals for general strip zeros. In the contour formula the test at a zero is the entire Laplace transform at $`z`$. The zeros can be off the line. At these two shifts the prime terms vanish as in the earlier arithmetic check. For the exact move one can first use positive smooth averages on $`[-h,h]`$ in the established prime explicit formula and let $`h\downarrow0`$. Compact supports stay inside positive log heights. The zero weights have an inverse-square summable majorant on each fixed window, and the physical prime score is continuous. This gives the pointwise formula, including the base-field term below. [Explicit formula and coefficient estimates](https://github.com/leecrossley/theta-host-split-inert-race/blob/theta-race-waiting-bound-v1.3.0/results/theta-race-waiting-bound/PROOF.md#the-transform-and-the-prime-powers) · [Even moments](#four-moments-for-the-missing-zeros)

### A bound for the high ordinates

For either group of eight primitive factors with conductors at most $`Q`$, use $`Q=11632`$ for the odd family and $`Q=16\cdot547\cdot727`$ for the base. If $`v\ge8590`$, then

```math
S_Q(v)=\sum_{|\gamma|\gt v}\frac1{1+\gamma^2}
\le \frac{80[\log(2Qv)+1]}{v-2}.
```

Here is a direct count. At $`s=2+iu`$ the real part of the Hadamard logarithmic derivative of a completed primitive factor gives a sum of $`(2-\beta)/[(2-\beta)^2+(u-\gamma)^2]`$. Its value is bounded by $`\tfrac12\log(q(|u|+2))+4.5`$. Indeed the Euler logarithmic derivative costs less than 2, and the two elementary factors for zeta cost at most $`3/2`$. For the gamma derivative use the series for $`\Re\psi(v_0+iw)`$, with $`1\le v_0\le3/2`$, $`w=u/2`$ and $`m=\max(1,\lceil|w|\rceil)`$. The first $`m`$ summands cost at most $`1+\log m`$. The tail is at most $`\sum_{n\ge m}[\tfrac12/(n+1)^2+w^2/(n+1)^3]\le1`$. This proves the stated allowance. Reflection $`\rho\mapsto1-\overline\rho`$ in the functional equation supplies the Hadamard constant just as at $`u=0`$. [Total reciprocal-square mass at zero](https://github.com/leecrossley/theta-host-split-inert-race/blob/theta-race-waiting-bound-v1.3.0/results/theta-race-waiting-bound/PROOF.md#total-mass-of-the-omitted-zeros)

At the midpoint of a unit interval each of its zeros contributes more than $`1/5`$. For a bin $`n\le|\gamma|\lt n+1`$, sum over both signs and eight factors. Since $`\log(Q(n+3))\gt9`$ for the two specified choices of $`Q`$ here, the count is at most $`80\log(Q(n+3))`$. For $`n\ge\lfloor v\rfloor`$ we have $`n+3\le2n`$. Integral comparison for the decreasing function $`\log(2Qu)/u^2`$ gives the displayed bound. Counting with multiplicity also covers characters with shared ordinates.

### What a displacement of an odd zero can change

Assume that the region in the theorem is empty. Platt’s finite theorem puts all zeros with $`|t|\le8590`$ on the line for these eight conductors. At greater ordinates up to $`B_0`$, reflection within a factor gives opposite real shifts $`\pm d`$ at the same ordinate, $`0\le d\le\epsilon`$. These two summands have the same oscillating factor $`e^{itx}`$. The following estimates refer to the average per member of the pair. [Platt, Theorem 7.1](https://arxiv.org/pdf/1305.3087#page=14)

For $`|t|\ge8590`$, the explicit kernels give $`|J(it)|\le50/(1+t^2)`$ and $`|\mathcal G(it)|\le40/(1+t^2)`$. On radii $`1/4`$ about the segment between $`it-\epsilon`$ and $`it+\epsilon`$ the real argument has absolute value below $`3/4`$ and the imaginary argument has absolute value at least $`|t|-1/4`$. In the explicit formulas the denominator factors $`z\pm1/2`$ have moduli at least this last number. The hyperbolic factors are bounded by the corresponding coshes of the real argument. This gives $`|\mathcal G|\le40/(1+t^2)`$ and $`|J|\le150/(1+t^2)`$ on the disks, hence both derivatives cost less than $`2000/(1+t^2)`$ on the segment. The bound 40 for $`\mathcal G`$ also holds throughout the critical strip at these ordinates.

After subtracting $`J(it)`$ the two kernel errors cost at most $`2000\epsilon e^{\epsilon x}/(1+t^2)`$ for $`0\le x\le X`$. The common coefficient multiplying the unperturbed kernel is $`\cosh(dx)-1`$. Thus the displacement through $`B_0`$ costs at most

```math
[2000\epsilon e^{\epsilon X}+50(\cosh(\epsilon X)-1)]S_{11632}(8590).
```

For $`B_0\lt|t|\le T`$ put $`p=49/50`$. Absence of zeros in the stated region gives $`e^{dx}\le e^{\epsilon X}(|t|/B_0)^p`$. The contribution to the displacement is at most $`(1000e^{\epsilon X}+50)D_p`$. An upper bound for the weighted sum over this band is

```math
D_p=\frac{160}{B_0}
 \left[\frac{A}{1-p}+\frac{p}{(1-p)^2}\right],
\qquad A=\log(2QB_0)+1,\quad Q=11632.
```

Indeed $`S_Q(v)\le160[\log(2Qv)+1]/v`$ there. Integrating the nonnegative weights by layers bounds their full weighted sum by $`S_Q(B_0)+(p/B_0^p)\int_{B_0}^\infty v^{p-1}S_Q(v)\,dv`$. For $`|t|\gt T`$ use the unrestricted strip bound $`(1000e^{X/2}+50)S_Q(T)\lt10^{-8}`$. Adding the three ranges gives $`|\sum J(z)e^{zx}-U_0(x)|\lt0.03`$ throughout the interval of interest.

For $`s=0,a`$ apply the reflected-pair estimate with 40 in place of 50 and $`a`$ in place of $`X`$. For all $`|t|\gt B_0`$ the unrestricted cost is at most $`40(e^{a/2}+1)S_Q(B_0)`$. Hence $`|V_*(s)-V_0(s)|\lt10^{-5}`$ at both shifts. [Outward-rounded error budgets](early_400.py)

### Controlling base-field zeros without a critical-line hypothesis

The base contribution for a nontrivial zero is $`-J((\rho-1)/2)e^{(\rho-1)x/2}/2`$. The pole at one vanishes because $`J(0)=0`$. Put $`Q_b=6362704`$ and $`T_b=e^{40}`$. The ordinates above $`T_b`$ in absolute value cost at most $`2000S_{Q_b}(T_b)`$, allowing for the halved argument and the factor $`1/2`$. Each base character has conductor dividing $`Q_b`$. Lifting it to modulus $`Q_b`$ only adds Euler factors with their zeros on $`\Re s=0`$.

McCurley's zero-free region for the full product modulo $`Q_b`$ shows that every zero with $`|\gamma|\le T_b`$, apart from a possible simple real exceptional zero, has $`(1-\beta)/2\ge \delta=[20\log(Q_bT_b)]^{-1}`$. His constant is $`9.645908801\lt10`$. The earlier total reciprocal-square mass for the sixteen factors is less than 1,000. Thus the nonexceptional base cost through $`T_b`$ is at most $`2\cdot10^6 e^{-\delta\ell}`$. [McCurley, explicit zero-free regions](https://www.sciencedirect.com/science/article/pii/0022314X84900891) · [McCurley, statement of Theorem 1.1](https://www.ams.org/mcom/1984-42-165/S0025-5718-1984-0726004-6/S0025-5718-1984-0726004-6.pdf)

For a real exceptional zero put $`d=(1-\beta)/2`$. The signed kernel of $`J`$ has total variation at most 6 and support in $`[-2,2]`$. Differentiation gives $`|J'(-v)|\le12e\lt40`$ for $`0\le v\le1/2`$. Since $`J(0)=0`$, this zero costs at most $`20d e^{-dx}\le20/(e x)`$. Even allowing one such zero per base factor gives $`160/(e\ell)`$. No lower bound for $`1-\beta`$ is needed. The higher prime-power remainder is at most $`10^6 e^{-x/10}`$ and trivial zeros cost at most $`128 e^{-(x-2)/2}/(1-e^{-(x-2)})`$. For each omitted rational prime at most $`1+4/\log2`$ exponents can meet a fixed log support of width four. The local logarithmic coefficients of either eight-factor product cost at most $`8\log p`$ per power. With the two level normalizations, all omitted-prime corrections are bounded by $`10^6 e^{-x/3}`$ at these heights. Together the base zeros, omitted primes, higher powers and trivial zeros cost less than $`0.001`$. [Arithmetic estimates and omitted primes](https://github.com/leecrossley/theta-host-split-inert-race/blob/theta-race-waiting-bound-v1.3.0/results/theta-race-waiting-bound/PROOF.md#bounding-square-contributions-without-the-base-field-zero-hypothesis)

### Two checks on the whole interval

Retain the 13,946 verified positive brackets and their conjugates. For one bracket containing $`t`$ the retained contribution is $`-4G(t)\sin(at/2)\sin(tx)/a`$. Subtract the retained weights from $`V_0(0)`$ and $`[V_0(0)-V_0(a)]/2`$ to get the nonnegative omitted masses $`m_0,m_1`$. They are the sums of $`G(t)`$ and $`G(t)\sin^2(at/2)`$ respectively over omitted members of the projected spectrum. The two exact gamma integrals above bound each mass from above after allowing $`10^{-5}`$ for its projection. Weighted Cauchy–Schwarz bounds the absolute omitted move by $`2\sqrt{m_0m_1}/a\lt0.831`$. This use of positivity concerns real ordinates after projection. [Gamma integral calculation](cell_check.py) · [True zero checks](thousand_zeros.py)

The direct calculation passes the full ordinate bracket and each closed time cell into the sine. It bisects 512 initial cells as needed. On all 1,062 leaves the retained sum is below 1.1. The Fourier calculation separately encloses 32,769 signed nodes, the degree-16 Taylor remainder, the interpolation between adjacent nodes and the displacement of each bracket from its midpoint. It also bounds the retained sum below 1.1 at every real point of $`[\ell,\ell+W]`$. The two computations use the same bracket table, kernels and ball arithmetic library. [Both interval checks](early_400.py) · [Fourier implementation](thousand_window.py) · [Check transcript](evidence/check.txt)

Adding every allowance gives $`f(x)\lt1.1+0.831+0.03+0.001=1.962`$. The endpoints $`y\gt H_0`$ and $`2y\le2e^WH_0`$ imply $`\ell\lt x=\log(2\sqrt2\,y)\le\ell+W`$. This proves the theorem by contraposition. The local Lean proofs check the deadline comparison and the numerical implication at the requested move size. The explicit formula, zero counts and kernel estimates are written proofs. [Lean deductions](PairedWindows.lean) · [Numerical verifier](early_400.py)

### A second witness below log height 1,620,543

A wider search identified the earlier closed interval $`[1620542.278,1620547.593]`$, at a different stretch of the same race. Put $`H_1=e^{1620542.278}/(2\sqrt2)`$. If any doubling after $`H_1`$ that finishes by $`2e^{5.315}H_1`$ gains at least $`1.944\log2`$, an odd-factor zero satisfies

```math
8590\lt|\gamma|\le e^{810330},\qquad
\beta\gt\frac12+4\cdot10^{-8}
 +\frac{49}{50\cdot1620548}\log^+\!\left(\frac{|\gamma|}{10^{12}}\right).
```

This is an unconditional implication. In particular the forced zero has $`\beta\gt0.50000004`$ and has $`\beta\gt3/4`$ if $`|\gamma|\ge e^{420000}`$. The announced primitive theorem therefore rules out the requested gain at $`H_1`$, even by the factor $`2e^{5.315}=406.72902625769\ldots`$. The starting log centre is over 289 times smaller than the previously published $`469417229.340`$ witness for factor 400. [Two-window verifier](early_400.py) · [Earlier witnesses](#statement)

Here the retained true-frequency sum has an interior peak which makes the two-moment allowance too large. Use the cubic majorant from the four-moment proof. With $`X=1620548`$ and $`\epsilon=4\cdot10^{-8}`$ the same paired-zero, weighted-band and unrestricted-tail estimates give an odd displacement below $`0.022`$. At the four even shifts $`0,a,2a,3a`$, insert $`3a`$ for the largest shift in the even-moment bounds above. Each exact moment then differs from the projected moment by less than $`2\cdot10^{-5}`$. The only interior prime contribution is at $`2a`$ from the prime power 4, as computed in the earlier cubic proof. Expanding the power-reduction identities, the total projected majorant before subtracting the retained pairs is

```math
\frac{71V_0(0)-42V_0(a)-6V_0(2a)-V_0(3a)}{54a}.
```

The error allowance for this quantity is only $`120(2\cdot10^{-5})/(54a)`$. After subtracting the sign-certified pairs, the omitted projected mass is below $`0.751`$. One calculation uses the closed gamma terms. The other encloses the two original real gamma integrals at $`0,a`$ and two additional integrals at $`ka`$, $`k=2,3`$. For these last two the interval is $`[(k-1)a,(k+1)a]`$, away from the singularity at zero. The integrand is $`e^{-t/2}g(ka-t)/(2\sinh t)`$, with the prime term at $`k=2`$ added separately. The computation encloses all 32,768 equal subintervals of each integral. [Four cosine moments](#four-moments-for-the-missing-zeros) · [Gamma integrals and both tail checks](early_400.py)

The retained sum is below $`1.17`$ on all 2,150 closed cells of the direct cover. The Fourier grid proves the same inequality at every point after bracket, Taylor and interpolation costs. The base and all other arithmetic terms still cost less than $`0.001`$ at this earlier $`\ell`$. Thus the normalized move is strictly below $`1.17+0.751+0.022+0.001=1.944`$ throughout. The width $`5.315`$ exceeds $`\log200`$. Both witnesses explicitly exclude a large gain. The separate four-digit example above excludes a large loss. [Exact margin and endpoint deductions](PairedWindows.lean) · [Replay transcript](evidence/check.txt)
