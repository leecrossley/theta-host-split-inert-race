# Two one-way windows in the dyadic 547/727 prime race

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

For arbitrarily large $`H`$ every eligible doubling by $`2e^6H`$ gains more than $`a/10`$. For arbitrarily large $`H`$ they all lose more than $`a/10`$. Hence no factor up to and including $`2e^6`$ can guarantee either move of magnitude greater than $`1.9735a`$ after every height.

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
