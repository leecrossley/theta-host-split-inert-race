# A fixed multiplicative deadline

## Statement

Put $`K=\mathbf Q(\zeta_{16})^+`$, $`F=K(\sqrt{547\cdot727})`$, $`E=K(\sqrt{-547},\sqrt{-727})`$, $`a=\log2`$, and $`B=2\sqrt2(\sqrt2-1)^2`$. For each base prime ideal $`v`$ let $`\epsilon_v=1`$ if it splits in $`E/F`$, and $`-1`$ if it is inert. Omit the places above $`2,547,727,343057911569,343070636977`$. The product of the normalized first-step theta trace jumps above a good prime gives $`\epsilon_v`$. Write

```math
W(u)=\max(0,\min(2-2/u,4/u-1)),\qquad
S(y)=\frac1{B\sqrt y}\sum_v-\epsilon_v\log Nv\,W(Nv/y).
```

Assume the generalised Riemann hypothesis (GRH) for the primitive Dirichlet $L$-functions in the 547/727 race and its base field, including the Riemann zeta function.

For every real $`H\ge1`$, there exist $`y_+,y_-\gt H`$ with $`2y_\pm\lt e^{175}H`$ and

```math
S(2y_+)-S(y_+)\gt 1.9735\log2,\qquad
S(2y_-)-S(y_-)\lt -1.9735\log2 .
```

The two locations can occur in either order. All logarithms are natural. [Field and trace identity](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-short-shifts/RaceTrace.lean#L83)

## Global input and transform

By the assumed GRH, every zero of a primitive Dirichlet $`L`$-function in the open critical strip has real part $`1/2`$. It includes the primitive character of modulus 1. The pinned formal statement of `LFunction_eq_zero_re_of_isPrimitive` uses Mathlib's `DirichletCharacter.LFunction`, which agrees with the ordinary Dirichlet series on $`\Re s\gt 1`$. This paper uses that global input as a hypothesis. Platt's finite-height result alone leaves the higher spectrum uncontrolled for this fixed window. [Pinned statement](https://github.com/leanprover-community/riemann-hypothesis/blob/371a1a26d5e17f604428e3a8ce6ce341cf517434/src/Statements.lean#L14) · [Mathlib identification](https://github.com/leanprover-community/mathlib4/blob/0df444a360eaa60ab8c11dca51a86af692955474/Mathlib/NumberTheory/LSeries/DirichletContinuation.lean#L75)

Write $`\chi_p=(\cdot/p)`$ and take the even quartic character $`\xi`$ of conductor 16 with $`\xi(5)=i`$. At primitive conductors the eight factors of $`L=\zeta_E/\zeta_F`$ have characters $`\chi_{547}\xi^j,\chi_{727}\xi^j`$, $`0\le j\lt 4`$. All are odd. The eight base factors of $`\zeta_F`$ have characters $`\xi^j,\chi_{547}\chi_{727}\xi^j`$, including ordinary zeta. Both lists are closed under conjugation.

Use the full prime sum for the moment. Put $`Q(x)=S(e^x/2)`$ and $`f(x)=[Q(x+a/2)-Q(x-a/2)]/a`$, the required change divided by $`a`$ for $`y=e^x/(2\sqrt2)`$. The bilateral transform on $`\Re z\gt 1/2`$ is $`J(z)D(1/2+z)`$ with

```math
G(z)=\frac{2\sqrt2\cosh(az)-3}{2(\sqrt2-1)^2(z^2-1/4)},\quad
J(z)=G(z)\frac{2\sinh(az/2)}a,\quad
D(s)=\frac{L'}L(s)-\frac{\zeta_F'}{\zeta_F}(2s)+D_R(s).
```

Direct Mellin integration gives $`G`$ with removable values filled in. Its kernel is nonnegative of mass 2 on $`[-a,a]`$. The kernel of $`J`$ has support $`[-j_0,j_0]`$, $`j_0=3a/2`$, and total variation at most $`4/a\lt 6`$. Euler coefficients give $`D(s)=-\sum_v\epsilon_v\log(Nv)(Nv)^{-s}`$. Prime squares cancel in the two logarithmic derivatives and residual powers have exponent at least three. The predecessor supplies the following two estimates, $`|J(u+it)|\le1000/(1+t^2)`$ for $`|u|\le1/2`$, and $`\sum_\rho(1+(\Im\rho)^2)^{-1}\lt 1000`$ across all sixteen factors with multiplicity. The latter follows from the completed logarithmic derivatives at 2 and the functional equation. [Transforms and absolute zero mass](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-waiting-bound/PROOF.md)

Let $`\nu`$ run over the zero ordinates of the eight odd factors with multiplicity. Then $`U(x)=\sum_\nu J(i\nu)e^{i\nu x}`$ converges absolutely and uniformly. For $`x\gt j_0`$ the explicit formula for $`f`$ contains $`U`$, the base terms $`-\tfrac12 J((\rho-1)/2)e^{(\rho-1)x/2}`$, the trivial-zero terms and $`f_R`$. The pole of ordinary zeta contributes zero because $`J(0)=0`$. To justify inversion, convolve first with four boxes of half-width $`\delta`$ and apply the predecessor's partial-fraction argument. Let $`\delta\downarrow0`$. The nontrivial sums have summable coefficient bounds and fixed shifted real parts 0 and $`-1/4`$. The trivial-zero series converges uniformly on compact subintervals of $`x\gt j_0`$. The original prime sum is continuous.

## Two inverse-square sums

The simple conductor-547 zero has ordinate $`3287178298086/10^{12}\lt \gamma\lt 3287178298087/10^{12}`$, with $`J(i\gamma)=iV`$ and $`V\gt 3.3428`$. The sign change, simplicity, absence of a coincident zero in the other factors, and nonvanishing of all eight factors at the centre are verified by intervals. Put $`A_\omega=\sum_{\nu\ne\omega}(\nu-\omega)^{-2}`$ for $`\omega=0,\gamma`$.

For an odd primitive character let $`\Lambda(s,\chi)=(q/\pi)^{(s+1)/2}\Gamma((s+1)/2)L(s,\chi)`$. Its genus-one product gives $`(\Lambda'/\Lambda)'(s)=-\sum_\rho(s-\rho)^{-2}`$. At $`s=1/2+i\omega`$ every summand is the positive inverse-square distance. At the selected zero first divide out its linear factor. For $`v_j=L^{(j)}(s,\chi)/j!`$, the regularized contribution of the $`L`$ part at that zero is $`2v_3/v_1-(v_2/v_1)^2`$. At an ordinary point it is $`2v_2/v_0-(v_1/v_0)^2`$. Add the second logarithmic derivative of the gamma factor. Rigorous complex balls give

```math
936\lt A_0\lt 937,\qquad 216\lt A_\gamma\lt 217 .
```

A second calculation uses first derivatives and positivity. At $`\sigma=1/1000`$ put $`P=\sigma^{-1}\Re\sum_\chi(\Lambda'/\Lambda)(1/2+\sigma+i\omega,\chi)-m/\sigma^2`$, where $`m`$ is the multiplicity at $`\omega`$. The functional equation and genus-one product give $`P=\sum_{\nu\ne\omega}p_\nu`$, $`p_\nu=(\sigma^2+(\nu-\omega)^2)^{-1}`$. Indeed, subtract the logarithmic derivatives at points reflected across the line. The difference is twice the real part and each partial fraction supplies $`2\sigma/(\sigma^2+(\nu-\omega)^2)`$. Since $`\sigma^2P\lt 1`$, $`A_\omega=\sum p_\nu/(1-\sigma^2p_\nu)\le P/(1-\sigma^2P)`$. This reproduces both upper bounds. [Curvature calculation](fixed_bound.py) · [Poisson calculation](poisson.py)

The mass of the $`G`$ kernel gives $`|J(it)|\le2|t|`$. Also $`|J(it)|\le7/2`$: use the first bound through $`7/4`$, interval evaluation of the trigonometric expression through 12, and $`2(3+2\sqrt2)/[2a(\sqrt2-1)^2(1/4+t^2)]`$ thereafter. To sharpen the central estimate, take five distinct positive ordinates, one in each interval $`(n,n+1)/10^9`$ for $`n=48536362,248572181,434216664,446011529,579615600`$. The checker reconfirms the crossings. Subtracting one copy of each positive and negative ordinate from the coarse bound gives

```math
\sum_\nu\frac{|J(i\nu)|}{\nu^2}
 \le \frac72 A_0-2\sum_{j=1}^5\frac{7/2-2\nu_j}{\nu_j^2}\lt 225.
```

For the selected test the corresponding sum over $`\nu\ne\gamma`$ is below $`(7/2)217`$.

## Fixed log window

Take $`X=\log H+115`$. Let $`\mu`$ be the law of $`x=X+v+w`$ for independent uniforms $`v,w`$ on $`[-30,30]`$. Its support is $`[\log H+55,\log H+175]`$. Its centred transform is $`\Phi(t)=(\sin(30t)/(30t))^2`$, with $`|\Phi(t)|\le1/(900t^2)`$. For $`M_\omega=\int f(x)e^{-i\omega x}\,d\mu`$, the odd-factor errors in $`M_0`$ and in $`M_\gamma-iV`$ are at most $`225/900`$ and $`(7/2)217/900`$.

For $`\lambda\gt 0`$, $`\int e^{-\lambda x}\,d\mu=e^{-\lambda(X-60)}[(1-e^{-60\lambda})/(60\lambda)]^2`$. Each nontrivial base zero has shifted real part $`-1/4`$. Their total averaged contribution is bounded by $`2\cdot10^6e^{-55/4}/225`$ for either moment, since $`1/(1+(\Im\rho/2)^2)\le4/(1+(\Im\rho)^2)`$ and the residue has a factor $`1/2`$. Each of the sixteen trivial-zero strings starts at shifted real part at most $`-1/2`$, with spacing at least 1. The kernel support and variation give a total below $`128e^{-53/2}/(1-e^{-53})`$.

At $`\sigma=7/20`$ the absolute coefficient mass of $`D_R=\sum_n c_n n^{-s}`$ satisfies

```math
\sum_n |c_n|n^{-\sigma}
 \le \frac{2}{1-2^{-\sigma}}\sum_v\frac{\log Nv}{(Nv)^{3\sigma}}
 \le \frac{16}{1-2^{-\sigma}}\left(-\frac{\zeta'}{\zeta}(3\sigma)\right)\lt 1600.
```

Here $`\sum_{v\mid p}\log Nv\le8\log p`$, $`Nv\ge p`$, and $`-\zeta'/\zeta(21/20)\lt 21`$. The two windows use $`n\le e^{x+j_0}`$ with weight per absolute coefficient at most $`[\sqrt2(e^{a/4}+e^{-a/4})/(Ba)]e^{-x/2}`$. Hence $`|f_R(x)|\lt 22000e^{-3x/20}`$, whose averaged error costs at most $`22000 e^{-165/20}/81`$. The three error bounds sum to less than $`1/10`$. The type primes split completely and all norms over $`2,547,727`$ are at most $`727^8`$. Every omitted norm is below $`e^{53}\lt e^{55-j_0}`$, so none enters the averaging window. [Type-prime splitting](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-short-shifts/SecondField.lean#L177)

The nonnegative tests $`1\mp\sin(\gamma x)`$ have masses at most $`1+1/90^2`$ because $`\gamma\gt 3`$. Their $`f`$ integrals are $`M_0\pm\Im M_\gamma`$. Subtracting both error budgets from $`V`$ leaves more than $`1.9735(1+1/90^2)`$. The first test forces a point with $`f\gt 1.9735`$, the second a point with $`f\lt -1.9735`$. Substituting $`y=e^x/(2\sqrt2)`$ gives both strict start and endpoint bounds. [Lean deductions](FixedFactor.lean) · [Expected output](evidence/check.txt)

The analytic estimates are written proofs. Lean checks the moment margins and the conversion to prime heights. The source verifier checks the named upstream theorem and external verification records. The full upstream kernel runs are Sutherland's.
