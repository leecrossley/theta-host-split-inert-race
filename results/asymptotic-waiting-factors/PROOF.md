# Asymptotic waiting factors for the 547/727 prime race

## What repeated misses imply

Let $`F=\mathbf Q(\zeta_{16})^+(\sqrt{547\cdot727})`$ and $`E=\mathbf Q(\zeta_{16})^+(\sqrt{-547},\sqrt{-727})`$. Use the width-two tent score $`S`$ with the same exclusions and trace sign as the earlier paper. A qualifying move of normalized magnitude $`t\gt0`$ is a gain or loss of more than $`t\log2`$ from $`y`$ to $`2y`$. A factor $`C`$ works after $`H`$ if each orientation has a start $`y\gt H`$ with $`2y\le CH`$. Failure at arbitrarily large heights means that no upper cutoff removes all failures. [Exact score and trace convention](https://github.com/leecrossley/theta-host-split-inert-race/blob/theta-race-waiting-bound-v1.2.0/results/theta-race-waiting-bound/PROOF.md)

Take the strip zeros of the eight primitive odd characters $`(\cdot/p)\xi^j`$, with $`p=547,727`$, $`0\le j\le3`$, and $`\xi`$ the even quartic character of conductor 16. Put $`\sigma_*=\sup_\rho\Re\rho`$ and $`\Theta=\sigma_*-1/2`$. Assume first that one of these factors has an off-line zero. Reflection then gives $`\Theta\gt0`$. Define

```math
S_*=\sum_{\gamma\in\Gamma_*}\frac1{\gamma},\qquad
\Gamma_*=\{\gamma>0:\sigma_*+i\gamma\ \mathrm{is\ a\ zero}\},
```

using each distinct ordinate once, even where several characters share it. An empty sum is zero.

**Theorem.** For each $`W\gt\pi S_*`$ and $`0\lt\eta\lt\Theta`$, every sufficiently large real $`H`$ has a gain and a loss of normalized magnitude greater than $`H^\eta`$ on two doublings with $`H\lt y`$ and $`2y\lt2e^W H`$. This implies the same conclusion at every fixed positive normalized threshold. In particular if $`\sigma_*`$ is unattained, any $`C\gt2`$ works eventually. The number 2 is the geometric infimum since a doubling beginning strictly after $`H`$ finishes strictly after $`2H`$.

### Spectral estimates and summable radii

For $`x=\log(2\sqrt2\,y)`$ write $`f(x)=(S(2y)-S(y))/\log2`$. In the absolutely convergent odd zero series the coefficient at $`z=\rho-1/2`$ is $`mJ(z)`$, where $`m`$ is the combined multiplicity and the move kernel satisfies

```math
J(z)=\frac{2\sqrt2\cosh(az)-3}{2(\sqrt2-1)^2(z^2-1/4)}
       \frac{2\sinh(az/2)}{a},\qquad a=\log2.
```

It is nonzero for $`0\lt\Re z\lt1/2`$. The tent numerator vanishes only on $`\Re z=\pm1/2`$, and the difference multiplier has its zeros on the imaginary axis. The family is closed under conjugation and reflection in the critical line. There are no real off-line zeros. For each strip zero the bound $`|J(z)|\le1000/(1+(\Im z)^2)`$ and the total zero-count bound in the previous proof give a finite sum of absolute coefficients. The signed prime-square correction, other powers and excluded primes are bounded at large heights. The trivial-zero terms tend to zero. Thus multiplying by $`e^{-dx}`$ for *any fixed* $`d\gt0`$ sends the arithmetic and trivial contributions to zero, uniformly on every compact window translated to infinity. The smoothed explicit formula may also be applied directly against each of the compact smooth tests below. [Transform, absolute mass and arithmetic correction](https://github.com/leecrossley/theta-host-split-inert-race/blob/theta-race-waiting-bound-v1.2.0/results/theta-race-waiting-bound/PROOF.md#the-tilted-prime-move)

For $`1/2\lt\sigma\le\sigma_*`$ let $`S(\sigma)`$ sum $`1/\gamma`$ over distinct **complex points** $`\beta+i\gamma`$ of these factors with $`\beta\ge\sigma`$ and $`\gamma\gt0`$. For fixed $`\sigma\gt1/2`$, Montgomery’s density bound has power $`3(1-\sigma)/(2-\sigma)\lt1`$ in the height. Summing its dyadic bound divided by the lower endpoint gives $`S(\sigma)\lt\infty`$. As $`\sigma\uparrow\sigma_*`$, these sums decrease to $`S_*`$ by domination by any one earlier summable series. If the supremum is attained take $`\sigma=\sigma_*`$. Otherwise take the real part of an actual zero sufficiently close to $`\sigma_*`$. In both cases choose $`g\gt0`$ with $`\sigma+ig`$ a zero and ensure

```math
d=\sigma-1/2>\eta,\qquad \pi S(\sigma)\lt W .
```

No equality or spacing between other ordinates is assumed. [Montgomery](https://www.mathnet.ru/eng/mat564) · [The density estimate in full](https://www.numdam.org/item/10.5802/jtnb.1162.pdf#page=11)

### A positive average with zero limiting mean

Tilt the signal to $`F_d(x)=e^{-dx}f(x)`$. The exponent of a zero $`\beta+i\gamma`$ is now $`z=(\beta-\sigma)+i\gamma`$. For one positive ordinate with $`\beta\ge\sigma`$, set $`b=\pi/(2\gamma)`$, $`u=\beta-\sigma`$. Give the points $`-b,b`$ weights proportional to $`e^{ub},e^{-ub}`$. Their probability has bilateral transform $`\cosh(b(z-u))/\cosh(bu)`$, zero at $`u\pm i\gamma`$. One factor kills both conjugates. Convolve one factor for each distinct positive complex point in the half-strip. Its radius costs just $`b`$. The sum of the radii is $`r=\pi S(\sigma)/2`$. Realizing the sum of independent draws constructs a compact probability $`P`$ on $`[-r,r]`$. Each target moment vanishes exactly by independence. This argument covers a finite or infinite list. In the unshifted real-frequency setting this is the diameter bound of Kozma and Oravecz for trigonometric polynomials. [Kozma and Oravecz, Theorem 1](https://arxiv.org/pdf/math/0601021)

Choose a small $`\epsilon\gt0`$ so that $`2r+2\epsilon\lt W`$. Smooth $`P`$ by a nonnegative smooth probability supported on $`[-\epsilon,\epsilon]`$ and positive inside. Denote the new density by $`K`$. All its target moments still vanish. It has integral one and support of diameter less than $`W`$. Somewhere it is strictly positive, hence there is a closed interval with nonempty interior on which $`K\ge c\gt0`$. This fact suffices even when the unsmoothed probability had a gap about its centre.

Translate $`K`$ to centre $`X`$. Its average of $`F_d`$ tends to zero as $`X\to\infty`$. Every term with $`\beta\ge\sigma`$ vanishes. Those with $`\beta\lt\sigma`$ tend individually to zero and are dominated by their summable absolute coefficients times a fixed constant. This is uniform for all sufficiently large real $`X`$. The same argument controls the arithmetic corrections. [Positive-filter precedent](https://github.com/leecrossley/theta-host-split-inert-race/blob/theta-race-waiting-bound-v1.2.0/results/theta-race-waiting-bound/PROOF.md#the-positive-two-point-filter)

### A signed test inside that average

We need a complex smooth density $`Q`$ supported strictly inside the interval on which $`K\ge c`$, whose bilateral transform vanishes at *every* normalized exponent with $`\beta\ge\sigma`$ except $`z_0=ig`$, and is nonzero at $`z_0`$. Its sign is free.

Choose a smaller open interval centred at some $`t_0`$ inside that positive interval. Take a cutoff $`G\gt2g`$ so large that the sum of the radii of the preceding two-point probabilities for $`\gamma\gt G`$ is smaller than both a quarter of the open interval's length and $`1/(8g)`$. Their convolution $`P_{\rm tail}`$ still kills all of those large ordinates. It has nonzero moment at $`ig`$ because $`\cos(gx)\ge\cos(1/8)\gt0`$ on its support. Choose a nonzero nonnegative smooth bump $`b_0`$ of integral one supported sufficiently near zero that $`\widehat b_0(ig)\ne0`$ for the same reason, and that its support plus the tail support fits strictly inside the chosen interval after translation by $`t_0`$.

Only finitely many distinct normalized exponents with $`|\gamma|\le G`$ remain. Include both signs of the ordinate. Form the complex polynomial

```math
p(z)=\prod_{\zeta\ne z_0}\frac{z-\zeta}{z_0-\zeta},
```

where the product ranges over those remaining exponents other than $`z_0`$. This includes $`-ig`$ and also a different real part at the *same* ordinate if present. Poles at the exact same complex point have already been combined. Apply $`p(-d/dx)`$ to the bump, convolve with $`P_{\rm tail}`$ and translate by $`t_0`$. Differentiation preserves the support and integration by parts multiplies the bilateral transform by $`p(z)`$. The resulting bounded smooth compact density $`Q`$ has exactly the required cancellations and $`\widehat Q(ig)\ne0`$. Extra zeros elsewhere are harmless.

With $`Q`$ translated to the same $`X`$, its average of $`F_d`$ is

```math
A e^{igX}+o(1),\qquad A=mJ(d+ig)\widehat Q(ig)\ne0 .
```

Indeed only the selected positive ordinate survives on the chosen half-strip. Smaller real parts and the arithmetic correction tend to zero by the same summable majorant, now times the fixed total variation of $`Q`$. This estimate requires no rate of approach of the smaller real parts to $`\sigma`$.

Fix $`\lambda\gt0`$ small enough that $`\lambda|Q(x)|\le K(x)/2`$ everywhere. For each $`X`$ choose a unit complex phase $`w_X`$ that makes the real part of the signed average positive and equal to its modulus. Use the two densities

```math
K(t)\ \pm\ \lambda\Re(w_X Q(t)).
```

Both are nonnegative, have mass at most $`3/2`$, and use the original support of diameter less than $`W`$. Their signed integrals after translation approach $`\pm\lambda|A|`$ in the sense of an error tending to zero. Consequently each orientation has a point where $`\pm F_d`$ exceeds a fixed positive constant. Place the common compact support strictly inside $`(\log(2\sqrt2 H),\log(2\sqrt2 H)+W)`$. This is possible for every $`H`$, and the translation tends to infinity with $`H`$. At each witness $`e^{dx}\ge(2\sqrt2 H)^d`$, which eventually beats $`H^\eta`$ after multiplication by the fixed constant. The starts are strictly after $`H`$ and the completions strictly before $`2e^W H`$. [Algebraic and integral checks](check_windows.py)

### The numerical consequence for 400

If some zero is off the line and 400 fails arbitrarily late at any fixed threshold, the theorem implies $`S_*\ge\log(200)/\pi`$. Platt’s verification through ordinate 8590 holds for all eight conductors, the largest being 11632. All rightmost ordinates therefore exceed 8590. With finitely many such ordinates their number $`N`$ satisfies $`N\gt8590 S_*`$. The independent arithmetic checks give

```math
14487<8590\log(200)/\pi<14488 .
```

Thus $`N\ge14488`$, with the assertion automatic for infinitely many ordinates. Equivalently, in the off-line case at most 14,487 distinct positive rightmost ordinates guarantee that 400 works eventually. This is an inverse statement about *repeated* misses. The finite example in the earlier publication uses a critical-line premise for the pointwise comparison. [Platt, Theorem 7.1](https://arxiv.org/abs/1305.3087) · [Arithmetic transcript](evidence/check.txt)

### One positive ordinate

Suppose the rightmost line has just one positive ordinate $`g`$ with combined nonzero coefficient $`A_0`$. Dominated convergence gives the uniform asymptotic

```math
e^{-\Theta x} f(x)=2\Re(A_0e^{igx})+o(1).
```

The supremum assertion above applies for every $`W\gt\pi/g`$. Conversely for $`0\lt W\lt\pi/g`$ take a closed interval of width $`W`$ strictly inside a positive half-wave, so the cosine there has a uniform positive lower bound. Translate by integer multiples of the period. The entire move is then positive on arbitrarily late such intervals, precluding the required loss. Negative half-waves likewise preclude a gain. This proves that the infimum of all eventually valid factors is exactly $`2e^{\pi/g}`$ for every fixed positive amplitude. At $`C\le2`$ no doubling can complete at all. Since $`g\gt8590`$, the infimum is strictly less than 2.000732. The proof gives no assertion at the boundary factor itself. [Rational and ball bounds](check_windows.py)

### A short reversal detects an off-line zero

Fix $`C_0=1601/800=2.00125`$ and the threshold $`t=1.9735`$. Call $`H`$ successful if a normalized gain above $`t`$ and a normalized loss below $`-t`$ both start after $`H`$ and finish by $`C_0 H`$. There is no prescribed order. The following three assertions are equivalent.

1. At least one of the eight odd Dirichlet factors has a zero off the critical line.
2. There are successful arbitrarily large real $`H`$.
3. Successful heights have positive lower logarithmic density, meaning that the lower limit of $`(\log X)^{-1}\int_{[1,X]\cap A}dH/H`$ is positive for their set $`A`$.

In the off-line case, the positive-density conclusion extends to **every** $`C\gt2`$ and every fixed positive threshold. For any fixed $`0\lt\eta\lt\Theta`$ it still holds if the threshold is replaced by $`H^\eta`$. The density can depend on the factor and exponent. Furthermore, every $`C\gt2`$ succeeds at every sufficiently large height at the original threshold if and only if the supremum $`\sigma_*`$ is **unattained**. Thus there are three distinct possibilities for short windows.

**All eight factors on the line.** By the prime ideal theorem for the fixed base field, the square contribution to the normalized level tends to a constant. Taking a dyadic difference removes it. Higher powers and bad primes give $`o(1)`$. The actual move is therefore an absolutely convergent sine series $`U(x)`$ plus $`o(1)`$ as $`x\to\infty`$. This only uses the eight odd critical-line assertions. The prime ideal theorem for the square contribution is unconditional. Uniform continuity already gives, for every fixed positive threshold, a sufficiently small width on which two opposite qualifying moves are impossible at all sufficiently late heights. [Prime-square argument](https://github.com/leecrossley/theta-host-split-inert-race/blob/theta-race-waiting-bound-v1.2.0/results/theta-race-waiting-bound/PROOF.md#every-positive-excess-exponent)

For the stated width write $`a=\log2`$, $`r=\sqrt2`$. The earlier four-moment certificate bounds the omitted absolute sine-coefficient mass by $`0.75`$ after retaining 13,946 positive ordinates with their partners. Its cosine moments use just the eight odd factors. The independent Cauchy–Schwarz certificate gives $`0.83`$, which alone would not close the inequality below. The earlier removal of the finite-height square bias used extra zero information. Here the limiting square-bias cancellation just proved suffices. For a retained positive ordinate $`g`$ the absolute sine coefficient is

```math
B(g)=\frac{2(3-2r\cos(ag))|\sin(ag/2)|}
             {a(r-1)^2(g^2+1/4)} .
```

The verifier encloses $`\sum_g g B(g)\lt3805`$ on the original disjoint-within-character brackets. A second evaluation uses $`3-2r\cos v=(r-1)^2+4r\sin^2(v/2)`$ and the supremum of $`2|\sin(gh/2)|`$ on the allowed displacements. For $`0\le h\le\log(1601/1600)\lt1/1600`$, the first bound alone gives

```math
|U(x+h)-U(x)| \lt 3805/1600+2(0.75)\lt3.89\lt 2(1.9735).
```

The $`o(1)`$ error preserves this strict inequality for all sufficiently large $`x`$. Two starts in the same requested window have log separation at most $`\log(1601/1600)`$, so none of those $`H`$ can be successful. The ordinate file here is a lossless copy of the exact published bracket file. Its hash is checked after decompression. The earlier release replays the individual sign brackets and the tail calculation. [Original zero and tail certificate](https://github.com/leecrossley/theta-host-split-inert-race/blob/theta-race-factor-400-v1.1.0/results/theta-race-factor-400/PROOF.md#a-four-digit-deadline-and-a-shallow-loss) · [Variation calculation](check_reversals.py)

**An unattained supremum.** All sufficiently late $`H`$ succeed for every $`C\gt2`$ and the growing threshold by the theorem above, because $`S_*=0`$.

**An attained supremum off the line.** With $`d=\Theta`$, uniform dominated convergence gives $`e^{-dx}f(x)=U_*(x)+o(1)`$, where the real function $`U_*`$ is the absolutely convergent trigonometric series on the extremal ordinates. It has zero mean and a nonzero coefficient, since multiplicities at the same complex point add and $`J`$ has no zeros there. Orthogonality in the ordinary long-interval mean recovers that coefficient. Uniform almost periodicity follows by truncation to finite sums. More explicitly, the closure of the finite phase flow is a compact group. Finitely many orbit translates of any neighbourhood of its identity cover that group. Thus the times returning near the identity have bounded gaps. The summable tail turns these into uniform approximate periods of $`U_*`$. In particular a nonnegative, nonzero such function has strictly positive mean, by recurring positive intervals. Applying this to both signs shows that $`U_*`$ takes both positive and negative values.

The function $`U_*`$ cannot vanish identically on any open interval. Choose a positive extremal ordinate $`g`$. Place the complex probe constructed above inside any such interval. Its tail measure can be given arbitrarily small support because $`\sum1/\gamma\lt\infty`$. The finite differential multiplier removes every remaining exponent including $`-ig`$, while the factor at $`ig`$ is nonzero. Integrating the absolutely convergent series against that probe would give a nonzero value. This proves the assertion even for infinitely many extremal ordinates.

There are therefore points of opposite strict sign arbitrarily close together. Indeed on a finite interval containing both signs, a positive distance between the positive and negative sets would leave an open zero interval between them. For any $`W=\log(C/2)\gt0`$, put two such points strictly inside an interval of width $`W`$. Their signed values exceed some common $`c\gt0`$ in magnitude. There is slack at both window endpoints. Uniform approximate periods with error below $`c/2`$ recur with bounded gaps. Around each translated window, an interval of starting log heights of fixed positive length still contains the translated pair in its interior. Such intervals have positive lower density on the log axis. Select one in the middle of every sufficiently large block of a fixed length, leaving room at both ends. The $`o(1)`$ leaves both signed normalized values strictly bounded away from zero after tilting. Undoing the tilt gives both actual magnitudes above $`H^\eta`$ for every sufficiently late such interval. This proves the third assertion in the attained case.

Finally $`U_*`$ is strictly positive on some closed interval of nonzero length. Choose a shorter width inside it and leave slack at both ends. Uniform approximate periods produce log intervals of starting heights with bounded gaps and a common positive length on which the *entire* window remains positive. The tilted error tends to zero there. Thus for some $`C\gt2`$ a loss fails on a set of positive lower logarithmic density. A negative interval likewise gives missing gains. This proves the converse in the every-height characterization. All these assertions concern the actual prime-only race. No independence or simplicity assumptions on the zeros enter. [Compact interpolation](#a-signed-test-inside-that-average) · [The earlier growth criterion on sparse grids](https://github.com/leecrossley/theta-host-split-inert-race/releases/tag/compact-filter-oscillation-v1.0.0)

### Both signs in a window ending by 2.011H

**Unconditionally, a set of starting heights of positive lower logarithmic density has both qualifying moves before $`2.011H`$.** More precisely that set contains log intervals of some fixed positive length with bounded gaps at sufficiently large heights. The two actual doublings can start less than $`0.55\%`$ apart in relative height. The statement allows the other starting heights to miss. Combined with the criterion above, the same test at $`2.00125H`$ distinguishes the two zero-location cases.

If any odd factor has an off-line zero, this is the positive-density conclusion just proved, including the bounded gaps between intervals of successful log heights. Suppose all eight odd factors satisfy the critical-line assertion. Let $`V(x)=\sum_\nu G(\nu)e^{i\nu x}`$ be the even level profile from the earlier explicit formula. Its difference $`U(s)=[V(a/2+s)-V(a/2-s)]/a`$ is the odd dyadic profile. At $`s=273/100000`$, put $`x=a/2\pm s`$. Both $`x`$ are strictly between zero and $`a=\log2`$ and satisfy $`e^{x+a}\lt3`$. The only possible forward prime power is 2. Its two odd-conductor coefficients cancel and the other characters are ramified there. No reflected prime power is in support. Write $`r=\sqrt2`$, $`c=2(r-1)^2`$, $`g(x)=4\sinh((a-x)/2)/(rc)`$ on this interval and $`A_*=\sum_\chi\log(q_\chi/\pi)-8(\gamma_E+a)`$. The even explicit formula then gives

```math
V(x)=A_*g(x)-\frac8c\{r[T(x+a)+T(x-a)]-3T(x)\}.
```

With $`z=e^{-|v|}`$ the two equivalent expressions used for the gamma kernel are

```math
T(v)=\frac{(1+z)\log(1+z)+(1-z)\log(1-z)}{2\sqrt z}
=\sum_{k\ge0}\frac{\exp(-(2k+3/2)|v|)}{(2k+3/2)^2-1/4}.
```

The second evaluation bounds the omitted positive tail geometrically after 128 terms. Every argument has absolute value larger than $`0.34`$. Both calculations enclose $`-1.990\lt U(s)\lt-1.989`$ and oddness gives the opposite value at $`-s`$. The width $`2s=0.00546`$ is strictly less than $`\log(2.011/2)`$. Leaving slack at both ends, translate this pair by uniform approximate periods of $`U`$. Both actual prime-only moves retain their strict margins since their errors tend to zero. The same bounded-gap argument from the preceding section yields the asserted density. This is recurrence of the *pair* in one window. The computation at the origin is an evaluation of the limiting profile. [Earlier even explicit formula](https://github.com/leecrossley/theta-race-distant-interval/blob/527fffe7b7808457b47a292144fea96b556158a8/PROOF.md#a-positive-bound-for-all-omitted-ordinates) · [Both numerical evaluations](check_close_pair.py)

## The sharp amplitude at factor four

Keep the score, the five finite exclusions and the strict start convention above. Let $`a=\log2`$ and let $`B_C(H)`$ be the supremum of the common normalized magnitudes of a gain and a loss on two eligible doublings starting strictly after $`H`$ and finishing by $`CH`$. Equivalently take the smaller of the positive parts of the supremum of $`f`$ and of $`-f`$ over the eligible log interval, where $`f`$ is the actual normalized move. A strict threshold below $`B_C(H)`$ has two actual witnesses. Put $`W=\log(C/2)`$, the width available to their starts in log coordinates. We consider $`2\lt C\le4`$.

**Theorem.** If the eight odd factors have their strip zeros on the critical line, then

```math
\limsup_{H\to\infty}B_4(H)=
A=\frac{V(0)-V(a)}{a}
\in(269.3640010193614394,\ 269.3640010193614395).
```

For every threshold $`0\lt t\lt A`$, strict successes at factor four have positive lower logarithmic density. If any of the eight has an off-line zero, the upper limit is infinite for every $`C\gt2`$ and every fixed positive threshold has successes of positive lower logarithmic density. In particular $`269.364\log2`$ in both directions works before $`3.999999997H`$ on sets of positive lower logarithmic density **without any assumption**. Recurrence before $`3.999999989H`$ at that same threshold is equivalent to an off-line zero. This window comparison does not assert an outcome at the exact threshold $`A`$ itself.

### The level profile decreases through the first dyadic interval

Work first in the critical-line case. The absolutely convergent even profile is $`V(t)=\sum_\nu G(\nu)e^{i\nu t}`$, over all odd-family ordinates with multiplicity. The coefficients $`G(\nu)`$ are real, even and positive. The limiting move is $`U(t)=(V(t+a/2)-V(t-a/2))/a`$. In particular $`f(t)=U(t)+o(1)`$ uniformly on the late tail. The base-field prime-square term tends to its constant by the prime ideal theorem, and this constant cancels in the difference. Only the eight odd factors enter the critical-line premise. [Exact profile](#a-short-reversal-detects-an-off-line-zero) · [Prime-only arithmetic remainder](https://github.com/leecrossley/theta-host-split-inert-race/blob/compact-filter-oscillation-v1.1.0/results/compact-filter-oscillation/PROOF.md#what-survives-on-sparse-grids)

For $`s=\sqrt2`$, $`c=2(s-1)^2`$ and $`0\le x\le a`$, the elementary explicit formula is

```math
V(x)=K\,\frac{4\sinh((a-x)/2)}{sc}
-\frac8c\left[s\big(T(x+a)+T(a-x)\big)-3T(x)\right],
\qquad K=\sum_\chi\log(q_\chi/\pi)-8(\gamma_{\rm E}+a).
```

Here $`q_\chi`$ runs over the eight primitive conductors, and $`K\gt47`$. In the translated tent no prime power above 4 occurs. At 2 the two nonvanishing characters cancel, at 3 the quartic twist sum vanishes, and at 4 the tent vanishes on the possible boundary $`x=a`$. The same holds at the reflected endpoint. Thus the formula extends to the entire closed interval, including its cusps. For $`x\ge0`$ the positive gamma kernel is

```math
T(x)=\sum_{k\ge0}\frac{e^{-(2k+3/2)x}}{(2k+1)(2k+2)},\qquad T(0)=\log2 .
```

On $`(0,\infty)`$, $`T'\lt0`$ and $`T''\gt0`$ by termwise differentiation on compact subintervals. For $`0\lt x\lt a`$, the derivative of the bracket in the explicit formula is $`s(T'(x+a)-T'(a-x))-3T'(x)\gt0`$. The tent is strictly decreasing and $`K\gt0`$. Hence $`V`$ is strictly decreasing on $`[0,a]`$. This last conclusion uses the absence of a nonzero prime contribution throughout this particular range. [Gamma series and tent formula](https://github.com/leecrossley/theta-host-split-inert-race/blob/9fdf8c9dea230464fc8ab635d1bdf1659c3539fd/results/asymptotic-waiting-factors/PROOF.md#both-signs-in-a-window-ending-by-2011h) · [Underlying even explicit formula](https://github.com/leecrossley/theta-race-distant-interval/blob/527fffe7b7808457b47a292144fea96b556158a8/PROOF.md#a-positive-bound-for-all-omitted-ordinates)

### An upper bound for every choice of two starts

Consider any two log centres at distance $`d\le W`$. Half the absolute difference between their limiting moves is at most

```math
\frac{2}{a}\sum_\nu G(\nu)
  |\sin(a\nu/2)\sin(d\nu/2)|
\ \le\ \frac{\sqrt{(V(0)-V(a))(V(0)-V(d))}}{a}.
```

The first bound follows by taking absolute values of the two Fourier differences, with both signs of every ordinate included. Cauchy–Schwarz gives the second since $`2\sum_\nu G(\nu)\sin^2(h\nu/2)=V(0)-V(h)`$. The smaller of two opposite magnitudes is at most half their separation in value. Strict decrease of $`V`$ bounds the right side by its value at $`d=W`$. The uniform prime-only remainder tends to zero independently of the two chosen starts. This proves a whole-window upper bound and permits either order and any overlap of the doublings. [Formal scalar deductions](ExtremalWindows.lean)

For a lower bound use the symmetric log centres $`-W/2`$ and $`W/2`$ in the limiting profile. Evenness of $`V`$ gives values $`L(W)`$ and $`-L(W)`$, where

```math
L(W)=\frac{V((a-W)/2)-V((a+W)/2)}{a},\qquad
P(W)=\frac{\sqrt{(V(0)-V(a))(V(0)-V(W))}}{a}.
```

For every $`t\lt L(W)`$, continuity lets us shorten the distance between the centres slightly, keeping both strict margins and leaving space at both endpoints of a window of width $`W`$. The uniformly almost-periodic zero series has a relatively dense set of translations preserving both margins. Uniform continuity gives log intervals of placement of fixed positive length at bounded gaps. The actual prime-only error tends to zero there. Such a set has positive lower logarithmic density in $`H`$. Consequently

```math
L(W)\le\limsup_{H\to\infty} B_C(H)\le P(W).
```

At $`W=a`$ both expressions equal $`A`$. The symmetric profile witnesses are an up move followed by a down move. To enforce the strict start at the exact completion factor four we first shorten their separation as above. Near equality also forces the ratio of the later start to the earlier start to tend to two, since any fixed shorter separation has strictly smaller $`P`$. Equality of zero frequencies or rational relations among ordinates causes no difficulty in this argument.

### Two deadlines eight billionths apart

Put $`t=269.364`$, $`C_+=3.999999997`$ and $`C_-=3.999999989`$. The two outward-rounded computations give

```math
L(\log(C_+/2))\gt t+0.00000015,\qquad
P(\log(C_-/2))\lt t-0.00000008.
```

The second computation integrates the positive series using $`1/((2k+1)(2k+2))=\int_0^1 u^{2k}(1-u)\,du`$. After summation and the change of variable $`v=1-u`$, the integrand for $`T(x)`$ is

```math
e^{-3x/2}\frac{v}{(1-z+zv)(1+z-zv)},\qquad z=e^{-x}.
```

At $`x=0`$ it reduces to $`1/(2-v)`$. At positive $`x`$ the verifier integrates over dyadic subintervals with complex ball quadrature. The first computation uses the closed logarithmic expression for $`T`$, including the exact limiting value at zero. Both methods enclose the endpoint and the two strict comparisons. The coefficient positivity, monotonicity and the all-start upper bound are analytic deductions. [Numerical replay](fourfold_amplitude.py) · [Transcript](evidence/check.txt)

If an odd factor has an off-line zero, the earlier positive-density theorem applies at each of these factors and every fixed threshold. It also gives the infinite upper limit by applying it to arbitrarily large fixed thresholds. Combining the two cases proves the unconditional success at $`C_+`$ and the equivalence at $`C_-`$. In the critical-line case the infimum of factors with arbitrarily late successes at amplitude $`t`$ belongs to $`[C_-,C_+]`$, an interval of width $`8\cdot10^{-9}`$. In the off-line case that infimum is 2. To obtain a fixed deadline after **every** $`H\ge1`$ for the amplitude $`t`$, use the bounded gaps between the intervals of successful log heights in either case. Enlarge the factor once more so that all bounded initial heights can use one fixed later successful pair. This argument gives no numerical bound for that all-height factor. [Off-line and bounded-gap conclusions](#a-short-reversal-detects-an-off-line-zero)
