# Identifying the two omitted prime conductors

## The observation model

Fix $`K=\mathbf Q(\zeta_{16})^+`$. A candidate is an unordered pair of distinct primes $`u,v\equiv3\bmod4`$. Set $`F=K(\sqrt{uv})`$ and $`E=K(\sqrt{-u},\sqrt{-v})`$. For a fixed known finite set $`\mathcal B`$ of odd primes, require $`u,v\notin\mathcal B`$. Omit all places over $`2,u,v`$ and $`\mathcal B`$. At each other prime ideal of $`F`$ let $`\epsilon=1`$ for split and $`-1`$ for inert in $`E/F`$. The 547/727 race uses $`\mathcal B=\{343057911569,343070636977\}`$. Its sign is the product of the normalized first-step theta trace jumps. The comparison family consists of these algebraically defined splitting races. Only the target needs the specific admissible theta datum. [Concrete theta trace](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-short-shifts/RaceTrace.lean#L83) · [Field and character frame](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-short-shifts/PROOF.md)

Write $`A(y)=B\sqrt y\,S(y)`$ for the raw score, where $`B=2\sqrt2(\sqrt2-1)^2`$ is common to all candidates. The score at a positive real height is

```math
A(y)=\sum_{n\ge1}c_n W(n/y),\qquad c_n=-\sum_{Nw=n}\epsilon_w\log n,\qquad
W(t)=\max(0,\min(2-2/t,4/t-1)).
```

Two candidates have the same observed data if their scores agree at every specified height. Identification means that only the target unordered pair gives those data. The conventions deliberately leave all masses of a conductor prime out of its own race. The theorem uses this fact.

## The exact observation bands

For an odd prime $`p`$ let $`d(p)`$ be its order in $`(\mathbf Z/16\mathbf Z)^\times/\{\pm1\}`$. Put $`J(n)=(n/4,n)`$. For each target conductor $`p`$, form the following list of possible unramified norms,

```math
\mathcal N(p)=
\begin{cases}
\{p^{d(p)}\},&d(p)=2\text{ or }4,\\ \{p,p^2\},&d(p)=1.
\end{cases}
```

**Theorem.** A finite schedule of positive rational heights identifies the target pair exactly if it meets every open interval $`J(n)`$ for both conductors and all $`n\in\mathcal N(p)`$. If any required interval is missed, infinitely many other pairs give the same observations. This last assertion holds for every bounded set of real heights missing such an interval, including a continuum. The theorem is unconditional: it uses Dirichlet's theorem as stated in Mathlib and unique factorisation, and no Riemann hypothesis.

To prove the local alternatives, at a good rational prime $`t`$ let $`s=(t/u)`$ and $`z=(t/v)`$. The group over the rationals is the product of the quartic cyclotomic and two quadratic components, which have pairwise coprime conductors. Complex conjugation has components $`(1,-1,-1)`$. In the real subfield the residue degree $`f`$ is the least positive multiple of $`d(t)`$ with $`(sz)^f=1`$. The $`8/f`$ prime ideals above $`t`$ all have norm $`t^f`$ and sign $`s^f`$. Their combined mass is $`-8s^f\log t`$. Thus $`d=2,4`$ forces $`f=d`$ and mass $`-8\log t`$ irrespective of the two quadratic signs. At $`d=1`$ the norm is $`t`$ if $`s=z`$, with mass $`-8s\log t`$, and $`t^2`$ otherwise, with mass $`-8\log t`$. A good prime has precisely one such norm. Its weight is positive exactly at the heights in $`J(t^f)`$. [Eight-character field computation](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-finite-euler-readout/PROOF.md#which-euler-factor-is-recovered) · [All residue cases](horizon.py)

At rational $`y`$, $`A(y)/8`$ is a finite rational linear combination of logarithms of distinct rational primes. Its coefficients are unique. Indeed, in any equality clear the finitely many denominators, exponentiate and move the negative powers across. Unique factorization of the two equal positive integers gives equality of each prime exponent. Consequently any candidate in which a target conductor $`p`$ is good must disagree at a rational height inside $`J(p^f)`$. The target has zero coefficient of $`\log p`$, and the candidate has a nonzero one there. If every required interval is hit, each target conductor is forced to be omitted. It lies outside the fixed auxiliary set and it is odd, so both conductors are forced. [Prime-product injectivity in Lean](ConductorHorizon.lean) · [Rational sampling checks](sampling.py)

For necessity suppose $`J(p^f)`$ was missed. Keep the other conductor $`r`$ fixed and let $`Y`$ be any finite upper bound for the observations. Choose a prime replacement $`q`$ larger than $`4Y`$, both target primes and all auxiliary primes. Require $`q\equiv p\bmod16`$ and $`q\equiv p\bmod t`$ for every odd prime $`t\le4Y`$ with $`d(t)=1`$, except $`p,r`$ and the auxiliary primes. These residues are units. At a good prime with $`d(t)>1`$ the mass is independent of the quadratic signs. Both races omit $`r`$ and the auxiliary primes. If $`d(p)=1`$, add a residue condition modulo $`p`$. To select $`f=1`$ demand $`(p/q)=(p/r)`$. To select $`f=2`$ demand the opposite. Both signs are available because $`q\equiv p\equiv3\bmod4`$ and quadratic reciprocity gives $`(p/q)=-(q/p)`$. For $`d(p)\gt 1`$ the unique possible degree is forced already.

The Chinese remainder theorem now gives a reduced residue class and Dirichlet's theorem supplies infinitely many primes in it exceeding every prescribed bound. At each rational prime $`t\le4Y`$ other than $`p`$, reciprocity preserves $`(t/q)=(t/p)`$ wherever the symbols are needed. Contributions and exclusions therefore match. The prime $`q`$ itself is beyond the entire support. At $`p`$ the replacement has a contribution only on the missed open interval. This proves equality at every prescribed height, including interval endpoints where the weight is zero. The replacement fields are distinct since the ramified rational prime changes. [Dirichlet's theorem in Mathlib](https://github.com/leanprover-community/mathlib4/blob/0df444a360eaa60ab8c11dca51a86af692955474/Mathlib/NumberTheory/LSeries/PrimesInAP.lean#L463) · [Explicit finite examples](evidence/twins.tsv)

## The exact endpoint and the number of readings

Complete observations on $`0\lt y\le Y`$ determine a target pair $`u,v`$ precisely when

```math
Y\gt Y_*=\tfrac14\max\{u^{\max(2,d(u))},v^{\max(2,d(v))}\}.
```

If $`Y\le Y_*`$ the open interval for the largest possible norm has been missed. If $`Y\gt Y_*`$ a rational point in each requisite band may be taken at or below $`Y`$. For the target, $`d(547)=4`$, $`d(727)=2`$, $`547^4=89526025681`$ and $`727^2=528529`$. The two bands are disjoint. Hence every single rational reading admits infinitely many alternatives, and any pair of rational readings hitting the two bands suffices. In fact the missed-band argument gives the same single-reading obstruction even for an arbitrary real height. One exact reading per band is a statement about complete real values. The next construction provides a practical error tolerance.

## Six measured scores with an error bound

For either odd diagnostic norm $`n=p^d`$ take raw readings at $`(n-\theta)/4,n/4,(n+\theta)/4`$, where $`0\lt \theta\le1`$. In the tick variable $`k=4y`$, the unit kernel at integer mass $`m`$ is
```math
[(k-m)_+-\tfrac32(k-2m)_++\tfrac12(k-4m)_+]/m,\qquad x_+=\max(0,x).
```

On this symmetric interval a hinge with integer knot has second difference $`\theta`$ at the central knot and zero otherwise. At odd $`n`$ only the first knot of mass $`n`$ can coincide with the centre. Thus

```math
\frac n\theta\{A((n+\theta)/4)-2A(n/4)+A((n-\theta)/4)\}=c_n .
```

For the two target tests the alternatives are exactly $`0`$ and $`-8\log p`$. Three errors bounded in absolute value by $`\delta`$ induce error at most $`4n\delta/\theta`$ in the recovered mass. Hence $`\delta\lt \theta\log p/n`$ suffices. For any $`Y\gt Y_*`$ choose $`\theta\gt 0`$ small enough that the last reading fits. The actual computation uses $`\theta=1`$. It reaches starting height $`22381506420.5`$. [Finite inverse](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-finite-euler-readout/PROOF.md) · [Fractional-step and rounding deductions](ConductorHorizon.lean)

For each of the six ticks the archive records integer endpoints for $`10^{12} A(k/4)`$. The widest enclosure spans two units. Independent rounding of each raw value to this accuracy is sufficient. The three-point linear combinations from the archived intervals enclose zero, and their unnormalised masses each have absolute value below one. An unramified alternative has magnitude $`8\log p\gt 8`$. A separate calculation recovers all six enclosures from the archived sieve bands using rational arithmetic. At the smaller norm a further check factors the defining quartic and quadratics over the individual residue fields and sums the three prime-ideal scores. [Six raw enclosures](evidence/raw-scores.tsv) [Direct polynomial check](independent.py)

The sieve enumerates good degree-one contributions in the progressions $`1,15\bmod16`$ and enumerates all higher-degree contributions from primes up to the square root of the cap. It accumulates signed $`\log t`$ and $`(\log t)/t^f`$ in bands separated by the kernel knots, in units $`2^{-80}`$ and $`2^{-112}`$. For each positive prime write $`t=2^k(1+j/2^{16})(1+z)`$, with $`0\le z\lt 2^{-16}`$. The first two logarithms have certified lower integer table entries. For the last logarithm use four alternating Taylor terms. For powers of $`z`$, the recurrence computing $`2^{80}z^i`$ rounds downward at each multiplication. Its absolute error stays below two integer units since $`z\lt 1/2`$. After division by $`i`$ and truncation the error is less than three units per term, and for $`i=1`$ it is less than one. The fifth Taylor term is below one unit. Padding the correction by 12 below and 13 above therefore suffices. The table contributes at most $`k+1`$ additional units, yielding the implemented upper padding $`k+14`$. Rounded reciprocal bounds are obtained by integer division. At the enforced cap of 90 billion, even the sum over every integer of $`\log n`$ fits in the chosen 128-bit units. For the reciprocals use $`\sum_{n\le x}(\log n)/n\le(1+\log x)\log x`$. The comparison with separate primality tests also checks direct Arb logarithms by band through 10,000,000, with band edges on either side of the first two sieve segment boundaries. [Segmented sieve](sieve.cpp) [Log table and sample calculation](measure.py) [Integer bands](evidence/bands.tsv)

## What an arbitrarily high tail determines

Fix two distinct pairs in this family. Choose a reduced residue class modulo $`m=16\prod q`$, the product ranging over the distinct conductors of the two pairs. Require residue one modulo 16 and quadratic signs $`+1`$ at both members of the first pair. Choose $`-1`$ at one conductor belonging exclusively to the second pair and $`+1`$ at its other conductor. These requirements are compatible by the Chinese remainder theorem. Every good prime $`t`$ in this class has degree one with nonzero mass in the first race and degree two in the second.

For every sufficiently large $`y`$, the progression contains a prime $`y\lt t\le2y`$. For an explicit sufficient range take $`y`$ greater than all excluded primes, greater than 4 and at least
```math
\max\!\left\{x_0(m),\,\exp(\varphi(m)/40)\right\},\qquad
x_0(m)=\begin{cases}
8\cdot10^9,&m\le10^5,\\ \exp(0.03\sqrt m\,(\log m)^3),&m\gt 10^5.
\end{cases}
```

Theorem 1.2 of Bennett, Martin, O'Bryant and Rechnitzer bounds both progression theta errors by $`x/(160\log x)`$ once $`x\ge x_0(m)`$, unconditionally. Thus $`\vartheta(2y;m,a)-\vartheta(y;m,a)\ge y/\varphi(m)-3y/(160\log y)\gt 0`$. At such a prime $`t`$, the first race assigns a positive kernel weight to its degree-one term and $`t^2\gt 4y`$ places the other term outside the support. Prime-log independence gives unequal scores at **every rational** $`y`$ in the stated range. Hence exact scores on any unbounded sequence of rational heights distinguish any two fixed candidates. Their cutoff can depend on both candidates. For any bounded set of real heights missing either diagnostic band, the missed-band construction gives infinitely many pairs agreeing with the target throughout the observations. In particular, a bounded window starting at or above $`727^2=528529`$ misses $`J(727^2)`$. [Bennett, Martin, O'Bryant and Rechnitzer, Theorem 1.2](https://personal.math.ubc.ca/~bennett/BeMaObRe-IJM-2018.pdf#page=5)

## Reproducible examples

Three explicit prime replacements demonstrate each residue-degree case. For targets $`(3,7),(23,3),(31,3)`$ the example pairs first differ in their integer-norm masses at 81, 529 and 961 respectively. The program proves the replacement conductors prime, checks the reduced residue conditions and exhausts all smaller prime norms. In the first example the two diagnostic bands overlap. Its one-reading test at $`y=25`$ is also checked against 27 different pairs using exact rational coefficients. The complete results above are analytic and algebraic proofs for arbitrarily large conductors. Seven named local deductions, including the Dirichlet progression and prime-product steps, are replayed in Lean. [Alternative fields](horizon.py) · [Schedule tests](sampling.py) · [Replay transcript](evidence/check.txt)
