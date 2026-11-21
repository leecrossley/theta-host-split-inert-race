# The four possible weights in a two-field prime-race inverse

Three disjoint doublings can make the 547/727 and 571/599 prime races balance exactly. Fix the integer start at **13,471** and keep the other two starts on the prime bands containing **2,457.452762** and **61,460.713124**. Given the four target measurements, the weight at the integer start solves a quartic.

**Exactly one of four real roots works.** This holds throughout the closed box of radius 0.0001 about the balanced target (both midpoint levels one, both moves zero). The other three roots violate the signs required by the original prime sums. The conclusion covers the whole ordered pair of prime bands. For almost every target in this box the four roots have Galois group $`S_4`$ over the target coordinates and the real algebraic numbers with the finite prime logarithms adjoined. In particular, field operations and iterated square roots cannot recover even the admissible weight. The proof uses finite prime counts and has no zero-location hypothesis. [Proof](PROOF.md#the-four-conjugates-of-the-integer-start-weight) · [Certificate](evidence/inverse-monodromy.txt)

**An elliptic curve records each way to choose a pair.** On a specified one-parameter slice, forgetting which of two complementary root pairs was selected gives a degree-two isogeny. The two elliptic curves have $`j`$-invariants approximately 2,159.399671 and 1,372.295119. The first has two real components and the second has one. Two independent calculations enclose both invariants to 18 decimal places. [Equations of the two curves](PROOF.md#the-two-elliptic-resolvents) · [Branch calculation](inverse_branch_recheck.py) · [Prior work](PRIOR-WORK.md#the-quartic-inverse-and-its-resolvents)

## Three disjoint doublings for two theta prime races

For the 547/727 and 571/599 theta prime races at width two, **three common starts reproduce any four targets in the following box exactly**. Take a positive weighted average of the midpoint levels and the normalized dyadic moves at the three starts. Use exactly the same weights and heights for both fields.

| Quantity | 547/727 | 571/599 |
| --- | --- | --- |
| Midpoint level | $`[-9,11]`$ | $`[-7,9]`$ |
| Move divided by $`\log2`$ | $`[-3.5,3.5]`$ | $`[-2.2,2.2]`$ |

**Finite theorem.** One start is an integer chosen from the two target levels. All three starts exceed 1,155 and all three doublings finish below 701,816. The doublings are pairwise disjoint and the last one finishes before 110 times the first start. Each weight exceeds $`1/400`$. These assertions follow from direct prime-ideal counts with no zero-location assumption. Almost every target in the box needs all three starts, even if observations at arbitrary positive heights are allowed. [Proof](PROOF.md#three-common-starts-for-the-four-measurements) · [23 timing configurations](evidence/three-start-charts.tsv) · [Lookup by target levels](evidence/three-starts.txt)

**Beyond every height.** The same configurations recur on the log-height axis with bounded gaps and positive lower uniform density under the critical-line condition for the sixteen primitive Dirichlet factors of the two fields. The published primitive-character theorem supplies this condition. Every successful window realizes the whole target box with the same bound of 110 measured from the first of the chosen starts. The wait from a prescribed $`H`$ to a working window has a finite factor whose value is open. The two earlier level–move ellipses lie entirely inside this box. [Transfer to actual primes](PROOF.md#the-limiting-curves-and-late-physical-starts) · [Primitive theorem and independent checks](https://math.mit.edu/~drew/grh/)

**An exact example.** At starts approximately 2,457.452761630442, exactly 13,471, and approximately 61,460.713123523080, one set of three positive weights gives both levels exactly one and both moves exactly zero. A quartic isolates the solution to 18 decimal places. An independent interval contraction verifies the original four equations. Every target within 0.0001 of this example has a unique solution in a specified parameter box with the same integer start. For almost every such target the local inverse has algebraic degree exactly four over the observations and the field of prime logarithms with real algebraic coefficients. [Balanced example and degree proof](PROOF.md#an-exact-balanced-example-and-the-degree-of-its-local-inverse) · [18-place brackets](evidence/three-inverse.txt)

The two calculations check complete arcs and all closed target cells. One uses the jumps of the tent at prime norms. The other factors defining polynomials, evaluates all primitive Conrey characters, and bounds the error of each secant. The original 17 local Lean deductions include the fixed-point algebra for the three-start argument. The analytic and topological theorems are written proofs. [Numerical replay](three_starts.py) · [Independent verifier](three_recheck.py) · [Lean](JointOscillation.lean) · [Prior work](PRIOR-WORK.md#three-common-starts)

## The preceding fixed-window bounds

**After every real $`H\ge1`$, the same finite completion factor $`C`$ works for every target pair in $`[-5.85,5.85]^2`$.** At most two actual common doublings reproduce the target pair of normalized moves exactly, using the same nonnegative weights in both fields. Each start is strictly after $`H`$, and each doubling finishes before $`CH`$. The fields are 547/727 and 571/599, both at width two. The proof establishes existence of $`C`$. Its numerical value remains open. The proof is unconditional. [Fixed-window theorem](PROOF.md#one-fixed-window-for-all-directions)

On **one actual doubling** the first race can outgain the second by more than $`11.7\log2`$. Another doubling gives the reverse advantage before the same deadline. Almost every pair in the square needs two common starts. The four-coordinate result reproduces both midpoint levels and both moves with at most four common starts. Its full closed target is the product of the ellipses below. Almost every point in that four-dimensional target needs at least three starts. The three-start construction above settles the dyadic box in a finite interval and on arbitrarily late windows with its stated primitive Dirichlet input. The zero-location-independent four-start theorem continues to apply at arbitrary fixed widths. [Common weights](PROOF.md#resolving-collisions-at-a-common-ordinate) [Lower bounds on the number of readings](PROOF.md#how-many-observations-are-necessary)

| Field | Width-two semiaxes for level and move | Semiaxes at any fixed width greater than one |
|---|---|---|
| 547/727 | 9.9375 and 3.4686 | 6.6921 and 1.9735 |
| 571/599 | 7.9731 and 2.1866 | 6.1027 and 1.5894 |

Each level ellipse is centred at level one and move zero. The widths may be chosen independently, with a fixed factor for each choice. For all sufficiently large $`H`$ the bounded targets also fit before $`H(\log H)^{8+\varepsilon}`$ for every fixed $`\varepsilon\gt0`$. The proof gives no numerical onset. Neither result assumes global Dirichlet zero location. The proof checks even directions that change with $`H`$ and handles several growth rates at the same zero ordinate. [Logarithmic exponent](PROOF.md#a-numerical-exponent-independent-of-the-zeros) · [Reproducible checks](evidence/check.txt)

## Earlier explicit and finite-height bounds

Compare the change in the 547/727 race with the change in the 571/599 race **on the same doubling**. After every real $`1\le H\le e^{2000}`$, there is a doubling where the first field outgains the second by more than $`8\log2`$, and another where the second outgains the first by more than $`8\log2`$. Both starts exceed $`H`$. Both doublings finish before $`262144H`$. The result is unconditional. [Finite-height proof](PROOF.md#the-finite-common-windows) · [Verification transcript](evidence/check.txt)

For every $`H\ge1`$ the corresponding advantage exceeds $`7\log2`$ in each direction before $`e^{14004}H^{17001/16999}`$. For each fixed $`\delta\gt0`$ it exceeds $`11.7\log2`$ in each direction before $`H^{1+\delta}`$ for all sufficiently large $`H`$. The first two deadlines are explicit. The last one gives no numerical onset. Both statements are unconditional. [All-height proof](PROOF.md#a-common-test-at-every-height) · [Power intervals](PROOF.md#larger-contrasts-on-power-intervals)

There is also an exact averaging result. Write $`D_A,D_B`$ for the two changes divided by $`\log2`$ and $`L_A,L_B`$ for the two midpoint scores.

| Starts $`H`$ | What a common weighted average can match | Starts needed | Completion deadline |
|---|---|---|---|
| $`1\le H\le e^{2000}`$ | Any pair in $`[-4,4]^2`$ for the moves | At most 2 | $`262144H`$ |
| $`1\le H\le e^{2000}`$ | Any four coordinates with $`\lvert L_i-1\rvert\le2`$ and $`\lvert D_i\rvert\le4`$ | At most 4 | $`2^{80}H`$ |
| Every $`H\ge1`$ | Any pair in $`[-3.5,3.5]^2`$ for the moves | At most 2 | $`e^{14004}H^{17001/16999}`$ |
| Every $`H\ge1`$ | In each field, any point in the ellipse centred at $`(1,0)`$ with semiaxes 1.17 and 2.45 | At most 4 | $`e^{14004}H^{17001/16999}`$ |

The starts and nonnegative weights are the same for both fields. A projected advantage is attained on one actual doubling. Prescribing two or four coordinates allows several doublings with weights adding to one. The explicit deadlines in this table use width two and the five stated prime exclusions for each host. In the fixed-factor ellipse theorem the two widths may be chosen independently. The original arbitrary-width, arbitrarily late theorem follows below. [Common starts and weights](PROOF.md#why-these-are-actual-prime-observations) [Original release](https://github.com/leecrossley/theta-host-split-inert-race/releases/tag/joint-levels-and-moves-v1.0.0)

Run `bash verify-release.sh` in the result directory after installing the pinned arithmetic package. The release includes 3,218 genuine positive zero brackets, two evaluations of 1,555 common averages that use different formulas for the kernel and the trigonometric sums but share the zero table, costs and geometry, eight early prime-ideal counts and thirteen local Lean deductions in the earlier proofs. The written analytic and geometric steps identify how the certificates imply the result. [File hashes](SHA256SUMS) · [Prior work](PRIOR-WORK.md)

## The earlier unrestricted-width result

Either the 547/727 race or the 571/599 race can move further on the same doubling, arbitrarily far out. At width two the difference between their normalized moves exceeds 5.6552 in each direction. The difference between their levels exceeds 17.9106 in each direction. These bounds are unconditional. They use nine finite zero brackets and real-axis checks for sixteen Dirichlet characters. Landau’s theorem lets the argument use these finitely many zeros without assuming GRH or independence of the remaining zeros. [Proof](PROOF.md) · [Verification output](evidence/check.txt)

For each pair $`i=(p,q)`$ put $`K=\mathbf Q(\zeta_{16})^+`$, $`F_i=K(\sqrt{pq})`$ and $`E_i=K(\sqrt{-p},\sqrt{-q})`$. A prime ideal of $`F_i`$ contributes positively for inertia and negatively for splitting in $`E_i`$. Write $`S_{i,b}(y)`$ for the compactly weighted count divided by $`2\sqrt b(\sqrt b-1)^2\sqrt y`$. The two races share the same starting height $`y`$. Their widths $`b_i\gt 1`$ may be chosen independently. [547/727 trace identity](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-short-shifts/RaceTrace.lean) · [571/599 trace identity](https://github.com/leecrossley/theta-host-split-inert-race/blob/8fca7cad37185e3d7951f0dffff5e2442d1dfa24/HostAdmissible.lean)

```math
L_i(y)=S_{i,b_i}(\sqrt2\,y),\qquad
D_i(y)=\frac{S_{i,b_i}(2y)-S_{i,b_i}(y)}{\log 2}.
```

There is a stronger joint result. For each race draw an ellipse centred at $`(1,0)`$ in its level–move plane, with the following horizontal and vertical semiaxes.

| Field | Level, width 2 | Move, width 2 | Level, any fixed width | Move, any fixed width |
|---|---:|---:|---:|---:|
| 547/727 | 9.9375 | 3.4686 | 6.6921 | 1.9735 |
| 571/599 | 7.9731 | 2.1866 | 6.1027 | 1.5894 |

Choose a point on or inside each ellipse. Beyond any specified height, a weighted average of **at most four common starts** matches all four coordinates exactly. The weights are nonnegative, sum to one and are the same for both races. For the two moves alone, any point of the corresponding rectangle takes at most two common starts. The two individual difference bounds 5.6552 and 17.9106 above each hold at single heights. [Convex-hull argument](PROOF.md#common-starts-and-the-number-of-readings)

Finite counts give examples where both normalized moves have magnitude greater than 10. In field order 547/727 then 571/599, the signs at starts 17, 97, 162 and 658 are respectively $`(-,+)`$, $`(-,-)`$, $`(+,-)`$ and $`(+,+)`$. These counts omit the primes over 2 and the two conductors in each race. [Direct counts](finite_examples.py)

## Verification

Use Python 3.10 or later with the pinned python-flint package. On Linux x86-64 with Git, curl, zstd and a C toolchain available, the script downloads the pinned Lean toolchain and Mathlib. With an existing checkout, pass its path as the sole argument and put Lean 4.33.1 on `PATH`.

```sh
python3 -m pip install -r requirements.txt
bash verify-release.sh
```

The replay checks the real segments, the nine sign changes, the mixed frequencies by two independent evaluations, the ellipse bounds, the finite examples, the deadline spectrum, the early prime-ideal counts, the all-height common test, both methods for the common averages, the fixed-window checks and the six original Lean deductions, four earlier deadline deductions and three fixed-window deductions. Both evaluations use Arb ball arithmetic. The Mellin, positive-weight and convex-geometric arguments are written proofs. [Checksums](SHA256SUMS) · [Prior work](PRIOR-WORK.md)
