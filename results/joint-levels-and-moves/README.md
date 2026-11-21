# A deadline for either theta prime race to outgain the other

Compare the change in the 547/727 race with the change in the 571/599 race **on the same doubling**. After every real $`1\le H\le e^{2000}`$, there is a doubling where the first field outgains the second by more than $`8\log2`$, and another where the second outgains the first by more than $`8\log2`$. Both starts exceed $`H`$. Both doublings finish before $`262144H`$. The result is unconditional. [Finite-height proof](PROOF.md#the-finite-common-windows) · [Verification transcript](evidence/check.txt)

For every $`H\ge1`$ the corresponding advantage exceeds $`7\log2`$ in each direction before $`e^{14004}H^{17001/16999}`$. For each fixed $`\delta\gt0`$ it exceeds $`11.7\log2`$ in each direction before $`H^{1+\delta}`$ for all sufficiently large $`H`$. The first two deadlines are explicit. The last one gives no numerical onset. [All-height proof](PROOF.md#a-common-test-at-every-height) · [Power intervals](PROOF.md#larger-contrasts-on-power-intervals)

There is also an exact averaging result. Write $`D_A,D_B`$ for the two changes divided by $`\log2`$ and $`L_A,L_B`$ for the two midpoint scores.

| Starts $`H`$ | What a common weighted average can match | Starts needed | Completion deadline |
|---|---|---|---|
| $`1\le H\le e^{2000}`$ | Any pair in $`[-4,4]^2`$ for the moves | At most 2 | $`262144H`$ |
| $`1\le H\le e^{2000}`$ | Any four coordinates with $`\lvert L_i-1\rvert\le2`$ and $`\lvert D_i\rvert\le4`$ | At most 4 | $`2^{80}H`$ |
| Every $`H\ge1`$ | Any pair in $`[-3.5,3.5]^2`$ for the moves | At most 2 | $`e^{14004}H^{17001/16999}`$ |
| Every $`H\ge1`$ | In each field, any point in the ellipse centred at $`(1,0)`$ with semiaxes 1.17 and 2.45 | At most 4 | $`e^{14004}H^{17001/16999}`$ |

The starts and nonnegative weights are the same for both fields. A projected advantage is attained on one actual doubling. Prescribing two or four coordinates allows several doublings with weights adding to one. The newer deadlines fix both window widths at two and use all five rational omissions specified for each host. The original arbitrary-width, arbitrarily late theorem follows below. [Common starts and weights](PROOF.md#why-these-are-actual-prime-observations) · [Original release](https://github.com/leecrossley/theta-host-split-inert-race/releases/tag/joint-levels-and-moves-v1.0.0)

Run `bash verify-release.sh` in the result directory after installing the pinned arithmetic package. The release includes 3,218 genuine positive zero brackets, two independent evaluations of 1,555 common averages, eight early prime-ideal counts and ten local Lean deductions. The written analytic and geometric steps identify how the certificates imply the result. [File hashes](SHA256SUMS) · [Prior work](PRIOR-WORK.md)

## The earlier unrestricted-width result

Either the 547/727 race or the 571/599 race can move further on the same doubling, arbitrarily far out. At width two the difference between their normalized moves exceeds 5.6552 in each direction. The difference between their levels exceeds 17.9106 in each direction. These bounds are unconditional. They use nine finite zero brackets and real-axis checks for sixteen Dirichlet characters. Landau’s theorem lets the argument use these finitely many zeros without assuming GRH or independence of the remaining zeros. [Proof](PROOF.md) · [Verification output](evidence/check.txt)

For each pair $`i=(p,q)`$ put $`K=\mathbf Q(\zeta_{16})^+`$, $`F_i=K(\sqrt{pq})`$ and $`E_i=K(\sqrt{-p},\sqrt{-q})`$. A prime ideal of $`F_i`$ contributes positively for inertia and negatively for splitting in $`E_i`$. Write $`S_{i,b}(y)`$ for the compactly weighted count divided by $`2\sqrt b(\sqrt b-1)^2\sqrt y`$. The two races share the same starting height $`y`$. Their widths $`b_i\gt 1`$ may be chosen independently. [547/727 trace identity](https://github.com/leecrossley/theta-race-short-shifts/blob/483c31f2df2132ed2165620ea55e20e969228bfb/RaceTrace.lean) · [571/599 trace identity](https://github.com/leecrossley/theta-host-split-inert-race/blob/8fca7cad37185e3d7951f0dffff5e2442d1dfa24/HostAdmissible.lean)

```math
L_i(y)=S_{i,b_i}(\sqrt2\,y),\qquad
D_i(y)=\frac{S_{i,b_i}(2y)-S_{i,b_i}(y)}{\log 2}.
```

There is a stronger joint result. For each race draw an ellipse centred at $`(1,0)`$ in its level–move plane, with the following horizontal and vertical semiaxes.

| Field | Level, width 2 | Move, width 2 | Level, any fixed width | Move, any fixed width |
|---|---:|---:|---:|---:|
| 547/727 | 9.9375 | 3.4686 | 6.6921 | 1.9735 |
| 571/599 | 7.9731 | 2.1866 | 6.1027 | 1.5894 |

Choose a point on or inside each ellipse. Beyond any specified height, a weighted average of **at most four common starts** matches all four coordinates exactly. The weights are nonnegative, sum to one and are the same for both races. For the two moves alone, any point of the corresponding rectangle takes at most two common starts. The two individual difference bounds in the opening paragraph each hold at single heights. [Convex-hull argument](PROOF.md#common-starts-and-the-number-of-readings)

Finite counts give examples where both normalized moves have magnitude greater than 10. In field order 547/727 then 571/599, the signs at starts 17, 97, 162 and 658 are respectively $`(-,+)`$, $`(-,-)`$, $`(+,-)`$ and $`(+,+)`$. These counts omit the primes over 2 and the two conductors in each race. [Direct counts](finite_examples.py)

## Verification

Use Python 3.10 or later with the pinned python-flint package. On Linux x86-64 with Git, curl, zstd and a C toolchain available, the script downloads the pinned Lean toolchain and Mathlib. With an existing checkout, pass its path as the sole argument and put Lean 4.33.1 on `PATH`.

```sh
python3 -m pip install -r requirements.txt
bash verify-release.sh
```

The replay checks the real segments, the nine sign changes, the mixed frequencies by two independent evaluations, the ellipse bounds, the finite examples and the six original Lean deductions plus four deadline deductions. Both evaluations use Arb ball arithmetic. The Mellin, positive-weight and convex-geometric arguments are written proofs. [Checksums](SHA256SUMS) · [Prior work](PRIOR-WORK.md)
