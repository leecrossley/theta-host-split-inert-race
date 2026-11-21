# Comparing two theta prime races at the same height

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

The replay checks the real segments, the nine sign changes, the mixed frequencies by two independent evaluations, the ellipse bounds, the finite examples and six Lean deductions. Both evaluations use Arb ball arithmetic. The Mellin, positive-weight and convex-geometric arguments are written proofs. [Checksums](SHA256SUMS) · [Prior work](PRIOR-WORK.md)
