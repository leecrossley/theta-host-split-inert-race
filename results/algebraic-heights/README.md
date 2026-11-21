# Transcendental crossing heights in two theta prime races

The 547/727 race hits both exact thresholds $`+1.9735\log 2`$ and $`-1.9735\log 2`$ after every starting height $`H\ge1`$, within one fixed multiplicative window. **Every height where it hits either threshold is transcendental.** Through $`H=e^{12500}`$ both exact hits finish by $`400H`$. Every finite interior boundary of a missed $`400H`$ deadline is also transcendental. The arithmetic conclusions are unconditional. Locating the first misses numerically uses the separate zero conditions of the directional result, which GRH implies and Alpöge's announced Dirichlet GRH theorem supplies. [Crossings and deadlines](PROOF.md#exact-crossings-and-deadlines) · [Finite deadline](https://github.com/leecrossley/theta-host-split-inert-race/releases/tag/theta-race-waiting-bound-v1.3.0) · [Directional first misses](https://github.com/leecrossley/theta-host-split-inert-race/releases/tag/theta-race-factor-400-v1.3.0)

Take the width-two score at $`m`$ distinct positive algebraic multiples of one sufficiently large algebraic start on both the 547/727 and 571/599 fields. Divide each score by $`\log 2`$. The $`2m`$ numbers and 1 are linearly independent over the algebraic numbers. For mixed algebraic widths and dilations, at every sufficiently large algebraic start, a matrix gives **every algebraic linear relation** among the readings. Its columns contain the three changes of slope of a tent. [Rank theorem](PROOF.md#the-rank-at-algebraic-starts) · [Exact width refinement](PROOF.md#width-refinement-and-exceptional-real-starts)

**The six readings at $`y,\sqrt2y,2y`$ have full rank exactly when $`y\gt79\sqrt2/8`$.** This holds at every positive algebraic start. Below and at that endpoint the first two readings on 571/599 use at most the logarithm of 7. Above it the rank is seven including 1. Exact determinants cover the first 901 cells. Six short chains of certified primes reach $`10^6`$, and two long chains bridge the finite calculation to the infinite range. [Sharp threshold and proof](PROOF.md#the-exact-transition-at-the-three-scales) · [Exact minors](evidence/initial-minors.tsv)

For every $`L\ge e^{366}`$ and each sign, the ordinary primes between $`L`$ and $`4L/3`$ whose two local splitting traces project above magnitude $`1/32`$ in that direction have total log weight greater than $`L/864`$. This holds for every unit projection and gives a quantitative count on the infinite range. [Prime bound and normalization](PROOF.md#an-explicit-infinite-range-and-a-prime-bound)

For either race, two rational starts in the same guaranteed window can bracket any nonzero real algebraic target of magnitude at most 5.85. The unique weight that gives this target by averaging the two normalised moves is transcendental. For a combination of the races with real algebraic coefficients $`v_A,v_B`$, the guaranteed range is $`0\lt |t|\le 5.85(|v_A|+|v_B|)`$. For the whole square of paired moves three rational starts suffice, and almost every target needs three. The corresponding sharp count for four coordinates is five. [Common windows and target range](PROOF.md#how-many-rational-starts-does-a-vector-target-need)

If an algebraic deadline fails but every rational height escapes a strict counterexample, the square of the logarithm of 2 must be a rational function of finitely many odd prime logarithms with algebraic coefficients. The proof gives its degree-two or degree-four equation. [Misses that a rational search could overlook](PROOF.md#what-an-isolated-missed-deadline-would-force)

The proof uses Baker’s theorem on logarithms, the prime number theorem in fixed progressions and Bordignon’s explicit progression bounds. It works directly with the prime-ideal scores. The rank is linear rank over the algebraic numbers. [Baker](https://doi.org/10.1112/S0025579300003971) · [Prior work](PRIOR-WORK.md)

## Replay

From this directory run

```sh
python3 verify_rank.py
bash verify-release.sh
```

The arithmetic checkers use the Python 3.10+ standard library. They verify exact determinants, 5,140 large prime witnesses, the local splitting tests, 43 real Dirichlet functions over a full interval and 36 geometric refinements of the tent. The full replay checks the hashes, compiles six local Lean deductions and replays their proofs in the kernel. On Linux x86-64 the script downloads Lean 4.33.1 and Mathlib revision `0df444a360eaa60ab8c11dca51a86af692955474` using curl, Git and tar with zstd support. It checks the download with sha256sum or shasum. On other platforms put Lean 4.33.1 on PATH and pass the directory of a clean, built checkout of that Mathlib revision to `verify-release.sh`. The trace script checks the moment identities and their numerical constants. The written proof uses Baker’s theorem, the prime number theorem, Bordignon’s progression bounds and the cited window theorems to deduce the prime bound and rank formula. [Lean deductions](AlgebraicHeights.lean) [Transcript](evidence/check.txt) [File hashes](SHA256SUMS)
