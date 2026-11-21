# No compact prime weight erases either theta race

On the 547/727 and 571/599 fields, take any nonzero, continuous real weight supported in a fixed closed interval of positive prime-norm ratios. The weighted inert-minus-split prime count swings strictly above and below its square-term bias by a positive amount, arbitrarily far out. The amount can depend on the weight. The proof is **unconditional**. It uses Bauer's theorem on simple critical-line zeros and a real-axis check on the two fields. [Proof](PROOF.md) · [Bauer, Corollary 1](https://matwbn.icm.edu.pl/ksiazki/aa/aa93/aa9313.pdf#page=2)

**A consequence for forecasting.** Even with the previous fifteen dyadic scores, every fixed affine linear prediction of the next score misses by more than **1.394** at arbitrarily large starts on each field. The corresponding bound for predicting the next normalized move is **1.19**. The start can be restricted to $`n^k`$ for any fixed real $`k\gt0`$, including squares, cubes or hundredth powers. The coefficients and intercept are chosen once and held fixed.

| Coordinates, including the one predicted | Two-sided radius for levels | Two-sided radius for normalized moves |
|---:|---:|---:|
| 4 | 1.46 | 2.45 |
| 8 | 1.4 | 1.83 |
| 16 | 1.394 | 1.19 |
| 32 | 0.15 | 0.14 |
| 64 | 0.021 | 0.002 |
| 128 | 0.0059 | 0.00012 |
| 256 | 0.00136 | 0.00006 |

The table applies separately to both fields at width two. A move is $`(S(2y)-S(y))/\log2`$. For a vector of real coefficients of Euclidean norm one, the corresponding linear combination exceeds the stated radius in both directions around its bias on every tail. An arbitrary intercept can shift the two errors. The absolute-error assertion above holds for every intercept. Two independent calculations certify all 28 field, length and observable combinations. [Finite certificate](evidence/quadratures.tsv) · [Independent calculation](independent.py)

A second consequence identifies what a sparse observation can see. For a weight with a log kernel whose second derivative is a finite measure, the growth exponent in either direction is the supremum of the real parts of the surviving zeros, minus one half. This holds on irregular grids with gaps $`O(y^{1-\eta})`$ for any fixed $`\eta\gt0`$. For the ordinary width-two score on either field, an eventual upper bound alone, or an eventual lower bound alone, at the cubes is equivalent to the critical-line assertion for that field's eight odd Dirichlet factors. [Sparse sampling theorem](PROOF.md#what-survives-on-sparse-grids) · [Prediction bounds](PROOF.md#from-quadrature-to-prediction)

## Reproduction

Use Python 3.10 or later. On Linux x86-64 the replay downloads pinned Lean and Mathlib. On other systems, use Lean 4.33.1 on `PATH` and pass the existing pinned Mathlib directory to the script.

```sh
python3 -m pip install -r requirements.txt
bash verify-release.sh
```

The replay checks 710 critical-line brackets, the real segments for sixteen characters, 28 positive quadratures by two methods and six local Lean deductions. The analytic proofs of survival, oscillation and sampling are written arguments. [Verification transcript](evidence/check.txt) · [File hashes](SHA256SUMS) · [Prior work](PRIOR-WORK.md)
