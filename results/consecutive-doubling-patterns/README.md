# All 4,096 sign patterns on six consecutive doublings

Pick a gain or a loss for each of the 547/727 and 571/599 prime races. Do that six times. **Every one of the 4,096 sequences occurs starting at arbitrarily large powers of two.** Both races start at the same height, and the twelve prescribed changes each exceed $`4\log 2`$ in magnitude. The recurrence uses Alpöge's primitive Dirichlet critical-line theorem (Dirichlet GRH), applied to the sixteen odd characters. [Proof](PROOF.md) · [Upstream verification](https://math.mit.edu/~drew/grh/)

Write $`S_i(y)`$ for field $`i`$'s normalized width-two inert-minus-split score, and $`D_i(y)=(S_i(2y)-S_i(y))/\log 2`$. A word of length $`m`$ specifies the signs of both fields at $`y,2y,\ldots,2^{m-1}y`$. Each entry of the table holds at actual power-of-two starts.

| Consecutive doublings | Joint patterns realized | Strict lower bound for each normalized move in its prescribed direction |
|---:|---:|---:|
| 4 | 256 | 14 |
| 5 | 1,024 | 9 |
| 6 | 4,096 | 4 |

There is a common set of integer translations that preserves all 5,376 witnesses, with lower uniform density at least $`2^{-252000000}`$. Each pattern has a fixed integer offset. Adding the same successful translation to every offset gives late starting exponents for the whole table. Two independent calculations verify 61,440 signed component bounds using 14,000 brackets for genuine zeros, an explicit remainder for the unlisted spectrum, and exact integer phase arithmetic. [Dyadic witnesses](evidence/dyadic-starts.tsv) · [Grid certificate](PROOF.md#every-six-step-word-on-the-dyadic-grid) · [Independent replay](grid_independent.py)

The same proof gives a deadline in terms of the starting exponent for the whole catalogue. It also gives finite examples and bounds on longer words.

- Direct counts give all sixteen two-step patterns at tabulated positive rational starts below 47,000, with every normalized move exceeding 1.6 in magnitude. These finite examples use elementary splitting and prime counts. [Rational witnesses](evidence/witnesses.tsv) · [Direct count](patterns.py)
- At margin 9, the number of different sufficiently late patterns is bounded by $`2^{241000}`$ independently of word length. At 547/727, 57 successive gains of that size are impossible sufficiently far out, as are 57 successive losses. At 571/599 the corresponding obstruction holds at 56. [Word-count and run bounds](PROOF.md#what-longer-patterns-can-occur)

## Verification

Use Python 3.11 or later and a C++17 compiler. The Linux x86-64 replay downloads pinned Lean and Mathlib. For an existing pinned Mathlib checkout, put Lean 4.33.1 on `PATH` and pass the checkout path to the script.

```sh
python3 -m pip install -r requirements.txt
bash verify-release.sh
```

The archive includes the checked upstream statement and verification records, 14,000 certified positive zero brackets, all the finite counts, independent local polynomial and primitive-character calculations, the constants and eight Lean deductions. The explicit-formula, recurrence and density arguments are written proofs. [Verification transcript](evidence/check.txt) · [Checksums](SHA256SUMS) · [Prior work](PRIOR-WORK.md)
