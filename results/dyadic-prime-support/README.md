# How many dyadic sampling grids isolate one prime?

On either theta field, choose an algebraic starting phase $`1\le\alpha\lt2`$ and read the split–inert race at $`\alpha,2\alpha,4\alpha,\ldots`$. Allow a finite linear combination of readings with real algebraic coefficients. How many independently chosen grids does it take to recover exactly $`\log_2 p`$ for a prescribed prime $`p`$?

| Field | Exactly one grid | Exactly two grids | Exactly three grids |
| --- | --- | --- | --- |
| 547/727 | 31 | 7, 47, 113, 127, 223, 271 | Every other good odd prime |
| 571/599 | 7 | 3, 79, 113, 127, 239, 241 | Every other good odd prime |
| Both fields | 7, 31, 47 | 3, 79, 97, 113, 127, 223, 239, 241, 271 | Every other odd prime |

The grids may be chosen afresh for each prime. On 547/727 alone the largest prime accessible with two grids is 271. On 571/599 alone it is 241. **97 is the only prime that two grids can isolate jointly while neither field can isolate it with two grids on its own.** The two phases for 97 are $`3/2`$ and $`7/4`$. Even with both fields, 17 needs three. [Exact formula for 97](PROOF.md#two-grids-using-both-fields) · [Fourteen single-field formulas](evidence/minimum-grid-formulas.tsv)

Three **rational** phases always suffice. On an individual field a good prime is one of the primes retained in its trace data. Every odd prime is good on at least one field. If the collection of algebraic phases must be fixed before the prime is chosen, any finite collection can isolate only finitely many prime logs, even when both fields are available. All these statements are unconditional. [Definitions](PROOF.md#the-two-readings) · [Minimum number of grids](PROOF.md#how-many-grids-isolate-a-prescribed-prime) · [Fixed collections](PROOF.md#fixing-the-grids-in-advance)

The proof also closes a gap in the count of ordinary primes. For every real $`u\gt271/2`$ and every nonzero direction in the two-field trace plane, at least two primes in $`(u,2u)`$ have nonzero projection in that direction. The endpoint fails. For $`\log u\ge1095`$ there are more than $`u/(500\log(2u))`$ such primes. [Exact endpoint](PROOF.md#two-grids-using-both-fields) · [Analytic estimate](PROOF.md#an-unconditional-overlap-using-paired-zeros)

## Shifting one grid

Read the 547/727 and 571/599 split–inert prime races at $`\alpha,2\alpha,4\alpha,\ldots`$ and allow finite linear combinations with real algebraic coefficients. At the algebraic phase $`\alpha=3/2`$, **neither field can isolate any individual prime logarithm on its own. Together they isolate exactly 31 and 47.** At $`\alpha=\sqrt{7/3}`$ even the two fields together isolate none.

Here is the complete pooled answer at every real algebraic phase $`1\le\alpha\lt2`$. A prime $`p`$ in the table denotes recovery of $`\log_2 p`$. Every upper endpoint shown with $`\le`$ belongs to its row.

| Algebraic sampling phase | Odd primes isolated by the two fields |
| --- | --- |
| $`1\le\alpha\le79/64`$ | 7, 31, 47 |
| $`79/64\lt\alpha\le97/64`$ | 31, 47 |
| $`97/64\lt\alpha\le49/32`$ | None |
| $`49/32\lt\alpha\le31/16`$ | 7 |
| $`31/16\lt\alpha\lt2`$ | 7, 31, 47 |

The 547/727 field alone isolates 31 exactly on $`[1,47/32]\cup(31/16,2)`$. The 571/599 field alone isolates 7 exactly on $`[1,79/64]\cup(49/32,2)`$. No other odd prime is individually recoverable on either field at an algebraic phase. Every nonzero raw reading at any one algebraic phase is linearly independent of the others. [Complete phase proof](PROOF.md#moving-the-sampling-grid) · [Exact rank certificate](phase_rank.py)

For each odd prime $`p`$ there is a countable dense set of transcendental phases at which a rational multiple of the single 547/727 reading at $`256\alpha`$ equals $`\log_2 p`$. At every such phase the formal expansion of this reading still has at least two nonzero prime coefficients. The union of phases where any prime is exactly recoverable by any finite algebraic combination is still countable and has measure zero. The phase in the construction depends on the target prime. These statements and the algebraic phase table are unconditional. [Real phases](PROOF.md#what-changes-at-transcendental-phases)

## The original grid

Fix the 547/727 and 571/599 theta fields and read each split–inert prime race at the powers of two, with $`b=r=2`$. Which individual prime weights can a finite algebraic linear combination of those readings recover exactly?

| Readings available | Recoverable odd prime weights |
| --- | --- |
| 547/727 alone | $`\log_2 31`$ |
| 571/599 alone | $`\log_2 7`$ |
| Both fields | $`\log_2 7,\ \log_2 31,\ \log_2 47`$ |

These lists are complete at all heights. The coefficients can be any real algebraic numbers, and the combination can use arbitrarily many readings provided it has finite length. The logarithms are normalized by $`\log 2`$. All nonzero individual readings on the two fields, together with 1, are linearly independent over the real algebraic numbers. Through exponent $`n\ge3`$ the readings span dimension $`2n-3`$. [Definitions, recovery formulas and proof](PROOF.md)

The proof gives an exact minimum number of contributing primes for any prescribed last reading on either field. For two fields, the sharp asymptotic minimum has three possible leading constants. A quantitative step shows that every direction sees more than $`X/(144\log(2X))`$ degree-one primes between $`X`$ and $`2X`$ once $`\log X\ge2000`$. All claims are unconditional. [Minimum on one field](PROOF.md#minimum-support-on-one-field) · [Two-field asymptotic](PROOF.md#how-many-ordinary-primes-must-remain) · [Explicit prime estimate](PROOF.md#ordinary-primes-in-every-direction)

## Replay

From this directory on Linux x86-64 with Git, curl, zstd and Python 3.10 or newer, run

```sh
bash verify-release.sh
```

With Lean 4.33.1 on the path and the pinned Mathlib checkout already available, run

```sh
bash verify-release.sh /path/to/mathlib4
```

The replay verifies the file hashes and compares its output with the archived transcript. For the grid minima it checks fourteen single-field formulas, the five-term formula for 97, 1,077 witnesses covering real norm octaves, all 2,664 pooled exclusions, and two independent estimates for the large-height range. Thirteen local Lean deductions are replayed in the kernel. The analytic number theory and the application of the certificates are in the written proof.

The inherited checks cover the original dyadic table, the four field-specific trace exclusions, 720 norm witnesses, 13,911 labelled primes, twelve real-character bounds, the one-grid phase table and 278,616 systems checked twice on 21,432 prime witnesses. A C++ compiler and GMP accelerate modular powers. The Python standard library verifies the same congruences. [Expected output](evidence/check.txt) · [Lean source](DyadicPrimeSupport.lean) · [Prior work and scope](PRIOR-WORK.md) · [Citation](CITATION.cff)

With SymPy 1.14.0 installed, the independent polynomial-matrix replay runs with

```sh
python3 phase_smith_recheck.py
```

It recomputes the Smith invariant factors, real-root counts and endpoint ranks. [Independent transcript](evidence/phase-smith.txt)
