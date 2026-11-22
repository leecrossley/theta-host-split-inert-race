# Which primes survive a shift of the dyadic sampling grid?

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

Allowing all **real** phases gives a different conclusion. For each odd prime $`p`$ there is a countable dense set of transcendental phases at which a rational multiple of the single 547/727 reading at $`256\alpha`$ equals $`\log_2 p`$. This is a numerical equality of physical readings, not an isolated column: the reading's formal prime coefficients still involve at least two primes. The union of phases where any prime is exactly recoverable by any finite algebraic combination is still countable and has measure zero. The phase in the construction depends on the target prime. These statements and the algebraic phase table are unconditional. [Real phases](PROOF.md#what-changes-at-transcendental-phases)

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

The replay checks the manifest, exact rational ranks, the four field-specific trace exclusions, 720 single-field prime-norm witnesses, 13,911 labelled primes and the twelve real-character certificates. It also checks the constants in both estimates for the large-height range, the phase table, 21,432 labelled prime witnesses and all 278,616 systems by two modular eliminations. It then replays nine Lean deductions in the kernel and compares the output with the archived transcript. A C++ compiler and GMP speed up modular exponentiation. The standard-library Python fallback verifies the same congruences. The global number-theory estimates and their application are in the written proof. [Expected output](evidence/check.txt) · [Lean source](DyadicPrimeSupport.lean) · [Prior work and scope](PRIOR-WORK.md) · [Citation](CITATION.cff)

With SymPy 1.14.0 installed, the independent polynomial-matrix replay runs with

```sh
python3 phase_smith_recheck.py
```

It recomputes the Smith invariant factors, real-root counts and endpoint ranks. [Independent transcript](evidence/phase-smith.txt)
