# Exactly three primes can be isolated from two dyadic prime races

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

The replay checks the manifest, exact rational ranks, the four field-specific trace exclusions, 720 single-field prime-norm witnesses, 13,911 labelled primes and the twelve real-character certificates. It also checks the constants in both estimates for the large-height range. It then replays seven Lean deductions in the kernel and compares the output with the archived transcript. A C++ compiler and GMP speed up modular exponentiation. The standard-library Python fallback verifies the same congruences. The global number-theory estimates and their application are in the written proof. [Expected output](evidence/check.txt) [Lean source](DyadicPrimeSupport.lean) [Prior work and scope](PRIOR-WORK.md) [Citation](CITATION.cff)
