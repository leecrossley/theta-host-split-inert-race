# When the 547/727 prime race identifies its conductors

The exact cutoff is **22,381,506,420.25**. Infinitely many other pairs of prime quadratic conductors produce exactly the same scores at every real starting height up to and including that number. Once observations at every real starting height extend past it, the 547/727 pair is uniquely determined. The comparison family has two distinct prime conductors congruent to three modulo four, with no size bound. Each race omits its two ramified primes and the same known auxiliary primes. The result is unconditional. It uses Dirichlet’s theorem on primes in progressions and, for the high tail, Theorem 1.2 of Bennett, Martin, O’Bryant and Rechnitzer. [Theorem and proof](PROOF.md)

Two exact rational measurements suffice. One must be inside each of these open intervals. Their disjointness makes two the minimum. The proof gives the corresponding criterion for every pair in the comparison family.

| Conductor | Starts where its absence can be detected |
|---|---|
| 727 | $`(727^2/4,727^2)=(132132.25,528529)`$ |
| 547 | $`(547^4/4,547^4)=(22381506420.25,89526025681)`$ |

Six measurements give a version that tolerates explicit rounding errors. A segmented sieve evaluated the six actual scores through norm **89,526,025,682**. Each score is enclosed to within two units in the twelfth decimal place. The sieve omits 547 and 727 by construction, so the recovered masses at $`547^4`$ and $`727^2`$ are zero by design. Each enclosure is narrow enough to distinguish zero from the unramified alternative. The error tolerance per reading is about $`7.0\cdot10^{-11}`$ at $`547^4`$ and $`1.25\cdot10^{-5}`$ at $`727^2`$. [Measured values](evidence/raw-scores.tsv) [Replay output](evidence/check.txt)

Any two fixed, distinct conductor pairs give different scores at **every sufficiently large rational height**. The proof gives a bound depending on both pairs. For any bounded set of real observations that misses either diagnostic interval, infinitely many other pairs agree with the 547/727 race at every observed height. This includes every bounded window starting at or above 528,529, the upper endpoint of the 727 interval. [Quantifiers and proof](PROOF.md#what-an-arbitrarily-high-tail-determines) [Prior work](PRIOR-WORK.md)

## Replay

Run these commands in this result directory using Python 3.10+ and a C++17 compiler with 128-bit integers.

```sh
python -m pip install -r requirements.txt
c++ -std=c++17 -O3 sieve.cpp -o /tmp/theta-conductor-sieve
python measure.py /tmp/theta-conductor-sieve
python independent.py
python horizon.py
python sampling.py
```

On x86-64 Linux with Git, curl, tar and zstd, `bash verify-release.sh` checks the file manifest, compares fresh arithmetic output with the archived files and replays seven Lean deductions. With Lean 4.33.1 on PATH and a built Mathlib Git checkout at `0df444a360eaa60ab8c11dca51a86af692955474`, run `bash verify-release.sh /path/to/mathlib4`. This mode verifies the Git revision and tracked Mathlib source and uses the supplied compiled dependency cache. The full identification theorem has a written proof. [Lean deductions](ConductorHorizon.lean) [File hashes](SHA256SUMS)
