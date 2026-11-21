# When the 547/727 prime race identifies its conductors

The exact cutoff is **22,381,506,420.25**. Infinitely many other pairs of prime quadratic conductors produce exactly the same scores at every real starting height up to and including that number. As soon as observations extend past it, the 547/727 pair is uniquely determined. The comparison family has two distinct prime conductors congruent to three modulo four, with no size bound. Each race omits its two ramified primes and the same known auxiliary primes. [Theorem and proof](PROOF.md)

Two exact rational measurements suffice. One must be inside each of these open intervals. Their disjointness makes two the minimum. The proof gives the corresponding criterion for every pair in the comparison family.

| Conductor | Starts where its absence can be detected |
|---|---|
| 727 | $`(727^2/4,727^2)=(132132.25,528529)`$ |
| 547 | $`(547^4/4,547^4)=(22381506420.25,89526025681)`$ |

Six measurements give a version that tolerates explicit rounding errors. A segmented sieve evaluated the six actual scores through norm **89,526,025,682**. Each score is enclosed to within two units in the twelfth decimal place. Together they exclude the unramified alternative at both primes. [Measured values](evidence/raw-scores.tsv) · [Replay output](evidence/check.txt)

There is a different conclusion high in the race. Any two fixed, distinct conductor pairs give different scores at **every sufficiently large rational height**, with an explicit bound depending on both pairs. Yet for every bounded observation window lying above 89,526,025,681, infinitely many pairs mimic the 547/727 race over the entire window. [Quantifiers and proof](PROOF.md#what-an-arbitrarily-high-tail-determines) · [Prior work](PRIOR-WORK.md)

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

On x86-64 Linux with Git, curl, tar and zstd, `bash verify-release.sh` checks the file manifest, regenerates the arithmetic evidence and replays seven Lean deductions. With Lean 4.33.1 on PATH and a built Mathlib Git checkout at `0df444a360eaa60ab8c11dca51a86af692955474`, run `bash verify-release.sh /path/to/mathlib4`. This mode verifies the Git revision and tracked Mathlib source and uses the supplied compiled dependency cache. The full identification theorem has a written proof. [Lean deductions](ConductorHorizon.lean) · [File hashes](SHA256SUMS)
