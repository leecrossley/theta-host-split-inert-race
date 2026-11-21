# Both 547/727 prime-race moves finish by 400H through exp(2000)

After every real starting height $`1\le H\le e^{2000}`$, the 547/727 prime-ideal race gains more than $`1.9735\log2`$ on one doubling and loses more than $`1.9735\log2`$ on another. Both starts are strictly after $`H`$ and both ends are **before $`400H`$**. Through $`e^{1000}`$, **$`256H`$** suffices. These bounds are unconditional. They use Platt's published finite verification of Dirichlet zeros and bound every higher zero wherever it lies in the critical strip. [Finite-height proof](PROOF.md#both-moves-through-exp2000) · [Verification transcript](evidence/check.txt)

| Starting heights | Deadline for both completions | Magnitude in each direction |
|---|---|---|
| $`1\le H\le10^9`$ | At most $`18H`$ by direct counting | More than $`1.9735\log2`$ |
| $`1\le H\le e^{1000}`$ | Before $`256H`$ | More than $`1.9735\log2`$ |
| $`1\le H\le e^{2000}`$ | Before $`400H`$ | More than $`1.9735\log2`$ |
| Every $`H\ge1`$ | Before $`e^{219}H^{17001/16999}`$ | More than $`3.32\log2`$ |

For every fixed $`\delta\gt0`$, both $`3.32\log2`$ moves also finish before $`H^{1+\delta}`$ for all sufficiently large $`H`$. This last statement gives no numerical starting height as a function of $`\delta`$. Separately, the all-height bound in the table puts the two stronger moves before $`H^{1.001}`$ once $`\log H\ge250000`$. [Shorter all-height bound](PROOF.md#a-shorter-all-height-bound) · [Arbitrarily small excess exponent](PROOF.md#every-positive-excess-exponent)

The finite proof checks 526 signed averages for the first range and 1,017 for the second. Each average has positive weights, so it locates at least one actual prime-score move of the stated size. Overlapping the intervals covered by the witnesses handles every real starting height. The four-digit obstruction at a much greater height assumes Dirichlet GRH, which Alpöge's announced theorem supplies. [First witness table](evidence/finite-1000.tsv) · [Second witness table](evidence/finite-2000.tsv) · [Four-digit obstruction](https://github.com/leecrossley/theta-host-split-inert-race/releases/tag/theta-race-factor-400-v1.1.0)

Run `bash verify-release.sh` to replay the complete paper on x86-64 Linux, or supply a built checkout of the pinned Mathlib with Lean 4.33.1 on the executable path. The release includes ten local Lean deductions. The analytic proof, interval computations and finite arithmetic are identified separately below. [File hashes](SHA256SUMS) · [Prior work](PRIOR-WORK.md)

## Earlier waiting interval

The original all-height proof and its direct prime counts are retained below. Its archived edition has its own citation. [Original release](https://github.com/leecrossley/theta-race-waiting-bound/releases/tag/v1.0.1) · [Original archived paper](https://doi.org/10.5281/zenodo.24410806)

After **every** starting height $`H\ge1`$, the split–inert lead on the 547/727 field both gains and loses more than $`1.9735\log 2`$ over a doubling. Both moves finish before

```math
e^{162000} H^{3001/2999}.
```

The exponent is about $`1.000667`$. The prefactor is large. For $`1\le H\le1{,}000{,}000`$, a direct count gives a much smaller bound of **$`18H`$** for both endpoints. For $`100\le H\le1{,}000{,}000`$, it gives $`(6851/657)H\lt 10.43H`$. At the other end of the scale, both moves finish before $`H^{1.001}`$ whenever $`\log H\ge500{,}000{,}000`$. [Theorem and proof](PROOF.md) · [Finite calculation](finite_windows.py)

## What is being measured?

Take $`K=\mathbf Q(\zeta_{16})^+`$, $`F=K(\sqrt{547\cdot727})`$ and $`E=K(\sqrt{-547},\sqrt{-727})`$. Count prime ideals of $`F`$, assigning $`+\log(Nv)`$ to an inert prime in $`E/F`$ and $`-\log(Nv)`$ to a split prime. The product of the normalized first-step theta trace jumps above each good prime gives its splitting sign. Use the good-prime convention in the proof. [Published trace identity](https://github.com/leecrossley/theta-race-short-shifts/blob/483c31f2df2132ed2165620ea55e20e969228bfb/RaceTrace.lean#L83)

For the score $`S(y)`$ multiply each signed log by $`W(Nv/y)`$ and divide the sum by $`2\sqrt2(\sqrt2-1)^2\sqrt y`$. The weight is $`2(1-1/u)`$ on $`1\lt u\le2`$, $`4/u-1`$ on $`2\lt u\lt 4`$ and zero elsewhere. A move is the difference $`S(2y)-S(y)`$. The theorem supplies $`y_+,y_-\gt H`$ with opposite moves exceeding the threshold in magnitude. Either move can happen first.

| Starting height $`H`$ | Verified gain from $`y`$ to $`2y`$ | Verified loss from $`y`$ to $`2y`$ |
|---:|---:|---:|
| 1 | 5 to 10 | 9 to 18 |
| 100 | 155 to 310 | 107 to 214 |
| 1,000 | 1,038 to 2,076 | 1,479 to 2,958 |
| 1,000,000 | 2,474,870 to 4,949,740 | 1,220,778 to 2,441,556 |

At one million the two changes divided by $`\log 2`$ lie in $`(4,5)`$ and $`(-16,-15)`$. The count certifies all real starting heights in the finite range, including heights equal to a sampled witness. [Calculation output](evidence/check.txt) · [Independent direct sums](independent.py)

## Why the bound is unconditional

The earlier short-shift proof established both moves arbitrarily far out. One zero of the conductor-547 Dirichlet function, near $`1/2+3.287178298086\thinspace i`$, has a coefficient larger than $`3.342843`$ in the normalized difference. A compact average isolates it on each prescribed interval. Wide averaging boxes separate the nearby frequencies. Increasing the number of narrow boxes as the starting height grows controls every zero above Platt's finite verification, including any zero off the line. Checking the base field's zeta factors in two complex disks accounts for the prime squares. The argument requires neither a global Riemann hypothesis nor linear independence of zero ordinates. The finite guarantees through $`e^{2000}`$ and the $`e^{219}`$ endpoint add Ramaré and Rumely's explicit prime counts, Rosser and Schoenfeld's bounds, the Montgomery–Vaughan large sieve and the 6,000-bracket zero certificate, none of which assumes where the zeros lie. [Spectrum and explicit formula](PROOF.md#an-isolated-zero-and-its-coefficient) · [Finite-height theorem](https://arxiv.org/pdf/1305.3087#page=14) · [Prior-work comparison](PRIOR-WORK.md)

## Replay

Use Python 3.10 or newer and the pinned version of python-flint. Each numerical test uses interval or exact arithmetic.

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python verify.py
python square_factors.py
python independent.py
python finite_windows.py
python finite_spectrum.py
python power_window.py
python carrier_window.py
python finite_deadlines.py direct
python finite_deadlines.py split
```

On x86-64 Linux with Git, curl, tar and zstd installed, run `bash verify-release.sh` in that environment for the full release replay. With Lean 4.33.1 on PATH and a built Mathlib at commit `0df444a360eaa60ab8c11dca51a86af692955474`, use `bash verify-release.sh /path/to/mathlib4`. The script verifies the file manifest, replays the finite and analytic calculations and the Lean module and compares the output byte for byte. The explicit formula and its error estimates are written mathematics. Lean checks the interval arithmetic in log coordinates and the deduction of both witnesses from the two weighted tests. [Lean statements](WaitingBound.lean) · [Expected output](evidence/check.txt) · [Release hashes](SHA256SUMS)
