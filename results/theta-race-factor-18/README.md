# When the 547/727 prime race first misses an 18H deadline

The sharp factor for obtaining both a gain and a loss of more than $`1.9735\log 2`$ is **$`17.493385009687\ldots`$** over every real starting height from $`1`$ through $`e^{40.0507}/(2\sqrt2)`$, about $`8.75\times10^{16}`$. Every larger factor works throughout this range. Equality already fails at $`H=1`$. From $`10^9`$ to the new upper endpoint, $`16.5H`$ suffices. This result uses finite prime and zero checks and assumes no global zero hypothesis. [Proof](PROOF.md#the-first-missed-18h-deadline) · [Two checks of the positive averages](onset_spectrum.py) · [Replay output](evidence/check.txt)

For the first miss at $`18H`$, the proof gives the following finite test. Suppose all zeros of the eight specified Dirichlet functions with $`8590\lt|\Im\rho|\le5{,}000{,}000`$ lie on the critical line. The **first missed move is a loss**, and its starting height $`H_*`$ satisfies

```math
40.0507 < \log(2\sqrt2\,H_*) \le 40.0936.
```

The upper endpoint is about $`9.14\times10^{16}`$. The first failed start is transcendental and its doubling hits the negative threshold exactly. From **91,386,963,248,971,966 through 91,393,856,157,112,503** every integer start misses the loss, a run of 6,892,908,140,538. The same finite condition applies to this run. Alpöge's announced Dirichlet GRH theorem supplies this input. The standalone finite replay uses Platt’s published verification through 8,590 and leaves the extra range to five million open. A single qualifying loss from any integer in the displayed run would force an off-line zero in that finite range. [Finite upper test](PROOF.md#a-finite-zero-test-for-the-first-miss) · [Whole-interval programs](onset_window.py) · [Announcement](https://alpo.ge/remark.pdf) [Independent upstream check](https://math.mit.edu/~drew/grh/)

Under the same finite condition this first loss stays at exactly the same height for **every** factor above $`17.493385009687\ldots`$ through $`2e^{2.1973}=18.0013\ldots`$. At the sharp factor itself the first miss is at height one. [Why the starting height stays fixed](PROOF.md#one-first-failed-height-for-an-interval-of-factors)

The unconditional bound **$`19.95H`$ works for both signs through $`H=e^{41.33}/(2\sqrt2)`$**, about $`3.15\times10^{17}`$. Under the finite zero condition the sharp factor on this larger range is $`C_*=2H_{\rm back}/H_*`$, strictly between $`18.0013`$ and $`19.95`$. Here $`H_{\rm back}`$ is the boundary where strict losses resume. For every factor above the height-one optimum and at most $`C_*`$, the complete set of failures is the single closed interval $`[H_*,2H_{\rm back}/C]`$. At $`C_*`$ only one real starting height fails and every integer start succeeds. The sharp integer factor is transcendental and lies within $`2.28\times10^{-16}`$ of the real factor. [Every exceptional start](PROOF.md#the-next-loss-and-every-exceptional-start)

A further cross-check in the family checkout independently integrates the four gamma terms, enumerates the signed residue classes and reproduces the integer endpoints at 95 decimal places. Run `python3 independent_onset.py` from this result directory. [Cross-check](independent_onset.py) · [Output](evidence/independent-onset.txt)

## Exactly when 6H misses an earlier move


Between one million and one billion, **$`6H`$ misses a required move at exactly 5,372,773 integer starting heights**. The failures occupy three closed intervals of real starts. The first misses a gain and the next two miss a loss. Everywhere else both moves of more than $`1.9735\log2`$ finish in time. This uses direct prime arithmetic on the 547/727 field with $`b=r=2`$. [Theorem and proof](PROOF.md#all-failing-starts-from-one-million-to-one-billion)

| Missing move at $`6H`$ | First failing integer $`H`$ | Last failing integer $`H`$ | Count |
|---|---:|---:|---:|
| Gain | 6,288,627 | 6,768,844 | 480,218 |
| Loss | 20,077,386 | 20,356,800 | 279,415 |
| Loss | 562,481,499 | 567,094,638 | 4,613,140 |

A single factor $`C`$ guarantees both moves by $`CH`$ for all real starts $`10^6\le H\le10^9`$ exactly when $`C`$ exceeds $`6.458177331046\ldots`$. A gain sets that limit. The sharp factor for losses alone is $`6.083501163950\ldots`$. Equality fails in each case. With integer starting heights the corresponding factors are slightly smaller, $`6.458177149109\ldots`$ and $`6.083501148481\ldots`$. The proof also gives the full set of exceptions for every real factor $`C\ge6`$. [Exact endpoints](PROOF.md#all-failing-starts-from-one-million-to-one-billion) [Numerical replay](evidence/check.txt)

Two prime sums with independently bounded errors locate the six crossings. An interpolation bound covers every point between their boundary neighborhoods, including all intervening prime norms. [Integer sieve](sieve.cpp) · [Independent prime sums](near_roots.py) · [Interval cover](late_slabs.py)

## The sharp finite deadline from height one

A gain and a loss of more than $`1.9735\log 2`$ both finish by $`CH`$ after **every real start from 1 through one billion** precisely when $`C`$ exceeds $`17.493385009687\ldots`$. The infimum itself fails. At $`H=1`$ the first loss reaches the required magnitude with equality. Every larger factor gives it room to cross. These finite conclusions follow from prime arithmetic. [Exact statement and proof](PROOF.md#the-sharp-finite-deadline)

The answer changes abruptly when the required move grows. Write $`D(y)=(S(2y)-S(y))/\log 2`$ and $`\tau=-D(201007/2)=3.997034484766\ldots`$. Just below $`\tau`$ the sharp factor approaches $`17.848747144\ldots`$. **At $`\tau`$ it jumps to $`19.723571102\ldots`$.** The last qualifying loss before a long gap disappears when it just touches the threshold. In particular, $`18H`$ guarantees both strict moves of size $`t\log2`$ through one billion exactly for $`0\le t\lt\tau`$. [The disappearing loss](PROOF.md#the-loss-that-disappears)

| Required magnitude in units of $`\log 2`$ | Infimum of the completion factor |
|---|---:|
| $`1.9735`$ | $`17.493385009687\ldots`$ |
| $`\tau`$ | $`19.723571102\ldots`$ |
| $`4`$ | $`19.727081988\ldots`$ |
| $`4.5`$ | $`20.211893558\ldots`$ |
| $`5`$ | $`20.592958664\ldots`$ |

Two computations agree at all 31,776 knots and endpoints through $`y=214710`$. The billion-height cover uses the original sieve through prime norm 32,467,011,360. [Numerical certificate](sharp_deadline.py) · [Prime computations](prime_knots.py) · [Verification output](evidence/check.txt)

## Earlier factor 18 result

[Original factor 18 publication on Zenodo](https://doi.org/10.5281/zenodo.24410788)

From every real starting height $`1\le H\le1{,}000{,}000{,}000`$, a gain and a loss of more than $`1.9735\log2`$ both finish by $`18H`$. The sieve proves this finite statement unconditionally. Assume the generalised Riemann hypothesis (GRH) for the primitive Dirichlet $`L`$-functions in the 547/727 race and its base field, including the Riemann zeta function. Alpöge's announced Dirichlet GRH theorem supplies this input. **At arbitrarily large heights, $`18H`$ misses the loss. At other arbitrarily large heights, it misses the gain.** Even the slightly larger factor $`2e^{11/5}=18.0500\ldots`$ fails. [Theorem and proof](PROOF.md)

Each move is a doubling $`S(y)\to S(2y)`$ with $`y\gt H`$. The score measures weighted inert prime ideals minus weighted split prime ideals on the 547/727 theta-lift host, with both the width and the move fixed at two. The allowed starts lie in $`(H,CH/2]`$ when both endpoints must finish by $`CH`$. [Prime convention and weight](PROOF.md#two-statements-about-the-same-race)

## What the direct count shows

From one million through one billion, both moves actually finish by **$`7H`$**. The sieve covers norms through 32,467,011,360. It certifies the slopes at 625 integer starts. Checking the next strictly later witness at each breakpoint covers every real starting height in the range. At $`H=1{,}000{,}000{,}000`$ the two verified moves are

| Move | Start | End | Change divided by $`\log2`$ |
|---|---:|---:|---:|
| Gain | 1,016,187,975 | 2,032,375,950 | $`[33.1637,33.1639]`$ |
| Loss | 1,714,594,397 | 3,429,188,794 | $`[-2.7357,-2.7356]`$ |

The finite result uses direct prime-ideal arithmetic. [All slope enclosures](evidence/windows.tsv) · [Coverage verifier](finite.py)

## Why the pattern eventually breaks

At large heights the normalized change approaches a uniformly convergent sum over zeros. On one closed log interval of length $`2.2`$, that sum stays above $`-1.8`$ everywhere. Two computations certify the entire interval, including the point where the kernel crosses the origin. Its shape recurs arbitrarily far out. A factor of 18 gives potential doubling starts only $`\log 9=2.19722\ldots`$ units of log height, so a whole search interval can fit inside a recurring stretch with no qualifying loss. Reflection yields stretches with no qualifying gain. The first-miss theorem above now supplies an explicit start under a finite zero condition. [Recurrence argument](PROOF.md#transfer-the-entire-interval) [Verification output](evidence/check.txt)

The recurrence proof uses Alpöge's announced primitive Dirichlet GRH theorem, independently matched to its formal statement by Sutherland's checks. Of the original nine local Lean lemmas, `recurring_no_loss` assumes a bound on late translates. The lemma `move_fits` places the doubling inside the logarithmic interval. The other seven check interpolation, endpoints and deadlines in the finite calculation. The analytic recurrence and its application to gains and losses are written out in the proof. [Announcement](https://alpo.ge/remark.pdf) [Independent formal checks](https://math.mit.edu/~drew/grh/) [Prior work](PRIOR-WORK.md)

Under GRH for the primitive Dirichlet functions of the race and its base field, including the zeta function, the companion result proves that $`e^{175}`$ works at every height. [Fixed-factor proof](https://github.com/leecrossley/theta-race-fixed-factor/releases/tag/v1.0.1)

## Replay

Use Python 3.10 or newer inside a virtual environment and a C++17 compiler. Run

```sh
python -m pip install -r requirements.txt
python profile.py
python profile_series.py
python finite.py
python sharp_deadline.py
c++ -std=c++17 -O3 sieve.cpp -o /tmp/theta-factor-18-sieve
python finite.py /tmp/theta-factor-18-sieve
python late_slabs.py /tmp/theta-factor-18-sieve
python late_roots.py /tmp/theta-factor-18-sieve
python splitting_square_bound.py /tmp/theta-factor-18-sieve /tmp/theta-factor-18-bins
python splitting_square_recheck.py /tmp/theta-factor-18-bins
python onset_spectrum.py hardy
python onset_spectrum.py exterior
python onset_spectrum.py direct
python onset_spectrum.py split
python onset_window.py
python onset_window.py cells
python check_source.py
```

The first `finite.py` command checks coverage from the stored table. The second regenerates every enclosure. For the full x86-64 Linux replay with Git, curl, tar and zstd, use `bash verify-release.sh`. With Lean 4.33.1 on PATH and a built Mathlib at commit `0df444a360eaa60ab8c11dca51a86af692955474`, use `bash verify-release.sh /path/to/mathlib4`. An optional second argument supplies local copies of the five published source and audit files. The release replay also compiles the Lean implications, rechecks the module with the kernel and compares the full output against the receipt. [Lean deductions](FactorObstruction.lean) · [File hashes](SHA256SUMS)
