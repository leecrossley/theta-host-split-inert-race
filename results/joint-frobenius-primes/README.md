# The last missing splitting patterns in two theta fields

For $`14493/4\le L\lt2^{1585}`$, **all nine pairs of local splitting traces** on the 547/727 and 571/599 fields occur at primes in $`(L,4L/3]`$. The common completely split primes 2,879 and 4,831 give the last missed interval in this range. At $`L=14493/4`$, its right endpoint is 4,831. Under Dirichlet GRH the same threshold works for every larger real start. [Last gaps](PROOF.md#where-the-last-two-gaps-close)

In the same certified range, **all sixteen quadratic sign patterns** occur from $`L=26211/4`$. The last missing pattern has signs $`(-,+,+,+)`$ in the order 547, 727, 571, 599. Its consecutive primes are 6,047 and 8,737. Their signs give the two trace coordinates by averaging in pairs. Under Dirichlet GRH all sixteen patterns continue at every larger start. [Local traces](PROOF.md#which-primes-and-which-traces)

For $`C\ge4/3`$ the table identifies every miss in $`(L,CL]`$ for $`1\le L\lt2^{1585}`$. A row $`(p,q)`$ gives the interval $`[p,q/C)`$ when it is nonempty. The best factors over this range are 337 for the four corners, 431 for all nine trace pairs and 769 for all sixteen patterns. They are also optimal over all starts under Dirichlet GRH. [Gap table](evidence/wide-gaps.tsv) [Endpoint conventions](PROOF.md#the-complete-map-of-missed-deadlines)

The finite proof uses 66,016 certified primes. The prime number theorem in fixed progressions also gives every pattern at all sufficiently large starts. For an explicit overlap and quantitative bounds we use Dirichlet GRH. Alpöge's announced Dirichlet GRH theorem supplies this input. For $`t=\log L\ge366`$ and $`\Delta\ge L/t^{10}`$ the prime log weight of each pattern in $`(L,L+\Delta]`$ lies between $`\Delta/65`$ and $`\Delta/63`$. [Prime-weight estimate](PROOF.md#all-sixteen-patterns-in-shrinking-intervals) [Prior work](PRIOR-WORK.md) [Announcement](https://alpo.ge/remark.pdf) [Independent formal checks](https://math.mit.edu/~drew/grh/)

## Replay

Run the verification in this directory.

```sh
bash verify-release.sh
```

The replay uses Python 3.10 or later, a C++17 compiler and the pinned Lean toolchain. GMP accelerates the modular powers when its development files are available. The Python standard library provides the same exact test. On Linux x86-64 the script obtains the pinned Lean and Mathlib dependencies. With those already installed, run `bash verify-release.sh /path/to/mathlib` using that revision. [Verification transcript](evidence/check.txt) · [File hashes](SHA256SUMS) · [Citation record](CITATION.cff)

For a separate modular-arithmetic replay with Python, run

```sh
python3 verify_windows.py --python
```

Lean checks the Proth primality implication, the endpoints of a missed deadline, the passage from consecutive primes to real starts and the subtraction of endpoint mass estimates. The prime lists, character partial sums and numerical inequalities are checked by the replay programs. The local splitting formula and Grenié and Molteni's estimate under Dirichlet GRH are mathematical inputs to the written proof. [Lean deductions](FrobeniusWindows.lean) · [Prime checker](verify_windows.py) · [Analytic checker](verify_analysis.py)
