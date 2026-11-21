# The last missing splitting patterns in two theta fields

For every real $`L\ge14493/4=3623.25`$, **all nine pairs of local splitting traces** on the 547/727 and 571/599 fields occur at primes in $`(L,4L/3]`$. The threshold is exact. The two fields last lack a common completely split prime between 2,879 and 4,831. At the threshold, 4,831 lies exactly at the right endpoint. [Last gaps](PROOF.md#where-the-last-two-gaps-close)

A prime supplies four quadratic signs, two for each field. The nine trace pairs remember the averages of the two pairs of signs. **All sixteen sign patterns occur in every such interval from $`L=26211/4=6552.75`$.** This threshold is also exact. The last missing pattern has signs $`(-,+,+,+)`$ in the conductor order 547, 727, 571, 599. Its consecutive primes are 6,047 and 8,737. The witnesses split in the common real cyclotomic field and avoid the fixed host omissions. [Local traces](PROOF.md#which-primes-and-which-traces)

For any factor $`C\ge4/3`$, a finite table identifies **every real starting height** $`L\ge1`$ where a specified sign pattern or trace pair is missing from $`(L,CL]`$. Read a row $`(p,q)`$ as the deficient interval $`[p,q/C)`$ whenever it is nonempty. The best factors that work from every $`L\ge1`$ are 337 for the four corner trace pairs, 431 for all nine trace pairs and 769 for all sixteen patterns. [Complete gap table](evidence/wide-gaps.tsv) · [Endpoint conventions](PROOF.md#the-complete-map-of-missed-deadlines)

The proof is unconditional. It verifies 66,016 labelled primes and uses a classical progression estimate beyond their range. The same estimate gives every sign pattern more than $`L/198`$ in prime log weight in $`(L,4L/3]`$ from $`\log L=1097`$. Write $`t=\log L`$. For $`t\ge6004`$ and interval width $`\Delta\ge L/t^{10}`$ the log weight of each pattern lies between $`\Delta/67`$ and $`\Delta/61`$. The narrower bounds $`\Delta/65`$ and $`\Delta/63`$ hold when $`\Delta\ge3L/t^{10}`$ or, at the original width, when $`t\ge12959`$. [Tail estimate](PROOF.md#all-sixteen-patterns-in-shrinking-intervals) [Prior work](PRIOR-WORK.md)

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

Lean checks the Proth primality implication, the endpoints of a missed deadline, the passage from consecutive primes to real starts and the subtraction of endpoint mass estimates. The prime lists, character partial sums and numerical inequalities are checked by the replay programs. The local splitting formula and the published estimate in arithmetic progressions are mathematical inputs to the written proof. [Lean deductions](FrobeniusWindows.lean) · [Prime checker](verify_windows.py) · [Analytic checker](verify_analysis.py)
