# First opposite runs at powers of two in two prime races

At an exact power of two, the 547/727 prime race first gains twice while the 571/599 race loses twice at **17,179,869,184**, or $`2^{34}`$. Reverse the two fields and the first start is **4,096**, or $`2^{12}`$. Each of the four changes exceeds $`1.9735\log 2`$ in its prescribed direction. The 547/727 gains first occur at a starting height **4,194,304** times that for the 571/599 gains. Both statements follow from unconditional finite prime counts. [Theorem and proof](PROOF.md#opposite-runs-on-the-unshifted-dyadic-grid) [Two independent checks](evidence/dyadic-starts.txt)

This difference persists even if the required margin is reduced to zero and every move merely has to be strict. The first starts remain $`2^{12}`$ and $`2^{34}`$ for **every** common normalized margin below $`3.4143751084\ldots`$ and at least zero. The exact upper boundary is the 571/599 gain at $`2^{12}`$. At equality the early direction disappears from the whole tested range. The opposite direction survives up to the larger strict threshold $`6.756437069\ldots`$. [Exact endpoint conventions](PROOF.md#opposite-runs-on-the-unshifted-dyadic-grid)

At strict margin zero, **14 of the 16** two-step patterns occur at the 35 starts from 1 through $`2^{34}`$. Two finite sieves give every first start in this range for **every nonnegative margin**, including equality at the cutoffs. Sixteen distinct limiting moves govern the table. They and 1 are linearly independent over the real algebraic numbers. [All 23 record changes](evidence/dyadic-records.tsv) [Prime witnesses for independence](evidence/dyadic-pivots.tsv)

Run `bash verify-dyadic-starts.sh` from this directory to check the margins and their prime witnesses. Run `bash verify-dyadic-starts.sh --rebuild` to repeat both prime counts through norm **274,877,906,944** and compare their full integer output with the archive.

## One, two or three readings in a two-field prime race

Allowing three readings on 547/727 instead of two lowers the infimum of the common upper observation height by a factor of **136,174,087,494.48…**. With two readings on that field, the upper limit for observations must be strictly greater than **37,232,367,690,617,444,520.125**. With three, **273,417,420.125** is the sharp lower endpoint. Any larger height suffices in each case. The same access to the 571/599 field is available throughout. [Theorem and proof](PROOF.md#the-price-of-limiting-observations-on-one-field) · [Independent checks](evidence/observation-budget.txt)

The target is the one surviving linear combination of the early prime and integer gain thresholds. It can be written $`3\sqrt{46695}\,D_B(46695)-2\sqrt{46691}\,D_B(46691)`$, where $`D_B`$ is the normalized move on 571/599. Fix any lower cutoff $`46769\le H\lt46993`$, for example **46,800**. Both fields may be observed at algebraic heights from $`H`$ through the common upper endpoint. The table allows any finite number of 571/599 readings and algebraic linear combinations. It gives the smallest upper endpoint as an infimum.

| Available 547/727 readings | Infimum of the common upper endpoint | Is equality sufficient? |
| ---: | ---: | --- |
| 0 | Recovery impossible | |
| 1 | 37,592,788,570,517,198,340.125 | No |
| 2 | 37,232,367,690,617,444,520.125 | No |
| 3 or more | 273,417,420.125 | No |

Two readings on 547/727 can cancel one unwanted prime coefficient while leaving the desired coefficient. That choice introduces much larger rational primes whose effects on 571/599 begin only at one eighth of their **squared** norms on 547/727. A third reading isolates the original prime directly. The proof checks all possible one- and two-reading times, including the sharp endpoints. The theorem assumes no zero-location hypothesis. [Why the larger primes are necessary](PROOF.md#what-the-other-field-must-wait-for) · [306 auxiliary norms and 188 interval tests](observation_budget.py) · [Independent prime sieve](observation_budget_recheck.py)

## Prime and square starting heights

The first prime starting height at which both fields gain more than $`1.9735\log 2`$ twice in a row is **20,306,543**. The first perfect square is **20,313,049**, or $`4507^2`$. Both doublings use the same starting height on the 547/727 and 571/599 fields. Among unrestricted integer starts, the first is 20,306,534. Two finite covers exclude all earlier starts. They use independent prime sieves, different meshes and different error bounds for the exclusion. [Proof](PROOF.md#prime-and-square-starting-heights) · [Replay transcript](evidence/arithmetic-starts.txt)

Three early thresholds govern what happens as the required gain increases. For a lower cutoff below 223, either field alone recovers all three thresholds. The cutoffs beginning at 46,687 allow readings from both fields. If only readings from height 46,687 onward are available, two independent threshold values can still be recovered exactly. At 46,703 just one linear combination survives. At **46,993 none survives**, however far the observations continue. These are exact statements about finite linear combinations with real algebraic coefficients and readings at algebraic heights. The proof identifies the individual rational prime responsible for each loss. [Recovery theorem](PROOF.md#which-early-thresholds-can-later-readings-recover) · [Independent prime enumeration](recover_thresholds_recheck.py)

At every starting cutoff from 46,769 up to 46,993 (excluded), neither field alone can reconstruct any nonzero combination. Together they recover one if the permitted upper observation height is **strictly above 273,417,420.125**. Every larger endpoint suffices. Through 539,330,952.25, exactly three readings on 547/727 are necessary and sufficient when combined with readings on 571/599. The missing rational prime is 46,769, whose last weight on one field has vanished before its weight on the other field begins. [Why both fields and that height are necessary](PROOF.md#when-both-fields-are-needed)

| Allowed starting heights | Last successful early start as the margin approaches its threshold from below | Threshold, rounded down | First allowed start at the exact threshold |
| --- | ---: | ---: | ---: |
| Perfect squares | 46,656 | 1.645437640571873959 | 20,295,025 |
| Rational primes | 46,691 | 1.677966548153054396 | 20,288,251 |
| Positive integers | 46,695 | 1.682222888957865575 | 20,288,501 |

Across **three successive margin intervals**, the first square exceeds 434 times the first prime, then the first prime exceeds 434 times the first integer, then the first integer exceeds 434 times the infimum for real starts. The proof gives the first prime and square for **every** margin from 0 through 2. The late real crossings at all three tabulated thresholds are transcendental. These conclusions use finite prime sums and Baker’s theorem and assume no hypothesis on zeros. [First-start formulas](PROOF.md#prime-and-square-starting-heights) · [Prime-logarithm coefficients](PROOF.md#which-early-thresholds-can-later-readings-recover)

## First starts among real and integer heights

On the 547/727 and 571/599 fields, the first integer start for **two consecutive gains on both fields**, each larger than $`1.9735\log 2`$, is **20,306,534**. The second doubling finishes at **81,226,136**. Every other two-step joint sign pattern already occurs at an integer start by **1,469**, with all four moves larger than $`2\log 2`$ in their prescribed directions. Two finite prime counts exclude **every real start** up to the threshold crossing. The calculation is unconditional. [First-occurrence proof](PROOF.md#first-occurrence-on-two-consecutive-doublings) · [Certificate](evidence/first-catalogue.txt)

The required margin makes a measurable difference. At normalized margin $`1.6`$, four gains first become possible just above $`46{,}610.99018\ldots`$. The early maximum is $`1.682449287659191709\ldots`$, at a noninteger height. At that margin, strict gains are first possible just above $`20{,}288{,}514.06421\ldots`$. Integer starts lose the early opportunity at the smaller margin $`1.682222888957865575\ldots`$. Between these two margins, the first successful integer is **more than 434 times** the infimum for real starts. Both changes of regime are sharp. [Real and integer thresholds](PROOF.md#real-and-integer-first-starts)

| Required gain, in units of $`\log 2`$ | Infimum of real starts, rounded down to five decimals | First integer start |
| ---: | ---: | ---: |
| 0, with strict gains | 45,740.15111 | 45,741 |
| 1.6 | 46,610.99018 | 46,611 |
| 1.9735 | 20,306,533.59533 | 20,306,534 |
| 2 | 20,308,190.27109 | 20,308,191 |

For **every** normalized margin from 0 through 2, two prime-sum curves determine the infimum of successful real starts and the least successful integer. At the real infimum one gain equals the threshold. Near the late crossing at margin 1.9735, the first $`B`$ move exceeds 20. The second $`A`$ and $`B`$ moves exceed 48 and 38 respectively. At every algebraic margin in this range the exact real onset is transcendental. [Prime-sum branches and arithmetic type](PROOF.md#real-and-integer-first-starts)

At the integer threshold, the early four-gain window is only 0.333262779289365032… units wide. Three equally spaced samples per unit miss 0.02116621319… percent of possible clock phases. Four guarantee detection. If samples instead repeat every octave, the sharp count is 97,121 offsets per octave. [Optimal sampling proof](PROOF.md#how-many-samples-detect-the-early-gain) · [Independent checks](evidence/clock-phases.txt)

To replay just the first-occurrence theorem, install the Python dependencies below and run `bash verify-first-catalogue.sh`. It builds two C++17 sieves, checks both full interval covers and recomputes the early maximum at higher precision. [Direct sieve](word_sieve.cpp) · [Independent sieve](word_sieve_recheck.cpp) · [Whole-interval verification](word_cover.py)

## All 4,096 sign patterns on six consecutive doublings

Assume the generalised Riemann hypothesis (GRH) for the sixteen odd primitive Dirichlet functions in the 547/727 and 571/599 races. Alpöge’s announced Dirichlet GRH theorem supplies this input. [Announcement](https://alpo.ge/remark.pdf) [Independent formal checks](https://math.mit.edu/~drew/grh/)

Pick a gain or a loss for each of the 547/727 and 571/599 prime races. Do that six times. **Every one of the 4,096 sequences occurs starting at arbitrarily large powers of two.** Both races start at the same height, and the twelve prescribed changes each exceed $`4\log 2`$ in magnitude. The recurrence, the dyadic certificate and the bounds on longer words use GRH. [Proof](PROOF.md)

Write $`S_i(y)`$ for field $`i`$'s normalized width-two inert-minus-split score, and $`D_i(y)=(S_i(2y)-S_i(y))/\log 2`$. A word of length $`m`$ specifies the signs of both fields at $`y,2y,\ldots,2^{m-1}y`$. Each entry of the table holds at actual power-of-two starts.

| Consecutive doublings | Joint patterns realized | Strict lower bound for each normalized move in its prescribed direction |
|---:|---:|---:|
| 4 | 256 | 14 |
| 5 | 1,024 | 9 |
| 6 | 4,096 | 4 |

There is a common set of integer translations that preserves all 5,376 witnesses, with lower uniform density at least $`2^{-252000000}`$. Each pattern has a fixed integer offset. Adding the same successful translation to every offset gives late starting exponents for the whole table. Two independent calculations verify 61,440 signed component bounds using 14,000 brackets for genuine zeros, an explicit remainder for the unlisted spectrum, which uses GRH, and exact integer phase arithmetic. [Dyadic witnesses](evidence/dyadic-starts.tsv) · [Grid certificate](PROOF.md#every-six-step-word-on-the-dyadic-grid) · [Independent replay](grid_independent.py)

The same proof gives a deadline in terms of the starting exponent for the whole catalogue, beyond a threshold that the proof does not make explicit. It also gives finite examples and bounds on longer words.

- Direct counts give all sixteen two-step patterns at tabulated positive rational starts below 47,000, with every normalized move exceeding 1.6 in magnitude. These finite examples use elementary splitting and prime counts, and need no hypothesis. [Rational witnesses](evidence/witnesses.tsv) · [Direct count](patterns.py)
- Under GRH, at margin 9, the number of different sufficiently late patterns is bounded by $`2^{241000}`$ independently of word length. At 547/727, 57 successive gains of that size are impossible sufficiently far out, as are 57 successive losses. At 571/599 the corresponding obstruction holds at 56. [Word-count and run bounds](PROOF.md#what-longer-patterns-can-occur)

## Verification

Use Python 3.11 or later and a C++17 compiler. The Linux x86-64 replay downloads pinned Lean and Mathlib. For an existing pinned Mathlib checkout, put Lean 4.33.1 on `PATH` and pass the checkout path to the script.

```sh
python3 -m pip install -r requirements.txt
bash verify-release.sh
```

The archive includes the checked upstream statement and verification records, 14,000 certified positive zero brackets, all the finite counts, independent local polynomial and primitive-character calculations, the constants, the arithmetic-grid and rank certificates and twenty-two Lean deductions. The explicit-formula, recurrence and density arguments are written proofs. [Verification transcript](evidence/check.txt) · [Checksums](SHA256SUMS) · [Prior work](PRIOR-WORK.md)
