# How wide an average removes the slow modes of two theta prime races?

One positive average of the 547/727 and 571/599 prime races can cancel every nontrivial-zero contribution of ordinate at most one in absolute value. The shortest possible span on the log-height line is **between 103.945 and 103.949**. At the upper endpoint an explicit filter uses at most 34 shifts. At the lower endpoint every choice of nonnegative weights leaves a measurable response at one of the sixteen positive frequencies. All three width bounds are **unconditional**. [Proof](PROOF.md) · [Filter certificates](evidence/design.json)

| Frequencies to cancel | Positive ordinates in the entire band | Certified interval for the minimum log width | Pairs of shifts in the construction |
|---|---:|---:|---:|
| 547/727 | 9 | (88.436, 88.440) | 10 |
| 571/599 | 7 | (47.698, 47.702) | 8 |
| Both fields with the same weights | 16 | (103.945, 103.949) | 17 |

A positive average assigns total weight one to observations at shifted log heights. The two extreme shifts determine its width. One joint design saves more than **32.185 log units** compared with applying the two shortest separate averages successively. At the exact minimum every optimal average is supported on finitely many shifts. Smooth nonnegative weights can get arbitrarily close to that minimum and can achieve all three displayed upper widths. [Extremizers and smooth weights](PROOF.md#which-averages-attain-the-minimum)

For an average of log width at most 20, every fixed linear adjustment from the other field still leaves a level error with upper logarithmic RMS above **2.50** on 547/727 and **0.946** on 571/599. The corresponding bounds for the doubling change divided by $`\log2`$ are **0.1214** and **0.1162**. These bounds allow an arbitrary compactly supported signed adjustment and an intercept. [Prime-score consequence](PROOF.md#what-remains-visible-in-the-prime-scores)

## Replay

```sh
python -m pip install -r requirements.txt
python band_zeros.py
python euler_mass.py
python positive_widths.py
python independent.py
```

The replay checks 810 zero brackets, counts the complete low band in two ways, proves positivity on continuous intervals and solves for positive weights at the true zero ordinates. For the full numerical and Lean replay on x86-64 Linux, run `bash verify-release.sh`. With Lean 4.33.1 and the pinned built Mathlib checkout already installed, run `bash verify-release.sh /path/to/mathlib4`. Five local Lean lemmas cover the numerical deductions listed in the transcript. The zero-count identity, convex argument and connection to the actual prime scores have written proofs. [Transcript](evidence/check.txt) · [Lean deductions](AveragingWidths.lean) · [Prior work](PRIOR-WORK.md)
