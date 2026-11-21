# How much of one theta race can the other predict?

Give a fixed linear predictor access to both races. It can use the last 15 doublings in the target field and any fixed linear average of the other field's scores across a bounded range of log heights. The predictor chooses the range, the signed weights and the intercept. It may even use future values of the other field.

On **547/727**, the resulting logarithmic root-mean-square error is strictly above **13.338** for the next level and **20.198** for the next normalized doubling move. On **571/599**, the corresponding bounds are **13.313** and **19.995**. These are lower bounds for the limiting upper RMS. They also hold when the starting heights are restricted to cubes, or to any fixed positive power of the integers, using harmonic weights. The bounds are unconditional. [Proof](PROOF.md) · [Verification transcript](evidence/check.txt)

| Previous target values | 547/727 level | 547/727 move | 571/599 level | 571/599 move |
|---:|---:|---:|---:|---:|
| 0 | 14.419 | 26.294 | 14.586 | 26.690 |
| 1 | 14.124 | 24.533 | 14.304 | 25.005 |
| 3 | 14.087 | 22.590 | 14.245 | 22.959 |
| 7 | 13.966 | 21.535 | 14.220 | 21.613 |
| 15 | 13.338 | 20.198 | 13.313 | 19.995 |
| 31 | 11.311 | 16.265 | 12.686 | 18.542 |
| 63 | 4.967 | 7.754 | 5.732 | 8.325 |
| 127 | 1.500 | 2.180 | 1.750 | 2.500 |

A level is the width-two weighted inert-minus-split prime-ideal score. A normalized move is $`(S(2y)-S(y))/\log 2`$. Each table entry allows every compactly supported finite signed average of the other field's levels. With 15 previous target values the error is measured at the sixteenth dyadic observation. The starts $`y=n^k`$ select the first observation in each vector. [Definitions and sampling](PROOF.md#the-two-fields-and-the-error)

The certificate locates **710 distinct field-specific frequencies**. At each one, the target field has a critical-line zero and every one of the other field's eight factors is nonzero. Consequently any choice of weights on the other field leaves a measurable part of the error. Under the critical-line hypothesis for the sixteen factors, the proof also identifies the exact least mean-square error from the other race with the total energy at all such unshared frequencies. Past observations alone attain the same infimum. [Exact distance](PROOF.md#the-exact-distance-when-the-two-spectra-are-on-the-line) · [Scope and prior work](PRIOR-WORK.md)

## Reproduction

Use Python 3.10 or later and the pinned dependencies. The Linux x86-64 replay obtains Lean and Mathlib. Other platforms can supply the pinned Mathlib directory and the matching Lean toolchain on `PATH`.

```sh
python3 -m pip install -r requirements.txt
bash verify-release.sh
```

Both exclusion methods check 5,680 factor intervals. One RMS calculation solves interval linear systems. A second uses exact rational prediction recursions after accounting for all input uncertainty. Lean checks six algebraic deductions, and the analytic arguments are given in the proof. [Bracket data](evidence/brackets.tsv) · [RMS bounds](radii.py) · [Rational recursion](prediction_recursion.py) · [File hashes](SHA256SUMS)
