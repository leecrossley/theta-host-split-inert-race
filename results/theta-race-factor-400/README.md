# Even 1,020H can miss a move in the 547/727 prime race

At $`b=r=2`$, both a large gain and a small loss occur after the explicit height $`H=e^{3855250213.808}/(2\sqrt2)`$. No doubling that starts after $`H`$ and finishes by $`2e^{6.235}H`$ loses more than $`1.9735\log2`$. The factor is **1,020.60123689658…**. In particular, 1,000 fails. There are arbitrarily late intervals of the same width missing either prescribed direction. [Proof of the four-digit obstruction](PROOF.md#a-four-digit-deadline-and-a-shallow-loss)

This comparison with the actual prime-ideal scores uses Dirichlet GRH for the eight odd factors and the base-field factors, including zeta. Alpöge's announced primitive Dirichlet theorem supplies that input. Two finite interval calculations cover the entire log window. [Critical-line statement](https://github.com/leanprover-community/riemann-hypothesis/blob/371a1a26d5e17f604428e3a8ce6ce341cf517434/src/Main.lean#L12) · [Independent upstream verification](https://math.mit.edu/~drew/grh/) · [Certificate](evidence/check.txt)

| Log centre or interval | Normalized doubling move |
|---|---|
| $`[3855250213.808,3855250220.043]`$ | Everywhere above $`-1.95`$ |
| $`3855250215.346`$ | Greater than $`12`$ |
| $`3855250215.954`$ | Less than $`-0.04`$ |

The move is $`(S(2y)-S(y))/\log2`$ at log centre $`x=\log(2\sqrt2\,y)`$. A start must be strictly greater than $`H`$, and its doubling must finish by the stated deadline. The interval proof checks 13,946 genuine zero brackets, a 32,769-node Fourier calculation and a second cover of 2,926 closed cells. [Extended bracket table](evidence/thousand-brackets.tsv) · [Independent cell proof](thousand_cells.py) · [Prior work](PRIOR-WORK.md)

The failure also occupies a positive fraction of log heights. For each prescribed direction, the lower uniform log-height density of starts that miss it by a factor of 1,000 is at least $`2^{-306812}/50`$. A single set of integer shifts supplies paired obstructions in both directions. [Density proof](PROOF.md#how-often-the-thousand-fold-obstruction-returns)

## The earlier factor-400 certificate

The two earlier one-way windows and their original publication record remain available. [Original tagged release](https://github.com/leecrossley/theta-race-factor-400/releases/tag/v1.0.0) · [Original archived paper](https://doi.org/10.5281/zenodo.24411635)

**A factor of 400 fails in both directions at explicit starting heights.** In the 547/727 prime-ideal race with $`b=r=2`$, every eligible doubling gains at one height and every eligible doubling loses at another. The two certified log intervals have widths $`6`$ and $`45/8`$, both greater than $`\log200`$. [Theorem and proof](PROOF.md)

| Starting height $`H`$ | Deadline for $`2y`$, with $`y\gt H`$ | Change in $`S(2y)-S(y)`$ for every eligible start |
|---|---|---|
| $`e^{469417229.340}/(2\sqrt2)`$ | $`2e^6 H=806.857586\ldots H`$ | Gain of more than $`0.1\log2`$ |
| $`e^{484947092.225}/(2\sqrt2)`$ | $`2e^{45/8}H=554.544569\ldots H`$ | Loss of more than $`0.1\log2`$ |

There are also arbitrarily large starting heights with only gains by $`2e^6H`$, and others with only losses by that same deadline. Thus a factor guaranteeing a gain, or one guaranteeing a loss, of more than $`1.9735\log2`$ after **every** height must exceed $`2e^6`$. The previous obstruction was $`392.7397507\ldots`$. [Previous result](https://github.com/leecrossley/theta-race-distant-interval/blob/527fffe7b7808457b47a292144fea96b556158a8/PROOF.md) [Prior-work comparison](PRIOR-WORK.md)

The proof uses four cosine moments to bound every omitted zero and certifies all points of both intervals. A second calculation uses two real integrals, Cauchy–Schwarz and direct interval evaluations of the sines. The two comparisons with actual prime scores require Dirichlet GRH. [Zero brackets](evidence/zero-brackets.tsv) [Verification output](evidence/check.txt) [Upstream verification](https://math.mit.edu/~drew/grh/)

## Replay

In a Python 3.10+ virtual environment run

```sh
python -m pip install -r requirements.txt
python check_zeros.py
python paired_windows.py
python cell_check.py
python check_source.py
```

For a full replay on x86-64 Linux run `bash verify-release.sh` with Git, curl, tar and zstd installed. With Lean 4.33.1 on PATH and a clean, built Git checkout of Mathlib at commit `0df444a360eaa60ab8c11dca51a86af692955474`, run `bash verify-release.sh /path/to/mathlib4`. An optional second argument names the directory of the five upstream records. When the live audit logs change, the script checks archived copies against the published hashes. Set `THETA_WORKERS` to change the number of zero-check processes. The script checks file hashes, verifies the source records, repeats the numerical proofs and replays nine local Lean deductions. The explicit formula and its analytic consequences are written proofs. [Local Lean proofs](PairedWindows.lean) [File hashes](SHA256SUMS)
