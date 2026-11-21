# The first missed gain and loss at factor 400 in the 547/727 prime race

Assume the eight odd primitive Dirichlet factors have no zero in the later region of the first-loss theorem. GRH for these factors implies this exclusion, and Alpöge's announced Dirichlet GRH theorem supplies it. [Announcement](https://alpo.ge/remark.pdf) [Independent formal checks](https://math.mit.edu/~drew/grh/) The first missed $`400H`$ deadline is then a gain. The first missed loss occurs at more than 26 times the centred log height of that gain. The score is the weighted split–inert prime-ideal race with $`b=r=2`$. Each move must exceed $`1.9735\log2`$ in magnitude and start strictly after $`H`$. [Two directional bounds](PROOF.md#the-first-missed-loss) · [First gain](PROOF.md#the-first-missed-move-at-factor-400)

| First missed direction | Bounds on the starting log height $`h=\log(2\sqrt2 H)`$ |
| --- | --- |
| Gain | $`1620542.262999\lt h\le1620542.276`$ |
| Loss | $`43380912.47\lt h\le77157950.44`$ |

At the upper loss endpoint, a whole interval of width 5.55 misses even $`514.475H`$. At the first missed loss the boundary doubling has value exactly $`-1.9735\log2`$. The ratio comparison $`26\lt h_-/h_+\lt 48`$ is for the logarithms. [Closed interval calculations](directional_loss.py)

**Unconditionally, any earlier missed loss through centred log height 43,380,912.47 would force an odd zero of real part above 0.5000000012 with absolute ordinate greater than 8,590 and at most 3.5 billion.** A separate, larger finite region supplies the upper loss endpoint. These regions have not been exhaustively checked. For every $`H\le e^{12500}`$ both moves are already guaranteed unconditionally. [Finite conditions](PROOF.md#the-first-missed-loss) · [Unconditional finite range](https://github.com/leecrossley/theta-host-split-inert-race/releases/tag/theta-race-waiting-bound-v1.3.0)

Two integer sweeps agree on 159,453,395 negative samples. Fourteen individual samples and two triangular averages complete the cover of real starting heights. Even the full sequence of dyadic cosine moments leaves a specific loss candidate at log centre 43,380,915.948 undecided among positive weighted spectra. [Triangular averages](loss_moments.py) · [Which centres the moments determine](PROOF.md#what-all-the-dyadic-moments-can-decide) · [Prior work](PRIOR-WORK.md#directional-first-failures)

## The first missed gain

Assume that the eight odd primitive Dirichlet factors of the 547/727 race have no zero in the explicit region described below. GRH for these functions implies this, and Alpöge's announced Dirichlet GRH theorem supplies it. [Announcement](https://alpo.ge/remark.pdf) [Independent formal checks](https://math.mit.edu/~drew/grh/)

The first starting height $`H_*\ge1`$ at which the dyadic 547/727 race misses a required move by $`400H_*`$ satisfies

```math
1620542.262999\lt\log(2\sqrt2\,H_*)\le1620542.276 .
```

The upper and lower heights differ by less than **1.31%**. The first miss is a **gain**. A loss greater than $`1.9735\log2`$ still finishes on time. The doubling that starts exactly at the boundary has gain **equal to** $`1.9735\log2`$. The deadline asks for a strictly later start. These conclusions need no zero-location input beyond that assumption. [First-failure proof](PROOF.md#the-first-missed-move-at-factor-400) · [Upstream theorem and audit](https://math.mit.edu/~drew/grh/)

**Unconditionally, a failure of this localization forces a Dirichlet zero in an explicit region.** One of those eight factors must then have a zero $`\beta+i\gamma`$ with $`\beta\gt0.50000004`$ and $`8590\lt|\gamma|\le e^{810330}`$. At $`|\gamma|\ge e^{420000}`$ the real part must exceed $`3/4`$. The base-field terms are bounded unconditionally. The certificate checks 13,946 positive zero brackets, all below height 1,365, so it does not search that region. A separate verification of zero exclusion throughout the much larger region remains open. [Exact finite condition](PROOF.md#a-second-witness-below-log-height-1620543) · [Verification transcript](evidence/check.txt)

Two integer Fourier sweeps certify 5,911,490 gain samples and 5,911,639 loss samples each. Five additional gains bridge the gaps. Open intervals around the witnesses cover **every real starting height** from log height 12,500 up to the lower endpoint, and the earlier unconditional result covers all lower heights. Two whole-interval calculations prove the large gain is absent on the later stretch. [Positive cover](first_cover.py) · [Negative interval calculations](early_400.py) · [Prior work](PRIOR-WORK.md#the-first-failure)

| Closed interval of log centres | Certified deadline factor | Upper bound for the normalized gain under the finite condition for that interval |
|---|---|---|
| $`[1620542.278,1620547.593]`$ | $`2e^{5.315}=406.729026\ldots`$ | $`1.944`$ |
| $`[3895686.221,3895691.520]`$ | $`2e^{5.299}=400.273146\ldots`$ | $`1.962`$ |

The second interval extends a published window that reached $`392.73975\ldots H`$. The first uses four moments to resolve an interior peak. An additional 23 closed cells extend its left edge to the upper endpoint of the first-failure bracket. [Earlier window](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-distant-interval/PROOF.md) · [New finite alternatives](PROOF.md#an-earlier-height-for-factor-400)

## The four-digit obstruction

Assume the generalised Riemann hypothesis (GRH) for the primitive Dirichlet functions in the 547/727 race and its base field, including the Riemann zeta function. Alpöge's announced Dirichlet GRH theorem supplies this input. [Announcement](https://alpo.ge/remark.pdf) [Independent formal checks](https://math.mit.edu/~drew/grh/)

At $`b=r=2`$, both a large gain and a small loss occur after the explicit height $`H=e^{3855250213.808}/(2\sqrt2)`$. No doubling that starts after $`H`$ and finishes by $`2e^{6.235}H`$ loses more than $`1.9735\log2`$. The factor is **1,020.60123689658…**. In particular, 1,000 fails. There are arbitrarily late intervals of the same width missing either prescribed direction. [Proof of the four-digit obstruction](PROOF.md#a-four-digit-deadline-and-a-shallow-loss)

GRH enters both the bound for the omitted zeros and the comparison with the actual prime-ideal scores, for the eight odd factors and the base-field factors. The zero brackets, the finite sums over the selected zeros and the Lean arithmetic need no hypothesis. Two finite interval calculations cover the entire log window. [Critical-line statement](https://github.com/leanprover-community/riemann-hypothesis/blob/371a1a26d5e17f604428e3a8ce6ce341cf517434/src/Main.lean#L12) · [Independent upstream verification](https://math.mit.edu/~drew/grh/) · [Certificate](evidence/check.txt)

| Log centre or interval | Normalized doubling move |
|---|---|
| $`[3855250213.808,3855250220.043]`$ | Everywhere above $`-1.95`$ |
| $`3855250215.346`$ | Greater than $`12`$ |
| $`3855250215.954`$ | Less than $`-0.04`$ |

The move is $`(S(2y)-S(y))/\log2`$ at log centre $`x=\log(2\sqrt2\,y)`$. A start must be strictly greater than $`H`$, and its doubling must finish by the stated deadline. The interval proof checks 13,946 genuine zero brackets, a 32,769-node Fourier calculation and a second cover of 2,926 closed cells. [Extended bracket table](evidence/thousand-brackets.tsv) · [Second cell proof](thousand_cells.py) · [Prior work](PRIOR-WORK.md)

The failure also occupies a positive fraction of log heights. For each prescribed direction, the lower uniform log-height density of starts that miss it by a factor of 1,000 is at least $`2^{-306812}/50`$. A single set of integer shifts supplies paired obstructions in both directions. [Density proof](PROOF.md#how-often-the-thousand-fold-obstruction-returns)

## The earlier factor-400 certificate

The two earlier one-way windows and their original publication record remain available. [Original tagged release](https://doi.org/10.5281/zenodo.24411635) · [Original archived paper](https://doi.org/10.5281/zenodo.24411635)

**A factor of 400 fails in both directions at explicit starting heights.** In the 547/727 prime-ideal race with $`b=r=2`$, every eligible doubling gains at one height and every eligible doubling loses at another. The two certified log intervals have widths $`6`$ and $`45/8`$, both greater than $`\log200`$. [Theorem and proof](PROOF.md)

| Starting height $`H`$ | Deadline for $`2y`$, with $`y\gt H`$ | Change in $`S(2y)-S(y)`$ for every eligible start |
|---|---|---|
| $`e^{469417229.340}/(2\sqrt2)`$ | $`2e^6 H=806.857586\ldots H`$ | Gain of more than $`0.1\log2`$ |
| $`e^{484947092.225}/(2\sqrt2)`$ | $`2e^{45/8}H=554.544569\ldots H`$ | Loss of more than $`0.1\log2`$ |

There are also arbitrarily large starting heights with only gains by $`2e^6H`$, and others with only losses by that same deadline. Thus a factor guaranteeing a gain, or one guaranteeing a loss, of more than $`1.9735\log2`$ after **every** height must exceed $`2e^6`$. The previous obstruction was $`392.7397507\ldots`$. [Previous result](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-race-distant-interval/PROOF.md) [Prior-work comparison](PRIOR-WORK.md)

The proof uses four cosine moments to bound every omitted zero and certifies all points of both intervals. A second calculation uses two real integrals, Cauchy–Schwarz and direct interval evaluations of the sines. The omitted-zero bounds and the two comparisons with actual prime scores require Dirichlet GRH. [Zero brackets](evidence/zero-brackets.tsv) [Verification output](evidence/check.txt) [Upstream verification](https://math.mit.edu/~drew/grh/)

## Replay

Run `python3 directional_loss.py` for the new finite tests, both integer sweeps, triangular moments, direct cells and the ambiguity estimate. The full release also checks the original 13,946 positive ordinate brackets.

In a Python 3.10+ virtual environment run

```sh
python -m pip install -r requirements.txt
python check_zeros.py
python paired_windows.py
python cell_check.py
python thousand_zeros.py
python thousand_window.py
python thousand_cells.py
python early_400.py
python first_cover.py
python directional_loss.py
python check_source.py
```

`thousand_zeros.py` reuses the 4,000 brackets that `check_zeros.py` checks, so run it after `check_zeros.py`.

`first_cover.py` and `directional_loss.py` compile `first_grid.cpp`, so this replay also needs a C++17 compiler supporting signed 128-bit integers.

For a full replay on x86-64 Linux run `bash verify-release.sh` with Git, curl, tar, zstd and a C++17 compiler supporting signed 128-bit integers installed. With Lean 4.33.1 on PATH and a clean, built Git checkout of Mathlib at commit `0df444a360eaa60ab8c11dca51a86af692955474`, run `bash verify-release.sh /path/to/mathlib4`. The family tree supplies the five pinned upstream records in the [geometric-grid evidence directory](https://github.com/leecrossley/theta-host-split-inert-race/tree/b9f9d9178f689c5aa9840eff4e61bdd75fe7268f/results/geometric-grid-deadlines/evidence/upstream). An optional second argument selects another directory with the same five records. When the live audit logs change, the script checks archived copies against the published hashes. Set `THETA_WORKERS` to change the number of zero-check processes. The script checks file hashes, verifies the source records, repeats the numerical proofs and replays seventeen local Lean deductions. The explicit formula and its analytic consequences are written proofs. [Local Lean proofs](PairedWindows.lean) [File hashes](SHA256SUMS)
