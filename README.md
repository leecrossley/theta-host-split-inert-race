# A factor of 392 fails in the 547/727 prime race

**The exact waiting factor $`C=2e^{132/25}=392.739750703596\ldots`$ fails.** In the 547/727 race with $`b=r=2`$, there are arbitrarily large starting heights $`H`$ where every doubling that finishes by $`CH`$ misses the gain of more than $`1.9735\log 2`$. At other arbitrarily large heights the same deadline misses the loss. The previous lower bound on a universal factor was $`18.286280578160\ldots`$. [Theorem and proof](PROOF.md) · [Previous interval](https://github.com/leecrossley/theta-race-profile-interval/blob/264a0976afd3a59266400cf8040dbe700692dda8/PROOF.md)

One starting height for the missed gain is explicit. Take $`H_*=\exp(3895686.23)/(2\sqrt2)`$. The closed interval of log centres runs from $`3895686.23`$ to $`3895691.51`$. Its length is exactly $`5.28`$. The normalized limiting score stays below $`1.8`$ throughout. The prime-only score differs from it by less than $`10^{-20}`$ there. [Explicit height and recurrence](PROOF.md#from-the-profile-to-prime-ideals)

The certificate brackets 5,200 positive zeros of eight primitive Dirichlet functions. A compact test bounds the total contribution of every omitted zero, including any roots missed between tabulated ones. Two calculations check the full interval. One interpolates a finite trigonometric sum with a curvature bound. The other integrates the gamma term and encloses the sine terms directly over closed cells. The unconditional result uses the independently checked primitive Dirichlet zero theorem. [Zero brackets](evidence/zero-brackets.tsv) · [Verification output](evidence/check.txt) · [Upstream verification](https://math.mit.edu/~drew/grh/) · [Prior work](PRIOR-WORK.md)

## Replay

In a Python 3.10+ virtual environment run

```sh
python -m pip install -r requirements.txt
python check_zeros.py
python spectral_bound.py
python interval_check.py
python check_source.py
```

The zero check uses up to 12 worker processes by default. Set `THETA_WORKERS` to select a different number. For the full replay on x86-64 Linux, run `bash verify-release.sh` with Git, curl, tar and zstd installed. With Lean 4.33.1 on PATH and a built Git checkout of Mathlib at commit `0df444a360eaa60ab8c11dca51a86af692955474` with no local changes, run `bash verify-release.sh /path/to/mathlib4`. An optional second argument gives the directory containing the five pinned upstream records. The script checks the hashes, reproduces the numerical evidence, checks the upstream records, compiles five local Lean deductions and compares the output with the release receipt. The explicit formula and its analytic estimates are proved in writing. [Lean deductions](DistantInterval.lean) · [File hashes](SHA256SUMS)
