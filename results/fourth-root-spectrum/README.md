# Recovering the two theta prime races from zeros up to a fourth-root bound

For the 547/727 and 571/599 fields, every normalized doubling move at log centre $`x\ge12494`$ can be reconstructed from a finite weighted sum of odd Dirichlet zeros with **absolute error below 0.0008**, unconditionally. The sum uses each zero’s actual complex value. Each field uses just its eight odd primitive factors. The zero cutoff is

```math
T(x)=8\cdot10^6\sqrt{\lceil x\rceil+1}\,e^{x/4},
\qquad x=\log(2\sqrt2\,y).
```

For a start at height $`y`$ this grows as $`y^{1/4}\sqrt{\log y}`$. The reconstructed quantity is the original unsmoothed $`[S(2y)-S(y)]/\log2`$. The cutoff and weights are explicit. The theorem requires a complete finite list of complex zeros as input to a numerical reconstruction. The released calculations verify its error bounds. [Formula and proof](PROOF.md#the-fourth-root-formula) · [Numerical replay](evidence/check.txt)

The first missed $`400H`$ move on 547/727 is already localized to a height ratio below 1.0131 using a finite zero condition. **Its ordinate cutoff can now be reduced from $`e^{810330}`$ to $`e^{405161}`$.** At the other published $`400.273H`$ window the cutoff falls from $`e^{1947900}`$ to $`e^{973947}`$. If the first localization fails, one of the eight odd factors unconditionally has a zero in the smaller stated region. Whether that region is empty remains a separate input to a purely finite proof of the localization. [Smaller finite conditions](PROOF.md#the-first-missed-move-at-factor-400) · [First-failure certificate](https://github.com/leecrossley/theta-host-split-inert-race/releases/tag/theta-race-factor-400-v1.2.0)

On each smooth piece the prime tent is a sum of two exponentials with rates $`1/2`$ and $`-1/2`$. Symmetric averaging multiplies both by the same factor. Only prime norms near a breakpoint contribute to the resulting error. This allows strong damping of high zeros while recovering the pointwise unsmoothed move. [Prime-side argument](PROOF.md#what-the-symmetric-average-changes) · [Checks directly on prime ideals](prime_check.py) · [Prior work](PRIOR-WORK.md)

## Replay

In a Python 3.10+ virtual environment run

```sh
python -m pip install -r requirements.txt
python spectrum.py
python prime_check.py
```

Run `bash verify-release.sh` for the full replay on x86-64 Linux with Git, curl, tar and zstd installed. With Lean 4.33.1 on PATH and a clean, built Mathlib checkout at `0df444a360eaa60ab8c11dca51a86af692955474`, run `bash verify-release.sh /path/to/mathlib4`. The six local Lean deductions cover the branch cancellation, its rational form, the variance algebra, the cutoff and the scalar error implications. The explicit formula, prime estimates and convolution argument are written proofs. [Lean file](FourthRoot.lean) · [Manifest](SHA256SUMS)
