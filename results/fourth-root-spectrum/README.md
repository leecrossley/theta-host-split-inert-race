# How many zeros recover the two theta prime races?

The unsmoothed normalized doubling move on each of the 547/727 and 571/599 fields can be recovered **with absolute error below 0.0008 from the true Dirichlet zeros up to $`5{,}000{,}000\,e^{x/4}`$**, at every log centre $`x=\log(2\sqrt2\,y)\ge12494`$. This cutoff is proportional to the fourth root of the starting prime height $`y`$. Each field uses its eight odd primitive factors. The real and imaginary parts of every zero in range, with multiplicity, determine the finite sum. [Exact weights and error theorem](PROOF.md#a-pure-fourth-root-cutoff-with-polynomial-weights) · [Replay transcript](evidence/check.txt)

The weights use a quintic spline and its first two derivatives. The quadratic correction accounts for zeros off the critical line. No global zero-location hypothesis is needed. The original box average carried an extra factor of $`\sqrt{\log y}`$. The finite zero condition for the first missed $`400H`$ gain on 547/727 now runs through $`e^{405153}`$. Whether the entire stated region is empty remains a separate question. [First-miss application](PROOF.md#recovery-at-every-stated-height) · [Original box-average release](https://github.com/leecrossley/theta-host-split-inert-race/releases/tag/fourth-root-spectrum-v1.0.0)

At sufficiently large heights, even the ordinary partial sum with the original coefficients converges to the move with a cutoff smaller than $`y^{1/4}`$ by every fixed power of $`\log y`$. This statement has an asymptotic onset. There is also an exact criterion for saving a **fixed power of $`y`$** when bounding the discarded terms in absolute value. Such a saving exists precisely when the eight Dirichlet factors have a fixed zero-free strip next to real part one at all sufficiently high ordinates. [Direct truncation](PROOF.md#a-sharp-truncation-at-sufficiently-large-heights) · [Proof of the equivalence](PROOF.md#when-a-fixed-power-below-one-quarter-works)

The certificate checks 24 Fourier integrals and derivatives by outward quadrature, all numerical error budgets, and the breakpoint identity against exact rational prime sums at eight starts in each field. Ten local Lean deductions cover the branch cancellation, moment and density algebra and the final scalar implications. The analytic Fourier and zero-density arguments are written proofs. The release supplies error estimates and weights for a complete finite zero list. It does not contain that list at these heights. [Numerical checks](quintic.py) · [Prime sums](prime_check.py) · [Prior work](PRIOR-WORK.md)

## Replay

With Python 3.10 or later in a virtual environment run

```sh
python -m pip install -r requirements.txt
python spectrum.py
python quintic.py
python prime_check.py
```

For the full check on x86-64 Linux, including the pinned Lean toolchain and Mathlib, run `bash verify-release.sh`. With Lean 4.33.1 on PATH and a clean, built Mathlib at commit `0df444a360eaa60ab8c11dca51a86af692955474`, run `bash verify-release.sh /path/to/mathlib4`. [Lean deductions](FourthRoot.lean) · [File hashes](SHA256SUMS)
