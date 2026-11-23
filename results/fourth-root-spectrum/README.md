# How many zeros recover the two theta prime races?

The unsmoothed normalised doubling move on each of the 547/727 and 571/599 fields can be recovered **with absolute error below 0.0008 from the nontrivial zeros of the odd Dirichlet factors up to $`5{,}000{,}000\,e^{x/4}`$**, at every log centre $`x=\log(2\sqrt2\,y)\ge12494`$. This cutoff is proportional to the fourth root of the starting prime height $`y`$. Each field uses its eight odd primitive factors. The real and imaginary parts of every nontrivial zero in range, with multiplicity, determine the finite sum. [Exact weights and error theorem](PROOF.md#a-pure-fourth-root-cutoff-with-polynomial-weights) [Replay transcript](evidence/check.txt)

The weights use a quintic spline and its first two derivatives. The two derivative terms account for zeros off the critical line. The reconstruction is unconditional. The polynomial weights vanish beyond the cutoff. A uniform Taylor bound controls the omitted odd terms throughout the critical strip. McCurley’s zero-free region and the bound for a possible exceptional real zero control the base-field contribution. The original box average carried an extra factor of $`\sqrt{\log y}`$. The finite zero condition for the first missed $`400H`$ gain on 547/727 now runs through $`e^{405153}`$. The localisation assumes that this region contains no zero of any of the eight odd factors. The 13,946 brackets recorded in the factor-400 certificate all lie below height 1,365. The required exclusion begins above 8,590. [First-miss application](PROOF.md#recovery-at-every-stated-height) [Original box-average release](https://github.com/leecrossley/theta-host-split-inert-race/releases/tag/fourth-root-spectrum-v1.0.0)

At sufficiently large heights, even the ordinary partial sum with the original coefficients converges to the move with a cutoff smaller than $`y^{1/4}`$ by every fixed power of $`\log y`$. This statement has an asymptotic onset. There is also an exact criterion for saving a **fixed power of $`y`$** when bounding the discarded terms in absolute value. Such a saving exists precisely when the eight Dirichlet factors have a fixed zero-free strip next to real part one at all sufficiently high ordinates. [Direct truncation](PROOF.md#a-sharp-truncation-at-sufficiently-large-heights) · [Proof of the equivalence](PROOF.md#when-a-fixed-power-below-one-quarter-works)

The certificate checks 24 Fourier integrals and derivatives by outward quadrature, all numerical error budgets, and the breakpoint identity against exact rational prime sums at eight starts in each field. Ten local Lean deductions cover the branch cancellation, moment and density algebra and the final scalar implications. The written proof gives the explicit formula, the Fourier and zero-density arguments, the strip Taylor estimate, the prime-side breakpoint estimate, the zero-mass estimates, the base-field bounds and the transfer to the first missed gain. The release supplies error estimates and weights for a complete finite zero list. It does not contain that list at these heights. [Numerical checks](quintic.py) [Prime sums](prime_check.py) [Prior work](PRIOR-WORK.md)

## Replay

With Python 3.10 or later in a virtual environment run

```sh
python -m pip install -r requirements.txt
python spectrum.py
python quintic.py
python prime_check.py
```

For the full check on x86-64 Linux, including the pinned Lean toolchain and Mathlib, run `bash verify-release.sh`. With Lean 4.33.1 on PATH and a clean, built Mathlib at commit `0df444a360eaa60ab8c11dca51a86af692955474`, run `bash verify-release.sh /path/to/mathlib4`. [Lean deductions](FourthRoot.lean) · [File hashes](SHA256SUMS)
