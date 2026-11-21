# What repeated factor-400 misses would say about the zeros

In the 547/727 dyadic prime race, suppose at least one of the eight odd Dirichlet factors has a zero off the critical line. **If $`400H`$ misses a qualifying gain or loss at arbitrarily large $`H`$, the rightmost real part is attained at 14,488 or more distinct positive ordinates.** More precisely, their reciprocal sum is at least $`\log(200)/\pi`$. The move threshold is $`1.9735\log2`$ and the same conclusion holds at every fixed positive threshold. [Theorem and proof](PROOF.md#what-repeated-misses-imply)

If the supremum is never attained, **every $`C\gt2`$ eventually works**. If it is attained at just one positive ordinate $`g`$, the exact infimum of the eventual completion factors is $`2e^{\pi/g}\lt2.000732`$. Multiplicity at that ordinate makes no difference. The proof covers every sufficiently large real starting height and both directions. It leaves the boundary value at the infimum open. [The single-ordinate case](PROOF.md#one-positive-ordinate)

For the actual field, the unconditional conclusion is an alternative. Arbitrarily late misses at 400 force either all eight odd factors to satisfy their critical-line assertion, or the large rightmost set above. A single finite miss does not imply this alternative. The earlier explicit large-height counterexample assumes Dirichlet GRH, which Alpöge's announced theorem supplies. Conversely the direct finite-height theorem shows $`400H`$ works through $`e^{2000}`$. [Earlier obstruction](https://github.com/leecrossley/theta-host-split-inert-race/releases/tag/theta-race-factor-400-v1.1.0) · [Finite positive result](https://github.com/leecrossley/theta-host-split-inert-race/releases/tag/theta-race-waiting-bound-v1.1.1)

**A gain and a loss above the stated magnitude both finish by $`2.00125H`$ at arbitrarily large $`H`$ if and only if one of the eight odd factors has a zero off the critical line.** If all eight factors satisfy their critical-line assertion, every sufficiently large starting height misses at least one of the two moves at this deadline. In the off-line case successes have positive lower logarithmic density for every fixed $`C\gt2`$. Every such factor works at *all* sufficiently late starts exactly when the supremum of the zero real parts is unattained. [Short-reversal criterion](PROOF.md#a-short-reversal-detects-an-off-line-zero)

**Both moves do fit before $`2.011H`$ for a set of starting heights of positive lower logarithmic density, unconditionally.** Those successful log heights contain intervals of fixed positive length with bounded gaps. In the critical-line case two evaluations of the gamma kernel locate the opposite signs at $`\pm0.00273`$ in the limiting profile. [Close-pair proof](PROOF.md#both-signs-in-a-window-ending-by-2011h)

## Verification

Run `bash verify-release.sh` with Python 3.10 or newer and python-flint 0.9.0 installed. The script checks the manifest, compares rational and interval calculations, tests the signed moment construction, evaluates both variation bounds on 13,946 archived ordinate brackets, checks a second close pair by two gamma-kernel evaluations, and replays five local Lean deductions. On a platform with the pinned Lean and a built Mathlib checkout, pass its path as the sole argument. The density theorem, compactly supported interpolation and uniform prime-score estimates are written proofs. [Numerical and formal transcript](evidence/check.txt) · [Prior work](PRIOR-WORK.md) · [Source hashes](SHA256SUMS)

```sh
python3 -m pip install -r requirements.txt
bash verify-release.sh
```
