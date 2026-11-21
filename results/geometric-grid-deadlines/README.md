# Both prime-race moves on geometric grids up to ratio four

On the 547/727 field, **every block of 223 consecutive dyadic starting heights** above $`e^{100}`$ contains a gain and a loss of more than $`1.9735\log 2`$. On 571/599 the bound is 319. The first height in the block can be any real number. For example, on the first field check the doublings starting at $`Y,2Y,\ldots,2^{222}Y`$. Both signs must occur. These are prime-ideal counts at width and dilation two. The proof uses the pinned primitive Dirichlet critical-line theorem (Dirichlet GRH) for the two fields' character factors. [Proof](PROOF.md) · [Arithmetic trace](https://github.com/leecrossley/theta-race-short-shifts/blob/483c31f2df2132ed2165620ea55e20e969228bfb/RaceTrace.lean#L83) · [Independent upstream verification](https://math.mit.edu/~drew/grh/)

| Allowed starts, at any offset | What is guaranteed in each field |
|---|---|
| Ratio 2 | Each direction occupies lower uniform density at least 1/223 of the starting indices for 547/727 and 1/319 for 571/599. |
| Any ratio between 2 and 4, inclusive | Every block of 4,095 consecutive starts above $`e^{100}`$ contains both directions. |
| Any fixed ratio greater than 1 and at most 4 | Each direction occupies lower uniform density at least 1/4,095. After every real $`H\ge1`$ both doublings finish by $`2e^{100}4^{4095}H\lt10^{2510}H`$. |
| Ratios with no upper bound | For every block length there are arbitrarily large power-of-two ratios and late starting offsets with that many successive gains on both fields. The same ratios also permit long losses and the other two choices of joint direction. |

For starts restricted to the actual powers of two, both signs therefore occur at high integer exponents at the densities above. The unrestricted continuous-height deadline $`10^{18}H`$ for 547/727 is sharper when the observer can choose any real starting point. Its guarantee for arbitrary grids covered log steps up to $`10^{-6}`$. [Earlier deadline](https://github.com/leecrossley/theta-race-filtered-deadline/blob/2400dc5ac2a7969d41e1a4bde1eb44dd6c63b7a2/PROOF.md) · [Earlier grid theorem](https://github.com/leecrossley/theta-race-log-grid/blob/dbdfc627384c841b78bd33656564b2cf9cb18890/PROOF.md)

## Verification

Two independent interval calculations cover the entire real interval of grid ratios. A quotient of sines and a separate product of cosines give the sharper dyadic counts. The finite certificate includes both complete covers, the 14,000 inherited zero brackets and a separate sign replay. Six local Lean deductions cover the extraction of a point from a positive average, the grid stride, and the numerical implications. The filter construction and the analytic estimates are written proofs. [Checksums](SHA256SUMS) · [Verification output](evidence/check.txt) · [Prior work](PRIOR-WORK.md)

Use Python 3.11 or later. On Linux x86-64 the replay downloads pinned Lean and Mathlib. With an existing pinned checkout, put Lean 4.33.1 on `PATH` and give the checkout path to the script. The optional last command repeats the earlier 28,000 endpoint signs.

```sh
python3 -m pip install -r requirements.txt
bash verify-release.sh
python3 check_grid_zeros.py
```
