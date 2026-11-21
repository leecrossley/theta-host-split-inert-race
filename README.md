# A split–inert prime race on a theta-lift host

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24410872.svg)](https://doi.org/10.5281/zenodo.24410872)

Set $K=\mathbf Q(\zeta_{16}+\zeta_{16}^{-1})$, $E=K(\sqrt{-571},\sqrt{-599})$ and $F=K(\sqrt{571\cdot599})$. These are the degree-16 CM field and its degree-8 real subfield in Alpöge’s explicit example. The primes being counted are prime ideals of $F$. Each gets a sign from the first step of the theta parameter recursion over $E/F$. [Host construction](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Formalism/LocalRingCardinality/SpanPowerMembership/ExampleLayerExistsGenerators.lean#L207)

**Both sides win by a fixed margin at every fixed multiplicative width.** For $b>1$, give a prime of norm $n$ the weight $\log n$ times

$$
W_b(t)=
\begin{cases}
b(1-1/t) & 1 < t\le b,\\
b^2/t-1 & b < t < b^2,\\
0 & \text{otherwise},
\end{cases}
\qquad t=n/y .
$$

Let $A_b(y)$ be the weight of the inert primes minus the weight of the split primes. Put $B_b=2\sqrt b(\sqrt b-1)^2$. Then

$$
\liminf_{y\to\infty}\frac{A_b(y)}{B_b\sqrt y}< -5.1027,
\qquad
\limsup_{y\to\infty}\frac{A_b(y)}{B_b\sqrt y}> 7.1027 .
$$

The inequalities are unconditional. Each holds for every fixed real $b>1$, even arbitrarily close to 1. The constants are independent of $b$. Each direction occurs at arbitrarily large locations, with the locations allowed to depend on the width. Dropping any fixed finite set of primes gives the same inequalities.

For $b=2$, a direct count gives $7.512<A_2(50)/(B_2\sqrt{50})<7.513$ and $-38.475<A_2(5000)/(B_2\sqrt{5000})<-38.474$. These are finite examples of both leads. [Direct prime-ideal count](independent.py)

## Where the trace enters

At each prime $w$ of $E$ above a prime $v$ of $F$, subtract the twisted trace of the preceding parameter row from the trace of the next row. Divide this difference by the trace of the added untwisted segment. Multiply those normalized jumps over all $w$ above $v$. Away from the level, the answer is $+1$ when $v$ splits and $-1$ when $v$ is inert. A split prime has two factors in this product. An inert prime has one. The race counts each $v$ once with its norm in $F$.

The Lean proof constructs character data satisfying the upstream admissibility record on this particular host and proves the product identity for every input parameter family. The arithmetic bounds follow from a written Mellin-transform argument and numerical certificates. They apply to the resulting trace observable and to the ordinary split–inert count. [Concrete admissible datum and trace theorem](HostAdmissible.lean) · [Trace extraction](FiberTrace.lean) · [Full proof](PROOF.md)

## Proof outline

1. The eight Dirichlet characters odd on complex conjugation give the factors of $\zeta_E/\zeta_F$. The squares of prime ideals contribute the positive term $B_b$ after normalization by $\sqrt y$.
2. Four zeros have ordinates close to $0.122949$, $0.279415$, $0.319523$ and $0.378007$. Two computations enclose a zero in each interval and exclude zeros at all 36 mixed frequencies up to sign. A separate interval cover checks the entire real segment from $1/2$ to 1 for all eight factors.
3. Two nonnegative products of cosines select the four residues with opposite signs. Their combined contribution exceeds $6.1027 B_b$ at every width. Landau’s theorem lets the argument use these finitely many residues without assuming GRH or independence of the remaining zeros. [Derivation and certificate bounds](PROOF.md)

## Prior work and scope

The first-step trace identity, the character restriction over a quadratic fibre, the weight and the general use of finite sets of zeros in oscillation theorems all have precedents. The contribution claimed here is the field-specific theorem, its realization by admissible theta parameter data and the four-ordinate certificate that gives both constants for all widths. The comparison through 21 November 2026 found no earlier derivation of those bounds for this race. [Source-by-source comparison](PRIOR-WORK.md)

## Verification

The numerical checks need Python 3 and python-flint 0.9.0. The two spectrum computations use different formulas. Both use Arb ball arithmetic.

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python verify.py
python independent.py
python real_axis.py
```

For a fresh replay of the numerical and Lean checks on x86-64 Linux, install the requirements above and run the following with Git, curl, tar, zstd and a C toolchain available. The script downloads the pinned Lean, Mathlib, original sources and trace module. The fresh replay builds about 2,800 upstream modules and needs several gigabytes for the Mathlib cache and compiled dependencies. The four-path mode trusts the supplied compiled cache and uses modification times to decide which upstream modules to rebuild. Use an empty build directory to recompile the upstream sources. [Replay instructions](PROOF.md) [Verification output](evidence/check.txt) [File hashes](SHA256SUMS)

```sh
bash verify-release.sh
```
