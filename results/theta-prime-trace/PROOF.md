# Derivation and context: recovering a character sum from two theta-lift traces

This is the detail behind the [summary](README.md). The Lean files are the checked version.

## Derivation

At an unramified prime ideal $v$, let $T_k(v)$ be the sum of the local Hecke parameters after step $k$. Write $q=Nv$, $z_k(v)=(\chi_S\chi_T^{-1})(v)$, and $A_t(q)=\sum_{j=0}^{t-1}q^{(t-1)/2-j}$. Here the two characters are those of step $k+1$. Apply the first trace to `stepPar_eq`. A twist multiplies the old trace by $z_k(v)$ and the new segment contributes $\chi_S(v)A_t(q)$. Thus

$$
  T_{k+1}(v)-z_k(v)T_k(v)=\chi_S(v)A_t(q).
$$

For $t>0$, $A_t(q)>0$. Divide and sum over any finite set of good prime ideals with weight $(\log Nv)W(Nv/y)$. Two adjacent rows recover the weighted character sum exactly, with the character values retaining their signs and phases. This is an identity for supplied rows. A bound for the resulting sum requires bounds for their difference. [Segment definition](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AbelianVarietyTheory/Definitions/Foundations/NormEulerPolySubOneLe.lean#L120), [Our application](ThetaPrimeTrace.lean)

The weight $W$ in this specialization is defined by

$$
 W(t)=
 \begin{cases}
 2(1-1/t) & 1 < t\le2,\\
 4/t-1 & 2 < t < 4,\\
 0 & \text{otherwise},
 \end{cases}
 \qquad
 P(y)=\sum_{y < p < 4y}(\log p)W(p/y).
$$

## From the upstream argument to the declarations

The source archive was matched to SHA-256 `9631b4f020199e580ba86178cb56c6e9b09f690c2369a80dc4334aa93f49c7fc`. The source links below point to mirror revision `40e56bab026552e089d7afcfeff08d2defb9c005`. Alpöge describes repeated lifts as a substitute for multiplicity. The source expresses the idea through a choice of scale, an iteration, a formula for its parameters and analytic estimates. [Author's note](https://alpo.ge/remark.pdf), [Independent verification of the original release](https://math.mit.edu/~drew/grh/)

| Step | What the declaration does | Upstream source |
|---|---|---|
| Choose the length | `GlobalTheory.Mem.Num.mAt` takes the even ceiling of $\sqrt{x\delta_1}$, where $\delta_1=x^{-\varepsilon/8}$. Thus the length is within two of $x^{1/2-\varepsilon/16}$. `RankVal` records bounds for admissible data. | [mAt](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AbelianVarietyTheory/Definitions/AutomorphicMembership/EvenCeil.lean#L21), [delta1](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Formalism/LocalRingCardinality/DataParams%E2%82%82.lean#L45), [RankVal](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AbelianVarietyTheory/Definitions/GlobalTwisting/PairRowIffQdRows.lean#L28) |
| Perform a lift | The `step` constructor of `IsIterLift` integrates the old function against an admissible theta kernel. Its local Hecke parameters obey `GlobalTheory.ChainOnto.chainPar_succ`. | [Integral constructor](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AbelianVarietyTheory/BlockTheory/TwistedDiscreteSeries/IteratedLiftChains.lean#L91), [Parameter recursion](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AbelianVarietyTheory/Definitions/LatticeFrames/ChainParSucc.lean#L229) |
| Count what accumulates | `ParClosed.stepPar_eq` identifies one twist and one new segment. `ParClosed.chainPar_closed` proves by induction that after $k$ steps all $k$ segments are present, indexed by their characters. | [One step and induction](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AbelianVarietyTheory/Proofs/DecoratedLatticeProofs/ThetaLiftChainHeckeParameterTwistNextTargetChar.lean#L47) |
| Relate the parameters to actual lifts | `lFull_holds` rewrites an Euler product using `thetaLiftChainHeckeParameter_twist_next_target_char`, producing products over the stages and the positions in a segment. `scalars_of_given_all` uses this factorization along with `thetaLiftChainStep_exact_innerProduct` to obtain the step scalars. Both the parameters and the scalar estimates matter. | [Euler factorization](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AbelianVarietyTheory/Proofs/DecoratedLatticeProofs/TwistedThetaChainLValueQuotientFactorization.lean#L51), [Assembly of the scalars](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AbelianVarietyTheory/Proofs/DecoratedLatticeProofs/ScalarsOfGivenData.lean#L53) |
| Supply a prime-sum estimate | `every_character_of_literature_of_vWin` passes the theta integral and local matching theorems to `every_character_of_theta_integrals_local_matching`. Its proof uses `twisted_square_theta_integrals_eq` in the evaluation grid. `zeroFreeAll_of_topG` combines tied families with `level_gm_le_on_family_scales` and `formB`. The comparison supplies a one-sided bound $G_m(x)\le x^{1/2+(e_D+e_L+7)\varepsilon}$ on family windows. | [Theta and local inputs](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/ThetaCorrespondence/ThetaPairLifts/ThirdPlaceSeries/EveryCharacterOfLiteratureOfVWin.lean#L810), [Evaluation grid](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/ThetaCorrespondence/ThetaPairLifts/ThirdPlaceSeries/SiegelEisensteinSecondTermHypotheses.lean#L52), [Convergence of the two branches](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/CharacterBoxes/ZeroFreeAllOfTopG.lean#L26), [Inequality type](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Formalism/SharpRestriction/RamificationFrobenius/TotientGaloisClasses/TopWin.lean#L26) |
| Exclude zeros | The proof body of `UpperTheory.Comparison₂.formB` invokes `Deduction.zeroFree_of_rpow_window`. This calls `lemmaL`, which obtains `poleData` and invokes `landau_window`. At the top, `riemannHypothesis` in `Main.lean` uses `rhOfAll`, which applies `Formalism.VW.rh_of_every` to `allChi`. | [formB](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/PrimitiveCosets/Imports.lean#L8), [Power-window reduction](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/DedekindLFunctions/LandauTheorems/ZeroFreeOfRpowWindow.lean#L11), [lemmaL](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/DedekindLFunctions/LandauTheorems/lemmaL.lean#L13), [landau_window](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/DedekindLFunctions/LandauTheorems/landauWindow.lean#L83), [Top-level reduction](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Supplements/Kernels/MainTheorems.lean#L7) |

## A character to test

A concrete geometric target is $E=\mathbb Q(\sqrt{-3})$ with the quadratic character $\chi_4(Nv)$ at prime ideals $v\nmid6$. Prime ideals of norm 7 contribute with a minus sign and those of norm 13 with a plus sign. Both rational primes split, as witnessed by roots 2 and 3 of $X^2+X+1$ modulo 7 and 13 respectively. Away from 2 and 3, the Euler product for this character factors as $L(s,\chi_4)L(s,\chi_4\chi_{-3})$, where $\chi_{-3}$ is the quadratic character of $E$. At a split prime the two factors agree. At an inert prime their product has one factor in $p^{-2s}$, as required by the norm. This gives an independent analytic benchmark. The first test is whether admissible theta data realize these character values. The second is whether bounding the difference of the two rows improves the explicit bound for the weighted difference of the two residue classes obtained from that Euler product.

Later results settle the direct route. The Eisenstein field cannot be a tower host. Every upstream host has degree divisible by four over $\mathbb Q$, and `EisensteinPrimeRace.no_quadratic_host` proves that no quadratic field is rationally isomorphic to a host. The Euler product gives the benchmark directly. Widen the dyadic tent to $W_b$, supported on $(y,b^2y)$ with its peak at $by$. Above a critical $b$ between $4.49334142447806$ and $4.49334142447807$, the norm-3 class eventually leads if and only if GRH holds for the characters of discriminants $-4$ and $12$. The trace identity still applies to the prime ideals of an actual host. On the 571/599 host it gives each good prime ideal of the degree-8 real subfield its split or inert sign, and that race oscillates at every fixed width. Transferring the Eisenstein race to a host would need a separate theorem comparing their weighted prime ideals. [Degree obstruction](https://github.com/leecrossley/eisenstein-prime-window-transition/blob/c1e11263b1cf6e8de187742f246b32992455658f/PROOF.md#what-the-lean-proof-covers), [Eisenstein cutoff](https://github.com/leecrossley/eisenstein-prime-window-transition/blob/c1e11263b1cf6e8de187742f246b32992455658f/README.md), [Race on a host's own primes](https://github.com/leecrossley/theta-host-split-inert-race/blob/c60c7344e19f5c9ada820a0fb69473dddfc8dcee/PROOF.md#field-and-local-observable)

## Reproducing the Lean check

The recorded run used Lean 4.33.1 on arm64 macOS. It freshly elaborated this repository's modules, used compiled dependencies from the pinned sources, replayed each local module with `leanchecker`, and inspected the named declarations' axioms. The receipt records the module counts and axiom lists. Compiler commit `819816b2e0a3bf405af45ae5c7af2491d8f5bee6` and Mathlib commit `0df444a360eaa60ab8c11dca51a86af692955474` fix the environment. [Lean release](https://github.com/leanprover/lean4/releases/tag/v4.33.1), [Mathlib revision](https://github.com/leanprover-community/mathlib4/tree/0df444a360eaa60ab8c11dca51a86af692955474)

The source archive for the transitive upstream modules is pinned by SHA-256 `9631b4f020199e580ba86178cb56c6e9b09f690c2369a80dc4334aa93f49c7fc`. [Source release](https://github.com/CocoSgt/riemann-hypothesis-lean/releases/tag/v2026-11-20), [Readable source mirror](https://github.com/Deicyde/riemann-hypothesis/tree/40e56bab026552e089d7afcfeff08d2defb9c005)

For a fresh build on x86-64 Linux with Git, curl, tar, zstd and Python 3.9 or later, run

```sh
bash verify-linux.sh
```

With the pinned toolchain and Mathlib cache already installed, use

```sh
cd /path/to/pinned/mathlib4
lake env python3 /path/to/repository/verify.py /path/to/upstream/lean/src /path/to/build
```

The build directory can begin empty. Reusing a dependency cache requires the same sources and compiler. A successful run ends with the kernel and axiom summary recorded in [evidence/check.txt](evidence/check.txt).

The Linux entry point checks the SHA-256 manifest before downloading the pinned dependencies and compares the resulting kernel and axiom summary with the recorded run.
