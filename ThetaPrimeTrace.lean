import DyadicArithmetic
import Development.AbelianVarietyTheory.Proofs.DecoratedLatticeProofs.ThetaLiftChainHeckeParameterTwistNextTargetChar

open Formalism.BasicNotions GlobalTheory.ChainOnto GlobalTheory.ChainOnto.ParClosed
open NumberField IsDedekindDomain
attribute [local instance] Formalism.BasicNotions.Host.cmStar

noncomputable section
namespace ThetaTransfer

def firstTrace (b : Multiset ℂˣ) : ℂ := (b.map (fun z : ℂˣ => (z : ℂ))).sum

private theorem firstTrace_add (b b' : Multiset ℂˣ) :
    firstTrace (b + b') = firstTrace b + firstTrace b' := by
  simp [firstTrace]

private theorem firstTrace_twist (b : Multiset ℂˣ) (z : ℂˣ) :
    firstTrace (b.map (z * ·)) = (z : ℂ) * firstTrace b := by
  induction b using Multiset.induction_on with
  | empty => simp [firstTrace]
  | cons a b ih =>
    simpa [firstTrace, mul_add] using congrArg (fun v : ℂ => (z : ℂ) * (a : ℂ) + v) ih

variable {D : D7}

private theorem firstTrace_segment (w : HeightOneSpectrum (𝓞 D.host.E)) (z : ℂˣ) (t : ℕ) :
    firstTrace (GlobalTheory.LFun.seg w z t) = (z : ℂ) * firstTrace (GlobalTheory.LFun.seg w 1 t) := by
  rw [← firstTrace_twist]
  congr 1
  simp [GlobalTheory.LFun.seg, Multiset.map_map]

/-- Every exponent in the untwisted segment is real and its norm base is positive. -/
private theorem segment_trace_ne_zero (w : HeightOneSpectrum (𝓞 D.host.E)) :
    firstTrace (GlobalTheory.LFun.seg w 1 D.R.t) ≠ 0 := by
  let q : ℝ := Ideal.absNorm w.asIdeal
  let r (j : ℕ) : ℝ := q ^ (((D.R.t : ℝ) - 1) / 2 - j)
  have hq : 0 < q := by
    dsimp [q]
    exact_mod_cast (GlobalTheory.LFun.one_lt_absNorm w).bot_lt
  have term (j : ℕ) : GlobalTheory.LFun.npow w (((D.R.t : ℂ) - 1) / 2 - j) = (r j : ℂ) := by
    dsimp [r, q, GlobalTheory.LFun.npow]
    rw [Complex.ofReal_cpow (Nat.cast_nonneg _)]
    push_cast
    rfl
  have form : firstTrace (GlobalTheory.LFun.seg w 1 D.R.t) =
      ((∑ j ∈ Finset.range D.R.t, r j : ℝ) : ℂ) := by
    rw [Complex.ofReal_sum]
    simp [firstTrace, GlobalTheory.LFun.seg, Multiset.map_map, Function.comp_def, Finset.sum, term]
  have hpos (j : ℕ) : 0 < r j := Real.rpow_pos_of_pos hq _
  have ht : 0 < D.R.t := (Nat.zero_le _).trans_lt D.R.three_m_lt_t
  have pos : 0 < ∑ j ∈ Finset.range D.R.t, r j :=
    Finset.sum_pos' (fun j _ => (hpos j).le) ⟨0, Finset.mem_range.mpr ht, hpos 0⟩
  rw [form]
  exact_mod_cast pos.ne'

/-- The previous stage and every older segment cancel in this one-step difference. -/
theorem theta_step_trace (c : Com D) (B : GlobalTheory.LFun.Fam D.host.E) (k : ℕ)
    (w : HeightOneSpectrum (𝓞 D.host.E)) :
    firstTrace (chainPar c B (k+1) w) -
      ((vf (c.chS (k+1) * (c.chT (k+1))⁻¹) w : ℂˣ) : ℂ) * firstTrace (chainPar c B k w) =
      ((vf (c.chS (k+1)) w : ℂˣ) : ℂ) * firstTrace (GlobalTheory.LFun.seg w 1 D.R.t) := by
  rw [chainPar_succ, GlobalTheory.ChainOnto.ParClosed.stepPar_eq]
  change firstTrace (((chainPar c B k w).map (_ * ·)) + GlobalTheory.LFun.seg w _ _) - _ = _
  rw [firstTrace_add, firstTrace_twist, firstTrace_segment]
  ring

/-- A finite prime-ideal sum of one character recovered from consecutive tower rows. -/
theorem theta_prime_trace_extract (c : Com D) (B : GlobalTheory.LFun.Fam D.host.E) (k : ℕ)
    (places : Finset (HeightOneSpectrum (𝓞 D.host.E))) (weight : HeightOneSpectrum (𝓞 D.host.E) → ℂ) :
    (∑ w ∈ places, weight w * ((firstTrace (chainPar c B (k+1) w) -
      ((vf (c.chS (k+1) * (c.chT (k+1))⁻¹) w : ℂˣ) : ℂ) * firstTrace (chainPar c B k w)) /
      firstTrace (GlobalTheory.LFun.seg w 1 D.R.t))) =
      ∑ w ∈ places, weight w * ((vf (c.chS (k+1)) w : ℂˣ) : ℂ) := by
  apply Finset.sum_congr rfl
  intro w _
  rw [theta_step_trace, mul_div_cancel_right₀ _ (segment_trace_ne_zero w)]

/-- The finite set may be any collection of good prime ideals. -/
def primeIdealWeight {K : Type} [Field K] [NumberField K] (y : ℝ)
    (w : HeightOneSpectrum (𝓞 K)) : ℂ :=
  (Real.log (Ideal.absNorm w.asIdeal) * tent (y / (Ideal.absNorm w.asIdeal : ℝ)) : ℝ)

theorem theta_band_prime_trace (c : Com D) (B : GlobalTheory.LFun.Fam D.host.E) (k : ℕ)
    (places : Finset (HeightOneSpectrum (𝓞 D.host.E))) (y : ℝ) :
    (∑ w ∈ places, primeIdealWeight y w * ((firstTrace (chainPar c B (k+1) w) -
      ((vf (c.chS (k+1) * (c.chT (k+1))⁻¹) w : ℂˣ) : ℂ) * firstTrace (chainPar c B k w)) /
      firstTrace (GlobalTheory.LFun.seg w 1 D.R.t))) =
      ∑ w ∈ places, primeIdealWeight y w * ((vf (c.chS (k+1)) w : ℂˣ) : ℂ) :=
  theta_prime_trace_extract c B k places (primeIdealWeight y)

#print axioms ThetaTransfer.theta_step_trace
#print axioms ThetaTransfer.theta_prime_trace_extract
#print axioms ThetaTransfer.theta_band_prime_trace

end ThetaTransfer
