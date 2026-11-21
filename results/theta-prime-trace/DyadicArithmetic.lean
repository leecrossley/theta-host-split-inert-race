import Development.AutomorphicArithmetic.DedekindLFunctions.CharFamily

open scoped ArithmeticFunction.vonMangoldt
open Complex

noncomputable section
namespace ThetaTransfer

/-- Write t = y/n. Positive weights have y < n < 4y. -/
def tent (t : ℝ) : ℝ :=
  max (4 * t - 1) 0 - 3 * max (2 * t - 1) 0 + 2 * max (t - 1) 0

theorem tent_pieces (t : ℝ) :
    tent t = if 1 ≤ t then 0 else if 1 / 2 ≤ t then 2 * (1 - t)
      else if 1 / 4 ≤ t then 4 * t - 1 else 0 := by
  unfold tent
  split_ifs <;> grind [max_def]

theorem tent_nonneg (t : ℝ) : 0 ≤ tent t := by
  rw [tent_pieces]
  split_ifs <;> linarith

theorem tent_le_one (t : ℝ) : tent t ≤ 1 := by
  rw [tent_pieces]
  split_ifs <;> linarith

def primePowerSum (y : ℝ) : ℝ :=
  ∑ n ∈ Finset.range (⌊4 * y⌋₊ + 1), Λ n * tent (y / n)

def primeSum (y : ℝ) : ℝ :=
  ∑ n ∈ (Finset.range (⌊4 * y⌋₊ + 1)).filter Nat.Prime, Real.log n * tent (y / n)

theorem primeSum_le_primePowerSum (y : ℝ) : primeSum y ≤ primePowerSum y := by
  unfold primeSum primePowerSum
  calc
    _ = ∑ n ∈ (Finset.range (⌊4 * y⌋₊ + 1)).filter Nat.Prime, Λ n * tent (y / n) := by
      apply Finset.sum_congr rfl
      intro n hn
      rw [ArithmeticFunction.vonMangoldt_apply_prime (Finset.mem_filter.mp hn).2]
    _ ≤ _ := Finset.sum_le_sum_of_subset_of_nonneg (Finset.filter_subset _ _)
      (fun n _ _ => mul_nonneg (ArithmeticFunction.vonMangoldt_nonneg) (tent_nonneg _))

theorem primePowerGap_le (y : ℝ) :
    primePowerSum y - primeSum y ≤ Chebyshev.psi (4*y) - Chebyshev.theta (4*y) := by
  rw [primePowerSum, primeSum, Nat.range_succ_eq_Icc_zero, Chebyshev.psi_eq_sum_Icc,
    Chebyshev.theta_eq_sum_Icc, Finset.sum_filter, Finset.sum_filter,
    ← Finset.sum_sub_distrib, ← Finset.sum_sub_distrib]
  apply Finset.sum_le_sum
  intro n _
  by_cases h : Nat.Prime n
  · simp [h, ArithmeticFunction.vonMangoldt_apply_prime h]
  · simpa [h] using mul_le_of_le_one_right (ArithmeticFunction.vonMangoldt_nonneg (n := n))
      (tent_le_one (y / n))

private theorem ramp_extend {x z : ℝ} (hxz : x ≤ z) :
    (Deduction.ramp (fun _ => 1) x).re =
      ∑ n ∈ Finset.range (⌊z⌋₊ + 1), Λ n * max (x / n - 1) 0 := by
  simp [Deduction.ramp, Complex.re_sum]
  apply Finset.sum_subset_zero_on_sdiff
  · exact Finset.range_subset_range.mpr (by gcongr)
  · intro i hi
    have hi' := (Finset.mem_sdiff.mp hi).2
    simp only [Finset.mem_range, not_lt] at hi'
    have hi0 : 0 < i := by omega
    have hxi : x < (i : ℝ) :=
      (Nat.floor_lt' (Nat.ne_of_gt hi0)).mp (by omega)
    rw [max_eq_right (by
      have := (div_le_one (Nat.cast_pos.mpr hi0 : 0 < (i : ℝ))).mpr hxi.le
      linarith)]
    simp
  · intros
    rfl

theorem primePowerSum_eq_ramps {y : ℝ} (hy : 0 < y) :
    primePowerSum y = (Deduction.ramp (fun _ => 1) (4*y)).re
      - 3 * (Deduction.ramp (fun _ => 1) (2*y)).re
      + 2 * (Deduction.ramp (fun _ => 1) y).re := by
  rw [ramp_extend (z := 4*y) le_rfl, ramp_extend (x := 2*y) (z := 4*y) (by linarith),
    ramp_extend (x := y) (z := 4*y) (by linarith)]
  simp only [primePowerSum, Finset.mul_sum, ← Finset.sum_sub_distrib, ← Finset.sum_add_distrib]
  apply Finset.sum_congr rfl
  intro n _
  simp only [tent, mul_div_assoc]
  ring

theorem deficit_eq_centered_ramps {y : ℝ} (hy : 0 < y) :
    2 * Real.log 2 * y - primePowerSum y =
      -Deduction.G (4*y) + 3 * Deduction.G (2*y) - 2 * Deduction.G y := by
  rw [primePowerSum_eq_ramps hy]
  have h4 : Real.log 4 = 2 * Real.log 2 := by
    rw [show (4 : ℝ) = 2 ^ 2 by norm_num, Real.log_pow]
    norm_num
  unfold Deduction.G
  rw [Real.log_mul (by norm_num) hy.ne', Real.log_mul (by norm_num) hy.ne', h4]
  ring

end ThetaTransfer
