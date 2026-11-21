import Mathlib

open scoped BigOperators

namespace AveragingWidths

theorem kernel_floor (u t : ℝ) (hl : 1/2 ≤ u) (hu : u ≤ 2) (ht : t^2 ≤ 1) :
    2/5 ≤ u/(u^2+t^2) := by
  have pos : 0 < u := by linarith
  have denom : 0 < u^2+t^2 := by positivity
  apply (le_div_iff₀ denom).2
  nlinarith only [mul_nonneg (sub_nonneg.mpr hl) (sub_nonneg.mpr hu),ht]

theorem positive_average {n : ℕ} (w g : Fin n → ℝ) (r : ℝ)
    (hw : ∀ j, 0 ≤ w j) (mass : ∑ j, w j = 1) (lower : ∀ j, r ≤ g j) :
    r ≤ ∑ j, w j * g j := by
  calc
    r = (∑ j, w j) * r := by rw [mass, one_mul]
    _ = ∑ j, w j * r := Finset.sum_mul _ _ _
    _ ≤ ∑ j, w j * g j := Finset.sum_le_sum fun j _ =>
      mul_le_mul_of_nonneg_left (lower j) (hw j)

theorem dual_separation {n : ℕ} (d m : Fin n → ℝ) (eps lower : ℝ)
    (large : lower < ∑ j, d j * m j) (budget : (∑ j, |d j|) * eps ≤ lower) :
    ∃ j, eps < |m j| := by
  by_contra! short
  have bound : ∀ j, d j * m j ≤ |d j| * eps := by
    intro j
    calc
      d j * m j ≤ |d j * m j| := le_abs_self _
      _ = |d j| * |m j| := abs_mul _ _
      _ ≤ |d j| * eps := mul_le_mul_of_nonneg_left (short j) (abs_nonneg _)
  have small : ∑ j, d j * m j ≤ (∑ j, |d j|) * eps := by
    rw [Finset.sum_mul]
    exact Finset.sum_le_sum fun j _ => bound j
  exact (not_lt_of_ge (small.trans budget)) large

theorem inverse_margin (w approx defect margin : ℝ)
    (close : |w-approx| ≤ defect) (positive : margin+defect < approx) :
    margin < w := by
  have := (abs_le.mp close).1
  linarith

theorem shared_saving (a b common : ℝ)
    (ha : 88436/1000 < a) (hb : 47698/1000 < b)
    (hc : common ≤ 103949/1000) : 32185/1000 < a+b-common := by
  linarith

end AveragingWidths

#print axioms AveragingWidths.kernel_floor
#print axioms AveragingWidths.positive_average
#print axioms AveragingWidths.dual_separation
#print axioms AveragingWidths.inverse_margin
#print axioms AveragingWidths.shared_saving
