import Mathlib

open scoped BigOperators

namespace PairedDensity

/-- The kernel proof applies to either sign after multiplying by that sign. -/
theorem translated_margin {head motion rem err value : ℝ}
    (h : head - motion - rem - err ≤ value)
    (hh : 164/100 - 43/1000 - 24/1000 < head)
    (hm : motion ≤ 1/10) (hr : rem ≤ 144/100) (he : err ≤ 1/1000) :
    1/40 < value := by linarith

/-- The coordinatewise error bounds add even when the frequencies are dependent. -/
theorem sum_budget {ι : Type*} (s : Finset ι) (change budget : ι → ℝ)
    (h : ∀ j ∈ s, |change j| ≤ budget j) :
    |∑ j ∈ s, change j| ≤ ∑ j ∈ s, budget j := by
  exact (Finset.abs_sum_le_sum_abs _ _).trans (Finset.sum_le_sum h)

/-- Collision pairs have total mass at least the square of the mass divided by the box count. -/
theorem collision_bound {ι : Type*} (s : Finset ι) (mass : ι → ℝ) :
    (∑ i ∈ s, mass i)^2 ≤ (s.card : ℝ) * ∑ i ∈ s, (mass i)^2 :=
  sq_sum_le_card_mul_sum_sq

/-- Convert the ordered-pair estimate into a finite-interval measure bound. -/
theorem finite_fraction {T M good pairs : ℝ} (hT : 0 < T) (hM : 0 < M)
    (lower : T^2 ≤ M*pairs) (upper : pairs ≤ 2*T*good) : T/(2*M) ≤ good := by
  apply (div_le_iff₀ (by positivity : (0:ℝ) < 2*M)).2
  apply le_of_mul_le_mul_left (a := T) _ hT
  calc
    T*T = T^2 := by ring
    _ ≤ M*pairs := lower
    _ ≤ M*(2*T*good) := mul_le_mul_of_nonneg_left upper hM.le
    _ = T*(good*(2*M)) := by ring

/-- Exact margins for the coordinate allocation and the last included marginal reduction. -/
theorem partition_budget :
    (788089118617009:ℝ)/7881299347898368 < 1/10 ∧
    (1/10:ℝ) ≤ 788180326378929/7881299347898368 := by norm_num

end PairedDensity

#print axioms PairedDensity.translated_margin
#print axioms PairedDensity.sum_budget
#print axioms PairedDensity.collision_bound
#print axioms PairedDensity.finite_fraction
#print axioms PairedDensity.partition_budget
