import Mathlib

namespace FourthRoot

/-- Symmetric sampling cancels a level that is affine in the height
before square-root normalization. The scale ratio is positive. -/
theorem affine_cell (A B y v : ℝ) (hv : 0 < v) :
    (A+B*(y*v)) + v*(A+B*(y/v)) - (v+1)*(A+B*y) = 0 := by
  field_simp
  ring

/-- The exceptional rational kernel is exactly the two shifted weights
after the square-root factors have been cleared. -/
theorem corrected_pair (l h m c : ℝ) (hc : 0 < c) :
    2*c*((l/c+c*h)/2-(c+1/c)*m/2) = l+c^2*h-(c^2+1)*m := by
  field_simp

/-- Independent mean-zero box variances add. -/
theorem box_variance (s n : ℝ) (hn : 0 < n) :
    n*((s/n)^2/3) = s^2/(3*n) := by
  field_simp

/-- The chosen cutoff places each high-frequency box factor below a quarter. -/
theorem box_cutoff {n R T t : ℝ} (hn : 0 < n) (hR : 0 < R)
    (range : T ≤ t) (packed : 8*n ≤ T*R) : 2*n/(t*R) ≤ 1/4 := by
  have htR : 0 < t*R := by nlinarith [mul_le_mul_of_nonneg_right range hR.le]
  apply (div_le_iff₀ htR).2
  nlinarith [mul_le_mul_of_nonneg_right range hR.le]

/-- Prime-side approximation, omitted zeros and arithmetic terms combine
without any sign or independence condition. -/
theorem reconstruction {f filtered spectrum : ℝ}
    (physical : |f-filtered| ≤ 11/100000000)
    (analytic : |filtered-spectrum| ≤ 79/100000 + 1/10000000000) :
    |f-spectrum| < 8/10000 := by
  have h : |f-spectrum| ≤ |f-filtered| + |filtered-spectrum| :=
    abs_sub_le f filtered spectrum
  linarith

/-- A smaller finite exclusion transfers a previously proved strict margin. -/
theorem signed_witness {head value cut error : ℝ} (close : |value-head| ≤ error)
    (margin : cut+error < head) : cut < value := by
  linarith [(abs_le.mp close).1]

/-- The six-box inclusion-exclusion polynomial at its centre. -/
theorem quintic_normalization :
    ((6:ℝ)^5-6*4^5+15*2^5)/2112 = 1 ∧
    20*((6:ℝ)^3-6*4^3+15*2^3)/2112 = -(5/11) := by
  constructor <;> norm_num

/-- At the moment scale used by the polynomial formula, the variance
part of the prime-side error is already below seven parts in a hundred million. -/
theorem quintic_variance_budget :
    (10:ℝ)^5*(5/11)/(5000000/6)^2 < 7/100000000 := by
  norm_num

/-- A uniform sub-quarter gap follows directly from the density exponent. -/
theorem density_gap {s b d : ℝ} (hs : s ≤ b) (hb : b ≤ 1)
    (hd0 : 0 ≤ d) (hd1 : d ≤ 1) (density : d ≤ (12/5)*(1-s)) :
    (s-1/2)/(2-d) ≤ 1/4 - (1-b)/5 := by
  apply (div_le_iff₀ (by linarith : 0 < 2-d)).2
  nlinarith [mul_nonneg (sub_nonneg.mpr hb) hd0]

/-- Numerical saving left by the zero-free region after the cutoff is lowered. -/
theorem density_exponential_budget :
    (1:ℝ)/1000 < 1/600-1/2000 ∧ 1-(5/2)*(1/3:ℝ) = 1/6 := by
  constructor <;> norm_num

end FourthRoot

#print axioms FourthRoot.affine_cell
#print axioms FourthRoot.corrected_pair
#print axioms FourthRoot.box_variance
#print axioms FourthRoot.box_cutoff
#print axioms FourthRoot.reconstruction
#print axioms FourthRoot.signed_witness
#print axioms FourthRoot.quintic_normalization
#print axioms FourthRoot.quintic_variance_budget
#print axioms FourthRoot.density_gap
#print axioms FourthRoot.density_exponential_budget
