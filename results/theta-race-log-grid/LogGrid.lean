import Mathlib

namespace LogGrid

/-- Exact reciprocal mass used when summing the absolute-sine Fourier coefficients. -/
theorem reciprocal_step {k : ℝ} (hk : 1 ≤ k) :
    1/(4*k^2-1) = (1/(2*k-1)-1/(2*k+1))/2 := by
  have ha : 0 < 2*k-1 := by linarith
  have hb : 0 < 2*k+1 := by linarith
  have split (x y : ℝ) (hx : x ≠ 0) (hy : y ≠ 0) :
      (1/x-1/y)/2 = (y-x)/(2*(x*y)) := by
    field_simp
  rw [split _ _ ha.ne' hb.ne']
  have diff : 2*k+1-(2*k-1) = 2 := by ring
  have prod : (2*k-1)*(2*k+1) = 4*k^2-1 := by ring
  rw [diff,prod,mul_comm (2:ℝ),div_mul_eq_div_div_swap]
  norm_num

/-- The uniform profile modulus preserves the strict prime-only threshold. -/
theorem primary_margin :
    (4917/1000:ℝ)-285/100-(144/100)/(48*(3/32))^2-7/1000-141/10000 > 19735/10000 := by
  norm_num

theorem second_margin :
    (4917/1000:ℝ)-285/100-(152/100)/(48*(3/32))^2-7/1000-85/10000 > 19735/10000 := by
  norm_num

/-- A buffer around the sampling point protects the whole interval of starts. -/
theorem buffered {left right x y d : ℝ}
    (lo : left+d < x) (hi : x < right-d) (perturb : |y-x| ≤ d) :
    left < y ∧ y < right := by
  obtain ⟨hm,hp⟩ := abs_le.mp perturb
  constructor <;> linarith

/-- Every interval of one mesh width contains a point of any translated lattice. -/
theorem grid_hits (offset left step : ℝ) (hs : 0 < step) :
    ∃ n : ℤ, left ≤ offset + n*step ∧ offset+n*step < left+step := by
  let q := (left-offset)/step
  refine ⟨⌈q⌉, ?_, ?_⟩
  · have h := (div_le_iff₀ hs).mp (Int.le_ceil q : (left-offset)/step ≤ (⌈q⌉ : ℝ))
    linarith
  · have h := mul_lt_mul_of_pos_right (Int.ceil_lt_add_one q) hs
    change (⌈q⌉ : ℝ)*step < ((left-offset)/step+1)*step at h
    rw [add_mul,div_mul_cancel₀ _ hs.ne',one_mul] at h
    linarith

end LogGrid

#print axioms LogGrid.reciprocal_step
#print axioms LogGrid.primary_margin
#print axioms LogGrid.second_margin
#print axioms LogGrid.buffered
#print axioms LogGrid.grid_hits
