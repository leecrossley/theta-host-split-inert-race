import Mathlib

namespace WidthCusp

/-- The rescaled edge of a moving compact window is a triangular profile. -/
theorem triangle (x : ℝ) :
    max (x+1) 0 + max (x-1) 0 - 2*max x 0 = max (1-|x|) 0 := by
  rcases le_or_gt x (-1) with h | h
  · rw [abs_of_nonpos (by linarith), max_eq_right (by linarith : x+1 ≤ 0),
      max_eq_right (by linarith : x-1 ≤ 0), max_eq_right (by linarith : x ≤ 0),
      max_eq_right (by linarith : 1 - -x ≤ 0)]
    ring
  · rcases le_or_gt x 0 with k | k
    · rw [abs_of_nonpos k, max_eq_left (by linarith : 0 ≤ x+1),
        max_eq_right (by linarith : x-1 ≤ 0), max_eq_right k,
        max_eq_left (by linarith : 0 ≤ 1 - -x)]
      ring
    · rcases le_or_gt x 1 with l | l
      · rw [abs_of_nonneg k.le, max_eq_left (by linarith : 0 ≤ x+1),
          max_eq_right (by linarith : x-1 ≤ 0), max_eq_left k.le,
          max_eq_left (by linarith : 0 ≤ 1-x)]
        ring
      · rw [abs_of_nonneg k.le, max_eq_left (by linarith : 0 ≤ x+1),
          max_eq_left (by linarith : 0 ≤ x-1), max_eq_left k.le,
          max_eq_right (by linarith : 1-x ≤ 0)]
        ring

/-- Exact cancellation of the interior exponential solution. -/
theorem hyperbolic_stencil (x h : ℝ) :
    Real.sinh (x+h)+Real.sinh (x-h)=2*Real.cosh h*Real.sinh x := by
  rw [Real.sinh_add, Real.sinh_sub]
  ring

theorem extremal_weight (h t : ℝ) : 0 ≤ Real.cosh h - Real.cos t := by
  linarith [Real.one_le_cosh h, Real.cos_le_one t]

theorem subtract_background (q g c a v : ℝ) :
    (2*(c*(1+g-2*a)-(1+v-q))) - 2*(c-1)-2*(c*g-v) = 2*(q-2*c*a) := by
  ring

theorem single_power (r c logp tr sh : ℝ)
    (hr : r ≠ 0) (hc : c ≠ 0) :
    4*r/c * (logp/r*tr)*sh = 4*logp*sh/c*tr := by
  field_simp

end WidthCusp

#print axioms WidthCusp.triangle
#print axioms WidthCusp.hyperbolic_stencil
#print axioms WidthCusp.extremal_weight
#print axioms WidthCusp.subtract_background
#print axioms WidthCusp.single_power
