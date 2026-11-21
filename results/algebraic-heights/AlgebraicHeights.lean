import Mathlib

namespace AlgebraicHeights

/-- A protected logarithm stays outside the span after nonzero scaling. -/
theorem protected_log {K V : Type*} [Field K] [AddCommGroup V] [Module K V]
    (W : Submodule K V) {v : V} (hv : v ∉ W) {c : K} (hc : c ≠ 0) :
    c • v ∉ W := by
  intro h
  have back := W.smul_mem c⁻¹ h
  apply hv
  simpa [smul_smul, hc] using back

/-- Two distinct zeros force an affine piece to vanish identically. -/
theorem affine_zeros {a b u v : ℝ} (huv : u ≠ v)
    (hu : a*u+b = 0) (hv : a*v+b = 0) : a = 0 ∧ b = 0 := by
  have h : a*(u-v) = 0 := by linarith only [hu, hv]
  have ha := (mul_eq_zero.mp h).resolve_right (sub_ne_zero.mpr huv)
  exact ⟨ha, by simpa [ha] using hu⟩

/-- The coefficients for inserting the intermediate knots at b and b^3. -/
theorem square_refinement (b x : ℝ) :
    (x^2-1)*(x^2-b^2) = (x-1)*(x-b)*(b+(b+1)*x+x^2) := by ring

/-- Clearing the positive square-root coordinate yields a quadratic. -/
theorem threshold_cell {u : ℝ} (hu : 0 < u) (a b t : ℝ) :
    a*u+b/u = t ↔ a*u^2-t*u+b = 0 := by
  constructor
  · intro h
    have cleared := congrArg (fun s : ℝ => s*u) h
    field_simp [ne_of_gt hu] at cleared
    nlinarith [cleared]
  · intro h
    apply (mul_left_inj' (ne_of_gt hu)).mp
    field_simp [ne_of_gt hu]
    nlinarith [h]

/-- Sliding a strict winning interval against an allowed later-start interval. -/
theorem eligible {l u h d : ℝ} (hlu : l < u) (hhd : h < d) :
    (∃ y, h < y ∧ y ≤ d ∧ l < y ∧ y < u) ↔ h < u ∧ l < d := by
  constructor
  · rintro ⟨y, hhy, hyd, hly, hyu⟩
    exact ⟨lt_trans hhy hyu, lt_of_lt_of_le hly hyd⟩
  · rintro ⟨hhu, hld⟩
    have gap : max h l < min u d := lt_min (max_lt hhu hlu) (max_lt hhd hld)
    obtain ⟨y, hy, hy'⟩ := exists_between gap
    exact ⟨y, lt_of_le_of_lt (le_max_left _ _) hy,
      le_trans hy'.le (min_le_right _ _), lt_of_le_of_lt (le_max_right _ _) hy,
      lt_of_lt_of_le hy' (min_le_left _ _)⟩

/-- A common positive threshold for two different quadratic cells forces
the homogeneous relation between their coefficient forms. -/
theorem simultaneous_cells {a b c d w s : ℝ}
    (h₁ : a*w^2+b = s*w) (h₂ : c*w^2+d = s*w) :
    s^2*(a-c)*(d-b) = (a*d-c*b)^2 := by
  have gap : (a-c)*w^2 = d-b := by linear_combination h₁ - h₂
  have det : a*d-c*b = s*w*(a-c) := by linear_combination a*h₂ - c*h₁
  rw [← gap, det]
  ring

end AlgebraicHeights

#print axioms AlgebraicHeights.protected_log
#print axioms AlgebraicHeights.affine_zeros
#print axioms AlgebraicHeights.square_refinement
#print axioms AlgebraicHeights.threshold_cell
#print axioms AlgebraicHeights.eligible
#print axioms AlgebraicHeights.simultaneous_cells
