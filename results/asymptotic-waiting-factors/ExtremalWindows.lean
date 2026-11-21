import Mathlib

namespace ExtremalWindows

/-- One tilted two-point measure annihilates an exponent at a prescribed ordinate. -/
theorem atom_zero (b u g : ℝ) (h : b*g = Real.pi/2) :
    Complex.cosh ((b:ℂ)*(((u:ℂ)+(g:ℂ)*Complex.I)-u)) = 0 := by
  have hh : (b:ℂ)*(g:ℂ) = ((Real.pi/2 : ℝ):ℂ) := by exact_mod_cast h
  have hz : (b:ℂ)*(((u:ℂ)+(g:ℂ)*Complex.I)-u) =
      ((Real.pi/2 : ℝ):ℂ)*Complex.I := by rw [← hh]; ring
  rw [hz, Complex.cosh_mul_I, ← Complex.ofReal_cos, Real.cos_pi_div_two]
  simp

/-- The local signed probe leaves both probability densities nonnegative. -/
theorem positive_tests {k p : ℝ} (h : |p| ≤ k/2) :
    k/2 ≤ k+p ∧ k/2 ≤ k-p ∧ k+p ≤ 3*k/2 ∧ k-p ≤ 3*k/2 := by
  rcases abs_le.mp h with ⟨hl, hr⟩
  constructor
  · linarith
  constructor
  · linarith
  constructor <;> linarith

/-- A compact translated support fits strictly inside the available log window. -/
theorem support_in_window {r e W h t : ℝ} (hw : 2*r+2*e < W)
    (ht : |t| ≤ r+e) :
    h < h+W/2+t ∧ h+W/2+t < h+W := by
  rcases abs_le.mp ht with ⟨hl, hr⟩
  constructor <;> linarith

/-- A finite set of ordinates above the verified height supplies the stated count. -/
theorem count_from_reciprocals {s : Finset ℝ} (hne : s.Nonempty)
    (height : ∀ g ∈ s, (8590 : ℝ) < g) (mass : 14487 < 8590 * ∑ g ∈ s, 1/g) :
    14488 ≤ s.card := by
  have hlt : ∑ g ∈ s, 8590 * (1/g) < ∑ _g ∈ s, (1:ℝ) := by
    apply Finset.sum_lt_sum_of_nonempty hne
    intro g hg
    have hp := height g hg
    have hgpos : 0 < g := by linarith
    rw [mul_one_div]
    exact (div_lt_one hgpos).2 hp
  rw [← Finset.mul_sum] at hlt
  simp only [Finset.sum_const, nsmul_eq_mul, mul_one] at hlt
  have : (14487 : ℝ) < s.card := lt_trans mass hlt
  exact_mod_cast this

/-- The upper half-wave factor decreases with the ordinate. -/
theorem half_wave_factor {g : ℝ} (hg : 8590 < g) :
    2*Real.exp (Real.pi/g) < 2*Real.exp (Real.pi/8590) := by
  have h := div_lt_div_of_pos_left Real.pi_pos (by norm_num : (0:ℝ)<8590) hg
  exact mul_lt_mul_of_pos_left (Real.exp_lt_exp.mpr h) (by norm_num)

/-- A two-sign amplitude is at most half the range between its witnesses. -/
theorem balanced_pair (a b t : ℝ) (ha : t ≤ a) (hb : t ≤ -b) :
    2*t ≤ a-b := by linarith

/-- A decreasing covariance bounds the Cauchy--Schwarz product. -/
theorem decreasing_covariance {v₀ v₁ v₂ a m : ℝ}
    (horder : v₁ ≤ v₂) (hgap : v₁ ≤ v₀) (hscale : 0 < a)
    (h : (a*m)^2 ≤ (v₀-v₁)*(v₀-v₂)) :
    m ≤ (v₀-v₁)/a := by
  apply (le_div_iff₀ hscale).mpr
  have hb := mul_le_mul_of_nonneg_left (show v₀-v₂ ≤ v₀-v₁ by linarith) (sub_nonneg.mpr hgap)
  by_contra! hbad
  have hsq := mul_self_lt_mul_self (sub_nonneg.mpr hgap) hbad
  nlinarith only [h,hb,hsq]

end ExtremalWindows

#print axioms ExtremalWindows.atom_zero
#print axioms ExtremalWindows.positive_tests
#print axioms ExtremalWindows.support_in_window
#print axioms ExtremalWindows.count_from_reciprocals
#print axioms ExtremalWindows.half_wave_factor
#print axioms ExtremalWindows.balanced_pair
#print axioms ExtremalWindows.decreasing_covariance
