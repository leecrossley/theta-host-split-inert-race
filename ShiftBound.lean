import Mathlib

open scoped BigOperators

noncomputable section

namespace ShiftBound

/-- The Fourier multiplier of a centred logarithmic difference, including its derivative. -/
noncomputable def multiplier (g t : ℝ) : ℝ :=
  if t = 0 then g else 2 * Real.sin (g*t/2) / t

private abbrev cap : ℝ := 69315 / 100000
def lower : Fin 5 → ℝ := ![48536362/10^9, 248572181/10^9, 434216664/10^9, 446011529/10^9, 579615600/10^9]
def upper : Fin 5 → ℝ := ![48536363/10^9, 248572182/10^9, 434216665/10^9, 446011530/10^9, 579615601/10^9]
private def floorTerm (l u : ℝ) : ℝ :=
  l*(1-(u*cap)^2/24) / (2*(1/4+u^2))

private theorem endpoint {l g u t : ℝ} (hl : 0 ≤ l) (hlg : l ≤ g) (hgu : g ≤ u)
    (hu : u ≤ 1) (ht : 0 ≤ t) (ht₂ : t ≤ Real.log 2) :
    floorTerm l u ≤ multiplier g t / (2*(1/4+g^2)) := by
  have hg := hl.trans hlg
  have hu₀ := hg.trans hgu
  have htc : t ≤ cap := ht₂.trans (Real.log_two_lt_d9.le.trans (by norm_num [cap]))
  have hc : (0:ℝ) ≤ cap := by norm_num
  have hprod : g*t ≤ u*cap := mul_le_mul hgu htc ht hu₀
  have hsq := pow_le_pow_left₀ (mul_nonneg hg ht) hprod 2
  have huc : u*cap ≤ 1 := le_trans (mul_le_mul_of_nonneg_right hu hc) (by norm_num)
  have huc₂ := pow_le_pow_left₀ (mul_nonneg hu₀ hc) huc 2
  have hpos : 0 ≤ 1-(u*cap)^2/24 := sub_nonneg.mpr
    ((div_le_iff₀ (by norm_num : (0:ℝ) < 24)).2 (huc₂.trans (by norm_num)))
  have hnum : l*(1-(u*cap)^2/24) ≤ g*(1-(g*t)^2/24) :=
    mul_le_mul hlg (sub_le_sub_left (div_le_div_of_nonneg_right hsq (by norm_num)) 1) hpos hg
  have hsin : g*(1-(g*t)^2/24) ≤ multiplier g t := by
    by_cases h : t = 0
    · subst t; simp [multiplier]
    · rw [multiplier, if_neg h, le_div_iff₀ (ht.lt_of_ne (Ne.symm h))]
      have hTaylor := Real.sin_ge_sub_cube (x := g*t/2) (by positivity)
      nlinarith only [hTaylor]
  exact div_le_div₀ ((mul_nonneg hl hpos).trans (hnum.trans hsin)) (hnum.trans hsin) (by positivity)
    (by nlinarith [pow_le_pow_left₀ hg hgu 2])

/-- Every width has this lower coefficient after centring and normalizing the window.
Here a is the square root of the width and c is the cosine of the zero phase. -/
theorem width_floor {a g c : ℝ} (ha : 1 < a) (hc : c ≤ 1) :
    1 / (2*(1/4+g^2)) ≤ (a^2+1-2*a*c) / (2*(a-1)^2*(1/4+g^2)) := by
  have ha₀ : 0 < a-1 := by linarith
  have hden : 0 < (1/4:ℝ)+g^2 := by positivity
  apply (div_le_div_iff₀ (by positivity) (by positivity)).2
  have h := mul_nonneg (mul_nonneg (by linarith : 0 ≤ 4*a) (sub_nonneg.mpr hc)) hden.le
  nlinarith only [h]

/-- Five ordinates in the certified brackets force a uniform bound through the derivative endpoint. -/
theorem uniform_short_shift (g : Fin 5 → ℝ) (hg : ∀ i, lower i ≤ g i ∧ g i ≤ upper i)
    {t : ℝ} (ht : 0 ≤ t) (ht₂ : t ≤ Real.log 2) :
    19735/10000 < ∑ i, multiplier (g i) t / (2*(1/4+(g i)^2)) := by
  apply lt_of_lt_of_le (b := ∑ i, floorTerm (lower i) (upper i))
  · norm_num [Fin.sum_univ_succ, floorTerm, lower, upper, cap]
  · apply Finset.sum_le_sum
    intro i _
    apply endpoint _ (hg i).1 (hg i).2 _ ht ht₂
    · fin_cases i <;> norm_num [lower]
    · fin_cases i <;> norm_num [upper]

end ShiftBound
end
