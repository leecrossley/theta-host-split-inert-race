import Mathlib

namespace PairedWindows

/-- The cubic majorizes the absolute half-angle sine when evaluated at its square. -/
theorem cubic_majorant (z : ℝ) :
    z ≤ (16*(z^2)^3 - 48*(z^2)^2 + 75*z^2 + 11)/54 := by
  have ident : (16*(z^2)^3 - 48*(z^2)^2 + 75*z^2 + 11)/54 - z =
      (8/27:ℝ)*(z-1/2)^2*(z-1)^2*((z+3/2)^2+1/2) := by ring
  exact sub_nonneg.mp (ident ▸ (by positivity : 0 ≤
    (8/27:ℝ)*(z-1/2)^2*(z-1)^2*((z+3/2)^2+1/2)))

/-- Curvature, phase and all-ordinate errors preserve the signed profile margin. -/
theorem spectral_budget {u p shift interp rem : ℝ}
    (h : p-shift-interp-rem ≤ u) (hp : 164/100 < p) (hs : shift < 43/1000)
    (hi : interp < 24/1000) (hr : rem < 144/100) : 13/100 < u := by linarith

/-- The explicit windows give both requested deadline comparisons. -/
theorem deadlines :
    Real.log 200 < (45/8:ℝ) ∧ 554 < 2*Real.exp (45/8) ∧ 806 < 2*Real.exp 6 := by
  have h₁ := Real.sum_le_exp_of_nonneg (x := (45/8:ℝ)) (by norm_num) 22
  have h₂ := Real.sum_le_exp_of_nonneg (x := (6:ℝ)) (by norm_num) 22
  norm_num [Finset.sum_range_succ, Nat.factorial_succ] at h₁ h₂
  constructor
  · apply (Real.log_lt_iff_lt_exp (by norm_num : (0:ℝ) < 200)).2
    linarith
  · constructor <;> linarith

/-- If a signed profile and its approximation have these margins then all moves have the same sign. -/
theorem strict_move {f u : ℝ} (close : |f-u| ≤ 1/100)
    (profile : 13/100 < u) : 1/10 < f := by
  linarith [(abs_le.mp close).1]

/-- Oddness transfers a positive interval to a negative interval of the same length. -/
theorem reflection {U : ℝ → ℝ} {l r b : ℝ} (odd : ∀ t, U (-t) = -U t)
    (positive : ∀ t ∈ Set.Icc l r, b < U t) :
    ∀ t ∈ Set.Icc (-r) (-l), b < -U t := by
  intro t ht
  obtain ⟨hl,hr⟩ := ht
  have hm : -t ∈ Set.Icc l r := ⟨by linarith,by linarith⟩
  simpa only [odd t] using positive (-t) hm

/-- Every eligible doubling lies in the certified closed interval of log centres. -/
theorem move_fits {H y L : ℝ} (hH : 0 < H) (hy : H < y)
    (hend : 2*y ≤ (2*Real.exp L)*H) :
    Real.log y + (3/2)*Real.log 2 ∈
      Set.Icc (Real.log H+(3/2)*Real.log 2) (Real.log H+(3/2)*Real.log 2+L) := by
  have hp : 0 < y := lt_trans hH hy
  have hlog := Real.strictMonoOn_log (Set.mem_Ioi.mpr hH) (Set.mem_Ioi.mpr hp) hy
  have hupper := Real.log_le_log hp (show y ≤ Real.exp L*H by linarith)
  rw [Real.log_mul (Real.exp_ne_zero _) hH.ne',Real.log_exp] at hupper
  constructor <;> linarith

/-- The 6.235-unit window excludes the four-digit deadline. -/
theorem thousand_deadline : 1020 < 2 * Real.exp (1247/200 : ℝ) := by
  have h := Real.sum_le_exp_of_nonneg (x := (1247/200:ℝ)) (by norm_num) 24
  norm_num [Finset.sum_range_succ, Nat.factorial_succ] at h
  linarith

/-- Uniform spectral and prime-only errors preserve the loss threshold. -/
theorem shallow_loss {f u : ℝ} (floor : -(197:ℝ)/100 < u)
    (error : |f-u| ≤ 1/1000) : -(19735:ℝ)/10000 < f := by
  linarith [(abs_le.mp error).1]

/-- An interior small loss and a large gain coexist with a threshold-free interval. -/
theorem both_move_signs {lo hi : ℝ} (fall : lo < -(1:ℝ)/25)
    (rise : 12 < hi) : lo < 0 ∧ 19735/10000 < hi := by
  constructor <;> linarith

end PairedWindows

#print axioms PairedWindows.cubic_majorant
#print axioms PairedWindows.spectral_budget
#print axioms PairedWindows.deadlines
#print axioms PairedWindows.strict_move
#print axioms PairedWindows.reflection
#print axioms PairedWindows.move_fits

#print axioms PairedWindows.thousand_deadline
#print axioms PairedWindows.shallow_loss
#print axioms PairedWindows.both_move_signs
