import Mathlib

namespace FactorObstruction

/-- Uniformly recurring copies of a profile with a strict margin give whole intervals without a loss. -/
theorem recurring_no_loss {F U : ℝ → ℝ}
    (profile : ∀ t ∈ Set.Icc (1014/1000) (3214/1000), -(9/5:ℝ) < U t)
    (late : ∀ B : ℝ, ∃ T > B, ∀ t ∈ Set.Icc (1014/1000) (3214/1000), |F (T+t)-U t| ≤ 1/10) :
    ∀ B : ℝ, ∃ s > B, ∀ x ∈ Set.Icc s (s+11/5), -(19735/10000:ℝ) < F x := by
  intro B
  obtain ⟨T,hT,h⟩ := late (B-1014/1000)
  refine ⟨T+1014/1000,by linarith,?_⟩
  intro x hx
  rcases hx with ⟨hl,hr⟩
  have ht : x-T ∈ Set.Icc (1014/1000) (3214/1000) := by constructor <;> linarith
  have a := (abs_le.mp (h (x-T) ht)).1
  rw [show T+(x-T)=x by ring] at a
  linarith [profile (x-T) ht]

/-- The logarithmic width available to a doubling is log(C/2). -/
theorem move_fits {H y s : ℝ} (hH : 0 < H) (hy : H < y)
    (hstart : s = Real.log H + (3/2)*Real.log 2)
    (hend : 2*y ≤ (2*Real.exp (11/5))*H) :
    Real.log y + (3/2)*Real.log 2 ∈ Set.Icc s (s+11/5) := by
  have hp : 0 < y := lt_trans hH hy
  have hlog := Real.strictMonoOn_log (Set.mem_Ioi.mpr hH) (Set.mem_Ioi.mpr hp) hy
  have hz : y ≤ Real.exp (11/5)*H := by linarith
  have hupper := Real.log_le_log hp hz
  rw [Real.log_mul (Real.exp_ne_zero _) hH.ne',Real.log_exp] at hupper
  constructor <;> linarith

end FactorObstruction

#print axioms FactorObstruction.recurring_no_loss
#print axioms FactorObstruction.move_fits
