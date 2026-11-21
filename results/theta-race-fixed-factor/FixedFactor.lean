import Mathlib

namespace FixedFactor

/-- The two nonnegative tests exceed the required signed thresholds despite all error bounds. -/
theorem test_margins {m v pg pl : ℝ}
    (hm : |m| ≤ 225/900 + 1/10)
    (hv : 33428/10000 - (7/2)*217/900 - 1/10 ≤ v)
    (hg : pg ≤ 1 + 1/90^2) (hl : pl ≤ 1 + 1/90^2) :
    (19735/10000)*pg < m+v ∧ m-v < -(19735/10000)*pl := by
  rcases abs_le.mp hm with ⟨h₁,h₂⟩
  constructor <;> linarith

/-- Convert the fixed log window to a deadline for both ends of a doubling. -/
theorem endpoint {H x : ℝ} (hH : 0 < H)
    (hx : x ∈ Set.Icc (Real.log H + 55) (Real.log H + 175)) :
    H < Real.exp (x-3*Real.log 2/2) ∧
    2*Real.exp (x-3*Real.log 2/2) < Real.exp 175 * H := by
  rcases hx with ⟨hlo,hhi⟩
  have ha : 0 < Real.log 2 := Real.log_pos (by norm_num)
  have hb := Real.log_two_lt_d9
  constructor
  · calc H = Real.exp (Real.log H) := (Real.exp_log hH).symm
         _ < _ := Real.exp_lt_exp.mpr (by linarith)
  · calc 2*Real.exp (x-3*Real.log 2/2)
           = Real.exp (Real.log 2+(x-3*Real.log 2/2)) := by
               rw [Real.exp_add, Real.exp_log (by norm_num : (0:ℝ)<2)]
         _ < Real.exp (175+Real.log H) := Real.exp_lt_exp.mpr (by linarith)
         _ = _ := by rw [Real.exp_add,Real.exp_log hH]

end FixedFactor

#print axioms FixedFactor.test_margins
#print axioms FixedFactor.endpoint
