import Mathlib

namespace ProfileInterval

/-- A recurring closed superlevel interval with the certified margin gives arbitrarily late intervals. -/
theorem recurring_no_loss {F U : ℝ → ℝ}
    (profile : ∀ t ∈ Set.Icc (1008796429610/10^12) (3221800332553/10^12),
        -(19735/10000:ℝ)+1/10^11 < U t)
    (late : ∀ B : ℝ, ∃ T > B, ∀ t ∈ Set.Icc (1008796429610/10^12) (3221800332553/10^12),
        |F (T+t)-U t| ≤ 1/10^11) :
    ∀ B : ℝ, ∃ s > B, ∀ x ∈ Set.Icc s (s+2213003902943/10^12), -(19735/10000:ℝ) < F x := by
  intro B
  obtain ⟨T,hT,close⟩ := late (B-1008796429610/10^12)
  refine ⟨T+1008796429610/10^12,by linarith,?_⟩
  intro x hx
  rcases hx with ⟨hl,hr⟩
  have ht : x-T ∈ Set.Icc (1008796429610/10^12) (3221800332553/10^12) := by
    constructor <;> linarith
  have h := (abs_le.mp (close (x-T) ht)).1
  rw [show T+(x-T)=x by ring] at h
  linarith [profile (x-T) ht]

/-- A doubling ending by 2 exp(L) H has its log centre in an interval of length L. -/
theorem move_fits {H y L : ℝ} (hH : 0 < H) (hy : H < y)
    (hend : 2*y ≤ (2*Real.exp L)*H) :
    Real.log y + (3/2)*Real.log 2 ∈
      Set.Icc (Real.log H+(3/2)*Real.log 2) (Real.log H+(3/2)*Real.log 2+L) := by
  have hp : 0 < y := lt_trans hH hy
  have hlog := Real.strictMonoOn_log (Set.mem_Ioi.mpr hH) (Set.mem_Ioi.mpr hp) hy
  have hupper := Real.log_le_log hp (show y ≤ Real.exp L*H by linarith)
  rw [Real.log_mul (Real.exp_ne_zero _) hH.ne',Real.log_exp] at hupper
  constructor <;> linarith

/-- Interior and exterior certificates bound the length of the component containing the witnesses. -/
theorem component_length {U : ℝ → ℝ} {a b l r d q : ℝ}
    (hd : 0 < d) (hl : l ∈ Set.Ioo a b) (hr : r ∈ Set.Ioo a b)
    (closed : ∀ x ∈ Set.Icc a b, q ≤ U x)
    (outside : U (l-d) < q ∧ U (r+d) < q) :
    r-l < b-a ∧ b-a < r-l+2*d := by
  rcases hl with ⟨hal,hlb⟩
  rcases hr with ⟨har,hrb⟩
  have left : l-d < a := by
    by_contra! h
    have := closed (l-d) ⟨h,by linarith⟩
    linarith [outside.1]
  have right : b < r+d := by
    by_contra! h
    have := closed (r+d) ⟨by linarith,h⟩
    linarith [outside.2]
  constructor <;> linarith

end ProfileInterval

#print axioms ProfileInterval.recurring_no_loss
#print axioms ProfileInterval.move_fits
#print axioms ProfileInterval.component_length
