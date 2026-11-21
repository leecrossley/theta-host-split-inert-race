import Mathlib

namespace DistantInterval

/-- The four certified errors and ordinate-grid lower bound imply the profile bound. -/
theorem spectral_budget {u p move interp rem : ℝ}
    (hu : u ≤ -p + move + interp + rem) (hp : 2/100 < p)
    (hm : move < 3/100) (hi : interp < 67/1000) (hr : rem < 172/100) :
    u < 18/10 := by linarith

/-- At this explicit interval the prime-only error and the profile bound exclude a gain. -/
theorem explicit_no_gain {f u : ℝ} (h : |f-u| ≤ 1/(10:ℝ)^20) (hu : u < 18/10) :
    f < 19735/10000 := by
  have h' := (abs_le.mp h).2
  linarith

/-- An upper bound on the profile gives the upper bound for its negative on the reflected interval. -/
theorem reflected_interval {U : ℝ → ℝ} {l r b : ℝ} (odd : ∀ x, U (-x) = -U x)
    (upper : ∀ t ∈ Set.Icc l r, U t < b) :
    ∀ x ∈ Set.Icc (-r) (-l), -U x < b := by
  intro x hx
  obtain ⟨hx₁,hx₂⟩ := hx
  have hm : -x ∈ Set.Icc l r := ⟨by linarith,by linarith⟩
  simpa only [odd x] using upper (-x) hm

/-- A profile below 1.9 and recurrent approximations of error at most 0.05 give late intervals. -/
theorem late_upper {F U : ℝ → ℝ} {l r : ℝ}
    (profile : ∀ t ∈ Set.Icc l r, U t < 19/10)
    (late : ∀ B : ℝ, ∃ T > B, ∀ t ∈ Set.Icc l r, |F (T+t)-U t| ≤ 1/20) :
    ∀ B : ℝ, ∃ s > B, ∀ x ∈ Set.Icc s (s+(r-l)), F x < 19735/10000 := by
  intro B
  obtain ⟨T,hT,close⟩ := late (B-l)
  refine ⟨T+l,by linarith,?_⟩
  intro x hx
  obtain ⟨hl,hr⟩ := hx
  have ht : x-T ∈ Set.Icc l r := ⟨by linarith,by linarith⟩
  have h := (abs_le.mp (close (x-T) ht)).2
  rw [show T+(x-T)=x by ring] at h
  linarith [profile (x-T) ht]

/-- Every permitted dyadic move has its log centre inside the closed certified window. -/
theorem move_fits {H y L : ℝ} (hH : 0 < H) (hy : H < y)
    (hend : 2*y ≤ (2*Real.exp L)*H) :
    Real.log y + (3/2)*Real.log 2 ∈
      Set.Icc (Real.log H+(3/2)*Real.log 2) (Real.log H+(3/2)*Real.log 2+L) := by
  have hp : 0 < y := lt_trans hH hy
  have hlog := Real.strictMonoOn_log (Set.mem_Ioi.mpr hH) (Set.mem_Ioi.mpr hp) hy
  have hupper := Real.log_le_log hp (show y ≤ Real.exp L*H by linarith)
  rw [Real.log_mul (Real.exp_ne_zero _) hH.ne',Real.log_exp] at hupper
  constructor <;> linarith

end DistantInterval

#print axioms DistantInterval.spectral_budget
#print axioms DistantInterval.explicit_no_gain
#print axioms DistantInterval.reflected_interval
#print axioms DistantInterval.late_upper
#print axioms DistantInterval.move_fits
