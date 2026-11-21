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


/-- A cell interpolated with nonnegative weights of sum at most one inherits a negative lower bound. -/
theorem knot_lower {u v a b t : ℝ} (ht : 0 ≤ t) (hu : 0 ≤ u) (hv : 0 ≤ v)
    (hs : u + v ≤ 1) (ha : -t ≤ a) (hb : -t ≤ b) : -t ≤ u*a + v*b := by
  nlinarith [mul_nonneg hu (sub_nonneg.mpr ha), mul_nonneg hv (sub_nonneg.mpr hb),
    mul_nonneg ht (sub_nonneg.mpr hs)]

/-- The endpoint factor itself fails at the left endpoint of a closed interval without a strict loss. -/
theorem gap_endpoint {F : ℝ → ℝ} {a b t : ℝ} (ha : 0 < a)
    (hgap : ∀ y ∈ Set.Icc a b, -t ≤ F y) :
    ¬ ∃ y > a, 2*y ≤ (2*b/a)*a ∧ F y < -t := by
  rintro ⟨y, hy, he, hl⟩
  rw [div_mul_cancel₀ _ ha.ne'] at he
  exact (not_lt_of_ge (hgap y ⟨hy.le, by linarith⟩)) hl

/-- Any larger deadline reaches a strict loss just beyond the right endpoint of the gap. -/
theorem beyond_gap {F : ℝ → ℝ} {a b H C t : ℝ} (hC : 0 ≤ C) (ha : a ≤ H)
    (hb : H ≤ b) (hroom : 2*b < C*a)
    (right : ∀ e > 0, ∃ y > b, y < b+e ∧ F y < -t) :
    ∃ y > H, 2*y ≤ C*H ∧ F y < -t := by
  have room : 0 < C*H/2-b := by nlinarith [mul_nonneg hC (sub_nonneg.mpr ha)]
  obtain ⟨y, hy, he, hl⟩ := right _ room
  exact ⟨y, lt_of_le_of_lt hb hy, by linarith, hl⟩

/-- A touching minimum is a strict loss precisely below its amplitude. -/
theorem strict_bridge {w t τ : ℝ} (hw : w = -τ) : w < -t ↔ t < τ := by
  rw [hw]
  constructor <;> intro h <;> linarith

end FactorObstruction

#print axioms FactorObstruction.recurring_no_loss
#print axioms FactorObstruction.move_fits

#print axioms FactorObstruction.knot_lower
#print axioms FactorObstruction.gap_endpoint
#print axioms FactorObstruction.beyond_gap
#print axioms FactorObstruction.strict_bridge
