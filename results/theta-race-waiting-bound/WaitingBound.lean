import Mathlib

open MeasureTheory

namespace WaitingBound

/-- Logarithmic support of the compact probability measure, rounded outward. -/
theorem box_support {h X n R : ℝ}
    (hX : (1 - 1 / 3000) * X = h + 81000)
    (hn : n ≤ X + 1) (hR : R = 80000 + n / 3000) :
    h + 999 < X - R ∧ X + R < (3001 / 2999) * h + 162000 := by
  constructor <;> linarith

theorem power_cutoff {h : ℝ} (hh : 500000000 ≤ h) :
    (3001 / 2999) * h + 162000 < (1001 / 1000) * h := by
  linarith

/-- If a nonnegative weighted average exceeds a threshold, a point exceeds it. -/
theorem gain_of_test {α : Type*} [MeasurableSpace α] {μ : Measure α}
    {f p : α → ℝ} {t : ℝ} (hp : ∀ x, 0 ≤ p x)
    (hf : Integrable (fun x ↦ f x * p x) μ) (hw : Integrable p μ)
    (h : t * (∫ x, p x ∂μ) < ∫ x, f x * p x ∂μ) : ∃ x, t < f x := by
  by_contra! hle
  have hmono := integral_mono hf (hw.const_mul t) (fun x ↦
    mul_le_mul_of_nonneg_right (hle x) (hp x))
  rw [integral_const_mul] at hmono
  linarith

/-- Reversing the inequality gives a point with a sufficiently large loss. -/
theorem loss_of_test {α : Type*} [MeasurableSpace α] {μ : Measure α}
    {f p : α → ℝ} {t : ℝ} (hp : ∀ x, 0 ≤ p x)
    (hf : Integrable (fun x ↦ f x * p x) μ) (hw : Integrable p μ)
    (h : (∫ x, f x * p x ∂μ) < -t * (∫ x, p x ∂μ)) : ∃ x, f x < -t := by
  by_contra! hle
  have hmono := integral_mono (hw.const_mul (-t)) hf (fun x ↦
    mul_le_mul_of_nonneg_right (hle x) (hp x))
  rw [integral_const_mul] at hmono
  linarith

theorem moment_margins {u v massGain massLoss : ℝ}
    (hu : |u| ≤ 1/100) (hv : 333/100 ≤ v)
    (hg : massGain ≤ 1 + 1/2^40) (hl : massLoss ≤ 1 + 1/2^40) :
    (19735/10000) * massGain < u + v ∧ u - v < -(19735/10000) * massLoss := by
  rcases abs_le.mp hu with ⟨hu₁, hu₂⟩
  constructor <;> linarith

end WaitingBound

#print axioms WaitingBound.box_support
#print axioms WaitingBound.power_cutoff
#print axioms WaitingBound.gain_of_test
#print axioms WaitingBound.loss_of_test
#print axioms WaitingBound.moment_margins
