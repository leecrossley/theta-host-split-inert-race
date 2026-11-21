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


theorem sphere_support (h : ℝ) :
    let L := (h+162000)/16999
    let X := h+81000+L
    let R := 80000+L
    X-R = h+1000 ∧ X+R < (17001/16999)*h+161020 := by
  dsimp
  constructor
  · ring
  · linarith

theorem stop_band (X L : ℝ) (size : 17000*L = X+81000) (hL : 0 ≤ L) :
    (X+80000)/2-8544*L ≤ -500 := by linarith

theorem overlap_witness {left right centre radius width offset : ℝ}
    (early : left < centre-radius) (late : centre+radius-width < left)
    (in_window : centre-radius ≤ right ∧ right ≤ centre+radius) :
    left-offset < right-offset ∧ right-offset < left+width-offset := by
  constructor <;> linarith [in_window.1,in_window.2]

theorem square_move {first last cap : ℝ} (hfirst : 0 ≤ first ∧ first ≤ cap)
    (hlast : 0 ≤ last ∧ last ≤ cap) : |last-first| ≤ cap := by
  exact abs_le.mpr ⟨by linarith [hfirst.1,hfirst.2,hlast.1,hlast.2],
    by linarith [hfirst.1,hfirst.2,hlast.1,hlast.2]⟩

theorem carrier_support (h : ℝ) :
    let L := (h+300)/16999
    let X := h+110+L
    let R := 108+L
    X-R = h+2 ∧ X+R < (17001/16999)*h+219 := by
  dsimp
  constructor
  · ring
  · linarith

end WaitingBound

#print axioms WaitingBound.box_support
#print axioms WaitingBound.power_cutoff
#print axioms WaitingBound.gain_of_test
#print axioms WaitingBound.loss_of_test
#print axioms WaitingBound.moment_margins
#print axioms WaitingBound.sphere_support
#print axioms WaitingBound.stop_band
#print axioms WaitingBound.overlap_witness
#print axioms WaitingBound.square_move
#print axioms WaitingBound.carrier_support
