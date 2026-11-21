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


/-- The paired positive two-point law really kills a prescribed nonreal exponent. -/
theorem filter_zero (δ η d : ℝ) (h : d*η = Real.pi/2) :
    Complex.cosh ((d:ℂ)*(((δ:ℂ)+(η:ℂ)*Complex.I)-δ)) *
      Complex.cosh ((d:ℂ)*(((δ:ℂ)+(η:ℂ)*Complex.I)+δ)) = 0 := by
  have hh : (d:ℂ)*(η:ℂ) = ((Real.pi/2 : ℝ):ℂ) := by exact_mod_cast h
  have hz : (d:ℂ)*(((δ:ℂ)+(η:ℂ)*Complex.I)-δ) =
      ((Real.pi/2 : ℝ):ℂ)*Complex.I := by rw [← hh]; ring
  rw [hz, Complex.cosh_mul_I, ← Complex.ofReal_cos, Real.cos_pi_div_two]
  simp

/-- The retained part of an imaginary-frequency response is bounded by its phase variance. -/
theorem filter_retention (x v : ℝ) :
    1-x^2 ≤ 1-(Real.sin x)^2/(Real.cosh v)^2 := by
  have hc := Real.one_le_cosh v
  have hd : 1 ≤ (Real.cosh v)^2 := by nlinarith [sq_nonneg (Real.cosh v-1)]
  have hu : (Real.sin x)^2/(Real.cosh v)^2 ≤ (Real.sin x)^2 :=
    div_le_self (sq_nonneg _) hd
  linarith [show (Real.sin x)^2 ≤ x^2 from sq_le_sq.mpr Real.abs_sin_le_abs]

/-- Pointwise losses add as an upper bound on the deficit of a finite product. -/
theorem filter_product {ι : Type*} (s : Finset ι) (a : ι → ℝ)
    (ha : ∀ i, 0 ≤ a i) (hb : ∀ i, a i ≤ 1) :
    1-(∑ i ∈ s, a i) ≤ ∏ i ∈ s, (1-a i) := by
  classical
  induction s using Finset.induction_on with
  | empty => simp
  | @insert i s hi ih =>
    rw [Finset.sum_insert hi, Finset.prod_insert hi]
    have hn : 0 ≤ ∑ j ∈ s, a j := Finset.sum_nonneg fun j _ => ha j
    calc
      1-(a i+∑ j ∈ s, a j) ≤ (1-a i)*(1-∑ j ∈ s, a j) := by
        nlinarith [mul_nonneg (ha i) hn]
      _ ≤ _ := mul_le_mul_of_nonneg_left ih (sub_nonneg.mpr (hb i))

/-- Translating a compact filtered measure includes the full support cost on both sides. -/
theorem filter_support {h R r A : ℝ} (hR : R ≤ 108+1.000004+r)
    (hA : 221+2*r ≤ A) :
    let X := h+2+R
    X-R = h+2 ∧ X+R < h+A := by
  dsimp
  constructor
  · ring
  · linarith


/-- The sharp column majorant for the cosine density telescopes. -/
theorem cosine_column (n : ℕ) :
    (∑ j ∈ Finset.range n, (2 : ℝ) /
      (4 * ((j : ℝ) + 1)^2 - 1)) = 1 - 1 / (2*(n:ℝ)+1) := by
  induction n with
  | zero => norm_num
  | succ n ih =>
    rw [Finset.sum_range_succ, ih]
    push_cast
    have h1 : (0:ℝ) < 2 * (n:ℝ) + 1 := by positivity
    have h2 : (0:ℝ) < 2 * ((n:ℝ) + 1) + 1 := by positivity
    have h3 : (0:ℝ) < 4 * ((n:ℝ) + 1)^2 - 1 := by
      nlinarith [show (0:ℝ) ≤ (n:ℝ) by positivity, sq_nonneg (n:ℝ)]
    field_simp [ne_of_gt h1, ne_of_gt h2, ne_of_gt h3]
    ring

/-- The interpolation perturbation leaves the full support inside the stated endpoint. -/
theorem grouped_support (h L U R p q : ℝ)
    (hR : R ≤ 108 + 1.000004 + p + q)
    (hp : p ≤ 270*U) (hq : q ≤ 2*L+2) :
    h+2+2*R < h+225+540*U+4*L := by linarith

/-- Every eligible sphere centre leaves a strict margin at the high ordinates. -/
theorem finite_sphere_exponent {top x d L : ℝ} (hL : 2*8492*L = top+82)
    (hx : x < top+7) (hd : d ≤ 1/16) : x/2+d-8492*L < -37 := by
  linarith

/-- An unlisted real frequency in a covering cell overspends the Poisson budget. -/
theorem exterior_cell {g v b : ℝ} (h : |g-v| ≤ 3/2)
    (hb : b < 1/4) : b < 1/(1+(g-v)^2) := by
  have hsq : (g-v)^2 ≤ (3/2 : ℝ)^2 := sq_le_sq.mpr (by simpa only [abs_of_pos (by norm_num : (0:ℝ)<3/2)] using h)
  have hd : 0 < 1+(g-v)^2 := by positivity
  have hc : (4/13 : ℝ) ≤ 1/(1+(g-v)^2) := (le_div_iff₀ hd).2 (by nlinarith)
  linarith

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

#print axioms WaitingBound.filter_zero
#print axioms WaitingBound.filter_retention
#print axioms WaitingBound.filter_product
#print axioms WaitingBound.filter_support

#print axioms WaitingBound.cosine_column
#print axioms WaitingBound.grouped_support

#print axioms WaitingBound.finite_sphere_exponent
#print axioms WaitingBound.exterior_cell
