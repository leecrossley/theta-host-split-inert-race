import Mathlib

namespace GeometricGrid

theorem forced_point {I : Type*} [Fintype I] (w f : I → ℝ) (θ : ℝ)
    (hw : ∀ i, 0 ≤ w i) (hs : θ * ∑ i, w i < ∑ i, w i * f i) :
    ∃ i, θ < f i := by
  by_contra! hf
  have hsum : (∑ i, w i * f i) ≤ ∑ i, w i * θ :=
    Finset.sum_le_sum fun i _ => mul_le_mul_of_nonneg_left (hf i) (hw i)
  rw [← Finset.sum_mul] at hsum
  linarith [mul_comm θ (∑ i, w i)]

theorem error_transfer (p f : ℝ) (h : (1.985 : ℝ) < p)
    (he : |f-p| < (0.007 : ℝ)) : (1.9735 : ℝ) < f := by
  have he' := (abs_lt.mp he).1
  linarith

theorem stride {a h : ℝ} (hpos : 0 < h) (ha : 0 < a) (hh : h ≤ a) :
    ∃ k : ℕ, a ≤ k*h ∧ k*h ≤ 2*a := by
  refine ⟨⌈a/h⌉₊, ?_, ?_⟩
  · exact (div_le_iff₀ hpos).mp (Nat.le_ceil (a/h))
  · have hz : 0 ≤ a/h := (div_pos ha hpos).le
    have hc := Nat.ceil_lt_add_one hz
    have hx := mul_lt_mul_of_pos_right hc hpos
    rw [add_mul, div_mul_cancel₀ _ hpos.ne', one_mul] at hx
    linarith

theorem block_counts : 2*112-1=223 ∧ 2*160-1=319 ∧ 2*2048-1=4095 := by
  decide

theorem positive_product {u v : ℝ} (hu : 0 ≤ u) (hv : 0 ≤ v) : 0 ≤ u*v :=
  mul_nonneg hu hv

theorem alias_run (x y : ℝ) (hx : (14 : ℝ) < x) (he : |y-x| < 1) :
    (13 : ℝ) < y := by
  have h := (abs_lt.mp he).1
  linarith

end GeometricGrid

#print axioms GeometricGrid.forced_point
#print axioms GeometricGrid.error_transfer
#print axioms GeometricGrid.stride
#print axioms GeometricGrid.block_counts
#print axioms GeometricGrid.positive_product
#print axioms GeometricGrid.alias_run
