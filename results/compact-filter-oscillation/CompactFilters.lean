import Mathlib

namespace CompactFilters

theorem weighted_bound {I : Type*} [Fintype I] (w x h : I → ℝ) (b : ℝ)
    (hw : ∀ i, 0 ≤ w i) (hx : ∀ i, x i ^ 2 ≤ b^2 * h i) :
    (∑ i, w i * x i ^ 2) ≤ b^2 * ∑ i, w i * h i := by
  rw [Finset.mul_sum]
  exact Finset.sum_le_sum fun i _ => calc
    w i * x i ^ 2 ≤ w i * (b^2 * h i) := mul_le_mul_of_nonneg_left (hx i) (hw i)
    _ = b^2 * (w i * h i) := by ring

theorem quadrature_radius (floor cost radius norm observed : ℝ)
    (hc : 0 < cost) (hn : 0 < norm) (ho : 0 ≤ observed)
    (margin : radius^2 * cost < floor)
    (control : floor * norm^2 ≤ observed^2 * cost) : radius*norm < observed := by
  by_contra! h
  have square := mul_self_le_mul_self ho h
  have hs : observed^2 ≤ (radius*norm)^2 := by simpa [pow_two] using square
  have hc' := mul_le_mul_of_nonneg_right hs hc.le
  have hm := mul_lt_mul_of_pos_right margin (sq_pos_of_pos hn)
  nlinarith only [control,hc',hm]

theorem affine_error (lo hi b e : ℝ) (sep : 2*e < hi-lo) :
    e < |hi-b| ∨ e < |lo-b| := by
  by_contra! h
  obtain ⟨_,hh⟩ := abs_le.mp h.1
  obtain ⟨hl,_⟩ := abs_le.mp h.2
  linarith

theorem two_sided (a b c : ℝ) (h : 0 < c) (ha : c ≤ a) (hb : b ≤ -c) :
    b < 0 ∧ 0 < a := by constructor <;> linarith

theorem unrepresented (zeros killed : ℕ) (h : killed < zeros) :
    ¬ zeros ≤ killed := Nat.not_le_of_gt h

theorem frame_perturbation (mass error cost radius : ℝ)
    (h : radius^2 * cost < mass-2*error) :
    radius^2 * cost + 2*error < mass := by linarith

theorem positive_probe (signal leak mass radius : ℝ)
    (hm : 0 < mass) (h : radius * mass < signal - leak) :
    radius < (signal - leak) / mass := (lt_div_iff₀ hm).mpr h

theorem grid_rounding (r δ b ε M K s : ℝ)
    (hM : 0 ≤ M) (hK : 0 ≤ K) (he : 0 ≤ ε)
    (hδ : ε * (1 + K) < δ) (heb : ε ≤ b)
    (hsmall : M ≤ K → r + δ < s)
    (hlarge : K < M → r + δ + b * M < s) :
    r < s - ε * (1 + M) := by
  by_cases h : M ≤ K
  · have horder : 1 + M ≤ 1 + K := by linarith
    linarith [mul_le_mul_of_nonneg_left horder he, hsmall h]
  · have hprod := mul_le_mul_of_nonneg_right heb hM
    have hcut := mul_nonneg he hK
    nlinarith only [hlarge (lt_of_not_ge h), hprod, hcut, hδ]

theorem rounding_exponent (d k : ℝ) (hd : 0 < d) (hk : 0 < k) :
    ∃ e : ℝ, 0 < e ∧ d - k < e ∧ e < d := by
  refine ⟨d - min d k / 2, ?_, ?_, ?_⟩
  · linarith [min_le_left d k]
  · linarith [min_le_right d k]
  · linarith [lt_min hd hk]

end CompactFilters

#print axioms CompactFilters.weighted_bound
#print axioms CompactFilters.quadrature_radius
#print axioms CompactFilters.affine_error
#print axioms CompactFilters.two_sided
#print axioms CompactFilters.unrepresented
#print axioms CompactFilters.frame_perturbation
#print axioms CompactFilters.positive_probe
#print axioms CompactFilters.grid_rounding
#print axioms CompactFilters.rounding_exponent
