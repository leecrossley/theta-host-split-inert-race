import Mathlib

namespace DyadicPowers

/-- Telescope the selected profile before bounding any of its zero coefficients. -/
theorem telescope (V : ℕ → ℝ) (n : ℕ) :
    ∑ k ∈ Finset.range n, (V (k+1)-V k) = V n-V 0 := by
  induction n with
  | zero => simp
  | succ n ih => rw [Finset.sum_range_succ, ih]; ring

/-- Certified real weights leave room for every omitted term and the prime-only remainder. -/
theorem small_move {selected tail err change : ℝ}
    (hs : selected ≤ (44/7)*175/32768) (ht : tail ≤ 2*(649/1000))
    (he : err ≤ 1/10^19) (h : |change| ≤ selected+tail+err) :
    |change| < 4/3 := by linarith

/-- This covers both the binary arc and the independently used decimal arc. -/
theorem phase_margin {x e : ℝ} (hx : 17/100 ≤ x) (he : |e| ≤ 1/1000) :
    1/6 < x+e := by linarith [(abs_le.mp he).1]

/-- A single true ordinate can certify a two-sided oscillation via its Fourier coefficient. -/
theorem coefficient_margin {g s : ℝ} (hg : 7/5 < g) (hs : 1/2 ≤ s) :
    7/5 < 2*g*s := by
  have h : 0 ≤ 2*g := by linarith
  have h' := mul_le_mul_of_nonneg_left hs h
  linarith

/-- Two distinct conjugate pairs in Parseval's sum force the required supremum. -/
theorem two_mode_margin {p q L : ℝ} (hp : 7/5 < p) (hq : 7/5 < q)
    (hL : 0 ≤ L) (h : 2*(p^2+q^2) ≤ L^2) : 14/5 < L := by
  by_contra! hle
  nlinarith [sq_nonneg (p-7/5), sq_nonneg (q-7/5), mul_nonneg hL (sub_nonneg.mpr hle)]

/-- The uniform norm triangle inequality excludes a small gap between exceptional exponents. -/
theorem separated {first second difference threshold : ℝ} (h : difference ≤ first+second)
    (large : 14/5 < difference) (small : threshold < 7/5) :
    ¬ (first ≤ threshold ∧ second ≤ threshold) := by
  rintro ⟨h₁,h₂⟩
  linarith

/-- A repeated box among the first M+1 phase vectors gives a positive index at most M. -/
theorem box_gap {u v M : ℕ} (hu : u < v) (hv : v ≤ M) :
    0 < v-u ∧ v-u ≤ M := by omega

end DyadicPowers

#print axioms DyadicPowers.telescope
#print axioms DyadicPowers.small_move
#print axioms DyadicPowers.phase_margin
#print axioms DyadicPowers.coefficient_margin
#print axioms DyadicPowers.two_mode_margin
#print axioms DyadicPowers.separated
#print axioms DyadicPowers.box_gap
