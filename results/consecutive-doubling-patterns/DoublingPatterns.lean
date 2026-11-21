import Mathlib

namespace DoublingPatterns

open scoped BigOperators

theorem sign_preserved {s x e t m : ℝ} (hs : s=1 ∨ s= -1)
    (hx : m < s*x) (he : |e| < t) : m-t < s*(x+e) := by
  rw [abs_lt] at he
  rcases hs with h | h <;> subst s <;> linarith

theorem separated_signs {x y e : ℝ} (he : 0 ≤ e) (hx : e < |x|)
    (hy : e < |y|) (h : |x-y| < 2*e) : (0 < x ↔ 0 < y) := by
  rw [abs_lt] at h
  rw [lt_abs] at hx hy
  constructor <;> intro hp <;> rcases hx with hx | hx <;> rcases hy with hy | hy <;> linarith

theorem move_error {x y x' y' e a : ℝ} (ha : 0 < a)
    (hx : |x-x'| ≤ e) (hy : |y-y'| ≤ e) :
    |(y-x)/a-(y'-x')/a| ≤ 2*e/a := by
  rw [← sub_div, abs_div, abs_of_pos ha]
  apply div_le_div_of_nonneg_right _ ha.le
  rw [abs_le] at hx hy ⊢
  constructor <;> linarith

theorem telescoping (S : ℕ → ℝ) (n : ℕ) :
    ∑ j ∈ Finset.range n, (S (j+1)-S j) = S n-S 0 := by
  induction n with
  | zero => simp
  | succ n ih => rw [Finset.sum_range_succ, ih]; ring

theorem excluded_run {m e M : ℝ} (he : 2*e < m-2*M)
    {s t : ℝ} (hs : 1-M-e ≤ s) (ht : t ≤ 1+M+e) : t-s < m := by
  linarith

theorem word_count {W C : Type*} [Fintype W] [Fintype C]
    (cell : W → C) (h : Function.Injective cell) : Fintype.card W ≤ Fintype.card C :=
  Fintype.card_le_of_injective cell h

theorem reflect_difference (f : ℝ → ℝ) (he : Function.Even f) (x a : ℝ) :
    (f (-x-a+a)-f (-x-a))/a = -(f (x+a)-f x)/a := by
  rw [show -x-a+a = -x by ring, show -x-a = -(x+a) by ring, he x, he (x+a)]
  ring

theorem cubic_majorant (z : ℝ) :
    z ≤ (11 + z^2*(75+z^2*(16*z^2-48)))/54 := by
  have h : 0 ≤ (8/27:ℝ)*(z-1/2)^2*(z-1)^2*((z+3/2)^2+1/2) := by positivity
  nlinarith only [h]

end DoublingPatterns

#print axioms DoublingPatterns.sign_preserved
#print axioms DoublingPatterns.separated_signs
#print axioms DoublingPatterns.move_error
#print axioms DoublingPatterns.telescoping
#print axioms DoublingPatterns.excluded_run
#print axioms DoublingPatterns.word_count
#print axioms DoublingPatterns.reflect_difference
#print axioms DoublingPatterns.cubic_majorant
