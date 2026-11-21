import Mathlib

namespace ExcursionProportions

noncomputable def q (x : ℝ) : ℝ :=
  (16515360000 / 325057991399) * (x - 5/27) * (x + 8/5)^2 * (10838/3823 - x)

/-- The signed indicator minorant is nonpositive on the closed lower half-line. -/
theorem minorant_left {x : ℝ} (hx : x ≤ 5/27) : q x ≤ 0 := by
  have ha : 0 ≤ (16515360000 / 325057991399 : ℝ) := by norm_num
  have hb : 0 ≤ 10838/3823 - x := by linarith
  exact mul_nonpos_of_nonpos_of_nonneg
    (mul_nonpos_of_nonpos_of_nonneg (mul_nonpos_of_nonneg_of_nonpos ha (sub_nonpos.mpr hx))
      (sq_nonneg _)) hb

/-- The same quartic is bounded above by one on the whole real line. -/
theorem minorant_cap (x : ℝ) : q x ≤ 1 := by
  have pos : 0 ≤ 16515360000*x^2 + 67380496000*x + 91322607600 := by
    nlinarith [sq_nonneg (16515360000*x+33690248000)]
  have h := mul_nonneg (sq_nonneg (x-39/20)) pos
  unfold q
  nlinarith only [h]

theorem signed_transfer {p r e f : ℝ} (hp : 5 < p) (hr : |r| ≤ 149/100)
    (he : |e| ≤ 152/100+7/1000) (hf : f = p+r+e) : 19735/10000 < f := by
  have h₁ := (abs_le.mp hr).1
  have h₂ := (abs_le.mp he).1
  linarith

/-- The elementary moment inequality after applying Cauchy--Schwarz on a tail. -/
theorem fourth_moment_tail {m t s p k : ℝ} (hk : 0 < k) (ht : t ≤ m)
    (hm : m ≤ t+s) (hs : s^2 ≤ p*k) : (m-t)^2/k ≤ p := by
  apply (div_le_iff₀ hk).2
  have h := mul_nonneg (show 0 ≤ s-(m-t) by linarith) (show 0 ≤ s+(m-t) by linarith)
  nlinarith

theorem asymptotic_margin :
    ((10065/10000:ℝ)-25/729)^2/(2*(30034/10000))
      -(923/100000)/(149/100)^2 > 3/20 := by norm_num

end ExcursionProportions

#print axioms ExcursionProportions.minorant_left
#print axioms ExcursionProportions.minorant_cap
#print axioms ExcursionProportions.signed_transfer
#print axioms ExcursionProportions.fourth_moment_tail
#print axioms ExcursionProportions.asymptotic_margin
