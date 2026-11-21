import Mathlib

namespace IsolatedPeaks

/-- Triangle inequalities for square-root deficits imply a signed height drop. -/
theorem chord_bound {a b c : ℝ} (hc : 0 ≤ c) (h : a ≤ b+c)
    (ha : b ≤ a) : (a-b)^2 ≤ c^2 := by
  have hh := mul_nonneg (show 0 ≤ c-(a-b) by linarith) (show 0 ≤ c+(a-b) by linarith)
  nlinarith

/-- Two high scores cannot realize a displacement with deficit at least 56. -/
theorem two_high_scores {a b c : ℝ}
    (h : a ≤ b+c) (ha : 0 ≤ a) (h56 : 56 ≤ a^2)
    (hb : b^2 ≤ 13713/1000) (hc : c^2 ≤ 13713/1000) : False := by
  have h₁ := mul_nonneg (show 0 ≤ b+c-a by linarith) (show 0 ≤ b+c+a by linarith)
  nlinarith [sq_nonneg (b-c)]

/-- A displacement whose square-root deficit exceeds sqrt(56) forces the forty-unit drop. -/
theorem forty_below {a b : ℝ} (ha : 0 ≤ a) (hb : 0 ≤ b)
    (hs : 56 ≤ a^2) (hk : b^2 ≤ 713/1000) : 40+2/10^12 < a^2-2*a*b := by
  have hb' : b < 85/100 := by nlinarith
  have hlimit : 748/100 ≤ a := by nlinarith
  have hprod := mul_nonneg (show 0 ≤ a-748/100 by linarith)
                            (show 0 ≤ a+748/100-2*b by linarith)
  nlinarith

theorem twenty_below {a b : ℝ} (ha : 0 ≤ a) (hb : 0 ≤ b)
    (hs : 56 ≤ a^2) (hk : b^2 ≤ 5713/1000) : 20+2/10^12 < a^2-2*a*b := by
  have hlimit : 7483/1000 ≤ a := by nlinarith
  have hb' : b < 2391/1000 := by nlinarith
  have hprod := mul_nonneg (show 0 ≤ a-7483/1000 by linarith)
                            (show 0 ≤ a+7483/1000-2*b by linarith)
  nlinarith

/-- In a window of the stated width the forbidden annulus gives a cluster of small diameter. -/
theorem cluster_span {l x y : ℝ} (hx : l ≤ x) (hy : y ≤ l+10^9)
    (forbid : 1/2 ≤ y-x → y-x ≤ 10^9 → False) : y-x < 1/2 := by
  by_contra! h
  exact forbid h (by linarith)

end IsolatedPeaks

#print axioms IsolatedPeaks.chord_bound
#print axioms IsolatedPeaks.two_high_scores
#print axioms IsolatedPeaks.forty_below
#print axioms IsolatedPeaks.twenty_below
#print axioms IsolatedPeaks.cluster_span
