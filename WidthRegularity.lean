import Mathlib

namespace WidthRegularity

theorem atom_floor {a s : ℝ} (hs : 0 < s) : 0 < a^2+s/2 := by positivity

/-- The diagonal in a conjugate pair has half the original squared cosine coefficient. -/
theorem diagonal (a b c : ℝ) :
    a*a + (b/2)^2+(b/2)^2+(c/2)^2+(c/2)^2 = a^2+(b^2+c^2)/2 := by ring

/-- Aggregating integer multiplicities before Parseval can only increase a lower bound. -/
theorem square_weight (m : ℕ) (x : ℝ) : (m:ℝ)*x^2 ≤ (m:ℝ)^2*x^2 := by
  apply mul_le_mul_of_nonneg_right _ (sq_nonneg x)
  exact_mod_cast (show m ≤ m^2 by nlinarith [Nat.le_mul_self m])

theorem critical_shell (t : ℝ) (h : 0 < t) :
    t^3*(t/(t^2)^2) = 1 := by
  field_simp

/-- Differentiating the Wronskian eliminates the interior Green-kernel term. -/
theorem green_identity (p p' p'' g g' : ℝ) :
    p''*g+p'*g'-(p'*g'+p*(g/4))=(p''-p/4)*g := by ring

end WidthRegularity

#print axioms WidthRegularity.atom_floor
#print axioms WidthRegularity.diagonal
#print axioms WidthRegularity.square_weight
#print axioms WidthRegularity.critical_shell
#print axioms WidthRegularity.green_identity
