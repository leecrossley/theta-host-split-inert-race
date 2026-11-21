import Mathlib

namespace CurvatureEuler

/-- Nonnegative deficits bound a centered second difference, with equality at alignment. -/
theorem peak_extreme {ι : Type*} [Fintype ι] (w x : ι → ℝ)
    (hw : ∀ i, 0 ≤ w i) (hx : ∀ i, |x i| ≤ 1) :
    |∑ i, w i*x i| ≤ ∑ i, w i := by
  calc
    |∑ i, w i*x i| ≤ ∑ i, |w i*x i| := Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ i, w i := by
      apply Finset.sum_le_sum
      intro i _
      rw [abs_mul, abs_of_nonneg (hw i)]
      simpa using mul_le_mul_of_nonneg_left (hx i) (hw i)

/-- An exact three-point sample on either side of an exponential-spline knot. -/
theorem stencil_jump {l m r leftSlope rightSlope s c : ℝ} (hs : s ≠ 0)
    (hl : l = m*c-2*leftSlope*s) (hr : r = m*c+2*rightSlope*s) :
    (l-2*c*m+r)/(2*s) = rightSlope-leftSlope := by
  rw [hl,hr]
  field_simp
  ring

theorem fourth_roots (p m d : ℝ) :
    let u := p-m
    let v := p+m-2*d
    let w := p+m+2*d
    (w+v+2*u)/4 = p ∧ (w+v-2*u)/4 = m ∧ (w-v)/4 = d := by
  dsimp
  constructor
  · ring
  constructor <;> ring

theorem two_adic (u v : ℝ) :
    let h := v-3*u
    let k := 3*u-3*v
    (-h-k)/2=v ∧ ((-h-k)/2-h)/3=u := by
  dsimp
  constructor <;> ring

/-- A real value determines a unique target integer inside an open half-unit ball. -/
theorem rounding {x : ℝ} {m n : ℤ} (hm : |x-m| < 1/2) (hn : |x-n| < 1/2) :
    m = n := by
  obtain ⟨ml,mu⟩ := abs_lt.mp hm
  obtain ⟨nl,nu⟩ := abs_lt.mp hn
  by_contra h
  rcases lt_or_gt_of_ne h with hmn | hnm
  · have hh : m+1 ≤ n := by omega
    exact (show ¬ ((m:ℝ)+1 ≤ (n:ℝ)) by linarith) (by exact_mod_cast hh)
  · have hh : n+1 ≤ m := by omega
    exact (show ¬ ((n:ℝ)+1 ≤ (m:ℝ)) by linarith) (by exact_mod_cast hh)

end CurvatureEuler

#print axioms CurvatureEuler.peak_extreme
#print axioms CurvatureEuler.stencil_jump
#print axioms CurvatureEuler.fourth_roots
#print axioms CurvatureEuler.two_adic
#print axioms CurvatureEuler.rounding
