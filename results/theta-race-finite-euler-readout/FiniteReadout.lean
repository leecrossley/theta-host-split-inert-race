import Mathlib

namespace FiniteReadout

def hinge (k n : ℤ) : ℝ := max ((k:ℝ)-(n:ℝ)) 0

theorem hinge_stencil (k n : ℤ) :
    hinge (k+1) n - 2*hinge k n + hinge (k-1) n = if k=n then 1 else 0 := by
  unfold hinge
  push_cast
  rcases lt_trichotomy k n with h | h | h
  · have hh : ((k:ℝ)+1 ≤ n) := by exact_mod_cast (show k+1 ≤ n by omega)
    rw [if_neg (ne_of_lt h), max_eq_right (by linarith), max_eq_right (by linarith),
      max_eq_right (by linarith)]
    ring
  · subst n
    simp
  · have hh : ((n:ℝ)+1 ≤ k) := by exact_mod_cast (show n+1 ≤ k by omega)
    rw [if_neg (ne_of_gt h), max_eq_left (by linarith), max_eq_left (by linarith),
      max_eq_left (by linarith)]
    ring

noncomputable def kernel (k n : ℤ) : ℝ :=
  (hinge k n - (3/2)*hinge k (2*n) + (1/2)*hinge k (4*n))/(n:ℝ)

theorem kernel_stencil (k n : ℤ) :
    kernel (k+1) n -2*kernel k n +kernel (k-1) n =
    ((if k=n then 1 else 0) - (3/2)*(if k=2*n then 1 else 0)
    + (1/2)*(if k=4*n then 1 else 0))/(n:ℝ) := by
  rw [← hinge_stencil k n, ← hinge_stencil k (2*n), ← hinge_stencil k (4*n)]
  unfold kernel
  ring

theorem odd_isolation (k n : ℤ) (h : k % 2 = 1) :
    kernel (k+1) n -2*kernel k n +kernel (k-1) n =
      (if k=n then 1 else 0)/(n:ℝ) := by
  rw [kernel_stencil, if_neg (show k ≠ 2*n by omega), if_neg (show k ≠ 4*n by omega)]
  ring

theorem dyadic_inverse_step (j : ℕ) :
    (2:ℤ)^(j+3)-1 = 3*(2^(j+2)-1)-2*(2^(j+1)-1) := by
  ring

theorem noise_bound {x y z e : ℝ} (hx : |x| ≤ e) (hy : |y| ≤ e) (hz : |z| ≤ e) :
    |x-2*y+z| ≤ 4*e := by
  rw [abs_le] at *
  constructor <;> linarith

theorem signed_rounding {x : ℝ} {m n : ℤ} (hm : |x-8*m| < 4) (hn : |x-8*n| < 4) :
    m=n := by
  rw [abs_lt] at *
  by_contra h
  rcases lt_or_gt_of_ne h with h | h
  · have ht : (m:ℝ)+1 ≤ n := by exact_mod_cast (show m+1 ≤ n by omega)
    linarith
  · have ht : (n:ℝ)+1 ≤ m := by exact_mod_cast (show n+1 ≤ m by omega)
    linarith

end FiniteReadout

#print axioms FiniteReadout.hinge_stencil
#print axioms FiniteReadout.kernel_stencil
#print axioms FiniteReadout.odd_isolation
#print axioms FiniteReadout.dyadic_inverse_step
#print axioms FiniteReadout.noise_bound
#print axioms FiniteReadout.signed_rounding
