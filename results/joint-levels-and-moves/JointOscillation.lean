import Mathlib

namespace JointOscillation

open scoped BigOperators

theorem projected_norm {g : ℝ} (hg : 0 ≤ g) (u v k : ℝ) :
    Real.sqrt ((g*u)^2 + (g*k*v)^2) = g * Real.sqrt (u^2 + (k*v)^2) := by
  rw [show (g*u)^2 + (g*k*v)^2 = g^2*(u^2+(k*v)^2) by ring,
    Real.sqrt_mul (sq_nonneg g), Real.sqrt_sq hg]

theorem disk_support {x y : ℝ} (h : x^2+y^2 ≤ 1) (u v : ℝ) :
    u*x+v*y ≤ Real.sqrt (u^2+v^2) := by
  apply Real.le_sqrt_of_sq_le
  have h1 := mul_nonneg (add_nonneg (sq_nonneg u) (sq_nonneg v)) (sub_nonneg.mpr h)
  nlinarith only [h1, sq_nonneg (u*y-v*x)]

theorem ellipse_support {ι : Type*} [Fintype ι] (g k : ι → ℝ)
    (hg : ∀ i, 0 ≤ g i) {x y : ℝ} (h : x^2+y^2 ≤ 1) (u v : ℝ) :
    (∑ i, g i)*u*x + (∑ i, g i*k i)*v*y ≤
      ∑ i, g i*Real.sqrt (u^2+(k i*v)^2) := by
  calc
    _ = ∑ i, (g i*(u*x+(k i*v)*y)) := by
      simp only [Finset.sum_mul, mul_add, Finset.sum_add_distrib]
      congr 1 <;> apply Finset.sum_congr rfl <;> intro i hi <;> ring
    _ ≤ _ := Finset.sum_le_sum fun i _ =>
      mul_le_mul_of_nonneg_left (disk_support h u (k i*v)) (hg i)

theorem projection_gap {x y a b t : ℝ} (hx : a < x) (hy : b < y)
    (ht : t ≤ a+b) : t < x+y := by linarith

theorem one_sided_pair {x y t : ℝ} (h : t < x-y) : -t > y-x := by linarith

theorem ellipse_axes {a b c d x y : ℝ} (ha : 0 ≤ a) (hb : 0 ≤ b)
    (hac : a ≤ c) (hbd : b ≤ d) :
    (a*x)^2+(b*y)^2 ≤ (c*x)^2+(d*y)^2 := by
  nlinarith only [mul_nonneg (mul_nonneg (sub_nonneg.mpr hac) (by linarith : 0 ≤ c+a)) (sq_nonneg x),
    mul_nonneg (mul_nonneg (sub_nonneg.mpr hbd) (by linarith : 0 ≤ d+b)) (sq_nonneg y)]

end JointOscillation

#print axioms JointOscillation.projected_norm
#print axioms JointOscillation.disk_support
#print axioms JointOscillation.ellipse_support
#print axioms JointOscillation.projection_gap
#print axioms JointOscillation.one_sided_pair
#print axioms JointOscillation.ellipse_axes
