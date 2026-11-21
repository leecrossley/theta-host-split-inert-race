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

theorem deadline_support (h : ℝ) :
    let L := (h+14084)/16999
    let X := h+7002+L
    let R := 7000+L
    X-R = h+2 ∧ X+R < (17001/16999)*h+14004 := by
  dsimp
  constructor
  · ring
  · linarith

theorem mixed_stop_band {h : ℝ} (hh : 100 ≤ h) :
    let L := (h+14084)/16999
    let X := h+7002+L
    (X+7000)/2-8503*L ≤ -41 := by
  dsimp
  linarith

theorem contrast_margin {x y r : ℝ} (hx : r < x) (hy : y < -r) :
    2*r < x-y := by linarith

theorem cover_endpoint {u x radius b width : ℝ}
    (left : b < x-radius) (right : x+radius-width < b)
    (position : x-radius ≤ u ∧ u ≤ x+radius) : b < u ∧ u < b+width := by
  constructor <;> linarith [position.1,position.2]


/-- Distinct real growth rates remain separable at equally spaced log translations. -/
theorem translated_det {n : ℕ} (d : Fin n → ℝ) (hd : Function.Injective d) :
    (Matrix.vandermonde (fun j => Real.exp (d j))).det ≠ 0 := by
  apply Matrix.det_vandermonde_ne_zero_iff.mpr
  exact Real.exp_injective.comp hd

/-- In the nearly perpendicular case the growing target and the discarded error share the margin. -/
theorem near_perpendicular {support target error growth margin : ℝ}
    (hs : target + margin ≤ support) (he : error + growth < margin) :
    target + growth < support - error := by linarith

/-- In the complementary case the extracted exponential dominates both error and target. -/
theorem transverse {signal error target growth : ℝ}
    (h : target + growth + error < signal) :
    target + growth < signal - error := by linarith

end JointOscillation

#print axioms JointOscillation.projected_norm
#print axioms JointOscillation.disk_support
#print axioms JointOscillation.ellipse_support
#print axioms JointOscillation.projection_gap
#print axioms JointOscillation.one_sided_pair
#print axioms JointOscillation.ellipse_axes
#print axioms JointOscillation.deadline_support
#print axioms JointOscillation.mixed_stop_band
#print axioms JointOscillation.contrast_margin
#print axioms JointOscillation.cover_endpoint

#print axioms JointOscillation.translated_det
#print axioms JointOscillation.near_perpendicular
#print axioms JointOscillation.transverse
