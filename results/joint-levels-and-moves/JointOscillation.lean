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


/-- The fixed point of the simplex map with strict inward boundary displacements
has five positive weights and zero displacement in every coordinate. -/
theorem simplex_fixed_displacement (w d : Fin 5 → ℝ) (s : ℝ)
    (hw : ∀ i, 0 ≤ w i) (hm : ∑ i, w i = 1) (hd : ∑ i, d i = 0)
    (hs : s = ∑ i, max 0 (w i + d i))
    (hf : ∀ i, max 0 (w i + d i) = w i * s)
    (hb : ∀ i, w i = 0 → 0 < d i) :
    (∀ i, 0 < w i) ∧ ∀ i, d i = 0 := by
  have hmass : ∑ i, (w i + d i) = 1 := by
    rw [Finset.sum_add_distrib, hm, hd, add_zero]
  have hscale : 1 ≤ s := by
    rw [hs, ← hmass]
    exact Finset.sum_le_sum (fun i _ => le_max_right _ _)
  have hpos : ∀ i, 0 < w i := by
    intro i
    by_contra! hi
    have hz := le_antisymm hi (hw i)
    have hle := le_max_right 0 (w i + d i)
    rw [hf i, hz, zero_mul, zero_add] at hle
    exact (not_lt_of_ge hle) (hb i hz)
  have hunclip (i) : w i + d i = w i * s := by
    have hmax := hf i
    have hgt : 0 < w i * s := mul_pos (hpos i) (by linarith)
    by_cases hn : 0 ≤ w i + d i
    · rwa [max_eq_right hn] at hmax
    · rw [max_eq_left (le_of_not_ge hn)] at hmax
      linarith
  have hone : s = 1 := by
    calc
      s = ∑ i, w i * s := by rw [← Finset.sum_mul, hm, one_mul]
      _ = ∑ i, (w i + d i) := Finset.sum_congr rfl (fun i _ => (hunclip i).symm)
      _ = 1 := hmass
  refine ⟨hpos, fun i => ?_⟩
  have h := hunclip i
  rw [hone, mul_one] at h
  linarith

/-- An upper bound on all three face contributions also bounds their common average. -/
theorem three_point_face {u v w x y z m : ℝ} (hu : 0 ≤ u) (hv : 0 ≤ v)
    (hw : 0 ≤ w) (hm : u+v+w = 1) (hx : x ≤ m) (hy : y ≤ m) (hz : z ≤ m) :
    u*x+v*y+w*z ≤ m := by
  calc
    u*x+v*y+w*z ≤ u*m+v*m+w*m := by gcongr
    _ = m := by rw [← add_mul, ← add_mul, hm, one_mul]


/-- Squaring the positive mass equation eliminates both square roots. -/
theorem quartic_elimination {a b c : ℝ} (mass : a+b+c = 1) :
    ((1-a)^2-b^2-c^2)^2-4*b^2*c^2 = 0 := by
  have h : 1-a = b+c := by linarith
  rw [h]
  ring

/-- Once the anchor parameter is known, recovering the other weights needs no new radical. -/
theorem recover_weight {a b c : ℝ} (mass : a+b+c = 1) (remaining : 0 < 1-a) :
    ((1-a)^2+b^2-c^2)/(2*(1-a)) = b := by
  apply (div_eq_iff (ne_of_gt (mul_pos (by norm_num) remaining))).mpr
  have h : 1-a = b+c := by linarith
  rw [h]
  ring

/-- On the positive sheet the sign discarded by squaring is prescribed. -/
theorem positive_mass_sign {h b c : ℝ} (hm : h+b+c=1) (hb : 0<b) (hc : 0<c) :
    0 < (1-h)^2-b^2-c^2 := by
  have h : 1-h = b+c := by linarith
  rw [h]
  nlinarith [mul_pos hb hc]

/-- The group deduction for a connected four-sheeted cover with simple branch generators. -/
theorem four_sheet_group (S : Set (Equiv.Perm (Fin 4)))
    (hs : ∀ σ ∈ S, σ.IsSwap) [MulAction.IsPretransitive (Subgroup.closure S) (Fin 4)] :
    Subgroup.closure S = ⊤ :=
  closure_of_isSwap_of_isPretransitive hs

/-- The discriminant of the quadratic left after projection from a point of a cubic. -/
theorem projection_discriminant {a b c x : ℝ} (hx : a*x^2+b*x+c=0) :
    (2*a*x+b)^2=b^2-4*a*c := by
  nlinarith [show (2*a*x+b)^2-(b^2-4*a*c)=4*a*(a*x^2+b*x+c) by ring,
    show 4*a*(a*x^2+b*x+c)=0 by rw [hx, mul_zero]]

/-- A bounded branch on a continuous normalized pencil supplies a root on the target fibre. -/
theorem specialized_root {f : ℝ × ℝ → ℝ} (hf : Continuous f) {u v : ℕ → ℝ} {s : ℝ}
    (hu : Filter.Tendsto u Filter.atTop (nhds 0))
    (hv : Filter.Tendsto v Filter.atTop (nhds s))
    (hz : ∀ n, f (u n,v n)=0) : f (0,s)=0 := by
  have ht : Filter.Tendsto (fun n => f (u n,v n)) Filter.atTop (nhds (f (0,s))) :=
    hf.continuousAt.tendsto.comp (hu.prodMk_nhds hv)
  have he : (fun n => f (u n,v n)) = fun _ => 0 := funext hz
  rw [he] at ht
  exact tendsto_nhds_unique ht tendsto_const_nhds

/-- Distinct four-variable exponents of individual degree at most four stay distinct. -/
theorem four_digit_exponents {a b c d i j k l : ℕ}
    (ha : a<5) (hb : b<5) (hc : c<5) (hi : i<5) (hj : j<5) (hk : k<5)
    (h : a+5*b+25*c+125*d=i+5*j+25*k+125*l) :
    a=i ∧ b=j ∧ c=k ∧ d=l := by omega

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

#print axioms JointOscillation.simplex_fixed_displacement
#print axioms JointOscillation.three_point_face

#print axioms JointOscillation.quartic_elimination
#print axioms JointOscillation.recover_weight

#print axioms JointOscillation.positive_mass_sign
#print axioms JointOscillation.four_sheet_group
#print axioms JointOscillation.projection_discriminant

#print axioms JointOscillation.specialized_root
#print axioms JointOscillation.four_digit_exponents
