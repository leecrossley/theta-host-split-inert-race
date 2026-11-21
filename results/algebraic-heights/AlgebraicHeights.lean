import Mathlib

namespace AlgebraicHeights

/-- A protected logarithm stays outside the span after nonzero scaling. -/
theorem protected_log {K V : Type*} [Field K] [AddCommGroup V] [Module K V]
    (W : Submodule K V) {v : V} (hv : v ∉ W) {c : K} (hc : c ≠ 0) :
    c • v ∉ W := by
  intro h
  have back := W.smul_mem c⁻¹ h
  apply hv
  simpa [smul_smul, hc] using back

/-- Two distinct zeros force an affine piece to vanish identically. -/
theorem affine_zeros {a b u v : ℝ} (huv : u ≠ v)
    (hu : a*u+b = 0) (hv : a*v+b = 0) : a = 0 ∧ b = 0 := by
  have h : a*(u-v) = 0 := by linarith only [hu, hv]
  have ha := (mul_eq_zero.mp h).resolve_right (sub_ne_zero.mpr huv)
  exact ⟨ha, by simpa [ha] using hu⟩

/-- The coefficients for inserting the intermediate knots at b and b^3. -/
theorem square_refinement (b x : ℝ) :
    (x^2-1)*(x^2-b^2) = (x-1)*(x-b)*(b+(b+1)*x+x^2) := by ring

/-- Clearing the positive square-root coordinate yields a quadratic. -/
theorem threshold_cell {u : ℝ} (hu : 0 < u) (a b t : ℝ) :
    a*u+b/u = t ↔ a*u^2-t*u+b = 0 := by
  constructor
  · intro h
    have cleared := congrArg (fun s : ℝ => s*u) h
    field_simp [ne_of_gt hu] at cleared
    nlinarith [cleared]
  · intro h
    apply (mul_left_inj' (ne_of_gt hu)).mp
    field_simp [ne_of_gt hu]
    nlinarith [h]

/-- Sliding a strict winning interval against an allowed later-start interval. -/
theorem eligible {l u h d : ℝ} (hlu : l < u) (hhd : h < d) :
    (∃ y, h < y ∧ y ≤ d ∧ l < y ∧ y < u) ↔ h < u ∧ l < d := by
  constructor
  · rintro ⟨y, hhy, hyd, hly, hyu⟩
    exact ⟨lt_trans hhy hyu, lt_of_lt_of_le hly hyd⟩
  · rintro ⟨hhu, hld⟩
    have gap : max h l < min u d := lt_min (max_lt hhu hlu) (max_lt hhd hld)
    obtain ⟨y, hy, hy'⟩ := exists_between gap
    exact ⟨y, lt_of_le_of_lt (le_max_left _ _) hy,
      le_trans hy'.le (min_le_right _ _), lt_of_le_of_lt (le_max_right _ _) hy,
      lt_of_lt_of_le hy' (min_le_left _ _)⟩

/-- A common positive threshold for two different quadratic cells forces
the homogeneous relation between their coefficient forms. -/
theorem simultaneous_cells {a b c d w s : ℝ}
    (h₁ : a*w^2+b = s*w) (h₂ : c*w^2+d = s*w) :
    s^2*(a-c)*(d-b) = (a*d-c*b)^2 := by
  have gap : (a-c)*w^2 = d-b := by linear_combination h₁ - h₂
  have det : a*d-c*b = s*w*(a-c) := by linear_combination a*h₂ - c*h₁
  rw [← gap, det]
  ring


/-- The inverse-height discriminant after clearing the positive square root. -/
theorem inverse_discriminant {u v s y : ℝ}
    (h : (u*y+v)^2 = s^2*y) :
    (2*u^2*y+2*u*v-s^2)^2 = s^2*(s^2-4*u*v) := by
  linear_combination 4*u^2*h

/-- Two exact nonzero-target observations recover the same squared normalizer. -/
theorem common_normalizer {u v p q x y t s a : ℝ}
    (hx : (u*x+v)^2 = t^2*a^2*x) (hy : (p*y+q)^2 = s^2*a^2*y) :
    s^2*y*(u*x+v)^2 = t^2*x*(p*y+q)^2 := by
  rw [hx, hy]
  ring

/-- A repeated root of the cleared cell equation cannot change its sign. -/
theorem repeated_cell {u v t w : ℝ} (hw : u*w^2-t*w+v=0)
    (hd : 2*u*w=t) (z : ℝ) :
    u*z^2-t*z+v = u*(z-w)^2 := by
  linear_combination hw + z*hd - w*hd

/-- Squaring yields paired roots of the height equation in one prime cell. -/
theorem paired_roots {u v t x y : ℝ} (hxy : x ≠ y)
    (hx : (u*x+v)^2=t^2*x) (hy : (u*y+v)^2=t^2*y) :
    u^2*x*y = v^2 := by
  apply (mul_right_cancel₀ (sub_ne_zero.mpr hxy) : _)
  linear_combination y*hx - x*hy


/-- In log height the squared value minus four times the squared derivative
is constant on a prime cell. The second bracket is twice that derivative. -/
theorem cell_energy (u v w : ℝ) (hw : w ≠ 0) :
    (u*w+v/w)^2 - (u*w-v/w)^2 = 4*u*v := by
  field_simp
  ring

/-- The change of this constant at a knot with a continuous score value. -/
theorem energy_jump (value before after : ℝ) :
    (value^2-4*after^2) - (value^2-4*before^2) =
      -4*(after-before)*(after+before) := by ring

/-- A simple zero of a polynomial prevents it from being a polynomial square. -/
theorem simple_root_obstruction {K : Type*} [Field K]
    (f : Polynomial K) (z : K) (hz : f.eval z = 0)
    (hs : f.derivative.eval z ≠ 0) : ¬ ∃ g : Polynomial K, f = g^2 := by
  rintro ⟨g, rfl⟩
  have hg : g.eval z = 0 := by
    rw [Polynomial.eval_pow] at hz
    by_contra hn
    exact (pow_ne_zero 2 hn) hz
  apply hs
  simp [Polynomial.derivative_pow, hg]

/-- Relations among scaled samples give relations between logarithmic forms. -/
theorem ratio_relation {ι : Type*} [Fintype ι] (x z : ℂ)
    (u v target coeff : ι → ℂ)
    (h : ∀ i, u i * x + v i = target i * z)
    (hz : ∑ i, coeff i * target i = 0) :
    x * (∑ i, coeff i * u i) = -(∑ i, coeff i * v i) := by
  have hzero : (∑ i, coeff i * (u i * x + v i)) = 0 := by
    calc
      _ = ∑ i, (coeff i * target i) * z := by
        apply Finset.sum_congr rfl
        intro i hi
        rw [h i, mul_assoc]
      _ = (∑ i, coeff i * target i) * z := (Finset.sum_mul ..).symm
      _ = 0 := by rw [hz, zero_mul]
  simp_rw [mul_add, ← mul_assoc] at hzero
  rw [Finset.sum_add_distrib, ← Finset.sum_mul] at hzero
  rw [mul_comm]
  exact eq_neg_of_add_eq_zero_left hzero

/-- Two different scaled observations of one quadratic cell isolate its slope. -/
theorem recover_scaled_slope {p q u v a t s : ℝ}
    (hp : p^2*u+v=p*t*a) (hq : q^2*u+v=q*s*a) :
    (p^2-q^2)*u=(p*t-q*s)*a := by linear_combination hp - hq

/-- A factor-four-thirds gap puts the next outer prime beyond every earlier norm. -/
theorem fresh_prime_order {y z p : ℝ} (h : 4*y ≤ 3*z) (hp : 6*z < p) :
    8*y < p := by linarith

/-- Clearing a signed equality at algebraically related starts recovers its height. -/
theorem dilation_elimination {u v p q y : ℝ} (h : u*y+v=p*y+q) :
    (u-p)*y=q-v ∧ (u-p)*(u*y+v)=u*q-p*v := by
  constructor
  · linear_combination h
  · linear_combination u*h

/-- The squared common amplitude is obtained without adjoining the height's square root. -/
theorem dilation_amplitude {u v p q y w : ℝ}
    (h : u*y+v=p*y+q) (hw : w*y=(u*y+v)^2) :
    (u-p)*(q-v)*w = (u*q-p*v)^2 := by
  obtain ⟨hy,ht⟩ := dilation_elimination h
  calc
    _ = (u-p)^2*(w*y) := by rw [← hy]; ring
    _ = ((u-p)*(u*y+v))^2 := by rw [hw]; ring
    _ = _ := by rw [ht]

/-- A descending late branch and an ascending early branch sweep a long ratio interval. -/
theorem passage_length {v₁ v₂ h a b : ℝ} (hv : v₂ ≤ v₁) (hh : h ≤ v₂)
    (ha : 0 ≤ a) (hab : b ≤ a) :
    h*(a-b) ≤ v₁*a-v₂*b := by
  nlinarith [mul_nonneg (sub_nonneg.mpr hv) ha,
    mul_nonneg (sub_nonneg.mpr hh) (sub_nonneg.mpr hab)]

/-- Between two neighbours of determinant one, every fraction uses their summed denominator or more. -/
theorem farey_denominator {p q l r a b : ℤ}
    (hd : r*a-l*b=1) (ha : 0 ≤ a) (hb : 0 ≤ b)
    (hl : 0 < p*a-q*l) (hr : 0 < q*r-p*b) : a+b ≤ q := by
  have xl : 1 ≤ p*a-q*l := by omega
  have xr : 1 ≤ q*r-p*b := by omega
  have hz : b*(p*a-q*l)+a*(q*r-p*b)=q := by
    calc
      _ = q*(r*a-l*b) := by ring
      _ = q := by rw [hd, mul_one]
  nlinarith [mul_nonneg hb (sub_nonneg.mpr xl), mul_nonneg ha (sub_nonneg.mpr xr)]

end AlgebraicHeights

#print axioms AlgebraicHeights.protected_log
#print axioms AlgebraicHeights.affine_zeros
#print axioms AlgebraicHeights.square_refinement
#print axioms AlgebraicHeights.threshold_cell
#print axioms AlgebraicHeights.eligible
#print axioms AlgebraicHeights.simultaneous_cells
#print axioms AlgebraicHeights.inverse_discriminant
#print axioms AlgebraicHeights.common_normalizer
#print axioms AlgebraicHeights.repeated_cell
#print axioms AlgebraicHeights.paired_roots
#print axioms AlgebraicHeights.cell_energy
#print axioms AlgebraicHeights.energy_jump
#print axioms AlgebraicHeights.simple_root_obstruction
#print axioms AlgebraicHeights.ratio_relation
#print axioms AlgebraicHeights.recover_scaled_slope
#print axioms AlgebraicHeights.fresh_prime_order
#print axioms AlgebraicHeights.dilation_elimination
#print axioms AlgebraicHeights.dilation_amplitude
#print axioms AlgebraicHeights.passage_length
#print axioms AlgebraicHeights.farey_denominator
