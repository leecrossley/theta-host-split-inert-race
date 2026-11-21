import Mathlib

namespace FilteredDeadline

/-- An omitted zero in a covered interval would exceed the audited Poisson remainder. -/
theorem no_missing {d R : ℝ} (hd : |d| ≤ 2) (hR : R < 4/25) :
    ¬ (1 / ((1/2:ℝ)^2+d^2) ≤ R) := by
  have hsq := (sq_le_sq.mpr (show |d| ≤ |(2:ℝ)| by simpa using hd))
  have hden : 0 < (1/2:ℝ)^2+d^2 := by positivity
  intro h
  have h' := (div_le_iff₀ hden).mp h
  have h'' : R*((1/2:ℝ)^2+d^2) < (4/25)*((1/2:ℝ)^2+d^2) :=
    mul_lt_mul_of_pos_right hR hden
  nlinarith

/-- A phase-aligned carrier controls both signs after all error allowances. -/
theorem carrier_budget {signed carrier leakage missing err threshold : ℝ}
    (h : carrier-leakage-missing-err ≤ signed)
    (margin : threshold < carrier-leakage-missing-err) : threshold < signed :=
  margin.trans_le h

/-- A nonnegative unit-mass average above a threshold requires some point above it. -/
theorem witness_of_average {α : Type*} {f : α → ℝ} {support : Set α}
    {value threshold : ℝ} (large : threshold < value)
    (upper : (∀ x ∈ support, f x ≤ threshold) → value ≤ threshold) :
    ∃ x ∈ support, threshold < f x := by
  by_contra! hn
  exact (not_le.mpr large) (upper hn)

/-- Two support radii and one period fit inside the admissible interval of log centres. -/
theorem centre_room {left right radius period : ℝ} (h : 2*radius+period < right-left) :
    left+radius+period < right-radius := by linarith

/-- Any point inside a translated support remains between the outer endpoints. -/
theorem support_fits {x centre radius left right : ℝ}
    (h₁ : left+radius < centre) (h₂ : centre < right-radius)
    (hx : |x-centre| ≤ radius) : left < x ∧ x < right := by
  obtain ⟨lo,hi⟩ := abs_le.mp hx
  constructor <;> linarith

end FilteredDeadline

#print axioms FilteredDeadline.no_missing
#print axioms FilteredDeadline.carrier_budget
#print axioms FilteredDeadline.witness_of_average
#print axioms FilteredDeadline.centre_room
#print axioms FilteredDeadline.support_fits
