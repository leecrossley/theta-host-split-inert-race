import Mathlib

namespace ConductorHorizon

theorem hinge_fraction (k n : ℤ) {d : ℝ} (hd : 0 < d) (h1 : d ≤ 1) :
    max ((k:ℝ)+d-n) 0 -2*max ((k:ℝ)-n) 0 +max ((k:ℝ)-d-n) 0 =
      if k=n then d else 0 := by
  rcases lt_trichotomy k n with h | h | h
  · have hh : (k:ℝ)+1 ≤ n := by exact_mod_cast (show k+1 ≤ n by omega)
    rw [if_neg (ne_of_lt h), max_eq_right (by linarith), max_eq_right (by linarith),
        max_eq_right (by linarith)]
    ring
  · subst n
    rw [if_pos rfl, max_eq_left (by linarith), sub_self, max_self, max_eq_right (by linarith)]
    ring
  · have hh : (n:ℝ)+1 ≤ k := by exact_mod_cast (show n+1 ≤ k by omega)
    rw [if_neg (ne_of_gt h), max_eq_left (by linarith), max_eq_left (by linarith),
        max_eq_left (by linarith)]
    ring

noncomputable def kernel (x : ℝ) (n : ℤ) : ℝ :=
  (max (x-n) 0 - (3/2)*max (x-(2*n:ℤ)) 0 + (1/2)*max (x-(4*n:ℤ)) 0)/n

theorem centered_jump (k n : ℤ) (hk : k%2=1) {d : ℝ} (hd : 0<d) (h1 : d≤1) :
    kernel ((k:ℝ)+d) n-2*kernel k n+kernel ((k:ℝ)-d) n =
      (if k=n then d else 0)/(n:ℝ) := by
  have ha := hinge_fraction k n hd h1
  have hb := hinge_fraction k (2*n) hd h1
  have hc := hinge_fraction k (4*n) hd h1
  rw [if_neg (show k≠2*n by omega)] at hb
  rw [if_neg (show k≠4*n by omega)] at hc
  unfold kernel
  rw [← ha]
  linear_combination (norm := ring_nf) ((-3/2:ℝ)/(n:ℝ))*hb + ((1/2:ℝ)/(n:ℝ))*hc

theorem strict_cutoff {N Y : ℝ} (h : N/4 < Y) :
    ∃ d : ℝ, 0 < d ∧ d ≤ 1 ∧ (N+d)/4 ≤ Y := by
  refine ⟨min 1 (4*Y-N), lt_min (by norm_num) (by linarith), min_le_left _ _, ?_⟩
  linarith [min_le_right (1:ℝ) (4*Y-N)]

theorem missing_or_present {s : ℝ} (hs : s=1 ∨ s= -1) {f : ℕ}
    (hf : f=2 ∨ f=4) : s^f=1 := by
  rcases hs with h | h <;> subst s <;> rcases hf with h | h <;> subst f <;> norm_num

theorem round_boundary {c e : ℝ} (h : c=0 ∨ c= -8) (he : |e|<4)
    (hobs : |c+e|<4) : c=0 := by
  rcases h with h | h
  · exact h
  · rw [abs_lt] at he hobs
    linarith

theorem larger_conductor {m a : ℕ} (hm : m ≠ 0) (ha : a.Coprime m) (B : ℕ) :
    ∃ p > B, p.Prime ∧ p ≡ a [MOD m] :=
  Nat.forall_exists_prime_gt_and_modEq B hm ha

theorem prime_products_injective {u v : ℕ →₀ ℕ}
    (hu : ∀ p ∈ u.support, Nat.Prime p) (hv : ∀ p ∈ v.support, Nat.Prime p)
    (h : u.prod (· ^ ·) = v.prod (· ^ ·)) : u = v := by
  calc
    u = (u.prod (· ^ ·)).factorization := (Nat.prod_pow_factorization_eq_self hu).symm
    _ = (v.prod (· ^ ·)).factorization := congrArg Nat.factorization h
    _ = v := Nat.prod_pow_factorization_eq_self hv

end ConductorHorizon

#print axioms ConductorHorizon.hinge_fraction
#print axioms ConductorHorizon.centered_jump
#print axioms ConductorHorizon.strict_cutoff
#print axioms ConductorHorizon.missing_or_present
#print axioms ConductorHorizon.round_boundary
#print axioms ConductorHorizon.larger_conductor
#print axioms ConductorHorizon.prime_products_injective
