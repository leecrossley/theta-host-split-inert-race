import Formalism.LocalRingCardinality.SpanPowerMembership.Primes

open Formalism.BasicNotions.Examples₂
namespace SecondHost

/-- Horner evaluation on base-256 digits, least significant digit first. -/
private def digitPower {n : ℕ} (a : ZMod n) : List ℕ → ZMod n
  | [] => 1
  | d :: ds => (digitPower a ds) ^ 256 * a ^ d

private theorem digitPower_eq {n : ℕ} (a : ZMod n) (ds : List ℕ) :
    digitPower a ds = a ^ Nat.ofDigits 256 ds := by
  induction ds with
  | nil => simp [digitPower]
  | cons d ds ih =>
    rw [digitPower, ih, Nat.ofDigits_cons, pow_add, mul_comm 256, pow_mul, mul_comm]

private def reducedPower (n a e : ℕ) : ZMod n :=
  digitPower (a : ZMod n) (Nat.digits 256 e)

private theorem reducedPower_eq (n a e : ℕ) : reducedPower n a e = (a : ZMod n) ^ e := by
  rw [reducedPower, digitPower_eq, Nat.ofDigits_digits]

theorem typePrimeA : Nat.Prime 343057911569 := by
  refine lucas_list 343057911569 3 [2, 547, 727, 53917] 4 (by decide +kernel) ?_ ?_ ?_
  · intro p hp
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hp
    rcases hp with rfl | rfl | rfl | rfl <;> first | assumption | norm_num
  · rw [← reducedPower_eq]
    decide +kernel
  · intro p hp
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hp
    rcases hp with rfl | rfl | rfl | rfl <;> rw [← reducedPower_eq] <;> decide +kernel

theorem typePrimeB : Nat.Prime 343070636977 := by
  refine lucas_list 343070636977 5 [2, 3, 547, 727, 1997] 4 (by decide +kernel) ?_ ?_ ?_
  · intro p hp
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hp
    rcases hp with rfl | rfl | rfl | rfl | rfl <;> first | assumption | norm_num
  · rw [← reducedPower_eq]
    decide +kernel
  · intro p hp
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hp
    rcases hp with rfl | rfl | rfl | rfl | rfl <;> rw [← reducedPower_eq] <;> decide +kernel

set_option maxRecDepth 100000 in
theorem normPrime : Nat.Prime 9492694782028222020118106107504940049676325753802488976527064603902501412216959612918260301437431961394880032513 := by
  have hg := typePrimeA
  have hh := typePrimeB
  refine lucas_list 9492694782028222020118106107504940049676325753802488976527064603902501412216959612918260301437431961394880032513 3 [2, 41, 547, 727, 343057911569, 343070636977] 8 (by decide +kernel) ?_ ?_ ?_
  · intro p hp
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hp
    rcases hp with rfl | rfl | rfl | rfl | rfl | rfl <;> first | assumption | norm_num
  · rw [← reducedPower_eq]
    decide +kernel
  · intro p hp
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hp
    rcases hp with rfl | rfl | rfl | rfl | rfl | rfl <;> rw [← reducedPower_eq] <;> decide +kernel

/-- The two finite exceptional-prime sets for the type primes are disjoint. -/
theorem disjointExceptional : Disjoint (Formalism.CyclotomicNorms.BB 9 343057911569) (Formalism.CyclotomicNorms.BB 9 343070636977) := by
  rw [Set.disjoint_left]
  rintro l (rfl | ⟨d, h1, h9, h⟩) (h' | ⟨d', h1', h9', h'⟩)
  · exact absurd (Set.mem_singleton_iff.mp h') (by norm_num)
  · exact absurd (exc_dvd (by norm_num) h1' h9' h') (by decide +kernel)
  · rw [Set.mem_singleton_iff.mp h'] at h
    exact absurd (exc_dvd (by norm_num) h1 h9 h) (by decide +kernel)
  · have hg : Nat.gcd (prodPowSubOne 343057911569) (prodPowSubOne 343070636977) ∣ [2, 3, 5, 7, 13, 17, 19, 31, 547, 727].prod ^ 43 := by decide +kernel
    have hm := mem_of_dvd_prod_pow h.1 _ _ (by decide +kernel)
      ((Nat.dvd_gcd (exc_dvd (by norm_num) h1 h9 h) (exc_dvd (by norm_num) h1' h9' h')).trans hg)
    have key : ∀ l ∈ [2, 3, 5, 7, 13, 17, 19, 31, 547, 727], ∀ d ∈ List.range 10, 1 ≤ d → ∃ k ∈ List.range 20,
        ¬ l ^ (k + 1) ∣ 343057911569 ^ d - 1 ∧ l ^ (4 * k) ≤ 343057911569 ^ d := by decide +kernel
    obtain ⟨k, -, hk1, hk2⟩ := key l hm d (List.mem_range.mpr (by omega)) h1
    exact not_exc hk1 hk2 h


end SecondHost
