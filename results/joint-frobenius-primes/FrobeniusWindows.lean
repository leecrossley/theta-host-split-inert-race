import Mathlib

namespace FrobeniusWindows

/-- A Proth congruence gives a high power of two in the multiplicative group
of every prime divisor. The size condition then forces primality. -/
theorem proth {N k n : ℕ} (hN : 1 < N) (hodd : Odd N)
    (hsize : N ≤ (2^(n+1))^2) (a : ZMod N)
    (hcert : a^(k*2^n) = -1) : N.Prime := by
  by_contra composite
  let d := N.minFac
  have hd : d.Prime := Nat.minFac_prime (by omega)
  have hdvd : d ∣ N := Nat.minFac_dvd N
  have hne : d ≠ 2 := by
    intro h
    exact hodd.not_two_dvd_nat (h ▸ hdvd)
  have htwo := hd.two_le
  let _ : Fact d.Prime := ⟨hd⟩
  let _ : Fact (2 < d) := ⟨by omega⟩
  let f := ZMod.castHom hdvd (ZMod d)
  let z := (f a)^k
  have half : z^(2^n) = -1 := by
    have h := congrArg f hcert
    simpa only [map_pow, map_neg, map_one, z, pow_mul] using h
  have full : z^(2^(n+1)) = 1 := by
    rw [pow_succ, pow_mul, half]
    norm_num
  have exactOrder : orderOf z = 2^(n+1) :=
    orderOf_eq_prime_pow (by rw [half]; exact ZMod.neg_one_ne_one) full
  have nz : z ≠ 0 := by
    intro hz
    rw [hz, zero_pow (by positivity : 2^n ≠ 0)] at half
    exact (neg_ne_zero.mpr one_ne_zero) half.symm
  have ordDvd := ZMod.orderOf_dvd_card_sub_one nz
  rw [exactOrder] at ordDvd
  have small := Nat.le_of_dvd (by omega : 0 < d-1) ordDvd
  have strict : 2^(n+1) < d := by omega
  have sq := Nat.pow_lt_pow_left strict (by decide : 2 ≠ 0)
  have upper : d^2 ≤ N := Nat.minFac_sq_le_self (by omega) composite
  omega

/-- A next prime q is missing precisely up to the open right endpoint q/C. -/
theorem deficient_gap {L p q C : ℝ} (hc : 0 < C) (hp : p ≤ L) (hq : L < q) :
    (¬ (L < q ∧ q ≤ C*L)) ↔ p ≤ L ∧ L < q/C := by
  rw [lt_div_iff₀ hc]
  constructor
  · intro h
    exact ⟨hp, by simpa [mul_comm] using lt_of_not_ge (fun t => h ⟨hq, t⟩)⟩
  · rintro ⟨_, h⟩ ⟨_, t⟩
    linarith only [h, t, mul_comm C L]

/-- A consecutive-prime ratio bounds the waiting factor at every real start in the gap. -/
theorem cover_between {L p q : ℝ} (hp : p ≤ L) (hq : L < q) (gap : 3*q ≤ 4*p) :
    L < q ∧ q ≤ 4*L/3 := by
  exact ⟨hq, by linarith⟩

/-- Passing a two-sided endpoint error to the mass in a moving interval. -/
theorem interval_mass {A B x y e f : ℝ}
    (hA : |A-x/64| ≤ e) (hB : |B-y/64| ≤ f) :
    |(B-A)-(y-x)/64| ≤ e+f := by
  rw [abs_le] at *
  constructor <;> linarith [hA.1, hA.2, hB.1, hB.2]

end FrobeniusWindows

#print axioms FrobeniusWindows.proth
#print axioms FrobeniusWindows.deficient_gap
#print axioms FrobeniusWindows.cover_between
#print axioms FrobeniusWindows.interval_mass
