import HostPrimes
import Development.AbelianVarietyTheory.Proofs.InertiaEigenvalues.LocalGapNonemptyAndRecordNonempty
import Development.AbelianVarietyTheory.Proofs.DivisorProofs.IsTwistPlacesAboveK
import Development.AutomorphicArithmetic.WeightCongruences.IsRootOfIntegral

/- Concrete input to Alpöge's exampleRecordOfData constructor.
   The prescription for the rank, level and norm generator follows Examples₂. -/
noncomputable section

open NumberField Formalism.FrobeniusOrbits Formalism.TwoAdicCyclotomicDecomposition
open Formalism.BasicNotions Formalism.BasicNotions.Examples₂ Formalism.L2.Host
open scoped Pointwise

namespace SecondHost

theorem existsLayer : ∃ (L : Type) (_ : Field L) (_ : NumberField L) (ζ a b : L),
    IsLayer 4 547 727 ζ a b := Formalism.TwoAdicCyclotomicDecomposition.layerExists 4 547 727

def L : Type := Classical.choose existsLayer

instance layerField : Field L := Classical.choose (Classical.choose_spec existsLayer)

theorem existsNumberField : ∃ (_ : NumberField L) (ζ a b : L), IsLayer 4 547 727 ζ a b :=
  Classical.choose_spec (Classical.choose_spec existsLayer)

instance layerNumberField : NumberField L := Classical.choose existsNumberField

theorem existsGenerators : ∃ (ζ a b : L), IsLayer 4 547 727 ζ a b := by
  obtain ⟨_, h⟩ := existsNumberField
  exact h

def zeta : L := Classical.choose existsGenerators
theorem existsPair : ∃ (a b : L), IsLayer 4 547 727 zeta a b := Classical.choose_spec existsGenerators
def alpha : L := Classical.choose existsPair
theorem existsBeta : ∃ b : L, IsLayer 4 547 727 zeta alpha b := Classical.choose_spec existsPair
def beta : L := Classical.choose existsBeta
theorem layer : IsLayer 4 547 727 zeta alpha beta := Classical.choose_spec existsBeta

theorem hyp : HostHyp 547 727 := by
  refine ⟨by norm_num, by norm_num, by norm_num, by norm_num, ?_⟩
  norm_num

theorem typeNonempty : (cmType zeta alpha beta 727).Nonempty := by
  apply Set.nonempty_of_ncard_ne_zero
  rw [(Formalism.FrobeniusOrbits.Final.cmType 2 547 727 hyp.toHypOdd L zeta alpha beta layer).2.1]
  positivity

def axes : Formalism.HC.HostAx R0.toFrame :=
  Formalism.HC.D7.hostAx R0.toFrame exampleRanks_h_injective hyp (le_refl 2) layer

theorem axesTypeNonempty : axes.cmType.Nonempty := by
  exact typeNonempty

def host : Host R0 where
  toHostAx := axes
  flPS := Formalism.HC.D7.psQ R0.toFrame exampleRanks_h_injective hyp (le_refl 2) layer
  sqq := Formalism.HC.D7.sqq R0.toFrame exampleRanks_h_injective hyp (le_refl 2) layer
  sqq_sq := Formalism.HC.D7.sqq_sq R0.toFrame exampleRanks_h_injective hyp (le_refl 2) layer
  gal_sqq := Formalism.HC.D7.gal_sqq R0.toFrame hyp layer
  QS_exact := Formalism.HC.D7.exact R0.toFrame exampleRanks_h_injective hyp (le_refl 2) layer
  QS_le := (by norm_num : 2 ^ (2 + 2) ≤ 4 ^ 2)
  o_pow := by rw [exampleRanks_o]; norm_num
  q_ne_three := (by norm_num : (547 : ℕ) ≠ 3)
  qd_big := (by norm_num : 4 * 15 * 9 + 1 < 547 ∧ 4 * 15 * 9 + 1 < 727)
  splitBound := 2
  split_small := fun p hp h2 hle => absurd (le_antisymm hle hp.two_le) h2
  small_le := fun l hl h2 hd =>
    absurd ((Nat.prime_dvd_prime_iff_eq hl Nat.prime_two).mp
      (hl.dvd_of_dvd_pow (show l ∣ 2 ^ (2 + 2) from hd))) h2
  φ₀ := axesTypeNonempty.some
  φ₀_mem := axesTypeNonempty.some_mem
  𝔟 := axes.der.principal.choose
  𝔟_spec := ⟨axes.der.principal.choose_spec.1, axes.der.principal.choose_spec.2.1⟩
  𝔟_stable := axes.der.principal.choose_spec.2.2
  levelPrimes := {991, 1231, 343057911569, 343070636977, 9492694782028222020118106107504940049676325753802488976527064603902501412216959612918260301437431961394880032513}
  g := ![343057911569, 343070636977]
  g_prime := by intro i; fin_cases i; exacts [typePrimeA, typePrimeB]
  g_ne := by norm_num
  g_level := by intro i; fin_cases i <;> simp
  pde := 991
  pde' := 1231
  pde_prime := by norm_num
  pde'_prime := by norm_num
  pde_ne := by norm_num
  pde_level := by simp
  pde_mod := by norm_num
  pde_ne_g := by intro i; fin_cases i <;> simp
  level_coprimeS := by
    intro p hp
    show Nat.Coprime p (2 ^ (2 + 2))
    simp only [Finset.mem_insert, Finset.mem_singleton] at hp
    rcases hp with rfl | rfl | rfl | rfl | rfl <;> norm_num
  ccK := ⊤
  ccK_ne_bot := top_ne_bot
  ccK_primitive := fun m hm => top_primitive m hm
  ccK_split := fun P hP hle => absurd (top_le_iff.1 hle) hP.ne_top
  ccK_coprime := fun p _ => by rw [Ideal.absNorm_top]; exact Nat.coprime_one_left p
  ccK_coprimeS := by rw [Ideal.absNorm_top]; exact Nat.coprime_one_left _
  TT := ∅
  T2 := {0}
  u0 := 0
  Lw := 50
  two_le_Lw := by norm_num
  second_set_range := by intro s hs; simp at hs; subst hs; rw [exampleRanks_baseRank]; norm_num
  second_set_no_relation := by
    intro s₁ h₁ s₂ h₂ _ _ h12
    simp at h₁ h₂
    exact absurd (h₁.trans h₂.symm) h12
  CT2 := 0

def levelPrime : ℕ := 9492694782028222020118106107504940049676325753802488976527064603902501412216959612918260301437431961394880032513

theorem levelPrime_prime : levelPrime.Prime := normPrime

theorem levelPrime_mem : levelPrime ∈ host.levelPrimes := by
  simp [host, levelPrime]

def multiplier : ℕ := 3613495791090103094491079215867623653897637253161937616

theorem levelModulus : (host.QS * host.q * host.d * (host.g 0 * host.g 1) ^ 2 : ℕ) = 88134043685124465719294615021161552534088713491754576 := by
  rfl

theorem sqrtSquare : host.sqd ^ 2 = -727 := by
  rw [bE_sq]
  have h : host.d = 727 := rfl
  rw [h]
  norm_num

theorem sqrtIntegral : IsIntegral ℤ (host.sqd) := by
  refine ⟨Polynomial.X ^ 2 + Polynomial.C 727, by monicity!, ?_⟩
  simp [sqrtSquare]

def generator : host.E := 1 + (multiplier : host.E) * host.sqd

theorem generatorIntegral : IsIntegral ℤ generator := by
  unfold generator
  refine IsIntegral.add isIntegral_one (IsIntegral.mul ?_ sqrtIntegral)
  have h := isIntegral_algebraMap (R := ℤ) (A := host.E) (x := (multiplier : ℤ))
  simpa using h

theorem generatorMem : generator ∈ host.K := by
  have hb : host.sqd ∈ host.K := IntermediateField.mem_adjoin_simple_self ℚ _
  have h1 : (1 : host.E) ∈ host.K := one_mem _
  have he : (multiplier : host.E) ∈ host.K := by simp
  exact add_mem h1 (mul_mem he hb)

theorem generatorNorm : generator * host.conj generator = (levelPrime : host.E) := by
  unfold generator
  rw [map_add, map_one, map_mul, map_natCast, conj_bE]
  have h := sqrtSquare
  unfold levelPrime multiplier
  push_cast
  linear_combination (-(3613495791090103094491079215867623653897637253161937616 : host.E) ^ 2) * h

def integerGenerator : 𝓞 host.E := ⟨generator, generatorIntegral⟩

def rootInteger : 𝓞 host.E := ⟨host.sqd, sqrtIntegral⟩

theorem generatorCongruence : integerGenerator - 1 ∈ Ideal.span
    {((host.QS * host.q * host.d * (host.g 0 * host.g 1) ^ 2 : ℕ) : 𝓞 host.E)} := by
  rw [levelModulus, Ideal.mem_span_singleton]
  refine ⟨41 * rootInteger, ?_⟩
  apply RingOfIntegers.ext
  have e1 : ((integerGenerator - 1 : 𝓞 host.E) : host.E) = (multiplier : host.E) * host.sqd := by
    show generator - 1 = _
    unfold generator
    ring
  rw [e1]
  push_cast
  show (multiplier : host.E) * host.sqd = (88134043685124465719294615021161552534088713491754576 : host.E) * (41 * host.sqd)
  unfold multiplier
  push_cast
  ring

instance factLevelPrime : Fact levelPrime.Prime := ⟨levelPrime_prime⟩

theorem typesSplit : ∀ k, (placesOver host.E (host.g k)).ncard = Module.finrank ℚ host.E := by
  refine g_split_of host (fun k => ?_)
  fin_cases k
  · show 343057911569 % (2 ^ (2 + 2) * 547 * 727) = 1; norm_num
  · show 343070636977 % (2 ^ (2 + 2) * 547 * 727) = 1; norm_num

theorem typeConditions : THyp host where
  cond := fun k => by
    fin_cases k
    · show 343057911569 % (2 ^ (2 + 2) * 547 * 727) = 1; norm_num
    · show 343070636977 % (2 ^ (2 + 2) * 547 * 727) = 1; norm_num
  bb := disjointExceptional
  jac := fun k => by
    fin_cases k
    · exact ⟨by show (343057911569 : ℕ) ≠ 2; norm_num, by show jacobiSym (-727) 343057911569 = 1; norm_num⟩
    · exact ⟨by show (343070636977 : ℕ) ≠ 2; norm_num, by show jacobiSym (-727) 343070636977 = 1; norm_num⟩

theorem levelPrime_distinct : levelPrime ≠ host.pde ∧ levelPrime ≠ host.pde' :=
  ⟨by show levelPrime ≠ 991; norm_num [levelPrime], by show levelPrime ≠ 1231; norm_num [levelPrime]⟩

theorem distinguishedSplit : (placesOver host.E host.pde).ncard = Module.finrank ℚ host.E ∧
    (placesOver host.E host.pde').ncard = Module.finrank ℚ host.E :=
  ⟨placesOver_ncard_of_split_inj host exampleRanks_h_injective 991 (hp := ⟨by norm_num⟩)
      (by show ¬ (991 : ℕ) ∣ 2 ^ (2 + 2) * 547 * 727; norm_num) (by show jacobiSym (-547) 991 = 1; norm_num)
      (by show jacobiSym (-727) 991 = 1; norm_num)
      (Or.inr (by show ((991 : ℕ) : ZMod 16) = -1; decide +kernel)),
    placesOver_ncard_of_split_inj host exampleRanks_h_injective 1231 (hp := ⟨by norm_num⟩)
      (by show ¬ (1231 : ℕ) ∣ 2 ^ (2 + 2) * 547 * 727; norm_num) (by show jacobiSym (-547) 1231 = 1; norm_num)
      (by show jacobiSym (-727) 1231 = 1; norm_num)
      (Or.inr (by show ((1231 : ℕ) : ZMod 16) = -1; decide +kernel))⟩

def makeRecord (G : LocalGap host) : D7 :=
  exampleRecordOfData host levelPrime typesSplit typeConditions levelPrime_mem levelPrime_distinct integerGenerator generatorMem generatorNorm generatorCongruence G

def record : D7 := makeRecord (GlobalTheory.MemLoc.localGap_nonempty host).some

theorem field_parameters : record.R.e = 2 ∧ record.host.q = 547 ∧ record.host.d = 727 :=
  ⟨rfl, rfl, rfl⟩

theorem field_degrees : Module.finrank ℚ record.host.E = 16 ∧
    Module.finrank ℚ record.host.F = 8 := by
  exact ⟨Formalism.FrobeniusOrbits.finrank_cmField (e := 2) hyp.toHypOdd layer,
    Formalism.FrobeniusOrbits.finrank_hostF (e := 2) hyp.toHypOdd layer⟩

theorem frame_exists : Nonempty (CharFrame record) := by
  obtain ⟨⟨P⟩, _⟩ := GlobalTheory.nonempty_piBar_and_omg record.R record.host
  have h (j : Fin record.R.m) : UpperTheory.MUQ
      (fun τ => record.host.emb (record.eTau τ : record.R.GG)) (fun τ => (record.hc τ).e j) := by
    obtain ⟨n, hn⟩ := hc0.e_int j
    exact UpperTheory.Mu.muq_of_const record.host record.eTau n (fun _ => hn)
  obtain ⟨C, _⟩ := GlobalTheory.BigOn.Inhab.exists_charFrame record P
    distinguishedSplit.1 distinguishedSplit.2 h
  exact ⟨C⟩

end SecondHost
end
