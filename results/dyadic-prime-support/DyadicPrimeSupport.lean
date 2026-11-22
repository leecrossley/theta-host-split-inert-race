import Mathlib

namespace DyadicPrimeSupport

abbrev Interior := Set.Ioo (0 : ℝ) 1

def line (u v x : ℝ) : ℝ := (1-x)*u+x*v

theorem two_zeros {u v x z : ℝ} (hxz : x ≠ z)
    (hx : line u v x = 0) (hz : line u v z = 0) : u = 0 ∧ v = 0 := by
  unfold line at hx hz
  have hd : (x-z)*(v-u) = 0 := by linear_combination hx-hz
  have huv : v = u := sub_eq_zero.mp ((mul_eq_zero.mp hd).resolve_left (sub_ne_zero.mpr hxz))
  rw [huv] at hx ⊢
  constructor <;> nlinarith only [hx]

theorem root_propagates {u v : ℝ} (t : Interior) (hv : v ≠ 0)
    (hx : line u v t = 0) : u ≠ 0 := by
  intro hu
  have htv : (t : ℝ)*v = 0 := by simpa [line, hu] using hx
  exact (mul_ne_zero (ne_of_gt t.2.1) hv) htv

theorem prescribe_root (v : ℝ) (t : Interior) :
    line (-((t : ℝ)*v)/(1-(t : ℝ))) v t = 0 := by
  have h : (1 : ℝ) - t ≠ 0 := sub_ne_zero.mpr (ne_of_gt t.2.2)
  unfold line
  field_simp [h]
  ring

def Silent : List (Finset Interior) → ℝ → Prop
  | [], _ => True
  | cell :: cells, v => ∃ u : ℝ, (∀ t ∈ cell, line u v t = 0) ∧ Silent cells u

def Passes : List (Finset Interior) → Prop
  | [] => True
  | cell :: cells => cell = ∅ ∨ (cell.card = 1 ∧ Passes cells)

theorem silent_zero (cells : List (Finset Interior)) : Silent cells 0 := by
  induction cells with
  | nil => trivial
  | cons cell cells ih => exact ⟨0, fun t _ => by simp [line], ih⟩

theorem invisible_iff (cells : List (Finset Interior)) (v : ℝ) (hv : v ≠ 0) :
    Silent cells v ↔ Passes cells := by
  classical
  induction cells generalizing v with
  | nil => exact Iff.rfl
  | cons cell cells ih =>
    change (∃ u : ℝ, (∀ t ∈ cell, line u v t = 0) ∧ Silent cells u) ↔
        cell = ∅ ∨ cell.card = 1 ∧ Passes cells
    by_cases he : cell = ∅
    · subst cell
      exact ⟨fun _ => Or.inl rfl, fun _ => ⟨0, by simp, silent_zero cells⟩⟩
    by_cases ho : cell.card = 1
    · obtain ⟨t, ht⟩ := Finset.card_eq_one.mp ho
      constructor
      · rintro ⟨u, hline, hcells⟩
        have hz := hline t (by simp [ht])
        exact Or.inr ⟨ho, (ih u (root_propagates t hv hz)).mp hcells⟩
      · intro good
        have htail := (good.resolve_left he).2
        let u := -((t : ℝ)*v)/(1-(t : ℝ))
        have hx : line u v t = 0 := prescribe_root v t
        refine ⟨u, ?_, (ih u (root_propagates t hv hx)).mpr htail⟩
        intro z hz
        have hzt : z = t := by simpa [ht] using hz
        simpa [hzt] using hx
    · have hn : 0 < cell.card := Finset.card_pos.mpr (Finset.nonempty_iff_ne_empty.mpr he)
      have hcard : 1 < cell.card := by omega
      obtain ⟨t, ht, z, hz, htz⟩ := Finset.one_lt_card.mp hcard
      constructor
      · rintro ⟨u, hline, _⟩
        have hne : (t : ℝ) ≠ (z : ℝ) := fun h => htz (Subtype.ext h)
        exact (hv (two_zeros hne (hline t ht) (hline z hz)).2).elim
      · intro good
        exact (ho (good.resolve_left he).1).elim

theorem opposite_linear_terms {a b c d x y z t : ℝ} (ha : a ≠ 0)
    (hxy : x ≠ y) (hzt : z ≠ t) (hs : 0 < x+y+z+t)
    (hx : a*x^2+b*x+c=0) (hy : a*y^2+b*y+c=0)
    (hz : a*z^2-b*z+d=0) (ht : a*t^2-b*t+d=0) : False := by
  have h₁ : (x-y)*(a*(x+y)+b) = 0 := by linear_combination hx-hy
  have h₂ : (z-t)*(a*(z+t)-b) = 0 := by linear_combination hz-ht
  have h₃ := (mul_eq_zero.mp h₁).resolve_left (sub_ne_zero.mpr hxy)
  have h₄ := (mul_eq_zero.mp h₂).resolve_left (sub_ne_zero.mpr hzt)
  have h₅ : a*(x+y+z+t) = 0 := by linear_combination h₃+h₄
  exact (mul_ne_zero ha (ne_of_gt hs)) h₅

theorem dyadic_sum (S : ℕ → ℝ) (n : ℕ) :
    ∑ j ∈ Finset.range (n+1), (S (j+1)-S j) = S (n+1)-S 0 := by
  induction n with
  | zero => simp
  | succ n ih => rw [Finset.sum_range_succ, ih]; ring

end DyadicPrimeSupport

#print axioms DyadicPrimeSupport.two_zeros
#print axioms DyadicPrimeSupport.root_propagates
#print axioms DyadicPrimeSupport.prescribe_root
#print axioms DyadicPrimeSupport.silent_zero
#print axioms DyadicPrimeSupport.invisible_iff
#print axioms DyadicPrimeSupport.opposite_linear_terms
#print axioms DyadicPrimeSupport.dyadic_sum
