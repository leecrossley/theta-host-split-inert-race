import Mathlib

namespace DoublingPatterns

open scoped BigOperators

theorem sign_preserved {s x e t m : ℝ} (hs : s=1 ∨ s= -1)
    (hx : m < s*x) (he : |e| < t) : m-t < s*(x+e) := by
  rw [abs_lt] at he
  rcases hs with h | h <;> subst s <;> linarith

theorem separated_signs {x y e : ℝ} (he : 0 ≤ e) (hx : e < |x|)
    (hy : e < |y|) (h : |x-y| < 2*e) : (0 < x ↔ 0 < y) := by
  rw [abs_lt] at h
  rw [lt_abs] at hx hy
  constructor <;> intro hp <;> rcases hx with hx | hx <;> rcases hy with hy | hy <;> linarith

theorem move_error {x y x' y' e a : ℝ} (ha : 0 < a)
    (hx : |x-x'| ≤ e) (hy : |y-y'| ≤ e) :
    |(y-x)/a-(y'-x')/a| ≤ 2*e/a := by
  rw [← sub_div, abs_div, abs_of_pos ha]
  apply div_le_div_of_nonneg_right _ ha.le
  rw [abs_le] at hx hy ⊢
  constructor <;> linarith

theorem telescoping (S : ℕ → ℝ) (n : ℕ) :
    ∑ j ∈ Finset.range n, (S (j+1)-S j) = S n-S 0 := by
  induction n with
  | zero => simp
  | succ n ih => rw [Finset.sum_range_succ, ih]; ring

theorem excluded_run {m e M : ℝ} (he : 2*e < m-2*M)
    {s t : ℝ} (hs : 1-M-e ≤ s) (ht : t ≤ 1+M+e) : t-s < m := by
  linarith

theorem word_count {W C : Type*} [Fintype W] [Fintype C]
    (cell : W → C) (h : Function.Injective cell) : Fintype.card W ≤ Fintype.card C :=
  Fintype.card_le_of_injective cell h

theorem reflect_difference (f : ℝ → ℝ) (he : Function.Even f) (x a : ℝ) :
    (f (-x-a+a)-f (-x-a))/a = -(f (x+a)-f x)/a := by
  rw [show -x-a+a = -x by ring, show -x-a = -(x+a) by ring, he x, he (x+a)]
  ring

theorem cubic_majorant (z : ℝ) :
    z ≤ (11 + z^2*(75+z^2*(16*z^2-48)))/54 := by
  have h : 0 ≤ (8/27:ℝ)*(z-1/2)^2*(z-1)^2*((z+3/2)^2+1/2) := by positivity
  nlinarith only [h]


theorem strict_peak_of_opposite_branches {I : Set ℝ} {f g : ℝ → ℝ} {c : ℝ}
    (hc : c ∈ I) (hf : StrictMonoOn f I) (hg : StrictAntiOn g I)
    (heq : f c = g c) {x : ℝ} (hx : x ∈ I) (hne : x ≠ c) :
    min (f x) (g x) < f c := by
  rcases lt_or_gt_of_ne hne with h | h
  · exact lt_of_le_of_lt (min_le_left _ _) (hf hx hc h)
  · rw [heq]
    exact lt_of_le_of_lt (min_le_right _ _) (hg hc hx h)

theorem catalogue_iff_last {ι α : Type*} {last : ι} {I : Set α} (occurs : ι → α → Prop)
    (early : ∀ i, i ≠ last → ∃ x ∈ I, occurs i x) :
    (∀ i, ∃ x ∈ I, occurs i x) ↔ ∃ x ∈ I, occurs last x := by
  constructor
  · exact fun h => h last
  · intro h i
    by_cases hi : i = last
    · simpa [hi] using h
    · exact early i hi

theorem cell_intersection {u v s t x : ℝ} (h : u ≠ s)
    (eq : u*x+v = s*x+t) : x = (t-v)/(u-s) := by
  apply (eq_div_iff (sub_ne_zero.mpr h)).2
  linarith only [eq, mul_comm x (u-s), sub_mul u s x]

theorem phase_count_needed {m : ℕ} (hm : 0 < m) (gap : Fin m → ℝ) {w : ℝ}
    (total : ∑ i, gap i = 1) (small : ∀ i, gap i < w) : 1 < (m : ℝ) * w := by
  let : Nonempty (Fin m) := ⟨⟨0, hm⟩⟩
  calc
    1 = ∑ i, gap i := total.symm
    _ < ∑ _i : Fin m, w := Finset.sum_lt_sum_of_nonempty Finset.univ_nonempty (fun i _ => small i)
    _ = (m : ℝ) * w := by simp


theorem successor_above_crossing {G : Set ℝ} {a s : ℝ} {success : ℝ → Prop}
    (hs : s ∈ G) (ha : a < s) (next : ∀ t ∈ G, a < t → s ≤ t)
    (works : a < s → success s) (before : ∀ t ∈ G, t < s → success t → a < t) :
    IsLeast {t | t ∈ G ∧ success t} s := by
  refine ⟨⟨hs, works ha⟩, ?_⟩
  intro t ht
  by_contra! h
  exact (not_lt_of_ge (next t ht.1 (before t ht.1 h ht.2))) h

theorem three_point_slope_jump {v left right step : ℝ} (hh : step ≠ 0) :
    ((v-left*step)-2*v+(v+right*step))/step = right-left := by
  apply (div_eq_iff hh).2
  ring

theorem recover_jump_off_grid {n t v left right : ℝ} (ht : t < n) :
    ((v+left*(t-n))-(v+right*(t-n)))/(n-t) = right-left := by
  apply (div_eq_iff (sub_ne_zero.mpr (ne_of_gt ht))).2
  ring

theorem last_surviving_combination {x y z : ℝ}
    (h1 : 31*x = 0) (h2 : 47*x+12*y+8*z = 0) :
    ∃ t : ℝ, x = 0 ∧ y = -2*t ∧ z = 3*t := by
  exact ⟨z/3, by linarith, by linarith, by ring⟩

theorem last_row_removes_combination {x y z : ℝ}
    (h1 : 31*x = 0) (h2 : 47*x+12*y+8*z = 0) (h3 : 337*x+302*y+298*z = 0) :
    x = 0 ∧ y = 0 ∧ z = 0 := by
  exact ⟨by linarith, by linarith, by linarith⟩


theorem two_sample_minor (x y u v : ℝ) :
    (x-u)*(y-v)-(x-v)*(y-u) = (x-y)*(u-v) := by ring

theorem two_backgrounds {x y u v a b : ℝ} (hxy : x ≠ y) (huv : u ≠ v)
    (h1 : a*(x-u)+b*(y-u) = 0) (h2 : a*(x-v)+b*(y-v) = 0) :
    a = 0 ∧ b = 0 := by
  have hz : (u-v)*(a+b) = 0 := by linear_combination h2-h1
  have hab : a+b = 0 := (mul_eq_zero.mp hz).resolve_left (sub_ne_zero.mpr huv)
  have hx : a*(x-y) = 0 := by linear_combination h1-(y-u)*hab
  have ha : a = 0 := (mul_eq_zero.mp hx).resolve_right (sub_ne_zero.mpr hxy)
  exact ⟨ha, by linarith⟩

theorem affine_extrapolation {x y : ℝ} (h : y ≠ x) (s u v : ℝ) :
    ((y-s)*(u*x+v)+(s-x)*(u*y+v))/(y-x) = u*s+v := by
  apply (div_eq_iff (sub_ne_zero.mpr h)).2
  ring

theorem first_strict_record {score : ℕ → ℝ} {q : ℝ} {n : ℕ}
    (hit : q < score n) (earlier : ∀ j < n, score j ≤ q) :
    IsLeast {j | q < score j} n := by
  refine ⟨hit, ?_⟩
  intro j hj
  by_contra! h
  exact (not_lt_of_ge (earlier j h)) hj

theorem last_prime_pivot {ι κ K : Type*} [LinearOrder ι] [Fintype ι] [Field K]
    (row : ι → κ → K) (pivot : ι → κ) (weight : ι → K)
    (diagonal : ∀ i, row i (pivot i) ≠ 0)
    (older : ∀ i j, j < i → row j (pivot i) = 0)
    (relation : ∀ p, ∑ i, weight i * row i p = 0) :
    ∀ i, weight i = 0 := by
  classical
  by_contra! h
  let used := Finset.univ.filter fun i => weight i ≠ 0
  have hn : used.Nonempty := by
    obtain ⟨i, hi⟩ := h
    exact ⟨i, Finset.mem_filter.mpr ⟨Finset.mem_univ i, hi⟩⟩
  let m := used.max' hn
  have hm : m ∈ used := used.max'_mem hn
  have hwm : weight m ≠ 0 := (Finset.mem_filter.mp hm).2
  have hsingle : ∑ i, weight i * row i (pivot m) = weight m * row m (pivot m) := by
    apply Finset.sum_eq_single m
    · intro j _ hj
      rcases lt_or_gt_of_ne hj with before | after
      · rw [older m j before, mul_zero]
      · have hz : weight j = 0 := by
          by_contra hnj
          have mem : j ∈ used := Finset.mem_filter.mpr ⟨Finset.mem_univ j, hnj⟩
          exact (not_lt_of_ge (Finset.le_max' used j mem)) after
        rw [hz, zero_mul]
    · simp
  exact (mul_ne_zero hwm (diagonal m)) (hsingle.symm.trans (relation (pivot m)))

end DoublingPatterns

#print axioms DoublingPatterns.sign_preserved
#print axioms DoublingPatterns.separated_signs
#print axioms DoublingPatterns.move_error
#print axioms DoublingPatterns.telescoping
#print axioms DoublingPatterns.excluded_run
#print axioms DoublingPatterns.word_count
#print axioms DoublingPatterns.reflect_difference
#print axioms DoublingPatterns.cubic_majorant
#print axioms DoublingPatterns.strict_peak_of_opposite_branches
#print axioms DoublingPatterns.catalogue_iff_last
#print axioms DoublingPatterns.cell_intersection
#print axioms DoublingPatterns.phase_count_needed
#print axioms DoublingPatterns.successor_above_crossing
#print axioms DoublingPatterns.three_point_slope_jump
#print axioms DoublingPatterns.recover_jump_off_grid
#print axioms DoublingPatterns.last_surviving_combination
#print axioms DoublingPatterns.last_row_removes_combination
#print axioms DoublingPatterns.two_sample_minor
#print axioms DoublingPatterns.two_backgrounds
#print axioms DoublingPatterns.affine_extrapolation
#print axioms DoublingPatterns.first_strict_record
#print axioms DoublingPatterns.last_prime_pivot
