import Mathlib

namespace CrossFieldPrediction

theorem private_frequency (target other multiplier : ℂ) (h : other = 0) :
    target - multiplier * other = target := by rw [h, mul_zero, sub_zero]

theorem dual_energy (energy inverseEntry radius : ℝ)
    (positive : 0 < inverseEntry) (normalised : 1 ≤ energy * inverseEntry)
    (certificate : radius^2 * inverseEntry < 1) : radius^2 < energy := by
  exact lt_of_mul_lt_mul_right (lt_of_lt_of_le certificate normalised) positive.le

theorem reflection_step (variance covariance k : ℝ) (variance_pos : 0 < variance)
    (stationary : covariance = k*variance) (other : ℝ) :
    variance*(1-k^2) ≤ variance - 2*other*covariance + other^2*variance := by
  rw [stationary]
  nlinarith only [mul_nonneg (sq_nonneg (other-k)) variance_pos.le]

theorem perturbation (trueEnergy centreEnergy cost normSq margin : ℝ)
    (err : |trueEnergy-centreEnergy| ≤ cost*normSq)
    (norm_nonneg : 0 ≤ normSq) (round_up : cost ≤ margin) :
    centreEnergy-margin*normSq ≤ trueEnergy := by
  have bound := (abs_le.mp err).1
  have scale := mul_le_mul_of_nonneg_right round_up norm_nonneg
  linarith

theorem two_energies (wanted hidden intercept : ℝ)
    (rest : 0 ≤ hidden) : wanted ≤ wanted+hidden+intercept^2 := by
  linarith [sq_nonneg intercept]

theorem strict_rms (observed bound : ℝ) (ho : 0 ≤ observed)
    (h : bound^2 < observed^2) : bound < observed := by
  by_contra! hn
  have hh := mul_self_le_mul_self ho hn
  nlinarith only [hh,h]

end CrossFieldPrediction

#print axioms CrossFieldPrediction.private_frequency
#print axioms CrossFieldPrediction.dual_energy
#print axioms CrossFieldPrediction.reflection_step
#print axioms CrossFieldPrediction.perturbation
#print axioms CrossFieldPrediction.two_energies
#print axioms CrossFieldPrediction.strict_rms
