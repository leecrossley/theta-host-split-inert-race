import ThetaPrimeTrace
import Development.AbelianVarietyTheory.Proofs.FoundationProofs.RamificationIdxEqOneOfIsSquare

open Formalism.BasicNotions GlobalTheory.ChainOnto GlobalTheory.ChainOnto.ParClosed
open NumberField IsDedekindDomain
open Formalism.L2.Host (aboveF)
open Formalism.MainNotions.SchurIrreducible (IsUnram)
attribute [local instance] Formalism.BasicNotions.Host.cmStar

noncomputable section
namespace HostRace

variable {D : D7}

def jump (c : Com D) (B : GlobalTheory.LFun.Fam D.host.E)
    (w : HeightOneSpectrum (𝓞 D.host.E)) : ℂ :=
  (ThetaTransfer.firstTrace (chainPar c B 1 w) -
    ((vf (c.chS 1 * (c.chT 1)⁻¹) w : ℂˣ) : ℂ) *
      ThetaTransfer.firstTrace (chainPar c B 0 w)) /
    ThetaTransfer.firstTrace (GlobalTheory.LFun.seg w 1 D.R.t)

theorem jump_eq (c : Com D) (B : GlobalTheory.LFun.Fam D.host.E)
    (w : HeightOneSpectrum (𝓞 D.host.E)) :
    jump c B w = ((vf (c.chS 1) w : ℂˣ) : ℂ) := by
  classical
  simpa [jump] using
    ThetaTransfer.theta_prime_trace_extract c B 0 {w} (fun _ => (1 : ℂ))

/-- At a good base prime the first-step trace product is its quadratic splitting sign. -/
theorem fiber_trace (c : Com D) (hc : c.Is) (B : GlobalTheory.LFun.Fam D.host.E)
    (v : HeightOneSpectrum (𝓞 D.host.F)) (hv : IsUnram (E := D.host.E) v) (hs : v ∉ c.S) :
    ∏ w ∈ aboveF (L := D.host.E) D.host.F v, jump c B w =
      (GlobalTheory.LFun.etaU D.host.E v : ℂˣ) := by
  classical
  have hm : 1 ∈ Finset.Icc 1 D.R.m := Finset.mem_Icc.mpr ⟨le_rfl, by have := D.R.two_le_m; omega⟩
  have hdeg := Algebra.IsQuadraticExtension.finrank_eq_two D.host.F D.host.E
  have hu : ∀ w : HeightOneSpectrum (𝓞 D.host.E), w.under (𝓞 D.host.F) = v →
      GlobalTheory.LFun.UnrAt (c.μ 1) w := by
    intro w hw
    apply hc.unr ⟨0, by have := D.R.two_le_m; omega⟩
    simpa [GlobalTheory.LFun.above, hw] using hs
  have he := GlobalTheory.LFun.ResEta.prod_val hdeg
    (show GlobalTheory.LFun.ResEta D.host.F D.host.E c.η 1 from hc.chars.2.2.1) hv
    (fun w _ => GlobalTheory.LFun.unrAt_of_unramEverywhere hc.chars.2.1 w)
  have hmu := GlobalTheory.LFun.ResEta.prod_val hdeg
    (show GlobalTheory.LFun.ResEta D.host.F D.host.E (c.μ 1) 0 from (hc.chars.2.2.2.2 1 hm).2) hv hu
  have hodd : c.chS 1 = c.μ 1 * c.η⁻¹ := by
    simp [Com.chS, Formalism.IndefR.Lift.chrStep]
  simp_rw [jump_eq, hodd]
  have h : ∏ w ∈ aboveF (L := D.host.E) D.host.F v, vf (c.μ 1 * c.η⁻¹) w =
      GlobalTheory.LFun.etaU D.host.E v := by
    simp only [vf, GlobalTheory.LFun.val_mul, GlobalTheory.LFun.val_inv,
      Finset.prod_mul_distrib, Finset.prod_inv_distrib, hmu, he, pow_zero, pow_one, one_mul]
    unfold GlobalTheory.LFun.etaU
    split <;> simp
  exact_mod_cast h

end HostRace

#print axioms HostRace.jump_eq
#print axioms HostRace.fiber_trace
