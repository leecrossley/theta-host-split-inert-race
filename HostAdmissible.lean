import FiberTrace
import HostFrame
import Development.AbelianVarietyTheory.Proofs.IwasawaDecomposition.IsUnitaryOfOrderDvd
import Development.AbelianVarietyTheory.Proofs.IwasawaDecomposition.ExistsIsInfFamFamily
import Development.AbelianVarietyTheory.Proofs.IwasawaDecomposition.ExistsCharOfTrace
import Development.AbelianVarietyTheory.Proofs.IwasawaDecomposition.ExistsTupChars
import Development.AbelianVarietyTheory.Proofs.IwasawaDecomposition.HasInfExpMulZero
import Development.AbelianVarietyTheory.Proofs.DivisorProofs.ExistsWeilRepresentationVacuumCosetVectorsGrid
import Development.AbelianVarietyTheory.Proofs.DivisorProofs.TypePlacesNotPairwiseConjugate
import Development.AbelianVarietyTheory.Proofs.DivisorProofs.CentralShiftExists
import Development.AbelianVarietyTheory.Proofs.DivisorProofs.UnderEqOrSmulUnderEq
import Development.AbelianVarietyTheory.Proofs.SatakeParameters.IsIdemOfAlgebraMap

noncomputable section

open NumberField IsDedekindDomain Formalism.BasicNotions Formalism.BasicNotions.TH
open Formalism.L2.Theta Formalism.L2.DIV Formalism.IndefR Formalism.IndefR.Links
open Formalism.IndefR.Mult Formalism.IndefR.Datum GlobalTheory.ChainOnto
attribute [local instance] Formalism.BasicNotions.Host.cmStar

namespace HostRace

/-- Assemble the untwisted level choices, keeping their original Hecke characters. -/
private theorem choices {D : D7} (C : Formalism.L2.Host.CharFrame D) :
    ∃ c : Com D, c.Is := by
  obtain ⟨L, _, hL⟩ := GlobalTheory.LevelChoices.exists_levCh_unit D C
  obtain ⟨θ, hθ⟩ := GlobalTheory.Member.exists_isInfFam_family D
  obtain ⟨ψ, hψ⟩ := GlobalTheory.StdAdelicChar.exists_isStdChar D
  let c : Com D := ⟨θ, ψ, L.mu0, L.eta, L.Pg, L.triv, L.chi⟩
  have hcover (k : Fin 2) (w : Pl D.host.E)
      (hk : ((D.host.g k : ℕ) : 𝓞 D.host.E) ∈ w.asIdeal) : w.under (𝓞 D.host.F) ∈ c.S := by
    classical
    apply Finset.mem_image.mpr
    rcases GlobalTheory.ConjOverQuad.div_or_conj_div D L.C.P k w hk with h | ⟨w', h, eq⟩
    · exact ⟨w, (hL.Pg w).2 ⟨k, h⟩, rfl⟩
    · refine ⟨w', (hL.Pg w').2 ⟨k, h⟩, ?_⟩
      ext x
      show algebraMap (𝓞 D.host.F) (𝓞 D.host.E) x ∈ w'.asIdeal ↔
        algebraMap (𝓞 D.host.F) (𝓞 D.host.E) x ∈ w.asIdeal
      rw [eq, Formalism.L2.Col.mem_conjIdeal]
      have fixed : Formalism.L2.Col.cO D.host.E (algebraMap (𝓞 D.host.F) (𝓞 D.host.E) x) =
          algebraMap _ _ x := RingOfIntegers.ext ((D.sp.cm.fixed _).2 ⟨x, rfl⟩)
      rw [fixed]
  have hres := GlobalTheory.Mem.EtaTilde.isEta_res D L.eta L.C.eta hL.eta L.C.isEta
  refine ⟨c, hθ, hψ, ⟨hL.eta.1, GlobalTheory.Mem.TupChar.unramEverywhere hL.eta L.C.isEta.unr,
    ⟨⟨_, hres⟩, fun η h a => ?_⟩, fun τ => ?_, fun j hj => ?_⟩,
    fun j τ => ?_, fun j w hw => ?_, fun j => ⟨_, hL.mu0 j, hL.chi j⟩, ?_⟩
  · rw [(GlobalTheory.IdeleNormChar.existsUnique_isEta D).unique h hres, pow_one, MonoidHom.comp_apply]
  · exact GlobalTheory.HeckeInfExp.hasInfExp_of_chrEq D _ _ _ hL.eta
      (GlobalTheory.InnerSteps.hasInfQ_tauE (a := fun _ => 1 / 2)
        ⟨fun _ => 1, fun _ => by norm_num, L.C.isEta.inf⟩) τ
  · obtain ⟨h1, h2⟩ := Finset.mem_Icc.mp hj
    obtain ⟨i, rfl⟩ : ∃ i, j = i + 1 := ⟨j - 1, by omega⟩
    let t : Fin D.R.m := ⟨i, by omega⟩
    refine ⟨(hL.mu0 t).1, ⟨_, hres⟩, fun η _ a => ?_⟩
    rw [pow_zero]
    exact GlobalTheory.GridWeil.triv_base (hL.mu0 t) (L.C.isMu0 t).triv a
  · have h := GlobalTheory.HeckeInfExp.hasInfExp_of_chrEq D _ _ _ (hL.mu0 j)
      (GlobalTheory.InnerSteps.hasInfQ_tauE (L.C.isMu0 j).inf) τ
    simpa only [Formalism.L2.Host.eExp, Equiv.symm_apply_apply] using h
  · exact (GlobalTheory.LFun.unrAt_iff_of_chrEq (hL.mu0 j) w).2
      ((L.C.isMu0 j).unr w fun k hk => hw (hcover k w hk))
  · exact GlobalTheory.CentralShift.exists_idemSum L.Pg
      (GlobalTheory.TypePlaces.not_isConjPl L.C.P L.Pg hL.Pg)

theorem admissible_exists : ∃ c : Com record, c.Is :=
  choices charFrame_exists.some

private theorem host_unram (D : D7) (v : HeightOneSpectrum (𝓞 D.host.F)) :
    Formalism.MainNotions.SchurIrreducible.IsUnram (E := D.host.E) v := by
  haveI := D.host.isCMField_hostE
  have hu := GlobalTheory.PlaceComp.unramified (E := D.host.E) D.sp.cm
  apply GlobalTheory.LFun.isUnram_iff.mpr
  intro w _
  haveI := w.isPrime
  haveI : Algebra.IsUnramifiedAt (𝓞 D.host.F) w.asIdeal :=
    Algebra.unramified_iff_forall.mp hu ⟨w.asIdeal, w.isPrime⟩
  exact w.asIdeal.ramificationIdx_eq_one (𝓞 D.host.F)

/-- One fixed admissible datum recovers every good fibre, for every input family. -/
theorem concrete_fiber_trace :
    ∃ c : Com record, c.Is ∧ ∀ (B : GlobalTheory.LFun.Fam record.host.E)
      (v : HeightOneSpectrum (𝓞 record.host.F)), v ∉ c.S →
      ∏ w ∈ Formalism.L2.Host.aboveF (L := record.host.E) record.host.F v, jump c B w =
        (GlobalTheory.LFun.etaU record.host.E v : ℂˣ) := by
  obtain ⟨c, hc⟩ := admissible_exists
  exact ⟨c, hc, fun B v hs => fiber_trace c hc B v (host_unram record v) hs⟩

end HostRace

#print axioms HostRace.admissible_exists
#print axioms HostRace.concrete_fiber_trace
