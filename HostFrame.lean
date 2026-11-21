import Development.AbelianVarietyTheory.Proofs.InertiaEigenvalues.LocalGapNonemptyAndRecordNonempty
import Development.AbelianVarietyTheory.Proofs.DivisorProofs.IsTwistPlacesAboveK
import Development.AutomorphicArithmetic.WeightCongruences.IsRootOfIntegral

noncomputable section

open Formalism.BasicNotions Formalism.BasicNotions.Examples₂ Formalism.L2.Host

namespace HostRace

def record : D7 := exampleRecord (GlobalTheory.MemLoc.localGap_nonempty H0).some

theorem charFrame_exists : Nonempty (CharFrame record) := by
  obtain ⟨⟨P⟩, _⟩ := GlobalTheory.nonempty_piBar_and_omg record.R record.host
  have hs := exampleHost_split_prime_split
  have hm (j : Fin record.R.m) : UpperTheory.MUQ
      (fun τ => record.host.emb (record.eTau τ : record.R.GG)) (fun τ => (record.hc τ).e j) := by
    obtain ⟨n, hn⟩ := hc0.e_int j
    exact UpperTheory.Mu.muq_of_const record.host record.eTau n (fun _ => hn)
  obtain ⟨C, _⟩ := GlobalTheory.BigOn.Inhab.exists_charFrame record P hs.1 hs.2 hm
  exact ⟨C⟩

end HostRace

#print axioms HostRace.charFrame_exists
