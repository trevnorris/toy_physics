import S10Audit.CAS.MinorLoci
import S10Audit.CAS.MinorBindings

namespace S10Audit.CAS.WLLoci
open S10Pilot
noncomputable section

def staticRank (k : Vec 3) : Prop := ((k 0 = 0 ∧ k 1 = 0 ∧ k 2 = 0))

theorem staticRank_geometry (k : Vec 3) : staticRank k ↔ (k = 0) := by
  classical
  locus_eval staticRank

theorem staticRank_minors (rho mu sigma z : ℝ) (k : Vec 3)
    (hm : mu ≠ 0) (hs : sigma ≠ 0) (_hs1 : sigma ≠ 1) :
    staticRank k ↔ (∀ i, WLMinors.staticRank rho mu sigma z k i = 0) := by
  rw [staticRank_geometry, WLMinors.staticRank_zero_iff rho mu sigma z k hs,
    CAS.staticRank_locus mu k hm]

theorem staticRank_all_minors (rho mu sigma : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : sigma ≠ 0) (hs1 : sigma ≠ 1) :
    staticRank k ↔ AllOrderedMinorsZero
      (rootMatrix rho mu sigma k 0) 2 := by
  rw [staticRank_minors rho mu sigma 0 k hm hs hs1,
    WLMinors.staticRank_complete rho mu sigma 0 k hr hs]

def staticTransverse (k : Vec 3) : Prop := ((k 0 = 0 ∧ k 1 = 0 ∧ k 2 = 0))

theorem staticTransverse_geometry (k : Vec 3) : staticTransverse k ↔ (k = 0) := by
  classical
  locus_eval staticTransverse

theorem staticTransverse_minors (rho mu sigma z : ℝ) (k : Vec 3)
    (hm : mu ≠ 0) (hs : sigma ≠ 0) (_hs1 : sigma ≠ 1) :
    staticTransverse k ↔ (∀ i, WLMinors.staticTransverse rho mu sigma z k i = 0) := by
  rw [staticTransverse_geometry, WLMinors.staticTransverse_zero_iff rho mu sigma z k hs,
    CAS.staticTransverse_locus mu k hm]

theorem staticTransverse_all_minors (rho mu sigma : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : sigma ≠ 0) (hs1 : sigma ≠ 1) :
    staticTransverse k ↔ AllOrderedMinorsZero
      (rootStack rho mu sigma k 0) 3 := by
  rw [staticTransverse_minors rho mu sigma 0 k hm hs hs1,
    WLMinors.staticTransverse_complete rho mu sigma 0 k hr hs]

def ordinaryRank (k : Vec 3) : Prop := ((k 1 = 0 ∧ k 2 = 0))

theorem ordinaryRank_geometry (k : Vec 3) : ordinaryRank k ↔ (k 1 = 0 ∧ k 2 = 0) := by
  classical
  locus_eval ordinaryRank

theorem ordinaryRank_minors (rho mu sigma z : ℝ) (k : Vec 3)
    (hm : mu ≠ 0) (hs : sigma ≠ 0) (_hs1 : sigma ≠ 1) :
    ordinaryRank k ↔ (∀ i, WLMinors.ordinaryRank rho mu sigma z k i = 0) := by
  rw [ordinaryRank_geometry, WLMinors.ordinaryRank_zero_iff rho mu sigma z k hs,
    CAS.ordinaryRank_locus mu sigma k hm _hs1]

theorem ordinaryRank_all_minors (rho mu sigma : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : sigma ≠ 0) (hs1 : sigma ≠ 1) :
    ordinaryRank k ↔ AllOrderedMinorsZero
      (rootMatrix rho mu sigma k 1) 2 := by
  rw [ordinaryRank_minors rho mu sigma 0 k hm hs hs1,
    WLMinors.ordinaryRank_complete rho mu sigma 0 k hr hs]

def ordinaryTransverse (k : Vec 3) : Prop := ((k 1 = 0 ∧ k 2 = 0))

theorem ordinaryTransverse_geometry (k : Vec 3) : ordinaryTransverse k ↔ (k 1 = 0 ∧ k 2 = 0) := by
  classical
  locus_eval ordinaryTransverse

theorem ordinaryTransverse_minors (rho mu sigma z : ℝ) (k : Vec 3)
    (hm : mu ≠ 0) (hs : sigma ≠ 0) (_hs1 : sigma ≠ 1) :
    ordinaryTransverse k ↔ (∀ i, WLMinors.ordinaryTransverse rho mu sigma z k i = 0) := by
  rw [ordinaryTransverse_geometry, WLMinors.ordinaryTransverse_zero_iff rho mu sigma z k hs,
    CAS.ordinaryTransverse_locus mu sigma k hm _hs1]

theorem ordinaryTransverse_all_minors (rho mu sigma : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : sigma ≠ 0) (hs1 : sigma ≠ 1) :
    ordinaryTransverse k ↔ AllOrderedMinorsZero
      (rootStack rho mu sigma k 1) 2 := by
  rw [ordinaryTransverse_minors rho mu sigma 0 k hm hs hs1,
    WLMinors.ordinaryTransverse_complete rho mu sigma 0 k hr hs]

def extraRank (k : Vec 3) : Prop := ((k 1 = 0 ∧ k 2 = 0))

theorem extraRank_geometry (k : Vec 3) : extraRank k ↔ (k 1 = 0 ∧ k 2 = 0) := by
  classical
  locus_eval extraRank

theorem extraRank_minors (rho mu sigma z : ℝ) (k : Vec 3)
    (hm : mu ≠ 0) (hs : sigma ≠ 0) (_hs1 : sigma ≠ 1) :
    extraRank k ↔ (∀ i, WLMinors.extraRank rho mu sigma z k i = 0) := by
  rw [extraRank_geometry, WLMinors.extraRank_zero_iff rho mu sigma z k hs,
    CAS.extraRank_locus mu sigma k hm hs _hs1]

theorem extraRank_all_minors (rho mu sigma : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : sigma ≠ 0) (hs1 : sigma ≠ 1) :
    extraRank k ↔ AllOrderedMinorsZero
      (rootMatrix rho mu sigma k 2) 2 := by
  rw [extraRank_minors rho mu sigma 0 k hm hs hs1,
    WLMinors.extraRank_complete rho mu sigma 0 k hr hs]

def extraTransverse (k : Vec 3) : Prop := ((k 0 = 0 ∧ (0 < k 1 ∨ k 1 < 0)) ∨ (k 0 = 0 ∧ (0 < k 2 ∨ k 2 < 0) ∧ k 1 = 0 ∧ (0 < k 2 ∨ k 2 < 0)) ∨ (k 1 = 0 ∧ k 2 = 0))

theorem extraTransverse_geometry (k : Vec 3) : extraTransverse k ↔ (k 0 = 0 ∨ (k 1 = 0 ∧ k 2 = 0)) := by
  classical
  locus_eval extraTransverse

theorem extraTransverse_minors (rho mu sigma z : ℝ) (k : Vec 3)
    (hm : mu ≠ 0) (hs : sigma ≠ 0) (_hs1 : sigma ≠ 1) :
    extraTransverse k ↔ (∀ i, WLMinors.extraTransverse rho mu sigma z k i = 0) := by
  rw [extraTransverse_geometry, WLMinors.extraTransverse_zero_iff rho mu sigma z k hs,
    CAS.extraTransverse_locus mu sigma k hm hs _hs1]

theorem extraTransverse_all_minors (rho mu sigma : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : sigma ≠ 0) (hs1 : sigma ≠ 1) :
    extraTransverse k ↔ AllOrderedMinorsZero
      (rootStack rho mu sigma k 2) 3 := by
  rw [extraTransverse_minors rho mu sigma 0 k hm hs hs1,
    WLMinors.extraTransverse_complete rho mu sigma 0 k hr hs]

def parallelPoint : Vec 3 := ![(27 / 1 : ℝ), (0 / 1 : ℝ), (0 / 1 : ℝ)]
theorem parallelPoint_nonzero : parallelPoint ≠ 0 := by
  intro h
  have hv := congrFun h 0
  norm_num [parallelPoint] at hv
theorem parallelPoint_target : extraTransverse parallelPoint := by
  rw [extraTransverse_geometry]; norm_num [parallelPoint, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]
theorem parallelPoint_separate : ordinaryRank parallelPoint ∧ parallelPoint 0 ≠ 0 := by
  rw [ordinaryRank_geometry]; norm_num [parallelPoint, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]

def perpendicularPoint : Vec 3 := ![(0 / 1 : ℝ), (27 / 1 : ℝ), (-1 / 2 : ℝ)]
theorem perpendicularPoint_nonzero : perpendicularPoint ≠ 0 := by
  intro h
  have hv := congrFun h 1
  norm_num [perpendicularPoint] at hv
theorem perpendicularPoint_target : extraTransverse perpendicularPoint := by
  rw [extraTransverse_geometry]; norm_num [perpendicularPoint, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]
theorem perpendicularPoint_separate : ¬ ordinaryRank perpendicularPoint ∧ perpendicularPoint 0 = 0 := by
  rw [ordinaryRank_geometry]; norm_num [perpendicularPoint, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]

end
end S10Audit.CAS.WLLoci
