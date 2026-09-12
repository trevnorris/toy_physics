import S10Audit.CAS.RecordSupport
import S10Audit.CAS.PYRecordsGeneric
import S10Audit.CAS.PYRecordsParallel
import S10Audit.CAS.PYRecordsPerpendicular
import S10Audit.CAS.WLRecordsGeneric
import S10Audit.CAS.WLRecordsParallel
import S10Audit.CAS.WLRecordsPerpendicular

set_option backward.isDefEq.respectTransparency false

namespace S10Audit.CAS.Records
open S10Pilot S10Anisotropic Polynomial
noncomputable section

def PY_generic_root0_sign : ReportedSign := .zero
theorem PY_generic_root0_sign_observed : PY_generic_root0_sign = .zero := rfl

theorem PY_generic_root0_sign_computed (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) (_hk : _k ≠ 0) :
    (classifiedSign 0).Holds (Expr.eval (values _rho _mu _sigma z _k) PYRootsGeneric.n5) := by
  have hr := _hd.1
  have hm := _hd.2.1
  have hs := _hd.2.2.1
  rw [PYRootsGeneric.q3_roots_distinct_cell0 _rho _mu _sigma z _k ]
  exact reference_root_sign _rho _mu _sigma _k 0 hr hm hs _hk

theorem PY_generic_root0_sign_sound (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) (_hk : _k ≠ 0) :
    PY_generic_root0_sign.Holds (Expr.eval (values _rho _mu _sigma z _k) PYRootsGeneric.n5) := by
  exact PY_generic_root0_sign_computed _rho _mu _sigma _lam z _k _hd _hk

def PY_generic_root1_sign : ReportedSign := .positive
theorem PY_generic_root1_sign_observed : PY_generic_root1_sign = .positive := rfl

theorem PY_generic_root1_sign_computed (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) (_hk : _k ≠ 0) :
    (classifiedSign 1).Holds (Expr.eval (values _rho _mu _sigma z _k) PYRootsGeneric.n17) := by
  have hr := _hd.1
  have hm := _hd.2.1
  have hs := _hd.2.2.1
  rw [PYRootsGeneric.q3_roots_distinct_cell1 _rho _mu _sigma z _k (by positivity)]
  exact reference_root_sign _rho _mu _sigma _k 1 hr hm hs _hk

theorem PY_generic_root1_sign_sound (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) (_hk : _k ≠ 0) :
    PY_generic_root1_sign.Holds (Expr.eval (values _rho _mu _sigma z _k) PYRootsGeneric.n17) := by
  exact PY_generic_root1_sign_computed _rho _mu _sigma _lam z _k _hd _hk

def PY_generic_root2_sign : ReportedSign := .undecided
theorem PY_generic_root2_sign_observed : PY_generic_root2_sign = .undecided := rfl

theorem PY_generic_root2_sign_computed (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) (_hk : _k ≠ 0) :
    (classifiedSign 2).Holds (Expr.eval (values _rho _mu _sigma z _k) PYRootsGeneric.n24) := by
  have hr := _hd.1
  have hm := _hd.2.1
  have hs := _hd.2.2.1
  rw [PYRootsGeneric.q3_roots_distinct_cell2 _rho _mu _sigma z _k (by positivity) (by positivity)]
  exact reference_root_sign _rho _mu _sigma _k 2 hr hm hs _hk

theorem PY_generic_root2_sign_sound (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) (_hk : _k ≠ 0) :
    PY_generic_root2_sign.Holds (Expr.eval (values _rho _mu _sigma z _k) PYRootsGeneric.n24) := by
  trivial

def PY_generic_spectrum_returned : Bool := true
theorem PY_generic_spectrum_returned_nonempty (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_spectrum_returned = true ↔ PYRootsGeneric.candidates _rho _mu _sigma z _k ≠ [] := by
  rw [PYRootsGeneric.candidates_reference _rho _mu _sigma z _k (ne_of_gt _hd.1) (ne_of_gt _hd.2.2.1)]
  simp [PY_generic_spectrum_returned]

theorem PY_generic_spectrum_returned_complete (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    Expr.eval (values _rho _mu _sigma z _k) PYRecordsGeneric.n25 = 0 ↔
      z ∈ PYRootsGeneric.candidates _rho _mu _sigma z _k := by
  rw [PYRecordsGeneric.q3_spectrum_solve_condition_operands_cell0 _rho _mu _sigma z _k]
  rw [rootPolynomial_complete _rho _mu _sigma z _k (ne_of_gt _hd.1) (ne_of_gt _hd.2.2.1),
    rootPolynomial_roots _rho _mu _sigma _k (ne_of_gt _hd.1) (ne_of_gt _hd.2.2.1),
    PYRootsGeneric.candidates_reference _rho _mu _sigma z _k (ne_of_gt _hd.1) (ne_of_gt _hd.2.2.1)]
  simp [referenceRoot]

def PY_parallel_root0_sign : ReportedSign := .zero
theorem PY_parallel_root0_sign_observed : PY_parallel_root0_sign = .zero := rfl

theorem PY_parallel_root0_sign_computed (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam)  :
    (classifiedSign 0).Holds (Expr.eval (values _rho _mu _sigma z _k) PYRootsParallel.n4) := by
  have hr := _hd.1
  have hm := _hd.2.1
  have hs := _hd.2.2.1
  rw [PYRootsParallel.q8_stratum1_q3_roots_distinct_cell0 _rho _mu _sigma z _k ]
  exact reference_root_sign _rho _mu _sigma PYLoci.parallelPoint 0 hr hm hs PYLoci.parallelPoint_nonzero

theorem PY_parallel_root0_sign_sound (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam)  :
    PY_parallel_root0_sign.Holds (Expr.eval (values _rho _mu _sigma z _k) PYRootsParallel.n4) := by
  exact PY_parallel_root0_sign_computed _rho _mu _sigma _lam z _k _hd

def PY_parallel_root1_sign : ReportedSign := .positive
theorem PY_parallel_root1_sign_observed : PY_parallel_root1_sign = .positive := rfl

theorem PY_parallel_root1_sign_computed (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam)  :
    (classifiedSign 1).Holds (Expr.eval (values _rho _mu _sigma z _k) PYRootsParallel.n3) := by
  have hr := _hd.1
  have hm := _hd.2.1
  have hs := _hd.2.2.1
  rw [PYRootsParallel.q8_stratum1_q3_roots_distinct_cell1 _rho _mu _sigma z _k (by positivity)]
  exact reference_root_sign _rho _mu _sigma PYLoci.parallelPoint 1 hr hm hs PYLoci.parallelPoint_nonzero

theorem PY_parallel_root1_sign_sound (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam)  :
    PY_parallel_root1_sign.Holds (Expr.eval (values _rho _mu _sigma z _k) PYRootsParallel.n3) := by
  exact PY_parallel_root1_sign_computed _rho _mu _sigma _lam z _k _hd

def PY_parallel_spectrum_returned : Bool := true
theorem PY_parallel_spectrum_returned_nonempty (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_parallel_spectrum_returned = true ↔ PYRootsParallel.candidates _rho _mu _sigma z _k ≠ [] := by
  rw [PYRootsParallel.candidates_reference _rho _mu _sigma z _k (ne_of_gt _hd.1) (ne_of_gt _hd.2.2.1)]
  simp [PY_parallel_spectrum_returned]

theorem PY_parallel_spectrum_returned_complete (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    Expr.eval (values _rho _mu _sigma z _k) PYRecordsParallel.n10 = 0 ↔
      z ∈ PYRootsParallel.candidates _rho _mu _sigma z _k := by
  rw [PYRecordsParallel.q8_stratum1_q3_spectrum_solve_condition_operands_cell0 _rho _mu _sigma z _k]
  exact PYRootsParallel.candidates_complete _rho _mu _sigma z _k (ne_of_gt _hd.1) (ne_of_gt _hd.2.1)
    _hd.2.2.1 _hd.2.2.2.1 z

def PY_perpendicular_root0_sign : ReportedSign := .zero
theorem PY_perpendicular_root0_sign_observed : PY_perpendicular_root0_sign = .zero := rfl

theorem PY_perpendicular_root0_sign_computed (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam)  :
    (classifiedSign 0).Holds (Expr.eval (values _rho _mu _sigma z _k) PYRootsPerpendicular.n4) := by
  have hr := _hd.1
  have hm := _hd.2.1
  have hs := _hd.2.2.1
  rw [PYRootsPerpendicular.q8_stratum2_q3_roots_distinct_cell0 _rho _mu _sigma z _k ]
  exact reference_root_sign _rho _mu _sigma PYLoci.perpendicularPoint 0 hr hm hs PYLoci.perpendicularPoint_nonzero

theorem PY_perpendicular_root0_sign_sound (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam)  :
    PY_perpendicular_root0_sign.Holds (Expr.eval (values _rho _mu _sigma z _k) PYRootsPerpendicular.n4) := by
  exact PY_perpendicular_root0_sign_computed _rho _mu _sigma _lam z _k _hd

def PY_perpendicular_root1_sign : ReportedSign := .positive
theorem PY_perpendicular_root1_sign_observed : PY_perpendicular_root1_sign = .positive := rfl

theorem PY_perpendicular_root1_sign_computed (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam)  :
    (classifiedSign 1).Holds (Expr.eval (values _rho _mu _sigma z _k) PYRootsPerpendicular.n3) := by
  have hr := _hd.1
  have hm := _hd.2.1
  have hs := _hd.2.2.1
  rw [PYRootsPerpendicular.q8_stratum2_q3_roots_distinct_cell1 _rho _mu _sigma z _k (by positivity)]
  exact reference_root_sign _rho _mu _sigma PYLoci.perpendicularPoint 1 hr hm hs PYLoci.perpendicularPoint_nonzero

theorem PY_perpendicular_root1_sign_sound (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam)  :
    PY_perpendicular_root1_sign.Holds (Expr.eval (values _rho _mu _sigma z _k) PYRootsPerpendicular.n3) := by
  exact PY_perpendicular_root1_sign_computed _rho _mu _sigma _lam z _k _hd

def PY_perpendicular_root2_sign : ReportedSign := .undecided
theorem PY_perpendicular_root2_sign_observed : PY_perpendicular_root2_sign = .undecided := rfl

theorem PY_perpendicular_root2_sign_computed (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam)  :
    (classifiedSign 2).Holds (Expr.eval (values _rho _mu _sigma z _k) PYRootsPerpendicular.n7) := by
  have hr := _hd.1
  have hm := _hd.2.1
  have hs := _hd.2.2.1
  rw [PYRootsPerpendicular.q8_stratum2_q3_roots_distinct_cell2 _rho _mu _sigma z _k (by positivity) (by positivity)]
  exact reference_root_sign _rho _mu _sigma PYLoci.perpendicularPoint 2 hr hm hs PYLoci.perpendicularPoint_nonzero

theorem PY_perpendicular_root2_sign_sound (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam)  :
    PY_perpendicular_root2_sign.Holds (Expr.eval (values _rho _mu _sigma z _k) PYRootsPerpendicular.n7) := by
  trivial

def PY_perpendicular_spectrum_returned : Bool := true
theorem PY_perpendicular_spectrum_returned_nonempty (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_perpendicular_spectrum_returned = true ↔ PYRootsPerpendicular.candidates _rho _mu _sigma z _k ≠ [] := by
  rw [PYRootsPerpendicular.candidates_reference _rho _mu _sigma z _k (ne_of_gt _hd.1) (ne_of_gt _hd.2.2.1)]
  simp [PY_perpendicular_spectrum_returned]

theorem PY_perpendicular_spectrum_returned_complete (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    Expr.eval (values _rho _mu _sigma z _k) PYRecordsPerpendicular.n11 = 0 ↔
      z ∈ PYRootsPerpendicular.candidates _rho _mu _sigma z _k := by
  rw [PYRecordsPerpendicular.q8_stratum2_q3_spectrum_solve_condition_operands_cell0 _rho _mu _sigma z _k]
  exact PYRootsPerpendicular.candidates_complete _rho _mu _sigma z _k (ne_of_gt _hd.1) (ne_of_gt _hd.2.1)
    _hd.2.2.1 _hd.2.2.2.1 z

def WL_generic_root0_sign : ReportedSign := .zero
theorem WL_generic_root0_sign_observed : WL_generic_root0_sign = .zero := rfl

theorem WL_generic_root0_sign_computed (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) (_hk : _k ≠ 0) :
    (classifiedSign 0).Holds (Expr.eval (values _rho _mu _sigma z _k) WLRootsGeneric.n5) := by
  have hr := _hd.1
  have hm := _hd.2.1
  have hs := _hd.2.2.1
  rw [WLRootsGeneric.q3_distinct_roots_cell0 _rho _mu _sigma z _k ]
  exact reference_root_sign _rho _mu _sigma _k 0 hr hm hs _hk

theorem WL_generic_root0_sign_sound (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) (_hk : _k ≠ 0) :
    WL_generic_root0_sign.Holds (Expr.eval (values _rho _mu _sigma z _k) WLRootsGeneric.n5) := by
  exact WL_generic_root0_sign_computed _rho _mu _sigma _lam z _k _hd _hk

def WL_generic_root0_conditions : List Prop := []
theorem WL_generic_root0_conditions_empty : WL_generic_root0_conditions = [] := rfl

theorem WL_generic_root0_conditions_root (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    Matrix.det ((1/2 : ℝ) • referenceMatrix _rho _mu _sigma (Expr.eval (values _rho _mu _sigma z _k) WLRootsGeneric.n5) _k) = 0 := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [WLRootsGeneric.q3_distinct_roots_cell0 _rho _mu _sigma z _k ]
  exact reference_root_isRoot _rho _mu _sigma _k 0 (ne_of_gt hr) (ne_of_gt hs)

def WL_generic_root1_sign : ReportedSign := .positive
theorem WL_generic_root1_sign_observed : WL_generic_root1_sign = .positive := rfl

theorem WL_generic_root1_sign_computed (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) (_hk : _k ≠ 0) :
    (classifiedSign 1).Holds (Expr.eval (values _rho _mu _sigma z _k) WLRootsGeneric.n17) := by
  have hr := _hd.1
  have hm := _hd.2.1
  have hs := _hd.2.2.1
  rw [WLRootsGeneric.q3_distinct_roots_cell1 _rho _mu _sigma z _k (by positivity)]
  exact reference_root_sign _rho _mu _sigma _k 1 hr hm hs _hk

theorem WL_generic_root1_sign_sound (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) (_hk : _k ≠ 0) :
    WL_generic_root1_sign.Holds (Expr.eval (values _rho _mu _sigma z _k) WLRootsGeneric.n17) := by
  exact WL_generic_root1_sign_computed _rho _mu _sigma _lam z _k _hd _hk

def WL_generic_root1_conditions : List Prop := []
theorem WL_generic_root1_conditions_empty : WL_generic_root1_conditions = [] := rfl

theorem WL_generic_root1_conditions_root (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    Matrix.det ((1/2 : ℝ) • referenceMatrix _rho _mu _sigma (Expr.eval (values _rho _mu _sigma z _k) WLRootsGeneric.n17) _k) = 0 := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [WLRootsGeneric.q3_distinct_roots_cell1 _rho _mu _sigma z _k (by positivity)]
  exact reference_root_isRoot _rho _mu _sigma _k 1 (ne_of_gt hr) (ne_of_gt hs)

def WL_generic_root2_sign : ReportedSign := .positive
theorem WL_generic_root2_sign_observed : WL_generic_root2_sign = .positive := rfl

theorem WL_generic_root2_sign_computed (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) (_hk : _k ≠ 0) :
    (classifiedSign 2).Holds (Expr.eval (values _rho _mu _sigma z _k) WLRootsGeneric.n24) := by
  have hr := _hd.1
  have hm := _hd.2.1
  have hs := _hd.2.2.1
  rw [WLRootsGeneric.q3_distinct_roots_cell2 _rho _mu _sigma z _k (by positivity) (by positivity)]
  exact reference_root_sign _rho _mu _sigma _k 2 hr hm hs _hk

theorem WL_generic_root2_sign_sound (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) (_hk : _k ≠ 0) :
    WL_generic_root2_sign.Holds (Expr.eval (values _rho _mu _sigma z _k) WLRootsGeneric.n24) := by
  exact WL_generic_root2_sign_computed _rho _mu _sigma _lam z _k _hd _hk

def WL_generic_root2_conditions : List Prop := []
theorem WL_generic_root2_conditions_empty : WL_generic_root2_conditions = [] := rfl

theorem WL_generic_root2_conditions_root (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    Matrix.det ((1/2 : ℝ) • referenceMatrix _rho _mu _sigma (Expr.eval (values _rho _mu _sigma z _k) WLRootsGeneric.n24) _k) = 0 := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [WLRootsGeneric.q3_distinct_roots_cell2 _rho _mu _sigma z _k (by positivity) (by positivity)]
  exact reference_root_isRoot _rho _mu _sigma _k 2 (ne_of_gt hr) (ne_of_gt hs)

def WL_parallel_root0_sign : ReportedSign := .zero
theorem WL_parallel_root0_sign_observed : WL_parallel_root0_sign = .zero := rfl

theorem WL_parallel_root0_sign_computed (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam)  :
    (classifiedSign 0).Holds (Expr.eval (values _rho _mu _sigma z _k) WLRootsParallel.n4) := by
  have hr := _hd.1
  have hm := _hd.2.1
  have hs := _hd.2.2.1
  rw [WLRootsParallel.stratum1_q3_distinct_roots_cell0 _rho _mu _sigma z _k ]
  exact reference_root_sign _rho _mu _sigma WLLoci.parallelPoint 0 hr hm hs WLLoci.parallelPoint_nonzero

theorem WL_parallel_root0_sign_sound (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam)  :
    WL_parallel_root0_sign.Holds (Expr.eval (values _rho _mu _sigma z _k) WLRootsParallel.n4) := by
  exact WL_parallel_root0_sign_computed _rho _mu _sigma _lam z _k _hd

def WL_parallel_root0_conditions : List Prop := []
theorem WL_parallel_root0_conditions_empty : WL_parallel_root0_conditions = [] := rfl

theorem WL_parallel_root0_conditions_root (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    Matrix.det ((1/2 : ℝ) • referenceMatrix _rho _mu _sigma (Expr.eval (values _rho _mu _sigma z _k) WLRootsParallel.n4) WLLoci.parallelPoint) = 0 := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [WLRootsParallel.stratum1_q3_distinct_roots_cell0 _rho _mu _sigma z _k ]
  exact reference_root_isRoot _rho _mu _sigma WLLoci.parallelPoint 0 (ne_of_gt hr) (ne_of_gt hs)

def WL_parallel_root1_sign : ReportedSign := .positive
theorem WL_parallel_root1_sign_observed : WL_parallel_root1_sign = .positive := rfl

theorem WL_parallel_root1_sign_computed (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam)  :
    (classifiedSign 1).Holds (Expr.eval (values _rho _mu _sigma z _k) WLRootsParallel.n7) := by
  have hr := _hd.1
  have hm := _hd.2.1
  have hs := _hd.2.2.1
  rw [WLRootsParallel.stratum1_q3_distinct_roots_cell1 _rho _mu _sigma z _k (by positivity)]
  exact reference_root_sign _rho _mu _sigma WLLoci.parallelPoint 1 hr hm hs WLLoci.parallelPoint_nonzero

theorem WL_parallel_root1_sign_sound (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam)  :
    WL_parallel_root1_sign.Holds (Expr.eval (values _rho _mu _sigma z _k) WLRootsParallel.n7) := by
  exact WL_parallel_root1_sign_computed _rho _mu _sigma _lam z _k _hd

def WL_parallel_root1_conditions : List Prop := []
theorem WL_parallel_root1_conditions_empty : WL_parallel_root1_conditions = [] := rfl

theorem WL_parallel_root1_conditions_root (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    Matrix.det ((1/2 : ℝ) • referenceMatrix _rho _mu _sigma (Expr.eval (values _rho _mu _sigma z _k) WLRootsParallel.n7) WLLoci.parallelPoint) = 0 := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [WLRootsParallel.stratum1_q3_distinct_roots_cell1 _rho _mu _sigma z _k (by positivity)]
  exact reference_root_isRoot _rho _mu _sigma WLLoci.parallelPoint 1 (ne_of_gt hr) (ne_of_gt hs)

def WL_perpendicular_root0_sign : ReportedSign := .zero
theorem WL_perpendicular_root0_sign_observed : WL_perpendicular_root0_sign = .zero := rfl

theorem WL_perpendicular_root0_sign_computed (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam)  :
    (classifiedSign 0).Holds (Expr.eval (values _rho _mu _sigma z _k) WLRootsPerpendicular.n4) := by
  have hr := _hd.1
  have hm := _hd.2.1
  have hs := _hd.2.2.1
  rw [WLRootsPerpendicular.stratum2_q3_distinct_roots_cell0 _rho _mu _sigma z _k ]
  exact reference_root_sign _rho _mu _sigma WLLoci.perpendicularPoint 0 hr hm hs WLLoci.perpendicularPoint_nonzero

theorem WL_perpendicular_root0_sign_sound (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam)  :
    WL_perpendicular_root0_sign.Holds (Expr.eval (values _rho _mu _sigma z _k) WLRootsPerpendicular.n4) := by
  exact WL_perpendicular_root0_sign_computed _rho _mu _sigma _lam z _k _hd

def WL_perpendicular_root0_conditions : List Prop := []
theorem WL_perpendicular_root0_conditions_empty : WL_perpendicular_root0_conditions = [] := rfl

theorem WL_perpendicular_root0_conditions_root (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    Matrix.det ((1/2 : ℝ) • referenceMatrix _rho _mu _sigma (Expr.eval (values _rho _mu _sigma z _k) WLRootsPerpendicular.n4) WLLoci.perpendicularPoint) = 0 := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [WLRootsPerpendicular.stratum2_q3_distinct_roots_cell0 _rho _mu _sigma z _k ]
  exact reference_root_isRoot _rho _mu _sigma WLLoci.perpendicularPoint 0 (ne_of_gt hr) (ne_of_gt hs)

def WL_perpendicular_root1_sign : ReportedSign := .positive
theorem WL_perpendicular_root1_sign_observed : WL_perpendicular_root1_sign = .positive := rfl

theorem WL_perpendicular_root1_sign_computed (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam)  :
    (classifiedSign 1).Holds (Expr.eval (values _rho _mu _sigma z _k) WLRootsPerpendicular.n9) := by
  have hr := _hd.1
  have hm := _hd.2.1
  have hs := _hd.2.2.1
  rw [WLRootsPerpendicular.stratum2_q3_distinct_roots_cell1 _rho _mu _sigma z _k (by positivity)]
  exact reference_root_sign _rho _mu _sigma WLLoci.perpendicularPoint 1 hr hm hs WLLoci.perpendicularPoint_nonzero

theorem WL_perpendicular_root1_sign_sound (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam)  :
    WL_perpendicular_root1_sign.Holds (Expr.eval (values _rho _mu _sigma z _k) WLRootsPerpendicular.n9) := by
  exact WL_perpendicular_root1_sign_computed _rho _mu _sigma _lam z _k _hd

def WL_perpendicular_root1_conditions : List Prop := []
theorem WL_perpendicular_root1_conditions_empty : WL_perpendicular_root1_conditions = [] := rfl

theorem WL_perpendicular_root1_conditions_root (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    Matrix.det ((1/2 : ℝ) • referenceMatrix _rho _mu _sigma (Expr.eval (values _rho _mu _sigma z _k) WLRootsPerpendicular.n9) WLLoci.perpendicularPoint) = 0 := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [WLRootsPerpendicular.stratum2_q3_distinct_roots_cell1 _rho _mu _sigma z _k (by positivity)]
  exact reference_root_isRoot _rho _mu _sigma WLLoci.perpendicularPoint 1 (ne_of_gt hr) (ne_of_gt hs)

def WL_perpendicular_root2_sign : ReportedSign := .positive
theorem WL_perpendicular_root2_sign_observed : WL_perpendicular_root2_sign = .positive := rfl

theorem WL_perpendicular_root2_sign_computed (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam)  :
    (classifiedSign 2).Holds (Expr.eval (values _rho _mu _sigma z _k) WLRootsPerpendicular.n12) := by
  have hr := _hd.1
  have hm := _hd.2.1
  have hs := _hd.2.2.1
  rw [WLRootsPerpendicular.stratum2_q3_distinct_roots_cell2 _rho _mu _sigma z _k (by positivity) (by positivity)]
  exact reference_root_sign _rho _mu _sigma WLLoci.perpendicularPoint 2 hr hm hs WLLoci.perpendicularPoint_nonzero

theorem WL_perpendicular_root2_sign_sound (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam)  :
    WL_perpendicular_root2_sign.Holds (Expr.eval (values _rho _mu _sigma z _k) WLRootsPerpendicular.n12) := by
  exact WL_perpendicular_root2_sign_computed _rho _mu _sigma _lam z _k _hd

def WL_perpendicular_root2_conditions : List Prop := []
theorem WL_perpendicular_root2_conditions_empty : WL_perpendicular_root2_conditions = [] := rfl

theorem WL_perpendicular_root2_conditions_root (_rho _mu _sigma _lam z : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    Matrix.det ((1/2 : ℝ) • referenceMatrix _rho _mu _sigma (Expr.eval (values _rho _mu _sigma z _k) WLRootsPerpendicular.n12) WLLoci.perpendicularPoint) = 0 := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [WLRootsPerpendicular.stratum2_q3_distinct_roots_cell2 _rho _mu _sigma z _k (by positivity) (by positivity)]
  exact reference_root_isRoot _rho _mu _sigma WLLoci.perpendicularPoint 2 (ne_of_gt hr) (ne_of_gt hs)

def PY_parallel_not_skipped : Bool := true
theorem PY_parallel_not_skipped_allowed :
    PY_parallel_not_skipped = true ∧ PYLoci.parallelPoint ≠ 0 ∧ PYLoci.extraTransverse PYLoci.parallelPoint :=
  ⟨rfl, PYLoci.parallelPoint_nonzero, PYLoci.parallelPoint_target⟩

def PY_perpendicular_not_skipped : Bool := true
theorem PY_perpendicular_not_skipped_allowed :
    PY_perpendicular_not_skipped = true ∧ PYLoci.perpendicularPoint ≠ 0 ∧ PYLoci.extraTransverse PYLoci.perpendicularPoint :=
  ⟨rfl, PYLoci.perpendicularPoint_nonzero, PYLoci.perpendicularPoint_target⟩

def PY_skipped_branch (_k : Vec 3) : Prop := (((_k 0 = 0) ∧ (_k 1 = 0)) ∧ (_k 2 = 0))
theorem PY_skipped_branch_geometry (_k : Vec 3) : PY_skipped_branch _k ↔ _k = 0 := by
  simp only [PY_skipped_branch, wavevector_zero_iff]
  tauto

def PY_skipped_test : Bool := false
theorem PY_skipped_test_decides : PY_skipped_test = true ↔
    ∃ _k : Vec 3, PY_skipped_branch _k ∧ 0 < normSq _k := by
  simp [PY_skipped_test, PY_skipped_branch_geometry, normSq_positive_iff]

theorem PY_skipped_reason_sound (_k : Vec 3) (h : PY_skipped_branch _k) : normSq _k = 0 := by
  rw [PY_skipped_branch_geometry] at h
  subst _k
  simp [normSq, dot]

def PY_skipped_aggregate_branch (_k : Vec 3) : Prop := (((_k 0 = 0) ∧ (_k 1 = 0)) ∧ (_k 2 = 0))
def PY_skipped_aggregate_test : Bool := false
theorem PY_skipped_aggregate_branch_primary (_k : Vec 3) :
    PY_skipped_aggregate_branch _k ↔ PY_skipped_branch _k := Iff.rfl
theorem PY_skipped_aggregate_test_decides : PY_skipped_aggregate_test = true ↔
    ∃ _k : Vec 3, PY_skipped_aggregate_branch _k ∧ 0 < normSq _k := by
  simp [PY_skipped_aggregate_test, PY_skipped_aggregate_branch_primary,
    PY_skipped_branch_geometry, normSq_positive_iff]
theorem PY_skipped_aggregate_reason_sound (_k : Vec 3) (h : PY_skipped_aggregate_branch _k) : normSq _k = 0 :=
  PY_skipped_reason_sound _k ((PY_skipped_aggregate_branch_primary _k).mp h)

end
end S10Audit.CAS.Records
