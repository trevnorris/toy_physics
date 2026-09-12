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

def PY_q8root1_locus0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := (((_k 0 = 0) ∧ (_k 1 = 0)) ∧ (_k 2 = 0))

theorem PY_q8root1_locus0_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : PY_q8root1_locus0 _rho _mu _sigma _lam _a _k ↔ Coincidence.PY_generic_locus0 _rho _mu _sigma _lam _a _k := by
  rfl

theorem PY_q8root1_locus0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : PY_q8root1_locus0 _rho _mu _sigma _lam _a _k ↔ coincidenceGeometry 0 _k := by
  exact (PY_q8root1_locus0_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.PY_generic_locus0_correct _rho _mu _sigma _lam _a _k _hd)

def PY_q8root1_locus1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := (((_k 0 = 0) ∧ (_k 1 = 0)) ∧ (_k 2 = 0))

theorem PY_q8root1_locus1_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : PY_q8root1_locus1 _rho _mu _sigma _lam _a _k ↔ Coincidence.PY_generic_locus1 _rho _mu _sigma _lam _a _k := by
  rfl

theorem PY_q8root1_locus1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : PY_q8root1_locus1 _rho _mu _sigma _lam _a _k ↔ coincidenceGeometry 1 _k := by
  exact (PY_q8root1_locus1_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.PY_generic_locus1_correct _rho _mu _sigma _lam _a _k _hd)

def PY_q8root1_locus2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((_k 1 = 0) ∧ (_k 2 = 0))

theorem PY_q8root1_locus2_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : PY_q8root1_locus2 _rho _mu _sigma _lam _a _k ↔ Coincidence.PY_generic_locus2 _rho _mu _sigma _lam _a _k := by
  rfl

theorem PY_q8root1_locus2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : PY_q8root1_locus2 _rho _mu _sigma _lam _a _k ↔ coincidenceGeometry 2 _k := by
  exact (PY_q8root1_locus2_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.PY_generic_locus2_correct _rho _mu _sigma _lam _a _k _hd)

def PY_q8root1_allowed0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := (((((_k 0 = 0) ∧ (_k 1 = 0)) ∧ (_k 2 = 0)) ∧ ((((_k 0 ^ 2) + (_k 1 ^ 2)) + (_k 2 ^ 2)) > 0)) ∧ (((((((((((((True ∧ (0 < (3 : ℝ))) ∧ (0 < _lam)) ∧ (0 < _mu)) ∧ (0 < _rho)) ∧ (0 < _sigma)) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ (_sigma ≠ 1)) ∧ (0 < (((_k 0 ^ 2) + (_k 1 ^ 2)) + (_k 2 ^ 2)))))

theorem PY_q8root1_allowed0_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : PY_q8root1_allowed0 _rho _mu _sigma _lam _a _k ↔ Coincidence.PY_generic_allowed0 _rho _mu _sigma _lam _a _k := by
  rfl

theorem PY_q8root1_allowed0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : PY_q8root1_allowed0 _rho _mu _sigma _lam _a _k ↔ allowedCoincidence 0 _k := by
  exact (PY_q8root1_allowed0_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.PY_generic_allowed0_correct _rho _mu _sigma _lam _a _k _hd)

def PY_q8root1_allowed1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := (((((_k 0 = 0) ∧ (_k 1 = 0)) ∧ (_k 2 = 0)) ∧ ((((_k 0 ^ 2) + (_k 1 ^ 2)) + (_k 2 ^ 2)) > 0)) ∧ (((((((((((((True ∧ (0 < (3 : ℝ))) ∧ (0 < _lam)) ∧ (0 < _mu)) ∧ (0 < _rho)) ∧ (0 < _sigma)) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ (_sigma ≠ 1)) ∧ (0 < (((_k 0 ^ 2) + (_k 1 ^ 2)) + (_k 2 ^ 2)))))

theorem PY_q8root1_allowed1_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : PY_q8root1_allowed1 _rho _mu _sigma _lam _a _k ↔ Coincidence.PY_generic_allowed1 _rho _mu _sigma _lam _a _k := by
  rfl

theorem PY_q8root1_allowed1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : PY_q8root1_allowed1 _rho _mu _sigma _lam _a _k ↔ allowedCoincidence 1 _k := by
  exact (PY_q8root1_allowed1_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.PY_generic_allowed1_correct _rho _mu _sigma _lam _a _k _hd)

def PY_q8root1_allowed2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((((_k 1 = 0) ∧ (_k 2 = 0)) ∧ ((((_k 0 ^ 2) + (_k 1 ^ 2)) + (_k 2 ^ 2)) > 0)) ∧ (((((((((((((True ∧ (0 < (3 : ℝ))) ∧ (0 < _lam)) ∧ (0 < _mu)) ∧ (0 < _rho)) ∧ (0 < _sigma)) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ (_sigma ≠ 1)) ∧ (0 < (((_k 0 ^ 2) + (_k 1 ^ 2)) + (_k 2 ^ 2)))))

theorem PY_q8root1_allowed2_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : PY_q8root1_allowed2 _rho _mu _sigma _lam _a _k ↔ Coincidence.PY_generic_allowed2 _rho _mu _sigma _lam _a _k := by
  rfl

theorem PY_q8root1_allowed2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : PY_q8root1_allowed2 _rho _mu _sigma _lam _a _k ↔ allowedCoincidence 2 _k := by
  exact (PY_q8root1_allowed2_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.PY_generic_allowed2_correct _rho _mu _sigma _lam _a _k _hd)

def PY_q8root1_status0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_q8root1_status0_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : PY_q8root1_status0 _rho _mu _sigma _lam _a _k ↔ Coincidence.PY_generic_status0 _rho _mu _sigma _lam _a _k := by
  rfl

theorem PY_q8root1_status0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : PY_q8root1_status0 _rho _mu _sigma _lam _a _k ↔ (∃ w : Vec 3, allowedCoincidence 0 w) := by
  exact (PY_q8root1_status0_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.PY_generic_status0_correct _rho _mu _sigma _lam _a _k _hd)

def PY_q8root1_status1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_q8root1_status1_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : PY_q8root1_status1 _rho _mu _sigma _lam _a _k ↔ Coincidence.PY_generic_status1 _rho _mu _sigma _lam _a _k := by
  rfl

theorem PY_q8root1_status1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : PY_q8root1_status1 _rho _mu _sigma _lam _a _k ↔ (∃ w : Vec 3, allowedCoincidence 1 w) := by
  exact (PY_q8root1_status1_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.PY_generic_status1_correct _rho _mu _sigma _lam _a _k _hd)

def PY_q8root1_status2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := True

theorem PY_q8root1_status2_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : PY_q8root1_status2 _rho _mu _sigma _lam _a _k ↔ Coincidence.PY_generic_status2 _rho _mu _sigma _lam _a _k := by
  rfl

theorem PY_q8root1_status2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : PY_q8root1_status2 _rho _mu _sigma _lam _a _k ↔ (∃ w : Vec 3, allowedCoincidence 2 w) := by
  exact (PY_q8root1_status2_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.PY_generic_status2_correct _rho _mu _sigma _lam _a _k _hd)

def PY_q8root2_locus0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := (((_k 0 = 0) ∧ (_k 1 = 0)) ∧ (_k 2 = 0))

theorem PY_q8root2_locus0_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : PY_q8root2_locus0 _rho _mu _sigma _lam _a _k ↔ Coincidence.PY_generic_locus0 _rho _mu _sigma _lam _a _k := by
  rfl

theorem PY_q8root2_locus0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : PY_q8root2_locus0 _rho _mu _sigma _lam _a _k ↔ coincidenceGeometry 0 _k := by
  exact (PY_q8root2_locus0_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.PY_generic_locus0_correct _rho _mu _sigma _lam _a _k _hd)

def PY_q8root2_locus1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := (((_k 0 = 0) ∧ (_k 1 = 0)) ∧ (_k 2 = 0))

theorem PY_q8root2_locus1_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : PY_q8root2_locus1 _rho _mu _sigma _lam _a _k ↔ Coincidence.PY_generic_locus1 _rho _mu _sigma _lam _a _k := by
  rfl

theorem PY_q8root2_locus1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : PY_q8root2_locus1 _rho _mu _sigma _lam _a _k ↔ coincidenceGeometry 1 _k := by
  exact (PY_q8root2_locus1_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.PY_generic_locus1_correct _rho _mu _sigma _lam _a _k _hd)

def PY_q8root2_locus2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((_k 1 = 0) ∧ (_k 2 = 0))

theorem PY_q8root2_locus2_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : PY_q8root2_locus2 _rho _mu _sigma _lam _a _k ↔ Coincidence.PY_generic_locus2 _rho _mu _sigma _lam _a _k := by
  rfl

theorem PY_q8root2_locus2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : PY_q8root2_locus2 _rho _mu _sigma _lam _a _k ↔ coincidenceGeometry 2 _k := by
  exact (PY_q8root2_locus2_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.PY_generic_locus2_correct _rho _mu _sigma _lam _a _k _hd)

def PY_q8root2_allowed0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := (((((_k 0 = 0) ∧ (_k 1 = 0)) ∧ (_k 2 = 0)) ∧ ((((_k 0 ^ 2) + (_k 1 ^ 2)) + (_k 2 ^ 2)) > 0)) ∧ (((((((((((((True ∧ (0 < (3 : ℝ))) ∧ (0 < _lam)) ∧ (0 < _mu)) ∧ (0 < _rho)) ∧ (0 < _sigma)) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ (_sigma ≠ 1)) ∧ (0 < (((_k 0 ^ 2) + (_k 1 ^ 2)) + (_k 2 ^ 2)))))

theorem PY_q8root2_allowed0_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : PY_q8root2_allowed0 _rho _mu _sigma _lam _a _k ↔ Coincidence.PY_generic_allowed0 _rho _mu _sigma _lam _a _k := by
  rfl

theorem PY_q8root2_allowed0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : PY_q8root2_allowed0 _rho _mu _sigma _lam _a _k ↔ allowedCoincidence 0 _k := by
  exact (PY_q8root2_allowed0_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.PY_generic_allowed0_correct _rho _mu _sigma _lam _a _k _hd)

def PY_q8root2_allowed1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := (((((_k 0 = 0) ∧ (_k 1 = 0)) ∧ (_k 2 = 0)) ∧ ((((_k 0 ^ 2) + (_k 1 ^ 2)) + (_k 2 ^ 2)) > 0)) ∧ (((((((((((((True ∧ (0 < (3 : ℝ))) ∧ (0 < _lam)) ∧ (0 < _mu)) ∧ (0 < _rho)) ∧ (0 < _sigma)) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ (_sigma ≠ 1)) ∧ (0 < (((_k 0 ^ 2) + (_k 1 ^ 2)) + (_k 2 ^ 2)))))

theorem PY_q8root2_allowed1_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : PY_q8root2_allowed1 _rho _mu _sigma _lam _a _k ↔ Coincidence.PY_generic_allowed1 _rho _mu _sigma _lam _a _k := by
  rfl

theorem PY_q8root2_allowed1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : PY_q8root2_allowed1 _rho _mu _sigma _lam _a _k ↔ allowedCoincidence 1 _k := by
  exact (PY_q8root2_allowed1_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.PY_generic_allowed1_correct _rho _mu _sigma _lam _a _k _hd)

def PY_q8root2_allowed2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((((_k 1 = 0) ∧ (_k 2 = 0)) ∧ ((((_k 0 ^ 2) + (_k 1 ^ 2)) + (_k 2 ^ 2)) > 0)) ∧ (((((((((((((True ∧ (0 < (3 : ℝ))) ∧ (0 < _lam)) ∧ (0 < _mu)) ∧ (0 < _rho)) ∧ (0 < _sigma)) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ (_sigma ≠ 1)) ∧ (0 < (((_k 0 ^ 2) + (_k 1 ^ 2)) + (_k 2 ^ 2)))))

theorem PY_q8root2_allowed2_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : PY_q8root2_allowed2 _rho _mu _sigma _lam _a _k ↔ Coincidence.PY_generic_allowed2 _rho _mu _sigma _lam _a _k := by
  rfl

theorem PY_q8root2_allowed2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : PY_q8root2_allowed2 _rho _mu _sigma _lam _a _k ↔ allowedCoincidence 2 _k := by
  exact (PY_q8root2_allowed2_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.PY_generic_allowed2_correct _rho _mu _sigma _lam _a _k _hd)

def PY_q8root2_status0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_q8root2_status0_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : PY_q8root2_status0 _rho _mu _sigma _lam _a _k ↔ Coincidence.PY_generic_status0 _rho _mu _sigma _lam _a _k := by
  rfl

theorem PY_q8root2_status0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : PY_q8root2_status0 _rho _mu _sigma _lam _a _k ↔ (∃ w : Vec 3, allowedCoincidence 0 w) := by
  exact (PY_q8root2_status0_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.PY_generic_status0_correct _rho _mu _sigma _lam _a _k _hd)

def PY_q8root2_status1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_q8root2_status1_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : PY_q8root2_status1 _rho _mu _sigma _lam _a _k ↔ Coincidence.PY_generic_status1 _rho _mu _sigma _lam _a _k := by
  rfl

theorem PY_q8root2_status1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : PY_q8root2_status1 _rho _mu _sigma _lam _a _k ↔ (∃ w : Vec 3, allowedCoincidence 1 w) := by
  exact (PY_q8root2_status1_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.PY_generic_status1_correct _rho _mu _sigma _lam _a _k _hd)

def PY_q8root2_status2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := True

theorem PY_q8root2_status2_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : PY_q8root2_status2 _rho _mu _sigma _lam _a _k ↔ Coincidence.PY_generic_status2 _rho _mu _sigma _lam _a _k := by
  rfl

theorem PY_q8root2_status2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : PY_q8root2_status2 _rho _mu _sigma _lam _a _k ↔ (∃ w : Vec 3, allowedCoincidence 2 w) := by
  exact (PY_q8root2_status2_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.PY_generic_status2_correct _rho _mu _sigma _lam _a _k _hd)

def PY_q8root3_locus0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := (((_k 0 = 0) ∧ (_k 1 = 0)) ∧ (_k 2 = 0))

theorem PY_q8root3_locus0_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : PY_q8root3_locus0 _rho _mu _sigma _lam _a _k ↔ Coincidence.PY_generic_locus0 _rho _mu _sigma _lam _a _k := by
  rfl

theorem PY_q8root3_locus0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : PY_q8root3_locus0 _rho _mu _sigma _lam _a _k ↔ coincidenceGeometry 0 _k := by
  exact (PY_q8root3_locus0_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.PY_generic_locus0_correct _rho _mu _sigma _lam _a _k _hd)

def PY_q8root3_locus1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := (((_k 0 = 0) ∧ (_k 1 = 0)) ∧ (_k 2 = 0))

theorem PY_q8root3_locus1_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : PY_q8root3_locus1 _rho _mu _sigma _lam _a _k ↔ Coincidence.PY_generic_locus1 _rho _mu _sigma _lam _a _k := by
  rfl

theorem PY_q8root3_locus1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : PY_q8root3_locus1 _rho _mu _sigma _lam _a _k ↔ coincidenceGeometry 1 _k := by
  exact (PY_q8root3_locus1_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.PY_generic_locus1_correct _rho _mu _sigma _lam _a _k _hd)

def PY_q8root3_locus2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((_k 1 = 0) ∧ (_k 2 = 0))

theorem PY_q8root3_locus2_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : PY_q8root3_locus2 _rho _mu _sigma _lam _a _k ↔ Coincidence.PY_generic_locus2 _rho _mu _sigma _lam _a _k := by
  rfl

theorem PY_q8root3_locus2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : PY_q8root3_locus2 _rho _mu _sigma _lam _a _k ↔ coincidenceGeometry 2 _k := by
  exact (PY_q8root3_locus2_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.PY_generic_locus2_correct _rho _mu _sigma _lam _a _k _hd)

def PY_q8root3_allowed0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := (((((_k 0 = 0) ∧ (_k 1 = 0)) ∧ (_k 2 = 0)) ∧ ((((_k 0 ^ 2) + (_k 1 ^ 2)) + (_k 2 ^ 2)) > 0)) ∧ (((((((((((((True ∧ (0 < (3 : ℝ))) ∧ (0 < _lam)) ∧ (0 < _mu)) ∧ (0 < _rho)) ∧ (0 < _sigma)) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ (_sigma ≠ 1)) ∧ (0 < (((_k 0 ^ 2) + (_k 1 ^ 2)) + (_k 2 ^ 2)))))

theorem PY_q8root3_allowed0_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : PY_q8root3_allowed0 _rho _mu _sigma _lam _a _k ↔ Coincidence.PY_generic_allowed0 _rho _mu _sigma _lam _a _k := by
  rfl

theorem PY_q8root3_allowed0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : PY_q8root3_allowed0 _rho _mu _sigma _lam _a _k ↔ allowedCoincidence 0 _k := by
  exact (PY_q8root3_allowed0_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.PY_generic_allowed0_correct _rho _mu _sigma _lam _a _k _hd)

def PY_q8root3_allowed1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := (((((_k 0 = 0) ∧ (_k 1 = 0)) ∧ (_k 2 = 0)) ∧ ((((_k 0 ^ 2) + (_k 1 ^ 2)) + (_k 2 ^ 2)) > 0)) ∧ (((((((((((((True ∧ (0 < (3 : ℝ))) ∧ (0 < _lam)) ∧ (0 < _mu)) ∧ (0 < _rho)) ∧ (0 < _sigma)) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ (_sigma ≠ 1)) ∧ (0 < (((_k 0 ^ 2) + (_k 1 ^ 2)) + (_k 2 ^ 2)))))

theorem PY_q8root3_allowed1_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : PY_q8root3_allowed1 _rho _mu _sigma _lam _a _k ↔ Coincidence.PY_generic_allowed1 _rho _mu _sigma _lam _a _k := by
  rfl

theorem PY_q8root3_allowed1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : PY_q8root3_allowed1 _rho _mu _sigma _lam _a _k ↔ allowedCoincidence 1 _k := by
  exact (PY_q8root3_allowed1_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.PY_generic_allowed1_correct _rho _mu _sigma _lam _a _k _hd)

def PY_q8root3_allowed2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((((_k 1 = 0) ∧ (_k 2 = 0)) ∧ ((((_k 0 ^ 2) + (_k 1 ^ 2)) + (_k 2 ^ 2)) > 0)) ∧ (((((((((((((True ∧ (0 < (3 : ℝ))) ∧ (0 < _lam)) ∧ (0 < _mu)) ∧ (0 < _rho)) ∧ (0 < _sigma)) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ (_sigma ≠ 1)) ∧ (0 < (((_k 0 ^ 2) + (_k 1 ^ 2)) + (_k 2 ^ 2)))))

theorem PY_q8root3_allowed2_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : PY_q8root3_allowed2 _rho _mu _sigma _lam _a _k ↔ Coincidence.PY_generic_allowed2 _rho _mu _sigma _lam _a _k := by
  rfl

theorem PY_q8root3_allowed2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : PY_q8root3_allowed2 _rho _mu _sigma _lam _a _k ↔ allowedCoincidence 2 _k := by
  exact (PY_q8root3_allowed2_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.PY_generic_allowed2_correct _rho _mu _sigma _lam _a _k _hd)

def PY_q8root3_status0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_q8root3_status0_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : PY_q8root3_status0 _rho _mu _sigma _lam _a _k ↔ Coincidence.PY_generic_status0 _rho _mu _sigma _lam _a _k := by
  rfl

theorem PY_q8root3_status0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : PY_q8root3_status0 _rho _mu _sigma _lam _a _k ↔ (∃ w : Vec 3, allowedCoincidence 0 w) := by
  exact (PY_q8root3_status0_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.PY_generic_status0_correct _rho _mu _sigma _lam _a _k _hd)

def PY_q8root3_status1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_q8root3_status1_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : PY_q8root3_status1 _rho _mu _sigma _lam _a _k ↔ Coincidence.PY_generic_status1 _rho _mu _sigma _lam _a _k := by
  rfl

theorem PY_q8root3_status1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : PY_q8root3_status1 _rho _mu _sigma _lam _a _k ↔ (∃ w : Vec 3, allowedCoincidence 1 w) := by
  exact (PY_q8root3_status1_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.PY_generic_status1_correct _rho _mu _sigma _lam _a _k _hd)

def PY_q8root3_status2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := True

theorem PY_q8root3_status2_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : PY_q8root3_status2 _rho _mu _sigma _lam _a _k ↔ Coincidence.PY_generic_status2 _rho _mu _sigma _lam _a _k := by
  rfl

theorem PY_q8root3_status2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : PY_q8root3_status2 _rho _mu _sigma _lam _a _k ↔ (∃ w : Vec 3, allowedCoincidence 2 w) := by
  exact (PY_q8root3_status2_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.PY_generic_status2_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate0_locus0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((((_k 0 = 0) ∧ ((_rho > 0) ∨ (_rho < 0))) ∧ ((_k 1 = 0) ∧ ((_rho > 0) ∨ (_rho < 0)))) ∧ ((_k 2 = 0) ∧ ((_rho > 0) ∨ (_rho < 0))))

theorem WL_generic_aggregate0_locus0_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate0_locus0 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_locus0 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate0_locus0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate0_locus0 _rho _mu _sigma _lam _a _k ↔ coincidenceGeometry 0 _k := by
  exact (WL_generic_aggregate0_locus0_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_locus0_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate0_allowed0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_generic_aggregate0_allowed0_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate0_allowed0 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_allowed0 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate0_allowed0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate0_allowed0 _rho _mu _sigma _lam _a _k ↔ allowedCoincidence 0 _k := by
  exact (WL_generic_aggregate0_allowed0_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_allowed0_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate0_outcome0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_generic_aggregate0_outcome0_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate0_outcome0 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_outcome0 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate0_outcome0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate0_outcome0 _rho _mu _sigma _lam _a _k ↔ (∃ w : Vec 3, allowedCoincidence 0 w) := by
  exact (WL_generic_aggregate0_outcome0_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_outcome0_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate0_status0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_generic_aggregate0_status0_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate0_status0 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_status0 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate0_status0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate0_status0 _rho _mu _sigma _lam _a _k ↔ (∃ w : Vec 3, allowedCoincidence 0 w) := by
  exact (WL_generic_aggregate0_status0_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_status0_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate0_equation0 (_rho _mu _sigma z : ℝ) (_k : Vec 3) : Prop :=
  Expr.eval (values _rho _mu _sigma z _k) WLRecordsGeneric.n15 = 0

theorem WL_generic_aggregate0_equation0_locus (_rho _mu _sigma _lam z : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate0_equation0 _rho _mu _sigma z _k ↔ WL_generic_aggregate0_locus0 _rho _mu _sigma _lam _a _k := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [WL_generic_aggregate0_equation0, WLRecordsGeneric.q3_root_coincidence_loci_cell2 _rho _mu _sigma z _k (by positivity) (by positivity),
    sub_eq_zero, Coincidence.WL_generic_pair0_roots _rho _mu _sigma _lam _k _hd,
    WL_generic_aggregate0_locus0_correct _rho _mu _sigma _lam _a _k _hd]

theorem WL_generic_aggregate0_status0_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate0_status0 _rho _mu _sigma _lam _a _k ↔ ∃ w : Vec 3, WL_generic_aggregate0_allowed0 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_generic_aggregate0_status0_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_aggregate0_allowed0_correct _rho _mu _sigma _lam _a _ _hd]

theorem WL_generic_aggregate0_outcome0_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate0_outcome0 _rho _mu _sigma _lam _a _k ↔ ∃ w : Vec 3, WL_generic_aggregate0_allowed0 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_generic_aggregate0_outcome0_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_aggregate0_allowed0_correct _rho _mu _sigma _lam _a _ _hd]

def WL_generic_aggregate0_locus1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((((((_k 2 = (-(Real.sqrt ((-(_k 1 ^ 2)) - ((_k 0 ^ 2) * _sigma))))) ∧ (((_k 0 < (-(Real.sqrt (-((_k 1 ^ 2) / _sigma))))) ∧ (_sigma < 0)) ∨ ((_sigma < 0) ∧ (_k 0 > (Real.sqrt (-((_k 1 ^ 2) / _sigma))))))) ∨ ((_k 2 = (Real.sqrt ((-(_k 1 ^ 2)) - ((_k 0 ^ 2) * _sigma)))) ∧ (((_k 0 < (-(Real.sqrt (-((_k 1 ^ 2) / _sigma))))) ∧ (_sigma < 0)) ∨ ((_sigma < 0) ∧ (_k 0 > (Real.sqrt (-((_k 1 ^ 2) / _sigma)))))))) ∨ (((_k 0 = (-(Real.sqrt (-((_k 1 ^ 2) / _sigma))))) ∧ (_sigma < 0)) ∧ ((_k 2 = 0) ∧ (_sigma < 0)))) ∨ (((_k 0 = (Real.sqrt (-((_k 1 ^ 2) / _sigma)))) ∧ (_sigma < 0)) ∧ ((_k 2 = 0) ∧ (_sigma < 0)))) ∨ ((((_k 0 = 0) ∧ (_sigma > 0)) ∧ ((_k 1 = 0) ∧ (_sigma > 0))) ∧ ((_k 2 = 0) ∧ (_sigma > 0))))

theorem WL_generic_aggregate0_locus1_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate0_locus1 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_locus1 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate0_locus1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate0_locus1 _rho _mu _sigma _lam _a _k ↔ coincidenceGeometry 1 _k := by
  exact (WL_generic_aggregate0_locus1_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_locus1_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate0_allowed1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_generic_aggregate0_allowed1_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate0_allowed1 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_allowed1 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate0_allowed1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate0_allowed1 _rho _mu _sigma _lam _a _k ↔ allowedCoincidence 1 _k := by
  exact (WL_generic_aggregate0_allowed1_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_allowed1_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate0_outcome1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_generic_aggregate0_outcome1_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate0_outcome1 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_outcome1 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate0_outcome1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate0_outcome1 _rho _mu _sigma _lam _a _k ↔ (∃ w : Vec 3, allowedCoincidence 1 w) := by
  exact (WL_generic_aggregate0_outcome1_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_outcome1_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate0_status1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_generic_aggregate0_status1_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate0_status1 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_status1 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate0_status1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate0_status1 _rho _mu _sigma _lam _a _k ↔ (∃ w : Vec 3, allowedCoincidence 1 w) := by
  exact (WL_generic_aggregate0_status1_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_status1_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate0_equation1 (_rho _mu _sigma z : ℝ) (_k : Vec 3) : Prop :=
  Expr.eval (values _rho _mu _sigma z _k) WLRecordsGeneric.n24 = 0

theorem WL_generic_aggregate0_equation1_locus (_rho _mu _sigma _lam z : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate0_equation1 _rho _mu _sigma z _k ↔ WL_generic_aggregate0_locus1 _rho _mu _sigma _lam _a _k := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [WL_generic_aggregate0_equation1, WLRecordsGeneric.q3_root_coincidence_loci_cell5 _rho _mu _sigma z _k (by positivity) (by positivity),
    sub_eq_zero, Coincidence.WL_generic_pair1_roots _rho _mu _sigma _lam _k _hd,
    WL_generic_aggregate0_locus1_correct _rho _mu _sigma _lam _a _k _hd]

theorem WL_generic_aggregate0_status1_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate0_status1 _rho _mu _sigma _lam _a _k ↔ ∃ w : Vec 3, WL_generic_aggregate0_allowed1 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_generic_aggregate0_status1_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_aggregate0_allowed1_correct _rho _mu _sigma _lam _a _ _hd]

theorem WL_generic_aggregate0_outcome1_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate0_outcome1 _rho _mu _sigma _lam _a _k ↔ ∃ w : Vec 3, WL_generic_aggregate0_allowed1 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_generic_aggregate0_outcome1_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_aggregate0_allowed1_correct _rho _mu _sigma _lam _a _ _hd]

def WL_generic_aggregate0_locus2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((_k 1 = 0) ∧ (_k 2 = 0))

theorem WL_generic_aggregate0_locus2_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate0_locus2 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_locus2 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate0_locus2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate0_locus2 _rho _mu _sigma _lam _a _k ↔ coincidenceGeometry 2 _k := by
  exact (WL_generic_aggregate0_locus2_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_locus2_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate0_allowed2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((((((True ∧ ((3 : ℝ) ≥ 1)) ∧ (_lam > 0)) ∧ (_mu > 0)) ∧ (_rho > 0)) ∧ (((0 < _sigma) ∧ (_sigma < 1)) ∨ (_sigma > 1))) ∧ ((((_k 0 < 0) ∧ (_k 1 = 0)) ∧ (_k 2 = 0)) ∨ (((_k 0 > 0) ∧ (_k 1 = 0)) ∧ (_k 2 = 0))))

theorem WL_generic_aggregate0_allowed2_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate0_allowed2 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_allowed2 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate0_allowed2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate0_allowed2 _rho _mu _sigma _lam _a _k ↔ allowedCoincidence 2 _k := by
  exact (WL_generic_aggregate0_allowed2_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_allowed2_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate0_outcome2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := True

theorem WL_generic_aggregate0_outcome2_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate0_outcome2 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_outcome2 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate0_outcome2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate0_outcome2 _rho _mu _sigma _lam _a _k ↔ (∃ w : Vec 3, allowedCoincidence 2 w) := by
  exact (WL_generic_aggregate0_outcome2_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_outcome2_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate0_status2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := True

theorem WL_generic_aggregate0_status2_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate0_status2 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_status2 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate0_status2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate0_status2 _rho _mu _sigma _lam _a _k ↔ (∃ w : Vec 3, allowedCoincidence 2 w) := by
  exact (WL_generic_aggregate0_status2_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_status2_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate0_equation2 (_rho _mu _sigma z : ℝ) (_k : Vec 3) : Prop :=
  Expr.eval (values _rho _mu _sigma z _k) WLRecordsGeneric.n32 = 0

theorem WL_generic_aggregate0_equation2_locus (_rho _mu _sigma _lam z : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate0_equation2 _rho _mu _sigma z _k ↔ WL_generic_aggregate0_locus2 _rho _mu _sigma _lam _a _k := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [WL_generic_aggregate0_equation2, WLRecordsGeneric.q3_root_coincidence_loci_cell8 _rho _mu _sigma z _k (by positivity) (by positivity),
    sub_eq_zero, Coincidence.WL_generic_pair2_roots _rho _mu _sigma _lam _k _hd,
    WL_generic_aggregate0_locus2_correct _rho _mu _sigma _lam _a _k _hd]

theorem WL_generic_aggregate0_status2_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate0_status2 _rho _mu _sigma _lam _a _k ↔ ∃ w : Vec 3, WL_generic_aggregate0_allowed2 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_generic_aggregate0_status2_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_aggregate0_allowed2_correct _rho _mu _sigma _lam _a _ _hd]

theorem WL_generic_aggregate0_outcome2_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate0_outcome2 _rho _mu _sigma _lam _a _k ↔ ∃ w : Vec 3, WL_generic_aggregate0_allowed2 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_generic_aggregate0_outcome2_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_aggregate0_allowed2_correct _rho _mu _sigma _lam _a _ _hd]

def WL_generic_aggregate0_pairs : List (ℕ × ℕ) := [(1, 2), (1, 3), (2, 3)]
theorem WL_generic_aggregate0_pairs_complete : WL_generic_aggregate0_pairs = [(1, 2), (1, 3), (2, 3)] := by decide

def WL_generic_aggregate1_locus0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((((_k 0 = 0) ∧ ((_rho > 0) ∨ (_rho < 0))) ∧ ((_k 1 = 0) ∧ ((_rho > 0) ∨ (_rho < 0)))) ∧ ((_k 2 = 0) ∧ ((_rho > 0) ∨ (_rho < 0))))

theorem WL_generic_aggregate1_locus0_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate1_locus0 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_locus0 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate1_locus0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate1_locus0 _rho _mu _sigma _lam _a _k ↔ coincidenceGeometry 0 _k := by
  exact (WL_generic_aggregate1_locus0_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_locus0_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate1_allowed0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_generic_aggregate1_allowed0_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate1_allowed0 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_allowed0 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate1_allowed0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate1_allowed0 _rho _mu _sigma _lam _a _k ↔ allowedCoincidence 0 _k := by
  exact (WL_generic_aggregate1_allowed0_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_allowed0_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate1_outcome0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_generic_aggregate1_outcome0_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate1_outcome0 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_outcome0 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate1_outcome0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate1_outcome0 _rho _mu _sigma _lam _a _k ↔ (∃ w : Vec 3, allowedCoincidence 0 w) := by
  exact (WL_generic_aggregate1_outcome0_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_outcome0_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate1_status0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_generic_aggregate1_status0_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate1_status0 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_status0 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate1_status0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate1_status0 _rho _mu _sigma _lam _a _k ↔ (∃ w : Vec 3, allowedCoincidence 0 w) := by
  exact (WL_generic_aggregate1_status0_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_status0_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate1_equation0 (_rho _mu _sigma z : ℝ) (_k : Vec 3) : Prop :=
  Expr.eval (values _rho _mu _sigma z _k) WLRecordsGeneric.n15 = 0

theorem WL_generic_aggregate1_equation0_locus (_rho _mu _sigma _lam z : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate1_equation0 _rho _mu _sigma z _k ↔ WL_generic_aggregate1_locus0 _rho _mu _sigma _lam _a _k := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [WL_generic_aggregate1_equation0, WLRecordsGeneric.root1_q8_root_coincidence_loci_cell2 _rho _mu _sigma z _k (by positivity) (by positivity),
    sub_eq_zero, Coincidence.WL_generic_pair0_roots _rho _mu _sigma _lam _k _hd,
    WL_generic_aggregate1_locus0_correct _rho _mu _sigma _lam _a _k _hd]

theorem WL_generic_aggregate1_status0_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate1_status0 _rho _mu _sigma _lam _a _k ↔ ∃ w : Vec 3, WL_generic_aggregate1_allowed0 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_generic_aggregate1_status0_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_aggregate1_allowed0_correct _rho _mu _sigma _lam _a _ _hd]

theorem WL_generic_aggregate1_outcome0_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate1_outcome0 _rho _mu _sigma _lam _a _k ↔ ∃ w : Vec 3, WL_generic_aggregate1_allowed0 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_generic_aggregate1_outcome0_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_aggregate1_allowed0_correct _rho _mu _sigma _lam _a _ _hd]

def WL_generic_aggregate1_locus1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((((((_k 2 = (-(Real.sqrt ((-(_k 1 ^ 2)) - ((_k 0 ^ 2) * _sigma))))) ∧ (((_k 0 < (-(Real.sqrt (-((_k 1 ^ 2) / _sigma))))) ∧ (_sigma < 0)) ∨ ((_sigma < 0) ∧ (_k 0 > (Real.sqrt (-((_k 1 ^ 2) / _sigma))))))) ∨ ((_k 2 = (Real.sqrt ((-(_k 1 ^ 2)) - ((_k 0 ^ 2) * _sigma)))) ∧ (((_k 0 < (-(Real.sqrt (-((_k 1 ^ 2) / _sigma))))) ∧ (_sigma < 0)) ∨ ((_sigma < 0) ∧ (_k 0 > (Real.sqrt (-((_k 1 ^ 2) / _sigma)))))))) ∨ (((_k 0 = (-(Real.sqrt (-((_k 1 ^ 2) / _sigma))))) ∧ (_sigma < 0)) ∧ ((_k 2 = 0) ∧ (_sigma < 0)))) ∨ (((_k 0 = (Real.sqrt (-((_k 1 ^ 2) / _sigma)))) ∧ (_sigma < 0)) ∧ ((_k 2 = 0) ∧ (_sigma < 0)))) ∨ ((((_k 0 = 0) ∧ (_sigma > 0)) ∧ ((_k 1 = 0) ∧ (_sigma > 0))) ∧ ((_k 2 = 0) ∧ (_sigma > 0))))

theorem WL_generic_aggregate1_locus1_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate1_locus1 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_locus1 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate1_locus1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate1_locus1 _rho _mu _sigma _lam _a _k ↔ coincidenceGeometry 1 _k := by
  exact (WL_generic_aggregate1_locus1_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_locus1_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate1_allowed1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_generic_aggregate1_allowed1_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate1_allowed1 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_allowed1 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate1_allowed1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate1_allowed1 _rho _mu _sigma _lam _a _k ↔ allowedCoincidence 1 _k := by
  exact (WL_generic_aggregate1_allowed1_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_allowed1_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate1_outcome1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_generic_aggregate1_outcome1_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate1_outcome1 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_outcome1 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate1_outcome1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate1_outcome1 _rho _mu _sigma _lam _a _k ↔ (∃ w : Vec 3, allowedCoincidence 1 w) := by
  exact (WL_generic_aggregate1_outcome1_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_outcome1_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate1_status1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_generic_aggregate1_status1_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate1_status1 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_status1 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate1_status1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate1_status1 _rho _mu _sigma _lam _a _k ↔ (∃ w : Vec 3, allowedCoincidence 1 w) := by
  exact (WL_generic_aggregate1_status1_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_status1_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate1_equation1 (_rho _mu _sigma z : ℝ) (_k : Vec 3) : Prop :=
  Expr.eval (values _rho _mu _sigma z _k) WLRecordsGeneric.n24 = 0

theorem WL_generic_aggregate1_equation1_locus (_rho _mu _sigma _lam z : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate1_equation1 _rho _mu _sigma z _k ↔ WL_generic_aggregate1_locus1 _rho _mu _sigma _lam _a _k := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [WL_generic_aggregate1_equation1, WLRecordsGeneric.root1_q8_root_coincidence_loci_cell5 _rho _mu _sigma z _k (by positivity) (by positivity),
    sub_eq_zero, Coincidence.WL_generic_pair1_roots _rho _mu _sigma _lam _k _hd,
    WL_generic_aggregate1_locus1_correct _rho _mu _sigma _lam _a _k _hd]

theorem WL_generic_aggregate1_status1_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate1_status1 _rho _mu _sigma _lam _a _k ↔ ∃ w : Vec 3, WL_generic_aggregate1_allowed1 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_generic_aggregate1_status1_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_aggregate1_allowed1_correct _rho _mu _sigma _lam _a _ _hd]

theorem WL_generic_aggregate1_outcome1_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate1_outcome1 _rho _mu _sigma _lam _a _k ↔ ∃ w : Vec 3, WL_generic_aggregate1_allowed1 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_generic_aggregate1_outcome1_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_aggregate1_allowed1_correct _rho _mu _sigma _lam _a _ _hd]

def WL_generic_aggregate1_locus2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((_k 1 = 0) ∧ (_k 2 = 0))

theorem WL_generic_aggregate1_locus2_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate1_locus2 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_locus2 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate1_locus2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate1_locus2 _rho _mu _sigma _lam _a _k ↔ coincidenceGeometry 2 _k := by
  exact (WL_generic_aggregate1_locus2_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_locus2_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate1_allowed2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((((((True ∧ ((3 : ℝ) ≥ 1)) ∧ (_lam > 0)) ∧ (_mu > 0)) ∧ (_rho > 0)) ∧ (((0 < _sigma) ∧ (_sigma < 1)) ∨ (_sigma > 1))) ∧ ((((_k 0 < 0) ∧ (_k 1 = 0)) ∧ (_k 2 = 0)) ∨ (((_k 0 > 0) ∧ (_k 1 = 0)) ∧ (_k 2 = 0))))

theorem WL_generic_aggregate1_allowed2_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate1_allowed2 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_allowed2 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate1_allowed2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate1_allowed2 _rho _mu _sigma _lam _a _k ↔ allowedCoincidence 2 _k := by
  exact (WL_generic_aggregate1_allowed2_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_allowed2_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate1_outcome2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := True

theorem WL_generic_aggregate1_outcome2_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate1_outcome2 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_outcome2 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate1_outcome2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate1_outcome2 _rho _mu _sigma _lam _a _k ↔ (∃ w : Vec 3, allowedCoincidence 2 w) := by
  exact (WL_generic_aggregate1_outcome2_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_outcome2_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate1_status2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := True

theorem WL_generic_aggregate1_status2_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate1_status2 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_status2 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate1_status2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate1_status2 _rho _mu _sigma _lam _a _k ↔ (∃ w : Vec 3, allowedCoincidence 2 w) := by
  exact (WL_generic_aggregate1_status2_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_status2_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate1_equation2 (_rho _mu _sigma z : ℝ) (_k : Vec 3) : Prop :=
  Expr.eval (values _rho _mu _sigma z _k) WLRecordsGeneric.n32 = 0

theorem WL_generic_aggregate1_equation2_locus (_rho _mu _sigma _lam z : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate1_equation2 _rho _mu _sigma z _k ↔ WL_generic_aggregate1_locus2 _rho _mu _sigma _lam _a _k := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [WL_generic_aggregate1_equation2, WLRecordsGeneric.root1_q8_root_coincidence_loci_cell8 _rho _mu _sigma z _k (by positivity) (by positivity),
    sub_eq_zero, Coincidence.WL_generic_pair2_roots _rho _mu _sigma _lam _k _hd,
    WL_generic_aggregate1_locus2_correct _rho _mu _sigma _lam _a _k _hd]

theorem WL_generic_aggregate1_status2_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate1_status2 _rho _mu _sigma _lam _a _k ↔ ∃ w : Vec 3, WL_generic_aggregate1_allowed2 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_generic_aggregate1_status2_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_aggregate1_allowed2_correct _rho _mu _sigma _lam _a _ _hd]

theorem WL_generic_aggregate1_outcome2_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate1_outcome2 _rho _mu _sigma _lam _a _k ↔ ∃ w : Vec 3, WL_generic_aggregate1_allowed2 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_generic_aggregate1_outcome2_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_aggregate1_allowed2_correct _rho _mu _sigma _lam _a _ _hd]

def WL_generic_aggregate1_pairs : List (ℕ × ℕ) := [(1, 2), (1, 3), (2, 3)]
theorem WL_generic_aggregate1_pairs_complete : WL_generic_aggregate1_pairs = [(1, 2), (1, 3), (2, 3)] := by decide

def WL_generic_aggregate2_locus0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((((_k 0 = 0) ∧ ((_rho > 0) ∨ (_rho < 0))) ∧ ((_k 1 = 0) ∧ ((_rho > 0) ∨ (_rho < 0)))) ∧ ((_k 2 = 0) ∧ ((_rho > 0) ∨ (_rho < 0))))

theorem WL_generic_aggregate2_locus0_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate2_locus0 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_locus0 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate2_locus0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate2_locus0 _rho _mu _sigma _lam _a _k ↔ coincidenceGeometry 0 _k := by
  exact (WL_generic_aggregate2_locus0_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_locus0_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate2_allowed0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_generic_aggregate2_allowed0_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate2_allowed0 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_allowed0 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate2_allowed0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate2_allowed0 _rho _mu _sigma _lam _a _k ↔ allowedCoincidence 0 _k := by
  exact (WL_generic_aggregate2_allowed0_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_allowed0_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate2_outcome0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_generic_aggregate2_outcome0_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate2_outcome0 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_outcome0 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate2_outcome0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate2_outcome0 _rho _mu _sigma _lam _a _k ↔ (∃ w : Vec 3, allowedCoincidence 0 w) := by
  exact (WL_generic_aggregate2_outcome0_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_outcome0_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate2_status0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_generic_aggregate2_status0_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate2_status0 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_status0 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate2_status0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate2_status0 _rho _mu _sigma _lam _a _k ↔ (∃ w : Vec 3, allowedCoincidence 0 w) := by
  exact (WL_generic_aggregate2_status0_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_status0_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate2_equation0 (_rho _mu _sigma z : ℝ) (_k : Vec 3) : Prop :=
  Expr.eval (values _rho _mu _sigma z _k) WLRecordsGeneric.n15 = 0

theorem WL_generic_aggregate2_equation0_locus (_rho _mu _sigma _lam z : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate2_equation0 _rho _mu _sigma z _k ↔ WL_generic_aggregate2_locus0 _rho _mu _sigma _lam _a _k := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [WL_generic_aggregate2_equation0, WLRecordsGeneric.root2_q8_root_coincidence_loci_cell2 _rho _mu _sigma z _k (by positivity) (by positivity),
    sub_eq_zero, Coincidence.WL_generic_pair0_roots _rho _mu _sigma _lam _k _hd,
    WL_generic_aggregate2_locus0_correct _rho _mu _sigma _lam _a _k _hd]

theorem WL_generic_aggregate2_status0_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate2_status0 _rho _mu _sigma _lam _a _k ↔ ∃ w : Vec 3, WL_generic_aggregate2_allowed0 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_generic_aggregate2_status0_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_aggregate2_allowed0_correct _rho _mu _sigma _lam _a _ _hd]

theorem WL_generic_aggregate2_outcome0_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate2_outcome0 _rho _mu _sigma _lam _a _k ↔ ∃ w : Vec 3, WL_generic_aggregate2_allowed0 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_generic_aggregate2_outcome0_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_aggregate2_allowed0_correct _rho _mu _sigma _lam _a _ _hd]

def WL_generic_aggregate2_locus1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((((((_k 2 = (-(Real.sqrt ((-(_k 1 ^ 2)) - ((_k 0 ^ 2) * _sigma))))) ∧ (((_k 0 < (-(Real.sqrt (-((_k 1 ^ 2) / _sigma))))) ∧ (_sigma < 0)) ∨ ((_sigma < 0) ∧ (_k 0 > (Real.sqrt (-((_k 1 ^ 2) / _sigma))))))) ∨ ((_k 2 = (Real.sqrt ((-(_k 1 ^ 2)) - ((_k 0 ^ 2) * _sigma)))) ∧ (((_k 0 < (-(Real.sqrt (-((_k 1 ^ 2) / _sigma))))) ∧ (_sigma < 0)) ∨ ((_sigma < 0) ∧ (_k 0 > (Real.sqrt (-((_k 1 ^ 2) / _sigma)))))))) ∨ (((_k 0 = (-(Real.sqrt (-((_k 1 ^ 2) / _sigma))))) ∧ (_sigma < 0)) ∧ ((_k 2 = 0) ∧ (_sigma < 0)))) ∨ (((_k 0 = (Real.sqrt (-((_k 1 ^ 2) / _sigma)))) ∧ (_sigma < 0)) ∧ ((_k 2 = 0) ∧ (_sigma < 0)))) ∨ ((((_k 0 = 0) ∧ (_sigma > 0)) ∧ ((_k 1 = 0) ∧ (_sigma > 0))) ∧ ((_k 2 = 0) ∧ (_sigma > 0))))

theorem WL_generic_aggregate2_locus1_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate2_locus1 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_locus1 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate2_locus1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate2_locus1 _rho _mu _sigma _lam _a _k ↔ coincidenceGeometry 1 _k := by
  exact (WL_generic_aggregate2_locus1_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_locus1_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate2_allowed1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_generic_aggregate2_allowed1_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate2_allowed1 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_allowed1 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate2_allowed1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate2_allowed1 _rho _mu _sigma _lam _a _k ↔ allowedCoincidence 1 _k := by
  exact (WL_generic_aggregate2_allowed1_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_allowed1_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate2_outcome1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_generic_aggregate2_outcome1_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate2_outcome1 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_outcome1 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate2_outcome1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate2_outcome1 _rho _mu _sigma _lam _a _k ↔ (∃ w : Vec 3, allowedCoincidence 1 w) := by
  exact (WL_generic_aggregate2_outcome1_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_outcome1_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate2_status1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_generic_aggregate2_status1_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate2_status1 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_status1 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate2_status1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate2_status1 _rho _mu _sigma _lam _a _k ↔ (∃ w : Vec 3, allowedCoincidence 1 w) := by
  exact (WL_generic_aggregate2_status1_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_status1_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate2_equation1 (_rho _mu _sigma z : ℝ) (_k : Vec 3) : Prop :=
  Expr.eval (values _rho _mu _sigma z _k) WLRecordsGeneric.n24 = 0

theorem WL_generic_aggregate2_equation1_locus (_rho _mu _sigma _lam z : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate2_equation1 _rho _mu _sigma z _k ↔ WL_generic_aggregate2_locus1 _rho _mu _sigma _lam _a _k := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [WL_generic_aggregate2_equation1, WLRecordsGeneric.root2_q8_root_coincidence_loci_cell5 _rho _mu _sigma z _k (by positivity) (by positivity),
    sub_eq_zero, Coincidence.WL_generic_pair1_roots _rho _mu _sigma _lam _k _hd,
    WL_generic_aggregate2_locus1_correct _rho _mu _sigma _lam _a _k _hd]

theorem WL_generic_aggregate2_status1_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate2_status1 _rho _mu _sigma _lam _a _k ↔ ∃ w : Vec 3, WL_generic_aggregate2_allowed1 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_generic_aggregate2_status1_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_aggregate2_allowed1_correct _rho _mu _sigma _lam _a _ _hd]

theorem WL_generic_aggregate2_outcome1_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate2_outcome1 _rho _mu _sigma _lam _a _k ↔ ∃ w : Vec 3, WL_generic_aggregate2_allowed1 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_generic_aggregate2_outcome1_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_aggregate2_allowed1_correct _rho _mu _sigma _lam _a _ _hd]

def WL_generic_aggregate2_locus2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((_k 1 = 0) ∧ (_k 2 = 0))

theorem WL_generic_aggregate2_locus2_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate2_locus2 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_locus2 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate2_locus2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate2_locus2 _rho _mu _sigma _lam _a _k ↔ coincidenceGeometry 2 _k := by
  exact (WL_generic_aggregate2_locus2_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_locus2_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate2_allowed2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((((((True ∧ ((3 : ℝ) ≥ 1)) ∧ (_lam > 0)) ∧ (_mu > 0)) ∧ (_rho > 0)) ∧ (((0 < _sigma) ∧ (_sigma < 1)) ∨ (_sigma > 1))) ∧ ((((_k 0 < 0) ∧ (_k 1 = 0)) ∧ (_k 2 = 0)) ∨ (((_k 0 > 0) ∧ (_k 1 = 0)) ∧ (_k 2 = 0))))

theorem WL_generic_aggregate2_allowed2_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate2_allowed2 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_allowed2 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate2_allowed2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate2_allowed2 _rho _mu _sigma _lam _a _k ↔ allowedCoincidence 2 _k := by
  exact (WL_generic_aggregate2_allowed2_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_allowed2_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate2_outcome2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := True

theorem WL_generic_aggregate2_outcome2_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate2_outcome2 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_outcome2 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate2_outcome2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate2_outcome2 _rho _mu _sigma _lam _a _k ↔ (∃ w : Vec 3, allowedCoincidence 2 w) := by
  exact (WL_generic_aggregate2_outcome2_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_outcome2_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate2_status2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := True

theorem WL_generic_aggregate2_status2_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate2_status2 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_status2 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate2_status2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate2_status2 _rho _mu _sigma _lam _a _k ↔ (∃ w : Vec 3, allowedCoincidence 2 w) := by
  exact (WL_generic_aggregate2_status2_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_status2_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate2_equation2 (_rho _mu _sigma z : ℝ) (_k : Vec 3) : Prop :=
  Expr.eval (values _rho _mu _sigma z _k) WLRecordsGeneric.n32 = 0

theorem WL_generic_aggregate2_equation2_locus (_rho _mu _sigma _lam z : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate2_equation2 _rho _mu _sigma z _k ↔ WL_generic_aggregate2_locus2 _rho _mu _sigma _lam _a _k := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [WL_generic_aggregate2_equation2, WLRecordsGeneric.root2_q8_root_coincidence_loci_cell8 _rho _mu _sigma z _k (by positivity) (by positivity),
    sub_eq_zero, Coincidence.WL_generic_pair2_roots _rho _mu _sigma _lam _k _hd,
    WL_generic_aggregate2_locus2_correct _rho _mu _sigma _lam _a _k _hd]

theorem WL_generic_aggregate2_status2_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate2_status2 _rho _mu _sigma _lam _a _k ↔ ∃ w : Vec 3, WL_generic_aggregate2_allowed2 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_generic_aggregate2_status2_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_aggregate2_allowed2_correct _rho _mu _sigma _lam _a _ _hd]

theorem WL_generic_aggregate2_outcome2_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate2_outcome2 _rho _mu _sigma _lam _a _k ↔ ∃ w : Vec 3, WL_generic_aggregate2_allowed2 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_generic_aggregate2_outcome2_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_aggregate2_allowed2_correct _rho _mu _sigma _lam _a _ _hd]

def WL_generic_aggregate2_pairs : List (ℕ × ℕ) := [(1, 2), (1, 3), (2, 3)]
theorem WL_generic_aggregate2_pairs_complete : WL_generic_aggregate2_pairs = [(1, 2), (1, 3), (2, 3)] := by decide

def WL_generic_aggregate3_locus0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((((_k 0 = 0) ∧ ((_rho > 0) ∨ (_rho < 0))) ∧ ((_k 1 = 0) ∧ ((_rho > 0) ∨ (_rho < 0)))) ∧ ((_k 2 = 0) ∧ ((_rho > 0) ∨ (_rho < 0))))

theorem WL_generic_aggregate3_locus0_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate3_locus0 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_locus0 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate3_locus0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate3_locus0 _rho _mu _sigma _lam _a _k ↔ coincidenceGeometry 0 _k := by
  exact (WL_generic_aggregate3_locus0_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_locus0_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate3_allowed0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_generic_aggregate3_allowed0_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate3_allowed0 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_allowed0 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate3_allowed0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate3_allowed0 _rho _mu _sigma _lam _a _k ↔ allowedCoincidence 0 _k := by
  exact (WL_generic_aggregate3_allowed0_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_allowed0_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate3_outcome0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_generic_aggregate3_outcome0_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate3_outcome0 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_outcome0 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate3_outcome0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate3_outcome0 _rho _mu _sigma _lam _a _k ↔ (∃ w : Vec 3, allowedCoincidence 0 w) := by
  exact (WL_generic_aggregate3_outcome0_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_outcome0_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate3_status0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_generic_aggregate3_status0_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate3_status0 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_status0 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate3_status0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate3_status0 _rho _mu _sigma _lam _a _k ↔ (∃ w : Vec 3, allowedCoincidence 0 w) := by
  exact (WL_generic_aggregate3_status0_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_status0_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate3_equation0 (_rho _mu _sigma z : ℝ) (_k : Vec 3) : Prop :=
  Expr.eval (values _rho _mu _sigma z _k) WLRecordsGeneric.n15 = 0

theorem WL_generic_aggregate3_equation0_locus (_rho _mu _sigma _lam z : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate3_equation0 _rho _mu _sigma z _k ↔ WL_generic_aggregate3_locus0 _rho _mu _sigma _lam _a _k := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [WL_generic_aggregate3_equation0, WLRecordsGeneric.root3_q8_root_coincidence_loci_cell2 _rho _mu _sigma z _k (by positivity) (by positivity),
    sub_eq_zero, Coincidence.WL_generic_pair0_roots _rho _mu _sigma _lam _k _hd,
    WL_generic_aggregate3_locus0_correct _rho _mu _sigma _lam _a _k _hd]

theorem WL_generic_aggregate3_status0_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate3_status0 _rho _mu _sigma _lam _a _k ↔ ∃ w : Vec 3, WL_generic_aggregate3_allowed0 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_generic_aggregate3_status0_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_aggregate3_allowed0_correct _rho _mu _sigma _lam _a _ _hd]

theorem WL_generic_aggregate3_outcome0_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate3_outcome0 _rho _mu _sigma _lam _a _k ↔ ∃ w : Vec 3, WL_generic_aggregate3_allowed0 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_generic_aggregate3_outcome0_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_aggregate3_allowed0_correct _rho _mu _sigma _lam _a _ _hd]

def WL_generic_aggregate3_locus1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((((((_k 2 = (-(Real.sqrt ((-(_k 1 ^ 2)) - ((_k 0 ^ 2) * _sigma))))) ∧ (((_k 0 < (-(Real.sqrt (-((_k 1 ^ 2) / _sigma))))) ∧ (_sigma < 0)) ∨ ((_sigma < 0) ∧ (_k 0 > (Real.sqrt (-((_k 1 ^ 2) / _sigma))))))) ∨ ((_k 2 = (Real.sqrt ((-(_k 1 ^ 2)) - ((_k 0 ^ 2) * _sigma)))) ∧ (((_k 0 < (-(Real.sqrt (-((_k 1 ^ 2) / _sigma))))) ∧ (_sigma < 0)) ∨ ((_sigma < 0) ∧ (_k 0 > (Real.sqrt (-((_k 1 ^ 2) / _sigma)))))))) ∨ (((_k 0 = (-(Real.sqrt (-((_k 1 ^ 2) / _sigma))))) ∧ (_sigma < 0)) ∧ ((_k 2 = 0) ∧ (_sigma < 0)))) ∨ (((_k 0 = (Real.sqrt (-((_k 1 ^ 2) / _sigma)))) ∧ (_sigma < 0)) ∧ ((_k 2 = 0) ∧ (_sigma < 0)))) ∨ ((((_k 0 = 0) ∧ (_sigma > 0)) ∧ ((_k 1 = 0) ∧ (_sigma > 0))) ∧ ((_k 2 = 0) ∧ (_sigma > 0))))

theorem WL_generic_aggregate3_locus1_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate3_locus1 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_locus1 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate3_locus1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate3_locus1 _rho _mu _sigma _lam _a _k ↔ coincidenceGeometry 1 _k := by
  exact (WL_generic_aggregate3_locus1_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_locus1_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate3_allowed1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_generic_aggregate3_allowed1_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate3_allowed1 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_allowed1 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate3_allowed1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate3_allowed1 _rho _mu _sigma _lam _a _k ↔ allowedCoincidence 1 _k := by
  exact (WL_generic_aggregate3_allowed1_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_allowed1_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate3_outcome1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_generic_aggregate3_outcome1_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate3_outcome1 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_outcome1 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate3_outcome1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate3_outcome1 _rho _mu _sigma _lam _a _k ↔ (∃ w : Vec 3, allowedCoincidence 1 w) := by
  exact (WL_generic_aggregate3_outcome1_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_outcome1_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate3_status1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_generic_aggregate3_status1_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate3_status1 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_status1 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate3_status1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate3_status1 _rho _mu _sigma _lam _a _k ↔ (∃ w : Vec 3, allowedCoincidence 1 w) := by
  exact (WL_generic_aggregate3_status1_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_status1_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate3_equation1 (_rho _mu _sigma z : ℝ) (_k : Vec 3) : Prop :=
  Expr.eval (values _rho _mu _sigma z _k) WLRecordsGeneric.n24 = 0

theorem WL_generic_aggregate3_equation1_locus (_rho _mu _sigma _lam z : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate3_equation1 _rho _mu _sigma z _k ↔ WL_generic_aggregate3_locus1 _rho _mu _sigma _lam _a _k := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [WL_generic_aggregate3_equation1, WLRecordsGeneric.root3_q8_root_coincidence_loci_cell5 _rho _mu _sigma z _k (by positivity) (by positivity),
    sub_eq_zero, Coincidence.WL_generic_pair1_roots _rho _mu _sigma _lam _k _hd,
    WL_generic_aggregate3_locus1_correct _rho _mu _sigma _lam _a _k _hd]

theorem WL_generic_aggregate3_status1_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate3_status1 _rho _mu _sigma _lam _a _k ↔ ∃ w : Vec 3, WL_generic_aggregate3_allowed1 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_generic_aggregate3_status1_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_aggregate3_allowed1_correct _rho _mu _sigma _lam _a _ _hd]

theorem WL_generic_aggregate3_outcome1_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate3_outcome1 _rho _mu _sigma _lam _a _k ↔ ∃ w : Vec 3, WL_generic_aggregate3_allowed1 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_generic_aggregate3_outcome1_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_aggregate3_allowed1_correct _rho _mu _sigma _lam _a _ _hd]

def WL_generic_aggregate3_locus2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((_k 1 = 0) ∧ (_k 2 = 0))

theorem WL_generic_aggregate3_locus2_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate3_locus2 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_locus2 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate3_locus2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate3_locus2 _rho _mu _sigma _lam _a _k ↔ coincidenceGeometry 2 _k := by
  exact (WL_generic_aggregate3_locus2_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_locus2_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate3_allowed2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((((((True ∧ ((3 : ℝ) ≥ 1)) ∧ (_lam > 0)) ∧ (_mu > 0)) ∧ (_rho > 0)) ∧ (((0 < _sigma) ∧ (_sigma < 1)) ∨ (_sigma > 1))) ∧ ((((_k 0 < 0) ∧ (_k 1 = 0)) ∧ (_k 2 = 0)) ∨ (((_k 0 > 0) ∧ (_k 1 = 0)) ∧ (_k 2 = 0))))

theorem WL_generic_aggregate3_allowed2_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate3_allowed2 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_allowed2 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate3_allowed2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate3_allowed2 _rho _mu _sigma _lam _a _k ↔ allowedCoincidence 2 _k := by
  exact (WL_generic_aggregate3_allowed2_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_allowed2_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate3_outcome2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := True

theorem WL_generic_aggregate3_outcome2_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate3_outcome2 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_outcome2 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate3_outcome2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate3_outcome2 _rho _mu _sigma _lam _a _k ↔ (∃ w : Vec 3, allowedCoincidence 2 w) := by
  exact (WL_generic_aggregate3_outcome2_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_outcome2_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate3_status2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := True

theorem WL_generic_aggregate3_status2_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_generic_aggregate3_status2 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_generic_status2 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_generic_aggregate3_status2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_generic_aggregate3_status2 _rho _mu _sigma _lam _a _k ↔ (∃ w : Vec 3, allowedCoincidence 2 w) := by
  exact (WL_generic_aggregate3_status2_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_generic_status2_correct _rho _mu _sigma _lam _a _k _hd)

def WL_generic_aggregate3_equation2 (_rho _mu _sigma z : ℝ) (_k : Vec 3) : Prop :=
  Expr.eval (values _rho _mu _sigma z _k) WLRecordsGeneric.n32 = 0

theorem WL_generic_aggregate3_equation2_locus (_rho _mu _sigma _lam z : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate3_equation2 _rho _mu _sigma z _k ↔ WL_generic_aggregate3_locus2 _rho _mu _sigma _lam _a _k := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [WL_generic_aggregate3_equation2, WLRecordsGeneric.root3_q8_root_coincidence_loci_cell8 _rho _mu _sigma z _k (by positivity) (by positivity),
    sub_eq_zero, Coincidence.WL_generic_pair2_roots _rho _mu _sigma _lam _k _hd,
    WL_generic_aggregate3_locus2_correct _rho _mu _sigma _lam _a _k _hd]

theorem WL_generic_aggregate3_status2_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate3_status2 _rho _mu _sigma _lam _a _k ↔ ∃ w : Vec 3, WL_generic_aggregate3_allowed2 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_generic_aggregate3_status2_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_aggregate3_allowed2_correct _rho _mu _sigma _lam _a _ _hd]

theorem WL_generic_aggregate3_outcome2_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_aggregate3_outcome2 _rho _mu _sigma _lam _a _k ↔ ∃ w : Vec 3, WL_generic_aggregate3_allowed2 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_generic_aggregate3_outcome2_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_aggregate3_allowed2_correct _rho _mu _sigma _lam _a _ _hd]

def WL_generic_aggregate3_pairs : List (ℕ × ℕ) := [(1, 2), (1, 3), (2, 3)]
theorem WL_generic_aggregate3_pairs_complete : WL_generic_aggregate3_pairs = [(1, 2), (1, 3), (2, 3)] := by decide

def WL_parallel_aggregate0_locus0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((_mu = 0) ∧ ((_rho > 0) ∨ (_rho < 0)))

theorem WL_parallel_aggregate0_locus0_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_parallel_aggregate0_locus0 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_parallel_locus0 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_parallel_aggregate0_locus0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_parallel_aggregate0_locus0 _rho _mu _sigma _lam _a _k ↔ False := by
  exact (WL_parallel_aggregate0_locus0_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_parallel_locus0_correct _rho _mu _sigma _lam _a _k _hd)

def WL_parallel_aggregate0_allowed0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_parallel_aggregate0_allowed0_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_parallel_aggregate0_allowed0 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_parallel_allowed0 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_parallel_aggregate0_allowed0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_parallel_aggregate0_allowed0 _rho _mu _sigma _lam _a _k ↔ False := by
  exact (WL_parallel_aggregate0_allowed0_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_parallel_allowed0_correct _rho _mu _sigma _lam _a _k _hd)

def WL_parallel_aggregate0_outcome0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_parallel_aggregate0_outcome0_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_parallel_aggregate0_outcome0 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_parallel_outcome0 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_parallel_aggregate0_outcome0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_parallel_aggregate0_outcome0 _rho _mu _sigma _lam _a _k ↔ False := by
  exact (WL_parallel_aggregate0_outcome0_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_parallel_outcome0_correct _rho _mu _sigma _lam _a _k _hd)

def WL_parallel_aggregate0_status0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_parallel_aggregate0_status0_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_parallel_aggregate0_status0 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_parallel_status0 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_parallel_aggregate0_status0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_parallel_aggregate0_status0 _rho _mu _sigma _lam _a _k ↔ False := by
  exact (WL_parallel_aggregate0_status0_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_parallel_status0_correct _rho _mu _sigma _lam _a _k _hd)

def WL_parallel_aggregate0_equation0 (_rho _mu _sigma z : ℝ) (_k : Vec 3) : Prop :=
  Expr.eval (values _rho _mu _sigma z _k) WLRecordsParallel.n8 = 0

theorem WL_parallel_aggregate0_equation0_locus (_rho _mu _sigma _lam z : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_parallel_aggregate0_equation0 _rho _mu _sigma z _k ↔ WL_parallel_aggregate0_locus0 _rho _mu _sigma _lam _a _k := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [WL_parallel_aggregate0_equation0, WLRecordsParallel.stratum1_q3_stratum_root_coincidence_records_cell2 _rho _mu _sigma z _k (by positivity) (by positivity),
    sub_eq_zero, Coincidence.WL_parallel_pair0_roots _rho _mu _sigma _lam _k _hd,
    WL_parallel_aggregate0_locus0_correct _rho _mu _sigma _lam _a _k _hd]

theorem WL_parallel_aggregate0_status0_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_parallel_aggregate0_status0 _rho _mu _sigma _lam _a _k ↔ ∃ w : Vec 3, WL_parallel_aggregate0_allowed0 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_parallel_aggregate0_status0_correct _rho _mu _sigma _lam _a _k _hd,
    WL_parallel_aggregate0_allowed0_correct _rho _mu _sigma _lam _a _ _hd]
  simp

theorem WL_parallel_aggregate0_outcome0_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_parallel_aggregate0_outcome0 _rho _mu _sigma _lam _a _k ↔ ∃ w : Vec 3, WL_parallel_aggregate0_allowed0 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_parallel_aggregate0_outcome0_correct _rho _mu _sigma _lam _a _k _hd,
    WL_parallel_aggregate0_allowed0_correct _rho _mu _sigma _lam _a _ _hd]
  simp

def WL_parallel_aggregate0_pairs : List (ℕ × ℕ) := [(1, 2)]
theorem WL_parallel_aggregate0_pairs_complete : WL_parallel_aggregate0_pairs = [(1, 2)] := by decide

def WL_perpendicular_aggregate0_locus0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((_mu = 0) ∧ ((_rho > 0) ∨ (_rho < 0)))

theorem WL_perpendicular_aggregate0_locus0_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_perpendicular_aggregate0_locus0 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_perpendicular_locus0 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_perpendicular_aggregate0_locus0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_perpendicular_aggregate0_locus0 _rho _mu _sigma _lam _a _k ↔ False := by
  exact (WL_perpendicular_aggregate0_locus0_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_perpendicular_locus0_correct _rho _mu _sigma _lam _a _k _hd)

def WL_perpendicular_aggregate0_allowed0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_perpendicular_aggregate0_allowed0_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_perpendicular_aggregate0_allowed0 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_perpendicular_allowed0 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_perpendicular_aggregate0_allowed0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_perpendicular_aggregate0_allowed0 _rho _mu _sigma _lam _a _k ↔ False := by
  exact (WL_perpendicular_aggregate0_allowed0_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_perpendicular_allowed0_correct _rho _mu _sigma _lam _a _k _hd)

def WL_perpendicular_aggregate0_outcome0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_perpendicular_aggregate0_outcome0_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_perpendicular_aggregate0_outcome0 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_perpendicular_outcome0 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_perpendicular_aggregate0_outcome0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_perpendicular_aggregate0_outcome0 _rho _mu _sigma _lam _a _k ↔ False := by
  exact (WL_perpendicular_aggregate0_outcome0_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_perpendicular_outcome0_correct _rho _mu _sigma _lam _a _k _hd)

def WL_perpendicular_aggregate0_status0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_perpendicular_aggregate0_status0_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_perpendicular_aggregate0_status0 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_perpendicular_status0 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_perpendicular_aggregate0_status0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_perpendicular_aggregate0_status0 _rho _mu _sigma _lam _a _k ↔ False := by
  exact (WL_perpendicular_aggregate0_status0_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_perpendicular_status0_correct _rho _mu _sigma _lam _a _k _hd)

def WL_perpendicular_aggregate0_equation0 (_rho _mu _sigma z : ℝ) (_k : Vec 3) : Prop :=
  Expr.eval (values _rho _mu _sigma z _k) WLRecordsPerpendicular.n10 = 0

theorem WL_perpendicular_aggregate0_equation0_locus (_rho _mu _sigma _lam z : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_perpendicular_aggregate0_equation0 _rho _mu _sigma z _k ↔ WL_perpendicular_aggregate0_locus0 _rho _mu _sigma _lam _a _k := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [WL_perpendicular_aggregate0_equation0, WLRecordsPerpendicular.stratum2_q3_stratum_root_coincidence_records_cell2 _rho _mu _sigma z _k (by positivity) (by positivity),
    sub_eq_zero, Coincidence.WL_perpendicular_pair0_roots _rho _mu _sigma _lam _k _hd,
    WL_perpendicular_aggregate0_locus0_correct _rho _mu _sigma _lam _a _k _hd]

theorem WL_perpendicular_aggregate0_status0_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_perpendicular_aggregate0_status0 _rho _mu _sigma _lam _a _k ↔ ∃ w : Vec 3, WL_perpendicular_aggregate0_allowed0 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_perpendicular_aggregate0_status0_correct _rho _mu _sigma _lam _a _k _hd,
    WL_perpendicular_aggregate0_allowed0_correct _rho _mu _sigma _lam _a _ _hd]
  simp

theorem WL_perpendicular_aggregate0_outcome0_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_perpendicular_aggregate0_outcome0 _rho _mu _sigma _lam _a _k ↔ ∃ w : Vec 3, WL_perpendicular_aggregate0_allowed0 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_perpendicular_aggregate0_outcome0_correct _rho _mu _sigma _lam _a _k _hd,
    WL_perpendicular_aggregate0_allowed0_correct _rho _mu _sigma _lam _a _ _hd]
  simp

def WL_perpendicular_aggregate0_locus1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((_mu = 0) ∧ (((((_rho > 0) ∧ (_sigma > 0)) ∨ ((_rho > 0) ∧ (_sigma < 0))) ∨ ((_rho < 0) ∧ (_sigma > 0))) ∨ ((_rho < 0) ∧ (_sigma < 0))))

theorem WL_perpendicular_aggregate0_locus1_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_perpendicular_aggregate0_locus1 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_perpendicular_locus1 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_perpendicular_aggregate0_locus1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_perpendicular_aggregate0_locus1 _rho _mu _sigma _lam _a _k ↔ False := by
  exact (WL_perpendicular_aggregate0_locus1_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_perpendicular_locus1_correct _rho _mu _sigma _lam _a _k _hd)

def WL_perpendicular_aggregate0_allowed1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_perpendicular_aggregate0_allowed1_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_perpendicular_aggregate0_allowed1 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_perpendicular_allowed1 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_perpendicular_aggregate0_allowed1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_perpendicular_aggregate0_allowed1 _rho _mu _sigma _lam _a _k ↔ False := by
  exact (WL_perpendicular_aggregate0_allowed1_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_perpendicular_allowed1_correct _rho _mu _sigma _lam _a _k _hd)

def WL_perpendicular_aggregate0_outcome1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_perpendicular_aggregate0_outcome1_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_perpendicular_aggregate0_outcome1 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_perpendicular_outcome1 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_perpendicular_aggregate0_outcome1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_perpendicular_aggregate0_outcome1 _rho _mu _sigma _lam _a _k ↔ False := by
  exact (WL_perpendicular_aggregate0_outcome1_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_perpendicular_outcome1_correct _rho _mu _sigma _lam _a _k _hd)

def WL_perpendicular_aggregate0_status1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_perpendicular_aggregate0_status1_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_perpendicular_aggregate0_status1 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_perpendicular_status1 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_perpendicular_aggregate0_status1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_perpendicular_aggregate0_status1 _rho _mu _sigma _lam _a _k ↔ False := by
  exact (WL_perpendicular_aggregate0_status1_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_perpendicular_status1_correct _rho _mu _sigma _lam _a _k _hd)

def WL_perpendicular_aggregate0_equation1 (_rho _mu _sigma z : ℝ) (_k : Vec 3) : Prop :=
  Expr.eval (values _rho _mu _sigma z _k) WLRecordsPerpendicular.n14 = 0

theorem WL_perpendicular_aggregate0_equation1_locus (_rho _mu _sigma _lam z : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_perpendicular_aggregate0_equation1 _rho _mu _sigma z _k ↔ WL_perpendicular_aggregate0_locus1 _rho _mu _sigma _lam _a _k := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [WL_perpendicular_aggregate0_equation1, WLRecordsPerpendicular.stratum2_q3_stratum_root_coincidence_records_cell5 _rho _mu _sigma z _k (by positivity) (by positivity),
    sub_eq_zero, Coincidence.WL_perpendicular_pair1_roots _rho _mu _sigma _lam _k _hd,
    WL_perpendicular_aggregate0_locus1_correct _rho _mu _sigma _lam _a _k _hd]

theorem WL_perpendicular_aggregate0_status1_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_perpendicular_aggregate0_status1 _rho _mu _sigma _lam _a _k ↔ ∃ w : Vec 3, WL_perpendicular_aggregate0_allowed1 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_perpendicular_aggregate0_status1_correct _rho _mu _sigma _lam _a _k _hd,
    WL_perpendicular_aggregate0_allowed1_correct _rho _mu _sigma _lam _a _ _hd]
  simp

theorem WL_perpendicular_aggregate0_outcome1_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_perpendicular_aggregate0_outcome1 _rho _mu _sigma _lam _a _k ↔ ∃ w : Vec 3, WL_perpendicular_aggregate0_allowed1 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_perpendicular_aggregate0_outcome1_correct _rho _mu _sigma _lam _a _k _hd,
    WL_perpendicular_aggregate0_allowed1_correct _rho _mu _sigma _lam _a _ _hd]
  simp

def WL_perpendicular_aggregate0_locus2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := (((_mu = 0) ∧ ((((((((0 < _sigma) ∧ (_sigma < 1)) ∧ (_rho > 0)) ∨ (((0 < _sigma) ∧ (_sigma < 1)) ∧ (_rho < 0))) ∨ ((_sigma > 1) ∧ (_rho > 0))) ∨ ((_sigma > 1) ∧ (_rho < 0))) ∨ ((_sigma < 0) ∧ (_rho > 0))) ∨ ((_sigma < 0) ∧ (_rho < 0)))) ∨ ((_sigma = 1) ∧ ((_rho > 0) ∨ (_rho < 0))))

theorem WL_perpendicular_aggregate0_locus2_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_perpendicular_aggregate0_locus2 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_perpendicular_locus2 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_perpendicular_aggregate0_locus2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_perpendicular_aggregate0_locus2 _rho _mu _sigma _lam _a _k ↔ False := by
  exact (WL_perpendicular_aggregate0_locus2_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_perpendicular_locus2_correct _rho _mu _sigma _lam _a _k _hd)

def WL_perpendicular_aggregate0_allowed2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_perpendicular_aggregate0_allowed2_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_perpendicular_aggregate0_allowed2 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_perpendicular_allowed2 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_perpendicular_aggregate0_allowed2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_perpendicular_aggregate0_allowed2 _rho _mu _sigma _lam _a _k ↔ False := by
  exact (WL_perpendicular_aggregate0_allowed2_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_perpendicular_allowed2_correct _rho _mu _sigma _lam _a _k _hd)

def WL_perpendicular_aggregate0_outcome2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_perpendicular_aggregate0_outcome2_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_perpendicular_aggregate0_outcome2 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_perpendicular_outcome2 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_perpendicular_aggregate0_outcome2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_perpendicular_aggregate0_outcome2 _rho _mu _sigma _lam _a _k ↔ False := by
  exact (WL_perpendicular_aggregate0_outcome2_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_perpendicular_outcome2_correct _rho _mu _sigma _lam _a _k _hd)

def WL_perpendicular_aggregate0_status2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_perpendicular_aggregate0_status2_primary (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : WL_perpendicular_aggregate0_status2 _rho _mu _sigma _lam _a _k ↔ Coincidence.WL_perpendicular_status2 _rho _mu _sigma _lam _a _k := by
  rfl

theorem WL_perpendicular_aggregate0_status2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) : WL_perpendicular_aggregate0_status2 _rho _mu _sigma _lam _a _k ↔ False := by
  exact (WL_perpendicular_aggregate0_status2_primary _rho _mu _sigma _lam _a _k).trans
    (Coincidence.WL_perpendicular_status2_correct _rho _mu _sigma _lam _a _k _hd)

def WL_perpendicular_aggregate0_equation2 (_rho _mu _sigma z : ℝ) (_k : Vec 3) : Prop :=
  Expr.eval (values _rho _mu _sigma z _k) WLRecordsPerpendicular.n19 = 0

theorem WL_perpendicular_aggregate0_equation2_locus (_rho _mu _sigma _lam z : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_perpendicular_aggregate0_equation2 _rho _mu _sigma z _k ↔ WL_perpendicular_aggregate0_locus2 _rho _mu _sigma _lam _a _k := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [WL_perpendicular_aggregate0_equation2, WLRecordsPerpendicular.stratum2_q3_stratum_root_coincidence_records_cell8 _rho _mu _sigma z _k (by positivity) (by positivity),
    sub_eq_zero, Coincidence.WL_perpendicular_pair2_roots _rho _mu _sigma _lam _k _hd,
    WL_perpendicular_aggregate0_locus2_correct _rho _mu _sigma _lam _a _k _hd]

theorem WL_perpendicular_aggregate0_status2_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_perpendicular_aggregate0_status2 _rho _mu _sigma _lam _a _k ↔ ∃ w : Vec 3, WL_perpendicular_aggregate0_allowed2 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_perpendicular_aggregate0_status2_correct _rho _mu _sigma _lam _a _k _hd,
    WL_perpendicular_aggregate0_allowed2_correct _rho _mu _sigma _lam _a _ _hd]
  simp

theorem WL_perpendicular_aggregate0_outcome2_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_perpendicular_aggregate0_outcome2 _rho _mu _sigma _lam _a _k ↔ ∃ w : Vec 3, WL_perpendicular_aggregate0_allowed2 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_perpendicular_aggregate0_outcome2_correct _rho _mu _sigma _lam _a _k _hd,
    WL_perpendicular_aggregate0_allowed2_correct _rho _mu _sigma _lam _a _ _hd]
  simp

def WL_perpendicular_aggregate0_pairs : List (ℕ × ℕ) := [(1, 2), (1, 3), (2, 3)]
theorem WL_perpendicular_aggregate0_pairs_complete : WL_perpendicular_aggregate0_pairs = [(1, 2), (1, 3), (2, 3)] := by decide

end
end S10Audit.CAS.Records
