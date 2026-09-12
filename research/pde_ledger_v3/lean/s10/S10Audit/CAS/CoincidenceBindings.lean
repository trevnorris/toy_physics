import S10Audit.CAS.CoincidenceSupport
import S10Audit.CAS.PYCoincidenceArithmeticGeneric
import S10Audit.CAS.PYCoincidenceArithmeticParallel
import S10Audit.CAS.PYCoincidenceArithmeticPerpendicular
import S10Audit.CAS.WLCoincidenceArithmeticGeneric
import S10Audit.CAS.WLCoincidenceArithmeticParallel
import S10Audit.CAS.WLCoincidenceArithmeticPerpendicular

set_option backward.isDefEq.respectTransparency false

namespace S10Audit.CAS.Coincidence
open S10Pilot S10Anisotropic
noncomputable section

theorem PY_generic_pair0_roots (_rho _mu _sigma _lam : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    referenceRoot _rho _mu _sigma _k 0 = referenceRoot _rho _mu _sigma _k 1 ↔
      coincidenceGeometry 0 _k := by
  rcases _hd with ⟨hr, hm, hs, hs1, _hl⟩
  have h := coincidence_pair_geometry _rho _mu _sigma _k 0 (ne_of_gt hr) (ne_of_gt hm) hs hs1
  exact h

theorem PY_generic_pair1_roots (_rho _mu _sigma _lam : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    referenceRoot _rho _mu _sigma _k 0 = referenceRoot _rho _mu _sigma _k 2 ↔
      coincidenceGeometry 1 _k := by
  rcases _hd with ⟨hr, hm, hs, hs1, _hl⟩
  have h := coincidence_pair_geometry _rho _mu _sigma _k 1 (ne_of_gt hr) (ne_of_gt hm) hs hs1
  exact h

theorem PY_generic_pair2_roots (_rho _mu _sigma _lam : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    referenceRoot _rho _mu _sigma _k 1 = referenceRoot _rho _mu _sigma _k 2 ↔
      coincidenceGeometry 2 _k := by
  rcases _hd with ⟨hr, hm, hs, hs1, _hl⟩
  have h := coincidence_pair_geometry _rho _mu _sigma _k 2 (ne_of_gt hr) (ne_of_gt hm) hs hs1
  exact h

def PY_generic_locus0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := (((_k 0 = 0) ∧ (_k 1 = 0)) ∧ (_k 2 = 0))

theorem PY_generic_locus0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_locus0 _rho _mu _sigma _lam _a _k ↔ (coincidenceGeometry 0 _k) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_generic_locus0

def PY_generic_locus1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := (((_k 0 = 0) ∧ (_k 1 = 0)) ∧ (_k 2 = 0))

theorem PY_generic_locus1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_locus1 _rho _mu _sigma _lam _a _k ↔ (coincidenceGeometry 1 _k) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_generic_locus1

def PY_generic_locus2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((_k 1 = 0) ∧ (_k 2 = 0))

theorem PY_generic_locus2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_locus2 _rho _mu _sigma _lam _a _k ↔ (coincidenceGeometry 2 _k) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_generic_locus2

def PY_generic_allowed0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := (((((_k 0 = 0) ∧ (_k 1 = 0)) ∧ (_k 2 = 0)) ∧ ((((_k 0 ^ 2) + (_k 1 ^ 2)) + (_k 2 ^ 2)) > 0)) ∧ (((((((((((((True ∧ (0 < (3 : ℝ))) ∧ (0 < _lam)) ∧ (0 < _mu)) ∧ (0 < _rho)) ∧ (0 < _sigma)) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ (_sigma ≠ 1)) ∧ (0 < (((_k 0 ^ 2) + (_k 1 ^ 2)) + (_k 2 ^ 2)))))

theorem PY_generic_allowed0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_allowed0 _rho _mu _sigma _lam _a _k ↔ (allowedCoincidence 0 _k) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_generic_allowed0

def PY_generic_allowed1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := (((((_k 0 = 0) ∧ (_k 1 = 0)) ∧ (_k 2 = 0)) ∧ ((((_k 0 ^ 2) + (_k 1 ^ 2)) + (_k 2 ^ 2)) > 0)) ∧ (((((((((((((True ∧ (0 < (3 : ℝ))) ∧ (0 < _lam)) ∧ (0 < _mu)) ∧ (0 < _rho)) ∧ (0 < _sigma)) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ (_sigma ≠ 1)) ∧ (0 < (((_k 0 ^ 2) + (_k 1 ^ 2)) + (_k 2 ^ 2)))))

theorem PY_generic_allowed1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_allowed1 _rho _mu _sigma _lam _a _k ↔ (allowedCoincidence 1 _k) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_generic_allowed1

def PY_generic_allowed2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((((_k 1 = 0) ∧ (_k 2 = 0)) ∧ ((((_k 0 ^ 2) + (_k 1 ^ 2)) + (_k 2 ^ 2)) > 0)) ∧ (((((((((((((True ∧ (0 < (3 : ℝ))) ∧ (0 < _lam)) ∧ (0 < _mu)) ∧ (0 < _rho)) ∧ (0 < _sigma)) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ (_sigma ≠ 1)) ∧ (0 < (((_k 0 ^ 2) + (_k 1 ^ 2)) + (_k 2 ^ 2)))))

theorem PY_generic_allowed2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_allowed2 _rho _mu _sigma _lam _a _k ↔ (allowedCoincidence 2 _k) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_generic_allowed2

def PY_generic_status0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_generic_status0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_status0 _rho _mu _sigma _lam _a _k ↔ ((∃ w : Vec 3, allowedCoincidence 0 w)) := by
  rw [allowedCoincidence_exists]
  norm_num [PY_generic_status0, Fin.ext_iff, Fin.coe_ofNat_eq_mod]

def PY_generic_status1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_generic_status1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_status1 _rho _mu _sigma _lam _a _k ↔ ((∃ w : Vec 3, allowedCoincidence 1 w)) := by
  rw [allowedCoincidence_exists]
  norm_num [PY_generic_status1, Fin.ext_iff, Fin.coe_ofNat_eq_mod]

def PY_generic_status2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := True

theorem PY_generic_status2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_status2 _rho _mu _sigma _lam _a _k ↔ ((∃ w : Vec 3, allowedCoincidence 2 w)) := by
  rw [allowedCoincidence_exists]
  norm_num [PY_generic_status2, Fin.ext_iff, Fin.coe_ofNat_eq_mod]

def PY_generic_witness0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_generic_witness0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_witness0 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_generic_witness0

theorem PY_generic_witness0_sound (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam)
    (hw : PY_generic_witness0 _rho _mu _sigma _lam _a _k) : allowedCoincidence 0 _k := by
  rw [PY_generic_witness0_correct _rho _mu _sigma _lam _a _k _hd] at hw
  exact hw.elim

def PY_generic_witness1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_generic_witness1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_witness1 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_generic_witness1

theorem PY_generic_witness1_sound (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam)
    (hw : PY_generic_witness1 _rho _mu _sigma _lam _a _k) : allowedCoincidence 1 _k := by
  rw [PY_generic_witness1_correct _rho _mu _sigma _lam _a _k _hd] at hw
  exact hw.elim

def PY_generic_witness2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := (((_k 0 = 1) ∧ (_k 1 = 0)) ∧ (_k 2 = 0))

theorem PY_generic_witness2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_witness2 _rho _mu _sigma _lam _a _k ↔ ((_k 0 = 1 ∧ _k 1 = 0 ∧ _k 2 = 0)) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_generic_witness2

theorem PY_generic_witness2_sound (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam)
    (hw : PY_generic_witness2 _rho _mu _sigma _lam _a _k) : allowedCoincidence 2 _k := by
  rw [PY_generic_witness2_correct _rho _mu _sigma _lam _a _k _hd] at hw
  rw [allowedCoincidence_geometry]
  rcases hw with ⟨h0, h1, h2⟩
  simp [h0, h1, h2]

theorem PY_generic_witness2_exists (_rho _mu _sigma _lam : ℝ) (_a : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    ∃ _k : Vec 3, PY_generic_witness2 _rho _mu _sigma _lam _a _k := by
  refine ⟨![1, 0, 0], ?_⟩
  rw [PY_generic_witness2_correct _rho _mu _sigma _lam _a _ _hd]
  norm_num [Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]

def PY_generic_branch_allowed0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := (((((_k 0 = 0) ∧ (_k 1 = 0)) ∧ (_k 2 = 0)) ∧ ((((_k 0 ^ 2) + (_k 1 ^ 2)) + (_k 2 ^ 2)) > 0)) ∧ (((((((((((((True ∧ (0 < (3 : ℝ))) ∧ (0 < _lam)) ∧ (0 < _mu)) ∧ (0 < _rho)) ∧ (0 < _sigma)) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ (_sigma ≠ 1)) ∧ (0 < (((_k 0 ^ 2) + (_k 1 ^ 2)) + (_k 2 ^ 2)))))

theorem PY_generic_branch_allowed0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_branch_allowed0 _rho _mu _sigma _lam _a _k ↔ (allowedCoincidence 0 _k) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_generic_branch_allowed0

def PY_generic_branch_allowed1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := (((((_k 0 = 0) ∧ (_k 1 = 0)) ∧ (_k 2 = 0)) ∧ ((((_k 0 ^ 2) + (_k 1 ^ 2)) + (_k 2 ^ 2)) > 0)) ∧ (((((((((((((True ∧ (0 < (3 : ℝ))) ∧ (0 < _lam)) ∧ (0 < _mu)) ∧ (0 < _rho)) ∧ (0 < _sigma)) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ (_sigma ≠ 1)) ∧ (0 < (((_k 0 ^ 2) + (_k 1 ^ 2)) + (_k 2 ^ 2)))))

theorem PY_generic_branch_allowed1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_branch_allowed1 _rho _mu _sigma _lam _a _k ↔ (allowedCoincidence 1 _k) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_generic_branch_allowed1

def PY_generic_branch_allowed2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((((_k 1 = 0) ∧ (_k 2 = 0)) ∧ ((((_k 0 ^ 2) + (_k 1 ^ 2)) + (_k 2 ^ 2)) > 0)) ∧ (((((((((((((True ∧ (0 < (3 : ℝ))) ∧ (0 < _lam)) ∧ (0 < _mu)) ∧ (0 < _rho)) ∧ (0 < _sigma)) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ (_sigma ≠ 1)) ∧ (0 < (((_k 0 ^ 2) + (_k 1 ^ 2)) + (_k 2 ^ 2)))))

theorem PY_generic_branch_allowed2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_branch_allowed2 _rho _mu _sigma _lam _a _k ↔ (allowedCoincidence 2 _k) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_generic_branch_allowed2

def PY_generic_branch_status0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_generic_branch_status0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_branch_status0 _rho _mu _sigma _lam _a _k ↔ ((∃ w : Vec 3, allowedCoincidence 0 w)) := by
  rw [allowedCoincidence_exists]
  norm_num [PY_generic_branch_status0, Fin.ext_iff, Fin.coe_ofNat_eq_mod]

def PY_generic_branch_status1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_generic_branch_status1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_branch_status1 _rho _mu _sigma _lam _a _k ↔ ((∃ w : Vec 3, allowedCoincidence 1 w)) := by
  rw [allowedCoincidence_exists]
  norm_num [PY_generic_branch_status1, Fin.ext_iff, Fin.coe_ofNat_eq_mod]

def PY_generic_branch_status2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := True

theorem PY_generic_branch_status2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_branch_status2 _rho _mu _sigma _lam _a _k ↔ ((∃ w : Vec 3, allowedCoincidence 2 w)) := by
  rw [allowedCoincidence_exists]
  norm_num [PY_generic_branch_status2, Fin.ext_iff, Fin.coe_ofNat_eq_mod]

def PY_generic_branch_witness0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_generic_branch_witness0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_branch_witness0 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_generic_branch_witness0

def PY_generic_branch_witness1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_generic_branch_witness1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_branch_witness1 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_generic_branch_witness1

def PY_generic_branch_witness2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := (((_k 0 = 1) ∧ (_k 1 = 0)) ∧ (_k 2 = 0))

theorem PY_generic_branch_witness2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_branch_witness2 _rho _mu _sigma _lam _a _k ↔ (_k 0 = 1 ∧ _k 1 = 0 ∧ _k 2 = 0) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_generic_branch_witness2

theorem PY_parallel_pair0_roots (_rho _mu _sigma _lam : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    referenceRoot _rho _mu _sigma PYLoci.parallelPoint 0 = referenceRoot _rho _mu _sigma PYLoci.parallelPoint 1 ↔
      False := by
  rcases _hd with ⟨hr, hm, hs, hs1, _hl⟩
  have h := coincidence_pair_geometry _rho _mu _sigma PYLoci.parallelPoint 0 (ne_of_gt hr) (ne_of_gt hm) hs hs1
  norm_num [pairLeft, pairRight, coincidenceGeometry, PYLoci.parallelPoint,
    Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons, Fin.ext_iff, Fin.coe_ofNat_eq_mod] at h
  norm_num [PYLoci.parallelPoint]
  exact h

def PY_parallel_locus0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_parallel_locus0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_parallel_locus0 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_parallel_locus0

def PY_parallel_allowed0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((False ∧ (1 > 0)) ∧ True)

theorem PY_parallel_allowed0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_parallel_allowed0 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_parallel_allowed0

def PY_parallel_status0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_parallel_status0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_parallel_status0 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_parallel_status0

def PY_parallel_witness0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_parallel_witness0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_parallel_witness0 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_parallel_witness0

theorem PY_parallel_witness0_sound (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam)
    (hw : PY_parallel_witness0 _rho _mu _sigma _lam _a _k) : False := by
  rw [PY_parallel_witness0_correct _rho _mu _sigma _lam _a _k _hd] at hw
  exact hw.elim

def PY_parallel_branch_allowed0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_parallel_branch_allowed0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_parallel_branch_allowed0 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_parallel_branch_allowed0

def PY_parallel_branch_status0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_parallel_branch_status0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_parallel_branch_status0 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_parallel_branch_status0

def PY_parallel_branch_witness0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_parallel_branch_witness0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_parallel_branch_witness0 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_parallel_branch_witness0

theorem PY_perpendicular_pair0_roots (_rho _mu _sigma _lam : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    referenceRoot _rho _mu _sigma PYLoci.perpendicularPoint 0 = referenceRoot _rho _mu _sigma PYLoci.perpendicularPoint 1 ↔
      False := by
  rcases _hd with ⟨hr, hm, hs, hs1, _hl⟩
  have h := coincidence_pair_geometry _rho _mu _sigma PYLoci.perpendicularPoint 0 (ne_of_gt hr) (ne_of_gt hm) hs hs1
  norm_num [pairLeft, pairRight, coincidenceGeometry, PYLoci.perpendicularPoint,
    Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons, Fin.ext_iff, Fin.coe_ofNat_eq_mod] at h
  norm_num [PYLoci.perpendicularPoint]
  exact h

theorem PY_perpendicular_pair1_roots (_rho _mu _sigma _lam : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    referenceRoot _rho _mu _sigma PYLoci.perpendicularPoint 0 = referenceRoot _rho _mu _sigma PYLoci.perpendicularPoint 2 ↔
      False := by
  rcases _hd with ⟨hr, hm, hs, hs1, _hl⟩
  have h := coincidence_pair_geometry _rho _mu _sigma PYLoci.perpendicularPoint 1 (ne_of_gt hr) (ne_of_gt hm) hs hs1
  norm_num [pairLeft, pairRight, coincidenceGeometry, PYLoci.perpendicularPoint,
    Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons, Fin.ext_iff, Fin.coe_ofNat_eq_mod] at h
  norm_num [PYLoci.perpendicularPoint]
  exact h

theorem PY_perpendicular_pair2_roots (_rho _mu _sigma _lam : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    referenceRoot _rho _mu _sigma PYLoci.perpendicularPoint 1 = referenceRoot _rho _mu _sigma PYLoci.perpendicularPoint 2 ↔
      False := by
  rcases _hd with ⟨hr, hm, hs, hs1, _hl⟩
  have h := coincidence_pair_geometry _rho _mu _sigma PYLoci.perpendicularPoint 2 (ne_of_gt hr) (ne_of_gt hm) hs hs1
  norm_num [pairLeft, pairRight, coincidenceGeometry, PYLoci.perpendicularPoint,
    Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons, Fin.ext_iff, Fin.coe_ofNat_eq_mod] at h
  norm_num [PYLoci.perpendicularPoint]
  exact h

def PY_perpendicular_locus0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_perpendicular_locus0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_perpendicular_locus0 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_perpendicular_locus0

def PY_perpendicular_locus1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_perpendicular_locus1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_perpendicular_locus1 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_perpendicular_locus1

def PY_perpendicular_locus2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_perpendicular_locus2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_perpendicular_locus2 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_perpendicular_locus2

def PY_perpendicular_allowed0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((False ∧ (1 > 0)) ∧ True)

theorem PY_perpendicular_allowed0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_perpendicular_allowed0 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_perpendicular_allowed0

def PY_perpendicular_allowed1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((False ∧ (1 > 0)) ∧ True)

theorem PY_perpendicular_allowed1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_perpendicular_allowed1 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_perpendicular_allowed1

def PY_perpendicular_allowed2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((False ∧ (1 > 0)) ∧ True)

theorem PY_perpendicular_allowed2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_perpendicular_allowed2 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_perpendicular_allowed2

def PY_perpendicular_status0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_perpendicular_status0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_perpendicular_status0 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_perpendicular_status0

def PY_perpendicular_status1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_perpendicular_status1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_perpendicular_status1 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_perpendicular_status1

def PY_perpendicular_status2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_perpendicular_status2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_perpendicular_status2 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_perpendicular_status2

def PY_perpendicular_witness0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_perpendicular_witness0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_perpendicular_witness0 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_perpendicular_witness0

theorem PY_perpendicular_witness0_sound (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam)
    (hw : PY_perpendicular_witness0 _rho _mu _sigma _lam _a _k) : False := by
  rw [PY_perpendicular_witness0_correct _rho _mu _sigma _lam _a _k _hd] at hw
  exact hw.elim

def PY_perpendicular_witness1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_perpendicular_witness1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_perpendicular_witness1 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_perpendicular_witness1

theorem PY_perpendicular_witness1_sound (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam)
    (hw : PY_perpendicular_witness1 _rho _mu _sigma _lam _a _k) : False := by
  rw [PY_perpendicular_witness1_correct _rho _mu _sigma _lam _a _k _hd] at hw
  exact hw.elim

def PY_perpendicular_witness2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_perpendicular_witness2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_perpendicular_witness2 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_perpendicular_witness2

theorem PY_perpendicular_witness2_sound (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam)
    (hw : PY_perpendicular_witness2 _rho _mu _sigma _lam _a _k) : False := by
  rw [PY_perpendicular_witness2_correct _rho _mu _sigma _lam _a _k _hd] at hw
  exact hw.elim

def PY_perpendicular_branch_allowed0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_perpendicular_branch_allowed0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_perpendicular_branch_allowed0 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_perpendicular_branch_allowed0

def PY_perpendicular_branch_allowed1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_perpendicular_branch_allowed1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_perpendicular_branch_allowed1 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_perpendicular_branch_allowed1

def PY_perpendicular_branch_allowed2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_perpendicular_branch_allowed2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_perpendicular_branch_allowed2 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_perpendicular_branch_allowed2

def PY_perpendicular_branch_status0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_perpendicular_branch_status0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_perpendicular_branch_status0 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_perpendicular_branch_status0

def PY_perpendicular_branch_status1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_perpendicular_branch_status1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_perpendicular_branch_status1 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_perpendicular_branch_status1

def PY_perpendicular_branch_status2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_perpendicular_branch_status2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_perpendicular_branch_status2 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_perpendicular_branch_status2

def PY_perpendicular_branch_witness0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_perpendicular_branch_witness0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_perpendicular_branch_witness0 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_perpendicular_branch_witness0

def PY_perpendicular_branch_witness1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_perpendicular_branch_witness1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_perpendicular_branch_witness1 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_perpendicular_branch_witness1

def PY_perpendicular_branch_witness2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem PY_perpendicular_branch_witness2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_perpendicular_branch_witness2 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval PY_perpendicular_branch_witness2

theorem WL_generic_pair0_roots (_rho _mu _sigma _lam : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    referenceRoot _rho _mu _sigma _k 0 = referenceRoot _rho _mu _sigma _k 1 ↔
      coincidenceGeometry 0 _k := by
  rcases _hd with ⟨hr, hm, hs, hs1, _hl⟩
  have h := coincidence_pair_geometry _rho _mu _sigma _k 0 (ne_of_gt hr) (ne_of_gt hm) hs hs1
  exact h

theorem WL_generic_pair1_roots (_rho _mu _sigma _lam : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    referenceRoot _rho _mu _sigma _k 0 = referenceRoot _rho _mu _sigma _k 2 ↔
      coincidenceGeometry 1 _k := by
  rcases _hd with ⟨hr, hm, hs, hs1, _hl⟩
  have h := coincidence_pair_geometry _rho _mu _sigma _k 1 (ne_of_gt hr) (ne_of_gt hm) hs hs1
  exact h

theorem WL_generic_pair2_roots (_rho _mu _sigma _lam : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    referenceRoot _rho _mu _sigma _k 1 = referenceRoot _rho _mu _sigma _k 2 ↔
      coincidenceGeometry 2 _k := by
  rcases _hd with ⟨hr, hm, hs, hs1, _hl⟩
  have h := coincidence_pair_geometry _rho _mu _sigma _k 2 (ne_of_gt hr) (ne_of_gt hm) hs hs1
  exact h

def WL_generic_premise0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((((((((((((((((_rho > 0) ∧ (_mu > 0)) ∧ (_lam > 0)) ∧ ((((_k 0 ^ 2) + (_k 1 ^ 2)) + (_k 2 ^ 2)) > 0)) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ ((3 : ℝ) > 0)) ∧ (_sigma > 0)) ∧ True) ∧ (_sigma ≠ 1))

theorem WL_generic_premise0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_premise0 _rho _mu _sigma _lam _a _k ↔ (_k ≠ 0) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval WL_generic_premise0

def WL_generic_locus0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((((_k 0 = 0) ∧ ((_rho > 0) ∨ (_rho < 0))) ∧ ((_k 1 = 0) ∧ ((_rho > 0) ∨ (_rho < 0)))) ∧ ((_k 2 = 0) ∧ ((_rho > 0) ∨ (_rho < 0))))

theorem WL_generic_locus0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_locus0 _rho _mu _sigma _lam _a _k ↔ (coincidenceGeometry 0 _k) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval WL_generic_locus0

def WL_generic_allowed0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_generic_allowed0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_allowed0 _rho _mu _sigma _lam _a _k ↔ (allowedCoincidence 0 _k) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval WL_generic_allowed0

def WL_generic_status0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_generic_status0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_status0 _rho _mu _sigma _lam _a _k ↔ ((∃ w : Vec 3, allowedCoincidence 0 w)) := by
  rw [allowedCoincidence_exists]
  norm_num [WL_generic_status0, Fin.ext_iff, Fin.coe_ofNat_eq_mod]

def WL_generic_outcome0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_generic_outcome0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_outcome0 _rho _mu _sigma _lam _a _k ↔ ((∃ w : Vec 3, allowedCoincidence 0 w)) := by
  rw [allowedCoincidence_exists]
  norm_num [WL_generic_outcome0, Fin.ext_iff, Fin.coe_ofNat_eq_mod]

def WL_generic_premise1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((((((((((((((((_rho > 0) ∧ (_mu > 0)) ∧ (_lam > 0)) ∧ ((((_k 0 ^ 2) + (_k 1 ^ 2)) + (_k 2 ^ 2)) > 0)) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ ((3 : ℝ) > 0)) ∧ (_sigma > 0)) ∧ True) ∧ (_sigma ≠ 1))

theorem WL_generic_premise1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_premise1 _rho _mu _sigma _lam _a _k ↔ (_k ≠ 0) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval WL_generic_premise1

def WL_generic_locus1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((((((_k 2 = (-(Real.sqrt ((-(_k 1 ^ 2)) - ((_k 0 ^ 2) * _sigma))))) ∧ (((_k 0 < (-(Real.sqrt (-((_k 1 ^ 2) / _sigma))))) ∧ (_sigma < 0)) ∨ ((_sigma < 0) ∧ (_k 0 > (Real.sqrt (-((_k 1 ^ 2) / _sigma))))))) ∨ ((_k 2 = (Real.sqrt ((-(_k 1 ^ 2)) - ((_k 0 ^ 2) * _sigma)))) ∧ (((_k 0 < (-(Real.sqrt (-((_k 1 ^ 2) / _sigma))))) ∧ (_sigma < 0)) ∨ ((_sigma < 0) ∧ (_k 0 > (Real.sqrt (-((_k 1 ^ 2) / _sigma)))))))) ∨ (((_k 0 = (-(Real.sqrt (-((_k 1 ^ 2) / _sigma))))) ∧ (_sigma < 0)) ∧ ((_k 2 = 0) ∧ (_sigma < 0)))) ∨ (((_k 0 = (Real.sqrt (-((_k 1 ^ 2) / _sigma)))) ∧ (_sigma < 0)) ∧ ((_k 2 = 0) ∧ (_sigma < 0)))) ∨ ((((_k 0 = 0) ∧ (_sigma > 0)) ∧ ((_k 1 = 0) ∧ (_sigma > 0))) ∧ ((_k 2 = 0) ∧ (_sigma > 0))))

theorem WL_generic_locus1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_locus1 _rho _mu _sigma _lam _a _k ↔ (coincidenceGeometry 1 _k) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval WL_generic_locus1

def WL_generic_allowed1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_generic_allowed1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_allowed1 _rho _mu _sigma _lam _a _k ↔ (allowedCoincidence 1 _k) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval WL_generic_allowed1

def WL_generic_status1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_generic_status1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_status1 _rho _mu _sigma _lam _a _k ↔ ((∃ w : Vec 3, allowedCoincidence 1 w)) := by
  rw [allowedCoincidence_exists]
  norm_num [WL_generic_status1, Fin.ext_iff, Fin.coe_ofNat_eq_mod]

def WL_generic_outcome1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_generic_outcome1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_outcome1 _rho _mu _sigma _lam _a _k ↔ ((∃ w : Vec 3, allowedCoincidence 1 w)) := by
  rw [allowedCoincidence_exists]
  norm_num [WL_generic_outcome1, Fin.ext_iff, Fin.coe_ofNat_eq_mod]

def WL_generic_premise2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((((((((((((((((_rho > 0) ∧ (_mu > 0)) ∧ (_lam > 0)) ∧ ((((_k 0 ^ 2) + (_k 1 ^ 2)) + (_k 2 ^ 2)) > 0)) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ ((3 : ℝ) > 0)) ∧ (_sigma > 0)) ∧ True) ∧ (_sigma ≠ 1))

theorem WL_generic_premise2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_premise2 _rho _mu _sigma _lam _a _k ↔ (_k ≠ 0) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval WL_generic_premise2

def WL_generic_locus2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((_k 1 = 0) ∧ (_k 2 = 0))

theorem WL_generic_locus2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_locus2 _rho _mu _sigma _lam _a _k ↔ (coincidenceGeometry 2 _k) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval WL_generic_locus2

def WL_generic_allowed2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((((((True ∧ ((3 : ℝ) ≥ 1)) ∧ (_lam > 0)) ∧ (_mu > 0)) ∧ (_rho > 0)) ∧ (((0 < _sigma) ∧ (_sigma < 1)) ∨ (_sigma > 1))) ∧ ((((_k 0 < 0) ∧ (_k 1 = 0)) ∧ (_k 2 = 0)) ∨ (((_k 0 > 0) ∧ (_k 1 = 0)) ∧ (_k 2 = 0))))

theorem WL_generic_allowed2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_allowed2 _rho _mu _sigma _lam _a _k ↔ (allowedCoincidence 2 _k) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval WL_generic_allowed2

def WL_generic_status2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := True

theorem WL_generic_status2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_status2 _rho _mu _sigma _lam _a _k ↔ ((∃ w : Vec 3, allowedCoincidence 2 w)) := by
  rw [allowedCoincidence_exists]
  norm_num [WL_generic_status2, Fin.ext_iff, Fin.coe_ofNat_eq_mod]

def WL_generic_outcome2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := True

theorem WL_generic_outcome2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_outcome2 _rho _mu _sigma _lam _a _k ↔ ((∃ w : Vec 3, allowedCoincidence 2 w)) := by
  rw [allowedCoincidence_exists]
  norm_num [WL_generic_outcome2, Fin.ext_iff, Fin.coe_ofNat_eq_mod]

theorem WL_parallel_pair0_roots (_rho _mu _sigma _lam : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    referenceRoot _rho _mu _sigma WLLoci.parallelPoint 0 = referenceRoot _rho _mu _sigma WLLoci.parallelPoint 1 ↔
      False := by
  rcases _hd with ⟨hr, hm, hs, hs1, _hl⟩
  have h := coincidence_pair_geometry _rho _mu _sigma WLLoci.parallelPoint 0 (ne_of_gt hr) (ne_of_gt hm) hs hs1
  norm_num [pairLeft, pairRight, coincidenceGeometry, WLLoci.parallelPoint,
    Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons, Fin.ext_iff, Fin.coe_ofNat_eq_mod] at h
  norm_num [WLLoci.parallelPoint]
  exact h

def WL_parallel_premise0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((((((((((((((_rho > 0) ∧ (_mu > 0)) ∧ (_lam > 0)) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ ((3 : ℝ) > 0)) ∧ (_sigma > 0)) ∧ True) ∧ (_sigma ≠ 1)) ∧ ((((((True ∧ (_rho > 0)) ∧ (_mu > 0)) ∧ (_lam > 0)) ∧ ((3 : ℝ) ≥ 1)) ∧ (_sigma > 1)) ∨ (((((True ∧ (_rho > 0)) ∧ (_mu > 0)) ∧ (_lam > 0)) ∧ ((3 : ℝ) ≥ 1)) ∧ ((0 < _sigma) ∧ (_sigma < 1))))) ∧ ((3 : ℝ) = 3))

theorem WL_parallel_premise0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_parallel_premise0 _rho _mu _sigma _lam _a _k ↔ (True) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval WL_parallel_premise0

def WL_parallel_locus0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((_mu = 0) ∧ ((_rho > 0) ∨ (_rho < 0)))

theorem WL_parallel_locus0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_parallel_locus0 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval WL_parallel_locus0

def WL_parallel_allowed0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_parallel_allowed0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_parallel_allowed0 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval WL_parallel_allowed0

def WL_parallel_status0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_parallel_status0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_parallel_status0 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval WL_parallel_status0

def WL_parallel_outcome0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_parallel_outcome0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_parallel_outcome0 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval WL_parallel_outcome0

theorem WL_perpendicular_pair0_roots (_rho _mu _sigma _lam : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    referenceRoot _rho _mu _sigma WLLoci.perpendicularPoint 0 = referenceRoot _rho _mu _sigma WLLoci.perpendicularPoint 1 ↔
      False := by
  rcases _hd with ⟨hr, hm, hs, hs1, _hl⟩
  have h := coincidence_pair_geometry _rho _mu _sigma WLLoci.perpendicularPoint 0 (ne_of_gt hr) (ne_of_gt hm) hs hs1
  norm_num [pairLeft, pairRight, coincidenceGeometry, WLLoci.perpendicularPoint,
    Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons, Fin.ext_iff, Fin.coe_ofNat_eq_mod] at h
  norm_num [WLLoci.perpendicularPoint]
  exact h

theorem WL_perpendicular_pair1_roots (_rho _mu _sigma _lam : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    referenceRoot _rho _mu _sigma WLLoci.perpendicularPoint 0 = referenceRoot _rho _mu _sigma WLLoci.perpendicularPoint 2 ↔
      False := by
  rcases _hd with ⟨hr, hm, hs, hs1, _hl⟩
  have h := coincidence_pair_geometry _rho _mu _sigma WLLoci.perpendicularPoint 1 (ne_of_gt hr) (ne_of_gt hm) hs hs1
  norm_num [pairLeft, pairRight, coincidenceGeometry, WLLoci.perpendicularPoint,
    Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons, Fin.ext_iff, Fin.coe_ofNat_eq_mod] at h
  norm_num [WLLoci.perpendicularPoint]
  exact h

theorem WL_perpendicular_pair2_roots (_rho _mu _sigma _lam : ℝ) (_k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    referenceRoot _rho _mu _sigma WLLoci.perpendicularPoint 1 = referenceRoot _rho _mu _sigma WLLoci.perpendicularPoint 2 ↔
      False := by
  rcases _hd with ⟨hr, hm, hs, hs1, _hl⟩
  have h := coincidence_pair_geometry _rho _mu _sigma WLLoci.perpendicularPoint 2 (ne_of_gt hr) (ne_of_gt hm) hs hs1
  norm_num [pairLeft, pairRight, coincidenceGeometry, WLLoci.perpendicularPoint,
    Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons, Fin.ext_iff, Fin.coe_ofNat_eq_mod] at h
  norm_num [WLLoci.perpendicularPoint]
  exact h

def WL_perpendicular_premise0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((((((((((((((_rho > 0) ∧ (_mu > 0)) ∧ (_lam > 0)) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ ((3 : ℝ) > 0)) ∧ (_sigma > 0)) ∧ True) ∧ (_sigma ≠ 1)) ∧ ((((((True ∧ (_rho > 0)) ∧ (_mu > 0)) ∧ (_lam > 0)) ∧ ((3 : ℝ) ≥ 1)) ∧ (_sigma > 1)) ∨ (((((True ∧ (_rho > 0)) ∧ (_mu > 0)) ∧ (_lam > 0)) ∧ ((3 : ℝ) ≥ 1)) ∧ ((0 < _sigma) ∧ (_sigma < 1))))) ∧ ((3 : ℝ) = 3))

theorem WL_perpendicular_premise0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_perpendicular_premise0 _rho _mu _sigma _lam _a _k ↔ (True) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval WL_perpendicular_premise0

def WL_perpendicular_locus0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((_mu = 0) ∧ ((_rho > 0) ∨ (_rho < 0)))

theorem WL_perpendicular_locus0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_perpendicular_locus0 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval WL_perpendicular_locus0

def WL_perpendicular_allowed0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_perpendicular_allowed0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_perpendicular_allowed0 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval WL_perpendicular_allowed0

def WL_perpendicular_status0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_perpendicular_status0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_perpendicular_status0 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval WL_perpendicular_status0

def WL_perpendicular_outcome0 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_perpendicular_outcome0_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_perpendicular_outcome0 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval WL_perpendicular_outcome0

def WL_perpendicular_premise1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((((((((((((((_rho > 0) ∧ (_mu > 0)) ∧ (_lam > 0)) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ ((3 : ℝ) > 0)) ∧ (_sigma > 0)) ∧ True) ∧ (_sigma ≠ 1)) ∧ ((((((True ∧ (_rho > 0)) ∧ (_mu > 0)) ∧ (_lam > 0)) ∧ ((3 : ℝ) ≥ 1)) ∧ (_sigma > 1)) ∨ (((((True ∧ (_rho > 0)) ∧ (_mu > 0)) ∧ (_lam > 0)) ∧ ((3 : ℝ) ≥ 1)) ∧ ((0 < _sigma) ∧ (_sigma < 1))))) ∧ ((3 : ℝ) = 3))

theorem WL_perpendicular_premise1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_perpendicular_premise1 _rho _mu _sigma _lam _a _k ↔ (True) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval WL_perpendicular_premise1

def WL_perpendicular_locus1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((_mu = 0) ∧ (((((_rho > 0) ∧ (_sigma > 0)) ∨ ((_rho > 0) ∧ (_sigma < 0))) ∨ ((_rho < 0) ∧ (_sigma > 0))) ∨ ((_rho < 0) ∧ (_sigma < 0))))

theorem WL_perpendicular_locus1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_perpendicular_locus1 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval WL_perpendicular_locus1

def WL_perpendicular_allowed1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_perpendicular_allowed1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_perpendicular_allowed1 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval WL_perpendicular_allowed1

def WL_perpendicular_status1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_perpendicular_status1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_perpendicular_status1 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval WL_perpendicular_status1

def WL_perpendicular_outcome1 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_perpendicular_outcome1_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_perpendicular_outcome1 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval WL_perpendicular_outcome1

def WL_perpendicular_premise2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := ((((((((((((((_rho > 0) ∧ (_mu > 0)) ∧ (_lam > 0)) ∧ True) ∧ True) ∧ True) ∧ True) ∧ True) ∧ ((3 : ℝ) > 0)) ∧ (_sigma > 0)) ∧ True) ∧ (_sigma ≠ 1)) ∧ ((((((True ∧ (_rho > 0)) ∧ (_mu > 0)) ∧ (_lam > 0)) ∧ ((3 : ℝ) ≥ 1)) ∧ (_sigma > 1)) ∨ (((((True ∧ (_rho > 0)) ∧ (_mu > 0)) ∧ (_lam > 0)) ∧ ((3 : ℝ) ≥ 1)) ∧ ((0 < _sigma) ∧ (_sigma < 1))))) ∧ ((3 : ℝ) = 3))

theorem WL_perpendicular_premise2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_perpendicular_premise2 _rho _mu _sigma _lam _a _k ↔ (True) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval WL_perpendicular_premise2

def WL_perpendicular_locus2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := (((_mu = 0) ∧ ((((((((0 < _sigma) ∧ (_sigma < 1)) ∧ (_rho > 0)) ∨ (((0 < _sigma) ∧ (_sigma < 1)) ∧ (_rho < 0))) ∨ ((_sigma > 1) ∧ (_rho > 0))) ∨ ((_sigma > 1) ∧ (_rho < 0))) ∨ ((_sigma < 0) ∧ (_rho > 0))) ∨ ((_sigma < 0) ∧ (_rho < 0)))) ∨ ((_sigma = 1) ∧ ((_rho > 0) ∨ (_rho < 0))))

theorem WL_perpendicular_locus2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_perpendicular_locus2 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval WL_perpendicular_locus2

def WL_perpendicular_allowed2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_perpendicular_allowed2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_perpendicular_allowed2 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval WL_perpendicular_allowed2

def WL_perpendicular_status2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_perpendicular_status2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_perpendicular_status2 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval WL_perpendicular_status2

def WL_perpendicular_outcome2 (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) : Prop := False

theorem WL_perpendicular_outcome2_correct (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_perpendicular_outcome2 _rho _mu _sigma _lam _a _k ↔ (False) := by
  rcases _hd with ⟨hr, hm, hs, hs1, hl⟩
  have hrn : ¬ _rho < 0 := not_lt_of_gt hr
  have hsn : ¬ _sigma < 0 := not_lt_of_gt hs
  have hmn : _mu ≠ 0 := ne_of_gt hm
  have hsplit := positive_sigma_split _sigma hs hs1
  coincidence_eval WL_perpendicular_outcome2

theorem PY_generic_equation0_locus (_rho _mu _sigma _lam z : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    Expr.eval (values _rho _mu _sigma z _k) PYCoincidenceArithmeticGeneric.n13 = 0 ↔
      PY_generic_locus0 _rho _mu _sigma _lam _a _k := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [PYCoincidenceArithmeticGeneric.q3_root_coincidence_equations_cell0 _rho _mu _sigma z _k (ne_of_gt hr) (ne_of_gt hs), sub_eq_zero,
    PY_generic_pair0_roots _rho _mu _sigma _lam _k _hd,
    PY_generic_locus0_correct _rho _mu _sigma _lam _a _k _hd]

theorem PY_generic_status0_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_status0 _rho _mu _sigma _lam _a _k ↔
      ∃ w : Vec 3, PY_generic_allowed0 _rho _mu _sigma _lam _a w := by
  simp_rw [PY_generic_status0_correct _rho _mu _sigma _lam _a _k _hd,
    PY_generic_allowed0_correct _rho _mu _sigma _lam _a _ _hd]

theorem PY_generic_equation1_locus (_rho _mu _sigma _lam z : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    Expr.eval (values _rho _mu _sigma z _k) PYCoincidenceArithmeticGeneric.n20 = 0 ↔
      PY_generic_locus1 _rho _mu _sigma _lam _a _k := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [PYCoincidenceArithmeticGeneric.q3_root_coincidence_equations_cell1 _rho _mu _sigma z _k (ne_of_gt hr) (ne_of_gt hs), sub_eq_zero,
    PY_generic_pair1_roots _rho _mu _sigma _lam _k _hd,
    PY_generic_locus1_correct _rho _mu _sigma _lam _a _k _hd]

theorem PY_generic_status1_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_status1 _rho _mu _sigma _lam _a _k ↔
      ∃ w : Vec 3, PY_generic_allowed1 _rho _mu _sigma _lam _a w := by
  simp_rw [PY_generic_status1_correct _rho _mu _sigma _lam _a _k _hd,
    PY_generic_allowed1_correct _rho _mu _sigma _lam _a _ _hd]

theorem PY_generic_equation2_locus (_rho _mu _sigma _lam z : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    Expr.eval (values _rho _mu _sigma z _k) PYCoincidenceArithmeticGeneric.n27 = 0 ↔
      PY_generic_locus2 _rho _mu _sigma _lam _a _k := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [PYCoincidenceArithmeticGeneric.q3_root_coincidence_equations_cell2 _rho _mu _sigma z _k (ne_of_gt hr) (ne_of_gt hs), sub_eq_zero,
    PY_generic_pair2_roots _rho _mu _sigma _lam _k _hd,
    PY_generic_locus2_correct _rho _mu _sigma _lam _a _k _hd]

theorem PY_generic_status2_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_status2 _rho _mu _sigma _lam _a _k ↔
      ∃ w : Vec 3, PY_generic_allowed2 _rho _mu _sigma _lam _a w := by
  simp_rw [PY_generic_status2_correct _rho _mu _sigma _lam _a _k _hd,
    PY_generic_allowed2_correct _rho _mu _sigma _lam _a _ _hd]

theorem PY_parallel_equation0_locus (_rho _mu _sigma _lam z : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    Expr.eval (values _rho _mu _sigma z _k) PYCoincidenceArithmeticParallel.n4 = 0 ↔
      PY_parallel_locus0 _rho _mu _sigma _lam _a _k := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [PYCoincidenceArithmeticParallel.q8_stratum1_q3_root_coincidence_equations_cell0 _rho _mu _sigma z _k (ne_of_gt hr) (ne_of_gt hs), sub_eq_zero,
    PY_parallel_pair0_roots _rho _mu _sigma _lam _k _hd,
    PY_parallel_locus0_correct _rho _mu _sigma _lam _a _k _hd]

theorem PY_parallel_status0_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_parallel_status0 _rho _mu _sigma _lam _a _k ↔
      ∃ w : Vec 3, PY_parallel_allowed0 _rho _mu _sigma _lam _a w := by
  simp_rw [PY_parallel_status0_correct _rho _mu _sigma _lam _a _k _hd,
    PY_parallel_allowed0_correct _rho _mu _sigma _lam _a _ _hd]
  simp

theorem PY_perpendicular_equation0_locus (_rho _mu _sigma _lam z : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    Expr.eval (values _rho _mu _sigma z _k) PYCoincidenceArithmeticPerpendicular.n4 = 0 ↔
      PY_perpendicular_locus0 _rho _mu _sigma _lam _a _k := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [PYCoincidenceArithmeticPerpendicular.q8_stratum2_q3_root_coincidence_equations_cell0 _rho _mu _sigma z _k (ne_of_gt hr) (ne_of_gt hs), sub_eq_zero,
    PY_perpendicular_pair0_roots _rho _mu _sigma _lam _k _hd,
    PY_perpendicular_locus0_correct _rho _mu _sigma _lam _a _k _hd]

theorem PY_perpendicular_status0_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_perpendicular_status0 _rho _mu _sigma _lam _a _k ↔
      ∃ w : Vec 3, PY_perpendicular_allowed0 _rho _mu _sigma _lam _a w := by
  simp_rw [PY_perpendicular_status0_correct _rho _mu _sigma _lam _a _k _hd,
    PY_perpendicular_allowed0_correct _rho _mu _sigma _lam _a _ _hd]
  simp

theorem PY_perpendicular_equation1_locus (_rho _mu _sigma _lam z : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    Expr.eval (values _rho _mu _sigma z _k) PYCoincidenceArithmeticPerpendicular.n7 = 0 ↔
      PY_perpendicular_locus1 _rho _mu _sigma _lam _a _k := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [PYCoincidenceArithmeticPerpendicular.q8_stratum2_q3_root_coincidence_equations_cell1 _rho _mu _sigma z _k (ne_of_gt hr) (ne_of_gt hs), sub_eq_zero,
    PY_perpendicular_pair1_roots _rho _mu _sigma _lam _k _hd,
    PY_perpendicular_locus1_correct _rho _mu _sigma _lam _a _k _hd]

theorem PY_perpendicular_status1_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_perpendicular_status1 _rho _mu _sigma _lam _a _k ↔
      ∃ w : Vec 3, PY_perpendicular_allowed1 _rho _mu _sigma _lam _a w := by
  simp_rw [PY_perpendicular_status1_correct _rho _mu _sigma _lam _a _k _hd,
    PY_perpendicular_allowed1_correct _rho _mu _sigma _lam _a _ _hd]
  simp

theorem PY_perpendicular_equation2_locus (_rho _mu _sigma _lam z : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    Expr.eval (values _rho _mu _sigma z _k) PYCoincidenceArithmeticPerpendicular.n11 = 0 ↔
      PY_perpendicular_locus2 _rho _mu _sigma _lam _a _k := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [PYCoincidenceArithmeticPerpendicular.q8_stratum2_q3_root_coincidence_equations_cell2 _rho _mu _sigma z _k (ne_of_gt hr) (ne_of_gt hs), sub_eq_zero,
    PY_perpendicular_pair2_roots _rho _mu _sigma _lam _k _hd,
    PY_perpendicular_locus2_correct _rho _mu _sigma _lam _a _k _hd]

theorem PY_perpendicular_status2_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_perpendicular_status2 _rho _mu _sigma _lam _a _k ↔
      ∃ w : Vec 3, PY_perpendicular_allowed2 _rho _mu _sigma _lam _a w := by
  simp_rw [PY_perpendicular_status2_correct _rho _mu _sigma _lam _a _k _hd,
    PY_perpendicular_allowed2_correct _rho _mu _sigma _lam _a _ _hd]
  simp

theorem WL_generic_equation0_locus (_rho _mu _sigma _lam z : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    Expr.eval (values _rho _mu _sigma z _k) WLCoincidenceArithmeticGeneric.n13 = 0 ↔
      WL_generic_locus0 _rho _mu _sigma _lam _a _k := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [WLCoincidenceArithmeticGeneric.root1_root2_q3_coincidence_operands_cell0 _rho _mu _sigma z _k (ne_of_gt hr) (ne_of_gt hs), sub_eq_zero,
    WL_generic_pair0_roots _rho _mu _sigma _lam _k _hd,
    WL_generic_locus0_correct _rho _mu _sigma _lam _a _k _hd]

theorem WL_generic_status0_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_status0 _rho _mu _sigma _lam _a _k ↔
      ∃ w : Vec 3, WL_generic_allowed0 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_generic_status0_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_allowed0_correct _rho _mu _sigma _lam _a _ _hd]

theorem WL_generic_outcome0_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_outcome0 _rho _mu _sigma _lam _a _k ↔
      ∃ w : Vec 3, WL_generic_allowed0 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_generic_outcome0_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_allowed0_correct _rho _mu _sigma _lam _a _ _hd]

theorem WL_generic_equation1_locus (_rho _mu _sigma _lam z : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    Expr.eval (values _rho _mu _sigma z _k) WLCoincidenceArithmeticGeneric.n21 = 0 ↔
      WL_generic_locus1 _rho _mu _sigma _lam _a _k := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [WLCoincidenceArithmeticGeneric.root1_root3_q3_coincidence_operands_cell0 _rho _mu _sigma z _k (ne_of_gt hr) (ne_of_gt hs), sub_eq_zero,
    WL_generic_pair1_roots _rho _mu _sigma _lam _k _hd,
    WL_generic_locus1_correct _rho _mu _sigma _lam _a _k _hd]

theorem WL_generic_status1_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_status1 _rho _mu _sigma _lam _a _k ↔
      ∃ w : Vec 3, WL_generic_allowed1 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_generic_status1_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_allowed1_correct _rho _mu _sigma _lam _a _ _hd]

theorem WL_generic_outcome1_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_outcome1 _rho _mu _sigma _lam _a _k ↔
      ∃ w : Vec 3, WL_generic_allowed1 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_generic_outcome1_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_allowed1_correct _rho _mu _sigma _lam _a _ _hd]

theorem WL_generic_equation2_locus (_rho _mu _sigma _lam z : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    Expr.eval (values _rho _mu _sigma z _k) WLCoincidenceArithmeticGeneric.n29 = 0 ↔
      WL_generic_locus2 _rho _mu _sigma _lam _a _k := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [WLCoincidenceArithmeticGeneric.root2_root3_q3_coincidence_operands_cell0 _rho _mu _sigma z _k (ne_of_gt hr) (ne_of_gt hs), sub_eq_zero,
    WL_generic_pair2_roots _rho _mu _sigma _lam _k _hd,
    WL_generic_locus2_correct _rho _mu _sigma _lam _a _k _hd]

theorem WL_generic_status2_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_status2 _rho _mu _sigma _lam _a _k ↔
      ∃ w : Vec 3, WL_generic_allowed2 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_generic_status2_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_allowed2_correct _rho _mu _sigma _lam _a _ _hd]

theorem WL_generic_outcome2_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_generic_outcome2 _rho _mu _sigma _lam _a _k ↔
      ∃ w : Vec 3, WL_generic_allowed2 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_generic_outcome2_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_allowed2_correct _rho _mu _sigma _lam _a _ _hd]

theorem WL_parallel_equation0_locus (_rho _mu _sigma _lam z : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    Expr.eval (values _rho _mu _sigma z _k) WLCoincidenceArithmeticParallel.n6 = 0 ↔
      WL_parallel_locus0 _rho _mu _sigma _lam _a _k := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [WLCoincidenceArithmeticParallel.stratum1_root1_root2_q3_coincidence_operands_cell0 _rho _mu _sigma z _k (ne_of_gt hr) (ne_of_gt hs), sub_eq_zero,
    WL_parallel_pair0_roots _rho _mu _sigma _lam _k _hd,
    WL_parallel_locus0_correct _rho _mu _sigma _lam _a _k _hd]

theorem WL_parallel_status0_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_parallel_status0 _rho _mu _sigma _lam _a _k ↔
      ∃ w : Vec 3, WL_parallel_allowed0 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_parallel_status0_correct _rho _mu _sigma _lam _a _k _hd,
    WL_parallel_allowed0_correct _rho _mu _sigma _lam _a _ _hd]
  simp

theorem WL_parallel_outcome0_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_parallel_outcome0 _rho _mu _sigma _lam _a _k ↔
      ∃ w : Vec 3, WL_parallel_allowed0 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_parallel_outcome0_correct _rho _mu _sigma _lam _a _k _hd,
    WL_parallel_allowed0_correct _rho _mu _sigma _lam _a _ _hd]
  simp

theorem WL_perpendicular_equation0_locus (_rho _mu _sigma _lam z : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    Expr.eval (values _rho _mu _sigma z _k) WLCoincidenceArithmeticPerpendicular.n8 = 0 ↔
      WL_perpendicular_locus0 _rho _mu _sigma _lam _a _k := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [WLCoincidenceArithmeticPerpendicular.stratum2_root1_root2_q3_coincidence_operands_cell0 _rho _mu _sigma z _k (ne_of_gt hr) (ne_of_gt hs), sub_eq_zero,
    WL_perpendicular_pair0_roots _rho _mu _sigma _lam _k _hd,
    WL_perpendicular_locus0_correct _rho _mu _sigma _lam _a _k _hd]

theorem WL_perpendicular_status0_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_perpendicular_status0 _rho _mu _sigma _lam _a _k ↔
      ∃ w : Vec 3, WL_perpendicular_allowed0 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_perpendicular_status0_correct _rho _mu _sigma _lam _a _k _hd,
    WL_perpendicular_allowed0_correct _rho _mu _sigma _lam _a _ _hd]
  simp

theorem WL_perpendicular_outcome0_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_perpendicular_outcome0 _rho _mu _sigma _lam _a _k ↔
      ∃ w : Vec 3, WL_perpendicular_allowed0 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_perpendicular_outcome0_correct _rho _mu _sigma _lam _a _k _hd,
    WL_perpendicular_allowed0_correct _rho _mu _sigma _lam _a _ _hd]
  simp

theorem WL_perpendicular_equation1_locus (_rho _mu _sigma _lam z : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    Expr.eval (values _rho _mu _sigma z _k) WLCoincidenceArithmeticPerpendicular.n11 = 0 ↔
      WL_perpendicular_locus1 _rho _mu _sigma _lam _a _k := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [WLCoincidenceArithmeticPerpendicular.stratum2_root1_root3_q3_coincidence_operands_cell0 _rho _mu _sigma z _k (ne_of_gt hr) (ne_of_gt hs), sub_eq_zero,
    WL_perpendicular_pair1_roots _rho _mu _sigma _lam _k _hd,
    WL_perpendicular_locus1_correct _rho _mu _sigma _lam _a _k _hd]

theorem WL_perpendicular_status1_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_perpendicular_status1 _rho _mu _sigma _lam _a _k ↔
      ∃ w : Vec 3, WL_perpendicular_allowed1 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_perpendicular_status1_correct _rho _mu _sigma _lam _a _k _hd,
    WL_perpendicular_allowed1_correct _rho _mu _sigma _lam _a _ _hd]
  simp

theorem WL_perpendicular_outcome1_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_perpendicular_outcome1 _rho _mu _sigma _lam _a _k ↔
      ∃ w : Vec 3, WL_perpendicular_allowed1 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_perpendicular_outcome1_correct _rho _mu _sigma _lam _a _k _hd,
    WL_perpendicular_allowed1_correct _rho _mu _sigma _lam _a _ _hd]
  simp

theorem WL_perpendicular_equation2_locus (_rho _mu _sigma _lam z : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    Expr.eval (values _rho _mu _sigma z _k) WLCoincidenceArithmeticPerpendicular.n17 = 0 ↔
      WL_perpendicular_locus2 _rho _mu _sigma _lam _a _k := by
  have hr := _hd.1
  have hs := _hd.2.2.1
  rw [WLCoincidenceArithmeticPerpendicular.stratum2_root2_root3_q3_coincidence_operands_cell0 _rho _mu _sigma z _k (ne_of_gt hr) (ne_of_gt hs), sub_eq_zero,
    WL_perpendicular_pair2_roots _rho _mu _sigma _lam _k _hd,
    WL_perpendicular_locus2_correct _rho _mu _sigma _lam _a _k _hd]

theorem WL_perpendicular_status2_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_perpendicular_status2 _rho _mu _sigma _lam _a _k ↔
      ∃ w : Vec 3, WL_perpendicular_allowed2 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_perpendicular_status2_correct _rho _mu _sigma _lam _a _k _hd,
    WL_perpendicular_allowed2_correct _rho _mu _sigma _lam _a _ _hd]
  simp

theorem WL_perpendicular_outcome2_decides (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    WL_perpendicular_outcome2 _rho _mu _sigma _lam _a _k ↔
      ∃ w : Vec 3, WL_perpendicular_allowed2 _rho _mu _sigma _lam _a w := by
  simp_rw [WL_perpendicular_outcome2_correct _rho _mu _sigma _lam _a _k _hd,
    WL_perpendicular_allowed2_correct _rho _mu _sigma _lam _a _ _hd]
  simp

theorem generic_locus0_cross_engine (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_locus0 _rho _mu _sigma _lam _a _k ↔ WL_generic_locus0 _rho _mu _sigma _lam _a _k := by
  rw [PY_generic_locus0_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_locus0_correct _rho _mu _sigma _lam _a _k _hd]

theorem generic_allowed0_cross_engine (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_allowed0 _rho _mu _sigma _lam _a _k ↔ WL_generic_allowed0 _rho _mu _sigma _lam _a _k := by
  rw [PY_generic_allowed0_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_allowed0_correct _rho _mu _sigma _lam _a _k _hd]

theorem generic_status0_cross_engine (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_status0 _rho _mu _sigma _lam _a _k ↔ WL_generic_status0 _rho _mu _sigma _lam _a _k := by
  rw [PY_generic_status0_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_status0_correct _rho _mu _sigma _lam _a _k _hd]

theorem generic_locus1_cross_engine (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_locus1 _rho _mu _sigma _lam _a _k ↔ WL_generic_locus1 _rho _mu _sigma _lam _a _k := by
  rw [PY_generic_locus1_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_locus1_correct _rho _mu _sigma _lam _a _k _hd]

theorem generic_allowed1_cross_engine (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_allowed1 _rho _mu _sigma _lam _a _k ↔ WL_generic_allowed1 _rho _mu _sigma _lam _a _k := by
  rw [PY_generic_allowed1_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_allowed1_correct _rho _mu _sigma _lam _a _k _hd]

theorem generic_status1_cross_engine (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_status1 _rho _mu _sigma _lam _a _k ↔ WL_generic_status1 _rho _mu _sigma _lam _a _k := by
  rw [PY_generic_status1_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_status1_correct _rho _mu _sigma _lam _a _k _hd]

theorem generic_locus2_cross_engine (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_locus2 _rho _mu _sigma _lam _a _k ↔ WL_generic_locus2 _rho _mu _sigma _lam _a _k := by
  rw [PY_generic_locus2_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_locus2_correct _rho _mu _sigma _lam _a _k _hd]

theorem generic_allowed2_cross_engine (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_allowed2 _rho _mu _sigma _lam _a _k ↔ WL_generic_allowed2 _rho _mu _sigma _lam _a _k := by
  rw [PY_generic_allowed2_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_allowed2_correct _rho _mu _sigma _lam _a _k _hd]

theorem generic_status2_cross_engine (_rho _mu _sigma _lam : ℝ) (_a _k : Vec 3) (_hd : coincidenceDomain _rho _mu _sigma _lam) :
    PY_generic_status2 _rho _mu _sigma _lam _a _k ↔ WL_generic_status2 _rho _mu _sigma _lam _a _k := by
  rw [PY_generic_status2_correct _rho _mu _sigma _lam _a _k _hd,
    WL_generic_status2_correct _rho _mu _sigma _lam _a _k _hd]

end
end S10Audit.CAS.Coincidence
