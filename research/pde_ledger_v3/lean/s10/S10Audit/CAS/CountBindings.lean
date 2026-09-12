import S10Audit.CAS.PYCounts
import S10Audit.CAS.WLCounts

namespace S10Audit.CAS
open S10Pilot
noncomputable section

namespace PYCounts

theorem parallel0_rank (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n0 =
      (((PYRerun.parallel0Matrix rho mu sigma).rank : ℕ) : ℝ) := by
  rw [q8_stratum1_root1_n2_rank_cell0 rho mu sigma z k,
    PYCountReference.parallel0_rank rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem parallel0_nullity (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((matrixNullity (PYRerun.parallel0Matrix rho mu sigma) : ℕ) : ℝ) := by
  rw [q8_stratum1_root1_n2_nullity_cell0 rho mu sigma z k,
    PYCountReference.parallel0_nullity rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem parallel0_stacked_rank (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n2 =
      (((PYRerun.parallel0Stack rho mu sigma).rank : ℕ) : ℝ) := by
  rw [q8_stratum1_root1_n3_stacked_rank_cell0 rho mu sigma z k,
    PYCountReference.parallel0_stacked_rank rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem parallel0_transverse_nullity (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n3 =
      ((matrixNullity (PYRerun.parallel0Stack rho mu sigma) : ℕ) : ℝ) := by
  rw [q8_stratum1_root1_n3_transverse_nullity_cell0 rho mu sigma z k,
    PYCountReference.parallel0_transverse_nullity rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem parallel0_difference (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((nullityDifference (PYRerun.parallel0Matrix rho mu sigma) (PYRerun.parallel0Stack rho mu sigma) : ℤ) : ℝ) := by
  rw [q8_stratum1_root1_n4_nullity_difference_cell0 rho mu sigma z k,
    PYCountReference.parallel0_difference rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem parallel0_basis_count (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((basisCount (PYRerun.parallel0Basis rho mu sigma) : ℕ) : ℝ) := by
  rw [q8_stratum1_root1_n7_basis_count_cell0 rho mu sigma z k,
    PYCountReference.parallel0_basis_count rho mu sigma]
  norm_num

theorem parallel0_basis_residual (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n3 =
      ((basisCountResidual (PYRerun.parallel0Basis rho mu sigma) (PYRerun.parallel0Matrix rho mu sigma) : ℤ) : ℝ) := by
  rw [q8_stratum1_root1_n7_basis_count_residual_cell0 rho mu sigma z k,
    PYCountReference.parallel0_basis_residual rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem parallel1_rank (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n1 =
      (((PYRerun.parallel1Matrix rho mu sigma).rank : ℕ) : ℝ) := by
  rw [q8_stratum1_root2_n2_rank_cell0 rho mu sigma z k,
    PYCountReference.parallel1_rank rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem parallel1_nullity (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n0 =
      ((matrixNullity (PYRerun.parallel1Matrix rho mu sigma) : ℕ) : ℝ) := by
  rw [q8_stratum1_root2_n2_nullity_cell0 rho mu sigma z k,
    PYCountReference.parallel1_nullity rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem parallel1_stacked_rank (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n1 =
      (((PYRerun.parallel1Stack rho mu sigma).rank : ℕ) : ℝ) := by
  rw [q8_stratum1_root2_n3_stacked_rank_cell0 rho mu sigma z k,
    PYCountReference.parallel1_stacked_rank rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem parallel1_transverse_nullity (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n0 =
      ((matrixNullity (PYRerun.parallel1Stack rho mu sigma) : ℕ) : ℝ) := by
  rw [q8_stratum1_root2_n3_transverse_nullity_cell0 rho mu sigma z k,
    PYCountReference.parallel1_transverse_nullity rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem parallel1_difference (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n3 =
      ((nullityDifference (PYRerun.parallel1Matrix rho mu sigma) (PYRerun.parallel1Stack rho mu sigma) : ℤ) : ℝ) := by
  rw [q8_stratum1_root2_n4_nullity_difference_cell0 rho mu sigma z k,
    PYCountReference.parallel1_difference rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem parallel1_basis_count (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n0 =
      ((basisCount (PYRerun.parallel1Basis rho mu sigma) : ℕ) : ℝ) := by
  rw [q8_stratum1_root2_n7_basis_count_cell0 rho mu sigma z k,
    PYCountReference.parallel1_basis_count rho mu sigma]
  norm_num

theorem parallel1_basis_residual (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n3 =
      ((basisCountResidual (PYRerun.parallel1Basis rho mu sigma) (PYRerun.parallel1Matrix rho mu sigma) : ℤ) : ℝ) := by
  rw [q8_stratum1_root2_n7_basis_count_residual_cell0 rho mu sigma z k,
    PYCountReference.parallel1_basis_residual rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular0_rank (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n0 =
      (((PYRerun.perpendicular0Matrix rho mu sigma).rank : ℕ) : ℝ) := by
  rw [q8_stratum2_root1_n2_rank_cell0 rho mu sigma z k,
    PYCountReference.perpendicular0_rank rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular0_nullity (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((matrixNullity (PYRerun.perpendicular0Matrix rho mu sigma) : ℕ) : ℝ) := by
  rw [q8_stratum2_root1_n2_nullity_cell0 rho mu sigma z k,
    PYCountReference.perpendicular0_nullity rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular0_stacked_rank (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n2 =
      (((PYRerun.perpendicular0Stack rho mu sigma).rank : ℕ) : ℝ) := by
  rw [q8_stratum2_root1_n3_stacked_rank_cell0 rho mu sigma z k,
    PYCountReference.perpendicular0_stacked_rank rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular0_transverse_nullity (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n3 =
      ((matrixNullity (PYRerun.perpendicular0Stack rho mu sigma) : ℕ) : ℝ) := by
  rw [q8_stratum2_root1_n3_transverse_nullity_cell0 rho mu sigma z k,
    PYCountReference.perpendicular0_transverse_nullity rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular0_difference (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((nullityDifference (PYRerun.perpendicular0Matrix rho mu sigma) (PYRerun.perpendicular0Stack rho mu sigma) : ℤ) : ℝ) := by
  rw [q8_stratum2_root1_n4_nullity_difference_cell0 rho mu sigma z k,
    PYCountReference.perpendicular0_difference rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular0_basis_count (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((basisCount (PYRerun.perpendicular0Basis rho mu sigma) : ℕ) : ℝ) := by
  rw [q8_stratum2_root1_n7_basis_count_cell0 rho mu sigma z k,
    PYCountReference.perpendicular0_basis_count rho mu sigma]
  norm_num

theorem perpendicular0_basis_residual (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n3 =
      ((basisCountResidual (PYRerun.perpendicular0Basis rho mu sigma) (PYRerun.perpendicular0Matrix rho mu sigma) : ℤ) : ℝ) := by
  rw [q8_stratum2_root1_n7_basis_count_residual_cell0 rho mu sigma z k,
    PYCountReference.perpendicular0_basis_residual rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular1_rank (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n0 =
      (((PYRerun.perpendicular1Matrix rho mu sigma).rank : ℕ) : ℝ) := by
  rw [q8_stratum2_root2_n2_rank_cell0 rho mu sigma z k,
    PYCountReference.perpendicular1_rank rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular1_nullity (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((matrixNullity (PYRerun.perpendicular1Matrix rho mu sigma) : ℕ) : ℝ) := by
  rw [q8_stratum2_root2_n2_nullity_cell0 rho mu sigma z k,
    PYCountReference.perpendicular1_nullity rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular1_stacked_rank (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n0 =
      (((PYRerun.perpendicular1Stack rho mu sigma).rank : ℕ) : ℝ) := by
  rw [q8_stratum2_root2_n3_stacked_rank_cell0 rho mu sigma z k,
    PYCountReference.perpendicular1_stacked_rank rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular1_transverse_nullity (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((matrixNullity (PYRerun.perpendicular1Stack rho mu sigma) : ℕ) : ℝ) := by
  rw [q8_stratum2_root2_n3_transverse_nullity_cell0 rho mu sigma z k,
    PYCountReference.perpendicular1_transverse_nullity rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular1_difference (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n3 =
      ((nullityDifference (PYRerun.perpendicular1Matrix rho mu sigma) (PYRerun.perpendicular1Stack rho mu sigma) : ℤ) : ℝ) := by
  rw [q8_stratum2_root2_n4_nullity_difference_cell0 rho mu sigma z k,
    PYCountReference.perpendicular1_difference rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular1_basis_count (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((basisCount (PYRerun.perpendicular1Basis rho mu sigma) : ℕ) : ℝ) := by
  rw [q8_stratum2_root2_n7_basis_count_cell0 rho mu sigma z k,
    PYCountReference.perpendicular1_basis_count rho mu sigma]
  norm_num

theorem perpendicular1_basis_residual (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n3 =
      ((basisCountResidual (PYRerun.perpendicular1Basis rho mu sigma) (PYRerun.perpendicular1Matrix rho mu sigma) : ℤ) : ℝ) := by
  rw [q8_stratum2_root2_n7_basis_count_residual_cell0 rho mu sigma z k,
    PYCountReference.perpendicular1_basis_residual rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular2_rank (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n0 =
      (((PYRerun.perpendicular2Matrix rho mu sigma).rank : ℕ) : ℝ) := by
  rw [q8_stratum2_root3_n2_rank_cell0 rho mu sigma z k,
    PYCountReference.perpendicular2_rank rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular2_nullity (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((matrixNullity (PYRerun.perpendicular2Matrix rho mu sigma) : ℕ) : ℝ) := by
  rw [q8_stratum2_root3_n2_nullity_cell0 rho mu sigma z k,
    PYCountReference.perpendicular2_nullity rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular2_stacked_rank (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n0 =
      (((PYRerun.perpendicular2Stack rho mu sigma).rank : ℕ) : ℝ) := by
  rw [q8_stratum2_root3_n3_stacked_rank_cell0 rho mu sigma z k,
    PYCountReference.perpendicular2_stacked_rank rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular2_transverse_nullity (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((matrixNullity (PYRerun.perpendicular2Stack rho mu sigma) : ℕ) : ℝ) := by
  rw [q8_stratum2_root3_n3_transverse_nullity_cell0 rho mu sigma z k,
    PYCountReference.perpendicular2_transverse_nullity rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular2_difference (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n3 =
      ((nullityDifference (PYRerun.perpendicular2Matrix rho mu sigma) (PYRerun.perpendicular2Stack rho mu sigma) : ℤ) : ℝ) := by
  rw [q8_stratum2_root3_n4_nullity_difference_cell0 rho mu sigma z k,
    PYCountReference.perpendicular2_difference rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular2_basis_count (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((basisCount (PYRerun.perpendicular2Basis rho mu sigma) : ℕ) : ℝ) := by
  rw [q8_stratum2_root3_n7_basis_count_cell0 rho mu sigma z k,
    PYCountReference.perpendicular2_basis_count rho mu sigma]
  norm_num

theorem perpendicular2_basis_residual (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n3 =
      ((basisCountResidual (PYRerun.perpendicular2Basis rho mu sigma) (PYRerun.perpendicular2Matrix rho mu sigma) : ℤ) : ℝ) := by
  rw [q8_stratum2_root3_n7_basis_count_residual_cell0 rho mu sigma z k,
    PYCountReference.perpendicular2_basis_residual rho mu sigma _hr _hm _hs _hs1]
  norm_num

end PYCounts

namespace WLCounts

theorem parallel0_rank (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n0 =
      (((WLRerun.parallel0Matrix rho mu sigma).rank : ℕ) : ℝ) := by
  rw [stratum1_root1_n2_rank_cell0 rho mu sigma z k,
    WLCountReference.parallel0_rank rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem parallel0_nullity (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((matrixNullity (WLRerun.parallel0Matrix rho mu sigma) : ℕ) : ℝ) := by
  rw [stratum1_root1_n2_nullity_cell0 rho mu sigma z k,
    WLCountReference.parallel0_nullity rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem parallel0_stacked_rank (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n2 =
      (((WLRerun.parallel0Stack rho mu sigma).rank : ℕ) : ℝ) := by
  rw [stratum1_root1_n3_stacked_rank_cell0 rho mu sigma z k,
    WLCountReference.parallel0_stacked_rank rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem parallel0_transverse_nullity (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n3 =
      ((matrixNullity (WLRerun.parallel0Stack rho mu sigma) : ℕ) : ℝ) := by
  rw [stratum1_root1_n3_transverse_nullity_cell0 rho mu sigma z k,
    WLCountReference.parallel0_transverse_nullity rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem parallel0_difference (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((nullityDifference (WLRerun.parallel0Matrix rho mu sigma) (WLRerun.parallel0Stack rho mu sigma) : ℤ) : ℝ) := by
  rw [stratum1_root1_n4_nullity_difference_cell0 rho mu sigma z k,
    WLCountReference.parallel0_difference rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem parallel0_basis_count (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((basisCount (WLRerun.parallel0Basis rho mu sigma) : ℕ) : ℝ) := by
  rw [stratum1_root1_n7_basis_count_cell0 rho mu sigma z k,
    WLCountReference.parallel0_basis_count rho mu sigma]
  norm_num

theorem parallel0_basis_residual (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n3 =
      ((basisCountResidual (WLRerun.parallel0Basis rho mu sigma) (WLRerun.parallel0Matrix rho mu sigma) : ℤ) : ℝ) := by
  rw [stratum1_root1_n7_count_residual_cell0 rho mu sigma z k,
    WLCountReference.parallel0_basis_residual rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem parallel1_rank (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n1 =
      (((WLRerun.parallel1Matrix rho mu sigma).rank : ℕ) : ℝ) := by
  rw [stratum1_root2_n2_rank_cell0 rho mu sigma z k,
    WLCountReference.parallel1_rank rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem parallel1_nullity (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n0 =
      ((matrixNullity (WLRerun.parallel1Matrix rho mu sigma) : ℕ) : ℝ) := by
  rw [stratum1_root2_n2_nullity_cell0 rho mu sigma z k,
    WLCountReference.parallel1_nullity rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem parallel1_stacked_rank (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n1 =
      (((WLRerun.parallel1Stack rho mu sigma).rank : ℕ) : ℝ) := by
  rw [stratum1_root2_n3_stacked_rank_cell0 rho mu sigma z k,
    WLCountReference.parallel1_stacked_rank rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem parallel1_transverse_nullity (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n0 =
      ((matrixNullity (WLRerun.parallel1Stack rho mu sigma) : ℕ) : ℝ) := by
  rw [stratum1_root2_n3_transverse_nullity_cell0 rho mu sigma z k,
    WLCountReference.parallel1_transverse_nullity rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem parallel1_difference (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n3 =
      ((nullityDifference (WLRerun.parallel1Matrix rho mu sigma) (WLRerun.parallel1Stack rho mu sigma) : ℤ) : ℝ) := by
  rw [stratum1_root2_n4_nullity_difference_cell0 rho mu sigma z k,
    WLCountReference.parallel1_difference rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem parallel1_basis_count (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n0 =
      ((basisCount (WLRerun.parallel1Basis rho mu sigma) : ℕ) : ℝ) := by
  rw [stratum1_root2_n7_basis_count_cell0 rho mu sigma z k,
    WLCountReference.parallel1_basis_count rho mu sigma]
  norm_num

theorem parallel1_basis_residual (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n3 =
      ((basisCountResidual (WLRerun.parallel1Basis rho mu sigma) (WLRerun.parallel1Matrix rho mu sigma) : ℤ) : ℝ) := by
  rw [stratum1_root2_n7_count_residual_cell0 rho mu sigma z k,
    WLCountReference.parallel1_basis_residual rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular0_rank (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n0 =
      (((WLRerun.perpendicular0Matrix rho mu sigma).rank : ℕ) : ℝ) := by
  rw [stratum2_root1_n2_rank_cell0 rho mu sigma z k,
    WLCountReference.perpendicular0_rank rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular0_nullity (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((matrixNullity (WLRerun.perpendicular0Matrix rho mu sigma) : ℕ) : ℝ) := by
  rw [stratum2_root1_n2_nullity_cell0 rho mu sigma z k,
    WLCountReference.perpendicular0_nullity rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular0_stacked_rank (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n2 =
      (((WLRerun.perpendicular0Stack rho mu sigma).rank : ℕ) : ℝ) := by
  rw [stratum2_root1_n3_stacked_rank_cell0 rho mu sigma z k,
    WLCountReference.perpendicular0_stacked_rank rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular0_transverse_nullity (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n3 =
      ((matrixNullity (WLRerun.perpendicular0Stack rho mu sigma) : ℕ) : ℝ) := by
  rw [stratum2_root1_n3_transverse_nullity_cell0 rho mu sigma z k,
    WLCountReference.perpendicular0_transverse_nullity rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular0_difference (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((nullityDifference (WLRerun.perpendicular0Matrix rho mu sigma) (WLRerun.perpendicular0Stack rho mu sigma) : ℤ) : ℝ) := by
  rw [stratum2_root1_n4_nullity_difference_cell0 rho mu sigma z k,
    WLCountReference.perpendicular0_difference rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular0_basis_count (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((basisCount (WLRerun.perpendicular0Basis rho mu sigma) : ℕ) : ℝ) := by
  rw [stratum2_root1_n7_basis_count_cell0 rho mu sigma z k,
    WLCountReference.perpendicular0_basis_count rho mu sigma]
  norm_num

theorem perpendicular0_basis_residual (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n3 =
      ((basisCountResidual (WLRerun.perpendicular0Basis rho mu sigma) (WLRerun.perpendicular0Matrix rho mu sigma) : ℤ) : ℝ) := by
  rw [stratum2_root1_n7_count_residual_cell0 rho mu sigma z k,
    WLCountReference.perpendicular0_basis_residual rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular1_rank (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n0 =
      (((WLRerun.perpendicular1Matrix rho mu sigma).rank : ℕ) : ℝ) := by
  rw [stratum2_root2_n2_rank_cell0 rho mu sigma z k,
    WLCountReference.perpendicular1_rank rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular1_nullity (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((matrixNullity (WLRerun.perpendicular1Matrix rho mu sigma) : ℕ) : ℝ) := by
  rw [stratum2_root2_n2_nullity_cell0 rho mu sigma z k,
    WLCountReference.perpendicular1_nullity rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular1_stacked_rank (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n0 =
      (((WLRerun.perpendicular1Stack rho mu sigma).rank : ℕ) : ℝ) := by
  rw [stratum2_root2_n3_stacked_rank_cell0 rho mu sigma z k,
    WLCountReference.perpendicular1_stacked_rank rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular1_transverse_nullity (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((matrixNullity (WLRerun.perpendicular1Stack rho mu sigma) : ℕ) : ℝ) := by
  rw [stratum2_root2_n3_transverse_nullity_cell0 rho mu sigma z k,
    WLCountReference.perpendicular1_transverse_nullity rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular1_difference (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n3 =
      ((nullityDifference (WLRerun.perpendicular1Matrix rho mu sigma) (WLRerun.perpendicular1Stack rho mu sigma) : ℤ) : ℝ) := by
  rw [stratum2_root2_n4_nullity_difference_cell0 rho mu sigma z k,
    WLCountReference.perpendicular1_difference rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular1_basis_count (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((basisCount (WLRerun.perpendicular1Basis rho mu sigma) : ℕ) : ℝ) := by
  rw [stratum2_root2_n7_basis_count_cell0 rho mu sigma z k,
    WLCountReference.perpendicular1_basis_count rho mu sigma]
  norm_num

theorem perpendicular1_basis_residual (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n3 =
      ((basisCountResidual (WLRerun.perpendicular1Basis rho mu sigma) (WLRerun.perpendicular1Matrix rho mu sigma) : ℤ) : ℝ) := by
  rw [stratum2_root2_n7_count_residual_cell0 rho mu sigma z k,
    WLCountReference.perpendicular1_basis_residual rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular2_rank (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n0 =
      (((WLRerun.perpendicular2Matrix rho mu sigma).rank : ℕ) : ℝ) := by
  rw [stratum2_root3_n2_rank_cell0 rho mu sigma z k,
    WLCountReference.perpendicular2_rank rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular2_nullity (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((matrixNullity (WLRerun.perpendicular2Matrix rho mu sigma) : ℕ) : ℝ) := by
  rw [stratum2_root3_n2_nullity_cell0 rho mu sigma z k,
    WLCountReference.perpendicular2_nullity rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular2_stacked_rank (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n0 =
      (((WLRerun.perpendicular2Stack rho mu sigma).rank : ℕ) : ℝ) := by
  rw [stratum2_root3_n3_stacked_rank_cell0 rho mu sigma z k,
    WLCountReference.perpendicular2_stacked_rank rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular2_transverse_nullity (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((matrixNullity (WLRerun.perpendicular2Stack rho mu sigma) : ℕ) : ℝ) := by
  rw [stratum2_root3_n3_transverse_nullity_cell0 rho mu sigma z k,
    WLCountReference.perpendicular2_transverse_nullity rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular2_difference (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n3 =
      ((nullityDifference (WLRerun.perpendicular2Matrix rho mu sigma) (WLRerun.perpendicular2Stack rho mu sigma) : ℤ) : ℝ) := by
  rw [stratum2_root3_n4_nullity_difference_cell0 rho mu sigma z k,
    WLCountReference.perpendicular2_difference rho mu sigma _hr _hm _hs _hs1]
  norm_num

theorem perpendicular2_basis_count (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((basisCount (WLRerun.perpendicular2Basis rho mu sigma) : ℕ) : ℝ) := by
  rw [stratum2_root3_n7_basis_count_cell0 rho mu sigma z k,
    WLCountReference.perpendicular2_basis_count rho mu sigma]
  norm_num

theorem perpendicular2_basis_residual (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma z k) n3 =
      ((basisCountResidual (WLRerun.perpendicular2Basis rho mu sigma) (WLRerun.perpendicular2Matrix rho mu sigma) : ℤ) : ℝ) := by
  rw [stratum2_root3_n7_count_residual_cell0 rho mu sigma z k,
    WLCountReference.perpendicular2_basis_residual rho mu sigma _hr _hm _hs _hs1]
  norm_num

end WLCounts

end
end S10Audit.CAS
