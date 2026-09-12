import S10Audit.CAS.PYGenericCounts
import S10Audit.CAS.WLGenericCounts

namespace S10Audit.CAS
open S10Pilot
noncomputable section

namespace PYGenericCounts

theorem root0_rank (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n0 =
      (((PYGenericCountReference.root0Matrix rho mu sigma z k).rank : ℕ) : ℝ) := by
  rw [root1_n2_rank_cell0 rho mu sigma z k,
    PYGenericCountReference.root0_rank rho mu sigma z k _hr _hm _h]
  norm_num

theorem root0_nullity (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((matrixNullity (PYGenericCountReference.root0Matrix rho mu sigma z k) : ℕ) : ℝ) := by
  rw [root1_n2_nullity_cell0 rho mu sigma z k,
    PYGenericCountReference.root0_nullity rho mu sigma z k _hr _hm _h]
  norm_num

theorem root0_stacked_rank (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n2 =
      (((PYGenericCountReference.root0Stack rho mu sigma z k).rank : ℕ) : ℝ) := by
  rw [root1_n3_stacked_rank_cell0 rho mu sigma z k,
    PYGenericCountReference.root0_stacked_rank rho mu sigma z k _hr _hm _h]
  norm_num

theorem root0_transverse_nullity (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n3 =
      ((matrixNullity (PYGenericCountReference.root0Stack rho mu sigma z k) : ℕ) : ℝ) := by
  rw [root1_n3_transverse_nullity_cell0 rho mu sigma z k,
    PYGenericCountReference.root0_transverse_nullity rho mu sigma z k _hr _hm _h]
  norm_num

theorem root0_difference (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((nullityDifference (PYGenericCountReference.root0Matrix rho mu sigma z k) (PYGenericCountReference.root0Stack rho mu sigma z k) : ℤ) : ℝ) := by
  rw [root1_n4_nullity_difference_cell0 rho mu sigma z k,
    PYGenericCountReference.root0_difference rho mu sigma z k _hr _hm _h]
  norm_num

theorem root0_basis_count (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((basisCount (PYGenericCountReference.root0Basis rho mu sigma z k) : ℕ) : ℝ) := by
  rw [root1_n7_basis_count_cell0 rho mu sigma z k,
    PYGenericCountReference.root0_basis_count rho mu sigma z k]
  norm_num

theorem root0_basis_residual (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n3 =
      ((basisCountResidual (PYGenericCountReference.root0Basis rho mu sigma z k) (PYGenericCountReference.root0Matrix rho mu sigma z k) : ℤ) : ℝ) := by
  rw [root1_n7_basis_count_residual_cell0 rho mu sigma z k,
    PYGenericCountReference.root0_basis_residual rho mu sigma z k _hr _hm _h]
  norm_num

theorem root1_rank (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n0 =
      (((PYGenericCountReference.root1Matrix rho mu sigma z k).rank : ℕ) : ℝ) := by
  rw [root2_n2_rank_cell0 rho mu sigma z k,
    PYGenericCountReference.root1_rank rho mu sigma z k _hr _hm _h]
  norm_num

theorem root1_nullity (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((matrixNullity (PYGenericCountReference.root1Matrix rho mu sigma z k) : ℕ) : ℝ) := by
  rw [root2_n2_nullity_cell0 rho mu sigma z k,
    PYGenericCountReference.root1_nullity rho mu sigma z k _hr _hm _h]
  norm_num

theorem root1_stacked_rank (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n0 =
      (((PYGenericCountReference.root1Stack rho mu sigma z k).rank : ℕ) : ℝ) := by
  rw [root2_n3_stacked_rank_cell0 rho mu sigma z k,
    PYGenericCountReference.root1_stacked_rank rho mu sigma z k _hr _hm _h]
  norm_num

theorem root1_transverse_nullity (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((matrixNullity (PYGenericCountReference.root1Stack rho mu sigma z k) : ℕ) : ℝ) := by
  rw [root2_n3_transverse_nullity_cell0 rho mu sigma z k,
    PYGenericCountReference.root1_transverse_nullity rho mu sigma z k _hr _hm _h]
  norm_num

theorem root1_difference (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n3 =
      ((nullityDifference (PYGenericCountReference.root1Matrix rho mu sigma z k) (PYGenericCountReference.root1Stack rho mu sigma z k) : ℤ) : ℝ) := by
  rw [root2_n4_nullity_difference_cell0 rho mu sigma z k,
    PYGenericCountReference.root1_difference rho mu sigma z k _hr _hm _h]
  norm_num

theorem root1_basis_count (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((basisCount (PYGenericCountReference.root1Basis rho mu sigma z k) : ℕ) : ℝ) := by
  rw [root2_n7_basis_count_cell0 rho mu sigma z k,
    PYGenericCountReference.root1_basis_count rho mu sigma z k]
  norm_num

theorem root1_basis_residual (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n3 =
      ((basisCountResidual (PYGenericCountReference.root1Basis rho mu sigma z k) (PYGenericCountReference.root1Matrix rho mu sigma z k) : ℤ) : ℝ) := by
  rw [root2_n7_basis_count_residual_cell0 rho mu sigma z k,
    PYGenericCountReference.root1_basis_residual rho mu sigma z k _hr _hm _h]
  norm_num

theorem root2_rank (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n0 =
      (((PYGenericCountReference.root2Matrix rho mu sigma z k).rank : ℕ) : ℝ) := by
  rw [root3_n2_rank_cell0 rho mu sigma z k,
    PYGenericCountReference.root2_rank rho mu sigma z k _hr _hm _h]
  norm_num

theorem root2_nullity (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((matrixNullity (PYGenericCountReference.root2Matrix rho mu sigma z k) : ℕ) : ℝ) := by
  rw [root3_n2_nullity_cell0 rho mu sigma z k,
    PYGenericCountReference.root2_nullity rho mu sigma z k _hr _hm _h]
  norm_num

theorem root2_stacked_rank (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n2 =
      (((PYGenericCountReference.root2Stack rho mu sigma z k).rank : ℕ) : ℝ) := by
  rw [root3_n3_stacked_rank_cell0 rho mu sigma z k,
    PYGenericCountReference.root2_stacked_rank rho mu sigma z k _hr _hm _h]
  norm_num

theorem root2_transverse_nullity (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n3 =
      ((matrixNullity (PYGenericCountReference.root2Stack rho mu sigma z k) : ℕ) : ℝ) := by
  rw [root3_n3_transverse_nullity_cell0 rho mu sigma z k,
    PYGenericCountReference.root2_transverse_nullity rho mu sigma z k _hr _hm _h]
  norm_num

theorem root2_difference (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((nullityDifference (PYGenericCountReference.root2Matrix rho mu sigma z k) (PYGenericCountReference.root2Stack rho mu sigma z k) : ℤ) : ℝ) := by
  rw [root3_n4_nullity_difference_cell0 rho mu sigma z k,
    PYGenericCountReference.root2_difference rho mu sigma z k _hr _hm _h]
  norm_num

theorem root2_basis_count (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((basisCount (PYGenericCountReference.root2Basis rho mu sigma z k) : ℕ) : ℝ) := by
  rw [root3_n7_basis_count_cell0 rho mu sigma z k,
    PYGenericCountReference.root2_basis_count rho mu sigma z k]
  norm_num

theorem root2_basis_residual (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n3 =
      ((basisCountResidual (PYGenericCountReference.root2Basis rho mu sigma z k) (PYGenericCountReference.root2Matrix rho mu sigma z k) : ℤ) : ℝ) := by
  rw [root3_n7_basis_count_residual_cell0 rho mu sigma z k,
    PYGenericCountReference.root2_basis_residual rho mu sigma z k _hr _hm _h]
  norm_num

end PYGenericCounts

namespace WLGenericCounts

theorem root0_rank (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n0 =
      (((WLGenericCountReference.root0Matrix rho mu sigma z k).rank : ℕ) : ℝ) := by
  rw [root1_n2_rank_cell0 rho mu sigma z k,
    WLGenericCountReference.root0_rank rho mu sigma z k _hr _hm _h]
  norm_num

theorem root0_nullity (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((matrixNullity (WLGenericCountReference.root0Matrix rho mu sigma z k) : ℕ) : ℝ) := by
  rw [root1_n2_nullity_cell0 rho mu sigma z k,
    WLGenericCountReference.root0_nullity rho mu sigma z k _hr _hm _h]
  norm_num

theorem root0_stacked_rank (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n2 =
      (((WLGenericCountReference.root0Stack rho mu sigma z k).rank : ℕ) : ℝ) := by
  rw [root1_n3_stacked_rank_cell0 rho mu sigma z k,
    WLGenericCountReference.root0_stacked_rank rho mu sigma z k _hr _hm _h]
  norm_num

theorem root0_transverse_nullity (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n3 =
      ((matrixNullity (WLGenericCountReference.root0Stack rho mu sigma z k) : ℕ) : ℝ) := by
  rw [root1_n3_transverse_nullity_cell0 rho mu sigma z k,
    WLGenericCountReference.root0_transverse_nullity rho mu sigma z k _hr _hm _h]
  norm_num

theorem root0_difference (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((nullityDifference (WLGenericCountReference.root0Matrix rho mu sigma z k) (WLGenericCountReference.root0Stack rho mu sigma z k) : ℤ) : ℝ) := by
  rw [root1_n4_nullity_difference_cell0 rho mu sigma z k,
    WLGenericCountReference.root0_difference rho mu sigma z k _hr _hm _h]
  norm_num

theorem root0_basis_count (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((basisCount (WLGenericCountReference.root0Basis rho mu sigma z k) : ℕ) : ℝ) := by
  rw [root1_n7_basis_count_cell0 rho mu sigma z k,
    WLGenericCountReference.root0_basis_count rho mu sigma z k]
  norm_num

theorem root0_basis_residual (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n3 =
      ((basisCountResidual (WLGenericCountReference.root0Basis rho mu sigma z k) (WLGenericCountReference.root0Matrix rho mu sigma z k) : ℤ) : ℝ) := by
  rw [root1_n7_count_residual_cell0 rho mu sigma z k,
    WLGenericCountReference.root0_basis_residual rho mu sigma z k _hr _hm _h]
  norm_num

theorem root1_rank (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n0 =
      (((WLGenericCountReference.root1Matrix rho mu sigma z k).rank : ℕ) : ℝ) := by
  rw [root2_n2_rank_cell0 rho mu sigma z k,
    WLGenericCountReference.root1_rank rho mu sigma z k _hr _hm _h]
  norm_num

theorem root1_nullity (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((matrixNullity (WLGenericCountReference.root1Matrix rho mu sigma z k) : ℕ) : ℝ) := by
  rw [root2_n2_nullity_cell0 rho mu sigma z k,
    WLGenericCountReference.root1_nullity rho mu sigma z k _hr _hm _h]
  norm_num

theorem root1_stacked_rank (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n0 =
      (((WLGenericCountReference.root1Stack rho mu sigma z k).rank : ℕ) : ℝ) := by
  rw [root2_n3_stacked_rank_cell0 rho mu sigma z k,
    WLGenericCountReference.root1_stacked_rank rho mu sigma z k _hr _hm _h]
  norm_num

theorem root1_transverse_nullity (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((matrixNullity (WLGenericCountReference.root1Stack rho mu sigma z k) : ℕ) : ℝ) := by
  rw [root2_n3_transverse_nullity_cell0 rho mu sigma z k,
    WLGenericCountReference.root1_transverse_nullity rho mu sigma z k _hr _hm _h]
  norm_num

theorem root1_difference (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n3 =
      ((nullityDifference (WLGenericCountReference.root1Matrix rho mu sigma z k) (WLGenericCountReference.root1Stack rho mu sigma z k) : ℤ) : ℝ) := by
  rw [root2_n4_nullity_difference_cell0 rho mu sigma z k,
    WLGenericCountReference.root1_difference rho mu sigma z k _hr _hm _h]
  norm_num

theorem root1_basis_count (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((basisCount (WLGenericCountReference.root1Basis rho mu sigma z k) : ℕ) : ℝ) := by
  rw [root2_n7_basis_count_cell0 rho mu sigma z k,
    WLGenericCountReference.root1_basis_count rho mu sigma z k]
  norm_num

theorem root1_basis_residual (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n3 =
      ((basisCountResidual (WLGenericCountReference.root1Basis rho mu sigma z k) (WLGenericCountReference.root1Matrix rho mu sigma z k) : ℤ) : ℝ) := by
  rw [root2_n7_count_residual_cell0 rho mu sigma z k,
    WLGenericCountReference.root1_basis_residual rho mu sigma z k _hr _hm _h]
  norm_num

theorem root2_rank (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n0 =
      (((WLGenericCountReference.root2Matrix rho mu sigma z k).rank : ℕ) : ℝ) := by
  rw [root3_n2_rank_cell0 rho mu sigma z k,
    WLGenericCountReference.root2_rank rho mu sigma z k _hr _hm _h]
  norm_num

theorem root2_nullity (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((matrixNullity (WLGenericCountReference.root2Matrix rho mu sigma z k) : ℕ) : ℝ) := by
  rw [root3_n2_nullity_cell0 rho mu sigma z k,
    WLGenericCountReference.root2_nullity rho mu sigma z k _hr _hm _h]
  norm_num

theorem root2_stacked_rank (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n2 =
      (((WLGenericCountReference.root2Stack rho mu sigma z k).rank : ℕ) : ℝ) := by
  rw [root3_n3_stacked_rank_cell0 rho mu sigma z k,
    WLGenericCountReference.root2_stacked_rank rho mu sigma z k _hr _hm _h]
  norm_num

theorem root2_transverse_nullity (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n3 =
      ((matrixNullity (WLGenericCountReference.root2Stack rho mu sigma z k) : ℕ) : ℝ) := by
  rw [root3_n3_transverse_nullity_cell0 rho mu sigma z k,
    WLGenericCountReference.root2_transverse_nullity rho mu sigma z k _hr _hm _h]
  norm_num

theorem root2_difference (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((nullityDifference (WLGenericCountReference.root2Matrix rho mu sigma z k) (WLGenericCountReference.root2Stack rho mu sigma z k) : ℤ) : ℝ) := by
  rw [root3_n4_nullity_difference_cell0 rho mu sigma z k,
    WLGenericCountReference.root2_difference rho mu sigma z k _hr _hm _h]
  norm_num

theorem root2_basis_count (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n1 =
      ((basisCount (WLGenericCountReference.root2Basis rho mu sigma z k) : ℕ) : ℝ) := by
  rw [root3_n7_basis_count_cell0 rho mu sigma z k,
    WLGenericCountReference.root2_basis_count rho mu sigma z k]
  norm_num

theorem root2_basis_residual (rho mu sigma z : ℝ) (k : Vec 3)
    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :
    Expr.eval (values rho mu sigma z k) n3 =
      ((basisCountResidual (WLGenericCountReference.root2Basis rho mu sigma z k) (WLGenericCountReference.root2Matrix rho mu sigma z k) : ℤ) : ℝ) := by
  rw [root3_n7_count_residual_cell0 rho mu sigma z k,
    WLGenericCountReference.root2_basis_residual rho mu sigma z k _hr _hm _h]
  norm_num

end WLGenericCounts

end
end S10Audit.CAS
