import S10Audit.CAS.CountSupport

set_option backward.isDefEq.respectTransparency false

namespace S10Audit.CAS
open S10Pilot
noncomputable section

namespace PYCountReference

theorem parallel0_kernel (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    LinearMap.ker (PYRerun.parallel0Matrix rho mu sigma).mulVecLin = rootMode sigma PYLoci.parallelPoint 0 := by
  ext a
  change (PYRerun.parallel0Matrix rho mu sigma).mulVec a = 0 ↔ a ∈ rootMode sigma PYLoci.parallelPoint 0
  rw [PYRerun.parallel0_basis_complete rho mu sigma a hr hm hs hs1,
    PYRerun.parallel0Basis_reference rho mu sigma hr (ne_of_gt hs),
    PYReference.parallel0_complete sigma hs hs1]

theorem parallel0_stack_kernel (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    LinearMap.ker (PYRerun.parallel0Stack rho mu sigma).mulVecLin = rootMode sigma PYLoci.parallelPoint 0 ⊓ transverseSpace PYLoci.parallelPoint := by
  have he : (PYRerun.parallel0Stack rho mu sigma) = appendConstraint (PYRerun.parallel0Matrix rho mu sigma) PYLoci.parallelPoint := by
    rw [PYRerun.parallel0Stack_reference rho mu sigma hr (ne_of_gt hs),
      PYRerun.parallel0Matrix_reference rho mu sigma hr (ne_of_gt hs)]
    rfl
  rw [he, appendConstraint_kernel, parallel0_kernel rho mu sigma hr hm hs hs1]

theorem parallel0_nullity (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : matrixNullity (PYRerun.parallel0Matrix rho mu sigma) = 1 := by
  rw [matrixNullity, parallel0_kernel rho mu sigma hr hm hs hs1]
  exact (PYReference.parallel0_dimensions sigma hs hs1).1

theorem parallel0_transverse_nullity (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : matrixNullity (PYRerun.parallel0Stack rho mu sigma) = 0 := by
  rw [matrixNullity, parallel0_stack_kernel rho mu sigma hr hm hs hs1]
  exact (PYReference.parallel0_dimensions sigma hs hs1).2

theorem parallel0_rank (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : (PYRerun.parallel0Matrix rho mu sigma).rank = 2 := by
  exact matrix_rank_of_nullity (PYRerun.parallel0Matrix rho mu sigma) (parallel0_nullity rho mu sigma hr hm hs hs1)

theorem parallel0_stacked_rank (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : (PYRerun.parallel0Stack rho mu sigma).rank = 3 := by
  exact matrix_rank_of_nullity (PYRerun.parallel0Stack rho mu sigma) (parallel0_transverse_nullity rho mu sigma hr hm hs hs1)

theorem parallel0_difference (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : nullityDifference (PYRerun.parallel0Matrix rho mu sigma) (PYRerun.parallel0Stack rho mu sigma) = (1 : ℤ) := by
  rw [nullityDifference, parallel0_nullity rho mu sigma hr hm hs hs1, parallel0_transverse_nullity rho mu sigma hr hm hs hs1]
  norm_num

theorem parallel0_basis_count (rho mu sigma : ℝ) : basisCount (PYRerun.parallel0Basis rho mu sigma) = (1 : ℕ) := by
  rfl

theorem parallel0_basis_residual (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : basisCountResidual (PYRerun.parallel0Basis rho mu sigma) (PYRerun.parallel0Matrix rho mu sigma) = (0 : ℤ) := by
  rw [basisCountResidual, parallel0_basis_count, parallel0_nullity rho mu sigma hr hm hs hs1]
  norm_num

theorem parallel0_physical_counts (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) (c : ℝ) (hc : c ≠ 0) :
    matrixNullity ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • PYLoci.parallelPoint) 0) (c • PYLoci.parallelPoint)) = 1 ∧ ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • PYLoci.parallelPoint) 0) (c • PYLoci.parallelPoint)).rank = 2 ∧
      matrixNullity (appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • PYLoci.parallelPoint) 0) (c • PYLoci.parallelPoint)) (c • PYLoci.parallelPoint)) = 0 ∧
      (appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • PYLoci.parallelPoint) 0) (c • PYLoci.parallelPoint)) (c • PYLoci.parallelPoint)).rank = 3 := by
  have hM := matrix_counts_of_kernel _ _
    (physical_modal_kernel rho mu sigma (referenceRoot rho mu sigma PYLoci.parallelPoint 0) c PYLoci.parallelPoint hc)
  have hS := matrix_counts_of_kernel _ _
    (physical_constraint_kernel rho mu sigma (referenceRoot rho mu sigma PYLoci.parallelPoint 0) c PYLoci.parallelPoint hc)
  have hB : ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma PYLoci.parallelPoint 0) PYLoci.parallelPoint) = (PYRerun.parallel0Matrix rho mu sigma) := (PYRerun.parallel0Matrix_reference rho mu sigma hr (ne_of_gt hs)).symm
  have hA : appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma PYLoci.parallelPoint 0) PYLoci.parallelPoint) PYLoci.parallelPoint = (PYRerun.parallel0Stack rho mu sigma) :=
    (PYRerun.parallel0Stack_reference rho mu sigma hr (ne_of_gt hs)).symm
  rw [hB] at hM
  rw [hA] at hS
  rw [referenceRoot_scale]
  exact ⟨hM.1.trans (parallel0_nullity rho mu sigma hr hm hs hs1), hM.2.trans (parallel0_rank rho mu sigma hr hm hs hs1),
    hS.1.trans (parallel0_transverse_nullity rho mu sigma hr hm hs hs1), hS.2.trans (parallel0_stacked_rank rho mu sigma hr hm hs hs1)⟩

theorem parallel1_kernel (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    LinearMap.ker (PYRerun.parallel1Matrix rho mu sigma).mulVecLin = rootMode sigma PYLoci.parallelPoint 1 := by
  ext a
  change (PYRerun.parallel1Matrix rho mu sigma).mulVec a = 0 ↔ a ∈ rootMode sigma PYLoci.parallelPoint 1
  rw [PYRerun.parallel1_basis_complete rho mu sigma a hr hm hs hs1,
    PYRerun.parallel1Basis_reference rho mu sigma hr (ne_of_gt hs),
    PYReference.parallel1_complete sigma hs hs1]

theorem parallel1_stack_kernel (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    LinearMap.ker (PYRerun.parallel1Stack rho mu sigma).mulVecLin = rootMode sigma PYLoci.parallelPoint 1 ⊓ transverseSpace PYLoci.parallelPoint := by
  have he : (PYRerun.parallel1Stack rho mu sigma) = appendConstraint (PYRerun.parallel1Matrix rho mu sigma) PYLoci.parallelPoint := by
    rw [PYRerun.parallel1Stack_reference rho mu sigma hr (ne_of_gt hs),
      PYRerun.parallel1Matrix_reference rho mu sigma hr (ne_of_gt hs)]
    rfl
  rw [he, appendConstraint_kernel, parallel1_kernel rho mu sigma hr hm hs hs1]

theorem parallel1_nullity (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : matrixNullity (PYRerun.parallel1Matrix rho mu sigma) = 2 := by
  rw [matrixNullity, parallel1_kernel rho mu sigma hr hm hs hs1]
  exact (PYReference.parallel1_dimensions sigma hs hs1).1

theorem parallel1_transverse_nullity (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : matrixNullity (PYRerun.parallel1Stack rho mu sigma) = 2 := by
  rw [matrixNullity, parallel1_stack_kernel rho mu sigma hr hm hs hs1]
  exact (PYReference.parallel1_dimensions sigma hs hs1).2

theorem parallel1_rank (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : (PYRerun.parallel1Matrix rho mu sigma).rank = 1 := by
  exact matrix_rank_of_nullity (PYRerun.parallel1Matrix rho mu sigma) (parallel1_nullity rho mu sigma hr hm hs hs1)

theorem parallel1_stacked_rank (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : (PYRerun.parallel1Stack rho mu sigma).rank = 1 := by
  exact matrix_rank_of_nullity (PYRerun.parallel1Stack rho mu sigma) (parallel1_transverse_nullity rho mu sigma hr hm hs hs1)

theorem parallel1_difference (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : nullityDifference (PYRerun.parallel1Matrix rho mu sigma) (PYRerun.parallel1Stack rho mu sigma) = (0 : ℤ) := by
  rw [nullityDifference, parallel1_nullity rho mu sigma hr hm hs hs1, parallel1_transverse_nullity rho mu sigma hr hm hs hs1]
  norm_num

theorem parallel1_basis_count (rho mu sigma : ℝ) : basisCount (PYRerun.parallel1Basis rho mu sigma) = (2 : ℕ) := by
  rfl

theorem parallel1_basis_residual (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : basisCountResidual (PYRerun.parallel1Basis rho mu sigma) (PYRerun.parallel1Matrix rho mu sigma) = (0 : ℤ) := by
  rw [basisCountResidual, parallel1_basis_count, parallel1_nullity rho mu sigma hr hm hs hs1]
  norm_num

theorem parallel1_physical_counts (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) (c : ℝ) (hc : c ≠ 0) :
    matrixNullity ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • PYLoci.parallelPoint) 1) (c • PYLoci.parallelPoint)) = 2 ∧ ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • PYLoci.parallelPoint) 1) (c • PYLoci.parallelPoint)).rank = 1 ∧
      matrixNullity (appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • PYLoci.parallelPoint) 1) (c • PYLoci.parallelPoint)) (c • PYLoci.parallelPoint)) = 2 ∧
      (appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • PYLoci.parallelPoint) 1) (c • PYLoci.parallelPoint)) (c • PYLoci.parallelPoint)).rank = 1 := by
  have hM := matrix_counts_of_kernel _ _
    (physical_modal_kernel rho mu sigma (referenceRoot rho mu sigma PYLoci.parallelPoint 1) c PYLoci.parallelPoint hc)
  have hS := matrix_counts_of_kernel _ _
    (physical_constraint_kernel rho mu sigma (referenceRoot rho mu sigma PYLoci.parallelPoint 1) c PYLoci.parallelPoint hc)
  have hB : ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma PYLoci.parallelPoint 1) PYLoci.parallelPoint) = (PYRerun.parallel1Matrix rho mu sigma) := (PYRerun.parallel1Matrix_reference rho mu sigma hr (ne_of_gt hs)).symm
  have hA : appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma PYLoci.parallelPoint 1) PYLoci.parallelPoint) PYLoci.parallelPoint = (PYRerun.parallel1Stack rho mu sigma) :=
    (PYRerun.parallel1Stack_reference rho mu sigma hr (ne_of_gt hs)).symm
  rw [hB] at hM
  rw [hA] at hS
  rw [referenceRoot_scale]
  exact ⟨hM.1.trans (parallel1_nullity rho mu sigma hr hm hs hs1), hM.2.trans (parallel1_rank rho mu sigma hr hm hs hs1),
    hS.1.trans (parallel1_transverse_nullity rho mu sigma hr hm hs hs1), hS.2.trans (parallel1_stacked_rank rho mu sigma hr hm hs hs1)⟩

theorem perpendicular0_kernel (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    LinearMap.ker (PYRerun.perpendicular0Matrix rho mu sigma).mulVecLin = rootMode sigma PYLoci.perpendicularPoint 0 := by
  ext a
  change (PYRerun.perpendicular0Matrix rho mu sigma).mulVec a = 0 ↔ a ∈ rootMode sigma PYLoci.perpendicularPoint 0
  rw [PYRerun.perpendicular0_basis_complete rho mu sigma a hr hm hs hs1,
    PYRerun.perpendicular0Basis_reference rho mu sigma hr (ne_of_gt hs),
    PYReference.perpendicular0_complete sigma hs hs1]

theorem perpendicular0_stack_kernel (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    LinearMap.ker (PYRerun.perpendicular0Stack rho mu sigma).mulVecLin = rootMode sigma PYLoci.perpendicularPoint 0 ⊓ transverseSpace PYLoci.perpendicularPoint := by
  have he : (PYRerun.perpendicular0Stack rho mu sigma) = appendConstraint (PYRerun.perpendicular0Matrix rho mu sigma) PYLoci.perpendicularPoint := by
    rw [PYRerun.perpendicular0Stack_reference rho mu sigma hr (ne_of_gt hs),
      PYRerun.perpendicular0Matrix_reference rho mu sigma hr (ne_of_gt hs)]
    rfl
  rw [he, appendConstraint_kernel, perpendicular0_kernel rho mu sigma hr hm hs hs1]

theorem perpendicular0_nullity (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : matrixNullity (PYRerun.perpendicular0Matrix rho mu sigma) = 1 := by
  rw [matrixNullity, perpendicular0_kernel rho mu sigma hr hm hs hs1]
  exact (PYReference.perpendicular0_dimensions sigma hs hs1).1

theorem perpendicular0_transverse_nullity (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : matrixNullity (PYRerun.perpendicular0Stack rho mu sigma) = 0 := by
  rw [matrixNullity, perpendicular0_stack_kernel rho mu sigma hr hm hs hs1]
  exact (PYReference.perpendicular0_dimensions sigma hs hs1).2

theorem perpendicular0_rank (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : (PYRerun.perpendicular0Matrix rho mu sigma).rank = 2 := by
  exact matrix_rank_of_nullity (PYRerun.perpendicular0Matrix rho mu sigma) (perpendicular0_nullity rho mu sigma hr hm hs hs1)

theorem perpendicular0_stacked_rank (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : (PYRerun.perpendicular0Stack rho mu sigma).rank = 3 := by
  exact matrix_rank_of_nullity (PYRerun.perpendicular0Stack rho mu sigma) (perpendicular0_transverse_nullity rho mu sigma hr hm hs hs1)

theorem perpendicular0_difference (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : nullityDifference (PYRerun.perpendicular0Matrix rho mu sigma) (PYRerun.perpendicular0Stack rho mu sigma) = (1 : ℤ) := by
  rw [nullityDifference, perpendicular0_nullity rho mu sigma hr hm hs hs1, perpendicular0_transverse_nullity rho mu sigma hr hm hs hs1]
  norm_num

theorem perpendicular0_basis_count (rho mu sigma : ℝ) : basisCount (PYRerun.perpendicular0Basis rho mu sigma) = (1 : ℕ) := by
  rfl

theorem perpendicular0_basis_residual (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : basisCountResidual (PYRerun.perpendicular0Basis rho mu sigma) (PYRerun.perpendicular0Matrix rho mu sigma) = (0 : ℤ) := by
  rw [basisCountResidual, perpendicular0_basis_count, perpendicular0_nullity rho mu sigma hr hm hs hs1]
  norm_num

theorem perpendicular0_physical_counts (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) (c : ℝ) (hc : c ≠ 0) :
    matrixNullity ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • PYLoci.perpendicularPoint) 0) (c • PYLoci.perpendicularPoint)) = 1 ∧ ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • PYLoci.perpendicularPoint) 0) (c • PYLoci.perpendicularPoint)).rank = 2 ∧
      matrixNullity (appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • PYLoci.perpendicularPoint) 0) (c • PYLoci.perpendicularPoint)) (c • PYLoci.perpendicularPoint)) = 0 ∧
      (appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • PYLoci.perpendicularPoint) 0) (c • PYLoci.perpendicularPoint)) (c • PYLoci.perpendicularPoint)).rank = 3 := by
  have hM := matrix_counts_of_kernel _ _
    (physical_modal_kernel rho mu sigma (referenceRoot rho mu sigma PYLoci.perpendicularPoint 0) c PYLoci.perpendicularPoint hc)
  have hS := matrix_counts_of_kernel _ _
    (physical_constraint_kernel rho mu sigma (referenceRoot rho mu sigma PYLoci.perpendicularPoint 0) c PYLoci.perpendicularPoint hc)
  have hB : ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma PYLoci.perpendicularPoint 0) PYLoci.perpendicularPoint) = (PYRerun.perpendicular0Matrix rho mu sigma) := (PYRerun.perpendicular0Matrix_reference rho mu sigma hr (ne_of_gt hs)).symm
  have hA : appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma PYLoci.perpendicularPoint 0) PYLoci.perpendicularPoint) PYLoci.perpendicularPoint = (PYRerun.perpendicular0Stack rho mu sigma) :=
    (PYRerun.perpendicular0Stack_reference rho mu sigma hr (ne_of_gt hs)).symm
  rw [hB] at hM
  rw [hA] at hS
  rw [referenceRoot_scale]
  exact ⟨hM.1.trans (perpendicular0_nullity rho mu sigma hr hm hs hs1), hM.2.trans (perpendicular0_rank rho mu sigma hr hm hs hs1),
    hS.1.trans (perpendicular0_transverse_nullity rho mu sigma hr hm hs hs1), hS.2.trans (perpendicular0_stacked_rank rho mu sigma hr hm hs hs1)⟩

theorem perpendicular1_kernel (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    LinearMap.ker (PYRerun.perpendicular1Matrix rho mu sigma).mulVecLin = rootMode sigma PYLoci.perpendicularPoint 1 := by
  ext a
  change (PYRerun.perpendicular1Matrix rho mu sigma).mulVec a = 0 ↔ a ∈ rootMode sigma PYLoci.perpendicularPoint 1
  rw [PYRerun.perpendicular1_basis_complete rho mu sigma a hr hm hs hs1,
    PYRerun.perpendicular1Basis_reference rho mu sigma hr (ne_of_gt hs),
    PYReference.perpendicular1_complete sigma hs hs1]

theorem perpendicular1_stack_kernel (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    LinearMap.ker (PYRerun.perpendicular1Stack rho mu sigma).mulVecLin = rootMode sigma PYLoci.perpendicularPoint 1 ⊓ transverseSpace PYLoci.perpendicularPoint := by
  have he : (PYRerun.perpendicular1Stack rho mu sigma) = appendConstraint (PYRerun.perpendicular1Matrix rho mu sigma) PYLoci.perpendicularPoint := by
    rw [PYRerun.perpendicular1Stack_reference rho mu sigma hr (ne_of_gt hs),
      PYRerun.perpendicular1Matrix_reference rho mu sigma hr (ne_of_gt hs)]
    rfl
  rw [he, appendConstraint_kernel, perpendicular1_kernel rho mu sigma hr hm hs hs1]

theorem perpendicular1_nullity (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : matrixNullity (PYRerun.perpendicular1Matrix rho mu sigma) = 1 := by
  rw [matrixNullity, perpendicular1_kernel rho mu sigma hr hm hs hs1]
  exact (PYReference.perpendicular1_dimensions sigma hs hs1).1

theorem perpendicular1_transverse_nullity (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : matrixNullity (PYRerun.perpendicular1Stack rho mu sigma) = 1 := by
  rw [matrixNullity, perpendicular1_stack_kernel rho mu sigma hr hm hs hs1]
  exact (PYReference.perpendicular1_dimensions sigma hs hs1).2

theorem perpendicular1_rank (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : (PYRerun.perpendicular1Matrix rho mu sigma).rank = 2 := by
  exact matrix_rank_of_nullity (PYRerun.perpendicular1Matrix rho mu sigma) (perpendicular1_nullity rho mu sigma hr hm hs hs1)

theorem perpendicular1_stacked_rank (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : (PYRerun.perpendicular1Stack rho mu sigma).rank = 2 := by
  exact matrix_rank_of_nullity (PYRerun.perpendicular1Stack rho mu sigma) (perpendicular1_transverse_nullity rho mu sigma hr hm hs hs1)

theorem perpendicular1_difference (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : nullityDifference (PYRerun.perpendicular1Matrix rho mu sigma) (PYRerun.perpendicular1Stack rho mu sigma) = (0 : ℤ) := by
  rw [nullityDifference, perpendicular1_nullity rho mu sigma hr hm hs hs1, perpendicular1_transverse_nullity rho mu sigma hr hm hs hs1]
  norm_num

theorem perpendicular1_basis_count (rho mu sigma : ℝ) : basisCount (PYRerun.perpendicular1Basis rho mu sigma) = (1 : ℕ) := by
  rfl

theorem perpendicular1_basis_residual (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : basisCountResidual (PYRerun.perpendicular1Basis rho mu sigma) (PYRerun.perpendicular1Matrix rho mu sigma) = (0 : ℤ) := by
  rw [basisCountResidual, perpendicular1_basis_count, perpendicular1_nullity rho mu sigma hr hm hs hs1]
  norm_num

theorem perpendicular1_physical_counts (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) (c : ℝ) (hc : c ≠ 0) :
    matrixNullity ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • PYLoci.perpendicularPoint) 1) (c • PYLoci.perpendicularPoint)) = 1 ∧ ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • PYLoci.perpendicularPoint) 1) (c • PYLoci.perpendicularPoint)).rank = 2 ∧
      matrixNullity (appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • PYLoci.perpendicularPoint) 1) (c • PYLoci.perpendicularPoint)) (c • PYLoci.perpendicularPoint)) = 1 ∧
      (appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • PYLoci.perpendicularPoint) 1) (c • PYLoci.perpendicularPoint)) (c • PYLoci.perpendicularPoint)).rank = 2 := by
  have hM := matrix_counts_of_kernel _ _
    (physical_modal_kernel rho mu sigma (referenceRoot rho mu sigma PYLoci.perpendicularPoint 1) c PYLoci.perpendicularPoint hc)
  have hS := matrix_counts_of_kernel _ _
    (physical_constraint_kernel rho mu sigma (referenceRoot rho mu sigma PYLoci.perpendicularPoint 1) c PYLoci.perpendicularPoint hc)
  have hB : ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma PYLoci.perpendicularPoint 1) PYLoci.perpendicularPoint) = (PYRerun.perpendicular1Matrix rho mu sigma) := (PYRerun.perpendicular1Matrix_reference rho mu sigma hr (ne_of_gt hs)).symm
  have hA : appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma PYLoci.perpendicularPoint 1) PYLoci.perpendicularPoint) PYLoci.perpendicularPoint = (PYRerun.perpendicular1Stack rho mu sigma) :=
    (PYRerun.perpendicular1Stack_reference rho mu sigma hr (ne_of_gt hs)).symm
  rw [hB] at hM
  rw [hA] at hS
  rw [referenceRoot_scale]
  exact ⟨hM.1.trans (perpendicular1_nullity rho mu sigma hr hm hs hs1), hM.2.trans (perpendicular1_rank rho mu sigma hr hm hs hs1),
    hS.1.trans (perpendicular1_transverse_nullity rho mu sigma hr hm hs hs1), hS.2.trans (perpendicular1_stacked_rank rho mu sigma hr hm hs hs1)⟩

theorem perpendicular2_kernel (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    LinearMap.ker (PYRerun.perpendicular2Matrix rho mu sigma).mulVecLin = rootMode sigma PYLoci.perpendicularPoint 2 := by
  ext a
  change (PYRerun.perpendicular2Matrix rho mu sigma).mulVec a = 0 ↔ a ∈ rootMode sigma PYLoci.perpendicularPoint 2
  rw [PYRerun.perpendicular2_basis_complete rho mu sigma a hr hm hs hs1,
    PYRerun.perpendicular2Basis_reference rho mu sigma hr (ne_of_gt hs),
    PYReference.perpendicular2_complete sigma hs hs1]

theorem perpendicular2_stack_kernel (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    LinearMap.ker (PYRerun.perpendicular2Stack rho mu sigma).mulVecLin = rootMode sigma PYLoci.perpendicularPoint 2 ⊓ transverseSpace PYLoci.perpendicularPoint := by
  have he : (PYRerun.perpendicular2Stack rho mu sigma) = appendConstraint (PYRerun.perpendicular2Matrix rho mu sigma) PYLoci.perpendicularPoint := by
    rw [PYRerun.perpendicular2Stack_reference rho mu sigma hr (ne_of_gt hs),
      PYRerun.perpendicular2Matrix_reference rho mu sigma hr (ne_of_gt hs)]
    rfl
  rw [he, appendConstraint_kernel, perpendicular2_kernel rho mu sigma hr hm hs hs1]

theorem perpendicular2_nullity (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : matrixNullity (PYRerun.perpendicular2Matrix rho mu sigma) = 1 := by
  rw [matrixNullity, perpendicular2_kernel rho mu sigma hr hm hs hs1]
  exact (PYReference.perpendicular2_dimensions sigma hs hs1).1

theorem perpendicular2_transverse_nullity (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : matrixNullity (PYRerun.perpendicular2Stack rho mu sigma) = 1 := by
  rw [matrixNullity, perpendicular2_stack_kernel rho mu sigma hr hm hs hs1]
  exact (PYReference.perpendicular2_dimensions sigma hs hs1).2

theorem perpendicular2_rank (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : (PYRerun.perpendicular2Matrix rho mu sigma).rank = 2 := by
  exact matrix_rank_of_nullity (PYRerun.perpendicular2Matrix rho mu sigma) (perpendicular2_nullity rho mu sigma hr hm hs hs1)

theorem perpendicular2_stacked_rank (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : (PYRerun.perpendicular2Stack rho mu sigma).rank = 2 := by
  exact matrix_rank_of_nullity (PYRerun.perpendicular2Stack rho mu sigma) (perpendicular2_transverse_nullity rho mu sigma hr hm hs hs1)

theorem perpendicular2_difference (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : nullityDifference (PYRerun.perpendicular2Matrix rho mu sigma) (PYRerun.perpendicular2Stack rho mu sigma) = (0 : ℤ) := by
  rw [nullityDifference, perpendicular2_nullity rho mu sigma hr hm hs hs1, perpendicular2_transverse_nullity rho mu sigma hr hm hs hs1]
  norm_num

theorem perpendicular2_basis_count (rho mu sigma : ℝ) : basisCount (PYRerun.perpendicular2Basis rho mu sigma) = (1 : ℕ) := by
  rfl

theorem perpendicular2_basis_residual (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : basisCountResidual (PYRerun.perpendicular2Basis rho mu sigma) (PYRerun.perpendicular2Matrix rho mu sigma) = (0 : ℤ) := by
  rw [basisCountResidual, perpendicular2_basis_count, perpendicular2_nullity rho mu sigma hr hm hs hs1]
  norm_num

theorem perpendicular2_physical_counts (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) (c : ℝ) (hc : c ≠ 0) :
    matrixNullity ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • PYLoci.perpendicularPoint) 2) (c • PYLoci.perpendicularPoint)) = 1 ∧ ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • PYLoci.perpendicularPoint) 2) (c • PYLoci.perpendicularPoint)).rank = 2 ∧
      matrixNullity (appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • PYLoci.perpendicularPoint) 2) (c • PYLoci.perpendicularPoint)) (c • PYLoci.perpendicularPoint)) = 1 ∧
      (appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • PYLoci.perpendicularPoint) 2) (c • PYLoci.perpendicularPoint)) (c • PYLoci.perpendicularPoint)).rank = 2 := by
  have hM := matrix_counts_of_kernel _ _
    (physical_modal_kernel rho mu sigma (referenceRoot rho mu sigma PYLoci.perpendicularPoint 2) c PYLoci.perpendicularPoint hc)
  have hS := matrix_counts_of_kernel _ _
    (physical_constraint_kernel rho mu sigma (referenceRoot rho mu sigma PYLoci.perpendicularPoint 2) c PYLoci.perpendicularPoint hc)
  have hB : ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma PYLoci.perpendicularPoint 2) PYLoci.perpendicularPoint) = (PYRerun.perpendicular2Matrix rho mu sigma) := (PYRerun.perpendicular2Matrix_reference rho mu sigma hr (ne_of_gt hs)).symm
  have hA : appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma PYLoci.perpendicularPoint 2) PYLoci.perpendicularPoint) PYLoci.perpendicularPoint = (PYRerun.perpendicular2Stack rho mu sigma) :=
    (PYRerun.perpendicular2Stack_reference rho mu sigma hr (ne_of_gt hs)).symm
  rw [hB] at hM
  rw [hA] at hS
  rw [referenceRoot_scale]
  exact ⟨hM.1.trans (perpendicular2_nullity rho mu sigma hr hm hs hs1), hM.2.trans (perpendicular2_rank rho mu sigma hr hm hs hs1),
    hS.1.trans (perpendicular2_transverse_nullity rho mu sigma hr hm hs hs1), hS.2.trans (perpendicular2_stacked_rank rho mu sigma hr hm hs hs1)⟩

end PYCountReference

namespace WLCountReference

theorem parallel0_kernel (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    LinearMap.ker (WLRerun.parallel0Matrix rho mu sigma).mulVecLin = rootMode sigma WLLoci.parallelPoint 0 := by
  ext a
  change (WLRerun.parallel0Matrix rho mu sigma).mulVec a = 0 ↔ a ∈ rootMode sigma WLLoci.parallelPoint 0
  rw [WLRerun.parallel0_basis_complete rho mu sigma a hr hm hs hs1,
    WLRerun.parallel0Basis_reference rho mu sigma hr (ne_of_gt hs),
    WLReference.parallel0_complete sigma hs hs1]

theorem parallel0_stack_kernel (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    LinearMap.ker (WLRerun.parallel0Stack rho mu sigma).mulVecLin = rootMode sigma WLLoci.parallelPoint 0 ⊓ transverseSpace WLLoci.parallelPoint := by
  have he : (WLRerun.parallel0Stack rho mu sigma) = appendConstraint (WLRerun.parallel0Matrix rho mu sigma) WLLoci.parallelPoint := by
    rw [WLRerun.parallel0Stack_reference rho mu sigma hr (ne_of_gt hs),
      WLRerun.parallel0Matrix_reference rho mu sigma hr (ne_of_gt hs)]
    rfl
  rw [he, appendConstraint_kernel, parallel0_kernel rho mu sigma hr hm hs hs1]

theorem parallel0_nullity (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : matrixNullity (WLRerun.parallel0Matrix rho mu sigma) = 1 := by
  rw [matrixNullity, parallel0_kernel rho mu sigma hr hm hs hs1]
  exact (WLReference.parallel0_dimensions sigma hs hs1).1

theorem parallel0_transverse_nullity (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : matrixNullity (WLRerun.parallel0Stack rho mu sigma) = 0 := by
  rw [matrixNullity, parallel0_stack_kernel rho mu sigma hr hm hs hs1]
  exact (WLReference.parallel0_dimensions sigma hs hs1).2

theorem parallel0_rank (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : (WLRerun.parallel0Matrix rho mu sigma).rank = 2 := by
  exact matrix_rank_of_nullity (WLRerun.parallel0Matrix rho mu sigma) (parallel0_nullity rho mu sigma hr hm hs hs1)

theorem parallel0_stacked_rank (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : (WLRerun.parallel0Stack rho mu sigma).rank = 3 := by
  exact matrix_rank_of_nullity (WLRerun.parallel0Stack rho mu sigma) (parallel0_transverse_nullity rho mu sigma hr hm hs hs1)

theorem parallel0_difference (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : nullityDifference (WLRerun.parallel0Matrix rho mu sigma) (WLRerun.parallel0Stack rho mu sigma) = (1 : ℤ) := by
  rw [nullityDifference, parallel0_nullity rho mu sigma hr hm hs hs1, parallel0_transverse_nullity rho mu sigma hr hm hs hs1]
  norm_num

theorem parallel0_basis_count (rho mu sigma : ℝ) : basisCount (WLRerun.parallel0Basis rho mu sigma) = (1 : ℕ) := by
  rfl

theorem parallel0_basis_residual (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : basisCountResidual (WLRerun.parallel0Basis rho mu sigma) (WLRerun.parallel0Matrix rho mu sigma) = (0 : ℤ) := by
  rw [basisCountResidual, parallel0_basis_count, parallel0_nullity rho mu sigma hr hm hs hs1]
  norm_num

theorem parallel0_physical_counts (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) (c : ℝ) (hc : c ≠ 0) :
    matrixNullity ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • WLLoci.parallelPoint) 0) (c • WLLoci.parallelPoint)) = 1 ∧ ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • WLLoci.parallelPoint) 0) (c • WLLoci.parallelPoint)).rank = 2 ∧
      matrixNullity (appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • WLLoci.parallelPoint) 0) (c • WLLoci.parallelPoint)) (c • WLLoci.parallelPoint)) = 0 ∧
      (appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • WLLoci.parallelPoint) 0) (c • WLLoci.parallelPoint)) (c • WLLoci.parallelPoint)).rank = 3 := by
  have hM := matrix_counts_of_kernel _ _
    (physical_modal_kernel rho mu sigma (referenceRoot rho mu sigma WLLoci.parallelPoint 0) c WLLoci.parallelPoint hc)
  have hS := matrix_counts_of_kernel _ _
    (physical_constraint_kernel rho mu sigma (referenceRoot rho mu sigma WLLoci.parallelPoint 0) c WLLoci.parallelPoint hc)
  have hB : ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma WLLoci.parallelPoint 0) WLLoci.parallelPoint) = (WLRerun.parallel0Matrix rho mu sigma) := (WLRerun.parallel0Matrix_reference rho mu sigma hr (ne_of_gt hs)).symm
  have hA : appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma WLLoci.parallelPoint 0) WLLoci.parallelPoint) WLLoci.parallelPoint = (WLRerun.parallel0Stack rho mu sigma) :=
    (WLRerun.parallel0Stack_reference rho mu sigma hr (ne_of_gt hs)).symm
  rw [hB] at hM
  rw [hA] at hS
  rw [referenceRoot_scale]
  exact ⟨hM.1.trans (parallel0_nullity rho mu sigma hr hm hs hs1), hM.2.trans (parallel0_rank rho mu sigma hr hm hs hs1),
    hS.1.trans (parallel0_transverse_nullity rho mu sigma hr hm hs hs1), hS.2.trans (parallel0_stacked_rank rho mu sigma hr hm hs hs1)⟩

theorem parallel1_kernel (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    LinearMap.ker (WLRerun.parallel1Matrix rho mu sigma).mulVecLin = rootMode sigma WLLoci.parallelPoint 1 := by
  ext a
  change (WLRerun.parallel1Matrix rho mu sigma).mulVec a = 0 ↔ a ∈ rootMode sigma WLLoci.parallelPoint 1
  rw [WLRerun.parallel1_basis_complete rho mu sigma a hr hm hs hs1,
    WLRerun.parallel1Basis_reference rho mu sigma hr (ne_of_gt hs),
    WLReference.parallel1_complete sigma hs hs1]

theorem parallel1_stack_kernel (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    LinearMap.ker (WLRerun.parallel1Stack rho mu sigma).mulVecLin = rootMode sigma WLLoci.parallelPoint 1 ⊓ transverseSpace WLLoci.parallelPoint := by
  have he : (WLRerun.parallel1Stack rho mu sigma) = appendConstraint (WLRerun.parallel1Matrix rho mu sigma) WLLoci.parallelPoint := by
    rw [WLRerun.parallel1Stack_reference rho mu sigma hr (ne_of_gt hs),
      WLRerun.parallel1Matrix_reference rho mu sigma hr (ne_of_gt hs)]
    rfl
  rw [he, appendConstraint_kernel, parallel1_kernel rho mu sigma hr hm hs hs1]

theorem parallel1_nullity (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : matrixNullity (WLRerun.parallel1Matrix rho mu sigma) = 2 := by
  rw [matrixNullity, parallel1_kernel rho mu sigma hr hm hs hs1]
  exact (WLReference.parallel1_dimensions sigma hs hs1).1

theorem parallel1_transverse_nullity (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : matrixNullity (WLRerun.parallel1Stack rho mu sigma) = 2 := by
  rw [matrixNullity, parallel1_stack_kernel rho mu sigma hr hm hs hs1]
  exact (WLReference.parallel1_dimensions sigma hs hs1).2

theorem parallel1_rank (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : (WLRerun.parallel1Matrix rho mu sigma).rank = 1 := by
  exact matrix_rank_of_nullity (WLRerun.parallel1Matrix rho mu sigma) (parallel1_nullity rho mu sigma hr hm hs hs1)

theorem parallel1_stacked_rank (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : (WLRerun.parallel1Stack rho mu sigma).rank = 1 := by
  exact matrix_rank_of_nullity (WLRerun.parallel1Stack rho mu sigma) (parallel1_transverse_nullity rho mu sigma hr hm hs hs1)

theorem parallel1_difference (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : nullityDifference (WLRerun.parallel1Matrix rho mu sigma) (WLRerun.parallel1Stack rho mu sigma) = (0 : ℤ) := by
  rw [nullityDifference, parallel1_nullity rho mu sigma hr hm hs hs1, parallel1_transverse_nullity rho mu sigma hr hm hs hs1]
  norm_num

theorem parallel1_basis_count (rho mu sigma : ℝ) : basisCount (WLRerun.parallel1Basis rho mu sigma) = (2 : ℕ) := by
  rfl

theorem parallel1_basis_residual (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : basisCountResidual (WLRerun.parallel1Basis rho mu sigma) (WLRerun.parallel1Matrix rho mu sigma) = (0 : ℤ) := by
  rw [basisCountResidual, parallel1_basis_count, parallel1_nullity rho mu sigma hr hm hs hs1]
  norm_num

theorem parallel1_physical_counts (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) (c : ℝ) (hc : c ≠ 0) :
    matrixNullity ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • WLLoci.parallelPoint) 1) (c • WLLoci.parallelPoint)) = 2 ∧ ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • WLLoci.parallelPoint) 1) (c • WLLoci.parallelPoint)).rank = 1 ∧
      matrixNullity (appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • WLLoci.parallelPoint) 1) (c • WLLoci.parallelPoint)) (c • WLLoci.parallelPoint)) = 2 ∧
      (appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • WLLoci.parallelPoint) 1) (c • WLLoci.parallelPoint)) (c • WLLoci.parallelPoint)).rank = 1 := by
  have hM := matrix_counts_of_kernel _ _
    (physical_modal_kernel rho mu sigma (referenceRoot rho mu sigma WLLoci.parallelPoint 1) c WLLoci.parallelPoint hc)
  have hS := matrix_counts_of_kernel _ _
    (physical_constraint_kernel rho mu sigma (referenceRoot rho mu sigma WLLoci.parallelPoint 1) c WLLoci.parallelPoint hc)
  have hB : ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma WLLoci.parallelPoint 1) WLLoci.parallelPoint) = (WLRerun.parallel1Matrix rho mu sigma) := (WLRerun.parallel1Matrix_reference rho mu sigma hr (ne_of_gt hs)).symm
  have hA : appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma WLLoci.parallelPoint 1) WLLoci.parallelPoint) WLLoci.parallelPoint = (WLRerun.parallel1Stack rho mu sigma) :=
    (WLRerun.parallel1Stack_reference rho mu sigma hr (ne_of_gt hs)).symm
  rw [hB] at hM
  rw [hA] at hS
  rw [referenceRoot_scale]
  exact ⟨hM.1.trans (parallel1_nullity rho mu sigma hr hm hs hs1), hM.2.trans (parallel1_rank rho mu sigma hr hm hs hs1),
    hS.1.trans (parallel1_transverse_nullity rho mu sigma hr hm hs hs1), hS.2.trans (parallel1_stacked_rank rho mu sigma hr hm hs hs1)⟩

theorem perpendicular0_kernel (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    LinearMap.ker (WLRerun.perpendicular0Matrix rho mu sigma).mulVecLin = rootMode sigma WLLoci.perpendicularPoint 0 := by
  ext a
  change (WLRerun.perpendicular0Matrix rho mu sigma).mulVec a = 0 ↔ a ∈ rootMode sigma WLLoci.perpendicularPoint 0
  rw [WLRerun.perpendicular0_basis_complete rho mu sigma a hr hm hs hs1,
    WLRerun.perpendicular0Basis_reference rho mu sigma hr (ne_of_gt hs),
    WLReference.perpendicular0_complete sigma hs hs1]

theorem perpendicular0_stack_kernel (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    LinearMap.ker (WLRerun.perpendicular0Stack rho mu sigma).mulVecLin = rootMode sigma WLLoci.perpendicularPoint 0 ⊓ transverseSpace WLLoci.perpendicularPoint := by
  have he : (WLRerun.perpendicular0Stack rho mu sigma) = appendConstraint (WLRerun.perpendicular0Matrix rho mu sigma) WLLoci.perpendicularPoint := by
    rw [WLRerun.perpendicular0Stack_reference rho mu sigma hr (ne_of_gt hs),
      WLRerun.perpendicular0Matrix_reference rho mu sigma hr (ne_of_gt hs)]
    rfl
  rw [he, appendConstraint_kernel, perpendicular0_kernel rho mu sigma hr hm hs hs1]

theorem perpendicular0_nullity (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : matrixNullity (WLRerun.perpendicular0Matrix rho mu sigma) = 1 := by
  rw [matrixNullity, perpendicular0_kernel rho mu sigma hr hm hs hs1]
  exact (WLReference.perpendicular0_dimensions sigma hs hs1).1

theorem perpendicular0_transverse_nullity (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : matrixNullity (WLRerun.perpendicular0Stack rho mu sigma) = 0 := by
  rw [matrixNullity, perpendicular0_stack_kernel rho mu sigma hr hm hs hs1]
  exact (WLReference.perpendicular0_dimensions sigma hs hs1).2

theorem perpendicular0_rank (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : (WLRerun.perpendicular0Matrix rho mu sigma).rank = 2 := by
  exact matrix_rank_of_nullity (WLRerun.perpendicular0Matrix rho mu sigma) (perpendicular0_nullity rho mu sigma hr hm hs hs1)

theorem perpendicular0_stacked_rank (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : (WLRerun.perpendicular0Stack rho mu sigma).rank = 3 := by
  exact matrix_rank_of_nullity (WLRerun.perpendicular0Stack rho mu sigma) (perpendicular0_transverse_nullity rho mu sigma hr hm hs hs1)

theorem perpendicular0_difference (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : nullityDifference (WLRerun.perpendicular0Matrix rho mu sigma) (WLRerun.perpendicular0Stack rho mu sigma) = (1 : ℤ) := by
  rw [nullityDifference, perpendicular0_nullity rho mu sigma hr hm hs hs1, perpendicular0_transverse_nullity rho mu sigma hr hm hs hs1]
  norm_num

theorem perpendicular0_basis_count (rho mu sigma : ℝ) : basisCount (WLRerun.perpendicular0Basis rho mu sigma) = (1 : ℕ) := by
  rfl

theorem perpendicular0_basis_residual (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : basisCountResidual (WLRerun.perpendicular0Basis rho mu sigma) (WLRerun.perpendicular0Matrix rho mu sigma) = (0 : ℤ) := by
  rw [basisCountResidual, perpendicular0_basis_count, perpendicular0_nullity rho mu sigma hr hm hs hs1]
  norm_num

theorem perpendicular0_physical_counts (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) (c : ℝ) (hc : c ≠ 0) :
    matrixNullity ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • WLLoci.perpendicularPoint) 0) (c • WLLoci.perpendicularPoint)) = 1 ∧ ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • WLLoci.perpendicularPoint) 0) (c • WLLoci.perpendicularPoint)).rank = 2 ∧
      matrixNullity (appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • WLLoci.perpendicularPoint) 0) (c • WLLoci.perpendicularPoint)) (c • WLLoci.perpendicularPoint)) = 0 ∧
      (appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • WLLoci.perpendicularPoint) 0) (c • WLLoci.perpendicularPoint)) (c • WLLoci.perpendicularPoint)).rank = 3 := by
  have hM := matrix_counts_of_kernel _ _
    (physical_modal_kernel rho mu sigma (referenceRoot rho mu sigma WLLoci.perpendicularPoint 0) c WLLoci.perpendicularPoint hc)
  have hS := matrix_counts_of_kernel _ _
    (physical_constraint_kernel rho mu sigma (referenceRoot rho mu sigma WLLoci.perpendicularPoint 0) c WLLoci.perpendicularPoint hc)
  have hB : ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma WLLoci.perpendicularPoint 0) WLLoci.perpendicularPoint) = (WLRerun.perpendicular0Matrix rho mu sigma) := (WLRerun.perpendicular0Matrix_reference rho mu sigma hr (ne_of_gt hs)).symm
  have hA : appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma WLLoci.perpendicularPoint 0) WLLoci.perpendicularPoint) WLLoci.perpendicularPoint = (WLRerun.perpendicular0Stack rho mu sigma) :=
    (WLRerun.perpendicular0Stack_reference rho mu sigma hr (ne_of_gt hs)).symm
  rw [hB] at hM
  rw [hA] at hS
  rw [referenceRoot_scale]
  exact ⟨hM.1.trans (perpendicular0_nullity rho mu sigma hr hm hs hs1), hM.2.trans (perpendicular0_rank rho mu sigma hr hm hs hs1),
    hS.1.trans (perpendicular0_transverse_nullity rho mu sigma hr hm hs hs1), hS.2.trans (perpendicular0_stacked_rank rho mu sigma hr hm hs hs1)⟩

theorem perpendicular1_kernel (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    LinearMap.ker (WLRerun.perpendicular1Matrix rho mu sigma).mulVecLin = rootMode sigma WLLoci.perpendicularPoint 1 := by
  ext a
  change (WLRerun.perpendicular1Matrix rho mu sigma).mulVec a = 0 ↔ a ∈ rootMode sigma WLLoci.perpendicularPoint 1
  rw [WLRerun.perpendicular1_basis_complete rho mu sigma a hr hm hs hs1,
    WLRerun.perpendicular1Basis_reference rho mu sigma hr (ne_of_gt hs),
    WLReference.perpendicular1_complete sigma hs hs1]

theorem perpendicular1_stack_kernel (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    LinearMap.ker (WLRerun.perpendicular1Stack rho mu sigma).mulVecLin = rootMode sigma WLLoci.perpendicularPoint 1 ⊓ transverseSpace WLLoci.perpendicularPoint := by
  have he : (WLRerun.perpendicular1Stack rho mu sigma) = appendConstraint (WLRerun.perpendicular1Matrix rho mu sigma) WLLoci.perpendicularPoint := by
    rw [WLRerun.perpendicular1Stack_reference rho mu sigma hr (ne_of_gt hs),
      WLRerun.perpendicular1Matrix_reference rho mu sigma hr (ne_of_gt hs)]
    rfl
  rw [he, appendConstraint_kernel, perpendicular1_kernel rho mu sigma hr hm hs hs1]

theorem perpendicular1_nullity (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : matrixNullity (WLRerun.perpendicular1Matrix rho mu sigma) = 1 := by
  rw [matrixNullity, perpendicular1_kernel rho mu sigma hr hm hs hs1]
  exact (WLReference.perpendicular1_dimensions sigma hs hs1).1

theorem perpendicular1_transverse_nullity (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : matrixNullity (WLRerun.perpendicular1Stack rho mu sigma) = 1 := by
  rw [matrixNullity, perpendicular1_stack_kernel rho mu sigma hr hm hs hs1]
  exact (WLReference.perpendicular1_dimensions sigma hs hs1).2

theorem perpendicular1_rank (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : (WLRerun.perpendicular1Matrix rho mu sigma).rank = 2 := by
  exact matrix_rank_of_nullity (WLRerun.perpendicular1Matrix rho mu sigma) (perpendicular1_nullity rho mu sigma hr hm hs hs1)

theorem perpendicular1_stacked_rank (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : (WLRerun.perpendicular1Stack rho mu sigma).rank = 2 := by
  exact matrix_rank_of_nullity (WLRerun.perpendicular1Stack rho mu sigma) (perpendicular1_transverse_nullity rho mu sigma hr hm hs hs1)

theorem perpendicular1_difference (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : nullityDifference (WLRerun.perpendicular1Matrix rho mu sigma) (WLRerun.perpendicular1Stack rho mu sigma) = (0 : ℤ) := by
  rw [nullityDifference, perpendicular1_nullity rho mu sigma hr hm hs hs1, perpendicular1_transverse_nullity rho mu sigma hr hm hs hs1]
  norm_num

theorem perpendicular1_basis_count (rho mu sigma : ℝ) : basisCount (WLRerun.perpendicular1Basis rho mu sigma) = (1 : ℕ) := by
  rfl

theorem perpendicular1_basis_residual (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : basisCountResidual (WLRerun.perpendicular1Basis rho mu sigma) (WLRerun.perpendicular1Matrix rho mu sigma) = (0 : ℤ) := by
  rw [basisCountResidual, perpendicular1_basis_count, perpendicular1_nullity rho mu sigma hr hm hs hs1]
  norm_num

theorem perpendicular1_physical_counts (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) (c : ℝ) (hc : c ≠ 0) :
    matrixNullity ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • WLLoci.perpendicularPoint) 1) (c • WLLoci.perpendicularPoint)) = 1 ∧ ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • WLLoci.perpendicularPoint) 1) (c • WLLoci.perpendicularPoint)).rank = 2 ∧
      matrixNullity (appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • WLLoci.perpendicularPoint) 1) (c • WLLoci.perpendicularPoint)) (c • WLLoci.perpendicularPoint)) = 1 ∧
      (appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • WLLoci.perpendicularPoint) 1) (c • WLLoci.perpendicularPoint)) (c • WLLoci.perpendicularPoint)).rank = 2 := by
  have hM := matrix_counts_of_kernel _ _
    (physical_modal_kernel rho mu sigma (referenceRoot rho mu sigma WLLoci.perpendicularPoint 1) c WLLoci.perpendicularPoint hc)
  have hS := matrix_counts_of_kernel _ _
    (physical_constraint_kernel rho mu sigma (referenceRoot rho mu sigma WLLoci.perpendicularPoint 1) c WLLoci.perpendicularPoint hc)
  have hB : ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma WLLoci.perpendicularPoint 1) WLLoci.perpendicularPoint) = (WLRerun.perpendicular1Matrix rho mu sigma) := (WLRerun.perpendicular1Matrix_reference rho mu sigma hr (ne_of_gt hs)).symm
  have hA : appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma WLLoci.perpendicularPoint 1) WLLoci.perpendicularPoint) WLLoci.perpendicularPoint = (WLRerun.perpendicular1Stack rho mu sigma) :=
    (WLRerun.perpendicular1Stack_reference rho mu sigma hr (ne_of_gt hs)).symm
  rw [hB] at hM
  rw [hA] at hS
  rw [referenceRoot_scale]
  exact ⟨hM.1.trans (perpendicular1_nullity rho mu sigma hr hm hs hs1), hM.2.trans (perpendicular1_rank rho mu sigma hr hm hs hs1),
    hS.1.trans (perpendicular1_transverse_nullity rho mu sigma hr hm hs hs1), hS.2.trans (perpendicular1_stacked_rank rho mu sigma hr hm hs hs1)⟩

theorem perpendicular2_kernel (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    LinearMap.ker (WLRerun.perpendicular2Matrix rho mu sigma).mulVecLin = rootMode sigma WLLoci.perpendicularPoint 2 := by
  ext a
  change (WLRerun.perpendicular2Matrix rho mu sigma).mulVec a = 0 ↔ a ∈ rootMode sigma WLLoci.perpendicularPoint 2
  rw [WLRerun.perpendicular2_basis_complete rho mu sigma a hr hm hs hs1,
    WLRerun.perpendicular2Basis_reference rho mu sigma hr (ne_of_gt hs),
    WLReference.perpendicular2_complete sigma hs hs1]

theorem perpendicular2_stack_kernel (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    LinearMap.ker (WLRerun.perpendicular2Stack rho mu sigma).mulVecLin = rootMode sigma WLLoci.perpendicularPoint 2 ⊓ transverseSpace WLLoci.perpendicularPoint := by
  have he : (WLRerun.perpendicular2Stack rho mu sigma) = appendConstraint (WLRerun.perpendicular2Matrix rho mu sigma) WLLoci.perpendicularPoint := by
    rw [WLRerun.perpendicular2Stack_reference rho mu sigma hr (ne_of_gt hs),
      WLRerun.perpendicular2Matrix_reference rho mu sigma hr (ne_of_gt hs)]
    rfl
  rw [he, appendConstraint_kernel, perpendicular2_kernel rho mu sigma hr hm hs hs1]

theorem perpendicular2_nullity (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : matrixNullity (WLRerun.perpendicular2Matrix rho mu sigma) = 1 := by
  rw [matrixNullity, perpendicular2_kernel rho mu sigma hr hm hs hs1]
  exact (WLReference.perpendicular2_dimensions sigma hs hs1).1

theorem perpendicular2_transverse_nullity (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : matrixNullity (WLRerun.perpendicular2Stack rho mu sigma) = 1 := by
  rw [matrixNullity, perpendicular2_stack_kernel rho mu sigma hr hm hs hs1]
  exact (WLReference.perpendicular2_dimensions sigma hs hs1).2

theorem perpendicular2_rank (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : (WLRerun.perpendicular2Matrix rho mu sigma).rank = 2 := by
  exact matrix_rank_of_nullity (WLRerun.perpendicular2Matrix rho mu sigma) (perpendicular2_nullity rho mu sigma hr hm hs hs1)

theorem perpendicular2_stacked_rank (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : (WLRerun.perpendicular2Stack rho mu sigma).rank = 2 := by
  exact matrix_rank_of_nullity (WLRerun.perpendicular2Stack rho mu sigma) (perpendicular2_transverse_nullity rho mu sigma hr hm hs hs1)

theorem perpendicular2_difference (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : nullityDifference (WLRerun.perpendicular2Matrix rho mu sigma) (WLRerun.perpendicular2Stack rho mu sigma) = (0 : ℤ) := by
  rw [nullityDifference, perpendicular2_nullity rho mu sigma hr hm hs hs1, perpendicular2_transverse_nullity rho mu sigma hr hm hs hs1]
  norm_num

theorem perpendicular2_basis_count (rho mu sigma : ℝ) : basisCount (WLRerun.perpendicular2Basis rho mu sigma) = (1 : ℕ) := by
  rfl

theorem perpendicular2_basis_residual (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) : basisCountResidual (WLRerun.perpendicular2Basis rho mu sigma) (WLRerun.perpendicular2Matrix rho mu sigma) = (0 : ℤ) := by
  rw [basisCountResidual, perpendicular2_basis_count, perpendicular2_nullity rho mu sigma hr hm hs hs1]
  norm_num

theorem perpendicular2_physical_counts (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) (c : ℝ) (hc : c ≠ 0) :
    matrixNullity ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • WLLoci.perpendicularPoint) 2) (c • WLLoci.perpendicularPoint)) = 1 ∧ ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • WLLoci.perpendicularPoint) 2) (c • WLLoci.perpendicularPoint)).rank = 2 ∧
      matrixNullity (appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • WLLoci.perpendicularPoint) 2) (c • WLLoci.perpendicularPoint)) (c • WLLoci.perpendicularPoint)) = 1 ∧
      (appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • WLLoci.perpendicularPoint) 2) (c • WLLoci.perpendicularPoint)) (c • WLLoci.perpendicularPoint)).rank = 2 := by
  have hM := matrix_counts_of_kernel _ _
    (physical_modal_kernel rho mu sigma (referenceRoot rho mu sigma WLLoci.perpendicularPoint 2) c WLLoci.perpendicularPoint hc)
  have hS := matrix_counts_of_kernel _ _
    (physical_constraint_kernel rho mu sigma (referenceRoot rho mu sigma WLLoci.perpendicularPoint 2) c WLLoci.perpendicularPoint hc)
  have hB : ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma WLLoci.perpendicularPoint 2) WLLoci.perpendicularPoint) = (WLRerun.perpendicular2Matrix rho mu sigma) := (WLRerun.perpendicular2Matrix_reference rho mu sigma hr (ne_of_gt hs)).symm
  have hA : appendConstraint ((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma WLLoci.perpendicularPoint 2) WLLoci.perpendicularPoint) WLLoci.perpendicularPoint = (WLRerun.perpendicular2Stack rho mu sigma) :=
    (WLRerun.perpendicular2Stack_reference rho mu sigma hr (ne_of_gt hs)).symm
  rw [hB] at hM
  rw [hA] at hS
  rw [referenceRoot_scale]
  exact ⟨hM.1.trans (perpendicular2_nullity rho mu sigma hr hm hs hs1), hM.2.trans (perpendicular2_rank rho mu sigma hr hm hs hs1),
    hS.1.trans (perpendicular2_transverse_nullity rho mu sigma hr hm hs hs1), hS.2.trans (perpendicular2_stacked_rank rho mu sigma hr hm hs hs1)⟩

end WLCountReference

end
end S10Audit.CAS
