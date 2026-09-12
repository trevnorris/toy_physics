import S10Audit.CAS.RootSupport
import S10Audit.CAS.PYRootsGeneric
import S10Audit.CAS.PYRootsParallel
import S10Audit.CAS.PYRootsPerpendicular
import S10Audit.CAS.WLRootsGeneric
import S10Audit.CAS.WLRootsParallel
import S10Audit.CAS.WLRootsPerpendicular

set_option backward.isDefEq.respectTransparency false

namespace S10Audit.CAS
open S10Pilot S10Anisotropic Polynomial
noncomputable section

namespace PYRootsGeneric

def candidatesTrees : List (Expr Symbol) := [n0, n17, n24]
def candidates (rho mu sigma _z : ℝ) (_k : Vec 3) : List ℝ := candidatesTrees.map (Expr.eval (values rho mu sigma _z _k))

theorem candidates_reference (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) : candidates rho mu sigma _z _k = [referenceRoot rho mu sigma _k 0, referenceRoot rho mu sigma _k 1, referenceRoot rho mu sigma _k 2] := by
  unfold candidates candidatesTrees
  simp only [List.map_cons, List.map_nil]
  rw [n0_eval rho mu sigma _z _k, q3_root_solutions_raw_cell1 rho mu sigma _z _k _hr, q3_root_solutions_raw_cell2 rho mu sigma _z _k _hr _hs]
  rfl

def distinctTrees : List (Expr Symbol) := [n5, n17, n24]
def distinct (rho mu sigma _z : ℝ) (_k : Vec 3) : List ℝ := distinctTrees.map (Expr.eval (values rho mu sigma _z _k))

theorem distinct_reference (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) : distinct rho mu sigma _z _k = [referenceRoot rho mu sigma _k 0, referenceRoot rho mu sigma _k 1, referenceRoot rho mu sigma _k 2] := by
  unfold distinct distinctTrees
  simp only [List.map_cons, List.map_nil]
  rw [q3_roots_distinct_cell0 rho mu sigma _z _k, q3_roots_distinct_cell1 rho mu sigma _z _k _hr, q3_roots_distinct_cell2 rho mu sigma _z _k _hr _hs]

theorem candidates_complete (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) (_hg : GenericChart sigma _k) (w : ℝ) :
    Matrix.det ((1/2 : ℝ) • referenceMatrix rho mu sigma w _k) = 0 ↔ w ∈ candidates rho mu sigma _z _k := by
  rw [rootPolynomial_complete rho mu sigma w _k _hr (ne_of_gt _hs),
    rootPolynomial_roots rho mu sigma _k _hr (ne_of_gt _hs), candidates_reference rho mu sigma _z _k _hr (ne_of_gt _hs)]
  simp [referenceRoot]

theorem distinct_complete (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) (_hg : GenericChart sigma _k) (w : ℝ) :
    Matrix.det ((1/2 : ℝ) • referenceMatrix rho mu sigma w _k) = 0 ↔ w ∈ distinct rho mu sigma _z _k := by
  rw [rootPolynomial_complete rho mu sigma w _k _hr (ne_of_gt _hs),
    rootPolynomial_roots rho mu sigma _k _hr (ne_of_gt _hs), distinct_reference rho mu sigma _z _k _hr (ne_of_gt _hs)]
  simp [referenceRoot]

theorem emitted_determinant_complete (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) (_hg : GenericChart sigma _k) (w : ℝ) :
    Expr.eval (values rho mu sigma w _k) PY.n161 = 0 ↔ w ∈ distinct rho mu sigma _z _k := by
  rw [PY.q3_determinant_cell0 rho mu sigma w _k]
  exact distinct_complete rho mu sigma _z _k _hr _hm _hs _hs1 _hg w

theorem distinct_nodup (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) (_hg : GenericChart sigma _k) : (distinct rho mu sigma _z _k).Nodup := by
  have hk : _k ≠ 0 := (genericChart_split sigma _k _hg).1
  have hq : perpSq 0 _k ≠ 0 := (genericChart_split sigma _k _hg).2
  rw [distinct_reference rho mu sigma _z _k _hr (ne_of_gt _hs)]
  exact split_root_nodup rho mu sigma _k _hr _hm _hs _hs1 hk hq

theorem multiplicities (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) (_hg : GenericChart sigma _k) :
    rootMultiplicity 0 (rootPolynomial rho mu sigma _k) = 1 ∧
      rootMultiplicity (referenceRoot rho mu sigma _k 1) (rootPolynomial rho mu sigma _k) = 1 ∧
      rootMultiplicity (referenceRoot rho mu sigma _k 2) (rootPolynomial rho mu sigma _k) = 1
    := by
  have hk : _k ≠ 0 := (genericChart_split sigma _k _hg).1
  have hq : perpSq 0 _k ≠ 0 := (genericChart_split sigma _k _hg).2
  exact split_root_multiplicities rho mu sigma _k _hr _hm _hs _hs1 hk hq

theorem candidate_multiset (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) (_hg : GenericChart sigma _k) :
    (rootPolynomial rho mu sigma _k).roots = (candidates rho mu sigma _z _k : Multiset ℝ) := by
  rw [rootPolynomial_roots rho mu sigma _k _hr (ne_of_gt _hs), candidates_reference rho mu sigma _z _k _hr (ne_of_gt _hs)]
  rfl

theorem distinct_count (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) (_hg : GenericChart sigma _k) :
    Expr.eval (values rho mu sigma _z _k) n25 = ((distinct rho mu sigma _z _k).toFinset.card : ℝ) := by
  classical
  rw [List.toFinset_card_of_nodup (distinct_nodup rho mu sigma _z _k _hr _hm _hs _hs1 _hg),
    q3_root_count_cell0 rho mu sigma _z _k, distinct_reference rho mu sigma _z _k _hr (ne_of_gt _hs)]
  norm_num

end PYRootsGeneric

namespace PYRootsParallel

def candidatesTrees : List (Expr Symbol) := [n0, n3]
def candidates (rho mu sigma _z : ℝ) (_k : Vec 3) : List ℝ := candidatesTrees.map (Expr.eval (values rho mu sigma _z _k))

theorem candidates_reference (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) : candidates rho mu sigma _z _k = [referenceRoot rho mu sigma PYLoci.parallelPoint 0, referenceRoot rho mu sigma PYLoci.parallelPoint 1] := by
  unfold candidates candidatesTrees
  simp only [List.map_cons, List.map_nil]
  rw [n0_eval rho mu sigma _z _k, q8_stratum1_q3_root_solutions_raw_cell1 rho mu sigma _z _k _hr]
  rfl

def distinctTrees : List (Expr Symbol) := [n4, n3]
def distinct (rho mu sigma _z : ℝ) (_k : Vec 3) : List ℝ := distinctTrees.map (Expr.eval (values rho mu sigma _z _k))

theorem distinct_reference (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) : distinct rho mu sigma _z _k = [referenceRoot rho mu sigma PYLoci.parallelPoint 0, referenceRoot rho mu sigma PYLoci.parallelPoint 1] := by
  unfold distinct distinctTrees
  simp only [List.map_cons, List.map_nil]
  rw [q8_stratum1_q3_roots_distinct_cell0 rho mu sigma _z _k, q8_stratum1_q3_roots_distinct_cell1 rho mu sigma _z _k _hr]

theorem candidates_complete (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) (w : ℝ) :
    Matrix.det ((1/2 : ℝ) • referenceMatrix rho mu sigma w PYLoci.parallelPoint) = 0 ↔ w ∈ candidates rho mu sigma _z _k := by
  have hq : perpSq 0 PYLoci.parallelPoint = 0 := by
    norm_num [PYLoci.parallelPoint, perpSq, normSq, dot, Fin.sum_univ_three, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]
  have he := (root_coincidence rho mu sigma PYLoci.parallelPoint _hr _hm (ne_of_gt _hs) _hs1).mpr hq
  rw [rootPolynomial_complete rho mu sigma w PYLoci.parallelPoint _hr (ne_of_gt _hs),
    rootPolynomial_roots rho mu sigma PYLoci.parallelPoint _hr (ne_of_gt _hs), candidates_reference rho mu sigma _z _k _hr (ne_of_gt _hs)]
  rw [he]
  simp [referenceRoot]

theorem distinct_complete (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) (w : ℝ) :
    Matrix.det ((1/2 : ℝ) • referenceMatrix rho mu sigma w PYLoci.parallelPoint) = 0 ↔ w ∈ distinct rho mu sigma _z _k := by
  have hq : perpSq 0 PYLoci.parallelPoint = 0 := by
    norm_num [PYLoci.parallelPoint, perpSq, normSq, dot, Fin.sum_univ_three, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]
  have he := (root_coincidence rho mu sigma PYLoci.parallelPoint _hr _hm (ne_of_gt _hs) _hs1).mpr hq
  rw [rootPolynomial_complete rho mu sigma w PYLoci.parallelPoint _hr (ne_of_gt _hs),
    rootPolynomial_roots rho mu sigma PYLoci.parallelPoint _hr (ne_of_gt _hs), distinct_reference rho mu sigma _z _k _hr (ne_of_gt _hs)]
  rw [he]
  simp [referenceRoot]

theorem emitted_determinant_complete (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) (w : ℝ) :
    Expr.eval (values rho mu sigma w _k) PYRerun.n23 = 0 ↔ w ∈ distinct rho mu sigma _z _k := by
  rw [PYRerun.q8_stratum1_q3_determinant_cell0 rho mu sigma w _k]
  exact distinct_complete rho mu sigma _z _k _hr _hm _hs _hs1 w

theorem distinct_nodup (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) : (distinct rho mu sigma _z _k).Nodup := by
  have hk : PYLoci.parallelPoint ≠ 0 := PYLoci.parallelPoint_nonzero
  have hq : perpSq 0 PYLoci.parallelPoint = 0 := by
    norm_num [PYLoci.parallelPoint, perpSq, normSq, dot, Fin.sum_univ_three, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]
  rw [distinct_reference rho mu sigma _z _k _hr (ne_of_gt _hs)]
  exact parallel_root_nodup rho mu sigma PYLoci.parallelPoint _hr _hm _hs hk

theorem multiplicities (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    rootMultiplicity 0 (rootPolynomial rho mu sigma PYLoci.parallelPoint) = 1 ∧
      rootMultiplicity (referenceRoot rho mu sigma PYLoci.parallelPoint 1) (rootPolynomial rho mu sigma PYLoci.parallelPoint) = 2
    := by
  have hk : PYLoci.parallelPoint ≠ 0 := PYLoci.parallelPoint_nonzero
  have hq : perpSq 0 PYLoci.parallelPoint = 0 := by
    norm_num [PYLoci.parallelPoint, perpSq, normSq, dot, Fin.sum_univ_three, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]
  exact parallel_root_multiplicities rho mu sigma PYLoci.parallelPoint _hr _hm _hs _hs1 hk hq

theorem distinct_count (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma _z _k) n5 = ((distinct rho mu sigma _z _k).toFinset.card : ℝ) := by
  classical
  rw [List.toFinset_card_of_nodup (distinct_nodup rho mu sigma _z _k _hr _hm _hs _hs1),
    q8_stratum1_q3_root_count_cell0 rho mu sigma _z _k, distinct_reference rho mu sigma _z _k _hr (ne_of_gt _hs)]
  norm_num

end PYRootsParallel

namespace PYRootsPerpendicular

def candidatesTrees : List (Expr Symbol) := [n0, n3, n7]
def candidates (rho mu sigma _z : ℝ) (_k : Vec 3) : List ℝ := candidatesTrees.map (Expr.eval (values rho mu sigma _z _k))

theorem candidates_reference (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) : candidates rho mu sigma _z _k = [referenceRoot rho mu sigma PYLoci.perpendicularPoint 0, referenceRoot rho mu sigma PYLoci.perpendicularPoint 1, referenceRoot rho mu sigma PYLoci.perpendicularPoint 2] := by
  unfold candidates candidatesTrees
  simp only [List.map_cons, List.map_nil]
  rw [n0_eval rho mu sigma _z _k, q8_stratum2_q3_root_solutions_raw_cell1 rho mu sigma _z _k _hr, q8_stratum2_q3_root_solutions_raw_cell2 rho mu sigma _z _k _hr _hs]
  rfl

def distinctTrees : List (Expr Symbol) := [n4, n3, n7]
def distinct (rho mu sigma _z : ℝ) (_k : Vec 3) : List ℝ := distinctTrees.map (Expr.eval (values rho mu sigma _z _k))

theorem distinct_reference (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) : distinct rho mu sigma _z _k = [referenceRoot rho mu sigma PYLoci.perpendicularPoint 0, referenceRoot rho mu sigma PYLoci.perpendicularPoint 1, referenceRoot rho mu sigma PYLoci.perpendicularPoint 2] := by
  unfold distinct distinctTrees
  simp only [List.map_cons, List.map_nil]
  rw [q8_stratum2_q3_roots_distinct_cell0 rho mu sigma _z _k, q8_stratum2_q3_roots_distinct_cell1 rho mu sigma _z _k _hr, q8_stratum2_q3_roots_distinct_cell2 rho mu sigma _z _k _hr _hs]

theorem candidates_complete (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) (w : ℝ) :
    Matrix.det ((1/2 : ℝ) • referenceMatrix rho mu sigma w PYLoci.perpendicularPoint) = 0 ↔ w ∈ candidates rho mu sigma _z _k := by
  rw [rootPolynomial_complete rho mu sigma w PYLoci.perpendicularPoint _hr (ne_of_gt _hs),
    rootPolynomial_roots rho mu sigma PYLoci.perpendicularPoint _hr (ne_of_gt _hs), candidates_reference rho mu sigma _z _k _hr (ne_of_gt _hs)]
  simp [referenceRoot]

theorem distinct_complete (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) (w : ℝ) :
    Matrix.det ((1/2 : ℝ) • referenceMatrix rho mu sigma w PYLoci.perpendicularPoint) = 0 ↔ w ∈ distinct rho mu sigma _z _k := by
  rw [rootPolynomial_complete rho mu sigma w PYLoci.perpendicularPoint _hr (ne_of_gt _hs),
    rootPolynomial_roots rho mu sigma PYLoci.perpendicularPoint _hr (ne_of_gt _hs), distinct_reference rho mu sigma _z _k _hr (ne_of_gt _hs)]
  simp [referenceRoot]

theorem emitted_determinant_complete (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) (w : ℝ) :
    Expr.eval (values rho mu sigma w _k) PYRerun.n32 = 0 ↔ w ∈ distinct rho mu sigma _z _k := by
  rw [PYRerun.q8_stratum2_q3_determinant_cell0 rho mu sigma w _k]
  exact distinct_complete rho mu sigma _z _k _hr _hm _hs _hs1 w

theorem distinct_nodup (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) : (distinct rho mu sigma _z _k).Nodup := by
  have hk : PYLoci.perpendicularPoint ≠ 0 := PYLoci.perpendicularPoint_nonzero
  have hq : perpSq 0 PYLoci.perpendicularPoint ≠ 0 := by
    norm_num [PYLoci.perpendicularPoint, perpSq, normSq, dot, Fin.sum_univ_three, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]
  rw [distinct_reference rho mu sigma _z _k _hr (ne_of_gt _hs)]
  exact split_root_nodup rho mu sigma PYLoci.perpendicularPoint _hr _hm _hs _hs1 hk hq

theorem multiplicities (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    rootMultiplicity 0 (rootPolynomial rho mu sigma PYLoci.perpendicularPoint) = 1 ∧
      rootMultiplicity (referenceRoot rho mu sigma PYLoci.perpendicularPoint 1) (rootPolynomial rho mu sigma PYLoci.perpendicularPoint) = 1 ∧
      rootMultiplicity (referenceRoot rho mu sigma PYLoci.perpendicularPoint 2) (rootPolynomial rho mu sigma PYLoci.perpendicularPoint) = 1
    := by
  have hk : PYLoci.perpendicularPoint ≠ 0 := PYLoci.perpendicularPoint_nonzero
  have hq : perpSq 0 PYLoci.perpendicularPoint ≠ 0 := by
    norm_num [PYLoci.perpendicularPoint, perpSq, normSq, dot, Fin.sum_univ_three, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]
  exact split_root_multiplicities rho mu sigma PYLoci.perpendicularPoint _hr _hm _hs _hs1 hk hq

theorem candidate_multiset (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    (rootPolynomial rho mu sigma PYLoci.perpendicularPoint).roots = (candidates rho mu sigma _z _k : Multiset ℝ) := by
  rw [rootPolynomial_roots rho mu sigma PYLoci.perpendicularPoint _hr (ne_of_gt _hs), candidates_reference rho mu sigma _z _k _hr (ne_of_gt _hs)]
  rfl

theorem distinct_count (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma _z _k) n8 = ((distinct rho mu sigma _z _k).toFinset.card : ℝ) := by
  classical
  rw [List.toFinset_card_of_nodup (distinct_nodup rho mu sigma _z _k _hr _hm _hs _hs1),
    q8_stratum2_q3_root_count_cell0 rho mu sigma _z _k, distinct_reference rho mu sigma _z _k _hr (ne_of_gt _hs)]
  norm_num

end PYRootsPerpendicular

namespace WLRootsGeneric

def candidatesTrees : List (Expr Symbol) := [n0, n17, n24]
def candidates (rho mu sigma _z : ℝ) (_k : Vec 3) : List ℝ := candidatesTrees.map (Expr.eval (values rho mu sigma _z _k))

theorem candidates_reference (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) : candidates rho mu sigma _z _k = [referenceRoot rho mu sigma _k 0, referenceRoot rho mu sigma _k 1, referenceRoot rho mu sigma _k 2] := by
  unfold candidates candidatesTrees
  simp only [List.map_cons, List.map_nil]
  rw [n0_eval rho mu sigma _z _k, q3_solutions_cell1 rho mu sigma _z _k _hr, q3_solutions_cell2 rho mu sigma _z _k _hr _hs]
  rfl

def distinctTrees : List (Expr Symbol) := [n5, n17, n24]
def distinct (rho mu sigma _z : ℝ) (_k : Vec 3) : List ℝ := distinctTrees.map (Expr.eval (values rho mu sigma _z _k))

theorem distinct_reference (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) : distinct rho mu sigma _z _k = [referenceRoot rho mu sigma _k 0, referenceRoot rho mu sigma _k 1, referenceRoot rho mu sigma _k 2] := by
  unfold distinct distinctTrees
  simp only [List.map_cons, List.map_nil]
  rw [q3_distinct_roots_cell0 rho mu sigma _z _k, q3_distinct_roots_cell1 rho mu sigma _z _k _hr, q3_distinct_roots_cell2 rho mu sigma _z _k _hr _hs]

theorem candidates_complete (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) (_hg : GenericChart sigma _k) (w : ℝ) :
    Matrix.det ((1/2 : ℝ) • referenceMatrix rho mu sigma w _k) = 0 ↔ w ∈ candidates rho mu sigma _z _k := by
  rw [rootPolynomial_complete rho mu sigma w _k _hr (ne_of_gt _hs),
    rootPolynomial_roots rho mu sigma _k _hr (ne_of_gt _hs), candidates_reference rho mu sigma _z _k _hr (ne_of_gt _hs)]
  simp [referenceRoot]

theorem distinct_complete (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) (_hg : GenericChart sigma _k) (w : ℝ) :
    Matrix.det ((1/2 : ℝ) • referenceMatrix rho mu sigma w _k) = 0 ↔ w ∈ distinct rho mu sigma _z _k := by
  rw [rootPolynomial_complete rho mu sigma w _k _hr (ne_of_gt _hs),
    rootPolynomial_roots rho mu sigma _k _hr (ne_of_gt _hs), distinct_reference rho mu sigma _z _k _hr (ne_of_gt _hs)]
  simp [referenceRoot]

theorem emitted_determinant_complete (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) (_hg : GenericChart sigma _k) (w : ℝ) :
    Expr.eval (values rho mu sigma w _k) WL.n77 = 0 ↔ w ∈ distinct rho mu sigma _z _k := by
  rw [WL.q3_determinant_cell0 rho mu sigma w _k]
  exact distinct_complete rho mu sigma _z _k _hr _hm _hs _hs1 _hg w

theorem distinct_nodup (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) (_hg : GenericChart sigma _k) : (distinct rho mu sigma _z _k).Nodup := by
  have hk : _k ≠ 0 := (genericChart_split sigma _k _hg).1
  have hq : perpSq 0 _k ≠ 0 := (genericChart_split sigma _k _hg).2
  rw [distinct_reference rho mu sigma _z _k _hr (ne_of_gt _hs)]
  exact split_root_nodup rho mu sigma _k _hr _hm _hs _hs1 hk hq

theorem multiplicities (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) (_hg : GenericChart sigma _k) :
    rootMultiplicity 0 (rootPolynomial rho mu sigma _k) = 1 ∧
      rootMultiplicity (referenceRoot rho mu sigma _k 1) (rootPolynomial rho mu sigma _k) = 1 ∧
      rootMultiplicity (referenceRoot rho mu sigma _k 2) (rootPolynomial rho mu sigma _k) = 1
    := by
  have hk : _k ≠ 0 := (genericChart_split sigma _k _hg).1
  have hq : perpSq 0 _k ≠ 0 := (genericChart_split sigma _k _hg).2
  exact split_root_multiplicities rho mu sigma _k _hr _hm _hs _hs1 hk hq

theorem candidate_multiset (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) (_hg : GenericChart sigma _k) :
    (rootPolynomial rho mu sigma _k).roots = (candidates rho mu sigma _z _k : Multiset ℝ) := by
  rw [rootPolynomial_roots rho mu sigma _k _hr (ne_of_gt _hs), candidates_reference rho mu sigma _z _k _hr (ne_of_gt _hs)]
  rfl

theorem distinct_count (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) (_hg : GenericChart sigma _k) :
    Expr.eval (values rho mu sigma _z _k) n25 = ((distinct rho mu sigma _z _k).toFinset.card : ℝ) := by
  classical
  rw [List.toFinset_card_of_nodup (distinct_nodup rho mu sigma _z _k _hr _hm _hs _hs1 _hg),
    q3_root_count_cell0 rho mu sigma _z _k, distinct_reference rho mu sigma _z _k _hr (ne_of_gt _hs)]
  norm_num

def filteredTrees : List (Expr Symbol) := candidatesTrees.filter frequencyFree
def filtered (rho mu sigma _z : ℝ) (_k : Vec 3) : List ℝ := filteredTrees.map (Expr.eval (values rho mu sigma _z _k))

theorem filter_keeps_all : filteredTrees = candidatesTrees := by rfl

theorem discarded_empty : candidatesTrees.filter (fun e => !frequencyFree e) = [] := by rfl

theorem filtered_reference (rho mu sigma _z : ℝ) (_k : Vec 3) : filtered rho mu sigma _z _k = candidates rho mu sigma _z _k := by
  rw [filtered, filter_keeps_all]; rfl

theorem candidate_count (rho mu sigma _z : ℝ) (_k : Vec 3) :
    Expr.eval (values rho mu sigma _z _k) n25 = ((candidates rho mu sigma _z _k).length : ℝ) := by
  rw [q3_root_candidate_count_before_filter_cell0 rho mu sigma _z _k]
  rfl

theorem filtered_count (rho mu sigma _z : ℝ) (_k : Vec 3) :
    Expr.eval (values rho mu sigma _z _k) n25 = ((filtered rho mu sigma _z _k).length : ℝ) := by
  rw [q3_root_candidate_count_after_filter_cell0 rho mu sigma _z _k]
  rw [filtered_reference]
  rfl

theorem list_counts (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) (_hg : GenericChart sigma _k) :
    Expr.eval (values rho mu sigma _z _k) n25 = ((filtered rho mu sigma _z _k).length : ℝ) ∧
      Expr.eval (values rho mu sigma _z _k) n25 = ((distinct rho mu sigma _z _k).toFinset.card : ℝ) := by
  classical
  constructor
  · rw [q3_root_list_counts_cell0 rho mu sigma _z _k, filtered_reference]; rfl
  · rw [List.toFinset_card_of_nodup (distinct_nodup rho mu sigma _z _k _hr _hm _hs _hs1 _hg),
      q3_root_list_counts_cell1 rho mu sigma _z _k, distinct_reference rho mu sigma _z _k _hr (ne_of_gt _hs)]
    norm_num

end WLRootsGeneric

namespace WLRootsParallel

def candidatesTrees : List (Expr Symbol) := [n0, n7, n7]
def candidates (rho mu sigma _z : ℝ) (_k : Vec 3) : List ℝ := candidatesTrees.map (Expr.eval (values rho mu sigma _z _k))

theorem candidates_reference (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) : candidates rho mu sigma _z _k = [referenceRoot rho mu sigma WLLoci.parallelPoint 0, referenceRoot rho mu sigma WLLoci.parallelPoint 1, referenceRoot rho mu sigma WLLoci.parallelPoint 1] := by
  unfold candidates candidatesTrees
  simp only [List.map_cons, List.map_nil]
  rw [n0_eval rho mu sigma _z _k, stratum1_q3_solutions_cell2 rho mu sigma _z _k _hr]
  rfl

def distinctTrees : List (Expr Symbol) := [n4, n7]
def distinct (rho mu sigma _z : ℝ) (_k : Vec 3) : List ℝ := distinctTrees.map (Expr.eval (values rho mu sigma _z _k))

theorem distinct_reference (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) : distinct rho mu sigma _z _k = [referenceRoot rho mu sigma WLLoci.parallelPoint 0, referenceRoot rho mu sigma WLLoci.parallelPoint 1] := by
  unfold distinct distinctTrees
  simp only [List.map_cons, List.map_nil]
  rw [stratum1_q3_distinct_roots_cell0 rho mu sigma _z _k, stratum1_q3_distinct_roots_cell1 rho mu sigma _z _k _hr]

theorem candidates_complete (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) (w : ℝ) :
    Matrix.det ((1/2 : ℝ) • referenceMatrix rho mu sigma w WLLoci.parallelPoint) = 0 ↔ w ∈ candidates rho mu sigma _z _k := by
  have hq : perpSq 0 WLLoci.parallelPoint = 0 := by
    norm_num [WLLoci.parallelPoint, perpSq, normSq, dot, Fin.sum_univ_three, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]
  have he := (root_coincidence rho mu sigma WLLoci.parallelPoint _hr _hm (ne_of_gt _hs) _hs1).mpr hq
  rw [rootPolynomial_complete rho mu sigma w WLLoci.parallelPoint _hr (ne_of_gt _hs),
    rootPolynomial_roots rho mu sigma WLLoci.parallelPoint _hr (ne_of_gt _hs), candidates_reference rho mu sigma _z _k _hr (ne_of_gt _hs)]
  rw [he]
  simp [referenceRoot]

theorem distinct_complete (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) (w : ℝ) :
    Matrix.det ((1/2 : ℝ) • referenceMatrix rho mu sigma w WLLoci.parallelPoint) = 0 ↔ w ∈ distinct rho mu sigma _z _k := by
  have hq : perpSq 0 WLLoci.parallelPoint = 0 := by
    norm_num [WLLoci.parallelPoint, perpSq, normSq, dot, Fin.sum_univ_three, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]
  have he := (root_coincidence rho mu sigma WLLoci.parallelPoint _hr _hm (ne_of_gt _hs) _hs1).mpr hq
  rw [rootPolynomial_complete rho mu sigma w WLLoci.parallelPoint _hr (ne_of_gt _hs),
    rootPolynomial_roots rho mu sigma WLLoci.parallelPoint _hr (ne_of_gt _hs), distinct_reference rho mu sigma _z _k _hr (ne_of_gt _hs)]
  rw [he]
  simp [referenceRoot]

theorem emitted_determinant_complete (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) (w : ℝ) :
    Expr.eval (values rho mu sigma w _k) WLRerun.n12 = 0 ↔ w ∈ distinct rho mu sigma _z _k := by
  rw [WLRerun.stratum1_q3_determinant_cell0 rho mu sigma w _k]
  exact distinct_complete rho mu sigma _z _k _hr _hm _hs _hs1 w

theorem distinct_nodup (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) : (distinct rho mu sigma _z _k).Nodup := by
  have hk : WLLoci.parallelPoint ≠ 0 := WLLoci.parallelPoint_nonzero
  have hq : perpSq 0 WLLoci.parallelPoint = 0 := by
    norm_num [WLLoci.parallelPoint, perpSq, normSq, dot, Fin.sum_univ_three, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]
  rw [distinct_reference rho mu sigma _z _k _hr (ne_of_gt _hs)]
  exact parallel_root_nodup rho mu sigma WLLoci.parallelPoint _hr _hm _hs hk

theorem multiplicities (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    rootMultiplicity 0 (rootPolynomial rho mu sigma WLLoci.parallelPoint) = 1 ∧
      rootMultiplicity (referenceRoot rho mu sigma WLLoci.parallelPoint 1) (rootPolynomial rho mu sigma WLLoci.parallelPoint) = 2
    := by
  have hk : WLLoci.parallelPoint ≠ 0 := WLLoci.parallelPoint_nonzero
  have hq : perpSq 0 WLLoci.parallelPoint = 0 := by
    norm_num [WLLoci.parallelPoint, perpSq, normSq, dot, Fin.sum_univ_three, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]
  exact parallel_root_multiplicities rho mu sigma WLLoci.parallelPoint _hr _hm _hs _hs1 hk hq

theorem candidate_multiset (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    (rootPolynomial rho mu sigma WLLoci.parallelPoint).roots = (candidates rho mu sigma _z _k : Multiset ℝ) := by
  have hq : perpSq 0 WLLoci.parallelPoint = 0 := by
    norm_num [WLLoci.parallelPoint, perpSq, normSq, dot, Fin.sum_univ_three, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]
  have he := (root_coincidence rho mu sigma WLLoci.parallelPoint _hr _hm (ne_of_gt _hs) _hs1).mpr hq
  rw [rootPolynomial_roots rho mu sigma WLLoci.parallelPoint _hr (ne_of_gt _hs), candidates_reference rho mu sigma _z _k _hr (ne_of_gt _hs)]
  rw [he]
  rfl

theorem distinct_count (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma _z _k) n8 = ((distinct rho mu sigma _z _k).toFinset.card : ℝ) := by
  classical
  rw [List.toFinset_card_of_nodup (distinct_nodup rho mu sigma _z _k _hr _hm _hs _hs1),
    stratum1_q3_root_count_cell0 rho mu sigma _z _k, distinct_reference rho mu sigma _z _k _hr (ne_of_gt _hs)]
  norm_num

def filteredTrees : List (Expr Symbol) := candidatesTrees.filter frequencyFree
def filtered (rho mu sigma _z : ℝ) (_k : Vec 3) : List ℝ := filteredTrees.map (Expr.eval (values rho mu sigma _z _k))

theorem filter_keeps_all : filteredTrees = candidatesTrees := by rfl

theorem discarded_empty : candidatesTrees.filter (fun e => !frequencyFree e) = [] := by rfl

theorem filtered_reference (rho mu sigma _z : ℝ) (_k : Vec 3) : filtered rho mu sigma _z _k = candidates rho mu sigma _z _k := by
  rw [filtered, filter_keeps_all]; rfl

theorem candidate_count (rho mu sigma _z : ℝ) (_k : Vec 3) :
    Expr.eval (values rho mu sigma _z _k) n9 = ((candidates rho mu sigma _z _k).length : ℝ) := by
  rw [stratum1_q3_root_candidate_count_before_filter_cell0 rho mu sigma _z _k]
  rfl

theorem filtered_count (rho mu sigma _z : ℝ) (_k : Vec 3) :
    Expr.eval (values rho mu sigma _z _k) n9 = ((filtered rho mu sigma _z _k).length : ℝ) := by
  rw [stratum1_q3_root_candidate_count_after_filter_cell0 rho mu sigma _z _k]
  rw [filtered_reference]
  rfl

theorem list_counts (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma _z _k) n9 = ((filtered rho mu sigma _z _k).length : ℝ) ∧
      Expr.eval (values rho mu sigma _z _k) n8 = ((distinct rho mu sigma _z _k).toFinset.card : ℝ) := by
  classical
  constructor
  · rw [stratum1_q3_root_list_counts_cell0 rho mu sigma _z _k, filtered_reference]; rfl
  · rw [List.toFinset_card_of_nodup (distinct_nodup rho mu sigma _z _k _hr _hm _hs _hs1),
      stratum1_q3_root_list_counts_cell1 rho mu sigma _z _k, distinct_reference rho mu sigma _z _k _hr (ne_of_gt _hs)]
    norm_num

end WLRootsParallel

namespace WLRootsPerpendicular

def candidatesTrees : List (Expr Symbol) := [n0, n9, n12]
def candidates (rho mu sigma _z : ℝ) (_k : Vec 3) : List ℝ := candidatesTrees.map (Expr.eval (values rho mu sigma _z _k))

theorem candidates_reference (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) : candidates rho mu sigma _z _k = [referenceRoot rho mu sigma WLLoci.perpendicularPoint 0, referenceRoot rho mu sigma WLLoci.perpendicularPoint 1, referenceRoot rho mu sigma WLLoci.perpendicularPoint 2] := by
  unfold candidates candidatesTrees
  simp only [List.map_cons, List.map_nil]
  rw [n0_eval rho mu sigma _z _k, stratum2_q3_solutions_cell1 rho mu sigma _z _k _hr, stratum2_q3_solutions_cell2 rho mu sigma _z _k _hr _hs]
  rfl

def distinctTrees : List (Expr Symbol) := [n4, n9, n12]
def distinct (rho mu sigma _z : ℝ) (_k : Vec 3) : List ℝ := distinctTrees.map (Expr.eval (values rho mu sigma _z _k))

theorem distinct_reference (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) : distinct rho mu sigma _z _k = [referenceRoot rho mu sigma WLLoci.perpendicularPoint 0, referenceRoot rho mu sigma WLLoci.perpendicularPoint 1, referenceRoot rho mu sigma WLLoci.perpendicularPoint 2] := by
  unfold distinct distinctTrees
  simp only [List.map_cons, List.map_nil]
  rw [stratum2_q3_distinct_roots_cell0 rho mu sigma _z _k, stratum2_q3_distinct_roots_cell1 rho mu sigma _z _k _hr, stratum2_q3_distinct_roots_cell2 rho mu sigma _z _k _hr _hs]

theorem candidates_complete (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) (w : ℝ) :
    Matrix.det ((1/2 : ℝ) • referenceMatrix rho mu sigma w WLLoci.perpendicularPoint) = 0 ↔ w ∈ candidates rho mu sigma _z _k := by
  rw [rootPolynomial_complete rho mu sigma w WLLoci.perpendicularPoint _hr (ne_of_gt _hs),
    rootPolynomial_roots rho mu sigma WLLoci.perpendicularPoint _hr (ne_of_gt _hs), candidates_reference rho mu sigma _z _k _hr (ne_of_gt _hs)]
  simp [referenceRoot]

theorem distinct_complete (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) (w : ℝ) :
    Matrix.det ((1/2 : ℝ) • referenceMatrix rho mu sigma w WLLoci.perpendicularPoint) = 0 ↔ w ∈ distinct rho mu sigma _z _k := by
  rw [rootPolynomial_complete rho mu sigma w WLLoci.perpendicularPoint _hr (ne_of_gt _hs),
    rootPolynomial_roots rho mu sigma WLLoci.perpendicularPoint _hr (ne_of_gt _hs), distinct_reference rho mu sigma _z _k _hr (ne_of_gt _hs)]
  simp [referenceRoot]

theorem emitted_determinant_complete (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) (w : ℝ) :
    Expr.eval (values rho mu sigma w _k) WLRerun.n28 = 0 ↔ w ∈ distinct rho mu sigma _z _k := by
  rw [WLRerun.stratum2_q3_determinant_cell0 rho mu sigma w _k]
  exact distinct_complete rho mu sigma _z _k _hr _hm _hs _hs1 w

theorem distinct_nodup (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) : (distinct rho mu sigma _z _k).Nodup := by
  have hk : WLLoci.perpendicularPoint ≠ 0 := WLLoci.perpendicularPoint_nonzero
  have hq : perpSq 0 WLLoci.perpendicularPoint ≠ 0 := by
    norm_num [WLLoci.perpendicularPoint, perpSq, normSq, dot, Fin.sum_univ_three, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]
  rw [distinct_reference rho mu sigma _z _k _hr (ne_of_gt _hs)]
  exact split_root_nodup rho mu sigma WLLoci.perpendicularPoint _hr _hm _hs _hs1 hk hq

theorem multiplicities (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    rootMultiplicity 0 (rootPolynomial rho mu sigma WLLoci.perpendicularPoint) = 1 ∧
      rootMultiplicity (referenceRoot rho mu sigma WLLoci.perpendicularPoint 1) (rootPolynomial rho mu sigma WLLoci.perpendicularPoint) = 1 ∧
      rootMultiplicity (referenceRoot rho mu sigma WLLoci.perpendicularPoint 2) (rootPolynomial rho mu sigma WLLoci.perpendicularPoint) = 1
    := by
  have hk : WLLoci.perpendicularPoint ≠ 0 := WLLoci.perpendicularPoint_nonzero
  have hq : perpSq 0 WLLoci.perpendicularPoint ≠ 0 := by
    norm_num [WLLoci.perpendicularPoint, perpSq, normSq, dot, Fin.sum_univ_three, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]
  exact split_root_multiplicities rho mu sigma WLLoci.perpendicularPoint _hr _hm _hs _hs1 hk hq

theorem candidate_multiset (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    (rootPolynomial rho mu sigma WLLoci.perpendicularPoint).roots = (candidates rho mu sigma _z _k : Multiset ℝ) := by
  rw [rootPolynomial_roots rho mu sigma WLLoci.perpendicularPoint _hr (ne_of_gt _hs), candidates_reference rho mu sigma _z _k _hr (ne_of_gt _hs)]
  rfl

theorem distinct_count (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma _z _k) n13 = ((distinct rho mu sigma _z _k).toFinset.card : ℝ) := by
  classical
  rw [List.toFinset_card_of_nodup (distinct_nodup rho mu sigma _z _k _hr _hm _hs _hs1),
    stratum2_q3_root_count_cell0 rho mu sigma _z _k, distinct_reference rho mu sigma _z _k _hr (ne_of_gt _hs)]
  norm_num

def filteredTrees : List (Expr Symbol) := candidatesTrees.filter frequencyFree
def filtered (rho mu sigma _z : ℝ) (_k : Vec 3) : List ℝ := filteredTrees.map (Expr.eval (values rho mu sigma _z _k))

theorem filter_keeps_all : filteredTrees = candidatesTrees := by rfl

theorem discarded_empty : candidatesTrees.filter (fun e => !frequencyFree e) = [] := by rfl

theorem filtered_reference (rho mu sigma _z : ℝ) (_k : Vec 3) : filtered rho mu sigma _z _k = candidates rho mu sigma _z _k := by
  rw [filtered, filter_keeps_all]; rfl

theorem candidate_count (rho mu sigma _z : ℝ) (_k : Vec 3) :
    Expr.eval (values rho mu sigma _z _k) n13 = ((candidates rho mu sigma _z _k).length : ℝ) := by
  rw [stratum2_q3_root_candidate_count_before_filter_cell0 rho mu sigma _z _k]
  rfl

theorem filtered_count (rho mu sigma _z : ℝ) (_k : Vec 3) :
    Expr.eval (values rho mu sigma _z _k) n13 = ((filtered rho mu sigma _z _k).length : ℝ) := by
  rw [stratum2_q3_root_candidate_count_after_filter_cell0 rho mu sigma _z _k]
  rw [filtered_reference]
  rfl

theorem list_counts (rho mu sigma _z : ℝ) (_k : Vec 3) (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :
    Expr.eval (values rho mu sigma _z _k) n13 = ((filtered rho mu sigma _z _k).length : ℝ) ∧
      Expr.eval (values rho mu sigma _z _k) n13 = ((distinct rho mu sigma _z _k).toFinset.card : ℝ) := by
  classical
  constructor
  · rw [stratum2_q3_root_list_counts_cell0 rho mu sigma _z _k, filtered_reference]; rfl
  · rw [List.toFinset_card_of_nodup (distinct_nodup rho mu sigma _z _k _hr _hm _hs _hs1),
      stratum2_q3_root_list_counts_cell1 rho mu sigma _z _k, distinct_reference rho mu sigma _z _k _hr (ne_of_gt _hs)]
    norm_num

end WLRootsPerpendicular

end
end S10Audit.CAS
