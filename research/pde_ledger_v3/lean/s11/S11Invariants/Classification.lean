import S11Invariants.Rotation

/-! Exhaustive D2 classification. Finite test rotations are used only for
necessary constraints; sufficiency is proved for every proper rotation. -/
namespace S11Invariants
noncomputable section

/-- Coefficients of trace-square, skew-square, traceless-symmetric norm,
and the trace/skew pairing, in that order. -/
def invariantForm (v : Fin 4 → ℝ) : Quad :=
  v 0 • monomial 0 0 + v 1 • monomial 1 1 +
  v 2 • (monomial 2 2 + monomial 3 3) + v 3 • monomial 0 1

theorem invariantForm_apply (v : Fin 4 → ℝ) (G : Mat) :
    invariantForm v G = v 0 * coordinates G 0 ^ 2 + v 1 * coordinates G 1 ^ 2 +
      v 2 * (coordinates G 2 ^ 2 + coordinates G 3 ^ 2) +
      v 3 * coordinates G 0 * coordinates G 1 := by
  simp [invariantForm, monomial, coordinateMap, pow_two]
  ring

theorem invariantForm_SO (v : Fin 4 → ℝ) : SOInvariant (invariantForm v) := by
  intro R hR G
  obtain ⟨a,b,hn,rfl⟩ := proper_rotation hR
  simp only [invariantForm_apply, coordinates_rotation, Matrix.cons_val]
  rw [spin_two_norm, hn]
  ring

/-- Necessary conditions extracted from two admissible rotations. -/
theorem invariant_coefficients {Q : Quad} (hQ : SOInvariant Q) (c : Coefficients)
    (hc : ∀ G, Q G = polynomial c (coordinates G)) :
    c 2 = 0 ∧ c 3 = 0 ∧ c 5 = 0 ∧ c 6 = 0 ∧ c 7 = c 9 ∧ c 8 = 0 := by
  have quarter (v : Coordinates) := hQ (rotation 0 1) (rotation_proper (by norm_num)) (decode v)
  simp only [hc, coordinates_rotation, coordinates_decode] at quarter
  have h2 := quarter ![1,0,1,0]
  have h3 := quarter ![1,0,0,1]
  have h5 := quarter ![0,1,1,0]
  have h6 := quarter ![0,1,0,1]
  norm_num [polynomial, Matrix.cons_val_two, Matrix.cons_val_three] at h2 h3 h5 h6
  have hc2 : c 2 = 0 := by linarith
  have hc3 : c 3 = 0 := by linarith
  have hc5 : c 5 = 0 := by linarith
  have hc6 : c 6 = 0 := by linarith
  have rational (v : Coordinates) :=
    hQ (rotation (3/5) (4/5)) (rotation_proper (by norm_num)) (decode v)
  simp only [hc, coordinates_rotation, coordinates_decode] at rational
  have hx := rational ![0,0,1,0]
  have hxy := rational ![0,0,1,1]
  norm_num [polynomial, Matrix.cons_val_two, Matrix.cons_val_three] at hx hxy
  exact ⟨hc2,hc3,hc5,hc6,by linarith,by linarith⟩

theorem SO_classification (Q : Quad) :
    SOInvariant Q ↔ ∃ v : Fin 4 → ℝ, Q = invariantForm v := by
  constructor
  · intro hQ
    obtain ⟨c,hc⟩ := quadratic_representation Q
    obtain ⟨h2,h3,h5,h6,h79,h8⟩ := invariant_coefficients hQ c hc
    refine ⟨![c 0,c 4,c 7,c 1], ?_⟩
    ext G
    rw [hc, invariantForm_apply]
    simp only [polynomial, Matrix.cons_val, h2,h3,h5,h6,h79,h8]
    ring
  · rintro ⟨v,rfl⟩
    exact invariantForm_SO v

theorem invariantForm_injective : Function.Injective invariantForm := by
  intro v w h
  have hv (z : Coordinates) := congrArg (fun Q : Quad => Q (decode z)) h
  simp only [invariantForm_apply, coordinates_decode] at hv
  have h0 := hv ![1,0,0,0]
  have h1 := hv ![0,1,0,0]
  have h2 := hv ![0,0,1,0]
  have h3 := hv ![1,1,0,0]
  norm_num [Matrix.cons_val_two, Matrix.cons_val_three] at h0 h1 h2 h3
  ext i
  fin_cases i
  · exact h0
  · exact h1
  · exact h2
  · change v 3 = w 3
    linarith

theorem invariantForm_O (v : Fin 4 → ℝ) : OInvariant (invariantForm v) ↔ v 3 = 0 := by
  rw [orthogonal_invariance_iff]
  constructor
  · rintro ⟨_,hf⟩
    have h := hf (decode ![1,1,0,0])
    norm_num [invariantForm_apply, coordinates_reflection, coordinates_decode,
      Matrix.cons_val_two, Matrix.cons_val_three] at h
    linarith
  · intro h
    refine ⟨invariantForm_SO v, fun G => ?_⟩
    simp only [invariantForm_apply, coordinates_reflection, Matrix.cons_val, h]
    ring

theorem invariantForm_odd (v : Fin 4 → ℝ) :
    ReflectionOdd (invariantForm v) ↔ v 0 = 0 ∧ v 1 = 0 ∧ v 2 = 0 := by
  constructor
  · intro hf
    have h0 := hf (decode ![1,0,0,0])
    have h1 := hf (decode ![0,1,0,0])
    have h2 := hf (decode ![0,0,1,0])
    norm_num [invariantForm_apply, coordinates_reflection, coordinates_decode,
      Matrix.cons_val_two, Matrix.cons_val_three] at h0 h1 h2
    exact ⟨by linarith,by linarith,by linarith⟩
  · rintro ⟨h0,h1,h2⟩ G
    simp only [invariantForm_apply, coordinates_reflection, Matrix.cons_val, h0,h1,h2]
    ring

theorem O_classification (Q : Quad) :
    OInvariant Q ↔ ∃ v : Fin 3 → ℝ, Q = invariantForm ![v 0,v 1,v 2,0] := by
  constructor
  · intro hQ
    obtain ⟨v,rfl⟩ := (SO_classification Q).mp ((orthogonal_invariance_iff Q).mp hQ).1
    have hv := (invariantForm_O v).mp hQ
    refine ⟨![v 0,v 1,v 2], ?_⟩
    congr 1
    ext i
    fin_cases i <;> simp [Matrix.cons_val_two, hv]
  · rintro ⟨v,rfl⟩
    exact (invariantForm_O _).mpr (by simp [Matrix.cons_val_three])

theorem odd_classification (Q : Quad) :
    SOInvariant Q ∧ ReflectionOdd Q ↔ ∃ d : ℝ, Q = invariantForm ![0,0,0,d] := by
  constructor
  · rintro ⟨hSO,hodd⟩
    obtain ⟨v,rfl⟩ := (SO_classification Q).mp hSO
    obtain ⟨h0,h1,h2⟩ := (invariantForm_odd v).mp hodd
    refine ⟨v 3, ?_⟩
    congr 1
    ext i
    fin_cases i <;> simp [h0,h1,h2]
  · rintro ⟨d,rfl⟩
    exact ⟨invariantForm_SO _, (invariantForm_odd _).mpr (by norm_num [Matrix.cons_val_two])⟩

end
end S11Invariants
