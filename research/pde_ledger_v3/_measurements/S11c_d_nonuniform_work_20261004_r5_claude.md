NEEDS REVISION FOR THIS WAVE-ORDER AND SAVED-GRADE CENTRE WORK METHOD

The architecture holds up against the stored operands and the native sources I read. The revisions are narrow. They are pre-registration and join items that would otherwise let a build mislabel an order or a representation. Nothing here is a result, a support law or a leakage value.

## Checked against source (hand reading; nothing executed)

**Operands.**
- **Normal.** The LAB_HELD `face_normal` VALUE[0] is `(-σw1_di/2 ×3, s)` for both signs and both DELTA_W and ZETA_C.
  - Its increments are `ε`-linear.
  - The DELTA_W fourth component carries `ε·σ(e·w1')` terms.
- **Measure.** `face_measure_shape_deriv` is `(1, ε·σ·W0 e·w1'/4, sum)` for DELTA_W and `(1, ±ε·σ·w1'·ζ'/2, …)` for ZETA_C. This matches the plan.
- **b bundle.** `BODY_FORCE.U` is literal `(0,0,0)`, and `PER_FACE_TRACTION[s]=(0,0,0,s·e_force/W0)`. I confirmed this from both the displayed values and the constructor at `S11c_b…py:4056–4070`.

**Face coefficients.** I re-derived them:
- `F_U = A·e_force·g/W0`
- `F_e = e_force`
- `F_c = 0`

They hold under the flat map and equal-and-opposite faces, as the plan says.

**Rate map.** With `n=(-g/2,s)`, `hdot = s(V+g·u/2) = ζ_t + sW0e_t/2 + s∇W·u/2`, as the plan says.

**Source-zero joins.**
- `physical_trace_fields` does set background pressure and velocity to zero.
- Background affinity is zero for a separate reason. `traction_raw` and `virtual_work_cases` build `exact_mu = parameter·μ_θ/ρ_br`, which vanishes at parameter 0. That is not `physical_trace_fields`.

**Wave-order ledger.** It is correct, including A1·t0·r1 at ε², A1·t1 as a force density, and `R0_c` staying unknown.

## Required revisions

1. **c2 area: pre-register the raw and retained values separately.**
   - `traction_pairing` computes `area=sqrt(sum(n*n for n in normal))`, and `normal` is the 4-vector `face_normal[...][0]`. The plan's `sum_i` should say it runs over all four components. Then `area = sqrt(1+Σg²/4)`, which is `D`.
   - The final `retained_shape` keeps only `η≤1, σ≤1`. In `shape_coefficients`, the Mul rule drops `σ·σ` and Pow expands about the truncated base. From that source reading, the retained c2 area at 00/10/01/11 is therefore likely 1, with `g²/8` as σ² debt.
   - The plan's phrase "not an identically unit area" is true only of the raw expression. Name both explicitly:
     - raw `D`
     - retained coefficient (expected 1, to be confirmed against actual records, not asserted)
     - named σ² area debt
   - The prompt forbids inferring "area ≡ 1" from the untruncated source. The same applies in the other direction: do not infer the retained value from the exact formula.

2. **Add a normalization join for `hdot`.**
   - The source `face_velocity_raw` is `ε·shape(normal_exact·v_exact)`, which uses the unit normal.
   - c2 then computes `(V − n_saved,par·u)/n_saved,4` with the non-unit saved normal, whose norm is `D`.
   - Since `D·n_exact = (-g/2, s)`, the exact identity gives `n_saved·r = D·V`, not `V`.
   - The difference is σ² (`V(1−D)`). It affects every pullback coefficient `F_U`, `F_e` and `F_c`, and must be recorded as discarded-order debt.
   - Pre-register that `n_saved[0]` coincides both with `D·n_exact` and with the σ¹ truncation of `n_exact`. The values alone cannot say which representation was intended.

3. **Make E_W grading explicit.**
   - Saved E_W has `−W0³κ_Wσ∇²w1/L` (σ¹) and `−2W0³ηκ_Wσ w1∇²w1/L` (ησ¹).
   - The published header is only `((0,0,0),)`, so it understates both.
   - `e_force` is built from `first_shape_series` (source line 4056). The face assignment is therefore a function of an already-truncated object.
   - The U candidate `A·e_force·g/W0` is σ² and ησ², outside the stored cap, against a literal-zero body U that is silent at those grades. That is order debt, not a retained mismatch.
   - Carry the full complement:
     - the retained-product coefficients
     - the unknown σ²-and-higher e_force terms
     - the D factor on A

4. **Handle two format hazards in the published residual and operator entries.**
   - **Face labels.** The residual display's `PER_FACE_TRACTION` labels are `0` and `0`, not `±1`, because they are label differences. Join faces by operator/support position and sign from the operator and support entries, never by the residual's labels.
   - **Zero-literal dimensions.** The operator's U slot and traction components 1–3 carry placeholder dimension `0` for literal zeros. The support carries `[-2,-2,1]` for those slots. A dimension check on those slots is vacuous, so do not count it as unit agreement.

5. **Cite the source of the affinity zero correctly.** The plan's "source-zero background pressure/affinity" over-attributes to `physical_trace_fields`, which zeroes only pressure and velocity. Cite the `parameter·μ_θ` construction for the affinity part.

## Optional wording

- Say that the plan is `plan.txt`, which stands in for plan2a.
- State the index range on `g_i` and `n_i` explicitly (three spatial entries, four with the fourth component).

## Reading coverage and uncertainty

I read:
- the review prompt, `guide.txt` and `plan.txt`
- the full `face_normal`, `face_measure`, operator and support displays
- the first 40 lines of the residual display
- the native views for c2 pairing, final retention, trace composition and `physical_trace_fields`
- the views for truncation, publication, area, normal, velocity, traction, virtual work and the case payload
- the background pairing and variation directives and the centre/thickness fields
- the b and c2 source excerpts for `shape_coefficients`, `multigrade` and `e_force`
- the first 60 lines of `packet-index.json`

I did not read:
- the `face_velocity`, traction, virtual-work or closure displays
- the chunk `.jsonl` files
- the JSON registries
- the receipt hashes
- the S11c-a build-face source or the 2b centre-row sources

The retained area of 1 is my reading of `shape_coefficients`, not a stored value. Section 2b (centre drive) and the `Q_port` and threshold material were not assessed beyond the plan's text.