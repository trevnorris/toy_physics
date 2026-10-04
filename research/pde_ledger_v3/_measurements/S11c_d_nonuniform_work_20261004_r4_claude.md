NEEDS REVISION FOR THIS EXPLICIT FORCE-PAIRING AND CENTRE WORK METHOD

The architecture is sound. The revisions below are narrow, but they change what the first build may treat as an operand, so they are substantive. This is a hand reading of source text only. Nothing was executed.

## What holds up

- **Force pairing.**
  - The views confirm that `BODY.E_W = e_force` and the face fourth component is `expand(face*e_force/W0)` (`native-background-force-representations.txt`). The published operator values agree: face +1 is `E_W/W0`, face −1 is `−E_W/W0`.
  - The directive's "ordered pairing (bulk-DOF body force + per-face traction)" (`native-background-pairing-directive.txt`) describes a bundle format. It does not say the two are independent powers.
  - The support symbols `f_hold_e_W_0` and `t_hold_*_0_4` are free, with no equality supplied.
  - The plan keeps all of this visible. It emits separate operands and residuals, refuses to sum them, and allows `ORIGINAL_FORCE_PAIRING_UNRESOLVED`. That is the correct treatment, and it does not reinterpret the directive.
- **Rate map.**
  - The plan's `hdot_s` matches `build_face_source`: `zeta_t + s*u_t·(scales*grad_W)/2` for LAB_HELD, with `DELTA_W` giving `s*W0*e_W_t/2`.
  - The `F_face,U/e/c` coefficients follow from `t_s·(u_t, hdot_s)`.
  - The c2 reconstruction `(V − n_par·u_t)/normal[3]` is the exact inverse of the registry `face_velocity = ε·shape(n·v)`.
- **Legacy residual.** `object_difference` returns `sp.Equivalent(...)` for `Boolean` operands. Treating the residual as literal provenance is right.
- **Centre.**
  - The `R0_c` / `UNSUPPLIED` rule is right.
  - So is keeping ε and ε² separate from defect λ².
  - So is keeping the physical `ZETA_C` column and `CENTER_FACE_GENERALIZED_ROW` separate. `face_generalized_force_rows` takes the centre row from the `(DELTA_W physical, ZETA_C virtual)` case, as the plan says.
- **Controls and stop rules.** The controls are sensitivity checks with an "untested, not passed" rule, and unexplained mismatches stop the run.

## Substantive revisions

1. **Tag every pullback operand with its ε-order and name its source.** The plan writes `t_s` and `A_s` without an ε-order tag.
   - `traction_raw` returns `epsilon*shape(exact)`, and `face_velocity_raw` returns `epsilon*shape(...)`. Both are ε¹ increments. Only `face_normal_raw` and `face_measure_raw` carry an ε⁰ entry (`background`).
   - The a-registry therefore has no ε⁰ traction. Any ε⁰ face traction is either `admissibility_operator_operand`, which has only a fourth component, or the background pressure embedded in the `virtual_work_shape_deriv` density.
   - The plan's phrase "background registry traction" and its `TRACTION_MEASURE_UNSUPPLIED` clause should say this. Otherwise a worker could treat an ε¹ increment as `R0_c`'s traction.
   - The expected default outcome for the b operator's face traction should be stated as `TRACTION_MEASURE_UNSUPPLIED`, and so `R0_c` as `UNSUPPLIED`, unless a source join appears.
   - Because `A_s` is ε⁰ (measure at parameter 0), also pre-register the ε² pieces: `A^(1)·t^(0)·rate^(1)`, `t^(1)·A^(1)`, and the rate-map correction. The R0 ledger currently only names `R0*zeta_t^(1)` and `R0*zeta_t^(2)`.

2. **Pre-register the c2 area fact.**
   - `normal_exact` is a unit 4-normal by construction: `(−s·∇h, s)/√(1+|∇h|²)`. The c2 expression `sqrt(Σ normal_i²)` is therefore identically 1.
   - c2 also takes `face_normal[...][0]`, the parameter-0 entry. The a measure is `measure_exact = √(1+|∇h|²)`. At the nonuniform background, `∇h0 = s∇W/2` is nonzero, so the two differ at order |∇W|²/4.
   - The "compare actual normalization" step should record this explicitly: area ≡ 1 versus `D_s`. It should also record that c2 drops the measure factor. Otherwise it reads as a numeric comparison to be discovered.
   - That difference is 20/02-grade debt and must stay unprojected, as the plan already requires.
   - From the source, the registry traction `−(P+Λ_X χ)·n_exact` with a unit normal is per true area. The plan's hedge about "per true or projected area" can be stated as determined for the a-registry. It remains open only for the b admissibility traction.

3. **State the expected U-row and centre consequences of the face representation.**
   - Applying the plan's own `F_face,U = Σ A(t_par + s·t_4·∇W/2)` to `t_4 = s·e_force/W0` gives `A·e_force·∇W/W0`. The operator's body `U` is literal `(0,0,0)`. So the face representation cannot equal the body representation at non-flat grade. The first build should expect and record that mismatch, not treat flat agreement as the pairing.
   - The vanishing `F_face,c = Σ t_4,s` follows from the antisymmetric assignment `t_s ∝ s`. The energy variation fixes only `Σ_s s·t_s·W0/2 = e_force`. The symmetric (centre) part is not determined by `e_force`. The plan already says a vanishing face piece does not prove the total. It should add that this zero is a property of the source assignment, not independent evidence.

4. **Clarify the support control.** `F_support,c` for the bundle depends on free `t_hold` symbols. The "centre-support control" should mutate the nonzero operator face component (live `e_force`) or a symbol coefficient. It should not report "resolved face piece" movement for a quantity that is symbolic.

## Optional wording

- Say whether `normal[3]` means the fourth component or index 3 of the tuple, and whether `[0]` means the background entry.
- Say once that the cancellation of `V_s`'s slope terms is a hand check against saved tables. By my hand algebra the tangential slope terms cancel in `V = n0·(u_t, hdot)`, giving `(s·ζ_t + W0·e_t/2)/D_s`. That is not verified here.

## Can the screen decide the question?

Yes, within limits. It can expose a work obstruction: the U-row mismatch, the area ≡ 1 versus `D_s` difference, or an unjoined traction area convention. It can also give a conditional consistent pairing. It needs no agency, centre constraint, grade or area normalization to do so. `ORIGINAL_FORCE_PAIRING_UNRESOLVED` and `UNSUPPLIED` are legitimate results. No real user choice has yet been shown. A choice would be needed only after a concrete residual is identified, such as the support-symbol compatibility residual. Q_port, physical total transverse survival including reflection, the 20/02 and second-order end-current duties, and the threshold remainder remain separate and unaddressed by this method.

## Reading coverage and uncertainty

I read in full:
- `guide.txt`, `plan.txt` and `review-prompt.md`.
- The three force/pairing/variation directive views.
- `native-c2-full-traction-pairing.txt`.
- `native-full-virtual-work.txt`, `native-full-face-source.txt` and `native-normal-shape.txt`.
- `native-true-area-shape.txt`, `native-outward-velocity.txt` and `native-traction-shape.txt`.
- `native-background-support-definition.txt`, `native-legacy-residual-definition.txt` and `native-centre-and-thickness-fields.txt`.
- `native-face-work-row-definition.txt`.
- The published operator-operand and support-operand JSON (all four cases, literal display and constructor).

I read only the header of `published-admissibility_residual.json`. Grep confirms two `Equivalent`/`Not` matches in it.

I did not read:
- The bodies of the `published-a-*` registries. Their single long lines were not displayable, so I relied on the source functions.
- `packet-index.json`, the pin joins, the 400 cells, the kinetic/energy/chemical views, or the `source/` JSONs.
- The `original/` sources, which I only confirmed exist.

The `native-slab-work.json` view has no centre token. The centre row is present in `native-face-work-row-definition.txt` and the original b sources, so I treated it as supplied.

The claims about ε-order, unit normal and area follow from reading function text. I did not check them against the stored tables, so they could be wrong where the registries differ from the functions.

This is a method verdict only. It is not a worker, a result, a support law or a leakage acceptance.