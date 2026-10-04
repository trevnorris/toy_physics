CLEAR FOR THIS BACKGROUND FORCE-PAIRING BUILD AND AMENDED METHOD

This is a source and method read only. I ran nothing, built no symbolic objects, and computed no hashes. Nothing here is a runtime result, and the verdict does not promote any physical status.

**What I checked against the actual operands**
- **Operator and support trees.** Both are selected by (LAB_HELD, RHO4_CONSTANT).
  - The operator E_W is `-W0³κσ(1+2η·w1)·Σw_didi/L`, which is exactly what `scientific_work` asserts at `S11c_d_background_force_pairing.py:95`.
  - The operator face fourth components are `±e/W0` (b-operator source lines 4056–4073, line 97), with literal-zero tangential slots.
  - The operator DIMENSION indexing (`dims[0][1][0][1][j]`, `dims[1][1][j][1][k]`) matches the displayed nesting. Literal-zero U and tangential slots are correctly vacuous. The nonzero THETA, E_W and face-4 units agree between operator and support: (-1,-2,1), (-1,-2,1) and (-2,-2,1).
- **Saved geometry.**
  - The saved normal background is `(-σ w_di/2, s)` and the measure background is 1. Both increments are linear in epsilon, and the code checks that per component.
  - Saved V is `W0·eW_t/2 + s·ζ_t`, which equals `n_saved·r` for the plan's rate `r`.
  - c2 (`S11c_c2…:957-963`) uses `(V − n_par·u_t)/n4` and `area = sqrt(Σn²)`, exactly as the worker's `hSaved` and `rawarea` do.
- **Exact versus saved geometry.** The exact source normal is `-face·∇h/D, face/D` with measure D (`S11c_a…:844-851`). This confirms the saved values are background-truncated. The mixed debt `h_mixed − h_saved = s(V_exact − V_saved)` is algebraically right and is zero at the four retained coefficients.
- **Pullback.** The coefficient map `F_U = A·Σ(t_par + s·t4·g/2)`, `F_e = A·Σ s·W0·t4/2`, `F_c = A·Σt4` follows from `t·r`. The source-assignment candidate `(A·e·g/W0, A·e, 0)` follows from the face fourth components. The U difference first appears at σ², so the retained comparison against literal body U=0 is order debt only. E_W differs by `(A−1)e`, which is σ³.
- **Support and centre.** Free support symbols are mapped through the same formal pullback and never set. The face-centre residual and the body-face compatibility are emitted, not zeroed.
- **Wave ledger.** The order-1 and order-2 coefficients are correct. The `A1·t1` term is stated separately from the force-density slot.
- **Controls.** All three have nonzero polynomial movement against the real operands. A zero movement would refuse.
- **Mechanics.**
  - Helper census, namespace injection (`sp`, `Str`), `zero`/`stage` persisting inputs and returns before `require`, and the exact-invocation and gate joins are consistent.
  - The symbol-name uniqueness check should pass: assumption sets are uniform across the views I grepped.
  - The test stderr shows 8 tests OK. They cover label transport, hash receipts and the literal export-line join only.

**Non-blocking statements to keep explicit**
1. The "exact" normal is `n_saved/|n_saved|`. The worker never reads the geometry source. I confirmed the match at `S11c_a…:849-851`, but the code only pins the file hash. Treat it as a candidate, not a source-joined certificate.
2. Pressure/velocity-zero and chemical-zero provenance (`physical_trace_fields` versus `traction_raw/virtual_work_cases`) is documentary. The worker emits no join for it, so the plan's wording "pins those origins" means hash pins only. No acoustic-zero certificate is claimed or produced.
3. The saved operator has a nonzero σ¹ THETA body force with `kappa_theta_W`. It is not paired, and `f_hold_theta_0` is not compared. If it is recorded in the output at all, record it as UNSUPPLIED. It is not in `support-and-centre-unresolved` (line 164) now, and adding it is a one-key change.
4. The three controls use A=1. They show coefficient sensitivity only, not physical power.
5. The pullback-identity zeros are near-tautological given zero tangential face slots. They show consistent transcription, not new physics.

**Shared premises**
- The g/2 slope convention and s²=1.
- Positive W0 and L.
- The (η, σ) ≤ 1 retained rectangle, matching the source MULTIGRADE (1,0,1).
- The face-minus-body U and E gaps are order debt only.
- Higher physical E and body coefficients stay UNSUPPLIED.

**Physical obligations that remain after a pass**
- TRACTION_MEASURE_UNSUPPLIED.
- Total R0_c UNSUPPLIED.
- ORIGINAL_FORCE_PAIRING_UNRESOLVED.
- The incremental centre drive, nonuniform work, physical 20/02, end-current and threshold remainder.
- Reflection, which is survival.

**Not inspected**
- Hashes: I did not recompute any, so I did not verify the pinned worker, launcher or method digests. I also did not verify that `plan.txt` equals the `methodPath` bytes.
- Full operator and support text: I read the operator display only through the RHO4 case and the start of RHOBR. I read the support JSON in full but not every RHOBR and material-advected duplicate.
- The traction, virtual_work and closure registries: I used them only for symbol-assumption census.
- Guard and supervisor internals, and the watcher.
- The legacy residual contents. That is intentional.