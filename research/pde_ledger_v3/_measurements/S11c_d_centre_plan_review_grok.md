I'll start by reading the packet index and the proposed plan, then check them against the governing physics and source materials.Next I’ll read the plan and the governing physics excerpts so I can form an independent view of the centre degree of freedom.I’ll inspect the native sources and saved-operand indexes next, focusing on centre rows, face work, and closure slots.I’ll keep going through `build_operator`, the saved operand indexes, and the amendment’s slab-to-face rule.I still need the rest of `build_operator`, amendment §3.3, and the closure/centre-row bookkeeping.PLAN: NEEDS REVISION

The independent centre displacement \(\zeta_c\) remains a face degree of freedom in a §1a, b §1a and c1 §1a. The supplied kinetic energy carries \(\partial_t u\) and \(\partial_t\delta W\) only. The saved b `KINETIC` rows match that: in-plane \(u_{tt}\) and \(e_{W,tt}\), with no \(\zeta_c\) inertia. The b constraint fold solves \(\delta_v\theta\) and extracts \(U\) and \(E_W\) rows only. Neither construction removes \(\zeta_c\). c2 consumes face velocity under the module constant `REPRESENTATION = 'DELTA_W'` (`S11c_c2_selfenergy_fold_sympy_audit.py` 40–61, 212). That is a supplied identification, not a derived invariant-sector result.

The centre generalized force is the coefficient of \(\delta_v\zeta_c\) in the ZETA_C virtual-work source (`face_generalized_force_rows`, b-native 81–105). The thickness face load that c2 closes is the **normalized** `FACE_VIRTUAL_WORK` \(E_W\) row (`build_operator` 634–651, 840–847). Those two rows are different pairings. The closed c2 per-face \(E_W\) origins are the image of each face’s `delta_p` / `d_w_delta_p` slots after reference-trace substitution, `retained_shape`, and `physical_fields` (`build_face` 230–234; `build_case` 290–304). They are the right reuse objects for a necessary centre-load diagnostic on the consumed thickness-source pressure. They are not a centre equation, a uniqueness proof, or a solved five-field residual.

## Blocking findings

### 1. Closed reuse must be the per-face signed combination; a single scalar on `FACE_SUM` is the wrong pairing

**Source.** a §1a: \(V_s = s\partial_t\zeta_s\), \(\zeta_s=\zeta_c+s\,\delta W/2\), outward \(n_s=s\hat w\). Saved `PY_S11CA_FACE_VELOCITY`:

- `DELTA_W`, both faces, LAB_HELD: \(\tfrac12 W_0 e_{W,t}\varepsilon\) (same sign).
- `ZETA_C`, LAB_HELD: \(+\varepsilon\zeta_{c,t}\) on face \(+1\), \(-\varepsilon\zeta_{c,t}\) on face \(-1\).
- MATERIAL_ADVECTED `DELTA_W` again equal on both faces; MATERIAL `ZETA_C` opposite in \(\zeta_{c,t}\), with a shared \(\sigma_W u_t\cdot\nabla w_1\) piece in both direction records.

c2 `main` (1190–1229) stores

`FACE_SUM = (F_{+}-(-F_{-}))/2`, `FACE_DIFFERENCE = (F_{+}-F_{-})/2`.

Saved `FACE_DIFFERENCE` / `E_W` is the literal `Integer(0)` for LAB_HELD/`RHO4_CONSTANT`, LAB_HELD/`RHOBR_CONSTANT`, and MATERIAL_ADVECTED/`RHO4_CONSTANT`. Saved `FACE_SUM` / `E_W` is the large nonzero closed object.

b `face_generalized_force_rows` differentiates thickness work in \(\delta_v e_W\) and centre work in \(\delta_v\zeta_c\). Those virtual displacements have opposite outward-face parity.

**Plan.** Step 2 extracts coefficients per face. Steps 3–4 then require “a scalar multiplier” / “this multiplier” that “maps the actual normalized pressure/jet face load to the centre load” and that “can pass through the original closure”. Step 5 builds the diagnostic “from the computed multipliers and the **saved** c2 parity combinations”.

**Effect.** On these operands the physically indicated closed object, once per-face slot coefficients exist, is

\[
\lambda_+ F_+ + \lambda_- F_- = (\lambda_++\lambda_-)\,\mathrm{FACE\_SUM} + (\lambda_+-\lambda_-)\,\mathrm{FACE\_DIFFERENCE}.
\]

A common two-face scalar applied to the summed thickness load, or a single \(\lambda\) times `FACE_SUM`, assigns the even closed \(E_W\) image (the only nonzero saved parity block) to the \(\delta_v\zeta_c\) pairing. That is the wrong reuse on the actual saved velocities, centre row, and parity returns.

The plan must split two tests:

1. Per-face slot scalars, then the signed combination above, using the saved sum and difference with c2’s conventions.
2. A single scalar for both faces, as a separate stronger claim, with a stop if it fails.

A saved `FACE_DIFFERENCE` of `Integer(0)` is a closed evenness of the thickness \(E_W\) origins. It is not the centre relation and is not a uniqueness proof. The plan already says a vanishing load is not a dynamics proof; it still needs the combination algebra so that evenness is not rewritten as a nonzero centre load on `FACE_SUM`.

This is a mandatory correction. Production must not launch until it is in the plan.

### 2. Full-row unjoined remainder must not veto reuse of the closed slot image

**Source.** The complete saved `CENTER_FACE_GENERALIZED_ROW` is 1752 bytes and is only the four slots `delta_p_±`, `d_w_delta_p_±` (same hash on all four cases). `FACE_VIRTUAL_WORK` / `ROWS` / `E_W` / `EXPANDED` is 192685 bytes (`RHO4`) or 119211 bytes (`RHOBR`). c2 `build_case` isolates per-face origins by

`(e - e.subs({slot:0 for slot in sub})).subs(sub)`

so the saved closed \(E_W\) origins are already the closed pressure/jet-slot image of the open operator.

**Plan.** Step 3: check “both slot relations and the full pressure-dependent rows”; “never drop an unjoined term”; “If a safe relation is unavailable, stop”.

**Effect.** If “safe relation” means proportionality of the entire normalized thickness row to the centre row, the job stops at the open remainder (affinity, \(\mu_\theta\), \(\Lambda_A/\Lambda_V\), velocity content, and any other non-slot terms). That remainder never enters the saved closed face \(E_W\) origins. The stop would then refuse the one reuse the inputs support: applying a successful **slot** scalar to those origins.

Required rule:

- Persist the open reconstruction residual, including every term outside the four slots.
- The closed-reuse premise is a scalar (per face) on the pressure/jet-slot content, including any pressure-dependent remainder **in those slots**.
- Non-slot open remainder is reported; it does not by itself fail the join to `CLOSED_SLAB_OPERATOR_TERM_ORIGINS` / `E_W`.
- Failure of the slot map, or a scalar that cannot pass through closure / `physical_fields` / `retained_shape` without depending on wave fields, coordinates, profile fields, or \((\varepsilon,\eta,\sigma_W)\), remains a legitimate stop.

This is a mandatory correction. Without it the stated diagnostic cannot use the saved closed returns.

## Nonblocking observations

1. **Normalization convention.** `mechanical_work_row_normalization` multiplies thickness `U` and `E_W` face rows before they are stored as `FACE_VIRTUAL_WORK` (`build_operator` 634–637, 840–847). `CENTER_FACE_GENERALIZED_ROW` in `FACE_GENERALIZED_FORCE_ROWS` / `LOCAL_SLAB_FACE_GENERALIZED_FORCE_ROWS` is the pre-multiplier \(\delta_v\zeta_c\) coefficient (705–715). The plan correctly forbids replacing normalized thickness with the unnormalized coefficient. It should name the centre operand as pre-normalization so implementers do not apply `ACTION_TO_STORED_ROW_MULTIPLIER` a second time. Units of \(\lambda\) then map stored-row \(E_W\) \([-1,-2,1]\) onto the native centre pairing.

2. **“After the native weak restriction.”** `CLOSED_SLAB_OPERATOR_TERM_ORIGINS` is `model['faces']` after close, `retained_shape`, and `physical_fields`. `extract()` is applied only to `CLOSED_COUPLING_KERNEL_TERM_ORIGINS` (`build_case` 305–311; `main` 444–451). The named \(E_W\) origins are the right objects for this centre pairing. The prose should call them post-close / pre-extract. Do not add a new `extract()` step.

3. **Pressure versus reference trace.** `build_face` writes `delta_p_s \leftarrow REFERENCE_PRESSURE` and `d_w_delta_p_s \leftarrow NORMAL_JET`, where `REFERENCE_PRESSURE` comes from `REFERENCE_VALUE_SOLVE` on physical `PRESSURE` and the jet (230–234). Open-row coefficients of `delta_p_±` therefore multiply the reference continuation. Reusing the already-substituted closed \(E_W\) origins preserves that identification. FOLD `PRESSURE` / `REFERENCE_PRESSURE` / `NORMAL_JET` are provenance; they are not a license to rerun the kernel/trace map or to rebind physical pressure into the open slots.

4. **`MULTIGRADE`.** FOLD-map grades are the whole-record aggregate and include higher shape powers (saved tuple through \((0,4,2)\)). Closed `TERM_ORIGINS` grades are the truncated support. The plan already separates those. Use the closed-row grades for the diagnostic; do not invent per-field grades.

5. **Shared MATERIAL advection.** MATERIAL `DELTA_W` and `ZETA_C` velocity records each carry \(\tfrac12\sigma_W u_t\cdot\nabla w_1\). Adding those records doubles that term. c2 already consumed `DELTA_W` only. Keep that identification; do not construct a new \(V_s\).

6. **Amendment citation.** The slab-to-face / centre-elimination dependency is Option B §2 (lines 182–188 of `accepted-scattering-form-scope.md`). §3.3 is the open-thickness coverage rule. The plan’s physics is the §2 map question. Fix the pointer.

7. **Controls.** Reversing one pressure-slot contribution and exchanging one face’s jet slot, on imported open operands only, is a routing check. Centre-row jet coefficients in the saved 1752-byte literal have the same structural factors on \(+\) and \(-\); the jet-exchange control may be degenerate there. Report that. Do not mutate completed c2 values.

8. **Claim limits.** The plan’s distinction between an arbitrary-field necessary load, a solved five-field residual, and a centre reduction is adequate. Keep it in the report. A vanishing closed centre load still leaves missing centre kinetics, a prescribed drive, and initial/boundary data.

9. **Identical open centre hashes.** The four-case `CENTER_FACE_GENERALIZED_ROW` literals share hash `bc159295…`. LAB and MATERIAL `FACE_VIRTUAL_WORK` \(E_W\) share a hash at fixed density. Closed c2 \(E_W\) origins still differ by anchoring and density. Run all four cases. Do not treat the shared open hash as a derived invariance.

## Source coverage

Read: `INDEX.md`, `centre-compatibility-plan.md`, `accepted-scattering-form-scope.md` §1–§3.3, `b-shared-physics.md`, `a-physics-context.md`, `c1-physics-context.md`, `c2-physics-context.md`, `a-native-sources.md`, `b-native-sources.md` (`face_generalized_force_rows`, `mechanical_work_row_normalization`, `kinetic_balance_from_energy`, `constraint_fold_from_source`, `build_operator` through origins), `c2-native-sources.md` (`REPRESENTATION`, `build_face`, `build_case`, `main` parity emit), `saved-operands.json` (centre rows, kinetic, constraint, face velocities), `closure-inputs.json` (FOLD `REFERENCE_PRESSURE` / `NORMAL_JET` / `PRESSURE` / `DENSITY_BINDING` / `IDENTIFICATIONS`), `reuse-inputs.json` (normalized face work, closed per-face and parity \(E_W\)), `normalized-face-literal-terms.json` (lexical slot mentions), `navigation-manifest.json`.

Large closed scientific values were used as addresses, sizes, hashes, and the copied small literals (`Integer(0)` difference, 1752-byte centre row, velocity values). No CAS, no restoration of the 1.7 MB trees, no independent coefficient extraction.

## Source joins and stopping conditions

Joins are adequate once findings 1–2 are written in: b unnormalized centre row; b normalized thickness `FACE_VIRTUAL_WORK` \(E_W\); a `DELTA_W` vs `ZETA_C` velocities without wholesale addition; c2 FOLD maps as slot provenance (reference pressure, global-\(w\) jet, density binding, `DELTA_W` identifications); c2 per-face and half-sum/difference \(E_W\) after close, without replaying kernel, Taylor, or `extract()`. Four anchoring×density cases, live RHO4 vs frozen RHOBR, and independent \((\varepsilon,\eta,\sigma_W)\) bookkeeping are in scope.

Stopping conditions are adequate for a failed **slot** map, a scalar that cannot pass through closure, a degenerate control, an undecided simplification, and the 900 s / 2 GiB cap, provided FOLD kernel payloads stay provenance-only and closed trees are not expanded. Using the saved parity combination in finding 1 keeps the algebraic job inside that cap. A failed reuse premise is an informative stop. A vanishing closed load is an informative necessary-condition result. Neither is a centre solver.

## Minimum safe next action

Revise the plan to (i) require the per-face signed closed combination and separate the common-scalar test, and (ii) tie closed reuse to the pressure/jet-slot map while persisting non-slot remainder without a full-row veto. Re-gate that text. Do not launch `_scratch/s11c/s11c-d-centre-compatibility-20260926/production` or treat this review as clearance of scientific results or upstream debts.
