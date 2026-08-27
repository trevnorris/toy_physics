# Build brief — S11c-a T7 cross-engine comparator (rev 4 = build from verified references + reviewed patches)

Three prose-spec rounds (rev-1/2/3) each drew new per-family engineering gaps from independent legs, because
the instrument is 39 bespoke per-family schemas that cannot be fully pre-enumerated in prose. rule 15:
change the author. The physics-bearing folds are now settled (three leg rounds); the working reference
implementations already handle the physics/projection families correctly; the per-family extraction is
EMPIRICALLY verifiable (a family that fails to join shows in the accounting). So this brief **re-expresses the
verified references, applies the reviewed patches, and delegates per-family extraction to you — with mandatory
accounting so no silent 0-join or false-agreement can hide.** The two build legs (fresh Claude + Grok) verify
the working instrument, not this prose.

## Object
`research/pde_ledger_v3/scripts/S11c_a_cross_engine_comparator.py`: for every one of the 39 covered tag
families, key each case on a full axis-typed key, apply the closed name/CAS-form map, and print `operand_A`
(PY), `operand_B` (WL), and `A − B` per joined case, plus per-family accounting
`{join, py_only, wl_only, duplicate_key, parse_failed, axis_set_mismatch}`. It computes and prints; decides
nothing. ⛔ No family carries a zero/nonzero target (rule 5); no "expected agreement" prior anywhere.

## Inputs (read-only, committed)
- PY: `research/pde_ledger_v3/scripts/out/S11c_a_interface_geometry_sympy_audit.out` (47 one-line
  `PY_S11CA_<Q>: <srepr>`; ⛔ do not run the engine; ⛔ `S11c_a_exports.py` is a LEDGER, not the tag stream;
  compare symbolically — PY `srepr` ordering is nondeterministic but math-invariant).
- WL: `research/pde_ledger_v3/mathematica/out/S11c_a_interface_geometry_mathematica_audit.out` (40 tags,
  multi-line payloads).
- Spec: `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md` (§1a, §1b, §3a, §3b face sum, §3c, §5–8).

## Base to re-express (verified sound by the legs for the families they cover)
`~/.s11_build/S11c_a_reconcile_fixed.py`, `S11c_a_cov_all.py`, `S11c_a_run_projection.py` (projection
`Dwin`/integral bridge), `S11c_a_scratch_loader.py` (`load_py`/`load_wl` with multi-line WL reassembly,
`split_top`/`arrow`/`py_cases`, `canon_key`). Re-express these into the committed comparator. Loaders: PY
single-line regex; WL multi-line reassembly (⛔ not S11b `read_transcript`, which does not reassemble).

## APPLY EVERY correction in these three reviewed findings files (they cite the exact payload leaves + patches)
- `~/.s11_build/comp_dir_rev3_codex_findings.md` — axis-order table, `FAMILY_FLATTENER`, `BoundIntegral`,
  VIRTUAL_WORK sum, FACE_SHIFT density mismatch, `MU_HEAD_BY_BRANCH`, `Inactive[Equal]`, residual reuse.
- `~/.s11_build/comp_dir_rev3_grok.txt` — CONTROL_FORM/HOMOGENEITY key/leaf schemas, `normalCoordinate→w`,
  `BOTH_FACES`, the current-fold (B) verdict.
- `~/.s11_build/rev3_orchestrator_derisk_findings.md` — control/bookkeeping object-nested structures.
Where those findings disagree or leave ambiguity, the SETTLED DECISIONS below govern.

## SETTLED DECISIONS (physics-bearing; three leg rounds; these override any ambiguity)
1. **Closed name/CAS map (physics-bearing).** PARAM/FIELD/FIELD_FACE/PROFILE renames; jet decode;
   `canon_jet_name`; `pynorm` (`d_w_<f>`↔`<f>_dw`); `XFACEX` face-detect; `waveOrder→1`, `virtualOrder→1`;
   width fold `deltaWidth→W_0·e_W` (§1a definitional, WL→PY, width field only); rational
   `expand`→`cancel(together)`. **ADD to the map:** the `normalCoordinate→w` binder rename;
   `Inactive[Equal]→Equal` head fold; the `Dwin` window bridge (window AST only) + `Inactive[Integrate]→
   sp.Integral`.
2. **Integral canonicalization = capture-safe `BoundIntegral`** (retain binder + `(lo,hi)`; substitute
   capture-avoidingly to a common binder; pull a factor OUT only when it is free of the binder; combine
   integrands ONLY over identical canonical `(lo,hi)`). ⛔ Not the opaque `ITG(f|_{w→W_INT})` that drops
   bounds.
3. **Perturbation-current fold = (A) ONLY.** A literal AST-head/symbol-base rename preserving every argument,
   nesting, and arity (`currentWPerturbation(args)→delta_j_bulk_4(args)`, `currentXPerturbation{1,2,3}(args)→
   delta_j_bulk_{1,2,3}(args)`, and the jet symbols). ⛔ Ban `AppliedUndef→bare Symbol` for any current name
   (the FIELD strip that hid #3). ⛔⛔ **DROP the (B) "held-context" diagnostic entirely** — do not fold, drop,
   or annotate away WL's `x,t` on the current. PY `delta_j_bulk_4(w)` and WL `delta_j_bulk_4(x1,x2,x3,{w,t})`
   remain distinct; the residual (and the `delta_rho_4D_bulk_t` density-time term) are SURFACED raw. Whether
   WL's `x,t` are inert spectators is a from-spec adjudication done AFTER this run — NOT in the comparator.
4. **MUMAP → per-branch REGISTRY, arg-preserving** (`MU_HEAD_BY_BRANCH={LAB_HELD:mu_theta_L,
   MATERIAL_ADVECTED:mu_theta_M}`; rename WL `muThetaOperand(*args)` to the branch head, args preserved),
   consulted ONLY for a case whose key pins the branch. ⛔ Never the global `mu_theta_L/M→mu_theta` collapse.
   PY's bare `mu_theta_L/M` stay structurally distinct.
5. **CONORMAL_DERIV compared raw** (no §3c Taylor fold; a phantom). Its residual is reported for the step
   record (prior verdict-A adjudication on record).
6. **Reuse S11b's typed recursive `residual`** (+ `_convert_parsed_containers`) for tuple/matrix/association/
   relation/text operands (scalar `A−B` is undefined there). Keep banning `compare_records`, `render_*`,
   `classify_*`, `main` (classification-first + `FINAL_OPERATIONAL_STATUS`).

## KEYING — full axis-typed key; report accounting (delegated; use the legs' schemas)
Build the key over a FIXED axis order `(OBJECT, BRANCH, DENSITY, FACE, DOF, VDOF, FIELD, ORIGIN, DIRECTION)`,
each token typed (not positional-guessed): reject two values for one axis; face ∈ `{PLUS, MINUS, BOTH_FACES}`
(PY `±1`/`BOTH_FACES` ↔ WL `FACE_*`); direction ∈ `{1,2,3}` (⛔ not ±1); WL `DOF_X`→`(DOF,X)`,
`VIRTUAL_DOF_X`→`(VDOF,X)`; PY positional DOF/VDOF per an explicit per-family schema table (⛔ never inferred
from token spelling). Emit the accounting per family; an unpaired/duplicate/parse-fail/axis-mismatch case is
emitted, never silently dropped.

## PER-FAMILY EXTRACTION — DISCOVER each family's structure from its payload, then extract (accounting mandatory)
The reference covers the simple shape-derivs + projection; the rest are bespoke — inspect the actual payload,
write the extractor, and REPORT. Known structures (verify against the payload):
- Simple shape-derivs (FACE_VELOCITY, RELATIVE_FLUX, TRACTION, EVOLUTION_MASS_BALANCE, CLOSURE): PY `VALUE`
  ε¹ vs WL `SHAPE_DERIVATIVE.EXPRESSION`. Triples (FACE_NORMAL, FACE_MEASURE, CONORMAL): PY `VALUE[2]` (total).
- KINEMATIC: PY `VALUE.(OPERAND_A,OPERAND_B,RESIDUAL)` vs WL `OPERAND_A_SHAPE_DERIVATIVE.EXPRESSION` etc.
  (descend to the `.EXPRESSION` leaf, not the parent association).
- VIRTUAL_CONSTRAINT: PY `VALUE` vs WL `NORMALIZED_VIRTUAL_MASS_VARIATION.EXPRESSION` (nested leaf).
- VIRTUAL_WORK: `operand_A = coeff_eps(PY.VALUE[UPPER] + PY.VALUE[LOWER])`; `operand_B =
  WL.SHAPE_DERIVATIVE.EXPRESSION` (§3b face sum). Emit UPPER/LOWER as PY context, not two comparisons.
- Projection (SHAPE_DERIV/STATIC/DYNAMIC/RESIDUAL): PY `VALUE` ε¹ (carrying `Subs(Derivative(O_window…))` +
  `sp.Integral`) vs WL `SHAPE_DERIVATIVE.EXPRESSION` | `EXPRESSION`; apply the Dwin + BoundIntegral machinery.
- Origins (EVOLUTION, PROJECTION): partition-sum. Evolution: PY named partition map vs WL per-origin cases.
  Projection: PY 2-tuple (DYNAMIC, STATIC) vs WL `DYNAMIC_SHAPE_DERIVATIVE`/`STATIC_SHAPE_DERIVATIVE` — emit
  TWO cross-engine sums.
- Controls object-nested (REP_INVARIANCE, CONTROL_INDEPENDENCE, UNIFORM_LIMIT, HOMOGENEITY): PY is
  `object → {case → operand}` (py_cases returns 0) — use `FAMILY_FLATTENER` (object outer, then base-family
  case decode). HOMOGENEITY: object rename (`RELATIVE_FLUX_LAW↔RELATIVE_FLUX`, `EVOLUTION_BALANCE↔EVOLUTION`);
  BASE/CONTROL compare PY `(term,dim)[1]` vs WL `SOURCE_TERM_DIMENSIONS`; RESIDUAL vs WL `RESIDUAL`.
- CONTROL_FORM: own object-dependent schema; keys `(OBJECT,BRANCH,DENSITY,FACE,DOF,DIRECTION)` (PY 528 vs WL
  960); the operand is NOT `VALUE`-wrapped (value directly on the outer tuple).
- FACE_SHIFT: flatten PY `VALUE` (`[0]→PRESSURE, [1][0:4]→VELOCITY, [2]→DENSITY, [3][0:4]→CURRENT`);
  ⭐ RETAIN the DENSITY axis and emit `axis_set_mismatch: WL missing DENSITY` (WL has 80 no-density cases; PY
  160 with density) — ⛔ do not dedupe PY or broadcast WL.
- Computed bookkeeping: DIMENSIONS (integer dim vectors); BACKGROUND_DENSITY_MAP (`PY.VALUE[1]↔WL.
  SIGMA_E_ZERO.EXPRESSION`, `PY.VALUE[2]↔WL.GRADIENT_SIGMA_E_ZERO.EXPRESSION`).
- Supplied bookkeeping (BACKGROUND_STATE, ADMISSIBILITY_PREMISE, FACE_MAP_*): compare as-emitted (state
  supplied/unfalsifiable); ⛔ do not broadcast a branch-agnostic premise across the other's branch cases —
  report a structure mismatch instead of manufacturing a match.
- `_LOCAL_` tags excluded (§7); emit the per-engine local inventory so the exclusion is visible.

## Definition of done (the build legs will check these empirically)
Every one of the 39 families prints its accounting line, and each has `join>0` OR a documented
`axis_set_mismatch`/`py_only`/`wl_only` with the reason. No family silently extracts 0. The script prints
operand A, B, and A−B before any guard; asserts nothing at runtime (put synthetic-fixture assertions in a
SEPARATE test file that never runs on measured payloads); exits 0 on any disagreement.

## Builder report (≤ 30 lines)
Per-family accounting summary; any family it could not extract (with the payload-shape reason); which folds
it added; runtime. State that §§1–3 + the admissibility premise + the supplied bookkeeping are
supplied/unfalsifiable and that no residual target was given.
