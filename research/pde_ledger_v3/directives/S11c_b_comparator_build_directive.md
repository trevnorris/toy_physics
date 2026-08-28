# Build brief — S11c-b T7 cross-engine comparator (delegated build; adapt the verified S11c-a instrument)

The comparator is per-family bespoke schemas that cannot be pre-enumerated in prose (rule 15 — the S11c-a
comparator took three prose rounds before this delegated form worked). So this brief **re-expresses the
verified S11c-a comparator as the mechanical base, states the S11c-b-specific reconciliation folds, and
delegates per-family extraction to you — with mandatory accounting so no silent 0-join or false-agreement can
hide.** The two build legs (fresh Claude + Grok) verify the working instrument, not this prose.

## Object
`research/pde_ledger_v3/scripts/S11c_b_cross_engine_comparator.py`: for every emitted `S11CB_*` tag family,
key each case on a full axis-typed key, apply the closed name/CAS-form folds below, and print `operand_A`
(PY), `operand_B` (WL), and `A − B` per joined case, plus per-family accounting
`{join, py_only, wl_only, duplicate_key, parse_failed, axis_set_mismatch}`. It computes and prints; **decides
nothing** (rule 2). ⛔ No family carries a zero/nonzero target (rule 5); ⛔ no "expected agreement" prior; ⛔
no `PASS`/`FAIL`/`VERDICT`/`FINAL_STATUS`. Exit 0 on any disagreement; nonzero only on operational failure.

## Inputs (read-only, committed)
- PY: `research/pde_ledger_v3/scripts/out/S11c_b_brane_operator_sympy_audit.out` (one-line
  `PY_S11CB_<Q>: <srepr>`; ⛔ do not run the engine; `S11c_b_exports.py` is a LEDGER, not the tag stream;
  compare symbolically — PY `srepr` ordering is math-invariant but not textually stable).
- WL: `research/pde_ledger_v3/mathematica/out/S11c_b_brane_operator_mathematica_audit.out` (multi-line
  `WL_S11CB_<Q>` payloads, `<| … |>` associations with `Inactive[…]`).
- Spec: `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md` (§1a sectors, §1c energy, §3b/§3c/§3d, §5).

## Mechanical base to re-express (verified sound)
`research/pde_ledger_v3/scripts/S11c_a_cross_engine_comparator.py` — reuse its `load_py` (single-line srepr
regex), `load_wl` (MULTI-LINE association reassembly — ⛔ not S11b's non-reassembling reader), `split_top`,
`arrow`, `wl_assoc_pairs`, `wl_field`, `make_key`/axis typing, the typed recursive `residual`
(+ container conversion for tuple/matrix/association/relation/text; scalar `A−B` undefined there), the
per-family `Accounting`, and the `BoundIntegral` integral canon (retain binder + `(lo,hi)`; capture-safe;
pull a factor out only if free of the binder; combine integrands only over identical canonical `(lo,hi)`).
⛔ Do NOT reuse any classification/verdict/`main`-status machinery.

## S11c-b axis-typed keying
Fixed axis order `(OBJECT, BRANCH, DENSITY, SECTOR, SOURCE, DIRECTION, DOF)`, each token TYPED (reject two
values for one axis, never positional-guess):
- BRANCH ∈ `{LAB_HELD, MATERIAL_ADVECTED}`; DENSITY ∈ `{RHO4_CONSTANT, RHOBR_CONSTANT}`.
- OBJECT: for the UNIFORM_LIMIT / CONTROL_FORM / CONTROL_INDEPENDENCE / REP_INVARIANCE families the leading
  token is the object under control (`SLAB_OPERATOR`/`COUPLING_KERNEL`/`ADMISSIBILITY`/`ENERGY_BASIS`) — key it.
- SOURCE ∈ `{W_BG, MU_R_BG}` (the CONTROL_FORM per-profile ablation, spec §5b); DIRECTION ∈ `{1,2,3}`.
- SECTOR: the coupling-kernel block label `{TRANSVERSE_TO_THICKNESS(≡TRANSVERSE_TO_THETA_EW_UL),
  THICKNESS_TO_TRANSVERSE(≡THETA_EW_UL_TO_TRANSVERSE)}` — reconcile the PY/WL spellings as an axis rename.
Emit the accounting per family; ⛔ an unpaired/duplicate/parse-fail/axis-mismatch case is emitted, never
silently dropped.

## S11c-b reconciliation folds (physics-bearing — get these right)
1. **Inherit the S11c-a name/CAS map** (PARAM/FIELD/PROFILE renames, jet decode `canon_jet_name`, `pynorm`
   `d_w_<f>`↔`<f>_dw`, `waveOrder/virtualOrder→1`, `Inactive[Equal]→Equal`, `Inactive[Integrate]→sp.Integral`,
   rational `expand→cancel(together)`). Add the `backgroundOrder`/`wave`/`order` bookkeeper renames if the WL
   payloads carry them.
2. **μ_θ per-branch REGISTRY, arg-preserving** (`{LAB_HELD:mu_theta_L, MATERIAL_ADVECTED:mu_theta_M}`; rename
   WL `muThetaOperand(*args)` to the branch head, args preserved), consulted ONLY where the key pins the
   branch. ⛔ Never the global `mu_theta_L/M→mu_theta` collapse.
3. **⭐ The energy basis is a non-unique QUOTIENT — compare AS-IS and FLAG; ⛔ never blanket-collapse.** The
   two engines may pick different basis representatives up to the non-unique quotient (modulo total in-plane
   divergences). ⛔ Do NOT name, assume, or fold any specific representative identity — a variable-coefficient
   IBP generates first-background-jet terms that are PHYSICS, so a pre-named collapse could mask them (spec
   §1d/D5). Compare `ENERGY_BASIS_VARIABLE`/`_COUNT`/`_NEW_INVARIANTS` operands raw; a representative
   difference is a
   documented `axis_set_mismatch`/nonzero residual SURFACED for post-run adjudication, ⛔ NOT reconciled by a
   name-map collapse in the comparator. (This is the pinned schema:
   `steps/S11c_a_interface_shape_derivatives.md` "never blanket-collapse".)
4. **⭐ The coupling-kernel ADJOINTNESS residual is a density MODULO compact-support in-plane IBP.** Before
   any zero-test on the adjointness operand, reduce modulo a total in-plane divergence (an adjoint operator
   yields a total-divergence density, ∫=0). Provide an explicit `modulo_total_divergence` canonicalizer (or
   surface the raw density AND its divergence-reduced form). ⛔ Do not read a nonzero adjoint density as
   "non-adjoint" — the real operator is genuinely non-self-adjoint via the dissipative `Λ_X` face response;
   SURFACE it, adjudicate post-run. The kernel BULK blocks (`TRANSVERSE_TO_THICKNESS` etc.) compare directly.
5. **CONTROL residuals compared raw** (UNIFORM_LIMIT, CONTROL_FORM, CONTROL_INDEPENDENCE, HOMOGENEITY,
   REP_INVARIANCE): these carry their own `RESIDUAL` operands; compare A/B/A−B, no target.
6. **Supplied/bookkeeping** (any `_LOCAL_`, ADMISSIBILITY_SUPPORT_OPERAND = symbolic `f_hold_0`/`t_hold_*_0`
   premise): compare as-emitted; ⛔ do not broadcast a branch-agnostic premise across the other's cases.

## Per-family extraction — DISCOVER each family's nested VALUE structure, then extract (accounting mandatory)
Inspect the actual payload and write the extractor. Known PY shapes (verify against the payload; WL mirrors
under `<|…|>` with `.EXPRESSION`/`.WEAK_PAIRING` leaves):
- `ENERGY_BASIS_VARIABLE/_COUNT/_NEW_INVARIANTS/_OMISSIONS`: per BRANCH → `VALUE` (the basis term list / the
  Integer count / the `(SOURCE,NAME,expr)` new-invariant list). Fold 3 governs.
- `SLAB_OPERATOR` (+ `_TERM_ORIGINS`), `MU_THETA_OPERATOR`: per `(BRANCH,DENSITY)` → `VALUE` → nested
  sub-object map (`U_BODY_BALANCE`/`THETA_*`/`E_W_*` → `{LOCAL, …}` → expression). Key the sub-object; compare
  the leaf expressions.
- `COUPLING_KERNEL` (+ `_TERM_ORIGINS`): per `(BRANCH,DENSITY)` → `VALUE` → `{SECTOR_LABELS, BULK_BLOCKS →
  {TRANSVERSE_TO_THICKNESS, THICKNESS_TO_TRANSVERSE}, ADJOINTNESS_*}`. Key the SECTOR; compare the block
  expressions (bulk directly, adjointness under fold 4).
- `ADMISSIBILITY_OPERATOR_OPERAND`/`_SUPPORT_OPERAND`/`_RESIDUAL`: per `(BRANCH,DENSITY)` → `VALUE` →
  `{BODY_FORCE → {U,(THETA),(E_W)}, per-face traction}`. Compare the body-force + traction leaves.
- `UNIFORM_LIMIT_*`: per `(OBJECT,BRANCH,DENSITY)` → `VALUE` (S11c-b operand vs S11B operand vs residual);
  the S11B operand is the engine's own uniform re-derivation — compare per object.
- `CONTROL_FORM_*`: per `(OBJECT,BRANCH,DENSITY,SOURCE,DIRECTION)` → `VALUE` (base vs ablated vs residual).
- `REP_INVARIANCE_*`/`CONTROL_INDEPENDENCE_*`: per `(OBJECT,BRANCH[,…])` → operands + residual.
- `DIMENSIONS`/`HOMOGENEITY_*`: integer `[L,T,M]` vectors; compare component-wise.
- `_LOCAL_*` excluded; emit each engine's local-tag inventory so the exclusion is visible.

## Definition of done (the build legs check these empirically)
Every emitted S11CB family prints its accounting line, each with `join>0` OR a documented
`axis_set_mismatch`/`py_only`/`wl_only` + reason. ⛔ No family silently extracts 0. Prints operand A, B, A−B
before any guard; asserts nothing on measured payloads (synthetic-fixture asserts go in a SEPARATE test file);
exits 0 on disagreement. A `RUN_ACCOUNTING` summary line reports `families`, `families_with_join`,
`families_with_unpaired`, `parse_failed`, `duplicate_key`.

## Builder report (≤30 lines)
Per-family accounting summary; any family you could not extract (with the payload-shape reason); which folds
you added; runtime. State that §§1–3 + the admissibility support premise + supplied bookkeeping are
supplied/unfalsifiable and that no residual target was given.
