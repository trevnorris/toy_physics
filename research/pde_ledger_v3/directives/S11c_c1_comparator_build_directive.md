# Build brief — S11c-c1 T7 cross-engine comparator (delegated build; adapt the verified S11c-a/S11c-b instrument)

The comparator is per-family bespoke schemas that cannot be pre-enumerated in prose (rule 15 — the S11c-a
comparator took three prose rounds before this delegated form worked). So this brief **re-expresses the
verified S11c-a comparator (and its closest sibling, the S11c-b comparator) as the mechanical base, fixes the
axis-typed keying and the emission contract, and delegates per-family extraction to you — with mandatory
accounting so no silent 0-join or false-agreement can hide.** The two build legs (fresh Claude + Grok) verify
the working instrument, not this prose.

⭐ **AUTHORITY (read first; ⛔ do not let this brief override it).** The frozen T7 contract for this sub-step is
`research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md:580-587` (N8, inherited verbatim). It is the
governing statement of what this comparator is: *"join by object name with the axis-typed keys of the S11c-a
reconciliation schema `steps/S11c_a_interface_shape_derivatives.md:233-253`, pair residual operands, reject a
native boolean as a residual operand, **three-valued**, repoint ablation. It computes and prints, deciding
nothing (rule 2). ⛔ No representational fold is pre-registered into the comparator; any cross-engine
representative difference is a computed residual adjudicated after the run. ⚠ A branch or regime the two engines
select or key differently is a computed residual to adjudicate, not a fold to pre-register (§1b: branch
selection is physics)."* The T7 source that defines "three-valued" is
`research/pde_ledger_v3/directives/S11_C17_C18_spec_repair_decisions_v2.md:42,53-60` — undecided/uncomparable is
an explicit **coverage finding the instrument emits**, not a post-run interpretation. Where this brief and that
contract appear to differ, the contract wins; report the discrepancy.

## Object
`research/pde_ledger_v3/scripts/S11c_c1_cross_engine_comparator.py`: for every emitted **shared** `S11CC1_*`
tag family (the 50 non-`LOCAL_` families), key each case on a full axis-typed key, apply the closed name/CAS-form
spelling folds below, and print `operand_A` (PY), `operand_B` (WL), and the typed recursive `A − B` per joined
case, plus per-family accounting `{join, py_only, wl_only, duplicate_key, parse_failed, axis_set_mismatch}`. It
computes and prints; **decides nothing** (rule 2).

⭐ **Three-valued (N8) — PRESERVE it in the instrument; do NOT suppress it.** The inherited `residual()`
(`S11b_cross_engine_comparator.py:682`) already returns one of: a computed residual **expression** (decidable —
zero or nonzero); a `BooleanNotResidualable` (a native boolean is NOT a residual operand — `:108`); a
`ResidualFailure` / `UndecidedResidual` (uncomparable / undecided leaf — `:114`,`:121`); or a `Mismatch`
(structural disagreement — `:100`). ⭐ **Print whichever it returns as the `A_minus_B` operand** — an
undecided/not-residualable outcome is a **coverage finding the instrument reports**, distinct from a parse
failure. ⛔ What is banned is a **per-case VERDICT/status token** the script decides:
`PASS`/`FAIL`/`VERDICT`/`FINAL_STATUS`/`AGREE`/`DISAGREE`/`STATUS`. The AGREE/DISAGREE/UNDECIDED characterization
is OUR post-run reading of the printed three-valued residual objects — ⛔ not a token you emit, and ⛔ not a
reason to strip the residual objects themselves. Exit 0 on any disagreement; nonzero only on operational failure
(missing input, ungrammatical stream). ⛔ No family carries a zero/nonzero target (rule 5); ⛔ no "expected
agreement" prior.

## Inputs (read-only, committed; ⛔ `datalad get` both before reading — they are git-annex pointers)
- PY: `research/pde_ledger_v3/scripts/out/S11c_c1_bulk_closure_sympy_audit.out` (one-line
  `PY_S11CC1_<Q>: <srepr>`; ⛔ do not run the SymPy engine; `S11c_c1_exports.py` is a LEDGER delta, NOT the tag
  stream; compare symbolically — PY `srepr` ordering is math-invariant but not textually stable). 63 tags.
- WL: `research/pde_ledger_v3/mathematica/out/S11c_c1_bulk_closure_mathematica_audit.out` (multi-line
  `WL_S11CC1_<Q>` payloads, `<| … |>` associations with `Inactive[…]`, `OperatorSum`/`OperatorComposition`/
  `FourierMultiplier`/`ScalarMultiplier`/`BoundaryMeasure`/`Conjugate`). 51 tags.
- The c1 **PY engine construction** `research/pde_ledger_v3/scripts/S11c_c1_bulk_closure_sympy_audit.py` is the
  authority for what each PY positional key slot MEANS (its axis constants and `key = (...)` lines) — read it to
  type PY keys; ⛔ do not guess a positional slot.
- Spec: `research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md` (§1 supplied inputs; §3a DtN; §3b face
  response + dissipation; §5 controls; §7 tag grammar and the keyed-CAS-map rule).

**Census (structural, not a target).** The two streams emit an **identical 50-family** non-`LOCAL_` set; the
`LOCAL_*` families differ by design (PY 13 export/fold-bookkeeping locals, WL only `LOCAL_TAG_NAMES`) and are
**excluded from the join** — emit each engine's local-tag inventory so the exclusion is visible (the S11c-b
`emit_local_inventory` shape). Iterate the join over the shared non-`LOCAL_` families; a family present in one
stream only is an accounting row, never a silent skip.

## Mechanical base to re-express (verified sound — reuse, do not re-invent)
`research/pde_ledger_v3/scripts/S11c_a_cross_engine_comparator.py` **and** the closest sibling
`research/pde_ledger_v3/scripts/S11c_b_cross_engine_comparator.py` (study it and
`test_S11c_b_cross_engine_comparator.py`). Reuse:
- `load_py` (single-line srepr regex, duplicate-tag rejection) and `load_wl` (MULTI-LINE `<| … |>` reassembly
  through the next tag line — ⛔ not S11b's non-reassembling reader).
- `split_top`, `arrow`, `wl_assoc_pairs`, `wl_field`, `py_tuple_args`, `py_pair`, `py_top_pairs`, `py_text`,
  `py_key_tokens`, `py_field` (the delimiter-aware structural readers).
- the typed recursive `residual` (imported via S11c-b from `S11b_cross_engine_comparator`) with its container
  conversion for tuple/matrix/association/relation/text and its three-valued returns above. ⚠ **`residual` DOES
  subtract scalar `sp.Expr` operands** (`S11b_cross_engine_comparator.py:~823`, `sp.factor(sp.cancel(sp.together(
  A−B)))`); a residual is undefined ONLY for a native boolean (→ `BooleanNotResidualable`) or a structurally
  mismatched / unsupported leaf (→ `Mismatch`) — reach and subtract every scalar leaf.
- the per-family `Accounting`; the `BoundIntegral` integral canon (retain binder + `(lo,hi)`; capture-safe;
  pull a factor out only if free of the binder; combine integrands only over identical canonical `(lo,hi)`).
- the S11c-b **lazy `build_case`/`materialize`/`release_case`** memory discipline (these payloads are ~82–91 MB;
  several single tags hold multi-megabyte srepr — do NOT materialize all cases before the first compare).
- the S11c-b **`RawCAS`/`raw_control_case`** outer-signature pattern — ⚠ see the STRICT scope in fold 5.
- the S11c-b **`SymbolicDifference`** oversized-leaf fallback.
- `mechanical_lower_camel` (`S11b_cross_engine_comparator.py:147`): the S10 snake_case→lowerCamel spelling rule,
  applied **INJECTIVELY** — ⚠ the FIRST piece keeps its case (`rho_m→rhoM`, `Lambda_A_0→LambdaA0`,
  `epsilon_shape→epsilonShape`). ⛔ **Check injectivity FIRST**: if two distinct reserved names collide under
  the spelling map, that is a finding to emit, never a silent merge.

⚠ **You must DEFINE this comparator's OWN axis vocabulary + typed `make_key` + `decode_py_key`/`decode_wl_key`
(the S11c-b PATTERN of redefining them for its own axis set).** ⛔ Do NOT literally reuse S11c-a's geometry
vocabulary (`FACES={PLUS,MINUS,BOTH_FACES}`, the `1→PLUS` rewrite) or its `decode_wl_key` — c1 emits
`FACE_1`/`FACE_-1`, `DIRECTION_1`, `PARITY_DELTA_W`, `OUTPUT_PROPAGATING`, which S11c-a's decoders reject with
`AxisError`. Strip a WL `FACE_`/`DIRECTION_`/`PARITY_` prefix to its typed value **without equating the axes**.
⛔ Do NOT reuse any classification / verdict / `main`-status machinery.

## S11c-c1 axis-typed keying — TYPE each token from the ENGINE'S construction; ⛔ never positional-guess, never merge distinct axes
The PY engine's axis constants (`S11c_c1_bulk_closure_sympy_audit.py:119-124`) are the closed vocabulary:
`ANCHORINGS=("LAB_HELD","MATERIAL_ADVECTED")` (**BRANCH**), `FACES=(1,-1)` (**FACE**), `DIRECTIONS=(1,2,3)`
(**DIRECTION**), `PARITIES=("THICKNESS","CENTRE_SHIFT")` (**PARITY**), `REGIMES=("PROPAGATING","EVANESCENT",
"GRAZING")` (**REGIME**, appearing as an ordered `(REGIME_OUT, REGIME_IN)` pair). Plus **DENSITY** ∈
`{RHO4_CONSTANT, RHOBR_CONSTANT}`, **OBJECT** (the object-under-control leading token, e.g. `DTN`), **MUTATION**
(`BASE`/`SIGNFLIP_INPUT`/`MOMENTUMFREEZE_INPUT`/`MOMENTUMFREEZE_OUTPUT`), and the closed **SCENARIO**/**LIMIT**/
**PORT** axes for the labelled cases (`SCENARIO`=`REAL_OMEGA_PROPAGATING_IMPERMEABLE_LAMBDA_X0_ZERO` etc.;
`LIMIT`=`OMEGA_TAU_TO_ZERO` etc.; `PORT`=the dissipation channel `A` etc.). ⛔ Each such token gets its own
named typed axis — ⛔ do NOT leave "own axis vs OBJECT sub-token" to the builder, and ⛔ do NOT drop it.

⭐ **The PY integer slot is NOT one axis — TYPE it per family from the PY constructor:**
- `Integer ∈ {1,-1}` is **FACE** in the `key=(anchoring,face[,density])` families — DtN core
  (`:1070`,`:1185`,`:1247`), face response, coeffs, energy, noninvertibility, `DTN_FLAT_SYMBOL`. The value `-1`
  ⇒ FACE (DIRECTION is never negative).
- `Integer ∈ {1,2,3}` is **DIRECTION** in `CONTROL_FORM_*` (`key=(object,anchoring,density,direction)`, `:1628`).
- ⚠ In `CONTROL_FORM`, the KEY's trailing integer is DIRECTION **and** the VALUE's inner map is keyed by FACE
  (`Integer(1)/Integer(-1)`) — these are DIFFERENT axes; ⛔ do not flatten the inner FACE onto the key's
  DIRECTION (a duplicate-axis `AxisError`) and ⛔ do not drop it (it hides per-face leaves) — key the inner map
  on FACE under the DIRECTION-keyed case.
- ⚠ **Do NOT assume `FACE_n ≡ DIRECTION_n`.** WL labels them separately (`FACE_1` in `DTN_KERNEL`;
  `DIRECTION_1` in `CONTROL_FORM`). Type each to its own axis; a legitimate join is PY-FACE↔WL-`FACE_1` within
  the FACE axis (strip the `FACE_` prefix), never FACE↔DIRECTION.

Emit the accounting per family; ⛔ an unpaired / duplicate / parse-fail / axis-mismatch case is emitted with its
reason, never silently dropped (reuse the S11c-b `axis_mismatch_detail`).

## S11c-c1 reconciliation folds — SPELLING ONLY; ⛔ every physics-bearing reconciliation is SEALED for post-run
1. **Inherit the S11c-a + S11c-b name/CAS map** (PARAM/FIELD/PROFILE renames, jet decode `canon_jet_name`,
   `Inactive[Equal]→HeldEqual`, `Inactive[Integrate]→sp.Integral`, rational `expand→cancel(together)`, the
   `__ONE__` bookkeeper renames). ⚠ **`Inactive[Greater]` and `Inactive[FourierTransform]` are NOT implemented
   in the base — keep them UNEVALUATED (held heads).** ⛔ Do NOT mirror `Inactive[Equal]` onto `Inactive[Greater]`
   so it evaluates to a native boolean (that manufactures a `BooleanNotResidualable`, or a false 0 if subtracted);
   WL `DTN_BY_REGIME_PAIR` carries `Inactive[Greater][…,0]` predicates — preserve them.
   ⭐ **Add a c1 head/symbol spelling ONLY when it is BOTH (a)** exactly `mechanical_lower_camel(<py>)` for a
   `<py>` that occurs as a real `Symbol('<py>')` in the PY stream, **AND (b)** a BARE symbol on both sides (not
   an applied head carrying arguments). The bare-symbol inverses verified present as real PY `Symbol`s are:
   `rhoM`←`rho_m`, `W0`←`W_0`, `LW`←`L_W`, `sigmaW`←`sigma_W`, `etaBg`←`eta_bg`, `epsilonShape`←`epsilon_shape`,
   `LambdaA0`←`Lambda_A_0`, `LambdaV0`←`Lambda_V_0`, `LambdaX0`←`Lambda_X_0`, `tauA`←`tau_A`, `tauV`←`tau_V`
   (re-verify each against the PY stream). ⛔ **Do NOT add** `qOut`, `kOne`/`kTwo`/`kThree`, `kPrimeOne`/…,
   `w1Profile…`/`w1Jet…`, `BoundaryMeasure`, `momentumSwapOutput`, `faceVelocityInput` — they are either applied
   heads carrying momentum/position arguments (adding them is the forbidden `AppliedUndef→Symbol` arg-strip) or
   have no real PY reserved name (a physics identification, not a spelling). ⚠ **`qOut[…]`, `w1Profile…[…]`,
   `w1Jet…[…]` STAY applied heads**, even though `mechanical_lower_camel("q_out")=="qOut"` — the head spelling
   matching is exactly the trap.
   ⚠ **The bulk sound speed is NOT a mechanical inverse.** PY has bare `Symbol('c_s0')`; WL writes `cS0`. But
   `mechanical_lower_camel("c_s0")=="cs0"` (lowercase), NOT `cS0` — so `c_s0↔cS0` is a **reviewed** spelling
   normalization, not a mechanical one. Confirm it is already in the inherited S11c-a `PARAM` map; ⛔ if it is
   not, do NOT invent a non-mechanical rename — let the `c_s0`/`cS0` difference surface as a residual and flag
   it for adjudication.
2. **μ_θ is OPAQUE and SEALED — no registry surgery.** ⚠ The S11c-b `mu_theta_L/M` / `muThetaOperand` registry
   does **not** exist in these streams (grep: absent both sides). PY emits face-specific composites
   `s11cc1_mu_theta_lab_held_{plus,minus}`, `s11cc1_mu_theta_material_advected_{plus,minus}`, `mu_theta_drive`;
   WL emits the single opaque bare symbol `muTheta`. ⛔ Do NOT build a registry mapping the PY composites onto
   `muTheta` — whether the PY composite and the WL opaque coincide is a **computed residual adjudicated after
   the run** (the μ representative fold is reviewed REGISTRY data, ⛔ never comparator operand surgery —
   `steps/S11c_a_interface_shape_derivatives.md:239-243,250-253`). Emit both as parsed (`muTheta` stays a bare
   opaque symbol; the PY composites stay as-emitted).
3. **⛔⛔ SEALED — do NOT reconcile these; surface the raw operands and let the typed residual stand for post-run
   adjudication (rule 5 / rule 6 / contract `:585-587`):**
   - **The two-momentum legs.** PY freezes them into opaque symbols `s11cc1_q_out_output`/`s11cc1_q_out_input`
     (and `s11cc1_k_output_1..3` / `s11cc1_k_input_1..3`); WL keeps them live: `qOut[omega,{kOne,kTwo,kThree}]` /
     `qOut[omega,{kPrimeOne,kPrimeTwo,kPrimeThree}]`. ⛔ Do NOT map one representation onto the other; ⛔ do NOT
     strip `qOut[...]`'s arguments to a bare symbol; ⛔ do NOT map `kOne`/`kPrimeOne` onto the PY `k_output`/
     `k_input` legs. Emit both; the residual is the measurement.
   - **The ω real-assumption.** PY carries `Symbol('omega', real=True)` in most objects but plain
     `Symbol('omega')` / `Symbol('q_out')` in some limit operands; WL carries plain `omega`. ⛔ Do NOT normalize
     `Symbol('omega')` and `Symbol('omega', real=True)` to a common assumption anywhere — emit the operands
     exactly as parsed; the assumption difference is a computed residual adjudicated post-run.
   - **The background density field (rule 17).** PY carries a BARE `Symbol('rho_br_bg_rho4_constant')`; WL
     carries the APPLIED field `rhoBrBgRho4Constant[xOne,xTwo,xThree]` (a function of position). ⛔ Do NOT
     arg-strip the WL field to the PY bare symbol — a varying background density frozen to a constant is exactly
     the rule-17 defect; surface both operands as parsed and let the residual stand for post-run adjudication
     (is PY's bare symbol a legitimate constant-density representation, or a freeze the WL field exposes?).
   - **The DtN operator ALGEBRA.** PY composes noncommutative `Symbol(..., commutative=False)` factors; WL uses
     `OperatorSum`/`OperatorComposition`/`FourierMultiplier`/`ScalarMultiplier`. Different operator
     representations — surface raw (see fold 5); ⛔ do NOT invent an operator-algebra bridge.
   - **Any parity / regime / face key the two engines select or shape differently** (§1b) — including
     `PROPAGATING` vs WL `OUTPUT_PROPAGATING`/`INPUT_PROPAGATING`, and PY `THICKNESS`/`CENTRE_SHIFT` (PARITY)
     vs WL `PARITY_DELTA_W`. ⛔ Type the literal token; ⛔ do NOT declare `OUTPUT_PROPAGATING` a "spelling" of
     `PROPAGATING`; surface the shape difference as accounting/residual, adjudicate post-run.
4. **CONTROL / limit / energy residual families compared as emitted, reaching every scalar leaf** (`CONTROL_
   BRANCH_*`, `CONTROL_FORM_*`, `CONTROL_INDEPENDENCE_*`, `UNIFORM_LIMIT_*`, `ZERO_JET_*`, `REP_INVARIANCE_*`,
   `HOMOGENEITY_*`, `ENERGY_*`): extract their `BASE`/`ABLATED`/`CORRUPTED`/operand + `RESIDUAL` leaves and
   compare A / B / A−B, **no target**. ⛔ NEVER blanket-collapse `X(args)→X` to make a control agree — hand-code
   the extractor and FLAG.
5. **⛔ `raw_control_case` (outer-signature) is permitted for EXACTLY two objects: `DTN_OPERATOR` and
   `NONINVERTIBILITY_CONDITION`'s `OPERATOR` leaf** (genuinely different operator algebras, fold 3). ⛔ For every
   other family you must **descend until no mechanically-addressable paired scalar leaf remains** — a control's
   inner scalar kernel (e.g. `CONTROL_FORM_BASE_OPERAND`'s `Integer(1)→Mul(I,omega,rho_m,…)` leaf; WL `OBJECT`
   leaf) MUST be reached and residualled, never stopped at an outer `ASSOCIATION_KEYS`/`TUPLE_ARITY` signature.
   An outer-signature comparison where a scalar leaf was reachable is a false-agreement path.
6. **Supplied / bookkeeping** (`_LOCAL_*`, and any supplied §1 premise): excluded from the join / compared
   as-emitted; ⛔ do not broadcast a branch-agnostic premise across the other engine's cases.

## Per-family extraction — ⛔ DISCOVER each family's ACTUAL PY and WL container from THIS payload (do NOT copy a sibling field name unverified)
⚠ **The two engines' VALUE containers frequently DIFFER, and several are NOT the S11c-b shapes.** For every
family, sample the real PY srepr head AND the real WL association head from these `.out`, then write the
extractor to reach the paired scalar leaves. ⛔ Do not reuse a sibling `wl_field(...)`/`decode_*` schema without
confirming the field/arity exists in THIS stream. Known traps and shapes (verify each against the payload):
- **`DIMENSIONS`** — WL field is **`OBJECT_DIMENSIONS`** (+ `DERIVATION_EQUATIONS`), NOT `NAMED_OBJECTS` (which
  the S11c-b `extract_dimensions` reads and which is **absent** here → a silent 0-extract if copied). PY is
  `(name, ImmutableDenseMatrix([L,T,M]))` tuples. Compare the `[L,T,M]` vectors component-wise AND the
  derivation-equation coverage.
- **`DTN_TERM_ORIGINS`** — ⛔ NOT the S11c-b `extract_term_origins` `(BRANCH,DENSITY)` schema. PY opens
  `(FLAT_OUTGOING_SOLUTION, …)`; WL is `LAB_HELD → ORIGINS → {BOUNDARY_SHIFT, CONORMAL_TILT, …}`. Discover it.
- **`CONTROL_FORM_*`** — ⛔ NOT the S11c-b `extract_control` 5-token `(OBJECT,BRANCH,DENSITY,SOURCE,DIRECTION)`
  schema (c1 has NO `SOURCE`; the key is 4-token `(OBJECT,BRANCH,DENSITY,DIRECTION)` → an `AxisError` on every
  case if copied). Type it 4-token; the VALUE inner map is FACE-keyed (axis note above).
- **`FACE_RESPONSE`** — PY is `(OPAQUE_MU_THETA_OPERATOR, CASES)` with cases → `RESOLVENT`/`RESOLVENT_DEFINITION`
  (noncommutative operator), NOT a `PRESSURE` map; WL keys `LAB_HELD|FACE_1|RHO4_CONSTANT → {PRESSURE, …}` AND
  carries separate `PARITY|BRANCH|DENSITY` view cases (don't omit them). `FACE_RESPONSE_COEFFS` — PY leaf is
  `PRESSURE → (tuple of Muls)`; WL is `PRESSURE_FLAT → {FACE_VELOCITY_COEFFICIENT, MU_THETA_COEFFICIENT, …}`.
- **`ENERGY_*`** — PY VALUE is `BASELINE`/`TRACTION_SIGN_REVERSED` (+ a `SCENARIO` token
  `REAL_OMEGA_PROPAGATING_IMPERMEABLE_LAMBDA_X0_ZERO`); WL is `OPERAND_A`/`OPERAND_B`/`FACE_TRACTION_PAIRING`
  with `Inactive[Integrate]` leaves (bound-integral canon). Their shells differ — reach the integrand leaves.
- **`PERMEABLE_PORT_HERMITIAN` / `PERMEABLE_DISSIPATION_VS_OMEGA_TAU`** — PY VALUE carries `FULL_BARE_DTN`,
  `BLOCK_HERMITIAN_FORM`, `STATUS_TOKEN` (no `PORT_MATRIX`); WL carries `PORT_MATRIX`/`HERMITIAN_FORM` under the
  `…|PARITY_DELTA_W` key. Match the Hermitian-form matrix leaves; the parity/regime key SHAPE differs (seal 3).
- **`DTN_KERNEL`** — PY VALUE opens `{FLAT_DIAGONAL, FIRST_SHAPE, ASSEMBLED, …}` (with `DiracDelta`); WL is
  `{OUTPUT_LEG, INPUT_LEG, MOMENTUM_TRANSFER, KERNEL.EXPRESSION}`. Compare the kernel-expression leaves (two-
  momentum seal 3 keeps `s11cc1_q_out_*` vs `qOut[...]` unreconciled).
- **`DTN_FLAT_SYMBOL`** — PY is a scalar under `(BRANCH,FACE)`; WL carries `C1_INDEPENDENT_DERIVATION`/
  `RADIATION_BRANCH` (dispersion, root candidates, selected root). ⛔ Do NOT force a join onto one representative
  — surface the shape difference; compare the scalar where both expose one.
- **`DEGENERATE_LOCI_*`** — PY is a bare `Tuple(Equality(...))`/`(var, solution)`; WL carries `VARIABLES`/
  `EQUATIONS`/`DOMAIN`/`SOLUTION` (WL-only keys). Reach the equation/solution leaves under the integral/rational
  canon.
- **`HOMOGENEITY_*`** — PY is `HEIGHT_NORMAL_NORMAL`/`TILT_DIRECTION_*` `[L,T,M]` vectors; WL is
  `SOURCE_DIMENSIONS`/`TERM_DIMENSIONS`. Compare the integer vectors component-wise where the containers align;
  surface the container mismatch otherwise.
- **DtN core** (`DTN_OPERATOR` raw per fold 5; `DTN_BY_PARITY` matrix; `DTN_BY_REGIME_PAIR` per regime-pair with
  the `OUTPUT_CONDITION`/`INPUT_CONDITION` held predicates; `DTN_HERMITIAN_PART`/`DTN_REACTIVE_PART`/
  `DTN_INERTIAL_LOADING`/`DTN_GRAZING_BEHAVIOUR`/`DTN_RIGID_SHIFT_OPERAND`/`_RESIDUAL`) — discover per payload.
- **Controls** (`CONTROL_BRANCH_*` per `(MUTATION,OBJECT,BRANCH)`; `CONTROL_INDEPENDENCE_*`; `REP_INVARIANCE_*`
  `EULERIAN`/`HANZAWA` operands + residual; `UNIFORM_LIMIT_*`/`ZERO_JET_*` `S11B`-operand vs `S11CC1`-operand vs
  residual) — extract the operand/residual leaves (fold 4).
- `_LOCAL_*` excluded; emit each engine's local-tag inventory.

## Tests (SEPARATE file `research/pde_ledger_v3/scripts/test_S11c_c1_cross_engine_comparator.py`; synthetic only)
Model on `test_S11c_b_cross_engine_comparator.py`; ⛔ do not load or run either measured engine. Include:
- typed-key construction rejects a duplicate axis and an unknown/untyped axis; **FACE and DIRECTION stay
  distinct** (a FACE case and a DIRECTION case with the same integer do NOT join); **`Integer(-1)` is rejected
  as a DIRECTION** value (FACE only);
- WL multi-line association reassembly;
- **injectivity of the spelling map** (two distinct reserved names do not collapse to one lowerCamel spelling);
- **μ_θ is never globally collapsed / registry-folded** (the PY composites and WL `muTheta` remain distinct
  operands);
- **three-valued residual is preserved**: a native boolean → `BooleanNotResidualable` (rejected as an operand,
  not silently subtracted); an undecided/uncomparable leaf → an emitted `UndecidedResidual`/`ResidualFailure`
  coverage outcome that is **distinct from a parse failure**; a **nested boolean beside an algebraic sibling** →
  only the boolean leaf is rejected while the sibling is still residualled;
- **the repoint ablation** (contract `:583-584`; model on `S10_cross_engine_comparator_repoint_ablation.py:38-47`
  and `S11b_comparator_build_directive.md:91-93`): substitute a **DIFFERENT** synthetic object's payload under a
  previously-paired NAME → the residual must MOVE (and flip away from zero). ⛔ A symbol/spelling rename is NOT a
  repoint and does not satisfy this;
- **`raw_control_case` is refused when a paired scalar leaf is reachable** (a control record with a scalar
  kernel is residualled at the leaf, not stopped at an outer signature);
- disagreement is a measurement: a synthetic PY≠WL pair prints `operand_A`/`operand_B`/`A_minus_B` +
  `ACCOUNTING` + `RUN_ACCOUNTING` and exits 0.

## Definition of done (the build legs check these empirically)
Every emitted shared `S11CC1` family prints its `ACCOUNTING` line. ⛔ **No family silently extracts 0**; every
family either joins cases OR emits documented `axis_set_mismatch`/`py_only`/`wl_only` rows with a reason. ⛔ **Do
NOT use `join>0` as an exit criterion** (it pressures manufacturing a join, rule 5) — and one `raw_control_case`
or generic case per family is NOT coverage: report a per-family **extracted-leaf count** so an under-measured
family is visible. Prints `operand_A`, `operand_B`, `A_minus_B` before any guard; asserts nothing on measured
payloads (synthetic-fixture asserts live in the test file only); exits 0 on disagreement, nonzero only on
operational failure. A `RUN_ACCOUNTING` summary reports `families`, `families_with_join`,
`families_with_unpaired`, `parse_failed`, `duplicate_key`, `runtime_seconds`; a `MEASUREMENT_SCOPE` line records
that §§1–2 + the supplied substrate/`μ_θ` operand + supplied bookkeeping are supplied/unfalsifiable and
`residual_target=none`. `LOCAL_INVENTORY` lines make the `_LOCAL_` exclusion visible for both engines.

⚠ **Memory / runtime.** These payloads are large; per S11c-b the FULL cross-engine residual may exceed this box
(30 GB). Keep the lazy materialize/release discipline; if the full run OOMs, that is an operational deferral
(recorded out-of-band), NOT a reason to narrow the comparison — report it, do not silently drop families.

## Builder report (≤30 lines)
Per-family accounting + extracted-leaf-count summary; any family you could not extract (with the payload-shape
reason); which spelling folds you added and the `mechanical_lower_camel` inverse + real-PY-`Symbol` check that
justifies each; the injectivity check result; the `raw_control_case` whitelist you used; runtime and peak RSS.
State that §§1–2 + the supplied substrate + the `μ_θ` operand + supplied bookkeeping are supplied/unfalsifiable
and that **no residual target was given** (rule 5).
