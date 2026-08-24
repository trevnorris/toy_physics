# S11c-a T7 comparator — SHARED-SCHEMA RE-EMIT plan (post-compact roadmap)

## Decision (user, 2026-08-24)
Build the S11c-a cross-engine comparator via **shared-schema re-emit**, NOT via comparator-side key
alignment. Both directive legs (Codex + Grok) found the comparator-side alignment (the directive
`S11c_a_comparator_build_directive.md`, now **superseded**) **UNSAFE**: the two engines serialize case keys
incompatibly, so a safe comparator would only *partially* cross-check — roughly **half the 39 objects,
including primary physics (FACE_NORMAL, traction, the projection cluster) and every control family, would
fall to KEY DISAGREE and get only orchestrator hand-adjudication** (the weak "typed conclusion, no
computation" instrument the method exists to remove — rule 2). The re-emit gives a FULL mechanical
cross-check of every object, a trivial trustworthy comparator, and fixes the same spec-§7 gap for
S11c-b…e. Cost accepted: ~2–3 build/review cycles. ⇒ `[[feedback_one_engine_fix_is_a_spec_question]]`,
`[[feedback_decompose_before_building_gates]]` (2 failed alignment designs ⇒ wrong shape).

## State at decision
- SymPy engine CLOSED+CLEAN `9b6438fa`; blind WL engine `ddecdbc2` (both committed + independently verified).
- PY tag stream (fresh run of 9b6438fa): `~/.s11_build/S11c_a_sympy_engine.out` (39 non-local + 8 local; NOT
  committed). WL committed `.out`: `research/pde_ledger_v3/mathematica/out/S11c_a_interface_geometry_mathematica_audit.out`.
- Superseded (uncommitted) comparator artifacts kept for reference only: `S11c_a_comparator_build_directive.md`
  + its `_measurements/` twin + `_legs/S11c_a_comparator_directive_review_prompt.md`. The two leg logs
  (`~/.s11_build/S11c_a_comparator_directive_{codex,grok}.log`) contain the full divergence evidence below.
- ⚠ Another session works untracked in repo `memory/` — commit EXPLICIT paths only, never `git add -A`.
- ⚠ Running the SymPy engine rewrites committed `S11c_a_exports.py` in `Dummy dummy_index=` counters only —
  `git checkout` it after any run.

## The divergence map (measured — both legs, verified at source)
Name-join is clean: **39 non-local stems, exact 1-to-1** (`grep -oE '^(PY|WL)_S11CA_[A-Z0-9_]+' | sort -u`).
Everything below is a SERIALIZATION divergence the shared §7 schema must eliminate.

**A. Case-key encoding.**
- PY = positional **Tuple of atoms**, e.g. `Tuple(Str('LAB_HELD'), Integer(1), Str('DELTA_W'))`; faces are
  `Integer 1/-1`; ⚠ integer `1` is BOTH face-plus and direction-1 (value-ambiguous). Order is per-object
  (sympy src: FACE_NORMAL=(branch,face,dof); controls add direction; VIRTUAL_WORK adds no virtual-dof;
  FACE_MAP=(1/-1); form-control=(quantity,branch,face,dof,direction) e.g. `("FACE_NORMAL","LAB_HELD",1,"DELTA_W",1)`).
- WL = single **pipe-joined string**, e.g. `"LAB_HELD|FACE_PLUS|DOF_DELTA_W"`. Grok: **49 tokens outside the
  {branch,face,dof,density} vocabulary** across 1190 keys — `DIRECTION_1/2/3`, `FIELD_PRESSURE/BULK_DENSITY/
  NORMAL_CURRENT/CURRENT_X1…`, `VIRTUAL_DOF_DELTA_W/ZETA_C`, `ORIGIN_*`, quantity/object names
  (`FACE_NORMAL`, `SIGMA_E_ZERO`, `TIME_DERIVATIVE_OF_PROJECTED_DENSITY`, …). ⚠ FACE_MAP outer keys are
  `PLUS`/`MINUS` (NOT `FACE_PLUS`/`FACE_MINUS`); PY FACE_MAP keys are `1`/`-1`.

**B. Field-envelope + field names.**
- PY inner: `VALUE`, `MULTIGRADE`, `DIMENSION_L_T_M` (flat at case level). KINEMATIC nests `OPERAND_A/OPERAND_B/
  RESIDUAL` *inside* `VALUE`.
- WL inner: `EXPRESSION`, `MULTIGRADE_EPSILON_ETA_SIGMAW`, `DIMENSION_L_T_M`. ⚠ WL nests the record under
  `SHAPE_DERIVATIVE` and beside `EXACT_SOURCE` for several objects (FACE_NORMAL, measure, velocity, traction);
  KINEMATIC siblings are `BOUND_SOURCE_LAW`, `OPERAND_A_SHAPE_DERIVATIVE`, `OPERAND_B_SHAPE_DERIVATIVE`,
  `RESIDUAL_A_MINUS_B`.
- Safe alias where flat: `VALUE↔EXPRESSION`, `MULTIGRADE↔MULTIGRADE_EPSILON_ETA_SIGMAW`; `DIMENSION_L_T_M` same.

**C. CAS representation (inherent — a trivial comparator still bridges these syntactically).**
- Held integral: WL `Inactive[Integrate][g,{normalCoordinate,-∞,∞}]` (1120) vs PY `Integral(g,(w,-∞,∞))`. Also
  `Inactive[Equal]` (1233), `Inactive[Greater]` (8). ⚠ `parse_mathematica("Equal[0,0]")` evaluates to native
  `True` — build relationals UNEVALUATED.
- Window: BOTH 2-arg — WL `windowFunction[gplus,gminus]` (660 apps, arity histogram {2:660}); PY
  `O_window(G_+,G_-)` via `Subs`/`Derivative` with `Dummy('window_argument_plus/minus', dummy_index=…)`
  (index is per-run noise → canonicalize by NAME, keep distinct). Head rename `O_window↔windowFunction` on the
  ORIGINAL head (mechanical_lower_camel gives `oWindow`, not `windowFunction`).
- Nested derivatives: WL `Derivative[o…,{a,b}][head][args]` (22656) with list-valued variable slots; the S11b
  `WL_DERIVATIVE` regex matches only flat integer orders + identifier args → these currently UNCOMPARED.
- Integer Association keys `<|1->…, -1->…|>` in BACKGROUND_STATE (18 hits) — the S11b `_parse_association_key`
  rejects integer keys → whole tag UNCOMPARED unless the parser is extended.

## The shared §7 schema to PIN (settle exactly in the amendment + 2 legs — do NOT prejudge here)
Pin, for every emitted keyed object, ONE of each:
1. **Case-key format** carrying ALL axes as explicit `(axis-label, value)` pairs (so it promotes to an
   Association via the existing S11b tuple-of-pairs machinery, or is emitted directly as an Association):
   axes = {quantity? , branch, face, dof, density, direction, field, origin, virtual_dof} as needed per
   object; ONE face convention (e.g. `FACE_PLUS/FACE_MINUS` everywhere including FACE_MAP), ONE direction
   token, no bare ambiguous integers. Both engines emit the SAME axis labels + token spellings.
2. **Flat field-envelope**: case-record = `{EXPRESSION (the computed result), MULTIGRADE_EPSILON_ETA_SIGMAW,
   DIMENSION_L_T_M}`, result under `EXPRESSION` at case level (WL flattens its `SHAPE_DERIVATIVE` nesting;
   PY renames `VALUE→EXPRESSION`). Any engine-only diagnostic (WL `EXACT_SOURCE`) is either dropped from the
   shared record or placed in a clearly-optional slot the comparator excludes.
3. The §1–§3 physics is UNTOUCHED — this amends only §7 (tag/emit grammar).
⚠ The CAS-representation differences in (C) are NOT pinned by §7 (they are the two CAS languages); the trivial
comparator still bridges those syntactically (Inactive→SymPy heads unevaluated, Dummy-by-name distinct,
nested-derivative lossless AST, integer Assoc keys) — carry the fold list Codex/Grok already produced.

## Workflow (post-compact)
0. **Complete schema enumeration** — one systematic pass: per object, both engines' exact key axes + field
   envelope. Grounds the amendment (the map above is representative, not exhaustive).
1. **§7 amendment** (orchestrator-authored; reopens closed spec `2926c71c`): pin the shared key + envelope +
   field-name schema. **2 legs (Codex + Grok) + fold once.** (rule 7 — physics-bearing: both engines read it.)
2. **Emit-layer patches** — one per engine (Codex builds each from a directive; each directive its own 2
   legs). ⛔ EMIT/serialization ONLY. **Mechanical physics-preservation gate:** extract each object's payload
   VALUE before/after and byte-compare modulo the relabeling — the physics VALUES must be UNCHANGED; only
   keys/field-names/envelope move. A value that moves is a build failure.
3. **Re-run both** (WL ~23 min ~9–10 GB, watch orphan kernel; PY ~3 min) → new `.out`s. Commit the re-emitted
   engines + the WL committed `.out`; regenerate the PY `.out`.
4. **Re-review the emit changes** — 2 legs per engine (LIGHT: confirm VALUES byte-unchanged + new schema
   correct + all 39 tags present; not a fresh physics review).
5. **Trivial comparator** — join `PY_S11CA_<Q>` ↔ `WL_S11CA_<Q>`, promote both keyed maps to Associations,
   residual per matching key, four verdicts (frozen T7 contract), + the (C) CAS-syntactic bridges only.
   Build (Codex) + 2 legs (Agent + Grok) with synthetic fixtures (rule 5). Then freeze + run against the real
   pair → the cross-engine result.
6. **S11c-a step record** (carry: T-c′ definitional-identity residual; §5a vacuous for T-d LAB_HELD virtual
   work; WL dead code `widthGradientAnsatz` 797-800) → family card. ⭐ Pin the shared §7 schema FORWARD for
   S11c-b…e so this does not recur.

## Governing reminders
Codex `-c model_reasoning_effort=xhigh`; `--sandbox danger-full-access` for Mathematica; ONE grok/user;
background via `run_in_background`, ⛔ never `&`; logs OUTSIDE repo, prompts UNDER project; rule-2 twin +
leak-gate every directive; commit reviewed baseline before overwriting; no commit before both legs report.
