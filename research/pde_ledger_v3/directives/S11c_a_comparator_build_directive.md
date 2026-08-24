# S11c-a cross-engine comparator — build directive

## Authority and boundary
Write `research/pde_ledger_v3/scripts/S11c_a_cross_engine_comparator.py`. It joins the two S11c-a engines'
emitted tag streams on object name and reports, per shared object, whether the payloads AGREE, DISAGREE, or
are UNDECIDED — plus UNCOMPARED (operational) rows and unpaired-tag coverage. `CLAUDE.md` binds. It reads two
`.out` transcripts given as argv; that read and its printed report are its only effects.

⭐ **Mechanical precedent — REUSE, do not re-derive:** `research/pde_ledger_v3/scripts/S11b_cross_engine_comparator.py`
(the frozen T7 contract, proven over a build + two fix rounds). Its transcript reader, `(PY|WL)_<QUANTITY>`
join, SymPy `srepr`/`sympify` parse, Wolfram `parse_mathematica` parse with `<|…|>`-Association handling,
`Derivative`→`partialDerivative` and `Inequality` preprocessing, INJECTIVE `mechanical_lower_camel`
transliteration with the `(kind,target)` collision guard, recursive structural `residual`
(`factor(cancel(together(py−wl)))` under a per-leaf time budget), native-boolean rejection, the four-verdict
`classify_residual` with precedence **DISAGREE > UNCOMPARED > UNDECIDED > AGREE**, `_LOCAL_` exclusion, and
the exit-code policy — are the working shape. **REUSE them unchanged** except where §"S11c-a delta" RETARGETS
or ADDS. ⛔ Do not re-implement or weaken any of them.

## The frozen T7 contract (semantics FIXED before this runs against the real engines)
Carried over verbatim from the S11b comparator (`CLAUDE.md` rule 5):
- Join by object name: pair `PY_S11CA_<QUANTITY>` with `WL_S11CA_<QUANTITY>` by exact non-local object name.
- Residual the paired payloads after INJECTIVE transliteration and the syntactic canonicalizations below;
  classify **AGREE** (every residual leaf structurally `0`) / **DISAGREE** (any nonzero residual leaf, or a
  structure/key mismatch) / **UNDECIDED** (a lone engine-emitted authoritative `STATUS_TOKEN`/coverage token,
  no severer sibling — the S11b `_is_authoritative_uncertain_leaf` gate) / **UNCOMPARED** (operational: parse
  failure, residual failure, budget-exceeded, transliteration collision, native-boolean leaf — a flagged
  NON-result, ⛔ never AGREE). Precedence DISAGREE > UNCOMPARED > UNDECIDED > AGREE.
- ⛔⛔ **Reject a native boolean as a residual OPERAND** (Python `bool`, `sp.true/false`, Wolfram `True/False`)
  at the LEAF, excluded from residual, ⛔ never AGREE; every sibling leaf still residuals.
- ⛔⛔ **DISAGREEMENTS ARE THE EXPECTED, VALUABLE OUTPUT.** The comparator REPORTS cross-engine disagreement;
  it ⛔ never eliminates it. ⛔ No canonicalization, alias, tolerance, or key-normalization below may turn a
  genuine payload-VALUE difference into AGREE. A nonzero residual, a structure/key mismatch, or an
  unreconciled representation is a finding for the orchestrator to adjudicate, ⛔ not a build failure — there
  is nothing to "fix" toward.
- No terminal verdict (`G10`): the tool prints per-object classifications, per-category lists, and a
  `FINAL_OPERATIONAL_STATUS: PASS|FAIL` that is **operational only**. Exit `0` = ran cleanly (any mix of
  AGREE/DISAGREE/UNDECIDED); `1` = an operational failure is present (shared-name UNCOMPARED, duplicate rows,
  format issue); `2` = malformed inputs (`ComparatorInputError`).

## S11c-a delta over the S11b comparator (each RETARGETS or ADDS; ⛔ preserves the contract)

**⛔⛔ GOVERNING RULE for every item below — bridges align REPRESENTATION, never VALUES.** Each
canonicalization/alias must be (i) a fixed table or purely syntactic rule, ⛔ never derived from the real
payloads; (ii) injective / value-preserving — it may rename a symbol, a function head, a key, a bound
variable, or a held-integral head, but it ⛔ may never move, drop, reorder, or coerce a payload VALUE, and
⛔ may never map two structurally different objects onto each other. Any case, key, field, or head that does
NOT match after the SAFE rules below is reported as a DISAGREE/UNCOMPARED with the specific mismatch, and the
matching siblings are STILL residualled — ⛔ a mismatch on one field never silently suppresses the physics
comparison of the others.

1. **Retarget the tag namespace `S11B_` → `S11CA_`.** The three hardcoded sites in the precedent:
   `TAG_LINE` regex, `is_local_name`, `_is_dimension_name`. (⚠ There are NO `S11CA_DIM_*` tags — the
   dimension special-path is inert here; L,T,M vectors are ordinary nested `DIMENSION_L_T_M` leaves and
   residual as integer tuples. Retarget the prefix anyway for correctness.)

2. **Outer case-key alignment (the keyed maps) — by dimension-classified token, order-independent.** The two
   engines encode a case key differently: **WL** = a single pipe-joined string, e.g.
   `"LAB_HELD|FACE_PLUS|DOF_DELTA_W"`; **PY** = a `Tuple` of atoms, e.g.
   `Tuple(Str('LAB_HELD'), Integer(1), Str('DELTA_W'))` (faces are integers `1`/`-1`). Canonicalize BOTH to a
   set of `(dimension, value)` pairs by classifying each token by its VALUE against these disjoint,
   fixed vocabularies, then join cases whose canonical `(dimension→value)` maps are EQUAL:
   - branch ∈ {`LAB_HELD`, `MATERIAL_ADVECTED`} (identical both engines)
   - face: `FACE_PLUS` ⇔ `1`, `FACE_MINUS` ⇔ `-1`
   - dof ∈ {`DELTA_W`, `ZETA_C`} (WL prefixes `DOF_`; strip it)
   - density ∈ {`RHO4_CONSTANT`, `RHOBR_CONSTANT`} (identical)
   Because the four vocabularies are disjoint, classification is unambiguous and order-independent. A case
   present in one engine with NO equal-map counterpart in the other is a **KEY DISAGREE** (surfaced with the
   unmatched canonical key), ⛔ never dropped and ⛔ never force-paired positionally. ⛔ Do not assume both
   engines key an object over the same dimensions — if the canonical key SETS differ, that is the finding.

3. **Inner field-name alias — a FIXED two-entry table on Association KEYS only.** Apply on both sides before
   the key-set comparison: PY `VALUE` → `EXPRESSION`; PY `MULTIGRADE` → `MULTIGRADE_EPSILON_ETA_SIGMAW`.
   (`DIMENSION_L_T_M` is already identical.) ⛔ This table is applied to KEYS only, ⛔ never to values, and is
   the WHOLE table — ⛔ do not add, infer, or object-condition any other alias. Fields that remain unmatched
   after this table — WL-only `EXACT_SOURCE` (the exact pre-jet source form, on a few objects), PY-only
   `OPERAND_A`/`OPERAND_B`/`RESIDUAL` (on `KINEMATIC_BALANCE`), or any other — are reported as the record's
   KEY difference (extra-on-WL / extra-on-PY, listed), while every MATCHING field (`EXPRESSION`,
   `MULTIGRADE_EPSILON_ETA_SIGMAW`, `DIMENSION_L_T_M`) is STILL residualled and its AGREE/DISAGREE reported.

4. **Held projection integrals — canonicalize head and bound variable, ⛔ never the integrand.** WL emits
   `Inactive[Integrate][integrand, {normalCoordinate, -Infinity, Infinity}]`; PY emits
   `Integral(integrand, Tuple(Symbol('w'), NegativeInfinity, Infinity))`. Canonicalize both to one held
   integral form over one canonical bound symbol (alpha-rename WL `normalCoordinate` and PY `w` to a single
   shared bound symbol; map `Inactive[Integrate]` → the same held head `sp.Integral` carries). Likewise
   `Inactive[Equal]` → `Eq`, `Inactive[Greater]` → `Gt`, and any other `Inactive[<rel>]` → its SymPy
   relational. ⛔ These touch only the HEAD and the BOUND-VARIABLE NAME; the integrand is residualled
   unchanged, so a real integrand difference still DISAGREEs. ⚠ `parse_mathematica` does NOT map
   `Inactive[Integrate]` itself — add it as a preprocessing substitution alongside the existing
   `Derivative`/`Inequality` ones.

5. **Nested-order derivatives.** WL carries derivatives of the form `Derivative[o1, o2, …, {a, b}][head][args]`
   (a multi-index `{a,b}` over a list-valued variable slot). The precedent's `WL_DERIVATIVE` regex matches
   only flat integer-order lists and silently leaves these to raw `parse_mathematica`. Extend the derivative
   preprocessing to canonicalize this nested-order form to the SAME `partialDerivative` head the flat form and
   the PY `Derivative` side land on, so equal derivatives residual to `0`. ⛔ A derivative of a different
   order must still DISAGREE.

6. **SymPy `Dummy` canonicalization — by name, ⛔ never by index.** PY payloads carry
   `Dummy('window_argument_plus', dummy_index=N)` / `…_minus` whose `dummy_index` is a per-run counter
   (non-deterministic between runs). Canonicalize each `Dummy(name, dummy_index=…)` to a stable form keyed by
   `name` only (e.g. a canonical `Symbol`/`Dummy` per distinct name), so the volatile index never enters the
   residual. ⛔ Do NOT collapse two DISTINCT names to one. (Without this, every PY projection payload is
   spuriously non-comparable.)

7. **Window function head — a name rename ONLY; ⛔ NEVER reconcile its parameterization.** PY names the window
   `O_window` and applies it to TWO arguments `O_window(G_+, G_-)` (the spec §3c `𝒪(G_+,G_-)`, via
   `Subs`/`Derivative` of an abstract function); WL names it `windowFunction` and (as emitted) applies it to a
   single composite argument. Treating `O_window`↔`windowFunction` as a transliterated head name is
   permitted. ⛔⛔ **But do NOT bridge the ARITY or the argument parameterization** — ⛔ never rewrite a
   two-argument window into a one-argument one (or vice-versa), and ⛔ never introduce an argument map to make
   them coincide. If, after the head rename, the two window objects differ in arity/arguments, that is a
   **structural DISAGREE** on the projection object for the ORCHESTRATOR to adjudicate (real physics
   difference vs. an unreconcilable representation choice), ⛔ never reconciled toward AGREE. ⚠ This is the
   one place a masked bridge would defeat the whole cross-engine test.

## Acceptance — executable, value-free (rule 5); FROZEN before the real run
⛔ The comparator carries NO measurement from the real S11c-a engines: no real `.out` path is a build input,
no real count/convention/contested value. Correctness is established ONLY by SYNTHETIC fixtures the build
authors — fabricated PY/WL tag-stream pairs with PLACEHOLDER symbols/values, shaped like the real grammar
(pipe-string vs tuple case keys with `FACE_PLUS`/`1`, `VALUE`/`EXPRESSION` fields, `Inactive[Integrate]` vs
`Integral`, nested `Derivative[…,{a,b}]`, `Dummy(dummy_index=…)`, `O_window`/`windowFunction`, `_LOCAL_`).
Carry over EVERY S11b fixture (sign-only DISAGREE; equal-integer AGREE and unequal-integer DISAGREE; native
boolean beside an algebraic sibling; UNDECIDED token vs a resolvable payload; parse/residual failure →
UNCOMPARED; tuple↔Association STRUCTURE DISAGREE and extra-key KEY DISAGREE; transliteration collision;
`_LOCAL_` handling; the repoint ablation) and ADD one per delta above, each with an EXACT expected per-name
outcome and each of which a DEFECTIVE comparator FAILS:
- **Case-key alignment**: a WL `"LAB_HELD|FACE_PLUS|DOF_DELTA_W"` case and its PY
  `Tuple(Str('LAB_HELD'),Integer(1),Str('DELTA_W'))` counterpart with EQUAL values → AGREE; identical keys but
  differing values → DISAGREE; a WL case whose canonical key has no PY counterpart → KEY DISAGREE (⛔ not
  dropped, ⛔ not force-paired).
- **Field alias**: `VALUE`↔`EXPRESSION` with equal values → AGREE; with DIFFERING values → DISAGREE (the
  alias must NOT mask a value difference); a record with a WL-only `EXACT_SOURCE` beside a matching
  `EXPRESSION` → the `EXPRESSION` still residuals (AGREE/DISAGREE on its own), and the record is flagged
  KEY-different for the extra field.
- **Held integral**: `Inactive[Integrate][g,{normalCoordinate,-∞,∞}]` vs `Integral(g,(w,-∞,∞))` with the same
  `g` (modulo the bound-var name) → AGREE; a DIFFERENT integrand → DISAGREE; a bound-variable-name-only
  difference → AGREE.
- **Nested derivative**: `Derivative[0,0,0,{0,1}][f][…]` vs the PY equivalent → AGREE; a different order →
  DISAGREE.
- **Dummy**: two payloads identical except `dummy_index` → AGREE; two DISTINCT dummy NAMES → not collapsed
  (they residual as distinct, so a payload that differs only by swapping the two window arguments DISAGREEs).
- **Window arity**: `O_window(a, b)` (2-arg) vs `windowFunction(c)` (1-arg) under the head rename →
  structural DISAGREE (⛔ NOT AGREE, ⛔ not reconciled); `O_window(a,b)` vs `windowFunction(a,b)` with equal
  args → AGREE (head rename is legitimate).
Every fixture + its literal stdout goes to named absolute build-scratch paths, reported in §report. ⛔ The
build does not run the comparator against the real S11c-a `.out` pair and ⛔ does not tune any bridge to a
real payload; the orchestrator runs the frozen comparator against the real pair after review.

## The three script clauses (verbatim, non-negotiable)
1. PRINT computed objects (residuals, classifications, unpaired lists); ⛔ do NOT state a physics conclusion.
2. PRINT the residual; ⛔ do NOT assert it zero.
3. Interpretation (whether a projection DISAGREE is real physics or an unreconcilable window representation,
   whether an extra `EXACT_SOURCE` field matters) belongs to the STEP RECORD, ⛔ not the comparator. Every
   residual is REACHED BY COMPUTATION from the two parsed payloads; every control (a fixture) re-enters at the
   transcript, ⛔ never at a residual.

## Report (§report) — under 25 lines
Deliverable path; each fixture + its stdout path; which delta items were implemented and where; any real
tag-shape the parser cannot yet ingest (`NOT_ESTABLISHED` + what is missing); and confirm no bridge/residual
was tuned to a real payload. ⛔ Do not state whether the real engines agree — that is the orchestrator's
frozen run after review.
