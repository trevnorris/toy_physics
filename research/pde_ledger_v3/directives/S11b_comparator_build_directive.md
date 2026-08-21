# S11b cross-engine comparator — build directive (folded once after two directive legs)

## Authority and boundary
Write `research/pde_ledger_v3/scripts/S11b_cross_engine_comparator.py`. It joins the two engines' emitted
tag streams on object name and reports, per shared object, whether the payloads AGREE, DISAGREE, or are
UNDECIDED — plus UNCOMPARED (operational) rows and unpaired-tag coverage. `CLAUDE.md` binds.

⭐ **Mechanical precedent — REUSE, do not re-derive:** `research/pde_ledger_v3/scripts/S10_cross_engine_comparator.py`
(proven over four fix rounds) and its repoint ablation `S10_cross_engine_comparator_repoint_ablation.py`.
Its transcript reader, `(PY|WL)_<QUANTITY>` join, Wolfram `InputForm` parse, SymPy parse, symbol
transliteration, matrix/relational/tuple `residual`, and three-valued classification are the working shape.
Reuse them. This build ADDS the rules in §"S11b delta" and REPLACES S10's places noted there.

## The frozen T7 contract (this comparator's semantics are FIXED before it is run against the real engines)
From `directives/S11_C17_C18_spec_repair_decisions_v2.md` (T7/T5) and `CLAUDE.md` rule 5:
- Join by object name: pair `PY_S11B_<QUANTITY>` with `WL_S11B_<QUANTITY>`.
- Residual the paired payloads after INJECTIVE symbol transliteration; classify AGREE (residual structurally
  0) / DISAGREE (nonzero or structural/key mismatch) / UNDECIDED (engine-emitted status token only) /
  UNCOMPARED (operational: parse or residual failure — a flagged NON-result, ⛔ never AGREE).
- ⛔⛔ **Reject a native boolean as a residual OPERAND** — the decisive T7 rule (see delta #2).
- ⛔⛔ **DISAGREEMENTS ARE THE EXPECTED, VALUABLE OUTPUT.** The comparator REPORTS cross-engine disagreement;
  it ⛔ never eliminates it. ⛔ No normalization, tolerance, key-invention, or boolean/undecided/uncompared
  bucket may turn a genuine payload difference into AGREE. A nonzero residual is a finding for the
  orchestrator, ⛔ not a build failure — there is nothing to "fix" toward.

## S11b delta over S10 (each item REPLACES or ADDS to the precedent)
1. **UNDECIDED vs UNCOMPARED — do not conflate (REPLACES S10's failure bucket).** UNDECIDED is reserved for
   an object whose payload IS an engine-emitted authoritative status token (a `STATUS_TOKEN`/coverage token
   the engine itself printed). A parse failure, an unresolved radical the residual cannot reduce, a
   structural mismatch the residual cannot process — these are **UNCOMPARED** (operational, flagged,
   reported), ⛔ never UNDECIDED and ⛔ never AGREE. ⚠ A broken parser must be VISIBLE as UNCOMPARED, never
   able to hide a difference behind UNDECIDED.
2. **Native-boolean rejection at the LEAF, not the object (REPLACES S10's `py==wl` scoring).** ⛔ Never let
   `py == wl` decide any leaf: `False == 0` and `True == 1` are Python-true and would score a
   boolean-vs-number pair as AGREE. A LEAF that is a native boolean (Python `bool`, SymPy `S.true`/`S.false`,
   or a Wolfram `True`/`False` token from either parser) is reported as its own `BOOLEAN_NOT_RESIDUALABLE`
   leaf outcome and ⛔ excluded from residual scoring — but ⛔ ONLY that leaf. ⛔ Do NOT skip the whole
   record: every OTHER leaf/sibling in the same Association/Tuple is still residualled and can still
   DISAGREE. Two `True`s are two non-residualable booleans, ⛔ not an agreement. Ordinary integers `0`/`1`
   remain fully residualable.
3. **Structure/key matching — residual only where structures MATCH; ⛔ never invent a cross-structure map.**
   The engines do NOT always emit an object with parallel internal structure (measured shapes differ:
   positional `Tuple` on one side, keyed Association on the other). Rule: residual leaf-by-leaf only when
   both payloads have the SAME shape — tuple↔tuple of equal length, Association↔Association with EQUAL KEY
   SETS (residual per matching key). Any of {tuple↔Association, unequal tuple length, unequal Association
   key sets, extra/missing key} is a **STRUCTURE/KEY DISAGREE**, reported with the mismatch. ⛔ Do NOT invent
   a key↔position table from the real payloads to force a comparison (that would let you tune toward
   agreement). A structural mismatch is an honest DISAGREE for the orchestrator to adjudicate.
4. **Wolfram `<|…|>` Associations must be preprocessed.** `parse_mathematica` RAISES on `<|…|>`. Preprocess
   into a keyed structure (key string -> parsed value) before comparison; residual per key under rule #3.
5. **Dimension tags: compare the VECTORS, ⛔ do NOT normalize.** Both engines emit `<ENGINE>_S11B_DIM_<NAME>`
   tags carrying an `L,T,M` integer vector, in different shells (a bare matrix/column on one side, an
   Association with a vector field on the other). Extract the integer vector from each shell and residual the
   two VECTORS. If they differ, report a **DIM DISAGREE** with the difference vector. ⛔ There is NO
   dimension-of-`W_0` tag and ⛔ NO W₀ normalization in this comparator: ⛔ do not rescale, absorb, or excuse
   a dimension difference — report the difference vector and let the orchestrator adjudicate off-comparator
   whether it is a known convention difference. ⚠ A wrong dimension is the class only the second engine
   catches; dimension vectors are first-class comparisons, ⛔ never skipped.
6. **Injective transliteration (ADDS a guard S10 lacks).** The symbol transliteration must be INJECTIVE: if
   two distinct source symbols would map to the same target (e.g. `a_b` and an existing `aB`), ⛔ do not
   silently collapse them — raise/flag it (`TRANSLITERATION_COLLISION`, an UNCOMPARED-class outcome for that
   object). A collision that silently collapses would make a real difference residual to zero.
7. **`_LOCAL_` tags are removed from the join entirely (§10).** A tag whose `<QUANTITY>` begins with the
   exact `LOCAL_` infix (`PY_S11B_LOCAL_…` / `WL_S11B_LOCAL_…`) is ⛔ never joined, residualled, or placed in
   the unpaired-disagreement list; list each engine's `_LOCAL_` set separately. A non-local tag present in
   only one engine IS an unpaired-coverage row that MUST be reported.

## Acceptance — executable, value-free (rule 5); FROZEN before the real run
⛔ The comparator carries NO measurement from the real engines: no real `.out` path is an input to the build,
no real count, no observed convention, no contested value. Correctness is established ONLY by SYNTHETIC
fixtures the build authors — fabricated PY/WL tag-stream pairs with PLACEHOLDER symbols and values, shaped
like the §10 grammar (positional tuples, `<|…|>` Associations, `DIM_` vectors, `_LOCAL_` tags). Each rule
above gets a fixture that a defective comparator FAILS. At minimum, and each with an EXACT expected per-name
outcome:
- **Sign-only DISAGREE**: two scalar payloads differing only in a sign → DISAGREE (⛔ not AGREE).
- **Bare-integer rows**: an equal integer pair → AGREE; an UNEQUAL integer pair (e.g. a count) → DISAGREE.
  ⛔ Both required — a comparator that mis-handles bare integers must fail here.
- **DIM vectors**: equal vectors (in the two different shells) → AGREE; vectors differing by a nonzero
  component → DIM DISAGREE with the printed difference. ⛔ No normalization makes the unequal case AGREE.
- **Nested boolean in an Association**: a record with `TEST_OBJECT -> True` beside an ALGEBRAIC sibling leaf.
  Cases: `True` vs `True`, `True` vs `1`, `False` vs `0`, `S.true` vs Python `True` — every boolean leaf is
  `BOOLEAN_NOT_RESIDUALABLE` (⛔ never AGREE), AND when the sibling algebraic leaf differs the record still
  DISAGREEs. Prove a plain `0` vs `0` sibling stays residualable.
- **Engine UNDECIDED token** vs a resolvable payload → the token row is UNDECIDED (third outcome, ⛔ not
  AGREE, ⛔ not forced to DISAGREE); a parse/residual FAILURE on a different fixture → UNCOMPARED (⛔ not
  UNDECIDED, ⛔ not AGREE).
- **Structure/key**: tuple↔Association for the same name → STRUCTURE DISAGREE; Association↔Association with an
  extra key on one side → KEY DISAGREE; ⛔ neither AGREE.
- **Transliteration collision**: a fixture where two source symbols collide under the map →
  TRANSLITERATION_COLLISION, ⛔ never a spurious AGREE.
- **Repoint ablation** (model on `S10_cross_engine_comparator_repoint_ablation.py`): substitute a DIFFERENT
  synthetic object's PY payload under a previously-agreeing name → the row must move (residual changes) AND
  flip to DISAGREE. ⛔ A symbol rename is NOT a repoint and does not satisfy this.
- **`_LOCAL_`**: an exact `LOCAL_` tag present in both engines must appear ONLY in the per-engine local lists
  (never COMPARED/DISAGREE); a PY-only and a WL-only NON-local tag must each appear in the unpaired list; a
  near-miss name containing `LOCAL` elsewhere (not the infix) is joined normally.
Every fixture + its literal stdout goes to named absolute build-scratch paths, reported in §report. ⛔ The
build does not run the comparator against the real S11b `.out` pair and ⛔ does not tune any residual to a
real payload; the orchestrator runs the frozen comparator against the real pair after review.

## The three script clauses (verbatim, non-negotiable)
1. PRINT computed objects (residuals, classifications, unpaired lists); ⛔ do NOT state a physics conclusion.
2. PRINT the residual; ⛔ do NOT assert it zero.
3. Interpretation (which disagreement matters, whether a dimension difference is a convention) belongs to the
   STEP RECORD, ⛔ not the comparator. Every residual is REACHED BY COMPUTATION from the two parsed payloads;
   every control (a fixture) re-enters at the transcript, ⛔ never at a residual.

## Report (§13) — under 25 lines
Deliverable path; each fixture + its stdout path; which delta items were implemented and where; any tag-shape
the parser cannot yet ingest (`NOT_ESTABLISHED` + what is missing); and confirm no residual/normalization was
tuned to a real payload. ⛔ Do not state whether the real engines agree — that is the orchestrator's frozen run.
