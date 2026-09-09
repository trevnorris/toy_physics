# Build brief — S11c-c2 N6 cross-engine comparator (delegated build; adapt the verified S11c-c1 instrument)

The N6 comparator is per-family bespoke schemas that cannot be pre-enumerated in prose (the S11c-a/c1 comparators
each took several rounds before this delegated form worked). So this brief **re-expresses the verified S11c-c1
comparator (and its S11c-a/S11c-b base) as the mechanical substrate, fixes the axis-typed keying and the
emission contract, states the N6-specific SEALS, and delegates per-family extraction to you — with mandatory
accounting so no silent 0-join or false-agreement can hide.** The two re-review legs (a fresh Claude/Opus agent
+ Grok) verify the working instrument, not this prose.

⭐ **AUTHORITY (read first; ⛔ do not let this brief override it).** The frozen `T7`/`N8` contract for this
sub-step is `research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md` **§N8 (the "The comparator (`N8`,
frozen `T7` contract)" paragraph, ~line 477)** + **§5c (the two-route N6 control, ~line 303)** + **§1b
(cross-engine AGREE vs UNDECIDED, ~line 76)**. It governs what this comparator is: *join the two blind engines'
emitted objects by name, pair residual operands, be **three-valued**, reject a native boolean as a residual
operand, PRINT and decide nothing (rule 2); surface the §3d representation questions; the reconcile is the
**staged representational bridge**, ⛔ never a blanket collapse; the giant families + the full per-family
symbolic residual remain **deferred** (`DEFERRED_HEAVY_RUNS.md`) — c2 **names, does not pre-adjudicate**,
whatever it cannot close on this box.* Where this brief and that contract appear to differ, the contract wins;
report the discrepancy. The definition of "three-valued" is `S11_C17_C18_spec_repair_decisions_v2.md:42,53-60`
(undecided/uncomparable is an explicit **coverage finding the instrument emits**, not a post-run interpretation).

⚠ **N6-SCOPING (this instrument is NOT the full c2 self-energy comparator).** §N8's *"surface the §3d representation
questions"* duty AND its *"giant families deferred ≥64 GB"* both concern the c2 **self-energy-increment** comparator;
the N6 streams carry NEITHER the §3d families (`t_s`, DtN whole-form, flat-symbol — verified ABSENT from the N6
emitters) NOR the ≥64 GB giants (N6's largest operand is ~80 MB). This N6 instrument implements N8's GENERIC clauses
(join-by-name, three-valued, no-native-boolean, no-blanket-collapse, name-don't-adjudicate-what-you-cannot-close) +
the **§5c two-route N6 physics** + a REQUIRED **PIT/blindness SEAL** clause (the engines' deliberately-different prime
sets — a required addition, ⛔ NOT a spec-currency error), deferring per-object only over a stated budget.

## Object
`research/pde_ledger_v3/scripts/S11c_c2_N6_cross_engine_comparator.py`: for every **cross-engine-shared** N6 tag
family (see the JOIN SET below), key each case on a full axis-typed key `(ANCHORING, DENSITY[, inner axes])`,
apply the closed name/CAS-form spelling folds, and print `operand_A` (SymPy), `operand_B` (WL), and the typed
`A − B` **structural + symbolic** residual per joined case, plus per-family accounting
`{join, sympy_only, wl_only, duplicate_key, parse_failed, axis_set_mismatch, unmatched_key, deferred_oversize,
zero_extract_failures}`. It computes and prints; **decides nothing** (rule 2). A separate test file lives at
`research/pde_ledger_v3/scripts/test_S11c_c2_N6_cross_engine_comparator.py` (synthetic fixtures only).

⭐ **Three-valued (N8) — PRESERVE it; do NOT suppress it.** Reuse the inherited `residual()` machinery
(`S11b_cross_engine_comparator.py`, imported through the S11c-c1 comparator) which returns one of: a computed
residual **expression** (decidable); a `BooleanNotResidualable` (a native boolean is NOT a residual operand); a
`ResidualFailure`/`UndecidedResidual` (uncomparable/undecided leaf); or a `Mismatch` (structural disagreement).
⭐ Print whichever it returns as the `A_minus_B` operand — an undecided/not-residualable/deferred outcome is a
**coverage finding the instrument reports**, distinct from a parse failure. ⛔ BANNED: any per-case
VERDICT/status token the script decides (`PASS`/`FAIL`/`VERDICT`/`AGREE`/`DISAGREE`/`STATUS`/`COVARIANT`). The
AGREE/UNDECIDED reading is OUR post-run adjudication of the printed residual objects — ⛔ not a token you emit.
Exit 0 on any disagreement; nonzero only on operational failure (missing input, ungrammatical stream). ⛔ No
family carries a zero/nonzero target (rule 5); ⛔ no "expected agreement" or "expected covariance" prior.

## ⛔⛔ THE N6 CROSS-ENGINE SEAL + the on-box channels (READ BEFORE DESIGNING THE JOIN)
The two engines are **deliberately blind** (the WL engine derives its own δ; `S9_export_chain_rebuild_directive.md:16-18`).
⭐ **SEAL means: surface each engine's value, ⛔ NEVER difference/equate it across engines.** SEAL is reserved for
EXACTLY these two:
1. ⛔⛔ **PIT residue tables + prime sets + sample points.** The two engines use **different prime sets** (the third
   prime differs) and **different draw counts + sample generators** (no shared sample manifest). A GF(p) residue at
   one engine's points mod its primes is meaningless subtracted from the other's. ⛔ NEVER difference
   `PROBE_NUMERATORS`/`numerator_denominator`; ⛔ never re-sample to "align" (that builds a THIRD engine, not a join).
   Each engine's PIT is its OWN internal zero-witness.
2. ⛔ **Digests.** WL integer `INPUT/OUTPUT_FINGERPRINT` vs SymPy sha256 hex (`source_sha256`,
   `builder_fingerprints`, `operand_fingerprints`) — different algorithms; surface each, ⛔ never difference/equate.

⚠⚠ **SEAL ≠ "skip the residual."** For a PHYSICS-BEARING representation difference — above all the blind WL carrier
`C_E` vs the imported-slab SymPy `C_E` — you STILL surface both operands and PRINT the three-valued symbolic residual;
you just ⛔ do not FOLD it to identity (the representational reconcile is the T7's post-run job,
[[feedback_reconcile_representational_bridge]]). ⛔ Do NOT mark a carrier/source operand "sealed" and skip its
symbolic residual — that discards the one real cross-engine channel for `C_E`.

⭐ **The TWO real cross-engine channels (what the comparator computes):**
- **(A) STRUCTURAL — compared ON MATCHED KEYS ONLY.** For each joined family compare the **dimensions** (`[L,T,M]` +
  support/consistency) and the **nonzero-SUPPORT**. ⚠⚠ The column bases are **DIFFERENT VOCABULARIES, not a
  permutation** ⇒ a typed inner-axis DECODER + pre-registered SCHEMA-normalization (see THE SCHEMA BRIDGE below) is
  MANDATORY; a post-run "figure out the columns" is a silent-0-join / false-count trap. ⭐ **Support is ONE-SIDED
  evidence, ⛔ not a two-valued boolean:** type each observation `NONZERO_WITNESSED` (a positive witness) vs
  `NO_NONZERO_FOUND` (⛔ NOT "structurally zero" — it propagates as **UNDECIDED**, both engines' samplers being
  one-sided); ⛔ never subtract native support booleans across engines.
- **(B) SYMBOLIC — the PRIMARY channel, TRACTABLE for essentially all N6 objects.** Both sides carry a symbolic form
  (WL `ARITHMETIC → Inactive[…]`; SymPy `<TAG>_NODES` reconstructed via its `ARITHMETIC_DAG`). ⚠ N6's LARGEST WL
  operand is ~80 MB (`CLOSURE_GUARD_RESIDUAL`) — the ordinary-tractable range on this 30 GB box, ⛔ NOT the ≥64 GB
  S11c-b bulk residual / c1 giants (those are a DIFFERENT deferral — ⛔ do NOT copy the 64 GB figure onto N6). So
  parse+translate+difference the symbolic operands (three-valued). A single object that exceeds a **stated per-object
  budget** (name explicit RSS + wallclock ceilings in the DoD) becomes a `deferred_oversize` accounting row (recorded
  to `DEFERRED_HEAVY_RUNS.md`), ⛔ never silently dropped, ⛔ never forced. Deferral is the RARE exception for N6, not
  the plan.

⭐ **SymPy-only families have NO WL sibling — `sympy_only` accounting rows, ⛔ NEVER a disagreement, ⛔ NEVER a residual
target.** The complete SymPy-only set = **every diagnostic `S11CC2_N6_*` EXCEPT the six guards** (`SLOT_GUARD_{NATIVE,
CARRIER,RESIDUAL}`, `CLOSURE_GUARD_{NATIVE,CARRIER,RESIDUAL}`, which DO have WL siblings) **plus every
`REP_INVARIANCE_*`, every `CONTROL_INDEPENDENCE_*`, and `PREMISES`** (WL's `N6_` namespace is thin: guards + `LOCAL_*`
+ `N6RC_DIMENSIONS` only). Emit each as a `sympy_only` accounting row whose reason is **"no WL sibling"** — ⛔ state
NOTHING about its expected value, sign, or that it is a superseded object (rule 5). They are SymPy-internal witnesses,
⛔ not part of the join.

## ⭐⭐ THE SCHEMA BRIDGE — pre-register the MECHANICAL key normalization; ⛔ NEVER pre-register a PHYSICS equality
The two engines key the same objects with DIFFERENT axis vocabularies (verified from the constructors). A blind
builder MUST have a typed inner-axis decoder built from BOTH constructors, plus a pre-registered SCHEMA-normalization
map — otherwise it silent-0-joins, compares nonzero COUNTS (false structural agreement), or invents a physics pairing.
⭐ **Pre-register ONLY the mechanical name/CAS identities** (deterministic emitter-schema relations, the SPELLING
analogue): grade `{1,η,σ} ↔ (η,σ)`; face `"SUM" ↔ 0`; axis-ORDER normalization to one common tuple order; the
`_NODES` literal decoding (see per-family extraction). Build `decode_py_key`/`decode_wl_key` from the two
constructors — SymPy `scripts/S11c_c2_N6_{reconcile,diagnostic}_sympy.py` (`columns`/`BLOCKS`/face-slot); WL
`mathematica/S11c_c2_N6_mathematica_audit.wl` (`covNames`/carrier/weak construction) — each token its OWN typed axis.
⛔⛔ **Do NOT pre-register a PHYSICS equality:** ⛔ NOT block-name equality (SymPy `(THICKNESS_TO_TRANSVERSE, THETA)`
vs WL `THETA_FROM_TRANSVERSE`), ⛔ NOT kernel-family equality (SymPy `0/6/9/12` vs WL `LOCAL_BARE`/`FOURIER_KOUT_Y`/…),
⛔ NOT WL-per-component (`COMPONENT_AXES {1,2,3}`) ↔ SymPy-formal-basis-scalar equality, ⛔ NOT WL-wave-column ↔
SymPy-wave-summed equality. Those are physics identifications the comparator must NOT assume — a wrong one
manufactures agreement ([[feedback_handcode_comparison_never_blanket_collapse]]).
⭐ **Emit a PAIRING TABLE per family** `{matched, sympy_only_column, wl_only_slot}`: compare support+symbolic ONLY on
the mechanically-matched keys; every unmatched key is a COMPUTED RESIDUAL / accounting row adjudicated post-run,
⛔ never a count, ⛔ never a silent drop. (⚠ The two decision legs SPLIT here — resolved to the conservative
name/CAS-only pre-registration + pairing table, ⛔ NOT the broader block/kernel-family pre-registration, per rule 5.)

## JOIN SET — the cross-engine-shared N6 families (both engines emit them)
Discover the exact per-family container from THIS payload (below), but the join set is:
- **Reconcile (SymPy `S11CC2_N6RC_*` ↔ WL `WL_S11CC2_N6RC_*`):** `R_N6`, `SPLIT_SUM`, `SPLIT_CHECK`,
  `CARRIER_EULERIAN`, `CARRIER_MATERIAL`, `CARRIER_CHANNEL`, `CARRIER_BRIDGE_RESIDUAL`, `SOURCE_EULERIAN`,
  `SOURCE_MATERIAL`, `SOURCE_CHANNEL`, `SOURCE_BRIDGE_RESIDUAL`, `CROSS_CHANNEL`, `EULERIAN_OPERAND`,
  `MATERIAL_OPERAND`, `DIMENSIONS`, `FROZEN_RELATIONS`, `ADVECTION_ABSENCE`. (`*_NODES` are the SymPy symbolic
  siblings — pair them to the WL `ARITHMETIC` field of the same object, not to a separate WL tag.)
- **Covariance (SymPy `S11CC2_N6COV_*` ↔ WL `WL_S11CC2_N6COV_*`):** `R_COV`, `R_COV_BASELINE`, `R_COV_INCREMENT`,
  **`R_COV_CONTROL_DELTA`**, `SOURCE_ACTUAL`, `SOURCE_PREDICTED`, `SOURCE_BASELINE`, **`SOURCE_CONTROL_DELTA`**,
  `FROZEN_PHI` (symbolic+support); **`PHI_DOMAIN_CENSUS`**, **`ACTUAL_CONTROL_PARAMETERS`** (STRUCTURAL join — the
  substitution-domain coverage guard + the live control settings; seal any digest subfields). ⚠ Both engines emit ALL
  of these (WL `covNames`; the SymPy covariance stream) — they are the covariance independence control + its residual
  + the domain census, ⛔ NOT bookkeeping. (SymPy `PRIMES`/`PIT_PROVENANCE` are dedicated tags; WL's live inside
  `WL_S11CC2_N6_LOCAL_PROBE` — surface each per-engine, ⛔ do NOT join/difference; `PROVENANCE` surfaced-not-joined,
  digest sealed.)
- **Guards (both):** `SLOT_GUARD_{NATIVE,CARRIER,RESIDUAL}`, `CLOSURE_GUARD_{NATIVE,CARRIER,RESIDUAL}` — structural
  (⚠ SymPy guards carry a `probe` axis, WL embeds the route in the slot key — see the loader duplicate rule).
Join by object name after stripping the `WL_` prefix and mapping `S11CC2_` ↔ `WL_S11CC2_`; ⛔ a family present in
one stream only is an accounting row, never a silent skip.

## Inputs (read-only, committed; ⛔ `datalad get` the WL `.out` before reading — it is a git-annex pointer)
- **WL** `research/pde_ledger_v3/mathematica/out/S11c_c2_N6_mathematica_audit.out` (~400 MB, **40 tags, one per
  line**; each `WL_S11CC2_<NAME> = <| … |>` payload keyed INTERNALLY by the 4 case-tuples
  `{"LAB_HELD"|"MATERIAL_ADVECTED","RHO4_CONSTANT"|"RHOBR_CONSTANT"}`; ⛔ do not run the WL engine). ⚠ Read it
  LAZILY (line-addressed, per case-key) — do NOT materialize all 40 tags before the first compare.
- **SymPy** (one JSON object per line, `{"object":"S11CC2_<NAME>", "anchoring":"LAB_HELD"|"MATERIAL_ADVECTED",
  "density":"RHO4_CONSTANT"|"RHOBR_CONSTANT", "data":{…}, …}`; ⚠ the case is TOP-LEVEL `anchoring`+`density`
  fields, ⛔ NOT a `CASE` array; each line is ONE case, all 4 cases concatenated per file):
  - `research/pde_ledger_v3/scripts/out/S11c_c2_N6_covariance_sympy.out` (`S11CC2_N6COV_*`)
  - `research/pde_ledger_v3/scripts/out/S11c_c2_N6_reconcile_sympy.out` (`S11CC2_N6RC_*`)
  - `research/pde_ledger_v3/scripts/out/S11c_c2_N6_diagnostic_sympy.out` (`S11CC2_N6_*` + SymPy-only set)
  ⛔ do not run the SymPy engines; ⛔ the `*_exports.py` are LEDGER deltas, NOT tag streams.
- **The SymPy engine sources** `scripts/S11c_c2_N6_{covariance,reconcile,diagnostic}_sympy.py` are the authority
  for what each `data.columns` slot MEANS and the `_NODES` DAG schema — read them to type keys + reconstruct the
  DAG; ⛔ do not guess a slot.
- **Spec** `directives/S11c_c2_SHARED_PHYSICS.md` (§1/§1b supplied inputs + AGREE/UNDECIDED; §5c the two N6 routes;
  §N8/T7 the comparator contract).

## Mechanical base to re-express (verified sound — reuse, do not re-invent)
`research/pde_ledger_v3/scripts/S11c_c1_cross_engine_comparator.py` **and** its base
`S11b_cross_engine_comparator.py` / `S11c_a_cross_engine_comparator.py` (study them + the c1 test). Reuse:
- ⛔ **c1's `load_wl` does NOT apply** — it recognizes the colon grammar `WL_S11CC1_<Q>:` and splits on `:`. N6 WL is
  a **single-line** `WL_S11CC2_<NAME> = <| … |>`, one tag per line, ` = ` delimiter. Write a **new N6 WL loader**
  (prefix `WL_S11CC2_`, ` = ` split, one tag/line, lazy line-addressing). Write a **new `load_py_jsonl`** for the
  SymPy streams keyed on `object` + TOP-LEVEL `anchoring` + `density` (⛔ there is NO `CASE` field). ⚠ **Duplicate
  identity must include the `probe` axis** — valid SymPy guard rows repeat the same `object`/`anchoring`/`density`
  with `probe="EULERIAN"` and `probe="MATERIAL"` (two legitimate routes); keying duplicates on
  object/anchoring/density alone WRONGLY rejects them. Map the SymPy guard `probe` label to the WL embedded route axis.
- the WL structural readers (`split_top`, `wl_assoc_pairs`, `wl_field`, `preprocess_wl`) and `parse_wl_value` /
  `canonical_basic` (Wolfram `Inactive[…]` → SymPy) — this is your WL symbolic translator.
- ⭐ **NEW: a SymPy `_NODES` DAG reader — reconstruct via the family's `ARITHMETIC_DAG`.** ⚠ `<TAG>_NODES` is NOT
  self-contained: its `root_nodes` are copies of the ROOT records only; every non-leaf `arg` is `{"ref": i}` into the
  shared **`ARITHMETIC_DAG.nodes`** list (the family's `*_ARITHMETIC_DAG` tag). Reconstruct as
  `ARITHMETIC_DAG.nodes[root_id]` following `ref` indices; op vocabulary = `number|variable|add|mul|pow`. ⚠⚠ **Literals
  are NOT all plain names** — applied functions/derivatives are `(sp.srepr(e), "formal_jet")` string literals (plus
  `global`/local-point/`algebraic_i`/integral-guard variants), so **"no srepr here" is FALSE**: decode a `formal_jet`
  literal via an `srepr`/`sympify` round-trip, ⛔ never treat it as an opaque name (that manufactures a symbolic
  disagreement). Build it from the emitters + the node constructors, ⛔ not by guessing. ⚠ Diagnostic GUARDS carry PIT
  + NO `_NODES` (symbolic is WL-only there); a missing `ARITHMETIC_DAG` for a joined family is `parse_failed`, ⛔ never
  a zero residual.
- the typed recursive `residual` with its three-valued returns; `mechanical_lower_camel` applied **INJECTIVELY**
  (check injectivity FIRST — a collision is a finding, never a silent merge).
- the S11c-b/c1 **lazy `build_case`/`materialize`/`release_case`** memory discipline (WL payloads are large).
- the per-family `Accounting`; the `BoundIntegral`/rational canon.
⚠ **DEFINE this comparator's OWN axis vocabulary + typed `make_key`/`decode_*`** — ⛔ do NOT literally reuse the
c1 geometry vocabulary (FACE/DIRECTION/PARITY/REGIME): N6 keys on `ANCHORING ∈ {LAB_HELD, MATERIAL_ADVECTED}`,
`DENSITY ∈ {RHO4_CONSTANT, RHOBR_CONSTANT}`, plus each family's inner container axes (channel/column/slot). Type
each token to its own axis; ⛔ never positional-guess, ⛔ never merge distinct axes. ⛔ Do not reuse any
classification / verdict / `main`-status machinery.

## Reconciliation folds — pre-register MECHANICAL name/CAS spelling ONLY; ⛔ physics-bearing rep differences get a COMPUTED residual (DO-NOT-FOLD), never a pre-registered collapse
1. **Inherit the S11c-a/b/c1 name/CAS map** (`Inactive[Equal]→HeldEqual`, `Inactive[Integrate]→sp.Integral`,
   rational `expand→cancel(together)`, the bare-symbol `mechanical_lower_camel` inverses re-verified against the
   real SymPy symbols). ⚠ Keep any `Inactive[Greater]`/`Inactive[FourierTransform]` UNEVALUATED (held heads); ⛔
   do not evaluate a held predicate to a native boolean.
2. **⛔ DO-NOT-FOLD (surface both operands + PRINT the three-valued residual; ⛔ never collapse to identity) — the
   physics-bearing representation differences (rule 5/6 / §N8/T7 / [[feedback_reconcile_representational_bridge]]):**
   - **The blind-carrier representation.** WL `CARRIER_EULERIAN`/`C_E` (blind graph-geometry) vs SymPy imported-slab
     `C_E`: SAME object, DIFFERENT representation — surface both AND **compute the symbolic residual** (⛔ do NOT skip
     it); whether they coincide is the T7 reconcile's post-run job, but the residual is COMPUTED here.
   - **`FROZEN_PHI`** (WL Wolfram substitution map vs SymPy `substitution_map` expr-STRINGS) — surface both + compute
     the map residual (translate where tractable, else `deferred_oversize`).
   - **Any anchoring/density/representation the two engines key or shape differently** (§1b) — type the literal
     token; ⛔ do not declare one a "spelling" of the other; surface the shape difference as accounting/residual +
     the pairing table.
   ⛔⛔ **The TRUE SEALS (⛔ NEVER differenced at all) are ONLY** the PIT residue tables + prime sets + sample points,
   and the digests (the SEAL section above). ⛔ Do NOT add the carrier / `FROZEN_PHI` / rep-differences to the true-seal
   list — those get a computed residual.
3. **CONTROL / guard families compared as emitted, reaching every addressable leaf** (`SLOT_GUARD_*`,
   `CLOSURE_GUARD_*`, `SPLIT_CHECK`, `SPLIT_SUM`): extract their operand + residual leaves and compare A / B /
   A−B, **no target**. ⛔ NEVER blanket-collapse `X(args)→X` to make a control agree — hand-code the extractor and
   FLAG.
4. **Supplied / bookkeeping / SymPy-only** (`_LOCAL_*`, `PROVENANCE`, `PREMISES`, the SymPy-only diagnostic set,
   any supplied §1 premise): excluded from the join / emitted as accounting; ⛔ do not broadcast a premise across
   the other engine's cases.

## Per-family extraction — ⛔ DISCOVER each family's ACTUAL container from THIS payload (recon traps below)
For every joined family, sample the real WL association head AND the real SymPy `data`/`_NODES` shape from these
`.out`, then write the extractor to reach the paired leaves. Known traps (verify each against the payload):
- **`DIMENSIONS` — THREE different containers; ⛔ do NOT read it from one place.** (i) SymPy per-column `[L,T,M]` is
  the **TOP-LEVEL `dimension` array** on EVERY numeric emission, index-aligned with `data.columns` (⛔ NOT inside
  `data`, ⛔ NOT the standalone `N6RC_DIMENSIONS`) — zip `dimension[i]` with `data.columns[i]`; there is NO
  `N6COV_DIMENSIONS`. (ii) WL per-leaf `[L,T,M]` lives in each leaf's own `DIMENSIONS` field (beside its `ARITHMETIC`).
  (iii) the standalone `N6RC_DIMENSIONS` tag is a THIRD, broader nested census (WL `OBJECTS → object → slot →
  {…SUPPORT…, ADDITION_CONSISTENCY}`; SymPy `data.{eulerian,material} → "('E_W',…)" →
  {computed,consistent,target,unknown,zero}`) — an ADDITIONAL guard family, ⛔ NOT the sole dimensions source. Compare
  the `[L,T,M]` vectors on matched columns; FLAG container mismatches.
- **`R_N6` / `R_COV` / residual families** — WL: `ARITHMETIC → Inactive[…]` (symbolic) + `PROBE_NUMERATORS`/
  `PROBE_DENOMINATORS` SparseArrays (PIT — SEALED); SymPy: `<TAG>_NODES` via `ARITHMETIC_DAG` (symbolic) +
  `data.numerator_denominator` + `data.columns` + `data.nonzero_modular_numerator` (PIT — SEALED). Compare the SYMBOLIC
  (translate, three-valued, defer only over budget) + the STRUCTURAL support ON MATCHED KEYS ONLY (the pre-registered
  schema bridge + pairing table; support typed `NONZERO_WITNESSED`/`NO_NONZERO_FOUND`). ⛔ never difference the raw
  PIT residues.
- **Carrier/source operands** (`CARRIER_EULERIAN/MATERIAL`, `SOURCE_*`, `EULERIAN/MATERIAL_OPERAND`, `*_CHANNEL`,
  `CROSS_CHANNEL`) — symbolic (WL `ARITHMETIC` vs SymPy `_NODES`) + support on matched keys. ⭐ The carrier-rep
  difference gets a COMPUTED three-valued residual (DO-NOT-FOLD), ⛔ it is NOT a true seal.
- **`FROZEN_PHI`** — WL `MAP → {var → Wolfram expr}`; SymPy `data.substitution_map = [[var, expr-string], …]`.
  Translate where tractable; else `deferred_oversize`.
- **`PRIMES`/`PIT_PROVENANCE`** — SymPy dedicated tags; WL inside `N6_LOCAL_PROBE`. Surface each per-engine (this
  is where the different prime sets are VISIBLE); ⛔ do not join/difference.
- `_LOCAL_*` excluded; emit each engine's local-tag inventory so the exclusion is visible.

## Tests (SEPARATE file; synthetic only; ⛔ do not load or run either measured engine)
Model on `test_S11c_c1_cross_engine_comparator.py`. Require a **source-derived synthetic fixture for every joined
SCHEMA**. Include:
- typed-key construction rejects a duplicate axis and an unknown/untyped axis; ANCHORING and DENSITY stay distinct;
- **WL SINGLE-LINE ` = ` grammar** parse AND SymPy `_NODES` round-trip through `ARITHMETIC_DAG` refs — INCLUDING a real
  `(srepr, "formal_jet")` literal that must round-trip (⛔ not become an opaque name); a missing `ARITHMETIC_DAG` for a
  joined family → `parse_failed`, ⛔ not a zero residual;
- **the `probe`-axis duplicate rule**: two SymPy guard rows with the same object/anchoring/density but
  `probe`=EULERIAN/MATERIAL are BOTH kept, ⛔ not rejected as duplicates;
- injectivity of the mechanical name/CAS map; the schema bridge pre-registers ONLY name/CAS identities (grade,
  `SUM↔0`) — a synthetic block/kernel-family NAME mismatch surfaces as `sympy_only_column`/`wl_only_slot` in the
  pairing table, ⛔ NOT auto-joined;
- **the PIT SEAL is enforced**: a synthetic pair of PIT residue tables with different primes/sample-counts is NEVER
  differenced — each engine's support surfaced + a `pit_sealed` marker, ⛔ not an A−B of the residues;
- **support is one-sided**: an all-zero SymPy `nonzero_modular_numerator` → `NO_NONZERO_FOUND` propagating as
  **UNDECIDED** (⛔ not "structurally zero", ⛔ not subtracted against the other engine's boolean); a nonzero column →
  `NONZERO_WITNESSED`;
- **carrier-rep computes a residual**: a synthetic blind-vs-imported `C_E` pair prints `A_minus_B` (⛔ NOT skipped as a seal);
- **three-valued residual preserved**: native boolean → `BooleanNotResidualable`; undecided/oversize leaf → an emitted
  `UndecidedResidual`/`deferred_oversize` DISTINCT from a parse failure;
- **the repoint ablation** (model on `S10_cross_engine_comparator_repoint_ablation.py`): a DIFFERENT synthetic object's
  payload under a previously-paired NAME → the symbolic/structural residual MOVES (flips away from agreement). ⛔ A
  symbol/spelling rename is NOT a repoint;
- **a SymPy-only family → `sympy_only`, never joined** (a synthetic `MU_RECONSTRUCTION_RESIDUAL` with no WL sibling →
  accounting row, ⛔ no manufactured WL operand);
- ⭐ **zero-extract from two NONEMPTY declared-shared containers is an OPERATIONAL EXTRACTOR FAILURE** — a synthetic
  fixture where both engines emit a joined family but the extractor reaches 0 leaves must trip the failure path, ⛔ NOT
  read as agreement or disagreement;
- disagreement is a measurement: a synthetic SymPy≠WL symbolic pair prints `operand_A`/`operand_B`/`A_minus_B` +
  `ACCOUNTING` and exits 0.

## Definition of done (the re-review legs check these empirically)
Every joined family prints its `ACCOUNTING` line. ⛔ **No family silently extracts 0**; every family either joins
cases OR emits documented `axis_set_mismatch`/`sympy_only`/`wl_only`/`deferred_oversize`/`unmatched_key` rows with a
reason. ⛔ **Do NOT use `join>0` as an exit criterion** (rule 5); report a per-family **extracted-leaf count**.
⭐⭐ **Zero extracted leaves from two NONEMPTY, declared-shared containers is an OPERATIONAL EXTRACTOR FAILURE
(nonzero exit) — ⛔ NOT a physics disagreement and ⛔ NOT agreement.** The coverage rule must distinguish a cautious
comparator from a nonfunctional one; the source-derived synthetic fixtures (test file) make a broken extractor fail a
TEST, not slip through the run. Prints `operand_A`, `operand_B`, `A_minus_B` before any guard; asserts nothing on
measured payloads; exits 0 on disagreement, nonzero on operational failure. ⭐ **State the per-object budget
explicitly** (name RSS + wallclock ceilings): a single object over budget → `deferred_oversize`; N6's largest operand
is ~80 MB, so deferral is RARE — an **all-`deferred_oversize` run is an operational failure, ⛔ not a pass**. A
`RUN_ACCOUNTING` summary reports `families`, `families_with_join`, `families_with_unpaired`, `parse_failed`,
`deferred_oversize`, `zero_extract_failures`, `runtime_seconds`; a `MEASUREMENT_SCOPE` line records §1–2 + the supplied
substrate as supplied/unfalsifiable, the PIT residues + digests + SymPy-only families as surfaced-not-joined, and
`residual_target=none`. `LOCAL_INVENTORY` lines make the `_LOCAL_` exclusion visible for both engines. ⚠ **This N6
instrument does NOT carry the main-c2 §3d representation questions** (`t_s` scalar-vs-4-vector, DtN whole-form,
flat-symbol) — those families are ABSENT from the N6 streams (they belong to the c2 self-energy comparator); ⛔ do NOT
require the builder to surface them here.

⚠ **Memory / runtime.** WL payloads are large (~400 MB stream); keep the lazy materialize/release discipline. If a
symbolic difference OOMs or exceeds the stated per-object budget, that is a `deferred_oversize` finding (recorded
out-of-band), ⛔ NOT a reason to narrow the comparison or drop a family — report it.

## Builder report (≤30 lines)
Per-family accounting + extracted-leaf-count summary; any family you deferred (with the size/time reason); the
SymPy `_NODES` DAG op-vocabulary you reconstructed and where you verified it; which spelling folds you added and
the injectivity result; how you enforced the PIT/digest/SymPy-only SEALS; the repoint-ablation result; runtime +
peak RSS. State that §1–2 + the supplied substrate are supplied/unfalsifiable and that **no residual target was
given** (rule 5), and that the PIT residues + digests are surfaced-not-differenced (blindness).
