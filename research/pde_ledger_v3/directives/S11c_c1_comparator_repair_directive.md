# Repair brief — S11c-c1 T7 cross-engine comparator (baseline 7141e6ad; two surgical fixes)

The comparator at `research/pde_ledger_v3/scripts/S11c_c1_cross_engine_comparator.py` was BUILT and
re-reviewed SOUND (two legs CLEAN — no false agreement, no hidden disagreement, every seal load-bearing). This
repair applies exactly two fixes surfaced by the re-review; ⛔ it changes NO soundness behaviour, NO seal, and
NO measurement semantics. Baseline commit **7141e6ad**. It computes and PRINTS, deciding nothing (rule 2) — that
is unchanged. ⛔ No verdict/status token; ⛔ no residual target on any MEASURED stream (rule 5).

## R1 — add the `cS0 ← c_s0` mechanical spelling fold (a legitimate inverse the build wrongly excluded)
The build's directive claimed `mechanical_lower_camel("c_s0") == "cs0"`; that is FALSE. Verify yourself:
`python3 -c "import S11b_cross_engine_comparator as m; print(m.mechanical_lower_camel('c_s0'))"` → **`cS0`**.
So WL `cS0` is the EXACT mechanical inverse of the real PY `Symbol('c_s0')` (present in the PY stream ~44k
times) — the SAME kind as the 11 bare-symbol folds already in the c1 map (`rhoM←rho_m` … `tauV←tau_V`).
- **Add `cS0 ← c_s0` to the c1 bare-symbol spelling map** (`C1_BARE_SYMBOL` / whatever the map is named),
  alongside the existing 11, so WL `cS0` normalizes to PY `c_s0`.
- ⛔ It must obey the SAME activation gating as the other 11 (fire ONLY when the real PY `Symbol('c_s0')` is
  present — the existing `verify_active_spelling_map` / active-fold mechanism), and it must keep the injectivity
  check passing (`checked_mechanical_symbol_map`, collisions=0). `c_s0`/`cS0` are BARE symbols on both sides —
  ⛔ this is NOT an applied head and must NOT arg-strip anything.
- Rationale (do NOT encode as a target): `cS0` vs `c_s0` is a pure spelling difference for the one bulk sound
  speed; leaving it unfolded injects a spurious `c_s0 − cS0` remainder into any joined leaf carrying both. It is
  currently latent (fires in no joined residual on the feasible families) but WL `cS0` is dense in the deferred
  PERMEABLE_* families, so the fold must be in place before that heavy run.
- (The spec-level source of the error — the `cS0` paragraph in the build directive — has ALREADY been corrected
  by the orchestrator; ⛔ you do NOT edit any directive. Your job is the code fix only.)

## R2 — parse the WL multi-range held-integral energy leaves so the operand is DISPLAYED, not `PARSE_FAILED`
⚠ **CORRECTED FORM (both decision legs — the earlier "range-less 1-arg" description was WRONG).** The WL parser
raises `InputError: Inactive[Integrate] does not have integrand and range` because the inherited
`S11c_a_cross_engine_comparator.preprocess_wl` (`:671-682`) accepts **exactly ONE** range, but every measured c1
WL integral is a **4-argument, TRIPLE-range** volume integral:
`Inactive[Integrate][<integrand>, {xOne,-Infinity,Infinity}, {xTwo,-Infinity,Infinity}, {xThree,-Infinity,Infinity}]`,
and the nested form is `Inactive[Limit][Inactive[Integrate][…4-arg…], outwardDistance -> Infinity]`. There are
**100** such integrals across **THREE** ENERGY families (verify: `ENERGY_RESIDUAL` 80, `ENERGY_BULK_FARFIELD_FLUX_OPERAND`
16, `ENERGY_FACE_TRACTION_OPERAND` 4). They are in UNPAIRED (join=0) leaves (PY ENERGY keys carry a `SCENARIO`
axis, WL do not) and are surfaced honestly as `<PARSE_FAILED>` — so this is a **display/completeness** fix, ⛔ NOT
a soundness fix and ⛔ NOT a residual change (it changes the parse-failure COUNTS in accounting, nothing else).
- Parse the **4-arg (integrand + three ranges)** `Inactive[Integrate]` and the outer `Inactive[Limit][…,
  outwardDistance -> Infinity]` as **held (unevaluated) heads** — the same pattern the comparator already uses
  for `Inactive[Greater]→HeldInactiveGreater` (a `HeldInactiveIntegrate` / `HeldInactiveLimit` `AppliedUndef`
  carrying its raw arguments: the integrand, the three range tuples, and the limit rule). ⛔ Do NOT evaluate the
  integral/limit; ⛔ do NOT invent or drop a range; keep the operand as an as-emitted held object.
- ⛔⛔ **c1 PREPROCESSOR ONLY — do NOT widen `S11c_a_cross_engine_comparator.preprocess_wl`.** The inherited
  **1-range** `Inactive[Integrate][f,{x,a,b}] → BOUND_INTEGRAL` canonicalization MUST stay intact
  (`test_bound_integral_keeps_binder_and_bounds`, test file `:223-234`), and `Inactive[Equal]→HeldEqual` MUST
  stay intact. Hold ONLY the unsupported multi-range Integrate and its nested Limit, in the c1 preprocessor,
  before/around the existing c1 `HeldInactiveGreater` rewrite — ⛔ never a blanket "hold every `Inactive` head".
- ⚠ The comparator's job is only to display what the frozen WL engine emitted, not to repair it. Held-parsing
  is the faithful representation; ⛔ do not reconcile it against the PY operand by any fold.

## What must NOT change (verify against the diff)
- The protected core is BYTE-IDENTICAL except the spelling-map addition (R1) and the WL held-parse handling
  (R2). ⛔ No change to: `make_key`/axis typing, the 5 SEALS (two-momentum legs, μ_θ, ω-assumption, background
  density field, regime/parity — all must stay unreconciled and load-bearing), the `raw_control_case` whitelist
  `{DTN_OPERATOR.WHOLE_OBJECT, NONINVERTIBILITY_CONDITION.OPERATOR}`, the three-valued residual objects, the
  per-family accounting, `RUN_ACCOUNTING`/`MEASUREMENT_SCOPE`/`LOCAL_INVENTORY`, exit-0-on-disagreement.
- ⛔ Do NOT add any other spelling fold, do NOT touch `qOut`/`kOne`/`muTheta`/`rhoBrBgRho4Constant`/`w1Profile`/
  `w1Jet` (they STAY applied/sealed), do NOT re-scope `raw_control_case`.

## Tests (extend `test_S11c_c1_cross_engine_comparator.py`; synthetic only; ⛔ do not run either measured engine)
Keep every existing test green (the sibling combined suite included) EXCEPT the one that R1 intentionally
changes: **`test_noninherited_sound_speed_spelling_is_exposed` (test file `:214-219`) asserts
`parse_wl_value("cS0") == Symbol('cS0')`, which R1 now folds to `c_s0` — REPLACE that test** with the R1
activation tests below (⛔ do not keep it, and ⛔ do not weaken the activation gate to preserve it). The per-family
accounting SCHEMA is unchanged; only the observed parse-failure COUNTS change (expected from R2). Add:
- **R1 activation + injectivity:** ⚠ the run-time gate is ALL-OR-NONE (`REQUIRED_C1_PY_SYMBOLS <= names` at
  `:2127`), so the positive fixture's PY stream must contain the **complete 12-symbol** required vocabulary (the
  11 existing + `c_s0`), not just `c_s0`. With that vocabulary present, a WL leaf `cS0` folds so a PY `c_s0` /
  WL `cS0` pair residuals to `Integer(0)`; ⛔ a WL leaf `cS0` with the required PY vocabulary ABSENT does NOT
  fold (stays `Symbol('cS0')`); the injectivity check reports 0 collisions with the 12th entry.
- **R1 no arg-strip regression:** the existing seal tests (qOut/rhoBrBg/w1* stay applied) remain green — adding
  `cS0` opened no applied-head strip.
- **R2 held parse (the REAL form):** a synthetic WL 4-arg `Inactive[Integrate][f[xOne,xTwo,xThree],
  {xOne,-Infinity,Infinity},{xTwo,-Infinity,Infinity},{xThree,-Infinity,Infinity}]` and the nested
  `Inactive[Limit][<that>, outwardDistance -> Infinity]` parse to a HELD head (an `AppliedUndef`, NOT an
  evaluated integral/limit, NOT a native value, NOT a `PARSE_FAILED`), carrying the integrand, all three range
  tuples, and the limit rule.
- **R2 preserves the inherited 1-range path:** a synthetic 1-range `Inactive[Integrate][f[w],{w,0,2}]` still
  canonicalizes to the inherited `BOUND_INTEGRAL` (binder + bounds retained) — R2 did not widen or break it.

## Builder report (≤20 lines)
The literal `mechanical_lower_camel("c_s0")` output; the diff scope (which lines changed for R1 vs R2, and that
R2 is c1-preprocessor-only with S11c-a's `preprocess_wl` untouched); the injectivity result with 12 folds; the
test result (all green with `test_noninherited_sound_speed_spelling_is_exposed` REPLACED + the new
activation/held-parse/BoundIntegral-preservation tests); confirmation the protected core + the 5 seals + the raw
whitelist are byte-identical (name the regions). ⛔ No residual target given.
