# S11c-c2 N6 ablation-harness directive — output-contract + WL-budget revision, re-review record

## Artifact
`directives/S11c_c2_N6_ablation_harness_directive.md`, **revision-1** (output-contract → compact PIT+digest;
WL budget → single-case pin + `drawCount=4`). Built on the knife-list cleared at commit `096d27a9`.
Revision-1 authored by Codex-sol (`gpt-5.6-sol` xhigh); the knife DESIGN is frozen from `096d27a9`.

## Role
Directive / build spec. Reviewed by reading + mechanical checks; the four harnesses do not exist yet
(executable ablation is the later build review — E2/L-R14). Spec ⇒ **review-until-clear**.

## Legs (Codex-authored directive → fresh Claude + Grok)
- **Fresh Opus agent** (`general-purpose`, model opus) — report: `_legs/S11c_c2_N6_rereview_agent_report.md`.
- **Grok** (`grok-4.6` high) — report: `_legs/S11c_c2_N6_rereview_grok_report.md`.
- Identical rendered prompt: `_legs/S11c_c2_N6_output_contract_rereview_prompt.md`; PRE→POST diff handed to
  both: `_legs/S11c_c2_N6_rereview.diff`.
- (Leg model = Opus per user 2026-09-08; the initial Sonnet agent was stopped and relaunched on Opus. The
  Opus leg did NOT trip the prior `[bio]` filter and produced a deep review.)

## Verdicts
- **Opus: SOUND** (Q1–Q5), with a retracted finding (nearly flagged `emitted_object_sha256` as a fictional
  engine field, then verified it is harness-computed), one non-blocking Q3 sampler caveat, two clarity nits.
- **Grok: NOT-SOUND** — three findings (K_junk inert under the pin; WL digest against nonexistent engine
  fields; `R_N6`/`R_cov` "direct payload" leak path).

The legs DISAGREED on the overall verdict (M1 — the disagreement is the measurement). Adjudicated by
orchestrator verification against the engine source (G4). **Grok's NOT-SOUND is correct**; all three findings
are real. Opus's deep pass was valuable but its SOUND verdict missed findings 1 and 3.

## Findings — orchestrator verification (each confirmed against the engine)

**F1 — WL K_junk is INERT under the `LAB_HELD/RHOBR_CONSTANT` pin (design defect in revision-1). CONFIRMED.**
- `audit.wl:21` `actualJunkCase = {"MATERIAL_ADVECTED", "RHO4_CONSTANT"};`
- `audit.wl:845` `activeJunk = If[case === actualJunkCase, actualJunkCoefficient, 0];` — the ONLY `case ===`
  gate in the file. Under the revision's `LAB_HELD/RHOBR_CONSTANT` pin, `activeJunk = 0` ⇒ the K_junk knife
  (`actualJunkCoefficient 0→1`) does nothing: WL K_junk FORM ≡ IDENTITY ≡ ×2 → bite invisible.
- K_carrier (`materialNormalKnife` @860) and K_split_route (`sourceChannel` @892) are NOT case-gated.
- ⇒ The WL pin is FORCED to `{"MATERIAL_ADVECTED","RHO4_CONSTANT"}` (= actualJunkCase) — the unique case
  where all three WL knives are simultaneously live. The "match the SymPy pin" rationale was wrong; the WL
  harness is an independent cert. (SymPy covariance `ACTUAL_JUNK` is passed unconditionally at
  covariance_sympy.py:187-188 — NOT case-gated — so the three SymPy harnesses stay on `LAB_HELD/RHOBR`.)

**F2 — `R_N6`/`R_cov`/`SPLIT_CHECK`/guard residuals are HEAVY PIT objects, not "small scalar direct
payloads" (leak path + self-contradiction). CONFIRMED.**
- WL wraps every numeric object, incl. `R_N6` (`addNumeric["RC","R_N6",...]` @897), in
  `"ARITHMETIC" -> (circuitExpression /@ values)` @804 — a heavy symbolic circuit.
- SymPy `add('R_N6', R)` (reconcile:270), `('R_COV', rcov)` (covariance:201) are residual TABLES passed to
  `n.pit`. "18/288", "0" are their PIT-projected fingerprints, not the raw tags.
- Revision-1 clause 1's "these remain direct payloads" is false, contradicts its own coverage note (routes
  them through `n.pit`), and is a path back to the ~285 MiB blowup. ⇒ Remove the carve-out; route them
  through the compact PIT fingerprint + digest like every other certified object.

**F3 — the per-object digest is HARNESS-computed, not an engine field (builder ambiguity). CONFIRMED.**
- `emitted_object_sha256`/`RUN_PROVENANCE` = 0 hits in the WL engine. They are the harness driver's records:
  `ablation_harness_wl.py:190` `emitted_hashes[name] = sha256(<full emitted payload text>)`. Because it
  hashes the FULL pre-projection payload (ARITHMETIC + PROBE_NUMERATORS), it DOES witness a same-support FORM
  change — the load-bearing backstop for a fingerprint that does not move. ⇒ Clause 1 must state the digest
  is a harness-computed SHA-256 of the full emitted payload before compaction (WL: emitted payload text;
  SymPy: `n.sha` over the parsed emitted PIT payload — raw objects are DAG Nodes, not `sympify`-able), hashed
  then discarded; digest only in the transcript.

**Non-blocking (recorded for the executable BUILD review):** Opus Q3 — the joint-rejection PIT sampler
(audit.wl:774–791) can shift the accepted sample set when a knife changes singularity structure, so a
DEAD/one-sided control's DIGEST may move even when its object is unchanged: a possible false POSITIVE on DEAD
controls (never a false negative; the nonzero fingerprint stays robust). Interpretation caveat for the build
review, not a directive defect.

## Disposition
Revision-1 = **NOT-SOUND**, 3 findings unresolved at this commit (preservation baseline, ⛔ not acceptance).
FOLD → revision-2 delegated to Codex-sol (keeps the directive Codex-authored; re-review stays fresh-Opus +
Grok): fix F1 (repin WL to actualJunkCase), F2 (remove R_N6/R_cov carve-out), F3 (digest = harness full-
payload SHA), + record the Q3 caveat. Then re-review-until-clear.
