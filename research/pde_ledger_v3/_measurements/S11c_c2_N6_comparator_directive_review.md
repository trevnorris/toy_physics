# S11c-c2 N6 cross-engine comparator BUILD DIRECTIVE — decision-leg gate + fold record

## Artifact / role
`directives/S11c_c2_N6_comparator_build_directive.md` — orchestrator-written build directive (decision-list class)
for the N6 cross-engine measurement instrument `scripts/S11c_c2_N6_cross_engine_comparator.py`. Gate = the
pre-builder decision review (G2 / [[feedback_directive_design_review]]): ⛔ no builder until 2 decision legs.

## Handed inputs (all committed)
WL N6 `.out` (`ae73b884`); 3 SymPy N6 `.out`s (`7c0790ab`, baselines reproduced: R_N6=18/288 nonzero in 3 cases,
R_cov=0 ×4, SPLIT_CHECK=0 ×4); c2 `SHARED_PHYSICS.md` §N8/T7 + §5c + §1b; the c1 comparator + base
(`S11c_c1/S11b/S11c_a_cross_engine_comparator.py`).

## Legs (orchestrator-written → Codex-sol + Grok; identical prompt; on sight)
- **Codex-sol** (`gpt-5.6-sol` xhigh) — report `directives/_legs/S11c_c2_N6_comparator_directive_review_codex.md`.
- **Grok** (`grok-4.6` high) — report `directives/_legs/S11c_c2_N6_comparator_directive_review_grok.md`.
- Prompt `directives/_legs/S11c_c2_N6_comparator_directive_review_prompt.md`. Both ran on the real payloads.

## Verdict: BOTH NOT-SOUND — strongly CONVERGENT (independently found the same core defects). Verified + folded ONE PASS.
Both legs affirmed the load-bearing design decisions as SOUND: the PIT SEAL (the engines' deliberately-different
prime sets + sample points make PIT residues non-cross-differenceable — a required clause, ⛔ not a spec-currency
error); three-valued / no-verdict; the authority trio (§N8/T7 + §5c + §1b); the reuse of `residual`/`canonical_basic`/
lazy cases; the NEW `_NODES` reader necessity; the named container traps; SymPy-only exclusion.

### Findings folded (all verified against the real sources, G4)
- **A** (both) — join set dropped 4 SHARED covariance families: `R_COV_CONTROL_DELTA`, `SOURCE_CONTROL_DELTA`
  (symbolic+support), `PHI_DOMAIN_CENSUS`, `ACTUAL_CONTROL_PARAMETERS` (structural). ADDED.
- **B** (both) — the column↔slot bridge must be a PRE-REGISTERED mechanical schema-normalization + typed decoder,
  ⛔ not "post-run adjudication". ⚠ LEG-SPLIT: Codex wanted block-name + kernel-family pre-registered; Grok wanted
  ONLY name/CAS identities + a pairing table. **Adjudicated to Grok (rule 5):** pre-register ONLY grade `{1,η,σ}↔(η,σ)`
  / face `SUM↔0` / axis-order / `_NODES` literal decoding; ⛔ NEVER pre-register block/kernel-family/component-vs-formal
  equality (a physics identification = the forbidden blanket-collapse [[feedback_handcode_comparison_never_blanket_collapse]]);
  emit a `{matched, sympy_only_column, wl_only_slot}` PAIRING TABLE, compare only matched keys, unmatched = residual.
  Added a dedicated "THE SCHEMA BRIDGE" section.
- **C** (both) — loader grammar: c1 `load_wl` (colon-split multiline) does NOT parse N6 (single-line ` = <|…|>`); +
  duplicate identity must include the `probe` axis (SymPy guards repeat obj/anchoring/density with probe=EULERIAN/MATERIAL).
  Removed the `CASE` field from the loader contract (already fixed in Inputs). REWRITTEN.
- **D** (both) — per-column `[L,T,M]` is the SymPy TOP-LEVEL `dimension` array (zip with `data.columns`), not
  `N6RC_DIMENSIONS`; three distinct dimension containers named.
- **E** (Codex) — nonzero-support is ONE-SIDED evidence: type `NONZERO_WITNESSED` vs `NO_NONZERO_FOUND`(=UNDECIDED);
  ⛔ never subtract native support booleans [[feedback_instrument_claims_on_invariants]].
- **F** (both) — `_NODES` reconstruct via `ARITHMETIC_DAG.nodes[ref]`; literals include `(srepr,"formal_jet")` etc.
  ("no srepr here" was FALSE) → round-trip; missing `ARITHMETIC_DAG` = parse_failed, not zero.
- **G** (both) — RULE-5 LEAK: dropped "REP_INVARIANCE_RESIDUAL nonzero is EXPECTED"; sympy_only reason = "no WL
  sibling" only; completed the sympy_only enumeration. DoD: zero-extract from two NONEMPTY shared containers =
  OPERATIONAL FAILURE (not agreement/disagreement); source-derived synthetic fixture per joined schema.
- **H** (both) — §N8/T7's §3d-surfacing + ≥64 GB deferral concern the c2 self-energy comparator; the N6 streams carry
  NEITHER (no t_s/DtN/flat-symbol families; largest N6 operand ~80 MB). Added N6-SCOPING clause + explicit per-object budget.
- **I** (Grok F5) — SEAL was overloaded: carrier-rep (blind C_E vs imported C_E) must COMPUTE a three-valued residual
  (DO-NOT-FOLD), ⛔ not be a true seal. Reserved SEAL strictly for PIT residues/primes/samples + digests.

### Gate discipline
Decision-list class ⇒ ONE decision-leg pass, fold once, go (G2 / c1 comparator precedent `84686a54`). ⛔ Not
iterated-to-green. Revision 2 is value-free (re-leak-gated) + internally consistent. The SCRIPT (astra-built) then
gets review-until-clear (fresh Opus + Grok); those legs are instructed to SCRUTINIZE the B schema-bridge adjudication.
Commit the reviewed folded directive before the build.
