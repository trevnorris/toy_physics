# Decision-review — S11c-c2 N6 cross-engine comparator BUILD DIRECTIVE

## Artifact (orchestrator-written build directive — review the DESIGN, ⛔ not a script yet)
`research/pde_ledger_v3/directives/S11c_c2_N6_comparator_build_directive.md`

This is a PRE-BUILD decision review of a build directive for a cross-engine measurement instrument. There is no
comparator script yet — ⛔ do NOT ablate a fictional script. Review whether the DIRECTIVE's decisions are sound,
complete, and value-free, so that a builder handed ONLY this directive produces a correct comparator. Substantiate
every finding against the real files below; ⛔ a prose assertion with no file/line evidence is discarded.

## What to check (derive your own view from the sources FIRST, then read the directive)
1. **The N6 CROSS-ENGINE SEAL is the load-bearing decision — is it correct?** The directive claims the two engines'
   PIT residues are NOT cross-differenceable because the engines use different prime sets AND different sample
   points (a deliberate blindness feature), so the cross-engine channels are SYMBOLIC (translate WL `Inactive[]` ↔
   SymPy `_NODES` DAG) + STRUCTURAL (dimensions + nonzero-support via a column↔slot bridge), with PIT residues,
   digests, and SymPy-only families SEALED (surfaced, never differenced). **Verify the premise empirically**:
   - the WL primes/sample points vs the SymPy primes/sample points (WL `mathematica/out/S11c_c2_N6_mathematica_audit.out`
     `N6_LOCAL_PROBE` tag; SymPy `scripts/out/S11c_c2_N6_{covariance,reconcile,diagnostic}_sympy.out` `PRIMES`/`PIT_PROVENANCE`);
   - that the symbolic operands genuinely exist on BOTH sides (WL `ARITHMETIC` field; SymPy `_NODES` DAG);
   - that the digests are incompatible schemes (WL integer fingerprints vs SymPy sha256).
   Is SEALING the PIT the RIGHT call, or is there a legitimate cross-engine numeric channel the directive wrongly
   discards? Is the "column↔slot bridge is a computed residual adjudicated post-run, not a pre-registered collapse"
   framing correct, or does it hide a needed structural comparison?
2. **Is the JOIN SET complete and correct?** Cross-check the directive's join set (reconcile `N6RC_*`, covariance
   `N6COV_*`, guards) against the ACTUAL tags each engine emits. Is any cross-engine-shared family MISSING from the
   join? Is any SymPy-only family (no WL sibling — `MU_RECONSTRUCTION_RESIDUAL`, `CARRIER_RECONSTRUCTION_RESIDUAL`,
   `REP_INVARIANCE_RESIDUAL`, `PREMISES`) wrongly placed in the join, or a genuinely-shared family wrongly excluded?
   Confirm WL's `N6_` namespace really is thin (only guards + `LOCAL_*` + `N6RC_DIMENSIONS`).
3. **Three-valued / no-verdict / no-target / rule-5 discipline.** Does the directive preserve the three-valued
   residual and ban per-case verdict tokens? Does it leak ANY expected value or agreement prior (a residual target,
   a baseline count, a prime value, a "should be zero/covariant")? Quote any leak.
4. **Authority application.** Does the directive correctly point at + apply c2 `SHARED_PHYSICS.md` §N8/T7 (~line
   477), §5c (~line 303), §1b (~line 76)? ⚠ SPEC-CURRENCY: §N8/T7 was written for the self-energy-increment
   comparator; is it valid to apply it to the N6 symbolic/structural test given the PIT blindness — or does the N6
   comparator need a distinct contract clause? Is the "giants deferred ≥64 GB, name-don't-adjudicate" scoping a
   faithful reading of §N8/T7, or an over/under-scope?
5. **Reuse claims + the NEW `_NODES` reader.** Are the reuse claims about the S11c-c1/S11c-b comparator base
   accurate (does `S11c_c1_cross_engine_comparator.py` actually provide `load_wl`, `parse_wl_value`/`canonical_basic`
   as a WL→SymPy translator, the three-valued `residual`, lazy materialize/release)? Is the directive right that a
   NEW SymPy `_NODES` op/args DAG reader is needed (no srepr on the SymPy side)? Sample the real `_NODES` payload +
   its emitter (`emit(name+'_NODES', …)` in the reconcile/covariance scripts) to confirm the DAG schema is
   reconstructable as specified.
6. **Container traps.** Are the per-family extraction traps (DIMENSIONS shape mismatch; `ARITHMETIC`+PIT dual form;
   `FROZEN_PHI` Wolfram-expr vs expr-string; `PRIMES`/`PIT_PROVENANCE` dedicated-vs-embedded) accurate against the
   real payloads? Any trap MISSING that would cause a silent 0-extract or false agreement?
7. **Over/under-specification.** Is anything under-specified such that a blind builder would guess wrong (rule 6 —
   under-specification has cost this ledger more than contamination)? Is anything over-specified into the builder's
   answer? Is the DoD able to FAIL (does it catch a silent 0-join, a manufactured join, a PIT-seal violation)?

## Sources (read these; derive independently)
- The directive under review (above).
- `directives/S11c_c2_SHARED_PHYSICS.md` §N8/T7, §5c, §1b (the governing contract — it WINS over the directive).
- The two engines' outputs: WL `mathematica/out/S11c_c2_N6_mathematica_audit.out` (⚠ 400 MB — targeted grep/awk
  only, NEVER cat); SymPy `scripts/out/S11c_c2_N6_{covariance,reconcile,diagnostic}_sympy.out`.
- The engine sources `scripts/S11c_c2_N6_{covariance,reconcile,diagnostic}_sympy.py` and
  `mathematica/S11c_c2_N6_mathematica_audit.wl` (for the tag/DAG/column schemas).
- The precedent comparator `scripts/S11c_c1_cross_engine_comparator.py` + its directive
  `directives/S11c_c1_comparator_build_directive.md` (the delegated-build form this mirrors).

## Physics filter
Report a finding only if it catches a way the DIRECTIVE would produce a WRONG or false-agreeing comparator (a
missed cross-engine channel, a wrong seal, a manufactured join, a leaked value, an under-specification a blind
builder guesses wrong, a spec-currency error). ⛔ Do not report style. State each finding as: the decision, the
file/line evidence, why it is wrong, and the minimal fix.

## Method
This is a DECISION review (one pass). Read the sources, form your own view of what the N6 cross-engine test CAN
and CANNOT establish given the two engines' actual outputs, THEN read the directive and report where it diverges.
⛔ Do not build or run the comparator. End with an overall SOUND / NOT-SOUND on the directive + the ranked findings.
