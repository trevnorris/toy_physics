# Decision-leg review: S11c-c2 N6 reconcile COLLAPSE instrument build directive

## Your role
You are one of two INDEPENDENT decision legs gating a builder (astra) that will write a CAS collapse-test
instrument. No builder launches until both legs report. Review the DIRECTIVE + its physics-bearing content (the
frozen bridge dictionary and the collapse-witness design). ⛔ The instrument does not exist yet — do not ablate a
nonexistent script; review the DESIGN by reading, and derive the maps against the real engines.

## Artifact under review
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_reconcile_collapse_build_directive.md`

## Read these (form your own view from the engines first)
- The framing this instrument implements (already vetted + folded):
  `research/pde_ledger_v3/_measurements/S11c_c2_N6_reconcile_question.md` and its vet record
  `…/S11c_c2_N6_reconcile_question_vet.md`.
- The two engines (the maps must be grounded in BOTH):
  `scripts/S11c_c2_N6_covariance_sympy.py`, `scripts/S11c_c2_N6_reconcile_sympy.py`,
  `mathematica/S11c_c2_N6_mathematica_audit.wl`, and the shared module `scripts/S11c_c2_N6_diagnostic_sympy.py`.
- The comparator whose extraction is reused: `scripts/S11c_c2_N6_cross_engine_comparator.py`.
- The c1 reconcile precedent (staged representational bridge): `_measurements/S11c_c1_comparator_reconcile.md`.
- The per-engine resolution to stay consistent with: `_measurements/S11c_c2_N6_RESOLVED.md`.

## What to check — report a finding only if it changes what the instrument would COMPUTE or CLAIM
1. **The frozen bridge dictionary — each of the 8 admitted maps.** For each: is it PRIMITIVE (a symbol/structure
   convention, a jet identity, a profile substitution), or is it a disguised WHOLE-OBJECT equality that would
   make the collapse VACUOUS (identifying the two engines' carriers/μ/sources/Φ outright)? Is it independently
   JUSTIFIED and GROUNDED in BOTH engines at the cited lines (verify the citations)? A map that is really the
   §3 question in disguise is a MUST finding (rule 5 / blanket-collapse).
2. **Completeness vs over-reach.** Is any NEEDED primitive MISSING (a convention the surfaced operands plainly
   require to be comparable, e.g. a units/normalization the engines emit differently)? Is any admitted map
   unnecessary or over-broad? Is the EXCLUDED list right (whole-object equality; LAB_HELD↔MATERIAL_ADVECTED;
   σ_W binding/σ_W→0; default on-shell/Fourier)?
3. **Retained order / no proxy.** Does the design force coefficientwise-per-grade `(η^i σ_W^j), i,j∈{0,1}` with
   no grade-combining, no σ_W binding, no σ_W→0? Confirm both engines retain that independent rectangle (cite).
4. **Collapse witness.** Is it genuinely three-valued (collapsed / residual_remains / undecided), does it PRINT
   the residual before the witness, is the fresh PIT independent of both engines' primes/points, is it bounded?
   ⛔ Is there any path where it asserts or emits a verdict?
5. **The 4 mandatory controls.** By design, will each BITE — blanket-collapse forces all-collapsed; one-sided
   dictionary corruption moves only the leaves using that map and drives ≥1 residual nonzero; dropped-primitive
   leaves a residual; grade-combining is impossible/flagged? Is any control tautological or unable to fail?
6. **BASELINE + census.** Is `SOURCE_BASELINE` correctly a nominal-control duplicate (not a 3rd discriminator)?
   Are the 8 census leaves audited as production-control + domain-coverage equivalence (not dismissed as
   metadata, not pre-decided)?
7. **Value-free / leak.** Does the directive leak an expected collapse outcome, count, or pass condition
   anywhere? (The prior-run nonzero counts are facts of the committed comparator run, not targets — is that line
   honored?)
8. **Fence + DoD + extraction reuse.** Is the astra fence (build→verify→run→report→STOP; no self-review; no
   outside edits) explicit? Is reusing the comparator's extraction (not re-implementing the join/schema bridge)
   correct and sufficient? Is the DoD able to detect a non-built or blanket-collapsing instrument?

## Method / evidence discipline
- Read the engines first; derive each candidate map yourself and confirm it is primitive + grounded. ⛔ A prose
  assertion is discarded — cite file+line for every claim about what a map is or what an engine computes.
- ⛔ You are not adjudicating whether the residuals collapse (that is the instrument's job). You are checking the
  directive would make the instrument ask the right question with primitive maps, at retained order, deciding
  nothing.
- End with a one-line verdict: **DIRECTIVE SOUND** (clear to build) or **DIRECTIVE NOT-SOUND**, listing the
  MUST-level findings first.
