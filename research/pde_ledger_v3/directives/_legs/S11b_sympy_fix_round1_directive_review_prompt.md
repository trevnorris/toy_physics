# Directive review — S11b SymPy engine fix-round-1

You are one of two independent legs (the other is a different engine) reviewing an **orchestrator-written
repair directive** BEFORE any builder runs it. This is a document review: does the directive specify the
right repair, completely, without leaking a value, without over-specifying, without breeding new defects,
and without disturbing physics that is already correct?

## Artifact (the thing under review)
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11b_sympy_engine_fix_round1_directive.md`
and its rule-2 twin
`/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S11b_sympy_engine_fix_round1_directive.md`.

## What you are handed (verify against these, do not take the directive on faith)
- The round-1 engine to be repaired: `research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py`
  (baseline `7dd89076`) and its product `research/pde_ledger_v3/scripts/S11b_exports.py`.
- The physics authority `research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md` (esp. §6 energy
  accounting / convention checks, and B7 dimensions) — the directive must NOT contradict it.
- The wiring authority `research/pde_ledger_v3/directives/S11b_sympy_build_directive.md`, the `F9`/`D`
  contracts (`S11_export_chain_decisions_v2.md`, `S9_REWRITE_PLAN.md`), and `CLAUDE.md`
  (rules 2, 3, 5, 11, 12, 15) + `.claude/skills/build/SKILL.md` (the three script clauses; the anti-tautology
  corollaries 1 & 3; the "hardening a check adds no-physics lines" measurement).

## The two defects the directive repairs (verify each is real, against the engine)
- **Fix 1 (F4):** the load-bearing impedance is a typed literal `Z_GENERAL` (line 138); the bulk solve that
  computes the same object (`z_plus`/`z_minus`, ~612–620, emitted as `Z_IMPERMEABLE`) is discarded; the
  assembly consumes the typed object (default `z_value=Z_GENERAL`, ~478). Read these lines and confirm.
- **Fix 2:** ~12 EXPORTED ledger rows are tautological self-checks (residual 0/−1 by construction). Read the
  flagged lines in the directive's table and confirm each is built from fresh placeholders / a self-identity
  and does not read the assembled `MODEL`; confirm the sibling `convention_check_conservative` (~943–953)
  IS genuine (the pattern the directive points to).

## What to check — report a finding only if it would change the repair or its outcome
1. **Accuracy.** Is any defect mis-described, any line reference wrong, or any flagged row actually genuine
   (not tautological)? Quote the engine line and the directive line.
2. **Completeness.** Is any tautological / typed-payload defect of the SAME class present in the engine but
   MISSING from the directive's list? (Scan for other self-checks that would be byte-identical under a form
   ablation of what they claim to police, or other typed load-bearing objects.) Name it with a line number.
3. **Leak (rule 5).** Does the directive state, or let a builder infer, any expected physics VALUE, sign,
   ratio, or count? The repair must be specifiable as an emission PROPERTY only.
4. **Over-specification (rule 3).** Does the directive dictate a code recipe where it should name the object
   and its required property? (Fix 1 should require "computed from the acoustic field; ablation moves the
   impedance-dependent tags," not a line-by-line rewrite.)
5. **Breeding risk (rule 15 / step-run).** Could a builder satisfy this directive by ADDING no-physics lines
   or new tautologies? Does the directive's "make genuine, or make §10-local/delete; ⛔ do not add checks"
   framing actually prevent that, and is the genuine-vs-local test (moves under a FORM ablation) decisive?
6. **Scope preservation.** Are the out-of-scope exclusions correct and safe (energy-basis count 11-vs-10;
   the B5 longitudinal fate)? Does the directive correctly forbid changing any physics VALUE and the export
   wiring (F9/bindings/digests/D3/_RELATIONALS/freeze/carry-forward)? Could Fix 1 (recomputing the impedance)
   inadvertently change an exported VALUE, and does the directive guard against that?
7. **Spec fidelity.** Does making the Fix-2 checks "genuine against the assembled MODEL" match what §6/B7
   actually mandate, or would it manufacture a check §6 does not ask for?

## Method
Read the sources first, form your own view of each defect, THEN read the directive and judge it. For any
claim of fact about the engine, quote the exact line. ⛔ A prose assertion with no engine citation is
discarded. If you believe a flagged row is genuine (not tautological), show the code path that makes it
able to fail. If you claim a defect is missing, give its tag and line.

## Output
Findings most-severe first, each with: the directive location, the engine/spec citation, and the concrete
fix. If the directive is correct on a given axis, say so per axis (1–7). End with a one-line verdict: is the
directive safe to build from as-is, or must something be fixed first?
