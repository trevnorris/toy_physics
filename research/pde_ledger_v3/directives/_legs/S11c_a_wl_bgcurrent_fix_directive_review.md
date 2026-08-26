# Independent review of a build DIRECTIVE (decision list) — S11c-a WL background-current fix

## Artifact under review
`research/pde_ledger_v3/directives/S11c_a_wl_bgcurrent_fix_directive.md` (+ its measurements twin
`research/pde_ledger_v3/directives/_measurements/S11c_a_wl_bgcurrent_fix_directive.md`).

This directive will be handed to a builder (Codex) to modify the Wolfram engine. You are reviewing the
DECISION LIST itself, before any build runs. Your job is to find defects in the directive that would cause a
faithful builder to compute the wrong thing, break something, or leak an answer — NOT to perform the fix.

## Read the source of truth FIRST, form your own view, then read the directive
- Spec: `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md` — §1 (the drain `v_bulk_normal_0`), §1b
  (`j = ρ_4D v_bulk`, continuity), §2d (`𝔅⁰`, `V_s⁰=J_s⁰=0`), §3c (the background current vanishes; "none may
  be introduced as a free premise").
- The WL engine to be modified: `research/pde_ledger_v3/mathematica/S11c_a_interface_geometry_mathematica_audit.wl`
  — the background current defs (438-446), `bulkVelocityZero` (115-119), `rhoBulkZero` (434-435),
  `traceFieldInventory` (612-618), the projection assembly (`projectionTermsSource`, ~798-822), and every
  consumer of `traceFieldInventory` (1106, 1869, 2082).
- The sibling SymPy engine (for the property being mirrored): `…/scripts/S11c_a_interface_geometry_sympy_audit.py`
  lines 659-664.

## What to check — report any defect you find
1. **Scope correctness (both directions).** Does the directive's change site + consumer list cover EVERY object
   that consumes the background current, or does it miss one? Conversely, does it risk changing something it
   should NOT — the perturbation current, the (legitimately nonzero) background density `rhoBulkZero`, the
   background velocity construction, or the pure face-geometry objects? Is there any object in the engine that
   legitimately NEEDS a nonzero background current (i.e. would the fix wrongly zero a real term somewhere)?
   Derive from the spec whether `j⁰=0` is genuinely scope-wide (T-c…T-i) or whether some object is out of that
   scope.
2. **Leaked answer / rule 5.** Does the directive leak an expected OUTPUT value, a cross-engine target, or a
   "fix until it matches" acceptance? The construction `j⁰ = ρ_4D⁰·v_bulk⁰` with the engine's existing zero
   background velocity is a supplied premise — is that legitimately a SUPPLIED premise, or does it function as
   a leaked answer? State your reasoning.
3. **Property vs recipe / exception.** Is the load-bearing requirement (P1) stated as a PROPERTY that must hold
   ("no free-premise background current; it is `ρ⁰·v⁰` reached by construction"), or as a brittle recipe / a
   named exception that could breed a regression?
4. **Acceptance is value-free AND able-to-fail.** Can the stated acceptance actually FAIL if the builder does
   the wrong thing? Is it verifiable by provenance rather than by grepping for deleted symbol names (the
   directive forbids self-certifying by name — is that enough)? Is there a way a builder could satisfy the
   letter of the acceptance while leaving the defect (e.g. a substitution rule zeroing the old symbols post-hoc
   — the directive forbids this; is the forbiddance sufficient)?
5. **No reliance on the continuity-cancellation.** The directive says the surviving-term-reduces-to-zero result
   (from the physics consults) must NOT be a premise of the build. Confirm the directive does not smuggle it in
   as something the builder must reproduce or rely on.
6. **The three script clauses (rule 2).** Are the print-not-assert / operands+residual / interpretation-in-record
   obligations present and correct for this build?

## Method
- This is a DOCUMENT review. Quote both the spec/engine source and the directive for every finding.
- Physics filter: report a defect only if it would cause wrong physics, a broken build, a leaked answer, or a
  missed/over-broad scope. Do not report style.
- You do NOT need to run the engine. If you make a claim about what the engine computes, back it by reading the
  raw source at named lines (or a small CAS check with literal stdout at a named /tmp path).

## What to return
A list of defects (or "no defect on axis N" per axis above), each with: the directive quote, the spec/engine
quote that decides it, and the concrete failure a builder would produce. If you believe the directive is
correct and complete on an axis, say so explicitly — a leg that finds nothing is weak evidence, so state what
you checked.
