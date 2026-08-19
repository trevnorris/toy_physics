# Independent review — the S11b unified decision list

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11b_unified_decisions.md`
Its rule-2 substantiation (the command + literal output behind each factual claim):
`/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S11b_unified_decisions.md`

## What this is
An orchestrator-written DECISION LIST (the "G-series") for rewriting step S11b ("the interface coupling
law") as ONE unified step on the project's export-chain scheme, subsuming two historical execution
stages: S11b-A (the bulk's response to moving faces) and S11b-B (the homogeneous interface assembly). It
is NOT the shared spec and NOT a build directive — it settles the structure, naming, provenance, and
scope decisions the spec and build must honour. Your job: find where a decision is **wrong in a way that
would corrupt the build**, because both engines will read the resulting spec, so a wrong shared decision
makes both engines agree on the same wrong thing.

## Sources of truth — read these FIRST, form your own view, THEN read the decision list
- S11b-A spec + record: `research/pde_ledger_v3/directives/S11bA_SHARED_PHYSICS.md`,
  `research/pde_ledger_v3/steps/S11bA_interface_response.md`
- S11b-B spec + record: `research/pde_ledger_v3/directives/S11bB_SHARED_PHYSICS.md`,
  `research/pde_ledger_v3/steps/S11bB_interface_assembly.md`
- The handoff: `research/pde_ledger_v3/steps/S11b_HANDOFF.md`
- What S11 defers to S11b: `research/pde_ledger_v3/steps/S11_stray_longitudinal.md` (ownership boundary
  near the end)
- The export-chain rules the list inherits: `research/pde_ledger_v3/directives/S11_export_chain_decisions_v2.md`
  (F1–F9); `research/pde_ledger_v3/directives/S9_export_chain_rebuild_directive.md` (the blind-Wolfram control)
- The frozen comparator contract: `research/pde_ledger_v3/directives/S11_C17_C18_spec_repair_decisions_v2.md` (T7)

## Method (DOCUMENT review)
Read the sources first and form an independent view of what S11b-A and S11b-B actually established and how
their objects relate. THEN read the decision list. For every finding, quote BOTH sides (source and
decision list) with `file:line`. A prose claim with no `file:line` citation is discarded.

⛔ Do NOT trust the decision list's own citations — open each cited line yourself and confirm it says what
the decision claims. (Two citation errors were already found and fixed during authoring; assume more may
remain.)

## The specific questions to settle — independently (read/derive; do not take the list's word)
1. **G2/G4/G5:** Is the "B-led, A-as-reduction" framing correct? Is B's `Λ_A⁰` genuinely a DIFFERENT
   coefficient from A's `Λ_p⁰` (affinity vs raw pressure), related only on a `μ_s=0` slice — or are they
   the same object (which would make "keep distinct" wrong)? Is A's single `τ` genuinely a restricted
   slice of B's three independent times `τ_A, τ_V, τ_X`? Cite the spec lines.
2. **G3:** Are S11's `v₀` and S11b's `v₀` genuinely DIFFERENT physical quantities (in-plane brane vs
   normal bulk drain), or the same quantity S11 merely pins to zero (which would make "name apart, not a
   premise override" wrong)? Read both specs and decide. Also test the claimed F9 hazard: would an
   object comparison of two bare `Symbol('v_0')`s actually score EQUAL and silently merge them? If you
   can, verify with a tiny SymPy script and report its literal output.
3. **G6:** Is `Λ_X⁰` genuinely a SUPPLIED/prescribed constitutive input in B (so the unified step derives
   its consequences, not its existence) — or does B actually DERIVE the channel's existence (which would
   change the provenance decision)?
4. **G14:** Are the items listed as "deferred to C" genuinely out of scope for a UNIFORM-background step,
   or does any of them actually have to be answered here?
5. Any place the list **restates** an inherited obligation (F1–F9, the three script clauses) in weaker
   words instead of pointing at it — a known failure mode. Flag it.
6. Any place a decision **presupposes the FORM of an answer** the engines are meant to compute (a fixed
   power law, sign, or functional form) — flag it; it leaks the answer to the builder.
7. Any leaked acceptance VALUE the builder could iterate toward (the list is supposed to WITHHOLD the
   `μ_s=0` reduction's value and state it only as an obligation-to-compute).

## Physics filter
Report a finding only if it catches a way a DECISION is wrong — a misnamed object, a mis-assigned
provenance (import / derive / supply / reduction-slice), a wrong scope boundary, a presupposed answer
form, a leaked value, or a restated-not-pointed obligation. ⛔ Do NOT report style, phrasing, or "you
could add more."

## Hygiene
Read-only. ⛔ Do not modify the working tree. If you run any check (e.g. a script to test the F9
symbol-equality claim), save the script AND its literal stdout to a named absolute path and report that
path. A prose "I checked and it's fine" is discarded.

## Output
A list of findings, each: the G-item, the defect, both-sided `file:line` quotes, and why it would corrupt
the build. If nothing survives the filter for an item, say so plainly. End with the single most important
thing you would change before this list is folded.
