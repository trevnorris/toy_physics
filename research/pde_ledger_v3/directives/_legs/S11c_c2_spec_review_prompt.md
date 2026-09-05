# Independent physics review — S11c-c2 SHARED_PHYSICS (the self-energy fold), DECISION-LIST gate

You are one of two independent legs reviewing an **orchestrator-written physics spec** *before any builder is
launched* (CLAUDE.md rule 7). This spec settles the structure, physics, provenance, and scope that the S11c-c2 build
directive and its two blind engines (SymPy + blind Wolfram) will honour. A defect here propagates into both engines
and the comparator, so it is the one artifact the builder trusts. Your job is to find every way the **decisions are
wrong or under-specified** — not to rubber-stamp.

## Artifact under review
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md`

## What S11c-c2 is
c1 solved the curved two-face outgoing bulk problem and exported the **closed permeable face response**
`(δp_s,J_s,t_s)(V_s,μ_θ)`, the nonlocal DtN operator, its flat symbol, and its two-momentum kernel. c2 **folds** that
closed response into the S11c-b variable-coefficient slab operator `S11CB_SLAB_OPERATOR` (whose θ-row and mechanical
rows still carry `δp_s` and the response kernels **symbolically**) and **re-extracts** the off-diagonal
transverse↔`{θ,e_W,u_L}` coupling from the **closed** full operator, yielding the coupled **nonlocal self-energy
operator**. It is a CODE build; this is the pre-build decision gate.

## What you are handed (read the SOURCES OF TRUTH first, form your own view, THEN read the spec)
Read in this order (a method instruction, not a blindness control):
1. The parent decision list `directives/S11c_decisions.md` (the N-series; the ratified S11c family scope).
2. The c1 step record `steps/S11c_c1_curved_bulk_closure.md` (the reviewed close: what is cross-engine AGREE vs
   UNDECIDED vs DEFERRED). The reconcile record `directives/_measurements/S11c_c1_comparator_reconcile.md`.
3. The S11c-b step record `steps/S11c_b_variable_coefficient_operator.md` (the slab operator; the DEFERRED
   cross-engine residual; the two whole-row sign conventions + #90's two flags that are cross-engine-UNVALIDATED).
4. The c1 spec `directives/S11c_c1_SHARED_PHYSICS.md` §1d (the Λ-channel placement; the line "the routing of `J_s`
   into the mass row and `t_s` into the mechanical rows is S11c-c2's, not c1's") and §3b (the response, the operator
   inverse, and the energy audit note that the traction-vs-slab pairing "is S11c-c2's, after the fold").
5. The S11c-b spec `directives/S11c_b_SHARED_PHYSICS.md` §3b (the divergence-form slab operator; the θ-row =
   sourced mass-evolution restored after variation) and §3c (the off-diagonal coupling by weak variational
   restriction — "extract both weak blocks from the §3b operator itself, not a parallel route").
6. The two real export files `scripts/S11c_b_exports.py` and `scripts/S11c_c1_exports.py`, and the fold module
   `scripts/ledger_fold.py` (`load_model(base, *deltas)` — line 102).
7. **Then** the spec under review.

## Required method (DOCUMENT branch + computation where a claim is computational)
Read the sources, form your own independent view of what the fold requires, and only then read the spec. Quote both
sides for every finding.

⛔ **A prose derivation is worth nothing.** Where a claim is computational, **write your own SymPy script BEFORE
relying on the spec's version, and save both the script and its literal stdout to named absolute paths** (under
`/var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/` or `/tmp`); without them your derivation claims
are discarded. Do NOT spawn Mathematica kernels — this is a spec review; the fold checks below are light symbolic
SymPy. Copy anything you run to `/tmp` and never modify the working tree.

## Physics filter
Report a finding only if it catches a way the c2 **decisions** would make the physics wrong, under-specified, or
un-reviewable — not "the spec could be worded better." Rank findings most-load-bearing first, each marked
MUST-FIX / SHOULD-FIX / NIT, with the exact spec line and the source line it contradicts.

## Contested questions you must SETTLE (with computation where marked ⚙)
1. **⚙ The non-commutation is the physics.** §2 claims `extract(close(SLAB)) ≠ close(extract(SLAB))` and that the
   self-energy IS that difference, fixed by the counterexample `R_x=x+p, p=αy`. Verify with your own minimal SymPy
   example that a weak-restriction "extract" and a substitution "close" genuinely fail to commute for an operator of
   the *shape* here (a diagonal transverse block + an off-diagonal block, closure threading a term from one sector
   into another). Is the spec's claim that the correct object is `extract(close(·))` (close FIRST) sound, or is there
   a case where extract-first is correct?
2. **⚙ The substitution-increment cancellation.** §3c/§7 claim the self-energy emitted as
   `extract(close(SLAB)) − extract(SLAB)` makes c2's cross-engine residual **independent of S11c-b's deferred,
   un-validated slab/coupling residual and its two sign conventions**, because those terms appear identically in
   `close(·)` and `open(·)` within each engine and cancel before the cross-engine diff. Verify this: (a) is the weak
   restriction linear, so `extract(close)−extract(open)=extract(close−open)`? (b) does the S11c-b base truly cancel,
   or does the closure MULTIPLY S11c-b terms (so they do NOT cancel and the deferred residual leaks into c2)? This is
   the spec's central design claim — stress it hard.
3. **The Λ-channel routing.** §1d/§3a route the closed `J_s` into the mass/θ row and the closed `t_s` (carrying
   `Λ_X`) into the mechanical rows. Is this the correct reading of c1 §1d + the S11c-b θ-row
   (`evolution_mass_balance − Σ closure_shape_deriv`, #90) + `closure_shape_deriv`'s content? Is `δp_s`'s entry into
   the mass row via `𝒜_s=μ_s−δp_s/ρ_m` handled correctly?
4. **Completeness of the three re-adjudications (§3d).** The spec carries exactly three c1-UNDECIDED items into c2 as
   mandatory re-adjudications: background density (field-vs-field, rule 17), `t_s` (4-vector vs scalar), and the DtN
   whole-form. Are these the COMPLETE set of c1-UNDECIDED/DEFERRED items the fold makes load-bearing? In particular:
   (a) c1 §3b lines 328–330 say the traction-vs-SLAB-row **energy/dissipation pairing** is "S11c-c2's, after the
   fold" — is the ABSENCE of any c2 energy/dissipation object a **scope defect**, or correctly deferred to S11c-e?
   Ground the answer in the parent decision list + the c1/S11c-b step records, not opinion. (b) Does the off-diagonal
   **flat-resolvent leg-labeling** (c1 UNDECIDED) enter c2's fold? (c) Does c2 silently inherit any c1 DEFERRED giant?
5. **Export topology.** §7 sets c2's base as `load_model(scripts/S11c_b_exports.py, scripts/S11c_c1_exports.py)` and
   the consume-set from BOTH parents. Verify against `ledger_fold.py` that a base+one-delta fold binds the c1 delta
   correctly; check for an F9 false-equal / key-collision hazard between the two parents' keys; confirm the guard
   (`check_consumer`/`assert_lookups_equal_manifest`/`assert_delta_is_minimal`) works for a two-parent base. Is
   leaving the exact `IMPORT_KEYS` root set to the build directive (rather than freezing it here) correct, or a gap?
6. **Rule 5 leak check.** Does the spec state, anywhere, an **expected value / sign / coefficient / count / form** of
   a c2 OUTPUT (as opposed to a SUPPLIED inherited fact)? The controls say "must be nonzero", "must move", "must
   vanish in the uniform limit", "must reduce in the zero-DtN limit" — are any of these leaking a computed answer, or
   are they legitimate structural properties (decoupling `N6`, non-commutation `§2`, limit reductions)?
7. **Rule 17 integrity.** Is the §3d.1/§5d density re-adjudication a genuine field-vs-field test, or does it smuggle
   the vacuous `∇ρ→0` / `∇W→0` uniform limit (forbidden by `N6`) as if it were a corruption?
8. **Rule 3 / rule 2 corollary 3.** Any decision that specifies a RECIPE where it should name an OBJECT (manufacturing
   a derivation-path question)? Any residual that is tautological (zero by construction for any input, e.g. the §3b
   adjointness residual or the §3c increment presented as a "check")? The spec claims the §3c increment is an export
   representation, NOT a check — is that framing honoured, or does it read as a tautological residual somewhere?

## Output
For each of the 8 questions: your verdict + the computation/quote that settles it. Then the ranked findings list.
End with the single most important thing that must change before a builder is launched, or state that the decisions
are sound to build against after any folds you list. Two legs finding real defects is the expected, productive
outcome — do not soften to reach agreement.
