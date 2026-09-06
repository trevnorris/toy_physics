# Decision review — the S11c-c2 EXPORT REPAIR build directive (publication-only)

You are one of two independent decision legs on an **orchestrator-written build directive**. This is a
**decision-list review: ONE pass** — check the requested decisions and their supporting evidence against the real
files, report anything that would make the builder build the **wrong** thing or **overstep**. ⛔ Do NOT run a
fictional-script ablation (the target script already exists and its physics is separately reviewed-clear; you are
reviewing the *directive's decisions*, not re-reviewing the physics). ⛔ Do NOT modify the working tree.

## Artifact under review
`research/pde_ledger_v3/directives/S11c_c2_export_repair_directive.md` — a publication-only repair of the export
serialization in `scripts/S11c_c2_selfenergy_fold_sympy_audit.py` (the 60 MB `scripts/S11c_c2_exports.py` must
shrink to carry only what a later step binds, in a transparent-compact form).

## What you are handed (read any)
- the directive above;
- the target script `scripts/S11c_c2_selfenergy_fold_sympy_audit.py` (esp. `EXPORT_ROOTS`@48, the `export_key`
  map in `run()`@~894, `publish`@807-853, `emit`@746);
- the export architecture `directives/export_ledger_bind_closure_design.md` (D1 membership = what a LATER step
  binds; D2 own-rows delta; D3 guard) and the fold machinery `scripts/ledger_fold.py`
  (`check_consumer`, `assert_lookups_equal_manifest`, `assert_delta_is_minimal`);
- the sub-step scopes `directives/S11c_decisions.md` (esp. the S11c-c "exports" column and the S11c-d row: what
  the profile-conditioned spectrum/scattering step will CONSUME);
- the physics authority `directives/S11c_c2_SHARED_PHYSICS.md` (§2 the three named objects, §3c the increment as
  a comparator/export representation, §7 the chain-output/export schema);
- the physics review adjudication `_measurements/S11c_c2_physics_review_adjudication.md` (context: physics is
  sound; the increment is not a downstream binding);
- the real parent exports `scripts/S11c_b_exports.py`, `scripts/S11c_c1_exports.py`.

## Settle these decisions (ground each in the real files, cite file:line)
1. **Scope.** Are the three named edit sites (`EXPORT_ROOTS`, the `export_key` map, `publish`) the **correct and
   complete** set to achieve the repair *without touching construction code*? Is there any way the required change
   forces an edit to physics/emit code (which the directive forbids)? Is the "STOP and report if it crosses into
   construction" fence correctly placed?
2. **R1 membership.** Is exporting **BOTH** `s11cc2ClosedSlabOperator` and `s11cc2ClosedCouplingKernel` (all 4
   cases) and dropping `s11cc2SelfEnergyIncrement` to EMIT-only **correct** per D1 (what S11c-d binds)? Does
   S11c-d need the full closed slab operator (its closure-modified diagonal) as well as the coupling kernel, or
   only one? Is the increment genuinely a non-binding comparator representation (§3c)? Would this membership
   under-export (break d) or over-export?
3. **R2 representation.** Is "transparent-compact (factor/cancel/collect/CSE), ⛔ not expand, ⛔ not opaque, +
   per-(object,case) semantic `expand(decoded − emitted)==0` hard-stop, upgrading the srepr-`==` roundtrip"
   **sound and sufficient**? Will a factored form remain importable + differentiable by a downstream `load_model`
   consumer? Is the semantic zero-test the right check (vs the current srepr equality)? Any way a compact form
   silently loses information the equivalence check wouldn't catch?
4. **R3 hygiene.** Is dropping/shrinking the `display = sp.sstr(value)` duplicate **safe** — does any consumer
   (`ledger_fold.load_model`, `check_consumer`, or any downstream step) bind or require `display`? Confirm against
   `ledger_fold.py` and how the parent exports' rows are consumed. Is there any OTHER non-binding oversized field
   that should also go, or any field the directive risks dropping that a consumer DOES need?
5. **R4 guards.** Are the preserved guards (`check_consumer`, `assert_delta_is_minimal`, F9 collision,
   `BUILD_INPUT_DIGESTS`, `_restore` roundtrip) the right and complete set for a correct minimal delta?
6. **Leak / overstep.** Does the directive leak any expected value (a target size, a specific factored form, a
   pre-registered residual)? Does it correctly bound the builder to build→verify→report and forbid the
   self-review/skill-reading sprawl that the previous build committed? Any instruction that would make the builder
   do the wrong thing?

## Output
For each of 1–6: your finding + evidence (file:line in the real files). Separate MUST-FIX (would build the wrong
thing or overstep) from nits. End with: **is this directive correct to hand to the builder as-is, or the exact
list to fold first?** (One pass — a concrete fold list, not iterate-to-green.)
