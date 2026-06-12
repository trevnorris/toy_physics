# Step 7 / branch-realization — parent-status freeze decision (Claude+Codex consult)

**Date:** 2026-06-12
**Shape:** read-only Claude+Codex consult (`codex-chat -s read-only`, session `019ebd26-2291-7d50-bbba-8523bd693bb2`; transcript `research/pde_ledger/redteam_adversarial/codex_logs/step7_parent_freeze_consult.txt`, prompt `..._consult_prompt.txt`). Both engines read the governing docs + ledger independently and **converged on all five points**. This is a MATH/METHODOLOGY decision (how to formalize the parent action so a solve can run), delegated to Claude+Codex by the user 2026-06-12; see `[[claude-codex-resolve-math]]`. No conceptual escalation is triggered (the one trip-wire — point 4 — came back clear).

## The question

At the Step-7 handoff to layer-3 (Stage-1 numerical branch realization), the branch-realization plan's §8 prereq #1 ("parent operator frozen") is NOT met: the strict parent action is `S_ψ[ψ,A;Σ] + S_EM[A]`, with `Σ/R` only a confinement-coupling argument, not an autonomous dynamical throat field; the quadratic wall closure `S_eta^(2)` is "effective unless promoted." `ACTUAL_BRANCH_PROTOCOL_V1.md` frames the fork: **declare strict `S_Σ[R]` promotion (Path A) OR explicitly label the wall an effective closure (Path B).**

## Findings (both engines converged)

1. **Parent operator is genuinely unfrozen** — no later freeze supersedes the V2 status correction. `notes/moving_throat_pde_program_compact.md:196,217`; executable guard `physical_nonlinear_model.py:126,135` returns `physical_export_permitted=False`; `NONLINEAR_PROTOCOL_V2.md:14`.

2. **Path B is a scientifically valid falsification test of the declared EFFECTIVE reduced branch — NOT circular.** Stage 001 introduces `S_eta^(2)` as an ansatz with *fixed* constitutive functions, not coefficients refit stage-by-stage (`stage001_geometry_lift.md:194,214`); V2-16 forbids residuals feeding back into branch data (`stage_v2_16_branch_freeze_no_refit_derivation.md:236,247`); and decisively, the prior V2-23 reduced branch ran **target-blind, came out open/stable, and FAILED the targets rather than being rescued** (`stage_v2_23_minimal_branch_solver_derivation.md:24,320`). A construction that bakes in the answer cannot fail target-blind — so B tests something real. Its value: the V2-23 negative was *noise-dominated*; a properly §5-validated B-run converts that untrustworthy negative into a **trustworthy verdict** (exactly the brief §9 deliverable).

3. **Path A is NOT derivable from the existing ledger.** The wall constitutive functions (`mu_eta`, `T_w`, `T_Omega`, `K_eta`) are red-team-classified `free_choice`/posited (`provenance/_synthesis/batch_01/fit_stage001_...mu_eta.yaml:11`, `...k_eta.yaml:19`). A promoted nonlinear `S_Σ[R]` would be a **fresh parent-level posit** unless independently derived — reopening a full layer-2-style audit cycle for the throat sector (action, variation, boundary terms, coupled residuals, Jacobian/manufactured tests, stability gates, no-refit provenance). Large scope; no compute can start until it's done.

4. **Conceptual flag — CLEAR (no user gate).** Path B, scoped to Stage 1, does **not** foreclose any physical claim. Stage 1 (2D, gravity-side coefficients) captures the P2 response but not the spin-½ chain anyway (`branch_realization_execution_plan.md:39`); the half-integer quantizer is conditional on an autonomous self-reproducing GNLS soliton, a step `lepton_work.md:1424,1444` records as "not yet file-grounded." The autonomous-field question becomes a SEPARATE later pre-registered test (= Path A as its own program). **The only thing that WOULD make this conceptual is declaring effective closure as a *permanent, program-wide abandonment* of ever completing the autonomous-field/eigenmode lane — which we are NOT doing.**

5. **Non-rescue.** B-first is clean only if the pre-registration says *in advance* "this is the effective-closure branch." A B miss falsifies that branch. It may NOT then be followed by "promote `S_Σ` and call it the same Stage-1 attempt." Path A later is allowed only as a new, separately pre-registered model with its own audit trail.

## Decision

**Path B first, correctly scoped.** Run the Stage-1 test as the **"effective-closure branch-realization test"** — NOT a permanent program-wide commitment. Path A (`S_Σ[R]` promotion) is preserved as a distinct future parent-completion program with its own derivation + audit + pre-registration.

Critical scoping correction vs. the first framing put to the user: the word **"permanent"** is what would have made B conceptual. The decision is B **scoped to this Stage-1 pre-registered test**, leaving A open later. With that scoping it is purely methodology — Claude+Codex's call.

`parent_action_status` to be frozen as `effective_closure` in the Stage-1 pre-registration.

## Next steps (flow from this decision)

- **Produce the Stage-1 pre-registration record** (prereq #3, currently absent) — freeze, per `ACTUAL_BRANCH_PROTOCOL_V1.md` + brief §3.4: `parent_action_status=effective_closure`, the branch, support-placement coordinate, D/N / open-impedance boundary + outgoing-port convention, wall constitutive packet, material inputs, port list, extraction formulas, and the **primary** observable. Immutable once frozen (non-rescue).
- Then the solver-build sequencing in `branch_realization_execution_plan.md` §7 (GPE benchmark → manufactured-solution tests → coarse Newton solve → convergence study → PML → conservation → error budget → WP3 tangent → verdict).
- **Compute-spend gate:** brief §3.4 — no compute until the pre-registration is frozen; spending GPU money is the user's call.
