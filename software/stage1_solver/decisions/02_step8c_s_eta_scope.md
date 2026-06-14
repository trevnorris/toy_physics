# Step 8c — `S_η^(ψ,A)` reachability + scope decision (Claude+Codex consult)

**Date:** 2026-06-14
**Shape:** read-only Claude+Codex consult (`codex exec --sandbox read-only`, model gpt-5.5 / reasoning xhigh).
Prompt: `software/stage1_solver/_scratch/step8c_s_eta_consult_prompt.md`; transcript:
`software/stage1_solver/_scratch/step8c_s_eta_consult_codex.log`. Both engines read the governing docs + ledger
independently. This is a MATH/METHODOLOGY determination delegated to Claude+Codex (see `[[claude-codex-resolve-math]]`);
the user pre-confirmed the route ("proceed surrogate + defer physical", 2026-06-14) and the **conceptual flag came back
CLEAR**, so no further user gate is triggered for the scope itself. (GATE A free_choice freeze + GATE B export-guard
flip remain separate USER gates, untouched here.)

## The question

§7 step 8 (WP3 P2 tangent) was decomposed 8a/8b/8c. 8c = current-carrying conservation (deferred from step 6) +
low-frequency response coefficients. The crux predicted at decomposition time: the **matter/gauge→wall source**
`S_η^(ψ)+S_η^(A)` on the RHS of the §4.7 linearized wall equation
(`μ_η∂_t²η − ∂_w(T_w∂_wη) − T_Ω Δ_{S²}η + K_η η = S_η^(ψ)+S_η^(A)+f_ext`). 8a/8b implemented the LHS wall operator and the
wall→matter direction (`δV_conf`, §3.3); the return direction `S_η^(ψ,A)` was deferred (driven by a target-blind
surrogate `f_ext` instead). Is `S_η^(ψ,A)` derivable in the frozen effective-closure scope, and what does its status
block?

## Findings (both engines converged; corrections noted)

1. **`S_η^(ψ,A)` status — OPEN (confirmed).** No complete, certified kernel exists. The only explicit matter form is
   schematic: `S_η^(ψ) ~ −(V_wall'/ℓc)·δρ + ...` (`research/pde_ledger/notes/stages/moving_throat_pde_stage001_geometry_lift.md:349`).
   The gauge source `S_η^(A)` is **named, not formulated** (`notes/moving_throat_pde_program_compact.md:1441`). §4.6
   (`...compact.md:1377-1381`) lists "the full coupled matter/gauge renormalization of these reduced lanes" + "the
   outgoing odd response and quadrupole-normalization branch" as **still open**. Completing it = fresh derivation work.

2. **Physical WP3 d-ln targets are BLOCKED (confirmed).** `d ln R_tr = d ln R_target = d ln ε_η = 0` (compact §13.4)
   are the induced-P₂₂ self-consistent response; the plan states the low-frequency coefficients are determined by
   feedback closure, not static readout (`branch_realization_execution_plan.md:147,151`). The §F extraction formulas are
   algebraic *after the solve* (`stage1_preregistration.md:81`), NOT a substitute for the physical tangent source. A
   surrogate-`f_ext` response is **methodology smoke, not a physical WP3 pass/fail**.

3. **Reachability split — CORRECTED.** My first cut over-protected the verdict by calling the whole "WP1 card" safe.
   The accurate split:
   - **Truly S_η-independent (ω⁰ / source-map, WP1-extractable):** `R_norm` (PRIMARY), `chi_Q`, `N_Q = chi_Q⁻¹`
     (`stage1_preregistration.md:141,203`; `compact.md:6440`). **These carry the Stage-1 falsification load — intact.**
   - **Low-frequency-response-gated (NOT safe-WP1):** `R_pole`, `P_2`, `P_4` depend on `D2/D4/N2/N4`, the ω²/ω⁴
     response-bundle coefficients (`compact.md:2677,2702,2710`). Their physical extraction needs the field→coefficient
     extraction map (TODO) on the realized branch's response operator; whether the bundle is fully determined within
     effective-closure scope or also needs the S_η return coupling is itself part of the open "full coupled matter/gauge
     renormalization" — so they are at minimum extraction-map-gated and potentially S_η-gated. **Not extracted in 8c.**
   - **Induced self-consistent (explicitly S_η-blocked):** `d ln R_tr / R_target / ε_η`. **Deferred to Path A.**

4. **Deferral is HONEST (confirmed).** The leading matter conjugate-variation term is structurally derivable, but the
   full return kernel is not frozen; including only the leading term would dress a partial source as complete. Treating
   full `S_η^(ψ,A)` closure as Path-A-class is correct — Codex's nuance: it **may be a smaller subprogram** than the full
   nonlinear `S_Σ[R]` promotion, but it still reopens derivation+audit and is out of effective-closure Stage-1 scope.
   (`stage1_preregistration.md:15`; `NONLINEAR_PROTOCOL_V2.md:14`.)

5. **Conceptual flag — CLEAR.** Reporting the physical WP3 d-ln targets as deferred/blocked (not pass/fail) does NOT
   abandon a physical claim beyond the already-frozen effective_closure scope (Path B Stage-1-scoped, Path A preserved:
   `branch_realization_parent_status_decision.md:18,24`). The primary `R_norm` test still carries the load.

6. **8c surrogate scope — SOUND, with two guardrails (Codex corrections, folded into the directive):**
   - **Current-carrying conservation:** appropriate because the 8b drive is complex/phase-carrying
     (`p2_driven_absorber.py:236`), BUT it must be a **forced/CAP balance** (`divergence = injected − absorbed`), not a
     closed `div J = 0` (step-6 already models source/sink balance, `conservation_diagnostics.py:463`). **Do NOT reuse
     `_matter_number_current` on the complex phasor lanes — its formula is for REAL branch fields
     (`coupled_branch.py:210`).** Build a **linearized phasor current around the WP1 background** and report source/CAP
     residuals.
   - **Low-frequency surrogate fit:** sound iff the fitted functionals are **target-blind, predeclared / fixed before
     runs, independent of the forcing overlap**, and labelled CAP/operator-methodology coefficients only (8b fitted
     nothing: `step8b_driven_absorber.md:148`).

## Decision

**Proceed with 8c as surrogate-only** (user-confirmed): linearized current/charge/energy **balance** with explicit
drive/CAP source accounting (phasor current around WP1), + low-frequency fits of **predeclared target-blind** response
functionals (methodology validation). **Defer the physical WP3 d-ln targets to Path A** (open `S_η^(ψ,A)`); record that
`R_pole/P_2/P_4` are low-freq-response-gated (not safe-WP1) so the step-9 verdict does NOT over-credit them; `R_norm` +
`chi_Q` + `N_Q` remain the intact primary falsification load. **Biggest trap = LABEL DRIFT** — surrogate coefficients,
closed-conservation language, or partial-`S_η` structure must never masquerade as physical WP3 target extraction.

## Flows into

- Directive `software/stage1_solver/directives/step8c_conservation_response.md` (guardrails folded in).
- The step-9 Stage-1 verdict: WP3 d-ln targets reported deferred/blocked (Path A); `R_pole/P_2/P_4` reported
  response-gated/not-extracted; `R_norm`/`chi_Q`/`N_Q` are the live primary test.
