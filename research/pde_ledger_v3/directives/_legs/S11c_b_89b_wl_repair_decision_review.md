# Decision review — S11c-b #89b WL repair directive (orchestrator-written)

You are one of two independent decision-review legs on an ORCHESTRATOR-WRITTEN repair directive, BEFORE the
builder (Codex) launches. Your job is to find where the directive is wrong, incomplete, or leaks an answer —
so a build round is not wasted. This is a document review; derive from the engine source, not from prose.

## Artifacts
- The repair directive to review: `research/pde_ledger_v3/directives/S11c_b_89b_wl_operator_repair.md`
- The engine it governs (working tree, uncommitted): `research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl`
- The build-review that produced the three findings (context): `directives/_measurements/S11c_b_89b_wl_operator_build_review.md`

## What to check (verify each against the actual engine code, cite line numbers)
1. **Defect A completeness — the dangerous one.** The directive requires every emitted operator-family object
   AND every operand compared against it to be the activate-then-reduce (derivative-normal, jet-retaining) form,
   while `KERNEL_SOURCE_OPERATOR` stays the un-reduced live operator. Independently enumerate EVERY reader of
   `["OPERATOR"]` / `["ORIGINS"]` and every per-depth / ablation / limit / cross-model operand in the engine, and
   decide: would applying the fix as described leave any comparison comparing UNLIKE forms (one activated, one
   still reduced-from-un-activated), or silently break any consumer? Name any site the directive fails to cover.
   Also confirm the directive is right that the §3c weak split consumes `KERNEL_SOURCE_OPERATOR` (un-reduced) and
   NOT the emitted `OPERATOR` — if that is wrong, the whole fix is wrong.
2. **Defect A correctness.** Is "activate the outer constrained divergences before the final background reduction"
   actually the correct object (the one that retains the mixed/higher background jets), and is the claim that
   `MU_THETA` / faces / divergence-form-source already use that route true in the code? Is the discarded
   `activatedOperator` really the right object?
3. **Defect B.** Confirm the `Times[factors__]`/`Plus[terms__]` patterns genuinely fail to attach (Flat/OneIdentity)
   and that the proposed property (walk composites to primitives; residual MOVES under a primitive-atom dimension
   mutation) is a decisive, non-tautological acceptance.
4. **Defect C.** Confirm the `corrupted = If[branch === First[branches], …, base]` change makes the residual
   `base − base = 0` for non-first branches, and that the proposed fix (genuine corruption, or an explicit
   `VALIDATED_ON_REPRESENTATIVE_BRANCH` marker) removes the tautology without hiding a real failure.
5. **Answer leakage (rule 5).** Does the directive leak an expected value — a coefficient, a jet count, a target
   the builder could iterate toward? The acceptance must be a residual between two independent routes / a "moves
   under mutation" property, never "equals <value>". Flag any leak.
6. **Scope / safety.** Does the directive correctly forbid touching the kernel path, the un-freeze mechanism, the
   tractability activation, the deferral gate, and every other emit? Is anything it asks for likely to regress a
   currently-correct object?

## Method + output
- Ground every claim in a line number from the engine. A prose assertion with no code citation is discarded.
- You may (optionally) confirm a mechanism with a tiny reduced-scale check (`basisRepresentativeIndices={16}`,
  `timeout 600`, one kernel, /tmp only, kill leftover kernels by exact pid) — but a document-level code reading
  is sufficient for most items.
- Report a numbered list: `severity — directive section / file:line — the gap or error — why it matters — the
  code evidence`. Severity ∈ {MUST-FIX-BEFORE-BUILD, SHOULD-FIX, NIT}. Then
  `VERDICT: N issues (M must-fix)` or `VERDICT: DIRECTIVE READY`. If you clear it, name the specific consumer
  sites you enumerated for defect A so the completeness claim is auditable.
