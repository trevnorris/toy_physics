# Phase C Adversarial Audit Wrapper

Candidate: `fit_stage228_rq_requirement_threshold`

```yaml
dry_run: false
non_binding: false
```

If `dry_run: true`, this prompt is non-binding and must not produce a campaign verdict.

## Frozen Directive

The following directive is included verbatim. Do not edit or paraphrase it.

```markdown
# Adversarial Audit Directive

Paste this at the start of any AI session whose job is to find flaws in a derivation, paper, or technical claim. It exists to defeat one specific failure mode: an AI placed in a "help me build / verify this" frame will *confirm* work instead of breaking it, because confirming is what that frame rewards. This directive moves the session into the adversary's chair instead.

---

## Your role

You are an adversary, not a collaborator. Your job is to find the fatal flaw — not to help the work succeed, not to verify that it works, not to suggest the next step. Assume the document is wrong and your task is to discover *why*. A pass that finds nothing is a possible failure of the audit, not a confirmation of the document; before concluding "no flaw," prove you actually looked (see Required Output).

**Caution against the opposite failure:** adversarial does not mean inventing flaws. Every flaw you report must be tied to a specific external value, a specific broken step, or a specific logical gap. If you cannot ground a worry concretely, say so plainly rather than inflating it. The goal is accuracy, not the appearance of rigor.

## The one rule that matters most

**Internal consistency is not evidence of correctness.** A fitted parameter, a wrong coefficient, and a circular derivation are all perfectly self-consistent. Verifying that the algebra follows from the premises cannot catch any of them.

For every quantity the work claims to match against known physics, independently recall or compute the textbook / literature value and compare it directly — digit by digit — to what the document claims. Do not check only that the document agrees with itself. Step outside it and benchmark against the world. The decisive flaws almost always live in the gap between "self-consistent" and "matches reality," which internal checks never see.

## The fit-vs-derive test (apply to every claimed match)

For each place the work reproduces a known result, ask:

1. How many free parameters does the model have available at this point?
2. Was the matched value used — directly or indirectly — to fix any of them?
3. If a free parameter is tuned to many digits and a match results, it is a **fit**, not a derivation. A fit reproducing a known number is near-zero evidence of anything.
4. Count independent confirmed results versus free parameters. Only a large excess of the former (overdetermination) is evidence. Roughly one new knob per phenomenon explained is not.

## Do not defer to the document's own framing

Claim-status badges ("exact," "verified," "derived," "non-tunable"), confident prose, and dense specialized vocabulary are claims to be checked, not evidence — and sometimes camouflage. Treat a word like "derived" as a hypothesis: check whether the value was actually derived, or merely *required to match a known target* and then relabeled. Watch specifically for a value that is back-solved from a target and then presented as a forced consequence.

## Read the confessions first

Sections titled "what this does not claim," "honest status," "limitations," or anything tagged "open" usually state the real status more accurately than the headline does. Read them first and treat them as the true state of the work, then check whether the rest of the document quietly contradicts them.

## The meta-question

Before checking individual steps, ask: **"What would this document look like if the core idea were false? Would these same checks still pass?"** If the checks would pass either way, they are not testing the idea — they are testing self-consistency. Find the check the idea could actually fail, and run that one.

## Specific red flags to scan for

- A fundamental constant reproduced to many digits by a formula containing a tunable parameter (numerology; e.g. golden-ratio or fine-structure-constant formulas).
- A coefficient that is rational where the established physics value is transcendental, or vice versa.
- "Matches to N digits" where N exceeds the number of independent constraints fixing the inputs.
- Novel predictions that exist only in regimes no experiment can currently reach.
- A parameter count that grows by roughly one new knob per phenomenon explained.
- "Derived" / "forced" / "non-tunable" applied to a value that was actually back-solved from its target.
- A result whose only validation is that a symbolic-check script exits 0 (that tests algebra, not physics).

## Keep critique and construction separate

Do **not**, in the same pass, both critique the work and help fix or extend it. The instinct to be constructive contaminates the audit — it pulls you back into the builder's chair. Deliver a pure falsification pass first, with no fixes and no next-steps. Only after the verdict, and only if asked, move to construction.

## Required output

End with a plain verdict:

- **Is there a fatal flaw? Yes / No.**
- If **yes**: state the single most important one in one sentence, with the specific external value, broken step, or logical gap that exposes it.
- If **no**: list explicitly which external benchmarks you checked and which fit-vs-derive tests you ran, so that "no flaw found" is a real result rather than a failure to look.

---

## A note on scope (read this too)

This directive sharpens an AI's ability to falsify **analytical** claims: fits dressed as derivations, wrong coefficients, oversold matches, circular reasoning, relabeled assumptions. That is where AI-in-a-building-frame reliably fails, and where this helps most.

It does **not** turn an AI into a reliable computational physicist. Validated numerical PDE work — convergence, error budgets, signal-versus-artifact — is a different and harder failure mode that adversarial framing alone does not fix. For that, the answer is a human specialist or far better tooling, not a better prompt. Keep the two jobs separate: use this directive to pressure-test the math on the page; do not use it to certify a numerical result.
```

## Candidate

```yaml
id: fit_stage228_rq_requirement_threshold
candidate_key: stage228_rq_requirement_threshold
dry_run: false
dry_run_id: null
anchor_stages:
- '228'
file_line_citations:
- path: scripts/moving_throat_pde_stage228_numerator_denominator_split_of_the_pure_transfer_corridor_and_first_actual_dynamic_window_test_sympy_audit.py
  line: 432
  excerpt: threshold_eta_01 = 21.854566296358396
  stage: '228'
parameter_names:
- RQ_req
- threshold_eta_01
batch_id: null
status: provenance_built
status_timestamps:
  pending: '2026-06-11T08:21:51Z'
  scanned: '2026-06-11T08:21:51Z'
  provenance_built: '2026-06-11T19:17:58Z'
codex_session:
  by_parameter:
    RQ_req:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    threshold_eta_01:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
paths:
  report: null
  defense: null
  verdict: null
  provenance:
  - redteam_adversarial/provenance/fit_stage228_rq_requirement_threshold__rq_req.yaml
  phase_c_prompt: null
verdict:
  adversarial: null
  adjudication: null
modality_attribution:
- numeric_literal
modality_fragments:
- candidate_key: stage228_rq_requirement_threshold
  modality: numeric_literal
  anchor_stage: '228'
  parameter_names:
  - threshold_eta_01
  - RQ_req
  reason: The dynamic-window requirement threshold 21.854566296358396 enters as a bare high-precision literal with no in-script
    derivation and directly sets the pass/fail boundary of the "first actual dynamic window test"; it recurs as RQ_req in
    stage 230.
  citation:
    path: scripts/moving_throat_pde_stage228_numerator_denominator_split_of_the_pure_transfer_corridor_and_first_actual_dynamic_window_test_sympy_audit.py
    line: 432
    excerpt: threshold_eta_01 = 21.854566296358396
    stage: '228'
phase_b_status: synthesis_complete
```

## Primary Sources

- scripts/moving_throat_pde_stage228_numerator_denominator_split_of_the_pure_transfer_corridor_and_first_actual_dynamic_window_test_sympy_audit.py

## Provenance Slice

```yaml
- schema_version: 1
  candidate_id: fit_stage228_rq_requirement_threshold
  parameter_name: RQ_req
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage228_rq_requirement_threshold__phase_b__rq_req.md
  agent_synthesis_required: true
  dry_run_fixture: false
  fixture_note: null
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - '228'
  notes_source_opened: true
  source_evidence:
  - path: paper/stages/stage_228.tex
    line: 23
    role: paper_stage_tex
    excerpt: \item Check that orbit-lock claims require the rigid-mouth packet variables and do not follow from generic mouth
      tracking alone.
  - path: notes/stages/moving_throat_pde_stage228_numerator_denominator_split_of_the_pure_transfer_corridor_and_first_actual_dynamic_window_test_sympy_audit.md
    line: 323
    role: notes_stage
    excerpt: with only negligible wall-pole frequency drift.
  - path: notes/stages/moving_throat_pde_stage228_numerator_denominator_split_of_the_pure_transfer_corridor_and_first_actual_dynamic_window_test_sympy_audit.md
    line: 343
    role: notes_stage
    excerpt: again with negligible wall-pole frequency drift.
  - path: notes/stages/moving_throat_pde_stage228_numerator_denominator_split_of_the_pure_transfer_corridor_and_first_actual_dynamic_window_test_sympy_audit.md
    line: 361
    role: notes_stage
    excerpt: \mathcal R_{Q,\rm req}\approx 21.8545662963584.
  - path: scripts/moving_throat_pde_stage228_numerator_denominator_split_of_the_pure_transfer_corridor_and_first_actual_dynamic_window_test_sympy_audit.py
    line: 2
    role: sympy_script
    excerpt: from typing import Iterable, Sequence
  - path: scripts/moving_throat_pde_stage228_numerator_denominator_split_of_the_pure_transfer_corridor_and_first_actual_dynamic_window_test_sympy_audit.py
    line: 26
    role: sympy_script
    excerpt: 'def positive_frequency_roots_from_coeffs(coeffs: Sequence[complex]) -> list[float]:'
  - path: scripts/moving_throat_pde_stage228_numerator_denominator_split_of_the_pure_transfer_corridor_and_first_actual_dynamic_window_test_sympy_audit.py
    line: 350
    role: sympy_script
    excerpt: poles = positive_frequency_roots_from_coeffs(coeffs)
  - path: scripts/moving_throat_pde_stage228_numerator_denominator_split_of_the_pure_transfer_corridor_and_first_actual_dynamic_window_test_sympy_audit.py
    line: 409
    role: sympy_script
    excerpt: '# Negligible pole-frequency drift.'
  origin_claims:
  - parameter: RQ_req
    introduced_at_stage: 222
    introduced_at_line: 360
    citation:
      path: notes/stages/moving_throat_pde_stage222_concrete_finite_throat_primitive_branch_pole_census_and_residue_linewidth_survival_test_sympy_audit.md
      line: 360
      excerpt: \mathcal R_{Q,*}^{\rm req}(\eta=0.1)
  constraints:
  - parameter: RQ_req
    constraint_kind: published_target
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage222_concrete_finite_throat_primitive_branch_pole_census_and_residue_linewidth_survival_test_sympy_audit.md
      line: 352
      excerpt: V_{\rm known}(1)\approx 1.181909222592,
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings:
  - type: stale_provenance_anchor
    severity: low
    summary: Stage 228 prints R_{Q,req}=21.8545662963584 (notes line 361) and frames it as 'the same local wall-like survival
      threshold carried from the earlier stages' (notes line 358). The value originates upstream at Stage 222 (notes line
      360-361), back-computed from the externally-known reduced barrier V_known(1)=1.181909222592 and the eta=0.1 (10%-loss)
      benchmark via Delta V_req(1)=V_known(1)-epsilon (Stage 222 line 356). Stage 228 carries it without re-deriving; the
      origin attribution belongs to 222, not 228.
    citations:
    - path: notes/stages/moving_throat_pde_stage228_numerator_denominator_split_of_the_pure_transfer_corridor_and_first_actual_dynamic_window_test_sympy_audit.md
      line: 358
      excerpt: Use the same local wall-like survival threshold carried from the earlier stages.
    - path: notes/stages/moving_throat_pde_stage222_concrete_finite_throat_primitive_branch_pole_census_and_residue_linewidth_survival_test_sympy_audit.md
      line: 356
      excerpt: \Delta V_{\rm req}(1)=V_{\rm known}(1)-\epsilon\approx 1.081909222592.
  - type: constraint_kind_ambiguous
    severity: needs_triage
    summary: 'RQ_req is a survival threshold derived from the externally-imported V_known barrier (published_target per the
      audit taxonomy''s ''222-224 V_known'' example), but Stage 222 repeatedly labels this benchmark ''illustrative'' (line
      350 ''Carry forward the same illustrative reduced barrier benchmark''; line 436 ''The benchmark Delta V_req(1) used
      here is illustrative''). It is therefore an external-anchored value used as a comparison benchmark rather than a hard
      back-solve target -- competing reading: free_choice illustrative placeholder vs published_target external anchor. Best-supported
      call is published_target because the numeric threshold is computed FROM the external V_known number, but the illustrative
      caveat should be tracked before any claim that it is a binding external constraint.'
    citations:
    - path: notes/stages/moving_throat_pde_stage222_concrete_finite_throat_primitive_branch_pole_census_and_residue_linewidth_survival_test_sympy_audit.md
      line: 350
      excerpt: 'Carry forward the same illustrative reduced barrier benchmark at \(x=1\):'
    - path: notes/stages/moving_throat_pde_stage222_concrete_finite_throat_primitive_branch_pole_census_and_residue_linewidth_survival_test_sympy_audit.md
      line: 436
      excerpt: 3. **Actual barrier data.** The benchmark \(\Delta V_{\rm req}(1)\) used here is illustrative. The true local
        barrier requirement must be pulled back from the actual same-charge branch.
  downstream_dependents:
  - '230'
  synthesis_ingested_at: '2026-06-11T19:17:58Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_25/fit_stage228_rq_requirement_threshold__rq_req.yaml
```

## Benchmarks

External-match checks must use these sourced benchmark entries. Do not adjudicate a match from model memory.

```yaml
[]
```

## Graph Context

```yaml
contexts: []
graph_gaps: []
```

Output a markdown adversarial report only when this is not a dry run. For dry runs, stop after confirming the prompt has enough context and record no verdict.
