# Phase C Adversarial Audit Wrapper

Candidate: `fit_stage_223_mathrm`

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
id: fit_stage_223_mathrm
candidate_key: stage_223_mathrm
dry_run: false
dry_run_id: null
anchor_stages:
- '223'
file_line_citations:
- path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
  line: 130
  role: notes_stage
  stage: '223'
  excerpt: In this stage the target is kept symbolic as \(P_{0,\mathrm{target}}\), because the first question on the primitive
    branch is compatibility rather than full calibration.
parameter_names:
- mathrm
batch_id: null
status: provenance_built
status_timestamps:
  pending: '2026-06-11T08:11:49Z'
  scanned: '2026-06-11T08:11:49Z'
  provenance_built: '2026-06-11T19:17:45Z'
codex_session:
  by_parameter:
    mathrm:
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
  - redteam_adversarial/provenance/fit_stage_223_mathrm__mathrm.yaml
  phase_c_prompt: null
verdict:
  adversarial: null
  adjudication: null
modality_attribution:
- claim_label
- numeric_literal
modality_fragments:
- candidate_key: stage_223_mathrm
  modality: numeric_literal
  anchor_stage: '223'
  parameter_names:
  - mathrm
  reason: target-related numerical literal or closed-form coefficient
  citation:
    path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 130
    role: notes_stage
    stage: '223'
    excerpt: In this stage the target is kept symbolic as \(P_{0,\mathrm{target}}\), because the first question on the primitive
      branch is compatibility rather than full calibration.
- candidate_key: stage_223_mathrm
  modality: claim_label
  anchor_stage: '223'
  parameter_names:
  - mathrm
  reason: claim label or status wording near target-related parameter
  citation:
    path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 130
    role: notes_stage
    stage: '223'
    excerpt: In this stage the target is kept symbolic as \(P_{0,\mathrm{target}}\), because the first question on the primitive
      branch is compatibility rather than full calibration.
phase_b_status: synthesis_complete
```

## Primary Sources

- notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md

## Provenance Slice

```yaml
- schema_version: 1
  candidate_id: fit_stage_223_mathrm
  parameter_name: mathrm
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage_223_mathrm__phase_b__mathrm.md
  agent_synthesis_required: true
  dry_run_fixture: false
  fixture_note: null
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - '223'
  notes_source_opened: true
  source_evidence:
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 123
    role: notes_stage
    excerpt: P_0=P_{0,\mathrm{target}},
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 125
    role: notes_stage
    excerpt: where for the fully calibrated moving-throat branch
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 127
    role: notes_stage
    excerpt: P_{0,\mathrm{target}}=\frac{54Gc_s^5}{5a^5c^5\,\hat m_0^{\,2}}.
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 130
    role: notes_stage
    excerpt: In this stage the target is kept symbolic as \(P_{0,\mathrm{target}}\), because the first question on the primitive
      branch is compatibility rather than full calibration.
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 138
    role: notes_stage
    excerpt: K_{\mathrm{pole}}
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 144
    role: notes_stage
    excerpt: K_{\mathrm{norm}}
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 146
    role: notes_stage
    excerpt: \frac{N_0}{P_{0,\mathrm{target}}}+B_0+Z_0.
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 151
    role: notes_stage
    excerpt: K_{\mathrm{pole}}=K_{\mathrm{norm}},
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 156
    role: notes_stage
    excerpt: \frac{N_0}{P_{0,\mathrm{target}}}
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 164
    role: notes_stage
    excerpt: P_{0,\mathrm{target,compat}}
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 172
    role: notes_stage
    excerpt: First, this is **not** a generic fit condition. It is the exact consequence of trying to hit the same isotropic
      target surface from the conservative one-pole side and the outgoing-normalization side.
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 174
    role: notes_stage
    excerpt: Second, the compatibility equation does **not** determine every coupling. It tells us whether the primitive branch
      wants a normalization target that is even compatible with its own conservative one-pole structure.
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 182
    role: notes_stage
    excerpt: P_{0,\mathrm{target,compat}}
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 190
    role: notes_stage
    excerpt: \frac{P^2/\Delta^2}{P_{0,\mathrm{target}}}
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 252
    role: notes_stage
    excerpt: P_{0,\mathrm{target,compat}}
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 259
    role: notes_stage
    excerpt: K_{\mathrm{compat}}
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 265
    role: notes_stage
    excerpt: D_{0,\mathrm{compat}}
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 267
    role: notes_stage
    excerpt: K_{\mathrm{compat}}-B_0-Z_0
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 279
    role: notes_stage
    excerpt: P_{0,\mathrm{target,compat}}\approx 0.00207.
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 287
    role: notes_stage
    excerpt: Using \(K=K_{\mathrm{compat}}\), the conservative pole census is
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 329
    role: notes_stage
    excerpt: 'Carry forward the same illustrative local barrier benchmark from Stage 222 at \(x=1\):'
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 331
    role: notes_stage
    excerpt: V_{\mathrm{known}}(1)\approx 1.181909222592,
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 335
    role: notes_stage
    excerpt: \Delta V_{\mathrm{req}}(1)\approx 1.081909222592.
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 339
    role: notes_stage
    excerpt: \mathcal R_{Q,*}^{\mathrm{req}}(\eta=0.1)
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 343
    role: notes_stage
    excerpt: \mathcal R_{Q,*}^{\mathrm{req}}(\eta=0.3)
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 347
    role: notes_stage
    excerpt: 'Now scan \(\lambda_W\) **along the exact compatibility surface**, i.e. always resetting \(K\) to \(K_{\mathrm{compat}}(\lambda_W)\).
      The resulting branch-compatible target and wall-like residue/linewidth figures are:'
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 349
    role: notes_stage
    excerpt: '| \(\lambda_W\) | \(P_{0,\mathrm{target,compat}}\) | \(K_{\mathrm{compat}}\) | lower wall \(\mathcal R_Q\) |
      upper wall \(\mathcal R_Q\) |'
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 359
    role: notes_stage
    excerpt: '- increasing \(\lambda_W\) raises the branch-compatible static target \(P_{0,\mathrm{target,compat}}\),'
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 360
    role: notes_stage
    excerpt: '- the same move lowers the required wall stiffness \(K_{\mathrm{compat}}\),'
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 367
    role: notes_stage
    excerpt: At the stricter \(10\%\)-loss benchmark,
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 371
    role: notes_stage
    excerpt: P_{0,\mathrm{target,compat}}
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 376
    role: notes_stage
    excerpt: P_{0,\mathrm{target,compat}}
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 380
    role: notes_stage
    excerpt: At the looser \(30\%\)-loss benchmark,
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 384
    role: notes_stage
    excerpt: P_{0,\mathrm{target,compat}}
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 389
    role: notes_stage
    excerpt: P_{0,\mathrm{target,compat}}
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 411
    role: notes_stage
    excerpt: '### 7.2 The dynamic corridor is not killed automatically by that calibration'
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 413
    role: notes_stage
    excerpt: On the concrete sample slice, once the branch is moved to the compatibility wall stiffness, both wall-like poles
      still clear the stricter \(10\%\)-loss benchmark.
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 459
    role: notes_stage
    excerpt: \frac{N_0}{P_{0,\mathrm{target,compat}}}
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 463
    role: notes_stage
    excerpt: 4. the primitive specialization of \(P_{0,\mathrm{target,compat}}\),
  - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
    line: 466
    role: notes_stage
    excerpt: P_{0,\mathrm{target,compat}},
  origin_claims:
  - parameter: mathrm
    introduced_at_stage: 22
    introduced_at_line: 270
    citation:
      path: notes/stages/moving_throat_pde_stage022_grouped_p2_normalization_bridge.md
      line: 270
      excerpt: '`mhat_0^2 * P_0 = 54 G c_s^5 / (5 a^5 c^5)`.'
  constraints:
  - parameter: mathrm
    constraint_kind: published_target
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
      line: 127
      excerpt: P_{0,\mathrm{target}}=\frac{54Gc_s^5}{5a^5c^5\,\hat m_0^{\,2}}.
  graph_context:
    contexts: []
    graph_gaps:
    - source: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
      attempted_sources:
      - notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
      - research/pde_ledger/notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
      - paper/stages/stage_223.tex
      - research/pde_ledger/paper/stages/stage_223.tex
      graph_gap: true
      reason: no atlas node tied to this source path
  provenance_findings:
  - type: graph_gap
    severity: needs_triage
    summary: The atlas graph wrapper returned no source nodes for one or more primary sources.
    gaps:
    - source: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
      attempted_sources:
      - notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
      - research/pde_ledger/notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
      - paper/stages/stage_223.tex
      - research/pde_ledger/paper/stages/stage_223.tex
      graph_gap: true
      reason: no atlas node tied to this source path
  - type: provenance_gap
    severity: needs_triage
    summary: '''mathrm'' is not a real parameter — it is a LaTeX-command fragment the scanner extracted from \mathrm{...}
      inside P_{0,\mathrm{target}}. The actual referent is P0_target (the 54/5 GR quadrupole external normalization target).
      Classified per its real referent (published_target); the candidate''s parameter_name should be treated as a scanner
      false-positive for P0_target.'
    citations:
    - path: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
      line: 123
      excerpt: P_0=P_{0,\mathrm{target}},
    - path: notes/stages/moving_throat_pde_stage022_grouped_p2_normalization_bridge.md
      line: 256
      excerpt: The universal GR quadrupole target is
  downstream_dependents:
  - '224'
  - '225'
  synthesis_ingested_at: '2026-06-11T19:17:45Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_24/fit_stage_223_mathrm__mathrm.yaml
```

## Benchmarks

External-match checks must use these sourced benchmark entries. Do not adjudicate a match from model memory.

```yaml
[]
```

## Graph Context

```yaml
contexts: []
graph_gaps:
- source: notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
  attempted_sources:
  - notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
  - research/pde_ledger/notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md
  - paper/stages/stage_223.tex
  - research/pde_ledger/paper/stages/stage_223.tex
  graph_gap: true
  reason: no atlas node tied to this source path
```

Output a markdown adversarial report only when this is not a dry run. For dry runs, stop after confirming the prompt has enough context and record no verdict.
