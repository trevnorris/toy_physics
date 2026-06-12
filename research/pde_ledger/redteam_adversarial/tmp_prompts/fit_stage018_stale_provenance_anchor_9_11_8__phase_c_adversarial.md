# Phase C Adversarial Audit Wrapper

Candidate: `fit_stage018_stale_provenance_anchor_9_11_8`

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
id: fit_stage018_stale_provenance_anchor_9_11_8
candidate_key: stage018_stale_provenance_anchor_9_11_8
dry_run: false
dry_run_id: null
anchor_stages:
- 018
file_line_citations:
- path: notes/CHECKPOINT_CONSTANT_PROVENANCE.md
  line: 376
  excerpt: '  carried forward with source anchor from the Stage-18 selected-branch'
  stage: 018
parameter_names:
- coeff_11
- coeff_8
- coeff_9
batch_id: null
status: provenance_built
status_timestamps:
  pending: '2026-06-11T08:21:51Z'
  scanned: '2026-06-11T08:21:51Z'
  provenance_built: '2026-06-11T18:25:24Z'
codex_session:
  by_parameter:
    coeff_11:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    coeff_8:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    coeff_9:
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
  - redteam_adversarial/provenance/fit_stage018_stale_provenance_anchor_9_11_8__coeff_11.yaml
  phase_c_prompt: null
verdict:
  adversarial: null
  adjudication: null
modality_attribution:
- existing_provenance
modality_fragments:
- candidate_key: stage018_stale_provenance_anchor_9_11_8
  modality: existing_provenance
  anchor_stage: 018
  parameter_names:
  - coeff_9
  - coeff_11
  - coeff_8
  reason: The stage-036 provenance entry anchors the carried 9/11/8 coefficients to "Stage-18 selected-branch D/N-normal-form",
    but current stage 018 is the parent-throat action bundle master (M_Sigma/K_Sigma integrals) and the D/N selected-branch
    normal form lives at stage 035 — a stale pre-renumber source anchor, i.e. the recorded provenance is self-contradictory
    and the coefficients' true origin is currently untraceable from this record.
  citation:
    path: notes/CHECKPOINT_CONSTANT_PROVENANCE.md
    line: 376
    excerpt: '  carried forward with source anchor from the Stage-18 selected-branch'
    stage: 018
phase_b_status: synthesis_complete
```

## Primary Sources

- notes/CHECKPOINT_CONSTANT_PROVENANCE.md

## Provenance Slice

```yaml
- schema_version: 1
  candidate_id: fit_stage018_stale_provenance_anchor_9_11_8
  parameter_name: coeff_11
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage018_stale_provenance_anchor_9_11_8__phase_b__coeff_11.md
  agent_synthesis_required: true
  dry_run_fixture: false
  fixture_note: null
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - 018
  notes_source_opened: true
  source_evidence: []
  origin_claims: []
  constraints: []
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings:
  - type: stale_provenance_anchor
    severity: high
    summary: CHECKPOINT_CONSTANT_PROVENANCE records the carried coefficients 9/11/8 as anchored to a 'Stage-18 selected-branch
      D/N-normal-form', but current stage 018 is the parent-throat action bundle master (M_Sigma/K_Sigma branch integrals)
      with no D/N selected-branch normal form; the unpadded 'Stage-18' is a pre-renumber anchor (the D/N selected-branch normal
      form lives near stage 035). The recorded source anchor for coeff_11 therefore points at the wrong stage.
    citations:
    - path: notes/CHECKPOINT_CONSTANT_PROVENANCE.md
      line: 376
      excerpt: '  carried forward with source anchor from the Stage-18 selected-branch'
    - path: notes/CHECKPOINT_CONSTANT_PROVENANCE.md
      line: 377
      excerpt: '  D/N-normal-form coefficients and the `8 / (pi^2 A)` dimensionless scaling'
    - path: paper/stages/stage_018.tex
      line: 1
      excerpt: '\section[Stage 018]{Stage 018: Parent throat action bundle master}'
    - path: scripts/moving_throat_pde_stage018_parent_throat_action_bundle_master_sympy_audit.py
      line: 2
      excerpt: '"""Master-note audit for step_16_parent_throat_action_bundle_master_notes.md."""'
  - type: provenance_gap
    severity: needs_triage
    summary: No file under notes/ establishes the origin or constraint of coeff_11 (the carried '9/11/8' coefficient) at stage
      018. Stage 018 has no per-stage notes/stages/ file; its listed note source (notes/em_projected/step_16_parent_throat_action_bundle_master_notes.md
      per STAGE_PROVENANCE_INDEX line 28) does not exist in the repo (notes/em_projected/ is absent). The only notes/ mention
      (CHECKPOINT_CONSTANT_PROVENANCE.md:376) is the stale anchor itself and gives no clean value genealogy, so constraint_kind
      cannot be assigned.
    citations:
    - path: notes/STAGE_PROVENANCE_INDEX.md
      line: 28
      excerpt: '| 018 | `paper/stages/stage_018.tex` | `notes/em_projected/step_16_parent_throat_action_bundle_master_notes.md`
        |'
    - path: notes/CHECKPOINT_CONSTANT_PROVENANCE.md
      line: 375
      excerpt: '- `9`, `11`, `8`'
  downstream_dependents:
  - '030'
  synthesis_ingested_at: '2026-06-11T18:25:24Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_02/fit_stage018_stale_provenance_anchor_9_11_8__coeff_11.yaml
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
