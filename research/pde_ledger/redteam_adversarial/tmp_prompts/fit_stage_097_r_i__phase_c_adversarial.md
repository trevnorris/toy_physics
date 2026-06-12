# Phase C Adversarial Audit Wrapper

Candidate: `fit_stage_097_r_i`

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
id: fit_stage_097_r_i
candidate_key: stage_097_r_i
dry_run: false
dry_run_id: null
anchor_stages:
- 097
file_line_citations:
- path: redteam/pass2/reports/stage_097.md
  line: 103
  role: pass2_stage_report
  stage: 097
  excerpt: Verdict is **clean**. I read the card, notes, and the part-IV appendix rows before opening the scripts; the script's
    verified claim (the K2/K4/Gamma5 reductions, the GR target, and the four-way R_i = N_Q - 1 collapse) matches the paper's
    stated deliverables exactly, with every constant agreeing. The R_i checks are non-tautological (target coefficients are
    built from the unsubstituted target K0, so a wrong power of Omega would survive), the positivity assumptions are physically
    warranted and required for the half-integer powers, the Mathematica engine independently re-derives the moments via series
    extraction, and both engines agree on every value with zero residuals. Output `.txt` files are newer than their scripts
    (see Value Reconciliation). Attacks on tautology, transliteration, symbol domain, parity, and constants all failed to
    find a defect.
parameter_names:
- R_i
batch_id: null
status: provenance_built
status_timestamps:
  pending: '2026-06-11T08:11:49Z'
  scanned: '2026-06-11T08:11:49Z'
  provenance_built: '2026-06-11T18:51:37Z'
codex_session:
  by_parameter:
    R_i:
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
  - redteam_adversarial/provenance/fit_stage_097_r_i__r_i.yaml
  phase_c_prompt: null
verdict:
  adversarial: null
  adjudication: null
modality_attribution:
- existing_provenance
modality_fragments:
- candidate_key: stage_097_r_i
  modality: existing_provenance
  anchor_stage: 097
  parameter_names:
  - R_i
  reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
  citation:
    path: redteam/pass2/reports/stage_097.md
    line: 103
    role: pass2_stage_report
    stage: 097
    excerpt: Verdict is **clean**. I read the card, notes, and the part-IV appendix rows before opening the scripts; the script's
      verified claim (the K2/K4/Gamma5 reductions, the GR target, and the four-way R_i = N_Q - 1 collapse) matches the paper's
      stated deliverables exactly, with every constant agreeing. The R_i checks are non-tautological (target coefficients
      are built from the unsubstituted target K0, so a wrong power of Omega would survive), the positivity assumptions are
      physically warranted and required for the half-integer powers, the Mathematica engine independently re-derives the moments
      via series extraction, and both engines agree on every value with zero residuals. Output `.txt` files are newer than
      their scripts (see Value Reconciliation). Attacks on tautology, transliteration, symbol domain, parity, and constants
      all failed to find a defect.
phase_b_status: synthesis_complete
```

## Primary Sources

- redteam/pass2/reports/stage_097.md

## Provenance Slice

```yaml
- schema_version: 1
  candidate_id: fit_stage_097_r_i
  parameter_name: R_i
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage_097_r_i__phase_b__r_i.md
  agent_synthesis_required: true
  dry_run_fixture: false
  fixture_note: null
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - 097
  notes_source_opened: true
  source_evidence: []
  origin_claims:
  - parameter: R_i
    introduced_at_stage: 97
    introduced_at_line: 104
    citation:
      path: notes/stages/moving_throat_pde_stage097_single_normalization_defect.md
      line: 104
      excerpt: '`R_0 := Kbar_0 / Kbar_0^target - 1,`'
  constraints:
  - parameter: R_i
    constraint_kind: published_target
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage097_single_normalization_defect.md
      line: 114
      excerpt: '`R_0 = R_2 = R_4 = R_5 = N_Q - 1`.'
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings:
  - type: provenance_gap
    severity: low
    summary: The mechanical bundle for fit_stage_097_r_i had EMPTY source_evidence. The parameter R_i resolves to the residual
      defects R_0/R_2/R_4/R_5 := Kbar_i/Kbar_i^target - 1, documented at stage097 notes lines 104-114. Each R_i measures distance
      from a GR target (published_target); the result that all four collapse to the single N_Q-1 is internal_consistency.
      Origin and constraint supplied by opening the per-stage note despite the empty bundle.
    citations:
    - path: notes/stages/moving_throat_pde_stage097_single_normalization_defect.md
      line: 110
      excerpt: '`R_5 := Gammabar_5 / (2G/(5c^5)) - 1,`'
  downstream_dependents:
  - 099
  - '100'
  synthesis_ingested_at: '2026-06-11T18:51:37Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_10/fit_stage_097_r_i__r_i.yaml
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
