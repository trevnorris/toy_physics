# Phase C Adversarial Audit Wrapper

Candidate: `fit_stage022_invariant_product_must_match_universal`

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
id: fit_stage022_invariant_product_must_match_universal
candidate_key: stage022_invariant_product_must_match_universal
dry_run: false
dry_run_id: null
anchor_stages:
- '022'
file_line_citations:
- path: paper/stages/stage_022.tex
  line: 7
  excerpt: \stagefield{Purpose}{Stage~022 translates the Stages~004--021 projection-first Maxwell packet, through its matched
    one-port outgoing normal form, into the grouped real \(P_2\) response language.  It converts conservative operator moments
    \(D_{An}\) into normalized response moments \(u_n^{(A)}\), tracks the outgoing-transfer prefactor \(P_{An}\), and isolates
    the exact invariant product that must match the universal quadrupole normalization.}
  stage: '022'
parameter_names:
- N_Q_universal
- P_An
- u_n_A
batch_id: null
status: provenance_built
status_timestamps:
  pending: '2026-06-11T08:21:51Z'
  scanned: '2026-06-11T08:21:51Z'
  provenance_built: '2026-06-11T18:25:24Z'
codex_session:
  by_parameter:
    N_Q_universal:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    P_An:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    u_n_A:
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
  - redteam_adversarial/provenance/fit_stage022_invariant_product_must_match_universal__n_q_universal.yaml
  phase_c_prompt: null
verdict:
  adversarial: null
  adjudication: null
modality_attribution:
- claim_label
modality_fragments:
- candidate_key: stage022_invariant_product_must_match_universal
  modality: claim_label
  anchor_stage: '022'
  parameter_names:
  - P_An
  - u_n_A
  - N_Q_universal
  reason: Declares an "exact invariant product that must match the universal quadrupole normalization"; the must-match phrasing
    defines a target-matching slot even though the claimstatus marks the equality Open.
  citation:
    path: paper/stages/stage_022.tex
    line: 7
    excerpt: \stagefield{Purpose}{Stage~022 translates the Stages~004--021 projection-first Maxwell packet, through its matched
      one-port outgoing normal form, into the grouped real \(P_2\) response language.  It converts conservative operator moments
      \(D_{An}\) into normalized response moments \(u_n^{(A)}\), tracks the outgoing-transfer prefactor \(P_{An}\), and isolates
      the exact invariant product that must match the universal quadrupole normalization.}
    stage: '022'
phase_b_status: synthesis_complete
```

## Primary Sources

- paper/stages/stage_022.tex

## Provenance Slice

```yaml
- schema_version: 1
  candidate_id: fit_stage022_invariant_product_must_match_universal
  parameter_name: N_Q_universal
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage022_invariant_product_must_match_universal__phase_b__n_q_universal.md
  agent_synthesis_required: true
  dry_run_fixture: false
  fixture_note: null
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - '022'
  notes_source_opened: true
  source_evidence:
  - path: paper/stages/stage_022.tex
    line: 5
    role: paper_stage_tex
    excerpt: \claimstatus{The bridge formulas are \StatusExactClosure{} inside the grouped response normalization convention.  The
      final equality to the universal normalization is a \StatusOpen{} branch-realization test.}
  - path: paper/stages/stage_022.tex
    line: 7
    role: paper_stage_tex
    excerpt: \stagefield{Purpose}{Stage~022 translates the Stages~004--021 projection-first Maxwell packet, through its matched
      one-port outgoing normal form, into the grouped real \(P_2\) response language.  It converts conservative operator moments
      \(D_{An}\) into normalized response moments \(u_n^{(A)}\), tracks the outgoing-transfer prefactor \(P_{An}\), and isolates
      the exact invariant product that must match the universal quadrupole normalization.}
  - path: paper/stages/stage_022.tex
    line: 9
    role: paper_stage_tex
    excerpt: \stagefield{Inputs}{The inputs are one conservative grouped-lane operator \(D_A^{\rm cons}(\omega)\), one outgoing
      transfer factor \(N_A(\omega)\) exported by the matched Stage~021 one-port normal form, the compact outgoing \(l=2\)
      fingerprint, and the grouped trace/anomaly map.}
  - path: paper/stages/stage_022.tex
    line: 86
    role: paper_stage_tex
    excerpt: \paragraph{4. Multiplication by the compact outgoing fingerprint.}
  - path: paper/stages/stage_022.tex
    line: 117
    role: paper_stage_tex
    excerpt: \paragraph{5. Universal normalization product.}
  - path: notes/stages/moving_throat_pde_stage022_grouped_p2_normalization_bridge.md
    line: 16
    role: notes_stage
    excerpt: '> lift the one-lane Stage-021 result to the grouped real `20/21/22` bundle, convert the Stage-003 and Stages-004--021
      conservative operator moments into the normalized grouped-response moments used by the 2.5PN package, and isolate the
      exact normalization product that still has to hit the universal quadrupole target.'
  - path: notes/stages/moving_throat_pde_stage022_grouped_p2_normalization_bridge.md
    line: 256
    role: notes_stage
    excerpt: The universal GR quadrupole target is
  - path: notes/stages/moving_throat_pde_stage022_grouped_p2_normalization_bridge.md
    line: 309
    role: notes_stage
    excerpt: At the universal point-particle 2.5PN level, the only internal quantity that matters is the isotropic static
      outgoing prefactor `P_0 = N_0 / D_0` (or `mhat_0^2 P_0` in invariant form).
  - path: notes/stages/moving_throat_pde_stage022_grouped_p2_normalization_bridge.md
    line: 334
    role: notes_stage
    excerpt: 4. and test whether the invariant normalization product `mhat_0^2 P_0` lands on the universal target.
  origin_claims:
  - parameter: N_Q_universal
    introduced_at_stage: 18
    introduced_at_line: 270
    citation:
      path: notes/stages/moving_throat_pde_stage022_grouped_p2_normalization_bridge.md
      line: 270
      excerpt: '`mhat_0^2 * P_0 = 54 G c_s^5 / (5 a^5 c^5)`.'
  constraints:
  - parameter: N_Q_universal
    constraint_kind: published_target
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage022_grouped_p2_normalization_bridge.md
      line: 256
      excerpt: The universal GR quadrupole target is
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings:
  - type: derive_vs_fit_mismatch
    severity: high
    summary: 'The invariant normalization product mhat_0^2 P_0 (= N_Q_universal) is REQUIRED to match the universal GR quadrupole
      value: the note defines gamma_GR = 2G/(5c^5) as ''The universal GR quadrupole target'' and back-solves the product to
      mhat_0^2 P_0 = 54Gc_s^5/(5a^5c^5). The whole point of the candidate name (''invariant_product_must_match_universal'')
      is the fit-to-external-target. This is a published_target constraint. It is honestly disclosed in the card claimstatus
      as a StatusOpen branch-realization test (stage_022.tex:5), so not overclaimed, but it must not be classed internal_consistency.'
    citations:
    - path: notes/stages/moving_throat_pde_stage022_grouped_p2_normalization_bridge.md
      line: 264
      excerpt: Including the source-map factor `mhat_0`, the invariant normalization condition is therefore
    - path: paper/stages/stage_022.tex
      line: 5
      excerpt: The final equality to the universal normalization is a \StatusOpen{} branch-realization test.
  downstream_dependents:
  - 019
  synthesis_ingested_at: '2026-06-11T18:25:24Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_02/fit_stage022_invariant_product_must_match_universal__n_q_universal.yaml
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
