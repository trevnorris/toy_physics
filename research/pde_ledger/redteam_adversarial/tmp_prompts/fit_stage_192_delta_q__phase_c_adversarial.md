# Phase C Adversarial Audit Wrapper

Candidate: `fit_stage_192_delta_q`

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
id: fit_stage_192_delta_q
candidate_key: stage_192_delta_q
dry_run: false
dry_run_id: null
anchor_stages:
- '192'
file_line_citations:
- path: paper/stages/stage_192.tex
  line: 13
  role: paper_stage_tex
  stage: '192'
  excerpt: \stagefield{Derivation ledger}{The derivation states the exact isotropic target surface, attaches the outgoing
    \(l=2\) fingerprint and source-map factor, proves higher odd terms are irrelevant at the current order, and collapses
    Packet A to \(\Delta_Q=\chi_Q-1\).}
parameter_names:
- Delta_Q
- chi_Q
- matched_fingerprint_value
batch_id: null
status: provenance_built
status_timestamps:
  pending: '2026-06-11T08:11:49Z'
  scanned: '2026-06-11T08:11:49Z'
  provenance_built: '2026-06-11T19:18:45Z'
codex_session:
  by_parameter:
    Delta_Q:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    chi_Q:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    matched_fingerprint_value:
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
  - redteam_adversarial/provenance/fit_stage_192_delta_q__delta_q.yaml
  phase_c_prompt: null
verdict:
  adversarial: null
  adjudication: null
modality_attribution:
- claim_label
- numeric_literal
modality_fragments:
- candidate_key: stage_192_delta_q
  modality: numeric_literal
  anchor_stage: '192'
  parameter_names:
  - Delta_Q
  - chi_Q
  - matched_fingerprint_value
  reason: target-related numerical literal or closed-form coefficient
  citation:
    path: paper/stages/stage_192.tex
    line: 13
    role: paper_stage_tex
    stage: '192'
    excerpt: \stagefield{Derivation ledger}{The derivation states the exact isotropic target surface, attaches the outgoing
      \(l=2\) fingerprint and source-map factor, proves higher odd terms are irrelevant at the current order, and collapses
      Packet A to \(\Delta_Q=\chi_Q-1\).}
- candidate_key: stage_192_delta_q
  modality: claim_label
  anchor_stage: '192'
  parameter_names:
  - Delta_Q
  - chi_Q
  - matched_fingerprint_value
  reason: claim label or status wording near target-related parameter
  citation:
    path: paper/stages/stage_192.tex
    line: 13
    role: paper_stage_tex
    stage: '192'
    excerpt: \stagefield{Derivation ledger}{The derivation states the exact isotropic target surface, attaches the outgoing
      \(l=2\) fingerprint and source-map factor, proves higher odd terms are irrelevant at the current order, and collapses
      Packet A to \(\Delta_Q=\chi_Q-1\).}
phase_b_status: synthesis_complete
```

## Primary Sources

- paper/stages/stage_192.tex

## Provenance Slice

```yaml
- schema_version: 1
  candidate_id: fit_stage_192_delta_q
  parameter_name: Delta_Q
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage_192_delta_q__phase_b__delta_q.md
  agent_synthesis_required: true
  dry_run_fixture: false
  fixture_note: null
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - '192'
  notes_source_opened: true
  source_evidence:
  - path: paper/stages/stage_192.tex
    line: 9
    role: paper_stage_tex
    excerpt: \stagefield{Inputs}{This stage imports the isotropic grouped-\(P_2\) one-pole front end, the compact outgoing
      \(l=2\) DtN fingerprint, and the natural point-particle source-map reduction.  It is the Packet-A branch-side finish
      line.}
  - path: paper/stages/stage_192.tex
    line: 13
    role: paper_stage_tex
    excerpt: \stagefield{Derivation ledger}{The derivation states the exact isotropic target surface, attaches the outgoing
      \(l=2\) fingerprint and source-map factor, proves higher odd terms are irrelevant at the current order, and collapses
      Packet A to \(\Delta_Q=\chi_Q-1\).}
  - path: paper/stages/stage_192.tex
    line: 20
    role: paper_stage_tex
    excerpt: \item Check the isotropic one-pole target surface before using \(\Delta_Q=\chi_Q-1\).
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 19
    role: notes_stage
    excerpt: \Delta_{\rm orbit}=(q_{\rm tr},q_{\rm nt},q_\eta).
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 36
    role: notes_stage
    excerpt: (\Delta_T,\Delta_{K_\eta},\Delta_\mu),
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 48
    role: notes_stage
    excerpt: (\Delta_T,\Delta_{K_\eta},\Delta_\mu),
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 52
    role: notes_stage
    excerpt: Q_{\rm quot}\Delta\mathbf x=0
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 54
    role: notes_stage
    excerpt: M_*\Delta\mathbf x=0.
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 66
    role: notes_stage
    excerpt: \Delta\mathbf x:=
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 68
    role: notes_stage
    excerpt: \Delta_\lambda\\
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 69
    role: notes_stage
    excerpt: \Delta_c\\
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 70
    role: notes_stage
    excerpt: \Delta_\gamma\\
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 71
    role: notes_stage
    excerpt: \Delta_U\\
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 72
    role: notes_stage
    excerpt: \Delta_{K_\eta}\\
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 73
    role: notes_stage
    excerpt: \Delta_W\\
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 74
    role: notes_stage
    excerpt: \Delta_\mu\\
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 75
    role: notes_stage
    excerpt: \Delta_T
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 93
    role: notes_stage
    excerpt: =M_*\,\Delta\mathbf x,
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 101
    role: notes_stage
    excerpt: 0 & 1+\delta_{U,*} & 1+\delta_{U,*} & -(2+\chi_{0,*}+\delta_{U,*}) & 0 & 0 & 0 & 1+\chi_{0,*}\\[4pt]
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 112
    role: notes_stage
    excerpt: (\Delta_\lambda,\Delta_c,\Delta_\gamma,\Delta_U,\Delta_W),
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 116
    role: notes_stage
    excerpt: (\Delta_T,\Delta_{K_\eta},\Delta_\mu).
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 210
    role: notes_stage
    excerpt: \boxed{\Delta\mathbf x_{\rm quot}:=S_{(T,K_\eta,\mu)}\,\mathbf q.}
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 240
    role: notes_stage
    excerpt: \Delta\mathbf x
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 242
    role: notes_stage
    excerpt: \underbrace{O_{\rm orb}\Delta\mathbf x}_{\Delta\mathbf x_{\rm orbit}}
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 244
    role: notes_stage
    excerpt: \underbrace{Q_{\rm quot}\Delta\mathbf x}_{\Delta\mathbf x_{\rm fail}}.
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 249
    role: notes_stage
    excerpt: \boxed{\Delta\mathbf x_{\rm orbit}\in\ker M_*,}
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 251
    role: notes_stage
    excerpt: \boxed{M_*\Delta\mathbf x_{\rm fail}=M_*\Delta\mathbf x=\mathbf q.}
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 262
    role: notes_stage
    excerpt: \mathbf q=M_*\Delta\mathbf x=
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 271
    role: notes_stage
    excerpt: \boxed{\Delta\mathbf x_{\rm fail}=Q_{\rm quot}\Delta\mathbf x=S_{(T,K_\eta,\mu)}\mathbf q.}
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 275
    role: notes_stage
    excerpt: \boxed{(\Delta_T)_{\rm fail}=\frac{q_{\rm tr}}{1+\chi_{0,*}},}
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 278
    role: notes_stage
    excerpt: \boxed{(\Delta_{K_\eta})_{\rm fail}=-q_\eta,}
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 281
    role: notes_stage
    excerpt: \boxed{(\Delta_\mu)_{\rm fail}=\frac{F_*}{1+\chi_{0,*}}q_{\rm tr}+q_{\rm nt}-q_\eta.}
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 283
    role: notes_stage
    excerpt: All five free similarity directions vanish in `\(\Delta\mathbf x_{\rm fail}\)`.
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 288
    role: notes_stage
    excerpt: '> \((\Delta_T,\Delta_{K_\eta},\Delta_\mu)\).'
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 298
    role: notes_stage
    excerpt: \Delta\mathbf x_{\rm orbit}=O_{\rm orb}\Delta\mathbf x,
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 302
    role: notes_stage
    excerpt: (\Delta_\lambda,\Delta_c,\Delta_\gamma,\Delta_U,\Delta_W)_{\rm orbit}
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 304
    role: notes_stage
    excerpt: (\Delta_\lambda,\Delta_c,\Delta_\gamma,\Delta_U,\Delta_W).
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 310
    role: notes_stage
    excerpt: \boxed{\alpha_*:=\frac{1+\delta_{U,*}}{1+\chi_{0,*}}.}
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 314
    role: notes_stage
    excerpt: \boxed{(\Delta_{K_\eta})_{\rm orbit}=2\Delta_c-\Delta_U,}
  - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
    line: 317
    role: notes_stage
    excerpt: \boxed{(\Delta_T)_{\rm orbit}=\Delta_U-\alpha_*(\Delta_\gamma+\Delta_c-\Delta_U),}
  origin_claims:
  - parameter: Delta_Q
    introduced_at_stage: 195
    introduced_at_line: 185
    citation:
      path: notes/stages/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch.md
      line: 185
      excerpt: \boxed{\Delta_Q:=\chi_Q-1.}
  constraints:
  - parameter: Delta_Q
    constraint_kind: internal_consistency
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch.md
      line: 185
      excerpt: \boxed{\Delta_Q:=\chi_Q-1.}
  graph_context:
    contexts: []
    graph_gaps:
    - source: paper/stages/stage_192.tex
      attempted_sources:
      - paper/stages/stage_192.tex
      - research/pde_ledger/paper/stages/stage_192.tex
      graph_gap: true
      reason: no atlas node tied to this source path
  provenance_findings:
  - type: graph_gap
    severity: needs_triage
    summary: The atlas graph wrapper returned no source nodes for one or more primary sources.
    gaps:
    - source: paper/stages/stage_192.tex
      attempted_sources:
      - paper/stages/stage_192.tex
      - research/pde_ledger/paper/stages/stage_192.tex
      graph_gap: true
      reason: no atlas node tied to this source path
  - type: paper_card_overclaim
    severity: high
    summary: Delta_Q does NOT appear in the stage-192 per-stage notes (192 is the orbit/quotient projector calculus; it has
      Delta_orbit, Delta_T, Delta_{K_eta}, Delta_mu but no Delta_Q, no chi_Q). The Delta_Q=chi_Q-1 attribution to stage 192
      comes ONLY from the card's shared boilerplate (stage_192.tex line 13 'collapses Packet A to Delta_Q=chi_Q-1' and Checks
      line 20), an identical block reused across cards 192-195. Delta_Q is genuinely DEFINED at stage 195 line 185 (boxed
      Delta_Q:=chi_Q-1), built on the chi_Q fixed at stage 194. The binding notes citation cannot be met at stage 192; recorded
      at the true stage-195 origin.
    citations:
    - path: paper/stages/stage_192.tex
      line: 13
      excerpt: '...and collapses Packet A to \(\Delta_Q=\chi_Q-1\).}'
    - path: notes/stages/moving_throat_pde_stage192_orbit_quotient_projectors.md
      line: 1
      excerpt: '# Moving-Throat PDE — Stage 192: Exact Orbit/Quotient Projectors and the Microscopic Orbit-Lock Split'
  downstream_dependents:
  - '196'
  - '197'
  synthesis_ingested_at: '2026-06-11T19:18:45Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_21/fit_stage_192_delta_q__delta_q.yaml
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
- source: paper/stages/stage_192.tex
  attempted_sources:
  - paper/stages/stage_192.tex
  - research/pde_ledger/paper/stages/stage_192.tex
  graph_gap: true
  reason: no atlas node tied to this source path
```

Output a markdown adversarial report only when this is not a dry run. For dry runs, stop after confirming the prompt has enough context and record no verdict.
