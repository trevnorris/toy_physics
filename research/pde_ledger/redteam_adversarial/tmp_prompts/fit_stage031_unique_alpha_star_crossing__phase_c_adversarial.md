# Phase C Adversarial Audit Wrapper

Candidate: `fit_stage031_unique_alpha_star_crossing`

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
id: fit_stage031_unique_alpha_star_crossing
candidate_key: stage031_unique_alpha_star_crossing
dry_run: false
dry_run_id: null
anchor_stages:
- '031'
file_line_citations:
- path: paper/stages/stage_031.tex
  line: 60
  excerpt: has a unique \(\alpha_*\in(0,\alpha_{\rm crit})\) such that \(P_{0,-}(\alpha_*)=P_{\rm target}\).  For the quadrupole
    target,
  stage: '031'
parameter_names:
- P_target
- alpha_star
batch_id: null
status: provenance_built
status_timestamps:
  pending: '2026-06-11T08:21:51Z'
  scanned: '2026-06-11T08:21:51Z'
  provenance_built: '2026-06-11T18:25:37Z'
codex_session:
  by_parameter:
    P_target:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    alpha_star:
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
  - redteam_adversarial/provenance/fit_stage031_unique_alpha_star_crossing__alpha_star.yaml
  phase_c_prompt: null
verdict:
  adversarial: null
  adjudication: null
modality_attribution:
- claim_label
modality_fragments:
- candidate_key: stage031_unique_alpha_star_crossing
  modality: claim_label
  anchor_stage: '031'
  parameter_names:
  - alpha_star
  - P_target
  reason: The "unique stable crossing theorem" guarantees an alpha_* exists matching any admissible target; uniqueness of
    the crossing is presented as a derivation result while it simultaneously certifies that alpha can always be tuned to the
    GR target.
  citation:
    path: paper/stages/stage_031.tex
    line: 60
    excerpt: has a unique \(\alpha_*\in(0,\alpha_{\rm crit})\) such that \(P_{0,-}(\alpha_*)=P_{\rm target}\).  For the quadrupole
      target,
    stage: '031'
phase_b_status: synthesis_complete
```

## Primary Sources

- paper/stages/stage_031.tex

## Provenance Slice

```yaml
- schema_version: 1
  candidate_id: fit_stage031_unique_alpha_star_crossing
  parameter_name: alpha_star
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage031_unique_alpha_star_crossing__phase_b__alpha_star.md
  agent_synthesis_required: true
  dry_run_fixture: false
  fixture_note: null
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - '031'
  notes_source_opened: true
  source_evidence:
  - path: paper/stages/stage_031.tex
    line: 9
    role: paper_stage_tex
    excerpt: \stagefield{Inputs}{The stage uses \(\lambda_-(\alpha_0)\), \(s_-(\alpha_0)\), \(P_{0,-}=\beta_0s_-/\lambda_-\),
      and the stable interval \(0\leq\alpha_0<\alpha_{\rm crit}\).}
  - path: paper/stages/stage_031.tex
    line: 18
    role: paper_stage_tex
    excerpt: \label{eq:app-stage031-ds-dalpha}
  - path: paper/stages/stage_031.tex
    line: 20
    role: paper_stage_tex
    excerpt: \frac{ds_-}{d\alpha_0}=
  - path: paper/stages/stage_031.tex
    line: 23
    role: paper_stage_tex
    excerpt: Using \(d\lambda_-/d\alpha_0=-s_-\), one obtains
  - path: paper/stages/stage_031.tex
    line: 25
    role: paper_stage_tex
    excerpt: \label{eq:app-stage031-dP-dalpha}
  - path: paper/stages/stage_031.tex
    line: 27
    role: paper_stage_tex
    excerpt: \frac{dP_{0,-}}{d\alpha_0}
  - path: paper/stages/stage_031.tex
    line: 28
    role: paper_stage_tex
    excerpt: =\beta_0\frac{(ds_-/d\alpha_0)\lambda_-+s_-^2}{\lambda_-^2}>0
  - path: paper/stages/stage_031.tex
    line: 41
    role: paper_stage_tex
    excerpt: \label{eq:app-stage031-alpha-crit}
  - path: paper/stages/stage_031.tex
    line: 43
    role: paper_stage_tex
    excerpt: \alpha_{\rm crit}=\frac{AB}{B\kappa_0^2+A\kappa_1^2}}.
  - path: paper/stages/stage_031.tex
    line: 45
    role: paper_stage_tex
    excerpt: As \(\alpha_0\to\alpha_{\rm crit}^-\), \(\lambda_-\to0^+\) while \(s_-\) remains positive, so
  - path: paper/stages/stage_031.tex
    line: 60
    role: paper_stage_tex
    excerpt: has a unique \(\alpha_*\in(0,\alpha_{\rm crit})\) such that \(P_{0,-}(\alpha_*)=P_{\rm target}\).  For the quadrupole
      target,
  - path: paper/stages/stage_031.tex
    line: 65
    role: paper_stage_tex
    excerpt: \stagefield{Output}{Stage~031 outputs the derivative identities \eqref{eq:app-stage031-ds-dalpha}--\eqref{eq:app-stage031-dP-dalpha},
      the refined threshold \eqref{eq:app-stage031-alpha-crit}, and the unique stable crossing theorem.}
  - path: paper/stages/stage_031.tex
    line: 68
    role: paper_stage_tex
    excerpt: \item The positivity of \(ds_-/d\alpha_0\) follows from \((\Delta K_{\rm ax})^2\Pi_\kappa>0\) and \(R_\lambda>0\).
  - path: notes/stages/moving_throat_pde_stage031_selected_branch_reachability.md
    line: 15
    role: notes_stage
    excerpt: It starts at the flat-branch value, diverges at the softening threshold, and therefore hits any target above
      its starting value exactly once before instability.
  - path: notes/stages/moving_throat_pde_stage031_selected_branch_reachability.md
    line: 32
    role: notes_stage
    excerpt: '`lambda_- = ( A + B - alpha_0 sigma - R ) / 2`,'
  - path: notes/stages/moving_throat_pde_stage031_selected_branch_reachability.md
    line: 34
    role: notes_stage
    excerpt: '`R = sqrt( (Delta K_ax + alpha_0 delta_kappa)^2 + 4 alpha_0^2 KappaProd )`.'
  - path: notes/stages/moving_throat_pde_stage031_selected_branch_reachability.md
    line: 38
    role: notes_stage
    excerpt: '`d s_- / d alpha_0 = 2 Delta K_ax^2 KappaProd / R^3`.'
  - path: notes/stages/moving_throat_pde_stage031_selected_branch_reachability.md
    line: 63
    role: notes_stage
    excerpt: '`d lambda_- / d alpha_0 = - s_-`,'
  - path: notes/stages/moving_throat_pde_stage031_selected_branch_reachability.md
    line: 67
    role: notes_stage
    excerpt: '`d P_{0,-} / d alpha_0`'
  - path: notes/stages/moving_throat_pde_stage031_selected_branch_reachability.md
    line: 68
    role: notes_stage
    excerpt: '`= beta_0 [ (d s_- / d alpha_0) lambda_- + s_-^2 ] / lambda_-^2`.'
  - path: notes/stages/moving_throat_pde_stage031_selected_branch_reachability.md
    line: 74
    role: notes_stage
    excerpt: '- `d s_- / d alpha_0 > 0`,'
  - path: notes/stages/moving_throat_pde_stage031_selected_branch_reachability.md
    line: 79
    role: notes_stage
    excerpt: '`d P_{0,-} / d alpha_0 > 0`.'
  - path: notes/stages/moving_throat_pde_stage031_selected_branch_reachability.md
    line: 88
    role: notes_stage
    excerpt: '## 3. Starting value at zero loading'
  - path: notes/stages/moving_throat_pde_stage031_selected_branch_reachability.md
    line: 90
    role: notes_stage
    excerpt: At `alpha_0 = 0`, the lower wall mode is just the flat `u_0` branch, so
  - path: notes/stages/moving_throat_pde_stage031_selected_branch_reachability.md
    line: 96
    role: notes_stage
    excerpt: Therefore the selected prefactor starts at
  - path: notes/stages/moving_throat_pde_stage031_selected_branch_reachability.md
    line: 100
    role: notes_stage
    excerpt: This is the exact stable-side starting value against which the universal target must be compared.
  - path: notes/stages/moving_throat_pde_stage031_selected_branch_reachability.md
    line: 108
    role: notes_stage
    excerpt: '`lambda_- lambda_+ = A B - alpha_0 ( B kappa_0^2 + A kappa_1^2 )`.'
  - path: notes/stages/moving_throat_pde_stage031_selected_branch_reachability.md
    line: 112
    role: notes_stage
    excerpt: '`alpha_crit = A B / ( B kappa_0^2 + A kappa_1^2 )`.'
  - path: notes/stages/moving_throat_pde_stage031_selected_branch_reachability.md
    line: 116
    role: notes_stage
    excerpt: '`lambda_-(alpha_crit) = 0`,'
  - path: notes/stages/moving_throat_pde_stage031_selected_branch_reachability.md
    line: 128
    role: notes_stage
    excerpt: as `alpha_0 -> alpha_crit^-` from the stable side.
  - path: notes/stages/moving_throat_pde_stage031_selected_branch_reachability.md
    line: 140
    role: notes_stage
    excerpt: '- `P_{0,-}(alpha_0)` is continuous on `0 <= alpha_0 < alpha_crit`,'
  - path: notes/stages/moving_throat_pde_stage031_selected_branch_reachability.md
    line: 142
    role: notes_stage
    excerpt: '- starts at `P_{0,-}(0)`,'
  - path: notes/stages/moving_throat_pde_stage031_selected_branch_reachability.md
    line: 143
    role: notes_stage
    excerpt: '- and diverges at `alpha_crit^-`.'
  - path: notes/stages/moving_throat_pde_stage031_selected_branch_reachability.md
    line: 147
    role: notes_stage
    excerpt: '> for every target value `P_target > P_{0,-}(0)`, there exists a **unique** stable-side loading `alpha_*` with
      `0 < alpha_* < alpha_crit` such that `P_{0,-}(alpha_*) = P_target`.'
  - path: notes/stages/moving_throat_pde_stage031_selected_branch_reachability.md
    line: 176
    role: notes_stage
    excerpt: 1. compute the physical stable-branch data `(Xi_0, beta_0, alpha_0)` from the coupled moving-throat operator,
  - path: notes/stages/moving_throat_pde_stage031_selected_branch_reachability.md
    line: 178
    role: notes_stage
    excerpt: '`P_{0,-}(alpha_0) = beta_0 (v.e_-)^2 / lambda_-`,'
  - path: notes/stages/moving_throat_pde_stage031_selected_branch_reachability.md
    line: 180
    role: notes_stage
    excerpt: 4. and check whether the resulting stable-side crossing sits on the natural passive/outgoing branch with `alpha_0
      < alpha_crit`.
  - path: scripts/moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.py
    line: 6
    role: sympy_script
    excerpt: 1. Exact derivative of the selected overlap s_-(alpha).
  - path: scripts/moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.py
    line: 7
    role: sympy_script
    excerpt: 2. Exact derivative formula for the selected prefactor P0_-(alpha).
  - path: scripts/moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.py
    line: 8
    role: sympy_script
    excerpt: 3. Initial-value formulas at alpha = 0.
  origin_claims:
  - parameter: alpha_star
    introduced_at_stage: 31
    introduced_at_line: 147
    citation:
      path: notes/stages/moving_throat_pde_stage031_selected_branch_reachability.md
      line: 147
      excerpt: '> for every target value `P_target > P_{0,-}(0)`, there exists a **unique** stable-side loading `alpha_*`
        with `0 < alpha_* < alpha_crit` such that `P_{0,-}(alpha_*) = P_target`.'
  constraints:
  - parameter: alpha_star
    constraint_kind: published_target
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage031_selected_branch_reachability.md
      line: 151
      excerpt: '`P_target = N_Q^(target) / mhat_-^2`'
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings:
  - type: derive_vs_fit_mismatch
    severity: needs_triage
    summary: 'alpha_star is the unique stable-side loading defined by P_{0,-}(alpha_*) = P_target, where P_target = N_Q^target/mhat_-^2
      and N_Q^target = 54/5 G c_s^5/(5 a^5 c^5) is the external GR quadrupole number. So alpha_star''s VALUE is whatever loading
      makes the selected prefactor hit the external target: it is back-solved from the published GR number via the monotone
      crossing. The stage''s content (existence + uniqueness of the crossing) is an internal_consistency theorem about the
      monotone P_{0,-}(alpha_0), but the alpha_star VALUE the candidate names is fixed by the external target -> published_target.
      Note the value is not numerically pinned in-stage (depends on the physical branch data Xi_0, beta_0, A, B not yet computed);
      it is the target-matching root.'
    citations:
    - path: notes/stages/moving_throat_pde_stage031_selected_branch_reachability.md
      line: 155
      excerpt: '`N_Q^(target) = 54 G c_s^5 / (5 a^5 c^5)`,'
    - path: notes/stages/moving_throat_pde_stage031_selected_branch_reachability.md
      line: 160
      excerpt: The only remaining question is whether the physical moving-throat branch puts the system in the right region
        of parameter space so that the crossing occurs on the natural passive/outgoing branch before other approximations
        fail.
  downstream_dependents:
  - '033'
  synthesis_ingested_at: '2026-06-11T18:25:37Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_03/fit_stage031_unique_alpha_star_crossing__alpha_star.yaml
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
