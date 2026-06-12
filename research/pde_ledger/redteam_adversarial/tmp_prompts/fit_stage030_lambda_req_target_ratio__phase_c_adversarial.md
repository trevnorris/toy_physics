# Phase C Adversarial Audit Wrapper

Candidate: `fit_stage030_lambda_req_target_ratio`

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
id: fit_stage030_lambda_req_target_ratio
candidate_key: stage030_lambda_req_target_ratio
dry_run: false
dry_run_id: null
anchor_stages:
- '030'
file_line_citations:
- path: scripts/moving_throat_pde_stage030_selected_mode_normalization_sympy_audit.py
  line: 169
  excerpt: lambda_req = sp.simplify(mhat**2 * beta0 * s_minus_closed / NQ_target)
  stage: '030'
parameter_names:
- lambda_req
batch_id: null
status: provenance_built
status_timestamps:
  pending: '2026-06-11T08:21:51Z'
  scanned: '2026-06-11T08:21:51Z'
  provenance_built: '2026-06-11T18:25:37Z'
codex_session:
  by_parameter:
    lambda_req:
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
  - redteam_adversarial/provenance/fit_stage030_lambda_req_target_ratio__lambda_req.yaml
  phase_c_prompt: null
verdict:
  adversarial: null
  adjudication: null
modality_attribution:
- existing_provenance
modality_fragments:
- candidate_key: stage030_lambda_req_target_ratio
  modality: existing_provenance
  anchor_stage: '030'
  parameter_names:
  - lambda_req
  reason: The required selected eigenvalue is defined as a ratio against NQ_target and only printed, never asserted — the
    script's own comment (lines 171-175) notes the gate identity collapses to 0 by construction — so the value is posited-to-meet-target
    with no independent in-stage check.
  citation:
    path: scripts/moving_throat_pde_stage030_selected_mode_normalization_sympy_audit.py
    line: 169
    excerpt: lambda_req = sp.simplify(mhat**2 * beta0 * s_minus_closed / NQ_target)
    stage: '030'
phase_b_status: synthesis_complete
```

## Primary Sources

- scripts/moving_throat_pde_stage030_selected_mode_normalization_sympy_audit.py

## Provenance Slice

```yaml
- schema_version: 1
  candidate_id: fit_stage030_lambda_req_target_ratio
  parameter_name: lambda_req
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage030_lambda_req_target_ratio__phase_b__lambda_req.md
  agent_synthesis_required: true
  dry_run_fixture: false
  fixture_note: null
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - '030'
  notes_source_opened: true
  source_evidence:
  - path: paper/stages/stage_030.tex
    line: 9
    role: paper_stage_tex
    excerpt: \stagefield{Inputs}{The inputs are the static selected matrix from Stage~029, the lower eigenvalue \(\lambda_-\),
      the selected eigenvector \(e_-\), and the overlap \(s_-=(v\cdot e_-)^2\).}
  - path: paper/stages/stage_030.tex
    line: 33
    role: paper_stage_tex
    excerpt: \label{eq:app-stage030-lambda-minus}
  - path: paper/stages/stage_030.tex
    line: 35
    role: paper_stage_tex
    excerpt: \lambda_-=\frac{A+B-\alpha_0\sigma-
  - path: paper/stages/stage_030.tex
    line: 40
    role: paper_stage_tex
    excerpt: R_\lambda=\sqrt{(\Delta K_{\rm ax}+\alpha_0\delta_\kappa)^2+4\alpha_0^2\Pi_\kappa}.
  - path: paper/stages/stage_030.tex
    line: 48
    role: paper_stage_tex
    excerpt: \sigma+\frac{(\Delta K_{\rm ax}+\alpha_0\delta_\kappa)\delta_\kappa+4\alpha_0\Pi_\kappa}{R_\lambda}
  - path: paper/stages/stage_030.tex
    line: 55
    role: paper_stage_tex
    excerpt: s_-=-\frac{d\lambda_-}{d\alpha_0}}.
  - path: paper/stages/stage_030.tex
    line: 63
    role: paper_stage_tex
    excerpt: D_{-0}=\lambda_-.
  - path: paper/stages/stage_030.tex
    line: 82
    role: paper_stage_tex
    excerpt: \Gamma_{5,-}=\frac{a^5}{27c_s^5}\frac{\beta_0s_-}{\lambda_-}.
  - path: paper/stages/stage_030.tex
    line: 88
    role: paper_stage_tex
    excerpt: P_{0,-}=\frac{\beta_0s_-}{\lambda_-}
  - path: paper/stages/stage_030.tex
    line: 89
    role: paper_stage_tex
    excerpt: =-\beta_0\frac{d\ln\lambda_-}{d\alpha_0}}.
  - path: paper/stages/stage_030.tex
    line: 98
    role: paper_stage_tex
    excerpt: \stagefield{Output}{Stage~030 outputs the lower eigenvalue \eqref{eq:app-stage030-lambda-minus}, the selected
      overlap \eqref{eq:app-stage030-s-minus}, the Hellmann--Feynman identity \eqref{eq:app-stage030-HF}, the selected prefactor
      \eqref{eq:app-stage030-P0minus}, and the target \eqref{eq:app-stage030-selected-target}.}
  - path: paper/stages/stage_030.tex
    line: 101
    role: paper_stage_tex
    excerpt: \item Differentiating the eigenvalue with respect to \(\alpha_0\) gives \(d\lambda_-/d\alpha_0=-e_-^Tvv^Te_-=-s_-\).
  - path: paper/stages/stage_030.tex
    line: 102
    role: paper_stage_tex
    excerpt: \item Normalizing \(D_-/D_{-0}\) gives the response coefficients as quotients by \(\lambda_-\).
  - path: paper/stages/stage_030.tex
    line: 103
    role: paper_stage_tex
    excerpt: \item The outgoing coefficient contains only \(\beta_0s_-/\lambda_-\), matching the selected projection of Stage~029.
  - path: paper/stages/stage_030.tex
    line: 106
    role: paper_stage_tex
    excerpt: \stagefield{Downstream use}{Stage~031 proves that \(P_{0,-}\) is monotone on the stable branch.  Stage~032 replaces
      the abstract \(\widehat m_-\) by the natural D/N source-map factor.  Stage~034 uses \(\lambda_-\) and \(s_-\) to build
      the softening-depth normal form.}
  - path: notes/stages/moving_throat_pde_stage030_selected_mode_normalization.md
    line: 18
    role: notes_stage
    excerpt: If the conservative selected lower wall eigenvalue is written as `lambda_-(omega)`, then the normalized selected-mode
      response has the low-frequency form
  - path: notes/stages/moving_throat_pde_stage030_selected_mode_normalization.md
    line: 28
    role: notes_stage
    excerpt: '`P_{0,-} = beta_0 (v.e_-)^2 / lambda_-(0)`'
  - path: notes/stages/moving_throat_pde_stage030_selected_mode_normalization.md
    line: 32
    role: notes_stage
    excerpt: '`P_{0,-} = - beta_0 d(ln lambda_-)/d alpha_0`.'
  - path: notes/stages/moving_throat_pde_stage030_selected_mode_normalization.md
    line: 61
    role: notes_stage
    excerpt: '`lambda_-`'
  - path: notes/stages/moving_throat_pde_stage030_selected_mode_normalization.md
    line: 70
    role: notes_stage
    excerpt: '`lambda_+ = ( A + B - alpha_0 sigma + R ) / 2`.'
  - path: notes/stages/moving_throat_pde_stage030_selected_mode_normalization.md
    line: 74
    role: notes_stage
    excerpt: '`D_{-0} = lambda_-`'
  - path: notes/stages/moving_throat_pde_stage030_selected_mode_normalization.md
    line: 76
    role: notes_stage
    excerpt: is the selected conservative wall stiffness at zero frequency.
  - path: notes/stages/moving_throat_pde_stage030_selected_mode_normalization.md
    line: 84
    role: notes_stage
    excerpt: '`(v.e_-)^2 = - d lambda_- / d alpha_0`.'
  - path: notes/stages/moving_throat_pde_stage030_selected_mode_normalization.md
    line: 109
    role: notes_stage
    excerpt: Write the selected wall operator in the same low-frequency form used in Stage 022,
  - path: notes/stages/moving_throat_pde_stage030_selected_mode_normalization.md
    line: 132
    role: notes_stage
    excerpt: '`           = beta_5 s_- / lambda_-`.'
  - path: notes/stages/moving_throat_pde_stage030_selected_mode_normalization.md
    line: 152
    role: notes_stage
    excerpt: '`P_{0,-} = beta_0 s_- / lambda_-`.'
  - path: notes/stages/moving_throat_pde_stage030_selected_mode_normalization.md
    line: 158
    role: notes_stage
    excerpt: '`P_{0,-} = - beta_0 d(ln lambda_-)/d alpha_0`.'
  - path: notes/stages/moving_throat_pde_stage030_selected_mode_normalization.md
    line: 183
    role: notes_stage
    excerpt: So the selected branch is required to hit exactly the same normalization stack as the isotropic Stage-022/025
      branch, but now with the selected-mode prefactor in place of the old isotropic lane prefactor.
  - path: notes/stages/moving_throat_pde_stage030_selected_mode_normalization.md
    line: 187
    role: notes_stage
    excerpt: '`lambda_- = mhat_-^2 beta_0 s_- / N_Q^(target)`,'
  - path: notes/stages/moving_throat_pde_stage030_selected_mode_normalization.md
    line: 199
    role: notes_stage
    excerpt: '`lambda_- = beta_0 s_- * 5 a^5 c^5 / (54 G c_s^5)`.'
  - path: notes/stages/moving_throat_pde_stage030_selected_mode_normalization.md
    line: 211
    role: notes_stage
    excerpt: '`lambda_- lambda_+ = A B - alpha_0 ( B kappa_0^2 + A kappa_1^2 )`.'
  - path: notes/stages/moving_throat_pde_stage030_selected_mode_normalization.md
    line: 220
    role: notes_stage
    excerpt: If `lambda_-` is too large, the selected outgoing coefficient is too small.
  - path: notes/stages/moving_throat_pde_stage030_selected_mode_normalization.md
    line: 221
    role: notes_stage
    excerpt: If `lambda_-` is too small, the branch is close to instability.
  - path: notes/stages/moving_throat_pde_stage030_selected_mode_normalization.md
    line: 225
    role: notes_stage
    excerpt: '> what conservative selected eigenvalue `lambda_-` and selected overlap `s_-` does the physical moving-throat
      branch produce, and does their ratio land on the universal target before the wall mode softens?'
  - path: notes/stages/moving_throat_pde_stage030_selected_mode_normalization.md
    line: 236
    role: notes_stage
    excerpt: '`P_{0,-} = beta_0 (v.e_-)^2 / lambda_-`'
  - path: notes/stages/moving_throat_pde_stage030_selected_mode_normalization.md
    line: 240
    role: notes_stage
    excerpt: '`P_{0,-} = - beta_0 d(ln lambda_-)/d alpha_0`.'
  - path: scripts/moving_throat_pde_stage030_selected_mode_normalization_sympy_audit.py
    line: 9
    role: sympy_script
    excerpt: 3. Exact Hellmann–Feynman selected overlap s_- = -(d lambda_- / d alpha).
  - path: scripts/moving_throat_pde_stage030_selected_mode_normalization_sympy_audit.py
    line: 97
    role: sympy_script
    excerpt: print("lambda_- =")
  - path: scripts/moving_throat_pde_stage030_selected_mode_normalization_sympy_audit.py
    line: 99
    role: sympy_script
    excerpt: print("lambda_+ =")
  - path: scripts/moving_throat_pde_stage030_selected_mode_normalization_sympy_audit.py
    line: 141
    role: sympy_script
    excerpt: '"P0_sel + beta0*d(log lambda_-)/d alpha",'
  origin_claims:
  - parameter: lambda_req
    introduced_at_stage: 30
    introduced_at_line: 187
    citation:
      path: notes/stages/moving_throat_pde_stage030_selected_mode_normalization.md
      line: 187
      excerpt: '`lambda_- = mhat_-^2 beta_0 s_- / N_Q^(target)`,'
  constraints:
  - parameter: lambda_req
    constraint_kind: published_target
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage030_selected_mode_normalization.md
      line: 185
      excerpt: 'Equivalently, the target becomes a direct conservative spectral condition:'
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings:
  - type: derive_vs_fit_mismatch
    severity: needs_triage
    summary: 'lambda_req is the conservative selected eigenvalue REQUIRED so the selected-mode prefactor hits the external
      GR target: lambda_- = mhat_-^2 beta_0 s_- / N_Q^target with N_Q^target = 54/5 G c_s^5/(5 a^5 c^5) (the external Einstein-quadrupole
      number imported at stage 022). At leading order lambda_- = beta_0 s_- * 5 a^5 c^5/(54 G c_s^5). So the ''required eigenvalue''
      is back-solved from the published GR target, not produced by the spectral problem; the notes frame it as the value the
      physical branch must REACH (''does their ratio land on the universal target''). Classification published_target. The
      card flags this as StatusOpen (whether the physical branch can reach it), consistent with a target rather than a derived
      output.'
    citations:
    - path: notes/stages/moving_throat_pde_stage030_selected_mode_normalization.md
      line: 191
      excerpt: '`N_Q^(target) = 54 G c_s^5 / (5 a^5 c^5)`.'
    - path: notes/stages/moving_throat_pde_stage030_selected_mode_normalization.md
      line: 225
      excerpt: '> what conservative selected eigenvalue `lambda_-` and selected overlap `s_-` does the physical moving-throat
        branch produce, and does their ratio land on the universal target before the wall mode softens?'
  downstream_dependents:
  - '031'
  - '033'
  synthesis_ingested_at: '2026-06-11T18:25:37Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_03/fit_stage030_lambda_req_target_ratio__lambda_req.yaml
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
