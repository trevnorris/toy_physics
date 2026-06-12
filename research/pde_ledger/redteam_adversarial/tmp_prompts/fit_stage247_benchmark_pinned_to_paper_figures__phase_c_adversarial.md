# Phase C Adversarial Audit Wrapper

Candidate: `fit_stage247_benchmark_pinned_to_paper_figures`

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
id: fit_stage247_benchmark_pinned_to_paper_figures
candidate_key: stage247_benchmark_pinned_to_paper_figures
dry_run: false
dry_run_id: null
anchor_stages:
- '247'
file_line_citations:
- path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
  line: 239
  excerpt: '    # F5: pin independently derived benchmark quantities to the paper figures.'
  stage: '247'
parameter_names:
- Lvar_soft
- lambda_L
batch_id: null
status: provenance_built
status_timestamps:
  pending: '2026-06-11T08:21:51Z'
  scanned: '2026-06-11T08:21:51Z'
  provenance_built: '2026-06-11T19:37:43Z'
codex_session:
  by_parameter:
    Lvar_soft:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    lambda_L:
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
  - redteam_adversarial/provenance/fit_stage247_benchmark_pinned_to_paper_figures__lambda_l.yaml
  phase_c_prompt: null
verdict:
  adversarial: null
  adjudication: null
modality_attribution:
- claim_label
modality_fragments:
- candidate_key: stage247_benchmark_pinned_to_paper_figures
  modality: claim_label
  anchor_stage: '247'
  parameter_names:
  - Lvar_soft
  - lambda_L
  reason: Script comment asserts the benchmark quantities are "independently derived" while the assertions literally pin them
    to paper-stated decimals (e.g. Lvar = 20.01677473), a derived-status claim coupled to value matching.
  citation:
    path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
    line: 239
    excerpt: '    # F5: pin independently derived benchmark quantities to the paper figures.'
    stage: '247'
phase_b_status: synthesis_complete
```

## Primary Sources

- scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py

## Provenance Slice

```yaml
- schema_version: 1
  candidate_id: fit_stage247_benchmark_pinned_to_paper_figures
  parameter_name: lambda_L
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage247_benchmark_pinned_to_paper_figures__phase_b__lambda_l.md
  agent_synthesis_required: true
  dry_run_fixture: false
  fixture_note: null
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - '247'
  notes_source_opened: true
  source_evidence:
  - path: paper/stages/stage_247.tex
    line: 58
    role: paper_stage_tex
    excerpt: \mathfrak L(r):=\Lambda(r)\varrho(r),
  - path: paper/stages/stage_247.tex
    line: 60
    role: paper_stage_tex
    excerpt: S_{\rm leak}=\frac{8\sqrt2\eta_{\rm leak}\mu_w\rho_0}{\pi^{5/2}\lambda^3}\mathfrak L,
  - path: paper/stages/stage_247.tex
    line: 66
    role: paper_stage_tex
    excerpt: \frac{512\eta_{\rm leak}^{2}\mu_wq\rho_0}{\pi^4\lambda^2}\mathfrak L^2,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 191
    role: notes_stage
    excerpt: \mathfrak L(r):=\Lambda(r)\varrho(r).
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 198
    role: notes_stage
    excerpt: \frac{8\sqrt2\,\eta_{\rm leak}\mu_w\rho_0}{\pi^{5/2}\lambda^3}
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 206
    role: notes_stage
    excerpt: \frac{512\,\eta_{\rm leak}^2\mu_w q\rho_0}{\pi^4\lambda^2}
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 296
    role: notes_stage
    excerpt: \lambda_L S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 298
    role: notes_stage
    excerpt: \lambda_W \mathcal W_w^{\rm sess}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 308
    role: notes_stage
    excerpt: '- `\lambda_L\ge 0` converts the open-subsystem leakage source into stationary barrier units,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 309
    role: notes_stage
    excerpt: '- `\lambda_W\ge 0` converts the reduced Session-I work scalar into stationary barrier units,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 320
    role: notes_stage
    excerpt: \lambda_L S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 322
    role: notes_stage
    excerpt: \lambda_W \mathcal W_w^{\rm sess}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 438
    role: notes_stage
    excerpt: \frac{512\,\eta_{\rm leak}^2\mu_w q\rho_0}{\pi^4\lambda^2}
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 443
    role: notes_stage
    excerpt: \lambda=1,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 456
    role: notes_stage
    excerpt: \mathfrak L(r_{\rm soft})=\Lambda(r_{\rm soft})\varrho(r_{\rm soft})
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 491
    role: notes_stage
    excerpt: \lambda_W=1.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 495
    role: notes_stage
    excerpt: \lambda_L,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 500
    role: notes_stage
    excerpt: \lambda_L
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 518
    role: notes_stage
    excerpt: \lambda_L=0.26971918.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 571
    role: notes_stage
    excerpt: \lambda_L=0.26971918.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 604
    role: notes_stage
    excerpt: \lambda_L S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 606
    role: notes_stage
    excerpt: \lambda_W \mathcal W_w^{\rm sess}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 618
    role: notes_stage
    excerpt: \lambda_L S_{\rm leak}+\lambda_W\mathcal W_w^{\rm sess}+\Delta E_{UV}+\mathcal M_\sigma,
  - path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
    line: 13
    role: sympy_script
    excerpt: 6. and a Session-I one-point benchmark decomposition at r_soft = 0.18.
  - path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
    line: 106
    role: sympy_script
    excerpt: 'Lvar = sp.symbols("Lvar", positive=True, real=True)  # shorthand for Lambda(r)*varrho(r)'
  - path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
    line: 173
    role: sympy_script
    excerpt: lambda_L, lambda_W = sp.symbols("lambda_L lambda_W", nonnegative=True, real=True)
  - path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
    line: 175
    role: sympy_script
    excerpt: V_eff = V_short - lambda_L * S_leak - lambda_W * W_sess - E_UV - M_sigma
  - path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
    line: 177
    role: sympy_script
    excerpt: lowering_expected = lambda_L * S_leak + lambda_W * W_sess + E_UV + M_sigma
  - path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
    line: 184
    role: sympy_script
    excerpt: '# 6. Session-I one-point benchmark at r_soft = 0.18.'
  - path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
    line: 186
    role: sympy_script
    excerpt: section("6. Session-I one-point benchmark")
  - path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
    line: 232
    role: sympy_script
    excerpt: '# F2: the inverted Lvar reproduces the recorded benchmark work scalar, and'
  - path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
    line: 233
    role: sympy_script
    excerpt: '#     matches the paper-stated Lvar(r_soft) = 20.01677473.'
  - path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
    line: 238
    role: sympy_script
    excerpt: lambda_L_soft = sp.N((Vshort_num - Wsess_obs - UVdrop_obs - M_sigma_num - Veff_obs) / S_soft, 16)
  - path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
    line: 239
    role: sympy_script
    excerpt: '# F5: pin independently derived benchmark quantities to the paper figures.'
  - path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
    line: 243
    role: sympy_script
    excerpt: assert abs(float(lambda_L_soft) - 0.26971918) < 1e-6
  - path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
    line: 244
    role: sympy_script
    excerpt: lambda_L_paper = sp.Float("0.26971918")
  - path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
    line: 245
    role: sympy_script
    excerpt: '# F4: lowered potential is below the baseline on the benchmark slice (lowering theorem).'
  - path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
    line: 246
    role: sympy_script
    excerpt: Veff_session = sp.N(Vshort_num - lambda_L_paper * S_soft - Wsess_obs - UVdrop_obs - M_sigma_num, 16)
  - path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
    line: 248
    role: sympy_script
    excerpt: '# F5: forward benchmark decomposition with the paper''s lambda_L (falsifiable closure).'
  - path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
    line: 249
    role: sympy_script
    excerpt: Veff_forward = sp.N(Vshort_num - lambda_L_paper * S_soft - Wsess_obs - UVdrop_obs - M_sigma_num, 16)
  origin_claims:
  - parameter: lambda_L
    introduced_at_stage: 247
    introduced_at_line: 518
    citation:
      path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
      line: 518
      excerpt: \lambda_L=0.26971918.
  constraints:
  - parameter: lambda_L
    constraint_kind: published_target
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
      line: 497
      excerpt: 'which is fixed exactly by the recorded softening point:'
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings:
  - type: derive_vs_fit_mismatch
    severity: needs_triage
    summary: This candidate framing ('benchmark pinned to paper figures') describes lambda_L=0.26971918 as 'independently
      derived benchmark quantities pinned to the paper figures' in the audit script (line 239), but lambda_L is itself back-solved
      from the recorded target V_eff^sess(r_soft)=1.74701126; the 'pin' asserts the back-solved value matches the recorded
      benchmark, which is the fit, not an independent derivation. The notes are honest about this (line 497), so the mismatch
      is between the script comment's 'independently derived' wording and the actual back-solve, not between notes and mechanism.
    citations:
    - path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
      line: 239
      excerpt: '# F5: pin independently derived benchmark quantities to the paper figures.'
    - path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
      line: 238
      excerpt: lambda_L_soft = sp.N((Vshort_num - Wsess_obs - UVdrop_obs - M_sigma_num - Veff_obs) / S_soft, 16)
    - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
      line: 497
      excerpt: 'which is fixed exactly by the recorded softening point:'
  downstream_dependents:
  - '248'
  synthesis_ingested_at: '2026-06-11T19:37:43Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_27/fit_stage247_benchmark_pinned_to_paper_figures__lambda_l.yaml
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
