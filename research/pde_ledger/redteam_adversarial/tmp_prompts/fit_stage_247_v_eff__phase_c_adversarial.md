# Phase C Adversarial Audit Wrapper

Candidate: `fit_stage_247_v_eff`

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
id: fit_stage_247_v_eff
candidate_key: stage_247_v_eff
dry_run: false
dry_run_id: null
anchor_stages:
- '247'
file_line_citations:
- path: redteam/pass2/reports/stage_247.md
  line: 51
  role: pass2_stage_report
  stage: '247'
  excerpt: '| (5) inequality `V_eff<=V_short` | py:247 (benchmark slice) / wl (M8 numeric forward) | match (numeric on branch)
    |'
parameter_names:
- V_eff
- V_short
batch_id: null
status: provenance_built
status_timestamps:
  pending: '2026-06-11T08:11:49Z'
  scanned: '2026-06-11T08:11:49Z'
  provenance_built: '2026-06-11T19:37:44Z'
codex_session:
  by_parameter:
    V_eff:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    V_short:
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
  - redteam_adversarial/provenance/fit_stage_247_v_eff__v_eff.yaml
  phase_c_prompt: null
verdict:
  adversarial: null
  adjudication: null
modality_attribution:
- existing_provenance
modality_fragments:
- candidate_key: stage_247_v_eff
  modality: existing_provenance
  anchor_stage: '247'
  parameter_names:
  - V_eff
  - V_short
  reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
  citation:
    path: redteam/pass2/reports/stage_247.md
    line: 51
    role: pass2_stage_report
    stage: '247'
    excerpt: '| (5) inequality `V_eff<=V_short` | py:247 (benchmark slice) / wl (M8 numeric forward) | match (numeric on branch)
      |'
phase_b_status: synthesis_complete
```

## Primary Sources

- redteam/pass2/reports/stage_247.md

## Provenance Slice

```yaml
- schema_version: 1
  candidate_id: fit_stage_247_v_eff
  parameter_name: V_eff
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage_247_v_eff__phase_b__v_eff.md
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
    line: 20
    role: paper_stage_tex
    excerpt: With primitive short-range source amplitudes \(\beta_Q,\beta_U,\beta_W\), the product-family coefficients are
  - path: paper/stages/stage_247.tex
    line: 76
    role: paper_stage_tex
    excerpt: The stationary compiler is \cref{eq:app-part08-intro-veff230}; the lowering identity is \cref{eq:app-part08-main-lowering-identity}.
  - path: paper/stages/stage_247.tex
    line: 93
    role: paper_stage_tex
    excerpt: This shows that, once the one-port baseline is fixed, the recorded stationary lowering can be decomposed into
      work, \(U/V\) drain, compensated source response, and one remaining leakage-to-barrier embedding coefficient.  It does
      not change the short-range kernel class.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 83
    role: notes_stage
    excerpt: \frac{V_{\rm eff}(r_{\rm soft})}{V_{\rm Coul}(r_{\rm soft})}\approx 0.31446203,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 131
    role: notes_stage
    excerpt: the exact one-port product-family coefficients are
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 292
    role: notes_stage
    excerpt: V_{\rm eff}^{(247)}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 318
    role: notes_stage
    excerpt: V_{\rm short}^{(1p)}(r)-V_{\rm eff}^{(247)}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 333
    role: notes_stage
    excerpt: V_{\rm eff}^{(247)}(r)\le V_{\rm short}^{(1p)}(r).
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 372
    role: notes_stage
    excerpt: V_{\rm eff}(r_{\rm soft})=1.74701126,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 378
    role: notes_stage
    excerpt: \frac{V_{\rm eff}(r_{\rm soft})}{V_{\rm Coul}(r_{\rm soft})}=0.31446203.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 493
    role: notes_stage
    excerpt: Then the Stage-247 compiler requires only one additional leakage embedding coefficient,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 510
    role: notes_stage
    excerpt: V_{\rm eff}^{\rm sess}(r_{\rm soft})}
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 545
    role: notes_stage
    excerpt: '> 5. and one remaining leakage-to-barrier embedding coefficient.'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 586
    role: notes_stage
    excerpt: 4. It reduces the remaining stationary embedding freedom, on the benchmark slice, to a single leakage-to-barrier
      coefficient once the work channel is taken in the Session-I normalization.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 600
    role: notes_stage
    excerpt: V_{\rm eff}^{(247)}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 616
    role: notes_stage
    excerpt: V_{\rm short}^{(1p)}-V_{\rm eff}^{(247)}
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 638
    role: notes_stage
    excerpt: 1. take `V_{\rm eff}^{(247)}(r)` as the lowered stationary front end,
  - path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
    line: 175
    role: sympy_script
    excerpt: V_eff = V_short - lambda_L * S_leak - lambda_W * W_sess - E_UV - M_sigma
  - path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
    line: 176
    role: sympy_script
    excerpt: lowering_gap = sp.expand(V_short - V_eff)
  - path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
    line: 213
    role: sympy_script
    excerpt: Veff_obs = sp.Float("1.74701126")
  - path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
    line: 238
    role: sympy_script
    excerpt: lambda_L_soft = sp.N((Vshort_num - Wsess_obs - UVdrop_obs - M_sigma_num - Veff_obs) / S_soft, 16)
  - path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
    line: 246
    role: sympy_script
    excerpt: Veff_session = sp.N(Vshort_num - lambda_L_paper * S_soft - Wsess_obs - UVdrop_obs - M_sigma_num, 16)
  - path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
    line: 247
    role: sympy_script
    excerpt: assert float(Vshort_num - Veff_session) >= 0
  - path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
    line: 249
    role: sympy_script
    excerpt: Veff_forward = sp.N(Vshort_num - lambda_L_paper * S_soft - Wsess_obs - UVdrop_obs - M_sigma_num, 16)
  - path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
    line: 250
    role: sympy_script
    excerpt: assert abs(float(Veff_forward) - float(Veff_obs)) < 1e-6
  - path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
    line: 264
    role: sympy_script
    excerpt: print("V_eff rebuilt         =", Vrebuild_soft)
  origin_claims:
  - parameter: V_eff
    introduced_at_stage: 247
    introduced_at_line: 368
    citation:
      path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
      line: 368
      excerpt: r_{\rm soft}=0.18,
  constraints:
  - parameter: V_eff
    constraint_kind: published_target
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
      line: 367
      excerpt: Appendix B.2 of the barrier-session write-up records the strongest stationary softening point at
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings:
  - type: constraint_kind_ambiguous
    severity: low
    summary: 'Two distinct senses of V_eff: (a) the compiler FORMULA V_eff^(247)(r) = V_short - lambda_L S_leak - lambda_W
      W_sess - DeltaE_UV - M_sigma (notes line 292) is an internally-defined construction; (b) the benchmark NUMBER V_eff(r_soft)=1.74701126
      (notes line 372) is a recorded figure imported from the Session-I run (Appendix B.2, notes line 367), used as the target
      that fixes lambda_L. The candidate''s value (1.74701126) is the recorded benchmark figure, so the value-level classification
      is published_target (an imported recorded target, not derived). The formula itself would be internal_consistency. Recorded
      as ambiguous because the genuine fit content is the back-solve of lambda_L to this V_eff target (see fit_stage247_lambda_*).'
    citations:
    - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
      line: 372
      excerpt: V_{\rm eff}(r_{\rm soft})=1.74701126,
    - path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
      line: 213
      excerpt: Veff_obs = sp.Float("1.74701126")
  downstream_dependents:
  - '248'
  synthesis_ingested_at: '2026-06-11T19:37:44Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_27/fit_stage_247_v_eff__v_eff.yaml
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
