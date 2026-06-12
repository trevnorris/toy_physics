# Phase C Adversarial Audit Wrapper

Candidate: `fit_stage250_goldilocks_inputs`

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
id: fit_stage250_goldilocks_inputs
candidate_key: stage250_goldilocks_inputs
dry_run: false
dry_run_id: null
anchor_stages:
- '250'
file_line_citations:
- path: research/pde_ledger/scripts/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.py
  line: 176
  excerpt: chi_num = 21.73204372
  stage: '250'
parameter_names:
- chi_peak
- chi_raw
- g_UV
- m_s_ratio_1836
batch_id: null
status: provenance_built
status_timestamps:
  pending: '2026-06-11T08:21:51Z'
  scanned: '2026-06-11T08:21:51Z'
  provenance_built: '2026-06-11T19:37:57Z'
codex_session:
  by_parameter:
    chi_peak:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    chi_raw:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    g_UV:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    m_s_ratio_1836:
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
  - redteam_adversarial/provenance/fit_stage250_goldilocks_inputs__chi_peak.yaml
  - redteam_adversarial/provenance/fit_stage250_goldilocks_inputs__chi_raw.yaml
  - redteam_adversarial/provenance/fit_stage250_goldilocks_inputs__g_uv.yaml
  - redteam_adversarial/provenance/fit_stage250_goldilocks_inputs__m_s_ratio_1836.yaml
  phase_c_prompt: null
verdict:
  adversarial: null
  adjudication: null
modality_attribution:
- existing_provenance
modality_fragments:
- candidate_key: stage250_goldilocks_inputs
  modality: existing_provenance
  anchor_stage: '250'
  parameter_names:
  - g_UV
  - chi_peak
  - chi_raw
  - m_s_ratio_1836
  reason: The Goldilocks-window benchmark hangs on declared inputs g_UV=0.95 (a suspiciously round coupling), chi_peak=21.73204372
    and chi_raw=50.74399964 (no recorded origin), and m=1836.15267343 - the CODATA proton/electron mass ratio silently imported
    as a toy-model literal; none of the named provenance sources derive these, and the safe-window margin v_safe/v_crit=1.259
    is sensitive to all of them.
  citation:
    path: research/pde_ledger/scripts/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.py
    line: 176
    excerpt: chi_num = 21.73204372
    stage: '250'
phase_b_status: synthesis_complete
```

## Primary Sources

- research/pde_ledger/scripts/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.py

## Provenance Slice

```yaml
- schema_version: 1
  candidate_id: fit_stage250_goldilocks_inputs
  parameter_name: chi_peak
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage250_goldilocks_inputs__phase_b__chi_peak.md
  agent_synthesis_required: true
  dry_run_fixture: false
  fixture_note: null
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - '250'
  notes_source_opened: true
  source_evidence:
  - path: paper/stages/stage_250.tex
    line: 20
    role: paper_stage_tex
    excerpt: \lambda_{\rm eff}\sqrt{\frac{m_s}{2(E-V_{\rm peak})}}.
  - path: paper/stages/stage_250.tex
    line: 27
    role: paper_stage_tex
    excerpt: \mu_\eta\ddot V=g_{UV}\chi_{\rm peak}V.
  - path: paper/stages/stage_250.tex
    line: 33
    role: paper_stage_tex
    excerpt: \sqrt{\frac{\mu_\eta}{g_{UV}\chi_{\rm peak}}}.
  - path: paper/stages/stage_250.tex
    line: 35
    role: paper_stage_tex
    excerpt: The survival ratio is
  - path: paper/stages/stage_250.tex
    line: 37
    role: paper_stage_tex
    excerpt: \label{eq:app-part08-stage250-survival-ratio}
  - path: paper/stages/stage_250.tex
    line: 41
    role: paper_stage_tex
    excerpt: \sqrt{\frac{m_sg_{UV}\chi_{\rm peak}}{2\mu_\eta(E-V_{\rm peak})}}.
  - path: paper/stages/stage_250.tex
    line: 48
    role: paper_stage_tex
    excerpt: E_{\rm edge}=V_{\rm peak}+\frac{g_{UV}\chi_{\rm peak}\lambda_{\rm eff}^{2}}{2\alpha},
  - path: paper/stages/stage_250.tex
    line: 53
    role: paper_stage_tex
    excerpt: (\lambda_{\rm eff},\chi_{\rm peak},g_{UV},\mu_\eta/m_s).
  - path: paper/stages/stage_250.tex
    line: 64
    role: paper_stage_tex
    excerpt: E_{\rm edge}-V_{\rm peak}\propto\lambda_{\rm eff}^{2}\chi_{\rm peak}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 31
    role: notes_stage
    excerpt: (r_{\rm peak},V_{\rm peak}),
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 54
    role: notes_stage
    excerpt: 4. the stability ratio
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 69
    role: notes_stage
    excerpt: '- Stage 248 supplied the dynamic crossing branch, the peak energy \(V_{\rm peak}\), and the turning-point trigger
      width \(\lambda_{\rm th}\).'
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 131
    role: notes_stage
    excerpt: and used the peak-over-barrier kinetic speed
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 134
    role: notes_stage
    excerpt: v_{\rm bar}(E)=\sqrt{\frac{2(E-V_{\rm peak})}{m_s}}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 144
    role: notes_stage
    excerpt: \lambda_{\rm eff}\sqrt{\frac{m_s}{2(E-V_{\rm peak})}}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 152
    role: notes_stage
    excerpt: For fixed \((m_s,\lambda_{\rm eff},V_{\rm peak})\),
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 156
    role: notes_stage
    excerpt: -\frac{\lambda_{\rm eff}\sqrt{m_s}}{2\sqrt{2}\,(E-V_{\rm peak})^{3/2}}<0.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 167
    role: notes_stage
    excerpt: \mu_\eta\,\ddot V = g_{UV}\,\chi_{\rm peak}\,V,
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 175
    role: notes_stage
    excerpt: \chi_{\rm peak}:=\max\chi_\lambda(r)
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 181
    role: notes_stage
    excerpt: \Gamma_{\rm coll}=\sqrt{\frac{g_{UV}\chi_{\rm peak}}{\mu_\eta}},
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 190
    role: notes_stage
    excerpt: \sqrt{\frac{\mu_\eta}{g_{UV}\chi_{\rm peak}}}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 198
    role: notes_stage
    excerpt: '## 3. Stability ratio and lower safe edge'
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 200
    role: notes_stage
    excerpt: Define the characteristic survival ratio
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 212
    role: notes_stage
    excerpt: \sqrt{\frac{m_s g_{UV}\chi_{\rm peak}}{2\mu_\eta\,(E-V_{\rm peak})}}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 231
    role: notes_stage
    excerpt: V_{\rm peak}
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 234
    role: notes_stage
    excerpt: \frac{g_{UV}\chi_{\rm peak}\lambda_{\rm eff}^2}{2}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 282
    role: notes_stage
    excerpt: \frac{2(V_{\rm peak}-V(r_0))}{m_s},
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 291
    role: notes_stage
    excerpt: \frac{\lambda_{\rm eff}^2 g_{UV}\chi_{\rm peak}}{\mu_\eta}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 300
    role: notes_stage
    excerpt: E-V_{\rm peak}=\frac{m_s}{2}\bigl(v_0^2-v_{\rm crit,new}^2\bigr),
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 310
    role: notes_stage
    excerpt: The corresponding survival ratio is
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 315
    role: notes_stage
    excerpt: \frac{\lambda_{\rm eff}\sqrt{g_{UV}\chi_{\rm peak}/\mu_\eta}}{\sqrt{v_0^2-v_{\rm crit,new}^2}}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 334
    role: notes_stage
    excerpt: Substituting into the stability ratio gives
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 339
    role: notes_stage
    excerpt: \frac{\sqrt{2}\,\lambda_{\rm eff}\sqrt{g_{UV}\chi_{\rm peak}}}{2\sqrt{\alpha}}
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 340
    role: notes_stage
    excerpt: \,(E-V_{\rm peak})^{-1/2}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 348
    role: notes_stage
    excerpt: V_{\rm peak}
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 350
    role: notes_stage
    excerpt: \frac{g_{UV}\chi_{\rm peak}\lambda_{\rm eff}^2}{2\alpha}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 365
    role: notes_stage
    excerpt: \lambda_{\rm eff},\qquad \chi_{\rm peak},\qquad g_{UV},\qquad \alpha.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 386
    role: notes_stage
    excerpt: \chi_\lambda(r)=\lambda\,\bigl|\partial_r\ln V(r)\bigr|,
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 390
    role: notes_stage
    excerpt: \chi_{\rm peak}=\max \chi_\lambda(r),
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 397
    role: notes_stage
    excerpt: \bigl(\lambda_{\rm eff},\chi_{\rm peak}\bigr).
  origin_claims:
  - parameter: chi_peak
    introduced_at_stage: 250
    introduced_at_line: 417
    citation:
      path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
      line: 417
      excerpt: ' \chi_{\rm peak}=21.73204372,'
  constraints:
  - parameter: chi_peak
    constraint_kind: free_choice
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
      line: 388
      excerpt: \chi_\lambda(r)=\lambda\,\bigl|\partial_r\ln V(r)\bigr|,
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings:
  - type: provenance_gap
    severity: needs_triage
    summary: chi_peak = 21.73204372 (max log-barrier gradient with trigger width) enters at stage 250 as a Session-II/III
      benchmark readback off the Stage-247/248 barrier profile. The notes give the defining formula but the numeric value
      is a hardcoded benchmark input, not evaluated in-ledger; the number traces to an external session run.
    citations:
    - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
      line: 406
      excerpt: Use the Session-II/III reduced benchmark values
    - path: scripts/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.py
      line: 176
      excerpt: chi_num = 21.73204372
  downstream_dependents:
  - '251'
  synthesis_ingested_at: '2026-06-11T19:37:57Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_28/fit_stage250_goldilocks_inputs__chi_peak.yaml
- schema_version: 1
  candidate_id: fit_stage250_goldilocks_inputs
  parameter_name: chi_raw
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage250_goldilocks_inputs__phase_b__chi_raw.md
  agent_synthesis_required: true
  dry_run_fixture: false
  fixture_note: null
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - '250'
  notes_source_opened: true
  source_evidence:
  - path: paper/stages/stage_250.tex
    line: 20
    role: paper_stage_tex
    excerpt: \lambda_{\rm eff}\sqrt{\frac{m_s}{2(E-V_{\rm peak})}}.
  - path: paper/stages/stage_250.tex
    line: 27
    role: paper_stage_tex
    excerpt: \mu_\eta\ddot V=g_{UV}\chi_{\rm peak}V.
  - path: paper/stages/stage_250.tex
    line: 33
    role: paper_stage_tex
    excerpt: \sqrt{\frac{\mu_\eta}{g_{UV}\chi_{\rm peak}}}.
  - path: paper/stages/stage_250.tex
    line: 35
    role: paper_stage_tex
    excerpt: The survival ratio is
  - path: paper/stages/stage_250.tex
    line: 37
    role: paper_stage_tex
    excerpt: \label{eq:app-part08-stage250-survival-ratio}
  - path: paper/stages/stage_250.tex
    line: 41
    role: paper_stage_tex
    excerpt: \sqrt{\frac{m_sg_{UV}\chi_{\rm peak}}{2\mu_\eta(E-V_{\rm peak})}}.
  - path: paper/stages/stage_250.tex
    line: 48
    role: paper_stage_tex
    excerpt: E_{\rm edge}=V_{\rm peak}+\frac{g_{UV}\chi_{\rm peak}\lambda_{\rm eff}^{2}}{2\alpha},
  - path: paper/stages/stage_250.tex
    line: 53
    role: paper_stage_tex
    excerpt: (\lambda_{\rm eff},\chi_{\rm peak},g_{UV},\mu_\eta/m_s).
  - path: paper/stages/stage_250.tex
    line: 64
    role: paper_stage_tex
    excerpt: E_{\rm edge}-V_{\rm peak}\propto\lambda_{\rm eff}^{2}\chi_{\rm peak}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 31
    role: notes_stage
    excerpt: (r_{\rm peak},V_{\rm peak}),
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 54
    role: notes_stage
    excerpt: 4. the stability ratio
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 69
    role: notes_stage
    excerpt: '- Stage 248 supplied the dynamic crossing branch, the peak energy \(V_{\rm peak}\), and the turning-point trigger
      width \(\lambda_{\rm th}\).'
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 131
    role: notes_stage
    excerpt: and used the peak-over-barrier kinetic speed
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 134
    role: notes_stage
    excerpt: v_{\rm bar}(E)=\sqrt{\frac{2(E-V_{\rm peak})}{m_s}}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 144
    role: notes_stage
    excerpt: \lambda_{\rm eff}\sqrt{\frac{m_s}{2(E-V_{\rm peak})}}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 152
    role: notes_stage
    excerpt: For fixed \((m_s,\lambda_{\rm eff},V_{\rm peak})\),
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 156
    role: notes_stage
    excerpt: -\frac{\lambda_{\rm eff}\sqrt{m_s}}{2\sqrt{2}\,(E-V_{\rm peak})^{3/2}}<0.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 167
    role: notes_stage
    excerpt: \mu_\eta\,\ddot V = g_{UV}\,\chi_{\rm peak}\,V,
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 175
    role: notes_stage
    excerpt: \chi_{\rm peak}:=\max\chi_\lambda(r)
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 181
    role: notes_stage
    excerpt: \Gamma_{\rm coll}=\sqrt{\frac{g_{UV}\chi_{\rm peak}}{\mu_\eta}},
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 190
    role: notes_stage
    excerpt: \sqrt{\frac{\mu_\eta}{g_{UV}\chi_{\rm peak}}}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 198
    role: notes_stage
    excerpt: '## 3. Stability ratio and lower safe edge'
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 200
    role: notes_stage
    excerpt: Define the characteristic survival ratio
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 212
    role: notes_stage
    excerpt: \sqrt{\frac{m_s g_{UV}\chi_{\rm peak}}{2\mu_\eta\,(E-V_{\rm peak})}}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 231
    role: notes_stage
    excerpt: V_{\rm peak}
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 234
    role: notes_stage
    excerpt: \frac{g_{UV}\chi_{\rm peak}\lambda_{\rm eff}^2}{2}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 282
    role: notes_stage
    excerpt: \frac{2(V_{\rm peak}-V(r_0))}{m_s},
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 291
    role: notes_stage
    excerpt: \frac{\lambda_{\rm eff}^2 g_{UV}\chi_{\rm peak}}{\mu_\eta}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 300
    role: notes_stage
    excerpt: E-V_{\rm peak}=\frac{m_s}{2}\bigl(v_0^2-v_{\rm crit,new}^2\bigr),
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 310
    role: notes_stage
    excerpt: The corresponding survival ratio is
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 315
    role: notes_stage
    excerpt: \frac{\lambda_{\rm eff}\sqrt{g_{UV}\chi_{\rm peak}/\mu_\eta}}{\sqrt{v_0^2-v_{\rm crit,new}^2}}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 334
    role: notes_stage
    excerpt: Substituting into the stability ratio gives
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 339
    role: notes_stage
    excerpt: \frac{\sqrt{2}\,\lambda_{\rm eff}\sqrt{g_{UV}\chi_{\rm peak}}}{2\sqrt{\alpha}}
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 340
    role: notes_stage
    excerpt: \,(E-V_{\rm peak})^{-1/2}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 348
    role: notes_stage
    excerpt: V_{\rm peak}
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 350
    role: notes_stage
    excerpt: \frac{g_{UV}\chi_{\rm peak}\lambda_{\rm eff}^2}{2\alpha}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 365
    role: notes_stage
    excerpt: \lambda_{\rm eff},\qquad \chi_{\rm peak},\qquad g_{UV},\qquad \alpha.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 386
    role: notes_stage
    excerpt: \chi_\lambda(r)=\lambda\,\bigl|\partial_r\ln V(r)\bigr|,
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 390
    role: notes_stage
    excerpt: \chi_{\rm peak}=\max \chi_\lambda(r),
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 397
    role: notes_stage
    excerpt: \bigl(\lambda_{\rm eff},\chi_{\rm peak}\bigr).
  origin_claims:
  - parameter: chi_raw
    introduced_at_stage: 250
    introduced_at_line: 526
    citation:
      path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
      line: 526
      excerpt: \chi_{\rm peak}^{\rm raw}\approx 50.74399964.
  constraints:
  - parameter: chi_raw
    constraint_kind: free_choice
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
      line: 518
      excerpt: Session III also repeated the same calculation using the raw model width
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings:
  - type: provenance_gap
    severity: low
    summary: chi_raw = chi_peak^raw = 50.74399964 is the steepest log-gradient evaluated with raw width lambda=1 (sensitivity
      variant). A Session-III readback off the barrier profile; the variant condition is documented but the numeric value
      is a hardcoded benchmark input, not derived in-ledger.
    citations:
    - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
      line: 521
      excerpt: \lambda=1
    - path: scripts/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.py
      line: 218
      excerpt: chi_raw = 50.74399964
  downstream_dependents: []
  synthesis_ingested_at: '2026-06-11T19:37:57Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_28/fit_stage250_goldilocks_inputs__chi_raw.yaml
- schema_version: 1
  candidate_id: fit_stage250_goldilocks_inputs
  parameter_name: g_UV
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage250_goldilocks_inputs__phase_b__g_uv.md
  agent_synthesis_required: true
  dry_run_fixture: false
  fixture_note: null
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - '250'
  notes_source_opened: true
  source_evidence:
  - path: paper/stages/stage_250.tex
    line: 20
    role: paper_stage_tex
    excerpt: \lambda_{\rm eff}\sqrt{\frac{m_s}{2(E-V_{\rm peak})}}.
  - path: paper/stages/stage_250.tex
    line: 27
    role: paper_stage_tex
    excerpt: \mu_\eta\ddot V=g_{UV}\chi_{\rm peak}V.
  - path: paper/stages/stage_250.tex
    line: 33
    role: paper_stage_tex
    excerpt: \sqrt{\frac{\mu_\eta}{g_{UV}\chi_{\rm peak}}}.
  - path: paper/stages/stage_250.tex
    line: 35
    role: paper_stage_tex
    excerpt: The survival ratio is
  - path: paper/stages/stage_250.tex
    line: 37
    role: paper_stage_tex
    excerpt: \label{eq:app-part08-stage250-survival-ratio}
  - path: paper/stages/stage_250.tex
    line: 41
    role: paper_stage_tex
    excerpt: \sqrt{\frac{m_sg_{UV}\chi_{\rm peak}}{2\mu_\eta(E-V_{\rm peak})}}.
  - path: paper/stages/stage_250.tex
    line: 48
    role: paper_stage_tex
    excerpt: E_{\rm edge}=V_{\rm peak}+\frac{g_{UV}\chi_{\rm peak}\lambda_{\rm eff}^{2}}{2\alpha},
  - path: paper/stages/stage_250.tex
    line: 53
    role: paper_stage_tex
    excerpt: (\lambda_{\rm eff},\chi_{\rm peak},g_{UV},\mu_\eta/m_s).
  - path: paper/stages/stage_250.tex
    line: 64
    role: paper_stage_tex
    excerpt: E_{\rm edge}-V_{\rm peak}\propto\lambda_{\rm eff}^{2}\chi_{\rm peak}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 31
    role: notes_stage
    excerpt: (r_{\rm peak},V_{\rm peak}),
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 54
    role: notes_stage
    excerpt: 4. the stability ratio
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 69
    role: notes_stage
    excerpt: '- Stage 248 supplied the dynamic crossing branch, the peak energy \(V_{\rm peak}\), and the turning-point trigger
      width \(\lambda_{\rm th}\).'
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 131
    role: notes_stage
    excerpt: and used the peak-over-barrier kinetic speed
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 134
    role: notes_stage
    excerpt: v_{\rm bar}(E)=\sqrt{\frac{2(E-V_{\rm peak})}{m_s}}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 144
    role: notes_stage
    excerpt: \lambda_{\rm eff}\sqrt{\frac{m_s}{2(E-V_{\rm peak})}}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 152
    role: notes_stage
    excerpt: For fixed \((m_s,\lambda_{\rm eff},V_{\rm peak})\),
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 156
    role: notes_stage
    excerpt: -\frac{\lambda_{\rm eff}\sqrt{m_s}}{2\sqrt{2}\,(E-V_{\rm peak})^{3/2}}<0.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 167
    role: notes_stage
    excerpt: \mu_\eta\,\ddot V = g_{UV}\,\chi_{\rm peak}\,V,
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 175
    role: notes_stage
    excerpt: \chi_{\rm peak}:=\max\chi_\lambda(r)
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 181
    role: notes_stage
    excerpt: \Gamma_{\rm coll}=\sqrt{\frac{g_{UV}\chi_{\rm peak}}{\mu_\eta}},
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 190
    role: notes_stage
    excerpt: \sqrt{\frac{\mu_\eta}{g_{UV}\chi_{\rm peak}}}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 198
    role: notes_stage
    excerpt: '## 3. Stability ratio and lower safe edge'
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 200
    role: notes_stage
    excerpt: Define the characteristic survival ratio
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 212
    role: notes_stage
    excerpt: \sqrt{\frac{m_s g_{UV}\chi_{\rm peak}}{2\mu_\eta\,(E-V_{\rm peak})}}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 231
    role: notes_stage
    excerpt: V_{\rm peak}
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 234
    role: notes_stage
    excerpt: \frac{g_{UV}\chi_{\rm peak}\lambda_{\rm eff}^2}{2}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 282
    role: notes_stage
    excerpt: \frac{2(V_{\rm peak}-V(r_0))}{m_s},
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 291
    role: notes_stage
    excerpt: \frac{\lambda_{\rm eff}^2 g_{UV}\chi_{\rm peak}}{\mu_\eta}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 300
    role: notes_stage
    excerpt: E-V_{\rm peak}=\frac{m_s}{2}\bigl(v_0^2-v_{\rm crit,new}^2\bigr),
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 310
    role: notes_stage
    excerpt: The corresponding survival ratio is
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 315
    role: notes_stage
    excerpt: \frac{\lambda_{\rm eff}\sqrt{g_{UV}\chi_{\rm peak}/\mu_\eta}}{\sqrt{v_0^2-v_{\rm crit,new}^2}}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 334
    role: notes_stage
    excerpt: Substituting into the stability ratio gives
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 339
    role: notes_stage
    excerpt: \frac{\sqrt{2}\,\lambda_{\rm eff}\sqrt{g_{UV}\chi_{\rm peak}}}{2\sqrt{\alpha}}
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 340
    role: notes_stage
    excerpt: \,(E-V_{\rm peak})^{-1/2}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 348
    role: notes_stage
    excerpt: V_{\rm peak}
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 350
    role: notes_stage
    excerpt: \frac{g_{UV}\chi_{\rm peak}\lambda_{\rm eff}^2}{2\alpha}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 365
    role: notes_stage
    excerpt: \lambda_{\rm eff},\qquad \chi_{\rm peak},\qquad g_{UV},\qquad \alpha.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 386
    role: notes_stage
    excerpt: \chi_\lambda(r)=\lambda\,\bigl|\partial_r\ln V(r)\bigr|,
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 390
    role: notes_stage
    excerpt: \chi_{\rm peak}=\max \chi_\lambda(r),
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 397
    role: notes_stage
    excerpt: \bigl(\lambda_{\rm eff},\chi_{\rm peak}\bigr).
  origin_claims:
  - parameter: g_UV
    introduced_at_stage: 250
    introduced_at_line: 415
    citation:
      path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
      line: 415
      excerpt: ' g_{UV}=0.95,'
  constraints:
  - parameter: g_UV
    constraint_kind: free_choice
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
      line: 172
      excerpt: '- \(g_{UV}\) is the non-rigid transfer-to-dressing coupling carried from Stage 245,'
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings:
  - type: stale_provenance_anchor
    severity: low
    summary: g_UV=0.95 is described as 'carried from Stage 245', but the numeric value 0.95 appears nowhere in the Stage 245
      notes (only the symbol chi_UV := g_UV*chi_lambda). The number first enters the ledger here at stage 250 as a posited
      near-unity benchmark coupling (free_choice). The 'carried from Stage 245' attribution applies to the symbol, not the
      value.
    citations:
    - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
      line: 172
      excerpt: '- \(g_{UV}\) is the non-rigid transfer-to-dressing coupling carried from Stage 245,'
    - path: notes/stages/moving_throat_pde_stage245_nonrigid_mouth_dressing_packet_and_uv_drain_compiler_sympy_audit.md
      line: 134
      excerpt: \chi_{UV}(r):=g_{UV}\,\chi_\lambda(r),
  downstream_dependents:
  - '251'
  - '253'
  synthesis_ingested_at: '2026-06-11T19:37:57Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_28/fit_stage250_goldilocks_inputs__g_uv.yaml
- schema_version: 1
  candidate_id: fit_stage250_goldilocks_inputs
  parameter_name: m_s_ratio_1836
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage250_goldilocks_inputs__phase_b__m_s_ratio_1836.md
  agent_synthesis_required: true
  dry_run_fixture: false
  fixture_note: null
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - '250'
  notes_source_opened: true
  source_evidence:
  - path: paper/stages/stage_250.tex
    line: 20
    role: paper_stage_tex
    excerpt: \lambda_{\rm eff}\sqrt{\frac{m_s}{2(E-V_{\rm peak})}}.
  - path: paper/stages/stage_250.tex
    line: 27
    role: paper_stage_tex
    excerpt: \mu_\eta\ddot V=g_{UV}\chi_{\rm peak}V.
  - path: paper/stages/stage_250.tex
    line: 33
    role: paper_stage_tex
    excerpt: \sqrt{\frac{\mu_\eta}{g_{UV}\chi_{\rm peak}}}.
  - path: paper/stages/stage_250.tex
    line: 35
    role: paper_stage_tex
    excerpt: The survival ratio is
  - path: paper/stages/stage_250.tex
    line: 37
    role: paper_stage_tex
    excerpt: \label{eq:app-part08-stage250-survival-ratio}
  - path: paper/stages/stage_250.tex
    line: 41
    role: paper_stage_tex
    excerpt: \sqrt{\frac{m_sg_{UV}\chi_{\rm peak}}{2\mu_\eta(E-V_{\rm peak})}}.
  - path: paper/stages/stage_250.tex
    line: 48
    role: paper_stage_tex
    excerpt: E_{\rm edge}=V_{\rm peak}+\frac{g_{UV}\chi_{\rm peak}\lambda_{\rm eff}^{2}}{2\alpha},
  - path: paper/stages/stage_250.tex
    line: 53
    role: paper_stage_tex
    excerpt: (\lambda_{\rm eff},\chi_{\rm peak},g_{UV},\mu_\eta/m_s).
  - path: paper/stages/stage_250.tex
    line: 64
    role: paper_stage_tex
    excerpt: E_{\rm edge}-V_{\rm peak}\propto\lambda_{\rm eff}^{2}\chi_{\rm peak}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 31
    role: notes_stage
    excerpt: (r_{\rm peak},V_{\rm peak}),
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 54
    role: notes_stage
    excerpt: 4. the stability ratio
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 69
    role: notes_stage
    excerpt: '- Stage 248 supplied the dynamic crossing branch, the peak energy \(V_{\rm peak}\), and the turning-point trigger
      width \(\lambda_{\rm th}\).'
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 131
    role: notes_stage
    excerpt: and used the peak-over-barrier kinetic speed
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 134
    role: notes_stage
    excerpt: v_{\rm bar}(E)=\sqrt{\frac{2(E-V_{\rm peak})}{m_s}}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 144
    role: notes_stage
    excerpt: \lambda_{\rm eff}\sqrt{\frac{m_s}{2(E-V_{\rm peak})}}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 152
    role: notes_stage
    excerpt: For fixed \((m_s,\lambda_{\rm eff},V_{\rm peak})\),
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 156
    role: notes_stage
    excerpt: -\frac{\lambda_{\rm eff}\sqrt{m_s}}{2\sqrt{2}\,(E-V_{\rm peak})^{3/2}}<0.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 167
    role: notes_stage
    excerpt: \mu_\eta\,\ddot V = g_{UV}\,\chi_{\rm peak}\,V,
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 175
    role: notes_stage
    excerpt: \chi_{\rm peak}:=\max\chi_\lambda(r)
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 181
    role: notes_stage
    excerpt: \Gamma_{\rm coll}=\sqrt{\frac{g_{UV}\chi_{\rm peak}}{\mu_\eta}},
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 190
    role: notes_stage
    excerpt: \sqrt{\frac{\mu_\eta}{g_{UV}\chi_{\rm peak}}}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 198
    role: notes_stage
    excerpt: '## 3. Stability ratio and lower safe edge'
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 200
    role: notes_stage
    excerpt: Define the characteristic survival ratio
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 212
    role: notes_stage
    excerpt: \sqrt{\frac{m_s g_{UV}\chi_{\rm peak}}{2\mu_\eta\,(E-V_{\rm peak})}}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 231
    role: notes_stage
    excerpt: V_{\rm peak}
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 234
    role: notes_stage
    excerpt: \frac{g_{UV}\chi_{\rm peak}\lambda_{\rm eff}^2}{2}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 282
    role: notes_stage
    excerpt: \frac{2(V_{\rm peak}-V(r_0))}{m_s},
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 291
    role: notes_stage
    excerpt: \frac{\lambda_{\rm eff}^2 g_{UV}\chi_{\rm peak}}{\mu_\eta}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 300
    role: notes_stage
    excerpt: E-V_{\rm peak}=\frac{m_s}{2}\bigl(v_0^2-v_{\rm crit,new}^2\bigr),
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 310
    role: notes_stage
    excerpt: The corresponding survival ratio is
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 315
    role: notes_stage
    excerpt: \frac{\lambda_{\rm eff}\sqrt{g_{UV}\chi_{\rm peak}/\mu_\eta}}{\sqrt{v_0^2-v_{\rm crit,new}^2}}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 334
    role: notes_stage
    excerpt: Substituting into the stability ratio gives
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 339
    role: notes_stage
    excerpt: \frac{\sqrt{2}\,\lambda_{\rm eff}\sqrt{g_{UV}\chi_{\rm peak}}}{2\sqrt{\alpha}}
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 340
    role: notes_stage
    excerpt: \,(E-V_{\rm peak})^{-1/2}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 348
    role: notes_stage
    excerpt: V_{\rm peak}
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 350
    role: notes_stage
    excerpt: \frac{g_{UV}\chi_{\rm peak}\lambda_{\rm eff}^2}{2\alpha}.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 365
    role: notes_stage
    excerpt: \lambda_{\rm eff},\qquad \chi_{\rm peak},\qquad g_{UV},\qquad \alpha.
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 386
    role: notes_stage
    excerpt: \chi_\lambda(r)=\lambda\,\bigl|\partial_r\ln V(r)\bigr|,
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 390
    role: notes_stage
    excerpt: \chi_{\rm peak}=\max \chi_\lambda(r),
  - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
    line: 397
    role: notes_stage
    excerpt: \bigl(\lambda_{\rm eff},\chi_{\rm peak}\bigr).
  origin_claims:
  - parameter: m_s_ratio_1836
    introduced_at_stage: 250
    introduced_at_line: 419
    citation:
      path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
      line: 419
      excerpt: ' m_s=\mu_\eta=1836.15267343.'
  constraints:
  - parameter: m_s_ratio_1836
    constraint_kind: published_target
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
      line: 423
      excerpt: '### 7.1 Proton-proxy classical threshold speed'
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings:
  - type: published_target
    severity: high
    summary: '1836.15267343 is the CODATA proton-to-electron mass ratio, imported into the ledger as the ''proton-proxy''
      benchmark mass (set for both m_s and mu_eta). This is the load-bearing external-target import flagged for this batch:
      a known published constant fitted in as a benchmark, with the script labelling it an ''Imported benchmark value''. The
      notes are honest that it is a ''proton-proxy'', so there is no derive-vs-fit overclaim -- but it IS an external published
      number, not a ledger-derived quantity.'
    citations:
    - path: notes/stages/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.md
      line: 419
      excerpt: ' m_s=\mu_\eta=1836.15267343.'
    - path: scripts/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.py
      line: 171
      excerpt: '# Imported benchmark values from the barrier-session write-up.'
  downstream_dependents:
  - '251'
  - '253'
  synthesis_ingested_at: '2026-06-11T19:37:57Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_28/fit_stage250_goldilocks_inputs__m_s_ratio_1836.yaml
```

## Benchmarks

External-match checks must use these sourced benchmark entries. Do not adjudicate a match from model memory.

```yaml
- family_id: fam_0458_m_s_ratio_1836
  claim: Ledger imports 1836.15267343 (stage 250 proton-proxy benchmark) as a match to the CODATA proton-to-electron mass
    ratio m_p/m_e.
  value: '1836.152 673 43(11)  [2018 CODATA: standard uncertainty 0.000 000 11, relative u_r = 6.0e-11]. Current 2022 CODATA
    value: 1836.152 673 426(32), relative u_r = 1.7e-11.'
  source_type: CODATA
  source_citation: NIST CODATA, proton-electron mass ratio (m_p/m_e), https://physics.nist.gov/cgi-bin/cuu/Value?mpsme ; CODATA
    Internationally recommended values of the fundamental physical constants (2018 and 2022 adjustments).
  obtained_note: 'WebFetch of https://physics.nist.gov/cgi-bin/cuu/Value?mpsme on 2026-06-11 returned the current 2022 value
    1836.152 673 426(32), u_r=1.7e-11. WebSearch (NIST/CODATA 2018 sources incl. physics.nist.gov/cuu/pdf/wallet_2018.pdf,
    all_2018.pdf) confirmed the 2018 value 1836.152 673 43(11), u=0.000 000 11, u_r=6.0e-11. AGREEMENT: the ledger''s 1836.15267343
    matches the 2018 CODATA value EXACTLY to all printed digits; it is the 2018 release (NOT 2014, NOT 2022). The 2022 value
    shifted to 1836.152 673 426(32), still consistent within uncertainty but with more digits, so the ledger number is specifically
    the 2018 CODATA recommended value.'
  id: bench_5d8d048b10f9
  ingested_at: '2026-06-11T19:49:46Z'
  ingested_from: redteam_adversarial/provenance/_benchmarks_sourced.yaml
```

## Graph Context

```yaml
contexts: []
graph_gaps: []
```

Output a markdown adversarial report only when this is not a dry run. For dry runs, stop after confirming the prompt has enough context and record no verdict.
