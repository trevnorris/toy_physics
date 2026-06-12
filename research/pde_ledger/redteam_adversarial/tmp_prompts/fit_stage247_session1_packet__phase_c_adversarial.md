# Phase C Adversarial Audit Wrapper

Candidate: `fit_stage247_session1_packet`

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
id: fit_stage247_session1_packet
candidate_key: stage247_session1_packet
dry_run: false
dry_run_id: null
anchor_stages:
- '247'
file_line_citations:
- path: research/pde_ledger/scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
  line: 213
  excerpt: Veff_obs = sp.Float("1.74701126")
  stage: '247'
parameter_names:
- G_W
- Rmix
- UVdrop_obs
- Veff_obs
- Wsess_obs
- beta_Q
- beta_U
- beta_W
- eta_leak
- mu_w
- xi_R
batch_id: null
status: provenance_built
status_timestamps:
  pending: '2026-06-11T08:21:51Z'
  scanned: '2026-06-11T08:21:51Z'
  provenance_built: '2026-06-11T19:37:44Z'
codex_session:
  by_parameter:
    G_W:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    Rmix:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    UVdrop_obs:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    Veff_obs:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    Wsess_obs:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    beta_Q:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    beta_U:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    beta_W:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    eta_leak:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    mu_w:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    xi_R:
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
  - redteam_adversarial/provenance/fit_stage247_session1_packet__g_w.yaml
  - redteam_adversarial/provenance/fit_stage247_session1_packet__rmix.yaml
  - redteam_adversarial/provenance/fit_stage247_session1_packet__uvdrop_obs.yaml
  - redteam_adversarial/provenance/fit_stage247_session1_packet__veff_obs.yaml
  - redteam_adversarial/provenance/fit_stage247_session1_packet__wsess_obs.yaml
  - redteam_adversarial/provenance/fit_stage247_session1_packet__beta_q.yaml
  - redteam_adversarial/provenance/fit_stage247_session1_packet__beta_u.yaml
  - redteam_adversarial/provenance/fit_stage247_session1_packet__beta_w.yaml
  - redteam_adversarial/provenance/fit_stage247_session1_packet__eta_leak.yaml
  - redteam_adversarial/provenance/fit_stage247_session1_packet__mu_w.yaml
  - redteam_adversarial/provenance/fit_stage247_session1_packet__xi_r.yaml
  phase_c_prompt: null
verdict:
  adversarial: null
  adjudication: null
modality_attribution:
- existing_provenance
modality_fragments:
- candidate_key: stage247_session1_packet
  modality: existing_provenance
  anchor_stage: '247'
  parameter_names:
  - Veff_obs
  - Wsess_obs
  - UVdrop_obs
  - Rmix
  - G_W
  - beta_Q
  - beta_U
  - beta_W
  - xi_R
  - eta_leak
  - mu_w
  reason: The Stage-247 benchmark consists of a 23-entry hand-set substitution dict (Kstar=4.0, OmU2=9.0, OmW2=16.0, GW=1.25,
    Rmix=1.35, betaQ=0.03, betaU=0.15, betaW=0.20, xi_R=0.9, eta_leak=0.03, mu_w=0.8, ...) plus three "observed" scalars (Veff_obs/Wsess_obs/UVdrop_obs)
    with no in-ledger derivation in any named provenance source - the entire one-point decomposition could be a fitted slice.
  citation:
    path: research/pde_ledger/scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
    line: 213
    excerpt: Veff_obs = sp.Float("1.74701126")
    stage: '247'
phase_b_status: synthesis_complete
```

## Primary Sources

- research/pde_ledger/scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py

## Provenance Slice

```yaml
- schema_version: 1
  candidate_id: fit_stage247_session1_packet
  parameter_name: G_W
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage247_session1_packet__phase_b__g_w.md
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
    line: 4
    role: paper_stage_tex
    excerpt: '\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py}.  Mathematica
      audit: none yet.}'
  - path: paper/stages/stage_247.tex
    line: 12
    role: paper_stage_tex
    excerpt: Q=G_U^2\Omega_W^2+2G_UG_WR_{\rm mix}+G_W^2\Omega_U^2,
  - path: paper/stages/stage_247.tex
    line: 16
    role: paper_stage_tex
    excerpt: P=\Omega_U^2G_W+R_{\rm mix}G_U,
  - path: paper/stages/stage_247.tex
    line: 20
    role: paper_stage_tex
    excerpt: With primitive short-range source amplitudes \(\beta_Q,\beta_U,\beta_W\), the product-family coefficients are
  - path: paper/stages/stage_247.tex
    line: 23
    role: paper_stage_tex
    excerpt: \mathcal C_6=\chi_{qq}\beta_Q^2,
  - path: paper/stages/stage_247.tex
    line: 25
    role: paper_stage_tex
    excerpt: \mathcal C_4=\chi_{qU}\beta_Q\beta_U+\chi_{qW}\beta_Q\beta_W,
  - path: paper/stages/stage_247.tex
    line: 29
    role: paper_stage_tex
    excerpt: \mathcal C_2=\chi_{UU}\beta_U^2+2\chi_{UW}\beta_U\beta_W+\chi_{WW}\beta_W^2.
  - path: paper/stages/stage_247.tex
    line: 60
    role: paper_stage_tex
    excerpt: S_{\rm leak}=\frac{8\sqrt2\eta_{\rm leak}\mu_w\rho_0}{\pi^{5/2}\lambda^3}\mathfrak L,
  - path: paper/stages/stage_247.tex
    line: 66
    role: paper_stage_tex
    excerpt: \frac{512\eta_{\rm leak}^{2}\mu_wq\rho_0}{\pi^4\lambda^2}\mathfrak L^2,
  - path: paper/stages/stage_247.tex
    line: 68
    role: paper_stage_tex
    excerpt: \Delta E_{UV}=\eta_{UV}\mathcal D_{UV},
  - path: paper/stages/stage_247.tex
    line: 72
    role: paper_stage_tex
    excerpt: \mathcal M_\sigma(r)=\xi_R\bigl[\mathcal R_\infty-\mathcal R(r)\bigr],
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
    line: 1
    role: notes_stage
    excerpt: '# Moving-Throat PDE — Stage 247: Relaxed Stationary Barrier Compiler from One-Port Short-Range, Leakage, `U/V`
      Drain, and Compensated Source Packets'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 17
    role: notes_stage
    excerpt: 2. the exact Stage-244 leakage / work compiler,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 33
    role: notes_stage
    excerpt: (S_{\rm leak},\,\mathcal W_w),
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 48
    role: notes_stage
    excerpt: 'It has four jobs:'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 51
    role: notes_stage
    excerpt: 2. import the Stage-244 leakage / work packet into barrier units,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 61
    role: notes_stage
    excerpt: '- **Stage 244**, which attached leakage and scalar-photon work to the selected-support packet,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 63
    role: notes_stage
    excerpt: '- **Stage 246**, which converted the sign-changing source branch into explicit mouth/source observables,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 98
    role: notes_stage
    excerpt: Q=G_U^2\Omega_W^2+2G_UG_WR_{\rm mix}+G_W^2\Omega_U^2,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 101
    role: notes_stage
    excerpt: P=\Omega_U^2G_W+R_{\rm mix}G_U,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 116
    role: notes_stage
    excerpt: \chi_{UU}=\frac{K_*\Omega_W^2-G_W^2}{\Delta D_0},
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 118
    role: notes_stage
    excerpt: \chi_{UW}=\frac{K_*R_{\rm mix}+G_UG_W}{\Delta D_0},
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 124
    role: notes_stage
    excerpt: P_U=G_U\Omega_W^2+R_{\rm mix}G_W.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 129
    role: notes_stage
    excerpt: \beta_Q,\qquad \beta_U,\qquad \beta_W,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 133
    role: notes_stage
    excerpt: \mathcal C_6=\chi_{qq}\beta_Q^2,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 136
    role: notes_stage
    excerpt: \mathcal C_4=\chi_{qU}\beta_Q\beta_U+\chi_{qW}\beta_Q\beta_W,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 139
    role: notes_stage
    excerpt: \mathcal C_2=\chi_{UU}\beta_U^2+2\chi_{UW}\beta_U\beta_W+\chi_{WW}\beta_W^2.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 187
    role: notes_stage
    excerpt: '### 2.1 Stage-244 leakage and reduced Session-I work scalar'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 196
    role: notes_stage
    excerpt: S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 198
    role: notes_stage
    excerpt: \frac{8\sqrt2\,\eta_{\rm leak}\mu_w\rho_0}{\pi^{5/2}\lambda^3}
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 206
    role: notes_stage
    excerpt: \frac{512\,\eta_{\rm leak}^2\mu_w q\rho_0}{\pi^4\lambda^2}
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 225
    role: notes_stage
    excerpt: \Delta E_{UV}(r):=\eta_{UV}\mathcal D_{UV}(r),
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 228
    role: notes_stage
    excerpt: where `\eta_{UV}` is the same Session-I barrier-lowering weight carried in the stationary script.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 247
    role: notes_stage
    excerpt: \mathcal M_\sigma(r):=\xi_R\bigl[\mathcal R_\infty-\mathcal R(r)\bigr],
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 296
    role: notes_stage
    excerpt: \lambda_L S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 308
    role: notes_stage
    excerpt: '- `\lambda_L\ge 0` converts the open-subsystem leakage source into stationary barrier units,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 320
    role: notes_stage
    excerpt: \lambda_L S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 344
    role: notes_stage
    excerpt: '- `(S_{\rm leak},\mathcal W_w^{\rm sess})` are support-side / open-system objects,'
  origin_claims:
  - parameter: G_W
    introduced_at_stage: 247
    introduced_at_line: 391
    citation:
      path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
      line: 391
      excerpt: G_W=1.25,
  constraints:
  - parameter: G_W
    constraint_kind: free_choice
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
      line: 381
      excerpt: Using the recorded one-port parameters
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings: []
  downstream_dependents:
  - '248'
  synthesis_ingested_at: '2026-06-11T19:37:44Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_27/fit_stage247_session1_packet__g_w.yaml
- schema_version: 1
  candidate_id: fit_stage247_session1_packet
  parameter_name: Rmix
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage247_session1_packet__phase_b__rmix.md
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
    line: 4
    role: paper_stage_tex
    excerpt: '\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py}.  Mathematica
      audit: none yet.}'
  - path: paper/stages/stage_247.tex
    line: 12
    role: paper_stage_tex
    excerpt: Q=G_U^2\Omega_W^2+2G_UG_WR_{\rm mix}+G_W^2\Omega_U^2,
  - path: paper/stages/stage_247.tex
    line: 16
    role: paper_stage_tex
    excerpt: P=\Omega_U^2G_W+R_{\rm mix}G_U,
  - path: paper/stages/stage_247.tex
    line: 20
    role: paper_stage_tex
    excerpt: With primitive short-range source amplitudes \(\beta_Q,\beta_U,\beta_W\), the product-family coefficients are
  - path: paper/stages/stage_247.tex
    line: 23
    role: paper_stage_tex
    excerpt: \mathcal C_6=\chi_{qq}\beta_Q^2,
  - path: paper/stages/stage_247.tex
    line: 25
    role: paper_stage_tex
    excerpt: \mathcal C_4=\chi_{qU}\beta_Q\beta_U+\chi_{qW}\beta_Q\beta_W,
  - path: paper/stages/stage_247.tex
    line: 29
    role: paper_stage_tex
    excerpt: \mathcal C_2=\chi_{UU}\beta_U^2+2\chi_{UW}\beta_U\beta_W+\chi_{WW}\beta_W^2.
  - path: paper/stages/stage_247.tex
    line: 60
    role: paper_stage_tex
    excerpt: S_{\rm leak}=\frac{8\sqrt2\eta_{\rm leak}\mu_w\rho_0}{\pi^{5/2}\lambda^3}\mathfrak L,
  - path: paper/stages/stage_247.tex
    line: 66
    role: paper_stage_tex
    excerpt: \frac{512\eta_{\rm leak}^{2}\mu_wq\rho_0}{\pi^4\lambda^2}\mathfrak L^2,
  - path: paper/stages/stage_247.tex
    line: 68
    role: paper_stage_tex
    excerpt: \Delta E_{UV}=\eta_{UV}\mathcal D_{UV},
  - path: paper/stages/stage_247.tex
    line: 72
    role: paper_stage_tex
    excerpt: \mathcal M_\sigma(r)=\xi_R\bigl[\mathcal R_\infty-\mathcal R(r)\bigr],
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
    line: 1
    role: notes_stage
    excerpt: '# Moving-Throat PDE — Stage 247: Relaxed Stationary Barrier Compiler from One-Port Short-Range, Leakage, `U/V`
      Drain, and Compensated Source Packets'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 17
    role: notes_stage
    excerpt: 2. the exact Stage-244 leakage / work compiler,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 33
    role: notes_stage
    excerpt: (S_{\rm leak},\,\mathcal W_w),
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 48
    role: notes_stage
    excerpt: 'It has four jobs:'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 51
    role: notes_stage
    excerpt: 2. import the Stage-244 leakage / work packet into barrier units,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 61
    role: notes_stage
    excerpt: '- **Stage 244**, which attached leakage and scalar-photon work to the selected-support packet,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 63
    role: notes_stage
    excerpt: '- **Stage 246**, which converted the sign-changing source branch into explicit mouth/source observables,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 98
    role: notes_stage
    excerpt: Q=G_U^2\Omega_W^2+2G_UG_WR_{\rm mix}+G_W^2\Omega_U^2,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 101
    role: notes_stage
    excerpt: P=\Omega_U^2G_W+R_{\rm mix}G_U,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 116
    role: notes_stage
    excerpt: \chi_{UU}=\frac{K_*\Omega_W^2-G_W^2}{\Delta D_0},
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 118
    role: notes_stage
    excerpt: \chi_{UW}=\frac{K_*R_{\rm mix}+G_UG_W}{\Delta D_0},
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 124
    role: notes_stage
    excerpt: P_U=G_U\Omega_W^2+R_{\rm mix}G_W.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 129
    role: notes_stage
    excerpt: \beta_Q,\qquad \beta_U,\qquad \beta_W,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 133
    role: notes_stage
    excerpt: \mathcal C_6=\chi_{qq}\beta_Q^2,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 136
    role: notes_stage
    excerpt: \mathcal C_4=\chi_{qU}\beta_Q\beta_U+\chi_{qW}\beta_Q\beta_W,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 139
    role: notes_stage
    excerpt: \mathcal C_2=\chi_{UU}\beta_U^2+2\chi_{UW}\beta_U\beta_W+\chi_{WW}\beta_W^2.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 187
    role: notes_stage
    excerpt: '### 2.1 Stage-244 leakage and reduced Session-I work scalar'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 196
    role: notes_stage
    excerpt: S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 198
    role: notes_stage
    excerpt: \frac{8\sqrt2\,\eta_{\rm leak}\mu_w\rho_0}{\pi^{5/2}\lambda^3}
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 206
    role: notes_stage
    excerpt: \frac{512\,\eta_{\rm leak}^2\mu_w q\rho_0}{\pi^4\lambda^2}
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 225
    role: notes_stage
    excerpt: \Delta E_{UV}(r):=\eta_{UV}\mathcal D_{UV}(r),
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 228
    role: notes_stage
    excerpt: where `\eta_{UV}` is the same Session-I barrier-lowering weight carried in the stationary script.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 247
    role: notes_stage
    excerpt: \mathcal M_\sigma(r):=\xi_R\bigl[\mathcal R_\infty-\mathcal R(r)\bigr],
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 296
    role: notes_stage
    excerpt: \lambda_L S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 308
    role: notes_stage
    excerpt: '- `\lambda_L\ge 0` converts the open-subsystem leakage source into stationary barrier units,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 320
    role: notes_stage
    excerpt: \lambda_L S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 344
    role: notes_stage
    excerpt: '- `(S_{\rm leak},\mathcal W_w^{\rm sess})` are support-side / open-system objects,'
  origin_claims:
  - parameter: Rmix
    introduced_at_stage: 247
    introduced_at_line: 393
    citation:
      path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
      line: 393
      excerpt: R_{\rm mix}=1.35,
  constraints:
  - parameter: Rmix
    constraint_kind: free_choice
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
      line: 381
      excerpt: Using the recorded one-port parameters
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings: []
  downstream_dependents:
  - '248'
  synthesis_ingested_at: '2026-06-11T19:37:44Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_27/fit_stage247_session1_packet__rmix.yaml
- schema_version: 1
  candidate_id: fit_stage247_session1_packet
  parameter_name: UVdrop_obs
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage247_session1_packet__phase_b__uvdrop_obs.md
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
    line: 4
    role: paper_stage_tex
    excerpt: '\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py}.  Mathematica
      audit: none yet.}'
  - path: paper/stages/stage_247.tex
    line: 12
    role: paper_stage_tex
    excerpt: Q=G_U^2\Omega_W^2+2G_UG_WR_{\rm mix}+G_W^2\Omega_U^2,
  - path: paper/stages/stage_247.tex
    line: 16
    role: paper_stage_tex
    excerpt: P=\Omega_U^2G_W+R_{\rm mix}G_U,
  - path: paper/stages/stage_247.tex
    line: 20
    role: paper_stage_tex
    excerpt: With primitive short-range source amplitudes \(\beta_Q,\beta_U,\beta_W\), the product-family coefficients are
  - path: paper/stages/stage_247.tex
    line: 23
    role: paper_stage_tex
    excerpt: \mathcal C_6=\chi_{qq}\beta_Q^2,
  - path: paper/stages/stage_247.tex
    line: 25
    role: paper_stage_tex
    excerpt: \mathcal C_4=\chi_{qU}\beta_Q\beta_U+\chi_{qW}\beta_Q\beta_W,
  - path: paper/stages/stage_247.tex
    line: 29
    role: paper_stage_tex
    excerpt: \mathcal C_2=\chi_{UU}\beta_U^2+2\chi_{UW}\beta_U\beta_W+\chi_{WW}\beta_W^2.
  - path: paper/stages/stage_247.tex
    line: 60
    role: paper_stage_tex
    excerpt: S_{\rm leak}=\frac{8\sqrt2\eta_{\rm leak}\mu_w\rho_0}{\pi^{5/2}\lambda^3}\mathfrak L,
  - path: paper/stages/stage_247.tex
    line: 66
    role: paper_stage_tex
    excerpt: \frac{512\eta_{\rm leak}^{2}\mu_wq\rho_0}{\pi^4\lambda^2}\mathfrak L^2,
  - path: paper/stages/stage_247.tex
    line: 68
    role: paper_stage_tex
    excerpt: \Delta E_{UV}=\eta_{UV}\mathcal D_{UV},
  - path: paper/stages/stage_247.tex
    line: 72
    role: paper_stage_tex
    excerpt: \mathcal M_\sigma(r)=\xi_R\bigl[\mathcal R_\infty-\mathcal R(r)\bigr],
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
    line: 1
    role: notes_stage
    excerpt: '# Moving-Throat PDE — Stage 247: Relaxed Stationary Barrier Compiler from One-Port Short-Range, Leakage, `U/V`
      Drain, and Compensated Source Packets'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 17
    role: notes_stage
    excerpt: 2. the exact Stage-244 leakage / work compiler,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 33
    role: notes_stage
    excerpt: (S_{\rm leak},\,\mathcal W_w),
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 48
    role: notes_stage
    excerpt: 'It has four jobs:'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 51
    role: notes_stage
    excerpt: 2. import the Stage-244 leakage / work packet into barrier units,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 61
    role: notes_stage
    excerpt: '- **Stage 244**, which attached leakage and scalar-photon work to the selected-support packet,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 63
    role: notes_stage
    excerpt: '- **Stage 246**, which converted the sign-changing source branch into explicit mouth/source observables,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 98
    role: notes_stage
    excerpt: Q=G_U^2\Omega_W^2+2G_UG_WR_{\rm mix}+G_W^2\Omega_U^2,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 101
    role: notes_stage
    excerpt: P=\Omega_U^2G_W+R_{\rm mix}G_U,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 116
    role: notes_stage
    excerpt: \chi_{UU}=\frac{K_*\Omega_W^2-G_W^2}{\Delta D_0},
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 118
    role: notes_stage
    excerpt: \chi_{UW}=\frac{K_*R_{\rm mix}+G_UG_W}{\Delta D_0},
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 124
    role: notes_stage
    excerpt: P_U=G_U\Omega_W^2+R_{\rm mix}G_W.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 129
    role: notes_stage
    excerpt: \beta_Q,\qquad \beta_U,\qquad \beta_W,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 133
    role: notes_stage
    excerpt: \mathcal C_6=\chi_{qq}\beta_Q^2,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 136
    role: notes_stage
    excerpt: \mathcal C_4=\chi_{qU}\beta_Q\beta_U+\chi_{qW}\beta_Q\beta_W,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 139
    role: notes_stage
    excerpt: \mathcal C_2=\chi_{UU}\beta_U^2+2\chi_{UW}\beta_U\beta_W+\chi_{WW}\beta_W^2.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 187
    role: notes_stage
    excerpt: '### 2.1 Stage-244 leakage and reduced Session-I work scalar'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 196
    role: notes_stage
    excerpt: S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 198
    role: notes_stage
    excerpt: \frac{8\sqrt2\,\eta_{\rm leak}\mu_w\rho_0}{\pi^{5/2}\lambda^3}
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 206
    role: notes_stage
    excerpt: \frac{512\,\eta_{\rm leak}^2\mu_w q\rho_0}{\pi^4\lambda^2}
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 225
    role: notes_stage
    excerpt: \Delta E_{UV}(r):=\eta_{UV}\mathcal D_{UV}(r),
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 228
    role: notes_stage
    excerpt: where `\eta_{UV}` is the same Session-I barrier-lowering weight carried in the stationary script.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 247
    role: notes_stage
    excerpt: \mathcal M_\sigma(r):=\xi_R\bigl[\mathcal R_\infty-\mathcal R(r)\bigr],
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 296
    role: notes_stage
    excerpt: \lambda_L S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 308
    role: notes_stage
    excerpt: '- `\lambda_L\ge 0` converts the open-subsystem leakage source into stationary barrier units,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 320
    role: notes_stage
    excerpt: \lambda_L S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 344
    role: notes_stage
    excerpt: '- `(S_{\rm leak},\mathcal W_w^{\rm sess})` are support-side / open-system objects,'
  origin_claims:
  - parameter: UVdrop_obs
    introduced_at_stage: 247
    introduced_at_line: 431
    citation:
      path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
      line: 431
      excerpt: \Delta E_{UV}(r_{\rm soft})=0.21064278.
  constraints:
  - parameter: UVdrop_obs
    constraint_kind: published_target
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
      line: 427
      excerpt: The Session-I table also records
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings:
  - type: published_target
    severity: needs_triage
    summary: UVdrop_obs (DeltaE_UV(r_soft))=0.21064278 is a recorded Session-I benchmark figure (notes line 427-431; script
      line 215 UVdrop_obs = sp.Float('0.21064278')). It is read directly from the recorded stationary run rather than recomputed
      from the Stage-245 drain formula on this slice, and enters the lambda_L back-solve numerator (script line 238). A hardcoded
      recorded target, not a stage-247 derivation.
    citations:
    - path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
      line: 215
      excerpt: UVdrop_obs = sp.Float("0.21064278")
    - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
      line: 431
      excerpt: \Delta E_{UV}(r_{\rm soft})=0.21064278.
  downstream_dependents:
  - '248'
  synthesis_ingested_at: '2026-06-11T19:37:44Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_27/fit_stage247_session1_packet__uvdrop_obs.yaml
- schema_version: 1
  candidate_id: fit_stage247_session1_packet
  parameter_name: Veff_obs
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage247_session1_packet__phase_b__veff_obs.md
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
    line: 4
    role: paper_stage_tex
    excerpt: '\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py}.  Mathematica
      audit: none yet.}'
  - path: paper/stages/stage_247.tex
    line: 12
    role: paper_stage_tex
    excerpt: Q=G_U^2\Omega_W^2+2G_UG_WR_{\rm mix}+G_W^2\Omega_U^2,
  - path: paper/stages/stage_247.tex
    line: 16
    role: paper_stage_tex
    excerpt: P=\Omega_U^2G_W+R_{\rm mix}G_U,
  - path: paper/stages/stage_247.tex
    line: 20
    role: paper_stage_tex
    excerpt: With primitive short-range source amplitudes \(\beta_Q,\beta_U,\beta_W\), the product-family coefficients are
  - path: paper/stages/stage_247.tex
    line: 23
    role: paper_stage_tex
    excerpt: \mathcal C_6=\chi_{qq}\beta_Q^2,
  - path: paper/stages/stage_247.tex
    line: 25
    role: paper_stage_tex
    excerpt: \mathcal C_4=\chi_{qU}\beta_Q\beta_U+\chi_{qW}\beta_Q\beta_W,
  - path: paper/stages/stage_247.tex
    line: 29
    role: paper_stage_tex
    excerpt: \mathcal C_2=\chi_{UU}\beta_U^2+2\chi_{UW}\beta_U\beta_W+\chi_{WW}\beta_W^2.
  - path: paper/stages/stage_247.tex
    line: 60
    role: paper_stage_tex
    excerpt: S_{\rm leak}=\frac{8\sqrt2\eta_{\rm leak}\mu_w\rho_0}{\pi^{5/2}\lambda^3}\mathfrak L,
  - path: paper/stages/stage_247.tex
    line: 66
    role: paper_stage_tex
    excerpt: \frac{512\eta_{\rm leak}^{2}\mu_wq\rho_0}{\pi^4\lambda^2}\mathfrak L^2,
  - path: paper/stages/stage_247.tex
    line: 68
    role: paper_stage_tex
    excerpt: \Delta E_{UV}=\eta_{UV}\mathcal D_{UV},
  - path: paper/stages/stage_247.tex
    line: 72
    role: paper_stage_tex
    excerpt: \mathcal M_\sigma(r)=\xi_R\bigl[\mathcal R_\infty-\mathcal R(r)\bigr],
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
    line: 1
    role: notes_stage
    excerpt: '# Moving-Throat PDE — Stage 247: Relaxed Stationary Barrier Compiler from One-Port Short-Range, Leakage, `U/V`
      Drain, and Compensated Source Packets'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 17
    role: notes_stage
    excerpt: 2. the exact Stage-244 leakage / work compiler,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 33
    role: notes_stage
    excerpt: (S_{\rm leak},\,\mathcal W_w),
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 48
    role: notes_stage
    excerpt: 'It has four jobs:'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 51
    role: notes_stage
    excerpt: 2. import the Stage-244 leakage / work packet into barrier units,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 61
    role: notes_stage
    excerpt: '- **Stage 244**, which attached leakage and scalar-photon work to the selected-support packet,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 63
    role: notes_stage
    excerpt: '- **Stage 246**, which converted the sign-changing source branch into explicit mouth/source observables,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 98
    role: notes_stage
    excerpt: Q=G_U^2\Omega_W^2+2G_UG_WR_{\rm mix}+G_W^2\Omega_U^2,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 101
    role: notes_stage
    excerpt: P=\Omega_U^2G_W+R_{\rm mix}G_U,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 116
    role: notes_stage
    excerpt: \chi_{UU}=\frac{K_*\Omega_W^2-G_W^2}{\Delta D_0},
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 118
    role: notes_stage
    excerpt: \chi_{UW}=\frac{K_*R_{\rm mix}+G_UG_W}{\Delta D_0},
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 124
    role: notes_stage
    excerpt: P_U=G_U\Omega_W^2+R_{\rm mix}G_W.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 129
    role: notes_stage
    excerpt: \beta_Q,\qquad \beta_U,\qquad \beta_W,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 133
    role: notes_stage
    excerpt: \mathcal C_6=\chi_{qq}\beta_Q^2,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 136
    role: notes_stage
    excerpt: \mathcal C_4=\chi_{qU}\beta_Q\beta_U+\chi_{qW}\beta_Q\beta_W,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 139
    role: notes_stage
    excerpt: \mathcal C_2=\chi_{UU}\beta_U^2+2\chi_{UW}\beta_U\beta_W+\chi_{WW}\beta_W^2.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 187
    role: notes_stage
    excerpt: '### 2.1 Stage-244 leakage and reduced Session-I work scalar'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 196
    role: notes_stage
    excerpt: S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 198
    role: notes_stage
    excerpt: \frac{8\sqrt2\,\eta_{\rm leak}\mu_w\rho_0}{\pi^{5/2}\lambda^3}
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 206
    role: notes_stage
    excerpt: \frac{512\,\eta_{\rm leak}^2\mu_w q\rho_0}{\pi^4\lambda^2}
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 225
    role: notes_stage
    excerpt: \Delta E_{UV}(r):=\eta_{UV}\mathcal D_{UV}(r),
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 228
    role: notes_stage
    excerpt: where `\eta_{UV}` is the same Session-I barrier-lowering weight carried in the stationary script.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 247
    role: notes_stage
    excerpt: \mathcal M_\sigma(r):=\xi_R\bigl[\mathcal R_\infty-\mathcal R(r)\bigr],
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 296
    role: notes_stage
    excerpt: \lambda_L S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 308
    role: notes_stage
    excerpt: '- `\lambda_L\ge 0` converts the open-subsystem leakage source into stationary barrier units,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 320
    role: notes_stage
    excerpt: \lambda_L S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 344
    role: notes_stage
    excerpt: '- `(S_{\rm leak},\mathcal W_w^{\rm sess})` are support-side / open-system objects,'
  origin_claims:
  - parameter: Veff_obs
    introduced_at_stage: 247
    introduced_at_line: 372
    citation:
      path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
      line: 372
      excerpt: V_{\rm eff}(r_{\rm soft})=1.74701126,
  constraints:
  - parameter: Veff_obs
    constraint_kind: published_target
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
      line: 367
      excerpt: Appendix B.2 of the barrier-session write-up records the strongest stationary softening point at
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings:
  - type: published_target
    severity: high
    summary: 'Veff_obs=1.74701126 is the recorded Session-I strongest-softening-point figure imported from Appendix B.2 of
      the barrier write-up (notes line 367-372; script line 213 Veff_obs = sp.Float(''1.74701126'')). It is THE target of
      the stage: lambda_L is back-solved so the compiler reproduces it (script line 238), and the forward ''falsifiable closure''
      (line 250) checks the same number, making that check tautological for lambda_L. This recorded benchmark figure is the
      load-bearing fit target for the whole stage-247 benchmark; high severity because it is the anchor the entire decomposition
      is tuned to.'
    citations:
    - path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
      line: 213
      excerpt: Veff_obs = sp.Float("1.74701126")
    - path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
      line: 238
      excerpt: lambda_L_soft = sp.N((Vshort_num - Wsess_obs - UVdrop_obs - M_sigma_num - Veff_obs) / S_soft, 16)
    - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
      line: 367
      excerpt: Appendix B.2 of the barrier-session write-up records the strongest stationary softening point at
  downstream_dependents:
  - '248'
  synthesis_ingested_at: '2026-06-11T19:37:44Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_27/fit_stage247_session1_packet__veff_obs.yaml
- schema_version: 1
  candidate_id: fit_stage247_session1_packet
  parameter_name: Wsess_obs
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage247_session1_packet__phase_b__wsess_obs.md
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
    line: 4
    role: paper_stage_tex
    excerpt: '\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py}.  Mathematica
      audit: none yet.}'
  - path: paper/stages/stage_247.tex
    line: 12
    role: paper_stage_tex
    excerpt: Q=G_U^2\Omega_W^2+2G_UG_WR_{\rm mix}+G_W^2\Omega_U^2,
  - path: paper/stages/stage_247.tex
    line: 16
    role: paper_stage_tex
    excerpt: P=\Omega_U^2G_W+R_{\rm mix}G_U,
  - path: paper/stages/stage_247.tex
    line: 20
    role: paper_stage_tex
    excerpt: With primitive short-range source amplitudes \(\beta_Q,\beta_U,\beta_W\), the product-family coefficients are
  - path: paper/stages/stage_247.tex
    line: 23
    role: paper_stage_tex
    excerpt: \mathcal C_6=\chi_{qq}\beta_Q^2,
  - path: paper/stages/stage_247.tex
    line: 25
    role: paper_stage_tex
    excerpt: \mathcal C_4=\chi_{qU}\beta_Q\beta_U+\chi_{qW}\beta_Q\beta_W,
  - path: paper/stages/stage_247.tex
    line: 29
    role: paper_stage_tex
    excerpt: \mathcal C_2=\chi_{UU}\beta_U^2+2\chi_{UW}\beta_U\beta_W+\chi_{WW}\beta_W^2.
  - path: paper/stages/stage_247.tex
    line: 60
    role: paper_stage_tex
    excerpt: S_{\rm leak}=\frac{8\sqrt2\eta_{\rm leak}\mu_w\rho_0}{\pi^{5/2}\lambda^3}\mathfrak L,
  - path: paper/stages/stage_247.tex
    line: 66
    role: paper_stage_tex
    excerpt: \frac{512\eta_{\rm leak}^{2}\mu_wq\rho_0}{\pi^4\lambda^2}\mathfrak L^2,
  - path: paper/stages/stage_247.tex
    line: 68
    role: paper_stage_tex
    excerpt: \Delta E_{UV}=\eta_{UV}\mathcal D_{UV},
  - path: paper/stages/stage_247.tex
    line: 72
    role: paper_stage_tex
    excerpt: \mathcal M_\sigma(r)=\xi_R\bigl[\mathcal R_\infty-\mathcal R(r)\bigr],
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
    line: 1
    role: notes_stage
    excerpt: '# Moving-Throat PDE — Stage 247: Relaxed Stationary Barrier Compiler from One-Port Short-Range, Leakage, `U/V`
      Drain, and Compensated Source Packets'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 17
    role: notes_stage
    excerpt: 2. the exact Stage-244 leakage / work compiler,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 33
    role: notes_stage
    excerpt: (S_{\rm leak},\,\mathcal W_w),
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 48
    role: notes_stage
    excerpt: 'It has four jobs:'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 51
    role: notes_stage
    excerpt: 2. import the Stage-244 leakage / work packet into barrier units,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 61
    role: notes_stage
    excerpt: '- **Stage 244**, which attached leakage and scalar-photon work to the selected-support packet,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 63
    role: notes_stage
    excerpt: '- **Stage 246**, which converted the sign-changing source branch into explicit mouth/source observables,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 98
    role: notes_stage
    excerpt: Q=G_U^2\Omega_W^2+2G_UG_WR_{\rm mix}+G_W^2\Omega_U^2,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 101
    role: notes_stage
    excerpt: P=\Omega_U^2G_W+R_{\rm mix}G_U,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 116
    role: notes_stage
    excerpt: \chi_{UU}=\frac{K_*\Omega_W^2-G_W^2}{\Delta D_0},
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 118
    role: notes_stage
    excerpt: \chi_{UW}=\frac{K_*R_{\rm mix}+G_UG_W}{\Delta D_0},
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 124
    role: notes_stage
    excerpt: P_U=G_U\Omega_W^2+R_{\rm mix}G_W.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 129
    role: notes_stage
    excerpt: \beta_Q,\qquad \beta_U,\qquad \beta_W,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 133
    role: notes_stage
    excerpt: \mathcal C_6=\chi_{qq}\beta_Q^2,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 136
    role: notes_stage
    excerpt: \mathcal C_4=\chi_{qU}\beta_Q\beta_U+\chi_{qW}\beta_Q\beta_W,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 139
    role: notes_stage
    excerpt: \mathcal C_2=\chi_{UU}\beta_U^2+2\chi_{UW}\beta_U\beta_W+\chi_{WW}\beta_W^2.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 187
    role: notes_stage
    excerpt: '### 2.1 Stage-244 leakage and reduced Session-I work scalar'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 196
    role: notes_stage
    excerpt: S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 198
    role: notes_stage
    excerpt: \frac{8\sqrt2\,\eta_{\rm leak}\mu_w\rho_0}{\pi^{5/2}\lambda^3}
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 206
    role: notes_stage
    excerpt: \frac{512\,\eta_{\rm leak}^2\mu_w q\rho_0}{\pi^4\lambda^2}
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 225
    role: notes_stage
    excerpt: \Delta E_{UV}(r):=\eta_{UV}\mathcal D_{UV}(r),
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 228
    role: notes_stage
    excerpt: where `\eta_{UV}` is the same Session-I barrier-lowering weight carried in the stationary script.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 247
    role: notes_stage
    excerpt: \mathcal M_\sigma(r):=\xi_R\bigl[\mathcal R_\infty-\mathcal R(r)\bigr],
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 296
    role: notes_stage
    excerpt: \lambda_L S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 308
    role: notes_stage
    excerpt: '- `\lambda_L\ge 0` converts the open-subsystem leakage source into stationary barrier units,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 320
    role: notes_stage
    excerpt: \lambda_L S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 344
    role: notes_stage
    excerpt: '- `(S_{\rm leak},\mathcal W_w^{\rm sess})` are support-side / open-system objects,'
  origin_claims:
  - parameter: Wsess_obs
    introduced_at_stage: 247
    introduced_at_line: 429
    citation:
      path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
      line: 429
      excerpt: \mathcal W_w^{\rm sess}(r_{\rm soft})=1.51632107,
  constraints:
  - parameter: Wsess_obs
    constraint_kind: published_target
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
      line: 427
      excerpt: The Session-I table also records
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings:
  - type: published_target
    severity: needs_triage
    summary: Wsess_obs=1.51632107 is a recorded Session-I benchmark figure read from the tabulated stationary run (notes line
      427 'The Session-I table also records'; script line 214 Wsess_obs = sp.Float('1.51632107')). It is a hardcoded recorded
      target, not a value derived within stage 247 from the work-law formula at this slice; the script's F2 check inverts
      it to Lvar=20.01677473 and re-substitutes (a consistency check against the recorded number, not an independent derivation).
      It is one of the recorded targets that, together with Veff_obs, fixes lambda_L.
    citations:
    - path: scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py
      line: 214
      excerpt: Wsess_obs = sp.Float("1.51632107")
    - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
      line: 427
      excerpt: The Session-I table also records
  downstream_dependents:
  - '248'
  synthesis_ingested_at: '2026-06-11T19:37:44Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_27/fit_stage247_session1_packet__wsess_obs.yaml
- schema_version: 1
  candidate_id: fit_stage247_session1_packet
  parameter_name: beta_Q
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage247_session1_packet__phase_b__beta_q.md
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
    line: 4
    role: paper_stage_tex
    excerpt: '\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py}.  Mathematica
      audit: none yet.}'
  - path: paper/stages/stage_247.tex
    line: 12
    role: paper_stage_tex
    excerpt: Q=G_U^2\Omega_W^2+2G_UG_WR_{\rm mix}+G_W^2\Omega_U^2,
  - path: paper/stages/stage_247.tex
    line: 16
    role: paper_stage_tex
    excerpt: P=\Omega_U^2G_W+R_{\rm mix}G_U,
  - path: paper/stages/stage_247.tex
    line: 20
    role: paper_stage_tex
    excerpt: With primitive short-range source amplitudes \(\beta_Q,\beta_U,\beta_W\), the product-family coefficients are
  - path: paper/stages/stage_247.tex
    line: 23
    role: paper_stage_tex
    excerpt: \mathcal C_6=\chi_{qq}\beta_Q^2,
  - path: paper/stages/stage_247.tex
    line: 25
    role: paper_stage_tex
    excerpt: \mathcal C_4=\chi_{qU}\beta_Q\beta_U+\chi_{qW}\beta_Q\beta_W,
  - path: paper/stages/stage_247.tex
    line: 29
    role: paper_stage_tex
    excerpt: \mathcal C_2=\chi_{UU}\beta_U^2+2\chi_{UW}\beta_U\beta_W+\chi_{WW}\beta_W^2.
  - path: paper/stages/stage_247.tex
    line: 60
    role: paper_stage_tex
    excerpt: S_{\rm leak}=\frac{8\sqrt2\eta_{\rm leak}\mu_w\rho_0}{\pi^{5/2}\lambda^3}\mathfrak L,
  - path: paper/stages/stage_247.tex
    line: 66
    role: paper_stage_tex
    excerpt: \frac{512\eta_{\rm leak}^{2}\mu_wq\rho_0}{\pi^4\lambda^2}\mathfrak L^2,
  - path: paper/stages/stage_247.tex
    line: 68
    role: paper_stage_tex
    excerpt: \Delta E_{UV}=\eta_{UV}\mathcal D_{UV},
  - path: paper/stages/stage_247.tex
    line: 72
    role: paper_stage_tex
    excerpt: \mathcal M_\sigma(r)=\xi_R\bigl[\mathcal R_\infty-\mathcal R(r)\bigr],
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
    line: 1
    role: notes_stage
    excerpt: '# Moving-Throat PDE — Stage 247: Relaxed Stationary Barrier Compiler from One-Port Short-Range, Leakage, `U/V`
      Drain, and Compensated Source Packets'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 17
    role: notes_stage
    excerpt: 2. the exact Stage-244 leakage / work compiler,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 33
    role: notes_stage
    excerpt: (S_{\rm leak},\,\mathcal W_w),
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 48
    role: notes_stage
    excerpt: 'It has four jobs:'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 51
    role: notes_stage
    excerpt: 2. import the Stage-244 leakage / work packet into barrier units,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 61
    role: notes_stage
    excerpt: '- **Stage 244**, which attached leakage and scalar-photon work to the selected-support packet,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 63
    role: notes_stage
    excerpt: '- **Stage 246**, which converted the sign-changing source branch into explicit mouth/source observables,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 98
    role: notes_stage
    excerpt: Q=G_U^2\Omega_W^2+2G_UG_WR_{\rm mix}+G_W^2\Omega_U^2,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 101
    role: notes_stage
    excerpt: P=\Omega_U^2G_W+R_{\rm mix}G_U,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 116
    role: notes_stage
    excerpt: \chi_{UU}=\frac{K_*\Omega_W^2-G_W^2}{\Delta D_0},
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 118
    role: notes_stage
    excerpt: \chi_{UW}=\frac{K_*R_{\rm mix}+G_UG_W}{\Delta D_0},
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 124
    role: notes_stage
    excerpt: P_U=G_U\Omega_W^2+R_{\rm mix}G_W.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 129
    role: notes_stage
    excerpt: \beta_Q,\qquad \beta_U,\qquad \beta_W,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 133
    role: notes_stage
    excerpt: \mathcal C_6=\chi_{qq}\beta_Q^2,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 136
    role: notes_stage
    excerpt: \mathcal C_4=\chi_{qU}\beta_Q\beta_U+\chi_{qW}\beta_Q\beta_W,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 139
    role: notes_stage
    excerpt: \mathcal C_2=\chi_{UU}\beta_U^2+2\chi_{UW}\beta_U\beta_W+\chi_{WW}\beta_W^2.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 187
    role: notes_stage
    excerpt: '### 2.1 Stage-244 leakage and reduced Session-I work scalar'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 196
    role: notes_stage
    excerpt: S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 198
    role: notes_stage
    excerpt: \frac{8\sqrt2\,\eta_{\rm leak}\mu_w\rho_0}{\pi^{5/2}\lambda^3}
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 206
    role: notes_stage
    excerpt: \frac{512\,\eta_{\rm leak}^2\mu_w q\rho_0}{\pi^4\lambda^2}
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 225
    role: notes_stage
    excerpt: \Delta E_{UV}(r):=\eta_{UV}\mathcal D_{UV}(r),
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 228
    role: notes_stage
    excerpt: where `\eta_{UV}` is the same Session-I barrier-lowering weight carried in the stationary script.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 247
    role: notes_stage
    excerpt: \mathcal M_\sigma(r):=\xi_R\bigl[\mathcal R_\infty-\mathcal R(r)\bigr],
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 296
    role: notes_stage
    excerpt: \lambda_L S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 308
    role: notes_stage
    excerpt: '- `\lambda_L\ge 0` converts the open-subsystem leakage source into stationary barrier units,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 320
    role: notes_stage
    excerpt: \lambda_L S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 344
    role: notes_stage
    excerpt: '- `(S_{\rm leak},\mathcal W_w^{\rm sess})` are support-side / open-system objects,'
  origin_claims:
  - parameter: beta_Q
    introduced_at_stage: 247
    introduced_at_line: 396
    citation:
      path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
      line: 396
      excerpt: \beta_Q=0.03,
  constraints:
  - parameter: beta_Q
    constraint_kind: free_choice
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
      line: 381
      excerpt: Using the recorded one-port parameters
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings: []
  downstream_dependents:
  - '248'
  synthesis_ingested_at: '2026-06-11T19:37:44Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_27/fit_stage247_session1_packet__beta_q.yaml
- schema_version: 1
  candidate_id: fit_stage247_session1_packet
  parameter_name: beta_U
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage247_session1_packet__phase_b__beta_u.md
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
    line: 4
    role: paper_stage_tex
    excerpt: '\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py}.  Mathematica
      audit: none yet.}'
  - path: paper/stages/stage_247.tex
    line: 12
    role: paper_stage_tex
    excerpt: Q=G_U^2\Omega_W^2+2G_UG_WR_{\rm mix}+G_W^2\Omega_U^2,
  - path: paper/stages/stage_247.tex
    line: 16
    role: paper_stage_tex
    excerpt: P=\Omega_U^2G_W+R_{\rm mix}G_U,
  - path: paper/stages/stage_247.tex
    line: 20
    role: paper_stage_tex
    excerpt: With primitive short-range source amplitudes \(\beta_Q,\beta_U,\beta_W\), the product-family coefficients are
  - path: paper/stages/stage_247.tex
    line: 23
    role: paper_stage_tex
    excerpt: \mathcal C_6=\chi_{qq}\beta_Q^2,
  - path: paper/stages/stage_247.tex
    line: 25
    role: paper_stage_tex
    excerpt: \mathcal C_4=\chi_{qU}\beta_Q\beta_U+\chi_{qW}\beta_Q\beta_W,
  - path: paper/stages/stage_247.tex
    line: 29
    role: paper_stage_tex
    excerpt: \mathcal C_2=\chi_{UU}\beta_U^2+2\chi_{UW}\beta_U\beta_W+\chi_{WW}\beta_W^2.
  - path: paper/stages/stage_247.tex
    line: 60
    role: paper_stage_tex
    excerpt: S_{\rm leak}=\frac{8\sqrt2\eta_{\rm leak}\mu_w\rho_0}{\pi^{5/2}\lambda^3}\mathfrak L,
  - path: paper/stages/stage_247.tex
    line: 66
    role: paper_stage_tex
    excerpt: \frac{512\eta_{\rm leak}^{2}\mu_wq\rho_0}{\pi^4\lambda^2}\mathfrak L^2,
  - path: paper/stages/stage_247.tex
    line: 68
    role: paper_stage_tex
    excerpt: \Delta E_{UV}=\eta_{UV}\mathcal D_{UV},
  - path: paper/stages/stage_247.tex
    line: 72
    role: paper_stage_tex
    excerpt: \mathcal M_\sigma(r)=\xi_R\bigl[\mathcal R_\infty-\mathcal R(r)\bigr],
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
    line: 1
    role: notes_stage
    excerpt: '# Moving-Throat PDE — Stage 247: Relaxed Stationary Barrier Compiler from One-Port Short-Range, Leakage, `U/V`
      Drain, and Compensated Source Packets'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 17
    role: notes_stage
    excerpt: 2. the exact Stage-244 leakage / work compiler,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 33
    role: notes_stage
    excerpt: (S_{\rm leak},\,\mathcal W_w),
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 48
    role: notes_stage
    excerpt: 'It has four jobs:'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 51
    role: notes_stage
    excerpt: 2. import the Stage-244 leakage / work packet into barrier units,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 61
    role: notes_stage
    excerpt: '- **Stage 244**, which attached leakage and scalar-photon work to the selected-support packet,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 63
    role: notes_stage
    excerpt: '- **Stage 246**, which converted the sign-changing source branch into explicit mouth/source observables,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 98
    role: notes_stage
    excerpt: Q=G_U^2\Omega_W^2+2G_UG_WR_{\rm mix}+G_W^2\Omega_U^2,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 101
    role: notes_stage
    excerpt: P=\Omega_U^2G_W+R_{\rm mix}G_U,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 116
    role: notes_stage
    excerpt: \chi_{UU}=\frac{K_*\Omega_W^2-G_W^2}{\Delta D_0},
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 118
    role: notes_stage
    excerpt: \chi_{UW}=\frac{K_*R_{\rm mix}+G_UG_W}{\Delta D_0},
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 124
    role: notes_stage
    excerpt: P_U=G_U\Omega_W^2+R_{\rm mix}G_W.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 129
    role: notes_stage
    excerpt: \beta_Q,\qquad \beta_U,\qquad \beta_W,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 133
    role: notes_stage
    excerpt: \mathcal C_6=\chi_{qq}\beta_Q^2,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 136
    role: notes_stage
    excerpt: \mathcal C_4=\chi_{qU}\beta_Q\beta_U+\chi_{qW}\beta_Q\beta_W,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 139
    role: notes_stage
    excerpt: \mathcal C_2=\chi_{UU}\beta_U^2+2\chi_{UW}\beta_U\beta_W+\chi_{WW}\beta_W^2.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 187
    role: notes_stage
    excerpt: '### 2.1 Stage-244 leakage and reduced Session-I work scalar'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 196
    role: notes_stage
    excerpt: S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 198
    role: notes_stage
    excerpt: \frac{8\sqrt2\,\eta_{\rm leak}\mu_w\rho_0}{\pi^{5/2}\lambda^3}
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 206
    role: notes_stage
    excerpt: \frac{512\,\eta_{\rm leak}^2\mu_w q\rho_0}{\pi^4\lambda^2}
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 225
    role: notes_stage
    excerpt: \Delta E_{UV}(r):=\eta_{UV}\mathcal D_{UV}(r),
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 228
    role: notes_stage
    excerpt: where `\eta_{UV}` is the same Session-I barrier-lowering weight carried in the stationary script.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 247
    role: notes_stage
    excerpt: \mathcal M_\sigma(r):=\xi_R\bigl[\mathcal R_\infty-\mathcal R(r)\bigr],
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 296
    role: notes_stage
    excerpt: \lambda_L S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 308
    role: notes_stage
    excerpt: '- `\lambda_L\ge 0` converts the open-subsystem leakage source into stationary barrier units,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 320
    role: notes_stage
    excerpt: \lambda_L S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 344
    role: notes_stage
    excerpt: '- `(S_{\rm leak},\mathcal W_w^{\rm sess})` are support-side / open-system objects,'
  origin_claims:
  - parameter: beta_U
    introduced_at_stage: 247
    introduced_at_line: 398
    citation:
      path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
      line: 398
      excerpt: \beta_{U0}=0.15,
  constraints:
  - parameter: beta_U
    constraint_kind: free_choice
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
      line: 381
      excerpt: Using the recorded one-port parameters
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings: []
  downstream_dependents:
  - '248'
  synthesis_ingested_at: '2026-06-11T19:37:44Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_27/fit_stage247_session1_packet__beta_u.yaml
- schema_version: 1
  candidate_id: fit_stage247_session1_packet
  parameter_name: beta_W
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage247_session1_packet__phase_b__beta_w.md
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
    line: 4
    role: paper_stage_tex
    excerpt: '\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py}.  Mathematica
      audit: none yet.}'
  - path: paper/stages/stage_247.tex
    line: 12
    role: paper_stage_tex
    excerpt: Q=G_U^2\Omega_W^2+2G_UG_WR_{\rm mix}+G_W^2\Omega_U^2,
  - path: paper/stages/stage_247.tex
    line: 16
    role: paper_stage_tex
    excerpt: P=\Omega_U^2G_W+R_{\rm mix}G_U,
  - path: paper/stages/stage_247.tex
    line: 20
    role: paper_stage_tex
    excerpt: With primitive short-range source amplitudes \(\beta_Q,\beta_U,\beta_W\), the product-family coefficients are
  - path: paper/stages/stage_247.tex
    line: 23
    role: paper_stage_tex
    excerpt: \mathcal C_6=\chi_{qq}\beta_Q^2,
  - path: paper/stages/stage_247.tex
    line: 25
    role: paper_stage_tex
    excerpt: \mathcal C_4=\chi_{qU}\beta_Q\beta_U+\chi_{qW}\beta_Q\beta_W,
  - path: paper/stages/stage_247.tex
    line: 29
    role: paper_stage_tex
    excerpt: \mathcal C_2=\chi_{UU}\beta_U^2+2\chi_{UW}\beta_U\beta_W+\chi_{WW}\beta_W^2.
  - path: paper/stages/stage_247.tex
    line: 60
    role: paper_stage_tex
    excerpt: S_{\rm leak}=\frac{8\sqrt2\eta_{\rm leak}\mu_w\rho_0}{\pi^{5/2}\lambda^3}\mathfrak L,
  - path: paper/stages/stage_247.tex
    line: 66
    role: paper_stage_tex
    excerpt: \frac{512\eta_{\rm leak}^{2}\mu_wq\rho_0}{\pi^4\lambda^2}\mathfrak L^2,
  - path: paper/stages/stage_247.tex
    line: 68
    role: paper_stage_tex
    excerpt: \Delta E_{UV}=\eta_{UV}\mathcal D_{UV},
  - path: paper/stages/stage_247.tex
    line: 72
    role: paper_stage_tex
    excerpt: \mathcal M_\sigma(r)=\xi_R\bigl[\mathcal R_\infty-\mathcal R(r)\bigr],
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
    line: 1
    role: notes_stage
    excerpt: '# Moving-Throat PDE — Stage 247: Relaxed Stationary Barrier Compiler from One-Port Short-Range, Leakage, `U/V`
      Drain, and Compensated Source Packets'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 17
    role: notes_stage
    excerpt: 2. the exact Stage-244 leakage / work compiler,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 33
    role: notes_stage
    excerpt: (S_{\rm leak},\,\mathcal W_w),
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 48
    role: notes_stage
    excerpt: 'It has four jobs:'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 51
    role: notes_stage
    excerpt: 2. import the Stage-244 leakage / work packet into barrier units,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 61
    role: notes_stage
    excerpt: '- **Stage 244**, which attached leakage and scalar-photon work to the selected-support packet,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 63
    role: notes_stage
    excerpt: '- **Stage 246**, which converted the sign-changing source branch into explicit mouth/source observables,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 98
    role: notes_stage
    excerpt: Q=G_U^2\Omega_W^2+2G_UG_WR_{\rm mix}+G_W^2\Omega_U^2,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 101
    role: notes_stage
    excerpt: P=\Omega_U^2G_W+R_{\rm mix}G_U,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 116
    role: notes_stage
    excerpt: \chi_{UU}=\frac{K_*\Omega_W^2-G_W^2}{\Delta D_0},
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 118
    role: notes_stage
    excerpt: \chi_{UW}=\frac{K_*R_{\rm mix}+G_UG_W}{\Delta D_0},
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 124
    role: notes_stage
    excerpt: P_U=G_U\Omega_W^2+R_{\rm mix}G_W.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 129
    role: notes_stage
    excerpt: \beta_Q,\qquad \beta_U,\qquad \beta_W,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 133
    role: notes_stage
    excerpt: \mathcal C_6=\chi_{qq}\beta_Q^2,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 136
    role: notes_stage
    excerpt: \mathcal C_4=\chi_{qU}\beta_Q\beta_U+\chi_{qW}\beta_Q\beta_W,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 139
    role: notes_stage
    excerpt: \mathcal C_2=\chi_{UU}\beta_U^2+2\chi_{UW}\beta_U\beta_W+\chi_{WW}\beta_W^2.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 187
    role: notes_stage
    excerpt: '### 2.1 Stage-244 leakage and reduced Session-I work scalar'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 196
    role: notes_stage
    excerpt: S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 198
    role: notes_stage
    excerpt: \frac{8\sqrt2\,\eta_{\rm leak}\mu_w\rho_0}{\pi^{5/2}\lambda^3}
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 206
    role: notes_stage
    excerpt: \frac{512\,\eta_{\rm leak}^2\mu_w q\rho_0}{\pi^4\lambda^2}
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 225
    role: notes_stage
    excerpt: \Delta E_{UV}(r):=\eta_{UV}\mathcal D_{UV}(r),
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 228
    role: notes_stage
    excerpt: where `\eta_{UV}` is the same Session-I barrier-lowering weight carried in the stationary script.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 247
    role: notes_stage
    excerpt: \mathcal M_\sigma(r):=\xi_R\bigl[\mathcal R_\infty-\mathcal R(r)\bigr],
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 296
    role: notes_stage
    excerpt: \lambda_L S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 308
    role: notes_stage
    excerpt: '- `\lambda_L\ge 0` converts the open-subsystem leakage source into stationary barrier units,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 320
    role: notes_stage
    excerpt: \lambda_L S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 344
    role: notes_stage
    excerpt: '- `(S_{\rm leak},\mathcal W_w^{\rm sess})` are support-side / open-system objects,'
  origin_claims:
  - parameter: beta_W
    introduced_at_stage: 247
    introduced_at_line: 400
    citation:
      path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
      line: 400
      excerpt: \beta_{W0}=0.20,
  constraints:
  - parameter: beta_W
    constraint_kind: free_choice
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
      line: 381
      excerpt: Using the recorded one-port parameters
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings: []
  downstream_dependents:
  - '248'
  synthesis_ingested_at: '2026-06-11T19:37:44Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_27/fit_stage247_session1_packet__beta_w.yaml
- schema_version: 1
  candidate_id: fit_stage247_session1_packet
  parameter_name: eta_leak
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage247_session1_packet__phase_b__eta_leak.md
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
    line: 4
    role: paper_stage_tex
    excerpt: '\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py}.  Mathematica
      audit: none yet.}'
  - path: paper/stages/stage_247.tex
    line: 12
    role: paper_stage_tex
    excerpt: Q=G_U^2\Omega_W^2+2G_UG_WR_{\rm mix}+G_W^2\Omega_U^2,
  - path: paper/stages/stage_247.tex
    line: 16
    role: paper_stage_tex
    excerpt: P=\Omega_U^2G_W+R_{\rm mix}G_U,
  - path: paper/stages/stage_247.tex
    line: 20
    role: paper_stage_tex
    excerpt: With primitive short-range source amplitudes \(\beta_Q,\beta_U,\beta_W\), the product-family coefficients are
  - path: paper/stages/stage_247.tex
    line: 23
    role: paper_stage_tex
    excerpt: \mathcal C_6=\chi_{qq}\beta_Q^2,
  - path: paper/stages/stage_247.tex
    line: 25
    role: paper_stage_tex
    excerpt: \mathcal C_4=\chi_{qU}\beta_Q\beta_U+\chi_{qW}\beta_Q\beta_W,
  - path: paper/stages/stage_247.tex
    line: 29
    role: paper_stage_tex
    excerpt: \mathcal C_2=\chi_{UU}\beta_U^2+2\chi_{UW}\beta_U\beta_W+\chi_{WW}\beta_W^2.
  - path: paper/stages/stage_247.tex
    line: 60
    role: paper_stage_tex
    excerpt: S_{\rm leak}=\frac{8\sqrt2\eta_{\rm leak}\mu_w\rho_0}{\pi^{5/2}\lambda^3}\mathfrak L,
  - path: paper/stages/stage_247.tex
    line: 66
    role: paper_stage_tex
    excerpt: \frac{512\eta_{\rm leak}^{2}\mu_wq\rho_0}{\pi^4\lambda^2}\mathfrak L^2,
  - path: paper/stages/stage_247.tex
    line: 68
    role: paper_stage_tex
    excerpt: \Delta E_{UV}=\eta_{UV}\mathcal D_{UV},
  - path: paper/stages/stage_247.tex
    line: 72
    role: paper_stage_tex
    excerpt: \mathcal M_\sigma(r)=\xi_R\bigl[\mathcal R_\infty-\mathcal R(r)\bigr],
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
    line: 1
    role: notes_stage
    excerpt: '# Moving-Throat PDE — Stage 247: Relaxed Stationary Barrier Compiler from One-Port Short-Range, Leakage, `U/V`
      Drain, and Compensated Source Packets'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 17
    role: notes_stage
    excerpt: 2. the exact Stage-244 leakage / work compiler,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 33
    role: notes_stage
    excerpt: (S_{\rm leak},\,\mathcal W_w),
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 48
    role: notes_stage
    excerpt: 'It has four jobs:'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 51
    role: notes_stage
    excerpt: 2. import the Stage-244 leakage / work packet into barrier units,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 61
    role: notes_stage
    excerpt: '- **Stage 244**, which attached leakage and scalar-photon work to the selected-support packet,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 63
    role: notes_stage
    excerpt: '- **Stage 246**, which converted the sign-changing source branch into explicit mouth/source observables,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 98
    role: notes_stage
    excerpt: Q=G_U^2\Omega_W^2+2G_UG_WR_{\rm mix}+G_W^2\Omega_U^2,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 101
    role: notes_stage
    excerpt: P=\Omega_U^2G_W+R_{\rm mix}G_U,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 116
    role: notes_stage
    excerpt: \chi_{UU}=\frac{K_*\Omega_W^2-G_W^2}{\Delta D_0},
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 118
    role: notes_stage
    excerpt: \chi_{UW}=\frac{K_*R_{\rm mix}+G_UG_W}{\Delta D_0},
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 124
    role: notes_stage
    excerpt: P_U=G_U\Omega_W^2+R_{\rm mix}G_W.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 129
    role: notes_stage
    excerpt: \beta_Q,\qquad \beta_U,\qquad \beta_W,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 133
    role: notes_stage
    excerpt: \mathcal C_6=\chi_{qq}\beta_Q^2,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 136
    role: notes_stage
    excerpt: \mathcal C_4=\chi_{qU}\beta_Q\beta_U+\chi_{qW}\beta_Q\beta_W,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 139
    role: notes_stage
    excerpt: \mathcal C_2=\chi_{UU}\beta_U^2+2\chi_{UW}\beta_U\beta_W+\chi_{WW}\beta_W^2.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 187
    role: notes_stage
    excerpt: '### 2.1 Stage-244 leakage and reduced Session-I work scalar'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 196
    role: notes_stage
    excerpt: S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 198
    role: notes_stage
    excerpt: \frac{8\sqrt2\,\eta_{\rm leak}\mu_w\rho_0}{\pi^{5/2}\lambda^3}
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 206
    role: notes_stage
    excerpt: \frac{512\,\eta_{\rm leak}^2\mu_w q\rho_0}{\pi^4\lambda^2}
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 225
    role: notes_stage
    excerpt: \Delta E_{UV}(r):=\eta_{UV}\mathcal D_{UV}(r),
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 228
    role: notes_stage
    excerpt: where `\eta_{UV}` is the same Session-I barrier-lowering weight carried in the stationary script.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 247
    role: notes_stage
    excerpt: \mathcal M_\sigma(r):=\xi_R\bigl[\mathcal R_\infty-\mathcal R(r)\bigr],
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 296
    role: notes_stage
    excerpt: \lambda_L S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 308
    role: notes_stage
    excerpt: '- `\lambda_L\ge 0` converts the open-subsystem leakage source into stationary barrier units,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 320
    role: notes_stage
    excerpt: \lambda_L S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 344
    role: notes_stage
    excerpt: '- `(S_{\rm leak},\mathcal W_w^{\rm sess})` are support-side / open-system objects,'
  origin_claims:
  - parameter: eta_leak
    introduced_at_stage: 247
    introduced_at_line: 445
    citation:
      path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
      line: 445
      excerpt: \eta_{\rm leak}=0.03,
  constraints:
  - parameter: eta_leak
    constraint_kind: free_choice
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
      line: 441
      excerpt: with the Session-I parameters
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings: []
  downstream_dependents:
  - '248'
  synthesis_ingested_at: '2026-06-11T19:37:44Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_27/fit_stage247_session1_packet__eta_leak.yaml
- schema_version: 1
  candidate_id: fit_stage247_session1_packet
  parameter_name: mu_w
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage247_session1_packet__phase_b__mu_w.md
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
    line: 4
    role: paper_stage_tex
    excerpt: '\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py}.  Mathematica
      audit: none yet.}'
  - path: paper/stages/stage_247.tex
    line: 12
    role: paper_stage_tex
    excerpt: Q=G_U^2\Omega_W^2+2G_UG_WR_{\rm mix}+G_W^2\Omega_U^2,
  - path: paper/stages/stage_247.tex
    line: 16
    role: paper_stage_tex
    excerpt: P=\Omega_U^2G_W+R_{\rm mix}G_U,
  - path: paper/stages/stage_247.tex
    line: 20
    role: paper_stage_tex
    excerpt: With primitive short-range source amplitudes \(\beta_Q,\beta_U,\beta_W\), the product-family coefficients are
  - path: paper/stages/stage_247.tex
    line: 23
    role: paper_stage_tex
    excerpt: \mathcal C_6=\chi_{qq}\beta_Q^2,
  - path: paper/stages/stage_247.tex
    line: 25
    role: paper_stage_tex
    excerpt: \mathcal C_4=\chi_{qU}\beta_Q\beta_U+\chi_{qW}\beta_Q\beta_W,
  - path: paper/stages/stage_247.tex
    line: 29
    role: paper_stage_tex
    excerpt: \mathcal C_2=\chi_{UU}\beta_U^2+2\chi_{UW}\beta_U\beta_W+\chi_{WW}\beta_W^2.
  - path: paper/stages/stage_247.tex
    line: 60
    role: paper_stage_tex
    excerpt: S_{\rm leak}=\frac{8\sqrt2\eta_{\rm leak}\mu_w\rho_0}{\pi^{5/2}\lambda^3}\mathfrak L,
  - path: paper/stages/stage_247.tex
    line: 66
    role: paper_stage_tex
    excerpt: \frac{512\eta_{\rm leak}^{2}\mu_wq\rho_0}{\pi^4\lambda^2}\mathfrak L^2,
  - path: paper/stages/stage_247.tex
    line: 68
    role: paper_stage_tex
    excerpt: \Delta E_{UV}=\eta_{UV}\mathcal D_{UV},
  - path: paper/stages/stage_247.tex
    line: 72
    role: paper_stage_tex
    excerpt: \mathcal M_\sigma(r)=\xi_R\bigl[\mathcal R_\infty-\mathcal R(r)\bigr],
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
    line: 1
    role: notes_stage
    excerpt: '# Moving-Throat PDE — Stage 247: Relaxed Stationary Barrier Compiler from One-Port Short-Range, Leakage, `U/V`
      Drain, and Compensated Source Packets'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 17
    role: notes_stage
    excerpt: 2. the exact Stage-244 leakage / work compiler,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 33
    role: notes_stage
    excerpt: (S_{\rm leak},\,\mathcal W_w),
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 48
    role: notes_stage
    excerpt: 'It has four jobs:'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 51
    role: notes_stage
    excerpt: 2. import the Stage-244 leakage / work packet into barrier units,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 61
    role: notes_stage
    excerpt: '- **Stage 244**, which attached leakage and scalar-photon work to the selected-support packet,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 63
    role: notes_stage
    excerpt: '- **Stage 246**, which converted the sign-changing source branch into explicit mouth/source observables,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 98
    role: notes_stage
    excerpt: Q=G_U^2\Omega_W^2+2G_UG_WR_{\rm mix}+G_W^2\Omega_U^2,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 101
    role: notes_stage
    excerpt: P=\Omega_U^2G_W+R_{\rm mix}G_U,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 116
    role: notes_stage
    excerpt: \chi_{UU}=\frac{K_*\Omega_W^2-G_W^2}{\Delta D_0},
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 118
    role: notes_stage
    excerpt: \chi_{UW}=\frac{K_*R_{\rm mix}+G_UG_W}{\Delta D_0},
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 124
    role: notes_stage
    excerpt: P_U=G_U\Omega_W^2+R_{\rm mix}G_W.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 129
    role: notes_stage
    excerpt: \beta_Q,\qquad \beta_U,\qquad \beta_W,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 133
    role: notes_stage
    excerpt: \mathcal C_6=\chi_{qq}\beta_Q^2,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 136
    role: notes_stage
    excerpt: \mathcal C_4=\chi_{qU}\beta_Q\beta_U+\chi_{qW}\beta_Q\beta_W,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 139
    role: notes_stage
    excerpt: \mathcal C_2=\chi_{UU}\beta_U^2+2\chi_{UW}\beta_U\beta_W+\chi_{WW}\beta_W^2.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 187
    role: notes_stage
    excerpt: '### 2.1 Stage-244 leakage and reduced Session-I work scalar'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 196
    role: notes_stage
    excerpt: S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 198
    role: notes_stage
    excerpt: \frac{8\sqrt2\,\eta_{\rm leak}\mu_w\rho_0}{\pi^{5/2}\lambda^3}
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 206
    role: notes_stage
    excerpt: \frac{512\,\eta_{\rm leak}^2\mu_w q\rho_0}{\pi^4\lambda^2}
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 225
    role: notes_stage
    excerpt: \Delta E_{UV}(r):=\eta_{UV}\mathcal D_{UV}(r),
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 228
    role: notes_stage
    excerpt: where `\eta_{UV}` is the same Session-I barrier-lowering weight carried in the stationary script.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 247
    role: notes_stage
    excerpt: \mathcal M_\sigma(r):=\xi_R\bigl[\mathcal R_\infty-\mathcal R(r)\bigr],
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 296
    role: notes_stage
    excerpt: \lambda_L S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 308
    role: notes_stage
    excerpt: '- `\lambda_L\ge 0` converts the open-subsystem leakage source into stationary barrier units,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 320
    role: notes_stage
    excerpt: \lambda_L S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 344
    role: notes_stage
    excerpt: '- `(S_{\rm leak},\mathcal W_w^{\rm sess})` are support-side / open-system objects,'
  origin_claims:
  - parameter: mu_w
    introduced_at_stage: 247
    introduced_at_line: 447
    citation:
      path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
      line: 447
      excerpt: \mu_w=0.8,
  constraints:
  - parameter: mu_w
    constraint_kind: free_choice
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
      line: 441
      excerpt: with the Session-I parameters
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings: []
  downstream_dependents:
  - '248'
  synthesis_ingested_at: '2026-06-11T19:37:44Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_27/fit_stage247_session1_packet__mu_w.yaml
- schema_version: 1
  candidate_id: fit_stage247_session1_packet
  parameter_name: xi_R
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage247_session1_packet__phase_b__xi_r.md
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
    line: 4
    role: paper_stage_tex
    excerpt: '\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py}.  Mathematica
      audit: none yet.}'
  - path: paper/stages/stage_247.tex
    line: 12
    role: paper_stage_tex
    excerpt: Q=G_U^2\Omega_W^2+2G_UG_WR_{\rm mix}+G_W^2\Omega_U^2,
  - path: paper/stages/stage_247.tex
    line: 16
    role: paper_stage_tex
    excerpt: P=\Omega_U^2G_W+R_{\rm mix}G_U,
  - path: paper/stages/stage_247.tex
    line: 20
    role: paper_stage_tex
    excerpt: With primitive short-range source amplitudes \(\beta_Q,\beta_U,\beta_W\), the product-family coefficients are
  - path: paper/stages/stage_247.tex
    line: 23
    role: paper_stage_tex
    excerpt: \mathcal C_6=\chi_{qq}\beta_Q^2,
  - path: paper/stages/stage_247.tex
    line: 25
    role: paper_stage_tex
    excerpt: \mathcal C_4=\chi_{qU}\beta_Q\beta_U+\chi_{qW}\beta_Q\beta_W,
  - path: paper/stages/stage_247.tex
    line: 29
    role: paper_stage_tex
    excerpt: \mathcal C_2=\chi_{UU}\beta_U^2+2\chi_{UW}\beta_U\beta_W+\chi_{WW}\beta_W^2.
  - path: paper/stages/stage_247.tex
    line: 60
    role: paper_stage_tex
    excerpt: S_{\rm leak}=\frac{8\sqrt2\eta_{\rm leak}\mu_w\rho_0}{\pi^{5/2}\lambda^3}\mathfrak L,
  - path: paper/stages/stage_247.tex
    line: 66
    role: paper_stage_tex
    excerpt: \frac{512\eta_{\rm leak}^{2}\mu_wq\rho_0}{\pi^4\lambda^2}\mathfrak L^2,
  - path: paper/stages/stage_247.tex
    line: 68
    role: paper_stage_tex
    excerpt: \Delta E_{UV}=\eta_{UV}\mathcal D_{UV},
  - path: paper/stages/stage_247.tex
    line: 72
    role: paper_stage_tex
    excerpt: \mathcal M_\sigma(r)=\xi_R\bigl[\mathcal R_\infty-\mathcal R(r)\bigr],
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
    line: 1
    role: notes_stage
    excerpt: '# Moving-Throat PDE — Stage 247: Relaxed Stationary Barrier Compiler from One-Port Short-Range, Leakage, `U/V`
      Drain, and Compensated Source Packets'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 17
    role: notes_stage
    excerpt: 2. the exact Stage-244 leakage / work compiler,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 33
    role: notes_stage
    excerpt: (S_{\rm leak},\,\mathcal W_w),
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 48
    role: notes_stage
    excerpt: 'It has four jobs:'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 51
    role: notes_stage
    excerpt: 2. import the Stage-244 leakage / work packet into barrier units,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 61
    role: notes_stage
    excerpt: '- **Stage 244**, which attached leakage and scalar-photon work to the selected-support packet,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 63
    role: notes_stage
    excerpt: '- **Stage 246**, which converted the sign-changing source branch into explicit mouth/source observables,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 98
    role: notes_stage
    excerpt: Q=G_U^2\Omega_W^2+2G_UG_WR_{\rm mix}+G_W^2\Omega_U^2,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 101
    role: notes_stage
    excerpt: P=\Omega_U^2G_W+R_{\rm mix}G_U,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 116
    role: notes_stage
    excerpt: \chi_{UU}=\frac{K_*\Omega_W^2-G_W^2}{\Delta D_0},
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 118
    role: notes_stage
    excerpt: \chi_{UW}=\frac{K_*R_{\rm mix}+G_UG_W}{\Delta D_0},
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 124
    role: notes_stage
    excerpt: P_U=G_U\Omega_W^2+R_{\rm mix}G_W.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 129
    role: notes_stage
    excerpt: \beta_Q,\qquad \beta_U,\qquad \beta_W,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 133
    role: notes_stage
    excerpt: \mathcal C_6=\chi_{qq}\beta_Q^2,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 136
    role: notes_stage
    excerpt: \mathcal C_4=\chi_{qU}\beta_Q\beta_U+\chi_{qW}\beta_Q\beta_W,
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 139
    role: notes_stage
    excerpt: \mathcal C_2=\chi_{UU}\beta_U^2+2\chi_{UW}\beta_U\beta_W+\chi_{WW}\beta_W^2.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 187
    role: notes_stage
    excerpt: '### 2.1 Stage-244 leakage and reduced Session-I work scalar'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 196
    role: notes_stage
    excerpt: S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 198
    role: notes_stage
    excerpt: \frac{8\sqrt2\,\eta_{\rm leak}\mu_w\rho_0}{\pi^{5/2}\lambda^3}
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 206
    role: notes_stage
    excerpt: \frac{512\,\eta_{\rm leak}^2\mu_w q\rho_0}{\pi^4\lambda^2}
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 225
    role: notes_stage
    excerpt: \Delta E_{UV}(r):=\eta_{UV}\mathcal D_{UV}(r),
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 228
    role: notes_stage
    excerpt: where `\eta_{UV}` is the same Session-I barrier-lowering weight carried in the stationary script.
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 247
    role: notes_stage
    excerpt: \mathcal M_\sigma(r):=\xi_R\bigl[\mathcal R_\infty-\mathcal R(r)\bigr],
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 296
    role: notes_stage
    excerpt: \lambda_L S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 308
    role: notes_stage
    excerpt: '- `\lambda_L\ge 0` converts the open-subsystem leakage source into stationary barrier units,'
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 320
    role: notes_stage
    excerpt: \lambda_L S_{\rm leak}(r)
  - path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
    line: 344
    role: notes_stage
    excerpt: '- `(S_{\rm leak},\mathcal W_w^{\rm sess})` are support-side / open-system objects,'
  origin_claims:
  - parameter: xi_R
    introduced_at_stage: 247
    introduced_at_line: 476
    citation:
      path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
      line: 476
      excerpt: \xi_R=0.9,
  constraints:
  - parameter: xi_R
    constraint_kind: free_choice
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md
      line: 467
      excerpt: For the compensated source packet, using
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings: []
  downstream_dependents:
  - '248'
  synthesis_ingested_at: '2026-06-11T19:37:44Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_27/fit_stage247_session1_packet__xi_r.yaml
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
