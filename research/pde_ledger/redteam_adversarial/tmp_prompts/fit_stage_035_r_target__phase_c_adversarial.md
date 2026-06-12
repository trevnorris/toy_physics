# Phase C Adversarial Audit Wrapper

Candidate: `fit_stage_035_r_target`

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
id: fit_stage_035_r_target
candidate_key: stage_035_r_target
dry_run: false
dry_run_id: null
anchor_stages:
- '035'
file_line_citations:
- path: notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
  line: 211
  role: notes_stage
  stage: '035'
  excerpt: 1. the unique stable normalization locus `xi_req` is determined by `F(xi,delta)=R_target`,
parameter_names:
- R_target
- xi_req
batch_id: null
status: provenance_built
status_timestamps:
  pending: '2026-06-11T08:11:49Z'
  scanned: '2026-06-11T08:11:49Z'
  provenance_built: '2026-06-11T18:25:49Z'
codex_session:
  by_parameter:
    R_target:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    xi_req:
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
  - redteam_adversarial/provenance/fit_stage_035_r_target__r_target.yaml
  phase_c_prompt: null
verdict:
  adversarial: null
  adjudication: null
modality_attribution:
- claim_label
- numeric_literal
modality_fragments:
- candidate_key: stage_035_r_target
  modality: numeric_literal
  anchor_stage: '035'
  parameter_names:
  - R_target
  - xi_req
  reason: target-related numerical literal or closed-form coefficient
  citation:
    path: notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
    line: 211
    role: notes_stage
    stage: '035'
    excerpt: 1. the unique stable normalization locus `xi_req` is determined by `F(xi,delta)=R_target`,
- candidate_key: stage_035_r_target
  modality: claim_label
  anchor_stage: '035'
  parameter_names:
  - R_target
  - xi_req
  reason: claim label or status wording near target-related parameter
  citation:
    path: notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
    line: 211
    role: notes_stage
    stage: '035'
    excerpt: 1. the unique stable normalization locus `xi_req` is determined by `F(xi,delta)=R_target`,
phase_b_status: synthesis_complete
```

## Primary Sources

- notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md

## Provenance Slice

```yaml
- schema_version: 1
  candidate_id: fit_stage_035_r_target
  parameter_name: R_target
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage_035_r_target__phase_b__r_target.md
  agent_synthesis_required: true
  dry_run_fixture: false
  fixture_note: null
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - '035'
  notes_source_opened: true
  source_evidence:
  - path: paper/stages/stage_035.tex
    line: 7
    role: paper_stage_tex
    excerpt: \stagefield{Purpose}{Stage~035 non-dimensionalizes the softening-depth target.  The selected normalization product
      factors into an onset value times a universal D/N shape function \(F(\xi,\delta)\).  Strict monotonicity of \(F\) upgrades
      the Stage~033 onset inequality into a unique stable normalization locus.}
  - path: paper/stages/stage_035.tex
    line: 50
    role: paper_stage_tex
    excerpt: The dimensionless target ratio is
  - path: paper/stages/stage_035.tex
    line: 52
    role: paper_stage_tex
    excerpt: \label{eq:app-stage035-Rtarget}
  - path: paper/stages/stage_035.tex
    line: 54
    role: paper_stage_tex
    excerpt: R_{\rm target}:=\frac{N_Q^{\rm target}}{N_-(0)}
  - path: paper/stages/stage_035.tex
    line: 55
    role: paper_stage_tex
    excerpt: =\frac{N_Q^{\rm target}A}{\beta_0\kappa_0^2}}.
  - path: paper/stages/stage_035.tex
    line: 60
    role: paper_stage_tex
    excerpt: \boxed{F(\xi,\delta)=R_{\rm target}}.
  - path: paper/stages/stage_035.tex
    line: 92
    role: paper_stage_tex
    excerpt: R_{\rm target}<1\Rightarrow\text{no stable hit},
  - path: paper/stages/stage_035.tex
    line: 94
    role: paper_stage_tex
    excerpt: R_{\rm target}=1\Rightarrow\xi_{\rm req}=0,
  - path: paper/stages/stage_035.tex
    line: 96
    role: paper_stage_tex
    excerpt: R_{\rm target}>1\Rightarrow\exists!\,\xi_{\rm req}\in(0,1)}.
  - path: paper/stages/stage_035.tex
    line: 120
    role: paper_stage_tex
    excerpt: If \(R_{\rm target}=1+\varepsilon_R\), then
  - path: paper/stages/stage_035.tex
    line: 125
    role: paper_stage_tex
    excerpt: For large target ratio,
  - path: paper/stages/stage_035.tex
    line: 128
    role: paper_stage_tex
    excerpt: 1-\xi_{\rm req}\simeq\frac{C(\delta)}{R_{\rm target}}.
  - path: paper/stages/stage_035.tex
    line: 131
    role: paper_stage_tex
    excerpt: \stagefield{Output}{Stage~035 outputs the shape function \eqref{eq:app-stage035-F}, the target ratio \eqref{eq:app-stage035-Rtarget},
      the monotonicity derivative \eqref{eq:app-stage035-F-derivative}, the existence/uniqueness theorem \eqref{eq:app-stage035-existence},
      and the required total loading \eqref{eq:app-stage035-alpha-req}.}
  - path: paper/stages/stage_035.tex
    line: 140
    role: paper_stage_tex
    excerpt: \stagefield{Downstream use}{Stage~036 pairs \(F\) with the support-feasibility function \(G\).  Part~III later
      computes continuum placement data that determine \(R_{\rm target}\), \(M_{\rm mix}\), and \(\delta\).}
  - path: notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
    line: 13
    role: notes_stage
    excerpt: '`F(xi,delta) = R_target`,'
  - path: notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
    line: 65
    role: notes_stage
    excerpt: So the full selected-branch target becomes
  - path: notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
    line: 67
    role: notes_stage
    excerpt: '`F(xi,delta) = R_target`,'
  - path: notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
    line: 71
    role: notes_stage
    excerpt: '`R_target := N_Q^(target) / N_-(0)`'
  - path: notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
    line: 72
    role: notes_stage
    excerpt: '`          = N_Q^(target) A / ( beta_0 kappa_0^2 )`.'
  - path: notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
    line: 75
    role: notes_stage
    excerpt: All dependence on the microscopic transfer factor and the overall wall scale has collapsed into the single dimensionless
      target ratio `R_target`.
  - path: notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
    line: 108
    role: notes_stage
    excerpt: '- if `R_target < 1`, the target is impossible on the stable branch,'
  - path: notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
    line: 109
    role: notes_stage
    excerpt: '- if `R_target = 1`, the only hit is the unloaded onset point `xi = 0`,'
  - path: notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
    line: 110
    role: notes_stage
    excerpt: '- if `R_target > 1`, there is one and only one stable selected-branch softening depth `xi_req` that hits the
      target.'
  - path: notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
    line: 129
    role: notes_stage
    excerpt: So once the unique `xi_req` solving `F(xi,delta)=R_target` is known, the unique required total directional loading
      follows immediately from the formula above.
  - path: notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
    line: 171
    role: notes_stage
    excerpt: So if the target is only slightly above onset,
  - path: notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
    line: 173
    role: notes_stage
    excerpt: '`R_target = 1 + eps_R`,'
  - path: notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
    line: 185
    role: notes_stage
    excerpt: For very large target ratio `R_target`, the unique solution lies close to softening and obeys
  - path: notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
    line: 187
    role: notes_stage
    excerpt: '`1 - xi_req ~= C(delta) / R_target`,'
  - path: notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
    line: 201
    role: notes_stage
    excerpt: The completed moving-throat PDE has to determine the microscopic quantities entering
  - path: notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
    line: 203
    role: notes_stage
    excerpt: '`R_target = N_Q^(target) A / ( beta_0 kappa_0^2 )`,'
  - path: notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
    line: 211
    role: notes_stage
    excerpt: 1. the unique stable normalization locus `xi_req` is determined by `F(xi,delta)=R_target`,
  - path: scripts/moving_throat_pde_stage035_dimensionless_normalization_locus_sympy_audit.py
    line: 54
    role: sympy_script
    excerpt: F_target = sp.simplify((9 * delta + 11 * xi) ** 4 / (81 * (1 - xi) * (9 * delta**2 + 18 * delta * xi + 11 * xi**2)
      ** 2))
  - path: scripts/moving_throat_pde_stage035_dimensionless_normalization_locus_sympy_audit.py
    line: 58
    role: sympy_script
    excerpt: expect_zero("F - closed D/N form", F - F_target)
  - path: scripts/moving_throat_pde_stage035_dimensionless_normalization_locus_sympy_audit.py
    line: 60
    role: sympy_script
    excerpt: R_target = sp.simplify(NQ * A / (beta0 * kappa0_sq))
  - path: scripts/moving_throat_pde_stage035_dimensionless_normalization_locus_sympy_audit.py
    line: 61
    role: sympy_script
    excerpt: print("R_target =", R_target)
  - path: scripts/moving_throat_pde_stage035_dimensionless_normalization_locus_sympy_audit.py
    line: 62
    role: sympy_script
    excerpt: 'print("Normalization locus: F(xi,delta) = R_target")'
  - path: scripts/moving_throat_pde_stage035_dimensionless_normalization_locus_sympy_audit.py
    line: 65
    role: sympy_script
    excerpt: dF = sp.simplify(sp.diff(F_target, xi))
  - path: scripts/moving_throat_pde_stage035_dimensionless_normalization_locus_sympy_audit.py
    line: 66
    role: sympy_script
    excerpt: dF_target = sp.simplify(
  - path: scripts/moving_throat_pde_stage035_dimensionless_normalization_locus_sympy_audit.py
    line: 73
    role: sympy_script
    excerpt: expect_zero("dF/dxi - manifestly positive form", dF - dF_target)
  - path: scripts/moving_throat_pde_stage035_dimensionless_normalization_locus_sympy_audit.py
    line: 74
    role: sympy_script
    excerpt: expect_zero("F(0,delta) - 1", sp.simplify(F_target.subs(xi, 0) - 1))
  origin_claims:
  - parameter: R_target
    introduced_at_stage: 29
    introduced_at_line: 71
    citation:
      path: notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
      line: 71
      excerpt: '`R_target := N_Q^(target) / N_-(0)`'
  constraints:
  - parameter: R_target
    constraint_kind: published_target
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
      line: 72
      excerpt: '`          = N_Q^(target) A / ( beta_0 kappa_0^2 )`.'
  graph_context:
    contexts: []
    graph_gaps:
    - source: notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
      attempted_sources:
      - notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
      - research/pde_ledger/notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
      - paper/stages/stage_035.tex
      - research/pde_ledger/paper/stages/stage_035.tex
      graph_gap: true
      reason: no atlas node tied to this source path
  provenance_findings:
  - type: graph_gap
    severity: needs_triage
    summary: The atlas graph wrapper returned no source nodes for one or more primary sources.
    gaps:
    - source: notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
      attempted_sources:
      - notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
      - research/pde_ledger/notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
      - paper/stages/stage_035.tex
      - research/pde_ledger/paper/stages/stage_035.tex
      graph_gap: true
      reason: no atlas node tied to this source path
  - type: constraint_kind_ambiguous
    severity: needs_triage
    summary: 'R_target := N_Q^target / N_-(0) is a ratio of a ledger-internal onset value N_-(0) (= beta_0 kappa_0^2/A, internal_consistency)
      to the EXTERNAL 2.5PN GR target N_Q^target = 54 G c_s^5/(5 a^5 c^5). At stage 035 R_target is carried symbolically (notes
      line 75: ''All dependence ... has collapsed into the single dimensionless target ratio R_target''; line 201-203: ''The
      completed moving-throat PDE has to determine the microscopic quantities entering R_target''), so it is not yet a numeric
      back-solve. Classified published_target because its DEFINING role is to package the external radiative demand against
      an internal scale; the competing reading is internal_consistency for the ratio''s algebra. The external value is committed
      downstream (Stage 038 Lambda, script line 90/140).'
    citations:
    - path: notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
      line: 75
      excerpt: All dependence on the microscopic transfer factor and the overall wall scale has collapsed into the single
        dimensionless target ratio `R_target`.
    - path: scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py
      line: 65
      excerpt: (54 * G * c_s**5 / (5 * a**5 * c_light**5)) * A / ((8 / sp.pi**2) * beta0)
  downstream_dependents:
  - '030'
  - 038
  - 039
  synthesis_ingested_at: '2026-06-11T18:25:49Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_04/fit_stage_035_r_target__r_target.yaml
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
- source: notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
  attempted_sources:
  - notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
  - research/pde_ledger/notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
  - paper/stages/stage_035.tex
  - research/pde_ledger/paper/stages/stage_035.tex
  graph_gap: true
  reason: no atlas node tied to this source path
```

Output a markdown adversarial report only when this is not a dry run. For dry runs, stop after confirming the prompt has enough context and record no verdict.
