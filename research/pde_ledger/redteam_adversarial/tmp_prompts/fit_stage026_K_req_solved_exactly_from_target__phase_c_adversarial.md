# Phase C Adversarial Audit Wrapper

Candidate: `fit_stage026_K_req_solved_exactly_from_target`

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
id: fit_stage026_K_req_solved_exactly_from_target
candidate_key: stage026_K_req_solved_exactly_from_target
dry_run: false
dry_run_id: null
anchor_stages:
- '026'
file_line_citations:
- path: notes/stages/moving_throat_pde_stage026_concrete_axial_overlaps.md
  line: 212
  excerpt: 'It can be solved exactly for the required wall stiffness:'
  stage: '026'
parameter_names:
- K_req
batch_id: null
status: provenance_built
status_timestamps:
  pending: '2026-06-11T08:21:51Z'
  scanned: '2026-06-11T08:21:51Z'
  provenance_built: '2026-06-11T18:25:36Z'
codex_session:
  by_parameter:
    K_req:
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
  - redteam_adversarial/provenance/fit_stage026_k_req_solved_exactly_from_target__k_req.yaml
  phase_c_prompt: null
verdict:
  adversarial: null
  adjudication: null
modality_attribution:
- claim_label
modality_fragments:
- candidate_key: stage026_K_req_solved_exactly_from_target
  modality: claim_label
  anchor_stage: '026'
  parameter_names:
  - K_req
  reason: The wall stiffness K_req is back-solved from the GR quadrupole target 54*G*c_s^5/(5*a^5*c^5) and labeled an exact
    solution; this is the clearest fit-shaped move in the band, a parameter defined by the value needed to match the target.
  citation:
    path: notes/stages/moving_throat_pde_stage026_concrete_axial_overlaps.md
    line: 212
    excerpt: 'It can be solved exactly for the required wall stiffness:'
    stage: '026'
phase_b_status: synthesis_complete
```

## Primary Sources

- notes/stages/moving_throat_pde_stage026_concrete_axial_overlaps.md

## Provenance Slice

```yaml
- schema_version: 1
  candidate_id: fit_stage026_K_req_solved_exactly_from_target
  parameter_name: K_req
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage026_k_req_solved_exactly_from_target__phase_b__k_req.md
  agent_synthesis_required: true
  dry_run_fixture: false
  fixture_note: null
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - '026'
  notes_source_opened: true
  source_evidence:
  - path: paper/stages/stage_026.tex
    line: 111
    role: paper_stage_tex
    excerpt: Solving for the required wall stiffness gives
  - path: paper/stages/stage_026.tex
    line: 113
    role: paper_stage_tex
    excerpt: \label{eq:app-stage026-Kreq}
  - path: paper/stages/stage_026.tex
    line: 115
    role: paper_stage_tex
    excerpt: K_{\rm req}
  - path: paper/stages/stage_026.tex
    line: 129
    role: paper_stage_tex
    excerpt: K_\eta+6T_\Omega=K_{\rm req}.
  - path: paper/stages/stage_026.tex
    line: 132
    role: paper_stage_tex
    excerpt: \stagefield{Output}{Stage~026 outputs the finite-throat D/N overlap \eqref{eq:app-stage026-kappa0}, the concrete
      couplings \eqref{eq:app-stage026-couplings}, the finite-throat normalization law \eqref{eq:app-stage026-normalization},
      and the required stiffness formula \eqref{eq:app-stage026-Kreq}.}
  - path: paper/stages/stage_026.tex
    line: 138
    role: paper_stage_tex
    excerpt: \item Solving \eqref{eq:app-stage026-normalization} for \(K\) gives \eqref{eq:app-stage026-Kreq}.
  - path: notes/stages/moving_throat_pde_stage026_concrete_axial_overlaps.md
    line: 28
    role: notes_stage
    excerpt: '- the brane-like internal gauge frequency,'
  - path: notes/stages/moving_throat_pde_stage026_concrete_axial_overlaps.md
    line: 29
    role: notes_stage
    excerpt: '- the mixed-channel frequency,'
  - path: notes/stages/moving_throat_pde_stage026_concrete_axial_overlaps.md
    line: 67
    role: notes_stage
    excerpt: 'For the trapped support channel and the mixed `A_w / F_(mu w) / J^w` channel, use the exact finite-throat D/N
      ladder already isolated in the frozen-wall benchmark:'
  - path: notes/stages/moving_throat_pde_stage026_concrete_axial_overlaps.md
    line: 85
    role: notes_stage
    excerpt: The associated support and mixed frequencies are taken to be
  - path: notes/stages/moving_throat_pde_stage026_concrete_axial_overlaps.md
    line: 108
    role: notes_stage
    excerpt: and keep `Omega_U` as the reduced internal restoring frequency of the brane-like zero mode `u_0` (so its frequency
      need not come from an axial gradient term).
  - path: notes/stages/moving_throat_pde_stage026_concrete_axial_overlaps.md
    line: 212
    role: notes_stage
    excerpt: 'It can be solved exactly for the required wall stiffness:'
  - path: notes/stages/moving_throat_pde_stage026_concrete_axial_overlaps.md
    line: 214
    role: notes_stage
    excerpt: '`K_req = kappa^2 lambda_B^2 / varpi^2 + Q / Delta`'
  - path: notes/stages/moving_throat_pde_stage026_concrete_axial_overlaps.md
    line: 238
    role: notes_stage
    excerpt: '`K_eta + 6 T_Omega = K_req`.'
  - path: notes/stages/moving_throat_pde_stage026_concrete_axial_overlaps.md
    line: 244
    role: notes_stage
    excerpt: It is to determine whether the actual quadrupolar wall stiffness of the throat satisfies the explicit algebraic
      relation above after support and internal gauge loading are included.
  - path: scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py
    line: 17
    role: sympy_script
    excerpt: • and exact solution of the normalization equation for the required wall
  - path: scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py
    line: 18
    role: sympy_script
    excerpt: stiffness K_req on the branch.
  - path: scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py
    line: 181
    role: sympy_script
    excerpt: '# IV. Exact branch-level normalization test and required wall stiffness'
  - path: scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py
    line: 195
    role: sympy_script
    excerpt: '# Solve exactly for the required wall stiffness and verify by back-substitution'
  - path: scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py
    line: 197
    role: sympy_script
    excerpt: K_req = sp.solve(sp.Eq(residual, 0), K)[0]
  - path: scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py
    line: 199
    role: sympy_script
    excerpt: subbanner("IV.2 — Exact required wall stiffness")
  - path: scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py
    line: 200
    role: sympy_script
    excerpt: print("K_req =")
  - path: scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py
    line: 201
    role: sympy_script
    excerpt: sp.pprint(sp.simplify(K_req))
  - path: scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py
    line: 203
    role: sympy_script
    excerpt: '# Independent structural check: the paper''s eq:app-stage026-Kreq states the'
  - path: scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py
    line: 204
    role: sympy_script
    excerpt: '# three-term decomposition K_req = B0 + Q/Delta + mhat^2 * kappa^2 * (...)^2'
  - path: scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py
    line: 206
    role: sympy_script
    excerpt: K_req_paper = (
  - path: scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py
    line: 212
    role: sympy_script
    excerpt: expect_zero("K_req - K_req_paper", sp.simplify(K_req - K_req_paper))
  - path: scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py
    line: 213
    role: sympy_script
    excerpt: expect_zero("residual @ K_req", residual.subs(K, K_req))
  origin_claims:
  - parameter: K_req
    introduced_at_stage: 26
    introduced_at_line: 214
    citation:
      path: notes/stages/moving_throat_pde_stage026_concrete_axial_overlaps.md
      line: 214
      excerpt: '`K_req = kappa^2 lambda_B^2 / varpi^2 + Q / Delta`'
  constraints:
  - parameter: K_req
    constraint_kind: published_target
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage026_concrete_axial_overlaps.md
      line: 218
      excerpt: 'So on this branch the GR quadrupole target is equivalent to one concrete statement:'
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings:
  - type: derive_vs_fit_mismatch
    severity: needs_triage
    summary: 'K_req is solved exactly so the one-mode prefactor equals the external GR quadrupole target 54/5 G c_s^5/(5 a^5
      c^5): the script computes K_req = solve(mhat^2 N0/D0 - target = 0, K), and the closed form K_req = B0 + Q/Delta + mhat^2
      kappa^2 (...)^2 / (target Delta^2) carries the external target in its denominator. So K_req''s value is defined by back-solving
      to the published GR number, not derived from ledger-internal stiffness physics. The algebra is exact, but the constraint
      that fixes the VALUE is the external target. The card frames this as ''solving for the required wall stiffness'' (line
      111), which is accurate, but the genealogy is published_target, not an independent derivation.'
    citations:
    - path: scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py
      line: 197
      excerpt: K_req = sp.solve(sp.Eq(residual, 0), K)[0]
    - path: scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py
      line: 210
      excerpt: / (target * Delta**2)
    - path: notes/stages/moving_throat_pde_stage026_concrete_axial_overlaps.md
      line: 7
      excerpt: '`mhat_rad^2 P^2 / [ Delta ( K Delta - Delta C^2 / varpi^2 - Q ) ] = 54 G c_s^5 / (5 a^5 c^5)`.'
  downstream_dependents:
  - '027'
  synthesis_ingested_at: '2026-06-11T18:25:36Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_03/fit_stage026_k_req_solved_exactly_from_target__k_req.yaml
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
