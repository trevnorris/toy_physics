# Phase C Adversarial Audit Wrapper

Candidate: `fit_stage106_odd_exactly_gr_target`

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
id: fit_stage106_odd_exactly_gr_target
candidate_key: stage106_odd_exactly_gr_target
dry_run: false
dry_run_id: null
anchor_stages:
- '106'
file_line_citations:
- path: notes/stages/moving_throat_pde_stage106_canonical_outgoing_reduced_closure.md
  line: 53
  excerpt: 'Equivalently, the normalized odd coefficient is exactly the GR/Burke–Thorne target:'
  stage: '106'
parameter_names:
- Gamma_5_normalized
batch_id: null
status: provenance_built
status_timestamps:
  pending: '2026-06-11T08:21:51Z'
  scanned: '2026-06-11T08:21:51Z'
  provenance_built: '2026-06-11T18:54:50Z'
codex_session:
  by_parameter:
    Gamma_5_normalized:
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
  - redteam_adversarial/provenance/fit_stage106_odd_exactly_gr_target__gamma_5_normalized.yaml
  phase_c_prompt: null
verdict:
  adversarial: null
  adjudication: null
modality_attribution:
- claim_label
modality_fragments:
- candidate_key: stage106_odd_exactly_gr_target
  modality: claim_label
  anchor_stage: '106'
  parameter_names:
  - Gamma_5_normalized
  reason: Exact-agreement-with-GR claim for the normalized odd coefficient; agreement achieved via the canonical-branch assumption
    could be back-filled from the target.
  citation:
    path: notes/stages/moving_throat_pde_stage106_canonical_outgoing_reduced_closure.md
    line: 53
    excerpt: 'Equivalently, the normalized odd coefficient is exactly the GR/Burke–Thorne target:'
    stage: '106'
phase_b_status: synthesis_complete
```

## Primary Sources

- notes/stages/moving_throat_pde_stage106_canonical_outgoing_reduced_closure.md

## Provenance Slice

```yaml
- schema_version: 1
  candidate_id: fit_stage106_odd_exactly_gr_target
  parameter_name: Gamma_5_normalized
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage106_odd_exactly_gr_target__phase_b__gamma_5_normalized.md
  agent_synthesis_required: true
  dry_run_fixture: false
  fixture_note: null
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - '106'
  notes_source_opened: true
  source_evidence:
  - path: paper/stages/stage_106.tex
    line: 24
    role: paper_stage_tex
    excerpt: \item Check the outgoing \(l=2\) DtN fingerprint against the normalized \(z=\omega a/c_s\) expansion.
  - path: notes/stages/moving_throat_pde_stage106_canonical_outgoing_reduced_closure.md
    line: 40
    role: notes_stage
    excerpt: 'With \(N_Q=1\), the canonical invariant low-frequency quadrupole coefficients are fixed to their target values:'
  - path: notes/stages/moving_throat_pde_stage106_canonical_outgoing_reduced_closure.md
    line: 49
    role: notes_stage
    excerpt: \overline \Gamma_5=\frac{2G}{5c^5}.
  - path: notes/stages/moving_throat_pde_stage106_canonical_outgoing_reduced_closure.md
    line: 53
    role: notes_stage
    excerpt: 'Equivalently, the normalized odd coefficient is exactly the GR/Burke–Thorne target:'
  - path: notes/stages/moving_throat_pde_stage106_canonical_outgoing_reduced_closure.md
    line: 56
    role: notes_stage
    excerpt: \gamma_{\rm quad}^{\rm eff}
  - path: notes/stages/moving_throat_pde_stage106_canonical_outgoing_reduced_closure.md
    line: 58
    role: notes_stage
    excerpt: \hat m_0^{\,2}\Gamma_5
  - path: notes/stages/moving_throat_pde_stage106_canonical_outgoing_reduced_closure.md
    line: 75
    role: notes_stage
    excerpt: Therefore the reduced nonspinning point-particle 2.5PN theorem is **closed on the canonical outgoing DtN branch**,
      conditional on that branch realization and the strict point-particle limit.
  - path: notes/stages/moving_throat_pde_stage106_canonical_outgoing_reduced_closure.md
    line: 79
    role: notes_stage
    excerpt: 'What remains open is no longer a reduced PN bookkeeping problem. It is the deeper PDE-side branch-selection
      theorem:'
  - path: scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py
    line: 6
    role: sympy_script
    excerpt: coefficient" and (iii) "outgoing l=2 DtN fingerprint against the normalized
  - path: scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py
    line: 53
    role: sympy_script
    excerpt: Gamma5_target = sp.simplify(2 * G / (5 * c**5))
  - path: scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py
    line: 58
    role: sympy_script
    excerpt: Gamma5 = sp.simplify(NQ_general * Gamma5_target)
  - path: scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py
    line: 63
    role: sympy_script
    excerpt: print("Gamma5 =", Gamma5)
  - path: scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py
    line: 70
    role: sympy_script
    excerpt: '"target identity Gamma5_target - 9 sqrt(K2_target^5 / K0_target^3)",'
  - path: scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py
    line: 71
    role: sympy_script
    excerpt: Gamma5_target - 9 * sp.sqrt(K2_target**5 / K0_target**3),
  - path: scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py
    line: 74
    role: sympy_script
    excerpt: gamma_eff_canonical = sp.simplify((m0hat**2 * Gamma5).subs(chi_Q, 1))
  - path: scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py
    line: 76
    role: sympy_script
    excerpt: '"canonical gamma_eff - target",'
  - path: scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py
    line: 77
    role: sympy_script
    excerpt: gamma_eff_canonical - Gamma5_target,
  - path: scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py
    line: 81
    role: sympy_script
    excerpt: '# N_Q = 1/chi_Q, so gamma_eff = m0hat^2 * N_Q * Gamma5_target collapses to'
  - path: scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py
    line: 82
    role: sympy_script
    excerpt: '# Gamma5_target / chi_Q. Expanding around chi_Q = 1 + Delta_Q gives a real'
  - path: scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py
    line: 83
    role: sympy_script
    excerpt: '# first-order slope of -Gamma5_target. This exercises the algebra non-trivially'
  - path: scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py
    line: 86
    role: sympy_script
    excerpt: gamma_eff_off = (m0hat**2 * Gamma5).subs([(m0hat, 1), (chi_Q, 1 + Delta_Q)])
  - path: scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py
    line: 87
    role: sympy_script
    excerpt: gamma_eff_series = sp.series(gamma_eff_off, Delta_Q, 0, 2).removeO()
  - path: scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py
    line: 88
    role: sympy_script
    excerpt: linear_coeff = sp.simplify(gamma_eff_series.coeff(Delta_Q, 1))
  - path: scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py
    line: 93
    role: sympy_script
    excerpt: zeroth_coeff = sp.simplify(gamma_eff_series.coeff(Delta_Q, 0))
  - path: scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py
    line: 95
    role: sympy_script
    excerpt: '"Delta_Q zeroth-order coefficient equals Gamma5_target",'
  - path: scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py
    line: 96
    role: sympy_script
    excerpt: zeroth_coeff - Gamma5_target,
  - path: scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py
    line: 102
    role: sympy_script
    excerpt: print("  The effective odd coefficient m0hat^2 Gamma5 then reproduces 2G/(5 c^5) exactly.")
  origin_claims:
  - parameter: Gamma_5_normalized
    introduced_at_stage: 106
    introduced_at_line: 56
    citation:
      path: notes/stages/moving_throat_pde_stage106_canonical_outgoing_reduced_closure.md
      line: 53
      excerpt: 'Equivalently, the normalized odd coefficient is exactly the GR/Burke–Thorne target:'
  constraints:
  - parameter: Gamma_5_normalized
    constraint_kind: published_target
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage106_canonical_outgoing_reduced_closure.md
      line: 53
      excerpt: 'Equivalently, the normalized odd coefficient is exactly the GR/Burke–Thorne target:'
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings:
  - type: published_target
    severity: low
    summary: The normalized odd coefficient gamma_quad^eff = m_hat0^2 Gamma_5 = 2G/5c^5 is set equal to the external GR/Burke-Thorne
      target. The chi_Q=1 -> N_Q=1 closure (internal, from 104/105) is what makes the EFFECTIVE coefficient land on this published
      number; the target value 2G/5c^5 itself is imported, not derived. The notes label it 'exactly the GR/Burke-Thorne target'
      and the closure is conditional (StatusOpen on actual branch realization), so the match is disclosed honestly -> clean
      published_target, no derive_vs_fit_mismatch.
    citations:
    - path: notes/stages/moving_throat_pde_stage106_canonical_outgoing_reduced_closure.md
      line: 60
      excerpt: \frac{2G}{5c^5}
    - path: scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py
      line: 102
      excerpt: print("  The effective odd coefficient m0hat^2 Gamma5 then reproduces 2G/(5 c^5) exactly.")
  downstream_dependents:
  - '107'
  synthesis_ingested_at: '2026-06-11T18:54:50Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_11/fit_stage106_odd_exactly_gr_target__gamma_5_normalized.yaml
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
