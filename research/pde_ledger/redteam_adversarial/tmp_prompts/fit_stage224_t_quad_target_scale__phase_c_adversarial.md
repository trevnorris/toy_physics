# Phase C Adversarial Audit Wrapper

Candidate: `fit_stage224_t_quad_target_scale`

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
id: fit_stage224_t_quad_target_scale
candidate_key: stage224_t_quad_target_scale
dry_run: false
dry_run_id: null
anchor_stages:
- '224'
file_line_citations:
- path: paper/stages/stage_224.tex
  line: 9
  excerpt: \stagefield{Inputs}{Imports the normalized prefactor packet $(\Delta_{\rm norm},a_{P_0},b_{P_0})$, the weak-axisymmetric
    signature $(1,1/2,-1)$, the line $b_{P_0}=3a_{P_0}$, and the target scale $T_{\rm quad}=54Gc_s^5/(5a^5c^5)$.}
  stage: '224'
parameter_names:
- T_quad
batch_id: null
status: provenance_built
status_timestamps:
  pending: '2026-06-11T08:21:51Z'
  scanned: '2026-06-11T08:21:51Z'
  provenance_built: '2026-06-11T19:17:46Z'
codex_session:
  by_parameter:
    T_quad:
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
  - redteam_adversarial/provenance/fit_stage224_t_quad_target_scale__t_quad.yaml
  phase_c_prompt: null
verdict:
  adversarial: null
  adjudication: null
modality_attribution:
- claim_label
modality_fragments:
- candidate_key: stage224_t_quad_target_scale
  modality: claim_label
  anchor_stage: '224'
  parameter_names:
  - T_quad
  reason: T_quad is presented as a fixed target scale with the same 54/5 coefficient family as P_0_target; the kill-test verdict
    is directly sensitive to this matched value, so its claimed fixedness is a high-leverage potential insertion point.
  citation:
    path: paper/stages/stage_224.tex
    line: 9
    excerpt: \stagefield{Inputs}{Imports the normalized prefactor packet $(\Delta_{\rm norm},a_{P_0},b_{P_0})$, the weak-axisymmetric
      signature $(1,1/2,-1)$, the line $b_{P_0}=3a_{P_0}$, and the target scale $T_{\rm quad}=54Gc_s^5/(5a^5c^5)$.}
    stage: '224'
phase_b_status: synthesis_complete
```

## Primary Sources

- paper/stages/stage_224.tex

## Provenance Slice

```yaml
- schema_version: 1
  candidate_id: fit_stage224_t_quad_target_scale
  parameter_name: T_quad
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage224_t_quad_target_scale__phase_b__t_quad.md
  agent_synthesis_required: true
  dry_run_fixture: false
  fixture_note: null
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - '224'
  notes_source_opened: true
  source_evidence:
  - path: paper/stages/stage_224.tex
    line: 9
    role: paper_stage_tex
    excerpt: \stagefield{Inputs}{Imports the normalized prefactor packet $(\Delta_{\rm norm},a_{P_0},b_{P_0})$, the weak-axisymmetric
      signature $(1,1/2,-1)$, the line $b_{P_0}=3a_{P_0}$, and the target scale $T_{\rm quad}=54Gc_s^5/(5a^5c^5)$.}
  - path: paper/stages/stage_224.tex
    line: 15
    role: paper_stage_tex
    excerpt: '\stagefield{Output}{Actual-branch kill test: the branch survives only if the transported quantity $(\Delta_{\rm
      norm}+T_{\rm quad})(1+|\varepsilon\Xi_1|)/\hat m_0^2$ stays below the critical prefactor budget.}'
  - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
    line: 34
    role: notes_stage
    excerpt: P_0^{(20)},\qquad P_0^{(21)},\qquad P_0^{(22)},
  - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
    line: 41
    role: notes_stage
    excerpt: \qquad
  - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
    line: 77
    role: notes_stage
    excerpt: \frac{\Delta_{\rm norm}+T_{\rm quad}}{\hat m_0^{\,2}},
  - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
    line: 78
    role: notes_stage
    excerpt: \qquad
  - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
    line: 79
    role: notes_stage
    excerpt: T_{\rm quad}:=\frac{54Gc_s^5}{5a^5c^5}.
  - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
    line: 104
    role: notes_stage
    excerpt: \qquad
  - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
    line: 111
    role: notes_stage
    excerpt: \qquad
  - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
    line: 124
    role: notes_stage
    excerpt: \qquad
  - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
    line: 126
    role: notes_stage
    excerpt: \qquad
  - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
    line: 139
    role: notes_stage
    excerpt: \frac{\Delta_{\rm norm}+T_{\rm quad}}{\hat m_0^{\,2}}+4a_{P_0},
  - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
    line: 140
    role: notes_stage
    excerpt: \frac{\Delta_{\rm norm}+T_{\rm quad}}{\hat m_0^{\,2}}-a_{P_0}+b_{P_0},
  - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
    line: 141
    role: notes_stage
    excerpt: \frac{\Delta_{\rm norm}+T_{\rm quad}}{\hat m_0^{\,2}}-a_{P_0}-b_{P_0}
  - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
    line: 170
    role: notes_stage
    excerpt: \hat m_0^{\,2}P_{\rm crit}-T_{\rm quad}.
  - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
    line: 177
    role: notes_stage
    excerpt: If the real branch already hits the universal quadrupole normalization exactly,
  - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
    line: 186
    role: notes_stage
    excerpt: \frac{T_{\rm quad}}{P_{\rm crit}}.
  - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
    line: 198
    role: notes_stage
    excerpt: \qquad
  - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
    line: 200
    role: notes_stage
    excerpt: \qquad
  - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
    line: 214
    role: notes_stage
    excerpt: \qquad
  - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
    line: 216
    role: notes_stage
    excerpt: \qquad
  - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
    line: 223
    role: notes_stage
    excerpt: \qquad
  - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
    line: 240
    role: notes_stage
    excerpt: \qquad
  - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
    line: 242
    role: notes_stage
    excerpt: \qquad
  - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
    line: 264
    role: notes_stage
    excerpt: \frac{\Delta_{\rm norm}+T_{\rm quad}}{\hat m_0^{\,2}}\bigl(1+|\epsilon\Xi_1|\bigr)
  - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
    line: 289
    role: notes_stage
    excerpt: \frac{T_{\rm quad}(1+|\epsilon\Xi_1|)}{P_{\rm crit}}.
  - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
    line: 306
    role: notes_stage
    excerpt: \qquad
  - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
    line: 317
    role: notes_stage
    excerpt: \qquad
  - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
    line: 325
    role: notes_stage
    excerpt: \qquad
  - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
    line: 336
    role: notes_stage
    excerpt: \qquad
  - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
    line: 344
    role: notes_stage
    excerpt: \qquad
  - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
    line: 404
    role: notes_stage
    excerpt: \bar P_0=\frac{\Delta_{\rm norm}+T_{\rm quad}}{\hat m_0^{\,2}},
  - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
    line: 408
    role: notes_stage
    excerpt: \Delta_{\rm norm}\le \hat m_0^{\,2}P_{\rm crit}-T_{\rm quad},
  - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
    line: 416
    role: notes_stage
    excerpt: \lambda_{20}=1,\quad \lambda_{21}=\frac12,\quad \lambda_{22}=-1,
  - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
    line: 421
    role: notes_stage
    excerpt: \qquad
  - path: scripts/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.py
    line: 21
    role: sympy_script
    excerpt: Delta_norm, T_quad, mhat0 = sp.symbols("Delta_norm T_quad mhat0", positive=True)
  - path: scripts/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.py
    line: 25
    role: sympy_script
    excerpt: Pbar = (Delta_norm + T_quad) / mhat0**2
  - path: scripts/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.py
    line: 50
    role: sympy_script
    excerpt: expected_rhs = sp.simplify(mhat0**2 * Pcrit - T_quad)
  - path: scripts/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.py
    line: 54
    role: sympy_script
    excerpt: calibrated_lower_bound = sp.simplify(T_quad / Pcrit)
  - path: scripts/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.py
    line: 122
    role: sympy_script
    excerpt: calibrated_weak_bound = sp.simplify(T_quad * (1 + sp.Abs(eps * Xi1)) / Pcrit)
  origin_claims:
  - parameter: T_quad
    introduced_at_stage: 22
    introduced_at_line: 270
    citation:
      path: notes/stages/moving_throat_pde_stage022_grouped_p2_normalization_bridge.md
      line: 270
      excerpt: '`mhat_0^2 * P_0 = 54 G c_s^5 / (5 a^5 c^5)`.'
  constraints:
  - parameter: T_quad
    constraint_kind: published_target
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
      line: 79
      excerpt: T_{\rm quad}:=\frac{54Gc_s^5}{5a^5c^5}.
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings:
  - type: derive_vs_fit_mismatch
    severity: low
    summary: T_quad=54Gc_s^5/(5a^5c^5) is exactly the external universal GR quadrupole normalization scale (back-solved from
      gamma_GR=2G/5c^5 at stage 022), used at stage 224 as the calibration target inside Delta_norm=mhat_0^2 P0bar - T_quad.
      It is a published_target (fit to an external GR number), not an internal-consistency derivation; the stage 224 card
      and notes label it the 'target scale', consistent with that.
    citations:
    - path: notes/stages/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.md
      line: 71
      excerpt: \frac{54Gc_s^5}{5a^5c^5}.
    - path: notes/stages/moving_throat_pde_stage022_grouped_p2_normalization_bridge.md
      line: 256
      excerpt: The universal GR quadrupole target is
  downstream_dependents:
  - '225'
  synthesis_ingested_at: '2026-06-11T19:17:46Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_24/fit_stage224_t_quad_target_scale__t_quad.yaml
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
