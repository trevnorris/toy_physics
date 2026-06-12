# Phase C Adversarial Audit Wrapper

Candidate: `fit_stage248_v0_cross_benchmark_speed`

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
id: fit_stage248_v0_cross_benchmark_speed
candidate_key: stage248_v0_cross_benchmark_speed
dry_run: false
dry_run_id: null
anchor_stages:
- '248'
file_line_citations:
- path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
  line: 513
  excerpt: The Session-II report also gave an explicit above-threshold demonstration speed
  stage: '248'
parameter_names:
- v0_cross
batch_id: null
status: provenance_built
status_timestamps:
  pending: '2026-06-11T08:21:51Z'
  scanned: '2026-06-11T08:21:51Z'
  provenance_built: '2026-06-11T19:37:44Z'
codex_session:
  by_parameter:
    v0_cross:
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
  - redteam_adversarial/provenance/fit_stage248_v0_cross_benchmark_speed__v0_cross.yaml
  phase_c_prompt: null
verdict:
  adversarial: null
  adjudication: null
modality_attribution:
- claim_label
modality_fragments:
- candidate_key: stage248_v0_cross_benchmark_speed
  modality: claim_label
  anchor_stage: '248'
  parameter_names:
  - v0_cross
  reason: v0_cross = 2.59221845 is imported from a session report and then shown to lie "inside the exact window"; the exact-window
    framing certifies an inherited demonstration value whose original selection criteria are not re-derived.
  citation:
    path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 513
    excerpt: The Session-II report also gave an explicit above-threshold demonstration speed
    stage: '248'
phase_b_status: synthesis_complete
```

## Primary Sources

- notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md

## Provenance Slice

```yaml
- schema_version: 1
  candidate_id: fit_stage248_v0_cross_benchmark_speed
  parameter_name: v0_cross
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage248_v0_cross_benchmark_speed__phase_b__v0_cross.md
  agent_synthesis_required: true
  dry_run_fixture: false
  fixture_note: null
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - '248'
  notes_source_opened: true
  source_evidence:
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 16
    role: notes_stage
    excerpt: 4. and the declared Session-II benchmark specialization.
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 55
    role: notes_stage
    excerpt: and the declared Session-II benchmark specialization.
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 58
    role: notes_stage
    excerpt: Session-II benchmark numbers explicitly confined to the benchmark-only
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 68
    role: notes_stage
    excerpt: '- The same-charge barrier audit had already shown that the linear dynamic mixed bundle does **not** generate
      a new conservative kernel class at large distance. Linear time dependence only makes the already-admissible short-range
      families frequency dependent, while the first outgoing correction is phase lag / pumping rather than direct conservative
      barrier lowering.'
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 144
    role: notes_stage
    excerpt: The minimum launch speed required to reach the top of the reduced barrier is determined by
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 180
    role: notes_stage
    excerpt: Let the lowered branch be single-peaked on \([r_{\rm contact},r_0]\), so that crossing the peak is sufficient
      to reach contact. Then if
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 323
    role: notes_stage
    excerpt: So once the lowered branch gives \(I_{\rm new}<I_{\rm Coul}\), the enhancement is completely fixed.
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 362
    role: notes_stage
    excerpt: This normal form is not needed for the Session-II benchmark itself, but it is the right local compiler for later
      near-threshold analysis.
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 396
    role: notes_stage
    excerpt: So the same turning-point event chain also determines the first dynamic confinement-width trigger.
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 405
    role: notes_stage
    excerpt: '## 6. Session-II benchmark specialization'
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 515
    role: notes_stage
    excerpt: v_{0,\rm cross}=2.59221845.
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 525
    role: notes_stage
    excerpt: E_{\rm Coul}(v_{0,\rm cross})
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 527
    role: notes_stage
    excerpt: \frac12 v_{0,\rm cross}^2+\frac15,
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 531
    role: notes_stage
    excerpt: r_{\rm turn,Coul}(v_{0,\rm cross})
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 533
    role: notes_stage
    excerpt: \frac{1}{E_{\rm Coul}(v_{0,\rm cross})}
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 538
    role: notes_stage
    excerpt: 'So the benchmark does exactly what Stage 248 needs it to do: it exhibits a concrete speed range where the lowered
      branch reaches the chosen contact scale and the Coulomb comparison does not.'
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 555
    role: notes_stage
    excerpt: That is within the numerical tolerance implied by the reported Session-II discrete reference value \(0.30222297\),
      so the closed-form Coulomb side of the compiler is consistent with the dynamic benchmark.
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 566
    role: notes_stage
    excerpt: 2. the barrier peak determines a finite launch threshold
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 586
    role: notes_stage
    excerpt: The event-chain improvement still sits on a short-range/open-system front end, and the linear dynamic mixed bundle
      still contributes phase-lag / pumping rather than a new conservative kernel class.
  - path: scripts/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.py
    line: 249
    role: sympy_script
    excerpt: vcross_num = 2.59221845
  - path: scripts/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.py
    line: 274
    role: sympy_script
    excerpt: Ecoul_cross_num = 0.5 * vcross_num**2 + 1.0 / r0_num
  - path: scripts/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.py
    line: 275
    role: sympy_script
    excerpt: rcoul_turn_cross_num = 1.0 / Ecoul_cross_num
  - path: scripts/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.py
    line: 297
    role: sympy_script
    excerpt: print("v_cross(session)      =", vcross_num)
  - path: scripts/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.py
    line: 298
    role: sympy_script
    excerpt: print("Coulomb turnback at v_cross =", rcoul_turn_cross_num)
  - path: scripts/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.py
    line: 308
    role: sympy_script
    excerpt: assert vcrit_num < vcross_num < vcontact_num
  - path: scripts/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.py
    line: 309
    role: sympy_script
    excerpt: assert abs(rcoul_turn_cross_num - rcoul_turn_report) < 3e-6
  origin_claims:
  - parameter: v0_cross
    introduced_at_stage: 248
    introduced_at_line: 515
    citation:
      path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
      line: 515
      excerpt: v_{0,\rm cross}=2.59221845.
  constraints:
  - parameter: v0_cross
    constraint_kind: published_target
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
      line: 513
      excerpt: The Session-II report also gave an explicit above-threshold demonstration speed
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings:
  - type: published_target
    severity: needs_triage
    summary: v0_cross=2.59221845 is a recorded figure imported from the Session-II report ('The Session-II report also gave
      an explicit above-threshold demonstration speed v0_cross=2.59221845', notes line 513-515; hardcoded in script line 249
      vcross_num = 2.59221845). Stage 248 uses it only to demonstrate it falls inside the independently-computed classical
      window 2.54139063 < v0 < 3.27278339 (notes line 517-519) and to reproduce the reported Coulomb turning point 0.28091705.
      It is a recorded benchmark target the compiler is checked against, not a stage-248 derivation.
    citations:
    - path: scripts/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.py
      line: 249
      excerpt: vcross_num = 2.59221845
    - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
      line: 513
      excerpt: The Session-II report also gave an explicit above-threshold demonstration speed
  downstream_dependents: []
  synthesis_ingested_at: '2026-06-11T19:37:44Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_27/fit_stage248_v0_cross_benchmark_speed__v0_cross.yaml
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
