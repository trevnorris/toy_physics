# Phase C Adversarial Audit Wrapper

Candidate: `fit_stage248_xi_turn_lambda_th_carried_hardcodes`

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
id: fit_stage248_xi_turn_lambda_th_carried_hardcodes
candidate_key: stage248_xi_turn_lambda_th_carried_hardcodes
dry_run: false
dry_run_id: null
anchor_stages:
- '248'
file_line_citations:
- path: research/pde_ledger/scripts/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.py
  line: 247
  excerpt: Xi_turn_num = 0.34437471
  stage: '248'
parameter_names:
- Xi_turn
- lambda_th
batch_id: null
status: provenance_built
status_timestamps:
  pending: '2026-06-11T08:39:39Z'
  scanned: '2026-06-11T08:39:39Z'
  provenance_built: '2026-06-11T19:37:44Z'
codex_session:
  by_parameter:
    Xi_turn:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    lambda_th:
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
  - redteam_adversarial/provenance/fit_stage248_xi_turn_lambda_th_carried_hardcodes__xi_turn.yaml
  phase_c_prompt: null
verdict:
  adversarial: null
  adjudication: null
modality_attribution:
- completeness_critic
modality_fragments:
- candidate_key: stage248_xi_turn_lambda_th_carried_hardcodes
  modality: completeness_critic
  anchor_stage: '248'
  parameter_names:
  - Xi_turn
  - lambda_th
  reason: Xi_turn=0.34437471 and lambda_th=0.42826825 (line 248) are bare benchmark hardcodes re-typed independently in the
    stage249 script (line 172) and re-quoted as lambda_eff in the stage250/253 notes (e.g. stage250 line 412), yet neither
    appears in any fragment — the b8 scans flagged the neighboring Session-IV literals (t_cross, tcollapse0, gamma_crit) but
    missed this carried pair on the same event chain.
  citation:
    path: research/pde_ledger/scripts/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.py
    line: 247
    excerpt: Xi_turn_num = 0.34437471
    stage: '248'
phase_b_status: synthesis_complete
```

## Primary Sources

- research/pde_ledger/scripts/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.py

## Provenance Slice

```yaml
- schema_version: 1
  candidate_id: fit_stage248_xi_turn_lambda_th_carried_hardcodes
  parameter_name: Xi_turn
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage248_xi_turn_lambda_th_carried_hardcodes__phase_b__xi_turn.md
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
  - path: paper/stages/stage_248.tex
    line: 4
    role: paper_stage_tex
    excerpt: '\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.py}.  Mathematica
      audit: \StageFile{mathematica/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_mathematica_audit.wl}.}'
  - path: paper/stages/stage_248.tex
    line: 32
    role: paper_stage_tex
    excerpt: For subbarrier energy \(V(r_0)<E<V_{\rm peak}\), turning points solve
  - path: paper/stages/stage_248.tex
    line: 34
    role: paper_stage_tex
    excerpt: \label{eq:app-part08-stage248-turning-points}
  - path: paper/stages/stage_248.tex
    line: 43
    role: paper_stage_tex
    excerpt: \label{eq:app-part08-stage248-turning-transport}
  - path: paper/stages/stage_248.tex
    line: 46
    role: paper_stage_tex
    excerpt: The WKB action and transmission ratio are \cref{eq:app-part08-main-wkb-action,eq:app-part08-main-transmission-ratio}.  The
      pure-Coulomb outer turning point is exact,
  - path: paper/stages/stage_248.tex
    line: 48
    role: paper_stage_tex
    excerpt: \label{eq:app-part08-stage248-coul-turn}
  - path: paper/stages/stage_248.tex
    line: 49
    role: paper_stage_tex
    excerpt: r_{\rm turn,Coul}(E)=\frac1E,
  - path: paper/stages/stage_248.tex
    line: 85
    role: paper_stage_tex
    excerpt: \Xi_{\rm turn}(E)=\Xi_1(r_+(E)),
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 1
    role: notes_stage
    excerpt: '# Moving-Throat PDE — Stage 248: Dynamic Event-Chain Compiler from the Relaxed Stationary Barrier Front End,
      Turning Points, Threshold Speed, and WKB'
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 15
    role: notes_stage
    excerpt: 3. the standard turning-point / WKB reduction for a smooth single-peak barrier,
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 40
    role: notes_stage
    excerpt: 'What was still missing was the **dynamic event chain** that turns that stationary front end into the objects
      actually used by the Session-II scattering test:'
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 44
    role: notes_stage
    excerpt: 3. the subbarrier turning points,
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 46
    role: notes_stage
    excerpt: 5. the dynamic turning-point diagnostics carried by the same lowered branch,
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 47
    role: notes_stage
    excerpt: 6. and the exact window condition that says when the lowered branch reaches contact while the pure Coulomb reference
      still turns back.
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 52
    role: notes_stage
    excerpt: '- `scripts/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.py`'
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 56
    role: notes_stage
    excerpt: '- `mathematica/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_mathematica_audit.wl`'
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 82
    role: notes_stage
    excerpt: '- “the outer turning point moves inward,”'
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 84
    role: notes_stage
    excerpt: '- and “the turning-point branch carries a definite dynamic gradient trigger.”'
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 192
    role: notes_stage
    excerpt: So the lowered branch reaches the contact scale while the pure Coulomb branch still turns back outside it.
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 198
    role: notes_stage
    excerpt: '## 3. Turning-point compiler and subbarrier WKB event chain'
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 205
    role: notes_stage
    excerpt: '### 3.1 Turning points'
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 207
    role: notes_stage
    excerpt: The classically forbidden interval is bounded by the turning points
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 217
    role: notes_stage
    excerpt: Here \(r_+(E)\) is the outer turning point and \(r_-(E)\) is the inner turning point.
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 219
    role: notes_stage
    excerpt: Differentiating the turning-point equation gives the exact transport law
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 227
    role: notes_stage
    excerpt: '- \(V''(r_+(E))<0\), hence \(dr_+/dE<0\): the outer turning point moves inward as energy rises;'
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 228
    role: notes_stage
    excerpt: '- \(V''(r_-(E))>0\), hence \(dr_-/dE>0\): the inner turning point moves outward as energy rises.'
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 262
    role: notes_stage
    excerpt: Differentiating under the integral sign and using the vanishing boundary terms at the turning points gives
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 281
    role: notes_stage
    excerpt: 'the outer turning point is exact:'
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 284
    role: notes_stage
    excerpt: r_{\rm turn,Coul}(E)=\frac1E.
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 339
    role: notes_stage
    excerpt: Then the turning points are
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 366
    role: notes_stage
    excerpt: '## 5. Dynamic turning-point diagnostics carried forward'
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 368
    role: notes_stage
    excerpt: Even before the magnetic/helicity branch is added, the turning-point event chain naturally carries two reduced
      diagnostics that are useful later.
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 370
    role: notes_stage
    excerpt: '### 5.1 Dynamic barrier scalar at the outer turning point'
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 375
    role: notes_stage
    excerpt: \Xi_{\rm turn}(E):=\Xi_1\bigl(r_+(E)\bigr).
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 386
    role: notes_stage
    excerpt: Solving the threshold condition \(\chi_\lambda=1\) at the outer turning point gives
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 396
    role: notes_stage
    excerpt: So the same turning-point event chain also determines the first dynamic confinement-width trigger.
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 428
    role: notes_stage
    excerpt: r_{\rm turn,new}=0.39096144,
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 434
    role: notes_stage
    excerpt: with the additional turning-point diagnostics
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 436
    role: notes_stage
    excerpt: \Xi_{\rm turn}=0.34437471,
  - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
    line: 468
    role: notes_stage
    excerpt: '### 6.2 Subbarrier launch speed and Coulomb turning point'
  origin_claims:
  - parameter: Xi_turn
    introduced_at_stage: 248
    introduced_at_line: 436
    citation:
      path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
      line: 436
      excerpt: \Xi_{\rm turn}=0.34437471,
  constraints:
  - parameter: Xi_turn
    constraint_kind: published_target
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
      line: 419
      excerpt: The reported lowered-branch dynamic observables were
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings:
  - type: published_target
    severity: needs_triage
    summary: Xi_turn=0.34437471 is a 'reported lowered-branch dynamic observable' carried in as a hardcode (notes line 419-436;
      script line 247 Xi_turn_num = 0.34437471). The symbolic definition Xi_turn(E)=Xi_1(r_+(E)) is verified only structurally
      (script line 200 'carried tag sampled at r_+(E)') and the numeric value is checked only for positivity (line 311 assert
      Xi_turn_num > 0); the script never recomputes 0.34437471 from the lowered branch. So the NUMBER is an imported recorded
      figure (published_target), not derived at stage 248. The companion lambda_th=0.42826825 (notes line 438, script line
      248) is the same kind of carried hardcode.
    citations:
    - path: scripts/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.py
      line: 247
      excerpt: Xi_turn_num = 0.34437471
    - path: scripts/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.py
      line: 311
      excerpt: assert Xi_turn_num > 0
    - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
      line: 419
      excerpt: The reported lowered-branch dynamic observables were
  - type: paper_card_overclaim
    severity: low
    summary: 'The paper card presents Xi_turn(E)=Xi_1(r_+(E)) and lambda_th(E)=|E/V''(r_+(E))| as ''dynamic diagnostics carried
      forward'' (stage_248.tex line 82-89) without disclosing that their benchmark numbers are imported Session-II report
      figures rather than recomputed; the notes (line 419 ''The reported ... observables'') and script (positivity-only check)
      are clearer that these are carried hardcodes. Low severity: the card frames them as carried tags, but a reader could
      assume the numbers were derived at stage 248.'
    citations:
    - path: paper/stages/stage_248.tex
      line: 89
      excerpt: These tags are later used in the helicity and Goldilocks compilers; they do not alter the scalar event chain.
    - path: notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md
      line: 436
      excerpt: \Xi_{\rm turn}=0.34437471,
  downstream_dependents: []
  synthesis_ingested_at: '2026-06-11T19:37:44Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_27/fit_stage248_xi_turn_lambda_th_carried_hardcodes__xi_turn.yaml
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
