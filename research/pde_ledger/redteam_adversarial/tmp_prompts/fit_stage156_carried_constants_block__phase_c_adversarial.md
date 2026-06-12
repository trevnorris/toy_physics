# Phase C Adversarial Audit Wrapper

Candidate: `fit_stage156_carried_constants_block`

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
id: fit_stage156_carried_constants_block
candidate_key: stage156_carried_constants_block
dry_run: false
dry_run_id: null
anchor_stages:
- '156'
file_line_citations:
- path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
  line: 28
  excerpt: T_hat_star = 0.901484054174204
  stage: '156'
parameter_names:
- Pi_star
- Sigma_0_star
- T_hat_star
- g_star
- r_F1
batch_id: null
status: provenance_built
status_timestamps:
  pending: '2026-06-11T08:21:51Z'
  scanned: '2026-06-11T08:21:51Z'
  provenance_built: '2026-06-11T18:52:53Z'
codex_session:
  by_parameter:
    Pi_star:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    Sigma_0_star:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    T_hat_star:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    g_star:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    r_F1:
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
  - redteam_adversarial/provenance/fit_stage156_carried_constants_block__pi_star.yaml
  - redteam_adversarial/provenance/fit_stage156_carried_constants_block__sigma0_star.yaml
  - redteam_adversarial/provenance/fit_stage156_carried_constants_block__t_hat_star.yaml
  - redteam_adversarial/provenance/fit_stage156_carried_constants_block__g_star.yaml
  - redteam_adversarial/provenance/fit_stage156_carried_constants_block__r_f1.yaml
  phase_c_prompt: null
verdict:
  adversarial: null
  adjudication: null
modality_attribution:
- numeric_literal
modality_fragments:
- candidate_key: stage156_carried_constants_block
  modality: numeric_literal
  anchor_stage: '156'
  parameter_names:
  - T_hat_star
  - Sigma_0_star
  - Pi_star
  - g_star
  - r_F1
  reason: The restored canonical point (Sigma_0_can ~ 4.65103, T_hat_can ~ 1.44671) advertised by the card is found by bisection
    seeded entirely from this carried five-float block, so the "canonical" restoration inherits the provenance of typed literals
    rather than re-derived inputs.
  citation:
    path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 28
    excerpt: T_hat_star = 0.901484054174204
    stage: '156'
phase_b_status: synthesis_complete
```

## Primary Sources

- scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py

## Provenance Slice

```yaml
- schema_version: 1
  candidate_id: fit_stage156_carried_constants_block
  parameter_name: Pi_star
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage156_carried_constants_block__phase_b__pi_star.md
  agent_synthesis_required: true
  dry_run_fixture: false
  fixture_note: null
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - '156'
  notes_source_opened: true
  source_evidence:
  - path: paper/stages/stage_156.tex
    line: 16
    role: paper_stage_tex
    excerpt: Exact compensation is restored at \(\Sigma_0^{\rm can}\approx4.65103\), \(\widehat T_{m,{\rm can}}\approx1.44671\).
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 19
    role: notes_stage
    excerpt: Scanning the exact self-consistent Family-1 map shows that the fixed-point moment
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 25
    role: notes_stage
    excerpt: has a unique positive root on that analyzed interval.
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 27
    role: notes_stage
    excerpt: That numerically located root is
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 37
    role: notes_stage
    excerpt: \Sigma_0=\frac{20}{9}\widehat T_m^2,
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 42
    role: notes_stage
    excerpt: \widehat T_{m,\rm can}
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 55
    role: notes_stage
    excerpt: At that renormalized traction, the exact self-consistent fixed point satisfies
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 86
    role: notes_stage
    excerpt: (\Pi_{\rm corr},\widehat T_{m,\rm corr})
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 91
    role: notes_stage
    excerpt: (\Pi_1,\widehat T_{m,1})
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 98
    role: notes_stage
    excerpt: (\Pi_{\rm can},\widehat T_{m,\rm can})
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 106
    role: notes_stage
    excerpt: (\Pi_*,\widehat T_{m,*})
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 127
    role: notes_stage
    excerpt: \frac{\widehat T_{m,\rm can}}{\widehat T_{m,*}}-1
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 9
    role: sympy_script
    excerpt: 1. Solve the fixed point as a function of Sigma0.
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 10
    role: sympy_script
    excerpt: 2. Solve the unique root g_fp(Sigma0) = g_*.
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 11
    role: sympy_script
    excerpt: 3. Report the renormalized canonical Sigma0, T_hat, S_can, Pi_can.
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 24
    role: sympy_script
    excerpt: rF1 = 1.77799353547498
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 25
    role: sympy_script
    excerpt: g_star = 0.758035078944663
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 26
    role: sympy_script
    excerpt: Pi_star = 1.50882951349316
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 27
    role: sympy_script
    excerpt: Sigma0_star = 1.80594111095636
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 28
    role: sympy_script
    excerpt: T_hat_star = 0.901484054174204
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 60
    role: sympy_script
    excerpt: return ((gv - rF1) ** 2) / (1 + rF1**2)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 62
    role: sympy_script
    excerpt: 'def next_sigma(sig: np.ndarray, Sigma0: float) -> np.ndarray:'
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 63
    role: sympy_script
    excerpt: ph = Sigma0 * (Ts(sig) - R(sig) * Tq(sig))
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 66
    role: sympy_script
    excerpt: 'def solve_fixed_point(Sigma0: float, maxit: int = 500) -> np.ndarray:'
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 67
    role: sympy_script
    excerpt: sig = normalize(Pi_star * np.exp(-Pi_star * x))
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 69
    role: sympy_script
    excerpt: sig_new = next_sigma(sig, Sigma0)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 75
    role: sympy_script
    excerpt: 'def g_fp(Sigma0: float) -> float:'
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 76
    role: sympy_script
    excerpt: sig = solve_fixed_point(Sigma0)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 103
    role: sympy_script
    excerpt: Sigma0_can = bisect(3.0, 6.0, g_star, 55)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 104
    role: sympy_script
    excerpt: sig_can = solve_fixed_point(Sigma0_can)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 108
    role: sympy_script
    excerpt: Pi_can = Sigma0_can * (1 - R_can * S_can)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 109
    role: sympy_script
    excerpt: T_hat_can = math.sqrt(9 * Sigma0_can / 20)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 111
    role: sympy_script
    excerpt: print("Sigma0_can =", Sigma0_can)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 116
    role: sympy_script
    excerpt: print("T_hat_can  =", T_hat_can)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 119
    role: sympy_script
    excerpt: print("Sigma0 ratio - 1 =", Sigma0_can / Sigma0_star - 1)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 120
    role: sympy_script
    excerpt: print("Pi ratio - 1     =", Pi_can / Pi_star - 1)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 121
    role: sympy_script
    excerpt: print("T_hat ratio - 1  =", T_hat_can / T_hat_star - 1)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 123
    role: sympy_script
    excerpt: 'if abs(g_can - g_star) > 1e-10:'
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 127
    role: sympy_script
    excerpt: 'if abs(Sigma0_can - 4.651033550168867) > 1e-10:'
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 128
    role: sympy_script
    excerpt: raise AssertionError("Sigma0_can deviates from notes value 4.651033550168867.")
  origin_claims:
  - parameter: Pi_star
    introduced_at_stage: 156
    introduced_at_line: 106
    citation:
      path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
      line: 106
      excerpt: (\Pi_*,\widehat T_{m,*})
  constraints:
  - parameter: Pi_star
    constraint_kind: internal_consistency
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
      line: 104
      excerpt: Relative to the original canonical point
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings:
  - type: provenance_gap
    severity: low
    summary: Pi_star (the ORIGINAL frozen canonical mouth bias ~1.5088) is only quoted as a comparison reference in stage156
      notes; its first-entry derivation lives upstream (frozen Family-1 fixed point, stages 134/155). Stage156 carries it
      as a baseline for the renormalization ratios, not as a value it constrains.
    citations:
    - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
      line: 104
      excerpt: Relative to the original canonical point
    - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
      line: 26
      excerpt: Pi_star = 1.50882951349316
  downstream_dependents:
  - '157'
  - '158'
  synthesis_ingested_at: '2026-06-11T18:52:53Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_17/fit_stage156_carried_constants_block__pi_star.yaml
- schema_version: 1
  candidate_id: fit_stage156_carried_constants_block
  parameter_name: Sigma0_star
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage156_carried_constants_block__phase_b__sigma0_star.md
  agent_synthesis_required: true
  dry_run_fixture: false
  fixture_note: null
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - '156'
  notes_source_opened: true
  source_evidence:
  - path: paper/stages/stage_156.tex
    line: 16
    role: paper_stage_tex
    excerpt: Exact compensation is restored at \(\Sigma_0^{\rm can}\approx4.65103\), \(\widehat T_{m,{\rm can}}\approx1.44671\).
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 19
    role: notes_stage
    excerpt: Scanning the exact self-consistent Family-1 map shows that the fixed-point moment
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 25
    role: notes_stage
    excerpt: has a unique positive root on that analyzed interval.
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 27
    role: notes_stage
    excerpt: That numerically located root is
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 37
    role: notes_stage
    excerpt: \Sigma_0=\frac{20}{9}\widehat T_m^2,
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 42
    role: notes_stage
    excerpt: \widehat T_{m,\rm can}
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 55
    role: notes_stage
    excerpt: At that renormalized traction, the exact self-consistent fixed point satisfies
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 86
    role: notes_stage
    excerpt: (\Pi_{\rm corr},\widehat T_{m,\rm corr})
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 91
    role: notes_stage
    excerpt: (\Pi_1,\widehat T_{m,1})
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 98
    role: notes_stage
    excerpt: (\Pi_{\rm can},\widehat T_{m,\rm can})
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 106
    role: notes_stage
    excerpt: (\Pi_*,\widehat T_{m,*})
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 127
    role: notes_stage
    excerpt: \frac{\widehat T_{m,\rm can}}{\widehat T_{m,*}}-1
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 9
    role: sympy_script
    excerpt: 1. Solve the fixed point as a function of Sigma0.
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 10
    role: sympy_script
    excerpt: 2. Solve the unique root g_fp(Sigma0) = g_*.
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 11
    role: sympy_script
    excerpt: 3. Report the renormalized canonical Sigma0, T_hat, S_can, Pi_can.
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 24
    role: sympy_script
    excerpt: rF1 = 1.77799353547498
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 25
    role: sympy_script
    excerpt: g_star = 0.758035078944663
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 26
    role: sympy_script
    excerpt: Pi_star = 1.50882951349316
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 27
    role: sympy_script
    excerpt: Sigma0_star = 1.80594111095636
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 28
    role: sympy_script
    excerpt: T_hat_star = 0.901484054174204
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 60
    role: sympy_script
    excerpt: return ((gv - rF1) ** 2) / (1 + rF1**2)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 62
    role: sympy_script
    excerpt: 'def next_sigma(sig: np.ndarray, Sigma0: float) -> np.ndarray:'
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 63
    role: sympy_script
    excerpt: ph = Sigma0 * (Ts(sig) - R(sig) * Tq(sig))
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 66
    role: sympy_script
    excerpt: 'def solve_fixed_point(Sigma0: float, maxit: int = 500) -> np.ndarray:'
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 67
    role: sympy_script
    excerpt: sig = normalize(Pi_star * np.exp(-Pi_star * x))
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 69
    role: sympy_script
    excerpt: sig_new = next_sigma(sig, Sigma0)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 75
    role: sympy_script
    excerpt: 'def g_fp(Sigma0: float) -> float:'
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 76
    role: sympy_script
    excerpt: sig = solve_fixed_point(Sigma0)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 103
    role: sympy_script
    excerpt: Sigma0_can = bisect(3.0, 6.0, g_star, 55)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 104
    role: sympy_script
    excerpt: sig_can = solve_fixed_point(Sigma0_can)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 108
    role: sympy_script
    excerpt: Pi_can = Sigma0_can * (1 - R_can * S_can)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 109
    role: sympy_script
    excerpt: T_hat_can = math.sqrt(9 * Sigma0_can / 20)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 111
    role: sympy_script
    excerpt: print("Sigma0_can =", Sigma0_can)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 116
    role: sympy_script
    excerpt: print("T_hat_can  =", T_hat_can)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 119
    role: sympy_script
    excerpt: print("Sigma0 ratio - 1 =", Sigma0_can / Sigma0_star - 1)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 120
    role: sympy_script
    excerpt: print("Pi ratio - 1     =", Pi_can / Pi_star - 1)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 121
    role: sympy_script
    excerpt: print("T_hat ratio - 1  =", T_hat_can / T_hat_star - 1)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 123
    role: sympy_script
    excerpt: 'if abs(g_can - g_star) > 1e-10:'
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 127
    role: sympy_script
    excerpt: 'if abs(Sigma0_can - 4.651033550168867) > 1e-10:'
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 128
    role: sympy_script
    excerpt: raise AssertionError("Sigma0_can deviates from notes value 4.651033550168867.")
  origin_claims:
  - parameter: Sigma0_star
    introduced_at_stage: 156
    introduced_at_line: 113
    citation:
      path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
      line: 113
      excerpt: \frac{\Sigma_0^{\rm can}}{\Sigma_0^*}-1
  constraints:
  - parameter: Sigma0_star
    constraint_kind: internal_consistency
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
      line: 104
      excerpt: Relative to the original canonical point
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings:
  - type: provenance_gap
    severity: low
    summary: Sigma0_star (~1.80594, the ORIGINAL frozen-traction canonical normalization) appears in stage156 only inside
      the renormalization-cost ratio Sigma0_can/Sigma0_star - 1; it first ENTERED the ledger upstream as the frozen Family-1
      fixed point (stage 155 notes line 9). Stage156 carries it, does not constrain it. Its own derivation is the frozen fixed-point
      solve, internal to the ledger.
    citations:
    - path: notes/stages/moving_throat_pde_stage155_frozen_traction_fixedpoint.md
      line: 9
      excerpt: \approx 1.80594111095636
    - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
      line: 27
      excerpt: Sigma0_star = 1.80594111095636
  downstream_dependents:
  - '157'
  - '158'
  synthesis_ingested_at: '2026-06-11T18:52:53Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_17/fit_stage156_carried_constants_block__sigma0_star.yaml
- schema_version: 1
  candidate_id: fit_stage156_carried_constants_block
  parameter_name: T_hat_star
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage156_carried_constants_block__phase_b__t_hat_star.md
  agent_synthesis_required: true
  dry_run_fixture: false
  fixture_note: null
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - '156'
  notes_source_opened: true
  source_evidence:
  - path: paper/stages/stage_156.tex
    line: 16
    role: paper_stage_tex
    excerpt: Exact compensation is restored at \(\Sigma_0^{\rm can}\approx4.65103\), \(\widehat T_{m,{\rm can}}\approx1.44671\).
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 19
    role: notes_stage
    excerpt: Scanning the exact self-consistent Family-1 map shows that the fixed-point moment
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 25
    role: notes_stage
    excerpt: has a unique positive root on that analyzed interval.
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 27
    role: notes_stage
    excerpt: That numerically located root is
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 37
    role: notes_stage
    excerpt: \Sigma_0=\frac{20}{9}\widehat T_m^2,
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 42
    role: notes_stage
    excerpt: \widehat T_{m,\rm can}
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 55
    role: notes_stage
    excerpt: At that renormalized traction, the exact self-consistent fixed point satisfies
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 86
    role: notes_stage
    excerpt: (\Pi_{\rm corr},\widehat T_{m,\rm corr})
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 91
    role: notes_stage
    excerpt: (\Pi_1,\widehat T_{m,1})
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 98
    role: notes_stage
    excerpt: (\Pi_{\rm can},\widehat T_{m,\rm can})
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 106
    role: notes_stage
    excerpt: (\Pi_*,\widehat T_{m,*})
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 127
    role: notes_stage
    excerpt: \frac{\widehat T_{m,\rm can}}{\widehat T_{m,*}}-1
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 9
    role: sympy_script
    excerpt: 1. Solve the fixed point as a function of Sigma0.
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 10
    role: sympy_script
    excerpt: 2. Solve the unique root g_fp(Sigma0) = g_*.
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 11
    role: sympy_script
    excerpt: 3. Report the renormalized canonical Sigma0, T_hat, S_can, Pi_can.
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 24
    role: sympy_script
    excerpt: rF1 = 1.77799353547498
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 25
    role: sympy_script
    excerpt: g_star = 0.758035078944663
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 26
    role: sympy_script
    excerpt: Pi_star = 1.50882951349316
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 27
    role: sympy_script
    excerpt: Sigma0_star = 1.80594111095636
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 28
    role: sympy_script
    excerpt: T_hat_star = 0.901484054174204
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 60
    role: sympy_script
    excerpt: return ((gv - rF1) ** 2) / (1 + rF1**2)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 62
    role: sympy_script
    excerpt: 'def next_sigma(sig: np.ndarray, Sigma0: float) -> np.ndarray:'
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 63
    role: sympy_script
    excerpt: ph = Sigma0 * (Ts(sig) - R(sig) * Tq(sig))
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 66
    role: sympy_script
    excerpt: 'def solve_fixed_point(Sigma0: float, maxit: int = 500) -> np.ndarray:'
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 67
    role: sympy_script
    excerpt: sig = normalize(Pi_star * np.exp(-Pi_star * x))
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 69
    role: sympy_script
    excerpt: sig_new = next_sigma(sig, Sigma0)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 75
    role: sympy_script
    excerpt: 'def g_fp(Sigma0: float) -> float:'
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 76
    role: sympy_script
    excerpt: sig = solve_fixed_point(Sigma0)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 103
    role: sympy_script
    excerpt: Sigma0_can = bisect(3.0, 6.0, g_star, 55)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 104
    role: sympy_script
    excerpt: sig_can = solve_fixed_point(Sigma0_can)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 108
    role: sympy_script
    excerpt: Pi_can = Sigma0_can * (1 - R_can * S_can)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 109
    role: sympy_script
    excerpt: T_hat_can = math.sqrt(9 * Sigma0_can / 20)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 111
    role: sympy_script
    excerpt: print("Sigma0_can =", Sigma0_can)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 116
    role: sympy_script
    excerpt: print("T_hat_can  =", T_hat_can)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 119
    role: sympy_script
    excerpt: print("Sigma0 ratio - 1 =", Sigma0_can / Sigma0_star - 1)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 120
    role: sympy_script
    excerpt: print("Pi ratio - 1     =", Pi_can / Pi_star - 1)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 121
    role: sympy_script
    excerpt: print("T_hat ratio - 1  =", T_hat_can / T_hat_star - 1)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 123
    role: sympy_script
    excerpt: 'if abs(g_can - g_star) > 1e-10:'
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 127
    role: sympy_script
    excerpt: 'if abs(Sigma0_can - 4.651033550168867) > 1e-10:'
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 128
    role: sympy_script
    excerpt: raise AssertionError("Sigma0_can deviates from notes value 4.651033550168867.")
  origin_claims:
  - parameter: T_hat_star
    introduced_at_stage: 156
    introduced_at_line: 127
    citation:
      path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
      line: 127
      excerpt: \frac{\widehat T_{m,\rm can}}{\widehat T_{m,*}}-1
  constraints:
  - parameter: T_hat_star
    constraint_kind: internal_consistency
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
      line: 104
      excerpt: Relative to the original canonical point
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings:
  - type: provenance_gap
    severity: low
    summary: T_hat_star (~0.90148, original frozen canonical mouth traction) appears in stage156 only inside the traction-renormalization
      ratio T_hat_can/T_hat_star - 1. It is tied to Sigma0_star via the exact closure Sigma0 = (20/9) T_hat^2; first entered
      upstream at the frozen fixed point (stage 155). Carried, not constrained, at stage156.
    citations:
    - path: notes/stages/moving_throat_pde_stage155_frozen_traction_fixedpoint.md
      line: 13
      excerpt: \approx 0.901484054174204
    - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
      line: 28
      excerpt: T_hat_star = 0.901484054174204
  downstream_dependents:
  - '157'
  - '158'
  synthesis_ingested_at: '2026-06-11T18:52:53Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_17/fit_stage156_carried_constants_block__t_hat_star.yaml
- schema_version: 1
  candidate_id: fit_stage156_carried_constants_block
  parameter_name: g_star
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage156_carried_constants_block__phase_b__g_star.md
  agent_synthesis_required: true
  dry_run_fixture: false
  fixture_note: null
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - '156'
  notes_source_opened: true
  source_evidence:
  - path: paper/stages/stage_156.tex
    line: 16
    role: paper_stage_tex
    excerpt: Exact compensation is restored at \(\Sigma_0^{\rm can}\approx4.65103\), \(\widehat T_{m,{\rm can}}\approx1.44671\).
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 19
    role: notes_stage
    excerpt: Scanning the exact self-consistent Family-1 map shows that the fixed-point moment
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 25
    role: notes_stage
    excerpt: has a unique positive root on that analyzed interval.
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 27
    role: notes_stage
    excerpt: That numerically located root is
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 37
    role: notes_stage
    excerpt: \Sigma_0=\frac{20}{9}\widehat T_m^2,
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 42
    role: notes_stage
    excerpt: \widehat T_{m,\rm can}
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 55
    role: notes_stage
    excerpt: At that renormalized traction, the exact self-consistent fixed point satisfies
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 86
    role: notes_stage
    excerpt: (\Pi_{\rm corr},\widehat T_{m,\rm corr})
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 91
    role: notes_stage
    excerpt: (\Pi_1,\widehat T_{m,1})
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 98
    role: notes_stage
    excerpt: (\Pi_{\rm can},\widehat T_{m,\rm can})
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 106
    role: notes_stage
    excerpt: (\Pi_*,\widehat T_{m,*})
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 127
    role: notes_stage
    excerpt: \frac{\widehat T_{m,\rm can}}{\widehat T_{m,*}}-1
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 9
    role: sympy_script
    excerpt: 1. Solve the fixed point as a function of Sigma0.
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 10
    role: sympy_script
    excerpt: 2. Solve the unique root g_fp(Sigma0) = g_*.
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 11
    role: sympy_script
    excerpt: 3. Report the renormalized canonical Sigma0, T_hat, S_can, Pi_can.
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 24
    role: sympy_script
    excerpt: rF1 = 1.77799353547498
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 25
    role: sympy_script
    excerpt: g_star = 0.758035078944663
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 26
    role: sympy_script
    excerpt: Pi_star = 1.50882951349316
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 27
    role: sympy_script
    excerpt: Sigma0_star = 1.80594111095636
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 28
    role: sympy_script
    excerpt: T_hat_star = 0.901484054174204
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 60
    role: sympy_script
    excerpt: return ((gv - rF1) ** 2) / (1 + rF1**2)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 62
    role: sympy_script
    excerpt: 'def next_sigma(sig: np.ndarray, Sigma0: float) -> np.ndarray:'
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 63
    role: sympy_script
    excerpt: ph = Sigma0 * (Ts(sig) - R(sig) * Tq(sig))
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 66
    role: sympy_script
    excerpt: 'def solve_fixed_point(Sigma0: float, maxit: int = 500) -> np.ndarray:'
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 67
    role: sympy_script
    excerpt: sig = normalize(Pi_star * np.exp(-Pi_star * x))
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 69
    role: sympy_script
    excerpt: sig_new = next_sigma(sig, Sigma0)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 75
    role: sympy_script
    excerpt: 'def g_fp(Sigma0: float) -> float:'
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 76
    role: sympy_script
    excerpt: sig = solve_fixed_point(Sigma0)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 103
    role: sympy_script
    excerpt: Sigma0_can = bisect(3.0, 6.0, g_star, 55)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 104
    role: sympy_script
    excerpt: sig_can = solve_fixed_point(Sigma0_can)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 108
    role: sympy_script
    excerpt: Pi_can = Sigma0_can * (1 - R_can * S_can)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 109
    role: sympy_script
    excerpt: T_hat_can = math.sqrt(9 * Sigma0_can / 20)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 111
    role: sympy_script
    excerpt: print("Sigma0_can =", Sigma0_can)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 116
    role: sympy_script
    excerpt: print("T_hat_can  =", T_hat_can)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 119
    role: sympy_script
    excerpt: print("Sigma0 ratio - 1 =", Sigma0_can / Sigma0_star - 1)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 120
    role: sympy_script
    excerpt: print("Pi ratio - 1     =", Pi_can / Pi_star - 1)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 121
    role: sympy_script
    excerpt: print("T_hat ratio - 1  =", T_hat_can / T_hat_star - 1)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 123
    role: sympy_script
    excerpt: 'if abs(g_can - g_star) > 1e-10:'
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 127
    role: sympy_script
    excerpt: 'if abs(Sigma0_can - 4.651033550168867) > 1e-10:'
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 128
    role: sympy_script
    excerpt: raise AssertionError("Sigma0_can deviates from notes value 4.651033550168867.")
  origin_claims:
  - parameter: g_star
    introduced_at_stage: 122
    introduced_at_line: 63
    citation:
      path: notes/stages/moving_throat_pde_stage122_mouth_source_compensation_test.md
      line: 63
      excerpt: \mathfrak g_-^{F1}\approx 0.758035078944663,
  constraints:
  - parameter: g_star
    constraint_kind: internal_consistency
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage122_mouth_source_compensation_test.md
      line: 51
      excerpt: gives the two exact compensated mouth-coupling values
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings:
  - type: provenance_gap
    severity: low
    summary: g_star = g_*_can ~0.758035 is the compensation target the stage156 root-solve drives g_fp(Sigma0) toward; stage156
      notes (line 58) restate it as the restored canonical moment but do NOT establish its origin. The originating derivation
      is the Stage 122 minus-root of the exact compensation family at r=r_F1 (solved internal-consistency condition), which
      is itself downstream of the geometric r_F1 (=> 37/20). Carried at stage156.
    citations:
    - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
      line: 58
      excerpt: \mathfrak g_{\rm can}=\mathfrak g_*
    - path: notes/stages/moving_throat_pde_stage122_mouth_source_compensation_test.md
      line: 34
      excerpt: Stage 119 already showed that the compensated canonical outgoing branch requires
  downstream_dependents:
  - '156'
  - '157'
  - '158'
  - '159'
  synthesis_ingested_at: '2026-06-11T18:52:53Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_17/fit_stage156_carried_constants_block__g_star.yaml
- schema_version: 1
  candidate_id: fit_stage156_carried_constants_block
  parameter_name: r_F1
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage156_carried_constants_block__phase_b__r_f1.md
  agent_synthesis_required: true
  dry_run_fixture: false
  fixture_note: null
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - '156'
  notes_source_opened: true
  source_evidence:
  - path: paper/stages/stage_156.tex
    line: 16
    role: paper_stage_tex
    excerpt: Exact compensation is restored at \(\Sigma_0^{\rm can}\approx4.65103\), \(\widehat T_{m,{\rm can}}\approx1.44671\).
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 19
    role: notes_stage
    excerpt: Scanning the exact self-consistent Family-1 map shows that the fixed-point moment
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 25
    role: notes_stage
    excerpt: has a unique positive root on that analyzed interval.
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 27
    role: notes_stage
    excerpt: That numerically located root is
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 37
    role: notes_stage
    excerpt: \Sigma_0=\frac{20}{9}\widehat T_m^2,
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 42
    role: notes_stage
    excerpt: \widehat T_{m,\rm can}
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 55
    role: notes_stage
    excerpt: At that renormalized traction, the exact self-consistent fixed point satisfies
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 86
    role: notes_stage
    excerpt: (\Pi_{\rm corr},\widehat T_{m,\rm corr})
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 91
    role: notes_stage
    excerpt: (\Pi_1,\widehat T_{m,1})
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 98
    role: notes_stage
    excerpt: (\Pi_{\rm can},\widehat T_{m,\rm can})
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 106
    role: notes_stage
    excerpt: (\Pi_*,\widehat T_{m,*})
  - path: notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md
    line: 127
    role: notes_stage
    excerpt: \frac{\widehat T_{m,\rm can}}{\widehat T_{m,*}}-1
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 9
    role: sympy_script
    excerpt: 1. Solve the fixed point as a function of Sigma0.
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 10
    role: sympy_script
    excerpt: 2. Solve the unique root g_fp(Sigma0) = g_*.
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 11
    role: sympy_script
    excerpt: 3. Report the renormalized canonical Sigma0, T_hat, S_can, Pi_can.
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 24
    role: sympy_script
    excerpt: rF1 = 1.77799353547498
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 25
    role: sympy_script
    excerpt: g_star = 0.758035078944663
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 26
    role: sympy_script
    excerpt: Pi_star = 1.50882951349316
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 27
    role: sympy_script
    excerpt: Sigma0_star = 1.80594111095636
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 28
    role: sympy_script
    excerpt: T_hat_star = 0.901484054174204
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 60
    role: sympy_script
    excerpt: return ((gv - rF1) ** 2) / (1 + rF1**2)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 62
    role: sympy_script
    excerpt: 'def next_sigma(sig: np.ndarray, Sigma0: float) -> np.ndarray:'
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 63
    role: sympy_script
    excerpt: ph = Sigma0 * (Ts(sig) - R(sig) * Tq(sig))
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 66
    role: sympy_script
    excerpt: 'def solve_fixed_point(Sigma0: float, maxit: int = 500) -> np.ndarray:'
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 67
    role: sympy_script
    excerpt: sig = normalize(Pi_star * np.exp(-Pi_star * x))
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 69
    role: sympy_script
    excerpt: sig_new = next_sigma(sig, Sigma0)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 75
    role: sympy_script
    excerpt: 'def g_fp(Sigma0: float) -> float:'
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 76
    role: sympy_script
    excerpt: sig = solve_fixed_point(Sigma0)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 103
    role: sympy_script
    excerpt: Sigma0_can = bisect(3.0, 6.0, g_star, 55)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 104
    role: sympy_script
    excerpt: sig_can = solve_fixed_point(Sigma0_can)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 108
    role: sympy_script
    excerpt: Pi_can = Sigma0_can * (1 - R_can * S_can)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 109
    role: sympy_script
    excerpt: T_hat_can = math.sqrt(9 * Sigma0_can / 20)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 111
    role: sympy_script
    excerpt: print("Sigma0_can =", Sigma0_can)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 116
    role: sympy_script
    excerpt: print("T_hat_can  =", T_hat_can)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 119
    role: sympy_script
    excerpt: print("Sigma0 ratio - 1 =", Sigma0_can / Sigma0_star - 1)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 120
    role: sympy_script
    excerpt: print("Pi ratio - 1     =", Pi_can / Pi_star - 1)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 121
    role: sympy_script
    excerpt: print("T_hat ratio - 1  =", T_hat_can / T_hat_star - 1)
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 123
    role: sympy_script
    excerpt: 'if abs(g_can - g_star) > 1e-10:'
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 127
    role: sympy_script
    excerpt: 'if abs(Sigma0_can - 4.651033550168867) > 1e-10:'
  - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
    line: 128
    role: sympy_script
    excerpt: raise AssertionError("Sigma0_can deviates from notes value 4.651033550168867.")
  origin_claims:
  - parameter: r_F1
    introduced_at_stage: 121
    introduced_at_line: 70
    citation:
      path: notes/stages/moving_throat_pde_stage121_geometric_r_selection.md
      line: 70
      excerpt: \approx 1.77799353547498.
  constraints:
  - parameter: r_F1
    constraint_kind: published_target
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage121_geometric_r_selection.md
      line: 60
      excerpt: \frac{L}{a}=\frac{37}{20},
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings:
  - type: derive_vs_fit_mismatch
    severity: needs_triage
    summary: r_F1 ~1.77799 is presented in the stage156 audit script as a bare carried literal (rF1 = 1.77799353547498) and
      the notes treat it as the fixed Family-1 core parameter, but its true origin (Stage 121) is r_F1 = sqrt((12/pi^2)(37/20)^2
      - 1), i.e. a deterministic geometric function of the POSITED preferred aspect ratio L/a = 37/20. The 37/20 is a flagged
      external/carried fit target. r_F1's value is therefore only as derived as 37/20 is; classifying it (or anything built
      on it) as internal_consistency would overstate. The originating choice should be traced to L/a=37/20, not called a free
      geometric derivation.
    citations:
    - path: notes/stages/moving_throat_pde_stage121_geometric_r_selection.md
      line: 60
      excerpt: \frac{L}{a}=\frac{37}{20},
    - path: notes/stages/moving_throat_pde_stage121_geometric_r_selection.md
      line: 67
      excerpt: \sqrt{\frac{12}{\pi^2}\left(\frac{37}{20}\right)^2-1}
    - path: scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py
      line: 24
      excerpt: rF1 = 1.77799353547498
  downstream_dependents:
  - '156'
  - '157'
  - '158'
  - '159'
  - '160'
  - '161'
  synthesis_ingested_at: '2026-06-11T18:52:53Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_17/fit_stage156_carried_constants_block__r_f1.yaml
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
