# Phase C Adversarial Audit Wrapper

Candidate: `fit_stage089_pe_window_upstream_misattribution`

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
id: fit_stage089_pe_window_upstream_misattribution
candidate_key: stage089_pe_window_upstream_misattribution
dry_run: false
dry_run_id: null
anchor_stages:
- 089
file_line_citations:
- path: research/pde_ledger/scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py
  line: 61
  excerpt: '# Source: scripts/output/moving_throat_pde_stage082_*_sympy_audit.txt and the'
  stage: 089
parameter_names:
- Pe_fail_chi
- Pe_suff_chi
batch_id: null
status: provenance_built
status_timestamps:
  pending: '2026-06-11T08:21:51Z'
  scanned: '2026-06-11T08:21:51Z'
  provenance_built: '2026-06-11T18:26:39Z'
codex_session:
  by_parameter:
    Pe_fail_chi:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    Pe_suff_chi:
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
  - redteam_adversarial/provenance/fit_stage089_pe_window_upstream_misattribution__pe_fail_chi.yaml
  phase_c_prompt: null
verdict:
  adversarial: null
  adjudication: null
modality_attribution:
- existing_provenance
modality_fragments:
- candidate_key: stage089_pe_window_upstream_misattribution
  modality: existing_provenance
  anchor_stage: 089
  parameter_names:
  - Pe_suff_chi
  - Pe_fail_chi
  reason: The hardcoded literals Pe_suff_chi = 96.5285247264386 and Pe_fail_chi = 11220.5441626259 (pinned at lines 65-66,
    "not rederived here") cite the stage 082 SymPy output as upstream source, but that output file contains neither value;
    they actually originate in stage 078's output (moving_throat_pde_stage078_family1_branch_verdict_sympy_audit.txt:5-6),
    so the recorded provenance chain is mis-aimed and the literals are anchored to a file that cannot corroborate them.
  citation:
    path: research/pde_ledger/scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py
    line: 61
    excerpt: '# Source: scripts/output/moving_throat_pde_stage082_*_sympy_audit.txt and the'
    stage: 089
phase_b_status: synthesis_complete
```

## Primary Sources

- research/pde_ledger/scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py

## Provenance Slice

```yaml
- schema_version: 1
  candidate_id: fit_stage089_pe_window_upstream_misattribution
  parameter_name: Pe_fail_chi
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage089_pe_window_upstream_misattribution__phase_b__pe_fail_chi.md
  agent_synthesis_required: true
  dry_run_fixture: false
  fixture_note: null
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - 089
  notes_source_opened: true
  source_evidence:
  - path: notes/stages/moving_throat_pde_stage089_family1_minimal_isotropic_verdict.md
    line: 30
    role: notes_stage
    excerpt: '`rho_suff^(chi) ≈ 3.46622291347846,`'
  - path: notes/stages/moving_throat_pde_stage089_family1_minimal_isotropic_verdict.md
    line: 32
    role: notes_stage
    excerpt: '`rho_fail^(chi) ≈ 3.46752913273870,`'
  - path: notes/stages/moving_throat_pde_stage089_family1_minimal_isotropic_verdict.md
    line: 42
    role: notes_stage
    excerpt: '`Delta_suff := rho_suff^(chi) - 4/3`'
  - path: notes/stages/moving_throat_pde_stage089_family1_minimal_isotropic_verdict.md
    line: 45
    role: notes_stage
    excerpt: '`Delta_fail := rho_fail^(chi) - 4/3`'
  - path: scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py
    line: 59
    role: sympy_script
    excerpt: '# CARRY-FORWARD: Pe_suff_chi and Pe_fail_chi are the loading-ratio Pe values'
  - path: scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py
    line: 60
    role: sympy_script
    excerpt: '# that produce the upstream rho_suff^(chi) and rho_fail^(chi) thresholds.'
  - path: scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py
    line: 65
    role: sympy_script
    excerpt: Pe_suff_chi = sp.Float("96.5285247264386")
  - path: scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py
    line: 66
    role: sympy_script
    excerpt: Pe_fail_chi = sp.Float("11220.5441626259")
  - path: scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py
    line: 71
    role: sympy_script
    excerpt: zeta_suff = sp.simplify(zeta_F1.subs(Pe, Pe_suff_chi))
  - path: scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py
    line: 72
    role: sympy_script
    excerpt: zeta_fail = sp.simplify(zeta_F1.subs(Pe, Pe_fail_chi))
  - path: scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py
    line: 74
    role: sympy_script
    excerpt: rho_fail = sp.simplify(Q.subs(zeta, zeta_fail))
  - path: scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py
    line: 92
    role: sympy_script
    excerpt: expect_close("rho_fail vs Stage-082 quote", rho_fail, sp.Float("3.46752913273870", 30), sp.Float("1e-12", 30))
  - path: scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py
    line: 96
    role: sympy_script
    excerpt: Delta_fail = sp.N(rho_fail - rho_min, 25)
  - path: scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py
    line: 103
    role: sympy_script
    excerpt: print("rho_fail  =", sp.N(rho_fail, 25))
  - path: scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py
    line: 111
    role: sympy_script
    excerpt: print("Delta_fail =", Delta_fail)
  - path: scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py
    line: 116
    role: sympy_script
    excerpt: 'if not (rho_min < rho_suff < rho_fail < rho_max):'
  - path: scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py
    line: 117
    role: sympy_script
    excerpt: raise AssertionError("Family-1 loading-ratio ordering failed.")
  - path: scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py
    line: 128
    role: sympy_script
    excerpt: '# success margin below is the can-fail assertion establishing the boxed Output'
  origin_claims:
  - parameter: Pe_fail_chi
    introduced_at_stage: 78
    introduced_at_line: 42
    citation:
      path: notes/stages/moving_throat_pde_stage078_family1_branch_verdict.md
      line: 42
      excerpt: '`Pe_req >= Pe_fail^(chi) := Theta_w^(chi) / 3.62605617972939e-4`'
  constraints:
  - parameter: Pe_fail_chi
    constraint_kind: internal_consistency
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage080_family1_zeta_thresholds.md
      line: 47
      excerpt: '`Pe_fail^(chi) = 11220.5441626259 lambda_mu^2,`'
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings:
  - type: stale_provenance_anchor
    severity: high
    summary: 'MISATTRIBUTION CONFIRMED. Stage 089''s script anchors Pe_fail_chi=11220.5441626259 to ''Source: scripts/output/moving_throat_pde_stage082_*_sympy_audit.txt''
      and labels the rho check ''Stage-082 quote'', but Stage 082 (master_quadrupole_residual) does NOT contain this value
      anywhere in its script or output. The value first enters the ledger at Stage 078 (defined as Theta_w^(chi)/3.62605617972939e-4
      = 11220.544...) and is carried via Stage 080 (family1_zeta_thresholds, whose own script comments it as ''Stage-078 numerical
      Pe thresholds''). The 089 ''Stage-082'' anchor is a stale/wrong upstream pointer — the candidate-id ''pe_window_upstream_misattribution''
      is borne out.'
    citations:
    - path: scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py
      line: 62
      excerpt: '# Source: scripts/output/moving_throat_pde_stage082_*_sympy_audit.txt and the'
    - path: scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py
      line: 92
      excerpt: expect_close("rho_fail vs Stage-082 quote", rho_fail, sp.Float("3.46752913273870", 30), sp.Float("1e-12", 30))
    - path: scripts/moving_throat_pde_stage080_family1_zeta_thresholds_sympy_audit.py
      line: 37
      excerpt: Pe_fail_chi = sp.Float('11220.5441626259') * lam**2
  downstream_dependents:
  - 089
  synthesis_ingested_at: '2026-06-11T18:26:39Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_09/fit_stage089_pe_window_upstream_misattribution__pe_fail_chi.yaml
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
