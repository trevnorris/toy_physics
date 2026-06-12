# Phase C Adversarial Audit Wrapper

Candidate: `fit_stage025_mhat_rad_tuned_to_target`

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
id: fit_stage025_mhat_rad_tuned_to_target
candidate_key: stage025_mhat_rad_tuned_to_target
dry_run: false
dry_run_id: null
anchor_stages:
- '025'
file_line_citations:
- path: paper/stages/stage_025.tex
  line: 65
  excerpt: \widehat m_{\rm rad}^{\,2}
  stage: '025'
parameter_names:
- mhat_rad
batch_id: null
status: provenance_built
status_timestamps:
  pending: '2026-06-11T08:21:51Z'
  scanned: '2026-06-11T08:21:51Z'
  provenance_built: '2026-06-11T18:25:36Z'
codex_session:
  by_parameter:
    mhat_rad:
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
  - redteam_adversarial/provenance/fit_stage025_mhat_rad_tuned_to_target__mhat_rad.yaml
  phase_c_prompt: null
verdict:
  adversarial: null
  adjudication: null
modality_attribution:
- existing_provenance
modality_fragments:
- candidate_key: stage025_mhat_rad_tuned_to_target
  modality: existing_provenance
  anchor_stage: '025'
  parameter_names:
  - mhat_rad
  reason: mhat_rad^2 is a free normalization knob fixed by equating the boxed P_0 product to N_Q^target (eq:app-stage025-normalization-target;
    script sample witness mhat^2 = 162/5) — tuned to meet the external target rather than independently derived, which the
    pass-2 report records only as MATCH.
  citation:
    path: paper/stages/stage_025.tex
    line: 65
    excerpt: \widehat m_{\rm rad}^{\,2}
    stage: '025'
phase_b_status: synthesis_complete
```

## Primary Sources

- paper/stages/stage_025.tex

## Provenance Slice

```yaml
- schema_version: 1
  candidate_id: fit_stage025_mhat_rad_tuned_to_target
  parameter_name: mhat_rad
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage025_mhat_rad_tuned_to_target__phase_b__mhat_rad.md
  agent_synthesis_required: true
  dry_run_fixture: false
  fixture_note: null
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - '025'
  notes_source_opened: true
  source_evidence:
  - path: paper/stages/stage_025.tex
    line: 7
    role: paper_stage_tex
    excerpt: '\stagefield{Purpose}{Stage~025 strips the isotropic branch to the smallest radial/axial module that can carry
      the outgoing quadrupole normalization: one stable BdG support mode, one brane-like internal gauge coordinate, and one
      mixed port-active coordinate.  The stage turns the formal ratio \(N_0/D_0\) into an explicit scalar expression.}'
  - path: paper/stages/stage_025.tex
    line: 16
    role: paper_stage_tex
    excerpt: Define the radial/axial couplings
  - path: paper/stages/stage_025.tex
    line: 65
    role: paper_stage_tex
    excerpt: \widehat m_{\rm rad}^{\,2}
  - path: notes/stages/moving_throat_pde_stage025_minimal_isotropic_normalization.md
    line: 8
    role: notes_stage
    excerpt: '`mhat_ang = 1`'
  - path: notes/stages/moving_throat_pde_stage025_minimal_isotropic_normalization.md
    line: 12
    role: notes_stage
    excerpt: That means the next honest simplification is not angular. It is **radial/axial**.
  - path: notes/stages/moving_throat_pde_stage025_minimal_isotropic_normalization.md
    line: 26
    role: notes_stage
    excerpt: '`mhat_rad^2 P^2 / [ Delta (K Delta - Delta C^2 / varpi^2 - Q) ] = 54 G c_s^5 / (5 a^5 c^5)`.'
  - path: notes/stages/moving_throat_pde_stage025_minimal_isotropic_normalization.md
    line: 29
    role: notes_stage
    excerpt: It is one exact scalar equation in the radial/axial overlap amplitudes.
  - path: notes/stages/moving_throat_pde_stage025_minimal_isotropic_normalization.md
    line: 33
    role: notes_stage
    excerpt: '## 1. Minimal isotropic radial/axial data'
  - path: notes/stages/moving_throat_pde_stage025_minimal_isotropic_normalization.md
    line: 41
    role: notes_stage
    excerpt: Define the common radial/axial overlap amplitudes
  - path: notes/stages/moving_throat_pde_stage025_minimal_isotropic_normalization.md
    line: 95
    role: notes_stage
    excerpt: 'Because Stage 024 already proved `mhat_ang = 1`, the remaining source normalization is purely radial/axial:'
  - path: notes/stages/moving_throat_pde_stage025_minimal_isotropic_normalization.md
    line: 97
    role: notes_stage
    excerpt: '`mhat_0 = mhat_rad`.'
  - path: notes/stages/moving_throat_pde_stage025_minimal_isotropic_normalization.md
    line: 101
    role: notes_stage
    excerpt: '`mhat_rad^2 P^2 / [ Delta ( K Delta - Delta C^2 / varpi^2 - Q ) ] = 54 G c_s^5 / (5 a^5 c^5)`.'
  - path: notes/stages/moving_throat_pde_stage025_minimal_isotropic_normalization.md
    line: 195
    role: notes_stage
    excerpt: On the minimal isotropic branch it needs only to determine the radial/axial amplitudes entering
  - path: notes/stages/moving_throat_pde_stage025_minimal_isotropic_normalization.md
    line: 205
    role: notes_stage
    excerpt: and `mhat_rad`.
  - path: notes/stages/moving_throat_pde_stage025_minimal_isotropic_normalization.md
    line: 209
    role: notes_stage
    excerpt: '`mhat_rad^2 P^2 / [ Delta (K Delta - Delta X - Q) ] ?= 54 G c_s^5 / (5 a^5 c^5)`.'
  - path: notes/stages/moving_throat_pde_stage025_minimal_isotropic_normalization.md
    line: 216
    role: notes_stage
    excerpt: '- if not, the failure is not mysterious — it is in one of the radial/axial quantities `X`, `Delta`, `Q`, `P`,
      or `mhat_rad`.'
  - path: notes/stages/moving_throat_pde_stage025_minimal_isotropic_normalization.md
    line: 232
    role: notes_stage
    excerpt: '### Radial/axial layer'
  - path: notes/stages/moving_throat_pde_stage025_minimal_isotropic_normalization.md
    line: 236
    role: notes_stage
    excerpt: '`mhat_rad^2 P^2 / [ Delta (K Delta - Delta C^2 / varpi^2 - Q) ] = 54 G c_s^5 / (5 a^5 c^5)`.'
  - path: scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py
    line: 20
    role: sympy_script
    excerpt: problem into one explicit scalar equation in radial/axial overlap amplitudes.
  - path: scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py
    line: 58
    role: sympy_script
    excerpt: mhat = sp.symbols("mhat", positive=True, real=True)
  - path: scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py
    line: 122
    role: sympy_script
    excerpt: '# imply mhat^2*P0 = gamma_GR/Gamma5_port = 54*G*c_s^5/(5*a^5*c^5).'
  - path: scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py
    line: 124
    role: sympy_script
    excerpt: equation_residual = sp.simplify(mhat**2 * P0_compact - target)
  - path: scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py
    line: 129
    role: sympy_script
    excerpt: '# Solvability check: mhat^2 = target / P0_compact must be positive on the'
  - path: scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py
    line: 132
    role: sympy_script
    excerpt: mhat_sq = sp.simplify(target / P0_compact)
  - path: scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py
    line: 135
    role: sympy_script
    excerpt: mhat_sq_at_sample = sp.nsimplify(mhat_sq.subs(sample_with_units))
  - path: scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py
    line: 136
    role: sympy_script
    excerpt: print(f"mhat^2 on sample = {mhat_sq_at_sample}")
  - path: scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py
    line: 137
    role: sympy_script
    excerpt: 'if mhat_sq_at_sample <= 0:'
  - path: scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py
    line: 138
    role: sympy_script
    excerpt: 'raise AssertionError(f"mhat^2 on sample is not positive: {mhat_sq_at_sample}")'
  - path: scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py
    line: 146
    role: sympy_script
    excerpt: residual_at_Pzero = sp.simplify((mhat**2 * P0_compact - target).subs(P_zero_sub))
  - path: scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py
    line: 147
    role: sympy_script
    excerpt: print(f"(mhat^2*P0 - target) at P=0 = {residual_at_Pzero}")
  - path: scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py
    line: 148
    role: sympy_script
    excerpt: expect_zero("(mhat^2*P0 - target) at P=0 equals -target", residual_at_Pzero + target)
  origin_claims:
  - parameter: mhat_rad
    introduced_at_stage: 25
    introduced_at_line: 97
    citation:
      path: notes/stages/moving_throat_pde_stage025_minimal_isotropic_normalization.md
      line: 97
      excerpt: '`mhat_0 = mhat_rad`.'
  constraints:
  - parameter: mhat_rad
    constraint_kind: published_target
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage025_minimal_isotropic_normalization.md
      line: 101
      excerpt: '`mhat_rad^2 P^2 / [ Delta ( K Delta - Delta C^2 / varpi^2 - Q ) ] = 54 G c_s^5 / (5 a^5 c^5)`.'
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings:
  - type: derive_vs_fit_mismatch
    severity: needs_triage
    summary: mhat_rad is the free amplitude solved so the one-mode prefactor hits the external GR target 54/5 G c_s^5/(5 a^5
      c^5); the audit script back-solves mhat^2 = target/P0_compact. The card frames the stage as an 'exact closure' but mhat_rad
      itself is the tuning knob fitted to the external GR quadrupole coefficient (gamma_GR = 2 G/(5 c^5)), not a ledger-internal
      derivation. Whether the branch realizes the required value is flagged StatusOpen in the card, consistent with a fit,
      not a derivation.
    citations:
    - path: scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py
      line: 132
      excerpt: mhat_sq = sp.simplify(target / P0_compact)
    - path: notes/stages/moving_throat_pde_stage022_grouped_p2_normalization_bridge.md
      line: 258
      excerpt: '`gamma_GR = 2 G / (5 c^5)`.'
    - path: paper/stages/stage_025.tex
      line: 5
      excerpt: Whether the actual branch realizes the formula at the required value is \StatusOpen{}.
  downstream_dependents:
  - '026'
  - '030'
  - '031'
  - '033'
  synthesis_ingested_at: '2026-06-11T18:25:36Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_03/fit_stage025_mhat_rad_tuned_to_target__mhat_rad.yaml
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
