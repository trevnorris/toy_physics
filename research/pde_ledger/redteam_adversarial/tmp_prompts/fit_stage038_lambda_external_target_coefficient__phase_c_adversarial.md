# Phase C Adversarial Audit Wrapper

Candidate: `fit_stage038_lambda_external_target_coefficient`

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
id: fit_stage038_lambda_external_target_coefficient
candidate_key: stage038_lambda_external_target_coefficient
dry_run: false
dry_run_id: null
anchor_stages:
- 038
file_line_citations:
- path: research/pde_ledger/paper/stages/stage_038.tex
  line: 28
  excerpt: \Lambda=\frac{27\pi^2Gc_s^5K_W^{\rm eff}}{20a^5c^5\mu_W}.
  stage: 038
parameter_names:
- G
- K_W_eff
- Lambda
- c_s
- mu_W
batch_id: null
status: provenance_built
status_timestamps:
  pending: '2026-06-11T08:21:51Z'
  scanned: '2026-06-11T08:21:51Z'
  provenance_built: '2026-06-11T18:25:49Z'
codex_session:
  by_parameter:
    G:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    K_W_eff:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    Lambda:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    c_s:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    mu_W:
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
  - redteam_adversarial/provenance/fit_stage038_lambda_external_target_coefficient__lambda.yaml
  phase_c_prompt: null
verdict:
  adversarial: null
  adjudication: null
modality_attribution:
- graph
modality_fragments:
- candidate_key: stage038_lambda_external_target_coefficient
  modality: graph
  anchor_stage: 038
  parameter_names:
  - Lambda
  - G
  - c_s
  - K_W_eff
  - mu_W
  reason: The 27/20 (product-law 54/5) GR-match coefficient enters the placement map as graph node EQ_FULL_BUNDLE_TARGET_SURFACE
    (expression mhat0^2 N0/D0 = 54Gc_s^5/(5a^5c^5)) whose only upstream support EQ_P0_TARGET carries status open_actual_branch_data,
    so the external target scale the band normalizes against could be back-filled to match GR rather than derived from branch
    data.
  citation:
    path: research/pde_ledger/paper/stages/stage_038.tex
    line: 28
    excerpt: \Lambda=\frac{27\pi^2Gc_s^5K_W^{\rm eff}}{20a^5c^5\mu_W}.
    stage: 038
phase_b_status: synthesis_complete
```

## Primary Sources

- research/pde_ledger/paper/stages/stage_038.tex

## Provenance Slice

```yaml
- schema_version: 1
  candidate_id: fit_stage038_lambda_external_target_coefficient
  parameter_name: Lambda
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage038_lambda_external_target_coefficient__phase_b__lambda.md
  agent_synthesis_required: true
  dry_run_fixture: false
  fixture_note: null
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - 038
  notes_source_opened: true
  source_evidence:
  - path: paper/stages/stage_038.tex
    line: 28
    role: paper_stage_tex
    excerpt: \Lambda=\frac{27\pi^2Gc_s^5K_W^{\rm eff}}{20a^5c^5\mu_W}.
  - path: paper/stages/stage_038.tex
    line: 41
    role: paper_stage_tex
    excerpt: R_{\rm target}=\frac{\Lambda(1-\epsilon_\eta)(1-\epsilon_W)^2}{Z_W(1+\rho)^2}}.
  - path: paper/stages/stage_038.tex
    line: 48
    role: paper_stage_tex
    excerpt: =\frac{8\Lambda(1-\epsilon_W)}{\pi^2}
  - path: paper/stages/stage_038.tex
    line: 51
    role: paper_stage_tex
    excerpt: Thus \((\epsilon_W,\Lambda)\) fix the product scale, \(\epsilon_\eta\) fixes the geometry lane through \(\delta\),
      and \((\epsilon_\eta,Z_W,\rho)\) redistribute the point along a fixed product curve.
  - path: notes/stages/moving_throat_pde_stage038_dimensionless_continuum_placement.md
    line: 28
    role: notes_stage
    excerpt: '`R_target = Lambda (1 - eps_eta) (1 - eps_W)^2 / [ Z_W (1 + rho)^2 ],`'
  - path: notes/stages/moving_throat_pde_stage038_dimensionless_continuum_placement.md
    line: 32
    role: notes_stage
    excerpt: '`Lambda := 27 pi^2 G c_s^5 K_W^(eff) / (20 a^5 c^5 mu_W).`'
  - path: notes/stages/moving_throat_pde_stage038_dimensionless_continuum_placement.md
    line: 38
    role: notes_stage
    excerpt: '`R_target M_mix = 8 Lambda (1 - eps_W) / pi^2`'
  - path: notes/stages/moving_throat_pde_stage038_dimensionless_continuum_placement.md
    line: 41
    role: notes_stage
    excerpt: That means three apparently independent microscopic lanes only redistribute the defect **along** a fixed product
      curve, while the mixed-sector stability ratio `eps_W` and the radiative demand scale `Lambda` set the product itself.
  - path: notes/stages/moving_throat_pde_stage038_dimensionless_continuum_placement.md
    line: 69
    role: notes_stage
    excerpt: '`Lambda := 27 pi^2 G c_s^5 K_W^(eff) / (20 a^5 c^5 mu_W).`'
  - path: notes/stages/moving_throat_pde_stage038_dimensionless_continuum_placement.md
    line: 91
    role: notes_stage
    excerpt: '`R_target = Lambda (1 - eps_eta) (1 - eps_W)^2 / [ Z_W (1 + rho)^2 ].`'
  - path: notes/stages/moving_throat_pde_stage038_dimensionless_continuum_placement.md
    line: 110
    role: notes_stage
    excerpt: '`= 8 Lambda (1 - eps_W) / pi^2`'
  - path: notes/stages/moving_throat_pde_stage038_dimensionless_continuum_placement.md
    line: 122
    role: notes_stage
    excerpt: 1. the pair `(eps_W, Lambda)` sets the product scale,
  - path: notes/stages/moving_throat_pde_stage038_dimensionless_continuum_placement.md
    line: 139
    role: notes_stage
    excerpt: '`d R_target / d eps_eta = - Lambda (1 - eps_W)^2 / [ Z_W (1 + rho)^2 ] < 0.`'
  - path: notes/stages/moving_throat_pde_stage038_dimensionless_continuum_placement.md
    line: 147
    role: notes_stage
    excerpt: '`d R_target / d eps_W = - 2 Lambda (1 - eps_eta) (1 - eps_W) / [ Z_W (1 + rho)^2 ] < 0.`'
  - path: notes/stages/moving_throat_pde_stage038_dimensionless_continuum_placement.md
    line: 156
    role: notes_stage
    excerpt: '`d R_target / d Z_W = - Lambda (1 - eps_eta) (1 - eps_W)^2 / [ Z_W^2 (1 + rho)^2 ] < 0.`'
  - path: notes/stages/moving_throat_pde_stage038_dimensionless_continuum_placement.md
    line: 164
    role: notes_stage
    excerpt: '`d R_target / d rho = - 2 Lambda (1 - eps_eta) (1 - eps_W)^2 / [ Z_W (1 + rho)^3 ].`'
  - path: notes/stages/moving_throat_pde_stage038_dimensionless_continuum_placement.md
    line: 183
    role: notes_stage
    excerpt: '`R_target M_mix = 8 Lambda (1 - eps_W) / pi^2`'
  - path: notes/stages/moving_throat_pde_stage038_dimensionless_continuum_placement.md
    line: 217
    role: notes_stage
    excerpt: '`R_target = Lambda (1 - eps_eta) (1 - eps_W)^2 / [ Z_W (1 + rho)^2 ].`'
  - path: notes/stages/moving_throat_pde_stage038_dimensionless_continuum_placement.md
    line: 228
    role: notes_stage
    excerpt: '> compute the dimensionless kernel ratios `(eps_eta, eps_W, rho, Z_W, delta_0, Lambda)` from the completed moving-throat
      PDE and check whether the resulting point lies inside the exact Stage-035/036 admissible region.'
  - path: scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py
    line: 80
    role: sympy_script
    excerpt: eps_eta, eps_W, rho, Z_W, delta0, Lambda = sp.symbols(
  - path: scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py
    line: 81
    role: sympy_script
    excerpt: '"eps_eta eps_W rho Z_W delta0 Lambda", positive=True, real=True'
  - path: scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py
    line: 90
    role: sympy_script
    excerpt: 'G: 20 * Lambda * a**5 * c_light**5 * mu_W / (27 * sp.pi**2 * c_s**5 * KWeff),'
  - path: scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py
    line: 100
    role: sympy_script
    excerpt: 'G: 20 * Lambda * a**5 * c_light**5 * mu_W / (27 * sp.pi**2 * c_s**5 * KWeff),'
  - path: scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py
    line: 115
    role: sympy_script
    excerpt: '"R_target - Lambda (1-eps_eta)(1-eps_W)^2 / [Z_W (1+rho)^2]",'
  - path: scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py
    line: 116
    role: sympy_script
    excerpt: R_dimless - Lambda * (1 - eps_eta) * (1 - eps_W) ** 2 / (Z_W * (1 + rho) ** 2),
  - path: scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py
    line: 135
    role: sympy_script
    excerpt: expect_zero("R_target M_mix - 8 Lambda (1-eps_W)/pi^2", product - 8 * Lambda * (1 - eps_W) / sp.pi**2)
  - path: scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py
    line: 141
    role: sympy_script
    excerpt: expect_zero("8 Lambda (1-eps_W)/pi^2 - NQ * KWeff(1-eps_W)/mu_W",
  - path: scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py
    line: 142
    role: sympy_script
    excerpt: 8 * Lambda * (1 - eps_W) / sp.pi**2
  - path: scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py
    line: 155
    role: sympy_script
    excerpt: R_expr = Lambda * (1 - eps_eta) * (1 - eps_W) ** 2 / (Z_W * (1 + rho) ** 2)
  - path: scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py
    line: 169
    role: sympy_script
    excerpt: expect_zero("d R / d Z_W factorization", dR_dZ + Lambda * (1 - eps_eta) * (1 - eps_W) ** 2 / (Z_W**2 * (1 + rho)
      ** 2))
  - path: scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py
    line: 171
    role: sympy_script
    excerpt: expect_zero("d R / d eps_eta factorization", dR_deps_eta + Lambda * (1 - eps_W) ** 2 / (Z_W * (1 + rho) ** 2))
  - path: scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py
    line: 173
    role: sympy_script
    excerpt: expect_zero("d R / d eps_W factorization", dR_deps_W + 2 * Lambda * (1 - eps_eta) * (1 - eps_W) / (Z_W * (1 +
      rho) ** 2))
  - path: scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py
    line: 175
    role: sympy_script
    excerpt: expect_zero("d R / d rho factorization", dR_drho + 2 * Lambda * (1 - eps_eta) * (1 - eps_W) ** 2 / (Z_W * (1
      + rho) ** 3))
  - path: scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py
    line: 178
    role: sympy_script
    excerpt: '# (0 < eps_eta < 1, 0 < eps_W < 1, 1 + rho > 0, Z_W > 0, Lambda > 0, delta0 > 0)'
  - path: scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py
    line: 192
    role: sympy_script
    excerpt: dR_dZ * Z_W**2 * (1 + rho)**2 / (Lambda * (1 - eps_eta) * (1 - eps_W)**2) + 1,
  - path: scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py
    line: 200
    role: sympy_script
    excerpt: dR_deps_eta * Z_W * (1 + rho)**2 / (Lambda * (1 - eps_W)**2) + 1,
  - path: scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py
    line: 208
    role: sympy_script
    excerpt: dR_deps_W * Z_W * (1 + rho)**2 / (2 * Lambda * (1 - eps_eta) * (1 - eps_W)) + 1,
  - path: scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py
    line: 216
    role: sympy_script
    excerpt: dR_drho * Z_W * (1 + rho)**3 / (2 * Lambda * (1 - eps_eta) * (1 - eps_W)**2) + 1,
  - path: scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py
    line: 226
    role: sympy_script
    excerpt: print("together with the product relation R_target M_mix = 8 Lambda (1-eps_W)/pi^2.")
  origin_claims:
  - parameter: Lambda
    introduced_at_stage: 038
    introduced_at_line: 69
    citation:
      path: notes/stages/moving_throat_pde_stage038_dimensionless_continuum_placement.md
      line: 69
      excerpt: '`Lambda := 27 pi^2 G c_s^5 K_W^(eff) / (20 a^5 c^5 mu_W).`'
  constraints:
  - parameter: Lambda
    constraint_kind: published_target
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage038_dimensionless_continuum_placement.md
      line: 67
      excerpt: and the radiative demand scale
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings:
  - type: published_target
    severity: high
    summary: 'Lambda is the radiative demand scale that PACKAGES the external 2.5PN GR quadrupole target into the dimensionless
      placement map. Its 27/20 coefficient is exactly (pi^2/8) x (54/5): the audit substitutes G = 20 Lambda a^5 c^5 mu_W/(27
      pi^2 c_s^5 K_W^eff) (line 90) so that the product law R_target M_mix = 8 Lambda(1-eps_W)/pi^2 equals N_Q^target K_W^eff(1-eps_W)/mu_W
      with N_Q^target = 54 G c_s^5/(5 a^5 c^5) hardcoded (lines 140-143). This is a FIT to an external published number (the
      5-coefficient 2.5PN energy-flux / 54-numerator quadrupole result), exactly the load-bearing case. NOTE: this is correctly
      DISCLOSED — the notes openly call Lambda ''the radiative demand scale'' (lines 67-69) and never claim it is internally
      derived, so there is no derive_vs_fit_mismatch; the candidate id itself names it ''external_target_coefficient''.'
    citations:
    - path: scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py
      line: 90
      excerpt: 'G: 20 * Lambda * a**5 * c_light**5 * mu_W / (27 * sp.pi**2 * c_s**5 * KWeff),'
    - path: scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py
      line: 140
      excerpt: NQ = 54 * G * c_s**5 / (5 * a**5 * c_light**5)
    - path: notes/stages/moving_throat_pde_stage038_dimensionless_continuum_placement.md
      line: 39
      excerpt: '`= 54 G c_s^5 K_W^(eff) (1 - eps_W) / (5 a^5 c^5 mu_W).`'
  downstream_dependents:
  - '032'
  - 039
  synthesis_ingested_at: '2026-06-11T18:25:49Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_04/fit_stage038_lambda_external_target_coefficient__lambda.yaml
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
