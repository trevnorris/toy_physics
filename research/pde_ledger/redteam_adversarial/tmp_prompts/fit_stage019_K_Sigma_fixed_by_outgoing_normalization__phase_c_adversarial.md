# Phase C Adversarial Audit Wrapper

Candidate: `fit_stage019_K_Sigma_fixed_by_outgoing_normalization`

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
id: fit_stage019_K_Sigma_fixed_by_outgoing_normalization
candidate_key: stage019_K_Sigma_fixed_by_outgoing_normalization
dry_run: false
dry_run_id: null
anchor_stages:
- 019
file_line_citations:
- path: paper/stages/stage_019.tex
  line: 27
  excerpt: The outgoing normalization fixes
  stage: 019
parameter_names:
- K_Sigma
- P_0_target
batch_id: null
status: provenance_built
status_timestamps:
  pending: '2026-06-11T08:21:51Z'
  scanned: '2026-06-11T08:21:51Z'
  provenance_built: '2026-06-11T18:25:24Z'
codex_session:
  by_parameter:
    K_Sigma:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    P_0_target:
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
  - redteam_adversarial/provenance/fit_stage019_k_sigma_fixed_by_outgoing_normalization__k_sigma.yaml
  phase_c_prompt: null
verdict:
  adversarial: null
  adjudication: null
modality_attribution:
- claim_label
modality_fragments:
- candidate_key: stage019_K_Sigma_fixed_by_outgoing_normalization
  modality: claim_label
  anchor_stage: 019
  parameter_names:
  - K_Sigma
  - P_0_target
  reason: Claims K_Sigma is "fixed" by an outgoing normalization equation containing P_{0,target}; the value is back-solved
    from a target prefactor, i.e. the canonical shape of a fit presented as fixed.
  citation:
    path: paper/stages/stage_019.tex
    line: 27
    excerpt: The outgoing normalization fixes
    stage: 019
phase_b_status: synthesis_complete
```

## Primary Sources

- paper/stages/stage_019.tex

## Provenance Slice

```yaml
- schema_version: 1
  candidate_id: fit_stage019_K_Sigma_fixed_by_outgoing_normalization
  parameter_name: K_Sigma
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage019_k_sigma_fixed_by_outgoing_normalization__phase_b__k_sigma.md
  agent_synthesis_required: true
  dry_run_fixture: false
  fixture_note: null
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - 019
  notes_source_opened: true
  source_evidence:
  - path: paper/stages/stage_019.tex
    line: 5
    role: paper_stage_tex
    excerpt: \claimstatus{This stage is \StatusExact{} for the displayed algebraic identities inside the declared parent/projection
      data, with \StatusReduced{} or \StatusExactClosure{} inherited where the local text invokes the matched reduction, mouth-kernel
      closure, parent-action packet, or one-port normal form.}
  - path: paper/stages/stage_019.tex
    line: 17
    role: paper_stage_tex
    excerpt: D_0=K_\Sigma-B_0-Z_0,\qquad
  - path: paper/stages/stage_019.tex
    line: 18
    role: paper_stage_tex
    excerpt: D_2=-(M_\Sigma+B_2+Z_2),\qquad
  - path: paper/stages/stage_019.tex
    line: 24
    role: paper_stage_tex
    excerpt: K_\Sigma=B_0+Z_0+
  - path: paper/stages/stage_019.tex
    line: 25
    role: paper_stage_tex
    excerpt: \frac{3(M_\Sigma+B_2+Z_2)^2}{B_4+Z_4}.
  - path: paper/stages/stage_019.tex
    line: 30
    role: paper_stage_tex
    excerpt: K_\Sigma=B_0+Z_0+\frac{N_0}{P_{0,\rm target}},
  - path: paper/stages/stage_019.tex
    line: 40
    role: paper_stage_tex
    excerpt: =\frac{3(M_\Sigma+B_2+Z_2)^2}{B_4+Z_4}.
  - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
    line: 5
    role: sympy_script
    excerpt: (KSigma, MSigma) to the isotropic grouped-P2 bundle moments and target surface.
  - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
    line: 25
    role: sympy_script
    excerpt: KSigma, MSigma = sp.symbols('KSigma MSigma', real=True, nonzero=True)
  - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
    line: 31
    role: sympy_script
    excerpt: D0 = KSigma - B0 - Z0
  - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
    line: 32
    role: sympy_script
    excerpt: D2 = -(MSigma + B2 + Z2)
  - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
    line: 38
    role: sympy_script
    excerpt: one_pole_numerator = D0 * (B4 + Z4) - 3 * (MSigma + B2 + Z2)**2
  - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
    line: 49
    role: sympy_script
    excerpt: K_from_one_pole = sp.simplify(sp.solve(sp.Eq(one_pole_defect, 0), KSigma)[0])
  - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
    line: 50
    role: sympy_script
    excerpt: K_from_norm = sp.simplify(sp.solve(sp.Eq(P0, P0_target), KSigma)[0])
  - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
    line: 51
    role: sympy_script
    excerpt: assert_zero('K from one-pole closed form', K_from_one_pole - (B0 + Z0 + 3*(MSigma + B2 + Z2)**2/(B4 + Z4)))
  - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
    line: 56
    role: sympy_script
    excerpt: assert_zero('compatibility equation', compatibility - (3*(MSigma + B2 + Z2)**2/(B4 + Z4) - N0/P0_target))
  - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
    line: 59
    role: sympy_script
    excerpt: N2_const_closed = 2 * N0 * (B2 + MSigma + Z2) / (B0 - KSigma + Z0)
  - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
    line: 65
    role: sympy_script
    excerpt: + 2 * B2 * MSigma
  - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
    line: 67
    role: sympy_script
    excerpt: '- 2 * B4 * KSigma'
  - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
    line: 69
    role: sympy_script
    excerpt: '- 2 * KSigma * Z4'
  - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
    line: 70
    role: sympy_script
    excerpt: + MSigma**2
  - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
    line: 71
    role: sympy_script
    excerpt: + 2 * MSigma * Z2
  - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
    line: 74
    role: sympy_script
    excerpt: ) / (B0 - KSigma + Z0) ** 2
  - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
    line: 83
    role: sympy_script
    excerpt: 'N4_const.subs({KSigma: K_from_one_pole})'
  - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
    line: 100
    role: sympy_script
    excerpt: N4_md_one_pole = -5 * (MSigma + B2 + Z2)**2 * N0 / (KSigma - B0 - Z0)**2
  - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
    line: 101
    role: sympy_script
    excerpt: assert_zero('N4 one-pole md form', N4_const.subs(KSigma, K_from_one_pole) - N4_md_one_pole.subs(KSigma, K_from_one_pole))
  - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
    line: 103
    role: sympy_script
    excerpt: M_roots = sp.solve(sp.Eq(one_pole_numerator, 0), MSigma)
  - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
    line: 104
    role: sympy_script
    excerpt: M_root_positive = sp.sqrt((KSigma - B0 - Z0) * (B4 + Z4) / 3) - (B2 + Z2)
  - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
    line: 105
    role: sympy_script
    excerpt: M_root_negative = -sp.sqrt((KSigma - B0 - Z0) * (B4 + Z4) / 3) - (B2 + Z2)
  - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
    line: 106
    role: sympy_script
    excerpt: root_gap = sp.sqrt((KSigma - B0 - Z0) * (B4 + Z4) / 3)
  - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
    line: 108
    role: sympy_script
    excerpt: '''one-pole numerator factorization in MSigma'','
  - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
    line: 109
    role: sympy_script
    excerpt: one_pole_numerator + 3 * (MSigma - M_root_positive) * (MSigma - M_root_negative),
  - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
    line: 116
    role: sympy_script
    excerpt: u2_root_positive = sp.simplify(u2.subs(MSigma, M_root_positive))
  - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
    line: 117
    role: sympy_script
    excerpt: u2_root_negative = sp.simplify(u2.subs(MSigma, M_root_negative))
  - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
    line: 129
    role: sympy_script
    excerpt: 'KSigma: sp.Integer(13),'
  - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
    line: 141
    role: sympy_script
    excerpt: 'KSigma: sp.Integer(4),'
  - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
    line: 153
    role: sympy_script
    excerpt: 'KSigma: sp.Integer(105),'
  - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
    line: 182
    role: sympy_script
    excerpt: MSigma_example = sp.integrate(mu_eta * beta**2, (w, -sp.oo, sp.oo))
  - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
    line: 183
    role: sympy_script
    excerpt: KSigma_example = sp.integrate(T_w * sp.diff(beta, w)**2 + (K_eta + 6*T_omega) * beta**2, (w, -sp.oo, sp.oo))
  - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
    line: 184
    role: sympy_script
    excerpt: assert_zero('concrete wall inertia integral', MSigma_example - sp.sqrt(sp.pi))
  origin_claims:
  - parameter: K_Sigma
    introduced_at_stage: 019
    introduced_at_line: 132
    citation:
      path: notes/stages/moving_throat_pde_stage023_full_grouped_bundle.md
      line: 210
      excerpt: '`D_(A,0) = K_A - B_(A,0) - Z_(A,0)`,'
  constraints:
  - parameter: K_Sigma
    constraint_kind: published_target
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage023_full_grouped_bundle.md
      line: 308
      excerpt: '`mhat_0^2 * N_0 / ( K - B_0 - Z_0 ) = 54 G c_s^5 / (5 a^5 c^5)`.'
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings:
  - type: derive_vs_fit_mismatch
    severity: high
    summary: 'The outgoing-normalization branch fixes K_Sigma by solving P0=N0/D0 against P0_target=54Gc_s^5/(5a^5c^5 mhat0^2)
      (stage_019 script: K_from_norm=solve(P0==P0_target,KSigma)). That P0_target is the GR quadrupole value 54/5, itself
      back-solved from the external target gamma_GR=2G/(5c^5) at stage 022. So on this branch K_Sigma is constrained by matching
      an EXTERNAL published GR coefficient, not by ledger-internal consistency. The stage card honestly labels K_Sigma''s
      two closed forms as a ''compatibility'' between the one-pole surface and the outgoing normalization (it does not claim
      K_Sigma is independently derived), but the outgoing-normalization route is a fit to a published target and should not
      be conflated with the genuinely internal one-pole route.'
    citations:
    - path: notes/stages/moving_throat_pde_stage022_grouped_p2_normalization_bridge.md
      line: 256
      excerpt: The universal GR quadrupole target is
    - path: notes/stages/moving_throat_pde_stage022_grouped_p2_normalization_bridge.md
      line: 258
      excerpt: '`gamma_GR = 2 G / (5 c^5)`.'
    - path: scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py
      line: 50
      excerpt: K_from_norm = sp.simplify(sp.solve(sp.Eq(P0, P0_target), KSigma)[0])
  - type: provenance_gap
    severity: needs_triage
    summary: Stage 019 has no per-stage notes/stages/ file and its STAGE_PROVENANCE_INDEX note source (notes/em_projected/step_17...)
      is absent from the repo, so the published_target classification of the outgoing-normalization K_Sigma branch is supported
      via the stage 022/023 notes that carry the same GR-derived 54/5 target, not via a stage-019 note.
    citations:
    - path: notes/STAGE_PROVENANCE_INDEX.md
      line: 29
      excerpt: '| 019 | `paper/stages/stage_019.tex` | `notes/em_projected/step_17_parent_throat_action_isotropic_bundle_notes.md`
        |'
    - path: paper/stages/stage_019.tex
      line: 30
      excerpt: K_\Sigma=B_0+Z_0+\frac{N_0}{P_{0,\rm target}},
  downstream_dependents:
  - '016'
  - 018
  - 019
  synthesis_ingested_at: '2026-06-11T18:25:24Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_02/fit_stage019_k_sigma_fixed_by_outgoing_normalization__k_sigma.yaml
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
