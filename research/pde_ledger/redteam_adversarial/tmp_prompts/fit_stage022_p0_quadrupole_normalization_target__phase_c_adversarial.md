# Phase C Adversarial Audit Wrapper

Candidate: `fit_stage022_p0_quadrupole_normalization_target`

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
id: fit_stage022_p0_quadrupole_normalization_target
candidate_key: stage022_p0_quadrupole_normalization_target
dry_run: false
dry_run_id: null
anchor_stages:
- '022'
file_line_citations:
- path: research/pde_ledger/paper/stages/stage_022.tex
  line: 133
  excerpt: \widehat m_0^{\,2}P_0=\frac{54Gc_s^5}{5a^5c^5} (graph node EQ_P0_TARGET, status open_actual_branch_data, anchors
    OPEN_QUAD_NORMALIZATION)
  stage: '022'
parameter_names:
- P_0
- gamma_GR
- mhat_0
batch_id: null
status: provenance_built
status_timestamps:
  pending: '2026-06-11T08:21:51Z'
  scanned: '2026-06-11T08:21:51Z'
  provenance_built: '2026-06-11T18:25:24Z'
codex_session:
  by_parameter:
    P_0:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    gamma_GR:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    mhat_0:
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
  - redteam_adversarial/provenance/fit_stage022_p0_quadrupole_normalization_target__p0_target.yaml
  phase_c_prompt: null
verdict:
  adversarial: null
  adjudication: null
modality_attribution:
- graph
modality_fragments:
- candidate_key: stage022_p0_quadrupole_normalization_target
  modality: graph
  anchor_stage: '022'
  parameter_names:
  - mhat_0
  - P_0
  - gamma_GR
  reason: This is an explicit external-match commitment to the GR quadrupole constant 2G/(5c^5) that the graph itself classifies
    as an open gate (OPEN_QUAD_NORMALIZATION) with no derivation edge producing mhat_0^2 P_0; any later branch data chosen
    to satisfy it would be a back-filled fit rather than a derivation.
  citation:
    path: research/pde_ledger/paper/stages/stage_022.tex
    line: 133
    excerpt: \widehat m_0^{\,2}P_0=\frac{54Gc_s^5}{5a^5c^5} (graph node EQ_P0_TARGET, status open_actual_branch_data, anchors
      OPEN_QUAD_NORMALIZATION)
    stage: '022'
phase_b_status: synthesis_complete
```

## Primary Sources

- research/pde_ledger/paper/stages/stage_022.tex

## Provenance Slice

```yaml
- schema_version: 1
  candidate_id: fit_stage022_p0_quadrupole_normalization_target
  parameter_name: P0_target
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage022_p0_quadrupole_normalization_target__phase_b__p0_target.md
  agent_synthesis_required: true
  dry_run_fixture: false
  fixture_note: null
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - '022'
  notes_source_opened: true
  source_evidence:
  - path: paper/stages/stage_022.tex
    line: 118
    role: paper_stage_tex
    excerpt: The GR quadrupole target is
  - path: paper/stages/stage_022.tex
    line: 131
    role: paper_stage_tex
    excerpt: \label{eq:app-stage022-p0-target}
  - path: paper/stages/stage_022.tex
    line: 135
    role: paper_stage_tex
    excerpt: On the leading natural source-map branch \(\widehat m_0=1+O(a^2/r^2)\), this reduces to a direct target for \(P_0\).
  - path: paper/stages/stage_022.tex
    line: 137
    role: paper_stage_tex
    excerpt: \stagefield{Output}{Stage~022 outputs the normalized response formulas \eqref{eq:app-stage022-u2u4}, the grouped
      inverse map \eqref{eq:app-stage022-grouped-inverse}, the prefactor coefficients \eqref{eq:app-stage022-pref-coeffs},
      the odd coefficient \eqref{eq:app-stage022-k-gamma}, and the invariant normalization test \eqref{eq:app-stage022-p0-target}.}
  - path: paper/stages/stage_022.tex
    line: 143
    role: paper_stage_tex
    excerpt: \item The final equality \eqref{eq:app-stage022-p0-target} is a branch test, not a proof of branch realization.
  - path: paper/stages/stage_022.tex
    line: 146
    role: paper_stage_tex
    excerpt: \stagefield{Downstream use}{Stage~023 uses the Stage~022 bridge to define \(P_0=N_0/D_0\) for the full grouped
      bundle.  The 2.5PN and 4PN bridge sections use \eqref{eq:app-stage022-p0-target} as the shared outgoing-normalization
      condition.  Later same-charge and weak-axisymmetric stages use the ratio \(P_1/P_0\) as the transported prefactor slope.}
  - path: notes/stages/moving_throat_pde_stage022_grouped_p2_normalization_bridge.md
    line: 16
    role: notes_stage
    excerpt: '> lift the one-lane Stage-021 result to the grouped real `20/21/22` bundle, convert the Stage-003 and Stages-004--021
      conservative operator moments into the normalized grouped-response moments used by the 2.5PN package, and isolate the
      exact normalization product that still has to hit the universal quadrupole target.'
  - path: notes/stages/moving_throat_pde_stage022_grouped_p2_normalization_bridge.md
    line: 256
    role: notes_stage
    excerpt: The universal GR quadrupole target is
  - path: notes/stages/moving_throat_pde_stage022_grouped_p2_normalization_bridge.md
    line: 279
    role: notes_stage
    excerpt: On the natural source-map branch `mhat_0 = 1 + O(a^2/r^2)`, this reduces to the familiar 2.5PN target
  - path: notes/stages/moving_throat_pde_stage022_grouped_p2_normalization_bridge.md
    line: 311
    role: notes_stage
    excerpt: That is a much sharper theorem target than “solve the whole quadrupole branch.”
  - path: notes/stages/moving_throat_pde_stage022_grouped_p2_normalization_bridge.md
    line: 318
    role: notes_stage
    excerpt: 2. and land the invariant product `mhat_0^2 P_0` on the target value above.
  - path: notes/stages/moving_throat_pde_stage022_grouped_p2_normalization_bridge.md
    line: 334
    role: notes_stage
    excerpt: 4. and test whether the invariant normalization product `mhat_0^2 P_0` lands on the universal target.
  - path: scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py
    line: 19
    role: sympy_script
    excerpt: • the invariant quadrupole normalization product that must hit the 2.5PN target.
  - path: scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py
    line: 249
    role: sympy_script
    excerpt: '# V. Quadrupole normalization product and targets'
  - path: scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py
    line: 252
    role: sympy_script
    excerpt: 'def quadrupole_normalization_target() -> None:'
  - path: scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py
    line: 287
    role: sympy_script
    excerpt: NQ_target = sp.simplify(sp.solve(sp.Eq(mhat**2 * P0 * Gamma5_port, gamma_GR), P0)[0])
  - path: scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py
    line: 295
    role: sympy_script
    excerpt: sp.pprint(NQ_target)
  - path: scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py
    line: 298
    role: sympy_script
    excerpt: '"mhat=1 target",'
  - path: scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py
    line: 299
    role: sympy_script
    excerpt: sp.simplify(NQ_target.subs(mhat, 1) - 54 * G * c_s**5 / (5 * a**5 * c**5)),
  - path: scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py
    line: 303
    role: sympy_script
    excerpt: K2_target = sp.simplify(NQ_target * A_stage4)
  - path: scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py
    line: 304
    role: sympy_script
    excerpt: K4_target = sp.simplify(NQ_target * B_stage4)
  - path: scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py
    line: 306
    role: sympy_script
    excerpt: print("K0_target =")
  - path: scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py
    line: 307
    role: sympy_script
    excerpt: sp.pprint(NQ_target)
  - path: scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py
    line: 308
    role: sympy_script
    excerpt: print("K2_target =")
  - path: scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py
    line: 309
    role: sympy_script
    excerpt: sp.pprint(K2_target)
  - path: scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py
    line: 310
    role: sympy_script
    excerpt: print("K4_target =")
  - path: scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py
    line: 311
    role: sympy_script
    excerpt: sp.pprint(K4_target)
  - path: scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py
    line: 314
    role: sympy_script
    excerpt: '"mhat=1 K2 target",'
  - path: scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py
    line: 315
    role: sympy_script
    excerpt: sp.simplify(K2_target.subs(mhat, 1) - 6 * G * c_s**3 / (5 * a**3 * c**5)),
  - path: scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py
    line: 318
    role: sympy_script
    excerpt: '"mhat=1 K4 target",'
  - path: scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py
    line: 319
    role: sympy_script
    excerpt: sp.simplify(K4_target.subs(mhat, 1) - 8 * G * c_s / (15 * a * c**5)),
  - path: scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py
    line: 334
    role: sympy_script
    excerpt: quadrupole_normalization_target()
  - path: scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py
    line: 353
    role: sympy_script
    excerpt: print("  • the 2.5PN target is now the invariant normalization product")
  origin_claims:
  - parameter: P0_target
    introduced_at_stage: 18
    introduced_at_line: 270
    citation:
      path: notes/stages/moving_throat_pde_stage022_grouped_p2_normalization_bridge.md
      line: 270
      excerpt: '`mhat_0^2 * P_0 = 54 G c_s^5 / (5 a^5 c^5)`.'
  constraints:
  - parameter: P0_target
    constraint_kind: published_target
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage022_grouped_p2_normalization_bridge.md
      line: 256
      excerpt: The universal GR quadrupole target is
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings:
  - type: derive_vs_fit_mismatch
    severity: high
    summary: 'Stage 022 is the ORIGIN of P0_target: the note back-solves the GR target gamma_GR=2G/(5c^5) through the fingerprint
      coefficient a^5/(27c_s^5) to get P_0 = 54Gc_s^5/(5a^5c^5) on the mhat_0=1 branch. The ''54/5'' is therefore a published_target
      (fit to the external GR quadrupole normalization), not a ledger-internal derivation. Disclosure is honest (card line
      143 calls it ''a branch test, not a proof of branch realization''), so no paper_card_overclaim; but the constraint is
      published_target.'
    citations:
    - path: notes/stages/moving_throat_pde_stage022_grouped_p2_normalization_bridge.md
      line: 281
      excerpt: '`P_0 = 54 G c_s^5 / (5 a^5 c^5)`'
    - path: paper/stages/stage_022.tex
      line: 133
      excerpt: \widehat m_0^{\,2}P_0=\frac{54Gc_s^5}{5a^5c^5}}
  downstream_dependents:
  - '016'
  - 019
  - '025'
  - '191'
  synthesis_ingested_at: '2026-06-11T18:25:24Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_02/fit_stage022_p0_quadrupole_normalization_target__p0_target.yaml
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
