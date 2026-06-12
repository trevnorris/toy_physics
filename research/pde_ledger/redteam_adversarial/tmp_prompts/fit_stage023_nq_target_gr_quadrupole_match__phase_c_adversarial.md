# Phase C Adversarial Audit Wrapper

Candidate: `fit_stage023_nq_target_gr_quadrupole_match`

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
id: fit_stage023_nq_target_gr_quadrupole_match
candidate_key: stage023_nq_target_gr_quadrupole_match
dry_run: false
dry_run_id: null
anchor_stages:
- '023'
file_line_citations:
- path: paper/stages/stage_023.tex
  line: 231
  excerpt: =\frac{54Gc_s^5}{5a^5c^5}}.
  stage: '023'
parameter_names:
- N_Q_target
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
    N_Q_target:
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
  - redteam_adversarial/provenance/fit_stage023_nq_target_gr_quadrupole_match__n_q.yaml
  phase_c_prompt: null
verdict:
  adversarial: null
  adjudication: null
modality_attribution:
- numeric_literal
modality_fragments:
- candidate_key: stage023_nq_target_gr_quadrupole_match
  modality: numeric_literal
  anchor_stage: '023'
  parameter_names:
  - N_Q_target
  - mhat_0
  - gamma_GR
  reason: The boxed isotropic normalization test imposes the external GR quadrupole coefficient 2G/(5c^5) (times the port
    factor) as a target that mhat_0 is required to satisfy, i.e. the normalization is fixed by matching to GR rather than
    derived from the toy model.
  citation:
    path: paper/stages/stage_023.tex
    line: 231
    excerpt: =\frac{54Gc_s^5}{5a^5c^5}}.
    stage: '023'
phase_b_status: synthesis_complete
```

## Primary Sources

- paper/stages/stage_023.tex

## Provenance Slice

```yaml
- schema_version: 1
  candidate_id: fit_stage023_nq_target_gr_quadrupole_match
  parameter_name: N_Q
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage023_nq_target_gr_quadrupole_match__phase_b__n_q.md
  agent_synthesis_required: true
  dry_run_fixture: false
  fixture_note: null
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - '023'
  notes_source_opened: true
  source_evidence:
  - path: paper/stages/stage_023.tex
    line: 7
    role: paper_stage_tex
    excerpt: \stagefield{Purpose}{Stage~023 writes the entire grouped real \(20/21/22\) bundle in one place.  It combines
      the wall baseline, BdG support moments, projection-first Maxwell/matched mixed conservative moments, and outgoing-transfer
      moments; introduces exact weighted projectors for the grouped lanes; states the isotropic normalization test in full-bundle
      variables; and derives the first-order anisotropy transport laws.}
  - path: paper/stages/stage_023.tex
    line: 9
    role: paper_stage_tex
    excerpt: \stagefield{Inputs}{The inputs are the Stage~003 BdG moments, the Stages~004--021 projected Maxwell-to-\(Z_n,N_n\)
      packet and matched one-port conservative/transfer moments, the Stage~022 response/prefactor conversion, and the grouped
      metric \(G_{\rm grp}=\operatorname{diag}(1,2,2)\).}
  origin_claims:
  - parameter: N_Q
    introduced_at_stage: 18
    introduced_at_line: 270
    citation:
      path: notes/stages/moving_throat_pde_stage022_grouped_p2_normalization_bridge.md
      line: 270
      excerpt: '`mhat_0^2 * P_0 = 54 G c_s^5 / (5 a^5 c^5)`.'
  constraints:
  - parameter: N_Q
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
    summary: 'The stage-023 isotropic normalization test N_Q (mhat_0^2 N_0/(K-B_0-Z_0) = 54Gc_s^5/(5a^5c^5)) is the full-bundle
      restatement of the GR quadrupole match: the right-hand side 54/5 descends from gamma_GR=2G/(5c^5) back-solved at stage
      022. Stage 023 frames it as a TEST against ''the universal target'' and asks whether the true branch ''lands on'' it
      (note line 314) — i.e. a fit to an external published number, classified published_target. The card claimstatus marks
      actual branch values StatusOpen (stage_023.tex:5), so this is honestly disclosed, not overclaimed. The 54/5 target value''s
      origin is stage 022, not 023.'
    citations:
    - path: notes/stages/moving_throat_pde_stage023_full_grouped_bundle.md
      line: 314
      excerpt: '> does the true moving-throat branch produce coupled bundle data `(K, M, B_n, Z_n, N_n)` whose isotropic static
        ratio `N_0/(K-B_0-Z_0)` lands on the universal target?'
    - path: paper/stages/stage_023.tex
      line: 231
      excerpt: =\frac{54Gc_s^5}{5a^5c^5}}
  downstream_dependents:
  - '020'
  - '025'
  - '191'
  synthesis_ingested_at: '2026-06-11T18:25:24Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_02/fit_stage023_nq_target_gr_quadrupole_match__n_q.yaml
```

## Benchmarks

External-match checks must use these sourced benchmark entries. Do not adjudicate a match from model memory.

```yaml
- family_id: fam_0174_n_q
  claim: Ledger back-solves to gamma_GR = 2 G / (5 c^5), the standard GR mass-quadrupole gravitational radiation-reaction
    / luminosity normalization (the dimensionless 2/5; quadrupole luminosity coefficient G/(5 c^5)).
  value: 'Quadrupole gravitational-wave luminosity L = dE/dt = (G / 5 c^5) < d^3 Q_ij/dt^3 . d^3 Q_ij/dt^3 > using the reduced
    (traceless) mass quadrupole Q_ij. Coefficient = G/(5 c^5); equivalently G/(45 c^5) when written with the non-reduced quadrupole.
    Circular-binary luminosity carries 32/5: L = (32/5)(G/c^5) mu^2 a^4 omega^6. The radiation-reaction normalization is the
    2/5 factor.'
  source_type: textbook
  source_citation: 'P. C. Peters, ''Gravitational Radiation and the Motion of Two Point Masses'', Phys. Rev. 136, B1224 (1964),
    DOI 10.1103/PhysRev.136.B1224 (ADS bibcode 1964PhRv..136.1224P); see also M. Maggiore, ''Gravitational Waves, Vol. 1:
    Theory and Experiments'', Oxford Univ. Press (2007), ISBN 978-0-19-857074-5, Ch. 3 (quadrupole formula L = (G/5c^5)<dddQ_ij
    dddQ_ij>).'
  obtained_note: 'WebSearch on 2026-06-11 (queries on the quadrupole luminosity formula and Maggiore Vol.1) returned the explicit
    formula L_GW = (1/5)(G/c^5)<d^3 Q_ij/dt^3 d^3 Q_ij/dt^3> with the reduced/traceless quadrupole, confirming the G/(5 c^5)
    coefficient and hence the 2/5; the equivalent non-reduced form carries G/(45 c^5). The 32/5 circular-binary coefficient
    was confirmed via the Peters & Mathews search (P0 = (32/5 c^5)... for a circular orbit). The primary paper Peters 1964
    Phys. Rev. 136, B1224 was confirmed by WebSearch (Semantic Scholar + NASA ADS 1964PhRv..136.1224P + OSTI 4688457) with
    DOI 10.1103/PhysRev.136.B1224; the APS full text (link.aps.org) returned HTTP 403 (paywall) so the coefficient itself
    was confirmed from the open lecture/textbook search extractions, not the paywalled APS page. AGREEMENT: the ledger''s
    gamma_GR = 2G/(5c^5) uses exactly the standard 2/5 / G/(5c^5) quadrupole radiation-reaction normalization.'
  id: bench_62250fff83d4
  ingested_at: '2026-06-11T19:49:46Z'
  ingested_from: redteam_adversarial/provenance/_benchmarks_sourced.yaml
```

## Graph Context

```yaml
contexts: []
graph_gaps: []
```

Output a markdown adversarial report only when this is not a dry run. For dry runs, stop after confirming the prompt has enough context and record no verdict.
