# Phase C Adversarial Audit Wrapper

Candidate: `fit_stage_100_gamma5_t`

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
id: fit_stage_100_gamma5_t
candidate_key: stage_100_gamma5_t
dry_run: false
dry_run_id: null
anchor_stages:
- '100'
file_line_citations:
- path: redteam/pass2/reports/stage_100.md
  line: 65
  role: pass2_stage_report
  stage: '100'
  excerpt: 'Stage 100 is a retarded 2.5PN factorization-ledger step. Its bottom-line claim (card line 16, boxed `\stagefield{Derivation
    ledger}`): "Full odd normalization factorizes as \(\widehat m_0^{\,2}\chi_Q N_Q=1\)." The card states the computation
    "isolates the reduced product \(\widehat m_0^{\,2}\chi_QN_Q=1\) and the canonical condition \(\chi_Q=1\)." The notes supply
    the supporting derivation: write the actual retarded branch as `Yhat_Q^ret = 3/4 + (1/4)/(1 - omega^2/Omega^2 - i chi_Q
    sigma_can omega^5) + O(omega^6)` with `sigma_can = 9/(8 Omega^5)`; extract the low-frequency tuple `Kbar_2 = Kbar_0/(4
    Omega^2)`, `Kbar_4 = Kbar_0/(4 Omega^4)`, `Gammabar_5 = chi_Q*9 Kbar_0/(32 Omega^5)`; with targets `K0_t = 64 G Omega^5/(45
    c^5)` and `Gamma5_t = 2 G/(5 c^5)` and `N_Q := K0/K0_t`, derive the even ratios `K2/K2_t = K4/K4_t = N_Q` and the odd
    ratio `Gamma5/Gamma5_t = chi_Q N_Q`; finally impose the observable condition `mhat_0^2 Gamma5 = Gamma5_t` to force `mhat_0^2
    chi_Q N_Q = 1`. The card''s three `\stagefield{Checks}` items list (i) factor-separation of the product, (ii) higher odd
    terms beginning beyond the 2.5PN coefficient, and (iii) the DtN l=2 fingerprint vs the z-expansion; the scripts explicitly
    delegate (ii) to stage 102 and (iii)/the `chi_Q=1` pin to stage 097 (carry-forward annotations), and own (i)/the closure
    here.'
- path: redteam/pass2/reports/stage_100.md
  line: 69
  role: pass2_stage_report
  stage: '100'
  excerpt: Both scripts build `Yhat_Q^ret` from the one-pole renormalized form, series-expand to O(omega^5), and read off
    K2, K4, Gamma5 from the coefficients (Gamma5 = imaginary part of the omega^5 coefficient times K0). They then construct
    the GR targets and N_Q = K0/K0_t and assert the three ratios `K2/K2_t - N_Q == 0`, `K4/K4_t - N_Q == 0`, `Gamma5/Gamma5_t
    - chi_Q*N_Q == 0`. They then impose the observable closure `mhat_0^2 Gamma5 = Gamma5_t`, form `closure_ratio = (mhat_0^2
    Gamma5 - Gamma5_t)/Gamma5_t`, and assert `closure_ratio - (mhat_0^2 chi_Q N_Q - 1) == 0`, i.e. that the closure residual
    factorizes as the headline `mhat_0^2 chi_Q N_Q - 1`. `chi_Q` is carried as a free real symbol (not pinned to 1) by design.
parameter_names:
- Gamma5_t
- Gammabar_5
- K0_t
- K2_t
- K4_t
- Kbar_0
- Kbar_2
- Kbar_4
- Yhat_Q
- chi_Q
- closure_ratio
- matched_fingerprint_value
- mhat_0
batch_id: null
status: provenance_built
status_timestamps:
  pending: '2026-06-11T08:11:49Z'
  scanned: '2026-06-11T08:11:49Z'
  provenance_built: '2026-06-11T18:51:38Z'
codex_session:
  by_parameter:
    Gamma5_t:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    Gammabar_5:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    K0_t:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    K2_t:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    K4_t:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    Kbar_0:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    Kbar_2:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    Kbar_4:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    Yhat_Q:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    chi_Q:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    closure_ratio:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    matched_fingerprint_value:
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
  - redteam_adversarial/provenance/fit_stage_100_gamma5_t__gamma5_t.yaml
  phase_c_prompt: null
verdict:
  adversarial: null
  adjudication: null
modality_attribution:
- existing_provenance
modality_fragments:
- candidate_key: stage_100_gamma5_t
  modality: existing_provenance
  anchor_stage: '100'
  parameter_names:
  - Gamma5_t
  - Gammabar_5
  - K0_t
  - K2_t
  - K4_t
  - Kbar_0
  - Kbar_2
  - Kbar_4
  - Yhat_Q
  - chi_Q
  - matched_fingerprint_value
  - mhat_0
  reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
  citation:
    path: redteam/pass2/reports/stage_100.md
    line: 65
    role: pass2_stage_report
    stage: '100'
    excerpt: 'Stage 100 is a retarded 2.5PN factorization-ledger step. Its bottom-line claim (card line 16, boxed `\stagefield{Derivation
      ledger}`): "Full odd normalization factorizes as \(\widehat m_0^{\,2}\chi_Q N_Q=1\)." The card states the computation
      "isolates the reduced product \(\widehat m_0^{\,2}\chi_QN_Q=1\) and the canonical condition \(\chi_Q=1\)." The notes
      supply the supporting derivation: write the actual retarded branch as `Yhat_Q^ret = 3/4 + (1/4)/(1 - omega^2/Omega^2
      - i chi_Q sigma_can omega^5) + O(omega^6)` with `sigma_can = 9/(8 Omega^5)`; extract the low-frequency tuple `Kbar_2
      = Kbar_0/(4 Omega^2)`, `Kbar_4 = Kbar_0/(4 Omega^4)`, `Gammabar_5 = chi_Q*9 Kbar_0/(32 Omega^5)`; with targets `K0_t
      = 64 G Omega^5/(45 c^5)` and `Gamma5_t = 2 G/(5 c^5)` and `N_Q := K0/K0_t`, derive the even ratios `K2/K2_t = K4/K4_t
      = N_Q` and the odd ratio `Gamma5/Gamma5_t = chi_Q N_Q`; finally impose the observable condition `mhat_0^2 Gamma5 = Gamma5_t`
      to force `mhat_0^2 chi_Q N_Q = 1`. The card''s three `\stagefield{Checks}` items list (i) factor-separation of the product,
      (ii) higher odd terms beginning beyond the 2.5PN coefficient, and (iii) the DtN l=2 fingerprint vs the z-expansion;
      the scripts explicitly delegate (ii) to stage 102 and (iii)/the `chi_Q=1` pin to stage 097 (carry-forward annotations),
      and own (i)/the closure here.'
- candidate_key: stage_100_gamma5_t
  modality: existing_provenance
  anchor_stage: '100'
  parameter_names:
  - Gamma5_t
  - K0_t
  - K2_t
  - K4_t
  - Yhat_Q
  - chi_Q
  - closure_ratio
  - mhat_0
  reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
  citation:
    path: redteam/pass2/reports/stage_100.md
    line: 69
    role: pass2_stage_report
    stage: '100'
    excerpt: Both scripts build `Yhat_Q^ret` from the one-pole renormalized form, series-expand to O(omega^5), and read off
      K2, K4, Gamma5 from the coefficients (Gamma5 = imaginary part of the omega^5 coefficient times K0). They then construct
      the GR targets and N_Q = K0/K0_t and assert the three ratios `K2/K2_t - N_Q == 0`, `K4/K4_t - N_Q == 0`, `Gamma5/Gamma5_t
      - chi_Q*N_Q == 0`. They then impose the observable closure `mhat_0^2 Gamma5 = Gamma5_t`, form `closure_ratio = (mhat_0^2
      Gamma5 - Gamma5_t)/Gamma5_t`, and assert `closure_ratio - (mhat_0^2 chi_Q N_Q - 1) == 0`, i.e. that the closure residual
      factorizes as the headline `mhat_0^2 chi_Q N_Q - 1`. `chi_Q` is carried as a free real symbol (not pinned to 1) by design.
phase_b_status: synthesis_complete
```

## Primary Sources

- redteam/pass2/reports/stage_100.md

## Provenance Slice

```yaml
- schema_version: 1
  candidate_id: fit_stage_100_gamma5_t
  parameter_name: Gamma5_t
  dry_run: false
  dry_run_id: null
  non_binding: false
  generated_content_kind: mechanical_evidence_bundle
  synthesis_status: complete
  agent_prompt_path: redteam_adversarial/tmp_prompts/fit_stage_100_gamma5_t__phase_b__gamma5_t.md
  agent_synthesis_required: true
  dry_run_fixture: false
  fixture_note: null
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - '100'
  notes_source_opened: true
  source_evidence:
  - path: scripts/moving_throat_pde_stage100_outgoing_normalization_factorization_sympy_audit.py
    line: 40
    role: sympy_script
    excerpt: Gamma5 = sp.simplify(sp.im(Yser.coeff(omega, 5)) * K0)
  - path: scripts/moving_throat_pde_stage100_outgoing_normalization_factorization_sympy_audit.py
    line: 45
    role: sympy_script
    excerpt: Gamma5_t = sp.simplify(2*G/(5*c**5))
  - path: scripts/moving_throat_pde_stage100_outgoing_normalization_factorization_sympy_audit.py
    line: 51
    role: sympy_script
    excerpt: print('Gamma5 =', Gamma5)
  - path: scripts/moving_throat_pde_stage100_outgoing_normalization_factorization_sympy_audit.py
    line: 55
    role: sympy_script
    excerpt: print('Gamma5/Gamma5_target - chiQ*NQ =', sp.simplify(Gamma5/Gamma5_t - chiQ*NQ_derived))
  - path: scripts/moving_throat_pde_stage100_outgoing_normalization_factorization_sympy_audit.py
    line: 58
    role: sympy_script
    excerpt: assert sp.simplify(Gamma5/Gamma5_t - chiQ*NQ_derived) == 0
  - path: scripts/moving_throat_pde_stage100_outgoing_normalization_factorization_sympy_audit.py
    line: 67
    role: sympy_script
    excerpt: closure_residual = sp.simplify(mhat0**2 * Gamma5 - Gamma5_t)
  - path: scripts/moving_throat_pde_stage100_outgoing_normalization_factorization_sympy_audit.py
    line: 68
    role: sympy_script
    excerpt: closure_ratio = sp.simplify(closure_residual / Gamma5_t)
  - path: scripts/moving_throat_pde_stage100_outgoing_normalization_factorization_sympy_audit.py
    line: 70
    role: sympy_script
    excerpt: print('closure_residual / Gamma5_target =', closure_ratio)
  origin_claims:
  - parameter: Gamma5_t
    introduced_at_stage: 100
    introduced_at_line: 69
    citation:
      path: notes/stages/moving_throat_pde_stage100_outgoing_normalization_factorization.md
      line: 69
      excerpt: '`Gammabar_5^target = 2 G/(5 c^5)`.'
  constraints:
  - parameter: Gamma5_t
    constraint_kind: published_target
    evidence_citation:
      path: notes/stages/moving_throat_pde_stage100_outgoing_normalization_factorization.md
      line: 68
      excerpt: The 2.5PN target odd coefficient is
  graph_context:
    contexts: []
    graph_gaps: []
  provenance_findings:
  - type: paper_card_overclaim
    severity: low
    summary: Gamma5_t = Gammabar_5^target = 2G/(5c^5) is the standard published GR 2.5PN (Burke-Thorne) radiation-reaction
      coefficient — an external target, correctly used here as the matching benchmark. Notes label it 'The 2.5PN target odd
      coefficient' (line 68). The script comment (line 62) describes the source-mapped factorization as matching 'the point-particle
      2.5PN coefficient ... the GR target'; this is an external anchor, not a ledger-internal derivation. Logged low to record
      the external-target nature; no actual overclaim since the card frames N_Q/chi_Q closure as Open.
    citations:
    - path: scripts/moving_throat_pde_stage100_outgoing_normalization_factorization_sympy_audit.py
      line: 45
      excerpt: Gamma5_t = sp.simplify(2*G/(5*c**5))
    - path: scripts/moving_throat_pde_stage100_outgoing_normalization_factorization_sympy_audit.py
      line: 62
      excerpt: '# (the point-particle 2.5PN coefficient with source-map factor mhat_0^2 matches'
  downstream_dependents:
  - '105'
  - '106'
  synthesis_ingested_at: '2026-06-11T18:51:38Z'
  agent_synthesis_path: redteam_adversarial/provenance/_synthesis/batch_10/fit_stage_100_gamma5_t__gamma5_t.yaml
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
