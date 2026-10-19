# S9b flow-only oracle-check decision list, review pass 3: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S9b_flow_oracle_check_decision_list_review_p3_lookups.md`
beside this file (generator `_scratch/s9b_oracle/gen/dl_review_p3_lookups.sh`).

## The pass
- **Artifact.** The list after the pass-2 repair of comparison 3, sha256 `c356ab4b…` (baseline
  `directives/_legs/S9b_flow_oracle_check_decision_list_review_baseline_p3.sha256`, OK at review). The previous
  version is `3e7ec3dc`.
- **Prompt and legs.** The prompt is `directives/_legs/S9b_flow_oracle_check_decision_list_review_prompt_p3.md`
  (sha256 `168c7cdd…`): pass 2's prompt, with the repair and pass 2's findings added. The legs were fresh Codex
  (`gpt-6.1-sol`, xhigh) and Grok (`grok-4.7`, effort high), launched identically, with scratch under
  `/var/projects/_review_scratch/s9b-oracle-dl-p3-{codex,grok}/` (hard-linked at
  `_scratch/s9b_oracle/dl_review_p3_{codex,grok}_evidence/`). Reports:
  `_scratch/s9b_oracle/dl_review_p3_codex_evidence/review.md` and `_scratch/s9b_oracle/dl_review_p3_grok.txt`.
  **Both returned SOUND.**

## Findings
None from either leg.

**Sound, both legs, each with its own retrieval scripts and stdout.**
- **Earlier findings.** All five are resolved:
  - the WL exchange address is `Cases/*/ImpliedExchange`, three per flow-only tag, with no `Exchange` field;
  - comparison 3 names PY's `PY_LOCAL_S9B_VECTOR_DISPERSION` and WL's `Dispersion` field, separates the spatial
    `Metric` and the constructed `FERMAT_OBJECT`, and treats the unbound dispersion symbols as open
    identifications;
  - the full line-element residual stays `UNMATCHED`;
  - §5's Visser conditions are in the map;
  - the oracle equation numbers match the manuscripts.
- **Question at the retained order.** Comparison 1 keeps the integration constant live and the retained branch.
  Comparison 2 keeps `ρ_br` live, substitutes no oracle density, and leaves Hamilton–Lisle `UNMATCHED`.
  Comparison 3 says it checks a supplied input. The claim limit stays a conditional profile correspondence.
- **Leakage and recipe.** No expected value, sign or passing residual. Every open identification is printed with
  its own residual, so no identification can be chosen by its residual.
- **Scope.** The out-of-scope items and their reasons hold. The user decision is quoted from its source.
- **Pointers.** The transcript blobs at `503736ff`, the record at `ff3bac63`, the scoping survey at `36f59ff1` and
  the disposition at `4760e6e2` match.

**Clear.** No finding is outstanding that changes what the check computes or what may be claimed from it. The
decision-list pass was folded once (pass 1). Its object content has been reviewed until clear (passes 2 and 3).
The list is cleared at this version, and the builder may start from it.
