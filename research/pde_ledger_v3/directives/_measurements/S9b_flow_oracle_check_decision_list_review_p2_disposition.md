# S9b flow-only oracle-check decision list, review pass 2: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S9b_flow_oracle_check_decision_list_review_p2_lookups.md`
beside this file (generator `_scratch/s9b_oracle/gen/dl_review_p2_lookups.sh`).

## The pass
- **Artifact.** The list after its one fold, sha256 `719e098e…` (baseline
  `directives/_legs/S9b_flow_oracle_check_decision_list_review_baseline_p2.sha256`, OK at review). The previous
  version is `504fa75c`. The fold changed what the check computes, so its content is reviewed until clear (spec row).
- **Prompt and legs.** The prompt is `directives/_legs/S9b_flow_oracle_check_decision_list_review_prompt_p2.md`
  (sha256 `92cc1a16…`): pass 1's prompt, with the revision and pass 1's findings added. The legs were fresh Codex
  (`gpt-6.1-sol`, xhigh) and Grok (`grok-4.7`, effort high), with scratch under
  `/var/projects/_review_scratch/s9b-oracle-dl-p2-{codex,grok}/`. Codex returned **SOUND**; Grok returned **NEEDS
  REVISION** with two findings.

## Findings

| # | Leg | Finding | Disposition |
|---|---|---|---|
| 1 | Grok 1 | **Comparison 3's S9b operand is not the dispersion the engines print.** The list (37–39) names `WL_LOCAL_S9B_SUPPLIED_INPUTS` "(its metric)" and `WL_LOCAL_S9B_FERMAT_OBJECT` as the printed dispersion. WL prints the dispersion as that tag's `Dispersion` field, `-(cSquared*(kr^2/aMetric + kphi^2/r^2)) + (omega - kr*vLocal)^2 == 0`. Its `Metric` field is the spatial metric. `FERMAT_OBJECT` is a separate, constructed tag. PY prints the dispersion as a polynomial in `s9b_chi`, `s9b_A` and `s9b_U`, and no transcript line binds those symbols, or WL's `aMetric` and `vLocal`, to `V`, `cγ` or `g`. | **Confirmed.** The lookups show the `Dispersion` and `Metric` fields, the Fermat tag's own fields, zero bindings of `aMetric` or `vLocal`, and `s9b_chi`, `s9b_A` and `s9b_U` only on the supplied dispersion and advection tags. The flow-only restriction tag sets only the two profiles to zero. This is material my fold added. |
| 2 | Grok 2 | **Comparison 3's oracle equation numbers are wrong.** Hamilton–Lisle's line element is (1); (2) and (4) are profiles. Visser's acoustic metric is (4). §2 holds continuity (5) and Euler (6) and no line element. §7's line element is (53), and §8's velocity and line element are (54) and (55). In §10, (85) is the metric with its factor removed, (86) the eikonal equation and (87) the ray form. | **Confirmed.** The lookups find each equation number at its line in the manuscripts' text. The list's "§2" for the acoustic line element and "(85)–(86)" for "eikonal equations" are wrong. This is material my fold added. |

**Sound, both legs.**
- Pass 1's findings 1 and 3 are resolved.
- The flow-profile and exchange comparisons ask the right question at the retained order, and comparison 1's claim
  limit holds.
- `ρ_br` and the integration constant stay live, and no oracle density is substituted.
- No expected value is stated, and the scope holds.

**Not clear.** Both findings are in comparison 3, the material the fold added, and they change its operands. The
reviewed version is preserved before the repair. The repair corrects the operand addresses from the lookups above,
and the repaired comparison 3 gets its own review.
