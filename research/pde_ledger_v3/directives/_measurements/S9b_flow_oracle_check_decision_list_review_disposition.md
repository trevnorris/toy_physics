# S9b flow-only oracle-check decision list, review pass 1: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S9b_flow_oracle_check_decision_list_review_lookups.md` beside
this file (generator `_scratch/s9b_oracle/gen/dl_review_lookups.sh`).

## The pass
- **Artifact.** `directives/S9b_flow_oracle_check_decision_list.md`, orchestrator-written, sha256 `98c20b82…`
  (baseline `directives/_legs/S9b_flow_oracle_check_decision_list_review_baseline.sha256`, OK at review).
- **Gate.** This is a pre-builder decision list: one Codex + Grok pass, folded once (G2). It also names the object the
  check computes, and that content is reviewed until clear (spec row).
- **Prompt and legs.** The prompt is `directives/_legs/S9b_flow_oracle_check_decision_list_review_prompt.md`
  (sha256 `baaf2cf0…`). The legs were Codex (`gpt-6.1-sol`, xhigh) and Grok (`grok-4.7`, effort high), launched
  identically, with scratch under `/var/projects/_review_scratch/s9b-oracle-dl-{codex,grok}/`. Reports:
  `_scratch/s9b_oracle/dl_review_codex_evidence/review.md` and `_scratch/s9b_oracle/dl_review_grok.txt`. **Both
  returned NEEDS REVISION.**

## Findings

| # | Leg | Finding | Disposition |
|---|---|---|---|
| 1 | Codex 2, Grok 1 | **The WL flow-only exchange field is misnamed.** The list names WL's `Exchange` (28–29). In both `WL_S9B_R_FLOW_ONLY_B_IMPLIED_JN_*` tags the field is `ImpliedExchange`, inside `Cases`. The only `"Exchange"` key is on `WL_LOCAL_S9B_MASS_BALANCE`, whose velocity is the unsubstituted profile. | **Confirmed.** The lookups count 3 `"ImpliedExchange" ->` and 0 `"Exchange" ->` in each flow-only tag, and 1 `"Exchange" ->` in `WL_LOCAL_S9B_MASS_BALANCE`. A builder following the list would take the wrong operand. |
| 2 | Codex 1, Grok 2 | **The line-element residual has no S9b operand, and the optical comparison is undefined.** The list asks for Visser's acoustic line element and its difference from an engine's printed object (49–67). The flow-only tags print no metric or line element. S9b's optical object is the supplied leading-eikonal dispersion (record lines 68–73), printed as `PY_LOCAL_S9B_VECTOR_DISPERSION`, `WL_LOCAL_S9B_SUPPLIED_INPUTS`' metric and `WL_LOCAL_S9B_FERMAT_OBJECT`. Visser's eikonal equation "is blatantly insensitive to any overall multiplicative prefactor (conformal factor)" (§10, (86)). Visser states two backgrounds. §7 sets `v=√(2GM/r)`, `c=1` and uses continuity, which fixes a density law, with a line element conformal to Painlevé–Gullstrand, (53). §8 is incompressible flow, (55). | **Confirmed.** The lookups count no `metric` or `lineelement` token in any flow-only tag, find the three optical tags, and show the dispersion and Visser's sentence. A metric residual would compare a representative the conformal factor makes arbitrary. Importing §7's density law would freeze live `ρ_br` (M3). The list must name each comparison's finite operands on both sides. An optical comparison concerns the characteristic cones with the conformal factor explicit. It compares S9b's **supplied** dispersion, not an S9b result, and must say so. |
| 3 | Grok 3 | **The identification map drops conditions §5 names for Visser.** The list adds only continuity, exchange and density for Visser. §5 also names background force accounting; barotropic, inviscid, irrotational flow in three dimensions; the sink examples' extra density and equation-of-state conditions; no shear-bearing sheet or off-sheet exchange; and scalar acoustics matching a declared limit. | **Confirmed.** The lookups show §5's Visser row. Each condition is a map entry that cites a source or prints `UNMATCHED`. A velocity residual must not be read as the force balance (Euler) or the acoustic comparison. |

**Sound, both legs.**
- The baseline hash matches.
- The user decision is quoted, the record applies no oracle, and the list states no result.
- Hamilton–Lisle's equations (2) and (4) are correctly identified, and its five §5 conditions are in the map.
- The PY flow-only fields exist as named, and both engines' counting forms are quoted as printed.
- The integration constant is kept live, and so is `ρ_br`.
- No value, sign or passing residual is stated.
- The out-of-scope items and their reasons hold.

## Fold (once)
The list is folded once with all three findings. Since the fold changes what the check computes, the revised object
content gets its own review before any builder starts (spec row). The reviewed version is committed before the fold
overwrites it.
