# Independent review — S11b WL engine fix round 2 (F-WL-3b/3c controls)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11b_interface_coupling_law_mathematica_audit.wl`
(just fixed; a NARROW change to the grazing + causality tails only). It re-derives the S11b interface
coupling law and emits one tag per named object; ⛔ no VERDICT.

## What to check — narrow
Fix round 2 (directive `research/pde_ledger_v3/directives/S11b_wl_engine_fix_round2_directive.md`) claims to:
(a) make `GRAZING_MODE_CLASSIFICATION` classify the ACTUAL leading `q⁰` sound-cone block (not a `chi5`-tuned
degenerate block), so its rank/nullspace now MOVE under a form change to that block; (b) DROP the two dead
`A−A` "unrelated" controls (the grazing `thresholdUnrelatedClassification` alias and the causality
`unrelatedCausalityAggregation`), keeping the genuine controls (grazing one-sided leading-block ablation;
causality removed-record presence). Verify exactly these, plus NO regression.

## Required method — SCRIPT; FOCUSED ablation
⚠ Long Mathematica runs get spuriously killed here. Keep each kernel run SHORT/TRUNCATED — reach only the
grazing/causality tags you need; ⛔ do not run the full engine when a truncated run suffices.
- **F-WL-3c**: on a /tmp copy, FORM-ablate the leading `q⁰` sound-cone block (a sign+off-diagonal row flip)
  and confirm the emitted `GRAZING_MODE_CLASSIFICATION` rank/nullspace now MOVE (they previously stayed fixed
  via `chi5`-tuning). Confirm the classification consumes the actual block (search for `classifyLeadingBlock[
  block]`, not `[stratumBlock]`), and that NO `A−A` unrelated-classification residual is emitted. ⛔ Do not
  state the rank value.
- **F-WL-3b**: confirm `CAUSALITY_CHECK` emits NO `A−A` unrelated aggregation, and that the removed-record
  presence control still FAILS on a dropped record.
- **No regression**: F-WL-1 (`ZPERM_SLICE_MAP` static `−(lambdaA0/rhoM)`-form, ω/τ-free), F-WL-2, F-WL-3a, and
  every Scope object (impedance/regimes, added mass, `ZPERM_SLICE`, transverse, breathing, longitudinal
  dispersion, roots) must be byte-identical to the committed baseline `ec89f9df` output
  `research/pde_ledger_v3/mathematica/out/S11b_interface_coupling_law_mathematica_audit.out`.
Save every ablation script + literal stdout to named absolute paths and report them.

## Kernel discipline (binds both legs identically)
⛔ Copy the .wl to /tmp, ablate the COPY, never the working tree. ⛔ `timeout 600` on EVERY kernel; ⛔ never
raise it; ⛔ one kernel at a time (2-seat licence); ⛔ kill on RSS>6 GB, `free -h` first if you do. Save
scripts+stdout to named paths.

## Report — under ~20 lines
Numbered findings (tag, line, ablation+stdout path, why); end with a verdict: sound to commit, or repair again?
