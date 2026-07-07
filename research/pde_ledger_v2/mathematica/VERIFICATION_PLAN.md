# Mathematica Verification Plan

## Goal

Use `mathematica/` as an independent symbolic verification layer for rebuilt
ledger stages.

This is not a line-by-line translation project. A Mathematica stage audit should
derive the same claim by an independent route whenever Mathematica can
reasonably check the claim.

## Operating Principles

1. Independent derivation, not transliteration.
2. Explicit assumptions.
3. No tautological checks.
4. Saved outputs are part of the audit trail.
5. Empty ledgers are clean no-ops for the runners.

## Proposed Layout

```text
mathematica/
  VERIFICATION_PLAN.md
  run_all_audits.sh
  run_one_audit.sh
  output/
    _summary.txt
    ledger_stageNNN_*.txt
  ledger_stageNNN_*_mathematica_audit.wl
```

## Script Contract

Each future Mathematica stage script should:

- run noninteractively with `math -script mathematica/<script>.wl`;
- avoid hidden kernel session state;
- print a stage banner and the key identities being checked;
- use Mathematica primitives appropriate to the claim;
- exit nonzero on failure.

## Runner Contract

`run_all_audits.sh` discovers `ledger_stage*_mathematica_audit.wl`, writes
outputs under `mathematica/output/`, and exits 0 with a clear "nothing to do"
message when no stage scripts exist.
