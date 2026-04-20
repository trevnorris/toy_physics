# Mathematica Mirror Policy

This document defines how the PDE ledger should talk about Mathematica
coverage.

Snapshot date: `2026-04-21`

## Rule

`Mathematica audit present` is an execution-coverage fact, not an independence
claim.

Unless a stage is explicitly listed below as an `independent mirror`, the
Mathematica file should be treated as:

- secondary execution coverage,
- useful for CAS replay and environment drift detection,
- but not by itself evidence of a genuinely independent second derivation.

This policy exists because much of the Mathematica corpus was generated as a
port of the SymPy logic. Those mirrors are still useful, but they should not be
described as independent corroboration.

## Current Independent-Mirror Set

These stages now have intentionally non-port Mathematica routes or materially
different verification structure from the SymPy side:

- `003`
  acceptable dual-CAS checkpoint coverage after the shared stress-harness
  hardening; no active independence defect remains
- `005`
  re-anchors the outgoing `l=2` coefficients through the Stage-4 exact
  fingerprint before solving the normalization product
- `042`
  uses a constructive `FindRoot` saturation route instead of symbolic replay
- `050`
  derives stationarity from the self-dual `C^2(r)=C^2(pi/r)` symmetry equation
- `052`
  closes the ordered three-zone regime algebra rather than width-only replay
- `072`
  rebuilds the Family-1 verdict from the Stage-62/63/69 formulas
- `073`
  acceptable as a narrow status-boundary replay because the checkpoint claim is
  itself an explicit carried-data verdict
- `168`
  reconstructs primitive microscopic ratios before assembling the carried
  packet
- `183`
  derives the Packet-B compiler from primitive monomials/orbit data
- `186`
  verifies the graph-composed scalar-closure / crossing route
- `201`
  rebuilds the actual support-five splice/budget ledger
- `222`
  uses the carried Stage-219/221 formulas for blind directions and orbit-lock
- `225`
  verifies orbit-lock through the direct-observable compiler
- `226`
  rebuilds the short-range kernel from the declared primitive profiles
- `231`
  has an exact symbolic route plus independent numerical stress on the
  event-chain benchmark family
- `236`
  has an exact symbolic route plus independent numerical stress on the
  material-threshold screening family

## Default Disposition For All Other Mathematica Files

If a stage is not in the list above:

- count it as `Mathematica present`,
- do not describe it as `independent dual-CAS support`,
- and rely on the stage review note or checkpoint trust audit to say whether the
  current mirror is good enough for the actual claim.

## Practical Use

- `STAGE_VERIFICATION_COVERAGE.md`
  should keep reporting raw Mathematica presence counts, but it must now say
  explicitly that those are not independence counts.
- `CHECKPOINT_TRUST_AUDIT.md`
  can still call a stage `strong` if the theorem path is exact and the current
  mirror quality is appropriate for the stated claim.
- Future widening work should upgrade only the load-bearing subset instead of
  trying to make every Mathematica file fully independent.
