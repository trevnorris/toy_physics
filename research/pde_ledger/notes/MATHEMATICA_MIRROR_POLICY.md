# Mathematica Mirror Policy

This document defines how the PDE ledger should talk about Mathematica
coverage.

Snapshot date: `2026-05-21`

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

Stages whose red-team audit has been completed are being upgraded to native
Mathematica derivations: `EulerEquations`, `VariationalD`, `SphericalHarmonicY`,
`Coefficient`/`Series`, `ThreeJSymbol`-composed Gaunt, and so on. Those new
files land under `mathematica/` (not `scripts/`) and are listed below in the
Independent-Mirror Set. The directory convention (`.py` lives in `scripts/`,
`.wl` lives in `mathematica/`) is enforced by the red-team workflow.

## Current Independent-Mirror Set

These stages now have intentionally non-port Mathematica routes or materially
different verification structure from the SymPy side:

- `001`
  red-team batch I.1 upgraded to native `EulerEquations`/`VariationalD` plus
  `SphericalHarmonicY` for the angular-Laplacian eigenvalue check
- `002`
  red-team batch I.1 replaced transliterated extraction with native
  `SphericalHarmonicY`, `Coefficient`-based M/K extraction, and
  `EulerEquations`; 5x5 multiplet matrix checks added
- `003`
  red-team batch I.1 restructured through `DiagonalMatrix`-valued `Series`
  and a single 4x4 overlap-matrix check; also patched a multi-line `lRed = ...`
  continuation defect that had captured only kinetic terms (downstream
  unaffected, flowed through `mMat/kMat/cMat/oMat`)
- `004`
  red-team batch I.1 created native mirror for M1-M6: density-level IBP via
  combined-integrand boundary-term identity, vector Bianchi signs, Gaussian
  normalization, matched-kernel overlap, delta-source ratio
- `005`
  red-team batch I.1 created native mirror for M1-M5 using independent test
  profiles for projection-by-parts and regulator limits
- `006`
  red-team batch I.1 created native mirror via `LeviCivitaTensor[3]`/`Sum`
  for Faraday/Ampere signs, plus mediator-parity checks (antisymmetric Z
  kills the projected leak; corrects auditor's prior wrong-parity claim)
- `007`
  red-team batch I.1 created native mirror covering 11 Gaussian/regulator
  overlap claims (M1-M11) with independent integrand routes
- `008`
  red-team batch I.1 created native mirror for M1-M7 including a
  Lorentzian-Gaussian non-matched-profile numeric check; SymPy companion
  gained an independent observer-kernel test with `sigma != lambda`
- `009`
  red-team batch I.1 created native mirror with 11 manifest items M1a-M5b
  including near-throat mouth-Gaussian asymptotic via `Series` at infinity;
  SymPy erfc closed form now derived rather than typed
- `010`
  red-team batch I.1 created native mirror for the full 17-claim manifest
  using `Series`/`Coefficient`, `Solve` with uniqueness checks, and
  `ThreeJSymbol`-composed Gaunt
- `011`
  red-team batch I.1 created native mirror for 11 manifest items via
  `Series`+`Coefficient` extraction (not SymPy `(expr_lin - expr_base)/eps`)
  and `ThreeJSymbol` directly
- `012`
  red-team batch I.1 created native mirror via `Series`/`Coefficient` for
  primitive-bridge expansions; SymPy companion gained explicit negative-control
  assertions replacing earlier tautological checks
- `013`
  red-team batch I.2 created native mirror (previously no `.wl`) covering
  M1-M6 mouth-Taylor primitive expansion via `Series` on the master
  primitives `(Q - H_port ell^2)/Delta(ell)` and `(P - G_w ell^2)^2/Delta(ell)^2`
- `014`
  red-team batch I.2 created native mirror (previously no `.wl`) covering
  M1-M10 gate-bridge claims via formal `D[..., ell]` Taylor lift, with two
  negative-existence solves, two Jacobian non-vanishing determinants, and a
  sign-flip mutation guard
- `015`
  red-team batch I.2 created native mirror (previously no `.wl`) using
  `Series`+`Coefficient` (not `D[..., {eps,2}]/2`), `ThreeJSymbol` directly
  for Gaunt overlaps, and closed-form Gaussian wall-overlap evaluation;
  M1-M9 covered with wall-only gate Jacobian determinant `1/27`
- `016`
  red-team batch I.2 created native mirror (previously no `.wl`) with
  dependent-symbol `R[t,w,u,v]` declaration, independent `D[L, {eps,2}]/2`
  quadratic expansion, symbolic IBP product-rule check, `SphericalHarmonicY`
  for Y20 eigenvalue/norm/stiffness, and direct sphere integration
- `017`
  red-team batch I.2 created native mirror (previously no `.wl`) using
  `Integrate[Sin[theta] * SphericalHarmonicY * SphericalHarmonicY *
  SphericalHarmonicY, ...]` for full Wigner independence; 12 manifest claims
  M1-M12 across 23 labeled checks
- `018`
  red-team batch I.2 created native mirror (previously no `.wl`) covering
  M1-M8 with `Series`/`Coefficient` re-derivation of `u2/u4` from the pole
  expansion (not a transliteration); SymPy companion gained closed-form
  `expected_dK/expected_dM` substitutions and an `Xi1_from_expected` block
- `019`
  red-team batch I.2 created native mirror (previously no `.wl`) covering
  M1-M12 (one-pole defect identity, both closed-form `KSigma` solutions,
  `N2_const/N4_const`, Jacobian determinant `D0^3`, mutation guards)
- `020`
  red-team batch I.2 created native mirror (previously no `.wl`) defining its
  own `GauntIntegral` from first principles via `ThreeJSymbol`, with
  CamelCase naming, `Module`/`SetAttributes`, and `Solve` for the
  weak-axisymmetric packet; SymPy companion lost the `m=0` Gaunt
  short-circuit, so the m=0 lane now exercises the Wigner machinery
- `021`
  red-team batch I.2 replaced the manual EL derivation with
  `Needs["VariationalMethods``"]` + `EulerEquations[lRed, {qFun[t], aFun[t],
  wFun[t]}, t]`, patched a recurrence of the I.1 stage 003 multi-line
  `lRed = ...` continuation defect via parenthesization, and rewrote
  Sections II.2/III/V to use `LinearSolve`, an analytic-derivative route,
  and `SphericalHankelH1[2, z]` instead of mirroring SymPy
- `022`
  red-team batch I.2 switched Sections I/II/IV from `Series` extraction to
  `Solve[coeffEqs, ...]` on `Expand[ansatz*denom - num]` and replaced
  hand-typed `j2 + I*y2` with `SphericalHankelH1[2, z]`; still re-anchors
  the outgoing `l=2` coefficients through the Stage-021 exact fingerprint
  before solving the normalization product
- `023`
  red-team batch I.2 added two algebraically-distinct cross-checks
  (numerical substitution for `Z_n`/`N_n` and direct small-z Bessel
  expansion via `Series` applied to `j2 + I*y2`); SymPy companion replaced
  solver-roundtrip substitutions with hand-typed closed-form comparisons
  `N2_target_closed = 2 D2 N0/D0` and `N4_target_closed = N0(2 D0 D4 + D2^2)/D0^2`
- `059`
  uses a constructive `FindRoot` saturation route instead of symbolic replay
- `067`
  derives stationarity from the self-dual `C^2(r)=C^2(pi/r)` symmetry equation
- `069`
  closes the ordered three-zone regime algebra rather than width-only replay
- `089`
  rebuilds the Family-1 verdict from the Stage-62/63/69 formulas
- `090`
  acceptable as a narrow status-boundary replay because the checkpoint claim is
  itself an explicit carried-data verdict
- `185`
  reconstructs primitive microscopic ratios before assembling the carried
  packet
- `200`
  derives the Packet-B compiler from primitive monomials/orbit data
- `203`
  verifies the graph-composed scalar-closure / crossing route
- `218`
  rebuilds the actual support-five splice/budget ledger
- `239`
  uses the carried Stage 236/238 formulas for blind directions and orbit-lock
- `242`
  verifies orbit-lock through the direct-observable compiler
- `243`
  rebuilds the short-range kernel from the declared primitive profiles
- `248`
  has an exact symbolic route plus independent numerical stress on the
  event-chain benchmark family
- `253`
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
