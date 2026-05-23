# Mathematica Mirror Policy

This document defines how the PDE ledger should talk about Mathematica
coverage.

Snapshot date: `2026-05-22` (batch III.2 close)

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

Red-team batches must explicitly screen each `.wl` for line-by-line
correspondence to its `.py` counterpart as the first audit step. If the `.wl`
mirrors SymPy primitives one-for-one -- shared local variable names, identical
section ordering, hand-typed polynomial answers, `pairings`-style recursion,
or `_expected` self-substitutions copied from the SymPy companion -- it must
be rewritten to a Mathematica-native primitive (`Integrate`, `Eigenvalues`,
`LinearSolve`, `Series`+`Coefficient`, `Solve`, `Factor`/`Apart`,
`EulerEquations`, `SphericalHarmonicY`, `ThreeJSymbol`, `SphericalHankelH1`)
before the batch closes. Batch II.1 found `mathematica_transliteration` on
every single one of its 13 stages; batch III.1 found it on 10 of 12 (with
the remaining two — 042 and 048 — passing transliteration screening only
because their closed-form-identity structure or independent `Solve`/
`Series` route already broke the line-by-line correspondence); batch III.2
found it on 6 of 12 (049, 051, 052, 058, 059, 060), the lower share due
to a higher concentration of pure tautology / hardcoded findings in the
other 5 dirty stages (050, 053, 054, 055, 057) plus one clean-on-first-read
stage (056) whose `Limit[..., {Pe, Infinity}]` and `Series.removeO` paths
already broke line-by-line correspondence. Treat transliteration as the
default expectation, not an exceptional finding.

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
- `024`
  red-team batch II.1 replaced Wick-pair `pairings` recursion with
  `Integrate[..., {theta, 0, Pi}, {phi, 0, 2 Pi}]` over Cartesian `n[i]`
  components for the angular moments; removed SymPy-named shorthands
  (`deltaPair`/`sPair`/`qPair`/`hPair`/`pPair`) and added an `xLane[lam_]`
  parameterizer in place of pre-substituted lane forms
- `025`
  red-team batch II.1 replaced transliterated algebra with `Factor`,
  `Apart[Together[..., k]]`, `Limit`-based derivatives, and `Reduce`
  solvability checks; SymPy companion added numerical sample-point checks
  on `P0/Delta/D0/mhat^2/dP0` and reanchored the `54/5` overlap target to
  Stage 023's exact derivation
- `026`
  red-team batch II.1 added two algebraically-distinct routes for the
  overlap law (indefinite-integral + boundary-evaluation path vs typed
  analytic short form); deleted `_expected` self-substitution rebuilds
  and the `K_req_expected` solver round-trip
- `027`
  red-team batch II.1 built `gEta = -tW*D[chi,{s,2}] + (kEta+6*tOmega)*chi`
  and evaluated `kGeo = FullSimplify[Integrate[chi*gEta, {s,0,l}], ...]`
  instead of hard-coding the answer; output prints `kGeo` in `Cos[2*theta]`
  canonical form rather than `sin^2(theta)`, proving the integral was
  actually evaluated
- `028`
  red-team batch II.1 replaced typed eigenvalue answers with
  `Eigenvalues[kEff]` sum/product checks and `Solve[detEff == 0, alpha]`
- `029`
  red-team batch II.1 used sequential elimination (phi -> W -> U) for the
  Schur block and `Eigensystem[keffAl]` for `kappa_sel` instead of
  transliterated linear algebra
- `030`
  red-team batch II.1 used `Eigenvalues[mMat]` on the explicit 2x2 wall
  block instead of hand-typed `lam_+`/`lam_-`
- `031`
  red-team batch II.1 derived `radcrit` from `T0^2*R^2.subs(alpha,
  alpha_crit)` instead of a hand-typed 9-term polynomial; SymPy companion
  replaced abstract `sp.Function("S")`/`sp.Function("L")` derivations with
  the physical `s_-`/`lam_-` symbols
- `032`
  red-team batch II.1 used `LinearSolve[kInt . y == bMat^T . z, ...]` to
  derive the Schur matrix from coefficients, with a `delta_kappa^2 +
  4*Kprod = sigma^2` identity check covering the previously-uncovered
  interior region
- `033`
  red-team batch II.1 added a numerical cross-check at two rational rule
  sets in addition to the symbolic route, and replaced the Mathematica
  `k0Onset` hardcoded form with `Solve[n0Mic == NQ, K0]`
- `034`
  red-team batch II.1 replaced the linear-solve self-check with
  `solve(gB+alpha_mix==alpha_x, lam)` in the original lambda variable,
  then substituted `lam = A-x`
- `035`
  red-team batch II.1 preserved target literals as the claim under test
  but switched `expectZero` LHSs from `fTarget`/`alphaReqTarget` to
  engine-derived `f`/`alphaReq` so wrong coefficients in either literal
  would surface
- `036`
  red-team batch II.1 replaced hand-typed `dGTarget`, `gMaxTarget`,
  `gSeriesTarget` with derived alternatives (polynomial form, `Limit`
  form, coefficient extraction); added a `disc + 72*delta^2 == 0`
  discriminant check; substantive symbolic kappa-based `F`-`R_target`
  identity inserted in both engines
- `037`
  red-team batch III.1 removed hand-supplied `xiTerm`/`alphaTerm`/
  `sigmaExpected` and `aExpected`/`deltaExpected` literals; the
  Mathematica mirror now reconstructs `xi` and `alpha` from two entries
  of `sigmaWall` and cross-checks the third as a substantive identity,
  and derives `A`/`delta` closed forms via `Together` numerator-and-
  denominator extraction against the Schur closed form
- `038`
  red-team batch III.1 dropped the pre-baked `(cEtaU*cUW + cEtaW*kU)^2
  -> zW kEtaEff kWEff kU^2 (1+rho)^2` substitution rule in `applyDimless`;
  added nine non-tautological sign assertions in both engines (each
  multiplies the symbolically-differentiated derivative by a manifestly
  positive template under the stated transfer-branch assumption and
  verifies the +/-1 residual)
- `039`
  red-team batch III.1 restructured the `.wl` so `deltaSplit`,
  `epsWSplit`, and `dDir` are derived from their own algebra rather
  than postulated, with the SymPy-side postulate moved to the RHS of
  a new `derived matches postulated` check (unlocks engine-disagreement
  detection); replaced the `z1/z0 - (kappa1/kappa0)*R_U` tautology with
  a structurally explicit kappa-rho residual and added flat-U baseline
  substitution checks
- `040`
  red-team batch III.1 added a genuine perturbed-matrix eigenvector
  residual check (both rows reduce to 0 against `M - alpha z z^T`),
  replaced the `series`-vs-`diff` tautology with a two-path cross-check
  for `H_F` (via `F_U` vs via `F_general`), and now derives `alphaReq`
  via `Solve[Det[...] == 0, alpha]` and the eigenvector via `NullSpace`
  instead of postulating the closed form
- `041`
  red-team batch III.1 made the source-tied `n_src` check
  non-tautological by deriving it from the general `n_expected` via
  `q -> t R_U, r -> t, t^2 -> lambda0` substitution in both engines,
  so the assertion now verifies the substitution actually reduces to
  the hand-written form
- `042`
  red-team batch III.1 verified clean as-is: rank-2 selected-mode
  Mathematica mirror is structurally parallel to SymPy but cross-checks
  identities through `FullSimplify`/`Together` canonical paths rather
  than copied algebra; not flagged as transliteration because the
  claims are pure closed-form identities and the agreement is genuine
- `043`
  red-team batch III.1 added five independent algebraic paths in the
  Mathematica script (`Det` forms for `dPhi`/`dPhiZ`, residue-ratio for
  `rPhi`, endpoint limits for `v.D_U.v`, `Series` expansion for
  mismatch); replaced tautological `A_phi^eff` and `M_supp` self-
  comparisons with genuine minimal-overlap and split-vs-minimal ratio
  anchors plus mu-independence derivatives
- `044`
  red-team batch III.1 added an independent `Solve` route for `xiPhys`
  in Mathematica, replaced tautological renames with non-tautological
  coefficient extraction from `branch_eq`, deleted an algebraically
  redundant tracking total-loading assertion, and added a literal slice
  at `Rphi=2` to constrain the bivariate dependence
- `045`
  red-team batch III.1 added a polynomial-extraction route
  (`coupling_density` -> `coeff(...)` -> `g_X_ext`) with four
  `extracted - reference` firewall assertions and an enumerated
  `channels` list giving `M_tr_channel_sum`; the `.wl` self-comparison
  on `mTrReq` was replaced with `Solve[collapsedNum == 0, mTrSym]`, and
  the branch numerator is now derived via `Series[..., {rPhi, rU, 0}] //
  Normal` rather than direct substitution
- `046`
  red-team batch III.1 removed hand-typed `pR`/`p1`/`p2`/`*Expected`
  literals; the `.wl` now uses direct `Together[D[...]]`, `Reduce[
  ForAll[...]]` sign claims, and `PolynomialQuotientRemainder` factor
  checks; both engines gained non-tautological boundary and three-point
  sign-sample assertions operating on `G_tr - G_flat` / `F_flat - F_tr`
- `047`
  red-team batch III.1 closed the `rho_0 - chi_0` and `sigma_0 - chi_0`
  tautologies (the `lamW/lamW` and `lamphi/lamphi` factors had
  cancelled by construction) and rewrote `mSupp`/`sEnhance` in the
  Mathematica mirror through independent algebraic routes; added a
  `PASS: S from ratio agrees with closed-form S` cross-engine identity
- `048`
  red-team batch III.1 verified clean as-is: independently `Solve`s for
  `zeta_req` and adds two limit-coefficient checks (softening, pole)
  absent from the SymPy script; not a transliteration
- `049`
  red-team batch III.2 deleted the `uniformDnOverlap` helper; the
  Mathematica overlap derivation now goes through
  `Integrate[chiN, {s, 0, l}]` (integer assumption) with `i0` obtained by
  `overlapFormula /. n -> 0`. Tautological `k_n` definitional check replaced
  with a non-trivial Neumann-boundary residual `cos(k_n L) == 0`
- `050`
  red-team batch III.2 replaced the `(2n+1)^2 / (2n+1)^2` impossibility-bound
  cancellation with a genuine admissibility-numerator residual; added a
  derivative-sign check on `zeta_n^(twin)(x)`; added the missing factored
  form for `S_n^(max) - S_n^(twin)`; introduced an explicit
  `ConditionalExpression` strip in `expectZero` so post-`Solve` zero residuals
  collapse cleanly under `$Assumptions`
- `051`
  red-team batch III.2 (CHECKPOINT) replaced the `M_mix(Z_W^(twin,req))`
  inverse-roundtrip with a forward-map comparison via the Stage 047/030
  `Z_W = pi^2 (1-eps_eta)(1-eps) M_mix / [8 (1+chi0)^2]` relation; routed
  the Mathematica side through independent `Factor[Together[...]]`
  canonicalization and `Solve[{gTr == 2 mMix, xi > 0}, xi, Reals]` for
  `xi_(2x)`. The `Pi_tr(xi->1-)` infinity check switched from `pi1 =!= Infinity`
  to `1/pi1 == 0` to handle Mathematica's non-deterministic `Limit` output
- `052`
  red-team batch III.2 broke the SymPy/Mathematica `applyDimless`-style
  tautological renames for `zeta_req`, `dZdPi`, `Delta_zeta`, and the
  softening fraction; the Mathematica mirror now derives each via independent
  `Solve`/`Together`/log-derivative routes
- `053`
  red-team batch III.2 replaced the hand-typed `2/pi - 1/2` linear-coefficient
  literals with `series_small.coeff(alpha, 1)` (SymPy) and
  `Coefficient[seriesSmall, alpha, 1]` (Mathematica), so the small-alpha
  coefficient assertion now depends on the integrated `Omega_alpha`
- `054`
  red-team batch III.2 removed two Mathematica hardcoded literals:
  `bExpr = a Tan[k ell]` is now obtained via `Solve` of the Neumann condition,
  and the `x floor` is obtained by inverting `aKMax == zetaReq` via `Solve`.
  SymPy was already clean
- `055`
  red-team batch III.2 re-anchored the `KX/KW equivalence` check to
  `(1/AK).subs(y, 0).subs(x, x_floor)` rather than the hand-typed `1 - x/4`,
  so the assertion now depends on `A_K` itself
- `056`
  red-team batch III.2 verified clean: SymPy and Mathematica derive
  `Omega_Pe`, the twin/finite-throat limits, the covariance identity, the
  small-Pe coefficient, and the large-Pe `-pi^3/8` correction through
  genuinely different mechanisms (`Series.removeO` vs `Limit[pe^2(Omega-pi/2)]`).
  Structural similarity is high but not transliteration
- `057`
  red-team batch III.2 replaced the `y_req identity` self-subtraction with a
  round-trip substitution of `y_req_sq` into the defining equation
  `zeta_req = Omega^2(kappa + pi^2/4)/(kappa + y^2)`; both engines now exercise
  this non-tautologically
- `058`
  red-team batch III.2 re-derived `fc`/`fs`/`delta` through independent
  `Integrate[]` calls in Mathematica (no SymPy antiderivative ansatz import);
  added bracket-gap closed-form + positivity sweep + `Delta_inf as Pe -> oo`
  limit checks across both engines; the constant-term analyticity identity
  was augmented with a genuine non-vanishing `Pe^1` coefficient assertion
- `059`
  red-team batch III.2 (previously: constructive `FindRoot` saturation route)
  swapped the Mathematica linear-coefficient path from `Series`/`Coefficient`
  to `Limit[D[Omega, pe], pe -> 0]`; restructured the circular saturation test
  to use an independent `zeta_req_probe = 2/5` and recover `Pe_star` via a
  fresh `Solve`. Substantive `(4-pi)/pi` claim preserved
- `060`
  red-team batch III.2 replaced the failing `sp.solve` with the explicit
  `Csol = a/(exp(a*L) - 1)` closed form plus Jacobian-aware rescaling
  assertions; swapped the tautological `Pe identification` for a
  `Solve[gamma]`-derived rate check; replaced the Onsager dissipation
  cancellation with a genuine positivity check (`sp.ask` / `Reduce[ForAll[...]]`)
  in both engines; added a `K_X = 0` support-BVP solve that confirms
  `Delta = Lambda L^2 sigma_0 / (2 T_X)` in the `K_m -> infty` limit.
  `material_change: true` was flagged by the verifier — second-pass should
  spot-check downstream Xi_micro consumers
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
