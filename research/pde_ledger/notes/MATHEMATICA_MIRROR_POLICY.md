# Mathematica Mirror Policy

This document defines how the PDE ledger should talk about Mathematica
coverage.

Snapshot date: `2026-05-26` (batch III.3 v2 close)

III.3 v2 re-audit (2026-05-26) caught *new* `mathematica_transliteration` patterns on 4 of 7 dirty stages (062 polytrope/parent-action chain port, 064 variational-closure chain, 068 P_res/Solve choreography, 070 Hw→Nphiphi→Xi chain). Down sharply from III.3 v1's 9/12 transliteration share — the v2 audit's paper-grounded reading caught remaining gaps. The III.3 standout was a **Mathematica-Integrate limitation**: Codex iter1 on stage 064 prescribed integrating `hFun[y] * chiSigmaFun[y]^2` (where `chiSigmaFun = gPhi*chiPhi/hFun`) and asserting equality with `gPhi^2 * Integrate[chiPhi^2/hFun]`. Algebraically the integrands are identical (`gPhi^2 chiPhi^2/hFun`), but `FullSimplify[Integrate[gPhi^2 * c^2/h] - gPhi^2 * Integrate[c^2/h]]` does NOT reduce to 0 — Mathematica cannot pull a constant factor outside `Integrate[...]` when the integrand contains unspecified symbolic functions. The exec-mathematica failed with surface form `-gPhi^2*Integrate[c^2/h] + Integrate[gPhi^2*c^2/h]`. **NEW pitfall #9 candidate: verify integrand equality first, NOT integral equality, when symbolic integrand factors are involved.** Orchestrator hot-fix: `FullSimplify[hFun*chiSigmaFun^2 - gPhi^2*chiPhi^2/hFun] == 0` (integrand-level check) then define `thetaGeneral = lambdaGeneral = gPhi^2 * i1Integral` directly. Integrand equality implies integral equality. Promote pitfall #9 to `codex.md` if recurs.

Pitfall #8 (heavy BVP `dsolve`/`DSolve`) **was promoted** from candidate to documented in `codex.md` "Common cross-engine pitfalls" item #1 before III.3 launched, as preemptive defense. No BVP verifications were prescribed in III.3, so pitfall #8 did not recur.

III.2 v2 re-audit (2026-05-26) caught a *new* `mathematica_transliteration` pattern on only 2 of 7 dirty stages (050 derivative/x_max/ceiling SymPy-targets ported byte-for-byte; 057 `aKX/xSub/aKKappa` chain), the lowest transliteration share seen yet in v2. The III.2 standout was a **performance pitfall**: Codex iter2 on stage 058 added a full `sp.dsolve` / `DSolve` symbolic BVP solve + boundary-condition `sp.solve`/`Solve` + `sp.simplify(phi_drop - Delta)` to verify F2 (Green-kernel construction). The sympy version hung at 100% CPU for 7+ hours before the orchestrator killed it; the Mathematica `DSolve` mirror would have had the same pathology. **pitfall #8 (heavy-machinery BVP checks via `dsolve` are not worth the symbolic cost)** — the same Green-function identity is verifiable in seconds via `Delta = integral(K * Sigma)` (kernel-integral identity). Orchestrator hot-fix replaced both engines' dsolve blocks: sympy uses a numerical sweep `Delta = integral(K * Sigma_Pe)` on 4 concrete `(α, η, Pe)` tuples; Mathematica relies on its pre-existing line 84 `delta independent integral matches combination form` check (which compares `Integrate[kernel*sigmaPe]` Green-function side vs the `Ic`/`Is` combination closed form). Also added `Pe == α` singularity guards in the sympy monotonicity/IVT sweeps because `Delta` has a removable 0/0 at `Pe = α` and `subs()` doesn't take the limit. Promoted to `codex.md` 2026-05-26.

III.1 v2 re-audit (2026-05-26) caught *new* `mathematica_transliteration` findings in stage 043 (full Sections 1-5 port) and stage 045 (coefficient-chain port). Stage 043's fix introduced a separate pitfall (NEW pitfall #7 candidate): **primitive-vs-derived substitution**. The directive's F2-Insertion2 numeric-point sign anchor prescribed substituting `sigma_0=0, rho_0=1` symbolically, but the Mathematica expression's primary symbols are the primitive couplings `gB, gU, gS, gR, gW, kU` — `sigma_0` and `rho_0` are derived quantities (`sigma_0 = gU·gS/(kU·gB)`, `rho_0 = gU·gR/(kU·gW)`). Symbolic substitution on derived names doesn't reduce the primitive expression. Iter2 fix: substitute on primitives that realize the derived value — `gS → 0` for sigma_0=0, `gW → gU·gR/kU` for rho_0=1. Defense: when a directive prescribes `subs(<derived>, value)`, lift to the primitive symbols that realize the derived value. Promote to `codex.md` if recurs.

II.1 v2 re-audit (2026-05-26) caught a *new* `mathematica_transliteration` finding in stage 024 (Sections III/V were line-by-line ports despite v1 II.1 having flagged the file as already rewritten). The remediation pass for 024 also fixed a separate **performance pitfall**: heavy `Table[tripleOverlap[...], {i,1,5}, {j,1,5}]` with 6-fold inner sums + FullSimplify hangs >18min when global-symbol context leaks from earlier sections; the fix is `ClearAll[<symbol-list>]` reset at the top of Section IV plus memoization of `i4`/`i6` sphere integrals via `i6[args] := i6[args] = Integrate[...]`. This pattern (heavy FullSimplify after large symbol additions in earlier sections) should be screened pre-emptively in future Mathematica-heavy stages. II.1 v2 also caught transliteration in 031 (full PART I-II port — replaced with `Eigenvalues`/`Eigenvectors`-based independent derivation) and 032 (Stage 15.4-15.5 port — replaced with `Eigensystem`-based independent path).

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
already broke line-by-line correspondence; batch III.3 found it on 9 of 12
(062-065, 067-070, 072) — share back up because the stage cluster 062-068
was particularly mirror-heavy on asymptotic leading-order forms, with
stages 061 and 066 passing transliteration screening on first read; batch
III.4 found it on 7 of 12 (073, 075, 076, 078, 079, 082, 083) — slightly
below III.3's share because the Family-1 numerology cluster (075-084)
distributed its 40 findings across a broader category mix (12
`hardcoded_result`, 14 `tautological_check`, 7 each of `math_translit`
and `insufficient_verification`), but `material_change: false` on every
stage. Treat transliteration as the default expectation, not an
exceptional finding.

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
- `062`
  red-team batch III.3 replaced hand-built `sigma_star` and EOS-monomial
  assertions with a real parent-action Gaussian-elimination `sp.solve` for
  `sigma_star(rho)` plus an `n=5 EOS identity` wrong-exponent probe printing
  a nonzero residual `K*rho^3*(5 - 6*rho)`; Mathematica side now derives
  `sigma_star` via `Series`/`SeriesCoefficient` rather than mirroring the
  SymPy choreography
- `063`
  red-team batch III.3 replaced all 8 SymPy + 8 Mathematica tautologies with
  `sp.solve` / `Reduce` derivations and added a Cauchy-saturation check that
  would fail under plausible bugs (e.g. `N_ss`-vs-`N_pp` swap in `G_max`);
  Mathematica side switched to `Reduce[... && gphiSq>0, ..., Reals]` to
  break the SymPy mirror
- `064`
  red-team batch III.3 instantiated `chi_phi(y)` / `H(y)` profiles and
  derived the `I1`/`I2` integral overlap reductions rather than carrying
  them as `const_subs` literals; one orchestrator-acceptable codex
  deviation (an extra `Npp -> I1*Hw` substitution required for F1's
  independently-verified integral identity to canonicalize) was retained
- `065`
  red-team batch III.3 anchored docstring claims (1)-(3) and (6) with
  concrete-Gaussian shell integrals `J2=0` (by parity), `I1` polynomial
  coefficient expansion, and `gphi`'s `1/ell` scaling from `V_conf`
  differentiation; `K_X`-cancellation identities remain as the surviving
  substantive checks
- `067`
  derives stationarity from the self-dual `C^2(r)=C^2(pi/r)` symmetry equation;
  red-team batch III.3 additionally integrated the transverse norms `2 w_f`
  and `w_g sqrt(pi/2)` explicitly in SymPy (via `.rewrite(sp.cosh)`-workaround
  `Integrate[sech^2]` and `Integrate[Gaussian^2]`) and made two
  "duality implies stationarity" blocks non-tautological while preserving
  the numerical duality/monotonicity checks
- `068`
  red-team batch III.3 derived `Wfail_res := Pe_req/(C2*Delta_inf)` and the
  `Wfail_match` cross-relation via `Solve` from explicit resonance-corrected
  premises rather than postulating them; `M3/M4` cross-relations now carry
  the load-bearing assertions; `material_change: true` flagged because the
  derivation route changed even though the symbolic content of the derived
  expressions matches the prior postulated forms
- `069`
  closes the ordered three-zone regime algebra rather than width-only replay;
  red-team batch III.3 (CHECKPOINT) additionally replaced
  `Cres2`/`Wfail_res`/`delta_fail` definitional identities with a
  parameterized `W_match` generator + monotonicity check (SymPy) and a
  `Cres2Prim` primitive + `Pres = 1/Cres2` derivation + `PresGap` via
  `Solve` (Mathematica); upstream `Pres`/`Wfail_match` carry-forward
  annotated with provenance comment
- `070`
  red-team batch III.3 inlined `1/Hw -> rhoW/(m*cSw^2)` in the Mathematica
  mirror and removed mirrored intermediates `J_1`, `gphi`, `I_1` while
  retaining the 3 `expectZero` calls; assembled forms now built from
  primitives and closed forms typed directly
- `071`
  red-team batch III.3 replaced the `eta - L/ell` tautology with a pin
  `K_m = pi a^2 hbar^2/(3 m rho_w)` and an `eta`-reconstruction from the
  closed-form `K_m` in both engines
- `072`
  red-team batch III.3 added 4 ratio-limit checks per engine comparing the
  full `Delta_0`/`Delta_inf` closed forms to the shell- and
  compression-dominated leading-order forms; SymPy collapses both shell
  ratios to `1` directly while Mathematica's `DeltaInf` shell ratio surfaces
  as the algebraically-equivalent surd `2/Sqrt[5] + 1/(5 + 2 Sqrt[5])` —
  divergent presentations confirm cross-engine independence
- `073`
  red-team batch III.4 fixed a Mathematica precedence bug (`eta /. (len/ell) -> 37 - 37`
  parsed as `(len/ell) -> 0` because `Rule` has lower precedence than `Plus`); both
  engines now build `eta`/`Lambda_ell` from symbolic `K_m` / `L/ell` with a symbolic
  `Lambda_ell - L/ell` identity check before numerical specialization
- `074`
  red-team batch III.4 replaced the `chi_lock = Lambda_ell/2` tautology with the
  physical substitution chain `chi_def = m_psi*c_s*L/hbar -> subs(c_s) -> subs(L/ell)`
  in SymPy; Mathematica side was already the more substantive engine here, deriving
  `chi_s` from `m c_s L / hbar` with the healing-length substitution rather than
  mirroring SymPy
- `075`
  red-team batch III.4 replaced the `100*Theta_w == 100*Theta_w` round-trip with a
  free-symbol `Delta_0`/`Delta_inf` algebraic identity check (alpha_sym/eta_sym must
  be proved by `simplify`/`FullSimplify` rather than being self-consistent literals);
  cross-engine independence achieved via the symbolic-route F2 leg
- `076`
  red-team batch III.4 derived `U` from `Integrate[P/rho^2]` rather than typing
  `K*rho^5/4`; routed `mu_star` and `c_sw` through `sp.solve`/`Solve` rather than
  postulating; paired the n=5 enthalpy identity with a non-tautological n=3
  fail-check (the assertion holds for `n=5` and fails for `n=3` because `P` itself is
  index-parameterized as `K*rho^n_poly`); renamed `thetaHealTarget` to
  `thetaHealReduced` to break line-by-line `.wl`/`.py` correspondence
- `077`
  red-team batch III.4 added the symbolic `1 - alpha_r * S(xi_*)^2 = 0` identity
  check in both engines (proves `xi_* = atanh(2/sqrt(alpha_r) - 1)` is the cut point
  rather than postulating it); removed the tautological `xi_* numeric check` from
  the `.wl`; added SymPy `expect_close` per-value Jensen-floor tolerance checks
- `078`
  red-team batch III.4 replaced literal-decimal `Theta_chi_coeff`/`Theta_J_coeff`
  copies with symbolic `Sinh`/`Cosh` closed forms in Mathematica and high-precision
  `ToExpression["...`40"]` loads; `expectApprox` targets now computed in-script;
  three branch-verdict ordering inequalities added; one acceptable codex deviation
  (removed a spurious `100` factor in the directive)
- `079`
  red-team batch III.4 replaced literal decimal `zeta0`/`zetaInf` copies with
  `Limit[aF1*omega^2, pe -> 0/Infinity]` (high-precision `expectApprox` fallback per
  directive); F3 added a `D[omega, pe]` symbolic-derivative slope check returning
  `(4-Pi)/(2*Pi)` independently from SymPy's `series`
- `080`
  red-team batch III.4 added four `zetaTarget*` independent-path numeric
  computations in Mathematica (replacing literal SymPy-output copies) and four
  ordering inequalities (chi-pair, J-pair, J<=chi suff, J<=chi fail); SymPy added
  four `zeta numeric check` lines exercising Stage-61 Pe constants
- `081`
  red-team batch III.4 derived `qq = piOfZeta / cMix` via `Solve[zeta == zetaExpr,
  piTr]` rather than postulating; orchestrator retrofitted the standard
  `ConditionalExpression[e_, _] :> e` strip (Solve introduced the wrapper, which
  broke the `=== 0` test and downstream `zeta -> {0, 1, ...}` substitutions);
  five magic-number self-checks replaced with residuals against `1 + zeta_*`
  functional form
- `082`
  red-team batch III.4 routed Mathematica `zetaReq` through `Solve` of `qMap`
  (with `ConditionalExpression` strip) rather than mirroring the SymPy formula;
  two `Xi_F1` tautological self-checks demoted to `print`; two new derivative
  assertions (`dR_quad/dzeta_phys + 1 = 0` and `dR_quad/dPi_tr - dzeta_req/dPi_tr`)
  added in both engines
- `083`
  red-team batch III.4 added `delta0Residual`/`deltaInfResidual`/`omegaResidual`
  cross-engine identity checks, `y_F1` defining-equation residual gate, and a
  `dzetaDpe` monotonicity sign-check; SOURCE-ANCHOR comment blocks for
  `Theta_chi_coeff`/`Theta_J_coeff`/`136900` literals; one benign deviation
  (`nsolve(prec=80)` bump for numerical stability)
- `084`
  red-team batch III.4 replaced two `(37^2 - 1369)`-style tautologies with the
  cross-route consistency check `(xiF1FromUpsilon /. upsilonW -> 100*thetaW) -
  xiF1FromTheta`; two hardcoded float-difference checks replaced with four
  `expectZero[If[TrueQ[...], 0, 1]]` ordering inequalities; added
  `FindRoot`/`Limit[zetaPhys, Pe -> Infinity]` block returning `2.467529229456...`
  matching `zetaMaxF1` to ~14 digits
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
