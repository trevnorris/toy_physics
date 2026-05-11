# Moving-Throat PDE — Numerical Stress Plan

## 2026-04-21 Status Update

- `045-048`: hardened in both CAS layers. The numerical check now evaluates
  `R_target(zeta)` at distinct `zeta` values instead of comparing a shared
  scalar placeholder to itself, and the support-load branch checks are now
  phrased as real threshold / pre-threshold assertions rather than a misleading
  signed residual print.
- `106`: hardened in both CAS layers. The old reciprocal-then-multiply
  tautology is gone. The stress layer now checks the exact canonical
  coefficient against the point-particle truncation, verifies the
  `O((a/r)^4)` remainder, and uses symmetric `chi_Q` perturbations to confirm
  the expected monotone and convex branch response.
- `154`: hardened in both CAS layers. The tautological per-profile
  `R(g)`-shift identity was removed; the harness now checks the loaded
  compensation moment `g_*` against its exact closed form and keeps the
  genuinely nonlinear fixed-point / transport checks.
- `155-156` plain harness: retired. The exploratory plain harness was weaker
  than the fixed-point harness and had drifted out of maintenance. Both plain
  entry points now delegate to the authoritative fixed-point stress harness so
  the stale logic cannot silently pass as an independent check. The old
  `stage155_156_config.json` file has been removed so the repo no longer keeps
  a dead scan grid with the known root embedded in it.

## Purpose

The symbolic audits are now strong in two CAS systems. The remaining numerical
hardening goal is different:

- test branch assumptions near their boundaries,
- probe degenerate or nearly singular regimes,
- compare exact formulas against direct residual substitution,
- and make sure the same branch is being selected numerically that the notes
  describe analytically.

This plan defines the first adversarial numerical pass.

## Priority Order

1. `045-048`
2. `106`
3. `125`
4. `142-157`
5. `185-187`
6. `003-004`, `058`, `170`

That order follows risk, not chronology.

## Common Test Rules

- Every sample point must record the assumptions it is intended to satisfy.
- Every sample point must be tagged as either:
  `interior`, `near-boundary`, `near-degenerate`, or `stress`.
- Every comparison should save:
  input parameters,
  exact/reference formulas used,
  numerical residuals,
  and a pass/fail tolerance.
- SymPy/Python and Mathematica should use at least one shared sample table for
  the same stages so branch disagreements are visible immediately.

## Test Families

## 1. Coherent support chain: Stages `045-048`

### Purpose

Stress the constructive coherent branch and the support-compensation theorem.

### Sample regions

- Interior coherent branch:
  `chi_0 in {0.25, 0.5, 1.0}`, `delta_U in {0.25, 1.0, 4.0}`
- Near tracking-flat boundary:
  `delta_U -> 0+`
- Near weak branch:
  `chi_0 -> 0+`
- Near support blow-up:
  `zeta -> (1/eps)^-`
- Near softening:
  `xi -> 1^-`

### Required checks

- Stage `045`: direct numerical equality `rho0 = sigma0` and
  `R_tr = R_U = R_phi`.
- Stage `046`: sign of `dG_tr/dR`, `dF_tr/dR`, and residual ordering
  `E_tr > E_flat`.
- Stage `047`: independence of `R_target` from `zeta` by direct numerical
  variation.
- Stage `048`: verify `zeta_req < zeta_crit < 1/eps` and that the target is
  reached before softening.

### Failure modes to watch

- hidden loss of monotonicity near `R -> 0`,
- numerical branch switching near `xi -> 1`,
- apparent compensation failure caused only by evaluating outside `0 < eps < 1`.

## 2. Canonical outgoing closure: Stage `106`

### Purpose

Check that the canonical branch is numerically distinct from nearby
non-canonical choices and that the GR-target coefficients are not an artifact
of an overconstrained symbolic normalization.

### Sample regions

- Point-particle interior:
  small `a/r`, e.g. `1e-3`, `1e-4`
- Transitional regime:
  moderate `a/r`, e.g. `0.05`, `0.1`
- Canonical branch perturbations:
  `chi_Q = 1 + delta_chi`, `mhat_0 = 1 + delta_m`

### Required checks

- Compare canonical `N_Q = 1` evaluation with nearby perturbed branch choices.
- Measure stability of the extracted `GammaBar_5` under decreasing `a/r`.
- Verify that direct numerical evaluation of the reduced closure approaches the
  symbolic target coefficients.

### Failure modes to watch

- accidental numerical insensitivity masking a real branch choice,
- apparent agreement only at symbolic normalization points.

## 3. Positive-source mouth theorem: Stage `125`

### Purpose

Stress the positive-source restriction against near-extremal and sharply
localized positive sources.

### Sample regions

- Smooth broad positive sources.
- Sharply localized positive sources near the endpoint.
- Sources with cosine moment near `0`, near `1`, and near the lower branch
  value `g_-^F1`.

### Required checks

- Numerically confirm `g[sigma] in [0,1]`.
- Verify exclusion of `g_+^F1` under positive normalization.
- Test robustness of the lower-branch admissibility under extreme but still
  positive localizations.

### Failure modes to watch

- numerical quadrature artifacts near endpoint-localized sources,
- accidental drift outside positivity due to basis truncation.

## 4. Mouth/core closure and co-evolution: Stages `142-157`

### Purpose

Stress the explicit Family-1 branch near its singular surface and along the
co-evolving fixed-point continuation.

### Sample regions

- Canonical neighborhood:
  `Pi` near `Pi_*`
- Large-bias regime:
  `Pi >> 1`
- Near singular denominator:
  `1 - R_q S_q -> 0+`
- Frozen-traction slice:
  `Sigma0 = Sigma0^*`
- Renormalized compensation bracket around the Stage-`156` root

### Required checks

- Stage `142`: direct residual of the self-consistent mouth/core equations.
- Stage `143`: confirm `g_Pi < 1` for finite `Pi` and singular growth of
  `Tmhat`.
- Stage `144`: numerical exclusion of upper/equal-normalized branches inside the
  explicit family.
- Stages `146-147`: compare first-order predictions against finite but small
  `eps` solves.
- Stages `150-152`: compare covariance-based first-order corrections with direct
  numerical profile integration.
- Stages `155-156`: solve the fixed-point map from multiple seeds and record
  whether all interior seeds converge to the same branch on the analyzed window.

### Failure modes to watch

- hidden branch hopping near `1 - R_q S_q = 0`,
- first-order formulas breaking earlier than expected,
- apparent uniqueness in `155-156` caused only by a narrow seed choice.

## 5. Microscopic invariant/quotient chain: Stages `185-187`

### Purpose

Stress whether the invariant coordinates remain numerically well-conditioned
near the boundary of the positive coherent sector.

### Sample regions

- Interior positive coherent states.
- Weak coherent edge:
  `chi0,* -> 0+`
- Weak split edge:
  `deltaU,* -> 0+`
- Small monomial drift samples tangent to `G_*`.
- Small monomial drift samples transverse to `G_*`.

### Required checks

- Stage `185`: compare direct observable defect with monomial-drift residual.
- Stage `186`: numerically verify that tangent directions in `T_id G_*` lie in
  `ker(M_*)`.
- Stage `187`: pick pairs of positive microscopic states with identical
  invariant triples and verify that the solved finite log-ratio laws map one to
  the other.

### Failure modes to watch

- rank loss of the selected minor near `chi0,* -> 0+`,
- numerically ill-conditioned monomial comparisons due to logarithmic scaling,
- false positives from comparing nearly equal large products instead of logs.

## 6. Foundational spot checks: Stages `003-004`, `058`, `170`

### Purpose

These are not the highest branch-risk stages, but they are structural choke
points and deserve a small adversarial pass.

### Required checks

- Stage `003`: sample separated poles with decreasing gap
  `varpi^2 - Omega_eta^2 -> 0+` and verify perturbative drift only before the
  resonance boundary.
- Stage `004`: compare the outgoing odd term with direct numerical evaluation of
  the reduced transfer factor on small `omega`.
- Stage `058`: solve `Pe = Xi Delta(Pe; kappa, eta)` numerically for small,
  moderate, and larger `Xi`, then compare with the bracket
  `Xi Delta_0 <= Pe_* <= Xi Delta_inf`.
- Stage `170`: numerically reconstruct `(K_A, G_A)` from grouped lane data and
  verify the hidden-even consistency law.

## Artifact Layout

- SymPy/Python:
  `scripts/numerical/`
- Mathematica:
  `research/pde_ledger/mathematica/numerical/`
- Shared recorded outputs:
  `notes/moving_throat/review/numerical/`

Suggested file naming:

- `stage045_048_tracking_stress.py`
- `stage106_canonical_outgoing_stress.py`
- `stage125_positive_source_stress.py`
- `stage142_140_mouth_core_stress.py`
- `stage185_187_orbit_stress.py`

and corresponding `.wl` files under the Mathematica numerical directory.

## Acceptance Criteria

- Each priority block has at least one shared sample table used by both SymPy
  and Mathematica.
- Every test reports whether it is checking:
  branch selection,
  residual substitution,
  asymptotic accuracy,
  or uniqueness/conditioning.
- Any disagreement is triaged immediately into:
  note-statement issue,
  script issue,
  branch/domain issue,
  or numerical-conditioning issue.

## Recommended First Implementation Slice

1. Build the shared sample-table format.
2. Implement the `045-048` stress harness in Python and Mathematica.
3. Implement the `155-156` multi-seed fixed-point stress check.
4. Implement the `185-187` invariant/quotient conditioning check.

That gives the best coverage of the remaining conceptual-risk surface with the
least additional symbolic work.
