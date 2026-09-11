# S9 pilot result — 2026-09-11

Lean now checks the chain from the supplied quadratic action, through integrated
stationarity and the local Euler–Lagrange PDE, to the three-dimensional
plane-wave mode census. The integrated-action/local-equation gap from the first
pilot is closed for smooth backgrounds on R^4 and smooth, compactly supported
variations. This reproduces the selected S9 results within that setting; it does
not establish every claim in the S9 record.

The subsequent arbitrary-dimensional S10 extension and its checked D=3 agreement
with this pilot are reported separately in [the S10 result](../s10/RESULT.md).

## The checked statement

Supply constant real coefficients `rho > 0`, `mu > 0`, a spatial wavevector
`k ≠ 0` in three dimensions, and the S9 density

```text
L(J) = rho/2 |J_time|² - mu/2 |curl(J)|².
```

For a smooth background `u` and any smooth, compactly supported vector field `h`,
define its action change using ordinary spacetime Lebesgue measure:

```text
ΔS[u,h](s) = ∫ [L(∂(u+s h)) - L(∂u)] d⁴x.
```

Lean proves that this integrand is integrable for every real `s`, computes the
derivative, and proves

```text
d/ds ΔS[u,h](0) = ∫ EL[u]·h d⁴x,
(∀ compact smooth h, d/ds ΔS[u,h](0) = 0) ↔ (∀ x, EL[u](x) = 0).
```

The action change integrates the pointwise density difference. It therefore
accommodates a nondecaying plane wave without assuming finite total action or
subtracting divergent integrals. When the background density is integrable,
`finiteAction_hasDerivAt` also proves the derivative formula for the total action
`S[u+s h]` itself. Positivity of the coefficients is unnecessary for these
variational identities; it enters the mode census below.

For the real field `u(t,x) = a cos(k·x - omega t)`, the resulting expression is

```text
EL[u] = cos(k·x - omega t) M(omega,k)a,
M(omega,k)a = rho omega² a - mu (|k|² a - k(k·a)).
```

The combined theorem `s9_variational_certificate` proves that a positive
frequency exists with `omega² = (mu/rho)|k|²`. At that frequency, stationarity
against **all** smooth compact test fields is equivalent to `k·a = 0`, and this
amplitude space has dimension two. At zero frequency, the stationary amplitude
space is `span{k}`, of dimension one.

`propagating_variational_mode_iff` establishes completeness within the ansatz:
every stationary plane wave with nonzero real frequency and nonzero amplitude
must be transverse and satisfy that dispersion relation. `coneValue_scaling`
proves squared-frequency homogeneity under `k ↦ s k`. The longitudinal direction
is present with zero restoring stiffness; it has not been removed as a degree of
freedom.

## How the proof connects to the action

| File | What its proofs establish |
|---|---|
| [Action.lean](S9Pilot/Action.lean) | Actual derivatives of the supplied density and its modal restriction; modal stationarity iff `M a = 0`. |
| [PlaneWave.lean](S9Pilot/PlaneWave.lean) | Coordinate derivatives of the trigonometric fields; the pointwise density; momentum obtained by differentiating `L`; the local Euler–Lagrange expression and its plane-wave evaluation. |
| [Analytic.lean](S9Pilot/Analytic.lean) | Smoothness, compact support of test derivatives, integrability, and coordinate integration by parts. |
| [Variation.lean](S9Pilot/Variation.lean) | Integrability and differentiation of the action change; integrated first variation; stationarity iff the pointwise local PDE. |
| [FiniteAction.lean](S9Pilot/FiniteAction.lean) | Agreement with differences of total actions and the total-action derivative when the background density is integrable. |
| [Spectrum.lean](S9Pilot/Spectrum.lean) | Transverse and longitudinal classifications, nonzero-frequency completeness, positivity, homogeneity, and dimensions two and one. |
| [Certificate.lean](S9Pilot/Certificate.lean) | The local-PDE census and concrete solution/nonsolution witnesses from the first pilot. |
| [VariationalCertificate.lean](S9Pilot/VariationalCertificate.lean) | The census through integrated stationarity, including positive and negative variational witnesses. |

The position-space route defines momentum using a derivative of `lagrangian`,
then differentiates it in spacetime. Equality to the proposed modal matrix is
proved. The integrated proof uses an exact quadratic expansion to differentiate
the actual integral. All finite sums and integral operations have integrability
proofs. Integration by parts follows from compact support. The fundamental lemma
for smooth compact test functions gives almost-everywhere vanishing; continuity
and full support of Lebesgue measure give pointwise vanishing. Neither boundary
cancellation nor the field equation is introduced as an axiom.

These proofs use Mathlib's calculus, measure theory, and linear algebra. Physlib
is installed and its existing wave theorem remains an installation check; the
custom S9 derivation does not reuse Physlib's physical field definitions. No CAS
output or numerical fingerprint is a proof premise.

## Verification

`lake build S9Pilot` passed with warnings treated as errors for the pilot modules.
All 21 audited declarations, including both combined certificates and the
integrated-action/PDE equivalence, report exactly

```text
[propext, Classical.choice, Quot.sound]
```

Their dependency chains contain no `sorryAx` or custom physics axiom. See
[VERIFICATION.txt](VERIFICATION.txt) for captured results and source hashes, and
[README.md](README.md) for pinned dependency versions.

For `rho = mu = omega = 1`, `k = (0,0,1)`, Lean proves that amplitude `(1,0,0)` is
stationary against every admissible test variation. For amplitude `(0,0,1)`, it
proves the existence of a smooth compact variation with nonzero first variation.

Two deliberate changes were checked in ignored scratch copies while preserving
the original claimed formulas and proofs:

- Changing only the density's shear sign from minus to plus made Lean reject the
  first-variation identity (exit 1, `action_sign_mutation.lean:52`).
- Reversing only the density difference in `relativeAction` made Lean reject the
  integrated quadratic-expansion proof (exit 1,
  `relative_action_sign_mutation.lean:74`).

The canonical source retained its original signs. These checks show that the
existing proofs fail under those changes; they do not establish that every
possible mistranscription would be detected.

## Remaining boundary

The formalized setting is D=3, flat whole spacetime R^4, constant material
coefficients, smooth backgrounds, and smooth compact variations. Smoothness is
a sufficient regularity assumption, not a claimed minimal one. Finite-domain
boundary conditions, interface jump laws, curved domains, and weak solutions
have not been formalized here.

The original S9 library does not include arbitrary-dimensional mode counts or
dimensional units; the later S10 libraries supply those extensions, with exact
D=3 baseline agreement. Remaining exclusions across this checkpoint include
completeness of general PDE solutions through Fourier superposition; the GNLS
no-transverse-mode argument; confinement and nonlinear physics; and derivation
of the curl-only action or its material constants from a microscopic medium.
The supplied density and its physical interpretation remain premises. The CAS
scripts, exports, and ledger prose have not themselves been certified or modified
by this pilot. No discrepancy with the selected S9 results was found. The
[coverage map](COVERAGE.md) distinguishes the completed pilot from the remaining
S9 ledger obligations.
