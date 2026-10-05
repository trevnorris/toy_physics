# S10 baseline formalization result — 2026-09-11

The arbitrary-dimensional baseline proof passes Lean. It derives the real
plane-wave operator from the supplied antisymmetric-derivative action, proves
its complete amplitude classification, and connects this classification to
stationarity against arbitrary smooth compact variations. Explicit equality
theorems recover the existing S9 pilot at D=3.

## Supplied action and checked conclusion

The input is the baseline action in
[S10_SHARED_PHYSICS.md](../../directives/S10_SHARED_PHYSICS.md): a D-component
in-plane displacement on flat spacetime R^(D+1), with constant real coefficients,

```text
A_ij(J) = J_(i+1),j - J_(j+1),i,
S_curl(J) = (1/2) Σ_i Σ_j A_ij(J)²,
L(J) = (rho/2) Σ_j J_time,j² - (mu/2) S_curl(J).
```

Here row zero of J is time, and the other rows are spatial derivatives. Both the
factor 1/2 in the double sum and the factor mu/2 in the density are retained.
The reduced modal expression is a theorem proved from this input, not an input
to the action definition.

For real `u(t,x) = a cos(k·x - omega t)`, Lean proves

```text
EL[u] = cos(k·x - omega t) M(omega,k)a,
M(omega,k)a = rho omega² a - mu (|k|² a - k(k·a)).
```

For every positive integer D, `rho > 0`, `mu > 0`, and `k ≠ 0`, the combined
theorem `s10_variational_certificate` supplies a positive cone frequency with

```text
omega² = (mu/rho)|k|².
```

At that frequency, the entire stationary plane-wave amplitude space is `k·a=0`
and has dimension D-1. At zero frequency, it is `span{k}` and has dimension one.
Stationarity here quantifies over all smooth compact vector variations, including
variations outside the plane-wave ansatz. A separate equivalence covers every
nonzero real frequency and nonzero amplitude: stationarity holds exactly when
the amplitude is transverse and the frequency satisfies the cone relation.

These statements hold for every nonzero wavevector, including coordinate axes;
there is no additional generic-direction or exceptional-stratum hypothesis for
this baseline action. `coneValue_scaling` also proves homogeneity under real
rescaling of the wavevector.

## What this strengthens and what the edge case says

The existing [S10 record](../../steps/S10_two_transverse_photons.md) reports a
finite sweep at D=2,3,4,5, with a genericity qualification on the recorded
computational evidence. The new theorem proves the baseline amplitude count for
all D and all nonzero wavevectors under the stated premises. In particular, it
proves positive cone value at D=5 without relying on a comparator sign verdict.
It does not repair or certify the CAS comparator, exports, or their provenance.

| Spatial dimension | Transverse amplitude dimension on the cone | Longitudinal amplitude dimension at zero frequency |
|---|---:|---:|
| 1 | 0 | 1 |
| 2 | 1 | 1 |
| 3 | 2 | 1 |
| 4 | 3 | 1 |
| 5 | 4 | 1 |
| D ≥ 2 | D-1 | 1 |

`exists_propagating_mode` proves that the positive-frequency amplitude space
contains a nonzero vector for every D ≥ 2. In D=1, the antisymmetric stiffness
vanishes identically, and `dimension_one_no_propagating_mode` excludes every
nonzero amplitude at every nonzero real frequency. Thus the positive value of
the cone formula by itself does not establish a propagating mode in D=1. This
edge case is consistent with the D-1 count and lies outside the old D=2..5 sweep.

## Integrated variation and phase averaging

The analytic setting is smooth backgrounds on all of R^(D+1), with smooth,
compactly supported variations and Lebesgue measure. Lean proves integrability
of the pointwise action change

```text
ΔS[u,h](s) = ∫ [L(∂(u+s h)) - L(∂u)] d^(D+1)x,
```

then its derivative, integration by parts, and stationarity iff the pointwise
Euler–Lagrange PDE. As in S9, no finite total action is assumed for an infinite
plane wave. If the background density is integrable, a separate theorem proves
the derivative of the total action itself. The fundamental lemma for compact
smooth tests and continuity supply the converse from stationarity to the PDE.

The phase-average route uses the actual integral from phi=0 to 2π. It proves

```text
<L> = (1/2) modalAction,
d/ds <L(a+s b)> at s=0 = (1/2) dot(M a, b).
```

The factor is explicit, and stationarity of the phase average is proved
equivalent to modal stationarity. No time period divided by omega is used, so
this averaging identity also holds at zero frequency. The phase average and
local PDE share the supplied action; agreement is not independent evidence for
its physical origin.

## Source map and S9 agreement

| File | Role |
|---|---|
| [Action.lean](S10Pilot/Action.lean) | Supplied double-sum stiffness; density variation; computed contraction; modal action and operator. |
| [PlaneWave.lean](S10Pilot/PlaneWave.lean) | Actual coordinate derivatives, momentum obtained from the density, and local PDE evaluation. |
| [Spectrum.lean](S10Pilot/Spectrum.lean) | Transverse and longitudinal spaces, rank-nullity count, all-frequency classification, positivity, and scaling. |
| [Certificate.lean](S10Pilot/Certificate.lean) | Combined local-PDE census. |
| [Analytic.lean](S10Pilot/Analytic.lean) | Smoothness, support, integrability, and coordinate integration by parts. |
| [Variation.lean](S10Pilot/Variation.lean) | Integrated action change and stationarity/PDE equivalence. |
| [FiniteAction.lean](S10Pilot/FiniteAction.lean) | Compatibility with finite total action. |
| [VariationalCertificate.lean](S10Pilot/VariationalCertificate.lean) | Combined census through integrated stationarity. |
| [PhaseAverage.lean](S10Pilot/PhaseAverage.lean) | Actual phase integral, its derivative, and equivalent stationarity. |
| [Specialization.lean](S10Pilot/Specialization.lean) | Exact D=3 equalities with S9. |
| [EdgeCases.lean](S10Pilot/EdgeCases.lean) | Existence for D ≥ 2, absence in D=1, and concrete D=5 positive/negative witnesses. |

The specialization proves equality of the stiffness, Lagrangian, modal action,
operator, cone value, transverse and longitudinal spaces, real plane wave,
local PDE, action change, and stationarity predicate with their S9 counterparts.
The S9 proof modules retain their prior contents. Both libraries share the
existing pinned Lean/Mathlib/Physlib environment; no new installation was needed.
The new derivation uses Mathlib directly, as did the S9 pilot.

## Verification

`lake build` checks both pilots, with warnings treated as errors for their
modules. The root [S10Pilot.lean](S10Pilot.lean) audits 34 declarations; their
dependency chains use only the standard logical axioms `propext`,
`Classical.choice`, and `Quot.sound` (or a subset). No audited chain contains
`sorryAx` or a custom physics axiom. The 21 original S9 audit declarations also
pass. See [VERIFICATION.txt](VERIFICATION.txt) for the captured audit and
source hashes, and [README.md](README.md) for the dependency pins and commands.

The concrete D=5 witness uses `rho=2`, `mu=8`, `omega=2`, and a wavevector along
the fifth axis. Its amplitude along the first axis is stationary against all
compact test fields. At the same frequency and wavevector, a longitudinal
amplitude has a compact test variation with nonzero first variation.

A deliberate normalization mutation in an ignored scratch copy replaced the
stiffness double sum's factor 1/2 by 1 while retaining the original claimed
formulas and proofs. Lean exited 1, rejecting the quadratic increment identity
at line 83 and the reduced modal-action identity at line 140. The canonical
source retained the prescribed normalization. This demonstrates sensitivity to
that particular input error; it is not a claim to detect every possible mismatch
between intended physics and a formal statement.

## Remaining boundary and next target

This report covers the S10 baseline. The subsequent full-gradient and
divergence-only controls are recorded in [CONTROLS_RESULT.md](CONTROLS_RESULT.md).
Anisotropic inertia and its special directions, sign reversal, and the remaining
controls have not been formalized. The next bounded target is anisotropic
inertia, which can split frequencies and change exact transversality.

The physical selection of D=3, supplied field content, curl-only action,
positive material constants, and unstrained nondissipative background remain
premises. There is no derivation of microscopic elasticity, out-of-plane
separation, boundary/interface laws, weak solutions, general Fourier/PDE
completeness, or nonlinear/confinement physics. Dimensional-unit analysis and
the GNLS obstruction remain outstanding from S9. No claim that all of S9 or
S10 is certified follows from this increment. No discrepancy with the selected
baseline results was found.
