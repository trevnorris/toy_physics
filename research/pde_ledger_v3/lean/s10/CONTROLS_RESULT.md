# S10 stiffness-form controls — 2026-09-11

The full-gradient and divergence-only controls are formalized in arbitrary
spatial dimension, from their supplied densities through integrated stationarity
and the local PDE to complete real plane-wave amplitude classifications. The
result confirms the S10 record's distinction: the transverse count is not unique
to curl-only stiffness, while the longitudinal sector distinguishes the forms.

## Inputs and the proved comparison

The two action inputs are the equations in
[S10_SHARED_PHYSICS.md](../../directives/S10_SHARED_PHYSICS.md), section 7:

```text
L_full(J) = (rho/2) Σ_j J_time,j² - (mu/2) Σ_i Σ_j J_spatial_i,j²,
L_div(J)  = (rho/2) Σ_j J_time,j² - (mu/2) (Σ_i J_spatial_i,i)².
```

Each density is defined directly. The stiffness forms are not introduced by
editing an already-derived modal matrix, and full-gradient stiffness is not
identified pointwise with curl stiffness plus divergence stiffness. One concrete
D=2 jet in the proof has baseline density 0, full-gradient density -1, and
divergence-only density -2 at rho=mu=1, so the inputs are demonstrably distinct.

Assume constant rho>0 and mu>0, a nonzero real spatial wavevector k, and real
`u(t,x)=a cos(k·x-omega t)`. There is one common positive cone parameter

```text
omega² = (mu/rho)|k|².
```

The certificate compares stationarity against **all** smooth compact vector
variations for the three action choices at that same frequency. The dimensions
below hold for every D≥1 and every nonzero wavevector, including coordinate axes.

| Action | Amplitude space on the cone | Total dimension | Exactly-transverse dimension | Amplitude space at zero frequency |
|---|---|---:|---:|---|
| Curl-only baseline | k-perpendicular space | D-1 | D-1 | span{k}, dimension 1 |
| Full-gradient | All D-component amplitudes | D | D-1 | {0}, dimension 0 |
| Divergence-only | span{k} | 1 | 0 | k-perpendicular space, dimension D-1 |

Exactly-transverse dimension means the dimension of the intersection of the
modal kernel with `k·a=0`. In particular, FULLGRAD's additional propagating
direction does not increase its exactly-transverse count. DIVONLY exchanges the
baseline's propagating and static sectors. Zero-frequency directions are counted
as surviving degrees of freedom, not removed.

The all-frequency theorems also exclude extra nonzero real-frequency branches:
for any nonzero amplitude, FULLGRAD is stationary exactly on this cone; for any
nonzero frequency and nonzero amplitude, DIVONLY is stationary exactly when the
amplitude is longitudinal and the frequency lies on the same cone.

At D=1 the two controls have equal stiffness densities and each supports a
longitudinal propagating amplitude. The baseline has zero transverse dimension
there, as proved in the preceding pilot. This coincidence of the controls in
one dimension is proved explicitly rather than excluded by an assumption.

## What Lean derives

The density derivatives produce the spatial momenta

```text
p_full,ij = -mu J_spatial_i,j,
p_div,ij  = -mu delta_ij Σ_m J_spatial_m,m,
```

and the same temporal momentum `rho J_time,j`. Actual coordinate derivatives of
the real plane wave then yield

```text
EL_full[u] = cos(phase) [(rho omega² - mu|k|²) a],
EL_div[u]  = cos(phase) [rho omega² a - mu k(k·a)].
```

The modal variation independently identifies these operators with the
amplitude derivatives of the supplied densities' restrictions. The phase average
uses an actual integral over phi=0..2π and keeps its factor 1/2: its first
variation is `(1/2) dot(M a,b)`. No time period involving division by omega is
introduced.

For either form, the analytic proof defines

```text
ΔS[u,h](s) = ∫ [L(∂(u+s h)) - L(∂u)] d^(D+1)x.
```

It proves this difference is integrable, differentiates its exact quadratic
expansion, integrates by parts with justified integrability, and proves
stationarity iff the pointwise local Euler–Lagrange expression vanishes.
Compact smooth tests, the fundamental lemma, and continuity provide the converse.
A separate theorem gives the derivative of total action when its background
density is integrable. Nondecaying plane waves therefore require no fictitious
finite total action assumption.

## Proof organization

S10 sources now live in `lean/s10/`; S9 sources remain in `lean/s9/`. The shared
configuration and dependency cache live at `lean/`. Library source directories
are explicit in `lakefile.toml`, and the dependency revisions are unchanged.
See the [shared build instructions](../README.md).

| File | Role |
|---|---|
| [Action.lean](S10Controls/Action.lean) | The two supplied densities; actual variations; reduced modal actions and operators. |
| [PlaneWave.lean](S10Controls/PlaneWave.lean) | Derived momenta, coordinate PDE, and real plane-wave evaluation. |
| [Analytic.lean](S10Controls/Analytic.lean) | Smoothness of the new momenta and Euler–Lagrange expressions. |
| [Variation.lean](S10Controls/Variation.lean) | Integrated first variation and stationarity/PDE equivalence for both actions. |
| [FiniteAction.lean](S10Controls/FiniteAction.lean) | Compatibility with finite total action. |
| [PhaseAverage.lean](S10Controls/PhaseAverage.lean) | Exact phase average, derivative, and stationarity equivalence. |
| [Spectrum.lean](S10Controls/Spectrum.lean) | Complete amplitude kernels, all-frequency classification, and disjointness of longitudinal/transverse spaces. |
| [Certificate.lean](S10Controls/Certificate.lean) | Integrated comparison, total/transverse dimensions, and positive/negative variational witnesses. |
| [Checks.lean](S10Controls/Checks.lean) | Concrete density distinction and D=1 coincidence. |

Coordinates, test fields, support/integration lemmas, and transverse/longitudinal
space dimensions are reused from S10Pilot. Each control's momentum, variation,
PDE, and modal operator are proved from its own density. No equality to the
baseline action is presumed.

## Verification

The control audit entry point [S10Controls.lean](S10Controls.lean) checks 32
stiffness-form declarations, many quantified over both forms, plus 16 declarations
for the [coefficient and sign controls](SCALAR_RESULT.md). Their dependency chains
contain only the standard logical axioms and no `sorryAx` or custom physics
axioms. The existing 21 S9 and 34 S10 baseline audit declarations also pass.
Warnings are errors for all pilot libraries. Captured output and source hashes
are in [CONTROLS_VERIFICATION.txt](CONTROLS_VERIFICATION.txt).

The variational witnesses establish two differences for nonzero amplitudes:
longitudinal cone waves pass both controls and fail the curl-only baseline;
transverse cone waves pass FULLGRAD, while DIVONLY admits a compact test
variation with nonzero first variation on that field.

Two isolated mutations exercise the action definitions: replacing FULLGRAD's
sum of squared derivatives by the squared divergence, and replacing DIVONLY's
square of the diagonal sum by the sum of diagonal squares. Lean rejected the unchanged quadratic-increment and stiffness-pair identities
for both mutations (exit 1). The canonical sources retain the supplied
definitions. These checks demonstrate sensitivity to those particular errors,
not to every possible mismatch between intended physics and a formal statement.

## Remaining boundary

The setting is flat whole spacetime, constant isotropic inertia, positive
coefficients, smooth backgrounds, and smooth compact variations. This increment
does not select physical D=3, justify a microscopic medium, introduce boundary
or interface laws, establish weak-solution or general Fourier/PDE completeness,
or certify the existing CAS scripts, comparator, or exports.

Anisotropic inertia now has a [separate proof report](ANISOTROPIC_RESULT.md).
Its split roots and special wavevector directions are analyzed there;
the present no-exception statements apply
to the three isotropic-inertia forms above. The sign control is covered in the
separate scalar-control report. Dimensional units,
the scalar-GNLS obstruction, confinement, and nonlinear physics remain outside
this increment. No discrepancy with the selected stiffness-control results was
found.
