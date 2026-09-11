# S10 one-axis anisotropic inertia — 2026-09-11

The one-axis inertia control is formalized from its supplied density through
integrated variation to a complete real plane-wave amplitude classification.
The proofs confirm the generic S10 counts and its parallel-wavevector stratum.
They also establish an additional perpendicular-wavevector stratum: the two
frequencies remain distinct, but every propagating amplitude is exactly transverse.

## Supplied action and setting

The action input is `XFORM_ANISO` in section 7 of
[S10_SHARED_PHYSICS.md](../../directives/S10_SHARED_PHYSICS.md):

```text
L(J) = (rho/2) Σ_i weight_i J_time,i² - (mu/2) S_curl(J),
weight_e = sigma,   weight_i = 1 for i != e,
S_curl(J) = (1/2) Σ_i Σ_j (J_spatial_i,j - J_spatial_j,i)².
```

`sigma` is the specification's `s_rho`. The distinguished index `e : Fin D`
corresponds to its first spatial component; allowing any distinguished index
changes only the labeling. The density uses the weighted sum directly. Lean
proves its identity with the baseline density plus
`(rho/2)(sigma-1) J_time,e²`. Only one inertia coefficient changes.

Assume constant `rho>0`, `mu>0`, `sigma>0`, `sigma != 1`, nonzero real k,
and real `u(t,x)=a cos(k·x-omega t)`. Spacetime is flat R^(D+1), the background
is smooth, and variations are smooth vector fields with compact support.
The physical action, field content, rest background, and dimension remain inputs.

## Derived operator and complete spectrum

The actual density derivative and the independently differentiated plane-wave
Euler–Lagrange expression give the same operator:

```text
M a = rho omega² [a + (sigma-1) a_e E] - mu [K a - k(k·a)],
K = |k|²,   E = unit vector along e.
```

The EL expression is `cos(phase) M a`. For nonzero frequency, stationarity
implies the weighted transverse condition

```text
k·a + (sigma-1) k_e a_e = 0.
```

This condition need not imply Euclidean transversality `k·a=0`.
Define

```text
p = k_e,
q = K - p² = |k - p E|²,
A = q + sigma p²,
w = A E - sigma p k.
```

Lean proves `q>=0`, `q=0` iff `k=p E`, `w_e=q`, and
`k·w=(1-sigma)pq`. If `q != 0`, the complete nonzero-frequency classification is:

| Candidate squared frequency | Full amplitude space | Dimension |
|---|---|---:|
| `(mu/rho) K` | `a_e=0` and `k·a=0` | D-2 |
| `(mu/rho) A/sigma = (mu/rho)(p²+q/sigma)` | span{w} | 1 |

Both frequency values are positive and distinct under these assumptions.
The proof treats `z=rho omega²/mu` as an arbitrary real parameter before
specializing to real frequencies: every nonzero real root with a nonzero
amplitude is positive. The certificate supplies positive real frequencies and
identifies their entire amplitude spaces with integrated action stationarity.
Every other nonzero real frequency has no nonzero stationary amplitude.

At zero frequency the full amplitude space is always span{k}, of dimension 1,
and its exactly-transverse intersection has dimension 0. The static longitudinal
degree of freedom survives the inertia change.

## All directional cases

N2 below is the dimension of the full amplitude kernel. N3 is the dimension of
its intersection with `k·a=0`, not a count inferred from a displayed basis.

| Direction of k | Ordinary frequency: N2 / N3 | Extra frequency: N2 / N3 | Total propagating N2 / N3 |
|---|---|---|---|
| Oblique: `p != 0`, `q > 0` | D-2 / D-2 | 1 / 0 | D-1 / D-2 |
| Perpendicular to E: `p=0`, `q>0` | D-2 / D-2 | 1 / 1 | D-1 / D-1 |
| Parallel to E: `q=0`, `k != 0` | One merged frequency, D-1 / D-1 | Same frequency, counted once | D-1 / D-1 |

These cases exhaust nonzero real wavevectors. Lean proves frequency coincidence
occurs exactly when `q=0`; for `q != 0`, it proves the extra polarization is
exactly transverse exactly when `p=0`. Thus a change in N3 can occur without
a frequency coincidence or a change in N2.

At D=3 the oblique count is `(1/1)+(1/0)`, and at D=4 it is
`(2/2)+(1/0)`, agreeing with the generic S10 record. In the perpendicular case
the last summand is instead `1/1`. The total mode count remains D-1.

The formulas cover arbitrary dimension with a distinguished index. At D=2
and `q>0`, the ordinary candidate frequency has zero nullity and is not an
actual mode; the extra mode exists. At D=1, no nonzero real-frequency mode
exists. Existence theorems give nonzero extra amplitudes whenever `q>0`, and
nonzero ordinary amplitudes when additionally D>=3. The exact `sigma=1`
limit recovers the baseline density and integrated stationary fields.

## Additional stratum relative to the ledger record

[S10_two_transverse_photons.md](../../steps/S10_two_transverse_photons.md)
lists the parallel stratum `k_2=...=k_D=0`, `k_1 != 0`. Its generic oblique
counts remain correct. The formal result additionally identifies
`k_1=0`, `k != 0`, with the same allowed `sigma>0`, `sigma != 1` assumptions,
as an N3 stratum. No inertia coefficient is returned to its isotropic value.

A concrete checked example has D=3, `rho=mu=1`, `sigma=2`, and `k=(0,1,0)`.
Its ordinary squared frequency is 1, its extra squared frequency is 1/2,
and `a=(1,0,0)` is a nonzero stationary extra amplitude with `k·a=0`.
Lean proves the extra kernel's transverse dimension is 1. For comparison,
`k=(1,1,0)` has a stationary extra amplitude `(1,-2,0)` at squared frequency
3/2, with `k·a=-1`.

The subsequent [focused CAS repair](../../_measurements/S10_anisotropic_strata_report.md)
adds stacked-matrix rank-drop discovery to both engines and reproduces these
counts at explicit points in D=3 and D=4. Its dedicated comparison passes for
ten roots across the two directional strata. The original broad comparator and
downstream export snapshots retain their separate, narrower evidence status.

## Variational chain and files

The relative action integrates the pointwise density change under a compact
variation. The proofs establish its integrability, differentiate its exact
quadratic expansion, justify integration by parts, and use the fundamental
lemma plus continuity to obtain pointwise EL vanishing from all compact tests.
Stationarity is therefore against arbitrary admissible test fields, including
those outside the plane-wave ansatz. Finite total action is addressed separately
with an explicit background integrability hypothesis. Nondecaying plane waves
do not require finite total action.

The phase average is an actual integral over 0..2π and is one half of the modal
action; its amplitude derivative retains the same factor. No division by omega
is used to define a time period.

| File | Role |
|---|---|
| [Action.lean](S10Anisotropic/Action.lean) | Direct weighted density, exact variations, and modal operator. |
| [PlaneWave.lean](S10Anisotropic/PlaneWave.lean) | Derived momenta, coordinate EL expression, plane-wave evaluation. |
| [Analytic.lean](S10Anisotropic/Analytic.lean) | Smoothness of the anisotropic momenta and EL expression. |
| [Variation.lean](S10Anisotropic/Variation.lean) | Integrated first variation and stationary-field/PDE equivalence. |
| [FiniteAction.lean](S10Anisotropic/FiniteAction.lean) | Compatibility with finite total action. |
| [PhaseAverage.lean](S10Anisotropic/PhaseAverage.lean) | Exact phase integral and its derivative. |
| [Geometry.lean](S10Anisotropic/Geometry.lean) | Axis decomposition, rank-nullity, positive/distinct frequency values. |
| [Spectrum.lean](S10Anisotropic/Spectrum.lean) | Exhaustive split, parallel, and static amplitude classification. |
| [Census.lean](S10Anisotropic/Census.lean) | Kernel dimensions and transverse intersections on all directional cases. |
| [Certificate.lean](S10Anisotropic/Certificate.lean) | Integrated classification, positive roots, frequency coincidence, isotropic limit. |
| [Checks.lean](S10Anisotropic/Checks.lean) | Nonzero mode existence, concrete direction checks, and D=1 exclusion. |

Coordinate, test-field, integration, and baseline linear-algebra helpers are
reused from S10Pilot. The existing S9, S10 baseline, and stiffness-control Lean
sources are unchanged. The shared environment has one additional library target,
`S10Anisotropic`, under `lean/s10/`; dependency pins are unchanged.

## Verification and boundary

The [root audit](S10Anisotropic.lean) checks 48 selected declarations, including
the two later [scaling results](S10Anisotropic/Scaling.lean). All depend
only on the standard logical axioms `propext`, `Classical.choice`, and
`Quot.sound`; there are no admissions or custom physics axioms. The full build
also checks 202 declarations in S9, the S10 baseline, the other controls,
and the later [Q6/Q7 audit](Q6_Q7_RESULT.md), for 250 in total.
Warnings are errors.
Captured output and source hashes are in
[ANISOTROPIC_VERIFICATION.txt](ANISOTROPIC_VERIFICATION.txt).

An isolated mutation replaced the one-axis weighted sum with uniform rescaling
of every inertia coefficient. Lean rejected the unchanged kinetic identity
(exit 1). A concrete canonical check separately verifies that a unit velocity
along the untouched second axis has kinetic norm 1, while one along the first
axis has kinetic norm 2 when sigma=2. These checks test this particular modeling
error; they do not certify every possible transcription of the intended physics.

This increment does not select physical D=3, derive the action from a microscopic
model, certify the CAS implementations or their stratum enumeration, establish
complex spectral or general Fourier/PDE completeness, or introduce boundaries,
weak solutions, dissipation, confinement, or nonlinear dynamics. Sign controls
now have a [separate proof](SCALAR_RESULT.md); action and root units are covered
by the later [Q6/Q7 audit](Q6_Q7_RESULT.md). The scalar-GNLS obstruction remains
separate work.
