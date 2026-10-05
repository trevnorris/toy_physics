# S9 compact statement and source connection

This record covers the bounded contract in [COVERAGE.md](COVERAGE.md).
Author: Codex. Both independent reviews are CLEAR for the fixed v1 packet.
See [FIDELITY_REVIEW.md](FIDELITY_REVIEW.md) for source pins and dispositions.

## Supplied action and native engines

The original [SymPy constructor](../../scripts/S9_light_requires_shear_sympy_audit.py)
uses `construct_curl_action(rho_br * identity3, mu_R)` with default
`stiffness_sign=-1`. The original [Wolfram constructor](../../mathematica/S9_light_requires_shear_mathematica_audit.wl)
uses `mainLagrangian = rhoBr velocityVector.velocityVector/2 -
muR curlVector.curlVector/2`, with `velocityVector = D[fieldVector,t]` and
`curlVector = Curl[fieldVector,{x,y,z}]`.

The exact parameter map sends `rho_br`/`rhoBr` to the Lean parameter `rho` and
`mu_R`/`muR` to the Lean parameter `mu`; these are positive real constants. The order of
coordinates is `(t,x,y,z)`, fields are `(u1,u2,u3)`, and
`J j i = partial_j u_i`. Thus both construct precisely

`L(J) = rho/2 sum_i J(0,i)^2 − mu/2 [(J(2,2)−J(3,1))^2 +
(J(3,0)−J(1,2))^2 + (J(1,1)−J(2,0))^2]`.

This is `S9Pilot.lagrangian` and `jetCurl`, including normalization. There is
no unspecified scalar multiple or polarization-basis change. The finite relative
action is the integral of `L(du + epsilon dh) − L(du)`, with smooth fields and
compact smooth test variations on all real spacetime. The original CAS local
calculation does not implement this analytic integration theorem.

Both engines use `exp(i(k·x − omega t))` for route A. Lean uses its real cosine
quadrature; for this real constant-coefficient second-order operator both give
the same real modal matrix after removing the phase factor:

`M(omega,k) = rho omega² I − mu (|k|² I − k kᵀ)`.

Their Euler–Lagrange sign is `dL/du − sum_j partial_j(dL/dJ_j)`; Lean derives
this sign from the actual integrated relative action. SymPy route B takes the
mixed Hessian of the opposite-phase pair. Wolfram takes that paired kernel,
forms `aᵀ K a / 2`, then differentiates twice. Since this kernel is symmetric,
both yield the same `M`, without a phase-average factor. This S9 paired-wave
normalization should not be confused with S10's real-cosine period average.

The compact source check (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S9_lean_source_check.py`) executes
only selected original SymPy constructor and route definitions. It compares the
action in independent jet symbols and the two native modal matrices with these
handwritten references, by exact symbolic zero residuals. It checks the
determinant `rho s (rho s − mu |k|²)^2` for `s=omega²`, and rejects a deliberately
reversed shear sign. Wolfram is connected by inspected native source definitions.
The instrument checks five literal anchors: coordinates, velocity, curl, MAIN
Lagrangian and wave phase. Those anchor checks do not cover every Wolfram line,
including the Euler–Lagrange sign and route-B normalization. The whole Wolfram
file is separately hash-pinned in the reviewed packet, including those lines.

The instrument records Lean source hashes for provenance; it does not parse
Lean or compare those hashes to an accepted revision. Both reviewers checked
the handwritten references and the source interpretation, and the author
recomputed the packet and live source hashes before closure. A later changed
source requires a new correspondence assessment; rerunning this instrument
alone does not renew fidelity clearance. This is **not** a fresh Wolfram
execution, a Lean-certified translator, or a new comparator/PIT claim.

## Exhaustive mathematical coverage

For positive `rho,mu` and nonzero real `k`, the cone value is strictly positive.
The existing S9 theorems give the entire stationary amplitude space:

| Frequency case | Kernel | Dimension |
|---|---|---|
| `omega = 0` | `longitudinalSpace k = span{k}` | 1 |
| `omega ≠ 0` and `omega² = (mu/rho)|k|²` | `transverseSpace k = ker(dot k)` | 2 |
| `omega ≠ 0` and `omega² ≠ (mu/rho)|k|²` | zero | 0 |

The last case follows by contradiction from `propagating_mode_iff` for any
alleged nonzero amplitude. These cases are disjoint and exhaustive by equality
case splits, including both signs of frequency. The zero vector is admitted in
every kernel; a nonzero mode must satisfy the stated cases. The determinant's
multiplicities agree here, but they are not substitutes for the proved full
subspace dimensions. CAS census checks must cover the zero root and cone root
separately and verify the full kernel/transverse intersection, not one witness.

Reuse `s9_variational_certificate`, `propagating_variational_mode_iff`,
`zero_frequency_iff` and the two `*_finrank` results. The exact S10 specialization
in [Specialization.lean](../s10/S10Pilot/Specialization.lean) identifies the
action, modal operator, local PDE, relative action and stationarity. Its existing
scalar coefficient, sign and full-subspace controls apply through that identity;
S10 reviews did not themselves review the original S9 engine constructors.

## Scalar phase claim and its limit

The supplied law is the dimensionless-phase convention
`v = (hbar/m) grad(theta)` recorded in `docs/model_map.md` and S9's inputs.
Here its spatial domain is D=3, as in this S9 pilot. A positive, uniform density
and a smooth single-valued phase provide the intended local Madelung regime.
No density dynamics are supplied or proved by this module.

For arbitrary real `theta0, A, B, omega` and `epsilon`, set
`theta = theta0 + epsilon [A cos(phi) + B sin(phi)]`,
`phi = k·x − omega t`. [Madelung.lean](S9Pilot/Madelung.lean) differentiates
the supplied law in space and then in `epsilon` at zero. The resulting cosine
and sine velocity amplitudes are respectively `(hbar/m) B k` and
`−(hbar/m) A k`. Both belong to `span{k}` for every choice of phase quadratures.
For nonzero `hbar,m` the cosine-amplitude range is the entire longitudinal span;
for nonzero `k` it has dimension 1 by the existing span theorem. This counts
polarization directions, **not** dynamical branches. It is compatible with two
real quadrature coefficients for a single spatial polarization.

For `k ≠ 0`, `span{k} ∩ transverseSpace k = {0}`: transverse velocity quadratures
must therefore vanish. At `k = 0` the velocity is identically zero, handled
explicitly rather than applying the nonzero-wavevector count. An admissible
unit-coefficient example produces a nonzero longitudinal velocity. No claim is
made that zero velocity forces the scalar density perturbation to vanish.

The result holds before imposing a scalar evolution equation and therefore
restricts any solutions inside this ansatz. It does not derive the GNLS
equation, its Bogoliubov dispersion, exactly one dynamical branch, a global
Helmholtz theorem, emergent collective-mode classifications, bulk confinement,
or the full broader P2 sentence that a scalar theory has no spin-1 excitations.
The supplied constitutive law and physical interpretation remain premises.

## Control map

The closure instrument leaves every canonical proof unchanged. Its source
mutations use the existing mathematical proofs; its concrete statement mutants
have separately compiled passing counterparts. Review the full source and
diagnostics in the result record (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S9_lean_contract_checks.json`).
The final suite passes: nine mathematical rejections, five explicit positive
controls, the sequential S9 rebuild and 32 axiom audits. Both independent reviews
are CLEAR; see the [review dispositions](FIDELITY_REVIEW.md).

| Load-bearing claim | Deliberate mutation / passing evidence |
|---|---|
| Action shear sign and first variation | Change only the supplied density's shear sign; `lagrangian_variation` must reject the resulting polynomial identity. Original `Action.lean` compiles. The native SymPy constructor check independently detects the opposite sign in both the action and modal operator. |
| Sign of the integrated action change | Reverse both terms in the claimed `relativeAction_expansion` while retaining the actual action and proof; the resulting equation must fail. Original `Variation.lean` compiles. |
| Exclusion of additional nonzero-frequency roots | Reuse the action/operator shear-sign controls above and the S10 coefficient/sign controls through exact specialization. The exact native determinant check independently identifies the root set. Mathematically, `propagating_mode_iff` excludes every nonzero off-cone amplitude. There is no separate off-cone fixture; this is shared control coverage of the same operator and cone, not an additional mutation. |
| Full transverse/static dimensions | At the admissible wavevector `(0,0,1)`, assert dimensions 3 instead of 2 and 2 instead of 1. Correct counterparts use the unchanged general dimension theorems and compile. The mutants reduce to false arithmetic. |
| Nonzero-wavevector restriction on the static census | At `k=0`, the vector `(0,0,1)` is stationary but is not in `span{0}`. The passing conjunction establishes this excluded case; falsely putting that vector in `span{0}` must fail. |
| Actual spatial gradient, not a supplied parallel amplitude | Replace the spatial derivative by the time derivative in the velocity definition; the actual phase-velocity identity must fail. The canonical Madelung module compiles. |
| Both real velocity quadratures and their relative sign | Reverse the sine velocity amplitude's sign; `linearVelocity_eq` must fail. Both quadratures are derived from the supplied phase in the passing module. |
| No nonzero transverse velocity polarization | The unit-coefficient cosine amplitude at `(0,0,1)` is not transverse. Its passing check and the false transverse-membership counterpart test this consequence nonvacuously. |
| Nonzero prefactor in the complete longitudinal amplitude range | At `hbar=0`, the cosine amplitude is zero. The passing check admits this limit; falsely asserting nonzero amplitude there must fail. Physically positive `hbar,m` and a nonzero-velocity witness ensure the intended family is nonempty. |

The initial instrument attempts (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S9_lean_contract_initial_attempt.json`)
are retained transparently. Their canonical rebuild and all 32 selected axiom
audits passed. The instrument refused to count a `simp`-only failure for a
reversed action definition and an `omega`-only failure for a wrong general
count. The final controls instead expose a wrong-sign integral expansion and
false concrete counts. These were changes to negative-test design, not repairs
to the mathematical definitions or theorems. Two positive-fixture elaboration
issues (a finite-vector component simplification and pointwise equality of a
zero operator) also stopped their runs before the corresponding negative test
could count. The fixtures were corrected; canonical proofs stayed unchanged.
Verified builds are reused only
when all canonical source hashes and build commands match.
