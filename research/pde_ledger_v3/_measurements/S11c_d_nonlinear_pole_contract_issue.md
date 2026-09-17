# S11c-d nonlinear-pole specification issue and proposed repair

Status: **confirmed specification issue; repair awaits user approval**. The
running wider-box three-momentum quadrature may continue against its frozen
inputs. No numerical result, governing authority, physics engine, export or
Lean file is changed by this assessment.

The Lean session supplied the counterexample and proposed an analytic
tail/Abel/stability proof increment. Its `lean/s11/S11C_ANALYTIC_ASSESSMENT.md`
is an assessment/proposal, not a completed proof or numerical error certificate.

## Exact evidence

`S11c_d_nonlinear_pole_contract_probe.py` computes six rational matrix-pencil
examples without importing the physics engine. Its saved JSON contains the
actual inverses, residues, nullspaces' derivative pairings, ranks and residuals.
All examples use dimensionless synthetic data, not the S11c-d physical input.

| Pencil / enclosed roots | Integral of L^-1 L' / (2 pi i) | Idempotency residual | Meaning |
| --- | --- | --- | --- |
| z / zero | 1 | 0 | Simple-pole control |
| z^2 / zero | 2 | 2 | Counterexample to unrestricted projector definition |
| diag(z, 2z, z+1) / zero | diag(1,1,0) | 0 | Full two-dimensional semisimple control |
| [[z,-1],[0,z]] / zero | identity | 0 | A defective affine pencil can still have a genuine Riesz projector |
| z^2-1 / +1 only | 1 | 0 | Isolated simple-root control |
| z^2-1 / both roots | 2 | 2 | Summing local nonlinear modal projectors need not give a projector |

For z^2, the inverse has a second-order pole but its residue is zero. A zero
residue therefore cannot certify absence of a pole or all singular response.
The derivative pairing at zero is singular. Each applicable semisimple modal
formula in the probe agrees exactly with the directly computed residue and
local contour integral; every inverse residual is zero.

## Affected scope and dependency disposition

The shared physics authority, section 3b (lines 495–523 at the recorded hash),
calls the logarithmic-derivative contour integral a Riesz projector for every
isolated pole, extends that prescription to multiple poles without sufficient
hypotheses, and describes its perturbative preservation as Riesz-rank and
pole-motion control. The counterexample directly contradicts that unrestricted
prescription. The program brief and build directive repeat the output contract.

The existing `EndModeFrequencyData` and `EndResolventAudit.pole` compute full
left/right nullspaces and form their local residue/projector only when the
derivative pairing has full rank. They also emit idempotency residuals. This is
the relevant semisimple construction, rather than the failed unrestricted
extension. Their numerical tolerance and domain limitations remain. This
inspection does not revalidate every accepted end-mode result.

The end-resolvent contours concern constant-end normal momentum, not the
profile-dependent frequency poles of section 3b. The latter construction remains
`POLES_RIESZ_OVERLAP` in the engine's outstanding list. No completed profile
bound-pole result is being withdrawn by this issue.

The running helper loads saved action operands and evaluates
`BoundedSourceFourierQuadrature.ThreeMomentum`, using the native group/assembly
path. It does not invoke the pole constructors. Its authority, engine and helper
hashes match their immutable production preflight. No repeated quadrature is
indicated by this pole-contract finding. The authority is nevertheless pinned
as a whole file: editing it during production would break the strict source
join. Leave it unchanged until the run is validated; later changes require an
explicit provenance/dependency join, never replacement hashes alone.

## Proposed repair, not yet adopted

1. State fixed domain/codomain spaces, a local analytic sheet, isolation and
   regularity hypotheses. In the nonlocal setting provide the required analytic
   Fredholm and finite-dimensional singular-part premises rather than inherit
   finite-matrix conclusions automatically.
2. For a single isolated semisimple zero, use full nullspace bases V and W,
   D = W L'(z*) V nonsingular, R = V D^-1 W. Then P = R L'(z*) is an idempotent
   local field-space projection. Keep the simple-pole derivative normalization
   and compute all residuals. This is consistent with the residue theorem in
   [Schumacher, Corollary 7.5 and Proposition 7.6](https://arxiv.org/html/2412.15985v1#S7).
3. For higher-order/defective poles retain the full Laurent principal part and
   root-chain/partial-multiplicity data. If a genuine Riesz projector is required,
   specify a suitable linearization or realization and its physical input/output
   maps. Do not call the unrestricted logarithmic-derivative integral a projector,
   or silently discard multiple poles. See also
   [Beyn's nonlinear contour construction](https://arxiv.org/abs/1003.1580).
4. Distinguish source-to-field residue/principal-part maps from field-space
   projectors when defining overlap. Include the actual forcing injection and
   observation maps; do not treat these maps as interchangeable just because
   their finite representations are square matrices.
5. Restrict contour perturbation claims to what their hypotheses prove:
   contour invertibility and the appropriate algebraic count under an admissible
   analytic homotopy; Riesz rank only in a justified spectral realization.
   Separate pole splitting, semisimplicity and quantitative displacement bounds.
6. Align the shared authority, repeated build/output requirements and future
   SymPy/Wolfram implementations; preserve historically accepted outputs with
   their original scope. Use the exact examples above plus basis, forcing-map
   and assumption controls for the repair. Audit actual consumers before any
   recomputation. No S11c-d Wolfram implementation was found in `mathematica/`.

The user requested confirmation for a newly discovered upstream physical repair.
Only the diagnostic and proposed scope are checkpointed now. Independent
completion validation/publication of the already running action quadrature can
proceed; adoption of this pole-contract repair waits for that confirmation.

## Analytic proof handoff

Keep the proposed tail/Abel plus stability-transfer theorem first. Require its
contract to name the trial/test/outgoing spaces, cover the complete operator,
separate constant/step distributional pieces from localized profile tails, and
retain both source and observation-map errors. Gaussian witness bounds alone
do not bound scattering states or an outgoing inverse. Weak Abel convergence
must not be silently promoted to operator-norm convergence.

For an operator error epsilon and an inverse bound kappa in compatible fixed
spaces, the application must supply actual constants and verify kappa*epsilon
< 1 before using the stability estimate. Profile tail estimates by themselves
do not supply kappa. The variable-coefficient extension must preserve coefficient
gradients and interface/boundary terms. These are requirements for the proposed
proof/application connection, not claims that the new proofs already exist.

The builder report's retained solver/export contract suffix remains byte-identical
(SHA256 `f01024512c9a028e7e524e5a95b5c455be42e125d10c36f94e20493109b013c2`).
