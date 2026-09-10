# S11c-d transverse sign diagnosis

Historical evidence recorded before the user-authorized inertia repair. The
equations, source locations, hashes, and run measurements below describe that
snapshot. Its transcript is preserved at
[`S11c_d_transverse_sign_diagnostic.before_inertia_repair.out`](../scripts/out/S11c_d_transverse_sign_diagnostic.before_inertia_repair.out).
The diagnostic command now reads the current exports. See the
[repair record](S11c_inertia_repair_report.md) for the repair and regenerated results.

2026-09-10. User-authorized focused investigation of the missing transverse
incident channel. This is a diagnostic, not a completed S11c-d build or a
cross-engine adjudication. Production engines and ledger exports are unchanged.

## Finding

The SymPy S11c-b slab assembly has a relative inertia/stiffness sign error in
the tested conservative sector. It is present before the c2 closure and before
the d Fourier reduction. The c2 closed rows preserve it on all tested transverse
ansatzes. It is not removed by an overall row sign or by choosing the opposite
harmonic time convention.

Let `K_T = mu_R + mu_S` denote the transverse stiffness computed from **S11c-b's
own exported energy**. For `u_p = epsilon q(x_d,t)`, `p != d`, the exported
equation and the independently varied supplied action give:

| Construction | Equation divided by epsilon |
|---|---|
| Exported S11c-b; exported c2 on the same ansatz | `-rho_br q_tt - K_T q_xx` |
| Euler variation of the supplied `T - U` | `+rho_br q_tt - K_T q_xx` |
| Exported minus action-derived | `-2 rho_br q_tt` |
| Exported plus action-derived | `-2 K_T q_xx` |

The stiffness coefficients agree while the inertia coefficients have opposite
signs. Aligning the rows by their stiffness coefficient leaves the same nonzero
inertia residual. The stored equation instead agrees, up to its overall sign,
with variation of `T + U` on this slice. No irreversible response is placed in
either diagnostic action: this is the uniform unforced transverse subsystem.

For either `exp(i k x - i omega t)` or `exp(i k x + i omega t)`, differentiation
of the ansatz yields:

| Construction | Fourier equation | Computed omega squared |
|---|---|---|
| Stored | `K_T k^2 + rho_br omega^2` | `-K_T k^2/rho_br` |
| Supplied action | `K_T k^2 - rho_br omega^2` | `+K_T k^2/rho_br` |

Thus positive `K_T` and positive density give exponential transverse growth or
decay in the stored conservative dynamics, where the supplied positive transverse
energy gives an oscillatory mode. No positivity assumption on the entire energy
basis is needed: the admitted pure-curl subcase `mu_S = 0`, `mu_R > 0` already
supplies the counterexample.

## Independent energy check

The diagnostic reads the retained stored-energy terms from `energy_basis_variable`
and applies the same transverse ansatz to them. Kinetic energy is constructed
from the supplied `rho_br |u_t|^2/2`. For the real standing wave
`q(x,t) = Q(t) cos(k x)`, integration over one dimensionless spatial period gives
the mean energy per volume:

    E = epsilon^2 [rho_br Q_t^2 + K_T k^2 Q^2]/4.

Solving each measured equation for `Q_tt` and substituting it into `dE/dt` gives:

    stored equation:         epsilon^2 K_T k^2 Q Q_t
    action-derived equation: 0

This uses the supplied `T + U` for the energy and independent differentiation of
the action for the equation. It does not define an energy exchange from the same
operator and then declare balance. The c2 theta and thickness responses vanish
on these transverse ansatzes, so no exchange with those fields accounts for this
unforced energy change in the tested slice.

## Trace to source

1. `scripts/S11c_b_brane_operator_sympy_audit.py:2290` constructs the stored-energy
   Euler derivative as `u_local - divergence(u_flux)` and **subtracts**
   `epsilon rhobr u_tt` at line 2325. The thickness construction makes the
   corresponding subtraction at line 2363. The acceleration symbols are declared
   as second time derivatives at lines 287–291.
2. The later constraint fold strips those negative kinetic terms by adding them
   back at lines 2894 and 2899, folds the material constraint on the internal
   variation, then again subtracts inertia at lines 2946 and 2970. The transverse
   result is already affected before the constraint reaction.
3. c2 reads the exported b rows at
   `scripts/S11c_c2_selfenergy_fold_sympy_audit.py:541`, substitutes the face
   closure at line 552, and maps jets to field derivatives at lines 557–558.
   `wave_jet` differentiates twice in time at line 170. The measured
   closed-minus-open transverse residual is zero for each tested case/axis.
4. The preceding d preflight independently used `EdgeReduction.value` and
   `ConstantEndPencil.strong_matrix` on the closed rows. Its displacement and
   weak curl-potential results had the same factor, with zero lift residual.
   This investigation locates that factor upstream of the reduction.

The actual b source SHA-256 is
`2a1fef275636dd210c8864ab1762859f63a4dfb054008fa21620e54b14e47fd2`, matching the
generator digest pinned by `S11c_b_exports.py`. This is not explained by a stale
export generated from a different b source. Git history places the raw negative
inertia insertion in `b17587de`; the later constraint-fold handling was added in
`82f53828`. These identify code history, not an independent validation.

## Thickness corroboration and the S11b comparison

The diagnostic applies the supplied uniform no-flux, zero-spatial-momentum
constraint `theta = -e_W` to the energy before varying. It removes prescribed
pressure and reciprocal-response forcing from the open operands and compares
the result with the corresponding exported b thickness row. Each case gives:

    stored - action-derived = -2 epsilon mu_W W_0^2 e_W,tt.

This corroborates that the same inertia assembly affects both mechanical rows.
It does not validate or modify pressure-work, closure-fold or nonlocal-response
signs in the full problem.

S11b's exported transverse equation has stiffness `mu_R + mu_S/2`, whereas b's
exported energy and equation use `mu_R + mu_S`. Both stiffness operands and their
difference `mu_S k^2/2` are emitted. That normalization change is an additional
unresolved inheritance issue; it is not silently removed by a coefficient
rescaling. For the sign check, the explicit `mu_S = 0` control makes the energy
convention common. After exposing the symbol-identity alignment, S11b minus the
action-derived Fourier equation is zero; stored b minus S11b is
`2 rho_br omega^2`.

## Relation to the recorded debt

`steps/S11c_b_variable_coefficient_operator.md:113` already lists the kinetic
`-K/+K` difference as unadjudicated. The d physics packet distinguishes it from
response-slot signs at `directives/S11c_d_SHARED_PHYSICS.md:164`. An increment or
response-slot comparison cannot establish the correctness of an unchanged
diagonal inertial term. In d's full operator, it affects channel existence.

The present result settles the tested SymPy conservative relative sign using
the supplied energy and action. No Wolfram engine or comparator was run, and
the c2 carrier/source/Phi families, face signs, and other recorded debts remain
unadjudicated. Uniform-regression code exists at
`scripts/S11c_b_brane_operator_sympy_audit.py:4963`; this investigation makes no
claim about which historical review or run inspected its residual.

## Repair boundary

S11c-d should remain paused while the upstream mechanical assembly is repaired.
A downstream sign flip or an artificial negative-stiffness sample would conceal
the defect without restoring the inherited model.

Derive the kinetic Euler term from the supplied `T` and combine it consistently
with the constrained stored-energy variation and prescribed external work.
Audit the raw assembly, kinetic removal before the constraint fold, and
reinsertion afterward together: changing only the raw minus signs leaves the
add-back logic inconsistent. Both u and e_W are affected. Keep the material
virtual constraint and sourced mass equation intact. Resolve the separate
`mu_S` normalization through an explicit energy-operand comparison.

After a source repair, regenerate the affected b export and downstream
dependents, including c2 and d, with input digests refreshed. Check c1's actual
dependencies rather than assuming it is unaffected. The orchestrator owns
review and cross-engine follow-up. The diagnostic added here provides small,
directly executable conservative regressions for that work.

## Reproduce and inspect

From `research/pde_ledger_v3`:

```bash
python -u scripts/S11c_d_transverse_sign_diagnostic.py > scripts/out/S11c_d_transverse_sign_diagnostic.out
```

The default path freshly loads the real three-parent fold and checks its direct
lookup manifest. `--development-cache` is explicitly marked exploratory and was
not used for the final artifact. Physical emissions carry computed grades and
restored units; residual values are printed, not asserted. A structural
dimension guard runs after its residual/unknown-unit records are printed.

The computation map is in `PY_S11CD_SIGN_EMISSION_LINES`. Action construction is
at lines 179–181; open/closed actions at 182–185; row-sign comparisons at
196–206; both Fourier conventions and the frequency solve at 208–214; S11b
comparison at 220–232; standing-wave energy at 233–249; thickness at 250–261.

The final fresh run completed successfully in 94.523 s with peak RSS
1,713,536 KiB and empty stderr. Its 433,736-byte output contains 403 unique
emission tags, including 394 physical records. The cache flag is false. All
recorded input digests still match the files, and the printed dimensional
constraint and unresolved-unit collections are both empty.

Parsing the completed output verifies 24 transverse restrictions (four cases,
six ordered propagation/polarization pairs). Every closed-minus-open vector
and closed scalar row vanishes; every stored-minus-action row has the stated
inertial discrepancy. Both time-phase conventions give the same pair of
opposite-sign frequency-squared solutions in all 48 evaluations. All four
standing-wave energy checks and all four thickness checks reproduce the
results above. The four pure-curl comparisons with S11b vanish for the
action-derived equation, while the stored equation differs by
`2 rho_br omega^2`.

SHA-256 of the final diagnostic:
`24586e3c97034cdd50230b5e6c404cc2135fb4a6554735ef0f50d7d0412ad6aa`.
SHA-256 of its completed output:
`961dc00bae81b77e1f4833fbba5d5e43d3b566c1add06db42c7d6d4fa3d025fc`.

No production engine, upstream export, or existing d output was changed, and
no commit, review leg, comparator, or Wolfram run was launched. The full d
build and export remain incomplete pending the upstream repair described above.
