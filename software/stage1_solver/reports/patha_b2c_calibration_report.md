# Path-A B2c Calibration

Overall B2c status: `root_bounded_by_resolved_floor`
Machine bundle: `software/stage1_solver/runs/patha_b2c_calibration/bundles/patha_b2c_calibration_bundle.json`
Bundle content hash: `bcf83db1850c54de4cc01290bfa4aee8e7b1e8e3a769f8e2d2e8e2cb546683bd`

## Verdict

Measured, fit-independent verdict: MISS. At every converged tau, the directly measured `P0` is far below the target `54/5=10.8` on the stable side; this is measured, not extrapolated.

| tau | P0 | orders_below_target | N0/10.8 | K |
| --- | --- | --- | --- | --- |
| 1 | 2.795e-09 | 9.59 | 1.998e-09 | 7.721e+00 |
| 0.03 | 2.138e-06 | 6.70 | 4.576e-08 | 2.316e-01 |
| 0.02925 | 1.461e-06 | 6.87 | 3.051e-08 | 2.258e-01 |
| 0.029125 | 1.224e-06 | 6.95 | 2.545e-08 | 2.249e-01 |

The measured deficit is 6.7-9.6 decimal orders, i.e. the ~7-to-9 order deficit is the headline result. Algebraically `P0=N0/D0`, so any stable root must satisfy `D0(tau*)=N0(tau*)/10.8`. In the measured records `N0/10.8` is only `1.998e-09` to `4.576e-08` while `K` is `2.249e-01` to `7.721e+00`; reaching the target therefore requires `D0 << K`, a knife-edge cancellation, regardless of the precise fitted tau*.
The `tau=1` structural win, meaning no fragile cancellation, does not survive calibration. The GR anchor can only be approached by riding the `D0 -> 0` stability boundary, so this is not a natural calibration.
The rigorous edge statement is only the upper bound `tau* < tau_floor ~= 0.0291`. Drift-aware estimates are fit-dependent: frozen deepest coefficients give `tau* ~ 3.43e-05`, the wide `tau=1` to `tau=0.03` power-law pair gives `tau* ~ 6.97e-04`, and the deepest close-pair fit is `above_floor_degenerate_not_a_valid_bound` (`tau* ~ 3.59e-02`). Thus tau* is not pinned by this run.
No stable-side R_norm sign bracket exists in the real converged re-solve set. The deepest full bundle is tau=2.912500000000e-02 with R_norm=-1.079999877620e+01; any root, if present on the stable branch, is only bounded above by that residual-gated floor and is not confirmed.
Numerical-floor evidence: No D0 or Jacobian-smallest-mode measurement exists at the failing tau because that solve did not converge. The deepest residual-gated full bundle has D0=O(0.1) and K/(B0+Z0) far from one; the lower near-floor failures are smooth line-search/continuation stalls with residuals rising monotonically as tau drops. This is consistent with a numerical (continuation/conditioning) floor, not a diagnosed physical marginal-stability edge.
The real dual-engine root finder was attempted on the real converged sample set, but no stable-side sign bracket exists, so there is no `root_agreement` to report.
Numerics status: the requested `1e-3`-`1e-2` band was not reached; the deepest converged tau is `0.029125`, about `2.9x` above the top of that band. Warm-start machinery is functional but converted no new converged point here: the deepest full evaluation converges cold, and the floor descent from `0.03` to the current floor is only about `2.9%`. This is a known limitation, not a material improvement; pinning tau* would require a small-tau preconditioner/linear-solver effort, deferred and not required for the MISS verdict.

No background, BdG, or Maxwell coefficient bundle is extrapolated below the recorded floor.

## Tau Floor

Lowest converged closed-background probe/evaluation: `tau_floor=2.912500000000e-02` with residual `3.950424206506e-10` from `full_evaluation`.

Converged floor/probe rows:

| tau | source | residual_linf | message |
| --- | --- | --- | --- |
| 2.912500000000e-02 | full_evaluation | 3.950424206506e-10 | closed B2a continuation completed |
| 2.912500000000e-02 | background_probe | 3.276117224726e-11 | closed B2a continuation completed |
| 2.925000000000e-02 | full_evaluation | 3.596428710395e-14 | closed B2a continuation completed |
| 2.925000000000e-02 | background_probe | 5.372991981889e-12 | closed B2a continuation completed |
| 2.950000000000e-02 | background_probe | 3.631851763775e-13 | closed B2a continuation completed |
| 3.000000000000e-02 | full_evaluation | 3.387974345515e-10 | closed B2a continuation completed |
| 3.000000000000e-02 | background_probe | 3.387974345515e-10 | closed B2a continuation completed |
| 1.000000000000e+00 | validated_seed | 2.749412086889e-10 | closed B2a continuation completed |

Failed or dropped tau rows:

| tau | source | residual_linf | message |
| --- | --- | --- | --- |
| 5.055913803509e-06 | background_probe | 7.523479070340e-04 | continuation stopped at eos_K=0.08: line search failed to reduce residual |
| 1.000000000000e-04 | background_probe | - | external timeout after 600 seconds |
| 1.000000000000e-02 | background_probe | 4.370505768220e-03 | continuation stopped at eos_K=0.08: line search failed to reduce residual |
| 2.000000000000e-02 | background_probe | 2.736135207282e-03 | continuation stopped at eos_K=0.08: line search failed to reduce residual |
| 2.800000000000e-02 | background_probe | 4.918899520096e-04 | continuation stopped at eos_K=0.18: line search failed to reduce residual |
| 2.900000000000e-02 | background_probe | 5.355903474474e-05 | continuation stopped at eos_K=0.18: line search failed to reduce residual |
| 2.906250000000e-02 | background_probe | 2.423840866221e-05 | continuation stopped at eos_K=0.18: line search failed to reduce residual |

## Real Resolved Evaluations

| tau | kind | R_norm | D0 | stable | background_residual_linf |
| --- | --- | --- | --- | --- | --- |
| 2.912500000000e-02 | resolved_deepest | -1.079999877620e+01 | 2.246076541249e-01 | true | 3.950424206506e-10 |
| 2.925000000000e-02 | resolved_deep | -1.079999853903e+01 | 2.255324765600e-01 | true | 3.596428710395e-14 |
| 3.000000000000e-02 | floor_resolved | -1.079999786247e+01 | 2.312081800702e-01 | true | 3.387974345515e-10 |
| 1.000000000000e+00 | seed_validated | -1.079999999721e+01 | 7.720905323816e+00 | true | 2.749412086889e-10 |

## Frozen Negative Control

Frozen-background root from the validated `tau=1` coefficients: `tau_frozen ~ 5.06e-06`; the frozen Schur edge is also `~5.06e-06`.
That prior is below the resolved floor: `True`. It is recorded only as the R6 negative control and is not used by the final status logic.
Kappa provenance note: local analytic formulas use continuum `kappa_hat=7.7209353105`; direct `K/tau` in the B2a-seeded records is `7.7209443582` (relative difference `1.17e-06`, discretization scale).
Measured resolved coefficient drift from `tau=1.000000000000e+00` to `tau=2.912500000000e-02` has max relative movement `1.173945078456e+01` across `{B,Z,N}`.

| coefficient | relative_drift |
| --- | --- |
| B0 | 5.784060921736e+00 |
| B2 | 5.426484583413e+00 |
| B4 | 5.174704939556e+00 |
| Z0 | 7.398111543305e+00 |
| Z2 | 6.938270275575e+00 |
| Z4 | 6.285883685425e+00 |
| N0 | 1.173945078456e+01 |
| N2 | 9.589633978695e+00 |
| N4 | 7.500655369108e+00 |


## Edge Diagnostic

Deepest full bundle: `tau=2.912500000000e-02`, `R_norm=-1.079999877620e+01`, `D0=2.246076541249e-01`, `P0=1.223800242913e-06`, `K/(B0+Z0)=8.490550997143e+02`, closed residual `3.950424206506e-10`.
Nearest failed probe below the floor: `tau=2.906250000000e-02`, residual `2.423840866221e-05`, message `continuation stopped at eos_K=0.18: line search failed to reduce residual`.
Deepest frozen-coefficient local edge estimate, reported only as one fit in the spread block: `tau* ~ 3.43e-05`, `D0(tau*) ~ 2.55e-08`, `K/(B0+Z0) ~ 1`, cancellation digits `~4.0`. This is fit-dependent and does not pin the edge.

Fit-spread estimates recorded in the machine bundle:

| fit | tau_estimate | status |
| --- | --- | --- |
| frozen_deepest_coefficients | 3.43e-05 | below_floor_fit_dependent |
| wide_pair_power_law_tau_1_to_0p03 | 6.97e-04 | below_floor_fit_dependent |
| close_pair_power_law_deepest_two | 3.59e-02 | above_floor_degenerate_not_a_valid_bound |

## Deepest Spot Convergence

BdG modal sweep at deepest τ: exported K=`30`, all-positive modes `100`, max truncation `1.453063862381e-05` against tolerance `1.000000000000e-04`, passed `True`.
Maxwell grid refinement `47x17` to `51x19` max relative Z/N movement `7.599119219021e-03`.

| coefficient | relative_change |
| --- | --- |
| Z0 | 5.540174219140e-04 |
| Z2 | 2.097497451521e-03 |
| Z4 | 7.599119219021e-03 |
| N0 | 8.117504215265e-04 |
| N2 | 1.196619699708e-03 |
| N4 | 7.239415178523e-03 |

## Naturalness At Bound

At the deepest residual-gated upper bound `tau=2.912500000000e-02`: `|ln tau|=3.536158367096e+00`, `K/(B0+Z0)=8.490550997143e+02`, `D0=2.246076541249e-01`, cancellation fraction `9.988222201358e-01` (`0.001` digits).
Because no real root is confirmed, this is bound/floor naturalness, not naturalness at a calibrated root.

## Held-Out Quantities

No calibrated `tau*` is confirmed in this run, so `R_pole/P2/P4` are not emitted as calibrated held-out predictions. The deepest resolved full bundle is reported below as a diagnostic only.

| tau | R_norm | D0 | R_pole | P2 | P4 |
| --- | --- | --- | --- | --- | --- |
| 2.912500000000e-02 | -1.079999877620e+01 | 2.246076541249e-01 | -2.546642778538e-03 | 3.773260192248e-08 | 1.291491506575e-07 |

## Error Bars

Recorded B2a/B2b budgets were propagated by coefficient perturbation at the deepest residual-gated bundle. `K` is exact for the frozen Hooke family; no error bar here converts the bound into a confirmed root, and no local tau* error bar is reported because tau* is not pinned by this run.

| quantity | value | abs_error |
| --- | --- | --- |
| D0 | 2.246076541249e-01 | 2.975552970676e-06 |
| R_norm | -1.079999877620e+01 | 5.000971716021e-09 |
| R_pole | -2.546642778538e-03 | 3.279572176983e-08 |
| P2 | 3.773260192248e-08 | 2.026711073392e-09 |
| P4 | 1.291491506575e-07 | 2.416546895602e-09 |
Recorded B2b cross-engine max relative spread (reported, not hidden): `4.368649426740e-02`.

## Dual Engine

B2c's new assembly was checked two ways: B1 `lane_extract` plus `observable_residuals`, and an independent SymPy expansion of the squared response-amplitude series through `x=omega^2` order two. Root finding, when bracketed, uses Brent on the primary model and independent bisection on the secondary model.

| tau | kind | passed | max_abs | max_rel |
| --- | --- | --- | --- | --- |
| 2.912500000000e-02 | resolved_deepest | true | 4.336808689942e-19 | 7.015095236760e-16 |
| 2.925000000000e-02 | resolved_deep | true | 1.776356839400e-15 | 8.019161492366e-16 |
| 3.000000000000e-02 | floor_resolved | true | 4.336808689942e-19 | 6.087010429971e-16 |
| 1.000000000000e+00 | seed_validated | true | 1.776356839400e-15 | 3.038514360124e-15 |

## Gates And Negative Controls

| gate | failed_as_expected | wrong_answer |
| --- | --- | --- |
| frozen_vs_resolved_injected_drift | true | frozen-background root finder that ignores tau-dependent {B,Z,N} drift |
| confirmation_gate_mislocated_tau | true | mis-located tau accepted without a fresh |R_norm| gate |
| stable_side_filter | true | unstable-side D0<=0 root accepted |

## Reproduce

Commands used by this route, each intended to stay under `timeout 600`:

```bash
env PYTHONPATH=software/stage1_solver/src python -m stage1_solver.patha_b2c_calibration seed-tau1
timeout 600 env PYTHONPATH=software/stage1_solver/src python -m stage1_solver.patha_b2c_calibration probe-background --tau 0.03 --background-grid 16,16
timeout 600 env PYTHONPATH=software/stage1_solver/src python -m stage1_solver.patha_b2c_calibration probe-background --tau 0.0295 --warm-start-background software/stage1_solver/runs/patha_b2c_calibration/backgrounds/patha_b2c_background_tau_0p03.json --warm-start-final-only --warm-start-wall-predictor
timeout 600 env PYTHONPATH=software/stage1_solver/src python -m stage1_solver.patha_b2c_calibration probe-background --tau 0.02925 --warm-start-background software/stage1_solver/runs/patha_b2c_calibration/backgrounds/patha_b2c_background_tau_0p0295.json --warm-start-final-only --warm-start-wall-predictor
timeout 600 env PYTHONPATH=software/stage1_solver/src python -m stage1_solver.patha_b2c_calibration probe-background --tau 0.029125 --warm-start-background software/stage1_solver/runs/patha_b2c_calibration/backgrounds/patha_b2c_background_tau_0p02925.json --warm-start-final-only --warm-start-wall-predictor
timeout 600 env PYTHONPATH=software/stage1_solver/src python -m stage1_solver.patha_b2c_calibration probe-background --tau 0.0290625 --warm-start-background software/stage1_solver/runs/patha_b2c_calibration/backgrounds/patha_b2c_background_tau_0p029125.json --warm-start-final-only --warm-start-wall-predictor
timeout 600 env PYTHONPATH=software/stage1_solver/src python -m stage1_solver.patha_b2c_calibration evaluate --tau 0.029125 --kind resolved_deepest --background-bundle software/stage1_solver/runs/patha_b2c_calibration/backgrounds/patha_b2c_background_tau_0p029125.json
env PYTHONPATH=software/stage1_solver/src python -m stage1_solver.patha_b2c_calibration calibrate
timeout 600 env PYTHONPATH=software/stage1_solver/src pytest -q software/stage1_solver/tests/test_patha_b2c_calibration.py
```

