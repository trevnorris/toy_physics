# W3/W4 round 4 report

Scope: no engine, step record, committed `.out`, or `checks_S*.yaml` changed; the prohibited S9 result record was not read.

## R1 shifted-registry printer
Driver: `/var/projects/toy_physics/research/pde_ledger_v3/reduction/w4_shifted_registry_printer.py`; intentional exit 1. Relevant literal stdout:
```text
S11_NUMERIC_REFERENCE_COMPARISON S11-py Q.brane.B_comp: engine-operand=[2 - Q.brane.D_brane,-2,1]; declared-reference-operand=[0,-2,1]; residual=[2 - Q.brane.D_brane,0,0]; declared-bindings=Q.brane.D_brane=3; residual-at-declared-bindings=[-1,0,0]
S11_NUMERIC_REFERENCE_COMPARISON S11-py Q.brane.mu_R: engine-operand=[2 - Q.brane.D_brane,-2,1]; declared-reference-operand=[0,-2,1]; residual=[2 - Q.brane.D_brane,0,0]; declared-bindings=Q.brane.D_brane=3; residual-at-declared-bindings=[-1,0,0]
S11_NUMERIC_REFERENCE_COMPARISON S11-py Q.brane.rho_br: engine-operand=[-Q.brane.D_brane,0,1]; declared-reference-operand=[-2,0,1]; residual=[2 - Q.brane.D_brane,0,0]; declared-bindings=Q.brane.D_brane=3; residual-at-declared-bindings=[-1,0,0]
ENGINE_DIMENSION_PIN: FAIL
```

## R2 weaker implementations
Runner: `/var/projects/toy_physics/research/pde_ledger_v3/reduction/w4_pin_completeness_runs.py`; wrapper exit 0. Literal stdout:
```text
PIN_WEAKER_TEST_GUARD drop-covered-quantity-check: FAIL_AS_REQUIRED
PIN_WEAKER_TEST_GUARD drop-population-count-check: FAIL_AS_REQUIRED
PIN_WEAKER_TEST_GUARD drop-errors-check: FAIL_AS_REQUIRED
PIN_WEAKER_TEST_GUARD drop-unmapped-symbol-raise: FAIL_AS_REQUIRED
W4_PIN_COMPLETENESS_WEAKER_IMPLEMENTATIONS: PASS
```
Each block also contains pytest's failed assertion: empty configured population, S11-only subset covering all required QIDs, nonempty error set, and missing unmapped-symbol exception, respectively.

## Documentation/method changes
README now records the pin's fixed-transcript scope, common-mode blind class, S11-only `B_comp` coverage, and the `mu_R`/`B_comp` attribution limit. The historical-spelling grep was removed from `w3_duplicate_pin_ablation.py`; its computed baseline/mutation exit flip remains.

## Regression evidence
- Loader, gate, acceptance, able-to-fail aggregate/cases/demos, pin/law checks, and mutation runners produced their guarded outcomes.
- Committed checks: S9 `exit=2, agree=12/12, registry=1/1`; S10 `exit=2, agree=517, coverage=545/690, registry=1/1`; established operational findings unchanged.
- Suites: current `138 passed in 144.92s`; pristine `HEAD` archive `87 passed in 47.18s`.
- Fence discovered S10-Python and S11-Python. Both exited 0, stderr was empty, and stdout was byte-identical: `937aed06…4886`, `5ed934e5…b009`. No output diff was emitted.
