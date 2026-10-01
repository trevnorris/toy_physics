# Central omega=3 development benchmark: completed, numerical loss unresolved

**Scope:** strict rest bulk (`v_bulk_normal_0=0`), selected transverse phase-speed ratio approximately 0.12, LAB_HELD/RHO4_CONSTANT, saved tangents `(1/5,1/10)`, one fixed profile. This is not the calibrated, draining model. The bare S9 coefficient ratio is separately 0.1. No neighboring frequency or alternative anchoring was run.

All ten planned central finite solves finished. The result is **UNRESOLVED**, not a positive leakage measurement, a no-leakage result, or a physical upper bound. The prescribed sign control flags an apparent gain in one incident column at the baseline full contrast. Refinement and domain changes also move the small signed imbalances substantially. The user requested stopping after this benchmark; no retry, new validator, neighboring-frequency run, or ratio-1 run is launched.

## Actual signed current results

The quantity is `D_T = 1 - (reflected + transmitted transverse current)/incident transverse current`, with full polarized cross terms. Positive means a deficit; negative means apparent gain. Entries below are parts per million of incident current (ppm = 1e-6). The four columns are incident channels, not four physical cases; their detailed mode/current ordering remains in the saved end maps. Columns 1 and 4 stay within 0.000012 ppm throughout these full-contrast settings.

| Full-contrast setting | Column 2 deficit (ppm) | Column 3 deficit (ppm) |
| --- | ---: | ---: |
| Baseline | -1.773911 | +1.633609 |
| Refined quadrature/basis/cutoff | -1.034869 | +0.942457 |
| Larger source/domain interval | -0.608245 | +0.554788 |
| Larger interval, regulator halved | -0.608218 | +0.554792 |

The declared empirical resolution floor is **1 ppm**. No column satisfies the prescribed positive-resolution criterion (above three times its comparison envelope and within the sensitivity requirement). Baseline column 2 is below -1 ppm, and its refined value also slightly exceeds that negative threshold. The final two settings place both small imbalances within the floor. A sign flag is not evidence for a physical active/gain channel; these results do not identify its cause.

The unmodified worker status is `CENTRAL_FINITE_PILOT_UNRESOLVED_STOP_NEIGHBORS`. In this report, its legacy `NO_DEFICIT_RESOLVED` entries are displayed as **DEFICIT_NOT_NUMERICALLY_RESOLVED**, as requested by the independent review. The original JSON remains immutable. `balance_result.json` separately reports `D < -E_num` for every setting, retains all raw statuses, and preserves incident/reflected/transmitted currents, condition numbers and residuals.

## Controls and their limits

- **Uniform background:** all four numerical settings have maximum absolute deficit at most 1.229e-11 (0.00001229 ppm), far below the 1 ppm floor. This is a good numerical uniform control, not a subtraction of absolute regulator absorption.
- **Contrast scaling:** baseline column 2 at quarter/half/full contrast is -0.506472/-0.994062/-1.773911 ppm; column 3 is +0.492952/+0.944425/+1.633609 ppm. The prescribed report leaves scaling unresolved for every incident column. These numbers do not establish a positive contrast-squared leakage law.
- **Refinement:** the full-contrast changes in columns 2/3 are 0.739042/0.691152 ppm. **Domain:** the next changes are 0.426624/0.387669 ppm. They are substantial relative to the reported small signal.
- **Regulator:** the last changes are 2.753e-11/3.864e-12 of incident current. Small sensitivity between these two finite regulator values is not an absolute regulator-error bound.
- **Envelope:** the code uses adjacent-setting movements with a 1e-6 floor. Across the whole sequence, the total spreads in columns 2/3 are 1.165694/1.078821 ppm. Report both facts; do not promote the floor to a rigorous error bound or change the acceptance rule after seeing these values.
- **Linear solves and partition:** reported ranks are full (1575 baseline, 2365 refined systems); scaled residuals are below 5.0e-16. Condition numbers span approximately 1.05e6–5.58e6. The runtime independent LU/SVD comparison had to pass before each saved solve returned; its complete operands/checks are retained. The largest saved normalized mixed-current contribution is 1.870e-19, so the recorded partition control does not explain the ppm differences. Tiny algebraic residuals do not demonstrate continuum or boundary accuracy.
- **Matrix construction:** all 80 native rows were assembled in each setting, grouped 40/30/10 by one/two/three momentum variables. Published matrix/action residuals are below 4e-15; the largest Fourier-measure residual is about 3.2e-12. These check assembly consistency, not physical accuracy.

The Gaussian refinement and selected independent middle-leg comparison were consumed from their completed saved returns. They are not newly repeated here. The two fresh replacement reviewers both literally returned **NEEDS REVISION** for the original-tolerance route. Their requested actual Gaussian matrix-route join, momentum-zero check at the adopted order, and contained saved-operand row/trial checks remain unperformed. Source-dispatch tests and the already implemented single-leg orders do not replace those checks. The coarse branch was unused. There is no paired independent method/result clearance; no reliable loss interpretation is claimed while those findings remain open. Their literal reports and prior disposition are preserved.

## Completion and preservation

The coordinator, resource guard and normalization supervisor finished normally with empty scientific stderr and byte-identical stdout/checks. That operational success is separate from the numerical unresolved status. Inspection verified all 161 input hashes, 74 source snapshots, 617 opaque SQLite blobs, 53 complete journal operations (20 byte-identical restored returns and 33 new operations), 90 partial-matrix checkpoints and 320 row-matrix blobs. Launch/gate/worker/manifest hashes and the guard/supervisor command chain join. Complete operation returns, intermediate matrices, final result and all prior startup/failure evidence remain under the original scratch root. Inspection uses JSON/source/hash/opaque bytes only, with no unpickling or new scientific calculation. Detailed integrity/resource counts are in `S11c_d_numerical_radiating_balance_checkpoint.json`.

Worker wall time was 54,424.11 seconds (**15 h 7 min**), including the user-requested approximately 33 min 28 s suspension. Elapsed time excluding that pause was about **14 h 34 min**; this is not a CPU-time measurement. The ten solve calls themselves totaled 86.10 seconds. Most time went into matrix quadrature, including approximately 48.9 million nodes across the four three-momentum groups. The result database is **8.47 GB**. The cgroup reached its 2 GiB cap and recorded 100,419 memory.max events, with zero OOM events, kills or swap. Worker peak RSS was approximately 1.60 GiB; minimum sampled host-available memory was 15.19 GB. All no-deadline, zero-swap, process/thread/CPU and host-memory protections remained in place. These are not “zero cap events.”

## Ratio-1 handoff: smallest viable next decision, no launch

The source-only feasibility assessment remains applicable: the exact-equal-speed case is not a supported parameter toggle. Match the **actual full transverse mode** at the reference end, rather than setting the bare S9 coefficient ratio and assuming the retained dispersion follows. The full-contrast RIGHT end has a different saved phase speed, so a single constant bulk speed need not match both ends. Changing bulk speed versus changing brane stiffness/density are different physical choices; calibration must specify which is intended.

Before paying for another finite matrix assembly, the smallest useful prerequisite is a **uniform reference-mode grazing calculation**: at matched speed the incident mode has `q_out=0`, precisely where the current selected-depth and joint-sheet checks refuse the state. Establish actual finite source limits, acoustic drive, current normalization and branch correspondence there. A zero transverse face drive might remove some apparent singularities, but that must be shown on the actual source expressions. Do not remove the guards, insert a tiny nonzero q, or substitute a near-equal-speed result for equality.

If those limits support a well-defined incoming/end map, reuse the original symbolic source records and the numerical quadrature/solver/control machinery. Rebind the material-dependent sources and compute the changed boundary, current, endpoint, integral and matrix evidence; old numerical arrays are not results at the new speed. Start with the matched uniform control before a profile solve, and retain the same sign/scaling/refinement/domain/regulator requirements. The present unresolved ppm signal and open operator-route checks cannot be inherited as cleared validation.

The source-only assessment cannot price the grazing treatment reliably. The measured central calculation supplies a realistic reference for the numerical part: roughly 14.6 hours excluding suspension on one CPU for this ten-solve suite, mostly quadrature. A speed change can alter cost, conditioning and singular behavior; that number is not a ratio-1 prediction. No multi-day method build, review submission or scientific run is authorized by this handoff. A successful rest-bulk ratio-1 benchmark would still omit drain flow and would not settle the calibrated draining model.

All prior accepted science, failures, review debt and incident history remain. Full Green/FORM/A11/A12, physical loss and analog-light calibration remain open. Leave Lean, S11_lean and the protected builder suffix unchanged. **Stop here for the user's next decision.**
