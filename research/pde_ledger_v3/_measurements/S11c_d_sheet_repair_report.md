# S11c-d: Fourier connection and fixed-frequency sheet repair

Historical record of checkpoint `f9e28f5`. The canonical main transcript has
since been regenerated for the
rectangular mode-jet checkpoint (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11c_d_rectangular_mode_jet_report.md`).
The sizes, hashes and verification results below describe the preserved
sheet-repair run; its run record and inventory remain unchanged.

2026-09-10. **The original counterexample is resolved within the computed Fourier domain.** The repair is in d; no upstream producer/export was changed. The complete d program remains unfinished. Four-case regeneration completed with exit 0, and the checked main transcript is published. This builder record and the instruments are unreviewed.

## Construction that resolves the scope question

The Fourier probe (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/scripts/S11c_d_fourier_sheet_probe.py`) consumes a newly generated reference pencil, its completed one-case transcript, the explicit physical input, and the actual c2 branch bindings captured by the scope lookup. It checks the run's source/input/cache/transcript hashes. The pre-repair engine run completed with exit 0, empty stderr, 801.55 seconds, and peak RSS 1,720,624 KiB.

From the computed reduced radical relation, the probe obtains the even quadratic in normal momentum and its roots. At this input, the branch points are `k_n = ±0.2 i` in inverse-length reference units. It scales momentum by the computed `a = 0.2`, producing the dimensionless radical square `1+p²`. The real-axis seed is evaluated from the exported c2 branch definition, with the conversion to the d radical obtained from the computed quadratic coefficient.

The probe then evaluates a Gamma integral representation of the normalized inverse radical, the Gaussian Fourier integral, and the resulting spatial kernel. SymPy returns that kernel as a Meijer-G expression. Integrating the kernel back gives a literal zero squared-radical residual. Computing the absolute-transform convergence condition gives `−1 < Im(p) < 1`, hence the open strip `|Im(k_n)| < 0.2`. The measured counterexample lies inside this strip. Its Fourier continuation is therefore fixed by this kernel and seed; choosing an unrelated radical half-plane is unnecessary.

Independent numerical integration of that computed spatial kernel, inverted and scaled by the exported seed, returns

`q = −0.08065640747826197 + 6.393830415309923 i`

in inverse-time reference units. The refined kernel calculation's sum with the old engine candidate is approximately `5.53e-16 − 1.92e-16 i`; the candidate was on the opposite branch. A separate implicit-ODE transport agrees at the endpoint to about `1.14e-13`. Its maximum sampled equation residual is `2.48e-10`, which is reported separately from the endpoint agreement. Increasing quadrature precision and cutoff changes the computed result toward the recorded endpoint; the quadrature routine's error estimate excludes the truncated tail and is not reported as a bound on the entire calculation.

Two small `(omega,k_n)` rectangle tests, with `omega = 1 ± 0.01 i`, compare frequency-first and momentum-first transport. Their endpoint differences are approximately `9.30e-16` and `1.26e-15`. These are local compatibility checks with S11b's real-momentum frequency continuation, not a global pole-sheet construction.

Thus no additional physical premise is needed to decide this original counterexample. This conclusion concerns the bulk radical carrier. It does not certify the full nonlocal resolvent, its other poles, or the unresolved upstream cross-engine response/sign debts.

## Repair and checks

[BulkSheetPath](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:1674) replaces the complex-root half-plane predicate. At the currently supported positive real frequency, it attaches a vertical path from `Re(k_n)` to the target. It computes branch points and path clearance, transports the real-axis seed with two adaptive step fractions, and emits the endpoints, residuals, and root-match data. Steps are limited by the distance to the computed branch points. A branch-point intersection or unresolved numerical match produces `UNRESOLVED`.

The flag records membership in this explicit fixed-frequency continuation chart. It is not a frequency-pole classification. Candidates on the other branch or outside the resolved path domain remain in the output. An unresolved label cannot pass the incoming/outgoing flag gate. Frequency slopes remain diagnostics; the unfinished bulk current is still required for flux normalization.

The cached path controls (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/scripts/S11c_d_sheet_path_controls.py`) recomputed labels for 198 candidates in nine one-case packets: six algebraic PIT packets and three explicit-input packets. Differential-equation transport completed for every resolved path. Maximum transport difference was `6.04e-12`, or `4.98e-12` relative. Label transitions were: 34 false→true, 28 true→false, and 28 becoming unresolved. These counts describe candidate records, not physical open-channel counts. All 18 deliberate branch-locus path probes were unresolved.

The build directive also requests `reduction/derived_or_declared.py` and `reduction/engine_output_checks.py`; neither file is present in this repository. The recorded compilation, import guards, emitted dimensional constraints and transcript inventory are the available builder checks, not execution of those missing tools.

The updated original probe (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/scripts/S11c_d_sheet_continuation_probe.py`) calls the repaired selector on the original counterexample. It emits `selectedPhysicalSheet = false`, a relative opposite-branch difference of `1.46e-16`, and empty dimensional constraints. The Fourier and cached-control instruments also emit empty dimensional constraints. These cached checks use the pinned pre-repair pencil as a construction operand; they are not presented as a fresh integrated execution of the repaired engine.

## Artifacts and remaining scope

Exact commands, run manifests, payload hashes, and exit/stderr records are in the run record (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11c_d_sheet_repair_runs.json`). Completed stdout is published at its expected paths:

- Main four-case/PIT [transcript](../scripts/out/S11c_d_mixing_scattering_sympy_audit.out): 55,823,604 bytes, from the repaired engine.
- Channel preflight (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/scripts/out/S11c_d_channel_reentry_preflight.out`): 15,882,954 bytes, from the freshly pinned **pre-repair** engine.
- Fourier probe (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/scripts/out/S11c_d_fourier_sheet_probe.out`): 36,365 bytes.
- Original counterexample probe (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/scripts/out/S11c_d_sheet_continuation_probe.out`): 13,844 bytes.
- Cached path controls (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/scripts/out/S11c_d_sheet_path_controls.out`): 313,768 bytes.

The full run used the engine's default four-case/PIT path; the explicit physical input is covered separately above. It completed in 2,495.66 seconds with peak RSS 1,721,064 KiB and empty stderr. The literal output inventory (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11c_d_sheet_output_inventory.json`) records 9,894 distinct tags, all 12 reference/end symbols and 24 mode packets, unchanged emitted source pins, one completion marker, and three empty dimensional records. All 528 candidates carry branch-path records and matching packet metadata: 224 labels are true, 240 false and 64 unresolved. These are algebraic PIT candidate records, not an open-channel census. The inventory exited 0 with empty stderr.

Publication replaced directory entries without writing through annex symlinks; the old annex content remains preserved by checkpoint `55e298b`. The previous main annex object's hash was checked again after publication and is unchanged. The user subsequently requested this repair checkpoint. Its storage policy sends the five `.out` payloads through DataLad/git-annex and keeps scripts, reports and JSON in ordinary Git. This preservation save does not confer review clearance; no push was requested.

Replaying the Fourier probe requires the pre-repair engine source pinned in its manifest (checkpoint `55e298b`), together with the current Fourier instrument and the regenerated baseline cache/transcript. The probe deliberately rejects a changed producer source. The cached selector checks then use the repaired engine with that baseline pencil. Scratch caches are reproducible intermediates, not committed artifacts; the run record preserves their hashes and producer commands.

The full Phase E verification matrix in the original plan is not complete: both frequency signs, global cut/winding paths, threshold families, boundary mutations and complete flux checks still require construction. The results above establish the bounded fixed-positive-frequency repair.

The generic sheet TODO remains: complex-frequency pole paths, cut-bank/continuum treatment, and full spectral coverage have not been completed by this repair. The other outstanding spectrum, current, scattering, pole, bookkeeping, control, and export constructions remain. No placeholder d export was created. All spec §1 operands remain supplied and unfalsifiable in this build; the separate shear-normalization and c2 cross-engine debts remain open. No review legs, comparator, Wolfram engine, or downstream stage ran.
