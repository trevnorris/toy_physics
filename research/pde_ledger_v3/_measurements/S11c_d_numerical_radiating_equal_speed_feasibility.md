# Equal-speed rest-bulk benchmark: source-only feasibility

2026-09-30 local. User authorized resuming the existing central-balance-v2 process, reporting that central benchmark with its agreed controls, and stopping. The separate equal-speed question authorizes source/metadata assessment only. No material-advected case, neighboring frequency, new physical input, build, or equal-speed run is launched.

**Conclusion:** a parameter change is straightforward to specify, but a trustworthy exact-equal-speed result is not a cheap rebind of the completed numerical arrays or a supported configuration of the present boundary instrument. The obstacle is the incoming mode meeting the bulk branch locus. This is a specific scope limitation, not evidence that the physical solution diverges or cannot exist. The cost of treating the coincidence is not established from metadata, so no quick-run or day estimate is promised.

## Specify what is being matched

The bare S9/R4 coefficient ratio is 0.1 for mu_R=1, rho_br=1 and c_s0=10. The actual saved full LEFT transverse branch at omega=3 instead has phase speed 1.2247448713915874 in reference units, with k_normal=2.439262183530097 and |k_parallel|^2=1/20. At the full-contrast RIGHT end the saved phase speed is 1.2186666955535794. These are fixed-point readouts, not a proof of a global isotropic dispersion or the calibrated light identification.

Setting c_s0=1 would match the bare coefficient formula but would not, merely by doing so, match the actual saved full transverse branch. A bulk speed near 1.2247448714 is the reference-side modal matching target suggested by the saved data if the brane inputs are retained. Its exact source-bound value and the re-bound modes would need to be established; a rounded decimal is not an exact degeneracy. Alternatively changing brane stiffness/density while leaving c_s0=10 changes different physics. Neither is uniquely prescribed by the phrase “set the ratio to 1.” A single bulk speed need not match both contrast-dependent ends. For a prospective comparison, reference-side actual transverse matching is the clearest candidate, not an authorization to change it here.

## What changes at coincidence

For a fixed-point transverse phase speed c_T, the rest-bulk relation evaluated on that branch is

    q_on^2 = omega^2/c_s0^2 - (k_normal^2 + |k_parallel|^2)
           = omega^2(1/c_s0^2 - 1/c_T^2).

At c_s0=c_T it gives q_on=0: the incident normal momentum lies on a branch endpoint b=sqrt(omega^2/c_s0^2-|k_parallel|^2). This is mode/branch coincidence at nonzero incident momentum, not necessarily the separate b=0 threshold at which the Fourier radiating interval first appears. For reference-side matching suggested above, the branch endpoint moves from 0.2 toward the present incident momentum 2.4392621835. The equality itself does not establish a nonzero bulk drive, radiated power, or infinite response.

The current code intentionally excludes this case:

- `S11c_d_numerical_radiating_boundary.py:193–205` retains the saved material mapping and changes only frequency/contrast. The saved seeds, reference correspondence and continuation results therefore cannot be adopted unchanged under a new speed.
- The same source at lines 358–364 requires selected physical depth momentum to have a resolved positive real or imaginary part. q=0 fails that selection.
- `JointBulkSheetPath` in `scripts/S11c_d_mixing_scattering_sympy_audit.py:5038` derives radical transport by dividing by the q derivative of the wave relation. Its `trace` rejects a zero seed or a branch-locus intersection (`:5109–5111`). The existing joint-sheet proof cannot pass through exact coincidence by lowering a tolerance.
- `S11c_d_numerical_radiating_integration_v3.py:122–132` binds the existing material-frozen source records and explicitly requires b^2=1/25 at omega=3. Its sin/cosh quadrature can in principle reuse its design at a different positive b, but its old nodes, coefficients, endpoint evidence and Gaussian comparisons are no longer results for those new inputs.

Permeability can make some combinations of the raw impedance Z0=rho_m*omega/q and the closure factor 1/(1+lambda_A Z0/rho_m^2) finite at grazing. That does not by itself establish every end lift, current form, selected basis or forced source action at coincidence. A decoupled transverse mode may have zero acoustic drive; a regular limiting expression must be established from the actual source rather than multiplying an undefined bulk quantity by a reported zero. Some apparent singularities may therefore be removable, but this assessment has not computed those limits.

The smallest real missing ingredient is an applicable grazing treatment of the selected incoming/end states and their physical current/face maps, with its source limits and branch correspondence. If that is regular, the numerical integration/assembly machinery may carry over. If it is not, it is a boundary-method change. Approaching unity without reaching it avoids exact coincidence but answers a different question and needs its own conditioning/distance-to-threshold checks; no such run is proposed here.

## Reuse and cost boundary

Reusable design and original operands: original symbolic pencils/source records, row and field units, profile, finite basis, transformed quadrature machinery, full polarized current formula, finite solver and uniform/scaling/domain/regulator/refinement/sign-control structure. The source-generation code retains original operands; restoring and rebinding them would be a new guarded operation, not a producer replay. No broad symbolic reconstruction is implied.

Not reusable as new-speed numerical evidence: material-bound frequency records, completed all-candidate paths, selected end/current/face maps, radical/denominator/endpoint certificates, middle and Gaussian integrals, affected finite matrices and solve results. Actual source basis arrays independent of the changed parameter might be reusable after operand/setting identity checks; blanket array reuse is not justified. See `S11c_d_frequency_source.py:119–137,166–181` for the original-to-material-bound packet construction.

The existing machinery makes a future implementation smaller than starting over, but almost all physics-dependent numerical work would still change. Exact coincidence is presently a method stop, so metadata does not support calling this an inexpensive additional run. Nothing in this assessment reopens the symbolic omega=1 track or authorizes a new campaign.

## Interpretation and current stop

Current result label: **rest-bulk (v_bulk_normal_0=0), c_gamma/c_s approximately 0.12 for the selected transverse branch, LAB_HELD/RHO4_CONSTANT development benchmark — not the calibrated, draining model**. Keep the fixed tangential direction, fixed profile, finite boundary/source-domain and unresolved retained-order physical-loss interpretation visible. The bare S9 ratio remains separately 0.1.

Leakage in the present benchmark would be evidence about this benchmark. No source result proves that its value must increase as the ratio approaches one: coupling, closure, phase matching and interference also change. Conversely, a deficit below the achieved numerical envelope here is not a calibrated-model no-loss result. A stationary profile and the absence of bulk flow remain independent assumptions.

The same PID4097233 was resumed with SIGCONT after checking start ticks85113744, exact command/cgroup, 2GiB/zero-swap/task containment, systemd infinity/Restart=no and the existing waiting completion hook. No source/gate/helper was edited and no new scientific process was launched. See `S11c_d_numerical_radiating_flow_calibration_resume.json`.

After this central job completes, inspect its actual numerical controls, integrity/resources and outstanding review qualifications, report the result together with this assessment, and STOP. Do not launch omega2.3/4, MATERIAL_ADVECTED, an equal-speed run or another validator automatically. Preserve partial work if it fails. The exact prior review verdicts remain as recorded; this source-only assessment does not clear them.
