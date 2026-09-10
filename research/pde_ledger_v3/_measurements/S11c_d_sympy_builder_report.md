# S11c-d checkpoint: paused for sheet continuation

2026-09-10. **Stopped at the user's requested change-of-approach boundary.** The repaired inertia supports propagating transverse modes at the explicit preflight input. A subsequent check found a separate defect in the prototype's complex-momentum sheet labels. Sheet continuation must be repaired before completing the bulk current or scattering. The four-case regeneration was interrupted; the published main `.out` and its annex object are unchanged. No export or commit was made.

All spec §1 inputs remain SUPPLIED and unfalsifiable. The S11b/b shear-normalization discrepancy and c2 cross-engine operand debt remain open. This finding concerns the d continuation rule and does not adjudicate upstream response/face signs. No review leg, comparator, Wolfram engine, or downstream stage was launched.

## Finding and next step

At engine line 1951, `PHYSICAL_BULK_SHEET` is assigned by a half-plane test on the radical. At the reference `LAB_HELD × RHO4_CONSTANT` input, candidate 8 has `k_n = 0.607303653886 + 0.008491689257 i`. The engine labels `q = 0.080656407478 − 6.393830415310 i` physical. Continuation from real `k_n` on the outgoing branch instead gives its opposite, `q = −0.080656407478 + 6.393830415310 i`.

The probe derives the radical relation from the cached **computed reduced pencil**, then transports its root along the straight momentum path. Refinements at 16, 64, and 256 steps agree. The opposite-root residual is `1.39e-17`; the largest radical-equation residual is `1.43e-14`; the minimum distance to a computed branch point is `0.636783449147` in inverse-length reference units. The engine's `q` is in inverse-time reference units and equals `c_s0` times the bulk normal wavenumber.

This is a bounded counterexample, not a generic sheet construction. Small equation and normalization residuals can occur on either sheet. **Do not use the current complex-root sheet flags to declare channels or classify poles.** The real-momentum transverse preflight is unaffected.

Next: implement continuation from a stated physical/outgoing reference, retaining branch points, paths, domain failures and unresolved cases; rerun the end-channel census; then resume bulk-current/scattering construction. This changes the immediate work order, so the repair was not started after the stop.

## Preserved additions

| Construction | Engine source anchor and coverage |
|---|---|
| Import/Fourier wiring | Existing three-parent fold and one-dimensional reduction preserved; inherited energy root added with lookup witness. |
| Explicit physical input | `ChannelInput`, line 1600: independent step/bump profiles, computed limits and derivative-integral/jump residuals, rational L/T/M input, physical homotopy. Separate from PIT. |
| Finite end solves | `solve_input`, line 1841: full reference/end quotient pencils, candidate roots and local jets. Complex sheet labels have the defect above; bulk continuous spectrum remains unconstructed. |
| Frequency normalization | `solve_sample`, line 1858: full nullspace pairing including the radical chain rule, normalized left bases where defined, and independent nearby-frequency checks. |
| Closed-field check | `field_lift`, line 1400: differentiate the existing sector ansatz, then test lifted modes in the five-field equations. |
| Slab current | `UniformSlabCurrent`, line 1415: reduce inherited energy before variation; extract the material constraint from the zero-transfer mass row; integrate normal boundary work; emit retained current and discarded grades. |

The changed-end constraint's thickness weight depends on `W_0/W_bg`; its first-order expansion is computed from the mass row. The slab pairing on physical right fields remains a **partial** current operand. The nonlocal bulk contribution, full biorthogonal current, derivative identity and flux normalization are unfinished.

## Checks and artifact status

The integrated one-case input run completed with exit 0, empty stderr, 703.57 s wall time and 14,310,873 bytes. Its transverse roots are doubly degenerate: `k_n = ±0.785281265959` at reference/left and `±0.789514618822` at right. Measured transverse frequency-normalization residuals reach `5.11e-16`; nearby-frequency differences are below `1.7e-13`. Its preserved transcript, `scripts/out/S11c_d_channel_reentry_preflight.out`, predates the latest slab-current refinement.

The latest focused slab-current check at reference/right completed in 31.62 s with peak RSS 83,792 KiB. Variation, boundary, and retained-constraint residuals are zero; dimension constraints are empty. Both scripts compile. A complete integrated run of the latest source has **not** finished; the interrupted full capture was not published.

The diagnostic is `scripts/S11c_d_sheet_continuation_probe.py` (relation at line 42; continuation at line 68), with completed output `scripts/out/S11c_d_sheet_continuation_probe.out`. It finished with empty stderr and dimension constraints. `S11c_d_sheet_continuation_diagnosis.json` records the measured data, hashes, focused-current checks and unchanged main-transcript hash. Input is `S11c_d_channel_preflight_input.json`.

Reproduction: run the engine with `--case LAB_HELD__RHO4_CONSTANT --dev-symbol-cache /tmp/s11cd-sheet-cache --channel-input-file _measurements/S11c_d_channel_preflight_input.json`, capturing stdout in a fresh scratch file. Run the probe with that reference cache via `--symbol-cache`, that stdout via `--transcript`, the same JSON via `--input`, and `--candidate-index 8`. The current complex sheet labels remain under repair.

The inertia-repair manifest and `S11c_inertia_d_checks.json` remain historical records of commit `a74da30`, including the still-published main transcript. They do not describe the current uncommitted source. See [the repair record](S11c_inertia_repair_report.md).

## Remaining program

The 12 live TODOs are: all-carrier Fourier round trips; full end spectra; mixed-grade jets; generic sheet continuation; closed bulk current/flux normalization; complete two-ended scattering; poles/Riesz/overlap; survival; flux bookkeeping; weak coefficients; §5 controls; and export. Frequency pairing is implemented for finite algebraic pencils; general-domain coverage remains outstanding.

`scripts/S11c_d_exports.py` is absent because its required roots have not been computed. No empty pole set, empty S-matrix, or placeholder export substitutes for an unexecuted solve.

## Retained user-approved solver/export contract


1. Preserve `EdgeReduction`, the positional three-parent fold, and the exact
   direct-lookup manifest. All numerical assembly must consume the computed
   reduced rows, including the full nonlocal terms and full coupling vertex.
2. Accept explicit, independent dimensionless profile functions w(xi), m(xi),
   their derivatives, asymptotic limits and tail information. A selected smooth
   step with an independently adjustable localized modulus bump is a numerical
   instance; it does not replace the interface class. Store the profile formula
   and digest in every case record. The current preflight input is recorded in `S11c_d_channel_preflight_input.json`.
3. Require a complete parameter map in a declared L/T/M unit frame, real
   continuum frequency and tangential momentum, small contrast, and
   sigma_W = eta_bg W_0/L_W for evaluations on the physical homotopy. Retain
   independent eta/sigma grades in the symbolic calculation. Test actual
   reference/end channel availability before attempting flux normalization.
   Do not manufacture an incident channel by assigning a sector label.
4. Compute both-end modes, left/right normalization, the S11b-derived current,
   and the variable-profile matching problem. Re-expand the continuum response
   to the retained rectangle; do not present a finite-contrast numerical
   solution as a higher-order continuum prediction. Retain evanescent matching
   modes, channel degeneracies and domain failures explicitly.
5. Numerical pole searches have an explicit profile, parameter map, sheet,
   bounded search region and isolating contours. Evaluate the retained operator
   without the continuum re-expansion. Record boundary/quadrature resolution,
   domain size, precision, root residuals, contour-count evidence, and changes
   under refinement. A bounded search does not establish a global pole set.
   An unsuccessful or inconclusive search is unresolved, not an empty pole set.
   Compute residues/projectors and sheet/decay/width/closure tests only for
   actually resolved candidates; emit spectral overlap, not capture probability.
6. Separate transparent symbolic expressions from evaluated numerical records.
   Symbolic operator/continuum/weak-coefficient exports remain differentiable
   SymPy expressions, compacted with algebraic equivalence checks. Numerical
   mode and pole datasets retain their input bindings, domain, convergence
   evidence, dimensions and truncated-model status. They are not stand-ins for
   a generic symbolic profile-dependent root function. This is the export
   distinction motivating the user-approved contract; the downstream consumer
   will need to bind the appropriate representation explicitly.
7. Fingerprints summarize already constructed/evaluated objects. They do not
   evaluate nonlocal integrals or replace a spectral solve. Algebraic PIT and
   physical numerical evaluation are separate records. All completed roots use
   fresh lowerCamel write-keys and the existing bind-closure/minimal-delta guards.
8. Finish the one-case path, then implement controls and bookkeeping, then run
   all four cases once and write the complete export. No review legs,
   comparator, Wolfram engine, downstream stage, or commit belongs to this lane.
