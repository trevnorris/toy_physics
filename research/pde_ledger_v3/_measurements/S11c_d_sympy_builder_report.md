# S11c-d builder checkpoint: fixed-frequency sheet repair

2026-09-10. The original complex-momentum sheet error is repaired in d. A computed Fourier kernel connects the real-axis outgoing seed to the counterexample inside its absolute-convergence strip; this finding requires no upstream producer/export change. The repaired four-case run completed and its canonical transcript is published. This is an unreviewed builder checkpoint, not completion of the S11c-d program.

## Computed repair

[BulkSheetPath](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:1674) replaces the radical half-plane predicate with root transport from real normal momentum, at positive real frequency. It computes branch points, path clearance, two adaptive refinements, radical residuals and root-match differences. Paths intersecting a branch point or failing numerical resolution emit `UNRESOLVED`; their candidates remain recorded. Incoming/outgoing flags cannot accept unresolved sheet labels. Their frequency slopes remain diagnostics until the full current is constructed.

The [Fourier probe](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_fourier_sheet_probe.py) derives the spatial kernel from the computed reduced radical and actual c2 branch bindings. At the explicit input it obtains the strip `|Im(k_n)| < 0.2` in inverse-length reference units. For the original candidate, its refined kernel integral gives `q = −0.08065640747826197 + 6.393830415309923 i` in inverse-time reference units: the opposite of the old physical-sheet label, with opposite-root difference `5.53e-16 − 1.92e-16 i`. The repaired selector labels the old root false. This is a fixed-frequency continuation chart, not a construction of complex-frequency pole sheets or the full nonlocal resolvent.

## Checks and artifact status

The fresh pre-repair one-case input run completed in 801.55 s, with peak RSS 1,720,624 KiB, exit 0 and empty stderr. Its pinned pencil feeds the Fourier probe and cached checks. The independent ODE endpoint differs from the Fourier endpoint by about `1.14e-13`; its maximum sampled equation residual is separately `2.48e-10`. Two local frequency/momentum path-order checks differ by at most `1.26e-15`. Quadrature precision and cutoff are recorded; the quadrature error estimate excludes the truncated tail.

Cached controls cover 198 candidate records in nine packets (six PIT, three explicit input). All 170 resolved paths complete independent ODE transport, with maximum relative difference `4.98e-12`. All 18 deliberately intersecting branch-locus paths remain unresolved. These are candidate counts, not open-channel counts. The three diagnostic instruments exit 0 and emit empty dimensional constraints. All six engine/instrument Python files compile. The repaired four-case/PIT run completed in 2,495.66 s with peak RSS 1,721,064 KiB, exit 0 and empty stderr. Its 55,823,604-byte transcript contains all 12 reference/end symbols and 24 mode packets (528 candidate records), with no missing path records or mode metadata, no duplicate tags, unchanged source pins and three empty dimensional records. It retains 64 unresolved sheet labels. The explicit physical input is covered by the separate baseline and cached checks above.

Commands, hashes and source/cache linkage are in [the run record](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_sheet_repair_runs.json); scalar results, canonical diagnostic output paths and scope are in [the repair report](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_sheet_repair_report.md). The [output inventory](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_sheet_output_inventory.json) records the full-run checks. The directive’s two named `reduction/` checker scripts are absent; they were not executed. The old annex payloads remain preserved at checkpoint `55e298b`. Storage policy for the user-requested preservation checkpoint: the five published `.out` payloads go through DataLad/git-annex; scripts, reports and JSON stay in ordinary Git. The checkpoint conveys no review clearance. No review leg, comparator, Wolfram engine or downstream stage ran.

## Remaining program and inherited limits

The positional three-parent fold, Fourier reduction, full quotient-pencil candidates, frequency normalization, field lift and partial slab current are preserved. The slab current still lacks the nonlocal bulk contribution, complete biorthogonal current, derivative identity and flux normalization. All spec §1 inputs remain SUPPLIED and unfalsifiable here; the separate shear-normalization and c2 cross-engine operand/sign debts remain open.

The 12 live TODOs remain: all-carrier Fourier round trips; full end spectra; mixed-grade jets; generic sheet continuation; closed bulk current/flux normalization; complete two-ended scattering; poles/Riesz/overlap; survival; flux bookkeeping; weak coefficients; §5 controls; and own-row export. The generic continuation TODO retains complex-frequency paths, cut-bank/continuum treatment and full spectral coverage. `scripts/S11c_d_exports.py` is absent because its required roots are uncomputed; no placeholder export or empty pole set was substituted.

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
