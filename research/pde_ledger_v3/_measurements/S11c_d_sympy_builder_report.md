# S11c-d builder checkpoint: regular mixed mode jets

2026-09-10. The engine now computes regular end-mode jets through the retained
`(eta,sigma_W)` rectangle. The full four-case run and output inventory completed,
and the regenerated transcript is published at its canonical path. This is an
unreviewed runnable checkpoint; the complete S11c-d engine and export remain
unfinished.

## Computed extension

[RectangularModeJets](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:1762)
derives the implicit radical chain rule and pencil Taylor coefficients from
the computed reduced operands. Its
[augmented solve](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:1826)
computes right and adjoint-left invariant-pair coefficients at grades `10`,
`01`, and `11`, retaining matrix normal momenta for degenerate clusters. The
[radical and classifier construction](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:1869)
computes the radical jet, overlap-inverse series and oblique classifier projector.
These are regular cluster jets, not Riesz data. Singular Jacobians retain explicit
status. The [emission](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:2359)
prints literal equation/gauge residuals and fingerprints the heavy objects,
with grade and restored dimension metadata. The prior first-grade route remains
computed for operand/difference checks.

## Checks and artifacts

Three sequential focused checks on one case passed with empty stderr and empty
dimensional constraints. The production check uses a second PIT sample. The
original-coordinate check computes absent sigma dependence at that constant end.
A mathematical parameter pullback exercises nonzero mixed coefficients; its
right/left coefficient residual maxima are `1.99967e-14` and `1.09836e-13` in the
declared numerical coefficient frame. Independent direct-pencil mixed residuals
fall `0.593519 → 0.297380 → 0.148845` when the step is halved twice. This is not a
physical profile ablation or a Section 5 control.

The full run took `3728.55` seconds, peak RSS `1721664` KiB, exit 0 and empty
stderr. Its
[78,707,339-byte transcript](/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_d_mixing_scattering_sympy_audit.out)
contains all 12 reference/end symbols and 24 mode packets, with 528 defined
rectangular jets (336 scalar and 192 two-dimensional nullspaces). All three
dimensional records are empty; tag, source-pin and residual-metadata coverage
checks passed. These are PIT candidates, not a physical open-channel census.
All 12 stored symbol payloads and all 528 root/nullity/sheet records match the
repaired checkpoint literally, including its 64 unresolved sheet labels.

Full-run grade-10 residual maxima remain recorded by physical dimension: for
example, `9.35918e-14` at `[L,T,M]=[-1,-2,1]` on the right and `5.85929e-14` at
`[0,0,0]` on the left. All grade-01 and grade-11 equation residuals are literal
zeros. The nonzero pullback check above exercises the mixed computation.
Detailed scalar results and per-object scope are in the
[construction report](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_rectangular_mode_jet_report.md),
with exact commands, hashes and producer/cache linkage in the
[run record](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_rectangular_mode_jet_runs.json)
and literal full-run counts in the
[inventory](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_mode_jet_output_inventory.json).

The engine and both new instruments compile. The directive-named
`reduction/derived_or_declared.py` and `reduction/engine_output_checks.py`
are absent and were not executed. Publication preserved the
previous annex payload; its hash was checked after replacing the main path.
The user subsequently requested a local preservation checkpoint. Its four
new/regenerated `.out` payloads (81,982,175 bytes) use DataLad/git-annex;
sources, reports and JSON remain ordinary Git. This save conveys no review
clearance. No review leg, comparator, Wolfram engine or downstream stage ran.

## Remaining program and inherited limits

`MIXED_GRADE_MODE_JETS` is removed from the live TODO; 11 constructions remain:
all-carrier inverse-Fourier round trips; full end-spectrum coverage; generic
sheet continuation; closed nonlocal bulk current and flux normalization;
complete two-ended scattering; poles/Riesz/overlap; survival; flux bookkeeping;
weak coefficients; Section 5 controls; and own-row export. The regular jets do
not settle singular clusters, generic individual branches or full spectral
coverage. `scripts/S11c_d_exports.py` remains absent because its required roots
are uncomputed; no placeholder export or empty pole set was supplied.

The positional parent fold, Fourier reduction, field lift, frequency
normalization and partial slab current are preserved. The full nonlocal current
and derivative/flux identities remain unfinished. Section 1 inputs remain
SUPPLIED and unfalsifiable here; the separate shear-normalization and c2
cross-engine operand/sign debts remain open. The prior fixed-frequency repair
is recorded in the
[historical repair report](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_sheet_repair_report.md).

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
