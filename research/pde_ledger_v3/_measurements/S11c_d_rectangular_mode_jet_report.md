# S11c-d regular rectangular mode jets

2026-09-10. This construction extends the repaired S11c-d engine at checkpoint
`f9e28f5f28dcb4b66d954d883b2977075598ca76`. It computes regular end-mode
coefficients through the retained `(eta,sigma_W)` rectangle. The four-case
regeneration and transcript inventory completed with exit 0. The canonical
main transcript is published. This is an unreviewed, runnable builder checkpoint;
the complete S11c-d program and its export remain unfinished.

## Computation

[RectangularModeJets](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:1762)
derives the implicit bulk-root chain rule and pencil Taylor coefficients from
the actual reduced pencil and its radical relation. Its invariant-pair ansatz
uses a matrix normal momentum for each computed nullspace, preserving a
degenerate cluster. Coefficient convolution retains both multiplication orders
in the mixed term. The augmented equation and gauge Jacobian is constructed by
coefficient probes, then solved successively at grades `10`, `01`, and `11`.
Singular Jacobians retain an explicit domain status.

The [right/adjoint-left construction](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:1869)
also computes the radical jet, overlap-inverse series, and oblique classifier
projector. These are regular cluster jets and nullspace classifiers, not a
Riesz projector or a construction of generic individual branches. The
[production emission](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:2359)
records literal equation/gauge residuals, restores physical dimensions, and
fingerprints the mode, root, radical and projector objects. The previous
first-grade route is still computed; both operands and their differences are
fingerprinted. Numeric coefficients are evaluated once per fingerprint sample
when repeated; the sampled projection and digest definitions are unchanged.

The action, imported parents, Fourier reduction, candidate roots and repaired
fixed-frequency sheet construction were not edited. Frequency normalization
and the partial slab current remain as before; full bulk-current construction
is still required.

## Focused checks

The [check instrument](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_rectangular_mode_jet_check.py:32)
validates the preceding four-case transcript/cache hashes and parent-source
hashes. Its producer is checkpoint `f9e28f5`; its consumer is the new engine.
All three checks use only `RIGHT_LAB_HELD_RHO4_CONSTANT`, run sequentially,
exit 0 with empty stderr, and emit empty dimensional constraints. Engine and
instrument hashes are unchanged across the checks.

Completed focused stdout is available at the canonical
[production](/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_d_rectangular_mode_jets_production.out),
[original-coordinate](/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_d_rectangular_mode_jets_native.out),
and [pullback](/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_d_rectangular_mode_jets_pullback.out)
paths. Exact commands, hashes and producer linkage are in the
[run record](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_rectangular_mode_jet_runs.json).

| Check | Wall seconds | Peak RSS, KiB | Result |
|---|---:|---:|---|
| Production packet, PIT sample 2 | 101.9093 | 122652 | 22 regular cluster jets; 14 scalar and 8 two-dimensional nullspaces |
| Original coordinates, PIT sample 0 | 50.7461 | 122524 | 22 regular cluster jets; computed mixed root/radical coefficients zero |
| Parameter-coordinate pullback, PIT sample 0 | 56.7918 | 122628 | 22 regular cluster jets; nonzero mixed coefficients |

In the original coordinates the computed constant-end pencil has no sigma
dependence. The check therefore also pulls back that computed pencil by
`eta = alpha + beta`, `sigma = beta`, with fresh dimensionless coordinates.
This is a mathematical coordinate check, not a physical profile ablation or a
Section 5 control. It exercises the mixed construction using the same pencil.

The original-coordinate maximum right/left coefficient residuals in the
declared numerical coefficient frame are `2.22322e-15` and `9.95245e-15`.
For the pullback they are `1.99967e-14` and `1.09836e-13`.
The pullback's largest mixed normal-momentum and radical coefficients are
`10.1127088` in inverse-length units and `4.52645849` in inverse-time units.
Direct evaluations of the original, undifferentiated pencil at four parameter
corners give mixed residuals `0.593518735`, `0.297379811`, and `0.148845381`
at steps `1e-3`, `5e-4`, and `2.5e-4`, respectively, in the declared numerical
coefficient frame. Their reduction under refinement checks the local expansion;
it does not establish global spectral accuracy.
This pullback also does not test arbitrary independent noncommuting parameter
directions or singular clusters.

The production packet separately records residuals with restored units. Its
largest right grade-10 coefficients by `[L,T,M]` dimension are `7.10738e-15`
at `[-3,-2,1]`, `3.55271e-15` at `[-3,-1,1]`, and `7.32411e-15` at
`[-1,-2,1]`. Left maxima are `1.58882e-14` at `[-2,0,0]` and `1.91319e-14`
at `[0,0,0]`. Its grade-01 and grade-11 equation residuals are literal zeros.

## Four-case run and publication

The single full four-case/default-PIT run completed in `3728.5501` seconds
with peak RSS `1721664` KiB, exit 0 and empty stderr. Source hashes were unchanged
throughout. The
[main transcript](/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_d_mixing_scattering_sympy_audit.out)
is `78707339` bytes, SHA-256
`a5f512c5af0e671377469c705002fa99a41095aa0f304f1c6a163d60002422ef`.

The [inventory](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_mode_jet_output_inventory.json)
completed in `84.2827` seconds with peak RSS `166172` KiB, exit 0 and empty
stderr. It records `23766` distinct tags, all 12 reference/end symbols and
24 mode packets, no duplicate tags or coverage gaps, unchanged emitted source
pins, one completion marker and three empty dimensional records. All 528
candidate records have defined rectangular jets: 336 scalar nullspaces and
192 two-dimensional nullspaces. These are algebraic PIT candidate counts,
not an open-channel census or an explicit physical-input run.

All 12 stored pencil-symbol payloads and all 528 stored `(k,q,nullity,sheet)`
records are literally identical to the repaired checkpoint's transcript.
Sheet counts remain 224 true, 240 false and 64 unresolved. This addition does
not expand the fixed-frequency sheet chart or resolve those domain limits.

The full-run right grade-10 equation residual maxima are `1.50990e-14` at
`[L,T,M]=[-3,-2,1]`, `7.11430e-14` at `[-3,-1,1]`, and `9.35918e-14` at
`[-1,-2,1]`. Left maxima are `2.38323e-14` at `[-2,0,0]` and `5.85929e-14`
at `[0,0,0]`. Grade-01 and grade-11 equation residuals are literal zeros in
every dimension group. The largest numerical-frame fingerprint projections
of the first-grade route differences are `5.03993e-16` (right), `5.82768e-16`
(left), and `1.81139e-16` (root); the corresponding projector-idempotence
projection maximum is `7.66989e-15`. These projections summarize computed
tensors and are not tensor-norm bounds. Each group contains 528 objects,
with three PIT projections per object and matching metadata.

Publication replaced the main directory entry without writing through its
annex symlink. The previous annex object's SHA-256 was rechecked after
publication and is unchanged. The user subsequently requested a local
preservation checkpoint. Its three focused transcripts and new main transcript
(81,982,175 bytes total) use DataLad/git-annex; scripts, reports and JSON remain
ordinary Git material. The preservation save conveys no review clearance.
The historical sheet-repair run record and inventory are retained unchanged.
Scratch symbol caches are reproducible intermediates; their commands, hashes
and producer/consumer source pins are retained in the run record.

## Scope

This closes the regular mixed-grade mode-jet construction. Full end-spectrum
coverage and singular-domain treatment remain under the existing spectrum
TODO. Generic sheet continuation, the full nonlocal current and flux
normalization, scattering, poles/Riesz/overlap, survival, flux bookkeeping,
weak coefficients, Section 5 controls and the own-row export remain open,
along with all-carrier Fourier round trips. No placeholder export or pole set
is supplied. Section 1 inputs remain supplied and unfalsifiable in this build;
the separate shear-normalization and c2 cross-engine debts remain open.
No upstream repair, review leg, comparator, Wolfram engine or downstream stage
is part of this checkpoint. The two directive-named reduction
checker scripts are absent and have not been executed.
