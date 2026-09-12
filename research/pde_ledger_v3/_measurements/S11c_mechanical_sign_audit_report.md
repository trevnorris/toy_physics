# S11c mechanical-load sign audit

Checkpoint `03491aca` committed the nine current-diagnostic files: eight ordinary
Git files and one DataLad/git-annex payload. The subsequent bounded audit locates
the mechanical sign disagreement at S11c-b's assembly of prescribed external
virtual work. No producer physics, export, or authority was changed.

## Measured source-to-row trace

The [instrument](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_mechanical_sign_audit.py)
replays S11c-a's native virtual-work construction and S11c-b's native force
extraction on the actual consumed operands. It covers both anchorings and both
density representatives, both faces, all four mechanical components, nonuniform
profiles, and the live chemical functional derivative. Reconstruction checks
use symbolic rational algebra, with original reciprocal bases retained.

All four geometry replays and force extractions agree with the exports. The
actual assembled face increment equals the extracted external generalized force
`Q`. Independently, the stored-energy stiffness coefficient and a functional time
variation of the supplied kinetic energy both give the action-to-row multiplier
`−1`, in every case. Applying that same multiplier to external virtual work
requires the load `−Q` in the stored left-hand-side convention. The complete
assembled increment plus this action-normalized load is exactly zero in all four
cases. This is a relative load sign; reversing the entire mechanical row would
also reverse its already anchored conservative terms.

Of 104 scalar physical comparison residuals, 94 are zero. The ten nonzero entries
are confined to the assembled-minus-action-load comparisons: the thickness entry
for each LAB_HELD case and all four mechanical entries for each MATERIAL_ADVECTED
case. Each of these ten entries has three nonzero carrier-PIT evaluations. These
are algebraic fingerprints of computed expressions, not spectral samples or a
claim that the load is nonzero at every parameter value.

The earliest attachment is
[S11c-b's row assembly](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py:2991):
its extraction computes the external-work coefficients, and its mechanical rows
add them with a plus sign. S11c-a's negative traction/work convention agrees with
the supplied law. The actual inherited S11b thickness row retains stiffness
coefficient `W_0²` and Fourier inertia coefficient `−W_0² ω²`.
[c2's closure](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:531)
substitutes pressure and its normal jets into the inherited rows; it does not
correct their relative sign.

The source trace also identifies a limitation in
[c2's power control](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:771).
It defines kinetic-and-stored power by subtracting face power from assembled slab
power. Substitution into its residual cancels slab power, leaving face power
minus traction power. That checks the force pairing but does not independently
anchor its placement against the kinetic/stored variation. No c2 control or
review engine was run or imported as S11c-d's current in this audit.

## Evidence and next boundary

The [transcript](/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_mechanical_sign_audit.out)
has 144 unique write-keys and is 882,058 bytes. The successful four-case run took
295.30 seconds and 1,659,452 KiB peak child RSS, with exit zero and empty stderr.
The [checked inventory](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_mechanical_sign_inventory.json)
has no metadata gaps, nonfinite operands, unresolved dimensions, or source/payload
pin issues. Every physical leaf carries restored dimensions and perturbation
order. Exit zero records completed emission, including the nonzero comparisons.
The [run record](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_mechanical_sign_runs.json)
includes the native inventory, source snapshots, the development dimension fix,
and the interrupted attempt at global simplification of the nonzero residual.

The source-hash join connects this audit to the prior reduced reference current
preflight: its b/c1/c2 inputs and current engine are unchanged. That earlier
preflight's five mechanical coefficients remain negatives of its reconstructed
face load; its mass-row join remains zero. This audit does not extend that
preflight to other ends or supply new spectral/exceptional-locus coverage.

The [repair plan](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_mechanical_sign_repair_plan.md)
covers b's load assembly/provenance, an independent c2 power control, and serial
regeneration through c1/c2 into d. No current evidence calls for a physical repair
to S11b or a. c1's direct-input manifest excludes the affected mechanical and
coupling roots; unchanged physical c1 values are a prediction to verify after
repair, not an already completed comparison.

Work stops before the upstream repair under the user's standing instruction.
The main d transcript and retained eight-point contract are byte-identical;
all ten broad TODOs and the separate shear/cross-engine debts remain open.
The new audit is uncommitted. No S10/Lean edits, full upstream regeneration,
review leg, comparator, Wolfram run, downstream task, incomplete export, or push.
