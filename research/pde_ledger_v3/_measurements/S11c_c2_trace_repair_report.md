# S11c-c2 pressure trace repair checkpoint

The user authorized the c2 repair identified at `34a0b38b`. The native
`reference_pressure_kernels` construction now extracts the affine trace from
inherited `face_shift`, derives its Fourier multiplication kernel from the
existing profile definition, and solves the two/three-leg trace operator for
reference pressure. `build_face` differentiates that reference continuation
for the normal jet and solves the same inherited affine trace for the pressure
slot before slab closure. The c1 pressure remains the physical-face operand.

The full retained eta/sigma rectangle includes the ordered product of the
height map with the first-shape pressure response. Output, input and middle
momenta remain distinct. No observed d residual is subtracted from a result,
and no pressure jet is set to zero. The b/c1 exports, d native reduction,
physics authorities and S10/Lean work are unchanged.

## Completed trace checks

All four anchoring/density cases and both faces were evaluated symbolically.
The [map inventory](S11c_c2_trace_repair_maps_checkpoint.json) records **208
exact-zero residual scalars**, 456 computed objects/fresh keys and 872 checked
metadata paths. This includes the actual inherited source-row pressure/jet
coefficient factorizations and the two/three-leg trace compositions. The
matrix inverse compositions are arithmetic consistency checks; the source
coefficient comparisons use the actual imported slab/traction/closure rows.
These results retain their original rational domains and establish neither
global spectral coverage nor constant rank on exceptional loci.

The published [map transcript](../scripts/out/S11c_c2_trace_repair_maps.out)
is 1,575,712 bytes, SHA256
`c224b35e1376a4ca708e86c512ee8cf929c267980a1db0fbfd12b25b3afff57a`.
The run took 52.61 seconds with 182,304 KiB peak RSS; input-cache creation is
separate. Frozen sources preserve the exact producer used for that result.

## One-case closure and replay

The complete LAB_HELD/RHO4_CONSTANT c2 model and native traction/power operands
have been computed and saved under `/tmp/s11c-trace-repair-20260913/case_one`.
The focused transcript first stopped on an explicit inverse wave-amplitude
factor in the power metadata. The diagnostic now handles exact Laurent
monomials, retains the reciprocal/excluded-zero operands, and does not Taylor
replace an unexpanded rational denominator. Model replay restores dimensions
through native declarations and recomputes the open coupling kernel against
its saved operand. Partial attempts remain separate from complete publications.

The completed one-case diagnostic has **59 exact-zero residual scalars**,
including the canonical power residual and six open-kernel replay entries.
Its 124 objects, 259 metadata paths and 124 keys were reconstructed by the
validator. Replay took 147.43 seconds and 214,500 KiB peak RSS, excluding the
original model/power construction. The [case transcript](../scripts/out/S11c_c2_trace_repair_case.out)
has 2,095,956 bytes and SHA256
`7d518735828d2913a6a0e93e71fdb0a5d6268a50c778c2d45e001f75df464142`.
See the [case inventory](S11c_c2_trace_repair_case_checkpoint.json).

## Full c2 regeneration

The four-case native producer completed with exit zero and empty stderr in
4,461.77 seconds, at 2,649,512 KiB peak RSS. Its 150 tags include the complete
controls. Export verification records 44 exact-zero expression differences
and 387 passing structural/metadata/reciprocal checks. The delta still has 70
rows: only `s11cc2ClosedSlabOperator` and `s11cc2ClosedCouplingKernel` changed;
no keys were added or removed. Direct lookups still equal the import manifest.
The [stage inventory](S11c_c2_trace_repair_stage_inventory.json) found no
execution or provenance failures.

The generated export is 30,672,491 bytes, SHA256
`97160450e8b1d9aa654cebb2ea0d2074df24186f337fad2ee850bb5b200e7579`.
The complete native transcript is 530,910,368 bytes, SHA256
`9e3a2da7d56cb7124e7b96ea47d0a459c48aecf5224d8d8431804316e9c6c289`.
The native transcript was atomically published to the standard c2 output path.
All four canonical power residuals are zero; all 24 diagnostic records have
finite arithmetic outside the explicitly unbounded integration endpoints.
The power check took 130.82 seconds and 661,700 KiB peak RSS. See the
[power inventory](S11c_c2_trace_repair_power_checkpoint.json) and
[power transcript](../scripts/out/S11c_c2_trace_repair_power.out) (273,533 bytes,
SHA256 `7aee3a1d8fe8fea2068f9803119ff1f1f56c8d314b92e44aaee69d30dc8bdd5c`). These c2 upstream transcripts retain the established full-object
format; new focused diagnostics and d heavy-object packets remain bounded.

The development-case reference/LEFT/RIGHT reduced operators have been rebuilt
through the unchanged native d reduction in 418.63 seconds, at 2,448,116 KiB
peak RSS, with exit zero and no unresolved dimension constraints. The scoped
producer computes no spectrum. The fresh RIGHT source construction completed
in 822.16 seconds at 153,048 KiB peak RSS, with exit zero and empty stderr.
All 33 top-level source residual scalars, including the former ten first-contrast
mass/mechanical discrepancies, are zero; all 127 exact arithmetic identities
are zero. Validation reconstructs 85 source objects, 1,651 cancellation objects,
368 native metadata paths and 3,484 tags. This uses the development case and
supplied endpoint limits while retaining all material symbols.

The [RIGHT transcript](../scripts/out/S11c_d_end_current_source_right_trace_repair.out)
is 8,789,042 bytes, SHA256
`97737027919eb8e217e0fd79df10481ce6862d7384f9b10932f2cdac4f3a8d1e`.
See the [source inventory](S11c_d_end_current_source_right_trace_repair_checkpoint.json).
The fresh LEFT run completed in 127.13 seconds at 123,724 KiB peak RSS, with
33 zero top-level source scalars and 82 zero arithmetic identities. Its 85
source objects agree with the preserved runtime baseline: 84 are structurally
identical and the remaining complete 52-entry parameter map agrees by symbolic
key. No algebraic difference remained. The [LEFT inventory](S11c_d_end_current_source_left_trace_repair_checkpoint.json)
records 1,066 cancellation objects, 368 native metadata paths and 2,314 tags.
The [LEFT transcript](../scripts/out/S11c_d_end_current_source_left_trace_repair.out)
is 2,625,497 bytes, SHA256
`57d28cc10b83167bf6b48ec70d0d600f1a7d9d5f6a1271a0803ad70ef6fcdd97`.

Both source joins are now established at the retained orders for the
development case. This clears the original source discrepancy; physical
current normalization still requires the subsequent two-frequency/end-mode
construction. Full d regeneration has completed under the repaired c2 input
in 10,396.23 seconds at 2,475,936 KiB peak RSS, with exit zero, empty stderr
and stable input/source hashes. All 12 end-symbol caches were saved. Its
83,972,570-byte transcript has SHA256
`f8b16e08920a9f5856bc26a4b2f37060673abbe87b9820e74014e7542d265fa3`.
All nine output inventories completed with exit zero and empty stderr.
The [full-run inventory](S11c_c2_trace_repair_d_full_checks.json) records 24
spectrum packets and 432 root/lift candidates: 336 scalar and 96 complete
two-dimensional spaces. Every packet has nine isolated radical roots of the
degree-11 polynomial, with zero degree/count residual. The sheet records retain
480 transported paths and 24 individually recorded branch-locus intersections.
There are 432 normal-momentum residues and 1,728 contours.

The threshold inventory retains 24 left and 24 right chain spaces, 96
individual chains, 336 bulk paths and 144 normal paths at the specified
targets. Earlier exceptional-census unresolved labels are followed by these
separate chain records; no global/complex-parameter coverage is inferred.
All 96 carrier roundtrips, the literal inverse/remainder/source-image/branch
residuals and the 222 integral/eight row reconstruction projections are zero.
The lossless codec reproduces 229,636 records and all 229,634 source-index
assignments from its 314,000,831-byte expanded view without a payload mismatch.
The established raw double-precision bank cancellation remains visible; the
worst bank inverse precision-refinement coefficient is 4.852e-25 in its
emitted unit frame. Numerical residuals and original domains remain in the
inventories.

The [main d transcript](../scripts/out/S11c_d_mixing_scattering_sympy_audit.out)
is atomically published. Its previous annex payload is unchanged; see the
[publication record](S11c_c2_trace_repair_d_publication.json). The
[scoped/full provenance joins](S11c_c2_trace_repair_d_source_joins.json) match
the weak and strong pencils, curl, weak units and energy basis at reference,
LEFT and RIGHT, with 229 shared dimension bindings agreeing at each end.
These are structural joins of the same native constructors, not an independent
physics derivation.
Both end-frequency packets have also been recomputed and validated: 36
candidates, 44 full basis directions and 90 exact-zero reconstruction/certificate
scalars. The largest algebraic coefficient-frame residual norm is 2.103e-13;
finite-frequency difference checks remain nonzero and refine under a halved
step. See the [frequency refresh](S11c_d_end_frequency_report.md). The new
reference-current source check completed in 87.34 seconds at 200,672 KiB
peak RSS, with exit zero and empty stderr. All 30 residual scalars and the
mechanical-orientation anchor residual are zero. Its 146 tags have complete
metadata, no nonfinite objects, no unknown dimensions and no dimension
constraints. The [current inventory](S11c_c2_trace_repair_d_current_checks.json)
and [publication checkpoint](S11c_c2_trace_repair_d_current_checkpoint.json)
pin the fresh [reference transcript](../scripts/out/S11c_d_nonlocal_current_reference_trace_repair.out):
1,495,327 bytes, SHA256
`34f5cd6bf8cf20c9a0106300a5b017a0b4de1b0759e56e54afd0d4ddf4236209`.
Older focused two-frequency/modal/adjoint packets retain their historical
provenance until the next current extension.

## Repair completion and next implementation

The full producer was run with
`_measurements/S11c_c2_trace_repair_run_stage.py --run-directory <fresh-dir> c2`.
The four-case power validation and c2 transcript publication are complete.
`S11c_c2_trace_repair_d_ends.py` has rebuilt the development-case reduced
reference/LEFT/RIGHT operators. Both symbolic source joins and LEFT baseline
regression are complete. The full d spectral regeneration has completed; its
output inventories, main publication and both end-frequency refreshes are
complete, as is the reference-current source refresh. This repair checkpoint
has no running computation or machine-capacity deferral. All complete outputs
are under scripts/out; partial attempts remain separately identified in scratch.

The next construction is the [both-end two-frequency physical current](S11c_d_both_end_current_plan.md),
including actual endpoint sources, full mode spaces, explicit adjoint maps and
outward flux orientation. The current native builders still restrict that
pairing to the reference state; extending their interface and operands is the
next implementation. Variable-profile matching follows. The ten broad d TODOs
remain, including complete scattering and section 3b profile-frequency poles.

The b/c1 exports, native d code and retained solver/export contract are unchanged.
The [final integrity inventory](S11c_c2_trace_repair_final_checks.json) records
12 published artifact/hash matches, stable producer sources, 14 parsed Python
files, 26 parsed JSON files and a clean tracked-source whitespace check.
No d export, review/comparator/Wolfram/downstream run, commit or push was made.
Generated outputs await DataLad/git-annex at a requested commit checkpoint.
