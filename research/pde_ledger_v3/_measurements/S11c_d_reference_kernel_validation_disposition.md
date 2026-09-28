# Saved-output validation stop and proposed correction

**Resolved on 2026-09-27:** the user replied “I approve” to the prepared
correction. The single authorized corrected run passed all 17 structural
checks and three numerical residual/mutation checks in 4.394 supervisor
seconds, peaking at 95,674,368 bytes. Both runs and all unchanged constructor
artifacts remain preserved. The [acceptance report](S11c_d_reference_kernel_report.md)
and [checkpoint](S11c_d_reference_kernel_checkpoint.json) give the final limited
disposition. The account below preserves the original diagnosis and proposal.

The first guarded constructor completed in 19.314 seconds, with 109,998,080
bytes peak whole-job memory, zero swap and no memory-limit events. All 115
operation receipts and 497 result files were preserved. The three recorded
left/right inverse residuals are below `9e-80`, and the independent LU
differences are below `8e-81` in the saved reference-unit frame. Strict
constructor stderr is empty; checks/stdout, inputs and posthashes match.
These are verified runtime and recorded-result facts, not final acceptance.

The separate saved-output validator stopped after 5.073 guarded seconds
(82,608,128 peak bytes). Four structural checks failed: inverse scalar factors,
density scalar factors, and left/right residual expression layout. Thirteen
other structural checks passed, including the source binding, 50 saved zero
certificate scalars, determinant domain and unit orientation. The validator
stopped before its numerical residual rechecks. Both the failure and original
validator remain unchanged under the `validation` run and the original helper.

## Source-based diagnosis

The validator rebuilt comparison expressions with `evaluate=False` and then
required exact expression-tree equality to restored expressions. The installed
SymPy `core/basic.py` defines `Basic.__getnewargs__` as `return self.args`;
its pickle reduction does not retain the `evaluate=False` construction flag
for Add/Mul. Reading the saved inverse's pickle opcodes without restoring it
also shows ordinary `NEWOBJ` construction for these nodes. Reconstruction can
therefore canonicalize the saved Add/Mul/Pow expressions. The failed checks
test preservation of an unevaluated layout that the codec does not promise.

This is a concrete validator defect and the likely explanation of all four
failures. It does not establish that every saved expression is correct. The
corrected validator must pass before the ingredient is accepted. No saved
scientific file will be rewritten to recover the original layout.

## Prepared bounded correction; execution awaits permission

`S11c_d_reference_kernel_validate_v2.py` changes only the comparison-side
Add/Mul/Pow constructions for those four checks to their default canonical
form. The determinant domain, source and unit joins, all thresholds, saved
operand checks and numerical residual/mutation validation remain unchanged.
The numerical phase reads saved probe matrices and performs residual products;
it neither reconstructs the symbolic inverse nor repeats LU solves, old
reductions, modes or scattering calculations.

The fixed proposal is one run in the fresh
`_scratch/s11c/s11c-d-reference-kernel-20260926/validation-v2` directory,
using the existing resource guard, normalization supervisor and hook-first
launcher: at most 900 seconds, 2 GiB, zero swap, one CPU, nice 15, 32 tasks,
one native thread. Inputs are the unchanged constructor outputs, frozen by
the same full-file inventory. The preparation receipt pins the corrected
validator, launcher, manifest, diagnosis source and completion message.

No corrected run has launched. The explicit no-retry boundary requires a
user decision before this second validation attempt. The full outgoing Green
prescription, FORM and A11/A12 remain unresolved regardless of this decision.
