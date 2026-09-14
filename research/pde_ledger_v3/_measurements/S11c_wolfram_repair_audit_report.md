# S11c Mathematica audit of the four SymPy repairs

Update: the user subsequently approved the pressure-trace repair. Its focused
and native four-case checks passed; see the [repair report](S11c_wolfram_pressure_trace_repair_report.md)
and [native checkpoint](S11c_wolfram_pressure_trace_native_checkpoint.json).
The original audit disposition and counterexample below are retained as history.

**A pressure-trace defect is present in the native Mathematica c2 N6 engine.**
The user-requested audit is checkpointed; producer repair awaits approval under
the standing instruction to stop when another upstream repair is needed.

| Repair tested | Mathematica disposition and scope |
| --- | --- |
| Conservative inertia sign | ABSENT_ON_AUDITED_DOMAIN: b's raw, frozen and live constructors use the positive kinetic variation relative to their stored-energy rows. |
| Mechanical-load normalization | ABSENT_ON_AUDITED_DOMAIN: b's prescribed physical face force enters the stored row with the opposite sign, with the correct virtual-work and power pairing. |
| Reference/physical pressure trace | PRESENT: c2 N6 assigns its physical-face pressure response to a reference-pressure slot and shifts it again. |
| Thickness kinetic coordinate | ABSENT_ON_AUDITED_DOMAIN: b uses the physical displacement `W0 eW`, including on a symbolic nonuniform background. |

The [checkpoint](S11c_wolfram_repair_audit_checkpoint.json) pins exact native
definitions, supplied authorities, pre-existing outputs and executed scripts.
The mechanical audit computes **200 zero residuals**: 24 inertia components,
64 constrained stored-energy probe components, 16 face-work comparisons,
80 generalized load components and 16 power comparisons. It covers both
coordinate routes, anchorings, density representatives and individual faces.
Negative-kinetic and wrong-coordinate action mutations produce 24 and six
nonzero components; load sign/scaling mutations produce 112. The stored-energy
probe tests variation/assembly orientation; it does not regenerate the N15
basis or establish full closed-c2 energy balance.

The trace audit executes native c1's closure solve and c2's response, pressure
image assignments, normal continuation and affine face evaluation. An independent
exact wave solves the displaced boundary problem. All 24 retained c1/response/
independent-affine control residuals vanish. The native fold has a nonzero
first-order height error on both faces and both density representatives. The
common map is independent of anchoring on this constant-height locus; its eight
case labels are not eight independent generic-domain proofs.

For the recorded lossless witness, `omega=13`, `|k|=12`, `q=5`, `c=1`,
`rhoM=5`, velocity drive `11`, and displacement height `0.01`, the physical
pressure is 143. The extra native shift contributes **143 i / 20** in the restored
pressure unit `L^-2 T^-2 M`. The dispersion residual is zero and all three
matching denominators are `5 i`. The independent reference-pressure continuation
reconstructs the physical pressure through the retained order. This counterexample
establishes presence; the nonconstant-profile ordered eta/sigma inverse remains
repair work, and no global sheet, threshold or pole coverage is inferred.

Final runs used `math -script <snapshotted-path>` serially, completing in 18.69 s
and 19.13 s with exit zero and empty stderr. Published evidence comprises
[393 mechanical records](../mathematica/out/S11c_wolfram_mechanical_repair_audit.out)
and [85 trace records](../mathematica/out/S11c_wolfram_pressure_trace_audit.out),
504,957 and 224,106 bytes. Each physical packet retains raw/retained/discarded
operands with dimensions and grades. Final trace metadata uses explicit
coefficient extraction to exclude beyond-order rational remainders; the earlier
diagnostic and all source snapshots remain in repository `_scratch/s11c/`.

No native producer, old transcript, SymPy export, accepted normalization result,
authority, S10/Lean file or retained solver/export contract changed. These focused
Wolfram files are audit instruments, not a full c2 self-energy or d engine.
No review/comparator, automation or push was run. No calculation remains active.
The next action is approval of the [c2 repair plan](S11c_wolfram_pressure_trace_repair_plan.md).
