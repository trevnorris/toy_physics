# S11c-d sheet repair: domain and dependency scope

Continuation update: the [Fourier construction and repair record](S11c_d_sheet_repair_report.md) resolves the original counterexample in a computed open Fourier strip. The earlier gate below is preserved as history. The fixed-positive-frequency d selector has been repaired; global pole-sheet coverage remains outstanding.

2026-09-10. **Historical Phase A state, before the construction linked above: continuation contract unresolved; solver repair and Phase B regeneration had not started.** The preservation checkpoint is `55e298b70a8b289c8e6b16ebd0b627edb20d7294`. This is a builder investigation record, not physics or review clearance.

## Finding

The earlier bounded counterexample identifies a disagreement between d's half-plane sheet flag and transport from a real-momentum outgoing seed. The source/export trace does not establish an earlier producer defect. It also does not yet justify assigning a **global** physical-sheet label to d's complex-normal-momentum candidates.

The missing connection is between the supplied outgoing real-Fourier operator and the complex-`k_n` domain used to classify end modes. S11b supplies a particular continuation in complex **frequency** at real in-plane momentum. Applying that rule when the momentum itself becomes complex needs a construction of the spatial contour/domain and its relation to the frequency prescription. It cannot be justified solely by the radical equation, a convention string, or the diagnostic's chosen straight path.

This is an unresolved construction/authority boundary. **The lookup does not prove that the physics is underdetermined or that a new physical premise is necessary.** The needed domain might follow from the already supplied outgoing resolvent. That implication has not been constructed and checked. Conversely, choosing cuts or path classes that change which candidates are called physical would require justification before becoming an engine rule.

## Evidence and limits

The read-only command

```text
python3 _measurements/S11c_d_sheet_scope_lookup.py > _measurements/S11c_d_sheet_scope_lookup.json 2> /tmp/s11cd-sheet-scope-lookup.stderr
```

ran from `/var/projects/toy_physics/research/pde_ledger_v3`, exited 0, and wrote 257,164 bytes of literal JSON stdout with empty stderr. The [lookup script](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_sheet_scope_lookup.py) retrieves existing payloads and source passages; it performs no spectral solve, continuation, substitution, or algebraic validity test. Its [stdout](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_sheet_scope_lookup.json) includes the command, checkpoint, source/export SHA-256 pins, exact excerpts, and stored branch/root objects. Large stored root expressions remain in that evidence file, not this report.

| Boundary | Observed object | What this establishes |
|---|---|---|
| Frozen fold into d | b: 2,441 rows; c1: 44; c2: 70; overwrite audit `[]` | The inspected three-parent fold is available. This is not a semantic-equivalence or physics test. |
| c2 closed slab and coupling exports | Both rows have four cases, each with all five payload slots and three `COMPUTED_BRANCH_BINDINGS` | The stored input/output/middle-leg definitions are present. No missing branch slot is demonstrated. |
| c2 branch expressions | `Piecewise` definitions; `omega` and all three momentum groups have stored `real=True` assumptions | These are real-axis representations. Their presence does not implement a complex-momentum path or prove correct continuation outside that domain. |
| d reduction | Source lines 715–732 preserve separate momentum groups and collect the two rows' branch equations | The source explicitly carries branch operands. This turn did not rerun the reduced-row reconstruction or verify its algebra. |
| d channel input and algebraic solve | Positive-real-frequency input gate; `analytic` replaces the positive-frequency real-axis presentation with an algebraic radical; line 1951 assigns the half-plane flag | The new complex-momentum classification occurs inside d. The old flag has no seed/path argument. |
| S11b stored spectrum | `roots` contains a `ConditionSet` and a false closed-form indicator; `sheet_of_each_root` contains conditional roots, the bulk relation, an upper-rim convention string, and opposite conditional roots | This is not a numerically resolved collection of sheet-certified roots that d can use as an oracle. It is not evidence that S11b selected the wrong sheet for a computed root. |

The [prior diagnosis](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_sheet_continuation_diagnosis.json) remains the numerical evidence for the bounded discrepancy. Its source hash matches the engine inspected here. No fresh cache/transcript pair or second transport method was run in this turn; those are Phase B work after the gate below. The old main `.out` remains a historical run and does not validate the checkpoint's latest source.

## Authority boundary

- [S11b §1/§1b](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md:84): the harmonic convention specifies real in-plane `k`; lines 131–139 define the complex-frequency path from upper-rim values. Lines 141–155 prohibit reselection by real-axis decay/radiation conditions at complex frequency and require cut/coalescence reporting.
- [S11c-c1 §1b](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md:115): inherits this frequency continuation and distinguishes bulk-normal `q_out` from the curved problem's ansatz momentum.
- [S11c-d §2](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:354): requires the full end modes with the outgoing/physical-sheet prescription. Lines 410–413 place the outgoing-resolvent convention inside the engine's construction. [§3b](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:495) separately requires sheet, normalizability, width, and closure tests for poles.

These passages supply the spectral objects and inherited boundary condition. They do not themselves give an explicit joint complex-frequency/normal-momentum contour contract. Absence of an explicit recipe is not a defect: the builder must first determine what follows from the supplied operator and boundary condition. The existing half-plane predicate does not perform that construction.

## Gate and next bounded task

Stop before replacing the predicate or regenerating parent exports. The [accepted plan's Phase A gate](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_sheet_continuation_repair_plan.md:31) requires a supported prescription or an explicit report of the unresolved domain. The user's instruction is to stop when the work needs a change of approach. The next task must establish the contour connection before implementing a global classification rule.

The precise question is:

> For d's reduced outgoing operator at fixed real `(omega,k_parallel)`, which continuation domain/path class in complex `k_n` is connected to the real-Fourier contour, and how must that contour be transported when `omega` is continued by S11b §1b?

Resolve that question by constructing the continuation from the supplied Fourier/outgoing-resolvent operands. Distinguish a locally continued real-axis germ from a global sheet assignment, record branch-point/cut encounters and path dependence, and retain normalizability as a separate test. If the supplied boundary condition determines the needed domain, document the derivation and proceed with the existing repair plan. If more than one inequivalent admissible choice remains and changes classification, take the explicit alternatives to the physics author/orchestrator; do not silently amend shared physics.

After that resolution, Phase B can use a fresh source-pinned one-case cache/transcript, independent numerical transport methods, and paired producer/export/consumer operands to locate the earliest failing boundary. Only then choose the regeneration scope. **No upstream restart is justified by this lookup alone; upstream branch correctness also remains unverified.**

## Preserved state

No engine, physics authority, parent export, or canonical transcript was changed. No four-case regeneration, review leg, comparator, Wolfram engine, downstream step, or additional commit ran. The new scope evidence and this record are uncommitted. The separate shear-normalization and cross-engine operand/sign debts remain open. S11c-d's 12 outstanding constructions and absent export remain unchanged.
