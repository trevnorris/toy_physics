# S11c-d sheet-continuation repair plan

Execution update: [the Fourier construction and repair](S11c_d_sheet_repair_report.md) resolves the original counterexample without an upstream producer change. The d selector now records root transport and unresolved paths. Four-case regeneration completed with exit 0 and the checked main transcript is published. Generic pole-sheet work and the wider Phase E verification matrix remain in the original program.

2026-09-10. Planning record following checkpoint `55e298b70a8b289c8e6b16ebd0b627edb20d7294`.
**Historical planning state, before the construction linked above:** implementation had not started. The original plan below organizes diagnosis and a conditional repair; it does not amend or clear a physics premise.

Earlier stop, 2026-09-10 (superseded by the execution update above): Phase A source/export lookup is recorded in [the scope report](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_sheet_continuation_scope.md). Work stopped at its unresolved continuation-domain gate before solver edits or Phase B regeneration. The next bounded task is to establish the complex-momentum contour connection to the supplied outgoing resolvent.

The immediate objective is to determine where sheet information first becomes incorrect or insufficient, repair that boundary, and establish whether S11c-d can resume. Reopening earlier steps requires evidence from their own declared domains. Completing this repair will not complete the remaining S11c-d scattering/pole program.

## 1. Preserved starting point

The user-requested checkpoint was saved before this plan using DataLad. Five source/report/JSON files are ordinary Git objects; the two new `.out` files are git-annex symlinks. Their combined payload is 14,319,129 bytes; their Git pointer blobs total 276 bytes. The worktree was clean immediately after the checkpoint. It is an unreviewed preservation commit, with the defect explicitly open.

| Artifact | Checkpoint status |
|---|---|
| `scripts/S11c_d_mixing_scattering_sympy_audit.py` | Latest unfinished engine; SHA-256 `263e3b191699e5e7df6c695bfd20617e1d1afa8dc6529aeee75016582d36ed4b`. |
| `scripts/out/S11c_d_channel_reentry_preflight.out` | Completed one-case run, 14,310,873 bytes; predates the latest slab-current refinement. |
| `scripts/out/S11c_d_sheet_continuation_probe.out` | Completed bounded diagnostic, 8,256 bytes. |
| `scripts/out/S11c_d_mixing_scattering_sympy_audit.out` | Earlier full run, 41,302,466 bytes; deliberately unchanged after the interrupted regeneration. It does not validate the latest source. |
| `scripts/S11c_d_exports.py` | Absent; required export roots remain uncomputed. |

The [diagnosis record](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_sheet_continuation_diagnosis.json) contains transcript hashes, input case, continuation refinements, and current-check measurements. The builder report stored in checkpoint `55e298b` records the preceding stop; its “no commit” and “uncommitted source” wording describes the state before that checkpoint. The working [builder report](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_sympy_builder_report.md) now records the repair.

Preserve the prior inertia repair at `a74da30afefffab345af910bbc78db6d07ec7ac2`. The separate shear-normalization discrepancy and c2 cross-engine operand/sign debts remain open; this plan does not bundle them into the sheet repair.

## 2. What is established, and what still needs a decision

At checkpoint `55e298b`, the prototype's `PHYSICAL_BULK_SHEET` assignment in `FullPencilModes.solve_sample` used a radical half-plane test. At the recorded reference `LAB_HELD × RHO4_CONSTANT` input, the candidate at

`k_n = 0.607303653886 + 0.008491689257 i`

is labeled physical with `q = 0.080656407478 − 6.393830415310 i`. Transport from its real-momentum outgoing seed along the diagnostic's straight path produces the opposite radical. The 16/64/256-step results agree; the recorded opposite-root residual is `1.39e-17`, maximum radical-equation residual is `1.42e-14`, and minimum distance to a computed branch point is `0.636783449147` in inverse-length reference units. Here `k_n` has dimension `[L^-1]`, while the engine radical `q = c_s0 q_out` has dimension `[T^-1]`.

This establishes disagreement between the flag and that continuation path. It does not establish a global physical sheet in joint complex frequency and momentum, nor an upstream defect. Small pencil, nullspace, and frequency-normalization residuals can hold on either radical sheet.

[S11b §1b](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md:113) supplies continuation in complex **frequency**, with real in-plane momentum: start from the upper-rim real-axis value, continue along the specified frequency rays, and report cuts/coalescence without reselecting by spatial decay. S11c-d also evaluates complex normal **momentum**. The first investigation must determine what the cleared authorities imply for that extension, including the intended paths and domain.

**Decision gate:** if the required complex-momentum or joint-domain prescription needs an additional physical premise, stop and present the precise missing convention, relevant authority passages, and consequences to the user/orchestrator. Do not promote the diagnostic's straight path into a global prescription, rewrite shared physics, or silently change the spectral problem.

## 3. Phase A — trace the branch contract and dependency graph

Read the applicable passages of the governing [program brief](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_sympy_build_PROGRAM_BRIEF.md), [d shared physics](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md), and [build directive](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_sympy_build_directive.md), together with inherited b/c1/c2 conventions. Record, for each consumer, which variables are real, which may become complex, the seed/cut/path data supplied, and where those data are represented in executable objects.

Trace backward from the failing d consumer, using the following source map as leads rather than evidence that each producer is defective:

| Boundary | What to inspect |
|---|---|
| d mode construction | `FullPencilModes.analytic` (line 1778) removes the positive-frequency real-axis `Piecewise` presentation and introduces an algebraic radical; `solve_sample` (line 1858) enumerates candidates and assigns flags. Determine which branch information survives each operation. |
| d Fourier reduction | `EdgeReduction` (line 700) collects `COMPUTED_BRANCH_BINDINGS` from both closed roots and reduces the input, output, and middle momentum legs. Check all five payload slots, symbol identity, substitutions, dimensions, and restored grades. Preserve this reduction and the three-parent fold. |
| c2 closed export | `outgoing_spectral` (line 437) computes real-axis `Piecewise` bindings; `emit` adds them to closed payloads (line 756). Compare producer objects with imported bindings on their declared domain before examining complex extensions. |
| c1 closure | It retains separate bulk-normal momentum legs and emits one-sided branch-sign and momentum-freeze controls (`task_branch`, line 2006). Determine whether the relevant branch information is inherited, symbolic, or evaluated. |
| b and a inheritance | b imports `S11c_a_exports.py`; a imports `S11b_exports.py`. Trace the actual rows used by c1/c2/d, distinguishing inherited rows from newly derived geometry or slab rows. |
| S11b branch authority/producer | Check the real-axis bulk ansatz/flux derivation and the active `S11b_interface_coupling_law_sympy_audit.py` producer. `SHEET_OF_EACH_ROOT` (line 1955) includes an upper-rim continuation label; determine whether the associated root data actually implement the supplied continuation or retain unresolved algebraic roots. A label alone is not verification. |

The current load graph is a chain through `S11b → a → b → c1`, with c2 folding **b and c1**, and d folding **b, c1, and c2**. Source-hash pins may introduce additional refresh obligations. Record both value dependencies and provenance dependencies; do not assume every earlier stage needs a new physical calculation.

**Phase A exit artifact:** a compact scope record with source/export hashes, exact row/slot lookup witnesses, domain contracts, and one of: a supported prescription to implement, or the explicit missing-premise stop described above. No upstream correctness claim is earned by source inspection alone.

## 4. Phase B — reproduce and locate the earliest failing boundary

Develop on the existing explicit `LAB_HELD__RHO4_CONSTANT` input. Produce a fresh cache/transcript pair with a manifest tying together the engine, imported exports, case, profile/parameter JSON, unit frame, cache, and command. The preserved diagnostic hashes its inputs but does not establish every cache-to-source/case relationship; close that provenance gap before using it to assign upstream responsibility.

Reproduce the counterexample using the radical relation extracted from the computed reduced pencil. Keep the old candidate index only as a historical locator. Match fresh candidates by their spectral values and, where needed, eigenspace/projector data; list ordering is not a mode identity.

Construct two numerical continuation routes with separate numerical failure modes: an adaptive algebraic root tracker and an implicit differential-equation transport derived from the same computed radical relation. Their shared relation and seed are disclosed; this checks transport, not an independent derivation of the physics. Derive/check the real-axis seed from the supplied bulk ansatz and computed outward flux/decay conditions. Do not force the observed opposite sign as a target answer.

For each suspected boundary, evaluate paired producer/export/consumer operands at the **same** bindings and in the **same declared domain**. Emit branch values, path/seed identifiers, radical and operator residuals, and uncertainties. Check both closed roots and all three momentum legs. Compare physical `q_out` only after computing the conversion from the d radical, including its `c_s0` factor.

A real-axis binding that works on its declared real Fourier domain is not shown defective merely because its consumer needs a complex extension. Conversely, satisfying the squared radical equation, agreeing at PIT points, or carrying a convention string does not verify analytic continuation.

**Phase B exit artifact:** a table identifying the earliest demonstrated failure, the reproducing input/path and literal residuals, and upstream boundaries that are verified, outside the tested domain, or unresolved. Preserve unsuccessful/inconclusive results explicitly.

## 5. Phase C — choose the smallest supported repair

| Finding | Repair and regeneration scope |
|---|---|
| Upstream operands agree within their supplied domains; d discards or fails to transport branch data | Repair d's analytic representation/continuation and its dependent consumers. Regenerate d diagnostics and the main four-case transcript after the focused checks. Leave upstream exports unchanged. |
| c2 loses/misrepresents required branch data during closure or export | Repair c2 at the computed operand/binding boundary. Regenerate its affected transcript/export and check both closed roots, bind closure, and actual consumers; then regenerate d. Preserve b/c1 values unless the trace identifies an additional dependency. |
| c1 or an earlier active producer violates its own supplied branch convention | Repair the earliest demonstrated producer. Recompute affected descendants in dependency order through the actual a/b/c1/c2/d graph. Record physical value changes separately from inherited-row or source-pin refreshes. Recheck the existing inertia repair wherever b-derived data are regenerated. |
| Authorities do not settle the extension needed by d, or evidence remains inconclusive | Stop at a runnable diagnostic checkpoint and report the unresolved premise/domain or numerical obstruction. Do not assign a physical flag by a fallback sign test. |

The scope record must be concrete before edits: failing function/row, affected consumers, anticipated regenerated files, and checks that establish the retained upstream boundary. If investigation expands the work into a different physics premise, another open debt, or a larger project reorganization, notify the user and stop before that expansion, as requested.

## 6. Phase D — implement branch transport consistently

Subject to the domain decision, carry a computed branch record with each evaluated spectral object: real-axis seed, continuation variables and path, cut/branch-point encounters, refinement/error data, and resolved/opposite/unresolved status relative to the declared prescription. Compute winding information where needed by that prescription. Never replace unresolved status with a global real/imaginary half-plane test.

Use that same branch record in the full reference/end pencils, radical chain-rule derivatives, nearby-frequency checks, lifted closed-field residuals, and local mode jets. Audit independent re-evaluations of the pencil for accidental principal-root reselection. Retain candidates on other sheets with their labels; do not delete a root to obtain a desired channel count.

Keep sheet membership separate from spatial decay/normalizability, frequency growth/decay, and incoming/outgoing flux. A spectral slope or a corrected sheet flag is not the unfinished full bulk current. Degenerate eigenspaces require subspace tracking; threshold or unresolved rank cases must retain their domain limitations.

Extend the existing engine and probe; preserve imported key ownership, the five-slot payload census, and the d one-dimensional reduction. All new physical expressions must be reached from the action/ansatz or computed imported/reduced operands. Emit computed objects and residuals before structural guards, with restored `[L,T,M]` dimensions and `(eps,eta,sigma_W)`/lambda order. New write-keys must be fresh, injective lowerCamel keys. Interpretation belongs in the report.

## 7. Phase E — focused verification, then complete regeneration

Use this coverage matrix to avoid treating the one counterexample as generic verification:

| Domain/check | Required evidence |
|---|---|
| Real frequency and real momentum | Both frequency signs, propagating and evanescent bulk regimes, and both face orientations; seed values and computed outward flux/decay operands. The current input has a closed propagating bulk cone, so additional declared physical inputs are needed. |
| Complex frequency, real momentum | Upper- and lower-half-plane paths from the supplied upper rim, prescribed cuts, opposite-sheet sensitivity, and explicit threshold/coalescence handling. Do not reapply spatial decay at complex frequency. |
| Complex normal momentum | The preserved counterexample, approved seed/path families, step and precision refinement, and comparison of the two transport methods. Test contractible path deformations away from branch points and separately record paths with different winding. |
| Branch points and cuts | Near-threshold conditioning, exact threshold limitations, cut crossings, and paths rejected or unresolved by the chosen domain contract. No silent reseeding. |
| Mode consistency | Radical equation, complete reduced-pencil/nullspace and lifted-field residuals, frequency normalization and finite differences with transported branches, and reference/left/right mode matching including degenerate subspaces. |
| Boundary liveness | Input/output/middle-leg substitutions and one-sided branch-sign or momentum-freeze mutations entering at actual operands, with both operands and computed residuals. No `A−A` controls. |
| Prior repair retention | Recompute the relevant kinetic-variation and propagating transverse-channel checks from the inherited energy/operator. Historical roots are comparison data, not values to hard-code. Keep the separate shear and sign debts visible. |

Report numerical precision, scale-aware residuals, refinement changes, and domain exclusions. Establish acceptance criteria from conditioning and independent error estimates; agreement between two calls to the same selector is insufficient. Finite coverage does not prove a global sheet construction.

Iterate cheaply on one case. Run at most one memory-heavy CAS process at a time and measure the actual process's wall time and peak RSS. Preserve cancellation/exit status and stderr. Once the implementation and focused checks are stable, run each affected upstream producer's required full case/task coverage in dependency order, then run all four d `(alpha,rho)` cases once against the refreshed parents. If a run fails, leave its incomplete capture unpublished and report what remains.

Capture stdout to a fresh scratch file, never through an existing annex symlink. Validate completion, source pins, case/tag coverage, dimensions/grades, and relevant consumer/export guards before replacing canonical files. Publish by replacing the directory entry with the completed payload; preserve prior annex objects. Save `.out` files with DataLad/git-annex and scripts, exports, JSON, and reports with ordinary Git according to the README/attributes. Any later commit must describe the actual repair and its remaining limitations; this plan does not confer review clearance.

Keep heavy output carrier-first with numeric-PIT fingerprints and SHA digests, at the existing tens-of-MB scale. Keep physical numerical continuation evidence separate from algebraic PIT. For changed upstream exports, run the existing minimal-delta, bind-closure, semantic-equivalence, and consumer checks appropriate to those producers.

## 8. Deliverables and return to S11c-d

The repair should leave:

1. The upstream scope/domain decision and a compact before/after boundary record, including unresolved domains and any required authority clarification.
2. The repaired engine/probe and only the upstream producers actually implicated by that record.
3. Completed diagnostic `.out` files, input/cache/source manifests, and refreshed canonical transcripts/exports for affected producers with validated provenance.
4. A concise updated d builder report linking scalar evidence, listing which constructions were completed, and stating exactly what remains.

Do not create `S11c_d_exports.py` merely to mark this repair finished. Its required scattering, pole, survival, and weak-coefficient roots must first be computed under the original program. Never substitute an empty S-matrix or pole set for an unexecuted construction.

Return to the original d program only when the continuation contract is settled for the needed domain, the original discrepancy is explained by measured before/after data, the affected producer chain and all four d cases have completed, and sheet state is used consistently by the mode consumers. A bounded diagnostic repair alone does not empty `GENERIC_DOMAIN_SHEET_CONTINUATION`; retain that TODO until its promised coverage exists.

Then continue the remaining spectrum work and bulk-current/flux construction before scattering and poles, in the brief's dependency order. The current 12 outstanding constructions remain the baseline until computation earns their removal. The lane remains **build → run → report → STOP**: no review legs, comparator, Wolfram engine, downstream stage, or `/build`/`/review-legs` skill invocation is part of this repair plan.

Historical endpoint of this original planning record: work stopped after planning. The execution update at the top and linked repair report describe the subsequent implementation and regeneration.
