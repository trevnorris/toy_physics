## MUST-FIX

The directive is not safe to hand to the builder as-is because R2 conflates two distinct checks:

1. Exact serialization round-trip: compact in-memory value ↔ `_restore(srepr(compact))`.
2. Semantic preservation: restored compact value ↔ original emitted value.

The first must be retained, and the second added separately. The semantic comparison also needs explicit recursive shape/key/label checks because the exported objects are nested tuples, not scalar expressions.

## 1. Scope — correct and complete

Finding: PASS.

The three edit sites are sufficient:

- `EXPORT_ROOTS` controls publication membership ([script:48](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:48)).
- The `export_key` map is evaluated only after `emit(...)`, so deleting the increment mapping leaves its stdout emission untouched ([script:880](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:880), [script:893](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:893)).
- All row construction, closure selection, serialization, import checking, and file replacement live inside `publish` ([script:807](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:807), [script:842](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:842)).
- The final call passes already-computed emitted payloads into `publish`; compacting there need not alter construction ([script:952](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:952)).
- `cas` merely converts containers/strings to SymPy structures; it does not force expansion ([script:64](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:64)).

No physics, `emit`, or construction edit is forced. The STOP fence is correctly placed.

## 2. R1 membership — correct

Finding: PASS, with one documentary caveat.

Both roots are justified:

- S11c-d owns scattering amplitudes, resonances, and local/Bloch/WKB spectrum, not merely the Born mixing calculation ([S11c_decisions.md:47](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_decisions.md:47), [S11c_decisions.md:83](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_decisions.md:83)).
- Closure modifies the full operator, explicitly including its transverse diagonal block ([shared physics:137](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md:137)). That diagonal information is absent from the off-diagonal coupling kernel and is required for the closed spectral equation/resolvent.
- The closed kernel is separately the weakly restricted off-diagonal object needed directly for linear mixing and scattering ([shared physics:187](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md:187)).
- Asking d to derive either missing object again would repeat c2’s close or extract construction.

The increment is correctly EMIT-only. It is `extract(close(SLAB)) − extract(SLAB)` and is explicitly a comparison/export representation, not a residual construction ([shared physics:203](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md:203)). Its load-bearing use is the stdout comparator ([shared physics:394](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md:394)); the adjudication independently directs exactly this membership ([adjudication:106](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_physics_review_adjudication.md:106)).

All four cases are populated by the two nested anchoring/density loops ([script:874](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:874)).

Caveat/nit: no concrete S11c-d `IMPORT_KEYS` manifest exists yet. The preceding build directive acknowledges that and declares both prospective binds explicitly ([build directive:197](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_sympy_build_directive.md:197)). The repair directive could repeat that caveat, but it does not change the correct build membership.

## 3. R2 representation — concept correct, guard specification MUST-FIX

Finding: ordinary transparent factoring is sound; the stated verification is not sufficient as written.

A restored ordinary SymPy expression remains importable and differentiable:

- `_restore` evaluates the serialized SymPy expression into a real object ([ledger_fold.py:51](/var/projects/toy_physics/research/pde_ledger_v3/scripts/ledger_fold.py:51)).
- `load_model` imports and returns row values without wrapping or string-level indirection ([ledger_fold.py:102](/var/projects/toy_physics/research/pde_ledger_v3/scripts/ledger_fold.py:102)).

The fault is “upgrade the roundtrip”:

- The current guard compares the compact in-memory row with the value reconstructed from its serialization ([script:842](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:842), [script:848](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:848)). That exact structural round-trip remains valuable after compaction and must not be replaced.
- A new, separate semantic comparison must compare the restored compact root against the original `objects` payload passed to `publish`.

The roots are nested case/payload/channel tuples, not subtractable scalars. The directive must require:

- identical case-key sets;
- identical tuple arities, mapping keys, `Str` labels, and matrix shapes;
- exact equality for non-compacted payload metadata;
- leafwise `sp.expand(decoded_leaf - emitted_leaf) == 0` for every algebraic leaf.

This matters because the existing `difference` helper uses `zip` without checking tuple lengths ([script:100](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:100)). Reusing it naïvely could ignore a dropped trailing channel and report zero.

CSE is safe only if its temporary symbols are fully decoded before the value enters `LEDGER`; no CSE temporary or hold/wrapper may survive in the runtime value.

Also keep the explicit `expand(diff)==0` test conservative. A substitute based only on `cancel(together(diff))` can accept transformations that cancel denominator factors and change singular loci—material to d’s resonance calculation. If such a rational canonicalizer is allowed, denominator/singularity preservation needs an additional guard.

## 4. R3 hygiene — consumer-safe

Finding: PASS operationally; minor schema wording nit.

Neither `load_model` nor the fold guards bind `display`:

- Export ingestion requires only string keys and mapping rows ([ledger_fold.py:82](/var/projects/toy_physics/research/pde_ledger_v3/scripts/ledger_fold.py:82)).
- Closure follows `value` and `dimension_key` only ([ledger_fold.py:225](/var/projects/toy_physics/research/pde_ledger_v3/scripts/ledger_fold.py:225)).
- Minimality examines only row keys ([ledger_fold.py:323](/var/projects/toy_physics/research/pde_ledger_v3/scripts/ledger_fold.py:323)).
- The actual c2 parent consumer binds `row['value']`, plus `step` for the two provenance-sensitive c1 rows ([script:177](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:177)).

Thus removing the multi-megabyte duplicate does not break a fold consumer.

Nit: prefer requiring a short bounded `display` rather than offering omission. The architecture notes that dropping it outright supersedes the retained human-readable schema clause ([design:137](/var/projects/toy_physics/research/pde_ledger_v3/directives/export_ledger_bind_closure_design.md:137)). A short identifier such as the root name preserves schema with negligible size.

There is no other comparable non-binding oversized row field. `VALUE` dominates; the multigrade, dimension, branch-binding, and Fourier-profile payloads are generated evidence ([script:754](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:754)) and should remain unchanged.

## 5. R4 guards — correct set, conditional on fixing R2

Finding: PASS after the R2 split.

The existing guards are appropriate:

- F9 collision hard-stop before folding candidates ([script:820](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:820)).
- Recursive closure and own-delta selection ([script:823](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:823)).
- Exact minimality assertion ([script:825](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:825)).
- Five required build-input digests ([script:827](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:827), [shared physics:386](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md:386)).
- Parent lookup/manifest equality remains preserved outside the edit sites ([script:858](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:858)).

But `_restore` exactness and emitted-value semantic equivalence must remain two guards, not one weakened replacement. Verification should run against the generated module/fold, not only pre-serialization candidates.

## 6. Leak / overstep — correctly bounded

Finding: PASS; one artifact-side-effect nit.

There is no leaked physics target:

- “60 MB” is an observed storage baseline, not an expected physics result.
- No required post-repair byte count or specific factored expression is supplied.
- Semantic zero is a serialization invariant, not a pre-registered physics residual.
- The directive explicitly ends at build → verify → report and forbids reviews, comparators, downstream work, skill-reading, and self-review artifacts ([repair directive:71](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_export_repair_directive.md:71)).

Nit: the full run inherently writes `_measurements/S11c_c2_sympy_guard_evidence.json` at the checkpoint calls ([script:864](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:864), [script:952](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:952)). The directive should acknowledge that expected side effect so it is not mistaken for another builder-created review artifact.

## Final disposition

Do not hand it to the builder as-is. Fold exactly this blocking correction first:

1. Retain the exact compact-in-memory ↔ `_restore(srepr(compact))` round-trip for every delta row.
2. Add a separate per-root/per-case emitted ↔ restored-compact comparison with strict recursive container/key/label/shape equality and leafwise `sp.expand(diff)==0`.
3. Require CSE temporaries to be fully decoded and forbid any equivalence method that can erase denominator/singularity information without a separate singular-locus check.

Nits worth folding concurrently: mandate a short `display` instead of permitting omission, repeat the “no d manifest yet” caveat, and acknowledge the checkpoint JSON side effect.
