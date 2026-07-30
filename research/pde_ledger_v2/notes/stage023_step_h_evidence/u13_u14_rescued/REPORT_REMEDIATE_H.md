# stage023 step-(h) remediation report

Work was performed from `research/pde_ledger_v2/scripts` on
`ledger-v2-rebuild`. No commit, checkout, restore, or stash was performed.
No value in `SOURCED_DIMS` or `EXPECTED_DIMS` was changed, and no physics was
changed. `ledger_dimensions.py`, the `.wl`, the `.out`, and the three named
validator/generator scripts were not edited.

## Custody and baseline

The current working bytes were copied before editing. The original and backup
`git hash-object` values matched:

```text
6ecc9fd00a1879c62fc1710d2f2c5c41ce741324  scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py
6ecc9fd00a1879c62fc1710d2f2c5c41ce741324  _scratch/stage023_h/backups/ledger_stage023_nullspace_underdetermination_sympy_audit.py.pre_h
8700fa0bd583fbd1c6d47cc82bbf81af803b94a1  scripts/output/ledger_stage023_nullspace_underdetermination_sympy_audit.txt
8700fa0bd583fbd1c6d47cc82bbf81af803b94a1  _scratch/stage023_h/backups/ledger_stage023_nullspace_underdetermination_sympy_audit.txt.pre_h
8273de4ae8e3ff666bac6c220c4ede0fafb4444b  scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.dimensions.txt
8273de4ae8e3ff666bac6c220c4ede0fafb4444b  _scratch/stage023_h/backups/ledger_stage023_nullspace_underdetermination_sympy_audit.dimensions.txt.pre_h
```

The canonical-table generator necessarily writes `CANONICAL_DIMENSIONS.md`.
That already-dirty file was also copied before acceptance command 4. Its
pre-command, post-command, and backup object hash was the same:

```text
75568d995214e7e0d7551593587e4c7b6961edb7
```

Thus the required generator execution did not change the other party's
Markdown bytes; no restore was needed.

## R1 — reachability of the two added branches

### Determination before editing

A complete top-level stage run was executed under a call/line profiler. It
exited 0. The argument-type and line-hit evidence was:

```text
RUN_EXIT=0
ARG_TYPES|fmt|Add=2|Half=1|Mul=5|NegativeOne=1|One=1|Symbol=48|Zero=4
ARG_TYPES|assert_no_float|Add=59|Integer=5|Mul=45|NegativeOne=2|One=8|Zero=449|int=20
LINE_HITS|line=74|hits=62
LINE_HITS|line=75|hits=0
LINE_HITS|line=83|hits=588
LINE_HITS|line=84|hits=0
LINE_HITS|line=85|hits=0
LINE_HITS|line=86|hits=0
```

The `isinstance` tests ran, but no `Dimension` reached either function and
neither branch body ran.

### `fmt`: remove the unreachable special case

The `Dimension` special case was redundant: the existing generic path calls
`compact`, which returns a non-SymPy object unchanged, and `sp.sstr` renders
the `Dimension` through its existing string representation. Removing the
branch reduces dead stage-local code without changing the shared module.

Executed before/after evidence for the exact case that the removed branch
appeared to cover:

```text
BEFORE_FMT=(1/2,0,-1)
AFTER_FMT=(1/2,0,-1)
AFTER_STR=(1/2,0,-1)
IDENTICAL=True
```

Diff:

```diff
 def fmt(expr: Any) -> str:
     if isinstance(expr, bool):
         return "True" if expr else "False"
     if isinstance(expr, str):
         return expr
-    if isinstance(expr, Dimension):
-        return str(expr)
     try:
         return sp.sstr(compact(expr))
```

### `assert_no_float`: repair the call path

The dimension branch performs a meaningful check for the newly emitted data,
so it was retained and made reachable at the sidecar boundary. All 29 records
are now traversed before they are returned to `emit_dimension_sidecar`.
This is stage-local; no `ledger_dimensions.py` change was required.

Diff:

```diff
 def dimension_records(dimension: dict[str, Any]) -> dict[str, Dimension]:
     computed = dimension["computed"]
-    return {
+    records = {
         ...
     }
+    assert_no_float("dimension_records", records)
+    return records
```

A complete instrumented `main()` after the repair returned 0, preserved the
tally, and reached the `Dimension` branch once for each emitted record:

```text
MAIN_RETURN=0
COUNTS|pass=111|fail=0
ARG_TYPES|fmt|Add=2|Half=1|Mul=5|NegativeOne=1|One=1|Symbol=48|Zero=4
ARG_TYPES|assert_no_float|Add=59|Dimension=29|Half=2|Integer=14|Mul=45|NegativeOne=13|One=25|Rational=2|Zero=495|dict=1|int=20
```

The executed mutation used an otherwise-valid `Dimension` whose
`sourced_dims.a.L` backing exponent was replaced with `sp.Float("1.0")`.
The guard fired at the intended labelled location:

```text
MUTATION=sp.Float('1.0') at sourced_dims.a.L
GUARD_FIRED=dimension_records.sourced_dims.a.L: Float atom(s) in exact audit expression: {1.00000000000000}
```

No pre-existing vacuous, entailed, or shadowed assertion was identified or
changed during this scoped remediation.

## R2 — sidecar failure and terminal status ordering

### Before

`ledger_dimensions.emit_dimension_sidecar` was replaced at runtime with a
function raising `OSError("injected sidecar emission fault")`. The unedited
stage exited 1 and printed the contradiction:

```text
AUDIT_STATUS=PASS
PHYSICS_VERDICT=FAIL_UNDERDETERMINED_NOT_PREDICTIVE (2/2, COMPLETING)
UNCAUGHT exception: OSError('injected sidecar emission fault')
TALLY sympy: 111 pass + 1 fail = 112 checks
OVERALL FAIL
```

### Change

Sidecar construction/emission now precedes `print_verdict_labels`, which owns
the terminal audit and physics status banners:

```diff
     }
     run_verdict_seam_and_scope(all_flags)
-    print_verdict_labels()
     emit_dimension_sidecar(
         __file__,
         dimension_records(run_gate(Mutation())["dimension"]),
     )
+    print_verdict_labels()
     return 0
```

### After

The same injected fault exited 1 and emitted:

```text
UNCAUGHT exception: OSError('injected sidecar emission fault')
TALLY sympy: 111 pass + 1 fail = 112 checks
OVERALL FAIL
```

There was no `AUDIT_STATUS` or `PHYSICS_VERDICT` line in the post-change fault
stdout. The normal successful stdout remains unchanged by the reorder because
sidecar emission itself prints nothing.

## R3 — committed Python transcript

The current transcript body differed from the fresh final run on nine lines.
The transcript was regenerated by replacing exactly those nine body lines:
the engine description, seven tuple renderings, and the sidecar-I/O
description.

```diff
-Engine: genuine SymPy symbolic Jacobian Matrix.rank; exact, float-free, standalone, zero file I/O.
+Engine: genuine SymPy symbolic Jacobian Matrix.rank; exact, float-free, standalone; audit results on stdout and labelled dimensions in a deterministic sidecar.
-  [A0]=(0, 1, -1) expected=(0, 1, -1)
+  [A0]=(0,1,-1) expected=(0,1,-1)
-  [A1]=(1, 1, -1) expected=(1, 1, -1)
+  [A1]=(1,1,-1) expected=(1,1,-1)
-  [T0]=(0, 0, 0) expected=(0, 0, 0)
+  [T0]=(0,0,0) expected=(0,0,0)
-  [T1]=(0, 0, 0) expected=(0, 0, 0)
+  [T1]=(0,0,0) expected=(0,0,0)
-  [epsilon0]=(0, 0, 0) expected=(0, 0, 0)
+  [epsilon0]=(0,0,0) expected=(0,0,0)
-  [epsilon1]=(0, 0, 0) expected=(0, 0, 0)
+  [epsilon1]=(0,0,0) expected=(0,0,0)
-  [P0_physical]=(0, 0, 0) expected=(0, 0, 0)
+  [P0_physical]=(0,0,0) expected=(0,0,0)
-  ZERO file I/O: standalone print-only audit; no scratch YAML, report, note, card, LaTeX, or registration write.
+  DIMENSION sidecar I/O: no scratch YAML, report, note, card, LaTeX, or registration write.
```

`git diff --numstat` reports `9 9` for the transcript. Separate `cmp` checks
against the pre-remediation backup returned 0 for both the first eight lines
and the last two lines, proving the committed wrapper shape was preserved.
The measured body command and its `cmp` against the fresh final stdout both
returned 0.

The producer also refreshed the already-present sidecar's source binding from
`10f125307d0a3a0bd09b19d32e0f4f66ce6512a4fc417f58ac2fca11c670edc4`
to
`ac218fbfaec569dc00370cfa332b6cd13ecf08d3a5ae26721c190b73fb4e0c66`;
the pinned module digest remained
`7cb8b7c93d189320f3ff857a665162e3290d1628e4511b2b387dd9a4f08e5bea`.

## Acceptance

All commands were run from `research/pde_ledger_v2/scripts` with the required
600-second cap where applicable.

1. `timeout 600 python3 ledger_stage023_nullspace_underdetermination_sympy_audit.py`

   ```text
   TALLY sympy: 111 pass + 0 fail = 111 checks
   OVERALL PASS
   EXIT_CODE: 0
   ```

2. `timeout 600 python3 compare_dimension_artifacts.py 023`

   ```text
   ARTIFACT_NAME_SET|stage=stage023|py=29|wl=29|shared=29|py_only=0|wl_only=0|source_coverage=not_checked
   RESULT|stage=stage023|status=PASS|mismatches=0
   EXIT_CODE: 0
   ```

3. `timeout 600 python3 check_ledger_dimensions_pin.py`

   ```text
   MODULE_PIN_OK|module=scripts/ledger_dimensions.py|authority=scripts/ledger_dimensions.accepted.sha256|accepted_sha256=7cb8b7c93d189320f3ff857a665162e3290d1628e4511b2b387dd9a4f08e5bea|consumer_artifacts=(none)
   EXIT_CODE: 0
   ```

4. `timeout 600 python3 generate_canonical_dimension_table.py`

   ```text
   WROTE|path=CANONICAL_DIMENSIONS.md|quantities=122|candidate_groups=6|needs_adjudication=3
   EXIT_CODE: 0
   ```

5. `tail -n +9 output/ledger_stage023_nullspace_underdetermination_sympy_audit.txt | head -n -2`

   ```text
   EXIT_CODE: 0
   ```

   `cmp` of that output against the stdout captured by acceptance command 1:

   ```text
   EXIT_CODE: 0
   ```

## Final repository state

`git status --short`:

```text
 M docs/development_pipeline.md
 M research/pde_ledger_v2/CANONICAL_DIMENSIONS.md
 M research/pde_ledger_v2/manifests/DIMENSION_REWRITE.md
 M research/pde_ledger_v2/scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py
 M research/pde_ledger_v2/scripts/output/ledger_stage023_nullspace_underdetermination_sympy_audit.txt
?? research/pde_ledger_v2/scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.dimensions.txt
```

The first three entries were pre-existing, unrelated work. Final
`git hash-object` values for every in-scope file modified by this remediation:

```text
0e01ceb07da93a24a37cee2aa6d94dad4e6b25f2  research/pde_ledger_v2/scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py
f263c1e9c73a7181123952c0b9c2773b1b8f5929  research/pde_ledger_v2/scripts/output/ledger_stage023_nullspace_underdetermination_sympy_audit.txt
e2de91a8241068cbbaff5bb99b9802a21c4a66c9  research/pde_ledger_v2/scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.dimensions.txt
```

---

# Round 2 — R4 re-decision and R5 card correction

Work remained on `ledger-v2-rebuild`. No commit, checkout, restore, or stash
was performed. No value in `SOURCED_DIMS` or `EXPECTED_DIMS`, no physics,
`.wl`, `.out`, `ledger_dimensions.py`, validator/generator script, note, or
manifest was edited.

## Round-2 custody

Before editing, the current bytes were copied with `cp` to
`/tmp/stage023_h_round2.vYHoVr`. Each original and backup pair had the same
`git hash-object`:

```text
0e01ceb07da93a24a37cee2aa6d94dad4e6b25f2
  scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py
  /tmp/stage023_h_round2.vYHoVr/ledger_stage023_nullspace_underdetermination_sympy_audit.py.pre_round2
e2de91a8241068cbbaff5bb99b9802a21c4a66c9
  scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.dimensions.txt
  /tmp/stage023_h_round2.vYHoVr/ledger_stage023_nullspace_underdetermination_sympy_audit.dimensions.txt.pre_round2
f263c1e9c73a7181123952c0b9c2773b1b8f5929
  scripts/output/ledger_stage023_nullspace_underdetermination_sympy_audit.txt
  /tmp/stage023_h_round2.vYHoVr/ledger_stage023_nullspace_underdetermination_sympy_audit.txt.pre_round2
effda1bda6e20b4e84b5684958e6406280ad4156
  paper/stages/stage_023.tex
  /tmp/stage023_h_round2.vYHoVr/stage_023.tex.pre_round2
adc3c614e96b6b8ba368dfbf32839701742e7697
  _scratch/stage023_h/REPORT_REMEDIATE_H.md
  /tmp/stage023_h_round2.vYHoVr/REPORT_REMEDIATE_H.md.pre_round2
75568d995214e7e0d7551593587e4c7b6961edb7
  CANONICAL_DIMENSIONS.md
  /tmp/stage023_h_round2.vYHoVr/CANONICAL_DIMENSIONS.md.pre_round2
```

No restore was needed. The acceptance generator left
`CANONICAL_DIMENSIONS.md` byte-identical to its backup, and the committed
Python transcript likewise remained byte-identical to its backup.

## R4 — option (b)

### Decision

Option **(b)** was applied: restore the `Dimension` branch in `fmt` and remove
the round-1 call that activated `assert_no_float` in `dimension_records`.

The activated check cannot reject a `Dimension` the stage or shared module
can construct. Keeping it would present an assertion as protection even
though `_exact` has already converted every accepted float input to a
`Rational`. Its only demonstrated failure requires post-construction
replacement of the private `_exponents` slot, and neither the stage nor the
module has such a path. Its failure surface is also misleading. Repairing
that handler for an impossible guard would expand this representation
refactor into the class of pre-existing assertion repair forbidden by
`DIMENSION_REWRITE.md` §10. Restoring the two helper shapes also returns
stage023 to the same shape as stages 012, 013, 016, and 018.

### Executed constructor/coercion evidence

The following probe used the live `ledger_dimensions.py`:

```text
Dim(0.5,0,0) -> (1/2,0,0)|types=Half,Zero,Zero|float_atoms=set()
Dim(sp.Float('1.0'),0,0) -> (1,0,0)|types=One,Zero,Zero|float_atoms=set()
from_mapping({L:0.1}) -> (3602879701896397/36028797018963968,0,0)|types=Rational,Zero,Zero|float_atoms=set()
Dim(1,0,0)**0.5 -> (1/2,0,0)|types=Half,Zero,Zero|float_atoms=set()
mul -> (1,1,0)|types=One,One,Zero|float_atoms=set()
truediv -> (1,-1,0)|types=One,NegativeOne,Zero|float_atoms=set()
mapping_item_overwrite=TypeError: 'mappingproxy' object does not support item assignment
PROBE_EXIT=0
```

Calling stage023's guard on the four reported cases was silent:

```text
GUARD_SILENT|Dim(0.5,0,0)|(1/2,0,0)
GUARD_SILENT|Dim(sp.Float('1.0'),0,0)|(1,0,0)
GUARD_SILENT|from_mapping({L:0.1})|(3602879701896397/36028797018963968,0,0)
GUARD_SILENT|Dim(1,0,0)**0.5|(1/2,0,0)
```

Replacing `_exponents` after construction remained the only demonstrated way
to fire it:

```text
PRIVATE_OVERWRITE_GUARD_FIRED|corrupt.L: Float atom(s) in exact audit expression: {1.00000000000000}
GUARD_PROBE_EXIT=0
```

Static source and AST execution found the module's two explicit
`Dimension(...)` calls inside `DimensionBasis.__call__` and
`DimensionBasis.from_mapping`. `Dimension.__init__` normalizes again through
`_exact`; multiplication and division route through `_combine` to
`from_mapping`, and powers route directly to `from_mapping`. The private-slot
write scan reported:

```text
AST_PRIVATE_EXPONENT_STORES|ledger_dimensions.py|[106]
AST_SETATTR_CALLS|ledger_dimensions.py|[]
AST_PRIVATE_EXPONENT_STORES|ledger_stage023_nullspace_underdetermination_sympy_audit.py|[]
AST_SETATTR_CALLS|ledger_stage023_nullspace_underdetermination_sympy_audit.py|[]
SOURCE_AUDIT_EXIT=0
```

Line 106 is the constructor's initial
`self._exponents = MappingProxyType(normalized)`, not an overwrite. Thus no
stage or module path overwrites an already-constructed exponent mapping.

The sibling-shape probe reported:

```text
ledger_stage012_dtn_pole_ladder_robin_sympy_audit.py|fmt_Dimension_branch=True|assert_no_float_Dimension_branch=True|
ledger_stage013_breathing_harmonic_mk_projection_sympy_audit.py|fmt_Dimension_branch=True|assert_no_float_Dimension_branch=True|
ledger_stage016_l2_so3_covariance_sympy_audit.py|fmt_Dimension_branch=True|assert_no_float_Dimension_branch=True|
ledger_stage018_dtn_hankel_fingerprint_sympy_audit.py|fmt_Dimension_branch=True|assert_no_float_Dimension_branch=True|
SIBLING_SHAPE_EXIT=0
```

### Executed failure-surface evidence

For a real `__main__` execution, `_exact` was fault-injected at runtime to
produce `sp.Float` exponents. This preserved the normal 111 checks and made
the activated guard raise `AuditFailure` after them. The handler output was:

```text
=======
Tallies
=======
TALLY sympy: 111 pass + 0 fail = 111 checks
OVERALL FAIL
INJECTED_PROCESS_EXIT=1
HANDLER_PROBE_WRAPPER_EXIT=0
```

An `rg` scan found no `dimension_records`, `Float atom`, or `UNCAUGHT` line in
that output. This reproduces the review claim: the exception is not counted
and the offending record is not named.

### Round-2 R4 diff

```diff
 def fmt(expr: Any) -> str:
     if isinstance(expr, bool):
         return "True" if expr else "False"
     if isinstance(expr, str):
         return expr
+    if isinstance(expr, Dimension):
+        return str(expr)
     try:
         return sp.sstr(compact(expr))

 def dimension_records(dimension: dict[str, Any]) -> dict[str, Dimension]:
     computed = dimension["computed"]
-    records = {
+    return {
         ...
     }
-    assert_no_float("dimension_records", records)
-    return records
```

The post-change execution installed a raising sentinel in place of
`assert_no_float`, then called `dimension_records`. It proves the assertion is
not reached at the sidecar boundary:

```text
POST_STATE_RECORD_COUNT=29
POST_STATE_ASSERT_NO_FLOAT_CALLS=0
POST_STATE_FMT_DIMENSION=(1,0,-1)
POST_STATE_PROBE_EXIT=0
POST_STATE_SHAPE|fmt|Dimension_branch=True
POST_STATE_SHAPE|assert_no_float|Dimension_branch=True
POST_STATE_SHAPE_EXIT=0
```

### Report-only module defect: float typos

`ledger_dimensions.py` was not edited. The defect is that `_exact` calls
`sp.Rational(value)` for Python `float` and `sp.Float`, silently converting a
decimal typo to an exact value that is generally not the intended decimal
rational. `0.1` becoming
`3602879701896397/36028797018963968` is the sharp example.

The module-level repair would reject inexact numeric input before the
`sp.Rational` call, for example by raising `TypeError` when `value` is a
Python `float` or `sp.Float`, while continuing to accept integers,
`sp.Rational`, and deliberately exact strings. This is an API/strictness
change, not a stage-local refactor.

The live tree has seven importers: the six previously converted stages 004,
011, 012, 013, 016, and 018, plus the stage023 conversion now in the tree. A
module change would require an intentional update to
`scripts/ledger_dimensions.accepted.sha256`, regeneration of all seven Python
dimension sidecars because each embeds `ledger_dimensions_sha256`, and
re-execution/re-baselining of their comparison and canonical-table gates. A
search found no decimal-float dimension constructor call in those seven
stage files, so no current declaration edit was demonstrated, but the shared
contract and all seven provenance bindings would still change.

### Catalogued pre-existing assertion defects (not repaired)

- `scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py:82-85`:
  the retained sibling-shape `Dimension` arm of `assert_no_float` is entailed
  for every `Dimension` constructible through the current module because
  `_exact` has already removed every `Float`.
- `scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py:1099-1106`:
  an `AuditFailure` reaching the `__main__` handler bypasses `_record_fail`
  and its message, producing the observed zero-fail tally beside
  `OVERALL FAIL`.

Neither item was repaired.

## R5 — card I/O claim

Only the verification sentence in `paper/stages/stage_023.tex` changed:

```diff
-Dual-engine, independent symbolic routes, zero file I/O --- the source pair's scratch-YAML handoff is severed. The \texttt{.wl} is
+Dual-engine, independent symbolic routes --- the source pair's scratch-YAML handoff is severed; the \texttt{.wl} remains zero
+file I/O, while the \texttt{.py} writes its deterministic labelled-dimension sidecar. The \texttt{.wl} is
```

The sentence keeps the card's register, preserves the true severed-YAML
claim, and distinguishes the zero-I/O Wolfram engine from the Python
engine's one deterministic sidecar write. There is no second pair-level
zero-file-I/O claim elsewhere in the card.

The stage-artifact search found:

- `notes/stages/ledger_stage023_nullspace_underdetermination.md:804` still
  says both engines perform zero file I/O. It is now stale and was not edited
  because Markdown notes are prohibited in this round.
- `notes/stage023_pathA34_nullspace_underdetermination_source_map.md:294-300`
  contains the old reshape plan that ends in pair-level zero file I/O; it is
  likewise stale and prohibited. Its `.wl`-scoped statement at `:211-212`
  remains true.
- The Wolfram-only claims in
  `mathematica/ledger_stage023_nullspace_underdetermination_mathematica_audit.wl:3`,
  `:771`, and `:788`, and the corresponding `.out:6/:228` and committed
  transcript `:14/:206`, remain true because that engine writes no file.
- The Python source at `:4-5`, `:1055`, and `:1074` and its committed
  transcript at `:14` and `:199` already describe the dimension sidecar.

No other artifact was changed for R5.

## Round-2 generated sidecar and transcript

Acceptance command 1 refreshed the sidecar source binding and nothing else
in that artifact:

```diff
-source_sha256=ac218fbfaec569dc00370cfa332b6cd13ecf08d3a5ae26721c190b73fb4e0c66
+source_sha256=46f90b7c937d7d229727821db0ba5c080a93240ff506afad41ba156de7dde5c1
```

The latter equals the live stage source's `sha256sum`. The module digest
remained
`7cb8b7c93d189320f3ff857a665162e3290d1628e4511b2b387dd9a4f08e5bea`.
R4 does not change stdout, so the transcript was not regenerated.

The exact acceptance-5 comparison evidence was:

```text
ACCEPTANCE_5_EXIT=0
  207 12921 /tmp/stage023_round2_reference_body.txt
  207 12921 /tmp/stage023_round2_live.txt
2b8feafd7d663a188615a21450d00d4f645d3dd6acfd8ab621d4e0500df68dd1  /tmp/stage023_round2_reference_body.txt
2b8feafd7d663a188615a21450d00d4f645d3dd6acfd8ab621d4e0500df68dd1  /tmp/stage023_round2_live.txt
```

## Round-2 acceptance

All five commands were run from `research/pde_ledger_v2/scripts` under
`timeout 600`.

1. `timeout 600 python3 ledger_stage023_nullspace_underdetermination_sympy_audit.py`

   ```text
   TALLY sympy: 111 pass + 0 fail = 111 checks
   OVERALL PASS
   ACCEPTANCE_1_EXIT=0
   ```

2. `timeout 600 python3 compare_dimension_artifacts.py 023`

   ```text
   ARTIFACT_NAME_SET|stage=stage023|py=29|wl=29|shared=29|py_only=0|wl_only=0|source_coverage=not_checked
   RESULT|stage=stage023|status=PASS|mismatches=0
   ACCEPTANCE_2_EXIT=0
   ```

3. `timeout 600 python3 check_ledger_dimensions_pin.py`

   ```text
   MODULE_PIN_OK|module=scripts/ledger_dimensions.py|authority=scripts/ledger_dimensions.accepted.sha256|accepted_sha256=7cb8b7c93d189320f3ff857a665162e3290d1628e4511b2b387dd9a4f08e5bea|consumer_artifacts=(none)
   ACCEPTANCE_3_EXIT=0
   ```

4. `timeout 600 python3 generate_canonical_dimension_table.py`

   ```text
   WROTE|path=CANONICAL_DIMENSIONS.md|quantities=122|candidate_groups=6|needs_adjudication=3
   ACCEPTANCE_4_EXIT=0
   ```

5. `timeout 600 bash -c 'tail -n +9 output/ledger_stage023_nullspace_underdetermination_sympy_audit.txt | head -n -2'`,
   compared under `timeout 600` with the fresh capture from command 1:

   ```text
   ACCEPTANCE_5_EXIT=0
   ```

## Round-2 final repository state

`git status --short`:

```text
 M docs/development_pipeline.md
 M research/pde_ledger_v2/CANONICAL_DIMENSIONS.md
 M research/pde_ledger_v2/manifests/DIMENSION_REWRITE.md
 M research/pde_ledger_v2/notes/stages/ledger_stage023_nullspace_underdetermination.md
 M research/pde_ledger_v2/paper/stages/stage_023.tex
 M research/pde_ledger_v2/scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py
 M research/pde_ledger_v2/scripts/output/ledger_stage023_nullspace_underdetermination_sympy_audit.txt
?? research/pde_ledger_v2/scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.dimensions.txt
```

The development-pipeline, canonical-table, manifest, note, and Python
transcript entries were already dirty at round-2 entry. The generator was
required to run, but the canonical table's object hash stayed unchanged. The
round-2 final `git hash-object` values for the three non-report files whose
bytes this round changed are:

```text
26995ae877f092ba447cac7edd8300bfbe776439  research/pde_ledger_v2/scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py
c6c1a1d06129d41f74348dbf01457485d148323a  research/pde_ledger_v2/scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.dimensions.txt
26409de7ed605a5285166189c76a6c34ee5287a7  research/pde_ledger_v2/paper/stages/stage_023.tex
```

For proof that the two required-but-byte-preserving artifacts were not
modified by round 2:

```text
f263c1e9c73a7181123952c0b9c2773b1b8f5929  Python transcript and pre-round-2 backup
75568d995214e7e0d7551593587e4c7b6961edb7  CANONICAL_DIMENSIONS.md and pre-round-2 backup
```
