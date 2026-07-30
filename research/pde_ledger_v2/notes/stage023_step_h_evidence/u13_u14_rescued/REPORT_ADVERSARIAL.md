# stage023 dimension conversion — adversarial-with-ablation leg (fresh agent, 2026-07-28)

Artifact under review: `research/pde_ledger_v2/scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py`
(uncommitted working-tree change) + its sidecar `…_sympy_audit.dimensions.txt` + the by-label comparison
`python3 scripts/compare_dimension_artifacts.py 023` against the committed
`mathematica/out/ledger_stage023_nullspace_underdetermination_mathematica_audit.out`.

**Every mutation was made on a disposable copy I created** at
`research/pde_ledger_v2/_scratch/stage023_h/adversarial/`. No file under `scripts/` or `mathematica/`
was written. No `git checkout`/`restore`/`stash` was used. Live-file digests are in §7.

## 0. Method and controls

```bash
S=research/pde_ledger_v2/_scratch/stage023_h/adversarial
mkdir -p $S/pde_ledger_v2
cp -r research/pde_ledger_v2/scripts research/pde_ledger_v2/mathematica $S/pde_ledger_v2/
```
Driver `$S/run_ablation.sh <tag> <patch.py>`: fresh `cp -r` of `scripts/` into `$S/runs/<tag>/scripts`
+ symlink to `mathematica/`, apply the patch **on disk at the real path**, `rm -rf __pycache__`, run the
stage, then run the comparator from that copy. Every run therefore starts from a byte-clean tree, so
"restore between targets" is structural rather than asserted. Patches: `$S/patches/*.py`; raw captures:
`$S/runs/<tag>/{stage.out,cmp.out}`. Python 3.10.12; one stage run ~27 s.

**Baseline (unmutated copy).**
```
$ python3 ledger_stage023_nullspace_underdetermination_sympy_audit.py   # exit 0
TALLY sympy: 111 pass + 0 fail = 111 checks / OVERALL PASS
$ python3 compare_dimension_artifacts.py 023                            # exit 0
ARTIFACT_NAME_SET|stage=stage023|py=29|wl=29|shared=29|py_only=0|wl_only=0|source_coverage=not_checked
RESULT|stage=stage023|status=PASS|mismatches=0
```

**Positive control on the Mathematica side (my own addition).** I re-ran the committed `.wl` from my copy
(`timeout 600 math -script …_mathematica_audit.wl`, exit 0) and diffed the whole transcript against the
committed `.out`: **`FULL_DIFF_LINES=0`** — byte-identical, all 29 `DIM` records included. So the committed
transcript is genuine and current for this stage. That matters for how §2/target-9 should be read: the
*process* gap is real, the *artifact* is presently faithful.

## 1. What the 29 emitted records actually are

`dimension_records()` emits **22 `sourced_dims.*` + 7 `computed_dims.*`**.

- The 22 `sourced_dims.*` are the literal `SOURCED_DIMS` table verbatim. Nothing computes them.
- The 7 `computed_dims.*` are genuinely produced by `dim_of` walking the Mul/Pow/Add tree — but **5 of the
  7** (`T0`, `T1`, `epsilon0`, `epsilon1`, `P0_physical`) have expected value `ZERO_DIM`, i.e. five of the
  seven teeth assert *dimensionless == dimensionless*. `dim_of` also returns `ZERO_DIM` for any purely
  numeric expression, so a computation that silently degenerated to a number still passes those five.
- The Mathematica side declares the **same literal table** (`baseDims` in the `.wl`) and computes the same 7
  through its own `dimOf`. The comparator therefore establishes **agreement between two hand-written
  transcriptions**, not derivation.

## 2. The ablation target list (orchestrator-owned)

Exit codes below are literal. "Fires" = the first `FAIL …` line printed by the stage.

### Target 1 — the 7 `[name] matches its sourced expected dimension` assertions

**(a) Isolated form — is each tooth individually able-to-fail?** Patch inserts, right after
`dim = baseline["dimension"]`, a local perturbation of that name's *expected* entry only
(`dim["expected"][<name>] = Dim(9,9,9)`), leaving `dimensional_ok` as computed.

| name | stage exit | first FAIL | checks reached | comparator |
|---|---|---|---|---|
| A0 | 1 | `[A0] matches its sourced expected dimension` | 43P+1F | exit 1, stale-sidecar |
| A1 | 1 | `[A1] …` | 44P+1F | exit 1 |
| T0 | 1 | `[T0] …` | 45P+1F | exit 1 |
| T1 | 1 | `[T1] …` | 46P+1F | exit 1 |
| epsilon0 | 1 | `[epsilon0] …` | 47P+1F | exit 1 |
| epsilon1 | 1 | `[epsilon1] …` | 48P+1F | exit 1 |
| P0_physical | 1 | `[P0_physical] …` | 49P+1F | exit 1 |

**All 7 are able-to-fail at their own assert.** OK (`$S/runs/t1_*`)

**(b) Realistic form — corrupt the source instead.** `EXPECTED_DIMS["A0"] = Dim(9,9,9)` at module level:
stage exit 1, but the assertion that fires is **`selector control reaches CROSS_L_RESIDUAL_PREDICTION`**
(22P+1F, ~60 checks *earlier*), not the target-1 tooth. `dimensional` is rung 3 of the `base_verdict`
ladder, so any genuine dimensional corruption flips the verdict before `run_dimensional_gate()` runs.
**Different-assertion outcome.** (`$S/runs/t1_real_A0`)

### Target 2 — `baseline dimensional_ok is True` — CANNOT FAIL AS SHIPPED

`dimensional_ok = (not errors) and all(computed == expected)`.
- The `all(...)` leg is entailed by the 7 preceding teeth: if it were false, one of them fired first.
- The `errors` leg is unreachable at this assert: if `dim_of` raised for name *X*, `computed` lacks *X* and
  line 909's `dim['computed'][name]` raises `KeyError` first.

Empirically: injecting an undeclared symbol into `P0_physical` (a real `DimError`) gives stage exit 1
firing **`selector control reaches CROSS_L_RESIDUAL_PREDICTION`** at 22P+1F — neither target 2 nor the
`KeyError`, because the verdict ladder intercepts even earlier (`$S/runs/t2_real_err`). Forcing the flag
directly (`dim["dimensional_ok"] = False`) does fire it (exit 1, 50P+1F, `$S/runs/t2_iso`) — but that
mutation has no source-level pre-image. **Target 2 is a restatement of target 1, not an independent tooth.**

### Target 3 — `corrupt sourced [M0]+=L reaches FAIL_DIMENSIONAL` — able-to-fail

Neutralising the mutation (`dims[M0] * Dim(1,0,0)` -> `* Dim(0,0,0)`): stage exit 1, fires **at its own
assert**, 51P+1F. (`$S/runs/t3`)

### Target 4 — `corrupt free q_free dimension reaches NO_FAIL locally` — VACUOUS
### Target 5 — `q_free appears in no checked expression` — VACUOUS

`q_free` occurs in exactly four places in the whole file: its `sp.Symbol` declaration (l.172), its
`SOURCED_DIMS` entry (l.482), the mutation branch (l.509), the mention-probe (l.535) — plus its emitted
record (l.563). **It appears in no physics expression anywhere.** Consequences:

- **Neutralising the mutation entirely** (`dims[q_free] = Dim(7,0,0)` -> `Dim(0,0,0)`, i.e. corrupting
  nothing at all): **stage exit 0, 111 pass + 0 fail, comparator exit 0, `mismatches=0`, sidecar
  byte-identical.** The control does not detect the absence of the corruption it is named for.
  (`$S/runs/t4_neutral`)
- Making `q_free` enter one checked expression (`"T0": transfers["T0"] * (1 + q_free)`) *does* fire target 4
  (exit 1, 52P+1F) — so the assertions are able-to-fail *in principle*, and only in a counterfactual where
  the stage is written differently. (`$S/runs/t45_enter`)

Target 5 is the premise; target 4 is its consequence. Neither constrains the dimension gate. The pair
establishes only that an unused symbol is unused.

### Target 6 — the `3f corrupt sourced dimension` failure-ablation block

**The description does not match the artifact: `check_failure_ablation` contains FIVE `expect_*` calls, not
four** (`.py:773-788`). Per-sub-assertion:

1. `… dynamic rerun with mutation reaches FAIL_DIMENSIONAL` — **shadowed.** For the 3f block this asserts
   the same fact as the line-915 tooth three lines earlier, so it can never fire first: my target-3
   mutation fired at 915. Removing 915 *and* neutralising the mutation makes it fire (exit 1, 53P+1F) — so
   it is able-to-fail, but **redundant as shipped**. (`$S/runs/t6_3f_shadow`)
2. `… dynamic rerun without mutation returns earned native-nullspace FAIL` — able-to-fail. Setting
   `without_context = run_gate(mutated)` fires it (exit 1, 34P+1F). (`$S/runs/t6_same_trace`)
3. `… two computed verdicts differ` and 4. `… own expected token fires` — **entailed by 1 AND 2 and cannot
   fail independently.** Given `trace[0]==expected` (assert 1) and `trace[1]==FAIL_UNDERDETERMINED` with
   `expected != FAIL_UNDERDETERMINED` (assert 2), `trace[0]!=trace[1]` and
   `trace[0]==expected and trace[1]!=expected` both follow. The t6_same_trace experiment confirms the
   ordering: with a degenerate trace, assert 2 fires, never 3 or 4.
5. `… independently rerun neutralized mutation does not fire` — able-to-fail on its *identity* leg
   (blanking `inert_rerun_token` fires it, exit 1, 37P+1F, `$S/runs/t6_token`). Its *substantive* leg
   (`not neutralized_mutation_fires`) is by construction: the neutral mutation is physics-identical to the
   baseline, so it cannot reach the expected FAIL.

Net: of five sub-assertions, **one is shadowed, two are entailed, two are genuinely able-to-fail**.

### Target 7 — the `EXPECTED_DIMS[name] if <expr> == 0 else dim_of(...)` branch (`.py:525`)

Replacing it with an unconditional `dim_of(expr, dims)`: **stage exit 1 at 90P+1F, firing
`3d perfect_return dynamic rerun with mutation reaches FAIL_OVERCANCEL`** — and the 7 baseline dimension
teeth passed on the way there. (`$S/runs/t7_nobranch`)

Reading: the branch is **not** exercised by the baseline (so the 7 emitted `computed_dims` are genuinely
computed, not shortcut), but it *is* load-bearing under the `perfect_return` mutation, where `Z->0` makes
`epsilon0`, `epsilon1`, `A0`, `A1` identically zero. In that branch the dimension of those four is **set
equal to the expected value by construction**, keeping `dimensional_ok` True so the ladder can report
`FAIL_OVERCANCEL` instead of `FAIL_DIMENSIONAL`. It is a deliberate, commented rule ("an identically zero
amplitude still occupies its independently sourced ledger slot"), and it is sound for the emitted artifact
— but it means the `dimensional` rung of the verdict ladder is *unfalsifiable* for any mutation that zeroes
an amplitude. Worth stating in the note; not a defect in the 29 records.

### Target 8 — the exactness/float guards on the dimension path — BOTH VACUOUS

- **The `Dimension` branch added to `assert_no_float` (`.py:83-86`) is dead code.** Replacing its body with
  `raise AuditFailure("ABLATION: … WAS reached")`: **stage exit 0, 111 pass + 0 fail, comparator exit 0.**
  It is never executed in a full run — `expect_bool` collapses the `Dimension` comparison to a `bool`
  before any residual is built, so no `Dimension` ever reaches `assert_no_float`. (`$S/runs/t8_deadcode`)
- **It could not fire even if it were reached.** `ledger_dimensions._exact()` is `sp.Rational(value)`,
  which converts any `float` to an exact `Rational`. A `Dimension` cannot hold a `Float` atom.
- **Consequence — an inexact author typo is invisible.** `R0: Dim(0.1, 1, -1)`: **stage exit 0, 111
  pass + 0 fail**, and the sidecar emits `exponents={3602879701896397/36028797018963968, 1, -1}`. Nothing
  in the stage objects; only the `.wl` value comparison catches it. (`$S/runs/i6_floatlit`)
- **And genuine `Float` exponents survive the comparator.** Patching the shared module so `_exact` returns
  `sp.Float(value)`: stage exit 0, 111 pass + 0 fail; sidecar reads `exponents={1.00000000000000, 0.0, 0.0}`.
  The comparator rejects it *only* via the module pin (`CONTROL_FAILURE: MODULE_PIN_MISMATCH`, exit 1).
  After `python3 check_ledger_dimensions_pin.py --accept`, the same float-bearing sidecar gives
  **`RESULT|stage=stage023|status=PASS|mismatches=0`, exit 0** — because `Fraction("1.00000000000000")`
  parses cleanly. **There is no exactness check anywhere on the sidecar text**; the pin is a
  change-detector, not an exactness gate. (`$S/runs/t8_float`)

### Target 9 — the sidecar's source-hash binding and the shared-module pin — mostly real

- **Post-run source edit, no re-run.** Green run, then append one comment line to the `.py`: comparator
  exit 1, `FAIL: stale Python dimension sidecar …: asserted source_sha256=10f1253…, computed
  source_sha256=2c19c80…`. Restoring the byte returns it to `PASS`.
- **Failing run leaves no fresh sidecar.** In all 13 exit-1 ablations the sidecar was not rewritten and the
  comparator reported staleness — the two controls compose correctly.
- **Header digests are mandatory.** Stripping `|source_sha256=…` -> `FAIL: … has no source_sha256
  assertion`, exit 1. Stripping both digests -> same.
- **Module pin is real** (see target 8) — but note it fires *before* the freshness check and masks it.
- **Asymmetry:** `require_fresh_python_sidecar` is applied to the Python sidecar only. The `.out` has **no
  `source_sha256`, no `wl_sha256`, no freshness check, and no manifest entry** (there is no
  `manifests/stages/stage023.json`; a repo-wide grep finds no artifact binding this `.out` to its `.wl`).
  The comparator prints `source_coverage=not_checked` in its own headline. My re-run of the `.wl` shows the
  transcript is currently faithful, so this is an **unpinned control, not a broken one**.

### Target 10 — targets I add, and description mismatches

**Description mismatches to record:** target 6 says "four sub-assertions"; there are five (see Target 6).
All other targets exist as described.

**Added targets, each executed:**

- **10a — An emitted record need not come from the computation its name claims.** Replacing
  `"computed_dims.T0": computed["T0"]` with the literal `ZERO_DIM`: **stage exit 0, 111 pass, comparator
  exit 0, `mismatches=0`.** Repointing `"computed_dims.epsilon0"` at `computed["T1"]`: identical result.
  Because 5 of the 7 computed records and 4 of the 22 sourced records share the all-zero triple, **any
  substitution within a value-class is invisible to both engines and to the comparator.**
  (`$S/runs/i1_hardcode`, `$S/runs/i2_swap`)
- **10b — 20 of the 22 `sourced_dims` records are pinned by nothing in the Python.** The stage's
  dimensional teeth constrain only *relations*, so whole gauge directions are free. One joint mutation
  changing **16 declarations at once** — `a=(2,0,0)`, `c_s=(2,0,-3)`, `omega=(0,0,-3)`, `D0=(0,-4,-6)`,
  `K0c=K_eta=T_Omega=Z0_ret=Z1_ret=(3,3,3)`, `Omega_U=Omega_W=(0,1,0)`, `R_mix=(0,2,0)`, `g_U=g_W=(0,0,0)`,
  `R0=(11,11,11)`, `R1=(12,12,12)` — gives **stage exit 0, 111 pass + 0 fail**, and **all 7 `computed_dims`
  records emerge byte-identical to the pristine ones**. The comparator reports `mismatches=16`. Adding the
  four dimensionally inert symbols (`eta_null`, `gain0`, `gain1`, `q_free`) gives 20. Only `M0` and `D1`
  are pinned by any Python assertion — and only *relative* to the declared `EXPECTED_DIMS["A0"]`/`["A1"]`,
  which are themselves declarations. **The sole control on 20 of the 22 sourced values, and on the absolute
  scale of the other two, is the `.wl`'s identical literal `baseDims` table.** (`$S/runs/i7_joint`, plus
  `i3_gauge` 5-way, `i5_port` 3-way, `i4_R0` 1-way, all green at 111 pass)
- **10c — The comparator agrees, it does not verify.** Same wrong value planted in *both* engines
  (`sourced_dims.R0 = (5,5,5)` in the `.py` and in my copy of the `.out`): **stage exit 0, comparator exit
  0, `py=29|wl=29|shared=29|mismatches=0`, `status=PASS`.** (`$S/runs/h4_agreed`)
- **10d — There is still no coverage floor.** Deleting 28 of the 29 records from `dimension_records()`
  *and* from my copy of the `.out`: **comparator exit 0,
  `ARTIFACT_NAME_SET|stage=stage023|py=1|wl=1|shared=1|py_only=0|wl_only=0` … `status=PASS`.** The
  `compared=0` hole recorded in `_scratch/COMPARATOR_FLOOR_FIX.md` was closed (`len(shared_names)==0` ->
  FAIL), but `compared=1 of 29` still passes: there is no minimum count and no required-name manifest.
  Nothing in the stage pins how many records it must emit. (`$S/runs/h5_floor`)

## 3. Ranked findings

1. **20 of the 22 emitted `sourced_dims` values are unearned by the Python** (10b). A 16-way joint
   redefinition leaves the stage at 111/111 with all 7 computed dims unchanged. The only control is the
   `.wl`'s identical literal table — a transcription cross-check, not two independent derivations.
2. **A record's value need not come from the computation its name claims** (10a). Hardcoding or
   cross-wiring any of the 25 records that share a triple with a sibling is invisible to both engines and
   the comparator.
3. **Targets 4 and 5 are vacuous**: `q_free` enters no expression, so removing the corruption entirely
   leaves the run green and the sidecar byte-identical.
4. **Target 2 cannot fail independently** of target 1; and any *source-level* dimensional corruption fires
   `selector control reaches CROSS_L_RESIDUAL_PREDICTION` ~60 checks earlier, so the dimensional gate's own
   assertions never do the detecting in practice.
5. **Both float/exactness guards on the dimension path are vacuous** (target 8): the `Dimension` branch of
   `assert_no_float` is dead code, cannot fire by construction, and float-valued exponents pass the
   comparator's `Fraction` parse.
6. **Two of five sub-assertions in `check_failure_ablation` are entailed; a third is shadowed** for the 3f
   block by the line-915 tooth.
7. **The comparator establishes agreement, not correctness** (10c), has no coverage floor above 1 record
   (10d), and prints `source_coverage=not_checked` itself.
8. **The Mathematica side is unpinned** — no `source_sha256`, no manifest, no freshness check (target 9).
   Mitigating: I re-ran the `.wl` and the committed `.out` reproduces byte-identically.
9. **The `==0 -> EXPECTED_DIMS` branch makes the `dimensional` verdict rung unfalsifiable** for any
   amplitude-zeroing mutation (target 7). Deliberate and documented; does not taint the emitted 29.

**Not found / UNDETERMINED.** I could not construct any mutation that changes an emitted `computed_dims.*`
value while leaving both the stage and the comparator green — every attempt either fired a named assertion
or produced a comparator `MISMATCH`. I could not defeat the sidecar `source_sha256` binding, the
header-digest requirement, or the module pin; all three are genuinely able-to-fail. I did not attempt to
verify that the *values* in the shared literal table are physically correct against stage008 / stage017 /
pathA_34 — that is the fidelity leg's territory, and it is where finding 1 lands the burden.

## 4. Agreement with the other instrument on this stage

I reached §2 and §3 before reading `_scratch/stage023_h/ABLATION_SUMMARY.md` (the orchestrator's
supplementary 51-row ablation) and `_scratch/COMPARATOR_FLOOR_FIX.md`. **No disagreements.** Independently
reproduced, before reading: the six stage-invisible declarations (`R0`, `R1`, `eta_null`, `gain0`, `gain1`,
`q_free`); that corrupted declarations fire `selector control reaches CROSS_L_RESIDUAL_PREDICTION` rather
than the dimensional teeth; that mis-binding an emitted record leaves the stage at 111 pass; and that
within-triple mis-bindings are invisible to all three instruments.

**What I add beyond it:** its axis 1 mutates one declaration at a time and therefore classifies 16 as
"caught by the stage". Finding 10b shows that catch is only of *relative* inconsistency — under a joint
perturbation 14 of those 16 are free as well, so the stage-unpinned count is **20 of 22, not 6**. Also new
here: the vacuity of both float guards, the entailed/shadowed sub-assertions, the `==0` branch's effect on
the `dimensional` rung, the agreed-wrong-value experiment, the surviving `compared=1 of 29` floor gap, and
the `.wl` re-run positive control.

## 5. Suggested repairs (not applied — this is a review instrument)

- Replace the `q_free` pair with a control whose corrupted symbol actually enters a checked expression.
- Add an explicit emitted-record manifest (name -> required count) so `compared=1 of 29` fails, and assert
  the count inside the stage.
- Give the `.out` a `wl_sha256`/manifest binding symmetric with the Python `source_sha256`.
- Either exercise `assert_no_float` on `Dimension` objects (and reject `float` inputs in `_exact` rather
  than silently rationalising them), or delete the dead branch and stop claiming a float guard.
- Reorder or duplicate the dimensional teeth ahead of the verdict-ladder assertions, so a corrupted
  declaration fires the dimensional assertion rather than the selector-control one.
- Drop sub-assertions 3 and 4 of `check_failure_ablation`, or restate them so they are not entailed.

## 6. Artifacts

`$S/run_ablation.sh` · `$S/patches/*.py` (24 patches) · `$S/runs/<tag>/{stage.out,cmp.out,scripts/}` for
tags `t1_{A0,A1,T0,T1,epsilon0,epsilon1,P0_physical}`, `t1_real_A0`, `t2_iso`, `t2_real_err`, `t3`,
`t4_neutral`, `t45_enter`, `t6_same_trace`, `t6_token`, `t6_3f_shadow`, `t7_nobranch`, `t8_deadcode`,
`t8_float`, `i1_hardcode`, `i2_swap`, `i3_gauge`, `i4_R0`, `i5_port`, `i6_floatlit`, `i7_joint`,
`h_tests`, `h4_agreed`, `h5_floor`.

## 7. Live-tree integrity

`git status --porcelain` outside `_scratch/` is unchanged from the start of this session (the four
pre-existing modifications plus the untracked sidecar). `diff -r --exclude=__pycache__` between
`research/pde_ledger_v2/scripts` and my copy: identical; between `research/pde_ledger_v2/mathematica` and
my copy: identical.

```
10f125307d0a3a0bd09b19d32e0f4f66ce6512a4fc417f58ac2fca11c670edc4  research/pde_ledger_v2/scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py
7cb8b7c93d189320f3ff857a665162e3290d1628e4511b2b387dd9a4f08e5bea  research/pde_ledger_v2/scripts/ledger_dimensions.py
87ae582d14802e9c23512be2a77bde99a180603e7df147a2849eeb3d600ebce1  research/pde_ledger_v2/scripts/compare_dimension_artifacts.py
cba3a9de579073d60fec8d11268cb5a626017349591290f0d5d274b7f637f096  research/pde_ledger_v2/scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.dimensions.txt
fae1d5cdc4c364857780854fd26d469edc4dfd4a61ab03583f3dd51d67212759  research/pde_ledger_v2/mathematica/out/ledger_stage023_nullspace_underdetermination_mathematica_audit.out
```
