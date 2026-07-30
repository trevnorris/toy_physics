# Prior-art survey — would an off-the-shelf Python mutation-testing tool replace the hand-rolled ablation harness?

Read-only survey. Date: 2026-07-28. Nothing was installed; no project file was modified except this one.

---

## 0. What the replacement would have to do (established from the artifacts, not assumed)

### 0.1 The target and its oracle

Target: `/var/projects/toy_physics/research/pde_ledger_v2/scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py`
— 1114 lines, standalone, **no pytest**. It is executed as a script (`if __name__ == "__main__":`) and terminates
with `raise SystemExit(0)` / `raise SystemExit(1)` (`:1096-1114`). It prints `PASS  <name>` / `FAIL  <name>` lines
(`_record_pass` / `_record_fail`, `:109-142`) and a `TALLY sympy: N pass + M fail` line.

Side effect: on the success path only (`:1089-1092`, i.e. *after* every check, just before `return 0`) it calls

```python
emit_dimension_sidecar(__file__, dimension_records(run_gate(Mutation())["dimension"]))
```

which writes `<stem>.dimensions.txt` **beside the source** (`ledger_dimensions.py:212-221`,
`destination = source.with_suffix(".dimensions.txt")`) with a header that embeds

```
source_sha256            = sha256(the .py file on disk, read at runtime)
ledger_dimensions_sha256 = sha256(ledger_dimensions.py)
```

External checker: `/var/projects/toy_physics/research/pde_ledger_v2/scripts/compare_dimension_artifacts.py 023`.
It (a) pins `ledger_dimensions.py`'s bytes against the separately-held authority file
`scripts/ledger_dimensions.accepted.sha256` via `check_ledger_dimensions_pin` (`:15`, `:251-256`);
(b) **re-computes** `sha256` of the on-disk `.py` and requires it to equal the digest the sidecar asserted
(`require_fresh_stage_source`, `:204-213`); (c) compares the 29 `DIM|` records against a frozen
`mathematica/out/ledger_stage023_*_mathematica_audit.out` and prints `MISMATCH <name>: py=…; wl=…`.
It locates the sidecar by globbing **`SCRIPT_DIR`** (`parse_stage`, `:60-64`) — i.e. the directory containing
`compare_dimension_artifacts.py` itself.

So the "test suite" is: `bash -c 'python3 <script> && python3 compare_dimension_artifacts.py 023'`
plus an out-of-band requirement that the multiset of `PASS` lines be unchanged.

### 0.2 The two ablation experiments actually run

**(a) `ablation/results.tsv` — value ablation, 22 rows.**
Columns: `var  stage_exit  sidecar_changed  records_moved  cmp_exit  mismatch_names`.
One row per entry of the `SOURCED_DIMS` dict literal (`:460-483`), e.g. `R0: Dim(0, 1, -1)`.
The recorded outcome is *not* killed/survived — it is a five-column record including "did the named emitted
record move" and "which name mismatched". Six rows (`R0 R1 eta_null gain0 gain1 q_free`) give
`stage_exit=0` and a genuine `MISMATCH sourced_dims.<name>`; the other sixteen give `stage_exit=1`.

**(b) `ablation/bind.tsv` — binding ablation, 29 rows.**
Columns: `record  stage_exit  cmp_exit  caught_by_name`. One row per entry of `dimension_records()`
(`:541-571`), each rebound to a *different symbol*. `b2_cmp.txt` shows the shape:
`MISMATCH sourced_dims.Z1_ret: py={L=1, M=0, T=0}; wl={L=0, M=1, T=-2}` — `{1,0,0}` is the dimension of `a`,
i.e. `"sourced_dims.Z1_ret": SOURCED_DIMS[Z1ret]` was rebound to `SOURCED_DIMS[a]`.
**This is an identifier-for-identifier swap, not a literal perturbation.**

### 0.3 ⭐ The failure mode that actually bit: the side-effect file

`cmp_a.txt` (an artifact of the real run) reads:

```
COMPARISON_SKIPPED: Python dimension sidecar freshness failed
FAIL: stale Python dimension sidecar …: asserted source_sha256=d4f093bc…, computed source_sha256=eb29e632… from …py
RESULT|stage=stage023|status=FAIL|mismatches=not_checked
```

Sixteen of the twenty-two `cmp_*.txt` files are 673 bytes — the *freshness-failure* text
(`cmp_a cmp_c_s cmp_D0 cmp_D1 cmp_gU cmp_gW cmp_K0c cmp_Keta cmp_M0 cmp_omega cmp_OmegaU cmp_OmegaW
cmp_Rmix cmp_TOmega cmp_Z0ret cmp_Z1ret`). Only six are 315–320 bytes — the *real `MISMATCH`* text.

⇒ **For 16 of 22 mutants the comparator's verdict was produced by the residual `.dimensions.txt` file, not by a
dimension comparison at all.** The sidecar is a persistent side effect that survives restoration of the source
(the script never reaches its `emit_` call when a check fires), so unless it is explicitly reset between mutants
the comparator reports "stale sidecar" instead of the finding you were looking for. Any tool adopted here must
reset an *emitted artifact*, not just the source file.

Snapshot/restore in the hand harness: `ablation/GREEN.py` (46 389 B) and `ablation/GREEN.sidecar` are the
byte-exact baselines; `docs/development_pipeline.md:19` mandates *"Every one restores by `cp` and proves the
restore with `git hash-object`"* and ⛔ *"Never `git checkout`/`restore`/`stash` to undo an ablation on
uncommitted work"*.

### 0.4 The process constraints a tool must not violate

- `docs/development_pipeline.md:330-332` — *"⛔ ablation target list owned by the **orchestrator** (per-tooth
  over every able-to-fail check / emitted record — never chosen by the builder)"*.
- `docs/development_pipeline.md:19` — *"⛔ The ablation target list is the **orchestrator's**, never the
  builder's"*.
- Phase 4(c), `:77` — *"Ablation is per-tooth over every able-to-fail check or emitted record."*

A tool that **selects its own mutation sites** (or samples them randomly) does not satisfy this clause even if
it is technically capable; the human target list is a *control*, not a convenience.

---

## 1. Availability — what is installed here

| tool | `import` result |
|---|---|
| `cosmic_ray` | `ModuleNotFoundError: No module named 'cosmic_ray'` |
| `mutmut` | `ModuleNotFoundError: No module named 'mutmut'` |
| `mutpy` | `ModuleNotFoundError: No module named 'mutpy'` |
| `mutatest` | not installed (absent from `pip list`) |
| `universalmutator` | not installed (absent from `pip list`) |

`python3 -m pip list | grep -i mut` returns exactly one line: `mutagen 1.47.0` — that is the **audio-metadata**
library, unrelated to mutation testing.

Local interpreter: **Python 3.10.12**, sympy 1.14.0.

PyPI **is** reachable from this box: a `pip download cosmic-ray --no-deps` resolved and fetched
`cosmic_ray-8.4.6-py3-none-any.whl` before failing on the deliberately-invalid output path. Installing was
therefore possible and was **not** done, per the survey's read-only constraint.

**Verified locally:** module availability; Python version; PyPI reachability; and the `__pycache__` behaviour in
§5 below.
**Not verified locally:** every behavioural claim about every tool. All of §2–§6 is from published
documentation and from reading each project's source on GitHub, never from execution.

---

## 2. cosmic-ray

Latest **8.4.6**, released 2026-04-02; declares Python ≥3.9 and lists 3.9–3.13 — compatible with the local 3.10.
Actively maintained (8.4.4 in Feb 2026).

### 2.1 Test-runner fit — **YES, with one wrinkle**

Config is TOML. The canonical example from the tutorial:

```toml
[cosmic-ray]
module-path = "mod.py"
timeout = 10.0
excluded-modules = []
test-command = "python -m unittest test_mod.py"

[cosmic-ray.distributor]
name = "local"
```

`test-command` is an **arbitrary command line** and pass/fail is **purely the exit code**. From
`src/cosmic_ray/testing.py`: the command is executed with `subprocess.run()`, exit status 0 ⇒
`TestOutcome.SURVIVED`, non-zero ⇒ `TestOutcome.KILLED`, `subprocess.TimeoutExpired` ⇒ `KILLED`, a launch
exception ⇒ `INCOMPETENT`. It sets `PYTHONDONTWRITEBYTECODE=1`.

⚠ **Wrinkle, and it is citable:** the command is split with **`shlex.split()`, not run through a shell**. So
`a && b`, pipes and redirects are *not* interpreted. The two-leg oracle must be wrapped:

```toml
test-command = "bash -c 'cd /var/projects/toy_physics/research/pde_ledger_v2 && ./_scratch/ablation_oracle.sh'"
```

That wrapper is also where the `PASS`-multiset comparison and the sidecar reset have to live.

### 2.2 Targeting — **NO include-list; filters are subtractive only**

The three shipped filters (`pyproject.toml` console-scripts: `cr-filter-operators`, `cr-filter-pragma`,
`cr-filter-git`) all **skip** work items:

- `cr-filter-operators` — regex over *operator names*, `[cosmic-ray.filters.operators-filter] exclude-operators = [...]`.
- `cr-filter-pragma` — skips any mutation on a line carrying `# pragma: no mutate`.
- `cr-filter-git` — skips any mutation on a line not changed relative to `[cosmic-ray.filters.git-filter] branch`.

There is **no "mutate only these lines" input.** Restricting to the ~22 dict entries via `cr-filter-pragma`
means annotating everything else — a denylist over a 1114-line file, which is the shape this project has a
standing rule against.

The supported route is a **custom filter**. The docs say only *"these filters are nothing more than simple
programs that modify a session in some way; it should be straightforward to write your own"* — they do **not**
specify an API. The mechanics are however plainly available: the session is a **plain SQLite file**, and
`src/cosmic_ray/work_item.py` shows `WorkItem(job_id, mutations: tuple[MutationSpec])` where each
`MutationSpec` carries `module_path`, `start_pos=(line,col)`, `end_pos=(line,col)`, `operator_name`,
`occurrence`, `operator_args`. An include-filter that marks everything outside a hand-written line set as
`WorkerOutcome.SKIPPED` is therefore writable — but it is code *you* write and maintain.

**Scale of the knock-out.** Measured on the target file by AST walk: 172 numeric literals, 52 `Compare` nodes,
62 boolean literals, 124 `BinOp` nodes, 1794 `Name` nodes. `NumberReplacer` alone emits **two** work items per
numeric literal (`OFFSETS = [+1, -1]`, deterministic, `eval(node.value) + OFFSETS[index]`). `cosmic-ray init`
would therefore produce order-10³ work items where the human wants 22 + 29 = 51. The filter must delete ~95 %.

### 2.3 ⛔ The bind ablation is **not expressible** by any shipped operator

The core provider (`src/cosmic_ray/operators/provider.py`) registers exactly:
`binary_operator_replacement.operators()`, `comparison_operator_replacement.operators()`,
`unary_operator_replacement.operators()`, `AddNot`, `ReplaceTrueWithFalse`, `ReplaceFalseWithTrue`,
`ReplaceAndWithOr`, `ReplaceOrWithAnd`, `ReplaceBreakWithContinue`, `ReplaceContinueWithBreak`,
`ExceptionReplacer`, `NumberReplacer`, `RemoveDecorator`, `ZeroIterationForLoop`, `VariableReplacer`,
`VariableInserter`, plus `NoOp` (name-accessible, excluded from iteration).

`VariableReplacer` sounds like the one you want and **is not**. From
`src/cosmic_ray/operators/variable_replacer.py`: it takes `Argument("cause_variable", …)` /
`Argument("effect_variable", …)` and, on a `Leaf` matching `cause_variable`, substitutes

```python
Number(start_pos=node.start_pos, value=str(randint(-100, 100)))
```

— a **random integer**, not another identifier, and **non-deterministic across runs**. Both properties are
disqualifying here: `SOURCED_DIMS[Z1ret] → SOURCED_DIMS[a]` cannot be produced, and a random replacement value
cannot be re-run to reproduce a banked TSV row.

Writing a custom operator is possible — `pyproject.toml` declares
`[project.entry-points."cosmic_ray.operator_providers"] core = "cosmic_ray.operators.provider:OperatorProvider"`,
loaded through stevedore, so third-party providers are a real extension point — but registering one means
**packaging and installing a distribution**, not dropping a file in `scripts/`. The base class
(`operators/operator.py`) requires `mutation_positions(node)`, `mutate(node, index)` and the `examples()`
classmethod.

### 2.4 The four failure modes

| | cosmic-ray |
|---|---|
| **snapshot / restore of the source** | ✅ **and architecturally correct here.** `src/cosmic_ray/mutating.py`'s `use_mutation` is *"A context manager that applies a mutation for the duration of a with-block. This applies a mutation to a file on disk, and after the with-block it put the unmutated code back in place."* Because it mutates **the real file at the real path**, `__file__`, the sidecar's destination path, and the `sha256(source.read_bytes())` freshness pin all behave exactly as they do in the hand harness. This is the single biggest compatibility fact in the survey. |
| **per-mutant result capture** | ⚠ partial. `WorkResult` has `worker_outcome`, `test_outcome`, `output` (captured worker output) and `diff`, all persisted in `session.sqlite`. So the wrapper can print the five TSV columns to stdout and they *are* stored — but `cr-report`/`cr-html`/`cr-rate` render survival rates, not your columns. Reconstructing `results.tsv` means a bespoke SQLite reader. |
| **resume after interrupt** | ✅ first class. *"Cosmic Ray can safely stop in the middle of a (potentially very long) session and be restarted. Since the session knows which work is already completed, it can continue where it left off."* |
| **⭐ side-effect files the script WRITES** | ❌ **not handled.** `use_mutation` restores *only* `module_path` — "it put the unmutated **code** back in place". There is **no pre-mutant / post-mutant hook** in the config schema or the CLI. The emitted `.dimensions.txt` is untouched, so cosmic-ray reproduces exactly the contamination of §0.3 unless the `test-command` wrapper itself does `rm -f …dimensions.txt` (or `cp GREEN.sidecar`) as its first action. Workable, but the tool contributes nothing. |

### 2.5 Verdict — **PARTIAL FIT**

Genuinely good: arbitrary-exit-code test command, on-disk mutate-and-restore at the real path (which is what
keeps the sha256 freshness pin meaningful), durable resumable session, per-mutant stdout retained.

Missing, each requiring code you write:
1. include-list targeting (custom filter over `session.sqlite`);
2. the identifier-rebinding operator for `bind.tsv` (custom operator + stevedore entry-point plugin);
3. side-effect (`.dimensions.txt`) reset between mutants (in the wrapper);
4. multi-column result extraction (bespoke reader over `WorkResult.output`);
5. `bash -c` wrapping because `test-command` is `shlex.split`, not shelled.

And it still does not satisfy `development_pipeline.md:330-332` on its own: `cosmic-ray init` enumerates the
targets, the orchestrator does not.

---

## 3. mutmut

Latest **3.6.0**, 2026-06-06; requires Python ≥3.10. Actively maintained.

### 3.1 Test-runner fit — ⛔ **NO. pytest only.**

`HISTORY.rst` 3.0.0: *"Execution model switched to mutation schemata, which enabled parallel execution"*, new
terminal UI, and the framework note **"Pytest only"**. 3.6.0 renamed `paths_to_mutate` → `source_paths` and
deprecated `tests_dir` in favour of `pytest_add_cli_args_test_selection`.

The whole documented config surface is pytest-shaped: `source_paths`, `also_copy`, `max_stack_depth`,
`only_mutate`, `do_not_mutate`, `mutate_only_covered_lines`, `type_check_command`, `debug`,
`use_setproctitle`, `do_not_mutate_patterns`, `cache_invalidation_files`, `cache_invalidation_exclude`,
`on_dependency_change`, `use_git_change_detection`, `timeout_constant`, `timeout_multiplier`,
`pytest_add_cli_args`, `pytest_add_cli_args_test_selection`. **There is no `runner` / `test-command` key.**

GitHub issue [boxed/mutmut#373](https://github.com/boxed/mutmut/issues/373) (2025-03-24), *"How to specify
runner, create JUnit XML in mutmut 3?"*, states: *"`--no-progress` (which is not crucial), `--runner`, and
`junitxml` are all gone from mutmut 3, and there doesn't seem to be configuration options to achieve the same
thing."*
⚠ **Caveat on evidence strength:** that is the *reporter's* claim. I could **not** retrieve a maintainer reply
on that issue. The independent corroboration is `HISTORY.rst`'s own "Pytest only" and the absence of any runner
key from the documented config surface.

The target has no pytest suite, and cannot trivially become one: its module scope ends in
`raise SystemExit(1)` on failure and the whole audit runs at import time.

### 3.2 ⛔ It relocates the source, which breaks the sidecar path

`ARCHITECTURE.rst`: *"We start by copying `source_paths` to `mutants/` and then mutate the `*.py` files in
there."* The docs add that `also_copy` exists because *"to run the full test suite some files are often needed
above the tests and the source"*.

Consequence: inside a mutant run `__file__` is `mutants/…/ledger_stage023_*.py`, so
`emit_dimension_sidecar(__file__, …)` writes `mutants/…/ledger_stage023_*.dimensions.txt`, while
`compare_dimension_artifacts.py`'s `one_matching_file(SCRIPT_DIR, "ledger_stage023_*_sympy_audit.dimensions.txt")`
globs whichever directory *it* was copied to. The frozen `mathematica/out/*.out` and
`ledger_dimensions.accepted.sha256` would each need `also_copy` entries. Recoverable in principle, but every
path assumption in the two checkers has to be re-established inside a tool-managed tree.

### 3.3 ⛔⛔ Mutation schemata destroy the sha256 pins — the decisive blocker

`ARCHITECTURE.rst`: *"The mutated files contains the original code and the mutants. With the
`MUTANT_UNDER_TEST` environment variable, we can specify (among other things) which mutant should be enabled.
If a mutant is not enabled, it will run the original code."*

So **the bytes on disk are identical for every mutant**; only an environment variable differs. Therefore:

- `emit_dimension_sidecar`'s `source_sha256 = sha256(source.read_bytes())` is a **constant across all mutants**,
  and equals neither the committed source's digest nor anything the frozen artifacts know. The freshness pin
  (`require_fresh_stage_source`) stops discriminating entirely — it degenerates into "the schemata file equals
  itself".
- If `ledger_dimensions.py` falls inside `source_paths`, its schemata-rewritten bytes will never equal the
  digest in `ledger_dimensions.accepted.sha256`, so `check_ledger_dimensions_pin` raises `ModulePinError` and
  `compare()` returns 1 with `COMPARISON_SKIPPED: ledger-dimensions module pin failed` for **every** mutant —
  a uniform false kill that tells you nothing.

This is a head-on collision with the project's architecture: mutmut's design **assumes the source bytes are not
themselves an observable**, and here they are the control.

### 3.4 Targeting and the other failure modes

- **Targeting:** `only_mutate` / `do_not_mutate` are **glob patterns over files** (`"src/api/*"`), not lines.
  Finer control is `do_not_mutate_patterns` (regex over expressions), `# pragma: no mutate` /
  `# pragma: no mutate block` / `# pragma: no mutate start/end` comments, `max_stack_depth`, and
  `mutmut run "my_module.my_function*"` (Unix-glob over *mutant names*, i.e. function scope). The
  function-glob is the closest thing to named targeting, but `SOURCED_DIMS` is a **module-level dict literal**,
  not inside a function, so that handle does not reach it. Everything else is a denylist.
- **Snapshot/restore:** it never touches your working tree — mutation happens in `mutants/`. Good hygiene, bad
  fit (§3.2/§3.3). `mutmut apply <mutant>` writes to the real file and the docs warn *"You should **REALLY**
  have the file you mutate under source code control and committed before you apply a mutant!"*
- **Per-mutant capture:** killed/survived plus 3.5.0's mutation statistics. No user-defined columns.
- **Resume:** ✅ *"You can stop the mutation run at any time and mutmut will restart where you left off."*
- **⭐ Side-effect files:** not addressed anywhere in the docs. The `mutants/` tree is per-run, not per-mutant,
  so an artifact written by mutant *n* is visible to mutant *n+1* — the same contamination, now inside a
  directory you do not control the lifecycle of.

### 3.5 Verdict — **POOR FIT**

Three independent blockers, any one of which is sufficient: pytest-only with no runner hook (§3.1); source
relocated into `mutants/` breaking `__file__`-derived paths (§3.2); and schemata making the on-disk bytes
identical for every mutant, which silently disables both sha256 controls (§3.3). §3.3 is the one that matters
most, because it does not error out — it makes a control **pass vacuously**.

---

## 4. MutPy

Latest **0.6.1, released 2019-11-17**. No release in ~7 years; effectively unmaintained.

- **Test-runner fit — ⛔ NO.** CLI is `mut.py --target <module> --unit-test <tests>`. `--runner` accepts
  **`unittest` (default) or `pytest` (experimental)** and nothing else. There is no arbitrary-shell-command
  option documented.
- **Mutation mechanism — ⛔ fatal here.** *"It applies mutation on AST level"* — in memory, injecting the
  mutated module through the import machinery. The audit script is executed as a **separate process**
  (`python3 <script>`); a subprocess would import nothing from MutPy's patched module registry and would run
  the **unmutated** file. Even if you wrapped the script in a `unittest` case, module-scope
  `raise SystemExit(1)` at import time and the `__main__` guard would have to be restructured.
- **Targeting — ⛔ NO.** Only `--percentage` (random mutation *sampling*) and `--coverage` (mutate covered code
  only). No line, name, or object selection.
- **Snapshot/restore:** N/A — nothing is written to disk, which also means nothing about the `.dimensions.txt`
  side effect is handled.
- **Per-mutant capture:** `--report` (YAML), `--report-html`, `--show-mutants`. Fixed schema.
- **Resume:** none documented.

**Verdict — POOR FIT**, and additionally an unmaintained dependency. Not a candidate.

---

## 5. mutatest (also surveyed)

Latest documented 3.1.0.

- **Test-runner fit — ✅ on paper.** `--testcmds/-t` takes an **arbitrary shell command** (docs example:
  `"pytest -m 'not slow'"`).
- **Mutation mechanism — ⛔ fatal here, and this is the hardest evidence in the survey.** The docs state
  *"Mutatest will manipulate local `__pycache__` files"* — it mutates **compiled bytecode caches**, not source.
  CPython **never writes or reads a bytecode cache for the module executed as `__main__`**.
  **Verified locally on this box:** running `python3 pcachetest/s.py` produced **no** `__pycache__` directory;
  subsequently `import s` produced one. Since the audit script is run *as a script*, a `__pycache__` mutation
  would be a **complete no-op** — mutatest would report every mutant as "survived" while nothing was ever
  actually mutated. That is the worst possible failure mode for this project: a silent, uniform false negative.
- **Targeting — ⛔ NO, and actively wrong-shaped.** `--nlocations/-n` selects **`n` random code locations**
  (default 10) and `--rseed` merely makes the randomness repeatable. `--only/-y` / `--skip/-k` select mutation
  *categories* (`"aa"`, `"bn"`, `"ix"`, …); `--exclude/-e` excludes whole `.py` files. There is no
  named-target mechanism. Random site selection is directly contrary to
  `development_pipeline.md:330-332`.
- **Resume:** none documented.

**Verdict — POOR FIT.**

---

## 6. universalmutator (also surveyed) — the closest thing in shape

Language-agnostic, regexp-rule-driven mutant generator (agroce/universalmutator).

- **Mutants are separate FILES on disk**: `mutate <file> --mutantDir mutants`.
- **Test command is genuinely arbitrary shell**, including multi-stage build/run/diff pipelines. The README's
  own example is
  `analyze_mutants src/main.rs "cargo clean; cargo build; rm mandel.png; target/debug/mandelbrot …; diff mandel.png origmandel.png" --mutantDir mutants`
  — note it already contains a `rm` of a **generated artifact before each trial**, i.e. the tool's usage model
  contemplates exactly the side-effect reset of §0.3, by putting it in the command.
- **Because mutation rules are regexps in a rule file**, an identifier-for-identifier rebinding
  (`SOURCED_DIMS[Z1ret]` → `SOURCED_DIMS[a]`) *is* expressible — the thing cosmic-ray's operator set cannot do
  without a plugin.
- **Results**: `killed.txt` / `not-killed.txt` — two flat filename lists. No per-mutant columns, no
  `stage_exit` / `sidecar_changed` / `records_moved` / `mismatch_names`.
- **Line restriction:** no `--lines`-style flag found in the README. ⚠ **Unverified** — I could not read the
  full `mutate --help` output without installing.
- **Resume:** none documented.

**Verdict — PARTIAL FIT (weakest evidence base of the five).** Worth a bounded local trial before cosmic-ray,
because the two hardest requirements (arbitrary shell oracle; a rebinding mutation) are both native. But it
gives you nothing on resume or on structured per-mutant results, and it is a much smaller, less-maintained
project than cosmic-ray.

---

## 7. Summary table

| requirement | cosmic-ray | mutmut | MutPy | mutatest | universalmutator |
|---|---|---|---|---|---|
| arbitrary shell command as the oracle | ✅ (`shlex.split`, so wrap in `bash -c`) | ⛔ pytest only since 3.0 | ⛔ unittest/pytest only | ✅ `--testcmds` | ✅ native |
| mutates the **real file at the real path** (keeps `__file__`, sidecar path, `sha256` pin honest) | ✅ | ⛔ copies to `mutants/` + schemata | ⛔ in-memory AST | ⛔ `__pycache__` only — **no-op for `__main__`** | ✅ separate mutant files |
| human-chosen include-list of named targets | ⛔ subtractive filters only; custom filter writable | ⛔ file globs / pragmas / function globs | ⛔ random `--percentage` | ⛔ random `--nlocations` | ⚠ unverified |
| identifier-rebinding mutation (`bind.tsv`) | ⛔ no operator; `VariableReplacer` substitutes `randint(-100,100)` | ⛔ | ⛔ | ⛔ | ✅ regexp rules |
| deterministic / re-runnable | ✅ `NumberReplacer` `OFFSETS=[+1,-1]`; ⛔ `VariableReplacer` random | ✅ | ✅ | ⚠ random sites, `--rseed` only | ✅ |
| restores the source between mutants | ✅ `use_mutation` | n/a (isolated tree) | n/a | n/a | ✅ |
| ⭐ resets **emitted side-effect files** | ⛔ no hook — wrapper must | ⛔ | n/a | n/a | ⚠ by convention, in the command |
| per-mutant multi-column result capture | ⚠ `WorkResult.output` in `session.sqlite`; bespoke reader needed | ⛔ killed/survived | ⛔ fixed YAML/HTML | ⛔ | ⛔ two flat lists |
| resume after interrupt | ✅ | ✅ | ⛔ | ⛔ | ⛔ |
| maintained | ✅ 8.4.6, 2026-04 | ✅ 3.6.0, 2026-06 | ⛔ 2019 | ⚠ | ⚠ |
| **verdict** | **PARTIAL** | **POOR** | **POOR** | **POOR** | **PARTIAL (weak evidence)** |

---

## 8. Migration cost, and lock-in to the current approach

### 8.1 Lock-in: essentially none in code, real in process

There **is no harness code to be locked into**. `_scratch/stage023_orch/` contains no driver script — the
ablations were performed by hand (as `docs/development_pipeline.md:19` explicitly authorises: *"A temporary
mutation whose only purpose is to make a gate fire, then be reverted, is part of the review process, and Claude
may perform it"*). The banked artefacts are `GREEN.py`, `GREEN.sidecar`, 22 + N `stdout_*.txt`, 22 `cmp_*.txt`,
and the two TSVs. Nothing imports them; nothing would break.

What **is** binding is process, and it is the part no tool satisfies:
the orchestrator owns the target list, per-tooth over every able-to-fail check and every emitted record
(`:330-332`, `:77`, `:19`). Any tool whose `init` step enumerates the sites has moved that ownership to the
tool.

### 8.2 Concrete cost of adopting cosmic-ray (the only serious candidate)

1. `pip install cosmic-ray` (reachable; 8.4.6 supports 3.10). Decide whether it belongs in the project's
   environment at all.
2. Write `_scratch/ablation_oracle.sh`: reset the sidecar (`cp GREEN.sidecar …dimensions.txt` or `rm -f`), run
   the audit script capturing stdout, diff the `PASS` multiset against `pass_base.txt`, run
   `compare_dimension_artifacts.py 023`, emit the five TSV columns on stdout, exit 0 only if everything holds.
   *This script is the entire oracle and cosmic-ray contributes nothing to it.*
3. Write a **custom include-filter** over `session.sqlite` driven by a human-written target list — because
   `init` will emit order-10³ work items against the wanted ~51 (§2.2).
4. Write and **package** a **custom operator** for the `bind.tsv` identifier rebinding, registered through
   `[project.entry-points."cosmic_ray.operator_providers"]` (§2.3). This is a distributable plugin, not a
   script in `scripts/`.
5. Write a **bespoke `session.sqlite` reader** to reconstruct `results.tsv`/`bind.tsv` from
   `WorkResult.output`, since `cr-report` reports survival rates (§2.4).
6. Re-establish the pipeline's `cp` + `git hash-object` restore proof on top of cosmic-ray's own restore, since
   the process rule names the mechanism, not just the outcome.

Items 2–5 are more code than the thing being replaced.

### 8.3 The honest alternative

A project-owned driver of roughly 100 lines does all of it and satisfies the process clause **by construction**,
because its input *is* the orchestrator's target list:

> read a hand-written table of `(name, line, old_text, new_text)`; for each row —
> `cp GREEN.py <source>`; apply the one textual edit; `cp GREEN.sidecar <sidecar>` (or delete it);
> run the audit script, record `stage_exit` and the `PASS` multiset and whether the sidecar changed and which
> records moved; run `compare_dimension_artifacts.py 023`, record `cmp_exit` and the `MISMATCH` names;
> `cp GREEN.py <source>` and prove the restore with `git hash-object`; append the TSV row.

That is a faithful re-derivation of both existing TSVs, it is resumable by skipping rows already present, and
it makes the ⭐ side-effect reset a first-class step rather than something smuggled into a wrapper.

### 8.4 The one genuinely transferable idea

Even on a NO verdict, cosmic-ray's **durable-session** design is worth copying: a persistent per-mutant record
(the two TSVs already are one) that lets an interrupted run resume by skipping completed rows, and that stores
each mutant's *captured stdout* alongside its verdict. The existing `stdout_*.txt` / `cmp_*.txt` pairs are that
idea implemented by hand; formalising the resume-by-skip is nearly free.

---

## 9. What I could not verify

- **No tool was executed.** Nothing was installed, per the survey's constraint. Every behavioural claim in
  §2–§6 rests on published documentation plus reading each project's source on GitHub. The two exceptions —
  module availability and the `__pycache__`/`__main__` behaviour — are marked as locally verified.
- **cosmic-ray's complete operator inventory.** I read `provider.py`'s registry, but the three
  `*_operator_replacement.operators()` factories expand to further classes I did not enumerate. The claim
  *"no shipped operator replaces one identifier with another identifier"* rests on `provider.py`'s explicit
  list plus `variable_replacer.py`'s source; a `cosmic-ray operators` run would settle it definitively.
- **The custom-filter API.** The docs say only *"it should be straightforward to write your own"*; there is no
  documented `FilterApp` contract. The route via the plain SQLite session and `MutationSpec.start_pos` is
  inferred from `work_item.py`, not from a documented recipe.
- **mutmut issue #373 has no visible maintainer reply.** The "no `--runner` in mutmut 3" statement is the
  reporter's; corroboration is `HISTORY.rst`'s "Pytest only" and the absence of any runner key from the config
  surface. I did not find an explicit maintainer denial or a documented workaround.
- **mutmut against a `SystemExit`-at-module-scope script.** No maintainer statement exists; §3.2/§3.3 are
  derived from `ARCHITECTURE.rst`'s description of the copy + schemata model, not from a run.
- **mutatest's source** was not read; the `__pycache__` mechanism is quoted from its own docs. The
  `__main__`-is-never-cached half of that argument *was* verified locally and is the load-bearing half.
- **universalmutator** is the weakest section: README-level only. Whether it can restrict to given line numbers,
  and exactly how it swaps a mutant file in and restores, are both **open**.
- **The precise semantics of the `sidecar_changed` / `records_moved` columns** in `results.tsv` were inferred
  from the artifacts (`cmp_*.txt` file sizes partition cleanly into 673-byte freshness failures and
  315–320-byte real mismatches); the harness that produced them was manual and left no driver script to read.
