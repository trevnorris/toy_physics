# Install and reproduce the Lean checks

These instructions cover the completed S9, S10 and bounded S11 Lean contracts.
Use the checkout containing the proofs and their `_measurements/*_contract_checks.json`
records. The independently reviewed D5 density classification is included in `all`;
the completed D5 bulk contract has a separate recorded suite and remains
outside this portable catalog.

The supported command-line path is Linux, or Linux inside WSL2 on Windows.
The tools also have macOS releases, but this repository's setup/replay path has
not been validated on macOS or native Windows. An editor is optional; VS Code
with the Lean extension can use the same environment.

## 1. System prerequisites

On Ubuntu/Debian or Ubuntu under WSL2:

```sh
sudo apt-get update
sudo apt-get install git curl ca-certificates zstd build-essential python3 python3-venv python3-pip
```

Use Python **3.10 or newer**. The Python environment is pinned to **SymPy 1.14.0**
and **mpmath 1.3.0** in `requirements.txt`; the replay tooling itself otherwise
uses the standard library. Git, curl, tar, zstd and a C compiler must be on PATH.
Install equivalent packages with your distribution's package manager if needed.

You need access to this repository. From a checkout:

```sh
cd research/pde_ledger_v3/lean
bash setup.sh --check
bash setup.sh
export PATH="${ELAN_HOME:-$HOME/.elan}/bin:$PATH"
.venv/bin/python verify.py --doctor
```

To select a different installed Python, run `PYTHON=python3.12 bash setup.sh`.
`--check` does not download or install anything. Plain `setup.sh`:

1. Installs [elan](https://github.com/leanprover/elan) when it is missing, using
   its official HTTPS installer, without modifying shell startup files or
   selecting a global default toolchain.
2. Installs the exact Lean release in `lean-toolchain`, including Lake.
3. Creates `.venv/` and installs `requirements.txt` into it.
4. Retrieves cached dependencies using the committed `lake-manifest.json`,
   then builds the three Physlib imports used by this ledger.
5. Checks the environment, pinned package revisions, direct imports and
   recorded control sources with `verify.py --doctor`.

Setup requires network access to Lean/elan releases, dependency repositories,
the Mathlib cache and PyPI. It writes `setup.log`, the project `.venv/` and
`.lake/` caches, plus elan's toolchain installation. Dependency caches are large;
allow several GB of disk space and extra space for fresh replay objects/logs.
Initial downloads can take considerably longer than subsequent runs.

The committed pins are Lean **4.33.0**, Physlib
`8b2b23701c07409f882c49e0bd290a214e44450c` and Mathlib
`db584cd6d46c92f209a44c0f1c829460d327499d`. Do not use `lake update` to repair an
installation: that can change the dependency selection. Restore missing pin
files from the same checkout instead. The doctor rejects dirty tracked package
sources or revisions different from the manifest.

## 2. Run proofs and controls

For an initial small run, use the completed analytic-error contract:

```sh
.venv/bin/python verify.py analytic-error
```

For all completed contracts, including S11:

```sh
.venv/bin/python verify.py all
```

Select one or several contracts, inspect the plan, or run only one layer:

```sh
.venv/bin/python verify.py --list
.venv/bin/python verify.py all --plan
.venv/bin/python verify.py s9 d2
.venv/bin/python verify.py d4-bulk --build-only
.venv/bin/python verify.py all --native-only
```

The available names are `s9`, `s10`, `homogeneous`, `d2`, `d2-dynamics`, `d3`,
`d3-bulk`, `d4`, `d4-odd`, `analytic-error`, `variable`, `poles`, `d4-bulk`, and `d5`.
`s10` includes the retained, already committed CAS bridge proofs; this runner
does not generate or expand that bridge. The complete run is substantial.
`--plan` validates control inputs and shows the work without launching Lean.

Each run prints its report location, then logs silently under
`research/pde_ledger_v3/_scratch/lean_portable/<run>/`. The final `PASS` message
and exit status zero indicate success for the selected mode. Inspect
`report.json`, `logs/`, and `snapshot/` on failure. A partial run remains failed;
it is never presented as a successful mutation test. Each invocation starts
with fresh local objects; there is currently no portable `--reuse-build` option.

Local proof builds and controls are sequential, with one Lean thread,
`-M4096`, strict warnings, and a default **600-second limit per process**.
The allocator limit is not an operating-system RSS cap. Leave additional RAM
for Python, the OS and any concurrent calculation jobs. On a slower machine,
an explicit `--timeout 1200` increases the process limit; a timeout still never
counts as mathematical rejection. Dependency installation uses Lake/cache
tools and is separate from these one-worker proof checks.

## 3. What is reproduced

The portable runner:

- Copies the selected current sources and test inputs into a separate snapshot.
- Rebuilds every local transitive Lean import in dependency order, using only
  the installed external package caches. It cannot read old local ledger
  `.olean` files through its `LEAN_PATH`.
- Checks printed axiom lists, including empty and multiline lists, against
  `propext`, `Classical.choice`, and `Quot.sound` only.
- Replays the exact control sources recorded in the completed contract reports.
  Older source mutations are reconstructed only when their replacement applies
  exactly once and the resulting source hash matches. Passing controls must
  compile; failed controls must reach the intended mathematical diagnostic.
- Checks D3/D4/D5 generated certificates with their existing `--check` mode and
  runs the selected compact SymPy source checks in the snapshot. Original
  reports are never overwritten. For the variable-coefficient check, the two
  upstream compact source checks run first in that same snapshot.
- Records current input, package, direct-import and output-object hashes. It
  detects input drift during the run. Use a stable checkout when reproducing.

Recent paired controls require exactly one error showing `False` in the named
control theorem. Some older source mutations also break downstream tactics or
produce a linter diagnostic; those secondary failures are not their evidence.
The runner retains the reviewed named mathematical failure. The S9
`phase_wrong_coordinate` case has an explicit exception limited to its exact
reviewed source and displayed false coordinate identity, not arbitrary rewrite
failures. Unit regressions cover these distinctions and resource-failure handling.

The reports already in `_measurements/` remain historical review evidence.
Their original author-oriented scripts can require historical object hashes,
absolute paths or preservation manifests from that workspace. **Use `verify.py`
for fresh-checkout replay**, rather than editing those guards or replacing their
records. New portable results are new execution evidence, not a claim that
compiler outputs are byte-identical across machines or that independent fidelity
review was repeated. See `FORMALIZATION_POLICY.md` for that distinction.

No Wolfram license, Claude/Grok account, GPU, CUDA installation or production
CAS run is required for these commands. The compact source checks execute
selected existing SymPy helpers; Wolfram anchors are source inspection only.
The pole check also reads the existing synthetic repair evidence as historical
input. Production CAS, comparator/PIT pipelines and independent review have
their own requirements and are not reproduced by this runner. Stored external
package caches remain part of the pinned dependency baseline; this is not a
rebuild of Mathlib from scratch.

## Troubleshooting and maintenance

- **`elan`, `lake` or `lean` not found:** export the elan PATH shown above.
  The setup script does this for itself; a parent shell needs its own export.
- **Missing `.venv/bin/python` or SymPy mismatch:** rerun setup with Python 3.10+
  and use `.venv/bin/python` for the checks.
- **Missing external import/cache or interrupted download:** inspect `setup.log`
  and rerun setup with network access. Do not remove the dependency pins.
- **Missing tracked input/annex content:** restore the named file from the
  checkout. The portable checks need the selected source and recorded compact
  evidence files, not every production export. Do not substitute a newer export
  into a historical packet to make a hash check pass.
- **Source drift or dirty dependency:** finish the edit or switch to a stable
  checkout, then start a new replay. Existing records remain unchanged.
- **Long run:** inspect its local logs when needed. There is no recurring
  model polling or automatic external transfer in this runner.

To test the replay instrument itself without a Lean build:

```sh
.venv/bin/python -m unittest -v test_verify test_setup
```

Ordinary editor-oriented Lake builds are also available after setup:

```sh
lake build S9Pilot
lake build S11D4Bulk
```

Those commands write the shared local `.lake/build/` directory and do not run
the contract mutation/native checks. Bare `lake build` uses the existing default
targets, which cover **S9/S10 only**. Use the explicit S11 target or the portable
`all` command when you want S11 coverage. Do not rebuild shared targets while
another job is guarding their object hashes.

## Validation of this installation path

Validated on Linux with Python 3.10.12 and the pinned dependencies on
2026-09-17. The compact record is [INSTALL_VALIDATION.json](INSTALL_VALIDATION.json).

| Check | Observed result |
|---|---|
| Tooling regressions | Ten passed, including bootstrap control flow, missing prerequisites, mutation diagnostic filtering and timeout child cleanup. |
| Compact native replay | All twelve available native suites passed in isolated snapshots. |
| Analytic-error replay | Five fresh local objects; 41 axiom audits; twelve mathematical rejections and sixteen positive controls. |
| S9 replay | Ten fresh local objects; 32 axiom audits; nine mathematical rejections and five positive controls. |
| Preservation | All 383 historical source/evidence files, 51 protected objects and two native inputs in the D5 preservation manifest remained unchanged. |

The two proof smoke runs used new local object directories and existing pinned
external caches. The bootstrap tests replaced network/install commands with
test stubs; elan, Lean and Python packages were not redownloaded onto a new OS.
The full `verify.py all` proof run, macOS and native Windows were not exercised
as part of this tooling change. All completed contract control inputs were
validated, and their historical diagnostic formats were covered by the tooling
tests. This validation adds no mathematical claim or fidelity-review clearance.

## D5 catalog integration

The completed D5 density classification (`ec27ecf0`) is now selectable with
`.venv/bin/python verify.py d5` and included in `all`. It has 45 local modules,
51 selected axiom audits and 30 paired/positive control executions. The runner
also checks the unchanged D5 generator and compact native span instrument.
The fresh D5 replay passed on Linux on 2026-09-18: all 45 objects, 51 axiom
audits, thirteen intended mathematical rejections, seventeen positive controls,
and the generator/native checks. Ten tooling regressions also pass. The run
used one worker and took approximately 73 minutes on this host; other hardware
may differ. Live inputs, logs, objects, fifteen clean dependency pins and five
direct Mathlib source/object pairs were checked after completion.

[INSTALL_D5_VALIDATION.json](INSTALL_D5_VALIDATION.json) records this new
execution evidence. All 499 protected historical source/evidence files and 232
shared ledger objects remained unchanged. The dated installation validation
above remains historical and is not rewritten. This replay used the existing
Python 3.10.12 environment with the pinned SymPy/mpmath versions and installed
external Lean caches; it did not repeat a fresh OS/bootstrap download test or
run the entire `all` proof catalog. The now-expanded catalog plans 232 unique
local modules and 349 control executions; those plan counts are not a claim of
a new full-catalog execution.

D5 bulk subsequently completed verification and both independent reviews; see
[s11/D5_BULK_FIDELITY_REVIEW.md](s11/D5_BULK_FIDELITY_REVIEW.md). It has not yet
been registered in the portable catalog. This closure updates status prose only;
`INSTALL_D5_VALIDATION.json` retains the hashes of the documents at replay time.
Runner, setup, tests and execution evidence are unchanged.
