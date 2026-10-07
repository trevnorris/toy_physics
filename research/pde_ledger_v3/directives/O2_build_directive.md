# O2 build directive: two engines and their ablation harnesses (orchestrator decision list)

**Author:** Claude (orchestrator), 2026-10-07. **Status:** folded once after one Codex + Grok pass (`CLAUDE.md` G2;
both NEEDS REVISION: no uniform-stiffness binding, refusals only where the spec allows, and a full named knife list).
There are two builds, one per engine. Each builder reads Part 1 and its own part.

## Part 1. Shared decisions (both engines)

1. **Physics authority.** `CLAUDE.md` and `AGENTS.md` bind. `research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md`
   (accepted `4680e251`) is the sole physics authority and the sole physics input, and it wins every conflict with
   this directive. Construct and print the objects of its §9. Its §10 "Deferred to the build" items are delegated
   to the builder, except the script-control tests, which item 8 decides. The builder chooses, implements and
   reports each deferred choice. No choice may turn an OPEN input into a chosen response. Add no expected value
   and no acceptance criterion.
2. **The three clauses, verbatim** (`.claude/skills/build/SKILL.md`):
   > **1. The script may PRINT computed objects. It may NOT state conclusions.** **2. PRINT the residual; do
   > NOT assert it.** Compute → emit → *then* assert. **3. Interpretation belongs to the STEP RECORD.**

   The structural rule, verbatim:
   > **The ONLY place the physical symbols may be combined by hand is in CONSTRUCTING THE ACTION and the
   > ANSATZ. Every other expression involving them must be REACHED BY COMPUTATION. Every control re-enters
   > the chain at the ACTION, ⛔ never at a result.**

   In this step, the spec's supplied equations and named OPEN operands play the role of the action, and its
   radial profiles (§1) are the ansatz. The component objects are reached by computation from them.
3. **Names.** A name that matches an upstream `LEDGER` row, or a symbol the spec uses, must denote the same
   physical object. Otherwise rename it. `w` is the bulk-normal coordinate. No symbol gets a numeric value.
4. **Tags.** Use the grammar `<ENGINE>_O2_<QUANTITY>`, where `<ENGINE>` is `PY` or `WL`. Both engines produce
   parallel tag sets, with one tag per named object, and one line per tag: `TAG: <payload>`, re-parseable.
   A name names the object, never its value, sign or shape. Engine-local tags carry `_LOCAL_` after the engine
   prefix. Each engine emits one tag listing its `_LOCAL_` names. This follows `S11b_SHARED_PHYSICS.md` §10,
   with the step tag `O2`.
5. **No verdicts.** No `VERDICT`, `PASS` or `FAIL`. A boolean-valued test is emitted as the CAS object the test
   returned. Emission never depends on a payload's value. An OPEN action the spec says to print is printed as
   an OPEN action, never refused or left absent. `NOT_ESTABLISHED`, with what is missing, is valid only where the
   spec says to report a missing restriction.
6. **No comparator here.** The cross-engine comparator is a separate artifact, built after both engines.
7. **Running.** Run every CAS demonstration through `scripts/s11c_guarded_run.py` with `--memory-gib 4` and the
   guard's default task limit. There are no time limits (`AGENTS.md`). If the guard kills a run or a kernel hits
   the task limit, stop and report it; a higher limit needs the user's authorization. ⛔ Never answer a kill by
   narrowing or cheapening the requested object. Run Wolfram scripts as `math -script <file>`. Demonstration
   output goes under your own scratch paths, never under `scripts/out/` or `mathematica/out/`. The production
   runs and their `out/*.out` files belong to the orchestrator after review.
8. **Ablation harness** (`docs/development_pipeline.md` §4). Each engine gets a committed harness that runs the
   live engine unchanged as the baseline. For each knife below, it runs a copy with exactly that one mutation at
   the named construction site. It prints the baseline payload, the corrupted payload and their difference for
   **every** §9 tag, including tags whose difference is zero. It prints no verdict. A knife whose construction an
   engine does not contain is reported, not invented. The knives:
   - **K1, mass law:** the live `ρ_br` in the supplied mass law is replaced by a constant.
   - **K2, graph velocity:** the bulk-direction component of the material velocity is dropped, so the supplied
     graph no longer enters it.
   - **K3, carried velocity:** premise 3's carried velocity `V` is replaced by an independent velocity.
   - **K4, metric:** the induced metric `g_ij` is replaced by `δ_ij`.
   - **K5, reference history:** the dependence of the full live stress on `ℛ_ref/strain^live` is removed.
   - **K6a–c, radial profiles:** in the material momentum storage and transport construction, each of `V_r`,
     `ρ_br` and `ξ_w` is replaced by a constant where it is differentiated (three separate knives).
   - **K7, tilted face normal:** the face normal used for bulk-normal loading is replaced by the untilted normal.
   - **K8, material transport:** the spatial transport part of the material momentum construction is removed,
     leaving storage only.
   - **K9, power pairing:** in the mechanical face/support power pairing, the paired velocity is replaced by an
     independent velocity.
   - **K10, energy transport:** the live energy transport `𝒥_E^live` is removed from the energy accounting.
   - **K11, relaxation power:** `𝒫_ref/relax^live` is removed from the energy accounting.
   - **K12, drive:** a separate body-force entry is inserted into the momentum construction.
   - **K13, carried bulk component:** the OPEN bulk-direction carried action is replaced by `j_n` times the graph
     velocity's bulk-direction component.

   The engine itself contains no control with an expected outcome.
9. **Stop.** Build, run demonstrations, report, then stop. Call, spawn or launch no other AI, reviewer or
   agent. Commit nothing; push nothing.
10. **Report** (in the final message):
    - the tags emitted;
    - each "Deferred to the build" item and how it was handled;
    - each knife, and where it acts;
    - every `NOT_ESTABLISHED`;
    - every guard kill;
    - every conflict or gap the spec says to stop and report.

## Part 2. SymPy engine

- **Write** `research/pde_ledger_v3/scripts/O2_live_balance_sympy_audit.py` and its harness
  `research/pde_ledger_v3/scripts/O2_live_balance_sympy_ablation.py`. The engine's only writes are its stdout
  tag stream and the delta `research/pde_ledger_v3/scripts/O2_exports.py`.
- **Import** through the fold: `load_model('scripts/S11c_b_exports.py', 'scripts/S11c_c1_exports.py',
  'scripts/S11c_c2_exports.py')` from `scripts/ledger_fold.py`. ⛔ Do not fold `scripts/S9b_exports.py`.
- **Bindings.** Bind a fold row only where the spec identifies it with a supplied object, obtained by
  computation. No fold row is bound as the value of a live profile. The uniform transverse anchor is a
  historical-domain record (spec §8.4), not the live `μ_⊥(r)`. ⛔ Every radial profile stays live.
- **Export**, following `directives/export_ledger_bind_closure_design.md` D1–D3, as the c1 and c2 engines do:
  - The delta holds only this step's own rows.
  - The outgoing bind-set is the constructed `ℬ_hold^live` component objects and `ℬ_E^steady`.
  - `IMPORT_KEYS` is the exact manifest of the fold rows this engine binds.
  - On an `F9` collision, `F9b` writes the bare key when equality is proved. Otherwise `F9c` writes `o2_<key>`
    and leaves the imported row as it stands. The report gives the comparison's three-valued outcome.
  - Keep `F3`, the `D3` round-trip, the `_RELATIONALS` reviver and the freeze.
  - `BUILD_INPUT_DIGESTS` pins every executable input consumed: this engine's source, the spec, the three fold
    files and `ledger_fold.py`.
  - Publish the delta only if every §9 object was printed.
- **Mechanical precedent, not authority:** `scripts/S11c_c2_selfenergy_fold_sympy_audit.py` for the fold, the
  manifest and the delta; `scripts/S11b_interface_coupling_law_sympy_audit.py` for the emission shape.

## Part 3. Wolfram engine (blind)

- **Write** `research/pde_ledger_v3/mathematica/O2_live_balance_mathematica_audit.wl` and its harness
  `research/pde_ledger_v3/mathematica/O2_live_balance_mathematica_ablation.wl`. The engine's only write is its
  stdout tag stream.
- **This is the blind engine.** It reads no file at run time, imports nothing, and re-derives everything from
  the spec's equations. It imports no `LEDGER` and writes none.
- **Handoff.** The orchestrator runs this build in a fresh worktree of the commit that holds this directive.
  No O2 SymPy engine exists there.
- **Mechanical precedent, not authority:** `mathematica/S11b_interface_coupling_law_mathematica_audit.wl`, for
  its emit, naming and flush shape only.
- **At most one kernel at a time.**
- **Executable checks, with no expected values:**
  1. Copy the finished `.wl` alone into an empty scratch directory. Run it there and in the repository. Both
     runs exit 0, and their streams are non-empty and byte-identical. Afterwards the scratch directory holds
     only the copied file.
  2. With stdout redirected, the capture file is observed growing while the kernel is alive. If every task
     finishes too fast to observe this, report that.
  3. `git status --porcelain research/pde_ledger_v3/mathematica/out/` prints nothing.
