# S9b build directive: two engines (orchestrator decision list)

**Author:** Claude (orchestrator), 2026-10-06. **Status:** folded once after one Codex + Grok pass
(`CLAUDE.md` G2). There are two builds, one per engine. Each builder reads Part 1 and its own part.

## Part 1. Shared decisions (both engines)

1. **Physics authority.** `CLAUDE.md` and `AGENTS.md` bind. `research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md`
   is the sole physics authority and the sole physics input, and it wins every conflict with this directive.
   Implement:
   - Parts A, B and C;
   - branch existence;
   - the order counting;
   - the supplied references;
   - the items under "Deferred to the build".

   The spec assigns the deferred items to this directive. This directive delegates them to the builder (user
   decision, 2026-10-06). The builder chooses, implements, and reports each choice (item 9). Add no expected
   value and no acceptance criterion.
2. **The three clauses, verbatim** (`.claude/skills/build/SKILL.md`):
   > **1. The script may PRINT computed objects. It may NOT state conclusions.** **2. PRINT the residual; do
   > NOT assert it.** Compute → emit → *then* assert. **3. Interpretation belongs to the STEP RECORD.**

   The structural rule, verbatim:
   > **The ONLY place the physical symbols may be combined by hand is in CONSTRUCTING THE ACTION and the
   > ANSATZ. Every other expression involving them must be REACHED BY COMPUTATION. Every control re-enters
   > the chain at the ACTION, ⛔ never at a result.**

   In this step, these play the role of the action: the supplied dispersion relation, the supplied mass
   balance, the supplied Part C responses and the supplied references. The profiles are the ansatz.
3. **Symbols that originate in this step:** `n`, `K`, `m`, `GM`, `s`, `f`, `j_n` and `ρ₀`. None has an
   upstream row or a numeric value. `ρ₀` is not `LEDGER['rho_m']`, which is a mass density. `c₀` is the far
   value of `c_γ`, so `c₀² = μ_⊥/ρ_br` far from the mass holds as a relation. `c₀` is not `c_s0`.
4. **Tags.** Use the grammar `<ENGINE>_S9B_<QUANTITY>`, where `<ENGINE>` is `PY` or `WL`. Both engines produce
   parallel tag sets, with one tag per named object, and one line per tag: `TAG: <payload>`, re-parseable.
   A name names the object, never its value, sign or shape. Engine-local tags carry `_LOCAL_` after the engine
   prefix. Each engine emits one tag listing its `_LOCAL_` names. This follows `S11b_SHARED_PHYSICS.md` §10,
   with the step tag `S9B`.
5. **No verdicts.** No `VERDICT`, `PASS` or `FAIL`. A boolean-valued test is emitted as the CAS object the test
   returned. Emission never depends on a payload's value. A refusal, `NOT_ESTABLISHED` with what is missing,
   is a valid output.
6. **No comparator here.** The cross-engine comparator is a separate artifact, built after both engines
   (precedent: `S11b_sympy_build_directive.md`, deviation (a)).
7. **Running.** Run every CAS demonstration through `scripts/s11c_guarded_run.py` with `--memory-gib 4`. There
   are no time limits (`AGENTS.md`). If the guard kills a run, record it and report it. ⛔ Never answer a kill by
   narrowing or cheapening the requested object. Demonstration output goes under your own scratch paths, never
   under `scripts/out/` or `mathematica/out/`. The production runs and their `out/*.out` files belong to the
   orchestrator after review.
8. **Stop.** Build, run demonstrations, report, then stop. Call, spawn or launch no other AI, reviewer or
   agent. Commit nothing; push nothing.
9. **Report** (in the final message):
   - the tags emitted;
   - each "Deferred to the build" item and how it was handled;
   - every `NOT_ESTABLISHED`;
   - every guard kill;
   - every conflict or gap, under the spec's "Stop and report".

## Part 2. SymPy engine

- **Write** `research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py`. Its only writes are its stdout
  tag stream and the delta `research/pde_ledger_v3/scripts/S9b_exports.py`.
- **Import** through the fold: `load_model('scripts/S11c_b_exports.py', 'scripts/S11c_c1_exports.py',
  'scripts/S11c_c2_exports.py')` from `scripts/ledger_fold.py`, as `scripts/S11c_d_mixing_scattering_sympy_audit.py`
  does.
- **Bindings.** These are the far-field values only. ⛔ Every radial profile stays live.
  - Far-field `ρ_br` → `LEDGER['rho_br']`.
  - `c_s0` → `LEDGER['c_s0']`.
  - Far-field `μ_⊥` is the transverse stiffness of `LEDGER['transverse_dispersion']`, the stored form of the
    spec's Anchor branch, obtained from that row by computation. ⛔ It is not `LEDGER['mu_R']` alone, because
    `steps/S11b_interface_coupling_law.md:82–85` records that this engine's basis differs from S9's. `mu_S`
    stays live.
- **Export**, following `directives/export_ledger_bind_closure_design.md` D1–D3, as the c1 and c2 engines do:
  - The delta holds only this step's own rows.
  - The outgoing bind-set is the Part B and Part C conditions and the `j_n` each implies.
  - `IMPORT_KEYS` is the exact manifest of the fold rows this engine binds.
  - On an `F9` collision, `F9b` writes the bare key when equality is proved. Otherwise `F9c` writes `s9b_<key>`
    and leaves the imported row as it stands. The report gives the comparison's three-valued outcome.
  - Keep `F3`, the `D3` round-trip, the `_RELATIONALS` reviver and the freeze, as in
    `S11b_sympy_build_directive.md`.
  - The conditions are relations.
  - `BUILD_INPUT_DIGESTS` pins every executable input consumed: this engine's source, the spec, the three fold
    files and `ledger_fold.py`.
  - Publish the delta only if Parts A, B, C and branch existence all completed (`F6`, first branch).
- **Mechanical precedent, not authority:** `scripts/S11c_c2_selfenergy_fold_sympy_audit.py` for the fold, the
  manifest and the delta; `scripts/S11b_interface_coupling_law_sympy_audit.py` for the emission shape.

## Part 3. Wolfram engine (blind)

- **Write** `research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl`. Its only write is its
  stdout tag stream.
- **This is the blind engine.** It reads no file at run time, imports nothing, and re-derives everything from
  the spec's equations. It imports no `LEDGER` and writes none.
- **Handoff.** The orchestrator runs this build in a fresh worktree of the commit that holds this directive.
  No S9b SymPy engine exists there.
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
