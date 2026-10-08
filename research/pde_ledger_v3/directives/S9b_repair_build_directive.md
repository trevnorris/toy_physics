# S9b repair build directive: Parts A–C in two engines (orchestrator decision list)

**Author:** Claude (orchestrator), 2026-10-08. **Status:** folded once after one Codex + Grok pass (`CLAUDE.md` G2;
both NEEDS REVISION, five findings, all accepted; dispositions in
`directives/_measurements/S9b_repair_build_directive_review_disposition.md`). Item 10 changes what the engines print,
so it was also reviewed as physics-bearing content until clear (CLAUDE.md scope precedence): three rounds, cleared by
both legs (`directives/_measurements/S9b_repair_build_directive_item10_r{1,2}_review_disposition.md`). There are two
builds, one per engine. Each builder reads Part 1 and its own part.

**Scope.** This implements D6 and D7 item 3 of `directives/S9b_repair_decision_list.md` (`79055918`) for Parts A–C
only. **Part D is held** (user decision, 2026-10-08) until the brane-material premise behind P1 is settled. This build
constructs and prints no Part D object. That narrows D7 item 3, which listed Part D. Not in this build: the
comparator, the production runs and the record.

## Part 1. Shared decisions (both engines)

1. **Physics authority.** `CLAUDE.md` and `AGENTS.md` bind. `research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md`
   (Parts A–C accepted at `97580010`) is the sole physics authority and the sole physics input, and it wins every conflict
   with this directive. Construct and print:
   - its "Governing object", order counting, observables, branch-existence conditions and supplied references;
   - its Deliverables Parts A, B and C.

   Its "Part D" deliverable and its "Part D inputs" section are not built here. Add no expected value and no
   acceptance criterion.
2. **The three clauses, verbatim** (`.claude/skills/build/SKILL.md`):
   > **1. The script may PRINT computed objects. It may NOT state conclusions.** **2. PRINT the residual; do
   > NOT assert it.** Compute → emit → *then* assert. **3. Interpretation belongs to the STEP RECORD.**

   The structural rule, verbatim:
   > **The ONLY place the physical symbols may be combined by hand is in CONSTRUCTING THE ACTION and the
   > ANSATZ. Every other expression involving them must be REACHED BY COMPUTATION. Every control re-enters
   > the chain at the ACTION, ⛔ never at a result.**

   In this step, these play the role of the action: the supplied dispersion relation, the supplied mass balance,
   the supplied Part C responses and the supplied references. The radial profiles are the ansatz.
3. **Names and bindings (D6 B1).**
   - A name that matches an upstream `LEDGER` row, or a symbol the spec uses, must denote the same physical
     object. Otherwise rename it. A bare-symbol match is not evidence that two objects are the same.
   - Bind a fold row only where the spec identifies it with a supplied object. Obtain the binding by computation
     and report it with its spec line. No fold row is bound as the value of a live profile.
   - No export changes an upstream row's status or role.
   - Symbols that originate in this step: `n`, `K`, `m`, `GM`, `s`, `f`, `j_n`, `ρ₀`, and item 10's flux `Φ`. None
     has a numeric value.
     `ρ₀` is not `LEDGER['rho_m']`, and `c₀` is not `c_s0`.
4. **Tags (D6 B7).** Use the grammar `<ENGINE>_S9B_<QUANTITY>`, where `<ENGINE>` is `PY` or `WL`. Both engines emit
   the same set of `<QUANTITY>` names, with one tag per named object and one line per tag: `TAG: <payload>`,
   re-parseable. A name names the object, never its value, sign or shape. Engine-local tags carry `_LOCAL_` after
   the engine prefix, and each engine emits one tag listing its `_LOCAL_` names. This follows
   `S11b_SHARED_PHYSICS.md` §10, with the step tag `S9B`.
5. **No verdicts (D6 B5).** No `VERDICT`, `PASS` or `FAIL`. A boolean-valued test is emitted as the CAS object the
   test returned. Emission never depends on a payload's value. Outside the branch-existence conditions the spec
   requires `NOT_ESTABLISHED`. That output, and the branch type, are produced from the computed conditions, never
   typed.
6. **The every-`b` conditions are reduced (D6 B2).**
   - Each Part B and Part C condition is solved relative to `GM` for every far-zone `b`. The solving is case by
     case, over every branch the reduction produces, and each case carries its domain.
   - An unevaluated quantified set, `ConditionSet` or `ForAll` is a restatement, not a reduction.
   - For each reduced condition, print the `j_n` it implies through `∇·(ρ_br V) = −j_n`, with `ρ_br` symbolic and
     the sign of `V` symbolic.
   - The method is the builder's choice (item 13). If an engine cannot reduce a condition, item 12 applies.
7. **The radar coefficient (D6 B3).** The coefficient of `ln(1/b²)` is obtained by computation from Part A's
   computed round-trip excess time, in the regime `Z_E, Z_R ≫ b`. It is not a rule applied by hand. Part B's
   radar `γ` uses that coefficient.
8. **Path dependence (D6 B4).** Whether the nonreciprocal part depends on the path or only on the endpoints is
   computed on the nonreciprocal object Part A computed, not on a placeholder.
9. **Part C domains (D6 B6).** Each Part C row's domain predicates are built from that row's own response. They are
   not built from another row, or from the amplitude before the response is applied.
10. **Forward case and single-mechanism restrictions (user requests, 2026-10-08).**
    - **Forward case: no far-zone loss.** This is the user's premise (2026-10-06/07/08): brane material leaves the
      brane only in throats, so `j_n ≡ 0` in the far zone, and the same brane mass crosses every sphere around the
      body.
      - `Φ` is a new live symbol, defined by one equation: `Φ ≡ ∮ ρ_br V·n_out dA`, the brane mass flux through a
        far-zone coordinate sphere with outward unit normal `n_out`. It has no value and no sign, and its relation
        to `GM` is open.
      - Solve the supplied mass balance for `V` under this premise, with `ρ_br` live, and express the integration
        constant through `Φ`.
      - Substitute that `V` into Part A's observables, first with `δ ≡ 0`, `ξ_w ≡ 0`, then with `δ` and `ξ_w`
        live.
      - First print the branch-existence and path-traversal conditions after the same substitution. Gate the
        forward observables on them as Part A does.
      - Print the deflection, the round trip's `ln(1/b²)` coefficient (item 7), both effective `γ`s, their
        difference and both Part B residuals against the references, each as a function of `b`. Also print the round-trip excess time, both
        one-way excess times and their nonreciprocal part, with `b`, `Z_E` and `Z_R` kept.
      - Print no matching condition for this case. Label every object with the premise.
    - For the deflection and for the radar `ln(1/b²)` coefficient, print each restriction below of:
      - Part A's computed observable;
      - Part B's effective `γ` and residual;
      - the difference between the two effective `γ`s;
      - Part B's condition, reduced as in item 6, with its implied `j_n`.
    - The restrictions:
      - **flow only:** `δ ≡ 0`, `ξ_w ≡ 0`, `V` live;
      - **speed only:** `V ≡ 0`, `ξ_w ≡ 0`, `δ` live;
      - **tilt only:** `δ ≡ 0`, `V ≡ 0`, `ξ_w` live.
    - For Part C's fixed-ratio and power responses, also print each condition restricted to `V ≡ 0`, `ξ_w ≡ 0`
      (bulk density alone), with its domain and the `j_n` it implies through the supplied mass balance, with
      `ρ_br` symbolic.
    - Each restriction is obtained by substituting into the computed general object, then reducing. Its tag and
      label name the restriction and the profiles it sets to zero, never an outcome. Restrictions add no premise.
    - Everything in this item is additional labelled output. This directive supplies the forward case's premise for
      that label only. Building this item is not a stop event under item 12, and it does not replace the live `j_n`
      in Parts A–C.
11. **Running.**
    - Run every CAS demonstration through `scripts/s11c_guarded_run.py` with `--pool s9b --memory-gib 8`. Wolfram
      runs also take `--tasks-max 64`. Both limits rest on the user's authorizations of 2026-10-06 and 2026-10-07.
    - If the guard refuses admission, kills a run, or a kernel hits the task limit, stop and report it. A higher
      limit needs the user. ⛔ Never answer a kill by narrowing or cheapening the requested object.
    - There are no time limits (`AGENTS.md`).
    - Run Wolfram scripts as `math -script <file>`.
    - Demonstration output goes under your own scratch paths, never under `scripts/out/` or `mathematica/out/`. The
      production runs and their `out/*.out` files belong to the orchestrator after review.
12. **Stop and report**, without choosing, when any of these happens (item 10's labelled outputs are not such events):
    - a second method failure;
    - a sub-problem the spec does not name;
    - a premise the spec does not supply.
13. **Deferred to the build** (implementation, not new physics). The builder chooses, implements and reports each
    of these. No choice may turn a live profile into a chosen one or add an expected value.
    - The symbolic handling of the every-`b` requirement, within item 6.
    - Component calculus on the supplied coordinate mass law.
    - The representation of general live radial profiles. A family restriction follows the spec's "Profiles" item.
14. **Ablation harness** (`docs/development_pipeline.md` §4).
    - Each engine gets a committed harness that runs the live engine unchanged as the baseline. For each knife
      below, it runs a copy with exactly that one mutation at the named construction site.
    - It prints the baseline payload, the corrupted payload and their difference for **every** tag, including tags
      whose difference is zero. It prints no verdict.
    - A knife whose construction an engine does not contain is reported, not invented. The engine itself contains
      no control with an expected outcome.

    The knives:
    - **K1, advection:** in the supplied kinetic symbol, `ω − V^i k_i` is replaced by `ω`.
    - **K2, metric:** `g^{ij}` is replaced by `δ^{ij}`.
    - **K3, local speed:** `c_γ(x)` in the dispersion relation is replaced by `c₀`.
    - **K4a–c, profile freezes:** each of the Parts A–C radial profiles `δ`, `V` and `ξ_w` is replaced by a constant
      where it is differentiated in the ray construction (three separate knives).
    - **K5, return leg:** the round trip's return leg propagates in the outgoing leg's direction.
    - **K6, radar source (a named coefficient knife):** the object from which the `ln(1/b²)` coefficient is
      extracted is replaced by one one-way excess time. It is named for a channel no FORM knife sees: which computed
      object the coefficient is taken from (`docs/development_pipeline.md` §4). K1–K3 are the FORM knives on the
      round-trip construction that feeds the coefficient.
    - **K7, quantifier:** the every-`b` reduction is carried out at one fixed far-zone `b`.
    - **K8, mass balance:** `ρ_br` in the supplied mass balance is replaced by a constant.
    - **K9a–b, Part C responses:** in the fixed-ratio response, the local `c_s(x)` is replaced by `c_s0`; in the
      power response, `ρ(x)` is replaced by `ρ₀` (two separate knives).
    - **K10, branch gate:** the computed branch-existence condition is replaced by its logical negation before the
      observables are gated on it.
    - **K11, non-radial flow:** at the advection construction, the advected velocity gains an azimuthal component,
      so it is no longer radial. It acts on the construction from which item 8 computes path dependence.
15. **Handoff.** Each builder works in its own fresh repository, exported from the commit that holds this
    directive, with no git history. The other engine's S9b files are absent from it. Builders are Codex
    `gpt-6-astra` at high effort, a fresh session per engine.
16. **Stop.** Build, run demonstrations, report, then stop. Call, spawn or launch no other AI, reviewer or agent.
    Commit nothing; push nothing.
17. **Report** (in the final message):
    - the tags emitted;
    - each deferred item (item 13) and how it was handled;
    - each binding and its spec line;
    - each family restriction and what it leaves out;
    - each knife, and where it acts;
    - every `NOT_ESTABLISHED`;
    - every guard refusal or kill;
    - every stop-and-report event (item 12).

## Part 2. SymPy engine

- **Repair and extend** `research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py`. **Write** its harness
  `research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py`. The engine's only writes are its stdout tag
  stream and the delta `research/pde_ledger_v3/scripts/S9b_exports.py`.
- **Import** through the fold: `load_model('scripts/S11c_b_exports.py', 'scripts/S11c_c1_exports.py',
  'scripts/S11c_c2_exports.py')` from `scripts/ledger_fold.py`. ⛔ Fold neither `scripts/S9b_exports.py` (D6
  Exclusion) nor `scripts/O2_exports.py`.
- **Export**, following `directives/export_ledger_bind_closure_design.md` D1–D3, as the c1 and c2 engines do:
  - The delta holds only this step's own rows.
  - The outgoing bind-set is the reduced Part B and Part C conditions and the `j_n` each implies.
  - `IMPORT_KEYS` is the exact manifest of the fold rows this engine binds.
  - On an `F9` collision, `F9b` writes the bare key only when the two rows are proved to be the same physical
    object (item 3). Otherwise `F9c` writes `s9b_<key>` and leaves the imported row as it stands. The report gives
    the comparison's three-valued outcome.
  - Keep `F3`, the `D3` round-trip, the `_RELATIONALS` reviver and the freeze.
  - `BUILD_INPUT_DIGESTS` pins every executable input consumed: this engine's source, the spec, the three fold
    files and `ledger_fold.py`.
  - Publish the delta only if Parts A, B and C, branch existence and item 10 all completed (`F6`, first branch).
- **Mechanical precedent, not authority:** `scripts/S11c_c2_selfenergy_fold_sympy_audit.py` for the fold, the
  manifest and the delta; `scripts/S11b_interface_coupling_law_sympy_audit.py` for the emission shape.

## Part 3. Wolfram engine (blind)

- **Repair and extend** `research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl`. **Write** its
  harness `research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl`. The engine's only write is
  its stdout tag stream.
- **This is the blind engine.** It reads no file at run time, imports nothing, and re-derives everything from the
  spec's equations. It imports no `LEDGER` and writes none.
- **Mechanical precedent, not authority:** `mathematica/S11b_interface_coupling_law_mathematica_audit.wl`, for its
  emit, naming and flush shape only.
- **At most one kernel at a time.**
- **Executable checks, with no expected values:**
  1. Copy the finished `.wl` alone into an empty scratch directory. Run it there and in the repository. Both runs
     exit 0, and their streams are non-empty and byte-identical. Afterwards the scratch directory holds only the
     copied file.
  2. With stdout redirected, the capture file is observed growing while the kernel is alive. If every task finishes
     too fast to observe this, report that.
  3. `git status --porcelain research/pde_ledger_v3/mathematica/out/` prints nothing.
