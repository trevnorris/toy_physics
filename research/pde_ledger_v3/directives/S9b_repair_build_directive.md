# S9b repair build directive: Parts A–C in two engines (orchestrator decision list)

**Author:** Claude (orchestrator), 2026-10-08. **Status:** folded once after one Codex + Grok pass (`CLAUDE.md` G2;
both NEEDS REVISION, five findings, all accepted; dispositions in
`directives/_measurements/S9b_repair_build_directive_review_disposition.md`). Item 10 changes what the engines print,
so it was also reviewed as physics-bearing content until clear (CLAUDE.md scope precedence): three rounds, cleared by
both legs (`directives/_measurements/S9b_repair_build_directive_item10_r{1,2}_review_disposition.md`). There are two
builds, one per engine. Each builder reads Part 1 and its own part.

**Amendment 1** (user, 2026-10-08). The forward case in item 10 now prints Part B's condition and its Part C
rewrites, where it printed none before, and the memory limit in item 11 is raised. Dispositions:
`directives/_measurements/S9b_repair_build_directive_amend1_review_disposition.md`.

**Amendment 2** (user, 2026-10-08). Item 14 evaluates K11's copy compactly. One two-leg pass, folded once (G2).
Dispositions: `directives/_measurements/S9b_repair_build_directive_amend2_review_disposition.md`.

**Amendment 3** (user, 2026-10-09), for the build's repair round 1 after review r0
(`directives/_measurements/S9b_repair_build_r0_review_disposition.md`):
- Item 4 now supplies the shared vocabulary.
- Item 13 requires general profiles.
- Parts 2 and 3 list each engine's repairs.

Reviewed until clear. Dispositions: `directives/_measurements/S9b_repair_build_directive_amend3_review_disposition.md`.

**Amendment 4** (orchestrator, 2026-10-09; a repair under the user's standing approval). The spec's amendment 1 makes
the radar comparison object the round trip's logarithmic slope `𝒮_RT` (spec observable 2). Items 4, 7, 10 and 14 and
Part 3 now name it. Reviewed until clear with the spec amendment. Dispositions:
`directives/_measurements/S9b_SHARED_PHYSICS_amend1_review_disposition.md`.

**Amendment 5** (orchestrator, 2026-10-09; a repair under the user's standing approval), for the build's repair
round 2 after review r1 (`directives/_measurements/S9b_repair_build_r1_review_disposition.md`). Items 5, 6 and 9
change, with the spec's amendment 2. Parts 2 and 3 point at them. Reviewed until clear with the spec amendment.
Dispositions: `directives/_measurements/S9b_SHARED_PHYSICS_amend2_review_disposition.md`.

**Amendment 6** (orchestrator, 2026-10-10; a repair under the user's standing approval), for the build's repair
round 3 after review r2 (`directives/_measurements/S9b_repair_build_r2_review_disposition.md`). Item 14 states what
a difference shows for a payload that carries a domain. Parts 2 and 3 point at the items each engine's repair must
meet. Reviewed until clear. Dispositions:
`directives/_measurements/S9b_repair_build_directive_amend6_review_disposition.md`.

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

   **Shared vocabulary (amendment 3).** Each object below is emitted by both engines under exactly this
   `<QUANTITY>` name. Every other emitted object is engine-local (`_LOCAL_`).

   Tokens:
   - `G<a><b><c>`: one retained optical monomial of the spec's "Order" item, where `a`, `b` and `c` are its
     exponents of `δ`, `V/c₀` and `(∂ξ_w)²`. There are twelve grades, `G000` to `G121`.
   - `<OBS>`: `DEFLECTION` (Δθ), or `RADAR` (the round trip's logarithmic slope `𝒮_RT`, item 7).
   - `<RESP>`: the Part C responses `CONSTANT` (`c_γ ≡ c₀`), `FIXED_RATIO` and `POWER`.

   Names for the spec's objects:

   | Name | Object |
   |---|---|
   | `BRANCH_EXISTENCE`, `PATH_TRAVERSAL`, `BRANCH_TYPE` | the two conditions and the branch type of the spec's "Branch existence" |
   | `A_DEFLECTION_G<abc>` | Δθ, at that grade |
   | `A_ROUND_TRIP_G<abc>` | the round-trip excess time, at that grade |
   | `A_ONE_WAY_ER_G<abc>`, `A_ONE_WAY_RE_G<abc>` | the one-way excess time from emitter to reflector, and from reflector to emitter, at that grade |
   | `A_NONRECIPROCAL_G<abc>` | their nonreciprocal part, at that grade: half of (emitter-to-reflector time minus reflector-to-emitter time) |
   | `A_NONRECIPROCAL_PATH_DEPENDENCE` | item 8's object |
   | `A_RADAR_LOG_G<abc>` | item 7's slope, for each grade in Part B's comparison sum (`G100`, `G010`, `G020`, `G001`) |
   | `B_GAMMA_<OBS>`, `B_GAMMA_DIFFERENCE` | Part B's effective `γ`s, and their difference: the deflection `γ` minus the radar `γ` |
   | `B_RESIDUAL_<OBS>` | Part B's residuals against the references |
   | `B_CONDITION_<OBS>`, `B_IMPLIED_JN_<OBS>` | Part B's every-`b` condition, and the `j_n` it implies |
   | `C_<RESP>_CONDITION_<OBS>`, `C_<RESP>_IMPLIED_JN_<OBS>` | Part C's rewritten condition, and the `j_n` it implies |
   | `C_<RESP>_N_DEPENDENCE` | where `n` enters that row's conditions |

   Item 10's objects are named by a prefix on the name of the object they restrict:
   - `R_FLOW_ONLY_`, `R_SPEED_ONLY_` and `R_TILT_ONLY_` prefix `A_DEFLECTION_G<abc>`, `A_RADAR_LOG_G<abc>`,
     `B_GAMMA_<OBS>`, `B_GAMMA_DIFFERENCE`, `B_RESIDUAL_<OBS>`, `B_CONDITION_<OBS>` and `B_IMPLIED_JN_<OBS>`.
   - `R_BULK_ONLY_` prefixes `C_FIXED_RATIO_CONDITION_<OBS>`, `C_FIXED_RATIO_IMPLIED_JN_<OBS>`,
     `C_POWER_CONDITION_<OBS>` and `C_POWER_IMPLIED_JN_<OBS>`.
   - `F_FLOW_` (the forward case's first stage, `δ ≡ 0`, `ξ_w ≡ 0`) and `F_LIVE_` (its second stage) prefix
     `BRANCH_EXISTENCE`, `PATH_TRAVERSAL`, `BRANCH_TYPE`, every graded `A_` name (`A_…_G<abc>`),
     `B_GAMMA_<OBS>`, `B_GAMMA_DIFFERENCE`, `B_RESIDUAL_<OBS>` and `B_CONDITION_<OBS>`. `F_LIVE_` also prefixes
     `C_<RESP>_CONDITION_<OBS>`.
   - `F_MASS_SOLUTION` is the forward case's solved `V`, expressed through `Φ` and `ρ_br`.
5. **No verdicts (D6 B5).** No `VERDICT`, `PASS` or `FAIL`. A boolean-valued test is emitted as the CAS object the
   test returned. Emission never depends on a payload's value. Outside the branch-existence conditions the spec
   requires `NOT_ESTABLISHED`. That output, and the branch type, are produced from the computed conditions, never
   typed. The branch-existence condition is computed through `c_γ²` without presupposing its sign (spec, "Branch
   existence"), so each branch type the spec lists can be produced from it.
6. **The every-`b` conditions are reduced (D6 B2).**
   - Each Part B and Part C condition is solved relative to `GM` for every far-zone `b`. The solving is case by
     case, over every branch the reduction produces, and each case carries its domain.
   - An unevaluated quantified set, `ConditionSet` or `ForAll` is a restatement, not a reduction.
   - For each reduced condition, print the `j_n` it implies through `∇·(ρ_br V) = −j_n`, with `ρ_br` symbolic and
     the sign of `V` symbolic. The implied `j_n` is a consequence of its own condition, not the supplied balance
     with `V` left free:
     - On every branch of the reduction, `V` is the velocity that condition determines, given the other live
       profiles and `GM`, and `j_n` is the balance evaluated on that `V`.
     - Every branch carries its domain, including any integration constants. A branch excluded by a supplied
       premise names that premise.
     - Where a condition leaves `V` undetermined, the print says so and names what stays free.
   - The `γ` difference is reduced by the engine's CAS, so its printed form shows its value on each stratum. It is
     not printed as an uncombined difference of two objects.
   - The method is the builder's choice (item 13). If an engine cannot reduce a condition, item 12 applies.
7. **The radar slope (D6 B3).** The spec's logarithmic slope `𝒮_RT` (observable 2) is obtained by computation from
   Part A's computed round-trip excess time. It is not a rule applied by hand. Part B's radar `γ` uses that slope.
8. **Path dependence (D6 B4).** Whether the nonreciprocal part depends on the path or only on the endpoints is
   computed on the nonreciprocal object Part A computed, not on a placeholder.
9. **Part C domains (D6 B6).** Each Part C row's domain predicates are built from that row's own response. They are
   not built from another row, or from the amplitude before the response is applied. Each predicate follows from the
   row's response being real and positive, or from the row's supplied inputs. The first-order counting in `f` is an
   order, not a numeric bound on any amplitude.
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
      - Print the deflection, the round trip's logarithmic slope (item 7), both effective `γ`s, their
        difference and both Part B residuals against the references, each as a function of `b`. Also print the round-trip excess time, both
        one-way excess times and their nonreciprocal part, with `b`, `Z_E` and `Z_R` kept.
      - At both stages, also print Part B's condition for each observable under the same substitution, reduced as
        in item 6. Its free inputs are `Φ` and the profiles live at that stage: `ρ_br`, plus `δ` and `ξ_w` at the
        second stage. `j_n ≡ 0` is the case's premise, not a free input, so item 6's implied-`j_n` print does not
        apply. The condition is a computed object. It is not a supplied relation between `Φ` and `GM`.
      - Also print the second-stage condition rewritten under each of Part C's three responses, as Part C rewrites
        Part B's conditions. Reduce each as in item 6 and build its domain as in item 9. `Φ`, `ρ_br`, `f`, `ξ_w` and
        the row's `n` or `s` stay symbolic. These are computed objects too, and the implied-`j_n` print does not
        apply to them.
      - Label every object with the premise.
    - For the deflection and for the radar slope (item 7), print each restriction below of:
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
    - Run every CAS demonstration through `scripts/s11c_guarded_run.py` with `--pool s9b --memory-gib 16`. Wolfram
      runs also take `--tasks-max 64`. These limits rest on the user's authorizations of 2026-10-06 and 2026-10-07,
      and on the user's raise to 16 GiB on 2026-10-08, after the SymPy harness was killed at 8 GiB.
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
    - The representation of general live radial profiles.

    **General profiles (amendment 3; the user's choice, 2026-10-09).**
    - In every shared-vocabulary object (item 4), each profile the object leaves live is a general radial
      function. The profiles are `δ`, `V`, `ξ_w`, and, where they enter, `ρ_br` and `f`. Item 10's restrictions
      and forward stages, and Part C's responses, apply their substitutions first.
    - Part A's observables are functionals of the live profiles. Each every-`b` condition is reduced, as item 6
      requires, to a condition on them.
    - No profile family is used for a shared-vocabulary object. An engine may also print family-restricted
      objects, as engine-local tags whose labels name the family.
    - If an engine cannot reduce a Parts A–C condition for general profiles, item 12 applies. If it cannot reduce
      an item 10 condition, that is not a stop event (item 10): it omits the object and reports it under item 17.
14. **Ablation harness** (`docs/development_pipeline.md` §4).
    - Each engine gets a committed harness that runs the live engine unchanged as the baseline. For each knife
      below, it runs a copy with exactly that one mutation at the named construction site.
    - It prints the baseline payload, the corrupted payload and their difference for **every** tag, including tags
      whose difference is zero. It prints no verdict.
    - For a payload that carries a domain, the difference shows separately whether the value changed and whether
      the domain changed. The domain is the condition a payload attaches to its computed object, such as the spec's
      branch-existence conditions or the condition that a case of a reduction carries. The value is the computed
      object itself, including when that object is a relation, a predicate such as a closedness test, or a
      condition. Each part's difference is computed by the engine, so an unchanged part prints as zero. This holds
      for K11's compact differences and method residual too. K11's sample domain is coverage, not a domain in this
      sense.
    - A knife whose construction an engine does not contain is reported, not invented. The engine itself contains
      no control with an expected outcome.

    The knives:
    - **K1, advection:** in the supplied kinetic symbol, `ω − V^i k_i` is replaced by `ω`.
    - **K2, metric:** `g^{ij}` is replaced by `δ^{ij}`.
    - **K3, local speed:** `c_γ(x)` in the dispersion relation is replaced by `c₀`.
    - **K4a–c, profile freezes:** each of the Parts A–C radial profiles `δ`, `V` and `ξ_w` is replaced by a constant
      where it is differentiated in the ray construction (three separate knives).
    - **K5, return leg:** the round trip's return leg propagates in the outgoing leg's direction.
    - **K6, radar source (a named coefficient knife):** the object whose slope is taken (item 7) is replaced
      by one one-way excess time. It is named for a channel no FORM knife sees: which computed
      object the slope is taken from (`docs/development_pipeline.md` §4). K1–K3 are the FORM knives on the
      round-trip construction that feeds the slope.
    - **K7, quantifier:** the every-`b` reduction is carried out at one fixed far-zone `b`.
    - **K8, mass balance:** `ρ_br` in the supplied mass balance is replaced by a constant.
    - **K9a–b, Part C responses:** in the fixed-ratio response, the local `c_s(x)` is replaced by `c_s0`; in the
      power response, `ρ(x)` is replaced by `ρ₀` (two separate knives).
    - **K10, branch gate:** the computed branch-existence condition is replaced by its logical negation before the
      observables are gated on it.
    - **K11, non-radial flow:** at the advection construction, the advected velocity gains an azimuthal component,
      so it is no longer radial. It acts on the construction from which item 8 computes path dependence.

    **K11 evaluation (amendment 2; the user's choice, 2026-10-08).** The guard killed K11's full copy at 8 GiB and
    again at 16 GiB, both times after the same 31 tags. For K11, the corrupted payload is a compact evaluation by a
    method of the builder's choice, and the full K11 copy is not run. Item 11's stop-and-report rule applies if the
    guard kills the compact run.
    - The method calls the engine's own constructions, or is mechanically extracted from them. The K11 mutation is
      its only source change, and it stays in force in the evaluated copy.
    - The same method, with the same configuration, is applied to an unmutated copy. The configuration includes
      every sample point and seed, every truncation order, and the arithmetic's precision.
    - For every tag, the harness prints:
      - the full baseline payload;
      - both compact evaluations and their difference;
      - the method's residual: the unmutated compact evaluation minus the full baseline payload reduced by the same
        configuration.
    - Each difference and residual is printed with its coverage: the retained grades it keeps, and its truncation
      order or sample domain. If the arithmetic is not exact, each also carries its precision and an error bound.
    - A tag the method does not reach is printed as not evaluated.
    - The engine, its baseline run and the other knives are unchanged.
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
    - every stop-and-report event (item 12);
    - for K11, the compact method, its coverage, and every tag printed as not evaluated;
    - every shared-vocabulary name (item 4) the engine does not emit, and why.

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
- **Repair round 1 (amendment 3).** Besides items 4 and 13, these must be true after the repair:
  - Item 8's nonreciprocal one-form is computed from this engine's own dispersion relation, with the advected
    velocity entering it as a vector field. K11's mutation then enters the dispersion.
  - Each item 10 restriction is substituted into its gates as well as into its observables.
- **Repair round 2 (amendment 5).** Items 5, 6 and 9 as amended, and the spec's amendment 2, hold after the repair.
- **Repair round 3 (amendment 6).** The branch domain is the set where both conditions in the spec's "Branch
  existence" paragraph hold. On every stratum inside it, including `GM = 0`, each `γ` object (`B_GAMMA_DEFLECTION`,
  `B_GAMMA_RADAR`, `B_GAMMA_DIFFERENCE`), and each of its restriction and forward copies, prints what that object's
  relation determines there, computed by the engine. `NOT_ESTABLISHED` is printed only outside the branch domain
  (item 5).

## Part 3. Wolfram engine (blind)

- **Repair and extend** `research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl`. **Write** its
  harness `research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl`. The engine's only write is
  its stdout tag stream.
- **This is the blind engine.** It reads no file at run time, imports nothing, and re-derives everything from the
  spec's equations. It imports no `LEDGER` and writes none.
- **Mechanical precedent, not authority:** `mathematica/S11b_interface_coupling_law_mathematica_audit.wl`, for its
  emit, naming and flush shape only.
- **At most one kernel at a time.**
- **Repair round 1 (amendment 3).** Besides items 4 and 13, these must be true after the repair:
  - Every along-ray gate quantifies over a constructed domain, the radii the ray traverses. No gate quantifies
    over an undefined head.
  - The radar `γ` is solved on every stratum of the radar slope (item 7). The `γ` difference, and each
    restriction and forward copy of the radar `γ`, use that solution.
- **Repair round 2 (amendment 5).** Items 5, 6 and 9 as amended, and the spec's amendment 2, hold after the repair.
- **Repair round 3 (amendment 6).** Item 14 as amended holds for the harness.
- **Executable checks, with no expected values:**
  1. Copy the finished `.wl` alone into an empty scratch directory. Run it there and in the repository. Both runs
     exit 0, and their streams are non-empty and byte-identical. Afterwards the scratch directory holds only the
     copied file.
  2. With stdout redirected, the capture file is observed growing while the kernel is alive. If every task finishes
     too fast to observe this, report that.
  3. `git status --porcelain research/pde_ledger_v3/mathematica/out/` prints nothing.
