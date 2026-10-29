# S15 reading build, review round 1: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S15_reading_build_review_r1_lookups.md`, beside this file
(generator `_scratch/s15/gen/s15_reading_r1_lookups.sh`). "M-n" cites a section of that file. Every computed result
quoted here is a leg's saved stdout; the legs' scripts were not re-run. The legs' closed forms and conditions stay in
their scratch files; this file and the lookups cite only their residuals, probes on toy expressions, and existence
counts.

## The round

- **Artifact.** Repair 1 of the S15 reading build, uncommitted on `5f0fc03a` (the round 0 baseline). The original
  author (the same Codex session) repaired both engines and the run record for V1 and V2 of
  `S15_reading_build_review_r0_disposition.md`. The build is:
  - `scripts/S15_separation_departure.py` (`1a6a839d…`) and `mathematica/S15_separation_departure.wl` (`407a40d4…`);
  - their transcripts, about 376 MB (`c88240f7…`) and 484 MB (`9495fdb1…`), written by the guarded runs
    `reading_repair_sympy_final_04` and `reading_repair_wolfram_final_01`, both exit 0;
  - `directives/_measurements/S15_separation_departure_run_record.md` (`990e38d0…`).

  All five matched the r1 baseline after the round. The working tree was unchanged apart from this file and the
  lookups, and no file outside `_scratch/` is newer than the prompt (M0).
- **Prompt and legs.** `_scratch/s15/reading_review/s15_reading_build_review_prompt_r1.md`, sha256 `dc3e4cde…` (M0):
  the round 0 prompt, with a repair item 0 (V1's and V2's "What must be true", the labels, the claims, dependence
  ablations, and whether items 1, 2, 5 and 6 changed). Items 3, 4 and 7–12 applied in full to the repaired code. Both
  legs received it identically, in parallel:
  - a fresh Claude agent (Opus): `/var/projects/_review_scratch/s15-reading-r1-claude/REVIEW.md`;
  - Grok (`grok-4.7`, effort high): `/var/projects/_review_scratch/s15-reading-r1-grok/REVIEW.md`, exit 0.

  The orchestrator recorded both reports' hashes before reading them, and they still match (M0).
- **Verdicts.** Both NEEDS REVISION. Each leg reached W1 on its own; W2 is Claude's alone.

## Findings

**W1. Item 4 is not computed for any closed form that depends on the separation, and neither engine's attempt can
return it.** (Claude 1, Grok 1.) **CONFIRMED.** V2 is not repaired.
- **The contract** (M1): Part B item 4, "the conditions under which that object vanishes for every `D` on the domain,
  exact, with every special exponent value treated. … printed as equations and inequalities in the named quantities;
  `Q` on the family under each condition." V2: "Part B item 4 holds for every closed form that item 3 prints,
  including those that depend on the separation." V1: "An object is labelled as not obtained only where that engine
  attempted it exactly and the attempt did not return it."
- **The engines** (M2):
  - SymPy abbreviates the `d`-dependent kernels before the solver (`scripts/S15_separation_departure.py:776`) and
    calls `solve_undetermined_coeffs` on each coefficient (`:782`). Its return is printed and collected in
    `attempted` (`:785–787`), which nothing else reads. Then `complete=(not reduced.has(d,…))` (`:788`), and for a
    remainder that holds `d`, `else: solutions=[]` (`:820`). Whatever the solver returns, no condition is printed.
  - Wolfram poses `SolveAlways` on `Implies[dist>rout[t], …]` (`mathematica/S15_separation_departure.wl:492`, `:499`).
    A returned list is printed as candidate conditions and kept as `answer` (`:501–507`). Then
    `If[… || !FreeQ[answer,dist] || ListQ[answer], <NOT OBTAINED>]` (`:510`) files both an unevaluated return and a
    returned list as not obtained.
- **The transcripts** (M3), the same counts in both engines:

  | | Item 3 closed form | Item 3 general field | Item 4 computed | Item 4 not obtained | Item 4 not applicable |
  |---|---|---|---|---|---|
  | SymPy | 144 | 192 | 96 | 48 | 192 |
  | Wolfram | 144 | 192 | 96 | 48 | 192 |

  - Each engine's `EXACT_OBJECT_GAP` index has 48 entries, all item 4. These are the drain-1 and mixed groups of the
    carried current, centre and radius transport, and region momentum and rate (both legs).
  - SymPy prints 258 `scalar_identity_solver_return` lines. 162 end in `None`, all in groups with drain 1.
  - Wolfram prints 48 `identity_solver_return` lines, each an unevaluated `SolveAlways[Implies[…]]`, and no
    `identity_candidate_conditions` line.
- **The attempts cannot return the object** (M4). These are Claude's probes on toy expressions:
  - SymPy's solver returns solutions for a coefficient rational in `d`, and the engine's gate on the same operand
    gives `False`. With a `d`-power abbreviated to a kernel function, the solver returns `None`.
  - Wolfram's `SolveAlways` on the `Implies` form returns unevaluated. On the plain equation it returns a list.
- **The conditions are reachable exactly** (M4). Each leg reached them by its own exact route and checked them
  against a second route. Representative objects:
  - **Claude:** the carried current, centre and radius transport, and region momentum and rate of the drain-1
    group (`d02_derive.py`, six exponent sets). They were checked against direct quadrature (`d03_quadrature.py`), with
    every residual at most `1.21e-13`.
  - **Grok:** the carried current, radius transport and region momentum (`d02_conditions.py`, `d03_vanishing.py`).
    They were checked on the engines' own printed objects in both engines (`d06`, `d07`, `d08`, `d08b`):
    - 24 residuals are exactly `0`;
    - four Wolfram residuals are `0` at 65–68 digits;
    - three are at most `2.2e-137`.
  - Grok's scope note: the inner surface, the other components and the mixed groups carry the same gate and were not
    given separate residuals.
- **The labels.** Each item-4 label says what was attempted. But each attempt was built so that it could not return
  conditions: SymPy discards any return once `d` remains, and Wolfram discards a returned list. The labels describe
  gates, not failed attempts, so they do not discharge item 4.
- **Effect.**
  - The build prints no condition, and no `Q` under a condition, for any object that depends on `d`. The reading
    Part B asks for cannot be made from either transcript.
  - Both engines decline the same 48 objects, so on item 4 they cannot disagree.
- **What must be true.**
  - Part B item 4 holds for each of the 48 objects in each engine. The conditions are computed exactly from the
    printed closed form, with every special exponent value treated. They are printed as equations and inequalities
    in the named quantities, and `Q` on the family is printed under each.
  - A return that gives conditions is not filed as not obtained. A label appears only where an exact attempt that
    could have returned the conditions did not, and the label says what was attempted, beside what was printed.
  - Each engine reaches its conditions by its own construction.
  - No item-1, item-2, item-5 or item-6 object changes, and no item-3 object changes except as W2 requires.

**W2. Wolfram's closed form for the region rate of the drain-1 group divides by zero at exponent values its branches
do not treat.** (Claude 2.) **CONFIRMED.** Grok did not compare the region rate (its scope note).
- **The contract** (M1):
  - Part B item 3: "with its dependence on `d` explicit and exact".
  - Part B: "No limit or expansion is taken in place of an exact object."
  - V1: "exact on the domain … every special exponent value treated".
  - V2: "No status line or run-record sentence states a property of a printed object that no computation in that
    engine produced."
- **The transcript** (M5). The tags are `reading/region/component_{1,2,3}/rate_labels_{drain1}/item3/object` and
  `…/Q_closed_operand` (lines 21458–21462, 22853–22857 and 24223–24227).
  - Each object line contains the divisors `(2 + aa[1] + bb[1])` and `(3 + aa[1] + bb[1])` 20 times each, and
    `(4 + aa[1] + bb[1])` and `(5 + aa[1] + bb[1])` 10 times each.
  - Each `Q_closed_operand` line contains them 4, 4, 2 and 2 times.
  - The only equalities in those lines are branch conditions on `aa[1]` and `bb[1]`, in two written forms.
  - Claude located the divisors inside Piecewise default branches (`e05_locate_zero.stdout`).
- **The legs' computation** (M5, Claude):
  - `e10_scan_strata.py` tests every divisor that depends only on the exponents, across all 96 item-3 lines of the 48
    objects, in both engines:
    - Wolfram: `SUMMARY WL: strata tested 606; not evaluable 24`;
    - SymPy: `SUMMARY PY: strata tested 84; not evaluable 0`.
  - The 24 are these six tags at the four exponent strata where a printed divisor vanishes, each a
    `ZeroDivisionError`.
  - The evaluator reads Wolfram's printed text faithfully: it agrees with Wolfram's own `N[]` to `6.483e-16` over 286
    values (`e06`).
  - Approaching each stratum, the printed text tends to the reference, with residuals of order `1e-9` at an offset of
    `1e-8` (`e09`). The singularity is removable, so on the stratum itself the printed object is `0/0`.
  - At the same points, SymPy's corresponding objects match the reference (`e08`).
- **The status line** (M5). All three components' status reads "Closed form from generating-polynomial angular
  moments, distance substitution and finite-interval Appell integration, including logarithmic branches." It is the
  same, byte for byte, in Claude's W1 ablation. That ablation drops the exceptional branch of the primitive and moves
  96 item-3 tags. The sentence is typed, not derived from what was computed.
- **Effect.** At those exponent values, Wolfram's region rate of the drain-1 group is not exact on the domain. Reading
  it there would take a limit, which Part B excludes, and the comparator would meet an object one engine cannot
  evaluate.
- **What must be true.**
  - Every item-3 closed form each engine prints is exact at every real exponent value on the domain, including those
    where a printed divisor vanishes, or carries a not-obtained label there that says what was attempted.
  - Each status line states only what that engine computed for that object.

## Examined and not filed

- **Reruns** (M8):
  - SymPy: both legs' unmodified reruns (`8115416b…`) differ from the production transcript only in the 12 pressure
    `whole_minus_expansion` tags. Each residual against production is `0` (Claude `c04`, Grok `d09`). Every item-3 and
    item-4 tag is byte-identical. Carried to the reading comparator's decision list (note b): compare SymPy by value.
  - Wolfram: Claude's rerun (`80265713…`) differs only in `WORKING_ON` progress lines, which carry `Module`
    temporaries (`c06`). Carried to the reading comparator's decision list as note c: Wolfram's progress lines vary
    between runs; compare object tags, not progress lines.
  - The run record does not claim that a rerun reproduces either transcript's bytes.
- **Blind residuals** (Claude, controls):
  - SymPy's `power_primitive_derivative/residual` still prints `0` when P2 drops the log stratum of the power
    primitive, so it cannot see a dropped stratum.
  - `kernel_definition_reconstruction/residual` (258 × `0`) is zero by construction.
  - The run record cites neither (M8). Carried to the S15 record decision list: cite neither as evidence of strata or
    of a check.
- **Centre transport** (Grok). No sampled exponent cleared it. Grok states that sampling is not a proof and files
  nothing. Claude's route gives its conditions (W1).
- **Process** (M0):
  - Claude's `d01_derive.py` stayed in `sympy.series` from 00:32 to 05:32 (unit `s11c-guard-a954bd5b63c5`). The leg
    reported it, did not stop it, and finished with `d02`. After both reports were in and their hashes recorded, the
    orchestrator stopped that unit while its launcher was alive. The pool released it (registry `{}`), and the
    orchestrator then stopped the agent.
  - Claude's `e02` exceeded its 16 GiB cap after 13401 s; its 286 printed values were used.
  - Grok's `d01` exceeded its memory after flushing its surface forms; Grok did not use its integration.
  - Grok withdrew `d05_numeric.stdout` itself, because its substitution missed the transcript's symbols.
  - Grok did not rerun Wolfram while Claude's kernel was running; Claude's rerun covers it.

## What the legs checked and found sound

- **V1, item 3** (M5).
  - Both engines print 144 closed forms and no item-3 gap.
  - Claude checked 15 representative objects at 27 parameter points against direct quadrature, covering every
    exponent stratum either engine branches on:
    - SymPy: `compared 330; max scaled residual 3.258e-08`;
    - Wolfram: `compared 286; max scaled residual 1.749e-15`.
    W2 is the exception.
  - Grok's residuals on the engines' own objects are in W1.
- **Item 4 where `d` is absent.** It is computed in both engines from factor zero sets (SymPy) and `Reduce` (Wolfram),
  not restated, and `Q` is printed under each condition (both legs).
- **Items 1, 2, 5 and 6 are unchanged** (M6).
  - Wolfram: 7056 tags byte-identical (Grok).
  - SymPy: 6372 byte-identical. The 12 pressure `whole_minus_expansion` tags differ in arrangement, and each
    residual of integrands is `0` (Claude `c02`, Grok `d07`).
- **V1 and V2** restate Part B items 3 and 4 and the label rule. They add nothing and drop nothing (both legs).
- **Run record.** Its commands, guard outcomes, peaks, sha256 values and counts match the files, and it does not
  claim the 48 were solved (both legs).
- **Emission and independence.**
  - No verdict tokens.
  - Profiles are substituted only on the item-3 and item-4 path.
  - Existence conditions say sufficient.
  - The Wolfram engine reads nothing from the SymPy side and is a separate construction (both legs).
- **Ablations** (M7). All were on copies, under the guard, and each moved its target:
  - **Claude:**
    - P1 (SymPy): `Q`'s closed operand multiplied by an exponent factor. Item-4 conditions follow: 1340 item-4 tags
      move and 54 are new; items 1, 5 and 6 are unchanged.
    - P3 (SymPy, one-sided): `X_1` shifted along a permutation of `D`. Items 2, 3 and 4 move; items 1, 5 and 6 do
      not, nor do the conditions where `d` is absent.
    - W1 (Wolfram FORM): the exceptional branch dropped. 96 item-3 tags and 342 item-4 tags move; items 1, 2, 5 and
      6 are unchanged.
  - **Grok:**
    - SymPy FORM: the separation chart `N d → N d²`; plus a perturbation of a drain-2 item-3 value before the
      condition solve. The item-4 status moves.
    - Wolfram FORM: the family exponent `a_k → a_k + 1`. The drain-2 and drain-1 item-3 objects move.

## Next

1. Preserve this reviewed baseline, with the transcripts annexed to GIN.
2. Repair W1 and W2. The repair brief carries the "What must be true" lines above, and none of the legs' closed forms
   or conditions.
3. The next review is round 3 of this build, so the orchestrator asks the user how to finish it (full, scoped or
   stop).
