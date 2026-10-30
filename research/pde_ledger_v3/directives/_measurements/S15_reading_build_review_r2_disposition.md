# S15 reading build, review round 2 (round 3 of the build): dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S15_reading_build_review_r2_lookups.md`, beside this file
(generator `_scratch/s15/gen/s15_reading_r2_lookups.sh`). "M-n" cites a section of that file. Every computed result
quoted here is a leg's saved stdout; the legs' scripts were not re-run. The legs' closed forms and conditions stay in
their scratch files. This file and the lookups cite only their residuals, the strata they did not finish, and
existence counts.

## The round

- **Artifact.** Repair 2 of the S15 reading build, uncommitted on `6c0fc66b` (the round 1 baseline). The original
  author (the same Codex session) repaired both engines and the run record for W1 and W2 of
  `S15_reading_build_review_r1_disposition.md`. The build is:
  - `scripts/S15_separation_departure.py` (`ccd4324e…`) and `mathematica/S15_separation_departure.wl` (`ff784cac…`);
  - their transcripts, about 1.2 GB (`138391ad…`) and 1.36 GB (`b133fdcb…`), written by the guarded runs
    `reading_w1_w2_sympy_final_01` (exit 0, 12545 s) and `reading_w1_w2_wolfram_final_02` (exit 0, 1609 s);
  - `directives/_measurements/S15_separation_departure_run_record.md` (`fe980d36…`).

  All five matched the r2 baseline after the round, and the working tree was unchanged apart from this file and the
  lookups (M0). One stray file was in the repository root and was moved before the lookups ran (Process).
- **Prompt and legs.** `_scratch/s15/reading_review/s15_reading_build_review_prompt_r2.md`, sha256 `be51c098…` (M0).
  It is the round 1 prompt with a repair item 0: W1's and W2's "What must be true", the labels, the claims,
  dependence ablations, and whether items 1, 2, 5 and 6 and the item-3 objects outside W2 changed. Items 3, 4 and 7–12
  applied in full to the repaired code. Both legs received it identically, in parallel:
  - a fresh Claude agent (Opus): `/var/projects/_review_scratch/s15-reading-r2-claude/REVIEW.md`;
  - Grok (`grok-4.7`, effort high): `/var/projects/_review_scratch/s15-reading-r2-grok/REVIEW.md`, exit 0.

  The orchestrator recorded each report's hash before reading either report, and both still match (M0).
- **Verdicts.** Both NEEDS REVISION, with one finding, which each leg reached on its own.

## Findings

**X1. For the 48 objects that depend on the separation, item 4 prints the object's own vanishing in another form,
not the conditions.** (Claude F1, Grok 1.) **CONFIRMED.** W1 is not repaired.
- **The contract** (M1). Part B item 4: "the conditions under which that object vanishes for every `D` on the domain,
  exact, with every special exponent value treated. A condition is on the exponents and, where they enter, on the
  amplitudes, the radii and their rates, and the velocities. It is printed as equations and inequalities in the named
  quantities; `Q` on the family under each condition." W1's "What must be true" and V2's say the same for these 48
  objects.
- **The engines** (M2):
  - SymPy sends each direction coefficient of a remainder that holds `d` to `finite_identity_condition`
    (`scripts/S15_separation_departure.py:1093–1102`). That function prints
    `Eq(Sum(Abs(reference**power*scalar)**2,(power,0,bound)),0)`, where `scalar` is an entry of an aliased ordered
    product of transfer matrices (`:293–299`, `:303`). The object's condition is the conjunction of these equations
    with `R > r_in > 0` (`:1103–1104`). `Q_under_condition` is an Euler template filled with the closed `Q` at
    `d = 2R` (`:1105–1108`).
  - Wolfram has the same construction: an `Inactive[Sum]` of squared product entries `== 0` (`.wl:610–611`), the
    conjunction (`:662`), and `Q_under_condition` from a restriction template at the same reference (`:663–666`).
  - Neither engine takes the step from these equations to values or relations of the named quantities.
  - That the finite equations are equivalent to the identity is stated in a docstring (`.py:30–36`), a comment
    (`.wl:613–623`) and a printed prose certificate in each transcript (M3). It is not computed.
  - The construction is new in repair 2. The matching lines go from 0 at `6c0fc66b` to 9 (SymPy) and 5 (Wolfram)
    (M2).
- **The transcripts** (M3), the same in both engines:
  - **Which objects.** 48 objects carry the finite-identity status. They are three components each of 16 groups: the
    drain-1 and mixed groups of the carried current, the centre transport and the radius transport at both surfaces,
    and of the region momentum and rate.
  - **Their lines.** Each has one `condition` line and one `Q_under_condition` line. None has a numbered or union
    condition.
  - **The conditions.** All 48 condition lines contain the jet sum, and every equation in them is a jet sum. None
    contains an exponent, an amplitude or a derivative symbol (0 of 48 for each of three patterns, in each engine).
    The sums' upper bounds reach 825513 (SymPy) and 2249810 (Wolfram).
  - **Labels.** `NOT OBTAINED` and `EXACT_OBJECT_GAP` occur 0 times in each transcript.
  - **Where `d` is absent,** the same engines print solved conditions in the amplitudes and velocities (M3).
- **The printed equation is the object's vanishing** (M4, Claude). Claude evaluated each engine's printed transfer,
  projection and initial tables. The terms `c_n` equal the Taylor coefficients at `2R` of `H` times the printed
  direction coefficient:
  - SymPy: residual/scale `5.6e-42` at `n = 0` (`e09`, 42 lines at two exponent points). Where the ratio reaches
    `1.0`, the scalar is `0` in the transfer route and of order `1e-47` in every direct term, so the scale there is
    round-off.
  - Wolfram: agreement at 20 digits for `n ≤ 6` (`e10e`).

  The sum is therefore zero exactly when the object vanishes for every `D`. It is the question in another form.
- **`Q` under the condition** (M3). There is one line per object: `Q` at `2R`, carried by the Euler template. It is
  the same expression whichever condition holds, so it does not give `Q` under each condition. Its operand is the
  unrestricted closed form at `d = 2R`, which still carries the general-exponent Piecewise.
- **The conditions are reachable** (M4). Each leg reached explicit conditions for representative objects by exact
  routes of its own, and checked them by a second route:
  - **Claude:** `d01_derive.py` and `d01c_momentum_rate.py`, with direct quadrature in `d02b_quadrature.py`. The
    largest scaled residuals run from `2.9e-07` to `4.9e-15`, each within that route's convergence estimate. One
    stratum of the momentum and rate did not finish within 900 s, and Claude reports it as not obtained.
  - **Grok:** `derive_item34.py`, `derive_special.py` and `derive_conditions.py`, using two routes plus quadrature: 16
    residuals `abs=0.0` and one `2.9e-33`. Grok's last run was cut off (Process); what it flushed is used.
- **Claims** (M5). The run record says:
  - "48 separation-dependent objects with finite-identity conditions" (l.15, l.30);
  - "Requested-object gaps indexed by this final engine: none" (l.18, l.33);
  - "Each engine solves its computed geometric scalar equation … to obtain Q under each condition" (l.39).

  The engines print those equations and do not solve them. The certificate lines state an equivalence that no
  computation in the engine produced (E1).
- **Effect.**
  - The reading Part B needs is in neither transcript: for which profiles each of these objects vanishes, and so
    where `Q` falls off as the inverse square.
  - Both engines use one construction: the same bound formula, the `2R` reference and the Euler template. On item 4
    they cannot disagree, and the comparator would compare two restatements.
- **Author (G4).** This defect class has now been filed in three rounds, each time on item 4 of these objects:
  - round 0, V2: the object `== 0`;
  - round 1, W1: solver returns discarded;
  - this round: the object's jet.

  This round's form sits in material repair 2 wrote. The repair therefore goes to a new author.
- **What must be true.**
  - For each of the 48 objects, in each engine, the conditions under which the item-3 object vanishes for every `D`
    on the domain:
    - computed exactly from the printed closed form;
    - complete, with every special exponent value treated;
    - printed as equations and inequalities in the named quantities: values of, or relations among, the exponents;
      and equations or inequalities in the amplitudes, the radii and their rates, and the velocities.
  - `Q` on the family under each condition: one per condition, computed on that condition.
  - An equation that sets the object, its coefficients, its derivatives or its jet to zero restates the question. It
    may be printed as an intermediate object, and it does not replace the conditions.
  - A not-obtained label appears only where an exact attempt that could have returned a condition did not. The label
    says what was attempted and stands beside what was printed.
  - No printed line or run-record sentence states a property that no computation in that engine produced: an
    equivalence, a bound, completeness, or "solves".
  - Each engine reaches its conditions by its own construction.
  - No object of items 1, 2, 3, 5 or 6 changes by value.

## Examined and not filed

- **Reruns** (M8):
  - **SymPy.** Claude's rerun (`claude-rerun_py2`, exit 0, 12708 s, `93d3eba1…`) has every item-4 tag identical
    (14040). 21 tags differ in bytes:
    - 6 item-2 pressure `whole_minus_expansion` tags;
    - 15 item-3 tags.

    All 21 are equal by value: 12 parse equal, and 9 are zero by both of Claude's routes (`x11b`). The differences are
    print order, and `polylog(0, r_in/d)` against its rational form. This extends comparator note b: compare SymPy by
    value, allowing `polylog(0, x)` and `x/(1−x)`.
  - **Wolfram.** Claude's rerun (`claude-rerun_wl`, exit 0, 1513 s) differs in 1848 item-4 tags. Each is identical
    after the generated symbols are renumbered (`x09`). This extends comparator note c: Wolfram's generated names
    (`Unique` aliases, `Module` temporaries) vary between runs, so compare after renaming them.
  - The run record does not claim that a rerun reproduces either transcript's bytes.
- **Grok's controls.** Grok ran no rerun, FORM ablation or one-sided corruption. Its derivation job held its only
  guarded slot, and Grok left it running (Process). Claude's controls cover item 12 for this round (M7).
- **The multiplicity bound itself** was not tested by either leg. It does not bear on X1, which holds whether or not
  the bound is right.

## Process

- **Grok's orphaned run (13:2x on 10-29).** Grok exited at 13:22 with its guarded job `grok-derive-04` (unit
  `s11c-guard-fcb350549c72`) still running, and that run's launcher had died. The guard never reclaims such a
  reservation by itself, and this one blocked Claude's admissions. The orchestrator handled it as follows:
  - It confirmed that no Grok process was alive and that the owner pid was dead.
  - It stopped the unit with `systemctl --user stop`, which left it inactive with no queued job.
  - It saved the registry as `_scratch/s15/reading_review/r2_registry_before_reconcile_grok-derive-04.json` (`29f5bb25…`;
    two records, `claude-d01` and `grok-derive-04`; M0).
  - After a dry run, it removed only Grok's record under the admission lock.

  `grok-derive-04` therefore has no `outcome.json` (M0). Carried to the next leg prompt: a one-shot leg's guarded run
  can outlive the leg, so check for it after each such leg ends.
- **The Claude leg's deadlock, and an amendment (19:49 on 10-29).**
  - **The deadlock.** `claude-d01` (20 GiB, no output since 12:08) and `claude-e10` (25 GiB, no output since 18:34)
    held enough of the pool that the leg's SymPy rerun could not be admitted. The prompt said not to stop guarded
    runs.
  - **The amendment.** The orchestrator sent the leg an operational note by `SendMessage`. Its text is verbatim in
    `_scratch/s15/reading_review/r2_claude_amendment_sendmessage.json` (M0). It carried no findings and no other
    leg's output. It let the leg stop its own runs that it judged stalled, with `systemctl --user stop` while the
    launcher was alive, then confirm that the unit was inactive and its registry entry gone. It did not let the leg
    edit the registry or touch other runs, and the decision stayed with the leg.
  - **Leg symmetry.** Grok had finished six hours earlier, so the amendment changed no leg's method relative to the
    other.
  - **What the leg stopped.** `claude-d01`, `claude-e10`, `claude-e10b`, `claude-d01b` and `claude-py_depB`. Each
    outcome records exit 1 (M0). Their replacements finished, and the registry is `{}` (M0).
  - **Carried to the next leg prompt:** let a leg stop its own stalled run in this order; size `--memory-gib` to the
    job, so that one leg's jobs do not block each other; never wrap the guard launcher in a shell `timeout`.
- **The Claude leg's other disclosures:**
  - `claude-d01` ran at 20 GiB beside the leg's other guarded jobs (M0).
  - `x04` first ran outside the guard. Its guarded rerun (`claude-x04g`) is byte-identical (M6).
  - `claude-x11` died on the leg's own import error and was rerun as `claude-x11b`.
- **A stray file in the repository root.** `x02_PYFORMA.stdout` was untracked and empty, with mtime 00:47:36. The
  leg's own `x02_PYFORMA.stdout`, in its scratch directory, was written 11 s later. The orchestrator moved the stray
  into `…/s15-reading-r2-claude/stray_repo_root/`, keeping its mtime, before the lookups ran (M0).

## What the legs checked and found sound

- **W2 is repaired** (M6).
  - **Points.** Claude evaluated both engines' item-3 object and `Q` operand for all 48 objects at 1,083 points. These
    include every divisor locus and Piecewise condition either engine prints, and every integer exponent stratum
    from −6 to 2.
  - **Results.** The largest scaled residual is `1.361e-03` in both engines, with nothing unevaluable; round 1's
    Wolfram had 30 unevaluable points. No row is beyond 10 times the reference's own convergence.
  - **The W2 divisors.** They occur 0 times in the six lines W2 named. Grok's audit finds every printed locus answer
    `False` or empty.
  - **Status lines.** The status lines of those tags now describe only the construction.
- **Items 1, 2, 5 and 6 are unchanged from round 1** (M6):
  - SymPy: 672, 3348, 4704 and 3246 tags are the same. The 12 item-2 pressure `whole_minus_expansion` tags differ in
    text only (`expand(r1 − r2) = 0`).
  - Wolfram: 672, 3360, 5376 and 3242 tags are the same.
  - Grok's hash comparison agrees.
- **Item 4 where `d` is absent** is unchanged. Wolfram's differing tags are identical after renaming the generated
  symbols (Claude `x03`).
- **W1 and W2** restate Part B items 3 and 4 and the label rule, and add nothing (both legs).
- **Items 7, 8 and 10** (Claude):
  - The conditions carry `R > r_in > 0`, and the reference `2R` is on the domain `d > R`.
  - No law, coefficient or normalization is introduced.
  - The `.wl` reads nothing from the SymPy side. The shared item-4 method is not transcribed line by line.
- **Controls** (M7, Claude). Every ablation ran on a copy, under the guard, and each moved its target:
  - **WL FORM** (`wl_formA`, `R` scaled with the separation): 51 item-3 tags and the item-4 tags of the same 33 objects
    move; items 1, 2, 5 and 6 are unchanged.
  - **PY FORM** (`py_formA`, the power primitive's special branch dropped, on a three-object subset): item 3 and item 4
    move for both separation-dependent objects.
  - **PY FORM** (`py_famA`, `R` scaled): item 3 and item 4 move for all three objects.
  - **Item 4 depends on item 3** in both engines (`wl_depB`, `py_depC`): item 4 moves exactly where item 3 does.
  - **WL one-sided balance corruption** (`wl_balC`, the outer product made isotropic): construction, reading `Q`,
    items 2, 3, 4 and 5 move, and item 6 does not.

## Next

1. Preserve this reviewed baseline, with the transcripts annexed to GIN.
2. Repair X1 by a new author (G4): a fresh Codex session. The brief carries Parts A and B, the engines as preserved,
   and the "What must be true" lines above. It carries none of the legs' closed forms or conditions.
3. Then round 4 of this build, in round 3's scoped shape, with the carried leg-prompt notes.

The premise was checked before the author change.
- The question is the step's own reading: which profiles give an inverse-square fall-off. The user chose this reading
  build on 10-26.
- It rests on no freeze: exponents stay symbolic, and amplitudes, radii and velocities stay live.
- Both legs reached its answer for representative objects by exact routes. What has failed three times is the
  construction, not the question.
