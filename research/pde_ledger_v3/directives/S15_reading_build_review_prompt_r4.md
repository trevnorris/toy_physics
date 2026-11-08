# Independent review: the S15 reading build, repair round 4 (Codex-written)

## Artifact
In `/var/projects/toy_physics`, uncommitted. The version reviewed in round 3 is committed at `90441b89` ("Preserve the
S15 reading build repair 3 as reviewed in round 3 (not accepted)"); round 2's at `4c4e2558`.
- To see the repair, run
  `git -C /var/projects/toy_physics diff 90441b89 -- research/pde_ledger_v3/scripts/S15_separation_departure.py research/pde_ledger_v3/mathematica/S15_separation_departure.wl research/pde_ledger_v3/mathematica/S15_separation_departure_run.py research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md`
  (`S15_separation_departure_run.py` is new in this repair).
- The round-3 transcripts are git-annex files at `90441b89`, present locally:
  `/var/projects/toy_physics/.git/annex/objects/gZ/F4/MD5E-s427868598--72ecbc788c6af3346eccc375e8857b7a.out/MD5E-s427868598--72ecbc788c6af3346eccc375e8857b7a.out`
  (SymPy) and
  `/var/projects/toy_physics/.git/annex/objects/WQ/zZ/MD5E-s472539680--3ee0dcc1d578da6add2b4c721297b136.out/MD5E-s472539680--3ee0dcc1d578da6add2b4c721297b136.out`
  (Wolfram). Read them in place; the md5 in each name is its content's md5.

The files:
- `research/pde_ledger_v3/scripts/S15_separation_departure.py`: the SymPy engine.
- `research/pde_ledger_v3/mathematica/S15_separation_departure.wl` and
  `research/pde_ledger_v3/mathematica/S15_separation_departure_run.py`: the Wolfram engine and its runner.
- `research/pde_ledger_v3/scripts/out/S15_separation_departure.out` and
  `research/pde_ledger_v3/mathematica/out/S15_separation_departure.out`: the two final transcripts.
- `research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md`, with
  `S15_separation_departure_attempt_history.md` and `.json` beside it: the builder's run record.
- The guard records of this repair's runs, `_scratch/s11c/s15/repair4_*` (`_scratch/s15/guard` links there). The final
  runs are the `repair4_decision_final_*` directories and any retry runs after them.
- The builder's working directory, `_scratch/s15/reading_build_work/repair3/` (the final runs' files are in
  `final_decision/`), and its session logs `_scratch/s15/reading_build/repair4*.log`. The last message of
  `_scratch/s15/reading_build/repair4o.log` is its final report.

Codex wrote all of these. Its contract is the one round 3 reviewed: Part A, the "Supplied equations", "The drains,
their fields and the region", "What the build computes" and "Model point" sections of
`research/pde_ledger_v3/directives/S15_decision_list.md` at `91248b67`; Part B, the "Supplied (unfalsifiable in this
build)", "The separation and the family (supplied definitions)", "The profile family", "What the build computes" and
"Model point" sections of `research/pde_ledger_v3/directives/S15_reading_decision_list.md` at `7f52b4d2`. The round-3
findings are Y1 and Y2 in `research/pde_ledger_v3/directives/_measurements/S15_reading_build_review_r3_disposition.md`.

The repair was directed by a chain of notes, each stating what must be true:
`_scratch/s15/reading_build/resume_s15_reading_repair_4.md` and `…/resume_s15_reading_repair_4d.md`, then
`research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md`, `4g`, `4h`, `4i`, `4l`, `4m`, `4n` and `4o`. Three
of them carry user decisions that change what the build computes:
- `4l` supplies the distance from drain 1 over the region and its series in `z = r/d`, as an identity.
- `4n`: six region objects that carry both drains are not attempted in the final runs; each carries a not-obtained
  label of a third kind, set by that decision. Every other item-4 object gets its final runs, and an item-4 attempt is
  labelled only for a 6-hour stop or for running out of memory as `4n` defines it.
- `4o`: SymPy attempts no item-4 object; each SymPy item-4 object carries a fourth-kind label. SymPy still prints items
  1, 2, 3, 5 and 6. Wolfram's item 4 goes on under `4n`. Item 4 therefore has one engine in this build.

The two engines are meant to build the objects independently: the Wolfram script imports nothing from the SymPy script
or its output and re-derives every object. Comparing the two engines is a later, separate step, so do not judge
whether they agree. Decide whether each engine computes the objects the contract and the notes ask for, at the
generality they ask for, and prints them faithfully. The run record is the builder's account; check it against the
artifacts, and do not take it as evidence.

## What to check
0. **The repair.**
   - **Wolfram item 4.** For each object Y1 names, other than the six of `4n`: does Wolfram print the conditions,
     computed exactly from its printed item-3 closed form, under which the object vanishes for every `D` on the domain,
     complete, with every special exponent value treated, as equations and inequalities in the named quantities, with
     `Q` on the family under each condition, one per condition, computed on that condition? Or does it carry a label of
     a kind the notes allow? Check a representative set against your own computation, by an exact route and an
     independent second route. Is the supplied series of `4l` used as an identity, never truncated in place of the
     object?
   - **Labels.** Does every label meet its rule? A 6-hour label: was the attempt still running 6 hours after it
     started? An out-of-memory label: did the attempt end for lack of memory in a run whose `--memory-gib` was the
     largest the guard admits when no other guarded job is running, with no other guarded job running at any time
     during the attempt? Check the guard records, `resource-samples.jsonl`, the registry and `journalctl --user`, not
     the run record. Is any label set on an object whose attempt could have returned its conditions? Is any return that
     gives conditions filed as not obtained? Does each label say what was attempted, beside what was printed?
   - **The six, and SymPy's item 4.** Do the six carry the third-kind label in Wolfram and nothing else of item 4?
     Does every SymPy item-4 object carry the fourth-kind label and no conditions? Does any printed line or run-record
     sentence present a label as a computed result, or the stopped runs' item-4 output as a result of this build?
   - **Restatement, claims, construction, dependence.** As in round 3: no printed condition is the object, its
     coefficients, its derivatives or its jet set to zero and restated; no status line or run-record sentence states a
     property no computation produced; each engine reaches its conditions by its own construction; on copies, show by
     computation that Wolfram's item-4 conditions depend on its item-3 closed form.
   - **The final transcripts.** How was each produced: which runs, from which source (sha256), and how the parts of
     Wolfram's final transcript were joined when its final run was continued and objects were retried? Does every
     object in each transcript come from the final code? Does each transcript end complete?
   - **Anything else.** Did the repair change anything else that changes what is computed or what may be claimed?
     Compare every item-1, item-2, item-3, item-5 and item-6 object, and in Wolfram every item-4 object where `d` is
     absent, with the round-3 transcripts, by value where the bytes differ (SymPy's printed form and Wolfram's
     generated names vary between runs).
   - **The run record.** Does it carry what `4n` and `4o` ask: the history of every item-4 run on each of the six, the
     stopped SymPy runs as history, and, for each engine, every labelled object and the kind of its label?

   Items 1, 2, 3, 5 and 6 below were reviewed in rounds 0 to 3: for them, check only whether the repair changed any
   object, and review any change against them. Items 4 and 7–12 apply in full to the code the repair changed and to the
   objects it prints.
1. **The objects `Q` (Part B item 1).**
   - Does each engine construct every object Part A's items 1 and 2 name, at each surface and for each Cartesian
     component, with the grouping by label, the whole pressure, the expansion with its order and existence condition,
     and the whole minus the expansion?
   - Are the three named sums built as Part B defines them, at each surface and over both surfaces: the rest part
     (carried, whole pressure, viscous, other), the transport (centre and radius), and the crossing (rest plus
     transport)?
   - Are these the same objects the accepted S15 step printed? The accepted engines and transcripts are listed below.
2. **The separation derivative and the two parts (item 2).**
   - Is `𝒟Q` the derivative Part B defines: drain 1's centre `X_1` moved along `D = X_1 − X_2`, and every changing
     quantity `F` moved by its change `F′`, at one instant?
   - Is every quantity Part B lists as changing given a change that is a general live object with the quantity's own
     arguments, including `Ẋ_1`? Is every quantity Part B lists as held (`X_2`, `Ẋ_2`, `r_in`, `R` and their rates, and
     the form of `p_br`) held, and kept in the objects rather than set to zero?
   - Is `E_Q^geo = (𝒟Q)|_{F′=0} + 2Q` and `E_Q^resp = E_Q − E_Q^geo`, as objects, whatever form a term is written in?
   - Does each departure object carry a condition under which its separation derivative exists, labelled as such where
     it is printed? Is any condition presented as one without which the derivative fails to exist, or as complete,
     without that being shown? Is any object printed where its derivative need not exist on the condition it carries?
3. **The profile family (item 3).**
   - Are the profiles `u_k = A_k(t) r^{a_k}` and `ρ_k = B_k(t) r^{b_k}`, with the amplitudes live functions of `t` and
     the exponents symbols, substituted in item 3 and item 4 only?
   - Where the integrand holds no general field besides the family's amplitudes, the region's quantities and the
     velocities, is the integral evaluated in closed form, exactly, with its dependence on `d` explicit and every
     special exponent value treated? Where it holds a general field, is the object printed with the family
     substituted and left unevaluated?
   - The run record lists objects an engine labels as not obtained in closed form. For each, is the label honest: does
     that engine's route in fact fail on it? Is the object reachable in closed form by an exact route? Show by
     computation either way. Is each object an engine prints as a closed form exact on the domain?
4. **The conditions (item 4).** For each item-3 object in closed form:
   - Is each printed condition a condition on the exponents and, where they enter, on the amplitudes, the radii and
     their rates, and the velocities, under which that object vanishes for every `D` on the domain? Is it exact, with
     every special exponent value treated, printed as equations and inequalities in the named quantities?
   - Is it reached by computation, or is it the object set to zero and restated?
   - Is the set of conditions complete for that object? Is `Q` on the family printed under each one?
5. **The velocity parts (item 5).** Is `Q − Q|_{Ẋ_k=0}` printed for each drain `k`, with `Q` itself unchanged?
6. **Dimensional homogeneity (item 6).** Is it computed for items 2 and 5 from a stated dimension assignment, with each
   change carrying the dimension of its quantity, or is it typed?
7. **Domain and exactness.** Is each object stated on `|X_1 − X_2| > R > r_in > 0`, with both cores excluded, under
   S15's pressure conditions and, for each departure object, its own condition? Is any limit or expansion taken in
   place of an exact object? Where an engine does not obtain an object exactly, does it print a label saying so beside
   what it did print?
8. **The supplied inputs.** Are Part A's balance, its fields, the split and the surfaces written as Part A states them,
   and the family and the power-law profiles as Part B defines them? These, and each field's argument list, are the
   only places the physical symbols may be combined by hand. Is any other expression involving them typed rather than
   reached by computation? Is a force law, a coefficient, a normalization or a mass introduced, or a printed object
   identified with a pull or with the force on a throat? Is anything that varies held fixed, set to zero, or given a
   profile outside the places Part B defines (CLAUDE.md M3)?
9. **Emission.**
   - For each printed object, which line computed it?
   - Are operands and residuals printed before any guard?
   - Does a script state a conclusion, or emit a `PASS`/`FAIL`/`AGREE`/`MATCH`/`VERDICT`/`STATUS` token or a prose
     result?
   - Does whether an object is printed depend on its value?
   - Does a tag name state a value, a sign or the shape of an answer?
   - Where a residual is printed, were its two operands produced by independent routes, or is it zero by
     construction? Where an object has no second route, is that said rather than dressed as a check?
10. **Independence of the engines.** Does the Wolfram script import or read anything from the SymPy side? Does it
    follow the SymPy script line by line, as a transcription rather than an independent construction?
11. **Reproducibility.** Does the run record name, for each engine, the command, the guard outcome, the transcript and
    its sha256? Does a guarded rerun of each final script reproduce its transcript?
12. **Controls.** Ablate each engine on copies, and report whether the printed output moves for each, quoting the
    literal diff:
    - at least one FORM ablation of a load-bearing object in each engine. Change its structure, not a coefficient.
      Examples: hold a changing quantity fixed along the family, move `X_2` instead of `X_1`, scale `R` with the
      separation, keep a change term inside the geometric part, drop a special exponent branch, or replace `p_br` by
      a linear law;
    - one-sided corruption of the supplied balance or of the family in one engine's copy.


## What you are handed
- The artifacts above, and the round-3 version at `90441b89` with its transcripts.
- The round 0 to round 3 dispositions and lookups:
  `research/pde_ledger_v3/directives/_measurements/S15_reading_build_review_r{0,1,2,3}_{disposition,lookups}.md`.
- The repair notes listed above, with their fact-check files in `research/pde_ledger_v3/directives/_measurements/`.
- The two decision lists: `research/pde_ledger_v3/directives/S15_reading_decision_list.md` (cleared at `7f52b4d2`),
  with its review dispositions beside it in `research/pde_ledger_v3/directives/_measurements/`
  (`S15_reading_decision_list_review_p{1,2,3}_disposition.md`), and
  `research/pde_ledger_v3/directives/S15_decision_list.md` (cleared at `91248b67`), with its dispositions in the same
  directory (`S15_decision_list_review_p{1,2,3}_disposition.md`).
- The accepted S15 step: `research/pde_ledger_v3/scripts/S15_two_drain_region.py`,
  `research/pde_ledger_v3/mathematica/S15_two_drain_region.wl`, and their transcripts
  `research/pde_ledger_v3/scripts/out/S15_two_drain_region.out` and
  `research/pde_ledger_v3/mathematica/out/S15_two_drain_region.out`.
- The supplied premises' sources: `research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md` lines 340–378 (P5 and P6),
  `research/pde_ledger_v3/steps/S9b_light_bending.md` lines 667–668, and
  `research/pde_ledger_v3/directives/S16_decision_list.md` lines 62–64, 76–86 and 89–90.
- The guard `scripts/s11c_guarded_run.py`, its registry `~/.local/state/s11c-guard/resource-pool/reservations.json`,
  and the user journal.
- `CLAUDE.md`, sections M1–M3, E1 and E2.

## Required method
- Your scratch directory is `/var/projects/_review_scratch/s15-reading-r4-<your leg name, claude or grok>/`. Create
  it, work only there, and change into it before you run any script, so that every file a script writes lands there.
  Write nothing in the repository tree, and do not work under `/tmp`.
- **Derive first.** Before you open either engine's source or the final transcripts, write your own CAS script that
  derives, from the supplied balance and definitions with every field and every change live, the item-3 closed forms
  of a representative set of the objects Y1 names, with their special exponent values, and their item-4 conditions,
  with `Q` under each. Include at least one carried group, one centre-transport group and one radius-transport group
  that contain drain 1, at a surface, and the region momentum or rate of one single-drain group. Where you do not
  obtain an object or a condition exactly, say what you attempted. Save the script and its literal stdout in your
  scratch directory. A claim without its script and literal stdout will be discarded; a prose derivation is worth
  nothing.
- Save every recomputation and ablation script, and its literal stdout, in your scratch directory, and report the
  paths.
- **Quote both sides** for every finding: the contract's or note's text and the engine's line.

## Running CAS
- Copy the artifacts to a new directory in your scratch directory, and ablate the copies. Never modify the working
  tree.
- Run every CAS job through the memory guard, in pooled mode:
  `python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s15-reading-review-r4 --memory-gib <GiB> --log-directory /var/projects/toy_physics/_scratch/s11c/s15-reading-review-r4-guard/<a new directory name per run> -- <command>`.
  Size `--memory-gib` to what the job needs. The final runs' memory use is in their guard records
  (`resource-samples.jsonl`). Never wrap the guard command in a shell `timeout`.
- Run Wolfram scripts as `math -script <file>`. The guard caps each job at 32 tasks, and a Wolfram kernel that starts
  more threads than that deadlocks with every thread idle. In any Wolfram script you write or run, set the kernel's
  thread limits first, as `research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/S9b_flow_oracle_check.wl` does
  (lines 9–10). Run one kernel at a time: before each Wolfram run, check `pgrep -x WolframKernel`, and if a kernel is
  already running, wait until it exits.
- The SymPy final run took about 30 minutes under the guard. Wolfram's final item-4 attempts took from seconds to
  several hours each, and some ended for lack of memory (their guard records and the runner's `FINISHED_KERNEL`
  lines). A full Wolfram rerun takes many hours; rerun the objects your findings need.
- If the guard refuses admission, for any reason, wait a few minutes and retry. Keep waiting until it admits you; do not
  finish without running the jobs your findings need. If the guard kills a run, report it.
- If one of your own guarded runs stops printing for a long stretch and you judge it stalled, you may stop it: run
  `systemctl --user stop <its unit>` while its guard launcher is still alive, then confirm that the unit is inactive and
  that its entry is gone from the registry. Never edit the registry, and never stop a run that is not yours. Report
  each run you stop.
- A guarded run must not outlive you. Before you exit, stop any of your guarded runs still active, as above, and confirm
  that the registry holds none of yours.
- Nobody will answer a question while you run. To wait for a guarded job, poll it (its process and its `outcome.json`)
  with a blocking shell command until it ends. Do not rely on queued notices or detached watchers. A turn that ends
  without a tool call ends the run.
- There is no time limit.
- Do not read other review outputs.

## Filter
Report a finding only if it changes what is computed or what may be claimed from the engines' output. Do not report
machinery that changes neither, and do not report that a script would be wrong on a different input.

## Output
- **Verdict:** CLEAR or NEEDS REVISION.
- **Findings,** numbered. For each: the location, what is wrong, the evidence (script path and literal stdout), and
  what must be true instead.
- A short list of what you checked and found sound.

Write the report to `REVIEW.md` in your scratch directory. Keep it under 1,500 words, not counting command output.

## Bounds
Write your report and exit. Do not spawn agents, and do not build `run_all`, `watch` or supervisor orchestration.

Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
