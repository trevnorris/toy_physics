# S9b repair build, review round 11 (K10 only): dispositions (orchestrator)

Round 11 repairs K10 only, in both harnesses, with a new author per engine (the user's 2026-10-18 decision, recorded in
`S9b_repair_build_r10_review_disposition.md`). The two engines' harnesses were reviewed separately: the Wolfram
builder finished about five hours before the SymPy builder's full run did, and the SymPy repair changes only that
harness's K10 attempts. Each section below has its own prompt, baseline, legs and lookups.

## How the round was run (both engines)
- **Builders.** A fresh Codex session per engine on item 15's builder model (`gpt-6-astra`, high effort), new to the
  harness. This ends the item-15 model deviation recorded since round 8. Briefs:
  `_scratch/s9b_build/repair_build/brief_{py,wl}_r11.md`.
- **Exports.** Each builder worked in a standalone checkout made from the round-10 preservation commit `26149ede`:
  `/var/projects/toy_physics-s9b-r11-py` and `/var/projects/toy_physics-s9b-r11-wl`. Each kept only its own engine's
  S9b files.
- **Export incident (orchestrator error).** My prune of the other engine's S9b files matched `*s9b*` against absolute
  paths, and each checkout's own directory name contains `s9b`, so the prune deleted every file except the kept ones.
  Both builders edited their harness and then stopped, because the guard was missing. I restored the 10,220 other
  paths without overwriting the builders' edits (list `_scratch/s9b_build/repair_build/restore_list_26149ede.txt`),
  committed the restore in each checkout (`c2b0ea1` SymPy, `36bd36c` Wolfram), and resumed both builders with
  `repair_build/resume_r11_restored.md`. No guarded run had started before the restore.

## Wolfram

### The round
- **Artifacts.** The ablation harness, its checks script, the rewritten build report, and 22 new
  `S9b_repair_wolfram_round11_*` measurement files. 30 files are in the baseline
  `_scratch/s9b_build/s9b_repair_build_review_baseline_r11wl.sha256`, all OK at review. The baseline also holds the
  engine, the SymPy engine and round-10 SymPy harness (for the cross-engine question only), the spec and the directive.
- **Contract.** The build directive at `36ccec3e` (amendments 1–12) and the round-10 disposition's "What must be true
  after the repair". The spec is at `34a994d5`.
- **Run history.**
  - The resumed builder started a detached driver (`_scratch/s9b_wl_repair11/run_final.py`, PID 195465) and a
    self-resume watcher (PID 195466), then ended its turn. I stopped the watcher, so that its resumption stayed in my
    hands, and resumed the session myself when the driver exited (`repair_build/resume_wl_r11c.md`, which bars
    detached drivers and watchers).
  - The driver's session ran checks (`s11c-guard-9bfe4f206cb8`), the full harness (`s11c-guard-d420f0fc3eba`), and
    the engine in the repository and alone (`s11c-guard-4d39d728ca3c`, `s11c-guard-489323f085ab`), all exit 0, on
    harness `4dc92515…` and checks `5688ac1c…`. The resumed builder then changed the harness and checks script
    (final `76659e38…`, `1fd77df8…`; the Claude leg found that run's K10 whole-object evaluations wrong), and reran
    checks and the full harness (`s11c-guard-8bbcce25edf3`, exit 0, 317 s). The earlier checks and full runs were
    renamed `*-prior-code`. The engine runs were not repeated; the engine is byte-identical to `86097a44`, so they
    are runs of the final engine.
  - The guard refused admission 10 times, each only because other jobs held memory: 2 in the driver's session and 8
    in the rerun. All runs were at 8 GiB, one kernel at a time.
- **Prompt and legs.** The prompt is `_scratch/s9b_build/s9b_repair_build_review_prompt_r11wl.md` (sha256
  `95bd5360…`). The legs were a fresh Claude agent (opus) and Grok (`grok-4.7`, effort high), run in parallel with one
  kernel at a time, with scratch under `/var/projects/_review_scratch/s9b-r11wl-{claude,grok}/`. Their reports are
  `_scratch/s9b_build/s9b_repair_build_review_r11wl_claude.md` and `…_r11wl_grok.txt`. Grok returned **CLEAR**;
  Claude returned **NEEDS REVISION**, with one finding about the record.
- **Lookups.** `S9b_repair_build_r11wl_review_lookups.md` beside this file (generator
  `_scratch/s9b_build/gen/s9b_build_r11wl_lookups.sh`).

### Findings

| # | Leg | Finding | Disposition |
|---|---|---|---|
| 1 | Claude 1 | **The measurements file and report leave out the driver session's superseded runs.** The measurements file says at line 3 that the packet "uses only the final round-11 sources and their current executions". Neither file mentions the driver session's checks and full runs on harness `4dc92515…`. Yet the refusal count (10) includes that session's two refusals, and the engine runs presented as final come from that session. That superseded full run printed wrong K10 whole-object evaluations while its `WholeEmittedObjectsEvaluated` flag read True; the flag tests only for unevaluated heads. Suggested fix: a paragraph in both files on that session's runs. No rerun. | **Facts confirmed; no repair.** The lookups show the two superseded runs, their hashes and the absence of any mention in either file. The record's own statements are true. Item 17 asks for every guard refusal or kill, and all 10 refusals are reported. The brief asked that nothing in the files describe a superseded run. The engine runs are runs of the byte-identical final engine. No witness, decision, mark or K11 pair rests on a superseded run; both legs confirmed that the packaged stream is the final full run. The superseded session is recorded above, with the leg's evidence (`…_r11wl_claude_evidence/out/superseded_run_evidence.stdout`). So nothing outstanding changes what the harness computes or what may be claimed from it. The flag's meaning (no unevaluated heads, not correct values) is noted; both legs evaluated the final K10 objects themselves and they match. |

**Sound, both legs, each with its own scripts and stdout.**
- **K10's mutation.** The K10 copy negates only the along-ray `Inactive[ForAll]` inside `propagationGate` (engine
  lines 262–264), in `rayGate` at line 266, before the observables are gated. This is the same condition the SymPy
  K10 copy negates (`real_domain`, SymPy engine line 696, gate copy at line 770). The ray construction and the
  traversal condition stay un-negated. `BRANCH_EXISTENCE`, `PATH_TRAVERSAL` and `BRANCH_TYPE` are byte-identical
  across the copies.
- **K10's witness.** It is at `WL_S9B_A_NONRECIPROCAL_PATH_DEPENDENCE`, domain part, with evaluations
  `{False, True}`. Each leg reproduced it with its own Wolfram script at the printed inputs: domain decisions
  corrupted False and baseline True, the two emitted objects differ, and traversal is True in both copies.
- **K10's site set.** 3,850 parts over 317 tags (both legs). The Claude leg found the site tags equal to the set of
  tags that differ between the copies.
- **Form ablations.** Claude's ablation replaced the negation with a logically neutral restructuring, and the
  witness became "undecided" (`{True, True}`). Grok's replaced it with a pointwise negation, and a witness remained;
  after path and fresh-symbol normalization only the manifest text differed outside K10.
- **Nothing else moved.** Every other knife (K1–K9b, K11), the identity marks and the baseline match round 10's
  output except for the checkout path in `"Source"` strings. The marks change only where K10's copy no longer
  differs: the three `BRANCH_EXISTENCE` tags are now marked identical. The engine is byte-identical to `86097a44`.
  The harness and checks diffs touch only K10.
- **The checks script** only prints, including its new K10 records.
- **No verdict, target or expected value** is printed by any artifact.

**Clear.** Under G4 nothing outstanding changes what the Wolfram harness computes or what may be claimed from it. The
Wolfram harness, checks script, report and round-11 measurement files are accepted as reviewed.

## SymPy

### The round
- **Artifacts.** The ablation harness, the tag inventory, the run measurements and the rewritten build report. They
  are in the baseline `_scratch/s9b_build/s9b_repair_build_review_baseline_r11py.sha256` with the engine, the export
  delta, the spec and the directive: 8 files, all OK at review in the builder's checkout
  (`/var/projects/toy_physics-s9b-r11-py`, where the legs reviewed them) and, after the copy, in this repository.
- **Contract.** As for Wolfram; the SymPy half of the round-10 disposition's "What must be true after the repair".
- **Run history.** The builder ran its probes and checks, then the full harness `s9b-round11-full-3` (exit 0,
  20,247 s, 18 GiB) and checks `s9b-round11-checks2-1` (exit 0), both started after its last edit. The guard refused
  admission 21 times, each only because other jobs held memory. The full run reached its 18 GiB `memory.max`
  (14,250 throttle events, no OOM kill) and finished.
- **Prompt and legs.** The prompt is `_scratch/s9b_build/s9b_repair_build_review_prompt_r11py.md` (sha256
  `6f45c0f9…`). The legs were a fresh Claude agent (opus) and Grok (`grok-4.7`, effort high), with scratch under
  `/var/projects/_review_scratch/s9b-r11py-{claude,grok}/`. Their reports are
  `_scratch/s9b_build/s9b_repair_build_review_r11py_claude.md` and `…_r11py_grok.txt`. Grok returned **CLEAR**;
  Claude returned **NEEDS REVISION**, with one finding about the record.
- **Lookups.** `S9b_repair_build_r11py_review_lookups.md` beside this file (generator
  `_scratch/s9b_build/gen/s9b_build_r11py_lookups.sh`).

### Findings

| # | Leg | Finding | Disposition |
|---|---|---|---|
| 1 | Claude 1 | **The measurements and the report omit three admitted runs.** `s9b-round11-probe-19` and `s9b-round11-gate1-1` are exploratory K10 runs of the builder's scratch scripts; `s9b-round11-checks1-1` is a checks run of a superseded harness (`53ebc382…`). All three exited 0. The measurements' `guard_runs` name only `full` and `checks2`, and the report says "**21** memory-contention admission refusals … no other refusal, guard kill or item 12 event." No record calls a superseded run final, and nothing printed comes from these runs. Suggested fix: list the three with command, log, source hash and exit code, labelled superseded or exploratory. No rerun. | **Facts confirmed; no repair.** The lookups show the three runs and that neither file names them. As for Wolfram, the record's own statements are true. Item 17 asks for every guard refusal or kill, and all 21 refusals are reported; none of the three runs was refused or killed. The brief asked that nothing in the report describe a superseded run, and that the measurement files come from the final runs. No witness, decision, mark, fixture or inventory row rests on these runs: both legs matched the final output's source pins to the baseline. The three runs are recorded here, with the leg's evidence (`…_r11py_claude_evidence/out/run_records.stdout`). So nothing outstanding changes what the harness or inventory computes, or what may be claimed from them. |

**Sound, both legs, each with its own scripts and stdout.**
- **K10's inputs and root sets.** The witness is at `PY_S9B_A_DEFLECTION_G000`, the `branch_complement` domain
  part, engine line 771, at `c_0=1, Z_R=Z_E=b=12` with profiles `V=0, xi=0, rho_br=1, mu_perp=1`. Every free
  input and live profile of that part is assigned; `j_n` and `Phi` do not occur in it. Each leg decided the engine's
  own ray `ConditionSet`s at those inputs with its own code: each ray set is `{(12, 12)}`, with the negative
  turning-point root kept and rejected by the base set. Both got baseline False and corrupted True, as printed.
- **K10's mutation** is unchanged: `observable_gate = real_domain` → `sp.Not(real_domain)` at engine line 770.
- **Form ablations.** Keeping only the first root of each solution set (Claude), or keeping a true candidate only
  under an added sign condition (Grok), removed the witness or failed the root-retention fixture after printing its
  value. Leaving evaluations undone (Claude) made the attempt "undecided": an unresolved evaluation furnishes no
  witness.
- **Nothing else moved.** Against round 10, 5,433 of 5,462 common lines are identical. The 16 non-K10 attempt and
  witness lines that differ are equal after sorting their point tuples (set printing order). The other differences
  are this round's new fixture fields. Every fixture result is unchanged and true. The engine and the export delta
  are byte-identical to `86097a44`.
- **The inventory** names every printed tag with the run, source hash and markers, and copies no payload.
- **No verdict, target or expected value** is printed by any artifact.

**Clear.** Under G4 nothing outstanding changes what the SymPy harness or inventory computes or what may be claimed
from them. The SymPy harness, inventory, run measurements and report are accepted as reviewed.

## Both engines
K10 now ablates the same check in both engines, the along-ray branch-existence condition, and each engine prints a
K10 witness that both of its legs reproduced. With round 11, every knife in both harnesses has a reviewed outcome,
and the S9b harnesses are accepted. Next on the S9b path is the S9b record.
