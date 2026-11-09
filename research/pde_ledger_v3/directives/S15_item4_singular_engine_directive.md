# S15 reading build, item 4: a Singular engine (builder directive)

Directive from the orchestrator. You are a new builder on this repository. `AGENTS.md` and `CLAUDE.md` apply. Paths
are relative to the repository root, `/var/projects/toy_physics`, unless they start with `/`. The lookups behind this
directive are in `research/pde_ledger_v3/directives/_measurements/S15_item4_singular_engine_directive_lookups.md`.

## Why

- Item 4 of `research/pde_ledger_v3/directives/S15_reading_decision_list.md` (its section "What the build computes")
  asks, for each object `Q` whose item-3 object is in closed form, for the conditions under which that object vanishes
  for every `D` on the domain, exact and complete, and for `Q` on the family under each condition.
- In the reading build as accepted at `cc62dc9f`, item 4 has one engine, the Wolfram engine. By a user decision
  recorded in `research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md`, the SymPy engine attempts no item-4
  object.
- **User decision (2026-11-09).** The orchestrator recommended that Singular become item 4's second engine before S15
  closes, so that item 4 has two engines that can disagree. The user replied: "Go with your recommendation. Sounds like
  Singular is a solid choice for this".
- This build adds that engine. A later, separate comparator compares its output with the Wolfram engine's. Do not
  compare the two engines here.

## Supplied (unfalsifiable in this build)

- Everything the reading decision list supplies: its sections "Supplied (unfalsifiable in this build)", "The
  separation and the family (supplied definitions)" and "The profile family", with the sections of
  `research/pde_ledger_v3/directives/S15_decision_list.md` they name.
- The section "Supplied (unfalsifiable in this trial)" of `research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md`:
  the distance from drain 1 over the region, and its series in `z = r/d`. The series is an identity on the whole region
  and may be used as one. A finite number of its terms is not the object.

## The objects

The objects `Q` are those the reading decision list's "What the build computes" names, at each surface and over the
region, for each Cartesian component and each group by label. This build takes each `Q` whose item-3 integrand on the
profile family holds no general field besides the family's amplitudes, the region's quantities and the velocities. For
each such `Q`, the engine prints:

1. its item-3 object, `E_Q^geo` on the profile family, exact, with its dependence on `d` explicit;
2. its item-4 objects: the conditions under which that object vanishes for every `D` on the domain, exact and complete,
   with every special exponent value treated, printed as equations and inequalities in the named quantities; and `Q`
   on the family under each condition.

Each object's name says its surface or the region, its piece or group, its component and its labels, so that a reader
can match it with the reading decision list's objects.

## What must be true

- The engine constructs every object from the supplied equations and definitions. It imports nothing from the SymPy or
  Wolfram engines of S15 or of the reading build, reads none of their files or outputs, and is not a transcription of
  either.
- The steps that decide when an object vanishes for every `D` run in Singular: splitting into the cases that question
  needs, special exponent values among them, and finding the conditions in each case. Other code may form the objects,
  split them into the cases that forming them needs, write them as Singular's input and read Singular's output. That
  code finds no condition and makes no other split, and every split it makes is printed. Python 3 with SymPy is
  available for it. Singular is `/var/projects/_tools/singular/bin/Singular`, version 4.4.1, with its libraries in
  `/var/projects/_tools/singular/share/singular/LIB/`.
- The construction is yours.
- Every object is exact on the domain under the conditions it carries. No limit or expansion is taken in place of an
  exact object.
- Any mathematical fact the construction relies on and does not compute is printed as a stated premise, beside each
  step that relies on it.
- The script prints computed objects and never states conclusions. Before each step it prints a line that starts with
  `WORKING_ON:` and names the step. A name may name an object, never its value or a comparison outcome. Whether an
  object is printed does not depend on its value.
- Where the engine does not obtain an object's conditions exactly and completely, it prints a label saying so, beside
  what it did print, with what was attempted. A branch that the construction sets aside without determining its
  conditions is printed as not obtained. An object with any branch printed as not obtained has not had its conditions
  returned.
- An attempt ends when it returns its conditions, when it finishes with a branch printed as not obtained, when it runs
  out of memory at its cap, or on an error in your code. Nothing stops it for time or for a period without output. If
  an attempt ends on an error in your code, you may fix the code and run again. The run record reports every attempt.
- The build computes no force law, coefficient, normalization or mass, and identifies no printed object with a pull or
  with the force on a throat. No acceptance criterion references an expected value.

## Where things go

- The engine's source goes in `research/pde_ledger_v3/singular/` (new). Its final transcript goes to
  `research/pde_ledger_v3/singular/out/S15_separation_departure_item4.out`. Transcripts there are not committed
  (`.gitignore`).
- The run record goes to `research/pde_ledger_v3/directives/_measurements/S15_separation_departure_singular_run_record.md`,
  and the run manifest to `research/pde_ledger_v3/directives/_measurements/S15_separation_departure_singular_run_manifest.json`.
  The manifest holds what `CLAUDE.md`'s "Repository infrastructure" section lists for a run manifest: the command; the
  sha256 of every script and input the run read; the start and end time, the exit code and the peak memory; the
  transcript's path, size and sha256; and, verbatim, every transcript line the run record relies on.
- No file outside `_scratch/` cites a path under `_scratch/`. Where the run record relies on a guard record, the
  manifest quotes what it relies on.
- `README.md` gains a section, beside "Mathematica Notes", that says how Singular is installed here and how to run this
  engine.

## How to run it

- Write only in `research/pde_ledger_v3/singular/`, the two `_measurements/` files above and `README.md`, and, for your
  working files and the guard's run logs, in new directories under `_scratch/s15/singular_engine/` and
  `_scratch/s11c/s15_singular/` (the guard takes its log directory only under `_scratch/s11c/`). Change nothing else
  in the repository or anywhere else, and run no `git` command that writes. The guard writes its own registry
  (`~/.local/state/s11c-guard/resource-pool/reservations.json`) when it admits and ends a job; change nothing in it
  yourself.
- Run every engine attempt through the guard, `scripts/s11c_guarded_run.py`, with `--pool s15-singular`. When no other
  guarded job is running, use the largest `--memory-gib` the guard admits. If the guard refuses admission, wait a few
  minutes and try again. Poll each job until it ends.

## The report, and where to stop

When the final run has ended, write a short final report, and stop. Give:
- your construction: which steps run in Singular and which run elsewhere;
- for each attempt: the run, its start and end, its `--memory-gib`, its outcome, its wall time and its peak memory;
- every object that carries a not-obtained label, and what was attempted for it;
- every file you wrote outside `_scratch/`, with its sha256.

This task ends at your report. Commit nothing, and start no other agent or review.

Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
