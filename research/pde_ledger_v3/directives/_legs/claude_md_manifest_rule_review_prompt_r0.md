# Independent review: the run-manifest rule (engine transcripts no longer committed)

## Artifact
An uncommitted amendment, written by the orchestrator, to four files in `/var/projects/toy_physics`: `CLAUDE.md`, `AGENTS.md`, `docs/development_pipeline.md` and `.gitignore`. See it with

    git -C /var/projects/toy_physics diff eceafc57 -- CLAUDE.md AGENTS.md docs/development_pipeline.md .gitignore

Its sha256 baseline:

    7f0c71140a0046c0f79c382379de61c8af47152fd9302c982f4486a3c9fbd3a5  CLAUDE.md
    ce029cc5249f06bc193b55565de7618d77ebdf4630a74ef98771e71cc7ac8f9b  AGENTS.md
    471a56c8b32cf5b21cf972025783bcb2515398cfaa8d5756ac1d23f90d8899b5  docs/development_pipeline.md
    e15360bdbe0451ab3d12f83bc46605903f4534cf3286de9330700a6bf31a9185  .gitignore

The working tree holds other uncommitted changes that belong to a build in progress: `research/pde_ledger_v3/scripts/S15_separation_departure.py`, `research/pde_ledger_v3/mathematica/S15_separation_departure.wl`, the S15 `.out` files, `research/pde_ledger_v3/mathematica/S15_separation_departure_run.py` and `.claude/skills/prep-for-compact/SKILL.md`. They are not part of this artifact. Do not modify them, and do not run anything that reads them heavily; a Wolfram job is running.

`CLAUDE.md` is the project's governing rule file. `AGENTS.md` instructs Codex builders.

## The decision it implements
The user asked on 2026-11-02: "Also, how much have you put in git annex? We should be smart and not commit stuff users can regenerate." The orchestrator reported the measured annex contents and asked two questions. The user answered:

- "Going forward, what should we commit for each engine run's transcript?" The user chose **"Manifest only (Recommended)"**, whose option text read: "Commit a small record per run: the command, script hash, runtime, the transcript's size and hash, and the lines the step record quotes. Full transcripts stay local while the step is open, then get deleted; rerun to regenerate. Needs an edit to CLAUDE.md's annex rule, reviewed like any rule change."
- The second question was about deleting old local-only annex versions. It has already been carried out and is not part of this artifact.

The same day the user also stated: "Nothing should be committed in _scratch. Nothing should reference _scratch. That's literally why it exists. So if there's something you need from it, move it out."

## What to check
Report a finding only if it changes what gets committed, what evidence exists, or what may be claimed. Wording preferences are not findings.

1. **Fidelity.** Does the amendment carry out the decision above, no more and no less? Check what the manifest holds, that new transcripts are not committed, that transcripts stay local while a step is open and are deleted when it closes, and how the transcripts committed earlier are treated.
2. **Other controls.** Does the change weaken any other control in `CLAUDE.md`? In particular:
   - E1: a claim carries its command and literal output;
   - G3 and G4: the reviewed baseline is preserved before a repair, so a finding can be reproduced;
   - M1: the cross-engine comparison, and the comparator's inputs;
   - E2: ablation evidence;
   - the observable gate record;
   - the rule that a regression check must run real consumers.

   Where a control loses evidence it needs, say what must be true for it to hold. Do not propose wording.
3. **Operability.** Can each consumer still get its inputs: a comparator, review legs comparing one round with the next, a record that quotes lines, a later step importing a step's exports, and a fresh clone? Are the `.gitignore` rules right? They must not hide files that are already tracked, and they must not repeat the failure the `.gitignore` comment records, where a blanket `*.out` rule hid committed evidence. Test them with `git check-ignore -v` and `git ls-files`.
4. **Leftovers.** List any remaining text that still tells someone to commit or annex a new engine transcript, with file and line. Check `CLAUDE.md`, `AGENTS.md`, `docs/development_pipeline.md`, `.claude/skills/*/SKILL.md` and `STATUS.md`.

## Method
This is a document review. Ground every claim in a command and its literal output (`git diff`, `git grep`, `git check-ignore`, `git ls-files`, `sed -n`). Save each command and its output in your scratch directory, and cite them in your report. A claim with no command behind it will be discarded.

## Scratch and report
Your scratch directory is `/var/projects/_review_scratch/manifest-rule-r0-<leg>/`, where `<leg>` is the name you are given (`codex` or `grok`). It already exists. Work only there; never in `/tmp`. Write your report to `REVIEW.md` in that directory. Give each finding a short title, the file and line, the evidence command and its output, and what must be true. End with a verdict: CLEAR (no finding) or NOT CLEAR (findings listed).

## Bounds
Read-only on the repository: no edits, no commits, no `git add`, no `git rm`, no `datalad` or `git annex` commands that change anything. Write your report and exit. Do not spawn agents or build supervisors.

Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
