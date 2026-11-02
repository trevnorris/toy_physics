# Independent review, round 1 (scoped): the run-manifest rule

## Artifact
An uncommitted amendment, written by the orchestrator, to four files in `/var/projects/toy_physics`: `CLAUDE.md`, `AGENTS.md`, `docs/development_pipeline.md` and `.gitignore`. See it with

    git -C /var/projects/toy_physics diff eceafc57 -- CLAUDE.md AGENTS.md docs/development_pipeline.md .gitignore

Its sha256 baseline for this round:

    1c1ce0e6b2869bf48978c5b8816a4f6da74b8b98c38764ce07d158ace4605711  CLAUDE.md
    ce029cc5249f06bc193b55565de7618d77ebdf4630a74ef98771e71cc7ac8f9b  AGENTS.md
    471a56c8b32cf5b21cf972025783bcb2515398cfaa8d5756ac1d23f90d8899b5  docs/development_pipeline.md
    e15360bdbe0451ab3d12f83bc46605903f4534cf3286de9330700a6bf31a9185  .gitignore

The working tree holds other uncommitted changes, from a build in progress:
- the S15 `.py`, `.wl` and `.out` files;
- `research/pde_ledger_v3/mathematica/S15_separation_departure_run.py`;
- `.claude/skills/prep-for-compact/SKILL.md`;
- this review's own record files under `research/pde_ledger_v3/directives/_measurements/` and `research/pde_ledger_v3/directives/_legs/`, named `claude_md_manifest_rule_*`.

The S15 files and the skill file are not part of the artifact; do not modify them, and do not read them heavily, because a Wolfram job is running.

## The decision it implements
The user, 2026-11-02: "We should be smart and not commit stuff users can regenerate." Asked "Going forward, what should we commit for each engine run's transcript?", the user chose **"Manifest only (Recommended)"**: "Commit a small record per run: the command, script hash, runtime, the transcript's size and hash, and the lines the step record quotes. Full transcripts stay local while the step is open, then get deleted; rerun to regenerate. Needs an edit to CLAUDE.md's annex rule, reviewed like any rule change."

## Round 0 and the repair
Round 0 found one defect, M1. The finding, the lookups that confirmed it, what must be true, and the repair are in
`research/pde_ledger_v3/directives/_measurements/claude_md_manifest_rule_review.md` (sha256 `6184e7dffa19e24384a65a85b60b49f1b69b6ce623f29811d734fca8ec3f89e2`).
Since round 0, the only change to the four files is the repair: lines 284–288 of `CLAUDE.md`. Round 0's `CLAUDE.md` had sha256 `7f0c71140a0046c0f79c382379de61c8af47152fd9302c982f4486a3c9fbd3a5`; the other three files are unchanged.

## What to check (this round is scoped to the repair)
Report a finding only if it changes what gets committed, what evidence exists, or what may be claimed. Wording preferences are not findings.
1. Does the statement under "What must be true" for M1 now hold? Any `.out` other than a new engine transcript in the two v3 `out/` directories, or a preserved copy under them, can be staged and committed as ordinary evidence. Check the repaired text against the tracked `.out` files outside those two directories (`git ls-files '*.out'`) and against `.gitignore` (`git check-ignore -v`).
2. Does the repair keep the ban on new engine transcripts and preserved copies, and leave the `git rm --cached` sentence meaning what it meant?
3. Did the repair introduce any new defect in the lines it touched?

## Method
This is a document review. Ground every claim in a command and its literal output (`git diff`, `git grep`, `git check-ignore`, `git ls-files`, `sed -n`). Save each command and its output in your scratch directory and cite them in your report. A claim with no command behind it will be discarded.

## Scratch and report
Your scratch directory is `/var/projects/_review_scratch/manifest-rule-r1-<leg>/`, where `<leg>` is the name you are given. It already exists. Work only there; never in `/tmp`. Write your report to `REVIEW.md` in that directory. Give each finding a short title, the file and line, the evidence command and its output, and what must be true. End with a verdict: CLEAR (no finding) or NOT CLEAR (findings listed).

## Bounds
Read-only on the repository: no edits, no commits, no `git add`, no `git rm`, no `datalad` or `git annex` commands that change anything. Write your report and exit. Do not spawn agents or build supervisors.

Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
