# Independent review: run-manifest rule

Leg: `grok`. Read-only on `/var/projects/toy_physics`. Commands and outputs are in `measurements/` beside this file.

**Verdict: NOT CLEAR.** One finding.

## Artifact

Working tree of `CLAUDE.md`, `AGENTS.md`, `docs/development_pipeline.md`, and `.gitignore` against `eceafc57233261403a75fc44f70c4e02549bb4ec`. The other uncommitted paths (S15 scripts, S15 transcripts, `S15_separation_departure_run.py`, `prep-for-compact/SKILL.md`) were not reviewed.

Command: `sha256sum` of the four paths, and `git diff eceafc57 --stat` of the same four.

```
7f0c71140a0046c0f79c382379de61c8af47152fd9302c982f4486a3c9fbd3a5  CLAUDE.md
ce029cc5249f06bc193b55565de7618d77ebdf4630a74ef98771e71cc7ac8f9b  AGENTS.md
471a56c8b32cf5b21cf972025783bcb2515398cfaa8d5756ac1d23f90d8899b5  docs/development_pipeline.md
e15360bdbe0451ab3d12f83bc46605903f4534cf3286de9330700a6bf31a9185  .gitignore
```

```
 .gitignore                   | 12 +++++++--
 AGENTS.md                    |  2 +-
 CLAUDE.md                    | 63 ++++++++++++++++++++++++++++----------------
 docs/development_pipeline.md |  5 ++--
 4 files changed, 54 insertions(+), 28 deletions(-)
```

Full diff: `measurements/01_diff_four_files.txt`.

## Finding — Staging ban covers every `.out`, not the two engine directories

`CLAUDE.md:284` (sentence continues at line 285):

```
- ⛔ Never stage a `.out`, and ⛔ never `git add -f` one. When a run replaces a transcript that is still tracked,
  the commit that adds the new manifest removes the old file from the index (`git rm --cached`); history keeps it.
```

The user's decision, and the rest of this section, cover new engine transcripts at `research/pde_ledger_v3/scripts/out/*.out` and `research/pde_ledger_v3/mathematica/out/*.out`. Line 284 forbids staging any file named `.out`. `CLAUDE.md` is the file the orchestrator follows when committing. `AGENTS.md:31` states the narrower rule (`research/pde_ledger_v3/*/out/*.out`). The governing sentence is the wider one, so the wider one is what a commit follows.

The same diff warns that a blanket `.out` rule drops committed evidence. `.gitignore:19-27`:

```
# ⛔ NOT a blanket `*.out`. A blanket rule silently dropped Mathematica arbiter/audit transcripts, which
# were committed dual-engine evidence (Phase 4(a)). LaTeX hyperref aux is ignored only where LaTeX runs, and
# the v3 engine transcripts only in their two out/ directories (below). ...
# 2026-09-01..11-02 stay tracked (git-annex + GIN); this rule only keeps new ones untracked.
```

The ignore rules do not hide those other files. Command: `git check-ignore -v --no-index` on `research/pde_ledger_v2/mathematica/out/ledger_stage001_solid_angle_second_moment_primitives_mathematica_audit.out`, `research/pde_ledger_v3/_measurements/S11c_c2_adjudication_verify_EG.out`, `research/pde_ledger_v3/directives/_legs/S11c_c1_grok_leg_derivations.out`, a `research/pde_ledger/scripts/output/.../*.out`, and `research/4d/scripts/euler.out`. Output is empty and the exit status is 1 (`measurements/19_ignore_scope_proof.txt`). They are not ignored. Line 284 still forbids staging them.

Command: `git ls-files '*.out'` grouped by parent, excluding paths directly under the two v3 engine `out/` directories (`measurements/19_ignore_scope_proof.txt` and `measurements/21_out_split.txt`):

```
      8 research/4d/scripts
     18 research/pde_ledger/scripts
     44 research/pde_ledger_v2/mathematica
      3 research/pde_ledger_v2/notes
      1 research/pde_ledger_v2/wlab3b.out
      1 research/pde_ledger_v2/wlab_arity.out
     16 research/pde_ledger_v3/directives
      4 research/pde_ledger_v3/_measurements
```

The 16 under `directives` are 14 `directives/_legs/*.out` and 2 `directives/_measurements/**/*.out`. The 4 under `research/pde_ledger_v3/_measurements/` are the E1 stdout files `S11c_c2_adjudication_verify_EG.out`, `S11c_c2_adjudication_verify_F.out`, `S11c_c2_export_repair_adjudication_verify.out`, and `S11c_c2_export_repair_rereview_claude_recompute_artifact.out`.

What must be true: a `.out` that is not a new transcript in those two `out/` directories, or a preserved copy under their `preserved/` directories, can be staged. That includes a v2 audit transcript, a `directives/_legs/*.out`, and a `_measurements/*.out` that holds literal stdout. E1 (`CLAUDE.md:109-112`) puts a claim's literal output in a `_measurements/` file. Four of those files are named `.out`. Under line 284 they cannot be staged, so a later update of that evidence does not enter the tree.

The following sentence, `git rm --cached` when a run replaces a transcript that is still tracked, uses the section's own meaning of transcript (the two globs at `CLAUDE.md:264-265`). The defect is the unscoped staging ban, not that sentence.

## Fidelity of the decision otherwise

The decision was: for each new engine run, commit a manifest (command, script hash, runtime, transcript size and hash, lines the step record quotes); do not commit the transcript; keep it local while the step is open; delete it when the step closes; regenerate by rerun. Legacy transcripts already committed stay annex content.

`CLAUDE.md:269-283` and `docs/development_pipeline.md:196-199` do that. The manifest also requires the sha256 of every input the run read, the exit code, the recorded peak memory, the transcript path, and lines a review disposition or a comparator claim relies on (`CLAUDE.md:271-274`). Those are additional fields inside the manifest. They do not commit a transcript. Tracked transcripts from 2026-09-01 through 2026-11-02 stay in git-annex, and a fresh clone still needs `datalad get` (`CLAUDE.md:286-290`). A replaced tracked transcript is removed from the index with history kept (`CLAUDE.md:284-285`).

## Other controls

E1. For an engine run, the committed literal output is the manifest's verbatim lines (`CLAUDE.md:22-23`, `CLAUDE.md:109-112`, `CLAUDE.md:274-276`). That holds for a manifest that is not itself blocked by the finding above. A committed claim quotes those lines.

G3. Unchanged. No commit before both reports (`CLAUDE.md:164-169`).

G4. The reviewed baseline commit for a transcript is its manifest, plus a local copy at `out/preserved/<name>.<round>.out` until the step closes (`CLAUDE.md:182-183`, `CLAUDE.md:277-280`). Command: `git check-ignore -v` matches both preserved paths to `.gitignore:30` and `:31` (`measurements/19_ignore_scope_proof.txt`). The copy is on disk and is not a committed blob. Round-to-round comparison while the step is open reads that file. After close, the copies are deleted and the committed manifest remains. What must be true, and what the text already requires: every line a finding relies on is in the manifest that was committed before the overwrite.

M1 and the comparator. Comparators read the engine `.out` streams, not the manifest. `research/pde_ledger_v3/scripts/S10_cross_engine_comparator.py:82-83` rejects an input whose suffix is not `.out`; lines 127 and 939 require exactly two `.out` files. `S11c_a_cross_engine_comparator.py:2042-2043` describes those paths as the committed PY and WL `.out` tag streams. `S15_cross_engine_comparator.py:1441` requires `--py` and `--wl`. While the step is open, `CLAUDE.md:277-278` leaves those files at the `out/` paths the comparator opens. On a fresh clone, a legacy tracked transcript is an annex symlink until `datalad get` (`CLAUDE.md:289-290`). A transcript created under this rule is absent until the engines are rerun; the Wolfram half needs a Mathematica licence (`CLAUDE.md:282-283`). Later steps import `scripts/*_exports.py`, which this rule still leaves as plain git (`CLAUDE.md:293`). Command: `git check-ignore -v --no-index -- research/pde_ledger_v3/scripts/S10_exports.py` exits 1 (`measurements/19_ignore_scope_proof.txt`). `git ls-files` shows the export modules directly under `research/pde_ledger_v3/scripts/`, not under `out/` (`measurements/14_ls_files_out_exports.txt`).

E2. `docs/development_pipeline.md:196-201` still commits the ablation harness, its invocation, its digests, and the literal `(baseline, corrupted, diff)` transcript through `_measurements/`. The parenthetical there says the full engine transcript is not committed. The ablation record holds if that `_measurements/` file can be staged. A file in that directory named `.out` is the finding above.

Observable gate record. `CLAUDE.md:63-68` is unchanged. A stdout path in that record can be the manifest, which is the file the tree holds.

Consumer check. `docs/development_pipeline.md:329-333` says a named comparator is re-run, and line 342 says a structured artifact is verified by running its consumer. This amendment does not edit those lines. The consumer can run while the `.out` streams are on disk. `CLAUDE.md:281-282` checks a regeneration against quoted lines and says a byte-hash mismatch alone is not a finding. That matches the existing drift guard at `docs/development_pipeline.md:192-195` (tagged payloads, not full stdout). What must be true: lines the comparator claim relies on are in the manifest, and the comparator itself is run on the `.out` streams when that consumer check is required. The manifest is not an input the comparator accepts.

## Operability of `.gitignore`

Patterns at `.gitignore:28-31` match new transcripts and preserved copies, and they do not exclude tracked files from the index.

Command: `git check-ignore -v` (default) on a new `scripts/out/*.out`, a new `mathematica/out/*.out`, and both `preserved/...r1.out` paths (`measurements/19_ignore_scope_proof.txt`):

```
.gitignore:28:research/pde_ledger_v3/scripts/out/*.out	research/pde_ledger_v3/scripts/out/NEW_untracked_example.out
.gitignore:29:research/pde_ledger_v3/mathematica/out/*.out	research/pde_ledger_v3/mathematica/out/NEW_untracked_example.out
.gitignore:30:research/pde_ledger_v3/scripts/out/preserved/	research/pde_ledger_v3/scripts/out/preserved/S15_separation_departure.r1.out
.gitignore:31:research/pde_ledger_v3/mathematica/out/preserved/	research/pde_ledger_v3/mathematica/out/preserved/S15_separation_departure.r1.out
```

Command: `git check-ignore -v` on the two tracked `S15_separation_departure.out` paths exits 1 with no output. `git ls-files` still lists 34 `scripts/out/*.out` and 16 `mathematica/out/*.out`. `git status --short` still shows `T research/pde_ledger_v3/mathematica/out/S15_separation_departure.out`. With `--no-index`, those tracked engine paths match `.gitignore:28` and `:29` (`measurements/16_check_ignore.txt`). The pattern applies to new files. It does not drop the tracked ones.

Not ignored (same `check-ignore --no-index` run, exit 1, empty): `S10_exports.py`, `scripts/out/O2_cross_engine_comparator_catalog.json`, `scripts/out/O2_cross_engine_comparator_accounting.jsonl`, `scripts/out/not-an-out.txt`, `lean/out/example.out`, `scripts/out/subdir/nested.out`. A manifest path `research/pde_ledger_v3/_measurements/S15_run_manifest.md` is not ignored (`measurements/20_manifest_ignore_and_leftover_scan.txt`, exit 1). `research/pde_ledger_v3/_measurements/.gitignore` only re-includes two `*.tar.gz` snapshots. Root `.gitattributes` still annexes only the two `out/*.out` globs, so a `_measurements/` file and an `*_exports.py` stay plain git.

`find research/pde_ledger_v3 -type d -name out` returns only `mathematica/out` and `scripts/out`.

## Leftovers

No remaining text in `CLAUDE.md`, `AGENTS.md`, `docs/development_pipeline.md`, `.claude/skills/*/SKILL.md`, or `STATUS.md` tells someone to commit or annex a new engine transcript.

Command: `rg -n 'datalad save|annexed automatically|Generating/updating|Commit/annex' ` on those trees. No matches (`measurements/20_manifest_ignore_and_leftover_scan.txt`). The removed `AGENTS.md` sentence "Use `datalad save`. The large `out/*.out` files are annexed automatically." is in the diff (`measurements/01_diff_four_files.txt`).

These sentences are about Codex session logs, and they say to keep them out of the tree:

- `.claude/skills/build/SKILL.md:35-39` redirects `codex exec` to a log outside the repo. Line 160: "transcripts still go outside the repo only to keep the tree clean."
- `docs/development_pipeline.md:380-383`: "Raw transcripts go outside the repository as tree hygiene."
- `.claude/skills/review-legs/SKILL.md:171`: "Raw transcripts still go outside the repository."

`STATUS.md` has no `datalad`, `annex`, or engine-transcript hit (`measurements/08_rg_leftovers_scoped.txt`).

## Verdict

NOT CLEAR.

1. `CLAUDE.md:284` forbids staging every `.out`. The decision, `AGENTS.md:31`, and `.gitignore:28-31` cover new v3 engine transcripts only. Ninety-five other tracked `.out` files, including v2 audit transcripts and `research/pde_ledger_v3/_measurements/*.out`, stay addable under gitignore and become unstageable under the governing sentence.
