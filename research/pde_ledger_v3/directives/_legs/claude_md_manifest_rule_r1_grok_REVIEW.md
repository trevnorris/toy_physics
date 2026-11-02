# Independent review, round 1 (scoped): the run-manifest rule

Leg: `grok`. Read-only on `/var/projects/toy_physics`. Command output is in this directory. The S15 build files and `prep-for-compact/SKILL.md` were not reviewed.

The M1 repair holds. No finding.

## Baseline

Command: `sha256sum` of the four paths and of the round-0 record (`01_baseline_and_repair_text.txt`).

```
1c1ce0e6b2869bf48978c5b8816a4f6da74b8b98c38764ce07d158ace4605711  CLAUDE.md
ce029cc5249f06bc193b55565de7618d77ebdf4630a74ef98771e71cc7ac8f9b  AGENTS.md
471a56c8b32cf5b21cf972025783bcb2515398cfaa8d5756ac1d23f90d8899b5  docs/development_pipeline.md
e15360bdbe0451ab3d12f83bc46605903f4534cf3286de9330700a6bf31a9185  .gitignore
6184e7dffa19e24384a65a85b60b49f1b69b6ce623f29811d734fca8ec3f89e2  claude_md_manifest_rule_review.md
```

These are the hashes named for this round. `AGENTS.md`, `docs/development_pipeline.md`, and `.gitignore` are the round-0 hashes. `CLAUDE.md` is the new hash.

Command: replace the current bullet at `CLAUDE.md:284-288` with the round-0 bullet quoted in the round-0 record and hash the result (`12_recon_hash.txt`).

```
current_block_found True
current_block_count 1
recon_sha256 7f0c71140a0046c0f79c382379de61c8af47152fd9302c982f4486a3c9fbd3a5
match True
current_sha256 1c1ce0e6b2869bf48978c5b8816a4f6da74b8b98c38764ce07d158ace4605711
```

`7f0c7114…` is the round-0 `CLAUDE.md` hash. The only byte change since round 0 is that one bullet.

Round-0 bullet, from `directives/_legs/claude_md_manifest_rule_r0_grok_REVIEW.md` lines 35-36:

```
- ⛔ Never stage a `.out`, and ⛔ never `git add -f` one. When a run replaces a transcript that is still tracked,
  the commit that adds the new manifest removes the old file from the index (`git rm --cached`); history keeps it.
```

Current bullet, command `nl -ba CLAUDE.md | sed -n '284,288p'` (`01_baseline_and_repair_text.txt`):

```
284  - ⛔ Never stage an engine transcript, meaning a `.out` in one of those two `out/` directories or a preserved copy
285    under them, and ⛔ never `git add -f` one. Any other `.out` (a `_measurements/` stdout file, a leg's derivation
286    output, a v2 audit transcript) is ordinary evidence and is committed as before. When a run replaces a
287    transcript that is still tracked, the commit that adds the new manifest removes the old file from the index
288    (`git rm --cached`); history keeps it.
```

"Those two" are the directories named at `CLAUDE.md:264-265` (`sed -n '262,296p'`, same file):

```
An engine run's full transcript (`research/pde_ledger_v3/scripts/out/*.out`,
`research/pde_ledger_v3/mathematica/out/*.out`)
```

## Check 1 — M1 "what must be true"

Required: any `.out` other than a new engine transcript in those two directories, or a preserved copy under them, can be staged and committed as ordinary evidence.

The repaired bullet says that. The parenthetical names a `_measurements/` stdout file, a leg's derivation output, and a v2 audit transcript.

Command: `git ls-files '*.out'`, then a path classification (`02_ls_files_out.txt`, `03_classify.txt`).

```
tracked_out_total 145
  inside-direct 50
  outside 95
tracked_preserved 0
tracked_paper_component 0
```

The 95 outside paths are the same buckets as round 0:

```
     8 research/4d/scripts
    18 research/pde_ledger/scripts
    44 research/pde_ledger_v2/mathematica
     3 research/pde_ledger_v2/notes
     1 research/pde_ledger_v2/wlab3b.out
     1 research/pde_ledger_v2/wlab_arity.out
     4 research/pde_ledger_v3/_measurements
    16 research/pde_ledger_v3/directives
```

The 16 under `directives` are 14 `directives/_legs/*.out` and 2 `directives/_measurements/reconciliation_probes/*.out` (`03_classify.txt`).

Command: `git check-ignore -v --stdin` on the 95 outside paths (`10_check_ignore_v_outside.txt`).

```
EXIT 1
```

Command: `git check-ignore -v --no-index --stdin` on the same 95 (`11_check_ignore_v_no_index_outside.txt`).

```
EXIT 1
```

Command: `git check-ignore -v --non-matching --stdin` on the same 95 (`14_outside_nonmatching.txt`).

```
outside_nonmatching_lines 95
outside_matched 0
outside_unmatched_marker 95
SAMPLE_UNMATCHED_FIRST ::	research/4d/scripts/brane_reduction.out
SAMPLE_UNMATCHED_LAST ::	research/pde_ledger_v3/directives/_measurements/reconciliation_probes/wl_provenance_probe.out
```

No ignore rule matches a tracked `.out` outside the two directories, including with `--no-index`, so a new file at one of those paths is addable with `git add`.

Command: `git check-ignore -v --non-matching --stdin` on the evidence classes the bullet names, plus a new file in each class (`13_check_ignore_v_nonmatching_probes.txt`). Unmatched lines are `::`:

```
::	research/pde_ledger_v3/_measurements/S11c_c2_adjudication_verify_EG.out
::	research/pde_ledger_v3/_measurements/NEW_stdout.out
::	research/pde_ledger_v3/directives/_legs/S11c_c1_grok_leg_derivations.out
::	research/pde_ledger_v3/directives/_legs/NEW_leg.out
::	research/pde_ledger_v3/directives/_measurements/reconciliation_probes/tower_depth_probe.out
::	research/pde_ledger_v2/mathematica/out/ledger_stage009_flat_slab_return_residual_mathematica_audit.out
::	research/pde_ledger_v2/mathematica/out/NEW_audit.out
::	research/pde_ledger_v2/wlab3b.out
::	research/pde_ledger_v2/notes/stage023_step_h_evidence/u13_u14_rescued/u13_expE.stage.out
::	research/4d/scripts/NEW.out
::	research/pde_ledger/scripts/NEW.out
```

Command: `git check-attr -a` on one file from each evidence class (`cmd_gitattributes.txt`). Each reports `annex.largefiles: nothing`. Root `.gitattributes` lines 6-8 set `annex.largefiles=nothing` for everything and `anything` only for the two `out/*.out` globs. Those evidence files stay plain git, which is how they were committed before.

The statement holds for every tracked `.out` outside the two directories and for a new file in each class the bullet names.

`git check-ignore -v` also matches one `.out` class that is neither an engine transcript nor a preserved copy (`13_check_ignore_v_nonmatching_probes.txt`):

```
.gitignore:24:**/paper/*.out	research/pde_ledger_v3/paper/example.out
```

`git ls-files '*.out'` lists zero paths with a `paper` component (`03_classify.txt`). `.gitignore:19-24` calls that pattern LaTeX hyperref aux, ignored where LaTeX runs. `.gitignore`'s hash is unchanged since round 0, so this pattern is not part of the repair. `AGENTS.md:29` still says never `git add -f`. That path does not change which evidence is committed, and it is not a finding.

## Check 2 — the ban, and `git rm --cached`

The bullet forbids staging an engine transcript, defined there as a `.out` in the two directories or a preserved copy under them, and forbids `git add -f` of that object.

Command: `git check-ignore -v --non-matching` on a new file in each directory and a preserved copy under each (`13_check_ignore_v_nonmatching_probes.txt`):

```
.gitignore:28:research/pde_ledger_v3/scripts/out/*.out	research/pde_ledger_v3/scripts/out/NEW_untracked_example.out
.gitignore:29:research/pde_ledger_v3/mathematica/out/*.out	research/pde_ledger_v3/mathematica/out/NEW_untracked_example.out
.gitignore:30:research/pde_ledger_v3/scripts/out/preserved/	research/pde_ledger_v3/scripts/out/preserved/S15_separation_departure.r1.out
.gitignore:31:research/pde_ledger_v3/mathematica/out/preserved/	research/pde_ledger_v3/mathematica/out/preserved/S15_separation_departure.r1.out
```

Command: `git check-ignore -v --no-index --stdin` on all 145 tracked `*.out` paths (`07_check_ignore_v_no_index_all_tracked.txt`, `08_no_index_split.txt`).

```
inside_tracked 50
inside_matched_by_no_index 50
inside_unmatched 0
    34 .gitignore:28:research/pde_ledger_v3/scripts/out/*.out
    16 .gitignore:29:research/pde_ledger_v3/mathematica/out/*.out
tracked_preserved 0
```

Command: `git check-ignore -v --stdin` on the same 145, without `--no-index` (`06_check_ignore_v_all_tracked.txt`):

```
EXIT 1
```

A new `.out` directly in either directory is ignored. Each of the 50 already tracked transcripts matches the same pattern once the index is set aside, and the default check leaves them tracked. Zero tracked paths sit under `preserved/`. A preserved copy at `out/preserved/<name>.<round>.out` is ignored. The prose ban and the ignore rules cover the same two directories and the same preserved copies.

Command: collapse whitespace and compare the `git rm --cached` clause in the current bullet with the round-0 bullet (`12_recon_hash.txt`).

```
rm_clause_in_current True
rm_clause_in_round0 True
```

The clause is: "When a run replaces a transcript that is still tracked, the commit that adds the new manifest removes the old file from the index (`git rm --cached`); history keeps it." The line break moved when the bullet grew. The words are the same. The trigger, the commit that performs `git rm --cached`, and history keeping the file are unchanged.

## Check 3 — the lines the repair touched

The new words are the definition of what may not be staged and the sentence that the other `.out` files are ordinary evidence committed as before. Check 1 shows those other files are the 95 tracked evidence `.out` files, and that a new file in each named class is not ignored. Check 2 shows the ban still covers a new transcript in either directory and a preserved copy under either. The `git rm --cached` clause is the round-0 clause.

`AGENTS.md:31` still names `research/pde_ledger_v3/*/out/*.out`. A five-component match of that shape hits 50 of the tracked `*.out` paths and 0 of the 95 outside paths (`15_agents_glob_shape.txt`):

```
agents_glob_shape_matches_tracked 50
agents_glob_shape_matches_outside 0
preserved_component 0
```

`AGENTS.md:29` still forbids `git add -f` for every ignored path. Both lines are outside the repair; their file hash is the round-0 hash.

No defect in the touched lines changes what is committed, what evidence is in the tree, or what a record can quote.

## Verdict

CLEAR (no finding).
