# Measurements — measurements-gate repair, review round 1 (generated 2026-10-10 13:08)

Generator: `_scratch/hook_fix/gen_lookups_r1.sh` (sha256sum/grep/sed/cat only). The reviewed r1 scripts are
frozen as `_scratch/hook_fix/r1_require_measurements_gate.sh` and `_scratch/hook_fix/r1_require_measurements.sh`.

## The reviewed version
```
$ sha256sum _scratch/hook_fix/r1_require_measurements_gate.sh _scratch/hook_fix/r1_require_measurements.sh
333263c145e26782cbca41e15ba905d96b923abb7de9049c1526b4b4f0a15a3b  _scratch/hook_fix/r1_require_measurements_gate.sh
658cdcd09ef1af81c69a6c401cbab01151e01541b18b61cf593eb939247ab66d  _scratch/hook_fix/r1_require_measurements.sh
```

```
$ cat _scratch/hook_fix/review_baseline_r1.sha256
658cdcd09ef1af81c69a6c401cbab01151e01541b18b61cf593eb939247ab66d  .claude/hooks/require_measurements.sh
333263c145e26782cbca41e15ba905d96b923abb7de9049c1526b4b4f0a15a3b  .claude/hooks/require_measurements_gate.sh
1cc6be119e6ec0c48596d57aef527d6703cdcfee354dafd6b2018c1b3c43c208  .claude/settings.json
```

```
$ sha256sum _scratch/hook_fix/review_prompt_r1.md
09278d8fdd4751326883950169b6a659b98c6caacafe8900abcb8bdaaf886f63  _scratch/hook_fix/review_prompt_r1.md
```

## Both verdicts
```
$ grep -n -o 'Verdict: [A-Z][A-Z ]*' _scratch/hook_fix/review_r1_codex_final.txt _scratch/hook_fix/review_r1_grok.txt
_scratch/hook_fix/review_r1_codex_final.txt:1:Verdict: NEEDS REVISION
_scratch/hook_fix/review_r1_grok.txt:3:Verdict: NEEDS REVISION
```

## A. Commits reached by a ref the gate does not check are skipped (Codex 1, Grok 2)
```
$ grep -n 'not --all\|refs/heads/git-annex) continue\|HEAD|refs/heads/\*)' _scratch/hook_fix/r1_require_measurements_gate.sh
58:    refs/heads/git-annex) continue ;;
59:    HEAD|refs/heads/*) ;;
68:new_commits=$(git rev-list "${tips[@]}" --not --all 2>&1) \
```

```
$ cat _scratch/hook_fix/r1_codex_evidence/26_stash_branch.excerpt.stdout
cwd=/tmp/measurements-gate-r1-d7_levce/repos/stash_branch
$ git stash push -m 'directive work in progress'
stdout:
Saved working directory and index state On main: directive work in progress

stderr:
<empty>

exit=0
elapsed_seconds=0.019587
cwd=/tmp/measurements-gate-r1-d7_levce/repos/stash_branch
$ git branch recovered 'stash^2'
stdout:
<empty>

stderr:
<empty>

exit=0
elapsed_seconds=0.007513
cwd=/tmp/measurements-gate-r1-d7_levce/repos/stash_branch
$ git switch recovered
stdout:
<empty>

stderr:
Switched to branch 'recovered'

exit=0
elapsed_seconds=0.001891
cwd=/tmp/measurements-gate-r1-d7_levce/repos/stash_branch
$ git show --format=fuller --name-status HEAD
stdout:
commit 4075835c9d5b411652f011413b91acdc9fe69955
Author:     Hook Reviewer <hook-review@example.invalid>
AuthorDate: Sat Oct 10 12:00:00 2026 +0000
Commit:     Hook Reviewer <hook-review@example.invalid>
CommitDate: Sat Oct 10 12:00:00 2026 +0000

    index on main: 83be10d seed

A	research/pde_ledger_v3/directives/test.md

stderr:
<empty>

exit=0
elapsed_seconds=0.001359
cwd=/tmp/measurements-gate-r1-d7_levce/repos/stash_branch
$ git cat-file -e HEAD:research/pde_ledger_v3/directives/_measurements/test.md
stdout:
<empty>

stderr:
fatal: Not a valid object name HEAD:research/pde_ledger_v3/directives/_measurements/test.md

exit=128
elapsed_seconds=0.001239
```

```
$ grep -n 'TAG EXIT\|BRANCH EXIT\|MAIN MOVED\|PULL EXIT\|PULLED DIRECTIVE' _scratch/hook_fix/r1_grok_evidence/stdout/followup.txt _scratch/hook_fix/r1_grok_evidence/stdout/*.txt
_scratch/hook_fix/r1_grok_evidence/stdout/run_tests.txt:885:TAG EXIT:0
_scratch/hook_fix/r1_grok_evidence/stdout/run_tests.txt:887:BRANCH EXIT:0
_scratch/hook_fix/r1_grok_evidence/stdout/run_tests.txt:889:MAIN MOVED TO BAD COMMIT
_scratch/hook_fix/r1_grok_evidence/stdout/run_tests.txt:920:PULL EXIT:0
_scratch/hook_fix/r1_grok_evidence/stdout/run_tests.txt:921:PULLED DIRECTIVE IS ON MAIN
```

## B. Rename detection keeps only the destination (Codex 2, Grok 3)
```
$ grep -n 'how=(-r -M --root)' _scratch/hook_fix/r1_require_measurements_gate.sh
79:    how=(-r -M --root)
```

```
$ cat _scratch/hook_fix/r1_codex_evidence/04_rename_out.excerpt.stdout
cwd=/tmp/measurements-gate-r1-d7_levce/repos/rename_out
$ git mv research/pde_ledger_v3/directives/test.md archived.md
stdout:
<empty>

stderr:
<empty>

exit=0
elapsed_seconds=0.001394
cwd=/tmp/measurements-gate-r1-d7_levce/repos/rename_out
$ git commit -m 'move directive out'
stdout:
[main 88d6263] move directive out
 1 file changed, 0 insertions(+), 0 deletions(-)
 rename research/pde_ledger_v3/directives/test.md => archived.md (100%)

stderr:
<empty>

exit=0
elapsed_seconds=0.017315
cwd=/tmp/measurements-gate-r1-d7_levce/repos/rename_out
$ git show --format=fuller --name-status HEAD
stdout:
commit 88d62635c9a91d056f8832d884cb0fd419af2d91
Author:     Hook Reviewer <hook-review@example.invalid>
AuthorDate: Sat Oct 10 12:00:00 2026 +0000
Commit:     Hook Reviewer <hook-review@example.invalid>
CommitDate: Sat Oct 10 12:00:00 2026 +0000

    move directive out

R100	research/pde_ledger_v3/directives/test.md	archived.md

stderr:
<empty>

exit=0
elapsed_seconds=0.001430
cwd=/tmp/measurements-gate-r1-d7_levce/repos/rename_out
$ git diff-tree -r -M --root --no-commit-id --name-only HEAD
stdout:
archived.md

stderr:
<empty>

exit=0
elapsed_seconds=0.001360
cwd=/tmp/measurements-gate-r1-d7_levce/repos/rename_out
$ git diff-tree -r --no-renames --root --no-commit-id --name-status HEAD
stdout:
A	archived.md
D	research/pde_ledger_v3/directives/test.md

stderr:
<empty>

exit=0
elapsed_seconds=0.001337
```

## C. NUL records are turned into lines and matched in the ambient locale (Codex 3, 4; Grok 1)
```
$ grep -n "tr '\\\\0' '\\\\n'\|grep -E \"\$DIR_RE\"" _scratch/hook_fix/r1_require_measurements_gate.sh
86:  files=$(git diff-tree "${how[@]}" --no-commit-id --name-only -z "$c" | tr '\0' '\n') \
88:  kept=$(git diff-tree "${how[@]}" --diff-filter=d --no-commit-id --name-only -z "$c" | tr '\0' '\n') \
91:  docs=$(grep -E "$DIR_RE" <<< "$files" || true)
```

```
$ cat _scratch/hook_fix/r1_codex_evidence/05_newline_directive.excerpt.stdout _scratch/hook_fix/r1_codex_evidence/06_fake_counterpart.excerpt.stdout
write '/tmp/measurements-gate-r1-d7_levce/repos/newline_doc/research/pde_ledger_v3/directives/a\nb.md' contents='claim\n'
cwd=/tmp/measurements-gate-r1-d7_levce/repos/newline_doc
$ git add -A
stdout:
<empty>

stderr:
<empty>

exit=0
elapsed_seconds=0.001555
cwd=/tmp/measurements-gate-r1-d7_levce/repos/newline_doc
$ git commit -m 'newline directive without measurements'
stdout:
[main 0794997] newline directive without measurements
 1 file changed, 1 insertion(+)
 create mode 100644 "research/pde_ledger_v3/directives/a\nb.md"

stderr:
<empty>

exit=0
elapsed_seconds=0.017540
cwd=/tmp/measurements-gate-r1-d7_levce/repos/newline_doc
$ git show --format=fuller --name-status HEAD
stdout:
commit 07949970e78dba4aaffde52f31e7efe46297c8e7
Author:     Hook Reviewer <hook-review@example.invalid>
AuthorDate: Sat Oct 10 12:00:00 2026 +0000
Commit:     Hook Reviewer <hook-review@example.invalid>
CommitDate: Sat Oct 10 12:00:00 2026 +0000

    newline directive without measurements

A	"research/pde_ledger_v3/directives/a\nb.md"

stderr:
<empty>

exit=0
elapsed_seconds=0.001425
cwd=/tmp/measurements-gate-r1-d7_levce/repos/newline_doc
$ git ls-tree -r --name-only HEAD
stdout:
.claude/hooks/require_measurements.sh
.claude/hooks/require_measurements_gate.sh
"research/pde_ledger_v3/directives/a\nb.md"
seed.txt

stderr:
<empty>

exit=0
elapsed_seconds=0.001289
write '/tmp/measurements-gate-r1-d7_levce/repos/fake_counterpart/research/pde_ledger_v3/directives/test.md' contents='claim\n'
write '/tmp/measurements-gate-r1-d7_levce/repos/fake_counterpart/junk\nresearch/pde_ledger_v3/directives/_measurements/test.md' contents='fake\n'
cwd=/tmp/measurements-gate-r1-d7_levce/repos/fake_counterpart
$ git add -A
stdout:
<empty>

stderr:
<empty>

exit=0
elapsed_seconds=0.001776
cwd=/tmp/measurements-gate-r1-d7_levce/repos/fake_counterpart
$ git commit -m 'fake counterpart in newline path'
stdout:
[main 6b96759] fake counterpart in newline path
 2 files changed, 2 insertions(+)
 create mode 100644 "junk\nresearch/pde_ledger_v3/directives/_measurements/test.md"
 create mode 100644 research/pde_ledger_v3/directives/test.md

stderr:
<empty>

exit=0
elapsed_seconds=0.022008
cwd=/tmp/measurements-gate-r1-d7_levce/repos/fake_counterpart
$ git show --format=fuller --name-status HEAD
stdout:
commit 6b96759ae8c3e83642fccf25004d0864e229aff9
Author:     Hook Reviewer <hook-review@example.invalid>
AuthorDate: Sat Oct 10 12:00:00 2026 +0000
Commit:     Hook Reviewer <hook-review@example.invalid>
CommitDate: Sat Oct 10 12:00:00 2026 +0000

    fake counterpart in newline path

A	"junk\nresearch/pde_ledger_v3/directives/_measurements/test.md"
A	research/pde_ledger_v3/directives/test.md

stderr:
<empty>

exit=0
elapsed_seconds=0.001491
cwd=/tmp/measurements-gate-r1-d7_levce/repos/fake_counterpart
$ git cat-file -e HEAD:research/pde_ledger_v3/directives/_measurements/test.md
stdout:
<empty>

stderr:
fatal: Not a valid object name HEAD:research/pde_ledger_v3/directives/_measurements/test.md

exit=128
elapsed_seconds=0.001204
```

```
$ cat _scratch/hook_fix/r1_codex_evidence/20_newline_false_positive.excerpt.stdout
write '/tmp/measurements-gate-r1-d7_levce/repos/newline_false_positive/unrelated\nresearch/pde_ledger_v3/directives/test.md' contents='ordinary file\n'
cwd=/tmp/measurements-gate-r1-d7_levce/repos/newline_false_positive
$ git add -A
stdout:
<empty>

stderr:
<empty>

exit=0
elapsed_seconds=0.001744
cwd=/tmp/measurements-gate-r1-d7_levce/repos/newline_false_positive
$ git commit -m 'ordinary newline path'
stdout:
<empty>

stderr:
BLOCKED — CLAUDE.md E1 (orchestrator half).

These commits change directives with no matching _measurements/ file in the same commit:

  7cc3db1 ordinary newline path
    research/pde_ledger_v3/directives/test.md  ->  research/pde_ledger_v3/directives/_measurements/test.md

A claim about an artifact carries the command that produced it. Run the commands,
write them and their LITERAL output to the path above, stage it, and commit again.
Regenerate the file from the commands; do not transcribe.

Measured 2026-08-12, the cost of skipping this: four export-chain designs and eight
review legs died on a question one `len(LEDGER)` answered.

If this document genuinely asserts nothing about any artifact, say so in the commit
message with:  no-measurements: <reason>
fatal: ref updates aborted by hook

exit=128
elapsed_seconds=0.023522
cwd=/tmp/measurements-gate-r1-d7_levce/repos/newline_false_positive
$ git diff --cached --name-status
stdout:
A	"unrelated\nresearch/pde_ledger_v3/directives/test.md"

stderr:
<empty>

exit=0
elapsed_seconds=0.001417
```

```
$ cat _scratch/hook_fix/r1_codex_evidence/28_invalid_utf8.excerpt.stdout
wrote byte path='/tmp/measurements-gate-r1-d7_levce/repos/invalid_utf8/research/pde_ledger_v3/directives/invalid_\udcff.md'
cwd=/tmp/measurements-gate-r1-d7_levce/repos/invalid_utf8
$ git add -A
stdout:
<empty>

stderr:
<empty>

exit=0
elapsed_seconds=0.001486
cwd=/tmp/measurements-gate-r1-d7_levce/repos/invalid_utf8
$ git commit -m 'invalid UTF-8 filename'
stdout:
[main c4e174a] invalid UTF-8 filename
 1 file changed, 1 insertion(+)
 create mode 100644 "research/pde_ledger_v3/directives/invalid_\377.md"

stderr:
<empty>

exit=0
elapsed_seconds=0.017180
```

## D. The hook-path rule matches any substring (Codex 5)
```
$ grep -n 'grep -qiE' _scratch/hook_fix/r1_require_measurements.sh
42:if grep -qiE 'hookspath|GIT_CONFIG' <<< "$cmd"; then
```

```
$ grep -n -A4 'hookspathology' _scratch/hook_fix/r1_codex_evidence/10_pre_installation.stdout
130:stdin='{"tool_input": {"command": "printf hookspathology"}}'
131-environment overrides={'CLAUDE_PROJECT_DIR': '/tmp/measurements-gate-r1-d7_levce/repos/pre_installation'}
132-stdout:
133-<empty>
134-
```

## E. The gate's comment on what PreToolUse refuses (Codex 6, Grok 4)
```
$ grep -n 'refuses commands that touch the hook path' _scratch/hook_fix/r1_require_measurements_gate.sh
32:# transcript. The PreToolUse hook refuses commands that touch the hook path.
```

```
$ grep -n 'Removing the hook inside the same command' _scratch/hook_fix/r1_require_measurements.sh
25:# (linked worktrees share its hook). Removing the hook inside the same command that
```

```
$ cat _scratch/hook_fix/r1_codex_evidence/19_comment_and_overmatch.excerpt.stdout
cwd=/tmp/measurements-gate-r1-d7_levce/repos/comment
$ bash /tmp/measurements-gate-r1-d7_levce/repos/comment/.claude/hooks/require_measurements.sh
stdin='{"tool_input": {"command": "rm .git/hooks/reference-transaction"}}'
environment overrides={'CLAUDE_PROJECT_DIR': '/tmp/measurements-gate-r1-d7_levce/repos/comment'}
stdout:
<empty>

stderr:
<empty>

exit=0
elapsed_seconds=0.028378
cwd=/tmp/measurements-gate-r1-d7_levce/repos/comment
$ bash /tmp/measurements-gate-r1-d7_levce/repos/comment/.claude/hooks/require_measurements.sh
stdin='{"tool_input": {"command": "chmod -x .claude/hooks/require_measurements_gate.sh"}}'
environment overrides={'CLAUDE_PROJECT_DIR': '/tmp/measurements-gate-r1-d7_levce/repos/comment'}
stdout:
<empty>

stderr:
<empty>

exit=0
elapsed_seconds=0.028491
cwd=/tmp/measurements-gate-r1-d7_levce/repos/comment
$ bash /tmp/measurements-gate-r1-d7_levce/repos/comment/.claude/hooks/require_measurements.sh
stdin='{"tool_input": {"command": "printf \\"%s\\\\n\\" hookspathology"}}'
environment overrides={'CLAUDE_PROJECT_DIR': '/tmp/measurements-gate-r1-d7_levce/repos/comment'}
stdout:
<empty>

stderr:
BLOCKED — require_measurements: this command mentions core.hooksPath or a GIT_CONFIG*
variable. Either can point git at another hook directory and switch off the
measurements gate (git's reference-transaction hook), so neither is allowed from
here. Restate the command without it.

exit=2
elapsed_seconds=0.028236
cwd=/tmp/measurements-gate-r1-d7_levce/repos/comment
$ bash /tmp/measurements-gate-r1-d7_levce/repos/comment/.claude/hooks/require_measurements.sh
stdin='{"tool_input": {"command": "mkdir -p my-hookspath-backup"}}'
environment overrides={'CLAUDE_PROJECT_DIR': '/tmp/measurements-gate-r1-d7_levce/repos/comment'}
stdout:
<empty>

stderr:
BLOCKED — require_measurements: this command mentions core.hooksPath or a GIT_CONFIG*
variable. Either can point git at another hook directory and switch off the
measurements gate (git's reference-transaction hook), so neither is allowed from
here. Restate the command without it.

exit=2
elapsed_seconds=0.028629
cwd=/tmp/measurements-gate-r1-d7_levce/repos/comment
$ bash -c 'printf "%s\n" hookspathology'
stdout:
hookspathology

stderr:
<empty>

exit=0
elapsed_seconds=0.001237
cwd=/tmp/measurements-gate-r1-d7_levce/repos/comment
$ git config alias.record commit
stdout:
<empty>

stderr:
<empty>

exit=0
elapsed_seconds=0.001401
write '/tmp/measurements-gate-r1-d7_levce/repos/comment/research/pde_ledger_v3/directives/test.md' contents='claim\n'
cwd=/tmp/measurements-gate-r1-d7_levce/repos/comment
$ git add -A
stdout:
<empty>

stderr:
<empty>

exit=0
elapsed_seconds=0.001615
cwd=/tmp/measurements-gate-r1-d7_levce/repos/comment
$ git record --no-verify -m 'unmeasured alias'
stdout:
<empty>

stderr:
BLOCKED — CLAUDE.md E1 (orchestrator half).

These commits change directives with no matching _measurements/ file in the same commit:

  1b2a7b8 unmeasured alias
    research/pde_ledger_v3/directives/test.md  ->  research/pde_ledger_v3/directives/_measurements/test.md

A claim about an artifact carries the command that produced it. Run the commands,
write them and their LITERAL output to the path above, stage it, and commit again.
Regenerate the file from the commands; do not transcribe.

Measured 2026-08-12, the cost of skipping this: four export-chain designs and eight
review legs died on a question one `len(LEDGER)` answered.

If this document genuinely asserts nothing about any artifact, say so in the commit
message with:  no-measurements: <reason>
fatal: ref updates aborted by hook

exit=128
elapsed_seconds=0.029001
cwd=/tmp/measurements-gate-r1-d7_levce/repos/comment
$ /tmp/measurements-gate-r1-d7_levce/commit_from_script.sh
stdout:
<empty>

stderr:
BLOCKED — CLAUDE.md E1 (orchestrator half).

These commits change directives with no matching _measurements/ file in the same commit:

  702f7a4 unmeasured script
    research/pde_ledger_v3/directives/test.md  ->  research/pde_ledger_v3/directives/_measurements/test.md

A claim about an artifact carries the command that produced it. Run the commands,
write them and their LITERAL output to the path above, stage it, and commit again.
Regenerate the file from the commands; do not transcribe.

Measured 2026-08-12, the cost of skipping this: four export-chain designs and eight
review legs died on a question one `len(LEDGER)` answered.

If this document genuinely asserts nothing about any artifact, say so in the commit
message with:  no-measurements: <reason>
fatal: ref updates aborted by hook

exit=128
elapsed_seconds=0.024235
```

