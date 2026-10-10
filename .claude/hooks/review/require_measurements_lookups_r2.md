# Measurements — measurements-gate repair, review round 2 (generated 2026-10-10 13:57)

Generator: `_scratch/hook_fix/gen_lookups_r2.sh` (sha256sum/grep/sed/cat/git show only). The reviewed r2
scripts are frozen as `_scratch/hook_fix/r2_require_measurements_gate.sh`, `_scratch/hook_fix/r2_install_measurements_gate.sh` and `_scratch/hook_fix/r2_require_measurements.sh`.

## The reviewed version
```
$ sha256sum _scratch/hook_fix/r2_require_measurements_gate.sh _scratch/hook_fix/r2_install_measurements_gate.sh _scratch/hook_fix/r2_require_measurements.sh
0ae77684bfd7fbfefc24003961b22dfe9b1dd40e6f88cb7febd8d91b69b91a5f  _scratch/hook_fix/r2_require_measurements_gate.sh
7e3bb54e05938426464311518d27e41bc82c1cf4dc86459d534273480713c113  _scratch/hook_fix/r2_install_measurements_gate.sh
fae06f491d4245db0f032657a12e8869a44dcb7c2c78d2d983013a2aeb86ab6c  _scratch/hook_fix/r2_require_measurements.sh
```

```
$ cat _scratch/hook_fix/review_baseline_r2.sha256
fae06f491d4245db0f032657a12e8869a44dcb7c2c78d2d983013a2aeb86ab6c  .claude/hooks/require_measurements.sh
0ae77684bfd7fbfefc24003961b22dfe9b1dd40e6f88cb7febd8d91b69b91a5f  .claude/hooks/require_measurements_gate.sh
7e3bb54e05938426464311518d27e41bc82c1cf4dc86459d534273480713c113  .claude/hooks/install_measurements_gate.sh
1cc6be119e6ec0c48596d57aef527d6703cdcfee354dafd6b2018c1b3c43c208  .claude/settings.json
```

```
$ sha256sum _scratch/hook_fix/review_prompt_r2.md
035d0a0d38895353df8cb3f62dc4b3f411372d0fa26651599422b8368a6f292a  _scratch/hook_fix/review_prompt_r2.md
```

## Both verdicts
```
$ grep -n -o 'Verdict: [A-Z][A-Z ]*' _scratch/hook_fix/review_r2_codex_final.txt _scratch/hook_fix/review_r2_grok.txt
_scratch/hook_fix/review_r2_codex_final.txt:1:Verdict: NEEDS REVISION
_scratch/hook_fix/review_r2_grok.txt:1:Verdict: NEEDS REVISION
```

## Prompt reference commit (Codex scope note)
```
$ grep -n '8ca68b2a' _scratch/hook_fix/review_prompt_r2.md
14:The version before this repair began is `git show 8ca68b2a:.claude/hooks/require_measurements.sh`. Earlier rounds,
39:5. **Semantics.** The gated paths and the counterpart path should match the version at `8ca68b2a`. The intended
```

```
$ git show 8ca68b2a:.claude/hooks/require_measurements.sh | grep -n 'DIR_RE='
46:DIR_RE='^research/pde_ledger_v3/directives/[^/]+\.md$'
```

```
$ git show 293e7e53:.claude/hooks/require_measurements.sh | grep -n 'DIR_RE='
50:DIR_RE='^research/pde_ledger_v3/directives/([^/]+\.md|_legs/[^/]*brief[^/]*\.md)$'
```

```
$ git log --format='%h %ad %s' --date=iso -- .claude/hooks/require_measurements.sh
ce377dde 2026-10-10 13:13:06 -0600 Preserve the measurements-gate repair as reviewed in round 1 (not accepted)
db1d7033 2026-10-10 12:21:50 -0600 Preserve the measurements-gate repair as reviewed in round 0 (not accepted)
293e7e53 2026-08-12 11:09:47 -0600 The directive is closed: one census I skipped for four rounds explained three findings
8ca68b2a 2026-08-12 10:21:22 -0600 A gate that refuses, not a reminder that decays
```

## A. Replace refs and grafts hide the stored commit (Codex 4, Grok 1)
```
$ grep -n 'git rev-list\|git diff-tree\|git log\|git cat-file\|GIT_NO_REPLACE\|GRAFT' _scratch/hook_fix/r2_require_measurements_gate.sh
74:  [ "$(git cat-file -t "$new" 2>/dev/null)" = commit ] || continue
94:new_commits=$(git rev-list --stdin --ignore-missing < "$tmp" 2>&1) \
100:  parents=$(git rev-list --no-walk --parents "$c" 2>&1) || refuse "could not read $c: $parents"
110:  git diff-tree "${how[@]}" --no-commit-id --name-only -z "$c" > "$tmp" || refuse "could not diff $c"
121:  msg=$(git log -1 --format=%B "$c") || refuse "could not read the message of $c"
128:    if [ -z "${in_commit[$m]:-}" ] || ! git cat-file -e "$c:$m" 2>/dev/null; then
133:  [ -n "$missing" ] && report+="  $(git log -1 --format='%h %s' "$c")"$'\n'"$missing"
```

```
$ grep -n 'diff-tree with replace\|diff-tree without replace\|UPDATE_EXIT\|^main=\|raw diff' -A1 _scratch/hook_fix/r2_grok_evidence/stdout/fu-replace.out
5:--- diff-tree with replace ---
6-3e3ec9889c13411275cdab3f4d52b65d50bcfd0a
--
11:--- diff-tree without replace ---
12-3e3ec9889c13411275cdab3f4d52b65d50bcfd0a
--
15:UPDATE_EXIT:0
16:main=3e3ec9889c13411275cdab3f4d52b65d50bcfd0a
17:--- raw diff of main vs parent ---
18-research/pde_ledger_v3/directives/claim.md
```

```
$ grep -n 'CONTROL_EXIT' _scratch/hook_fix/r2_grok_evidence/stdout/fu-replace-control.out
18:CONTROL_EXIT:128
```

```
$ grep -n 'BRANCH_EXIT\|claim.md' _scratch/hook_fix/r2_grok_evidence/stdout/fu-graft-silent.out
1:BRANCH_EXIT:0
5:research/pde_ledger_v3/directives/claim.md
```

```
$ grep -n 'replace\|update-ref\|exit=' _scratch/hook_fix/r2_codex_evidence/replacement.stdout
1:$ cd /tmp/measurements-r2-review-wm6j4rma/replacement && git init -q
2:exit=0
3:$ cd /tmp/measurements-r2-review-wm6j4rma/replacement && bash .claude/hooks/install_measurements_gate.sh
4:baseline written: /tmp/measurements-r2-review-wm6j4rma/replacement/.git/measurements-gate-baseline (1 commits)
5:hook linked: /tmp/measurements-r2-review-wm6j4rma/replacement/.git/hooks/reference-transaction -> /tmp/measurements-r2-review-wm6j4rma/replacement/.claude/hooks/require_measurements_gate.sh
6:exit=0
7:$ cd /tmp/measurements-r2-review-wm6j4rma/replacement && git show --format=fuller --stat --no-renames 5b521c6630948e4345f94ad6f5142a632176eb16
18:exit=0
19:$ cd /tmp/measurements-r2-review-wm6j4rma/replacement && git replace 5b521c6630948e4345f94ad6f5142a632176eb16 194889915064b56cb2d218b2e07fe78e4bd29129
20:exit=0
21:$ cd /tmp/measurements-r2-review-wm6j4rma/replacement && git update-ref refs/heads/masked 5b521c6630948e4345f94ad6f5142a632176eb16
22:exit=0
23:$ cd /tmp/measurements-r2-review-wm6j4rma/replacement && git replace -d 5b521c6630948e4345f94ad6f5142a632176eb16
24:Deleted replace ref '5b521c6630948e4345f94ad6f5142a632176eb16'
25:exit=0
26:$ cd /tmp/measurements-r2-review-wm6j4rma/replacement && git show --format=fuller --stat --no-renames refs/heads/masked
37:exit=0
```

## B. rev-list stderr is read as commit ids (Codex 8, Grok 2)
```
$ sed -n '94,95p;100p' _scratch/hook_fix/r2_require_measurements_gate.sh
new_commits=$(git rev-list --stdin --ignore-missing < "$tmp" 2>&1) \
  || refuse "could not list the new commits: $new_commits"
  parents=$(git rev-list --no-walk --parents "$c" 2>&1) || refuse "could not read $c: $parents"
```

```
$ cat _scratch/hook_fix/r2_grok_evidence/stdout/fu-graft-ordinary.out
hint: Support for <GIT_DIR>/info/grafts is deprecated
hint: and will be removed in a future Git version.
hint: 
hint: Please use "git replace --convert-graft-file"
hint: to convert the grafts into replace refs.
hint: 
hint: Turn this message off by running
hint: "git config advice.graftFileDeprecated false"
BLOCKED — require_measurements: could not read hint:: fatal: invalid object name 'hint'.
fatal: ref updates aborted by hook
COMMIT_EXIT:128
```

```
$ grep -n 'hint\|BLOCKED\|exit=' _scratch/hook_fix/r2_codex_evidence/grafts_checks.stdout
2:exit=0
6:exit=0
8:hint: Support for <GIT_DIR>/info/grafts is deprecated
9:hint: and will be removed in a future Git version.
10:hint: 
11:hint: Please use "git replace --convert-graft-file"
12:hint: to convert the grafts into replace refs.
13:hint: 
14:hint: Turn this message off by running
15:hint: "git config advice.graftFileDeprecated false"
16:BLOCKED — require_measurements: could not read hint:: fatal: invalid object name 'hint'.
18:exit=128
20:exit=0
24:exit=0
26:exit=0
28:exit=0
40:exit=0
```

## C. symbolic-ref moves a branch or HEAD without the hook (Codex 1, Grok 3)
```
$ sed -n '20,21p' _scratch/hook_fix/r2_require_measurements_gate.sh
# message cleanup removed. Git runs this hook before every reference update, and
# --no-verify does not skip it. So it judges each commit itself: its own diff
```

```
$ sed -n '40,49p' _scratch/hook_fix/r2_codex_evidence/refs.stdout
$ cd /tmp/measurements-r2-review-wm6j4rma/refs && git symbolic-ref refs/heads/symbolic refs/tags/unsafe
exit=0
$ cd /tmp/measurements-r2-review-wm6j4rma/refs && git rev-parse refs/heads/symbolic
9503e9fca95da5b0fef66ca27685965502aa7c7e
exit=0
$ cd /tmp/measurements-r2-review-wm6j4rma/refs && git symbolic-ref HEAD refs/heads/symbolic
exit=0
$ cd /tmp/measurements-r2-review-wm6j4rma/refs && git rev-parse HEAD
9503e9fca95da5b0fef66ca27685965502aa7c7e
exit=0
```

```
$ grep -n 'SYM_EXIT\|rev-parse main\|symbolic-ref main\|claim.md' _scratch/hook_fix/r2_grok_evidence/stdout/fu-sym.out
2:SYM_EXIT:0
8:research/pde_ledger_v3/directives/claim.md
10:MAIN_SYM_EXIT:0
11:rev-parse main: 3e3ec9889c13411275cdab3f4d52b65d50bcfd0a
12:symbolic-ref main: refs/tags/parked
```

## D. git-annex and synced/* are not judged, and can become HEAD (Codex 2)
```
$ sed -n '5,8p;68p' _scratch/hook_fix/r2_require_measurements_gate.sh
# commands that produced the claim. This refuses to let a branch (or a detached HEAD)
# move onto a commit that adds, changes or deletes a directives/*.md, or a _legs/
# brief, unless that commit adds or changes its _measurements/ counterpart (deleting
# the counterpart does not count), or its message carries the reason below.
    refs/heads/git-annex|refs/heads/synced/*) continue ;;
```

```
$ sed -n '7,11p' _scratch/hook_fix/r2_codex_evidence/corrected_refs.stdout
$ cd /tmp/measurements-r2-review-wm6j4rma/excluded_heads_corrected && git update-ref refs/heads/synced/unsafe da55b668e070165a52394a603521d935cf3b83c3
exit=0
$ cd /tmp/measurements-r2-review-wm6j4rma/excluded_heads_corrected && git checkout -f synced/unsafe
Switched to branch 'synced/unsafe'
exit=0
```

## E. Worktree-qualified HEAD updates are not judged (Codex 3)
```
$ sed -n '69,70p' _scratch/hook_fix/r2_require_measurements_gate.sh
    HEAD|refs/heads/*) ;;
    *) continue ;;
```

```
$ sed -n '11,12p' _scratch/hook_fix/r2_codex_evidence/worktree_refs.stdout
$ cd /tmp/measurements-r2-review-wm6j4rma/worktree_refs && git update-ref worktrees/linked_head/HEAD 8ae719ad1cd5685e972651d523dfa8edfdc958de
exit=0
```

## F. The hatch is read from rendered log output (Codex 5)
```
$ sed -n '120,121p' _scratch/hook_fix/r2_require_measurements_gate.sh
  # Escape hatch: an explicit, recorded reason in the message as git stored it.
  msg=$(git log -1 --format=%B "$c") || refuse "could not read the message of $c"
```

```
$ grep -n 'gpg.program\|showSignature\|gpg: no-measurements\|update-ref\|exit=' _scratch/hook_fix/r2_codex_evidence/signature_output.stdout
2:exit=0
6:exit=0
10:exit=0
11:$ cd /tmp/measurements-r2-review-wm6j4rma/signature_output && git config gpg.program /tmp/measurements-r2-review-wm6j4rma/signature_output/verify-test.sh
12:exit=0
13:$ cd /tmp/measurements-r2-review-wm6j4rma/signature_output && git config log.showSignature true
14:exit=0
25:exit=0
27:gpg: no-measurements: external verifier diagnostic
30:exit=0
31:$ cd /tmp/measurements-r2-review-wm6j4rma/signature_output && git update-ref refs/heads/signed 4bbb1e1fad861cf5d932c8727c163c6b6e191376
32:exit=0
```

```
$ grep -a -n 'logOutputEncoding\|BLOCKED\|exit=' _scratch/hook_fix/r2_codex_evidence/raw_messages.stdout
2:exit=0
6:exit=0
7:$ cd /tmp/measurements-r2-review-wm6j4rma/ebcdic_message && git config i18n.logOutputEncoding IBM1047
8:exit=0
16:exit=0
19:exit=0
22:BLOCKED — CLAUDE.md E1 (orchestrator half).
39:exit=128
41:exit=0
45:exit=0
49:exit=0
58:exit=0
60:BLOCKED — CLAUDE.md E1 (orchestrator half).
77:exit=128
81:exit=0
89:exit=0
91:BLOCKED — CLAUDE.md E1 (orchestrator half).
108:exit=128
```

## G. A historical checkout leaves the gate link dangling (Codex 6)
```
$ sed -n '10,11p;16,19p;41p' _scratch/hook_fix/r2_install_measurements_gate.sh
# 2. Links git's reference-transaction hook to require_measurements_gate.sh. Linked
#    worktrees share the main checkout's hooks, so run this once per clone.
main=${common%/.git}                 # the main checkout; linked worktrees use its gate
hooks=$(git rev-parse --path-format=absolute --git-path hooks)
baseline="$common/measurements-gate-baseline"
gate="$main/.claude/hooks/require_measurements_gate.sh"
  ln -s "$gate" "$link"
```

```
$ cat _scratch/hook_fix/r2_codex_evidence/historical_checkout.stdout
$ cd /tmp/measurements-r2-review-wm6j4rma/historical_checkout && git init -q
exit=0
$ cd /tmp/measurements-r2-review-wm6j4rma/historical_checkout && bash .claude/hooks/install_measurements_gate.sh
baseline written: /tmp/measurements-r2-review-wm6j4rma/historical_checkout/.git/measurements-gate-baseline (1 commits)
hook linked: /tmp/measurements-r2-review-wm6j4rma/historical_checkout/.git/hooks/reference-transaction -> /tmp/measurements-r2-review-wm6j4rma/historical_checkout/.claude/hooks/require_measurements_gate.sh
exit=0
$ cd /tmp/measurements-r2-review-wm6j4rma/historical_checkout && git checkout --detach 06f04572410083893e85032e028d58d1bba0b24e
HEAD is now at 06f0457 history before new hooks
exit=0
$ cd /tmp/measurements-r2-review-wm6j4rma/historical_checkout && ls -l .git/hooks/reference-transaction
lrwxrwxrwx 1 trevnorris trevnorris 99 Oct 10 13:19 .git/hooks/reference-transaction -> /tmp/measurements-r2-review-wm6j4rma/historical_checkout/.claude/hooks/require_measurements_gate.sh
exit=0
gate target exists=False
$ cd /tmp/measurements-r2-review-wm6j4rma/historical_checkout && git commit -m 'unmeasured after historical checkout'
[detached HEAD 9575c4e] unmeasured after historical checkout
 1 file changed, 1 insertion(+)
 create mode 100644 research/pde_ledger_v3/directives/after-checkout.md
exit=0
$ cd /tmp/measurements-r2-review-wm6j4rma/historical_checkout && git show --format=fuller --stat --no-renames HEAD
commit 9575c4e65fc19f71bfaa091e8bc223ae1bbff3e2
Author:     Review <review@example.invalid>
AuthorDate: Sat Oct 10 13:19:28 2026 -0600
Commit:     Review <review@example.invalid>
CommitDate: Sat Oct 10 13:19:28 2026 -0600

    unmeasured after historical checkout

 research/pde_ledger_v3/directives/after-checkout.md | 1 +
 1 file changed, 1 insertion(+)
exit=0
```

```
$ cat _scratch/hook_fix/r2_codex_evidence/pretool_old_checkout.stdout
$ cd /tmp/measurements-r2-review-wm6j4rma/historical_checkout && bash -c '"$CLAUDE_PROJECT_DIR/.claude/hooks/require_measurements.sh"'
stdin=b'{"tool_input": {"command": "git status"}}'
bash: line 1: /tmp/measurements-r2-review-wm6j4rma/historical_checkout/.claude/hooks/require_measurements.sh: No such file or directory
exit=127
```

## H. The baseline omits detached HEADs (Codex 7, Grok 4)
```
$ sed -n '6,9p;27,29p' _scratch/hook_fix/r2_install_measurements_gate.sh
# 1. Records the baseline, once: the commit every ref points at now. The gate treats
#    those commits, and everything they reach, as already judged, so history from
#    before the gate is not re-judged. An existing baseline is never overwritten,
#    because refreshing it would grandfather whatever arrived since.
  git for-each-ref --format='%(objectname)' \
    | while read -r o; do git rev-parse -q --verify "$o^{commit}" 2>/dev/null || true; done \
    | sort -u > "$tmp"
```

```
$ grep -n 'detach\|BLOCKED\|aborted\|exit=' _scratch/hook_fix/r2_codex_evidence/corrected_detached.stdout
1:$ cd /tmp/measurements-r2-review-wm6j4rma/detached_baseline_corrected && git init -q
2:exit=0
3:$ cd /tmp/measurements-r2-review-wm6j4rma/detached_baseline_corrected && git checkout -f --detach da55b668e070165a52394a603521d935cf3b83c3
5:exit=0
6:$ cd /tmp/measurements-r2-review-wm6j4rma/detached_baseline_corrected && git rev-parse HEAD
8:exit=0
9:$ cd /tmp/measurements-r2-review-wm6j4rma/detached_baseline_corrected && bash .claude/hooks/install_measurements_gate.sh
10:baseline written: /tmp/measurements-r2-review-wm6j4rma/detached_baseline_corrected/.git/measurements-gate-baseline (1 commits)
11:hook linked: /tmp/measurements-r2-review-wm6j4rma/detached_baseline_corrected/.git/hooks/reference-transaction -> /tmp/measurements-r2-review-wm6j4rma/detached_baseline_corrected/.claude/hooks/require_measurements_gate.sh
12:exit=0
13:$ cd /tmp/measurements-r2-review-wm6j4rma/detached_baseline_corrected && cat .git/measurements-gate-baseline
15:exit=0
16:$ cd /tmp/measurements-r2-review-wm6j4rma/detached_baseline_corrected && git checkout master
28:exit=0
29:$ cd /tmp/measurements-r2-review-wm6j4rma/detached_baseline_corrected && git checkout --detach da55b668e070165a52394a603521d935cf3b83c3
30:BLOCKED — CLAUDE.md E1 (orchestrator half).
46:fatal: ref updates aborted by hook
47:exit=128
48:$ cd /tmp/measurements-r2-review-wm6j4rma/detached_baseline_corrected && git worktree add --detach /tmp/measurements-r2-review-wm6j4rma/old_worktree_corrected da55b668e070165a52394a603521d935cf3b83c3
49:Preparing worktree (detached HEAD da55b66)
50:BLOCKED — CLAUDE.md E1 (orchestrator half).
66:fatal: ref updates aborted by hook
67:exit=128
68:$ cd /tmp/measurements-r2-review-wm6j4rma/detached_worktree_baseline && git init -q
69:exit=0
70:$ cd /tmp/measurements-r2-review-wm6j4rma/detached_worktree_baseline && git worktree add --detach /tmp/measurements-r2-review-wm6j4rma/preexisting_detached_worktree da55b668e070165a52394a603521d935cf3b83c3
71:Preparing worktree (detached HEAD da55b66)
73:exit=0
74:$ cd /tmp/measurements-r2-review-wm6j4rma/preexisting_detached_worktree && bash .claude/hooks/install_measurements_gate.sh
75:baseline written: /tmp/measurements-r2-review-wm6j4rma/detached_worktree_baseline/.git/measurements-gate-baseline (1 commits)
76:hook linked: /tmp/measurements-r2-review-wm6j4rma/detached_worktree_baseline/.git/hooks/reference-transaction -> /tmp/measurements-r2-review-wm6j4rma/detached_worktree_baseline/.claude/hooks/require_measurements_gate.sh
77:exit=0
78:$ cd /tmp/measurements-r2-review-wm6j4rma/detached_worktree_baseline && cat .git/measurements-gate-baseline
80:exit=0
81:$ cd /tmp/measurements-r2-review-wm6j4rma/detached_worktree_baseline && git branch old-worktree-history da55b668e070165a52394a603521d935cf3b83c3
82:BLOCKED — CLAUDE.md E1 (orchestrator half).
98:fatal: ref updates aborted by hook
99:exit=128
```

```
$ cat _scratch/hook_fix/r2_grok_evidence/stdout/fu-detach-back.out
PRE=ddce3cb0e4baaa3f9dd3b4ad44a3b664590f7520
baseline:
f29b313ce13c70b5ddc596435ed43ab956b0eeba
PRE NOT in baseline
CHECKOUT_MAIN_EXIT:0
now main f29b313ce13c70b5ddc596435ed43ab956b0eeba
BLOCKED — CLAUDE.md E1 (orchestrator half).

These commits change directives with no matching _measurements/ file in the same commit:

  ddce3cb detached pre-install directive
    research/pde_ledger_v3/directives/pre.md  ->  research/pde_ledger_v3/directives/_measurements/pre.md

A claim about an artifact carries the command that produced it. Run the commands,
write them and their LITERAL output to the path above, stage it, and commit again.
Regenerate the file from the commands; do not transcribe.

Measured 2026-08-12, the cost of skipping this: four export-chain designs and eight
review legs died on a question one `len(LEDGER)` answered.

If this document genuinely asserts nothing about any artifact, say so in the commit
message with:  no-measurements: <reason>
fatal: ref updates aborted by hook
CHECKOUT_BACK_EXIT:128
HEAD=f29b313ce13c70b5ddc596435ed43ab956b0eeba
```

```
$ grep -n 'NOT in baseline\|BACK' _scratch/hook_fix/r2_grok_evidence/stdout/fu2-wt-back.out
2:WTPRE NOT in baseline
22:BACK:128
```

## I. A dangling baseline symlink is overwritten (Codex 9)
```
$ sed -n '8p;23p;30p' _scratch/hook_fix/r2_install_measurements_gate.sh
#    before the gate is not re-judged. An existing baseline is never overwritten,
if [ -e "$baseline" ]; then
  mv "$tmp" "$baseline"
```

```
$ grep -n 'measurements-gate-baseline' _scratch/hook_fix/r2_codex_evidence/installer.stdout
4:baseline written: /tmp/measurements-r2-review-wm6j4rma/installer_unrelated/.git/measurements-gate-baseline (1 commits)
11:baseline written: /tmp/measurements-r2-review-wm6j4rma/installer_repeat/.git/measurements-gate-baseline (1 commits)
15:baseline kept: /tmp/measurements-r2-review-wm6j4rma/installer_repeat/.git/measurements-gate-baseline (1 commits)
21:$ cd /tmp/measurements-r2-review-wm6j4rma/installer_dangling && ls -l .git/measurements-gate-baseline
22:lrwxrwxrwx 1 trevnorris trevnorris 16 Oct 10 13:16 .git/measurements-gate-baseline -> missing-baseline
25:baseline written: /tmp/measurements-r2-review-wm6j4rma/installer_dangling/.git/measurements-gate-baseline (1 commits)
28:$ cd /tmp/measurements-r2-review-wm6j4rma/installer_dangling && ls -l .git/measurements-gate-baseline
29:-rw------- 1 trevnorris trevnorris 41 Oct 10 13:16 .git/measurements-gate-baseline
```

