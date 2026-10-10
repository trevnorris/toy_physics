# Measurements — measurements-gate repair, review round 0 (generated 2026-10-10 12:22)

Generator: `_scratch/hook_fix/gen_lookups_r0.sh` (sha256sum/grep/sed/cat only). The reviewed r0 scripts are
frozen from `db1d7033` as `_scratch/hook_fix/r0_require_measurements.sh` and `_scratch/hook_fix/r0_require_measurements_commit_msg.sh`.

## The reviewed version
```
$ sha256sum _scratch/hook_fix/r0_require_measurements.sh _scratch/hook_fix/r0_require_measurements_commit_msg.sh
be3f46bb1e2b4621831132f93d6744134e5554397823b405fc425ab594e80b8a  _scratch/hook_fix/r0_require_measurements.sh
c9d2a31ab66b9030edcdea1f28eae092564eef24b041f8f54c1564b85f96ca91  _scratch/hook_fix/r0_require_measurements_commit_msg.sh
```

```
$ cat _scratch/hook_fix/review_baseline_r0.sha256
be3f46bb1e2b4621831132f93d6744134e5554397823b405fc425ab594e80b8a  .claude/hooks/require_measurements.sh
c9d2a31ab66b9030edcdea1f28eae092564eef24b041f8f54c1564b85f96ca91  .claude/hooks/require_measurements_commit_msg.sh
1cc6be119e6ec0c48596d57aef527d6703cdcfee354dafd6b2018c1b3c43c208  .claude/settings.json
```

```
$ sha256sum _scratch/hook_fix/review_prompt_r0.md
74ea34b9a51ed1bfa9eebacda9775dfaaca82d56fc7699194dcfa6d2efa80af9  _scratch/hook_fix/review_prompt_r0.md
```

## Both verdicts
```
$ grep -n -o 'Verdict: [A-Z][A-Z ]*' _scratch/hook_fix/review_r0_codex_final.txt _scratch/hook_fix/review_r0_grok.txt
_scratch/hook_fix/review_r0_codex_final.txt:1:Verdict: NEEDS REVISION
_scratch/hook_fix/review_r0_grok.txt:1:Verdict: NEEDS REVISION
```

## The commit-msg gate compares with HEAD and greps the uncleaned message (Codex 1, 2)
```
$ grep -n 'git diff --cached --name-only\|no-measurements:\[\[' _scratch/hook_fix/r0_require_measurements_commit_msg.sh
42:if ! staged=$(git diff --cached --name-only 2>&1); then
53:   && grep -qiE 'no-measurements:[[:space:]]*[^[:space:]]' "$msg_file"; then
```

```
$ grep -n -A12 'git commit --amend --no-edit' _scratch/hook_fix/r0_codex_evidence/partial_amend_loses_measurement.stdout
48:STDIN: {"tool_input": {"command": "git commit --amend --no-edit"}}
49-STDOUT:
50-STDERR:
51-EXIT: 0
52:COMMAND: /bin/bash -c 'git commit --amend --no-edit'
53-CWD: /tmp/measurements-gate-review.jfOaBi7z/partial_amend_loses_measurement
54-STDOUT:
55-[master 013e616] paired claim
56- Date: Sat Oct 10 11:52:09 2026 -0600
57- 1 file changed, 1 insertion(+)
58- create mode 100644 research/pde_ledger_v3/directives/claim.md
59-STDERR:
60-EXIT: 0
61-COMMAND: git show --format= --name-status HEAD
62-CWD: /tmp/measurements-gate-review.jfOaBi7z/partial_amend_loses_measurement
63-STDOUT:
64-A	research/pde_ledger_v3/directives/claim.md
```

```
$ grep -n -A6 'cleanup=strip' _scratch/hook_fix/r0_codex_evidence/message_cleanup_strip.stdout
37:STDIN: {"tool_input": {"command": "git commit --cleanup=strip -F message.txt"}}
38-STDOUT:
39-STDERR:
40-EXIT: 0
41:COMMAND: /bin/bash -c 'git commit --cleanup=strip -F message.txt'
42-CWD: /tmp/measurements-gate-review.jfOaBi7z/message_cleanup_strip
43-STDOUT:
44-[master b654641] A claim
45- 1 file changed, 1 insertion(+)
46- create mode 100644 research/pde_ledger_v3/directives/claim.md
47-STDERR:
```

## Paths are read without -z, so git quotes a non-ASCII name (Codex 3)
```
$ grep -n 'name-only' _scratch/hook_fix/r0_require_measurements_commit_msg.sh
42:if ! staged=$(git diff --cached --name-only 2>&1); then
```

```
$ grep -n 'caf' _scratch/hook_fix/r0_codex_evidence/unicode_basename.stdout
27:WRITE /tmp/measurements-gate-review.jfOaBi7z/unicode_basename/research/pde_ledger_v3/directives/café.md: 'A claim about an artifact.\n'
28:COMMAND: git add -- 'research/pde_ledger_v3/directives/café.md'
45: create mode 100644 "research/pde_ledger_v3/directives/caf\303\251.md"
58:A	"research/pde_ledger_v3/directives/caf\303\251.md"
```

## Only commands with a commit token are seen; the git gate is commit-msg only (Codex 4)
```
$ grep -n 'cherry-pick' _scratch/hook_fix/r0_codex_evidence/cherry_pick_edit.stdout
52:STDIN: {"tool_input": {"command": "GIT_EDITOR=/tmp/measurements-gate-review.jfOaBi7z/cherry_pick_message_editor.sh git cherry-pick -e source"}}
56:COMMAND: /bin/bash -c 'GIT_EDITOR=/tmp/measurements-gate-review.jfOaBi7z/cherry_pick_message_editor.sh git cherry-pick -e source'
```

## The parser over-refuses (Codex 5, 6, 7; Grok 3)
```
$ grep -n 'texts.append(c)\|ARG_LONG = \|name in ARG_LONG\|hookspath\|GIT_CONFIG' _scratch/hook_fix/r0_require_measurements.sh
9:#     or GIT_CONFIG_* in the environment), or
81:ARG_LONG = {"--message", "--file", "--reuse-message", "--reedit-message",
95:            if name in ARG_LONG and "=" not in t:
122:                if any("hookspath" in b.lower() for b in before):
124:                if any(b.startswith("GIT_CONFIG") for b in before):
125:                    why = why or "GIT_CONFIG_* in the environment"
131:    texts.append(c)                         # a misread heredoc must not hide a command
```

```
$ grep -n -B1 -A1 'BLOCKED' _scratch/hook_fix/r0_codex_evidence/false_positive_data.stdout
32-STDERR:
33:BLOCKED — require_measurements: this commit turns git's hooks off (-n).
34-The measurements gate is git's commit-msg hook; a commit that skips it is ungated.
--
51-STDERR:
52:BLOCKED — require_measurements: this commit turns git's hooks off (-n).
53-The measurements gate is git's commit-msg hook; a commit that skips it is ungated.
--
70-STDERR:
71:BLOCKED — require_measurements: this command may run `git commit`, and it could not
72-be tokenised (an unbalanced quote outside any heredoc), so the hook cannot tell
```

```
$ grep -n -B1 -A1 'BLOCKED' _scratch/hook_fix/r0_codex_evidence/abbreviated_message_argument.stdout _scratch/hook_fix/r0_codex_evidence/verify_reenabled.stdout
_scratch/hook_fix/r0_codex_evidence/abbreviated_message_argument.stdout-44-STDERR:
_scratch/hook_fix/r0_codex_evidence/abbreviated_message_argument.stdout:45:BLOCKED — require_measurements: this commit turns git's hooks off (-n).
_scratch/hook_fix/r0_codex_evidence/abbreviated_message_argument.stdout-46-The measurements gate is git's commit-msg hook; a commit that skips it is ungated.
--
_scratch/hook_fix/r0_codex_evidence/verify_reenabled.stdout-44-STDERR:
_scratch/hook_fix/r0_codex_evidence/verify_reenabled.stdout:45:BLOCKED — require_measurements: this commit turns git's hooks off (--no-verify).
_scratch/hook_fix/r0_codex_evidence/verify_reenabled.stdout-46-The measurements gate is git's commit-msg hook; a commit that skips it is ungated.
```

```
$ grep -n -B1 -A1 'BLOCKED' _scratch/hook_fix/r0_codex_evidence/safe_hookspath_override.stdout _scratch/hook_fix/r0_codex_evidence/safe_env_config.stdout
_scratch/hook_fix/r0_codex_evidence/safe_hookspath_override.stdout-44-STDERR:
_scratch/hook_fix/r0_codex_evidence/safe_hookspath_override.stdout:45:BLOCKED — require_measurements: this commit turns git's hooks off (core.hooksPath).
_scratch/hook_fix/r0_codex_evidence/safe_hookspath_override.stdout-46-The measurements gate is git's commit-msg hook; a commit that skips it is ungated.
--
_scratch/hook_fix/r0_codex_evidence/safe_env_config.stdout-44-STDERR:
_scratch/hook_fix/r0_codex_evidence/safe_env_config.stdout:45:BLOCKED — require_measurements: this commit turns git's hooks off (GIT_CONFIG_* in the environment).
_scratch/hook_fix/r0_codex_evidence/safe_env_config.stdout-46-The measurements gate is git's commit-msg hook; a commit that skips it is ungated.
```

## The parser under-refuses (Codex 8; Grok 2, 4)
```
$ grep -n 'line continuations' _scratch/hook_fix/r0_require_measurements.sh
61:    s = s.replace("\\\n", " ")              # line continuations
```

```
$ grep -n -A3 'mit --no-verify' _scratch/hook_fix/r0_codex_evidence/line_continuation_bypass.stdout
36:STDIN: {"tool_input": {"command": "git com\\\nmit --no-verify -m claim"}}
37-STDOUT:
38-STDERR:
39-EXIT: 0
--
41:mit --no-verify -m claim'
42-CWD: /tmp/measurements-gate-review.jfOaBi7z/line_continuation_bypass
43-STDOUT:
44-[master 04d5e6d] claim
```

```
$ grep -n -A3 'core.hooksPath /dev/null &&' _scratch/hook_fix/r0_codex_evidence/prior_config_bypass.stdout
36:STDIN: {"tool_input": {"command": "git config core.hooksPath /dev/null && git commit -m claim"}}
37-STDOUT:
38-STDERR:
39-EXIT: 0
40:COMMAND: /bin/bash -c 'git config core.hooksPath /dev/null && git commit -m claim'
41-CWD: /tmp/measurements-gate-review.jfOaBi7z/prior_config_bypass
42-STDOUT:
43-[master d32b9da] claim
```

```
$ grep -n -B2 -A3 'masked_by_config\|masked_by_env\|slipped' _scratch/hook_fix/r0_grok_evidence/stdout.txt _scratch/hook_fix/r0_grok_evidence/followup_stdout.txt _scratch/hook_fix/r0_grok_evidence/followup2_stdout.txt
_scratch/hook_fix/r0_grok_evidence/stdout.txt-461-pretool_exit=0 stderr=
_scratch/hook_fix/r0_grok_evidence/stdout.txt-462-bash_exit=0
_scratch/hook_fix/r0_grok_evidence/stdout.txt:463:[master 238bcba] masked_by_config| 1 file changed, 1 insertion(+)| create mode 100644 research/pde_ledger_v3/directives/ev.md|
_scratch/hook_fix/r0_grok_evidence/stdout.txt:464:HEAD=238bcba masked_by_config
_scratch/hook_fix/r0_grok_evidence/stdout.txt-465--- I5d export GIT_CONFIG then commit --
_scratch/hook_fix/r0_grok_evidence/stdout.txt-466-pretool_exit=0 stderr=
_scratch/hook_fix/r0_grok_evidence/stdout.txt-467-bash_exit=0
_scratch/hook_fix/r0_grok_evidence/stdout.txt:468:[master c2355f3] masked_by_env| 1 file changed, 1 insertion(+)| create mode 100644 research/pde_ledger_v3/directives/ev.md|
_scratch/hook_fix/r0_grok_evidence/stdout.txt:469:HEAD=c2355f3 masked_by_env
_scratch/hook_fix/r0_grok_evidence/stdout.txt-470-
_scratch/hook_fix/r0_grok_evidence/stdout.txt-471-===== W worktrees =====
_scratch/hook_fix/r0_grok_evidence/stdout.txt-472-main hooks=/tmp/hook_review_meas/wtmain/.git/hooks
--
_scratch/hook_fix/r0_grok_evidence/followup_stdout.txt-39--- pretool on the continued flag --
_scratch/hook_fix/r0_grok_evidence/followup_stdout.txt-40-CMD=git commit --no-\
_scratch/hook_fix/r0_grok_evidence/followup_stdout.txt:41:verify -m slipped
_scratch/hook_fix/r0_grok_evidence/followup_stdout.txt-42-pretool_exit=0
_scratch/hook_fix/r0_grok_evidence/followup_stdout.txt-43-stderr=
_scratch/hook_fix/r0_grok_evidence/followup_stdout.txt-44--- actual commit --
_scratch/hook_fix/r0_grok_evidence/followup_stdout.txt-45-bash_exit=0
_scratch/hook_fix/r0_grok_evidence/followup_stdout.txt:46:[master b6ba0e7] slipped| 1 file changed, 1 insertion(+)| create mode 100644 research/pde_ledger_v3/directives/lc.md|
_scratch/hook_fix/r0_grok_evidence/followup_stdout.txt:47:HEAD=b6ba0e7 slipped
_scratch/hook_fix/r0_grok_evidence/followup_stdout.txt-48-files:
_scratch/hook_fix/r0_grok_evidence/followup_stdout.txt-49-  README
_scratch/hook_fix/r0_grok_evidence/followup_stdout.txt-50-  research/pde_ledger_v3/directives/lc.md
--
_scratch/hook_fix/r0_grok_evidence/followup_stdout.txt-81-===== confirm I5 trees still =====
_scratch/hook_fix/r0_grok_evidence/followup_stdout.txt-82--- repoI log --
_scratch/hook_fix/r0_grok_evidence/followup_stdout.txt:83:c2355f3 masked_by_env
_scratch/hook_fix/r0_grok_evidence/followup_stdout.txt-84-e4d1ed8 init
_scratch/hook_fix/r0_grok_evidence/followup_stdout.txt-85--- 238bcba --
_scratch/hook_fix/r0_grok_evidence/followup_stdout.txt-86-research/pde_ledger_v3/directives/ev.md
```

## The installation check is tied to CLAUDE_PROJECT_DIR (Codex 9; Grok 1)
```
$ grep -n 'proj=\|want=' _scratch/hook_fix/r0_require_measurements.sh
174:proj=${CLAUDE_PROJECT_DIR:-/var/projects/toy_physics}
175:want="$proj/.claude/hooks/require_measurements_commit_msg.sh"
```

## The comment on scripts and aliases (Grok 5)
```
$ grep -n 'git still runs the commit-msg hook' _scratch/hook_fix/r0_require_measurements.sh
19:# seen here; git still runs the commit-msg hook for the first two.
```

```
$ grep -n -B2 -A2 'from script skip\|via alias skip' _scratch/hook_fix/r0_grok_evidence/stdout.txt _scratch/hook_fix/r0_grok_evidence/followup_stdout.txt _scratch/hook_fix/r0_grok_evidence/followup2_stdout.txt
_scratch/hook_fix/r0_grok_evidence/stdout.txt-518-pretool_exit=0 stderr=
_scratch/hook_fix/r0_grok_evidence/stdout.txt-519-bash_exit=0
_scratch/hook_fix/r0_grok_evidence/stdout.txt:520:[master 713a5e7] from script skip| 1 file changed, 1 insertion(+)| create mode 100644 research/pde_ledger_v3/directives/s.md|
_scratch/hook_fix/r0_grok_evidence/stdout.txt:521:HEAD=713a5e7 from script skip
_scratch/hook_fix/r0_grok_evidence/stdout.txt-522--- S3 alias ci=commit, hook installed --
_scratch/hook_fix/r0_grok_evidence/stdout.txt-523-pretool_exit=0 stderr=
--
_scratch/hook_fix/r0_grok_evidence/stdout.txt-528-pretool_exit=0 stderr=
_scratch/hook_fix/r0_grok_evidence/stdout.txt-529-alias_exit=0
_scratch/hook_fix/r0_grok_evidence/stdout.txt:530:[master ca95a9d] via alias skip| 1 file changed, 1 insertion(+)| create mode 100644 research/pde_ledger_v3/directives/s.md|
_scratch/hook_fix/r0_grok_evidence/stdout.txt:531:HEAD=ca95a9d via alias skip
_scratch/hook_fix/r0_grok_evidence/stdout.txt-532--- S5 commit-tree --
_scratch/hook_fix/r0_grok_evidence/stdout.txt-533-commit-tree=14071df5b4a4afe1aa6f1533d00e362203c96e4f
```

