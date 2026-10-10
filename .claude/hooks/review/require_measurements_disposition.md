# The measurements commit gate: review dispositions (orchestrator)

**Artifact.** The control that enforces `CLAUDE.md` E1's grounding duty for directives: a commit that changes a
`research/pde_ledger_v3/directives/*.md` (or a `_legs/` brief) ships `directives/_measurements/<same basename>` or
carries `no-measurements: <reason>`. It is orchestrator-written. It changes a check, so it is reviewed by Codex and
Grok until clear, and both reports come before any commit (`CLAUDE.md` G1, G3, G4; "Other" row of the review table).

**Why it was repaired (2026-10-10, user: "Sure fix it").** The version committed at `8ca68b2a` ran as a Claude Code
`PreToolUse` hook. It read the index before the Bash command ran, and it skipped any command Python's `shlex` could
not tokenise. Directive commits `ede8aa21` and `34a994d5` each staged and committed in one command, so neither was
gated. Each message also held one apostrophe.

## Round 0

**Reviewed version.** `PreToolUse` hook `require_measurements.sh` (sha256 `be3f46bb…`), rewritten to parse the
command for hook-disabling options and to check installation, plus a new git `commit-msg` hook
`require_measurements_commit_msg.sh` (sha256 `c9d2a31a…`) doing the check. Baseline
`_scratch/hook_fix/review_baseline_r0.sha256`; preserved unaccepted at `db1d7033`. Identical prompt
`_scratch/hook_fix/review_prompt_r0.md` (sha256 `74ea34b9…`) to both legs; both reported before I adjudicated.
- **Codex (`gpt-6.1-sol`, xhigh): NEEDS REVISION, 9 findings.** `_scratch/hook_fix/review_r0_codex_final.txt`; full
  report and reproductions copied to `_scratch/hook_fix/r0_codex_evidence/`.
- **Grok (`grok-4.7`): NEEDS REVISION, 5 findings.** `_scratch/hook_fix/review_r0_grok.txt`; scripts and stdout copied
  to `_scratch/hook_fix/r0_grok_evidence/`.

**Verification.** Mechanical lookups in `require_measurements_lookups_r0.md` (generator
`_scratch/hook_fix/gen_lookups_r0.sh`): each finding's cited lines in the frozen r0 scripts, and the leg's own
reproduction output.

| # | Finding | Disposition |
|---|---|---|
| Codex 1 | `--amend` is compared with the amended commit, not its parent. Deleting the measurement and amending leaves a directive-only commit. | **ACCEPT.** `git diff --cached` (r0 L42); the leg's commit `013e616` records only the directive. |
| Codex 2 | The hatch is read from the uncleaned message; `--cleanup=strip` removes a reason on a comment line after the hook accepted it. | **ACCEPT.** r0 L53 greps the message file; the leg's commit `b654641` has no reason. |
| Codex 3 | Paths are read without `-z`; git quotes `café.md`, the anchored pattern misses it. | **ACCEPT.** r0 L42; the leg's output shows the quoted path committed. |
| Codex 4 | `cherry-pick -e` can drop a source commit's reason and create an ungated commit. | **ACCEPT.** The r0 git gate is `commit-msg` only; the leg's cherry-pick committed. |
| Codex 5 | The parser refuses harmless data (`printf` arguments, heredoc text) and a valid `END-MESSAGE` heredoc. | **ACCEPT.** r0 L131 rescans the raw text; the leg's three refusals. |
| Codex 6, Grok 3 | Abbreviated argument-taking long options (`--mes '-n'`) are not consumed, and a later `--verify` is ignored. | **ACCEPT.** r0 L81, L95; the leg's refusals of commands git ran with the hook. |
| Codex 7 | A hook path that keeps the gate, or an unrelated `GIT_CONFIG_*` setting, is refused and labelled "turns hooks off". | **ACCEPT the label; keep the refusal.** The label is wrong. The refusal is kept, worded accurately: the hook cannot see before the command runs what such a setting resolves to, and the project has no use for either. |
| Codex 8, Grok 2, Grok 4 | `com\`+newline+`mit`, `--no-\`+newline+`verify`, and a hook path set by an earlier command in the same call (`git config …`, `export GIT_CONFIG_…`) commit ungated. | **ACCEPT.** r0 L61 turns a continuation into a space; the legs' commits `04d5e6d`, `b6ba0e7`, `d32b9da`, `238bcba`, `c2355f3`. |
| Codex 9, Grok 1 | The installation check is tied to `CLAUDE_PROJECT_DIR`: a linked worktree's inherited hook is refused, and a worktree-specific hook path reached by `git -C` is not seen. | **ACCEPT.** r0 L174–175. |
| Grok 5 | The comment says git still runs `commit-msg` for scripts and aliases; one that passes `--no-verify` does not. | **ACCEPT.** r0 L19; the leg's commits `713a5e7`, `ca95a9d`. |

**Repair (round 1).** Most defects were in the command parser, which I wrote in round 0. Repairing it would mean more
parsing. The repair removes both jobs the parser did:
- **The gate moves to git's `reference-transaction` hook** (`require_measurements_gate.sh`). Git runs it before every
  ref update, and `--no-verify` does not skip it. It judges each commit that no ref reached before the update: its
  own diff against its own parents, with `-z` paths, and its message as git recorded it. That covers amend,
  cherry-pick, rebase, `commit-tree` + `update-ref`, scripts and aliases (Codex 1–4, 8; Grok 2, 4, 5).
- **The `PreToolUse` hook stops parsing.** It refuses a command whose text mentions a hook path or `GIT_CONFIG*`, and
  it refuses commands mentioning `git` or `datalad` while the hook git will run is not the gate (Codex 5–7, 9;
  Grok 1, 3).
- **Rule changes, recorded:** the reason after `no-measurements:` must be non-empty; a counterpart counts only if the
  commit adds or changes it (deleting it does not); a merge is judged by the paths that differ from every parent.

The `commit-msg` hook is removed.

## Round 1

**Reviewed version.** `require_measurements_gate.sh` (sha256 `333263c1…`) and `require_measurements.sh` (sha256
`658cdcd0…`); baseline `_scratch/hook_fix/review_baseline_r1.sha256`; preserved unaccepted at `ce377dde`. Identical
prompt `_scratch/hook_fix/review_prompt_r1.md` (sha256 `09278d8f…`) to fresh Codex and Grok legs; both reported
before I adjudicated.
- **Codex (`gpt-6.1-sol`, xhigh): NEEDS REVISION, 6 findings.** `_scratch/hook_fix/review_r1_codex_final.txt`;
  report and reproductions copied to `_scratch/hook_fix/r1_codex_evidence/`.
- **Grok (`grok-4.7`): NEEDS REVISION, 4 findings,** with a list of what held (ordinary operations, amend,
  cherry-pick, rebase, `am`, `--no-verify`, aliases, scripts, worktrees, `datalad save`, message cleanup).
  `_scratch/hook_fix/review_r1_grok.txt`; scripts and stdout copied to `_scratch/hook_fix/r1_grok_evidence/`.

**Verification.** Mechanical lookups in `require_measurements_lookups_r1.md` (generator
`_scratch/hook_fix/gen_lookups_r1.sh`).

| # | Finding | Disposition |
|---|---|---|
| A (Codex 1, Grok 2) | `rev-list … --not --all` skips any commit some ref reaches, including refs the gate never checks. A commit reached first through the stash (`stash^2`), a tag, a remote-tracking ref (a pull from a clone without the hook) or `refs/heads/git-annex` lands on a branch unchecked. | **ACCEPT.** r1 L68; the legs' branch at `stash^2`, tag-then-branch and pull. Grok's suggested fix (subtract local branches only) would re-judge old history: 89 commits reachable only from the tags `s10-as-built`, `s11-as-built`, `wip-2026-08-05-unreviewed` and `archive/pre-cleanup-2026-10-04` fail the rule (mechanical count, this session), so checking one out would be refused. Repair: a baseline of every ref at installation is grandfathered; local branches (not `git-annex`, `synced/*`) and the moving refs' old values are subtracted; anything else is checked when it first reaches HEAD or a branch. |
| B (Codex 2, Grok 3) | `-M --name-only` keeps only a rename's destination, so moving a directive out of `directives/` is not seen. | **ACCEPT.** r1 L79. Repair: no rename detection; a rename is a deletion plus an addition, and both sides are judged. A paired rename inside `directives/` now needs the reason line (recorded rule change). |
| C (Codex 3, 4; Grok 1) | NUL records are turned into lines and matched in the ambient locale: a newline in a name splits it (a bypass, a fake counterpart, and a false refusal), and an invalid UTF-8 name is not matched. | **ACCEPT.** r1 L86–91. Repair: NUL-delimited records read into arrays, matched with `LC_ALL=C`, counterpart looked up by exact key and confirmed to exist in the commit. |
| D (Codex 5) | The text rule matches the substring `hookspath` anywhere (`printf hookspathology`). | **ACCEPT.** r1 L42. Repair: match the setting name `core.hooksPath` and words beginning `GIT_CONFIG`, at word boundaries. |
| E (Codex 6, Grok 4) | The gate's comment says the PreToolUse hook refuses commands that touch the hook path; it does not see `rm .git/hooks/reference-transaction`. | **ACCEPT.** r1 L32. Repair: the comment states what the PreToolUse hook checks and what it does not see. |

**Repair (round 2).** `require_measurements_gate.sh` with the changes above; a new `install_measurements_gate.sh`
that writes the baseline once (never overwritten) and links the hook; `require_measurements.sh` with the narrowed
text rule and a baseline check. A, C and D are in material I wrote in round 1. If round 2 finds defects in the
material changed here, the next revision changes author (G4).

## Round 2

**Reviewed version.** `require_measurements_gate.sh` (sha256 `0ae77684…`), `install_measurements_gate.sh` (sha256
`7e3bb54e…`) and `require_measurements.sh` (sha256 `fae06f49…`); baseline `_scratch/hook_fix/review_baseline_r2.sha256`;
frozen as `_scratch/hook_fix/r2_*.sh`. Identical prompt `_scratch/hook_fix/review_prompt_r2.md` (sha256 `035d0a0d…`)
to fresh Codex and Grok legs; both reported before I adjudicated.
- **Codex (`gpt-6.1-sol`, xhigh): NEEDS REVISION, 9 findings.** `_scratch/hook_fix/review_r2_codex_final.txt`;
  scripts and stdout copied to `_scratch/hook_fix/r2_codex_evidence/`.
- **Grok (`grok-4.7`): NEEDS REVISION, 4 findings,** with a list of what held. `_scratch/hook_fix/review_r2_grok.txt`;
  scripts and stdout copied to `_scratch/hook_fix/r2_grok_evidence/`.

**Verification.** Mechanical lookups in `require_measurements_lookups_r2.md` (generator
`_scratch/hook_fix/gen_lookups_r2.sh`): each finding's cited lines in the frozen r2 scripts, and the legs' own
reproduction output.

**Prompt correction.** The r2 prompt named `8ca68b2a` as the version before the repair. The last committed version
before the repair is `293e7e53`, which already gates `_legs/` briefs (lookups, "Prompt reference commit"). Codex
treated the briefs as intended; the semantics are unchanged. Later prompts name `293e7e53`.

| # | Finding | Disposition |
|---|---|---|
| A (Codex 4, Grok 1) | Replace refs and a grafts file make the gate's git commands judge a substitute commit, not the one the ref stores. | **ACCEPT.** r2 L74–133 run `rev-list`, `diff-tree`, `log` and `cat-file` with replacement and grafts in force; the legs' `update-ref` and `branch` of a directive commit exit 0. Present since r1. |
| B (Codex 8, Grok 2) | `2>&1` makes git's advice text into commit ids; a grafts file refuses an ordinary commit. | **ACCEPT.** r2 L94, L100; the legs' `BLOCKED … could not read hint:` and exit 128. Present since r1. |
| C (Codex 1, Grok 3) | `git symbolic-ref` points a branch or HEAD at an unjudged commit without running the hook; the comment says git runs it before every reference update. | **ACCEPT.** r2 L20; the legs' `symbolic-ref` exits 0 and the branch resolves to the directive commit. Present since r1. |
| D (Codex 2) | `refs/heads/git-annex` and `refs/heads/synced/*` are not judged, and checkout makes them HEAD. | **ACCEPT.** r2 L68; the leg's `update-ref refs/heads/synced/unsafe` then `checkout`, both exit 0. The `synced/*` exclusion is new in r2. |
| E (Codex 3) | A worktree-qualified HEAD (`worktrees/<name>/HEAD`) is not judged. | **ACCEPT.** r2 L69–70; the leg's `update-ref worktrees/linked_head/HEAD` exits 0. Present since r1. |
| F (Codex 5) | The hatch is read from `git log` output, which configuration can add to (signature output) or re-encode. | **ACCEPT.** r2 L121; the leg's `gpg: no-measurements:` line admitted a commit whose message has no reason, and `i18n.logOutputEncoding` refused a recorded reason. Present since r1. |
| G (Codex 6) | The installed hook links into the checked-out tree, so checking out history from before the gate leaves the link dangling and git skips it; the `PreToolUse` script is then also absent (exit 127). | **ACCEPT.** Installer L16–19, L41; the leg's commit after a historical checkout exits 0 with a directive. The link design is from r1. |
| H (Codex 7, Grok 4) | The baseline omits detached HEADs, so returning to a pre-install detached commit is refused. | **ACCEPT.** Installer L27–29; the legs' checkout back exits 128. New in r2. |
| I (Codex 9) | A dangling symlink at the baseline path is overwritten, contrary to the comment. | **ACCEPT.** Installer L8, L23, L30; the leg's `ls -l` before and after. New in r2. |

**Repair (round 3): change of author (G4).** D, H and I are defects in material written in round 2, and I have
folded three times. Codex authors round 3 from a brief that states what must be true and points at the round-2
reports. The brief is reviewed by two legs before Codex starts (G2). Codex's version is reviewed by a fresh Claude
agent and Grok (G1) until clear.
