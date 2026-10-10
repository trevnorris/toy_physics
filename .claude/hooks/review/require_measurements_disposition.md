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
