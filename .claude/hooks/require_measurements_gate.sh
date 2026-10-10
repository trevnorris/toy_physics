#!/bin/bash
# git reference-transaction hook — CLAUDE.md E1, orchestrator half. This is the gate itself.
#
# A document under directives/ that makes a claim about an artifact must ship the
# commands that produced the claim. This refuses to let a branch (or a detached HEAD)
# move to a new commit that adds, changes or deletes a directives/*.md, or a _legs/
# brief, unless the same commit adds or changes its _measurements/ counterpart
# (deleting the counterpart does not count), or its message carries the reason below.
#
# Installed as .git/hooks/reference-transaction, a symlink to this file:
#     ln -s ../../.claude/hooks/require_measurements_gate.sh .git/hooks/reference-transaction
# Linked worktrees share it. The PreToolUse hook require_measurements.sh refuses git
# and datalad commands while it is not installed.
#
# Why here. Measured 2026-10-10: the first version ran before the Bash command and
# read the index as it stood then, so `git add D && git commit ...` was never gated,
# nor was a command shlex could not tokenise (ede8aa21 and 34a994d5 went through).
# A commit-msg hook then fixed those but missed --amend (it compared with the
# amended commit, not its parent), cherry-pick, --no-verify, and a reason that
# message cleanup removed. Git runs this hook before every reference update, and
# --no-verify does not skip it. So it judges each new commit itself: its own diff
# against its own parents, and its message as recorded.
#
# Which commits. Only commits that no ref reached before this update: a new commit,
# an amended one, a cherry-pick, a rebase step, a commit-tree made reachable. Moving
# a branch to a commit some ref already reaches (checkout, reset, fast-forward) checks
# nothing. Replaying an old commit that has neither counterpart nor reason is refused.
# Remote-tracking refs, tags, the stash and refs/heads/git-annex are not checked.
#
# It is not tamper-proof against the agent; it is tamper-EVIDENT. Removing the hook,
# or pointing core.hooksPath elsewhere, switches it off, and that is visible in the
# transcript. The PreToolUse hook refuses commands that touch the hook path.
#
# Escape hatch, deliberately permanent and public: put
#     no-measurements: <reason>
# in the commit message. A directive that genuinely asserts nothing about an
# artifact is legitimate — saying so in the record is the price. The reason must
# be on the same line and must not be empty.

set -uo pipefail

if [ "${1:-}" != prepared ]; then
  cat >/dev/null
  exit 0
fi

# Briefs are gated too. Measured 2026-08-12: a fix brief demanded an end-state the
# declared population provably could not deliver, carried no measurement, and cost a
# review round. A brief that asserts what is achievable is making a claim about an
# artifact, and rule 2 binds it the same way.
DIR_RE='^research/pde_ledger_v3/directives/([^/]+\.md|_legs/[^/]*brief[^/]*\.md)$'

refuse() { printf 'BLOCKED — require_measurements: %s\n' "$1" >&2; exit 1; }

tips=()
while read -r old new ref; do
  case "$ref" in
    refs/heads/git-annex) continue ;;
    HEAD|refs/heads/*) ;;
    *) continue ;;
  esac
  case "$new" in *[!0]*) ;; *) continue ;; esac            # all zeros: a deletion
  [ "$(git cat-file -t "$new" 2>/dev/null)" = commit ] || continue
  tips+=("$new")
done
[ ${#tips[@]} -eq 0 ] && exit 0

new_commits=$(git rev-list "${tips[@]}" --not --all 2>&1) \
  || refuse "could not list the new commits: $new_commits"
[ -z "$new_commits" ] && exit 0

report=""
for c in $new_commits; do
  parents=$(git rev-list --no-walk --parents "$c" 2>&1) || refuse "could not read $c: $parents"
  set -- $parents
  if [ $# -le 2 ]; then
    # Root or ordinary commit: its diff against its parent (or the empty tree).
    # Renames are detected, so a renamed directive is gated under its new name.
    how=(-r -M --root)
  else
    # Merge: only paths that differ from every parent, the merge's own change.
    how=(-r -c)
  fi
  # Every changed path is a candidate directive; a counterpart counts only if the
  # commit adds or changes it — deleting the measurements file does not ship it.
  files=$(git diff-tree "${how[@]}" --no-commit-id --name-only -z "$c" | tr '\0' '\n') \
    || refuse "could not diff $c"
  kept=$(git diff-tree "${how[@]}" --diff-filter=d --no-commit-id --name-only -z "$c" | tr '\0' '\n') \
    || refuse "could not diff $c"

  docs=$(grep -E "$DIR_RE" <<< "$files" || true)
  [ -z "$docs" ] && continue

  # Escape hatch: an explicit, recorded reason in the message as git stored it.
  # (Here-strings, not pipes: under pipefail an early-exiting grep -q can fail the pipe.)
  msg=$(git log -1 --format=%B "$c") || refuse "could not read the message of $c"
  grep -qiE 'no-measurements:[[:space:]]*[^[:space:]]' <<< "$msg" && continue

  missing=""
  while IFS= read -r doc; do
    [ -z "$doc" ] && continue
    m="research/pde_ledger_v3/directives/_measurements/$(basename "$doc")"
    grep -qxF "$m" <<< "$kept" || missing="$missing    $doc  ->  $m"$'\n'
  done <<< "$docs"
  [ -n "$missing" ] && report="$report  $(git log -1 --format='%h %s' "$c")"$'\n'"$missing"
done

[ -z "$report" ] && exit 0

cat >&2 <<EOF
BLOCKED — CLAUDE.md E1 (orchestrator half).

These commits change directives with no matching _measurements/ file in the same commit:

$report
A claim about an artifact carries the command that produced it. Run the commands,
write them and their LITERAL output to the path above, stage it, and commit again.
Regenerate the file from the commands; do not transcribe.

Measured 2026-08-12, the cost of skipping this: four export-chain designs and eight
review legs died on a question one \`len(LEDGER)\` answered.

If this document genuinely asserts nothing about any artifact, say so in the commit
message with:  no-measurements: <reason>
EOF
exit 1
