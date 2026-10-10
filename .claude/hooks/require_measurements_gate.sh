#!/bin/bash
# git reference-transaction hook — CLAUDE.md E1, orchestrator half. This is the gate itself.
#
# A document under directives/ that makes a claim about an artifact must ship the
# commands that produced the claim. This refuses to let a branch (or a detached HEAD)
# move onto a commit that adds, changes or deletes a directives/*.md, or a _legs/
# brief, unless that commit adds or changes its _measurements/ counterpart (deleting
# the counterpart does not count), or its message carries the reason below.
#
# Installed by .claude/hooks/install_measurements_gate.sh, which links this file as
# .git/hooks/reference-transaction (linked worktrees share it) and records the
# baseline described below. The PreToolUse hook require_measurements.sh refuses git
# and datalad commands while either is missing.
#
# Why here. Measured 2026-10-10: the first version ran before the Bash command and
# read the index as it stood then, so `git add D && git commit ...` was never gated,
# nor was a command shlex could not tokenise (ede8aa21 and 34a994d5 went through).
# A commit-msg hook then fixed those but missed --amend (it compared with the
# amended commit, not its parent), cherry-pick, --no-verify, and a reason that
# message cleanup removed. Git runs this hook before every reference update, and
# --no-verify does not skip it. So it judges each commit itself: its own diff
# against its own parents, and its message as recorded.
#
# Which commits. When HEAD or a branch moves, every commit it would then reach is
# checked, except commits reached by: the baseline (every ref that existed when the
# gate was installed — history from before the gate, including tags whose old
# commits would fail, is not re-judged), any local branch other than
# refs/heads/git-annex and refs/heads/synced/*, and the moving refs' old values.
# Those were checked when they reached a branch, or predate the gate. So a commit
# that arrives later through a tag, a remote-tracking ref, the stash or the
# git-annex branch is checked when it first lands on HEAD or a branch.
# Updates of other refs (tags, remotes, the stash, refs/heads/git-annex) are not
# checked themselves.
#
# It is not tamper-proof against the agent; it is tamper-EVIDENT. Removing or
# replacing the hook, editing the baseline, or pointing git at another hook
# directory switches it off, and that is visible in the transcript. The PreToolUse
# hook refuses command text that names core.hooksPath or a GIT_CONFIG* variable, and
# git or datalad commands while the gate or its baseline is missing; it does not see
# a hook removed or replaced inside the same command.
#
# Escape hatch, deliberately permanent and public: put
#     no-measurements: <reason>
# in the commit message. A directive that genuinely asserts nothing about an
# artifact is legitimate — saying so in the record is the price. The reason must
# be on the same line and must not be empty.

set -uo pipefail
export LC_ALL=C                      # paths are bytes; match them as bytes

if [ "${1:-}" != prepared ]; then
  cat >/dev/null
  exit 0
fi

# Briefs are gated too. Measured 2026-08-12: a fix brief demanded an end-state the
# declared population provably could not deliver, carried no measurement, and cost a
# review round. A brief that asserts what is achievable is making a claim about an
# artifact, and rule 2 binds it the same way.
DIR_RE='^research/pde_ledger_v3/directives/([^/]+\.md|_legs/[^/]*brief[^/]*\.md)$'
MEAS=research/pde_ledger_v3/directives/_measurements

refuse() { printf 'BLOCKED — require_measurements: %s\n' "$1" >&2; exit 1; }

tips=() olds=()
while read -r old new ref; do
  case "$ref" in
    refs/heads/git-annex|refs/heads/synced/*) continue ;;
    HEAD|refs/heads/*) ;;
    *) continue ;;
  esac
  case "$old" in *[!0]*) olds+=("$old") ;; esac
  case "$new" in *[!0]*) ;; *) continue ;; esac            # all zeros: a deletion
  [ "$(git cat-file -t "$new" 2>/dev/null)" = commit ] || continue
  tips+=("$new")
done
[ ${#tips[@]} -eq 0 ] && exit 0

baseline="$(git rev-parse --path-format=absolute --git-common-dir 2>/dev/null)/measurements-gate-baseline"
[ -r "$baseline" ] || refuse "no baseline at $baseline; run: bash .claude/hooks/install_measurements_gate.sh"

tmp=$(mktemp) || refuse "could not create a temporary file"
trap 'rm -f "$tmp"' EXIT

# Commits the moving refs would reach that nothing checked or grandfathered reaches.
# --ignore-missing: a baseline commit that gc has since removed is simply skipped.
{
  printf '%s\n' "${tips[@]}"
  sed 's/^/^/' "$baseline"
  [ ${#olds[@]} -gt 0 ] && printf '^%s\n' "${olds[@]}"
  git for-each-ref --format='%(refname) ^%(objectname)' refs/heads \
    | grep -v '^refs/heads/git-annex \|^refs/heads/synced/' | sed 's/.* //'
} > "$tmp"
new_commits=$(git rev-list --stdin --ignore-missing < "$tmp" 2>&1) \
  || refuse "could not list the new commits: $new_commits"
[ -z "$new_commits" ] && exit 0

report=""
for c in $new_commits; do
  parents=$(git rev-list --no-walk --parents "$c" 2>&1) || refuse "could not read $c: $parents"
  set -- $parents
  if [ $# -le 2 ]; then
    # Root or ordinary commit: its diff against its parent (or the empty tree).
    # No rename detection: a directive renamed away is a deletion of the directive.
    how=(-r --no-renames --root)
  else
    # Merge: only paths that differ from every parent, the merge's own change.
    how=(-r --no-renames -c)
  fi
  git diff-tree "${how[@]}" --no-commit-id --name-only -z "$c" > "$tmp" || refuse "could not diff $c"
  changed=() docs=()
  mapfile -d '' -t changed < "$tmp"
  declare -A in_commit=()
  for p in "${changed[@]}"; do
    in_commit["$p"]=1
    [[ $p =~ $DIR_RE ]] && docs+=("$p")
  done
  [ ${#docs[@]} -eq 0 ] && { unset in_commit; continue; }

  # Escape hatch: an explicit, recorded reason in the message as git stored it.
  msg=$(git log -1 --format=%B "$c") || refuse "could not read the message of $c"
  if grep -qiE 'no-measurements:[[:space:]]*[^[:space:]]' <<< "$msg"; then unset in_commit; continue; fi

  missing=""
  for doc in "${docs[@]}"; do
    m="$MEAS/${doc##*/}"
    # The counterpart counts if this commit changed it and it exists afterwards.
    if [ -z "${in_commit[$m]:-}" ] || ! git cat-file -e "$c:$m" 2>/dev/null; then
      missing+="$(printf '    %q  ->  %q' "$doc" "$m")"$'\n'
    fi
  done
  unset in_commit
  [ -n "$missing" ] && report+="  $(git log -1 --format='%h %s' "$c")"$'\n'"$missing"
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
