#!/bin/bash
# git commit-msg hook — CLAUDE.md E1, orchestrator half. This is the gate itself.
#
# A document under directives/ that makes a claim about an artifact must ship the
# commands that produced the claim. This refuses the commit when a directives/*.md
# (or a _legs/ brief) is in the commit without its _measurements/ counterpart.
#
# Installed as .git/hooks/commit-msg, a symlink to this file:
#     ln -s ../../.claude/hooks/require_measurements_commit_msg.sh .git/hooks/commit-msg
# The PreToolUse hook require_measurements.sh refuses a commit while it is not
# installed, or when the command turns hooks off.
#
# Why it runs inside git. Measured 2026-10-10: the first version ran before the
# Bash command and read the index as it stood then. A command that staged a
# directive and committed it in one line (`git add D && git commit ...`) was never
# gated, and neither was one whose heredoc message held an odd apostrophe (the
# command did not tokenise, and the hook skipped it). Two directive commits went
# through that way, ede8aa21 and 34a994d5. Git runs this hook with the final
# index and the final message, however the command is phrased.
#
# It is not tamper-proof against the agent; it is tamper-EVIDENT. Evading it, or
# using the escape hatch, is visible in the transcript and in the commit message.
#
# Escape hatch, deliberately permanent and public: put
#     no-measurements: <reason>
# in the commit message. A directive that genuinely asserts nothing about an
# artifact is legitimate — saying so in the record is the price. The reason must
# be on the same line and must not be empty.

set -uo pipefail

msg_file=${1:-}

# Briefs are gated too. Measured 2026-08-12: a fix brief demanded an end-state the
# declared population provably could not deliver, carried no measurement, and cost a
# review round. A brief that asserts what is achievable is making a claim about an
# artifact, and rule 2 binds it the same way.
DIR_RE='^research/pde_ledger_v3/directives/([^/]+\.md|_legs/[^/]*brief[^/]*\.md)$'

# The commit's contents. Under `git commit -a` or `git commit <paths>`, git points
# GIT_INDEX_FILE at the index it will commit, and `git diff --cached` reads it.
if ! staged=$(git diff --cached --name-only 2>&1); then
  printf 'BLOCKED — require_measurements: could not list the files in this commit:\n%s\n' "$staged" >&2
  exit 1
fi
[ -z "$staged" ] && exit 0

docs=$(printf '%s\n' "$staged" | grep -E "$DIR_RE" || true)
[ -z "$docs" ] && exit 0

# Escape hatch: an explicit, recorded reason in the commit message.
if [ -n "$msg_file" ] && [ -r "$msg_file" ] \
   && grep -qiE 'no-measurements:[[:space:]]*[^[:space:]]' "$msg_file"; then
  exit 0
fi

missing=""
while IFS= read -r doc; do
  [ -z "$doc" ] && continue
  base=$(basename "$doc")
  m="research/pde_ledger_v3/directives/_measurements/$base"
  printf '%s\n' "$staged" | grep -qxF "$m" || missing="$missing  $doc  ->  $m"$'\n'
done <<< "$docs"

[ -z "$missing" ] && exit 0

cat >&2 <<EOF
BLOCKED — CLAUDE.md E1 (orchestrator half).

These directives are in this commit with no matching _measurements/ file beside them:

$missing
A claim about an artifact carries the command that produced it. Run the commands,
write them and their LITERAL output to the path above, stage it, and commit again.
Regenerate the file from the commands; do not transcribe.

Measured 2026-08-12, the cost of skipping this: four export-chain designs and eight
review legs died on a question one \`len(LEDGER)\` answered.

If this document genuinely asserts nothing about any artifact, say so in the commit
message with:  no-measurements: <reason>
EOF
exit 1
