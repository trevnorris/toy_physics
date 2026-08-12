#!/bin/bash
# PreToolUse(Bash) gate — CLAUDE.md rule 2, orchestrator half.
#
# A document under directives/ that makes a claim about an artifact must ship the
# commands that produced the claim. This refuses the commit when a directives/*.md
# is staged without its _measurements/ counterpart.
#
# It is not tamper-proof against the agent; it is tamper-EVIDENT. Evading it, or
# using the escape hatch, is visible in the transcript and in the commit message.
#
# Escape hatch, deliberately permanent and public: put
#     no-measurements: <reason>
# in the commit message. A directive that genuinely asserts nothing about an
# artifact is legitimate — saying so in the record is the price.

set -uo pipefail

payload=$(cat)

# Extract the command, and decide whether it actually INVOKES `git commit`.
# Substring matching is wrong: a command that merely mentions the string inside a
# quoted argument is not a commit. Measured 2026-08-12 — the first version of this
# hook blocked its own test harness that way. shlex tokenises respecting quotes,
# so only adjacent bare tokens count.
read -r cmd is_commit <<<"$(printf '%s' "$payload" | python3 -c '
import json, shlex, sys
try:
    c = json.load(sys.stdin).get("tool_input", {}).get("command", "")
except Exception:
    c = ""
hit = False
try:
    t = shlex.split(c, comments=False)
    hit = any(t[i] in ("git", "/usr/bin/git") and t[i+1] == "commit" for i in range(len(t)-1))
except ValueError:
    hit = False          # unparseable => do not gate, do not crash the tool call
print("x", "1" if hit else "0")
')"
[ "${is_commit:-0}" = "1" ] || exit 0
cmd=$(printf '%s' "$payload" | python3 -c 'import json,sys
try: print(json.load(sys.stdin).get("tool_input",{}).get("command",""))
except Exception: print("")' 2>/dev/null)

cd "${CLAUDE_PROJECT_DIR:-/var/projects/toy_physics}" || exit 0

DIR_RE='^research/pde_ledger_v3/directives/[^/]+\.md$'
staged=$(git diff --cached --name-only 2>/dev/null)
[ -z "$staged" ] && exit 0

docs=$(printf '%s\n' "$staged" | grep -E "$DIR_RE" || true)
[ -z "$docs" ] && exit 0

# Escape hatch: an explicit, recorded reason in the commit message.
# The token is searched in the whole command AND in any -F message file. `sed` was
# used here first and was wrong: it is line-wise, so a multi-line -m message only
# ever yielded its first line and the token on line 3 was invisible. Caught by the
# test matrix, not by reading.
printf '%s' "$cmd" | grep -qi 'no-measurements:' && exit 0
for f in $(printf '%s' "$cmd" | tr ' ' '\n' | grep -v '^-' | grep -E '\.(txt|md)$'); do
  [ -r "$f" ] && grep -qi 'no-measurements:' "$f" && exit 0
done

missing=""
while IFS= read -r doc; do
  [ -z "$doc" ] && continue
  base=$(basename "$doc")
  m="research/pde_ledger_v3/directives/_measurements/$base"
  printf '%s\n' "$staged" | grep -qxF "$m" || missing="$missing  $doc  ->  $m"$'\n'
done <<< "$docs"

[ -z "$missing" ] && exit 0

cat >&2 <<EOF
BLOCKED — CLAUDE.md rule 2 (orchestrator half).

These directives are staged with no matching _measurements/ file staged beside them:

$missing
A claim about an artifact carries the command that produced it. Run the commands,
write them and their LITERAL output to the path above, stage it, and commit again.
Regenerate the file from the commands; do not transcribe.

Measured 2026-08-12, the cost of skipping this: four export-chain designs and eight
review legs died on a question one \`len(LEDGER)\` answered.

If this document genuinely asserts nothing about any artifact, say so in the commit
message with:  no-measurements: <reason>
EOF
exit 2
