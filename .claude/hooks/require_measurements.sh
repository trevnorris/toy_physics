#!/bin/bash
# PreToolUse(Bash) — keeps the measurements gate installed. CLAUDE.md E1, orchestrator half.
#
# The gate itself is git's reference-transaction hook, require_measurements_gate.sh.
# Git runs it for every commit however the command is written, so this hook does
# not try to recognise commits. It refuses a Bash command when either holds:
#   - the command text mentions core.hooksPath or a GIT_CONFIG* variable, anywhere,
#     in any form. Either can point git at another hook directory and switch the
#     gate off, and this hook cannot see what such a setting will resolve to before
#     the command runs. The project has no use for either, so it refuses them all,
#     including settings that would leave the gate in place.
#   - the command text mentions git or datalad, and the hook git will run in the
#     project's repository is not the gate (missing, replaced, not executable, or
#     masked by a configured core.hooksPath).
#
# Measured 2026-10-10: the first version parsed the command to decide whether it
# committed, read the index before the command ran, and skipped what it could not
# tokenise; directive commits went through ungated (ede8aa21, 34a994d5). Its
# replacement parsed shell and git options to find hook-disabling flags, and the
# review found that parsing both over- and under-refused. Neither job is needed
# now that the gate runs inside git.
#
# It is not tamper-proof against the agent; it is tamper-EVIDENT. The installation
# check covers the project's repository as git resolves it from CLAUDE_PROJECT_DIR
# (linked worktrees share its hook). Removing the hook inside the same command that
# then commits, or committing in a repository configured beforehand by other means,
# is not seen here.

set -uo pipefail

payload=$(cat)
cmd=$(printf '%s' "$payload" | python3 -c '
import json, sys
raw = sys.stdin.read()
try:
    c = json.loads(raw).get("tool_input", {}).get("command", "")
    print(c if isinstance(c, str) else raw, end="")
except Exception:
    print(raw, end="")             # unreadable payload: check its whole text
' 2>/dev/null) || cmd=$payload

if grep -qiE 'hookspath|GIT_CONFIG' <<< "$cmd"; then
  cat >&2 <<'EOF'
BLOCKED — require_measurements: this command mentions core.hooksPath or a GIT_CONFIG*
variable. Either can point git at another hook directory and switch off the
measurements gate (git's reference-transaction hook), so neither is allowed from
here. Restate the command without it.
EOF
  exit 2
fi

grep -qE '(^|[^[:alnum:]_.-])(git|datalad)([^[:alnum:]_]|$)' <<< "$cmd" || exit 0

proj=${CLAUDE_PROJECT_DIR:-/var/projects/toy_physics}
gate=require_measurements_gate.sh
hooks=$(git -C "$proj" rev-parse --path-format=absolute --git-path hooks 2>/dev/null) || hooks=""
common=$(git -C "$proj" rev-parse --path-format=absolute --git-common-dir 2>/dev/null) || common=""
have=""
[ -n "$hooks" ] && have=$(readlink -f "$hooks/reference-transaction" 2>/dev/null)
ok=0
for want in "$proj/.claude/hooks/$gate" "${common%/.git}/.claude/hooks/$gate"; do
  [ -n "$have" ] && [ -x "$have" ] && [ "$have" = "$(readlink -f "$want" 2>/dev/null)" ] && ok=1
done
[ $ok = 1 ] && exit 0

cat >&2 <<EOF
BLOCKED — require_measurements: git's reference-transaction hook is not the measurements gate.

The gate runs inside git, so a commit made without it is ungated. Expected
  ${hooks:-<git hooks directory>}/reference-transaction  ->  .claude/hooks/$gate
found
  ${have:-nothing}
Install it from the main checkout's root with:
  ln -s ../../.claude/hooks/$gate .git/hooks/reference-transaction
If git reads hooks from somewhere else (a configured hook path), remove that setting.
EOF
exit 2
