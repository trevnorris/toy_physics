#!/bin/bash
# PreToolUse(Bash) — keeps the measurements gate installed. CLAUDE.md E1, orchestrator half.
#
# The gate itself is git's reference-transaction hook, require_measurements_gate.sh.
# Git runs it for every commit however the command is written, so this hook does
# not try to recognise commits. It refuses a Bash command when either holds:
#   - the command text contains the setting name core.hooksPath (any case), or a
#     word beginning GIT_CONFIG. Either can point git at another hook directory and
#     switch the gate off, and this hook cannot see before the command runs what
#     such a setting will resolve to. The project has no use for either, so it
#     refuses them all, including settings that would leave the gate in place.
#   - the command text contains the word git or datalad, and either the hook git
#     will run in the project's repository is not the gate (missing, replaced, not
#     executable, or masked by a configured hook path) or the gate's baseline file
#     is missing.
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
# (linked worktrees share its hook). A hook removed or replaced inside the same
# command that then commits, or a repository configured beforehand by other means,
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

if grep -qiE '(^|[^[:alnum:]_])core\.hookspath([^[:alnum:]_]|$)' <<< "$cmd" \
   || grep -qE '(^|[^[:alnum:]_])GIT_CONFIG' <<< "$cmd"; then
  cat >&2 <<'EOF'
BLOCKED — require_measurements: this command names core.hooksPath or a GIT_CONFIG*
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
[ -n "$common" ] && [ -r "$common/measurements-gate-baseline" ] || ok=0
[ $ok = 1 ] && exit 0

cat >&2 <<EOF
BLOCKED — require_measurements: the measurements gate is not installed.

The gate runs inside git, so a commit made without it is ungated. Expected
  ${hooks:-<git hooks directory>}/reference-transaction  ->  .claude/hooks/$gate
  ${common:-<git common directory>}/measurements-gate-baseline
found
  hook: ${have:-nothing}
  baseline: $([ -n "$common" ] && [ -r "$common/measurements-gate-baseline" ] && echo present || echo missing)
Install it from the repository with:
  bash .claude/hooks/install_measurements_gate.sh
If git reads hooks from somewhere else (a configured hook path), remove that setting.
EOF
exit 2
