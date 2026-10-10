#!/bin/bash
# PreToolUse(Bash) — keeps the measurements gate in force. CLAUDE.md E1, orchestrator half.
#
# The gate itself is git's commit-msg hook, require_measurements_commit_msg.sh: it
# sees the commit's final contents and message. This hook runs before any Bash
# command that may commit, and refuses it when
#   - that git hook is not installed where git will run it, or
#   - the command turns hooks off (--no-verify or -n on commit, core.hooksPath,
#     or GIT_CONFIG_* in the environment), or
#   - the command cannot be tokenised, so neither can be ruled out.
#
# Measured 2026-10-10: the previous version read the index before the command
# ran, and skipped any command shlex could not tokenise. A directive staged and
# committed in one command, or committed with an odd apostrophe in a heredoc
# message, went through ungated (ede8aa21, 34a994d5).
#
# It is not tamper-proof against the agent; it is tamper-EVIDENT. A commit made
# from inside a script, through an alias, or with plumbing (commit-tree) is not
# seen here; git still runs the commit-msg hook for the first two.

set -uo pipefail

payload=$(cat)

# Decide whether the command may invoke `git commit`, and whether it turns hooks
# off. Substring matching is wrong: a command that merely mentions the string
# inside a quoted argument is not a commit (measured 2026-08-12 — the first version
# of this hook blocked its own test harness that way). So the command is tokenised
# respecting quotes, heredoc bodies are removed first (they are data, and an
# apostrophe in one breaks tokenising), and only bare tokens count.
verdict=$(printf '%s' "$payload" | python3 -c '
import json, os, re, shlex, sys

raw = sys.stdin.read()
try:
    c = json.loads(raw).get("tool_input", {}).get("command", "")
except Exception:
    # Not a payload we can read. Refuse only if it could be a commit.
    print("unparsed" if "commit" in raw else "pass")
    sys.exit(0)
if not isinstance(c, str):
    c = ""

HEREDOC = re.compile(r"(?<!<)<<(?!<)(-?)\s*\\?([\x27\x22]?)([A-Za-z_][A-Za-z0-9_]*)\2")

def strip_heredocs(s):
    out, pending = [], []
    for line in s.split("\n"):
        if pending:
            delim, dash = pending[0]
            if (line.lstrip("\t") if dash else line) == delim:
                pending.pop(0)
            continue
        out.append(line)
        pending += [(m.group(3), m.group(1) == "-") for m in HEREDOC.finditer(line)]
    return "\n".join(out)

OPS = set("();<>|&\n")

def simple_commands(s):
    s = s.replace("\\\n", " ")              # line continuations
    lex = shlex.shlex(s, posix=True, punctuation_chars="();<>|&\n")
    lex.whitespace = " \t\r"
    lex.whitespace_split = True
    cmds, cur = [], []
    for t in lex:                           # raises ValueError if unbalanced
        if t and set(t) <= OPS:
            if cur:
                cmds.append(cur)
            cur = []
        else:
            cur.append(t)
    if cur:
        cmds.append(cur)
    return cmds

# Short options of `git commit` that take an argument: the rest of the cluster,
# or else the next token, is that argument.
ARG_SHORT = set("mFCct")
OPT_ARG_SHORT = set("Su")                  # optional argument, attached only
ARG_LONG = {"--message", "--file", "--reuse-message", "--reedit-message",
            "--template", "--author", "--date", "--fixup", "--squash",
            "--cleanup", "--trailer", "--pathspec-from-file"}

def commit_bypass(after):
    i = 0
    while i < len(after):
        t = after[i]
        if t == "--":
            break
        if t.startswith("--"):
            name = t.split("=", 1)[0]
            if name.startswith("--no-veri"):
                return "--no-verify"
            if name in ARG_LONG and "=" not in t:
                i += 1
        elif t.startswith("-") and len(t) > 1:
            for k, ch in enumerate(t[1:]):
                if ch == "n":
                    return "-n"
                if ch in ARG_SHORT:
                    if k == len(t) - 2:
                        i += 1
                    break
                if ch in OPT_ARG_SHORT:
                    break
        i += 1
    return None

def scan(s):
    """Return (is_commit, bypass_reason)."""
    found, why = False, None
    for cmd in simple_commands(s):
        for gi, t in enumerate(cmd):
            if os.path.basename(t) != "git":
                continue
            for ci in range(gi + 1, len(cmd)):
                if cmd[ci] != "commit":
                    continue
                found = True
                before = cmd[:ci]
                if any("hookspath" in b.lower() for b in before):
                    why = why or "core.hooksPath"
                if any(b.startswith("GIT_CONFIG") for b in before):
                    why = why or "GIT_CONFIG_* in the environment"
                why = why or commit_bypass(cmd[ci + 1:])
    return found, why

texts = [strip_heredocs(c)]
if texts[0] != c:
    texts.append(c)                         # a misread heredoc must not hide a command
parsed, found, why = 0, False, None
for s in texts:
    try:
        f, w = scan(s)
    except ValueError:
        continue
    parsed += 1
    found = found or f
    why = why or w

if parsed == 0:
    looks = re.search(r"\bgit\b", c) and re.search(r"\bcommit\b", c)
    print("unparsed" if looks else "pass")
elif why:
    print("bypass " + why)
else:
    print("commit" if found else "pass")
')

case "$verdict" in
  pass) exit 0 ;;
  unparsed)
    cat >&2 <<'EOF'
BLOCKED — require_measurements: this command may run `git commit`, and it could not
be tokenised (an unbalanced quote outside any heredoc), so the hook cannot tell
whether it turns hooks off. Restate the command so it tokenises.
EOF
    exit 2 ;;
  bypass\ *)
    cat >&2 <<EOF
BLOCKED — require_measurements: this commit turns git's hooks off (${verdict#bypass }).
The measurements gate is git's commit-msg hook; a commit that skips it is ungated.
Commit without it. If a directive in the commit genuinely asserts nothing about any
artifact, say so in the message with:  no-measurements: <reason>
EOF
    exit 2 ;;
  commit) ;;
  *)
    echo "BLOCKED — require_measurements: unexpected verdict '$verdict' from the tokeniser." >&2
    exit 2 ;;
esac

proj=${CLAUDE_PROJECT_DIR:-/var/projects/toy_physics}
want="$proj/.claude/hooks/require_measurements_commit_msg.sh"
hooks=$(git -C "$proj" rev-parse --path-format=absolute --git-path hooks 2>/dev/null) || hooks=""
have=""
[ -n "$hooks" ] && have=$(readlink -f "$hooks/commit-msg" 2>/dev/null)
if [ -n "$have" ] && [ "$have" = "$(readlink -f "$want")" ] && [ -x "$want" ]; then
  exit 0
fi

cat >&2 <<EOF
BLOCKED — require_measurements: git's commit-msg hook is not the measurements gate.

The gate runs inside git, so a commit made without it is ungated. Expected
  ${hooks:-<git hooks directory>}/commit-msg  ->  $want
found
  ${have:-nothing}
Install it from the repository root with:
  ln -s ../../.claude/hooks/require_measurements_commit_msg.sh .git/hooks/commit-msg
(If core.hooksPath is set, git reads hooks from there instead; unset it.)
EOF
exit 2
