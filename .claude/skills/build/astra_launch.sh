#!/bin/bash
# Durable astra (codex) build launcher with a grok/claude KILL-WATCHDOG.
#
# WHY: a builder (codex/astra) given danger-full-access can spawn the `claude` (Claude Code)
# CLI and `grok` as its own "review" subprocesses — invalid (whatever writes does not review);
# it has happened twice (2026-09-05, 2026-09-07). ⚠ The COST concern is SPECIFICALLY `claude`:
# a codex-spawned Claude Code CLI is a NEW full-cost Claude session. `grok` CLI calls carry
# NO API cost here — grok is killed only because the self-review itself is invalid, not for cost.
# A prose "do not review" in the directive is only a REQUEST (CLAUDE.md S1); this enforces it
# structurally: a watchdog kills any claude-code / grok process in ASTRA'S OWN SESSION on sight.
# Session-scoped ⇒ it NEVER touches the orchestrator's own review legs or any other session on
# this shared box. astra (codex) and Mathematica (wolfram) are never killed.
#
# Usage:
#   setsid bash astra_launch.sh <DIRECTIVE_ABS> <SCRATCH_DIR_ABS> [<codex-model>] >/dev/null 2>&1 < /dev/null &
# Produces in SCRATCH_DIR: astra.log (+ EXIT=), astra.DONE, watchdog.log.
# The orchestrator MUST leak-gate the directive for process-words BEFORE calling this
# (see the /build skill) — the watchdog is the backstop, not the only control.
set -u
DIRECTIVE="$1"; SCR="$2"; MODEL="${3:-gpt-6-astra}"
REPO="/var/projects/toy_physics"
mkdir -p "$SCR"
WLOG="$SCR/watchdog.log"; : > "$WLOG"

# --- kill-watchdog: grok + claude-code in THIS session only ---
# FAIL-SAFE: only run if we are our own fresh setsid session leader (SID == our PID).
# If not (someone ran this without setsid, sharing the orchestrator's session), ABORT —
# never risk killing the orchestrator's own claude/grok.
watchdog() {
  local sid pid comm cl
  sid=$(ps -o sess= -p $$ | tr -d ' ')
  if [ "$sid" != "$$" ]; then
    echo "$(date +%T) WATCHDOG ABORTED: not a fresh setsid session (SID=$sid != PID=$$) — refusing to kill" >> "$WLOG"
    return
  fi
  while :; do
    for pid in $(pgrep -s "$sid" 2>/dev/null); do
      [ "$pid" = "$$" ] && continue
      comm=$(cat "/proc/$pid/comm" 2>/dev/null) || continue
      cl=$(tr '\0' ' ' < "/proc/$pid/cmdline" 2>/dev/null)
      case "$cl" in *codex*) continue ;; esac                 # astra itself — never kill
      case "$comm" in *[Ww]olfram*|MathKernel|wolframscript) continue ;; esac  # Mathematica — never kill
      if [ "$comm" = "grok" ] || { case "$cl" in *claude-code*|*anthropic-ai/claude*|*/claude\ *|claude\ *) true ;; *) false ;; esac; }; then
        kill -KILL "$pid" 2>/dev/null && echo "$(date +%T) WATCHDOG KILLED $pid ($comm): ${cl:0:100}" >> "$WLOG"
      fi
    done
    sleep 3
  done
}
watchdog & WD=$!

cd "$REPO"
codex exec -m "$MODEL" -c model_reasoning_effort=high --sandbox danger-full-access "$(<"$DIRECTIVE")" > "$SCR/astra.log" 2>&1 < /dev/null
echo "EXIT=$?" >> "$SCR/astra.log"
kill "$WD" 2>/dev/null
# final sweep: any grok/claude the watchdog's last cycle missed (same fresh-session fail-safe)
sid=$(ps -o sess= -p $$ | tr -d ' ')
[ "$sid" = "$$" ] || sid="__no_session__"   # abort sweep if not our own setsid session
for pid in $(pgrep -s "$sid" 2>/dev/null); do
  comm=$(cat "/proc/$pid/comm" 2>/dev/null); cl=$(tr '\0' ' ' < "/proc/$pid/cmdline" 2>/dev/null)
  case "$cl" in *codex*) continue ;; esac
  { [ "$comm" = "grok" ] || case "$cl" in *claude-code*|*anthropic-ai/claude*) true ;; *) false ;; esac; } && kill -KILL "$pid" 2>/dev/null
done
touch "$SCR/astra.DONE"
