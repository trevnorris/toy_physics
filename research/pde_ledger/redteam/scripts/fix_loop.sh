#!/bin/bash
# fix_loop.sh — autonomous codex fix-pass over a list of directive_ready stages.
#
# Designed to run unattended (e.g., overnight). Halts cleanly on:
#   - codex non-zero exit
#   - any ## Blocked: block appearing in the directive (codex couldn't apply something)
#   - post-codex script re-exec failure (sanity check)
#   - stage in unexpected initial state
#
# A halt writes a marker file with the failing stage + reason; subsequent re-runs
# of this script will skip already-completed stages (codex_applied / verified) and
# resume from the failed one once the user fixes the underlying issue.
#
# Logs are tee'd to redteam/fix_batch_<BATCH>.log so the human can read what
# happened in the morning without needing to chase notifications.
#
# Usage:
#   bash redteam/scripts/fix_loop.sh I.1 002 003 004 005 006 007 008 009 010 011 012
#   bash redteam/scripts/fix_loop.sh I.1  # auto-discovers directive_ready stages in batch
#
# First arg is the batch ID (used for log file name). Remaining args are stage
# numbers in order. If only the batch ID is given, the script asks the manifest
# for all directive_ready stages in that batch.

set -uo pipefail

BATCH="${1:-}"
shift || true

if [[ -z "$BATCH" ]]; then
  echo "usage: fix_loop.sh <BATCH-ID> [stages...]" >&2
  exit 2
fi

# Resolve project root from this script's location (.../redteam/scripts/fix_loop.sh).
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$HERE/../.." && pwd)"
cd "$PROJECT_ROOT"

RT=/var/projects/toy_physics/.claude/skills/redteam-audit/lib/redteam.sh
MP=/var/projects/toy_physics/.claude/skills/redteam-audit/lib/manifest.py
LOG_BASENAME="redteam/fix_batch_${BATCH//\//_}.log"
MARKER="redteam/fix_batch_${BATCH//\//_}.marker"

mkdir -p redteam
# Preserve log across re-launches (helpful when resuming after a halt). The
# marker file is the canonical "result of last run" signal.
if [[ -f "$LOG_BASENAME" ]]; then
  printf '\n========== RESUME at %s ==========\n' "$(date -Iseconds)" >> "$LOG_BASENAME"
fi
rm -f "$MARKER"

log() {
  printf '[%s] %s\n' "$(date -Iseconds)" "$*" | tee -a "$LOG_BASENAME"
}

halt() {
  local stage="$1"; shift
  local reason="$*"
  log "HALT at stage $stage: $reason"
  printf 'HALTED\nstage: %s\nreason: %s\ntimestamp: %s\n' "$stage" "$reason" "$(date -Iseconds)" > "$MARKER"
  exit 1
}

success() {
  log "ALL_DONE: every stage in batch $BATCH reached codex_applied or beyond"
  printf 'ALL_DONE\nbatch: %s\ntimestamp: %s\n' "$BATCH" "$(date -Iseconds)" > "$MARKER"
  exit 0
}

# If no explicit stage list, ask the manifest for directive_ready stages in this batch.
if [[ "$#" -eq 0 ]]; then
  STAGES=$(bash "$RT" ls --state directive_ready 2>/dev/null | awk -v b="$BATCH" '$2==b{print $1}')
  if [[ -z "$STAGES" ]]; then
    log "No directive_ready stages found in batch $BATCH. Nothing to do."
    success
  fi
else
  STAGES="$*"
fi

log "Starting fix loop for batch $BATCH"
log "Stages to process: $STAGES"

# --- per-stage processing -------------------------------------------------

# Run sanity exec + verifier-prompt prep for a stage that's already in
# codex_applied (i.e., codex already ran but we didn't run the orchestrator's
# prep). Idempotent: safe to re-run; halts cleanly if scripts no longer pass.
ensure_prep_for_codex_applied() {
  local n="$1"
  local sympy_path math_path
  sympy_path=$(bash "$RT" paths "$n" | awk '/^sympy:/{print $2}')
  math_path=$(bash "$RT" paths "$n" | awk '/^mathematica:/{print $2}')

  if [[ -n "$sympy_path" && "$sympy_path" != "null" ]]; then
    bash "$RT" exec-sympy "$n" >> "$LOG_BASENAME" 2>&1
    local src=$?
    if [[ $src -ne 0 ]]; then
      halt "$n" "codex_applied stage's sympy script no longer exits 0 (exit $src)"
    fi
  fi
  if [[ -n "$math_path" && "$math_path" != "null" ]]; then
    bash "$RT" exec-mathematica "$n" >> "$LOG_BASENAME" 2>&1
    local mrc=$?
    if [[ $mrc -ne 0 ]]; then
      halt "$n" "codex_applied stage's mathematica script no longer exits 0 (exit $mrc)"
    fi
  fi

  mkdir -p redteam/tmp_prompts
  bash "$RT" render-verify-prompt "$n" > "redteam/tmp_prompts/verify_prompt_${n}.md"
  bash "$RT" capture-diff "$n" >> "$LOG_BASENAME"
  log "Stage $n: prep complete (sanity exec passed, verifier prompt rendered)"
  return 0
}

process_stage() {
  local n="$1"

  log "--- Stage $n: begin ---"

  local cur
  cur=$(bash "$RT" stage-info "$n" | awk '/^status:/{print $2}')

  case "$cur" in
    directive_ready)
      : # proceed
      ;;
    codex_applied)
      log "Stage $n already codex_applied — running sanity exec + verifier prep (idempotent)"
      # fall through to a streamlined codex_applied path below
      ensure_prep_for_codex_applied "$n"
      return $?
      ;;
    verified)
      log "Stage $n already verified — skipping entirely"
      return 0
      ;;
    needs_rework)
      # Resume-friendly: a previous needs_rework can be re-entered, but we
      # need the orchestrator (human or main session) to write a delta
      # directive first. Halt for review.
      halt "$n" "stage is in needs_rework — requires delta directive (manual)"
      ;;
    *)
      halt "$n" "unexpected pre-status: $cur"
      ;;
  esac

  # Confirm directive file exists
  local directive="redteam/directives/stage_${n}.md"
  if [[ ! -f "$directive" ]]; then
    halt "$n" "missing directive file: $directive"
  fi

  # set-status to fixing
  bash "$RT" set-status "$n" fixing 2>&1 | tee -a "$LOG_BASENAME"

  # Get current iteration count for the codex-invoke call (iter = current + 1)
  local iter
  iter=$(bash "$RT" stage-info "$n" | awk '/^iteration_count:/{print $2+1}')
  log "Stage $n: codex-invoke iter=$iter directive=$directive"

  # Run codex (per codex.md preamble, codex itself iterates against math -script
  # / python3 until exit 0 or hits the 5-iter cap with Blocked)
  bash "$RT" codex-invoke "$n" "$directive" "$iter" >> "$LOG_BASENAME" 2>&1
  local rc=$?

  if [[ $rc -ne 0 ]]; then
    halt "$n" "codex wrapper exited non-zero ($rc)"
  fi

  # Check for Blocked findings (codex couldn't make progress on at least one)
  if grep -q '^## Blocked:' "$directive" 2>/dev/null; then
    local blocked_count
    blocked_count=$(grep -c '^## Blocked:' "$directive")
    halt "$n" "codex marked $blocked_count finding(s) as Blocked — human review needed"
  fi

  # Re-scan after codex so the manifest knows about any file codex created
  # (e.g. a new .wl for a missing_mathematica finding). Without this, the
  # downstream paths lookup returns null and the sanity exec fails.
  bash "$RT" scan "$n" >> "$LOG_BASENAME" 2>&1

  # Sanity exec: independently confirm scripts run clean. (Codex claims it
  # iterated to exit 0, but trust-but-verify costs us seconds and catches the
  # case where codex's iteration logic short-circuited.)
  local sympy_path math_path
  sympy_path=$(bash "$RT" paths "$n" | awk '/^sympy:/{print $2}')
  math_path=$(bash "$RT" paths "$n" | awk '/^mathematica:/{print $2}')

  if [[ -n "$sympy_path" && "$sympy_path" != "null" ]]; then
    log "Stage $n: sanity exec sympy"
    bash "$RT" exec-sympy "$n" >> "$LOG_BASENAME" 2>&1
    local src=$?
    if [[ $src -ne 0 ]]; then
      halt "$n" "post-codex sympy exec failed (exit $src) — codex claimed success but script doesn't run"
    fi
  fi

  if [[ -n "$math_path" && "$math_path" != "null" ]]; then
    log "Stage $n: sanity exec mathematica"
    bash "$RT" exec-mathematica "$n" >> "$LOG_BASENAME" 2>&1
    local mrc=$?
    if [[ $mrc -ne 0 ]]; then
      halt "$n" "post-codex mathematica exec failed (exit $mrc) — codex claimed success but script doesn't run"
    fi
  fi

  bash "$RT" set-status "$n" codex_applied 2>&1 | tee -a "$LOG_BASENAME"
  log "Stage $n: codex_applied ✓"

  # Pre-render the verifier prompt so the orchestrator's main session can spawn
  # the verifier agent without an extra Bash round-trip. (Tucked under
  # tmp_prompts/ which is gitignored.)
  mkdir -p redteam/tmp_prompts
  bash "$RT" render-verify-prompt "$n" > "redteam/tmp_prompts/verify_prompt_${n}.md"
  bash "$RT" capture-diff "$n" >> "$LOG_BASENAME"
  log "Stage $n: verifier prompt rendered, diff captured"

  return 0
}

for n in $STAGES; do
  process_stage "$n"
done

success
