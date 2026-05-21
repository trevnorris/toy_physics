#!/usr/bin/env bash
#
# redteam — CLI entry point for the redteam-audit skill.
#
# Reads .redteam-config.yaml from the current working directory. All paths
# are resolved relative to that file's location.
#
# Subcommands:
#   bootstrap [PROJECT-DIR]    Generate a starter .redteam-config.yaml by scanning the
#                              project tree. Default project-dir is current working dir.
#                              Pass --force to overwrite an existing config, --batch-size N
#                              to control auto-batching, --out PATH to write elsewhere.
#   detect [PROJECT-DIR]       Print the raw detection report (YAML) without writing a config.
#                              Useful for debugging or feeding into a wrapper tool.
#   init                       Seed redteam/MANIFEST.yaml and BATCHES.md from config.
#   status [--batch ID|--state S]   Print manifest summary.
#   scan <stage|batch:ID|all>  Re-scan filesystem for script presence.
#   ls --state S               List stages in state S, one per line (stage<TAB>batch<TAB>status).
#   batch-info <ID>            Print batch metadata + stage list.
#   stage-info <NNN>           Dump a single stage's manifest entry as YAML.
#   set-status <NNN> <STATE>   Manually set a stage's status (use sparingly).
#   set-iter <NNN> <N>         Set iteration count for a stage.
#   mark-stale-downstream <NNN>  Mark every stage > NNN as upstream_stale=true.
#   exec-sympy <NNN>           Run a stage's SymPy script via configured runner.
#   exec-mathematica <NNN>     Run a stage's Mathematica script via configured runner.
#   codex-invoke <NNN> <DIRECTIVE-PATH> [ITER]   Invoke codex against the directive.
#                              Resumes the unit's session if one exists; otherwise starts
#                              a new session. Captures log to redteam/codex_logs/<NNN>_iter<N>.txt
#                              and stores the session id back in the manifest.
#   codex-reset <NNN>          Clear the unit's codex_session field (forces next codex-invoke
#                              to start fresh). Use after upstream changes invalidate context.
#   render-audit-prompt <NNN>  Output prompts/auditor.md with all {PLACEHOLDER} substitutions
#                              filled in for the given unit. The orchestrator pipes this to
#                              the Agent tool as the audit prompt.
#   render-verify-prompt <NNN> Same for prompts/verifier.md (report/directive/log/diff paths).
#   capture-diff <NNN>         git diff HEAD of the unit's script files; writes to
#                              redteam/exec_logs/stage_<NNN>_diff.patch; prints the path.
#   next-batch                 Print the earliest batch with any audit-eligible unit
#                              (status in {pending, upstream_stale}).
#   blocked                    List all stages in any blocked_* state.
#   paths <NNN>                Print the resolved file paths for a stage (YAML).
#   render-batches             Regenerate BATCHES.md from the manifest.
#   help                       Show this help.
#
# Valid stage statuses (state machine):
#   pending, auditing, audited, directive_ready, fixing, codex_applied,
#   verifying, verified, needs_rework, blocked_unfixable,
#   blocked_critical_downstream, upstream_stale

set -euo pipefail

SKILL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LIB_DIR="$SKILL_DIR/lib"
YAML="python3 $LIB_DIR/_yaml.py"
SCAN="python3 $LIB_DIR/scan.py"
MANIFEST_PY="python3 $LIB_DIR/manifest.py"
DETECT_PY="python3 $LIB_DIR/detect.py"
BOOTSTRAP_PY="python3 $LIB_DIR/bootstrap.py"
RENDER_PY="python3 $LIB_DIR/render_prompt.py"

# --- locate config -----------------------------------------------------------
# bootstrap / detect / help work without a config; everything else requires one.

find_config() {
  local dir="$PWD"
  while [[ "$dir" != "/" ]]; do
    if [[ -f "$dir/.redteam-config.yaml" ]]; then
      echo "$dir/.redteam-config.yaml"
      return 0
    fi
    dir="$(dirname "$dir")"
  done
  return 1
}

ensure_config_loaded() {
  if [[ -n "${CONFIG:-}" ]]; then return 0; fi
  CONFIG="$(find_config)" || {
    echo "error: no .redteam-config.yaml found in $PWD or any parent" >&2
    echo "       run '$(basename "${BASH_SOURCE[0]}") bootstrap' first to generate one" >&2
    exit 1
  }
  PROJECT_ROOT="$(dirname "$CONFIG")"
  REDTEAM_DIR_REL="$($YAML get "$CONFIG" project.redteam_dir)"
  REDTEAM_DIR="$PROJECT_ROOT/$REDTEAM_DIR_REL"
  MANIFEST="$REDTEAM_DIR/MANIFEST.yaml"
  BATCHES_MD="$REDTEAM_DIR/BATCHES.md"
  STAGE_PAD="$($YAML get "$CONFIG" stage_pad)"
}

pad_stage() {
  printf "%0${STAGE_PAD}d" "$1"
}

# --- batch / stage helpers ---------------------------------------------------

batch_for_stage() {
  local n=$1
  $YAML batches "$CONFIG" | while IFS=$'\t' read -r id start end label; do
    if (( n >= start && n <= end )); then
      echo "$id"
      return 0
    fi
  done
}

all_stages() {
  local start end excl
  start="$($YAML get "$CONFIG" stages.start)"
  end="$($YAML get "$CONFIG" stages.end)"
  local exclude
  exclude="$($YAML get-list "$CONFIG" stages.exclude 2>/dev/null || true)"
  for ((n=start; n<=end; n++)); do
    if echo "$exclude" | grep -qx "$n"; then continue; fi
    echo "$n"
  done
}

is_checkpoint() {
  local n=$1
  local cps
  cps="$($YAML get-list "$CONFIG" checkpoints 2>/dev/null || true)"
  echo "$cps" | grep -qx "$n"
}

is_status_only_candidate() {
  local n=$1
  local cands
  cands="$($YAML get-list "$CONFIG" status_only_candidates 2>/dev/null || true)"
  echo "$cands" | grep -qx "$n"
}

# --- manifest mutation -------------------------------------------------------

ensure_manifest() {
  mkdir -p "$REDTEAM_DIR"/{reports,directives,verifications,batches,exec_logs,codex_logs}
  $YAML init "$MANIFEST"
}

manifest_set_stage_field() {
  # set-status, set-iter, etc. all funnel through here.
  local stage_str=$1
  local field=$2
  local value=$3
  $YAML set "$MANIFEST" "stages.$stage_str.$field" "$value"
}

manifest_get_stage_field() {
  local stage_str=$1
  local field=$2
  $YAML get "$MANIFEST" "stages.$stage_str.$field" 2>/dev/null || echo ""
}

# --- subcommands -------------------------------------------------------------

cmd_bootstrap() {
  # First positional arg may be a project dir; rest are passed through to bootstrap.py.
  local target="."
  local args=()
  if [[ $# -gt 0 && "$1" != --* ]]; then
    target="$1"
    shift
  fi
  args+=("$(realpath "$target")")
  args+=("$@")
  $BOOTSTRAP_PY "${args[@]}"
}

cmd_detect() {
  local target="${1:-.}"
  $DETECT_PY "$(realpath "$target")"
}

cmd_init() {
  ensure_config_loaded
  ensure_manifest

  # Seed manifest header if absent.
  if ! $YAML has "$MANIFEST" project_name; then
    local proj
    proj="$($YAML get "$CONFIG" project.name)"
    $YAML set "$MANIFEST" project_name "$proj" --raw
    $YAML set "$MANIFEST" schema_version 2
    $YAML set "$MANIFEST" created_at "$(date -Iseconds)" --raw
  fi

  local n stage_str batch_id scan_out
  local created=0 updated=0
  while read -r n; do
    stage_str="$(pad_stage "$n")"
    batch_id="$(batch_for_stage "$n")"
    # Only seed if not already present (idempotent).
    if $YAML has "$MANIFEST" "stages.$stage_str"; then
      updated=$((updated+1))
      continue
    fi
    scan_out="$($SCAN "$CONFIG" "$n")"
    # Build the stage entry as YAML and pipe to set-yaml.
    {
      echo "stage: $n"
      echo "stage_str: \"$stage_str\""
      echo "batch_id: \"$batch_id\""
      echo "status: pending"
      echo "iteration_count: 0"
      echo "upstream_stale: false"
      echo "is_checkpoint: $(is_checkpoint "$n" && echo true || echo false)"
      echo "is_status_only_candidate: $(is_status_only_candidate "$n" && echo true || echo false)"
      echo "last_audit_date: null"
      echo "last_audit_findings: 0"
      echo "last_verify_date: null"
      echo "files:"
      echo "$scan_out" | sed 's/^/  /' | sed '/^  stage:/d;/^  stage_str:/d'
    } | $YAML set-yaml "$MANIFEST" "stages.$stage_str"
    created=$((created+1))
  done < <(all_stages)

  echo "init: created=$created, already-present=$updated"
  cmd_render_batches
}

cmd_scan() {
  ensure_config_loaded
  local target=${1:-all}
  local stages=()
  if [[ "$target" == "all" ]]; then
    while read -r n; do stages+=("$n"); done < <(all_stages)
  elif [[ "$target" =~ ^batch: ]]; then
    local bid="${target#batch:}"
    while read -r s; do stages+=("$((10#$s))"); done < <($YAML stages-in-batch "$CONFIG" "$bid")
  else
    stages+=("$((10#$target))")
  fi
  local n stage_str scan_out
  for n in "${stages[@]}"; do
    stage_str="$(pad_stage "$n")"
    scan_out="$($SCAN "$CONFIG" "$n")"
    # Strip the stage/stage_str header since the manifest entry already has them.
    local files_yaml
    files_yaml=$(echo "$scan_out" | sed '/^stage:/d;/^stage_str:/d')
    echo "$files_yaml" | $YAML set-yaml "$MANIFEST" "stages.$stage_str.files"
  done
  echo "scan: ${#stages[@]} stage(s) updated"
}

cmd_status() {
  ensure_config_loaded
  local mode="summary"
  local arg=""
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --batch) mode="batch"; arg="$2"; shift 2 ;;
      --state) mode="state"; arg="$2"; shift 2 ;;
      *) echo "unknown arg: $1" >&2; exit 2 ;;
    esac
  done
  case "$mode" in
    summary) $MANIFEST_PY "$CONFIG" summary ;;
    batch)   $MANIFEST_PY "$CONFIG" batch-status "$arg" ;;
    state)   $MANIFEST_PY "$CONFIG" state-list "$arg" ;;
  esac
}

cmd_ls() {
  ensure_config_loaded
  local target=""
  if [[ "${1:-}" == "--state" ]]; then target="$2"; fi
  $MANIFEST_PY "$CONFIG" state-list "$target"
}

cmd_batch_info() {
  ensure_config_loaded
  local bid=$1
  $YAML batches "$CONFIG" | while IFS=$'\t' read -r id start end label; do
    if [[ "$id" == "$bid" ]]; then
      echo "id: $id"
      echo "range: [$start, $end]"
      echo "label: $label"
      echo "stages:"
      $YAML stages-in-batch "$CONFIG" "$bid" | sed 's/^/  - /'
    fi
  done
}

cmd_stage_info() {
  ensure_config_loaded
  local stage_str
  stage_str="$(pad_stage "$((10#$1))")"
  $YAML get "$MANIFEST" "stages.$stage_str"
}

cmd_set_status() {
  ensure_config_loaded
  local stage_str=$(pad_stage "$((10#$1))")
  local new=$2
  manifest_set_stage_field "$stage_str" status "$new"
  manifest_set_stage_field "$stage_str" last_status_change "$(date -Iseconds)"
  # Keep upstream_stale flag synced with status. The flag is a denormalization
  # of "status == upstream_stale" — without this sync, re-auditing a stale
  # unit and moving it back to verified would leave the flag stuck at true.
  if [[ "$new" == "upstream_stale" ]]; then
    manifest_set_stage_field "$stage_str" upstream_stale true
  else
    manifest_set_stage_field "$stage_str" upstream_stale false
  fi
  echo "set-status: $stage_str → $new"
}

cmd_set_iter() {
  ensure_config_loaded
  local stage_str=$(pad_stage "$((10#$1))")
  manifest_set_stage_field "$stage_str" iteration_count "$2"
}

cmd_mark_stale_downstream() {
  ensure_config_loaded
  $MANIFEST_PY "$CONFIG" mark-stale-downstream "$1"
}

cmd_paths() {
  ensure_config_loaded
  local stage_str
  stage_str="$(pad_stage "$((10#$1))")"
  $YAML get "$MANIFEST" "stages.$stage_str.files"
}

cmd_exec_script() {
  ensure_config_loaded
  local engine=$1   # sympy | mathematica
  local stage_str
  stage_str="$(pad_stage "$((10#$2))")"
  local script_path
  script_path="$(_path_field "$stage_str" "files.$engine.path" "")"
  if [[ -z "$script_path" ]]; then
    echo "error: no $engine script for stage $stage_str" >&2
    return 1
  fi
  local runner
  runner="$($YAML get "$CONFIG" "runners.$engine")"
  if [[ -z "$runner" || "$runner" == "null" ]]; then
    echo "error: no runner configured for engine '$engine'" >&2
    return 1
  fi
  # Build argv by splitting the runner on whitespace and substituting {script}.
  # Avoids eval / shell-quoting hazards from paths or runner strings.
  local -a runner_tokens=()
  read -ra runner_tokens <<< "$runner"
  local -a argv=()
  local tok
  for tok in "${runner_tokens[@]}"; do
    argv+=("${tok//\{script\}/$script_path}")
  done

  local log="$REDTEAM_DIR/exec_logs/stage_${stage_str}_${engine}.log"
  mkdir -p "$(dirname "$log")"
  echo "exec: ${argv[*]}"
  {
    echo "# stage: $stage_str"
    echo "# engine: $engine"
    echo "# argv: ${argv[*]@Q}"
    echo "# date: $(date -Iseconds)"
    echo "---"
  } > "$log"
  set +e
  ( cd "$PROJECT_ROOT" && "${argv[@]}" ) >> "$log" 2>&1
  local rc=$?
  set -e
  echo "# exit_code: $rc" >> "$log"
  manifest_set_stage_field "$stage_str" "last_${engine}_exit" "$rc"
  manifest_set_stage_field "$stage_str" "last_${engine}_run" "$(date -Iseconds)"
  echo "exit: $rc  log: $log"
  return $rc
}

cmd_codex_invoke() {
  ensure_config_loaded
  if [[ $# -lt 2 ]]; then
    echo "usage: codex-invoke <NNN> <DIRECTIVE-PATH> [ITER]" >&2
    exit 2
  fi
  local unit
  unit="$(pad_stage "$((10#$1))")"
  local directive_path=$2
  local iter=${3:-1}

  if [[ ! -f "$directive_path" ]]; then
    echo "error: directive not found: $directive_path" >&2
    exit 1
  fi

  mkdir -p "$REDTEAM_DIR/codex_logs"
  local log="$REDTEAM_DIR/codex_logs/${unit}_iter${iter}.txt"

  local wrapper
  wrapper="$($YAML get "$CONFIG" codex.chat_wrapper)"
  wrapper="${wrapper/#\~/$HOME}"
  if [[ ! -x "$wrapper" ]]; then
    echo "error: codex wrapper not executable: $wrapper" >&2
    exit 1
  fi
  local sandbox
  sandbox="$($YAML get "$CONFIG" codex.sandbox)"

  local session
  session="$(manifest_get_stage_field "$unit" codex_session)"

  echo "codex: unit=$unit iter=$iter session=${session:-NEW} log=$log"

  set +e
  if [[ -n "$session" && "$session" != "null" ]]; then
    # Resume — codex already has its system prompt cached in the session.
    "$wrapper" --resume "$session" -C "$PROJECT_ROOT" < "$directive_path" > "$log" 2>&1
  else
    # New session — prepend codex.md preamble (system prompt) to the directive.
    local combined
    combined="$(mktemp)"
    sed -e "s|{UNIT_ID}|$unit|g" \
        -e "s|{DIRECTIVE_PATH}|$directive_path|g" \
        "$SKILL_DIR/prompts/codex.md" > "$combined"
    {
      echo ""
      echo "---"
      echo ""
      echo "# Active directive"
      echo ""
      cat "$directive_path"
    } >> "$combined"
    "$wrapper" -s "$sandbox" -C "$PROJECT_ROOT" < "$combined" > "$log" 2>&1
    rm -f "$combined"
  fi
  local rc=$?
  set -e

  # Parse session id from wrapper output (appended as `codex_session_id: <uuid>`).
  local new_session
  new_session=$(grep -oP 'codex_session_id: \K[0-9a-f-]+' "$log" | tail -1 || true)

  if [[ -n "$new_session" ]]; then
    manifest_set_stage_field "$unit" codex_session "$new_session"
  elif [[ -n "$session" && "$session" != "null" && $rc -ne 0 ]]; then
    # Resume failed; session likely expired. Clear it so the next call starts fresh.
    echo "warning: codex --resume failed (session may have expired); clearing session id" >&2
    manifest_set_stage_field "$unit" codex_session ""
  fi

  # Append log path to codex_log_paths list (atomic, via manifest.py).
  local rel_log="${log#$PROJECT_ROOT/}"
  echo "$rel_log" | $MANIFEST_PY "$CONFIG" append-stage-list "$unit" codex_log_paths

  manifest_set_stage_field "$unit" iteration_count "$iter"
  manifest_set_stage_field "$unit" last_codex_run "$(date -Iseconds)"
  manifest_set_stage_field "$unit" last_codex_exit "$rc"

  echo "codex: exit=$rc"
  return $rc
}

cmd_codex_reset() {
  ensure_config_loaded
  local unit
  unit="$(pad_stage "$((10#$1))")"
  manifest_set_stage_field "$unit" codex_session ""
  echo "codex-reset: session cleared for unit $unit"
}

# Internal: print a manifest field value with a sane fallback string if missing.
_field_or() {
  local unit=$1 field=$2 fallback=$3
  local v
  v="$(manifest_get_stage_field "$unit" "$field")"
  if [[ -z "$v" || "$v" == "null" ]]; then
    echo "$fallback"
  else
    echo "$v"
  fi
}

# Internal: resolve a manifest path (relative since schema v2) to absolute.
# Pre-v2 manifests may have absolute paths already; if so, return as-is.
# Fallback strings like "(missing)" pass through untouched.
_resolve_path() {
  local p=$1
  if [[ -z "$p" || "$p" == "null" || "$p" == "(missing)" ]]; then
    echo "$p"
    return
  fi
  if [[ "$p" == /* ]]; then
    echo "$p"
  else
    echo "$PROJECT_ROOT/$p"
  fi
}

# Internal: fetch a path field and resolve to absolute in one step.
_path_field() {
  local unit=$1 field=$2 fallback=$3
  _resolve_path "$(_field_or "$unit" "$field" "$fallback")"
}

cmd_render_audit_prompt() {
  ensure_config_loaded
  if [[ $# -lt 1 ]]; then
    echo "usage: render-audit-prompt <NNN>" >&2
    exit 2
  fi
  local unit
  unit="$(pad_stage "$((10#$1))")"

  local batch_id is_checkpoint is_status_only
  batch_id="$(_field_or "$unit" batch_id "")"
  is_checkpoint="$(_field_or "$unit" is_checkpoint false)"
  is_status_only="$(_field_or "$unit" is_status_only_candidate false)"

  local sympy_path math_path sympy_out math_out
  sympy_path="$(_path_field "$unit" files.sympy.path "(missing)")"
  math_path="$(_path_field "$unit" files.mathematica.path "(missing)")"
  sympy_out="$(_path_field "$unit" files.sympy_output.path "(missing)")"
  math_out="$(_path_field "$unit" files.mathematica_output.path "(missing)")"

  local report_path="$REDTEAM_DIR/reports/stage_${unit}.md"
  local directive_path="$REDTEAM_DIR/directives/stage_${unit}.md"

  $RENDER_PY "$SKILL_DIR/prompts/auditor.md" \
    "UNIT_ID=$unit" \
    "STAGE_STR=$unit" \
    "BATCH_ID=$batch_id" \
    "IS_CHECKPOINT=$is_checkpoint" \
    "IS_STATUS_ONLY_CANDIDATE=$is_status_only" \
    "SYMPY_PATH=$sympy_path" \
    "MATH_PATH=$math_path" \
    "SYMPY_OUT_PATH=$sympy_out" \
    "MATH_OUT_PATH=$math_out" \
    "REPORT_PATH=$report_path" \
    "DIRECTIVE_PATH=$directive_path"
}

cmd_render_verify_prompt() {
  ensure_config_loaded
  if [[ $# -lt 1 ]]; then
    echo "usage: render-verify-prompt <NNN>" >&2
    exit 2
  fi
  local unit
  unit="$(pad_stage "$((10#$1))")"

  local batch_id
  batch_id="$(_field_or "$unit" batch_id "")"

  local report_path="$REDTEAM_DIR/reports/stage_${unit}.md"
  local directive_path="$REDTEAM_DIR/directives/stage_${unit}.md"
  local sympy_log="$REDTEAM_DIR/exec_logs/stage_${unit}_sympy.log"
  local math_log="$REDTEAM_DIR/exec_logs/stage_${unit}_mathematica.log"
  local diff_path="$REDTEAM_DIR/exec_logs/stage_${unit}_diff.patch"
  local verification_path="$REDTEAM_DIR/verifications/stage_${unit}.md"

  $RENDER_PY "$SKILL_DIR/prompts/verifier.md" \
    "UNIT_ID=$unit" \
    "STAGE_STR=$unit" \
    "BATCH_ID=$batch_id" \
    "REPORT_PATH=$report_path" \
    "DIRECTIVE_PATH=$directive_path" \
    "SYMPY_LOG_PATH=$sympy_log" \
    "MATH_LOG_PATH=$math_log" \
    "DIFF_PATH=$diff_path" \
    "VERIFICATION_PATH=$verification_path"
}

cmd_capture_diff() {
  ensure_config_loaded
  if [[ $# -lt 1 ]]; then
    echo "usage: capture-diff <NNN>" >&2
    exit 2
  fi
  local unit
  unit="$(pad_stage "$((10#$1))")"

  local diff_path="$REDTEAM_DIR/exec_logs/stage_${unit}_diff.patch"
  mkdir -p "$(dirname "$diff_path")"

  local sympy_path math_path
  sympy_path="$(_path_field "$unit" files.sympy.path "")"
  math_path="$(_path_field "$unit" files.mathematica.path "")"

  local paths=()
  [[ -n "$sympy_path" ]] && paths+=("$sympy_path")
  [[ -n "$math_path"  ]] && paths+=("$math_path")

  if [[ ${#paths[@]} -eq 0 ]]; then
    : > "$diff_path"
    echo "$diff_path"
    return 0
  fi

  # Diff vs HEAD captures whatever has changed since the last commit.
  # Assumes the working tree was clean before codex ran; if not, the verifier
  # sees pre-existing changes mixed in and the user is responsible for that.
  ( cd "$PROJECT_ROOT" && git diff HEAD -- "${paths[@]}" ) > "$diff_path" 2>/dev/null || :
  echo "$diff_path"
}

cmd_next_batch() {
  ensure_config_loaded
  $MANIFEST_PY "$CONFIG" next-batch
}

cmd_blocked() {
  ensure_config_loaded
  $MANIFEST_PY "$CONFIG" blocked
}

cmd_render_batches() {
  ensure_config_loaded
  $MANIFEST_PY "$CONFIG" render-batches
}

cmd_help() {
  sed -n '3,48p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
}

# --- dispatch ----------------------------------------------------------------

main() {
  local sub="${1:-help}"; shift || true
  case "$sub" in
    bootstrap)              cmd_bootstrap "$@" ;;
    detect)                 cmd_detect "$@" ;;
    init)                   cmd_init "$@" ;;
    scan)                   cmd_scan "${1:-all}" ;;
    status)                 cmd_status "$@" ;;
    ls)                     cmd_ls "$@" ;;
    batch-info)             cmd_batch_info "$@" ;;
    stage-info)             cmd_stage_info "$@" ;;
    set-status)             cmd_set_status "$@" ;;
    set-iter)               cmd_set_iter "$@" ;;
    mark-stale-downstream)  cmd_mark_stale_downstream "$@" ;;
    exec-sympy)             cmd_exec_script sympy "$@" ;;
    exec-mathematica)       cmd_exec_script mathematica "$@" ;;
    codex-invoke)           cmd_codex_invoke "$@" ;;
    codex-reset)            cmd_codex_reset "$@" ;;
    render-audit-prompt)    cmd_render_audit_prompt "$@" ;;
    render-verify-prompt)   cmd_render_verify_prompt "$@" ;;
    capture-diff)           cmd_capture_diff "$@" ;;
    next-batch)             cmd_next_batch ;;
    blocked)                cmd_blocked ;;
    paths)                  cmd_paths "$@" ;;
    render-batches)         cmd_render_batches ;;
    help|--help|-h)         cmd_help ;;
    *) echo "unknown subcommand: $sub" >&2; cmd_help; exit 2 ;;
  esac
}

main "$@"
