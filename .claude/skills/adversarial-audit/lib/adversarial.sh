#!/usr/bin/env bash
#
# adversarial-audit — CLI entry point for the layer-2 fit-vs-derive audit skill.
#
# Subcommands:
#   init                         Create the adversarial report tree and manifest.
#   status                       Print manifest summary.
#   summary                      Alias for status.
#   candidate-info <ID>          Dump one candidate entry.
#   render-phase-a-prompts       Render blind Phase A modality prompts only.
#   phase-a-scan --stages ...    Run the blind Phase A modalities and union them.
#   phase-b-build <ID>           Build one parameter-value provenance slice.
#   phase-c-render <ID>          Render the adversarial prompt for one candidate.
#   set-status <ID> <STATUS>     Advance one candidate status under the manifest lock.
#   codex-defense <ID> <PARAM>   Invoke Codex defense, resuming per parameter.
#   graph <cmd> ...              Wrapper around graph/query_graph.py.
#   dry-run --stages 003 104 105 Run the non-binding A -> B -> C dry-run.
#   purge-dry-run <ID|all>       Remove dry-run entries/artifacts.
#   help                         Show this help.

set -euo pipefail

SKILL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LIB_DIR="$SKILL_DIR/lib"
CORE="$LIB_DIR/core.py"
PY_TIMEOUT=(timeout 600 python3)

find_config() {
  local dir="$PWD"
  while [[ "$dir" != "/" ]]; do
    if [[ -f "$dir/.redteam-config.yaml" ]]; then
      echo "$dir/.redteam-config.yaml"
      return 0
    fi
    if [[ -f "$dir/research/pde_ledger/.redteam-config.yaml" ]]; then
      echo "$dir/research/pde_ledger/.redteam-config.yaml"
      return 0
    fi
    dir="$(dirname "$dir")"
  done
  return 1
}

ensure_config_loaded() {
  if [[ -n "${CONFIG:-}" ]]; then return 0; fi
  CONFIG="$(find_config)" || {
    echo "error: no .redteam-config.yaml with adversarial section found from $PWD" >&2
    exit 1
  }
  PROJECT_ROOT="$(dirname "$CONFIG")"
  ARTIFACT_ROOT="$("${PY_TIMEOUT[@]}" "$CORE" "$CONFIG" artifact-root)"
  MANIFEST_LOCK="$ARTIFACT_ROOT/.manifest.lock"
}

_manifest_locked() {
  ensure_config_loaded
  mkdir -p "$ARTIFACT_ROOT"
  flock "$MANIFEST_LOCK" env ADVERSARIAL_MANIFEST_LOCKED=1 "$@"
}

run_core() {
  ensure_config_loaded
  "${PY_TIMEOUT[@]}" "$CORE" "$CONFIG" "$@"
}

cmd_init() {
  ensure_config_loaded
  _manifest_locked "${PY_TIMEOUT[@]}" "$CORE" "$CONFIG" init
}

cmd_status() {
  run_core status
}

cmd_candidate_info() {
  run_core candidate-info "$@"
}

cmd_render_phase_a_prompts() {
  run_core render-phase-a-prompts "$@"
}

cmd_phase_a_scan() {
  ensure_config_loaded
  _manifest_locked "${PY_TIMEOUT[@]}" "$CORE" "$CONFIG" phase-a-scan "$@"
}

cmd_phase_b_build() {
  ensure_config_loaded
  _manifest_locked "${PY_TIMEOUT[@]}" "$CORE" "$CONFIG" phase-b-build "$@"
}

cmd_phase_c_render() {
  ensure_config_loaded
  _manifest_locked "${PY_TIMEOUT[@]}" "$CORE" "$CONFIG" phase-c-render "$@"
}

cmd_set_status() {
  ensure_config_loaded
  _manifest_locked "${PY_TIMEOUT[@]}" "$CORE" "$CONFIG" set-status "$@"
}

cmd_dry_run() {
  ensure_config_loaded
  _manifest_locked "${PY_TIMEOUT[@]}" "$CORE" "$CONFIG" dry-run "$@"
}

cmd_purge_dry_run() {
  ensure_config_loaded
  _manifest_locked "${PY_TIMEOUT[@]}" "$CORE" "$CONFIG" purge-dry-run "$@"
}

cmd_graph() {
  ensure_config_loaded
  local query_graph graph_path
  query_graph="$(run_core config-path query_graph)"
  graph_path="$(run_core config-path atlas_graph)"
  "${PY_TIMEOUT[@]}" "$query_graph" --graph "$graph_path" "$@"
}

cmd_codex_defense() {
  ensure_config_loaded
  if [[ $# -lt 2 ]]; then
    echo "usage: codex-defense <CANDIDATE-ID> <PARAMETER> [ITER]" >&2
    exit 2
  fi
  local candidate_id="$1"
  local parameter="$2"
  local iter="${3:-1}"
  local prompt_path log_path wrapper sandbox session rc new_session

  prompt_path="$(run_core render-defense-prompt "$candidate_id" "$parameter")"
  log_path="$(run_core codex-log-path "$candidate_id" "$parameter" "$iter")"
  wrapper="$(run_core config-path codex.chat_wrapper)"
  sandbox="$(run_core config-value codex.sandbox)"
  session="$(run_core parameter-session "$candidate_id" "$parameter")"

  if [[ ! -x "$wrapper" ]]; then
    echo "error: codex wrapper not executable: $wrapper" >&2
    exit 1
  fi
  mkdir -p "$(dirname "$log_path")"
  echo "codex-defense: candidate=$candidate_id parameter=$parameter iter=$iter session=${session:-NEW} log=$log_path"

  set +e
  if [[ -n "$session" && "$session" != "null" ]]; then
    timeout 600 "$wrapper" --resume "$session" -C "$PROJECT_ROOT" < "$prompt_path" > "$log_path" 2>&1
  else
    timeout 600 "$wrapper" -s "$sandbox" -C "$PROJECT_ROOT" < "$prompt_path" > "$log_path" 2>&1
  fi
  rc=$?
  set -e

  new_session="$(grep -oP 'codex_session_id: \K[0-9a-f-]+' "$log_path" | tail -1 || true)"
  _manifest_locked "${PY_TIMEOUT[@]}" "$CORE" "$CONFIG" record-codex-defense \
    "$candidate_id" "$parameter" "$iter" "$rc" "$log_path" "$new_session"
  echo "codex-defense: exit=$rc"
  return "$rc"
}

cmd_help() {
  sed -n '3,28p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
}

main() {
  local sub="${1:-help}"; shift || true
  case "$sub" in
    init)              cmd_init "$@" ;;
    status|summary)    cmd_status "$@" ;;
    candidate-info)    cmd_candidate_info "$@" ;;
    render-phase-a-prompts) cmd_render_phase_a_prompts "$@" ;;
    phase-a-scan)      cmd_phase_a_scan "$@" ;;
    phase-b-build)     cmd_phase_b_build "$@" ;;
    phase-c-render)    cmd_phase_c_render "$@" ;;
    set-status)        cmd_set_status "$@" ;;
    codex-defense)     cmd_codex_defense "$@" ;;
    graph)             cmd_graph "$@" ;;
    dry-run)           cmd_dry_run "$@" ;;
    purge-dry-run)     cmd_purge_dry_run "$@" ;;
    help|--help|-h)    cmd_help ;;
    *) echo "unknown subcommand: $sub" >&2; cmd_help; exit 2 ;;
  esac
}

main "$@"
