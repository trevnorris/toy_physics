#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
runner="$repo_root/software/em_charge_attribute/u2_stage0_runner.py"

if [[ "${1:-}" == "--shell-probe" ]]; then
  exec /usr/bin/python3 -I "$runner" --repo "$repo_root" --self-probe
fi

if [[ "${1:-}" == "--production-shell-probe" ]]; then
  exec /usr/bin/python3 -I "$runner" --repo "$repo_root" --production-self-probe
fi

exec /usr/bin/python3 -I "$runner" --repo "$repo_root" "$@"
