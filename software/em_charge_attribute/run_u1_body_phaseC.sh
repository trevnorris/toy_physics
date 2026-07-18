#!/usr/bin/env bash

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

if [[ "${1:-}" == "--shell-probe" ]]; then
  exec /usr/bin/python3 -I "$script_dir/u1_phaseC_stage0_runner.py" --self-probe
fi

if [[ "${1:-}" == "--production-shell-probe" ]]; then
  exec /usr/bin/python3 -I "$script_dir/u1_phaseC_stage0_runner.py" --production-self-probe
fi

exec /usr/bin/python3 -I "$script_dir/u1_phaseC_stage0_runner.py" "$@"
