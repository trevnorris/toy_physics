#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"

# Install the project-selected version without changing elan's global default.
elan toolchain install "$(cat lean-toolchain)"
if [[ ! -f lake-manifest.json ]]; then
  lake update
fi
# Fetch precompiled Mathlib dependencies; compile only the Physlib imports used here.
lake exe cache get
lake build
