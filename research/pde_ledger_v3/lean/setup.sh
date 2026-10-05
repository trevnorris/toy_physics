#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"

if [[ ${1:-} == --help ]]; then
  cat <<'EOF'
Usage: bash setup.sh [--check]
Install the pinned Lean/Python dependencies; see INSTALL.md.
--check checks prerequisites only, without downloads or installation.
Setup does not run the proof suites. Afterwards use:
  .venv/bin/python verify.py all
PYTHON selects Python 3.10+ (default: python3).
EOF
  exit 0
fi
if [[ $# -gt 1 || ( $# -eq 1 && $1 != --check ) ]]; then
  echo 'Unknown option; use --help.' >&2
  exit 2
fi
for cmd in git curl zstd tar cc; do
  command -v "$cmd" >/dev/null || { echo "Missing $cmd; see INSTALL.md." >&2; exit 1; }
done
pde_python=${PYTHON:-python3}
"$pde_python" -c 'import sys, venv; sys.exit(0 if sys.version_info >= (3,10) else "Python 3.10+ required")'
for pin in lean-toolchain lake-manifest.json lakefile.toml requirements.txt; do
  [[ -s $pin ]] || { echo "Missing $pin; restore it from this checkout (do not run lake update)." >&2; exit 1; }
done
if [[ ${1:-} == --check ]]; then
  echo 'Prerequisites available. No downloads or installation performed.'
  exit 0
fi

# Keep command output in a durable log; surface failure and its location.
exec 3>&1
exec >>setup.log 2>&1
trap 'rc=$?; echo "Setup failed (exit $rc). See $(pwd)/setup.log" >&3; exit "$rc"' ERR
echo "Setup started: $(date -u +%FT%TZ)"
export PATH="${ELAN_HOME:-$HOME/.elan}/bin:$PATH"
if ! command -v elan >/dev/null; then
  pde_installer=$(mktemp)
  trap 'rm -f -- "${pde_installer:-}"' EXIT
  curl --proto '=https' --tlsv1.2 -fsSL https://elan.lean-lang.org/elan-init.sh -o "$pde_installer"
  sh "$pde_installer" -y --default-toolchain none --no-modify-path
fi
# Never change the global default or resolve newer dependency versions.
elan toolchain install "$(cat lean-toolchain)"
"$pde_python" -m venv .venv
.venv/bin/python -m pip install --disable-pip-version-check -r requirements.txt
export LAKE_CACHE_DIR="$PWD/.lake/cache"
lake exe cache get
# Cache retrieval supplies Mathlib. Build only the three Physlib imports used
# by this ledger, not Physlib's entire library or any local proof target.
lake build +Physlib.ClassicalMechanics.WaveEquation.Basic:olean \
  +Physlib.Mathematics.LeviCivita.Basic:olean +Physlib.Units.Dimension:olean
.venv/bin/python verify.py --doctor
echo 'Setup complete. Run .venv/bin/python verify.py all (see INSTALL.md).' >&3
