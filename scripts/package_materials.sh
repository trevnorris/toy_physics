#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"
research_root="$repo_root/research"

if [[ ! -d "$research_root" ]]; then
  echo "ERROR: research directory not found at $research_root" >&2
  exit 1
fi

if ! command -v zip >/dev/null 2>&1; then
  echo "ERROR: zip is not installed." >&2
  exit 1
fi

mapfile -d '' paper_dirs < <(find "$research_root" -mindepth 2 -maxdepth 2 -type d -name paper -print0 | sort -z)

if [[ "${#paper_dirs[@]}" -eq 0 ]]; then
  echo "No research/*/paper directories found under $research_root" >&2
  exit 1
fi

for paper_dir in "${paper_dirs[@]}"; do
  bundle_dir="$(dirname "$paper_dir")"
  bundle_name="$(basename "$bundle_dir")"
  zip_path="$bundle_dir/materials.zip"
  tmp_dir="$(mktemp -d "${TMPDIR:-/tmp}/${bundle_name}.materials.XXXXXX")"
  tmp_zip="$tmp_dir/materials.zip"

  echo "==> $bundle_name"
  echo "    cleaning paper build artifacts"
  find "$paper_dir" -maxdepth 1 -type f \( -name '*.aux' -o -name '*.log' -o -name '*.out' \) -delete

  echo "    building $(realpath --relative-to="$repo_root" "$zip_path")"
  (
    cd "$bundle_dir"
    zip -q -r "$tmp_zip" . \
      -x 'paper/*.pdf' \
      -x 'materials.zip' \
      -x './materials.zip' \
      -x 'redteam/*' \
      -x 'redteam' \
      -x '.redteam-config.yaml'
  )

  mv -f "$tmp_zip" "$zip_path"
  rmdir "$tmp_dir"
done

echo "Done."
