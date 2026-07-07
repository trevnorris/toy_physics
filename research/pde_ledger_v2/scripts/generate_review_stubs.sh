#!/usr/bin/env bash
#
# Generate review stub files for ledger stage notes that do not already have
# one. Safe to re-run: never overwrites existing files.
#
# Usage: bash scripts/generate_review_stubs.sh

set -euo pipefail

LEDGER_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
NOTES_DIR="$LEDGER_ROOT/notes/stages"
REVIEW_DIR="$NOTES_DIR/review"
SCRIPTS_DIR="$LEDGER_ROOT/scripts"

mkdir -p "$REVIEW_DIR"

get_batch() {
  echo "NEW|Unassigned"
}

created=0
skipped=0

for notes_file in "$NOTES_DIR"/ledger_stage???_*.md; do
  [ -f "$notes_file" ] || continue

  basename=$(basename "$notes_file")
  num=$(echo "$basename" | grep -oE 'stage[0-9]{3}' | grep -oE '[0-9]{3}')
  topic=$(echo "$basename" | sed 's|.*stage[0-9][0-9][0-9]_||; s|\.md$||; s|_| |g')

  review_file="$REVIEW_DIR/stage_${num}_review.md"

  if [ -f "$review_file" ]; then
    skipped=$((skipped + 1))
    continue
  fi

  script_file=$(ls "$SCRIPTS_DIR"/ledger_stage${num}_*_sympy_audit*.py 2>/dev/null | head -1 || true)
  if [ -n "$script_file" ]; then
    script_line="\`scripts/$(basename "$script_file")\`"
  else
    script_line="None"
  fi

  batch_info=$(get_batch)
  batch_num="${batch_info%%|*}"
  batch_label="${batch_info##*|}"

  {
    echo "# Review: Stage ${num} - ${topic^}"
    echo ""
    echo "**Batch:** ${batch_num} - ${batch_label}"
    echo "**Status:** Pending"
    echo ""
    echo "## Files Under Review"
    echo ""
    echo "- **Notes:** \`notes/stages/${basename}\`"
    echo "- **Script:** ${script_line}"
    echo ""
    echo "## Review Checklist"
    echo ""
    echo "- [ ] Equation-level correctness"
    echo "- [ ] Logical flow from prior stage(s)"
    echo "- [ ] Assumptions stated and justified"
    echo "- [ ] Notation consistent with prior stages"
    echo "- [ ] SymPy script faithfully implements notes"
    echo "- [ ] Script runs without error"
    echo "- [ ] Script output matches notes claims"
    echo "- [ ] No missing edge cases or branches"
    echo ""
    echo "## Agent Reviews"
    echo ""
    echo "<!-- Agents: append reviews below this line. -->"
  } > "$review_file"

  created=$((created + 1))
  echo "  Created: stage_${num}_review.md"
done

echo ""
echo "Done. Created: ${created}, Skipped (already exist): ${skipped}"
