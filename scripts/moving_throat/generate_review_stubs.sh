#!/usr/bin/env bash
#
# Generate review stub files for any moving-throat stages that don't
# already have one.  Safe to re-run — never overwrites existing files.
#
# Usage:  bash scripts/moving_throat/generate_review_stubs.sh

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
NOTES_DIR="$REPO_ROOT/notes/moving_throat"
REVIEW_DIR="$NOTES_DIR/review"
SCRIPTS_DIR="$REPO_ROOT/scripts/moving_throat"

mkdir -p "$REVIEW_DIR"

get_batch() {
  local n=$1
  if (( n <= 9 ));   then echo "1|Geometry Lift & Coupling"
  elif (( n <= 19 )); then echo "2|Wall Profiles & Loading"
  elif (( n == 20 )); then echo "3|Continuum Kernel Extraction [CP]"
  elif (( n <= 29 )); then echo "4|Kernel Continuation"
  elif (( n == 30 )); then echo "5|Coherent Kernel Map [CP]"
  elif (( n <= 39 )); then echo "6|Support & Threshold Analysis"
  elif (( n == 40 )); then echo "7|Physical Parameter Placement [CP]"
  elif (( n <= 52 )); then echo "8|Operator & Gain Analysis"
  elif (( n <= 66 )); then echo "9|Wall Branch & Family-1 Geometry"
  elif (( n == 67 )); then echo "10|Full Reduced PDE Write-Up [CP]"
  elif (( n <= 72 )); then echo "11|Loading Ratio & Isotropic Verdict"
  elif (( n <= 82 )); then echo "12|Geometry Lane"
  elif (( n <= 89 )); then echo "13|Outgoing DtN"
  elif (( n <= 103 )); then echo "14|General DtN & Core Extraction"
  elif (( n <= 119 )); then echo "15|Positive Source & Mouth Dynamics"
  elif (( n <= 128 )); then echo "16|Core-to-Mouth Gain"
  elif (( n <= 140 )); then echo "17|Rigidity & Corrections"
  elif (( n <= 152 )); then echo "18|Linear Defect Transport & Final"
  else echo "NEW|Unassigned — update REVIEW_PLAN.md"
  fi
}

# Checkpoint stages
declare -A IS_CP
for cp in 020 030 040 067 082 089 100 130 150 152; do
  IS_CP[$cp]=1
done

created=0
skipped=0

for notes_file in "$NOTES_DIR"/moving_throat_pde_stage???_*.md; do
  [ -f "$notes_file" ] || continue

  basename=$(basename "$notes_file")
  num=$(echo "$basename" | grep -oP 'stage\K\d{3}')
  topic=$(echo "$basename" | sed 's|.*stage[0-9]*_||; s|\.md$||; s|_| |g')

  review_file="$REVIEW_DIR/stage_${num}_review.md"

  if [ -f "$review_file" ]; then
    skipped=$((skipped + 1))
    continue
  fi

  # Check for script
  script_file=$(ls "$SCRIPTS_DIR"/moving_throat_pde_stage${num}_*_sympy_audit*.py 2>/dev/null | head -1 || true)
  if [ -n "$script_file" ]; then
    script_line="\`scripts/moving_throat/$(basename "$script_file")\`"
  else
    script_line="None (status/consolidation stage)"
  fi

  # Batch info
  batch_info=$(get_batch $((10#$num)))
  batch_num="${batch_info%%|*}"
  batch_label="${batch_info##*|}"

  # Checkpoint note
  cp_note=""
  if [[ -n "${IS_CP[$num]:-}" ]]; then
    cp_note=$'\n**This is a CHECKPOINT stage.** Also verify cross-stage consistency (Protocol C).\n'
  fi

  cat > "$review_file" << ENDOFTEMPLATE
# Review: Stage ${num} — ${topic^}

**Batch:** ${batch_num} — ${batch_label}
**Status:** Pending
${cp_note}
## Files Under Review

- **Notes:** \`notes/moving_throat/${basename}\`
- **Script:** ${script_line}

## Review Checklist

- [ ] Equation-level correctness (signs, factors, indices, limits)
- [ ] Logical flow from prior stage(s)
- [ ] Assumptions stated and justified
- [ ] Notation consistent with prior stages
- [ ] Physical interpretation sensible
- [ ] SymPy script faithfully implements notes
- [ ] Script runs without error
- [ ] Script output matches notes claims
- [ ] No missing edge cases or branches

## Agent Reviews

<!-- Agents: append your review below this line using the template:

### Agent: [Model Name] — [Date]
**Verdict:** [PASS | MINOR | ISSUE | BLOCK]

**Notes Derivation Review:**

**Script Review:**

**Issues Found:**

**Questions:**

---
-->
ENDOFTEMPLATE

  created=$((created + 1))
  echo "  Created: stage_${num}_review.md"
done

echo ""
echo "Done. Created: ${created}, Skipped (already exist): ${skipped}"
