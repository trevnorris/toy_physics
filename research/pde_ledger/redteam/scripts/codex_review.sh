#!/usr/bin/env bash
# Read-only adversarial Codex review of the orchestrator-authored edits for ONE stage.
# Codex reads the script(s) + saved outputs + directive + paper card and emits a
# PASS/FINDINGS verdict to redteam/codex_reviews/stage_<NNN>.md.
#
#   - read-only sandbox  -> Codex physically cannot edit any file
#   - no execution       -> no `math -script`, so reviews are PARALLEL-SAFE
#
# The full codex transcript is kept at stage_<NNN>.md.raw for traceability; the
# clean verdict report is extracted to stage_<NNN>.md.
#
# Usage: codex_review.sh <stage>
set -euo pipefail

STAGE="${1:?usage: codex_review.sh <stage>}"
ROOT=/var/projects/toy_physics/research/pde_ledger
WRAPPER="$HOME/.claude/hooks/codex-chat/codex-chat"
PREAMBLE="$ROOT/redteam/prompts/codex_review.md"
OUTDIR="$ROOT/redteam/codex_reviews"
mkdir -p "$OUTDIR"
OUT="$OUTDIR/stage_${STAGE}.md"
RAW="$OUT.raw"

cd "$ROOT"

# Resolve files for this stage (gracefully tolerate a missing engine).
PY=$(ls scripts/moving_throat_pde_stage${STAGE}_*_sympy_audit.py 2>/dev/null | head -1 || true)
WL=$(ls mathematica/moving_throat_pde_stage${STAGE}_*_mathematica_audit.wl 2>/dev/null | head -1 || true)
PYO=$(ls scripts/output/moving_throat_pde_stage${STAGE}_*_sympy_audit.txt 2>/dev/null | head -1 || true)
WLO=$(ls mathematica/output/moving_throat_pde_stage${STAGE}_*_mathematica_audit.txt 2>/dev/null | head -1 || true)
DIRECTIVE="redteam/directives/stage_${STAGE}.md"
CARD="paper/stages/stage_${STAGE}.tex"

if [[ -z "$PY" && -z "$WL" ]]; then
  echo "codex_review: stage $STAGE has no .py or .wl script — nothing to review" >&2
  exit 2
fi

# Build either/or path lines (no concatenation when set).
sympy_line="${PY:+$ROOT/$PY}";   : "${sympy_line:=(absent)}"
syout_line="${PYO:+$ROOT/$PYO}"; : "${syout_line:=(absent)}"
wl_line="${WL:+$ROOT/$WL}";      : "${wl_line:=(absent)}"
wlout_line="${WLO:+$ROOT/$WLO}"; : "${wlout_line:=(absent)}"

files_block=$(cat <<EOF

---

## Files for this review (stage ${STAGE})

Read each path that is listed (skip any marked "(absent)"):

- SymPy script:        ${sympy_line}
- SymPy saved output:  ${syout_line}
- Mathematica script:  ${wl_line}
- Mathematica output:  ${wlout_line}
- Directive (what the fix had to do): $ROOT/$DIRECTIVE
- Paper stage card (what to prove):   $ROOT/$CARD

Write your verdict report for stage ${STAGE} now, following the format above exactly.
EOF
)

prompt="$(cat "$PREAMBLE")
${files_block}"

echo "codex_review: stage=$STAGE  py=${PY:+yes}${PY:-no} wl=${WL:+yes}${WL:-no} -> $OUT" >&2

# read-only sandbox, ephemeral (no persistent session), prompt via stdin.
set +e
printf '%s' "$prompt" | "$WRAPPER" -s read-only -C "$ROOT" --ephemeral > "$RAW" 2>&1
rc=$?
set -e

# Extract the final clean verdict report (the LAST `---\nstage:` frontmatter block)
# from the full codex transcript.
python3 - "$RAW" "$OUT" <<'PYEXTRACT'
import sys
raw, out = sys.argv[1], sys.argv[2]
lines = open(raw, encoding='utf-8', errors='replace').read().splitlines()
start = None
for i in range(len(lines) - 1):
    if lines[i].strip() == '---' and lines[i + 1].startswith('stage:'):
        start = i
if start is None:
    open(out, 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
    sys.exit(0)
end = len(lines)
for j in range(start, len(lines)):
    if lines[j].startswith('codex_session_id:'):
        end = j
        break
block = lines[start:end]
def junk(s):
    s = s.strip()
    return s in ('', '---') or s.startswith('tokens used') or s.replace(',', '').isdigit()
while block and junk(block[-1]):
    block.pop()
open(out, 'w', encoding='utf-8').write('\n'.join(block) + '\n')
PYEXTRACT

verdict=$(grep -iE '^verdict:' "$OUT" | head -1 || true)
echo "codex_review: stage=$STAGE exit=$rc ${verdict:-(no verdict line parsed)}" >&2
exit "$rc"
