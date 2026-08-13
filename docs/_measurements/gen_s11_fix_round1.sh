#!/bin/bash
cd /var/projects/toy_physics || exit 1
S=research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py
OUT=research/pde_ledger_v3/directives/_measurements/S11_engine_fix_round1_brief.md
{
echo "# Measurements behind \`S11_engine_fix_round1_brief.md\`"
echo
echo "Commands and their literal output. CLAUDE.md rule 2. Regenerate; do not transcribe."
echo
echo '## 1. The isolated D5 reproduction and its verdict'
echo '```'
echo '$ ulimit -v 3000000; python3 ~/.s11_build/repro_d5.py'
cat "$HOME/.s11_build/d5_repro_verdict.out"
echo '```'
echo
echo '## 2. The defective zero test, both call sites'
echo '```'
echo "\$ sed -n '868,870p;920p' $S"
sed -n '868,870p;920p' "$S"
echo
echo "\$ grep -n 'matrix_rank(\|iszerofunc' $S"
grep -n 'matrix_rank(\|iszerofunc' "$S"
echo '```'
echo
echo '## 3. Diagnostics are only reachable after the whole loop'
echo '```'
echo "\$ grep -n 'except Exception' $S"
grep -n 'except Exception' "$S"
echo
echo "\$ sed -n '1666,1670p;1678,1688p' $S"
sed -n '1666,1670p;1678,1688p' "$S"
echo '```'
echo
echo '## 4. MAIN D5 never completed — 188 of D4 249 suffixes absent'
echo '```'
echo '$ comm -23 <(D4 suffixes) <(D5 suffixes) | wc -l ; distinct counts ; last D5 tag'
F=$HOME/.s11_build/s11_stalled_FINAL.out
comm -23 <(grep -o '^PY_S11_MAIN_D4_[A-Z0-9_]*' "$F" | sed 's/^PY_S11_MAIN_D4_//' | sort -u) \
         <(grep -o '^PY_S11_MAIN_D5_[A-Z0-9_]*' "$F" | sed 's/^PY_S11_MAIN_D5_//' | sort -u) | wc -l
echo "D4 distinct: $(grep -o '^PY_S11_MAIN_D4_[A-Z0-9_]*' "$F" | sort -u | wc -l)"
echo "D5 distinct: $(grep -o '^PY_S11_MAIN_D5_[A-Z0-9_]*' "$F" | sort -u | wc -l)"
echo "last D5 tag: $(grep -o '^PY_S11_MAIN_D5_[A-Z0-9_]*' "$F" | tail -1)"
echo 'absent at D5, incl. the two that fill main_dim_data:'
comm -23 <(grep -o '^PY_S11_MAIN_D4_[A-Z0-9_]*' "$F" | sed 's/^PY_S11_MAIN_D4_//' | sort -u) \
         <(grep -o '^PY_S11_MAIN_D5_[A-Z0-9_]*' "$F" | sed 's/^PY_S11_MAIN_D5_//' | sort -u) \
  | grep -E '^(DIM_COEFFICIENTS|COEFFICIENT_ORDERING)$'
echo '```'
echo
echo '## 5. The script has no out/ writing'
echo '```'
echo "\$ grep -c \"out/\|OUT_DIR\" $S"
grep -c "out/\|OUT_DIR" "$S"
echo "committed out/ file: $(wc -l < research/pde_ledger_v3/scripts/out/S11_stray_longitudinal_sympy_audit.out) lines (truncated casualty)"
echo '```'
} > "$OUT"
wc -l "$OUT"
