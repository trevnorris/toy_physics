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

# appended after the two legs — what they overturned
{
echo
echo '## 6. BOTH LEGS OVERTURNED ITEM 1 (verdicts: ~/.s11_build/fix1_leg_{codex,grok}.log)'
echo
echo 'Codex, measuring the matrix that reaches rank and the stronger zero predicate:'
echo '```'
echo 'MB_SHAPE (5, 5)   MB_ENTRY_OPS max 22   DET_OPS 36   DET_STRLEN 190'
echo 'ROOT 1 PRE_MAX_OPS 10  PRE_TOTAL_OPS 150  PRE_MAX_STRLEN 49'
echo 'ROOT 1 ZERO_ENTRIES_STRUCT 0   ZERO_ENTRIES_SIMPLIFY 0'
echo "RANK_EXC MemoryError  SECONDS 63.294   # with iszerofunc=factor(cancel(e))==0"
echo "ZERO_CALLS {'n': 48, 'alg_zero': 9, 'structural_false_alg_zero': 0, 'max_ops': 34661}"
echo '```'
echo '=> zero missed structural zeros; a STRONGER predicate still OOMs; entries swell 22 -> 34661 ops'
echo '   during reduction. The zero test is not the cause.'
echo
echo 'Grok, on the next wall and the dead site:'
echo '```'
echo 'all_minors(m_r, 4) at D=5 R1: ~26s per 4x4 det, 25 minors, timeout at 120s   (:906-915, :1289)'
echo 'nullspace_basis (:918-920) has NO callers; live path generic_nullspace_vectors (:924, :1245)'
echo 'MAIN D2-D4 ranks already agree with an independent rank route'
echo '```'
echo
echo 'Handler severity, where the legs DISAGREED (a finding, per rule 6):'
echo '```'
echo ':955  silent, no ISSUES   (both)        :1450 silent, no ISSUES  (both)'
echo ':1668 swallows MemoryError (both)       :756  wrong-payload risk (Codex) vs soft (Grok)'
echo '```'
} >> research/pde_ledger_v3/directives/_measurements/S11_engine_fix_round1_brief.md
