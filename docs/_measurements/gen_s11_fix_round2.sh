#!/bin/bash
cd /var/projects/toy_physics || exit 1
S=research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py
OUT=research/pde_ledger_v3/directives/_measurements/S11_engine_fix_round2_brief.md
{
echo "# Measurements behind \`S11_engine_fix_round2_brief.md\`"
echo
echo "Commands and literal output. CLAUDE.md rule 2. Regenerate; do not transcribe."
echo "⛔ NO S11 RESULT APPEARS HERE — rule 5. The legs' reports carry values; this file carries only"
echo "the defects and their locations, because a builder reads it."
echo
echo '## 1. The path that cannot resolve (BLOCKING)'
echo '```'
echo "\$ grep -n 'REPO_ROOT' $S | head -3"
grep -n 'REPO_ROOT' "$S" | head -3
echo
echo '$ ls /var/projects/toy_physics/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md'
ls /var/projects/toy_physics/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md 2>&1
echo '$ ls research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md'
ls research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md
echo '```'
echo
echo '## 2. The publish sits inside the cell try; the end-of-run publish is uncaught'
echo '```'
echo "\$ grep -n 'write_exports(\|publish_completed_main\|except Exception as exc' $S | sed -n '1,14p'"
grep -n 'write_exports(\|publish_completed_main\|except Exception as exc' "$S" | sed -n '1,14p'
echo '```'
echo
echo '## 3. The undecided zero test and the discarded minor (REGRESSION)'
echo '```'
echo "\$ sed -n '895,901p;996,999p' $S"
sed -n '895,901p;996,999p' "$S"
echo '```'
echo
echo '## 4. skipped = declared - completed, evaluated mid-loop'
echo '```'
echo "\$ sed -n '1712,1714p' $S"
sed -n '1712,1714p' "$S"
echo '(leg probe: 18 pairs listed as skipped at MAIN-completion, none attempted;'
echo " value_kind COMPUTED_OBJECT. Probe: ~/.s11_build/leg_probes/probe_publish.py)"
echo '```'
echo
echo '## 5. The size gates, and what the SPEC requires instead'
echo '```'
echo "\$ sed -n '747,750p;769,773p' $S"
sed -n '747,750p;769,773p' "$S"
echo
echo '$ sed -n "245p;285p" research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md'
sed -n '245p;285p' research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md
echo
echo '$ grep -c UNAVAILABLE_EXACT_SYSTEM_SIZE  (token occurrences in the engine)'
grep -c 'UNAVAILABLE_EXACT_SYSTEM_SIZE' "$S"
echo '(leg measured 9/27/36/36 emissions at MAIN D2/D3/D4/D5, and that a refused'
echo ' system returns in well under a second. ⛔ The returned VALUE is withheld: rule 5.)'
echo '```'
echo
echo '## 6. The section 10 quota'
echo '```'
echo "\$ grep -n 'ISSUES\[:20\]\|ISSUES.append' $S | wc -l ; grep -n 'ISSUES\[' $S"
grep -n 'ISSUES.append' "$S" | wc -l
grep -n 'ISSUES\[' "$S"
echo '(leg measured 72 entries from MAIN alone: D2=6 D3=18 D4=24 D5=24; MAIN runs first)'
echo '```'
echo
echo '## 7-8. The third wall, and the dead function'
echo '```'
echo "\$ sed -n '774,777p' $S"
sed -n '774,777p' "$S"
echo
echo "\$ grep -n 'nullspace_basis\|generic_nullspace_vectors' $S"
grep -n 'nullspace_basis\|generic_nullspace_vectors' "$S"
echo '```'
} > "$OUT"
wc -l "$OUT"
