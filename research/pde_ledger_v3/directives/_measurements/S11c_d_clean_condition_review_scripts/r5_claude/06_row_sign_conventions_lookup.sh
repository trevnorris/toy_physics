#!/bin/bash
# Mechanical lookup: are the step record's whole-row sign conventions (dated at the S11c-b close) still
# a description of the current SymPy engine?  Literal git/grep output only.
cd /var/projects/toy_physics
E=research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py
W=research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl
R=research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md
echo "== step record last commit"; git log -1 --format='%h %ad %s' --date=short -- $R
echo "== step record :112-114"; sed -n 112,114p $R
echo "== SymPy engine commits after the record"; git log --format='%h %ad %s' --date=short bcb9f7d7..HEAD -- $E
echo "== WL engine commits after the record"; git log --format='%h %ad %s' --date=short bcb9f7d7..HEAD -- $W
echo "== SymPy kinetic term in the U body row at the record's commit bcb9f7d7"
git show bcb9f7d7:$E | grep -n -E 'epsilon \* rhobr \* u_tt|u_kinetic' | head -8
echo "== SymPy kinetic term in the U body row at HEAD"
grep -n -E 'epsilon \* rhobr \* u_tt|\+ \(u_kinetic\[a\] if include_kinetic' $E | head -8
echo "== SymPy kinetic_balance_from_energy at HEAD (sign source)"; sed -n 2317,2336p $E
echo "== SymPy face multiplier at the record's commit (count of 'face_multiplier')"; git show bcb9f7d7:$E | grep -c face_multiplier
echo "== SymPy face multiplier at HEAD"; grep -n 'face_multiplier' $E
echo "== stored ACTION_TO_STORED_ROW_MULTIPLIER values in the committed exports"
grep -o "ACTION_TO_STORED_ROW_MULTIPLIER'), [^)]\{0,20\}" research/pde_ledger_v3/scripts/S11c_b_exports.py | sort | uniq -c
grep -o "ACTION_TO_STORED_ROW_MULTIPLIER[^)]\{0,30\}" research/pde_ledger_v3/scripts/S11c_b_exports.py | sort | uniq -c
echo "== WL kinetic and face-row sign lines"; grep -n -E 'kineticU *=|uRows *=' $W
