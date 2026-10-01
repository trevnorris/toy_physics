#!/bin/bash
# Mechanical citation lookups (verbatim retrieval / literal-match counts).
cd /var/projects/toy_physics
N=docs/native_light_em_and_vortex_throat_interpretation.md
echo "== 'director' in native_light :112-120"; sed -n 112,120p $N | grep -c -i director
echo "== 'director' anywhere in native_light (line numbers)"; grep -n -i -o "director[a-z]*" $N | head
echo "== 'C13' in V3_STEP_PLAN.md (line numbers)"; grep -n -o "C13[^|]\{0,60\}" research/pde_ledger_v3/V3_STEP_PLAN.md | head -5
echo "== V3_STEP_PLAN :1116-1126 mentions C13?"; sed -n 1116,1126p research/pde_ledger_v3/V3_STEP_PLAN.md | grep -c C13
echo "== S11c-a engine :76-105 COORDINATE declarations"; sed -n 76,105p research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py | grep -n 'COORDINATE'
echo "== S11c-b engine :552-599 vs :568-577 declared names"; sed -n 552,599p research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py | grep -o -E '"(theta_t|zeta_c|zeta_c_t|delta_v_zeta_c)"|f"zeta_c_d\{i\}"|f"delta_v_u_\{component\}(_d\{component\}|_d\{direction\})?"|"delta_v_theta", "delta_v_e_W"'
echo "== 'w' symbol declaration"; grep -n 'inherited("w"' research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py
echo "== S11c-a virtual_vertical"; sed -n 853,860p research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py
echo "== S11c-a spec J0 = 0 statement"; sed -n 368,370p research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md
echo "== S11b drain line :104"; sed -n 100,106p research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md
