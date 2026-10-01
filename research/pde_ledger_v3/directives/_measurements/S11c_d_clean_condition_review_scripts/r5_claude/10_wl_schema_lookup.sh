#!/bin/bash
# Literal-match counts: do the existing S11c-b engines carry the code paths that B's join key requires?
cd /var/projects/toy_physics
W=research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl
E=research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py
for k in U_BODY_BALANCE THETA_BALANCE E_W_BALANCE ADVECTIVE_MASS_OPERAND FACE_FLUX_BOUNDARY_OPERANDS FACE_GENERALIZED_FORCE_ROWS MU_THETA_FACE_BINDING projection_term_origins evolution_term_origins; do
  printf "%-30s PY %4s  WL %4s\n" "$k" "$(grep -c "$k" $E)" "$(grep -c "$k" $W)"
done
echo "== WL schema keys near :1139-1206"; grep -n -o -E '"(U_MOMENTUM_ROWS|THICKNESS_ROW|MASS_EVOLUTION_ROW|CENTER_FACE_GENERALIZED_ROW)"' $W | head
echo "== B join-key and blind-WL lines"
sed -n 273,275p research/pde_ledger_v3/directives/S11c_d_zinvariant_operator_blocks_directive.md
sed -n 321,322p research/pde_ledger_v3/directives/S11c_d_zinvariant_operator_blocks_directive.md
