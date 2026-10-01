#!/usr/bin/env bash
set -euo pipefail

engine=/tmp/s11cd_clean_review_codex/S11c_b_brane_operator_sympy_audit.copy.py
wl=/tmp/s11cd_clean_review_codex/S11c_b_brane_operator_mathematica_audit.copy.wl
exports=/tmp/s11cd_clean_review_codex/S11c_b_exports.copy.py

echo 'ENGINE_DIRECTION_AND_CASE_DECLARATIONS'
rg -n '^BRANCHES =|^DENSITY_REPS =|^DIRECTIONS =|grad_W = tuple|grad_mu = tuple' "$engine"

echo 'WOLFRAM_THREE_COORDINATE_BACKGROUND'
rg -n -m 8 'spatialCoordinates =|widthBackground\[xOne, xTwo, xThree\]|modulusBackground\[xOne, xTwo, xThree\]' "$wl"

echo 'OPERATOR_ROW_KEYS'
rg -n -m 20 'operator\["(U_BODY_BALANCE|THETA_BALANCE|E_W_BALANCE|FACE_FLUX_BOUNDARY_OPERANDS|FACE_GENERALIZED_FORCE_ROWS|MU_THETA_FACE_BINDING)"\]' "$engine"

echo 'FACE_AND_BULK_OPERAND_DECLARATIONS'
rg -n -m 20 'delta_p_\{face_name\}|bulk velocity perturbation|bulk-current perturbation|CENTER_FACE_GENERALIZED_ROW' "$engine"

echo 'ACTUAL_EXPORT_Z_DIRECTION_CARRIERS'
rg -o -m 12 'u_3_t\*w1_profile_d3|delta_v_bulk_plus_3\*w1_profile_d3|u_3\*w1_profile_d3d[123]|zeta_c_d3|e_W_d3' "$exports" | sort -u

echo 'PRIMARY_BULK_PHI_DYNAMIC_ROW_OCCURRENCES'
count=$(rg -c 'operator\["(PHI|BULK_PHI|DELTA_P)_BALANCE"\]' "$engine" || true)
echo "${count:-0}"
