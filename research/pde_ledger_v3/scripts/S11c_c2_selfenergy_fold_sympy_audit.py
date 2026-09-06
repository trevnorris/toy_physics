#!/usr/bin/env python3
"""S11c-c2 SymPy fold; authority: directives/S11c_c2_SHARED_PHYSICS.md.

All inherited physics is supplied.  The integration signs below represent
operators on arbitrary compact-support fields, with explicit computed kernels;
they are not stand-ins for a missing substitution or weak restriction.
The kernel inverse is evaluated in the rectangular retained-grade kernel algebra. See the
builder report for the scope of that representation and the normal extension.
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from collections.abc import Mapping
import hashlib
import inspect
import json
import os
import re
import resource
import sys
import time
import sympy as sp
from sympy.core.function import AppliedUndef
from sympy.core.symbol import Str

HERE = Path(__file__).resolve()
ROOT = HERE.parent.parent
sys.path.insert(0, str(HERE.parent))
from ledger_fold import load_model, check_consumer, assert_lookups_equal_manifest, assert_delta_is_minimal

IMPORT_KEYS = (
    'slab_operator', 'slab_operator_term_origins', 'mu_theta_operator',
    'closure_shape_deriv', 'face_normal', 'conormal_deriv',
    'face_measure_shape_deriv', 'face_velocity', 'relative_flux',
    'kinematic_balance', 'traction', 'face_shift', 'background_density_map',
    's11c_c1_face_response', 's11c_c1_face_response_coeffs',
    'dtn_operator', 'dtn_flat_symbol', 'dtn_kernel', 'coupling_kernel',
    'rho_m', 'rho_br', 'W_0', 'epsilon_shape', 'eta_bg', 'sigma_W',
    'W_bg', 'L_W', 'mu_R', 'omega', 'c_s0',
    'Lambda_A_0', 'Lambda_V_0', 'Lambda_X_0',
    'rho_br_bg_rho4_constant',
)
ANCHORINGS = ('LAB_HELD', 'MATERIAL_ADVECTED')
DENSITIES = ('RHO4_CONSTANT', 'RHOBR_CONSTANT')
FACES = (1, -1)
REPRESENTATION = 'DELTA_W'
EXPORT_ROOTS = frozenset((
    's11cc2ClosedSlabOperator', 's11cc2ClosedCouplingKernel',
))
NEW_DIMENSIONS = {}
COMPUTATION_LINES = {}
EMITTED_KEYS = set()
OUTPUT_ROWS = {}
COMPUTED_BINDINGS = {}

# Dimensional declaration schema of the supplied jet ansatz.  Populated below
# with the inherited S11c-b/c1 declaration dimensions (no result values).
DIMENSION_SCHEMA = {'B_rho_3': [-1, -2, 1], 'C': [-2, -2, 1], 'G_W_u': [-1, -2, 1], 'G_theta_u': [-1, -2, 1], 'L_W': [1, 0, 0], 'Lambda_A_0': [-5, 1, 1], 'Lambda_V_0': [-4, 0, 1], 'Lambda_X_0': [-4, 0, 1], 'W_0': [1, 0, 0], 'W_bg': [1, 0, 0], 'W_bg_d1': [0, 0, 0], 'W_bg_d2': [0, 0, 0], 'W_bg_d3': [0, 0, 0], 'c_s0': [1, -1, 0], 'd_w_delta_j_bulk_minus_1': [-4, -1, 1], 'd_w_delta_j_bulk_minus_2': [-4, -1, 1], 'd_w_delta_j_bulk_minus_3': [-4, -1, 1], 'd_w_delta_j_bulk_minus_4': [-4, -1, 1], 'd_w_delta_j_bulk_plus_1': [-4, -1, 1], 'd_w_delta_j_bulk_plus_2': [-4, -1, 1], 'd_w_delta_j_bulk_plus_3': [-4, -1, 1], 'd_w_delta_j_bulk_plus_4': [-4, -1, 1], 'd_w_delta_p_minus': [-3, -2, 1], 'd_w_delta_p_plus': [-3, -2, 1], 'd_w_delta_rho_4D_face_minus': [-5, 0, 1], 'd_w_delta_rho_4D_face_plus': [-5, 0, 1], 'd_w_delta_v_bulk_minus_1': [0, -1, 0], 'd_w_delta_v_bulk_minus_2': [0, -1, 0], 'd_w_delta_v_bulk_minus_3': [0, -1, 0], 'd_w_delta_v_bulk_minus_4': [0, -1, 0], 'd_w_delta_v_bulk_plus_1': [0, -1, 0], 'd_w_delta_v_bulk_plus_2': [0, -1, 0], 'd_w_delta_v_bulk_plus_3': [0, -1, 0], 'd_w_delta_v_bulk_plus_4': [0, -1, 0], 'd_w_trace_grad_f_1': [-2, 0, 0], 'd_w_trace_grad_f_2': [-2, 0, 0], 'd_w_trace_grad_f_3': [-2, 0, 0], 'd_w_trace_grad_f_4': [-2, 0, 0], 'delta_j_bulk_1': [-3, -1, 1], 'delta_j_bulk_2': [-3, -1, 1], 'delta_j_bulk_3': [-3, -1, 1], 'delta_j_bulk_4': [-3, -1, 1], 'delta_p_minus': [-2, -2, 1], 'delta_p_plus': [-2, -2, 1], 'delta_rho_4D_bulk_t': [-4, -1, 1], 'delta_rho_4D_face_minus': [-4, 0, 1], 'delta_rho_4D_face_plus': [-4, 0, 1], 'delta_v_bulk_minus_1': [1, -1, 0], 'delta_v_bulk_minus_2': [1, -1, 0], 'delta_v_bulk_minus_3': [1, -1, 0], 'delta_v_bulk_minus_4': [1, -1, 0], 'delta_v_bulk_plus_1': [1, -1, 0], 'delta_v_bulk_plus_2': [1, -1, 0], 'delta_v_bulk_plus_3': [1, -1, 0], 'delta_v_bulk_plus_4': [1, -1, 0], 'delta_v_e_W': [0, 0, 0], 'delta_v_theta': [0, 0, 0], 'delta_v_u_1': [1, 0, 0], 'delta_v_u_1_d1': [0, 0, 0], 'delta_v_u_1_d2': [0, 0, 0], 'delta_v_u_1_d3': [0, 0, 0], 'delta_v_u_2': [1, 0, 0], 'delta_v_u_2_d1': [0, 0, 0], 'delta_v_u_2_d2': [0, 0, 0], 'delta_v_u_2_d3': [0, 0, 0], 'delta_v_u_3': [1, 0, 0], 'delta_v_u_3_d1': [0, 0, 0], 'delta_v_u_3_d2': [0, 0, 0], 'delta_v_u_3_d3': [0, 0, 0], 'delta_v_zeta_c': [1, 0, 0], 'e_W': [0, 0, 0], 'e_W_bg': [0, 0, 0], 'e_W_d1': [-1, 0, 0], 'e_W_d1d1': [-2, 0, 0], 'e_W_d1d1d1': [-3, 0, 0], 'e_W_d1d1d1d1': [-4, 0, 0], 'e_W_d1d1d1d2': [-4, 0, 0], 'e_W_d1d1d1d3': [-4, 0, 0], 'e_W_d1d1d2': [-3, 0, 0], 'e_W_d1d1d2d2': [-4, 0, 0], 'e_W_d1d1d2d3': [-4, 0, 0], 'e_W_d1d1d3': [-3, 0, 0], 'e_W_d1d1d3d3': [-4, 0, 0], 'e_W_d1d2': [-2, 0, 0], 'e_W_d1d2d2': [-3, 0, 0], 'e_W_d1d2d2d2': [-4, 0, 0], 'e_W_d1d2d2d3': [-4, 0, 0], 'e_W_d1d2d3': [-3, 0, 0], 'e_W_d1d2d3d3': [-4, 0, 0], 'e_W_d1d3': [-2, 0, 0], 'e_W_d1d3d3': [-3, 0, 0], 'e_W_d1d3d3d3': [-4, 0, 0], 'e_W_d2': [-1, 0, 0], 'e_W_d2d2': [-2, 0, 0], 'e_W_d2d2d2': [-3, 0, 0], 'e_W_d2d2d2d2': [-4, 0, 0], 'e_W_d2d2d2d3': [-4, 0, 0], 'e_W_d2d2d3': [-3, 0, 0], 'e_W_d2d2d3d3': [-4, 0, 0], 'e_W_d2d3': [-2, 0, 0], 'e_W_d2d3d3': [-3, 0, 0], 'e_W_d2d3d3d3': [-4, 0, 0], 'e_W_d3': [-1, 0, 0], 'e_W_d3d3': [-2, 0, 0], 'e_W_d3d3d3': [-3, 0, 0], 'e_W_d3d3d3d3': [-4, 0, 0], 'e_W_probe': [0, 0, 0], 'e_W_probe_d1': [-1, 0, 0], 'e_W_probe_d1d1': [-2, 0, 0], 'e_W_probe_d1d1d1': [-3, 0, 0], 'e_W_probe_d1d1d1d1': [-4, 0, 0], 'e_W_probe_d1d1d1d2': [-4, 0, 0], 'e_W_probe_d1d1d1d3': [-4, 0, 0], 'e_W_probe_d1d1d2': [-3, 0, 0], 'e_W_probe_d1d1d2d2': [-4, 0, 0], 'e_W_probe_d1d1d2d3': [-4, 0, 0], 'e_W_probe_d1d1d3': [-3, 0, 0], 'e_W_probe_d1d1d3d3': [-4, 0, 0], 'e_W_probe_d1d2': [-2, 0, 0], 'e_W_probe_d1d2d2': [-3, 0, 0], 'e_W_probe_d1d2d2d2': [-4, 0, 0], 'e_W_probe_d1d2d2d3': [-4, 0, 0], 'e_W_probe_d1d2d3': [-3, 0, 0], 'e_W_probe_d1d2d3d3': [-4, 0, 0], 'e_W_probe_d1d3': [-2, 0, 0], 'e_W_probe_d1d3d3': [-3, 0, 0], 'e_W_probe_d1d3d3d3': [-4, 0, 0], 'e_W_probe_d2': [-1, 0, 0], 'e_W_probe_d2d2': [-2, 0, 0], 'e_W_probe_d2d2d2': [-3, 0, 0], 'e_W_probe_d2d2d2d2': [-4, 0, 0], 'e_W_probe_d2d2d2d3': [-4, 0, 0], 'e_W_probe_d2d2d3': [-3, 0, 0], 'e_W_probe_d2d2d3d3': [-4, 0, 0], 'e_W_probe_d2d3': [-2, 0, 0], 'e_W_probe_d2d3d3': [-3, 0, 0], 'e_W_probe_d2d3d3d3': [-4, 0, 0], 'e_W_probe_d3': [-1, 0, 0], 'e_W_probe_d3d3': [-2, 0, 0], 'e_W_probe_d3d3d3': [-3, 0, 0], 'e_W_probe_d3d3d3d3': [-4, 0, 0], 'e_W_probe_t': [0, -1, 0], 'e_W_t': [0, -1, 0], 'e_W_tt': [0, -2, 0], 'epsilon_shape': [0, 0, 0], 'eta_bg': [0, 0, 0], 'f_hold_e_W_0': [-1, -2, 1], 'f_hold_theta_0': [-1, -2, 1], 'f_hold_u_1_0': [-2, -2, 1], 'f_hold_u_2_0': [-2, -2, 1], 'f_hold_u_3_0': [-2, -2, 1], 'grad_theta_1': [-1, 0, 0], 'grad_theta_2': [-1, 0, 0], 'grad_theta_3': [-1, 0, 0], 'k': [-1, 0, 0], 'k_W': [-3, -2, 1], 'kappa_W': [-3, -2, 1], 'kappa_theta': [1, -2, 1], 'kappa_theta_W': [1, -2, 1], 'm1_profile': [0, 0, 0], 'm1_profile_d1': [0, 0, 0], 'm1_profile_d1d1': [0, 0, 0], 'm1_profile_d1d1d1': [0, 0, 0], 'm1_profile_d1d1d1d1': [0, 0, 0], 'm1_profile_d1d1d1d1d1': [0, 0, 0], 'm1_profile_d1d1d1d1d2': [0, 0, 0], 'm1_profile_d1d1d1d1d3': [0, 0, 0], 'm1_profile_d1d1d1d2': [0, 0, 0], 'm1_profile_d1d1d1d2d2': [0, 0, 0], 'm1_profile_d1d1d1d2d3': [0, 0, 0], 'm1_profile_d1d1d1d3': [0, 0, 0], 'm1_profile_d1d1d1d3d3': [0, 0, 0], 'm1_profile_d1d1d2': [0, 0, 0], 'm1_profile_d1d1d2d2': [0, 0, 0], 'm1_profile_d1d1d2d2d2': [0, 0, 0], 'm1_profile_d1d1d2d2d3': [0, 0, 0], 'm1_profile_d1d1d2d3': [0, 0, 0], 'm1_profile_d1d1d2d3d3': [0, 0, 0], 'm1_profile_d1d1d3': [0, 0, 0], 'm1_profile_d1d1d3d3': [0, 0, 0], 'm1_profile_d1d1d3d3d3': [0, 0, 0], 'm1_profile_d1d2': [0, 0, 0], 'm1_profile_d1d2d2': [0, 0, 0], 'm1_profile_d1d2d2d2': [0, 0, 0], 'm1_profile_d1d2d2d2d2': [0, 0, 0], 'm1_profile_d1d2d2d2d3': [0, 0, 0], 'm1_profile_d1d2d2d3': [0, 0, 0], 'm1_profile_d1d2d2d3d3': [0, 0, 0], 'm1_profile_d1d2d3': [0, 0, 0], 'm1_profile_d1d2d3d3': [0, 0, 0], 'm1_profile_d1d2d3d3d3': [0, 0, 0], 'm1_profile_d1d3': [0, 0, 0], 'm1_profile_d1d3d3': [0, 0, 0], 'm1_profile_d1d3d3d3': [0, 0, 0], 'm1_profile_d1d3d3d3d3': [0, 0, 0], 'm1_profile_d2': [0, 0, 0], 'm1_profile_d2d2': [0, 0, 0], 'm1_profile_d2d2d2': [0, 0, 0], 'm1_profile_d2d2d2d2': [0, 0, 0], 'm1_profile_d2d2d2d2d2': [0, 0, 0], 'm1_profile_d2d2d2d2d3': [0, 0, 0], 'm1_profile_d2d2d2d3': [0, 0, 0], 'm1_profile_d2d2d2d3d3': [0, 0, 0], 'm1_profile_d2d2d3': [0, 0, 0], 'm1_profile_d2d2d3d3': [0, 0, 0], 'm1_profile_d2d2d3d3d3': [0, 0, 0], 'm1_profile_d2d3': [0, 0, 0], 'm1_profile_d2d3d3': [0, 0, 0], 'm1_profile_d2d3d3d3': [0, 0, 0], 'm1_profile_d2d3d3d3d3': [0, 0, 0], 'm1_profile_d3': [0, 0, 0], 'm1_profile_d3d3': [0, 0, 0], 'm1_profile_d3d3d3': [0, 0, 0], 'm1_profile_d3d3d3d3': [0, 0, 0], 'm1_profile_d3d3d3d3d3': [0, 0, 0], 'mu_R': [-1, -2, 1], 'mu_R_bg': [-1, -2, 1], 'mu_R_bg_d1': [-2, -2, 1], 'mu_R_bg_d2': [-2, -2, 1], 'mu_R_bg_d3': [-2, -2, 1], 'mu_S': [-1, -2, 1], 'mu_W': [-3, 0, 1], 'mu_theta_L': [-1, -2, 1], 'mu_theta_M': [-1, -2, 1], 'omega': [0, -1, 0], 'phi_L_d1d1': [0, 0, 0], 'phi_L_d1d1d1': [-1, 0, 0], 'phi_L_d1d1d1d1': [-2, 0, 0], 'phi_L_d1d1d1d1d1': [-3, 0, 0], 'phi_L_d1d1d1d1d2': [-3, 0, 0], 'phi_L_d1d1d1d1d3': [-3, 0, 0], 'phi_L_d1d1d1d2': [-2, 0, 0], 'phi_L_d1d1d1d2d2': [-3, 0, 0], 'phi_L_d1d1d1d2d3': [-3, 0, 0], 'phi_L_d1d1d1d3': [-2, 0, 0], 'phi_L_d1d1d1d3d3': [-3, 0, 0], 'phi_L_d1d1d2': [-1, 0, 0], 'phi_L_d1d1d2d2': [-2, 0, 0], 'phi_L_d1d1d2d2d2': [-3, 0, 0], 'phi_L_d1d1d2d2d3': [-3, 0, 0], 'phi_L_d1d1d2d3': [-2, 0, 0], 'phi_L_d1d1d2d3d3': [-3, 0, 0], 'phi_L_d1d1d3': [-1, 0, 0], 'phi_L_d1d1d3d3': [-2, 0, 0], 'phi_L_d1d1d3d3d3': [-3, 0, 0], 'phi_L_d1d2': [0, 0, 0], 'phi_L_d1d2d2': [-1, 0, 0], 'phi_L_d1d2d2d2': [-2, 0, 0], 'phi_L_d1d2d2d2d2': [-3, 0, 0], 'phi_L_d1d2d2d2d3': [-3, 0, 0], 'phi_L_d1d2d2d3': [-2, 0, 0], 'phi_L_d1d2d2d3d3': [-3, 0, 0], 'phi_L_d1d2d3': [-1, 0, 0], 'phi_L_d1d2d3d3': [-2, 0, 0], 'phi_L_d1d2d3d3d3': [-3, 0, 0], 'phi_L_d1d3': [0, 0, 0], 'phi_L_d1d3d3': [-1, 0, 0], 'phi_L_d1d3d3d3': [-2, 0, 0], 'phi_L_d1d3d3d3d3': [-3, 0, 0], 'phi_L_d2d2': [0, 0, 0], 'phi_L_d2d2d2': [-1, 0, 0], 'phi_L_d2d2d2d2': [-2, 0, 0], 'phi_L_d2d2d2d2d2': [-3, 0, 0], 'phi_L_d2d2d2d2d3': [-3, 0, 0], 'phi_L_d2d2d2d3': [-2, 0, 0], 'phi_L_d2d2d2d3d3': [-3, 0, 0], 'phi_L_d2d2d3': [-1, 0, 0], 'phi_L_d2d2d3d3': [-2, 0, 0], 'phi_L_d2d2d3d3d3': [-3, 0, 0], 'phi_L_d2d3': [0, 0, 0], 'phi_L_d2d3d3': [-1, 0, 0], 'phi_L_d2d3d3d3': [-2, 0, 0], 'phi_L_d2d3d3d3d3': [-3, 0, 0], 'phi_L_d3d3': [0, 0, 0], 'phi_L_d3d3d3': [-1, 0, 0], 'phi_L_d3d3d3d3': [-2, 0, 0], 'phi_L_d3d3d3d3d3': [-3, 0, 0], 'rho_4D_bg_rho4_constant': [-4, 0, 1], 'rho_4D_bg_rhobr_constant': [-4, 0, 1], 'rho_br': [-3, 0, 1], 'rho_br_bg_rho4_constant': [-3, 0, 1], 'rho_br_bg_rhobr_constant': [-3, 0, 1], 'rho_m': [-4, 0, 1], 's11cc1_V_lab_held_minus': [1, -1, 0], 's11cc1_V_lab_held_plus': [1, -1, 0], 's11cc1_V_material_advected_minus': [1, -1, 0], 's11cc1_V_material_advected_plus': [1, -1, 0], 's11cc1_div_height_grad_operator': [0, 0, 0], 's11cc1_dtn_operator_lab_held_minus': [0, 0, 0], 's11cc1_dtn_operator_lab_held_plus': [0, 0, 0], 's11cc1_dtn_operator_material_advected_minus': [0, 0, 0], 's11cc1_dtn_operator_material_advected_plus': [0, 0, 0], 's11cc1_first_shape_impedance_operator': [0, 0, 0], 's11cc1_flat_impedance_operator': [0, 0, 0], 's11cc1_flat_normal_dtn_inverse': [0, 0, 0], 's11cc1_flat_normal_dtn_operator': [0, 0, 0], 's11cc1_height_multiplication_operator': [0, 0, 0], 's11cc1_identity_operator': [0, 0, 0], 's11cc1_k_input_1': [-1, 0, 0], 's11cc1_k_input_2': [-1, 0, 0], 's11cc1_k_input_3': [-1, 0, 0], 's11cc1_k_output_1': [-1, 0, 0], 's11cc1_k_output_2': [-1, 0, 0], 's11cc1_k_output_3': [-1, 0, 0], 's11cc1_mu_theta_lab_held_minus': [-1, -2, 1], 's11cc1_mu_theta_lab_held_plus': [-1, -2, 1], 's11cc1_mu_theta_material_advected_minus': [-1, -2, 1], 's11cc1_mu_theta_material_advected_plus': [-1, -2, 1], 's11cc1_q_out_input': [-1, 0, 0], 's11cc1_q_out_output': [-1, 0, 0], 's11cc1_response_resolvent_lab_held_minus_rho4_constant': [0, 0, 0], 's11cc1_response_resolvent_lab_held_minus_rhobr_constant': [0, 0, 0], 's11cc1_response_resolvent_lab_held_plus_rho4_constant': [0, 0, 0], 's11cc1_response_resolvent_lab_held_plus_rhobr_constant': [0, 0, 0], 's11cc1_response_resolvent_material_advected_minus_rho4_constant': [0, 0, 0], 's11cc1_response_resolvent_material_advected_minus_rhobr_constant': [0, 0, 0], 's11cc1_response_resolvent_material_advected_plus_rho4_constant': [0, 0, 0], 's11cc1_response_resolvent_material_advected_plus_rhobr_constant': [0, 0, 0], 's11cc1_w1_profile_hat_transfer': [3, 0, 0], 's11cc1_w1_profile_jet_hat_1': [3, 0, 0], 's11cc1_w1_profile_jet_hat_2': [3, 0, 0], 's11cc1_w1_profile_jet_hat_3': [3, 0, 0], 'sigma_W': [0, 0, 0], 't_hold_minus_0_1': [-2, -2, 1], 't_hold_minus_0_2': [-2, -2, 1], 't_hold_minus_0_3': [-2, -2, 1], 't_hold_minus_0_4': [-2, -2, 1], 't_hold_plus_0_1': [-2, -2, 1], 't_hold_plus_0_2': [-2, -2, 1], 't_hold_plus_0_3': [-2, -2, 1], 't_hold_plus_0_4': [-2, -2, 1], 'tau_A': [0, 1, 0], 'tau_V': [0, 1, 0], 'tau_X': [0, 1, 0], 'theta': [0, 0, 0], 'theta_d1d1': [-2, 0, 0], 'theta_d1d1d1': [-3, 0, 0], 'theta_d1d1d1d1': [-4, 0, 0], 'theta_d1d1d1d2': [-4, 0, 0], 'theta_d1d1d1d3': [-4, 0, 0], 'theta_d1d1d2': [-3, 0, 0], 'theta_d1d1d2d2': [-4, 0, 0], 'theta_d1d1d2d3': [-4, 0, 0], 'theta_d1d1d3': [-3, 0, 0], 'theta_d1d1d3d3': [-4, 0, 0], 'theta_d1d2': [-2, 0, 0], 'theta_d1d2d2': [-3, 0, 0], 'theta_d1d2d2d2': [-4, 0, 0], 'theta_d1d2d2d3': [-4, 0, 0], 'theta_d1d2d3': [-3, 0, 0], 'theta_d1d2d3d3': [-4, 0, 0], 'theta_d1d3': [-2, 0, 0], 'theta_d1d3d3': [-3, 0, 0], 'theta_d1d3d3d3': [-4, 0, 0], 'theta_d2d2': [-2, 0, 0], 'theta_d2d2d2': [-3, 0, 0], 'theta_d2d2d2d2': [-4, 0, 0], 'theta_d2d2d2d3': [-4, 0, 0], 'theta_d2d2d3': [-3, 0, 0], 'theta_d2d2d3d3': [-4, 0, 0], 'theta_d2d3': [-2, 0, 0], 'theta_d2d3d3': [-3, 0, 0], 'theta_d2d3d3d3': [-4, 0, 0], 'theta_d3d3': [-2, 0, 0], 'theta_d3d3d3': [-3, 0, 0], 'theta_d3d3d3d3': [-4, 0, 0], 'theta_probe': [0, 0, 0], 'theta_probe_d1': [-1, 0, 0], 'theta_probe_d1d1': [-2, 0, 0], 'theta_probe_d1d1d1': [-3, 0, 0], 'theta_probe_d1d1d1d1': [-4, 0, 0], 'theta_probe_d1d1d1d2': [-4, 0, 0], 'theta_probe_d1d1d1d3': [-4, 0, 0], 'theta_probe_d1d1d2': [-3, 0, 0], 'theta_probe_d1d1d2d2': [-4, 0, 0], 'theta_probe_d1d1d2d3': [-4, 0, 0], 'theta_probe_d1d1d3': [-3, 0, 0], 'theta_probe_d1d1d3d3': [-4, 0, 0], 'theta_probe_d1d2': [-2, 0, 0], 'theta_probe_d1d2d2': [-3, 0, 0], 'theta_probe_d1d2d2d2': [-4, 0, 0], 'theta_probe_d1d2d2d3': [-4, 0, 0], 'theta_probe_d1d2d3': [-3, 0, 0], 'theta_probe_d1d2d3d3': [-4, 0, 0], 'theta_probe_d1d3': [-2, 0, 0], 'theta_probe_d1d3d3': [-3, 0, 0], 'theta_probe_d1d3d3d3': [-4, 0, 0], 'theta_probe_d2': [-1, 0, 0], 'theta_probe_d2d2': [-2, 0, 0], 'theta_probe_d2d2d2': [-3, 0, 0], 'theta_probe_d2d2d2d2': [-4, 0, 0], 'theta_probe_d2d2d2d3': [-4, 0, 0], 'theta_probe_d2d2d3': [-3, 0, 0], 'theta_probe_d2d2d3d3': [-4, 0, 0], 'theta_probe_d2d3': [-2, 0, 0], 'theta_probe_d2d3d3': [-3, 0, 0], 'theta_probe_d2d3d3d3': [-4, 0, 0], 'theta_probe_d3': [-1, 0, 0], 'theta_probe_d3d3': [-2, 0, 0], 'theta_probe_d3d3d3': [-3, 0, 0], 'theta_probe_d3d3d3d3': [-4, 0, 0], 'theta_probe_t': [0, -1, 0], 'theta_t': [0, -1, 0], 'trace_grad_f_1': [-1, 0, 0], 'trace_grad_f_2': [-1, 0, 0], 'trace_grad_f_3': [-1, 0, 0], 'trace_grad_f_4': [-1, 0, 0], 'u_1': [1, 0, 0], 'u_1_d1': [0, 0, 0], 'u_1_d1d1': [-1, 0, 0], 'u_1_d1d1d1': [-2, 0, 0], 'u_1_d1d1d1d1': [-3, 0, 0], 'u_1_d1d1d1d2': [-3, 0, 0], 'u_1_d1d1d1d3': [-3, 0, 0], 'u_1_d1d1d2': [-2, 0, 0], 'u_1_d1d1d2d2': [-3, 0, 0], 'u_1_d1d1d2d3': [-3, 0, 0], 'u_1_d1d1d3': [-2, 0, 0], 'u_1_d1d1d3d3': [-3, 0, 0], 'u_1_d1d2': [-1, 0, 0], 'u_1_d1d2d2': [-2, 0, 0], 'u_1_d1d2d2d2': [-3, 0, 0], 'u_1_d1d2d2d3': [-3, 0, 0], 'u_1_d1d2d3': [-2, 0, 0], 'u_1_d1d2d3d3': [-3, 0, 0], 'u_1_d1d3': [-1, 0, 0], 'u_1_d1d3d3': [-2, 0, 0], 'u_1_d1d3d3d3': [-3, 0, 0], 'u_1_d2': [0, 0, 0], 'u_1_d2d2': [-1, 0, 0], 'u_1_d2d2d2': [-2, 0, 0], 'u_1_d2d2d2d2': [-3, 0, 0], 'u_1_d2d2d2d3': [-3, 0, 0], 'u_1_d2d2d3': [-2, 0, 0], 'u_1_d2d2d3d3': [-3, 0, 0], 'u_1_d2d3': [-1, 0, 0], 'u_1_d2d3d3': [-2, 0, 0], 'u_1_d2d3d3d3': [-3, 0, 0], 'u_1_d3': [0, 0, 0], 'u_1_d3d3': [-1, 0, 0], 'u_1_d3d3d3': [-2, 0, 0], 'u_1_d3d3d3d3': [-3, 0, 0], 'u_1_t': [1, -1, 0], 'u_1_t_d1': [0, -1, 0], 'u_1_t_d2': [0, -1, 0], 'u_1_t_d3': [0, -1, 0], 'u_1_tt': [1, -2, 0], 'u_2': [1, 0, 0], 'u_2_d1': [0, 0, 0], 'u_2_d1d1': [-1, 0, 0], 'u_2_d1d1d1': [-2, 0, 0], 'u_2_d1d1d1d1': [-3, 0, 0], 'u_2_d1d1d1d2': [-3, 0, 0], 'u_2_d1d1d1d3': [-3, 0, 0], 'u_2_d1d1d2': [-2, 0, 0], 'u_2_d1d1d2d2': [-3, 0, 0], 'u_2_d1d1d2d3': [-3, 0, 0], 'u_2_d1d1d3': [-2, 0, 0], 'u_2_d1d1d3d3': [-3, 0, 0], 'u_2_d1d2': [-1, 0, 0], 'u_2_d1d2d2': [-2, 0, 0], 'u_2_d1d2d2d2': [-3, 0, 0], 'u_2_d1d2d2d3': [-3, 0, 0], 'u_2_d1d2d3': [-2, 0, 0], 'u_2_d1d2d3d3': [-3, 0, 0], 'u_2_d1d3': [-1, 0, 0], 'u_2_d1d3d3': [-2, 0, 0], 'u_2_d1d3d3d3': [-3, 0, 0], 'u_2_d2': [0, 0, 0], 'u_2_d2d2': [-1, 0, 0], 'u_2_d2d2d2': [-2, 0, 0], 'u_2_d2d2d2d2': [-3, 0, 0], 'u_2_d2d2d2d3': [-3, 0, 0], 'u_2_d2d2d3': [-2, 0, 0], 'u_2_d2d2d3d3': [-3, 0, 0], 'u_2_d2d3': [-1, 0, 0], 'u_2_d2d3d3': [-2, 0, 0], 'u_2_d2d3d3d3': [-3, 0, 0], 'u_2_d3': [0, 0, 0], 'u_2_d3d3': [-1, 0, 0], 'u_2_d3d3d3': [-2, 0, 0], 'u_2_d3d3d3d3': [-3, 0, 0], 'u_2_t': [1, -1, 0], 'u_2_t_d1': [0, -1, 0], 'u_2_t_d2': [0, -1, 0], 'u_2_t_d3': [0, -1, 0], 'u_2_tt': [1, -2, 0], 'u_3': [1, 0, 0], 'u_3_d1': [0, 0, 0], 'u_3_d1d1': [-1, 0, 0], 'u_3_d1d1d1': [-2, 0, 0], 'u_3_d1d1d1d1': [-3, 0, 0], 'u_3_d1d1d1d2': [-3, 0, 0], 'u_3_d1d1d1d3': [-3, 0, 0], 'u_3_d1d1d2': [-2, 0, 0], 'u_3_d1d1d2d2': [-3, 0, 0], 'u_3_d1d1d2d3': [-3, 0, 0], 'u_3_d1d1d3': [-2, 0, 0], 'u_3_d1d1d3d3': [-3, 0, 0], 'u_3_d1d2': [-1, 0, 0], 'u_3_d1d2d2': [-2, 0, 0], 'u_3_d1d2d2d2': [-3, 0, 0], 'u_3_d1d2d2d3': [-3, 0, 0], 'u_3_d1d2d3': [-2, 0, 0], 'u_3_d1d2d3d3': [-3, 0, 0], 'u_3_d1d3': [-1, 0, 0], 'u_3_d1d3d3': [-2, 0, 0], 'u_3_d1d3d3d3': [-3, 0, 0], 'u_3_d2': [0, 0, 0], 'u_3_d2d2': [-1, 0, 0], 'u_3_d2d2d2': [-2, 0, 0], 'u_3_d2d2d2d2': [-3, 0, 0], 'u_3_d2d2d2d3': [-3, 0, 0], 'u_3_d2d2d3': [-2, 0, 0], 'u_3_d2d2d3d3': [-3, 0, 0], 'u_3_d2d3': [-1, 0, 0], 'u_3_d2d3d3': [-2, 0, 0], 'u_3_d2d3d3d3': [-3, 0, 0], 'u_3_d3': [0, 0, 0], 'u_3_d3d3': [-1, 0, 0], 'u_3_d3d3d3': [-2, 0, 0], 'u_3_d3d3d3d3': [-3, 0, 0], 'u_3_t': [1, -1, 0], 'u_3_t_d1': [0, -1, 0], 'u_3_t_d2': [0, -1, 0], 'u_3_t_d3': [0, -1, 0], 'u_3_tt': [1, -2, 0], 'u_L_1': [1, 0, 0], 'u_L_1_d1': [0, 0, 0], 'u_L_1_d1d1': [-1, 0, 0], 'u_L_1_d1d1d1': [-2, 0, 0], 'u_L_1_d1d1d1d1': [-3, 0, 0], 'u_L_1_d1d1d1d2': [-3, 0, 0], 'u_L_1_d1d1d1d3': [-3, 0, 0], 'u_L_1_d1d1d2': [-2, 0, 0], 'u_L_1_d1d1d2d2': [-3, 0, 0], 'u_L_1_d1d1d2d3': [-3, 0, 0], 'u_L_1_d1d1d3': [-2, 0, 0], 'u_L_1_d1d1d3d3': [-3, 0, 0], 'u_L_1_d1d2': [-1, 0, 0], 'u_L_1_d1d2d2': [-2, 0, 0], 'u_L_1_d1d2d2d2': [-3, 0, 0], 'u_L_1_d1d2d2d3': [-3, 0, 0], 'u_L_1_d1d2d3': [-2, 0, 0], 'u_L_1_d1d2d3d3': [-3, 0, 0], 'u_L_1_d1d3': [-1, 0, 0], 'u_L_1_d1d3d3': [-2, 0, 0], 'u_L_1_d1d3d3d3': [-3, 0, 0], 'u_L_1_d2': [0, 0, 0], 'u_L_1_d2d2': [-1, 0, 0], 'u_L_1_d2d2d2': [-2, 0, 0], 'u_L_1_d2d2d2d2': [-3, 0, 0], 'u_L_1_d2d2d2d3': [-3, 0, 0], 'u_L_1_d2d2d3': [-2, 0, 0], 'u_L_1_d2d2d3d3': [-3, 0, 0], 'u_L_1_d2d3': [-1, 0, 0], 'u_L_1_d2d3d3': [-2, 0, 0], 'u_L_1_d2d3d3d3': [-3, 0, 0], 'u_L_1_d3': [0, 0, 0], 'u_L_1_d3d3': [-1, 0, 0], 'u_L_1_d3d3d3': [-2, 0, 0], 'u_L_1_d3d3d3d3': [-3, 0, 0], 'u_L_1_t': [1, -1, 0], 'u_L_1_t_d1': [0, -1, 0], 'u_L_1_t_d2': [0, -1, 0], 'u_L_1_t_d3': [0, -1, 0], 'u_L_2': [1, 0, 0], 'u_L_2_d1': [0, 0, 0], 'u_L_2_d1d1': [-1, 0, 0], 'u_L_2_d1d1d1': [-2, 0, 0], 'u_L_2_d1d1d1d1': [-3, 0, 0], 'u_L_2_d1d1d1d2': [-3, 0, 0], 'u_L_2_d1d1d1d3': [-3, 0, 0], 'u_L_2_d1d1d2': [-2, 0, 0], 'u_L_2_d1d1d2d2': [-3, 0, 0], 'u_L_2_d1d1d2d3': [-3, 0, 0], 'u_L_2_d1d1d3': [-2, 0, 0], 'u_L_2_d1d1d3d3': [-3, 0, 0], 'u_L_2_d1d2': [-1, 0, 0], 'u_L_2_d1d2d2': [-2, 0, 0], 'u_L_2_d1d2d2d2': [-3, 0, 0], 'u_L_2_d1d2d2d3': [-3, 0, 0], 'u_L_2_d1d2d3': [-2, 0, 0], 'u_L_2_d1d2d3d3': [-3, 0, 0], 'u_L_2_d1d3': [-1, 0, 0], 'u_L_2_d1d3d3': [-2, 0, 0], 'u_L_2_d1d3d3d3': [-3, 0, 0], 'u_L_2_d2': [0, 0, 0], 'u_L_2_d2d2': [-1, 0, 0], 'u_L_2_d2d2d2': [-2, 0, 0], 'u_L_2_d2d2d2d2': [-3, 0, 0], 'u_L_2_d2d2d2d3': [-3, 0, 0], 'u_L_2_d2d2d3': [-2, 0, 0], 'u_L_2_d2d2d3d3': [-3, 0, 0], 'u_L_2_d2d3': [-1, 0, 0], 'u_L_2_d2d3d3': [-2, 0, 0], 'u_L_2_d2d3d3d3': [-3, 0, 0], 'u_L_2_d3': [0, 0, 0], 'u_L_2_d3d3': [-1, 0, 0], 'u_L_2_d3d3d3': [-2, 0, 0], 'u_L_2_d3d3d3d3': [-3, 0, 0], 'u_L_2_t': [1, -1, 0], 'u_L_2_t_d1': [0, -1, 0], 'u_L_2_t_d2': [0, -1, 0], 'u_L_2_t_d3': [0, -1, 0], 'u_L_3': [1, 0, 0], 'u_L_3_d1': [0, 0, 0], 'u_L_3_d1d1': [-1, 0, 0], 'u_L_3_d1d1d1': [-2, 0, 0], 'u_L_3_d1d1d1d1': [-3, 0, 0], 'u_L_3_d1d1d1d2': [-3, 0, 0], 'u_L_3_d1d1d1d3': [-3, 0, 0], 'u_L_3_d1d1d2': [-2, 0, 0], 'u_L_3_d1d1d2d2': [-3, 0, 0], 'u_L_3_d1d1d2d3': [-3, 0, 0], 'u_L_3_d1d1d3': [-2, 0, 0], 'u_L_3_d1d1d3d3': [-3, 0, 0], 'u_L_3_d1d2': [-1, 0, 0], 'u_L_3_d1d2d2': [-2, 0, 0], 'u_L_3_d1d2d2d2': [-3, 0, 0], 'u_L_3_d1d2d2d3': [-3, 0, 0], 'u_L_3_d1d2d3': [-2, 0, 0], 'u_L_3_d1d2d3d3': [-3, 0, 0], 'u_L_3_d1d3': [-1, 0, 0], 'u_L_3_d1d3d3': [-2, 0, 0], 'u_L_3_d1d3d3d3': [-3, 0, 0], 'u_L_3_d2': [0, 0, 0], 'u_L_3_d2d2': [-1, 0, 0], 'u_L_3_d2d2d2': [-2, 0, 0], 'u_L_3_d2d2d2d2': [-3, 0, 0], 'u_L_3_d2d2d2d3': [-3, 0, 0], 'u_L_3_d2d2d3': [-2, 0, 0], 'u_L_3_d2d2d3d3': [-3, 0, 0], 'u_L_3_d2d3': [-1, 0, 0], 'u_L_3_d2d3d3': [-2, 0, 0], 'u_L_3_d2d3d3d3': [-3, 0, 0], 'u_L_3_d3': [0, 0, 0], 'u_L_3_d3d3': [-1, 0, 0], 'u_L_3_d3d3d3': [-2, 0, 0], 'u_L_3_d3d3d3d3': [-3, 0, 0], 'u_L_3_t': [1, -1, 0], 'u_L_3_t_d1': [0, -1, 0], 'u_L_3_t_d2': [0, -1, 0], 'u_L_3_t_d3': [0, -1, 0], 'u_T_1': [1, 0, 0], 'u_T_1_d1': [0, 0, 0], 'u_T_1_d1d1': [-1, 0, 0], 'u_T_1_d1d1d1': [-2, 0, 0], 'u_T_1_d1d1d1d1': [-3, 0, 0], 'u_T_1_d1d1d1d2': [-3, 0, 0], 'u_T_1_d1d1d1d3': [-3, 0, 0], 'u_T_1_d1d1d2': [-2, 0, 0], 'u_T_1_d1d1d2d2': [-3, 0, 0], 'u_T_1_d1d1d2d3': [-3, 0, 0], 'u_T_1_d1d1d3': [-2, 0, 0], 'u_T_1_d1d1d3d3': [-3, 0, 0], 'u_T_1_d1d2': [-1, 0, 0], 'u_T_1_d1d2d2': [-2, 0, 0], 'u_T_1_d1d2d2d2': [-3, 0, 0], 'u_T_1_d1d2d2d3': [-3, 0, 0], 'u_T_1_d1d2d3': [-2, 0, 0], 'u_T_1_d1d2d3d3': [-3, 0, 0], 'u_T_1_d1d3': [-1, 0, 0], 'u_T_1_d1d3d3': [-2, 0, 0], 'u_T_1_d1d3d3d3': [-3, 0, 0], 'u_T_1_d2': [0, 0, 0], 'u_T_1_d2d2': [-1, 0, 0], 'u_T_1_d2d2d2': [-2, 0, 0], 'u_T_1_d2d2d2d2': [-3, 0, 0], 'u_T_1_d2d2d2d3': [-3, 0, 0], 'u_T_1_d2d2d3': [-2, 0, 0], 'u_T_1_d2d2d3d3': [-3, 0, 0], 'u_T_1_d2d3': [-1, 0, 0], 'u_T_1_d2d3d3': [-2, 0, 0], 'u_T_1_d2d3d3d3': [-3, 0, 0], 'u_T_1_d3': [0, 0, 0], 'u_T_1_d3d3': [-1, 0, 0], 'u_T_1_d3d3d3': [-2, 0, 0], 'u_T_1_d3d3d3d3': [-3, 0, 0], 'u_T_1_t': [1, -1, 0], 'u_T_1_t_d1': [0, -1, 0], 'u_T_1_t_d2': [0, -1, 0], 'u_T_1_t_d3': [0, -1, 0], 'u_T_2': [1, 0, 0], 'u_T_2_d1': [0, 0, 0], 'u_T_2_d1d1': [-1, 0, 0], 'u_T_2_d1d1d1': [-2, 0, 0], 'u_T_2_d1d1d1d1': [-3, 0, 0], 'u_T_2_d1d1d1d2': [-3, 0, 0], 'u_T_2_d1d1d1d3': [-3, 0, 0], 'u_T_2_d1d1d2': [-2, 0, 0], 'u_T_2_d1d1d2d2': [-3, 0, 0], 'u_T_2_d1d1d2d3': [-3, 0, 0], 'u_T_2_d1d1d3': [-2, 0, 0], 'u_T_2_d1d1d3d3': [-3, 0, 0], 'u_T_2_d1d2': [-1, 0, 0], 'u_T_2_d1d2d2': [-2, 0, 0], 'u_T_2_d1d2d2d2': [-3, 0, 0], 'u_T_2_d1d2d2d3': [-3, 0, 0], 'u_T_2_d1d2d3': [-2, 0, 0], 'u_T_2_d1d2d3d3': [-3, 0, 0], 'u_T_2_d1d3': [-1, 0, 0], 'u_T_2_d1d3d3': [-2, 0, 0], 'u_T_2_d1d3d3d3': [-3, 0, 0], 'u_T_2_d2': [0, 0, 0], 'u_T_2_d2d2': [-1, 0, 0], 'u_T_2_d2d2d2': [-2, 0, 0], 'u_T_2_d2d2d2d2': [-3, 0, 0], 'u_T_2_d2d2d2d3': [-3, 0, 0], 'u_T_2_d2d2d3': [-2, 0, 0], 'u_T_2_d2d2d3d3': [-3, 0, 0], 'u_T_2_d2d3': [-1, 0, 0], 'u_T_2_d2d3d3': [-2, 0, 0], 'u_T_2_d2d3d3d3': [-3, 0, 0], 'u_T_2_d3': [0, 0, 0], 'u_T_2_d3d3': [-1, 0, 0], 'u_T_2_d3d3d3': [-2, 0, 0], 'u_T_2_d3d3d3d3': [-3, 0, 0], 'u_T_2_t': [1, -1, 0], 'u_T_2_t_d1': [0, -1, 0], 'u_T_2_t_d2': [0, -1, 0], 'u_T_2_t_d3': [0, -1, 0], 'u_T_3': [1, 0, 0], 'u_T_3_d1': [0, 0, 0], 'u_T_3_d1d1': [-1, 0, 0], 'u_T_3_d1d1d1': [-2, 0, 0], 'u_T_3_d1d1d1d1': [-3, 0, 0], 'u_T_3_d1d1d1d2': [-3, 0, 0], 'u_T_3_d1d1d1d3': [-3, 0, 0], 'u_T_3_d1d1d2': [-2, 0, 0], 'u_T_3_d1d1d2d2': [-3, 0, 0], 'u_T_3_d1d1d2d3': [-3, 0, 0], 'u_T_3_d1d1d3': [-2, 0, 0], 'u_T_3_d1d1d3d3': [-3, 0, 0], 'u_T_3_d1d2': [-1, 0, 0], 'u_T_3_d1d2d2': [-2, 0, 0], 'u_T_3_d1d2d2d2': [-3, 0, 0], 'u_T_3_d1d2d2d3': [-3, 0, 0], 'u_T_3_d1d2d3': [-2, 0, 0], 'u_T_3_d1d2d3d3': [-3, 0, 0], 'u_T_3_d1d3': [-1, 0, 0], 'u_T_3_d1d3d3': [-2, 0, 0], 'u_T_3_d1d3d3d3': [-3, 0, 0], 'u_T_3_d2': [0, 0, 0], 'u_T_3_d2d2': [-1, 0, 0], 'u_T_3_d2d2d2': [-2, 0, 0], 'u_T_3_d2d2d2d2': [-3, 0, 0], 'u_T_3_d2d2d2d3': [-3, 0, 0], 'u_T_3_d2d2d3': [-2, 0, 0], 'u_T_3_d2d2d3d3': [-3, 0, 0], 'u_T_3_d2d3': [-1, 0, 0], 'u_T_3_d2d3d3': [-2, 0, 0], 'u_T_3_d2d3d3d3': [-3, 0, 0], 'u_T_3_d3': [0, 0, 0], 'u_T_3_d3d3': [-1, 0, 0], 'u_T_3_d3d3d3': [-2, 0, 0], 'u_T_3_d3d3d3d3': [-3, 0, 0], 'u_T_3_t': [1, -1, 0], 'u_T_3_t_d1': [0, -1, 0], 'u_T_3_t_d2': [0, -1, 0], 'u_T_3_t_d3': [0, -1, 0], 'v_bulk_normal_0': [1, -1, 0], 'w1_profile': [0, 0, 0], 'w1_profile_d1': [0, 0, 0], 'w1_profile_d1d1': [0, 0, 0], 'w1_profile_d1d1d1': [0, 0, 0], 'w1_profile_d1d1d1d1': [0, 0, 0], 'w1_profile_d1d1d1d1d1': [0, 0, 0], 'w1_profile_d1d1d1d1d2': [0, 0, 0], 'w1_profile_d1d1d1d1d3': [0, 0, 0], 'w1_profile_d1d1d1d2': [0, 0, 0], 'w1_profile_d1d1d1d2d2': [0, 0, 0], 'w1_profile_d1d1d1d2d3': [0, 0, 0], 'w1_profile_d1d1d1d3': [0, 0, 0], 'w1_profile_d1d1d1d3d3': [0, 0, 0], 'w1_profile_d1d1d2': [0, 0, 0], 'w1_profile_d1d1d2d2': [0, 0, 0], 'w1_profile_d1d1d2d2d2': [0, 0, 0], 'w1_profile_d1d1d2d2d3': [0, 0, 0], 'w1_profile_d1d1d2d3': [0, 0, 0], 'w1_profile_d1d1d2d3d3': [0, 0, 0], 'w1_profile_d1d1d3': [0, 0, 0], 'w1_profile_d1d1d3d3': [0, 0, 0], 'w1_profile_d1d1d3d3d3': [0, 0, 0], 'w1_profile_d1d2': [0, 0, 0], 'w1_profile_d1d2d2': [0, 0, 0], 'w1_profile_d1d2d2d2': [0, 0, 0], 'w1_profile_d1d2d2d2d2': [0, 0, 0], 'w1_profile_d1d2d2d2d3': [0, 0, 0], 'w1_profile_d1d2d2d3': [0, 0, 0], 'w1_profile_d1d2d2d3d3': [0, 0, 0], 'w1_profile_d1d2d3': [0, 0, 0], 'w1_profile_d1d2d3d3': [0, 0, 0], 'w1_profile_d1d2d3d3d3': [0, 0, 0], 'w1_profile_d1d3': [0, 0, 0], 'w1_profile_d1d3d3': [0, 0, 0], 'w1_profile_d1d3d3d3': [0, 0, 0], 'w1_profile_d1d3d3d3d3': [0, 0, 0], 'w1_profile_d2': [0, 0, 0], 'w1_profile_d2d2': [0, 0, 0], 'w1_profile_d2d2d2': [0, 0, 0], 'w1_profile_d2d2d2d2': [0, 0, 0], 'w1_profile_d2d2d2d2d2': [0, 0, 0], 'w1_profile_d2d2d2d2d3': [0, 0, 0], 'w1_profile_d2d2d2d3': [0, 0, 0], 'w1_profile_d2d2d2d3d3': [0, 0, 0], 'w1_profile_d2d2d3': [0, 0, 0], 'w1_profile_d2d2d3d3': [0, 0, 0], 'w1_profile_d2d2d3d3d3': [0, 0, 0], 'w1_profile_d2d3': [0, 0, 0], 'w1_profile_d2d3d3': [0, 0, 0], 'w1_profile_d2d3d3d3': [0, 0, 0], 'w1_profile_d2d3d3d3d3': [0, 0, 0], 'w1_profile_d3': [0, 0, 0], 'w1_profile_d3d3': [0, 0, 0], 'w1_profile_d3d3d3': [0, 0, 0], 'w1_profile_d3d3d3d3': [0, 0, 0], 'w1_profile_d3d3d3d3d3': [0, 0, 0], 'zeta_c': [1, 0, 0], 'zeta_c_d1': [0, 0, 0], 'zeta_c_d2': [0, 0, 0], 'zeta_c_d3': [0, 0, 0], 'zeta_c_t': [1, -1, 0]}
WAVE_NAMES = frozenset(['e_W', 'e_W_d1', 'e_W_d1d1', 'e_W_d1d1d1', 'e_W_d1d1d1d1', 'e_W_d1d1d1d2', 'e_W_d1d1d1d3', 'e_W_d1d1d2', 'e_W_d1d1d2d2', 'e_W_d1d1d2d3', 'e_W_d1d1d3', 'e_W_d1d1d3d3', 'e_W_d1d2', 'e_W_d1d2d2', 'e_W_d1d2d2d2', 'e_W_d1d2d2d3', 'e_W_d1d2d3', 'e_W_d1d2d3d3', 'e_W_d1d3', 'e_W_d1d3d3', 'e_W_d1d3d3d3', 'e_W_d2', 'e_W_d2d2', 'e_W_d2d2d2', 'e_W_d2d2d2d2', 'e_W_d2d2d2d3', 'e_W_d2d2d3', 'e_W_d2d2d3d3', 'e_W_d2d3', 'e_W_d2d3d3', 'e_W_d2d3d3d3', 'e_W_d3', 'e_W_d3d3', 'e_W_d3d3d3', 'e_W_d3d3d3d3', 'e_W_t', 'e_W_tt', 'grad_theta_1', 'grad_theta_2', 'grad_theta_3', 'theta', 'theta_d1d1', 'theta_d1d1d1', 'theta_d1d1d1d1', 'theta_d1d1d1d2', 'theta_d1d1d1d3', 'theta_d1d1d2', 'theta_d1d1d2d2', 'theta_d1d1d2d3', 'theta_d1d1d3', 'theta_d1d1d3d3', 'theta_d1d2', 'theta_d1d2d2', 'theta_d1d2d2d2', 'theta_d1d2d2d3', 'theta_d1d2d3', 'theta_d1d2d3d3', 'theta_d1d3', 'theta_d1d3d3', 'theta_d1d3d3d3', 'theta_d2d2', 'theta_d2d2d2', 'theta_d2d2d2d2', 'theta_d2d2d2d3', 'theta_d2d2d3', 'theta_d2d2d3d3', 'theta_d2d3', 'theta_d2d3d3', 'theta_d2d3d3d3', 'theta_d3d3', 'theta_d3d3d3', 'theta_d3d3d3d3', 'theta_t', 'u_1', 'u_1_d1', 'u_1_d1d1', 'u_1_d1d1d1', 'u_1_d1d1d1d1', 'u_1_d1d1d1d2', 'u_1_d1d1d1d3', 'u_1_d1d1d2', 'u_1_d1d1d2d2', 'u_1_d1d1d2d3', 'u_1_d1d1d3', 'u_1_d1d1d3d3', 'u_1_d1d2', 'u_1_d1d2d2', 'u_1_d1d2d2d2', 'u_1_d1d2d2d3', 'u_1_d1d2d3', 'u_1_d1d2d3d3', 'u_1_d1d3', 'u_1_d1d3d3', 'u_1_d1d3d3d3', 'u_1_d2', 'u_1_d2d2', 'u_1_d2d2d2', 'u_1_d2d2d2d2', 'u_1_d2d2d2d3', 'u_1_d2d2d3', 'u_1_d2d2d3d3', 'u_1_d2d3', 'u_1_d2d3d3', 'u_1_d2d3d3d3', 'u_1_d3', 'u_1_d3d3', 'u_1_d3d3d3', 'u_1_d3d3d3d3', 'u_1_t', 'u_1_t_d1', 'u_1_t_d2', 'u_1_t_d3', 'u_1_tt', 'u_2', 'u_2_d1', 'u_2_d1d1', 'u_2_d1d1d1', 'u_2_d1d1d1d1', 'u_2_d1d1d1d2', 'u_2_d1d1d1d3', 'u_2_d1d1d2', 'u_2_d1d1d2d2', 'u_2_d1d1d2d3', 'u_2_d1d1d3', 'u_2_d1d1d3d3', 'u_2_d1d2', 'u_2_d1d2d2', 'u_2_d1d2d2d2', 'u_2_d1d2d2d3', 'u_2_d1d2d3', 'u_2_d1d2d3d3', 'u_2_d1d3', 'u_2_d1d3d3', 'u_2_d1d3d3d3', 'u_2_d2', 'u_2_d2d2', 'u_2_d2d2d2', 'u_2_d2d2d2d2', 'u_2_d2d2d2d3', 'u_2_d2d2d3', 'u_2_d2d2d3d3', 'u_2_d2d3', 'u_2_d2d3d3', 'u_2_d2d3d3d3', 'u_2_d3', 'u_2_d3d3', 'u_2_d3d3d3', 'u_2_d3d3d3d3', 'u_2_t', 'u_2_t_d1', 'u_2_t_d2', 'u_2_t_d3', 'u_2_tt', 'u_3', 'u_3_d1', 'u_3_d1d1', 'u_3_d1d1d1', 'u_3_d1d1d1d1', 'u_3_d1d1d1d2', 'u_3_d1d1d1d3', 'u_3_d1d1d2', 'u_3_d1d1d2d2', 'u_3_d1d1d2d3', 'u_3_d1d1d3', 'u_3_d1d1d3d3', 'u_3_d1d2', 'u_3_d1d2d2', 'u_3_d1d2d2d2', 'u_3_d1d2d2d3', 'u_3_d1d2d3', 'u_3_d1d2d3d3', 'u_3_d1d3', 'u_3_d1d3d3', 'u_3_d1d3d3d3', 'u_3_d2', 'u_3_d2d2', 'u_3_d2d2d2', 'u_3_d2d2d2d2', 'u_3_d2d2d2d3', 'u_3_d2d2d3', 'u_3_d2d2d3d3', 'u_3_d2d3', 'u_3_d2d3d3', 'u_3_d2d3d3d3', 'u_3_d3', 'u_3_d3d3', 'u_3_d3d3d3', 'u_3_d3d3d3d3', 'u_3_t', 'u_3_t_d1', 'u_3_t_d2', 'u_3_t_d3', 'u_3_tt'])


def cas(value):
    if isinstance(value, str):
        return Str(value)
    if isinstance(value, bool):
        return sp.sympify(value)
    if isinstance(value, Mapping):
        return sp.Tuple(*(sp.Tuple(cas(k), cas(v)) for k, v in value.items()))
    if isinstance(value, (list, tuple, set, frozenset)):
        return sp.Tuple(*(cas(x) for x in value))
    return sp.sympify(value)


def named(value, key):
    return next(v for k, v in value if str(k) == key)


def axes(value):
    return tuple(int(x) if isinstance(x, sp.Integer) else str(x) for x in value)


def cases(value):
    return {axes(k): named(v, 'VALUE') for k, v in value}


def tree(value, fn):
    if isinstance(value, dict):
        return {k: tree(v, fn) for k, v in value.items()}
    if isinstance(value, sp.MatrixBase):
        return value.applyfunc(fn)
    if isinstance(value, (sp.Tuple, tuple, list)):
        return sp.Tuple(*(tree(v, fn) for v in value))
    if isinstance(value, Str):
        return value
    return fn(value)


def difference(a, b):
    if isinstance(a, dict):
        return {k: difference(a[k], b[k]) for k in a}
    if isinstance(a, sp.Tuple):
        return sp.Tuple(*(difference(x, y) for x, y in zip(a, b)))
    expression=a-b
    integrals=sorted(expression.atoms(sp.Integral),key=sp.default_sort_key)
    temporary={v:sp.Dummy('integralCoefficient') for v in integrals}
    expanded=sp.expand(expression.xreplace(temporary))
    return expanded.xreplace({v:k for k,v in temporary.items()})


def atom_named(value, name):
    matches = {a for a in cas(value).atoms(sp.Symbol) if a.name == name}
    if len(matches) != 1:
        raise ValueError((name, len(matches)))
    return matches.pop()


def coordinate(name, dimension):
    result = sp.Symbol(name, real=True)
    NEW_DIMENSIONS[result] = tuple(dimension)
    return result


X = tuple(coordinate(f's11cc2X{i}', (1, 0, 0)) for i in range(1, 4))
Y = tuple(coordinate(f's11cc2Y{i}', (1, 0, 0)) for i in range(1, 4))
TIME = coordinate('s11cc2Time', (0, 1, 0))
NORMAL = coordinate('s11cc2NormalCoordinate', (1, 0, 0))
MIDDLE = tuple(coordinate(f's11cc2MiddleMomentum{i}',(-1,0,0)) for i in range(1,4))
MIDDLE_Q = sp.Symbol('s11cc2MiddleNormalMomentum')
NEW_DIMENSIONS[MIDDLE_Q]=(-1,0,0)


def field(name, point, dimension):
    f = sp.Function(re.sub(r'_([a-zA-Z0-9])', lambda m:m[1].upper(), name))
    NEW_DIMENSIONS[f] = tuple(dimension)
    return f(*point, TIME)


def wave_jet(atom, point=X, sector=None):
    """The trial-field ansatz, with every spatial and temporal jet mapped.

    Solenoidal trials are curls; longitudinal trials are gradients.  This is
    the compact-support weak-restriction convention, without a Fourier projector.
    """
    name = str(atom)
    if name.startswith('grad_theta_'):
        name = 'theta_d' + name.rsplit('_', 1)[1]
    match = re.fullmatch(r'(u_[123]|theta|e_W)((?:_t{1,2})?(?:_?d[123])*)', name)
    if not match:
        raise ValueError(('wave jet', name))
    base, suffix = match.groups()
    if sector == 'TRANSVERSE' and base.startswith('u_'):
        i = int(base[-1]) - 1
        j, k = (i + 1) % 3, (i + 2) % 3
        potentials = [field(f's11cc2TrialA{n}', point, (2, 0, 0)) for n in range(3)]
        value = sp.diff(potentials[k], point[j]) - sp.diff(potentials[j], point[k])
    elif sector == 'LONGITUDINAL' and base.startswith('u_'):
        value = sp.diff(field('s11cc2TrialPhi', point, (2, 0, 0)), point[int(base[-1])-1])
    elif sector == 'THETA' and base == 'theta':
        value = field('s11cc2TrialTheta', point, (0, 0, 0))
    elif sector == 'E_W' and base == 'e_W':
        value = field('s11cc2TrialE', point, (0, 0, 0))
    elif sector is not None:
        value = sp.S.Zero
    else:
        value = field('s11cc2Field' + base, point, (1, 0, 0) if base.startswith('u_') else (0, 0, 0))
    for direction in re.findall(r'd([123])', suffix):
        value = sp.diff(value, point[int(direction)-1])
    if '_tt' in suffix:
        value = sp.diff(value, TIME, 2)
    elif '_t' in suffix:
        value = sp.diff(value, TIME)
    return value


def bind_inputs(fold):
    values = {key: fold[key]['value'] for key in IMPORT_KEYS}
    # The prefixed rows are a distinct producer, a condition the lookup guard
    # alone cannot establish.
    producer = tuple(fold[key]['step'] for key in ('s11c_c1_face_response', 's11c_c1_face_response_coeffs'))
    if producer != ('S11c-c1', 'S11c-c1'):
        raise ValueError(producer)
    return Inputs(values)


class Inputs:
    def __init__(self, values):
        self.values = values
        self.atoms = {a.name: a for value in values.values() for a in value.atoms(sp.Symbol)}
        self.eps = values['epsilon_shape']
        self.eta = values['eta_bg']
        self.sigma = values['sigma_W']
        self.slab = cases(values['slab_operator'])
        self.origins = cases(values['slab_operator_term_origins'])
        self.mu = cases(values['mu_theta_operator'])
        self.response = cases(named(values['s11c_c1_face_response'], 'CASES'))
        self.kernel = cases(values['dtn_kernel'])
        self.density = cases(values['background_density_map'])
        self.geometry = {key: cases(values[key]) for key in (
            'face_normal', 'face_velocity', 'traction', 'face_shift',
            'conormal_deriv', 'face_measure_shape_deriv', 'relative_flux',
            'kinematic_balance', 'closure_shape_deriv')}
        self.open_kernel = cases(values['coupling_kernel'])
        self.profiles = {}
        for eq in cas(self.density).atoms(sp.Equality):
            if isinstance(eq.lhs, sp.Symbol) and eq.lhs != self.sigma:
                self.profiles[eq.lhs] = eq.rhs
        self.wave = {a for a in self.atoms.values() if a.name in WAVE_NAMES}
        self.dimension_inference = self.infer_dimensions()

    def infer_dimensions(self):
        """Infer the inherited generated-coefficient units from supplied rows.

        These coefficients are created during the predecessor build, so its
        module-initialisation dimensional declarations do not include them.
        Every occurrence is retained as a dimensional consistency operand.
        """
        constraints=[]
        for case in self.slab:
            operands=[(self.mu[case][1],(-1,-2,1))]
            rows=expanded_rows(self.slab[case])
            operands += [(e,(-2,-2,1)) for e in rows['U']]
            operands += [(rows['THETA'],(-3,-1,1)),(rows['E_W'],(-1,-2,1))]
            for expression,target in operands:
                for term in sp.Add.make_args(sp.expand(expression)):
                    missing=[a for a in term.free_symbols if a.name not in DIMENSION_SCHEMA and a not in NEW_DIMENSIONS]
                    if len(missing)==1 and missing[0].name.startswith('gamma_'):
                        a=missing[0]
                        exponent=term.as_powers_dict().get(a,0)
                        if exponent:
                            other=dimension(term/a**exponent)
                            if sp.nan not in other:
                                DIMENSION_SCHEMA[a.name]=tuple((t-d)/exponent for t,d in zip(target,other))
                                dimension.cache_clear()
                                constraints.append(sp.Tuple(a,sp.ImmutableMatrix(DIMENSION_SCHEMA[a.name])))
        return sp.Tuple(*constraints)

    def a(self, name):
        return self.atoms[name]

    def physical_fields(self, value, point=X):
        replacements = {a: wave_jet(a, point) for a in value.atoms(sp.Symbol) if a in self.wave}
        return value.xreplace(replacements)

    def at_source(self, expression):
        """Translate the full source to its integration point, including density.

        Background jets are independent supplied coefficient fields here.  No
        field is replaced by its constant representative during translation.
        """
        replacements = {}
        for a in expression.atoms(sp.Symbol):
            if a in self.wave:
                replacements[a] = wave_jet(a, Y)
            elif re.match(r'(?:W_bg|mu_R_bg|[wm]1_profile)(?:$|_d)', a.name):
                base=a.name.split('_d')[0]
                name=re.sub(r'_([a-zA-Z0-9])',lambda m:m[1].upper(),'s11cc2Coefficient'+base)
                function=sp.Function(name)
                NEW_DIMENSIONS[function]=tuple(DIMENSION_SCHEMA[base])
                indices=[int(n)-1 for n in re.findall(r'd([123])',a.name)]
                value=sp.diff(function(*Y),*(Y[i] for i in indices)) if indices else function(*Y)
                if base in ('w1_profile','m1_profile'):
                    value*=self.values['L_W']**len(indices)
                replacements[a]=value
        return expression.xreplace(replacements)

    @lru_cache(maxsize=16384)
    def dx(self, expression, direction):
        if isinstance(expression,sp.Add):
            return sp.Add(*(self.dx(arg,direction) for arg in expression.args))
        if isinstance(expression,sp.Mul):
            return sp.Add(*(sp.Mul(*(self.dx(arg,direction) if i==j else arg for j,arg in enumerate(expression.args)))
                            for i in range(len(expression.args))))
        if isinstance(expression,sp.Integral):
            integrand=sp.diff(expression.function,X[direction])
            return sp.S.Zero if integrand==0 else sp.Integral(integrand,*expression.limits)
        result = sp.diff(expression, X[direction])
        for a in expression.free_symbols:
            name = a.name
            target = None
            if name in ('W_bg', 'mu_R_bg'):
                target = self.a(name + f'_d{direction+1}')
            elif name in ('w1_profile', 'm1_profile'):
                scale = self.values['W_0'] if name == 'w1_profile' else self.values['mu_R']
                bg = 'W_bg' if name == 'w1_profile' else 'mu_R_bg'
                target = self.a(bg + f'_d{direction+1}') / (scale * self.eta)
            elif re.match(r'[wm]1_profile_d', name):
                base = name.split('_d')[0]
                indices = sorted(re.findall(r'd([123])', name) + [str(direction+1)])
                new_name = base + '_' + ''.join('d' + i for i in indices)
                target_atom = self.atoms.get(new_name)
                if target_atom is None:
                    target_atom=coordinate('s11cc2'+''.join(w.title() for w in new_name.split('_')),(0,0,0))
                target = target_atom / self.values['L_W']
            elif re.fullmatch(r'(W_bg|mu_R_bg)_d[123]', name):
                bg, first = name.split('_d')
                profile = 'w1_profile' if bg == 'W_bg' else 'm1_profile'
                indices = sorted((first, str(direction+1)))
                p = self.a(profile + '_' + ''.join('d' + i for i in indices))
                scale = 1 if bg == 'W_bg' else self.values['mu_R'] / self.values['W_0']
                target = scale * self.sigma * p / self.values['L_W']
            if target is not None:
                result += sp.diff(expression, a) * target
        return result


@lru_cache(maxsize=262144)
def restricted(expression, sector):
    wave_functions=tuple(sp.Function('s11cc2Field'+s) for s in ('u1','u2','u3','theta','eW'))
    if not expression.has(*wave_functions):
        return expression
    if isinstance(expression,AppliedUndef):
        base=expression.func.__name__[len('s11cc2Field'):]
        base={'u1':'u_1','u2':'u_2','u3':'u_3','eW':'e_W'}.get(base,base)
        return wave_jet(sp.Symbol(base),expression.args[:3],sector)
    if isinstance(expression,sp.Derivative):
        return sp.diff(restricted(expression.expr,sector),*expression.variables)
    if isinstance(expression,sp.Integral):
        integrand=restricted(expression.function,sector)
        return sp.S.Zero if integrand==0 else sp.Integral(integrand,*expression.limits)
    value=expression.func(*(restricted(arg,sector) for arg in expression.args))
    return sp.expand(value,deep=False) if isinstance(value,sp.Add) else value


def extract(rows, inputs):
    u, theta, e = rows['U'], rows['THETA'], rows['E_W']
    ut = [restricted(r, 'TRANSVERSE') for r in u]
    test_theta = field('s11cc2TestTheta', X, (0, 0, 0))
    test_e = field('s11cc2TestE', X, (0, 0, 0))
    test_phi = field('s11cc2TestPhi', X, (2, 0, 0))
    test_a = [field(f's11cc2TestA{i}', X, (2, 0, 0)) for i in range(3)]
    forward = {
        'THETA': test_theta * restricted(theta, 'TRANSVERSE'),
        'E_W': test_e * restricted(e, 'TRANSVERSE'),
        'DIV_U': -test_phi * sum(inputs.dx(ut[i], i) for i in range(3)),
    }
    reverse = {}
    for target in ('THETA', 'E_W', 'LONGITUDINAL'):
        vector = [restricted(r, target) for r in u]
        curl = [inputs.dx(vector[(i+2)%3], (i+1)%3) - inputs.dx(vector[(i+1)%3], (i+2)%3) for i in range(3)]
        reverse['DIV_U' if target == 'LONGITUDINAL' else target] = sum(a*b for a, b in zip(test_a, curl))
    return {'TRANSVERSE_TO_THICKNESS': forward, 'THICKNESS_TO_TRANSVERSE': reverse}


def matrix_evaluate(expression, replacements):
    """Evaluate the inherited noncommutative inverse definition in kernel algebra."""
    if expression in replacements:
        return replacements[expression]
    if isinstance(expression, sp.Add):
        result = sp.zeros(2)
        for arg in expression.args:
            v = matrix_evaluate(arg, replacements)
            result += v if isinstance(v, sp.MatrixBase) else v * sp.eye(2)
        return result
    if isinstance(expression, sp.Mul):
        result = sp.eye(2)
        for arg in expression.args:
            result = result * matrix_evaluate(arg, replacements)
        return result
    if expression.is_commutative:
        return expression
    raise TypeError(expression)


def kernel_bridge(inputs, anchoring, face, density, overrides):
    response = inputs.response[(anchoring, face, density)]
    z = atom_named(response, f's11cc1_dtn_operator_{anchoring.lower()}_{"plus" if face == 1 else "minus"}')
    resolvent = named(response, 'RESOLVENT')
    definition = named(response, 'RESOLVENT_DEFINITION')[1]
    identity = atom_named(definition, 's11cc1_identity_operator')
    raw = inputs.kernel[(anchoring, face)]
    diagonal = named(raw, 'FLAT_DIAGONAL')
    deltas = diagonal.atoms(sp.DiracDelta)
    z0out = diagonal.xreplace({d: sp.S.One for d in deltas})
    DIMENSION_SCHEMA[z.name] = dimension(z0out)
    dimension.cache_clear()
    qo = inputs.a('s11cc1_q_out_output')
    qi = inputs.a('s11cc1_q_out_input')
    kout = [inputs.a(f's11cc1_k_output_{i}') for i in range(1, 4)]
    kin = [inputs.a(f's11cc1_k_input_{i}') for i in range(1, 4)]
    swap = dict(zip(kout, kin)) | {qo: qi}
    z0in = z0out.xreplace(swap)
    z1 = named(raw, 'FIRST_SHAPE')
    if overrides.get('FLIP_FACE_SLOPE', False) and face == 1:
        jet=inputs.a('s11cc1_w1_profile_jet_hat_1')
        normal,corrupted_normal,scale=normal_slope_control(inputs,anchoring,face)
        z1=z1.subs(jet,scale*jet)
    # On-shell reduction retains both normal-momentum legs.
    dispersion = inputs.values['omega']**2 / inputs.values['c_s0']**2
    z1 = sp.expand(z1).subs(qi**2, dispersion - sum(k*k for k in kin))
    z1 = z1.subs(qo**2, dispersion - sum(k*k for k in kout))
    z_matrix = sp.Matrix([[z0out, z1], [0, z0in]]).subs({k:v for k,v in overrides.items() if isinstance(k,sp.Basic)})
    if overrides.get('ZERO_DTN', False):
        z_matrix = z_matrix * 0
    inverse_operand = matrix_evaluate(definition.subs({k:v for k,v in overrides.items() if isinstance(k,sp.Basic)}), {identity: sp.eye(2), z: z_matrix})
    inverse = inverse_operand.upper_triangular_solve(sp.eye(2))
    response_matrix = inverse * z_matrix
    # Second scattering is required at the mixed retained grade.  A three-leg
    # triangular representation evaluates the ordered operator product, with
    # its middle leg integrated below.  There is no single-momentum division.
    transfer=fourier_profiles(inputs,z_matrix[0,1],tuple(kout),tuple(kin))
    left_map=dict(zip(kin,MIDDLE))|{qi:MIDDLE_Q}
    right_map=dict(zip(kout,MIDDLE))|{qo:MIDDLE_Q}
    z_middle=z0out.xreplace(dict(zip(kout,MIDDLE))|{qo:MIDDLE_Q}).subs({k:v for k,v in overrides.items() if isinstance(k,sp.Basic)})
    if overrides.get('ZERO_DTN', False):
        z_middle = sp.S.Zero
    z_three=sp.Matrix([[z_matrix[0,0],transfer.xreplace(left_map),0],
                       [0,z_middle,transfer.xreplace(right_map)],
                       [0,0,z_matrix[1,1]]])
    # Extract the coefficient of Z from the inherited inverse operand.
    coefficient=sp.expand(definition).coeff(z)
    coefficient=coefficient.subs({k:v for k,v in overrides.items() if isinstance(k,sp.Basic)})
    three_inverse=(sp.eye(3)+coefficient*z_three).upper_triangular_solve(sp.eye(3))
    second=(three_inverse*z_three)[0,2]
    second=sp.diff(second,inputs.eta,inputs.sigma).subs({inputs.eta:0,inputs.sigma:0})*inputs.eta*inputs.sigma
    return response, z, resolvent, z_matrix, inverse_operand, inverse, response_matrix, tuple(kout), tuple(kin), qo, second


def normal_slope_control(inputs,anchoring,face):
    normal=inputs.geometry['face_normal'][(anchoring,face,REPRESENTATION)][0]
    corrupted=sp.MutableDenseMatrix(normal)
    corrupted[0]=-corrupted[0]
    scale=sp.cancel(corrupted[0]/normal[0])
    return normal,corrupted.as_immutable(),scale


def fourier_profiles(inputs,expression,kout,kin):
    spectral={}
    for name in ('s11cc1_w1_profile_hat_transfer',*[f's11cc1_w1_profile_jet_hat_{i}' for i in range(1,4)]):
        function=sp.Function('s11cc2Fourier'+''.join(w.title() for w in name.removeprefix('s11cc1_').split('_')))
        NEW_DIMENSIONS[function]=(3,0,0)
        spectral[inputs.a(name)]=function(*(a-b for a,b in zip(kout,kin)))
    return expression.xreplace(spectral)


@lru_cache(maxsize=16)
def outgoing_spectral(inputs,kout,kin):
    spectral = {}
    for q, momenta in ((inputs.a('s11cc1_q_out_output'),kout),(inputs.a('s11cc1_q_out_input'),kin),(MIDDLE_Q,MIDDLE)):
        equation = q**2 + sum(k*k for k in momenta) - inputs.values['omega']**2/inputs.values['c_s0']**2
        roots = sp.solve(equation, q)
        root = roots[-1]
        propagating = inputs.values['omega']**2/inputs.values['c_s0']**2 - sum(k*k for k in momenta)
        solved=sp.Piecewise((sp.sign(inputs.values['omega'])*root, propagating >= 0),(root,True))
        function=sp.Function('s11cc2OutgoingNormalMomentum')
        NEW_DIMENSIONS[function]=(-1,0,0)
        representative=function(*momenta)
        COMPUTED_BINDINGS[representative]=solved
        spectral[q]=representative
    return spectral


def kernel_apply(inputs, diagonal, off_diagonal, source, kout, kin, second=sp.S.Zero):
    spectral=outgoing_spectral(inputs,tuple(kout),tuple(kin))
    diagonal = diagonal.xreplace(spectral)
    off_diagonal = fourier_profiles(inputs,off_diagonal,kout,kin).xreplace(spectral)
    second=second.xreplace(spectral)
    local_source = inputs.at_source(source)
    phase0 = sp.exp(sp.I * sum(k * (x-y) for k,x,y in zip(kout,X,Y)))
    phase1 = sp.exp(sp.I * (sum(k*x for k,x in zip(kout,X)) - sum(k*y for k,y in zip(kin,Y))))
    limits0 = tuple((v, -sp.oo, sp.oo) for v in (*kout, *Y))
    limits1 = tuple((v, -sp.oo, sp.oo) for v in (*kout, *kin, *Y))
    # c1 uses DiracDelta(k-k') without a (2*pi)^3 coefficient; the transform
    # convention has an unnormalised forward transform and normalised inverse.
    p0 = integral(phase0 * diagonal * local_source / (2*sp.pi)**3, *limits0)
    p1 = integral(phase1 * off_diagonal * local_source / (2*sp.pi)**3, *limits1)
    p2 = integral(phase1 * second * local_source / (2*sp.pi)**3,
                     *limits1,*((v,-sp.oo,sp.oo) for v in MIDDLE))
    return p0 + p1 + p2


def integral(integrand,*limits):
    return sp.S.Zero if integrand==0 else sp.Integral(integrand,*limits)


def build_face(inputs, anchoring, face, density, overrides=None, mu_override=None, velocity_override=None):
    overrides = {} if overrides is None else overrides
    numeric_overrides = {k:v for k,v in overrides.items() if isinstance(k,sp.Basic)}
    response, z, r, zm, inverse_operand, inverse, pmat, ko, ki, qo, second = kernel_bridge(inputs, anchoring, face, density, overrides)
    # Stage 0: field-vs-field density binding precedes DELTA_P composition.
    rho_live = inputs.density[(density,)][1]
    density_map = {inputs.values['rho_br_bg_rho4_constant']: rho_live}
    delta_p_source = named(response, 'DELTA_P').subs(density_map, simultaneous=True)
    # Stage 1: strip the imported operator product by computation to obtain its
    # source, while evaluating that product with the two-leg matrix above.
    source = delta_p_source.subs({r: sp.S.One, z: sp.S.One}, simultaneous=True) / inputs.eps
    label = 'plus' if face == 1 else 'minus'
    c1_mu = atom_named(source, f's11cc1_mu_theta_{anchoring.lower()}_{label}')
    c1_v = atom_named(source, f's11cc1_V_{anchoring.lower()}_{label}')
    # Stage 2: amplitude-to-amplitude identification.  The imported slab inputs
    # themselves carry epsilon; their amplitudes enter the c1 source amplitude.
    mu_amplitude = inputs.mu[(anchoring,density)][1] / inputs.eps
    velocity = inputs.geometry['face_velocity'][(anchoring,face,REPRESENTATION)] / inputs.eps
    if mu_override is not None:
        mu_amplitude = mu_override(mu_amplitude)
    if velocity_override is not None:
        velocity = velocity_override(velocity)
    stage2 = {c1_mu: mu_amplitude, c1_v: velocity}
    composed_source = sp.expand(source.subs(stage2, simultaneous=True).subs(numeric_overrides).xreplace(inputs.profiles))
    pressure = kernel_apply(inputs, pmat[0,0], pmat[0,1], composed_source, ko, ki,second)
    # Stage 3: normal jet from the outgoing continuation ansatz, differentiated
    # before evaluating at the reference face.  The output leg belongs to w.
    reference = face * inputs.values['W_0'] / 2
    extension = sp.exp(sp.I * face * qo * (NORMAL-reference))
    jet_diagonal = sp.diff(extension * pmat[0,0], NORMAL).subs(NORMAL,reference)
    jet_transfer = sp.diff(extension * pmat[0,1], NORMAL).subs(NORMAL,reference)
    jet_second=sp.diff(extension*second,NORMAL).subs(NORMAL,reference)
    normal_jet = kernel_apply(inputs, jet_diagonal, jet_transfer, composed_source, ko, ki,jet_second)
    pressure_slot = inputs.a('delta_p_' + label)
    jet_slot = inputs.a('d_w_delta_p_' + label)
    replacements = {pressure_slot: pressure, jet_slot: normal_jet}
    return replacements, {
        'DENSITY_BINDING': sp.Tuple(*[sp.Tuple(k,v) for k,v in density_map.items()]),
        'DELTA_P_SOURCE': delta_p_source,
        'IDENTIFICATIONS': sp.Tuple(*[sp.Tuple(k,v) for k,v in stage2.items()]),
        'DTN_KERNEL_MATRIX': zm,
        'RESOLVENT_INVERSE_OPERAND': inverse_operand,
        'RESOLVENT_KERNEL_MATRIX': inverse,
        'PRESSURE_KERNEL_MATRIX': pmat,
        'PRESSURE_SECOND_SCATTERING_KERNEL':second,
        'PRESSURE': pressure,
        'NORMAL_JET': normal_jet,
    }


def expanded_rows(operator):
    return {name: named(named(operator, row), 'EXPANDED') for name,row in (
        ('U','U_BODY_BALANCE'),('THETA','THETA_BALANCE'),('E_W','E_W_BALANCE'))}


def build_case(inputs, anchoring, density, *, overrides=None, routing=None, mu_override=None, velocity_override=None, input_map=None):
    """All controls enter at these imported operands, before close/extract."""
    progress_path=ROOT/'_measurements/S11c_c2_sympy_progress.json'
    progress={'anchoring':anchoring,'density':density,'routing':routing,
              'overrides':str(overrides),'mu_form':mu_override is not None,
              'velocity_control':velocity_override is not None,'started':time.time()}
    progress_path.write_text(json.dumps(progress,indent=2))
    overrides = {} if overrides is None else overrides
    numeric_overrides = {k:v for k,v in overrides.items() if isinstance(k,sp.Basic)}
    operator = inputs.slab[(anchoring,density)]
    imported = expanded_rows(operator)
    pre_map = {} if input_map is None else input_map
    source_rows = tree(imported, lambda e:e.subs(numeric_overrides).xreplace(pre_map))
    substitutions, maps, face_rows = {}, {}, {}
    for face in FACES:
        sub, fmap = build_face(inputs,anchoring,face,density,overrides,mu_override,velocity_override)
        substitutions.update(sub)
        maps[face] = fmap
        face_rows[face] = tree(source_rows, lambda e:(e-e.subs({slot:0 for slot in sub}, simultaneous=True)).subs(sub, simultaneous=True))
    # Full row closure precedes the weak restriction.  The open coupling root
    # does not occur anywhere in this construction or in the increment.
    closed = {name: tree(value, lambda e:e.subs(substitutions, simultaneous=True))
              if not (routing == 'MECHANICAL_ONLY' and name == 'THETA') else value
              for name,value in source_rows.items()}
    closed=tree(closed,lambda e:retained_shape(e,inputs))
    source_rows=tree(source_rows,lambda e:retained_shape(e,inputs))
    opened = tree(source_rows, inputs.physical_fields)
    closed = tree(closed, inputs.physical_fields)
    face_rows = {s: tree(v, lambda e:inputs.physical_fields(retained_shape(e,inputs))) for s,v in face_rows.items()}
    closed_kernel = tree(extract(closed, inputs),lambda e:retained_shape(e,inputs))
    open_kernel = tree(extract(opened, inputs),lambda e:retained_shape(e,inputs))
    increment = difference(closed_kernel, open_kernel)
    progress['finished']=time.time()
    progress_path.write_text(json.dumps(progress,indent=2))
    return {'closed':closed, 'open':opened, 'closed_kernel':closed_kernel,
            'open_kernel':open_kernel, 'increment':increment, 'faces':face_rows,
            'maps':maps, 'substitutions':substitutions, 'operator':operator}


@lru_cache(maxsize=8192)
def retained_shape(expression,inputs):
    expression=expression.xreplace(inputs.profiles)
    coefficients=shape_coefficients(expression,inputs.eta,inputs.sigma)
    return sp.Add(*(v*inputs.eta**a*inputs.sigma**b for (a,b),v in coefficients.items() if a<=1 and b<=1))


@lru_cache(maxsize=131072)
def shape_coefficients(expression,eta,sigma):
    """Evaluated rectangular Taylor algebra, commuting with explicit integrals.

    This avoids repeatedly alpha-renaming all bound momenta in SymPy's generic
    Integral differentiation.  Each kernel coefficient is still computed.
    """
    if not expression.has(eta,sigma):
        return {(0,0):expression}
    if expression==eta:
        return {(1,0):sp.S.One}
    if expression==sigma:
        return {(0,1):sp.S.One}
    if isinstance(expression,sp.Add):
        result={}
        for arg in expression.args:
            for grade,v in shape_coefficients(arg,eta,sigma).items():
                result[grade]=result.get(grade,0)+v
        return result
    if isinstance(expression,sp.Mul):
        result={(0,0):sp.S.One}
        for arg in expression.args:
            updated={}
            for (a,b),v in result.items():
                for (c,d),w in shape_coefficients(arg,eta,sigma).items():
                    if a+c<=1 and b+d<=1:
                        g=(a+c,b+d)
                        updated[g]=updated.get(g,0)+v*w
            result=updated
        return result
    if isinstance(expression,sp.Pow):
        if expression.base in (eta,sigma) and expression.exp.is_Integer:
            return {(int(expression.exp),0) if expression.base==eta else (0,int(expression.exp)):sp.S.One}
        base=shape_coefficients(expression.base,eta,sigma)
        a=sp.Dummy('seriesA');b=sp.Dummy('seriesB');c=sp.Dummy('seriesC');d=sp.Dummy('seriesD')
        h=sp.Dummy('seriesH');j=sp.Dummy('seriesJ')
        # Universal coefficient ansatz; the four coefficients are differentiated,
        # then bound to the recursively computed physical coefficients.
        ansatz=(a+h*b+j*c+h*j*d)**expression.exp
        substitutions={a:base.get((0,0),0),b:base.get((1,0),0),c:base.get((0,1),0),d:base.get((1,1),0)}
        result={}
        for r,s in ((0,0),(1,0),(0,1),(1,1)):
            derivative=sp.diff(ansatz,*([h]*r+[j]*s)) if r+s else ansatz
            result[(r,s)]=derivative.subs({h:0,j:0}).xreplace(substitutions)
        return result
    if isinstance(expression,sp.Integral):
        return {g:sp.S.Zero if v==0 else sp.Integral(v,*expression.limits)
                for g,v in shape_coefficients(expression.function,eta,sigma).items()}
    if isinstance(expression,sp.Derivative):
        return {g:sp.diff(v,*expression.variables) for g,v in shape_coefficients(expression.expr,eta,sigma).items()}
    if isinstance(expression,sp.Piecewise):
        branches=[(shape_coefficients(arg.expr,eta,sigma),arg.cond) for arg in expression.args]
        support=set().union(*(set(branch) for branch,cond in branches))
        return {g:sp.Piecewise(*[(branch.get(g,0),cond) for branch,cond in branches]) for g in support}
    raise TypeError(('shape algebra',expression.func))


@lru_cache(maxsize=131072)
def dimension(expression):
    zero = (sp.S.Zero,)*3
    if expression in NEW_DIMENSIONS:
        return NEW_DIMENSIONS[expression]
    if expression.is_Number or expression.is_number is True or isinstance(expression, (Str, sp.logic.boolalg.BooleanAtom)):
        return zero
    if isinstance(expression, sp.Symbol):
        return tuple(DIMENSION_SCHEMA.get(expression.name, (sp.nan,)*3))
    if isinstance(expression, AppliedUndef):
        return NEW_DIMENSIONS.get(expression.func, (sp.nan,)*3)
    if isinstance(expression, sp.Add):
        ds = [dimension(a) for a in expression.args if a != 0]
        return ds[0] if ds and all(d == ds[0] for d in ds) else (sp.nan,)*3
    if isinstance(expression, sp.Mul):
        ds = [dimension(a) for a in expression.args]
        return tuple(sum(d[i] for d in ds) for i in range(3))
    if isinstance(expression, sp.Pow) and expression.exp.is_number:
        return tuple(expression.exp*d for d in dimension(expression.base))
    if isinstance(expression, sp.Derivative):
        d = dimension(expression.expr)
        for v,n in expression.variable_count:
            d = tuple(a-n*b for a,b in zip(d,dimension(v)))
        return d
    if isinstance(expression, sp.Integral):
        ds = [dimension(expression.function)] + [dimension(l[0]) for l in expression.limits]
        return tuple(sum(d[i] for d in ds) for i in range(3))
    if expression.func == sp.exp:
        d = dimension(expression.args[0])
        return zero if d == zero else (sp.nan,)*3
    if expression.func == sp.DiracDelta:
        return tuple(-d for d in dimension(expression.args[0]))
    if isinstance(expression,sp.Piecewise):
        ds=[dimension(arg.expr) for arg in expression.args]
        return ds[0] if all(d==ds[0] for d in ds) else (sp.nan,)*3
    if expression.func in (sp.sign,sp.Abs):
        return zero if expression.func==sp.sign else dimension(expression.args[0])
    if isinstance(expression,sp.Equality):
        left,right=dimension(expression.lhs),dimension(expression.rhs)
        return left if left==right else (sp.nan,)*3
    return (sp.nan,)*3


def dimensions(value):
    if isinstance(value,dict):
        return {k:dimensions(v) for k,v in value.items()}
    if isinstance(value,(sp.MatrixBase,sp.Tuple,tuple,list)):
        return sp.Tuple(*(dimensions(v) for v in value))
    if isinstance(value,Str):
        return sp.zeros(3,1).as_immutable()
    return sp.ImmutableMatrix(dimension(value))


@lru_cache(maxsize=131072)
def grades(expression, epsilon, eta, sigma):
    """Structural multigrade support, including arbitrary-profile integrands.

    Denominators with shape dependence are reported by their literal degree
    rather than silently assigning a Taylor truncation to an unexpanded inverse.
    """
    z = frozenset(((0,0,0),))
    for i,knob in enumerate((epsilon,eta,sigma)):
        if expression == knob:
            g=[0,0,0];g[i]=1
            return frozenset((tuple(g),))
    if isinstance(expression, (sp.Integral, sp.Derivative)):
        return grades(expression.function if isinstance(expression,sp.Integral) else expression.expr,epsilon,eta,sigma)
    if isinstance(expression, sp.Add):
        return frozenset().union(*(grades(a,epsilon,eta,sigma) for a in expression.args))
    if isinstance(expression, sp.Mul):
        result=z
        for a in expression.args:
            result=frozenset(tuple(x+y for x,y in zip(g,h)) for g in result for h in grades(a,epsilon,eta,sigma))
        return result
    if isinstance(expression,sp.Pow) and expression.exp.is_Integer:
        support=grades(expression.base,epsilon,eta,sigma)
        if expression.exp < 0 and (0,0,0) in support and len(support)>1:
            # Analytic denominator: support of the retained rectangular Taylor
            # representative.  This does not assert a grade of the final object.
            taylor=sp.series(expression,eta,0,2).removeO()
            taylor=sp.series(taylor,sigma,0,2).removeO()
            return grades(sp.expand(taylor),epsilon,eta,sigma)
        return frozenset(tuple(expression.exp*x for x in g) for g in support)
    return z


def grade_object(value, inputs):
    result=set()
    def add(expression):
        expression=expression.xreplace(inputs.profiles)
        result.update(grades(expression,inputs.eps,inputs.eta,inputs.sigma))
        return expression
    tree(value, add)
    return sp.Tuple(*(sp.Tuple(*g) for g in sorted(result)))


@lru_cache(maxsize=4)
def profile_bindings(inputs):
    ko=tuple(inputs.a(f's11cc1_k_output_{i}') for i in range(1,4))
    ki=tuple(inputs.a(f's11cc1_k_input_{i}') for i in range(1,4))
    phase=sp.exp(-sp.I*sum((a-b)*y for a,b,y in zip(ko,ki,Y)))
    equations=[]
    pairs=[('w1_profile','s11cc1_w1_profile_hat_transfer')]+[(f'w1_profile_d{i}',f's11cc1_w1_profile_jet_hat_{i}') for i in range(1,4)]
    for local,transformed in pairs:
        local_field=inputs.at_source(inputs.a(local))
        transformed_field=fourier_profiles(inputs,inputs.a(transformed),ko,ki)
        definition=sp.Integral(phase*local_field/(2*sp.pi)**3,*((y,-sp.oo,sp.oo) for y in Y))
        equations.append(sp.Eq(transformed_field,definition,evaluate=False))
    return sp.Tuple(*equations)


def emit(quantity, value, inputs, case=(), *, key=None):
    tag = 'PY_S11CC2_' + quantity
    if case:
        tag += '_' + '_'.join(str(c).replace('-','MINUS') for c in case)
    write_key = key or 's11cc2' + ''.join(w.title().replace('_','') for w in tag.removeprefix('PY_S11CC2_').split('_'))
    if write_key in EMITTED_KEYS:
        raise ValueError(('duplicate write-key',write_key))
    EMITTED_KEYS.add(write_key)
    body={'VALUE':value,'MULTIGRADE':grade_object(value,inputs),
          'DIMENSION_L_T_M':dimensions(value)}
    if quantity in ('CLOSED_SLAB_OPERATOR','CLOSED_COUPLING_KERNEL','SELF_ENERGY_INCREMENT'):
        body['COMPUTED_BRANCH_BINDINGS']=sp.Tuple(*(sp.Eq(k,v,evaluate=False) for k,v in COMPUTED_BINDINGS.items()))
        body['FOURIER_PROFILE_BINDINGS']=profile_bindings(inputs)
    payload = cas(body)
    print(tag + ': ' + sp.srepr(payload), flush=True)
    COMPUTATION_LINES[tag] = {'emit_line':inspect.currentframe().f_back.f_lineno,'write_key':write_key}
    return payload


def control(inputs, case, baseline, name, **kwargs):
    changed=build_case(inputs,*case,**kwargs)
    emit(name+'_OPERAND',changed['increment'],inputs,case)
    emit(name+'_RESIDUAL',difference(changed['increment'],baseline['increment']),inputs,case)
    return changed


def traction_pairing(inputs, case, model, *, flip=False):
    anchoring,density=case
    covectors={}
    power=sp.S.Zero
    for face in FACES:
        value=inputs.geometry['traction'][(anchoring,face,REPRESENTATION,density)]
        if flip:
            value=-value
        traction=inputs.physical_fields(value.subs(model['substitutions'],simultaneous=True))
        # Compose the held-fixed constitutive operand in the native covector.
        mu_atom=inputs.a('mu_theta_L' if anchoring=='LAB_HELD' else 'mu_theta_M')
        mu=inputs.physical_fields(inputs.mu[case][1]/inputs.eps)
        traction=traction.subs(mu_atom,mu)
        covectors[face]=traction
        normal=inputs.geometry['face_normal'][(anchoring,face,REPRESENTATION)][0]
        v=inputs.physical_fields(inputs.geometry['face_velocity'][(anchoring,face,REPRESENTATION)]/inputs.eps)
        tangential=[wave_jet(inputs.a(f'u_{i}_t')) for i in range(1,4)]
        height_velocity=(v-sum(normal[i]*tangential[i] for i in range(3)))/normal[3]
        virtual_velocity=sp.Matrix([*tangential,height_velocity])
        area=sp.sqrt(sum(n*n for n in normal))
        power += (traction.dot(virtual_velocity))*area
    face_rows=named(model['operator'],'FACE_GENERALIZED_FORCE_ROWS')
    independent_force={k:named(face_rows,k) for k in ('U','E_W')}
    independent_force=tree(independent_force,lambda e:inputs.physical_fields(e.subs(model['substitutions'],simultaneous=True)))
    velocities=[wave_jet(inputs.a(f'u_{i}_t')) for i in range(1,4)]
    e_velocity=wave_jet(inputs.a('e_W_t'))
    slab_power=sum(r*v for r,v in zip(model['closed']['U'],velocities))+model['closed']['E_W']*e_velocity
    face_power=sum(r*v for r,v in zip(independent_force['U'],velocities))+independent_force['E_W']*e_velocity
    kinetic_stored=slab_power-face_power
    pairing=tree({'SLAB_POWER':slab_power,'KINETIC_STORED_POWER':kinetic_stored,
                  'TRACTION_POWER':power,'FACE_GENERALIZED_POWER':face_power},lambda e:retained_shape(e,inputs))
    residual=pairing['SLAB_POWER']-pairing['KINETIC_STORED_POWER']-pairing['TRACTION_POWER']
    return tree(covectors,lambda e:retained_shape(e,inputs)),pairing,residual


def publication_compact(value):
    """Factor VALUE entries while preserving reciprocal and Integral boundaries.

    Temporary atoms shield denominators and calculus objects from factor_terms;
    they are substituted back immediately, never serialized or used as CSE.
    """
    @lru_cache(maxsize=None)
    def compact(expression):
        if isinstance(expression, sp.MatrixBase):
            return expression.applyfunc(compact)
        if isinstance(expression, sp.Tuple):
            return sp.Tuple(*(compact(v) for v in expression))
        if isinstance(expression, Str) or expression.is_Atom:
            return expression
        if isinstance(expression, sp.Integral):
            return sp.Integral(compact(expression.function), *expression.limits)
        if (expression.is_Pow and expression.exp.is_negative) or isinstance(
                expression, (sp.Function, sp.Derivative)):
            return expression
        protected = {}

        @lru_cache(maxsize=None)
        def shield(node):
            if isinstance(node, sp.Integral) or isinstance(node, (sp.Function, sp.Derivative)) or (
                    node.is_Pow and node.exp.is_negative):
                token = sp.Dummy(commutative=node.is_commutative)
                protected[token] = compact(node)
                return token
            if not node.args:
                return node
            return node.func(*(shield(a) for a in node.args))

        shielded = shield(expression)
        grouped = sp.collect(shielded, list(protected), exact=True)
        result = sp.factor_terms(grouped, fraction=False).xreplace(protected)
        if result.atoms(sp.Dummy) - expression.atoms(sp.Dummy):
            raise ValueError('publication temporary escaped')
        # Factoring can expose an identically zero coefficient and erase the
        # reciprocal it multiplied. Keep that original subexpression intact.
        before_poles = {p for p in expression.atoms(sp.Pow) if p.exp.is_negative}
        after_poles = {p for p in result.atoms(sp.Pow) if p.exp.is_negative}
        return result if before_poles == after_poles else expression

    return compact(value)


def publish(inputs, fold, objects):
    if set(objects) != EXPORT_ROOTS:
        raise ValueError('publication root membership')
    compact_objects = {}
    for key, value in objects.items():
        compact_objects[key] = sp.Tuple(*(sp.Tuple(case, sp.Tuple(*(
            sp.Tuple(label, publication_compact(item) if str(label) == 'VALUE' else item)
            for label, item in payload))) for case, payload in value))
    candidates={k:{'value':v,'display':k,'value_kind':'COMPUTED_OBJECT',
                   'class':'DERIVED','step':'S11c-c2','route':'F9A_ABSENT'} for k,v in compact_objects.items()}
    all_atoms=set().union(*(v.atoms(sp.Symbol,AppliedUndef) for v in objects.values()))
    for atom in all_atoms:
        declaration=atom.func if isinstance(atom,AppliedUndef) else atom
        if declaration in NEW_DIMENSIONS:
            name=declaration.__name__ if isinstance(atom,AppliedUndef) else str(atom)
            candidates[name]={'value':declaration,'display':str(declaration),'value_kind':'COMPUTED_OBJECT',
                              'class':'COORDINATE','step':'S11c-c2','route':'F9A_ABSENT',
                              'dimension_key':name+'Dimension'}
            candidates[name+'Dimension']={'value':sp.ImmutableMatrix(NEW_DIMENSIONS[declaration]),
                'display':str(NEW_DIMENSIONS[declaration]),'value_kind':'COMPUTED_OBJECT','class':'STRUCTURAL','step':'S11c-c2','route':'F9A_ABSENT'}
    collisions=set(candidates)&set(fold)
    if collisions:
        raise ValueError(('F9 collision',sorted(collisions)))
    print('EXPORT_F9_COLLISIONS = ' + repr(sorted(collisions)), flush=True)
    combined=dict(fold)|candidates
    closure=check_consumer(combined,EXPORT_ROOTS)['closure']-set(fold)
    delta={k:candidates[k] for k in sorted(closure)}
    guard=assert_delta_is_minimal(delta,closure)
    digests={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (
        HERE,ROOT/'scripts/S11c_b_exports.py',ROOT/'scripts/S11c_c1_exports.py',
        ROOT/'directives/S11c_c2_SHARED_PHYSICS.md',ROOT/'scripts/ledger_fold.py')}
    lines=['# Generated own-rows delta; S11c-c2.','from types import MappingProxyType',
        'import sympy as sp','from sympy.core.symbol import Str',
        'from sympy.functions.elementary.piecewise import ExprCondPair',
        'def _restore(s):',"    return eval(s, {'__builtins__': {}, **vars(sp), 'Str': Str, 'ExprCondPair': ExprCondPair})",
        'IMPORT_KEYS = '+repr(IMPORT_KEYS),'BUILD_INPUT_DIGESTS = MappingProxyType('+repr(digests)+')','_LEDGER = {']
    for k,row in delta.items():
        fields=[]
        for field_name,v in row.items():
            rendered='_restore('+repr(sp.srepr(v))+')' if field_name=='value' else repr(v)
            fields.append(repr(field_name)+': '+rendered)
        lines.append(repr(k)+': {'+', '.join(fields)+'},')
    lines += ['}','LEDGER = MappingProxyType({k: MappingProxyType(v) for k,v in _LEDGER.items()})','del _LEDGER','']
    path=ROOT/'scripts/S11c_c2_exports.py'
    temporary=path.with_suffix('.py.tmp')
    temporary.write_text('\n'.join(lines))
    namespace={}
    exec(compile(temporary.read_text(),str(path),'exec'),namespace)
    restored=namespace['LEDGER']
    comparison=sp.Tuple(*(sp.sympify(delta[k]['value']==restored[k]['value']) for k in delta))
    emit('EXPORT_ROUNDTRIP',comparison,inputs)
    if any(delta[k]['value']!=restored[k]['value'] for k in delta):
        raise ValueError('serialization roundtrip')
    # Separate semantic guard: the generated module's restored values versus
    # the original emitted payload, with strict containers before arithmetic.
    evidence = []

    def record(path, check, result):
        line = f'EXPORT_SEMANTIC {path} {check} = {result}'
        print(line, flush=True)
        evidence.append(line)
        if result is not True:
            raise ValueError(('emitted/compact mismatch', path, check, result))

    def entries(value, path):
        record(path, 'mapping_tuple', isinstance(value, sp.Tuple))
        record(path, 'pair_arities', all(isinstance(p, sp.Tuple) and len(p) == 2 for p in value))
        result = dict(value)
        record(path, 'unique_keys', len(result) == len(value))
        return result

    def poles(value):
        # Exact bases AND exponents: stronger than a pole-set comparison.
        return {p for p in value.atoms(sp.Pow) if p.exp.is_negative}

    def semantic(original, decoded, path):
        if isinstance(original, Mapping):
            record(path, 'mapping_type', type(decoded) is type(original))
            record(path, 'mapping_keys', set(original) == set(decoded))
            for key in original:
                semantic(original[key], decoded[key], f'{path}/{key}')
        elif isinstance(original, sp.MatrixBase):
            record(path, 'matrix_type', type(decoded) is type(original))
            record(path, 'matrix_shape', decoded.shape == original.shape)
            for i in range(original.rows):
                for j in range(original.cols):
                    semantic(original[i, j], decoded[i, j], f'{path}[{i},{j}]')
        elif isinstance(original, (sp.Tuple, tuple, list)):
            record(path, 'tuple_type', type(decoded) is type(original))
            record(path, 'tuple_arity', len(decoded) == len(original))
            if original and all(isinstance(p, sp.Tuple) and len(p) == 2 and
                                isinstance(p[0], Str) for p in original):
                left, right = entries(original, path + '/emitted'), entries(decoded, path + '/restored')
                record(path, 'mapping_keys', set(left) == set(right))
                for key in left:
                    semantic(left[key], right[key], f'{path}/{key}')
            else:
                for i in range(len(original)):
                    semantic(original[i], decoded[i], f'{path}[{i}]')
        elif isinstance(original, Str):
            record(path, 'Str_label', type(decoded) is Str and original == decoded)
        else:
            record(path, 'algebraic_leaf', isinstance(original, sp.Expr) and isinstance(decoded, sp.Expr))
            record(path, 'reciprocal_powers_unchanged', poles(original) == poles(decoded))
            # Expand integrands with their unchanged reciprocal/calculus atoms
            # protected, then protect identical Integral atoms as in difference.
            # These local check-only dummies never enter a published value.
            atoms = (poles(original) | poles(decoded) |
                     original.atoms(sp.Function, sp.Derivative) |
                     decoded.atoms(sp.Function, sp.Derivative))
            protected = {v: sp.Dummy('publicationCoefficient') for v in atoms}
            unprotect = {v: k for k, v in protected.items()}

            @lru_cache(maxsize=None)
            def normalize_integrals(value):
                replacements = {}
                for integral in value.atoms(sp.Integral):
                    integrand = sp.expand(normalize_integrals(integral.function).xreplace(protected))
                    integrand = integrand.xreplace(unprotect)
                    replacements[integral] = (sp.S.Zero if integrand == 0 else
                                              sp.Integral(integrand, *integral.limits))
                return value.xreplace(replacements)

            expression = normalize_integrals(decoded) - normalize_integrals(original)
            integrals = {v: sp.Dummy('publicationIntegral') for v in expression.atoms(sp.Integral)}
            residual = sp.expand(expression.xreplace(integrals))
            print(f'EXPORT_SEMANTIC {path} expanded_difference = {residual}', flush=True)
            evidence.append(f'EXPORT_SEMANTIC {path} expanded_difference = {residual}')
            record(path, 'expanded_difference_is_zero', residual == 0)

    sizes = {}
    expected_cases = {cas((a, d)) for a in ANCHORINGS for d in DENSITIES}
    for key in sorted(EXPORT_ROOTS):
        original_cases = entries(objects[key], key + '/emitted')
        decoded_cases = entries(restored[key]['value'], key + '/restored')
        record(key, 'case_keys', set(original_cases) == set(decoded_cases) == expected_cases)
        sizes[key] = {}
        for case in original_cases:
            path_name = key + '/' + '/'.join(map(str, case))
            original = entries(original_cases[case], path_name + '/emitted')
            decoded = entries(decoded_cases[case], path_name + '/restored')
            record(path_name, 'payload_keys', set(original) == set(decoded))
            for label in original:
                if str(label) == 'VALUE':
                    semantic(original[label], decoded[label], path_name + '/VALUE')
                    counts = {'emitted_srepr_bytes': len(sp.srepr(original[label]).encode()),
                              'compact_srepr_bytes': len(sp.srepr(decoded[label]).encode())}
                    sizes[key]['/'.join(map(str, case))] = counts
                    print('EXPORT_VALUE_BYTES ' + path_name + ' = ' + repr(counts), flush=True)
                else:
                    record(path_name + '/' + str(label), 'metadata_exact', original[label] == decoded[label])
    restored_closure = check_consumer(dict(fold) | dict(restored), EXPORT_ROOTS)['closure'] - set(fold)
    restored_guard = assert_delta_is_minimal(restored, restored_closure)
    print('EXPORT_RESTORED_CLOSURE = ' + repr(sorted(restored_closure)), flush=True)
    print('EXPORT_RESTORED_MINIMALITY = ' + repr({k: sorted(v) for k, v in restored_guard.items()}), flush=True)
    record('delta', 'closure_unchanged', restored_closure == closure)
    temporary.replace(path)
    return {'keys':sorted(delta),'own_closure':sorted(closure),'guard':{k:sorted(v) for k,v in guard.items()},'digests':digests,
            'publication_semantic': evidence, 'publication_value_bytes': sizes}


def run():
    started=time.monotonic()
    fold,audit=load_model(str(ROOT/'scripts/S11c_b_exports.py'),str(ROOT/'scripts/S11c_c1_exports.py'))
    closure_audit=check_consumer(fold,IMPORT_KEYS)
    lookup_audit=assert_lookups_equal_manifest(bind_inputs,fold,IMPORT_KEYS)
    inputs=lookup_audit['result']
    measurements={'fold':audit,'fold_rows':len(fold),'import_keys':sorted(IMPORT_KEYS),
        'lookups':sorted(lookup_audit['lookups']),'closure':sorted(closure_audit['closure'])}
    checkpoint=ROOT/('_measurements/S11c_c2_sympy_triage_guard_evidence.json' if os.environ.get('S11CC2_PACKAGE')=='TRIAGE' else '_measurements/S11c_c2_sympy_guard_evidence.json')
    checkpoint.write_text(json.dumps(measurements,indent=2))
    emit('DIMENSION_COEFFICIENT_BINDINGS',inputs.dimension_inference,inputs)
    if os.environ.get('S11CC2_PACKAGE')=='TRIAGE':
        case=(ANCHORINGS[0],DENSITIES[0])
        model=build_case(inputs,*case)
        for name,key in (('CLOSED_SLAB_OPERATOR','closed'),('CLOSED_COUPLING_KERNEL','closed_kernel'),('SELF_ENERGY_INCREMENT','increment')):
            emit(name,model[key],inputs,case)
        return
    exports={k:{} for k in EXPORT_ROOTS}
    baselines={}
    for anchoring in ANCHORINGS:
        for density in DENSITIES:
            case=(anchoring,density)
            model=build_case(inputs,*case)
            baselines[case]=model
            for quantity,obj in (
                ('CLOSED_SLAB_OPERATOR',model['closed']),
                ('CLOSED_SLAB_OPERATOR_TERM_ORIGINS',model['faces']),
                ('CLOSED_SLAB_OPERATOR_PARITY_BLOCKS',{
                    'FACE_SUM':tree(difference(model['faces'][1],tree(model['faces'][-1],lambda e:-e)),lambda e:e/2),
                    'FACE_DIFFERENCE':tree(difference(model['faces'][1],model['faces'][-1]),lambda e:e/2)}),
                ('CLOSED_COUPLING_KERNEL',model['closed_kernel']),
                ('CLOSED_COUPLING_KERNEL_TERM_ORIGINS',{s:extract(v,inputs) for s,v in model['faces'].items()}),
                ('SELF_ENERGY_CLOSED_EXTRACT_OPERAND',model['closed_kernel']),
                ('SELF_ENERGY_OPEN_EXTRACT_OPERAND',model['open_kernel']),
                ('SELF_ENERGY_INCREMENT',model['increment']),
                ('FOLD_SYMBOL_MAP',model['maps']),
            ):
                payload=emit(quantity,obj,inputs,case)
                export_key={'CLOSED_SLAB_OPERATOR':'s11cc2ClosedSlabOperator',
                            'CLOSED_COUPLING_KERNEL':'s11cc2ClosedCouplingKernel'}.get(quantity)
                if export_key:
                    exports[export_key][case]=payload
            # Regression-only bind of the predecessor's already-extracted kernel.
            regression=named(inputs.open_kernel[case],'COMPLETE_OPERATOR_BLOCKS')
            closed_regression=tree(regression,lambda e:inputs.physical_fields(e.subs(model['substitutions'],simultaneous=True)))
            # The predecessor uses jet trial/test symbols; align those explicitly
            # in regression_coordinates before comparing to the function ansatz.
            aligned=regression_coordinates(closed_regression,inputs)
            emit('ORDERING_EXTRACT_FIRST_OPERAND',aligned,inputs,case)
            emit('ORDERING_COMMUTATOR',difference(model['closed_kernel'],aligned),inputs,case)
            covectors,pairing,residual=traction_pairing(inputs,case,model)
            emit('TRACTION_MECHANICAL_CONTRIB',covectors,inputs,case)
            emit('TRACTION_SLAB_POWER_PAIRING',pairing,inputs,case)
            emit('TRACTION_SLAB_POWER_PAIRING_RESIDUAL',residual,inputs,case)
            _,corrupt_pairing,corrupt_residual=traction_pairing(inputs,case,model,flip=True)
            emit('TRACTION_SIGN_OPERAND',corrupt_pairing,inputs,case)
            emit('TRACTION_SIGN_RESIDUAL',corrupt_residual-residual,inputs,case)
            whole_atoms=inputs.values['dtn_operator'].atoms(sp.Symbol)
            used=cas(model['increment']).atoms(sp.Symbol)&whole_atoms
            emit('DTN_WHOLEFORM_DEPENDENCE',sp.Tuple(*sorted((a for a in used if not a.is_commutative),key=str)),inputs,case)
            flat=cases(inputs.values['dtn_flat_symbol'])[(anchoring,1)]
            bridge=model['maps'][1]['DTN_KERNEL_MATRIX'][0,0]
            emit('FLAT_SYMBOL_USAGE',{'REGRESSION':flat,'KERNEL_DIAGONAL':bridge,'RESIDUAL':bridge-flat},inputs,case)
            control(inputs,case,model,'ROUTING_MECHANICAL_ONLY',routing='MECHANICAL_ONLY')
            control(inputs,case,model,'ROUTING_TRACTION_CHANNEL',overrides={inputs.values['Lambda_X_0']:0})
            control(inputs,case,model,'ZERO_DTN',overrides={'ZERO_DTN':True})
            control(inputs,case,model,'LAMBDA_A_LIMIT',overrides={inputs.values['Lambda_A_0']:0})
            control(inputs,case,model,'IMPERMEABLE_LIMIT',overrides={inputs.values['Lambda_A_0']:0,inputs.values['Lambda_V_0']:0})
            control(inputs,case,model,'UNIFORM_LIMIT',overrides={inputs.eta:0,inputs.sigma:0,inputs.values['W_bg']:inputs.values['W_0']})
            control(inputs,case,model,'MU_R_FORM',mu_override=lambda e:modulus_form(e,inputs))
            measurements['last_case']=case
            checkpoint.write_text(json.dumps(measurements,indent=2))
    for anchoring in ANCHORINGS:
        emit('DENSITY_LIVE_MINUS_FROZEN',difference(baselines[(anchoring,DENSITIES[0])]['increment'],baselines[(anchoring,DENSITIES[1])]['increment']),inputs,(anchoring,))
    for density in DENSITIES:
        lab=baselines[(ANCHORINGS[0],density)]
        material=baselines[(ANCHORINGS[1],density)]
        emit('ANCHORING_L_MINUS_M',difference(lab['increment'],material['increment']),inputs,(density,))
    measurements['export']=publish(inputs,fold,{k:cas(v) for k,v in exports.items()})
    measurements['computation_lines']=COMPUTATION_LINES
    measurements['elapsed_seconds']=time.monotonic()-started
    measurements['peak_rss_kib']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    checkpoint.write_text(json.dumps(measurements,indent=2))


def modulus_form(expression,inputs):
    # Change the profile form in the composed operand, including every jet.
    m=inputs.a('m1_profile')
    replacements={m:m**2}
    coordinates=sp.symbols('profileXi1:4')
    profile=sp.Function('profileAnsatz')(*coordinates)
    for a in expression.atoms(sp.Symbol):
        if a.name.startswith('m1_profile_d'):
            indices=[int(i)-1 for i in re.findall(r'd([123])',a.name)]
            v=sp.diff(profile**2,*(coordinates[i] for i in indices))
            jet_map={profile:m}
            for jet in v.atoms(sp.Derivative):
                suffix=''.join('d'+str(coordinates.index(q)+1) for q in jet.variables)
                jet_map[jet]=inputs.a('m1_profile_'+suffix)
            replacements[a]=v.xreplace(jet_map)
    for a,rhs in inputs.profiles.items():
        if a.name.startswith('mu_R_bg'):
            replacements[a]=rhs.xreplace(replacements)
    return expression.xreplace(replacements)


def representation_pullback(value,inputs,density):
    rho_br=inputs.density[(density,)][1]
    rho4=rho_br/inputs.values['W_bg']
    theta_shift=sum(wave_jet(inputs.a(f'u_{i+1}'))*inputs.dx(rho4,i) for i in range(3))/rho4
    e_shift=sum(wave_jet(inputs.a(f'u_{i+1}'))*inputs.a(f'W_bg_d{i+1}') for i in range(3))/inputs.values['W_0']
    theta_shift=sp.cancel(theta_shift.xreplace(inputs.profiles))
    e_shift=e_shift.xreplace(inputs.profiles)
    def pull(expression):
        replacements={}
        for atom in expression.atoms(AppliedUndef):
            point=atom.args[:3]
            translate=dict(zip(X,point))
            shift_theta=theta_shift if tuple(point)==X else inputs.at_source(theta_shift)
            shift_e=e_shift if tuple(point)==X else inputs.at_source(e_shift)
            if atom.func.__name__=='s11cc2Fieldtheta':
                replacements[atom]=atom+shift_theta.xreplace(translate)
            elif atom.func.__name__=='s11cc2FieldeW':
                replacements[atom]=atom+shift_e.xreplace(translate)
        @lru_cache(maxsize=65536)
        def compose(value):
            if not value.has(sp.Function('s11cc2Fieldtheta'),sp.Function('s11cc2FieldeW')):
                return value
            if value in replacements:
                return replacements[value]
            if isinstance(value,sp.Derivative):
                result=compose(value.expr)
                for variable in value.variables:
                    result=inputs.dx(result,X.index(variable)) if variable in X else sp.diff(result,variable)
                return result
            if isinstance(value,sp.Integral):
                return integral(compose(value.function),*value.limits)
            if not value.args:
                return value
            return value.func(*(compose(arg) for arg in value.args))
        return compose(expression)
    return tree(value,pull)


def regression_coordinates(value,inputs):
    tests={'v_theta_s11cb':field('s11cc2TestTheta',X,(0,0,0)),
           'v_e_W_s11cb':field('s11cc2TestE',X,(0,0,0)),
           'psi_L_s11cb':field('s11cc2TestPhi',X,(2,0,0))}
    tests.update({f'A_T_s11cb_{i+1}':field(f's11cc2TestA{i}',X,(2,0,0)) for i in range(3)})
    def transform(expression):
        replacements={}
        for a in expression.atoms(sp.Symbol):
            name=a.name
            if name in tests:
                replacements[a]=tests[name]
            elif name.startswith('u_T_'):
                replacements[a]=wave_jet(sp.Symbol(name.replace('u_T_','u_',1)),X,'TRANSVERSE')
            elif name.startswith('u_L_'):
                replacements[a]=wave_jet(sp.Symbol(name.replace('u_L_','u_',1)),X,'LONGITUDINAL')
            elif name.startswith('theta_probe'):
                replacements[a]=wave_jet(sp.Symbol(name.replace('theta_probe','theta',1)),X,'THETA')
            elif name.startswith('e_W_probe'):
                replacements[a]=wave_jet(sp.Symbol(name.replace('e_W_probe','e_W',1)),X,'E_W')
            elif name.startswith('phi_L_d'):
                p=field('s11cc2TrialPhi',X,(2,0,0))
                for i in re.findall(r'd([123])',name):
                    p=sp.diff(p,X[int(i)-1])
                replacements[a]=p
        return expression.xreplace(replacements)
    v=tree(value,transform)
    result={str(k):{str(n):e for n,e in row} for k,row in v}
    for direction,row in result.items():
        for target,value in row.items():
            sector='TRANSVERSE' if direction=='TRANSVERSE_TO_THICKNESS' else ('LONGITUDINAL' if target=='DIV_U' else target)
            row[target]=retained_shape(restricted(value,sector),inputs)
    return result


if __name__ == '__main__':
    run()
