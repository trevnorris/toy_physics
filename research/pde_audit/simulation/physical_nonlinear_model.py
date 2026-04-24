#!/usr/bin/env python3
"""Physical nonlinear model status for the PDE audit simulation bundle.

This module intentionally does not export candidate packets.  The ledger does
not yet contain a completed strict parent-level nonlinear moving-throat branch
equation set.  What is available is an effective wall/interface closure plus
reduced BdG and Maxwell/mixed blocks.  The verifier records that boundary so the
release bundle cannot silently promote the manufactured/exporter lanes into a
physical branch realization.
"""

from __future__ import annotations

import hashlib
import json
from typing import Any, Dict, List


def stable_hash(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":"), default=float).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def build_status() -> Dict[str, Any]:
    equations: List[Dict[str, Any]] = [
        {
            "name": "promoted_confinement_variation",
            "source": "research/pde_ledger/notes/stages/moving_throat_pde_stage001_geometry_lift.md",
            "status": "available",
            "expression": "delta V_conf = - V_wall'(Sigma_0/ell_c) eta / ell_c",
            "scope": "gives a wall force/source from the promoted confinement coupling",
            "evidence_patterns": [
                "delta V_conf = - (V_wall'(Sigma_0/ell_c)/ell_c) eta",
                "S_eta^(psi) ~ - (V_wall'(Sigma_0/ell_c)/ell_c) delta rho",
            ],
        },
        {
            "name": "effective_quadratic_wall_operator",
            "source": "research/pde_ledger/notes/stages/moving_throat_pde_stage001_geometry_lift.md",
            "status": "effective_closure",
            "expression": "mu_eta q_tt - d_w(T_w q_w) + (K_eta + l(l+1) T_Omega) q = source_lm",
            "scope": "valid if S_eta^(2) is included/promoted as an effective wall action",
            "evidence_patterns": [
                "Minimal distributed geometry action (new ansatz)",
                "S_eta^(2) = (1/2) int dt dw dOmega sqrt(gamma_0)",
                "mu_eta q_{lm,tt} - partial_w( T_w partial_w q_{lm} )",
            ],
        },
        {
            "name": "axisymmetric_collective_reduction",
            "source": "research/pde_ledger/notes/stages/moving_throat_pde_stage002_breathing_reduction.md",
            "status": "available",
            "expression": "M_AB Qddot^B + K_AB Q^B = 0",
            "scope": "lowest-mode conservative reduction of the distributed wall action",
            "evidence_patterns": [
                "M_{AB}\\,\\ddot Q^B+K_{AB}\\,Q^B=0.",
                "lowest-mode truncation of the new distributed wall theory",
            ],
        },
        {
            "name": "stable_bdg_schur_kernel",
            "source": "research/pde_ledger/notes/stages/moving_throat_pde_stage003_bdg_coupling.md",
            "status": "reduced_branch_conditional",
            "expression": "D_eff(omega)=K-omega^2 M-C(Omega^2-omega^2 I)^(-1)C^T",
            "scope": "requires a stable positive-energy BdG support sector on the branch",
            "evidence_patterns": [
                "stable normal-mode reduction",
                "controlled reduced-sector move, not a full theorem of the unsolved PDE",
                "D^{\\rm eff}_0(\\omega)",
            ],
        },
        {
            "name": "maxwell_mixed_self_energy",
            "source": "research/pde_ledger/notes/stages/moving_throat_pde_stage004_maxwell_mixed_response.md",
            "status": "reduced_representative",
            "expression": "Sigma=(g_A^2 W + 2 g_A g_W R + g_W^2 A)/(A W - R^2)",
            "scope": "one-lane localized Maxwell/mixed representative, not a unique full radiation theorem",
            "evidence_patterns": [
                "Exact conservative Maxwell/mixed self-energy",
                "`= [ g_{A,l}^2 W_l",
                "positive conservative transfer factor",
            ],
        },
    ]
    blockers = [
        {
            "id": "strict_parent_dynamic_wall_not_promoted",
            "source": "research/pde_audit/notes/stage_v2_01_parent_wall_action_derivation.md",
            "reason": "The confinement-only parent action gives a force term for eta, not an autonomous moving-throat PDE. A promoted S_eta/S_Sigma action is required before a strict physical exporter is allowed.",
            "evidence_patterns": [
                "confinement-only parent action gives a force/source term",
                "does **not** give an autonomous moving-throat PDE",
                "STRICT_PARENT_DYNAMIC_WALL: FAIL unless S_eta/S_Sigma is included in S_total.",
            ],
        },
        {
            "id": "full_stationary_branch_equations_not_frozen",
            "source": "research/pde_ledger/notes/stages/moving_throat_pde_stage001_geometry_lift.md",
            "reason": "Stage 001 explicitly introduces the distributed wall action as a new ansatz and does not solve the full nonlinear moving-throat branch.",
            "evidence_patterns": [
                "It does **not** try to solve the full nonlinear moving-throat problem.",
                "Minimal distributed geometry action (new ansatz)",
            ],
        },
        {
            "id": "bdg_support_spectrum_branch_conditional",
            "source": "research/pde_ledger/notes/stages/moving_throat_pde_stage003_bdg_coupling.md",
            "reason": "The Schur kernel assumes stable positive-energy support modes; the full GNLS/BdG branch spectrum is not solved here.",
            "evidence_patterns": [
                "stable normal-mode reduction",
                "controlled reduced-sector move, not a full theorem of the unsolved PDE",
                "stable support spectrum and overlap matrices",
            ],
        },
        {
            "id": "maxwell_mixed_outgoing_not_unique_parent_theorem",
            "source": "research/pde_ledger/notes/stages/moving_throat_pde_stage004_maxwell_mixed_response.md",
            "reason": "The outgoing mixed port is a reduced representative, not a complete microscopic radiation theorem.",
            "evidence_patterns": [
                "remaining gap is the passive/outgoing normalization on the actual moving-throat branch",
                "does not yet prove that it equals the GR value on the actual moving-throat branch",
                "Do **not** jump yet to the full nonlinear PDE.",
            ],
        },
    ]
    status = {
        "schema": "pde_audit_physical_nonlinear_model_status/v1",
        "strict_parent_dynamic_wall_pass": False,
        "effective_wall_closure_available": True,
        "physical_export_permitted": False,
        "packets_emitted": False,
        "target_residuals_computed": False,
        "equation_inventory": equations,
        "blocking_reasons": blockers,
        "required_before_physical_export": [
            "Promote and freeze S_eta/S_Sigma in the declared parent action, or explicitly downgrade the lane to effective closure.",
            "Specify the coupled nonlinear stationary residual equations for R/eta, stable BdG modes, localized Maxwell modes, and mixed/outgoing ports.",
            "Add residual/Jacobian checks for the coupled physical residual.",
            "Add mesh-refinement and continuation checks on the coupled physical branch.",
            "Only then emit V2-22B-compatible physical_nonlinear_packets.",
        ],
    }
    status["status_hash"] = stable_hash(status)
    return status
