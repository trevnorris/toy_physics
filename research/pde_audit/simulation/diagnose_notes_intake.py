#!/usr/bin/env python3
"""Intake unincorporated derivation notes into the simulation handoff.

This diagnostic is deliberately non-generative.  It reads selected notes,
checks source anchors, and records what those notes imply for the next
simulation step.  It does not produce solver packets, rank candidates, or use
target residuals to mutate any branch family.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Dict, List, Mapping


SIM_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SIM_DIR.parents[2]
OUTPUT_DIR = SIM_DIR / "output"


ANCHORS = [
    {
        "id": "fivepn_computational_handoff",
        "source": "notes/5pn/5pn_stage354_355_computational_handoff.md",
        "role": "Defines the remaining actual-branch computation after reduced algebra.",
        "patterns": [
            "The reduced theorem side is now essentially closed.",
            "The current chain has already reduced the problem to four finish-line conditions on the actual PDE-selected branch:",
            "`N_Q = 1`",
            "What is still missing is the **actual PDE-selected orbit-lock point**",
            "Do **not** manually project it back onto the canonical outgoing branch just to make the theorem pass.",
            "Because `zeta` drops out of the orbit packet exactly",
        ],
    },
    {
        "id": "fivepn_final_packet_shape",
        "source": "notes/5pn/5pn_final_wrap_up.md",
        "role": "Freezes Packet A/Packet B as the next actual-branch payloads.",
        "patterns": [
            "Compute **Packet A** from the actual moving-throat branch:",
            "together with `mhat_0`.",
            "Compute **Packet B** from the actual branch:",
            "The symbolic stage is finished.",
        ],
    },
    {
        "id": "family1_support_non_bottleneck",
        "source": "notes/5pn/5pn_stage350_351_family1_branch_location.md",
        "role": "Shows the explicit support/source branch is not the active bottleneck.",
        "patterns": [
            "So the support/source side is no longer just",
            "on the canonical compact passive/outgoing branch one still has",
            "What is **still not numerically present in the files** is a completed PDE-selected point",
        ],
    },
    {
        "id": "barrier_actual_branch_finish_line",
        "source": "notes/summaries/barrier_audit_summary.md",
        "role": "Separates reduced corridor survival from final actual-branch verdict.",
        "patterns": [
            "support/source enhancement is already strong enough;",
            "the remaining unresolved kill test is the actual coherent placement / orbit packet together with outgoing normalization.",
            "So the final answer depends on the actual moving-throat branch data, not on another broad unresolved coefficient family.",
        ],
    },
    {
        "id": "strict_parent_wall_status",
        "source": "notes/moving_throat_pde_program_compact.md",
        "role": "Preserves the strict-parent versus effective-wall firewall.",
        "patterns": [
            "the wall PDE in Sections `3-4` is an effective closure.",
            "Wall PDE strict parent status: requires + S_Sigma[Sigma]",
            "Do not treat",
        ],
    },
    {
        "id": "no_refit_status_firewall",
        "source": "notes/pde_audit_full.md",
        "role": "Requires pre-target branch freeze and forbids post-hoc rescues.",
        "patterns": [
            "Define the parent action, gauge convention, wall/interface action, open-exit boundary protocol, projection/source map, support-profile family, number of modes/ports, stability gates, and extraction formulas **before** evaluating any target residual.",
            "A branch may be accepted only if the frozen pre-target packet produces stable coefficients satisfying the target residuals.",
            "A branch may not be rescued by changing mode support",
        ],
    },
    {
        "id": "uv_dressing_candidate",
        "source": "research/pde_ledger/notes/stages/moving_throat_pde_stage228_nonrigid_mouth_dressing_packet_and_uv_drain_compiler_sympy_audit.md",
        "role": "Promotable reduced U/V dressing lane for a future physical exporter.",
        "patterns": [
            "Reduced non-rigid mouth/dressing free energy",
            "The stationarity equations are linear:",
            "Then the exact stationary packet is",
            "So the dressing activation is a genuinely non-rigid effect rather than a rephrasing of the old transfer leg.",
        ],
    },
    {
        "id": "microscopic_export_kernel_candidate",
        "source": "research/pde_ledger/notes/stages/moving_throat_pde_stage234_microscopic_damping_export_kernel_replacing_the_phenomenological_v_leg_envelope_law_sympy_audit.md",
        "role": "Promotable odd passive/export kernel, but not an even N2/N4 shape controller by itself.",
        "patterns": [
            "Selected mixed quadrupole outlet",
            "First microscopic export kernel",
            "The first microscopic law is **super-Ohmic**",
            "the first microscopic active-leg equation is",
        ],
    },
    {
        "id": "constant_prefactor_target_equations",
        "source": "research/pde_ledger/notes/stages/moving_throat_pde_stage006_full_grouped_bundle.md",
        "role": "Defines the exact N2/N4 correlation target for outgoing moment shape.",
        "patterns": [
            "`P_2 = 0`,",
            "`P_4 = 0`.",
            "`N_2 = 2 D_2 N_0 / D_0`,",
            "outgoing-transfer moments `N_2,N_4` to be correlated",
        ],
    },
    {
        "id": "outgoing_support_blind_split",
        "source": "notes/5pn/5pn_stage205_normalized_packet_bridge.md",
        "role": "Shows support tuning alone cannot control the outgoing transfer moments.",
        "patterns": [
            "Exact support-blind theorem for the outgoing-transfer bundle",
            "is **exactly support-blind** in the explicit BdG pair",
            "the **outgoing transfer side** is support-blind",
            "the **conservative wall operator** is not.",
        ],
    },
    {
        "id": "leakage_work_support_side",
        "source": "research/pde_ledger/notes/stages/moving_throat_pde_stage227_selected_branch_leakage_and_scalar_photon_work_compiler_sympy_audit.md",
        "role": "Promotable leakage/work lane, explicitly support-side rather than orbit-side.",
        "patterns": [
            "proves that the whole leakage/work lane is **support-side**, not orbit-side.",
            "This proves that Stage 227 is a **selected-support compiler**, not an orbit-lock compiler.",
            "The leakage/work lane depends only on the support-side variables",
        ],
    },
    {
        "id": "atomic_p22_finite_throat_shape_source",
        "source": "notes/atomic_p22_bridge/paper/atomic_p22_bridge.tex",
        "role": "Concrete finite-throat P0/P2 mouth source and P22 shape-response physics.",
        "patterns": [
            "the first shape-sensitive atomic load is not the raw Coulomb depth itself",
            "the first shape-sensitive external loading on a defect is governed by the partner-field Hessian",
            "which is precisely the reduced atomic realization of the already inherited \\(P_0\\oplus P_2\\) support hierarchy.",
            "pure Coulomb loading excites the quadrupolar $P_2$ response channel",
            "The finite-throat kernel therefore matches the point-source tide in the far field but softens to zero at overlap",
        ],
    },
    {
        "id": "lepton_radiative_p0_flux_hook",
        "source": "notes/lepton_mass_notes.md",
        "role": "Concrete open/radiative scalar mouth-flux hook for future P0 normalization work.",
        "patterns": [
            "The exact mouth operator on the D/N branch is",
            "Exact quantum-geometric hammer stress",
            "a nonzero DC mouth source exists only on an explicitly open / radiative branch.",
            "\\dot M_{\\rm mouth}^{\\rm(dc,rad)}",
        ],
    },
    {
        "id": "lepton_scalar_hammer_p22_veto",
        "source": "notes/lepton_mass_notes.md",
        "role": "Rules out using scalar P0 hammering as direct P22/N2/N4 moment-shape control.",
        "patterns": [
            "The first is a symmetry veto: the exact scalar mouth hammer does **not** linearly drive the area-preserving \\(P_{22}\\) ellipse.",
            "The scalar \\(P_0\\) hammer does zero linear work on the \\(P_{22}\\) mouth mode",
            "\\boxed{P_0\\text{ hammer} \\not\\to P_{22}\\text{ oscillator at linear order.}}",
        ],
    },
    {
        "id": "atom_work_intrinsic_p22_bracing",
        "source": "notes/atom_work.md",
        "role": "Concrete intrinsic P22 bracing source conditional on half-flux/mixed-sector closure.",
        "patterns": [
            "Intrinsic core stress as a permanent quadrupole bracer",
            "So the atomic tide no longer has to create \\(P_{22}\\) from nothing.",
            "So the intrinsic \\(P_{22}\\) source is no longer a placeholder.",
            "\\boxed{h_\\alpha>0.}",
        ],
    },
    {
        "id": "simulation_coefficient_map",
        "source": "research/pde_audit/notes/stage_v2_19_isotropic_full_bundle_target_surface_derivation.md",
        "role": "Maps simulation A/C/D0 and outgoing prefactor residuals to reduced D/N moments.",
        "patterns": [
            "## 1. Isotropic full-bundle definitions",
            "The outgoing prefactor is defined by",
            "so the exact one-pole surface is",
            "The constant-prefactor branch requires",
            "On the one-pole surface `D_0C=3A^2`, this simplifies to",
        ],
    },
]


def sha256_json(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":"), default=float).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def find_line(text: str, pattern: str) -> int | None:
    index = text.find(pattern)
    if index < 0:
        return None
    return text.count("\n", 0, index) + 1


def check_anchor(anchor: Mapping[str, Any]) -> Dict[str, Any]:
    rel_source = str(anchor["source"])
    path = PROJECT_ROOT / rel_source
    if not path.exists():
        return {
            "id": anchor["id"],
            "source": rel_source,
            "role": anchor["role"],
            "exists": False,
            "pass": False,
            "issues": [f"missing source: {rel_source}"],
            "patterns": [],
        }

    text = path.read_text(encoding="utf-8")
    pattern_rows: List[Dict[str, Any]] = []
    issues: List[str] = []
    for pattern in anchor["patterns"]:
        line = find_line(text, str(pattern))
        found = line is not None
        if not found:
            issues.append(f"missing pattern in {rel_source}: {pattern!r}")
        pattern_rows.append({
            "pattern": pattern,
            "found": found,
            "line": line,
        })

    return {
        "id": anchor["id"],
        "source": rel_source,
        "role": anchor["role"],
        "exists": True,
        "source_sha256": sha256_text(text),
        "pass": not issues,
        "issues": issues,
        "patterns": pattern_rows,
    }


def build_report(anchor_checks: List[Dict[str, Any]]) -> Dict[str, Any]:
    failed = [item["id"] for item in anchor_checks if not item["pass"]]
    report: Dict[str, Any] = {
        "schema": "pde_audit_simulation_notes_intake_report/v1",
        "post_hoc_only": True,
        "candidate_generation_mutated": False,
        "target_residuals_used_to_generate_candidates": False,
        "target_blind_hit_claimed": False,
        "anchor_count": len(anchor_checks),
        "failed_anchor_count": len(failed),
        "failed_anchors": failed,
        "anchors": anchor_checks,
        "conclusion": {
            "does_more_derivation_help": True,
            "how_it_helps": (
                "The notes narrow the next simulation task to an actual moving-throat "
                "branch extraction with frozen parent/action, support, outgoing, and "
                "orbit-lock packets."
            ),
            "does_it_license_retuning_current_candidates": False,
            "retuning_reason": (
                "The notes repeatedly require pre-target freeze and warn against "
                "projecting realized branches back onto canonical outgoing or support "
                "surfaces after residuals are known."
            ),
            "support_source_bottleneck_active": False,
            "actual_branch_packet_required": True,
            "strict_parent_or_effective_closure_must_be_declared": True,
            "primary_next_artifact": "actual_branch_protocol_v1",
            "missing_physics_interpretation": (
                "The current exporter lacks the physical branch equations needed to "
                "select Packet A/Packet B and outgoing moment shape.  That is missing "
                "branch-realization physics in the simulation layer, not a proof that "
                "the symbolic target surface is inconsistent."
            ),
            "one_pole_support_alone_is_not_enough": True,
            "uniform_outgoing_scale_alone_is_not_enough": True,
            "outgoing_moment_shape_control_required": True,
            "atom_lepton_notes_are_useful": True,
            "atom_lepton_notes_are_ready_exporter": False,
            "atom_lepton_gap": (
                "The atom/lepton notes provide finite-throat P0/P2 forcing, "
                "scalar radiative P0 flux hooks, and intrinsic P22 bracing. "
                "They do not yet provide a calibrated actual-branch map into "
                "D0, C, P0, N2, and N4."
            ),
        },
        "actual_branch_finish_line": {
            "finish_line_conditions": [
                "dln_R_tr = 0",
                "dln_R_target = 0",
                "dln_epsilon_eta = 0",
                "N_Q = 1",
            ],
            "packet_A_fields": [
                "D_A0",
                "D_A2",
                "D_A4",
                "N_A0",
                "N_A2",
                "N_A4",
                "mhat_0",
            ],
            "packet_B_alternatives": [
                ["m_T", "m_K", "m_mu"],
                ["R_tr", "R_nt", "R_eta"],
                ["q_tr", "q_nt", "q_eta"],
            ],
            "constant_prefactor_shape_equations": [
                "N_2 = 2 D_2 N_0 / D_0",
                "N_4 = (2 D_0 (D_2 N_2 + D_4 N_0) - 3 D_2^2 N_0) / D_0^2",
            ],
            "simulation_aliases": {
                "A": "-D_2 = M + B2 + Z2",
                "C": "-D_4 = B4 + Z4",
                "D0": "K - B0 - Z0",
                "one_pole_surface": "D0*C/(3*A^2) = 1",
                "P0": "N0/D0",
                "P2": "(D0*N2 + 2*A*N0)/D0^2",
                "P4": "(D0^2*N4 + 2*D0*(A*N2 + C*N0) + 3*A^2*N0)/D0^3",
            },
        },
        "promotable_reduced_lanes": [
            {
                "name": "leakage_work",
                "source_anchor": "leakage_work_support_side",
                "recommended_use": "support-side source/work accounting only",
                "caveat": "not an orbit-lock compiler",
            },
            {
                "name": "uv_dressing",
                "source_anchor": "uv_dressing_candidate",
                "recommended_use": "non-rigid mouth/dressing branch component",
                "caveat": "must be frozen before target residuals and coupled to a real branch exporter",
            },
            {
                "name": "microscopic_export_kernel",
                "source_anchor": "microscopic_export_kernel_candidate",
                "recommended_use": "odd passive/export term in the active V equation",
                "caveat": "does not by itself supply independent even N2/N4 outgoing moment-shape control",
            },
            {
                "name": "atomic_finite_throat_p22",
                "source_anchor": "atomic_p22_finite_throat_shape_source",
                "recommended_use": "replace point-source forcing with finite-throat P0/P2 and P22 response sources",
                "caveat": "static shape/source physics, not a calibrated outgoing N2/N4 transfer law",
            },
            {
                "name": "radiative_scalar_p0_flux",
                "source_anchor": "lepton_radiative_p0_flux_hook",
                "recommended_use": "candidate scalar P0 normalization hook on an explicitly open/radiative branch",
                "caveat": "constitutive response coefficient and coupling to exporter P0 remain unresolved",
            },
            {
                "name": "intrinsic_p22_bracing",
                "source_anchor": "atom_work_intrinsic_p22_bracing",
                "recommended_use": "candidate intrinsic P22 mouth-shape source for future branch equations",
                "caveat": "conditional on half-flux/mixed-sector closure and still not an outgoing N2/N4 law",
            },
        ],
        "recommended_work_packages": [
            {
                "id": "WP0",
                "name": "freeze_protocol",
                "requirement": "Declare parent action or effective closure, gauge, boundary class, support basis, ports, stability gates, and extraction formulas before target evaluation.",
            },
            {
                "id": "WP1",
                "name": "stationary_isotropic_branch",
                "requirement": "Solve or continue the completed stationary branch without target feedback.",
            },
            {
                "id": "WP2",
                "name": "packet_A_and_outgoing_export",
                "requirement": "Extract D/N moments, mhat_0, N_Q or chi_Q, and source hashes from the realized branch.",
            },
            {
                "id": "WP3",
                "name": "weak_axisymmetric_tangent",
                "requirement": "Extract the first tangent needed for orbit-lock diagnostics on the same branch.",
            },
            {
                "id": "WP4",
                "name": "orbit_lock_packet",
                "requirement": "Evaluate dln_R_tr, dln_R_target, and dln_epsilon_eta without mixing in support enhancement.",
            },
            {
                "id": "WP5",
                "name": "frozen_verdict",
                "requirement": "Run the existing target-blind guard and post-hoc evaluator against the frozen actual-branch packet.",
            },
        ],
    }
    report["pass"] = bool(
        not failed
        and report["post_hoc_only"]
        and not report["candidate_generation_mutated"]
        and not report["target_residuals_used_to_generate_candidates"]
        and not report["target_blind_hit_claimed"]
    )
    report["report_hash"] = sha256_json(report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description="Build an evidence-backed intake of unincorporated derivation notes")
    parser.add_argument("--output-dir", default=str(OUTPUT_DIR), help="Simulation output directory")
    parser.add_argument("--report-prefix", default="notes_intake_report", help="Output basename without extension")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    anchor_checks = [check_anchor(anchor) for anchor in ANCHORS]
    report = build_report(anchor_checks)
    write_json(output_dir / f"{args.report_prefix}.json", report)

    conclusion = report["conclusion"]
    lines = [
        "PDE audit notes-intake report",
        "=" * 34,
        f"pass: {report['pass']}",
        f"anchor_count: {report['anchor_count']}",
        f"failed_anchor_count: {report['failed_anchor_count']}",
        f"does_more_derivation_help: {conclusion['does_more_derivation_help']}",
        f"actual_branch_packet_required: {conclusion['actual_branch_packet_required']}",
        f"support_source_bottleneck_active: {conclusion['support_source_bottleneck_active']}",
        f"retune_current_candidates_allowed: {conclusion['does_it_license_retuning_current_candidates']}",
        f"outgoing_moment_shape_control_required: {conclusion['outgoing_moment_shape_control_required']}",
        f"primary_next_artifact: {conclusion['primary_next_artifact']}",
        f"report_hash: {report['report_hash']}",
    ]
    for anchor in anchor_checks:
        lines.append(f"{'PASS' if anchor['pass'] else 'FAIL'}  {anchor['id']}")
    (output_dir / f"{args.report_prefix}.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0 if report["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
