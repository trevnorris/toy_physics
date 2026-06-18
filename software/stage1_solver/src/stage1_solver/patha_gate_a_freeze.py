"""Path-A GATE-A freeze stamp for decision-11.

This module only writes the target-blind freeze sheet. It intentionally does
not run a closed solve, extraction, calibration, or observable evaluation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
from typing import Any

from .backend import library_versions
from .config import source_revision


FROZEN_ROOT = Path("software/stage1_solver/frozen/pathA_gate_a")
SCHEMA = "pathA_gate_a_decision_11_freeze/v1"


def canonical_json(data: Any) -> str:
    return json.dumps(data, sort_keys=True, separators=(",", ":"))


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def repo_root() -> Path:
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
        if out:
            return Path(out)
    except (OSError, subprocess.CalledProcessError):
        pass
    return Path(__file__).resolve().parents[4]


def _source_hashes(root: Path) -> dict[str, str]:
    rel_paths = [
        "software/stage1_solver/src/stage1_solver/patha_static_balance.py",
        "software/stage1_solver/src/stage1_solver/patha_gate_a_freeze.py",
        "software/stage1_solver/src/stage1_solver/patha_closed_newton.py",
        "software/stage1_solver/src/stage1_solver/patha_closed_validation.py",
        "software/stage1_solver/src/stage1_solver/coupled_branch.py",
        "software/stage1_solver/src/stage1_solver/config.py",
        "research/pde_audit/scripts/stage_v2_21_branch_extraction_fixture.py",
        "research/pde_audit/scripts/stage_v2_22a_profile_to_coefficient_adapter.py",
    ]
    hashes: dict[str, str] = {}
    for rel in rel_paths:
        path = root / rel
        hashes[rel] = sha256_file(path) if path.exists() else "missing"
    return hashes


def build_freeze_sheet() -> dict[str, Any]:
    root = repo_root()
    return {
        "schema": SCHEMA,
        "parent_action_status": "path_A_promoted_S_Sigma",
        "branch_identity": {
            "decision_record": "software/stage1_solver/decisions/11_pathA_gate_a_freeze_sheet.md",
            "gate": "GATE-A",
            "target_blind": True,
            "no_post_residual_refit": True,
            "firewalled_simulation_tree": "not_touched",
            "physical_export_permitted_guard": "not_touched",
        },
        "family": {
            "id": "homogeneous_isotropic_hooke_v1",
            "parameters": ["tau", "a", "w_min", "w_max"],
            "calibrated_parameters": ["tau"],
            "pinned_parameters": {"a": "geometry.a", "w_min": "0", "w_max": "L"},
            "parameter_domains": {"tau": "tau>0", "a": "a>0", "w": "w_min<=w<=w_max"},
            "forms": [
                "mu(R,w)=tau",
                "T_w(R,w)=tau",
                "T_Omega(R,w)=tau/a^2",
                "U(R,w)=1/2*(tau/a^2)*(R-a)^2",
            ],
            "constitutive_derivatives": [
                "U_R(R,w)=(tau/a^2)*(R-a)",
                "U_RR(R,w)=tau/a^2",
                "T_w_R(R,w)=0",
                "T_w_RR(R,w)=0",
                "mu_R(R,w)=0",
                "T_Omega_R(R,w)=0",
            ],
            "ties": ["mu_Sigma=tau", "T_w=tau", "T_Omega=tau/a^2"],
            "g": "0",
            "domains": {"w": "[0,L]", "R": "R>0"},
            "units": {"mu_Sigma": "1", "T_w": "1", "T_Omega": "L^-2", "U": "1"},
            "R_star": "a",
            "admissibility": [
                "tau>0",
                "a>0",
                "U is bounded below with a single well at R=a",
            ],
        },
        "geometry": {
            "a": 1,
            "L": 1.85,
            "L_exact": "37/20",
            "boundary_class": {
                "mouth": "Dirichlet R0(0)=a",
                "exit": "natural zero-traction T_w*R0_prime(L)=0",
                "Y_L": "0",
                "hard_cap": "none",
                "positivity": "R0(w)>0",
            },
            "modal_boundary_conditions": [
                "eta(0)=0",
                "T_w*eta_prime(L)=0",
            ],
        },
        "source_port": {
            "mhat0": 1,
            "S_port": 1,
            "gauge": "canonical real STF grouped P2 basis; H=Z convention",
            "constants": {"G": 1, "c": 1, "c_s": 1},
        },
        "calibration_objective": {
            "anchor_equivalence": "R_norm=0 <=> P0=N0/D0=54/5 under mhat0^2*S_port=1",
            "anchor_constant": "54*G*c_s^5/(5*a^5*c^5)",
            "root_finder": "deterministic one-dimensional root-find on tau",
            "branch_selection": "stable-side D0>0",
            "naturalness_reporting": [
                "report tau_star",
                "report |ln tau_star|",
                "report K/(B0+Z0) cancellation ratio and digit count",
                "no hard exclusion from naturalness band",
            ],
            "mode_selection": {
                "channel": "l=2",
                "branch": "lowest-positive L2 eigenmode",
                "normalization": "integral chi^2 dw=1",
                "orientation": "integral chi dw>0",
            },
            "tolerances": {
                "closed_newton_linf_atol": "2.0e-9",
                "closed_newton_linf_rtol": "2.0e-9",
                "step_atol": "1.0e-12",
                "step_rtol": "1.0e-11",
                "gmres_rtol": "2.0e-9",
                "gmres_atol": "1.0e-11",
                "gmres_restart": "160",
                "gmres_maxiter": "12",
                "line_search_max_iters": "18",
                "finite_difference_jvp_epsilon": "1.0e-5",
            },
            "mesh_grid_convergence_ladder": {
                "chunk_1b_smoke_grid": [6, 6],
                "chunk_1c_validation_levels": [[4, 4], [8, 8], [16, 16]],
                "refinement_ratio": 2,
                "min_levels": 3,
                "expected_order": 2.0,
                "max_level_seconds": 590.0,
            },
        },
        "extraction_formulas": [
            "D0=K-B0-Z0",
            "D2=-(M+B2+Z2)",
            "D4=-(B4+Z4)",
            "P0=N0/D0",
            "P2=(D0*N2-2*D2*N0)/D0^2",
            "P4=(D0^2*N4-2*D0*(D2*N2+D4*N0)+3*D2^2*N0)/D0^3",
            "R_norm=mhat0^2*S_port*N0/D0-54*G*c_s^5/(5*a^5*c^5)",
            "R_pole=D0*(B4+Z4)-3*(M+B2+Z2)^2",
            "K=tau*kappahat",
        ],
        "source_revision": {
            "git_head": source_revision(),
            "library_versions": library_versions(),
            "source_file_sha256": _source_hashes(root),
        },
    }


def write_freeze() -> tuple[str, Path]:
    root = repo_root()
    sheet = build_freeze_sheet()
    text = canonical_json(sheet)
    freeze_hash = sha256_text(text)
    out_dir = root / FROZEN_ROOT / freeze_hash
    out_dir.mkdir(parents=True, exist_ok=True)
    freeze_path = out_dir / "freeze_sheet.json"
    freeze_path.write_text(text, encoding="utf-8")
    (out_dir / "freeze_hash.txt").write_text(freeze_hash + "\n", encoding="utf-8")
    return freeze_hash, freeze_path


def load_freeze(freeze_hash: str) -> tuple[dict[str, Any], Path]:
    root = repo_root()
    freeze_path = root / FROZEN_ROOT / freeze_hash / "freeze_sheet.json"
    if not freeze_path.exists():
        raise FileNotFoundError(f"freeze sheet not found: {freeze_path}")
    text = freeze_path.read_text(encoding="utf-8")
    actual = sha256_text(text)
    if actual != freeze_hash:
        raise ValueError(f"freeze hash mismatch: dir={freeze_hash} content={actual}")
    return json.loads(text), freeze_path


def scan_freeze_sheet(sheet: dict[str, Any]) -> dict[str, Any]:
    """Check for forbidden computed-output fields while allowing spec strings."""

    forbidden_numeric_symbols = {"R_norm", "D0", "P0", "N0"}
    forbidden_field_terms = ("pass", "fail", "target_value")
    allowed_field_names = {"target_blind", "no_post_residual_refit"}
    issues: list[str] = []

    def walk(value: Any, path: str) -> None:
        if isinstance(value, dict):
            for key, child in value.items():
                key_text = str(key)
                key_lower = key_text.lower()
                child_path = f"{path}.{key_text}" if path else key_text
                if key_text in forbidden_numeric_symbols and isinstance(child, (int, float)):
                    issues.append(f"{child_path} stores a numeric computed coefficient")
                if key_text not in allowed_field_names:
                    for term in forbidden_field_terms:
                        if term in key_lower:
                            issues.append(f"{child_path} uses forbidden field term {term!r}")
                    if "residual" in key_lower:
                        issues.append(f"{child_path} uses forbidden field term 'residual'")
                walk(child, child_path)
        elif isinstance(value, list):
            for index, child in enumerate(value):
                walk(child, f"{path}[{index}]")

    walk(sheet, "")
    return {"passed": not issues, "issues": issues}


def main() -> int:
    parser = argparse.ArgumentParser(description="Path-A GATE-A freeze stamp helper")
    sub = parser.add_subparsers(dest="cmd", required=True)
    sub.add_parser("freeze")
    scan = sub.add_parser("scan")
    scan.add_argument("freeze_hash")
    args = parser.parse_args()

    if args.cmd == "freeze":
        freeze_hash, freeze_path = write_freeze()
        print(f"candidate_freeze_hash: {freeze_hash}")
        print(f"freeze_sheet: {freeze_path}")
        return 0
    if args.cmd == "scan":
        sheet, freeze_path = load_freeze(args.freeze_hash)
        result = scan_freeze_sheet(sheet)
        print(json.dumps({"freeze_sheet": str(freeze_path), **result}, indent=2, sort_keys=True))
        return 0 if result["passed"] else 1
    raise AssertionError(args.cmd)


if __name__ == "__main__":
    raise SystemExit(main())
