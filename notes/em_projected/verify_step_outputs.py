#!/usr/bin/env python3
"""Run the SymPy audit scripts and check a small set of live output invariants.

This is a lightweight regression harness for the step_*_sympy.py bundle. It is
not a full markdown parser; it checks the parts of the live output that the
notes rely on most heavily.
"""
from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parent
DEFAULT_OUT_DIR = pathlib.Path("/tmp/em_projected_verify")
EXPECTATIONS: dict[str, list[str]] = {
    "step_04_projection_reduction_comparison_sympy.py": [
        "STATUS: PASS",
        "observer-dependent effective laws unless extra closure conditions are imposed.",
    ],
    "step_03_projected_maxwell_vector_sympy.py": [
        "STATUS: PASS",
        "Concrete checks: Faraday/Ampere sign maps, explicit boundary discharge, E-vs-D distinction, and nonzero leakage all pass.",
    ],
    "step_08_projected_maxwell_push_bundle_master_sympy.py": [
        "STATUS: PASS",
        "z0 cancellation from compatibility",
    ],
    "step_11_projected_maxwell_mouth_taylor_master_sympy.py": [
        "STATUS: PASS",
        "Checked one-sided Taylor projection, bottleneck dependencies, G_W transport entry, and mechanism sieve.",
    ],
    "step_12_projected_maxwell_mouth_taylor_gate_bridge_sympy.py": [
        "STATUS: PASS",
        "G_W' first tunes the constant-prefactor transport through delta P2 and delta P4.",
    ],
    "step_13_parent_throat_action_master_sympy.py": [
        "STATUS: PASS",
        "zero-support/zero-mixed specialization of the full even gates",
    ],
    "step_14_parent_throat_action_candidate_sympy.py": [
        "STATUS: PASS",
        "Y20 modal EL is recovered from the genuine S^2 reduction.",
    ],
    "step_15_parent_throat_action_weak_axisym_sympy.py": [
        "STATUS: PASS",
        "Therefore pure wall anisotropy closes the even gates only on the trivial branch delta_K = delta_M = 0.",
    ],
    "step_16_parent_throat_action_bundle_master_sympy.py": [
        "STATUS: PASS",
        "Checked isotropic compatibility, exact weak-axisymmetric wall-slope solve, and residual Xi1.",
    ],
    "step_17_parent_throat_action_isotropic_bundle_sympy.py": [
        "STATUS: PASS",
        "only the positive root gives u2>0; the negative root gives u2<0.",
    ],
    "step_19_parent_throat_action_actual_branch_export_sympy.py": [
        "STATUS: PASS",
        "branch_id = v2_local_parent_background_galerkin_demo",
        "Parent-background Galerkin port export on basis n=(0,2,4):",
        "Operator-adapted generalized eigenvalues:",
        "4-mode Z0, Z2, Z4 = 0.5083076399936368 -0.4143227047792866 0.3408544970748318",
        "5-mode mixed residual = sqrt(pi)*(-890747 + 629856*sqrt(2))/1259712",
        "adapted 5-mode Z0, Z2, Z4 = 0.5083084120357607 -0.4143247482933136 0.340857769373979",
        "residual packet labels = (R_pole, R_norm, R_P2, R_P4)",
        "R_pole = -13.13467188936131",
        "adapted 5-mode isotropic residues = (-13.134428427110006, -10.337190829817821, 0.37009832625683897, 0.8889265110475489)",
        "one-pole residue = -13.13467188936131",
        "static residue   = -10.33723418002391",
        "P2 residue       = 0.3700576386561055",
        "P4 residue       = 0.8888349036886998",
    ],
    "step_20_parent_throat_action_branch_family_scan_sympy.py": [
        "STATUS: PASS",
        "branch_id = v2_local_parent_background_reduced_family_scan",
        "residual packet (R_pole, R_norm, R_P2, R_P4) = (-13.134593938872369, -10.33719584868593, 0.3700984456976848, 0.8889149882257383)",
        "singular values = [71.044210206226, 1.376398462469, 0.14131229158]",
        "irreducible linearized norm = 10.108287522271972",
        "best actual residual norm = 15.828401363791917",
    ],
    "step_21_parent_throat_action_outgoing_family_scan_sympy.py": [
        "STATUS: PASS",
        "branch_id = v2_local_parent_background_outgoing_family_scan",
        "singular values = [71.044213096589, 1.715640711598, 0.141855721241, 0.108888020712]",
        "irreducible linearized norm = 4.4642745655749085e-13",
        "best actual residual norm = 11.328790887276268",
        "best actual residual packet = (3.421090838716161, -10.79989076029757, -6.610561955921492e-05, 4.642628358125823e-05)",
    ],
    "step_22_parent_throat_action_static_normalized_slice_sympy.py": [
        "STATUS: PASS",
        "branch_id = v2_local_parent_background_static_normalized_slice",
        "root scale = 0.091265380859375",
        "mhat0_static = 206.8317356734255",
        "packet norm = 0.0002055057029418912",
        "not a realized-branch success verdict",
    ],
    "step_23_parent_throat_action_normalization_frontier_sympy.py": [
        "STATUS: PASS",
        "branch_id = v2_local_parent_background_normalization_frontier_scan",
        "sample_scales = [0.0, 0.03, 0.05, 0.08, 0.088, 0.09, 0.092]",
        "first sampled point with Q_iso < 1 = (0.09, 194.6081703105869, 0.4513177752288337)",
        "best sampled Q_iso point = (0.092, 214.2709506975902, 0.26705543084121786)",
        "Q_iso first drops below 1 only once mhat0_req is already about 1.95e2.",
    ],
    "step_24_parent_throat_action_outgoing_amplitude_frontier_sympy.py": [
        "STATUS: PASS",
        "branch_id = v2_local_parent_background_outgoing_amplitude_frontier_scan",
        "sample_lambda_out = [1, 5, 20, 50, 100, 200, 500, 1000, 2000]",
        "mhat0_req <= 5   -> (0.092, 2000.0, 4.7912441136331765, 0.4394839373049669)",
        "Q_iso <= 1.0 -> (0.09, 2000.0, 4.351570977913287, 0.618690285150578)",
        "Q_iso <= 0.5 -> (0.092, 2000.0, 4.7912441136331765, 0.4394839373049669)",
        "static normalization no longer looks structurally fatal if this outgoing-amplitude family is physically admissible.",
    ],
    "step_25_parent_throat_action_outgoing_amplitude_admissibility_sympy.py": [
        "STATUS: PASS",
        "branch_id = v2_local_parent_background_outgoing_amplitude_admissibility",
        "X = mhat0^2 N_Q = P0_target/P0_base",
        "chi_Q = P0_base/P0_target",
        "N_Q   = P0_target/(P0_base*mhat0**2)",
        "elimination Jacobian determinant = P0_target",
        "inferred chi_Q = 2.6404494712422554e-05",
        "outgoing finish-line defect N_Q - 1 = 1999.0",
        "lambda_out mutation residual (1999 instead of 2000) = -0.000500000000000056",
        "step-24 improvement does not rescue outgoing admissibility.",
    ],
    "step_26_parent_throat_action_natural_source_map_outgoing_burden_sympy.py": [
        "STATUS: PASS",
        "branch_id = v2_local_parent_background_natural_source_burden",
        "natural-source N_Q burden = lambda_out*mhat0_req**2",
        "minimal Robin rho(N_Q) = 3 - 3*N_Q_nat",
        "natural-source N_Q burden = 37872.3399516344",
        "natural-source N_Q burden = 45912.04031284913",
        "minimal Robin rho = -113614.01985490319",
        "minimal Robin rho = -137733.1209385474",
        "At fixed frozen scale it only redistributes that burden between lambda_out and mhat0_req.",
    ],
    "step_27_parent_throat_action_transfer_amplitude_interpretation_sympy.py": [
        "STATUS: PASS",
        "branch_id = v2_local_parent_background_transfer_amplitude_interpretation",
        "exact transfer-shape amplitude law = T_eff^2 -> lambda_out * T_eff^2",
        "canonical branch-B chi_Q = 1",
        "sigma(lambda_out) = 1 - lambda_out",
        "hybrid branch-B sigma = -19.0",
        "hybrid branch-B sigma = -1999.0",
        "Stage95 branch B gives a concrete exact outlet family where chi_Q stays canonical while the overall outgoing amplitude scales by 1 - sigma.",
    ],
    "step_28_parent_throat_action_hybrid_amplitude_budget_frontier_sympy.py": [
        "STATUS: PASS",
        "branch_id = v2_local_parent_background_hybrid_amplitude_budget_frontier",
        "sigma_budgets = [4, 19, 49, 199, 1999]",
        "|sigma| <= 4 -> (0.092, 1.0, 0.0, 214.2709506975902, 0.26705543084121786)",
        "|sigma| <= 19 -> (0.092, 1.0, 0.0, 214.2709506975902, 0.26705543084121786)",
        "|sigma| <= 4 -> (0.09, 5.0, -4.0, 87.03141955826572, 0.4513189656743675)",
        "|sigma| <= 19 -> (0.09, 20.0, -19.0, 43.51570977913286, 0.45133756597801233)",
        "|sigma| <= 1999 -> (0.09, 2000.0, -1999.0, 4.351570977913287, 0.618690285150578)",
        "The best-Q frontier measures how much isotropic defect reduction is available at each amplitude budget.",
    ],
    "step_29_parent_throat_action_moderate_branchb_sector_sympy.py": [
        "STATUS: PASS",
        "branch_id = v2_local_parent_background_moderate_branchb_sector",
        "sigma_budgets = [19, 49]",
        "q_sector_cap = 1.5",
        "(0.092, 20.0, -19.0, 47.912441136331765, 0.2670781822626949, 2.995732273553991)",
        "(0.092, 50.0, -49.0, 30.302488449910886, 0.26719789465663646, 3.912023005428146)",
        "q increase = 2.2751421477",
        "normalization drop factor = 4.47213595499958",
        "q increase = 0.0001424638154186",
        "normalization drop factor = 7.071067811865475",
        "it is the moderate sigma = -19 to -49 corridor.",
    ],
    "step_30_parent_throat_action_branchb_parameter_burden_sympy.py": [
        "STATUS: PASS",
        "branch_id = v2_local_parent_background_branchb_parameter_burden",
        "rho = 4*sigma",
        "kappa = 1/3",
        "canonical gamma = 1/9",
        "branch-B rho = -76.0",
        "branch-B rho = -196.0",
        "natural-source Robin burden ratio = 1494.9213138803052",
        "natural-source Robin burden ratio = 702.7200047885071",
        "it is quantitatively much cheaper at the outlet-parameter level.",
    ],
    "step_31_parent_throat_action_electron_fast_falsification_sympy.py": [
        "STATUS: PASS",
        "branch_id = v2_local_parent_background_electron_fast_falsification",
        "P0_base = 0.00023523241237827255",
        "lambda_out = 50 -> P0 = 0.011761620618913627 mhat0 = 30.302488449910886 mhat0^2 P0 = 10.8",
        "S_port * (mhat0^2 lambda_out P0_base) = 54 G c_s^5 / (5 r_e^5 c^5)",
        "c_s = 9.159491211330665e-05",
        "Omega_Q = 48756308603.324455",
        "Omega_Q / omega_Compton = 6.280244857660478e-11",
        "required K = 1.1054545016851398e+52",
        "required S_port = 1.0235689830417961e+51",
        "lambda load guard = 50.00000000000001",
        "the naive direct electron identification is falsified quickly",
    ],
    "step_32_parent_throat_action_dimensional_port_map_sympy.py": [
        "STATUS: PASS",
        "branch_id = v2_local_parent_background_dimensional_port_map",
        "n = 5",
        "alpha_opt = 2",
        "kappa_add = 1/2",
        "alpha^2 = 3/4",
        "K_vec = 2/pi**2",
        "S_port(P0_red,mhat0,c_s,a) = 54*G*c_s**5/(5*P0_red*a**5*c**5*mhat_0**2)",
        "On the patched target sheet mhat0^2 P0_red = 54/5:",
        "S_port = G*c_s**5/(a**5*c**5)",
        "S_port = 1  =>  c_s = a*c/G**(1/5)",
        "lambda_out = 1 - sigma",
        "gamma_quad_eff = 2*G/(5*c**5)",
    ],
    "step_33_parent_throat_action_runtime_monitor_falsifier_sympy.py": [
        "STATUS: PASS",
        "branch_id = v2_local_parent_background_runtime_monitor_falsifier",
        "S_rho     = -[W J_w] + ∫ W' J_w dw",
        "R_cont = -S_rho + divJ + dtrho",
        "R_Pois_lin   = -S_rho + lap_phi*rho_0",
        "Newton plateau  Q_r = 4*pi*A",
        "Yukawa flux     Q_r = 4*pi*A*(mu*r + 1)*exp(-mu*r)",
        "1/r^2 impostor Q_r = 8*pi*A/r",
        "mu_eff^2(Newton) = 0",
        "mu_eff^2(Yukawa) = mu**2",
        "mu_eff^2(1/r^2 impostor) = 2/r**2",
        "alpha_n = 2",
        "R_opt = (-2*Phi_eff + c_probe**2*(1 - N_probe))/Phi_eff",
        "Concrete-flow monitor exhibit:",
        "periodic max |R_Pois_exact| = 0.0",
        "Bad-optics alpha_fit tail = 1.4",
        "The decisive hard falsifiers are Yukawa-like screening in the exterior field and the wrong weak-field optical coefficient.",
    ],
    "step_34_parent_throat_action_cfd_runtime_postprocessor_sympy.py": [
        "STATUS: PASS",
        "branch_id = v2_local_parent_background_cfd_runtime_postprocessor",
        "max |R_cont| = 1.4210854715202004e-14",
        "max |R_Pois_exact| = 0.0",
        "alpha_fit tail = 2.0000000000000004",
        "Newton Q_r tail cv = 0.00048311949183278685",
        "Newton mu_eff^2 tail median = -0.003597967054167384",
        "Yukawa Q_r tail cv = 0.25427545643672567",
        "Yukawa mu_eff^2 tail median = 1.954292666547206",
        "monitor CLI = python cfd_runtime_monitor_postprocess.py snapshot.npz --output-json summary.json",
        "The postprocessor now computes S_rho, brane continuity residuals, exact/linearized Poisson residuals, exterior Q_r plateaus, mu_eff^2 tails, and alpha_fit tails directly from a CFD snapshot.",
    ],
    "step_35_parent_throat_action_failfast_classifier_sympy.py": [
        "STATUS: PASS",
        "branch_id = v2_local_parent_background_failfast_classifier",
        "continuity_rel_max = 0.05",
        "poisson_rel_max = 0.05",
        "q_tail_cv_max = 0.05",
        "mu_eff2_tail_abs_max = 0.25",
        "alpha_fit_tail_error_max = 0.1",
        "alpha_fit_tail_std_max = 0.1",
        "Newton-like exterior = PASS",
        "Yukawa exterior = FAIL",
        "Bad optics exterior = FAIL",
        "Projection-broken snapshot = FAIL",
        "Missing optics snapshot = INCOMPLETE",
        "Near-zero source snapshot = INCOMPLETE",
        "Yukawa: exterior Q_r plateau is too noisy: 0.254275; effective exterior mass scale is too large: 1.95429",
        "Bad optics: alpha_fit mean is not close to 2: 1.4",
        "Projection break: continuity residual too large: 0.3; exact Poisson residual too large: 0.3; exterior Q_r plateau is too noisy: 2.33516",
        "Near-zero source: source scale is near zero; projection residuals are not load-bearing",
        "classify    = python cfd_runtime_failfast.py summary.json --output-json verdict.json",
    ],
    "step_36_parent_throat_action_snapshot_adapter_sympy.py": [
        "STATUS: PASS",
        "branch_id = v2_local_parent_background_snapshot_adapter",
        "4D wavefunction snapshot -> full_4d monitor schema",
        "3D monopole state dump   -> projected_3d monitor schema",
        "source_schema = full_4d",
        "rel jx error = 5.470539562816678e-05",
        "rel jy error = 5.470539562816674e-05",
        "rel jz error = 0.00020415416182304036",
        "max |R_cont| = 1.8245167713208398e-16",
        "source_schema = projected_3d",
        "max |R_cont| = 1.4210854715202004e-14",
        "max |R_Pois_exact| = 0.0",
        "reconstructed source rel error = 0.0",
        "missing W without --allow-uniform-W = rejected",
        "reconstructed monopole source without lambda = rejected",
        "wavefunction -> python cfd_snapshot_adapters.py wavefunction-4d raw_wave.npz runtime_snapshot.npz",
        "wavefunction fallback W -> add --allow-uniform-W explicitly",
        "monopole     -> python cfd_snapshot_adapters.py monopole-3d raw_state.npz runtime_snapshot.npz",
    ],
    "step_37_parent_throat_action_monopole_jsonl_fastscreen_sympy.py": [
        "STATUS: PASS",
        "branch_id = v2_local_parent_background_monopole_jsonl_fastscreen",
        "dP_slope target = -1 +/- 0.35",
        "geff_slope target = -2 +/- 0.35",
        "mach_max <= 0.6",
        "min_fit_points = 8",
        "good log = PASS",
        "bad log = FAIL",
        "weak log = INCOMPLETE",
        "threshold-boundary log = PASS",
        "just-outside-threshold log = FAIL",
        "malformed JSONL log = INCOMPLETE",
        "bad log: dP slope misses -1 target: -0.32; g_eff slope misses -2 target: -1.11; flow is too compressible for weak-field screen: mach_max=0.84",
        "weak log warnings: insufficient dP fit points: 4; insufficient geff fit points: 5",
        "malformed log warnings: malformed JSONL line 2",
        "python single_throat_monopole_jsonl_fastscreen.py monopole.log --output-json monopole_verdict.json",
    ],
}


def run_script(path: pathlib.Path, out_dir: pathlib.Path) -> str:
    proc = subprocess.run(
        [sys.executable, str(path)],
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    output = proc.stdout
    (out_dir / f"{path.stem}.out").write_text(output)
    if proc.returncode != 0:
        raise RuntimeError(f"{path.name} exited with code {proc.returncode}\n{output}")
    if "STATUS: PASS" not in output:
        raise RuntimeError(f"{path.name} did not report STATUS: PASS\n{output}")
    return output


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("scripts", nargs="*")
    args = parser.parse_args()

    out_dir = pathlib.Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.scripts:
        targets = [ROOT / script for script in args.scripts]
    else:
        targets = sorted(ROOT.glob("step_*_sympy.py"))

    for path in targets:
        output = run_script(path, out_dir)
        for needle in EXPECTATIONS.get(path.name, []):
            if needle not in output:
                raise RuntimeError(f"{path.name} missing expected output fragment:\n{needle}")
        print(f"OK {path.name}")

    print(f"Wrote logs to {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
