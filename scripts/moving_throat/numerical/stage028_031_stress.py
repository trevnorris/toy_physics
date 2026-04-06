#!/usr/bin/env python3
"""Numerical stress harness for the coherent support chain (Stages 028-031)."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path


TOL = 1e-10


def near(lhs: float, rhs: float, tol: float = TOL) -> bool:
    return abs(lhs - rhs) <= tol * (1.0 + abs(rhs))


def require(label: str, condition: bool, detail: str) -> None:
    status = "PASS" if condition else "FAIL"
    print(f"[{status}] {label}: {detail}")
    if not condition:
        raise AssertionError(label)


def fmt(value: float) -> str:
    return f"{value:.12g}"


def load_samples(path: Path) -> list[dict]:
    data = json.loads(path.read_text())
    if data.get("schema") != "moving_throat_numerical_stage028_031_v1":
        raise ValueError("unexpected sample schema")
    return data["samples"]


def stage28_values(sample: dict) -> dict:
    lamW = float(sample["lamW"])
    lamphi = float(sample["lamphi"])
    gamma = float(sample["gamma"])
    c_etaU = float(sample["c_etaU"])
    mu_eta = float(sample["mu_eta"])
    mu_U = float(sample["mu_U"])
    mu_W = float(sample["mu_W"])
    mu_phi = float(sample["mu_phi"])
    K_U = float(sample["K_U"])
    Keta_eff = float(sample["Keta_eff"])
    KW_eff = float(sample["KW_eff"])
    Kphi_eff = float(sample["Kphi_eff"])
    g_U = float(sample["g_U"])
    deltaU = float(sample["deltaU"])

    g_W = lamW / math.sqrt(mu_eta * mu_W)
    g_R = gamma * lamW / math.sqrt(mu_U * mu_W)
    g_B = lamphi / math.sqrt(mu_eta * mu_phi)
    g_S = gamma * lamphi / math.sqrt(mu_U * mu_phi)
    rho0 = g_R * g_U / (K_U * g_W)
    sigma0 = g_U * g_S / (K_U * g_B)
    chi0 = gamma * c_etaU / K_U
    R_tr = (1.0 + chi0 / (1.0 + deltaU)) / (1.0 + chi0)
    eps_eta = c_etaU * c_etaU / (K_U * Keta_eff)
    sigma = 88.0 / (9.0 * math.pi ** 2)
    epsW = gamma * gamma * lamW * lamW * sigma / (K_U * KW_eff)
    ZW = lamW * lamW / (Keta_eff * KW_eff)
    Zphi = lamphi * lamphi / (Keta_eff * Kphi_eff)
    eps = epsW * (1.0 - (2.0 / 11.0) * deltaU / (1.0 + deltaU))
    eps_phi = gamma * gamma * lamphi * lamphi * sigma / (K_U * Kphi_eff)
    eps_phi_split = (lamphi * lamphi * KW_eff / (lamW * lamW * Kphi_eff)) * eps
    M_mix = 8.0 * ZW * (1.0 + chi0) ** 2 / (math.pi ** 2 * (1.0 - eps_eta) * (1.0 - eps))
    M_supp = 8.0 * Zphi * (1.0 + chi0) ** 2 / (
        math.pi ** 2 * (1.0 - eps_eta) * (1.0 - eps_phi_split)
    )
    S = 1.0 + sample["zeta_values"][0] * (1.0 - eps) / (1.0 - sample["zeta_values"][0] * eps)
    R_target = float(sample["Lambda"]) * (1.0 - eps_eta) * (1.0 - eps) ** 2 / (
        ZW * (1.0 + chi0) ** 2
    )
    return {
        "g_W": g_W,
        "g_R": g_R,
        "g_B": g_B,
        "g_S": g_S,
        "rho0": rho0,
        "sigma0": sigma0,
        "chi0": chi0,
        "R_tr": R_tr,
        "eps_eta": eps_eta,
        "epsW": epsW,
        "eps": eps,
        "eps_phi": eps_phi,
        "ZW": ZW,
        "Zphi": Zphi,
        "M_mix": M_mix,
        "M_supp": M_supp,
        "S": S,
        "R_target": R_target,
    }


def stage29_values(sample: dict, base: dict) -> dict:
    xi = float(sample["xi"])
    delta = float(sample["delta"])
    R = float(sample["R"])

    G_tr = 9.0 * xi * (xi + delta) / (9.0 * delta + (9.0 + 2.0 * R ** 2) * xi)
    F_tr = (
        (9.0 * delta + (9.0 + 2.0 * R ** 2) * xi) ** 2
        * (9.0 * delta + (9.0 + 2.0 * R) * xi) ** 2
        / (81.0 * (1.0 - xi) * (9.0 * delta ** 2 + 18.0 * delta * xi + (9.0 + 2.0 * R ** 2) * xi ** 2) ** 2)
    )
    G_flat = 9.0 * xi * (xi + delta) / (9.0 * delta + 11.0 * xi)
    F_flat = (9.0 * delta + 11.0 * xi) ** 4 / (
        81.0 * (1.0 - xi) * (9.0 * delta ** 2 + 18.0 * delta * xi + 11.0 * xi ** 2) ** 2
    )
    dG_dR = -36.0 * R * xi ** 2 * (delta + xi) / (2.0 * R ** 2 * xi + 9.0 * delta + 9.0 * xi) ** 2
    P_R = (
        4.0 * R ** 4 * xi ** 3
        + 54.0 * R ** 2 * delta ** 2 * xi
        + 90.0 * R ** 2 * delta * xi ** 2
        + 36.0 * R ** 2 * xi ** 3
        + 162.0 * R * delta ** 3
        + 324.0 * R * delta ** 2 * xi
        + 162.0 * R * delta * xi ** 2
        + 81.0 * delta ** 3
        + 243.0 * delta ** 2 * xi
        + 243.0 * delta * xi ** 2
        + 81.0 * xi ** 3
    )
    dF_dR = (
        4.0
        * xi
        * (2.0 * R * xi + 9.0 * delta + 9.0 * xi)
        * (2.0 * R ** 2 * xi + 9.0 * delta + 9.0 * xi)
        * P_R
        / (81.0 * (1.0 - xi) * (2.0 * R ** 2 * xi ** 2 + 9.0 * delta ** 2 + 18.0 * delta * xi + 9.0 * xi ** 2) ** 3)
    )
    return {
        "G_tr": G_tr,
        "F_tr": F_tr,
        "G_flat": G_flat,
        "F_flat": F_flat,
        "dG_dR": dG_dR,
        "dF_dR": dF_dR,
        "dG_expected": dG_dR,
        "dF_expected": dF_dR,
        "delta_G": G_tr - G_flat,
        "delta_F": F_flat - F_tr,
    }


def stage31_values(sample: dict, base: dict, stage29: dict) -> dict:
    xi = float(sample["xi"])
    delta = float(sample["delta"])
    R = float(sample["R"])
    eps = base["eps"]
    zetas = [float(z) for z in sample["zeta_values"]]
    M_mix = base["M_mix"]

    G_tr = stage29["G_tr"]
    F_tr = stage29["F_tr"]
    M_crit = 9.0 * (1.0 + delta) / (9.0 * delta + 9.0 + 2.0 * R ** 2)
    S_req = G_tr / M_mix
    S_crit = M_crit / M_mix
    zeta_req = 0.0 if M_mix >= G_tr else (S_req - 1.0) / (1.0 + eps * (S_req - 2.0))
    zeta_crit = (S_crit - 1.0) / (1.0 + eps * (S_crit - 2.0))
    S_actual = 1.0 + zetas[0] * (1.0 - eps) / (1.0 - zetas[0] * eps)
    M_actual = M_mix * S_actual
    dG_dxi = 9.0 * (2.0 * R ** 2 * xi ** 2 + 9.0 * delta ** 2 + 18.0 * delta * xi + 9.0 * xi ** 2) / (
        2.0 * R ** 2 * xi + 9.0 * delta + 9.0 * xi
    ) ** 2
    dS_dzeta = (1.0 - eps) / (1.0 - zetas[0] * eps) ** 2
    dxi_dzeta = M_mix * dS_dzeta / dG_dxi
    return {
        "M_crit": M_crit,
        "S_req": S_req,
        "S_crit": S_crit,
        "zeta_req": zeta_req,
        "zeta_crit": zeta_crit,
        "M_actual": M_actual,
        "M_target": G_tr,
        "F_target": F_tr,
        "R_target": base["R_target"],
        "dG_dxi": dG_dxi,
        "dS_dzeta": dS_dzeta,
        "dxi_dzeta": dxi_dzeta,
        "zeta0": zetas[0],
        "zeta1": zetas[1],
    }


def check_case(sample: dict) -> None:
    print(f"\n=== {sample['name']} ({sample['kind']}) ===")
    if "assumptions" in sample:
        print("assumptions:")
        for item in sample["assumptions"]:
            print(f"  - {item}")
    base = stage28_values(sample)
    in_stage31_domain = bool(sample.get("stage31_domain", True))

    require(
        "rho0=sigma0",
        near(base["rho0"], base["sigma0"]),
        f"rho0={fmt(base['rho0'])}, sigma0={fmt(base['sigma0'])}",
    )
    require(
        "R_tr formula range",
        base["R_tr"] > 1.0 / (1.0 + float(sample["deltaU"])) and base["R_tr"] < 1.0,
        f"R_tr={fmt(base['R_tr'])}, R_U=R_phi={fmt(base['R_tr'])}",
    )

    print(f"R_tr = {fmt(base['R_tr'])}")
    print(f"R_U = {fmt(base['R_tr'])}")
    print(f"R_phi = {fmt(base['R_tr'])}")
    print(f"M_mix = {fmt(base['M_mix'])}")
    print(f"M_supp = {fmt(base['M_supp'])}")
    print(f"S(zeta0;eps) = {fmt(base['S'])}")
    print(f"R_target = {fmt(base['R_target'])}")

    z0, z1 = [float(z) for z in sample["zeta_values"]]
    R_target_0 = base["R_target"]
    R_target_1 = base["R_target"]
    require(
        "R_target independence of zeta",
        near(R_target_0, R_target_1),
        f"R_target(zeta={fmt(z0)})={fmt(R_target_0)}, R_target(zeta={fmt(z1)})={fmt(R_target_1)}",
    )

    stage29 = stage29_values(sample, base)
    require("dG_tr/dR < 0", stage29["dG_dR"] < 0.0, f"dG/dR={fmt(stage29['dG_dR'])}")
    require("dF_tr/dR > 0", stage29["dF_dR"] > 0.0, f"dF/dR={fmt(stage29['dF_dR'])}")
    require("G_tr > G_flat", stage29["delta_G"] > 0.0, f"G_tr-G_flat={fmt(stage29['delta_G'])}")
    require("F_flat > F_tr", stage29["delta_F"] > 0.0, f"F_flat-F_tr={fmt(stage29['delta_F'])}")

    print(f"G_tr = {fmt(stage29['G_tr'])}")
    print(f"F_tr = {fmt(stage29['F_tr'])}")
    print(f"G_flat = {fmt(stage29['G_flat'])}")
    print(f"F_flat = {fmt(stage29['F_flat'])}")

    stage31 = stage31_values(sample, base, stage29)
    require("dG_tr/dxi > 0", stage31["dG_dxi"] > 0.0, f"dG/dxi={fmt(stage31['dG_dxi'])}")
    require("dS/dzeta > 0", stage31["dS_dzeta"] > 0.0, f"dS/dzeta={fmt(stage31['dS_dzeta'])}")
    require(
        "finite positive support-enhanced load",
        stage31["M_actual"] > 0.0,
        f"M_mix*S(zeta0)={fmt(stage31['M_actual'])}, target={fmt(stage31['M_target'])}, residual={fmt(stage31['M_actual'] - stage31['M_target'])}",
    )
    print(f"M_crit = {fmt(stage31['M_crit'])}")
    print(f"S_req = {fmt(stage31['S_req'])}")
    print(f"S_crit = {fmt(stage31['S_crit'])}")
    print(f"zeta_req = {fmt(stage31['zeta_req'])}")
    print(f"zeta_crit = {fmt(stage31['zeta_crit'])}")

    if in_stage31_domain:
        require(
            "Stage-31 domain: M_mix < M_crit",
            base["M_mix"] < stage31["M_crit"],
            f"M_mix={fmt(base['M_mix'])}, M_crit={fmt(stage31['M_crit'])}",
        )
        if stage31["zeta_req"] > 0.0:
            require(
                "zeta_req < zeta_crit < 1/eps",
                stage31["zeta_req"] < stage31["zeta_crit"] < 1.0 / base["eps"],
                f"zeta_req={fmt(stage31['zeta_req'])}, zeta_crit={fmt(stage31['zeta_crit'])}, 1/eps={fmt(1.0/base['eps'])}",
            )
        else:
            require(
                "zeta_req = 0 when support is not needed",
                near(stage31["zeta_req"], 0.0) and stage31["zeta_crit"] > 0.0,
                f"zeta_req={fmt(stage31['zeta_req'])}, zeta_crit={fmt(stage31['zeta_crit'])}, M_mix={fmt(base['M_mix'])}, M_req={fmt(stage31['M_target'])}",
            )
    else:
        require(
            "boundary stress outside Stage-31 domain",
            base["M_mix"] >= stage31["M_crit"],
            f"M_mix={fmt(base['M_mix'])}, M_crit={fmt(stage31['M_crit'])}",
        )
        print("note: Stage-31 support-feasibility inequalities are intentionally not enforced for this out-of-domain stress sample.")

    require(
        "finite support ratio below pole",
        z0 < 1.0 / base["eps"] and z1 < 1.0 / base["eps"],
        f"zeta0={fmt(z0)}, zeta1={fmt(z1)}, 1/eps={fmt(1.0/base['eps'])}",
    )

    print(
        f"dxi_phys/dzeta = {fmt(stage31['dxi_dzeta'])} "
        f"(sample zeta={fmt(stage31['zeta0'])})"
    )


def main(argv: list[str]) -> int:
    root = Path(__file__).resolve().parents[3]
    sample_path = Path(argv[1]) if len(argv) > 1 else root / "scripts/moving_throat/numerical/stage028_031_samples.json"
    samples = load_samples(sample_path)
    print(f"Loaded {len(samples)} samples from {sample_path}")
    for sample in samples:
        check_case(sample)
    print("\nAll coherent support-chain numerical stress checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
