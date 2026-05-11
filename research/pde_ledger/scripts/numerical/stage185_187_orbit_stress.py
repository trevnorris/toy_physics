#!/usr/bin/env python3
"""Numerical stress harness for the invariant/orbit closure chain (Stages 168-170)."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np


ORDER = [
    "lambda_W",
    "c_etaU",
    "gamma",
    "K_U",
    "K_eta_eff",
    "K_W_eff",
    "mu_W",
    "T_U",
]

FREE_ORDER = ["Delta_lambda", "Delta_c", "Delta_gamma", "Delta_U", "Delta_W"]
FULL_ORDER = [
    "Delta_lambda",
    "Delta_c",
    "Delta_gamma",
    "Delta_U",
    "Delta_eta",
    "Delta_W",
    "Delta_mu",
    "Delta_T",
]


def fmt(value: float) -> str:
    return f"{value:.12g}"


def require(label: str, condition: bool, detail: str) -> None:
    status = "PASS" if condition else "FAIL"
    print(f"[{status}] {label}: {detail}")
    if not condition:
        raise AssertionError(label)


def load_config(path: Path) -> dict:
    data = json.loads(path.read_text())
    if data.get("schema") != "moving_throat_numerical_stage185_187_v1":
        raise ValueError("unexpected sample schema")
    return data


def vector_from_mapping(mapping: dict, order: list[str]) -> np.ndarray:
    return np.array([float(mapping[key]) for key in order], dtype=float)


def to_mapping(vector: np.ndarray, order: list[str]) -> dict[str, float]:
    return {key: float(value) for key, value in zip(order, vector, strict=True)}


def near(lhs: float, rhs: float, tol: float) -> bool:
    return abs(lhs - rhs) <= tol * (1.0 + abs(rhs))


class OrbitHarness:
    def __init__(self, config: dict) -> None:
        tol = config["tolerances"]
        self.stage185_step = float(config["stage185_step"])
        self.stage186_step = float(config["stage186_step"])
        self.monomial_linear_tol = float(tol["monomial_linear_tol"])
        self.observable_linear_tol = float(tol["observable_linear_tol"])
        self.tangent_zero_tol = float(tol["tangent_zero_tol"])
        self.kernel_tol = float(tol["kernel_tol"])
        self.invariant_tol = float(tol["invariant_tol"])
        self.transverse_floor = float(tol["transverse_floor"])

    def state_vector(self, case: dict) -> np.ndarray:
        return vector_from_mapping(case["reference_state"], ORDER)

    def derived(self, state: np.ndarray) -> dict[str, float]:
        lam, ceta, gamma, k_u, k_eta, k_w, mu, t_u = state
        chi = gamma * ceta / k_u
        delta = t_u / k_u
        eps_w = gamma * gamma * lam * lam / (k_u * k_w)
        eps = eps_w * (1.0 - (2.0 / 11.0) * delta / (1.0 + delta))
        zratio = lam * lam * mu / (k_eta * k_w * k_w)
        eps_eta = ceta * ceta / (k_u * k_eta)
        e_star = 2.0 * eps_w / (1.0 - eps) * (11.0 + 9.0 * delta) / (11.0 * (1.0 + delta))
        f_star = 2.0 * chi / (1.0 + delta) + 4.0 * eps_w * delta / (
            11.0 * (1.0 - eps) * (1.0 + delta) * (1.0 + delta)
        )
        c_tr = chi * delta / ((1.0 + chi) * (1.0 + delta) * (1.0 + chi + delta))
        a_tr = 2.0 * chi / ((1.0 + chi) * (1.0 + delta))
        b_star = 2.0 * (1.0 + chi + delta) / delta
        c_star = ((1.0 + chi) * (1.0 + delta) * (1.0 + chi + delta)) / (chi * delta)
        r_tr = (1.0 + chi / (1.0 + delta)) / (1.0 + chi)
        t_sq = zratio * (1.0 + chi) * (1.0 + chi) / ((1.0 - eps) * (1.0 - eps))
        t_star = r_tr ** (-c_star)
        n_star = t_sq * (r_tr ** b_star)
        return {
            "chi": chi,
            "delta": delta,
            "eps_w": eps_w,
            "eps": eps,
            "zratio": zratio,
            "eps_eta": eps_eta,
            "E": e_star,
            "F": f_star,
            "C_tr": c_tr,
            "A_tr": a_tr,
            "B_star": b_star,
            "C_star": c_star,
            "R_tr": r_tr,
            "T_sq": t_sq,
            "T_star": t_star,
            "N_star": n_star,
        }

    def monomial_invariants(self, state: np.ndarray, ref: dict[str, float]) -> dict[str, float]:
        d = self.derived(state)
        ctr = (d["chi"] ** (1.0 + ref["delta"])) * (d["delta"] ** (1.0 + ref["chi"]))
        cnt = d["zratio"] * (d["eps_w"] ** ref["E"]) * (d["delta"] ** (-ref["F"]))
        return {"Ctr": ctr, "Cnt": cnt, "eps_eta": d["eps_eta"]}

    def observable_invariants(self, state: np.ndarray, ref: dict[str, float]) -> dict[str, float]:
        d = self.derived(state)
        return {
            "T_star": d["R_tr"] ** (-ref["C_star"]),
            "N_star": d["T_sq"] * (d["R_tr"] ** ref["B_star"]),
            "eps_eta": d["eps_eta"],
        }

    def matrix(self, ref: dict[str, float]) -> np.ndarray:
        chi = ref["chi"]
        delta = ref["delta"]
        e_star = ref["E"]
        f_star = ref["F"]
        return np.array(
            [
                [0.0, 1.0 + delta, 1.0 + delta, -(2.0 + chi + delta), 0.0, 0.0, 0.0, 1.0 + chi],
                [2.0 * (1.0 + e_star), 0.0, 2.0 * e_star, f_star - e_star, -1.0, -(2.0 + e_star), 1.0, -f_star],
                [0.0, 2.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0],
            ],
            dtype=float,
        )

    def selected_minor(self, ref: dict[str, float]) -> np.ndarray:
        chi = ref["chi"]
        f_star = ref["F"]
        return np.array(
            [
                [0.0, 0.0, 1.0 + chi],
                [-1.0, 1.0, -f_star],
                [-1.0, 0.0, 0.0],
            ],
            dtype=float,
        )

    def full_tangent_from_free(self, ref: dict[str, float], free: dict) -> np.ndarray:
        lam, ceta, gamma, k_u, k_w = [float(free[key]) for key in FREE_ORDER]
        chi = ref["chi"]
        delta = ref["delta"]
        e_star = ref["E"]
        f_star = ref["F"]
        deta = 2.0 * ceta - k_u
        dt = k_u - ((1.0 + delta) / (1.0 + chi)) * (gamma + ceta - k_u)
        dmu = (
            deta
            + 2.0 * k_w
            - 2.0 * lam
            - e_star * (2.0 * gamma + 2.0 * lam - k_u - k_w)
            - f_star * ((1.0 + delta) / (1.0 + chi)) * (gamma + ceta - k_u)
        )
        return np.array([lam, ceta, gamma, k_u, deta, k_w, dmu, dt], dtype=float)

    def linear_monomial_drifts(self, ref: dict[str, float], drift: np.ndarray) -> dict[str, float]:
        lam, ceta, gamma, k_u, k_eta, k_w, mu, t_u = drift
        sigma_chi = gamma + ceta - k_u
        sigma_delta = t_u - k_u
        sigma_eps = 2.0 * gamma + 2.0 * lam - k_u - k_w
        sigma_z = 2.0 * lam + mu - k_eta - 2.0 * k_w
        sigma_eta = 2.0 * ceta - k_u - k_eta
        sigma_tr = (1.0 + ref["chi"]) * sigma_delta + (1.0 + ref["delta"]) * sigma_chi
        sigma_nt = sigma_z + ref["E"] * sigma_eps - ref["F"] * sigma_delta
        return {
            "Sigma_tr": sigma_tr,
            "Sigma_nt": sigma_nt,
            "Sigma_eta": sigma_eta,
        }

    def linear_observable_defect(self, ref: dict[str, float], drifts: dict[str, float]) -> dict[str, float]:
        theta = -ref["C_tr"] * drifts["Sigma_tr"]
        xi = ref["A_tr"] * drifts["Sigma_tr"] + drifts["Sigma_nt"]
        r_plus_xi = -(ref["eps_eta"] / (1.0 - ref["eps_eta"])) * drifts["Sigma_eta"]
        return {"Theta_1": theta, "Xi_1": xi, "R1_plus_Xi1": r_plus_xi}

    def stepped_state(self, base: np.ndarray, drift: np.ndarray, step: float) -> np.ndarray:
        return base * np.exp(step * drift)

    def stage185_compare(self, case: dict, ref_state: np.ndarray, ref: dict[str, float]) -> None:
        print("\n-- Stage 185: observable vs microscopic drifts --")
        ref_obs = self.observable_invariants(ref_state, ref)
        ref_mon = self.monomial_invariants(ref_state, ref)

        tangent_drifts = [self.full_tangent_from_free(ref, item) for item in case["tangent_free_vectors"]]
        transverse_drifts = [vector_from_mapping(item, FULL_ORDER) for item in case["transverse_full_vectors"]]
        all_vectors = [("tangent", item["name"], drift) for item, drift in zip(case["tangent_free_vectors"], tangent_drifts, strict=True)]
        all_vectors += [("transverse", item["name"], drift) for item, drift in zip(case["transverse_full_vectors"], transverse_drifts, strict=True)]

        for kind, name, drift in all_vectors:
            step_state = self.stepped_state(ref_state, drift, self.stage185_step)
            obs = self.observable_invariants(step_state, ref)
            mon = self.monomial_invariants(step_state, ref)
            obs_drifts = {
                "Sigma_tr": math.log(obs["T_star"] / ref_obs["T_star"]) / self.stage185_step,
                "Sigma_nt": math.log(obs["N_star"] / ref_obs["N_star"]) / self.stage185_step,
                "Sigma_eta": math.log(obs["eps_eta"] / ref_obs["eps_eta"]) / self.stage185_step,
            }
            mon_drifts = {
                "Sigma_tr": math.log(mon["Ctr"] / ref_mon["Ctr"]) / self.stage185_step,
                "Sigma_nt": math.log(mon["Cnt"] / ref_mon["Cnt"]) / self.stage185_step,
                "Sigma_eta": math.log(mon["eps_eta"] / ref_mon["eps_eta"]) / self.stage185_step,
            }
            lin = self.linear_monomial_drifts(ref, drift)
            print(
                f"{kind} {name}: "
                f"obs=({fmt(obs_drifts['Sigma_tr'])}, {fmt(obs_drifts['Sigma_nt'])}, {fmt(obs_drifts['Sigma_eta'])}) "
                f"mon=({fmt(mon_drifts['Sigma_tr'])}, {fmt(mon_drifts['Sigma_nt'])}, {fmt(mon_drifts['Sigma_eta'])})"
            )
            for key in ("Sigma_tr", "Sigma_nt", "Sigma_eta"):
                require(
                    f"{name} monomial drift matches linear {key}",
                    near(mon_drifts[key], lin[key], self.monomial_linear_tol),
                    f"direct={fmt(mon_drifts[key])}, linear={fmt(lin[key])}",
                )
                require(
                    f"{name} observable drift matches monomial {key}",
                    near(obs_drifts[key], mon_drifts[key], self.observable_linear_tol),
                    f"observable={fmt(obs_drifts[key])}, monomial={fmt(mon_drifts[key])}",
                )
            defects = self.linear_observable_defect(ref, mon_drifts)
            print(
                f"  defects: Theta_1={fmt(defects['Theta_1'])}, "
                f"Xi_1={fmt(defects['Xi_1'])}, R_1+Xi_1={fmt(defects['R1_plus_Xi1'])}"
            )
            if kind == "transverse":
                max_drift = max(abs(mon_drifts[key]) for key in ("Sigma_tr", "Sigma_nt", "Sigma_eta"))
                require(
                    f"{name} is genuinely transverse",
                    max_drift > self.transverse_floor,
                    f"max monomial drift={fmt(max_drift)}",
                )

    def stage186_check(self, case: dict, ref_state: np.ndarray, ref: dict[str, float]) -> None:
        print("\n-- Stage 186: tangent-space kernel checks --")
        matrix = self.matrix(ref)
        minor = self.selected_minor(ref)
        det_minor = float(np.linalg.det(minor))
        svals = np.linalg.svd(matrix, compute_uv=False)
        minor_svals = np.linalg.svd(minor, compute_uv=False)
        cond_matrix = float(svals[0] / svals[-1])
        cond_minor = float(minor_svals[0] / minor_svals[-1])
        print(f"det selected minor = {fmt(det_minor)}")
        print(f"cond(M_*) = {fmt(cond_matrix)}")
        print(f"cond(selected minor) = {fmt(cond_minor)}")
        require(
            "selected minor determinant = 1 + chi_*",
            near(det_minor, 1.0 + ref["chi"], self.kernel_tol),
            f"det={fmt(det_minor)}, expected={fmt(1.0 + ref['chi'])}",
        )

        ref_mon = self.monomial_invariants(ref_state, ref)
        for item in case["tangent_free_vectors"]:
            drift = self.full_tangent_from_free(ref, item)
            residual = matrix @ drift
            max_residual = float(np.max(np.abs(residual)))
            require(
                f"{item['name']} lies in ker(M_*)",
                max_residual <= self.kernel_tol,
                f"max residual={fmt(max_residual)}",
            )
            step_state = self.stepped_state(ref_state, drift, self.stage186_step)
            mon = self.monomial_invariants(step_state, ref)
            mon_drifts = [
                math.log(mon["Ctr"] / ref_mon["Ctr"]) / self.stage186_step,
                math.log(mon["Cnt"] / ref_mon["Cnt"]) / self.stage186_step,
                math.log(mon["eps_eta"] / ref_mon["eps_eta"]) / self.stage186_step,
            ]
            max_drift = max(abs(value) for value in mon_drifts)
            require(
                f"{item['name']} preserves the invariant triple to first order",
                max_drift <= self.tangent_zero_tol,
                f"max drift={fmt(max_drift)}",
            )

    def stage187_check(self, case: dict, ref_state: np.ndarray, ref: dict[str, float]) -> None:
        print("\n-- Stage 187: finite orbit/quotient checks --")
        ref_mon = self.monomial_invariants(ref_state, ref)
        for item in case["finite_orbit_pairs"]:
            delta = self.full_tangent_from_free(ref, item)
            paired_state = self.stepped_state(ref_state, delta, 1.0)
            recovered = np.log(paired_state / ref_state)
            max_recover = float(np.max(np.abs(recovered - delta)))
            require(
                f"{item['name']} recovers its finite log-ratio vector",
                max_recover <= self.invariant_tol,
                f"max recovery error={fmt(max_recover)}",
            )

            lam, ceta, gamma, k_u, _, k_w, _, _ = delta
            deta_expected = 2.0 * ceta - k_u
            dt_expected = k_u - ((1.0 + ref["delta"]) / (1.0 + ref["chi"])) * (gamma + ceta - k_u)
            dmu_expected = (
                deta_expected
                + 2.0 * k_w
                - 2.0 * lam
                - ref["E"] * (2.0 * gamma + 2.0 * lam - k_u - k_w)
                - ref["F"] * ((1.0 + ref["delta"]) / (1.0 + ref["chi"])) * (gamma + ceta - k_u)
            )
            require(
                f"{item['name']} Delta_eta finite law",
                near(delta[4], deta_expected, self.invariant_tol),
                f"direct={fmt(delta[4])}, expected={fmt(deta_expected)}",
            )
            require(
                f"{item['name']} Delta_T finite law",
                near(delta[7], dt_expected, self.invariant_tol),
                f"direct={fmt(delta[7])}, expected={fmt(dt_expected)}",
            )
            require(
                f"{item['name']} Delta_mu finite law",
                near(delta[6], dmu_expected, self.invariant_tol),
                f"direct={fmt(delta[6])}, expected={fmt(dmu_expected)}",
            )

            pair_mon = self.monomial_invariants(paired_state, ref)
            ratio_logs = {
                "Ctr": math.log(pair_mon["Ctr"] / ref_mon["Ctr"]),
                "Cnt": math.log(pair_mon["Cnt"] / ref_mon["Cnt"]),
                "eps_eta": math.log(pair_mon["eps_eta"] / ref_mon["eps_eta"]),
            }
            for key, value in ratio_logs.items():
                require(
                    f"{item['name']} preserves {key}",
                    abs(value) <= self.invariant_tol,
                    f"log ratio={fmt(value)}",
                )

    def run_case(self, case: dict) -> None:
        print(f"\n=== {case['name']} ===")
        ref_state = self.state_vector(case)
        ref = self.derived(ref_state)
        print(
            "reference coefficients: "
            f"chi0*={fmt(ref['chi'])}, deltaU*={fmt(ref['delta'])}, "
            f"epsilon_W*={fmt(ref['eps_w'])}, epsilon*={fmt(ref['eps'])}, "
            f"E*={fmt(ref['E'])}, F*={fmt(ref['F'])}, epsilon_eta*={fmt(ref['eps_eta'])}"
        )
        self.stage185_compare(case, ref_state, ref)
        self.stage186_check(case, ref_state, ref)
        self.stage187_check(case, ref_state, ref)


def main(argv: list[str]) -> int:
    default_config = Path(__file__).resolve().with_name("stage185_187_orbit_samples.json")
    config_path = (
        Path(argv[1])
        if len(argv) > 1
        else default_config
    )
    config = load_config(config_path)
    harness = OrbitHarness(config)
    print(f"Loaded config from {config_path}")
    for case in config["cases"]:
        harness.run_case(case)
    print("\nAll stage-185/187 orbit-conditioning stress checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
