#!/usr/bin/env python3
"""Finite-deformation stress harness for the Stage 129-130 rigidity bridge."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import mpmath as mp


def fmt(value: float | mp.mpf) -> str:
    return mp.nstr(value, 12)


def near(lhs: float | mp.mpf, rhs: float | mp.mpf, tol: float) -> bool:
    return abs(lhs - rhs) <= tol * (1.0 + abs(rhs))


def require(label: str, condition: bool, detail: str) -> None:
    status = "PASS" if condition else "FAIL"
    print(f"[{status}] {label}: {detail}")
    if not condition:
        raise AssertionError(label)


def load_config(path: Path) -> dict:
    data = json.loads(path.read_text())
    if data.get("schema") != "moving_throat_numerical_stage129_130_v1":
        raise ValueError("unexpected sample schema")
    return data


class RigidityHarness:
    def __init__(self, config: dict) -> None:
        mp.mp.dps = 80
        self.config = config
        const = config["constants"]
        tol = config["tolerances"]
        self.value_tol = float(tol["value_tol"])
        self.integral_tol = float(tol["integral_tol"])
        self.compensation_residual_tol = mp.mpf(str(tol["compensation_residual_tol"]))

        self.Pi_star_expected = mp.mpf(str(const["Pi_star_expected"]))
        self.g_star_expected = mp.mpf(str(const["g_star_expected"]))
        self.S_star_expected = mp.mpf(str(const["S_star_expected"]))
        self.gprime_star_expected = mp.mpf(str(const["gprime_star_expected"]))
        self.Sprime_star_expected = mp.mpf(str(const["Sprime_star_expected"]))
        self.T_hat_star_expected = mp.mpf(str(const["T_hat_star_expected"]))
        self.A_T_expected = mp.mpf(str(const["A_T_expected"]))
        self.B_T_expected = mp.mpf(str(const["B_T_expected"]))
        self.lambda_pi_zero_expected = mp.mpf(str(const["lambda_pi_zero_expected"]))
        self.lambda_t_zero_expected = mp.mpf(str(const["lambda_T_zero_expected"]))

        self.kappa = mp.pi / 2
        self.rF1 = mp.sqrt(4107 - 100 * mp.pi**2) / (10 * mp.pi)
        self.g_minus = self.rF1 - mp.sqrt(1 + self.rF1**2) / 2

        self.Pi_star = mp.findroot(lambda p: self.g_formula(p) - self.g_minus, mp.mpf("1.5"))
        self.g_star = self.g_formula(self.Pi_star)
        self.S_star = self.S_formula(self.Pi_star)
        self.gprime_star = mp.diff(self.g_formula, self.Pi_star)
        self.Sprime_star = mp.diff(self.S_formula, self.Pi_star)
        self.T_hat_star = mp.sqrt(mp.mpf(9) * (self.Pi_star / (1 - self.S_star / 4)) / 20)
        self.A_T = -(
            mp.mpf(9)
            / (40 * self.T_hat_star)
            * (
                1 / (self.gprime_star * (1 - self.S_star / 4))
                + self.Pi_star * self.Sprime_star / (4 * self.gprime_star * (1 - self.S_star / 4) ** 2)
            )
        )
        self.B_T = (mp.mpf(9) / (40 * self.T_hat_star)) * self.Pi_star / (4 * (1 - self.S_star / 4) ** 2)

        self.g_uniform = 2 / mp.pi
        self.S_uniform = 2 * mp.tanh(mp.pi / 2) / mp.pi
        self.g_derivative = mp.pi / 4
        self.S_derivative = (1 + mp.sinh(mp.pi / 2)) / (2 * mp.cosh(mp.pi / 2))

    def g_formula(self, p: mp.mpf) -> mp.mpf:
        return 2 * p * (2 * p * mp.e**p + mp.pi) / ((4 * p**2 + mp.pi**2) * (mp.e**p - 1))

    def S_formula(self, p: mp.mpf) -> mp.mpf:
        return p * (
            self.kappa * mp.tanh(self.kappa) + p * (mp.e**(-p) * mp.sech(self.kappa) - 1)
        ) / ((1 - mp.e**(-p)) * (self.kappa**2 - p**2))

    def sigma_exp(self, x: mp.mpf, p: mp.mpf) -> mp.mpf:
        return p * mp.e ** (-p * x) / (1 - mp.e**(-p))

    def varsigma_lambda(self, x: mp.mpf, lam: mp.mpf) -> mp.mpf:
        return (1 - lam) * 1 + lam * (mp.pi / 2) * mp.cos(mp.pi * x / 2)

    def g_lambda_formula(self, lam: mp.mpf) -> mp.mpf:
        return (1 - lam) * self.g_uniform + lam * self.g_derivative

    def S_lambda_formula(self, lam: mp.mpf) -> mp.mpf:
        return (1 - lam) * self.S_uniform + lam * self.S_derivative

    def sigma_total(self, x: mp.mpf, p: mp.mpf, eps: mp.mpf, lam: mp.mpf) -> mp.mpf:
        return (1 - eps) * self.sigma_exp(x, p) + eps * self.varsigma_lambda(x, lam)

    def nint(self, fn) -> mp.mpf:
        return mp.quad(fn, [0, 1])

    def check_base_constants(self) -> None:
        print("=== Stage 129-130: base constants ===")
        require(
            "Pi_* matches audited value",
            near(self.Pi_star, self.Pi_star_expected, self.value_tol),
            f"Pi_*={fmt(self.Pi_star)}, expected={fmt(self.Pi_star_expected)}",
        )
        require(
            "g_* matches audited value",
            near(self.g_star, self.g_star_expected, self.value_tol),
            f"g_*={fmt(self.g_star)}, expected={fmt(self.g_star_expected)}",
        )
        require(
            "S_* matches audited value",
            near(self.S_star, self.S_star_expected, self.value_tol),
            f"S_*={fmt(self.S_star)}, expected={fmt(self.S_star_expected)}",
        )
        require(
            "g_*' matches audited value",
            near(self.gprime_star, self.gprime_star_expected, self.value_tol),
            f"g_*'={fmt(self.gprime_star)}, expected={fmt(self.gprime_star_expected)}",
        )
        require(
            "S_*' matches audited value",
            near(self.Sprime_star, self.Sprime_star_expected, self.value_tol),
            f"S_*'={fmt(self.Sprime_star)}, expected={fmt(self.Sprime_star_expected)}",
        )
        require(
            "T_* matches audited value",
            near(self.T_hat_star, self.T_hat_star_expected, self.value_tol),
            f"T_*={fmt(self.T_hat_star)}, expected={fmt(self.T_hat_star_expected)}",
        )
        require(
            "A_T matches audited value",
            near(self.A_T, self.A_T_expected, self.value_tol),
            f"A_T={fmt(self.A_T)}, expected={fmt(self.A_T_expected)}",
        )
        require(
            "B_T matches audited value",
            near(self.B_T, self.B_T_expected, self.value_tol),
            f"B_T={fmt(self.B_T)}, expected={fmt(self.B_T_expected)}",
        )

    def run_case(self, case: dict) -> None:
        name = case["name"]
        lam = mp.mpf(str(case["lambda"]))
        print(f"\n=== {name} ({case['kind']}) ===")
        print(f"lambda = {fmt(lam)}")
        if "assumptions" in case:
            print("assumptions:")
            for item in case["assumptions"]:
                print(f"  - {item}")

        if name == "bias_neutral_direction":
            require(
                "bias-neutral lambda matches Stage 131",
                near(lam, self.lambda_pi_zero_expected, self.value_tol),
                f"lambda={fmt(lam)}, expected={fmt(self.lambda_pi_zero_expected)}",
            )
        if name == "traction_neutral_direction":
            require(
                "traction-neutral lambda matches Stage 131",
                near(lam, self.lambda_t_zero_expected, self.value_tol),
                f"lambda={fmt(lam)}, expected={fmt(self.lambda_t_zero_expected)}",
            )

        norm_var = self.nint(lambda x: self.varsigma_lambda(x, lam))
        g_var_int = self.nint(lambda x: self.varsigma_lambda(x, lam) * mp.cos(mp.pi * x / 2))
        S_var_int = self.nint(lambda x: self.varsigma_lambda(x, lam) * mp.cosh((mp.pi / 2) * (1 - x)) / mp.cosh(mp.pi / 2))
        g_var = self.g_lambda_formula(lam)
        S_var = self.S_lambda_formula(lam)
        require("varsigma_lambda normalization", near(norm_var, 1, self.integral_tol), f"integral={fmt(norm_var)}")
        require("varsigma_lambda is nonnegative at x=1", self.varsigma_lambda(mp.mpf("1"), lam) >= -1e-30, f"varsigma(1)={fmt(self.varsigma_lambda(mp.mpf('1'), lam))}")
        require("g_lambda integral matches formula", near(g_var_int, g_var, self.integral_tol), f"integral={fmt(g_var_int)}, formula={fmt(g_var)}")
        require("S_lambda integral matches formula", near(S_var_int, S_var, self.integral_tol), f"integral={fmt(S_var_int)}, formula={fmt(S_var)}")

        dPi_per_eps = -(g_var - self.g_star) / self.gprime_star
        dS_per_eps = (S_var - self.S_star) - (self.Sprime_star / self.gprime_star) * (g_var - self.g_star)
        dT_per_eps = self.A_T * (g_var - self.g_star) + self.B_T * (S_var - self.S_star)

        print(
            f"linear dPi/eps={fmt(dPi_per_eps)}, "
            f"linear dS/eps={fmt(dS_per_eps)}, linear dT/eps={fmt(dT_per_eps)}"
        )

        pi_errors: list[mp.mpf] = []
        s_errors: list[mp.mpf] = []
        t_errors: list[mp.mpf] = []
        quadratic_coeffs: list[mp.mpf] = []

        for sample in case["samples"]:
            eps = mp.mpf(str(sample["epsilon"]))
            predicted_dPi = eps * dPi_per_eps
            predicted_dS = eps * dS_per_eps
            predicted_dT = eps * dT_per_eps

            p_exact = mp.findroot(
                lambda p: (1 - eps) * self.g_formula(p) + eps * g_var - self.g_star,
                self.Pi_star + predicted_dPi,
            )
            sigma_norm = self.nint(lambda x: self.sigma_total(x, p_exact, eps, lam))
            g_int = self.nint(lambda x: self.sigma_total(x, p_exact, eps, lam) * mp.cos(mp.pi * x / 2))
            S_int = self.nint(
                lambda x: self.sigma_total(x, p_exact, eps, lam) * mp.cosh((mp.pi / 2) * (1 - x)) / mp.cosh(mp.pi / 2)
            )
            compensation_residual = (1 - eps) * self.g_formula(p_exact) + eps * g_var - self.g_star
            t_exact = mp.sqrt(mp.mpf(9) * p_exact / (20 * (1 - S_int / 4)))

            delta_pi_exact = p_exact - self.Pi_star
            delta_s_exact = S_int - self.S_star
            delta_t_exact = t_exact - self.T_hat_star

            print(f"\nepsilon = {fmt(eps)}")
            print(
                f"  Pi_exact={fmt(p_exact)}, deltaPi_exact={fmt(delta_pi_exact)}, "
                f"deltaPi_linear={fmt(predicted_dPi)}"
            )
            print(
                f"  deltaS_exact={fmt(delta_s_exact)}, deltaS_linear={fmt(predicted_dS)}, "
                f"deltaT_exact={fmt(delta_t_exact)}, deltaT_linear={fmt(predicted_dT)}"
            )

            require("sigma_total normalization", near(sigma_norm, 1, self.integral_tol), f"integral={fmt(sigma_norm)}")
            require("sigma_total stays nonnegative at x=1", self.sigma_total(mp.mpf('1'), p_exact, eps, lam) >= -1e-30, f"sigma_total(1)={fmt(self.sigma_total(mp.mpf('1'), p_exact, eps, lam))}")
            require("g_total integral matches g_*", near(g_int, self.g_star, self.integral_tol), f"g_int={fmt(g_int)}, g_*={fmt(self.g_star)}")
            require("S_total integral matches affine formula", near(S_int, (1 - eps) * self.S_formula(p_exact) + eps * S_var, self.integral_tol), f"S_int={fmt(S_int)}")
            require("exact compensation residual", abs(compensation_residual) <= self.compensation_residual_tol, f"residual={fmt(compensation_residual)}")

            if "deltaPi_abs_tol" in sample:
                require(
                    f"epsilon={fmt(eps)} bias shift stays neutral",
                    abs(delta_pi_exact) <= mp.mpf(str(sample["deltaPi_abs_tol"])),
                    f"deltaPi_exact={fmt(delta_pi_exact)}",
                )
            else:
                pi_rel = abs(delta_pi_exact - predicted_dPi) / abs(predicted_dPi)
                pi_errors.append(pi_rel)
                require(
                    f"epsilon={fmt(eps)} deltaPi stays within first-order envelope",
                    pi_rel <= mp.mpf(str(sample["deltaPi_rel_envelope"])),
                    f"relative error={fmt(pi_rel)}, envelope={fmt(sample['deltaPi_rel_envelope'])}",
                )

            s_rel = abs(delta_s_exact - predicted_dS) / abs(predicted_dS)
            s_errors.append(s_rel)
            require(
                f"epsilon={fmt(eps)} deltaS stays within first-order envelope",
                s_rel <= mp.mpf(str(sample["deltaS_rel_envelope"])),
                f"relative error={fmt(s_rel)}, envelope={fmt(sample['deltaS_rel_envelope'])}",
            )

            if "deltaT_quadratic_coeff_max" in sample:
                quad_coeff = abs(delta_t_exact) / (eps**2)
                quadratic_coeffs.append(quad_coeff)
                require(
                    f"epsilon={fmt(eps)} deltaT remains quadratic-small",
                    quad_coeff <= mp.mpf(str(sample["deltaT_quadratic_coeff_max"])),
                    f"|deltaT|/eps^2={fmt(quad_coeff)}, max={fmt(sample['deltaT_quadratic_coeff_max'])}",
                )
            else:
                t_rel = abs(delta_t_exact - predicted_dT) / abs(predicted_dT)
                t_errors.append(t_rel)
                require(
                    f"epsilon={fmt(eps)} deltaT stays within first-order envelope",
                    t_rel <= mp.mpf(str(sample["deltaT_rel_envelope"])),
                    f"relative error={fmt(t_rel)}, envelope={fmt(sample['deltaT_rel_envelope'])}",
                )

        if len(pi_errors) > 1:
            require(
                "deltaPi relative error grows with epsilon",
                all(b > a for a, b in zip(pi_errors[:-1], pi_errors[1:], strict=True)),
                f"errors={[fmt(err) for err in pi_errors]}",
            )
        if len(s_errors) > 1:
            require(
                "deltaS relative error grows with epsilon",
                all(b > a for a, b in zip(s_errors[:-1], s_errors[1:], strict=True)),
                f"errors={[fmt(err) for err in s_errors]}",
            )
        if len(t_errors) > 1:
            require(
                "deltaT relative error grows with epsilon",
                all(b > a for a, b in zip(t_errors[:-1], t_errors[1:], strict=True)),
                f"errors={[fmt(err) for err in t_errors]}",
            )
        if len(quadratic_coeffs) > 1:
            require(
                "quadratic traction coefficient grows mildly with epsilon",
                all(b >= a for a, b in zip(quadratic_coeffs[:-1], quadratic_coeffs[1:], strict=True)),
                f"coeffs={[fmt(coeff) for coeff in quadratic_coeffs]}",
            )


def main(argv: list[str]) -> int:
    root = Path(__file__).resolve().parents[3]
    config_path = Path(argv[1]) if len(argv) > 1 else root / "scripts/moving_throat/numerical/stage129_130_rigidity_samples.json"
    config = load_config(config_path)
    harness = RigidityHarness(config)
    print(f"Loaded config from {config_path}")
    harness.check_base_constants()
    for case in config["cases"]:
        harness.run_case(case)
    print("\nAll stage-129/130 rigidity stress checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
