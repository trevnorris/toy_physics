#!/usr/bin/env python3
"""Stress harness for the full-profile mouth correction chain (Stages 133-135)."""

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
    if data.get("schema") != "moving_throat_numerical_stage133_135_v1":
        raise ValueError("unexpected sample schema")
    return data


class ProfileCorrectionHarness:
    def __init__(self, config: dict) -> None:
        mp.mp.dps = 80
        const = config["constants"]
        tol = config["tolerances"]
        self.value_tol = float(tol["value_tol"])
        self.derivative_tol = float(tol["derivative_tol"])
        self.curvature_tol = float(tol["curvature_tol"])

        self.Pi_star_expected = mp.mpf(str(const["Pi_star_expected"]))
        self.Sigma_m_star_expected = mp.mpf(str(const["Sigma_m_star_expected"]))
        self.g_star_expected = mp.mpf(str(const["g_star_expected"]))
        self.S_star_expected = mp.mpf(str(const["S_star_expected"]))
        self.gprime_star_expected = mp.mpf(str(const["gprime_star_expected"]))
        self.A_T_expected = mp.mpf(str(const["A_T_expected"]))
        self.B_T_expected = mp.mpf(str(const["B_T_expected"]))
        self.T_hat_star_expected = mp.mpf(str(const["T_hat_star_expected"]))
        self.cov_cR_expected = mp.mpf(str(const["cov_cR_expected"]))
        self.cov_KR_expected = mp.mpf(str(const["cov_KR_expected"]))
        self.delta_g_expected = mp.mpf(str(const["delta_g_expected"]))
        self.delta_S_expected = mp.mpf(str(const["delta_S_expected"]))
        self.deltaPi_expected = mp.mpf(str(const["deltaPi_expected"]))
        self.deltaT_expected = mp.mpf(str(const["deltaT_expected"]))
        self.g1_expected = mp.mpf(str(const["g1_expected"]))
        self.S1_expected = mp.mpf(str(const["S1_expected"]))
        self.Pi1_expected = mp.mpf(str(const["Pi1_expected"]))
        self.T1_expected = mp.mpf(str(const["T1_expected"]))

        self.fd_steps = [mp.mpf(str(x)) for x in config["finite_difference_steps"]]
        self.curvature_steps = [mp.mpf(str(x)) for x in config["curvature_steps"]]
        self.lambda_samples = config["lambda_samples"]

        self.kappa = mp.pi / 2
        self.rF1 = mp.sqrt(4107 - 100 * mp.pi**2) / (10 * mp.pi)
        self.gMinus = self.rF1 - mp.sqrt(1 + self.rF1**2) / 2

        self.Pi_star = mp.findroot(lambda p: self.g_formula(p) - self.gMinus, mp.mpf("1.5"))
        self.S_star = self.s_formula(self.Pi_star)
        self.g_star = self.g_formula(self.Pi_star)
        self.gprime_star = mp.diff(self.g_formula, self.Pi_star)
        self.sprime_star = mp.diff(self.s_formula, self.Pi_star)
        self.Sigma_m_star = self.Pi_star / (4 - self.S_star)
        self.T_hat_star = mp.sqrt(mp.mpf(9) * (self.Pi_star / (1 - self.S_star / 4)) / 20)
        self.A_T = -(
            mp.mpf(9)
            / (40 * self.T_hat_star)
            * (
                1 / (self.gprime_star * (1 - self.S_star / 4))
                + self.Pi_star * self.sprime_star / (4 * self.gprime_star * (1 - self.S_star / 4) ** 2)
            )
        )
        self.B_T = (mp.mpf(9) / (40 * self.T_hat_star)) * self.Pi_star / (4 * (1 - self.S_star / 4) ** 2)

    def g_formula(self, p: mp.mpf) -> mp.mpf:
        return 2 * p * (2 * p * mp.e**p + mp.pi) / ((4 * p**2 + mp.pi**2) * (mp.e**p - 1))

    def s_formula(self, p: mp.mpf) -> mp.mpf:
        return p * (
            self.kappa * mp.tanh(self.kappa) + p * (mp.e**(-p) * mp.sech(self.kappa) - 1)
        ) / ((1 - mp.e**(-p)) * (self.kappa**2 - p**2))

    def t_s(self, x: mp.mpf) -> mp.mpf:
        return (1 - mp.e**(-self.Pi_star * x)) / (self.Pi_star * (1 - mp.e**(-self.Pi_star))) - x * mp.e**(-self.Pi_star) / (
            1 - mp.e**(-self.Pi_star)
        )

    def t_q(self, x: mp.mpf) -> mp.mpf:
        c_q = self.Pi_star / ((1 - mp.e**(-self.Pi_star)) * (self.kappa**2 - self.Pi_star**2))
        a_q = c_q * (self.kappa * mp.sinh(self.kappa) + self.Pi_star * mp.e**(-self.Pi_star)) / (
            self.kappa * mp.cosh(self.kappa)
        )
        return a_q * mp.sinh(self.kappa * x) - c_q * mp.cosh(self.kappa * x) + c_q * mp.e**(-self.Pi_star * x)

    def sigma_star(self, x: mp.mpf) -> mp.mpf:
        return self.Pi_star * mp.e**(-self.Pi_star * x) / (1 - mp.e**(-self.Pi_star))

    def r_star(self, x: mp.mpf) -> mp.mpf:
        return self.Sigma_m_star * (4 * self.t_s(x) - self.t_q(x)) - self.Pi_star * x

    def c_kernel(self, x: mp.mpf) -> mp.mpf:
        return mp.cos(mp.pi * x / 2)

    def s_kernel(self, x: mp.mpf) -> mp.mpf:
        return mp.cosh((mp.pi / 2) * (1 - x)) / mp.cosh(mp.pi / 2)

    def nint(self, fn) -> mp.mpf:
        return mp.quad(fn, [0, 1])

    def avg(self, fn) -> mp.mpf:
        return self.nint(lambda t: self.sigma_star(t) * fn(t))

    def nonlinear_sigma(self, lam: mp.mpf):
        def unnormalized(x: mp.mpf) -> mp.mpf:
            return mp.e ** (-self.Pi_star * x - lam * self.r_star(x))

        z = self.nint(unnormalized)
        return lambda x: unnormalized(x) / z

    def check_stage133(self) -> None:
        print("=== Stage 133: full-profile residual geometry ===")
        require(
            "Pi_* matches expected canonical bias",
            near(self.Pi_star, self.Pi_star_expected, self.value_tol),
            f"Pi_*={fmt(self.Pi_star)}, expected={fmt(self.Pi_star_expected)}",
        )
        require(
            "Sigma_m* matches expected canonical gain",
            near(self.Sigma_m_star, self.Sigma_m_star_expected, self.value_tol),
            f"Sigma_m*={fmt(self.Sigma_m_star)}, expected={fmt(self.Sigma_m_star_expected)}",
        )
        require(
            "g_* matches expected canonical overlap",
            near(self.g_star, self.g_star_expected, self.value_tol),
            f"g_*={fmt(self.g_star)}, expected={fmt(self.g_star_expected)}",
        )
        require(
            "S_* matches expected canonical response",
            near(self.S_star, self.S_star_expected, self.value_tol),
            f"S_*={fmt(self.S_star)}, expected={fmt(self.S_star_expected)}",
        )

        require("R_*(0)=0", abs(self.r_star(mp.mpf("0"))) <= self.value_tol, f"R(0)={fmt(self.r_star(mp.mpf('0')))}")
        for step in self.fd_steps:
            forward = (self.r_star(step) - self.r_star(mp.mpf("0"))) / step
            print(f"forward derivative at h={fmt(step)} -> {fmt(forward)}")
            require(
                f"R_*'(0) small at h={fmt(step)}",
                abs(forward) <= self.derivative_tol,
                f"forward diff={fmt(forward)}",
            )
            require(
                f"R_*(h) is negative at h={fmt(step)}",
                self.r_star(step) < 0,
                f"R(h)={fmt(self.r_star(step))}",
            )

        target_r2 = -3 * self.Sigma_m_star * self.Pi_star / (1 - mp.e**(-self.Pi_star))
        for step in self.curvature_steps:
            second = (self.r_star(step) - 2 * self.r_star(mp.mpf("0")) + self.r_star(-step)) / (step**2)
            print(f"central second derivative at h={fmt(step)} -> {fmt(second)}")
            require(
                f"R_*''(0) matches target at h={fmt(step)}",
                abs(second - target_r2) <= self.curvature_tol,
                f"finite diff={fmt(second)}, target={fmt(target_r2)}",
            )
            require(
                f"R_*''(0) stays negative at h={fmt(step)}",
                second < 0,
                f"finite diff={fmt(second)}",
            )

    def check_stage135(self) -> None:
        print("\n=== Stage 135: covariance correction and nonlinear stress ===")
        avg_r = self.avg(self.r_star)
        avg_c = self.avg(self.c_kernel)
        avg_s = self.avg(self.s_kernel)

        cov_cR = self.nint(lambda t: self.sigma_star(t) * (self.c_kernel(t) - avg_c) * (self.r_star(t) - avg_r))
        cov_sR = self.nint(lambda t: self.sigma_star(t) * (self.s_kernel(t) - avg_s) * (self.r_star(t) - avg_r))

        delta_g = -cov_cR
        delta_s = -cov_sR
        delta_pi = cov_cR / self.gprime_star
        delta_t = self.A_T * delta_g + self.B_T * delta_s

        require(
            "Cov(c,R_*) matches Stage 135",
            near(cov_cR, self.cov_cR_expected, self.value_tol),
            f"cov_cR={fmt(cov_cR)}, expected={fmt(self.cov_cR_expected)}",
        )
        require(
            "Cov(K_q,R_*) matches Stage 135",
            near(cov_sR, self.cov_KR_expected, self.value_tol),
            f"cov_KR={fmt(cov_sR)}, expected={fmt(self.cov_KR_expected)}",
        )
        require(
            "delta g_act matches Stage 135",
            near(delta_g, self.delta_g_expected, self.value_tol),
            f"delta_g={fmt(delta_g)}, expected={fmt(self.delta_g_expected)}",
        )
        require(
            "delta S_act matches Stage 135",
            near(delta_s, self.delta_S_expected, self.value_tol),
            f"delta_S={fmt(delta_s)}, expected={fmt(self.delta_S_expected)}",
        )
        require(
            "delta Pi_act matches Stage 135",
            near(delta_pi, self.deltaPi_expected, self.value_tol),
            f"deltaPi={fmt(delta_pi)}, expected={fmt(self.deltaPi_expected)}",
        )
        require(
            "delta T_act matches Stage 135",
            near(delta_t, self.deltaT_expected, self.value_tol),
            f"deltaT={fmt(delta_t)}, expected={fmt(self.deltaT_expected)}",
        )

        pi_errors = []
        t_errors = []
        for case in self.lambda_samples:
            lam = mp.mpf(str(case["lambda"]))
            sigma = self.nonlinear_sigma(lam)
            g_lam = self.nint(lambda t: sigma(t) * self.c_kernel(t))
            s_lam = self.nint(lambda t: sigma(t) * self.s_kernel(t))
            delta_g_actual = g_lam - self.g_star
            delta_s_actual = s_lam - self.S_star
            delta_pi_actual = -delta_g_actual / self.gprime_star
            delta_t_actual = self.A_T * delta_g_actual + self.B_T * delta_s_actual

            delta_g_linear = lam * delta_g
            delta_s_linear = lam * delta_s
            delta_pi_linear = lam * delta_pi
            delta_t_linear = lam * delta_t

            pi_rel = abs(delta_pi_actual - delta_pi_linear) / abs(delta_pi_linear)
            t_rel = abs(delta_t_actual - delta_t_linear) / abs(delta_t_linear)
            pi_errors.append(pi_rel)
            t_errors.append(t_rel)

            print(f"\n{case['name']} ({case['kind']}): lambda={fmt(lam)}")
            print(
                f"  delta_g_actual={fmt(delta_g_actual)}, delta_g_linear={fmt(delta_g_linear)}, "
                f"delta_S_actual={fmt(delta_s_actual)}, delta_S_linear={fmt(delta_s_linear)}"
            )
            print(
                f"  deltaPi_actual={fmt(delta_pi_actual)}, deltaPi_linear={fmt(delta_pi_linear)}, "
                f"deltaT_actual={fmt(delta_t_actual)}, deltaT_linear={fmt(delta_t_linear)}"
            )

            require(
                f"{case['name']} keeps the correction direction",
                delta_g_actual < 0 and delta_s_actual < 0 and delta_pi_actual > 0 and delta_t_actual > 0,
                f"delta_g={fmt(delta_g_actual)}, delta_S={fmt(delta_s_actual)}, deltaPi={fmt(delta_pi_actual)}, deltaT={fmt(delta_t_actual)}",
            )
            require(
                f"{case['name']} deltaPi stays within the linear envelope",
                pi_rel <= float(case["relative_envelope"]),
                f"relative error={fmt(pi_rel)}, envelope={fmt(case['relative_envelope'])}",
            )
            require(
                f"{case['name']} deltaT stays within the linear envelope",
                t_rel <= float(case["relative_envelope"]),
                f"relative error={fmt(t_rel)}, envelope={fmt(case['relative_envelope'])}",
            )

            if lam == 1:
                require(
                    "lambda=1 g_1 matches Stage 135",
                    near(g_lam, self.g1_expected, self.value_tol),
                    f"g_1={fmt(g_lam)}, expected={fmt(self.g1_expected)}",
                )
                require(
                    "lambda=1 S_1 matches Stage 135",
                    near(s_lam, self.S1_expected, self.value_tol),
                    f"S_1={fmt(s_lam)}, expected={fmt(self.S1_expected)}",
                )
                require(
                    "lambda=1 Pi_1 matches Stage 135",
                    near(self.Pi_star + delta_pi_actual, self.Pi1_expected, self.value_tol),
                    f"Pi_1={fmt(self.Pi_star + delta_pi_actual)}, expected={fmt(self.Pi1_expected)}",
                )
                require(
                    "lambda=1 T_1 matches Stage 135",
                    near(self.T_hat_star + delta_t_actual, self.T1_expected, self.value_tol),
                    f"T_1={fmt(self.T_hat_star + delta_t_actual)}, expected={fmt(self.T1_expected)}",
                )

        require(
            "deltaPi linearization error grows with lambda",
            all(b > a for a, b in zip(pi_errors[:-1], pi_errors[1:], strict=True)),
            f"errors={[fmt(err) for err in pi_errors]}",
        )
        require(
            "deltaT linearization error grows with lambda",
            all(b > a for a, b in zip(t_errors[:-1], t_errors[1:], strict=True)),
            f"errors={[fmt(err) for err in t_errors]}",
        )


def main(argv: list[str]) -> int:
    root = Path(__file__).resolve().parents[3]
    config_path = Path(argv[1]) if len(argv) > 1 else root / "scripts/moving_throat/numerical/stage133_135_profile_correction_samples.json"
    config = load_config(config_path)
    h = ProfileCorrectionHarness(config)
    print(f"Loaded config from {config_path}")
    h.check_stage133()
    h.check_stage135()
    print("\nAll stage-133/135 profile-correction stress checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
