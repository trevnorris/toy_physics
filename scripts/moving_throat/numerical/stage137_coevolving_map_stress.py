#!/usr/bin/env python3
"""Numerical bridge harness for the Stage 137 co-evolving core-mouth map."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np


def fmt(value: float) -> str:
    return f"{value:.12g}"


def near(lhs: float, rhs: float, tol: float) -> bool:
    return abs(lhs - rhs) <= tol * (1.0 + abs(rhs))


def require(label: str, condition: bool, detail: str) -> None:
    status = "PASS" if condition else "FAIL"
    print(f"[{status}] {label}: {detail}")
    if not condition:
        raise AssertionError(label)


def load_config(path: Path) -> dict:
    data = json.loads(path.read_text())
    if data.get("schema") != "moving_throat_numerical_stage137_v1":
        raise ValueError("unexpected sample schema")
    return data


class CoevolvingMapHarness:
    def __init__(self, config: dict) -> None:
        const = config["constants"]
        tol = config["tolerances"]

        self.rF1 = float(const["rF1"])
        self.g_star = float(const["g_star"])
        self.S_star = float(const["S_star"])
        self.Pi_star = float(const["Pi_star"])
        self.Sigma0_star = float(const["Sigma0_star"])
        self.T_hat_star = float(const["T_hat_star"])
        self.Sigma0_can = float(const["Sigma0_can"])
        self.S_can = float(const["S_can"])
        self.Pi_can = float(const["Pi_can"])
        self.T_hat_can = float(const["T_hat_can"])
        self.g_fp_star = float(const["g_fp_star"])
        self.S_fp_star = float(const["S_fp_star"])
        self.R_fp_star = float(const["R_fp_star"])
        self.Pi_fp_star = float(const["Pi_fp_star"])

        self.normalization_tol = float(tol["normalization_tol"])
        self.positivity_tol = float(tol["positivity_tol"])
        self.moment_tol = float(tol["moment_tol"])
        self.fixed_point_residual_tol = float(tol["fixed_point_residual_tol"])
        self.kernel_slope_tol = float(tol["kernel_slope_tol"])
        self.phi_slope_tol = float(tol["phi_slope_tol"])

        self.grid_points = int(config["grid_points"])
        self.max_iterations = int(config["max_iterations"])
        self.fixed_point_tol = float(config["fixed_point_tol"])
        self.bisection_iterations = int(config["bisection_iterations"])
        self.fd_steps = [float(x) for x in config["fd_steps"]]
        self.representative_profiles = config["representative_profiles"]
        self.transport_cases = config["transport_cases"]

        self.x = np.linspace(0.0, 1.0, self.grid_points)
        dx = self.x[1] - self.x[0]
        self.w = np.full(self.grid_points, dx)
        self.w[0] = dx / 2.0
        self.w[-1] = dx / 2.0
        self.kappa = math.pi / 2.0
        self.c_grid = np.cos(math.pi * self.x / 2.0)
        self.Kq_grid = np.cosh(self.kappa * (1.0 - self.x)) / math.cosh(self.kappa)
        self.sqrt_term = math.sqrt(1.0 + self.rF1**2)
        self._canonical_profile = self.normalize(self.Pi_star * np.exp(-self.Pi_star * self.x))

    def normalize(self, raw: np.ndarray) -> np.ndarray:
        return raw / float(np.sum(raw * self.w))

    def profile_min(self, sigma: np.ndarray) -> float:
        return float(np.min(sigma))

    def g(self, sigma: np.ndarray) -> float:
        return float(np.sum(self.c_grid * sigma * self.w))

    def S(self, sigma: np.ndarray) -> float:
        return float(np.sum(self.Kq_grid * sigma * self.w))

    def R(self, sigma: np.ndarray) -> float:
        gv = self.g(sigma)
        return ((gv - self.rF1) ** 2) / (1.0 + self.rF1**2)

    def sigma_var(self, lam: float) -> np.ndarray:
        raw = (1.0 - lam) * np.ones_like(self.x) + lam * (math.pi / 2.0) * np.cos(math.pi * self.x / 2.0)
        return self.normalize(raw)

    def seed_profile(self, kind: str, lam: float | None = None) -> np.ndarray:
        if kind == "canonical_exponential":
            return self._canonical_profile.copy()
        if kind == "uniform":
            return self.normalize(np.ones_like(self.x))
        if kind == "derivative_match":
            return self.normalize((math.pi / 2.0) * np.cos(math.pi * self.x / 2.0))
        if kind == "convex_blend":
            if lam is None:
                raise ValueError("convex_blend requires lambda")
            return self.sigma_var(lam)
        raise ValueError(f"unknown seed kind: {kind}")

    def Ts(self, sigma: np.ndarray) -> np.ndarray:
        cum_sig = np.cumsum(sigma * self.w)
        cum_y = np.cumsum(self.x * sigma * self.w)
        return cum_y + self.x * (cum_sig[-1] - cum_sig)

    def Tq(self, sigma: np.ndarray) -> np.ndarray:
        a_term = np.cumsum(np.sinh(self.kappa * self.x) * sigma * self.w)
        b_term = np.cumsum((np.cosh(self.kappa * (1.0 - self.x)) * sigma * self.w)[::-1])[::-1]
        return (
            np.cosh(self.kappa * (1.0 - self.x)) * a_term
            + np.sinh(self.kappa * self.x) * b_term
        ) / (self.kappa * math.cosh(self.kappa))

    def phi_grid(self, sigma: np.ndarray, sigma0: float) -> np.ndarray:
        return sigma0 * (self.Ts(sigma) - self.R(sigma) * self.Tq(sigma))

    def next_sigma(self, sigma: np.ndarray, sigma0: float) -> np.ndarray:
        phi = self.phi_grid(sigma, sigma0)
        phi_shift = phi - float(np.min(phi))
        return self.normalize(np.exp(-phi_shift))

    def solve_fixed_point(self, sigma0: float, initial: np.ndarray) -> tuple[np.ndarray, int, float]:
        sigma = initial.copy()
        err = math.inf
        for niter in range(self.max_iterations):
            sigma_new = self.next_sigma(sigma, sigma0)
            err = float(np.max(np.abs(sigma_new - sigma)))
            sigma = sigma_new
            if err < self.fixed_point_tol:
                return sigma, niter + 1, err
        raise RuntimeError(f"fixed-point solve did not converge for sigma0={sigma0}")

    def ts_at(self, sigma: np.ndarray, x_value: float) -> float:
        return float(np.sum(np.minimum(x_value, self.x) * sigma * self.w))

    def tq_at(self, sigma: np.ndarray, x_value: float) -> float:
        xmin = np.minimum(x_value, self.x)
        xmax = np.maximum(x_value, self.x)
        kernel = np.sinh(self.kappa * xmin) * np.cosh(self.kappa * (1.0 - xmax)) / (self.kappa * math.cosh(self.kappa))
        return float(np.sum(kernel * sigma * self.w))

    def phi_at(self, sigma: np.ndarray, sigma0: float, x_value: float) -> float:
        return sigma0 * (self.ts_at(sigma, x_value) - self.R(sigma) * self.tq_at(sigma, x_value))

    def forward_slope(self, fn, h: float) -> float:
        return (4.0 * fn(h) - fn(2.0 * h)) / (2.0 * h)

    def check_slopes(self, sigma: np.ndarray, sigma0: float, prefix: str) -> float:
        S_value = self.S(sigma)
        R_value = self.R(sigma)
        pi_expected = sigma0 * (1.0 - R_value * S_value)
        phi_slopes: list[float] = []

        for h in self.fd_steps:
            ts_slope = self.forward_slope(lambda z: self.ts_at(sigma, z), h)
            tq_slope = self.forward_slope(lambda z: self.tq_at(sigma, z), h)
            phi_slope = self.forward_slope(lambda z: self.phi_at(sigma, sigma0, z), h)
            phi_slopes.append(phi_slope)

            require(
                f"{prefix} T_s'(0) at h={fmt(h)}",
                near(ts_slope, 1.0, self.kernel_slope_tol),
                f"slope={fmt(ts_slope)}",
            )
            require(
                f"{prefix} T_q'(0)=S at h={fmt(h)}",
                near(tq_slope, S_value, self.kernel_slope_tol),
                f"slope={fmt(tq_slope)}, S={fmt(S_value)}",
            )
            require(
                f"{prefix} Phi'(0)=Sigma0(1-RS) at h={fmt(h)}",
                near(phi_slope, pi_expected, self.phi_slope_tol),
                f"slope={fmt(phi_slope)}, expected={fmt(pi_expected)}",
            )

        avg_phi = float(sum(phi_slopes) / len(phi_slopes))
        require(
            f"{prefix} forward-slope self-consistency",
            near(avg_phi, pi_expected, self.phi_slope_tol),
            f"avg slope={fmt(avg_phi)}, expected={fmt(pi_expected)}",
        )
        return avg_phi

    def exact_shift_value(self, g_value: float) -> float:
        dg = g_value - self.g_star
        return 0.25 - dg / self.sqrt_term + dg * dg / (1.0 + self.rF1**2)

    def representative_profile(self, case: dict) -> np.ndarray:
        kind = case["kind"]
        if kind == "explicit_canonical":
            return self.seed_profile("canonical_exponential")
        if kind == "derivative_match":
            return self.seed_profile("derivative_match")
        if kind == "convex_blend":
            return self.seed_profile("convex_blend", float(case["lambda"]))
        if kind == "fixed_point":
            seed_kind = case.get("seed_kind", "canonical_exponential")
            initial = self.seed_profile(seed_kind)
            sigma, _, _ = self.solve_fixed_point(float(case["sigma0"]), initial)
            return sigma
        raise ValueError(f"unknown representative profile kind: {kind}")

    def check_representative_profiles(self) -> None:
        print("=== Stage 137: representative profile checks ===")
        for case in self.representative_profiles:
            sigma0 = float(case["sigma0"])
            sigma = self.representative_profile(case)
            g_value = self.g(sigma)
            S_value = self.S(sigma)
            R_value = self.R(sigma)
            norm = float(np.sum(sigma * self.w))
            min_sigma = self.profile_min(sigma)
            exact_shift = self.exact_shift_value(g_value)
            pi_direct = self.check_slopes(sigma, sigma0, case["name"])
            print(f"\n=== {case['name']} ({case['kind']}) ===")
            print(f"sigma0 = {fmt(sigma0)}")
            if "assumptions" in case:
                print("assumptions:")
                for item in case["assumptions"]:
                    print(f"  - {item}")
            print(
                f"g={fmt(g_value)}, S={fmt(S_value)}, R={fmt(R_value)}, "
                f"Pi_direct={fmt(pi_direct)}, T_hat={fmt(math.sqrt(9.0 * sigma0 / 20.0))}"
            )

            require(
                f"{case['name']} normalization",
                abs(norm - 1.0) <= self.normalization_tol,
                f"integral={fmt(norm)}",
            )
            require(
                f"{case['name']} positivity",
                min_sigma >= -self.positivity_tol,
                f"min sigma={fmt(min_sigma)}",
            )
            require(
                f"{case['name']} exact R(g) shift law",
                near(R_value, exact_shift, self.moment_tol),
                f"R={fmt(R_value)}, exact shift={fmt(exact_shift)}",
            )

            if "expected_g" in case:
                require(
                    f"{case['name']} expected g",
                    near(g_value, float(case["expected_g"]), self.moment_tol),
                    f"g={fmt(g_value)}, expected={fmt(float(case['expected_g']))}",
                )
            if "expected_S" in case:
                require(
                    f"{case['name']} expected S",
                    near(S_value, float(case["expected_S"]), self.moment_tol),
                    f"S={fmt(S_value)}, expected={fmt(float(case['expected_S']))}",
                )
            if "expected_R" in case:
                require(
                    f"{case['name']} expected R",
                    near(R_value, float(case["expected_R"]), self.moment_tol),
                    f"R={fmt(R_value)}, expected={fmt(float(case['expected_R']))}",
                )
            if "expected_Pi" in case:
                require(
                    f"{case['name']} expected Pi",
                    near(pi_direct, float(case["expected_Pi"]), self.phi_slope_tol),
                    f"Pi={fmt(pi_direct)}, expected={fmt(float(case['expected_Pi']))}",
                )
            if case.get("expect_compensated"):
                require(
                    f"{case['name']} hits g_*",
                    near(g_value, self.g_star, self.moment_tol),
                    f"g={fmt(g_value)}, g_*={fmt(self.g_star)}",
                )
                require(
                    f"{case['name']} hits R=1/4",
                    near(R_value, 0.25, self.moment_tol),
                    f"R={fmt(R_value)}",
                )
            if case["kind"] == "fixed_point":
                residual = float(np.max(np.abs(self.next_sigma(sigma, sigma0) - sigma)))
                require(
                    f"{case['name']} fixed-point residual",
                    residual <= self.fixed_point_residual_tol,
                    f"L_inf residual={fmt(residual)}",
                )

    def check_transport(self) -> None:
        print("\n=== Stage 137: first-order transport from the canonical branch ===")
        sigma_star = self.seed_profile("canonical_exponential")
        Pi_star_direct = self.check_slopes(sigma_star, self.Sigma0_star, "canonical_transport_base")
        require(
            "canonical transport base Pi matches Pi_*",
            near(Pi_star_direct, self.Pi_star, self.phi_slope_tol),
            f"Pi_direct={fmt(Pi_star_direct)}, Pi_*={fmt(self.Pi_star)}",
        )

        for case in self.transport_cases:
            lam = float(case["lambda"])
            sigma0 = float(case["sigma0"])
            sigma_var = self.sigma_var(lam)
            g_var = self.g(sigma_var)
            S_var = self.S(sigma_var)
            dg_rate = g_var - self.g_star
            dS_rate = S_var - self.S_star
            dR_rate = -dg_rate / self.sqrt_term
            dPi_rate = -sigma0 * (0.25 * dS_rate + self.S_star * dR_rate)

            print(f"\n=== {case['name']} ({case['kind']}) ===")
            print(f"lambda = {fmt(lam)}")
            print(
                f"direction rates: dg/dε={fmt(dg_rate)}, dS/dε={fmt(dS_rate)}, "
                f"dR/dε={fmt(dR_rate)}, dPi/dε={fmt(dPi_rate)}"
            )

            for sample in case["samples"]:
                eps = float(sample["epsilon"])
                sigma_eps = self.normalize((1.0 - eps) * sigma_star + eps * sigma_var)
                norm = float(np.sum(sigma_eps * self.w))
                min_sigma = self.profile_min(sigma_eps)
                g_eps = self.g(sigma_eps)
                S_eps = self.S(sigma_eps)
                R_eps = self.R(sigma_eps)
                exact_shift = self.exact_shift_value(g_eps)
                Pi_eps = self.check_slopes(sigma_eps, sigma0, f"{case['name']} ε={fmt(eps)}")
                deltaPi_exact = Pi_eps - Pi_star_direct
                deltaPi_linear = eps * dPi_rate
                deltaPi_rel = abs(deltaPi_exact - deltaPi_linear) / abs(deltaPi_linear)

                print(f"epsilon = {fmt(eps)}")
                print(
                    f"  g={fmt(g_eps)}, S={fmt(S_eps)}, R={fmt(R_eps)}, "
                    f"deltaPi_exact={fmt(deltaPi_exact)}, deltaPi_linear={fmt(deltaPi_linear)}"
                )

                require(
                    f"{case['name']} ε={fmt(eps)} normalization",
                    abs(norm - 1.0) <= self.normalization_tol,
                    f"integral={fmt(norm)}",
                )
                require(
                    f"{case['name']} ε={fmt(eps)} positivity",
                    min_sigma >= -self.positivity_tol,
                    f"min sigma={fmt(min_sigma)}",
                )
                require(
                    f"{case['name']} ε={fmt(eps)} exact R(g) shift law",
                    near(R_eps, exact_shift, self.moment_tol),
                    f"R={fmt(R_eps)}, exact shift={fmt(exact_shift)}",
                )
                require(
                    f"{case['name']} ε={fmt(eps)} first-order deltaPi envelope",
                    deltaPi_rel <= float(sample["deltaPi_rel_envelope"]),
                    f"relative error={fmt(deltaPi_rel)}, envelope={fmt(float(sample['deltaPi_rel_envelope']))}",
                )


def main() -> int:
    config_path = Path(__file__).with_name("stage137_coevolving_map_samples.json")
    harness = CoevolvingMapHarness(load_config(config_path))
    harness.check_representative_profiles()
    harness.check_transport()
    print("\nStage 137 co-evolving map stress harness passed.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover - terminal harness
        print(f"\n[FAIL] {exc}")
        raise
