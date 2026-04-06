#!/usr/bin/env python3
"""Numerical stress harness for the explicit Family-1 mouth branch (Stages 125-127)."""

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
    if data.get("schema") != "moving_throat_numerical_stage125_127_v1":
        raise ValueError("unexpected sample schema")
    return data


class MouthBranchHarness:
    def __init__(self, config: dict) -> None:
        const = config["constants"]
        self.rF1 = float(const["rF1"])
        self.g_minus_expected = float(const["g_minus_expected"])
        self.g_plus_expected = float(const["g_plus_expected"])
        self.Pi_star_expected = float(const["Pi_star_expected"])
        self.Pi_match_expected = float(const["Pi_match_expected"])
        self.T_hat_star_expected = float(const["T_hat_star_expected"])
        self.T_hat_match_expected = float(const["T_hat_match_expected"])
        self.R_inf_expected = float(const["R_inf_expected"])
        self.denom_inf_expected = float(const["denom_inf_expected"])
        self.that_sqrtPi_limit_expected = float(const["that_sqrtPi_limit_expected"])

        self.quadrature_tol = float(config["tolerances"]["quadrature_tol"])
        self.normalization_tol = float(config["tolerances"]["normalization_tol"])
        self.residual_tol = float(config["tolerances"]["residual_tol"])
        self.root_tol = float(config["tolerances"]["root_tol"])
        self.ratio_tol = float(config["tolerances"]["ratio_tol"])
        self.bisection_iterations = int(config["bisection_iterations"])

        n = int(config["quadrature_points"])
        self.x = np.linspace(0.0, 1.0, n)
        self.c_grid = np.cos(math.pi * self.x / 2.0)
        self.Kq_grid = np.cosh((math.pi / 2.0) * (1.0 - self.x)) / math.cosh(math.pi / 2.0)
        self.sample_points = config["sample_points"]
        self.monotone_scan_points = [float(x) for x in config["monotone_scan_points"]]
        self.root_checks = config["root_checks"]
        self.finite_exclusion_points = [float(x) for x in config["finite_exclusion_points"]]
        self.asymptotic_samples = [float(x) for x in config["asymptotic_samples"]]

    def trapz(self, values: np.ndarray) -> float:
        return float(np.trapz(values, self.x))

    def sigma_profile(self, pi_value: float) -> np.ndarray:
        raw = pi_value * np.exp(-pi_value * self.x) / (1.0 - math.exp(-pi_value))
        return raw / self.trapz(raw)

    def g_formula(self, pi_value: float) -> float:
        return 2.0 * pi_value * (2.0 * pi_value * math.exp(pi_value) + math.pi) / (
            (4.0 * pi_value * pi_value + math.pi * math.pi) * (math.exp(pi_value) - 1.0)
        )

    def S_formula(self, pi_value: float) -> float:
        numerator = pi_value * (
            (math.pi / 2.0) * math.tanh(math.pi / 2.0)
            + pi_value * (math.exp(-pi_value) / math.cosh(math.pi / 2.0) - 1.0)
        )
        denominator = (1.0 - math.exp(-pi_value)) * (((math.pi / 2.0) ** 2) - pi_value * pi_value)
        return numerator / denominator

    def R_from_g(self, g_value: float) -> float:
        return ((g_value - self.rF1) ** 2) / (1.0 + self.rF1**2)

    def sigma0_formula(self, pi_value: float) -> float:
        gv = self.g_formula(pi_value)
        sv = self.S_formula(pi_value)
        return pi_value / (1.0 - self.R_from_g(gv) * sv)

    def that_formula(self, pi_value: float) -> float:
        return math.sqrt(9.0 * self.sigma0_formula(pi_value) / 20.0)

    def bisection_root(self, target: float, lo: float, hi: float) -> float:
        flo = self.g_formula(lo) - target
        fhi = self.g_formula(hi) - target
        require(
            f"root bracket sign change on [{fmt(lo)}, {fmt(hi)}]",
            flo == 0.0 or fhi == 0.0 or flo * fhi < 0.0,
            f"f(lo)={fmt(flo)}, f(hi)={fmt(fhi)}",
        )
        left = lo
        right = hi
        for _ in range(self.bisection_iterations):
            mid = 0.5 * (left + right)
            fm = self.g_formula(mid) - target
            if flo == 0.0:
                return left
            if fhi == 0.0:
                return right
            if flo * fm <= 0.0:
                right = mid
                fhi = fm
            else:
                left = mid
                flo = fm
        return 0.5 * (left + right)

    def stage125_samples(self) -> None:
        print("=== Stages 125-126: explicit positive-mouth branch samples ===")
        for case in self.sample_points:
            pi_value = float(case["Pi"])
            sigma = self.sigma_profile(pi_value)
            norm = self.trapz(sigma)
            min_sigma = float(np.min(sigma))
            g_quad = self.trapz(sigma * self.c_grid)
            S_quad = self.trapz(sigma * self.Kq_grid)
            g_exact = self.g_formula(pi_value)
            S_exact = self.S_formula(pi_value)
            R_quad = self.R_from_g(g_quad)
            sigma0 = self.sigma0_formula(pi_value)
            denom_exact = 1.0 - self.R_from_g(g_exact) * S_exact
            that = self.that_formula(pi_value)
            residual = pi_value - sigma0 * (1.0 - R_quad * S_quad)

            print(f"\n{case['name']} ({case['kind']}): Pi={fmt(pi_value)}")
            if "assumptions" in case:
                print("assumptions:")
                for item in case["assumptions"]:
                    print(f"  - {item}")
            print(
                f"  g_quad={fmt(g_quad)}, g_exact={fmt(g_exact)}, "
                f"S_quad={fmt(S_quad)}, S_exact={fmt(S_exact)}"
            )
            print(
                f"  R={fmt(R_quad)}, denom={fmt(denom_exact)}, "
                f"Sigma0={fmt(sigma0)}, T_hat={fmt(that)}"
            )

            require(
                f"{case['name']} positive source normalization",
                abs(norm - 1.0) <= self.normalization_tol,
                f"integral={fmt(norm)}",
            )
            require(
                f"{case['name']} source stays nonnegative",
                min_sigma >= -1e-15,
                f"min sigma={fmt(min_sigma)}",
            )
            require(
                f"{case['name']} g quadrature matches formula",
                near(g_quad, g_exact, self.quadrature_tol),
                f"quadrature={fmt(g_quad)}, formula={fmt(g_exact)}",
            )
            require(
                f"{case['name']} S quadrature matches formula",
                near(S_quad, S_exact, self.quadrature_tol),
                f"quadrature={fmt(S_quad)}, formula={fmt(S_exact)}",
            )
            require(
                f"{case['name']} self-consistent residual",
                abs(residual) <= self.residual_tol,
                f"residual={fmt(residual)}",
            )
            require(
                f"{case['name']} stays in 2/pi < g < 1",
                2.0 / math.pi < g_exact < 1.0,
                f"2/pi={fmt(2.0 / math.pi)}, g={fmt(g_exact)}",
            )
            require(
                f"{case['name']} denominator stays positive",
                denom_exact > 0.8,
                f"1-R_q S_q={fmt(denom_exact)}",
            )

        require(
            "canonical T_hat matches audited value",
            near(self.that_formula(self.Pi_star_expected), self.T_hat_star_expected, self.root_tol),
            f"T_hat(Pi_*)={fmt(self.that_formula(self.Pi_star_expected))}, expected={fmt(self.T_hat_star_expected)}",
        )
        require(
            "derivative-match T_hat matches audited value",
            near(self.that_formula(self.Pi_match_expected), self.T_hat_match_expected, self.root_tol),
            f"T_hat(Pi_match)={fmt(self.that_formula(self.Pi_match_expected))}, expected={fmt(self.T_hat_match_expected)}",
        )

    def stage127_branch_selection(self) -> None:
        print("\n=== Stage 127: branch selection and finite-bias exclusion ===")
        scan_values = [self.g_formula(pi_value) for pi_value in self.monotone_scan_points]
        for pi_value, g_value in zip(self.monotone_scan_points, scan_values, strict=True):
            print(f"scan Pi={fmt(pi_value)} -> g(Pi)={fmt(g_value)}")

        diffs = [b - a for a, b in zip(scan_values[:-1], scan_values[1:], strict=True)]
        require(
            "g(Pi) increases on the explicit branch scan",
            all(diff > 0.0 for diff in diffs),
            f"min diff={fmt(min(diffs))}",
        )

        roots: dict[str, float] = {}
        for block in self.root_checks:
            target = math.pi / 4.0 if block["target_key"] == "pi_over_4" else getattr(self, block["target_key"])
            block_roots = []
            for lo, hi in block["brackets"]:
                root = self.bisection_root(target, float(lo), float(hi))
                block_roots.append(root)
                print(f"{block['name']} root from [{fmt(float(lo))}, {fmt(float(hi))}] = {fmt(root)}")
            ref_root = block_roots[0]
            for idx, root in enumerate(block_roots[1:], start=1):
                require(
                    f"{block['name']} bracket {idx} agrees with bracket 0",
                    near(root, ref_root, self.root_tol),
                    f"root={fmt(root)}, ref={fmt(ref_root)}",
                )
            require(
                f"{block['name']} root matches expected value",
                near(ref_root, float(block["expected_pi"]), self.root_tol),
                f"root={fmt(ref_root)}, expected={fmt(float(block['expected_pi']))}",
            )
            roots[block["name"]] = ref_root

        require(
            "derivative-match point lies above canonical point",
            roots["derivative_match_branch"] > roots["lower_compensated_branch"],
            f"Pi_match={fmt(roots['derivative_match_branch'])}, Pi_*={fmt(roots['lower_compensated_branch'])}",
        )

        exclusion_values = [self.g_formula(pi_value) for pi_value in self.finite_exclusion_points]
        require(
            "finite explicit branch never reaches g=1",
            all(g_value < 1.0 for g_value in exclusion_values),
            f"max g(Pi) on finite window={fmt(max(exclusion_values))}",
        )
        require(
            "finite explicit branch never reaches g_+",
            all(g_value < self.g_plus_expected for g_value in exclusion_values),
            f"max g(Pi)={fmt(max(exclusion_values))}, g_+={fmt(self.g_plus_expected)}",
        )
        require(
            "upper compensated branch remains outside the positive-source range",
            self.g_plus_expected > 1.0,
            f"g_+={fmt(self.g_plus_expected)}",
        )

    def stage126_asymptotics(self) -> None:
        print("\n=== Stage 126: singular-limit asymptotics ===")
        denom_errors = []
        ratio_errors = []
        g_tail = []
        for pi_value in self.asymptotic_samples:
            g_value = self.g_formula(pi_value)
            S_value = self.S_formula(pi_value)
            R_value = self.R_from_g(g_value)
            denom = 1.0 - R_value * S_value
            ratio = self.that_formula(pi_value) / math.sqrt(pi_value)
            g_tail.append(g_value)
            denom_errors.append(abs(denom - self.denom_inf_expected))
            ratio_errors.append(abs(ratio - self.that_sqrtPi_limit_expected))
            print(
                f"Pi={fmt(pi_value)}: 1-g={fmt(1.0 - g_value)}, "
                f"denom={fmt(denom)}, T_hat/sqrt(Pi)={fmt(ratio)}"
            )
            require(
                f"Pi={fmt(pi_value)} keeps g(Pi)<1",
                g_value < 1.0,
                f"g(Pi)={fmt(g_value)}",
            )
            require(
                f"Pi={fmt(pi_value)} asymptotic T_hat ratio stays near limit",
                abs(ratio - self.that_sqrtPi_limit_expected) <= self.ratio_tol,
                f"ratio={fmt(ratio)}, limit={fmt(self.that_sqrtPi_limit_expected)}",
            )

        require(
            "1-g(Pi) decreases toward the singular limit",
            all(b < a for a, b in zip([1.0 - g for g in g_tail[:-1]], [1.0 - g for g in g_tail[1:]], strict=True)),
            f"tail gaps={[fmt(1.0 - g) for g in g_tail]}",
        )
        require(
            "denominator approaches positive large-Pi limit",
            all(b < a for a, b in zip(denom_errors[:-1], denom_errors[1:], strict=True)),
            f"errors={[fmt(err) for err in denom_errors]}",
        )
        require(
            "T_hat/sqrt(Pi) approaches the predicted asymptotic coefficient",
            all(b < a for a, b in zip(ratio_errors[:-1], ratio_errors[1:], strict=True)),
            f"errors={[fmt(err) for err in ratio_errors]}",
        )
        require(
            "large-Pi denominator limit matches Stage 126",
            near(1.0 - self.R_inf_expected, self.denom_inf_expected, self.root_tol),
            f"1-R_inf={fmt(1.0 - self.R_inf_expected)}, expected={fmt(self.denom_inf_expected)}",
        )


def main(argv: list[str]) -> int:
    root = Path(__file__).resolve().parents[3]
    config_path = Path(argv[1]) if len(argv) > 1 else root / "scripts/moving_throat/numerical/stage125_127_mouth_branch_samples.json"
    config = load_config(config_path)
    h = MouthBranchHarness(config)
    print(f"Loaded config from {config_path}")
    h.stage125_samples()
    h.stage127_branch_selection()
    h.stage126_asymptotics()
    print("\nAll stage-125/127 mouth-branch stress checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
