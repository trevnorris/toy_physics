#!/usr/bin/env python3
"""Multi-seed numerical stress harness for the co-evolving fixed-point chain."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np


def near(lhs: float, rhs: float, tol: float) -> bool:
    return abs(lhs - rhs) <= tol * (1.0 + abs(rhs))


def require(label: str, condition: bool, detail: str) -> None:
    status = "PASS" if condition else "FAIL"
    print(f"[{status}] {label}: {detail}")
    if not condition:
        raise AssertionError(label)


def fmt(value: float) -> str:
    return f"{value:.12g}"


def load_config(path: Path) -> dict:
    data = json.loads(path.read_text())
    if data.get("schema") != "moving_throat_numerical_stage138_139_v1":
        raise ValueError("unexpected sample schema")
    return data


class FixedPointHarness:
    def __init__(self, config: dict) -> None:
        const = config["constants"]
        self.rF1 = float(const["rF1"])
        self.g_star = float(const["g_star"])
        self.Pi_star = float(const["Pi_star"])
        self.Sigma0_star = float(const["Sigma0_star"])
        self.T_hat_star = float(const["T_hat_star"])
        self.Sigma0_can_expected = float(const["Sigma0_can_expected"])
        self.S_can_expected = float(const["S_can_expected"])
        self.Pi_can_expected = float(const["Pi_can_expected"])
        self.T_hat_can_expected = float(const["T_hat_can_expected"])
        self.grid_points = int(config["grid_points"])
        self.max_iterations = int(config["max_iterations"])
        self.fixed_point_tol = float(config["fixed_point_tol"])
        self.profile_tol = float(config["profile_tolerance"])
        self.moment_tol = float(config["moment_tolerance"])
        self.root_tol = float(config["root_tolerance"])
        self.bisection_iterations = int(config["bisection_iterations"])
        self.seeds = config["seeds"]

        self.x = np.linspace(0.0, 1.0, self.grid_points)
        dx = self.x[1] - self.x[0]
        self.w = np.full(self.grid_points, dx)
        self.w[0] = dx / 2.0
        self.w[-1] = dx / 2.0
        self.kappa = math.pi / 2.0
        self.c_grid = np.cos(math.pi * self.x / 2.0)
        self.Kq_grid = np.cosh(self.kappa * (1.0 - self.x)) / math.cosh(self.kappa)
        self._canonical_seed = self.normalize(self.Pi_star * np.exp(-self.Pi_star * self.x))
        self._g_cache: dict[float, float] = {}

    def normalize(self, sig: np.ndarray) -> np.ndarray:
        return sig / np.sum(sig * self.w)

    def seed_profile(self, seed: dict) -> np.ndarray:
        kind = seed["kind"]
        if kind == "canonical_exponential":
            raw = self.Pi_star * np.exp(-self.Pi_star * self.x)
        elif kind == "uniform":
            raw = np.ones_like(self.x)
        elif kind == "derivative_match":
            raw = (math.pi / 2.0) * np.cos(math.pi * self.x / 2.0)
        elif kind == "convex_blend":
            lam = float(seed["lambda"])
            raw = (1.0 - lam) * np.ones_like(self.x) + lam * (math.pi / 2.0) * np.cos(math.pi * self.x / 2.0)
        else:
            raise ValueError(f"unknown seed kind: {kind}")
        return self.normalize(raw)

    def Ts(self, sig: np.ndarray) -> np.ndarray:
        cum_sig = np.cumsum(sig * self.w)
        cum_y = np.cumsum(self.x * sig * self.w)
        return cum_y + self.x * (cum_sig[-1] - cum_sig)

    def Tq(self, sig: np.ndarray) -> np.ndarray:
        a_term = np.cumsum(np.sinh(self.kappa * self.x) * sig * self.w)
        b_term = np.cumsum((np.cosh(self.kappa * (1.0 - self.x)) * sig * self.w)[::-1])[::-1]
        return (
            np.cosh(self.kappa * (1.0 - self.x)) * a_term
            + np.sinh(self.kappa * self.x) * b_term
        ) / (self.kappa * math.cosh(self.kappa))

    def g(self, sig: np.ndarray) -> float:
        return float(np.sum(self.c_grid * sig * self.w))

    def S(self, sig: np.ndarray) -> float:
        return float(np.sum(self.Kq_grid * sig * self.w))

    def R(self, sig: np.ndarray) -> float:
        gv = self.g(sig)
        return ((gv - self.rF1) ** 2) / (1.0 + self.rF1 ** 2)

    def Phi(self, sig: np.ndarray, sigma0: float) -> np.ndarray:
        return sigma0 * (self.Ts(sig) - self.R(sig) * self.Tq(sig))

    def next_sigma(self, sig: np.ndarray, sigma0: float) -> np.ndarray:
        ph = self.Phi(sig, sigma0)
        ph_shift = ph - np.min(ph)
        return self.normalize(np.exp(-ph_shift))

    def solve_fixed_point(self, sigma0: float, initial: np.ndarray | None = None) -> tuple[np.ndarray, int, float]:
        sig = self._canonical_seed.copy() if initial is None else initial.copy()
        err = math.inf
        for n in range(self.max_iterations):
            sig_new = self.next_sigma(sig, sigma0)
            err = float(np.max(np.abs(sig_new - sig)))
            sig = sig_new
            if err < self.fixed_point_tol:
                return sig, n + 1, err
        raise RuntimeError(f"fixed-point solve did not converge for sigma0={sigma0}")

    def fixed_point_summary(self, sigma0: float, initial: np.ndarray) -> dict:
        sig, niter, err = self.solve_fixed_point(sigma0, initial)
        g_fp = self.g(sig)
        S_fp = self.S(sig)
        R_fp = self.R(sig)
        Pi_fp = sigma0 * (1.0 - R_fp * S_fp)
        T_hat = math.sqrt(9.0 * sigma0 / 20.0)
        return {
            "sigma": sig,
            "iterations": niter,
            "error": err,
            "g": g_fp,
            "S": S_fp,
            "R": R_fp,
            "Pi": Pi_fp,
            "T_hat": T_hat,
        }

    def g_fp_default(self, sigma0: float) -> float:
        key = round(float(sigma0), 14)
        if key not in self._g_cache:
            summary = self.fixed_point_summary(float(sigma0), self._canonical_seed)
            self._g_cache[key] = summary["g"]
        return self._g_cache[key]

    def bracket_root(self, lo: float, hi: float) -> tuple[float, float]:
        flo = self.g_fp_default(lo) - self.g_star
        fhi = self.g_fp_default(hi) - self.g_star
        require(
            "root bracket sign change",
            flo == 0.0 or fhi == 0.0 or flo * fhi < 0.0,
            f"[{fmt(lo)}, {fmt(hi)}] -> flo={fmt(flo)}, fhi={fmt(fhi)}",
        )
        return flo, fhi

    def bisect_root(self, lo: float, hi: float) -> float:
        flo, fhi = self.bracket_root(lo, hi)
        left = lo
        right = hi
        for _ in range(self.bisection_iterations):
            mid = 0.5 * (left + right)
            fm = self.g_fp_default(mid) - self.g_star
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


def check_seed_consistency(h: FixedPointHarness, case: dict) -> None:
    sigma0 = float(case["sigma0"])
    print(f"\n=== {case['name']} ({case['kind']}) ===")
    print(f"sigma0 = {fmt(sigma0)}")
    if "assumptions" in case:
        print("assumptions:")
        for item in case["assumptions"]:
            print(f"  - {item}")

    summaries = []
    for seed in h.seeds:
        initial = h.seed_profile(seed)
        summary = h.fixed_point_summary(sigma0, initial)
        summaries.append((seed["name"], summary))
        print(
            f"seed {seed['name']}: iterations={summary['iterations']}, "
            f"max_residual={fmt(summary['error'])}, g={fmt(summary['g'])}, "
            f"S={fmt(summary['S'])}, R={fmt(summary['R'])}, Pi={fmt(summary['Pi'])}"
        )

    ref_name, ref = summaries[0]
    for name, summary in summaries[1:]:
        require(
            f"profile agreement vs {ref_name} for {name}",
            float(np.max(np.abs(summary["sigma"] - ref["sigma"]))) <= h.profile_tol,
            f"L_inf={fmt(float(np.max(np.abs(summary['sigma'] - ref['sigma']))))}",
        )
        require(
            f"g agreement vs {ref_name} for {name}",
            near(summary["g"], ref["g"], h.moment_tol),
            f"g={fmt(summary['g'])}, ref={fmt(ref['g'])}",
        )
        require(
            f"S agreement vs {ref_name} for {name}",
            near(summary["S"], ref["S"], h.moment_tol),
            f"S={fmt(summary['S'])}, ref={fmt(ref['S'])}",
        )
        require(
            f"R agreement vs {ref_name} for {name}",
            near(summary["R"], ref["R"], h.moment_tol),
            f"R={fmt(summary['R'])}, ref={fmt(ref['R'])}",
        )
        require(
            f"Pi agreement vs {ref_name} for {name}",
            near(summary["Pi"], ref["Pi"], h.moment_tol),
            f"Pi={fmt(summary['Pi'])}, ref={fmt(ref['Pi'])}",
        )

    if near(sigma0, h.Sigma0_star, h.root_tol):
        dg = ref["g"] - h.g_star
        pred = -dg / math.sqrt(1.0 + h.rF1 ** 2) + dg * dg / (1.0 + h.rF1 ** 2)
        direct = ref["R"] - 0.25
        require(
            "Stage-138 transport-law consistency",
            near(pred, direct, h.moment_tol),
            f"pred={fmt(pred)}, direct={fmt(direct)}",
        )


def check_stage139(h: FixedPointHarness, config: dict) -> None:
    print("\n=== stage139_scan_and_root ===")
    scan_points = [float(x) for x in config["stage139_scan_points"]]
    scan_values = [h.g_fp_default(sigma0) for sigma0 in scan_points]
    for sigma0, value in zip(scan_points, scan_values, strict=True):
        print(f"scan sigma0={fmt(sigma0)} -> g_fp={fmt(value)}")
    diffs = [b - a for a, b in zip(scan_values[:-1], scan_values[1:], strict=True)]
    require(
        "monotone increase on analyzed scan grid",
        all(diff > 0.0 for diff in diffs),
        f"min diff={fmt(min(diffs))}",
    )

    roots = []
    for lo, hi in config["stage139_root_brackets"]:
        lo = float(lo)
        hi = float(hi)
        root = h.bisect_root(lo, hi)
        roots.append(root)
        print(f"root from bracket [{fmt(lo)}, {fmt(hi)}] = {fmt(root)}")

    ref_root = roots[0]
    for idx, root in enumerate(roots[1:], start=1):
        require(
            f"root consistency vs bracket 0 for bracket {idx}",
            near(root, ref_root, h.root_tol),
            f"root={fmt(root)}, ref={fmt(ref_root)}",
        )

    final = h.fixed_point_summary(ref_root, h._canonical_seed)
    require("g(root) = g_star", near(final["g"], h.g_star, h.moment_tol), f"g={fmt(final['g'])}, g_star={fmt(h.g_star)}")
    require("R(root) = 1/4", near(final["R"], 0.25, h.moment_tol), f"R={fmt(final['R'])}")
    require(
        "Sigma0_can matches expected audit value",
        near(ref_root, h.Sigma0_can_expected, h.root_tol),
        f"root={fmt(ref_root)}, expected={fmt(h.Sigma0_can_expected)}",
    )
    require(
        "S_can matches expected audit value",
        near(final["S"], h.S_can_expected, h.moment_tol),
        f"S={fmt(final['S'])}, expected={fmt(h.S_can_expected)}",
    )
    require(
        "Pi_can matches expected audit value",
        near(final["Pi"], h.Pi_can_expected, h.moment_tol),
        f"Pi={fmt(final['Pi'])}, expected={fmt(h.Pi_can_expected)}",
    )
    require(
        "T_hat_can matches expected audit value",
        near(final["T_hat"], h.T_hat_can_expected, h.moment_tol),
        f"T_hat={fmt(final['T_hat'])}, expected={fmt(h.T_hat_can_expected)}",
    )

    print(f"Sigma0_can = {fmt(ref_root)}")
    print(f"g_can      = {fmt(final['g'])}")
    print(f"S_can      = {fmt(final['S'])}")
    print(f"R_can      = {fmt(final['R'])}")
    print(f"Pi_can     = {fmt(final['Pi'])}")
    print(f"T_hat_can  = {fmt(final['T_hat'])}")

    print("\nRelative shifts from original canonical point:")
    print(f"Sigma0 ratio - 1 = {fmt(ref_root / h.Sigma0_star - 1.0)}")
    print(f"Pi ratio - 1     = {fmt(final['Pi'] / h.Pi_star - 1.0)}")
    print(f"T_hat ratio - 1  = {fmt(final['T_hat'] / h.T_hat_star - 1.0)}")

    print("\n=== stage139_negative_controls ===")
    for control in config["stage139_negative_controls"]:
        sigma0 = float(control["sigma0"])
        summary = h.fixed_point_summary(sigma0, h._canonical_seed)
        side = control["side"]
        print(
            f"{control['name']}: sigma0={fmt(sigma0)}, g={fmt(summary['g'])}, "
            f"R={fmt(summary['R'])}, Pi={fmt(summary['Pi'])}"
        )
        require(
            f"{control['name']} stays away from the compensated moment",
            abs(summary["g"] - h.g_star) >= float(control["min_g_gap"]),
            f"|g-g_star|={fmt(abs(summary['g'] - h.g_star))}, min gap={fmt(float(control['min_g_gap']))}",
        )
        require(
            f"{control['name']} stays away from R=1/4",
            abs(summary["R"] - 0.25) >= float(control["min_R_gap"]),
            f"|R-1/4|={fmt(abs(summary['R'] - 0.25))}, min gap={fmt(float(control['min_R_gap']))}",
        )
        if side == "below":
            require(
                f"{control['name']} remains on the pre-compensation side",
                summary["g"] < h.g_star and summary["R"] > 0.25,
                f"g={fmt(summary['g'])}, g_star={fmt(h.g_star)}, R={fmt(summary['R'])}",
            )
        elif side == "above":
            require(
                f"{control['name']} remains on the post-compensation side",
                summary["g"] > h.g_star and summary["R"] < 0.25,
                f"g={fmt(summary['g'])}, g_star={fmt(h.g_star)}, R={fmt(summary['R'])}",
            )
        else:
            raise ValueError(f"unknown negative-control side: {side}")


def main(argv: list[str]) -> int:
    root = Path(__file__).resolve().parents[2]
    config_path = Path(argv[1]) if len(argv) > 1 else root / "scripts/numerical/stage138_139_fixedpoint_samples.json"
    config = load_config(config_path)
    h = FixedPointHarness(config)
    print(f"Loaded config from {config_path}")
    for case in config["seed_check_sigmas"]:
        check_seed_consistency(h, case)
    check_stage139(h, config)
    print("\nAll stage-138/139 fixed-point stress checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
