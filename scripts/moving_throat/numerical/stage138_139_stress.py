#!/usr/bin/env python3
"""Numerical stress harness for the Family-1 co-evolving mouth/core map (Stages 138-139)."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path


def near(lhs: float, rhs: float, tol: float) -> bool:
    return abs(lhs - rhs) <= tol * (1.0 + abs(rhs))


def require(label: str, condition: bool, detail: str) -> None:
    status = "PASS" if condition else "FAIL"
    print(f"[{status}] {label}: {detail}")
    if not condition:
        raise AssertionError(label)


def fmt(value: float) -> str:
    return f"{value:.12g}"


def trapezoid_weights(n: int) -> list[float]:
    if n < 2:
        raise ValueError("grid needs at least 2 points")
    h = 1.0 / (n - 1)
    w = [h] * n
    w[0] = h / 2.0
    w[-1] = h / 2.0
    return w


def load_config(path: Path) -> dict:
    data = json.loads(path.read_text())
    if data.get("schema") != "moving_throat_numerical_stage138_139_v1":
        raise ValueError("unexpected config schema")
    return data


def grid_points(n: int) -> list[float]:
    return [i / (n - 1) for i in range(n)]


def exp_profile(pi: float, xs: list[float], weights: list[float]) -> list[float]:
    vals = [math.exp(-pi * x) for x in xs]
    z = sum(v * w for v, w in zip(vals, weights))
    return [v / z for v in vals]


def uniform_profile(xs: list[float], weights: list[float]) -> list[float]:
    vals = [1.0 for _ in xs]
    z = sum(v * w for v, w in zip(vals, weights))
    return [v / z for v in vals]


def endpoint_profile(xs: list[float], weights: list[float], mix: float) -> list[float]:
    # Positive, normalized source with extra weight near x=0.
    vals = [math.exp(-mix * 6.0 * x) + 0.15 for x in xs]
    z = sum(v * w for v, w in zip(vals, weights))
    return [v / z for v in vals]


def derivative_like_profile(xs: list[float], weights: list[float]) -> list[float]:
    vals = [math.cos(math.pi * x / 2.0) for x in xs]
    z = sum(v * w for v, w in zip(vals, weights))
    return [v / z for v in vals]


def seed_profile(seed: dict, xs: list[float], weights: list[float]) -> list[float]:
    name = seed["name"]
    if name == "canonical_exponential":
        return exp_profile(float(seed["Pi"]), xs, weights)
    if name == "uniform_broad":
        return uniform_profile(xs, weights)
    if name == "derivative_like":
        return derivative_like_profile(xs, weights)
    if name == "endpoint_biased":
        return endpoint_profile(xs, weights, float(seed["mix"]))
    raise ValueError(f"unknown seed {name}")


def moments(profile: list[float], xs: list[float], weights: list[float]) -> dict[str, float]:
    c = [math.cos(math.pi * x / 2.0) for x in xs]
    s = [math.cosh(math.pi * (1.0 - x) / 2.0) / math.cosh(math.pi / 2.0) for x in xs]
    g = sum(p * ci * w for p, ci, w in zip(profile, c, weights))
    S = sum(p * si * w for p, si, w in zip(profile, s, weights))
    return {"g": g, "S": S}


def kernels(xs: list[float]) -> tuple[list[list[float]], list[list[float]]]:
    n = len(xs)
    Kt = [[0.0] * n for _ in range(n)]
    Kq = [[0.0] * n for _ in range(n)]
    norm = (math.pi / 2.0) * math.cosh(math.pi / 2.0)
    for i, x in enumerate(xs):
        for j, y in enumerate(xs):
            mn = x if x < y else y
            mx = x if x > y else y
            Kt[i][j] = mn
            Kq[i][j] = (
                math.sinh(math.pi * mn / 2.0)
                * math.cosh(math.pi * (1.0 - mx) / 2.0)
                / norm
            )
    return Kt, Kq


def apply_kernel(kernel: list[list[float]], profile: list[float], weights: list[float]) -> list[float]:
    return [sum(kij * pj * w for kij, pj, w in zip(row, profile, weights)) for row in kernel]


def update_map(
    profile: list[float],
    sigma0: float,
    r_f1: float,
    xs: list[float],
    weights: list[float],
    Kt: list[list[float]],
    Kq: list[list[float]],
) -> dict:
    mom = moments(profile, xs, weights)
    g = mom["g"]
    S = mom["S"]
    R = (g - r_f1) ** 2 / (1.0 + r_f1 ** 2)
    Ts = apply_kernel(Kt, profile, weights)
    Tq = apply_kernel(Kq, profile, weights)
    phi = [sigma0 * (ts - R * tq) for ts, tq in zip(Ts, Tq)]
    vals = [math.exp(-p) for p in phi]
    z = sum(v * w for v, w in zip(vals, weights))
    new_profile = [v / z for v in vals]
    new_mom = moments(new_profile, xs, weights)
    new_R = (new_mom["g"] - r_f1) ** 2 / (1.0 + r_f1 ** 2)
    return {
        "profile": new_profile,
        "g": new_mom["g"],
        "S": new_mom["S"],
        "R": new_R,
        "phi0": phi[0],
        "phi1": phi[-1],
    }


def fixed_point(
    seed: list[float],
    sigma0: float,
    r_f1: float,
    xs: list[float],
    weights: list[float],
    Kt: list[list[float]],
    Kq: list[list[float]],
    tol: float,
    max_iter: int,
    relax: float,
) -> dict:
    prof = seed[:]
    last = None
    for it in range(1, max_iter + 1):
        upd = update_map(prof, sigma0, r_f1, xs, weights, Kt, Kq)
        prof_next = [(1.0 - relax) * p + relax * q for p, q in zip(prof, upd["profile"])]
        z = sum(p * w for p, w in zip(prof_next, weights))
        prof_next = [p / z for p in prof_next]
        mom = moments(prof_next, xs, weights)
        R = (mom["g"] - r_f1) ** 2 / (1.0 + r_f1 ** 2)
        change = max(abs(a - b) for a, b in zip(prof_next, prof))
        if last is not None:
            delta_mom = max(abs(mom["g"] - last["g"]), abs(mom["S"] - last["S"]), abs(R - last["R"]))
        else:
            delta_mom = float("inf")
        prof = prof_next
        last = {"g": mom["g"], "S": mom["S"], "R": R}
        if change < tol and delta_mom < tol:
            return {"profile": prof, "iterations": it, **last}
    raise RuntimeError("fixed point did not converge")


def bracket_root(samples: list[tuple[float, float]], target: float) -> tuple[float, float]:
    low = high = None
    for sigma0, gfp in samples:
        if gfp < target:
            low = (sigma0, gfp)
        if gfp > target and low is not None:
            high = (sigma0, gfp)
            break
    if low is None or high is None:
        raise RuntimeError("could not bracket target root from samples")
    return low[0], high[0]


def solve_root(
    sigma_lo: float,
    sigma_hi: float,
    target: float,
    seed_profiles: list[list[float]],
    r_f1: float,
    xs: list[float],
    weights: list[float],
    Kt: list[list[float]],
    Kq: list[list[float]],
    tol: float,
    max_iter: int,
    relax: float,
) -> dict:
    lo = sigma_lo
    hi = sigma_hi
    flow = None
    fhi = None
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        fp_vals = []
        for seed in seed_profiles:
            fp = fixed_point(seed, mid, r_f1, xs, weights, Kt, Kq, tol, max_iter, relax)
            fp_vals.append(fp["g"])
        gmid = sum(fp_vals) / len(fp_vals)
        if flow is None:
            flow = sum(
                fixed_point(seed, lo, r_f1, xs, weights, Kt, Kq, tol, max_iter, relax)["g"]
                for seed in seed_profiles
            ) / len(seed_profiles)
            fhi = sum(
                fixed_point(seed, hi, r_f1, xs, weights, Kt, Kq, tol, max_iter, relax)["g"]
                for seed in seed_profiles
            ) / len(seed_profiles)
        if (flow - target) * (gmid - target) <= 0:
            hi = mid
            fhi = gmid
        else:
            lo = mid
            flow = gmid
        if abs(hi - lo) < tol:
            break
    return {"sigma0": 0.5 * (lo + hi), "g_mid": gmid, "low": lo, "high": hi}


def main(argv: list[str]) -> int:
    root = Path(__file__).resolve().parents[3]
    config_path = Path(argv[1]) if len(argv) > 1 else root / "scripts/moving_throat/numerical/stage138_139_config.json"
    cfg = load_config(config_path)
    xs = grid_points(int(cfg["grid"]["points"]))
    weights = trapezoid_weights(len(xs))
    Kt, Kq = kernels(xs)
    const = cfg["constants"]
    r_f1 = float(const["r_F1"])
    g_star = float(const["g_star"])
    sigma0_star = float(cfg["frozen_traction"]["Sigma0"])
    stage138_ref = cfg["stage138_reference"]
    stage139_ref = cfg["stage139_reference"]
    tol = float(cfg["tolerances"]["fixed_point_tol"])
    mom_tol = float(cfg["tolerances"]["moment_tol"])
    max_iter = int(cfg["tolerances"]["max_iter"])
    relax = float(cfg["tolerances"]["relaxation"])

    print(f"Loaded config from {config_path}")
    print(f"Grid points: {len(xs)}")

    frozen_seeds = cfg["frozen_traction"]["seeds"]
    frozen_results = []
    print("\n=== Frozen traction multi-seed solve ===")
    for seed in frozen_seeds:
        profile = seed_profile(seed, xs, weights)
        fp = fixed_point(profile, sigma0_star, r_f1, xs, weights, Kt, Kq, tol, max_iter, relax)
        frozen_results.append(fp)
        print(
            f"[PASS] seed {seed['name']}: iterations={fp['iterations']}, "
            f"g_fp={fmt(fp['g'])}, S_fp={fmt(fp['S'])}, R_fp={fmt(fp['R'])}"
        )

    g_vals = [fp["g"] for fp in frozen_results]
    S_vals = [fp["S"] for fp in frozen_results]
    R_vals = [fp["R"] for fp in frozen_results]
    g_avg = sum(g_vals) / len(g_vals)
    S_avg = sum(S_vals) / len(S_vals)
    R_avg = sum(R_vals) / len(R_vals)
    for val in g_vals:
        require("frozen-traction g consistency", near(val, g_avg, mom_tol), f"g_fp={fmt(val)}, mean={fmt(g_avg)}")
    for val in S_vals:
        require("frozen-traction S consistency", near(val, S_avg, mom_tol), f"S_fp={fmt(val)}, mean={fmt(S_avg)}")
    for val in R_vals:
        require("frozen-traction R consistency", near(val, R_avg, mom_tol), f"R_fp={fmt(val)}, mean={fmt(R_avg)}")
    require(
        "frozen traction stays below compensated target",
        g_avg < g_star,
        f"g_fp(mean)={fmt(g_avg)}, target={fmt(g_star)}",
    )
    require(
        "frozen traction remains above compensated R threshold",
        R_avg > 0.25,
        f"R_fp(mean)={fmt(R_avg)}",
    )
    require(
        "frozen traction matches Stage 138 reference g",
        near(g_avg, float(stage138_ref["g_fp"]), 1e-5),
        f"g_fp(mean)={fmt(g_avg)}, ref={fmt(float(stage138_ref['g_fp']))}",
    )
    require(
        "frozen traction matches Stage 138 reference S",
        near(S_avg, float(stage138_ref["S_fp"]), 1e-5),
        f"S_fp(mean)={fmt(S_avg)}, ref={fmt(float(stage138_ref['S_fp']))}",
    )
    require(
        "frozen traction matches Stage 138 reference R",
        near(R_avg, float(stage138_ref["R_fp"]), 1e-5),
        f"R_fp(mean)={fmt(R_avg)}, ref={fmt(float(stage138_ref['R_fp']))}",
    )

    print(f"frozen-trction mean g = {fmt(g_avg)}")
    print(f"frozen-traction mean S = {fmt(S_avg)}")
    print(f"frozen-traction mean R = {fmt(R_avg)}")

    print("\n=== Root search samples ===")
    root_seeds = cfg["root_search"]["seeds"]
    seed_profiles = [seed_profile(seed, xs, weights) for seed in root_seeds]
    sample_values = []
    for sigma0 in cfg["root_search"]["samples"]:
        per_seed = []
        for seed, profile in zip(root_seeds, seed_profiles):
            fp = fixed_point(profile, float(sigma0), r_f1, xs, weights, Kt, Kq, tol, max_iter, relax)
            per_seed.append(fp["g"])
        g_mean = sum(per_seed) / len(per_seed)
        sample_values.append((float(sigma0), g_mean))
        print(f"sigma0={fmt(float(sigma0))} -> g_fp(mean)={fmt(g_mean)}")
        for val in per_seed:
            require(
                f"seed agreement at Sigma0={fmt(float(sigma0))}",
                near(val, g_mean, mom_tol),
                f"seed g_fp={fmt(val)}, mean={fmt(g_mean)}",
            )
    for (_, left_val), (_, right_val) in zip(sample_values, sample_values[1:]):
        require(
            "sample monotonicity on analyzed interval",
            right_val > left_val,
            f"g_left={fmt(left_val)}, g_right={fmt(right_val)}",
        )
    lo, hi = bracket_root(sample_values, g_star)
    require("root bracket order", lo < hi, f"lo={fmt(lo)}, hi={fmt(hi)}")
    print(f"bracket from samples: [{fmt(lo)}, {fmt(hi)}]")

    # Solve the root by bisection on the averaged fixed-point moment.
    target = float(cfg["root_search"]["target_g"])
    left = lo
    right = hi
    left_val = next(val for sigma0, val in sample_values if sigma0 == left)
    right_val = next(val for sigma0, val in sample_values if sigma0 == right)
    for _ in range(60):
        mid = 0.5 * (left + right)
        mid_vals = [
            fixed_point(profile, mid, r_f1, xs, weights, Kt, Kq, tol, max_iter, relax)["g"]
            for profile in seed_profiles
        ]
        mid_val = sum(mid_vals) / len(mid_vals)
        if (left_val - target) * (mid_val - target) <= 0:
            right = mid
            right_val = mid_val
        else:
            left = mid
            left_val = mid_val
        if abs(right - left) < float(cfg["tolerances"]["root_tol"]):
            break

    sigma0_can = 0.5 * (left + right)
    root_fp = [
        fixed_point(profile, sigma0_can, r_f1, xs, weights, Kt, Kq, tol, max_iter, relax)
        for profile in seed_profiles
    ]
    root_g = sum(fp["g"] for fp in root_fp) / len(root_fp)
    root_S = sum(fp["S"] for fp in root_fp) / len(root_fp)
    root_R = sum(fp["R"] for fp in root_fp) / len(root_fp)
    for fp in root_fp:
        require("root-seed consistency", near(fp["g"], root_g, mom_tol), f"g_fp={fmt(fp['g'])}, mean={fmt(root_g)}")

    print("\n=== Renormalized canonical solve ===")
    print(f"sigma0_can = {fmt(sigma0_can)}")
    print(f"g_can(mean) = {fmt(root_g)}")
    print(f"S_can(mean) = {fmt(root_S)}")
    print(f"R_can(mean) = {fmt(root_R)}")
    require("root moment close to target", near(root_g, target, 1e-7), f"g_can={fmt(root_g)}, target={fmt(target)}")
    require("compensated branch recovered", near(root_R, 0.25, 1e-7), f"R_can={fmt(root_R)}")
    require(
        "root traction matches Stage 139 reference",
        near(sigma0_can, float(stage139_ref["Sigma0_can"]), 1e-4),
        f"sigma0_can={fmt(sigma0_can)}, ref={fmt(float(stage139_ref['Sigma0_can']))}",
    )
    require(
        "root mouth traction matches Stage 139 reference",
        near(math.sqrt(9.0 * sigma0_can / 20.0), float(stage139_ref["Tm_can"]), 1e-4),
        f"Tm_can={fmt(math.sqrt(9.0 * sigma0_can / 20.0))}, ref={fmt(float(stage139_ref['Tm_can']))}",
    )

    tm_can = math.sqrt(9.0 * sigma0_can / 20.0)
    print(f"Tm_can = {fmt(tm_can)}")
    print(f"Pi_can estimated = {fmt(sigma0_can * (1.0 - root_R * root_S))}")

    print("\nAll Family-1 co-evolution stress checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
