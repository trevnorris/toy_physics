#!/usr/bin/env python3
# locality_scan.py (patched v3.2)

"""
python locality_scan.py \
  --mode match \
  --eta_ref 0.012 \
  --eta_bounds 0.002,0.25 \
  --rtol 0.03 --max_iter 14 \
  --sigmas 1.2 \
  --yg_list 0,1,2,3,4,5,6,7,8,9,10,11 \
  --yg_ref 0 \
  --seeds 1,2,3 \
  --jitters 0.02,0.03 \
  --tmax 12 --dt 0.005 \
  --csv_runs wp2_loc_match_dense_runs.csv \
  --csv_agg  wp2_loc_match_dense_agg.csv \
  --fig_out  wp2_loc_match_dense.png
"""

from __future__ import annotations
import argparse, csv, math
from dataclasses import dataclass, asdict
from typing import Any, Dict, List, Tuple

import numpy as np
import matplotlib.pyplot as plt

import wp2_onset_scan as base


def _parse_list(s: str, cast=float) -> List:
    return [cast(x) for x in s.split(",") if x.strip() != ""]


def _safe_div(a: float, b: float) -> float:
    if b == 0.0 or (isinstance(b, float) and math.isnan(b)):
        return float("nan")
    return a / b


def _trapz_safe(y, x) -> float:
    try:
        if y is None or x is None:
            return float("nan")
        y = np.asarray(y, dtype=float)
        x = np.asarray(x, dtype=float)
        if y.size < 2 or x.size < 2:
            return float("nan")
        if not np.all(np.isfinite(x)):
            return float("nan")
        return float(np.trapz(y, x))
    except Exception:
        return float("nan")


def _meanstd_finite(vals: np.ndarray) -> Tuple[float, float]:
    vals = np.asarray(vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return float("nan"), float("nan")
    return float(vals.mean()), float(vals.std(ddof=0))


def _fmt_onset(t_onset: float, Tmax: float) -> str:
    if hasattr(base, "fmt_onset"):
        return base.fmt_onset(t_onset, Tmax)
    if isinstance(t_onset, float) and math.isnan(t_onset):
        return f">{Tmax:.2f}"
    return f"{t_onset:.4f}"


@dataclass
class Row:
    mode: str        # fixed | match
    kind: str        # OFF | RUN | REF
    case: str
    eta_peak: float
    yg: float
    sigma: float
    seed: int
    jitter: float
    t_onset_curl: float
    t_onset_div: float
    jc_ratio: float
    jd_ratio: float
    Deta_int: float
    Deta_target: float
    Deta_err: float
    eff_curl: float
    eff_div: float
    Tmax: float


def _run_case(
    *,
    name: str,
    P: Dict[str, Any],
    eta_peak: float,
    yg: float,
    sigma: float,
    onset_cfg: Dict[str, Any],
    off_jc_late: float,
    off_jd_late: float,
    mode: str,
    kind: str,
    Deta_target: float = float("nan"),
) -> Row:
    out = base.run_sim_timeseries(name, **{**P, "eta_peak": float(eta_peak), "yg": float(yg), "sigma": float(sigma)})
    s = base.summarize(out, **onset_cfg)

    t_arr = out.get("t", None)
    Deta_t = out.get("Deta", None)
    Deta_int = _trapz_safe(Deta_t, t_arr)

    if kind == "OFF":
        jc_ratio = 1.0
        jd_ratio = 1.0
    else:
        jc_ratio = _safe_div(s["jc_late_mean"], off_jc_late)
        jd_ratio = _safe_div(s["jd_late_mean"], off_jd_late)

    eff_c = float("nan")
    eff_d = float("nan")
    if np.isfinite(Deta_int) and Deta_int > 0:
        if np.isfinite(jc_ratio):
            eff_c = (1.0 - jc_ratio) / Deta_int
        if np.isfinite(jd_ratio):
            eff_d = (1.0 - jd_ratio) / Deta_int

    Derr = float("nan")
    if np.isfinite(Deta_target) and np.isfinite(Deta_int):
        Derr = Deta_int - Deta_target

    return Row(
        mode=mode, kind=kind, case=name,
        eta_peak=float(eta_peak), yg=float(yg), sigma=float(sigma),
        seed=int(P["seed"]), jitter=float(P["jitter"]),
        t_onset_curl=float(s["tcurl"]), t_onset_div=float(s["tdiv"]),
        jc_ratio=float(jc_ratio), jd_ratio=float(jd_ratio),
        Deta_int=float(Deta_int),
        Deta_target=float(Deta_target),
        Deta_err=float(Derr),
        eff_curl=float(eff_c), eff_div=float(eff_d),
        Tmax=float(s["Tmax"]),
    )


def _match_eta_for_work(
    *,
    P: Dict[str, Any],
    yg: float,
    sigma: float,
    onset_cfg: Dict[str, Any],
    off_jc_late: float,
    off_jd_late: float,
    target: float,
    eta_min: float,
    eta_max: float,
    rtol: float,
    atol: float,
    max_iter: int,
    cache: Dict[Tuple[int, float, float, float, float], Row],
    name_prefix: str,
) -> Row:

    def eval_eta(ep: float) -> Row:
        key = (int(P["seed"]), float(P["jitter"]), float(yg), float(sigma), float(ep))
        if key in cache:
            return cache[key]
        row = _run_case(
            name=f"{name_prefix}_ep{ep:.6g}",
            P=P, eta_peak=ep, yg=yg, sigma=sigma,
            onset_cfg=onset_cfg,
            off_jc_late=off_jc_late, off_jd_late=off_jd_late,
            mode="match", kind="RUN",
            Deta_target=target,
        )
        cache[key] = row
        return row

    # evaluate bounds
    r_lo = eval_eta(eta_min)
    r_hi = eval_eta(eta_max)

    if not np.isfinite(r_lo.Deta_int) or not np.isfinite(r_hi.Deta_int):
        cand = [r for r in [r_lo, r_hi] if np.isfinite(r.Deta_int)]
        return cand[0] if cand else r_hi

    # if not bracketed, return closer bound
    if (r_lo.Deta_int - target) * (r_hi.Deta_int - target) > 0:
        return r_lo if abs(r_lo.Deta_int - target) < abs(r_hi.Deta_int - target) else r_hi

    lo, hi = eta_min, eta_max
    r_best = r_lo if abs(r_lo.Deta_int - target) < abs(r_hi.Deta_int - target) else r_hi
    tol = max(atol, rtol * target)

    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        r_mid = eval_eta(mid)
        if not np.isfinite(r_mid.Deta_int):
            hi = mid
            continue

        if abs(r_mid.Deta_int - target) < abs(r_best.Deta_int - target):
            r_best = r_mid
        if abs(r_mid.Deta_int - target) <= tol:
            return r_mid

        if (r_lo.Deta_int - target) * (r_mid.Deta_int - target) <= 0:
            hi = mid
            r_hi = r_mid
        else:
            lo = mid
            r_lo = r_mid

    return r_best


def _print_runs(rows: List[Row]):
    print("\n=== LOCALITY RUNS SUMMARY (bulk, current-based) ===")
    print("mode, kind, case, eta_peak, yg, sigma, seed, jitter, t_onset_curl, t_onset_div, "
          "jc_late/OFF, jd_late/OFF, Deta_int, Deta_target, Deta_err, eff_curl, eff_div")
    for r in rows:
        print(
            f"{r.mode}, {r.kind}, {r.case}, {r.eta_peak:.6g}, {r.yg:.3f}, {r.sigma:.3f}, {r.seed:d}, {r.jitter:g}, "
            f"{_fmt_onset(r.t_onset_curl, r.Tmax)}, {_fmt_onset(r.t_onset_div, r.Tmax)}, "
            f"{r.jc_ratio:.6e}, {r.jd_ratio:.6e}, {r.Deta_int:.6e}, {r.Deta_target:.6e}, {r.Deta_err:.6e}, "
            f"{r.eff_curl:.6e}, {r.eff_div:.6e}"
        )


def _aggregate(rows: List[Row]) -> List[Dict[str, Any]]:
    # Exclude OFF and REF from aggregation by design
    use = [r for r in rows if r.kind == "RUN"]
    groups: Dict[Tuple[str, float, float], List[Row]] = {}
    for r in use:
        groups.setdefault((r.mode, r.sigma, r.yg), []).append(r)

    out = []
    for (mode, sigma, yg), rs in sorted(groups.items(), key=lambda k: (k[0][0], k[0][1], k[0][2])):
        def arr(getter):
            return np.array([getter(x) for x in rs], dtype=float)

        jc_m, jc_s = _meanstd_finite(arr(lambda x: x.jc_ratio))
        jd_m, jd_s = _meanstd_finite(arr(lambda x: x.jd_ratio))
        De_m, De_s = _meanstd_finite(arr(lambda x: x.Deta_int))
        efc_m, efc_s = _meanstd_finite(arr(lambda x: x.eff_curl))
        efd_m, efd_s = _meanstd_finite(arr(lambda x: x.eff_div))
        ep_m, ep_s = _meanstd_finite(arr(lambda x: x.eta_peak))
        err_m, err_s = _meanstd_finite(arr(lambda x: x.Deta_err))

        out.append(dict(
            mode=mode, sigma=float(sigma), yg=float(yg), n=len(rs),
            eta_peak_mean=ep_m, eta_peak_std=ep_s,
            jc_ratio_mean=jc_m, jc_ratio_std=jc_s,
            jd_ratio_mean=jd_m, jd_ratio_std=jd_s,
            Deta_int_mean=De_m, Deta_int_std=De_s,
            Deta_err_mean=err_m, Deta_err_std=err_s,
            eff_curl_mean=efc_m, eff_curl_std=efc_s,
            eff_div_mean=efd_m, eff_div_std=efd_s,
        ))
    return out


def _print_agg(agg: List[Dict[str, Any]]):
    print("\n=== LOCALITY AGG SUMMARY (mean±std over seeds/jitter) ===")
    print("mode, sigma, yg, n, eta_peak_mean, eta_peak_std, jc_ratio_mean, jc_ratio_std, jd_ratio_mean, jd_ratio_std, "
          "Deta_int_mean, Deta_int_std, Deta_err_mean, Deta_err_std, eff_curl_mean, eff_curl_std, eff_div_mean, eff_div_std")
    for r in agg:
        print(
            f"{r['mode']}, {r['sigma']:.3f}, {r['yg']:.3f}, {r['n']:d}, "
            f"{r['eta_peak_mean']:.6e}, {r['eta_peak_std']:.6e}, "
            f"{r['jc_ratio_mean']:.6e}, {r['jc_ratio_std']:.6e}, "
            f"{r['jd_ratio_mean']:.6e}, {r['jd_ratio_std']:.6e}, "
            f"{r['Deta_int_mean']:.6e}, {r['Deta_int_std']:.6e}, "
            f"{r['Deta_err_mean']:.6e}, {r['Deta_err_std']:.6e}, "
            f"{r['eff_curl_mean']:.6e}, {r['eff_curl_std']:.6e}, "
            f"{r['eff_div_mean']:.6e}, {r['eff_div_std']:.6e}"
        )


def _write_csv(path: str, dict_rows: List[Dict[str, Any]]):
    if not dict_rows:
        return
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(dict_rows[0].keys()))
        w.writeheader()
        for r in dict_rows:
            w.writerow(r)


def _plot(rows: List[Row], agg: List[Dict[str, Any]], fig_out: str):
    modes = sorted(set(r["mode"] for r in agg))
    sigmas = sorted(set(r["sigma"] for r in agg))

    nrows = len(modes)
    fig, ax = plt.subplots(nrows, 3, figsize=(12, 3.4 * nrows), sharex=True)
    if nrows == 1:
        ax = np.array([ax])

    for i, mode in enumerate(modes):
        for sigma in sigmas:
            a = [r for r in agg if r["mode"] == mode and r["sigma"] == sigma]
            if not a:
                continue
            xs = [r["yg"] for r in a]
            yc = [r["jc_ratio_mean"] for r in a]
            yd = [r["jd_ratio_mean"] for r in a]
            ef = [r["eff_curl_mean"] for r in a]
            ax[i, 0].plot(xs, yc, marker="o", label=f"sigma={sigma:g}")
            ax[i, 1].plot(xs, yd, marker="o", label=f"sigma={sigma:g}")
            ax[i, 2].plot(xs, ef, marker="o", label=f"sigma={sigma:g}")

        ax[i, 0].set_title(f"{mode}: curl ratio vs yg")
        ax[i, 1].set_title(f"{mode}: div ratio vs yg")
        ax[i, 2].set_title(f"{mode}: eff_curl vs yg")
        for j in range(3):
            ax[i, j].grid(True, alpha=0.3)
            ax[i, j].set_xlabel("yg")

        ax[i, 0].set_ylabel("ratio vs OFF")
        ax[i, 2].set_ylabel("(1-ratio)/Deta_int")

        h, lab = ax[i, 0].get_legend_handles_labels()
        if h:
            ax[i, 0].legend()

    fig.suptitle("Locality summary (ratios and efficiency)")
    fig.tight_layout()
    fig.savefig(fig_out, dpi=160)

    # --- ratio vs work: plot per-run cloud + mean markers ---
    run_rows = [r for r in rows if r.kind == "RUN"]
    fig2, ax2 = plt.subplots(1, 2, figsize=(10, 4))
    for mode in modes:
        rr = [r for r in run_rows if r.mode == mode]
        x = np.array([r.Deta_int for r in rr], dtype=float)
        y1 = np.array([r.jc_ratio for r in rr], dtype=float)
        y2 = np.array([r.jd_ratio for r in rr], dtype=float)
        ax2[0].scatter(x, y1, label=f"{mode} runs", alpha=0.35)
        ax2[1].scatter(x, y2, label=f"{mode} runs", alpha=0.35)

        aa = [r for r in agg if r["mode"] == mode]
        xm = np.array([r["Deta_int_mean"] for r in aa], dtype=float)
        y1m = np.array([r["jc_ratio_mean"] for r in aa], dtype=float)
        y2m = np.array([r["jd_ratio_mean"] for r in aa], dtype=float)
        ax2[0].scatter(xm, y1m, marker="D", label=f"{mode} means")
        ax2[1].scatter(xm, y2m, marker="D", label=f"{mode} means")

    ax2[0].set_title("curl ratio vs Deta_int")
    ax2[1].set_title("div ratio vs Deta_int")
    ax2[0].set_xlabel("Deta_int")
    ax2[1].set_xlabel("Deta_int")
    ax2[0].set_ylabel("ratio vs OFF")
    ax2[1].set_ylabel("ratio vs OFF")
    ax2[0].grid(True, alpha=0.3)
    ax2[1].grid(True, alpha=0.3)
    ax2[0].legend()
    ax2[1].legend()
    fig2.tight_layout()
    fig2.savefig(fig_out.replace(".png", "_vs_work.png"), dpi=160)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["fixed", "match", "both"], default="match")

    ap.add_argument("--eta_peak", type=float, default=0.012)
    ap.add_argument("--eta_ref", type=float, default=None)
    ap.add_argument("--eta_bounds", type=str, default="0.002,0.2")
    ap.add_argument("--rtol", type=float, default=0.03)
    ap.add_argument("--atol", type=float, default=0.0)
    ap.add_argument("--max_iter", type=int, default=12)

    ap.add_argument("--sigmas", type=str, required=True)
    ap.add_argument("--yg_list", type=str, required=True)
    ap.add_argument("--yg_ref", type=float, default=None)

    ap.add_argument("--t0", type=float, default=3.0)
    ap.add_argument("--tau", type=float, default=3.0)
    ap.add_argument("--tmax", type=float, default=12.0)
    ap.add_argument("--dt", type=float, default=0.005)

    ap.add_argument("--seeds", type=str, default="1,2,3")
    ap.add_argument("--jitters", type=str, default="0.02,0.03")

    ap.add_argument("--sigma_off", type=float, default=1.2)

    ap.add_argument("--t_baseline_max", type=float, default=3.0)
    ap.add_argument("--t_late", type=float, default=10.0)
    ap.add_argument("--nsig", type=float, default=6.0)
    ap.add_argument("--hold", type=int, default=6)

    ap.add_argument("--csv_runs", type=str, default="wp2_loc_runs.csv")
    ap.add_argument("--csv_agg", type=str, default="wp2_loc_agg.csv")
    ap.add_argument("--fig_out", type=str, default="wp2_loc.png")
    args = ap.parse_args()

    sigmas = _parse_list(args.sigmas, float)
    ygs = _parse_list(args.yg_list, float)
    seeds = _parse_list(args.seeds, int)
    jitters = _parse_list(args.jitters, float)

    if args.yg_ref is None:
        args.yg_ref = float(ygs[0])
    eta_ref = args.eta_peak if args.eta_ref is None else args.eta_ref
    eta_min, eta_max = _parse_list(args.eta_bounds, float)

    Nx = 512 if getattr(base, "use_gpu", False) else 256
    nsteps = int(round(args.tmax / args.dt))

    P_common = dict(
        Nx=Nx, Ny=Nx,
        Lx=24.0, Ly=24.0,
        dt=args.dt, nsteps=nsteps,
        g=1.0, mu=1.0, V0=0.0,
        n_vort=24, y_sep=2.0, core=0.25,
        t0=args.t0, tau=args.tau,
        xg=0.0,
        gate_alpha=0.25,
        DIFF_LIMIT=0.35,
        diag_stride=10,
    )

    onset_cfg = dict(
        t_baseline_max=args.t_baseline_max,
        t_late=args.t_late,
        nsig=args.nsig,
        hold=args.hold,
    )

    rows: List[Row] = []
    cache: Dict[Tuple[int, float, float, float, float], Row] = {}

    for seed in seeds:
        for jitter in jitters:
            P = dict(P_common, seed=int(seed), jitter=float(jitter))

            # OFF baseline (needed for ratios)
            off = _run_case(
                name="OFF", P=P, eta_peak=0.0, yg=0.0, sigma=args.sigma_off,
                onset_cfg=onset_cfg, off_jc_late=1.0, off_jd_late=1.0,
                mode="fixed", kind="OFF"
            )
            rows.append(off)
            off_jc_late = base.summarize(base.run_sim_timeseries("OFF_tmp", **{**P, "eta_peak": 0.0, "yg": 0.0, "sigma": args.sigma_off}), **onset_cfg)["jc_late_mean"]
            off_jd_late = base.summarize(base.run_sim_timeseries("OFF_tmp2", **{**P, "eta_peak": 0.0, "yg": 0.0, "sigma": args.sigma_off}), **onset_cfg)["jd_late_mean"]

            # FIXED
            if args.mode in ("fixed", "both"):
                for sigma in sigmas:
                    for yg in ygs:
                        rows.append(_run_case(
                            name=f"fixed_yg{yg:g}_sig{sigma:g}",
                            P=P, eta_peak=args.eta_peak, yg=yg, sigma=sigma,
                            onset_cfg=onset_cfg,
                            off_jc_late=off_jc_late, off_jd_late=off_jd_late,
                            mode="fixed", kind="RUN"
                        ))

            # MATCH
            if args.mode in ("match", "both"):
                ref = _run_case(
                    name=f"ref_match_yg{args.yg_ref:g}_sig{sigmas[0]:g}",
                    P=P, eta_peak=float(eta_ref), yg=float(args.yg_ref), sigma=float(sigmas[0]),
                    onset_cfg=onset_cfg,
                    off_jc_late=off_jc_late, off_jd_late=off_jd_late,
                    mode="match", kind="REF"
                )
                rows.append(ref)
                target = ref.Deta_int

                for sigma in sigmas:
                    for yg in ygs:
                        # IMPORTANT: don’t run an extra bisection at yg_ref; use eta_ref there
                        if float(yg) == float(args.yg_ref):
                            rows.append(_run_case(
                                name=f"match_yg{yg:g}_sig{sigma:g}",
                                P=P, eta_peak=float(eta_ref), yg=yg, sigma=sigma,
                                onset_cfg=onset_cfg,
                                off_jc_late=off_jc_late, off_jd_late=off_jd_late,
                                mode="match", kind="RUN",
                                Deta_target=target
                            ))
                            continue

                        rows.append(_match_eta_for_work(
                            P=P, yg=yg, sigma=sigma,
                            onset_cfg=onset_cfg,
                            off_jc_late=off_jc_late, off_jd_late=off_jd_late,
                            target=target, eta_min=float(eta_min), eta_max=float(eta_max),
                            rtol=float(args.rtol), atol=float(args.atol),
                            max_iter=int(args.max_iter),
                            cache=cache,
                            name_prefix=f"match_yg{yg:g}_sig{sigma:g}",
                        ))

    _print_runs(rows)
    _write_csv(args.csv_runs, [asdict(r) for r in rows])
    print(f"\n[Wrote] {args.csv_runs}")

    agg = _aggregate(rows)
    _print_agg(agg)
    _write_csv(args.csv_agg, agg)
    print(f"\n[Wrote] {args.csv_agg}")

    _plot(rows, agg, args.fig_out)
    print(f"[Wrote] {args.fig_out} and {args.fig_out.replace('.png','_vs_work.png')}")


if __name__ == "__main__":
    main()
