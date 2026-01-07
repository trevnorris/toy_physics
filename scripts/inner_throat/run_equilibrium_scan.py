import json
import time
from pathlib import Path

import numpy as np

from paper7_4d_gnls_gpu import (
    select_backend,
    Params,
    Grid4D,
    build_vconf_F1,
    find_ground_state_imagtime,
    force_residuals,
)

# -----------------------
# Small table helpers
# -----------------------
def print_table(rows, cols, max_rows=50, title=None):
    if title:
        print("\n" + title)

    # Try pandas first (nice formatting), otherwise fallback to manual
    try:
        import pandas as pd  # type: ignore
        df = pd.DataFrame(rows)
        df = df[cols]
        if len(df) > max_rows:
            print(df.head(max_rows).to_string(index=False))
            print(f"... ({len(df) - max_rows} more rows)")
        else:
            print(df.to_string(index=False))
        return
    except Exception:
        pass

    # Manual fixed-width formatting
    widths = {c: max(len(c), *(len(f"{r.get(c,'')}") for r in rows)) for c in cols}
    header = "  ".join(c.ljust(widths[c]) for c in cols)
    sep = "  ".join("-" * widths[c] for c in cols)
    print(header)
    print(sep)

    for i, r in enumerate(rows[:max_rows]):
        line = "  ".join(f"{r.get(c,'')}".ljust(widths[c]) for c in cols)
        print(line)

    if len(rows) > max_rows:
        print(f"... ({len(rows) - max_rows} more rows)")


def main():
    # -----------------------
    # Backend + config
    # -----------------------
    xp, fft, asnumpy, use_gpu = select_backend(prefer_gpu=True)

    cfg_path = Path("paper7_defs.json")
    P = Params.from_json(str(cfg_path))

    # -----------------------
    # Grid (scale up once stable)
    # -----------------------
    grid = Grid4D(
        xp, fft,
        Nx=64, Ny=64, Nz=64, Nw=64,
        Lx=16.0, Ly=16.0, Lz=16.0, Lw=16.0,
        dtype_real=P.dtype_real
    )

    # -----------------------
    # Scan ranges
    # -----------------------
    a_list = np.linspace(0.4, 0.9, 7)
    L_list = np.linspace(2.0, 4.0, 9)

    # -----------------------
    # Solver knobs
    # -----------------------
    Ntarget = 1.0
    dtau = 1e-3
    max_steps = 800
    tol = 1e-7
    report_every = 50

    deriv_method = "analytic"
    if deriv_method == "analytic":
        da_fd = None
        dL_fd = None
    else:
        da_fd = 0.01
        dL_fd = 0.01

    # -----------------------
    # Run scan
    # -----------------------
    psi_prev = None
    rows = []
    t_start = time.time()

    for L in L_list:
        for a in a_list:
            V = build_vconf_F1(grid, float(a), float(L), P)

            psi = find_ground_state_imagtime(
                grid, P, V,
                Ntarget=Ntarget,
                dtau=dtau,
                max_steps=max_steps,
                tol=tol,
                report_every=report_every,
                psi0=psi_prev
            )

            Ra, RL, info = force_residuals(
                grid, P, psi,
                a=float(a), L=float(L),
                deriv_method=deriv_method,
                da=da_fd, dL=dL_fd
            )

            res_norm = float(np.sqrt(Ra * Ra + RL * RL))

            row = {
                "a": float(a),
                "L": float(L),
                "Ra": float(Ra),
                "RL": float(RL),
                "|R|": res_norm,
                "Ia": float(info["Ia"]),
                "IL": float(info["IL"]),
                "IL_plus": float(info["IL_plus"]),
                "IL_minus": float(info["IL_minus"]),
                "IL_split_rel": float(info["IL_split_rel"]),
                "dEda": float(info["dEda"]),
                "dEdL": float(info["dEdL"]),
            }
            rows.append(row)

            print(f"(a,L)=({a:.3f},{L:.3f})  Ra={Ra:+.3e}  RL={RL:+.3e}  |R|={res_norm:.3e}")

            psi_prev = psi  # warm-start

    elapsed = time.time() - t_start
    print(f"\n[Scan] Completed {len(rows)} points in {elapsed:.1f}s (gpu={use_gpu}).")

    # -----------------------
    # Sort + summarize
    # -----------------------
    rows_sorted = sorted(rows, key=lambda r: r["|R|"])
    best = rows_sorted[0]
    print("\n[Best] lowest |R| point:")
    print(json.dumps(best, indent=2))

    def edge_suggestion(best, a_list, L_list, pad_frac=0.25):
        a_min, a_max = float(np.min(a_list)), float(np.max(a_list))
        L_min, L_max = float(np.min(L_list)), float(np.max(L_list))
        a_rng = a_max - a_min
        L_rng = L_max - L_min

        on_min_a = abs(best["a"] - a_min) < 1e-12
        on_max_a = abs(best["a"] - a_max) < 1e-12
        on_min_L = abs(best["L"] - L_min) < 1e-12
        on_max_L = abs(best["L"] - L_max) < 1e-12

        if not (on_min_a or on_max_a or on_min_L or on_max_L):
            return None

        # Expand the box outward in the directions that hit the boundary.
        new_a_min = a_min - (pad_frac * a_rng) if on_min_a else a_min
        new_a_max = a_max + (pad_frac * a_rng) if on_max_a else a_max
        new_L_min = L_min - (pad_frac * L_rng) if on_min_L else L_min
        new_L_max = L_max + (pad_frac * L_rng) if on_max_L else L_max

        # Clamp at >0 for physical positivity
        new_a_min = max(1e-3, new_a_min)
        new_L_min = max(1e-3, new_L_min)

        return {
            "reason": f"best point is on scan boundary (a_min={on_min_a}, a_max={on_max_a}, L_min={on_min_L}, L_max={on_max_L})",
            "suggested_box": {"a_min": new_a_min, "a_max": new_a_max, "L_min": new_L_min, "L_max": new_L_max},
            "suggested_counts": {"Na": len(a_list), "NL": len(L_list)},
        }

    suggest = edge_suggestion(best, a_list, L_list, pad_frac=0.5)
    if suggest:
        print("\n[Suggestion] The best point is on the scan boundary.")
        print("  " + suggest["reason"])
        box = suggest["suggested_box"]
        print(f"  Next scan box: a in [{box['a_min']:.3g}, {box['a_max']:.3g}], L in [{box['L_min']:.3g}, {box['L_max']:.3g}]")


    # Print a nice table of the best points
    cols = ["a", "L", "Ra", "RL", "|R|", "Ia", "IL", "IL_split_rel"]
    print_table(rows_sorted, cols, max_rows=25, title="[Top 25] Sorted by |R| (best first)")

    # -----------------------
    # Save outputs (npy + csv + json)
    # -----------------------
    out_dir = Path(".")
    npy_path = out_dir / "paper7_scan_results.npy"
    csv_path = out_dir / "paper7_scan_results.csv"
    json_path = out_dir / "paper7_scan_results.json"

    # NPY: raw numeric array (still useful)
    arr = np.array([[r["a"], r["L"], r["Ra"], r["RL"], r["|R|"], r["Ia"], r["IL"], r["dEda"], r["dEdL"]] for r in rows],
                   dtype=float)
    np.save(npy_path, arr)

    # CSV: human-friendly
    import csv
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # JSON: includes metadata + rows
    payload = {
        "meta": {
            "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "gpu": bool(use_gpu),
            "grid": {"Nx": grid.Nx, "Ny": grid.Ny, "Nz": grid.Nz, "Nw": grid.Nw,
                     "Lx": grid.Lx, "Ly": grid.Ly, "Lz": grid.Lz, "Lw": grid.Lw},
            "solver": {"Ntarget": Ntarget, "dtau": dtau, "max_steps": max_steps, "tol": tol,
                       "report_every": report_every, "deriv_method": deriv_method,
                       "da_fd": da_fd, "dL_fd": dL_fd},
            "freezeSHA256": getattr(P, "freezeSHA256", ""),
        },
        "rows": rows_sorted,
    }
    with open(json_path, "w") as f:
        json.dump(payload, f, indent=2)

    print("\n[Saved]")
    print(f"  {npy_path}")
    print(f"  {csv_path}")
    print(f"  {json_path}")


if __name__ == "__main__":
    main()
