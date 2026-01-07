from pathlib import Path

from paper7_4d_gnls_gpu import (
    select_backend,
    Params,
    Grid4D,
    build_vconf_F1,
    find_ground_state_imagtime,
    force_residuals,
)


def print_table(rows, cols, title=None):
    if title:
        print("\n" + title)
    widths = {c: max(len(c), *(len(f"{r.get(c,'')}") for r in rows)) for c in cols}
    header = "  ".join(c.ljust(widths[c]) for c in cols)
    sep = "  ".join("-" * widths[c] for c in cols)
    print(header)
    print(sep)
    for r in rows:
        line = "  ".join(f"{r.get(c,'')}".ljust(widths[c]) for c in cols)
        print(line)


def main():
    xp, fft, asnumpy, use_gpu = select_backend(prefer_gpu=True)

    repo_root = Path(__file__).resolve().parents[2]
    cfg_path = repo_root / "paper7_defs.json"
    P = Params.from_json(str(cfg_path))

    grid = Grid4D(
        xp, fft,
        Nx=64, Ny=64, Nz=64, Nw=64,
        Lx=16.0, Ly=16.0, Lz=16.0, Lw=16.0,
        dtype_real=P.dtype_real,
    )

    a0 = 0.8166667
    L0 = 2.75

    V = build_vconf_F1(grid, a0, L0, P)
    psi = find_ground_state_imagtime(
        grid, P, V,
        Ntarget=1.0,
        dtau=1e-3,
        max_steps=1000,
        tol=1e-7,
        report_every=50,
        psi0=None,
    )

    rows = []

    Ra, RL, info = force_residuals(
        grid, P, psi,
        a=a0, L=L0,
        deriv_method="analytic",
    )
    rows.append({
        "method": "analytic",
        "dL": "-",
        "IL": f"{info['IL']:+.6e}",
        "-IL": f"{-info['IL']:+.6e}",
        "RL": f"{RL:+.6e}",
        "IL_split": f"{info['IL_split_rel']:.3e}",
        "Lhalf_err": f"{info.get('IL_Lhalf_err', 0.0):.3e}",
    })

    dL_list = [0.05, 0.1, 0.2, 0.3]
    da_fd = max(grid.dx, 0.5 * P.dr)

    for dL in dL_list:
        Ra, RL, info = force_residuals(
            grid, P, psi,
            a=a0, L=L0,
            deriv_method="fd64",
            da=da_fd,
            dL=dL,
            deriv_dtype="float64",
        )
        rows.append({
            "method": "fd64",
            "dL": f"{dL:.3f}",
            "IL": f"{info['IL']:+.6e}",
            "-IL": f"{-info['IL']:+.6e}",
            "RL": f"{RL:+.6e}",
            "IL_split": f"{info['IL_split_rel']:.3e}",
            "Lhalf_err": "-",
        })

    print_table(
        rows,
        cols=["method", "dL", "IL", "-IL", "RL", "IL_split", "Lhalf_err"],
        title=f"[diagnose_dVdL] (a={a0:.7f}, L={L0:.3f}, gpu={use_gpu})",
    )


if __name__ == "__main__":
    main()
