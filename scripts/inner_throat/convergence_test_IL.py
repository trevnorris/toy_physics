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

    a0 = 0.8166667
    L0 = 2.75
    grid_sizes = [48, 56, 64]

    rows = []
    for N in grid_sizes:
        grid = Grid4D(
            xp, fft,
            Nx=N, Ny=N, Nz=N, Nw=N,
            Lx=16.0, Ly=16.0, Lz=16.0, Lw=16.0,
            dtype_real=P.dtype_real,
        )

        V = build_vconf_F1(grid, a0, L0, P)
        psi = find_ground_state_imagtime(
            grid, P, V,
            Ntarget=1.0,
            dtau=1e-3,
            max_steps=1200,
            tol=1e-7,
            report_every=50,
            psi0=None,
        )

        Ra, RL, info = force_residuals(
            grid, P, psi,
            a=a0, L=L0,
            deriv_method="analytic",
        )

        rows.append({
            "N": f"{N}",
            "Ra": f"{Ra:+.3e}",
            "RL": f"{RL:+.3e}",
            "Ia": f"{info['Ia']:+.3e}",
            "IL": f"{info['IL']:+.3e}",
            "IL_split": f"{info['IL_split_rel']:.3e}",
        })

    print_table(
        rows,
        cols=["N", "Ra", "RL", "Ia", "IL", "IL_split"],
        title=f"[convergence_test_IL] (a={a0:.7f}, L={L0:.3f}, gpu={use_gpu})",
    )


if __name__ == "__main__":
    main()
