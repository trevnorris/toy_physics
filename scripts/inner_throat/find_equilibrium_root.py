import json
import math
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


def load_best_point(path: Path):
    if not path.exists():
        return None
    with open(path, "r") as f:
        data = json.load(f)
    rows = data.get("rows") if isinstance(data, dict) else None
    if not rows:
        return None
    return min(rows, key=lambda r: abs(r.get("|R|", 1e30)))


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

    scan_path = repo_root / "paper7_scan_results.json"
    best = load_best_point(scan_path)
    if best:
        a = float(best.get("a", 0.8))
        L = float(best.get("L", 2.75))
    else:
        a = 0.8166667
        L = 2.75

    Ntarget = 1.0
    dtau = 1e-3
    max_steps = 1200
    max_steps_probe = 300
    tol_imag = 1e-7
    report_every = 50
    tol_root = 1e-4
    max_iter = 6

    def solve_ground_state(a_val, L_val, psi0=None, max_steps_override=None):
        steps = max_steps if max_steps_override is None else max_steps_override
        V = build_vconf_F1(grid, a_val, L_val, P)
        psi = find_ground_state_imagtime(
            grid, P, V,
            Ntarget=Ntarget,
            dtau=dtau,
            max_steps=steps,
            tol=tol_imag,
            report_every=report_every,
            psi0=psi0,
        )
        return psi

    def eval_residuals(a_val, L_val, psi0=None, max_steps_override=None):
        psi = solve_ground_state(a_val, L_val, psi0=psi0, max_steps_override=max_steps_override)
        Ra, RL, info = force_residuals(
            grid, P, psi,
            a=a_val, L=L_val,
            deriv_method="analytic",
        )
        Rnorm = math.sqrt(Ra * Ra + RL * RL)
        return psi, Ra, RL, Rnorm, info

    psi, Ra, RL, Rnorm, info = eval_residuals(a, L, psi0=None)
    print(f"[Root] start (a={a:.6f}, L={L:.6f}) R=({Ra:+.3e}, {RL:+.3e}) |R|={Rnorm:.3e}")

    for it in range(1, max_iter + 1):
        da = max(grid.dx, 0.5 * P.dr)
        dL = max(grid.dw, 0.5 * P.dw)

        psi_ap, Ra_ap, RL_ap, _, _ = eval_residuals(a + da, L, psi0=psi, max_steps_override=max_steps_probe)
        psi_am, Ra_am, RL_am, _, _ = eval_residuals(a - da, L, psi0=psi, max_steps_override=max_steps_probe)
        psi_Lp, Ra_Lp, RL_Lp, _, _ = eval_residuals(a, L + dL, psi0=psi, max_steps_override=max_steps_probe)
        psi_Lm, Ra_Lm, RL_Lm, _, _ = eval_residuals(a, L - dL, psi0=psi, max_steps_override=max_steps_probe)

        J = np.array([
            [(Ra_ap - Ra_am) / (2.0 * da), (Ra_Lp - Ra_Lm) / (2.0 * dL)],
            [(RL_ap - RL_am) / (2.0 * da), (RL_Lp - RL_Lm) / (2.0 * dL)],
        ], dtype=float)
        Rvec = np.array([Ra, RL], dtype=float)

        try:
            delta = -np.linalg.solve(J, Rvec)
        except np.linalg.LinAlgError:
            print("[Root] Jacobian solve failed; stopping.")
            break

        accepted = False
        for scale in [1.0, 0.5, 0.25, 0.1, 0.05]:
            a_try = a + scale * float(delta[0])
            L_try = L + scale * float(delta[1])
            if a_try <= 0.0 or L_try <= 0.0:
                continue
            psi_try, Ra_try, RL_try, Rnorm_try, info_try = eval_residuals(a_try, L_try, psi0=psi)
            if Rnorm_try < Rnorm:
                a, L = a_try, L_try
                psi, Ra, RL, Rnorm, info = psi_try, Ra_try, RL_try, Rnorm_try, info_try
                accepted = True
                break

        if not accepted:
            print("[Root] line search failed to improve residual; stopping.")
            break

        print(f"[Root] iter={it} (a={a:.6f}, L={L:.6f}) R=({Ra:+.3e}, {RL:+.3e}) |R|={Rnorm:.3e}")
        if Rnorm < tol_root:
            break

    out = {
        "a": a,
        "L": L,
        "Ra": Ra,
        "RL": RL,
        "|R|": Rnorm,
        "Ia": info.get("Ia"),
        "IL": info.get("IL"),
        "IL_plus": info.get("IL_plus"),
        "IL_minus": info.get("IL_minus"),
        "IL_split_rel": info.get("IL_split_rel"),
        "meta": {
            "grid": {"N": grid.Nx, "L": grid.Lx},
            "solver": {
                "Ntarget": Ntarget,
                "dtau": dtau,
                "max_steps": max_steps,
                "tol_imag": tol_imag,
                "tol_root": tol_root,
            },
            "gpu": bool(use_gpu),
            "seed": "scan" if best else "default",
        },
    }

    out_path = repo_root / "paper7_root_result.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)

    print(f"[Root] saved {out_path}")


if __name__ == "__main__":
    main()
