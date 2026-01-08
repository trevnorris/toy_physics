"""
wp2_onset_scan_v1.py

- Uses current-based bulk diagnostics (no 1/rho):
    j = Im(conj(psi) * grad psi)
    curl j, div j
- Runs:
    OFF
    ON_center (eta_peak=0.08, yg=0.0, sigma=1.2)
    ON_displaced (eta_peak=0.08, yg=11.0, sigma=0.6)   # stronger "off-sheet" control
    plus a centered eta scan
- Outputs:
    1) Plain-text table (with ">Tmax" instead of NaN onset)
    2) Summary plot (onset times + late-time means vs eta)
"""

import numpy as np
import matplotlib.pyplot as plt

# --- GPU SETUP ---
try:
    import cupy as cp
    xp = cp
    use_gpu = True
    print("[System] CuPy detected. GPU Mode Engaged.")
except ImportError:
    xp = np
    use_gpu = False
    print("[System] CuPy not found. CPU Mode (Slower).")

def to_cpu(a):
    if use_gpu:
        return a.get()
    return a

# =========================
# Helpers
# =========================
def wrap_delta(d, L):
    return (d + 0.5 * L) % L - 0.5 * L

def make_kgrid(Nx, Ny, Lx, Ly, xp):
    kx = 2.0 * np.pi * xp.fft.fftfreq(Nx, d=Lx / Nx)
    ky = 2.0 * np.pi * xp.fft.fftfreq(Ny, d=Ly / Ny)
    KX, KY = xp.meshgrid(kx, ky, indexing="ij")
    K2 = KX**2 + KY**2
    kmax2 = float(to_cpu(xp.max(K2)))
    return KX, KY, K2, kmax2

def grad_spectral(f, KX, KY, xp):
    F = xp.fft.fft2(f)
    fx = xp.fft.ifft2(1j * KX * F)
    fy = xp.fft.ifft2(1j * KY * F)
    return fx, fy

def div_spectral(fx, fy, KX, KY, xp):
    Fx = xp.fft.fft2(fx)
    Fy = xp.fft.fft2(fy)
    return xp.fft.ifft2(1j * KX * Fx + 1j * KY * Fy)

def masked_mean(a, mask):
    num = xp.sum(a * mask)
    den = xp.sum(mask)
    if float(to_cpu(den)) <= 0.0:
        return xp.nan
    return num / den

def tukey_window(t, t0, tau, alpha=0.25):
    if (t < t0) or (t > t0 + tau):
        return 0.0
    s = (t - t0) / tau  # [0,1]
    a = float(alpha)
    if a <= 0.0:
        return 1.0
    if a >= 1.0:
        return 0.5 * (1.0 - np.cos(2.0 * np.pi * s))
    if s < 0.5 * a:
        return 0.5 * (1.0 - np.cos(2.0 * np.pi * s / a))
    if s <= 1.0 - 0.5 * a:
        return 1.0
    return 0.5 * (1.0 - np.cos(2.0 * np.pi * (1.0 - s) / a))

def eta_gaussian(X, Y, Lx, Ly, eta_peak, xg, yg, sigma, t, t0, tau, alpha=0.25):
    w = tukey_window(float(t), float(t0), float(tau), alpha=alpha)
    if w == 0.0 or eta_peak == 0.0:
        return 0.0 * X
    dx = wrap_delta(X - xg, Lx)
    dy = wrap_delta(Y - yg, Ly)
    return (eta_peak * w) * xp.exp(-(dx*dx + dy*dy) / (2.0 * sigma * sigma))

def imprint_vortex_rows(X, Y, Lx, Ly, *,
                        n_vort=24, y_sep=2.0, core=0.25, jitter=0.0, seed=0):
    rng = np.random.default_rng(seed)
    psi = xp.ones_like(X, dtype=xp.complex128)

    xs = np.linspace(-0.5*Lx, 0.5*Lx, n_vort, endpoint=False)
    xs = xs + jitter * rng.standard_normal(size=xs.shape)

    y1 = +0.5 * y_sep
    y2 = -0.5 * y_sep

    for x0 in xs:
        # + vortex
        dx = wrap_delta(X - x0, Lx)
        dy = wrap_delta(Y - y1, Ly)
        r2 = dx*dx + dy*dy
        phase = xp.arctan2(dy, dx)
        amp = xp.sqrt(r2 / (r2 + core*core))
        psi *= amp * xp.exp(1j * (+1.0) * phase)

        # - vortex
        dx = wrap_delta(X - x0, Lx)
        dy = wrap_delta(Y - y2, Ly)
        r2 = dx*dx + dy*dy
        phase = xp.arctan2(dy, dx)
        amp = xp.sqrt(r2 / (r2 + core*core))
        psi *= amp * xp.exp(1j * (-1.0) * phase)

    rho = xp.abs(psi)**2
    psi = psi / xp.sqrt(xp.mean(rho))
    return psi

# =========================
# Current-based diagnostics (no 1/rho)
# =========================
def current_from_psi(psi, KX, KY, xp):
    dpsi_dx, dpsi_dy = grad_spectral(psi, KX, KY, xp)
    rho = xp.abs(psi)**2
    jx = xp.imag(xp.conj(psi) * dpsi_dx)
    jy = xp.imag(xp.conj(psi) * dpsi_dy)
    return jx, jy, rho, dpsi_dx, dpsi_dy

def curl_and_div(fx, fy, KX, KY, xp):
    dfy_dx, _ = grad_spectral(fy, KX, KY, xp)
    _, dfx_dy = grad_spectral(fx, KX, KY, xp)
    curl_z = dfy_dx - dfx_dy
    div_v = xp.real(div_spectral(fx, fy, KX, KY, xp))
    return xp.real(curl_z), div_v

# =========================
# Run + record time series
# =========================
def run_sim_timeseries(case_name,
                       *,
                       Nx, Ny, Lx, Ly,
                       dt, nsteps,
                       g, mu, V0,
                       n_vort, y_sep, core, jitter, seed,
                       eta_peak, t0, tau, sigma, xg, yg,
                       gate_alpha, DIFF_LIMIT,
                       RHO_CORE=0.2, RHO_BULK=0.6,
                       diag_stride=10):

    x = xp.linspace(-0.5*Lx, 0.5*Lx, Nx, endpoint=False)
    y = xp.linspace(-0.5*Ly, 0.5*Ly, Ny, endpoint=False)
    X, Y = xp.meshgrid(x, y, indexing="ij")
    KX, KY, K2, kmax2 = make_kgrid(Nx, Ny, Lx, Ly, xp)

    psi = imprint_vortex_rows(X, Y, Lx, Ly,
                              n_vort=n_vort, y_sep=y_sep, core=core,
                              jitter=jitter, seed=seed)

    kin_phase = xp.exp(-1j * 0.5 * K2 * dt)

    t_hist = []
    jcurl2_bulk_hist = []
    jdiv2_bulk_hist = []
    Deta_hist = []
    eta_max_hist = []
    core_frac_hist = []

    t = 0.0
    for n in range(1, nsteps + 1):
        # NL half
        rho = xp.abs(psi)**2
        psi *= xp.exp(-1j * (V0 + g * rho - mu) * (0.5 * dt))

        # kinetic full
        psi_hat = xp.fft.fft2(psi)
        psi_hat *= kin_phase
        psi = xp.fft.ifft2(psi_hat)

        # NL half
        rho = xp.abs(psi)**2
        psi *= xp.exp(-1j * (V0 + g * rho - mu) * (0.5 * dt))

        # diffusion gate
        eta_max_val = 0.0
        Deta_val = 0.0
        if eta_peak != 0.0:
            eta = eta_gaussian(X, Y, Lx, Ly, eta_peak, xg, yg, sigma, t, t0, tau, alpha=gate_alpha)
            eta_max_val = float(to_cpu(xp.max(eta)))
            if eta_max_val > 0.0:
                dt_sub_max = DIFF_LIMIT / (eta_max_val * kmax2 + 1e-30)
                n_sub = int(np.ceil(dt / dt_sub_max))
                n_sub = max(1, n_sub)
                dt_sub = dt / n_sub
                for _ in range(n_sub):
                    dpsi_dx, dpsi_dy = grad_spectral(psi, KX, KY, xp)
                    fx = eta * dpsi_dx
                    fy = eta * dpsi_dy
                    psi = psi + dt_sub * div_spectral(fx, fy, KX, KY, xp)

                # gate work proxy at end of step
                dpsi_dx, dpsi_dy = grad_spectral(psi, KX, KY, xp)
                grad2 = (xp.abs(dpsi_dx)**2 + xp.abs(dpsi_dy)**2).astype(xp.float64)
                Deta_val = float(to_cpu(xp.mean(eta.astype(xp.float64) * grad2)))

        t = n * dt

        if (n % diag_stride) == 0:
            jx, jy, rho, _, _ = current_from_psi(psi, KX, KY, xp)
            curl_j, div_j = curl_and_div(jx, jy, KX, KY, xp)

            mask_bulk = (rho > RHO_BULK).astype(xp.float64)
            jcurl2 = (curl_j*curl_j).astype(xp.float64)
            jdiv2 = (div_j*div_j).astype(xp.float64)

            jcurl2_bulk = masked_mean(jcurl2, mask_bulk)
            jdiv2_bulk = masked_mean(jdiv2, mask_bulk)
            core_frac = xp.mean((rho < RHO_CORE).astype(xp.float64))

            t_hist.append(float(t))
            jcurl2_bulk_hist.append(float(to_cpu(jcurl2_bulk)))
            jdiv2_bulk_hist.append(float(to_cpu(jdiv2_bulk)))
            Deta_hist.append(float(Deta_val))
            eta_max_hist.append(float(abs(eta_max_val)))
            core_frac_hist.append(float(to_cpu(core_frac)))

    return {
        "case": case_name,
        "t": np.array(t_hist),
        "jcurl2_bulk": np.array(jcurl2_bulk_hist),
        "jdiv2_bulk": np.array(jdiv2_bulk_hist),
        "Deta": np.array(Deta_hist),
        "eta_max": np.array(eta_max_hist),
        "core_frac": np.array(core_frac_hist),
        "Tmax": float(dt * nsteps),
    }

# =========================
# Onset + summary metrics
# =========================
def onset_time(t, y, t_baseline_max=3.0, nsig=6.0, hold=6):
    """
    Onset = first time after t_baseline_max when y stays above (mu+nsig*sigma)
    for 'hold' consecutive samples. Returns (t_onset, threshold).
    """
    if len(t) == 0:
        return np.nan, np.nan
    base = y[t <= t_baseline_max]
    if len(base) < max(10, hold):
        return np.nan, np.nan
    mu = float(np.mean(base))
    sig = float(np.std(base) + 1e-30)
    thr = mu + nsig * sig

    idx0 = np.searchsorted(t, t_baseline_max)
    for i in range(idx0, len(t) - hold):
        if np.all(y[i:i+hold] > thr):
            return float(t[i]), float(thr)
    return np.nan, float(thr)

def late_mean_window(t, y, t0=10.0, window=2.0):
    """
    Mean over a fixed window [t0, t0+window].
    This makes comparisons fair even if runs have different Tmax.
    """
    if len(t) == 0:
        return np.nan
    t1 = t0 + window
    sel = y[(t >= t0) & (t <= t1)]
    if len(sel) == 0:
        return np.nan
    return float(np.mean(sel))

def summarize(out, t_baseline_max=3.0, t_late=10.0, nsig=6.0, hold=6):
    t = out["t"]
    jc = out["jcurl2_bulk"]
    jd = out["jdiv2_bulk"]
    tc, thr_c = onset_time(t, jc, t_baseline_max=t_baseline_max, nsig=nsig, hold=hold)
    td, thr_d = onset_time(t, jd, t_baseline_max=t_baseline_max, nsig=nsig, hold=hold)

    return {
        "tcurl": tc,
        "tdiv": td,
        "thr_curl": thr_c,
        "thr_div": thr_d,
        "jc_max": float(np.max(jc)) if len(jc) else np.nan,
        "jd_max": float(np.max(jd)) if len(jd) else np.nan,
        "jc_late_mean": late_mean_window(t, jc, t0=t_late, window=2.0),
        "jd_late_mean": late_mean_window(t, jd, t0=t_late, window=2.0),
        "Deta_max": float(np.max(out["Deta"])) if len(out["Deta"]) else np.nan,
        "Tmax": out["Tmax"],
    }

def fmt_onset(t_onset, Tmax):
    if np.isnan(t_onset):
        return f">{Tmax:.2f}"
    return f"{t_onset:.4f}"

# =========================
# Main experiment set
# =========================
if __name__ == "__main__":

    Nx = 512 if use_gpu else 256
    Ny = Nx

    # Shared parameters (match your previous runs)
    P = dict(
        Nx=Nx, Ny=Ny,
        Lx=24.0, Ly=24.0,
        dt=0.005, nsteps=2400,      # Tmax=12
        g=1.0, mu=1.0, V0=0.0,
        n_vort=24, y_sep=2.0, core=0.25, jitter=0.0, seed=1,
        t0=3.0, tau=3.0,
        xg=0.0,
        gate_alpha=0.25,
        DIFF_LIMIT=0.35,
        diag_stride=10,
    )

    # Onset metric settings
    onset_cfg = dict(t_baseline_max=3.0, t_late=10.0, nsig=6.0, hold=6)

    # Define runs:
    runs = []

    # OFF (baseline)
    runs.append(("OFF", dict(eta_peak=0.0, yg=0.0, sigma=1.2)))

    # ON centered (extend time a bit to strengthen the "no onset" claim)
    runs.append(("ON_center", dict(eta_peak=0.08, yg=0.0, sigma=1.2, nsteps=3600)))  # Tmax=18

    # ON displaced (control; also extend)
    runs.append(("ON_displaced", dict(eta_peak=0.08, yg=11.0, sigma=0.6, nsteps=3600)))  # Tmax=18

    # Centered eta scan (extra resolution near the knee ~ 0.01–0.02)
    eta_scan = [
        0.003, 0.005, 0.0075,
        0.010, 0.012, 0.015, 0.018,
        0.020, 0.0225, 0.025,
        0.030, 0.040, 0.060, 0.080
    ]
    # Centered eta scan (keep default Tmax unless you want to extend select points)
    for ep in eta_scan:
        runs.append((f"scan_ep{ep:.4f}", dict(eta_peak=ep, yg=0.0, sigma=1.2)))

    # Execute
    results = []
    for name, upd in runs:
        cfg = dict(P)         # start from defaults
        cfg.update(upd)       # allow overrides like nsteps, dt, etc.

        out = run_sim_timeseries(name, **cfg)
        s = summarize(out, **onset_cfg)
        results.append((name, upd["eta_peak"], upd["yg"], upd["sigma"], s))
        print(f"[Done] {name}")

    def safe_div(a, b):
        if b == 0 or np.isnan(b):
            return np.nan
        return a / b

    # Plain-text table
    Tmax = results[0][4]["Tmax"] if results else 0.0
    # Find OFF reference for normalization
    off_row = next((r for r in results if r[0] == "OFF"), None)
    off_jc_late = off_row[4]["jc_late_mean"] if off_row else np.nan
    off_jd_late = off_row[4]["jd_late_mean"] if off_row else np.nan

    print("\n=== SUMMARY (bulk, current-based) ===")
    print("case, eta_peak, yg, sigma, t_onset_curl, t_onset_div, jc_late_mean, jd_late_mean, "
          "jc_late/ OFF, jd_late/ OFF, jc_max, jd_max, Deta_max")
    for (name, ep, yg, sig, s) in results:
        jc_ratio = safe_div(s["jc_late_mean"], off_jc_late)
        jd_ratio = safe_div(s["jd_late_mean"], off_jd_late)

        print(
            f"{name}, {ep:.4f}, {yg:.3f}, {sig:.3f}, "
            f"{fmt_onset(s['tcurl'], s['Tmax'])}, {fmt_onset(s['tdiv'], s['Tmax'])}, "
            f"{s['jc_late_mean']:.6e}, {s['jd_late_mean']:.6e}, "
            f"{jc_ratio:.6e}, {jd_ratio:.6e}, "
            f"{s['jc_max']:.6e}, {s['jd_max']:.6e}, {s['Deta_max']:.6e}"
        )

    # =========================
    # Summary plot (eta scan only + special markers)
    # =========================
    # Pull scan points
    scan_rows = [(name, ep, yg, sig, s) for (name, ep, yg, sig, s) in results if name.startswith("scan_ep")]
    scan_rows = sorted(scan_rows, key=lambda r: r[1])  # sort by eta

    etas = np.array([r[1] for r in scan_rows], dtype=float)
    tcurl = np.array([r[4]["tcurl"] for r in scan_rows], dtype=float)
    tdiv  = np.array([r[4]["tdiv"]  for r in scan_rows], dtype=float)
    jc_late = np.array([r[4]["jc_late_mean"] for r in scan_rows], dtype=float)
    jd_late = np.array([r[4]["jd_late_mean"] for r in scan_rows], dtype=float)

    # Replace NaNs in onset times for plotting (put them slightly above Tmax)
    Tmax_plot = Tmax
    NO_ONSET_Y = Tmax_plot + 0.35
    tcurl_plot = np.where(np.isnan(tcurl), NO_ONSET_Y, tcurl)
    tdiv_plot  = np.where(np.isnan(tdiv),  NO_ONSET_Y, tdiv)

    # Special points
    def get_summary(name_key):
        for (name, ep, yg, sig, s) in results:
            if name == name_key:
                return ep, s
        return None, None

    ep_off, s_off = get_summary("OFF")
    ep_c, s_c = get_summary("ON_center")
    ep_d, s_d = get_summary("ON_displaced")
    off_jc = s_off["jc_late_mean"] if s_off is not None else None
    off_jd = s_off["jd_late_mean"] if s_off is not None else None

    fig, ax = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    ax = ax.ravel()

    # (1) Onset curl
    ax[0].plot(etas, tcurl_plot, marker="o", label="scan: t_onset(curl)")
    ax[0].axhline(Tmax_plot, linestyle="--", linewidth=1)
    ax[0].set_ylabel("t_onset (curl j)  (NaN shown as >Tmax)")
    ax[0].set_title("Onset time vs eta_peak (curl criterion)")

    # mark special
    # Add OFF as an eta=0 point for visual context (if available)
    if s_off is not None:
        ax[0].scatter([0.0], [s_off["tcurl"] if not np.isnan(s_off["tcurl"]) else NO_ONSET_Y],
                      marker="s", s=80, label="OFF")
        ax[1].scatter([0.0], [s_off["tdiv"] if not np.isnan(s_off["tdiv"]) else NO_ONSET_Y],
                      marker="s", s=80, label="OFF")
    # Mark ON_center / ON_displaced late means on the bottom panels
    if s_c is not None:
        ax[2].scatter([ep_c], [s_c["jc_late_mean"]], marker="^", s=80, label="ON_center")
        ax[3].scatter([ep_c], [s_c["jd_late_mean"]], marker="^", s=80, label="ON_center")

    if s_d is not None:
        ax[2].scatter([ep_d], [s_d["jc_late_mean"]], marker="x", s=80, label="ON_displaced")
        ax[3].scatter([ep_d], [s_d["jd_late_mean"]], marker="x", s=80, label="ON_displaced")

    ax[2].legend()
    ax[3].legend()

    # (2) Onset div
    ax[1].plot(etas, tdiv_plot, marker="o", label="scan: t_onset(div)")
    ax[1].axhline(Tmax_plot, linestyle="--", linewidth=1)
    ax[1].set_ylabel("t_onset (div j)  (NaN shown as >Tmax)")
    ax[1].set_title("Onset time vs eta_peak (div criterion)")
    ax[1].legend()

    # (3) Late mean curl (log scale helps)
    ax[2].plot(etas, jc_late, marker="o", label="scan: late mean curl")
    ax[2].set_yscale("log")
    ax[2].set_xlabel("eta_peak")
    ax[2].set_ylabel("mean_{t>=10} (curl j)^2 (bulk)")
    ax[2].set_title("Late-time bulk activity vs eta_peak (curl)")
    if off_jc is not None and not np.isnan(off_jc):
        ax[2].axhline(off_jc, linestyle="--", linewidth=1, label="OFF late mean")
        ax[2].scatter([0.0], [off_jc], marker="s", s=80)
    ax[2].legend()

    # (4) Late mean div
    ax[3].plot(etas, jd_late, marker="o", label="scan: late mean div")
    ax[3].set_yscale("log")
    ax[3].set_xlabel("eta_peak")
    ax[3].set_ylabel("mean_{t>=10} (div j)^2 (bulk)")
    ax[3].set_title("Late-time bulk activity vs eta_peak (div)")
    if off_jd is not None and not np.isnan(off_jd):
        ax[3].axhline(off_jd, linestyle="--", linewidth=1, label="OFF late mean")
        ax[3].scatter([0.0], [off_jd], marker="s", s=80)
    ax[3].legend()

    for a in ax[:2]:
        a.set_ylim(0.0, Tmax_plot + 0.8)

    plt.tight_layout()
    fig.savefig("wp2_onset_scan_v1_summary.png", dpi=200)
    print("[Wrote] wp2_onset_scan_v1_summary.png")
    plt.show()

    # Optional: write a compact CSV summary
    csv = "wp2_onset_scan_summary_v1.csv"
    with open(csv, "w", encoding="utf-8") as f:
        f.write("case,eta_peak,yg,sigma,t_onset_curl,t_onset_div,jc_max,jd_max,jc_late_mean,jd_late_mean,Deta_max\n")
        for (name, ep, yg, sig, s) in results:
            f.write(
                f"{name},{ep:.6f},{yg:.6f},{sig:.6f},"
                f"{s['tcurl']},{s['tdiv']},"
                f"{s['jc_max']},{s['jd_max']},"
                f"{s['jc_late_mean']},{s['jd_late_mean']},"
                f"{s['Deta_max']}\n"
            )
    print(f"\n[Wrote] {csv}")
