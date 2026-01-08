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
    """
    Smooth pulse with a plateau (Tukey window).
    alpha in (0,1): fraction of the interval used for cosine tapers.
    """
    if (t < t0) or (t > t0 + tau):
        return 0.0
    s = (t - t0) / tau  # in [0,1]
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
                        n_vort=24,
                        y_sep=2.0,
                        core=0.25,
                        jitter=0.0,
                        seed=0):
    """
    Double vortex-sheet analog:
      row at y= +y_sep/2: +1 vortices
      row at y= -y_sep/2: -1 vortices
    """
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

    # normalize mean density ~ 1
    rho = xp.abs(psi)**2
    psi = psi / xp.sqrt(xp.mean(rho))
    return psi

# =========================
# Current-based diagnostics (NO 1/rho)
# =========================
def current_from_psi(psi, KX, KY, xp):
    """
    j = Im( psi* grad psi )  (no division by rho)
    Returns (jx, jy, rho, dpsi_dx, dpsi_dy)
    """
    dpsi_dx, dpsi_dy = grad_spectral(psi, KX, KY, xp)
    rho = xp.abs(psi)**2
    jx = xp.imag(xp.conj(psi) * dpsi_dx)
    jy = xp.imag(xp.conj(psi) * dpsi_dy)
    return jx, jy, rho, dpsi_dx, dpsi_dy

def curl_and_div_of_vector(fx, fy, KX, KY, xp):
    """
    curl_z = d/dx fy - d/dy fx
    div    = d/dx fx + d/dy fy
    """
    dfy_dx, dfy_dy = grad_spectral(fy, KX, KY, xp)
    dfx_dx, dfx_dy = grad_spectral(fx, KX, KY, xp)
    curl_z = dfy_dx - dfx_dy
    div = dfx_dx + dfy_dy
    return xp.real(curl_z), xp.real(div)

# =========================
# Core integrator
# =========================
def run_sim(case_name,
            *,
            Nx=256, Ny=256,
            Lx=24.0, Ly=24.0,
            dt=0.005, nsteps=2400,     # total time = dt*nsteps
            g=1.0, mu=1.0, V0=0.0,
            # sheet IC params
            n_vort=24, y_sep=2.0, core=0.25, jitter=0.0, seed=0,
            # gate params
            eta_peak=0.0, t0=3.0, tau=3.0, sigma=1.2, xg=0.0, yg=0.0,
            gate_alpha=0.25,
            DIFF_LIMIT=0.35,
            # diagnostics thresholds
            RHO_CORE=0.2,
            RHO_BULK=0.6,
            diag_stride=10,
            # output controls
            WRITE_FILES=True,
            save_prefix="wp2_out_current",
            snap_times=(0.0, 3.0, 4.0, 5.0, 6.0, 8.0, 12.0)):

    # grid
    x = xp.linspace(-0.5*Lx, 0.5*Lx, Nx, endpoint=False)
    y = xp.linspace(-0.5*Ly, 0.5*Ly, Ny, endpoint=False)
    X, Y = xp.meshgrid(x, y, indexing="ij")
    KX, KY, K2, kmax2 = make_kgrid(Nx, Ny, Lx, Ly, xp)

    # initial condition
    psi = imprint_vortex_rows(X, Y, Lx, Ly,
                              n_vort=n_vort, y_sep=y_sep, core=core,
                              jitter=jitter, seed=seed)

    # kinetic phase
    kin_phase = xp.exp(-1j * 0.5 * K2 * dt)

    # diagnostics arrays
    t_hist = []
    core_frac_hist = []
    min_rho_hist = []

    # current-based
    jcurl2_all_hist = []
    jcurl2_bulk_hist = []
    jdiv2_all_hist = []
    jdiv2_bulk_hist = []
    j2_all_hist = []
    j2_bulk_hist = []
    max_j_hist = []

    # gate work / activity
    Deta_hist = []
    eta_max_hist = []

    # snapshots (rho + |j| + curl(j))
    snaps = {}

    def record_snapshot(t, psi, eta_for_diag):
        jx, jy, rho, dpsi_dx, dpsi_dy = current_from_psi(psi, KX, KY, xp)
        curl_j, _ = curl_and_div_of_vector(jx, jy, KX, KY, xp)
        jmag = xp.sqrt(jx*jx + jy*jy)
        # store CPU copies
        snaps[float(t)] = (to_cpu(rho), to_cpu(jmag), to_cpu(curl_j), float(eta_for_diag))

    total_time = dt * nsteps
    print(f"[{case_name}] Nx={Nx} Ny={Ny} dt={dt} nsteps={nsteps} T={total_time:.3f} eta_peak={eta_peak}")

    # snap at t=0 (eta=0 at t=0)
    record_snapshot(0.0, psi, 0.0)

    t = 0.0
    for n in range(1, nsteps + 1):

        # --- Strang split: NL half ---
        rho = xp.abs(psi)**2
        psi *= xp.exp(-1j * (V0 + g * rho - mu) * (0.5 * dt))

        # --- kinetic full ---
        psi_hat = xp.fft.fft2(psi)
        psi_hat *= kin_phase
        psi = xp.fft.ifft2(psi_hat)

        # --- NL half ---
        rho = xp.abs(psi)**2
        psi *= xp.exp(-1j * (V0 + g * rho - mu) * (0.5 * dt))

        # --- gated diffusion ---
        eta_max_val = 0.0
        eta_for_diag = 0.0
        if eta_peak != 0.0:
            eta = eta_gaussian(X, Y, Lx, Ly, eta_peak, xg, yg, sigma, t, t0, tau, alpha=gate_alpha)
            eta_max_val = float(to_cpu(xp.max(eta)))
            eta_for_diag = eta_max_val
            if eta_max_val > 0.0:
                dt_sub_max = DIFF_LIMIT / (eta_max_val * kmax2 + 1e-30)
                n_sub = int(np.ceil(dt / dt_sub_max))
                n_sub = max(1, n_sub)
                dt_sub = dt / n_sub
                for _ in range(n_sub):
                    dpsi_dx, dpsi_dy = grad_spectral(psi, KX, KY, xp)
                    fx = eta * dpsi_dx
                    fy = eta * dpsi_dy
                    div_term = div_spectral(fx, fy, KX, KY, xp)
                    psi = psi + dt_sub * div_term

        # time
        t = n * dt

        # snapshots near requested times
        for ts in snap_times:
            if (ts not in snaps) and (abs(t - ts) <= 0.5 * dt):
                record_snapshot(ts, psi, eta_for_diag)

        # diagnostics
        if (n % diag_stride) == 0:
            jx, jy, rho, dpsi_dx, dpsi_dy = current_from_psi(psi, KX, KY, xp)
            curl_j, div_j = curl_and_div_of_vector(jx, jy, KX, KY, xp)

            j2 = (jx*jx + jy*jy).astype(xp.float64)
            jcurl2 = (curl_j*curl_j).astype(xp.float64)
            jdiv2 = (div_j*div_j).astype(xp.float64)

            mask_bulk = (rho > RHO_BULK).astype(xp.float64)

            core_frac = xp.mean((rho < RHO_CORE).astype(xp.float64))
            min_rho = xp.min(rho)
            max_j = xp.max(xp.sqrt(j2))

            # gate dissipation power proxy: < eta |grad psi|^2 >
            # (uses eta_max_val to decide whether eta is active; if inactive, treat as 0)
            if eta_peak != 0.0 and eta_max_val > 0.0:
                # reuse eta at this step (already exists)
                grad2 = (xp.abs(dpsi_dx)**2 + xp.abs(dpsi_dy)**2).astype(xp.float64)
                Deta = xp.mean((eta.astype(xp.float64) * grad2))
            else:
                Deta = xp.array(0.0, dtype=xp.float64)

            # store
            t_hist.append(float(t))
            core_frac_hist.append(float(to_cpu(core_frac)))
            min_rho_hist.append(float(to_cpu(min_rho)))

            jcurl2_all_hist.append(float(to_cpu(xp.mean(jcurl2))))
            jcurl2_bulk_hist.append(float(to_cpu(masked_mean(jcurl2, mask_bulk))))
            jdiv2_all_hist.append(float(to_cpu(xp.mean(jdiv2))))
            jdiv2_bulk_hist.append(float(to_cpu(masked_mean(jdiv2, mask_bulk))))

            j2_all_hist.append(float(to_cpu(xp.mean(j2))))
            j2_bulk_hist.append(float(to_cpu(masked_mean(j2, mask_bulk))))
            max_j_hist.append(float(to_cpu(max_j)))

            Deta_hist.append(float(to_cpu(Deta)))
            eta_max_hist.append(float(abs(eta_max_val)))

    # final snapshot
    if float(total_time) not in snaps:
        record_snapshot(float(total_time), psi, eta_for_diag)

    out = {
        "case": case_name,
        "t": np.array(t_hist),
        "core_frac": np.array(core_frac_hist),
        "min_rho": np.array(min_rho_hist),
        "jcurl2_all": np.array(jcurl2_all_hist),
        "jcurl2_bulk": np.array(jcurl2_bulk_hist),
        "jdiv2_all": np.array(jdiv2_all_hist),
        "jdiv2_bulk": np.array(jdiv2_bulk_hist),
        "j2_all": np.array(j2_all_hist),
        "j2_bulk": np.array(j2_bulk_hist),
        "max_j": np.array(max_j_hist),
        "Deta": np.array(Deta_hist),
        "eta_max": np.array(eta_max_hist),
        "snaps": snaps,
        "params": dict(
            Nx=Nx, Ny=Ny, Lx=Lx, Ly=Ly, dt=dt, nsteps=nsteps,
            g=g, mu=mu, V0=V0,
            n_vort=n_vort, y_sep=y_sep, core=core, jitter=jitter, seed=seed,
            eta_peak=eta_peak, t0=t0, tau=tau, sigma=sigma, xg=xg, yg=yg,
            gate_alpha=gate_alpha, DIFF_LIMIT=DIFF_LIMIT,
            RHO_CORE=RHO_CORE, RHO_BULK=RHO_BULK, diag_stride=diag_stride
        )
    }

    if WRITE_FILES:
        csv_path = f"{save_prefix}_{case_name}_diag.csv"
        header = (
            "t,eta_max,Deta,core_frac,min_rho,max_j,"
            "j2_all,j2_bulk,jcurl2_all,jcurl2_bulk,jdiv2_all,jdiv2_bulk\n"
        )
        data = np.column_stack([
            out["t"], out["eta_max"], out["Deta"],
            out["core_frac"], out["min_rho"], out["max_j"],
            out["j2_all"], out["j2_bulk"],
            out["jcurl2_all"], out["jcurl2_bulk"],
            out["jdiv2_all"], out["jdiv2_bulk"],
        ])
        with open(csv_path, "w", encoding="utf-8") as f:
            f.write(header)
            np.savetxt(f, data, delimiter=",", fmt="%.8e")
        print(f"[{case_name}] Wrote {csv_path}")

    return out

def spike_report(out, label, field, k=5):
    t = out["t"]
    y = out[field]
    if len(y) == 0:
        print(f"[{label}] No data for {field}")
        return
    idx = np.argsort(y)[::-1][:k]
    print(f"\n[{label}] Top {k} spikes in {field}:")
    for i in idx:
        print(
            f"  t={t[i]:.4f}  {field}={y[i]:.6e}  "
            f"Deta={out['Deta'][i]:.3e}  eta_max={out['eta_max'][i]:.3e}  "
            f"min_rho={out['min_rho'][i]:.3e}  max_j={out['max_j'][i]:.3e}"
        )

# =========================
# Main: OFF vs ON
# =========================
if __name__ == "__main__":

    # Heavier defaults if GPU present
    Nx_default = 512 if use_gpu else 256
    Ny_default = Nx_default

    common = dict(
        Nx=Nx_default, Ny=Ny_default,
        Lx=24.0, Ly=24.0,
        dt=0.005,
        nsteps=2400,          # total time = 12
        g=1.0, mu=1.0, V0=0.0,
        n_vort=24,
        y_sep=2.0,
        core=0.25,
        jitter=0.00,
        seed=1,
        RHO_CORE=0.2,
        RHO_BULK=0.6,
        diag_stride=10,
        WRITE_FILES=True,
        save_prefix="wp2_current",
        snap_times=(0.0, 3.0, 4.0, 5.0, 6.0, 8.0, 12.0),
    )

    gate = dict(
        eta_peak=0.08,
        t0=3.0,
        tau=3.0,
        sigma=1.2,
        xg=0.0,
        yg=0.0,
        gate_alpha=0.25,
        DIFF_LIMIT=0.35,
    )

    print("\n[Run] Case A: Gate OFF")
    A = run_sim("gate_off",
                eta_peak=0.0, t0=gate["t0"], tau=gate["tau"], sigma=gate["sigma"],
                xg=gate["xg"], yg=gate["yg"], gate_alpha=gate["gate_alpha"],
                DIFF_LIMIT=gate["DIFF_LIMIT"],
                **common)

    print("\n[Run] Case B: Gate ON")
    B = run_sim("gate_on", **gate, **common)

    # Spike reports (these are the ones we care about now)
    spike_report(A, "OFF", "jcurl2_all", k=5)
    spike_report(B, "ON ", "jcurl2_all", k=5)
    spike_report(A, "OFF", "jcurl2_bulk", k=5)
    spike_report(B, "ON ", "jcurl2_bulk", k=5)

    spike_report(A, "OFF", "jdiv2_all", k=5)
    spike_report(B, "ON ", "jdiv2_all", k=5)
    spike_report(A, "OFF", "jdiv2_bulk", k=5)
    spike_report(B, "ON ", "jdiv2_bulk", k=5)

    spike_report(A, "OFF", "Deta", k=5)
    spike_report(B, "ON ", "Deta", k=5)

    # =========================
    # One figure with key curves
    # =========================
    fig, ax = plt.subplots(5, 1, figsize=(9, 14), sharex=True)

    ax[0].plot(A["t"], A["core_frac"], label="core frac (off)")
    ax[0].plot(B["t"], B["core_frac"], label="core frac (on)")
    ax[0].set_ylabel("core fraction (rho<RHO_CORE)")
    ax[0].legend()

    ax[1].plot(A["t"], A["jcurl2_all"], label="⟨(curl j)^2⟩ all (off)")
    ax[1].plot(B["t"], B["jcurl2_all"], label="⟨(curl j)^2⟩ all (on)")
    ax[1].plot(A["t"], A["jcurl2_bulk"], "--", label="⟨(curl j)^2⟩ bulk (off)")
    ax[1].plot(B["t"], B["jcurl2_bulk"], "--", label="⟨(curl j)^2⟩ bulk (on)")
    ax[1].set_ylabel("current-enstrophy proxy")
    ax[1].legend()

    ax[2].plot(A["t"], A["jdiv2_all"], label="⟨(div j)^2⟩ all (off)")
    ax[2].plot(B["t"], B["jdiv2_all"], label="⟨(div j)^2⟩ all (on)")
    ax[2].plot(A["t"], A["jdiv2_bulk"], "--", label="⟨(div j)^2⟩ bulk (off)")
    ax[2].plot(B["t"], B["jdiv2_bulk"], "--", label="⟨(div j)^2⟩ bulk (on)")
    ax[2].set_ylabel("current-compressibility proxy")
    ax[2].legend()

    ax[3].plot(A["t"], A["Deta"], label="D_eta (off)")
    ax[3].plot(B["t"], B["Deta"], label="D_eta (on)")
    ax[3].set_ylabel("gate work  ⟨η|∇ψ|²⟩")
    ax[3].legend()

    ax[4].plot(A["t"], A["min_rho"], label="min rho (off)")
    ax[4].plot(B["t"], B["min_rho"], label="min rho (on)")
    ax[4].plot(B["t"], B["eta_max"], "--", label="eta_max (on)")
    ax[4].set_ylabel("min rho / eta_max")
    ax[4].set_xlabel("time")
    ax[4].legend()

    plt.tight_layout()
    plt.show()

    print("\n[Done] Look at jcurl2_bulk and jdiv2_bulk: if ON shows a clear gate-window bump there, we can start onset scans.")
