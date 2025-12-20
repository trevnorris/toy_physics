#!/usr/bin/env python3
"""
quantum_vortex_reconnection_wp1_diagnostics.py

WP1: Topology-switch "gate" + diagnostics (sound + core topology proxy)

This is an updated version of the project's 3D GPE script that adds:
  1) A pulsed, localized diffusion ("eta gate") that can open topology change.
  2) A stability-safe implementation (automatic diffusion subcycling).
  3) Extra diagnostics:
       - core volume fraction (rho < RHO_TH)
       - number of large connected core components (n_big)
  4) Plots:
       - 3D core scatter (topology snapshot)
       - Compressible energy (sound) with gate overlay
       - Core fraction + component count with gate overlay

Notes:
- Works on CPU (NumPy) and optionally GPU (CuPy).
- Component counting uses scipy.ndimage.label; if SciPy isn't available,
  the script will still run and will plot core fraction only.

Run:
  python quantum_vortex_reconnection_wp1_diagnostics.py
"""

import time
import numpy as np
import matplotlib.pyplot as plt

# --- Optional SciPy for connected components ---
try:
    from scipy.ndimage import label as nd_label, gaussian_filter
    _HAS_SCIPY = True
except Exception:
    nd_label = None
    _HAS_SCIPY = False

# --- GPU SETUP ---
try:
    import cupy as cp
    xp = cp
    use_gpu = True
    print("[System] CuPy detected. GPU Mode Engaged.")
except Exception:
    xp = np
    use_gpu = False
    print("[System] CuPy not found. CPU Mode (Slower).")

# -----------------------
# CONFIGURATION (edit me)
# -----------------------
N = 128              # Grid size
L = 16.0             # Domain size
DT = 0.05            # Real-time step
STEPS_RELAX = 100    # Imaginary-time relaxation steps (quiet start)
STEPS_REAL = 400     # Real-time steps
PLOT_EVERY = 5       # Plot cadence (smaller = better timing resolution)

# --- Gate ("switch slot") parameters ---
# eta(x,t) = D_ETA * Wt(t) * Wx(x), where Wt is a smooth-ish trapezoid and Wx is a Gaussian.
D_ETA = 0.08         # peak eta inside gate (units L^2 / T)
PULSE_T0 = 8.0
PULSE_TAU = 8.0
PULSE_SIGMA = 2.0    # spatial localization width
PULSE_CENTER = (0.0, 0.0, 0.0)  # (x0,y0,z0)

# --- Diffusion stability ---
# Explicit diffusion is stable only if dt * eta * k_max^2 is sufficiently small.
# DIFF_LIMIT is the max allowed value for (dt_sub * eta_peak * k2_max).
DIFF_LIMIT = 0.4     # smaller = more stable, slower (0.2–0.6 are reasonable)

# --- Core/topology proxy settings ---
RHO_TH = 0.05        # core threshold for topology proxy
MIN_VOX = 20000      # minimum voxel count to call a component "big"
RHO_VIS = 0.02       # visualization threshold for scatter (rho < RHO_VIS)

# --- Plot downsample (scatter can be heavy) ---
SCATTER_STRIDE = 4   # bigger = fewer points


# -----------------------
# GRID + FOURIER SPACE
# -----------------------
x = xp.linspace(-L/2, L/2, N)
y = xp.linspace(-L/2, L/2, N)
z = xp.linspace(-L/2, L/2, N)
X, Y, Z = xp.meshgrid(x, y, z, indexing="ij")
if use_gpu:
    X_cpu = X.get()
    Y_cpu = Y.get()
    Z_cpu = Z.get()
else:
    X_cpu = X
    Y_cpu = Y
    Z_cpu = Z

dx = float((x[1] - x[0]).get() if use_gpu else (x[1] - x[0]))
k = xp.fft.fftfreq(N, d=dx/(2*xp.pi)) * 2*xp.pi
KX, KY, KZ = xp.meshgrid(k, k, k, indexing="ij")
K2 = KX**2 + KY**2 + KZ**2

# --- DEALIASING MASK (2/3 Rule) ---
k_max = float((k.max()).get() if use_gpu else k.max())
dealias_mask = (xp.abs(KX) < 2/3 * k_max) & (xp.abs(KY) < 2/3 * k_max) & (xp.abs(KZ) < 2/3 * k_max)

# Use the max k^2 inside the dealiased region for stability estimates
K2_eff_max = float((xp.max(K2 * dealias_mask)).get() if use_gpu else xp.max(K2 * dealias_mask))


# -----------------------
# INITIAL CONDITION: ORTHOGONAL VORTEX TUBES
# -----------------------
def vortex_tube_phase(X_grid, Y_grid, x0, y0):
    return xp.arctan2(Y_grid - y0, X_grid - x0)

def vortex_tube_density(X_grid, Y_grid, x0, y0):
    r = xp.sqrt((X_grid - x0)**2 + (Y_grid - y0)**2)
    return (r**2) / (r**2 + 1.0)  # simple Padé-ish core

# Tube 1 (along X): phase winding in (Y,Z)
phi1 = vortex_tube_phase(Y, Z, 0.0, -1.2)
rho1 = vortex_tube_density(Y, Z, 0.0, -1.2)

# Tube 2 (along Y): phase winding in (X,Z)
phi2 = vortex_tube_phase(X, Z, 0.0,  1.2)
rho2 = vortex_tube_density(X, Z, 0.0,  1.2)

psi_phase = phi1 + phi2
psi_rho = rho1 * rho2
psi = xp.sqrt(psi_rho) * xp.exp(1j * psi_phase)

# Save analytic phase for pinning during relaxation
initial_phase = xp.angle(psi)


# -----------------------
# DIAGNOSTICS
# -----------------------
def gradients_of_psi(psi_curr):
    """Return spectral gradients of psi."""
    psi_k = xp.fft.fftn(psi_curr)
    gx = xp.fft.ifftn(1j * KX * psi_k)
    gy = xp.fft.ifftn(1j * KY * psi_k)
    gz = xp.fft.ifftn(1j * KZ * psi_k)
    return gx, gy, gz

def mass_current_and_velocity(psi_curr, eps=1e-12):
    """
    j = Im(conj(psi) * grad psi)
    v ≈ j / rho   (superfluid velocity proxy)
    """
    gx, gy, gz = gradients_of_psi(psi_curr)
    rho = xp.abs(psi_curr)**2
    jx = xp.imag(xp.conj(psi_curr) * gx)
    jy = xp.imag(xp.conj(psi_curr) * gy)
    jz = xp.imag(xp.conj(psi_curr) * gz)
    vmag2 = (jx*jx + jy*jy + jz*jz) / (rho*rho + eps)
    return (jx, jy, jz), rho, vmag2

def get_compressible_energy(psi_curr):
    """
    Helmholtz-like projection:
      - compute current j
      - project j_k onto k (compressible part)
      - Ec ~ 0.5 * ∫ |j_c|^2
    """
    (jx, jy, jz), _, _ = mass_current_and_velocity(psi_curr)

    jk_x = xp.fft.fftn(jx)
    jk_y = xp.fft.fftn(jy)
    jk_z = xp.fft.fftn(jz)

    k_dot_j = KX*jk_x + KY*jk_y + KZ*jk_z

    K2_safe = K2.copy()
    K2_safe[0, 0, 0] = 1.0

    jc_x = KX * (k_dot_j / K2_safe)
    jc_y = KY * (k_dot_j / K2_safe)
    jc_z = KZ * (k_dot_j / K2_safe)
    jc_x[0, 0, 0] = 0.0
    jc_y[0, 0, 0] = 0.0
    jc_z[0, 0, 0] = 0.0

    Ec = 0.5 * xp.sum(xp.abs(jc_x)**2 + xp.abs(jc_y)**2 + xp.abs(jc_z)**2)
    Ec = float(Ec.get() if use_gpu else Ec)
    return Ec / (N**3)

def core_metrics_cpu(rho_cpu, rho_th=0.05, min_vox=20000):
    """
    CPU-only: compute core fraction and number of large connected core components.
    """
    if not _HAS_SCIPY:
        core = rho_cpu < rho_th
        core_frac = float(core.mean())
        return core_frac, None

    rho_s = gaussian_filter(rho_cpu, sigma=1.0)
    core = rho_s < rho_th
    core_frac = float(core.mean())

    structure = np.ones((3, 3, 3), dtype=np.int8)  # 26-connectivity
    lab, n = nd_label(core, structure=structure)
    if n == 0:
        return core_frac, 0

    sizes = np.bincount(lab.ravel())
    # sizes[0] is background
    big = (sizes[1:] >= int(min_vox))
    n_big = int(big.sum())
    return core_frac, n_big


# -----------------------
# GATE: eta(x,t) and diffusion step
# -----------------------
def gate_profile_time(t, t0, tau, ramp=0.15):
    """
    Smooth trapezoid:
      - ramps up over ramp*tau
      - holds
      - ramps down over ramp*tau
    Returns scalar in [0,1].
    """
    if tau <= 0:
        return 0.0
    t1 = t0
    t2 = t0 + ramp * tau
    t3 = t0 + (1 - ramp) * tau
    t4 = t0 + tau

    if t < t1 or t > t4:
        return 0.0
    if t2 <= t <= t3:
        return 1.0
    # ramp up
    if t1 <= t < t2:
        s = (t - t1) / max(t2 - t1, 1e-12)
        return 0.5 - 0.5*np.cos(np.pi*s)
    # ramp down
    if t3 < t <= t4:
        s = (t - t3) / max(t4 - t3, 1e-12)
        return 0.5 + 0.5*np.cos(np.pi*s)
    return 0.0

# Precompute spatial window Wx on the simulation grid
x0, y0, z0 = PULSE_CENTER
Wx = xp.exp(-((X - x0)**2 + (Y - y0)**2 + (Z - z0)**2) / (2 * (PULSE_SIGMA**2)))

def laplacian_field(field):
    """Spectral Laplacian: ∇^2 f = ifft(-k^2 * fft(f))"""
    fk = xp.fft.fftn(field)
    fk *= (-K2)
    fk *= dealias_mask
    return xp.fft.ifftn(fk)

def apply_diffusion_gate(psi_curr, eta_field, dt, eta_peak):
    """
    Apply explicit diffusion: psi <- psi + dt * eta(x)*∇^2 psi
    Uses subcycling for stability.

    Stability proxy:
      dt_sub * eta_peak * K2_eff_max <= DIFF_LIMIT
    """
    if eta_peak <= 0.0:
        return psi_curr, 0  # no substeps

    # choose number of substeps
    dt_sub_max = DIFF_LIMIT / (eta_peak * K2_eff_max + 1e-30)
    n_sub = int(np.ceil(dt / max(dt_sub_max, 1e-30)))
    n_sub = max(n_sub, 1)
    dt_sub = dt / n_sub

    psi_out = psi_curr
    for _ in range(n_sub):
        lap = laplacian_field(psi_out)
        psi_out = psi_out + (dt_sub * eta_field) * lap
    return psi_out, n_sub


# -----------------------
# HISTORY + PLOTTING
# -----------------------
history_t = []
history_gate = []
history_Ec = []
history_core_frac = []
history_nbig = []

# Scalars preview lines
print("[WP1 Scalars Preview]")
print("  Canonical IC: orthogonal vortex tubes")
print(f"  Gate timing: t0={PULSE_T0}, tau={PULSE_TAU}")
print(f"  Gate peak eta ≈ {D_ETA:.4f} (units L^2/T)")

# ld, chi (delta is a user-defined thickness scale; we print with delta=1.0)
delta_assumed = 1.0
ld = np.sqrt(max(D_ETA, 0.0) * max(PULSE_TAU, 0.0))
chi = ld / delta_assumed if delta_assumed > 0 else float("nan")
print(f"  ld = sqrt(eta_peak*tau) ≈ {ld:.4f}")
print(f"  χ = ld/δ with δ=1.0 ≈ {chi:.4f}")

if not _HAS_SCIPY:
    print("[Note] SciPy not found; will plot core fraction but skip component counting (n_big).")

plt.ion()
fig = plt.figure(figsize=(11, 9))
plt.subplots_adjust(hspace=0.35)
ax_top = fig.add_subplot(311, projection="3d")
ax_mid = fig.add_subplot(312)
ax_bot = fig.add_subplot(313)

_u0_printed = False  # print U0/Rm once, near first time the gate is active


# -----------------------
# PHASE 1: QUIET START (Imaginary time relaxation with phase pinning)
# -----------------------
print("--- PHASE 1: QUIET START (Relaxation) ---")
print("Damping initial phonon shock while pinning vortices...")

for step in range(STEPS_RELAX):
    rho = xp.abs(psi)**2
    psi *= xp.exp(-rho * DT)  # imaginary-time "potential" decay

    # phase pinning
    amp = xp.abs(psi)
    psi = amp * xp.exp(1j * initial_phase)

    # imaginary-time kinetic decay
    psi_k = xp.fft.fftn(psi)
    psi_k *= xp.exp(-0.5 * K2 * DT)
    psi_k *= dealias_mask
    psi = xp.fft.ifftn(psi_k)

    # normalize (keep max density ~ 1)
    norm_factor = 1.0 / xp.max(xp.abs(psi))
    psi *= norm_factor

print("Relaxation Complete. System is 'Quiet'.")
print("--- PHASE 2: UNITARY EVOLUTION (The Experiment) ---")


# -----------------------
# PHASE 2: REAL TIME EVOLUTION
# -----------------------
t_start = time.time()

for step in range(STEPS_REAL):
    t = step * DT

    # --- Unitary split-step (GPE-like) ---
    rho = xp.abs(psi)**2
    psi *= xp.exp(-1j * rho * DT)  # nonlinear potential

    psi_k = xp.fft.fftn(psi)
    psi_k *= xp.exp(-0.5j * K2 * DT)  # kinetic
    psi_k *= dealias_mask
    psi = xp.fft.ifftn(psi_k)

    # --- Gate: diffusion "switch slot" ---
    g = gate_profile_time(t, PULSE_T0, PULSE_TAU)
    eta_peak = float(D_ETA * g)
    gate_on = (eta_peak > 0.0)

    if gate_on:
        eta_field = (D_ETA * g) * Wx
        psi, n_sub = apply_diffusion_gate(psi, eta_field, DT, eta_peak)
    else:
        n_sub = 0

    # --- Plot / diagnostics cadence ---
    if step % PLOT_EVERY == 0:
        # Sound metric
        Ec = get_compressible_energy(psi)

        # Core metrics (computed on CPU for robustness)
        rho_now = xp.abs(psi)**2
        rho_cpu = rho_now.get() if use_gpu else np.asarray(rho_now)

        core_frac, n_big = core_metrics_cpu(rho_cpu, rho_th=RHO_TH, min_vox=MIN_VOX)

        history_t.append(t)
        history_gate.append(float(g))
        history_Ec.append(Ec)
        history_core_frac.append(core_frac)
        history_nbig.append(-1 if n_big is None else int(n_big))

        # Print U0/Rm once when the gate is first active and we have a snapshot
        if (not _u0_printed) and gate_on:
            # Estimate RMS velocity in the gate region, weighted by Wx
            (_, _, _), _, vmag2 = mass_current_and_velocity(psi)
            # weighted average in the gate region (Wx as weight)
            w = Wx
            num = xp.sum(w * vmag2)
            den = xp.sum(w) + 1e-30
            u0 = float((xp.sqrt(num / den)).get() if use_gpu else xp.sqrt(num / den))

            # Rm ~ U0*delta/eta_peak, with delta=1
            rm = (u0 * delta_assumed) / max(eta_peak, 1e-30)

            print(f"  U0 (RMS in gate) ≈ {u0:.4f}")
            print(f"  Rm ~ U0*δ/eta_peak with δ=1.0: Rm ≈ {rm:.0f}")
            _u0_printed = True

        # --- Plotting ---
        ax_top.clear()
        mask = rho_cpu < RHO_VIS
        xs = X_cpu[mask][::SCATTER_STRIDE]
        ys = Y_cpu[mask][::SCATTER_STRIDE]
        zs = Z_cpu[mask][::SCATTER_STRIDE]
        ax_top.scatter(xs, ys, zs, c="k", s=1, alpha=0.35)
        ax_top.set_xlim(-L/2, L/2); ax_top.set_ylim(-L/2, L/2); ax_top.set_zlim(-L/2, L/2)
        ax_top.set_title(f"Step {step}   t={t:.2f}   Gate={'ON' if gate_on else 'OFF'}   substeps={n_sub}")

        ax_mid.clear()
        ax_mid.plot(history_t, history_Ec, linewidth=2)
        # Scale gate overlay to ~20% of max energy for visibility
        if len(history_Ec) > 0:
            scale = 0.2 * max(history_Ec)
        else:
            scale = 1.0
        ax_mid.plot(history_t, [scale * gg for gg in history_gate], "k--", linewidth=1.5)
        ax_mid.set_xlabel("Time")
        ax_mid.set_ylabel("Compressible Energy (Sound)")
        ax_mid.set_title("Acoustic Emission (look for change during/after gate)")
        ax_mid.grid(True)
        ax_mid.legend(["Compressible Energy", "Gate (scaled)"], loc="upper right")

        ax_bot.clear()
        ax_bot.plot(history_t, history_core_frac, linewidth=2)
        ax_bot.set_xlabel("Time")
        ax_bot.set_ylabel("Core fraction (rho < threshold)")
        ax_bot.grid(True)

        # Plot n_big if available
        if _HAS_SCIPY:
            ax2 = ax_bot.twinx()
            # Replace -1 with NaN (shouldn't happen if SciPy is present)
            n_big_series = [np.nan if v < 0 else v for v in history_nbig]
            ax2.plot(history_t, n_big_series, linestyle="--", linewidth=2)
            ax2.set_ylabel("n_big (connected core components)")
            ax_bot.legend(["core_frac"], loc="upper left")
            ax2.legend(["n_big"], loc="upper right")
        else:
            ax_bot.legend(["core_frac"], loc="upper left")

        # Gate overlay on bottom plot too
        g_scale2 = 0.25 * max(history_core_frac) if len(history_core_frac) else 0.1
        ax_bot.plot(history_t, [g_scale2 * gg for gg in history_gate], "k--", linewidth=1.2)

        plt.pause(0.001)

t_end = time.time()
print(f"[Done] Runtime: {t_end - t_start:.2f} s")

plt.ioff()
plt.show()
