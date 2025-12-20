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

# --- CONFIGURATION ---
N = 192
L = 16.0
DT = 0.05
STEPS_RELAX = 100
STEPS_REAL = 400
SAMPLE_EVERY = 10
W_CORE = 1.0
RHO_CORE_TH = 0.2
N_SHELLS = 40
DO_PLOT = True  # set True if you still want live plots

# --- GRID (periodic) ---
dx = L / N
x = xp.linspace(-L/2, L/2, N, endpoint=False)
y = xp.linspace(-L/2, L/2, N, endpoint=False)
z = xp.linspace(-L/2, L/2, N, endpoint=False)
X, Y, Z = xp.meshgrid(x, y, z, indexing='ij')

# --- WAVENUMBERS (correct 2π) ---
k = (2*np.pi) * xp.fft.fftfreq(N, d=dx)
KX, KY, KZ = xp.meshgrid(k, k, k, indexing='ij')
K2 = KX**2 + KY**2 + KZ**2
K2_safe = K2.copy()
K2_safe[0, 0, 0] = 1.0

# --- DEALIASING MASK (2/3 rule based on Nyquist) ---
k_ny = np.pi / dx
k_cut = (2.0/3.0) * k_ny
dealias_mask = (xp.abs(KX) <= k_cut) & (xp.abs(KY) <= k_cut) & (xp.abs(KZ) <= k_cut)

# --- RADIAL SHELLS ---
R_GRID = xp.sqrt(X**2 + Y**2 + Z**2)
R_BINS = np.linspace(0, L/2, N_SHELLS + 1)
R_BIN_CENTERS = 0.5 * (R_BINS[1:] + R_BINS[:-1])
R_BINS_XP = xp.asarray(R_BINS)

# --- INITIAL CONDITION: ORTHOGONAL TUBES ---
def vortex_phase(u, v, u0, v0):
    return xp.arctan2(v - v0, u - u0)

def vortex_density(u, v, u0, v0):
    r = xp.sqrt((u - u0)**2 + (v - v0)**2)
    return (r**2) / (r**2 + W_CORE**2)

# Tube along X => phase depends on (y,z)
phi1 = vortex_phase(Y, Z, 0.0, -1.2)
rho1 = vortex_density(Y, Z, 0.0, -1.2)

# Tube along Y => phase depends on (x,z)
phi2 = vortex_phase(X, Z, 0.0,  1.2)
rho2 = vortex_density(X, Z, 0.0,  1.2)

psi = xp.sqrt(rho1 * rho2) * xp.exp(1j * (phi1 + phi2))
initial_phase = xp.angle(psi)

# --- DIAGNOSTICS ---
def helmholtz_energies_and_core(psi_curr):
    """
    Returns:
      Ec, Ei : compressible/incompressible kinetic energies (not densities)
      mass   : ∫|psi|^2 dV
      core_cells, core_vol_phys
    """
    rho = xp.abs(psi_curr)**2

    psi_k = xp.fft.fftn(psi_curr)
    grad_x = xp.fft.ifftn(1j * KX * psi_k)
    grad_y = xp.fft.ifftn(1j * KY * psi_k)
    grad_z = xp.fft.ifftn(1j * KZ * psi_k)

    jx = xp.imag(xp.conj(psi_curr) * grad_x)
    jy = xp.imag(xp.conj(psi_curr) * grad_y)
    jz = xp.imag(xp.conj(psi_curr) * grad_z)

    # Project j -> compressible component in k-space:
    jkx = xp.fft.fftn(jx) * dealias_mask
    jky = xp.fft.fftn(jy) * dealias_mask
    jkz = xp.fft.fftn(jz) * dealias_mask

    k_dot_j = KX*jkx + KY*jky + KZ*jkz
    jc_kx = KX * (k_dot_j / K2_safe)
    jc_ky = KY * (k_dot_j / K2_safe)
    jc_kz = KZ * (k_dot_j / K2_safe)
    jc_kx[0,0,0] = 0; jc_ky[0,0,0] = 0; jc_kz[0,0,0] = 0

    # Back to real space
    jc_x = xp.fft.ifftn(jc_kx)
    jc_y = xp.fft.ifftn(jc_ky)
    jc_z = xp.fft.ifftn(jc_kz)

    # Incompressible part
    ji_x = jx - jc_x
    ji_y = jy - jc_y
    ji_z = jz - jc_z

    # Energies (real-space integral)
    Ec = 0.5 * dx**3 * xp.sum((xp.abs(jc_x)**2 + xp.abs(jc_y)**2 + xp.abs(jc_z)**2).real)
    Ei = 0.5 * dx**3 * xp.sum((xp.abs(ji_x)**2 + xp.abs(ji_y)**2 + xp.abs(ji_z)**2).real)

    mass = dx**3 * xp.sum(rho)

    core_cells = xp.sum(rho < RHO_CORE_TH)
    core_vol_phys = core_cells * dx**3

    return float(to_cpu(Ec)), float(to_cpu(Ei)), float(to_cpu(mass)), int(to_cpu(core_cells)), float(to_cpu(core_vol_phys))

def shell_average_delta_rho_sq(psi_curr):
    rho = xp.abs(psi_curr)**2
    delta = rho - 1.0
    val = (delta*delta).real

    r_flat = R_GRID.ravel()
    v_flat = val.ravel()

    hist_vals, _ = xp.histogram(r_flat, bins=R_BINS_XP, weights=v_flat)
    hist_cnts, _ = xp.histogram(r_flat, bins=R_BINS_XP)

    hist_cnts = xp.maximum(hist_cnts, 1)
    shell_avg = hist_vals / hist_cnts
    return to_cpu(shell_avg)

# --- QUIET START: IMPROVED ITP (still pinned) ---
# Use a more standard imaginary-time nonlinear factor with (1 - rho)
print("--- PHASE 1: RELAXATION (Quiet Start) ---")
for _ in range(STEPS_RELAX):
    rho = xp.abs(psi)**2
    psi *= xp.exp((1.0 - rho) * (DT/2.0))                 # nonlinear (imag time)
    psi = xp.abs(psi) * xp.exp(1j * initial_phase)        # pin phase
    psi_k = xp.fft.fftn(psi) * xp.exp(-0.5 * K2 * DT) * dealias_mask
    psi = xp.fft.ifftn(psi_k)
    rho = xp.abs(psi)**2
    psi *= xp.exp((1.0 - rho) * (DT/2.0))
    # normalize to mean density ~1: ∫|psi|^2 dV = L^3
    mass = dx**3 * xp.sum(rho)
    psi *= xp.sqrt((L**3) / mass)

# --- STORAGE ---
t_hist, Ec_hist, Ei_hist, mass_hist = [], [], [], []
core_cells_hist, core_vol_hist = [], []
Shell_Map = []

print("--- PHASE 2: UNITARY EVOLUTION (Measurement) ---")
if DO_PLOT:
    plt.ion()
    fig = plt.figure(figsize=(12, 10))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(212)

for step in range(STEPS_REAL):
    # Strang split (unitary)
    rho = xp.abs(psi)**2
    psi *= xp.exp(-1j * (rho - 1.0) * (DT/2.0))  # subtract 1 to remove trivial global phase
    psi_k = xp.fft.fftn(psi) * xp.exp(-0.5j * K2 * DT) * dealias_mask
    psi = xp.fft.ifftn(psi_k)
    rho = xp.abs(psi)**2
    psi *= xp.exp(-1j * (rho - 1.0) * (DT/2.0))

    if step % SAMPLE_EVERY == 0:
        t = step * DT
        shell_avg = shell_average_delta_rho_sq(psi)
        Ec, Ei, mass, core_cells, core_vol = helmholtz_energies_and_core(psi)

        t_hist.append(t)
        Ec_hist.append(Ec)
        Ei_hist.append(Ei)
        mass_hist.append(mass)
        core_cells_hist.append(core_cells)
        core_vol_hist.append(core_vol)
        Shell_Map.append(shell_avg)

        # concise text log every sample
        if len(t_hist) == 1:
            Ec0, Ei0, m0 = Ec, Ei, mass
        dEc = Ec - Ec0
        dEi = Ei - Ei0
        dm  = mass - m0
        print(f"[t={t:6.2f}] Ec={Ec: .4e} Ei={Ei: .4e}  dEc={dEc:+.2e} dEi={dEi:+.2e}  mass={mass: .4e} dm={dm:+.2e}  core_cells={core_cells}")

        if DO_PLOT:
            ax1.clear()
            ax1.plot(t_hist, Ec_hist, 'r-', lw=2)
            ax1.set_title("Compressible Energy Ec")
            ax1.grid(True)

            ax2.clear()
            ax2.plot(t_hist, core_cells_hist, 'k-', lw=2)
            ax2.set_title("Core cells (<rho_th)")
            ax2.grid(True)

            ax3.clear()
            map_data = np.array(Shell_Map)
            if np.max(map_data) > 0:
                map_data = map_data / np.max(map_data)
            ax3.imshow(map_data, aspect='auto', origin='lower',
                       extent=[0, L/2, 0, t_hist[-1]])
            ax3.set_title("<(delta rho)^2> shell average (normalized)")
            ax3.set_xlabel("r")
            ax3.set_ylabel("t")
            plt.pause(0.001)

if DO_PLOT:
    plt.ioff()
    plt.show()

# --- SUMMARY REPORT (TEXT) ---
t_arr = np.array(t_hist)
Ec_arr = np.array(Ec_hist)
Ei_arr = np.array(Ei_hist)
mass_arr = np.array(mass_hist)
core_arr = np.array(core_cells_hist)

# Event time estimate: maximum positive slope of Ec
dEc_dt = np.gradient(Ec_arr, t_arr)
i_evt = int(np.argmax(dEc_dt))
t_evt = t_arr[i_evt]

# Pre/post windows
w = max(2, len(t_arr)//10)
i0 = max(0, i_evt - w)
i1 = min(len(t_arr), i_evt + w)

Ec_pre = Ec_arr[i0:i_evt].mean() if i_evt > i0 else Ec_arr[i_evt]
Ec_post = Ec_arr[i_evt:i1].mean() if i1 > i_evt else Ec_arr[i_evt]
dEc_step = Ec_post - Ec_pre

core_pre = core_arr[i0:i_evt].mean() if i_evt > i0 else core_arr[i_evt]
core_post = core_arr[i_evt:i1].mean() if i1 > i_evt else core_arr[i_evt]
dcore = core_pre - core_post

# Wave speed estimate: track radius of peak shell signal vs time after event
Shell = np.array(Shell_Map)  # shape (Nt, Nr)
peak_idx = np.argmax(Shell, axis=1)
r_peak = R_BIN_CENTERS[peak_idx]

mask = (t_arr >= t_evt) & (r_peak > 0.0) & (r_peak < L/2)
t_fit = t_arr[mask]
r_fit = r_peak[mask]

v_fit = np.nan
if len(t_fit) >= 5:
    coef = np.polyfit(t_fit, r_fit, 1)
    v_fit = coef[0]

print("\n--- SIMULATION SUMMARY ---")
print(f"Grid: {N}^3, L={L}, dx={dx:.6f}, dt={DT}, samples={len(t_arr)}")
print(f"Core threshold rho_th={RHO_CORE_TH}")
print(f"Event time (max dEc/dt): t_evt ≈ {t_evt:.3f}")
print(f"Ec step (post-pre over window): ΔEc ≈ {dEc_step:.4e}")
print(f"Core proxy drop (pre-post):     Δ(core_cells) ≈ {dcore:.4e}")
print(f"Mass drift: max|m-m0| ≈ {np.max(np.abs(mass_arr - mass_arr[0])):.3e}")
if not np.isnan(v_fit):
    print(f"Estimated pulse speed from shell peak: v ≈ {v_fit:.3f} (expected c_s ~ 1 in these units)")
else:
    print("Pulse speed fit: insufficient data / unstable peak tracking")

# Save data for reproducibility
np.savez("gpe_diagnostics.npz",
         t=t_arr, Ec=Ec_arr, Ei=Ei_arr, mass=mass_arr, core_cells=core_arr,
         shell_map=Shell, r_centers=R_BIN_CENTERS, t_evt=t_evt)
print("Saved: gpe_diagnostics.npz")
